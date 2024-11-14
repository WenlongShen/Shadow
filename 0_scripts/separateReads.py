#-*- coding: UTF-8 -*- 
#!/usr/bin/env python

import sys, getopt

def usage():
	print("hello!\n")

def processFile(inputFile, outputFile, paired):
	
	inFile = open(inputFile, "r")

	# unique: unique and in chroms (based on mapq)
	# multi: multi and in chroms
	# chrM: unique and multi in chrM
	# otherChr: unique and multi not in chroms & chrM
	# unmap: unmap or not in the same chrom
	outTypes = [".unique", ".multi", ".chrM", ".chrOthers", ".unmap"]
	outFile0 = open(outputFile+outTypes[0]+".sam", "w")
	outFile1 = open(outputFile+outTypes[1]+".sam", "w")
	outFile2 = open(outputFile+outTypes[2]+".sam", "w")
	outFile3 = open(outputFile+outTypes[3]+".sam", "w")
	outFile4 = open(outputFile+outTypes[4]+".sam", "w")
	outFile5 = open(outputFile+".mapping.stat", "w")

	outFile5.write("Total\tUnique\tMulti\tChrM\tChrOthers\tUnmap\n")
	count = [0,0,0,0,0,0]

	searchTerms = ["@", "XS", "AS", "chr", "_"]
	mapq_thresh = 30
	chroms = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", \
				"chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", \
				"chr21", "chr22", "chrX", "chrY"]
	chromM = ["chrM"]

	line1 = inFile.readline()
	while line1:
		if line1.find(searchTerms[0]) == 0:
			outFile0.write("%s" % line1)
			outFile1.write("%s" % line1)
			outFile2.write("%s" % line1)
			outFile3.write("%s" % line1)
			outFile4.write("%s" % line1)

		else:
			if paired == 1:
				line2 = inFile.readline()
				line1Arr = line1.split("\t")
				line2Arr = line2.split("\t")
				if line1Arr[0] != line2Arr[0]:
					print("unpaired error...\n")
					break
				count[0] += 1
				
				if (line1Arr[2] != line2Arr[2]) or line1Arr[2] == "*" or line2Arr[2] == "*":
					# unmap
					outFile4.write("%s%s" % (line1, line2))
					count[5] += 1
				else:
					inChroms1 = 0
					inChroms2 = 0
					for chrom in chroms:
						if line1Arr[2] == chrom:
							inChroms1 = 1
						if line2Arr[2] == chrom:
							inChroms2 = 1

					for chrom in chromM:
						if line1Arr[2] == chrom:
							inChroms1 = 2
						if line2Arr[2] == chrom:
							inChroms2 = 2

					if inChroms1 == 1 and inChroms2 == 1:
						if line1.find(searchTerms[2]) != -1 and line2.find(searchTerms[2]) != -1:
							# unique

							# if (line1.find(searchTerms[1]) == -1 \
							# 	and line1.find(searchTerms[2]) != -1) \
							# 	and (line2.find(searchTerms[1]) == -1 \
							# 	and line2.find(searchTerms[2]) != -1):
							if int(line1Arr[4]) > mapq_thresh and int(line2Arr[4]) > mapq_thresh:
								outFile0.write("%s%s" % (line1, line2))
								count[1] += 1

							# multi

							# elif (line1.find(searchTerms[1]) != -1 \
							# 	and line1.find(searchTerms[2]) != -1) \
							# 	or (line2.find(searchTerms[1]) != -1 \
							# 	and line2.find(searchTerms[2]) != -1):
							else:
								outFile1.write("%s%s" % (line1, line2))
								count[2] += 1

						else:
							outFile4.write("%s%s" % (line1, line2))
							count[5] += 1

					elif inChroms1 == 2 and inChroms2 == 2:
						# chrM
						outFile2.write("%s%s" % (line1, line2))
						count[3] += 1

					else:
						# otherChr
						outFile3.write("%s%s" % (line1, line2))
						count[4] += 1
					
			else:
				line1Arr = line1.split("\t")
				count[0] += 1

				if line1Arr[2] == "*":
					# unmap
					outFile4.write("%s" % (line1))
					count[5] += 1
				else:
					inChroms1 = 0
					for chrom in chroms:
						if line1Arr[2] == chrom:
							inChroms1 = 1

					for chrom in chromM:
						if line1Arr[2] == chrom:
							inChroms1 = 2

					if inChroms1 == 1:
						# unique
						if (line1.find(searchTerms[1]) == -1 \
							and line1.find(searchTerms[2]) != -1):
							outFile0.write("%s" % (line1))
							count[1] += 1

						# multi
						elif (line1.find(searchTerms[1]) != -1 \
							and line1.find(searchTerms[2]) != -1):
							outFile1.write("%s" % (line1))
							count[2] += 1

						else:
							outFile4.write("%s" % (line1))
							count[5] += 1

					elif inChroms1 == 2:
						# chrM
						outFile2.write("%s" % (line1))
						count[3] += 1

					else:
						# otherChr
						outFile3.write("%s" % (line1))
						count[4] += 1

		line1 = inFile.readline()

	if count[0] != count[1]+count[2]+count[3]+count[4]+count[5]:
		print("counting error...\n")
	outFile5.write("%d\t%d\t%d\t%d\t%d\t%d\n" % (count[0],count[1],count[2],count[3],count[4],count[5]))



if __name__ == "__main__":
	opts, args = getopt.getopt(sys.argv[1:], "hi:o:", ["pe", "se"])
	paired = 0
	inputFile = ""
	outputFile = ""
	for op, value in opts:
		if op == "-i":
			inputFile = value
		elif op == "-o":
			outputFile = value
		elif op == "--pe":
			paired = 1
		elif op == "--se":
			paired = 0
		elif op == "-h":
			usage()
			sys.exit()

	processFile(inputFile, outputFile, paired)
