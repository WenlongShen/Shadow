#-*- coding: UTF-8 -*- 
#!/usr/bin/env python

import sys, getopt

def usage():
	print("hello!\n")

def processFile(sampleName, inputStat, inputRmdup, inputPbc, outputFile):
	
	inFileStat  = open(inputStat, "r")
	inFileRmdup = open(inputRmdup, "r")
	inFilePbc   = open(inputPbc, "r")

	outFile = open(outputFile, "a")

	outFile.write("%s\t" % (sampleName))

	line = inFileStat.readline()
	line = inFileStat.readline().strip()
	lineArr = line.split("\t")
	outFile.write("%d\t%d\t" % (int(lineArr[0]), int(lineArr[1])+int(lineArr[2])))

	for i in range(0,7):
		line = inFileRmdup.readline()
	line = inFileRmdup.readline().strip()
	lineArr = line.split("\t")
	outFile.write("%d\t%d\t" % (int(lineArr[2]), int(lineArr[2])-int(lineArr[6])))

	line = inFilePbc.readline().strip()
	lineArr = line.split("\t")
	outFile.write("%f\t%f\n" % (float(lineArr[3]), float(lineArr[4])))


if __name__ == "__main__":
	opts, args = getopt.getopt(sys.argv[1:], "hn:o:", ["i1=", "i2=", "i3="])
	sampleName = ""
	inputStat = ""
	inputRmdup = ""
	inputPbc = ""
	outputFile = ""
	for op, value in opts:
		if op == "-n":
			sampleName = value
		elif op == "-o":
			outputFile = value
		elif op == "--i1":
			inputStat = value
		elif op == "--i2":
			inputRmdup = value
		elif op == "--i3":
			inputPbc = value
		elif op == "-h":
			usage()
			sys.exit()

	processFile(sampleName, inputStat, inputRmdup, inputPbc, outputFile)
