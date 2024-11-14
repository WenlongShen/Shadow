#!/bin/bash

# settings
source settings_atac.config


# QC and clean the raw sequence data
start_time=$(date +%s)

if [ ! -d ${cleandata_dir} ]; then
	mkdir -p ${cleandata_dir}
fi

for sample in ${samples[*]}; do
	${run_cmd} fastp --in1 ${rawdata_dir}/${sample}_R1.fastq.gz --in2 ${rawdata_dir}/${sample}_R2.fastq.gz \
					--out1 ${cleandata_dir}/${sample}.R1.clean.fastq.gz --out2 ${cleandata_dir}/${sample}.R2.clean.fastq.gz \
					-R ${sample} -h ${cleandata_dir}/${sample}.qc.html -j ${cleandata_dir}/${sample}.qc.json \
					--thread ${threads}
done

end_time=$(date +%s)
running_time=$((${end_time}-${start_time}))
echo -e "\nDONE QC and clean.\nTime taken to execute commands is ${running_time} seconds.\n"



# mapping
start_time=$(date +%s)

if [ ! -d ${bowtie2_dir} ]; then
	mkdir -p ${bowtie2_dir}
fi

for sample in ${samples[*]}; do
	if [ "${aligner}" = "bowtie2" ]; then
		${run_cmd} bowtie2 -x ${reference_bowtie2} --threads ${threads} \
						${bowtie2_para} \
						-1 ${cleandata_dir}/${sample}.R1.clean.fastq.gz \
						-2 ${cleandata_dir}/${sample}.R2.clean.fastq.gz \
						-S ${bowtie2_dir}/${sample}.sam \
						2> ${bowtie2_dir}/${sample}.mapping.log
	fi
done

end_time=$(date +%s)
running_time=$((${end_time}-${start_time}))
echo -e "\nDONE mapping.\nTime taken to execute commands is ${running_time} seconds.\n"



# separate reads
start_time=$(date +%s)

for sample in ${samples[*]}; do
	${run_cmd} python ${separateReads} -i ${bowtie2_dir}/${sample}.sam -o ${bowtie2_dir}/${sample} --pe
	${run_cmd} rm ${bowtie2_dir}/${sample}.sam
	${run_cmd} bedtools bamtobed -bedpe -i ${bowtie2_dir}/${sample}.unique.sam | \
		awk 'BEGIN{{OFS="\t"}} {{print $1,$2,$3,$9,$4,$5,$6,$10}}' | \
		bedtools sort -i | uniq -c | \
		awk 'BEGIN{{m0=0;m1=0;m2=0}} ($1==1){{m1=m1+1}} ($1==2){{m2=m2+1}} {{m0=m0+1}} END{{m1_m2=-1.0; if(m2>0) m1_m2=m1/m2; m1_m0=-1.0; if (m0>0) m1_m0=m1/m0; printf "%d\t%d\t%d\t%f\t%f\n", m0,m1,m2,m1_m0,m1_m2}}' > ${bowtie2_dir}/${sample}.mapping.pbc
done

end_time=$(date +%s)
running_time=$((${end_time}-${start_time}))
echo -e "\nDONE PBC counting.\nTime taken to execute commands is ${running_time} seconds.\n"



# remove duplicates
start_time=$(date +%s)

for sample in ${samples[*]}; do
	${run_cmd} picard SortSam \
						INPUT=${bowtie2_dir}/${sample}.unique.sam \
						OUTPUT=${bowtie2_dir}/${sample}.unique.sortByQuery.sam \
						SORT_ORDER=queryname
	${run_cmd} rm ${bowtie2_dir}/${sample}.unique.sam
	${run_cmd} picard MarkDuplicates \
						INPUT=${bowtie2_dir}/${sample}.unique.sortByQuery.sam \
						OUTPUT=${bowtie2_dir}/${sample}.unique.rmdup.sam \
						REMOVE_DUPLICATES=true \
						METRICS_FILE=${bowtie2_dir}/${sample}.unique.rmdup.metrics
	${run_cmd} rm ${bowtie2_dir}/${sample}.unique.sortByQuery.sam
done

end_time=$(date +%s)
running_time=$((${end_time}-${start_time}))
echo -e "\nDONE removing duplicate reads.\nTime taken to execute commands is ${running_time} seconds.\n"



# mapping summary
start_time=$(date +%s)

${run_cmd} echo -e "Samples\tTotal\tKeepChr\tUniqueMapping\tRmDup\tPBC1\tPBC2" > ${mapping_dir}/mappingSummary.txt
for sample in ${samples[*]}; do
	${run_cmd} python ${mappingSummary} -n ${sample} \
										--i1=${bowtie2_dir}/${sample}.mapping.stat \
										--i2=${bowtie2_dir}/${sample}.unique.rmdup.metrics \
										--i3=${bowtie2_dir}/${sample}.mapping.pbc \
										-o ${mapping_dir}/mappingSummary.txt
done

end_time=$(date +%s)
running_time=$((${end_time}-${start_time}))
echo -e "\nDONE mapping summary.\nTime taken to execute commands is ${running_time} seconds.\n"



# bam correlation
start_time=$(date +%s)

sample_qc_files=""
sample_qc_labels=""
for sample in ${samples[*]}; do
	${run_cmd} samtools view -@ ${threads} \
							-bS ${bowtie2_dir}/${sample}.unique.rmdup.sam \
							-o ${bowtie2_dir}/${sample}.unique.rmdup.bam
	${run_cmd} rm ${bowtie2_dir}/${sample}.unique.rmdup.sam
	${run_cmd} samtools sort -@ ${threads} ${bowtie2_dir}/${sample}.unique.rmdup.bam \
							-o ${bowtie2_dir}/${sample}.unique.rmdup.sorted.bam
	${run_cmd} rm ${bowtie2_dir}/${sample}.unique.rmdup.bam
	${run_cmd} samtools index ${bowtie2_dir}/${sample}.unique.rmdup.sorted.bam

	sample_qc_files=${sample_qc_files}"${bowtie2_dir}/${sample}.unique.rmdup.sorted.bam "
	sample_qc_labels=${sample_qc_labels}"${sample} "
done

${run_cmd} multiBamSummary bins \
				--bamfiles ${sample_qc_files} \
				--labels ${sample_qc_labels} \
				--numberOfProcessors ${threads} \
				-o ${bowtie2_dir}/multiBamSummary_result.npz \
				--outRawCounts ${bowtie2_dir}/multiBamSummary_readCounts.tab

${run_cmd} plotCorrelation -in ${bowtie2_dir}/multiBamSummary_result.npz \
				--corMethod pearson --skipZeros --removeOutliers \
				--plotTitle "Pearson Correlation" \
				--whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
				-o ${bowtie2_dir}/heatmap_PearsonCorr.png \
				--outFileCorMatrix ${bowtie2_dir}/PearsonCorr.tab

${run_cmd} plotCorrelation -in ${bowtie2_dir}/multiBamSummary_result.npz \
				--corMethod spearman --skipZeros --removeOutliers \
				--plotTitle "Spearman Correlation" \
				--whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
				-o ${bowtie2_dir}/heatmap_SpearmanCorr.png \
				--outFileCorMatrix ${bowtie2_dir}/SpearmanCorr.tab

${run_cmd} plotPCA -in ${bowtie2_dir}/multiBamSummary_result.npz \
				-o ${bowtie2_dir}/PCA.png \
				-T "PCA"

end_time=$(date +%s)
running_time=$((${end_time}-${start_time}))
echo -e "\nDONE calculating sample correlation.\nTime taken to execute commands is ${running_time} seconds.\n"



# tn5 shifting
start_time=$(date +%s)

if [ ! -d ${tn5shift_dir} ]; then
	mkdir -p ${tn5shift_dir}
fi

for sample in ${samples[*]}; do
	${run_cmd} alignmentSieve -b ${bowtie2_dir}/${sample}.unique.rmdup.sorted.bam \
							-o ${tn5shift_dir}/${sample}.unique.rmdup.sorted.tn5shift.bam \
							--ATACshift --numberOfProcessors ${threads}
	${run_cmd} samtools sort -@ ${threads} ${tn5shift_dir}/${sample}.unique.rmdup.sorted.tn5shift.bam \
							-o ${tn5shift_dir}/${sample}.unique.rmdup.sorted.tn5shift.sorted.bam
	${run_cmd} rm ${tn5shift_dir}/${sample}.unique.rmdup.sorted.tn5shift.bam
done

end_time=$(date +%s)
running_time=$((${end_time}-${start_time}))
echo -e "\nDONE tn5 shifting.\nTime taken to execute commands is ${running_time} seconds.\n"



# bam merging
start_time=$(date +%s)

if [ ! -d ${mergedBam_dir} ]; then
	mkdir -p ${mergedBam_dir}
fi

for sample in ${sample_names[*]}; do
	${run_cmd} samtools merge ${mergedBam_dir}/${sample}.bam ${tn5shift_dir}/${sample}_1.unique.rmdup.sorted.tn5shift.sorted.bam ${tn5shift_dir}/${sample}_2.unique.rmdup.sorted.tn5shift.sorted.bam -@ ${threads}
	${run_cmd} samtools index ${mergedBam_dir}/${sample}.bam -@ ${threads}
	${run_cmd} bamCoverage -b ${mergedBam_dir}/${sample}.bam -o ${mergedBam_dir}/${sample}.bigWig -of bigwig --binSize 1 --ignoreDuplicates --extendReads --normalizeUsing BPM --effectiveGenomeSize ${reference_genomesize} -p ${threads}
done

end_time=$(date +%s)
running_time=$((${end_time}-${start_time}))
echo -e "\nDONE bam merging.\nTime taken to execute commands is ${running_time} seconds.\n"



# peak calling
start_time=$(date +%s)

if [ ! -d ${macs3_dir} ]; then
	mkdir -p ${macs3_dir}
fi

for sample in ${samples[*]}; do
	${run_cmd} macs3 callpeak -t ${tn5shift_dir}/${sample}.unique.rmdup.sorted.tn5shift.sorted.bam \
						-n ${macs3_dir}/${sample} \
						${macs3_para}
done


if [ ! -d ${chipr_dir} ]; then
	mkdir -p ${chipr_dir}
fi

sample_chipr_files=""
for sample in ${cell_names[*]}; do
	sample_chipr_files=""
	for rep in ${rep_numbers[*]}; do
		sample_chipr_files=${sample_chipr_files}"${macs3_dir}/${sample}-I_${rep}_peaks.narrowPeak ${macs3_dir}/${sample}-M_${rep}_peaks.narrowPeak"
	done
	${run_cmd} chipr -i ${sample_chipr_files} -o ${chipr_dir}/${sample}.cache.txt ${chipr_para}
done

end_time=$(date +%s)
running_time=$((${end_time}-${start_time}))
echo -e "\nDONE peak calling.\nTime taken to execute commands is ${running_time} seconds.\n"

