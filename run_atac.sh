#!/bin/bash

# settings
source atac.settings.config


# QC and clean the raw sequence data
start_time=$(date +%s)

if [ ! -d ${cleandata_dir} ]; then
	mkdir -p ${cleandata_dir}
fi

for sample in ${samples[*]}; do
	${run_cmd} fastp --in1 ${rawdata_dir}/${sample}.R1.fq.gz --in2 ${rawdata_dir}/${sample}.R2.fq.gz \
					--out1 ${cleandata_dir}/${sample}.R1.clean.fq.gz --out2 ${cleandata_dir}/${sample}.R2.clean.fq.gz \
					-R ${sample} -h ${cleandata_dir}/${sample}.qc.html -j ${cleandata_dir}/${sample}.qc.json \
					--thread ${threads}
done

end_time=$(date +%s)
running_time=$((${end_time}-${start_time}))
echo -e "\nDONE QC and clean.\nTime taken to execute commands is ${running_time} seconds.\n"


# mapping
start_time=$(date +%s)

if [ ! -d ${mapping_dir} ]; then
	mkdir -p ${mapping_dir}
fi

for sample in ${samples[*]}; do
	if [ "${aligner}" = "bowtie2" ]; then
		${run_cmd} bowtie2 -x ${reference_bowtie2} --threads ${threads} \
						${bowtie2_para} \
						-1 ${cleandata_dir}/${sample}.R1.clean.fq.gz \
						-2 ${cleandata_dir}/${sample}.R2.clean.fq.gz \
						-S ${mapping_dir}/${sample}.sam
	fi
done

end_time=$(date +%s)
running_time=$((${end_time}-${start_time}))
echo -e "\nDONE mapping.\nTime taken to execute commands is ${running_time} seconds.\n"


# select unique
start_time=$(date +%s)

for sample in ${samples[*]}; do
	${run_cmd} python ${selectUnique} -i ${mapping_dir}/${sample}.sam -o ${mapping_dir}/${sample}.unique.sam -p
	rm ${mapping_dir}/${sample}.sam
done

end_time=$(date +%s)
running_time=$((${end_time}-${start_time}))
echo -e "\nDONE selecting unique reads.\nTime taken to execute commands is ${running_time} seconds.\n"


# remove duplicates
start_time=$(date +%s)

for sample in ${samples[*]}; do
	${run_cmd} picard SortSam \
						INPUT=${mapping_dir}/${sample}.unique.sam \
						OUTPUT=${mapping_dir}/${sample}.unique.sortByQuery.sam \
						SORT_ORDER=queryname
	rm ${mapping_dir}/${sample}.unique.sam
	${run_cmd} picard MarkDuplicates \
						INPUT=${mapping_dir}/${sample}.unique.sortByQuery.sam \
						OUTPUT=${mapping_dir}/${sample}.unique.rmDup.sam \
						REMOVE_DUPLICATES=true \
						METRICS_FILE=${mapping_dir}/${sample}.unique.rmDup.metrics
	rm ${mapping_dir}/${sample}.unique.sortByQuery.sam

	
done

end_time=$(date +%s)
running_time=$((${end_time}-${start_time}))
echo -e "\nDONE removing duplicate reads.\nTime taken to execute commands is ${running_time} seconds.\n"


# bam correlation
start_time=$(date +%s)

sample_qc_files=""
sample_qc_labels=""
for sample in ${samples[*]}; do
	${run_cmd} samtools view -@ ${threads} \
							-bS ${mapping_dir}/${sample}.unique.rmDup.sam \
							-o ${mapping_dir}/${sample}.unique.rmDup.bam
	rm ${mapping_dir}/${sample}.unique.rmDup.sam
	${run_cmd} samtools sort -@ ${threads} ${mapping_dir}/${sample}.unique.rmDup.bam \
							-o ${mapping_dir}/${sample}.unique.rmDup.sorted.bam
	rm ${mapping_dir}/${sample}.unique.rmDup.bam
	${run_cmd} samtools index ${mapping_dir}/${sample}.unique.rmDup.sorted.bam

	sample_qc_files=${sample_qc_files}"${mapping_dir}/${sample}.unique.rmDup.sorted.bam "
	sample_qc_labels=${sample_qc_labels}"${sample} "
done

${run_cmd} multiBamSummary bins \
				--bamfiles ${sample_qc_files} \
				--labels ${sample_qc_labels} \
				--numberOfProcessors ${threads} \
				-o ${mapping_dir}/multiBamSummary_result.npz \
				--outRawCounts ${mapping_dir}/multiBamSummary_readCounts.tab

${run_cmd} plotCorrelation -in ${mapping_dir}/multiBamSummary_result.npz \
				--corMethod pearson --skipZeros --removeOutliers \
				--plotTitle "Pearson Correlation" \
				--whatToPlot heatmap \
				-o ${mapping_dir}/heatmap_PearsonCorr.png \
				--outFileCorMatrix ${mapping_dir}/PearsonCorr.tab

${run_cmd} plotCorrelation -in ${mapping_dir}/multiBamSummary_result.npz \
				--corMethod spearman --skipZeros \
				--plotTitle "Spearman Correlation" \
				--whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
				-o ${mapping_dir}/heatmap_SpearmanCorr.png \
				--outFileCorMatrix ${mapping_dir}/SpearmanCorr.tab

${run_cmd} plotPCA -in ${mapping_dir}/multiBamSummary_result.npz \
				-o ${mapping_dir}/PCA.png \
				-T "PCA"

end_time=$(date +%s)
running_time=$((${end_time}-${start_time}))
echo -e "\nDONE calculating sample correlation.\nTime taken to execute commands is ${running_time} seconds.\n"


# peak calling
start_time=$(date +%s)

if [ ! -d ${atac_dir}/${macs2_dir} ]; then
	mkdir -p ${atac_dir}/${macs2_dir}
fi

if [ ! -d ${atac_dir}/${hmmratac_dir} ]; then
	mkdir -p ${atac_dir}/${hmmratac_dir}
fi

for sample in ${samples[*]}; do
	${run_cmd} macs2 callpeak -t ${mapping_dir}/${sample}.unique.rmDup.sorted.bam -f BAMPE \
						-n ${atac_dir}/${macs2_dir}/${sample} \
						${macs2_para}

	${run_cmd} java -jar ${hmmratac} \
							-b ${mapping_dir}/${sample}.unique.rmDup.sorted.bam \
							-i ${mapping_dir}/${sample}.unique.rmDup.sorted.bam.bai \
							-g ${reference_size} \
							-o ${atac_dir}/${hmmratac_dir}/${sample} \
							${hmmratac_para}
done

end_time=$(date +%s)
running_time=$((${end_time}-${start_time}))
echo -e "\nDONE peak calling.\nTime taken to execute commands is ${running_time} seconds.\n"

