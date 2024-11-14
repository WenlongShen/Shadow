#!/bin/bash

# settings
source settings_rna.config


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

sample_qc_files=""
sample_qc_labels=""
for sample in ${samples[*]}; do
	if [ "${aligner}" = "hisat2" ]; then
		${run_cmd} hisat2 -x ${reference_hisat2} -p ${threads} --dta \
						-1 ${cleandata_dir}/${sample}.R1.clean.fq.gz \
						-2 ${cleandata_dir}/${sample}.R2.clean.fq.gz \
						-S ${mapping_dir}/${sample}.sam
	fi

	${run_cmd} samtools view -@ ${threads} -bS ${mapping_dir}/${sample}.sam -o ${mapping_dir}/${sample}.bam
	rm ${mapping_dir}/${sample}.sam
	${run_cmd} samtools sort -@ ${threads} ${mapping_dir}/${sample}.bam -o ${mapping_dir}/${sample}.sorted.bam
	rm ${mapping_dir}/${sample}.bam
	${run_cmd} samtools index ${mapping_dir}/${sample}.sorted.bam

	sample_qc_files=${sample_qc_files}"${mapping_dir}/${sample}.sorted.bam "
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
echo -e "\nDONE mapping.\nTime taken to execute commands is ${running_time} seconds.\n"


# stringtie
start_time=$(date +%s)

if [ ! -d ${stringtie_dir} ]; then
	mkdir -p ${stringtie_dir}
fi

for sample in ${samples[*]}; do
	${run_cmd} stringtie ${mapping_dir}/${sample}.sorted.bam \
					-e -G ${reference_gtf} \
					-o ${stringtie_dir}/${sample}.stringtie.gtf \
					-A ${stringtie_dir}/${sample}.abund.txt \
					-p ${threads}
done

${run_cmd} python ${prepDE} -i ${stringtie_dir} -g ${stringtie_dir}/gene_count_matrix.csv -t ${stringtie_dir}/transcript_count_matrix.csv

end_time=$(date +%s)
running_time=$((${end_time}-${start_time}))
echo -e "\nDONE stringtie.\nTime taken to execute commands is ${running_time} seconds.\n"

