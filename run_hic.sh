#!/bin/bash

# settings
source hic.settings.config


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
					-U ${cleandata_dir}/${sample}.R1.clean.fq.gz \
					${bowtie2_para} \
					| samtools view -Shb - > ${mapping_dir}/${sample}.R1.bam
		${run_cmd} bowtie2 -x ${reference_bowtie2} --threads ${threads} \
					-U ${cleandata_dir}/${sample}.R2.clean.fq.gz \
					${bowtie2_para} \
					| samtools view -Shb - > ${mapping_dir}/${sample}.R2.bam
	fi
done

end_time=$(date +%s)
running_time=$((${end_time}-${start_time}))
echo -e "\nDONE mapping.\nTime taken to execute commands is ${running_time} seconds.\n"


# hicexplorer
start_time=$(date +%s)

if [ ! -d ${hic_dir} ]; then
	mkdir -p ${hic_dir}
fi
if [ ! -d ${hic_dir}/${hic_raw_dir} ]; then
	mkdir -p ${hic_dir}/${hic_raw_dir}
fi

${run_cmd} hicFindRestSite --fasta ${reference_fa} --searchPattern ${restriction_site} \
					-o ${hic_dir}/${hic_raw_dir}/${enzyme}_cut_site_${reference_version}.bed

end_time=$(date +%s)
running_time=$((${end_time}-${start_time}))
echo -e "\nDONE preprocessing.\nTime taken to execute commands is ${running_time} seconds.\n"


# create raw Hi-C matrices
start_time=$(date +%s)

for resolution in ${resolutions[*]}; do
	if [ ! -d ${hic_dir}/${hic_raw_dir}/${resolution} ]; then
		mkdir -p ${hic_dir}/${hic_raw_dir}/${resolution}
	fi

	for sample in ${samples[*]}; do
		${run_cmd} hicBuildMatrix --samFiles ${mapping_dir}/${sample}.R1.bam ${mapping_dir}/${sample}.R2.bam \
						--binSize ${resolution} \
						--restrictionSequence ${restriction_site} --danglingSequence ${dangling_end} \
						--restrictionCutFile ${hic_dir}/${hic_raw_dir}/${enzyme}_cut_site_${reference_version}.bed \
						--chromosomeSizes ${hic_dir}/${hic_raw_dir}/${reference_version}.chrom.sizes \
						--threads ${threads} --inputBufferSize 100000 \
						-o ${hic_dir}/${hic_raw_dir}/${resolution}/${sample}.${resolution}.cool \
						--outBam ${hic_dir}/${hic_raw_dir}/${resolution}/${sample}.${resolution}.bam \
						--QCfolder ${hic_dir}/${hic_raw_dir}/${resolution}/${sample}_QC
	done

	sample_qc_labels=""
	sample_qc_files=""
	sample_qc_matrices=""
	for sample in ${samples[*]}; do
		sample_qc_labels=${sample_qc_labels}"${sample} "
		sample_qc_files=${sample_qc_files}"${hic_dir}/${hic_raw_dir}/${resolution}/${sample}_QC/QC.log "
		sample_qc_matrices=${sample_qc_matrices}"${hic_dir}/${hic_raw_dir}/${resolution}/${sample}.${resolution}.cool "
	done
	${run_cmd} hicQC --logfiles ${sample_qc_files} --labels ${sample_qc_labels} \
		--outputFolder ${hic_dir}/${hic_raw_dir}/${resolution}/All_samples_QC
	${run_cmd} hicCorrelate -m ${sample_qc_matrices} \
				--method=pearson --log1p \
				--labels ${sample_qc_labels} \
				--outFileNameHeatmap ${hic_dir}/${hic_raw_dir}/${resolution}/${resolution}_pearson_heatmap \
				--outFileNameScatter ${hic_dir}/${hic_raw_dir}/${resolution}/${resolution}_pearson_scatterplot \
				--plotFileFormat pdf --threads ${threads}
	${run_cmd} hicCorrelate -m ${sample_qc_matrices} \
				--method=spearman --log1p \
				--labels ${sample_qc_labels} \
				--outFileNameHeatmap ${hic_dir}/${hic_raw_dir}/${resolution}/${resolution}_spearman_heatmap \
				--outFileNameScatter ${hic_dir}/${hic_raw_dir}/${resolution}/${resolution}_spearman_scatterplot \
				--plotFileFormat pdf --threads ${threads}

	# create merge (sum) matrices from replicates
	for sample in ${sample_names[*]}; do
		sample_matrices=""
		for rep in ${rep_numbers[*]}; do
			if [[ ${sample} = "B" ]]; then
				sample_matrices=${sample_matrices}"${hic_dir}/${hic_raw_dir}/${resolution}/${sample}-${rep}.${resolution}.cool "
			else
				if [[ ${rep} -ne "3" ]]; then
					sample_matrices=${sample_matrices}"${hic_dir}/${hic_raw_dir}/${resolution}/${sample}-${rep}.${resolution}.cool "
				fi
			fi
		done
		${run_cmd} hicSumMatrices -m ${sample_matrices} \
						-o ${hic_dir}/${hic_raw_dir}/${resolution}/${sample}.all.${resolution}.cool
	done

done

end_time=$(date +%s)
running_time=$((${end_time}-${start_time}))
echo -e "\nDONE creation of raw Hi-C matrices.\nTime taken to execute commands is ${running_time} seconds.\n"


# normalize and correct matrices
start_time=$(date +%s)

if [ ! -d ${hic_dir}/${hic_corrected_dir} ]; then
	mkdir -p ${hic_dir}/${hic_corrected_dir}
fi

for resolution in ${resolutions[*]}; do
	if [ ! -d ${hic_dir}/${hic_corrected_dir}/${resolution} ]; then
		mkdir -p ${hic_dir}/${hic_corrected_dir}/${resolution}
	fi

	sample_hic_raw=""
	sample_hic_normalized=""
	for sample in ${sample_names[*]}; do
		sample_hic_raw=${sample_hic_raw}"${hic_dir}/${hic_raw_dir}/${resolution}/${sample}.all.${resolution}.cool "
		sample_hic_normalized=${sample_hic_normalized}"${hic_dir}/${hic_corrected_dir}/${resolution}/${sample}.all.${resolution}.${normalize_method}.cool "
	done

	${run_cmd} hicNormalize -m ${sample_hic_raw} \
					-o ${sample_hic_normalized} \
					--normalize ${normalize_method}	\
					--setToZeroThreshold ${setToZeroThreshold} 	
done

end_time=$(date +%s)
running_time=$((${end_time}-${start_time}))
echo -e "\nDONE normaliztion of raw Hi-C matrices.\nTime taken to execute commands is ${running_time} seconds.\n"


# normalize and correct matrices
start_time=$(date +%s)

if [ ! -d ${hic_dir}/${hic_corrected_dir} ]; then
	mkdir -p ${hic_dir}/${hic_corrected_dir}
fi

for resolution in ${resolutions[*]}; do
	if [ ! -d ${hic_dir}/${hic_corrected_dir}/${resolution} ]; then
		mkdir -p ${hic_dir}/${hic_corrected_dir}/${resolution}
	fi

	for sample in ${sample_names[*]}; do
		${run_cmd} hicCorrectMatrix diagnostic_plot -m ${hic_dir}/${hic_corrected_dir}/${resolution}/${sample}.all.${resolution}.${normalize_method}.cool \
						-o ${hic_dir}/${hic_corrected_dir}/${resolution}/${sample}.all.${resolution}.${normalize_method}.diagnostic_plot.png
		${run_cmd} hicCorrectMatrix correct -m ${hic_dir}/${hic_corrected_dir}/${resolution}/${sample}.all.${resolution}.${normalize_method}.cool \
						--chromosomes ${chrom} \
						--correctionMethod ICE \
						--filterThreshold -1.5 3 \
						--outFileName ${hic_dir}/${hic_corrected_dir}/${resolution}/${sample}.all.${resolution}.${normalize_method}.corrected.cool
	done
done

end_time=$(date +%s)
running_time=$((${end_time}-${start_time}))
echo -e "\nDONE correction of normalized Hi-C matrices.\nTime taken to execute commands is ${running_time} seconds.\n"


# AB compartments analysis
start_time=$(date +%s)

if [ ! -d ${hic_dir}/${hic_domain_dir} ]; then
	mkdir -p ${hic_dir}/${hic_domain_dir}
fi
if [ ! -d ${hic_dir}/${hic_domain_dir}/${hic_ab_dir} ]; then
	mkdir -p ${hic_dir}/${hic_domain_dir}/${hic_ab_dir}
fi

for resolution in ${resolutions[*]}; do
	if [ ! -d ${hic_dir}/${hic_domain_dir}/${hic_ab_dir}/${resolution} ]; then
		mkdir -p ${hic_dir}/${hic_domain_dir}/${hic_ab_dir}/${resolution}
	fi
	for sample in ${sample_names[*]}; do
		${run_cmd} hicPCA -m ${hic_dir}/${hic_corrected_dir}/${resolution}/${sample}.all.${resolution}.${normalize_method}.corrected.cool \
				-o ${hic_dir}/${hic_domain_dir}/${hic_ab_dir}/${resolution}/${sample}.all.${resolution}.${normalize_method}.corrected.pca1.bw \
				${hic_dir}/${hic_domain_dir}/${hic_ab_dir}/${resolution}/${sample}.all.${resolution}.${normalize_method}.corrected.pca2.bw \
				-f bigwig
	done
done

end_time=$(date +%s)
running_time=$((${end_time}-${start_time}))
echo -e "\nDONE AB compartments calling.\nTime taken to execute commands is ${running_time} seconds.\n"


# tad calling
start_time=$(date +%s)

if [ ! -d ${hic_dir}/${hic_domain_dir}/${hic_tad_dir} ]; then
	mkdir -p ${hic_dir}/${hic_domain_dir}/${hic_tad_dir}
fi

for resolution in ${resolutions[*]}; do
	if [ ! -d ${hic_dir}/${hic_domain_dir}/${hic_tad_dir}/${resolution} ]; then
		mkdir -p ${hic_dir}/${hic_domain_dir}/${hic_tad_dir}/${resolution}
	fi
	for sample in ${sample_names[*]}; do
		${run_cmd} hicFindTADs -m ${hic_dir}/${hic_corrected_dir}/${resolution}/${sample}.all.${resolution}.${normalize_method}.corrected.cool \
				--outPrefix ${hic_dir}/${hic_domain_dir}/${hic_tad_dir}/${resolution}/${sample}.all.${resolution}.${normalize_method}.corrected.tad \
				--minDepth $[3*${resolution}] \
				--maxDepth $[10*${resolution}] \
				--step ${resolution} \
				--thresholdComparisons 0.05 \
				--delta 0.01 \
				--correctForMultipleTesting fdr \
				-p ${threads}
	done
done

end_time=$(date +%s)
running_time=$((${end_time}-${start_time}))
echo -e "\nDONE TADs calling.\nTime taken to execute commands is ${running_time} seconds.\n"


# loop detection
start_time=$(date +%s)

if [ ! -d ${hic_dir}/${hic_domain_dir}/${hic_loop_dir} ]; then
	mkdir -p ${hic_dir}/${hic_domain_dir}/${hic_loop_dir}
fi

for resolution in ${resolutions[*]}; do
	if [ ! -d ${hic_dir}/${hic_domain_dir}/${hic_loop_dir}/${resolution} ]; then
		mkdir -p ${hic_dir}/${hic_domain_dir}/${hic_loop_dir}/${resolution}
	fi
	for sample in ${sample_names[*]}; do
		${run_cmd} hicDetectLoops -m ${hic_dir}/${hic_corrected_dir}/${resolution}/${sample}.all.${resolution}.${normalize_method}.corrected.cool \
				-o ${hic_dir}/${hic_domain_dir}/${hic_loop_dir}/${resolution}/${sample}.all.${resolution}.${normalize_method}.corrected.loops.bedgraph \
				--pValuePreselection 0.05 --pValue 0.01 \
				--threads ${threads}
	done
done

end_time=$(date +%s)
running_time=$((${end_time}-${start_time}))
echo -e "\nDONE loops calling.\nTime taken to execute commands is ${running_time} seconds.\n"
