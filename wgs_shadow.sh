#!/bin/bash

# settings
source settings.config

samples=$(echo -e ${samples} | tr "," "\n")

# QC and clean data
if [ ! -d ${cleandata_dir} ]; then
	mkdir -p ${cleandata_dir}
fi

start_time=$(date +%s)

for sample in ${samples}; do
	${fastp} --in1 ${rawdata_dir}/${sample}.R1.fq.gz --in2 ${rawdata_dir}/${sample}.R2.fq.gz \
		--out1 ${cleandata_dir}/${sample}.R1.clean.fq.gz --out2 ${cleandata_dir}/${sample}.R2.clean.fq.gz \
		-R ${sample} -h ${cleandata_dir}/${sample}.qc.html -j ${cleandata_dir}/${sample}.qc.json
done

end_time=$(date +%s)
running_time=$((${end_time}-${start_time}))
echo -e "\nDONE QC and clean.\nTime taken to execute commands is ${running_time} seconds.\n"


# mapping
if [ ! -d ${mapping_dir} ]; then
	mkdir -p ${mapping_dir}
fi

start_time=$(date +%s)

for sample in ${samples}; do
	if [ "${aligner}" = "bwa" ]; then
		ID=${sample}_RG # Read Group, or Lane ID
		PL=ILLUMINA # ILLUMINA, SLX, SOLEXA, SOLID, 454, LS454, COMPLETE, PACBIO, IONTORRENT, CAPILLARY, HELICOS, UNKNOWN
		LB=${sample}_lib # Library
		SM=${sample} # Sample
		${bwa} mem -t ${threads} -M \
			-R "@RG\tID:$ID\tPL:$PL\tLB:$LB\tSM:$SM" \
			${reference_bwa} \
			${cleandata_dir}/${sample}.R1.clean.fq.gz ${cleandata_dir}/${sample}.R2.clean.fq.gz \
			| samtools view -bS - \
			| samtools sort -@ ${threads} -o ${mapping_dir}/${sample}.bam
	fi
done

end_time=$(date +%s)
running_time=$((${end_time}-${start_time}))
echo -e "\nDONE mapping.\nTime taken to execute commands is ${running_time} seconds.\n"


# mark/remove duplicates
start_time=$(date +%s)

for sample in ${samples}; do
	${gatk} MarkDuplicatesSpark \
		-I ${mapping_dir}/${sample}.bam \
		-O ${mapping_dir}/${sample}.markdup.bam \
		-M ${mapping_dir}/${sample}.markdup.metrics \
		--remove-all-duplicates false \
		--spark-verbosity WARN
done

end_time=$(date +%s)
running_time=$((${end_time}-${start_time}))
echo -e "\nDONE mark/remove duplicates.\nTime taken to execute commands is ${running_time} seconds.\n"


# fix tags
start_time=$(date +%s)

for sample in ${samples}; do
	${gatk} FixMateInformation \
		-I ${mapping_dir}/${sample}.markdup.bam \
		-O ${mapping_dir}/${sample}.markdup.fixmate.bam \
		--CREATE_INDEX true \
		--VERBOSITY WARNING
	${gatk} SetNmMdAndUqTags \
		-R ${reference_fa} \
		-I ${mapping_dir}/${sample}.markdup.fixmate.bam \
		-O ${mapping_dir}/${sample}.markdup.fixmate.fixtags.bam \
		--CREATE_INDEX true \
		--VERBOSITY WARNING
done

end_time=$(date +%s)
running_time=$((${end_time}-${start_time}))
echo -e "\nDONE fix tags.\nTime taken to execute commands is ${running_time} seconds.\n"


# BQSR
if [ ! -d ${bqsr_dir} ]; then
	mkdir -p ${bqsr_dir}
fi

start_time=$(date +%s)

for sample in ${samples}; do
	${gatk} BaseRecalibrator \
		-R ${reference_fa} \
		--known-sites ${dbsnp} \
		--known-sites ${G1000} \
		--known-sites ${mills} \
		-I ${mapping_dir}/${sample}.markdup.fixmate.fixtags.bam \
		-O ${bqsr_dir}/${sample}.bqsr.recal.table
	${gatk} ApplyBQSR \
		-R ${reference_fa} \
		-I ${mapping_dir}/${sample}.markdup.fixmate.fixtags.bam \
		--bqsr-recal-file ${bqsr_dir}/${sample}.bqsr.recal.table \
		-O ${bqsr_dir}/${sample}.bqsr.bam \
		--create-output-bam-index true
done

end_time=$(date +%s)
running_time=$((${end_time}-${start_time}))
echo -e "\nDONE BQSR.\nTime taken to execute commands is ${running_time} seconds.\n"


# call germline SNPs and indels
# HaplotypeCaller
if [ ! -d ${hc_dir} ]; then
	mkdir -p ${hc_dir}
fi

start_time=$(date +%s)

for sample in ${samples}; do
	${gatk} HaplotypeCaller \
		-R ${reference_fa} \
		-ERC GVCF \
		-I ${bqsr_dir}/${sample}.bqsr.bam \
		-O ${hc_dir}/${sample}.hc.g.vcf.gz
done

end_time=$(date +%s)
running_time=$((${end_time}-${start_time}))
echo -e "\nDONE HaplotypeCaller.\nTime taken to execute commands is ${running_time} seconds.\n"


# call germline SNPs and Indels
# joint genotyping
start_time=$(date +%s)

sample_gvcfs=""
for sample in ${samples}; do
	sample_gvcfs=${sample_gvcfs}"-V ${hc_dir}/${sample}.hc.g.vcf.gz "
done

${gatk} CombineGVCFs \
	-R ${reference_fa} \
	${sample_gvcfs} \
	-O ${hc_dir}/${project_name}.hc.g.vcf.gz 
${gatk} GenotypeGVCFs \
	-R ${reference_fa} \
	-V ${hc_dir}/${project_name}.hc.g.vcf.gz \
	-O ${hc_dir}/${project_name}.hc.vcf.gz 

end_time=$(date +%s)
running_time=$((${end_time}-${start_time}))
echo -e "\nDONE joint genotyping.\nTime taken to execute commands is ${running_time} seconds.\n"


# VQSR
if [ "${filtration_vqsr}" = "true" ]; then
	if [ ! -d ${vqsr_dir} ]; then
		mkdir -p ${vqsr_dir}
	fi

	# call germline SNPs
	if [ "${call_snps}" = "true" ]; then
		start_time=$(date +%s)

		${gatk} VariantRecalibrator \
			-R ${reference_fa} \
			${resource_hapmap} \
			${resource_omni} \
			${resource_1000G} \
			${resource_dbsnp} \
			${snps_vqsr_params} \
			-V ${hc_dir}/${project_name}.hc.vcf.gz \
			-O ${vqsr_dir}/${project_name}.hc.snps.recal.table \
			--tranches-file ${vqsr_dir}/${project_name}.hc.snps.recal.tranches \
			--rscript-file ${vqsr_dir}/${project_name}.hc.snps.recal.plots.R
		${gatk} ApplyVQSR \
			-R ${reference_fa} \
			-V ${hc_dir}/${project_name}.hc.vcf.gz \
			--recal-file ${vqsr_dir}/${project_name}.hc.snps.recal.table \
			--tranches-file ${vqsr_dir}/${project_name}.hc.snps.recal.tranches \
			${ts_filter_level} \
			-mode SNP \
			-O ${vqsr_dir}/${project_name}.hc.snps.vqsr.vcf.gz \
			--create-output-variant-index true

		# output snps
		${gatk} SelectVariants \
			-R ${reference_fa} \
			-select-type SNP \
			-V ${vqsr_dir}/${project_name}.hc.snps.vqsr.vcf.gz \
			-O ${vqsr_dir}/${project_name}.hc.snps.vqsr_filt.vcf.gz

		end_time=$(date +%s)
		running_time=$((${end_time}-${start_time}))
		echo -e "\nDONE VQSR snps.\nTime taken to execute commands is ${running_time} seconds.\n"
	fi

	# VQSR
	# call germline Indels
	if [ "${call_indels}" = "true" ]; then
		start_time=$(date +%s)

		${gatk} VariantRecalibrator \
			-R ${reference_fa} \
			${resource_mills} \
			${resource_dbsnp} \
			${indels_vqsr_params} \
			-V ${hc_dir}/${project_name}.hc.vcf.gz \
			-O ${vqsr_dir}/${project_name}.hc.indels.recal.table \
			--tranches-file ${vqsr_dir}/${project_name}.hc.indels.recal.tranches \
			--rscript-file ${vqsr_dir}/${project_name}.hc.indels.recal.plots.R
		${gatk} ApplyVQSR \
			-R ${reference_fa} \
			-V ${hc_dir}/${project_name}.hc.vcf.gz \
			--recal-file ${vqsr_dir}/${project_name}.hc.indels.recal.table \
			--tranches-file ${vqsr_dir}/${project_name}.hc.indels.recal.tranches \
			${ts_filter_level} \
			-mode INDEL \
			-O ${vqsr_dir}/${project_name}.hc.indels.vqsr.vcf.gz \
			--create-output-variant-index true

		# output indels
		${gatk} SelectVariants \
			-R ${reference_fa} \
			-select-type INDEL \
			-V ${vqsr_dir}/${project_name}.hc.indels.vqsr.vcf.gz \
			-O ${vqsr_dir}/${project_name}.hc.indels.vqsr_filt.vcf.gz

		end_time=$(date +%s)
		running_time=$((${end_time}-${start_time}))
		echo -e "\nDONE VQSR indels.\nTime taken to execute commands is ${running_time} seconds.\n"
	fi

	# combine variants
	start_time=$(date +%s)

	if [ ! -d ${final_vcf_dir} ]; then
		mkdir -p ${final_vcf_dir}
	fi

	project_vcfs=""
	if [ "${call_snps}" = "true" ]; then
		project_vcfs=${project_vcfs}"-I ${vqsr_dir}/${project_name}.hc.snps.vqsr_filt.vcf.gz "
	fi
	if [ "${call_indels}" = "true" ]; then
		project_vcfs=${project_vcfs}"-I ${vqsr_dir}/${project_name}.hc.indels.vqsr_filt.vcf.gz "
	fi
	${gatk} MergeVcfs \
		${project_vcfs} \
		-O ${final_vcf_dir}/${project_name}.hc.vqsr_filt.vcf.gz

	end_time=$(date +%s)
	running_time=$((${end_time}-${start_time}))
	echo -e "\nDONE combine variants.\nTime taken to execute commands is ${running_time} seconds.\n"

fi


# hard filter
if [ "${filtration_hard}" = "true" ]; then
	if [ ! -d ${hard_dir} ]; then
		mkdir -p ${hard_dir}
	fi

	# call germline SNPs
	if [ "${call_snps}" = "true" ]; then
		start_time=$(date +%s)

		${gatk} SelectVariants \
			-R ${reference_fa} \
			-select-type SNP \
			-V ${hc_dir}/${project_name}.hc.vcf.gz \
			-O ${hard_dir}/${project_name}.hc.snps.vcf.gz
		${gatk} VariantFiltration \
			-R ${reference_fa} \
			-filter "QD < 2.0" --filter-name "QD2" \
			-filter "QUAL < 30.0" --filter-name "QUAL30" \
			-filter "SOR > 3.0" --filter-name "SOR3" \
			-filter "FS > 60.0" --filter-name "FS60" \
			-filter "MQ < 40.0" --filter-name "MQ40" \
			-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
			-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
			-V ${hard_dir}/${project_name}.hc.snps.vcf.gz \
			-O ${hard_dir}/${project_name}.hc.snps.hard_filt.vcf.gz

		end_time=$(date +%s)
		running_time=$((${end_time}-${start_time}))
		echo -e "\nDONE hard filter snps.\nTime taken to execute commands is ${running_time} seconds.\n"
	fi

	# call germline Indels
	if [ "${call_indels}" = "true" ]; then
		start_time=$(date +%s)

		${gatk} SelectVariants \
			-R ${reference_fa} \
			-select-type INDEL \
			-V ${hc_dir}/${project_name}.hc.vcf.gz \
			-O ${hard_dir}/${project_name}.hc.indels.vcf.gz
		${gatk} VariantFiltration \
			-R ${reference_fa} \
			-filter "QD < 2.0" --filter-name "QD2" \
			-filter "QUAL < 30.0" --filter-name "QUAL30" \
			-filter "SOR > 3.0" --filter-name "SOR3" \
			-filter "FS > 200.0" --filter-name "FS200" \
			-filter "MQ < 40.0" --filter-name "MQ40" \
			-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
			-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
			-V ${hard_dir}/${project_name}.hc.indels.vcf.gz \
			-O ${hard_dir}/${project_name}.hc.indels.hard_filt.vcf.gz

		end_time=$(date +%s)
		running_time=$((${end_time}-${start_time}))
		echo -e "\nDONE hard filter indels.\nTime taken to execute commands is ${running_time} seconds.\n"
	fi

	# combine variants
	start_time=$(date +%s)

	if [ ! -d ${final_vcf_dir} ]; then
		mkdir -p ${final_vcf_dir}
	fi

	project_vcfs=""
	if [ "${call_snps}" = "true" ]; then
		project_vcfs=${project_vcfs}"-I ${hard_dir}/${project_name}.hc.snps.hard_filt.vcf.gz "
	fi
	if [ "${call_indels}" = "true" ]; then
		project_vcfs=${project_vcfs}"-I ${hard_dir}/${project_name}.hc.indels.hard_filt.vcf.gz "
	fi
	${gatk} MergeVcfs \
		${project_vcfs} \
		-O ${final_vcf_dir}/${project_name}.hc.hard_filt.vcf.gz

	end_time=$(date +%s)
	running_time=$((${end_time}-${start_time}))
	echo -e "\nDONE combine variants.\nTime taken to execute commands is ${running_time} seconds.\n"

fi


# annotate variants
if [ "${annotation}" = "funcotator" ]; then
	if [ ! -d ${anno_dir} ]; then
		mkdir -p ${anno_dir}
	fi

	start_time=$(date +%s)

	input_vcf=""
	if [ "${filtration_vqsr}" = "true" ]; then
		input_vcf=${final_vcf_dir}/${project_name}.hc.vqsr_filt.vcf.gz
	fi
	if [ "${filtration_hard}" = "true" ]; then
		input_vcf=${final_vcf_dir}/${project_name}.hc.hard_filt.vcf.gz
	fi
	output=`basename $input_vcf`
	output=${output%%.vcf.gz}

	${gatk} Funcotator \
		-R ${reference_fa} \
		-V ${input_vcf} \
		-O ${anno_dir}/${output}.anno_funcotator.vcf \
		--output-file-format VCF \
		--data-sources-path ${annotation_funcotator} \
		--ref-version ${reference_version}

	end_time=$(date +%s)
	running_time=$((${end_time}-${start_time}))
	echo -e "\nDONE annotation with Funcotator.\nTime taken to execute commands is ${running_time} seconds.\n"

fi

# annotate variants
if [ "${annotation}" = "funcotator" ]; then
	if [ ! -d ${vqsr_dir} ]; then
		mkdir -p ${anno_dir}
	fi

	start_time=$(date +%s)

	input_vcf=""
	if [ "${filtration_vqsr}" = "true" ]; then
		input_vcf=${final_vcf_dir}/${project_name}.hc.vqsr_filt.vcf.gz
	fi
	if [ "${filtration_hard}" = "true" ]; then
		input_vcf=${final_vcf_dir}/${project_name}.hc.hard_filt.vcf.gz
	fi
	output=`basename input_vcf`
	output=${output%%.vcf.gz}

	VEP --fasta $reference/Homo_sapiens_assembly38.fasta \
		--vcf --merged --fork ${threads} --hgvs --force_overwrite --everything \
		--offline --dir_cache ${annotation_vep} \
		-i ${input_vcf} \
		-o ${anno_dir}/${output}.anno_vep.vcf

	end_time=$(date +%s)
	running_time=$((${end_time}-${start_time}))
	echo -e "\nDONE annotation with VEP.\nTime taken to execute commands is ${running_time} seconds.\n"

fi

# annotate variants
if [ "${annotation}" = "vep" ]; then
	if [ ! -d ${anno_dir} ]; then
		mkdir -p ${anno_dir}
	fi

	start_time=$(date +%s)

	input_vcf=""
	if [ "${filtration_vqsr}" = "true" ]; then
		input_vcf=${final_vcf_dir}/${project_name}.hc.vqsr_filt.vcf.gz
	fi
	if [ "${filtration_hard}" = "true" ]; then
		input_vcf=${final_vcf_dir}/${project_name}.hc.hard_filt.vcf.gz
	fi
	output=`basename $input_vcf`
	output=${output%%.vcf.gz}

	${vep} --fasta ${reference_fa} \
		--vcf --fork ${threads} --hgvs --force_overwrite --everything \
		--offline --dir_cache ${annotation_vep} \
		-i ${input_vcf} \
		-o ${anno_dir}/${output}.anno_vep.vcf

	end_time=$(date +%s)
	running_time=$((${end_time}-${start_time}))
	echo -e "\nDONE annotation with VEP.\nTime taken to execute commands is ${running_time} seconds.\n"

fi

