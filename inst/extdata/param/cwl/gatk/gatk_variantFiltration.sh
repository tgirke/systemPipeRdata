#!/bin/bash -l

java=$1
input_vcf=$2
output_vcf=$3
ref=$4

gatk VariantFiltration \
    --java-options "${java}" \
    -V ${input_vcf} \
    -O ${output_vcf} \
    -R ${ref} \
  	-filter-name "QD_filter" -filter "QD < 2.0" \
  	-filter-name "FS_filter" -filter "FS > 60.0" \
  	-filter-name "MQ_filter" -filter "MQ < 40.0" \
  	-filter-name "SOR_filter" -filter "SOR > 4.0" \
  	-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
  	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
  	-genotype-filter-expression "DP < 10" \
  	-genotype-filter-name "DP_filter" \
  	-genotype-filter-expression "GQ < 10" \
  	-genotype-filter-name "GQ_filter"
