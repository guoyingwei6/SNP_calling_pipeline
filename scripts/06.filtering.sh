# 06.filtering.sh

export chr=${1}

# HardFilter
bcftools view \
-v snps \
-e 'QD < 2.0 || QUAL < 30.0 || SOR > 3.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' \
-m2 -M2 \
05.combined_vcf/chr${chr}.vcf.gz \
-Oz \
-o 06.filtered_vcf/chr${chr}.filter.vcf.gz

# think over the threshold
bcftools view -i 'F_MISSING < 0.1 & MAF > 0.05' 06.filtered_vcf/chr${chr}.filter.vcf.gz -Oz -o 06.filtered_vcf/chr${chr}.filtered.vcf.gz

# index
bcftools index 06.filtered_vcf/chr${chr}.filtered.vcf.gz