# 0405.Combine_Genotype_GVCFs.sh
export chr=${1}
export ref=/path/00.reference/GCF_002263795.1_ARS-UCD1.2_genomic.addY.rename.fna

# CombineGVCFs 
find 03.gvcf/ -name "*.chr${chr}.gvcf.gz" > 04.combined_gvcf/chr${chr}.gvcf.lst
gatk CombineGVCFs -R ${ref} -V 04.combined_gvcf/chr${chr}.gvcf.lst -O 04.combined_gvcf/chr${chr}.gvcf.gz
 
# GenotypeGVCFs
gatk GenotypeGVCFs -R ${ref} -V 04.combined_gvcf/chr${chr}.gvcf.gz -O 05.combined_vcf/chr${chr}.vcf.gz