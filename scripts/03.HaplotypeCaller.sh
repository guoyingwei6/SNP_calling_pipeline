# 03.HaplotypeCaller.sh
# HaplotypeCaller
export sm=${1}
export bam=${2}
export chr=${3}
export ref=/path/00.reference/GCF_002263795.1_ARS-UCD1.2_genomic.addY.rename.fna
 
gatk HaplotypeCaller -R ${ref} -I ${bam} -L ${chr} -ERC GVCF -O 03.gvcf/${sm}/${sm}.chr${chr}.gvcf.gz