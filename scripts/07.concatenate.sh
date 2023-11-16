# vcf lst
find 06.filtered_vcf/ -name chr*.filtered.vcf.gz | sort -V > 07.concatenate/vcf_lst

# concat
bcftools concat -f 07.concatenate/vcf_lst -Oz -o 07.concatenate/all_samples.all_chr.vcf.gz

# index 
bcftools index 07.concatenate/all_samples.all_chr.vcf.gz