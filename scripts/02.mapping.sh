# 02.mapping.sh
# define variable
export sm=${1} #sample name，作为该脚本外部的第1个参数
export fq1=${2} #fq1的绝对路径，作为该脚本外部的第2个参数
export fq2=${3} #fq2的绝对路径，作为该脚本外部的第3个参数
export thread=8 #线程数，根据集群情况自行修改
export ref=/path/00.reference/GCF_002263795.1_ARS-UCD1.2_genomic.addY.rename.fna  #参考基因组,记得改成绝对路径
 
# mapping and sorting
bwa mem -t ${thread} -R '@RG\tID:'${sm}'\tLB:'${sm}'\tPL:ILLUMINA\tSM:' ${sm} ${ref} ${fq1} ${fq2} | samtools sort -o 02.bam/${sm}/${sm}.sorted.bam - -@ ${thread}

# dedup
gatk MarkDuplicates --spark-runner LOCAL -I 02.bam/${sm}/${sm}.sorted.bam -O 02.bam/${sm}/${sm}.dedup.bam -M 02.bam/${sm}/${sm}.dup.metrics --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --REMOVE_DUPLICATES true

#stat
samtools flagstat 02.bam/${sm}/${sm}.dedup.bam > 02.bam/${sm}/${sm}.dedup.flagstat

samtools coverage 02.bam/${sm}/${sm}.dedup.bam > 02.bam/${sm}/${sm}.dedup.coverage