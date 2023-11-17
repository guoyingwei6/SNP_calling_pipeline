# SNP_calling_pipeline

## 00.Preparing

### 00.Notice
- **保留所有的脚本、提交任务的命令、每个步骤的日志文件。**
- **采用下文指定的流程、软件版本和参数，每个步骤的输出文件名后缀必须同示例保持一致。**
- **确保软件所在目录已存在于PTAH变量中，或自行添加绝对路径。**

### 01.Software Version
- bwa==0.7.17
- gatk==4.4.0.0
- samtools==1.18
- bcftools==1.18
- fastp==0.23.4

### 02.Creating Directories

```shell
mkdir log 00.reference 01.fastq 02.bam 03.gvcf 04.combined_gvcf 05.combined_vcf 06.filtered_vcf 07.concatenate
```

### 03.Reference Genome

0. 以NCBI上*Bos taurus*最新的ARS-UCD2.0版本基因组为参考基因组，并保存至目录`00.reference`中,对染色体进行重命名（常染色体，X，T，MT重命名，scafolds不用管），python脚本在scripts目录下：

```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/263/795/GCF_002263795.3_ARS-UCD2.0/GCF_002263795.3_ARS-UCD2.0_genomic.fna.gz
python3 scripts/ChangeChromosomeNameInFastaOrGff.py -i GCF_002263795.3_ARS-UCD2.0_genomic.fna.gz --fa -g -N chr_name -o GCF_002263795.3_ARS-UCD2.0_genomic.rename.fna
```

通过以下步骤为参考基因组建立索引：

1. Generate the BWA index using BWA. This will create the “.fasta.amb”, “.fasta.ann”, “.fasta.bwt”, “.fasta.pac” and “.fasta.sa” files.
```shell
bwa index GCF_002263795.1_ARS-UCD1.2_genomic.addY.rename.fna
```
 
2. Generate the FASTA file index using samtools. This will create the “.fasta.fai” file.
```shell
samtools faidx GCF_002263795.1_ARS-UCD1.2_genomic.addY.rename.fna
```
 
3. Generate the sequence dictionary using Picard. This will create the “.dict” file.
```shell
java -jar picard_2.18.17.jar CreateSequenceDictionary R=GCF_002263795.1_ARS-UCD1.2_genomic.addY.rename.fna O=GCF_002263795.1_ARS-UCD1.2_genomic.addY.rename.dict
```


> **NOTE：开始之前，检查所有的样本的fastq文件，是否都是1对，如果某些样本存在多个lane导致下机数据不止1对fastq文件，直接把该样本所有的fastq1合并在一起，所有的fastq2合并在一起，再进行接下来的操作。** 



## 01.Quality Control

采用`fastp`进行质控，默认参数即可，保留`report`：
下面代码框内容保存为`01.cleaning.sh`:

```shell
# 01.cleaning.sh
${sm}=$1
fastp -i ${sm}_1.fq.gz -I ${sm}_2.fq.gz -o ${sm}_1.clean.fq.gz -O ${sm}_2.clean.fq.gz -h ${sm}.html -j ${sm}.json -w 8
```

- 无论raw.fq.gz前缀是什么格式，clean.fq.gz统一命名为sample_1.clean.fq.gz及sample_2.clean.fq.gz
- -w 线程数自行调节

通过以下命令将上面脚本提交任务至计算节点（以slurm为例），每个样本一个任务（合理调整线程和内存）。：

```
# 所有样本名保存至文件lst中，每行1个
for i in `cat lst`; do sbatch -J ${sm} -p normal -N 1 -c 32 -o log/${sm}_cleaning.o -e log/${sm}_cleaning.e 01.cleaning.sh ${sm}; done
```

## 02.Mapping

开始之前，需要准备样本list文件，三列分别是：样本名，fq1的绝对路径，fq2的绝对路径，以tab分割。
**不要所有样本的bam放在同一个文件夹下面，要为每个样本创建单独的目录：**

```shell
cd 02.mapping
for i in `cut -f1 sample_lst`; do mkdir ${i}; done
cd ..
```

下面代码框内容保存为`02.mapping.sh`:

```shell
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

```
通过以下命令将上面脚本提交任务至计算节点（以slurm为例），每个样本一个任务（合理调整线程和内存）：

```
cat sample_lst | while read sm fq1 fq2; do sbatch -J ${sm} -p normal -N 1 -c 32 -o log/${sm}_mapping.o -e log/${sm}_mapping.e 02.mapping.sh ${sm} ${fq1} ${fq2}; done
```



## 03.Calling

### 01.HaplotypeCaller

开始之前，需要提前准备包含样本名和bam文件绝对路径的list文件，一共2列，以tab分割。
**不要所有样本的gvcf放在同一个文件夹下面，要为每个样本创建单独的目录：**

```shell
cd 03.gvcf
for i in `cut -f1 bam_lst`; do mkdir ${i}; done
cd ..
```

下面代码框内容保存为`03.HaplotypeCaller.sh`:
```shell
# 03.HaplotypeCaller.sh
# HaplotypeCaller
export sm=${1}
export bam=${2}
export chr=${3}
export ref=/path/00.reference/GCF_002263795.1_ARS-UCD1.2_genomic.addY.rename.fna
 
gatk HaplotypeCaller -R ${ref} -I ${bam} -L ${chr} -ERC GVCF -O 03.gvcf/${sm}/${sm}.chr${chr}.gvcf.gz
```

这一步要求分染色体进行calling，**需要29条常染色体加上X Y MT共32条染色体**，每个个体每条染色体分别递交任务，注意嵌套循环：
```
cat bam_lst | while read sm bam; do for chr in {1..29} X Y MT; do sbatch -J ${sm} -p normal -N 1 -c 32 -o log/${sm}_${chr}_calling.o -e log/${sm}_${chr}_calling.e 03.HaplotypeCaller.sh ${sm} ${bam} ${chr}; done; done
```


---


### 02.CombineGVCFs and GenotypeGVCFs
下面代码框内容保存为`04.05.Combine_Genotype_GVCFs.sh`:
```shell
# 04.05.Combine_Genotype_GVCFs.sh
export chr=${1}
export ref=/path/00.reference/GCF_002263795.1_ARS-UCD1.2_genomic.addY.rename.fna

# CombineGVCFs 
find 03.gvcf/ -name "*.chr${chr}.gvcf.gz" > 04.combined_gvcf/chr${chr}.gvcf.lst
gatk CombineGVCFs -R ${ref} -V 04.combined_gvcf/chr${chr}.gvcf.lst -O 04.combined_gvcf/chr${chr}.gvcf.gz
 
# GenotypeGVCFs
gatk GenotypeGVCFs -R ${ref} -V 04.combined_gvcf/chr${chr}.gvcf.gz -O 05.combined_vcf/chr${chr}.vcf.gz
```

按照上面的示例，每条染色体提交一次任务：
```
for i in {1..29} X Y MT; do 
```


### 03.filtering
下面代码框内容保存为`06.filtering.sh`:
```shell
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
```

### 04.concatenate

将全部个体的每条染色体合并在一起，生成全部样本全基因组的vcf：
下面代码框内容保存为`07.concatenate.sh`:

```shell
# vcf lst
find 06.filtered_vcf/ -name chr*.filtered.vcf.gz | sort -V > 07.concatenate/vcf_lst

# concat
bcftools concat -f 07.concatenate/vcf_lst -Oz -o 07.concatenate/all_samples.all_chr.vcf.gz

# index 
bcftools index 07.concatenate/all_samples.all_chr.vcf.gz
```

