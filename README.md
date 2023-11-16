# SNP_calling_pipeline

## 00.Preparing

### 00.Notice
- **保留所有的脚本、提交任务的命令、每个步骤的日志文件。**
- **采用下文指定的流程、软件版本和参数。**
- **确保软件所在目录已存在于PTAH变量中，或自行添加绝对路径。**

### 01.Software Version
- bwa==0.7.17
- samtools==1.18
- picard==3.1.1
- gatk==4.4.0.0
- fastp==0.23.4

### 02.Creating Directories

```shell
mkdir log 00.reference 01.fastq 02.mapping 03.calling
```

### 03.Reference Genome
> **请直接下载百度网盘链接中提供的参考基因组（UCD1.2加一条Y染色体），保存至目录`00.reference`中，并运行`md5sum -c md5.txt`查看文件是否下载完整。**  
>
> 参考基因组网盘链接: https://pan.baidu.com/s/1mvMld7s4PHXNZOtgJvS7qA?pwd=jazm 
> 

通过以下步骤为参考基因组建立索引：

1. Generate the BWA index using BWA. This will create the “.fasta.amb”, “.fasta.ann”, “.fasta.bwt”, “.fasta.pac” and “.fasta.sa” files.
```shell
bwa index GCF_002742125.1_Oar_rambouillet_v1.0_genomic.rename.AddY.fna
```
 
2. Generate the FASTA file index using samtools. This will create the “.fasta.fai” file.
```shell
samtools faidx GCF_002742125.1_Oar_rambouillet_v1.0_genomic.rename.AddY.fna
```
 
3. Generate the sequence dictionary using Picard. This will create the “.dict” file.
```shell
java -jar picard_2.18.17.jar CreateSequenceDictionary R=GCF_002742125.1_Oar_rambouillet_v1.0_genomic.rename.AddY.fna O=GCF_002742125.1_Oar_rambouillet_v1.0_genomic.rename.AddY.dict
```


> **NOTE：开始之前，检查所有的样本的fastq文件，是否都是1对，如果某些样本存在多个pane导致下机数据不止1对fastq文件，直接把该样本所有的fastq1合并在一起，所有的fastq2合并在一起，再进行接下来的操作。** 



## 01.Quality Control

采用`fastp`进行质控，默认参数即可，保留`report`：
下面代码框内容保存为01.cleaning.sh:

```shell
${sm}=$1
fastp -i ${sm}_1.fq.gz -I ${sm}_2.fq.gz -o ${sm}_1.clean.fq.gz -O ${sm}_2.clean.fq.gz -h ${sm}.html -j ${sm}.json -w 8
```

- 无论raw.fq.gz前缀是什么格式，clean.fq.gz统一命名为sample_1.clean.fq.gz及sample_1.clean.fq.gz
- -w 线程数自行调节

通过以下命令将上面脚本提交任务至计算节点（以slurm为例），每个样本一个任务：

```
# 所有样本名保存至文件lst中，每行1个
for i in `cat lst`; do sbatch -J ${sm} -p normal -N 1 -c 32 -o log/${sm}_cleaning.o -e log/${sm}_cleaning.e 01.cleaning.sh ${sm}; done
```

## 02.Mapping

开始之前，需要准备样本list文件，三列分别是：样本名，fq1的绝对路径，fq2的绝对路径，以tab分割。
**不要所有样本的bam放在同一个文件夹下面，为每个样本创建目录：**

```shell
cd 02.mapping
for i in `cut -f1 sample_lst`; do mkdir ${i}; done
cd ..
```

下面代码框内容保存为02.mapping.sh:

```shell
# define variable
export sm=${1} #sample name，作为该脚本外部的第1个参数
export fq1=${2} #fq1的绝对路径，作为该脚本外部的第2个参数
export fq2=${3} #fq2的绝对路径，作为该脚本外部的第3个参数
export thread=8 #线程数，根据集群情况自行修改
export ref=Oar4.0_add_CPY.fa  #参考基因组,记得改成绝对路径
 
# mapping and sorting
bwa mem -t ${thread} -R '@RG\tID:'${sm}'\tLB:'${sm}'\tPL:ILLUMINA\tSM:' ${sm} ${ref} ${fq1} ${fq2} | samtools sort -o 02.mapping/${sm}/${sm}.sorted.bam - -@ ${thread}

# dedup
gatk MarkDuplicates --spark-runner LOCAL -I 02.mapping/${sm}/${sm}.sorted.bam -O 02.mapping/${sm}/${sm}.dedup.bam -M 02.mapping/${sm}/${sm}.dup.metrics --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --REMOVE_DUPLICATES true

#stat
samtools flagstat 02.mapping/${sm}/${sm}.dedup.bam > 02.mapping/${sm}/${sm}.dedup.flagstat

samtools coverage 02.mapping/${sm}/${sm}.dedup.bam > 02.mapping/${sm}/${sm}.dedup.coverage

```
通过以下命令将上面脚本提交任务至计算节点（以slurm为例），每个样本一个任务：

```
cat sample_lst | while read sm fq1 fq2; do sbatch -J ${sm} -p normal -N 1 -c 32 -o log/${sm}_mapping.o -e log/${sm}_mapping.e 02.mapping.sh ${sm} ${fq1} ${fq2}; done
```



## 03.Calling

