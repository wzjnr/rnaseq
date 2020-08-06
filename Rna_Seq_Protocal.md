# Protocal

## 1.下载sra

## 2.sra转换成fastq

## 3.QC质控

## 4.建立对照基因组index

## 5.将fq比对到index上

## 6.将sam转换成bam并建立index

```bash
#view语句

-b output BAM  
# 该参数设置输出 BAM 格式，默认下输出是 SAM 格式文件
-h print header for the SAM output 
# 默认下输出的 sam 格式文件不带 header，该参数设定输出sam文件时带 header 信息 
-H print SAM header only (no alignments)  
# 仅仅输出文件的头文件
-S input is SAM  
# 默认下输入是 BAM 文件，若是输入是 SAM 文件，则最好加该参数，否则有时候会报错。 
-u uncompressed BAM output (force -b)  
# 该参数的使用需要有-b参数，能节约时间，但是需要更多磁盘空间。 
-c print only the count of matching records
# 仅输出匹配的统计记录 
-L FILE  only include reads overlapping this BED FILE [null] 
#  仅包括和bed文件存在overlap的reads
-o FILE  output file name [stdout] 
# 输出文件的名称
-F INT  only include reads with none of the FLAGS in INT present [0] 
# 过滤flag，仅输出指定FLAG值的序列
-q INT   only include reads with mapping quality >= INT [0]    
# 比对的最低质量值，一般认为20就为unique比对了，可以结合上述-bF参数使用使用提取特定的比对结果
-@ Number of additional threads to use [0]
# 指使用的线程数

# 将sam文件转换成bam文件
samtools view -bS abc.sam > abc.bam
# BAM转换为SAM
samtools view -h -o out.sam out.bam
# 提取比对到参考序列上的比对结果
samtools view -bF 4 abc.bam > abc.F.bam
# 提取paired reads中两条reads都比对到参考序列上的比对结果，只需要把两个4+8的值12作为过滤参数即可
samtools view -bF 12 abc.bam > abc.F12.bam
# 提取没有比对到参考序列上的比对结果
samtools view -bf 4 abc.bam > abc.f.bam
# 提取bam文件中比对到caffold1上的比对结果，并保存到sam文件格式
samtools view abc.bam scaffold1 > scaffold1.sam
# 提取scaffold1上能比对到30k到100k区域的比对结果
samtools view abc.bam scaffold1:30000-100000 $gt; scaffold1_30k-100k.sam
# 根据fasta文件，将 header 加入到 sam 或 bam 文件中
samtools view -T genome.fasta -h scaffold1.sam > scaffold1.h.sam
```

```bash
#sort语句

Usage: samtools sort [option] <in.bam> -o <out.prefix>  
-n Sort by read name
#设定排序方式按short reads的ID排序。默认下是按序列在fasta文件中的顺序（即header）和序列从左往右的位点排序。
-m INT     Set maximum memory per thread; suffix K/M/G recognized [768M]
# 设置每个线程的最大内存，单位可以是K/M/G，默认是 768M。对于处理大数据时，如果内存够用，则设置大点的值，以节约时间。
-t TAG     Sort by value of TAG. Uses position as secondary index (or read name if -n is set)
# 按照TAG值排序
-o FILE    Write final output to FILE rather than standard output 
# 输出到文件中，加文件名

#例子
#  tmp.bam 按照序列位置排序，并将结果输出到tmp.sort.bam  
samtools sort -n tmp.bam  -o tmp.sort.bam    
samtools view tmp.sort.bam 
```

```bash
#merge和cat语句
#merge将多个已经sort了的bam文件融合成一个bam文件。融合后的文件不需要则是已经sort过了的。而cat命令不需要将bam文件进行sort。

Usage: samtools merge [-nurlf] [-h inh.sam] [-b <bamlist.fofn>] <out.bam> <in1.bam> [<in2.bam> ... <inN.bam>]

Options: 
  -n         Input files are sorted by read name
# 输入文件是经过sort -n的
  -t TAG     Input files are sorted by TAG value
# 输入文件是经过sort -t的
  -r         Attach RG tag (inferred from file names)
# 添加上RG标签
  -u         Uncompressed BAM output
# 输出未压缩的bam
  -f         Overwrite the output BAM if exist
# 覆盖已经存在的bam
  -1         Compress level 1
# 1倍压缩
  -l INT     Compression level, from 0 to 9 [-1]
# 指定压缩倍数
  -R STR     Merge file in the specified region STR [all]
  -h FILE    Copy the header in FILE to <out.bam> [in1.bam]

$samtools cat
Usage: samtools cat [options] <in1.bam>  [... <inN.bam>]
       samtools cat [options] <in1.cram> [... <inN.cram>]

Options: -b FILE  list of input BAM/CRAM file names, one per line
         -h FILE  copy the header from FILE [default is 1st input file]
         -o FILE  output BAM/CRAM
```

```bash
#index语句

Usage: samtools index <in.bam> [out.index]

samtools index abc.sort.bam
```



```bash
##使用SAMtools进行操作脚本
#注意：由于未知原因，bioinfo环境操作缺少组件，需要先转换到samtools环境再操作

for i in {16..37}
do
    samtools view -Sb ~/ncbi/public/sra/darmor_sam/SRR78952${i}.sam > ~/ncbi/public/sra/darmor_sam/SRR78952${i}.bam
    samtools  sort -n ~/ncbi/public/sra/darmor_sam/SRR78952${i}.bam -o ~/ncbi/public/sra/darmor_sam/SRR78952${i}_sorted.bam
    samtools index ~/ncbi/public/sra/darmor_sam/SRR78952${i}_sorted.bam
done

##最终得到
#未排序的*.sam;
#排了序的*_sored.bam;
#以及index
```

