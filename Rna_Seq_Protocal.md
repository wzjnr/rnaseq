# protocal

## 1.下载sra

## 2.sra转换成fastq

## 3.QC质控

## 4.建立对照基因组index

## 5.将fq比对到index上

## 6.将sam转换成bam并建立index

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

