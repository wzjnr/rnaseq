# RNA-Seq expression analysis code and data

This repository serves as a collection of publicly available RNA-Seq datasets (moslty in maize), provides codes / pipeline to systematically analyze them as well as QC and summarized data output. The focus is on large-scale (multiple tissues / developmental stages, multiple inbred / hybrid lines) Illumina RNA-Seq experiments, but also inlcude experiments done with other sequencing technologies (3' RNA-Seq, etc.).

Raw sequencing reads were downloaded from [NCBI Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra), 

trimmed using [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) / [fastp](https://github.com/OpenGene/fastp) 

and mapped to the [maize B73 AGP_v4 genome](http://plants.ensembl.org/Zea_mays/Info/Index) using [Hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml) / [STAR](https://github.com/alexdobin/STAR).  

Uniquely mapped reads were assigned to and counted for the 46,117 reference gene models ([Ensembl Plants v37](http://plants.ensembl.org/Zea_mays/Info/Index)) using [featureCounts](http://bioinf.wehi.edu.au/featureCounts/).  

Raw read counts were then normalized using the [TMM normalization approach](https://bioconductor.org/packages/release/bioc/html/edgeR.html) to give CPMs (Counts Per Million reads)

and then further normalized by gene CDS lengths to give FPKM (Fragments Per Kilobase of exon per Million reads) values.  

Hierarchical clustering and principal component analysis were used to explore sample clustering pattern.

See [this table](/data/studies.tsv) for a list of collected datasets.

[Check here](output.md) for a walk through of output files.