
- [The standard procedures SNP Calling and QC (Done)](#the-standard-procedures-snp-calling-and-qc-done)
  - [1. Read mapping and dealing `.bam` file](#1-read-mapping-and-dealing-bam-file)
  - [2. 利用GATK4进行SNP Calling](#2-利用gatk4进行snp-calling)
  - [3. Genotype quality control](#3-genotype-quality-control)
- [GWAS and QTL analysis (coming soon)](#gwas-and-qtl-analysis-coming-soon)
- [Analysis of differentially expressed genes at population level (Done)](#analysis-of-differentially-expressed-genes-at-population-level-done)
- [Weighted Gene Correlation Network Analysis (WGCNA) (Coming soon)](#weighted-gene-correlation-network-analysis-wgcna-coming-soon)
- [eQTL mapping and eQTG identification (Coming soon)](#eqtl-mapping-and-eqtg-identification-coming-soon)
  - [1. Normalization of gene expression matrices](#1-normalization-of-gene-expression-matrices)
  - [2. Using peer to remove batch effect](#2-using-peer-to-remove-batch-effect)
  - [3. eQTL mapping](#3-eqtl-mapping)
  - [4. Calculate the interval of QTL](#4-calculate-the-interval-of-qtl)
  - [5. Identification of eGTG based on WGCNA](#5-identification-of-egtg-based-on-wgcna)
- [TWAS using FUSION (Coming soon)](#twas-using-fusion-coming-soon)
  - [Prepare input files](#prepare-input-files)
  - [Compute predictive models for each gene](#compute-predictive-models-for-each-gene)
  - [Perform the expression imputation](#perform-the-expression-imputation)
- [Fine mapping of genes using cTWAS (Coming soon)](#fine-mapping-of-genes-using-ctwas-coming-soon)
- [其它有用的命令 (In Progress)](#其它有用的命令-in-progress)
  - [PLINK v1.9软件的使用 (In Progress)](#plink-v19软件的使用-in-progress)
  - [利用blupADC快速实现基因型数据的格式转换](#利用blupadc快速实现基因型数据的格式转换)


# The standard procedures SNP Calling and QC (Done)
## 1. Read mapping and dealing `.bam` file
1.1 利用`samtools`和`BWA`软件构建参考基因组索引
```
samtools ref.fa
bwa index ref.fa
```
1.2 利用`BWA`的`mem`算法将clean reads比对到参考基因组
```
bwa mem -t 6 -R "@RG\tID:C1388\tPL:illumina\tPU:illumina\tLB:C1388\tSM:C1388" NtaSR1_genome_1.0.fa 01.RawData/C1388_1.fq.gz 01.RawData/C1388_2.fq.gz > 00.mapping/C1388.bam
```
1.3 将`.bam`文件按染色体坐标进行排序
```
/dt1/share/software/samtools-1.18/samtools sort -@ 2 -T tmp/C1388 -o 00.mapping/C1388.sort.bam 00.mapping/C1388.bam
```
1.4 利用`picard`移除除PCR duplicated reads
```
java -Xmx10g -jar /dt1/yinm/software/picard/picard.jar MarkDuplicates
		INPUT=00.mapping/C1075.sort.bam \
		OUTPUT=00.mapping/C1075.rmdup.bam \
		METRICS_FILE=00.mapping/C1075.dup.txt \
		REMOVE_DUPLICATES=true
```
1.5 构建`.bam`文件索引
```
/dt1/share/software/samtools-1.18/samtools index 00.mapping/C1075.rmdup.bam
```
1.6 利用`mosdepth`统计样本的有效测序深度
```
mosdepth -t 4 00.mapping/C1 00.mapping/C1.rmdup.bam
```

  
## 2. 利用GATK4进行SNP Calling
2.1 利用`HaplotypeCaller`将`.bam`文件转换成`.gvcf`
```
/dt1/share/software/gatk-4.4.0.0/gatk --java-options "-Xmx20g" HaplotypeCaller \
		-R /dt2/yinm/project/tobacco/00.ref_genome/NtaSR1_genome_1.0.fa \
		-I 00.mapping/C1075.rmdup.bam \
		-O 01.gvcf/C1075.g.vcf.gz \
		--emit-ref-confidence GVCF
```
2.2 利用`CombineGVCFs`合并`GVCF`（按染色体拆分并行）
```
/dt1/share/software/gatk-4.4.0.0/gatk CombineGVCFs \
		-R /dt2/yinm/project/tobacco/00.ref_genome/NtaSR1_genome_1.0.fa \
		-V 03.gvcf.list \
		-L Chr01  \
		-O 03.combineGVCF_split_Chr/03.Chr01_merge_220samples.gvcf.gz`
```
2.3 利用`GenotypeGVCFs`联合分型
```
/dt1/share/software/gatk-4.4.0.0/gatk GenotypeGVCFs \
		-R /dt2/yinm/project/tobacco/00.ref_genome/NtaSR1_genome_1.0.fa \
		-V 03.combineGVCF_split_Chr/03.Chr01_merge_220samples.gvcf.gz \
		-O 04.genotype/Chr01_raw_220samples.vcf.gz
```
2.4 拆分SNP和INDEL
```bash
#SNP:
/dt1/share/software/gatk-4.4.0.0/gatk  SelectVariants \
	-V 04.genotype/Chr01_raw_220samples.vcf.gz \
	--select-type-to-include SNP \
	-O 05.filter/Chr01_raw_220samples.SNP.vcf.gz

#INDEL:
/dt1/share/software/gatk-4.4.0.0/gatk  SelectVariants \
		-V 04.genotype/Chr01_raw_220samples.vcf.gz \
		--select-type-to-include INDEL \
		-O 05.filter/Chr01_raw_220samples.INDEL.vcf.gz
```

  
## 3. Genotype quality control
3.1 利用`VariantFiltration`对SNP和INDEL分别进行硬过滤
```bash
#SNP:
/dt1/share/software/gatk-4.4.0.0/gatk VariantFiltration \
		-V 05.filter/Chr01_raw_220samples.SNP.vcf.gz \
		-filter "QD < 2.0" --filter-name "QD2" \
		-filter "QUAL < 30.0" --filter-name "QUAL30" \
		-filter "SOR > 3.0" --filter-name "SOR3" \
		-filter "FS > 60.0" --filter-name "FS60" \
		-filter "MQ < 40.0" --filter-name "MQ40" \
		-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
		-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
		-O 05.filter/Chr01_220samples.SNP.hardfilter.vcf.gz

#INDEL
/dt1/share/software/gatk-4.4.0.0/gatk VariantFiltration \
		-V 05.filter/Chr01_raw_220samples.INDEL.vcf.gz \
		-filter "QD < 2.0" --filter-name "QD2" \
		-filter "QUAL < 30.0" --filter-name "QUAL30" \
		-filter "FS > 200.0" --filter-name "FS200" \
		-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
		-O 05.filter/Chr01_220samples.INDEL.hardfilter.vcf.gz
```
3.2 利用[filter_by_sample_depth.pl](https://github.com/Mou-Yin/Population_genomics_and_transcriptomics/blob/main/scripts/filter_by_sample_depth.pl)根据每个样本的平均深度过滤 VCF 文件中的 SNP 和 Indel 位点。具体来说，它检查每个样本在每个位点的深度（DP 值），如果深度异常（低于平均深度的三分之一或高于三倍），则将该样本的基因型设置为缺失（./.）
```perl
perl filter_by_sample_depth.pl 00.summary_sample_depth.txt 05.filter/Chr01_220samples.INDEL.hardfilter.vcf > 05.filter/Chr01_220samples.INDEL.DPfilter.vcf
perl filter_by_sample_depth.pl 00.summary_sample_depth.txt 05.filter/Chr01_220samples.SNP.hardfilter.vcf > 05.filter/Chr01_220samples.SNP.DPfilter.vcf
# 注意，该脚本同时过滤了多等位位点和没有“PASS”标签的位点

head -n 4 00.summary_sample_depth.txt                      
# C1075	11.35
# C1102	11.46
# C1103	12.81
# C1104	12.21
```

3.3 合并每条染色体的vcf
```
/dt1/share/software/gatk-4.4.0.0/gatk  GatherVcfs -I 05.filter/Chr01_220samples.INDEL.final.recode.vcf -I 05.filter/Chr01_220samples.SNP.final.recode.vcf ......(etc.) -O tobacco220.vcf
```

3.4 计算样本的杂合度。对于自交系群体，可以考虑过滤掉高杂合样本。相反，对于杂交系群体，可以考虑过滤掉高纯合的样本。
```bash
# 将vcf转化成plink的格式
plink --vcf tobacco220.vcf --make-bed --chr-set 24 --recode --out tobacco220
# 计算杂合率
plink --bfile tobacco220 --het --out tobacco220

# 计算结果 .het 文件，F列的值越大，表示样本纯和度越高，负值代表样本是高杂合的。
#   FID     IID       O(HOM)       E(HOM)        N(NM)            F
# C1075   C1075      2966558    2.137e+06      2969718       0.9962
# C1102   C1102      2973284    2.142e+06      2975577       0.9973
```


3.5 计算位点杂合率。对于自交系群体，过滤掉杂合基因型比例大于5%的位点
```bash
plink --vcf your_vcf_file --make-bed --chr-set 24 --recode --out tobacco220
plink --bfile vcf_file_prefix --chr-set 24 --hardy --out output_prefix


# 计算结果 .hwe 文件                               
# CHR	SNP	TEST	A1	A2	GENO	O(HET)	E(HET)	P
# 1	01_391	ALL(NP)	T	C	0/1/213	0.004673	0.004662	1
# 1	01_743	ALL(NP)	T	C	6/0/207	0	0.05475	1.894e-12
# 第七列是位点杂合率
```

3.6 过滤高杂合的样本，以及高缺失率和低频的位点
```bash
# 过滤高杂合度的样本，以及过滤次等位基因频率小于0.05、缺失率小于0.2的位点
# INDEL
vcftools --vcf your_vcf_file --keep Tobacco215.list --max-missing 0.8 --maf 0.05 --recode --recode-INFO-all --out output_prefix
# SNP
vcftools --vcf your_vcf_file --keep Tobacco215.list --max-missing 0.8 --maf 0.05 --recode --recode-INFO-all --out output_prefix

head -n 3 Tobacco215.list                     
#C1587
#C1154
#C727
```
> [!TIP]
> 在上述中，先用GATK的`VariantFiltration`程序对每个位点进行了硬过滤；然后基于每个样本的平均有效深度对每个样本的每个位点的深度进行检查，将测序深度异常的位点换成`./.`；之后我计算了样本的杂合度，发现有一个样本高杂合，因为我是自交系群体，因此过滤掉了这个样本；然后我针对每个位点，计算了群体中杂合基因型的比例，过滤比例大于5%的位点；同时过滤了位点的缺失率大于20%的位点，以及MAF小于0.05的位点。
>   
> 然而，在实际分析中，需要根据自己的分析目的进行样本和位点的过滤。不同的分析分析应该执行不同的过滤，特别是MAF。例如，计算核苷酸多样性应该保留稀有或低频的位点；相反，在推断群体结构、亲缘关系等，需要过滤稀有或低频的位点。
>   
> 此外，VCF的过滤顺序也很重要，通常应该先进行样本过滤，后执行位点过滤。所以，一般先根据样本缺失率、杂合度等进行计算和过滤，再根据位点缺失率、基因型频率、杂合率等进行过滤。


# GWAS and QTL analysis (coming soon)


# Analysis of differentially expressed genes at population level (Done)
群体水平的差异表达基因分析应该使用Wilcoxon检验


1. 读取所有样本的表达量矩阵(00.sample860_filter_expr_mat.txt)，该文件包含215份样本的4个生长发育时期的基因表达数据。
```
# 读取表达量文件
data <- fread("/dt2/yinm/project/tobacco/01.DEGs/00.sample860_filter_expr_mat.txt",header= T)
# 提取苗期的表达量数据
data_stage =  filter(data,Stage=="Seedling")
data_stage = as.data.frame(data_stage)
rownames(data_stage) <- data_stage$Ind1
# 提取表达量矩阵并转置
data1 = data_stage[,c(8:ncol(data_stage))]
data2 = t(data1) # 每一行一个基因，每一列一个样本
```

2. 读取样本分组信息
```
info <- read.table("/dt2/yinm/project/tobacco/00.sample_info/tobacco220_info.txt",header = T)
# tobacco220_info.txt
# Type	ID	Accessions	Name	Era	ID
# Fluecured	C1	KT001	Hongda	1989	C0001
# Fluecured	C2	KT002	K326	1991	C0002
# Fluecured	C3	KT003	K346	1997	C0003

FT <- info[info$Type == "Fluecured", 2]
CT <- info[info$Type == "Cigar", 2]
BT <- info[info$Type == "Burley", 2]
ST <- info[info$Type == "Suncured", 2]
OT <- info[info$Type == "Oriental", 2]

# 仅保留在表达矩阵中的样本，确保用于分析的样本确实存在于基因表达矩阵
FT <- intersect(FT, colnames(data2))
CT <- intersect(CT, colnames(data2))
BT <- intersect(BT, colnames(data2))
ST <- intersect(ST, colnames(data2))
OT <- intersect(OT, colnames(data2))
```

3. 提取这两个群体的基因表达矩阵，并过滤
```
# 设置要比较的两个群体
pop1 = FT
pop2 = CT

two_pop_expr = data2[,c(pop1,pop2)]
nrow(two_pop_expr)

# 过滤掉在少于10个样本中表达的基因，即要求每个基因至少在10个样本中的表达，TPM > 0.5则认为表达。
filtered_expr <- two_pop_expr[rowSums(two_pop_expr > 0.5) >= 10, ] 
nrow(filtered_expr)
```

4. 对表达量矩阵进行标准化。在比较两个群体的时候，应该提取两个群体的基因表达矩阵，然后进行标准化。
```
norm_expr = normalize.quantiles.robust(x = as.matrix(filtered_expr),use.log2 = FALSE,keep.names = T)
```

5. 定义函数，用于计算每个基因的平均表达量，log2(fold change)，和Wilcoxon检验的p值
```
calculate_metrics <- function(gene_expression) {
  # 提取群体1的基因表达
  exp_group1 <- as.numeric(gene_expression[pop1])
  exp_group2 <- as.numeric(gene_expression[pop2])
  # 提取群体2的基因表达
  mean_group1 <- mean(exp_group1, na.rm = TRUE)
  mean_group2 <- mean(exp_group2, na.rm = TRUE)
  # 计算 fold change
  log2_fold_change <- log2(mean_group1 + 0.1) - log2(mean_group2 + 0.1)
  ## 计算p值
  p_value <- tryCatch({
    wilcox.test(exp_group1, exp_group2)$p.value
  }, error = function(e) NA)
  # tryCatch({ ... }, error = function(e) NA): tryCatch() 是一个错误处理机制，用于捕获代码块中可能发生的错误
  return(c(mean_group1, mean_group2, log2_fold_change, p_value))
}
```

6. 进行计算，fdr检验，和DEGs提取
```
result <- data.frame(Gene = rownames(norm_expr))
# 计算
result[, c("Mean_FT", "Mean_CT", "log2FC", "pvalue")] <- 
  t(apply(norm_expr, 1, calculate_metrics)) #apply()函数用于对矩阵或数组的行或列应用函数 1代表对行进行操作

result$Mean_FT <- round(result$Mean_FT,digits = 3) # 保留3为小数
result$Mean_CT <- round(result$Mean_CT,digits = 3)

# fdr
result$padj <- p.adjust(result$pvalue, method = "fdr",  n = length(result$pvalue))

# 提取差异表达的基因
deg = filter(result, (log2FC >= 1 | log2FC <= -1) & padj < 0.01 )

# 保存结果
write.table(deg, "/dt2/yinm/project/tobacco/01.DEGs/02.DEG_resluts/DEGs_Seedling_pop1_vs_pop2.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(result, "/dt2/yinm/project/tobacco/01.DEGs/02.DEG_resluts/AllGene_Seedling_pop1_vs_pop2.txt", sep = "\t", row.names = FALSE, quote = FALSE)
```



# Weighted Gene Correlation Network Analysis (WGCNA) (Coming soon)

# eQTL mapping and eQTG identification (Coming soon)
## 1. Normalization of gene expression matrices

## 2. Using peer to remove batch effect

## 3. eQTL mapping

## 4. Calculate the interval of QTL

## 5. Identification of eGTG based on WGCNA



# TWAS using FUSION (Coming soon)

## Prepare input files

## Compute predictive models for each gene

## Perform the expression imputation



# Fine mapping of genes using cTWAS (Coming soon)





# 其它有用的命令 (In Progress)
## PLINK v1.9软件的使用 (In Progress)

ped/map, bed/bim/fam, tped/tfam, 和VCF格式互换
```bash
# vcf格式转plink格式
plink --vcf input_vcf_file --keep-allele-order --allow-extra-chr --recode --out out_file_prefix

# ped/map转换为tped/tfam
plink --file geno_file_prefix --recode --transpose --out out_file_prefix

# tped/tfam转换为ped/map
plink --tfile geno_file_prefix --recode --out out_file_prefix

# ped/map转换为bed/bim/fam
plink --file geno_file_prefix --make-bed --out out_file_prefix

# bed/bim/fam转换为ped/map
plink --bfile geno_file_prefix --recode --out out_file_prefix

# tped/tfam转换为bed/bim/fam
plink --tfile geno_file_prefix --make-bed --out out_file_prefix

# bed/bim/fam转换为tped/tfam
plink --bfile geno_file_prefix --recode --transpose --out out_file_prefix

## bed/bim/fam 转为 vcf
plink --bfile geno_file_prefix --export vcf --out out_file_prefix
```

使用`--keep-allele-order`参数保持REF和ALT的顺序
使用`--allow-extra-chr`允许染色体的ID是非chr1或1的形式

`PLINK`是为人类研究定制的软件，使用`--chr-set`允许处理超过22条的染色体
```
plink --vcf geno.vcf --make-bed --chr-set 24  --recode --out output_prefix
```
`GCTA`默认也只能处理22条染色体，需要加--autosome-num这个参数限制常染色体的数量.

plink计算PCA
```
plink -bfile geno_file_prefix --pca 10 --allow-extra-chr --threads 20 --out out_file_prefix
```


利用`--indep-pairwise`对基因型数据进行LD prune
```
plink --file W200k.final.convert --allow-extra-chr --keep-allele-order  --indep-pairwise 1000kb 1 0.5 --out W200k.final.prune
```

计算IBS

```
plink --allow-extra-chr --file W200k.final.convert.vcf --distance ibs 1-ibs allele-ct flat-missing square --out IBS
```


随机挑选100个SNP
```
plink --vcf W200k.final.1B.vcf --make-bed  --thin-count 1000 --allow-extra-chr  --recode --keep-allele-order --out W200k.final.1B.randomly1000SNP.vcf
```

转换成VCF
```
plink --bfile W200k.final.1B.randomly1000SNP.vcf --allow-extra-chr  --export vcf --out W200k.final.1B.randomly1000SNP
```

计算位点的杂合度
```
plink --bfile $in --hardy  --out $out
```

计算样本的杂合度
```
plink --file $in --het --out $out
```








## 利用blupADC快速实现基因型数据的格式转换

```R
format_result=geno_format(
    input_data_hmp=example_data_hmp,  # provided data variable
    output_data_type=c("Plink","BLUPF90","Numeric"),# output data format
	output_data_path=getwd(),   #output data path      
    output_data_name="blupADC", #output data name    
    return_result = TRUE,       #save result in R environment
    cpu_cores=1                 # number of cpu 
    )
```
blupADC还可以实现基因型数据的质控与填充、品种分析及基因型数据重复性检测、品种分析及基因型数据重复性检测、系谱可视化、亲缘关系矩阵的构建(A,G, H)、利用DMU或BLUPF90软件进行遗传评估[https://github.com/TXiang-lab/blupADC](https://github.com/TXiang-lab/blupADC)


