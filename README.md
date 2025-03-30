# **The pipelines and codes of population genomics and transcriptomics.**


- [**The pipelines and codes of population genomics and transcriptomics.**](#the-pipelines-and-codes-of-population-genomics-and-transcriptomics)
- [SNP检查和质量控制标准操作](#snp检查和质量控制标准操作)
  - [1. 重测序的比对和`.bam`文件处理](#1-重测序的比对和bam文件处理)
  - [2. 利用GATK4进行SNP Calling](#2-利用gatk4进行snp-calling)
  - [3. Genotype quality control](#3-genotype-quality-control)
- [其它有用的命令](#其它有用的命令)
  - [PLINK v1.9软件的使用](#plink-v19软件的使用)
  - [利用`blupADC`快速实现基因型数据的格式转换](#利用blupadc快速实现基因型数据的格式转换)


# SNP检查和质量控制标准操作
## 1. 重测序的比对和`.bam`文件处理 
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
3.2 利用[filter_by_sample_depth.pl](https://github.com/Mou-Yin/Population_genomics_and_transcriptomics/edit/main/filter_by_sample_depth.pl)根据每个样本的平均深度过滤 VCF 文件中的 SNP 和 Indel 位点。具体来说，它检查每个样本在每个位点的深度（DP 值），如果深度异常（低于平均深度的三分之一或高于三倍），则将该样本的基因型设置为缺失（./.）
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


计算结果 .hwe 文件                               
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






# 其它有用的命令
## PLINK v1.9软件的使用

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









## 利用`blupADC`快速实现基因型数据的格式转换

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


