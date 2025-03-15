# **Population transcriptomics of common tobacco**
I will continuously update the pipelines and codes of population genomics and transcriptomics, including
1. SNP calling and QC
2. PCA
3. Admixture
4. Constucting Neighbour-Joining tree
5. Selective sweeps
6. Inference of haplotye blocks
7. SNP and haploblock-based GWAS
8. QTL analysis
9. Quantification of genes or transcripts
10. Normalization of expression matrix
11. Elimination of RNA-seq batch effect
12. eQTL
13. TWAS
14. SMR
15. Other


  
## The standard procedures of SNP Calling and QC
### 1. Mapping of reads and dealing of .bam files 
1.1 利用samtools和BWA软件构建参考基因组索引
```
samtools ref.fa
bwa index ref.fa
```
1.2 利用BWA的mem算法将clean reads比对到参考基因组
```
bwa mem -t 6 -R "@RG\tID:C1388\tPL:illumina\tPU:illumina\tLB:C1388\tSM:C1388" NtaSR1_genome_1.0.fa 01.RawData/C1388_1.fq.gz 01.RawData/C1388_2.fq.gz > 00.mapping/C1388.bam
```
1.3 将bam文件按染色体坐标进行排序
```
/dt1/share/software/samtools-1.18/samtools sort -@ 2 -T tmp/C1388 -o 00.mapping/C1388.sort.bam 00.mapping/C1388.bam
```
1.4 利用picard移除除PCR duplicated reads
```
java -Xmx10g -jar /dt1/yinm/software/picard/picard.jar MarkDuplicates
		INPUT=00.mapping/C1075.sort.bam \
		OUTPUT=00.mapping/C1075.rmdup.bam \
		METRICS_FILE=00.mapping/C1075.dup.txt \
		REMOVE_DUPLICATES=true
```
1.5 构建bam文件索引
```
/dt1/share/software/samtools-1.18/samtools index 00.mapping/C1075.rmdup.bam
```
1.6 利用mosdepth统计样本的有效测序深度
```
mosdepth -t 4 00.mapping/C1 00.mapping/C1.rmdup.bam
```

  
### 2. SNP calling using GATK4
2.1 利用HaplotypeCaller将bam文件转化成gvcf
```
/dt1/share/software/gatk-4.4.0.0/gatk --java-options "-Xmx20g" HaplotypeCaller \
		-R /dt2/yinm/project/tobacco/00.ref_genome/NtaSR1_genome_1.0.fa \
		-I 00.mapping/C1075.rmdup.bam \
		-O 01.gvcf/C1075.g.vcf.gz \
		--emit-ref-confidence GVCF
```
2.2 利用CombineGVCFs合并GVCF（按染色体拆分并行）
```
/dt1/share/software/gatk-4.4.0.0/gatk CombineGVCFs \
		-R /dt2/yinm/project/tobacco/00.ref_genome/NtaSR1_genome_1.0.fa \
		-V 03.gvcf.list \
		-L Chr01  \
		-O 03.combineGVCF_split_Chr/03.Chr01_merge_220samples.gvcf.gz`
```
2.3 利用GenotypeGVCFs联合分型
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

  
## 3. Genotype Quality Control
3.1 利用VariantFiltration对SNP和INDEL分别进行硬过滤
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
3.2 根据每个样本的平均深度过滤 VCF 文件中的 SNP 和 Indel 位点。具体来说，它检查每个样本在每个位点的深度（DP 值），如果深度异常（低于平均深度的三分之一或高于三倍），则将该样本的基因型设置为缺失（./.）
```perl
perl filter_by_sample_depth.pl 00.summary_sample_depth.txt 05.filter/Chr01_220samples.INDEL.hardfilter.vcf > 05.filter/Chr01_220samples.INDEL.DPfilter.vcf
perl filter_by_sample_depth.pl 00.summary_sample_depth.txt 05.filter/Chr01_220samples.SNP.hardfilter.vcf > 05.filter/Chr01_220samples.SNP.DPfilter.vcf
# 该脚本同时过滤了多等位位点和没有“PASS”标签的位点

head -n 4 00.summary_sample_depth.txt                      
# C1075	11.35
# C1102	11.46
# C1103	12.81
# C1104	12.21

# filter_by_sample_depth.pl脚本的内容
############### filter_by_sample_depth.pl ###############
#!/usr/bin/perl
use strict;
use warnings;
# Get the path to the depth file
my ($depth_file,$vcf_file) =  @ARGV or die "Usage: $0 depth_file < input.vcf > output.vcf\n";
# Read the sample average depth information
open F1,"$depth_file";
my %sample_depth;
while (<F1>) {
    chomp;
    next if /^\s*$/; # Skip empty lines
    my ($sample, $avg_depth) = split;
    $sample_depth{$sample} = $avg_depth;
}
close F1;
# Process the VCF file
open F2,"$vcf_file";
my @samples;
while (<F2>) {
    chomp;
    if ( $_ =~ /^##/) {
        # VCF header lines, print as is
        print "$_\n";
    } elsif ($_ =~ /^#CHROM/) {
        # Extract sample names
        print "$_\n";
        my @fields = split(/\t/, $_);
        @samples = @fields[9..$#fields];
    } else {
        my @fields = split /\t/;
	next if $fields[6] ne "PASS";
	next if $fields[4] =~ /,/;
	my $format = $fields[8];
        my @format_fields = split(/:/, $format);
        my %format_index = map { $format_fields[$_] => $_ } 0..$#format_fields;
        # Check if DP field exists
        if (exists $format_index{'DP'}) {
            my $dp_index = $format_index{'DP'};
            for (my $i = 9; $i < @fields; $i++) {
                my $sample = $samples[$i - 9];
                my $avg_depth = $sample_depth{$sample};
                my $genotype = $fields[$i];
                next if $genotype =~ /\.\/\./; # Skip missing genotypes
                my @genotype_fields = split(/:/, $genotype);
                my $dp_value = $genotype_fields[$dp_index];
                $dp_value = 0 if $dp_value eq '.'; # Treat missing DP as 0
                # Check if DP value is within range
                if ($dp_value <= $avg_depth / 3 || $dp_value >= $avg_depth * 3) {
                    $fields[$i] = './.';
                }
            }
        } else {
            # If no DP field, set all genotypes to missing
            for (my $i = 9; $i < @fields; $i++) {
                $fields[$i] = './.';
            }
        }
        print join("\t", @fields), "\n";
    }
}
close F2;
######################################################################
```
3.3 合并每条染色体的vcf
```
/dt1/share/software/gatk-4.4.0.0/gatk  GatherVcfs -I 05.filter/Chr01_220samples.INDEL.final.recode.vcf -I 05.filter/Chr01_220samples.SNP.final.recode.vcf ......(etc.) -O tobacco220.vcf
```

3.4 计算样本的杂合度
```bash
# 将vcf转化成plink的格式
plink --vcf tobacco220.vcf --make-bed --chr-set 24 --recode --out tobacco220
# 计算杂合率
plink --bfile tobacco220 --het --out tobacco220

# 计算结果 .het 文件
#   FID     IID       O(HOM)       E(HOM)        N(NM)            F
# C1075   C1075      2966558    2.137e+06      2969718       0.9962
# C1102   C1102      2973284    2.142e+06      2975577       0.9973
```
F列的值越大，表示样本纯和度越高，负值代表样本是高杂合的。对于自交系群体，可以考虑过滤掉高杂合样本\

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

3.6过滤高杂合的样本，以及高缺失率和低频的位点
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
  
在实际分析中需要根据自己的分析目的进行样本和位点的过滤。由于我的群体是一个自交系群体，因此我过滤掉了杂合的样本（F < 0）。不同的分析分析应该执行不同的过滤，特别是MAF。例如，计算核苷酸多样性应该保留稀有或低频的位点。此外，在对VCF进行过滤的时候，过滤顺序应该是先进行样本过滤，后执行位点过滤。所以，一般先对样本的缺失率、杂合度等进行计算，再对位点的缺失率、基因型频率、杂合率等进行过滤。\
    
