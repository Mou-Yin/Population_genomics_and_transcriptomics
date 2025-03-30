
#!/usr/bin/perl
use strict;
use warnings;
#use IO::Uncompress::Gunzip qw(gunzip $GunzipError);  # 导入 Gunzip 模块，可读取.gz文件

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
