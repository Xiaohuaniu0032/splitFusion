use strict;
use warnings;
use File::Basename;

my ($infile) = @ARGV;
# infile is *.sr_mapq_uniqness_check.txt
# outfile is *.uniq_check.result.txt

my $d = dirname($infile);
my $name = (split /\./, basename($infile))[0];
#print "$name\n";
my $outfile = "$d/$name\.uniq_check.result.txt";
print "make output file: $outfile\n";

open O, ">$outfile" or die;
print O "fsGene\thgene_pass_n\thgene_nopass_n\ttgene_pass_n\ttgene_nopass_n\tall_nopass_n\tn\tnopass_pct\tuniqCheckResult\n";

my %fs_gene_uniq_check;
my %hgene_info;

open IN, "$infile" or die;
<IN>;
while (<IN>){
	chomp;
	my @arr = split /\t/;
	$fs_gene_uniq_check{$arr[0]}{$arr[3]}{$arr[-1]} += 1; # gene->hgene/tgene->PASS/NOPASS += 1
}
close IN;

foreach my $gene (keys %fs_gene_uniq_check){
	my ($hgene_pass_n,$hgene_nopass_n);
	my ($tgene_pass_n,$tgene_nopass_n);

	if (exists $fs_gene_uniq_check{$gene}{'hgene'}){
		if (exists $fs_gene_uniq_check{$gene}{'hgene'}{'PASS'}){
			$hgene_pass_n = $fs_gene_uniq_check{$gene}{'hgene'}{'PASS'};
		}else{
			$hgene_pass_n = 0;
		}
		
		if (exists $fs_gene_uniq_check{$gene}{'hgene'}{'NOPASS'}){
			$hgene_nopass_n = $fs_gene_uniq_check{$gene}{'hgene'}{'NOPASS'};
		}else{
			$hgene_nopass_n = 0;
		}
	}else{
		$hgene_pass_n = 0;
		$hgene_nopass_n = 0;
	}

	if (exists $fs_gene_uniq_check{$gene}{'tgene'}){
		if (exists $fs_gene_uniq_check{$gene}{'tgene'}{'PASS'}){
			$tgene_pass_n = $fs_gene_uniq_check{$gene}{'tgene'}{'PASS'};
		}else{
			$tgene_pass_n = 0;
		}
		
		if (exists $fs_gene_uniq_check{$gene}{'tgene'}{'NOPASS'}){
			$tgene_nopass_n = $fs_gene_uniq_check{$gene}{'tgene'}{'NOPASS'};
		}else{
			$tgene_nopass_n = 0;
		}
	}else{
		$tgene_pass_n = 0;
		$tgene_nopass_n = 0;
	}

	# 判断规则:hgene & tgene必须同时包含PASS,且NOPASS占比<0.85
	
	my $align_uniq;
	
	my $nopass_n = $hgene_nopass_n + $tgene_nopass_n;
	my $n = $hgene_pass_n + $hgene_nopass_n + $tgene_pass_n + $tgene_nopass_n;
	my $nopass_pct = sprintf "%.2f", $nopass_n / $n;
	if ($hgene_pass_n > 1 and $tgene_pass_n > 1 and $nopass_pct < 0.85){
		$align_uniq = 'Uniq';
	}else{
		$align_uniq = 'NonUniq';
	}
	
	print O "$gene\t$hgene_pass_n\t$hgene_nopass_n\t$tgene_pass_n\t$tgene_nopass_n\t$nopass_n\t$n\t$nopass_pct\t$align_uniq\n";
}

close O;


