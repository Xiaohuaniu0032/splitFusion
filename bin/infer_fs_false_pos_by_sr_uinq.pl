use strict;
use warnings;
use File::Basename;

my ($infile) = @ARGV;
# infile is *.sr_mapq_uniqness_check.txt
# outfile is *.uniq_check.result.txt

my $d = dirname($infile);
my $name = (split /\./, basename($infile))[0];
my $outfile = "$d/$name\.uniq_check.result.txt";
print "make output file: $outfile\n";

open O, ">$outfile" or die;
print O "sample\tfsGene\thgene_pass_n\thgene_nopass_n\ttgene_pass_n\ttgene_nopass_n\tpass_n(hgene+tgene)\tall_n(pass+nopass)\tpass_pct\tuniqness_check\n";

my %fs_gene_uniq_check;
my %hgene_info;

open IN, "$infile" or die;
<IN>;
while (<IN>){
	chomp;
	my @arr = split /\t/;
	# skip cigar_check == NOPASS
	if ($arr[6] eq 'NOPASS'){
		print "$_ skipped for cigar_check failed\n";
		next;
	}
	
	# skip mate chr == NOPASS
	if ($arr[8] eq 'NOPASS'){
		print "$_ skipped for mate_chr_check failed\n";
		next;
	}
	
	# skip sa_ori == NOPASS
	if ($arr[-5] eq 'NOPASS'){
		print "$_ skipped for sa_ori_check failed\n";
		next;
	}
	
	# the left is effective fusion-support reads
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

	
	# hgene_pass_n >= 2 &&  tgene_pass_n >= 2 (根据真阳性融合确定，真阳性融合基因hgene_pass_n和tgene_pass_n几乎都大于0，2是一个合理的阈值）
	# 同时，hgene_pass_n + tgene_pass_n / N > 0.01 （根据假阳性融合确定0.01）
	my $align_uniq;
	
	# some gene may only have hgene or tgene
	my $hg_n = $hgene_pass_n + $hgene_nopass_n;
	my $tg_n = $tgene_pass_n + $tgene_nopass_n;
	my $all_n = $hg_n + $tg_n;
	my $pass_n = $hgene_pass_n + $tgene_pass_n;
	
	my $pass_pct;
	if ($all_n == 0){
		$pass_pct = 0;
	}else{
		$pass_pct = sprintf "%.2f", ($hgene_pass_n + $tgene_pass_n) / $all_n;
	}
	
	if ($hgene_pass_n >= 2 and $tgene_pass_n >= 2 and $pass_pct > 0.01){
		$align_uniq = 'Uniq';
	}else{
		$align_uniq = 'NonUniq';
	}
	
	print O "$name\t$gene\t$hgene_pass_n\t$hgene_nopass_n\t$tgene_pass_n\t$tgene_nopass_n\t$pass_n\t$all_n\t$pass_pct\t$align_uniq\n";
}

close O;


