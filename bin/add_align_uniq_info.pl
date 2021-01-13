use strict;
use warnings;
use File::Basename;
use FindBin qw/$Bin/;

my ($svFile, $aln_uniq_result, $outdir) = @ARGV;

my $name = (split /\./, basename($svFile))[0];
my $outfile = "$outdir/$name\.svscan.fusion.xls";

open O, ">$outfile" or die;

# get align check result
my %aln_uniq;
open IN, "$aln_uniq_result" or die;
<IN>;
while (<IN>){
	chomp;
	my @arr = split /\t/;
	$aln_uniq{$arr[1]} = $arr[-1]; # gene=>Uniq/NonUniq
}
close IN;

# re-check those gene(s) if IfReport == YES
open IN, "$svFile" or die;
my $h = <IN>;
print O "$h";
while (<IN>){
	chomp;
	my @arr = split /\t/;
	if (exists $aln_uniq{$arr[1]}){
		if ($arr[7] eq 'Y' || $arr[8] eq 'Y'){
			# 对于常见融合（COSMIC & COCOKB数据库中）不检查
			print O "$_\n";
		}else{
			if ($arr[9] eq "YES"){
				# 检查需要报出的不常见融合
				my $aln_res = $aln_uniq{$arr[1]};
				if ($aln_res ne "Uniq"){
					print "check $arr[1] fusion gene's align uniqness: Failed\n";
					my $FilterReason = "split read align uniqness check Failed";
					$arr[9] = "NO";
					$arr[10] = $FilterReason;

					# output
					my $flag = 0;
					for my $val (@arr){
						$flag += 1;
						if ($flag < scalar(@arr)){
							print O "$val\t";
						}else{
							print O "$val\n";
						}
					}
				}else{
					print O "$_\n";
					print "check $arr[1] fusion gene's align uniqness: Successed\n";
				}	
			}else{
				print O "$_\n";
			}
		}
	}else{
		print O "$_\n";
	}
}
close IN;
close O;
