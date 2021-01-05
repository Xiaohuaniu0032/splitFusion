use strict;
use warnings;
use File::Basename;
use FindBin qw/$Bin/;

my ($svFile, $aln_uniq_result, $outdir) = @ARGV;

my $name = (split /\./, basename($svFile))[0];
my $outfile = "$outdir/$name\.svscan.fusion.xls";

open O, ">$outfile" or die;
#print O "\#SampleID\tFusionGene\tHeadGene\tTailGene\tFusionReadsNumber\tTotalReadsNumber\tFusionRate(%)\tCOSMIC\tOncoKB\tIfReport\tFilterReason\tGene1\tChr1\tJunctionPosition1\tStrand1\tTranscript1\tRegion1\tGene2\tChr2\tJunctionPosition2\tStrand2\tTranscript2\tRegion2\tFusionSequence\tsvType\tsvSize\tinsSeq\tinsLen\tFusionOriPattern\n";

# get align check result
my %aln_uniq;
open IN, "$aln_uniq_result" or die;
<IN>;
while (<IN>){
	chomp;
	my @arr = split /\t/;
	$aln_uniq{$arr[0]} = $arr[-1]; # gene=>Uniq/NonUniq
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
		if ($arr[9] eq "YES"){
			# check aln uniq
			my $aln_res = $aln_uniq{$arr[1]};
			#print "$arr[1]\t$aln_res\n";
			if ($aln_res ne "Uniq"){
				print "check $arr[1] fusion gene's align uniqness: Failed\n";
				my $FilterReason = "split read align uniqness check Failed";
				# update origin result
				$arr[9] = "NO";
				$arr[10] = $FilterReason;
			}else{
				print "check $arr[1] fusion gene's align uniqness: Successed\n";
			}
		}
	}
	
	my $flag = 0;
	for my $val (@arr){
		$flag += 1;
		if ($flag < scalar(@arr)){
			print O "$val\t";
		}else{
			print O "$val\n";
		}
	}
}
close IN;
close O;
			
				
