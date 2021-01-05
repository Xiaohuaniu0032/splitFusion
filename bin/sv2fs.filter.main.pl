use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use FindBin qw/$Bin/;

# 从SV结果中筛选可靠的fusion结果

# 1. sv_annotate_filter_ctDNA_v3.py
# 2. check_split_read_uniq.py
# 3. this script (make final fusion result)

# 2021/1/4

my ($svFile,$bam,$bed,$annotFile,$fsDBDIR,$blist,$samtools,$outdir,$py3,$mapq,$refFlat);
GetOptions(
	"sv:s"         => \$svFile,       # Need
	"bam:s"        => \$bam,          # Need
	"bed:s"        => \$bed,          # Need
	"annot:s"      => \$annotFile,    # Default: /data1/workdir/wangce/software/svscan/svdb/dna/anndb/anno.sort.rmchr.gz
	"fsDBDIR:s"    => \$fsDBDIR,      # Default: /data1/workdir/wangce/software/svscan/svdb/fusion
	"blist:s"      => \$blist,        # Default: /data1/workdir/wangce/software/svscan/svdb/fusion/fusion.blacklist.txt
	"samtools:s"   => \$samtools,     # Default: /data1/workdir/wangce/software/samtools-1.8/samtools
	"od:s"         => \$outdir,       # Need
	"py3:s"        => \$py3,          # Default: /home/fulongfei/miniconda3/bin/python3
	"mapq:i"       => \$mapq,         # Default: 10
	"refFlat:s"    => \$refFlat,      # Default: /home/fulongfei/workdir/git_repo/GeExCNV/public_db/refFlat.txt
	) or die;


# default value

if (not defined $annotFile){
	$annotFile = "/data1/workdir/wangce/software/svscan/svdb/dna/anndb/anno.sort.rmchr.gz";
}

if (not defined $fsDBDIR){
	$fsDBDIR = "/data1/workdir/wangce/software/svscan/svdb/fusion";
}

if (not defined $blist){
	$blist = "/data1/workdir/wangce/software/svscan/svdb/fusion/fusion.blacklist.txt";
}

if (not defined $samtools){
	$samtools = "/data1/workdir/wangce/software/samtools-1.8/samtools";
}

if (not defined $py3){
	$py3 = "/home/fulongfei/miniconda3/bin/python3";
}

if (not defined $refFlat){
	$refFlat = "/home/fulongfei/workdir/git_repo/GeExCNV/public_db/refFlat.txt";
}

if (not defined $mapq){
	$mapq = 10;
}

# print args into output
print "############ ARGS ############\n";
print "SV file is: $svFile\n";
print "BAM file is: $bam\n";
print "BED file is: $bed\n";
print "Annot file is: $annotFile\n";
print "Fusion DBDIR is: $fsDBDIR\n";
print "Black list is: $blist\n";
print "samtools is: $samtools\n";
print "Output DIR is: $outdir\n";
print "Python3 bin is: $py3\n";
print "MAPQ cutoff is: $mapq\n";
print "refFlat file is: $refFlat\n";
print "############ ARGS ############\n";

# main filter process

# outfile is *.svscan.fusion.xls
my $cmd = "$py3 $Bin/sv_annotate_filter_ctDNA_v3.py -sv $svFile -bam $bam -bed $bed -annot $annotFile -fsDBDIR $fsDBDIR -blist $blist -samtools $samtools -od $outdir";
system($cmd);


# outfile is *.sr_mapq_uniqness_check.txt
my $name = (split /\./, basename($svFile))[0];
my $fs = "$outdir/$name\.svscan.fusion.xls";
$cmd = "$py3 $Bin/check_split_read_uniq.py -bam $bam -f $fs -refFlag $refFlat -od $outdir";
system($cmd);

# outfile is *.uniq_check.result.txt
my $uniq_check_f = "$outdir/$name\.sr_mapq_uniqness_check.txt";
$cmd = "perl $Bin/infer_fs_false_pos_by_sr_uinq.pl $uniq_check_f";
system($cmd);

# outfile is *.svscan.fusion.xls
`mv $fs $outdir/$name\.svscan.fusion.temp.xls`;
my $fs_tmp = "$outdir/$name\.svscan.fusion.temp.xls";
my $aln_res = "$outdir/$name\.uniq_check.result.txt";
$cmd = "perl $Bin/add_align_uniq_info.pl $fs_tmp $aln_res $outdir";
system($cmd);

