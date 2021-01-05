use strict;
use warnings;
use Parallel::ForkManager;

my ($file, $process_num) = @ARGV;

my %commands;
open IN, "<", "$file" or die;
while (<IN>){
    chomp;
    $commands{$_}=1;
}
close IN;

my $pm = Parallel::ForkManager->new($process_num);

foreach my $cmd (keys %commands){
    my $pid = $pm->start and next;
    my $cmd2 = "sh $cmd";
    &run($cmd2);
    $pm->finish;
}

$pm->wait_all_children;

sub run{
    my ($cmd) = @_;
    print STDERR "$cmd has been runed locally by system()\n";
    system($cmd);
}
