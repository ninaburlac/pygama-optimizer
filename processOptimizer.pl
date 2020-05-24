#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use YAML::Tiny 'LoadFile';
use YAML::Tiny;
use Data::Dumper;
use POSIX qw(ceil floor);

use constant false => 0;
use constant true  => 1;

my $localdir = $ENV{PWD};

my $run = $ARGV[0];
my $filter = $ARGV[1];
my $window = $ARGV[2];
#my $date = $ARGV[2];

my @cmd = ( "!wind_command!", "!opt_command!");

my $myUser = $ENV{'USER'};

# save current umask
my $old_umask = umask;
umask 0000;
my $outdir = $localdir . "/_outs";
mkdir "$outdir", 0770 unless -d "$outdir";

my $wind_command = "";
if ($window){ 
    $wind_command = "python " . $localdir . "/optimizer.py -r " . $run . " -w -d " . $localdir . " > " . $outdir . "/optimizer_" . $run . "_wind.out";
}
my $opt_command = "python " . $localdir . "/optimizer.py -r " . $run . " -c " . $filter . " -g -p -f -t -d " . $localdir . " > " . $outdir . "/optimizer_" . $run . "_" . $filter . ".out";
    
my $scriptlog = $outdir . "/file_" . $run . "-" . $filter . ".out";
my $scripterr = $outdir . "/file_" . $run . "-" . $filter . ".err";
my $jobName = "op" . $run . "-" . $filter;

my $scripttobecopied = $localdir . "/script.template.sh";
open IN, $scripttobecopied or die "Can't read source file $scripttobecopied: $!\n";

my $currscript = $outdir . "/script_" . $run . "-" . $filter . ".sh";
open OUT, ">$currscript" or die "Can't write on file $currscript: $!\n";
while (<IN>) {
    s/$cmd[0]/$wind_command/g;
    s/$cmd[1]/$opt_command/g;
    print OUT $_;
}
close IN;
close OUT;

my $QUEUEcmd = "qsub -P short -N " . $jobName . " -e localhost:". $scripterr . " -o localhost:" . $scriptlog . " " . $currscript;
system($QUEUEcmd);
print $QUEUEcmd . "\n";


my $countJob = "qstat -u " . $myUser . " |  grep op" . $run . "-" . $filter . " | wc -l ";
print $countJob . "\n";
my $inQueue = 1;

while ($inQueue) {
    sleep 30;
    my $actualJob = `$countJob`;
    print $actualJob;
    if($actualJob==0) {
	$inQueue=0;
	
	umask $old_umask;
	
	my $from = 'optimizer';
	my $to = 'valerio.dandrea\@lngs.infn.it';
	my $subject = "optimizer-run" . $run;
	my $body = "Hi, \n this is an automatic message from pygama optimizer\n";
	my $sendMail = "echo \"". $body . "\" | mailx -s \"". $subject . "\" -r " . $from . " " .  $to;
	system($sendMail);
    }
}
