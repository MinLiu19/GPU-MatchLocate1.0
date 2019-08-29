#!/usr/bin/env perl
use warnings;
#The script is modified from Miao Zhang's M&L.
#Please merge all potential events into one file 'AllEvents', and select one best within INTD (e.g.,6 sec) using SelectFinal.
#Time is relative to the origin time of the data (I use the beginning time of one day).
@ARGV >= 2 || die "perl $0 Year Month Day Allevents\n";
my ($year,$month,$day,$eventfile) = @ARGV;
chomp($eventfile);
my $detect = "DetectedFinal.dat";

#keep one event within 6 sec.
my $D = "6.0";

system("../../bin/SelectFinal  -D$D $eventfile");

$out = $eventfile.".final";
open(FL,"< $out");
my @events = <FL>;
close(FL);

open(OUT,">$detect");
foreach $_(@events){
    chomp($_);
    if(substr($_,0,1) eq "#"){
        print OUT "#No.      Date        Time         Lat.      Lon.        Dep.     Mag.   Coef.      Times_MAD           Reference\n";
    }else{
        my ($num,$time,$lat,$lon,$dep,$mag,$coef,$r1,$ref) = split(" ",$_);
        my $hh = int($time/3600);
        my $min = int(($time-$hh*3600)/60);
        my $sec = $time - $hh*3600 - $min*60;
        printf OUT "%4d   %04d/%02d/%02d   %02d:%02d:%06.3f   %7.4f   %8.4f   %6.2f   %5.2f  %6.4f   %8.4f      %s\n",$num,$year,$month,$day,$hh,$min,$sec,$lat,$lon,$dep,$mag,$coef,$r1,$ref;
    }
}
unlink $out;
