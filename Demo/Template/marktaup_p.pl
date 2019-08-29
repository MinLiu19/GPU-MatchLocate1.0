#!/usr/bin/env perl
use strict;
use warnings;
#Usually, only direct waves (e.g., Pg, Sg) are used in a regional area.
my $dir = ".";
my @files = `ls -d $dir/2019*`;
#P wave
my $pha = "P,p";

#S wave
#my $pha = "S,s"; 
my $distmax = "200";

#my $model = "JMA"; If you have model JMA2001.
my $model = "prem";

my $INPUT = "INPUT";
if(-e $INPUT){}else{`mkdir $INPUT`;}

my $temp = "temp";
open(SAC,"|sac");
system("bash modify.sh");
for(my $i=0;$i<@files;$i++){
	chomp($files[$i]);
	my ($jk,$FILE) = split("\/",$files[$i]);
	my $NAME = "$INPUT/$FILE";
	my @sta = `ls $files[$i]/*`;
	open(FL,">$NAME");
foreach $_(@sta){
	chomp($_);
	my ($jk0,$gcarc,$dist,$knetwk,$kstnm,$kcmpnm,$evdp) = split(" ",`saclst gcarc dist knetwk kstnm kcmpnm evdp f $_`);chomp($evdp);
	if($dist <= $distmax){
	print"$_\n";
	system("taup_time -mod $model -h $evdp -km $dist -ph $pha > $temp"); #See TauP.
	open(TP,"<$temp");
	my @par = <TP>;
	close(TP);
	my $t1 = -100;
	shift(@par);shift(@par);shift(@par);shift(@par);shift(@par);
	my ($jk1,$jk2,$jk3,$p,$takeoff,$incident,$jk4,$jk5,$jk6);
    ($jk1,$jk2,$jk3,$t1,$p,$takeoff,$incident,$jk4,$jk5,$jk6) = split(" ",shift(@par));
	my $dt_dgc = $p; #Ray parameter, horizontal slowness
	$takeoff = $takeoff*3.1415926/180;
    $t1 = sprintf("%.2f",$t1); #We have to change it to number of points in the code and make sure it can be devided by your sampling interval exactly.
    my $dt_dh = - ($p/111.11)/(sin($takeoff)/cos($takeoff)); #depth slowness
    #my $dt_dh = 0.0;
	if($t1 < 0 ){print STDERR "stop!! Wrong with the arrival time!\n"; exit;}
	print SAC "rh $_\n";
	print SAC "ch t2 $t1\n";
	print SAC "wh\n";
	unlink $temp;
	my $station =  "$knetwk.$kstnm.$kcmpnm";
	my $D = "$dt_dgc/$dt_dh";
	printf FL "%10s %10.3f %10.4f/%.4e\n",$station,$t1,$dt_dgc,$dt_dh;
	}
}
	close(FL);
}
print SAC "q\n";
close(SAC);
