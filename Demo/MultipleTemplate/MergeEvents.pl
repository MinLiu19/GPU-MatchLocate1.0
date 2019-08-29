#/usr/bin/env perl 
use warnings;

$dir = "../";
$day = "$ARGV[0]";
$file = "Allevents";

open(FL,">$file");
chdir "$dir";
@events = glob "$day";
foreach $_(@events){
    chomp($_);
    open(JK,"<$_");
    @pars = <JK>;
    close(JK);
    $num = @pars-1;
    shift(@pars);
    for($i=0;$i<$num;$i++){
        chomp($pars[$i]);
        ($jk,$time,$lat,$lon,$dep,$mag,$coef,$mad,$ref) = split(" ",$pars[$i]);
        if($mad>=9){
	print FL "$time $lat $lon $dep $mag $coef $mad $ref\n";
    }}
}
close(FL);
