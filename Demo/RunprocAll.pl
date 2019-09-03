#!/usr/bin/env perl
use warnings;
# Modified by Min Liu 09/2018
# From Miao Zhang's M&L pacakge
# Event Format:
#     Date        evLat.    evLon.   evDep.  Mag. refLat.  refLon. refDep.
#20120902030045  37.799   139.998    7.8     2.1  37.799   139.998  7.8
#Template dir
$dir1 = "./Template";
#Slave continuous data direction
$dir2 = "./Trace";
#Beginning date
$year0 = "2019";
$month0 = "07";
$day0 = "04";
#Day length
$dleng = "1"; 

for($i=0;$i<$dleng;$i++){
	if($i == 0){
	$year = $year0; $month=$month0; $day = $day0;
	}else{
	($year,$month,$day) = &Timeadd($year0,$month0,$day0,1);}
	
	$year0 = $year; $month0 = $month; $day0 = $day;
	print"$year$month$day\n";
	if(length($month)==1){$month = "0".$month;} 
	if(length($day)==1){$day = "0".$day;} 
	$outfile ="$year$month$day";
	system("bash GPU_MatchLocate.sh $dir1 $dir2/$year$month$day");
	system("mv EventCase.out $outfile");
}

sub Timeadd{
   my($yyear,$mm,$dday,$adday) = @_;
   $dday = $dday + $adday;	
   if (($mm==1) || ($mm==3) || ($mm==5) || ($mm==7) || ($mm==8) || ($mm==10) || ($mm==12)){
      if ($dday >31) {
         $dday = 1;
         $mm = $mm + 1;
         if ($mm > 12) {
            $mm = 1;
            $yyear = $yyear + 1;
         }
      }
   }    
   if (($mm==4) || ($mm==6) || ($mm==9) || ($mm==11)){
      if ($dday >30) {
         $dday = 1;
         $mm = $mm + 1;
         if ($mm > 12) {
            $mm = 1;
            $yyear = $yyear + 1;
         }
      }
   }    
   if ($mm == 2) {
      if ((($yyear%4 == 0) && ($yyear%100 != 0)) || ($yyear%400 == 0)){
         if ($dday >29) {
            $dday = 1;
            $mm = $mm + 1;
         }
      }
      else{
        if ($dday >28) {
            $dday = 1;
            $mm = $mm + 1;
         }
      }
   }

   my @time = ($yyear,$mm,$dday);
   return(@time);
}
