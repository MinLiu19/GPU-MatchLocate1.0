#!/bin/bash

#############User-specified########
#P/S arrival header		  #
p_phase="t2"			  #
s_phase="t1"			  #
#Signal window			  #
t1=-1				  #
t2=5				  #
#Noise window			  #
t3=-3				  #
t4=-1				  #
###################################

temp_signal_e="signal.dat"
temp_noise_e="noise.dat"
S_travel="S_travel.dat"
SNR_wf="wf_SNR.dat"
weights="weights.dat"

cat ../catalog.dat|gawk '{print $1}'>temp.list
for template in `cat temp.list`
do
	cd $template
	saclst $s_phase f *.*.*|gawk '{print 1/$2}' > $S_travel
	saclst $s_phase f *.*.*|gawk '{print "../../../bin/sac_e "$2+"'"$t1"'",$2+"'"$t2"'",$1}' | sh > $temp_signal_e
	saclst $p_phase f *.*.*|gawk '{print "../../../bin/sac_e "$2+"'"$t3"'",$2+"'"$t4"'",$1}' | sh > $temp_noise_e
	paste $temp_signal_e $temp_noise_e|gawk '{printf "%s %12.2f\n",$1,$2/$4}'|gawk '{if($2<1){print $1"	1";} else {print $1"	"$2}}' > $SNR_wf
	sum_snr=`cat $SNR_wf|gawk '{SUM+=(log($2)/log(10))}END{print SUM}'`
	sum_s_travel=`cat $S_travel|gawk '{SUM+=$1}END{print SUM}'`
	paste $SNR_wf $S_travel|gawk '{print $1,(((log($2)/log(10))/"'"$sum_snr"'")+($3/"'"$sum_s_travel"'"))/2}'>$weights
	paste ../INPUT/$template $weights|gawk '{print $1"	"$2"	"$3"	"$5}'>tmp
	mv tmp ../INPUT/$template
	rm *dat 
	cd ..
done
rm temp.list
