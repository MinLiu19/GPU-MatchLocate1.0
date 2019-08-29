#!/bin/bash
# Event Format:
#     Date        evLat.    evLon.   evDep.  Mag. refLat.  refLon. refDep.
#20120902030045  37.799   139.998    7.8     2.1  37.799   139.998  7.8
dir1=$1
dir2=$2
#inputfile=$4

#Search area   (Maxlat/Maxlon/Maxdepth)
#"0/0/0" corresponds to the matched filter technique
R="0/0/0"

#Search step
I="0.001/0.001/1"

#Time Window/times MAD
#Only one detection will be kept in the time window(e.g.., 2 sec).
#User specify multiple times MAD as threshold
D="6.0/9";

#Template Window (windowlength/before/after)
#The cross-correlation window based on the marked t1 in your templates.
T="6.0/1/5"

#Station parameter
#1.Channel number of each template
#2.Template_dir dt_dD(horizontal slowness)/dt_dh(vertical slowness)
#3.Trace_dir
INPUT="INPUT.in"

#The number of Templates
N="28"

#The step of searching and delta (step/delta/segments)
G="1/0.01/1"

#Output style(0:no output stacked cc value 1:output the stacked cc value)
O="0"

starttime=`date`
printf "******Start Running GPU_Match&Locate on $starttime\n"
ls $dir2 >tmp1
for event_id in `cat catalog.dat|gawk '{print $1}'`
do
cat Template/INPUT/$event_id|gawk '{print $1}'>>tmp1
done
cat tmp1|sort|uniq -c|gawk '{if($1=="'$N'"+1)print $2}'>station.dat
cat station.dat|wc -l>$INPUT
for event in `cat catalog.dat|gawk '{print $1}'`
do
for station in `cat station.dat`
do
cat $dir1/INPUT/$event|gawk '{print "'$dir1/$event/'"$1"	"$3"	"$4}'|grep $station>>$INPUT
done
done
for station in `cat station.dat`
do
echo "$dir2/$station">>$INPUT
done

np=`echo $R/$I|gawk -F/ '{print (2*$1/$4+1)*(2*$2/$5+1)*(2*$3/$6+1)}'`
echo "There are $np potential locations!"
../bin/GPU_MatchLocate $INPUT -R$R -I$I -T$T  -D$D  -N$N -G$G -O$O
if [ $O == 1 ];then 
if [ ! -d ./stackdir ];then
mkdir stackdir
fi
mv *stack stackdir
fi
rm station.dat tmp1 
endtime=`date`
printf "******End Running GPU_Match&Locate on $endtime\n"
