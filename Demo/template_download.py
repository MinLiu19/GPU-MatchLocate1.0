#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 16:23:58 2019

@author: mliu
"""

from obspy.clients.fdsn import Client
from obspy import UTCDateTime, read
import os

#User-specied 
starttime = UTCDateTime("2019-07-04T22:00:00.00")
endtime = UTCDateTime("2019-07-04T23:00:00.00")
mindepth = -100
maxdepth = 100
minmagnitude=2.5
maxmagnitude=6
latitude=35.7080
longitude=-117.5037
minradius=0
maxradius=0.5
tmpfile="tmp.catalog"
catalog="catalog.dat"
stationfile = "./IRIS.sta"
out = open(catalog,"w")



#request catalog
client = Client("SCEDC")
cat = client.get_events(catalog="SCEDC",starttime=starttime,endtime=endtime,mindepth=mindepth,maxdepth=maxdepth,minmagnitude=minmagnitude,maxmagnitude=maxmagnitude,latitude=latitude,longitude=longitude,minradius=minradius,maxradius=maxradius)
cat.write(tmpfile,format="CNV")

with open(tmpfile,"r") as events:
        for event in events:
            if(event != "\n"):
                client = Client("SCEDC")
                event = event.strip("\n")
                ymd, hm, s, lat, lon , dep, mag,jk = event.split()
                year="20"+ymd[:2]
                month=ymd[2:4]
                day=ymd[4:6]
                hour=hm[:2]
                minutes=hm[2:4]
                sec,msec=s.split(".")
                if (len(sec)<2):
                    sec="0"+sec
                evla=lat[:7]
                evlo=lon[:8]
                evdp=dep
                evmag=mag
                out.write('{}{}{}{}{}{}.{} {} -{} {} {} {} -{} {}\n'.format(year, month,day,hour,minutes,sec,msec,evla,evlo,evdp,evmag,evla,evlo,evdp))
                #download template
                client = Client("IRIS")
                tb = UTCDateTime(year+"-"+month+"-"+day+"T"+hour+":"+minutes+":"+sec+"."+msec)-10
                te = UTCDateTime(year+"-"+month+"-"+day+"T"+hour+":"+minutes+":"+sec+"."+msec)+50
                templatedir = "./Template/" + year+month+day+hour+minutes+sec+"."+msec + "/"
                if not os.path.exists(templatedir):
                    os.makedirs(templatedir)
                with open(stationfile, "r") as f:
                        for station in f:
                            
                            stlo, stla, net, sta, chan, elev = station.split()
                            chane = chan[:2] + "E"
                            chann = chan[:2] + "N"
                            chanz = chan[:2] + "Z"
                            chan0 = [chane, chann, chanz]
                            for chan1 in chan0:
                                print(sta,chan1)
                                try:
                                    st = client.get_waveforms(network=net, station=sta, channel=chan1,starttime=tb,endtime=te, location = '--')
                                    st.detrend("demean")
                                    st.detrend("linear")
                                    st.filter('bandpass', freqmin=2, freqmax=8)
                                    st.write(filename=templatedir + net + '.' + sta + '.' + chan1,format="SAC")
                                    st1 = read(templatedir+ net + '.' + sta + '.' + chan1)
                                    st1[0].stats.sac.stla=stla
                                    st1[0].stats.sac.stlo=stlo
                                    st1[0].stats.sac.stel=elev
                                    st1[0].stats.sac.evla=evla
                                    st1[0].stats.sac.evlo="-"+evlo
                                    st1[0].stats.sac.evdp=evdp
                                    st1[0].stats.sac.mag=evmag				    
                                    st1[0].stats.sac.nzsec+=10
                                    st1[0].stats.sac.b=-10
                                    st1.write(filename=templatedir + net + '.' + sta + '.' + chan1,format="SAC")
                                except:
                                    print("doesn't exist:",sta, chan1)
os.remove(tmpfile)
