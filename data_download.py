#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 17 22:58:57 2019

@author: mliu
"""

from obspy import UTCDateTime, read
from obspy.clients.fdsn import Client
import os

source = "IRIS"
stationfile = "./IRIS.sta"
tracedir = "./Trace/20190704/"
tb = UTCDateTime("2019-07-04T00:00:00.00")
te = tb+86400
client = Client(source)

if not os.path.exists(tracedir):
    os.makedirs(tracedir)

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
                st.write(filename=tracedir + net + '.' + sta + '.' + chan1,format="SAC")
                st1 = read(tracedir + net + '.' + sta + '.' + chan1)
                st1[0].stats.sac.stla=stla
                st1[0].stats.sac.stlo=stlo
                st1[0].stats.sac.stel=elev
                st1.write(filename=tracedir + net + '.' + sta + '.' + chan1,format="SAC")
            except:
                print("doesn't exist:",sta, chan1)
