# GPU-MatchLocate1.0
An improved match and locate method

1.The GPU-M&L package contains following main files:

	bin: executable file
	src: source code
	Demo: an example

2.Revision History:

	06/08 2018 M. Liu Initial coding
	14/08 2019 M. Liu Version 1.0 released


3.References:

	1)Zhang, M., and L. Wen (2015), An effective method for small event detection: match and locate (M&L), Geophysical Journal International, 200(3),1523-1537.
	2)Beaucé, E., W. B. Frank, and A. Romanenko (2017), Fast matched filter (FMF): An efficient seismic matched‐filter search for both CPU and GPU architectures, Seismological Research Letters, 89(1), 165-172.
	3)Liu, M. H. Li, M. Zhang, and T. Wang (2019), Graphics Processing Unit-based Match&Locate (GPU-M&L): An improved Match and Locate technique and its application. (submitted)

4.Usage:

	Please see user guide

5.Demo:
        
	5.1 compile source code
	$cd ./src
	$make
	$cd sacCC
	$make

  	5.2 download continuous data (User-specified)
	$cd ../../Demo
	$python data_download.py

	5.3 download templates and create routine catalog (User-specified)
	$python template_download.py

	5.4 set common origin time
	$cd ./Trace
	$perl SACH_O.pl 20190704

	5.5 mark P/S wave arrival time and calculate slowness
	$cd ../Template
	$perl marktaup_p.pl
	$perl marktaup_s.pl

	5.6 calculate weighting factor for each trace
	$bash cal_weights.sh

	5.7 detect and locate events
	$cd ..
	$perl RunprocAll.pl

	5.8 generate new catalog
	$cd ./MultipleTemplate
	$perl MergeEvents.pl 20190704
	$perl SelectFinal.pl 2019 07 04 Allevents

	5.9 compare waveforms between template and detection
	$perl PlotEventWaveformJapan.pl DetectedFinal.dat 1

6.Author:

	Min Liu, China University of Geosciences (Beijing), liumin@cugb.edu.cn
