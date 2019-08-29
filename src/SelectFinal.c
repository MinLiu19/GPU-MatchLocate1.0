/******************************************************************************************************
//Miao Zhang
//zhmiao@mail.ustc.edu.cn
//Select events in each day finally.
*********************************************************************************************************
//only one event with highest CC value is selected in a specific time window by the program.
//modified by Min Liu (2018.07)
//liumin_cugb@163.com
 *******************************************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <omp.h>
#include <unistd.h>
#include "sac.h"

#define NPEVE 10000000	//max possible determined events for all locaiton



//define sturcture to store the potential events.
typedef struct event {
    double eventLat;
    double eventLon;
    double eventH;
    double event_Times_Mad;
    double eventTime;
    double eventCoef;
    double eventMag;
    char eventRef[128];
} EVENT;

//Globle parameter
float INTD;
char outfile[128];

//Subroutine
void DetermineEvent(EVENT *loc, long int);
int Time_compare_Time(const void *a,const void *b);
int Time_compare_Coef(const void *a,const void *b);

//Main function begin here.
int main(int argc, char **argv) {
    int i,k,kk,error;
    FILE *fp1;
    char inputfile[128];
    extern float INTD;
    long int gk;
    EVENT *loc, *loc1;
    extern char outfile[128];
    float HHN;

    error = 0;
    for(i=1; !error && i<argc; i++) {
        if(argv[i][0] == '-') {
            switch(argv[i][1]) {
            /*case 'H':
                sscanf(&argv[i][2],"%f/%f",&THRESH1,&THRESH2);
                break;
            case 'D':
                sscanf(&argv[i][2],"%f",&INTD);
                break;
            case 'N':
                sscanf(&argv[i][2],"%f/%f/%f",&CC,&NN1,&NN2);
                break;
            default:
                error = 1;
                break;*/
	    case 'D':
                sscanf(&argv[i][2],"%f",&INTD);
                break;
	    default:
                error = 1;
                break;
            }
        }
	}
  /*  if(argc < 3 || error == 1) {
        fprintf(stderr,"Usage: SelectFinal -H(CC1/CC2) -N(CC0/NN1/NN2) -D(INTD) Eventfile\n");
        fprintf(stderr,"-H: segmented CC thresholds.\n");
        fprintf(stderr,"-N: background mean CC, segmented SNR coefficients.\n");
        fprintf(stderr,"-D: keep one event within INTD sec.\n");
        fprintf(stderr,"Eventfile: potential events detected by all templates.\n");
        return -1;
    }*/

    //Read the input file
    for(i=2; i<argc; i++) {
        strcpy(inputfile,argv[i]);
        sprintf(outfile,"%s.final",inputfile);
	if((fp1 = fopen(inputfile,"r")) == NULL) {
            fprintf(stderr,"Unable to open file [Event files] %s\n",inputfile);
            exit(-1);
        }

        loc = (EVENT *)malloc(sizeof(EVENT)*NPEVE);
        loc1 = (EVENT *)malloc(sizeof(EVENT)*NPEVE);
        if(loc == NULL || loc1 == NULL) {
            fprintf(stderr,"There are no enough memory for loc and loc1\n");
            exit(-1);
        }


        k = 0;
        kk = 0;
        //fprintf(stderr,"Start read all potential events!\n");
        while((fscanf(fp1,"%lf %lf %lf %lf %lf %lf %lf %s", &loc1[k].eventTime, &loc1[k].eventLat,
                      &loc1[k].eventLon, &loc1[k].eventH, &loc1[k].eventMag, &loc1[k].eventCoef,
                      &loc1[k].event_Times_Mad, loc1[k].eventRef)) != EOF) {
                loc[k].eventTime = loc1[k].eventTime;
                loc[k].eventLat = loc1[k].eventLat;
                loc[k].eventLon = loc1[k].eventLon;
                loc[k].eventH = loc1[k].eventH;
                loc[k].eventMag = loc1[k].eventMag;
                loc[k].eventCoef = loc1[k].eventCoef;
                loc[k].event_Times_Mad = loc1[k].event_Times_Mad;
                strcpy(loc[k].eventRef,loc1[k].eventRef);
            k++;}
        //fprintf(stderr,"End read all potential events!\n");

        gk = k;
qsort(loc,gk,sizeof(EVENT),Time_compare_Time);
        //fprintf(stderr,"Start select events!\n");
        DetermineEvent(loc,gk);
        //fprintf(stderr,"End select events!\n");

        //Free memory
        free(loc1);
        free(loc);
    }
    return 0;
}

void DetermineEvent(EVENT *loc,long int gk) {
    EVENT *loc2;
    extern char outfile[128];
    extern float INTD;
    int i,j,k,nn;
    FILE *fp;

    loc2 = (EVENT *)malloc(sizeof(EVENT)*gk);
    if(loc2 == NULL) {
        fprintf(stderr,"Can't locate loc2\n");
        exit(-1);
    }

    fp = fopen(outfile,"w");
    if(fp == NULL) {
        fprintf(stderr,"Can't open output file in DetermineEvent\n");
        exit(-1);
    }
    fprintf(fp,"#Event     Time     Lat.      Lon.        Depth    Mag.    Coef.      Times_MAD     Reference\n");

    if(gk > 0) {
        loc2[0].eventTime = loc[0].eventTime;
        loc2[0].eventLat = loc[0].eventLat;
        loc2[0].eventLon = loc[0].eventLon;
        loc2[0].eventH = loc[0].eventH;
        loc2[0].eventMag = loc[0].eventMag;
        loc2[0].eventCoef = loc[0].eventCoef;
        loc2[0].event_Times_Mad = loc[0].event_Times_Mad;
        strcpy(loc2[0].eventRef,loc[0].eventRef);

        k = 0;
        for(i=1; i<gk; i++) {
                if(fabs(loc[i].eventTime - loc2[k].eventTime) < INTD) {
			if(loc[i].eventCoef>loc2[k].eventCoef){
				loc2[k].eventTime = loc[i].eventTime;
                		loc2[k].eventLat = loc[i].eventLat;
                		loc2[k].eventLon = loc[i].eventLon;
                		loc2[k].eventH = loc[i].eventH;
                		loc2[k].eventMag = loc[i].eventMag;
                		loc2[k].eventCoef = loc[i].eventCoef;
                		loc2[k].event_Times_Mad = loc[i].event_Times_Mad;
               			strcpy(loc2[k].eventRef,loc[i].eventRef);
			}
		}
		else{
			k++;
                	loc2[k].eventTime = loc[i].eventTime;
                	loc2[k].eventLat = loc[i].eventLat;
                	loc2[k].eventLon = loc[i].eventLon;
                	loc2[k].eventH = loc[i].eventH;
                	loc2[k].eventMag = loc[i].eventMag;
                	loc2[k].eventCoef = loc[i].eventCoef;
                	loc2[k].event_Times_Mad = loc[i].event_Times_Mad;
                	strcpy(loc2[k].eventRef,loc[i].eventRef);
		}

        }
	k++;
        qsort(loc2,k,sizeof(EVENT),Time_compare_Time);

        for(i=0; i<k; i++) {
            fprintf(fp,"%4d   %9.3lf   %7.4f   %8.4f   %6.2f   %5.2f    %6.4f    %.2f    %s\n",i+1,
                    loc2[i].eventTime,loc2[i].eventLat,loc2[i].eventLon,loc2[i].eventH,loc2[i].eventMag,
                    loc2[i].eventCoef,loc2[i].event_Times_Mad,loc2[i].eventRef);
        }
    }
    fclose(fp);
    free(loc2);
}

//Compare function
int Time_compare_Time(const void *a,const void *b) {
    double v1 = ((EVENT *)a)->eventTime;
    double v2 = ((EVENT *)b)->eventTime;
    return (v1<v2) ? -1 : (v1>v2) ? 1 : 0;
}

//Compare function
int Time_compare_Coef(const void *b,const void *a) {
    double v1 = ((EVENT *)a)->eventCoef;
    double v2 = ((EVENT *)b)->eventCoef;
    return (v1<v2) ? -1 : (v1>v2) ? 1 : 0;
}
