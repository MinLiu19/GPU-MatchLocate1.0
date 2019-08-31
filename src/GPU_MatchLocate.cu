/*
	GPU-MatchLocate : GPU-based Match and Locate [The source code of GPU-M&L is mainly modified from M&L (Zhang and Wen, 2015)] 
	
	References:
		1.Zhang, M., and L. Wen (2015), An effective method for small event detection: match and locate (M&L), Geophysical Journal International, 200(3),1523-1537.
		2.Beaucé, E., W. B. Frank, and A. Romanenko (2017), Fast matched filter (FMF): An efficient seismic matched‐filter search for both CPU and GPU architectures, Seismological Research Letters, 89(1), 165-172.
		3.Liu, M. H. Li, M. Zhang, and T. Wang (2019), Graphics Processing Unit-based Match&Locate (GPU-M&L): An improved Match and Locate technique and its application. (submitted)
	
	Usage:
		see the usage below.
	Author: Min Liu, China University of Geosciences (Beijing) (liumin@cugb.edu.cn)
	Revision History
	06/08 2018 M. Liu Initial coding
	14/08 2019 M. Liu Version 1.0 released
 */
#include <cuda.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <unistd.h>
#include <omp.h>
#include <cuda_runtime.h>
#include <time.h>
#include <sys/sysinfo.h>
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
//GPU SETTING
#define BLOCKSIZE 512
#define WARPSIZE 32
#define STABILITY_THRESHOLD 0.000001f
#define D2R  .017453292519943295769237
#define R2D    57.2957795130823208768
#define OUTFILE "EventCase.out"
#define MARKT 1
extern "C" {
#include "sac.h"
}

//define sturcture to store the potential events.-----------------------------------------------------
typedef struct event {
    double eventLat;
    double eventLon;
    double eventH;
    double eventTime;
    double eventCoef;
    float  event_times_MAD;
} EVENT;
//----------------------------------------------------------------------------------------------------
//Function declear -----------------------------------------------------------------------------------

void Checkdata(int, char**, int*, float*);
int Time_compare_F(const void *a,const void *b);
int Time_compare_D(const void *a,const void *b);
int compare (const void *a, const void *b);
void calc_moveout(double, double, double, int *, int, int, int);
void GCinit(const double *lat1,const double *lon1,const double *lat2,const double *lon2,double *GCarc);
void Select_event(double evla, double evlo, double evdp, float *cc_sums, EVENT **loc2, long int *KK,int n_corr_final,int n_grid, float threshold,int id, int i);
float Mag(float initime,double lat, double lon, double h, int template_id);
float CalculateMedian(float *arrValue, int max);
void DetermineEvent(EVENT *loc,long int NN, int template_id);
float Threshold_detection(float *cc_sum, int npts);
int Time_compare_Coef(const void *a,const void *b);
int Time_compare_Time(const void *a,const void *b);
extern "C" int read_sachead(const char *, SACHEAD *);
extern "C" float   *read_sac(const char *, SACHEAD *);


//----------------------------------------------------------------------------------------------------
//Global parameter------------------------------------------------------------------------------------
char **traces,**templates,**ref_template;
double delay, delta;
double *Dt_Dgc,*Dt_Dh,*tshift,*gstla,*gstlo;
double *evla0,*evlo0,*evdp0,*mag0;
double *temp_lon,*temp_lat,*temp_h;
int ntrace,times;
double tb;
double template_window,before,after;
double INTD;
int n_weights;
float median,MAD;
int event_id;
int out_style;
//----------------------------------------------------------------------------------------------------
//check the fault of GPU card

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true) {

    if (code != cudaSuccess) 
    {
        fprintf(stderr, "An error occured in the kernel: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}
//----------------------------------------------------------------------------------------------------
//kernel of slide_corss_correlation-------------------------------------------------------------------
__global__ void slide_cross_correlation(float *templates, 
				       float *sum_square_template, 
				       float *data,
				       size_t step,
				       size_t n_samples_template,
				       size_t n_samples_data,
				       size_t n_stations, 
				       size_t n_components,
				       int n_corr,
				       float *cc_mat){
	//each thread matches the template to one time in the data
	int idx, first_sample_block, first_sample_trace, last_sample_trace;
	int i, s, c; 
	int data_offset, templates_offset, sum_square_template_offset, cc_mat_offset;
 	float numerator, denominator, sum_square_data;
    	float data_sample;
    	int t_idx;
	//clock_t start = clock();
	int count_template = (n_samples_template / WARPSIZE + 1) * WARPSIZE;
	extern __shared__ float shared[];
	float *templates_s = &shared[0];
	float *data_s = &shared[count_template];
	
	// 1 block processes one channel to blockDim.x / step different positions in time
        idx = blockIdx.x/n_stations * blockDim.x;//the point of n_corr
        first_sample_block = idx * step;
	s = blockIdx.x % n_stations;

	for (c = 0; c < n_components; c++){
		cc_mat_offset = (first_sample_block / step + threadIdx.x) * n_stations * n_components + s * n_components + c;
		templates_offset = s * n_samples_template * n_components + c * n_samples_template;
            	sum_square_template_offset = s * n_components + c;
		first_sample_trace = first_sample_block;
		last_sample_trace = first_sample_trace + n_samples_template + threadIdx.x * step;
		data_offset = s * n_samples_data * n_components + c * n_samples_data + first_sample_trace;
//initialize sums
		sum_square_data = 0.0f;
		numerator = 0.0f;


//load template and data into shared memory
		t_idx = threadIdx.x;
		while(t_idx < n_samples_template) {
                	templates_s[t_idx] = templates[templates_offset + t_idx];
                	if ((first_sample_trace + t_idx) < n_samples_data) data_s[t_idx] = data[data_offset + t_idx];
                		t_idx += blockDim.x;
            		}
            	while(t_idx < (blockDim.x * step + n_samples_template)){
                	if ((first_sample_trace + t_idx) < n_samples_data) data_s[t_idx] = data[data_offset + t_idx];
                		t_idx += blockDim.x;
            		}
		__syncthreads(); // make sure the waveforms are read before keep going

// calculate correlation coefficient
		if (last_sample_trace < n_samples_data){
// if not, corresponds to an ill-defined CC with some samples out of the bounds
			for(i = 0; i < n_samples_template; i++) {
				data_sample = data_s[i + threadIdx.x * step];
                    		numerator += data_sample * templates_s[i];
                    		sum_square_data += data_sample * data_sample; 			
			}
			denominator = sum_square_data * sum_square_template[sum_square_template_offset];
			if (cc_mat_offset < n_corr * n_stations * n_components){
				// check that this thread is not ouf of the chunk's bounds
				
				 if (denominator > STABILITY_THRESHOLD) cc_mat[cc_mat_offset] = numerator * rsqrtf(sum_square_data * sum_square_template[sum_square_template_offset]);
				
				}
			}
		__syncthreads(); // wait for every thread to finish before leaving the kernel 
		}
	//*time = clock() - start;
	
}	
//----------------------------------------------------------------------------------------------------
//kernel of shifting and stacking---------------------------------------------------------------------------
__global__ void sum_cross_correlation(float *cc_mat, 
				      float *cc_sum,
				      int *moveout,
				      float *weights,
				      int n_stations,
				      int n_components, 
				      int n_corr, 
				      int n_grid,
				      int n_weights,
				      int id_segment,
				      int n_corr_final){
	int id,idx,idy, ch, moveout_offset;
	idx = blockIdx.x * blockDim.x + threadIdx.x;
	idy = blockIdx.y;
	id = blockIdx.x * blockDim.x + threadIdx.x +idy*n_corr;

	if (((idx + n_corr*id_segment) < n_corr_final) && (idx < n_corr) && idy<n_grid){
        // first condition: check if we are not outside cc_sum's length
        // second condition: check if we are not outside the chunk's size
		float *cc_mat_offset;
		cc_mat_offset = cc_mat + idx * n_stations * n_components + n_stations * n_components*id_segment*n_corr;
		for (ch = 0; ch < (n_stations * n_components); ch++){ 
			moveout_offset = moveout[ch+idy*n_stations*n_components]*n_stations*n_components ;		
			cc_sum[id] += (weights[ch]*cc_mat_offset[ch+moveout_offset]);
		}    
			
	}

}
//----------------------------------------------------------------------------------------------------
int main (int argc, char **argv){

/*------------------------------------------Variable definition-------------------------------------------------*/

//1.Host
	extern double tb,delay,template_window,after,before;
	extern int out_style;
	out_style=0;
	char inputfile[256];
	int i,ii,iii,j, k, m, n_samples_data, n_samples_template,*max_moveout;
	int error=0;
	float t_data;
	SACHEAD hd,hd0;
	float **ref,**obj,*ref_p,*obj_p,*sum_square_ref;
	FILE *fp;
	FILE *fp1, *fp2;
	int step,segment_size;	
	int n_stations, n_components, n_corr; 
	int n_templates;
	//float ppp;
	float *weights;
	int SEGMENTS;
//2.Search grid	
	double *lat0,*lon0,*h0,dlat,dlon,dh,lat,lon,h;
	double maxlat,maxlon,maxh;
	int nlat,nlon,ndep;
	int *moveouts;
	n_weights=0;

//3.Device
	float *templates_d = NULL;
        float *data_d = NULL;
	float *cc_mat_d = NULL;
	float *sum_square_templates_d = NULL;
	int *gpu_moveout=NULL;
	float *gpu_weights=NULL;
	int template_id;
	float *cc_sum_d;
        int n_corr_final;
        long int pn;
        long int *gk;
        float *cc_sums_t,*cc_sums_tt;
        float *templates_d_t;
        float *sum_square_templates_d_t;
        float *gpu_weights_t;
        int *gpu_moveout_t;
        int maxSharedMem;
        int count_template;
        int count_data;
        int sharedMem;
        float threshold;
        float *cc_sums_stack;
/*---------------------------------------------------------------------------------------------------------------------*/ 

/*--------------------------------------------Load parameter----------------------------------------------------------------*/
	for(i=0;!error && i<argc; i++){
		if(argv[i][0] == '-'){
			switch(argv[i][1]){
				case 'R':
					sscanf(&argv[i][2],"%lf/%lf/%lf",&maxlat,&maxlon,&maxh);
					break;
				case 'I':
                			sscanf(&argv[i][2],"%lf/%lf/%lf",&dlat,&dlon,&dh);
                			break;
            			case 'T':
                			sscanf(&argv[i][2],"%lf/%lf/%lf",&template_window,&before,&after);
                			break;
				case 'D':
               				sscanf(&argv[i][2],"%lf/%d",&INTD,&times);
                			break;	
				case 'N':
					sscanf(&argv[i][2],"%d",&n_templates);
					break;
				case 'G':
					sscanf(&argv[i][2],"%d/%lf/%d",&step,&delay,&SEGMENTS);
					break;
				case 'O':
					sscanf(&argv[i][2],"%d",&out_style);				
					break;
				default:
               				error = 1;
                			break;
			}		
		}	
	}
/*---------------------------------------------------------------------------------------------------------------------------*/

/*-----------------------------------------------------------Usage----------------------------------------------------------*/
	if(argc < 9 || error == 1) {
		fprintf(stderr, "Usage: GPU_MatchLocate -R(maxlat/maxlon/maxh) -I(dlat/dlon/dh) -T(template_window/before/after -D(INTD/Times) -N(n_templates) -G(step/delay) -O(0 or 1) INPUT.in\n");
		fprintf(stderr, "-R: searching area (e.g., 0.05/0.05/5.0).\n");
		fprintf(stderr, "-I: searching interval (e.g., 0.01/0.01/1.0).\n");
		fprintf(stderr, "-T: time length of the reference phase (e.g., 4.0/1.0/3.0).\n");
		fprintf(stderr, "-D: keep one event within INTD sec and set threshold of detection (e.g., 6/9).\n");
		fprintf(stderr, "-N: the number of template.\n");
		fprintf(stderr, "-G: moveout setting for step and delay (e.g., 1/0.01).\n");
		fprintf(stderr, "-O: if output the stack ccv(e.g. default 0:no output).\n");
		fprintf(stderr, "INPUT.in: template and continuous seismograms, horizontal slowness and vertical slowness of the reference phase and weighting factor of each component.\n");
		return -1;
	}
/*---------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------Read the INPUT.IN---------------------------------------------------*/
	strcpy(inputfile, argv[1]);
	if((fp1 = fopen(inputfile, "r")) == NULL) {
		fprintf(stderr,"Unable to open file [INPUT.in] %s\n", inputfile);
		exit(-1);
	}
	if((fp2 = fopen("catalog.dat", "r")) == NULL) {
		fprintf(stderr,"Unable to open file catalog.dat \n");
                exit(-1);
	}
	fscanf(fp1, "%d", &ntrace);
	if(ntrace <= 0){
                fprintf(stderr,"Please check your %s!! ntrace = %d <=zero!!\n", inputfile,ntrace);
                exit(-1);
        }
/*---------------------------------------------------------------------------------------------------------------------------*/

/*------------------------------------------------Check data's shape --------------------------------------------------------*/
	n_components=3;
	if(ntrace%n_components!=0){
		fprintf(stderr,"The shape of data is wrong, Please check the number of stations\n");
		exit(-1);
	}
	else n_stations=ntrace/n_components;
/*---------------------------------------------------------------------------------------------------------------------------*/

/*--------------------------------------------------Memory request from Host-------------------------------------------------*/
	templates = (char **)malloc(sizeof(char*)*ntrace*n_templates);
	traces = (char **)malloc(sizeof(char*)*ntrace);
	ref_template = (char **)malloc(sizeof(char*)*n_templates);
	for(i=0; i<ntrace*n_templates; i++){
		templates[i]=(char *)malloc(sizeof(char)*256);
	}
	for(i=0; i<ntrace; i++){
		traces[i]=(char *)malloc(sizeof(char)*256);
	}
	for(i=0; i<n_templates; i++){
		ref_template[i]=(char *)malloc(sizeof(char)*256);
	}
	Dt_Dgc    =  (double *)malloc(sizeof(double)*ntrace*n_templates);
    	Dt_Dh     =  (double *)malloc(sizeof(double)*ntrace*n_templates);
    	tshift    =  (double *)malloc(sizeof(double)*ntrace*n_templates);
    	gstla     =  (double *)malloc(sizeof(double)*ntrace*n_templates);
    	gstlo     =  (double *)malloc(sizeof(double)*ntrace*n_templates);
	evla0     =  (double *)malloc(sizeof(double)*n_templates);
	evlo0     =  (double *)malloc(sizeof(double)*n_templates);
	evdp0     =  (double *)malloc(sizeof(double)*n_templates);
	temp_lat  =  (double *)malloc(sizeof(double)*n_templates);
        temp_lon  =  (double *)malloc(sizeof(double)*n_templates);
        temp_h    =  (double *)malloc(sizeof(double)*n_templates);
	mag0      =  (double *)malloc(sizeof(double)*n_templates);
	weights   =  (float *)malloc(sizeof(float)*ntrace*n_templates);
/*---------------------------------------------------------------------------------------------------------------------------*/

/*---------------------------------------------Read catalog-------------------------------------------------------------*/
	for(i=0; i<n_templates; i++){
		fscanf(fp2,"%s %lf %lf %lf %lf %lf %lf %lf", ref_template[i], &temp_lat[i], &temp_lon[i], &temp_h[i], &mag0[i], &evla0[i], &evlo0[i], &evdp0[i]);
	}
/*----------------------------------------------------------------------------------------------------------------------*/

/*--------------------------------------Arrange the weight for each component-------------------------------------------*/
	for(i=0; i<ntrace*n_templates; i++){
		fscanf(fp1, "%s %lf/%lf %f", templates[i], &Dt_Dgc[i],&Dt_Dh[i],&weights[i]);
	}
	for(i=0; i<ntrace; i++ ){
                fscanf(fp1, "%s", traces[i]);
	}
 	fclose(fp2);	
	fclose(fp1);
/*-----------------------------------------------------------------------------------------------------------------------*/

/*--------------------------------------------------------Check data-----------------------------------------------------*/
	Checkdata(ntrace, traces, &n_samples_data, &t_data);
	n_samples_template=template_window/delay;	
	n_corr=(n_samples_data-n_samples_template)/step;
/*------------------------------------------------------------------------------------------------------------------------*/


/*-----------------------------------------Load template and continues data into host memory------------------------------*/
	//read ref and obj	
	ref = (float **)calloc(ntrace*n_templates,sizeof(float *));
	for(i=0; i<ntrace*n_templates; i++)ref[i] = (float *)calloc(n_samples_template,sizeof(float));
    	obj = (float **)calloc(ntrace,sizeof(float *));
    	for(i=0; i<ntrace; i++)obj[i] = (float *)calloc(n_samples_data,sizeof(float));	
		for(i=0; i<ntrace; i++){
			if((obj[i] = read_sac(traces[i],&hd))==NULL) {
            		fprintf(stderr,"Can't open traces file %s\n",traces[i]);
            		exit(-1);
        	}
		tb = hd.b;
	}
	for(i=0; i<ntrace*n_templates; i++){
		if((ref[i] = read_sac2(templates[i],&hd0,MARKT,-before,after))==NULL) {
            		fprintf(stderr,"Can't open templates file %s\n",templates[i]);
            		exit(-1);
		}
	    	tshift[i] =  hd0.t1;
	    	gstla[i]  =  hd0.stla;
            	gstlo[i]  =  hd0.stlo;
	}
	ref_p=(float *)calloc(n_samples_template*ntrace*n_templates,sizeof(float));
	obj_p=(float *)calloc(n_samples_data*ntrace,sizeof(float));
	sum_square_ref=(float *)calloc(ntrace*n_templates,sizeof(float));
/*-----------------------------------------------------------------------------------------------------------------------*/

/*-------------------------Flaten the tempalte and trace(station*component*npts)-----------------------------------------*/
	for (i=0;i<ntrace; i++){
		for(j=0;j<n_samples_data;j++){
			obj_p[i*n_samples_data+j]=obj[i][j];
		}
	}
	for (i=0;i<ntrace*n_templates; i++){
                for(j=0;j<n_samples_template;j++){
           		ref_p[i*n_samples_template+j]=ref[i][j];
			sum_square_ref[i]+=ref[i][j]*ref[i][j];
                }
        }
/*------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------Generate search grids--------------------------------------------------*/
	
	nlat = (int)(2*maxlat/dlat + 1);
	nlon = (int)(2*maxlon/dlon + 1);
	ndep = (int)(2*maxh/dh + 1);
	moveouts=(int *)malloc(sizeof(int)*nlat*nlon*ndep*ntrace*n_templates);
	lat0  =  (double *)malloc(sizeof(double)*n_templates);
    	lon0  =  (double *)malloc(sizeof(double)*n_templates);
	h0    =  (double *)malloc(sizeof(double)*n_templates);
	for(i=0; i<n_templates; i++){
		lat0[i] = evla0[i] - maxlat;
		lon0[i] = evlo0[i] - maxlon;
		h0[i]   = evdp0[i] - maxh;	
		for(m=0; m<nlat*nlon*ndep; m++){
			ii = (int)(m/(nlon*ndep));
			j = (int)((m - ii*nlon*ndep)/ndep);
			k = m - ii*nlon*ndep - j*ndep;
			lat = lat0[i] + ii*dlat;
			lon = lon0[i] + j*dlon;
			h = h0[i] + k*dh;
			calc_moveout(lat,lon,h,moveouts,m,i,nlat*nlon*ndep);
		}
	}
	max_moveout = (int *)malloc(sizeof(int)*n_templates);
	for(j=0;j<n_templates;j++){
		max_moveout[j]=0;	
		for(i=0; i<ntrace*nlat*nlon*ndep;i=i+3){
			if(max_moveout[j]<moveouts[i+j*ntrace*nlat*nlon*ndep])max_moveout[j]=moveouts[i+j*ntrace*nlat*nlon*ndep];
		}
	}
/*------------------------------------------------------------------------------------------------------------------------*/	
	
/*--------------------------Size of variables to create on the device (GPU)-----------------------------------------------*/
	size_t sizeof_templates = sizeof(float) * n_samples_template * n_stations * n_components * n_templates;
	size_t sizeof_data = sizeof(float) * n_samples_data * n_stations * n_components;
	size_t sizeof_sum_square_templates = sizeof(float) * n_stations * n_components * n_templates;
	size_t sizeof_cc_mat = sizeof(float) * (n_corr) * n_stations * n_components;
	size_t sizeof_weights = sizeof(float) * ntrace * n_templates;	
	size_t sizeof_moveout = sizeof(int) * ntrace*nlat*nlon*ndep*n_templates;
	size_t sizeof_total = sizeof_templates + sizeof_data + sizeof_sum_square_templates + sizeof_cc_mat + sizeof_weights + sizeof_moveout;
	//size_t sizeof_bandpass = sizeof_templates + sizeof_data + sizeof_sum_square_templates;
	//printf("size= %zu\n",sizeof_bandpass);		
	cudaSetDevice(0);
	cudaDeviceProp props;
	cudaGetDeviceProperties(&props, 0);	
	cudaDeviceSetCacheConfig(cudaFuncCachePreferShared);
/*------------------------------------------------------------------------------------------------------------------------*/

	
/*------------------------------------------------Check memory------------------------------------------------------------*/
//1.check host memory													
	size_t freeMem = 0;
	size_t totalMem = 0;
	cudaMemGetInfo(&freeMem, &totalMem);
        if (sizeof_total > freeMem) {
        	printf("%zu bytes are requested on GPU #0 whereas it has only %zu free bytes.\n", sizeof_total, freeMem);
            	printf("Reduce the number of templates processed in one batch.\n");
            	exit(0);
        }
//2.check device memory
	maxSharedMem = props.sharedMemPerBlock;
	count_template = (n_samples_template / WARPSIZE + 1) * WARPSIZE;
	count_data = ((n_samples_template + BLOCKSIZE * step) / WARPSIZE + 1) * WARPSIZE;
	sharedMem = (count_template + count_data) * sizeof(float);	
	if (sharedMem > maxSharedMem){
		printf("The maximum shared memory available (%i bytes) on this card is not enough to search!\
			please reduce the samples of template", maxSharedMem);
		exit(0);
		}
/*------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------- allocate GPU memory-----------------------------------------------*/
        cudaMalloc((void**)&templates_d, sizeof_templates);
        cudaMalloc((void**)&data_d, sizeof_data);
        cudaMalloc((void**)&cc_mat_d, sizeof_cc_mat);
        cudaMalloc((void**)&sum_square_templates_d, sizeof_sum_square_templates);
	cudaMalloc((void**)&gpu_moveout,sizeof_moveout);
	cudaMalloc((void**)&gpu_weights,sizeof_weights);
/*------------------------------------------------------------------------------------------------------------------------*/
/*-----------------------------------Transfer data from Host memory to Device memory--------------------------------------*/
	cudaMemcpy(gpu_weights, weights, sizeof_weights, cudaMemcpyHostToDevice);
	cudaMemcpy(gpu_moveout, moveouts, sizeof_moveout, cudaMemcpyHostToDevice);
	cudaMemcpy(templates_d, ref_p, sizeof_templates, cudaMemcpyHostToDevice);
        cudaMemcpy(data_d, obj_p, sizeof_data, cudaMemcpyHostToDevice);
        cudaMemcpy(sum_square_templates_d, sum_square_ref, sizeof_sum_square_templates, cudaMemcpyHostToDevice);
/*---------------------------------------Free useless Memory--------------------------------------------------------------*/
	free(obj_p);
        free(ref_p);
        for(i=0; i<ntrace; i++) {
                free(obj[i]);
        }
        free(obj);
        for(i=0; i<ntrace * n_templates; i++) {
                free(ref[i]);
        }
        free(ref);
/*------------------------------------------------------------------------------------------------------------------------*/
	fp = fopen(OUTFILE,"a");
    	if(fp == NULL) {
        	fprintf(stderr,"Can't open output file in DetermineEvent\n");
        	exit(-1);
    	}
    	fprintf(fp,"#Event   Time      Lat.      Lon.        Depth    Mag.    Coef.    Times_MAD	Reference\n");	
	fclose(fp);
	template_id=0;
/*-------------------------------------------Loop Templates---------------------------------------------------------------*/
	while(template_id < n_templates){
		printf("The template\t%s\tis working and\t%d%%\ttasks are accomplished!\n",ref_template[template_id],100*(template_id+1)/n_templates);
		n_weights=0;
		for(i=0;i<ntrace;i++){
			if(weights[i+template_id * n_stations * n_components]>0){
				n_weights += 1;
			}
		}
		cc_sum_d = NULL;
//1.SCC's length of each template
		n_corr_final = (n_samples_data-n_samples_template-max_moveout[template_id])/step;
		segment_size=n_corr_final/SEGMENTS + 0.5;
		size_t sizeof_cc_sum = sizeof(float) *segment_size*nlat*nlon*ndep;
//2.Check device memory again
		printf("%d %d\n",(n_samples_data-n_samples_template-max_moveout[template_id])/step,n_corr_final);
		cudaMemGetInfo(&freeMem, &totalMem);
        	if (sizeof_cc_sum > freeMem) {
            		printf("%zu MB are requested on GPU #0 whereas it has only %zu free MB on device memory. Please increase the SEGMENTS\n", (sizeof_cc_sum/1024)/1024, (freeMem/1024)/1024);
            		exit(0);
        	}		
		cudaMalloc((void**)&cc_sum_d, sizeof_cc_sum);
		EVENT *loc,**loc1;
		cc_sums_t = NULL;//store the results
		templates_d_t = NULL;
        	sum_square_templates_d_t = NULL;
        	gpu_weights_t = NULL;
		gpu_moveout_t = NULL;			
		templates_d_t = templates_d + template_id * n_samples_template * n_stations * n_components;
        	sum_square_templates_d_t = sum_square_templates_d +template_id * n_stations * n_components;	
		cc_sums_tt = NULL;//store the results
		gpu_moveout_t = gpu_moveout + template_id * n_stations * n_components * nlat*nlon*ndep;
		gpu_weights_t = gpu_weights + template_id * n_stations * n_components;
//3.Check Host memory again
		pn= (int)(t_data*2/SEGMENTS+1);
		//ppp = (2*nlat*nlon*ndep*pn*sizeof(EVENT)+2*nlat*nlon*ndep*segment_size*sizeof(float)+sizeof(long int)*nlat*nlon*ndep)/(1024*1024);
                struct sysinfo si;
                sysinfo(&si);
		//printf("%d\n",pn*sizeof(EVENT)/(1024*1024));
		//if(ppp>si.freeram/(1024*1024)){
              //		fprintf(stderr,"%.0lf MB. Memory are requested on Host, but there are only %d MB.\n",ppp,int(si.freeram/(1024*1024)));
		//	exit(0);
		//}
//4.Request memory
		cc_sums_t = (float *)calloc(nlat*nlon*ndep*segment_size,sizeof(float ));
		cc_sums_tt = (float *)calloc(nlat*nlon*ndep*segment_size,sizeof(float));
                gk=NULL;
		loc=NULL;
		loc1=NULL;
		gk = (long int *)malloc(sizeof(long int)*nlat*nlon*ndep);
                loc = (EVENT *)malloc(sizeof(EVENT)*nlat*nlon*ndep*pn);
		loc1 = (EVENT **)malloc(sizeof(EVENT)*nlat*nlon*ndep);
                if(loc == NULL || loc1 == NULL){
			printf("There are errors in Memory, Please check it! (NT: Increase the segments.)\n");
			exit(-1);
		}
		for(i=0; i<nlat*nlon*ndep; i++) {
                        loc1[i] = (EVENT *)malloc(sizeof(EVENT)*pn);
                }
//5.Calculate SCCs
		cudaMemset(cc_mat_d, 0, sizeof_cc_mat);
		dim3 BS(BLOCKSIZE);
		dim3 GS(ceilf(n_corr / (float)BS.x)* n_stations);
		slide_cross_correlation<<<GS, BS, sharedMem>>>(templates_d_t,
					    sum_square_templates_d_t,	
					    data_d,
					    step,
					    n_samples_template,
					    n_samples_data,
					    n_stations,
					    n_components,
					    n_corr,	
					    cc_mat_d);
			
		gpuErrchk(cudaPeekAtLastError());
        	gpuErrchk(cudaDeviceSynchronize());

//
//6.Loop SEGMENTS to stack SCCs
		for(i=0;i<SEGMENTS;i++){//2018.4.29
			dim3 BS_sum(BLOCKSIZE);
			dim3 GS_sum(ceilf(segment_size/ (int)BS_sum.x),nlat*nlon*ndep);
			cudaMemset(cc_sum_d, 0, sizeof_cc_sum);
			sum_cross_correlation<<<GS_sum, BS_sum>>>(cc_mat_d, 
				cc_sum_d,
				gpu_moveout_t,
				gpu_weights_t,
				n_stations,
				n_components,
				segment_size,
				nlat*nlon*ndep,
				n_weights,
				i,
				n_corr_final);
			gpuErrchk(cudaPeekAtLastError());
        		gpuErrchk(cudaDeviceSynchronize());
//7.Transfer stacked SCCs from device to host
			cudaMemcpy(cc_sums_tt, cc_sum_d, sizeof_cc_sum,	cudaMemcpyDeviceToHost);
			gpuErrchk(cudaPeekAtLastError());
        		gpuErrchk(cudaDeviceSynchronize());
			for(j=0; j<nlat*nlon*ndep; j++){
				for(ii=0;ii<segment_size;ii++){
					cc_sums_t[j*segment_size+ii]=cc_sums_tt[ii+j*segment_size];
				}
			}
//8.Calculate the threshold of detection
			cc_sums_stack = NULL;
			cc_sums_stack=cc_sums_t+nlat*nlon*ndep/2*segment_size;
			threshold = Threshold_detection(cc_sums_stack, segment_size);
//9.Loop grids to select events
			#pragma omp parallel for shared(loc1,gk,cc_sums_t,segment_size,i,template_id) firstprivate(nlat,nlon,ndep,dlat,dlon,dh) private(m,iii,j,k,lat,lon,h)
			for(m=0;m<nlat*nlon*ndep;m++){
				iii = (int)(m/(nlon*ndep));
				j = (int)((m - iii*nlon*ndep)/ndep);
				k = m - iii*nlon*ndep - j*ndep;
				lat = lat0[template_id] + iii*dlat;
				lon = lon0[template_id] + j*dlon;
				h   = h0[template_id] + k*dh;
				Select_event(lat,lon,h,cc_sums_t,loc1,&gk[m],segment_size,m,threshold, template_id, i);
				}
			#pragma omp barrier
//10.Merge selected events
			m = 0;
			for(iii=0; iii<nlat*nlon*ndep; iii++) {
				for(j=0; j<gk[iii]; j++) {
					loc[m].eventTime = loc1[iii][j].eventTime;
					loc[m].eventCoef = loc1[iii][j].eventCoef;
					loc[m].eventLat  = loc1[iii][j].eventLat;
					loc[m].eventLon  = loc1[iii][j].eventLon;
					loc[m].eventH    = loc1[iii][j].eventH;
					loc[m].event_times_MAD	 = loc1[iii][j].event_times_MAD;
					m++;
                                }
                        }
			
//11.Select again
			DetermineEvent(loc,m,template_id);
		}
		cudaDeviceSynchronize();
		template_id++;	
//12.Free local memory
        	cudaFree(cc_sum_d);
		free(cc_sums_t);
        	free(cc_sums_tt);
		free(gk);
		free(loc);
		for(i=0; i<nlat*nlon*ndep; i++) free(loc1[i]);
		free(loc1);
	}//loop end	
/*---------------------------------------------------------------------------------------------------------------*/

/*--------------------------------------------Free memory GPU----------------------------------------------------*/
//1.Device memory
	cudaFree(cc_sum_d);
	cudaFree(templates_d);
	cudaFree(data_d);
	cudaFree(sum_square_templates_d);
	cudaFree(gpu_moveout);
	cudaFree(cc_mat_d);	
	cudaFree(gpu_weights);
//2.Host Memory
	free(h0);
	free(lat0);
	free(lon0);	
	free(max_moveout);
 	free(evla0);
	free(evlo0);
	free(evdp0);
	free(mag0);
 	free(Dt_Dgc);
    	free(Dt_Dh);
    	free(tshift);
    	free(gstla);
    	free(gstlo);
	free(weights);
	free(moveouts);	
	free(sum_square_ref);
	for(i=0; i<ntrace * n_templates; i++) {
        	free(templates[i]);
    	}

	free(templates);
	for(i=0; i<ntrace; i++) {
       		free(traces[i]);
    	}
	free(traces);
	for(i=0; i<n_templates; i++){
		free(ref_template[i]);
	}
	free(ref_template);
	return 0;
/*---------------------------------------------------------------------------------------------------------------*/
}
//calculate the threshold of detection (this subroutine is mainly modified from Prof. Peng's CPU_WFCC)
float Threshold_detection(float *cc_sum, int npts){
	extern float median, MAD;	
	extern int n_weights,times;
	float *cc_sum2,*cc_sum1;
	cc_sum1=(float *)calloc(npts,sizeof(float ));
	for(int i=0;i<npts;i++)cc_sum1[i]=cc_sum[i];	
	cc_sum2=(float *)calloc(npts,sizeof(float ));	
	qsort(cc_sum1, npts, sizeof(float), compare);
	if(npts%2 == 0)median = (cc_sum1[npts/2]+cc_sum1[npts/2-1])/2;	
	if(npts%2 == 1)median = cc_sum1[(npts-1)/2];
	for(int i=0; i<npts; i++)cc_sum2[i] = fabsf(cc_sum1[i]-median);
	qsort(cc_sum2, npts, sizeof(float), compare);	
	if(npts%2 == 0)MAD = (cc_sum2[npts/2]+cc_sum2[npts/2-1])/2;	
	if(npts%2 == 1)MAD = cc_sum2[(npts-1)/2];
	free(cc_sum1);
	free(cc_sum2);
	return median+times*MAD;	
	}
//select again
void DetermineEvent(EVENT *loc,long int NN, int template_id) {
	EVENT *loc2;
    	int i,j,k,nn;
    	FILE *fp;
    	float mag;
    	extern double INTD;
    	extern int event_id;
    	extern char **ref_template;
    	loc2 = (EVENT *)malloc(sizeof(EVENT)*NN);
    	if(loc2 == NULL) {
        	fprintf(stderr,"Can't locate loc2 in DetermineEvent\n");
        	exit(-1);
    	}
    //Sort by coefficient (from large to small).
    	qsort(loc,NN,sizeof(EVENT),Time_compare_Coef);
    	fp = fopen(OUTFILE,"a");
    	if(fp == NULL) {
        	fprintf(stderr,"Can't open output file in DetermineEvent\n");
        	exit(-1);
    	}
    	if(NN > 0) {
        	loc2[0].eventTime = loc[0].eventTime;
        	loc2[0].eventCoef = loc[0].eventCoef;
        	loc2[0].eventLat = loc[0].eventLat;
        	loc2[0].eventLon = loc[0].eventLon;
        	loc2[0].eventH = loc[0].eventH;
		loc2[0].event_times_MAD = loc[0].event_times_MAD;
        	k = 1;
        	for(i=1; i<NN; i++) {
            		nn = 0;
            		for(j=0; j<k; j++) {
                		if(fabs(loc[i].eventTime - loc2[j].eventTime) > INTD) {
                    		nn++;
                	}
            	}
            		if(nn == k || loc[i].eventCoef == 1.00) {
                		loc2[k].eventTime = loc[i].eventTime;
                		loc2[k].eventCoef = loc[i].eventCoef;
                		loc2[k].eventLat = loc[i].eventLat;
                		loc2[k].eventLon = loc[i].eventLon;
                		loc2[k].eventH = loc[i].eventH;
				loc2[k].event_times_MAD = loc[i].event_times_MAD;
                		k++;
            		}

        	}	

//Sort by origin time (from small to large).
        	qsort(loc2,k,sizeof(EVENT),Time_compare_Time);
        	for(i=0; i<k; i++) {
            		mag = Mag(loc2[i].eventTime,loc2[i].eventLat,loc2[i].eventLon,loc2[i].eventH,template_id);
            		fprintf(fp,"%4d   %9.3lf   %7.4f   %8.4f   %6.2f   %5.2f    %6.4f    %.2f	%s\n",++event_id,
                    	loc2[i].eventTime,loc2[i].eventLat,loc2[i].eventLon,loc2[i].eventH,mag,loc2[i].eventCoef,loc2[i].event_times_MAD,ref_template[template_id]);
		}
    	}
    	fclose(fp);
    	free(loc2);
}
//mag
float Mag(float initime,double lat, double lon, double h, int template_id) {
    	double mag,t1,t2,Dt;
    	float *master,*ar,*ratio;
    	double GCarc,GCarc0;
    	extern int ntrace;
    	extern double *gstla,*gstlo,*Dt_Dgc,*Dt_Dh,*mag0;
    	extern double before,template_window,after;
    	int i,k,npts1,npts2;
    	float armax,mastermax,median;
    	SACHEAD hd0,hd1;
    	int tmark2 = -3;
    	ratio = (float *)malloc(ntrace*sizeof(float));
    	for(i=0; i<ntrace; i++) {
        	if( (master=read_sac2(templates[i+template_id*ntrace],&hd0,MARKT,-before,after)) == NULL) {
            		fprintf(stderr,"erro master in Mag\n");
            		exit(-1);
        	}
        	npts1 = hd0.npts;
        	mastermax = master[0];
        	for(k=0; k<npts1; k++) {
            		if(fabs(master[k]) > mastermax)mastermax = fabs(master[k]);
        	}
        	GCinit(&gstla[i+template_id*ntrace],&gstlo[i+template_id*ntrace],&lat,&lon,&GCarc);
        	GCinit(&gstla[i+template_id*ntrace],&gstlo[i+template_id*ntrace],&evla0[template_id],&evlo0[template_id],&GCarc0);
        	Dt = Dt_Dgc[i+template_id*ntrace]*(GCarc-GCarc0) + Dt_Dh[i+template_id*ntrace]*(h-evdp0[template_id]);
        	t1 = initime + tshift[i+template_id*ntrace] + Dt - before;
        	t2 = t1 + template_window;
        	if( (ar=read_sac2(traces[i],&hd1,tmark2,t1,t2)) == NULL) {
            		fprintf(stderr,"erro event in Mag\n");
           		// exit(-1);
        	}
        	npts2 = hd1.npts;
        	armax = ar[0];
        	for(k=0; k<npts2; k++) {
            		if(fabs(ar[k]) > armax)armax = fabs(ar[k]);
        	}
        	ratio[i] = armax/mastermax; 
	}
    	free(master);
    	free(ar);
    	median = CalculateMedian(ratio,ntrace);
    	mag = mag0[template_id] + log10(median);
    	free(ratio);
    	return mag;
}
//median
float CalculateMedian(float *arrValue, int max)
{
    	float median = 0;
    	float *value;
    	int i, j;
    	float temp;
    	value = (float *)malloc(max*sizeof(float));
    	for(i = 0; i < max; i++)
                value[i] = arrValue[i];
    	for(i = 0; i < max; i++){
        	for(j = 0; j < max - i - 1; j++){
            		if(value[j] > value[j + 1]){
                		temp = value[j];
                		value[j] = value[j + 1];
                		value[j + 1] = temp;
            		}
        	}
    	}
    	if( (max % 2) == 1){
        	median =  value[ (max + 1) / 2 - 1];
    	}
    	else{
        	median = (value[max / 2] + value[max / 2 - 1]) / 2;
    	}
    	free(value);
    	return median;
}
//
//select event
void Select_event(double evla, double evlo, double evdp, float *cc_sums, EVENT **loc2, long int *KK,int n_corr_final,int n_grid, float threshold, int id, int segment_size_id){
	SACHEAD hd;
	int j,k;
	float *cc_sum;
	extern char **ref_template;
	extern float median,MAD;
	extern int out_style;
	extern double tb;
	extern double delay,before;
	char outfile[128];
	cc_sum = cc_sums + n_grid * n_corr_final;
	hd = sachdr(delay,n_corr_final,tb+before);
	hd.evla = evla;
	hd.evlo = evlo;
	hd.evdp = evdp;
	if(out_style==1){	
		sprintf(outfile,"%s_%lf.2_%lf.2_%lf.3.stack",ref_template[id],evla,evlo,evdp);
		write_sac(outfile,hd,cc_sum);
	}
	
	k=0;
	for(j=0;j<hd.npts;j++){
		if(cc_sum[j]>threshold){
				loc2[n_grid][k].eventTime = (double)(hd.b) + delay*j + delay*segment_size_id*n_corr_final;
				loc2[n_grid][k].eventCoef = cc_sum[j];
				loc2[n_grid][k].eventLat = evla;
				loc2[n_grid][k].eventLon = evlo;
				loc2[n_grid][k].eventH = evdp;
				loc2[n_grid][k].event_times_MAD = (cc_sum[j] - median)/MAD;				
				k++;
			}
		}
	*KK=k;
}

//calculate each grid's moveout------------------------------------------------
void calc_moveout(double lat, double lon, double h, int *moveout, int offset, int offset_template, int total_grids){
	extern int ntrace;
	extern double *Dt_Dgc,*Dt_Dh,*tshift;
	extern double *temp_lat, *temp_lon, *temp_h, *gstla,*gstlo;
	int i;
	extern double delay;
	double GCarc,GCarc0;
	extern double *Dt_Dgc,*Dt_Dh;
	for(i=0;i<ntrace;i++){
		GCinit(&gstla[i+offset_template*ntrace],&gstlo[i+offset_template*ntrace],&lat,&lon,&GCarc);
		GCinit(&gstla[i+offset_template*ntrace],&gstlo[i+offset_template*ntrace],&temp_lat[offset_template],&temp_lon[offset_template],&GCarc0);
		moveout[offset*ntrace+i+offset_template*total_grids*ntrace] =(int)((tshift[i+offset_template*ntrace] + Dt_Dgc[i+offset_template*ntrace]*(GCarc-GCarc0) + Dt_Dh[i+offset_template*ntrace]*(h-temp_h[offset_template]))/delay+0.5);
	}
}
//compute great circle
void GCinit(const double *lat1,const double *lon1,const double *lat2,const double *lon2,double *GCarc)
{
    	double x1,yp1,z1,x2,y2,z2;
    	double the1,phe1,the2,phe2;
    	the1=(90.0-*lat1)*D2R;         /* convert to radius */
    	phe1=(*lon1)*D2R;
    	the2=(90.0-*lat2)*D2R;
    	phe2=(*lon2)*D2R;
    	x1=sin(the1)*cos(phe1);
    	yp1=sin(the1)*sin(phe1);
    	z1=cos(the1);
    	x2=sin(the2)*cos(phe2);
    	y2=sin(the2)*sin(phe2);
    	z2=cos(the2);
    	*GCarc=acos((x1)*(x2)+(yp1)*(y2)+(z1)*(z2));
    	if( fabs(*GCarc-M_PI) <=1.e-16) {
        	fprintf(stderr," The great circle is not determined! (GCinit)\n");
        	exit(1);
    	}
    	*GCarc *= R2D;
}
//---------------------------------------------------------------------------------

int Time_compare_Time(const void *a,const void *b) {
    	double v1 = ((EVENT *)a)->eventTime;
    	double v2 = ((EVENT *)b)->eventTime;
    	return (v1<v2) ? -1 : (v1>v2) ? 1 : 0;
}
int Time_compare_Coef(const void *b,const void *a) {
    	double v1 = ((EVENT *)a)->eventCoef;
    	double v2 = ((EVENT *)b)->eventCoef;
    	return (v1<v2) ? -1 : (v1>v2) ? 1 : 0;
}

int Time_compare_F(const void *a,const void *b) {
   	float v1 = *(float *)a;
   	float v2 = *(float *)b;
    	return (v1<v2) ? -1 : (v1>v2) ? 1 : 0;
}

int Time_compare_D(const void *a,const void *b) {
    	int v1 = *(int *)a;
    	int v2 = *(int *)b;
  	return (v1<v2) ? -1 : (v1>v2) ? 1 : 0;
}

int     compare (const void *a, const void *b) {
        const float *da = (const float *) a;
        const float *db = (const float *) b;
        return  (*da > *db) - (*da < *db);
}

// check the continuous data, you can ignore below part.
void Checkdata(int ntrace, char **traces, int *npp, float *tleng) {
   	float *nnb,*nne;
   	double *nnd;
  	int *nnl, i;
  	extern double delta;
  	SACHEAD hd;
  	nnb     =  (float *)malloc(sizeof(float)*ntrace);
    	nne     =  (float *)malloc(sizeof(float)*ntrace);
    	nnd     =  (double *)malloc(sizeof(double)*ntrace);
    	nnl     =  (int *)malloc(sizeof(int)*ntrace);
    	for(i=0; i<ntrace; i++) {
        	read_sachead(traces[i],&hd);
         	nnb[i] = hd.b;
         	nne[i] = hd.e;
         	nnd[i] = (double)(hd.delta);
         	nnl[i] = (int)(((hd.e - hd.b))/hd.delta+0.5);
         	*tleng = hd.e - hd.b;
         	if(hd.o != 0) {
            		fprintf(stderr,"Please set origin time ZERO [ch o 0] in your continuous data.\n");
            		exit(-1);
         	}
    	}
    	qsort(nnb,ntrace,sizeof(nnb[0]),Time_compare_F); // sort from small to large
    	qsort(nne,ntrace,sizeof(nne[0]),Time_compare_F);
    	qsort(nnd,ntrace,sizeof(nnd[0]),Time_compare_F);
    	qsort(nnl,ntrace,sizeof(nnl[0]),Time_compare_D);
    	if(fabs(nnb[ntrace-1] - nnb[0]) > 1.0e-3) {
        	fprintf(stderr,"Please Check Your Continuous Data!!! Not the same begin time.\n");
        	exit(-1);
    	}
    	if(fabs(nnd[ntrace-1] - nnd[0]) > 1.0e-5) {
        	fprintf(stderr,"Please Check Your Continuous Data!!! Not the same sampling invertal.\n");
        	exit(-1);
    	}
    	*npp = nnl[0];
    	if(fabs(nne[ntrace-1] - nne[0]) > 1.0e-3) {
        	fprintf(stderr,"Please Check Your Continuous Data!!! Not the same end time.\n");
        	fprintf(stderr,"Shortest end time: %.2f sec; Longest end time %.2f sec. The shortest one would be used.\n",nne[ntrace-1],nne[0]);
        	//exit(-1); //make sure the code don't stop, please reference above warning.
    }
    	delta = nnd[0];
    	free(nnb);
    	free(nne);
    	free(nnd);
    	free(nnl);
 	//All data length should be equal.
     	//If there is a gap in continuous seismogram, please set the gap value to be zero.
}
