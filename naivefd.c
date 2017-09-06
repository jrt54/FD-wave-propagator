#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "fd.h"
#include <time.h>
#include <unistd.h>
#include <getopt.h>

//Replace prev with expression for next
//For now, velocity is constant
//source injection is seperate
void Setup(float *prev, const int nx, const int ny, const int nz)
{
const int z_row=ny*nx;

//+pressure[idx+(-i-j+1)*z_row]//*2.0/(density[idx+(-i+1)*z_row]+density[idx-i*z_row])
#pragma omp parallel for collapse(2) 
for(int i=0; i<nx; ++i){
for(int j=0; j<ny; ++j){
for(int k=0; k<nz; ++k){
prev[i+j*nx+k*z_row]=0.0;
}}}
}

void Processargs(int argc, char** argv, int *fd_radius, int *nx, int *ny, int *nz, int *solver_option)
{

int c;
int digit_optind = 0;

           while (1) {
               int this_option_optind = optind ? optind : 1;
               int option_index = 0;
               static struct option long_options[] = {
                   {"radius",    required_argument, 0,  'o' },
                   {"solver",    required_argument, 0,  's' },
                   {"nx",         required_argument,                 0,  'x' },
                   {"ny",         required_argument,                 0,  'y' },
                   {"nz",         required_argument,                 0,  'z' }
               };

               c = getopt_long(argc, argv, "abc:d:012",
                        long_options, &option_index);
               if (c == -1)
                   break;

               switch (c) {

               case 'x':
                   *nx=atoi(optarg);
                   printf("option nx with value '%d'\n", *nx);
                   break;
               case 'y':
                   *ny=atoi(optarg);
                   printf("option ny with value '%d'\n", *ny);
                   break;

               case 'z':
                   *nz=atoi(optarg);
                   printf("option nz with value '%d'\n", *nz);
                   break;
               
               case 'o':
                   //printf("option space_order with value '%d'\n", *fd_radius);
                   *fd_radius=atoi(optarg);
                   printf("option space_order with value '%d'\n", *fd_radius);
                   break;

               case 's':
                   if (strcmp(optarg, "slow")==0) {
			*solver_option=0;
                   	printf("choosing to use slow, low-memory method \n");
                   	}
		   else if (strcmp(optarg, "pre")==0) {
			*solver_option=1;
                   	printf("choosing to precompute \n");
			}
                   else if (strcmp(optarg, "opt")==0) {
                   	printf("choosing to try new opt \n");
			*solver_option=2;
			}
                   break;
               case '?':
                   break;

               default:
                   printf("?? getopt returned character code 0%o ??\n", c);
               }
           }

           if (optind < argc) {
               printf("non-option ARGV-elements: ");
               while (optind < argc)
                   printf("%s ", argv[optind++]);
               printf("\n");
           }
}
int main(int argc, char** argv)
{

const float freq=15.0;
float vel=3000;
float Np=5.0;
float h=vel/(freq*Np); //grid spacing
const int nt=50;
//const float dt = .5*h/vel; //later insert CFL conditions
int nx = 100;
int ny = 100;
int nz = 100;
int fd_radius=8;//approximate derivative from x-r, ...x+r
int solver_option=2;
Processargs(argc, argv, &fd_radius, &nx, &ny, &nz, &solver_option);

const int total_size = nx*ny*(nz);
const int src_x = (int)nx/2;
const int src_y = (int)ny/2;
const int src_z = (int)nz/2;
const int src_idx = src_x + nx*src_y + nx*ny*src_z;

const float dt = h/(vel*3.0); //later insert CFL conditions
const int space_order=fd_radius*2+1;

float* full_coeff = (float*) malloc(space_order*sizeof(float));
float* grid_points = (float*) malloc(space_order*sizeof(float));
//constructing centered grid

for(int i=0; i<space_order; ++i){
    if (i<fd_radius)
	{
	grid_points[i]=(-fd_radius+.5+i)*h;
	}
    else if (i==fd_radius)
	{
	grid_points[i]=0.0;
	}
    else if (i>fd_radius)
	{
	grid_points[i]=(-fd_radius-.5+i)*h;
	}
}



const int diff_order=1; //taking first derivative
fd_coeff(full_coeff, 0.0, diff_order, grid_points, space_order);

float* coeff = (float*) malloc((fd_radius+1)*sizeof(float));
for(int i=0;i<=fd_radius;++i)
{
coeff[i]=full_coeff[fd_radius+i];
}

/*
float coeff[2]={0,1};

for(int i=0; i<=fd_radius;++i)
{
coeff[i]*=1.0/h;
}
*/

for(int i =0; i<=fd_radius; ++i){
printf("unscaled(by grid size) coeff[%d]=%f \n", i, coeff[i]*h);
}



float *source = (float*) malloc(sizeof(float)*nt);
Rickerwavelet(nt, dt, 15, source, vel);
for(int i=0; i<nt; ++i)
{
//printf("source[%d]: %f \n",i, source[i]);
}

float *density = (float*) malloc(sizeof(float)*total_size);
for(int i=0; i<total_size; ++i)
{
density[i]=1.0;
//printf("density[%d]: %f \n",i, density[i]);
}

float *pressure_prev = (float*) malloc(sizeof(float)*total_size);
//pressure_prev = memset(pressure_prev, 0, total_size*sizeof(pressure_prev[0]));
Setup(pressure_prev, nx, ny, nz);
/*for(int i=0; i<total_size; ++i)
{
printf("pressure_prev[%d]: %f \n",i, pressure[i]);
}*/


float *pressure_curr = (float*) malloc(sizeof(float)*total_size);
//pressure_curr = memset(pressure_curr, 0, total_size*sizeof(pressure_curr[0]));
Setup(pressure_curr, nx, ny, nz);

clock_t t;
time_t w;
t=clock();
w=time(NULL);
unsigned long long update_time=0;
if(solver_option==2){
printf("Using new opt option \n");
AllTimeStepOpt(pressure_prev, pressure_curr, density, nx, ny, nz, fd_radius, coeff, dt, vel, nt, source, src_idx, &update_time);
}
else if(solver_option==1){
printf("Using precompute option \n");
AllTimeStepPreCompute(pressure_prev, pressure_curr, density, nx, ny, nz, fd_radius, coeff, dt, vel, nt, source, src_idx, &update_time);
}
else if(solver_option==0){
printf("Using slowest, low-memory option \n");
AllTimeStep(pressure_prev, pressure_curr, density, nx, ny, nz, fd_radius, coeff, dt, vel, nt, source, src_idx);
}
t=clock()-t;
w=time(NULL)-w;
double time_taken=((double)t)/CLOCKS_PER_SEC; //in seconds
double clock_taken=((double)w); //in seconds
printf("timestepping took %f seconds clocktime \n", (float)time_taken);
printf("timestepping took %ld seconds walltime \n", (long)clock_taken);
printf("updating took %f seconds clocktime \n", (float)update_time/1000L);
printf("CLOCKS_PER_SEC is %lu \n", CLOCKS_PER_SEC);

{
int resolution=40;
int sizeofnx=sizeof(nx)*(sizeof(char)/sizeof(int)+1);//get size of level digit
int sizeofny=sizeof(ny)*(sizeof(char)/sizeof(int)+1);//get size of level digit
int sizeofnz=sizeof(nz)*(sizeof(char)/sizeof(int)+1);//get size of level digit
int sizeofnt=sizeof(nt)*(sizeof(char)/sizeof(int)+1);//get size of level digit
int sizeofopt=sizeof(solver_option)*(sizeof(char)/sizeof(int)+1);//get size of level digit
int sizeofnumbers=sizeofnx+sizeofny+sizeofnz+sizeofnt+sizeofopt;
char *plotname = malloc(sizeof(char)*(strlen("soln_nx=_ny=_nz=_nt=_opt=") +sizeofnumbers+ strlen(".png")));
//sprintf(plotname, "plots/soln_nx=%d_ny=%d_nz=%d_nt=%d_opt=%d.png", nx, ny, nz, nt, solver_option);
sprintf(plotname, "plots/soln_nx=%d_ny=%d_nz=%d_nt=%d.png", nx, ny, nz, nt);
PlotSolution(pressure_curr, nx, ny, nz, src_z, resolution, plotname);
}


return 0;

}
