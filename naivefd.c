#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "fd.h"
#include <time.h>


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


int main()
{

const float freq=15.0;
float vel=3000;
float Np=5.0;
float h=vel/(freq*Np); //grid spacing
const int nt=50;
//const float dt = .5*h/vel; //later insert CFL conditions
const float dt = h/(vel*3.0); //later insert CFL conditions
const int nx = 100;
const int ny = 100;
const int nz = 100;
const int total_size = nx*ny*(nz);
const int src_x = (int)nx/2;
const int src_y = (int)ny/2;
const int src_z = (int)nz/2;
const int src_idx = src_x + nx*src_y + nx*ny*src_z;


const int fd_radius=8;//approximate derivative from x-r, ...x+r
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
AllTimeStepOpt(pressure_prev, pressure_curr, density, nx, ny, nz, fd_radius, coeff, dt, vel, nt, source, src_idx, &update_time);
//AllTimeStepPreCompute(pressure_prev, pressure_curr, density, nx, ny, nz, fd_radius, coeff, dt, vel, nt, source, src_idx, &update_time);
//AllTimeStep(pressure_prev, pressure_curr, density, nx, ny, nz, fd_radius, coeff, dt, vel, nt, source, src_idx);
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
int sizeofnumbers=sizeofnx+sizeofny+sizeofnz+sizeofnt;
char *plotname = malloc(sizeof(char)*(strlen("soln_nx=_ny=_nz=_nt=") +sizeofnx+ strlen(".png")));
sprintf(plotname, "soln_nx=%d_ny=%d_nz=%d_nt=%d.png", nx, ny, nz, nt);
PlotSolution(pressure_curr, nx, ny, nz, src_z, resolution, plotname);
}


return 0;

}
