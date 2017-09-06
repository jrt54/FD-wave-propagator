#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "fd.h"
//#include <ctime.h>
#include <sys/time.h>
#define MIN(x, y) (((x) < (y)) ? (x): (y))
#define MAX(x, y) (((x) > (y)) ? (x): (y))


void fd_coeff(float *coeff, const float eval_point, const int order, float *points, const int num_points)
{
    float c1, c2, c3;
    float x_0=eval_point;
    float center=0;

    


//  float* d = (float*) malloc((order+1)*num_points*num_points*sizeof(float));
    float d[(order+1)*num_points*num_points];
    int m_idx = (order+1)*num_points;
    int n_idx = num_points;
    
    //array initializer 1
    /*
    memset(d, 0.f, sizeof(d));
    */

    //array initializer 2
    int sizeofd = (order+1)*(num_points)*(num_points)*sizeof(float);
    memset(d, 0.f, sizeofd);
    

    //array initializer 3
    /*
    for(int m=0; m <= order; ++m){
	for(int n=0; n< num_points; ++n){
	    for(int v=0; v<num_points;++v){
	    d[m*m_idx+n*n_idx+v]=0.f;
	    }}}
    */
    

    d[0]=1.f;
    c1 = 1.f;

    for(int n=1; n<=num_points-1;++n){
        c2=1.f;
	for(int v=0; v<=n-1; ++v){
            c3 = points[n] - points[v];
            c2 = c2*c3;
            for(int m=0; m<=MIN(n, order); ++m){
		d[m*m_idx+n*n_idx + v] = (points[n]-x_0)*d[m*m_idx + (n-1)*n_idx + v] - m*d[(m-1)*m_idx + (n-1)*n_idx + v];
		d[m*m_idx + n*n_idx + v] *= 1.f/c3;
            }
	}
	for(int m=0; m<= MIN(n, order); ++m){
            d[m*m_idx+n*n_idx+n] = m*d[(m-1)*m_idx+(n-1)*n_idx+(n-1)] - (points[n-1]-x_0)*d[m*m_idx+(n-1)*n_idx+n-1];
            d[m*m_idx+n*n_idx+n] *= c1/c2;
	}
        c1=c2;
    }

    for(int i=0; i<num_points; ++i){
	coeff[i] = d[order*m_idx+(num_points-1)*n_idx + i];
    }

//    free(d);
}

//#include "fd.h"
void Rickerwavelet(int nt, float dt, float freq, float *source, float velocity)
{
float velocity_sq = velocity*velocity;
float pi_freq_sq = M_PI*M_PI*freq*freq; 
float pi_freq = M_PI*freq; 
float shift = 1.0/freq;
for(int i=0; i<nt; ++i){
//float r = pi_freq*(i*dt-shift); 
//source[i] = velocity_sq*(1-2*r*r)*exp(-r*r);


float curr_t = -shift + i*dt;
float curr_t_sq=curr_t*curr_t;
source[i] = (1-2*pi_freq_sq*curr_t_sq)*exp(-pi_freq_sq*curr_t_sq);
//printf("source[%d]: %f \n",i, source[i]);
}

}

//double for loop calculation
void Calcoperator(float *pressure, float *density, float *val, int ix, int nx, int iy, int ny, int iz, int nz, const int fd_radius, float *coeff)
{
int z_row=nx*ny;//row-major order adjustment
int idx=ix + nx*iy + z_row*iz; //current row-major index
//indicies start at 0 b/c coeff[0] = 0 in centered, first-order derivative approximation
float sol=0;
for(int i=1; i<=fd_radius;++i){
for(int j=1; j<=fd_radius;++j){
//printf("i, j = %d %d \n", i, j);
sol+=coeff[i]*coeff[j]*(

//below is the dx compoent 
pressure[idx+i+j-1]*2.0/(density[idx+i-1]+density[idx+i]) //p_++/q(.5+)
-pressure[idx-i+j]*2.0/(density[idx-i+1]+density[idx-i])//p_-+/q(.5-)
-pressure[idx+i-j]*2.0/(density[idx+i-1]+density[idx+i]) //p_++/q(.5+)
+pressure[idx-i-j+1]*2.0/(density[idx-i+1]+density[idx-i])//p_-+/q(.5-)

//dy component
+pressure[idx+(i+j-1)*nx]*2.0/(density[idx+(i-1)*nx]+density[idx+i*nx]) 
-pressure[idx+(-i+j)*nx]*2.0/(density[idx+(-i+1)*nx]+density[idx-i*nx])
-pressure[idx+(i-j)*nx]*2.0/(density[idx+(i-1)*nx]+density[idx+i*nx]) 
+pressure[idx+(-i-j+1)*nx]*2.0/(density[idx+(-i+1)*nx]+density[idx-i*nx])

//dz
+pressure[idx+(i+j-1)*z_row]*2.0/(density[idx+(i-1)*z_row]+density[idx+i*z_row]) 
-pressure[idx+(-i+j)*z_row]*2.0/(density[idx+(-i+1)*z_row]+density[idx-i*z_row])
-pressure[idx+(i-j)*z_row]*2.0/(density[idx+(i-1)*z_row]+density[idx+i*z_row]) 
+pressure[idx+(-i-j+1)*z_row]*2.0/(density[idx+(-i+1)*z_row]+density[idx-i*z_row])
);

}}

*val=sol;
}


void TimeStep(float *prev, float *curr, float *density, const int nx, const int ny, const int nz, const int fd_radius, float *coeff, const float dt, const float vel)
{
const int z_row=ny*nx;
float scalar = vel*vel*dt*dt;



const int offset = 2*fd_radius-1;
#pragma omp parallel for collapse(3) 
for(int i=offset; i<nx-offset; ++i){
for(int j=offset; j<ny-offset; ++j){
for(int k=offset; k<nz-offset; ++k){
float operator;
int idx = i+j*nx+k*z_row;
//printf("x, y, z = %d %d %d \n", i, j, k);
Calcoperator(curr, density, &operator, i,nx, j, ny, k, nz, fd_radius, coeff);
//printf("Operator calced = %f  \n", operator);
prev[i+j*nx+k*z_row]=2*curr[i+j*nx+k*z_row]-prev[i+j*nx+k*z_row]+density[idx]*operator*scalar;


//prev[i+j*nx+k*z_row]=2*curr[i+j*nx+k*z_row]-prev[i+j*nx+k*z_row];
}}}



}


void AllTimeStep(float *prev, float *curr, float *density, const int nx, const int ny, const int nz, const int fd_radius, float *coeff, const float dt, const float vel, const int nt, float *source, const int src_idx)
{
const float scalar = dt*dt*vel*vel;
for(int timestep=0; timestep<nt; timestep+=2)
{
TimeStep(prev, curr, density, nx, ny, nz, fd_radius, coeff, dt, vel);
prev[src_idx] += source[timestep]*scalar;
//prev[src_idx] += source[timestep];

TimeStep(curr, prev, density, nx, ny, nz, fd_radius, coeff, dt, vel);
curr[src_idx] += source[timestep+1]*scalar;
//prev[src_idx] += source[timestep+1];

}

}


void TimeStepPrecompute(float *prev, float *curr, float *density, float *dx, float *dy, float *dz, const int nx, const int ny, const int nz, const int fd_radius, float *coeff, const float dt, const float vel, uint64* log)
{
const int z_row=ny*nx;

float scalar = vel*vel*dt*dt;

#pragma omp parallel for collapse(3) 
for(int i=fd_radius-1; i<nx-fd_radius-1; ++i){
for(int j=fd_radius-1; j<ny-fd_radius-1; ++j){
for(int k=fd_radius-1; k<nz-fd_radius-1; ++k){
int idx = i+j*nx+k*z_row;
float der_x=0;
float der_y=0;
float der_z=0;

for(int r=1; r<=fd_radius;++r){
//all these are implicitly referring to idx+.5
der_x += coeff[r]*(curr[idx+r]-curr[idx-(r-1)]);
der_y += coeff[r]*(curr[idx+r*nx]-curr[idx-(r-1)*nx]);
der_z += coeff[r]*(curr[idx+r*z_row]-curr[idx-(r-1)*z_row]);
}

//all these are implicitly referring to idx+.5
dx[idx]=der_x;
dy[idx]=der_y;
dz[idx]=der_z;

}}}


const int offset = 2*fd_radius-1;
#pragma omp parallel for collapse(3) 
for(int i=offset; i<nx-offset; ++i){
for(int j=offset; j<ny-offset; ++j){
for(int k=offset; k<nz-offset; ++k){
int idx = i+j*nx+k*z_row;
float operator=0;
for(int r=1; r<=fd_radius;++r){

//dx component
operator+=coeff[r]*dx[idx+(r-1)]*2.0/(density[idx+r]+density[idx+r-1]);
operator-=coeff[r]*dx[idx-r]*2.0/(density[idx-r] + density[idx-r+1]);

//dy component
operator+=coeff[r]*dy[idx+(r-1)*nx]*2.0/(density[idx+r*nx]+density[idx+(r-1)*nx]);
operator-=coeff[r]*dy[idx-r*nx]*2.0/(density[idx-r*nx] + density[idx-(r+1)*nx]);

//dz component
operator+=coeff[r]*dz[idx+(r-1)*z_row]*2.0/(density[idx+r*z_row]+density[idx+(r-1)*z_row]);
operator-=coeff[r]*dz[idx-r*z_row]*2.0/(density[idx-r*z_row] + density[idx-(r+1)*z_row]);
}

uint64 update_start=GetTimeMs64();
prev[idx]=2*curr[idx]-prev[idx]+density[idx]*operator*scalar;
*log+=((GetTimeMs64()-update_start)); //in seconds

}}}


}

void AllTimeStepPreCompute(float *prev, float *curr, float *density, const int nx, const int ny, const int nz, const int fd_radius, float *coeff, const float dt, const float vel, const int nt, float *source, const int src_idx, uint64* log)
{
const float scalar = dt*dt*vel*vel;

const int total_size = nx*ny*nz;
float *dx = (float*) malloc(sizeof(float)*nx*ny*nz);
dx = memset(dx, 0, total_size*sizeof(dx[0]));
float *dy = (float*) malloc(sizeof(float)*nx*ny*nz);
dy = memset(dy, 0, total_size*sizeof(dy[0]));
float *dz = (float*) malloc(sizeof(float)*nx*ny*nz);
dz = memset(dz, 0, total_size*sizeof(dz[0]));



for(int timestep=0; timestep<nt; timestep+=2)
{
TimeStepPrecompute(prev, curr, density, dx, dy, dz, nx, ny, nz, fd_radius, coeff, dt, vel, log);
prev[src_idx] += source[timestep]*scalar;
//prev[src_idx] += source[timestep];

TimeStepPrecompute(curr, prev, density, dx, dy, dz, nx, ny, nz, fd_radius, coeff, dt, vel, log);
curr[src_idx] += source[timestep+1]*scalar;
//prev[src_idx] += source[timestep+1];

}
}

unsigned long long GetTimeMs64()
{
#ifdef _WIN32
 /* Windows */
 FILETIME ft;
 LARGE_INTEGER li;

 /* Get the amount of 100 nano seconds intervals elapsed since January 1, 1601 (UTC) and copy it
  * to a LARGE_INTEGER structure. */
 GetSystemTimeAsFileTime(&ft);
 li.LowPart = ft.dwLowDateTime;
 li.HighPart = ft.dwHighDateTime;

 uint64 ret = li.QuadPart;
 ret -= 116444736000000000LL; /* Convert from file time to UNIX epoch time. */
 ret /= 10000; /* From 100 nano seconds (10^-7) to 1 millisecond (10^-3) intervals */

 return ret;
#else
 /* Linux */
 struct timeval tv;

 gettimeofday(&tv, NULL);

 unsigned long long ret = tv.tv_usec;
 /* Convert from micro seconds (10^-6) to milliseconds (10^-3) */
 ret /= 1000;

 /* Adds the seconds (10^0) after converting them to milliseconds (10^-3) */
 ret += (tv.tv_sec * 1000);

 return ret;
#endif
}

void TimeStepOpt(float *prev, float *curr, float *density, float *grad, const int nx, const int ny, const int nz, const int fd_radius, const float *coeff, const float dt, const float vel, unsigned long long *log)
{
const int z_row=ny*nx;
const int offset = 2*fd_radius-1;

const float scalar = vel*vel*dt*dt;


grad = memset(grad, 0, nx*ny*nz*sizeof(grad[0]));


//#pragma omp parallel for reduction(+:grad) 
#pragma omp parallel for collapse(3)// reduction(+:grad) 
for(int i=fd_radius-1; i<nx-fd_radius-1; ++i){
for(int j=fd_radius-1; j<ny-fd_radius-1; ++j){
for(int k=fd_radius-1; k<nz-fd_radius-1; ++k){
int idx = i+j*nx+k*z_row;
float der_x=0;
float der_y=0;
float der_z=0;

for(int r=1; r<=fd_radius;++r){
//all these are implicitly referring to idx+.5
der_x += coeff[r]*(curr[idx+r]-curr[idx-(r-1)]);
der_y += coeff[r]*(curr[idx+r*nx]-curr[idx-(r-1)*nx]);
der_z += coeff[r]*(curr[idx+r*z_row]-curr[idx-(r-1)*z_row]);
}

der_x *= 2.0/(density[idx+1]+density[idx]);
der_y *= 2.0/(density[idx+nx]+density[idx]);
der_z *= 2.0/(density[idx+z_row]+density[idx]);

//prev[idx]+=coeff[1]*der_x*2.0/(density[idx+1]+density[idx]); //can be done with first iterate of below for loop

unsigned long long update_start=GetTimeMs64();
	for(int r=0; r<=fd_radius-1;++r){
	//all these are implicitly referring to idx+.5

	#pragma omp atomic update
	grad[idx-r]+=coeff[r+1]*der_x;
	#pragma omp atomic update
	grad[idx-r*nx]+=coeff[r+1]*der_y;
	#pragma omp atomic update
	grad[idx-r*z_row]+=coeff[r+1]*der_z;

	#pragma omp atomic update
	grad[idx+(r+1)]-=coeff[r+1]*der_x;
	#pragma omp atomic update
	grad[idx+(r+1)*nx]-=coeff[r+1]*der_y;
	#pragma omp atomic update
	grad[idx+(r+1)*z_row]-=coeff[r+1]*der_z;
	}

//*log+=((double)(CLOCKS_PER_SEC))/CLOCKS_PER_SEC; //in seconds
//#pragma omp atomic update
*log+=((GetTimeMs64()-update_start)); //in seconds
}}}


#pragma omp parallel for collapse(3) 
for(int i=fd_radius-1; i<nx-fd_radius-1; ++i){
for(int j=fd_radius-1; j<ny-fd_radius-1; ++j){
for(int k=fd_radius-1; k<nz-fd_radius-1; ++k){
unsigned long long update_start=GetTimeMs64();
int idx = i+j*nx+k*z_row;
prev[idx]=2*curr[idx]-prev[idx]+scalar*density[idx]*grad[idx];
*log+=((GetTimeMs64()-update_start)); //in seconds
}}}



}


void AllTimeStepOpt(float *prev, float *curr, float *density, const int nx, const int ny, const int nz, const int fd_radius, float *coeff, const float dt, const float vel, const int nt, float *source, const int src_idx, unsigned long long *log)
{
float *grad = (float*) malloc(sizeof(float)*nx*ny*nz);
const float scalar = dt*dt*vel*vel;
for(int timestep=0; timestep<nt; timestep+=2)
{
TimeStepOpt(prev, curr, density, grad, nx, ny, nz, fd_radius, coeff, dt, vel, log);
prev[src_idx] += source[timestep]*scalar;
//prev[src_idx] += source[timestep];

TimeStepOpt(curr, prev, density, grad, nx, ny, nz, fd_radius, coeff, dt, vel, log);
curr[src_idx] += source[timestep+1]*scalar;
//prev[src_idx] += source[timestep+1];

}

}

//plot solution at z_slice
void PlotSolution(float *solution, const int nx, const int ny, const int nz, const int z_coord, int resolution, const char *FILENAME)
{
    int NUM_POINTS=resolution;
    
    //PetscScalar LLxcoord = mypatch->LLxcoord;    
    //PetscScalar LLycoord = mypatch->LLycoord;    
    //PetscScalar URxcoord = mypatch->URxcoord;    
    //PetscScalar URycoord = mypatch->URycoord;     
    int NUM_COMMANDS=5;
    //char * commandsForGnuplot[] = {"set title \"Errorplot\"", "plot 'data.temp'"};
    char * commandsForGnuplot[] = {
					"set terminal png",
					"set pm3d interpolate 30,30",
					"set view map",
					"set title 'Solutionplot'", 
    					"splot 'data.temp' using 1:2:3 with image \n"
				};
    /*double xvals[nx];
    for(int i=0; i< nx; ++i)
    {
    xvals[i] =  i 
    }	
    
    double yvals[ny];
    for(int i=0; i< NUM_POINTS; ++i)
    {
    yvals[i] = i;  
    }*/	

    FILE * temp = fopen("data.temp", "w");
    //Opens an interface that one can use to send commands as if they were typing into the
    //     gnuplot command line.  "The -persistent" keeps the plot open even after your
    //     C program terminates.
    //
    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
    //FILE * gnuplotPipe = fopen ("gnuplot -persistent", "w");
    fprintf(gnuplotPipe, "set output \"%s\"\n", FILENAME);  
    for (int i=0; i < nx; ++i)
    for (int j=0; j < ny; ++j)
    {
    fprintf(temp, "%d %d %lf \n", i, j, solution[i+j*nx+z_coord*nx*ny]); //Write the data to a temporary file
    //fprintf(temp, "%d %d %lf \n", i, j, solution[z_coord+j*nx+i*nx*ny]); //x slice
    }
    
    
    for (int i=0; i < NUM_COMMANDS; i++)
    {
    fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]); //Send commands to gnuplot one by one.
    }
    fclose(temp);
    pclose(gnuplotPipe);

}
//calculate deriv of input and store in dx array
//void Calcdxarray(float *input, float *dx)
//{

//}
