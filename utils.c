#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
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
//printf("fd_radius=%d \n", fd_radius); 
//indicies start at 0 b/c coeff[0] = 0 in centered, first-order derivative approximation
float sol=0;
//for(int i=1; i<=1;++i){
//for(int j=1; j<=1;++j){
for(int i=1; i<=fd_radius;++i){
for(int j=1; j<=fd_radius;++j){

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
//+pressure[idx+(-i-j+1)*z_row]//*2.0/(density[idx+(-i+1)*z_row]+density[idx-i*z_row])
/*
*/
);

}}

//+pressure[idx+(-i-j+1)*z_row]//*2.0/(density[idx+(-i+1)*z_row]+density[idx-i*z_row])
*val=sol;
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
//    fprintf(temp, "%lf %lf %lf \n", xvals[i], yvals[j], PatchSolution(*mypatch, xvals[i], yvals[j])); //Write the data to a temporary file
    fprintf(temp, "%d %d %lf \n", i, j, solution[i+j*nx+z_coord*nx*ny]); //Write the data to a temporary file
    }
    
    
    for (int i=0; i < NUM_COMMANDS; i++)
    {
    fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]); //Send commands to gnuplot one by one.
    }
    fclose(temp);
    pclose(gnuplotPipe);

printf("Ending coeff generator \n");
}

//calculate deriv of input and store in dx array
//void Calcdxarray(float *input, float *dx)
//{

//}
