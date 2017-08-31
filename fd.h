#include <stdio.h>
#include <math.h>
#include <stdlib.h>


void fd_coeff(float *coeff, const float eval_point, const int order, float *points, const int num_points);
void Rickerwavelet(int nt, float dt, float freq, float *source, float velocity);
void Calcoperator(float *pressure, float *density, float *val, int ix, int nx, int iy, int ny, int iz, int nz, const int space_order, float *coeff);
void PlotSolution(float *solution, const int nx, const int ny, const int nz, const int z_coord, int resolution, const char *FILENAME);
void TimeStep(float *prev, float *curr, float *density, const int nx, const int ny, const int nz, const int fd_radius, float *coeff, const float dt, const float vel);
void AllTimeStep(float *prev, float *curr, float *density, const int nx, const int ny, const int nz, const int fd_radius, float *coeff, const float dt, const float vel, const int nt, float *source, const int src_idx);
void TimeStepPrecompute(float *prev, float *curr, float *density, float *dx, float *dy, float *dz, const int nx, const int ny, const int nz, const int fd_radius, float *coeff, const float dt, const float vel);
void AllTimeStepPreCompute(float *prev, float *curr, float *density, const int nx, const int ny, const int nz, const int fd_radius, float *coeff, const float dt, const float vel, const int nt, float *source, const int src_idx);
void TimeStepOpt(float *prev, float *curr, float *density, const int nx, const int ny, const int nz, const int fd_radius, float *coeff, const float dt, const float vel);
void AllTimeStepOpt(float *prev, float *curr, float *density, const int nx, const int ny, const int nz, const int fd_radius, float *coeff, const float dt, const float vel, const int nt, float *source, const int src_idx);
