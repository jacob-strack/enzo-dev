//Star maker routine that places a star particle at given position at current time. No checks before forming. No mass removed from grid.
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "phys_constants.h"

int star_maker_force(int *nx, int *ny, int *nz, int *size, float *u, float *v, float *w, float *dt, FLOAT *t, FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,float *mp, int *ctype, float *tcp, float *tdp, int *nstar){
    int i, j, k, index; //defining indices for loops 
    printf("in star_maker_now");
    if(*t < ForceStarTime && *t + *dt > ForceStarTime){
    printf("in sf loop \n");
    *xp = ForceStarPosition[0]; 
    *yp = ForceStarPosition[1]; 
    *zp = ForceStarPosition[2]; 
    //find velocity at position to place particle 
    *up = 0; 
    *vp = 0; 
    *wp = 0; 
    *mp = 1.0; 
    *ctype = -2; 
    //setting creation time 
    *tcp = (float) *t; 
    *tdp = 1.0 + *tcp; //want particle to die fast 
    *nstar += 1;
    printf("nstar: %d \n",*nstar);
    }
    return SUCCESS;
}


