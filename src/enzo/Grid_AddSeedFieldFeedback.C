
#ifdef USE_MPI
#include "mpi.h"
#endif
#include <math.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "list.h"

#define NO_DEATH 0
#define KILL_STAR 1
#define KILL_ALL 2

float incomplete_gamma(float r);
int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);

int grid::AddSeedFieldFeedback(float tau){
                            //Loop over SN list
                            //Check for Feedback
                            //Inject into cells up to cutoff radius  
    if(ProcessorNumber != MyProcessorNumber){
        return SUCCESS;
    }
    if(UseMagneticSupernovaFeedback && HydroMethod == 6){
        snsf_source_terms S;
        int igrid, index;
        std::vector<SuperNova>::iterator P = this->MagneticSupernovaList.begin(); 
        FLOAT cell_center[3]; 
        FLOAT dx, dy,dz, dist_to_sn; 
        //Loop over grids 
        static float DensityUnits, LengthUnits, TemperatureUnits = 1,
                     TimeUnits, VelocityUnits, MassUnits;
        if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
                     &TimeUnits, &VelocityUnits, Time) == FAIL) {
          ENZO_FAIL("Error in GetUnits.");}
        float MU = DensityUnits * pow(LengthUnits,3);
        for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
            for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
              for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
                igrid = i+(j+k*GridDimension[1])*GridDimension[0];
	            index = i+ ElectricDims[2][0] *(j+ ElectricDims[2][1]*k);
                std::vector<SuperNova>::iterator P = this->MagneticSupernovaList.begin();
                cell_center[0] = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
                cell_center[1] = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
                cell_center[2] = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
                while(P != this->MagneticSupernovaList.end()){
                    //printf("pos %.5e %.5e %.5e \n", P->getPosition()[0], P->getPosition()[1], P->getPosition()[2]);
                    //printf("cell %.5e %.5e %.5e \n", cell_center[0], cell_center[1], cell_center[2]); 
                    dx = (P->getPosition()[0] - cell_center[0]);
                    dy = (P->getPosition()[1] - cell_center[1]);
                    dz = (P->getPosition()[2] - cell_center[2]);
                    dist_to_sn = sqrt(dx*dx + dy*dy + dz*dz);
                    float L = CellWidth[0][i];
                    //printf("dx %.5e dy %.5e dz %.5e \n", dx, dy, dz);
                    //printf("width %.5e \n", L);
                    tau = .2;
                    float sn_time_active = this->ReturnTime() - P->getTime();
                    //printf("ReturnTime: %.5e pGetTime %.5e \n", this->ReturnTime(), P->getTime());
                    //Find electric field z component, add time derivative of vector potential
                    printf("dist: %.5e, 5L: %.5e, 5*tau %.5e , sn_time_active %.5e time %.5e ratio %.5e \n", dist_to_sn, 5*L, 5*tau,sn_time_active, this->ReturnTime(), sn_time_active/(0.5*tau));
                    if(dist_to_sn < 3*L && sn_time_active < 5*tau && sn_time_active > 0){
                        printf("dist: %.5e, 3L: %.5e, 5*tau %.5e , sn_time_active %.5e time %.5e ratio %.5e \n", dist_to_sn, 3*L, 5*tau,sn_time_active, this->ReturnTime(), sn_time_active/(0.5*tau));
                        float B_0 = (4 /pow(3.141*pow(L,3), 1/2)) * sqrt((3.3e48) / (MU *  pow(LengthUnits,2)) * pow(TimeUnits,2));
                        printf("B_0: %.5e \n", B_0);
                        printf("1st: %.5e 2nd: %.5e \n", (4 /pow(3.141*pow(L,3), 1/2)),sqrt((3.3e48) / (MU *  pow(LengthUnits,2)) * pow(TimeUnits,2)));
                        printf("TimeUnits %.5e LengthUnits %.5e MassUnits %.5e DensityUnits %.5e \n", TimeUnits, LengthUnits, MU, DensityUnits);
                        printf("indices: %d %d %d \n", i,j,k);
                        //printf("ElectricField Added %.5e \n",B_0*pow(2,-1/4)*(L/tau)*exp(-dz*dz / (2*L*L))*exp(-(sn_time_active)/tau)*incomplete_gamma(sqrt(dx*dx + dy*dy)));
                        ElectricField[2][index] -= B_0*pow(2,-1/4)*(L/tau)*exp(-sn_time_active/tau)*exp(-dz*dz / (2*L*L))*incomplete_gamma(sqrt(dx*dx + dy*dy));//B_0*pow(2,-1/4)*(L/tau)*exp(-dz*dz / (2*L*L))*exp(-(sn_time_active)/tau)*incomplete_gamma(sqrt(dx*dx + dy*dy)); 
                    }
                      
                    P++;
                }
              } // End of k for-loop                                                                                                         
            } // End of j for-loop                                                                                                           
        } // End of i for-loop            
    }
    return SUCCESS;
}
float incomplete_gamma(float r){//function that computes the lower incomplete gamma function from linear interpolation
    float x_pts[21] = {0.0, 5e-6,1e-5,2.5e-5,5e-5,7.5e-5, 1e-4, 1.5e-4, 2e-4, 2.5e-4, 3e-4, 3.5e-4, 4e-4, 4.5e-4, 5e-4, 5.5e-4, 6e-4, 6.5e-4, 7e-4, 7.5e-4, 8e-4}; 
    float y_pts[21] = {0.0, 0.000140983, 0.000237103, 0.000471399, 0.00079279,0.00107454, 0.00133328, 0.00180709, 0.0022422, 0.00265062, 0.00303895,0.00341134, 0.00377059, 0.00411874, 0.00445731, 0.00478750, 0.00511023, 0.00542629, 0.00573630, 0.00604081,0.00634026};
    //float x_pts[10] = {.01, .02, .03, .04, .05, .06, .07, .08, .09, .1};
    //float y_pts[10] = {.04198, .07031, .09498, .1172, .1380, .1576, .1761, .1939,.2109,.2273};

    float y_interp;
    if(r > 0.0008)
        y_interp = 0.0;
    if(r < 0.0008){ 
    int index_l=0;
    int index_u = 0;
    for(int i = 0; i < 10; i++){
        if(x_pts[i] > x_pts[index_l] && x_pts[i] < r)
            index_l = i; 
        if(x_pts[10 - i] < x_pts[10 - index_u] && x_pts[10 - i] > r)
            index_u = 10 - i; 
    }
    y_interp =((y_pts[index_u] - y_pts[index_l]) / (x_pts[index_u] - x_pts[index_l])) * (r - x_pts[index_l]) + y_pts[index_l];
    printf("x_l %.5e x_h %.5e \n", x_pts[index_l], x_pts[index_u]);
    }
    printf("s: %.5e incgamma: %.5e \n", r,y_interp);
    return y_interp;}

