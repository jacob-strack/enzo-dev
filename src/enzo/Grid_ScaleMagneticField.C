//grid member function to rescale the magnetic field to a set strength, on average
//scale by same constant everywhere, don't want to introduce divergence

#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

int grid::ScaleMagneticField(float mean){
	if( ProcessorNumber != MyProcessorNumber)
		return SUCCESS; 
	
	int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum, B1Num, B2Num, B3Num;
	this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, 
					   TENum, B1Num, B2Num, B3Num);
	
	FLOAT B_sq = 0; //keep track of total magnetic field strength	
	int BNum[3] = {B1Num, B2Num, B3Num}; 
	int size = 1;
        	
	//loop over fields to get total B magnitude
	for(int dim = 0; dim < GridRank; dim++)
		size *= GridDimension[dim];
	
	for(int i = 0; i < size; i++)
		for(int dim = 0; dim < GridRank; dim++)
			B_sq += pow(BaryonField[BNum[dim]][i],2);

	//Get average magnetic field strength
	FLOAT B_st = pow(B_sq, 0.5) / size; 
	FLOAT scale = mean / B_st; 
	//Now loop again, using mean value to rescale
	for(int i = 0; i < size; i++){ 
		for(int dim = 0; dim < GridRank; dim++)
			BaryonField[BNum[dim]][i] *= scale; 
	}
	return SUCCESS;
}
