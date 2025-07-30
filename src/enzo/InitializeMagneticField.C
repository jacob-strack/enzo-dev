// Initialize Magnetic field by taking either curl of vector potential (electric field) or setting constant field


#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"
#include "phys_constants.h"

void InitializeMagneticField(LevelHierarchyEntry *LevelArray[], int TopologyMethod){
    if(TopologyMethod > 0){
	MHD_ProjectE=TRUE;
	MHD_ProjectB=FALSE;
	for(int mylevel = MaximumRefinementLevel; mylevel > 0; mylevel--){
		LevelHierarchyEntry *Temp = LevelArray[mylevel]; 
		while(Temp != NULL){
			if(Temp->GridData->ProjectSolutionToParentGrid(*LevelArray[mylevel-1]->GridData) == FAIL){
				fprintf(stderr, "Electric project failed. \n");
				//return FAIL;
			}	
			Temp = Temp->NextGridThisLevel;
		}	
	}
	for(int mylevel = 0; mylevel <= MaximumRefinementLevel; mylevel++){
		LevelHierarchyEntry *Temp = LevelArray[mylevel]; 
		while(Temp != NULL){
			int CurlStart[3] = {0,0,0};
			int CurlEnd[3];
			if(MyProcessorNumber == Temp->GridData->ReturnProcessorNumber()){
				for(int ind = 0; ind < 3; ind++)
					CurlEnd[ind] = Temp->GridData->GetGridDimension(ind) - 1;
				Temp->GridData->MHD_Curl(CurlStart, CurlEnd,0);
				//rescale magnetic field to given mean 
				FLOAT this_scale = 1e-16; //hard code shit for now, easy to generalize later
				Temp->GridData->ScaleMagneticField(this_scale); 
				Temp->GridData->CenterMagneticField();
			}
			Temp = Temp->NextGridThisLevel;
		}
	}
	}
	//Idea: Topology = 2, combine both methods. A constant background part and a Gaussian white noise potential part
	//At least it would take care of total field happening at some points. kind of a dumb way to fix a problem though
    
    else{
	    MHD_ProjectE=FALSE;
	    MHD_ProjectB=TRUE;
	    for (int level = MaximumRefinementLevel; level > 0; level--) {
	      LevelHierarchyEntry *Temp = LevelArray[level];
	      while (Temp != NULL) {
		if (Temp->GridData->ProjectSolutionToParentGrid(
					   *LevelArray[level-1]->GridData) == FAIL) {
		  fprintf(stderr, "Error in grid->ProjectSolutionToParentGrid.\n");
		  //return FAIL;
		}
		Temp = Temp->NextGridThisLevel;
	      }
	    }
    }
    MHD_ProjectE=TRUE;
    MHD_ProjectB=FALSE;
}

