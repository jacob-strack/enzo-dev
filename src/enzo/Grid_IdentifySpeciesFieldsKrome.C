/***********************************************************************
/
/  GRID CLASS (IDENTIFY THE MULTI-SPECIES FIELDS)
/
/  written by: KROME developers, based on the original ENZO routine
/  date:       2013
/
************************************************************************/
 
#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
 
/* function prototypes */
 
int FindField(int f, int farray[], int n);
 
 
int grid::IdentifySpeciesFieldsKrome(
 int &DeNum, int &HMNum, int &HINum, int &HeINum,
 int &H2INum, int &HIINum, int &HeIINum, int &H2IINum,
 int &HeIIINum
){

 DeNum = HMNum = HINum = HeINum = H2INum =
 HIINum = HeIINum = H2IINum = HeIIINum = 0;
 
  /* Find Fields for the KROME species */
 DeNum = FindField(ElectronDensity, FieldType, NumberOfBaryonFields);
 HMNum = FindField(HMDensity, FieldType, NumberOfBaryonFields);
 HINum = FindField(HIDensity, FieldType, NumberOfBaryonFields);
 HeINum = FindField(HeIDensity, FieldType, NumberOfBaryonFields);
 H2INum = FindField(H2IDensity, FieldType, NumberOfBaryonFields);
 HIINum = FindField(HIIDensity, FieldType, NumberOfBaryonFields);
 HeIINum = FindField(HeIIDensity, FieldType, NumberOfBaryonFields);
 H2IINum = FindField(H2IIDensity, FieldType, NumberOfBaryonFields);
 HeIIINum = FindField(HeIIIDensity, FieldType, NumberOfBaryonFields);


  /* Error if any not found. */
  if (DeNum<0 || HMNum<0 || HINum<0 || HeINum<0 ||
 H2INum<0 || HIINum<0 || HeIINum<0 || H2IINum<0 ||
 HeIIINum<0) {
    ENZO_VFAIL("De=%"ISYM", HM=%"ISYM", HI=%"ISYM", HeI=%"ISYM", H2I=%"ISYM", HII=%"ISYM", HeII=%"ISYM", H2II=%"ISYM", HeIII=%"ISYM"\n",
	    DeNum, HMNum, HINum, HeINum, H2INum, HIINum,
 HeIINum, H2IINum, HeIIINum)
  }
 
  return SUCCESS;
}
