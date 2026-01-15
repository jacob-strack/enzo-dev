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
 int &DeNum, int &HMNum, int &CMNum, int &OMNum,
 int &HINum, int &HeINum, int &H2INum, int &CINum,
 int &OINum, int &OHINum, int &COINum, int &CHINum,
 int &CH2INum, int &C2INum, int &HCOINum,
 int &H2OINum, int &O2INum, int &CO_TOTALINum,
 int &H2O_TOTALINum, int &HIINum, int &HeIINum,
 int &H2IINum, int &CIINum, int &OIINum, int &HOCIINum,
 int &HCOIINum, int &H3IINum, int &CHIINum,
 int &CH2IINum, int &COIINum, int &CH3IINum,
 int &OHIINum, int &H2OIINum, int &H3OIINum,
 int &O2IINum, int &HeIIINum
){

 DeNum = HMNum = CMNum = OMNum = HINum = HeINum =
 H2INum = CINum = OINum = OHINum = COINum =
 CHINum = CH2INum = C2INum = HCOINum = H2OINum =
 O2INum = CO_TOTALINum = H2O_TOTALINum = HIINum =
 HeIINum = H2IINum = CIINum = OIINum = HOCIINum =
 HCOIINum = H3IINum = CHIINum = CH2IINum =
 COIINum = CH3IINum = OHIINum = H2OIINum =
 H3OIINum = O2IINum = HeIIINum = 0;
 
  /* Find Fields for the KROME species */
 DeNum = FindField(ElectronDensity, FieldType, NumberOfBaryonFields);
 HMNum = FindField(HMDensity, FieldType, NumberOfBaryonFields);
 CMNum = FindField(CMDensity, FieldType, NumberOfBaryonFields);
 OMNum = FindField(OMDensity, FieldType, NumberOfBaryonFields);
 HINum = FindField(HIDensity, FieldType, NumberOfBaryonFields);
 HeINum = FindField(HeIDensity, FieldType, NumberOfBaryonFields);
 H2INum = FindField(H2IDensity, FieldType, NumberOfBaryonFields);
 CINum = FindField(CIDensity, FieldType, NumberOfBaryonFields);
 OINum = FindField(OIDensity, FieldType, NumberOfBaryonFields);
 OHINum = FindField(OHIDensity, FieldType, NumberOfBaryonFields);
 COINum = FindField(COIDensity, FieldType, NumberOfBaryonFields);
 CHINum = FindField(CHIDensity, FieldType, NumberOfBaryonFields);
 CH2INum = FindField(CH2IDensity, FieldType, NumberOfBaryonFields);
 C2INum = FindField(C2IDensity, FieldType, NumberOfBaryonFields);
 HCOINum = FindField(HCOIDensity, FieldType, NumberOfBaryonFields);
 H2OINum = FindField(H2OIDensity, FieldType, NumberOfBaryonFields);
 O2INum = FindField(O2IDensity, FieldType, NumberOfBaryonFields);
 CO_TOTALINum = FindField(CO_TOTALIDensity, FieldType, NumberOfBaryonFields);
 H2O_TOTALINum = FindField(H2O_TOTALIDensity, FieldType, NumberOfBaryonFields);
 HIINum = FindField(HIIDensity, FieldType, NumberOfBaryonFields);
 HeIINum = FindField(HeIIDensity, FieldType, NumberOfBaryonFields);
 H2IINum = FindField(H2IIDensity, FieldType, NumberOfBaryonFields);
 CIINum = FindField(CIIDensity, FieldType, NumberOfBaryonFields);
 OIINum = FindField(OIIDensity, FieldType, NumberOfBaryonFields);
 HOCIINum = FindField(HOCIIDensity, FieldType, NumberOfBaryonFields);
 HCOIINum = FindField(HCOIIDensity, FieldType, NumberOfBaryonFields);
 H3IINum = FindField(H3IIDensity, FieldType, NumberOfBaryonFields);
 CHIINum = FindField(CHIIDensity, FieldType, NumberOfBaryonFields);
 CH2IINum = FindField(CH2IIDensity, FieldType, NumberOfBaryonFields);
 COIINum = FindField(COIIDensity, FieldType, NumberOfBaryonFields);
 CH3IINum = FindField(CH3IIDensity, FieldType, NumberOfBaryonFields);
 OHIINum = FindField(OHIIDensity, FieldType, NumberOfBaryonFields);
 H2OIINum = FindField(H2OIIDensity, FieldType, NumberOfBaryonFields);
 H3OIINum = FindField(H3OIIDensity, FieldType, NumberOfBaryonFields);
 O2IINum = FindField(O2IIDensity, FieldType, NumberOfBaryonFields);
 HeIIINum = FindField(HeIIIDensity, FieldType, NumberOfBaryonFields);


  /* Error if any not found. */
  if (DeNum<0 || HMNum<0 || CMNum<0 || OMNum<0 ||
 HINum<0 || HeINum<0 || H2INum<0 || CINum<0 ||
 OINum<0 || OHINum<0 || COINum<0 || CHINum<0 ||
 CH2INum<0 || C2INum<0 || HCOINum<0 || H2OINum<0 ||
 O2INum<0 || CO_TOTALINum<0 || H2O_TOTALINum<0 ||
 HIINum<0 || HeIINum<0 || H2IINum<0 || CIINum<0 ||
 OIINum<0 || HOCIINum<0 || HCOIINum<0 || H3IINum<0 ||
 CHIINum<0 || CH2IINum<0 || COIINum<0 || CH3IINum<0 ||
 OHIINum<0 || H2OIINum<0 || H3OIINum<0 ||
 O2IINum<0 || HeIIINum<0) {
    ENZO_VFAIL("De=%"ISYM", HM=%"ISYM", CM=%"ISYM", OM=%"ISYM", HI=%"ISYM", HeI=%"ISYM", H2I=%"ISYM", CI=%"ISYM", OI=%"ISYM", OHI=%"ISYM", COI=%"ISYM", CHI=%"ISYM", CH2I=%"ISYM", C2I=%"ISYM", HCOI=%"ISYM", H2OI=%"ISYM", O2I=%"ISYM", CO_TOTALI=%"ISYM", H2O_TOTALI=%"ISYM", HII=%"ISYM", HeII=%"ISYM", H2II=%"ISYM", CII=%"ISYM", OII=%"ISYM", HOCII=%"ISYM", HCOII=%"ISYM", H3II=%"ISYM", CHII=%"ISYM", CH2II=%"ISYM", COII=%"ISYM", CH3II=%"ISYM", OHII=%"ISYM", H2OII=%"ISYM", H3OII=%"ISYM", O2II=%"ISYM", HeIII=%"ISYM"\n",
	    DeNum, HMNum, CMNum, OMNum, HINum, HeINum,
 H2INum, CINum, OINum, OHINum, COINum, CHINum,
 CH2INum, C2INum, HCOINum, H2OINum, O2INum,
 CO_TOTALINum, H2O_TOTALINum, HIINum, HeIINum,
 H2IINum, CIINum, OIINum, HOCIINum, HCOIINum,
 H3IINum, CHIINum, CH2IINum, COIINum, CH3IINum,
 OHIINum, H2OIINum, H3OIINum, O2IINum, HeIIINum)
  }
 
  return SUCCESS;
}
