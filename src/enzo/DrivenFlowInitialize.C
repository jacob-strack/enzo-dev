/***********************************************************************
/
/  INITIALIZE DRIVEN FLOW SIMULATION
/
/  written by: Wolfram Schmidt
/  date:       May, 2005
/  modified1:  April, 2007
/  modified2: Sep, 2014: updated to support Enzo 2.4   // P. Grete
/
/  PURPOSE: Initializes simulation of flow driven by stochastic forcing
/
/  RETURNS: SUCCESS or FAIL
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
#include "Hierarchy.h"
#include "TopGridData.h"
#include "phys_constants.h"

StochasticForcing Forcing;
void MHDCTSetupFieldLabels();
int ReadEquilibriumTable(char *name, FLOAT Time);
int DrivenFlowInitialize(FILE *fptr, FILE *Outfptr, 
             HierarchyEntry &TopGrid, TopGridData &MetaData, 
             int SetBaryonFields)
{
  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *ElectronName = "Electron_Density";
  char *HIName    = "HI_Density";
  char *HIIName   = "HII_Density";
  char *HeIName   = "HeI_Density";
  char *HeIIName  = "HeII_Density";
  char *HeIIIName = "HeIII_Density";
  char *HMName    = "HM_Density";
  char *H2IName   = "H2I_Density";
  char *H2IIName  = "H2II_Density";
  char *DIName    = "DI_Density";
  char *DIIName   = "DII_Density";
  char *HDIName   = "HDI_Density";
  char *MetalName   = "Metal_Density";
  char *MetalIaName = "MetalSNIa_Density";
  char *CMName = "CMDensity"; 
  char *OMName = "OMDensity"; 
  char *CIName = "CIDensity"; 
  char *OIName = "OIDensity"; 
  char *OHIName = "OHIDensity"; 
  char *COIName = "COIDensity"; 
  char *CHIName = "CHIDensity"; 
  char *CH2IName = "CH2IDensity"; 
  char *C2IName = "C2IDensity"; 
  char *HCOIName = "HCOIDensity"; 
  char *H2OIName = "H2O0IDensity"; 
  char *O2IName = "O2IDensity"; 
  char *CO_TOTALIName = "CO_TOTALIDensity"; 
  char *H2O_TOTALIName = "H2O_TOTALIDensity"; 
  char *CIIName = "CIIDensity"; 
  char *OIIName = "OIIDensity"; 
  char *HOCIIName = "HOCIIDensity";
  char *HCOIIName = "HCOIIDensity"; 
  char *H3IIName = "H3IIDensity"; 
  char *CHIIName = "CHIIDensity"; 
  char *CH2IIName = "CH2IIDensity"; 
  char *COIIName = "COIIDensity"; 
  char *CH3IIName = "CH3IIDensity"; 
  char *OHIIName = "OHIIDensity"; 
  char *H2OIIName = "H2OIIDensity"; 
  char *H3OIIName = "H3OIIDensity"; 
  char *O2IIName = "O2IIDensity"; 
  char *BxName = "Bx";
  char *ByName = "By";
  char *BzName = "Bz";
  char *PhiName = "Phi";
  char *StochAccel1Name = "x-acceleration";
  char *StochAccel2Name = "y-acceleration";
  char *StochAccel3Name = "z-acceleration";

  /* declarations */

  char line[MAX_LINE_LENGTH];
  int ret;

  /* set default parameters specifying the random force field */

  int DrivenFlowAlpha[3]       = {1, 1, 1};       // ratio of domain size to characteristic length
  int DrivenFlowSeed = 20150418;                  // seed of random number generator
  float DrivenFlowBandWidth[3] = {1.0, 1.0, 1.0}; // band width (1.0 = maximal)
  float DrivenFlowAutoCorrl[3] = {1.0, 1.0, 1.0}; // ratio auto-correlation to large-eddy turn-over time scale
  float DrivenFlowMach[3]      = {1.0, 1.0, 1.0}; // Mach number
  float DrivenFlowWeight       = 1.0;             // weight of solenoidal components

  /* set other default parameters */

  float DrivenFlowDensity     = 1.0; // initial mass density
  float DrivenFlowPressure    = 1.0; // initial pressure

  float DrivenFlowMagField           = 0.0; // initial magnetic field

  forcing_type DrivenFlowProfile;            //defined in typedefs.h
  float DrivenFlowDomainLength[3];
  float DrivenFlowVelocity[3];
  float SoundSpeed;

  /* read input from file */
  if (debug) printf("DrivenFlowInitialize: reading problem-specific parameters.\n");

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "DrivenFlowProfile = %"ISYM, &DrivenFlowProfile);
    
    ret += sscanf(line, "DrivenFlowAlpha = %"ISYM" %"ISYM" %"ISYM, 
          DrivenFlowAlpha, DrivenFlowAlpha+1, DrivenFlowAlpha+2);
    ret += sscanf(line, "DrivenFlowSeed = %"ISYM, &DrivenFlowSeed);

    ret += sscanf(line, "DrivenFlowBandWidth = %"FSYM"%"FSYM"%"FSYM, 
          DrivenFlowBandWidth, DrivenFlowBandWidth+1, DrivenFlowBandWidth+2);

    ret += sscanf(line, "DrivenFlowAutoCorrl = %"FSYM"%"FSYM"%"FSYM, 
          DrivenFlowAutoCorrl, DrivenFlowAutoCorrl+1, DrivenFlowAutoCorrl+2);

    ret += sscanf(line, "DrivenFlowMach = %"FSYM"%"FSYM"%"FSYM, 
          DrivenFlowMach, DrivenFlowMach+1, DrivenFlowMach+2);

    ret += sscanf(line, "DrivenFlowWeight = %"FSYM, &DrivenFlowWeight);

    ret += sscanf(line, "DrivenFlowDensity = %"FSYM, &DrivenFlowDensity);
    ret += sscanf(line, "DrivenFlowPressure = %"FSYM, &DrivenFlowPressure);
    ret += sscanf(line, "DrivenFlowMagField = %"FSYM, &DrivenFlowMagField);
    
    ret += sscanf(line, "use_krome = %"ISYM, &use_krome);
    /* if the line is suspicious, issue a warning */
    printf("\n use_krome %d \n", use_krome);
    if (ret == 0 && strstr(line, "=") && strstr(line, "DrivenFlow"))
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

  }

  ReadEquilibriumTable("equilibrium_table_60_030-Zsun.h5", MetaData.Time); 


  /* thermodynamic initial values */

  if (MultiSpecies == 0) {
    if (EquationOfState == 0)
      SoundSpeed = sqrt(Gamma * DrivenFlowPressure / DrivenFlowDensity);
    else
      SoundSpeed = IsothermalSoundSpeed;
    
  } else {
    fprintf(stderr,"DrivenFlowInitialize: Multispecies != 0 untested at this point.\n");
    if (EquationOfState == 0)
      SoundSpeed = sqrt(Gamma * DrivenFlowPressure / DrivenFlowDensity);
    else
      SoundSpeed = IsothermalSoundSpeed;
    //return FALSE;
  }

  if (SelfGravity) {
      fprintf(stderr,"DrivenFlowInitialize: SelfGravity untested at this point.\n");
      return FALSE;
  }

  if ((HydroMethod != MHD_RK) && (HydroMethod != HD_RK) && (HydroMethod != MHD_Li)) {
      fprintf(stderr,"DrivenFlowInitialize: Only support for MUSCL framework and MHDCT at this point.\n");
      return FALSE;
  }

  if (MetaData.TopGridRank != 3) {
      fprintf(stderr,"DrivenFlowInitialize: Only 3D tested at this point.\n");
      return FALSE;
  }


  // set proper internal unit for magnetic field
  DrivenFlowMagField /= sqrt(4*pi);
  
  /* compute characteristic velocity from speed of sound and Mach numbers */

  for (int dim = 0; dim < MetaData.TopGridRank; dim++) {
      DrivenFlowVelocity[dim] = DrivenFlowMach[dim] * SoundSpeed;
      DrivenFlowDomainLength[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
      if (debug)
          printf("dim = %"ISYM" vel = %"FSYM" len = %"FSYM"\n",
            dim,DrivenFlowVelocity[dim],DrivenFlowDomainLength[dim]);
  }

  /* Begin grid initialization */

  HierarchyEntry *CurrentGrid; // all level 0 grids on this processor first
  CurrentGrid = &TopGrid;
  while (CurrentGrid != NULL) {
      if (CurrentGrid->GridData->DrivenFlowInitializeGrid(DrivenFlowDensity,
              DrivenFlowPressure,DrivenFlowMagField,SetBaryonFields) == FAIL) {
          fprintf(stderr, "Error in DrivenFlowInitializeGrid.\n");
          return FAIL;
      }
      CurrentGrid = CurrentGrid->NextGridThisLevel;
  }
  
  if (SetBaryonFields) {
  /* create a stochasitc forcing object with the specified parameters */
  Forcing.Init(MetaData.TopGridRank,
           DrivenFlowProfile,
           DrivenFlowAlpha,
           DrivenFlowDomainLength,
           DrivenFlowBandWidth,
           DrivenFlowVelocity,
           DrivenFlowAutoCorrl,
           DrivenFlowWeight,
           DrivenFlowSeed);

  /* set up field names and units */

  int count = 0;
  DataLabel[count++] = DensName;
  
  DataLabel[count++] = Vel1Name;
  DataLabel[count++] = Vel2Name;
  DataLabel[count++] = Vel3Name;


  if(EquationOfState == 0)
    DataLabel[count++] = TEName;
  if (DualEnergyFormalism)
    DataLabel[count++] = GEName;

  if ( UseMHD ) {
    DataLabel[count++] = BxName;
    DataLabel[count++] = ByName;
    DataLabel[count++] = BzName;
  }
  if( HydroMethod == MHD_RK ){
    DataLabel[count++] = PhiName;
  }
 if (MultiSpecies) {
   DataLabel[count++] = ElectronName;
   DataLabel[count++] = HIName;
   DataLabel[count++] = HIIName;
   DataLabel[count++] = HeIName;
   DataLabel[count++] = HeIIName;
   DataLabel[count++] = HeIIIName;
   if (MultiSpecies > 1) {
     DataLabel[count++] = HMName;
     DataLabel[count++] = H2IName;
     DataLabel[count++] = H2IIName;
   }
   if (MultiSpecies > 2) {
     DataLabel[count++] = DIName;
     DataLabel[count++] = DIIName;
     DataLabel[count++] = HDIName;
   }
   DataLabel[count++] = CMName;  
   DataLabel[count++] = OMName; 
   DataLabel[count++] = CIName; 
   DataLabel[count++] = OIName; 
   DataLabel[count++] = OHIName; 
   DataLabel[count++] = COIName; 
   DataLabel[count++] = CHIName; 
   DataLabel[count++] = CH2IName; 
   DataLabel[count++] = C2IName; 
   DataLabel[count++] = HCOIName; 
   DataLabel[count++] = H2OIName; 
   DataLabel[count++] = O2IName; 
   DataLabel[count++] = CO_TOTALIName; 
   DataLabel[count++] = H2O_TOTALIName; 
   DataLabel[count++] = CIIName; 
   DataLabel[count++] = OIIName; 
   DataLabel[count++] = HOCIIName; 
   DataLabel[count++] = HCOIIName; 
   DataLabel[count++] = H3IIName; 
   DataLabel[count++] = CHIIName; 
   DataLabel[count++] = CH2IIName; 
   DataLabel[count++] = COIIName; 
   DataLabel[count++] = CH3IIName; 
   DataLabel[count++] = OHIIName; 
   DataLabel[count++] = H2OIIName; 
   DataLabel[count++] = H3OIIName; 
   DataLabel[count++] = O2IIName; 
 }
  MHDCTSetupFieldLabels();

    DataLabel[count++] = StochAccel1Name;
    DataLabel[count++] = StochAccel2Name;
    DataLabel[count++] = StochAccel3Name;

  /* Write parameters to parameter output file */

  if (debug) printf("DrivenFlowInitialize: writing parameters to output file.\n");

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "Profile    = %"ISYM"\n\n", DrivenFlowProfile);

    fprintf(Outfptr, "Alpha1     = %"ISYM"\n",   DrivenFlowAlpha[0]);
    fprintf(Outfptr, "Alpha2     = %"ISYM"\n",   DrivenFlowAlpha[1]);
    fprintf(Outfptr, "Alpha3     = %"ISYM"\n\n", DrivenFlowAlpha[2]);

    fprintf(Outfptr, "BandWidth1 = %"FSYM"\n",   DrivenFlowBandWidth[0]);
    fprintf(Outfptr, "BandWidth2 = %"FSYM"\n",   DrivenFlowBandWidth[1]);
    fprintf(Outfptr, "BandWidth3 = %"FSYM"\n\n", DrivenFlowBandWidth[2]);

    fprintf(Outfptr, "AutoCorrl1 = %"FSYM"\n",   DrivenFlowAutoCorrl[0]);
    fprintf(Outfptr, "AutoCorrl2 = %"FSYM"\n",   DrivenFlowAutoCorrl[1]);
    fprintf(Outfptr, "AutoCorrl3 = %"FSYM"\n\n", DrivenFlowAutoCorrl[2]);

    fprintf(Outfptr, "SolnWeight = %"FSYM"\n\n", DrivenFlowWeight);

    fprintf(Outfptr, "Density    = %"FSYM"\n",   DrivenFlowDensity);
    fprintf(Outfptr, "Pressure   = %"FSYM"\n\n", DrivenFlowPressure);

    fprintf(Outfptr, "Velocity1  = %"FSYM"\n",   DrivenFlowVelocity[0]);
    fprintf(Outfptr, "Velocity2  = %"FSYM"\n",   DrivenFlowVelocity[1]);
    fprintf(Outfptr, "Velocity3  = %"FSYM"\n\n", DrivenFlowVelocity[2]);
  }
  }
  return SUCCESS;
}
