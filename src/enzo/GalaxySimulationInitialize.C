/***********************************************************************
/
/  INITIALIZE A GALAXY SIMULATION
/
/  written by: Greg Bryan
/  date:       May, 1998
/  modified1:  Elizabeth Tasker, March 2004
/
/  PURPOSE:
/
/    Set up a number of spherical objects
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

// This routine intializes a new simulation based on the parameter file.
//

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

#include <string.h>
#include <stdio.h>
#include <math.h>

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
#define VCIRC_TABLE_LENGTH 10000
void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);

int GetUnits(float *DensityUnits, float *LengthUnits,
       float *TemperatureUnits, float *TimeUnits,
       float *VelocityUnits, double *MassUnits, FLOAT Time);
//int Enzo_Version;
//function prototypes that are for only the isogal build 
void MHDCTSetupFieldLabels();
float GetMagneticUnits(float DensityUnits, float LengthUnits, float TimeUnits);		
int ReadEquilibriumTable(char * name, FLOAT Time);

int nlines(const char* fname);
int InitializeParticles(grid *thisgrid_orig, HierarchyEntry &TopGrid, TopGridData &MetaData, FLOAT * Center);
int ReadParticlesFromFile(PINT *Number, int *Type, FLOAT *Position[],
			    float *Velocity[], float* Mass, const char* fname,
			    Eint32 particle_type, int &c, FLOAT dx, FLOAT * Center);
float InterpolateVcircTable(FLOAT radius, FLOAT * VCircRadius, float * VCircVelocity);
void ReadInVcircData(FLOAT * VCircRadius, float * VCircVelocity);

int GalaxySimulationInitialize(FILE *fptr, FILE *Outfptr, 
			  HierarchyEntry &TopGrid, TopGridData &MetaData, ExternalBoundary &Exterior, int SetBaryons)
{
  char *DensName    = "Density";
  char *TEName      = "TotalEnergy";
  char *GEName      = "GasEnergy";
  char *Vel1Name    = "x-velocity";
  char *Vel2Name    = "y-velocity";
  char *Vel3Name    = "z-velocity";
  char *CRName      = "CREnergyDensity";
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
  char *BxName      = "Bx";
  char *ByName      = "By";
  char *BzName      = "Bz";
  char *PhiName     = "Phi";
  char *PotName     = "PotentialField";
  char *Acceleration0Name = "Acceleration_x";
  char *Acceleration1Name = "Acceleration_y";
  char *Acceleration2Name = "Acceleration_z";
  /* declarations */

  char  line[MAX_LINE_LENGTH];
  int   dim, ret, level, disk, i;
  static int GalaxySimulationDebugHold = FALSE;
  /* make sure it is 3D */
  
  if (MetaData.TopGridRank != 3) {
    ENZO_VFAIL("Cannot do GalaxySimulation in %"ISYM" dimension(s)\n", MetaData.TopGridRank)
  }

  /* set default parameters */

  FLOAT GalaxySimulationGasMass,//this was float, I changed to FLOAT to make versions
    GalaxySimulationGalaxyMass,//since i dont think this should break anything
    GalaxySimulationCR,
    GalaxySimulationDiskTemperature,
    GalaxySimulationAngularMomentum[MAX_DIMENSION],
    GalaxySimulationUniformVelocity[MAX_DIMENSION],
    GalaxySimulationUniformDensity,
    GalaxySimulationUniformCR,
    GalaxySimulationUniformEnergy;

  FLOAT GalaxySimulationDiskRadius,
    GalaxySimulationDiskPosition[MAX_DIMENSION],
    GalaxySimulationDiskScaleHeightz,
    GalaxySimulationDiskScaleHeightR,
    GalaxySimulationTruncationRadius;


  FLOAT GalaxySimulationInitialTemperature,
    GalaxySimulationDarkMatterConcentrationParameter,
    GalaxySimulationInflowTime,
    GalaxySimulationInflowDensity;
	
  int GalaxySimulationGasHalo;
  FLOAT GalaxySimulationGasHaloScaleRadius,
	GalaxySimulationGasHaloDensity;

  int GalaxySimulationIterateRebuildHierarchy = TRUE; // IF you have a solution senstive AMR strategy, this should be true.
  int GalaxySimulationStaticHierarchyAfterInit = TRUE; // IF you have a solution senstive AMR strategy, this should be true.

  int   GalaxySimulationRefineAtStart,
    GalaxySimulationUseMetallicityField;
		
  float ZeroBField[3] = {0.0, 0.0, 0.0};
//definitions of variables used in the isogal build
	FLOAT GalaxySimulationEquilibrateChem, 
	GalaxySimulationDiskDensityCap,	     
	GalaxySimulationGasHaloDensity2,
	GalaxySimulationGasHaloTemperature,
	GalaxySimulationGasHaloAlpha,
	GalaxySimulationGasHaloZeta,
	GalaxySimulationGasHaloZeta2,
	GalaxySimulationGasHaloCoreEntropy,
	GalaxySimulationGasHaloRatio,
	GalaxySimulationGasHaloMetallicity,
	GalaxySimulationDiskMetallicityEnhancementFactor,
	GalaxySimulationGasHaloRotationScaleVelocity,
	GalaxySimulationGasHaloRotationScaleRadius,
	GalaxySimulationGasHaloRotationIndex;
	

  	char GalaxySimulationEquilibriumFile[MAX_LINE_LENGTH] = "equilibrium_table_50.h5"; 	       
  	int GalaxySimulationGasHaloRotation;		
  	FLOAT GalaxySimulationInitialBfield[3] = {0.0, 0.0, 0.0};
  	int GalaxySimulationInitialBfieldTopology= 0; //Uniform cartesiani
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];

  /* Default Values */
  //shared defaults
  GalaxySimulationRefineAtStart      = TRUE;
  GalaxySimulationUseMetallicityField  = FALSE;
  GalaxySimulationInitialTemperature = 1000.0;
  GalaxySimulationDiskRadius         = 0.2;      // CODE UNITS
  GalaxySimulationDiskTemperature    = 1.e4;     // [K]
  GalaxySimulationDiskScaleHeightz   = 325e-6;   // Mpc
  GalaxySimulationDiskScaleHeightR   = 3500e-6;  // Mpc
  GalaxySimulationTruncationRadius   = .026; // [ Mpc ]
  GalaxySimulationDarkMatterConcentrationParameter = 10;//messing with this to try to avoid negative or zero energy error 
  GalaxySimulationGasMass            = 4.0e10;
  GalaxySimulationGalaxyMass         = 1.0e12;
  GalaxySimulationDiskTemperature    = 1000.0;
  GalaxySimulationGasHalo            = 0; // uniform halo w/ densicm and UniformTemperature
  GalaxySimulationGasHaloScaleRadius = .001; // Mpc
  GalaxySimulationGasHaloDensity     = 1.8e-27; // cgs
  GalaxySimulationInflowTime         = -1;
  GalaxySimulationInflowDensity      = 0;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    GalaxySimulationDiskPosition[dim] = 0.5*(DomainLeftEdge[dim] +
					     DomainRightEdge[dim]);
    GalaxySimulationAngularMomentum[dim] = 0.0;
    GalaxySimulationUniformVelocity[dim] = 0.0;
  }
  GalaxySimulationUniformDensity = 1.0E-28;
  GalaxySimulationUniformEnergy = 1.0;
  GalaxySimulationCR = .01;
  GalaxySimulationUniformCR = .01;
  //default values of quantities only defined in isogal build
GalaxySimulationDiskDensityCap = 0.0; 
GalaxySimulationEquilibrateChem = 0; //bool
GalaxySimulationGasHaloDensity2    = 0.0; // cgs
GalaxySimulationGasHaloTemperature = 1.0e+6;  // Kelvin
GalaxySimulationGasHaloAlpha       = 2.0/3.0;  // unitless
GalaxySimulationGasHaloZeta        = 0;
GalaxySimulationGasHaloZeta2       = 0;
GalaxySimulationGasHaloCoreEntropy = 5.0;  // keV cm^2
GalaxySimulationGasHaloRatio       = 10; // ratio of cooling time to freefall time
GalaxySimulationGasHaloMetallicity = 0.1; // Zsun
GalaxySimulationGasHaloRotation    = 0; // off
GalaxySimulationGasHaloRotationScaleVelocity = 180.0; // km/s
GalaxySimulationGasHaloRotationScaleRadius   = 10.0; // kpc
GalaxySimulationGasHaloRotationIndex         = 0.0; // unitless
GalaxySimulationDiskMetallicityEnhancementFactor = 3.0; // w.r.t to halo metallicity
  
  /* read input from file */
//setup needed for chemistry in the isogal build
char *dummy = new char[MAX_LINE_LENGTH];
dummy[0] = 0;

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
    
    ret = 0;
    //read parameter file for common variables between builds
    ret += sscanf(line, "GalaxySimulationEnzoVersion		       = %"ISYM, &Enzo_Version); //added int value to specify which version of enzo to use 
    ret += sscanf(line, "GalaxySimulationRefineAtStart = %"ISYM,
		  &GalaxySimulationRefineAtStart);
    ret += sscanf(line, "GalaxySimulationIterateRebuildHierarchy = %"ISYM,
		  &GalaxySimulationIterateRebuildHierarchy);
    ret += sscanf(line, "GalaxySimulationStaticHierarchyAfterInit = %"ISYM,
		  &GalaxySimulationStaticHierarchyAfterInit);

    ret += sscanf(line, "GalaxySimulationUseMetallicityField = %"ISYM,
		  &GalaxySimulationUseMetallicityField);
    ret += sscanf(line, "GalaxySimulationInitialTemperature = %"FSYM,
		  &GalaxySimulationInitialTemperature);
    ret += sscanf(line, "GalaxySimulationUniformDensity = %"FSYM,
		  &GalaxySimulationUniformDensity);
    ret += sscanf(line, "GalaxySimulationUniformVelocity = %"FSYM" %"FSYM" %"FSYM,
                  &GalaxySimulationUniformVelocity[0], &GalaxySimulationUniformVelocity[1],
                  &GalaxySimulationUniformVelocity[2]);
    ret += sscanf(line, "GalaxySimulationDiskRadius = %"PSYM,
		  &GalaxySimulationDiskRadius);
    ret += sscanf(line, "GalaxySimulationGalaxyMass = %"FSYM,
		  &GalaxySimulationGalaxyMass);
    ret += sscanf(line, "GalaxySimulationGasMass = %"FSYM,
		  &GalaxySimulationGasMass);
    ret += sscanf(line, "GalaxySimulationCR = %"FSYM,
		  &GalaxySimulationCR);
    ret += sscanf(line, "GalaxySimulationUniformCR = %"FSYM,
		  &GalaxySimulationUniformCR);
    ret += sscanf(line, "GalaxySimulationDiskPosition = %"PSYM" %"PSYM" %"PSYM, 
		  &GalaxySimulationDiskPosition[0],
		  &GalaxySimulationDiskPosition[1],
		  &GalaxySimulationDiskPosition[2]);
    ret += sscanf(line, "GalaxySimulationDiskScaleHeightz = %"PSYM,
		  &GalaxySimulationDiskScaleHeightz);
    ret += sscanf(line, "GalaxySimulationDiskScaleHeightR = %"PSYM,
		  &GalaxySimulationDiskScaleHeightR);
    ret += sscanf(line, "GalaxySimulationTruncationRadius = %"PSYM,
		  &GalaxySimulationTruncationRadius);
    ret += sscanf(line, "GalaxySimulationDarkMatterConcentrationParameter = %"FSYM,
		  &GalaxySimulationDarkMatterConcentrationParameter);
    ret += sscanf(line, "GalaxySimulationDiskTemperature = %"FSYM,
		  &GalaxySimulationDiskTemperature);
    ret += sscanf(line, "GalaxySimulationGasHalo = %"ISYM,
		  &GalaxySimulationGasHalo);
    ret += sscanf(line, "GalaxySimulationGasHaloScaleRadius = %"FSYM,
		  &GalaxySimulationGasHaloScaleRadius);
    ret += sscanf(line, "GalaxySimulationGasHaloDensity = %"FSYM,
		  &GalaxySimulationGasHaloDensity);
    ret += sscanf(line, "GalaxySimulationInflowTime = %"FSYM,
		  &GalaxySimulationInflowTime);
    ret += sscanf(line, "GalaxySimulationInflowDensity = %"FSYM,
		  &GalaxySimulationInflowDensity);
    ret += sscanf(line, "GalaxySimulationAngularMomentum = %"FSYM" %"FSYM" %"FSYM,
		  &GalaxySimulationAngularMomentum[0],
		  &GalaxySimulationAngularMomentum[1],
		  &GalaxySimulationAngularMomentum[2]);
	    ret += sscanf(line, "GalaxySimulationDiskDensityCap = %"FSYM,
			  &GalaxySimulationDiskDensityCap);    
	    ret += sscanf(line, "GalaxySimulationEquilibrateChem = %"FSYM,
			  &GalaxySimulationEquilibrateChem);
	    if (sscanf(line, "GalaxySimulationEquilibriumFile = %s", dummy) == 1) {
	      strcpy(GalaxySimulationEquilibriumFile, dummy);
	      ret++;
    }
	    ret += sscanf(line, "GalaxySimulationGasHaloDensity2 = %"FSYM,
			  &GalaxySimulationGasHaloDensity2);
	    ret += sscanf(line, "GalaxySimulationGasHaloTemperature = %"FSYM,
			  &GalaxySimulationGasHaloTemperature);
	    ret += sscanf(line, "GalaxySimulationGasHaloAlpha = %"FSYM,
			  &GalaxySimulationGasHaloAlpha);
	    ret += sscanf(line, "GalaxySimulationGasHaloZeta = %"FSYM,
			  &GalaxySimulationGasHaloZeta);
	    ret += sscanf(line, "GalaxySimulationGasHaloZeta2 = %"FSYM,
			  &GalaxySimulationGasHaloZeta2);
	    ret += sscanf(line, "GalaxySimulationGasHaloCoreEntropy = %"FSYM,
			  &GalaxySimulationGasHaloCoreEntropy);
	    ret += sscanf(line, "GalaxySimulationGasHaloRatio = %"FSYM,
			  &GalaxySimulationGasHaloRatio);
	    ret += sscanf(line, "GalaxySimulationGasHaloMetallicity = %"FSYM,
			  &GalaxySimulationGasHaloMetallicity);
	    ret += sscanf(line, "GalaxySimulationGasHaloRotation = %"ISYM,
			  &GalaxySimulationGasHaloRotation);
	    ret += sscanf(line, "GalaxySimulationGasHaloRotationScaleVelocity = %"FSYM,
			  &GalaxySimulationGasHaloRotationScaleVelocity);
	    ret += sscanf(line, "GalaxySimulationGasHaloRotationScaleRadius = %"FSYM,
			  &GalaxySimulationGasHaloRotationScaleRadius);
	    ret += sscanf(line, "GalaxySimulationGasHaloRotationIndex = %"FSYM,
			  &GalaxySimulationGasHaloRotationIndex);
	    ret += sscanf(line, "GalaxySimulationDiskMetallicityEnhancementFactor = %"FSYM,
			  &GalaxySimulationDiskMetallicityEnhancementFactor);
	    ret += sscanf(line, "GalaxySimulationInitialBfield = %"FSYM" %"FSYM" %"FSYM, 
	      &GalaxySimulationInitialBfield[0],
	      &GalaxySimulationInitialBfield[1],
	      &GalaxySimulationInitialBfield[2]);
	    ret += sscanf(line, "GalaxySimulationInitialBfieldTopology = %"ISYM,
	      &GalaxySimulationInitialBfieldTopology);
	    ret += sscanf(line, "DiskGravityDarkMatterMass = %"FSYM, &DiskGravityDarkMatterMass);
	    ret += sscanf(line, "DiskGravityDarkMatterConcentration = %"FSYM, &DiskGravityDarkMatterConcentration); 
	    ret += sscanf(line, "HydroMethod = %"ISYM, &HydroMethod);
	    ret += sscanf(line, "RiemannSolver = %"ISYM, &RiemannSolver);
	    ret += sscanf(line, "ReconstructionMethod = %"ISYM, &ReconstructionMethod);
        ret += sscanf(line, "GalaxySimulationDebugHold = %"ISYM, &GalaxySimulationDebugHold);
    
    /* if the line is suspicious, issue a warning */
    if (ret == 0 && strstr(line, "=") && strstr(line, "GalaxySimulation") 
	&& line[0] != '#' && !strstr(line,"RPSWind") && !strstr(line,"PreWind"))
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

  } // end input from parameter file
  FLOAT VCircRadius[VCIRC_TABLE_LENGTH];
  float VCircVelocity[VCIRC_TABLE_LENGTH];
  ReadInVcircData(VCircRadius, VCircVelocity);
  
  /* fix wind values wrt units */
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits;
  double MassUnits;
  if (GetUnits(&DensityUnits, &LengthUnits,&TemperatureUnits, &TimeUnits,
               &VelocityUnits, &MassUnits, MetaData.Time) == FAIL){
    return FAIL;
  }
  if(Enzo_Version==2){
	  // If using DiskGravity, make two GalaxySimulation parameters consistent
	  if (DiskGravity > 0) {
	    GalaxySimulationGalaxyMass = DiskGravityDarkMatterMass;
	    GalaxySimulationDarkMatterConcentrationParameter = DiskGravityDarkMatterConcentration;//this line has a problem with the RHS being zero and making things blow up
	  }
	  if (GalaxySimulationEquilibrateChem){
	    ReadEquilibriumTable(GalaxySimulationEquilibriumFile, MetaData.Time);
      }

	  delete [] dummy;
	  if( UseMHD ){
	      float MagneticUnits = GetMagneticUnits(DensityUnits, LengthUnits, TimeUnits);
	      for( dim=0; dim<3; dim++ ){
		  GalaxySimulationInitialBfield[dim] /=MagneticUnits;
	      }
	  }

		
  }
  //convert things to code units
  GalaxySimulationRPSWindDensity = GalaxySimulationRPSWindDensity/DensityUnits;
  GalaxySimulationRPSWindPressure = GalaxySimulationRPSWindPressure/DensityUnits/LengthUnits/LengthUnits*TimeUnits*TimeUnits;
  GalaxySimulationRPSWindVelocity[0] = GalaxySimulationRPSWindVelocity[0]/LengthUnits*TimeUnits;
  GalaxySimulationRPSWindVelocity[1] = GalaxySimulationRPSWindVelocity[1]/LengthUnits*TimeUnits;
  GalaxySimulationRPSWindVelocity[2] = GalaxySimulationRPSWindVelocity[2]/LengthUnits*TimeUnits;
  GalaxySimulationRPSWindShockSpeed = GalaxySimulationRPSWindShockSpeed/LengthUnits*TimeUnits;
  GalaxySimulationRPSWindDelay = GalaxySimulationRPSWindDelay/TimeUnits;

  /* Align gaseous and stellar disks */
  if( DiskGravity > 0 ){
    for( i = 0 ; i < MAX_DIMENSION ; i++ )
      DiskGravityAngularMomentum[i] = GalaxySimulationAngularMomentum[i];
  } // end DiskGravity if

  /* set up grid */

  HierarchyEntry *CurrentGrid; // all level 0 grids on this processor first
  CurrentGrid = &TopGrid;
  while (CurrentGrid != NULL) {
	  CurrentGrid->GridData->GalaxySimulationInitializeGrid(GalaxySimulationDiskRadius,
					GalaxySimulationGalaxyMass, 
					GalaxySimulationGasMass,
					GalaxySimulationDiskPosition, 
					GalaxySimulationDiskScaleHeightz,
					GalaxySimulationDiskScaleHeightR,
					GalaxySimulationTruncationRadius, 
					GalaxySimulationDiskDensityCap,
					GalaxySimulationDarkMatterConcentrationParameter,
					GalaxySimulationDiskTemperature, 
					GalaxySimulationInitialTemperature,
					GalaxySimulationUniformDensity,
					GalaxySimulationEquilibrateChem,
					GalaxySimulationGasHalo,
					GalaxySimulationGasHaloScaleRadius,
					GalaxySimulationGasHaloDensity,
					GalaxySimulationGasHaloDensity2,
					GalaxySimulationGasHaloTemperature,
					GalaxySimulationGasHaloAlpha,
					GalaxySimulationGasHaloZeta,
					GalaxySimulationGasHaloZeta2,
					GalaxySimulationGasHaloCoreEntropy,
					GalaxySimulationGasHaloRatio,
					GalaxySimulationGasHaloMetallicity,
					GalaxySimulationGasHaloRotation,
					GalaxySimulationGasHaloRotationScaleVelocity,
					GalaxySimulationGasHaloRotationScaleRadius,
					GalaxySimulationGasHaloRotationIndex,
					GalaxySimulationDiskMetallicityEnhancementFactor,
					GalaxySimulationAngularMomentum,
					GalaxySimulationUniformVelocity,
					GalaxySimulationUseMetallicityField,
					GalaxySimulationInflowTime,
					GalaxySimulationInflowDensity,0,
					GalaxySimulationInitialBfield,
					GalaxySimulationInitialBfieldTopology,
	VCircRadius, VCircVelocity,
					GalaxySimulationCR,
          SetBaryons
							       );
    CurrentGrid = CurrentGrid->NextGridThisLevel;
  }

  if ( DiskGravity + PointSourceGravity == FALSE ){
      InitializeParticles(TopGrid.GridData, TopGrid, MetaData, GalaxySimulationDiskPosition);
  }
  
  /* Convert minimum initial overdensity for refinement to mass
     (unless MinimumMass itself was actually set). */
switch(Enzo_Version){
	case 1: //stock density minimum refinement mass routine 
	  if (MinimumMassForRefinement[0] == FLOAT_UNDEFINED) {
	    MinimumMassForRefinement[0] = MinimumOverDensityForRefinement[0];
	    for (int dim = 0; dim < MetaData.TopGridRank; dim++)
	      MinimumMassForRefinement[0] *=(DomainRightEdge[dim]-DomainLeftEdge[dim])/
		float(MetaData.TopGridDims[dim]);
	  }
	  break; 
	case 2: //isogal version
	  for (i = 0; i < MAX_FLAGGING_METHODS; i++)
	    if (MinimumMassForRefinement[i] == FLOAT_UNDEFINED) {
	      MinimumMassForRefinement[i] = MinimumOverDensityForRefinement[i];
	      for (dim = 0; dim < MetaData.TopGridRank; dim++)
		MinimumMassForRefinement[i] *=
		(DomainRightEdge[dim]-DomainLeftEdge[dim])/
		float(MetaData.TopGridDims[dim]);
	    }
	  break;
}
  /* If requested, refine the grid to the desired level. */
if(SetBaryons){
  if (GalaxySimulationRefineAtStart) {

    /* Declare, initialize and fill out the LevelArray. */

    LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
      LevelArray[level] = NULL;
    AddLevel(LevelArray, &TopGrid, 0);

    /* Add levels to the maximum depth or until no new levels are created,
       and re-initialize the level after it is created. */

    for (level = 0; level < MaximumRefinementLevel; level++) {
        //For refinement strategies that depend on the initial conditions on each level being resolved,
        //you'll need to iterate. For pre determined AMR structures, you do not.
        if ( GalaxySimulationIterateRebuildHierarchy || level == 0 )
            if (RebuildHierarchy(&MetaData, LevelArray, level) == FAIL) {
                fprintf(stderr, "Error in RebuildHierarchy.\n");
                return FAIL;
            }

        if (LevelArray[level+1] == NULL)
            break;
        if ( debug ) fprintf(stderr,"Build level %d\n",level+1);
        LevelHierarchyEntry *Temp = LevelArray[level+1];
        while (Temp != NULL) {
            if (Temp->GridData->GalaxySimulationInitializeGrid(GalaxySimulationDiskRadius,
                        GalaxySimulationGalaxyMass, 
                        GalaxySimulationGasMass,
                        GalaxySimulationDiskPosition, 
                        GalaxySimulationDiskScaleHeightz,
                        GalaxySimulationDiskScaleHeightR,
                        GalaxySimulationTruncationRadius, 
                        GalaxySimulationDiskDensityCap,
                        GalaxySimulationDarkMatterConcentrationParameter,
                        GalaxySimulationDiskTemperature, 
                        GalaxySimulationInitialTemperature,
                        GalaxySimulationUniformDensity,
                        GalaxySimulationEquilibrateChem,
                        GalaxySimulationGasHalo,
                        GalaxySimulationGasHaloScaleRadius,
                        GalaxySimulationGasHaloDensity,
                        GalaxySimulationGasHaloDensity2,
                        GalaxySimulationGasHaloTemperature,
                        GalaxySimulationGasHaloAlpha,
                        GalaxySimulationGasHaloZeta,
                        GalaxySimulationGasHaloZeta2,
                        GalaxySimulationGasHaloCoreEntropy,
                        GalaxySimulationGasHaloRatio,
                        GalaxySimulationGasHaloMetallicity,
                        GalaxySimulationGasHaloRotation,
                        GalaxySimulationGasHaloRotationScaleVelocity,
                        GalaxySimulationGasHaloRotationScaleRadius,
                        GalaxySimulationGasHaloRotationIndex,
                        GalaxySimulationDiskMetallicityEnhancementFactor,
                        GalaxySimulationAngularMomentum,
                        GalaxySimulationUniformVelocity,
                        GalaxySimulationUseMetallicityField,
                        GalaxySimulationInflowTime,
                        GalaxySimulationInflowDensity,level,
                        GalaxySimulationInitialBfield,
                        GalaxySimulationInitialBfieldTopology,
			VCircRadius, VCircVelocity,
                        GalaxySimulationCR, 
                        SetBaryons
                            )
                            == FAIL) {
                                ENZO_FAIL("Error in GalaxySimulationInitialize[Sub]Grid.");
                            }// end subgrid if
            else {
                Temp->GridData->_GalaxySimulationInitialization = 1;
            }
            Temp = Temp->NextGridThisLevel;
        } // end: loop over levels
    }//end of switch statement
    /* Loop back from the bottom, restoring the consistency among levels. */

    MHD_ProjectE=FALSE;
    MHD_ProjectB=TRUE;
    for (level = MaximumRefinementLevel; level > 0; level--) {
      LevelHierarchyEntry *Temp = LevelArray[level];
      while (Temp != NULL) {
	if (Temp->GridData->ProjectSolutionToParentGrid(
				   *LevelArray[level-1]->GridData) == FAIL) {
	  fprintf(stderr, "Error in grid->ProjectSolutionToParentGrid.\n");
	  return FAIL;
	}
	Temp = Temp->NextGridThisLevel;
      }
    }
    MHD_ProjectE=TRUE;
    MHD_ProjectB=FALSE;

    if ( GalaxySimulationStaticHierarchyAfterInit && 0)
        MetaData.StaticHierarchy=TRUE;//this breaks things for now

  } // end: if (GalaxySimulationRefineAtStart)
} 
  /* If Galaxy is Subject to ICM Wind, Initialize the exterior */

  if ( GalaxySimulationRPSWind > 0 ) {
    Exterior.Prepare(TopGrid.GridData);
	
    const int MAX_BNDRY_VARS = 6;
    float InflowValue[MAX_BNDRY_VARS], Dummy[MAX_BNDRY_VARS];
    InflowValue[0] = GalaxySimulationRPSWindDensity;
    InflowValue[1] = GalaxySimulationRPSWindPressure/(Gamma-1.0)/GalaxySimulationRPSWindDensity;
    if (HydroMethod != 2) {
      InflowValue[1] = InflowValue[1] + 0.5*(   pow(GalaxySimulationRPSWindVelocity[0],2)
	                                            + pow(GalaxySimulationRPSWindVelocity[1],2)
	                                            + pow(GalaxySimulationRPSWindVelocity[2],2));
    }
    InflowValue[2] = GalaxySimulationRPSWindVelocity[0];
    InflowValue[3] = GalaxySimulationRPSWindVelocity[1];
    InflowValue[4] = GalaxySimulationRPSWindVelocity[2];
    if (GalaxySimulationUseMetallicityField)
      InflowValue[5] = 1.0e-10;
  
    if (Exterior.InitializeExternalBoundaryFace(0, inflow, outflow, InflowValue,
						Dummy) == FAIL) {
      fprintf(stderr, "Error in InitializeExternalBoundaryFace.\n");
      return FAIL;
    }
    if (MetaData.TopGridRank > 1)
      Exterior.InitializeExternalBoundaryFace(1, inflow, outflow,
					      InflowValue, Dummy);
    if (MetaData.TopGridRank > 2)
      Exterior.InitializeExternalBoundaryFace(2, inflow, outflow,
					      InflowValue, Dummy);
	
    /* Set Global Variables for RPS Wind (see ExternalBoundary_SetGalaxySimulationBoundary.C)*/

    GalaxySimulationRPSWindDelay += TopGrid.GridData->ReturnTime();
    GalaxySimulationRPSWindTotalEnergy = InflowValue[1]; 	
    GalaxySimulationPreWindDensity     = GalaxySimulationUniformDensity/DensityUnits;
    GalaxySimulationPreWindTotalEnergy = GalaxySimulationInitialTemperature/TemperatureUnits/((Gamma-1.0)*0.6); 
    GalaxySimulationPreWindVelocity[0] = 0.0;
    GalaxySimulationPreWindVelocity[1] = 0.0;
    GalaxySimulationPreWindVelocity[2] = 0.0;
  }

  // If we used the Equilibrium Table, delete it
  if (GalaxySimulationEquilibrateChem && Enzo_Version == 2){
    if (MultiSpecies) {
      delete [] EquilibriumTable.HI;
      delete [] EquilibriumTable.HII;
      delete [] EquilibriumTable.HeI;
      delete [] EquilibriumTable.HeII;
      delete [] EquilibriumTable.HeIII;
      delete [] EquilibriumTable.de;
      if (MultiSpecies > 1) {
        delete [] EquilibriumTable.HM;
        delete [] EquilibriumTable.H2I;
        delete [] EquilibriumTable.H2II;
      }
      if (MultiSpecies > 2) {
        delete [] EquilibriumTable.DI;
        delete [] EquilibriumTable.DII;
        delete [] EquilibriumTable.HDI;
      }
    }
  }
if(SetBaryons){
 /* set up field names and units */
 int count = 0;
 DataLabel[count++] = DensName;
 DataLabel[count++] = TEName;
 if (DualEnergyFormalism)
     DataLabel[count++] = GEName;
 if ( WritePotential )
     DataLabel[count++] = PotName;
 if(WriteAcceleration){
     DataLabel[count++] = Acceleration0Name;
     DataLabel[count++] = Acceleration1Name;
     DataLabel[count++] = Acceleration2Name;
 }
 DataLabel[count++] = Vel1Name;
 if(MetaData.TopGridRank > 1)
     DataLabel[count++] = Vel2Name;
 if(MetaData.TopGridRank > 2)
   DataLabel[count++] = Vel3Name;
  if( UseMHD){
      DataLabel[count++] = BxName;
      DataLabel[count++] = ByName;
      DataLabel[count++] = BzName;
  }
  if (HydroMethod == MHD_RK){
      DataLabel[count++] = PhiName;
  }
 if(CRModel)
   DataLabel[count++] = CRName;
 //fields used in isogal build  was an if here with Enzo_Version 2 removed bc stock toggle fields were not being set right
if(Enzo_Version == 2){ 
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
 }
}
 if(GalaxySimulationUseMetallicityField)
   DataLabel[count++] = MetalName;
 if (StarMakerTypeIaSNe)
   DataLabel[count++] = MetalIaName;
for(int l = 0; l < count; l++)
 for (i = 0; i < count; i++)
   DataUnits[i] = NULL;
 MHDCTSetupFieldLabels();
 /* Write parameters to parameter output file */

 if (MyProcessorNumber == ROOT_PROCESSOR) {

   fprintf(Outfptr, "GalaxySimulationRefineAtStart      = %"ISYM"\n",
	   GalaxySimulationRefineAtStart);
   fprintf(Outfptr, "GalaxySimulationUseMetallicityField          = %"ISYM"\n",
	   GalaxySimulationUseMetallicityField);
   fprintf(Outfptr, "GalaxySimulationInitialTemperature = %"GOUTSYM"\n",
	   GalaxySimulationInitialTemperature);
   fprintf(Outfptr, "GalaxySimulationUniformDensity = %"GOUTSYM"\n",
     GalaxySimulationUniformDensity);
   fprintf(Outfptr, "GalaxySimulationUniformVelocity    = %"GOUTSYM" %"GOUTSYM" %"GOUTSYM"\n",
	   GalaxySimulationUniformVelocity[0], GalaxySimulationUniformVelocity[1],
	   GalaxySimulationUniformVelocity[2]);
   fprintf(Outfptr, "GalaxySimulationDiskRadius = %"GOUTSYM"\n",
	   GalaxySimulationDiskRadius);
   fprintf(Outfptr, "GalaxySimulationGalaxyMass = %"GOUTSYM"\n",
	   GalaxySimulationGalaxyMass);
   fprintf(Outfptr, "GalaxySimulationGasMass = %"GOUTSYM"\n",
	   GalaxySimulationGasMass);
   fprintf(Outfptr, "GalaxySimulationUniformCR = %"GOUTSYM"\n",
     GalaxySimulationUniformCR);
   fprintf(Outfptr, "GalaxySimulationCR = %"GOUTSYM"\n",
     GalaxySimulationCR);
   fprintf(Outfptr, "GalaxySimulationDiskScaleHeightz = %"GOUTSYM"\n",
	   GalaxySimulationDiskScaleHeightz);
   fprintf(Outfptr, "GalaxySimulationDiskScaleHeightR = %"GOUTSYM"\n",
	   GalaxySimulationDiskScaleHeightR);
   fprintf(Outfptr, "GalaxySimulationTruncationRadius = %"GOUTSYM"\n",
     GalaxySimulationTruncationRadius);
   fprintf(Outfptr, "GalaxySimulationDarkMatterConcentrationParameter = %"GOUTSYM"\n",
	   GalaxySimulationDarkMatterConcentrationParameter);
   fprintf(Outfptr, "GalaxySimulationDiskTemperature = %"GOUTSYM"\n",
	   GalaxySimulationDiskTemperature);
   fprintf(Outfptr, "GalaxySimulationGasHalo = %"ISYM"\n",
     GalaxySimulationGasHalo);
   fprintf(Outfptr, "GalaxySimulationGasHaloScaleRadius = %"GOUTSYM"\n",
     GalaxySimulationGasHaloScaleRadius);
   fprintf(Outfptr, "GalaxySimulationGasHaloDensity = %"GOUTSYM"\n",
     GalaxySimulationGasHaloDensity);
   fprintf(Outfptr, "GalaxySimulationInflowTime = %"GOUTSYM"\n",
	   GalaxySimulationInflowTime);
   fprintf(Outfptr, "GalaxySimulationInflowDensity = %"GOUTSYM"\n",
	   GalaxySimulationInflowDensity);
   fprintf(Outfptr, "GalaxySimulationDiskPosition = ");
   WriteListOfFloats(Outfptr, MetaData.TopGridRank, GalaxySimulationDiskPosition);
   fprintf(Outfptr, "GalaxySimulationAngularMomentum = ");
   WriteListOfFloats(Outfptr, MetaData.TopGridRank, GalaxySimulationAngularMomentum);
 }
}

  printf("IterateRebuild %d \n", GalaxySimulationIterateRebuildHierarchy);
#ifdef USE_MPI

 // BWO: this forces the synchronization of the various point source gravity
 // parameters between processors.  If this is not done, things go to pieces!

 MPI_Barrier(MPI_COMM_WORLD);
 MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
 MPI_Bcast(&PointSourceGravityConstant,1,DataType,ROOT_PROCESSOR, MPI_COMM_WORLD);
 MPI_Bcast(&PointSourceGravityCoreRadius,1,DataType,ROOT_PROCESSOR, MPI_COMM_WORLD);
 int dbfile = 0; 
while(GalaxySimulationDebugHold && dbfile == 0){
	FILE *file = fopen("GO", "r"); 
	if(file != NULL){dbfile = 1;}
}
#endif
 return SUCCESS;

}
int InitializeParticles(grid *thisgrid, HierarchyEntry &TopGrid, TopGridData &MetaData, FLOAT * Center){

    int GridRank, dim;
    int *Dims = new int[3]; //hard coded jerk.
    FLOAT Left[MAX_DIMENSION];
    FLOAT Right[MAX_DIMENSION];
    FLOAT CellWidth[MAX_DIMENSION];
    thisgrid->ReturnGridInfo(&GridRank, Dims, Left, Right);
    for ( dim=0;dim<GridRank;dim++){
        CellWidth[dim] = (Right[dim]-Left[dim])/(Dims[dim]-2*NumberOfGhostZones);
    }





    int nBulge, nDisk, nHalo, nParticles;
    nBulge = nlines("bulge.dat");
    if(debug) fprintf(stderr, "InitializeParticles: Number of Bulge Particles %"ISYM"\n", nBulge);
    nDisk = nlines("disk.dat");
    if(debug) fprintf(stderr, "InitializeParticles: Number of Disk Particles %"ISYM"\n", nDisk);
    nHalo = nlines("halo.dat");
    if(debug) fprintf(stderr, "InitializeParticles: Number of Halo Particles %"ISYM"\n", nHalo);
    nParticles = nBulge + nDisk + nHalo;
    if(debug) fprintf(stderr, "InitializeParticles: Total Number of Particles %"ISYM"\n", nParticles);

    // Initialize particle arrays
    PINT *Number = new PINT[nParticles];
    int *Type = new int[nParticles];
    FLOAT *Position[MAX_DIMENSION];
    float *Velocity[MAX_DIMENSION];
    for (int i = 0; i < GridRank; i++)
    {
      Position[i] = new FLOAT[nParticles];
      Velocity[i] = new float[nParticles];
    }
    float *Mass = new float[nParticles];
    float *Attribute[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];
    for (int i = 0; i < NumberOfParticleAttributes; i++)
    {
      Attribute[i] = new float[nParticles];
      for (int j = 0; j < nParticles; j++)
	Attribute[i][j] = FLOAT_UNDEFINED;
    }

    FLOAT dx = CellWidth[0];
    // Read them in and assign them as we go
    int count = 0;
    ReadParticlesFromFile(
      Number, Type, Position, Velocity, Mass,
      "bulge.dat", PARTICLE_TYPE_STAR, count, dx,Center);
    ReadParticlesFromFile(
      Number, Type, Position, Velocity, Mass,
      "disk.dat", PARTICLE_TYPE_STAR, count, dx,Center);
    ReadParticlesFromFile(
      Number, Type, Position, Velocity, Mass,
      "halo.dat", PARTICLE_TYPE_DARK_MATTER, count, dx,Center);

    thisgrid->SetNumberOfParticles(count);
    thisgrid->SetParticlePointers(Mass, Number, Type, Position,
				  Velocity, Attribute);
    thisgrid->ExtraFunction("In Particle Creation");
    MetaData.NumberOfParticles = count;

    return SUCCESS;
}
int ReadParticlesFromFile(PINT *Number, int *Type, FLOAT *Position[],
			    float *Velocity[], float* Mass, const char* fname,
			    Eint32 particle_type, int &c, FLOAT dx, FLOAT * Center){
    FILE *fptr;
    char line[MAX_LINE_LENGTH];
    int ret;
    FLOAT x, y, z;
    float vx, vy, vz;
    double mass;

    float DensityUnits=1, LengthUnits=1, VelocityUnits=1, TimeUnits=1,
      TemperatureUnits=1;
    double MassUnits=1;

    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		 &TimeUnits, &VelocityUnits, &MassUnits, 0) == FAIL) {
      ENZO_FAIL("Error in GetUnits.");
    }

    fptr = fopen(fname, "r");

    while(fgets(line, MAX_LINE_LENGTH, fptr) != NULL)
    {
      ret +=
	sscanf(line,
	       "%"PSYM" %"PSYM" %"PSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
	       &x, &y, &z, &vx, &vy, &vz, &mass);

      Position[0][c] = x * kpc_cm / LengthUnits + Center[0];
      Position[1][c] = y * kpc_cm / LengthUnits + Center[1];
      Position[2][c] = z * kpc_cm / LengthUnits + Center[2];

      Velocity[0][c] = vx * km_cm / VelocityUnits;
      Velocity[1][c] = vy * km_cm / VelocityUnits;
      Velocity[2][c] = vz * km_cm / VelocityUnits;

      // Particle masses are actually densities.
      Mass[c] = mass * 1e9 * SolarMass / MassUnits / dx / dx / dx;
      Type[c] = particle_type;
      Number[c] = c++;
    }

    fclose(fptr);

    return c;
}

int nlines(const char* fname) {

  FILE* fptr = fopen(fname, "r");
  if(fptr==NULL){fprintf(stderr, "nlines error: %s doesn't exist!\n", fname); fflush(stderr);}
  int ch, n = 0;

  do
  {
    ch = fgetc(fptr);
    if(ch == '\n')
      n++;
  } while (ch != EOF);

  fclose(fptr);
  if (debug) fprintf(stderr,"Read %"ISYM" lines \n", n);
  return n;
}

void ReadInVcircData(FLOAT * VCircRadius, float * VCircVelocity)
{
FILE *fptr;
char line[MAX_LINE_LENGTH];
int i=0, ret;
float vcirc;
FLOAT rad;

fptr = fopen("vcirc.dat" , "r");

if(fptr==NULL){fprintf(stderr, "vcirc error: vcirc.dat doesn't exist!\n"); fflush(stderr);}

while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL)
{
  ret += sscanf(line, "%"PSYM" %"FSYM, &rad, &vcirc);
  VCircRadius[i] = rad*kpc_cm; // 3.08567758e21 = kpc/cm
  VCircVelocity[i] = vcirc*1e5; // 1e5 = (km/s)/(cm/s)
  i += 1;
}

fclose(fptr);
} // ReadInVcircData

float InterpolateVcircTable(FLOAT radius, FLOAT * VCircRadius, float * VCircVelocity)
{
int i;

for (i = 0; i < VCIRC_TABLE_LENGTH; i++)
  if (radius < VCircRadius[i])
break;

if (i == 0)
  return (VCircVelocity[i]) * (radius - VCircRadius[0]) / VCircRadius[0];
else if (i == VCIRC_TABLE_LENGTH)
  ENZO_FAIL("Fell off the circular velocity interpolation table");

// we know the radius is between i and i-1
return VCircVelocity[i-1] +
  (VCircVelocity[i] - VCircVelocity[i-1]) *
  (radius - VCircRadius[i-1])  /
  (VCircRadius[i] - VCircRadius[i-1]);
}
