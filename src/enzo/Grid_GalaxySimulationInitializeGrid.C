/***********************************************************************
[6~/
/  GRID CLASS (INITIALIZE THE GRID FOR A GALAXY SIMULATION)
/
/  written by: Greg Bryan
/  date:       May, 1998
/  modified1:  Elizabeth Tasker, Feb, 2004
/  modified1:  Elizabeth Tasker, Oct, 2006 (tidied up)
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <random>
#include <math.h>
#include <assert.h>
#include "preincludes.h" 
#include "EnzoTiming.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"
#include "phys_constants.h"

#define Mpc (3.0856e24)         //Mpc [cm] 
#define SolarMass (1.989e33)    //Solar Mass [g]
#define GravConst (6.67e-8)     //Gravitational Constant [cm3g-1s-2]
#define pi (3.14159)
#define mh (1.67e-24)           //Mass of Hydrogen [g]
#define kboltz (1.381e-16)      //Boltzmann's Constant [ergK-1]
#define kboltzKeV (8.617e-8)    // keV per K
#define mu (0.6)
#define CM_PER_KM (1.0e5)
#define CM_PER_KPC (3.0856e21)
#define KEV_PER_ERG (6.242e8)
#define VCIRC_TABLE_LENGTH 10000

int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);

int GetUnits(float *DensityUnits, float *LengthUnits,
       float *TemperatureUnits, float *TimeUnits,
       float *VelocityUnits, double *MassUnits, FLOAT Time);

int nlines(const char* fname);
int InitializeParticles(grid *thisgrid_orig, FLOAT *Center); //, HierarchyEntry &TopGrid, TopGridData &MetaData, FLOAT * Center);
int ReadParticlesFromFile(PINT *Number, int *Type, FLOAT *Position[],
			    float *Velocity[], float* Mass, const char* fname,
			    Eint32 particle_type, int &c, FLOAT dx, FLOAT * Center);
float InterpolateVcircTable(FLOAT radius, FLOAT * VCircRadius, float * VCircVelocity);
void ReadInVcircData(FLOAT * VCircRadius, float * VCircVelocity);


int CosmologyGetUnits(float *DensityUnits, float *LengthUnits,
                      float *TemperatureUnits, float *TimeUnits,
                      float *VelocityUnits, FLOAT Time);

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

int FindField(int field, int farray[], int numfields);

/* Internal routines */
void setup_chem(float density, float temperature, int equilibrate,
		float& DEdest, float&  HIdest, float& HIIdest,
		float& HeIdest, float& HeIIdest, float& HeIIIdest,
		float& HMdest, float& H2Idest, float& H2IIdest,
		float& DIdest, float& DIIdest, float& HDIdest);
float gasvel(FLOAT radius, float DiskDensity, FLOAT ExpansionFactor, 
             float GalaxyMass, FLOAT ScaleHeightR, FLOAT ScaleHeightz, 
             float DMConcentration, FLOAT Time);
float gauss_mass(FLOAT r, FLOAT z, FLOAT xpos, FLOAT ypos, FLOAT zpos, 
                 FLOAT inv [3][3], float DiskDensity,
                 FLOAT ScaleHeightR, FLOAT ScaleHeightz,
                 FLOAT cellwidth);
void rot_to_disk(FLOAT xpos, FLOAT ypos, FLOAT zpos, FLOAT &xrot,
                 FLOAT &yrot, FLOAT &zrot, FLOAT inv [3][3]);
double bilinear_interp(double x, double y, 
                       double x1, double x2, double y1, double y2,
                       double f_x1y1, double f_x1y2, 
                       double f_x2y1, double f_x2y2);


float NFWDarkMatterMassEnclosed(FLOAT R);
double MassEnclosed_r(FLOAT *binned_mass, double rad); 

/* struct to carry around data required for circumgalactic media
   if we need to generate radial profiles of halo quantities via 
   numerical integration */
struct CGMdata {
  double *n_rad, *T_rad, *rad, *press;
  int nbins;
  double R_inner, R_outer, dr;

  CGMdata(int n_bins) {
    
    nbins = n_bins;
  
    n_rad = new double[nbins];
    T_rad = new double[nbins];
    rad   = new double[nbins];

    for(int i=0; i<nbins; i++) n_rad[i] = T_rad[i] = rad[i] = -1.0;
  }
  
  ~CGMdata() {
    if (n_rad) delete[] n_rad;
    if (T_rad) delete[] T_rad;
    if (rad) delete[] rad;
  }
};

/* Internal Routines for CGM setup */
float HaloGasDensity(FLOAT R, struct CGMdata&);
float HaloGasTemperature(FLOAT R, struct CGMdata&);

/* Internal Routines for DiskGravity Setup */
double DiskGravityCircularVelocity(double rsph, double rcyl, double z);
double DiskGravityStellarAccel(double rcyl, double z);
double DiskGravityBulgeAccel(double rsph);

static float DensityUnits, LengthUnits, TemperatureUnits = 1,
             TimeUnits, VelocityUnits, MassUnits;

double gScaleHeightR, gScaleHeightz, densicm, MgasScale, Picm,
       TruncRadius, SmoothRadius, SmoothLength,Ticm;

/* Global variables (within this file) for circumgalactic medium setup 
   (also used a bit for disk potential setup) */
int GalaxySimulationGasHalo, EquilibrateChem;
double GalaxySimulationGasHaloScaleRadius,
  GalaxySimulationGasHaloDensity, GalaxySimulationGasHaloDensity2,
  GalaxySimulationGasHaloTemperature, GalaxySimulationGasHaloAlpha,
  GalaxySimulationGasHaloZeta, GalaxySimulationGasHaloZeta2,
  GalaxySimulationGasHaloCoreEntropy, 
  GalaxySimulationGalaxyMass, GalaxySimulationDMConcentration,
  GalaxySimulationGasHaloMetallicity,
  GalaxySimulationDiskMetallicityEnhancementFactor,
  GalaxySimulationGasHaloRatio;

/* declarations for a bunch of functions needed to generate radial profiles
   of halo quantities via numerical integration - see the actual functions for
   descriptions. */
double sigmoid(double r, double r0, double k, double y0, double y_off);
double halo_S_of_r(double r);
double halo_S_of_r(double r, grid* Grid, FLOAT *binned_mass);
double halo_dSdr(double r, double n);
double halo_dn_dr(double r, double n);
double halo_dP_dr(double r, double P, grid* Grid, FLOAT *binned_mass);
double sigmoid_halo_dP_dr(double r, double P, double log_r0,
			  double k, double y0, double y_off);
double halo_g_of_r(double r);
double halo_mod_g_of_r(double r, FLOAT *binned_mass);
double halo_mod_DMmass_at_r(double r);
double halo_mod_DMmass_at_r(FLOAT *binned_mass, double r);
void halo_init(struct CGMdata& CGM_data, grid* Grid, FLOAT *binned_mass, float Rstop=-1, int GasHalo_override=0);

float InterpolateVcircTable(FLOAT radius, FLOAT * VCircRadius, float * VCircVelocity);
int grid::GalaxySimulationInitializeGrida(FLOAT DiskRadius,
           FLOAT GalaxyMass,
           FLOAT GasMass,
           FLOAT DiskPosition[MAX_DIMENSION], 
           FLOAT ScaleHeightz,
           FLOAT ScaleHeightR,
           FLOAT GalaxyTruncationRadius,
           FLOAT DiskDensityCap, 
           FLOAT DMConcentration,
           FLOAT DiskTemperature,
           FLOAT InitialTemperature,
           FLOAT UniformDensity,
           int   EquilChem,
	   int   GasHalo,
           FLOAT GasHaloScaleRadius,
           FLOAT GasHaloDensity,
           FLOAT GasHaloDensity2,
           FLOAT GasHaloTemperature,
           FLOAT GasHaloAlpha,
           FLOAT GasHaloZeta,
           FLOAT GasHaloZeta2,
           FLOAT GasHaloCoreEntropy,
	   FLOAT GasHaloRatio,
           FLOAT GasHaloMetallicity,
           int   UseHaloRotation,
           FLOAT RotationScaleVelocity,
           FLOAT RotationScaleRadius,
           FLOAT RotationPowerLawIndex,
           FLOAT DiskMetallicityEnhancementFactor,
           FLOAT AngularMomentum[MAX_DIMENSION],
           FLOAT UniformVelocity[MAX_DIMENSION], 
           int UseMetallicityField, 
           FLOAT GalaxySimulationInflowTime,
           FLOAT GalaxySimulationInflowDensity,
           int level,
           FLOAT GalaxySimulationInitialBfield[MAX_DIMENSION],
           int GalaxySimulationInitialBfieldTopology,
            FLOAT  VCircRadius[], float  VCircVelocity[], 
	   FLOAT radius_bins[100], FLOAT binned_mass[100],
           FLOAT GalaxySimulationCR,
           int SetBaryons, 
	   int bin_mass
          )
{
  /* declarations */

  int dim, i, j, k, m, field, disk, size, MetalNum, MetalIaNum, vel;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum,
    H2IINum, DINum, DIINum, HDINum, B1Num, B2Num, B3Num, PhiNum;
  int CMNum, OMNum, CINum, OINum, OHINum, COINum, CHINum, CH2INum, C2INum, HCOINum, H2OINum, O2INum, CO_TOTALINum, H2O_TOTALINum, CIINum, OIINum, HOCIINum, HCOIINum, H3IINum, CHIINum, CH2IINum, COIINum, CH3IINum, OHIINum, H2OIINum, H3OIINum, O2IINum;  
  float DiskDensity, DiskVelocityMag;
  int CRNum, DensNum;
  int MassEnclosedNum;
  int isDiskNum, DebugNum; 
  int MassEnclosedSet = 0;  
  /* global-scope variables for disk potential functions (would be better if not global) */

  gScaleHeightR = ScaleHeightR;
  gScaleHeightz = ScaleHeightz;
  densicm = UniformDensity;  // gas density if no halo is used
  MgasScale = GasMass;
  Ticm = InitialTemperature;  // gas temperature if no halo is used
  Picm = kboltz*UniformDensity*Ticm/(mu*mh);  // gas pressure if no halo is used
  TruncRadius = GalaxyTruncationRadius;
  SmoothRadius = TruncRadius*.02/.026;
  SmoothLength = TruncRadius - SmoothRadius;

  /* set all of the gas halo quantities to variables that are global
     within this file, so the calls to HaloGasDensity() and HaloGasTemperature()
     are unchanged. */
  EquilibrateChem = EquilChem;
  GalaxySimulationGasHalo = GasHalo;   // integer, >= 0
  GalaxySimulationGasHaloScaleRadius = GasHaloScaleRadius;  // in mpc
  GalaxySimulationGasHaloDensity = GasHaloDensity; // in grams/cm^3
  GalaxySimulationGasHaloDensity2 = GasHaloDensity2;
  GalaxySimulationGasHaloTemperature = GasHaloTemperature;  // in Kelvin
  GalaxySimulationGasHaloAlpha = GasHaloAlpha;  // power-law index; unitless
  GalaxySimulationGasHaloZeta = GasHaloZeta;
  GalaxySimulationGasHaloZeta2 = GasHaloZeta2;
  GalaxySimulationGasHaloCoreEntropy = GasHaloCoreEntropy;  // power-law index; unitless
  GalaxySimulationGasHaloRatio = GasHaloRatio; // ratio of cooling time to freefall time
  GalaxySimulationGalaxyMass = GalaxyMass;
  GalaxySimulationDMConcentration = DMConcentration;
  GalaxySimulationGasHaloMetallicity = GasHaloMetallicity; // Zsun
  GalaxySimulationDiskMetallicityEnhancementFactor = DiskMetallicityEnhancementFactor; // w.r.t to halo
 
  /* create fields */
  NumberOfBaryonFields = 0;
  DensNum = NumberOfBaryonFields;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism)
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  if (WritePotential)
	  FieldType[NumberOfBaryonFields++] = GravPotential;
  if (WriteAcceleration){
	FieldType[NumberOfBaryonFields++] = Acceleration0;
	FieldType[NumberOfBaryonFields++] = Acceleration1;
	FieldType[NumberOfBaryonFields++] = Acceleration2;
  }
  vel = NumberOfBaryonFields;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  if (GridRank > 1) 
    FieldType[NumberOfBaryonFields++] = Velocity2;
  if (GridRank > 2)
    FieldType[NumberOfBaryonFields++] = Velocity3;
  if (UseMHD) {
    FieldType[B1Num = NumberOfBaryonFields++] = Bfield1;
    FieldType[B2Num = NumberOfBaryonFields++] = Bfield2;
    FieldType[B3Num = NumberOfBaryonFields++] = Bfield3;
  }
  if(HydroMethod == MHD_RK ){
    FieldType[PhiNum = NumberOfBaryonFields++] = PhiField;
  }
  if (UsePoissonDivergenceCleaning) {
    printf("Using Div Cleaning");
    FieldType[NumberOfBaryonFields++] = Phi_pField;
  }

  /* If cosmic rays present, set up field */
  CRNum = NumberOfBaryonFields;
  if( CRModel )
    FieldType[NumberOfBaryonFields++] = CRDensity;

  if (MultiSpecies) {
    FieldType[DeNum    = NumberOfBaryonFields++] = ElectronDensity;
    FieldType[HINum    = NumberOfBaryonFields++] = HIDensity;
    FieldType[HIINum   = NumberOfBaryonFields++] = HIIDensity;
    FieldType[HeINum   = NumberOfBaryonFields++] = HeIDensity;
    FieldType[HeIINum  = NumberOfBaryonFields++] = HeIIDensity;
    FieldType[HeIIINum = NumberOfBaryonFields++] = HeIIIDensity;
    if (MultiSpecies > 1) {
      FieldType[HMNum    = NumberOfBaryonFields++] = HMDensity;
      FieldType[H2INum   = NumberOfBaryonFields++] = H2IDensity;
      FieldType[H2IINum  = NumberOfBaryonFields++] = H2IIDensity;
    }
    if (MultiSpecies > 2) {
      FieldType[DINum   = NumberOfBaryonFields++] = DIDensity;
      FieldType[DIINum  = NumberOfBaryonFields++] = DIIDensity;
      FieldType[HDINum  = NumberOfBaryonFields++] = HDIDensity;
    }
    
    FieldType[CMNum = NumberOfBaryonFields++] = CMDensity;
    FieldType[OMNum = NumberOfBaryonFields++] = OMDensity;
    FieldType[CINum = NumberOfBaryonFields++] = CIDensity; 
    FieldType[OINum = NumberOfBaryonFields++] = OIDensity; 
    FieldType[OHINum = NumberOfBaryonFields++] = OHIDensity; 
    FieldType[COINum = NumberOfBaryonFields++] = COIDensity; 
    FieldType[CHINum = NumberOfBaryonFields++] = CHIDensity; 
    FieldType[CH2INum = NumberOfBaryonFields++] = CH2IDensity; 
    FieldType[C2INum = NumberOfBaryonFields++] = C2IDensity; 
    FieldType[HCOINum = NumberOfBaryonFields++] = HCOIDensity; 
    FieldType[H2OINum = NumberOfBaryonFields++] = H2OIDensity;
    FieldType[O2INum = NumberOfBaryonFields++] = O2IDensity;
    FieldType[CO_TOTALINum = NumberOfBaryonFields++] = CO_TOTALIDensity; 
    FieldType[H2O_TOTALINum = NumberOfBaryonFields++] = H2O_TOTALIDensity;
    FieldType[CIINum = NumberOfBaryonFields++] = CIIDensity; 
    FieldType[OIINum = NumberOfBaryonFields++] = OIIDensity; 
    FieldType[HOCIINum = NumberOfBaryonFields++] = HOCIIDensity;
    FieldType[HCOIINum = NumberOfBaryonFields++] = HCOIIDensity;
    FieldType[H3IINum = NumberOfBaryonFields++] = H3IIDensity; 
    FieldType[CHIINum = NumberOfBaryonFields++] = CHIIDensity; 
    FieldType[CH2IINum = NumberOfBaryonFields++] = CH2IIDensity; 
    FieldType[COIINum = NumberOfBaryonFields++] = COIIDensity;
    FieldType[CH3IINum = NumberOfBaryonFields++] = CH3IIDensity; 
    FieldType[OHIINum = NumberOfBaryonFields++] = OHIIDensity; 
    FieldType[H2OIINum = NumberOfBaryonFields++] = H2OIIDensity; 
    FieldType[H3OIINum = NumberOfBaryonFields++] = H3OIIDensity; 
    FieldType[O2IINum = NumberOfBaryonFields++] = O2IIDensity;
  }
    FieldType[MassEnclosedNum = NumberOfBaryonFields++] = MassEnclosed;
    FieldType[isDiskNum = NumberOfBaryonFields++] = isDisk; 
    FieldType[DebugNum = NumberOfBaryonFields++] = DebugField;
  if (UseMetallicityField)
    FieldType[MetalNum = NumberOfBaryonFields++] = Metallicity; /* fake it with metals */
  if (StarMakerTypeIaSNe)
    FieldType[MetalIaNum = NumberOfBaryonFields++] = MetalSNIaDensity;

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber){ 
   	std::cout << "Processor leaving " << ProcessorNumber << std::endl; 
	return SUCCESS;
  }
  if(!SetBaryons)
    return SUCCESS; //leave if baryon flag not on
  /* Set various units. */

  float CriticalDensity = 1, BoxLength = 1;
  FLOAT a, dadt, ExpansionFactor = 1;
  if (ComovingCoordinates) {
    CosmologyComputeExpansionFactor(Time, &a, &dadt);
    ExpansionFactor = a/(1.0+InitialRedshift);
    CosmologyGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		      &TimeUnits, &VelocityUnits, Time);
    CriticalDensity = 2.78e11*POW(HubbleConstantNow, 2); // in Msolar/Mpc^3
    BoxLength = ComovingBoxSize*ExpansionFactor/HubbleConstantNow;  // in Mpc
  } else if( PointSourceGravity ){
    ENZO_FAIL("ERROR IN GALAXY SIM GRID INITIALIZE: non-cosmology units not supported for point source gravity");
  } else {
    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
                 &TimeUnits, &VelocityUnits, Time) == FAIL) {
      ENZO_FAIL("Error in GetUnits.");
    } // end get units error if  
  } // end units if/else
 
  /* correct background density if it's not given in code units */
  if( UniformDensity < 1.0E-10 ){
    UniformDensity /= DensityUnits;
    if( debug && MyProcessorNumber == ROOT_PROCESSOR ) 
      fprintf(stdout,"Converting GalaxySimulationUniformDensity = %"GSYM" from CGS to code units\n",UniformDensity);
  } // end uniform density if

  /* Set up inflow */
  if (GalaxySimulationInflowTime > 0.0){
    TimeActionType[0] = 2;
    TimeActionParameter[0] = GalaxySimulationInflowDensity*DensityUnits;
    TimeActionTime[0] = GalaxySimulationInflowTime*1e9/TimeUnits;
  }

  /* Scale gas halo rotation quantities to code units.
   * gas halo rotation variable are NOT global */
  RotationScaleVelocity *= CM_PER_KM; // km/s to cm/s
  RotationScaleVelocity /= LengthUnits/TimeUnits; // cm/s to code length/code time
  RotationScaleRadius *= CM_PER_KPC;  // kpc to cm
  RotationScaleRadius /= LengthUnits;  // cm to code length

  /*  initializes halo radius, density, temperature profiles 
      for circumgalactic medium if needed (i.e., for CGM profiles that
      require integration to get quantities we care about. 
      Assumes disk is located at the center of the domain. */

  double far_left, far_right, largest_rad;
  
  far_left = DomainLeftEdge[0]-DiskPosition[0];
  far_right = DomainRightEdge[0]-DiskPosition[0];
  
  for (int i=1; i<GridRank; ++i) {
    if (DomainLeftEdge[i] - DiskPosition[i]< far_left)
      far_left = DomainLeftEdge[i]-DiskPosition[i];
    if (DomainRightEdge[i] - DiskPosition[i] > far_right)
      far_right = DomainRightEdge[i]-DiskPosition[i];
  }

  largest_rad = sqrt(3) * (far_right - far_left) / 2.0 * LengthUnits;
  /* compute size of fields */
  size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* allocate fields */
  this->AllocateGrids();

  /* I'm commenting this out because the metal field should
     be set during grid initialization rather than just setting
     it as a constant color field. -- DWS */
  // /* set metals to small value */
  //  if (UseMetallicityField)
  //    for (i = 0; i < size; i++)
  //      BaryonField[MetalNum][i] = 1.0e-10;
 
  /* Loop over the mesh. */
  float density, disk_dens;
  FLOAT halo_vmag, disk_vel[MAX_DIMENSION], Velocity[MAX_DIMENSION];
  FLOAT temperature, disk_temp, init_temp, initial_metallicity;
  FLOAT r_sph, x, y = 0, z = 0;
  FLOAT this_x, this_y, this_z; 
  int n = 0, iter;
  double mean = 0.0; 
  double stdev = 1.0; 
  std::default_random_engine gen; 
  std::normal_distribution<double> distribution(mean,stdev);


  //The setting of density currently "works" in 3 passes: 
  //1. Set Disk and UniformDensity into cells 
  //2. Set down mass from particles 
  //3. Set down gas halo after calculating enclosed mass from previous steps


  //first pass over grids: this sets density from particles, disk, and UniformDensity background
  for (k = 0; k < GridDimension[2]; k++)
    for (j = 0; j < GridDimension[1]; j++)
      for (i = 0; i < GridDimension[0]; i++, n++) {
	      density = UniformDensity;
	      BaryonField[isDiskNum][n] = 0; 
	if (UseMetallicityField) {
	  /* Set a background metallicity value that will scale with density.
	     If the cell is in the disk, this wifll be increased by a factor
	     of 3.  This should really be a parameter that is read in -- DWS */ 
	  initial_metallicity = GalaxySimulationGasHaloMetallicity;
	}

	/* Compute position */

	x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - DiskPosition[0];
	if (GridRank > 1)
	  y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - DiskPosition[1];
	if (GridRank > 2)
	  z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k] - DiskPosition[2];
        double max_pos = 0.5; //this assumes that the disk is at the same coord in x,y,z. change if it's not. 	
	double min_pos = -0.5; 
	/*if(x > max_pos)
		x = 1 - x; 
	if(x < min_pos)
		x = 1 + x; 
	if (GridRank > 1){
	  if(y > max_pos)
		  y = 1 - y; 
	  if(y < min_pos)
		  y = 1 + y; 
	}
	if (GridRank > 2){
	  if(z > max_pos) 
		  z = 1 - z; 
	  if(z < min_pos)
		  z = 1 + z; 
	}*/

	for (dim = 0; dim < MAX_DIMENSION; dim++){
	  Velocity[dim] = 0;
	  disk_vel[dim] = 0;
  }

	/* Find distance from center. */
	r_sph = sqrt(POW(fabs(x), 2) +
		     POW(fabs(y), 2) +
		     POW(fabs(z), 2) );
	r_sph = max(r_sph, 0.1*CellWidth[0][0]);
    
	
	  //r_cyl = sqrt(POW(fabs(x-DiskPosition[0]), 2) +
	  //POW(fabs(y-DiskPosition[1]), 2) );
	

	FLOAT xpos, ypos, zpos, rsph, zheight, rcyl, theta; 
	float CellMass;
	FLOAT rp_hat[3];
	FLOAT yhat[3];

	/* Loop over dims if using Zeus (since vel's face-centered). */

	for (dim = 0; dim < 1+(HydroMethod == Zeus_Hydro ? GridRank : 0);
	     dim++) {

	  /* Compute position. */
	  xpos = x-(dim == 1 ? 0.5*CellWidth[0][0] : 0.0);
	  ypos = y-(dim == 2 ? 0.5*CellWidth[1][0] : 0.0);
	  zpos = z-(dim == 3 ? 0.5*CellWidth[2][0] : 0.0);

	  /* Compute z and r_perp (AngularMomentum is angular momentum 
	     and must have unit length). */    

	  /* magnitude of z = r.L in L direction */

	  zheight = AngularMomentum[0]*xpos + 
	    AngularMomentum[1]*ypos +
	    AngularMomentum[2]*zpos;

	  /* position in plane of disk */

	  rp_hat[0] = xpos - zheight*AngularMomentum[0];
	  rp_hat[1] = ypos - zheight*AngularMomentum[1];
	  rp_hat[2] = zpos - zheight*AngularMomentum[2];
	  rcyl = sqrt(rp_hat[0]*rp_hat[0] + rp_hat[1]*rp_hat[1] + rp_hat[2]*rp_hat[2]);

	  /* Normalize the vector r_perp = unit vector pointing along plane of disk */

	  rp_hat[0] = rp_hat[0]/rcyl;
	  rp_hat[1] = rp_hat[1]/rcyl;
	  rp_hat[2] = rp_hat[2]/rcyl;

	  /* If requested, calculate velocity for CGM halo.
	   * Will be replaced wtih disk velocity later if appropriate */
	  if (UseHaloRotation){
	    /* polar angle as measured from the angular momentum vector*/
	    theta = acos(zheight/r_sph);

	    halo_vmag = RotationScaleVelocity // code units
	      * POW(r_sph/RotationScaleRadius, 
		    RotationPowerLawIndex);

	    if (r_sph <= RotationScaleRadius)
	      halo_vmag = RotationScaleVelocity;

	    halo_vmag *= sin(theta)*sin(theta);
	  
	    /* Cylindrical velocity */
	    Velocity[0] = halo_vmag * (AngularMomentum[1]*rp_hat[2] -
				       AngularMomentum[2]*rp_hat[1]);
	    Velocity[1] = halo_vmag * (AngularMomentum[2]*rp_hat[0] -
				       AngularMomentum[0]*rp_hat[2]);
	    Velocity[2] = halo_vmag * (AngularMomentum[0]*rp_hat[1] -
				       AngularMomentum[1]*rp_hat[0]);
	  }

	  disk_dens = 0.0;
	  if (r_sph < DiskRadius) {

	    /* Beyond truncation radius */
	    if( fabs(rcyl*LengthUnits/Mpc) > TruncRadius ){
	      disk_dens = 0.0;
	      break;
	    }

	    /* Find another vector perpendicular to r_perp and AngularMomentum */

	    yhat[0] = AngularMomentum[1]*rp_hat[2] - AngularMomentum[2]*rp_hat[1];
	    yhat[1] = AngularMomentum[2]*rp_hat[0] - AngularMomentum[0]*rp_hat[2];
	    yhat[2] = AngularMomentum[0]*rp_hat[1] - AngularMomentum[1]*rp_hat[0];

	    /* generate rotation matrix */
	    FLOAT inv[3][3],temp;
	    int i,j;

	    // matrix of basis vectors in coordinate system defined by the galaxy
	    inv[0][0] = rp_hat[0];
	    inv[0][1] = yhat[0];
	    inv[0][2] = AngularMomentum[0];
        
	    inv[1][0] = rp_hat[1];
	    inv[1][1] = yhat[1];
	    inv[1][2] = AngularMomentum[1];
        
	    inv[2][0] = rp_hat[2];
	    inv[2][1] = yhat[2];
	    inv[2][2] = AngularMomentum[2];

	    // Matrix is orthogonal by construction so inverse = transpose
	    for (i=0;i<3;i++)
	      for (j=i+1;j<3;j++){
          temp = inv[i][j];
          inv[i][j] = inv[j][i];
          inv[j][i] = temp;
	      }

	    DiskDensity = (GasMass * SolarMass
			   / (8.0*pi*ScaleHeightz*Mpc*POW(ScaleHeightR*Mpc,2.0)))
	      / DensityUnits;   //Code units (rho_0) 
	    CellMass = gauss_mass(rcyl*LengthUnits, zheight*LengthUnits,
				  xpos*LengthUnits, ypos*LengthUnits,
				  zpos*LengthUnits, inv, 
				  DiskDensity*DensityUnits,
				  ScaleHeightR*Mpc, ScaleHeightz*Mpc, 
				  CellWidth[0][0]*LengthUnits);

	    disk_dens = CellMass/POW(CellWidth[0][0]*LengthUnits,3)/DensityUnits;

	    if ((disk_dens > DiskDensityCap) && (DiskDensityCap > 0))
	      disk_dens = DiskDensityCap;

	    /* Inside CGM */
	    if (disk_dens < UniformDensity){ //was just density, now the background UniformDensity bc the gas halo isn't laid down yet 
		    break;
	    }

	    //
	    // calculate velocity
	    //
	    if (PointSourceGravity > 0 )
	      DiskVelocityMag = gasvel(rcyl, DiskDensity, ExpansionFactor,
				       GalaxyMass, ScaleHeightR,
				       ScaleHeightz, DMConcentration, Time);

	    else if( DiskGravity > 0 )
	      DiskVelocityMag = DiskGravityCircularVelocity(r_sph*LengthUnits,
							    rcyl*LengthUnits,
							    zheight*LengthUnits)
		    /VelocityUnits;
        else{
            DiskVelocityMag = InterpolateVcircTable(r_sph*LengthUnits, VCircRadius, VCircVelocity)/VelocityUnits; 
        }
        
	    if (PointSourceGravity*DiskGravity != FALSE ) 
	      ENZO_FAIL("Cannot activate both PointSource and Disk gravity options for Isolated Galaxy");

	    /* Compute velocty: L x r_perp. */
	    if (dim == 0 || dim == 1)
	      disk_vel[0] = DiskVelocityMag*(AngularMomentum[1]*rp_hat[2] -
					     AngularMomentum[2]*rp_hat[1]);
	    if (dim == 0 || dim == 2)
	      disk_vel[1] = DiskVelocityMag*(AngularMomentum[2]*rp_hat[0] -
					     AngularMomentum[0]*rp_hat[2]);
	    if (dim == 0 || dim == 3)
	      disk_vel[2] = DiskVelocityMag*(AngularMomentum[0]*rp_hat[1] -
					     AngularMomentum[1]*rp_hat[0]);

	  } // end: if (r_sph < DiskRadius)

	  /* Replace CGM ("Halo") defaults with disk if dense enough; i.e.
	   * replace 'density', 'temperature', 'initial_metallicity', and
	   * 'Velocity' (which are currently set to CGM values) with their
	   * appropriate disk values */
	  if (disk_dens > density && fabs(rcyl*LengthUnits/Mpc) <= TruncRadius){
	    density = disk_dens;
	    temperature = DiskTemperature;
        
	    /* Here we're setting the disk to be X times more enriched -- DWS */
	    if( UseMetallicityField )
	      initial_metallicity *= GalaxySimulationDiskMetallicityEnhancementFactor;
          
	    /* Replace default/CGM velocity with disk velocity */
	    Velocity[0] = disk_vel[0];
	    Velocity[1] = disk_vel[1];
	    Velocity[2] = disk_vel[2];
	    BaryonField[isDiskNum][n] = 1; //flag we are in the disk for InitializeGridb 
	  }

	} // end: loop over dims 

	/* Set density. */
	BaryonField[0][n] = density;
	
	//set velocity 
	BaryonField[vel][n] = Velocity[0]; 
	BaryonField[vel+1][n] = Velocity[1]; 
	BaryonField[vel+2][n] = Velocity[2]; 

	if(bin_mass){
	double delta_r = 1.0 / 100; //code units
	int radius_bin = int(r_sph / delta_r);	
	binned_mass[radius_bin] += density * CellWidth[0][0] * CellWidth[1][0] * CellWidth[2][0]; }
      } //end: first pass over grid
	
  //if on the root grid, init particles now 
  if ( DiskGravity + PointSourceGravity == FALSE && level == 0){
      InitializeParticles(this, DiskPosition);
  }
  
  //now deposit particle masses
  if(level==0 && bin_mass){
  double MassUnitsDouble = double(DensityUnits)*POW(double(LengthUnits), 3.0);
  n = 0; 
  density = 0;
  double total_dm_mass = 0.0; //for debug, variable to check how much mass is in the DM particles 
  double total_init_star_mass = 0.0; //same, but for the initial star particles  
  int total_dm_particles = 0; 
  for(int p = 0; p < this->NumberOfParticles; p++){
  	int ind; 
	int x,y,z; 
	if(ParticleType[p] == PARTICLE_TYPE_DARK_MATTER){
		total_dm_mass += ParticleMass[p]*CellWidth[0][0]*CellWidth[1][0]*CellWidth[2][0]/SolarMass*MassUnitsDouble; 
		total_dm_particles++; 
	}
	if(ParticleType[p] != PARTICLE_TYPE_DARK_MATTER)
		total_init_star_mass += ParticleMass[p]*CellWidth[0][0]*CellWidth[1][0]*CellWidth[2][0]/SolarMass*MassUnitsDouble; 
	x = int((ParticlePosition[0][p] - GridLeftEdge[0]) / CellWidth[0][0]) + NumberOfGhostZones; 
	y = int((ParticlePosition[1][p] - GridLeftEdge[1]) / CellWidth[1][0]) + NumberOfGhostZones; 
	z = int((ParticlePosition[2][p] - GridLeftEdge[2]) / CellWidth[2][0]) + NumberOfGhostZones;
        if((ParticlePosition[0][p] < GridLeftEdge[0]) || (ParticlePosition[1][p] < GridLeftEdge[1]) || (ParticlePosition[2][p] < GridLeftEdge[2]) || (ParticlePosition[0][p] > GridRightEdge[0]) || (ParticlePosition[1][p] > GridRightEdge[1]) || (ParticlePosition[2][p] > GridRightEdge[2])){
		std::cout << "no" << std::endl;
		//continue;
	}
	ind = (z * GridDimension[1] + y) * GridDimension[0] + x; 
	density = ParticleMass[p]*CellWidth[0][0]*CellWidth[1][0]*CellWidth[2][0]; //code mass  
	float xpos, ypos, zpos;
    	xpos = CellLeftEdge[0][x] + 0.5*CellWidth[0][x] - DiskPosition[0];
	ypos = CellLeftEdge[1][y] + 0.5*CellWidth[1][y] - DiskPosition[1];
	zpos = CellLeftEdge[2][z] + 0.5*CellWidth[2][z] - DiskPosition[2];
        double max_pos = 0.5; //this assumes that the disk is at the same coord in x,y,z. change if it's not. 	
	double min_pos = -0.5; 
	if(xpos > max_pos)
		xpos = 1 - xpos; 
	if(xpos < min_pos)
		xpos = 1 + xpos; 
	if (GridRank > 1){
	  if(ypos > max_pos)
		  ypos = 1 - ypos; 
	  if(ypos < min_pos)
		  ypos = 1 + ypos; 
	}
	if (GridRank > 2){
	  if(zpos > max_pos) 
		  zpos = 1 - zpos; 
	  if(zpos < min_pos)
		  zpos = 1 + zpos; 
	}
    	r_sph = sqrt(POW(fabs(xpos), 2) +
		     POW(fabs(ypos), 2) +
		     POW(fabs(zpos), 2) );
	r_sph = max(r_sph, 0.1*CellWidth[0][0]);
  	BaryonField[MassEnclosedNum][ind] = ParticlePosition[0][p];
	double delta_r = 1.0 / 100; //code units
    	int ind_r = int(r_sph / delta_r);
	binned_mass[ind_r] += density; 	
  	}
  std::cout << "total num dm particles: " << total_dm_particles << std::endl;
  std::cout << "total dark matter mass: " << total_dm_mass << std::endl;
  std::cout << "total init star mass: " << total_init_star_mass << std::endl;
  }

 return SUCCESS; 
} 

int grid::GalaxySimulationInitializeGridb(FLOAT DiskRadius,
           FLOAT GalaxyMass,
           FLOAT GasMass,
           FLOAT DiskPosition[MAX_DIMENSION], 
           FLOAT ScaleHeightz,
           FLOAT ScaleHeightR,
           FLOAT GalaxyTruncationRadius,
           FLOAT DiskDensityCap, 
           FLOAT DMConcentration,
           FLOAT DiskTemperature,
           FLOAT InitialTemperature,
           FLOAT UniformDensity,
           int   EquilChem,
	   int   GasHalo,
           FLOAT GasHaloScaleRadius,
           FLOAT GasHaloDensity,
           FLOAT GasHaloDensity2,
           FLOAT GasHaloTemperature,
           FLOAT GasHaloAlpha,
           FLOAT GasHaloZeta,
           FLOAT GasHaloZeta2,
           FLOAT GasHaloCoreEntropy,
	   FLOAT GasHaloRatio,
           FLOAT GasHaloMetallicity,
           int   UseHaloRotation,
           FLOAT RotationScaleVelocity,
           FLOAT RotationScaleRadius,
           FLOAT RotationPowerLawIndex,
           FLOAT DiskMetallicityEnhancementFactor,
           FLOAT AngularMomentum[MAX_DIMENSION],
           FLOAT UniformVelocity[MAX_DIMENSION], 
           int UseMetallicityField, 
           FLOAT GalaxySimulationInflowTime,
           FLOAT GalaxySimulationInflowDensity,
           int level,
           FLOAT GalaxySimulationInitialBfield[MAX_DIMENSION],
           int GalaxySimulationInitialBfieldTopology,
            FLOAT  VCircRadius[], float  VCircVelocity[],
	   FLOAT radius_bins[100], FLOAT binned_mass[100],
           FLOAT GalaxySimulationCR,
           int SetBaryons, 
	   int bin_mass 
          )
{
  /* declarations */
  int dim, i, j, k, m, field, disk, size, MetalNum, MetalIaNum, vel;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum,
    H2IINum, DINum, DIINum, HDINum, B1Num, B2Num, B3Num, PhiNum;
  int CMNum, OMNum, CINum, OINum, OHINum, COINum, CHINum, CH2INum, C2INum, HCOINum, H2OINum, O2INum, CO_TOTALINum, H2O_TOTALINum, CIINum, OIINum, HOCIINum, HCOIINum, H3IINum, CHIINum, CH2IINum, COIINum, CH3IINum, OHIINum, H2OIINum, H3OIINum, O2IINum;  
  int CRNum, DensNum;
  int TENum, GENum, Vel1Num, Vel2Num, Vel3Num;  
  int MassEnclosedNum, isDiskNum, DebugNum;
  float density, disk_dens;
  FLOAT halo_vmag, disk_vel[MAX_DIMENSION], Velocity[MAX_DIMENSION];
  FLOAT temperature, disk_temp, init_temp, initial_metallicity;
  FLOAT r_sph, x, y = 0, z = 0;
  double mean = 0.0; 
  double stdev = 1.0; 
  std::default_random_engine gen; 
  std::normal_distribution<double> distribution(mean,stdev);
  
  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber){ 
   	std::cout << "Processor leaving " << ProcessorNumber << std::endl; 
	return SUCCESS;
  }
  if(!SetBaryons)
    return SUCCESS; //leave if baryon flag not on
  
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
  				       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }
  if (CRModel)
    if ((CRNum = FindField(CRDensity, FieldType, NumberOfBaryonFields)) < 0)
      ENZO_FAIL("Cannot Find Cosmic Rays");
  MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields);
  MetalIaNum = FindField(MetalSNIaDensity, FieldType, NumberOfBaryonFields); 
  DeNum = FindField(ElectronDensity, FieldType, NumberOfBaryonFields); 
  HINum = FindField(HIDensity, FieldType, NumberOfBaryonFields); 
  HIINum = FindField(HIIDensity, FieldType, NumberOfBaryonFields);
  HeINum = FindField(HeIDensity, FieldType, NumberOfBaryonFields); 
  HeIINum = FindField(HeIIDensity, FieldType, NumberOfBaryonFields); 
  HeIIINum = FindField(HeIIIDensity, FieldType, NumberOfBaryonFields);
  HMNum = FindField(HMDensity, FieldType, NumberOfBaryonFields); 
  H2INum = FindField(H2IDensity, FieldType, NumberOfBaryonFields); 
  H2IINum = FindField(H2IIDensity, FieldType, NumberOfBaryonFields); 
  DINum = FindField(DIDensity, FieldType, NumberOfBaryonFields); 
  DIINum = FindField(DIIDensity, FieldType, NumberOfBaryonFields); 
  HDINum = FindField(HDIDensity, FieldType, NumberOfBaryonFields);
  PhiNum = FindField(PhiField, FieldType, NumberOfBaryonFields); 
  CMNum = FindField(CMDensity, FieldType, NumberOfBaryonFields); 
  OMNum = FindField(OMDensity, FieldType, NumberOfBaryonFields); 
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
  MassEnclosedNum = FindField(MassEnclosed, FieldType, NumberOfBaryonFields);   
  isDiskNum = FindField(isDisk, FieldType, NumberOfBaryonFields);   
  DebugNum = FindField(DebugField, FieldType, NumberOfBaryonFields); 
  //halo init needs to get called after all the mass is placed
  //so that we can tell what the actual m_enclosed is at a given r 
  int n = 0;
  double far_left, far_right, largest_rad;
  
  largest_rad = sqrt(3*0.25) * LengthUnits;//periodic wrap means this should be true. unless you changed the box size in code units. 
  struct CGMdata CGM_data(8162);
  halo_init(CGM_data, this, binned_mass, -1.0);
  double MassUnitsDouble = double(DensityUnits)*POW(double(LengthUnits), 3.0);
  float rho_crit = 1.8788e-29*0.49; 
  float R200 = pow(3.0/(4.0*3.14159) * binned_mass[99]*MassUnitsDouble/(200.0*rho_crit), 1./3.);
  //begin second pass over grid
  for (k = 0; k < GridDimension[2]; k++)
    for (j = 0; j < GridDimension[1]; j++)
      for (i = 0; i < GridDimension[0]; i++, n++) {
	
	  if (UseMetallicityField) {
	  /* Set a background metallicity value that will scale with density.
	     If the cell is in the disk, this wifll be increased by a factor
	     of 3.  This should really be a parameter that is read in -- DWS */ 
	  initial_metallicity = GalaxySimulationGasHaloMetallicity;
	}

	/* Compute position */

	x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	if (GridRank > 1)
	  y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	if (GridRank > 2)
	  z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
	
	double xpos = x - DiskPosition[0];
        double max_pos = 0.5; //this assumes that the disk is at the same coord in x,y,z. change if it's not. 	
	double min_pos = -0.5; 
	if(xpos > max_pos)
		xpos = 1 - xpos; 
	if(xpos < min_pos) 
		xpos = 1 + xpos; 
	double ypos = y - DiskPosition[1]; 
	if(ypos > max_pos)
		ypos = 1 - ypos; 
	if(ypos < min_pos) 
		ypos = 1 + ypos; 
	double zpos = z - DiskPosition[2]; 
	if(zpos > max_pos)
		zpos = 1 - zpos; 
	if(zpos < min_pos) 
		zpos = 1 + zpos; 
	
	/* Find distance from center. */

	r_sph = sqrt(POW(fabs(xpos), 2) +
		     POW(fabs(ypos), 2) +
		     POW(fabs(zpos), 2) );
	r_sph = max(r_sph, 0.1*CellWidth[0][0]);
	density = 0.0; 
	//add halo if we're in the right spot
	if(r_sph*LengthUnits < R200){
		temperature = disk_temp = init_temp = HaloGasTemperature(r_sph, CGM_data);
		density += HaloGasDensity(r_sph, CGM_data)/DensityUnits; 
	}
	else{ 
		temperature = InitialTemperature; //background box temp if outside the CGM   
	}
	if (BaryonField[isDiskNum][n]){
	    temperature = DiskTemperature;//set disk temp as specified by parameter file
	}
	if (UseMetallicityField) {
	  BaryonField[MetalNum][n] = initial_metallicity 
	    * CoolData.SolarMetalFractionByMass 
	    * density;
	}

	/* This should probably be scaled with density in some way to be
	   a proper metallicity -- DWS (loop redundancy addressed by CEK) */
	if (StarMakerTypeIaSNe)
	  BaryonField[MetalIaNum][n] = 1.0e-10;
   

    //new chem densities; init to 0 they get populated with time evolution for now 
    if(CRModel)
    	BaryonField[CMNum][n] = 0.0*density; 
    if(MultiSpecies){
    BaryonField[OMNum][n] = 0.0*density; 
    BaryonField[HCOINum][n] = 0.0; 
    BaryonField[CO_TOTALINum][n] = 0.0; 
    BaryonField[H2O_TOTALINum][n] = 0.0; 
    BaryonField[HOCIINum][n] = 0.0; 
    BaryonField[H3IINum][n] = 0.0; 
    BaryonField[CHIINum][n] = 0.0; 
    BaryonField[CH2IINum][n] = 0.0; 
    BaryonField[COIINum][n] = 0.0; 
    BaryonField[OHIINum][n] = 0.0; 
    BaryonField[H2OIINum][n] = 0.0; 
    BaryonField[H3OIINum][n] = 0.0; 
    BaryonField[O2IINum][n] = 0.0; 
    }
    for(dim = 0; dim < 3; dim++){
	BaryonField[Vel1Num][n] += UniformVelocity[dim];
	BaryonField[Vel2Num][n] += UniformVelocity[dim];
	BaryonField[Vel3Num][n] += UniformVelocity[dim];
    }

    //add the halo density, it's already in code units
    	BaryonField[0][n] += density;

    /* Set energy (thermal and then total if necessary). */

	BaryonField[1][n] = temperature/TemperatureUnits / ((Gamma-1.0)*mu);

	if (DualEnergyFormalism)
	  BaryonField[2][n] = BaryonField[1][n];

	if (HydroMethod != Zeus_Hydro)
	  for (dim = 0; dim < GridRank; dim++)
	    BaryonField[1][n] += 0.5*POW(BaryonField[Vel1Num+dim][n], 2);

	if (BaryonField[1][n] <= 0.0)
	  printf("G_GSIC: negative or zero energy  n = %"ISYM"  temp = %"FSYM"   e = %"FSYM"\n",
		 n, temperature, BaryonField[1][n]);

	if ( UseMHD ){
	  switch ( GalaxySimulationInitialBfieldTopology ){
	  case 0: //uniform
	    for (dim = 0; dim < GridRank; dim++) {
	      if( UseMHDCT ){
		MagneticField[dim][n] = GalaxySimulationInitialBfield[dim];
	      }
	      BaryonField[B1Num+dim][n] = GalaxySimulationInitialBfield[dim];
	    }
	    break;
	  case 1: //white noise vector potential 
		{	
		//generate random values for vector potential
		//store vector potential values in ElectricField. Bc Grid_MHD_Curl does it, that's why. 
		for(int ind = 0; ind < GridRank; ind++){
			ElectricField[ind][n] = distribution(gen)*CellWidth[0][0]*1e-9; //the magnetic field has a cellwidth dependence from the curl, so undo it so the potential isn't level dependent 
			}
	  	}
	      	break;
          default:
	    ENZO_FAIL("undefined value of GalaxySimulationInitialBfieldTopology");
	  }
	  BaryonField[1][n] += 0.5*(BaryonField[B1Num][n]*BaryonField[B1Num][n]
				    +BaryonField[B2Num][n]*BaryonField[B2Num][n]
				    +BaryonField[B3Num][n]*BaryonField[B3Num][n])/
	    BaryonField[0][n];
	}//UseMHD
	if( CRModel )
	  BaryonField[CRNum][n] = BaryonField[DensNum][n] * GalaxySimulationCR;

	// Set multispecies fields!
	// this attempts to set them such that species conservation is maintained,
	// using the method in CosmologySimulationInitializeGrid.C
	if(MultiSpecies){
	  if (MultiSpecies == 3)
	    setup_chem(BaryonField[DensNum][n], temperature, EquilibrateChem,
		       BaryonField[DeNum][n], BaryonField[HINum][n], BaryonField[HIINum][n],
		       BaryonField[HeINum][n], BaryonField[HeIINum][n], BaryonField[HeIIINum][n],
		       BaryonField[HMNum][n], BaryonField[H2INum][n], BaryonField[H2IINum][n],
		       BaryonField[DINum][n], BaryonField[DIINum][n], BaryonField[HDINum][n]);
	  else if (MultiSpecies == 2) {
	    float temp;
	    setup_chem(BaryonField[DensNum][n], temperature, EquilibrateChem,
		       BaryonField[DeNum][n], BaryonField[HINum][n], BaryonField[HIINum][n],
		       BaryonField[HeINum][n], BaryonField[HeIINum][n], BaryonField[HeIIINum][n],
		       BaryonField[HMNum][n], BaryonField[H2INum][n], BaryonField[H2IINum][n],
		       temp, temp, temp);
	  }
	  else {
	    float temp;
	    setup_chem(BaryonField[DensNum][n], temperature, EquilibrateChem,
		       BaryonField[DeNum][n], BaryonField[HINum][n], BaryonField[HIINum][n],
		       BaryonField[HeINum][n], BaryonField[HeIINum][n], BaryonField[HeIIINum][n],
		       temp, temp, temp,
		       temp, temp, temp);
	  }
	} // if(MultiSpecies)
	
      } // end loop over grids

  return SUCCESS;

} // end Grid::GalaxySimulationInitializeGrid

/*
* Initialization routines
* Order:
*   NFW mass (NFWDarkMatterMassEnclosed)
*   cell mass (gauss_mass, rot_to_disk)
*   disk velocity (gas_vel OR ???)
*   chemistry (setup_chem, bilinear_interp)
*   CGM profile
*/


/* halo galaxy mass at a given radius, using user-defined global parameters for galaxy
   quantities and assuming that all halo mass is in an NFW halo.  This is not totally
   correct near the center of the halo, but since we're using it for the CGM initialization 
   and are dealing with radii that aren't particularly near the center of the halo, this 
   approximation is probably fine. 

   Input is the radius in CGS units; output is the enclosed mass at that radius in CGS units.
*/
double NFWDarkMatterMassEnclosed(double r){
  double M, C, R200, rho_0, Rs, M_within_r;
  double rho_crit = 1.8788e-29*0.49;
  
  // GSDarkMatterConcentration is the same as DiskGravityDarkMatterConcentration
  // if the latter is in use, same with GSGalaxyMass & DGDarkMatterMass

  M = GalaxySimulationGalaxyMass * SolarMass;  // halo total mass in CGS
  C = GalaxySimulationDMConcentration;  // concentration parameter for NFW halo
  
  R200 = POW(3.0/(4.0*pi)*M/(200.*rho_crit),1./3.);  // virial radius in CGS
  Rs = R200/C;  // scale radius of NFW halo in CGS
  rho_0 = 200.0*POW(C,3)/3.0/(log(1.0+C) - C/(1.0+C))*rho_crit;  // rho_0 for NFW halo in CGS
  // mass w/in radius R
  M_within_r = 4.0*pi*rho_0*POW(Rs,3.0)*(log((Rs+r)/Rs) - r/(Rs+r));
  return M_within_r;
  
}

//function that gets the mass enclosed at a given radius (input in CGS) returns in CGS 
double MassEnclosed_r(FLOAT *binned_mass, double rad){
	rad /= LengthUnits; //to code
    	double MassUnitsDouble = double(DensityUnits)*POW(double(LengthUnits), 3.0);
	double delta_r = 1.0 / 100; 
	int bin_index = rad / delta_r;
	double prev_menc = 0.0;
	double next_menc = 0.0; 
	double this_r_bin = delta_r * (bin_index + 0.5); 
	double prev_r_bin = delta_r * ((bin_index - 1) + 0.5); 
	double next_r_bin = delta_r * ((bin_index + 1) + 0.5); 
	double ans; 
	if(prev_r_bin < 0){
		ans = binned_mass[bin_index] / delta_r * (rad); 
		return ans*MassUnitsDouble; 	
	}
	if(next_r_bin >= 100){
		prev_menc = binned_mass[bin_index - 1];
	    	ans = (binned_mass[bin_index] - prev_menc) / delta_r * (rad - prev_r_bin) + prev_menc; 
		return ans*MassUnitsDouble; 	
	}
	prev_menc = binned_mass[bin_index - 1]; 
	next_menc = binned_mass[bin_index + 1]; 
	ans = 0.5 * (((binned_mass[bin_index] - prev_menc) / (delta_r) * (rad - prev_r_bin) + prev_menc) + ((next_menc - binned_mass[bin_index]) / delta_r * (next_r_bin - rad) + binned_mass[bin_index])); 
	return ans*MassUnitsDouble; //cgs 
}

// Computes the total mass in a given cell by integrating the density profile using 5-point Gaussian quadrature
float gauss_mass(FLOAT r, FLOAT z, FLOAT xpos, FLOAT ypos, FLOAT zpos, FLOAT inv [3][3], float DiskDensity, FLOAT ScaleHeightR, FLOAT ScaleHeightz, FLOAT cellwidth)
{
  
    FLOAT EvaluationPoints [5] = {-0.90617985,-0.53846931,0.0,0.53846931,0.90617985};
    FLOAT Weights [5] = {0.23692689,0.47862867,0.56888889,0.47862867,0.23692689};
    FLOAT xResult [5];
    FLOAT yResult [5];
    float Mass = 0;
    FLOAT xrot,yrot,zrot;
    int i,j,k;
    FLOAT rrot;
    FLOAT rcoord, zcoord;
    float  kludge2;

    for (i=0;i<5;i++) {

        xResult[i] = 0.0;
        for (j=0;j<5;j++) {
        
            yResult[j] = 0.0;
            for (k=0;k<5;k++) {

                rot_to_disk(xpos+EvaluationPoints[i]*cellwidth/2.0,ypos+EvaluationPoints[j]*cellwidth/2.0,zpos+EvaluationPoints[k]*cellwidth/2.0,xrot,yrot,zrot,inv);
                rrot = sqrt(POW(xrot,2)+POW(yrot,2));

                if( PointSourceGravity > 0 )
                    yResult[j] += cellwidth/2.0*Weights[k]*PEXP(-rrot/ScaleHeightR)/POW(cosh(zrot/(2.0*ScaleHeightz)),2);
                else if( DiskGravity > 0 ){
                    if( rrot/Mpc < SmoothRadius )
                        yResult[j] += cellwidth/2.0*Weights[k]/cosh(rrot/ScaleHeightR)/cosh(fabs(zrot)/ScaleHeightz);
                    else if( rrot/Mpc < TruncRadius )
                        yResult[j] += cellwidth/2.0*Weights[k]/cosh(rrot/ScaleHeightR)/cosh(fabs(zrot)/ScaleHeightz)
                                        *0.5*(1.0+cos(pi*(rrot-SmoothRadius*Mpc)/(SmoothLength*Mpc)));
                } // end disk gravity if
		else{
                    rcoord = sqrt((POW(xpos+EvaluationPoints[i]*cellwidth/2.0, 2.0) +
                                POW(ypos+EvaluationPoints[j]*cellwidth/2.0, 2.0) ) );
                    zcoord = fabs(zpos+EvaluationPoints[k]*cellwidth/2.0);
                    kludge2 = 
                        cellwidth/2.0 * Weights[k] * 
                        PEXP(-rcoord/ScaleHeightR) *
                        PEXP(-fabs(zcoord)/ScaleHeightz);
                    yResult[j] +=kludge2;

		}

            }
            xResult[i] += cellwidth/2.0*Weights[j]*yResult[j];
        }
        Mass += cellwidth/2.0*Weights[i]*xResult[i];
    }  
    Mass *= DiskDensity;
    return Mass;
}

//Finds coordinates in rotated coordinate system; used by gauss_mass
void rot_to_disk(FLOAT xpos, FLOAT ypos, FLOAT zpos, FLOAT &xrot, FLOAT &yrot, FLOAT &zrot, FLOAT inv [3][3])
{
  xrot = xpos*inv[0][0] + ypos*inv[0][1] + zpos*inv[0][2];
  yrot = xpos*inv[1][0] + ypos*inv[1][1] + zpos*inv[1][2];
  zrot = xpos*inv[2][0] + ypos*inv[2][1] + zpos*inv[2][2];
}



//
// Disk velocity with PointSourceGravity
//
float gasvel(FLOAT radius, float DiskDensity, FLOAT ExpansionFactor, float GalaxyMass, FLOAT ScaleHeightR, FLOAT ScaleHeightz, float DMConcentration, FLOAT Time)
{

 double OMEGA=OmegaLambdaNow+OmegaMatterNow;                 //Flat Universe

 double r = radius*LengthUnits/100;    // Radius [m]

 double M_200 = GalaxyMass*SolarMass/1000.0;      // Virial Mass [kg]

 double H = sqrt(HubbleConstantNow*100*HubbleConstantNow*100*(OmegaLambdaNow+OmegaMatterNow*POW(ExpansionFactor,-3)-(OMEGA-1.)*POW(ExpansionFactor,-2)));                                

 double r_200 = (1.63e-2*POW(GalaxyMass,1.0/3.0)*POW((OmegaLambdaNow+OmegaMatterNow*POW(ExpansionFactor, -3)-(OMEGA-1.0)*POW(ExpansionFactor,-2)),-1.0/3.0)*ExpansionFactor*POW(H,-2.0/3.0)*POW(100,2.0/3.0))*Mpc/1.0e5;
 //virial radius [m]: M_200/M_Solar = GalaxyMass

 double M_gas, M_DM, M_Tot, Acc, V_Circ;
 double f_C = log(1.0+DMConcentration)-DMConcentration/(1.0+DMConcentration);
 double r_s = r_200/DMConcentration;  //[m]

 // Mass of gas disk and DM at given radius

  M_gas=8.0*M_PI*ScaleHeightz*Mpc/100*ScaleHeightR*Mpc/100*ScaleHeightR*Mpc/100*DiskDensity*DensityUnits*1000*PEXP(-r/(ScaleHeightR*Mpc/100))*(PEXP(r/(ScaleHeightR*Mpc/100))-r/(ScaleHeightR*Mpc/100)-1.0);

     M_DM=(M_200/f_C)*(log(1.0+r/r_s)-(r/r_s)/(1.0+r/r_s));

     if (SelfGravity==1){
        M_Tot=M_DM+M_gas;
     }
     else{
        M_Tot=M_DM;
     }

  float DensityUnits=1, LengthUnits=1, VelocityUnits=1, TimeUnits=1,
    TemperatureUnits=1;
  double MassUnits=1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
         &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }

  double MassUnitsDouble=1.0;

  if(ComovingCoordinates)
    MassUnitsDouble = double(DensityUnits)*POW(double(LengthUnits), 3.0);

  // Set the point source gravity parameters.  This is the DM mass (in g)
  //   within rs.  The core radius to rs in cm.
  //
  // BWO 10 July 2009: Both of these values are now converted to code units, because 
  // otherwise the values go over 32-bit precision.  This is used in
  // Grid::ComputeAccelerationFieldExternal, and converted back to CGS where needed.
  //

  PointSourceGravityConstant = (M_200/f_C)*(log(1.0+1.0)-1.0/(1.0+1.0))*1000.0 / MassUnitsDouble;
  PointSourceGravityCoreRadius = r_s*100.0 / LengthUnits;

  // Force per unit mass on disk (i.e. acceleration) [ms-2]

  Acc=((GravConst/1000.0)*M_Tot)/(r*r);

 // Magnitude of Circular Velocity of disk 

  V_Circ = sqrt(r*Acc)*100;       //cms-1

  return (V_Circ/VelocityUnits);  //code units
}


//
// Disk velocity with DiskGravity. Arguments should be in cgs
//
double DiskGravityStellarAccel(double rcyl, double z){
    // return (cylindrical) radial component of acceleration ...? in cgs?
    // Potential from  Miyamoto & Nagai 75

    double accelcylR;

    accelcylR = GravConst * DiskGravityStellarDiskMass*SolarMass * rcyl
              / sqrt(POW(
                    POW(rcyl,2)
                  + POW(DiskGravityStellarDiskScaleHeightR*Mpc
                        + sqrt(POW(z,2)
                             + POW(DiskGravityStellarDiskScaleHeightz*Mpc,2)),
                        2),
                     3));

    return accelcylR;
}

double DiskGravityBulgeAccel(double rsph) { // cgs arguments
    double accelsph;

    accelsph = GravConst * DiskGravityStellarBulgeMass*SolarMass
             / POW(rsph + DiskGravityStellarBulgeR*Mpc,2);

    return accelsph;
}

double DiskGravityCircularVelocity(double rsph, double rcyl, double z) {
    double acc, velmag;
    
    acc = GravConst*NFWDarkMatterMassEnclosed(rsph)/POW(rsph,2)
        + DiskGravityStellarAccel(rcyl, z)
        + DiskGravityBulgeAccel(rsph);

    velmag = sqrt(acc*rcyl); // cgs
    return velmag;
}

/* Function for initializing chemistry */
void setup_chem(float density, float temperature, int equilibrate,
		float& DEdest, float& HIdest, float& HIIdest,
		float& HeIdest, float& HeIIdest, float& HeIIIdest,
		float& HMdest, float& H2Idest, float& H2IIdest,
		float& DIdest, float& DIIdest, float& HDIdest)
{
  if (equilibrate) {
    /*  What temperature and density bins does the cell fall between? 
     *  'density' is in code units; 'temperature' is K
     *  Table arrays contain fractions.
     *  densities are returned in code units
     */
        
    // Start by assuming values are larger than those in table;
    // set to dim_size-1 for highest available value
    bool interpolate = true;
    int dens_indx, temp_indx, iter;
    dens_indx = temp_indx = EquilibriumTable.dim_size-1;
    
    for (iter=0; iter < EquilibriumTable.dim_size; ++iter) {
      if (density < EquilibriumTable.density[iter]) {
	dens_indx = iter-1;
	break;
      }
    }
    
    for (iter=0; iter<EquilibriumTable.dim_size; ++iter) {
      if (temperature < EquilibriumTable.temperature[iter]) {
	temp_indx = iter-1;
	break;
      }
    }

    // Density or temperature lower than in table
    if (dens_indx == -1) {
      dens_indx = 0;
      interpolate = false;
    }
    if (temp_indx == -1){
      temp_indx = 0;
      interpolate = false;
    }

    // Density or temp higher than table; unchanged from inital value
    if (dens_indx == EquilibriumTable.dim_size-1 ||
	temp_indx == EquilibriumTable.dim_size-1)
      interpolate = false;

    if (interpolate) {
      HIdest = density*bilinear_interp(density, temperature, 
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.HI[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.HI[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.HI[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.HI[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);

      HIIdest =  density*bilinear_interp(density, temperature, 
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.HII[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.HII[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.HII[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.HII[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);

      HeIdest =  density*bilinear_interp(density, temperature, 
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.HeI[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.HeI[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.HeI[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.HeI[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);

      HeIIdest =  density*bilinear_interp(density, temperature, 
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.HeII[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.HeII[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.HeII[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.HeII[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);

      HeIIIdest =  density*bilinear_interp(density, temperature, 
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.HeIII[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.HeIII[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.HeIII[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.HeIII[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);
      
      DEdest =  density*bilinear_interp(density, temperature, 
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.de[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.de[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.de[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.de[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);

      if (MultiSpecies > 1) {
	HMdest =  density*bilinear_interp(density, temperature, 
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.HM[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.HM[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.HM[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.HM[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);

	H2Idest =  density*bilinear_interp(density, temperature,
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.H2I[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.H2I[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.H2I[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.H2I[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);

	H2IIdest =  density*bilinear_interp(density, temperature,
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.H2II[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.H2II[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.H2II[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.H2II[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);
      }
      if (MultiSpecies > 2) {
	DIdest =  density*bilinear_interp(density, temperature, 
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.DI[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.DI[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.DI[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.DI[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);

	DIIdest =  density*bilinear_interp(density, temperature, 
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.DII[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.DII[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.DII[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.DII[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);

	HDIdest =  density*bilinear_interp(density, temperature, 
  EquilibriumTable.density[dens_indx],
  EquilibriumTable.density[dens_indx+1],
  EquilibriumTable.temperature[temp_indx],
  EquilibriumTable.temperature[temp_indx+1],
  EquilibriumTable.HDI[EquilibriumTable.dim_size * temp_indx + dens_indx],
  EquilibriumTable.HDI[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx],
  EquilibriumTable.HDI[EquilibriumTable.dim_size * temp_indx + dens_indx+1],
  EquilibriumTable.HDI[EquilibriumTable.dim_size * (temp_indx+1) + dens_indx+1]);
      }
    } // end interpolate
    else { // don't interpolate; density and/or temp at edge of table
      HIdest =  EquilibriumTable.HI[EquilibriumTable.dim_size * temp_indx + dens_indx];

      HIIdest = EquilibriumTable.HII[EquilibriumTable.dim_size * temp_indx + dens_indx];

      HeIdest = EquilibriumTable.HeI[EquilibriumTable.dim_size * temp_indx + dens_indx];

      HeIIdest = EquilibriumTable.HeII[EquilibriumTable.dim_size * temp_indx + dens_indx];

      HeIIIdest = EquilibriumTable.HeIII[EquilibriumTable.dim_size * temp_indx + dens_indx];

      DEdest =  EquilibriumTable.de[EquilibriumTable.dim_size * temp_indx + dens_indx];

      if (MultiSpecies > 1) {
	HMdest =  EquilibriumTable.HM[EquilibriumTable.dim_size * temp_indx + dens_indx];

	H2Idest = EquilibriumTable.H2I[EquilibriumTable.dim_size * temp_indx + dens_indx];

	H2IIdest = EquilibriumTable.H2II[EquilibriumTable.dim_size * temp_indx + dens_indx];
      }
      if (MultiSpecies > 2) {
	DIdest =  EquilibriumTable.DI[EquilibriumTable.dim_size * temp_indx + dens_indx];

	DIIdest = EquilibriumTable.DII[EquilibriumTable.dim_size * temp_indx + dens_indx];

	HDIdest = EquilibriumTable.HDI[EquilibriumTable.dim_size * temp_indx + dens_indx];
      }
    } // end no interpolation
  } // end if equilibrate
  else {
    HIdest = TestProblemData.HI_Fraction * density *
      TestProblemData.HydrogenFractionByMass;

    HeIdest = TestProblemData.HeI_Fraction * density *
      (1.0-TestProblemData.HydrogenFractionByMass);
       
    HeIIdest = TestProblemData.HeII_Fraction * density *
      (1.0-TestProblemData.HydrogenFractionByMass);
       
    HeIIIdest = (1.0 - TestProblemData.HydrogenFractionByMass) *
      density - HeIdest - HeIIdest;

    if(MultiSpecies > 1){
      HMdest = TestProblemData.HM_Fraction *
	TestProblemData.HydrogenFractionByMass * density;
   
      H2Idest = 2 * TestProblemData.H2I_Fraction *
	TestProblemData.HydrogenFractionByMass * density;
   
      H2IIdest = 2 * TestProblemData.H2II_Fraction 
	* TestProblemData.HydrogenFractionByMass * density;
    }

    // HII density is calculated by subtracting off the various ionized fractions
    // from the total
    HIIdest = TestProblemData.HydrogenFractionByMass * density - HIdest;
    if (MultiSpecies > 1)
      HIIdest -= (HMdest + H2IIdest + H2Idest);

    // Electron "density" (remember, this is a factor of m_p/m_e scaled from the 'normal'
    // density for convenience) is calculated by summing up all of the ionized species.
    // The factors of 0.25 and 0.5 in front of HeII and HeIII are to fix the fact that we're
    // calculating mass density, not number density (because the BaryonField values are 4x as
    // heavy for helium for a single electron)
    DEdest = HIIdest + 0.25*HeIIdest + 0.5*HeIIIdest;
    
    if (MultiSpecies > 1)
      DEdest += 0.5*H2IIdest - HMdest;
       
    DEdest = max(DEdest, tiny_number);
       
    // Set deuterium species (assumed to be a negligible fraction of the total, so not
    // counted in the conservation)
    if(MultiSpecies > 2){
      DIdest = HIdest * TestProblemData.DeuteriumToHydrogenRatio;
      DIIdest = HIIdest *	TestProblemData.DeuteriumToHydrogenRatio;
      HDIdest = H2Idest *	0.75 * TestProblemData.DeuteriumToHydrogenRatio;
    }

  } // end not equilibrate
}
double bilinear_interp(double x, double y, 
                       double x1, double x2, double y1, double y2,
                       double f_x1y1, double f_x1y2, 
                       double f_x2y1, double f_x2y2) {
    double interp;

    interp = f_x1y1*(x2-x)*(y2-y) + f_x2y1*(x-x1)*(y2-y) 
           + f_x1y2*(x2-x)*(y-y1) + f_x2y2*(x-x1)*(y-y1);

    interp *= 1/( (x2-x1)*(y2-y1) );

    return interp;

}

/* -------------------- BEGINNING OF Routines used for initializing the circumgalactic medium -------------------- */
/* 
   Computes halo gas density values assuming a variety of user-specifiable models
   for the CGM, toggled by the variable GalaxySimulationGasHalo.  Depending on the
   specific model chosen, different global parameters are needed (as set near the beginning
   of Grid::GalaxySimluationInitializeGrid).  Halo types are:

   GalaxySimulationGasHalo = 0  -- "zero CGM" - sets to a very low density/temperature
   GalaxySimulationGasHalo = 1  -- assuming hydrostatic equilibrium of CGM given an NFW dark matter halo 
                                   and a temperature as a function of radius set by the virial theorem.
   GalaxySimulationGasHalo = 2  -- assumes density, temperature set according to T = Tvir and entropy
                                   as a power-law function of radius.
   GalaxySimulationGasHalo = 3  -- as #2, but the entropy distribution has a floor value, so S = S_f + S_0 (r/r_0)^alpha
   GalaxySimulationGasHalo = 4  -- assumes a hydrostatic equilibrium of CGM given an NFW dark matter halo
                                   and an entropy that is a power-law function of radius.
   GalaxySimulationGasHalo = 5  -- as #4, but the entropy distribution has a floor value, so S = S_f + S_0 (r/r_0)^alpha
   GalaxySimulationGasHalo = 6  -- as #4, but the entropy distribution follows that for a precipitation-regulated NFW halo
                                   in Voit 2019 (ApJ)
   GalaxySimulationGasHalo = 7  -- previously reserved
   GalaxySimulationGasHalo = 8  -- a density and temperature profile fit to the entropy profiles in Voit 2019 (ApJ)     

   Inputs:  R - spherical radius, code units

   Returns:  density, grams/cm^3

   Note: using global variables w/following units:

   GalaxySimulationGasHalo: integer, >= 0
   GalaxySimulationGasHaloScaleRadius, units of Mpc
   GalaxySimulationGasHaloDensity, units of grams/cm^3
   GalaxySimulationGasHaloTemperature, units of Kelvin
   GalaxySimulationGasHaloAlpha, power-law index; unitless
   GalaxySimulationGasHaloCoreEntropy, units of keV cm^2
   GalaxySimulationGasHaloMetallicity, units of Zsun
*/
float HaloGasDensity(FLOAT R, struct CGMdata& CGM_data){

  if(GalaxySimulationGasHalo < 1){
    /* "zero CGM" - sets a very low density */
   
    return densicm;

  } else if(GalaxySimulationGasHalo == 1){
    /* gets density assuming hydrostatic equilibrium using a temperature 
       as a function of radius given by virial theorem */
    
    double T0,haloDensity;
    T0 = HaloGasTemperature(GalaxySimulationGasHaloScaleRadius*Mpc/LengthUnits,
			    CGM_data);
    haloDensity = GalaxySimulationGasHaloDensity*(T0/HaloGasTemperature(R, CGM_data));
    haloDensity /= POW((R*LengthUnits/GalaxySimulationGasHaloScaleRadius/Mpc),3);
    return min(haloDensity,GalaxySimulationGasHaloDensity);
    
  } else if(GalaxySimulationGasHalo == 2){
    /* assumes entropy is a power-law function of radius and T = Tvir, so
       n(r) = n_0 * (r/r_0)**(-alpha/(gamma-1))
       where n_0 is (user-supplied) number density at (user-supplied) radius r_0,
       alpha is (user-supplied) power-law exponent,
       gamma is adiabatic index.
    */
    double scale_radius_cgs, this_radius_cgs, power_law_exponent;
    
    scale_radius_cgs = GalaxySimulationGasHaloScaleRadius*Mpc;
    this_radius_cgs = R*LengthUnits;
    power_law_exponent = -1.0*GalaxySimulationGasHaloAlpha/(Gamma-1.0);
    
    return GalaxySimulationGasHaloDensity*POW(this_radius_cgs/scale_radius_cgs, power_law_exponent);
    
  } else if(GalaxySimulationGasHalo == 3){
    /* assumes entropy is a  power-law function of radius and T = Tvir that has a core (i.e., minimum 
       entropy value), so:

       n(r) = (Tvir / (Score + S_0*(r/r_0)^alpha))^(1/(gamma-1))

       where n_0 is (user-supplied) number density at (user-supplied) radius r_0,
       Tvir is user-supplied temperature,
       Score is a user-supplied entropy,
       alpha is (user-supplied) power-law exponent,
       gamma is adiabatic index.
    */

    double scale_radius_cgs, this_radius_cgs, power_law_exponent, S_0, T_kev, n_0, this_number_density;

    // get radii in common set of units
    scale_radius_cgs = GalaxySimulationGasHaloScaleRadius*Mpc;
    this_radius_cgs = R*LengthUnits;

    // now get number density using expression above

    T_kev = GalaxySimulationGasHaloTemperature*kboltzKeV;  // halo temperature in keV
    n_0 = GalaxySimulationGasHaloDensity / (mu*mh);  // convert n_0 to electron number density 
    S_0 = T_kev / POW(n_0,Gamma-1.0);   // S_0 in units of kev cm^2

    // get number density at this radius giving the requested info
    this_number_density = T_kev/(GalaxySimulationGasHaloCoreEntropy + S_0*POW(this_radius_cgs/scale_radius_cgs, power_law_exponent));
    this_number_density = POW(this_number_density, 1.0/(Gamma-1.0));

    return this_number_density*mu*mh;  // return physical density
    
  } else if(GalaxySimulationGasHalo >= 4 && GalaxySimulationGasHalo <= 7){
    /* assumes entropy is a power-law function of radius OR a cored power-law function
       of radius and gas is in hydrostatic equilibrium w/the NFW halo.  */

    double this_radius_cgs, Rstart;
    int index;

    this_radius_cgs = R*LengthUnits;  // radius in CGS
    index = int((this_radius_cgs-CGM_data.R_inner)/CGM_data.dr + 1.0e-3);  // index in array of CGM values
    if(index<0) index=0;  // check our indices
    if(index>=CGM_data.nbins) index=CGM_data.nbins-1;
    return CGM_data.n_rad[index]*mu*mh;  // return physical density

  } else if(GalaxySimulationGasHalo == 8){
    /* Eqn 24 in the Appendix of Voit 2019; a fit to the theoretical density profile of a precipitation-regulated NFW halo.
       Equation gives the electron number density, but pull the same trick as methods 2 & 3 and assume n_e = n */
    double this_radius_kpc, this_number_density;
    this_radius_kpc = R*LengthUnits/CM_PER_KPC;
    
    this_number_density = POW( POW(this_radius_kpc,GalaxySimulationGasHaloZeta) / GalaxySimulationGasHaloDensity, 2);
    this_number_density += POW( POW(this_radius_kpc/100,GalaxySimulationGasHaloZeta2) / GalaxySimulationGasHaloDensity2, 2);
    this_number_density = POW(this_number_density, -0.5);

    return this_number_density*mu*mh;  // return physical density
    
  } else {
    ENZO_FAIL("Grid::GalaxySimulationInitializeGrid - invalid choice of GalaxySimulationGasHalo in HaloGasDensity().");
  }
  
} // end HaloGasDensity


/* 
   Computes halo gas temperature values assuming a variety of user-specifiable models
   for the CGM, toggled by the variable GalaxySimulationGasHalo.  The properties of each of the
   models are described immediately above this in the comments for the function HaloGasDensity(). 

   Inputs:  R - spherical radius, code units

   Returns:  Temperature, Kelvin
*/

float HaloGasTemperature(FLOAT R, struct CGMdata& CGM_data){

  if(GalaxySimulationGasHalo < 1){
    /* "zero CGM" - sets a very low temperature */
   
    return Ticm;

  } else if(GalaxySimulationGasHalo == 1){
    /* gets temperature as a function of radius given by virial theorem */
    printf("\n NFWDarkMatterMassEnclosed(R) = %e \n", NFWDarkMatterMassEnclosed(R));
    return GravConst*NFWDarkMatterMassEnclosed(R)*mu*mh/(3.0*kboltz*R*LengthUnits);
    
  } else if(GalaxySimulationGasHalo == 2){
    /* assumes entropy is a power-law function of radius and T = Tvir */

    return GalaxySimulationGasHaloTemperature;
    
  } else if(GalaxySimulationGasHalo == 3){

    /* assumes entropy is a cored power-law function of radius and T = Tvir */

    return GalaxySimulationGasHaloTemperature;

  } else if(GalaxySimulationGasHalo >= 4 && GalaxySimulationGasHalo <= 7){
    /* assumes entropy is a power-law function of radius and gas is in hydrostatic equilibrium */

    double this_radius_cgs;
    int index;
    this_radius_cgs = R*LengthUnits;  // radius in CGS
    index = int((this_radius_cgs-CGM_data.R_inner)/CGM_data.dr+1.0e-3);  // index in array of CGM values
    if(index<0) index=0;  // check our indices
    if(index>=CGM_data.nbins) index=CGM_data.nbins-1;
    return CGM_data.T_rad[index];  // return temperature in Kelvin

  } else if(GalaxySimulationGasHalo == 8){
    /* Theoretical temperature profile of a precipitation-regulated NFW halo, using fits to n(r) and S(r) */
    double this_radius_kpc, this_number_density, this_entropy;
    this_radius_kpc = R*LengthUnits/CM_PER_KPC;
    
    this_number_density = POW( POW(this_radius_kpc,GalaxySimulationGasHaloZeta) / GalaxySimulationGasHaloDensity, 2);
    this_number_density += POW( POW(this_radius_kpc/100,GalaxySimulationGasHaloZeta2) / GalaxySimulationGasHaloDensity2, 2);
    this_number_density = POW(this_number_density, -0.5);

    this_entropy = GalaxySimulationGasHaloCoreEntropy * POW(this_radius_kpc, GalaxySimulationGasHaloAlpha);

    return this_entropy * POW(this_number_density, Gamma-1.0) / kboltzKeV; // units of K
    
  } else {
    ENZO_FAIL("Grid::GalaxySimulationInitializeGrid - invalid choice of GalaxySimulationGasHalo in HaloGasTemperature().");
  }
  
}

/* Initializes arrays of number density, temperature, and radius for 
   choices of circumgalactic medium that require numerical integration based
   on user-defined parameters.  These quantities are stored in a global struct
   with arrays of size nbins for convenience (global within this file, at least).
   Rstop is the outer boundary of the integrtation in CGS; if negative, |Rstop|*R200 is used. 
   nbins defaults to 8192. */
 void halo_init(struct CGMdata& CGM_data, grid* Grid, FLOAT *binned_mass, float Rstop, int GasHalo_override){

  int halo_type=GalaxySimulationGasHalo;
  if (GasHalo_override) // not 0
    halo_type = GasHalo_override;
  
  if(halo_type < 4 || halo_type > 7) return;

  double k1, k2, k3, k4;
  double M, R200, rho_crit = 1.8788e-29*0.49;
  double Rstart;

  int index;
  
  double MassUnitsDouble = double(DensityUnits)*POW(double(LengthUnits), 3.0);
  M = binned_mass[99] * MassUnitsDouble;  // halo total mass in CGS
   
  R200 = pow(3.0/(4.0*3.14159)*M/(200.*rho_crit),1./3.);  // virial radius in CGS
  if (Rstop < 0)
    Rstop = fabs(Rstop)*R200;
  CGM_data.R_outer = Rstop;// integrate out to the virial radius of halo

  Rstart = 0.00*LengthUnits; // you could force a different start if you liked

  // stepsize for RK4 integration and radial bins
  CGM_data.R_inner = Rstart;
  CGM_data.dr = (CGM_data.R_outer - CGM_data.R_inner)/ double(CGM_data.nbins); 
  
  if (halo_type < 6){
    double T0, n0, r0, dr;
    double this_n, this_radius, temperature;
    
    // set some quantities based on user inputs; this defines our integration
    T0 = GalaxySimulationGasHaloTemperature;
    n0 = GalaxySimulationGasHaloDensity / (mu*mh);
    r0 = GalaxySimulationGasHaloScaleRadius*Mpc;
  
    // used for our numerical integration
    dr = CGM_data.dr;
    this_n = n0;
    this_radius = r0;

    // set the bin that we start at (otherwise it doesn't get set!)
    index = int((this_radius - CGM_data.R_inner)/dr + 1.0e-3);
    CGM_data.n_rad[index] = this_n;
    CGM_data.T_rad[index] = T0;
    CGM_data.rad[index] = this_radius;

    /* starting at the point where the user has defined the radius, density, and 
       temperature, use RK4 to integrate the number density outward to R_outer using the expression 
       for dn_dr in another function.  Calculate the temperature using the entropy at this radius. */
    while(this_radius <= CGM_data.R_outer){
    
      // calculate RK4 coefficients.
      k1 = halo_dn_dr(this_radius,          this_n);
      k2 = halo_dn_dr(this_radius + 0.5*dr, this_n + 0.5*dr*k1);
      k3 = halo_dn_dr(this_radius + 0.5*dr, this_n + 0.5*dr*k2);
      k4 = halo_dn_dr(this_radius + dr,     this_n + dr*k3);

      // update density and radius
      this_n += (1.0/6.0) * dr * (k1 + 2.0*k2 + 2.0*k3 + k4);
      this_radius += dr;  // new radius

      // calculate temperature at this radius using known entropy 
      temperature = halo_S_of_r(this_radius) * POW(this_n,Gamma-1.0);
      
      // store everything in the struct
      index = int((this_radius - CGM_data.R_inner)/dr + 1.0e-3);    
      CGM_data.n_rad[index] = this_n;
      CGM_data.T_rad[index] = temperature;
      CGM_data.rad[index] = this_radius;
    }
        
    /* now we do the same thing as above, but integrating inward to zero radius. */
    this_n = n0;
    this_radius = r0;
    dr *= -1.0;
    
    while(this_radius > CGM_data.R_inner){
      
      // calculate RK4 coefficients.
      k1 = halo_dn_dr(this_radius,          this_n);
      k2 = halo_dn_dr(this_radius + 0.5*dr, this_n + 0.5*dr*k1);
      k3 = halo_dn_dr(this_radius + 0.5*dr, this_n + 0.5*dr*k2);
      k4 = halo_dn_dr(this_radius + dr,     this_n + dr*k3);

      // update density and radius
      this_n += (1.0/6.0) * dr * (k1 + 2.0*k2 + 2.0*k3 + k4);
      this_radius += dr;  // new radius
      
      // calculate temperature at this radius using known entropy 
      temperature = halo_S_of_r(this_radius) * POW(this_n,Gamma-1.0);

      // store everything in the struct
      index = int((this_radius - CGM_data.R_inner)/(-1.0*dr) + 1.0e-3);

      if(index >= 0){
	CGM_data.n_rad[index] = this_n;
	CGM_data.T_rad[index] = temperature;
	CGM_data.rad[index] = this_radius;
      }
    }

  } else if (halo_type == 6) { 

    /* Integrate pressure assuming HSE & S(r) from Voit 2019, then convert to n & T,
       instead of integrating n(r) directly as with methods 4 & 5. This makes the boundary
       condition easier to handle.*/
    double dr, rmax, vcirc2_max;
    double this_press, this_ent, this_radius;//, this_n;
    double mu_ratio = 1.17/mu; // mu_e/mu
    double T_floor = 4e4; // IGM, was 4e4

    // boundary condition & quantities for integration
    dr = -1.0*CGM_data.dr;
    this_radius = R200;
    this_ent = halo_S_of_r(this_radius, Grid, binned_mass); // in erg*cm^2
    rmax = 2.163*R200/GalaxySimulationDMConcentration;
    vcirc2_max = GravConst * MassEnclosed_r(binned_mass,rmax)/rmax;
    this_press = mu_ratio*POW(0.25*mu*mh*vcirc2_max/POW(this_ent, 1./Gamma),
		     Gamma/(Gamma-1.));
    float UniformDensity = 1e-30; //cgs
    float UniformTemperature = 1000; //K 
    this_press = UniformDensity*kboltz*UniformTemperature; //backround pressure 
    // set the bin that we start at (otherwise it doesn't get set!)
    index = int((this_radius - CGM_data.R_inner)/(-1.0*dr)  + 1.0e-3);
    CGM_data.n_rad[index] = 2 * POW(this_press/(mu_ratio*this_ent), 1./Gamma); // n_e ~ n_i
    CGM_data.T_rad[index] = POW( POW(this_press/mu_ratio, Gamma-1.) * this_ent, 1./Gamma) / kboltz;
    CGM_data.rad[index] = this_radius;

    int start_index = index; 
    // integrate inward from R200    
    while(this_radius > CGM_data.R_inner){
      double p1 = log10(mu_ratio*(CGM_data.n_rad[index] / 2)*kboltz*CGM_data.T_rad[index]);   
      // calculate RK4 coefficients.
      if(this_radius + dr < 0) //leave if new radius is negative 
          break; 
      k1 = halo_dP_dr(this_radius,          this_press,             Grid, binned_mass);
      k2 = halo_dP_dr(this_radius + 0.5*dr, this_press + 0.5*dr*k1, Grid, binned_mass);
      k3 = halo_dP_dr(this_radius + 0.5*dr, this_press + 0.5*dr*k2, Grid, binned_mass);
      k4 = halo_dP_dr(this_radius + dr,     this_press + dr*k3,     Grid, binned_mass);
      // update radius, pressure, entropy
      this_radius += dr;  // new radius
      this_press += (1.0/6.0) * dr * (k1 + 2.0*k2 + 2.0*k3 + k4); // P @ new radius
      this_ent = halo_S_of_r(this_radius, Grid, binned_mass); // entropy @ new radius
      // store density and temperature in the struct
      index = int((this_radius - CGM_data.R_inner)/(-1.0*dr) + 1.0e-3);
      CGM_data.n_rad[index] = 2 * POW(this_press/(mu_ratio*this_ent), 1./Gamma);
      CGM_data.T_rad[index] = POW(POW(this_press/mu_ratio, Gamma-1.) * this_ent, 1./Gamma) / kboltz;
      CGM_data.rad[index] = this_radius;
    }

    std::cout << "inward integration completed" << std::endl;  
    // Reset to boundary state
    dr = CGM_data.dr;
    this_radius = R200;
    this_ent = halo_S_of_r(this_radius, Grid, binned_mass); // in erg*cm^2
    rmax = 2.163*R200/GalaxySimulationDMConcentration;
    vcirc2_max = GravConst * MassEnclosed_r(binned_mass, rmax)/rmax;
    this_press = mu_ratio*POW(0.25*mu*mh*vcirc2_max/POW(this_ent, 1./Gamma),
		     Gamma/(Gamma-1.));
    this_press = UniformDensity*kboltz*UniformTemperature; //backround pressure 

    // Construct sigmoid to transition temperature to a constant
    double this_temp, this_dens;
    double deriv, r0, y0, y_offset, k;

    index = int((this_radius - CGM_data.R_inner)/(1.0*dr) + 1.0e-3);
    this_dens = 2 * POW(this_press/(mu_ratio*this_ent), 1./Gamma);
    this_temp = POW( POW(this_press/mu_ratio, Gamma-1.) * this_ent, 1./Gamma) / kboltz;

    deriv = (log10(this_temp) - log10(CGM_data.T_rad[index-1]))
          / (log10(this_radius) - log10(this_radius-dr));

    r0 = log10(this_radius);
    y0 = 2.0 * log10( T_floor / this_temp );
    //assert (y0 < 0.0);
    y_offset = log10(this_temp) - y0/2.0;
    k = fabs(4.0/y0 * deriv);

    // Set constant dlog(P)/dlog(r)
    double prev_press, dlP_dlr, this_dPdr, press_vir;
    prev_press = mu_ratio * CGM_data.n_rad[index-2]/2.0 * kboltz*CGM_data.T_rad[index-2];
    press_vir = this_press; 
    
    dlP_dlr = (log10(this_press) - log10(prev_press))
            / (log10(this_radius) - log10(this_radius-dr));
    //assert (dlP_dlr < 0.0);
    
    while(this_radius <= CGM_data.R_outer){
      this_dPdr = this_press/this_radius * dlP_dlr;
      //k1 = halo_dP_dr(this_radius,          this_press,             Grid, binned_mass);
      //k2 = halo_dP_dr(this_radius + 0.5*dr, this_press + 0.5*dr*k1, Grid, binned_mass);
      //k3 = halo_dP_dr(this_radius + 0.5*dr, this_press + 0.5*dr*k2, Grid, binned_mass);
      //k4 = halo_dP_dr(this_radius + dr,     this_press + dr*k3,     Grid, binned_mass);
      // update radius, pressure, entropy
      this_radius += dr;  // new radius
      //this_press += (1.0/6.0) * dr * (k1 + 2.0*k2 + 2.0*k3 + k4); // P @ new radius
      this_ent = halo_S_of_r(this_radius, Grid, binned_mass); // entropy @ new radius
      //std::cout << "outward integration pressure " << this_press << std::endl;
      // update density and radius
      this_dens = -2.0 * this_dPdr/(1.22*mh*halo_mod_g_of_r(this_radius, binned_mass)); // n_e = n_i
      this_temp = POW(10, sigmoid(log10(this_radius), r0, k, y0, y_offset));
      this_press = POW(10, dlP_dlr*log10(this_radius/R200) + log10(press_vir));
      // store everything in the struct
      index = int((this_radius - CGM_data.R_inner)/dr + 1.0e-3);    
      if (index < CGM_data.nbins) {
	CGM_data.n_rad[index] = 2 * POW(this_press/(mu_ratio*this_ent), 1./Gamma);
	CGM_data.T_rad[index] = POW(POW(this_press/mu_ratio, Gamma-1.) * this_ent, 1./Gamma) / kboltz;
	CGM_data.T_rad[index] = this_temp; 
	CGM_data.rad[index] = this_radius;
      }
      else
	break;
    }
  }
    
  if (CGM_data.R_inner == 0) {
    // this integration acts a little squirrelly around r=0 because the mass values are garbage.  Cheap fix.
    CGM_data.rad[0]=CGM_data.rad[1];
    CGM_data.n_rad[0]=CGM_data.n_rad[1];
    CGM_data.T_rad[0]=CGM_data.T_rad[1];
  }
  
  return;
}

/* Halo entropy as a function of radius for the user-specified CGM types that require numerical
   integration. 

   Input is radius in CGS units.  output is entropy in CGS units (Kelvin cm^2) 
*/
double halo_S_of_r(double r){

  double Tvir, n0, r0, Smin, S0;

  // calculate a bunch of things based on user inputs
  Tvir = GalaxySimulationGasHaloTemperature;  // in Kelvin
  n0 = GalaxySimulationGasHaloDensity / (mu*mh);  // convert from density to electron number density (cm^-3)
  r0 = GalaxySimulationGasHaloScaleRadius*Mpc;  // scale radius in CGS
  Smin = GalaxySimulationGasHaloCoreEntropy/8.621738e-8;  // given in keV cm^2, converted to Kelvin cm^2
  S0 = Tvir / POW(n0,Gamma-1);  // entropy at scale radius, in units of Kelvin cm^2

  if(GalaxySimulationGasHalo == 4){

    return S0*POW(r/r0,GalaxySimulationGasHaloAlpha);  // has units of Kelvin cm^2

  } else if (GalaxySimulationGasHalo == 5){

    return Smin + S0*POW(r/r0,GalaxySimulationGasHaloAlpha);  // has units of Kelvin cm^2

  } else {
    ENZO_FAIL("halo_S_of_r: GalaxySimulationGasHalo set incorrectly.");
  }

}

/* More complex entropy profile from Voit 2019 that requires calculation of the cooling function.
   This one returns entropy in erg cm^2 instead of K cm^2 */
double halo_S_of_r(double r, grid* Grid, FLOAT *binned_mass){
  if (GalaxySimulationGasHalo == 6){

    double M, C, r_vir, r_max, rho_crit = 1.8788e-29*0.49;
    double vcirc2, vcirc2_max;
    double Tgrav, Tgrav_therm;
    double MassUnitsDouble = double(DensityUnits)*POW(double(LengthUnits), 3.0);
    M = binned_mass[99]*MassUnitsDouble;  // total mass in CGS
    C = GalaxySimulationDMConcentration;  // concentration parameter for NFW halo
    r_vir = POW(3.0/(4.0*3.14159)*M/(200.*rho_crit),1./3.);  // virial radius in CGS
    r_max = 2.163 * r_vir/C;
    
    vcirc2 = GravConst * MassEnclosed_r(binned_mass, r) / r;
    vcirc2_max = GravConst * MassEnclosed_r(binned_mass, r_max) / r_max;
    
    Tgrav = mu*mh * vcirc2 / kboltz; // 2x gravitational "temperature"
    Tgrav_therm = Tgrav / TemperatureUnits / ((Gamma-1.0)*mu); // code
  
    /* Calculate the cooling function Lambda using Grackle */
    double Lambda;
    double dens = mh/DensityUnits; // code
    double vx=0, vy=0, vz=0;
    double hi, hii, hei, heii, heiii, de, hm, h2i, h2ii, di, dii, hdi, metal; // species
    int dim=1;

    // setup_chem has densities in code, temperature in K
    setup_chem(dens, Tgrav, EquilibrateChem, de, hi, hii, hei, heii, heiii, hm, h2i, h2ii, di, dii, hdi);
    metal = GalaxySimulationGasHaloMetallicity * CoolData.SolarMetalFractionByMass * dens;

    // temporarily disable UV background; makes S(r) trend downward at large r instead of upward
    // because of low Tgrav
    int saved_UVB = grackle_data->UVbackground;
    //grackle_data->UVbackground = 0;
    Grid->GrackleCustomCoolRate(1, &dim, &Lambda,
				&dens, &Tgrav_therm,
				&vx, &vy, &vz,
				&hi, &hii,
				&hei, &heii, &heiii,
				&de,
				&hm, &h2i, &h2ii,
				&di, &dii, &hdi,
				&metal);
    grackle_data->UVbackground = saved_UVB;

    // to cgs
    Lambda = fabs(Lambda) * POW(mh,2) * POW(LengthUnits,2) / ( POW(TimeUnits,3) * DensityUnits);
    
    double n_e = DensityUnits * de * 0.000544617 / 9.109e-28; //undo the grackle normalization from setup_chem 
    //n_e = Density*de; //Is this factor needed?

    double n_hi = DensityUnits * hi / mh; 
    double n_hii = DensityUnits * hii / mh; 
    double n_hm = DensityUnits * hm / mh; 
 
    double m_he = 6.64e-24; 

    double n_hei = DensityUnits * hei / m_he; 
    double n_heii = DensityUnits * heii / m_he; 
    double n_heiii = DensityUnits * heiii / m_he; 

    double m_h2 = 2*mh; 

    double n_h2i = DensityUnits * h2i / m_h2; 
    double n_h2ii = DensityUnits * h2ii / m_h2; 

    double m_d = 3.345e-24; 

    double n_di = DensityUnits * di / m_d; 
    double n_dii = DensityUnits * dii / m_d; 

    double m_hd = 5.018e-24; 

    double n_hd = DensityUnits * hdi / m_hd; 

    //double n_metal = DensityUnits * metal / (3*mh); 

    double n_i = n_hii + n_heii + n_heiii + n_h2ii + n_dii + n_hm; 

    double n = n_hi + n_hii + n_hm + n_hei + n_heii + n_heiii + n_h2i + n_h2ii + n_di + n_dii + n_hd + n_e; 
    /* Calculate entropy S(r) in erg cm^2 */
    double S_precip = POW(2*mu*mh, 1./3.) * POW(r*Lambda*GalaxySimulationGasHaloRatio/3.0, 2./3.);
    //double S_precip = POW(2*mu*mh, 1./3.) * POW(20 * r * Lambda * n_i / (n * 3), 2./3.); 
    double S_nfw = 39. * vcirc2_max/1e10/4e4 * POW(r/r_vir, 1.1) / KEV_PER_ERG; // See Voit 19 Eqn 10 for assumptions

    // TODO blend with an entropy cap
    return (S_nfw + S_precip);
    
  } else {
    ENZO_FAIL("halo_S_of_r: GalaxySimulationGasHalo set incorrectly.");
  }  
    
}

/* dEntropy/dr as a function of radius for the user-specified CGM types that require numerical
   integration. 

   Input is radius in CGS units; output is entropy gradient in CGS units (Kelvin*cm) */
double halo_dSdr(double r, double n){

  double Tvir, alpha, n0, r0, Smin, S0;

  // calculate a bunch of things based on user inputs
  Tvir = GalaxySimulationGasHaloTemperature;  // in Kelvin
  n0 = GalaxySimulationGasHaloDensity / (mu*mh);  // convert from density to electron number density (cm^-3)
  r0 = GalaxySimulationGasHaloScaleRadius*Mpc;  // scale radius in CGS
  Smin = GalaxySimulationGasHaloCoreEntropy/8.621738e-8;  // given in keV cm^2, converted to Kelvin cm^2
  S0 = Tvir / POW(n0,Gamma-1);  // entropy at scale radius, in units of Kelvin cm^2

  if(GalaxySimulationGasHalo == 4 || GalaxySimulationGasHalo == 5){

    // has units of Kelvin*cm, same deriv for both halo types (since constant drops out)
    return S0*GalaxySimulationGasHaloAlpha*
      POW(r/r0,GalaxySimulationGasHaloAlpha-1.0)/r0;

  } else {
    ENZO_FAIL("halo_dSdr: GalaxySimulationGasHalo set incorrectly.");
  }
}

/* a flexible, contortable sigmoid function. Gets used for log(K) */
double sigmoid(double x, double x0, double k, double y0, double y_off) {
  return y0 / (1 + exp(-k*(x-x0))) + y_off;
}

/* dn/dr as a function of radius and halo electron number density.  This quantity is calculated
   by assuming that gravity and pressure are in hydrostatic equilibrium in a halo with a specified 
   entropy profile S(r).

   Input is radius in CGM units and electron number density in units 
   of particles per cm^-3.  Output is dn/dr in CGS units, so particles per cm^4. 
*/
double halo_dn_dr(double r, double n){
  
  return -1.0*( n*1.22*mh*halo_g_of_r(r) + kboltz*POW(n,Gamma)*halo_dSdr(r,n) ) /
    ( Gamma * kboltz * halo_S_of_r(r) * POW(n, Gamma-1));
}

// double halo_dn_dr(double r, double n, double T) {
//   return -1.0 *  * 1.22*mh*n / (kboltz*T);
// }

double halo_dP_dr(double r, double P, grid* Grid, FLOAT *binned_mass) {
    double ret =  -1.0 * halo_mod_g_of_r(r, binned_mass) * 1.22 * mh * POW( P/(1.1/mu) / halo_S_of_r(r,Grid, binned_mass),
						    1./Gamma );
    if(halo_mod_g_of_r(r, binned_mass) < 0)
	    ENZO_FAIL("negative g"); 
    if(ret > 0)
        ENZO_FAIL("positive dp/dr"); 
    if(isnan(ret) && !isnan(P)){
	std::cout << "NAN IN dP_dr" << std::endl; 
	std::cout << "P " << P << std::endl;
        std::cout << "POW " << POW(P/(1.1/mu) / halo_S_of_r(r,Grid,binned_mass), 1./Gamma) << std::endl;	
        std::cout << "halo s of r " << halo_S_of_r(r,Grid, binned_mass) << " r " << r << std::endl;
	std::cout << "g" << halo_mod_g_of_r(r, binned_mass) << std::endl; 
	std::cout.flush();
        ENZO_FAIL("nan in dp/dr");
    }
    return ret;
}

/* halo gravitational acceleration as a function of radius.

   Input is the radius in CGS units and returns the MAGNITUDE of the 
   acceleration in CGS units.  */
double halo_g_of_r(double r){
  return GravConst*NFWDarkMatterMassEnclosed(r)/(r*r); 
}

double halo_mod_g_of_r(double r, FLOAT *binned_mass){
  return GravConst*MassEnclosed_r(binned_mass, r)/(r*r);
}

/* A modified NFW halo, where small radii (r <= 2.163*Rs) have enclosed mass
   that satisfies v_circ(r) = v_circ,max = v_circ(2.163*Rs). Beyond this
   boundary radius, the enclosed mass is that of a standard NFW halo.
*/
double halo_mod_DMmass_at_r(double r){

  double M, C, R200, Rs, Rmax;
  double rho_crit = 1.8788e-29*0.49;
  
  M = GalaxySimulationGalaxyMass * SolarMass;  // halo total mass in CGS
  C = GalaxySimulationDMConcentration;  // concentration parameter for NFW halo
  
  R200 = POW(3.0/(4.0*3.14159)*M/(200.*rho_crit),1./3.);  // virial radius in CGS
  Rs = R200/C;  // scale radius of NFW halo in CGS
  Rmax = 2.163*Rs;
  
  if (r <= Rmax) {
    return r / Rmax * NFWDarkMatterMassEnclosed(Rmax);
  }
  else {
    return NFWDarkMatterMassEnclosed(r);
  }
}

double halo_mod_DMmass_at_r(float* binned_mass, double r){

  double M, C, R200, Rs, Rmax;
  double rho_crit = 1.8788e-29*0.49;
  
  M = GalaxySimulationGalaxyMass * SolarMass;  // halo total mass in CGS
  C = GalaxySimulationDMConcentration;  // concentration parameter for NFW halo
  
  R200 = POW(3.0/(4.0*3.14159)*M/(200.*rho_crit),1./3.);  // virial radius in CGS
  Rs = R200/C;  // scale radius of NFW halo in CGS
  Rmax = 2.163*Rs;
  
  if (r <= Rmax) {
    return r / Rmax * MassEnclosed_r(binned_mass, Rmax);
  }
  else {
    return MassEnclosed_r(binned_mass, r);
  }
}
/* -------------------- END OF Routines used for initializing the circumgalactic medium -------------------- */


   void grid::SetParticleAttributes(float *Attribute[]){
    if(ProcessorNumber != MyProcessorNumber)
	    return; //return if not on right proc
    //defining everything needed
    int x, y, z, index[this->NumberOfParticles];
    int DensNum, MetalNum, SNColourNum; 
    float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits, MassUnits;
    //get the mass units to convert for star particle lifetime
    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, Time);
    MassUnits = DensityUnits * pow(LengthUnits, 3); 
    //allocate field for total density (density field + mass from initial dark matter particles)
    //for Hopkins-like initial star maker particles w/o calling fortran routine
    int size = GridDimension[0]*GridDimension[1]*GridDimension[2]; 
    float *dens_tot = new float[size]; 
    DensNum = FindField(Density, FieldType, NumberOfBaryonFields); 
    MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields); 
    SNColourNum = FindField(SNColour, FieldType, NumberOfBaryonFields); 
    //first the total density is just the same as in BaryonField
    for(int i = 0; i < size; i++)
	   dens_tot[i] = BaryonField[DensNum][i]; 
    for(int i = 0; i < NumberOfParticles; i++){
	   if(ParticleType[i] != PARTICLE_TYPE_DARK_MATTER)
		   continue; //if it's not dark matter we don't care yet 
	   if(ParticlePosition[0][i] < GridLeftEdge[0] || ParticlePosition[0][i] > GridRightEdge[0] || ParticlePosition[1][i] < GridLeftEdge[1] || ParticlePosition[1][i] > GridRightEdge[1] || ParticlePosition[2][i] < GridLeftEdge[2] || ParticlePosition[2][i] > GridRightEdge[2])
		   continue; //particle off grid
	    //Get the correct index in total density field from particle's position 
	    int ind; 
	    x = ((ParticlePosition[0][i] - GridLeftEdge[0]) / CellWidth[0][0]) + NumberOfGhostZones; 
	    y = ((ParticlePosition[1][i] - GridLeftEdge[1]) / CellWidth[1][0]) + NumberOfGhostZones;
	    z = ((ParticlePosition[2][i] - GridLeftEdge[2]) / CellWidth[2][0]) + NumberOfGhostZones; 
	    ind = (z * GridDimension[1] + y) * GridDimension[0] + x;
	    //add mass of dark matter particle uniformly over cell
	    dens_tot[ind] += (ParticleMass[i] * SolarMass / MassUnits) / (CellWidth[0][0]*CellWidth[1][0]*CellWidth[2][0]); 
	    ParticleAttribute[0][i] = -1.0; //Birthtime negative makes flagging agora particles easy
	    ParticleAttribute[1][i] = 9000*Myr_s; //just some big lifetime. probably a better way to handle this
	    //leave metallicity undefined bc duh 
    } 
    for(int ind = 0; ind < this->NumberOfParticles; ind++){
	    if(ParticleType[ind] == PARTICLE_TYPE_DARK_MATTER)
		    continue; //now we only care if it's an actual star particle
	   if(ParticlePosition[0][ind] < GridLeftEdge[0] || ParticlePosition[0][ind] > GridRightEdge[0] || ParticlePosition[1][ind] < GridLeftEdge[1] || ParticlePosition[1][ind] > GridRightEdge[1] || ParticlePosition[2][ind] < GridLeftEdge[2] || ParticlePosition[2][ind] > GridRightEdge[2])
		   continue; //particle off grid
	    //get index in grid where star particle is located 
	    x = ((ParticlePosition[0][ind] - GridLeftEdge[0]) / CellWidth[0][0]) + NumberOfGhostZones; 
	    y = ((ParticlePosition[1][ind] - GridLeftEdge[1]) / CellWidth[1][0]) + NumberOfGhostZones;
	    z = ((ParticlePosition[2][ind] - GridLeftEdge[2]) / CellWidth[2][0]) + NumberOfGhostZones; 
	    index[ind] = (z * GridDimension[1] + y) * GridDimension[0] + x;
	    float p_metal_density;
	    //get metal density of cell
	    if(SNColourNum == -1)
	    	p_metal_density = BaryonField[MetalNum][index[ind]];
	    else
		p_metal_density = BaryonField[MetalNum][index[ind]] + BaryonField[SNColourNum][index[ind]];
	    if(DensNum != -1){ 
		//set birthtime, lifetime, and metallicity of star particle 
	    	float p_density = BaryonField[DensNum][index[ind]];
	    	float p_metal_frac = p_metal_density / p_density;  
	    	ParticleAttribute[0][ind] = -1.0; //Birthtime negative makes flagging agora particles easy
	    	ParticleAttribute[1][ind] = max(StarMakerMinimumDynamicalTime*yr_s, POW(3*PI / (32*GravitationalConstant*dens_tot[index[ind]]), 0.5)*TimeUnits);
	    	ParticleAttribute[2][ind] = p_metal_frac;  
	    }
    }
   delete[] dens_tot; //don't need this field for anything else (and if we did, it's going to be out of scope by this point)
   } 


int InitializeParticles(grid *thisgrid, FLOAT *Center){ //, HierarchyEntry &TopGrid, TopGridData &MetaData, FLOAT * Center){

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
    printf("ID: %d count %d \n", thisgrid->GetGridID(), count); 
    thisgrid->SetNumberOfParticles(count);
    thisgrid->SetParticlePointers(Mass, Number, Type, Position,
				  Velocity, Attribute);
    thisgrid->SetParticleAttributes(Attribute);
    thisgrid->ExtraFunction("In Particle Creation");
    //MetaData.NumberOfParticles = count;

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
    double total_mass = 0.0; 

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
      total_mass += mass*1e9; 
    }

    std::cout << fname << " " << total_mass << std::endl; 

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
else if (i == VCIRC_TABLE_LENGTH){
	std::cout << "Falling off table at radius " << radius << " max radius " << VCircRadius[0] << std::endl; 
  	ENZO_FAIL("Fell off the circular velocity interpolation table");
}

// we know the radius is between i and i-1
return VCircVelocity[i-1] +
  (VCircVelocity[i] - VCircVelocity[i-1]) *
  (radius - VCircRadius[i-1])  /
  (VCircRadius[i] - VCircRadius[i-1]);
}
