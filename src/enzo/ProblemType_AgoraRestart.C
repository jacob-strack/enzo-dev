/***********************************************************************
/
/  Agora isolated galaxy restart
/
/  written by: Nathan Goldbaum
/  date:       March, 2013
/
/  PURPOSE:
/  https://sites.google.com/site/projectagoraworkspace/metagroup1/group2
/  https://www.dropbox.com/sh/1xzt1rysy9v3a9l/AAAHZyjrfTz88aG12H0Q_Rqla
/
************************************************************************/

#ifdef NEW_PROBLEM_TYPES
#include <stdio.h>
#include <iostream>
#include "preincludes.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "ProblemType.h"
#include "EventHooks.h"
#include "phys_constants.h"


#define VCIRC_TABLE_LENGTH 10000
#define KEV_PER_ERG (6.242e8)

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
void mt_init(unsigned_int seed);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);
inline int nlines(const char* fname);
void halo_init(struct CGMdata& CGM_data, grid* Grid, TopGridData &MetaData, FLOAT *binned_mass, int halo_type, float C, float Rstop=-1, int GasHalo_override=0);
double halo_S_of_r_Agora(double r, grid* Grid, FLOAT *binned_mass, TopGridData &MetaData);
double MassEnclosed_r(FLOAT *binned_mass, double rad, grid* Grid); 
double halo_dP_dr_Agora(double r, double P, grid* Grid, FLOAT *binned_mass, TopGridData &MetaData);
double halo_mod_g_of_r(double r, FLOAT *binned_mass, grid* Grid);
float HaloGasDensity(FLOAT R, struct CGMdata& CGM_data, grid* Grid);
float HaloGasTemperature(FLOAT R, struct CGMdata& CGM_data, grid* Grid);
void setup_chem(float density, float temperature, int equilibrate,
		float& DEdest, float&  HIdest, float& HIIdest,
		float& HeIdest, float& HeIIdest, float& HeIIIdest,
		float& HMdest, float& H2Idest, float& H2IIdest,
		float& DIdest, float& DIIdest, float& HDIdest);
int ReadEquilibriumTable(char * name, FLOAT Time);


int nlines(const char* fname) {

  FILE* fptr = fopen(fname, "r");
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

class ProblemType_AgoraRestart;

class AgoraRestartGrid : private grid
{
  friend class ProblemType_AgoraRestart;
};

class ProblemType_AgoraRestart : public EnzoProblemType
{
private:
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  FLOAT CenterPosition[MAX_DIMENSION];
  float Bfield[MAX_DIMENSION];
  FLOAT ScaleLength;
  FLOAT ScaleHeight;
  float DiskMass;
  float GasFraction;
  float DiskTemperature;
  float DiskMetallicity;
  float HaloMass;
  float HaloTemperature;
  float HaloMetallicity;
  FLOAT VCircRadius[VCIRC_TABLE_LENGTH];
  float VCircVelocity[VCIRC_TABLE_LENGTH];
  int RefineAtStart;
  int AgoraRestartGasHalo = 1; 
  float AgoraRestartGasHaloRatio = 10.;
  float AgoraRestartDMConcentration = 10.; 

public:
  ProblemType_AgoraRestart() : EnzoProblemType()
  {
    if (MyProcessorNumber == 0)
      std::cout << "Creating problem type Agora Restart" << std::endl;
  }

  ~ProblemType_AgoraRestart() {}

  virtual int InitializeFromRestart(
    HierarchyEntry &TopGrid, TopGridData &MetaData)
  {
    return SUCCESS;
  }

  virtual int InitializeSimulation(
    FILE *fptr, FILE *Outfptr,
    HierarchyEntry &TopGrid, TopGridData &MetaData)
  {
    if(debug)
    {
      printf("Entering AgoraRestartInitialize\n");
      fflush(stdout);
    }

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
    char *MetalName = "Metal_Density";
    char *MetalSNIaName = "MetalSNIa_Density";
    char *MetalSNIIName = "MetalSNII_Density";
    char *BxName = "Bx";
    char *ByName = "By";
    char *BzName = "Bz";
    char *PhiName = "Phi";

    /* local declarations */

    char line[MAX_LINE_LENGTH];
    int  i, ret, level;

    /* make sure it is 3D */

    if (MetaData.TopGridRank != 3)
    {
      printf("Cannot do AcoraRestart in %"ISYM" dimension(s)\n",
	     MetaData.TopGridRank);
      ENZO_FAIL("Agora Restart simulations must be 3D!");
    }

    for (i=0; i < MAX_DIMENSION; i++)
    {
      this->CenterPosition[i] = 0.5;
      this->Bfield[i] = 0.;
    }

    // These come from Oscar's sample output.  The units are:
    // Velocity: km/s
    // Mass: 10^9 Msun
    // Length: kpc
    // Temperature: K
    this->ScaleLength         = .0343218;
    this->ScaleHeight         = .00343218;
    this->DiskMass            = 42.9661;
    this->GasFraction         = 0.2;
    this->DiskTemperature     = 1e4;
    this->DiskMetallicity     = 0.0;
    this->HaloMass            = 0.10000;
    this->HaloTemperature     = this->DiskTemperature;
    this->HaloMetallicity     = 0.0;
    this->RefineAtStart       = TRUE;

    // set this from global data (kind of a hack)
    TestProblemData.MultiSpecies = MultiSpecies;
    int AgoraRestartGasHalo = 0; 
    /* read input from file */
    while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL)
    {
      ret = 0;
      ret += sscanf(line, "AgoraRestartCenterPosition = %"PSYM" %"PSYM" %"PSYM,
		    CenterPosition, CenterPosition+1, CenterPosition+2);
      ret += sscanf(line, "AgoraRestartScaleLength = %"PSYM, &ScaleLength);
      ret += sscanf(line, "AgoraRestartScaleHeight = %"PSYM, &ScaleHeight);
      ret += sscanf(line, "AgoraRestartDiskMass = %"FSYM, &DiskMass);
      ret += sscanf(line, "AgoraRestartGasFraction = %"FSYM, &GasFraction);
      ret += sscanf(line, "AgoraRestartDiskTemperature = %"FSYM,
		    &DiskTemperature);
      ret += sscanf(line, "AgoraRestartDiskMetallicity = %"FSYM,
		    &DiskMetallicity);
      ret += sscanf(line, "AgoraRestartHaloMass = %"FSYM, &HaloMass);
      ret += sscanf(line, "AgoraRestartHaloTemperature = %"FSYM,
		    &HaloTemperature);
      ret += sscanf(line, "AgoraRestartHaloMetallicity = %"FSYM,
                    &HaloMetallicity);
      ret += sscanf(line, "AgoraRestartMagneticField = %"FSYM" %"FSYM" %"FSYM,
		    Bfield, Bfield+1, Bfield+2);

      ret += sscanf(line, "AgoraRestartRefineAtStart = %"ISYM,
		    &RefineAtStart);
      ret += sscanf(line, "AgoraRestartHydrogenFractionByMass = %"FSYM,
		    &TestProblemData.HydrogenFractionByMass);
      ret += sscanf(line, "AgoraRestartHeliumFractionByMass = %"FSYM,
		    &TestProblemData.HeliumFractionByMass);
      ret += sscanf(line, "AgoraRestartMetalFractionByMass = %"FSYM,
		    &TestProblemData.MetalFractionByMass);
      ret += sscanf(line, "AgoraRestartDeuteriumToHydrogenRatio = %"FSYM,
		    &TestProblemData.DeuteriumToHydrogenRatio);
      ret += sscanf(line, "AgoraRestartInitialHIFraction  = %"FSYM,
		    &TestProblemData.HI_Fraction);
      ret += sscanf(line, "AgoraRestartInitialHIIFraction  = %"FSYM,
		    &TestProblemData.HII_Fraction);
      ret += sscanf(line, "AgoraRestartInitialHeIFraction  = %"FSYM,
		    &TestProblemData.HeI_Fraction);
      ret += sscanf(line, "AgoraRestartInitialHeIIFraction  = %"FSYM,
		    &TestProblemData.HeII_Fraction);
      ret += sscanf(line, "AgoraRestartInitialHeIIIFraction  = %"FSYM,
		    &TestProblemData.HeIII_Fraction);
      ret += sscanf(line, "AgoraRestartInitialHMFraction  = %"FSYM,
		    &TestProblemData.HM_Fraction);
      ret += sscanf(line, "AgoraRestartInitialH2IFraction  = %"FSYM,
		    &TestProblemData.H2I_Fraction);
      ret += sscanf(line, "AgoraRestartInitialH2IIFraction  = %"FSYM,
		    &TestProblemData.H2II_Fraction);
      ret += sscanf(line, "AgoraRestartInitialDIFraction  = %"FSYM,
		    &TestProblemData.DI_Fraction);
      ret += sscanf(line, "AgoraRestartInitialDIIFraction  = %"FSYM,
		    &TestProblemData.DII_Fraction);
      ret += sscanf(line, "AgoraRestartInitialHDIFraction  = %"FSYM,
		    &TestProblemData.HDI_Fraction);
      ret += sscanf(line, "AgoraRestartUseMetallicityField  = %"ISYM,
		    &TestProblemData.UseMetallicityField);
      ret += sscanf(line, "AgoraRestartGasHalo  = %"ISYM,
		    &AgoraRestartGasHalo);


      if (ret == 0 && strstr(line, "=") &&
	  (strstr(line, "AgoraRestart") || strstr(line, "TestProblem")) &&
	  line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR)
	fprintf(stderr,
		"*** warning: the following parameter line from AgoraRestart was not interpreted:\n%s\n",
		line);

    } // end input from parameter file
    
    if(AgoraRestartGasHalo)
	    std::cout << "Gas Halo ON" << std::endl;
    
    // Read in circular velocity table

    this->ReadInVcircData();




    /* set up top grid */

    float dummy_density = 1.0;
    float dummy_gas_energy = 1.0; // Only used if DualEnergyFormalism is True
    float dummy_total_energy = 1.0;
    float dummy_velocity[3] = {0.0, 0.0, 0.0};
    float dummy_b_field[3] = {1e-20, 1e-20, 1e-20}; // Only set if HydroMethod = mhd_rk

    if (this->InitializeUniformGrid(
	  TopGrid.GridData, dummy_density, dummy_total_energy,
	  dummy_gas_energy, dummy_velocity, dummy_b_field) == FAIL)
    {
      ENZO_FAIL("Error in InitializeUniformGrid");
    }

    this->InitializeParticles(TopGrid.GridData, TopGrid, MetaData);
    
    //bin mass for hydrostatic halo 
    float binned_mass[100]; 
    for(int ind = 0; ind < 100; ind++)
	    binned_mass[ind] = 0.0; 
    
    AgoraRestartGrid *thisgrid =
      static_cast<AgoraRestartGrid *>(TopGrid.GridData);
    // loop through the particles and deposit mass 
    for(int p = 0; p < MetaData.NumberOfParticles; p++){
   	float x,y,z;
	float ppos[3]; 
	thisgrid->ReturnParticlePosition(p,ppos); 
	x = ppos[0] - CenterPosition[0];
	y = ppos[1] - CenterPosition[1];
	z = ppos[2] - CenterPosition[2];
	float ParticleMass = TopGrid.GridData->ReturnParticleMass(p);
	float p_dens = ParticleMass*TopGrid.GridData->GetCellWidth(0,0)*TopGrid.GridData->GetCellWidth(1,0)*TopGrid.GridData->GetCellWidth(2,0); //code mass 
	float r_sph = sqrt(POW(fabs(x),2) + POW(fabs(y),2) + POW(fabs(z),2)); 
	float delta_r = sqrt(3.0) / 100; 
	int ind_r = int(r_sph / delta_r); //code_length
	if(ind_r >= 100)
		ENZO_FAIL("Bad index"); 	
	binned_mass[ind_r] += p_dens; 
    }
    //Now make binned_mass a cumulative sum in radius 
    float total_mass_enc = 0.0;
    for(int i = 0; i < 100; i++){
	total_mass_enc += binned_mass[i]; 
	binned_mass[i] = total_mass_enc; 
    }
    
    ReadEquilibriumTable("equilibrium_table_60_030-Zsun.h5", MetaData.Time);
    this->InitializeGrida(TopGrid.GridData, TopGrid, MetaData); //setup baryons. needed for CGM setup

    //fill CGM data here to be used to add halo later. Adding here means one integration for entire domain. 
    struct CGMdata CGM_data(8192);
    halo_init(CGM_data, thisgrid, MetaData, binned_mass, 6, 10); 
    if(AgoraRestartGasHalo) 
    	this->InitializeGridb(TopGrid.GridData, TopGrid, MetaData, binned_mass, CGM_data);

    /* Convert minimum initial overdensity for refinement to mass
       (unless MinimumMass itself was actually set). */

    if (MinimumMassForRefinement[0] == FLOAT_UNDEFINED) {
      MinimumMassForRefinement[0] = MinimumOverDensityForRefinement[0];
      for (int dim = 0; dim < MetaData.TopGridRank; dim++)
	MinimumMassForRefinement[0] *=(DomainRightEdge[dim]-DomainLeftEdge[dim])/
	  float(MetaData.TopGridDims[dim]);
    }

    /* If requested, refine the grid to the desired level. */

    
    
    if (RefineAtStart)
    {
      /* Declare, initialize, and fill out the first level of the LevelArray. */
      LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
      for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
	LevelArray[level] = NULL;
      AddLevel(LevelArray, &TopGrid, 0);

      /* Add levels to the maximum depth or until no new levels are created,
	 and re-initialize the level after it is created. */
      for (level = 0; level < MaximumRefinementLevel; level++) {
	if (RebuildHierarchy(&MetaData, LevelArray, level) == FAIL) {
	  fprintf(stderr, "Error in RebuildHierarchy.\n");
	  return FAIL;
	}
	if (LevelArray[level+1] == NULL)
	  break;
	LevelHierarchyEntry *Temp = LevelArray[level+1];
	while (Temp != NULL) {
	  if (this->InitializeGrida(Temp->GridData, TopGrid, MetaData) == FAIL)
	  {
	    ENZO_FAIL("Error in AgoraRestart->InitializeGrida");
	  }
	  if(AgoraRestartGasHalo){
		  if (this->InitializeGridb(Temp->GridData, TopGrid, MetaData, binned_mass, CGM_data) == FAIL)
		  {
		    ENZO_FAIL("Error in AgoraRestart->InitializeGridb");
		  }
	  }
	  Temp = Temp->NextGridThisLevel;
	} // end: loop over grids on this level
      } // end: loop over levels
    }


  // If we used the Equilibrium Table, delete it
  if (1){
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

    /* set up field names and units */
    int count = 0;
    DataLabel[count++] = DensName;
    DataLabel[count++] = Vel1Name;
    if(MetaData.TopGridRank > 1)
      DataLabel[count++] = Vel2Name;
    if(MetaData.TopGridRank > 2)
      DataLabel[count++] = Vel3Name;
    DataLabel[count++] = TEName;
    if (DualEnergyFormalism)
      DataLabel[count++] = GEName;

    if (HydroMethod == MHD_RK) {
      DataLabel[count++] = (char*) BxName;
      DataLabel[count++] = (char*) ByName;
      DataLabel[count++] = (char*) BzName;
      DataLabel[count++] = (char*) PhiName;
    }

    if (MultiSpecies)
    {
      DataLabel[count++] = ElectronName;
      DataLabel[count++] = HIName;
      DataLabel[count++] = HIIName;
      DataLabel[count++] = HeIName;
      DataLabel[count++] = HeIIName;
      DataLabel[count++] = HeIIIName;
      if (MultiSpecies > 1)
      {
	DataLabel[count++] = HMName;
	DataLabel[count++] = H2IName;
	DataLabel[count++] = H2IIName;
      }
      if (MultiSpecies > 2)
      {
	DataLabel[count++] = DIName;
	DataLabel[count++] = DIIName;
	DataLabel[count++] = HDIName;
      }
    }
    if (TestProblemData.UseMetallicityField)
      DataLabel[count++] = MetalName;
    if (StarMakerTypeIaSNe)
        DataLabel[count++] = MetalSNIaName;
    if (StarMakerTypeIISNeMetalField)
        DataLabel[count++] = MetalSNIIName;
    for (i = 0; i < count; i++)
      DataUnits[i] = NULL;



    if (MyProcessorNumber == ROOT_PROCESSOR)
    {
      fprintf(Outfptr, "AgoraRestartCenterPosition          = %"
	      PSYM" %"PSYM" %"PSYM"\n",
	      CenterPosition[0], CenterPosition[1], CenterPosition[2]);
      fprintf(Outfptr, "AgoraRestartMagneticField           = %"FSYM" %"FSYM" %"FSYM,
		    Bfield[0], Bfield[1], Bfield[2]);
      fprintf(Outfptr, "AgoraRestartScaleLength             = %"PSYM"\n",
	      ScaleLength);
      fprintf(Outfptr, "AgoraRestartScaleHeight             = %"PSYM"\n",
	      ScaleHeight);
      fprintf(Outfptr, "AgoraRestartDiskMass                = %"FSYM"\n",
	      DiskMass);
      fprintf(Outfptr, "AgoraRestartGasFraction             = %"FSYM"\n",
	      GasFraction);
      fprintf(Outfptr, "AgoraRestartDiskTemperature         = %"FSYM"\n",
	      DiskTemperature);
      fprintf(Outfptr, "AgoraRestartHaloMass                = %"FSYM"\n",
	      HaloMass);
      fprintf(Outfptr, "AgoraRestartHaloTemperature         = %"FSYM"\n",
	      HaloTemperature);
      fprintf(Outfptr, "AgoraRestartRefineAtStart           = %"ISYM"\n",
	      RefineAtStart);
      fprintf(Outfptr, "AgoraRestartHydrogenFractionByMass = %"FSYM"\n",
	      TestProblemData.HydrogenFractionByMass);
      fprintf(Outfptr, "AgoraRestartHeliumFractionByMass = %"FSYM"\n",
	      TestProblemData.HeliumFractionByMass);
      fprintf(Outfptr, "AgoraRestartMetalFractionByMass = %"FSYM"\n",
	      TestProblemData.MetalFractionByMass);
      fprintf(Outfptr, "AgoraRestartInitialHIFraction  = %"FSYM"\n",
	      TestProblemData.HI_Fraction);
      fprintf(Outfptr, "AgoraRestartInitialHIIFraction  = %"FSYM"\n",
	      TestProblemData.HII_Fraction);
      fprintf(Outfptr, "AgoraRestartInitialHeIFraction  = %"FSYM"\n",
	      TestProblemData.HeI_Fraction);
      fprintf(Outfptr, "AgoraRestartInitialHeIIFraction  = %"FSYM"\n",
	      TestProblemData.HeII_Fraction);
      fprintf(Outfptr, "AgoraRestartInitialHeIIIIFraction  = %"FSYM"\n",
	      TestProblemData.HeIII_Fraction);
      fprintf(Outfptr, "AgoraRestartInitialHMFraction  = %"FSYM"\n",
	      TestProblemData.HM_Fraction);
      fprintf(Outfptr, "AgoraRestartInitialH2IFraction  = %"FSYM"\n",
	      TestProblemData.H2I_Fraction);
      fprintf(Outfptr, "AgoraRestartInitialH2IIFraction  = %"FSYM"\n",
	      TestProblemData.H2II_Fraction);
      fprintf(Outfptr, "AgoraRestartInitialDIFraction  = %"FSYM"\n",
	      TestProblemData.DI_Fraction);
      fprintf(Outfptr, "AgoraRestartInitialDIIFraction  = %"FSYM"\n",
	      TestProblemData.DII_Fraction);
      fprintf(Outfptr, "AgoraRestartInitialHDIFraction  = %"FSYM"\n",
	      TestProblemData.HDI_Fraction);
      fprintf(Outfptr, "AgoraRestartUseMetallicityField  = %"ISYM"\n",
	      TestProblemData.UseMetallicityField);
    }
 
    return SUCCESS;

  } // InitializeSimulation

  int InitializeGrida(grid *thisgrid_orig, HierarchyEntry &TopGrid, TopGridData &MetaData){
    if(debug)
      printf("Entering AgoraRestart InitializeGrida\n");

    AgoraRestartGrid *thisgrid =
      static_cast<AgoraRestartGrid *>(thisgrid_orig);

    if (thisgrid->ProcessorNumber != MyProcessorNumber)
      return SUCCESS;

    /* Get units */
    float DensityUnits=1, LengthUnits=1, VelocityUnits=1, TimeUnits=1,
      TemperatureUnits=1;
    double MassUnits=1;

    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		 &TimeUnits, &VelocityUnits, &MassUnits, thisgrid->Time) == FAIL) {
      ENZO_FAIL("Error in GetUnits.");
    }

    /* Identify physical quantities */
    int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num, PhiNum, MetalNum;

    int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;

    if (thisgrid->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
					     Vel3Num, TENum, B1Num, B2Num, B3Num, PhiNum) == FAIL) {
      fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
      ENZO_FAIL("");
    }

    if (TestProblemData.MultiSpecies)
      if (thisgrid->IdentifySpeciesFields(
	    DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
	    HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL)
	ENZO_FAIL("Error in grid->IdentifySpeciesFields.");

    int MetallicityField = FALSE;
    if ((MetalNum = FindField(
	   Metallicity, thisgrid->FieldType, thisgrid->NumberOfBaryonFields)
	  ) != -1)
      MetallicityField = TRUE;
    else
      MetalNum = 0;

    int dim, i, j, k, size, index=0;
    float RhoZero, DiskGasEnergy, DiskDensity, HaloGasEnergy, HaloDensity,
      BoxVolume, vcirc, mu;
    FLOAT x, y, z, radius, xy_radius, cellwidth;

    /* Compute size of this grid */
    size = 1;
    for (dim = 0; dim < thisgrid->GridRank; dim++)
      size *= thisgrid->GridDimension[dim];
    cellwidth = thisgrid->CellWidth[0][0];

    /* Compute the size of the box */
    BoxVolume = 1.;
    for (dim = 0; dim < TopGrid.GridData->GetGridRank(); dim++)
      BoxVolume *= (DomainRightEdge[dim] - DomainLeftEdge[dim]);

    /* Find the mean molecular weight */

    if (TestProblemData.MultiSpecies == FALSE)
      mu = Mu;
    else
    {
      // Atomic hydrogen
      mu = TestProblemData.HydrogenFractionByMass *
	(TestProblemData.HI_Fraction + 2.0*TestProblemData.HII_Fraction);

      // Helium
      mu += TestProblemData.HeliumFractionByMass / 4.0 *
	(TestProblemData.HeI_Fraction + 2.0*TestProblemData.HeII_Fraction +
	 3.0*TestProblemData.HeIII_Fraction);

      // Molecular hydrogen, ignore Deuterium
      if (TestProblemData.MultiSpecies > 1)
	mu += TestProblemData.HydrogenFractionByMass / 2.0 *
	  (TestProblemData.H2I_Fraction + 2.0*TestProblemData.H2II_Fraction);

      // Metals
      if (TestProblemData.UseMetallicityField)
	mu += TestProblemData.MetalFractionByMass / 16.0;

      mu = POW(mu, -1);

    }

    /* Find global physical properties */
    RhoZero = this->DiskMass * this->GasFraction / (4.*pi) /
      (POW((this->ScaleLength),2)*(this->ScaleHeight));

    HaloGasEnergy = this->HaloTemperature / mu / (Gamma - 1) /
      TemperatureUnits;

    HaloDensity = this->HaloMass / BoxVolume;

    DiskGasEnergy = this->DiskTemperature / mu / (Gamma - 1) /
      TemperatureUnits;

    /* Loop over the mesh. */
    float temperature; 

    for (k = 0; k < thisgrid->GridDimension[2]; k++)
    {
      for (j = 0; j < thisgrid->GridDimension[1]; j++)
      {
	for (i = 0; i < thisgrid->GridDimension[0]; i++, index++)
	{
	  /* Compute position */

	  x = (thisgrid->CellLeftEdge[0][i] + 0.5*thisgrid->CellWidth[0][i]) *
	    LengthUnits;
	  y = (thisgrid->CellLeftEdge[1][j] + 0.5*thisgrid->CellWidth[1][j]) *
	    LengthUnits;
	  z = (thisgrid->CellLeftEdge[2][k] + 0.5*thisgrid->CellWidth[2][k]) *
	    LengthUnits;

	  x -= this->CenterPosition[0]*LengthUnits;
	  y -= this->CenterPosition[1]*LengthUnits;
	  z -= this->CenterPosition[2]*LengthUnits;

	  radius = sqrt(POW(x, 2) +
			POW(y, 2) +
			POW(z, 2) );

	  xy_radius = sqrt(POW(x, 2) +
			   POW(y, 2) );

	  /* Find disk density, halo density and internal energy */

	  DiskDensity = gauss_mass(RhoZero, x/LengthUnits, y/LengthUnits,
				   z/LengthUnits, cellwidth) / POW(cellwidth, 3);
	  if ((HaloDensity*HaloTemperature > DiskDensity*DiskTemperature))
	  {
	    thisgrid->BaryonField[DensNum][index] = 1e-31/DensityUnits; //HaloDensity; //a low background density
	    thisgrid->BaryonField[TENum][index] = HaloGasEnergy;
	    if (DualEnergyFormalism)
	      thisgrid->BaryonField[GENum][index] = HaloGasEnergy;
	    temperature = HaloTemperature; //I guess? this is a background temperature for here. It can't be the pNFW halo model yet bc I need this temperature for chem init. 


	    thisgrid->BaryonField[Vel1Num][index] = 0; //no halo rotation
	    thisgrid->BaryonField[Vel2Num][index] = 0;
	    thisgrid->BaryonField[Vel3Num][index] = 0;

	    if (TestProblemData.UseMetallicityField)
	      thisgrid->BaryonField[MetalNum][index] = thisgrid->BaryonField[DensNum][index] *
		TestProblemData.MetalFractionByMass * HaloMetallicity;
	  }
	  else // Ok, we're in the disk
	  {
	    thisgrid->BaryonField[DensNum][index] = DiskDensity;
	    thisgrid->BaryonField[DensNum][index] = 1e-31/DensityUnits; //HaloDensity; //a low background density
	    vcirc = this->InterpolateVcircTable(xy_radius);
	    

	    thisgrid->BaryonField[Vel1Num][index] =
	      -vcirc*y/xy_radius/VelocityUnits;
	    thisgrid->BaryonField[Vel2Num][index] =
	      vcirc*x/xy_radius/VelocityUnits;
	    thisgrid->BaryonField[Vel3Num][index] = 0;

	    thisgrid->BaryonField[TENum][index] = DiskGasEnergy; 
	    if(HydroMethod != Zeus_Hydro) {
	      thisgrid->BaryonField[TENum][index] +=  0.5 *
		(POW(thisgrid->BaryonField[Vel1Num][index],2) +
		 POW(thisgrid->BaryonField[Vel2Num][index],2) +
		 POW(thisgrid->BaryonField[Vel3Num][index],2));
	    }
	    
	    if (DualEnergyFormalism)
	      {
		thisgrid->BaryonField[GENum][index] = DiskGasEnergy;
	      }

	    if (TestProblemData.UseMetallicityField) {
	      thisgrid->BaryonField[MetalNum][index] = thisgrid->BaryonField[DensNum][index] *
		TestProblemData.MetalFractionByMass * DiskMetallicity;

	    }
	    temperature = DiskTemperature;
	  }


      if (StarMakerTypeIaSNe) {
          int SNIaNum = FindField(MetalSNIaDensity , thisgrid->FieldType, thisgrid->NumberOfBaryonFields);
          if(SNIaNum != -1) {
              thisgrid->BaryonField[SNIaNum][index] = thisgrid->BaryonField[DensNum][index] * 1.0e-8;
          }
      }
      if (StarMakerTypeIISNeMetalField) {
          int SNIINum = FindField(MetalSNIIDensity , thisgrid->FieldType, thisgrid->NumberOfBaryonFields);
          if(SNIINum != -1) {
              thisgrid->BaryonField[SNIINum][index] = thisgrid->BaryonField[DensNum][index] * 1.0e-8;
          }
          else {
              ENZO_FAIL("Thought we would find a SNII field but did not.");
          }
      }
      if(1){ //init chem the way GalaxySimulation does with EquilibriumTable. 
	     //trying to be consistent with what is done in S(r) for gas halo. 
	  int EquilibrateChem = 1;

	  if (MultiSpecies == 3)
	    setup_chem(thisgrid->BaryonField[DensNum][index], temperature, EquilibrateChem,
		       thisgrid->BaryonField[DeNum][index], thisgrid->BaryonField[HINum][index], thisgrid->BaryonField[HIINum][index],
		       thisgrid->BaryonField[HeINum][index], thisgrid->BaryonField[HeIINum][index], thisgrid->BaryonField[HeIIINum][index],
		       thisgrid->BaryonField[HMNum][index], thisgrid->BaryonField[H2INum][index], thisgrid->BaryonField[H2IINum][index],
		       thisgrid->BaryonField[DINum][index], thisgrid->BaryonField[DIINum][index], thisgrid->BaryonField[HDINum][index]);
	  else if (MultiSpecies == 2) {
	    float temp;
	    setup_chem(thisgrid->BaryonField[DensNum][index], temperature, EquilibrateChem,
		       thisgrid->BaryonField[DeNum][index], thisgrid->BaryonField[HINum][index], thisgrid->BaryonField[HIINum][index],
		       thisgrid->BaryonField[HeINum][index], thisgrid->BaryonField[HeIINum][index], thisgrid->BaryonField[HeIIINum][index],
		       thisgrid->BaryonField[HMNum][index], thisgrid->BaryonField[H2INum][index], thisgrid->BaryonField[H2IINum][index],
		       temp, temp, temp);
	  }
	  else {
	    float temp;
	    setup_chem(thisgrid->BaryonField[DensNum][index], temperature, EquilibrateChem,
		       thisgrid->BaryonField[DeNum][index], thisgrid->BaryonField[HINum][index], thisgrid->BaryonField[HIINum][index],
		       thisgrid->BaryonField[HeINum][index], thisgrid->BaryonField[HeIINum][index], thisgrid->BaryonField[HeIIINum][index],
		       temp, temp, temp,
		       temp, temp, temp);
	  }

      }

	  if(TestProblemData.MultiSpecies && 0)
	  {
	    thisgrid->BaryonField[HINum][index] = TestProblemData.HI_Fraction *
	      TestProblemData.HydrogenFractionByMass * thisgrid->BaryonField[DensNum][index];

	    thisgrid->BaryonField[HIINum][index] = TestProblemData.HII_Fraction *
	      TestProblemData.HydrogenFractionByMass * thisgrid->BaryonField[DensNum][index];

	    thisgrid->BaryonField[HeINum][index] = TestProblemData.HeI_Fraction *
	      TestProblemData.HeliumFractionByMass * thisgrid->BaryonField[DensNum][index];

	    thisgrid->BaryonField[HeIINum][index] = TestProblemData.HeII_Fraction *
	      TestProblemData.HeliumFractionByMass * thisgrid->BaryonField[DensNum][index];

	    thisgrid->BaryonField[HeIIINum][index] = TestProblemData.HeIII_Fraction *
	      TestProblemData.HeliumFractionByMass * thisgrid->BaryonField[DensNum][index];

	    if(TestProblemData.MultiSpecies > 1){
	      thisgrid->BaryonField[HMNum][index] = TestProblemData.HM_Fraction *
		TestProblemData.HydrogenFractionByMass * thisgrid->BaryonField[DensNum][index];

	      thisgrid->BaryonField[H2INum][index] = 2 * TestProblemData.H2I_Fraction *
		TestProblemData.HydrogenFractionByMass * thisgrid->BaryonField[DensNum][index];

	      thisgrid->BaryonField[H2IINum][index] = 2 * TestProblemData.H2II_Fraction *
		TestProblemData.HydrogenFractionByMass * thisgrid->BaryonField[DensNum][index];
	    }

	    if (TestProblemData.MultiSpecies > 1)
	      thisgrid->BaryonField[HIINum][index] -=
		(thisgrid->BaryonField[HMNum][index] + thisgrid->BaryonField[H2IINum][index]
		 + thisgrid->BaryonField[H2INum][index]);

	    // Electron "density" (remember, this is a factor of m_p/m_e scaled
	    // from the 'normal' density for convenience) is calculated by
	    // summing up all of the ionized species.  The factors of 0.25 and
	    // 0.5 in front of HeII and HeIII are to fix the fact that we're
	    // calculating mass density, not number density (because the
	    // thisgrid->BaryonField values are 4x as heavy for helium for a single
	    // electron)
	    thisgrid->BaryonField[DeNum][index] = thisgrid->BaryonField[HIINum][index] +
	      0.25*thisgrid->BaryonField[HeIINum][index] +
	      0.5*thisgrid->BaryonField[HeIIINum][index];

	    if (TestProblemData.MultiSpecies > 1)
	      thisgrid->BaryonField[DeNum][index] += 0.5*thisgrid->BaryonField[H2IINum][index] -
		thisgrid->BaryonField[HMNum][index];

	    // Set deuterium species (assumed to be a negligible fraction of the
	    // total, so not counted in the conservation)
	    if(TestProblemData.MultiSpecies > 2){
	      thisgrid->BaryonField[DINum ][index] =
		CoolData.DeuteriumToHydrogenRatio * thisgrid->BaryonField[HINum][index];
	      thisgrid->BaryonField[DIINum][index] =
		CoolData.DeuteriumToHydrogenRatio * thisgrid->BaryonField[HIINum][index];
	      thisgrid->BaryonField[HDINum][index] = 0.75 *
		CoolData.DeuteriumToHydrogenRatio * thisgrid->BaryonField[H2INum][index];
	    }
	  } // if(TestProblemData.MultiSpecies)

	    if (HydroMethod == MHD_RK)
	      {
		thisgrid->BaryonField[B1Num][index] = Bfield[0];
		thisgrid->BaryonField[B2Num][index] = Bfield[1];
		thisgrid->BaryonField[B3Num][index] = Bfield[2];

		thisgrid->BaryonField[TENum][index] += 
		  0.5*(POW(thisgrid->BaryonField[B1Num][index], 2) +
		       POW(thisgrid->BaryonField[B2Num][index], 2) + 
		       POW(thisgrid->BaryonField[B3Num][index], 2))/thisgrid->BaryonField[DensNum][index];
	      }



	} // i
      } // j
    } // k

    return SUCCESS;

  }

  int InitializeGridb(grid *thisgrid_orig, HierarchyEntry &TopGrid,
		     TopGridData &MetaData, float* binned_mass, struct CGMdata& CGM_data)
  {
    AgoraRestartGrid *thisgrid =
      static_cast<AgoraRestartGrid *>(thisgrid_orig);

    if (thisgrid->ProcessorNumber != MyProcessorNumber)
      return SUCCESS;
    /* Get units */
    float DensityUnits=1, LengthUnits=1, VelocityUnits=1, TimeUnits=1,
      TemperatureUnits=1;
    double MassUnits=1;

    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		 &TimeUnits, &VelocityUnits, &MassUnits, thisgrid->Time) == FAIL) {
      ENZO_FAIL("Error in GetUnits.");
    }
    int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num, PhiNum, MetalNum;

    int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;

    if (thisgrid->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
					     Vel3Num, TENum, B1Num, B2Num, B3Num, PhiNum) == FAIL) {
      fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
      ENZO_FAIL("");
    }
    if (TestProblemData.MultiSpecies)
      if (thisgrid->IdentifySpeciesFields(
	    DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
	    HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL)
	ENZO_FAIL("Error in grid->IdentifySpeciesFields.");

    int MetallicityField = FALSE;
    if ((MetalNum = FindField(
	   Metallicity, thisgrid->FieldType, thisgrid->NumberOfBaryonFields)
	  ) != -1)
      MetallicityField = TRUE;
    else
      MetalNum = 0;

    //this function will esentially just lay down the halo. everything else should be taken care of by now.
    float ScaleLength         = .0343218;
    float ScaleHeight         = .00343218;
    float DiskMass            = 42.9661;
    float GasFraction         = 0.2;
    float DiskTemperature     = 1e4;
    float DiskMetallicity     = 0.0;
    float HaloMass            = 0.10000;
    float HaloTemperature     = DiskTemperature;
    float HaloMetallicity     = 0.0;
    float RhoZero = DiskMass * GasFraction / (4.*pi) /
      (POW((ScaleLength),2)*(ScaleHeight));
    float BoxVolume = 1.;
    for (int dim = 0; dim < TopGrid.GridData->GetGridRank(); dim++)
      BoxVolume *= (DomainRightEdge[dim] - DomainLeftEdge[dim]);
    float HaloDensity = HaloMass / BoxVolume;
    float DiskDensity = DiskMass / BoxVolume;
    float HaloGasEnergy = HaloTemperature / Mu / (Gamma - 1) /
      TemperatureUnits;
    int i,j,k,index = 0; 
    int size; 
    FLOAT x, y, z, radius, xy_radius, cellwidth;
    /* Compute size of this grid */
    size = 1;
    for (int dim = 0; dim < thisgrid->GridRank; dim++)
      size *= thisgrid->GridDimension[dim];
    cellwidth = thisgrid->CellWidth[0][0];
    for (k = 0; k < thisgrid->GridDimension[2]; k++)
    {
      for (j = 0; j < thisgrid->GridDimension[1]; j++)
      {
	for (i = 0; i < thisgrid->GridDimension[0]; i++, index++)
	{
	  /* Compute position */

	  x = (thisgrid->CellLeftEdge[0][i] + 0.5*thisgrid->CellWidth[0][i]) *
	    LengthUnits;
	  y = (thisgrid->CellLeftEdge[1][j] + 0.5*thisgrid->CellWidth[1][j]) *
	    LengthUnits;
	  z = (thisgrid->CellLeftEdge[2][k] + 0.5*thisgrid->CellWidth[2][k]) *
	    LengthUnits;

	  x -= this->CenterPosition[0]*LengthUnits;
	  y -= this->CenterPosition[1]*LengthUnits;
	  z -= this->CenterPosition[2]*LengthUnits;

	  radius = sqrt(POW(x, 2) +
			POW(y, 2) +
			POW(z, 2) );

	  xy_radius = sqrt(POW(x, 2) +
			   POW(y, 2) );
	  float HaloDensityPrecip = HaloGasDensity(radius, CGM_data, thisgrid);
	  float HaloTemperaturePrecip = HaloGasTemperature(radius, CGM_data, thisgrid);
	  if(HaloTemperaturePrecip < 0.0) continue; //we're above the stop radius, so nothing to do for CGM 
	  float HaloGasEnergyPrecip = HaloTemperaturePrecip / Mu / (Gamma - 1) /
      TemperatureUnits;
	  float HaloGasEnergy = HaloTemperature / Mu / (Gamma - 1) /
      TemperatureUnits;
	  //lay down the gas halo
          DiskDensity = gauss_mass(RhoZero, x/LengthUnits, y/LengthUnits,
				   z/LengthUnits, cellwidth) / POW(cellwidth, 3); 
	  if ( (HaloDensity*HaloTemperature > DiskDensity*DiskTemperature) )
		  {
		    thisgrid->BaryonField[DensNum][index] = HaloDensityPrecip;
		    thisgrid->BaryonField[TENum][index] -= HaloGasEnergy; //get rid of energy from overwritten placeholder density
		    thisgrid->BaryonField[TENum][index] += HaloGasEnergyPrecip;//but DO NOT get rid of the magnetic energy
		    if (DualEnergyFormalism)
		      thisgrid->BaryonField[GENum][index] = HaloGasEnergyPrecip;
		    if (TestProblemData.UseMetallicityField)
		      thisgrid->BaryonField[MetalNum][index] = thisgrid->BaryonField[DensNum][index] *
			TestProblemData.MetalFractionByMass * HaloMetallicity; 
		  }
	}
      }
    }

    return SUCCESS;

  } // InitializeGrid

  void InitializeParticles(grid *thisgrid_orig, HierarchyEntry &TopGrid,
			  TopGridData &MetaData)
  {
    AgoraRestartGrid *thisgrid =
      static_cast<AgoraRestartGrid *>(thisgrid_orig);

    mt_init(thisgrid->ID);

    if(debug)
      printf("Entering AgoraRestart InitializeParticles\n");

    // Determine the number of particles of each type
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
    for (int i = 0; i < thisgrid->GridRank; i++)
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

    FLOAT dx = thisgrid->CellWidth[0][0];

    // Read them in and assign them as we go
    int count = 0;
    this->ReadParticlesFromFile(
      Number, Type, Position, Velocity, Mass,
      "bulge.dat", PARTICLE_TYPE_STAR, count, dx);
    this->ReadParticlesFromFile(
      Number, Type, Position, Velocity, Mass,
      "disk.dat", PARTICLE_TYPE_STAR, count, dx);
    this->ReadParticlesFromFile(
      Number, Type, Position, Velocity, Mass,
      "halo.dat", PARTICLE_TYPE_DARK_MATTER, count, dx);

    thisgrid->SetNumberOfParticles(count);
    thisgrid->SetParticlePointers(Mass, Number, Type, Position,
				  Velocity, Attribute);
    MetaData.NumberOfParticles = count;
    if(debug) fprintf(stderr, "InitializeParticles: Set Number of Particles %"ISYM"\n", count);

  }
  
  void InitialParticlePositions(grid *thisgrid_orig, HierarchyEntry &TopGrid,
			  TopGridData &MetaData, float* ParticlePositions)
  {
    AgoraRestartGrid *thisgrid =
      static_cast<AgoraRestartGrid *>(thisgrid_orig);

    mt_init(thisgrid->ID);

    if(debug)
      printf("Entering AgoraRestart InitializeParticles\n");

    // Determine the number of particles of each type
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
    for (int i = 0; i < thisgrid->GridRank; i++)
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

    FLOAT dx = thisgrid->CellWidth[0][0];

    // Read them in and assign them as we go
    int count = 0;
    this->ReadParticlesFromFile(
      Number, Type, Position, Velocity, Mass,
      "bulge.dat", PARTICLE_TYPE_STAR, count, dx);
    this->ReadParticlesFromFile(
      Number, Type, Position, Velocity, Mass,
      "disk.dat", PARTICLE_TYPE_STAR, count, dx);
    this->ReadParticlesFromFile(
      Number, Type, Position, Velocity, Mass,
      "halo.dat", PARTICLE_TYPE_DARK_MATTER, count, dx);

    thisgrid->SetNumberOfParticles(count);
    thisgrid->SetParticlePointers(Mass, Number, Type, Position,
				  Velocity, Attribute);
    MetaData.NumberOfParticles = count;
    if(debug) fprintf(stderr, "InitializeParticles: Set Number of Particles %"ISYM"\n", count);

  }

  float gauss_mass(
    float RhoZero, FLOAT xpos, FLOAT ypos, FLOAT zpos, FLOAT cellwidth)
  {
    // Computes the total mass in a given cell by integrating the density
    // profile using 5-point Gaussian quadrature.
    // http://mathworld.wolfram.com/Legendre-GaussQuadrature.html
    FLOAT EvaluationPoints [5] = {-0.90617985,-0.53846931,0.0,0.53846931,0.90617985};
    FLOAT Weights [5] = {0.23692689,0.47862867,0.56888889,0.47862867,0.23692689};
    FLOAT xResult [5];
    FLOAT yResult [5];
    FLOAT r, z;
    float Mass = 0;
    int i,j,k;

    for (i=0;i<5;i++)
    {
      xResult[i] = 0.0;
      for (j=0;j<5;j++)
      {
	yResult[j] = 0.0;
	for (k=0;k<5;k++)
	{
	  r = sqrt((POW(xpos+EvaluationPoints[i]*cellwidth/2.0, 2.0) +
		    POW(ypos+EvaluationPoints[j]*cellwidth/2.0, 2.0) ) );
	  z = fabs(zpos+EvaluationPoints[k]*cellwidth/2.0);
	  yResult[j] +=
	    cellwidth/2.0 * Weights[k] * RhoZero *
	    PEXP(-r/this->ScaleLength) *
	    PEXP(-fabs(z)/this->ScaleHeight);
	}
	xResult[i] += cellwidth/2.0*Weights[j]*yResult[j];
      }
      Mass += cellwidth/2.0*Weights[i]*xResult[i];
    }
    return Mass;
  }

  void ReadInVcircData(void)
  {
    FILE *fptr;
    char line[MAX_LINE_LENGTH];
    int i=0, ret;
    float vcirc;
    FLOAT rad;

    fptr = fopen("vcirc.dat" , "r");

    while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL)
    {
      ret += sscanf(line, "%"PSYM" %"FSYM, &rad, &vcirc);
      this->VCircRadius[i] = rad*kpc_cm; // 3.08567758e21 = kpc/cm
      this->VCircVelocity[i] = vcirc*1e5; // 1e5 = (km/s)/(cm/s)
      i += 1;
    }

    fclose(fptr);
  } // ReadInVcircData

  float InterpolateVcircTable(FLOAT radius)
  {
    int i;

    for (i = 0; i < VCIRC_TABLE_LENGTH; i++)
      if (radius < this->VCircRadius[i])
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

  int ReadParticlesFromFile(PINT *Number, int *Type, FLOAT *Position[],
			    float *Velocity[], float* Mass, const char* fname,
			    Eint32 particle_type, int &c, FLOAT dx)
  {
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

      Position[0][c] = x * kpc_cm / LengthUnits + this->CenterPosition[0];
      Position[1][c] = y * kpc_cm / LengthUnits + this->CenterPosition[1];
      Position[2][c] = z * kpc_cm / LengthUnits + this->CenterPosition[2];

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
  } // ReadParticlesFromFile

}; // class declaration


//.. register:
namespace {
    EnzoProblemType_creator_concrete<ProblemType_AgoraRestart>
        agora_restart("AgoraRestart");
}


#endif



 void halo_init(struct CGMdata& CGM_data, grid* Grid, TopGridData &MetaData, FLOAT *binned_mass, int halo_type, float C, float Rstop, int GasHalo_override){

    /* Get units */
    float DensityUnits=1, LengthUnits=1, VelocityUnits=1, TimeUnits=1,
      TemperatureUnits=1;
    double MassUnits=1;

    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		 &TimeUnits, &VelocityUnits, &MassUnits, Grid->ReturnTime()) == FAIL) {
      ENZO_FAIL("Error in GetUnits.");
    }
  if (GasHalo_override) // not 0
    halo_type = GasHalo_override;
  
  double k1, k2, k3, k4;
  double M, R200, rho_crit = 1.8788e-29*0.49;
  double Rstart;

  int index;
  
  double MassUnitsDouble = double(DensityUnits)*POW(double(LengthUnits), 3.0);
  M = binned_mass[99] * MassUnitsDouble;  // DM halo total mass in CGS
  R200 = pow(3.0/(4.0*3.14159)*M/(200.*rho_crit),1./3.);  // virial radius in CGS
  if (Rstop < 0)
    Rstop = fabs(Rstop)*R200;
  CGM_data.R_outer = Rstop;// integrate out to the virial radius of halo

  Rstart = 0.01*LengthUnits; // you could force a different start if you liked

  // stepsize for RK4 integration and radial bins
  CGM_data.R_inner = Rstart;
  CGM_data.dr = (CGM_data.R_outer - CGM_data.R_inner)/ double(CGM_data.nbins); 
  

    /* Integrate pressure assuming HSE & S(r) from Voit 2019, then convert to n & T,
       instead of integrating n(r) directly as with methods 4 & 5. This makes the boundary
       condition easier to handle.*/
    double dr, rmax, vcirc2_max;
    double this_press, this_ent, this_radius;//, this_n;
    double mu_ratio = 1.17/Mu; // mu_e/mu
    double T_floor = 4e4; // IGM, was 4e4

    // boundary condition & quantities for integration
    dr = -1.0*CGM_data.dr;
    this_radius = R200;
    this_ent = halo_S_of_r_Agora(this_radius, Grid, binned_mass, MetaData); // in erg*cm^2
    rmax = 2.163*R200/C;
    vcirc2_max = GravConst * MassEnclosed_r(binned_mass,rmax,Grid)/rmax;
    this_press = mu_ratio*POW(0.25*Mu*mh*vcirc2_max/POW(this_ent, 1./Gamma),
		     Gamma/(Gamma-1.));
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
      k1 = halo_dP_dr_Agora(this_radius,          this_press,             Grid, binned_mass, MetaData);
      k2 = halo_dP_dr_Agora(this_radius + 0.5*dr, this_press + 0.5*dr*k1, Grid, binned_mass, MetaData);
      k3 = halo_dP_dr_Agora(this_radius + 0.5*dr, this_press + 0.5*dr*k2, Grid, binned_mass, MetaData);
      k4 = halo_dP_dr_Agora(this_radius + dr,     this_press + dr*k3,     Grid, binned_mass, MetaData);
      // update radius, pressure, entropy
      std::cout << "inward integration r " << this_radius << " " << dr << std::endl;
      this_radius += dr;  // new radius
      this_press += (1.0/6.0) * dr * (k1 + 2.0*k2 + 2.0*k3 + k4); // P @ new radius
      this_ent = halo_S_of_r_Agora(this_radius, Grid, binned_mass,MetaData); // entropy @ new radius
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
    this_ent = halo_S_of_r_Agora(this_radius, Grid, binned_mass,MetaData); // in erg*cm^2
    rmax = 2.163*R200/C;
    vcirc2_max = GravConst * MassEnclosed_r(binned_mass, rmax, Grid)/rmax;
    this_press = mu_ratio*POW(0.25*Mu*mh*vcirc2_max/POW(this_ent, 1./Gamma),
		     Gamma/(Gamma-1.));
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
    assert (y0 < 0.0);
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
      std::cout << "outward integration " << this_radius << " " << CGM_data.R_outer << std::endl;
      //this_dPdr = this_press/this_radius * dlP_dlr;
      k1 = halo_dP_dr_Agora(this_radius,          this_press,             Grid, binned_mass, MetaData);
      k2 = halo_dP_dr_Agora(this_radius + 0.5*dr, this_press + 0.5*dr*k1, Grid, binned_mass, MetaData);
      k3 = halo_dP_dr_Agora(this_radius + 0.5*dr, this_press + 0.5*dr*k2, Grid, binned_mass, MetaData);
      k4 = halo_dP_dr_Agora(this_radius + dr,     this_press + dr*k3,     Grid, binned_mass, MetaData);
      // update radius, pressure, entropy
      this_radius += dr;  // new radius
      this_press += (1.0/6.0) * dr * (k1 + 2.0*k2 + 2.0*k3 + k4); // P @ new radius
      this_ent = halo_S_of_r_Agora(this_radius, Grid, binned_mass,MetaData); // entropy @ new radius
      // update density and radius
      //this_dens = -2.0 * this_dPdr/(1.22*mh*halo_mod_g_of_r(this_radius, binned_mass)); // n_e = n_i
      //this_temp = POW(10, sigmoid(log10(this_radius), r0, k, y0, y_offset));
      //this_press = POW(10, dlP_dlr*log10(this_radius/R200) + log10(press_vir));
      // store everything in the struct
      index = int((this_radius - CGM_data.R_inner)/dr + 1.0e-3); 
      if (index < CGM_data.nbins) {
	CGM_data.n_rad[index] = 2 * POW(this_press/(mu_ratio*this_ent), 1./Gamma);
	CGM_data.T_rad[index] = POW(POW(this_press/mu_ratio, Gamma-1.) * this_ent, 1./Gamma) / kboltz;
	CGM_data.rad[index] = this_radius;
      }
      else
	break;
  }
    
  if (CGM_data.R_inner == 0) {
    // this integration acts a little squirrelly around r=0 because the mass values are garbage.  Cheap fix.
    CGM_data.rad[0]=CGM_data.rad[1];
    CGM_data.n_rad[0]=CGM_data.n_rad[1];
    CGM_data.T_rad[0]=CGM_data.T_rad[1];
  }
  
  return;
}

double halo_dP_dr_Agora(double r, double P, grid* Grid, FLOAT *binned_mass, TopGridData &MetaData) {
    double ret =  -1.0 * halo_mod_g_of_r(r, binned_mass,Grid) * 1.22 * mh * POW( P/(1.1/Mu) / halo_S_of_r_Agora(r,Grid, binned_mass,MetaData),
						    1./Gamma );
    if(halo_mod_g_of_r(r, binned_mass,Grid) < 0){
	    std::cout << halo_mod_g_of_r(r, binned_mass, Grid) << std::endl;
	    ENZO_FAIL("negative g"); 
    }
    if(ret > 0)
        ENZO_FAIL("positive dp/dr"); 
    if(isnan(ret) && !isnan(P)){
	std::cout << "NAN IN dP_dr" << std::endl;
    double delta_r = sqrt(3.0)/100.0; 
        std::cout << "halo s of r " << halo_S_of_r_Agora(r,Grid, binned_mass,MetaData) << " r " << r << std::endl;
	std::cout << "g" << halo_mod_g_of_r(r, binned_mass,Grid) << std::endl; 
	std::cout.flush();
        ENZO_FAIL("nan in dp/dr");
    }
    return ret;
}

/* More complex entropy profile from Voit 2019 that requires calculation of the cooling function.
   This one returns entropy in erg cm^2 instead of K cm^2 */
    double halo_S_of_r_Agora(double r, grid* Grid, FLOAT *binned_mass, TopGridData &MetaData){
    double M, C, r_vir, r_max, rho_crit = 1.8788e-29*0.49;
    double vcirc2, vcirc2_max;
    double Tgrav, Tgrav_therm;
    /* Get units */
    float DensityUnits=1, LengthUnits=1, VelocityUnits=1, TimeUnits=1,
      TemperatureUnits=1;
    double MassUnits=1;

    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		 &TimeUnits, &VelocityUnits, &MassUnits, Grid->ReturnTime()) == FAIL) {
      ENZO_FAIL("Error in GetUnits.");
    }
    double MassUnitsDouble = double(DensityUnits)*POW(double(LengthUnits), 3.0);
    M = binned_mass[99]*MassUnitsDouble;  // total mass in CGS
    C = 10;  // concentration parameter for NFW halo
    r_vir = POW(3.0/(4.0*3.14159)*M/(200.*rho_crit),1./3.);  // virial radius in CGS
    r_max = 2.163 * r_vir/C;
    
    vcirc2 = GravConst * MassEnclosed_r(binned_mass, r,Grid) / r;
    vcirc2_max = GravConst * MassEnclosed_r(binned_mass, r_max,Grid) / r_max;
    Tgrav = Mu*mh * vcirc2 / kboltz; // 2x gravitational "temperature"
    Tgrav_therm = Tgrav / TemperatureUnits / ((Gamma-1.0)*Mu); // code
  
    /* Calculate the cooling function Lambda using Grackle */
    double Lambda;
    double dens = mh/DensityUnits; // code
    double vx=0, vy=0, vz=0;
    double hi, hii, hei, heii, heiii, de, hm, h2i, h2ii, di, dii, hdi, metal; // species
    int dim=1;

    // setup_chem has densities in code, temperature in K
    setup_chem(dens, Tgrav, 1, de, hi, hii, hei, heii, heiii, hm, h2i, h2ii, di, dii, hdi);
    metal = 1e-6 * 0.02041 * dens;

    // temporarily disable UV background; makes S(r) trend downward at large r instead of upward
    // because of low Tgrav
    int saved_UVB = grackle_data->UVbackground;
    grackle_data->UVbackground = 0;
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
    double GasHaloRatio = 10;  
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
    double n_di, n_dii, n_hd; 
    if(MultiSpecies > 2){ 
        n_di = DensityUnits * di / m_d; 
        n_dii = DensityUnits * dii / m_d; 
        double m_hd = 5.018e-24; 
        n_hd = DensityUnits * hdi / m_hd; 
    }

    //double n_metal = DensityUnits * metal / (3*mh); 

    double n_i = n_hii + n_heii + n_heiii + n_h2ii + n_hm; 
    if(MultiSpecies > 2) 
        n_i += n_dii; 
    double n = n_hi + n_hii + n_hm + n_hei + n_heii + n_heiii + n_h2i + n_h2ii  + n_e; 
    if(MultiSpecies > 2) 
        n += n_di + n_dii + n_hd;
    /* Calculate entropy S(r) in erg cm^2 */
    double S_precip = POW(2*Mu*mh, 1./3.) * POW(r*Lambda*GasHaloRatio/3.0, 2./3.);
    //double S_precip = POW(2*mu*mh, 1./3.) * POW(20 * r * Lambda * n_i / (n * 3), 2./3.); 
    double S_nfw = 39. * vcirc2_max/1e10/4e4 * POW(r/r_vir, 1.1) / KEV_PER_ERG; // See Voit 19 Eqn 10 for assumptions
    // TODO blend with an entropy cap
    return (S_nfw + S_precip);
    
}

double MassEnclosed_r(FLOAT *binned_mass, double rad, grid* Grid){
	    /* Get units */
	    float DensityUnits=1, LengthUnits=1, VelocityUnits=1, TimeUnits=1,
	      TemperatureUnits=1;
	    double MassUnits=1;

	    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
			 &TimeUnits, &VelocityUnits, &MassUnits, Grid->ReturnTime()) == FAIL) {
	      ENZO_FAIL("Error in GetUnits.");
	    }
    	double MassUnitsDouble = double(DensityUnits)*POW(double(LengthUnits), 3.0);
	double delta_r = sqrt(3.0) / 100; 
	double prev_menc = 0.0;
	double next_menc = 0.0; 
	double ans; 
	rad /= LengthUnits; //to code
	int bin_index = rad / delta_r;
	double this_r_bin = delta_r * (bin_index + 0.5); 
	double prev_r_bin = delta_r * ((bin_index - 1) + 0.5); 
	double next_r_bin = delta_r * ((bin_index + 1) + 0.5); 
	if(prev_r_bin < 0){
		ans = binned_mass[bin_index] / delta_r * (rad); 
		return ans*MassUnitsDouble; 	
	}
	if(next_r_bin >= 100){
        if(debug)
            std::cout << "rbin too big!" << std::endl; 
		prev_menc = binned_mass[bin_index - 1];
	    	ans = (binned_mass[bin_index] - prev_menc) / delta_r * (rad - prev_r_bin) + prev_menc; 
		return ans*MassUnitsDouble; 	
	}
	prev_menc = binned_mass[bin_index - 1]; 
	next_menc = binned_mass[bin_index + 1]; 
	ans = 0.5 * (((binned_mass[bin_index] - prev_menc) / (delta_r) * (rad - prev_r_bin) + prev_menc) + ((next_menc - binned_mass[bin_index]) / delta_r * (next_r_bin - rad) + binned_mass[bin_index])); 
	return ans*MassUnitsDouble; //cgs 
}
double halo_mod_g_of_r(double r, FLOAT *binned_mass, grid* Grid){
  return GravConst*MassEnclosed_r(binned_mass, r, Grid)/(r*r);
}

float HaloGasDensity(FLOAT R, struct CGMdata& CGM_data, grid* Grid){
    /* assumes entropy is a power-law function of radius OR a cored power-law function
       of radius and gas is in hydrostatic equilibrium w/the NFW halo.  */

    double this_radius_cgs, Rstart;
    int index;
    /* Get units */
    float DensityUnits=1, LengthUnits=1, VelocityUnits=1, TimeUnits=1,
      TemperatureUnits=1;
    double MassUnits=1;

    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		 &TimeUnits, &VelocityUnits, &MassUnits, Grid->ReturnTime()) == FAIL) {
      ENZO_FAIL("Error in GetUnits.");
    }
    this_radius_cgs = R;  // radius in CGS
    index = int((this_radius_cgs-CGM_data.R_inner)/CGM_data.dr + 1.0e-3);  // index in array of CGM values
    if(index<0) index=0;  // check our indices
    //if(index>=CGM_data.nbins) index=CGM_data.nbins-1;
    if(index >= CGM_data.nbins) return 0.0; //outside of CGM
    return CGM_data.n_rad[index]*Mu*mh / DensityUnits;  // return physical density in code units
} // end HaloGasDensity

float HaloGasTemperature(FLOAT R, struct CGMdata& CGM_data, grid* Grid){
    /* assumes entropy is a power-law function of radius and gas is in hydrostatic equilibrium */

    /* Get units */
    float DensityUnits=1, LengthUnits=1, VelocityUnits=1, TimeUnits=1,
      TemperatureUnits=1;
    double MassUnits=1;

    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		 &TimeUnits, &VelocityUnits, &MassUnits, Grid->ReturnTime()) == FAIL) {
      ENZO_FAIL("Error in GetUnits.");
    }
    double this_radius_cgs;
    int index;
    this_radius_cgs = R; // radius in CGS
    index = int((this_radius_cgs-CGM_data.R_inner)/CGM_data.dr+1.0e-3);  // index in array of CGM values
    if(index<0) index=0;  // check our indices
    //if(index>=CGM_data.nbins) index=CGM_data.nbins-1;
    if(index >= CGM_data.nbins) return -1.0; //outside CGM, return negative temperature for easy flagging 
    return CGM_data.T_rad[index] / TemperatureUnits;  // return temperature in code units
  
}

