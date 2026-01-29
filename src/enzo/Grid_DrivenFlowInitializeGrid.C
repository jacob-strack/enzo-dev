/***********************************************************************
/
/  GRID CLASS: DrivenFlowInitializeGrid
/
/  written by: Wolfram Schmidt
/  date:       May, 2005
/  modified1:  Sep, 2014: modified to support enzo 2.4 // P. Grete
/
/  PURPOSE: Initializes grid for a driven flow simulation
/
/  RETURNS: FAIL or SUCCESS
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

int GetUnits(float *DensityUnits, float *LengthUnits,
              float *TemperatureUnits, float *TimeUnits,
              float *VelocityUnits, FLOAT *MassUnits, FLOAT Time);

int grid::DrivenFlowInitializeGrid(float DrivenFlowDensity,
                   float DrivenFlowPressure, float DrivenFlowMagField, int SetBaryonFields)
{
  /* create fields */
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum,
    H2IINum, DINum, DIINum, HDINum, B1Num, B2Num, B3Num, PhiNum;
  int CMNum, OMNum, CINum, OINum, OHINum, COINum, CHINum, CH2INum, C2INum, HCOINum, H2OINum, O2INum, CO_TOTALINum, H2O_TOTALINum, CIINum, OIINum, HOCIINum, HCOIINum, H3IINum, CHIINum, CH2IINum, COIINum, CH3IINum, OHIINum, H2OIINum, H3OIINum, O2IINum;  

  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;

  int vel = NumberOfBaryonFields;

  FieldType[NumberOfBaryonFields++] = Velocity1;
  FieldType[NumberOfBaryonFields++] = Velocity2;
  FieldType[NumberOfBaryonFields++] = Velocity3;

  if (EquationOfState == 0)
    FieldType[NumberOfBaryonFields++] = TotalEnergy;

  if (DualEnergyFormalism)
      FieldType[NumberOfBaryonFields++] = InternalEnergy;

  if (UseMHD) {
    FieldType[NumberOfBaryonFields++] = Bfield1;
    FieldType[NumberOfBaryonFields++] = Bfield2;
    FieldType[NumberOfBaryonFields++] = Bfield3;
  }
  if (HydroMethod == MHD_RK){
    FieldType[NumberOfBaryonFields++] = PhiField;
    
    if(UsePoissonDivergenceCleaning)
      FieldType[NumberOfBaryonFields++] = Phi_pField;
  }
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

  int accel = NumberOfBaryonFields;

  FieldType[NumberOfBaryonFields++] = DrivingField1;
  FieldType[NumberOfBaryonFields++] = DrivingField2;
  FieldType[NumberOfBaryonFields++] = DrivingField3;


  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (!SetBaryonFields)
    return SUCCESS;

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
        VelocityUnits;
  FLOAT MassUnits = 1;

  int size = 1;

  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];


  this->AllocateGrids();

  /* set density, total energy and velocity in problem dimension */

  float Energy = DrivenFlowPressure / ((Gamma-1.0) * DrivenFlowDensity);

  for (int i = 0; i < size; i++) {
    BaryonField[iden][i] = DrivenFlowDensity;
    BaryonField[ietot][i] = Energy;
    //new chem densities; init to 0 they get populated with time evolution for now 
    BaryonField[CMNum][i] = 0.0; 
    BaryonField[OMNum][i] = 0.0; 
    BaryonField[HCOINum][i] = 0.0; 
    BaryonField[CO_TOTALINum][i] = 0.0; 
    BaryonField[H2O_TOTALINum][i] = 0.0; 
    BaryonField[HOCIINum][i] = 0.0; 
    BaryonField[H3IINum][i] = 0.0; 
    BaryonField[CHIINum][i] = 0.0; 
    BaryonField[CH2IINum][i] = 0.0; 
    BaryonField[COIINum][i] = 0.0; 
    BaryonField[OHIINum][i] = 0.0; 
    BaryonField[H2OIINum][i] = 0.0; 
    BaryonField[H3OIINum][i] = 0.0; 
    BaryonField[O2IINum][i] = 0.0; 
  }

  if (DualEnergyFormalism) ///CF (internal energy = total energy, because v=0)
    for (int i = 0; i < size; i++) {
      BaryonField[ieint][i] = Energy;
    }
  printf("DrivenFlowInitializeGrid %"FSYM" %"FSYM" %"FSYM" %"FSYM"\n",
    DrivenFlowDensity,DrivenFlowPressure,Gamma,Energy);
  
  if (HydroMethod == MHD_RK) {
      for (int i = 0; i < size; i++) {
          BaryonField[iBx  ][i]  = DrivenFlowMagField;
      }
      Energy += 0.5 * pow(DrivenFlowMagField,2) / DrivenFlowDensity;
  }
  
  if ( UseMHDCT ){
    for ( int i = 0; i < MagneticSize[0]; i++){
       MagneticField[0][i] = DrivenFlowMagField;
    }
    Energy += 0.5 * pow(DrivenFlowMagField,2) / DrivenFlowDensity;
  this->CenterMagneticField();
  }
  
  if (EquationOfState == 0)
    for( int i = 0; i < size; i++)
      BaryonField[ietot][i] = Energy;

  return SUCCESS;
}
