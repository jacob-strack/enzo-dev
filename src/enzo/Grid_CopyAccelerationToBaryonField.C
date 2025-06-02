
/***********************************************************************
/
/  GRID CLASS (COPY acceleration FIELD TO BARYON FIELD)
/
/  written by: Alexei Kritsuk
/  date:       Aug 2001
/  modified1:  Robert Harkness
/  date:       June 2004
/  modified2: dcollins, september 2024
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/
 
// Copy the potential field to baryon field
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
int FindField(int field, int farray[], int numfields);
 
int grid::CopyAccelerationToBaryonField()
{

  if (WriteAcceleration == FALSE)
    return SUCCESS;

  // Return if this doesn't concern us
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  // Find the field
 
  int fieldx = FindField(Acceleration0, FieldType, NumberOfBaryonFields);
  int fieldy = FindField(Acceleration1, FieldType, NumberOfBaryonFields);
  int fieldz = FindField(Acceleration2, FieldType, NumberOfBaryonFields);
 
  // Check to make sure BaryonField "GravPotential" exists
 
  if (BaryonField[fieldx] == NULL) {
    ENZO_FAIL("gx field missing.\n");
  }
  if (BaryonField[fieldy] == NULL) {
    ENZO_FAIL("gy field missing.\n");
  }
  if (BaryonField[fieldz] == NULL) {
    ENZO_FAIL("gz field missing.\n");
  }
 
  // Check to make sure PotentialField exists
 
  if ( AccelerationField == NULL) {
    ENZO_FAIL("AccelerationField missing.\n");
  }

  int size = GridDimension[0]*GridDimension[1]*GridDimension[2];

  for ( int dim=0;dim<GridRank;dim++)
      for(int index=0;index<size;index++){
          BaryonField[fieldx+dim][index] = AccelerationField[dim][index];
      }
  return SUCCESS;
}


