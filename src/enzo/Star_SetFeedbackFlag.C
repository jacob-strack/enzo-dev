/***********************************************************************
/
/  FROM THE ATTRIBUTES, DETERMINE THE STAR TYPE AND FEEDBACK MODE
/
/  written by: John Wise
/  date:       November, 2005
/  modified1: Ji-hoon Kim
/             July, 2009
/
/  NOTES:  When the star particle is created, it is assigned the 
/          ParticleType from the grid.  Change this to represent the 
/          proper star_type (typedefs.h)
/
************************************************************************/
#include <stdlib.h>
#include <stdio.h>
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
#include "LevelHierarchy.h"
#include "phys_constants.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

void Star::SetFeedbackFlag(int flag)
{
  this->FeedbackFlag = flag;
  return;
}

#ifdef LARGE_INTS
void Star::SetFeedbackFlag(Eint32 flag)
{
  this->FeedbackFlag = flag;
  return;
}
#endif

int Star::SetFeedbackFlag(FLOAT Time)
{
  printf("In SetFeedbackFlag \n");
  const float TypeIILowerMass = 11, TypeIIUpperMass = 40.1;
  const float PISNLowerMass = 140, PISNUpperMass = 260;
  const float StarClusterSNeStart = 4.0;   // Myr after cluster is born
  const float StarClusterSNeEnd = 20.0; // Myr (lifetime of a 8 SolarMass star)

  int abs_type;
  float AgeInMyr;

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
    VelocityUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);

  abs_type = ABS(this->type);
  switch (abs_type) {

  case PopIII:
    if (this->type < 0) // birth
      this->FeedbackFlag = FORMATION;
    else if (Time > this->BirthTime + this->LifeTime) // endpoint
      if (((this->Mass >= PISNLowerMass && this->Mass <= PISNUpperMass) ||
	   (this->Mass >= TypeIILowerMass && this->Mass <= TypeIIUpperMass)) &&
	  PopIIISupernovaExplosions == TRUE)
	this->FeedbackFlag = SUPERNOVA;
      else
	this->FeedbackFlag = NO_FEEDBACK; // BH formation
    else // main sequence
      this->FeedbackFlag = NO_FEEDBACK;
    break;

  case SimpleSource:
    if (this->type < 0) // birth
      this->FeedbackFlag = FORMATION;
    
  case PopII:
    AgeInMyr = (Time - BirthTime) * TimeUnits / Myr_s;
    if (this->type > 0)
      if ((AgeInMyr > StarClusterSNeStart && AgeInMyr < StarClusterSNeEnd) ||
	  StarClusterUnresolvedModel)
	this->FeedbackFlag = CONT_SUPERNOVA;
      else
	this->FeedbackFlag = NO_FEEDBACK;
    else
      this->FeedbackFlag = FORMATION;
    break;

  case BlackHole:
    this->FeedbackFlag = NO_FEEDBACK;
    break;

  case PopIII_CF:
    if (this->type < 0) 
      this->FeedbackFlag = COLOR_FIELD;
    else 
      this->FeedbackFlag = NO_FEEDBACK;
    break;

  /* For MBH particle. Even with the NO_FEEDBACK flag, 
     the particle still can act as a Radiation Source if RadiativeTransfer = 1. */  
  case MBH:
    AgeInMyr = (Time - BirthTime) * TimeUnits / Myr_s;
    if (this->type > 0 && AgeInMyr > 0 && MBHFeedback > 0) {
      if (MBHFeedback == 1) 
	this->FeedbackFlag = MBH_THERMAL;
      if (MBHFeedback >= 2 && MBHFeedback <= 5) 
	this->FeedbackFlag = MBH_JETS;
    }
    else
      this->FeedbackFlag = NO_FEEDBACK; //It could still be a Radiation Source. 

#define NOT_SEDOV_TEST
#ifdef SEDOV_TEST
  //if (this->type > 0 && AgeInMyr > 0 && AgeInMyr < 0.001)  //for Sedov test (injecting for 1kyr)
    if (this->type > 0 && AgeInMyr > 0)                      //for Ostriker & McKee test (injecting continuously)
      this->FeedbackFlag = MBH_THERMAL;
    else
      this->FeedbackFlag = NO_FEEDBACK;      
#endif
    break;

  default:
    this->FeedbackFlag = NO_FEEDBACK;
    break;

   //this->type = abs_type;
 case PARTICLE_TYPE_STAR:
   if(AgoraICFeedback == FALSE && BirthTime < 0)
      this->FeedbackFlag = NO_FEEDBACK;
   printf("CORRECT PARTICLE TYPE \n");
   if(UseMagneticSupernovaFeedback)
     this->FeedbackFlag = SUPERNOVA_SEEDFIELD;
   break;

  }

  return SUCCESS;
}
