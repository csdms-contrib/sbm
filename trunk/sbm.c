#include <stdlib.h>		/* Sorted Bedforms Be Heir */
#include <stdio.h>
#include <math.h>		/* BTW, doubling Xmax and Ymax for periodic BC copy */
#include <time.h>		/* is now very unnecessary. Remove later */

#define INITIAL_DEPTH 21.25	/* meters !!!remember there is a 25*cell_height elevation to be subtracted to have the effective water depth */

#define VELOCITY_MEAN 0.1414	/* for use with GenRandomGaussian  x and y vels, m/s. 0.07071-> |v|=0.1 */
#define VELOCITY_SIGMA 0.0	/* for use with GenRandomGaussian */
/* WARNING! VELOCITY_MEAN NEEDS TO BE > 3*VELOCITY_SIGMA */

#define WAVEHEIGHT_MEAN 2.50	/*(meters) for use with GenRandomGaussian */
#define WAVEHEIGHT_SIGMA 0.0	/*        for use with GenRandomGaussian */
/* WARNING! WAVEHEIGHT_MEAN NEEDS TO BE > 3*WAVEHEIGHT_SIGMA */

#define AVE_PERCENT_COARSE 0.20	/* WARNING! The REAL average of % coarse is 1.5*AVE_PERCENT_COARSE */
#define	AggRateFine	0.0		/* aggradation rate of fine, in m/yr. It can be negative (erosion) and it will still work ok */
#define	AggRateCoarse	0.0		/* aggradation rate of coarse, in m/yr */
#define	SecondsPerYear	31536000.0		/* what it says */

#define CELL_WIDTH 5		/* m */
#define CELL_HEIGHT 0.05	/* m */

#define Xmax 100
#define Ymax 100
#define Zmax 100		/* this # should be related to depth
				   i.e., Zmax = DEPTH/CELLHEIGHT +
				   height of sediment accumulation */

#define AdjustTime 10		/* number of iterations sed trans across bdries */
				/* is given to adjust to change of direction */

#define PMAX_MEAN  0.5		/* for use with GenRandomGaussian */
#define PMAX_SIGMA 0.2		/* ditto */
#define NMAX_MEAN  0.5		/* ditto */
#define NMAX_SIGMA 0.2		/* ditto */



/* sediment transport constants */
/*#define Wf 0.0184     /* fall velocity of fine sediment */
/*#define Wc 0.0867     /* fall velocity of coarse sediment */
#define Wf 0.01			/* fall velocity of fine sediment */
#define Wc 0.10			/* fall velocity of coarse sediment */
#define Es 0.035		/* efficiency factor */
#define rho 1000		/* water density, kg/m^3 */
#define rhoS 2600		/* sediment density, kg/m^3 */
#define PI 3.1415927		/* PI */
#define g 9.8			/* gravity, m/s */
#define ALPH 0.01		/* O(0.01) constant, in eff. prof. ht. dep. on %coarse */
#define A 10			/* O(10) constant, in Cf dep. on %coarse */
#define ProfHtOverFine 0.001	/* 2*Uo*(1/w)*this -> effective prof. ht. when %coarse = 0 */
/* constants for new 7/03 sed trans XXX */
/*#define dfine 0.00022	    /*  diameter of fine sed, in m */
/*#define dcoarse 0.00075       /*  diameter of coarse sed, in m */
#define dfine 0.00015		/*  diameter of fine sed, in m */
#define dcoarse 0.001		/*  diameter of coarse sed, in m */
#define ConvertToEffGrSz  1.0	/* ConvertToEffGrSz times Se provides the effective d50 for determining ripple dims ... */
#define T 10.0			/* wave period, in s */
#define A2 10.0			/* O(10) const in ripple wavelength as func. of bed comp. */
#define A3 2.0			/* O(1) const in ripple aspect ratio as func. of bed comp. */
#define mu 0.001		/* water visc., kg/ms */
#define coeffslope  1.0		/* coefficient tuning sediment transport '2' */
#define nu 0.000001		/* kinematic visc., m2/s */
#define ezexpb -0.0		/* exponent to account for hiding in bedload */
#define ezexps 0.0		/* exponent to account for hiding in suspension */
#define karman 0.4		/* von karman constant */
#define hustar 2.0		/* von karman constant */


double timeStep = 100.0;
double flux4time = 1.0;
double ripfinemaxL = -99.99;
double ripfinemaxA = -99.99;
double ripfineminL = +99.99;
double ripfineminA = +99.99;
double ripcoarsemaxL = -99.99;
double ripcoarsemaxA = -99.99;
double ripcoarseminL = +99.99;
double ripcoarseminA = +99.99;

#define FORCING_DURATION 86400.0	/* 86400 s = 24hrs.; 
					   after FORCING_DURATION has passed, new
					   wave heights and current velocities are generated */
#define CURRENT_REVERSAL_PROB 2.0	/* probability that current will reverse itself
					   after one FORCING_DURATION */
#define POS_NEG_RATIO 2.0	/* current is POS_NEG_RATIO times as likely to switch
				   from neg to pos as from pos to neg */
double currentDirectionX = 1.0;
double currentDirectionY = 1.0;

double tmpBedslope;
double maxBedslope = 0.01;	/* for use in FindTimeStep(); */

double totalElapsedTime = 0.0;	/* in seconds */
double timeSinceForcingUpdate = 0.0;	/* in seconds */

/* integration weights (gaussian) */
double gaussw[5] =
  { 0.295524224714753, 0.269266719309996, 0.219086362515982,
0.149451349150581, 0.066671344308688 };
double gaussx[5] =
  { 0.148874338981631, 0.433395394129247, 0.679409568299024,
0.865063366688985, 0.973906528517172 };

/*file saving stuff*/
#define StartSavingAt 0.	/*units = FORCING_DURATIONSs */
#define FrameSpacing  1.0	/* units = FORCING_DURATIONs (presently 1 day); 
				   wave frequency ~= 864/day                    */
int NumFramesSaved = 0;		/* for keeping of track of when to save a new one */

//int maxRunTime = 30.;		/* measured in units of FORCING_DURATIONs; */
double maxRunTime = 30.;		/* measured in units of FORCING_DURATIONs; */
					  /* presently, 1 FORCING_DURATION = 1 day. */

struct cellStruct
/* this structure represents all charateristics of a cell in 2
   dimensions, X * Y.  the third dimension is represented by the
   arrays percentFull and percentCoarse, which are themselves arrays
   of Zmax size */
{
  double localFluxFineX, localFluxFineY, localFluxCoarseX, localFluxCoarseY;
  double excessFineSedOutX, excessFineSedOutY, excessCoarseSedOutX,
    excessCoarseSedOutY;
  double volumeIn;
  double volumeOut;
  double fineVolumeAdded;
  double coarseVolumeAdded;
  double percentFull[Zmax];
  double percentCoarse[Zmax];
  double volumeDeposited;
  double consistencyBelow;	/* this will represent the consistency */
  int activeZ;			/* active cell in the Z dimension, i.e. top of sed buildup */
  int hystcountfine;
  int hystcountcoarse;
  double RippleAspectRatioFineMix;
  double RippleLengthFineMix;
  double RippleHeightFineMix;
  double RippleAspectRatioCoarseMix;
  double RippleLengthCoarseMix;
  double RippleHeightCoarseMix;
  double depth;			/* indicates the depth of water */
  double sedimentHeight;	/*indicates sediment buildup */
};

struct cellStruct area[Xmax][Ymax];
/* this is the declaration of our main data structure.
   above components can be accessed by area.volumeIn, etc... i.e.
   name of struct period field;
   allocate double space for periodic boundary
   conditions */

/* int currentDirection;   1~X positive
			  2~X negative
			  3~Y positive
			  4~Y negative */

double excessOutofIterFineX[Ymax], localFluxOutofIterFineX[Ymax];
double excessOutofIterCoarseX[Ymax], localFluxOutofIterCoarseX[Ymax];
double excessOutofIterFineY[Xmax], localFluxOutofIterFineY[Xmax];
double excessOutofIterCoarseY[Xmax], localFluxOutofIterCoarseY[Xmax];

double currentVelocityX, currentVelocityY;
double waveHeight = WAVEHEIGHT_MEAN;
double waveHeightold = WAVEHEIGHT_MEAN;
int CELL_VOLUME = CELL_HEIGHT * CELL_WIDTH * CELL_WIDTH;


int DEBUG = 0, DEBUG2 = 0, DEBUG3 = 0, DEBUG4 = 0, DEBUG5 = 0;
	/* indicates debug mode, turn on (1) or off (0) printf statements */

/* forward declarations of functions */

/* helper functions */
double Raise (double b, double e);	/* implements ^ for doubles */
double GenRandomPercentage ();	/* generate a random coarseness */



/* initializers */
void InitConds ();		/* initialize conditions */
void InitCondsFromFile ();	/* initialize conditions from a file */
void PrintToFile ();		/* like it says */
FILE *SaveFile;
FILE *SaveForcing;
char savefilename[24];
char startFromFile = 'n';	/* user inputted */
char filename[24] = "none";	/* user inputted */
char readfilename[24] = "save.0004";

/*checkers*/
void FindAveBedHt ();

/* sed trans */

void DoIteration ();		/* function that makes things happen */
void DoIterationDummy ();	/* function that doesn't really make things happen */

void SedTransFine ();
void SedTransCoarse ();

void AdjustCells ();
void ZeroVars ();

double GenRandomGaussian (double, double);

double RunTimeClock ();
double ForcingClock ();
void UpdateRunTimeClock ();
void UpdateForcingClock ();
void WriteForcingConditions ();

void
initialize (void)
{
  /*  srand(500); */

  /*printf("shall we start from a file (y or n)? \n");
     scanf("%c", &startFromFile); */
  /*  if (startFromFile == 'y')
     {
     printf("filename to read? \n");
     scanf("%s", &readfilename);
     printf("and what is the number of iterations?");
     scanf("%d", &iterationCnt);     
     InitCondsFromFile();               
     } */
  if (startFromFile == 'y')
    InitCondsFromFile ();
  if (startFromFile == 'n')
    InitConds ();

  SaveForcing = fopen ("save.forcing", "w");

}

void
run (double until)
{
  int p, z;
  int icount = 0;

  if (until>0)
    maxRunTime = until;
fprintf (stderr, "%f\n", RunTimeClock ());
fprintf (stderr, "%f\n", maxRunTime * FORCING_DURATION);
  while ((RunTimeClock () < maxRunTime * FORCING_DURATION))
  {
/*   FORCING CONDITIONS ALONG A COASTLINE (CONSTRAINED CURRENT DIRECTION) */
     /**/ if (currentDirectionX > 0)
    {
      if (GenRandomPercentage () < CURRENT_REVERSAL_PROB)
      {
	currentDirectionX = -currentDirectionX;
	currentDirectionY = -currentDirectionY;
      }				/*  // if */
    }				/*  // if */
    else
    {
      if (GenRandomPercentage () < CURRENT_REVERSAL_PROB * POS_NEG_RATIO)
      {
	currentDirectionX = -currentDirectionX;
	currentDirectionY = -currentDirectionY;
      }				/*  // if */
    }				/*  // else */


/*    FORCING CONDITIONS AWAY FROM A COASTLINE (NO CONSTRAINED CURRENT DIRECTION */
/*
      if (GenRandomPercentage() > 0.5) currentDirectionX =  1.0;
      else                             currentDirectionX = -1.0;
      if (GenRandomPercentage() > 0.5) currentDirectionY =  1.0;
      else                             currentDirectionY = -1.0;
*/



/*    waveHeight = GenRandomGaussian(WAVEHEIGHT_MEAN, WAVEHEIGHT_SIGMA);
      currentVelocityX = -currentDirectionX * GenRandomGaussian(VELOCITY_MEAN, VELOCITY_SIGMA);
      currentVelocityY = -currentDirectionY * GenRandomGaussian(VELOCITY_MEAN, VELOCITY_SIGMA); */
    waveHeight = WAVEHEIGHT_MEAN;
    currentVelocityX = currentDirectionX * VELOCITY_MEAN;
    currentVelocityY = currentDirectionY * VELOCITY_MEAN;

    icount = icount + 1;
    if ((icount > 19) && (icount < 22))
    {
      waveHeight = 5.0;
      timeStep = 10.0;
    }
    else
    {
      timeStep = 100.0;
    }

    printf (" waveHeight= %f \n", waveHeight);
    printf (" currentVelocityX= %f \n", currentVelocityX);
    printf (" currentVelocityY= %f \n", currentVelocityY);


    for (p = 1; p <= AdjustTime; p++)
    {
      DoIterationDummy ();	/*Gives sed trans across boundaries time */
      /*to adjust to new direction */
    }				/* // for */

fprintf (stderr, "%f\n", ForcingClock ());
fprintf (stderr, "%f\n", FORCING_DURATION);
    while (ForcingClock () < FORCING_DURATION)
    {
      DoIteration ();
      z = area[50][50].activeZ;


      if (RunTimeClock () / FORCING_DURATION >
	  StartSavingAt + NumFramesSaved * FrameSpacing)
      {
	NumFramesSaved++;
	PrintToFile ();
	printf ("\nFrame %i saved\n", NumFramesSaved);
      }				/* // if */

      /* //   printf("."); */
      /*      printf("Time = %f hrs\n", RunTimeClock() / 3600.); */
      /* //   printf("%f %i \n", RunTimeClock()/FORCING_DURATION, NumFramesSaved); */
      FindAveBedHt ();
    }				/* // while */

    UpdateForcingClock ();

  }				/* // while */

}

void
finalize (void)
{
  printf ("Total elapsed time = %f hrs\n", RunTimeClock () / 3600.);
  fclose (SaveForcing);
}

int
new_main ()
{
  initialize ();
  run (maxRunTime);
  finalize ();

  return 0;
}

/* main looks at how many long you want to run, and then runs it */
old_main ()
{
  int p, z;
  int icount = 0;

  /*  srand(500); */

  /*printf("shall we start from a file (y or n)? \n");
     scanf("%c", &startFromFile); */
  /*  if (startFromFile == 'y')
     {
     printf("filename to read? \n");
     scanf("%s", &readfilename);
     printf("and what is the number of iterations?");
     scanf("%d", &iterationCnt);     
     InitCondsFromFile();               
     } */
  if (startFromFile == 'y')
    InitCondsFromFile ();
  if (startFromFile == 'n')
    InitConds ();

  SaveForcing = fopen ("save.forcing", "w");

/*    printf("Max run time = %f hrs", maxRunTime*FORCING_DURATION/3600.); */
  while ((RunTimeClock () < maxRunTime * FORCING_DURATION))
  {
fprintf (stderr ,"The time is %f\n", RunTimeClock ());
/*   FORCING CONDITIONS ALONG A COASTLINE (CONSTRAINED CURRENT DIRECTION) */
     /**/ if (currentDirectionX > 0)
    {
      if (GenRandomPercentage () < CURRENT_REVERSAL_PROB)
      {
	currentDirectionX = -currentDirectionX;
	currentDirectionY = -currentDirectionY;
      }				/*  // if */
    }				/*  // if */
    else
    {
      if (GenRandomPercentage () < CURRENT_REVERSAL_PROB * POS_NEG_RATIO)
      {
	currentDirectionX = -currentDirectionX;
	currentDirectionY = -currentDirectionY;
      }				/*  // if */
    }				/*  // else */


/*    FORCING CONDITIONS AWAY FROM A COASTLINE (NO CONSTRAINED CURRENT DIRECTION */
/*
      if (GenRandomPercentage() > 0.5) currentDirectionX =  1.0;
      else                             currentDirectionX = -1.0;
      if (GenRandomPercentage() > 0.5) currentDirectionY =  1.0;
      else                             currentDirectionY = -1.0;
*/



/*    waveHeight = GenRandomGaussian(WAVEHEIGHT_MEAN, WAVEHEIGHT_SIGMA);
      currentVelocityX = -currentDirectionX * GenRandomGaussian(VELOCITY_MEAN, VELOCITY_SIGMA);
      currentVelocityY = -currentDirectionY * GenRandomGaussian(VELOCITY_MEAN, VELOCITY_SIGMA); */
    waveHeight = WAVEHEIGHT_MEAN;
    currentVelocityX = currentDirectionX * VELOCITY_MEAN;
    currentVelocityY = currentDirectionY * VELOCITY_MEAN;

    icount = icount + 1;
    if ((icount > 19) && (icount < 22))
    {
      waveHeight = 5.0;
      timeStep = 10.0;
    }
    else
    {
      timeStep = 100.0;
    }

    printf (" waveHeight= %f \n", waveHeight);
    printf (" currentVelocityX= %f \n", currentVelocityX);
    printf (" currentVelocityY= %f \n", currentVelocityY);


    for (p = 1; p <= AdjustTime; p++)
    {
      DoIterationDummy ();	/*Gives sed trans across boundaries time */
      /*to adjust to new direction */
    }				/* // for */

    while (ForcingClock () < FORCING_DURATION)
    {
      DoIteration ();
      z = area[50][50].activeZ;


      if (RunTimeClock () / FORCING_DURATION >
	  StartSavingAt + NumFramesSaved * FrameSpacing)
      {
	NumFramesSaved++;
	PrintToFile ();
	printf ("\nFrame %i saved\n", NumFramesSaved);
      }				/* // if */

      /* //   printf("."); */
      /*      printf("Time = %f hrs\n", RunTimeClock() / 3600.); */
      /* //   printf("%f %i \n", RunTimeClock()/FORCING_DURATION, NumFramesSaved); */
      FindAveBedHt ();
    }				/* // while */

    UpdateForcingClock ();

  }				/* // while */

  printf ("Total elapsed time = %f hrs\n", RunTimeClock () / 3600.);
  fclose (SaveForcing);

}				/*  // main */




double
Raise (double b, double e)
/* this function acts like ^
   necessary because ^ does not allow doubles to be operands */
{
  if (b > 0)
    return pow (b, e);
  else
    return -pow (fabs (b), e);
}

double
GenRandomPercentage ()
/* this function generates a random number between 0 and 1 */
{
  int num = rand ();		/* rand returns a number from 1 to RAND_MAX */
  double realNum = (double) num / RAND_MAX;	/* smoosh it to b/w 0 and 1 */
  return realNum;
}

double
GenRandomGaussian (double mean, double sigma)	/* according to a Gaussian distribution */
{
  double rand, dist, step;
  double normalization;
  double totdist;
  double renorm;
  double nbins, excursion, signz;
  int i;

  nbins = 100;
  step = 3 * sigma / nbins;

  i = 0;
  totdist = 0.0;
  for (i = 0; i < nbins; i++)
  {
    totdist = totdist + exp (-(step * i) * (step * i) / 2 / sigma / sigma);
  }				/* // while normalization = normalization + step*exp(-(step*i - mean)*(step*i - mean)/2/sigma/sigma); */
  /* printf("ck1 = %f \n", totdist);       */

  rand = GenRandomPercentage ();
  i = 0;
  dist = 0.0;
  while (dist < rand)
  {
    dist =
      dist +
      (1 / totdist) * exp (-(step * i) * (step * i) / 2 / sigma / sigma);
    i++;
  }				/* // while */

  excursion = step * i;
  rand = GenRandomPercentage ();
  if (rand > 0.5)
    signz = 1;
  else
    signz = -1;

  /* printf("Time = %f \n", step*i);         */
  return mean + signz * excursion;
  /* //  return mean; */
}



void
InitConds ()
/* initializes area, %full, %coarse, buildup so far in each cell */
{
  int x, y, z;

  for (x = 0; x < Xmax; x++)
    for (y = 0; y < Ymax; y++)
      for (z = 0; z < Zmax; z++)
      {
	if ((z >= 0) && (z < 25))
	  /* near the floor, say its already full */
	{
	  area[x][y].percentCoarse[z] = AVE_PERCENT_COARSE +
	    AVE_PERCENT_COARSE * GenRandomPercentage ();
	  area[x][y].percentFull[z] = 1;
	}
	else if (z == 25)
	  /*  arbitrary, sometime in the middle it becomes uneven */
	{
	  area[x][y].percentCoarse[z] = AVE_PERCENT_COARSE +
	    AVE_PERCENT_COARSE * GenRandomPercentage ();
	  area[x][y].percentFull[z] = 0.5 /*GenRandomPercentage() */ ;
	  area[x][y].activeZ = 25;
	  area[x][y].hystcountfine = 0;
	  area[x][y].hystcountcoarse = 0;
	  area[x][y].RippleLengthFineMix = 0.;
	  area[x][y].RippleHeightFineMix = 0.;
	  area[x][y].RippleAspectRatioFineMix = 0.;
	  area[x][y].RippleLengthCoarseMix = 0.;
	  area[x][y].RippleHeightCoarseMix = 0.;
	  area[x][y].RippleAspectRatioCoarseMix = 0.;
	  /*if ( (x > Xmax/2 - 5 ) && 
	     (y > Ymax/2 - 5 )  )
	     {
	     area[x][y].percentCoarse[z] = 0.95;
	     area[x][y].percentFull[z] = 0.5;
	     } */
	}
	else
	  /* somewhere near the top it is empty, all water  */
	{
	  area[x][y].percentCoarse[z] = 0;
	  area[x][y].percentFull[z] = 0;
	}
	area[x][y].volumeIn = 0;
	area[x][y].volumeOut = 0;
	area[x][y].depth = INITIAL_DEPTH;
      }
  /* conceptual image:  cross section 
     |--------------------------------------------------------------|
     |                                                              |
     |       water, so %coarse = 0, % full = 0                      |
     |                                                              |
     |                                                              |
     |                                                              |
     |                                                              |
     |                                                              |
     z=25 |-----active Z, i.e., top layer of sediment--------------------|
     |  at z=25, random amount full, and random amount coarse       |
     |                                                              |
     |                                                              |
     |   all sediment, random amt of coarse, % full = 1             |
     |                                                              |
     |--------------------------------------------------------------|
   */


}


void
InitCondsFromFile ()
{
  int x, y, z;
  double tmp;


  if (DEBUG2)
    printf ("readfilename: %s \n", readfilename);
  SaveFile = fopen (readfilename, "r");
  for (x = 0; x < Xmax; x++)
    for (y = 0; y < Ymax; y++)
    {
      fscanf (SaveFile, "%d", &area[x][y].activeZ);
      /*      printf("x: %i, y: %i, area[x][y].activeZ: %d  \n ",x,y,area[x][y].activeZ);     */
    }


  for (x = 0; x < Xmax; x++)
    for (y = 0; y < Ymax; y++)
    {
      fscanf (SaveFile, "%lf", &area[x][y].percentFull[area[x][y].activeZ]);
      /*printf("x: %i, y: %i, area[x][y].percentFull[area[x][y].activeZ]: %f  \n ",x,y,area[x][y].percentFull[area[x][y].activeZ]); */
    }

  for (x = 0; x < Xmax; x++)
    for (y = 0; y < Ymax; y++)
    {
      fscanf (SaveFile, "%lf", &tmp);
      /*  EffectivePercentCoarse = (double)tmp;       */
    }


  for (x = 0; x < Xmax; x++)
    for (y = 0; y < Ymax; y++)
      for (z = 0; z <= area[x][y].activeZ; z++)
      {
	fscanf (SaveFile, "%lf", &area[x][y].percentCoarse[z]);
	/*      printf("x: %i, y: %i, area[x][y].percentCoarse[z]: %f  \n ",x,y,area[x][y].percentCoarse[z]); */
      }
  fclose (SaveFile);
/*		printf("ck3  \n ");	*/

  /*fill in above and below activeZ */
  for (x = 0; x < Xmax; x++)
    for (y = 0; y < (Ymax); y++)
    {
      for (z = area[x][y].activeZ + 1; z < Zmax; z++)
      {
	area[x][y].percentFull[z] = 0;
	area[x][y].percentCoarse[z] = 0;
      }
      for (z = 0; z < area[x][y].activeZ; z++)
	area[x][y].percentFull[z] = 1;
      area[x][y].volumeIn = 0;
      area[x][y].volumeOut = 0;
      area[x][y].depth = INITIAL_DEPTH;
    }

}


void
PrintToFile ()
{
  double EffectivePercentCoarse = -9.9;
  int x, y, z;
/*saving the usual files START*/
  sprintf (savefilename, "save.%04d", NumFramesSaved);
  SaveFile = fopen (savefilename, "w");
  if (!SaveFile)
  {
    printf ("problem opening output file\n");
    exit (1);
  }

  for (x = 0; x < Xmax; x++)
  {
    for (y = 0; y < Ymax; y++)
      fprintf (SaveFile, "%d ", area[x][y].activeZ);
    fprintf (SaveFile, "\n");
  }
  fprintf (SaveFile, "\n");
  for (x = 0; x < Xmax; x++)
  {
    for (y = 0; y < Ymax; y++)
      fprintf (SaveFile, "%12.5f",
	       area[x][y].percentFull[area[x][y].activeZ]);
    fprintf (SaveFile, "\n");
  }
  fprintf (SaveFile, "\n");
  for (x = 0; x < Xmax; x++)
  {
    for (y = 0; y < Ymax; y++)
    {
      z = area[x][y].activeZ;
      EffectivePercentCoarse = (area[x][y].percentFull[z] *
				area[x][y].percentCoarse[z]
				+ area[x][y].percentCoarse[z - 1]
				+ area[x][y].percentCoarse[z - 2]
				+ (1 - area[x][y].percentFull[z]) *
				area[x][y].percentCoarse[z - 3]) / 3;
      fprintf (SaveFile, "%12.5f", EffectivePercentCoarse);
    }
    fprintf (SaveFile, "\n");
  }
  fprintf (SaveFile, "\n");
  for (x = 0; x < Xmax; x++)
    for (y = 0; y < Ymax; y++)
    {
      for (z = 0; z <= area[x][y].activeZ; z++)
	fprintf (SaveFile, "%12.5f", area[x][y].percentCoarse[z]);
      fprintf (SaveFile, "\n");
    }
  fprintf (SaveFile, "\n");
  fclose (SaveFile);
/*saving the usual files END*/

/*saving the ripple geometry START*/
  sprintf (savefilename, "ripp.%04d", NumFramesSaved);
  SaveFile = fopen (savefilename, "w");

  fprintf (SaveFile, "%8.5f", ripfinemaxL);
  fprintf (SaveFile, "\n");
  fprintf (SaveFile, "%8.5f", ripfinemaxA);
  fprintf (SaveFile, "\n");
  fprintf (SaveFile, "%8.5f", ripfineminL);
  fprintf (SaveFile, "\n");
  fprintf (SaveFile, "%8.5f", ripfineminA);
  fprintf (SaveFile, "\n");
  fprintf (SaveFile, "%8.5f", ripcoarsemaxL);
  fprintf (SaveFile, "\n");
  fprintf (SaveFile, "%8.5f", ripcoarsemaxA);
  fprintf (SaveFile, "\n");
  fprintf (SaveFile, "%8.5f", ripcoarseminL);
  fprintf (SaveFile, "\n");
  fprintf (SaveFile, "%8.5f", ripcoarseminA);
  fprintf (SaveFile, "\n");
  fprintf (SaveFile, "\n");
  fclose (SaveFile);
/*saving the ripple geometry END*/

  WriteForcingConditions ();
}




void
FindAveBedHt ()
{
  int x, y, HtSum = 0;
  double HtAve;

  for (x = 0; x < Xmax; x++)
    for (y = 0; y < Ymax; y++)
      HtSum += area[x][y].activeZ +
	area[x][y].percentFull[area[x][y].activeZ];

  HtAve = (double) HtSum / (Xmax * Ymax);
/*  printf("Ave bed height: %f \n",HtAve); */

}


/* 4 sed. trans. proceedures:
   fine(positive, negative)
   coarse(positive, negative).
   Will need others for when can have pos x direction
   at same time as neg y, and v.v. (now both pos or both neg)
*/


void
SedTransFine ()
{
  double coeff = 16 * Es * rho / (3 * PI * Wf);
  double convertToVolume = CELL_WIDTH * timeStep / ((rhoS - rho) * g * 0.6);
  double effectiveProfileHeight = -99.99;
  double EffectivePercentCoarse = -99.99;
  double percentDeposited = -99.99;
  double Cf = -99.99;
  double localFlux = -99.99;
  double OrbitalVel = -99.99;
  double excessSed = -99.99;
  double excessSedInX = -99.99, excessSedInY = -99.99;
  double localFluxInX = -99.99, localFluxInY = -99.99;
  /* constants for new 7/03 set trans XXX */
  double ConvertToImmersedWt = (rhoS - rho) * g;	/* from volumetric flux, c*u */
  double OrbitalExcurs = -99.99;
  double d50, d50Ripples;
  double diffusivity = -99.99;
  double RippleLength = -99.99;
  double XX = -99.99;
  double RippleAspectRatiohyst = -99.99;
  double RippleHeightFineMixhyst = -99.99;
  double RippleLengthFineMixhyst = -99.99;
  double FrictionParam = -99.99;
  double ShieldsParam = -99.99;
  double ModifiedShieldsParam = -99.99;
  double BedConc = -99.99;
  double IntegratedConc = -99.99;
  double IntegratedConcx = -99.99;
  double IntegratedConcy = -99.99;
  double nikuk = -99.99;
  double fw1 = -99.99;
  double teta2aux = -99.99;
  double teta2 = -99.99;
  double km = -99.99;
  double FrictionVel = -99.99;
  double fw = -99.99;
  double zeta0 = -99.99;
  double zeta = -99.99;
  double zetamax = -99.99;
  double auxexp = -99.99;
  double auxslope = -99.99;
  double auxflufine1 = -99.99;
  double auxflufine2 = +99.99;
  double auxsort = -99.99;
  double auxripplel = -99.99;
  double auxripplemin = +99.99;
  double auxripplemax = -99.99;
  double auxflufine1EPC = -99.0;
  double auxflufine1BC = -99.0;
  double auxflufine1zeta = +99.0;
  double auxflufine1d50 = +99.0;
  double auxfinerippleL1 = +99.0;
  double auxfinerippleA1 = +99.0;
  double auxfineFV1 = -99.0;
  double auxflufine2EPC = +99.0;
  double auxflufine2BC = +99.0;
  double auxflufine2zeta = +99.0;
  double auxflufine2d50 = +99.0;
  double auxfinerippleL2 = +99.0;
  double auxfinerippleA2 = +99.0;
  double adfluxcoarse = -99.0;
  double slopefluxcoarse = -99.0;
  double adfluxfine = -99.0;
  double slopefluxfine = -99.0;
  double deltavr = -99.99;
  double kavr = -99.99;
  double VRvr = -99.99;
  double qbx = -99.99;
  double qby = -99.99;
  double qsx1 = -99.99;
  double qsy1 = -99.99;
  double qsx2 = -99.99;
  double qsy2 = -99.99;
  double Kdw = -99.99;
  double Kdw0 = -99.99;
  double Kaux = -99.99;
  double Se = 0.0;
  double Le = -99.99;
  double abl = -99.99;
  double deltabl = -99.99;
  double Ufcbl = -99.99;
  double kwobl = -99.99;
  double auxexp2 = -99.99;
  double stepbl = -99.99;
  double qbl = -99.99;
  double qobl = -99.99;
  double Ubl = -99.99;
  double Czbl = -99.99;
  double Uobl = -99.99;
  double Czobl = -99.99;
  double zetacenterofmass = -99.99;
  double CenterMassSum1 = -99.99;
  double CenterMassSum2 = -99.99;
  double RippleAspectRatio = -99.99;
  double Shieldscritical = -99.99;
  double ustar = -99.99;
  double vstar = -99.99;
  double Vbl = -99.99;
  double Vobl = -99.99;
  double qblu = -99.99;
  double qoblu = -99.99;
  double qblv = -99.99;
  double qoblv = -99.99;
  double CD = -99.99;
  double ShieldsRipple = -99.99;
  double coeffcrit = 1.0;
  double ssu = -99.99;
  double ssv = -99.99;
  double xrint = -99.99;
  double xmint = -99.99;
  double caux1 = -99.99;
  double caux2 = -99.99;
  double cauy1 = -99.99;
  double cauy2 = -99.99;
  double dxint = -99.99;
  double zeta1 = -99.99;


  /* convert from immersed-weight transport flux */
  int x;
  int y;
  int z;
  int jint;
  double uMagnitude, vMagnitude;
  double D;
  double bedslopeX = -99.99;
  double bedslopeY = -99.99;
  /* double bedslopeX,bedslopeY;  */

  int beginX, beginY, endX, endY;
  int incrementX, incrementY;

  uMagnitude = currentVelocityX;
  vMagnitude = currentVelocityY;

  if (uMagnitude > 0)
  {
    beginX = 0;
    endX = Xmax;
    incrementX = 1;
  }				/* // if */
  else
  {
    beginX = Xmax - 1;
    endX = -1;
    incrementX = -1;
    uMagnitude = -uMagnitude;
  }				/* // else */

  if (vMagnitude > 0)
  {
    beginY = 0;
    endY = Ymax;
    incrementY = 1;
  }				/* // if */
  else
  {
    beginY = Ymax - 1;
    endY = -1;
    incrementY = -1;
    vMagnitude = -vMagnitude;
  }				/* // else */

  for (x = beginX; x != endX; x = x + incrementX)
    for (y = beginY; y != endY; y = y + incrementY)
    {

      z = area[x][y].activeZ;
      /*D = area[x][y].depth; */
      D = area[x][y].depth - (z + area[x][y].percentFull[z]) * CELL_HEIGHT;
      coeffcrit = 1.0;

      if (x == beginX)
      {
	 /**/ excessSedInX = excessOutofIterFineX[y];
	localFluxInX = localFluxOutofIterFineX[y];	/*these give time delay */
      }
      else
      {
	 /**/ excessSedInX = area[x - incrementX][y].excessFineSedOutX;
	localFluxInX = area[x - incrementX][y].localFluxFineX;
      }
      if (y == beginY)
      {
	 /**/ excessSedInY = excessOutofIterFineY[x];
	localFluxInY = localFluxOutofIterFineY[x];
      }
      else
      {
	 /**/ excessSedInY = area[x][y - incrementY].excessFineSedOutY;
	localFluxInY = area[x][y - incrementY].localFluxFineY;
      }
       /**/ if (x == endX - incrementX)
      {
	bedslopeX = (area[beginX][y].activeZ
		     + area[beginX][y].percentFull[area[beginX][y].activeZ]
		     - (z +
			area[x][y].percentFull[z])) * (CELL_HEIGHT /
						       CELL_WIDTH);
      }
      else
      {
	bedslopeX = (area[x + incrementX][y].activeZ
		     + area[x +
			    incrementX][y].percentFull[area[x +
							    incrementX][y].
						       activeZ] - (z +
								   area[x][y].
								   percentFull
								   [z])) *
	  (CELL_HEIGHT / CELL_WIDTH);
      }
      if (y == endY - incrementY)
      {
	bedslopeY = (area[x][beginY].activeZ
		     + area[x][beginY].percentFull[area[x][beginY].activeZ]
		     - (z +
			area[x][y].percentFull[z])) * (CELL_HEIGHT /
						       CELL_WIDTH);
      }
       /**/
      else
      {
	bedslopeY = (area[x][y + incrementY].activeZ
		     + area[x][y +
			       incrementY].percentFull[area[x][y +
							       incrementY].
						       activeZ] - (z +
								   area[x][y].
								   percentFull
								   [z])) *
	  (CELL_HEIGHT / CELL_WIDTH);
      }

      tmpBedslope = sqrt (bedslopeX * bedslopeX + bedslopeY * bedslopeY);
      if (tmpBedslope > maxBedslope)
      {
	maxBedslope = tmpBedslope;
      }

      /*  Dispersion relation is solved using Dean and Dalrymple p. 72. Starting value for the loop is Kdw0 */
      Kdw0 =
	Raise (2 * PI / T,
	       2) * (1 / g) * (Raise (tanh (Raise (2 * PI / T, 2) * (D / g)),
				      -0.5));
      Kdw = Raise (2 * PI / T, 2) * (1 / g) * (1 / tanh (Kdw0 * D));
      Kaux = fabs (1 - Kdw0 / Kdw);
      while (Kaux >= 0.01)
      {
	Kdw0 = Kdw;
	Kdw = Raise (2 * PI / T, 2) * (1 / g) * (1 / tanh (Kdw0 * D));
	Kaux = fabs (1 - Kdw0 / Kdw);
      }
      OrbitalVel = (PI * waveHeight / T) * (1 / sinh (Kdw * D));


      /* ave %coarse over top 3*CELL_HEIGHT in a way that changes smoothly */
      EffectivePercentCoarse =
	(area[x][y].percentFull[z] * area[x][y].percentCoarse[z] +
	 area[x][y].percentCoarse[z - 1] + area[x][y].percentCoarse[z - 2] +
	 (1 - area[x][y].percentFull[z]) * area[x][y].percentCoarse[z -
								    3]) / 3;


      OrbitalExcurs = 0.5 * waveHeight / sinh (Kdw * D);
      d50 =
	dfine * (1 - EffectivePercentCoarse) +
	dcoarse * (EffectivePercentCoarse);
      /*Se=(1-EffectivePercentCoarse)*pow(dfine-d50,2) + EffectivePercentCoarse*pow(dcoarse-d50,2);
         Se=Raise(Se,0.5);
         d50Ripples = ConvertToEffGrSz*d50;
         d50Ripples = d50+ConvertToEffGrSz*Se; */
      d50Ripples =
	dfine * (1 - EffectivePercentCoarse) +
	ConvertToEffGrSz * dcoarse * (EffectivePercentCoarse);
      XX =
	4 * nu * Raise ((OrbitalExcurs * 2 * PI / (T * 1.0)),
			2) / (d50Ripples * Raise (1.65 * g * d50Ripples,
						  1.5));
      if (XX <= 2)
      {
	RippleAspectRatiohyst = 0.15 / Raise (XX, 0.11);
	RippleLengthFineMixhyst = 1.96 * OrbitalExcurs / Raise (XX, 0.28);
	RippleHeightFineMixhyst =
	  RippleAspectRatiohyst * RippleLengthFineMixhyst;
      }
      else
      {
	RippleAspectRatiohyst = 0.166 / Raise (XX, 0.24);
	RippleLengthFineMixhyst = 2.71 * OrbitalExcurs / Raise (XX, 0.75);
	RippleHeightFineMixhyst =
	  RippleAspectRatiohyst * RippleLengthFineMixhyst;
      }

      /* Bed concentration    
         Shieldscritical=Raise(dfine/d50,ezexps)*0.04; */
      Shieldscritical = 0.04;
      nikuk = 2.5 * d50;
      FrictionParam =
	exp (5.213 * Raise (nikuk / OrbitalExcurs, 0.194) - 5.977);
/*ATTENTION! d50 is used to evaluate ShieldsParam. The assumption is that the whole bed composition is key to putting sediment in suspension,
the diffusivity profile is instead related to the individual grain size under consideration
	ShieldsParam = (FrictionParam/(3.3*g*dfine) )*Raise(OrbitalExcurs*2*PI/T,2);	*/
      ShieldsParam =
	(FrictionParam / (3.3 * g * d50)) * Raise (OrbitalExcurs * 2 * PI / T,
						   2);
      /* Ripple hysterisis. Ripple dimensions are updated only if the shields parameter, 0.04, is exceeded    */
      ShieldsRipple =
	(FrictionParam / (3.3 * g * d50Ripples)) * Raise (OrbitalExcurs * 2 *
							  PI / T, 2);
      if (ShieldsRipple > Shieldscritical)
      {
	area[x][y].RippleAspectRatioFineMix = RippleAspectRatiohyst;
	area[x][y].RippleLengthFineMix = RippleLengthFineMixhyst;
	area[x][y].RippleHeightFineMix = RippleHeightFineMixhyst;
      }
      /*else */

      ModifiedShieldsParam =
	ShieldsParam / pow (1 - PI * area[x][y].RippleAspectRatioFineMix, 2);
      BedConc = 0.005 * Raise (ModifiedShieldsParam - Shieldscritical, 3);

/* velocity profile */
/* inside the boundary layer	Fredsoe and Deigaard, p. 59-60, deltabl from eq.2.45(p.25)*/
/* we are adapting eq.2.45 to account for the presence of ripples (use of km)	*/
/* concentration profile (Nielsen p. 258; Fredsoe and Deigaard p.305 and eq.10.24 for diffusivity )*/
/* we are assuming the current doesn't affect concentration profile and only advects sediment*/
/* we are neglecting the flux contribution to sediment fluxes below the ripple crest because:
1) we have noidea of the flow in there
2) we assume that contribution is included in the bedload prediction
	abl = OrbitalVel * T /(2*PI);
	deltabl = km * 0.09 * Raise(abl/km,0.82);
	kwobl = 30*deltabl/(exp(0.4*VRvr/Ufcbl));
	
	VRvr=Raise(uMagnitude*uMagnitude+vMagnitude*vMagnitude,0.5);
	Ufcbl = Raise(2/PI*fw*OrbitalVel*VRvr,0.5);*/

      auxexp2 =
	area[x][y].RippleLengthFineMix * OrbitalExcurs * 2 * PI / T *
	exp (1.5 - 4500 * dfine -
	     1.2 * log (OrbitalExcurs * 2 * PI / T / Wf));

      zeta0 =
	(28 * area[x][y].RippleHeightFineMix *
	 area[x][y].RippleAspectRatioFineMix + 2.5 * d50) / 30;
      Czbl = 0.0;
      Czobl = 0.0;
      ustar = 0.4 * uMagnitude / log (hustar / (zeta0));
      vstar = 0.4 * vMagnitude / log (hustar / (zeta0));
      qoblu = 0.0;
      qoblv = 0.0;
      qblu = 0.0;
      qblv = 0.0;
      stepbl = 0.025;

      if (BedConc > 0.0)
      {

/*		
 outside the boundary layer	Fredsoe and Deigaard, p. 60-61*/
	Czobl = 1.0;
	zeta = zeta0;
	zeta = area[x][y].RippleHeightFineMix;
	zeta1 = 1.0 * area[x][y].RippleHeightFineMix;
	Uobl = uMagnitude;
	Vobl = vMagnitude;
	while (Czobl > 0.000003)
	{
	  Czobl = BedConc * exp (-Wf * (zeta - zeta1) / auxexp2);
	  if (zeta < hustar)
	  {
	    Uobl = 1 / 0.4 * ustar * log (zeta / (zeta0));
	    Vobl = 1 / 0.4 * vstar * log (zeta / (zeta0));
	  }
	  else
	  {
	    Uobl = uMagnitude;
	    Vobl = vMagnitude;
	  }
	  qoblu = qoblu + Czobl * Uobl * stepbl;
	  qoblv = qoblv + Czobl * Vobl * stepbl;
	  zeta = zeta + stepbl;
	}

	IntegratedConcx = qoblv;
	IntegratedConcy = qoblv;
      }
      else
      {
	IntegratedConcx = 0.0;
	IntegratedConcy = 0.0;
	coeffcrit = 0.0;
      }
/*bedload (Van Rijn Principles of coastal morphology p.4.38) */
      qbx =
	coeffcrit * 9.1 * Raise (ModifiedShieldsParam - Shieldscritical,
				 1.78) * Raise (g * 1.65 * d50 * d50 * d50,
						0.5) * uMagnitude /
	fabs (uMagnitude);
      qby =
	coeffcrit * 9.1 * Raise (ModifiedShieldsParam - Shieldscritical,
				 1.78) * Raise (g * 1.65 * d50 * d50 * d50,
						0.5) * vMagnitude /
	fabs (vMagnitude);

/*	Cf = 0.01 * (1 + A * EffectivePercentCoarse); */

      CD = Raise (karman / log (hustar / (zeta0)), 2);
      qsx1 = (ConvertToImmersedWt * IntegratedConcx);
      qsx2 =
	coeffcrit * coeffslope * coeff * CD * (1 / (5 * Wf)) *
	Raise (OrbitalVel, 5) * bedslopeX;

      qsy1 = (ConvertToImmersedWt * IntegratedConcy);
      qsy2 =
	coeffcrit * coeffslope * coeff * CD * (1 / (5 * Wf)) *
	Raise (OrbitalVel, 5) * bedslopeY;

      area[x][y].localFluxFineX =
	(1 - EffectivePercentCoarse) * (qbx + qsx1 - qsx2);
      area[x][y].localFluxFineY =
	(1 - EffectivePercentCoarse) * (qby + qsy1 - qsy2);
      
/*	if (x==50 && y==50)
	{
		printf("qbx: %f qsx1: %f qsx2: %f \n ",qbx, qsx1, qsx2);
		printf("orbital vel: %f excursion: %f  \n ",OrbitalVel, OrbitalExcurs);
		printf("bedconc: %f coeffcrit: %f  int conc: %f \n ",BedConc, coeffcrit, IntegratedConcx);
		printf("ModifiedShieldsParam: %f Shieldscritical: %f   \n \n",ModifiedShieldsParam,Shieldscritical);
		
	}	*/
	if (EffectivePercentCoarse > auxflufine1)
      {
	auxflufine1EPC = area[x][y].activeZ;
	auxflufine1 = EffectivePercentCoarse;
	auxflufine1d50 = d50Ripples;
	auxflufine1BC = XX;
	auxflufine1zeta = d50;
	ripfinemaxL = area[x][y].RippleLengthFineMix;
	ripfinemaxA = area[x][y].RippleHeightFineMix;
	adfluxcoarse = qsx1;
	slopefluxcoarse = qsx2;
      }

      if (EffectivePercentCoarse < auxflufine2)
      {
	auxflufine2EPC = area[x][y].activeZ;
	auxflufine2 = EffectivePercentCoarse;
	auxflufine2d50 = d50Ripples;
	auxflufine2BC = XX;
	auxflufine2zeta = d50;
	ripfineminL = area[x][y].RippleLengthFineMix;
	ripfineminA = area[x][y].RippleHeightFineMix;
	adfluxfine = qsx1;
	slopefluxfine = qsx2;
      }
      /**/
	/*effectiveProfileHeight = 0.001*D/Wf;
	   effectiveProfileHeight= (2*OrbitalVel/Wf)*( ProfHtOverFine + ALPH * 1.00 );
	   effectiveProfileHeight = 0.35;
	   effectiveProfileHeight = (2*OrbitalVel/Wf)*( ProfHtOverFine +
	   ALPH * EffectivePercentCoarse ); */
	effectiveProfileHeight = auxexp2 / Wf;

      if (uMagnitude > vMagnitude)
	percentDeposited = Wf * CELL_WIDTH / (uMagnitude *
					      effectiveProfileHeight);
      if (uMagnitude <= vMagnitude)
	percentDeposited = Wf * CELL_WIDTH / (vMagnitude *
					      effectiveProfileHeight);

      if (percentDeposited < 0)
      {
	printf (" \n Woo! This can not be! (percent deposited < 0) \n");
	percentDeposited = 0.0000000001;
	printf (" xfine= %f \n", x);
	printf (" yfine= %f \n", y);
	printf (" %g\n %g\n %g\n %g \n %g \n",
		area[x][y].percentCoarse[z],
		area[x][y].percentCoarse[z - 1],
		area[x][y].percentCoarse[z - 2],
		area[x][y].percentCoarse[z - 3], area[x][y].percentFull[z]);
	printf (" um= %f \n", uMagnitude);
	printf (" vm= %f \n", vMagnitude);
	printf (" eph= %f \n", effectiveProfileHeight);
	printf (" ProfHtOverFine = %f \n", ProfHtOverFine);
	printf (" EffectivePercentCoarse= %f \n", EffectivePercentCoarse);
      }
      if (percentDeposited > 1)
      {
	percentDeposited = 1;
      }

      excessSed = (localFluxInX - area[x][y].localFluxFineX)
	+ (localFluxInY - area[x][y].localFluxFineY)
	+ (excessSedInX + excessSedInY);
      /*think about negative excessOut... */
      area[x][y].excessFineSedOutX =
	(1 - percentDeposited) * excessSed
	* (uMagnitude / (uMagnitude + vMagnitude));
      area[x][y].excessFineSedOutY =
	(1 - percentDeposited) * excessSed
	* (vMagnitude / (uMagnitude + vMagnitude));
      area[x][y].fineVolumeAdded =
	convertToVolume * percentDeposited * excessSed;


      if (DEBUG5 && area[x][y].volumeOut < 0)
      {
	printf ("\nBeep: vol. out < 0 ; %f\n  at %d,%d, localflux %f, \n \n ",
		area[x][y].volumeOut, x, y, localFlux);
	printf ("percentcoarse[z]=%f\n", area[x][y].percentCoarse[z]);
	printf ("percentcoarse[z-1]=%f\n", area[x][y].percentCoarse[z - 1]);
	printf ("percentcoarse[z-2]=%f\n", area[x][y].percentCoarse[z - 2]);
      }

      /*record how much susp. load to pass in the other side next iteration */
      if (x == endX - incrementX)
      {
	excessOutofIterFineX[y] = area[x][y].excessFineSedOutX;
	localFluxOutofIterFineX[y] = area[x][y].localFluxFineX;
      }
      if (y == endY - incrementY)
      {
	excessOutofIterFineY[x] = area[x][y].excessFineSedOutY;
	localFluxOutofIterFineY[x] = area[x][y].localFluxFineY;
      }

/*if (x==3 && y==78)
{ 
  printf(" ModifiedShieldsParam: %f ShieldsParam: %f \n",ModifiedShieldsParam,ShieldsParam);
  printf(" BedConc: %f diffusivity: %f \n",BedConc,auxexp2);
  printf(" d50: %f d50Ripples: %f \n",d50,d50Ripples);
  printf(" qsx1: %f qsx2: %f \n",qsx1,qsx2);
  printf(" qsy1: %f qbx: %f \n",qsy1,qbx);
  printf(" RippleL: %f RippleRatio: %f \n",area[x][y].RippleLengthFineMix,area[x][y].RippleAspectRatioFineMix);
  printf(" orb vel: %f EffPercentCoarse: %f \n",OrbitalVel,EffectivePercentCoarse);
  printf("x,y: %d %d vMagnitude: %g  effProfHt: %g \n",x,y, vMagnitude, effectiveProfileHeight);
  
}	
  printf("   convertToVol: %g percent deposited: %g excessSed %g\n",convertToVolume,
						percentDeposited,excessSed);
  printf("   localFluxInX: %g localFluxFineX(x,y): %g \n", localFluxInX, area[x][y].localFluxFineX);
  printf("   localFluxInY: %g localFluxFineY(x,y): %g \n", localFluxInY, area[x][y].localFluxFineY);
  printf("   excessSedInX: %g excessSedInY: %g \n",excessSedInX, excessSedInY);
  printf("   fineVolAdded: %g \n", area[x][y].fineVolumeAdded);*/
    }
  
    /*      printf("percent coarse in coarsest = %f \n",auxflufine1);
       printf("active z in coarsest = %f \n",auxflufine1EPC);
       printf("d50ripples coarsest = %f \n",auxflufine1d50);
       printf("XX in coarsest = %f \n",auxflufine1BC);
       printf("ripple L in coarsest = %f \n",ripfinemaxL);
       printf("ripple Ratio in coarsest = %f \n \n",ripfinemaxA/ripfinemaxL);
       printf("percent coarse in finest = %f \n",auxflufine2);
       printf("active z in finest = %f \n",auxflufine2EPC);
       printf("d50ripples finest = %f \n",auxflufine2d50);
       printf("XX in finest = %f \n",auxflufine2BC);
       printf("ripple L in finest = %f \n",ripfineminL);
       printf("ripple Ratio in finest = %f \n \n",ripfineminA/ripfineminL); */
}


void
SedTransCoarse ()
{
  double coeff = 16 * Es * rho / (3 * PI * Wc);
  double convertToVolume = CELL_WIDTH * timeStep / ((rhoS - rho) * g * 0.6);
  double effectiveProfileHeight = -99.99;
  double EffectivePercentCoarse = -99.99;
  double percentDeposited = -99.99;
  double Cf = -99.99;
  double localFlux = -99.99;
  double OrbitalVel = -99.99;
  double excessSed = -99.99;
  double excessSedInX = -99.99, excessSedInY = -99.99;
  double localFluxInX = -99.99, localFluxInY = -99.99;
  double ConvertToImmersedWt = (rhoS - rho) * g;	/* from volumetric flux, c*u */
  double OrbitalExcurs = -99.99;
  double d50, d50Ripples;
  double diffusivity = -99.99;
  double RippleLength = -99.99;
  double XX = -99.99;
  double BedConc = -99.99;
  double RippleAspectRatiohyst = -99.99;
  double RippleHeightCoarseMixhyst = -99.99;
  double RippleLengthCoarseMixhyst = -99.99;
  double FrictionParam = -99.99;
  double ShieldsParam = -99.99;
  double ModifiedShieldsParam = -99.99;
  double IntegratedConc = -99.99;
  double IntegratedConcx = -99.99;
  double IntegratedConcy = -99.99;
  double nikuk = -99.99;
  double fw1 = -99.99;
  double teta2aux = -99.99;
  double teta2 = -99.99;
  double km = -99.99;
  double FrictionVel = -99.99;
  double fw = -99.99;
  double zeta0 = -99.99;
  double zeta = -99.99;
  double zetamax = -99.99;
  double auxexp = -99.99;
  double auxflufine1 = -99.99;
  double auxflufine2 = +99.99;
  double auxsort = -99.99;
  double auxripplel = -99.99;
  double auxripplemin = +99.99;
  double auxripplemax = -99.99;
  double auxflucoarse1 = -99.0;
  double auxflucoarse1EPC = -99.0;
  double auxflucoarse1BC = -99.0;
  double auxflucoarse1zeta = -99.0;
  double auxflucoarse1zeta0 = -99.0;
  double auxflucoarse1d50 = -99.0;
  double auxcoarserippleL1 = -99.0;
  double auxcoarserippleA1 = -99.0;
  double auxflucoarseIX1 = -99.0;
  double auxcoarseFV1 = -99.0;
  double auxflucoarse2 = +99.0;
  double auxflucoarse2EPC = +99.0;
  double auxflucoarse2BC = +99.0;
  double auxflucoarse2zeta = +99.0;
  double auxflucoarse2zeta0 = +99.0;
  double auxflucoarse2d50 = +99.0;
  double auxcoarserippleL2 = +99.0;
  double auxcoarserippleA2 = +99.0;
  double auxflucoarseIX2 = -99.0;
  double deltavr = -99.99;
  double ksvr = -99.99;
  double kavr = -99.99;
  double VRvr = -99.99;
  double qbx = -99.99;
  double qby = -99.99;
  double qsx1 = -99.99;
  double qsy1 = -99.99;
  double qsx2 = -99.99;
  double qsy2 = -99.99;
  double Kdw = -99.99;
  double Kdw0 = -99.99;
  double Kaux = -99.99;
  double Se = -99.99;
  double Le = -99.99;
  double abl = -99.99;
  double deltabl = -99.99;
  double Ufcbl = -99.99;
  double kwobl = -99.99;
  double auxexp2 = -99.99;
  double stepbl = -99.99;
  double qbl = -99.99;
  double qobl = -99.99;
  double Ubl = -99.99;
  double Czbl = -99.99;
  double Uobl = -99.99;
  double Czobl = -99.99;
  double zetacenterofmass = -99.99;
  double CenterMassSum1 = -99.99;
  double CenterMassSum2 = -99.99;
  double RippleAspectRatio = -99.99;
  double Shieldscritical = -99.99;
  double ustar = -99.99;
  double vstar = -99.99;
  double Vbl = -99.99;
  double Vobl = -99.99;
  double qblu = -99.99;
  double qoblu = -99.99;
  double qblv = -99.99;
  double qoblv = -99.99;
  double CD = -99.99;
  double ShieldsRipple = -99.99;
  double coeffcrit = 1.0;
  double ssu = -99.99;
  double ssv = -99.99;
  double xrint = -99.99;
  double xmint = -99.99;
  double caux1 = -99.99;
  double caux2 = -99.99;
  double cauy1 = -99.99;
  double cauy2 = -99.99;
  double dxint = -99.99;
  double zeta1 = -99.99;

  /* convert from immersed-weight transport flux */
  int x;
  int y;
  int z;
  int jint;
  double uMagnitude, vMagnitude;
  double D;
  double bedslopeX, bedslopeY;

  int beginX, beginY, endX, endY;
  int incrementX, incrementY;

  uMagnitude = currentVelocityX;
  vMagnitude = currentVelocityY;

  if (uMagnitude > 0)
  {
    beginX = 0;
    endX = Xmax;
    incrementX = 1;
  }				/* // if */
  else
  {
    beginX = Xmax - 1;
    endX = -1;
    incrementX = -1;
    uMagnitude = -uMagnitude;
  }				/* // else */

  if (vMagnitude > 0)
  {
    beginY = 0;
    endY = Ymax;
    incrementY = 1;
  }				/* // if */
  else
  {
    beginY = Ymax - 1;
    endY = -1;
    incrementY = -1;
    vMagnitude = -vMagnitude;
  }				/* // else */

  for (x = beginX; x != endX; x = x + incrementX)
    for (y = beginY; y != endY; y = y + incrementY)
    {

      z = area[x][y].activeZ;
      D = area[x][y].depth - (z + area[x][y].percentFull[z]) * CELL_HEIGHT;
      coeffcrit = 1.0;

      if (x == beginX)
      {
	 /**/ excessSedInX = excessOutofIterCoarseX[y];
/*MUCH TO DO HERE: DEFINE THESE VARS, ZERO THEM, MAKE SURE THEY'RE SET...*/
	localFluxInX = localFluxOutofIterCoarseX[y];	/*these give time delay */
      }
      else
      {
	 /**/ excessSedInX = area[x - incrementX][y].excessCoarseSedOutX;
	localFluxInX = area[x - incrementX][y].localFluxCoarseX;
      }
      if (y == beginY)
      {
	 /**/ excessSedInY = excessOutofIterCoarseY[x];
	localFluxInY = localFluxOutofIterCoarseY[x];
      }
      else
      {
	 /**/ excessSedInY = area[x][y - incrementY].excessCoarseSedOutY;
	localFluxInY = area[x][y - incrementY].localFluxCoarseY;
      }


       /**/ if (x == endX - incrementX)
      {
	bedslopeX = (area[beginX][y].activeZ
		     + area[beginX][y].percentFull[area[beginX][y].activeZ]
		     - (z +
			area[x][y].percentFull[z])) * (CELL_HEIGHT /
						       CELL_WIDTH);
      }
      else
      {
	bedslopeX = (area[x + incrementX][y].activeZ
		     + area[x +
			    incrementX][y].percentFull[area[x +
							    incrementX][y].
						       activeZ] - (z +
								   area[x][y].
								   percentFull
								   [z])) *
	  (CELL_HEIGHT / CELL_WIDTH);
      }
      if (y == endY - incrementY)
      {
	bedslopeY = (area[x][beginY].activeZ
		     + area[x][beginY].percentFull[area[x][beginY].activeZ]
		     - (z +
			area[x][y].percentFull[z])) * (CELL_HEIGHT /
						       CELL_WIDTH);
      }
       /**/
      else
      {
	bedslopeY = (area[x][y + incrementY].activeZ
		     + area[x][y +
			       incrementY].percentFull[area[x][y +
							       incrementY].
						       activeZ] - (z +
								   area[x][y].
								   percentFull
								   [z])) *
	  (CELL_HEIGHT / CELL_WIDTH);
      }

      tmpBedslope = sqrt (bedslopeX * bedslopeX + bedslopeY * bedslopeY);
      if (tmpBedslope > maxBedslope)
      {
	maxBedslope = tmpBedslope;
      }

      /*  Dispersion relation is solved using Dean and Dalrymple p. 72. Starting value for the loop is Kdw0 */
      Kdw0 =
	Raise (2 * PI / T,
	       2) * (1 / g) * (Raise (tanh (Raise (2 * PI / T, 2) * (D / g)),
				      -0.5));
      Kdw = Raise (2 * PI / T, 2) * (1 / g) * (1 / tanh (Kdw0 * D));
      Kaux = fabs (1 - Kdw0 / Kdw);
      while (Kaux >= 0.01)
      {
	Kdw0 = Kdw;
	Kdw = Raise (2 * PI / T, 2) * (1 / g) * (1 / tanh (Kdw0 * D));
	Kaux = fabs (1 - Kdw0 / Kdw);
      }
      OrbitalVel = (PI * waveHeight / T) * (1 / sinh (Kdw * D));
      /* ave %coarse over top 3*CELL_HEIGHT in a way that changes smoothly */
      EffectivePercentCoarse =
	(area[x][y].percentFull[z] * area[x][y].percentCoarse[z] +
	 area[x][y].percentCoarse[z - 1] + area[x][y].percentCoarse[z - 2] +
	 (1 - area[x][y].percentFull[z]) * area[x][y].percentCoarse[z -
								    3]) / 3;


/* 	Cf = 0.01 * (1 + A * EffectivePercentCoarse);*/

/* entrainment */
      OrbitalExcurs = 0.5 * waveHeight / sinh (Kdw * D);
      d50 =
	dfine * (1 - EffectivePercentCoarse) +
	dcoarse * (EffectivePercentCoarse);
      /*      Se=(1-EffectivePercentCoarse)*pow(dfine-d50,2) + EffectivePercentCoarse*pow(dcoarse-d50,2);
         Se=Raise(Se,0.5);
         d50Ripples = d50+ConvertToEffGrSz*Se; */
      d50Ripples =
	dfine * (1 - EffectivePercentCoarse) +
	ConvertToEffGrSz * dcoarse * (EffectivePercentCoarse);
      XX =
	4 * nu * Raise ((OrbitalExcurs * 2 * PI / (T * 1.0)),
			2) / (d50Ripples * Raise (1.65 * g * d50Ripples,
						  1.5));
      if (XX <= 2)
      {
	RippleAspectRatiohyst = 0.15 / Raise (XX, 0.11);
	RippleLengthCoarseMixhyst = 1.96 * OrbitalExcurs / Raise (XX, 0.28);
	RippleHeightCoarseMixhyst =
	  RippleLengthCoarseMixhyst * RippleAspectRatiohyst;
      }
      else
      {
	RippleAspectRatiohyst = 0.166 / Raise (XX, 0.24);
	RippleLengthCoarseMixhyst = 2.71 * OrbitalExcurs / Raise (XX, 0.75);
	RippleHeightCoarseMixhyst =
	  RippleLengthCoarseMixhyst * RippleAspectRatiohyst;
      }


      /* Bed concentration    
         Shieldscritical=Raise(dcoarse/d50,ezexps)*0.04; */
      Shieldscritical = 0.04;
      nikuk = 2.5 * d50;
      FrictionParam =
	exp (5.213 * Raise (nikuk / OrbitalExcurs, 0.194) - 5.977);
      /* Ripple hysterisis. Ripple dimensions are updated only if the shields parameter, 0.04, is exceeded    
         ShieldsParam = ( FrictionParam/(3.3*g*dcoarse) )*Raise(OrbitalExcurs*2*PI/T,2); */
      ShieldsParam =
	(FrictionParam / (3.3 * g * d50)) * Raise (OrbitalExcurs * 2 * PI / T,
						   2);
      ShieldsRipple =
	(FrictionParam / (3.3 * g * d50Ripples)) * Raise (OrbitalExcurs * 2 *
							  PI / T, 2);
      if (ShieldsRipple > Shieldscritical)
      {
	area[x][y].RippleAspectRatioCoarseMix = RippleAspectRatiohyst;
	area[x][y].RippleLengthCoarseMix = RippleLengthCoarseMixhyst;
	area[x][y].RippleHeightCoarseMix = RippleHeightCoarseMixhyst;
      }
      ModifiedShieldsParam =
	ShieldsParam / pow (1 - PI * area[x][y].RippleAspectRatioCoarseMix,
			    2);
      BedConc = 0.005 * Raise (ModifiedShieldsParam - Shieldscritical, 3);



/* velocity profile */
/* inside the boundary layer	Fredsoe and Deigaard, p. 59-60, deltabl from eq.2.45(p.25)*/
/* we are adapting eq.2.45 to account for the presence of ripples (use of km)	*/
/* concentration profile (Nielsen p. 258; Fredsoe and Deigaard p.305 and eq.10.24 for diffusivity -above is instead eq. 8.13)*/
/* we are assuming the current doesn't affect concentration profile and only advects sediment*/
/* we are neglecting the flux contribution to sediment fluxes below the ripple crest because:
1) we have noidea of the flow in there
2) we assume that contribution is included in the bedload prediction
					abl = OrbitalVel * T /(2*PI);
					deltabl = km * 0.09 * Raise(abl/km,0.82);
					kwobl = 30*deltabl/(exp(0.4*VRvr/Ufcbl));
					VRvr=Raise(uMagnitude*uMagnitude+vMagnitude*vMagnitude,0.5);
					Ufcbl = Raise(2/PI*fw*OrbitalVel*VRvr,0.5);*/

      auxexp2 =
	area[x][y].RippleLengthCoarseMix * OrbitalExcurs * 2 * PI / T *
	exp (1.5 - 4500 * dcoarse -
	     1.2 * log (OrbitalExcurs * 2 * PI / T / Wc));
      zeta0 =
	(28 * area[x][y].RippleHeightCoarseMix *
	 area[x][y].RippleAspectRatioCoarseMix + 2.5 * d50) / 30;
      Czbl = 0.0;
      Czobl = 0.0;
      ustar = 0.4 * uMagnitude / log (hustar / (zeta0));
      vstar = 0.4 * vMagnitude / log (hustar / (zeta0));
      qoblu = 0.0;
      qoblv = 0.0;
      qblu = 0.0;
      qblv = 0.0;
      stepbl = 0.025;

      if (BedConc > 0.0)
      {
	Czobl = 1.0;
	zeta = zeta0;
	zeta = area[x][y].RippleHeightCoarseMix;
	zeta1 = 1.0 * area[x][y].RippleHeightCoarseMix;
	Uobl = uMagnitude;
	Vobl = vMagnitude;
	while (Czobl > 0.000003)
	{
	  Czobl = BedConc * exp (-Wc * (zeta - zeta1) / auxexp2);
	  if (zeta < hustar)
	  {
	    Uobl = 1 / 0.4 * ustar * log (zeta / (zeta0));
	    Vobl = 1 / 0.4 * vstar * log (zeta / (zeta0));
	  }
	  else
	  {
	    Uobl = uMagnitude;
	    Vobl = vMagnitude;
	  }
/*	if (x==50 && y==50)
	{
		printf("zeta: %f Czobl: %f Uobl: %f \n ",zeta, Czobl, Uobl);
	}	*/
	  qoblu = qoblu + Czobl * Uobl * stepbl;
	  qoblv = qoblv + Czobl * Vobl * stepbl;
	  zeta = zeta + stepbl;
	}

	IntegratedConcx = qoblv;
	IntegratedConcy = qoblv;
      }
      else
      {
	IntegratedConcx = 0.0;
	IntegratedConcy = 0.0;
	coeffcrit = 0.0;
	if (coeffcrit < 0.5)
	{
	  printf
	    (" !!!!!!!!!!!!!!!!!!!!!ModifiedShieldsParam: %f Shieldscritical: %f \n",
	     ModifiedShieldsParam, Shieldscritical);
	  printf (" x: %i y: %i \n", x, y);
	}
      }

/*bedload (Van Rijn Principles of coastal morphology p.4.38) */
      qbx =
	coeffcrit * 9.1 * Raise (ModifiedShieldsParam - Shieldscritical,
				 1.78) * Raise (g * 1.65 * d50 * d50 * d50,
						0.5) * uMagnitude /
	fabs (uMagnitude);
      qby =
	coeffcrit * 9.1 * Raise (ModifiedShieldsParam - Shieldscritical,
				 1.78) * Raise (g * 1.65 * d50 * d50 * d50,
						0.5) * vMagnitude /
	fabs (vMagnitude);



      CD = Raise (karman / log (hustar / (zeta0)), 2);
      qsx1 = (ConvertToImmersedWt * IntegratedConcx);
      qsx2 =
	coeffcrit * coeffslope * coeff * CD * (1 / (5 * Wc)) *
	Raise (OrbitalVel, 5) * bedslopeX;

      qsy1 = (ConvertToImmersedWt * IntegratedConcy);
      qsy2 =
	coeffcrit * coeffslope * coeff * CD * (1 / (5 * Wc)) *
	Raise (OrbitalVel, 5) * bedslopeY;

      area[x][y].localFluxCoarseX =
	(EffectivePercentCoarse) * (qbx + qsx1 - qsx2);
      area[x][y].localFluxCoarseY =
	(EffectivePercentCoarse) * (qby + qsy1 - qsy2);

/*	if (x==50 && y==50)
	{
		printf("CCCqbx: %f qsx1: %f qsx2: %f \n ",qbx, qsx1, qsx2);
		printf("zeta0: %f area[x][y].RippleHeightCoarseMix: %f  \n ",zeta0, area[x][y].RippleHeightCoarseMix);
		printf("bedconc: %f coeffcrit: %f  int conc: %f \n ",BedConc, coeffcrit, IntegratedConcx);
		printf("ModifiedShieldsParam: %f Shieldscritical: %f   \n \n",ModifiedShieldsParam,Shieldscritical);
		
	}	*/


      /*effectiveProfileHeight = 0.001*D/Wc;
         effectiveProfileHeight= (2*OrbitalVel/Wc)*( ProfHtOverFine + ALPH * 1.00 );
         effectiveProfileHeight = 0.15;
         effectiveProfileHeight = (2*OrbitalVel/Wc)*( ProfHtOverFine +
         ALPH * EffectivePercentCoarse ); */
      effectiveProfileHeight = auxexp2 / Wc;


      if (uMagnitude > vMagnitude)
	percentDeposited = Wc * CELL_WIDTH / (uMagnitude *
					      effectiveProfileHeight);
      if (uMagnitude <= vMagnitude)
	percentDeposited = Wc * CELL_WIDTH / (vMagnitude *
					      effectiveProfileHeight);

      if (percentDeposited < 0)
      {
	printf (" \n Woo! This can not be! (percent deposited < 0) \n");
	percentDeposited = 0.0000000001;
	printf (" xcoarse= %d \n", x);
	printf (" ycoarse= %d \n", y);
	printf (" %g\n %g\n %g\n %g \n %g \n",
		area[x][y].percentCoarse[z],
		area[x][y].percentCoarse[z - 1],
		area[x][y].percentCoarse[z - 2],
		area[x][y].percentCoarse[z - 3], area[x][y].percentFull[z]);
	printf (" um= %f \n", uMagnitude);
	printf (" vm= %f \n", vMagnitude);
	printf (" eph= %f \n", effectiveProfileHeight);
	printf (" ProfHtOverFine= %f \n", ProfHtOverFine);
	printf (" EffectivePercentCoarse= %f \n", EffectivePercentCoarse);

      }
      if (percentDeposited > 1)
      {
	percentDeposited = 1;
      }

      excessSed = (localFluxInX - area[x][y].localFluxCoarseX)
	+ (localFluxInY - area[x][y].localFluxCoarseY)
	+ (excessSedInX + excessSedInY);
      /*think about negative excessOut... */
      area[x][y].excessCoarseSedOutX =
	(1 - percentDeposited) * excessSed
	* (uMagnitude / (uMagnitude + vMagnitude));
      area[x][y].excessCoarseSedOutY =
	(1 - percentDeposited) * excessSed
	* (vMagnitude / (uMagnitude + vMagnitude));
      area[x][y].coarseVolumeAdded =
	convertToVolume * percentDeposited * excessSed;

      if (DEBUG5 && area[x][y].volumeOut < 0)
      {
	printf ("\nBeep: vol. out < 0 ; %f\n  at %d,%d, localflux %f, \n \n ",
		area[x][y].volumeOut, x, y, localFlux);
	printf ("percentcoarse[z]=%f\n", area[x][y].percentCoarse[z]);
	printf ("percentcoarse[z-1]=%f\n", area[x][y].percentCoarse[z - 1]);
	printf ("percentcoarse[z-2]=%f\n", area[x][y].percentCoarse[z - 2]);
      }

      /*record how much susp. load to pass in the other side next iteration */
      if (x == endX - incrementX)
      {
	excessOutofIterCoarseX[y] = area[x][y].excessCoarseSedOutX;
	localFluxOutofIterCoarseX[y] = area[x][y].localFluxCoarseX;
      }
      if (y == endY - incrementY)
      {
	excessOutofIterCoarseY[x] = area[x][y].excessCoarseSedOutY;
	localFluxOutofIterCoarseY[x] = area[x][y].localFluxCoarseY;
      }
/*if (x==3 && y==78)
{ 
  printf(" CCCqsx1: %f qsx2: %f \n",qsx1,qsx2);
  printf(" CCC_D: %f z: %i \n",D,z);
  printf(" CCC_zeta: %f z: %i \n",zeta,z);
  printf(" CCC_CD: %f bedslopeX: %f \n",CD,bedslopeX);
  printf(" CCCModifiedShieldsParam: %f Shieldscritical: %f \n",ModifiedShieldsParam,Shieldscritical);
  printf(" CCCXX: %f OrbitalExcurs: %f \n",XX,OrbitalExcurs);
  printf(" CCCBedConc: %f diffusivity: %f \n",BedConc,auxexp2);
  printf(" CCCRippleL: %f RippleRatio: %f \n",area[x][y].RippleLengthCoarseMix,area[x][y].RippleAspectRatioCoarseMix);
  printf(" CCCorb vel: %f EffPercentCoarse: %f \n",OrbitalVel,EffectivePercentCoarse);
  printf("CCCx,y: %d %d vMagnitude: %g  effProfHt: %g \n",x,y, vMagnitude, effectiveProfileHeight);	
  printf("   CCCconvertToVol: %g percent deposited: %g excessSed %g\n",convertToVolume,
						percentDeposited,excessSed);
  printf("   CCClocalFluxInX: %g localFluxcoarseX(x,y): %g \n", localFluxInX, area[x][y].localFluxCoarseX);
  printf("   CCClocalFluxInY: %g localFluxcoarseY(x,y): %g \n", localFluxInY, area[x][y].localFluxCoarseY);
  printf("   CCCexcessSedInX: %g excessSedInY: %g \n",excessSedInX, excessSedInY);
  printf("   CCCcoarseVolAdded: %g \n", area[x][y].coarseVolumeAdded);
  
}*/


    }
}




void
AdjustCells ()
/* this function adjusts the information in the cells based on the new
   data, as derived from sediment transport */
{
  int x, y, z;
  double newEntry, newFullness, underflow, overflow;
  double oldCoarse;
  int level;

  for (x = 0; x < Xmax; x++)
    for (y = 0; y < Ymax; y++)
    {
      z = area[x][y].activeZ;
      area[x][y].fineVolumeAdded +=
	(AggRateFine / SecondsPerYear) * timeStep * CELL_WIDTH * CELL_WIDTH;
      area[x][y].coarseVolumeAdded +=
	(AggRateCoarse / SecondsPerYear) * timeStep * CELL_WIDTH * CELL_WIDTH;
      newEntry =
	((area[x][y].fineVolumeAdded +
	  area[x][y].coarseVolumeAdded) / CELL_VOLUME);
      newFullness = newEntry + area[x][y].percentFull[z];

      if ((newFullness >= 0) && (newFullness <= 1))
      {
	oldCoarse = area[x][y].percentCoarse[z] *
	  area[x][y].percentFull[z] * CELL_VOLUME;

	area[x][y].percentCoarse[z] =
	  (area[x][y].coarseVolumeAdded +
	   oldCoarse) / (newFullness * CELL_VOLUME);

	if (area[x][y].percentCoarse[z] < 0)
	{
	  level = z;
	  do
	  {
	    area[x][y].percentCoarse[level - 1] =
	      area[x][y].percentCoarse[level - 1]
	      + area[x][y].percentCoarse[level] * newFullness;
	    area[x][y].percentCoarse[level] = 0;
	    level = level - 1;
	  }
	  while (area[x][y].percentCoarse[level] < 0 && level > 0);
	}


	else if (area[x][y].percentCoarse[z] > 1)
	{
	  level = z;
	  do
	  {
	    area[x][y].percentCoarse[level - 1] =
	      area[x][y].percentCoarse[level - 1]
	      + (area[x][y].percentCoarse[level] - 1) * newFullness;
	    area[x][y].percentCoarse[level] = 1;
	    level = level - 1;
	  }
	  while (area[x][y].percentCoarse[level] > 1 && level > 0);
	}

	area[x][y].percentFull[z] = newFullness;
      }

      else if (newFullness > 1)
      {
	overflow = newFullness - 1;

	if (DEBUG4)
	  printf ("OVERFLOW \n");

	oldCoarse = area[x][y].percentCoarse[z] *
	  area[x][y].percentFull[z] * CELL_VOLUME;

	area[x][y].percentCoarse[z] =
	  (area[x][y].coarseVolumeAdded +
	   oldCoarse) / (newFullness * CELL_VOLUME);

	if (area[x][y].percentCoarse[z] < 0)
	{
	  level = z;
	  do
	  {
	    area[x][y].percentCoarse[level - 1] =
	      area[x][y].percentCoarse[level - 1]
	      + area[x][y].percentCoarse[level] * newFullness;
	    area[x][y].percentCoarse[level] = 0;
	    level = level - 1;
	  }
	  while (area[x][y].percentCoarse[level] < 0 && level > 0);
	}


	else if (area[x][y].percentCoarse[z] > 1)
	{
	  level = z;
	  do
	  {
	    area[x][y].percentCoarse[level - 1] =
	      area[x][y].percentCoarse[level - 1]
	      + (area[x][y].percentCoarse[level] - 1) * newFullness;
	    area[x][y].percentCoarse[level] = 1;
	    level = level - 1;
	  }
	  while (area[x][y].percentCoarse[level] > 1 && level > 0);
	}
	area[x][y].percentFull[z] = 1;
	area[x][y].activeZ = area[x][y].activeZ + 1;	/* move up one */
	z = area[x][y].activeZ;	/* set z again */

	area[x][y].percentFull[z] = overflow;
	area[x][y].percentCoarse[z] = area[x][y].percentCoarse[z - 1];

      }
      else if (newFullness < 0)
      {
	underflow = newFullness + 1;	/* newFullness is negative */
	if (DEBUG4)
	  printf ("UNDERFLOW \n");

/* if ((x==50)&&(y==50))  
{           printf("newFullness= %f \n",newFullness);
            printf("finevol= %f \n",area[x][y].fineVolumeAdded);
            printf("corsvol= %f \n",area[x][y].coarseVolumeAdded);
} */

	oldCoarse =
	  area[x][y].percentCoarse[z] * area[x][y].percentFull[z] *
	  CELL_VOLUME + area[x][y].percentCoarse[z - 1] * CELL_VOLUME;

	area[x][y].percentCoarse[z] =
	  (area[x][y].coarseVolumeAdded +
	   oldCoarse) / (underflow * CELL_VOLUME);
	area[x][y].activeZ = area[x][y].activeZ - 1;	/* move down one */

	z = area[x][y].activeZ;	/* set z again */

	area[x][y].percentCoarse[z] = area[x][y].percentCoarse[z + 1];

	area[x][y].percentCoarse[z + 1] = 0;

	area[x][y].percentFull[z + 1] = 0;

	area[x][y].percentFull[z] = underflow;

	if (area[x][y].percentCoarse[z] < 0)
	{
	  level = z;
	  do
	  {
	    area[x][y].percentCoarse[level - 1] =
	      area[x][y].percentCoarse[level - 1]
	      + area[x][y].percentCoarse[level] * underflow;
	    area[x][y].percentCoarse[level] = 0;
	    level = level - 1;
	  }
	  while (area[x][y].percentCoarse[level] < 0 && level > 0);
	}


	else if (area[x][y].percentCoarse[z] > 1)
	{
	  level = z;
	  do
	  {
	    area[x][y].percentCoarse[level - 1] =
	      area[x][y].percentCoarse[level - 1]
	      + (area[x][y].percentCoarse[level] - 1) * underflow;
	    area[x][y].percentCoarse[level] = 1;
	    level = level - 1;
	  }
	  while (area[x][y].percentCoarse[level] > 1 && level > 0);
	}



      }

      /* adjust sediment height and water depth
         necessary? *//* yes, eventually water depth (for waves) */
    }
}


void
ZeroVars ()
/*resets a couple of things for a fresh interation*/
{
  int x, y;

  for (y = 0; y < 2 * Ymax; y++)
    for (x = 0; x < 2 * Xmax; x++)
    {
      area[x][y].fineVolumeAdded = area[x][y].coarseVolumeAdded = 0;
    }
}

void
DoIterationDummy ()
     /* controls iteration by looking at what direction the flow is
        and then calling the proper sediment transport function 
        DUMMY VERSION DOESN'T ADJUST CELLS, SO THE SED TRANS OUT OF ITERATION CAN ADJUST 
        TO NEW DIRECTION FOR A FEW ITERATIONS FIRST */
{

  SedTransFine ();
  SedTransCoarse ();

}



void
DoIteration ()
     /* controls iteration by looking at what direction the flow is
        and then calling the proper sediment transport function */
{

  SedTransFine ();
  SedTransCoarse ();
  AdjustCells ();
  UpdateForcingClock ();
  UpdateRunTimeClock ();
  /* NumFramesSaved++;
     PrintToFile();
     printf("\nFrame %i saved\n", NumFramesSaved); */


}


double
RunTimeClock ()
{
  return totalElapsedTime;	/* // units of seconds */
}				/*  //RunTimeClock */


double
ForcingClock ()
{
  return timeSinceForcingUpdate;	/* // units of seconds */
}				/*  // ForcingClock */


void
UpdateRunTimeClock ()
{
  totalElapsedTime = totalElapsedTime + timeStep;	/* // units of seconds */
}				/*  // UpdateRunTimeClock */


void
UpdateForcingClock ()
{
  if (timeSinceForcingUpdate >= FORCING_DURATION)
    timeSinceForcingUpdate = 0.0;
  else
    timeSinceForcingUpdate = timeSinceForcingUpdate + timeStep;
  /*  // again, units of seconds */
}				/*  // UpdateForcingClock */

void
WriteForcingConditions ()
{
  fprintf (SaveForcing, "%i %f %f %f %f\n",
	   NumFramesSaved, RunTimeClock (), waveHeight,
	   currentDirectionX * currentVelocityX,
	   currentDirectionY * currentVelocityY);
  return;
}				/*  // WriteForcingConditions */
