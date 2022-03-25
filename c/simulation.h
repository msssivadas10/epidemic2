#ifndef SIMUTILS_H
#define SIMUTILS_H 2

/* all simulation variables */

extern int Nparticles, Nsteps, Repeat, SimulationKey; 
extern float Boxsize, TimeStep;
extern int Subdivisions;

extern int RecoveryTime; 
extern float InfectionProbability; 
extern float InfectionRadius;

extern float ForceStrength, ForceSpeed;
extern float ForceStrengthStart, ForceStrengthStop;
extern int ForceStrengthStep;

extern float InfectionRadiusStart, InfectionRadiusStop;
extern int InfectionRadiusStep;

extern char OutputDir[256];

/* Read input parameters from a file. */
void readParameters(char *fname);

#endif // SIMUTILS_H