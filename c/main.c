#include <stdio.h>
#include <stdlib.h>
#include "epidemic.h"
#include "simulation.h"

int Nparticles;
float Boxsize;
int Nsteps;
int RecoveryTime;
float InfectionProbability;
float InfectionRadius;
int Repeat                 = 1;
int Subdivisions           = 10;
float TimeStep             = 0.1;
float ForceStrength        = 0.0, 
      ForceSpeed           = 0.0;
int SimulationKey          = 1;
float ForceStrengthStart   = 0.0, 
      ForceStrengthStop    = 100.0;
int ForceStrengthStep      = 21;
float InfectionRadiusStart = 0.0, 
      InfectionRadiusStop  = 0.2;
int InfectionRadiusStep    = 21;
char OutputDir[256];

/* Simulate the particles system to get averaged statistics */
void simulation1(void)
{
    /* setup simulation */
    setupSimulation(Nparticles, Boxsize, Nsteps, Subdivisions, TimeStep);

    setInfectionRadius(InfectionRadius);
    setInfectionProbability(InfectionProbability);
    setRecoveryTime(RecoveryTime);

    if (fabsf(ForceStrength) < 1e-08)
    {
        switchForce(0);
    }
    else
    {
        switchForce(1);
        setForceStrength(ForceStrength);
        setForceSpeed(ForceSpeed);
    }

    struct Stats *stats;
    stats = (struct Stats*) malloc(Nsteps * sizeof(struct Stats)); // store average stats 

    for (int rep = 0; rep < Repeat; rep++)
    {
        struct Stats *__stats;
        __stats = (struct Stats*) malloc(Nsteps * sizeof(struct Stats));

        __simulation(__stats);

        // if (rep == 0)
        // {
        //     for (int j = 0; j < nsteps; j++)
        //     {
        //         stats[j].s = 0.0; stats[j].i = 0.0; stats[j].r = 0.0; 
        //     }
        // }

        for (int j = 0; j < Nsteps; j++)
        {
            stats[j].s += (__stats[j].s - stats[j].s) / (float) (rep + 1);
            stats[j].i += (__stats[j].i - stats[j].i) / (float) (rep + 1);
            stats[j].r += (__stats[j].r - stats[j].r) / (float) (rep + 1);
        }
    }

    /* write output to file */
    FILE *file = fopen("output.txt", "w");
    if (file == NULL)
    {
        printf("Error opening file.\n");
        exit(0);
    }

    for (int t = 0; t < Nsteps; t++)
        fprintf(file, "%d\t%16.8f\t%16.8f\t%16.8f\n", t, stats[t].s, stats[t].i, stats[t].r);
    
    fclose(file);
}

/* Simulation to find maximum infection as function of force strength and infection radius. */
void simulation2()
{
    ;
}

int main()
{
    simulation1();

    return 0;
}
