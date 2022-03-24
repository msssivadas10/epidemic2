#include <stdio.h>
#include <stdlib.h>
#include "epidemic.h"

int npart;
float boxsize;
int Trec;
float pi;
float ri;
int nsteps;
int repeat;
int subdiv;
float dt;

void simulation()
{
    struct Stats *stats;
    stats = (struct Stats*) malloc(nsteps * sizeof(struct Stats)); // store average stats 

    for (int rep = 0; rep < repeat; rep++)
    {
        struct Stats *__stats;
        __stats = (struct Stats*) malloc(nsteps * sizeof(struct Stats));

        __simulation(__stats);

        // if (rep == 0)
        // {
        //     for (int j = 0; j < nsteps; j++)
        //     {
        //         stats[j].s = 0.0; stats[j].i = 0.0; stats[j].r = 0.0; 
        //     }
        // }

        for (int j = 0; j < nsteps; j++)
        {
            stats[j].s += (__stats[j].s - stats[j].s) / (float) (rep + 1);
            stats[j].i += (__stats[j].i - stats[j].i) / (float) (rep + 1);
            stats[j].r += (__stats[j].r - stats[j].r) / (float) (rep + 1);
        }
    }

    // write output to file
    FILE *file = fopen("output.txt", "w");
    if (file == NULL)
    {
        printf("Error opening file.\n");
        exit(0);
    }

    for (int t = 0; t < nsteps; t++)
        fprintf(file, "%d\t%16.8f\t%16.8f\t%16.8f\n", t, stats[t].s, stats[t].i, stats[t].r);
    
    fclose(file);
}

int main()
{

    npart   = 500;
    boxsize = 10.0;
    Trec    = 10;
    pi      = 0.1;
    ri      = 0.1;
    nsteps  = 100;
    repeat  = 10;
    subdiv  = 10;
    dt      = 0.1;

    setupSimulation(npart, boxsize, nsteps, subdiv, dt);

    setInfectionRadius(ri);
    setInfectionProbability(pi);
    setRecoveryTime(Trec);
    switchForce(0);

    simulation();

    return 0;
}
