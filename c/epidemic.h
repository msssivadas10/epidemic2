#ifndef EPIDEMIC_H
#define EPIDEMIC_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define min(x, y) ((x) < (y) ? (x) : (y))
#define max(x, y) ((x) < (y) : (y) : (x))
#define sqr(x) ((x) * (x))
#define hypot2(x, y) ((x) * (x) + (y) * (y))

typedef enum {EXPOSED, INFECTED, RECOVERED} status;

/* store the configuration of the particles */
struct Particles
{
    /* data */
    int instant;    // time instant: time = instant * dt
    float *x, *y;   // positions
    float *vx, *vy; // velocities
    status *state;  // infection state 
    int *t2rec;     // time left to recover
    int **grid;
};

/* Store the state at an instant */
struct Stats
{
    float s, i, r; // number of exposed, infected and recovered
};

/* Initialise a particle system (random walkers). */
void initialiseParticles(struct Particles *rw);

/* Free memory. */
void freeParticles(struct Particles *rw);

/* Put the particles on a hash-grid. */
void prepareGrid(struct Particles *rw);

/* Update the particles configuration. */
void updateParticles(struct Particles *rw);

/* Infect a random particle. */
void infectRandom(struct Particles *rw);

/* Get the particles nearest to a point, within a radius. */
void getNearest(struct Particles *rw, float x, float y, float r, int *nid, int *nn);

/* Get the counts at present instant. */
void getCurrentCounts(struct Particles *rw, struct Stats *stats);

/* Run a (normal) simulation once */
void simulation(struct Stats *stats);


#endif // EPIDEMIC_H