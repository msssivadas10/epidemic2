#include "epidemic.h"

/* input parameters */
int npart;                // number of particles / random walkers
int nsteps;               // number of time steps
float boxsize;            // size of the simulation box
float ri;                 // radius of influense
float pi;                 // transmission probability
int Trec;                 // recovery time period (in steps)

/* optional parameters */
float dt         = 0.1;   // time-step
int repeat       = 1;     // number of repetitions
int subdiv       = 10;    // number of subdivisions of the region
int forceOn      = 0;     // state of the force field
int movingForce  = 0;     //  force field is moving 
float forceSpeed = 1.0;   // speed of the force field
float forceValue = 1.0;   // strength of the force
float forceRad   = 0.05;  // force radius parameter
float __soft     = 0.1;

/* calculated/other parameters */
int ncells;               // number of cells = subdiv^2
float cellsize;           // size of a cell in the box
float force_x,  force_y ; // position of the central force
float force_vx, force_vy; // velocity of force field

void initialiseParticles(struct Particles *rw)
{
    rw->x     = (float*) malloc(npart * sizeof(float));
    rw->y     = (float*) malloc(npart * sizeof(float));

    rw->vx    = (float*) malloc(npart * sizeof(float));
    rw->vy    = (float*) malloc(npart * sizeof(float));

    rw->state = (status*) malloc(npart * sizeof(int));
    rw->t2rec = (int*) malloc(npart * sizeof(int));

    /* put the particles random position and velocity */

    struct timespec __ts;
    clock_gettime(CLOCK_REALTIME, &__ts);
    size_t __seed = (size_t) __ts.tv_sec * 1000000000L + (size_t) __ts.tv_nsec;
    initGenerator(__seed);

    for (int k = 0; k < npart; k++)
    {
        // random positions in [0, boxsize]
        rw->x[k] = randf() * boxsize;
        rw->y[k] = randf() * boxsize;

        // random velocity in [-1, 1]
        rw->vx[k] = 2.0 * randf() - 1.0;
        rw->vy[k] = 2.0 * randf() - 1.0;

        // all are initially not infected
        rw->state[k] = EXPOSED;
        rw->t2rec[k] = 0;
    }

    // initialise the force field
    force_x  = boxsize * 0.5; 
    force_y  = boxsize * 0.5; 
    if (movingForce)
    {
        force_vx = (2.0 * randf() - 1.0) * forceSpeed;
        force_vy = (2.0 * randf() - 1.0) * forceSpeed;
    }

    rw->instant = 0;
    rw->npart   = npart;

    // make a grid and place the particles on it
    ncells   = subdiv * subdiv;
    cellsize = boxsize / subdiv;
    rw->grid = (int**) malloc(ncells * sizeof(int*));
    for (int i = 0; i < ncells; i++)
    {
        rw->grid[i] = (int*) malloc((npart + 1) * sizeof(int));
    }

    prepareGrid(rw);
}

void freeParticles(struct Particles *rw)
{
    free(rw->x);
    free(rw->y);
    free(rw->vx);
    free(rw->vy);
    free(rw->state);
    free(rw->t2rec);
    free(rw->grid);
}

void prepareGrid(struct Particles *rw)
{
    // reset all the grid cells to empty
    for (int i = 0; i < ncells; i++)
        rw->grid[i][0] = 0;

    for (int k = 0; k < npart; k++)
    {
        int i = rw->x[k] / cellsize,
            j = rw->y[k] / cellsize;
        int c = i + j * subdiv;
        int n = rw->grid[c][0];

        rw->grid[c][n+1] = k;
        rw->grid[c][0]++;
    }
}

void getNearest(struct Particles *rw, float x, float y, float r, int *nid, int *nn)
{
    // get the cell contatining the point
    int ci   = x / cellsize,
        cj   = y / cellsize,
        step = ceilf(r / cellsize);
    float r2 = r * r;

    int __nn = 0;
    for (int i = ci - step; i < ci + step + 1; i++)
    {
        if (i < 0 || i >= subdiv)
            continue; // cell outside the box

        for (int j = cj - step; j < cj + step + 1; j++)
        {
            if (j < 0 || j >= subdiv)
                continue; // cell outside the box

            int cell = i + j * subdiv;    // cell to look
            int num  = rw->grid[cell][0]; // number of particles in this cell

            if (!num)
                continue; // no particles in this cell

            for (int k = 0; k < num; k++)
            {
                int   p  = rw->grid[cell][k+1];

                // check if this particle is near the point
                if ( hypot2(rw->x[p] - x, rw->y[p] - y) <= r2 )
                {
                    nid[__nn] = p;
                    __nn++;
                }
            }
        }
    }

    *nn = __nn;
}

void getAcceleration(struct Particles *rw, float *ax, float *ay)
{
    // get the acceleration due to force field
    for (int k = 0; k < npart; k++)
    {
        if (!forceOn)
        {
            ax[k] = 0.0; ay[k] = 0.0;
            continue;
        }
        float dx = rw->x[k] - force_x,
              dy = rw->y[k] - force_y;
        float rk = sqrtf(dx * dx + dy * dy + __soft);
        float ak = - forceValue * expf( -sqr( rk / boxsize ) / forceRad ) / rk;

        ax[k] = ak * dx; ay[k] = ak * dy;
    }

    // TODO: get acceleration due to interaction
}

void updateParticles(struct Particles *rw)
{
    /* transmit the infection */
    int *nearest = (int*) malloc(npart * sizeof(int)); // to store nearest neighbours
    int ncount = 0; // number of neighbours

    for (int k = 0; k < npart; k++)
    {
        if (rw->state[k] != INFECTED)
            continue; // particle has no infection
        getNearest(rw, rw->x[k], rw->y[k], ri, nearest, &ncount);
        if (!ncount)
            continue; // no nearest particles

        for (int p = 0; p < ncount; p++)
        {
            int q = nearest[p];

            if (rw->state[q] != EXPOSED)
                continue; // infect only exposed particles

            // infect with a probability
            if (randf() < pi)
            {
                rw->state[q] = INFECTED;
                rw->t2rec[q] = Trec;
            }
        }
    }

    free(nearest);

    /* update particle positions and velocity */
    float *ax, *ay; // acceleration
    ax = (float*) malloc(npart * sizeof(float));
    ay = (float*) malloc(npart * sizeof(float));
    if (forceOn)
        getAcceleration(rw, ax, ay);

    for (int k = 0; k < npart; k++)
    {
        float xk = rw->x[k] + rw->vx[k] * dt;
        float yk = rw->y[k] + rw->vy[k] * dt;
        
        if (forceOn)
        {
            // force is on: update velocity

            rw->vx[k] += ax[k] * dt;
            rw->vy[k] += ay[k] * dt;
        }

        // make particles reflect at boundary
        if (xk < 0.0)
        {
            xk = -xk; rw->vx[k] = -rw->vx[k];
        }
        if (xk > boxsize)
        {
            xk = 2 * boxsize - xk; rw->vx[k] = -rw->vx[k];
        }
        if (yk < 0.0)
        {
            yk = -yk; rw->vy[k] = -rw->vy[k];
        }
        if (yk > boxsize)
        {
            yk = 2 * boxsize - yk; rw->vy[k] = -rw->vy[k];
        }
        
        rw->x[k] = xk; rw->y[k] = yk;

        // update recovery status
        if (rw->state[k] == INFECTED)
        {
            rw->t2rec[k]--;
            
            if (rw->t2rec[k] < 0)
            {
                rw->t2rec[k] = 0;
                rw->state[k] = RECOVERED;
            }
        }
    }

    // update the force field
    if (movingForce)
    {
        force_x += force_vx * dt;
        force_y += force_vy * dt;
    }

    rw->instant++;

    // re-compute the grid for current configuration
    prepareGrid(rw);
}

void infectRandom(struct Particles *rw)
{
    int k = randi() % npart; // select a random particle
    
    rw->state[k] = INFECTED;
    rw->t2rec[k] = Trec;
}

void getCurrentCounts(struct Particles *rw, struct Stats *stats)
{
    stats->s = 0;
    stats->r = 0;
    stats->i = 0;
    for (int k = 0; k < npart; k++)
    {
        if (rw->state[k] == EXPOSED) 
            stats->s++;
        else if (rw->state[k] == INFECTED) 
            stats->i++;
        else if (rw->state[k] == RECOVERED) 
            stats->r++;
    }
}


/* set parameters (constants) */

void setParticleNumber(int __npart)
{
    if (__npart <= 0)
    {
        printf("Error: number of particles must be positive\n"); exit(0);
    }
    npart = __npart;
}

void setBoxsize(float __boxsize)
{
    if (__boxsize <= 1.0)
    {
        printf("Error: Boxsize must be atleast 1\n"); exit(0);
    }
    boxsize = __boxsize;
}

void setSteps(int __nsteps)
{
    if (__nsteps <= 0)
    {
        printf("Error: number of steps must be positive\n"); exit(0);
    }
    nsteps = __nsteps;
}

void setTimestep(float __dt)
{
    if (__dt < 0.0)
    {
        printf("Error: time step must be positive\n"); exit(0);
    }
    dt = __dt;
}

void setSubdivisions(int __subdiv)
{
    if (__subdiv < 1)
    {
        printf("Error: number of subdivisions must be atleast 1\n"); exit(0);
    }
    subdiv = __subdiv;
    ncells = subdiv * subdiv;
}

void setupSimulation(int __npart, float __boxsize, int __nsteps, int __subdiv, float __dt)
{
    setParticleNumber(__npart);
    setBoxsize(__boxsize);
    setSteps(__nsteps);
    setSubdivisions(__subdiv);
    setTimestep(__dt);
}

/* set parameters (variables) */

void setInfectionRadius(float __ri)
{
    if (__ri < 0.0)
    {
        printf("Error: infection radius must be positive\n"); exit(0);
    }
    ri = __ri * boxsize;
}

void setInfectionProbability(float __pi)
{
    if (__pi < 0.0 || __pi > 1.0)
    {
        printf("Error: infection probability must be in the range [0, 1]\n"); exit(0);
    }
    pi = __pi;
}

void setRecoveryTime(int __Trec)
{
    if (__Trec < 0)
    {
        printf("Error: recovery period must be positive\n"); exit(0);
    }
    Trec = __Trec;
}

void setForceStrength(float __value)
{
    forceValue = __value;
}

void setForceSpeed(float __value)
{
    forceSpeed = __value;
}

void setForceRadius(float __value)
{
    forceRad = __value;
}

void lockForce(int __value)
{
    movingForce = __value;
}

void switchForce(int __value)
{
    forceOn = __value;
}

void __simulation(struct Stats *stats)
{
    struct Particles rw;

    initialiseParticles(&rw);

    infectRandom(&rw);

    int step = 0;
    while (step < nsteps)
    {
        struct Stats __stats;
        getCurrentCounts(&rw, &__stats);
        stats[step] = __stats;

        if (__stats.i == 0.0)
        {
            // printf("No infected left: ending run at step %d\n", step);
            break;
        }

        updateParticles(&rw);

        step++;
    }

    for (int i = step; i < nsteps; i++)
    {
        stats[i].s = stats[step].s;
        stats[i].i = 0.0;
        stats[i].r = stats[step].r;
    }

    freeParticles(&rw);
}

