#include "mt19937.h"
#include "epidemic.h"

/* input parameters */
int npart;      // number of particles / random walkers
int nsteps;     // number of time steps
float boxsize;  // size of the simulation box
float ri;       // radius of influense
float pi;       // transmission probability
int Trec;     // recovery time period (in steps)

/* optional parameters */
float dt   = 0.1;   // time-step
int repeat = 1;     // number of repetitions
int subdiv = 10;    // number of subdivisions of the region

int writeCounts = 1; // write the counts as a file

/* calculated parameters */
int ncells;     // number of cells = subdiv^2
float cellsize; // size of a cell in the box


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

    rw->instant = 0;

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
    for (int k = 0; k < npart; k++)
    {
        float xk = rw->x[k] + rw->vx[k] * dt;
        float yk = rw->y[k] + rw->vy[k] * dt;

        // rw->vx[k] += ax[k] * dt;
        // rw->vy[k] += ay[k] * dt;

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

void simulation(struct Stats *stats)
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

void __simulation()
{
    struct Stats *stats;
    stats = (struct Stats*) malloc(nsteps * sizeof(struct Stats));

    for (int rep = 0; rep < repeat; rep++)
    {
        struct Stats *__stats;
        __stats = (struct Stats*) malloc(nsteps * sizeof(struct Stats));

        __simulation(__stats);

        for (int j = 0; j < nsteps; j++)
        {
            stats[j].s += __stats[j].s;
            stats[j].i += __stats[j].i;
            stats[j].r += __stats[j].r;
        }
    }

    // take averages
    for (int j = 0; j < nsteps; j++)
    {
        stats[j].s /= (float) repeat;
        stats[j].i /= (float) repeat;
        stats[j].r /= (float) repeat;
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

    npart = 500;
    boxsize = 10.0;
    Trec   = 10;
    pi     = 0.1;
    ri     = 1;
    nsteps = 50;
    repeat = 10;

    __simulation();

    return 0;
}