#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

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

char* strip(char *__str)
{
    // Ref.: <https://stackoverflow.com/questions/122616how-do-i-trim-leading-trailing-whitespace-in-a-standard-way>
    char *__end;

    // strip leading spaces
    while (isspace((unsigned char) *__str)) __str++;

    if (*__str == 0) return __str; // all are spaces

    // strip trailing spaces
    __end = __str + strlen(__str) - 1;
    while (__end > __str && isspace((unsigned char) *__end)) __end--;

    __end[1] = 0; // write the null charecter

    return __str;
}

void readParameters(char *fname)
{
    enum {INTEGER, FLOAT, STRING} __types;
    enum {UNSET, SET} __status;
    const int NKEYS = 64; // maximum number of keys

    /* variable table */
    char keys[NKEYS][32];  // variable keys
    int  type[NKEYS];      // variable types 
    int  status[NKEYS];    // status: set or unset
    void *address[NKEYS];  // variable addresses

    int nkeys = 0;

    strcpy(keys[nkeys], "Nparticles");
    address[nkeys] = &Nparticles;
    status[nkeys]  = UNSET;     
    type[nkeys++]  = INTEGER;

    strcpy(keys[nkeys], "Boxsize");
    address[nkeys] = &Boxsize;
    status[nkeys]  = UNSET;     
    type[nkeys++]  = FLOAT;

    strcpy(keys[nkeys], "Nsteps");
    address[nkeys] = &Nsteps;
    status[nkeys]  = UNSET;     
    type[nkeys++]  = INTEGER;

    strcpy(keys[nkeys], "RecoveryTime");
    address[nkeys] = &RecoveryTime;
    status[nkeys]  = UNSET;     
    type[nkeys++]  = INTEGER;

    strcpy(keys[nkeys], "InfectionProbability");
    address[nkeys] = &InfectionProbability;
    status[nkeys]  = UNSET;     
    type[nkeys++]  = FLOAT;

    strcpy(keys[nkeys], "InfectionRadius");
    address[nkeys] = &InfectionRadius;
    status[nkeys]  = UNSET;     
    type[nkeys++]  = FLOAT;

    strcpy(keys[nkeys], "Repeat");
    address[nkeys] = &Repeat;
    status[nkeys]  = SET;     
    type[nkeys++]  = INTEGER;

    strcpy(keys[nkeys], "Subdivisions");
    address[nkeys] = &Subdivisions;
    status[nkeys]  = SET;     
    type[nkeys++]  = INTEGER;

    strcpy(keys[nkeys], "TimeStep");
    address[nkeys] = &TimeStep;
    status[nkeys]  = SET;     
    type[nkeys++]  = FLOAT;

    strcpy(keys[nkeys], "ForceStrength");
    address[nkeys] = &ForceStrength;
    status[nkeys]  = SET;     
    type[nkeys++]  = FLOAT;

    strcpy(keys[nkeys], "ForceSpeed");
    address[nkeys] = &ForceSpeed;
    status[nkeys]  = SET;     
    type[nkeys++]  = FLOAT;

    strcpy(keys[nkeys], "SimulationKey");
    address[nkeys] = &SimulationKey;
    status[nkeys]  = SET;     
    type[nkeys++]  = INTEGER;

    strcpy(keys[nkeys], "ForceStrengthStart");
    address[nkeys] = &ForceStrengthStart;
    status[nkeys]  = SET;     
    type[nkeys++]  = FLOAT;

    strcpy(keys[nkeys], "ForceStrengthStop");
    address[nkeys] = &ForceStrengthStop;
    status[nkeys]  = SET;     
    type[nkeys++]  = FLOAT;

    strcpy(keys[nkeys], "ForceStrengthStep");
    address[nkeys] = &ForceStrengthStep;
    status[nkeys]  = SET;     
    type[nkeys++]  = INTEGER;

    strcpy(keys[nkeys], "InfectionRadiusStart");
    address[nkeys] = &InfectionRadiusStart;
    status[nkeys]  = SET;     
    type[nkeys++]  = FLOAT;

    strcpy(keys[nkeys], "InfectionRadiusStop");
    address[nkeys] = &InfectionRadiusStop;
    status[nkeys]  = SET;     
    type[nkeys++]  = FLOAT;

    strcpy(keys[nkeys], "InfectionRadiusStep");
    address[nkeys] = &InfectionRadiusStep;
    status[nkeys]  = SET;     
    type[nkeys++]  = INTEGER;

    strcpy(keys[nkeys], "OutputDir");
    address[nkeys] = &OutputDir;
    status[nkeys]  = UNSET;     
    type[nkeys++]  = STRING;


    /* read parameters from file */
    FILE *fd = fopen(fname, "r");
    if (fd == NULL)
    {
        printf("Error: error opening file '%s', no such file or directory.\n", fname);
        exit(0);
    }

    char buf[256];
    char buf1[32], buf2[32], buf3[32], *__buf1;

    while (!feof(fd))
    {
        if (fgets(buf, 256, fd) != NULL)
        {
            if ( sscanf(buf, "%[^=]=%[^;\n];%s", buf1, buf2, buf3) < 2 )
                continue;

            if (buf1[0] == ';')
                continue;
            

            __buf1 = strip(buf1);
            for (int j = 0; j < nkeys; j++)
            {
                if ( strcmp(__buf1, keys[j]) == 0 )
                {
                    if (type[j] == INTEGER)
                    {
                        *((int*) address[j]) = atoi(buf2);
                    }
                    else if (type[j] == FLOAT)
                    {
                        *((float*) address[j]) = (float) atof(buf2); 
                    }
                    else if (type[j] == STRING)
                    {
                        strcpy(address[j], buf2);
                    }
                    status[j] = SET;
                    break;
                }
            }
        }
    }
    
    fclose(fd);

    // check if all parameters are set
    for (int j = 0; j < nkeys; j++)
    {
        if (status[j] == UNSET)
        {
            printf("Error: value of '%s' is not set.\n", keys[j]);
            exit(0);
        }
    }
}

int main(int argc, char *argv[])
{

    if (argc < 2)
    {
        printf("Error: no parameter file is specified.\n");
        exit(0);
    }
    else if (argc > 2)
    {
        printf("Error: invalid argument '%s'.\n", argv[2]);
        exit(0);
    }

    readParameters(argv[1]);
    
    return 0;
}
