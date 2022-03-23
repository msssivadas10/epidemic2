#ifndef MT19937_H
#define MT19937_H 1

#include <stdio.h>
#include <stdlib.h>

#define size_t unsigned long

/* MT19937 parameters  */
#define MT19937_W 32L
#define MT19937_N 624L
#define MT19937_M 397L
#define MT19937_R 31L
#define MT19937_A 0x9908b0dfL
#define MT19937_U 11L
#define MT19937_D 0xffffffffL
#define MT19937_S 7L
#define MT19937_B 0x9d2c5680L
#define MT19937_T 15L
#define MT19937_C 0xefc60000L
#define MT19937_L 18L
#define MT19937_F 1812433253L

#define LOWER_W(x) ((x) & (0xffffffffL)) // get the lower w bits of x
#define MT19937_LM ((1L << MT19937_R) - 1L)
#define MT19937_UM LOWER_W(~MT19937_LM)

extern size_t* MT19937_mt;
extern size_t __index;

/* Initialise the MT19937 random number generator from a seed. */
void initGenerator(size_t seed);

/* Generate a random integer in the range [0, 2^32-1]. */
size_t randi(void);

/* Generate a random number between 0 and 1 (exclusive) */
float randf(void);

#endif