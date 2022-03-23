#include "mt19937.h"

size_t* MT19937_mt;
size_t __index = MT19937_N + 1L;

void initGenerator(size_t seed)
{
    MT19937_mt = (size_t*) malloc(MT19937_N * sizeof(size_t));
    
    MT19937_mt[0] = seed;
    __index       = MT19937_N;

    for (size_t i = 1; i < MT19937_N; i++)
    {
        MT19937_mt[i] = LOWER_W( MT19937_F * ( MT19937_mt[i-1] ^ ( MT19937_mt[i-1] >> (MT19937_W - 2L) ) ) + i );
    }
}

void __twist(void)
{
    size_t x, xA;

    for (size_t i = 0; i < MT19937_N; i++)
    {
        x  = ( MT19937_mt[i] & MT19937_UM ) + ( MT19937_mt[ (i+1) % MT19937_N ] & MT19937_LM );
        xA = x >> 1;
        if (x % 2 != 0L) 
            xA = xA ^ MT19937_A;
        
        MT19937_mt[i] = MT19937_mt[ (i + MT19937_M) % MT19937_N ] ^ xA;
    }
    __index = 0;
}

size_t __extract_number(void)
{
    size_t y;

    if (__index >= MT19937_N)
        __twist();

    y = MT19937_mt[__index];
    y = y ^ ( (y >> MT19937_U) & MT19937_D );
    y = y ^ ( (y << MT19937_S) & MT19937_B );
    y = y ^ ( (y << MT19937_T) & MT19937_C );
    y = y ^ (y >> MT19937_L);

    __index++;
    return LOWER_W(y);
}

size_t randi(void)
{
    return __extract_number();
}

float randf(void)
{
    return (float) randi() / (float) (1L << 32);
}
