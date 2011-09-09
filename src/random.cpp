/* **********************************
    Random Number Generator

    File: random.cpp
   ********************************** */

#include <math.h>
#include <time.h>

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

long long int idum;

double ran1(long long int *idum)
{
    int j;
    long long int k;
    static long long int iy=0;
    static long long int iv[NTAB];
    long double temp;
    
    if (*idum <= 0 || !iy) {
        if (-(*idum) < 1) 
            *idum = 1;
        else 
            *idum = -(*idum);
        for (j = NTAB + 7; j >= 0; j--) {
            k = (*idum)/IQ;
            *idum = IA*(*idum - k*IQ) - IR*k;
            if (*idum < 0) 
                *idum += IM;
            if (j < NTAB) 
                iv[j] = *idum;
        }
        iy=iv[0];
    }
    k = (*idum)/IQ;
    *idum = IA*(*idum - k*IQ) - IR*k;
    if (*idum < 0) 
        *idum += IM;
    j = iy/NDIV;
    iy = iv[j];
    iv[j] = *idum;
    temp = AM*iy;
    if (temp > RNMX)
        return RNMX;
    else 
        return temp;
}

void seed_ran(long long int seed) {
    idum = seed;
    ran1(&idum);
}

double ran() {
    return ran1(&idum);
}

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
