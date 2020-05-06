#include "par.h"

/** Rounds n up to a power of 2. */
int nextpow2(int n)
{
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    n++;
    return n;
}

/**
 * Computes the exclusive prefix sum of an array in-place.
 * @param data the data array
 * @param n the length of the array
 * @pre n is a power of 2
 * @remark code taken from asst2.pdf
 */
void exscan_int(int *data, int n)
{
    int twod, twod1, i;

    // Upsweep phase
    for (twod = 1; twod < n; twod <<= 1)
    {
        twod1 = twod << 1;
        #pragma omp parallel for
        for (i = 0; i < n; i += twod1)
            data[i + twod1 - 1] += data[i + twod - 1];
    }
    data[n - 1] = 0;

    // Downsweep phase
    for (twod = n >> 1; twod >= 1; twod >>= 1)
    {
        twod1 = twod << 1;
        #pragma omp parallel for
        for (i = 0; i < n; i += twod1)
        {
            int t = data[i + twod - 1];
            data[i + twod - 1] = data[i + twod1 - 1];
            data[i + twod1 - 1] += t;
        }
    }
}