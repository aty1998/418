/**
 * @file main.cpp
 * @brief Runs an image segmenter on std input and prints to std output.
 * @author Andrew Yang (atyang)
 */

#include "mst.hpp"

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include <map>

using namespace std;

/** RGB pixel type. */
struct Pixel
{
    int r, g, b;
};

/** Computes the L1 distance between two pixels. */
int pcmp(Pixel a, Pixel b)
{
    return abs(a.r - b.r) + abs(a.g - b.g) + abs(a.b - b.b);
}

/** Prints usage information. */
void usage(char* argv[])
{
    cerr << "Usage: " << argv[0];
    cerr << " [-t threads] [-c credits]" << endl;
}

int main(int argc, char *argv[])
{
    // Read command-line arguments
    int opt;
    int nthread = 1;
    int icred = 10000000;
    while ((opt = getopt(argc, argv, "ht:c:")) != -1)
    {
        switch (opt)
        {
            case 't':
                nthread = atoi(optarg);
                break;

            case 'c':
                icred = atoi(optarg);
                break;

            case 'h':
                usage(argv);
                exit(EXIT_SUCCESS);
            
            default:
                usage(argv);
                exit(EXIT_FAILURE);
        }
    }

    omp_set_num_threads(nthread);

    // Read in image data
    int nrow, ncol;
    cin >> nrow >> ncol;
    
    Pixel **data = new Pixel*[nrow];
    for (int r = 0; r < nrow; r++)
        data[r] = new Pixel[ncol];
    
    for (int r = 0; r < nrow; r++)
        for (int c = 0; c < ncol; c++)
            cin >> data[r][c].r >> data[r][c].g >> data[r][c].b;

    // Segment image
    SegOMP<Pixel> seg(nrow, ncol, data, pcmp, nthread, icred);
    seg.run();

    // Print segmented image
    cout << nrow << " " << ncol << endl;
    for (int r = 0; r < nrow; r++)
    {
        for (int c = 0; c < ncol; c++)
        {
            if (c > 0)
                cout << " ";
            cout << seg.label[r][c];
        }
        cout << endl;
    }

    // Print instrumentation information
    int nlabels = 5;
    string labels[] = {"joiner", "star", "contract", "prune", "segment"};
    for (int i = 0; i < nlabels; i++)
    {
        long long total = 0;
        vector<long long> times = seg.timer.times[labels[i]];
        for (size_t j = 0; j < times.size(); j++)
        {
            total += times[j];
        }
        cerr << labels[i] << " " << total << endl;
    }
    cerr << "total " << seg.timer.total << endl;
}