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

struct Pixel
{
    int r, g, b;
};

int pcmp(Pixel a, Pixel b)
{
    return abs(a.r - b.r) + abs(a.g - b.g) + abs(a.b - b.b);
}

int main(int argc, char *argv[])
{
    int opt;
    int nthread = 1;
    int icred = 10000000;

    while ((opt = getopt(argc, argv, "o:i:t:c:")) != -1)
    {
        switch (opt)
        {
            case 't':
                nthread = atoi(optarg);
                break;

            case 'c':
                icred = atoi(optarg);
                break;
            
            default:
                cerr << "Usage: " << argv[0] << " [-t threads]" << endl;
                exit(EXIT_FAILURE);
        }
    }

    omp_set_num_threads(nthread);

    int nrow, ncol;
    cin >> nrow >> ncol;
    
    Pixel **data = new Pixel*[nrow];
    for (int r = 0; r < nrow; r++)
        data[r] = new Pixel[ncol];
    
    for (int r = 0; r < nrow; r++)
        for (int c = 0; c < ncol; c++)
            cin >> data[r][c].r >> data[r][c].g >> data[r][c].b;

    SegOMP<Pixel> seg(nrow, ncol, data, pcmp, nthread, icred);
    seg.run();

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

    int nlabels = 4;
    string labels[] = {"joiner", "star", "contract", "prune"};
    
    for (int i = 0; i < nlabels; i++)
    {
        long long total = 0;
        vector<long long> times = seg.timer.times[labels[i]];
        for (size_t j = 0; j < times.size(); j++)
        {
            total += times[j];
        }

        cerr << labels[i] << ":" << endl;
        cerr << "    total (ms): " << total / 1000.0 << endl;
    }

    cerr << seg.timer.total / 1000 << endl;
}