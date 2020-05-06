#ifndef __MST_H__
#define __MST_H__

#include "par.hpp"
#include "instrument.hpp"

#include <omp.h>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <queue>

#define INF (1 << 20)

const int NDIR = 8;
const int DR[] = {0, -1, -1, -1, 0, 1, 1, 1};
const int DC[] = {1, 1, 0, -1, -1, -1, 0, 1};

/** Edge struct. */
struct Edge
{
    int u, v, w, ou, ov;
    bool operator <(const Edge& rhs) const
    {
        return w < rhs.w;
    }
};

/** Class that segments an image using Boruvka's MST algorithm. */
template <typename T>
class SegOMP
{
public:
    Timer timer;
    int nsegs;
    int **label;

    SegOMP(int nrow, int ncol, T **data, int cmp(T, T), int nthread, int icred);

    void run();

private:
    int nthread;
    int nrow, ncol;
    int nvert, nedge;

    int *offsets;
    Edge *edges;

    int *_offsets;
    int **_counts;
    int **_boffsets;
    Edge *_edges;

    int *credits;
    int *_credits;

    int *joiners;
    int *jprev;

    unsigned int *seeds;
    bool *heads;
    int *starmap;

    int *indexmap;
    
    int mstsize;
    Edge *mstedges;

    int cell(int r, int c) const;
    int row(int vid) const;
    int col(int vid) const;
    int bin_start(int size, int ti) const;
    int bin_size(int size) const;

    void find_joiners();
    void find_stars();
    void contract_stars();
    void prune_vertices();

    void segment();
};

using namespace std;

template <typename T>
int SegOMP<T>::cell(int r, int c) const
{
    return r * ncol + c;
}

template <typename T>
int SegOMP<T>::row(int vid) const
{
    return vid / ncol;
}

template <typename T>
int SegOMP<T>::col(int vid) const
{
    return vid % ncol;
}

template <typename T>
int SegOMP<T>::bin_size(int size) const
{
    return (size + nthread - 1) / nthread;
}

template <typename T>
int SegOMP<T>::bin_start(int size, int ti) const
{
    return min(bin_size(size) * ti, size);
}

template <typename T>
SegOMP<T>::SegOMP(int nrow, int ncol, T **data, int cmp(T, T), int nthread,
                  int icred) : timer()
{
    assert(nrow >= 0 && ncol >= 0 && nthread > 0);
    this->nrow = nrow;
    this->ncol = ncol;
    this->nthread = nthread;

    nvert = nrow * ncol;
    int asize = nextpow2(nvert + 1);
    offsets = new int[asize];

    // Determine offsets for adjacency array and total number of edges
    #pragma omp parallel for
    for (int r = 0; r < nrow; r++)
    {
        for (int c = 0; c < ncol; c++)
        {
            int vid = r * ncol + c;
            offsets[vid] = 0;
            for (int d = 0; d < NDIR; d++)
            {
                int nr = r + DR[d];
                int nc = c + DC[d];
                if (nr < 0 || nr >= nrow || nc < 0 || nc >= ncol)
                    continue;
                offsets[vid]++;
            }
        }
    }
    exscan_int(offsets, asize);
    nedge = offsets[nvert];
    edges = new Edge[nedge + 1];

    _offsets = new int[asize];
    _counts = new int*[nthread];
    _boffsets = new int*[nthread];
    for (int ti = 0; ti < nthread; ti++)
    {
        _counts[ti] = new int[asize];
        _boffsets[ti] = new int[asize];
    }
    _edges = new Edge[nedge + 1];

    credits = new int[asize];
    for (int vid = 0; vid < nvert; vid++)
        credits[vid] = icred;
    _credits = new int[asize];

    // Fill in adjacency array
    #pragma omp parallel for
    for (int r = 0; r < nrow; r++)
    {
        for (int c = 0; c < ncol; c++)
        {
            int vid = r * ncol + c;
            int eid = offsets[vid];
            T val = data[r][c];

            for (int d = 0; d < NDIR; d++)
            {
                int nr = r + DR[d];
                int nc = c + DC[d];
                if (nr < 0 || nr >= nrow || nc < 0 || nc >= ncol)
                    continue;
                int diff = cmp(val, data[nr][nc]);
                int ovid = cell(nr, nc);
                edges[eid++] = {vid, ovid, diff, vid, ovid};
            }

            assert(eid == offsets[vid + 1]);
        }
    }

    joiners = new int[nvert];
    jprev = new int[nthread];

    seeds = new unsigned int[nvert];
    for (int vid = 0; vid < nvert; vid++)
    {
        seeds[vid] = vid;
        rand_r(&seeds[vid]);
    }
    heads = new bool[nvert];
    starmap = new int[nvert];

    indexmap = new int[asize];

    mstsize = 0;
    mstedges = new Edge[nvert - 1];
    label = new int*[nrow];
    for (int r = 0; r < nrow; r++)
        label[r] = new int[ncol];
}

template <typename T>
void SegOMP<T>::run()
{
    while (nedge > 0)
    {
        find_joiners();
        find_stars();
        contract_stars();
        prune_vertices();
    }

#if DEBUG
    cout << endl;
    for (int ei = 0; ei < mstsize; ei++)
    {
        Edge e = mstedges[ei];
        cout << "MST edge " << e.ou << " " << e.ov << endl;
    }
#endif

    segment();
}

template <typename T>
void SegOMP<T>::find_joiners()
{
    timer.start("joiner");

    // Initialize joiners to infinite-weight edge
    edges[nedge].w = INF;
    #pragma omp parallel for
    for (int vid = 0; vid < nvert; vid++)
        joiners[vid] = nedge;

    // Populate joiners array by iterating through edges
    #pragma omp parallel for
    for (int ti = 0; ti < nthread; ti++)
    {
        // Start and end EIDs in this bin
        int bs = bin_start(nedge, ti);
        int be = bin_start(nedge, ti + 1);

        int us = edges[bs].u;
        bool conflict = bs > offsets[us];

        int u = us;
        jprev[ti] = nedge;

        for (int eid = bs; eid < be; eid++)
        {
            // Move onto next source vertex if done with neighbors
            if (eid == offsets[u + 1])
                u++;

            // Update joiners array if edge has lower weight
            if (u == us && conflict)
            {
                if (edges[eid] < edges[jprev[ti]])
                    jprev[ti] = eid;
            }
            else
            {
                if (edges[eid] < edges[joiners[u]])
                    joiners[u] = eid;
            }
        }
    }

    // Resolve conflicts by performing reductions for each vertex
    int bsize = bin_size(nedge);
    #pragma omp parallel for
    for (int vid = 0; vid < nvert; vid++)
    {
        int ts = offsets[vid] / bsize;
        int te = (offsets[vid + 1] - 1) / bsize;

        for (int ti = ts + 1; ti <= te; ti++)
        {
            if (edges[jprev[ti]] < edges[joiners[vid]])
                joiners[vid] = jprev[ti];
        }
    }
    
    timer.stop("joiner");

#if DEBUG
    for (int vid = 0; vid < nvert; vid++)
    {
        int j = joiners[vid];
        Edge e = edges[j];
        cout << "vid joiner: " << vid << " " << e.v << endl;
    }
#endif
}

template <typename T>
void SegOMP<T>::find_stars()
{
    timer.start("star");

    // Flip coins to decide star centers
    #pragma omp parallel for
    for (int vid = 0; vid < nvert; vid++)
        heads[vid] = rand_r(&seeds[vid]) % 2;

    // Mark (T, H) edges as satellites and map to star centers
    #pragma omp parallel for
    for (int vid = 0; vid < nvert; vid++)
    {
        int eid = joiners[vid];
        int v = edges[eid].v;
        bool satedge = !heads[vid] && heads[v];
        starmap[vid] = satedge ? v : vid;
        indexmap[vid] = satedge; 
    }
    
    // Compute prefix sum of satellite edge mark array
    exscan_int(indexmap, nextpow2(nvert + 1));

    // Add satellite edges into MST edges
    #pragma omp parallel for
    for (int vid = 0; vid < nvert; vid++)
    {
        if (starmap[vid] != vid)
            mstedges[indexmap[vid] + mstsize] = edges[joiners[vid]];
    }
    mstsize += indexmap[nvert];

    // Update vertex credits
    #pragma omp parallel for
    for (int vid = 0; vid < nvert; vid++)
    {
        int start = offsets[vid];
        int end = offsets[vid + 1];
        
        // New credits = min credits in star - total weight of star
        int min_credit = credits[vid];
        int weight_sum = 0;
        for (int eid = start; eid < end; eid++)
        {
            Edge e = edges[eid];
            if (heads[vid] && !heads[e.v])
            {
                weight_sum += e.w;
                if (credits[e.v] < min_credit)
                    min_credit = credits[e.v];
            }
        }
        credits[vid] = min_credit - weight_sum;
    }

    timer.stop("star");

#if DEBUG
    for (int vid = 0; vid < nvert; vid++)
    {
        if (starmap[vid] != vid)
        {
            int eid = joiners[vid];
            Edge e = edges[eid];
            cout << "MST edge " << e.ou << " " << e.ov << endl;
        }
    }

    for (int vid = 0; vid < nvert; vid++)
    {
        cout << "vid star: " << vid << " " << starmap[vid] << endl;
    }
#endif
}

template <typename T>
void SegOMP<T>::contract_stars()
{
    timer.start("contract");

    // Vertices that will remain after contraction
    #pragma omp parallel for
    for (int vid = 0; vid < nvert; vid++)
        indexmap[vid] = starmap[vid] == vid;
    
    // Map old VIDs to new VIDs
    exscan_int(indexmap, nextpow2(nvert + 1));
    int _nvert = indexmap[nvert];
    int asize = nextpow2(_nvert + 1);

    // Compute vertex degree counts within each bin
    #pragma omp parallel for
    for (int ti = 0; ti < nthread; ti++)
    {
        // Zero-out vertex degrees
        for (int nvid = 0; nvid < _nvert; nvid++)
            _counts[ti][nvid] = 0;

        int bs = bin_start(nedge, ti);
        int be = bin_start(nedge, ti + 1);

        // Add edges to appropriate degree counts
        for (int eid = bs; eid < be; eid++)
        {
            Edge e = edges[eid];
            int pu = starmap[e.u];
            int pv = starmap[e.v];
            if (pu != pv && e.w <= credits[pu] && e.w <= credits[pv])
                _counts[ti][indexmap[pu]]++;
        }
    }

    // Compute offsets within each bin and collect total vertex counts
    #pragma omp parallel for
    for (int nvid = 0; nvid < _nvert; nvid++)
    {
        int total = 0;
        for (int ti = 0; ti < nthread; ti++)
        {
            _boffsets[ti][nvid] = total;
            total += _counts[ti][nvid];
        }
        _offsets[nvid] = total;
    }

    // Compute new vertex offsets
    exscan_int(_offsets, asize);
    int _nedge = _offsets[_nvert];

    // Place edges into new places
    #pragma omp parallel for
    for (int ti = 0; ti < nthread; ti++)
    {
        // Zero-out vertex degrees
        for (int nvid = 0; nvid < _nvert; nvid++)
            _counts[ti][nvid] = 0;

        int bs = bin_start(nedge, ti);
        int be = bin_start(nedge, ti + 1);

        // Add edges to appropriate degree counts
        for (int eid = bs; eid < be; eid++)
        {
            Edge e = edges[eid];
            int pu = starmap[e.u];
            int pv = starmap[e.v];
            if (pu != pv && e.w <= credits[pu] && e.w <= credits[pv])
            {
                int nu = indexmap[pu];
                int nv = indexmap[pv];

                int off = _offsets[nu];
                int boff = _boffsets[ti][nu];
                int ei = _counts[ti][nu];
                _edges[off + boff + ei] = {nu, nv, e.w, e.ou, e.ov};
                _counts[ti][nu]++;
            }
        }
    }

    // Place credits into new positions
    #pragma omp parallel for
    for (int vid = 0; vid < nvert; vid++)
    {
        if (starmap[vid] == vid)
            _credits[indexmap[vid]] = credits[vid];
    }
    
    // Copy new credits into credits array
    #pragma omp parallel for
    for (int nvid = 0; nvid < _nvert; nvid++)
        credits[nvid] = _credits[nvid];

    // Copy new offsets into offsets array
    #pragma omp parallel for
    for (int nvid = 0; nvid <= _nvert; nvid++)
        offsets[nvid] = _offsets[nvid];
    nvert = _nvert;

    // Copy new edges into edges array
    #pragma omp parallel for
    for (int neid = 0; neid < _nedge; neid++)
        edges[neid] = _edges[neid];
    nedge = _nedge;
    
    timer.stop("contract");

#if DEBUG
    for (int eid = 0; eid < nedge; eid++)
    {
        int u = edges[eid].u;
        int v = edges[eid].v;
        int ou = edges[eid].ou;
        int ov = edges[eid].ov;
        int w = edges[eid].w;
        cout << "new edge: " << ou << " " << ov << " " << w;
        cout << " " << u << " " << v << endl; 
    }
#endif
}

template <typename T>
void SegOMP<T>::prune_vertices()
{
    timer.start("prune");

    // Vertices that have degree > 0
    #pragma omp parallel for
    for (int vid = 0; vid < nvert; vid++)
        indexmap[vid] = offsets[vid] != offsets[vid + 1];
    
    // Map old VIDs to new VIDs
    exscan_int(indexmap, nextpow2(nvert + 1));
    int _nvert = indexmap[nvert];
    
    // Remap endpoints of each edge
    #pragma omp parallel for
    for (int eid = 0; eid < nedge; eid++)
    {
        edges[eid].u = indexmap[edges[eid].u];
        edges[eid].v = indexmap[edges[eid].v];
    }

    // Remap credits and offsets
    #pragma omp parallel for
    for (int vid = 0; vid < nvert; vid++)
    {
        if (offsets[vid] != offsets[vid + 1])
        {
            _credits[indexmap[vid]] = credits[vid];
            _offsets[indexmap[vid]] = offsets[vid];
        }
    }

    // Copy new credits into credits array
    #pragma omp parallel for
    for (int nvid = 0; nvid < _nvert; nvid++)
        credits[nvid] = _credits[nvid];

    // Copy new offsets into offsets array
    #pragma omp parallel for
    for (int nvid = 0; nvid <= _nvert; nvid++)
        offsets[nvid] = _offsets[nvid];
    nvert = _nvert;
    
    timer.stop("prune");
}

template <typename T>
void SegOMP<T>::segment()
{
    timer.start("segment");

    int n = nrow * ncol;
    vector<vector<int> > AL(n);

    for (int ei = 0; ei < mstsize; ei++)
    {
        Edge e = mstedges[ei];
        int u = e.ou;
        int v = e.ov;
        AL[u].push_back(v);
        AL[v].push_back(u);
    }

    for (int vid = 0; vid < n; vid++)
        label[row(vid)][col(vid)] = -1;

    nsegs = 0;
    queue<int> Q;
    for (int src = 0; src < n; src++)
    {
        if (label[row(src)][col(src)] != -1)
            continue;
        
        Q.push(src);

        while (!Q.empty())
        {
            int u = Q.front();
            Q.pop();

            if (label[row(u)][col(u)] != -1)
                continue;
            label[row(u)][col(u)] = nsegs;

            for (int i = 0; i < AL[u].size(); i++)
            {
                int v = AL[u][i];
                if (label[row(v)][col(v)] == -1)
                    Q.push(v);
            }
        }

        nsegs++;
    }

    timer.stop("segment");

#if DEBUG
    for (int r = 0; r < nrow; r++)
    {
        for (int c = 0; c < ncol; c++)
        {
            cout << label[cell(r, c)];
        }
        cout << endl;
    }
#endif
}

#endif /* __MST_H__ */