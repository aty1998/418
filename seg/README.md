# Parallel Image Segmentation using Boruvka's MST Algorithm

The program `seg-omp` implements an OpenMP image segmenter that computes minimum
spanning forests with Boruvka's algorithm. The script `segment.py` is a wrapper
around `seg-omp` that converts images into a format the segmenter can use. This
is my class project submission for CMU's 15-418 course.

## Using seg-omp and segment.py

```bash
./seg-omp [-h] [-t nthread] [-c credit] < input.txt > output.txt
```

The `nthread` argument specifies how many OpenMP threads to use. The `credit`
argument is the initial credit for each vertex, used to limit the size of
strongly-connected components in the segmented image.

Due to environment constraints, `seg-omp` only operates with text I/O. The input
should contain R and C, the number of rows and columns in the image, separated
by a space. The next R lines will contain 3C space-separated integers, with
each consecutive triple of integers representing the RGB value of a pixel.

The output will contain R and C separated by a space. Given S is the number of
segments, the next R lines will contain C integers from 0 to S - 1 inclusive,
specifying which segment each pixel belongs to.

The script also prints instrumentation data to stderr.

```bash
python3 segment.py [-h] [-t threads] [-c credit] ifile
```

The `ifile` argument specifies the input image path. The script converts the
image into text format, feeds it into `seg-omp`, and converts the output to
another image. Note that this script will create directories `text` and `images`
in the working directory if they do not already exist.