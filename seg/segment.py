from PIL import Image
import numpy as np
import subprocess
import argparse
import os

parser = argparse.ArgumentParser(description="Segments an image.")
parser.add_argument("ifile", help="input file")
parser.add_argument("-c", "--credit", type=int, help="initial credit")
parser.add_argument("-t", "--threads", type=int, help="number of threads")
args = parser.parse_args()

if not os.path.exists("text"):
    os.mkdir("text")
if not os.path.exists("images"):
    os.mkdir("images")

icred = args.credit or 1000000000
nthread = args.threads or 1

name, _ = os.path.splitext(os.path.basename(args.ifile))

# Open input image and convert to RGB array
im = Image.open(args.ifile)
im.load()
data = np.asarray(im, dtype="int32")

# Print RGB array as text file
itpath = os.path.join("text", "%s-in.txt" % name)
with open(itpath, "w") as f:
    f.write("%d %d\n" % (data.shape[0], data.shape[1]))
    for r in range(data.shape[0]):
        f.write(" ".join(map(str, np.ravel(data[r]))))
        f.write("\n")

# Run segmenter program on text file
otpath = os.path.join("text", "%s-out.txt" % name)
cmd = "./seg-omp -t %d -c %d < %s > %s" % (nthread, icred, itpath, otpath)
subprocess.call(cmd, shell=True)

# Read in segmenter output from text file
with open(otpath, "r") as f:
    line = f.readline()
    nrow, ncol = map(int, line.split(" "))
    data = np.zeros((nrow, ncol))
    for r in range(nrow):
        line = f.readline()
        data[r] = np.array(list(map(int, line.split(" "))))

# Randomize segment colors and apply them
segs = np.max(data) + 1
R = np.random.randint(0, 256, (int(segs), 3))
A = np.zeros((nrow, ncol, 3))
for r in range(nrow):
    for c in range(ncol):
        s = int(data[r, c])
        A[r, c, :] = R[s, :]

# Save segmented image to file
im = Image.fromarray(A.astype(np.uint8))
im.save("images/%s-out.png" % name)