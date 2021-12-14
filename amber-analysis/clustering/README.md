# Clustering

These scripts pertain to clustering (typically for QM/MM).
See the [Comp Chem Website](https://emleddin.github.io/comp-chem-website/Analysisguide-clustering.html)
and the [LICHEM Tutorial](https://emleddin.github.io/LICHEM-tutorial/03-MD-cluster/index.html)
for further information.

## `original_clust.in`
This is an example cpptraj script for an initial pass at clustering.

## `frames_from_clusters.sh`
This script can be used to write a list of frames associated with an individual
cluster.

## `clust_framelist.py`
This will read in the file created using `frames_from_clusters.sh` and write
out a comma-separated list of frames to include in `subclust.in`.

## `subclust.in`
This is an example cpptraj script for a second round of clustering from a
cluster with parameters most like what you're seeking ("subclustering.")

## `write-from-subclust.in`
This is an example cpptraj script for writing out the individual frames
identified through subclustering.
