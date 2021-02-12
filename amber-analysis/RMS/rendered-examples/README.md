# Rendered Examples

This directory has examples of the plotting scripts for RMS using toy data.
The toy RMSD data were created using:
```python
import pandas as pd
import numpy as np

## Set number of frames for plotting (50000/500 is 100 frames)
num_frames = 50000
## Set Gaussian center for the first dataset
start_rmsd  = 3
## Set number for StdDev
my_stdev = 0.5
## Use an increment to add to center so lines aren't overlapping for additional
## datasets
increment  = 0.25
## Choose how many datasets to create
datasets   = 4

for i in range(datasets):
    ## Create a sample data frame with num_frames values
    d1 = pd.DataFrame(np.arange(1,num_frames+1))
    ## Name that column "Frame"
    d1.columns = ['Frame']
    ## Add a normalized (Gaussian) dataset centered on 3 with StDev my_stdev
    ## Add increment * interation (first is 0 because Python)
    d1['RMS'] = np.random.normal(loc=(start_rmsd+(increment*i)),scale=my_stdev,
     size=num_frames)
    ## Create a uniform dataset for StDev from 0 to 0.1
    d1['AvgSTDEV'] = np.random.uniform(low=0., high=0.1, size=num_frames)

    ## Save a file for Gnuplot examples
    with open(f'toy_rmsd_data_{i+1}.dat', 'w+') as f:
        f.write("  Frame     RMSD    AvgStDev\n")
        np.savetxt(f, d1.values, fmt='%10.5f')
```

The toy total hydrogen bond data were created using:
```python
import pandas as pd
import numpy as np

## Set number of frames for plotting (50000/500 is 100 frames)
num_frames = 50000
## Set Gaussian center for the first dataset
start_hb  = 185
## Set number for StdDev
my_stdev = 5
## Use an increment to add to center so lines aren't overlapping for additional
## datasets
increment  = 7.5
## Choose how many datasets to create
datasets   = 4

for i in range(datasets):
    ## Create a sample data frame with num_frames values
    d1 = pd.DataFrame(np.arange(1,num_frames+1))
    ## Name that column "Frame"
    d1.columns = ['Frame']
    ## Add a normalized (Gaussian) dataset centered on start_hb w/ StDev my_stdev
    ## Add increment * interation (first is 0 because Python)
    d1['RMS'] = np.random.normal(loc=(start_hb+(increment*i)),scale=my_stdev,
     size=num_frames)
    ## Create a uniform dataset for StDev from 0 to 0.1
    d1['AvgSTDEV'] = np.random.uniform(low=0., high=0.1, size=num_frames)

    ## Save a file for Gnuplot examples
    with open(f'toy_hb_data_{i+1}.dat', 'w+') as f:
        f.write("  Frame     HB     AvgStDev\n")
        np.savetxt(f, d1.values, fmt='%10.5f')
```

The toy RMSF data were created using:
```python
import pandas as pd
import numpy as np

## Set number of residues
num_res = 455
## Set Gaussian center for the first dataset
start_rmsf  = 1.5
## Set number for StdDev
my_stdev = 2
## Choose how many datasets to create
datasets   = 4

for i in range(datasets):
    ## Create a sample data frame with num_frames values
    d1 = pd.DataFrame(np.arange(1,num_res+1))
    ## Name that column "Frame"
    d1.columns = ['Frame']
    ## Add a normalized (Gaussian) dataset centered on start_rmsf w/ StDev my_stdev
    ## Force all to be positive numbers
    d1['RMS'] = np.random.normal(loc=start_rmsf,scale=my_stdev,size=num_res)
    d1['RMS'] = abs(d1['RMS'])
    ## Create a uniform dataset for StDev from 0 to 0.1
    d1['AvgSTDEV'] = np.random.uniform(low=0., high=0.1, size=num_res)

    ## Save a file for Gnuplot examples
    with open(f'toy_rmsf_data_{i+1}.dat', 'w+') as f:
        f.write("  Frame     RMSF    AvgStDev\n")
        np.savetxt(f, d1.values, fmt='%10.5f')
```

The EPS files generated with `gnuplot` can be converted to PNGs using
[ImageMagick](https://imagemagick.org/index.php).
```bash
convert avg-rmsd-etc-rmsd.eps -rotate 90 avg-rmsd-etc-rmsd.png
convert avg-rmsd-etc-hbond.eps -rotate 90 avg-rmsd-etc-hbond.png
convert avg-rmsd-etc-rmsf.eps -rotate 90 avg-rmsd-etc-rmsf.png

convert rmsd-etc-rmsd.eps -rotate 90 rmsd-etc-rmsd.png
convert rmsd-etc-hbond.eps -rotate 90 rmsd-etc-hbond.png
convert rmsd-etc-rmsf.eps -rotate 90 rmsd-etc-rmsf.png
```

The file names match the script that generated them with the prefix `test-` and
the final word reflecting which analysis was plotted in that figure.

## Results of `avg-rmsd-etc.gnu`
<img src="https://raw.github.com/emleddin/research-scripts/main/amber-analysis/RMS/rendered-examples/test-avg-rmsd-etc-rmsd.png?raw=true" alt="A line graph of RMSD with average standard deviation for WT and Mutant C" width="500"/>

<img src="https://raw.github.com/emleddin/research-scripts/main/amber-analysis/RMS/rendered-examples/test-avg-rmsd-etc-hbond.png?raw=true" alt="A line graph of the total number of hydrogen bonds with average standard deviation for WT and Mutant C" width="500"/>

<img src="https://raw.github.com/emleddin/research-scripts/main/amber-analysis/RMS/rendered-examples/test-avg-rmsd-etc-rmsf.png?raw=true" alt="A line graph of RMSF with average standard deviation for WT and Mutant C" width="500"/>

## Results of `moving-averages.py`
<img src="https://raw.github.com/emleddin/research-scripts/main/amber-analysis/RMS/rendered-examples/test-moving-avgs-rmsd.png?raw=true" alt="A line graph with the moving average over 1 ns of RMSD for WT, Mutant A, Mutant B, and Mutant C" width="500"/>

<img src="https://raw.github.com/emleddin/research-scripts/main/amber-analysis/RMS/rendered-examples/test-moving-avgs-hbond.png?raw=true" alt="A line graph with the moving average over 1 ns of the number of hydrogen bonds for WT, Mutant A, Mutant B, and Mutant C" width="500"/>

<img src="https://raw.github.com/emleddin/research-scripts/main/amber-analysis/RMS/rendered-examples/test-moving-avgs-rmsf.png?raw=true" alt="A line graph of RMSF for WT, Mutant A, Mutant B, and Mutant C" width="500"/>

## Results of `rmsd-etc.gnu`
<img src="https://raw.github.com/emleddin/research-scripts/main/amber-analysis/RMS/rendered-examples/test-rmsd-etc-rmsd.png?raw=true" alt="A line graph of RMSD with WT, Mutant A, Mutant B, and Mutant C" width="500"/>

<img src="https://raw.github.com/emleddin/research-scripts/main/amber-analysis/RMS/rendered-examples/test-rmsd-etc-hbond.png?raw=true" alt="A line graph of the total number of hydrogen bonds with WT, Mutant A, Mutant B, and Mutant C" width="500"/>

<img src="https://raw.github.com/emleddin/research-scripts/main/amber-analysis/RMS/rendered-examples/test-rmsd-etc-rmsf.png?raw=true" alt="A line graph of RMSF with WT, Mutant A, Mutant B, and Mutant C" width="500"/>
