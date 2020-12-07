# Rendered Examples

This directory has examples of the plotting scripts for EDA using toy data.
The toy data were created through commands like:
```python
import pandas as pd
import numpy as np

## Create a sample data frame with 450 values
d1 = pd.DataFrame(np.arange(1,451))
## Name that column "Residues"
d1.columns = ['Residues']
## Add a normalized (Gaussian) dataset centered on 0 with StDev 0.5
d1['DiffE'] = np.random.normal(loc=0,scale=0.5,size=450)
## Create a uniform dataset for StDev from 0 to 0.1
d1['AvgSTDEV'] = np.random.uniform(low=0., high=0.1, size=450)

## Save a file for Gnuplot examples
with open("toy_data.dat", 'w+') as f:
    f.write("  Residue   TotIntE    AvgStDev\n")
    np.savetxt(f, d1.values, fmt='%10.5f')
```

The file names match the script that generated them with the prefix `test-`.

## Result of `EDA-bar-comp.py`
<img src="https://raw.github.com/emleddin/research-scripts/main/amber-analysis/EDA/rendered-examples/test-EDA-bar-comp.png?raw=true" alt="A bar plot comparing two systems with average standard deviation" width="500"/>

## Result of `EDA-diffs.py`
<img src="https://raw.github.com/emleddin/research-scripts/main/amber-analysis/EDA/rendered-examples/test-EDA-diffs.png?raw=true" alt="A figure with 5 subplots of energy and average standard deviation" width="500"/>

## Result of `EDA-single-diff.py`
<img src="https://raw.github.com/emleddin/research-scripts/main/amber-analysis/EDA/rendered-examples/test-EDA-single-diff.png?raw=true" alt="A bar plot with average standard deviation" width="500"/>

## Result of `EDA-candlestick.gnu`
<img src="https://raw.github.com/emleddin/research-scripts/main/amber-analysis/EDA/rendered-examples/test-EDA-candlestick.png?raw=true" alt="A bar and whisker plot with average standard deviation" width="500"/>
