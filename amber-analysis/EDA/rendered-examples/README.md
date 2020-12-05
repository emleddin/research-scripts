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
```

The file names match the script that generated them with the prefix `test-`.
