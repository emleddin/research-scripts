# Normal Mode Analysis (NMA)
These scripts pertain to normal mode analysis (related to principle component
analysis) of protein complexes.


## `normalmodeplot.gnu`
A gnuplot script used to generate graphs of the first three normal modes for a
specific system.

## `eigenplots.py`
A python script to plot the contributions of the different normal modes and
determine which are important for the overall motion of the protein.
This information can then be used to determine the number of modes to plot
using `NMA_plot_mult.py`.

## `NMA_plot_mult.py`
A python script to plot the most relevant normal modes from the `.nmd` file
generated using `cpptraj` or `VMD`.
When using the `.nmd` from `cpptraj`, this replaces `normalmodeplot.gnu`.
