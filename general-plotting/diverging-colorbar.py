import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

## Create a Bottom - Middle - Top gradient
## Ex: Bottom (-1), Middle (0), Top (1)
## This can and will also be reversed!
## Ex: Purple - Orange - White gradient
##      Light var: Purple: #998EC3, White: #F7F7F7, Orange: #F1A340
##          153,142,195; 247,247,247; 241,163,64
##      Dark var: Purple: #5E3C99, White: #F7F7F7, Orange: #E66101
##          94,60,153; 247,247,247; 230,97,1
Bot_RGB = [ 94,  60, 153] # Dark Purple
Mid_RGB = [247, 247, 247] # White
Top_RGB = [230,  97,   1] # Dark Orange

## Names for Colorbars
## Use a string
Bot_to_Top_Name = "PurpleOrange"
Top_to_Bot_Name = "OrangePurple"

def NewCmap(Bot_RGB, Mid_RGB, Top_RGB, Bot_to_Top_Name, Top_to_Bot_Name):
    """Builds a new diverging colormap.
    Parameters
    ----------
    Bot_RGB : list (3 items)
        The RGB values for the bottom-most color for the colorbar.
        Black would be [0,0,0].
    Mid_RGB : list (3 items)
        The middle color, where things diverge.
    Top_RGB : list (3 items)
        The top-most color for the colorbar.
    Bot_to_Top_Name : str
        How to name the colorbar from bottom to top in the mpl map
    Top_to_Bot_Name : str
        How to name the colorbar from top to bottom in the mpl map
    Returns
    -------
    newcmap : matplotlib.colors.ListedColormap
        A diverging colormap built using the input colors.
    """
    ## Set up RGB colorspace
    N = 256
    ## Build the Bottom to Mid Gradient
    bvals = np.ones((N, 4))
    bvals[:, 0] = np.linspace(Bot_RGB[0]/N, Mid_RGB[0]/N, N)
    bvals[:, 1] = np.linspace(Bot_RGB[1]/N, Mid_RGB[1]/N, N)
    bvals[:, 2] = np.linspace(Bot_RGB[2]/N, Mid_RGB[2]/N, N)
    bottom_color = ListedColormap(bvals)
    #
    ## Build the Mid to Top Gradient
    tvals = np.ones((N, 4))
    tvals[:, 0] = np.linspace(Mid_RGB[0]/N, Top_RGB[0]/N, N)
    tvals[:, 1] = np.linspace(Mid_RGB[1]/N, Top_RGB[1]/N, N)
    tvals[:, 2] = np.linspace(Mid_RGB[2]/N, Top_RGB[2]/N, N)
    top_color = ListedColormap(tvals)
    #
    ## Save the Bottom - Top gradient
    newcolors = np.vstack((bottom_color(np.linspace(0,1,128)),
                              top_color(np.linspace(0,1,128))))
    newcmap = ListedColormap(newcolors, name=Bot_to_Top_Name)
    mpl.cm.register_cmap(cmap=newcmap)
    ## Save the reverse, Top - Bottom gradient
    revnewcolors = np.vstack((top_color(np.linspace(1,0,128)),
                              bottom_color(np.linspace(1,0,128))))
    rev_newcmap = ListedColormap(revnewcolors, name=Top_to_Bot_Name)
    mpl.cm.register_cmap(cmap=rev_newcmap)

# https://stackoverflow.com/questions/16595138/standalone-colorbar-matplotlib
def PrintColbar(outfile, orient, plot_cmap, bar_min=-1.0, bar_max=1.0):
    """Saves a figure with just the colorbar.
    Parameters
    ----------
    outfile : str
        The save name for the figure.
    orient : str
        Orientation of the colorbar, either 'horizontal' or 'vertical'.
    plot_cmap : str
        Name of the colormap to use for the colorbar.
    bar_min : float
        Minimum value for the colorbar.
    bar_max : float
        Maximum value for the colorbar.
    """
    ## Save a picture with the colorbar!
    fig = plt.figure()
    ## left, bottom, width, height
    if orient.lower() == 'vertical':
        ax = fig.add_axes([0.05, 0.80, 0.1, 0.9])
    elif orient.lower() == 'horizontal':
        ax = fig.add_axes([0.05, 0.80, 0.9, 0.1])
    #
    use_cmap = mpl.cm.get_cmap(plot_cmap)
    norm = mpl.colors.Normalize(vmin=bar_min, vmax=bar_max)
    cb = mpl.colorbar.ColorbarBase(ax, norm=norm, orientation=orient,
         cmap=use_cmap)
    #
    plt.savefig(outfile, bbox_inches='tight', dpi=1000)

#------------------------------ Run the code! --------------------------------#
## Set-up the colormap
NewCmap(Bot_RGB, Mid_RGB, Top_RGB, Bot_to_Top_Name, Top_to_Bot_Name)

## Save a figure of the forward
PrintColbar(outfile="cbar_"+Bot_to_Top_Name+".png", bar_min=-1.0, bar_max=1.0,
            orient="horizontal", plot_cmap=Bot_to_Top_Name)

## Save a figure of the reverse
PrintColbar(outfile="cbar_"+Top_to_Bot_Name+".png", bar_min=-1.0, bar_max=1.0,
            orient="vertical", plot_cmap=Top_to_Bot_Name)
