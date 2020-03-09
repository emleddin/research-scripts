import numpy as np
import prody as prd
import matplotlib.pyplot as plt

## Typically 3 modes will be enough
num_of_modes = 4

## Create a list of tuple (infile_name, outfile_name, system_name)
## These are the NMD file, the PNG file, and the tag for the system to determine
## the top X ticks through plot_ticks
in_out_sys_names = [
  ("WT_protein_system_r1_100.nmd", "WT_protein_system_r1_NMA.png", "WT"),
  ("WT_protein_system_r2_100.nmd", "WT_protein_system_r2_NMA.png", "WT"),
  ("WT_protein_system_r3_100.nmd", "WT_protein_system_r3_NMA.png", "WT"),
  ("MUT_A_system_r1_100.nmd", "MUT_A_system_r1_NMA.png", "MUTA"),
  ("MUT_A_system_r2_100.nmd", "MUT_A_system_r2_NMA.png", "MUTA"),
  ("MUT_A_system_r3_100.nmd", "MUT_A_system_r3_NMA.png", "MUTA"),
  ("MUT_B_system_r1_100.nmd", "MUT_B_system_r1_NMA.png", "MUTB"),
  ("MUT_B_system_r2_100.nmd", "MUT_B_system_r2_NMA.png", "MUTB"),
  ("MUT_B_system_r3_100.nmd", "MUT_B_system_r3_NMA.png", "MUTB"),
]

def plot_ticks(sys, NMA_data):
    """Sets top xticks. You NEED the 0 and NMA_data.numAtoms(), otherwise the
    scale will be turned off.
    **This is an example, you'll need to modify it for your system.**
    """
    if sys == "MUTA":
        labels_top = ["", "MUTA", "GS linker", "", "DNA", ""]
        places_top = [0, 141, 334, 346, 431, NMA_data.numAtoms()]
    elif sys == "MUTB":
        labels_top = ["", "GS linker", "", "MUTB", "DNA", ""]
        places_top = [0, 334, 346, 378, 431, NMA_data.numAtoms()]
    elif sys == "WT":
        labels_top = ["","GS linker", "", "DNA", ""]
        places_top = [0, 334, 346, 431, NMA_data.numAtoms()]
    else:
        labels_top = []
        places_top = []
    return labels_top, places_top


def NMA_plots(filename,outfile,sys):
        """Creates a plot of the most important modes for a system.
        Parameters
        ----------
        filename : str
            An NMD file.
        outfile: str
            Name of the output PNG.
        sys: str
            A name (e.g., WT, MUT A, MUT B, etc.) for the system. Replicates
            should have the same name.
        """
        NMA_data,Atom_Group = prd.parseNMD(filename)
        eigens = NMA_data.getEigvals()

        labels_top, places_top = plot_ticks(sys, NMA_data)

        scales=[]
        temp = open(filename)
        lines = temp.readlines()
        temp.close()
        for line in lines:
            if 'mode' in line[:5]:
                scales.append(float(line.split()[:3][-1]))

        ## Make an array of the number of atoms for plotting
        x_vals = np.arange(0, NMA_data.numAtoms(), 1)

        fig = plt.figure(figsize=(10,8),dpi=300)
        ax = fig.add_subplot(1,1,1)
        for i in range(num_of_modes):
            dataset = [np.linalg.norm(NMA_data.getEigvecs()[:,i][n:n+3])*scales[i]*eigens[i] for n in range(0, NMA_data.numEntries(), 3)]
            ax.bar(x_vals, dataset, width=1.0, label="Mode "+str(i+1))

        ax_top = ax.twiny()
        ax_top.set_xticks(places_top)
        ax_top.set_xticklabels(labels_top, fontdict=None, minor=False)

        ax.legend()
        ax.set_xlabel("Residue Number")
        ax.set_ylabel("PCA Square Fluctuations")
        ax.set_xlim([0,x_vals.size]) ## Remove white space at edge
        plt.tight_layout()
        fig.savefig(outfile,dpi=300)
        plt.close()

for filename,outfile,sys in in_out_sys_names:
    NMA_plots(filename,outfile,sys)
