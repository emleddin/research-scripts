import pandas as pd

in_cluster = "clust_num_0.txt"
out_cluster = "out_clust0.txt"

## Read the infile and write comma separated list
clust = pd.read_csv(in_cluster, header=None, delim_whitespace=True)
test = clust.loc[:,0].values.tolist()
test2 = ''.join(str(i)+"," for i in test)

## Write the outfile
f=open(out_cluster,"w+")
f.write(test2)
f.close()

