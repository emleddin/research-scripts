#!/bin/bash

nvt=5Y2S_OC_OCO_db_clustnum_v_time.dat

grep " 0" $nvt > clust_num_0.txt
grep " 1" $nvt > clust_num_1.txt
grep " 2" $nvt > clust_num_2.txt

