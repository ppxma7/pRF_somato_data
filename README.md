# pRF_somato_data

Processed maps, .mat files for population receptive field mapping performed using 7T-fMRI, in November 2021.

All data are processed and analysed in MATLAB 2019b.

The **data** folder contains processed `.mat` files. For each subject, there is a *baa* folder which has split up the pRF ROI, per model. The *basetips* folder contains the pRF fit maps, i.e. what the pRF fit outputs. Those files with *_big* include additional Brodmann areas.

To explode the pRFs, use `grabData.m`

To remake the figures as seen in the manuscript, use the scripts with the prefix `prfsummary_`

The **pRF_somato** folder contains the fitting code implemented in MATLAB and requires the mrTools software suite: [https://github.com/ppxma7/mrTools]()

