#! /usr/bin/env python
## Rename files ending with "Skim.root"

import os
import sys
import subprocess
import numpy as np
import ROOT as R

R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn

## Location of postprocessed input files
IN_DIR = '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/data/PNet_v2_2024_11_22/'

## Loop over datasets
#print(IN_DIR)
for dset in os.listdir(IN_DIR):
    #print(IN_DIR+dset+'/')
    ## Loop over processings / eras
    for proc in os.listdir(IN_DIR+dset+'/'):
        #print(IN_DIR+dset+'/'+proc+'/')
        ## Loop over crab date tags
        for crab in os.listdir(IN_DIR+dset+'/'+proc+'/'):
            if crab.endswith('.root'): continue
            #print(IN_DIR+dset+'/'+proc+'/'+crab+'/')
            ## Loop over sub-directories
            for subd in os.listdir(IN_DIR+dset+'/'+proc+'/'+crab+'/'):
                print(IN_DIR+dset+'/'+proc+'/'+crab+'/'+subd+'/')
                ## Loop over files
                for fname in os.listdir(IN_DIR+dset+'/'+proc+'/'+crab+'/'+subd+'/'):
                    if not fname.endswith('Skim.root'): continue
                    #print(fname)
                    strs = fname.split('_')
                    #print(strs)
                    new_name = 'PNet_v1_Skim_'+strs[2]+'.root'
                    #print(new_name)
                    os.rename(IN_DIR+dset+'/'+proc+'/'+crab+'/'+subd+'/'+fname,
                              IN_DIR+dset+'/'+proc+'/'+crab+'/'+subd+'/'+new_name)

print('\n*** All done!!! ***\n\n')
