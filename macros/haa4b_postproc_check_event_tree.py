#! /usr/bin/env python
## Check that Event trees are not empty

import os
import sys
import subprocess
import numpy as np
import ROOT as R

R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn

## Location of postprocessed input files
IN_DIR = '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/data/PNet_v2_2024_11_22/re_post/'

## Loop over datasets
#print(IN_DIR)
for dset in os.listdir(IN_DIR):
    #print(IN_DIR+dset+'/')
    ## Loop over processings / eras
    for proc in os.listdir(IN_DIR+dset+'/'):
        #print(IN_DIR+dset+'/'+proc+'/')
        ## Loop over crab date tags
        for crab in os.listdir(IN_DIR+dset+'/'+proc+'/'):
            #print(IN_DIR+dset+'/'+proc+'/'+crab+'/')
            ## Loop over sub-directories
            for subd in os.listdir(IN_DIR+dset+'/'+proc+'/'+crab+'/'):
                print(IN_DIR+dset+'/'+proc+'/'+crab+'/'+subd+'/')
                ## Loop over files
                for fname in os.listdir(IN_DIR+dset+'/'+proc+'/'+crab+'/'+subd+'/'):
                    in_file = IN_DIR+dset+'/'+proc+'/'+crab+'/'+subd+'/'+fname
                    #print(in_file)
                    chain = R.TChain('Events')
                    chain.Add( in_file )
                    nEntries = chain.GetEntries()
                    if nEntries <= 0:
                        print('%d entries in %s' % (nEntries, in_file))
                    del chain

print('\n*** All done!!! ***\n\n')
