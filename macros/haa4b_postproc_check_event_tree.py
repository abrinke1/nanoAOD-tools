#! /usr/bin/env python
## Check that Event trees are not empty, and count total events

import os
import sys
import subprocess
import numpy as np
import ROOT as R

R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn

## Location of postprocessed input files
IS_DATA = False
HADD_ONLY = False
EVENT_TREE = False  ## Use 'Event' tree even for MC (slower, counts only "passing" / saved events)
IN_DIR = '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/%s/PNet_v2_2024_11_22/' % ('data' if IS_DATA else 'MC')

## Loop over datasets
#print(IN_DIR)
for dset in os.listdir(IN_DIR):
    #if dset == 'JetHT' or dset == 'EGamma': continue
    if not 'DYJetsToLL_M-50_HT-200to400' in dset: continue
    #if (not 'HToAATo4B_Pt150_M-' in dset or not 'VBFH' in dset): continue
    #print(IN_DIR+dset+'/')
    ## Loop over processings / eras
    for proc in os.listdir(IN_DIR+dset+'/'):
        #if not (dset+'/'+proc == 'EGamma/r1_Run2018D'): continue
        print('\nStarting to look at '+dset+'/'+proc)
        hadd_sum = 0
        hadd_neg = 0
        hadd_files = 0
        hadd_empty = 0
        nano_sum = 0
        skim_sum = 0
        nano_files = 0
        skim_files = 0
        nano_empty = 0
        skim_empty = 0
        #print(IN_DIR+dset+'/'+proc+'/')
        ## Loop over crab date tags (and hadded files)
        for crab in os.listdir(IN_DIR+dset+'/'+proc+'/'):
            ## Count events in hadded files
            if crab.endswith('.root'):
                hadd_files += 1
                in_file_hadd = IN_DIR+dset+'/'+proc+'/'+crab
                ## Get number of entries from 'Events' tree for data
                if IS_DATA or EVENT_TREE:
                    chain_hadd = R.TChain('Events')
                    chain_hadd.Add( in_file_hadd )
                    try:    nEntries_hadd = chain_hadd.GetEntries()
                    except: nEntries_hadd = 0
                    hadd_sum += nEntries_hadd
                ## Get number of processed events from 'Runs' tree for MC (faster)
                else:
                    chain_hadd = R.TChain('Runs')
                    chain_hadd.Add( in_file_hadd )
                    nEntries_hadd = 0
                    nEntries_neg  = 0
                    try:
                        for iEntry in range(chain_hadd.GetEntries()):
                            chain_hadd.GetEntry(iEntry)
                            nEntries_hadd += chain_hadd.genEventCount
                            nEntries_neg  += chain_hadd.genEventCountNeg
                    except: nEntries_hadd += 0
                    hadd_sum += nEntries_hadd
                    hadd_neg += nEntries_neg

                if nEntries_hadd <= 0:
                    print('HADD: %d entries in %s' % (nEntries_hadd, in_file_hadd))
                    hadd_empty += 1
                del chain_hadd
                continue
            if HADD_ONLY: continue
            #print(IN_DIR+dset+'/'+proc+'/'+crab+'/')
            if crab.endswith('.sh'): continue
            ## Loop over sub-directories
            for subd in os.listdir(IN_DIR+dset+'/'+proc+'/'+crab+'/'):
                print(IN_DIR+dset+'/'+proc+'/'+crab+'/'+subd+'/')
                if subd.endswith('.sh'): continue
                ## Loop over files
                for fname in os.listdir(IN_DIR+dset+'/'+proc+'/'+crab+'/'+subd+'/'):
                    if not fname.endswith('.root'): continue
                    kind = 'SKIM' if 'Skim' in fname else 'NANO'
                    nano_files += (kind == 'NANO')
                    skim_files += (kind == 'SKIM')
                    in_file = IN_DIR+dset+'/'+proc+'/'+crab+'/'+subd+'/'+fname
                    #print(in_file)
                    ## Get number of entries from 'Events' tree for data
                    if IS_DATA or EVENT_TREE:
                        chain = R.TChain('Events')
                        chain.Add( in_file )
                        try:    nEntries = chain.GetEntries()
                        except: nEntries = 0
                    ## Get number of processed events from 'Runs' tree for MC (faster)
                    else:
                        chain = R.TChain('Runs')
                        chain.Add( in_file )
                        nEntries = 0
                        try:
                            for iEntry in range(chain.GetEntries()):
                                chain.GetEntry(iEntry)
                                nEntries += chain.genEventCount
                        except: nEntries += 0

                    nano_sum += (nEntries if kind == 'NANO' else 0)
                    skim_sum += (nEntries if kind == 'SKIM' else 0)
                    if nEntries <= 0:
                        print('%s: %d entries in %s' % (kind, nEntries, in_file))
                        nano_empty += (kind == 'NANO')
                        skim_empty += (kind == 'SKIM')
                    del chain
                ## End loop over files
            ## End loop over sub-directories
        ## End loop over crab date tags (and hadded files)
        print('\n*** Summary for '+dset+'/'+proc+' ***')
        print('HADD has %d events, %d / %d empty files' % (hadd_sum, hadd_empty, hadd_files))
        if not IS_DATA and not EVENT_TREE and hadd_neg > 0:
            print('         %d negative-weight events!' % hadd_neg)
        print('SKIM has %d events, %d / %d empty files' % (skim_sum, skim_empty, skim_files))
        print('NANO has %d events, %d / %d empty files' % (nano_sum, nano_empty, nano_files))

print('\n*** All done!!! ***\n\n')
