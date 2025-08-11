#! /usr/bin/env python
## Compute ParticleNet WZvsQCD thresholds with the same efficiency as official WvsQCD thresholds
## https://twiki.cern.ch/twiki/bin/viewauth/CMS/ParticleNetSFs#2016_RunBCDEF

import os
import sys
import subprocess
import numpy as np
import ROOT as R

R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn

MAX_EVT = -1  ## Maximum number of events to process per MC sample
PRT_EVT = 10000  ## Print every Nth event while processing
DEBUG   = False
NBINS   = 10000  ## Number of bins in histogram

## Location of postprocessed input files
IN_DIRS  = '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/YEAR/MC/PNet_v2_2024_11_22/'
IN_FILES = 'SUSY_WH_WToAll_HToAATo4B_Pt150_M-mA_TuneCP5_13TeV_madgraph_pythia8/r1/PNet_v1_Skim_Lp_weight.root'

YEARS  = ['2016APV','2016','2017','2018']
MASSES = ['12']+[str(15+5*i) for i in range(10)]
#YEARS  = ['2018']
#MASSES = ['30']
CUTS   = {'2016APV':0.974, '2016':0.974, '2017':0.978, '2018':0.980}

## Histograms for WvsQCD and WZvsQCD distributions
hst = {}

for year in YEARS:
    print('\n\n*** Beginning to look at year %s ***\n' % year)
    hst[year] = {}
    hst[year]['W']  = R.TH1D('h_WvsQCD_%s' % year, 'PNet WvsQCD %s' % year,  NBINS, 0, 1)
    hst[year]['WZ'] = R.TH1D('h_WZvsQCD_%s' % year,'PNet WZvsQCD %s' % year, NBINS, 0, 1)

    ## Combine input files into a single "chain"
    in_file_names = [IN_DIRS.replace('YEAR',year)+IN_FILES.replace('mA',mA) for mA in MASSES]
    chains = {}
    chains['Events'] = 0
    for i in range(len(in_file_names)):
        print('Adding file %s' % in_file_names[i])
        for key in chains.keys():
            if i == 0: chains[key] = R.TChain(key)
            chains[key].Add( in_file_names[i] )

    ## Loop through events, select, and count
    ch = chains['Events']  ## Shortcut expression
    nEntries = ch.GetEntries()
    print('\nEntering loop over %d events\n' % (nEntries))

    for iEvt in range(nEntries):

        if iEvt > MAX_EVT and MAX_EVT > 0: break
        if (iEvt % PRT_EVT) == 0: print('Looking at event #%d / %d' % (iEvt, nEntries))

        ch.GetEntry(iEvt)

        if ch.Haa4b_nFatX < 1:  continue  ## Good W/Z candidate AK8 must exist
        if ch.Haa4b_iFatH < 0:  continue  ## Good Haa4b candidate AK8 must exist
        if ch.Haa4b_isLep == 1: continue  ## Cannot have triggerable lepton
        if ch.Haa4b_FatH_tagHaa4b_v2a + \
           ch.Haa4b_FatH_tagHaa4b_v2b < 2.0*0.84: continue  ## Higgs signal must pass WP80

        for iJ in range(ch.nFatJet):
            if iJ == ch.Haa4b_iFatH: continue  ## Cannot be Higgs candidate
            if ch.FatJet_Haa4b_candX[iJ] != 1: continue  ## Must be W/Z candidate
            if abs(ch.FatJet_genPart_pdgId[iJ]) != 24: continue  ## Must be GEN-matched to true W boson
            hst[year]['W'] .Fill(min(max(ch.FatJet_particleNet_WvsQCD[iJ],  0.1/NBINS), 1.0-0.1/NBINS))
            hst[year]['WZ'].Fill(min(max(ch.FatJet_particleNet_WZvsQCD[iJ], 0.1/NBINS), 1.0-0.1/NBINS))
            break

    ## End loop: for iEvt in range(nEntries)
## End loop: for year in YEARS

print('\n*** Finished looking at all years! ***\n\n')

for year in YEARS:
    cutW = CUTS[year]
    iW = hst[year]['W'].Integral()
    nW = hst[year]['W'].Integral(int(cutW*NBINS)+1, NBINS)
    print('\nIn year %s, WvsQCD cut = %.3f, efficiency = %.2f%% (%d / %d)' % (year, cutW, 100.0*nW/iW, nW, iW))
    cutWZ = -1.0
    for iX in reversed(range(1, NBINS+1)):
        if hst[year]['WZ'].Integral(iX, NBINS) > nW:
            iWZ = hst[year]['WZ'].Integral()
            nWZ = hst[year]['WZ'].Integral(iX, NBINS)
            mWZ = hst[year]['WZ'].Integral(iX+1, NBINS)
            cutWZn = hst[year]['WZ'].GetBinLowEdge(iX)
            cutWZm = hst[year]['WZ'].GetBinLowEdge(iX+1)
            print('             WZvsQCD cut = %.5f, efficiency = %.3f%% (%d / %d)' % (cutWZn, 100.0*nWZ/iWZ, nWZ, iWZ))
            print('             WZvsQCD cut = %.5f, efficiency = %.3f%% (%d / %d)' % (cutWZm, 100.0*mWZ/iWZ, mWZ, iWZ))
            break
    ## End loop: for iX in reversed(range(1, NBINS+1))
## End loop: for year in YEARS

print('\n*** All done!!! ***\n\n')
