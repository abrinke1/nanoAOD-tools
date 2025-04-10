#! /usr/bin/env python
## Counts events passing different categories according to Haa4b postprocessor branches

import os
import sys
import subprocess
import numpy as np
import ROOT as R

R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn

MAX_EVT = -1     ## Maximum number of events to process per MC sample
PRT_EVT = 10000  ## Print every Nth event while processing
VERBOSE = False
DEBUG   = [] ## [luminosityBlock, event] to debug

# ## Location of postprocessed input files
# IN_DIR = '/eos/cms/store/user/abrinke1/NanoPost/hadd/'
# ## Haa4b selection categories, as of end-of-summer 2024
# CATS = ['gg0lHi','gg0lLo','VBFjjHi','VBFjjLo','ZvvHi','ZvvLo','Vjj','ttHad',
#        'WlvHi_m','WlvHi_e','WlvLo_m','WlvLo_e','Zll_m','Zll_e',
#        'ttbblv_m','ttbblv_e','ttblv_m','ttblv_e','ttll_m','ttll_e',
#        '2lSS','3l','other']
## MC samples to test
#SAMPS = ['GluGluH_M-15','GluGluH_M-25','GluGluH_M-40','GluGluH_M-55','VBFH_M-30',
#         'WH_M-20','ZH_M-20','ZH_M-40','TTH_M-15','TTH_M-55']

CATS = ['mu_sel','mu_sel_trig','mu_sel_trigFat','mu_sel_trigBtag',
        'mu_sel_0b','mu_sel_0b_trig','mu_sel_0b_trigFat','mu_sel_0b_trigBtag']

# SAMPS = ['SingleMuon']
# IN_DIR = '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/data/PNet_v2_2024_11_22/SingleMuon/r1_Run2018C/'
SAMPS = ['TTToSemiLeptonic']
IN_DIR = '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/PNet_v2_2024_11_22/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/r1/'


## Count passing events for each sample in each category
count = {}
for samp in SAMPS:

    print('\n\n*** Beginning to look at sample %s ***\n' % samp)

    count[samp] = {}
    for cat in CATS:
        count[samp][cat] = 0

    in_file_names = []
    for fn in os.listdir(IN_DIR):
        if samp == 'SingleMuon':
            if fn.endswith('7.root'):
                in_file_names.append(IN_DIR+fn)
        elif samp == 'TTToSemiLeptonic':
            if fn.endswith('9_8.root'):
                in_file_names.append(IN_DIR+fn)
        elif fn.startswith(samp+'_Skim'):
            in_file_names.append(IN_DIR+fn)

    chains = {}
    chains['Events'] = 0

    ## Combine input files into a single "chain"
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

        if len(DEBUG) > 0:
            if ch.luminosityBlock != DEBUG[0] or ch.event != DEBUG[1]:
                continue
            print('Starting LS = %d, event = %d' % (ch.luminosityBlock, ch.event))

        if VERBOSE: print('tagHaa4b_v1 = %.4f, cat_idx = %d, passFilters = %d' % (ch.Haa4b_FatH_tagHaa4b_v1,
                                                                                  ch.Haa4b_cat_idx,
                                                                                  ch.Haa4b_passFilters))
        #if ch.Haa4b_FatH_tagHaa4b_v1 < 0.975:    continue  ## Minimum Haa4b vs. QCD cut, WP60
        if ch.Haa4b_isHad + ch.Haa4b_isLep != 1: continue  ## AK8 Higgs candidate with no overlapping trigger lepton
        if ch.Haa4b_passFilters != 1:            continue  ## Must pass event filters
        if VERBOSE: print('Passed initial selection!')

        ## gg0l
        if 'gg0lHi' in CATS or 'gg0lLo' in CATS:
            if ch.Haa4b_cat_gg0l == 1 and ch.Haa4b_FatH_tagHaa4b_v1 > 0.992 and \
               (ch.Haa4b_trigFat + ch.Haa4b_trigBtag >= 1) and ch.Haa4b_nJetBtag == 0:
                if ch.Haa4b_FatH_pt > 400:
                    count[samp]['gg0lHi'] += 1
                elif ch.Haa4b_FatH_pt > 250:
                    count[samp]['gg0lLo'] += 1
        ## VBFjj
        if 'VBFjjHi' in CATS or 'VBFjjLo' in CATS:
            if ch.Haa4b_cat_VBFjj == 1 and ch.Haa4b_FatH_tagHaa4b_v1 > 0.992 and \
               (ch.Haa4b_trigFat + ch.Haa4b_trigBtag + ch.Haa4b_trigVBF) >= 1 and \
               ch.Haa4b_FatH_pt > 250 and ch.Haa4b_nJetBtag == 0:
                if ch.Haa4b_dijet_mass > 900 and ch.Haa4b_dijet_dEta > 3.0:
                    count[samp]['VBFjjHi'] += 1
                else:
                    count[samp]['VBFjjLo'] += 1
        ## Zvv
        if 'ZvvHi' in CATS or 'ZvvLo' in CATS:
            if ch.Haa4b_cat_Zvv == 1 and ch.Haa4b_trigMET == 1 and ch.Haa4b_FatH_pt > 200:
                if ch.MET_pt > 300:
                    count[samp]['ZvvHi'] += 1
                else:
                    count[samp]['ZvvLo'] += 1
        ## Vjj
        if 'Vjj' in CATS:
            if ch.Haa4b_cat_Vjj == 1 and (ch.Haa4b_trigFat + ch.Haa4b_trigBtag) >= 1 and \
               ch.Haa4b_FatH_pt > 250 and ch.Haa4b_FatX_tagWZ_max > 0.98:
                count[samp]['Vjj'] += 1
        ## ttHad
        if 'ttHad' in CATS:
            if ch.Haa4b_cat_ttHad == 1 and (ch.Haa4b_trigFat + ch.Haa4b_trigBtag) >= 1 and \
               ch.Haa4b_FatH_pt > 250 and ch.Haa4b_FatX_tagTop_max > 0.8 and ch.Haa4b_nJetBtag >= 1:
                count[samp]['ttHad'] += 1

        ## All lepton categories must have trigger lepton and fire trigger
        if (ch.Haa4b_isMu + ch.Haa4b_trigMu) < 2 and \
           (ch.Haa4b_isEle + ch.Haa4b_trigEle) < 2: continue
        flav = ('_m' if ch.Haa4b_isMu else '_e')

        ## Wlv
        if 'WlvHi'+flav in CATS or 'WlvLo'+flav in CATS:
            if ch.Haa4b_cat_Wlv == 1 and ch.Haa4b_nJetBtag == 0:
                pure = ('Hi' if ch.Haa4b_lepMET_pt > 300 else 'Lo')
                count[samp]['Wlv'+pure+flav] += 1
        ## Zll
        if 'Zll'+flav in CATS:
            if ch.Haa4b_cat_Zll == 1:
                count[samp]['Zll'+flav] += 1
        ## ttlv
        if 'ttblv'+flav in CATS or 'ttbblv'+flav in CATS:
            if ch.Haa4b_iLep1 >= 0 and ch.Haa4b_iLep2 < 0 and ch.Haa4b_nJetBtag >= 1:
                nbs = ('bb' if ch.Haa4b_nJetBtag >= 2 else 'b')
                count[samp]['tt'+nbs+'lv'+flav] += 1
        ## ttll
        if 'ttll'+flav in CATS:
            if ch.Haa4b_dilep_charge == 0 and ch.Haa4b_dilep_mass > 12 and ch.Haa4b_nJetBtag >= 1:
                if ch.Haa4b_cat_Zll == 0 and ch.Haa4b_cat_3l == 0:
                    count[samp]['ttll'+flav] += 1
        ## 2lSS, 3l, other
        if '2lSS' in CATS:
            if ch.Haa4b_cat_2lSS == 1:
                count[samp]['2lSS'] += 1
        if '2l' in CATS:
            if ch.Haa4b_cat_3l == 1:
                count[samp]['3l'] += 1
        if 'other' in CATS:
            if ch.Haa4b_cat_other == 1:
                count[samp]['other'] += 1

        ## SingleMuon control region for trigger studies
        if 'mu_sel' in CATS:
            if ch.Haa4b_isMu == 1 and ch.Haa4b_trigMu == 1:
                count[samp]['mu_sel'] += 1
                count[samp]['mu_sel_0b'] += (ch.Haa4b_nJetBtag == 0)
                if (ch.Haa4b_trigFat == 1 or ch.Haa4b_trigBtag == 1):
                    count[samp]['mu_sel_trig'] += 1
                    count[samp]['mu_sel_0b_trig'] += (ch.Haa4b_nJetBtag == 0)
                if ch.Haa4b_trigFat == 1:
                    count[samp]['mu_sel_trigFat'] += 1
                    count[samp]['mu_sel_0b_trigFat'] += (ch.Haa4b_nJetBtag == 0)
                if ch.Haa4b_trigBtag == 1:
                    count[samp]['mu_sel_trigBtag'] += 1
                    count[samp]['mu_sel_0b_trigBtag'] += (ch.Haa4b_nJetBtag == 0)

        if len(DEBUG) > 0:
            if ch.luminosityBlock != DEBUG[0] or ch.event != DEBUG[1]:
                print('Finished LS = %d, event = %d' % (ch.luminosityBlock, ch.event))
                break
    ## End loop: for iEvt in range(nEntries)

    print('\n*** Finished looking at sample %s ***\n\n' % samp)

## End loop: for samp in SAMPS

print('\n*** Finished looking at all samples! ***\n\n')

for samp in SAMPS:
    for cat in CATS:
        print('%s %s = %d' % (samp, cat, count[samp][cat]))

print('\n*** All done!!! ***\n\n')
