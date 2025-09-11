#! /usr/bin/env python
## Counts events passing leptonic categories

import os
import sys
import subprocess
import numpy as np
import ROOT as R

R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn

MAX_EVT = -1  ## Maximum number of events to process per MC sample
PRT_EVT = 10000    ## Print every Nth event while processing
VERBOSE = False
DEBUG   = [] ## [luminosityBlock, event] to debug
YEAR = '2018'

# SAMPS = ['SingleMuon']
# IN_DIR = '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/%s/data/PNet_v2_2024_11_22/SingleMuon/r1_Run%sC/' % (YEAR,YEAR)
# SAMPS = ['TTToSemiLeptonic']
# IN_DIR = '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/%s/MC/PNet_v2_2024_11_22/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/r1/' % YEAR
SAMPS = ['TTH_HToAATo4B_M-45']
IN_DIR = '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/%s/MC/PNet_v2_2024_11_22/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8/r1/' % YEAR

CUTS = ['all','nPV','noise',
        'fatH_presel','mu_veto','ele_veto',
        'bjet_veto','MET_veto','fatWZ_veto',
        'VBFjj','VBFjjLo','VBFjjHi',
        'VBFjjLoPtLo','VBFjjLoPtHi','VBFjjHiPtLo','VBFjjHiPtHi',
        'X4bSB','X4bSR','4GenB']
CATS = {}
CATS['VBFjj'] = [['all'],['nPV'],['noise'],['fatH_presel'],
                 ['trigJetHT','trigBTag'],['mu_veto'],['ele_veto'],
                 ['bjet_veto'],['fatWZ_veto'],['MET_veto'],
                 ['VBFjj'],['VBFjjLo','VBFjjHi'],
                 ['VBFjjLoPtLo','VBFjjLoPtHi','VBFjjHiPtLo','VBFjjHiPtHi'],
                 ['X4bSB'],['X4bSR'],['4GenB']]


## Function to iterate counters based on per-event flags
def fill_counts(counter, selected):
    for cat in CATS.keys():
        for iCuts in range(len(CATS[cat])):
            if len(CATS[cat][iCuts]) == 1:
                cutstr = CATS[cat][iCuts][0]
                key = '%d_%s' % (iCuts, cutstr)
                if selected[cutstr]:
                    counter[cat][key] = counter[cat][key] + 1
                else: ## Subsequent cuts depend on earlier cuts
                    break
            else:
                atleastone = False
                for jCut in range(len(CATS[cat][iCuts])):
                    cutstr = CATS[cat][iCuts][jCut]
                    key = '%d_%s' % (iCuts, cutstr)
                    if selected[cutstr]:
                        counter[cat][key] = counter[cat][key] + 1
                        atleastone = True
                if atleastone:
                    counter[cat]['%d_any' % iCuts] = counter[cat]['%d_any' % iCuts] + 1
                else: ## Subsequent cuts depend on earlier cuts
                    break
        ## End loop: for iCuts in range(len(CATS[cat]))
    ## End loop: for cat in CATS.keys()
    return counter
## End function: fill_counts(counter, selected)    


########################################################
##                  MAIN FUNCTION
## Count passing events for each sample in each category
########################################################
count = {}
for samp in SAMPS:

    print('\n\n*** Beginning to look at sample %s ***\n' % samp)

    count[samp] = {}
    for cat in CATS.keys():
        count[samp][cat] = {}
        for iCuts in range(len(CATS[cat])):
            if len(CATS[cat][iCuts]) == 1:
                count[samp][cat]['%d_%s' % (iCuts, CATS[cat][iCuts][0])] = 0
            else:
                count[samp][cat]['%d_any' % iCuts] = 0
                for jCut in range(len(CATS[cat][iCuts])):
                    count[samp][cat]['%d_%s' % (iCuts, CATS[cat][iCuts][jCut])] = 0
        ## End loop: for iCuts in range(len(CATS[cat]))
    ## End loop: for cat in CATS.keys()

    in_file_names = []
    for fn in os.listdir(IN_DIR):
        if samp == 'SingleMuon':
            if fn.endswith('7.root'):
                in_file_names.append(IN_DIR+fn)
        elif samp == 'TTToSemiLeptonic':
            if fn.endswith('9_8.root'):
                in_file_names.append(IN_DIR+fn)
        elif samp == 'TTH_HToAATo4B_M-45':
            if fn.endswith('.root'):
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
    
    ##################
    ## MAIN EVENT LOOP
    for iEvt in range(nEntries):

        if iEvt >= MAX_EVT and MAX_EVT > 0: break
        if (iEvt % PRT_EVT) == 0: print('Looking at event #%d / %d' % (iEvt, nEntries))

        ch.GetEntry(iEvt)

        if len(DEBUG) > 0:
            if ch.luminosityBlock != DEBUG[0] or ch.event != DEBUG[1]:
                continue
            print('Starting LS = %d, event = %d' % (ch.luminosityBlock, ch.event))

        ## Initialize per-event bools for each cut
        sel = {}
        for cut in CUTS:
           sel[cut] = False
        sel['all'] = True

        ## Table 5 in Ch. 4 of AN2023_047_v4
        sel['noise'] = (ch.Flag_goodVertices and ch.Flag_globalSuperTightHalo2016Filter and
                        ch.Flag_HBHENoiseFilter and ch.Flag_HBHENoiseIsoFilter and
                        ch.Flag_eeBadScFilter and ch.Flag_BadPFMuonFilter and
                        ch.Flag_BadPFMuonDzFilter and ch.Flag_EcalDeadCellTriggerPrimitiveFilter and
                        ('2016' in YEAR or ch.Flag_ecalBadCalibFilter))
        if not sel['noise']:
            count[samp] = fill_counts(count[samp], sel)
            continue

        sel['nPV'] = (ch.PV_npvsGood > 0)
        if not sel['nPV']:
            count[samp] = fill_counts(count[samp], sel)
            continue

        ## Higgs candidate AK8 jet selection from Table 1 in Ch. 4 of AN2023_047_v4
        xFatH = -99
        xFatH_X4b = -99
        iFatWZs = []
        for iFat in range(ch.nFatJet):
            if      ch.FatJet_pt[iFat]   <= 250: continue
            if  abs(ch.FatJet_eta[iFat]) >= 2.4: continue
            if     ch.FatJet_jetId[iFat] != 6:   continue
            if ch.FatJet_msoftdrop[iFat] <= 20:  continue
            if (YEAR == '2018' and ch.FatJet_particleNet_WZvsQCD[iFat] > 0.9873) or \
               (YEAR == '2017' and ch.FatJet_particleNet_WZvsQCD[iFat] > 0.9858) or \
               ('2016' in YEAR and ch.FatJet_particleNet_WZvsQCD[iFat] > 0.9843):
                iFatWZs.append(iFat)
            if ch.FatJet_particleNetMD_XbbvsQCD[iFat] <= 0.75: continue
            iFatH_X4b = 0.5*(ch.FatJet_PNet_X4b_v2a_Haa4b_score[iFat] + \
                             ch.FatJet_PNet_X4b_v2b_Haa4b_score[iFat])
            if iFatH_X4b > xFatH_X4b:
                xFatH     = iFat
                xFatH_X4b = iFatH_X4b
        ## End loop: for iFat in range(ch.nFatJet)
        if xFatH < 0:
            count[samp] = fill_counts(count[samp], sel)
            continue
        sel['fatH_presel'] = True
        vFatH = R.TLorentzVector()
        vFatH.SetPtEtaPhiM(ch.FatJet_pt[xFatH], ch.FatJet_eta[xFatH],
                           ch.FatJet_phi[xFatH], ch.FatJet_mass[xFatH])

        ## Trigger: Table 18 in Appendix C of AN2023_047_v4
        if YEAR == '2018':
            sel['trigJetHT'] = (ch.HLT_AK8PFJet500 or ch.HLT_PFHT1050 or ch.HLT_PFJet500 or \
                                ch.HLT_AK8PFHT800_TrimMass50 or ch.HLT_AK8PFJet400_TrimMass30)
            sel['trigBTag']  = (ch.HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4 or \
                                ch.HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71 or \
                                ch.HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2 or \
                                ch.HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1)
        elif YEAR == '2017':
            sel['trigJetHT'] = (ch.HLT_AK8PFJet500 or ch.HLT_PFHT1050 or ch.HLT_PFJet500 or \
                                ch.HLT_AK8PFHT800_TrimMass50 or ch.HLT_AK8PFJet400_TrimMass30 or \
                                ch.HLT_AK8PFJet360_TrimMass30)
            sel['trigBTag'] = (ch.HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2 or \
                               ch.HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2 or \
                               ch.HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5 or \
                               ch.HLT_AK8PFHT750_TrimMass50 or ch.HLT_AK8PFJet380_TrimMass30 or \
                               ch.HLT_DoublePFJets100MaxDeta1p6_DoubleCaloBTagCSV_p33 or \
                               ch.HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0)
        elif '2016' in YEAR:
            sel['trigJetHT'] = (ch.HLT_PFJet450 or ch.HLT_DiCentralPFJet430 or \
                                ch.HLT_PFHT650_WideJetMJJ900DEtaJJ1p5 or ch.HLT_PFHT750_4JetPt50 or \
                                ch.HLT_PFHT800 or ch.HLT_PFHT900 or ch.HLT_AK8PFJet360_TrimMass30 or \
                                ch.HLT_AK8PFJet450 or ch.HLT_AK8PFHT650_TrimR0p1PT0p03Mass50 or \
                                ch.HLT_AK8PFHT700_TrimR0p1PT0p03Mass50)
            sel['trigBTag']  = (ch.HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20 or \
                                ch.HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20 or \
                                ch.HLT_PFHT400_SixJet30_DoubleBTagCSV_p056 or \
                                ch.HLT_PFHT450_SixJet40_BTagCSV_p056 or \
                                ch.HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20 or \
                                ch.HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160 or \
                                ch.HLT_DoubleJetsC112_DoubleBTagCSV_p026_DoublePFJetsC172 or \
                                ch.HLT_DoubleJet90_Double30_TripleBTagCSV_p08 or \
                                ch.HLT_QuadJet45_TripleBTagCSV_p087)
        if not (sel['trigJetHT'] or sel['trigBTag']):
            count[samp] = fill_counts(count[samp], sel)
            continue

        ## Muon selection from Table 2 in Ch. 4 of AN2023_047_v4
        hasTrgMu = False
        for iMu in range(ch.nMuon):
            if ch.Muon_pt[iMu] <= 26.0 + 3.0*(YEAR == '2017'): continue
            if abs(ch.Muon_eta[iMu]) >= 2.4:  continue
            if ch.Muon_mediumId[iMu] == 0 and \
               (ch.Muon_pt[iMu] <= 200 or
                ch.Muon_highPtId[iMu] == 0):  continue
            if ch.Muon_miniPFRelIso_all[iMu] >= 0.10: continue
            if abs(ch.Muon_dxy[iMu])         >= 0.02: continue
            if abs(ch.Muon_dz[iMu])          >= 0.10: continue
            hasTrgMu = True
            break
        ## End loop: for iMu in range(ch.nMuon)
        sel['mu_veto'] = (not hasTrgMu)
        if not sel['mu_veto']:
            count[samp] = fill_counts(count[samp], sel)
            continue

        ## Electron selection from Table 3 in Ch. 4 of AN2023_047_v4
        hasTrgEle = False
        for iEle in range(ch.nElectron):
            if ch.Electron_pt[iEle] <= 30.0 + 5.0*(not '2016' in YEAR): continue
            if         abs(ch.Electron_eta[iEle]) >= 2.5: continue
            if ch.Electron_mvaFall17V2Iso_WP90[iEle] == 0: continue
            if ch.Electron_mvaFall17V2Iso_WP80[iEle] == 0 and \
               (ch.Electron_pt[iEle] <= 35 or
                ch.Electron_cutBased_HEEP[iEle] == 0): continue
            if abs(ch.Electron_dxy[iEle])     >= 0.02: continue
            if abs(ch.Electron_dz[iEle])      >= 0.10: continue
            hasTrgEle = True
            break
        ## End loop: for iEle in range(ch.nElectron)
        sel['ele_veto'] = (not hasTrgEle)
        if not sel['ele_veto']:
            count[samp] = fill_counts(count[samp], sel)
            continue

        ## AK4 jet selection from Table 4 in Ch. 4 of AN2023_047_v4
        nBJet = 0
        btagWPM = 0.2783 if YEAR == '2018' else (0.3040 if YEAR == '2017' else (0.2489 if 'post' in YEAR else 0.2598))
        sel['bjet_veto'] = True
        vJets = []
        for iJet in range(ch.nJet):
            if     ch.Jet_pt[iJet]   <= 30: continue
            if    ch.Jet_jetId[iJet] != 6:  continue
            if   ch.Jet_puId[iJet] < 4 and \
                    ch.Jet_pt[iJet]  <= 50: continue
            vJet = R.TLorentzVector()
            vJet.SetPtEtaPhiM(ch.Jet_pt[iJet], ch.Jet_eta[iJet],
                              ch.Jet_phi[iJet], ch.Jet_mass[iJet])
            if vJet.DeltaR(vFatH) <= 0.8: continue
            ## Check for b-tagged AK4 jets near Higgs AK8
            if       abs(ch.Jet_eta[iJet]) < 2.4 and \
                ch.Jet_btagDeepFlavB[iJet] > btagWPM:
                vBJet = R.TLorentzVector()
                if vBJet.DeltaR(vFatH) <= 1.2:
                    sel['bjet_veto'] = False
                continue
            ## Store non-b-tagged jets for VBFjj selection
            vJets.append(vJet)
        ## End loop: for iJet in range(ch.nJet)

        # vMET = R.TLorentzVector()
        # vMET.SetPtEtaPhiM(ch.MET_pt, 0, ch.MET_phi, 0)

        if not sel['bjet_veto']:
            count[samp] = fill_counts(count[samp], sel)
            continue
        sel['fatWZ_veto'] = (len(iFatWZs) == 0 or (len(iFatWZs) == 1 and iFatWZs[0] == xFatH))
        if not sel['fatWZ_veto']:
            count[samp] = fill_counts(count[samp], sel)
            continue
        sel['MET_veto'] = (ch.MET_pt <= 200)
        if not sel['MET_veto']:
            count[samp] = fill_counts(count[samp], sel)
            continue
        if len(vJets) >= 2:
            j1 = vJets[0]
            j2 = vJets[1]
        sel['VBFjj'] = (len(vJets) >= 2 and abs(j1.Eta() - j2.Eta()) > 2.2 and (j1+j2).M() > 450)
        if not sel['VBFjj']:
            count[samp] = fill_counts(count[samp], sel)
            continue

        sel['VBFjjHi'] = (abs(j1.Eta() - j2.Eta()) > 3.0 and (j1+j2).M() > 900)
        sel['VBFjjLo'] = (not sel['VBFjjHi'])
        for xSel in ['VBFjjLo','VBFjjHi']:
            sel[xSel+'PtLo'] = (sel[xSel] and vFatH.Pt()  < 400)
            sel[xSel+'PtHi'] = (sel[xSel] and vFatH.Pt() >= 400)

        sel['X4bSB'] = (xFatH_X4b > 0.84)
        sel['X4bSR'] = (xFatH_X4b > 0.96)
        sel['4GenB'] = ch.FatJet_nBHadrons[xFatH] >= 4
        count[samp] = fill_counts(count[samp], sel)

    ## End loop: for iEvt in range(nEntries)

    print('\n*** Finished looking at sample %s ***\n' % samp)

## End loop: for samp in SAMPS

print('\n*** Finished looking at all samples! ***\n\n')

for samp in SAMPS:
    for cat in CATS.keys():
        print('\nCut-flow for %s in %s:' % (cat, samp))
        for iCuts in range(len(CATS[cat])):
            if len(CATS[cat][iCuts]) == 1:
                cutstr = CATS[cat][iCuts][0]
                key = '%d_%s' % (iCuts, cutstr)
                print('%s : %d' % (key, count[samp][cat][key]))
            else:
                for jCut in range(len(CATS[cat][iCuts])):
                    cutstr = CATS[cat][iCuts][jCut]
                    key = '%d_%s' % (iCuts, cutstr)
                    print('  * %s : %d' % (key, count[samp][cat][key]))
                print('%d_any : %d' % (iCuts, count[samp][cat]['%d_any' % iCuts]))
        ## End loop: for iCuts in range(len(CATS[cat]))
    ## End loop: for cat in CATS.keys()
## End loop: for samp in SAMPS

print('\n*** All done!!! ***\n\n')
