#! /usr/bin/env python
## Counts events passing leptonic categories

import os
import sys
import subprocess
import numpy as np
import ROOT as R

R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn

MAX_EVT = -1  ## Maximum number of events to process per MC sample
PRT_EVT = 1000    ## Print every Nth event while processing
VERBOSE = False
DEBUG   = [] ## [luminosityBlock, event] to debug
YEAR = '2018'

# SAMPS = ['SingleMuon']
# IN_DIR = '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/%s/data/PNet_v2_2024_11_22/SingleMuon/r1_Run%sC/' % (YEAR,YEAR)
# SAMPS = ['TTToSemiLeptonic']
# IN_DIR = '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/%s/MC/PNet_v2_2024_11_22/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/r1/' % YEAR
SAMPS = ['TTH_HToAATo4B_M-45']
IN_DIR = '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/%s/MC/PNet_v2_2024_11_22/SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-45_TuneCP5_13TeV_madgraph_pythia8/r1/' % YEAR

CUTS = ['all','noise','trigMu','trigEle','1mu','1ele',
        'mumu','muele','elemu','eleele','lep_trig',
        'fatH_presel','fatH_lep_ovlp','1mu_sel','1ele_sel',
        'mumu_sel','muele_sel','elemu_sel','eleele_sel',
        'Zveto','Zsel','eq0b','eq1b','ge1b','ge2b',
        'ttbll','ttbmm','ttbme','ttbem','ttbee',
        'Zll','Zmm','Zee',
        'WlvLo','WlvHi','ttblv','ttbblv',
        'WmvLo','WmvHi','ttbmv','ttbbmv',
        'WevLo','WevHi','ttbev','ttbbev',
        'X4bSB','X4bSR','4GenB']
CATS = {}
CATS['ttbll'] = [['all'],['noise'],['trigMu','trigEle'],
                 ['mumu','muele','elemu','eleele'],['lep_trig'],
                 ['fatH_presel'],['fatH_lep_ovlp'],
                 ['mumu_sel','muele_sel','elemu_sel','eleele_sel'],
                 ['Zveto'],['ge1b'],['ttbll','ttbmm','ttbme','ttbem','ttbee'],
                 ['X4bSB'],['X4bSR'],['4GenB']]
CATS['Zll'] = [['all'],['noise'],['trigMu','trigEle'],
               ['mumu','eleele'],['lep_trig'],
               ['fatH_presel'],['fatH_lep_ovlp'],
               ['mumu_sel','eleele_sel'],
               ['Zsel'],['Zll','Zmm','Zee'],
               ['X4bSB'],['X4bSR'],['4GenB']]
CATS['1Lep'] = [['all'],['noise'],['trigMu','trigEle'],
                ['1mu','1ele'],['fatH_presel'],['fatH_lep_ovlp'],
                ['1mu_sel','1ele_sel'],['eq0b','eq1b','ge2b'],
                ['WlvLo','WmvLo','WevLo','WlvHi','WmvHi','WevHi',
                 'ttblv','ttbmv','ttbev','ttbblv','ttbbmv','ttbbev'],
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

        ## Tables 20 and 21 in Appendix C of AN2023_047_v4
        if YEAR == '2018':
            sel['trigMu'] = ( (ch.HLT_IsoMu24 or ch.HLT_Mu50 or ch.HLT_OldMu100 or ch.HLT_TkMu100) and \
                              (ch.L1_SingleMu22 or ch.L1_SingleMu25) )
            sel['trigEle'] = ( ch.HLT_Ele32_WPTight_Gsf or ch.HLT_Ele35_WPTight_Gsf_L1EGMT or
                               ch.HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 or
                               ch.HLT_Ele115_CaloIdVT_GsfTrkIdT )

        elif YEAR == '2017':
            sel['trigMu'] = ( (ch.HLT_IsoMu27 or ch.HLT_Mu50 or ch.HLT_OldMu100 or ch.HLT_TkMu100) and \
                              (ch.L1_SingleMu22 or ch.L1_SingleMu25) )
            sel['trigEle'] = ( ch.HLT_Ele32_WPTight_Gsf or ch.HLT_Ele35_WPTight_Gsf or
                               ch.HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 or
                               ch.HLT_Ele115_CaloIdVT_GsfTrkIdT or ch.HLT_Photon200 )
        elif '2016' in YEAR:
            sel['trigMu'] = ( (ch.HLT_IsoMu24 or ch.HLT_IsoTkMu24 or ch.HLT_Mu50 or ch.HLT_TkMu50) and \
                              (ch.L1_SingleMu22) )
            sel['trigEle'] = ( ch.HLT_Ele27_WPTight_Gsf or
                               ch.HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 or
                               ch.HLT_Ele115_CaloIdVT_GsfTrkIdT or ch.HLT_Photon175 )
        if not (sel['trigMu'] or sel['trigEle']):
            count[samp] = fill_counts(count[samp], sel)
            continue

        ## Muon selection from Table 2 in Ch. 4 of AN2023_047_v4
        iSelMu = []
        iTrgMu = []
        for iMu in range(ch.nMuon):
            if     ch.Muon_pt[iMu]   <= 10:   continue
            if abs(ch.Muon_eta[iMu]) >= 2.4:  continue
            if ch.Muon_mediumId[iMu] == 0 and \
               (ch.Muon_pt[iMu] <= 200 or
                ch.Muon_highPtId[iMu] == 0):  continue
            if ch.Muon_miniPFRelIso_all[iMu] >= 0.10: continue
            if abs(ch.Muon_dxy[iMu])         >= 0.02: continue
            if abs(ch.Muon_dz[iMu])          >= 0.10: continue
            iSelMu.append(iMu)
            if ch.Muon_pt[iMu] > 26.0 + 3.0*(YEAR == '2017'):
                vMu = R.TLorentzVector()
                vMu.SetPtEtaPhiM(ch.Muon_pt[iMu], ch.Muon_eta[iMu],
                                 ch.Muon_phi[iMu], ch.Muon_mass[iMu])
                for iTO in range(ch.nTrigObj):
                    if ch.TrigObj_id[iTO] != 13: continue
                    vTO = R.TLorentzVector()
                    vTO.SetPtEtaPhiM(ch.TrigObj_pt[iTO], ch.TrigObj_eta[iTO],
                                     ch.TrigObj_phi[iTO], 0)
                    if vTO.DeltaR(vMu) < 0.5:
                        iTrgMu.append(iMu)
                        break
        ## End loop: for iMu in range(ch.nMuon)

        ## Electron selection from Table 3 in Ch. 4 of AN2023_047_v4
        iSelEle = []
        iTrgEle = []
        for iEle in range(ch.nElectron):
            if             ch.Electron_pt[iEle]   <= 10:  continue
            if         abs(ch.Electron_eta[iEle]) >= 2.5: continue
            if ch.Electron_mvaFall17V2Iso_WPL[iEle] == 0: continue
            if ch.Electron_mvaFall17V2Iso_WP90[iEle] == 0 and \
               (ch.Electron_pt[iEle] <= 35 or
                ch.Electron_cutBased_HEEP[iEle] == 0): continue
            if abs(ch.Electron_dxy[iEle])     >= 0.02: continue
            if abs(ch.Electron_dz[iEle])      >= 0.10: continue
            iSelEle.append(iEle)
            ## TODO: add trigger matching logic
            if ch.Electron_pt[iEle] > 30.0 + 5.0*(not '2016' in YEAR) and \
               ch.Electron_mvaFall17V2Iso_WP90[iEle] == 1 and \
               (ch.Electron_mvaFall17V2Iso_WP80[iEle] == 1 or \
                (ch.Electron_pt[iEle] > 35 and \
                 ch.Electron_cutBased_HEEP[iEle] == 1)):
                vEle = R.TLorentzVector()
                vEle.SetPtEtaPhiM(ch.Electron_pt[iEle], ch.Electron_eta[iEle],
                                 ch.Electron_phi[iEle], ch.Electron_mass[iEle])
                for iTO in range(ch.nTrigObj):
                    if ch.TrigObj_id[iTO] != 11: continue
                    vTO = R.TLorentzVector()
                    vTO.SetPtEtaPhiM(ch.TrigObj_pt[iTO], ch.TrigObj_eta[iTO],
                                     ch.TrigObj_phi[iTO], 0)
                    if vTO.DeltaR(vEle) < 0.5:
                        iTrgEle.append(iEle)
                        break
        ## End loop: for iEle in range(ch.nElectron)

        sel['mumu']   = len(iSelMu) >= 2
        sel['muele']  = len(iSelMu) >= 1 and len(iSelEle) >= 1
        sel['elemu']  = sel['muele']
        sel['eleele'] = len(iSelEle) >= 2

        if len(iTrgMu) + len(iTrgEle) == 0:
            count[samp] = fill_counts(count[samp], sel)
            continue
        sel['lep_trig'] = True
        sel['1mu']  = (len(iTrgMu) == 1 and len(iTrgEle) == 0)
        sel['1ele'] = (len(iTrgMu) == 0 and len(iTrgEle) == 1)


        ## Higgs candidate AK8 jet selection from Table 1 in Ch. 4 of AN2023_047_v4
        xFatH = -99
        xFatH_X4b = -99
        for iFat in range(ch.nFatJet):
            if      ch.FatJet_pt[iFat]   <= 250: continue
            if  abs(ch.FatJet_eta[iFat]) >= 2.4: continue
            if     ch.FatJet_jetId[iFat] != 6:   continue
            if ch.FatJet_msoftdrop[iFat] <= 20:  continue
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

        ## Veto events where a trigger lepton overlaps the Higgs AK8 candidate
        nOvlp = 0
        jSelMu = []
        jSelEle = []
        vSelMu = []
        vSelEle = []
        vLep = None
        vFatH = R.TLorentzVector()
        vFatH.SetPtEtaPhiM(ch.FatJet_pt[xFatH], ch.FatJet_eta[xFatH],
                           ch.FatJet_phi[xFatH], ch.FatJet_mass[xFatH])
        for iMu in iSelMu:
            vMu = R.TLorentzVector()
            vMu.SetPtEtaPhiM(ch.Muon_pt[iMu], ch.Muon_eta[iMu],
                             ch.Muon_phi[iMu], ch.Muon_mass[iMu])
            if iMu in iTrgMu and vMu.DeltaR(vFatH) <= 0.8:
                nOvlp += 1
            elif vMu.DeltaR(vFatH) > 0.8:
                jSelMu.append(iMu)
                vSelMu.append(vMu)
                if not vLep and iMu in iTrgMu:
                    vLep = vMu

        for iEle in iTrgEle:
            vEle = R.TLorentzVector()
            vEle.SetPtEtaPhiM(ch.Electron_pt[iEle], ch.Electron_eta[iEle],
                              ch.Electron_phi[iEle], ch.Electron_mass[iEle])
            if iEle in iTrgEle and vEle.DeltaR(vFatH) <= 0.8:
                nOvlp += 1
            elif vEle.DeltaR(vFatH) > 0.8:
                jSelEle.append(iEle)
                vSelEle.append(vEle)
                if not vLep and iEle in iTrgEle:
                    vLep = vEle

        if nOvlp > 0:
            count[samp] = fill_counts(count[samp], sel)
            continue
        sel['fatH_lep_ovlp'] = True

        sel['1mu_sel']    = len(jSelMu) == 1 and len(jSelEle) == 0
        sel['1ele_sel']   = len(jSelMu) == 0 and len(jSelEle) == 1
        sel['mumu_sel']   = len(jSelMu) == 2 and len(jSelEle) == 0
        sel['muele_sel']  = len(jSelMu) == 1 and len(jSelEle) == 1 and len(iTrgMu) == 1
        sel['elemu_sel']  = len(jSelMu) == 1 and len(jSelEle) == 1 and len(iTrgMu) == 0
        sel['eleele_sel'] = len(jSelMu) == 0 and len(jSelEle) == 2
        sel_dilep_sel = (len(jSelMu) + len(jSelEle) == 2)

        sel['Zveto'] = ( (sel['mumu_sel'] and ch.Muon_charge[jSelMu[0]] + ch.Muon_charge[jSelMu[1]] == 0 and \
                          (vSelMu[0]+vSelMu[1]).M() > 12 and abs((vSelMu[0]+vSelMu[1]).M() - 91) > 10) or \
                         (sel['eleele_sel'] and ch.Electron_charge[jSelEle[0]] + ch.Electron_charge[jSelEle[1]] == 0 and \
                          (vSelEle[0]+vSelEle[1]).M() > 12 and abs((vSelEle[0]+vSelEle[1]).M() - 91) > 10) or \
                         ((sel['muele_sel'] or sel['elemu_sel']) and (vSelMu[0]+vSelEle[0]).M() > 12 and \
                          ch.Muon_charge[jSelMu[0]] + ch.Electron_charge[jSelEle[0]] == 0) )
        sel['Zsel'] = ( (sel['mumu_sel'] and ch.Muon_charge[jSelMu[0]] + ch.Muon_charge[jSelMu[1]] == 0 and \
                         abs((vSelMu[0]+vSelMu[1]).M() - 91) < 10) or \
                         (sel['eleele_sel'] and ch.Electron_charge[jSelEle[0]] + ch.Electron_charge[jSelEle[1]] == 0 and \
                          abs((vSelEle[0]+vSelEle[1]).M() - 91) < 10) )

        ## b-tagged AK4 jet selection from Table 4 in Ch. 4 of AN2023_047_v4
        nBJet = 0
        btagWPM = 0.2783 if YEAR == '2018' else (0.3040 if YEAR == '2017' else (0.2489 if 'post' in YEAR else 0.2598))
        for iJet in range(ch.nJet):
            if      ch.Jet_pt[iJet]   <= 30:  continue
            if  abs(ch.Jet_eta[iJet]) >= 2.4: continue
            if     ch.Jet_jetId[iJet] != 6:   continue
            if   ch.Jet_puId[iJet] < 4 and \
                    ch.Jet_pt[iJet]    <= 50: continue
            if ch.Jet_btagDeepFlavB[iJet] <= btagWPM: continue
            vJet = R.TLorentzVector()
            vJet.SetPtEtaPhiM(ch.Jet_pt[iJet], ch.Jet_eta[iJet],
                              ch.Jet_phi[iJet], ch.Jet_mass[iJet])
            if vJet.DeltaR(vFatH) <= 0.8: continue
            ovlpLep = False
            for vLep in (vSelMu+vSelEle):
                if vLep.DeltaR(vJet) <= 0.4:
                    ovlpLep = True
                    break
            if not ovlpLep:
                nBJet += 1
        ## End loop: for iJet in range(ch.nJet)
        sel['eq0b'] = (nBJet == 0)
        sel['eq1b'] = (nBJet == 1)
        sel['ge1b'] = (nBJet >= 1)
        sel['ge2b'] = (nBJet >= 2)

        vMET = R.TLorentzVector()
        vMET.SetPtEtaPhiM(ch.MET_pt, 0, ch.MET_phi, 0)
        
        if sel['1mu_sel'] or sel['1ele_sel']:
            if sel['eq0b'] and abs(vFatH.DeltaPhi(vLep+vMET)) > np.pi*0.75:
                sel['WlvLo'] = ((vLep+vMET).Pt() <  300)
                sel['WlvHi'] = ((vLep+vMET).Pt() >= 300)
            sel['ttblv']  = sel['eq1b']
            sel['ttbblv'] = sel['ge2b']
            for xCat in ['WlvLo','WlvHi','ttblv','ttbblv']:
                sel[xCat.replace('lv','mv')] = (sel[xCat] and sel['1mu_sel'])
                sel[xCat.replace('lv','ev')] = (sel[xCat] and sel['1ele_sel'])

        if sel_dilep_sel and sel['Zveto'] and sel['ge1b']:
            sel['ttbll'] = True
            sel['ttbmm'] = sel['mumu_sel']
            sel['ttbme'] = sel['muele_sel']
            sel['ttbem'] = sel['elemu_sel']
            sel['ttbee'] = sel['eleele_sel']
        if sel_dilep_sel and sel['Zsel']:
            sel['Zll'] = True
            sel['Zmm'] = sel['mumu_sel']
            sel['Zee'] = sel['eleele_sel']

        sel['X4bSB'] = (xFatH_X4b > 0.66)
        sel['X4bSR'] = (xFatH_X4b > 0.93)
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
