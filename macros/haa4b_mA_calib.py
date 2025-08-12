#! /usr/bin/env python
## Calibration study for mass(a) regression in terms of AK8 sub-jets
## Central NanoAOD branches, and custom branches from processor and postprocessor:
## https://cms-nanoaod-integration.web.cern.ch/autoDoc/NanoAODv9/2018UL/doc_SingleMuon_Run2018D-UL2018_MiniAODv2_NanoAODv9-v1.html
## https://gitlab.cern.ch/abrinke1/cmssw/-/blob/PNet_v2_2024_11_22_bkg/PhysicsTools/NanoAOD/python/jets_Hto4b_cff.py
## https://github.com/abrinke1/nanoAOD-tools/blob/nanoPostProc_SS/python/postprocessing/modules/haa4b/objectSelection.py
## https://github.com/abrinke1/nanoAOD-tools/blob/nanoPostProc_SS/python/postprocessing/modules/haa4b/genParticles.py

import os
import sys
import subprocess
import numpy as np
import ROOT as R

R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn

MAX_EVT = -1     ## Maximum number of events to process per MC sample
PRT_EVT = 10000  ## Print every Nth event while processing
DEBUG   = False

## Location of postprocessed input files
IN_DIRS  = '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/YEAR/MC/PNet_v2_2024_11_22/'
IN_FILES = 'SUSY_GluGluH_01J_HToAATo4B_Pt150_M-mA_TuneCP5_13TeV_madgraph_pythia8/r1/PNet_v1_Skim_Lp_weight.root'

#YEARS  = ['2016APV','2016','2017','2018']
YEARS  = ['2018']
#MASSES = ['12']+[str(15+5*i) for i in range(10)]
MASSES = ['12']+[str(15+5*i) for i in range(5)]  ## Through 35 GeV
#MASSES = ['15']

## Histograms for mass(a) distributions
hst = {}
## Count events passing different selection cuts
count = {}
## Color codes for histograms
color = {}
color['Med_Med'] = R.kBlack
color['Hi_Hi']   = R.kBlue
color['Lo_Lo']   = R.kRed
color['Hi_Med']  = R.kCyan-3
color['Med_Hi']  = R.kViolet
color['Lo_Med']  = R.kMagenta
color['Med_Lo']  = R.kGreen+1
color['Hi_Lo']   = R.kOrange+8
color['Lo_Hi']   = R.kTeal+7

if not os.path.exists('plots/png/haa4b_mA_calib'):
    os.makedirs('plots/png/haa4b_mA_calib')

for year in YEARS:
    print('\n\n*** Beginning to look at year %s ***\n' % year)
    hst[year] = {}
    count[year] = {}
    for mA in MASSES:
        for var in ['jms','jes','jrs','sjms','sjes','sjrs']:
            for idx in ['1','2']:
                ## 1D jet/sub-jet mass and energy and mass/pT scale histograms
                hst[year]['mA_%s_%s%s' % (mA, var, idx)] = R.TH1D('h_1D_mA_%s_%s%s_%s' % (mA, var, idx, year),
                                                                  'm(a) = %s, %s%s (%s)' % (mA, var, idx, year),
                                                                  60, 0, 3.0)
            for js1 in ['Hi','Med','Lo']:
                for js2 in ['Hi','Med','Lo']:
                    ## Categorize by overlapping AK4 jet ("j") or AK8 soft-drop sub-jet ("sj")
                    ## mass scale ("ms") or energy scale ("es") relative to GEN mass(a)
                    cat = '%s_%s1_%s_%s2_%s' % (mA, var, js1, var, js2)
                    hst[year]['mA_'+cat] = R.TH1D('h_1D_mA_'+cat+'_'+year,
                                                  'm(a) = %s, %s1 %s %s2 %s (%s)' % (mA, var, js1, var, js2, year),
                                                  int(mA)*5, 5, float(mA)+max(5, int(float(mA)/3)))
                    hst[year]['mA_'+cat].SetLineColor(color['%s_%s' % (js1, js2)])
    ## End loop: for js1 in ['Hi','Med','Lo']

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

    for key in ['None','iH','X4b','4B','iA','2j','2sj']:
        count[year][key] = 0
    for iEvt in range(nEntries):

        if iEvt >= MAX_EVT and MAX_EVT > 0: break
        if (iEvt % PRT_EVT) == 0: print('Looking at event #%d / %d' % (iEvt, nEntries))

        ch.GetEntry(iEvt)
        count[year]['None'] += 1

        ## Good Haa4b candidate AK8 must exist
        iH = ch.Haa4b_iFatH
        if iH < 0: continue
        if ch.FatJet_pt[iH] < 250.0: continue
        count[year]['iH'] += 1
        ## Higgs AK8 must pass X4b WP60
        if ch.FatJet_PNet_X4b_v2a_Haa4b_score[iH] + \
           ch.FatJet_PNet_X4b_v2b_Haa4b_score[iH] < 2.0*0.93: continue  ## Higgs signal must pass WP60
        count[year]['X4b'] += 1
        ## Higgs AK8 must contain 4 GEN-level b-hadrons and quarks
        if ch.FatJet_nBHadrons[iH] < 4: continue
        if ch.Haa4b_FatH_nBQuarks  < 4: continue
        count[year]['4B'] += 1
        ## Event must have good GEN "a" bosons
        iA1 = ch.GEN_a1_idx
        iA2 = ch.GEN_a2_idx
        if iA1 < 0 or iA2 < 0: continue
        count[year]['iA'] += 1
        mA1  = int(ch.GenPart_mass[iA1])
        mA2  = int(ch.GenPart_mass[iA2])
        ptA1 = ch.GenPart_pt[iA1]
        ptA2 = ch.GenPart_pt[iA2]
        assert(mA1 == mA2)
        assert(str(mA1) in MASSES)

        ## Look for two highest-pT AK4 jets overlapping Higgs AK8
        jmass = []
        jpt   = []
        for ij in range(ch.nJet):
            if ch.Jet_Haa4b_ovlp_iFat[ij] != iH: continue
            if ch.Jet_Haa4b_presel[ij]    != 1:  continue
            jmass.append(ch.Jet_mass[ij])
            jpt  .append(ch.Jet_pt[ij])
            if len(jmass) == 2: break

        ## Look for Higgs AK8 sub-jets
        sjmass = []
        sjpt   = []
        if ch.FatJet_subJetIdx1[iH] >= 0:
            sj1 = ch.FatJet_subJetIdx1[iH]
            sjmass.append(ch.SubJet_mass[sj1])
            sjpt  .append(ch.SubJet_pt[sj1])
        if ch.FatJet_subJetIdx2[iH] >= 0:
            sj2 = ch.FatJet_subJetIdx2[iH]
            sjmass.append(ch.SubJet_mass[sj2])
            sjpt  .append(ch.SubJet_pt[sj2])

        ## Fill histograms based on relation between AK4 jet RECO and GEN mass/pT
        if len(jmass) == 2:
            count[year]['2j'] += 1
            jms1r = jmass[0] / mA1
            jms2r = jmass[1] / mA2
            jes1r = jpt[0] / ptA1
            jes2r = jpt[1] / ptA2
            jrs1r = jms1r / jes1r
            jrs2r = jms2r / jes2r
            jms1 = 'Hi' if jms1r > 1.80 else 'Lo' if jms1r < 1.20 else 'Med'
            jms2 = 'Hi' if jms2r > 1.25 else 'Lo' if jms2r < 0.75 else 'Med'
            jes1 = 'Hi' if jes1r > 1.05 else 'Lo' if jes1r < 0.85 else 'Med'
            jes2 = 'Hi' if jes2r > 1.00 else 'Lo' if jes2r < 0.75 else 'Med'
            jrs1 = 'Hi' if jrs1r > 1.90 else 'Lo' if jrs1r < 1.30 else 'Med'
            jrs2 = 'Hi' if jrs2r > 1.35 else 'Lo' if jrs2r < 0.95 else 'Med'
            cat_jes = '%d_jes1_%s_jes2_%s' % (mA1, jes1, jes2)
            cat_jms = '%d_jms1_%s_jms2_%s' % (mA1, jms1, jms2)
            cat_jrs = '%d_jrs1_%s_jrs2_%s' % (mA1, jrs1, jrs2)
            hst[year]['mA_'+cat_jes].Fill(ch.FatJet_PNet_34massAa[iH])
            hst[year]['mA_'+cat_jms].Fill(ch.FatJet_PNet_34massAa[iH])
            hst[year]['mA_'+cat_jrs].Fill(ch.FatJet_PNet_34massAa[iH])
            hst[year]['mA_%s_jms1' % mA1].Fill(jms1r)
            hst[year]['mA_%s_jms2' % mA1].Fill(jms2r)
            hst[year]['mA_%s_jes1' % mA1].Fill(jes1r)
            hst[year]['mA_%s_jes2' % mA1].Fill(jes2r)
            hst[year]['mA_%s_jrs1' % mA1].Fill(jrs1r)
            hst[year]['mA_%s_jrs2' % mA1].Fill(jrs2r)

            
        ## Fill histograms based on relation between subjet RECO and GEN mass/pT
        if len(sjmass) == 2:
            count[year]['2sj'] += 1
            sjms1r = sjmass[0] / mA1
            sjms2r = sjmass[1] / mA2
            sjes1r = sjpt[0] / ptA1
            sjes2r = sjpt[1] / ptA2
            sjrs1r = sjms1r / sjes1r
            sjrs2r = sjms2r / sjes2r
            sjms1 = 'Hi' if sjms1r > 2.00 else 'Lo' if sjms1r < 1.20 else 'Med'
            sjms2 = 'Hi' if sjms2r > 1.35 else 'Lo' if sjms2r < 0.80 else 'Med'
            sjes1 = 'Hi' if sjes1r > 1.05 else 'Lo' if sjes1r < 0.85 else 'Med'
            sjes2 = 'Hi' if sjes2r > 1.05 else 'Lo' if sjes2r < 0.80 else 'Med'
            sjrs1 = 'Hi' if sjrs1r > 2.00 else 'Lo' if sjrs1r < 1.25 else 'Med'
            sjrs2 = 'Hi' if sjrs2r > 1.40 else 'Lo' if sjrs2r < 0.95 else 'Med'
            cat_sjes = '%d_sjes1_%s_sjes2_%s' % (mA1, sjes1, sjes2)
            cat_sjms = '%d_sjms1_%s_sjms2_%s' % (mA1, sjms1, sjms2)
            cat_sjrs = '%d_sjrs1_%s_sjrs2_%s' % (mA1, sjrs1, sjrs2)
            hst[year]['mA_'+cat_sjes].Fill(ch.FatJet_PNet_34massAa[iH])
            hst[year]['mA_'+cat_sjms].Fill(ch.FatJet_PNet_34massAa[iH])
            hst[year]['mA_'+cat_sjrs].Fill(ch.FatJet_PNet_34massAa[iH])
            hst[year]['mA_%s_sjms1' % mA1].Fill(sjms1r)
            hst[year]['mA_%s_sjms2' % mA1].Fill(sjms2r)
            hst[year]['mA_%s_sjes1' % mA1].Fill(sjes1r)
            hst[year]['mA_%s_sjes2' % mA1].Fill(sjes2r)
            hst[year]['mA_%s_sjrs1' % mA1].Fill(sjrs1r)
            hst[year]['mA_%s_sjrs2' % mA1].Fill(sjrs2r)


    ## End loop: for iEvt in range(nEntries)
## End loop: for year in YEARS

print('\n*** Finished looking at all years! ***\n\n')

for year in YEARS:
    nTot = count[year]['None']
    nPrev = nTot
    print('\n%s had %d total events' % (year, nTot))
    for key in ['iH','X4b','4B','iA','2j','2sj']:
        nX = count[year][key]
        print('  * %s cut: %d (%.2f%% lost, %.2f%% total)' % (key, nX, 100*(nPrev-nX)/nPrev, 100*nX/nTot))
        nPrev = nX if key != '2j' else nPrev
## End loop: for year in YEARS

print('\nWriting out to plots/haa4b_mA_calib.root')
outf = R.TFile('plots/haa4b_mA_calib.root', 'recreate')
outf.cd()
cans = {}
maxs = {}
for year in YEARS:
    cans[year] = {}
    maxs[year] = {}
    for key in hst[year].keys():
        hist = hst[year][key]
        keyA = key.replace('Hi','All').replace('Med','All').replace('Lo','All')
        hist.SetLineWidth(2)
        hist.Write()
        R.gStyle.SetOptStat(1110) ## Display rms, mean, and number of entries
        can = R.TCanvas(hist.GetName())
        can.cd()
        hist.Draw('histe')
        can.SaveAs('plots/png/haa4b_mA_calib/'+hist.GetName()+'.png')
        del can

        R.gStyle.SetOptStat(0) ## Don't display stat boxes
        if 'All' in keyA and not keyA in cans[year].keys():
            maxs[year][keyA] = hist.GetMaximum()
            cans[year][keyA] = R.TCanvas('%s_%s' % (keyA, year))
            cans[year][keyA].cd()
            hist.Draw('histe')
        elif 'All' in keyA:
            maxs[year][keyA] = max(maxs[year][keyA], hist.GetMaximum())
            cans[year][keyA].cd()
            hist.Draw('histesame')
    ## End loop: for key in hst[year].keys()

    R.gStyle.SetOptStat(0) ## Don't display stat boxes
    for keyA in cans[year].keys():
        cans[year][keyA].cd()
        leg = R.TLegend(0.12, 0.42, 0.48, 0.88)
        for js1 in ['Hi','Med','Lo']:
            for js2 in ['Hi','Med','Lo']:
                key = keyA.replace('1_All','1_'+js1).replace('2_All','2_'+js2)
                hist = hst[year][key]
                hist.SetTitle(hist.GetTitle().replace('1 '+js1,'1 All').replace('2 '+js2,'2 All'))
                hist.GetYaxis().SetRangeUser(0, 1.2*maxs[year][keyA])
                leg.AddEntry(hist, '%.2f #pm %.2f (%3s-%3s)' % (hist.GetMean(), hist.GetStdDev(), js1, js2), 'l')
        leg.Draw()
        cans[year][keyA].SaveAs('plots/png/haa4b_mA_calib/%s_%s.png' % (keyA, year))
        del leg

## End loop: for year in YEARS
        
outf.Write()
outf.Close()

print('\n*** All done!!! ***\n\n')

