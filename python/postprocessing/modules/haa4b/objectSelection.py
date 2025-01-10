## H --> aa --> 4b search with "boosted" decays to AK8 jets
## Select object candidates and "clean" between types of objects

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import ROOT
import os
import sys
import math
import random
ROOT.PyConfig.IgnoreCommandLineOptions = True

class Haa4bObjectSelectionProducer(Module):

    def __init__(self, IS_MC, YEAR):
        self.isMC     = IS_MC
        self.year_str = YEAR
        ## Latest object selection information:
        ## https://indico.cern.ch/event/1430644/#2-update-on-higgs-aa-4b-booste
        ## https://indico.cern.ch/event/1450818/#17-andrew-brinkerhoff
        ## https://indico.cern.ch/event/1424480/#17-andrew-brinkerhoff (HEM info)
        self.Fat_sel_ID = 6
        self.Fat_sel_pt = 170.0
        self.Fat_sel_eta = 2.4
        self.Fat_sel_msoft = 20.0
        self.Fat_candH_Xbb = 0.75
        ## Loose tagger WPs, see "ParticleNet top and W tag working points for UL" e-mail
        ## https://twiki.cern.ch/twiki/bin/viewauth/CMS/ParticleNetSFs#Working_Points
        self.Fat_candX_Vqq = 0.9
        self.Fat_candX_Tbqq = 0.5

        self.Mu_sel_pt = 10.0
        self.Mu_sel_eta = 2.4
        self.Mu_sel_iso = 0.10
        self.Mu_sel_dxy = 0.02
        self.Mu_sel_dz = 0.10
        #self.Mu_sel_ID = (mediumPromptId or (pt > 53 and highPtId))
        self.Mu_trig_pt = 26.0

        self.Ele_sel_pt = 10.0
        self.Ele_sel_eta = 2.5
        self.Ele_sel_crack = [1.44, 1.57]
        #self.Ele_sel_ID = (mvaFall17V2Iso_WPL and (mvaFall17V2Iso_WP90 or (pt > 35 and cutBased_HEEP)
        self.Ele_trig_pt = 35.0
        self.Ele_trig_dxy = 0.02
        self.Ele_trig_dz = 0.10
        #self.Ele_trig_ID = (mvaFall17V2Iso_WP90 and (mvaFall17V2Iso_WP80 or cutBased_HEEP))

        self.Jet_presel_pt = 15.0
        self.Jet_presel_ID = 6
        #self.Jet_presel_puId = (pt > 50 or puId >= 4)
        self.Jet_sel_pt = 30.0
        self.bJet_sel_eta = 2.4
        self.bJet_sel_flavB = 0.2783
        self.bJet_candH_dR = 1.4


    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        ## Critical to use "O" and not "b" or "B" for bools!!!
        ## https://root.cern/doc/master/TBranch_8cxx_source.html
        ## In lines 365 - 371, we see that the type for "O" is TLeafO (bool type),
        ## but "B" and "b" are TLeafB (8 bit integer), and for some reason PyROOT
        ## chokes on the bool type:
        ## https://root-forum.cern.ch/t/reading-boolean-data-from-ttree-using-pyroot/39178
        ## For further info see:
        ## https://docs.python.org/3/library/array.html
        ## python/postprocessing/framework/output.py _rootBranchType2PythonArray and L26 and L27

        ## Branches for object selection, per object
        self.out.branch("FatJet_HEM",             "O", lenVar="nFatJet")
        self.out.branch("FatJet_Haa4b_sel",       "O", lenVar="nFatJet")
        self.out.branch("FatJet_Haa4b_candH",     "O", lenVar="nFatJet")
        self.out.branch("FatJet_Haa4b_candX",     "O", lenVar="nFatJet")
        self.out.branch("FatJet_Haa4b_ovlp_iMu",  "I", lenVar="nFatJet")
        self.out.branch("FatJet_Haa4b_ovlp_iEle", "I", lenVar="nFatJet")
        if self.isMC:
            self.out.branch("FatJet_genPartIdx",      "I", lenVar="nFatJet")
            self.out.branch("FatJet_genPart_pdgId",   "I", lenVar="nFatJet")
            self.out.branch("FatJet_genPart_dR",      "F", lenVar="nFatJet")
            self.out.branch("FatJet_nBPartons",       "I", lenVar="nFatJet")
            self.out.branch("FatJet_nCPartons",       "I", lenVar="nFatJet")
        self.out.branch("Muon_Haa4b_sel",       "O", lenVar="nMuon")
        self.out.branch("Muon_Haa4b_trig",      "O", lenVar="nMuon")
        self.out.branch("Muon_Haa4b_ovlp_iFat", "I", lenVar="nMuon")
        self.out.branch("Electron_HEM",             "O", lenVar="nElectron")
        self.out.branch("Electron_Haa4b_sel",       "O", lenVar="nElectron")
        self.out.branch("Electron_Haa4b_trig",      "O", lenVar="nElectron")
        self.out.branch("Electron_Haa4b_ovlp_iFat", "I", lenVar="nElectron")
        self.out.branch("Jet_HEM",             "O", lenVar="nJet")
        self.out.branch("Jet_HEM_EM",          "O", lenVar="nJet")
        self.out.branch("Jet_Haa4b_presel",    "O", lenVar="nJet")
        self.out.branch("Jet_Haa4b_sel",       "O", lenVar="nJet")
        self.out.branch("Jet_Haa4b_btag",      "O", lenVar="nJet")
        self.out.branch("Jet_Haa4b_candH",     "O", lenVar="nJet")
        self.out.branch("Jet_Haa4b_ovlp_iFat", "I", lenVar="nJet")
        self.out.branch("Jet_Haa4b_ovlp_iMu",  "I", lenVar="nJet")
        self.out.branch("Jet_Haa4b_ovlp_iEle", "I", lenVar="nJet")
        if self.isMC:
            self.out.branch("Jet_genPartIdx",      "I", lenVar="nJet")
            self.out.branch("Jet_genPart_pdgId",   "I", lenVar="nJet")
            self.out.branch("Jet_genPart_dR",      "F", lenVar="nJet")

        ## Branches for counting selected objects, per event
        self.out.branch("Haa4b_nFatSel",      "I")
        self.out.branch("Haa4b_nFatH",        "I")
        self.out.branch("Haa4b_nFatHj",       "I")
        self.out.branch("Haa4b_nFatX",        "I")
        self.out.branch("Haa4b_nMuSel",       "I")
        self.out.branch("Haa4b_nMuTrig",      "I")
        self.out.branch("Haa4b_nEleSel",      "I")
        self.out.branch("Haa4b_nEleTrig",     "I")
        self.out.branch("Haa4b_nJetSel",      "I")
        self.out.branch("Haa4b_nJetSelNonH",  "I")
        self.out.branch("Haa4b_nJetBtag",     "I")
        self.out.branch("Haa4b_nJetBtagNonH", "I")

        ## Branches for preliminary category designations
        CATS = ['gg0l','VBFjj','Vjj','ttHad','Zvv','Zll','Wlv','ttlv','ttll','2lSS','3l','other']
        for cat in CATS:
            self.out.branch("Haa4b_cat_%s" % cat, "O")
            self.out.branch("Haa4b_HEM_%s" % cat, "O")
        self.out.branch("Haa4b_cat_idx", "I")
        self.out.branch("Haa4b_nCats",   "I")
        
        ## Branches with special or multi-object quantities
        self.out.branch("Haa4b_iFatH",         "I")
        self.out.branch("Haa4b_iFatHj_jet",    "I")
        self.out.branch("Haa4b_iFatX1",        "I")
        self.out.branch("Haa4b_iFatX2",        "I")
        self.out.branch("Haa4b_iLep1",         "I")
        self.out.branch("Haa4b_iLep2",         "I")
        self.out.branch("Haa4b_Lep1_pdgId",    "I")
        self.out.branch("Haa4b_Lep2_pdgId",    "I")
        self.out.branch("Haa4b_iJetSel1",      "I")
        self.out.branch("Haa4b_iJetSel2",      "I")
        self.out.branch("Haa4b_iJetSelNonH1",  "I")
        self.out.branch("Haa4b_iJetSelNonH2",  "I")
        self.out.branch("Haa4b_iJetBtag1",     "I")
        self.out.branch("Haa4b_iJetBtag2",     "I")
        self.out.branch("Haa4b_iJetBtagNonH1", "I")
        self.out.branch("Haa4b_iJetBtagNonH2", "I")

        self.out.branch("Haa4b_FatH_mass",  "F")
        self.out.branch("Haa4b_FatHj_mass", "F")
        self.out.branch("Haa4b_FatX1_mass", "F")
        self.out.branch("Haa4b_FatX2_mass", "F")
        
        self.out.branch("Haa4b_FatH_msoftdrop",  "F")
        self.out.branch("Haa4b_FatHj_msoftdrop", "F")
        self.out.branch("Haa4b_FatX1_msoftdrop", "F")
        self.out.branch("Haa4b_FatX2_msoftdrop", "F")

        self.out.branch("Haa4b_FatH_pt",  "F")
        self.out.branch("Haa4b_FatHj_pt", "F")
        self.out.branch("Haa4b_FatX1_pt", "F")
        self.out.branch("Haa4b_FatX2_pt", "F")

        self.out.branch("Haa4b_FatHj_dR",     "F")
        self.out.branch("Haa4b_FatHj_jet_pt", "F")
        
        self.out.branch("Haa4b_FatH_MET_dPhi",    "F")
        self.out.branch("Haa4b_FatH_lepMET_dPhi", "F")
        self.out.branch("Haa4b_FatH_dilep_dPhi",  "F")

        if self.isMC:
            self.out.branch("Haa4b_FatH_nBQuarks",   "I")
            self.out.branch("Haa4b_FatH_nCQuarks",   "I")
            self.out.branch("Haa4b_FatH_nLFQuarks",  "I")
            self.out.branch("Haa4b_FatHj_jet_nBQuarks",  "I")
            self.out.branch("Haa4b_FatHj_jet_nCQuarks",  "I")
            self.out.branch("Haa4b_FatHj_jet_nLFQuarks", "I")
            self.out.branch("Haa4b_FatX1_nBQuarks",  "I")
            self.out.branch("Haa4b_FatX1_nCQuarks",  "I")
            self.out.branch("Haa4b_FatX1_nLFQuarks", "I")
            self.out.branch("Haa4b_FatX2_nBQuarks",  "I")
            self.out.branch("Haa4b_FatX2_nCQuarks",  "I")
            self.out.branch("Haa4b_FatX2_nLFQuarks", "I")

        self.out.branch("Haa4b_FatH_tagHaa34b_v1",  "F")
        self.out.branch("Haa4b_FatH_tagHaa34b_v2a", "F")
        self.out.branch("Haa4b_FatH_tagHaa34b_v2b", "F")
        self.out.branch("Haa4b_FatH_tagHaa4b_v1",  "F")
        self.out.branch("Haa4b_FatH_tagHaa4b_v2a", "F")
        self.out.branch("Haa4b_FatH_tagHaa4b_v2b", "F")

        self.out.branch("Haa4b_FatX_tagWZ_max",  "F")
        self.out.branch("Haa4b_FatX_tagTop_max", "F")
        self.out.branch("Haa4b_FatX1_tagWZ",     "F")
        self.out.branch("Haa4b_FatX2_tagWZ",     "F")
        self.out.branch("Haa4b_FatX1_tagTop",    "F")
        self.out.branch("Haa4b_FatX2_tagTop",    "F")

        self.out.branch("Haa4b_dijet_mass",   "F")
        self.out.branch("Haa4b_dijet_pt",     "F")
        self.out.branch("Haa4b_dijet_dEta",   "F")
        self.out.branch("Haa4b_dijet_dPhi",   "F")
        self.out.branch("Haa4b_dilep_mass",   "F")
        self.out.branch("Haa4b_dilep_pt",     "F")
        self.out.branch("Haa4b_dilep_dR",     "F")
        self.out.branch("Haa4b_dilep_dPhi",   "F")
        self.out.branch("Haa4b_dilep_charge", "I")
        self.out.branch("Haa4b_dilep_SFOS",   "O")
        self.out.branch("Haa4b_lepMET_pt",    "F")
        self.out.branch("Haa4b_lepMET_MT",    "F")
        self.out.branch("Haa4b_lepMET_dPhi",  "F")

        ## Branches for categorization and triggers, per event
        self.out.branch("year",   "I")
        self.out.branch("isData", "O")
        self.out.branch("isMC",   "O")
        self.out.branch("Haa4b_isHad",      "O")
        self.out.branch("Haa4b_isLep",      "O")
        self.out.branch("Haa4b_isMu",       "O")
        self.out.branch("Haa4b_isEle",      "O")
        self.out.branch("Haa4b_trigFat",    "O")
        self.out.branch("Haa4b_trigBtag",   "O")
        self.out.branch("Haa4b_trigVBF",    "O")
        self.out.branch("Haa4b_trigMET",    "O")
        self.out.branch("Haa4b_trigMu",     "O")
        self.out.branch("Haa4b_trigEle",    "O")
        self.out.branch("Haa4b_passFilters","O")
        

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        year      = int(self.year_str)
        data_year = (2016*(event.run >= 264283) +
                     0001*(event.run >= 293953) +
                     0001*(event.run >= 312300) +
                     0001*(event.run >  327564))
        if (not self.isMC) and year != data_year:
            print('\nWeird error! Nominal year %d, but run = %d implies year = %d.' % (year, event.run, data_year))


        CATS = ['gg0l','VBFjj','Vjj','ttHad','Zvv','Zll','Wlv','ttlv','ttll','2lSS','3l','other']

        ################################
        ## Set vectors for GEN particles
        if self.isMC:
            iGP = {}
            iGP['H'] = event.GEN_H_idx
            iGP['X1'] = event.GEN_X1_idx
            iGP['X2'] = event.GEN_X2_idx
            iGP['a1'] = event.GEN_a1_idx
            iGP['a2'] = event.GEN_a2_idx
            iGP['b11'] = event.GEN_b11_idx
            iGP['b12'] = event.GEN_b12_idx
            iGP['b21'] = event.GEN_b21_idx
            iGP['b22'] = event.GEN_b22_idx
            iGP['d11'] = event.GEN_d11_idx
            iGP['d12'] = event.GEN_d12_idx
            iGP['d13'] = event.GEN_d13_idx
            iGP['d21'] = event.GEN_d21_idx
            iGP['d22'] = event.GEN_d22_idx
            iGP['d23'] = event.GEN_d23_idx
            GPks = iGP.keys()
            GPqs = [gpk for gpk in GPks if len(gpk) == 3]
            has_GPs = (sum([(iGP[gpk] >= 0) for gpk in GPks]) > 0)
            if has_GPs:
                GPs = Collection(event, "GenPart")
            vGP = {}
            idGP = {}
            for gpk in GPks:
                vGP[gpk] = ROOT.TLorentzVector()
                idGP[gpk] = -99
                if iGP[gpk] >= 0:
                    gp = GPs[iGP[gpk]]
                    vGP[gpk].SetPtEtaPhiM(gp.pt, gp.eta, gp.phi, gp.mass)
                    idGP[gpk] = gp.pdgId
        

        #######################
        ## Set bits for FatJets
        FatJets = Collection(event, "FatJet")
        nFatJets = len(FatJets)
        Fat_HEM_bits    = [0] * nFatJets
        Fat_sel_bits    = [0] * nFatJets
        Fat_candH_bits  = [0] * nFatJets
        Fat_candX_bits  = [0] * nFatJets
        Fat_ovlp_iMu   = [-99] * nFatJets
        Fat_ovlp_iEle  = [-99] * nFatJets
        Fat_genPartIdx = [-99] * nFatJets
        if self.isMC:
            Fat_GP_pdgId  = [-99] * nFatJets
            Fat_GP_dR     = [-99] * nFatJets
            Fat_nBPartons = [-99] * nFatJets
            Fat_nCPartons = [-99] * nFatJets

        max_Haa_score = -99.0
        max_WZ_tag  = -99.0
        max_top_tag = -99.0
        iCandH   = -99
        vCandH   = ROOT.TLorentzVector()
        vCandHsd = ROOT.TLorentzVector()
        iCandsX = []
        vCandsX = []
        iFat = -1
        for fj in FatJets:
            iFat += 1
            vFat = ROOT.TLorentzVector()
            vFat.SetPtEtaPhiM(fj.pt, fj.eta, fj.phi, fj.mass)
            ## Match FatJets to heavy GEN particles
            if self.isMC:
                iGF = fj.genJetAK8Idx
                Fat_nBPartons[iFat] = (event.GenJetAK8_nBPartons[iGF] if iGF >= 0 else -99)
                Fat_nCPartons[iFat] = (event.GenJetAK8_nCPartons[iGF] if iGF >= 0 else -99)
                min_GP_dR = 0.8
                for gpk in ['H','X1','X2']:
                    if iGP[gpk] >= 0 and vGP[gpk].DeltaR(vFat) < min_GP_dR:
                        min_GP_dR = vGP[gpk].DeltaR(vFat)
                        Fat_genPartIdx[iFat] = iGP[gpk]
                        Fat_GP_pdgId[iFat]   = idGP[gpk]
                        Fat_GP_dR[iFat]      = vGP[gpk].DeltaR(vFat)
            ## Flag FatJets in the HEM region
            if abs(fj.phi + 1.22) < 0.55 and fj.eta < -1.1:
                Fat_HEM_bits[iFat] = 1
            ## Common FatJet selection
            if fj.jetId     <  self.Fat_sel_ID:    continue
            if fj.pt        <= self.Fat_sel_pt:    continue
            if abs(fj.eta)  >= self.Fat_sel_eta:   continue
            if fj.msoftdrop <= self.Fat_sel_msoft: continue
            Fat_sel_bits[iFat] = 1
            ## Loose tagging cuts for W/Z/top candidates
            if fj.particleNet_WZvsQCD > self.Fat_candX_Vqq or fj.particleNet_TvsQCD > self.Fat_candX_Tbqq:
                Fat_candX_bits[iFat] = 1
                iCandsX.append(iFat)
                vCandX = ROOT.TLorentzVector()
                vCandX.SetPtEtaPhiM(fj.pt, fj.eta, fj.phi, fj.mass)
                vCandsX.append(vCandX)
            ## Xbb selection for Haa4b candidates
            if fj.particleNetMD_XbbvsQCD <= self.Fat_candH_Xbb: continue
            Haa_score = fj.PNet_X4b_v2a_Haa34b_score + fj.PNet_X4b_v2b_Haa34b_score
            if Haa_score > max_Haa_score:
                max_Haa_score = Haa_score
                iCandH = iFat
                vCandH.SetPtEtaPhiM(fj.pt, fj.eta, fj.phi, fj.mass)
                vCandHsd.SetPtEtaPhiM(fj.pt, fj.eta, fj.phi, fj.msoftdrop)
        ## End loop: for fj in FatJets

        ## Set the single Higgs candidate, exclude from W/Z/top candidates ("X")
        if iCandH >= 0:
            Fat_candH_bits[iCandH] = 1
            Fat_candX_bits[iCandH] = 0
            for iX in reversed(range(len(iCandsX))):
                if iCandsX[iX] == iCandH:
                    del iCandsX[iX]
                    del vCandsX[iX]

        self.out.fillBranch("FatJet_HEM",           Fat_HEM_bits)
        self.out.fillBranch("FatJet_Haa4b_sel",     Fat_sel_bits)
        self.out.fillBranch("FatJet_Haa4b_candH",   Fat_candH_bits)
        if self.isMC:
            self.out.fillBranch("FatJet_genPartIdx",    Fat_genPartIdx)
            self.out.fillBranch("FatJet_genPart_pdgId", Fat_GP_pdgId)
            self.out.fillBranch("FatJet_genPart_dR",    Fat_GP_dR)
            self.out.fillBranch("FatJet_nBPartons",     Fat_nBPartons)
            self.out.fillBranch("FatJet_nCPartons",     Fat_nCPartons)
        self.out.fillBranch("Haa4b_nFatSel", sum(Fat_sel_bits))
        self.out.fillBranch("Haa4b_nFatH",   sum(Fat_candH_bits))


        #####################
        ## Set bits for Muons
        Muons = Collection(event, "Muon")
        nMuons = len(Muons)
        Mu_sel_bits  = [0] * nMuons
        Mu_trig_bits = [0] * nMuons
        Mu_ovlp_iFat = [-99] * nMuons

        iSelMus = []
        vSelMus = []
        iTrigMus = []
        vTrigMus = []
        iMu = -1
        for mu in Muons:
            iMu += 1
            ## Check for overlap with candidate FatJets
            vMu = ROOT.TLorentzVector()
            vMu.SetPtEtaPhiM(mu.pt, mu.eta, mu.phi, 0.10566)
            if iCandH >= 0 and vMu.DeltaR(vCandH) < 0.8:
                Mu_ovlp_iFat[iMu] = iCandH
                if Fat_ovlp_iMu[iCandH] < 0:
                    Fat_ovlp_iMu[iCandH] = iMu
            for iX in range(len(iCandsX)):
                if vMu.DeltaR(vCandsX[iX]) < 0.8:
                    if Mu_ovlp_iFat[iMu] < 0:
                        Mu_ovlp_iFat[iMu] = iCandsX[iX]
                    if Fat_ovlp_iMu[iCandsX[iX]] < 0:
                        Fat_ovlp_iMu[iCandsX[iX]] = iMu
            ## Common muon selection
            if mu.pt               <= self.Mu_sel_pt:  continue
            if abs(mu.eta)         >= self.Mu_sel_eta: continue
            if mu.miniPFRelIso_all >= self.Mu_sel_iso: continue
            if abs(mu.dxy)         >= self.Mu_sel_dxy: continue
            if abs(mu.dz)          >= self.Mu_sel_dz:  continue
            if not ( mu.mediumPromptId >= 1 or
                     (mu.pt > 200 and mu.highPtId >= 1) ): continue
            Mu_sel_bits[iMu] = 1
            iSelMus.append(iMu)
            vSelMus.append(vMu)
            ## Triggerable muon selection
            if mu.pt <= self.Mu_trig_pt: continue
            Mu_trig_bits[iMu] = 1
            iTrigMus.append(iMu)
            vTrigMus.append(vMu)
        ## End loop: for mu in Muons

        self.out.fillBranch("Muon_Haa4b_sel",  Mu_sel_bits)
        self.out.fillBranch("Muon_Haa4b_trig", Mu_trig_bits)
        self.out.fillBranch("Muon_Haa4b_ovlp_iFat", Mu_ovlp_iFat)
        self.out.fillBranch("FatJet_Haa4b_ovlp_iMu", Fat_ovlp_iMu)
        self.out.fillBranch("Haa4b_nMuSel",  sum(Mu_sel_bits))
        self.out.fillBranch("Haa4b_nMuTrig", sum(Mu_trig_bits))


        #########################
        ## Set bits for Electrons
        Electrons = Collection(event, "Electron")
        nElectrons = len(Electrons)
        Ele_HEM_bits  = list(0 for ele in Electrons)
        Ele_sel_bits  = list(0 for ele in Electrons)
        Ele_trig_bits = list(0 for ele in Electrons)
        Ele_ovlp_iFat = list(-99 for ele in Electrons)

        iSelEles = []
        vSelEles = []
        iTrigEles = []
        vTrigEles = []
        iEle = -1
        for ele in Electrons:
            iEle += 1
            ## Flag Electrons in the HEM region
            if abs(ele.phi + 1.22) < 0.45 and ele.eta < -1.2:
                Ele_HEM_bits[iEle] = 1
            ## Check for overlap with candidate FatJets
            vEle = ROOT.TLorentzVector()
            vEle.SetPtEtaPhiM(ele.pt, ele.eta, ele.phi, 0.000511)
            if iCandH >= 0 and vEle.DeltaR(vCandH) < 0.8:
                Ele_ovlp_iFat[iEle] = iCandH
                if Fat_ovlp_iEle[iCandH] < 0:
                    Fat_ovlp_iEle[iCandH] = iEle
            for iX in range(len(iCandsX)):
                if vEle.DeltaR(vCandsX[iX]) < 0.8:
                    if Ele_ovlp_iFat[iEle] < 0:
                        Ele_ovlp_iFat[iEle] = iCandsX[iX]
                    if Fat_ovlp_iEle[iCandsX[iX]] < 0:
                        Fat_ovlp_iEle[iCandsX[iX]] = iEle
            ## Common electron selection
            if ele.pt           <= self.Ele_sel_pt:  continue
            if abs(ele.eta)     >= self.Ele_sel_eta: continue
            if (abs(ele.eta)    >= self.Ele_sel_crack[0] and
                abs(ele.eta)    <= self.Ele_sel_crack[1]): continue
            if ele.mvaFall17V2Iso_WPL < 1: continue
            if not ( ele.mvaFall17V2Iso_WP90 >= 1 or
                     (ele.pt > 35 and ele.cutBased_HEEP >= 1) ): continue
            Ele_sel_bits[iEle] = 1
            iSelEles.append(iEle)
            vSelEles.append(vEle)
            ## Triggerable electron selection
            if ele.pt       <= self.Ele_trig_pt:  continue
            if abs(ele.dxy) >= self.Ele_trig_dxy: continue
            if abs(ele.dz)  >= self.Ele_trig_dz:  continue
            if ele.mvaFall17V2Iso_WP90 < 1:   continue
            if not ( ele.mvaFall17V2Iso_WP80 >= 1 or
                     ele.cutBased_HEEP >= 1): continue
            Ele_trig_bits[iEle] = 1
            iTrigEles.append(iEle)
            vTrigEles.append(vEle)
        ## End loop: for ele in Electrons

        ## Get selected electrons in HEM region
        Ele_sel_HEM_bits = [Ele_sel_bits[ii]*Ele_HEM_bits[ii] for ii in range(len(Ele_sel_bits))]

        self.out.fillBranch("Electron_HEM",        Ele_HEM_bits)
        self.out.fillBranch("Electron_Haa4b_sel",  Ele_sel_bits)
        self.out.fillBranch("Electron_Haa4b_trig", Ele_trig_bits)
        self.out.fillBranch("Electron_Haa4b_ovlp_iFat", Ele_ovlp_iFat)
        self.out.fillBranch("FatJet_Haa4b_ovlp_iEle", Fat_ovlp_iEle)
        self.out.fillBranch("Haa4b_nEleSel",  sum(Ele_sel_bits))
        self.out.fillBranch("Haa4b_nEleTrig", sum(Ele_trig_bits))

        ## Events with at least one triggerable lepton belong to lepton categories
        nTrigLeps = len(iTrigMus) + len(iTrigEles)
        nSelLeps  = len(iSelMus)  + len(iSelEles)
        ## Find how many trigger and selected leptons overlap AK8 Higgs candidate
        nTrigLepsOvlpH = 0
        iSelLepsNonOvlpH  = []
        idSelLepsNonOvlpH = []
        vSelLepsNonOvlpH  = []
        ## Fill with triggerable leptons first, then secondary leptons
        for isTrig in [1,0]:
            for iMu in range(nMuons):
                if Mu_trig_bits[iMu] != isTrig: continue
                if iCandH >= 0 and Mu_ovlp_iFat[iMu] == iCandH:
                    nTrigLepsOvlpH += Mu_trig_bits[iMu]
                elif Mu_sel_bits[iMu] == 1:
                    iSelLepsNonOvlpH.append(iMu)
                    idSelLepsNonOvlpH.append(-13*Muons[iMu].charge)
                    vSelLepsNonOvlpH.append(vSelMus[iSelMus.index(iMu)])
            for iEle in range(nElectrons):
                if Ele_trig_bits[iEle] != isTrig: continue
                if iCandH >= 0 and Ele_ovlp_iFat[iEle] == iCandH:
                    nTrigLepsOvlpH += Ele_trig_bits[iEle]
                elif Ele_sel_bits[iEle] == 1:
                    iSelLepsNonOvlpH.append(iEle)
                    idSelLepsNonOvlpH.append(-11*Electrons[iEle].charge)
                    vSelLepsNonOvlpH.append(vSelEles[iSelEles.index(iEle)])
        ## End loop: for isTrig in [1,0]
        nSelLepsNonOvlpH = len(vSelLepsNonOvlpH)
        ## If FatJet overlaps trigger lepton, exclude from W/Z/top candidates ("X")
        for vTrigLep in vTrigMus+vTrigEles:
            for iX in reversed(range(len(iCandsX))):
                Fat_candX_bits[iCandsX[iX]] = 0
                if vCandsX[iX].DeltaR(vTrigLep) < 0.8:
                    del iCandsX[iX]
                    del vCandsX[iX]
        max_WZ_tag  = max([-99]+[FatJets[iX].particleNet_WZvsQCD for iX in iCandsX])
        max_top_tag = max([-99]+[FatJets[iX].particleNet_TvsQCD  for iX in iCandsX])
        self.out.fillBranch("FatJet_Haa4b_candX", Fat_candX_bits)
        self.out.fillBranch("Haa4b_nFatX",   sum(Fat_candX_bits))

        #####################
        ## Set bits for Jets
        Jets = Collection(event, "Jet")
        nJets = len(Jets)
        Jet_HEM_bits    = [0] * nJets
        Jet_HEM_EM_bits = [0] * nJets
        Jet_presel_bits = [0] * nJets
        Jet_sel_bits    = [0] * nJets
        Jet_btag_bits   = [0] * nJets
        Jet_candH_bits  = [0] * nJets
        Jet_ovlp_iFat = [-99] * nJets
        Jet_ovlp_iMu  = [-99] * nJets
        Jet_ovlp_iEle = [-99] * nJets
        if self.isMC:
            Jet_genPartIdx = [-99] * nJets
            Jet_GP_pdgId   = [-99] * nJets
            Jet_GP_dR      = [-99] * nJets

        vJets     = []
        iSelJets  = []
        vSelJets  = []
        iJetCandH = -99
        vJetCandH = ROOT.TLorentzVector()
        iJet = -1
        for jet in Jets:
            iJet += 1
            vJet = ROOT.TLorentzVector()
            vJet.SetPtEtaPhiM(jet.pt, jet.eta, jet.phi, jet.mass)
            vJets.append(vJet)
            ## Match Jets to decay products of heavy GEN particles
            if self.isMC:
                min_GP_dR = 0.4
                for gpq in GPqs:
                    if iGP[gpq] >= 0 and vGP[gpq].DeltaR(vJet) < min_GP_dR:
                        min_GP_dR = vGP[gpq].DeltaR(vJet)
                        Jet_genPartIdx[iJet] = iGP[gpq]
                        Jet_GP_pdgId[iJet]   = idGP[gpq]
                        Jet_GP_dR[iJet]      = vGP[gpq].DeltaR(vJet)
            ## Flag Jets in the HEM region
            if abs(jet.phi + 1.22) < 0.45 and jet.eta < -1.2 and jet.eta > -3.2:
                Jet_HEM_bits[iJet] = 1
            if abs(jet.phi + 1.22) < 0.55 and jet.eta < -1.1 and jet.eta > -3.2:
                if jet.pt*(jet.chEmEF + jet.neEmEF) > 10:
                    Jet_HEM_EM_bits[iJet] == 1
            ## Check for overlap with candidate FatJets and selected leptons
            if iCandH >= 0 and vJet.DeltaR(vCandH) < 0.8:
                Jet_ovlp_iFat[iJet] = iCandH
            else:
                for iX in range(len(iCandsX)):
                    if vJet.DeltaR(vCandsX[iX]) < 0.8:
                        Jet_ovlp_iFat[iJet] = iCandsX[iX]
                        break
            for sMu in range(len(iSelMus)):
                if vJet.DeltaR(vSelMus[sMu]) < 0.4:
                    Jet_ovlp_iMu[iJet] = iSelMus[sMu]
                    break
            for sEle in range(len(iSelEles)):
                if vJet.DeltaR(vSelEles[sEle]) < 0.4:
                    Jet_ovlp_iEle[iJet] = iSelEles[sEle]
                    break
            ## Common jet pre-selection
            if jet.pt    <= self.Jet_presel_pt: continue
            if jet.jetId <  self.Jet_presel_ID: continue
            if (jet.pt <= 50 and jet.puId < 4): continue
            Jet_presel_bits[iJet] = 1
            ## Reject all AK4 jets overlapping Higgs candidate AK8 jet
            if iCandH >= 0 and Jet_ovlp_iFat[iJet] == iCandH: continue
            ## If at least one triggerable lepton, reject jets overlapping selected leptons
            ## (Selcted leptons must not be overlapping Higgs candidate AK8 to "count")
            if nTrigLeps > 0:
                if Jet_ovlp_iMu[iJet]  >= 0 and (iCandH < 0 or Mu_ovlp_iFat [Jet_ovlp_iMu [iJet]] != iCandH): continue
                if Jet_ovlp_iEle[iJet] >= 0 and (iCandH < 0 or Ele_ovlp_iFat[Jet_ovlp_iEle[iJet]] != iCandH): continue
            ## "Regular" jet selection (i.e. not part of Haa4b candidate)
            if jet.pt > self.Jet_sel_pt:
                Jet_sel_bits[iJet] = 1
                iSelJets.append(iJet)
                vSelJets.append(vJet)
            ## b-tagged jet selection
            if abs(jet.eta)      >= self.bJet_sel_eta: continue
            if jet.btagDeepFlavB <= self.bJet_sel_flavB: continue
            if jet.pt > self.Jet_sel_pt:
                Jet_btag_bits[iJet] = 1
            ## b-tagged jet closest to Higgs candidate AK8 jet
            if iCandH >= 0 and vJet.DeltaR(vCandH) < self.bJet_candH_dR:
                if iJetCandH >= 0 and vJetCandH.DeltaR(vCandH) <= vJet.DeltaR(vCandH): continue
                Jet_candH_bits  = [0] * nJets
                Jet_candH_bits[iJet] = 1
                iJetCandH = iJet
                vJetCandH = vJet
        ## End loop: for jet in Jets

        ## Get b-tagged jets which are not associated with Higgs decay
        Jet_sel_nonH_bits  = [Jet_sel_bits[ii] *(1 - Jet_candH_bits[ii]) for ii in range(len(Jet_sel_bits))]
        Jet_btag_nonH_bits = [Jet_btag_bits[ii]*(1 - Jet_candH_bits[ii]) for ii in range(len(Jet_btag_bits))]
        ## Get selected and b-tagged jets in HEM region
        Jet_sel_HEM_bits  = [Jet_sel_bits [ii]*Jet_HEM_bits[ii] for ii in range(len(Jet_sel_bits))]
        Jet_btag_HEM_bits = [Jet_btag_bits[ii]*Jet_HEM_bits[ii] for ii in range(len(Jet_btag_bits))]
        ## Pick 2 highest-pT jets for di-jet pair, excluding jetCandH
        iSelJetsNonH = [ij for ij in iSelJets if ij != iJetCandH]
        vSelJetsNonH = [vj for vj in vSelJets if (iJetCandH < 0 or vj.DeltaR(vJetCandH) > 0.05)]
        if len(iSelJetsNonH) != len(vSelJetsNonH):
            print('\n\nWIERD ERROR! %d vs. %s sel non-H jets!' % (len(iSelJetsNonH), len(vSelJetsNonH)))
        ## Get indices for b-tagged jets
        iBtagJets     = [ij for ij in iSelJets  if Jet_btag_bits[ij] == 1]
        iBtagJetsNonH = [ij for ij in iBtagJets if ij != iJetCandH]

        self.out.fillBranch("Jet_HEM",             Jet_HEM_bits)
        self.out.fillBranch("Jet_HEM_EM",          Jet_HEM_EM_bits)
        self.out.fillBranch("Jet_Haa4b_presel",    Jet_presel_bits)
        self.out.fillBranch("Jet_Haa4b_sel",       Jet_sel_bits)
        self.out.fillBranch("Jet_Haa4b_btag",      Jet_btag_bits)
        self.out.fillBranch("Jet_Haa4b_candH",     Jet_candH_bits)
        self.out.fillBranch("Jet_Haa4b_ovlp_iFat", Jet_ovlp_iFat)
        self.out.fillBranch("Jet_Haa4b_ovlp_iMu",  Jet_ovlp_iMu)
        self.out.fillBranch("Jet_Haa4b_ovlp_iEle", Jet_ovlp_iEle)
        if self.isMC:
            self.out.fillBranch("Jet_genPartIdx",    Jet_genPartIdx)
            self.out.fillBranch("Jet_genPart_pdgId", Jet_GP_pdgId)
            self.out.fillBranch("Jet_genPart_dR",    Jet_GP_dR)
        self.out.fillBranch("Haa4b_nJetSel",       sum(Jet_sel_bits))
        self.out.fillBranch("Haa4b_nJetSelNonH",   sum(Jet_sel_nonH_bits))
        self.out.fillBranch("Haa4b_nJetBtag",      sum(Jet_btag_bits))
        self.out.fillBranch("Haa4b_nJetBtagNonH",  sum(Jet_btag_nonH_bits))
        self.out.fillBranch("Haa4b_nFatHj",        sum(Jet_candH_bits))


        ############################
        ## Compute useful quantities
        vMET = ROOT.TLorentzVector()
        vMET.SetPtEtaPhiM(event.MET_pt, 0, event.MET_phi, 0)

        is2j = (len(vSelJetsNonH) >= 2)
        dijet_mass = ((vSelJetsNonH[0]+vSelJetsNonH[1]).M()              if is2j else -99)
        dijet_pt   = ((vSelJetsNonH[0]+vSelJetsNonH[1]).Pt()             if is2j else -99)
        dijet_dEta = (abs(vSelJetsNonH[0].Eta() - vSelJetsNonH[1].Eta()) if is2j else -99)
        dijet_dPhi = (abs(vSelJetsNonH[0].DeltaPhi(vSelJetsNonH[1]))     if is2j else -99)
        ## Pick highest-pT muons then highest-pT electrons for di-lepton pair
        hasL = (nSelLepsNonOvlpH >= 1)
        is2L = (nSelLepsNonOvlpH >= 2)
        if hasL:
            vLep1 = vSelLepsNonOvlpH[0]
            idLep1 = idSelLepsNonOvlpH[0]
            vLep1T = ROOT.TLorentzVector()
            vLep1T.SetPtEtaPhiM(vLep1.Pt(), 0, vLep1.Phi(), 0)
        if is2L:
            vLep2 = vSelLepsNonOvlpH[1]
            idLep2 = idSelLepsNonOvlpH[1]
        dilep_mass   = ((vLep1+vLep2).M()           if is2L else -99)
        dilep_pt     = ((vLep1+vLep2).Pt()          if is2L else -99)
        dilep_dR     = (vLep1.DeltaR(vLep2)         if is2L else -99)
        dilep_dPhi   = (abs(vLep1.DeltaPhi(vLep2))  if is2L else -99)
        dilep_charge = (2*(1-(idLep1>0)-(idLep2>0)) if is2L else -99)
        dilep_SFOS   = ((idLep1 + idLep2 == 0)      if is2L else 0)
        ## Pick highest-pT lepton and MET for lep+MET calculations
        lepMET_pt   = ((vLep1T+vMET).Pt()         if hasL else -99)
        lepMET_MT   = ((vLep1T+vMET).M()          if hasL else -99)
        lepMET_dPhi = (abs(vLep1T.DeltaPhi(vMET)) if hasL else -99)
        ## Other dPhi quantities
        FatH_MET_dPhi    = (abs(vCandH.DeltaPhi(vMET))        if iCandH >= 0 else -99)
        FatH_lepMET_dPhi = (abs(vCandH.DeltaPhi(vMET+vLep1))  if iCandH >= 0 and hasL else -99)
        FatH_dilep_dPhi  = (abs(vCandH.DeltaPhi(vLep1+vLep2)) if iCandH >= 0 and is2L else -99)

        self.out.fillBranch("Haa4b_iFatH",         iCandH)
        self.out.fillBranch("Haa4b_iFatHj_jet",    iJetCandH)
        self.out.fillBranch("Haa4b_iFatX1",        iCandsX[0] if len(iCandsX) > 0 else -99)
        self.out.fillBranch("Haa4b_iFatX2",        iCandsX[1] if len(iCandsX) > 1 else -99)
        self.out.fillBranch("Haa4b_iLep1",         iSelLepsNonOvlpH[0] if len(iSelLepsNonOvlpH) > 0 else -99)
        self.out.fillBranch("Haa4b_iLep2",         iSelLepsNonOvlpH[1] if len(iSelLepsNonOvlpH) > 1 else -99)
        self.out.fillBranch("Haa4b_Lep1_pdgId",    idSelLepsNonOvlpH[0] if len(idSelLepsNonOvlpH) > 0 else -99)
        self.out.fillBranch("Haa4b_Lep2_pdgId",    idSelLepsNonOvlpH[1] if len(idSelLepsNonOvlpH) > 1 else -99)
        self.out.fillBranch("Haa4b_iJetSel1",      iSelJets[0] if len(iSelJets) > 0 else -99)
        self.out.fillBranch("Haa4b_iJetSel2",      iSelJets[1] if len(iSelJets) > 1 else -99)
        self.out.fillBranch("Haa4b_iJetSelNonH1",  iSelJetsNonH[0] if len(iSelJetsNonH) > 0 else -99)
        self.out.fillBranch("Haa4b_iJetSelNonH2",  iSelJetsNonH[1] if len(iSelJetsNonH) > 1 else -99)
        self.out.fillBranch("Haa4b_iJetBtag1",     iBtagJets[0] if len(iBtagJets) > 0 else -99)
        self.out.fillBranch("Haa4b_iJetBtag2",     iBtagJets[1] if len(iBtagJets) > 1 else -99)
        self.out.fillBranch("Haa4b_iJetBtagNonH1", iBtagJetsNonH[0] if len(iBtagJetsNonH) > 0 else -99)
        self.out.fillBranch("Haa4b_iJetBtagNonH2", iBtagJetsNonH[1] if len(iBtagJetsNonH) > 1 else -99)

        self.out.fillBranch("Haa4b_FatH_mass",  vCandH.M()             if iCandH >= 0 else -99)
        self.out.fillBranch("Haa4b_FatHj_mass", (vCandH+vJetCandH).M() if iJetCandH >= 0 else -99)
        self.out.fillBranch("Haa4b_FatX1_mass", vCandsX[0].M()         if len(iCandsX) > 0 else -99)
        self.out.fillBranch("Haa4b_FatX2_mass", vCandsX[1].M()         if len(iCandsX) > 1 else -99)

        self.out.fillBranch("Haa4b_FatH_msoftdrop",  vCandHsd.M()                  if iCandH >= 0 else -99)
        self.out.fillBranch("Haa4b_FatHj_msoftdrop", (vCandHsd+vJetCandH).M()      if iJetCandH >= 0 else -99)
        self.out.fillBranch("Haa4b_FatX1_msoftdrop", FatJets[iCandsX[0]].msoftdrop if len(iCandsX) > 0 else -99)
        self.out.fillBranch("Haa4b_FatX2_msoftdrop", FatJets[iCandsX[1]].msoftdrop if len(iCandsX) > 1 else -99)

        self.out.fillBranch("Haa4b_FatH_pt",  vCandH.Pt()             if iCandH >= 0 else -99)
        self.out.fillBranch("Haa4b_FatHj_pt", (vCandH+vJetCandH).Pt() if iJetCandH >= 0 else -99)
        self.out.fillBranch("Haa4b_FatX1_pt", vCandsX[0].Pt()         if len(iCandsX) > 0 else -99)
        self.out.fillBranch("Haa4b_FatX2_pt", vCandsX[1].Pt()         if len(iCandsX) > 1 else -99)

        if self.isMC:
            self.out.fillBranch("Haa4b_FatH_nBQuarks",   sum([(vGP[i].DeltaR(vCandH) < 0.8 and abs(idGP[i]) == 5) for i in GPqs]))
            self.out.fillBranch("Haa4b_FatH_nCQuarks",   sum([(vGP[i].DeltaR(vCandH) < 0.8 and abs(idGP[i]) == 4) for i in GPqs]))
            self.out.fillBranch("Haa4b_FatH_nLFQuarks",  sum([(vGP[i].DeltaR(vCandH) < 0.8 and abs(idGP[i]) <= 3) for i in GPqs]))
            self.out.fillBranch("Haa4b_FatHj_jet_nBQuarks",  sum([((vGP[i].DeltaR(vJetCandH) < 0.4 and abs(idGP[i]) == 5)
                                                               if iJetCandH >= 0 else 0) for i in GPqs]))
            self.out.fillBranch("Haa4b_FatHj_jet_nCQuarks",  sum([((vGP[i].DeltaR(vJetCandH) < 0.4 and abs(idGP[i]) == 4)
                                                               if iJetCandH >= 0 else 0) for i in GPqs]))
            self.out.fillBranch("Haa4b_FatHj_jet_nLFQuarks", sum([((vGP[i].DeltaR(vJetCandH) < 0.4 and abs(idGP[i]) <= 3)
                                                               if iJetCandH >= 0 else 0) for i in GPqs]))
            self.out.fillBranch("Haa4b_FatX1_nBQuarks",  sum([((vGP[i].DeltaR(vCandsX[0]) < 0.8 and abs(idGP[i]) == 5)
                                                           if len(iCandsX) > 0 else 0) for i in GPqs]))
            self.out.fillBranch("Haa4b_FatX1_nCQuarks",  sum([((vGP[i].DeltaR(vCandsX[0]) < 0.8 and abs(idGP[i]) == 4)
                                                           if len(iCandsX) > 0 else 0) for i in GPqs]))
            self.out.fillBranch("Haa4b_FatX1_nLFQuarks", sum([((vGP[i].DeltaR(vCandsX[0]) < 0.8 and abs(idGP[i]) <= 3)
                                                           if len(iCandsX) > 0 else 0) for i in GPqs]))
            self.out.fillBranch("Haa4b_FatX2_nBQuarks",  sum([((vGP[i].DeltaR(vCandsX[1]) < 0.8 and abs(idGP[i]) == 5)
                                                           if len(iCandsX) > 1 else 0) for i in GPqs]))
            self.out.fillBranch("Haa4b_FatX2_nCQuarks",  sum([((vGP[i].DeltaR(vCandsX[1]) < 0.8 and abs(idGP[i]) == 4)
                                                           if len(iCandsX) > 1 else 0) for i in GPqs]))
            self.out.fillBranch("Haa4b_FatX2_nLFQuarks", sum([((vGP[i].DeltaR(vCandsX[1]) < 0.8 and abs(idGP[i]) <= 3)
                                                           if len(iCandsX) > 1 else 0) for i in GPqs]))

        self.out.fillBranch("Haa4b_FatHj_dR", vCandH.DeltaR(vJetCandH) if iJetCandH >= 0 else -99)
        self.out.fillBranch("Haa4b_FatHj_jet_pt", vJetCandH.Pt()       if iJetCandH >= 0 else -99)
        
        self.out.fillBranch("Haa4b_FatH_MET_dPhi",    FatH_MET_dPhi)
        self.out.fillBranch("Haa4b_FatH_lepMET_dPhi", FatH_lepMET_dPhi)
        self.out.fillBranch("Haa4b_FatH_dilep_dPhi",  FatH_dilep_dPhi)

        self.out.fillBranch("Haa4b_FatH_tagHaa34b_v1",  FatJets[iCandH].PNet_X4b_v1_Haa34b_vs_QCD  if iCandH >= 0 else -99)
        self.out.fillBranch("Haa4b_FatH_tagHaa34b_v2a", FatJets[iCandH].PNet_X4b_v2a_Haa34b_score if iCandH >= 0 else -99)
        self.out.fillBranch("Haa4b_FatH_tagHaa34b_v2b", FatJets[iCandH].PNet_X4b_v2b_Haa34b_score if iCandH >= 0 else -99)
        self.out.fillBranch("Haa4b_FatH_tagHaa4b_v1",  FatJets[iCandH].PNet_X4b_v1_Haa4b_vs_QCD  if iCandH >= 0 else -99)
        self.out.fillBranch("Haa4b_FatH_tagHaa4b_v2a", FatJets[iCandH].PNet_X4b_v2a_Haa4b_score if iCandH >= 0 else -99)
        self.out.fillBranch("Haa4b_FatH_tagHaa4b_v2b", FatJets[iCandH].PNet_X4b_v2b_Haa4b_score if iCandH >= 0 else -99)

        self.out.fillBranch("Haa4b_FatX_tagWZ_max",  max_WZ_tag)
        self.out.fillBranch("Haa4b_FatX_tagTop_max", max_top_tag)
        self.out.fillBranch("Haa4b_FatX1_tagWZ",  FatJets[iCandsX[0]].particleNet_WZvsQCD if len(iCandsX) > 0 else -99)
        self.out.fillBranch("Haa4b_FatX2_tagWZ",  FatJets[iCandsX[1]].particleNet_WZvsQCD if len(iCandsX) > 1 else -99)
        self.out.fillBranch("Haa4b_FatX1_tagTop", FatJets[iCandsX[0]].particleNet_TvsQCD  if len(iCandsX) > 0 else -99)
        self.out.fillBranch("Haa4b_FatX2_tagTop", FatJets[iCandsX[1]].particleNet_TvsQCD  if len(iCandsX) > 1 else -99)

        self.out.fillBranch("Haa4b_dijet_mass",   dijet_mass)
        self.out.fillBranch("Haa4b_dijet_pt",     dijet_pt)
        self.out.fillBranch("Haa4b_dijet_dEta",   dijet_dEta)
        self.out.fillBranch("Haa4b_dijet_dPhi",   dijet_dPhi)
        self.out.fillBranch("Haa4b_dilep_mass",   dilep_mass)
        self.out.fillBranch("Haa4b_dilep_pt",     dilep_pt)
        self.out.fillBranch("Haa4b_dilep_dR",     dilep_dR)
        self.out.fillBranch("Haa4b_dilep_dPhi",   dilep_dPhi)
        self.out.fillBranch("Haa4b_dilep_charge", dilep_charge)
        self.out.fillBranch("Haa4b_dilep_SFOS",   dilep_SFOS)
        self.out.fillBranch("Haa4b_lepMET_pt",    lepMET_pt)
        self.out.fillBranch("Haa4b_lepMET_MT",    lepMET_MT)
        self.out.fillBranch("Haa4b_lepMET_dPhi",  lepMET_dPhi)


        ####################
        ## Define categories
        cat_bit = {}
        HEM_bit = {}
        HEM_candH = ((iCandH >= 0    and Fat_HEM_bits[iCandH] == 1) or \
                     (iJetCandH >= 0 and Jet_HEM_bits[iJetCandH] == 1))
        for cat in CATS:
            cat_bit[cat] = 0
            HEM_bit[cat] = HEM_candH
        ## 0-lepton categories: gg0l, VBFjj, Vjj, ttHad, Zvv
        if iCandH >= 0 and nTrigLeps == 0:
            ## For now, only make "exclusive" categories among high-MET and low-MET >=1b, 0b+dijet, and 0b
            if vMET.Pt() > 200:
                if FatH_MET_dPhi > 0.5*math.pi:
                    cat_bit['Zvv'] = 1
                for ii in range(len(vJets)):
                    if Jet_HEM_EM_bits[ii] == 1 and abs(vJets[ii].DeltaPhi(vMET)) < 0.3:
                        HEM_bit['Zvv'] = 1
            elif sum(Jet_btag_nonH_bits) > 0:
                cat_bit['ttHad'] = 1
                HEM_bit['ttHad'] = HEM_candH or (sum(Jet_btag_HEM_bits) > 0)
            elif dijet_mass > 450 and dijet_dEta > 2.2:
                ## https://indico.cern.ch/event/1450818/#34-mohamed-darwish
                cat_bit['VBFjj'] = 1
                HEM_bit['VBFjj'] = HEM_candH or (sum(Jet_sel_HEM_bits) > 0)
            else:
                cat_bit['gg0l'] = 1
            ## Hadronic categories based on 2nd AK8 jet are provisional, non-exclusive
            if len(iCandsX) > 0 and vMET.Pt() <= 200:
                cat_bit['Vjj']   = (max_WZ_tag  > self.Fat_candX_Vqq)
                cat_bit['ttHad'] = (max_top_tag > self.Fat_candX_Tbqq) or cat_bit['ttHad']
                HEM_bit['Vjj']   = HEM_candH or (sum([Fat_HEM_bits[iX] for iX in iCandsX]) > 0)
                HEM_bit['ttHad'] = HEM_candH or (sum([Fat_HEM_bits[iX] for iX in iCandsX]) > 0) or HEM_bit['ttHad']
        ## >=1-lepton categories: Zll, Wlv, ttlv, ttll, 2lSS, 3l
        if iCandH >= 0 and nTrigLeps >= 1 and nTrigLepsOvlpH == 0:
            for lepCat in ['3l','ttll','Zll','2lSS','ttlv','Wlv']:
                HEM_bit[lepCat] = HEM_candH or (sum(Ele_sel_HEM_bits) > 0)
            if nSelLepsNonOvlpH >= 3:
                cat_bit['3l'] = 1
            elif nSelLepsNonOvlpH == 2:
                if dilep_charge == 0:
                    if dilep_SFOS == 1 and abs(dilep_mass - 90.0) < 10.0:
                        cat_bit['Zll'] = 1
                    elif sum(Jet_btag_nonH_bits) > 0 and dilep_mass > 12.0:
                        cat_bit['ttll'] = 1
                else:
                    cat_bit['2lSS'] = 1
            elif nSelLepsNonOvlpH == 1:
                if sum(Jet_btag_nonH_bits) > 0:
                    cat_bit['ttlv'] = 1
                elif FatH_lepMET_dPhi > 0.75*math.pi:
                    cat_bit['Wlv'] = 1
        ## Figure how many good Higgs candidate events we're losing "through the cracks"
        nCats = sum([cat_bit[cat] for cat in CATS])
        cat_bit['other'] = iCandH >= 1 and nCats == 0
        ## Fill category bits
        for cat in CATS:
            self.out.fillBranch("Haa4b_cat_%s" % cat, cat_bit[cat])
            self.out.fillBranch("Haa4b_HEM_%s" % cat, HEM_bit[cat])
            
	cat_idx = (2*cat_bit['VBFjj'] + 5*cat_bit['Zvv'] + 6*cat_bit['Zll'] +
                   7*cat_bit['Wlv'] + 8*cat_bit['ttlv'] + 9*cat_bit['ttll'])
        if cat_idx == 0:
            if cat_bit['ttHad'] == 1 and (cat_bit['Vjj'] == 0 or max_top_tag > self.Fat_candX_Tbqq):
                cat_idx = 4
            elif cat_bit['Vjj'] == 1:
                cat_idx = 3
            elif cat_bit['gg0l'] == 1:
                cat_idx = 1
        self.out.fillBranch("Haa4b_cat_idx", cat_idx)
        self.out.fillBranch("Haa4b_nCats",   nCats)

        self.out.fillBranch("year",   year)
        self.out.fillBranch("isData", (not self.isMC))
        self.out.fillBranch("isMC",   self.isMC)
        self.out.fillBranch("Haa4b_isHad", (iCandH >= 0 and nTrigLeps == 0))
        self.out.fillBranch("Haa4b_isLep", (iCandH >= 0 and nTrigLeps >= 1 and nTrigLepsOvlpH == 0))
        self.out.fillBranch("Haa4b_isMu",  (iCandH >= 0 and len(iTrigMus) >= 1 and nTrigLepsOvlpH == 0))
        self.out.fillBranch("Haa4b_isEle", (iCandH >= 0 and len(iTrigEles) >= 1 and len(iTrigMus) == 0 and nTrigLepsOvlpH == 0))

        ## Logic for triggers only available part of the year (2018)
        trig = {}
        try: trig['HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71'] = \
           event.HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71
        except: trig['HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71'] = False
        try: trig['HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2'] = \
           event.HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2
        except: trig['HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2'] = False
        try: trig['HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1'] = \
           event.HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1
        except: trig['HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1'] = False
        try: trig['HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1'] = \
           event.HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1
        except: trig['HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1'] = False
        try: trig['HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4'] = \
           event.HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4
        except: trig['HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4'] = False
        ## End logic for triggers only available part of the year (2018)


        self.out.fillBranch("Haa4b_trigFat",  ( ((event.HLT_PFJet500 or event.HLT_AK8PFJet500 or event.HLT_AK8PFJet400_TrimMass30 or
                                                  trig['HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4']) and (event.L1_SingleJet180)) or
                                                ((event.HLT_AK8PFHT800_TrimMass50 or event.HLT_PFHT1050) and (event.L1_SingleJet180 or event.L1_HTT360er)) ) )
        self.out.fillBranch("Haa4b_trigBtag", ( (trig['HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71'] and
                                                 (event.L1_DoubleJet112er2p3_dEta_Max1p6 or event.L1_DoubleJet150er2p5)) or
                                                (event.HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5 and
                                                 (event.L1_HTT320er or event.L1_HTT360er or event.L1_HTT400er or event.L1_ETT2000 or
                                                  event.L1_HTT320er_QuadJet_70_55_40_40_er2p4 or
                                                  event.L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3)) ) )
        self.out.fillBranch("Haa4b_trigVBF", ( (trig['HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2'] or
                                                trig['HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1']) and
                                               (event.L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5 or
                                                event.L1_HTT320er or event.L1_SingleJet180) ) )
        self.out.fillBranch("Haa4b_trigMET", ( ((event.HLT_PFMET120_PFMHT120_IDTight_PFHT60 or event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60) and
                                                (event.L1_ETMHF90_HTT60er or event.L1_ETMHF100_HTT60er or event.L1_ETMHF110_HTT60er)) or
                                               ((trig['HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1'] or
                                                 event.HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned or event.HLT_PFMETTypeOne140_PFMHT140_IDTight) and
                                                (event.L1_ETMHF100 or event.L1_ETMHF110 or event.L1_ETMHF120 or event.L1_ETMHF130)) ) )
        self.out.fillBranch("Haa4b_trigMu",  ( (event.HLT_IsoMu24 or event.HLT_Mu50) and (event.L1_SingleMu22 or event.L1_SingleMu25) ) )
        self.out.fillBranch("Haa4b_trigEle", ( (event.HLT_Ele32_WPTight_Gsf or event.HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 or
                                                event.HLT_Ele115_CaloIdVT_GsfTrkIdT or event.HLT_Ele35_WPTight_Gsf_L1EGMT) ) )
        self.out.fillBranch("Haa4b_passFilters", (event.Flag_goodVertices and event.Flag_globalSuperTightHalo2016Filter and event.Flag_HBHENoiseFilter and
                                                  event.Flag_HBHENoiseIsoFilter and event.Flag_eeBadScFilter and event.Flag_BadPFMuonFilter and
                                                  event.Flag_BadPFMuonDzFilter and event.Flag_ecalBadCalibFilter and event.Flag_EcalDeadCellTriggerPrimitiveFilter) )


        ############
        ## All done!
        return True

Haa4bObjectSelectionBranches2018MC   = lambda: Haa4bObjectSelectionProducer(True,  '2018')
Haa4bObjectSelectionBranches2018Data = lambda: Haa4bObjectSelectionProducer(False, '2018')
Haa4bObjectSelectionBranches2017MC   = lambda: Haa4bObjectSelectionProducer(True,  '2017')
Haa4bObjectSelectionBranches2017Data = lambda: Haa4bObjectSelectionProducer(False, '2017')
Haa4bObjectSelectionBranches2016MC   = lambda: Haa4bObjectSelectionProducer(True,  '2016')
Haa4bObjectSelectionBranches2016Data = lambda: Haa4bObjectSelectionProducer(False, '2016')
