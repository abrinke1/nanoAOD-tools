## H --> aa --> 4b search with "boosted" decays to AK8 jets
## Add trigger and filter bits

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import ROOT
import os
import sys
import math
import random
ROOT.PyConfig.IgnoreCommandLineOptions = True

class Haa4bTriggerBitsProducer(Module):

    def __init__(self, IS_MC, YEAR):
        self.isMC     = IS_MC
        self.year_str = YEAR

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

        self.out.branch("Haa4b_trigFat",    "O")
        self.out.branch("Haa4b_trigBtag",   "O")
        self.out.branch("Haa4b_trigVBF",    "O")
        self.out.branch("Haa4b_trigMET",    "O")
        self.out.branch("Haa4b_trigMu",     "O")
        self.out.branch("Haa4b_trigEle",    "O")
        self.out.branch("Haa4b_passFilters","O")
    ## End function: beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree)

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

        ###############################
        ## 2018 triggers and filters ##
        ###############################
        if (year == 2018):

            ## Logic for triggers only available part of the year (2018)
            trig = {}
            try: trig['L1_HTT360er'] = \
               event.L1_HTT360er
            except: trig['L1_HTT360er'] = False
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


            self.out.fillBranch("Haa4b_trigFat",
                                ( ((event.HLT_PFJet500 or event.HLT_AK8PFJet500 or event.HLT_AK8PFJet400_TrimMass30 or
                                    trig['HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4']) and
                                   (event.L1_SingleJet180)) or
                                  ((event.HLT_AK8PFHT800_TrimMass50 or event.HLT_PFHT1050) and
                                   (event.L1_SingleJet180 or trig['L1_HTT360er'])) ) )
            self.out.fillBranch("Haa4b_trigBtag",
                                ( (trig['HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71'] and
                                   (event.L1_DoubleJet112er2p3_dEta_Max1p6 or event.L1_DoubleJet150er2p5)) or
                                  (event.HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5 and
                                   (event.L1_HTT320er or trig['L1_HTT360er'] or event.L1_HTT400er or event.L1_ETT2000 or
                                    event.L1_HTT320er_QuadJet_70_55_40_40_er2p4 or
                                    event.L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3)) ) )
            self.out.fillBranch("Haa4b_trigVBF",
                                ( (trig['HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2'] or
                                   trig['HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1']) and
                                  (event.L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5 or
                                   event.L1_HTT320er or event.L1_SingleJet180) ) )
            self.out.fillBranch("Haa4b_trigMET",
                                ( ((event.HLT_PFMET120_PFMHT120_IDTight_PFHT60 or
                                    event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60) and
                                   (event.L1_ETMHF90_HTT60er or event.L1_ETMHF100_HTT60er or event.L1_ETMHF110_HTT60er)) or
                                  ((trig['HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1'] or
                                    event.HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned or
                                    event.HLT_PFMETTypeOne140_PFMHT140_IDTight) and
                                   (event.L1_ETMHF100 or event.L1_ETMHF110 or event.L1_ETMHF120 or event.L1_ETMHF130)) ) )
            self.out.fillBranch("Haa4b_trigMu",
                                ( (event.HLT_IsoMu24 or event.HLT_Mu50) and
                                  (event.L1_SingleMu22 or event.L1_SingleMu25) ) )
            self.out.fillBranch("Haa4b_trigEle",
                                ( (event.HLT_Ele32_WPTight_Gsf or event.HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 or
                                   event.HLT_Ele115_CaloIdVT_GsfTrkIdT or event.HLT_Ele35_WPTight_Gsf_L1EGMT) ) )
            self.out.fillBranch("Haa4b_passFilters",
                                (event.Flag_goodVertices and event.Flag_globalSuperTightHalo2016Filter and
                                 event.Flag_HBHENoiseFilter and event.Flag_HBHENoiseIsoFilter and
                                 event.Flag_eeBadScFilter and event.Flag_BadPFMuonFilter and
                                 event.Flag_BadPFMuonDzFilter and event.Flag_ecalBadCalibFilter and
                                 event.Flag_EcalDeadCellTriggerPrimitiveFilter) )
        ## End conditional: if (year == 2018)

        ############
        ## All done!
        return True


Haa4bTriggerBitsBranches2018MC   = lambda: Haa4bTriggerBitsProducer(True,  '2018')
Haa4bTriggerBitsBranches2018Data = lambda: Haa4bTriggerBitsProducer(False, '2018')
Haa4bTriggerBitsBranches2017MC   = lambda: Haa4bTriggerBitsProducer(True,  '2017')
Haa4bTriggerBitsBranches2017Data = lambda: Haa4bTriggerBitsProducer(False, '2017')
Haa4bTriggerBitsBranches2016MC   = lambda: Haa4bTriggerBitsProducer(True,  '2016')
Haa4bTriggerBitsBranches2016Data = lambda: Haa4bTriggerBitsProducer(False, '2016')
