#!/usr/bin/env python
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetUncertainties import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.PrefireCorr import *
from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducerUL import *
from PhysicsTools.NanoAODTools.postprocessing.modules.haa4b.objectSelection import Haa4bObjectSelectionProducer
from PhysicsTools.NanoAODTools.postprocessing.modules.haa4b.triggerBits import Haa4bTriggerBitsProducer
from PhysicsTools.NanoAODTools.postprocessing.modules.haa4b.genParticles import Haa4bGenParticlesBranches
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from importlib import import_module
import os
import sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import re


runLocally = False
isMC = True
DataYear  = 'UL2018'
RunPeriod = 'D'
  
RunEraForMC = {
    'UL2016_preVFP': 'E',
    'UL2016': 'G',
    'UL2017': 'F',
    'UL2018': 'D',      
}

if not runLocally:
    # This takes care of converting the input files from CRAB
    from PhysicsTools.NanoAODTools.postprocessing.framework.crabhelper import inputFiles, runsAndLumis

    crabFiles = inputFiles()
    for iFile in range(0, len(crabFiles)):

        # '/store/mc/'
        dataType = re.search('/store/(?P<DataType>\w+)/', crabFiles[iFile]).group('DataType')
        isMC = True if dataType.lower() == 'mc' else False

        if not isMC:
            # '/store/data/Run2018A/'
            r_ = re.search('Run(?P<DataYear>\d{4})(?P<RunPeriod>[a-zA-Z])', crabFiles[iFile])
            DataYear  = r_.group('DataYear')
            RunPeriod = r_.group('RunPeriod')
        else:
            # '/store/mc/RunIISummer20UL18MiniAODv2/'
            r_ = re.search('RunIISummer20UL(?P<DataYear>\d{2})MiniAOD', crabFiles[iFile])
            DataYear  = '20%s' % (r_.group('DataYear'))
            if 'APV' in crabFiles[iFile]:
                DataYear += '_preVFP'

        if 'UL' in crabFiles[iFile]:
            DataYear = 'UL%s' % (DataYear)

        if isMC:
            RunPeriod = RunEraForMC[DataYear]

        print('Hto4b_postproc.py:: ', crabFiles[iFile], DataYear, RunPeriod)
        break

    print('Hto4b_postproc.py:: DataYear: ', DataYear, ', RunPeriod: ', RunPeriod)




jmeCorrectionsAK8 = createJMECorrector(
    isMC          = isMC, 
    dataYear      = DataYear, 
    runPeriod     = RunPeriod, 
    jesUncert     = "Total", # Options: Total, All
    jetType       = "AK8PFPuppi", # Options: AK8PFPuppi, AK4PFchs
    noGroom       = False,
    applyWJMS     = False, 
    applyMsdJMS   = False,
    applyMsdJMR   = False,
    applySmearing = True,
    isFastSim     = False,
    applyHEMfix   = ('18' in DataYear), # Additional "HEMIssue" uncertainty
    splitJER      = False,
    saveMETUncs   = ['T1', 'T1Smear'])

jmeCorrectionsAK4 = createJMECorrector(
    isMC          = isMC, 
    dataYear      = DataYear, 
    runPeriod     = RunPeriod, 
    jesUncert     = "Total", 
    jetType       = "AK4PFchs", # Options: AK8PFPuppi, AK4PFchs
    noGroom       = False,
    applySmearing = True,
    isFastSim     = False,
    applyHEMfix   = ('18' in DataYear), # Additional "HEMIssue" uncertainty
    splitJER      = False,
    saveMETUncs   = ['T1', 'T1Smear'])


modulesToRun = [
    jmeCorrectionsAK8(),
    jmeCorrectionsAK4()
]

if isMC:
    ## PU weights
    puWeights = None
    if   '16' in DataYear:
        puWeights = puWeight_UL2016
    elif '17' in DataYear:
        puWeights = puWeight_UL2017
    elif '18' in DataYear:
        puWeights = puWeight_UL2018
    else:
        print('Hto4b_postproc.py:: puWeights DataYear ', DataYear, ' not compatible *** ERROR ***')
        exit(0)
        

    ## Prefiring corrections
    '''
    l1PrefirCorr2017 = lambda: PrefCorr(
        jetroot="L1prefiring_jetpt_2017BtoF.root",
        jetmapname="L1prefiring_jetpt_2017BtoF",
        photonroot="L1prefiring_photonpt_2017BtoF.root",
        photonmapname="L1prefiring_photonpt_2017BtoF",
        branchnames=[
            "PrefireWeight", "PrefireWeight_Up", "PrefireWeight_Down"
        ]
    )
    '''

    # b-tag SF producer
    btagSF = lambda: btagSFProducer(
        era = DataYear, 
        algo='deepjet', 
        selectedWPs=['M'], #['M', 'shape_corr'],
        sfFileName=None, 
        verbose=0, 
        jesSystsForShape=["jes"]
    )

    modulesToRun.extend([
        puWeights(),
        btagSF(),
        Haa4bGenParticlesBranches(),
    ])

    modulesToRun.extend([
        Haa4bObjectSelectionProducer(isMC, DataYear.replace('UL', ''))
    ])

    if '18' in DataYear:
        modulesToRun.extend([
            Haa4bTriggerBitsProducer(isMC, DataYear.replace('UL', ''))
        ])

fnames = ["PNet_v1.root"] 
if runLocally:
    fnames = ["/afs/cern.ch/work/s/ssawant/private/htoaa/NanoAODProduction_wPNetHToAATo4B/CMSSW_10_6_30/src/PhysicsTools/NanoAOD/output/HtoAA_addHto4bPlus_ggH_HtoAA_MH-125_MA-32.5_Pt170_Eta2p4_Msoft10_Xbb0p6_skimFatCand_1k.root"]



p = PostProcessor(
    outputDir=".", 
    inputFiles=fnames, 
    #cut="",  
    modules=modulesToRun, 
    #provenance=True,
    #fwkJobReport=True,
    #jsonInput=runsAndLumis()
    )
p.run()

print("DONE Hto4b_postproc.py")
