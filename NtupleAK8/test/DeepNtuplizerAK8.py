import FWCore.ParameterSet.Config as cms

# ---------------------------------------------------------

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')

options.outputFile = 'output.root'
options.inputFiles = [
#'file:85820ACA-657B-BC44-AC74-AACD6D54B348.root'
#'root://cmsxrootd.fnal.gov//store/relval/CMSSW_10_3_0_pre3/RelValTTbar_13/MINIAODSIM/103X_upgrade2018_realistic_v4-v1/20000/D3C5A8CB-17C6-8548-B133-177A42C8B012.root'
#'root://cmsxrootd.fnal.gov//store/mc/RunIISummer17MiniAOD/BulkGravTohhTohbbhbb_narrow_M-1000_13TeV-madgraph/MINIAODSIM/92X_upgrade2017_realistic_v10-v2/90000/04DF3196-3B99-E711-AE12-008CFAC93CE8.root'
#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAODv2/GluGluToBulkGravitonToHHTo4C_M-1000_narrow_13TeV-madgraph-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/86032B22-7943-E811-AF6F-D0946626135C.root'
#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/GluGluHToBB_M125_13TeV_powheg_pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/20000/C8932584-5006-E811-9840-141877410512.root'
#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/GluGluToBulkGravitonToHHTo4C_M-1000_narrow_13TeV-madgraph-pythia8/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/90000/0820742A-B729-E811-AC15-24BE05C6E711.root', 'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/GluGluToBulkGravitonToHHTo4C_M-1000_narrow_13TeV-madgraph-pythia8/MINIAODSIM/PU2017_94X_mc2017_realistic_v11-v1/90000/3290B112-B129-E811-9DE9-008CFAF55422.root'
#'/store/mc/RunIIAutumn18MiniAOD/GluGluToBulkGravitonToHHTo4B_M-1000_narrow_TuneCP5_PSweights_13TeV-madgraph_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/20000/4FF5E1EB-773A-4B42-9BFE-30C27EFDF3F3.root'
#'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAODv2/GluGluToBulkGravitonToHHTo4C_M-1000_narrow_13TeV-madgraph-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/86032B22-7943-E811-AF6F'
#'root://cmsxrootd.fnal.gov//store/mc/RunIISummer16MiniAODv2/BulkGravTohhTohbbhbb_narrow_M-2500_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/0A83E4E2-34B6-E611-89A0-549F35AE4FA2.root',
#'file:/eos/user/a/anovak/022C3683-D4AB-E611-AC4D-3417EBE70078.root'  #include file: for local files, for catalogues /store..
'/store/mc/RunIIFall17MiniAODv2/GluGluHToCC_M125_13TeV_powheg_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/30000/72164088-CB67-E811-9D0D-008CFA197AC4.root'
]
options.maxEvents = -1#10000

options.register('inputScript', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, "input Script")
options.register('skipEvents', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "skip N events")
options.register('job', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "job number")
options.register('nJobs', 1, VarParsing.multiplicity.singleton, VarParsing.varType.int, "total jobs")
options.register('fjKeepFlavors', [], VarParsing.multiplicity.list, VarParsing.varType.int, "Types of fatjet to keep in this sample")
# options.register('gluonReduction', 0.0, VarParsing.multiplicity.singleton, VarParsing.varType.float, "gluon reduction")
options.register('inputDataset',
                 '',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Input dataset")

options.setupTags(tag='%d', ifCond='nJobs > 1', tagArg='job')
options.parseArguments()

#options.inputDataset='/RelValQCD_FlatPt_15_3000_13/CMSSW_10_1_0_pre2-100X_mcRun2_asymptotic_v2_FastSim-v1/MINIAODSIM'

# ---------------------------------------------------------

process = cms.Process("DNNFiller")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000       
if not options.inputScript:  # this is probably for testing
	process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.options = cms.untracked.PSet(
   allowUnscheduled = cms.untracked.bool(True),  
   wantSummary=cms.untracked.bool(False)
)

print ('Using output file ' + options.outputFile)

process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string(options.outputFile))

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))

process.source = cms.Source('PoolSource',
    fileNames=cms.untracked.vstring(options.inputFiles),
    skipEvents=cms.untracked.uint32(options.skipEvents)
)

#process.source.eventsToProcess = cms.untracked.VEventRange("1:3127296-1:3127297")


if options.inputScript:
    process.load(options.inputScript)

numberOfFiles = len(process.source.fileNames)
numberOfJobs = options.nJobs
jobNumber = options.job

process.source.fileNames = process.source.fileNames[jobNumber:numberOfFiles:numberOfJobs]
if options.nJobs > 1:
    print ("running over these files:")
    print (process.source.fileNames)

# ---------------------------------------------------------

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
# process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

# ---------------------------------------------------------

bTagInfos = [
	'pfImpactParameterTagInfos',
	'pfInclusiveSecondaryVertexFinderTagInfos',
	'pfDeepCSVTagInfos',
	'pfBoostedDoubleSVAK8TagInfos'
]

bTagDiscriminators = [
	'softPFMuonBJetTags',
	'softPFElectronBJetTags',
	'pfJetBProbabilityBJetTags',
	'pfJetProbabilityBJetTags',
	'pfCombinedInclusiveSecondaryVertexV2BJetTags',
	'pfDeepCSVJetTags:probudsg',
	'pfDeepCSVJetTags:probb',
	'pfDeepCSVJetTags:probc',
	'pfDeepCSVJetTags:probbb',
	 #'pfDeepCSVJetTags:probcc',

	'pfBoostedDoubleSecondaryVertexAK8BJetTags',


      'pfDeepDoubleBvLJetTags:probHbb',
      'pfDeepDoubleCvLJetTags:probHcc',
      'pfDeepDoubleCvBJetTags:probHcc',
      'pfMassIndependentDeepDoubleBvLJetTags:probHbb',
      'pfMassIndependentDeepDoubleCvLJetTags:probHcc',
      'pfMassIndependentDeepDoubleCvBJetTags:probHcc',

    "pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:bbvsLight",
    "pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ccvsLight",
    "pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZHccvsQCD",
    "pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZHbbvsQCD"
]

jetCorrectionsAK8 = ('AK8PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(
        process,
        labelName = "DeepFlavour",
        jetSource = cms.InputTag('slimmedJetsAK8'),
        jetCorrections = jetCorrectionsAK8,
        pfCandidates = cms.InputTag('packedPFCandidates'),
        pvSource = cms.InputTag("offlineSlimmedPrimaryVertices"),
        svSource = cms.InputTag('slimmedSecondaryVertices'),
        muSource = cms.InputTag('slimmedMuons'),
        elSource = cms.InputTag('slimmedElectrons'),
        btagInfos = bTagInfos,
        btagDiscriminators = bTagDiscriminators,
        explicitJTA = False
)

if hasattr(process,'updatedPatJetsTransientCorrectedDeepFlavour'):
	process.updatedPatJetsTransientCorrectedDeepFlavour.addTagInfos = cms.bool(True) 
	process.updatedPatJetsTransientCorrectedDeepFlavour.addBTagInfo = cms.bool(True)
else:
	raise ValueError('I could not find updatedPatJetsTransientCorrectedDeepFlavour to embed the tagInfos, please check the cfg')
# ---------------------------------------------------------

'''
from RecoJets.JetProducers.ak4GenJets_cfi import ak8GenJets
process.ak8GenJetsWithNu = ak4GenJets.clone(src='packedGenParticles')
 
 ## Filter out neutrinos from packed GenParticles
process.packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedGenParticles"), cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16"))
 ## Define GenJets
process.ak8GenJetsRecluster = ak8GenJets.clone(src='packedGenParticlesForJetsNoNu')

 
process.patGenJetMatchWithNu = cms.EDProducer("GenJetMatcher",  # cut on deltaR; pick best by deltaR           
    src         = cms.InputTag("selectedUpdatedPatJetsDeepFlavour"),      # RECO jets (any View<Jet> is ok) 
    matched     = cms.InputTag("ak8GenJetsWithNu"),        # GEN jets  (must be GenJetCollection)              
    mcPdgId     = cms.vint32(),                      # n/a   
    mcStatus    = cms.vint32(),                      # n/a   
    checkCharge = cms.bool(False),                   # n/a   
    maxDeltaR   = cms.double(0.8),                   # Minimum deltaR for the match   
    #maxDPtRel   = cms.double(3.0),                  # Minimum deltaPt/Pt for the match (not used in GenJetMatcher)                     
    resolveAmbiguities    = cms.bool(True),          # Forbid two RECO objects to match to the same GEN object 
    resolveByMatchQuality = cms.bool(False),         # False = just match input in order; True = pick lowest deltaR pair first          
)

process.patGenJetMatchRecluster = cms.EDProducer("GenJetMatcher",  # cut on deltaR; pick best by deltaR           
    src         = cms.InputTag("selectedUpdatedPatJetsDeepFlavour"),      # RECO jets (any View<Jet> is ok) 
    matched     = cms.InputTag("ak8GenJetsRecluster"),        # GEN jets  (must be GenJetCollection)              
    mcPdgId     = cms.vint32(),                      # n/a   
    mcStatus    = cms.vint32(),                      # n/a   
    checkCharge = cms.bool(False),                   # n/a   
    maxDeltaR   = cms.double(0.8),                   # Minimum deltaR for the match   
    #maxDPtRel   = cms.double(3.0),                  # Minimum deltaPt/Pt for the match (not used in GenJetMatcher)                     
    resolveAmbiguities    = cms.bool(True),          # Forbid two RECO objects to match to the same GEN object 
    resolveByMatchQuality = cms.bool(False),         # False = just match input in order; True = pick lowest deltaR pair first          
)

process.genJetSequence = cms.Sequence(process.packedGenParticlesForJetsNoNu*process.ak4GenJetsWithNu*process.ak4GenJetsRecluster*process.patGenJetMatchWithNu*process.patGenJetMatchRecluster)

# ---------------------------------------------------------

# Very Loose IVF SV collection
from PhysicsTools.PatAlgos.tools.helpers import loadWithPrefix
loadWithPrefix(process, 'RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff', "looseIVF")
process.looseIVFinclusiveCandidateVertexFinder.primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices")
process.looseIVFinclusiveCandidateVertexFinder.tracks = cms.InputTag("packedPFCandidates")
process.looseIVFinclusiveCandidateVertexFinder.vertexMinDLen2DSig = cms.double(0.)
process.looseIVFinclusiveCandidateVertexFinder.vertexMinDLenSig = cms.double(0.)
process.looseIVFinclusiveCandidateVertexFinder.fitterSigmacut = 20

process.looseIVFcandidateVertexArbitrator.primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices")
process.looseIVFcandidateVertexArbitrator.tracks = cms.InputTag("packedPFCandidates")
process.looseIVFcandidateVertexArbitrator.secondaryVertices = cms.InputTag("looseIVFcandidateVertexMerger")
process.looseIVFcandidateVertexArbitrator.fitterSigmacut = 20
'''

# ---------------------------------------------------------

# DeepNtuplizer
process.load("DeepNTuples.NtupleAK8.DeepNtuplizerAK8_cfi")
process.deepntuplizer.jets = cms.InputTag('selectedUpdatedPatJetsDeepFlavour')
process.deepntuplizer.bDiscriminators = bTagDiscriminators 
process.deepntuplizer.bDiscriminators.append('pfCombinedMVAV2BJetTags')
process.deepntuplizer.LooseSVs = cms.InputTag("looseIVFinclusiveCandidateSecondaryVertices")

process.deepntuplizer.fjKeepFlavors = cms.untracked.vuint32(options.fjKeepFlavors)
process.deepntuplizer.isQCDSample = ('/QCD_' in options.inputDataset or '/RelValQCD_' in options.inputDataset)

# process.deepntuplizer.gluonReduction = cms.double(options.gluonReduction)
# process.p = cms.Path(process.QGTagger + process.genJetSequence * process.deepntuplizer)
# process.p = cms.Path(process.deepntuplizer)

#Trick to make it work in >=9_1_X
process.tsk = cms.Task()
for mod in process.producers_().itervalues():
    process.tsk.add(mod)
for mod in process.filters_().itervalues():
    process.tsk.add(mod)

process.p = cms.Path(
#	process.QGTagger + process.genJetSequence*  
	process.deepntuplizer,
	process.tsk
	)
