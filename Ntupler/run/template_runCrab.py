import os
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = '_requestName_'
config.General.workArea = '_workArea_'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '_psetName_'
# config.JobType.sendExternalFolder = True

config.Data.inputDataset = '_inputDataset_'
config.JobType.pyCfgParams = _pyCfgParams_
config.JobType.pyCfgParams += ['inputDataset=%s' % config.Data.inputDataset]
# config.JobType.inputFiles = []
config.JobType.sendExternalFolder = True

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = _unitsPerJob_
#config.Data.splitting = 'Automatic'
#config.Data.unitsPerJob = 200
#config.Data.outLFNDirBase = '/store/user/%s/_outLFNDirBase_' % (getUsernameFromSiteDB())
config.Data.outLFNDirBase = '_outLFNDirBase_'
config.Data.publication = False
config.Data.allowNonValidInputDataset = True

# uncomment this part to use CRAB to submit to FNAL_LPC
# config.Data.ignoreLocality = True
# config.Site.whitelist = ['T3_US_FNALLPC']
# config.Site.ignoreGlobalBlacklist = True

# Run in DESY/RWTH
#config.Site.whitelist = ['T2_CH_CERN', 'T2_DE_DESY', 'T2_DE_RWTH']
config.Site.whitelist = ['T2_DE_RWTH']
#config.Site.ignoreGlobalBlacklist = True
#config.Site.blacklist = ['T3_US_UMiss']
config.User.voGroup = 'dcms'
config.User.voRole = 'NULL'

#config.Site.storageSite = '_storageSite_'
config.Site.storageSite = "T2_DE_RWTH"


# uncomment this part to use CRAB to submit to FNAL LPC
#config.Data.ignoreLocality = True
#config.Site.whitelist = ['T3_US_FNALLPC']
#config.Site.ignoreGlobalBlacklist = True

# uncomment this to run in DESY/RWTH
#config.Data.ignoreLocality = True
#config.Site.whitelist = ['T2_CH_CERN', 'T2_DE_DESY', 'T2_DE_RWTH']
#config.Site.ignoreGlobalBlacklist = True

config.Site.storageSite = '_storageSite_'
