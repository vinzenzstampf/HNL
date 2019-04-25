'''
######################################################################################################
## FROM heppy_crab.py 
######################################################################################################
# datasets to run as defined from cfg file
# number of jobs to run per dataset decided based on splitFactor and fineSplitFactor from cfg file
'''

import imp, os, json
from PhysicsTools.HeppyCore.framework.heppy_loop import _heppyGlobalOptions
_heppyGlobalOptions["isCrab"] = True
optjsonfile = open('options.json','w')
optjsonfile.write(json.dumps(_heppyGlobalOptions))
optjsonfile.close()

handle = open('../hn3l_cfg.py', 'r')
cfo = imp.load_source('../hn3l_cfg.py'.split('/')[-1].rstrip(".py"), '../hn3l_cfg.py', handle)
conf = cfo.config
handle.close()

os.system("tar czf python.tar.gz --dereference --directory $CMSSW_BASE python")
os.system("tar czf cmgdataset.tar.gz --directory $HOME .cmgdataset")
os.system("tar czf cafpython.tar.gz --directory /afs/cern.ch/cms/caf/ python")

os.environ["PROD_LABEL"]     = 'test_heppy_crab'
os.environ["CMG_VERSION"]    = 'CMSSW_9_4_6_patch1'
os.environ["USEAAA"]         = 'full'
os.environ["STAGEOUTREMDIR"] = 'test_heppyTrees'
os.environ["CFG_FILE"]       = '../hn3l_cfg.py'
os.environ["OUTSITE"]        = 'T3_CH_PSI'
os.environ["ONLYUNPACKED"]   = 'False'
os.environ["MAXNUMEVENTS"]   = str(10000)

from PhysicsTools.HeppyCore.framework.heppy_loop import split
#for comp in conf.components:
#    if getattr(comp,"useAAA",False):
#        raise RuntimeError, 'Components should have useAAA disabled in the cfg when running on crab - tune the behaviour of AAA in the crab submission instead!'
#    os.environ["DATASET"] = str(comp.name)
#    os.environ["NJOBS"] = str(len(split([comp])))
#    os.system("crab submit %s -c heppy_crab_config_env.py"%("--dryrun" if options.dryrun else ""))
comp = conf.components[0]
if getattr(comp,"useAAA",False):
    raise RuntimeError, 'Components should have useAAA disabled in the cfg when running on crab - tune the behaviour of AAA in the crab submission instead!'
os.environ["DATASET"] = str(comp.name)
os.environ["NJOBS"] = str(len(split([comp])))
#os.environ["ONLYUNPACKED"] = str(options.only_unpacked)

'''
######################################################################################################
## FROM heppy_crab_config.py 
######################################################################################################
here we set all crab options that are not fixed
values we'll be taken from environment variables set in launchall.py
fixed options will be taken from heppy_crab_config.py
'''

import os
from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName   = 'heppy_crab_fake_pset.py'
config.JobType.scriptExe  = 'heppy_crab_script.sh'
config.JobType.disableAutomaticOutputCollection = True
# config.JobType.sendPythonFolder = True  #doesn't work, not supported yet? do it by hand

config.JobType.inputFiles = ['FrameworkJobReport.xml','heppy_crab_script.py','cmgdataset.tar.gz', 'python.tar.gz', 'cafpython.tar.gz','options.json']
config.JobType.outputFiles = []

config.section_("Data")
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventBased'
config.Data.outLFNDirBase = '/store/user/' + os.environ["USER"]
config.Data.publication = False

config.section_("Site")

'''
######################################################################################################
## FROM heppy_crab_config_env.py 
######################################################################################################
import imp, os
file = open( "heppy_crab_config.py", 'r' )
cfg = imp.load_source( 'cfg', "heppy_crab_config.py", file)
config = cfg.config
'''

print "Will send dataset", os.environ["DATASET"], "with", os.environ["NJOBS"], "jobs"

config.General.requestName = os.environ["DATASET"] + "_" + os.environ["CMG_VERSION"] # task name
config.General.workArea = 'crab_' + os.environ["DATASET"] + "_" + os.environ["PROD_LABEL"] # crab dir name

# this will divide task in *exactly* NJOBS jobs (for this we need JobType.pluginName = 'PrivateMC' and Data.splitting = 'EventBased')
config.Data.unitsPerJob = 10
config.Data.totalUnits = config.Data.unitsPerJob * int(os.environ["NJOBS"])

config.JobType.inputFiles.append(os.environ["CFG_FILE"])
# arguments to pass to scriptExe. They have to be like "arg=value". 
config.JobType.scriptArgs = ["dataset="+os.environ["DATASET"], "total="+os.environ["NJOBS"], "useAAA="+os.environ["USEAAA"], "cfgfile="+os.environ["CFG_FILE"].split('/')[-1]]

#final output: /store/user/$USER/output_dir/cmg_version/production_label/dataset/$date_$time/0000/foo.bar
config.Data.outLFNDirBase += '/' + os.environ["STAGEOUTREMDIR"] + '/' + os.environ["CMG_VERSION"]
config.Data.outputPrimaryDataset =  os.environ["PROD_LABEL"]
config.Data.outputDatasetTag = os.environ["DATASET"]
config.Data.ignoreLocality = (os.environ["USEAAA"]!="local") # "full" or "eos"
if os.environ["ONLYUNPACKED"]!="True": config.JobType.outputFiles.append("heppyOutput.tgz")
if (os.environ["USEAAA"]!="local"): config.Site.whitelist = ["T2_CH_CSCS", "T2_IT_Legnaro", "T2_UK_London_IC", "T2_UK_SGrid_Bristol", "T2_DE_DESY", "T2_ES_CIEMAT", "T2_IT_Rome", "T2_AT_Vienna","T2_DE_RWTH","T2_FR_GRIF_IRFU", "T2_HU_Budapest", "T2_FR_IPHC", "T2_BE_IIHE", "T2_IT_Pisa", "T2_ES_IFCA", "T2_UK_London_Brunel", "T2_US_Purdue", "T2_UA_KIPT", "T2_US_MIT", "T2_US_Wisconsin", "T2_US_UCSD", "T2_US_Vanderbilt", "T2_US_Caltech"]

config.Site.storageSite = os.environ["OUTSITE"]
Dryrun = True
os.system("crab submit %s -c heppy_crab_config_env.py"%("--dryrun" if Dryrun else ""))

'''
######################################################################################################
## FROM heppy_crab.py 
######################################################################################################
'''
os.system("rm options.json")
os.system("rm python.tar.gz")
os.system("rm cmgdataset.tar.gz")
os.system("rm cafpython.tar.gz")
