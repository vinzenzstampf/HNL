import ROOT
import numpy as np
from DataFormats.FWLite import Events, Handle
from CMGTools.RootTools.samples.ComponentCreator import ComponentCreator
from PhysicsTools.Heppy.physicsutils.PileUpSummaryInfo import PileUpSummaryInfo
import multiprocessing as mtp
from pdb import set_trace

def FILL(events):
    for j, event in enumerate(events):
        if j%1000==0:
            print '\t\tprocessing the %d-th event of the %d-th batch' %(j, i+1)
        event.getByLabel(label, handle)
        puinfos = map(PileUpSummaryInfo, handle.product())
        for pu in puinfos:
            if pu.getBunchCrossing()==0:
                h_ti.Fill(pu.nTrueInteractions())
                break

pool = mtp.Pool(8)

ROOT.gROOT.SetBatch()        # don't pop up canvases

creator = ComponentCreator()

TTJets_amcat = creator.makeMCComponent(
    name    = 'TTJets_amcat', 
    dataset = '/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM',
    user    = 'CMS', 
    pattern = '.*root', 
    xSec    = 1, 
    useAAA  = True
)

outfile = ROOT.TFile.Open('pileup_TTJets_amcat.root', 'recreate')

# fill with true interactions
h_ti = ROOT.TH1F('pileup', 'pileup', 200, 0, 200)

handle  = Handle ('std::vector<PileupSummaryInfo>')
label = ("slimmedAddPileupInfo")

totevents = 0
nfiles = len(TTJets_amcat.files)
batch = 20
maxend = (nfiles - nfiles%batch) / batch + 1

for i in range(maxend):
#for i in range(3):
    
    begin  = batch * (i)
    end    = batch * (i+1)
    
    print 'running of %d-th batch of %d files out of %d total batches' %(i+1, batch, maxend)

    events = Events(TTJets_amcat.files[begin:end])

    set_trace()
    pool.Map(FILL)

    totevents += j
    print '\tprocessed %d events so far' %(events.size())

outfile.cd()
h_ti.Write()
outfile.Close()

