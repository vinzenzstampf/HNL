'''
This is the base analyzer going through data and trying to identify HNL->3L events.
'''

import ROOT
from itertools import product, combinations
from math import sqrt, pow
import PhysicsTools.HeppyCore.framework.config as cfg
from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from PhysicsTools.Heppy.physicsobjects.PhysicsObjects import Lepton
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaR2
from PhysicsTools.Heppy.physicsobjects.GenParticle   import GenParticle
from PhysicsTools.Heppy.physicsobjects.Photon        import Photon
from PhysicsTools.Heppy.physicsobjects.Tau           import Tau
from PhysicsTools.Heppy.physicsobjects.Muon          import Muon
from PhysicsTools.Heppy.physicsobjects.Electron      import Electron
()
        chi2 = 0.
        ndof = 0.
        bsvtx = ROOT.reco.Vertex(point, error, chi2, ndof, 2) # size? say 3? does it matter?
                                        
        event.recoSv.disp2DFromBS      = ROOT.VertexDistanceXY().distance(event.recoSv, bsvtx)
        event.recoSv.disp2DFromBS_sig  = event.recoSv.disp2DFromBS.significance()
        event.recoSv.prob              = ROOT.TMath.Prob(event.recoSv.chi2(), int(event.recoSv.ndof()))
        
        dilep_p4 = event.displaced_dilepton_reco_cand.lep1().p4() + event.displaced_dilepton_reco_cand.lep2().p4()

        perp = ROOT.math.XYZVector(dilep_p4.px(),
                                   dilep_p4.py(),
                                   0.)
        
        dxybs = ROOT.GlobalPoint(-1*((event.beamspot.x0() - event.recoSv.x()) + (event.recoSv.z() - event.beamspot.z0()) * event.beamspot.dxdz()), 
                                 -1*((event.beamspot.y0() - event.recoSv.y()) + (event.recoSv.z() - event.beamspot.z0()) * event.beamspot.dydz()),
                                  0)
        
        vperp = ROOT.math.XYZVector(dxybs.x(), dxybs.y(), 0.)
        
        cos = vperp.Dot(perp)/(vperp.R()*perp.R())
        
        event.recoSv.disp2DFromBS_cos = cos

        ########################################################################################
        # Define event.selectedLeptons, will be used by JetAnalyzer.py
        ########################################################################################

        # the selected 3 leptons must be leptons and not jets
        event.selectedLeptons = [event.the_3lep_cand.l0(),
                                 event.the_3lep_cand.l1(),
                                 event.the_3lep_cand.l2()]

        # plus any isolated electron or muon is also a good lepton rather than a jet
        event.selMuons     = [mu  for mu  in event.muons     if self.preselectPromptMuons    (mu , pt=10) and mu .relIsoR(R=0.3, dBetaFactor=0.5, allCharged=0)<0.15]
        event.selElectrons = [ele for ele in event.electrons if self.preselectPromptElectrons(ele, pt=10) and ele.relIsoR(R=0.3, dBetaFactor=0.5, allCharged=0)<0.15]
        # RM: what about taus?

        event.selectedLeptons += event.selMuons
        event.selectedLeptons += event.selElectrons

        ########################################################################################
        # Extra prompt and isolated lepton veto
        ########################################################################################        
        event.veto_mus   = [ele for ele in event.selMuons     if ele.physObj not in [event.the_3lep_cand.l0().physObj, event.the_3lep_cand.l1().physObj, event.the_3lep_cand.l2().physObj] ]
        event.veto_eles  = [mu  for mu  in event.selElectrons if mu .physObj not in [event.the_3lep_cand.l0().physObj, event.the_3lep_cand.l1().physObj, event.the_3lep_cand.l2().physObj] ]

        if len(event.veto_eles): event.veto_save_ele = sorted([ele for ele in event.veto_eles], key = lambda x : x.pt, reverse = True)[0] 
        if len(event.veto_mus ): event.veto_save_mu  = sorted([mu  for mu  in event.veto_mus ], key = lambda x : x.pt, reverse = True)[0] 

        ########################################################################################
        # charged PF isolation
        ########################################################################################        

        event.the_3lep_cand.abs_tot_iso03_rhoArea    = totIso(event, 'rhoArea', 0.3) 
        event.the_3lep_cand.abs_tot_iso04_rhoArea    = totIso(event, 'rhoArea', 0.4) 
        event.the_3lep_cand.abs_tot_iso05_rhoArea    = totIso(event, 'rhoArea', 0.5) 

        event.the_3lep_cand.rel_tot_iso03_rhoArea    = event.the_3lep_cand.abs_tot_iso03_rhoArea / event.the_3lep_cand.hnVisP4().pt()
        event.the_3lep_cand.rel_tot_iso04_rhoArea    = event.the_3lep_cand.abs_tot_iso04_rhoArea / event.the_3lep_cand.hnVisP4().pt()
        event.the_3lep_cand.rel_tot_iso05_rhoArea    = event.the_3lep_cand.abs_tot_iso05_rhoArea / event.the_3lep_cand.hnVisP4().pt()

        event.the_3lep_cand.abs_tot_iso03_deltaBeta  = totIso(event, 'dBeta', 0.3) 
        event.the_3lep_cand.abs_tot_iso04_deltaBeta  = totIso(event, 'dBeta', 0.4) 
        event.the_3lep_cand.abs_tot_iso05_deltaBeta  = totIso(event, 'dBeta', 0.5) 

        event.the_3lep_cand.rel_tot_iso03_deltaBeta  =  event.the_3lep_cand.abs_tot_iso03_deltaBeta /  event.the_3lep_cand.hnVisP4().pt()
        event.the_3lep_cand.rel_tot_iso04_deltaBeta  =  event.the_3lep_cand.abs_tot_iso04_deltaBeta /  event.the_3lep_cand.hnVisP4().pt()
        event.the_3lep_cand.rel_tot_iso05_deltaBeta  =  event.the_3lep_cand.abs_tot_iso05_deltaBeta /  event.the_3lep_cand.hnVisP4().pt()


        #####################################################################################
        # After passing all selections and we have an HNL candidate, pass a "true" boolean!
        #####################################################################################
        return True



def totIso(event, offset_mode, dRCone):
    ch_pu_iso = chargedHadronIso(event, dRCone, True)
    ch_pv_iso = chargedHadronIso(event, dRCone, False)
    neu_iso   = neutralHadronIso(event, dRCone)
    ph_iso    = photonIso(event, dRCone)
    if offset_mode == 'rhoArea': 
        eta = event.the_3lep_cand.hnVisP4().eta()
        offset = offset_rhoArea(event.rho, dRCone, eta)
    if offset_mode == 'dBeta': 
        offset = offset_dBeta(0.5, ch_pu_iso)
    tot_iso = ch_pv_iso + max(0., ph_iso + neu_iso - offset)
    # if dRCone == 0.3:
        # print '2M dr %.1f: ch_pv_iso: %.2f, neu_iso: %.2f, ph_iso: %.2f, ch_pu_iso: %.2f, l1+l2pt: %.2f, id: %i'%(dRCone, ch_pv_iso, neu_iso, ph_iso, ch_pu_iso, event.the_3lep_cand.l1().pt()+event.the_3lep_cand.l2().pt(), event.eventId)
    return tot_iso

def offset_rhoArea(rho, dRCone, eta):
    area = 0.0
    if abs(eta) < 0.8000: area = 0.0566
    if abs(eta) > 0.8000 and abs(eta) < 1.3000: area = 0.0562
    if abs(eta) > 1.3000 and abs(eta) < 2.0000: area = 0.0363
    if abs(eta) > 2.0000 and abs(eta) < 2.2000: area = 0.0119
    if abs(eta) > 2.2000 and abs(eta) < 2.4000: area = 0.0064
    if dRCone != 0.3: area *= ( (dRCone ** 2) / (0.3 **2) )
    # print 'area = {a}, offset = {o}'.format(a = area, o = area * rho) 
    return area * rho
    
def offset_dBeta(dBeta, ch_pu_iso):
    # print 'dbeta = {db}, offset = {o}'.format(db = dBeta, o = dBeta * ch_pu_iso)
    return ch_pu_iso * dBeta
        
def chargedHadronIso(event, dRCone, PU = False): 
    if PU == True:
        charged_pfs = [ipf for ipf in event.pfs if ( ipf.charge() != 0 and ipf.pt() > 0.5 )]
        charged_pfs = [ipf for ipf in charged_pfs if ipf.fromPV() <= 1]
    if PU == False:
        charged_pfs = [ipf for ipf in event.pfs if ( ipf.charge() != 0 and abs(ipf.pdgId()) != 11 and abs(ipf.pdgId()) != 13 and ipf.pt() > 0.5 )]
        charged_pfs = [ipf for ipf in charged_pfs if (ipf.fromPV() > 1 and abs(ipf.pdgId()) == 211) ]
    charged_pfs  = [ipf for ipf in charged_pfs  if deltaR(ipf, event.the_3lep_cand.hnVisP4()) < dRCone ]
    ch_iso       = sum([ipf.pt() for ipf in charged_pfs]) 
#    for i in charged_pfs: print 'charged hadron, pile up:\t', dRCone, PU, i.dz(), i.dxy(), i.pt(), i.pdgId()
    return ch_iso

def neutralHadronIso(event, dRCone): 
    neutral_pfs = [ipf for ipf in event.pfs if ( ipf.charge() == 0 and abs(ipf.pdgId()) != 11 and abs(ipf.pdgId()) != 13 and abs(ipf.pdgId()) != 22 )]
    neutral_pfs = [ipf for ipf in neutral_pfs if ipf.pt() > 0.5]
    neutral_pfs = [ipf for ipf in neutral_pfs if deltaR(ipf, event.the_3lep_cand.hnVisP4()) < dRCone ]
    neu_iso     = sum([ipf.pt() for ipf in neutral_pfs])
#    for i in neutral_pfs: print 'neutral hadron\t', dRCone, i.dz(), i.dxy(), i.pt(), i.pdgId()
    return neu_iso

def photonIso(event, dRCone): 
    photon_pfs = [ipf for ipf in event.pfs if abs(ipf.pdgId()) == 22]
    photon_pfs = [ipf for ipf in photon_pfs if ipf.pt() > 0.5]
    photon_pfs = [ipf for ipf in photon_pfs if deltaR(ipf, event.the_3lep_cand.hnVisP4()) < dRCone ] 
    ph_iso     = sum([ipf.pt() for ipf in photon_pfs])
#    for i in photon_pfs: print 'photon\t\t', dRCone, i.dz(), i.dxy(), i.pt(), i.pdgId()
    return ph_iso
