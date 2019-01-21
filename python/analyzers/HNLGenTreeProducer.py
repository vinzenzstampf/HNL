import ROOT
import numpy as np
from CMGTools.HNL.analyzers.TreeProducerBase import TreeProducerBase
from PhysicsTools.HeppyCore.utils.deltar import deltaR, bestMatch
from CMGTools.HNL.utils.utils import isAncestor, displacement2D, displacement3D, makeRecoVertex # utility functions
from pdb import set_trace

class HNLTreeProducer(TreeProducerBase):
    '''
    RM: add more info:
    - gen impact parameter  ==> how to do it at gen level?
    - gen met including pu (all neutrinos)
    - di-lepton isolation
    make this iherit from a common reco tree producer, then specialise by lepton flavour
    '''
    def declareVariables(self, setup):
        '''
        '''
        # event variables
        self.bookEvent(self.tree)
        self.var      (self.tree, 'n_cands')
        self.var      (self.tree, 'rho')
        self.var      (self.tree, 'nLeptons')
        self.var      (self.tree, 'nElectrons')
        self.var      (self.tree, 'nMuons')
        
        # reco variables
        self.bookHNL (self.tree, 'hnl')
        self.var     (self.tree, 'hnl_iso03_abs_rhoArea')
        self.var     (self.tree, 'hnl_iso04_abs_rhoArea')
        self.var     (self.tree, 'hnl_iso05_abs_rhoArea')
        self.var     (self.tree, 'hnl_iso03_abs_deltaBeta')
        self.var     (self.tree, 'hnl_iso04_abs_deltaBeta')
        self.var     (self.tree, 'hnl_iso05_abs_deltaBeta')
#        self.var     (self.tree, 'hnl_iso_abs_met')
        self.var     (self.tree, 'hnl_iso03_rel_rhoArea')
        self.var     (self.tree, 'hnl_iso04_rel_rhoArea')
        self.var     (self.tree, 'hnl_iso05_rel_rhoArea')
        self.var     (self.tree, 'hnl_iso03_rel_deltaBeta')
        self.var     (self.tree, 'hnl_iso04_rel_deltaBeta')
        self.var     (self.tree, 'hnl_iso05_rel_deltaBeta')
#        self.var     (self.tree, 'hnl_iso_rel_met')
        
        if   self.cfg_ana.promptLepType == 'ele':
            self.bookEle (self.tree, 'l0')
            self.var(self.tree, 'hlt_Ele27_WPTight_Gsf'        )
            self.var(self.tree, 'hlt_Ele32_WPTight_Gsf'        )
            self.var(self.tree, 'hlt_Ele35_WPTight_Gsf'        )
            self.var(self.tree, 'hlt_Ele115_CaloIdVT_GsfTrkIdT')
            self.var(self.tree, 'hlt_Ele135_CaloIdVT_GsfTrkIdT')

        elif self.cfg_ana.promptLepType == 'mu':
            self.bookMuon(self.tree, 'l0')
            self.var(self.tree, 'hlt_IsoMu24'                   )
            self.var(self.tree, 'hlt_IsoMu27'                   )
            self.var(self.tree, 'hlt_Mu50'                      )
        else:
             print 'ERROR: prompt lepton type non specified or missing! Exit'
             exit(0)

        if self.cfg_ana.L1L2LeptonType == 'mm':
            self.bookMuon(self.tree, 'l1' )
            self.bookMuon(self.tree, 'l2' )
        if self.cfg_ana.L1L2LeptonType == 'ee':
            self.bookEle(self.tree, 'l1'  )
            self.bookEle(self.tree, 'l2'  )
        if self.cfg_ana.L1L2LeptonType == 'em':
            self.bookEle(self.tree, 'l1'  )
            self.bookMuon(self.tree, 'l2' )
        
        # book the matched  gen particle
        self.bookSimpleGenParticle(self.tree, 'l0_gen_match')
        self.bookSimpleGenParticle(self.tree, 'l1_gen_match')
        self.bookSimpleGenParticle(self.tree, 'l2_gen_match')

        # relevant for signal: check if reco matched with gen, save a bool
        self.var(self.tree, 'l0_is_real')
        self.var(self.tree, 'l1_is_real')
        self.var(self.tree, 'l2_is_real')

        # reco primary vertex
        self.var(self.tree, 'pv_x')
        self.var(self.tree, 'pv_y')
        self.var(self.tree, 'pv_z')
        self.var(self.tree, 'pv_xe')
        self.var(self.tree, 'pv_ye')
        self.var(self.tree, 'pv_ze')

        # beamspot
        self.var(self.tree, 'bs_x')
        self.var(self.tree, 'bs_y')
        self.var(self.tree, 'bs_z')
        self.var(self.tree, 'bs_sigma_x')
        self.var(self.tree, 'bs_sigma_y')
        self.var(self.tree, 'bs_sigma_z')
        self.var(self.tree, 'bs_dxdz')
        self.var(self.tree, 'bs_dydz')
        
        # reco HN decay vertex (when present)
        self.bookVertex(self.tree, 'sv')

        # lepton vetoes
        self.var(self.tree, 'pass_e_veto')
        self.var(self.tree, 'pass_m_veto')

        # save vetoing lepton
        self.bookEle(self.tree, 'veto_ele')
        self.bookMuon(self.tree, 'veto_mu')
    
        # invariant masses with vetoing leptons
        self.var(self.tree, 'hnl_m_0Vele')
        self.var(self.tree, 'hnl_m_1Vele')
        self.var(self.tree, 'hnl_m_2Vele')
        self.var(self.tree, 'hnl_m_0Vmu')
        self.var(self.tree, 'hnl_m_1Vmu')
        self.var(self.tree, 'hnl_m_2Vmu')
        
        # gen level particles
        self.bookHNL     (self.tree, 'hnl_gen')
        self.bookParticle(self.tree, 'l0_gen' )
        self.bookParticle(self.tree, 'l1_gen' )
        self.bookParticle(self.tree, 'l2_gen' )
        self.bookParticle(self.tree, 'n_gen'  )

        # gen primary vertex
        self.var(self.tree, 'pv_gen_x')
        self.var(self.tree, 'pv_gen_y')
        self.var(self.tree, 'pv_gen_z')

        # gen HN decay vertex
        self.var(self.tree, 'sv_gen_x')
        self.var(self.tree, 'sv_gen_y')
        self.var(self.tree, 'sv_gen_z')

        # displacements
        self.var(self.tree, 'hnl_2d_gen_disp')
        self.var(self.tree, 'hnl_3d_gen_disp')

        self.var(self.tree, 'hnl_2d_disp')
        self.var(self.tree, 'hnl_3d_disp')

        self.var(self.tree, 'hnl_2d_disp_sig')
        self.var(self.tree, 'hnl_3d_disp_sig')

        # met information
        self.bookExtraMetInfo(self.tree)

        self.var(self.tree, 'Flag_goodVertices')
        self.var(self.tree, 'Flag_globalSuperTightHalo2016Filter')
        self.var(self.tree, 'Flag_HBHENoiseFilter')
        self.var(self.tree, 'Flag_HBHENoiseIsoFilter')
        self.var(self.tree, 'Flag_EcalDeadCellTriggerPrimitiveFilter')
        self.var(self.tree, 'Flag_BadPFMuonFilter')
        self.var(self.tree, 'Flag_BadChargedCandidateFilter')
        self.var(self.tree, 'Flag_eeBadScFilter')
        self.var(self.tree, 'Flag_ecalBadCalibFilter')
#        self.var(self.tree, 'Flag_any_met_filters')
        
        # jet information
        self.bookJet(self.tree, 'j1' , fill_extra=False)
        self.bookJet(self.tree, 'j2' , fill_extra=False)
        self.bookJet(self.tree, 'bj1', fill_extra=False)
        self.bookJet(self.tree, 'bj2', fill_extra=False)
        self.var(self.tree, 'htj' )
        self.var(self.tree, 'htbj')
        self.var(self.tree, 'nj'  )
        self.var(self.tree, 'nbj' )
        
        # LHE weight
        self.var(self.tree, 'lhe_weight')

    def process(self, event):
        '''
        '''
        self.readCollections(event.input)
        self.tree.reset()
        self.event = event
        # event variables 
        self.fillEvent(self.tree, event)
        self.fill     (self.tree, 'rho'    , event.rho)

        # reco HNL
        # output of MC analysis ONLY FOR SIGNAL
        if hasattr(event, 'the_hnl'):
            self.fillHNL     (self.tree, 'hnl_gen', event.the_hnl      )
            self.fillParticle(self.tree, 'l0_gen' , event.the_hnl.l0() )
            self.fillParticle(self.tree, 'l1_gen' , event.the_hnl.l1() )
            self.fillParticle(self.tree, 'l2_gen' , event.the_hnl.l2() )
            self.fillParticle(self.tree, 'n_gen'  , event.the_hnl.met())

        # true primary vertex
        if hasattr(event, 'the_hnl'):
            self.fill(self.tree, 'pv_gen_x', event.the_hn.vx())
            self.fill(self.tree, 'pv_gen_y', event.the_hn.vy())
            self.fill(self.tree, 'pv_gen_z', event.the_hn.vz())

        # true HN decay vertex
        if hasattr(event, 'the_hn'):
            self.fill(self.tree, 'sv_gen_x', event.the_hn.lep1.vertex().x()) # don't use the final lepton to get the vertex from!
            self.fill(self.tree, 'sv_gen_y', event.the_hn.lep1.vertex().y()) # don't use the final lepton to get the vertex from!
            self.fill(self.tree, 'sv_gen_z', event.the_hn.lep1.vertex().z()) # don't use the final lepton to get the vertex from!

            # displacements
            self.fill(self.tree, 'hnl_2d_gen_disp', displacement2D(event.the_hn.lep1, event.the_hn))
            self.fill(self.tree, 'hnl_3d_gen_disp', displacement3D(event.the_hn.lep1, event.the_hn))

        # HLT bits & matches
        trig_list = [trig.name for trig in event.trigger_infos if trig.fired]
        if self.cfg_ana.promptLepType == 'ele':
            self.fill(self.tree, 'hlt_Ele27_WPTight_Gsf'               , any('HLT_Ele27_WPTight_Gsf'                 in name for name in trig_list))
            self.fill(self.tree, 'hlt_Ele32_WPTight_Gsf'               , any('HLT_Ele32_WPTight_Gsf'                 in name for name in trig_list))
            self.fill(self.tree, 'hlt_Ele35_WPTight_Gsf'               , any('HLT_Ele35_WPTight_Gsf'                 in name for name in trig_list))
            self.fill(self.tree, 'hlt_Ele115_CaloIdVT_GsfTrkIdT'       , any('HLT_Ele115_CaloIdVT_GsfTrkIdT'         in name for name in trig_list))
            self.fill(self.tree, 'hlt_Ele135_CaloIdVT_GsfTrkIdT'       , any('HLT_Ele135_CaloIdVT_GsfTrkIdT'         in name for name in trig_list))
        if self.cfg_ana.promptLepType == 'mu':
            self.fill(self.tree, 'hlt_IsoMu24'                          , any('HLT_IsoMu24'                           in name for name in trig_list))
            self.fill(self.tree, 'hlt_IsoMu27'                          , any('HLT_IsoMu27'                           in name for name in trig_list))
            self.fill(self.tree, 'hlt_Mu50'                             , any('HLT_Mu50'                              in name for name in trig_list))

#               if len(event.cleanJets )>0: self.fillJet(self.tree, 'j1' , event.cleanJets [0], fill_extra=False)
        if len(event.cleanJets )>1: self.fillJet(self.tree, 'j2' , event.cleanJets [1], fill_extra=False)
        if len(event.cleanBJets)>0: self.fillJet(self.tree, 'bj1', event.cleanBJets[0], fill_extra=False)
        if len(event.cleanBJets)>1: self.fillJet(self.tree, 'bj2', event.cleanBJets[1], fill_extra=False)

        self.fill(self.tree, 'htj' , event.HT_cleanJets   )
        self.fill(self.tree, 'htbj', event.HT_bJets       )
        self.fill(self.tree, 'nj'  , len(event.cleanJets) )
        self.fill(self.tree, 'nbj' , len(event.cleanBJets))

        # LHE weight
        self.fill(self.tree, 'lhe_weight', np.sign(getattr(event, 'LHE_originalWeight', 1.)))
                
        self.fillTree(event)


