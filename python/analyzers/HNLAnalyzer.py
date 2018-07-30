'''
This is the main analyzer going through data and trying to identify HNL->3L events.
'''

import ROOT
from itertools import product, combinations
from math import sqrt, pow

from PhysicsTools.Heppy.analyzers.core.Analyzer      import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle    import AutoHandle
from PhysicsTools.Heppy.physicsobjects.GenParticle   import GenParticle
from PhysicsTools.Heppy.physicsobjects.Muon          import Muon
from PhysicsTools.Heppy.physicsobjects.Electron      import Electron
from PhysicsTools.Heppy.physicsobjects.PhysicsObject import PhysicsObject
from CMGTools.HNL.utils.utils                        import isAncestor, displacement2D, displacement3D, makeRecoVertex
from PhysicsTools.HeppyCore.utils.deltar             import deltaR, deltaPhi

from CMGTools.HNL.physicsobjects.DiLepton      import DiLepton
from CMGTools.HNL.physicsobjects.DisplacedMuon import DisplacedMuon
from pdb import set_trace

# load custom library to ROOT. This contains the kinematic vertex fitter class
ROOT.gSystem.Load('libCMGToolsHNL')
from ROOT import HNLKinematicVertexFitter as VertexFitter

class HNLAnalyzer(Analyzer):
    '''
    '''

    def declareHandles(self):
        super(HNLAnalyzer, self).declareHandles()
        
        self.handles['ele']    = AutoHandle(('slimmedElectrons', '','PAT'), 'std::vector<pat::Electron>')
        self.handles['muons']    = AutoHandle(('slimmedMuons'                 ,'','PAT' ), 'std::vector<pat::Muon>'                        )
        self.handles['dsamuons'] = AutoHandle(('displacedStandAloneMuons'     ,'','RECO'), 'std::vector<reco::Track>'                      )
        self.handles['dgmuons' ] = AutoHandle(('displacedGlobalMuons'         ,'','RECO'), 'std::vector<reco::Track>'                      )
        self.handles['pvs']      = AutoHandle(('offlineSlimmedPrimaryVertices','','PAT' ), 'std::vector<reco::Vertex>'                     )
        self.handles['svs']      = AutoHandle(('slimmedSecondaryVertices'     ,'','PAT' ), 'std::vector<reco::VertexCompositePtrCandidate>')
        self.handles['beamspot'] = AutoHandle(('offlineBeamSpot'              ,'','RECO'), 'reco::BeamSpot'                                )
        self.handles['met']      = AutoHandle(('slimmedMETs'                  ,'','PAT' ), 'std::vector<pat::MET>'                         )

    def assignVtx(self, particles, vtx):    
        for ip in particles:
            ip.associatedVertex = vtx

    def beginLoop(self, setup):
        super(HNLAnalyzer, self).beginLoop(setup)
        self.counters.addCounter('HNL')
        count = self.counters.counter('HNL')
        count.register('all events')
        count.register('good gen')
        count.register('> 0 di-muon')
        count.register('> 0 di-muon + vtx')
        count.register('dimuons')

        # initiate the VertexFitter
        self.vtxfit = VertexFitter()

        # create a std::vector<RecoChargedCandidate> to be passed to the fitter
        self.tofit = ROOT.std.vector('reco::RecoChargedCandidate')()

    def buildDisplacedMuons(self, collection):
        muons = [DisplacedMuon(mm, collection) for mm in collection]
        return muons

    def process(self, event):
        self.readCollections(event.input)
        self.counters.counter('HNL').inc('all events')

        #####################################################################################
        # produce collections and map our objects to convenient Heppy objects
        #####################################################################################

        event.ele         = map(Electron, self.handles['ele'].product())
        # make muon collections
        event.muons       = map(Muon, self.handles['muons'].product())
        event.dsamuons    = self.buildDisplacedMuons(self.handles['dsamuons'].product())
        event.dgmuons     = self.buildDisplacedMuons(self.handles['dgmuons' ].product())

        # make vertex objects 
        event.pvs         = self.handles['pvs'     ].product()
        event.svs         = self.handles['svs'     ].product()
        event.beamspot    = self.handles['beamspot'].product()

        # make met object
        event.met         = self.handles['met'].product().at(0)

        # assign to the leptons the primary vertex, will be needed to compute a few quantities
        # TODO! understand exactly to which extent it is reasonable to assign the PV to *all* leptons
        #        regardless whether they're displaced or not
        if len(event.pvs):
            myvtx = event.pvs[0]
        else:
            myvtx = event.beamspot
        
        self.assignVtx(event.muons, myvtx)
        self.assignVtx(event.ele, myvtx)

        event.n_ele = len(event.ele)

        # PROMPT ANALYSIS # SO FAR ONLY ELECTRONS

       # MUONS TODO
        mu_cand = []
        matchable_mu = [mu for mu in event.muons] 
        # selection
        mu_sel_eta = 2.4; mu_sel_pt = 3; mu_sel_vtx = 0.2 
        # match collections
        matchable_mu_sel_pt = [mu for mu in matchable_mu if (mu.pt() > mu_sel_pt)] 
        matchable_mu_sel_eta = [mu for mu in matchable_mu if (abs(mu.eta()) < mu_sel_eta)] 
        matchable_mu_sel_id = [mu for mu in matchable_mu if (mu.looseId() == True)] 
        # https://github.com/rmanzoni/cmgtools-lite/blob/825_HTT/H2TauTau/python/proto/analyzers/TauEleAnalyzer.py#L193
        matchable_mu_sel_vtx = [mu for mu in matchable_mu if abs(mu.dz()) < mu_sel_vtx] # TODO what about dxy component ?
        # https://github.com/rmanzoni/cmgtools-lite/blob/825_HTT/H2TauTau/python/proto/analyzers/TauEleAnalyzer.py#L104
        mu_cand = [mu for mu in matchable_mu if (mu in matchable_mu_sel_pt and mu in matchable_mu_sel_eta and mu in matchable_mu_sel_id and mu in matchable_mu_sel_vtx)]

        # ELECTRONS
        ele_cand = []
#        matchable_ele = [ele for ele in event.ele]
        # selection
#        ele_sel_eta = 2.5; ele_sel_pt = 3; ele_sel_vtx = 0.2 
        # match collections
#        matchable_ele_sel_pt = [ele for ele in matchable_ele if (ele.pt() > ele_sel_pt)] 
#        matchable_ele_sel_eta = [ele for ele in matchable_ele if (abs(ele.eta()) < ele_sel_eta)] 
#        matchable_ele_sel_id = [ele for ele in matchable_ele if (ele.mvaIDRun2('NonTrigSpring15MiniAOD', 'POG90') == True)] 
        # https://github.com/rmanzoni/cmgtools-lite/blob/825_HTT/H2TauTau/python/proto/analyzers/TauEleAnalyzer.py#L193
#        matchable_ele_sel_vtx = [ele for ele in matchable_ele if abs(ele.dz()) < ele_sel_vtx] # TODO what about dxy component ?
        
#        ele_cand = [ele for ele in matchable_ele if (ele in matchable_ele_sel_pt and ele in matchable_ele_sel_eta and ele in matchable_ele_sel_id and ele in matchable_ele_sel_vtx)]
        
        prompt_cand = ele_cand + mu_cand
        the_prompt_cand = None
        # EVALUATING THE PROMPT SELECTION: EFF / PUR
        event.prompt_ana_success = -99 # NO RECO FOUND
#        if not len(prompt_cand): return False
        if len(prompt_cand): 
        # selection: pick candidate with highest pt 
        # there must be something better; maybe if both are matched check some additional stuff
            the_prompt_cand = sorted(prompt_cand, key = lambda lep: lep.pt(), reverse = True)[0]
#            event.the_prompt_cand = the_prompt_cand # TODO WHEN TRIGGER STUFF IS DONE, MOVE THIS LINE AFTER TRIGGERS
            # REMOVING PROMPT LEPTON FROM MATCHES 
            if the_prompt_cand in mu_cand:
                event.muons.remove(the_prompt_cand)
                if hasattr(event.the_hnl.l0().bestmatch, 'physObj'):
                    if  the_prompt_cand.physObj == event.the_hnl.l0().bestmatch.physObj:
                        event.prompt_ana_success = 1
                else: event.prompt_ana_success = -13 # FAKE MUONS 

        # TRIGGER MATCHING
        # match only if the trigger fired
        event.fired_triggers = [info.name for info in getattr(event, 'trigger_infos', []) if info.fired]

        # trigger matching
        if hasattr(self.cfg_ana, 'trigger_match') and len(self.cfg_ana.trigger_match.keys())>0:
                                   
            for ele in the_prompt_cand:
                
#                triplet.hltmatched = [] # initialise to no match
                ele.hltmatched = [] # initialise to no match
                
#               triplet.trig_objs = OrderedDict()
#               triplet.trig_objs[1] = [] # initialise to no trigger objct matches
#               triplet.trig_objs[2] = [] # initialise to no trigger objct matches
#               triplet.trig_objs[3] = [] # initialise to no trigger objct matches
                ele.trig_objs = OrderedDict()
                ele.trig_objs[1] = [] # initialise to no trigger objct matches
    
#               triplet.trig_matched = OrderedDict()
#               triplet.trig_matched[1] = False # initialise to no match
#               triplet.trig_matched[2] = False # initialise to no match
#               triplet.trig_matched[3] = False # initialise to no match
                ele.trig_matched = OrderedDict()
                ele.trig_matched[1] = False # initialise to no match

#               triplet.best_trig_match = OrderedDict()
#               triplet.best_trig_match[1] = OrderedDict()
#               triplet.best_trig_match[2] = OrderedDict()
#               triplet.best_trig_match[3] = OrderedDict()
                ele.best_trig_match = OrderedDict()
                ele.best_trig_match[1] = OrderedDict()

                # add all matched objects to each muon
                for info in event.trigger_infos:
                                    
                    mykey = '_'.join(info.name.split('_')[:-1])

                    # start with simple matching
#                    these_objects1 = sorted([obj for obj in info.objects if deltaR(triplet.mu1(), obj)<0.15], key = lambda x : deltaR(x, triplet.mu1()))
#                    these_objects2 = sorted([obj for obj in info.objects if deltaR(triplet.mu2(), obj)<0.15], key = lambda x : deltaR(x, triplet.mu2()))
#                    these_objects3 = sorted([obj for obj in info.objects if deltaR(triplet.mu3(), obj)<0.15], key = lambda x : deltaR(x, triplet.mu3()))
                    these_objects = sorted([obj for obj in info.objects if deltaR(ele, obj)<0.15], key = lambda x : deltaR(x, ele))

#                   triplet.trig_objs[1] += these_objects1
#                   triplet.trig_objs[2] += these_objects2
#                   triplet.trig_objs[3] += these_objects3
                    ele.trig_objs[1] += these_objects

                    # get the set of trigger types from the cfg 
                    trigger_types_to_match = self.cfg_ana.trigger_match[mykey][1]
                    
                    # list of tuples of matched objects
                    good_matches = []

                    # initialise the matching to None
#                    triplet.best_trig_match[1][mykey] = None
#                    triplet.best_trig_match[2][mykey] = None
#                    triplet.best_trig_match[3][mykey] = None
                    ele.best_trig_match[1][mykey] = None

                    # investigate all the possible matches (eles, pairs or singlets)
                   # for to1, to2, to3 in product(these_objects1, these_objects2, these_objects3):
                    for t_o in these_objects:
                        # avoid double matches!
                        # if to1==to2 or to1==to3 or to2==to3:
                        #    continue

                        # intersect found trigger types to desired trigger types
                        itypes = Counter()
                        for ikey in trigger_types_to_match.keys():
                           # itypes[ikey] = sum([1 for iobj in [to1, to2, to3] if iobj.triggerObjectTypes()[0]==ikey])
                            itypes[ikey] = sum([1 for iobj in t_o if iobj.triggerObjectTypes()[0]==ikey])
                                            
                        # all the types to match are matched then assign the 
                        # corresponding trigger object to each ele
                        if itypes & trigger_types_to_match == trigger_types_to_match:
                            #good_matches.append((to1, to2, to3))
                            good_matches.append(t_o)
                    
                    
                    if len(good_matches):
                        # good_matches.sort(key = lambda x : deltaR(x[0], ele.mu1()) + deltaR(x[1], ele.mu2()) + deltaR(x[2], ele.mu3()))        
                        good_matches.sort(key = lambda x : deltaR(x, ele))        

                        # ONLY for HLT_DoubleMu3_Trk_Tau3mu
                        # it might happen that more than one combination of trk mu mu is found,
                        # make sure that the online 3-body mass cut is satisfied by the matched objects
#                       if mykey == 'HLT_DoubleMu3_Trk_Tau3mu':
#                          
#                          good_matches_tmp = []
#                          
#                          for im in good_matches:
#                              p4_1 = ROOT.TLorentzVector()
#                              p4_2 = ROOT.TLorentzVector()
#                              p4_3 = ROOT.TLorentzVector()
#
#                              p4_1.SetPtEtaPhiM(im[0].pt(), im[0].eta(), im[0].phi(), 0.10565999895334244)                        
#                              p4_2.SetPtEtaPhiM(im[1].pt(), im[1].eta(), im[1].phi(), 0.10565999895334244)
#                              p4_3.SetPtEtaPhiM(im[2].pt(), im[2].eta(), im[2].phi(), 0.10565999895334244)
#                                      
#                              totp4 = p4_1 + p4_2 + p4_3
#                              
#                              if totp4.M()>1.6 and totp4.M()<2.02:
#                                  good_matches_tmp.append(im)
#                          
#                          good_matches = good_matches_tmp
#                          
#                      ele.best_trig_match[1][mykey] = good_matches[0][0] if len(good_matches) and len(good_matches[0])>0 else None
#                      ele.best_trig_match[2][mykey] = good_matches[0][1] if len(good_matches) and len(good_matches[0])>1 else None
#                      ele.best_trig_match[3][mykey] = good_matches[0][2] if len(good_matches) and len(good_matches[0])>2 else None
                
                # iterate over the path:filters dictionary
                #     the filters MUST be sorted correctly: i.e. first filter in the dictionary 
                #     goes with the first muons and so on
                for k, vv in self.cfg_ana.trigger_match.iteritems():

                    if not any(k in name for name in event.fired_triggers):
                         continue
                    
                    v = vv[0]
                                                                 
                    for ii, filters in enumerate(v):
                        if not ele.best_trig_match[ii+1][k]:
                            continue
                        if set([filters]) & set(ele.best_trig_match[ii+1][k].filterLabels()):
                            ele.trig_matched[ii+1] = True                 
                    
                    ismatched = sum(ele.trig_matched.values())            
                                
                    if len(v) == ismatched:
                        ele.hltmatched.append(k)

#            seltau3mu = [triplet for triplet in seltau3mu if len(triplet.hltmatched)>0]
            the_prompt_cand = [ele for ele in the_prompt_cand if len(ele.hltmatched)>0]
            
#            if len(the_prompt_cand) == 0:
#                return False #TODO  UNCOMMENT THIS IN FINAL VERSION
#            self.counters.counter('Tau3Mu').inc('trigger matched')

#        event.seltau3mu = seltau3mu

#        event.tau3mu = self.bestTriplet(event.seltau3mu)                        
        event.the_prompt_cand = the_prompt_cand

#        return True #TODO  UNCOMMENT THIS IN FINAL VERSION

        
        self.counters.counter('HNL').inc('good gen')

        # some simple preselection
        event.muons    = [imu for imu in event.muons    if imu.pt()>3.]
        event.dsamuons = [imu for imu in event.dsamuons if imu.pt()>3.]
        event.dgmuons  = [imu for imu in event.dgmuons  if imu.pt()>3.]

        # create all the possible di-muon pairs out of the three different collections
        dimuons = combinations(event.muons + event.dsamuons + event.dgmuons, 2)
        
        dimuons = [(mu1, mu2) for mu1, mu2 in dimuons if deltaR(mu1, mu2)>0.01]
        
        self.counters.counter('HNL').inc('> 0 di-muon + vtx')
        
        ########################################################################################
        # select only dimuon pairs with mutual vertices (surviving the kinematic vertex fitter)
        ########################################################################################
        
#         print '\ngen truth'
#         print '\tl0\n', event.the_hnl.l0()
#         print '\tl1\n', event.the_hn.lep1, '\n', event.the_hnl.l1().bestmatch
#         print '\tl2\n', event.the_hn.lep2, '\n', event.the_hnl.l2().bestmatch
#         print '\tgen  sv x=%.2f y=%.2f z=%.2f' %(event.the_hn.lep1.vertex().x(), event.the_hn.lep1.vertex().y(), event.the_hn.lep1.vertex().z())
#         print '\treco sv x=%.2f y=%.2f z=%.2f' %(event.recoSv.x(), event.recoSv.y(), event.recoSv.z())
#         print '\tgen 2d displacement %.3f' %displacement2D(event.the_hn.lep1, event.the_hn)
#         print '\treco 2d displacement %.3f' %(-99. if not hasattr(event, 'recoSv') else event.recoSv.disp2DFromBS.value())
#         print '\t3d displacement %.3f' %displacement3D(event.the_hn.lep1, event.the_hn)
#         print '\tcos             %.6f' %event.the_hnl.cos()
#         print '\tllmass          %.3f' %event.the_hnl.mass12()

        dimuonsvtx = []
        for index, pair in enumerate(dimuons):
            if pair[0]==pair[1]:
                continue
            
            print pair[0]
            print pair[1]
            
            self.tofit.clear()
            for il in pair:
                if il.muonBestTrack().isNull():       # check that the track is valid, there are photons around too!
                    continue
                ic = ROOT.reco.RecoChargedCandidate() # instantiate a dummy RecoChargedCandidate
                ic.setCharge(il.charge())             # assign the correct charge
                ic.setTrack(il.muonBestTrack())
                # if the reco particle is a displaced thing, it does not have the p4() method, so let's build it 
                myp4 = ROOT.Math.LorentzVector('<ROOT::Math::PxPyPzE4D<double> >')(ic.track().px(), ic.track().py(), ic.track().pz(), sqrt(il.mass()**2 + ic.track().px()**2 + ic.track().py()**2 + ic.track().pz()**2))
                ic.setP4(myp4)                        # assign the correct p4
                self.tofit.push_back(ic)

            # further sanity check: two *distinct* tracks
            if self.tofit.size() == 2 and self.tofit[0].track() != self.tofit[1].track():
                svtree = self.vtxfit.Fit(self.tofit) # the actual vertex fitting!
                if not svtree.get().isEmpty() and svtree.get().isValid(): # check that the vertex is good
                    svtree.movePointerToTheTop()
                    sv = svtree.currentDecayVertex().get()
                    dimuonsvtx.append(DiLepton(pair, makeRecoVertex(sv, kinVtxTrkSize=2), myvtx, event.beamspot))
#                     print dimuonsvtx[-1]
            else:
                print 'FAILED!'
            
        self.counters.counter('HNL').inc('> 0 di-muon + vtx')

        
        ########################################################################################
        # candidate choice by different criteria
        ########################################################################################

#         event.hnl_minchi2     = None if not len(dimuonsvtx) else sorted(dimuonsvtx, key = lambda x : (x.isSS(),  x.chi2()                        ), reverse=False)[0]
#         event.hnl_maxpt       = None if not len(dimuonsvtx) else sorted(dimuonsvtx, key = lambda x : (x.isSS(), -x.pt()                          ), reverse=False)[0]
#         event.hnl_mindr       = None if not len(dimuonsvtx) else sorted(dimuonsvtx, key = lambda x : (x.isSS(),  x.dr()                          ), reverse=False)[0]
#         event.hnl_maxdphi     = None if not len(dimuonsvtx) else sorted(dimuonsvtx, key = lambda x : (x.isSS(), -x.dphi()                        ), reverse=False)[0]
#         event.hnl_mindeta     = None if not len(dimuonsvtx) else sorted(dimuonsvtx, key = lambda x : (x.isSS(),  x.deta()                        ), reverse=False)[0]
#         event.hnl_maxdisp2dbs = None if not len(dimuonsvtx) else sorted(dimuonsvtx, key = lambda x : (x.isSS(), -x.disp2DFromBS()                ), reverse=False)[0]
#         event.hnl_maxdisp2dpv = None if not len(dimuonsvtx) else sorted(dimuonsvtx, key = lambda x : (x.isSS(), -x.disp2DFromPV()                ), reverse=False)[0]
#         event.hnl_maxdisp3dpv = None if not len(dimuonsvtx) else sorted(dimuonsvtx, key = lambda x : (x.isSS(), -x.disp3DFromPV()                ), reverse=False)[0]
#         event.hnl_maxdls2dbs  = None if not len(dimuonsvtx) else sorted(dimuonsvtx, key = lambda x : (x.isSS(), -x.disp2DFromBSSignificance()    ), reverse=False)[0]
#         event.hnl_maxdls2dpv  = None if not len(dimuonsvtx) else sorted(dimuonsvtx, key = lambda x : (x.isSS(), -x.disp2DFromPVSignificance()    ), reverse=False)[0]
#         event.hnl_maxdls3dpv  = None if not len(dimuonsvtx) else sorted(dimuonsvtx, key = lambda x : (x.isSS(), -x.disp3DFromPVSignificance()    ), reverse=False)[0]
        event.hnl_maxcos      = None if not len(dimuonsvtx) else sorted(dimuonsvtx, key = lambda x : (x.isSS(), -x.cosTransversePointingAngleBS()), reverse=False)[0]

        event.dimuonsvtx = dimuonsvtx

#         print 'hnl_minchi2    ', event.hnl_minchi2    
#         print 'hnl_maxpt      ', event.hnl_maxpt      
#         print 'hnl_mindr      ', event.hnl_mindr      
#         print 'hnl_maxdphi    ', event.hnl_maxdphi    
#         print 'hnl_mindeta    ', event.hnl_mindeta    
#         print 'hnl_maxdisp2dbs', event.hnl_maxdisp2dbs
#         print 'hnl_maxdisp2dpv', event.hnl_maxdisp2dpv
#         print 'hnl_maxdisp3dpv', event.hnl_maxdisp3dpv
#         print 'hnl_maxdls2dbs ', event.hnl_maxdls2dbs 
#         print 'hnl_maxdls2dpv ', event.hnl_maxdls2dpv 
#         print 'hnl_maxdls3dpv ', event.hnl_maxdls3dpv 
#         print 'hnl_maxcos     ', event.hnl_maxcos     
 
#         import pdb ; pdb.set_trace()        


#         if not len(dimuonsvtx):
#             import pdb ; pdb.set_trace()        

#
#
#
#
#
#       
#        #####################################################################################
#        # MUCO
#        # Concatenate all Muon Reconstructions:
#        # Create an array of DisplacedMuon objects, 
#        # summarizing all sMu and dSAMus into a single array, 
#        # avoid redundancies with dR<0.2
#        #####################################################################################
#        # Merge Reco Muons
#        # Create an array of DisplacedMuon objects, summarizing all sMu and dSAMus into a single array, while avoiding redundancies through dR<0.2
#        dMus = []
#        dxy_cut = 10 # cut selection for sMu / dSAMu in cm
#        event.n_sMuOnly = 0
#        event.n_dSAMuOnly = 0
#        event.n_sMuRedundant = 0
#        event.n_dSAMuRedundant = 0
#        for smu in event.sMu:
#            matches = []
#            matches = [dsa for dsa in event.dSAMu if (deltaR(smu,dsa)<0.2 and abs((smu.pt()-dsa.pt())/smu.pt())<0.3)] 
#            if len(matches) == 0:
#                dmu = smu
#                dmu.reco = 1 # sMu = 1, dSAMu = 2
#                dmu.redundancy = 0
#                dMus.append(dmu)
#                event.n_sMuOnly += 1
#            else:
#                bestmatch = sorted(matches, key = lambda dsa: deltaR(smu,dsa), reverse = False)[0] 
##                 import pdb ; pdb.set_trace()
#                if smu.vx()**2 + smu.vy()**2 < dxy_cut**2:
##                 if smu.dxy() < dxy_cut:
#                    dmu = smu
#                    dmu.reco = 1 # sMu = 1, dSAMu = 2 
#                    dmu.redundancy = len(matches)
#                    dMus.append(dmu)
#                    event.n_sMuRedundant += 1
#                else:
#                    dmu = bestmatch
#                    dmu.reco = 2 # sMu = 1, dSAMu = 2
#                    dmu.redundancy = len(matches)
#                    dMus.append(dmu)
#                    event.n_dSAMuRedundant += 1
#
#                event.dSAMu.remove(bestmatch)    
#
#        for dsa in event.dSAMu:
#            dmu = dsa
#            dmu.reco = 2 # sMu = 1, dSAMu = 2
#            dmu.redundancy = 0 
#            dMus.append(dmu)
#            event.n_dSAMuOnly += 1
#       
#        event.n_dMu = len(dMus) # important to understand how well the "Merge Reco Muons" process went. 
#       
#        #####################################################################################
#        # Qualify the performance of the MUCO
#        #####################################################################################
#        event.flag_MUCOsuccess = False    
#        l1matched=False
#        l2matched=False
#        if len(dMus) > 1 and hasattr(event.the_hnl.l1().bestmatch, 'physObj') and hasattr(event.the_hnl.l2().bestmatch,'physObj'):
#            for dmu in dMus:
#                if dmu.physObj == event.the_hnl.l1().bestmatch.physObj:
#                    l1matched = True
#                if dmu.physObj == event.the_hnl.l2().bestmatch.physObj:
#                    l2matched = True
#        if l1matched and l2matched:
#            event.flag_MUCOsuccess = True
#
#        #####################################################################################
#        # select only events with good gen events
#        #####################################################################################
#        if not( abs(event.the_hnl.l1().pdgId())==13   and \
#                abs(event.the_hnl.l2().pdgId())==13   and \
#                abs(event.the_hnl.l1().eta())   < 2.4 and \
#                abs(event.the_hnl.l2().eta())   < 2.4 and \
#                abs(event.the_hnl.l0().eta())   < 2.4): 
#            return False
#
#        self.counters.counter('HNL').inc('good gen')
#
#        # #####################################################################################
#        # TODO: Preselection for the reco muons
#        # #####################################################################################
#
#
#
#
#        #####################################################################################
#        # collect all muon pairs
#        #####################################################################################
#        event.pairs = [pair for pair in combinations(dMus,2)] 
#        event.n_pairs = len(event.pairs)
#        event.flag_IsThereTHEDimuon = False
#
#        event.n_dimuon = 0
#        if len(event.pairs) > 0:
#            self.counters.counter('HNL').inc('pairs')
#
#            ########################################################################################
#            # select only dimuon pairs with mutual vertices (surviving the kinematic vertex fitter)
#            ########################################################################################
#            dimuons = []
#            for index, pair in enumerate(event.pairs):
#                if not pair[0]==pair[1]:
#                    self.tofit.clear()
#                    for il in pair:
#                        # if the reco particle is a displaced thing, it does not have the p4() method, so let's build it 
##                         myp4 = ROOT.Math.LorentzVector('<ROOT::Math::PxPyPzE4D<double> >')(il.px(), il.py(), il.pz(), sqrt(il.mass()**2 + il.px()**2 + il.py()**2 + il.pz()**2))
#                        ic = ROOT.reco.RecoChargedCandidate() # instantiate a dummy RecoChargedCandidate
#                        ic.setCharge(il.charge())           # assign the correct charge
##                         ic.setP4(myp4)                      # assign the correct p4
#                        
#                        if hasattr(il, 'isStandAloneMuon') and hasattr(il, 'isGlobalMuon'):
#                            if il.isStandAloneMuon() and not il.isGlobalMuon():
#                                ic.setTrack(il.standAloneMuon())
#                            else:
#                                ic.setTrack(il.standAloneMuon())
##                                 ic.setTrack(il.track())                                
#                        else:
#                            ic.setTrack(il.track())
#
#                        myp4 = ROOT.Math.LorentzVector('<ROOT::Math::PxPyPzE4D<double> >')(ic.track().px(), ic.track().py(), ic.track().pz(), sqrt(il.mass()**2 + ic.track().px()**2 + ic.track().py()**2 + ic.track().pz()**2))
#                        ic.setP4(myp4)                      # assign the correct p4
#
#                        # if il.reco == 1: # sMu = 1, dSAMu = 2
#                            # ic.setTrack(il.outerTrack())             # set the correct TrackRef
#                        # if il.reco == 2: # sMu = 1, dSAMu = 2
#                            # ic.setTrack(il.physObj.track())             # set the correct TrackRef
#                        if ic.track().isNonnull():          # check that the track is valid, there are photons around too!
#                            self.tofit.push_back(ic)
#                    if index==0: import pdb ; pdb.set_trace()
#                    # further sanity check: two *distinct* tracks
#                    if self.tofit.size() == 2 and self.tofit[0].track() != self.tofit[1].track():
#                        svtree = self.vtxfit.Fit(self.tofit) # the actual vertex fitting!
#                        if not svtree.get().isEmpty() and svtree.get().isValid(): # check that the vertex is good
#                            svtree.movePointerToTheTop()
#                            sv = svtree.currentDecayVertex().get()
#                            dimuons.append(DiMuon(pair, makeRecoVertex(sv, kinVtxTrkSize=2)))
#                    if index==0: import pdb ; pdb.set_trace()
#
#            #####################################################################################
#            # Check whether the correct dimuon is part of the collection dimuons
#            #####################################################################################
#            # if abs(event.the_hnl.l1().bestmatch.pt() - 20.650056)<0.001:
#                # set_trace()
#            if len(dimuons) > 0 and hasattr(event.the_hnl.l1().bestmatch, 'physObj') and hasattr(event.the_hnl.l2().bestmatch,'physObj'):
#                for dimu in dimuons:
#                    dMu1 = dimu.pair[0]
#                    dMu2 = dimu.pair[1]
#                    if (dMu1.physObj == event.the_hnl.l1().bestmatch.physObj or dMu1.physObj == event.the_hnl.l2().bestmatch.physObj) and \
#                       (dMu2.physObj == event.the_hnl.l1().bestmatch.physObj or dMu2.physObj == event.the_hnl.l2().bestmatch.physObj):
#                        event.flag_IsThereTHEDimuon = True
#
#            if event.flag_MUCOsuccess and not event.flag_IsThereTHEDimuon:
#                import pdb ; pdb.set_trace()
#
#            #####################################################################################
#            # select the best dimuon pairs 
#            #####################################################################################
#            if len(dimuons) > 0:
#                self.counters.counter('HNL').inc('dimuons')
#                
#                event.n_dimuon = len(dimuons)
#                 
#                # select the dimuon with lowest vertex fit chi2 as the HNL dimuon candidate
#                dimuonChi2 = sorted(dimuons, key = lambda x: (x.isSS(),x.chi2()), reverse = False)[0] 
#                event.dimuonChi2 = dimuonChi2
#                event.dMu1Chi2 = sorted(dimuonChi2.pair, key = lambda x: x.pt(), reverse = True)[0]
#                event.dMu2Chi2 = sorted(dimuonChi2.pair, key = lambda x: x.pt(), reverse = False)[0] 
#                
#                # select the dimuon with largest displacement
#                dimuonDxy = sorted(dimuons, key = lambda x: (x.isSS(),x.dxy()), reverse = True)[0] 
#                event.dimuonDxy = dimuonDxy
#                event.dMu1Dxy = sorted(dimuonDxy.pair, key = lambda x: x.pt(), reverse = True)[0]
#                event.dMu2Dxy = sorted(dimuonDxy.pair, key = lambda x: x.pt(), reverse = False)[0] 
#
#                # select leptons ito added momenta's pt
#                dimuonMaxPt = sorted(dimuons, key = lambda x: (x.pair[0].p4() + x.pair[1].p4()).pt(), reverse = True)[0] 
#                event.dimuonMaxPt = dimuonMaxPt
#                event.dMu1MaxPt = sorted(dimuonMaxPt.pair, key = lambda x: x.pt(), reverse = True)[0] 
#                event.dMu2MaxPt = sorted(dimuonMaxPt.pair, key = lambda x: x.pt(), reverse = False)[0]
#
#                # select closest leptons ito dr
#                dimuonMinDr12 = sorted(dimuons, key = lambda x: deltaR(x.pair[0],x.pair[1]), reverse = False)[0]
#                event.dimuonMinDr12 = dimuonMinDr12
#                event.dMu1MinDr12 = sorted(dimuonMinDr12.pair, key = lambda x: x.pt(), reverse = True)[0] 
#                event.dMu2MinDr12 = sorted(dimuonMinDr12.pair, key = lambda x: x.pt(), reverse = False)[0]
#
#                # select leptons farthest to l0 ito dr of added momenta (l1+l2)
#                # DEPENDENT ON GEN INFO
#                dimuonMaxDr0a12 = sorted(dimuons, key = lambda x: deltaR(x.pair[0].p4()+x.pair[1].p4(),event.the_hnl.l0().p4()), reverse = True)[0]
#                event.dimuonMaxDr0a12 = dimuonMaxDr0a12
#                event.dMu1MaxDr0a12 = sorted(dimuonMaxDr0a12.pair, key = lambda x: x.pt(), reverse = True)[0] 
#                event.dMu2MaxDr0a12 = sorted(dimuonMaxDr0a12.pair, key = lambda x: x.pt(), reverse = False)[0]
#
#
#            #####################################################################################
#            # TODO: Final Qualification and 'ok' to nominate the selection dimuon as HNL candidate
#            #####################################################################################
#


        return True
