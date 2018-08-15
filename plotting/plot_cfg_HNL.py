import copy
from collections import namedtuple
from operator import itemgetter

from numpy import array

from CMGTools.HNL.plotter.PlotConfigs import HistogramCfg, VariableCfg
from CMGTools.HNL.plotter.categories_TauMu import cat_Inc
from CMGTools.HNL.plotter.HistCreator import createHistograms, createTrees
from CMGTools.HNL.plotter.HistDrawer import HistDrawer
from CMGTools.HNL.plotter.Variables import taumu_vars, hnl_vars, getVars
from CMGTools.HNL.plotter.helper_methods import getVertexWeight
from CMGTools.HNL.samples.samples_mc_2017 import hnl_bkg
from pdb import set_trace
# from CMGTools.HNL.plotter.qcdEstimationMSSMltau import estimateQCDWMSSM, createQCDWHistograms
from CMGTools.HNL.plotter.defaultGroups import createDefaultGroups

from CMGTools.HNL.plotter.Samples import createSampleLists
from CMGTools.HNL.plotter.metrics import ams_hists

Cut = namedtuple('Cut', ['name', 'cut'])

int_lumi = 41000.0 #### FIXME 

binning_mssm = array([0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,190.,200.,225.,250.,275.,300.,325.,350.,400.,500.,700.,900.,1100.,1300.,1500.,1700.,1900.,2100.,2300.,2500.,2700.,2900.,3100.,3300.,3500.,3700.,3900.])

binning_mssm_btag = array([0.,20.,40.,60.,80.,100.,120.,140.,160.,180.,200.,250.,300.,350.,400.,500.,700.,900.,1100.,1300.,1500.,1700.,1900.,2100.,2300.,2500.,2700.,2900.,3100.,3300.,3500.,3700.,3900.])

binning_mva = array([0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.925, 0.95, 0.975, 0.985, 0.9925, 1.001])
binning_mva2 = array([0., 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.925, 0.95, 0.975, 1.001])

def prepareCuts():
    cuts = []
    inc_cut = '&&'.join([cat_Inc])
    # inc_cut += '&& l2_decayModeFinding'

    mt_cut = 'mt<40'
    inc_cut += '&& nbj<3'

    cuts.append(Cut('inclusive', inc_cut + '&&  l1_q != l2_q'))

    return cuts, mt_cut

def createSamples(analysis_dir, total_weight, qcd_from_same_sign, w_qcd_mssm_method, r_qcd_os_ss):
    hist_dict = {}
    sample_dict = {}
#    set_trace()
    samples_mc, samples_data, samples, all_samples, sampleDict = createSampleLists(analysis_dir=analysis_dir)

    sample_dict['all_samples'] = all_samples

    if qcd_from_same_sign and not w_qcd_mssm_method:
        samples_qcdfromss = [s for s in all_samples if s.name != 'QCD']
        samples_ss = copy.deepcopy(samples_qcdfromss)

        samples_ss = [s for s in samples_ss if not s.is_signal]

        for sample in samples_ss:
            if sample.name != 'data_obs':
                # Subtract background from data
                sample.scale = -1.

        qcd = HistogramCfg(name='QCD', var=None, cfgs=samples_ss, cut=None, total_scale=r_qcd_os_ss, lumi=int_lumi, weight=total_weight)

        samples_qcdfromss.append(qcd)
        sample_dict['samples_qcdfromss'] = samples_qcdfromss

    if w_qcd_mssm_method:
        sample_dict['samples_mssm_method'] = createQCDWHistograms(samples, hist_dict, int_lumi, weight=total_weight, r_qcd_os_ss=r_qcd_os_ss)

    return sample_dict, hist_dict

def createVariables():
    # Taken from Variables.py; can get subset with e.g. getVars(['mt', 'mvis'])
    # variables = taumu_vars
    # variables = getVars(['_norm_', 'mt', 'mvis', 'l1_pt', 'l2_pt', 'l1_eta', 'l2_eta', 'n_vertices', 'n_jets', 'n_bjets'])
    variables = hnl_vars

    return variables

def makePlots(variables, cuts, total_weight, sample_dict, hist_dict, qcd_from_same_sign, w_qcd_mssm_method, mt_cut, friend_func, dc_postfix, make_plots=True, create_trees=False):
    ams_dict = {}
    sample_names = set()
    for cut in cuts:
        cfg_main = HistogramCfg(name=cut.name, var=None, cfgs=sample_dict['all_samples'], cut=cut.cut, lumi=int_lumi, weight=total_weight)
    
        cfg_main.vars = variables

        plots = createHistograms(cfg_main, verbose=False, friend_func=friend_func)
        for variable in variables:
        # for plot in plots.itervalues():
            plot = plots[variable.name]
            createDefaultGroups(plot)
            if not w_qcd_mssm_method:
                plot.Group('W', ['W', 'W1Jets', 'W2Jets', 'W3Jets', 'W4Jets'])
            plot.Group('Electroweak', ['VV', 'W'])
            # plot.Group('Single t', ['T_tWch', 'TBar_tWch', 'TToLeptons_sch', 'TToLeptons_tch'])
#            plot.Group('ZTT', ['ZTT', 'ZJ'], style=plot.Hist('ZTT').style)
            if make_plots:
                HistDrawer.draw(plot, plot_dir='/afs/cern.ch/work/v/vstampf/plots/stacked/'+cut.name)
            if variable.name in ['mvis', 'svfit_transverse_mass', 'svfit_mass', 'mva', 'mva2div1', 'mva1', 'mva2', 'l2_nc_ratio']:
                plot.WriteDataCard(filename='datacard_{mode}_{var}.root'.format(mode=mode, var=variable.name), dir='mt_' + cut.name, mode='UPDATE', postfix=dc_postfix) #mt = mu-tau
#            for signal_hist in plot.SignalHists():
#                sample_names.add(signal_hist.name)
#                ams_dict[variable.name + '__' + cut.name + '__' + signal_hist.name + '_'] = ams_hists(signal_hist.weighted, plot.BGHist().weighted)

    print '\nOptimisation results:'
    all_vals = ams_dict.items()
    for sample_name in sample_names:
        vals = [v for v in all_vals if sample_name + '_' in v[0]]
        vals.sort(key=itemgetter(1))
        for key, item in vals:
            print item, key

        print '\nBy variable'
        for variable in variables:
            name = variable.name
            print '\nResults for variable', name
            for key, item in vals:
                if key.startswith(name + '__'):
                    print item, key

if __name__ == '__main__':
        
    data2016G = False

    friend_func = None
    
    qcd_from_same_sign = True
    w_qcd_mssm_method = True
    r_qcd_os_ss = 1.17

    run_central = True
    add_ttbar_sys = False
    add_tes_sys = False

    analysis_dir = '/eos/user/v/vstampf/ntuples/bkg_mc_prompt_e/' # input

    total_weight = ''
#    total_weight = 'weight * (1. - 0.0772790*(l2_gen_match == 5 && l2_decayMode==0) - 0.138582*(l2_gen_match == 5 && l2_decayMode==1) - 0.220793*(l2_gen_match == 5 && l2_decayMode==10) )' # Tau ID eff scale factor

    if data2016G:
        total_weight = '(' + total_weight + '*' + getVertexWeight(True) + ')'

    print total_weight

    cuts, mt_cut = prepareCuts()

    variables = createVariables()

    sample_dict, hist_dict = createSamples(analysis_dir, total_weight, qcd_from_same_sign=False, w_qcd_mssm_method=False, r_qcd_os_ss=None)
#    sample_dict_tes = {'all_samples':[s for s in sample_dict['all_samples'] if s.name in tes_samples]}
    makePlots(variables, cuts, total_weight, sample_dict, hist_dict={}, qcd_from_same_sign=False, w_qcd_mssm_method=False, mt_cut=mt_cut, friend_func=lambda f: f.replace('TESUp', 'TESUpMultiMVA'), dc_postfix='_CMS_scale_t_mt_13TeVUp', make_plots=True)



