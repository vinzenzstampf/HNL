import os
from glob import glob
from CMGTools.HNL.samples.signal import HN3L_M_2p5_V_0p0173205080757_e_onshell, HN3L_M_2p5_V_0p00707106781187_e_onshell, HN3L_M_2p5_V_0p00707106781187_mu_onshell, HN3L_M_2p5_V_0p004472135955_e_onshell

HN3L_M_2p5_V_0p004472135955_e_onshell.files = glob('/eos/user/v/vstampf/miniaod/HN3L_M2p5_V0p004_eos_pre2017_NLO/heavyNeutrino*.root')
HN3L_M_2p5_V_0p00707106781187_e_onshell.files = glob('/eos/user/v/vstampf/miniaod/HN3L_M2p5_V0p007_eos_pre2017_NLO/heavyNeutrino*.root')
HN3L_M_2p5_V_0p00707106781187_mu_onshell.files = glob('/eos/user/v/vstampf/miniaod/HN3L_M2p5_V0p007_mos_pre2017_NLO/heavyNeutrino*.root')
HN3L_M_2p5_V_0p0173205080757_e_onshell.files = glob('/eos/user/v/vstampf/miniaod/HN3L_M2p5_V0p017_eos_pre2017_NLO/heavyNeutrino*.root')

