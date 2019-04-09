import ROOT as rt


data_B_mmm    = 'root://cms-xrd-transit.cern.ch//store/user/dezhu/2_ntuples/HN3Lv2.0/mmm/data/Single_mu_2017B/HNLTreeProducer/tree.root'
data_B_mmm_FR = '/eos/user/v/vstampf/plots/DDE/friendTree_mmm_190403_17h_41m/fr_012_mmm_sfr_190403_17h_41m.root'

f_in    = rt.TFile.Open(data_B_mmm)
f_in_FR = rt.TFile(data_B_mmm_FR)

t = f_in.Get('tree')

'''
https://root.cern.ch/doc/master/classTTree.html#a011d362261b694ee7dd780bad21f030b
If the friend tree has the same name as the original tree, you can give it an alias in the context of the friendship:
tree.AddFriend("tree1 = tree","friendfile1.root")
tree.Draw("var:ft1.v1")
'''

t.AddFriend('tfr = tree', f_in_FR)

c_fr = rt.TCanvas('fr', 'fr'); c_fr.cd()
print '\n\tt.Scan( event, event, fover1minusf012 )' 
