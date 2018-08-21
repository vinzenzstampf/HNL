import ROOT as rt


f_out = rt.TFile('pileup_TTJets_total.root', 'recreate')

h_out = rt.TH1F()

for i in range(116):
   f_in = rt.TFile('pileup_TTJets_amcat_batch_%i.root'%i)
   h_in = f_in.Get('histoname')
   h_out.Add(h_in)
   f_in.Close()

f_out.cd()
h_out.Write()
f_out.Close()
   
