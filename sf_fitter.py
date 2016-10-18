#----------------------------------------------------------------------------------
#------------Efficiency fitter and event scale factor distribution plotter---------
#----------------------------------------------------------------------------------
# This tool fits the efficiencies from the tnp tool output, and creates a TProfile
# with the event scale factors for a 4l sample
# Thomas Gaskin, October 2016
#----------------------------------------------------------------------------------

import ROOT as rt
import numpy as np
from array import array
import math

#get user input for files
Eff_location ='T_and_P/egm_tnp_analysis/results/test/HZZ/egammaEffi.txt_SF2D.root'
Sample_location = 'qqZZ.root'

Eff = rt.TFile(Eff_location)
Sample = rt.TFile(Sample_location)
output_file =rt.TFile('SF_'+Sample_location, 'recreate')

Eta = [0.0, 0.8, 1.479, 2.0, 2.5]
ptbins = np.array([7, 15, 20, 30, 40, 50, 60, 80, 110, 150, 200], dtype=float)
linecolor=2
mc_color=12
data_color=12

mc_hists=[]
data_graphs=[]

#get the MC for each Eta Bin by projecting the TH2 from the file
for i in range (1, 5):
	mc_hists.append(Eff.Get('MC_eff_Eta_'+str(Eta[i-1])+'_'+str(Eta[i])))
	data_graphs.append(Eff.Get('Data_eff_Eta_'+str(Eta[i-1])+'_'+str(Eta[i])))

for i in range (0, 4):
	mc_hists[i].SetTitle(str(Eta[i])+'<#eta<'+str(Eta[i+1]))
#we get the data efficiencies from the TGraphs which are the output from the tnp tool. These TGraphs are arrays that are not sorted in x, 
#hence we need to sort them before filling the values into a new histogram

for data_graph in data_graphs:
	data_graph.Sort()

#we now create four new histograms, and fill them with the values from the TGraph
d1 = rt.TH1F('data_1', '0<#eta<0.8', 10, ptbins)
d2 = rt.TH1F('data_2', '0.8<#eta<1.479', 10, ptbins)
d3 = rt.TH1F('data_3', '1.479<#eta<2.0', 10,  ptbins)
d4 = rt.TH1F('data_4', '2.0<#eta<2.5', 10, ptbins)
data_hists =[d1, d2, d3, d4]
for i in range(0, 10):
	for j in range(0, 4):
		data_hists[j].SetBinContent(i+1, data_graphs[j].GetY()[i])
		data_hists[j].SetBinError(i+1, data_graphs[j].GetErrorY(i))

#fit histograms with exponential model
expo_fit = rt.TF1('expo_fit','[0]*([6]-[1]*exp([2]-[3]*x))*([4]+[7]*exp([8]-[5]*x))', 0, 10)
fitparameters = np.array([0.77, 0.3, 1, 0.01, 1, 0.01, 1, 1, 1])
for i in range (0, 9):
	expo_fit.SetParameter(i, fitparameters[i])
expo_fit.SetLineColor(linecolor)

#for the error bars, we fit the endpoints of the error bars. we create new histograms that have as data points the hi/lo error point
#we fit these hists, and draw them on top of the actual histograms on the canvas
mc_err_hi_Clones=[]
mc_err_lo_Clones=[]
d_err_hi_Clones=[]
d_err_lo_Clones=[]
for mc_hist in mc_hists:
	mc_err_hi_Clones.append(mc_hist.Clone())
	mc_err_lo_Clones.append(mc_hist.Clone())
for data_hist in data_hists:
	d_err_hi_Clones.append(data_hist.Clone())
	d_err_lo_Clones.append(data_hist.Clone())

#set bin content of the cloned hists as the hi/lo error bar points
for i in range(0, 4):
	for j in range (1, 11):
		mc_err_hi_Clones[i].SetBinContent(j, mc_hists[i].GetBinContent(j)+mc_hists[i].GetBinError(j))		#mc hi error
		mc_err_lo_Clones[i].SetBinContent(j, mc_hists[i].GetBinContent(j)-mc_hists[i].GetBinError(j))		#mc lo error
		d_err_hi_Clones[i].SetBinContent(j, data_hists[i].GetBinContent(j)+data_hists[i].GetBinError(j))	#data hi error
		d_err_lo_Clones[i].SetBinContent(j, data_hists[i].GetBinContent(j)-data_hists[i].GetBinError(j))	#data lo error

error_hists=mc_err_hi_Clones + mc_err_lo_Clones + d_err_hi_Clones + d_err_lo_Clones
mc_err_hists=mc_err_hi_Clones + mc_err_lo_Clones
data_err_hists=d_err_hi_Clones + d_err_lo_Clones
hists=data_hists+mc_hists

#fit hists
for hist in hists:
	hist.Fit('expo_fit')
	hist.GetYaxis().SetTitle('efficiencies')

expo_fit.SetLineStyle(2)
for hist in error_hists:
	hist.Fit('expo_fit')

mc_fits=[]
data_fits=[]
for i in range(0, 4):
	mc_fits.append(mc_hists[i].GetFunction('expo_fit'))
	data_fits.append(data_hists[i].GetFunction('expo_fit'))

#create the pull factors histograms. Here, the pull factors are the difference between analytic function and data point.
#the errors on the pull factors are the errors from data point
mc_1_stat = rt.TH1F('mc_1_stat', '', 10, ptbins)	
mc_2_stat = rt.TH1F('mc_2_stat', '', 10, ptbins)	
mc_3_stat = rt.TH1F('mc_3_stat', '', 10, ptbins)	
mc_4_stat = rt.TH1F('mc_4_stat', '', 10, ptbins)	

d1_stat = rt.TH1F('d1_stat', '', 10, ptbins)	
d2_stat = rt.TH1F('d2_stat', '', 10, ptbins)	
d3_stat = rt.TH1F('d3_stat', '', 10, ptbins)	
d4_stat = rt.TH1F('d4_stat', '', 10, ptbins)	

mc_stats = [mc_1_stat, mc_2_stat, mc_3_stat, mc_4_stat]
data_stats = [d1_stat, d2_stat, d3_stat, d4_stat]

for i in range(0, 4):
	for j in range(0, 10):
		funcvalue = mc_fits[i].Eval((ptbins[j+1]+ptbins[j])/2)	#evaluate the function at the same point as the data point, ie in the middle of the bin
		histvalue = mc_hists[i].GetBinContent(j+1)
		stat = (histvalue-funcvalue)
		mc_stats[i].Fill(ptbins[j], stat)
		mc_stats[i].SetBinError(j+1, mc_hists[i].GetBinError(j+1))

for i in range(0, 4):
	for j in range(0, 10):
		funcvalue = data_fits[i].Eval((ptbins[j+1]+ptbins[j])/2)
		histvalue = data_hists[i].GetBinContent(j+1)
		stat = (histvalue-funcvalue)
		data_stats[i].Fill(ptbins[j], stat)
		data_stats[i].SetBinError(j+1, data_hists[i].GetBinError(j+1))

#set colours etc and draw
for mc_hist in mc_hists:
	mc_hist.SetMarkerColor(mc_color)
	mc_hist.SetLineColor(mc_color)
	mc_hist.SetMarkerSize(1.4)
	mc_hist.SetMarkerStyle(22)
for data_hist in data_hists:
	data_hist.SetMarkerColor(data_color)
	data_hist.SetLineColor(data_color)
	data_hist.SetMarkerSize(1.2)
	data_hist.SetMarkerStyle(22)
for error_hist in mc_err_hists:
	error_hist.SetLineWidth(0)
	error_hist.SetMarkerStyle(8)
	error_hist.SetMarkerColor(mc_color)
	error_hist.SetMarkerSize(1)
for error_hist in data_err_hists:
	error_hist.SetLineWidth(0)
	error_hist.SetMarkerStyle(8)
	error_hist.SetMarkerColor(data_color)
	error_hist.SetMarkerSize(1)
for hist in hists:
	hist.SetStats(False)
	hist.GetYaxis().SetTitleSize(0.06)
	hist.GetYaxis().SetTitleOffset(0.5)
	hist.GetXaxis().SetLabelSize(0)
for hist in mc_stats+data_stats:
	hist.SetStats(False)
	hist.SetMarkerStyle(8)
	hist.SetMarkerColor(9)
	hist.GetXaxis().SetTitle('p_{T} [GeV]')
	hist.GetXaxis().SetTitleSize(0.075)
	hist.GetXaxis().SetLabelSize(0.075)
	hist.GetYaxis().SetLabelSize(0.075)
 	hist.GetYaxis().SetTitle('data point - function')
	hist.GetYaxis().SetTitleSize(0.1)
	hist.GetYaxis().SetTitleOffset(0.3)

#start drawing the canvasses. There are eight pads. The following are the positions for lower left and upper right corners:
#points are measured relative to the canvas size. The bottom left point of the canvas has coordinates (0,0), the top right corner has coordinates (1,1)
x1=0.02857143
x2=0.48571429
x3=0.51428571
x4=0.97142857
y1=0.02
y2=0.17633
y3=0.17833
y4=0.475
y5=0.495
y6=0.65133
y7=0.65333
y8=0.95

line = rt.TLine(7, 0, 200, 0)
line.SetLineStyle(2)
line.SetLineColor(2)

title = rt.TText()
title.SetTextColor(1)
title.SetTextSize(0.03)
title.SetTextAlign(11)

c1 = rt.TCanvas('MC_eff', 'MC efficiencies', 700, 1000)
title.DrawText(x1, 0.960, 'MC efficiencies')
c1.Draw()
pad1 = rt.TPad('pad1', '',x1,y7,x2,y8) 
pad2 = rt.TPad('pad2', '',x3,y7,x4,y8) 
pad3 = rt.TPad('pad3', '',x1,y5,x2,y6)
pad4 = rt.TPad('pad4', '',x3,y5,x4,y6)
pad5 = rt.TPad('pad5', '',x1,y3,x2,y4)
pad6 = rt.TPad('pad6', '',x3,y3,x4,y4)
pad7 = rt.TPad('pad7', '',x1,y1,x2,y2)
pad8 = rt.TPad('pad8', '',x3,y1,x4,y2)
pads = [pad1, pad2, pad3, pad4, pad5, pad6, pad7, pad8]

j=0
k=0
c1.cd()

for i in range(0, 8):
	pads[i].Draw()
	pads[i].cd()
	if(i==0 or i==1 or i==4 or i==5):	
		pads[i].SetBottomMargin(0.01)
		pads[i].SetTopMargin(0.01)
		mc_hists[j].Draw('expo_fit')
		mc_err_hi_Clones[j].Draw('expo_fit, same')
		mc_err_lo_Clones[j].Draw('expo_fit, same')
		j=j+1
	if(i==2 or i==3 or i==6 or i==7):
		pads[i].SetBottomMargin(0.2)
		pads[i].SetTopMargin(0.001)
		pads[i].SetTickx()
		mc_stats[k].Draw()
		line.Draw()
		k=k+1
	c1.cd()
c1.Write()

c2 = rt.TCanvas('Data_eff', 'Data efficiencies', 700, 1000)
title.DrawText(x1, 0.960, 'Data efficiencies, template fit')
c2.Draw()

j=0
k=0
c2.cd()

for i in range(0, 8):
	pads[i].Draw()
	pads[i].cd()
	if(i==0 or i==1 or i==4 or i==5):	
		pads[i].SetBottomMargin(0.01)
		pads[i].SetTopMargin(0.01)
		data_hists[j].Draw('expo_fit')
		d_err_hi_Clones[j].Draw('expo_fit, same')
		d_err_lo_Clones[j].Draw('expo_fit, same')
		j=j+1
	if(i==2 or i==3 or i==6 or i==7):
		pads[i].SetBottomMargin(0.2)
		pads[i].SetTopMargin(0.001)
		pads[i].SetTickx()
		data_stats[k].Draw()
		line.Draw()
		k=k+1
	c2.cd()
c2.Write()

c1.Close()
c2.Close()

#Draw Scalefactor map
p_Scalefactor = rt.TProfile('p_Scalefactor_analytic', 'Event scale factor', 1000, 0, 1000)

#loop over all events, selecting only electrons, retrieve scale factor from analytic function, and multiply
Tree = Sample.Get('ZZTree/candTree')
print 'Creating scale factor distribution ...'
for event in Tree:
	if (Tree.Z1Flav * Tree.Z2Flav == 121*121):
		sf=1.000000
		SF=[]
		for i in range(0, 4):
			lookupPt=event.LepPt[i]
			if lookupPt>200:
				lookupPt = 199
			Pt=[]
			Pt.append(lookupPt)
			if 0<=abs(event.LepEta[i])<0.8:
				SF.append(data_fits[0].Eval(lookupPt)/mc_hists[0].GetBinContent(np.digitize(Pt, ptbins)[0]))
			if 0.8<=abs(event.LepEta[i])<1.479:
				SF.append(data_fits[1].Eval(lookupPt)/mc_hists[1].GetBinContent(np.digitize(Pt, ptbins)[0]))
			if 1.479<=abs(event.LepEta[i])<2.0:
				SF.append(data_fits[2].Eval(lookupPt)/mc_hists[2].GetBinContent(np.digitize(Pt, ptbins)[0]))
			if 2.0<=abs(event.LepEta[i])<=2.5:
				SF.append(data_fits[3].Eval(lookupPt)/mc_hists[3].GetBinContent(np.digitize(Pt, ptbins)[0]))
		for i in range(0, 3):
			sf=sf*SF[i]
		p_Scalefactor.Fill(Tree.ZZMass, sf)

p_Scalefactor.GetXaxis().SetTitle('4e mass [GeV]')
p_Scalefactor.GetYaxis().SetTitle('event scale factor')
p_Scalefactor.Write()

SF_Map = Eff.Get('EGamma_SF2D')
SF_Map.GetXaxis().SetRangeUser(0, 2.5)
SF_Map.GetYaxis().SetTitle('p_{T} [GeV]')
SF_Map.GetYaxis().SetTitleSize(0.05)
SF_Map.GetYaxis().SetTitleOffset(1)
SF_Map.Write()

Eff.Close()
Sample.Close()
print 'Done. File saved as: ./SF_'+Sample_location





















