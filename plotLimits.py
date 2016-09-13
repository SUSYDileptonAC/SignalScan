#! /usr/bin/env python

import ROOT
import sys
sys.path.append('cfg/')
from frameworkStructure import pathes
sys.path.append(pathes.basePath)
from messageLogger import messageLogger as log

ROOT.gROOT.SetBatch(True)

from array import *

### Tool to make limit plots
def plot():
	import re
	import ratios
	from ROOT import TCanvas, TPad, TH1F, TH2F, TH1I, THStack, TLegend, TMath, TF1, TGraph, TGraph2D, TFile, TObjArray, gROOT
	from setTDRStyle import setTDRStyle
	import pickle
	from defs import sbottom_masses
	from math import sqrt
	
	### define canvas, pad and style
	canv = TCanvas("canv", "canv",800,800)
	plotPad = ROOT.TPad("plotPad","plotPad",0,0,1,1)
	style=setTDRStyle()	
	style.SetPadRightMargin(0.18)	
	plotPad.UseCurrentStyle()
	plotPad.Draw()	
	plotPad.cd()
	
	### The different contours used: observed +/- 1 sigma, expected +/- 1 sigma and the x-section limit
	exclusionContours = ["obsR","obsR_up","obsR_down","expR","expR_up","expR_down","obsXsecLimit"]
	Graphs = {}
	Histograms = {}
	
	printLumi = "2.3"
	
	### mass range
	m_n_min = 150
	m_n_max = 650
	m_b_min = 450
	m_b_max = 700
	
	m_b_stepsize = 25
	
	### Assume 15% theoretical uncertainty. In reasonable agreement with
	### https://twiki.cern.ch./twiki/bin/view/LHCPhysics/SUSYCrossSections13TeVstopsbottom
	### where the cross sections were taken from
	PDFAndScaleUncert = 0.15
	
	### store the masses and exclusion values for each point
	Exclusions = {}	
	
	for exclusionContour in exclusionContours:
		Exclusions[exclusionContour] = []
		
	masses_b = []
	masses_n = []
	cross_sections = []
	
	### loop over mass points	
	m_b = m_b_min
	while m_b <= m_b_max:
			
		m_n = m_n_min
		
		while m_n < m_b:
			
			if m_n < 300:
				m_n_stepsize = 25
			else:
				m_n_stepsize = 50
			
			### fetch the limit files
			limitFile = open("Limits/T6bbllslepton_%s_%s.result.txt"%(str(m_b),str(m_n)),"r")
			
			### store masses and cross sections
			masses_b.append(m_b)
			masses_n.append(m_n)
			M_SBOTTOM = "m_b_"+str(m_b)
			xSection = getattr(sbottom_masses, M_SBOTTOM).cross_section13TeV
			cross_sections.append(xSection)
			
			### fetch the different results
			for line in limitFile:
				
				if "CLs observed asymptotic" in line:
					observed_R = float(re.findall("\d+\.\d+",line)[0])

					Exclusions["obsR"].append(observed_R)
					Exclusions["obsR_up"].append(observed_R*(1+PDFAndScaleUncert))
					Exclusions["obsR_down"].append(observed_R*(1-PDFAndScaleUncert))
					
					Exclusions["obsXsecLimit"].append(observed_R*xSection)
										
				if "CLs expected asymptotic" in line:
					Exclusions["expR"].append(float(re.findall("\d+\.\d+",line)[0]))
					
				if "CLs expected m1sigma asymptotic" in line:
					Exclusions["expR_down"].append(float(re.findall("\d+\.\d+",line)[0]))
					
				if "CLs expected p1sigma asymptotic" in line:
					Exclusions["expR_up"].append(float(re.findall("\d+\.\d+",line)[0]))
				
				
			m_n += m_n_stepsize		
		m_b += m_b_stepsize
	
	### binning for the plot	
	bin_size =12.5
	nxbins = int(min(500,(m_b_max-m_b_min)/bin_size))
	nybins = int(min(500,(m_n_max-m_n_min)/bin_size))
	
	### create 2D graphs from the arrays of masses and r-values
	for exclusionContour in exclusionContours:
		Graphs[exclusionContour] = TGraph2D("Graph_%s"%(exclusionContour),exclusionContour, len(Exclusions[exclusionContour]), array('d',masses_b), array('d',masses_n), array('d',Exclusions[exclusionContour]))
		Graphs[exclusionContour].SetNpx(nxbins)
		Graphs[exclusionContour].SetNpy(nybins)
	
	### set contour to be plotted. Distinguish points with r < 1 (excluded)
	### and r > 1 (not excluded)	
	contours = array('d',[1.0])
	
	### clolor coded histogram for the observed x-section limit
	Histograms["obsXsecLimit"] = Graphs["obsXsecLimit"].GetHistogram()
	Histograms["obsXsecLimit"].SetTitle(";m_{#tilde{b}} [GeV]; m_{#tilde{#chi_{2}^{0}}} [GeV]")
	Histograms["obsXsecLimit"].GetZaxis().SetRangeUser(0.1,1.)
	
	### exclusion contours	
	Histograms["expR"] = Graphs["expR"].GetHistogram()	
	Histograms["expR"].SetContour(1,contours)
	Histograms["expR"].SetLineWidth(4)
	Histograms["expR"].SetLineStyle(2)
	Histograms["expR"].SetLineColor(2)
	Histograms["expR"].Smooth()
		
	Histograms["expR_up"] = Graphs["expR_up"].GetHistogram()	
	Histograms["expR_up"].SetContour(1,contours)
	Histograms["expR_up"].SetLineWidth(2)
	Histograms["expR_up"].SetLineStyle(2)
	Histograms["expR_up"].SetLineColor(2)
	Histograms["expR_up"].Smooth()
	
	Histograms["expR_down"] = Graphs["expR_down"].GetHistogram()	
	Histograms["expR_down"].SetContour(1,contours)
	Histograms["expR_down"].SetLineWidth(2)
	Histograms["expR_down"].SetLineStyle(2)
	Histograms["expR_down"].SetLineColor(2)
	Histograms["expR_down"].Smooth()
	
	Histograms["obsR"] = Graphs["obsR"].GetHistogram()	
	Histograms["obsR"].SetContour(1,contours)
	Histograms["obsR"].SetLineWidth(4)
	Histograms["obsR"].SetLineStyle(1)
	Histograms["obsR"].SetLineColor(1)
	Histograms["obsR"].Smooth()
	
	Histograms["obsR_up"] = Graphs["obsR_up"].GetHistogram()	
	Histograms["obsR_up"].SetContour(1,contours)
	Histograms["obsR_up"].SetLineWidth(2)
	Histograms["obsR_up"].SetLineStyle(1)
	Histograms["obsR_up"].SetLineColor(1)
	Histograms["obsR_up"].Smooth()
	
	Histograms["obsR_down"] = Graphs["obsR_down"].GetHistogram()	
	Histograms["obsR_down"].SetContour(1,contours)
	Histograms["obsR_down"].SetLineWidth(2)
	Histograms["obsR_down"].SetLineStyle(1)
	Histograms["obsR_down"].SetLineColor(1)
	Histograms["obsR_down"].Smooth()
	
	
	### labels
	latex = ROOT.TLatex()
	latex.SetTextSize(0.03)
	latex.SetNDC(True)
	latexLumi = ROOT.TLatex()
	latexLumi.SetTextFont(42)
	latexLumi.SetTextAlign(31)
	latexLumi.SetTextSize(0.04)
	latexLumi.SetNDC(True)
	latexCMS = ROOT.TLatex()
	latexCMS.SetTextFont(61)
	latexCMS.SetTextSize(0.05)
	latexCMS.SetNDC(True)
	latexCMSExtra = ROOT.TLatex()
	latexCMSExtra.SetTextFont(52)
	latexCMSExtra.SetTextSize(0.0375)
	latexCMSExtra.SetNDC(True)
	
	latexLegendHeader = ROOT.TLatex() 
	latexLegendHeader.SetTextSize(0.025)
	latexLegendHeader.SetNDC(True)
	
	### Triangular overlay to make the diagonal look nicer
	Overlay = ROOT.TGraph(0)
	Overlay.SetPoint(0, 450, 415)
	Overlay.SetPoint(1, 700, 640)
	Overlay.SetPoint(2, 450, 640)
	Overlay.SetPoint(3, 450, 415)
	Overlay.SetFillColor(0)
	
	### Somehow the legend content does not stick to the defined legend size
	### Thus we are using two legends, a larger one for the frame
	### and a smaller one with the actual content
	
	leg1 = ROOT.TLegend(0.18,0.72,0.83,0.92)
	leg1.SetBorderSize(1)
	leg1.SetLineWidth(2)
	
	leg = ROOT.TLegend(0.185,0.725,0.825,0.82)
	leg.SetTextSize(0.025)
	leg.SetBorderSize(0)
	leg.SetLineWidth(0)
	leg.AddEntry(Histograms["expR"], "Expected limit, #pm 1 #sigma_{exp.}","l")
	leg.AddEntry(Histograms["obsR"], "Observed limit, #pm 1 #sigma_{theory}","l")
	
	### set margins around the plot to improve the style
	plotPad.SetTopMargin(0.08)
	plotPad.SetBottomMargin(0.16)
	plotPad.SetLeftMargin(0.18)
	plotPad.SetRightMargin(0.17)
	
	### Draw the histograms
	plotPad.SetLogz()
	Histograms["obsXsecLimit"].GetYaxis().SetTitleOffset(1.3)
	Histograms["obsXsecLimit"].SetZTitle("95% CL upper limit on #sigma [pb]")
	Histograms["obsXsecLimit"].GetZaxis().SetLabelSize(0.035)
	Histograms["obsXsecLimit"].GetZaxis().SetTitleOffset(0.85)
	Histograms["obsXsecLimit"].Draw("COLZ")
	Histograms["expR_up"].Draw("SAMECONT3")
	Histograms["expR_down"].Draw("SAMECONT3")
	Histograms["expR"].Draw("SAMECONT3")
	Histograms["obsR_up"].Draw("SAMECONT3")
	Histograms["obsR_down"].Draw("SAMECONT3")
	Histograms["obsR"].Draw("SAMECONT3")
	Overlay.Draw("f")
	
	latexLumi.DrawLatex(0.83, 0.94, "%s fb^{-1} (13 TeV)"%(printLumi))
	latexCMS.DrawLatex(0.18,0.94,"CMS")
	latexCMSExtra.DrawLatex(0.285,0.94,"Private Work")
	plotPad.RedrawAxis()
	leg1.Draw("same")
	leg.Draw("same")
	latexLegendHeader.DrawLatex(0.2, 0.85, "#splitline{pp#rightarrow#tilde{b}#tilde{b}, #tilde{b}#rightarrow#tilde{#chi}_{2}^{0}b, #tilde{#chi}_{2}^{0}#rightarrow#tilde{l}l/Z#tilde{#chi}_{1}^{0}, #tilde{l}#rightarrow#tilde{#chi}_{1}^{0}l; m_{#tilde{#chi_{1}}^{0}}= 100 GeV}{m_{#tilde{l}} = 0.5(m_{#tilde{#chi}_{2}^{0}}+ m_{#tilde{#chi}_{1}^{0}}); NLO+NLL exclusion}")
	canv.Update()
	canv.Print("fig/LimitPlot.pdf")
	
	
	### Write observed x -section histogram and graphs with the exclusion
	### limits into a root file. In this way other people can use
	### the limits more easily
	Histograms["obsXsecLimit"].SetName("XSecUpperLimit")
	
	Graphs["expR"].SetName("ExpectedUpperLimit")
	Graphs["expR_up"].SetName("ExpectedUpperLimitUp")
	Graphs["expR_down"].SetName("ExpectedUpperLimitDown")
	Graphs["obsR"].SetName("ObservedUpperLimit")
	Graphs["obsR_up"].SetName("ObservedUpperLimitUp")
	Graphs["obsR_down"].SetName("ObservedUpperLimitDown")
	
	Graphs["expR"].SetTitle("Expected Upper Limit")
	Graphs["expR_up"].SetTitle("Expected Upper Limit + 1 #sigma")
	Graphs["expR_down"].SetTitle("Expected Upper Limit - 1 #sigma")
	Graphs["obsR"].SetTitle("Observed Upper Limit")
	Graphs["obsR_up"].SetTitle("Observed Upper Limit + 1 #sigma")
	Graphs["obsR_down"].SetTitle("Observed Upper Limit - 1 #sigma")
	
	f1 = TFile("fig/T6bbllslepton_XSecUpperLimit_and_ExclusionContours.root","RECREATE")
	Histograms["obsXsecLimit"].Write()
	Graphs["expR"].Write()
	Graphs["expR_up"].Write()
	Graphs["expR_down"].Write()
	Graphs["obsR"].Write()
	Graphs["obsR_up"].Write()
	Graphs["obsR_down"].Write()
	f1.Close()
	
	

	
					
			
	
	
	
if (__name__ == "__main__"):
    plot()
