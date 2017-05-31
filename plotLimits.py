#! /usr/bin/env python

import ROOT
import sys
sys.path.append("/.automount/home/home__home4/institut_1b/schomakers/FrameWork/frameWorkBase")
from messageLogger import messageLogger as log

from array import *

def plot():
	import re
	import ratios
	from ROOT import TCanvas, TPad, TH1F, TH2F, TH1I, THStack, TLegend, TMath, TF1, TGraph, TGraph2D, TFile, TObjArray, gROOT
	from setTDRStyle import setTDRStyle
	import pickle
	from defs import sbottom_masses
	from math import sqrt
	
	canv = TCanvas("canv", "canv",800,800)
	plotPad = ROOT.TPad("plotPad","plotPad",0,0,1,1)
	style=setTDRStyle()	
	style.SetPadRightMargin(0.18)	
	style.SetTitleYOffset(1.54)	
	plotPad.UseCurrentStyle()
	plotPad.Draw()	
	plotPad.cd()
	
	exclusionContours = ["obsR","obsR_up","obsR_down","expR","expR_up","expR_down","expR_2up","expR_2down","obsXsecLimit"]
	Graphs = {}
	Histograms = {}
	
	printLumi = "35.9"
	
	m_n_min = 150
	m_n_max = 1550
	#~ m_n_max = 1450
	m_b_min = 700
	m_b_max = 1600
	#~ m_b_max = 1450
	
	
	Exclusions = {}	
	
	for exclusionContour in exclusionContours:
		Exclusions[exclusionContour] = []
		
	masses_b = []
	masses_n = []
	cross_sections = []
		
	m_b = m_b_min
	while m_b <= m_b_max:
		if m_b < 800:
			stepsize = 25
		else:
			stepsize = 50
			
		m_n = m_n_min
		
		while m_n < m_b:
			
			if not ((m_b == 775 and m_n == 750) or (m_b == 800 and m_n == 150) or (m_b == 950 and m_n == 250) or (m_b == 950 and m_n == 300) or (m_b == 950 and m_n == 500) or (m_b == 950 and m_n == 550) or (m_b == 950 and m_n == 850) or (m_b == 950 and m_n == 900)):	
				#~ print "Limits/T6bbllslepton_%s_%s.result.txt"%(str(m_b),str(m_n))
				limitFile = open("Limits/T6bbllslepton_%s_%s.result.txt"%(str(m_b),str(m_n)),"r")
				#~ limitFile = open("LimitsMT2OnZ/T6bbllslepton_%s_%s.result.txt"%(str(m_b),str(m_n)),"r")
				masses_b.append(m_b)
				masses_n.append(m_n)
				M_SBOTTOM = "m_b_"+str(m_b)
				xSection = getattr(sbottom_masses, M_SBOTTOM).cross_section13TeV
				xsectionUncert = getattr(sbottom_masses, M_SBOTTOM).cross_sectionUncertainty
				cross_sections.append(xSection)
				
				for line in limitFile:
					
					if "CLs observed asymptotic" in line:
						observed_R = float(re.findall("\d+\.\d+",line)[0])
	
						Exclusions["obsR"].append(observed_R)
						Exclusions["obsR_up"].append(observed_R*(1+xsectionUncert))
						Exclusions["obsR_down"].append(observed_R*(1-xsectionUncert))
						
						Exclusions["obsXsecLimit"].append(observed_R*xSection)
											
					if "CLs expected asymptotic" in line:
						Exclusions["expR"].append(float(re.findall("\d+\.\d+",line)[0]))
						
					if "CLs expected m1sigma asymptotic" in line:
						Exclusions["expR_down"].append(float(re.findall("\d+\.\d+",line)[0]))
						
					if "CLs expected p1sigma asymptotic" in line:
						Exclusions["expR_up"].append(float(re.findall("\d+\.\d+",line)[0]))
					
					if "CLs expected m2sigma asymptotic" in line:
						Exclusions["expR_2down"].append(float(re.findall("\d+\.\d+",line)[0]))
						
					if "CLs expected p2sigma asymptotic" in line:
						Exclusions["expR_2up"].append(float(re.findall("\d+\.\d+",line)[0]))
				
			m_n += stepsize		
		m_b += stepsize
		
	bin_size =12.5
	nxbins = int(min(500,(m_b_max-m_b_min)/bin_size))
	nybins = int(min(500,(m_n_max-m_n_min)/bin_size))
	
	for exclusionContour in exclusionContours:
		Graphs[exclusionContour] = TGraph2D("Graph_%s"%(exclusionContour),exclusionContour, len(Exclusions[exclusionContour]), array('d',masses_b), array('d',masses_n), array('d',Exclusions[exclusionContour]))
		#~ Graphs[exclusionContour].SetNpx(nxbins)
		#~ Graphs[exclusionContour].SetNpy(nybins)
	
	dots = TGraph(len(masses_b), array('d',masses_b), array('d',masses_n))
	
	contours = array('d',[1.0])
	
	Histograms["obsXsecLimit"] = Graphs["obsXsecLimit"].GetHistogram()
	Histograms["obsXsecLimit"].SetTitle("Observed cross section limit;m_{#tilde{b}} [GeV]; m_{#tilde{#chi}_{2}^{0}} [GeV]")
	Histograms["obsXsecLimit"].GetZaxis().SetRangeUser(0.0002,0.05)
		
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
	
	Histograms["expR_2up"] = Graphs["expR_2up"].GetHistogram()	
	Histograms["expR_2up"].SetContour(1,contours)
	Histograms["expR_2up"].SetLineWidth(2)
	Histograms["expR_2up"].SetLineStyle(2)
	Histograms["expR_2up"].SetLineColor(2)
	Histograms["expR_2up"].Smooth()
	
	Histograms["expR_2down"] = Graphs["expR_2down"].GetHistogram()	
	Histograms["expR_2down"].SetContour(1,contours)
	Histograms["expR_2down"].SetLineWidth(2)
	Histograms["expR_2down"].SetLineStyle(2)
	Histograms["expR_2down"].SetLineColor(2)
	Histograms["expR_2down"].Smooth()
	
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
	
	Overlay = ROOT.TGraph(0)
	Overlay.SetPoint(0, 700, 625)
	Overlay.SetPoint(1, 1650, 1575)
	Overlay.SetPoint(2, 700, 1575)
	Overlay.SetPoint(3, 700, 625)
	Overlay.SetFillColor(0)
	
	leg1 = ROOT.TLegend(0.19,0.72,0.83,0.92)
	leg1.SetBorderSize(1)
	leg1.SetLineWidth(2)
	
	leg = ROOT.TLegend(0.195,0.725,0.825,0.82)
	leg.SetTextSize(0.025)
	leg.SetBorderSize(0)
	leg.SetLineWidth(0)
	leg.AddEntry(Histograms["expR"], "Expected limit, #pm 1 (2) #sigma_{exp.}","l")
	leg.AddEntry(Histograms["obsR"], "Observed limit, #pm 1 #sigma_{theory}","l")
	
	plotPad.SetTopMargin(0.08)
	plotPad.SetBottomMargin(0.13)
	plotPad.SetLeftMargin(0.19)
	plotPad.SetRightMargin(0.17)
	
	
	
	
	plotPad.DrawFrame(700,150,1600,1750,";m_{#tilde{b}} [GeV]; m_{#tilde{#chi}_{2}^{0}} [GeV]")
	
	plotPad.SetLogz()
	#~ Histograms["obsXsecLimit"].GetYaxis().SetTitleOffset(1.5)
	Histograms["obsXsecLimit"].SetZTitle("95% CL upper limit on #sigma [pb]")
	Histograms["obsXsecLimit"].GetZaxis().SetLabelSize(0.035)
	Histograms["obsXsecLimit"].GetZaxis().SetTitleOffset(1.0)
	Histograms["obsXsecLimit"].Draw("SAME COLZ")
	Histograms["expR_2up"].Draw("SAMECONT3")
	Histograms["expR_2down"].Draw("SAMECONT3")
	Histograms["expR_up"].Draw("SAMECONT3")
	Histograms["expR_down"].Draw("SAMECONT3")
	Histograms["expR"].Draw("SAMECONT3")
	Histograms["obsR_up"].Draw("SAMECONT3")
	Histograms["obsR_down"].Draw("SAMECONT3")
	Histograms["obsR"].Draw("SAMECONT3")
	Overlay.Draw("f")
	
	latexLumi.DrawLatex(0.83, 0.94, "%s fb^{-1} (13 TeV)"%(printLumi))
	#~ latexCMS.DrawLatex(0.18,0.94,"CMS Preliminary")
	latexCMS.DrawLatex(0.185,0.94,"CMS")
	#~ latexCMSExtra.DrawLatex(0.29,0.94,"Preliminary")
	plotPad.RedrawAxis()
	leg1.Draw("same")
	leg.Draw("same")
	latexLegendHeader.DrawLatex(0.2, 0.85, "#splitline{pp#rightarrow#tilde{b}#tilde{b}, #tilde{b}#rightarrow#tilde{#chi}_{2}^{0}b, #tilde{#chi}_{2}^{0}#rightarrow#tilde{l}l/Z#tilde{#chi}_{1}^{0}, #tilde{l}#rightarrow#tilde{#chi}_{1}^{0}l; m_{#tilde{#chi_{1}}^{0}}= 100 GeV}{m_{#tilde{l}} = 0.5(m_{#tilde{#chi}_{2}^{0}}+ m_{#tilde{#chi}_{1}^{0}}); NLO+NLL exclusion}")
	canv.Update()
	canv.Print("fig/LimitPlot.pdf")
	canv.Print("fig/LimitPlot.png")
	
	
	
	Histograms["obsXsecLimit"].SetName("XSecUpperLimit")
	Histograms["expR"].SetName("ExpectedUpperLimit")
	Histograms["expR_up"].SetName("ExpectedUpperLimitUp")
	Histograms["expR_down"].SetName("ExpectedUpperLimitDown")
	Histograms["expR_2up"].SetName("ExpectedUpperLimitUp2")
	Histograms["expR_2down"].SetName("ExpectedUpperLimitDown2")
	Histograms["obsR"].SetName("ObservedUpperLimit")
	Histograms["obsR_up"].SetName("ObservedUpperLimitUp")
	Histograms["obsR_down"].SetName("ObservedUpperLimitDown")
	
	Graphs["expR"].SetName("ExpectedUpperLimit")
	Graphs["expR_up"].SetName("ExpectedUpperLimitUp")
	Graphs["expR_down"].SetName("ExpectedUpperLimitDown")
	Graphs["expR_2up"].SetName("ExpectedUpperLimitUp2")
	Graphs["expR_2down"].SetName("ExpectedUpperLimitDown2")
	Graphs["obsR"].SetName("ObservedUpperLimit")
	Graphs["obsR_up"].SetName("ObservedUpperLimitUp")
	Graphs["obsR_down"].SetName("ObservedUpperLimitDown")
	
	Histograms["expR"].SetTitle("Expected Upper Limit;m_{#tilde{b}} [GeV]; m_{#tilde{#chi}_{2}^{0}} [GeV]")
	Histograms["expR_up"].SetTitle("Expected Upper Limit + 1 #sigma;m_{#tilde{b}} [GeV]; m_{#tilde{#chi}_{2}^{0}} [GeV]")
	Histograms["expR_down"].SetTitle("Expected Upper Limit - 1 #sigma;m_{#tilde{b}} [GeV]; m_{#tilde{#chi}_{2}^{0}} [GeV]")
	Histograms["expR_2up"].SetTitle("Expected Upper Limit + 2 #sigma;m_{#tilde{b}} [GeV]; m_{#tilde{#chi}_{2}^{0}} [GeV]")
	Histograms["expR_2down"].SetTitle("Expected Upper Limit - 2 #sigma;m_{#tilde{b}} [GeV]; m_{#tilde{#chi}_{2}^{0}} [GeV]")
	Histograms["obsR"].SetTitle("Observed Upper Limit;m_{#tilde{b}} [GeV]; m_{#tilde{#chi}_{2}^{0}} [GeV]")
	Histograms["obsR_up"].SetTitle("Observed Upper Limit + 1 #sigma;m_{#tilde{b}} [GeV]; m_{#tilde{#chi}_{2}^{0}} [GeV]")
	Histograms["obsR_down"].SetTitle("Observed Upper Limit - 1 #sigma;m_{#tilde{b}} [GeV]; m_{#tilde{#chi}_{2}^{0}} [GeV]")
	
	f1 = TFile("fig/T6bbllslepton_XSecUpperLimit_and_ExclusionContours.root","RECREATE")
	Histograms["obsXsecLimit"].Write()
	Histograms["expR"].Write()
	Histograms["expR_up"].Write()
	Histograms["expR_down"].Write()
	Histograms["expR_2up"].Write()
	Histograms["expR_2down"].Write()
	Histograms["obsR"].Write()
	Histograms["obsR_up"].Write()
	Histograms["obsR_down"].Write()
	f1.Close()
	
	

	
					
			
	
	
	
if (__name__ == "__main__"):
    plot()
