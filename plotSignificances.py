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
	
	Graphs = {}
	Histograms = {}
	
	printLumi = "35.9"
	
	m_n_min = 150
	m_n_max = 1550
	m_b_min = 700
	m_b_max = 1600
	
		
	masses_b = []
	masses_n = []
	significances = []
		
	m_b = m_b_min
	while m_b <= m_b_max:
		if m_b < 800:
			stepsize = 25
		else:
			stepsize = 50
			
		m_n = m_n_min
		
		while m_n < m_b:
			
			#~ if not ((m_b == 775 and m_n == 750) or (m_b == 800 and m_n == 150) or (m_b == 950 and m_n == 250) or (m_b == 950 and m_n == 300) or (m_b == 950 and m_n == 500) or (m_b == 950 and m_n == 550) or (m_b == 950 and m_n == 850) or (m_b == 950 and m_n == 900)):	
				
			#~ print "Limits/T6bbllslepton_%s_%s.result.txt"%(str(m_b),str(m_n))
			limitFile = open("Significances/T6bbllslepton_%s_%s.result.txt"%(str(m_b),str(m_n)),"r")
			masses_b.append(m_b)
			masses_n.append(m_n)
			M_SBOTTOM = "m_b_"+str(m_b)
			
			for line in limitFile:
				
				if "observed significance" in line:
					#~ significance = float(re.findall("\d+\.\d+",line)[0])
					significance = float(re.findall(r"[-+]?\d*\.\d+|\d+",line)[0])
					#~ print significance
					
					if significance < -3.5:
						significance = -3.49

					significances.append(significance)
					
			m_n += stepsize		
		m_b += stepsize
		
	bin_size =12.5
	nxbins = int(min(500,(m_b_max-m_b_min)/bin_size))
	nybins = int(min(500,(m_n_max-m_n_min)/bin_size))
	
	
	Graph = TGraph2D("Graph_significance","Local significance", len(significances), array('d',masses_b), array('d',masses_n), array('d',significances))
	Graph.SetNpx(nxbins)
	Graph.SetNpy(nybins)
	
	
	
	Histogram = Graph.GetHistogram()
	Histogram.SetTitle(";m_{#tilde{b}} [GeV]; m_{#tilde{#chi}_{2}^{0}} [GeV]")
	Histogram.GetZaxis().SetRangeUser(-3.5,3.5)
	
	
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
	latexLegendHeader.SetTextSize(0.032)
	latexLegendHeader.SetTextFont(42)
	latexLegendHeader.SetNDC(True)
	
	Overlay = ROOT.TGraph(0)
	Overlay.SetPoint(0, 700, 625)
	Overlay.SetPoint(1, 1650, 1575)
	Overlay.SetPoint(2, 700, 1575)
	Overlay.SetPoint(3, 700, 625)
	Overlay.SetFillColor(0)
	
	oneLine = ROOT.TLine(400, 375, 950, 900)
	oneLine.SetLineStyle(9)
	oneLine.SetLineWidth(2)
	
	leg1 = ROOT.TLegend(0.19,0.78,0.83,0.92)
	leg1.SetBorderSize(1)
	leg1.SetLineWidth(2)
	
	leg = ROOT.TLegend(0.195,0.785,0.825,0.82)
	leg.SetTextSize(0.025)
	leg.SetBorderSize(0)
	leg.SetLineWidth(0)
	
	plotPad.SetTopMargin(0.08)
	plotPad.SetBottomMargin(0.16)
	plotPad.SetLeftMargin(0.19)
	plotPad.SetRightMargin(0.17)
	
	plotPad.DrawFrame(700,150,1600,1750,";m_{#tilde{b}} [GeV]; m_{#tilde{#chi}_{2}^{0}} [GeV]")
	
	Histogram.SetZTitle("Observed Local Significance (#sigma)")
	Histogram.GetZaxis().SetLabelSize(0.045)
	Histogram.GetZaxis().SetTitleSize(0.045)
	Histogram.GetZaxis().SetTitleOffset(1.1)
	Histogram.Draw("SAME COLZ")
	Overlay.Draw("f")
	
	
	latexLumi.DrawLatex(0.83, 0.94, "%s fb^{-1} (13 TeV)"%(printLumi))
	latexCMS.DrawLatex(0.185,0.94,"CMS")
	latexCMSExtra.DrawLatex(0.289,0.94,"Preliminary")
	plotPad.RedrawAxis()
	leg1.Draw("same")
	leg.Draw("same")
	latexLegendHeader.DrawLatex(0.2, 0.84, "#splitline{pp #rightarrow #tilde{b} #tilde{b}, #tilde{b} #rightarrow #tilde{#chi}_{2}^{0} b, #tilde{#chi}_{2}^{0} #rightarrow #tilde{l} l / Z #tilde{#chi}_{1}^{0}, #tilde{l} #rightarrow #tilde{#chi}_{1}^{0} l}{m_{#tilde{#chi_{1}}^{0}} = 100 GeV, m_{#tilde{l}} = 0.5 (m_{#tilde{#chi}_{2}^{0}} + m_{#tilde{#chi}_{1}^{0}})}")
	canv.Update()
	canv.Print("fig/Significances.pdf")
	
	
	
	Histogram.SetName("ObsSignificance")
	
	Graph.SetName("ObsSignificance")
	
	Graph.SetTitle("Observed Local Significance in #sigma")
	
	f1 = TFile("fig/T6bbllslepton_Significance.root","RECREATE")
	Histogram.Write()
	f1.Close()
	

	
					
			
	
	
	
if (__name__ == "__main__"):
    plot()
