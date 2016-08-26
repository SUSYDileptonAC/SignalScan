#! /usr/bin/env python

import ROOT
import sys
sys.path.append('cfg/')
from frameworkStructure import pathes
sys.path.append(pathes.basePath)
from messageLogger import messageLogger as log

ROOT.gROOT.SetBatch(True)

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
	plotPad.UseCurrentStyle()
	plotPad.Draw()	
	plotPad.cd()
	
	Graphs = {}
	Histograms = {}
	
	printLumi = "12.9"
	
	m_n_min = 150
	m_n_max = 750
	m_b_min = 450
	m_b_max = 700
	
	PDFAndScaleUncert = 0.15
	
		
	masses_b = []
	masses_n = []
	significances = []
		
	m_b = m_b_min
	m_b_stepsize = 25
	while m_b <= m_b_max:			
		m_n = m_n_min
		
		while m_n < m_b:
			
			if m_n < 300:
				m_n_stepsize = 25
			else:
				m_n_stepsize = 50			
			
			#~ print "Limits/T6bbllslepton_%s_%s.result.txt"%(str(m_b),str(m_n))
			limitFile = open("Significances/T6bbllslepton_%s_%s.result.txt"%(str(m_b),str(m_n)),"r")
			masses_b.append(m_b)
			masses_n.append(m_n)
			M_SBOTTOM = "m_b_"+str(m_b)
			
			for line in limitFile:
				
				if "observed significance" in line:
					significance = float(re.findall("\d+\.\d+",line)[0])

					significances.append(significance)
					
			m_n += m_n_stepsize		
		m_b += m_bstepsize
		
	bin_size =12.5
	nxbins = int(min(500,(m_b_max-m_b_min)/bin_size))
	nybins = int(min(500,(m_n_max-m_n_min)/bin_size))
	
	masses_b.append(700)
	masses_n.append(750)
	
	significances.append(0.)
	
	
	Graph = TGraph2D("Graph_significance","Significance", len(significances), array('d',masses_b), array('d',masses_n), array('d',significances))
	Graph.SetNpx(nxbins)
	Graph.SetNpy(nybins)
	
	
	
	Histogram = Graph.GetHistogram()
	Histogram.SetTitle(";m_{#tilde{b}} [GeV]; m_{#tilde{#chi_{2}^{0}}} [GeV]")
	Histogram.GetZaxis().SetRangeUser(0.,2.)
	
	
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
	Overlay.SetPoint(0, 450, 315)
	Overlay.SetPoint(1, 700, 640)
	Overlay.SetPoint(2, 450, 640)
	Overlay.SetPoint(3, 450, 415)
	Overlay.SetFillColor(0)
	
	oneLine = ROOT.TLine(450, 425, 700, 650)
	oneLine.SetLineStyle(9)
	oneLine.SetLineWidth(2)
	
	leg1 = ROOT.TLegend(0.18,0.80,0.83,0.92)
	leg1.SetBorderSize(1)
	leg1.SetLineWidth(2)
	
	leg = ROOT.TLegend(0.185,0.805,0.825,0.82)
	leg.SetTextSize(0.025)
	leg.SetBorderSize(0)
	leg.SetLineWidth(0)
	
	plotPad.SetTopMargin(0.08)
	plotPad.SetBottomMargin(0.16)
	plotPad.SetLeftMargin(0.18)
	plotPad.SetRightMargin(0.17)
	
	Histogram.GetYaxis().SetTitleOffset(1.3)
	Histogram.SetZTitle("observed significance (#sigma)")
	Histogram.GetZaxis().SetLabelSize(0.035)
	Histogram.GetZaxis().SetTitleOffset(0.85)
	Histogram.Draw("COLZ")
	Overlay.Draw("f")
	
	
	latexLumi.DrawLatex(0.83, 0.94, "%s fb^{-1} (13 TeV)"%(printLumi))
	latexCMS.DrawLatex(0.18,0.94,"CMS")
	latexCMSExtra.DrawLatex(0.285,0.94,"Preliminary")
	plotPad.RedrawAxis()
	leg1.Draw("same")
	leg.Draw("same")
	latexLegendHeader.DrawLatex(0.2, 0.85, "#splitline{pp#rightarrow#tilde{b}#tilde{b}, #tilde{b}#rightarrow#tilde{#chi}_{2}^{0}b, #tilde{#chi}_{2}^{0}#rightarrow#tilde{l}l/Z#tilde{#chi}_{1}^{0}, #tilde{l}#rightarrow#tilde{#chi}_{1}^{0}l; m_{#tilde{#chi_{1}}^{0}}= 100 GeV}{m_{#tilde{l}} = 0.5(m_{#tilde{#chi}_{2}^{0}}+ m_{#tilde{#chi}_{1}^{0}})}")
	canv.Update()
	canv.Print("fig/Significances.pdf")
	
	
	
	Histogram.SetName("ObsSignificance")
	
	Graph.SetName("ObsSignificance")
	
	Graph.SetTitle("Observed significance in #sigma")
	
	f1 = TFile("fig/T6bbllslepton_Significance.root","RECREATE")
	Histogram.Write()
	f1.Close()
	

	
					
			
	
	
	
if (__name__ == "__main__"):
    plot()
