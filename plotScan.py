#! /usr/bin/env python

import ROOT
import sys
sys.path.append('cfg/')
from frameworkStructure import pathes
sys.path.append(pathes.basePath)
from messageLogger import messageLogger as log

from array import *

from optparse import OptionParser


def loadPickles(path):
	from glob import glob
	import pickle
	result = {}
	for pklPath in glob(path):
		pklFile = open(pklPath, "r")
		result.update(pickle.load(pklFile))
	return result
	
			   
def plot():
	import ratios
	from ROOT import TCanvas, TPad, TH1F, TH2F, TH1I, THStack, TLegend, TMath, TF1, TFile, TGraph2D
	from setTDRStyle import setTDRStyle
	import pickle
	from defs import sbottom_masses
	from math import sqrt
	
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
	latexCMS.SetTextSize(0.055)
	latexCMS.SetNDC(True)
	latexCMSExtra = ROOT.TLatex()
	latexCMSExtra.SetTextFont(52)
	latexCMSExtra.SetTextSize(0.03)
	latexCMSExtra.SetNDC(True)
	
	canv = TCanvas("canv", "canv",1024,768)
	plotPad = ROOT.TPad("plotPad","plotPad",0,0,1,1)
	setTDRStyle()
	style=setTDRStyle()	
	style.SetPadRightMargin(0.175)	
	plotPad.UseCurrentStyle()
	plotPad.Draw()	
	plotPad.cd()	
	
	nEvents = -1
	lumi = 12900.
	printLumi = "12.9"
	
	m_neutr_1_fix = False
	#~ m_neutr_1_fix = True

	
	
	path = "shelves"
	generalSignalLabel = "T6bbllslepton"
	
	observables = ["mll"]
	leptonCombinations = ["SF-OF"]
	
	regions = ["EdgeLegacy","LowMassLowNll","LowMassHighNll","HighMassLowNll","HighMassHighNll",]
	
	
	Graphs = {}
	Histograms = {}			
	uncertaintyArrays = {}
	#~ uncertaintySources = ["Yield","StatUncertainty","SystUncertainty","TotalUncertainty","Efficiency","ISRUncertainty","pileupUncertainty","JESUncertainty","LeptonUncertainty","LeptonFastSimUncertainty","LeptonFullSimUncertainty","BTagUncertaintyLight","BTagUncertaintyHeavy"]
	uncertaintySources = ["Yield","StatUncertainty","SystUncertainty","TotalUncertainty","Efficiency","ISRUncertainty","pileupUncertainty","JESUncertainty","LeptonFullSimUncertainty","LeptonFastSimUncertainty","BTagUncertaintyLight","BTagUncertaintyHeavy","MetUncertainty"]

	for region in regions:
		for uncertaintySource in uncertaintySources:
			uncertaintyArrays["%s_%s"%(uncertaintySource,region)] = []
			
	
	
	title = "Simplified Model Scan; m(#tilde{b}) [GeV]; m(#tilde{#chi}_{2}^{0}) [GeV]"

	masses_b = []
	masses_n = []
	
	m_n_min = 150
	m_n_max = 900
	m_b_min = 500
	m_b_max = 950
	
	bin_size =25
	nxbins = int(min(500,(m_b_max-m_b_min)/bin_size))
	nybins = int(min(500,(m_n_max-m_n_min)/bin_size))
	
	
	TriggerEffUncertainty = 0.05
	LumiUncertainty = 0.062
	ElectronUncertainty = 0.05
	MuonUncertainty = 0.03	
					
	m_b = m_b_min
	while m_b <= m_b_max:
		print m_b
		if m_b < 800:
			stepsize = 25
		else:
			stepsize = 50
		
		M_SBOTTOM = "m_b_"+str(m_b)
		m_sbottom = str(m_b)
		xsection = getattr(sbottom_masses, M_SBOTTOM).cross_section13TeV
		
		m_n = m_n_min
		
		while m_n < m_b:
			
			if not ((m_b == 775 and m_n == 750) or (m_b == 800 and m_n == 150) or (m_b == 950 and m_n == 900) or (m_b == 950 and m_n == 850) or (m_b == 950 and m_n == 550) or (m_b == 950 and m_n == 500) or (m_b == 950 and m_n == 300) or (m_b == 950 and m_n == 250)):
				#~ print m_n
				m_neutralino_2 = str(m_n)
				
				masses_b.append(m_b)
				masses_n.append(m_n)
			
				
	
				Pickles = {}
				Yields = {}
				MCEvents = {}
				JES = {}
				#~ Lepton = {}
				LeptonFastSim = {}
				LeptonFullSim = {}
				Pileup = {}
				ISR= {}
				BTag = {}
				Met = {}
				
				Uncertainties = {}		
				
				
				for region in regions:
	
					Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)] = loadPickles("%s/%s_msbottom_%s_mneutralino_%s_%s_EE.pkl"%(path,generalSignalLabel,m_sbottom,m_neutralino_2,region))
					Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)] = loadPickles("%s/%s_msbottom_%s_mneutralino_%s_%s_EMu.pkl"%(path,generalSignalLabel,m_sbottom,m_neutralino_2,region))
					Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)] = loadPickles("%s/%s_msbottom_%s_mneutralino_%s_%s_MuMu.pkl"%(path,generalSignalLabel,m_sbottom,m_neutralino_2,region))
					
					Yields["EE_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EEval"]
					Yields["EMu_%s"%region] = Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuval"]
					Yields["MuMu_%s"%region] =  Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuval"]
					Yields["SFOF_%s"%region] = Yields["EE_%s"%region] + Yields["MuMu_%s"%region] - Yields["EMu_%s"%region]
					Uncertainties["Yield_%s"%region] = max(Yields["SFOF_%s"%region],0)
					
					MCEvents["EE_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EEMCEvents"]
					MCEvents["EMu_%s"%region] = Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuMCEvents"]
					MCEvents["MuMu_%s"%region] =  Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuMCEvents"]
					MCEvents["SFOF_%s"%region] = MCEvents["EE_%s"%region] + MCEvents["MuMu_%s"%region] - MCEvents["EMu_%s"%region]
					MCEvents["SFOF_%s"%region] = max(MCEvents["SFOF_%s"%region],0)
					
					Uncertainties["Efficiency_%s"%region] = Yields["SFOF_%s"%region]/(xsection*lumi)
					
					if MCEvents["SFOF_%s"%region] > 100:
						Uncertainties["StatUncertainty_%s"%region] = sqrt(MCEvents["EE_%s"%region]+MCEvents["EMu_%s"%region]+MCEvents["MuMu_%s"%region])/MCEvents["SFOF_%s"%region]
					else:
						Uncertainties["StatUncertainty_%s"%region] = 0
					
					### JES Uncertainty
					JES["Mean_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EEJESMean"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuJESMean"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuJESMean"] 
					
					JES["JESUp_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EEJESUp"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuJESUp"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuJESUp"] 
					JES["JESDown_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EEJESDown"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuJESDown"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuJESDown"] 
					
					if JES["Mean_%s"%region] > 0 and MCEvents["SFOF_%s"%region] > 100:
						Uncertainties["JESUncertainty_%s"%region] = max(abs(JES["JESUp_%s"%region]-JES["Mean_%s"%region])/JES["Mean_%s"%region],abs(JES["JESDown_%s"%region]-JES["Mean_%s"%region])/JES["Mean_%s"%region])
					else:
						Uncertainties["JESUncertainty_%s"%region] = 0
					
					
					
					### Lepton FastSim Uncertainty
					
					LeptonFastSim["Mean_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EELeptonFastSimMean"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuLeptonFastSimMean"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuLeptonFastSimMean"] 
					LeptonFastSim["LeptonNoFastSimScaleFactors_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EELeptonNoFastSimScaleFactor"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuLeptonNoFastSimScaleFactor"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuLeptonNoFastSimScaleFactor"] 
					
					if LeptonFastSim["Mean_%s"%region] > 0 and MCEvents["SFOF_%s"%region] > 100:
						Uncertainties["LeptonFastSimUncertainty_%s"%region] = abs(LeptonFastSim["LeptonNoFastSimScaleFactors_%s"%region]-LeptonFastSim["Mean_%s"%region])/LeptonFastSim["Mean_%s"%region]
					else:
						Uncertainties["LeptonFastSimUncertainty_%s"%region] = 0
					
					### Lepton FullSim Uncertainty
					
					LeptonFullSim["Mean_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EELeptonFullSimMean"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuLeptonFullSimMean"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuLeptonFullSimMean"] 
					LeptonFullSim["LeptonNoFullSimScaleFactors_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EELeptonNoFullSimScaleFactor"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuLeptonNoFullSimScaleFactor"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuLeptonNoFullSimScaleFactor"] 
					
					if LeptonFullSim["Mean_%s"%region] > 0 and MCEvents["SFOF_%s"%region] > 100:
						Uncertainties["LeptonFullSimUncertainty_%s"%region] = abs(LeptonFullSim["LeptonNoFullSimScaleFactors_%s"%region]-LeptonFullSim["Mean_%s"%region])/LeptonFullSim["Mean_%s"%region]
						Uncertainties["LeptonFullSimUncertainty_%s"%region] = 0.07
					else:
						Uncertainties["LeptonFullSimUncertainty_%s"%region] = 0
						
					
					### Lepton Uncertainty
					#~ Lepton["LeptonMean_%s"%region] = Yields["SFOF_%s"%region]
					#~ Lepton["LeptonUp_%s"%region] = Yields["EE_%s"%region] * (1+ElectronUncertainty)**2 + Yields["MuMu_%s"%region] * (1+MuonUncertainty)**2 - Yields["EMu_%s"%region] * (1+ElectronUncertainty) * (1+MuonUncertainty)
					#~ Lepton["LeptonDown_%s"%region] = Yields["EE_%s"%region] * (1-ElectronUncertainty)**2 + Yields["MuMu_%s"%region] * (1-MuonUncertainty)**2 - Yields["EMu_%s"%region] * (1-ElectronUncertainty) * (1-MuonUncertainty)
					#~ 
					#~ if Lepton["LeptonMean_%s"%region] > 0:
						#~ Uncertainties["LeptonUncertainty_%s"%region] = max(abs(Lepton["LeptonUp_%s"%region]-Lepton["LeptonMean_%s"%region])/Lepton["LeptonMean_%s"%region],abs(Lepton["LeptonDown_%s"%region] -Lepton["LeptonMean_%s"%region])/Lepton["LeptonMean_%s"%region] )
					#~ else:
						#~ Uncertainties["LeptonUncertainty_%s"%region] = 0
					#~ 
					#~ Lepton["Mean_%s"%region] = Yields["SFOF_%s"%region]
					#~ Lepton["ElectronUncertainty_EE_%s"%region] = Yields["EE_%s"%region] * ElectronUncertainty
					#~ Lepton["ElectronUncertainty_EMu_%s"%region] = Yields["EMu_%s"%region] * ElectronUncertainty
					#~ Lepton["MuonUncertainty_EMu_%s"%region] = Yields["EMu_%s"%region] * MuonUncertainty
					#~ Lepton["MuonUncertainty_MuMu_%s"%region] = Yields["MuMu_%s"%region] * MuonUncertainty
					#~ 
					#~ if Lepton["Mean_%s"%region] > 0:
						#~ Uncertainties["LeptonUncertainty_%s"%region] = sqrt( 4 * Lepton["ElectronUncertainty_EE_%s"%region]**2 + 4 * Lepton["MuonUncertainty_MuMu_%s"%region]**2
																			#~ + Lepton["ElectronUncertainty_EMu_%s"%region]**2 + Lepton["MuonUncertainty_EMu_%s"%region]**2
																			#~ + 2 * Lepton["ElectronUncertainty_EMu_%s"%region] * Lepton["MuonUncertainty_EMu_%s"%region] 
																			#~ + 8 * Lepton["ElectronUncertainty_EE_%s"%region] * Lepton["MuonUncertainty_MuMu_%s"%region]
																			#~ - 4 * Lepton["ElectronUncertainty_EE_%s"%region] * sqrt(Lepton["ElectronUncertainty_EMu_%s"%region]**2 + Lepton["MuonUncertainty_EMu_%s"%region]**2 + 2* Lepton["MuonUncertainty_EMu_%s"%region] * Lepton["ElectronUncertainty_EMu_%s"%region])
																			#~ - 4 * Lepton["MuonUncertainty_MuMu_%s"%region] * sqrt(Lepton["ElectronUncertainty_EMu_%s"%region]**2 + Lepton["MuonUncertainty_EMu_%s"%region]**2 + 2* Lepton["MuonUncertainty_EMu_%s"%region] * Lepton["ElectronUncertainty_EMu_%s"%region])
																			#~ ) / Lepton["Mean_%s"%region]
																			#~ 
					#~ else:
						#~ Uncertainties["LeptonUncertainty_%s"%region] = 0
					
					###  Pileup Uncertainty
					Pileup["Mean_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EEPileupMean"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuPileupMean"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuPileupMean"] 
					
					Pileup["PileupUp_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EEPileupUp"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuPileupUp"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuPileupUp"] 
					Pileup["PileupDown_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EEPileupDown"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuPileupDown"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuPileupDown"] 
					
					if Pileup["Mean_%s"%region] > 0 and MCEvents["SFOF_%s"%region] > 100:
						Uncertainties["pileupUncertainty_%s"%region] = max(abs(Pileup["PileupUp_%s"%region]-Pileup["Mean_%s"%region])/Pileup["Mean_%s"%region],abs(Pileup["PileupDown_%s"%region] -Pileup["Mean_%s"%region])/Pileup["Mean_%s"%region] )
					else:
						Uncertainties["pileupUncertainty_%s"%region] = 0
					
					### ISR Uncertainty
					ISR["Mean_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EEISRMean"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuISRMean"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuISRMean"] 
					
					ISR["ISRUp_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EEISRUp"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuISRUp"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuISRUp"]  
					ISR["ISRDown_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EEISRDown"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuISRDown"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuISRDown"]  
					
					if ISR["Mean_%s"%region] > 0 and MCEvents["SFOF_%s"%region] > 100:
						Uncertainties["ISRUncertainty_%s"%region] = max(abs(ISR["ISRUp_%s"%region]-ISR["Mean_%s"%region])/ISR["Mean_%s"%region],abs(ISR["ISRDown_%s"%region]-ISR["Mean_%s"%region])/ISR["Mean_%s"%region])
					else:
						Uncertainties["ISRUncertainty_%s"%region] = 0
						
					if Yields["SFOF_%s"%region] > 0 and MCEvents["SFOF_%s"%region] > 100:
						Uncertainties["TriggerEffUncertainty_%s"%region] = sqrt((Yields["EE_%s"%region]*TriggerEffUncertainty)**2 + (Yields["MuMu_%s"%region]*TriggerEffUncertainty)**2 + (Yields["EMu_%s"%region]*TriggerEffUncertainty)**2)/Yields["SFOF_%s"%region]
					else:
						Uncertainties["TriggerEffUncertainty_%s"%region] = 0
							
					### BTag Uncertainty
					BTag["Mean_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EEbTagMean"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMubTagMean"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMubTagMean"] 
					
					BTag["BTagHeavy_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EEbTagHeavy"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMubTagHeavy"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMubTagHeavy"]  
					
					BTag["BTagLight_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EEbTagLight"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMubTagLight"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMubTagLight"]  
					
					if BTag["Mean_%s"%region] > 0 and MCEvents["SFOF_%s"%region] >100:
						Uncertainties["BTagUncertaintyHeavy_%s"%region] = abs(BTag["BTagHeavy_%s"%region]-BTag["Mean_%s"%region])/BTag["Mean_%s"%region]
						Uncertainties["BTagUncertaintyLight_%s"%region] = abs(BTag["BTagLight_%s"%region]-BTag["Mean_%s"%region])/BTag["Mean_%s"%region]
					else:
						Uncertainties["BTagUncertaintyHeavy_%s"%region] = 0
						Uncertainties["BTagUncertaintyLight_%s"%region] = 0
						
					### FastSim Met Uncertainty					
					Met["Met_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EEMet"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuMet"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuMet"]  
					Met["GenMet_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EEGenMet"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuGenMet"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuGenMet"]  
					
					if Yields["SFOF_%s"%region] > 0 and MCEvents["SFOF_%s"%region] >100:
						Uncertainties["MetUncertainty_%s"%region] = 0.5*abs(Met["Met_%s"%region]-Met["GenMet_%s"%region])/Yields["SFOF_%s"%region]
					else:
						Uncertainties["MetUncertainty_%s"%region] = 0
					
					#~ ### Trigger efficiency	
					#~ if Yields["SFOF_%s"%region] > 0:
						#~ Uncertainties["TriggerEffUncertainty_%s"%region] = sqrt((Yields["EE_%s"%region]*TriggerEffUncertainty)**2 + (Yields["MuMu_%s"%region]*TriggerEffUncertainty)**2 + (Yields["EMu_%s"%region]*TriggerEffUncertainty)**2)/Yields["SFOF_%s"%region]
					#~ else:
						#~ Uncertainties["TriggerEffUncertainty_%s"%region] = 0
					Uncertainties["TriggerEffUncertainty_%s"%region] = TriggerEffUncertainty
							
					### Total syst uncertainty
					#~ Uncertainties["SystUncertainty_%s"%region] = sqrt(Uncertainties["JESUncertainty_%s"%region]**2 + Uncertainties["LeptonFastSimUncertainty_%s"%region]**2 + Uncertainties["LeptonFullSimUncertainty_%s"%region]**2 + Uncertainties["pileupUncertainty_%s"%region]**2 + Uncertainties["ISRUncertainty_%s"%region]**2 + Uncertainties["BTagUncertaintyLight_%s"%region]**2 + Uncertainties["BTagUncertaintyHeavy_%s"%region]**2  + Uncertainties["TriggerEffUncertainty_%s"%region]**2 + LumiUncertainty**2)
					Uncertainties["SystUncertainty_%s"%region] = sqrt(Uncertainties["JESUncertainty_%s"%region]**2 + Uncertainties["LeptonFullSimUncertainty_%s"%region]**2 + Uncertainties["pileupUncertainty_%s"%region]**2 + Uncertainties["ISRUncertainty_%s"%region]**2 + Uncertainties["BTagUncertaintyLight_%s"%region]**2 + Uncertainties["BTagUncertaintyHeavy_%s"%region]**2  + Uncertainties["MetUncertainty_%s"%region]**2 + Uncertainties["TriggerEffUncertainty_%s"%region]**2 + LumiUncertainty**2)
	
					### Totalaluncertainty
					Uncertainties["TotalUncertainty_%s"%region] = sqrt(Uncertainties["SystUncertainty_%s"%region]**2 + Uncertainties["StatUncertainty_%s"%region]**2)
			
					for uncertainty in uncertaintySources:
						uncertaintyArrays["%s_%s"%(uncertainty,region)].append(Uncertainties["%s_%s"%(uncertainty,region)])
						#~ Histograms["%s_%s"%(uncertainty,region)].SetBinContent(Histograms["%s_%s"%(uncertainty,region)].GetXaxis().FindBin(m_b),Histograms["%s_%s"%(uncertainty,region)].GetYaxis().FindBin(m_n),Uncertainties["%s_%s"%(uncertainty,region)])
					
										
								
							
			m_n += stepsize		
		m_b += stepsize
				
	for region in regions:
		
		for uncertainty in uncertaintySources:
			Graphs["%s_%s"%(uncertainty,region)]=TGraph2D("Graph_%s_%s"%(uncertainty,region),"%s_%s"%(uncertainty,region), len(uncertaintyArrays["%s_%s"%(uncertainty,region)]), array('d',masses_b), array('d',masses_n), array('d',uncertaintyArrays["%s_%s"%(uncertainty,region)]))
			Graphs["%s_%s"%(uncertainty,region)].SetNpx(nxbins)
			Graphs["%s_%s"%(uncertainty,region)].SetNpy(nybins)
			Histograms["%s_%s"%(uncertainty,region)] = Graphs["%s_%s"%(uncertainty,region)].GetHistogram()
			Histograms["%s_%s"%(uncertainty,region)].SetTitle(";m_{#tilde{b}} [GeV]; m_{#tilde{#chi_{2}^{0}}} [GeV]")
		
		if region == "EdgeLegacy":
			region_label = "8 TeV legacy signal region"
		else:
			region_label = "Signal Region"
		
		region_label_2 = ""
		
		if "LowMass" in region:
			region_label_2 += "low mass"	
		elif "HighMass" in region:
			region_label_2 = "high mass"	
			
		if "LowNll" in region:
			region_label_2 += ", ttbar like"
		elif "HighNll" in region:
			region_label_2 += ", non-ttbar like"		


		plotPad.SetLogz()
		Histograms["Yield_%s"%region].SetZTitle("SF-OF yield")
		Histograms["Yield_%s"%region].Draw("colz")
		latexLumi.DrawLatex(0.85, 0.96, "%s fb^{-1} (13 TeV)"%(printLumi,))
		latexCMS.DrawLatex(0.19,0.89,"CMS")
		#~ latexCMSExtra.DrawLatex(0.19,0.85,"Private Work - Simulation")
		latexCMSExtra.DrawLatex(0.19,0.85,"Simulation")
		#~ latexCMSExtra.DrawLatex(0.19,0.85,"Unpublished")
		latex.DrawLatex(0.175, 0.75, "#splitline{Simplified Model}{#splitline{T6bbslepton, m(#tilde{#chi}_{1}^{0})=100 GeV}{#splitline{"+region_label+"}{"+region_label_2+"}}}")
		canv.Update()
		canv.Print("fig/Yields/T6bbllslepton_%s.pdf"%(region))
		
		plotPad.SetLogz(0)
		Histograms["StatUncertainty_%s"%region].SetZTitle("SF-OF rel. stat. uncertainty")
		Histograms["StatUncertainty_%s"%region].Draw("colz")
		latexLumi.DrawLatex(0.85, 0.96, "%s fb^{-1} (13 TeV)"%(printLumi,))
		latexCMS.DrawLatex(0.19,0.89,"CMS")
		#~ latexCMSExtra.DrawLatex(0.19,0.85,"Private Work - Simulation")
		latexCMSExtra.DrawLatex(0.19,0.85,"Simulation")
		#~ latexCMSExtra.DrawLatex(0.19,0.85,"Unpublished")
		latex.DrawLatex(0.175, 0.75, "#splitline{Simplified Model}{#splitline{T6bbslepton, m(#tilde{#chi}_{1}^{0})=100 GeV}{#splitline{"+region_label+"}{"+region_label_2+"}}}")
		canv.Update()
		canv.Print("fig/statUncertainties/T6bbllslepton_%s_stat_err.pdf"%(region))
		
		
		Histograms["SystUncertainty_%s"%region].SetZTitle("SF-OF rel. syst. uncertainty")
		Histograms["SystUncertainty_%s"%region].Draw("colz")
		latexLumi.DrawLatex(0.85, 0.96, "%s fb^{-1} (13 TeV)"%(printLumi,))
		latexCMS.DrawLatex(0.19,0.89,"CMS")
		#~ latexCMSExtra.DrawLatex(0.19,0.85,"Private Work - Simulation")
		latexCMSExtra.DrawLatex(0.19,0.85,"Simulation")
		#~ latexCMSExtra.DrawLatex(0.19,0.85,"Unpublished")
		latex.DrawLatex(0.175, 0.75, "#splitline{Simplified Model}{#splitline{T6bbslepton, m(#tilde{#chi}_{1}^{0})=100 GeV}{#splitline{"+region_label+"}{"+region_label_2+"}}}")
		canv.Update()
		canv.Print("fig/sysUncertainties/T6bbllslepton_%s_syst_err.pdf"%(region))
		
		
		Histograms["TotalUncertainty_%s"%region].SetZTitle("SF-OF total rel. uncertainty")
		Histograms["TotalUncertainty_%s"%region].Draw("colz")
		latexLumi.DrawLatex(0.85, 0.96, "%s fb^{-1} (13 TeV)"%(printLumi,))
		latexCMS.DrawLatex(0.19,0.89,"CMS")
		#~ latexCMSExtra.DrawLatex(0.19,0.85,"Private Work - Simulation")
		latexCMSExtra.DrawLatex(0.19,0.85,"Simulation")
		#~ latexCMSExtra.DrawLatex(0.19,0.85,"Unpublished")
		latex.DrawLatex(0.175, 0.75, "#splitline{Simplified Model}{#splitline{T6bbslepton, m(#tilde{#chi}_{1}^{0})=100 GeV}{#splitline{"+region_label+"}{"+region_label_2+"}}}")
		canv.Update()
		canv.Print("fig/totUncertainties/T6bbllslepton_%s_tot_err.pdf"%(region))
		
		Histograms["Efficiency_%s"%region].SetZTitle("acceptance #times efficiency")
		Histograms["Efficiency_%s"%region].Draw("colz")
		latexLumi.DrawLatex(0.85, 0.96, "%s fb^{-1} (13 TeV)"%(printLumi,))
		latexCMS.DrawLatex(0.19,0.89,"CMS")
		#~ latexCMSExtra.DrawLatex(0.19,0.85,"Private Work - Simulation")
		latexCMSExtra.DrawLatex(0.19,0.85,"Simulation")
		#~ latexCMSExtra.DrawLatex(0.19,0.85,"Unpublished")
		latex.DrawLatex(0.175, 0.75, "#splitline{Simplified Model}{#splitline{T6bbslepton, m(#tilde{#chi}_{1}^{0})=100 GeV}{#splitline{"+region_label+"}{"+region_label_2+"}}}")
		canv.Update()
		canv.Print("fig/Efficiencies/T6bbllslepton_%s_signalEfficiency.pdf"%(region))
		
		f1 = TFile("fig/T6bbllslepton_%s_Yields_Eff_Uncertainties.root"%(region),"RECREATE")
		Histograms["Yield_%s"%region].Write()
		Histograms["Efficiency_%s"%region].Write()
		Histograms["StatUncertainty_%s"%region].Write()
		Histograms["SystUncertainty_%s"%region].Write()
		Histograms["TotalUncertainty_%s"%region].Write()
		f1.Close()
		
		
			
		for uncertainty in uncertaintySources:
			if not ( uncertainty == "Yield" or uncertainty == "StatUncertainty" or uncertainty == "SystUncertainty" or uncertainty == "TotalUncertainty" or uncertainty == "Efficiency"):			
				Histograms["%s_%s"%(uncertainty,region)].SetZTitle("%s"%uncertainty)
				Histograms["%s_%s"%(uncertainty,region)].Draw("colz")
				latexLumi.DrawLatex(0.85, 0.96, "%s fb^{-1} (13 TeV)"%(printLumi,))
				latexCMS.DrawLatex(0.19,0.89,"CMS")
				#~ latexCMSExtra.DrawLatex(0.19,0.85,"Private Work - Simulation")
				latexCMSExtra.DrawLatex(0.19,0.85,"Simulation")
				#~ latexCMSExtra.DrawLatex(0.19,0.85,"Unpublished")
				latex.DrawLatex(0.175, 0.75, "#splitline{Simplified Model}{#splitline{T6bbslepton, m(#tilde{#chi}_{1}^{0})=100 GeV}{#splitline{"+region_label+"}{"+region_label_2+"}}}")
				canv.Update()
				canv.Print("fig/%s/T6bbllslepton_%s_%s.pdf"%(uncertainty,region,uncertainty))
				#~ canv.Print("fig/%s/T6bbllslepton_%s_%s.root"%(uncertainty,region,uncertainty))
													
					
			
											
									
				

# This method just waits for a button to be pressed
def waitForInput():
    raw_input("Press any key to continue!")
    return

# entry point
#-------------
if (__name__ == "__main__"):
    # use option parser to allow verbose mode
    parser = OptionParser()
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=False,
                                  help="talk about everything")
    #~ parser.add_option("-b", "--base_path", dest="base_path", default="/user/schomakers/AnalysisData/PAT/Histos/sw532v0474/cutsV22DileptonTriggerSignal/TTJets/processedTrees_Simulation",
                                  #~ help="path to the directory containing simulated events")
    #~ parser.add_option("-n", "--nEvents", dest="nEvents", default="-1",
                                  #~ help="number of events to read (default = -1 = all). use smaller numbers for tests")
    (opts, args) = parser.parse_args()
    if (opts.verbose):
        # print out all output
        log.outputLevel = 5
    else:
        # ignore output with "debug" level
        log.outputLevel = 4

    # start
    #~ plot(opts.base_path, int(opts.nEvents), opts.observable)
    plot()
 
