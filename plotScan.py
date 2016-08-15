#! /usr/bin/env python

import ROOT
import sys
sys.path.append('cfg/')
from frameworkStructure import pathes
sys.path.append(pathes.basePath)
from messageLogger import messageLogger as log

from array import *

from optparse import OptionParser

ROOT.gROOT.SetBatch(True)


def loadPickles(path):
	from glob import glob
	import pickle
	result = {}
	for pklPath in glob(path):
		pklFile = open(pklPath, "r")
		result.update(pickle.load(pklFile))
	return result
	
			   
def plot(systematics):
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
	lumi = 2300.
	printLumi = "2.3"
	

	
	if systematics:
		path = "shelves"
	else:
		path = "shelvesYields"
		
	generalSignalLabel = "T6bbllslepton"
	
	observables = ["mll"]
	leptonCombinations = ["SF-OF"]
	
	etaRegions = ["central","forward"]
	massRegions = ["lowMass","belowZ","onZ","aboveZ","highMass"]
	bTagRegions = ["noBTag","geOneBTag","inclusiveBTags"]
	bTagRegionsExclusive = ["noBTag","geOneBTag"]
	
	regions = []
	regionCombinations = []
	regionCombinationsInclusive = []
	
	for etaRegion in etaRegions:
		for massRegion in massRegions:
			for bTagRegion in bTagRegions:
				regions.append("%s_%s_%s"%(etaRegion,massRegion,bTagRegion))
	

	for etaRegion in etaRegions:	
		for bTagRegion in bTagRegions:
			regionCombinations.append("%s_%s"%(etaRegion,bTagRegion))	
			
	for massRegion in massRegions:
		for bTagRegion in bTagRegions:
			regionCombinations.append("%s_%s"%(massRegion,bTagRegion))	
			
	Graphs = {}
	Histograms = {}			
	uncertaintyArrays = {}		
	
	if systematics:			
		uncertaintySources = ["Yield","StatUncertainty","SystUncertainty","TotalUncertainty","Efficiency","ISRUncertainty","pileupUncertainty","JESUncertainty","LeptonFullSimUncertainty","LeptonFastSimUncertainty","BTagUncertaintyLight","BTagUncertaintyHeavy",]
	else:			
		uncertaintySources = ["Yield","StatUncertainty","Efficiency"]
	
	for uncertaintySource in uncertaintySources:	
		for region in regions:
			uncertaintyArrays["%s_%s"%(uncertaintySource,region)] = []
		for bTagRegion in bTagRegions:
			uncertaintyArrays["%s_%s"%(uncertaintySource,bTagRegion)] = []
		for regionCombination in regionCombinations:
			uncertaintyArrays["%s_%s"%(uncertaintySource,regionCombination)] = []

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
	LumiUncertainty = 0.027
	
	m_n_min = 150
	m_b_min = 400
	m_b_max = 950
	
	
	TriggerEffUncertainty = 0.05
	PDFUncertainty = 0.0
	LumiUncertainty = 0.046
	ElectronUncertainty = 0.05
	MuonUncertainty = 0.01	
					
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
			
			if not ((m_b == 775 and m_n == 750) or (m_b == 800 and m_n == 150)):
				#~ print m_n
				m_neutralino_2 = str(m_n)
				
				masses_b.append(m_b)
				masses_n.append(m_n)
			
				
	
				Pickles = {}
				Yields = {}
				MCEvents = {}
				if systematics:
					JES = {}
					LeptonFastSim = {}
					LeptonFullSim = {}
					Pileup = {}
					ISR= {}
					BTag = {}
				
				Uncertainties = {}		
				
				
				for region in regions:	
					Pickles["%s_%s"%(m_sbottom,m_neutralino_2)] = loadPickles("%s/%s_msbottom_%s_mneutralino_%s_%s.pkl"%(path,generalSignalLabel,m_sbottom,m_neutralino_2,region))
					
					Yields["EE_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["EE"]["EEval"]
					Yields["EMu_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuval"]
					Yields["MuMu_%s"%region] =  Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuval"]
					Yields["SFOF_%s"%region] = Yields["EE_%s"%region] + Yields["MuMu_%s"%region] - Yields["EMu_%s"%region]
					Uncertainties["Yield_%s"%region] = max(Yields["SFOF_%s"%region],0)
					
					MCEvents["EE_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["EE"]["EEMCEvents"]
					MCEvents["EMu_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuMCEvents"]
					MCEvents["MuMu_%s"%region] =  Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuMCEvents"]
					MCEvents["SFOF_%s"%region] = MCEvents["EE_%s"%region] + MCEvents["MuMu_%s"%region] - MCEvents["EMu_%s"%region]
					MCEvents["SFOF_%s"%region] = max(MCEvents["SFOF_%s"%region],0)
					
					if MCEvents["SFOF_%s"%region] > 0:
						Uncertainties["StatUncertainty_%s"%region] = sqrt(MCEvents["EE_%s"%region]+MCEvents["EMu_%s"%region]+MCEvents["MuMu_%s"%region])/MCEvents["SFOF_%s"%region]
					else:
						Uncertainties["StatUncertainty_%s"%region] = 0
					
					Uncertainties["Efficiency_%s"%region] = Yields["SFOF_%s"%region]/(xsection*lumi)
					
					if systematics:
					
						### JES Uncertainty
						JES["JESUp_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["EE"]["EEJESUp"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuJESUp"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuJESUp"] 
						JES["JESDown_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["EE"]["EEJESDown"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuJESDown"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuJESDown"] 			
						
						### Lepton FastSim Uncertainty					
						LeptonFastSim["NoLeptonFastSimScaleFactors_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["EE"]["EENoLeptonFastSimScaleFactor"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuNoLeptonFastSimScaleFactor"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuNoLeptonFastSimScaleFactor"] 
											
						### Lepton FullSim Uncertainty
						LeptonFullSim["NoLeptonFullSimScaleFactors_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["EE"]["EENoLeptonFullSimScaleFactor"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuNoLeptonFullSimScaleFactor"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuNoLeptonFullSimScaleFactor"] 
						
						###  Pileup Uncertainty
						Pileup["PileupUp_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["EE"]["EEPileupUp"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuPileupUp"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuPileupUp"] 
						Pileup["PileupDown_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["EE"]["EEPileupDown"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuPileupDown"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuPileupDown"] 
											
						### ISR Uncertainty
						ISR["ISRUp_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["EE"]["EEISRUp"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuISRUp"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuISRUp"]  
						ISR["ISRDown_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["EE"]["EEISRDown"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuISRDown"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuISRDown"]  
														
						### BTag Uncertainty
						BTag["BTagHeavy_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["EE"]["EEbTagHeavy"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMubTagHeavy"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["EMu"]["EMubTagHeavy"]  
						
						BTag["BTagLight_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["EE"]["EEbTagLight"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMubTagLight"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["EMu"]["EMubTagLight"]  
						
						if Yields["SFOF_%s"%region] > 0:
							Uncertainties["JESUncertainty_%s"%region] = max(abs(JES["JESUp_%s"%region]-Yields["SFOF_%s"%region])/Yields["SFOF_%s"%region],abs(JES["JESDown_%s"%region]-Yields["SFOF_%s"%region])/Yields["SFOF_%s"%region])
							Uncertainties["LeptonFastSimUncertainty_%s"%region] = abs(LeptonFastSim["NoLeptonFastSimScaleFactors_%s"%region]-Yields["SFOF_%s"%region])/Yields["SFOF_%s"%region]
							Uncertainties["LeptonFullSimUncertainty_%s"%region] = abs(LeptonFullSim["NoLeptonFullSimScaleFactors_%s"%region]-Yields["SFOF_%s"%region])/Yields["SFOF_%s"%region]
							Uncertainties["pileupUncertainty_%s"%region] = max(abs(Pileup["PileupUp_%s"%region]-Yields["SFOF_%s"%region])/Yields["SFOF_%s"%region],abs(Pileup["PileupDown_%s"%region] -Yields["SFOF_%s"%region])/Yields["SFOF_%s"%region] )
							Uncertainties["ISRUncertainty_%s"%region] = max(abs(ISR["ISRUp_%s"%region]-Yields["SFOF_%s"%region])/Yields["SFOF_%s"%region],abs(ISR["ISRDown_%s"%region]-Yields["SFOF_%s"%region])/Yields["SFOF_%s"%region])
							Uncertainties["BTagUncertaintyHeavy_%s"%region] = abs(BTag["BTagHeavy_%s"%region]-Yields["SFOF_%s"%region])/Yields["SFOF_%s"%region]
							Uncertainties["BTagUncertaintyLight_%s"%region] = abs(BTag["BTagLight_%s"%region]-Yields["SFOF_%s"%region])/Yields["SFOF_%s"%region]
						else:
							Uncertainties["StatUncertainty_%s"%region] = 0
							Uncertainties["JESUncertainty_%s"%region] = 0
							Uncertainties["LeptonFastSimUncertainty_%s"%region] = 0
							Uncertainties["LeptonFullSimUncertainty_%s"%region] = 0
							Uncertainties["pileupUncertainty_%s"%region] = 0
							Uncertainties["ISRUncertainty_%s"%region] = 0
							Uncertainties["BTagUncertaintyHeavy_%s"%region] = 0
							Uncertainties["BTagUncertaintyLight_%s"%region] = 0
								
						### Total syst uncertainty
						Uncertainties["SystUncertainty_%s"%region] = sqrt(Uncertainties["JESUncertainty_%s"%region]**2 + Uncertainties["LeptonFullSimUncertainty_%s"%region]**2 + Uncertainties["pileupUncertainty_%s"%region]**2 + Uncertainties["ISRUncertainty_%s"%region]**2 + Uncertainties["BTagUncertaintyLight_%s"%region]**2 + Uncertainties["BTagUncertaintyHeavy_%s"%region]**2  + TriggerEffUncertainty**2 + LumiUncertainty**2)
		
						### Totalaluncertainty
						Uncertainties["TotalUncertainty_%s"%region] = sqrt(Uncertainties["SystUncertainty_%s"%region]**2 + Uncertainties["StatUncertainty_%s"%region]**2)
				
					for uncertainty in uncertaintySources:
						uncertaintyArrays["%s_%s"%(uncertainty,region)].append(Uncertainties["%s_%s"%(uncertainty,region)])
						
										
							
						
				for regionCombination in regionCombinations:
					
					region1, region2 = regionCombination.split('_')
					
					Yields["EE_%s"%regionCombination] = 0	
					Yields["EMu_%s"%regionCombination] = 0	
					Yields["MuMu_%s"%regionCombination] = 0	
					Yields["SFOF_%s"%regionCombination] = 0
					
					Uncertainties["Yield_%s"%regionCombination] = 0
					
					MCEvents["EE_%s"%regionCombination] = 0
					MCEvents["EMu_%s"%regionCombination] = 0
					MCEvents["MuMu_%s"%regionCombination] = 0
					MCEvents["SFOF_%s"%regionCombination] = 0
					
					if systematics:
						
						JES["JESUp_%s"%regionCombination] =  0
						JES["JESDown_%s"%regionCombination] =  0
						
						LeptonFastSim["NoLeptonFastSimScaleFactors_%s"%regionCombination] =  0
						LeptonFullSim["NoLeptonFullSimScaleFactors_%s"%regionCombination] =  0
						
						Pileup["PileupUp_%s"%regionCombination] =  0
						Pileup["PileupDown_%s"%regionCombination] =  0
											
						ISR["ISRUp_%s"%regionCombination] = 0
						ISR["ISRDown_%s"%regionCombination] = 0
											
						BTag["BTagHeavy_%s"%regionCombination] = 0
						BTag["BTagLight_%s"%regionCombination] = 0
						
					### Add information of single regions into combinations
					for region in regions:
							
						if region1 in region and region2 in region:
								
							Yields["EE_%s"%regionCombination] += Yields["EE_%s"%region]
							Yields["EMu_%s"%regionCombination] += Yields["EMu_%s"%region]
							Yields["MuMu_%s"%regionCombination] += Yields["MuMu_%s"%region]
							Yields["SFOF_%s"%regionCombination] += Yields["SFOF_%s"%region]
							Uncertainties["Yield_%s"%regionCombination] += Uncertainties["Yield_%s"%region]
							
							MCEvents["EE_%s"%regionCombination] += MCEvents["EE_%s"%region]
							MCEvents["EMu_%s"%regionCombination] += MCEvents["EMu_%s"%region]
							MCEvents["MuMu_%s"%regionCombination] += MCEvents["MuMu_%s"%region]
							MCEvents["SFOF_%s"%regionCombination] += MCEvents["SFOF_%s"%region]
							
							if systematics:
								JES["JESUp_%s"%regionCombination] +=  JES["JESUp_%s"%region]
								JES["JESDown_%s"%regionCombination] +=  JES["JESDown_%s"%region]
								
								LeptonFastSim["NoLeptonFastSimScaleFactors_%s"%regionCombination] +=  LeptonFastSim["NoLeptonFastSimScaleFactors_%s"%region]
								
								LeptonFullSim["NoLeptonFullSimScaleFactors_%s"%regionCombination] +=  LeptonFullSim["NoLeptonFullSimScaleFactors_%s"%region]
													
								Pileup["PileupUp_%s"%regionCombination] +=  Pileup["PileupUp_%s"%region]
								Pileup["PileupDown_%s"%regionCombination] +=  Pileup["PileupDown_%s"%region]
												
								ISR["ISRUp_%s"%regionCombination] += ISR["ISRUp_%s"%region]
								ISR["ISRDown_%s"%regionCombination] += ISR["ISRDown_%s"%region]		
											
								BTag["BTagHeavy_%s"%regionCombination] += BTag["BTagHeavy_%s"%region]
								BTag["BTagLight_%s"%regionCombination] += BTag["BTagLight_%s"%region]		
					
					##Calculate uncertainties for combinations			
					Uncertainties["Efficiency_%s"%regionCombination] = Yields["SFOF_%s"%regionCombination]/(xsection*lumi)
						
					if MCEvents["SFOF_%s"%regionCombination] > 0:
						Uncertainties["StatUncertainty_%s"%regionCombination] = sqrt(MCEvents["EE_%s"%regionCombination]+MCEvents["EMu_%s"%regionCombination]+MCEvents["MuMu_%s"%regionCombination])/MCEvents["SFOF_%s"%regionCombination]
					else:
						Uncertainties["StatUncertainty_%s"%regionCombination] = 0
					
					if systematics:
						if Yields["SFOF_%s"%regionCombination] > 0:
							Uncertainties["JESUncertainty_%s"%regionCombination] = max(abs(JES["JESUp_%s"%regionCombination]-Yields["SFOF_%s"%regionCombination])/Yields["SFOF_%s"%regionCombination],abs(JES["JESDown_%s"%regionCombination]-Yields["SFOF_%s"%regionCombination])/Yields["SFOF_%s"%regionCombination])
							Uncertainties["LeptonFastSimUncertainty_%s"%regionCombination] = abs(LeptonFastSim["NoLeptonFastSimScaleFactors_%s"%regionCombination]-Yields["SFOF_%s"%regionCombination])/Yields["SFOF_%s"%regionCombination]
							Uncertainties["LeptonFullSimUncertainty_%s"%regionCombination] = abs(LeptonFullSim["NoLeptonFullSimScaleFactors_%s"%regionCombination]-Yields["SFOF_%s"%regionCombination])/Yields["SFOF_%s"%regionCombination]
							Uncertainties["pileupUncertainty_%s"%regionCombination] = max(abs(Pileup["PileupUp_%s"%regionCombination]-Yields["SFOF_%s"%regionCombination])/Yields["SFOF_%s"%regionCombination],abs(Pileup["PileupDown_%s"%regionCombination] -Yields["SFOF_%s"%regionCombination])/Yields["SFOF_%s"%regionCombination] )
							Uncertainties["ISRUncertainty_%s"%regionCombination] = max(abs(ISR["ISRUp_%s"%regionCombination]-Yields["SFOF_%s"%regionCombination])/Yields["SFOF_%s"%regionCombination],abs(ISR["ISRDown_%s"%regionCombination]-Yields["SFOF_%s"%regionCombination])/Yields["SFOF_%s"%regionCombination])
							Uncertainties["BTagUncertaintyHeavy_%s"%regionCombination] = abs(BTag["BTagHeavy_%s"%regionCombination]-Yields["SFOF_%s"%regionCombination])/Yields["SFOF_%s"%regionCombination]
							Uncertainties["BTagUncertaintyLight_%s"%regionCombination] = abs(BTag["BTagLight_%s"%regionCombination]-Yields["SFOF_%s"%regionCombination])/Yields["SFOF_%s"%regionCombination]
						else:
							Uncertainties["JESUncertainty_%s"%regionCombination] = 0
							Uncertainties["LeptonFastSimUncertainty_%s"%regionCombination] = 0
							Uncertainties["LeptonFullSimUncertainty_%s"%regionCombination] = 0
							Uncertainties["pileupUncertainty_%s"%regionCombination] = 0
							Uncertainties["ISRUncertainty_%s"%regionCombination] = 0
							Uncertainties["BTagUncertaintyHeavy_%s"%regionCombination] = 0
							Uncertainties["BTagUncertaintyLight_%s"%regionCombination] = 0
						
								
						### Total syst uncertainty
						Uncertainties["SystUncertainty_%s"%regionCombination] = sqrt(Uncertainties["JESUncertainty_%s"%regionCombination]**2 + Uncertainties["LeptonFastSimUncertainty_%s"%regionCombination]**2 + Uncertainties["LeptonFullSimUncertainty_%s"%regionCombination]**2 + Uncertainties["pileupUncertainty_%s"%regionCombination]**2 + Uncertainties["ISRUncertainty_%s"%regionCombination]**2 + Uncertainties["BTagUncertaintyLight_%s"%regionCombination]**2  + Uncertainties["BTagUncertaintyHeavy_%s"%regionCombination]**2  + TriggerEffUncertainty**2 + LumiUncertainty**2)
		
						### Total uncertainty
						Uncertainties["TotalUncertainty_%s"%regionCombination] = sqrt(Uncertainties["SystUncertainty_%s"%regionCombination]**2 + Uncertainties["StatUncertainty_%s"%regionCombination]**2)		
											
					for uncertainty in uncertaintySources:
						uncertaintyArrays["%s_%s"%(uncertainty,regionCombination)].append(Uncertainties["%s_%s"%(uncertainty,regionCombination)])
						
									
				for bTagRegion in bTagRegions:
					
	
					Yields["EE_%s"%bTagRegion] = 0	
					Yields["EMu_%s"%bTagRegion] = 0	
					Yields["MuMu_%s"%bTagRegion] = 0	
					Yields["SFOF_%s"%bTagRegion] = 0
					
					Uncertainties["Yield_%s"%bTagRegion] = 0
					
					MCEvents["EE_%s"%bTagRegion] = 0
					MCEvents["EMu_%s"%bTagRegion] = 0
					MCEvents["MuMu_%s"%bTagRegion] = 0
					MCEvents["SFOF_%s"%bTagRegion] = 0
					
					if systematics:
						JES["JESUp_%s"%bTagRegion] =  0
						JES["JESDown_%s"%bTagRegion] =  0
																
						LeptonFastSim["NoLeptonFastSimScaleFactors_%s"%bTagRegion] =  0
												
						LeptonFullSim["NoLeptonFullSimScaleFactors_%s"%bTagRegion] =  0
						
						Pileup["PileupUp_%s"%bTagRegion] =  0
						Pileup["PileupDown_%s"%bTagRegion] =  0
										
						ISR["ISRUp_%s"%bTagRegion] = 0
						ISR["ISRDown_%s"%bTagRegion] = 0
											
						BTag["BTagHeavy_%s"%bTagRegion] = 0
						BTag["BTagLight_%s"%bTagRegion] = 0
						
					### Add information of single regions into combinations
					for region in regions:						
						if bTagRegion in region:
								
							Yields["EE_%s"%bTagRegion] += Yields["EE_%s"%region]
							Yields["EMu_%s"%bTagRegion] += Yields["EMu_%s"%region]
							Yields["MuMu_%s"%bTagRegion] += Yields["MuMu_%s"%region]
							Yields["SFOF_%s"%bTagRegion] += Yields["SFOF_%s"%region]
							Uncertainties["Yield_%s"%bTagRegion] += Uncertainties["Yield_%s"%region]
							
							MCEvents["EE_%s"%bTagRegion] += MCEvents["EE_%s"%region]
							MCEvents["EMu_%s"%bTagRegion] += MCEvents["EMu_%s"%region]
							MCEvents["MuMu_%s"%bTagRegion] += MCEvents["MuMu_%s"%region]
							MCEvents["SFOF_%s"%bTagRegion] += MCEvents["SFOF_%s"%region]
							
							if systematics:
								JES["JESUp_%s"%bTagRegion] +=  JES["JESUp_%s"%region]
								JES["JESDown_%s"%bTagRegion] +=  JES["JESDown_%s"%region]
													
								LeptonFastSim["NoLeptonFastSimScaleFactors_%s"%bTagRegion] += LeptonFastSim["NoLeptonFastSimScaleFactors_%s"%region]
												
								LeptonFullSim["NoLeptonFullSimScaleFactors_%s"%bTagRegion] += LeptonFullSim["NoLeptonFullSimScaleFactors_%s"%region]
								
								Pileup["PileupUp_%s"%bTagRegion] +=  Pileup["PileupUp_%s"%region]
								Pileup["PileupDown_%s"%bTagRegion] +=  Pileup["PileupDown_%s"%region]
								
								ISR["ISRUp_%s"%bTagRegion] += ISR["ISRUp_%s"%region]
								ISR["ISRDown_%s"%bTagRegion] += ISR["ISRDown_%s"%region]		
						
								BTag["BTagHeavy_%s"%bTagRegion] += BTag["BTagHeavy_%s"%region]
								BTag["BTagLight_%s"%bTagRegion] += BTag["BTagLight_%s"%region]		
					
					##Calculate uncertaities for combinations			
					Uncertainties["Efficiency_%s"%bTagRegion] = Yields["SFOF_%s"%bTagRegion]/(xsection*lumi)
					if MCEvents["SFOF_%s"%bTagRegion] > 0:
						Uncertainties["StatUncertainty_%s"%bTagRegion] = sqrt(MCEvents["EE_%s"%bTagRegion]+MCEvents["EMu_%s"%bTagRegion]+MCEvents["MuMu_%s"%bTagRegion])/MCEvents["SFOF_%s"%bTagRegion]
					else:
						Uncertainties["StatUncertainty_%s"%bTagRegion] = 0
							
					
					if systematics:	
						if Yields["SFOF_%s"%bTagRegion] > 0:
							Uncertainties["JESUncertainty_%s"%bTagRegion] = max(abs(JES["JESUp_%s"%bTagRegion]-Yields["SFOF_%s"%bTagRegion])/Yields["SFOF_%s"%bTagRegion],abs(JES["JESDown_%s"%bTagRegion]-Yields["SFOF_%s"%bTagRegion])/Yields["SFOF_%s"%bTagRegion])
							Uncertainties["LeptonFastSimUncertainty_%s"%bTagRegion] = abs(LeptonFastSim["NoLeptonFastSimScaleFactors_%s"%bTagRegion]-Yields["SFOF_%s"%bTagRegion])/Yields["SFOF_%s"%bTagRegion]
							Uncertainties["LeptonFullSimUncertainty_%s"%bTagRegion] = abs(LeptonFullSim["NoLeptonFullSimScaleFactors_%s"%bTagRegion]-Yields["SFOF_%s"%bTagRegion])/Yields["SFOF_%s"%bTagRegion]
							Uncertainties["pileupUncertainty_%s"%bTagRegion] = max(abs(Pileup["PileupUp_%s"%bTagRegion]-Yields["SFOF_%s"%bTagRegion])/Yields["SFOF_%s"%bTagRegion],abs(Pileup["PileupDown_%s"%bTagRegion] -Yields["SFOF_%s"%bTagRegion])/Yields["SFOF_%s"%bTagRegion] )
							Uncertainties["ISRUncertainty_%s"%bTagRegion] = max(abs(ISR["ISRUp_%s"%bTagRegion]-Yields["SFOF_%s"%bTagRegion])/Yields["SFOF_%s"%bTagRegion],abs(ISR["ISRDown_%s"%bTagRegion]-Yields["SFOF_%s"%bTagRegion])/Yields["SFOF_%s"%bTagRegion])
							Uncertainties["BTagUncertaintyHeavy_%s"%bTagRegion] = abs(BTag["BTagHeavy_%s"%bTagRegion]-Yields["SFOF_%s"%bTagRegion])/Yields["SFOF_%s"%bTagRegion]
							Uncertainties["BTagUncertaintyLight_%s"%bTagRegion] = abs(BTag["BTagHeavy_%s"%bTagRegion]-Yields["SFOF_%s"%bTagRegion])/Yields["SFOF_%s"%bTagRegion]
						else:
							Uncertainties["StatUncertainty_%s"%bTagRegion] = 0
							Uncertainties["JESUncertainty_%s"%bTagRegion] = 0
							Uncertainties["LeptonFastSimUncertainty_%s"%bTagRegion] = 0
							Uncertainties["LeptonFullSimUncertainty_%s"%bTagRegion] = 0
							Uncertainties["pileupUncertainty_%s"%bTagRegion] = 0
							Uncertainties["ISRUncertainty_%s"%bTagRegion] = 0
							Uncertainties["BTagUncertaintyHeavy_%s"%bTagRegion] = 0
							Uncertainties["BTagUncertaintyLight_%s"%bTagRegion] = 0
					
						
								
						### Total syst uncertainty
						Uncertainties["SystUncertainty_%s"%bTagRegion] = sqrt(Uncertainties["JESUncertainty_%s"%bTagRegion]**2 + Uncertainties["LeptonFastSimUncertainty_%s"%bTagRegion]**2 + Uncertainties["LeptonFullSimUncertainty_%s"%bTagRegion]**2 + Uncertainties["pileupUncertainty_%s"%bTagRegion]**2 + Uncertainties["ISRUncertainty_%s"%bTagRegion]**2  + Uncertainties["BTagUncertaintyLight_%s"%bTagRegion]**2  + Uncertainties["BTagUncertaintyHeavy_%s"%bTagRegion]**2  + TriggerEffUncertainty**2 + LumiUncertainty**2)
		
						### Total uncertainty
						Uncertainties["TotalUncertainty_%s"%bTagRegion] = sqrt(Uncertainties["SystUncertainty_%s"%bTagRegion]**2 + Uncertainties["StatUncertainty_%s"%bTagRegion]**2)		
					
					for uncertainty in uncertaintySources:
						uncertaintyArrays["%s_%s"%(uncertainty,bTagRegion)].append(Uncertainties["%s_%s"%(uncertainty,bTagRegion)])
													
										
								
							
			m_n += stepsize		
		m_b += stepsize
				
	for region in regions:
		
		for uncertainty in uncertaintySources:
			Graphs["%s_%s"%(uncertainty,region)]=TGraph2D("Graph_%s_%s"%(uncertainty,region),"%s_%s"%(uncertainty,region), len(uncertaintyArrays["%s_%s"%(uncertainty,region)]), array('d',masses_b), array('d',masses_n), array('d',uncertaintyArrays["%s_%s"%(uncertainty,region)]))
			Graphs["%s_%s"%(uncertainty,region)].SetNpx(nxbins)
			Graphs["%s_%s"%(uncertainty,region)].SetNpy(nybins)
			Histograms["%s_%s"%(uncertainty,region)] = Graphs["%s_%s"%(uncertainty,region)].GetHistogram()
			Histograms["%s_%s"%(uncertainty,region)].SetTitle(";m_{#tilde{b}} [GeV]; m_{#tilde{#chi_{2}^{0}}} [GeV]")
		
		if "central" in region:
			region_label = "Central Signal Region"
		elif "forward" in region:
			region_label = "Forward Signal Region"
		
		region_label_2 = ""
		
		if "lowMass" in region:
			region_label_2 = "low mass"	
		elif "belowZ" in region:
			region_label_2 = "below Z"	
		elif "onZ" in region:
			region_label_2 = "on Z"	
		elif "aboveZ" in region:
			region_label_2 = "above Z"	
		elif "highMass" in region:
			region_label_2 = "high mass"	
			
		if "noBTag" in region:
			region_label_2 += ", no b-tag"
		elif "geOneBTag" in region:
			region_label_2 += ", #geq 1 b-tag"		


		plotPad.SetLogz()
		Histograms["Yield_%s"%region].SetZTitle("SF-OF yield")
		Histograms["Yield_%s"%region].Draw("colz")
		latexLumi.DrawLatex(0.85, 0.96, "%s fb^{-1} (13 TeV)"%(printLumi,))
		latexCMS.DrawLatex(0.19,0.89,"CMS")
		latexCMSExtra.DrawLatex(0.19,0.85,"Private Work - Simulation")
		latex.DrawLatex(0.175, 0.75, "#splitline{Simplified Model}{#splitline{T6bbslepton, m(#tilde{#chi}_{1}^{0})=100 GeV}{#splitline{"+region_label+"}{"+region_label_2+"}}}")
		canv.Update()
		canv.Print("fig/Yields/T6bbllslepton_%s.pdf"%(region))
		
		plotPad.SetLogz(0)
		Histograms["StatUncertainty_%s"%region].SetZTitle("SF-OF rel. stat. uncertainty")
		Histograms["StatUncertainty_%s"%region].Draw("colz")
		latexLumi.DrawLatex(0.85, 0.96, "%s fb^{-1} (13 TeV)"%(printLumi,))
		latexCMS.DrawLatex(0.19,0.89,"CMS")
		latexCMSExtra.DrawLatex(0.19,0.85,"Private Work - Simulation")
		latex.DrawLatex(0.175, 0.75, "#splitline{Simplified Model}{#splitline{T6bbslepton, m(#tilde{#chi}_{1}^{0})=100 GeV}{#splitline{"+region_label+"}{"+region_label_2+"}}}")
		canv.Update()
		canv.Print("fig/statUncertainties/T6bbllslepton_%s_stat_err.pdf"%(region))
		
		Histograms["Efficiency_%s"%region].SetZTitle("acceptance #times efficiency")
		Histograms["Efficiency_%s"%region].Draw("colz")
		latexLumi.DrawLatex(0.85, 0.96, "%s fb^{-1} (13 TeV)"%(printLumi,))
		latexCMS.DrawLatex(0.19,0.89,"CMS")
		latexCMSExtra.DrawLatex(0.19,0.85,"Private Work - Simulation")
		latex.DrawLatex(0.175, 0.75, "#splitline{Simplified Model}{#splitline{T6bbslepton, m(#tilde{#chi}_{1}^{0})=100 GeV}{#splitline{"+region_label+"}{"+region_label_2+"}}}")
		canv.Update()
		canv.Print("fig/Efficiencies/T6bbllslepton_%s_signalEfficiency.pdf"%(region))
		
		if systematics:
			Histograms["SystUncertainty_%s"%region].SetZTitle("SF-OF rel. syst. uncertainty")
			Histograms["SystUncertainty_%s"%region].Draw("colz")
			latexLumi.DrawLatex(0.85, 0.96, "%s fb^{-1} (13 TeV)"%(printLumi,))
			latexCMS.DrawLatex(0.19,0.89,"CMS")
			latexCMSExtra.DrawLatex(0.19,0.85,"Private Work - Simulation")
			latex.DrawLatex(0.175, 0.75, "#splitline{Simplified Model}{#splitline{T6bbslepton, m(#tilde{#chi}_{1}^{0})=100 GeV}{#splitline{"+region_label+"}{"+region_label_2+"}}}")
			canv.Update()
			canv.Print("fig/systUncertainties/T6bbllslepton_%s_syst_err.pdf"%(region))
			
			
			Histograms["TotalUncertainty_%s"%region].SetZTitle("SF-OF total rel. uncertainty")
			Histograms["TotalUncertainty_%s"%region].Draw("colz")
			latexLumi.DrawLatex(0.85, 0.96, "%s fb^{-1} (13 TeV)"%(printLumi,))
			latexCMS.DrawLatex(0.19,0.89,"CMS")
			latexCMSExtra.DrawLatex(0.19,0.85,"Private Work - Simulation")
			latex.DrawLatex(0.175, 0.75, "#splitline{Simplified Model}{#splitline{T6bbslepton, m(#tilde{#chi}_{1}^{0})=100 GeV}{#splitline{"+region_label+"}{"+region_label_2+"}}}")
			canv.Update()
			canv.Print("fig/totUncertainties/T6bbllslepton_%s_tot_err.pdf"%(region))
		
			
			for uncertainty in uncertaintySources:
				if not ( uncertainty == "Yield" or uncertainty == "StatUncertainty" or uncertainty == "SystUncertainty" or uncertainty == "TotalUncertainty" or uncertainty == "Efficiency"):			
					Histograms["%s_%s"%(uncertainty,region)].SetZTitle("%s"%uncertainty)
					Histograms["%s_%s"%(uncertainty,region)].Draw("colz")
					latexLumi.DrawLatex(0.85, 0.96, "%s fb^{-1} (13 TeV)"%(printLumi,))
					latexCMS.DrawLatex(0.19,0.89,"CMS")
					latexCMSExtra.DrawLatex(0.19,0.85,"Private Work - Simulation")
					latex.DrawLatex(0.175, 0.75, "#splitline{Simplified Model}{#splitline{T6bbslepton, m(#tilde{#chi}_{1}^{0})=100 GeV}{#splitline{"+region_label+"}{"+region_label_2+"}}}")
					canv.Update()
					canv.Print("fig/%s/T6bbllslepton_%s_%s.pdf"%(uncertainty,region,uncertainty))
				
				
	for regionCombination in regionCombinations:
		
		for uncertainty in uncertaintySources:
			Graphs["%s_%s"%(uncertainty,regionCombination)]=TGraph2D("Graph_%s_%s"%(uncertainty,regionCombination),"%s_%s"%(uncertainty,regionCombination), len(uncertaintyArrays["%s_%s"%(uncertainty,regionCombination)]), array('d',masses_b), array('d',masses_n), array('d',uncertaintyArrays["%s_%s"%(uncertainty,regionCombination)]))
			Graphs["%s_%s"%(uncertainty,regionCombination)].SetNpx(nxbins)
			Graphs["%s_%s"%(uncertainty,regionCombination)].SetNpy(nybins)
			Histograms["%s_%s"%(uncertainty,regionCombination)] = Graphs["%s_%s"%(uncertainty,regionCombination)].GetHistogram()
			Histograms["%s_%s"%(uncertainty,regionCombination)].SetTitle(";m_{#tilde{b}} [GeV]; m_{#tilde{#chi_{2}^{0}}} [GeV]")
		
		if "central" in regionCombination:
			regionCombination_label = "Central Signal Region"
		elif "forward" in regionCombination:
			regionCombination_label = "Forward Signal Region"
		
		regionCombination_label_2 = ""
		
		if "lowMass" in regionCombination:
			regionCombination_label_2 = "low mass"	
		elif "belowZ" in regionCombination:
			regionCombination_label_2 = "below Z"	
		elif "onZ" in regionCombination:
			regionCombination_label_2 = "on Z"	
		elif "aboveZ" in regionCombination:
			regionCombination_label_2 = "above Z"	
		elif "highMass" in regionCombination:
			regionCombination_label_2 = "high mass"	
			
		if "noBTag" in regionCombination:
			if not regionCombination_label_2 == "":
				regionCombination_label_2 += ", no b-tag"
			else:
				regionCombination_label_2 += "no b-tag"
		elif "geOneBTag" in regionCombination:
			if not regionCombination_label_2 == "":
				regionCombination_label_2 += ", #geq 1 b-tag"
			else:
				regionCombination_label_2 += "#geq 1 b-tag"


		plotPad.SetLogz()
		Histograms["Yield_%s"%regionCombination].SetZTitle("SF-OF yield")
		Histograms["Yield_%s"%regionCombination].Draw("colz")
		latexLumi.DrawLatex(0.85, 0.96, "%s fb^{-1} (13 TeV)"%(printLumi,))
		latexCMS.DrawLatex(0.19,0.89,"CMS")
		latexCMSExtra.DrawLatex(0.19,0.85,"Private Work - Simulation")
		latex.DrawLatex(0.175, 0.75, "#splitline{Simplified Model}{#splitline{T6bbslepton, m(#tilde{#chi}_{1}^{0})=100 GeV}{#splitline{"+regionCombination_label+"}{"+regionCombination_label_2+"}}}")
		canv.Update()
		canv.Print("fig/Yields/T6bbllslepton_%s.pdf"%(regionCombination))
		
		f1 = TFile("fig/Yields/TH2_T6bbllslepton_%s_signalYields.root"%(regionCombination),"RECREATE")
		Histograms["Yield_%s"%regionCombination].Write()
		f1.Close()
		
		plotPad.SetLogz(0)
		Histograms["StatUncertainty_%s"%regionCombination].SetZTitle("SF-OF rel. stat. uncertainty")
		Histograms["StatUncertainty_%s"%regionCombination].Draw("colz")
		latexLumi.DrawLatex(0.85, 0.96, "%s fb^{-1} (13 TeV)"%(printLumi,))
		latexCMS.DrawLatex(0.19,0.89,"CMS")
		latexCMSExtra.DrawLatex(0.19,0.85,"Private Work - Simulation")
		latex.DrawLatex(0.175, 0.75, "#splitline{Simplified Model}{#splitline{T6bbslepton, m(#tilde{#chi}_{1}^{0})=100 GeV}{#splitline{"+regionCombination_label+"}{"+regionCombination_label_2+"}}}")
		canv.Update()
		canv.Print("fig/statUncertainties/T6bbllslepton_%s_stat_err.pdf"%(regionCombination))
		
		Histograms["Efficiency_%s"%regionCombination].SetZTitle("acceptance #times efficiency")
		Histograms["Efficiency_%s"%regionCombination].Draw("colz")
		latexLumi.DrawLatex(0.85, 0.96, "%s fb^{-1} (13 TeV)"%(printLumi,))
		latexCMS.DrawLatex(0.19,0.89,"CMS")
		latexCMSExtra.DrawLatex(0.19,0.85,"Private Work - Simulation")
		latex.DrawLatex(0.175, 0.75, "#splitline{Simplified Model}{#splitline{T6bbslepton, m(#tilde{#chi}_{1}^{0})=100 GeV}{#splitline{"+regionCombination_label+"}{"+regionCombination_label_2+"}}}")
		canv.Update()
		canv.Print("fig/Efficiencies/T6bbllslepton_%s_signalEfficiency.pdf"%(regionCombination))
		
		if systematics:
		
			Histograms["SystUncertainty_%s"%regionCombination].SetZTitle("SF-OF rel. syst. uncertainty")
			Histograms["SystUncertainty_%s"%regionCombination].Draw("colz")
			latexLumi.DrawLatex(0.85, 0.96, "%s fb^{-1} (13 TeV)"%(printLumi,))
			latexCMS.DrawLatex(0.19,0.89,"CMS")
			latexCMSExtra.DrawLatex(0.19,0.85,"Private Work - Simulation")
			latex.DrawLatex(0.175, 0.75, "#splitline{Simplified Model}{#splitline{T6bbslepton, m(#tilde{#chi}_{1}^{0})=100 GeV}{#splitline{"+regionCombination_label+"}{"+regionCombination_label_2+"}}}")
			canv.Update()
			canv.Print("fig/systUncertainties/T6bbllslepton_%s_syst_err.pdf"%(regionCombination))
			
			
			Histograms["TotalUncertainty_%s"%regionCombination].SetZTitle("SF-OF total rel. uncertainty")
			Histograms["TotalUncertainty_%s"%regionCombination].Draw("colz")
			latexLumi.DrawLatex(0.85, 0.96, "%s fb^{-1} (13 TeV)"%(printLumi,))
			latexCMS.DrawLatex(0.19,0.89,"CMS")
			latexCMSExtra.DrawLatex(0.19,0.85,"Private Work - Simulation")
			latex.DrawLatex(0.175, 0.75, "#splitline{Simplified Model}{#splitline{T6bbslepton, m(#tilde{#chi}_{1}^{0})=100 GeV}{#splitline{"+regionCombination_label+"}{"+regionCombination_label_2+"}}}")
			canv.Update()
			canv.Print("fig/totUncertainties/T6bbllslepton_%s_tot_err.pdf"%(regionCombination))
			
				
			for uncertainty in uncertaintySources:
				if not ( uncertainty == "Yield" or uncertainty == "StatUncertainty" or uncertainty == "SystUncertainty" or uncertainty == "TotalUncertainty" or uncertainty == "Efficiency"):				
					Histograms["%s_%s"%(uncertainty,regionCombination)].SetZTitle("%s"%uncertainty)
					Histograms["%s_%s"%(uncertainty,regionCombination)].Draw("colz")
					latexLumi.DrawLatex(0.85, 0.96, "%s fb^{-1} (13 TeV)"%(printLumi,))
					latexCMS.DrawLatex(0.19,0.89,"CMS")
					latexCMSExtra.DrawLatex(0.19,0.85,"Private Work - Simulation")
					latex.DrawLatex(0.175, 0.75, "#splitline{Simplified Model}{#splitline{T6bbslepton, m(#tilde{#chi}_{1}^{0})=100 GeV}{#splitline{"+regionCombination_label+"}{"+regionCombination_label_2+"}}}")
					canv.Update()
					canv.Print("fig/%s/T6bbllslepton_%s_%s.pdf"%(uncertainty,regionCombination,uncertainty))
															
	for bTagRegion in bTagRegions:
		
		for uncertainty in uncertaintySources:
			Graphs["%s_%s"%(uncertainty,bTagRegion)]=TGraph2D("Graph_%s_%s"%(uncertainty,bTagRegion),"%s_%s"%(uncertainty,bTagRegion), len(uncertaintyArrays["%s_%s"%(uncertainty,bTagRegion)]), array('d',masses_b), array('d',masses_n), array('d',uncertaintyArrays["%s_%s"%(uncertainty,bTagRegion)]))
			Graphs["%s_%s"%(uncertainty,bTagRegion)].SetNpx(nxbins)
			Graphs["%s_%s"%(uncertainty,bTagRegion)].SetNpy(nybins)
			Histograms["%s_%s"%(uncertainty,bTagRegion)] = Graphs["%s_%s"%(uncertainty,bTagRegion)].GetHistogram()
			Histograms["%s_%s"%(uncertainty,bTagRegion)].SetTitle(";m_{#tilde{b}} [GeV]; m_{#tilde{#chi_{2}^{0}}} [GeV]")
		
		if "central" in bTagRegion:
			bTagRegion_label = "Central Signal Region"
		elif "forward" in bTagRegion:
			bTagRegion_label = "Forward Signal Region"
		else:
			bTagRegion_label = "Signal Region"
		
		bTagRegion_label_2 = ""
		
		if "lowMass" in bTagRegion:
			bTagRegion_label_2 = "low mass"	
		elif "belowZ" in bTagRegion:
			bTagRegion_label_2 = "below Z"	
		elif "onZ" in bTagRegion:
			bTagRegion_label_2 = "on Z"	
		elif "aboveZ" in bTagRegion:
			bTagRegion_label_2 = "above Z"	
		elif "highMass" in bTagRegion:
			bTagRegion_label_2 = "high mass"	
			
		if "noBTag" in bTagRegion:
			if not bTagRegion_label_2 == "":
				bTagRegion_label_2 += ", no b-tag"
			else:
				bTagRegion_label_2 += "no b-tag"
		elif "geOneBTag" in bTagRegion:
			if not bTagRegion_label_2 == "":
				bTagRegion_label_2 += ", #geq 1 b-tag"
			else:
				bTagRegion_label_2 += "#geq 1 b-tag"


		plotPad.SetLogz()
		Histograms["Yield_%s"%bTagRegion].SetZTitle("SF-OF yield")
		Histograms["Yield_%s"%bTagRegion].Draw("colz")
		latexLumi.DrawLatex(0.85, 0.96, "%s fb^{-1} (13 TeV)"%(printLumi,))
		latexCMS.DrawLatex(0.19,0.89,"CMS")
		latexCMSExtra.DrawLatex(0.19,0.85,"Private Work - Simulation")
		latex.DrawLatex(0.175, 0.75, "#splitline{Simplified Model}{#splitline{T6bbslepton, m(#tilde{#chi}_{1}^{0})=100 GeV}{#splitline{"+bTagRegion_label+"}{"+bTagRegion_label_2+"}}}")
		canv.Update()
		canv.Print("fig/Yields/T6bbllslepton_%s.pdf"%(bTagRegion))
		
		f1 = TFile("fig/Yields/TH2_T6bbllslepton_%s_signalYields.root"%(bTagRegion),"RECREATE")
		Histograms["Yield_%s"%bTagRegion].Write()
		f1.Close()
		
		plotPad.SetLogz(0)
		Histograms["StatUncertainty_%s"%bTagRegion].SetZTitle("SF-OF rel. stat. uncertainty")
		Histograms["StatUncertainty_%s"%bTagRegion].Draw("colz")
		latexLumi.DrawLatex(0.85, 0.96, "%s fb^{-1} (13 TeV)"%(printLumi,))
		latexCMS.DrawLatex(0.19,0.89,"CMS")
		latexCMSExtra.DrawLatex(0.19,0.85,"Private Work - Simulation")
		latex.DrawLatex(0.175, 0.75, "#splitline{Simplified Model}{#splitline{T6bbslepton, m(#tilde{#chi}_{1}^{0})=100 GeV}{#splitline{"+bTagRegion_label+"}{"+bTagRegion_label_2+"}}}")
		canv.Update()
		canv.Print("fig/statUncertainties/T6bbllslepton_%s_stat_err.pdf"%(bTagRegion))
		
		Histograms["Efficiency_%s"%bTagRegion].SetZTitle("acceptance #times efficiency")
		Histograms["Efficiency_%s"%bTagRegion].Draw("colz")
		latexLumi.DrawLatex(0.85, 0.96, "%s fb^{-1} (13 TeV)"%(printLumi,))
		latexCMS.DrawLatex(0.19,0.89,"CMS")
		latexCMSExtra.DrawLatex(0.19,0.85,"Private Work - Simulation")
		latex.DrawLatex(0.175, 0.75, "#splitline{Simplified Model}{#splitline{T6bbslepton, m(#tilde{#chi}_{1}^{0})=100 GeV}{#splitline{"+bTagRegion_label+"}{"+bTagRegion_label_2+"}}}")
		canv.Update()
		canv.Print("fig/Efficiencies/T6bbllslepton_%s_signalEfficiency.pdf"%(bTagRegion))
		
		if systematics:
		
			Histograms["SystUncertainty_%s"%bTagRegion].SetZTitle("SF-OF rel. syst. uncertainty")
			Histograms["SystUncertainty_%s"%bTagRegion].Draw("colz")
			latexLumi.DrawLatex(0.85, 0.96, "%s fb^{-1} (13 TeV)"%(printLumi,))
			latexCMS.DrawLatex(0.19,0.89,"CMS")
			latexCMSExtra.DrawLatex(0.19,0.85,"Private Work - Simulation")
			latex.DrawLatex(0.175, 0.75, "#splitline{Simplified Model}{#splitline{T6bbslepton, m(#tilde{#chi}_{1}^{0})=100 GeV}{#splitline{"+bTagRegion_label+"}{"+bTagRegion_label_2+"}}}")
			canv.Update()
			canv.Print("fig/systUncertainties/T6bbllslepton_%s_syst_err.pdf"%(bTagRegion))
			
			
			Histograms["TotalUncertainty_%s"%bTagRegion].SetZTitle("SF-OF total rel. uncertainty")
			Histograms["TotalUncertainty_%s"%bTagRegion].Draw("colz")
			latexLumi.DrawLatex(0.85, 0.96, "%s fb^{-1} (13 TeV)"%(printLumi,))
			latexCMS.DrawLatex(0.19,0.89,"CMS")
			latexCMSExtra.DrawLatex(0.19,0.85,"Private Work - Simulation")
			latex.DrawLatex(0.175, 0.75, "#splitline{Simplified Model}{#splitline{T6bbslepton, m(#tilde{#chi}_{1}^{0})=100 GeV}{#splitline{"+bTagRegion_label+"}{"+bTagRegion_label_2+"}}}")
			canv.Update()
			canv.Print("fig/totUncertainties/T6bbllslepton_%s_tot_err.pdf"%(bTagRegion))
			
				
			for uncertainty in uncertaintySources:
				if not ( uncertainty == "Yield" or uncertainty == "StatUncertainty" or uncertainty == "SystUncertainty" or uncertainty == "TotalUncertainty" or uncertainty == "Efficiency"):				
					Histograms["%s_%s"%(uncertainty,bTagRegion)].SetZTitle("%s"%uncertainty)
					Histograms["%s_%s"%(uncertainty,bTagRegion)].Draw("colz")
					latexLumi.DrawLatex(0.85, 0.96, "%s fb^{-1} (13 TeV)"%(printLumi,))
					latexCMS.DrawLatex(0.19,0.89,"CMS")
					latexCMSExtra.DrawLatex(0.19,0.85,"Private Work - Simulation")
					latex.DrawLatex(0.175, 0.75, "#splitline{Simplified Model}{#splitline{T6bbslepton, m(#tilde{#chi}_{1}^{0})=100 GeV}{#splitline{"+bTagRegion_label+"}{"+bTagRegion_label_2+"}}}")
					canv.Update()
					canv.Print("fig/%s/T6bbllslepton_%s_%s.pdf"%(uncertainty,bTagRegion,uncertainty))
								
								
								
					
			
										
					
			
											
									
				

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
    parser.add_option("-s", "--systematics", dest="systematics", action="store_true", default=False,
                                  help="plot maps for systematic uncertainties")
    (opts, args) = parser.parse_args()
    if (opts.verbose):
        # print out all output
        log.outputLevel = 5
    else:
        # ignore output with "debug" level
        log.outputLevel = 4

    # start
    plot(opts.systematics)
 
