#! /usr/bin/env python

import ROOT
ROOT.gROOT.SetBatch(True)
import sys
sys.path.append("/.automount/home/home__home4/institut_1b/schomakers/FrameWork/frameWorkBase")
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
	from ROOT import TCanvas, TPad, TH1F, TH2F,TH2D, TH1I, THStack, TLegend, TMath, TF1, TFile, TGraph2D
	from setTDRStyle import setTDRStyle
	import pickle
	from defs import sbottom_masses
	from math import sqrt
	
	from corrections import rSFOFDirect
	
	signalDenominatorFile = TFile("/.automount/home/home__home4/institut_1b/schomakers/FrameWork/SignalScan/T6bbllsleptonDenominatorHisto7.root")
	denominatorHisto = TH2F(signalDenominatorFile.Get("massScan"))
	ISRNormalizationHisto = TH2F(signalDenominatorFile .Get("ISRNormalization"))
	ISRNormalizationHistoUp = TH2F(signalDenominatorFile .Get("ISRNormalizationUp"))
	ISRNormalizationHistoDown = TH2F(signalDenominatorFile .Get("ISRNormalizationDown"))
	
	
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
	lumi = 35867.
	printLumi = "35.9"
	
	m_neutr_1_fix = False
	#~ m_neutr_1_fix = True

	
	
	path = "shelvesSystematics"
	generalSignalLabel = "T6bbllslepton"
	
	observables = ["mll"]
	leptonCombinations = ["SF-OF"]
	
	massRegions = ["20To60","60To86","86To96","96To150","150To200","200To300","300To400","Above400"]
	nLLRegions = ["lowNll","highNll"]
	MT2Regions = ["highMT2"]
	
	regionCombinations = []
	regions = []
	
	for massRegion in massRegions:
		for nLLRegion in nLLRegions:
			for MT2Region in MT2Regions:
				regions.append("%s_%s_%s"%(massRegion,nLLRegion,MT2Region))
				
	for nLLRegion in nLLRegions:
		regionCombinations.append(nLLRegion)	
			
	
	Graphs = {}
	Histograms = {}			
	uncertaintyArrays = {}
	#~ uncertaintySources = ["Yield","StatUncertainty","SystUncertainty","TotalUncertainty","Efficiency","ISRUncertainty","pileupUncertainty","JESUncertainty","LeptonUncertainty","LeptonFastSimUncertainty","LeptonFullSimUncertainty","BTagUncertaintyLight","BTagUncertaintyHeavy"]
	uncertaintySources = ["Yield","StatUncertainty","SystUncertainty","TotalUncertainty","Efficiency","EfficiencyUnscaled","ISRUncertainty","pileupUncertainty","JESUncertainty","LeptonFullSimUncertainty","LeptonFastSimUncertainty","BTagUncertaintyLight","BTagUncertaintyHeavy","MetUncertainty","ScaleUncertainty"]

	
	for uncertaintySource in uncertaintySources:
		uncertaintyArrays["%s_highNll"%(uncertaintySource)] = []
		uncertaintyArrays["%s_lowNll"%(uncertaintySource)] = []
		for region in regions:
			uncertaintyArrays["%s_%s"%(uncertaintySource,region)] = []
			
	
	
	title = "Simplified Model Scan; m(#tilde{b}) [GeV]; m(#tilde{#chi}_{2}^{0}) [GeV]"

	masses_b = []
	masses_n = []
	
	m_n_min = 150
	m_n_max = 1450
	m_b_min = 700
	m_b_max = 1500
	
	bin_size =25
	nxbins = int(min(500,(m_b_max-m_b_min)/bin_size))
	nybins = int(min(500,(m_n_max-m_n_min)/bin_size))
	
	
	TriggerEffUncertainty = 0.03
	LumiUncertainty = 0.026
	FastSimUncertainty = 0.04
					
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
			
			if not ((m_b == 775 and m_n == 750) or (m_b == 800 and m_n == 150) or (m_b == 950 and m_n == 250) or (m_b == 950 and m_n == 300) or (m_b == 950 and m_n == 500) or (m_b == 950 and m_n == 550) or (m_b == 950 and m_n == 850) or (m_b == 950 and m_n == 900)):	
							
				
				m_neutralino_2 = str(m_n)
				
				masses_b.append(m_b)
				masses_n.append(m_n)
				
				sampleName = "T6bbllslepton_msbottom_%s_mneutralino_%s"%(m_sbottom,m_neutralino_2)
				
				denominator = denominatorHisto.GetBinContent(denominatorHisto.GetXaxis().FindBin(int(sampleName.split("_")[2])),denominatorHisto.GetYaxis().FindBin(int(sampleName.split("_")[4])))
				ISRNormalization = ISRNormalizationHisto.GetBinContent(ISRNormalizationHisto.GetXaxis().FindBin(int(sampleName.split("_")[2])),denominatorHisto.GetYaxis().FindBin(int(sampleName.split("_")[4])))
					
				
	
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
				ScaleShift = {}
				
				Uncertainties = {}	
				
				Pickles["%s_%s"%(m_sbottom,m_neutralino_2)] = loadPickles("%s/%s.pkl"%(path,sampleName))	
				
				
				for region in regions:
					
					if region == "edgeLegacy":
						RSFOF = rSFOFDirect.central.val
					else:
						RSFOF = rSFOFDirect.inclusive.val
	
					
					Yields["EE_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_EE"]["Val"]
					Yields["EMu_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_EMu"]["Val"]*RSFOF
					Yields["MuMu_%s"%region] =  Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_MuMu"]["Val"]
					Yields["SFOF_%s"%region] = Yields["EE_%s"%region] + Yields["MuMu_%s"%region] - Yields["EMu_%s"%region]
					Uncertainties["Yield_%s"%region] = max(Yields["SFOF_%s"%region],0)
					
					MCEvents["EE_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_EE"]["MCEvents"]
					MCEvents["EMu_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_EMu"]["MCEvents"]
					MCEvents["MuMu_%s"%region] =  Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_MuMu"]["MCEvents"]
					MCEvents["SFOF_%s"%region] = MCEvents["EE_%s"%region] + MCEvents["MuMu_%s"%region] - MCEvents["EMu_%s"%region]
					MCEvents["SFOF_%s"%region] = max(MCEvents["SFOF_%s"%region],0)
					
					Uncertainties["Efficiency_%s"%region] = Yields["SFOF_%s"%region]/(xsection*lumi*ISRNormalization)
					Uncertainties["EfficiencyUnscaled_%s"%region] = MCEvents["SFOF_%s"%region]/denominator
					
					if MCEvents["SFOF_%s"%region] > 0:
						Uncertainties["StatUncertainty_%s"%region] = sqrt(MCEvents["EE_%s"%region]+MCEvents["EMu_%s"%region]+MCEvents["MuMu_%s"%region])/MCEvents["SFOF_%s"%region]
					else:
						Uncertainties["StatUncertainty_%s"%region] = 0
					
					### JES Uncertainty
					JES["Mean_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_EE"]["JESMean"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_MuMu"]["JESMean"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_EMu"]["JESMean"]*RSFOF 
					
					JES["JESUp_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_EE"]["JESUp"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_MuMu"]["JESUp"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_EMu"]["JESUp"]*RSFOF
					JES["JESDown_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_EE"]["JESDown"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_MuMu"]["JESDown"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_EMu"]["JESDown"]*RSFOF
					
					### Only consider bins with a sufficient number of MC events
					### otherwise nearly empty bins will dominate the plots
					if JES["Mean_%s"%region] > 0 and MCEvents["SFOF_%s"%region] > 100:
						Uncertainties["JESUncertainty_%s"%region] = max(abs(JES["JESUp_%s"%region]-JES["Mean_%s"%region])/JES["Mean_%s"%region],abs(JES["JESDown_%s"%region]-JES["Mean_%s"%region])/JES["Mean_%s"%region])
					else:
						Uncertainties["JESUncertainty_%s"%region] = 0
					
					### Lepton FastSim Uncertainty
					
					LeptonFastSim["Mean_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EE"%region]["LeptonFastSimMean"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_MuMu"%region]["LeptonFastSimMean"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EMu"%region]["LeptonFastSimMean"] 
					LeptonFastSim["MeanShifted_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EE"%region]["LeptonFastSimMean"] * (1+FastSimUncertainty) + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_MuMu"%region]["LeptonFastSimMean"]*(1+FastSimUncertainty)  - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EMu"%region]["LeptonFastSimMean"]*(1+FastSimUncertainty)  
					
					if LeptonFastSim["Mean_%s"%region] > 0:
						Uncertainties["LeptonFastSimUncertainty_%s"%region] = abs(LeptonFastSim["MeanShifted_%s"%region]-LeptonFastSim["Mean_%s"%region])/LeptonFastSim["Mean_%s"%region]
					else:
						Uncertainties["LeptonFastSimUncertainty_%s"%region] = 0
						
					### Lepton FullSim Uncertainty
					
					LeptonFullSim["Mean_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EE"%region]["LeptonFullSimMean"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_MuMu"%region]["LeptonFullSimMean"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EMu"%region]["LeptonFullSimMean"] 
					LeptonFullSim["MeanShifted_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EE"%region]["LeptonFullSimScaleFactorShifted"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_MuMu"%region]["LeptonFullSimScaleFactorShifted"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EMu"%region]["LeptonFullSimScaleFactorShifted"] 
					
					if LeptonFullSim["Mean_%s"%region] > 0:
						Uncertainties["LeptonFullSimUncertainty_%s"%region] = abs(LeptonFullSim["MeanShifted_%s"%region]-LeptonFullSim["Mean_%s"%region])/LeptonFullSim["Mean_%s"%region]
					else:
						Uncertainties["LeptonFullSimUncertainty_%s"%region] = 0
					
					
					###  Pileup Uncertainty
					Pileup["Mean_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_EE"]["PileupMean"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_MuMu"]["PileupMean"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_EMu"]["PileupMean"]*RSFOF
					
					Pileup["PileupHigh_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_EE"]["PileupHigh"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_MuMu"]["PileupHigh"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_EMu"]["PileupHigh"]*RSFOF
					Pileup["PileupLow_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_EE"]["PileupLow"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_MuMu"]["PileupLow"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_EMu"]["PileupLow"]*RSFOF 
					
					if Pileup["Mean_%s"%region] > 0 and MCEvents["SFOF_%s"%region] > 100:
						Uncertainties["pileupUncertainty_%s"%region] = max(abs(Pileup["PileupHigh_%s"%region]-Pileup["Mean_%s"%region])/Pileup["Mean_%s"%region],abs(Pileup["PileupLow_%s"%region] -Pileup["Mean_%s"%region])/Pileup["Mean_%s"%region] )
					else:
						Uncertainties["pileupUncertainty_%s"%region] = 0
					
					Uncertainties["pileupUncertainty_%s"%region] = 0.02
					
					### ISR Uncertainty
					ISR["Mean_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_EE"]["ISRMean"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_MuMu"]["ISRMean"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_EMu"]["ISRMean"]*RSFOF
					
					ISR["ISRUp_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_EE"]["ISRUp"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_MuMu"]["ISRUp"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_EMu"]["ISRUp"]*RSFOF
					ISR["ISRDown_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_EE"]["ISRDown"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_MuMu"]["ISRDown"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_EMu"]["ISRDown"]*RSFOF  
					
					if ISR["Mean_%s"%region] > 0 and MCEvents["SFOF_%s"%region] > 100:
						Uncertainties["ISRUncertainty_%s"%region] = max(abs(ISR["ISRUp_%s"%region]-ISR["Mean_%s"%region])/ISR["Mean_%s"%region],abs(ISR["ISRDown_%s"%region]-ISR["Mean_%s"%region])/ISR["Mean_%s"%region])
					else:
						Uncertainties["ISRUncertainty_%s"%region] = 0
						
					if Yields["SFOF_%s"%region] > 0 and MCEvents["SFOF_%s"%region] > 100:
						Uncertainties["TriggerEffUncertainty_%s"%region] = sqrt((Yields["EE_%s"%region]*TriggerEffUncertainty)**2 + (Yields["MuMu_%s"%region]*TriggerEffUncertainty)**2 + (Yields["EMu_%s"%region]*TriggerEffUncertainty*RSFOF)**2)/Yields["SFOF_%s"%region]
					else:
						Uncertainties["TriggerEffUncertainty_%s"%region] = 0
							
					### BTag Uncertainty
					BTag["Mean_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_EE"]["BTagMean"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_MuMu"]["BTagMean"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_EMu"]["BTagMean"]*RSFOF
					
					BTag["BTagHeavy_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_EE"]["BTagHeavy"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_MuMu"]["BTagHeavy"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_EMu"]["BTagHeavy"]*RSFOF  
					
					BTag["BTagLight_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_EE"]["BTagLight"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_MuMu"]["BTagLight"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_EMu"]["BTagLight"]*RSFOF
					
					if BTag["Mean_%s"%region] > 0 and MCEvents["SFOF_%s"%region] >100:
						Uncertainties["BTagUncertaintyHeavy_%s"%region] = abs(BTag["BTagHeavy_%s"%region]-BTag["Mean_%s"%region])/BTag["Mean_%s"%region]
						Uncertainties["BTagUncertaintyLight_%s"%region] = abs(BTag["BTagLight_%s"%region]-BTag["Mean_%s"%region])/BTag["Mean_%s"%region]
					else:
						Uncertainties["BTagUncertaintyHeavy_%s"%region] = 0
						Uncertainties["BTagUncertaintyLight_%s"%region] = 0
						
					### FastSim Met Uncertainty					
					Met["Met_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_EE"]["Met"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_MuMu"]["Met"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_EMu"]["Met"]*RSFOF
					Met["GenMet_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_EE"]["GenMet"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_MuMu"]["GenMet"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)][region+"_EMu"]["GenMet"]*RSFOF  
					
					if Yields["SFOF_%s"%region] > 0 and MCEvents["SFOF_%s"%region] >100:
						Uncertainties["MetUncertainty_%s"%region] = 0.5*abs(Met["Met_%s"%region]-Met["GenMet_%s"%region])/Yields["SFOF_%s"%region]
					else:
						Uncertainties["MetUncertainty_%s"%region] = 0
					
					### Scale uncertainty
					Uncertainties["ScaleUncertainty_%s"%region] = 0
					ScaleShift["Mean_%s"%region] = BTag["Mean_%s"%region]
					ScaleShift["MaxShift_%s"%region] = 0
					if BTag["Mean_%s"%region] > 0:
						
						for scaleIndex in range(1,9):
							ScaleShift["ScaleShift%s"%str(scaleIndex)] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EE"%region]["ScaleShifted%s"%str(scaleIndex)] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_MuMu"%region]["ScaleShifted%s"%str(scaleIndex)] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EMu"%region]["ScaleShifted%s"%str(scaleIndex)]
							Uncertainties["ScaleUncertainty_%s"%region] = max(Uncertainties["ScaleUncertainty_%s"%region], abs( (ScaleShift["ScaleShift%s"%str(scaleIndex)] - ScaleShift["Mean_%s"%region])/ScaleShift["Mean_%s"%region]))
							if abs( ScaleShift["ScaleShift%s"%str(scaleIndex)] - ScaleShift["Mean_%s"%region])/ScaleShift["Mean_%s"%region] > ScaleShift["MaxShift_%s"%region]:
								ScaleShift["MaxShift_%s"%region] = ScaleShift["ScaleShift%s"%str(scaleIndex)]
								
					
					#~ ### Trigger efficiency	
					#~ if Yields["SFOF_%s"%region] > 0:
						#~ Uncertainties["TriggerEffUncertainty_%s"%region] = sqrt((Yields["EE_%s"%region]*TriggerEffUncertainty)**2 + (Yields["MuMu_%s"%region]*TriggerEffUncertainty)**2 + (Yields["EMu_%s"%region]*TriggerEffUncertainty)**2)/Yields["SFOF_%s"%region]
					#~ else:
						#~ Uncertainties["TriggerEffUncertainty_%s"%region] = 0
					Uncertainties["TriggerEffUncertainty_%s"%region] = TriggerEffUncertainty
							
					### Total syst uncertainty
					Uncertainties["SystUncertainty_%s"%region] = sqrt(Uncertainties["ScaleUncertainty_%s"%region]**2 + Uncertainties["JESUncertainty_%s"%region]**2 + Uncertainties["LeptonFastSimUncertainty_%s"%region]**2 + Uncertainties["LeptonFullSimUncertainty_%s"%region]**2 + Uncertainties["pileupUncertainty_%s"%region]**2 + Uncertainties["ISRUncertainty_%s"%region]**2 + Uncertainties["BTagUncertaintyLight_%s"%region]**2 + Uncertainties["BTagUncertaintyHeavy_%s"%region]**2  + Uncertainties["TriggerEffUncertainty_%s"%region]**2 + LumiUncertainty**2)
					
					### Totalaluncertainty
					Uncertainties["TotalUncertainty_%s"%region] = sqrt(Uncertainties["SystUncertainty_%s"%region]**2 + Uncertainties["StatUncertainty_%s"%region]**2)
			
					for uncertainty in uncertaintySources:
						uncertaintyArrays["%s_%s"%(uncertainty,region)].append(Uncertainties["%s_%s"%(uncertainty,region)])
						
				for regionCombination in regionCombinations:
						
						Yields["EE_%s"%regionCombination] = 0	
						Yields["EMu_%s"%regionCombination] = 0	
						Yields["MuMu_%s"%regionCombination] = 0	
						Yields["SFOF_%s"%regionCombination] = 0
						
						Uncertainties["Yield_%s"%regionCombination] = 0
						
						MCEvents["EE_%s"%regionCombination] = 0
						MCEvents["EMu_%s"%regionCombination] = 0
						MCEvents["MuMu_%s"%regionCombination] = 0
						MCEvents["SFOF_%s"%regionCombination] = 0
						
						JES["Mean_%s"%regionCombination] =  0
						JES["JESUp_%s"%regionCombination] =  0
						JES["JESDown_%s"%regionCombination] =  0
						
						LeptonFastSim["Mean_%s"%regionCombination] =  0
						LeptonFastSim["MeanShifted_%s"%regionCombination] =  0
						
						LeptonFullSim["Mean_%s"%regionCombination] =  0
						LeptonFullSim["MeanShifted_%s"%regionCombination] =  0
						
						Pileup["Mean_%s"%regionCombination] = 0
						Pileup["PileupHigh_%s"%regionCombination] =  0
						Pileup["PileupLow_%s"%regionCombination] =  0
						
						ISR["Mean_%s"%regionCombination] = 0					
						ISR["ISRUp_%s"%regionCombination] = 0
						ISR["ISRDown_%s"%regionCombination] = 0
							
						BTag["Mean_%s"%regionCombination] = 0					
						BTag["BTagHeavy_%s"%regionCombination] = 0
						BTag["BTagLight_%s"%regionCombination] = 0
							
						Met["Met_%s"%regionCombination] = 0					
						Met["GenMet_%s"%regionCombination] = 0
						
						ScaleShift["Mean_%s"%regionCombination] = 0
						ScaleShift["MaxShift_%s"%regionCombination] = 0
						
						
					### Add information of single regions into combinations
						for region in regions:
							#~ etaBin, massBin, bTagBin = region.split('_')	
							#~ if etaBin in regionCombination or massBin in regionCombination or bTagBin in regionCombination:
								
							if regionCombination in region:
									
								Yields["EE_%s"%regionCombination] += Yields["EE_%s"%region]
								Yields["EMu_%s"%regionCombination] += Yields["EMu_%s"%region]
								Yields["MuMu_%s"%regionCombination] += Yields["MuMu_%s"%region]
								Yields["SFOF_%s"%regionCombination] += Yields["SFOF_%s"%region]
								Uncertainties["Yield_%s"%regionCombination] += Uncertainties["Yield_%s"%region]
								
								MCEvents["EE_%s"%regionCombination] += MCEvents["EE_%s"%region]
								MCEvents["EMu_%s"%regionCombination] += MCEvents["EMu_%s"%region]
								MCEvents["MuMu_%s"%regionCombination] += MCEvents["MuMu_%s"%region]
								MCEvents["SFOF_%s"%regionCombination] += MCEvents["SFOF_%s"%region]
								
								JES["Mean_%s"%regionCombination] +=  JES["Mean_%s"%region]
								JES["JESUp_%s"%regionCombination] +=  JES["JESUp_%s"%region]
								JES["JESDown_%s"%regionCombination] +=  JES["JESDown_%s"%region]
								
								LeptonFastSim["Mean_%s"%regionCombination] +=  LeptonFastSim["Mean_%s"%region]
								LeptonFastSim["MeanShifted_%s"%regionCombination] +=  LeptonFastSim["MeanShifted_%s"%region]
								
								LeptonFullSim["Mean_%s"%regionCombination] +=  LeptonFullSim["Mean_%s"%region]
								LeptonFullSim["MeanShifted_%s"%regionCombination] +=  LeptonFullSim["MeanShifted_%s"%region]
								
								Pileup["Mean_%s"%regionCombination] +=  Pileup["Mean_%s"%region]						
								Pileup["PileupHigh_%s"%regionCombination] +=  Pileup["PileupHigh_%s"%region]
								Pileup["PileupLow_%s"%regionCombination] +=  Pileup["PileupLow_%s"%region]
								
								ISR["Mean_%s"%regionCombination] += ISR["Mean_%s"%region]						
								ISR["ISRUp_%s"%regionCombination] += ISR["ISRUp_%s"%region]
								ISR["ISRDown_%s"%regionCombination] += ISR["ISRDown_%s"%region]		
						
								BTag["Mean_%s"%regionCombination] += BTag["Mean_%s"%region]						
								BTag["BTagHeavy_%s"%regionCombination] += BTag["BTagHeavy_%s"%region]
								BTag["BTagLight_%s"%regionCombination] += BTag["BTagLight_%s"%region]		
						
								Met["Met_%s"%regionCombination] += Met["Met_%s"%region]						
								Met["GenMet_%s"%regionCombination] += Met["GenMet_%s"%region]
								
								ScaleShift["Mean_%s"%regionCombination] += ScaleShift["Mean_%s"%region]
								ScaleShift["MaxShift_%s"%regionCombination] += ScaleShift["MaxShift_%s"%region]
						
						##Calculate uncertaities for combinations			
						Uncertainties["Efficiency_%s"%regionCombination] = Yields["SFOF_%s"%regionCombination]/(xsection*lumi*ISRNormalization)
						Uncertainties["EfficiencyUnscaled_%s"%regionCombination] = MCEvents["SFOF_%s"%regionCombination]/denominator
	
							
						if MCEvents["SFOF_%s"%regionCombination] > 0:
							Uncertainties["StatUncertainty_%s"%regionCombination] = sqrt(MCEvents["EE_%s"%regionCombination]+MCEvents["EMu_%s"%regionCombination]+MCEvents["MuMu_%s"%regionCombination])/MCEvents["SFOF_%s"%regionCombination]
						else:
							Uncertainties["StatUncertainty_%s"%regionCombination] = 0
						
						if JES["Mean_%s"%regionCombination] > 0 and MCEvents["SFOF_%s"%regionCombination] > 100:
						#~ if JES["Mean_%s"%regionCombination] > 0:
							Uncertainties["JESUncertainty_%s"%regionCombination] = max(abs(JES["JESUp_%s"%regionCombination]-JES["Mean_%s"%regionCombination])/JES["Mean_%s"%regionCombination],abs(JES["JESDown_%s"%regionCombination]-JES["Mean_%s"%regionCombination])/JES["Mean_%s"%regionCombination])
						else:
							Uncertainties["JESUncertainty_%s"%regionCombination] = 0
							
						if LeptonFastSim["Mean_%s"%regionCombination] > 0 and MCEvents["SFOF_%s"%regionCombination] > 100:
						#~ if LeptonFastSim["Mean_%s"%regionCombination] > 0:
							Uncertainties["LeptonFastSimUncertainty_%s"%regionCombination] = abs(LeptonFastSim["MeanShifted_%s"%regionCombination]-LeptonFastSim["Mean_%s"%regionCombination])/LeptonFastSim["Mean_%s"%regionCombination]
						else:
							Uncertainties["LeptonFastSimUncertainty_%s"%regionCombination] = 0
						
						if LeptonFullSim["Mean_%s"%regionCombination] > 0 and MCEvents["SFOF_%s"%regionCombination] > 100:
						#~ if LeptonFullSim["Mean_%s"%regionCombination] > 0:
							Uncertainties["LeptonFullSimUncertainty_%s"%regionCombination] = abs(LeptonFullSim["MeanShifted_%s"%regionCombination]-LeptonFullSim["Mean_%s"%regionCombination])/LeptonFullSim["Mean_%s"%regionCombination]
						else:
							Uncertainties["LeptonFullSimUncertainty_%s"%regionCombination] = 0
						
						if Pileup["Mean_%s"%regionCombination] > 0 and MCEvents["SFOF_%s"%regionCombination] > 100:
						#~ if Pileup["Mean_%s"%regionCombination] > 0:
							Uncertainties["pileupUncertainty_%s"%regionCombination] = max(abs(Pileup["PileupHigh_%s"%regionCombination]-Pileup["Mean_%s"%regionCombination])/Pileup["Mean_%s"%regionCombination],abs(Pileup["PileupLow_%s"%regionCombination] -Pileup["Mean_%s"%regionCombination])/Pileup["Mean_%s"%regionCombination] )
						else:
							Uncertainties["pileupUncertainty_%s"%regionCombination] = 0
							
						Uncertainties["pileupUncertainty_%s"%regionCombination] = 0.02
						
						if ISR["Mean_%s"%regionCombination] > 0 and MCEvents["SFOF_%s"%regionCombination] > 100:
						#~ if ISR["Mean_%s"%regionCombination] > 0:
							Uncertainties["ISRUncertainty_%s"%regionCombination] = max(abs(ISR["ISRUp_%s"%regionCombination]-ISR["Mean_%s"%regionCombination])/ISR["Mean_%s"%regionCombination],abs(ISR["ISRDown_%s"%regionCombination]-ISR["Mean_%s"%regionCombination])/ISR["Mean_%s"%regionCombination])
						else:
							Uncertainties["ISRUncertainty_%s"%regionCombination] = 0
							
						if BTag["Mean_%s"%regionCombination] > 0 and MCEvents["SFOF_%s"%regionCombination] > 100:
						#~ if BTag["Mean_%s"%regionCombination] > 0:
							Uncertainties["BTagUncertaintyHeavy_%s"%regionCombination] = abs(BTag["BTagHeavy_%s"%regionCombination]-BTag["Mean_%s"%regionCombination])/BTag["Mean_%s"%regionCombination]
							Uncertainties["BTagUncertaintyLight_%s"%regionCombination] = abs(BTag["BTagLight_%s"%regionCombination]-BTag["Mean_%s"%regionCombination])/BTag["Mean_%s"%regionCombination]
						else:
							Uncertainties["BTagUncertaintyHeavy_%s"%regionCombination] = 0
							Uncertainties["BTagUncertaintyLight_%s"%regionCombination] = 0
							
						if Met["Met_%s"%regionCombination] > 0  and MCEvents["SFOF_%s"%regionCombination] > 100:
						#~ if Met["Met_%s"%regionCombination] > 0:
							Uncertainties["MetUncertainty_%s"%regionCombination] = 0.5*abs(Met["Met_%s"%regionCombination]-Met["GenMet_%s"%regionCombination])/Yields["SFOF_%s"%regionCombination]
						else:
							Uncertainties["MetUncertainty_%s"%regionCombination] = 0
					
						if ScaleShift["Mean_%s"%regionCombination] > 0  and MCEvents["SFOF_%s"%regionCombination] > 100:
						#~ if ScaleShift["Mean_%s"%regionCombination] > 0:
							Uncertainties["ScaleUncertainty_%s"%regionCombination] = abs( (ScaleShift["MaxShift_%s"%regionCombination] - ScaleShift["Mean_%s"%regionCombination])/ScaleShift["Mean_%s"%regionCombination])
						else:
							Uncertainties["ScaleUncertainty_%s"%regionCombination] = 0
							
							
						#~ if Yields["SFOF_%s"%regionCombination] > 0:
							#~ Uncertainties["TriggerEffUncertainty_%s"%regionCombination] = sqrt((Yields["EE_%s"%regionCombination]*TriggerEffUncertainty)**2 + (Yields["MuMu_%s"%regionCombination]*TriggerEffUncertainty)**2 + (Yields["EMu_%s"%regionCombination]*TriggerEffUncertainty)**2)/Yields["SFOF_%s"%regionCombination]
						#~ else:
							#~ Uncertainties["TriggerEffUncertainty_%s"%regionCombination] = 0
						Uncertainties["TriggerEffUncertainty_%s"%regionCombination] = 0.05
								
						### Total syst uncertainty
						Uncertainties["SystUncertainty_%s"%regionCombination] = sqrt(Uncertainties["ScaleUncertainty_%s"%regionCombination]**2 + Uncertainties["JESUncertainty_%s"%regionCombination]**2 + Uncertainties["LeptonFastSimUncertainty_%s"%regionCombination]**2 + Uncertainties["LeptonFullSimUncertainty_%s"%regionCombination]**2 + Uncertainties["pileupUncertainty_%s"%regionCombination]**2 + Uncertainties["ISRUncertainty_%s"%regionCombination]**2 + Uncertainties["BTagUncertaintyLight_%s"%regionCombination]**2  + Uncertainties["BTagUncertaintyHeavy_%s"%regionCombination]**2  + Uncertainties["TriggerEffUncertainty_%s"%regionCombination]**2 + LumiUncertainty**2)
						
						### Total uncertainty
						Uncertainties["TotalUncertainty_%s"%regionCombination] = sqrt(Uncertainties["SystUncertainty_%s"%regionCombination]**2 + Uncertainties["StatUncertainty_%s"%regionCombination]**2)		
											
						for uncertainty in uncertaintySources:
							uncertaintyArrays["%s_%s"%(uncertainty,regionCombination)].append(Uncertainties["%s_%s"%(uncertainty,regionCombination)])					
								
			m_n += stepsize		
			
		m_b += stepsize
				
	for regionCombination in regionCombinations:
		
		for uncertainty in uncertaintySources:
			if regionCombination == "lowNll":
				regionName = "ttbar_like"
			else:
				regionName = "non_ttbar_like"
			Graphs["%s_%s"%(uncertainty,regionCombination)]=TGraph2D("%s_%s"%(uncertainty,regionName),"%s_%s"%(uncertainty,regionCombination), len(uncertaintyArrays["%s_%s"%(uncertainty,regionCombination)]), array('d',masses_b), array('d',masses_n), array('d',uncertaintyArrays["%s_%s"%(uncertainty,regionCombination)]))
			Graphs["%s_%s"%(uncertainty,regionCombination)].SetNpx(nxbins)
			Graphs["%s_%s"%(uncertainty,regionCombination)].SetNpy(nybins)
			Histograms["%s_%s"%(uncertainty,regionCombination)] = Graphs["%s_%s"%(uncertainty,regionCombination)].GetHistogram()
			Histograms["%s_%s"%(uncertainty,regionCombination)].SetTitle(";m_{#tilde{b}} [GeV]; m_{#tilde{#chi_{2}^{0}}} [GeV]")
		
		region_label = "Signal Region"
			
		if "lowNll" in regionCombination:
			region_label_2 = "ttbar like"
		elif "highNll" in region:
			region_label_2 = "non-ttbar like"		


		plotPad.SetLogz()
		Histograms["Yield_%s"%regionCombination].SetZTitle("SF-OF yield")
		Histograms["Yield_%s"%regionCombination].Draw("colz")
		latexLumi.DrawLatex(0.85, 0.96, "%s fb^{-1} (13 TeV)"%(printLumi,))
		latexCMS.DrawLatex(0.19,0.89,"CMS")
		#~ latexCMSExtra.DrawLatex(0.19,0.85,"Private Work - Simulation")
		latexCMSExtra.DrawLatex(0.19,0.85,"Simulation")
		#~ latexCMSExtra.DrawLatex(0.19,0.85,"Unpublished")
		latex.DrawLatex(0.175, 0.75, "#splitline{Simplified Model}{#splitline{T6bbslepton, m(#tilde{#chi}_{1}^{0})=100 GeV}{#splitline{"+region_label+"}{"+region_label_2+"}}}")
		canv.Update()
		canv.Print("fig/Yields/T6bbllslepton_%s.pdf"%(regionCombination))
		
		plotPad.SetLogz(0)
		Histograms["StatUncertainty_%s"%regionCombination].SetZTitle("SF-OF rel. stat. uncertainty")
		Histograms["StatUncertainty_%s"%regionCombination].Draw("colz")
		latexLumi.DrawLatex(0.85, 0.96, "%s fb^{-1} (13 TeV)"%(printLumi,))
		latexCMS.DrawLatex(0.19,0.89,"CMS")
		#~ latexCMSExtra.DrawLatex(0.19,0.85,"Private Work - Simulation")
		latexCMSExtra.DrawLatex(0.19,0.85,"Simulation")
		#~ latexCMSExtra.DrawLatex(0.19,0.85,"Unpublished")
		latex.DrawLatex(0.175, 0.75, "#splitline{Simplified Model}{#splitline{T6bbslepton, m(#tilde{#chi}_{1}^{0})=100 GeV}{#splitline{"+region_label+"}{"+region_label_2+"}}}")
		canv.Update()
		canv.Print("fig/statUncertainties/T6bbllslepton_%s_stat_err.pdf"%(regionCombination))
		
		
		Histograms["SystUncertainty_%s"%regionCombination].SetZTitle("SF-OF rel. syst. uncertainty")
		Histograms["SystUncertainty_%s"%regionCombination].Draw("colz")
		latexLumi.DrawLatex(0.85, 0.96, "%s fb^{-1} (13 TeV)"%(printLumi,))
		latexCMS.DrawLatex(0.19,0.89,"CMS")
		#~ latexCMSExtra.DrawLatex(0.19,0.85,"Private Work - Simulation")
		latexCMSExtra.DrawLatex(0.19,0.85,"Simulation")
		#~ latexCMSExtra.DrawLatex(0.19,0.85,"Unpublished")
		latex.DrawLatex(0.175, 0.75, "#splitline{Simplified Model}{#splitline{T6bbslepton, m(#tilde{#chi}_{1}^{0})=100 GeV}{#splitline{"+region_label+"}{"+region_label_2+"}}}")
		canv.Update()
		canv.Print("fig/sysUncertainties/T6bbllslepton_%s_syst_err.pdf"%(regionCombination))
		
		
		Histograms["TotalUncertainty_%s"%regionCombination].SetZTitle("SF-OF total rel. uncertainty")
		Histograms["TotalUncertainty_%s"%regionCombination].Draw("colz")
		latexLumi.DrawLatex(0.85, 0.96, "%s fb^{-1} (13 TeV)"%(printLumi,))
		latexCMS.DrawLatex(0.19,0.89,"CMS")
		#~ latexCMSExtra.DrawLatex(0.19,0.85,"Private Work - Simulation")
		latexCMSExtra.DrawLatex(0.19,0.85,"Simulation")
		#~ latexCMSExtra.DrawLatex(0.19,0.85,"Unpublished")
		latex.DrawLatex(0.175, 0.75, "#splitline{Simplified Model}{#splitline{T6bbslepton, m(#tilde{#chi}_{1}^{0})=100 GeV}{#splitline{"+region_label+"}{"+region_label_2+"}}}")
		canv.Update()
		canv.Print("fig/totUncertainties/T6bbllslepton_%s_tot_err.pdf"%(regionCombination))
		
		Histograms["Efficiency_%s"%regionCombination].SetZTitle("acceptance #times efficiency")
		Histograms["Efficiency_%s"%regionCombination].Draw("colz")
		latexLumi.DrawLatex(0.85, 0.96, "%s fb^{-1} (13 TeV)"%(printLumi,))
		latexCMS.DrawLatex(0.19,0.89,"CMS")
		#~ latexCMSExtra.DrawLatex(0.19,0.85,"Private Work - Simulation")
		latexCMSExtra.DrawLatex(0.19,0.85,"Simulation")
		#~ latexCMSExtra.DrawLatex(0.19,0.85,"Unpublished")
		latex.DrawLatex(0.175, 0.75, "#splitline{Simplified Model}{#splitline{T6bbslepton, m(#tilde{#chi}_{1}^{0})=100 GeV}{#splitline{"+region_label+"}{"+region_label_2+"}}}")
		canv.Update()
		canv.Print("fig/Efficiencies/T6bbllslepton_%s_signalEfficiency.pdf"%(regionCombination))
		
		Histograms["EfficiencyUnscaled_%s"%regionCombination].SetZTitle("acceptance #times efficiency")
		Histograms["EfficiencyUnscaled_%s"%regionCombination].Draw("colz")
		latexLumi.DrawLatex(0.85, 0.96, "%s fb^{-1} (13 TeV)"%(printLumi,))
		latexCMS.DrawLatex(0.19,0.89,"CMS")
		#~ latexCMSExtra.DrawLatex(0.19,0.85,"Private Work - Simulation")
		latexCMSExtra.DrawLatex(0.19,0.85,"Simulation")
		#~ latexCMSExtra.DrawLatex(0.19,0.85,"Unpublished")
		latex.DrawLatex(0.175, 0.75, "#splitline{Simplified Model}{#splitline{T6bbslepton, m(#tilde{#chi}_{1}^{0})=100 GeV}{#splitline{"+region_label+"}{"+region_label_2+"}}}")
		canv.Update()
		canv.Print("fig/Efficiencies/T6bbllslepton_%s_signalEfficiencyUnscaled.pdf"%(regionCombination))
		
		
		
			
		for uncertainty in uncertaintySources:
			if not ( uncertainty == "Yield" or uncertainty == "StatUncertainty" or uncertainty == "SystUncertainty" or uncertainty == "TotalUncertainty" or uncertainty == "Efficiency" or uncertainty == "EfficiencyUnscaled"):			
				Histograms["%s_%s"%(uncertainty,regionCombination)].SetZTitle("%s"%uncertainty)
				Histograms["%s_%s"%(uncertainty,regionCombination)].Draw("colz")
				latexLumi.DrawLatex(0.85, 0.96, "%s fb^{-1} (13 TeV)"%(printLumi,))
				latexCMS.DrawLatex(0.19,0.89,"CMS")
				#~ latexCMSExtra.DrawLatex(0.19,0.85,"Private Work - Simulation")
				latexCMSExtra.DrawLatex(0.19,0.85,"Simulation")
				#~ latexCMSExtra.DrawLatex(0.19,0.85,"Unpublished")
				latex.DrawLatex(0.175, 0.75, "#splitline{Simplified Model}{#splitline{T6bbslepton, m(#tilde{#chi}_{1}^{0})=100 GeV}{#splitline{"+region_label+"}{"+region_label_2+"}}}")
				canv.Update()
				canv.Print("fig/%s/T6bbllslepton_%s_%s.pdf"%(uncertainty,regionCombination,uncertainty))
	
	histoFile = TFile("fig/T6bbllslepton_XSecUpperLimit_and_ExclusionContours.root")
	
	histoXSection = TH2D(histoFile.Get("XSecUpperLimit"))
	graphExpectedUpperLimit = TGraph2D(histoFile.Get("ExpectedUpperLimit"))
	graphExpectedUpperLimitUp = TGraph2D(histoFile.Get("ExpectedUpperLimitUp"))
	graphExpectedUpperLimitDown = TGraph2D(histoFile.Get("ExpectedUpperLimitDown"))
	graphExpectedUpperLimitUp2 = TGraph2D(histoFile.Get("ExpectedUpperLimitUp2"))
	graphExpectedUpperLimitDown2 = TGraph2D(histoFile.Get("ExpectedUpperLimitDown2"))
	graphObservedUpperLimit = TGraph2D(histoFile.Get("ObservedUpperLimit"))
	graphObservedUpperLimitUp = TGraph2D(histoFile.Get("ObservedUpperLimitUp"))
	graphObservedUpperLimitDown = TGraph2D(histoFile.Get("ObservedUpperLimitDown"))
		
	f1 = TFile("fig/SummaryFile.root","RECREATE")
	for regionCombination in regionCombinations:
		Histograms["Yield_%s"%regionCombination].Write()
		Histograms["Efficiency_%s"%regionCombination].Write()
		Histograms["StatUncertainty_%s"%regionCombination].Write()
		Histograms["SystUncertainty_%s"%regionCombination].Write()
		Histograms["TotalUncertainty_%s"%regionCombination].Write()
		
	histoXSection.Write()
	graphExpectedUpperLimit.Write()
	graphExpectedUpperLimitUp.Write()
	graphExpectedUpperLimitDown.Write()
	graphExpectedUpperLimitUp2.Write()
	graphExpectedUpperLimitDown2.Write()
	graphObservedUpperLimit.Write()
	graphObservedUpperLimitUp.Write()
	graphObservedUpperLimitDown.Write()	
	
	f1.Close()			
	
	
	

							
					
			
											
									
				

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
 
