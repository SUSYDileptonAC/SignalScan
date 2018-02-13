#! /usr/bin/env python

import sys
sys.path.append("/.automount/home/home__home4/institut_1b/schomakers/FrameWork/frameWorkBase")

from helpers import readTrees, totalNumberOfGeneratedEvents,createHistoFromTree
from locations import locations

from math import sqrt
from array import array

import ROOT

from sys import argv
import pickle	
from ROOT import TCanvas, TPad, TH1F, TH1I, THStack, TLegend, TF1, TH2F, TH2D, TFile, TMath
import ratios
from defs import sbottom_masses
from corrections import triggerEffs, rSFOFDirect

from math import sqrt
from setTDRStyle import setTDRStyle

from ConfigParser import ConfigParser
config = ConfigParser()
config.read("/.automount/home/home__home4/institut_1b/schomakers/FrameWork/SubmitScripts/Input/Master80X_MC.ini")
attic = []


ROOT.gStyle.SetOptStat(0)

sleptonPoints = {
				"1":{"m_b":900,"m_n":150},
				"2":{"m_b":900,"m_n":500},
				"3":{"m_b":1000,"m_n":300},
				"4":{"m_b":1200,"m_n":200},
				"5":{"m_b":1200,"m_n":1000},
				}

mllCuts = {
			"noMllCut":"1>0",
			"20To60":"mll < 60 && mll > 20",
			"60To86":"mll < 86 && mll > 60",
			"86To96":"mll < 96 && mll > 86",
			"96To150":"mll < 150 && mll > 96",
			"150To200":"mll < 200 && mll > 150",
			"200To300":"mll > 200 && mll < 300",
			"300To400":"mll > 300 && mll < 400",
			"Above400":"mll > 400 ",
			}

NllCuts = {
			"inclusiveNll":"1 > 0",
			"lowNll":"nLL < 21",
			"highNll":"nLL > 21",
			}

signalRegionCuts = {
			"basicCut":"1>0",
			"mll20Cut":"p4.M() > 20",
			"pt25Cut":"p4.Pt() > 25 && p4.M() > 20",
			"nJets2Cut":"nJets >= 2 && pt > 25 && mll > 20",
			"deltaPhiCut":"abs(deltaPhiJetMet1) > 0.4  && abs(deltaPhiJetMet2) > 0.4 && nJets >= 2 && pt > 25 && mll > 20",
			"met150Cut":"met > 150 && abs(deltaPhiJetMet1) > 0.4  && abs(deltaPhiJetMet2) > 0.4 && nJets >= 2 && pt > 25 && mll > 20",
			"signalRegionCut":"met > 150 && MT2 > 80 && abs(deltaPhiJetMet1) > 0.4  && abs(deltaPhiJetMet2) > 0.4 && nJets >= 2 && pt > 25 && mll > 20",
			}			
			

def signalYields(tree,cuts,scalingLumi):
	histo = createHistoFromTree(tree, "bTagWeight", cuts, 40, 0, 40)
	histo.Scale(scalingLumi)
	yields = float(histo.Integral(0,-1))
	histo.Delete()
	return yields
			
def produceSignalEfficiency(tree,cuts,weight):
	histo = createHistoFromTree(tree, "bTagWeight", cuts, 40, 0, 40)
	result = histo.Integral(0,-1)/weight
	histo.Delete()
	return result
			
def MCEventYield(tree,cuts):
	histo = createHistoFromTree(tree, "bTagWeight", cuts, 40, 0, 40)
	result = histo.Integral(0,-1)
	histo.Delete()
	return result	

def producePileupUncertainty(tree,cuts):
	
	#~ cuts = cuts.replace("chargeProduct < 0","chargeProduct < 0 && nGenVertices >= 20")
	cuts = cuts.replace("chargeProduct < 0","chargeProduct < 0 && nVertices > 16")
	histoUp = createHistoFromTree(tree, "bTagWeight", cuts, 40, 0, 40)
	yields = [float(histoUp.Integral(0,-1))]
	histoUp.Delete()
	
	#~ cuts = cuts.replace("nGenVertices >= 20","nGenVertices < 20")
	cuts = cuts.replace("nVertices > 16","nVertices <= 16")
	histoDown = createHistoFromTree(tree, "bTagWeight", cuts, 40, 0, 40)
	yields.append(float(histoDown.Integral(0,-1)))
	histoDown.Delete()

	return yields[0],yields[1]	
	
def produceJESUncertainty(tree,cuts):

	cuts = cuts.replace("nJets > 1","nShiftedJetsJESUp > 1")	
	cuts = cuts.replace("met","metJESUp")	
	cuts = cuts.replace("nLL","nLLJESUp")
	cuts = cuts.replace("MT2","MT2JESUp")
	histoUp = createHistoFromTree(tree, "bTagWeight", cuts, 40, 0, 40)
	yields = [float(histoUp.Integral(0,-1))]
	histoUp.Delete()
	
	cuts = cuts.replace("nShiftedJetsJESUp","nShiftedJetsJESDown")	
	cuts = cuts.replace("metJESUp","metJESDown")	
	cuts = cuts.replace("MT2JESUp","MT2JESDown")	
	histoDown = createHistoFromTree(tree, "bTagWeight", cuts, 40, 0, 40)
	yields.append(float(histoDown.Integral(0,-1)))
	histoDown.Delete()
	
	return yields[0],yields[1]

			
def produceISRUncertainty(tree,cuts):
	
	cuts = cuts.replace("ISRCorrection","(ISRCorrection+ISRUncertainty)")
	histoUp = createHistoFromTree(tree, "bTagWeight", cuts, 40, 0, 40)
	yields = [float(histoUp.Integral(0,-1))]
	histoUp.Delete()

	cuts = cuts.replace("(ISRCorrection+ISRUncertainty)","(ISRCorrection-ISRUncertainty)")
	histoDown = createHistoFromTree(tree, "bTagWeight", cuts, 40, 0, 40)
	yields.append(float(histoDown.Integral(0,-1)))
	histoDown.Delete()
	
	return yields[0],yields[1]
			
def produceLeptonFastSimUncertainty(tree,cuts):
	
	result = 0.
	cuts = cuts.replace("leptonFastSimScaleFactor1*leptonFastSimScaleFactor2*","")
	histo = createHistoFromTree(tree, "bTagWeight", cuts, 40, 0, 40)
	result = histo.Integral(0,-1)
	histo.Delete()
	
	return result
	
def produceLeptonFullSimUncertainty(tree,cuts):
	
	result = 0.
	#~ cuts = cuts.replace("leptonFullSimScaleFactor1*leptonFullSimScaleFactor2*","")
	cuts = cuts.replace("leptonFullSimScaleFactor1*","leptonFullSimScaleFactor1*(1+leptonFullSimScaleFactorErr1)*")
	cuts = cuts.replace("leptonFullSimScaleFactor2*","leptonFullSimScaleFactor2*(1+leptonFullSimScaleFactorErr2)*")
	histo = createHistoFromTree(tree, "bTagWeight", cuts, 40, 0, 40)
	result = histo.Integral(0,-1)
	histo.Delete()
	
	return result
	
def produceLeptonFastSimUncertainty(tree,cuts):
	
	result = 0.
	#~ cuts = cuts.replace("leptonFastSimScaleFactor1*leptonFastSimScaleFactor2*","")
	cuts = cuts.replace("leptonFastSimScaleFactor1*","leptonFastSimScaleFactor1*(1+leptonFastSimScaleFactorErr1)*")
	cuts = cuts.replace("leptonFastSimScaleFactor2*","leptonFastSimScaleFactor2*(1+leptonFastSimScaleFactorErr2)*")
	histo = createHistoFromTree(tree, "bTagWeight", cuts, 40, 0, 40)
	result = histo.Integral(0,-1)
	histo.Delete()
	
	return result
	
def produceBTagUncertainty(tree,cuts):
	
	cuts = cuts.replace("bTagWeight","bTagWeight * (1 + bTagWeightErrHeavy)")
	histoUp = createHistoFromTree(tree, "bTagWeight", cuts, 40, 0, 40)
	yields = [float(histoUp.Integral(0,-1))]
	histoUp.Delete()
	
	cuts = cuts.replace("bTagWeight * (1 + bTagWeightErrHeavy)","bTagWeight * (1 + bTagWeightErrLight)")
	histoDown = createHistoFromTree(tree, "bTagWeight", cuts, 40, 0, 40)
	yields.append(float(histoDown.Integral(0,-1)))
	histoDown.Delete()
	
	return yields[0],yields[1]
	
	
def produceScaleShiftUncertainty(tree,cuts,weightIndex):
	
	
	cuts = cuts.replace("bTagWeight","bTagWeight * scaleWeight%s"%str(weightIndex+1))
	histo = createHistoFromTree(tree, "bTagWeight", cuts, 40, 0, 40)
	result = float(histo.Integral(0,-1))
	histo.Delete()	
	
	return result
	

		
if (__name__ == "__main__"):
	import ratios
	
	import pickle
	
	signalDenominatorFile = TFile("/.automount/home/home__home4/institut_1b/schomakers/FrameWork/SignalScan/T6bbllsleptonDenominatorHisto7.root")
	denominatorHisto = TH2F(signalDenominatorFile.Get("massScan"))
	denominatorHistoLowPU = TH2F(signalDenominatorFile.Get("massScanLowPU"))
	denominatorHistoHighPU = TH2F(signalDenominatorFile.Get("massScanHighPU"))
	ISRNormalizationHisto = TH2F(signalDenominatorFile.Get("ISRNormalization"))
	ISRNormalizationHistoUp = TH2F(signalDenominatorFile.Get("ISRNormalizationUp"))
	ISRNormalizationHistoDown = TH2F(signalDenominatorFile.Get("ISRNormalizationDown"))
	
	scaleNormalizationHistos = {}
	for i in range(1,9):
		scaleNormalizationHistos["ScaleWeight"+str(i)] = TH2F(signalDenominatorFile.Get("ScaleWeight"+str(i)))

	lumi = 35867.
	
	### more inclusive signal regions, do not yet divide in mll and 
	### likelihood regions
	
	#~ path = locations.dataSetPath	
	#~ CutRegions = ["basicCut","mll20Cut","pt25Cut"]
	#~ mllCutRegions = ["noMllCut"]
	#~ NllCutRegions = ["inclusiveNll"]
	
	### closer to signal region -> can be calculated with NLL trees
	### split in mll and likelihood bins
		
	path = locations.dataSetPathSignalNLL	
	CutRegions = ["met150Cut","nJets2Cut","signalRegionCut","deltaPhiCut"]
	mllCutRegions = ["noMllCut","20To60","60To86","86To96","96To150","150To200","200To300","300To400","Above400"]
	NllCutRegions = ["inclusiveNll","lowNll","highNll"]
	
	points = ["1","2","3","4","5"]
	
	### Factor to calculate the expected number of events with at least 2 leptons
	BR_factor =  0.25      + 0.5   + 0.25*(0.13+0.06*0.124+0.0011*(0.015+0.113+0.312))
	
	triggerEffUncertainty = 0.03
	lumiUncertainty = 0.026
	
	EETrees = readTrees(path, "EE")	
	MMTrees = readTrees(path, "MuMu")	
	
	EETriggerEff = triggerEffs.inclusive.effEE.val
	EMuTriggerEff = triggerEffs.inclusive.effEM.val
	MuMuTriggerEff = triggerEffs.inclusive.effMM.val
	RSFOF = rSFOFDirect.inclusive.val

	## loop over example points
	for point in points:
		m_b = sleptonPoints[point]["m_b"]
		m_sbottom = str(m_b)
		M_SBOTTOM = "m_b_"+str(m_b)
		m_n_2 =sleptonPoints[point]["m_n"]
		m_neutralino_2 = str(m_n_2)
		
		xsection = getattr(sbottom_masses, M_SBOTTOM).cross_section13TeV
		
		sampleName = "T6bbllslepton_msbottom_%s_mneutralino_%s"%(m_sbottom,m_neutralino_2)
		fileName = "T6bbllslepton_msbottom_%s_mneutralino_%s"%(m_sbottom,m_neutralino_2)
		
		denominator = denominatorHisto.GetBinContent(denominatorHisto.GetXaxis().FindBin(int(sampleName.split("_")[2])),denominatorHisto.GetYaxis().FindBin(int(sampleName.split("_")[4])))
		denominatorLowPU = denominatorHistoLowPU.GetBinContent(denominatorHistoLowPU.GetXaxis().FindBin(int(sampleName.split("_")[2])),denominatorHistoLowPU.GetYaxis().FindBin(int(sampleName.split("_")[4])))
		denominatorHighPU = denominatorHistoHighPU.GetBinContent(denominatorHistoHighPU.GetXaxis().FindBin(int(sampleName.split("_")[2])),denominatorHistoHighPU.GetYaxis().FindBin(int(sampleName.split("_")[4])))
					
		ISRNormalization = ISRNormalizationHisto.GetBinContent(ISRNormalizationHisto.GetXaxis().FindBin(int(sampleName.split("_")[2])),denominatorHisto.GetYaxis().FindBin(int(sampleName.split("_")[4])))
		ISRNormalizationUp = ISRNormalizationHistoUp.GetBinContent(ISRNormalizationHisto.GetXaxis().FindBin(int(sampleName.split("_")[2])),denominatorHisto.GetYaxis().FindBin(int(sampleName.split("_")[4])))
		ISRNormalizationDown = ISRNormalizationHistoDown.GetBinContent(ISRNormalizationHisto.GetXaxis().FindBin(int(sampleName.split("_")[2])),denominatorHisto.GetYaxis().FindBin(int(sampleName.split("_")[4])))
		
		scaleNormalizations = []
		
		for weightIndex in range(1,9):
			scaleNormalizations.append(scaleNormalizationHistos["ScaleWeight"+str(weightIndex)].GetBinContent(scaleNormalizationHistos["ScaleWeight"+str(weightIndex)].GetXaxis().FindBin(int(sampleName.split("_")[2])),scaleNormalizationHistos["ScaleWeight"+str(weightIndex)].GetYaxis().FindBin(int(sampleName.split("_")[4]))))
			
		totalEventYield = lumi*xsection*BR_factor
		totalEventYieldUncert = lumi*xsection*BR_factor/sqrt(denominator)
		print totalEventYield
		print totalEventYieldUncert
		
		scalingLumi = lumi*xsection/(denominator*ISRNormalization)
		weight = denominator*ISRNormalization
	
		### loop over regions
		### calculate yields and uncertainties
		
		for region in CutRegions:
			signalRegionCut = signalRegionCuts[region]
			 
			for mllCutRegion in mllCutRegions:
				mllCut = mllCuts[mllCutRegion]
				
				for NllCutRegion in NllCutRegions:
					NllCut = NllCuts[NllCutRegion]
			
			
					cuts = "bTagWeight*ISRCorrection*leptonFullSimScaleFactor1*leptonFullSimScaleFactor2*leptonFastSimScaleFactor1*leptonFastSimScaleFactor2*(met / caloMet < 5 && nBadMuonJets == 0 && nUnmatchedJets == 0 && deltaR > 0.1 && chargeProduct < 0 && ((pt1 > 25 && pt2 > 20) || (pt1 > 20 && pt2 > 25)) && abs(eta1) < 2.4 && abs(eta2) < 2.4 && %s && %s && %s)"%(signalRegionCut,mllCut,NllCut)
					unweightedCuts = "(met / caloMet < 5 && nBadMuonJets == 0 && nUnmatchedJets == 0 && chargeProduct < 0 && deltaR > 0.1 && ((pt1 > 25 && pt2 > 20) || (pt1 > 20 && pt2 > 25)) && abs(eta1) < 2.4 && abs(eta2) < 2.4 && %s && %s && %s)"%(signalRegionCut,mllCut,NllCut)	
						
					
					for sample, tree in EETrees.iteritems():
						if sample == sampleName:
							EEsignalYieldMet = signalYields(tree,cuts, EETriggerEff * scalingLumi)
							EEsignalEfficiencyMet = produceSignalEfficiency(tree, cuts, weight)
							EEMCEventsMet = MCEventYield(tree,unweightedCuts)
							EEISRUpMet,EEISRDownMet = produceISRUncertainty(tree,cuts)
							EEJESUpMet,EEJESDownMet= produceJESUncertainty(tree,cuts)
							EELeptonFastSimScaleFactorShiftedMet = produceLeptonFastSimUncertainty(tree,cuts)
							EELeptonFullSimScaleFactorShiftedMet = produceLeptonFullSimUncertainty(tree,cuts)
							EEPileupHighMet,EEPileupLowMet = producePileupUncertainty(tree,cuts)
							EEBTagHeavyMet,EEBTagLightMet = produceBTagUncertainty(tree,cuts)
							
							
							EEScaleShiftedMet = []
							for weightIndex in range(0,8):
								EEScaleShiftedMet.append(produceScaleShiftUncertainty(tree,cuts,weightIndex)/scaleNormalizations[weightIndex])
							
							if denominator > 0:
								EEISRUpMet = EEISRUpMet*ISRNormalization/ISRNormalizationUp
								EEISRDownMet = EEISRDownMet*ISRNormalization/ISRNormalizationDown
								
								EEPileupHighMet = EEPileupHighMet * denominator/denominatorHighPU
								EEPileupLowMet = EEPileupLowMet * denominator/denominatorLowPU
							
							
							cuts = cuts.replace("met","genMet")
							
							EEsignalYieldGenMet = signalYields(tree,cuts, EETriggerEff * scalingLumi)
							EEsignalEfficiencyGenMet = produceSignalEfficiency(tree, cuts, weight)
							EEMCEventsGenMet = MCEventYield(tree,unweightedCuts)
							EEISRUpGenMet,EEISRDownGenMet = produceISRUncertainty(tree,cuts)
							EEJESUpGenMet,EEJESDownGenMet= produceJESUncertainty(tree,cuts)
							EELeptonFastSimScaleFactorShiftedGenMet = produceLeptonFastSimUncertainty(tree,cuts)
							EELeptonFullSimScaleFactorShiftedGenMet = produceLeptonFullSimUncertainty(tree,cuts)
							EEPileupHighGenMet,EEPileupLowGenMet = producePileupUncertainty(tree,cuts)
							EEBTagHeavyGenMet,EEBTagLightGenMet = produceBTagUncertainty(tree,cuts)
							
							EEScaleShiftedGenMet = []
							for weightIndex in range(0,8):
								EEScaleShiftedGenMet.append(produceScaleShiftUncertainty(tree,cuts,weightIndex)/scaleNormalizations[weightIndex])
							
							if denominator > 0:
								EEISRUpGenMet = EEISRUpGenMet*ISRNormalization/ISRNormalizationUp
								EEISRDownGenMet = EEISRDownGenMet*ISRNormalization/ISRNormalizationDown
								
								EEPileupHighGenMet = EEPileupHighGenMet * denominator/denominatorHighPU
								EEPileupLowGenMet = EEPileupLowGenMet * denominator/denominatorLowPU						
							
							cuts = cuts.replace("genMet","met")
							
							EEsignalYield = 0.5* (EEsignalYieldMet + EEsignalYieldGenMet)
							EEsignalEfficiency = 0.5* (EEsignalEfficiencyMet + EEsignalEfficiencyGenMet)
							EEMCEvents = 0.5* (EEMCEventsMet + EEMCEventsGenMet)
							EEISRUp = 0.5* (EEISRUpMet + EEISRUpGenMet) 
							EEISRDown = 0.5* (EEISRDownMet + EEISRDownGenMet) 
							EEJESUp = 0.5* (EEJESUpMet + EEJESUpGenMet) 
							EEJESDown = 0.5* (EEJESDownMet + EEJESDownGenMet) 
							EELeptonFastSimScaleFactorShifted = 0.5* (EELeptonFastSimScaleFactorShiftedMet + EELeptonFastSimScaleFactorShiftedGenMet) 
							EELeptonFullSimScaleFactorShifted = 0.5* (EELeptonFullSimScaleFactorShiftedMet + EELeptonFullSimScaleFactorShiftedGenMet) 
							EEPileupHigh = 0.5* (EEPileupHighMet + EEPileupHighGenMet) 
							EEPileupLow = 0.5* (EEPileupLowMet + EEPileupLowGenMet) 
							EEBTagHeavy = 0.5* (EEBTagHeavyMet + EEBTagHeavyGenMet) 
							EEBTagLight = 0.5* (EEBTagLightMet + EEBTagLightGenMet)						
							
							
							EEScaleShifted = []
							for weightIndex in range(0,8):
								EEScaleShifted.append(0.5*(EEScaleShiftedMet[weightIndex]+EEScaleShiftedGenMet[weightIndex]))
							
							UnscaledYield = EEsignalYield / (EETriggerEff * scalingLumi)
							UnscaledYieldEE = UnscaledYield
							
							if UnscaledYield > 0:
								EEISRUncertainty = max(abs(EEISRUp-UnscaledYield)/UnscaledYield,abs(EEISRDown-UnscaledYield)/UnscaledYield) 
								EEJESUncertainty = max(abs(EEJESUp-UnscaledYield)/UnscaledYield,abs(EEJESDown-UnscaledYield)/UnscaledYield) 
								EEPileupUncertainty = max(abs(EEPileupHigh-UnscaledYield)/UnscaledYield,abs(EEPileupLow-UnscaledYield)/UnscaledYield) 
								EELeptonFullSimUncertainty = abs(EELeptonFullSimScaleFactorShifted-UnscaledYield)/UnscaledYield
								EELeptonFastSimUncertainty = abs(EELeptonFastSimScaleFactorShifted-UnscaledYield)/UnscaledYield
								EEBTagHeavyUncertainty = abs(EEBTagHeavy-UnscaledYield)/UnscaledYield
								EEBTagLightUncertainty = abs(EEBTagLight-UnscaledYield)/UnscaledYield
								EEMetUncertainty = abs(EEsignalYieldMet-EEsignalYield)/EEsignalYield
								EEStatUncertainty = sqrt(EEMCEvents)
								
								EEScaleUncertainty = 0
								for weightIndex in range(0,8):
									EEScaleUncertainty = max (EEScaleUncertainty, abs(EEScaleShifted[weightIndex]-UnscaledYield)/UnscaledYield)
							else:
								EEISRUncertainty = 0.
								EEJESUncertainty = 0. 
								EEPileupUncertainty = 0.
								EELeptonFullSimUncertainty = 0.
								EELeptonFastSimUncertainty = 0.
								EEBTagHeavyUncertainty = 0.
								EEBTagLightUncertainty = 0.
								EEMetUncertainty = 0.
								EEStatUncertainty = 0.
								EEScaleUncertainty = 0.
							
							
					
					for sample, tree in MMTrees.iteritems():
						if sample == sampleName:
							MuMusignalYieldMet = signalYields(tree,cuts, MuMuTriggerEff * scalingLumi)
							MuMusignalEfficiencyMet = produceSignalEfficiency(tree, cuts, weight)
							MuMuMCEventsMet = MCEventYield(tree,unweightedCuts)
							MuMuISRUpMet,MuMuISRDownMet = produceISRUncertainty(tree,cuts)
							MuMuJESUpMet,MuMuJESDownMet= produceJESUncertainty(tree,cuts)
							MuMuLeptonFastSimScaleFactorShiftedMet = produceLeptonFastSimUncertainty(tree,cuts)
							MuMuLeptonFullSimScaleFactorShiftedMet = produceLeptonFullSimUncertainty(tree,cuts)
							MuMuPileupHighMet,MuMuPileupLowMet = producePileupUncertainty(tree,cuts)
							MuMuBTagHeavyMet,MuMuBTagLightMet = produceBTagUncertainty(tree,cuts)
							
							MuMuScaleShiftedMet = []
							for weightIndex in range(0,8):
								MuMuScaleShiftedMet.append(produceScaleShiftUncertainty(tree,cuts,weightIndex)/scaleNormalizations[weightIndex])
							
							
							if ISRNormalization > 0:
								MuMuISRUpMet = MuMuISRUpMet*ISRNormalization/ISRNormalizationUp
								MuMuISRDownMet = MuMuISRDownMet*ISRNormalization/ISRNormalizationDown
								
								MuMuPileupHighMet = MuMuPileupHighMet * denominator/denominatorHighPU
								MuMuPileupLowMet = MuMuPileupLowMet * denominator/denominatorLowPU
							
							
							cuts = cuts.replace("met","genMet")
							
							MuMusignalYieldGenMet = signalYields(tree,cuts, MuMuTriggerEff * scalingLumi)
							MuMusignalEfficiencyGenMet = produceSignalEfficiency(tree, cuts, weight)
							MuMuMCEventsGenMet = MCEventYield(tree,unweightedCuts)
							MuMuISRUpGenMet,MuMuISRDownGenMet = produceISRUncertainty(tree,cuts)
							MuMuJESUpGenMet,MuMuJESDownGenMet= produceJESUncertainty(tree,cuts)
							MuMuLeptonFastSimScaleFactorShiftedGenMet = produceLeptonFastSimUncertainty(tree,cuts)
							MuMuLeptonFullSimScaleFactorShiftedGenMet = produceLeptonFullSimUncertainty(tree,cuts)
							MuMuPileupHighGenMet,MuMuPileupLowGenMet = producePileupUncertainty(tree,cuts)
							MuMuBTagHeavyGenMet,MuMuBTagLightGenMet = produceBTagUncertainty(tree,cuts)
							
							MuMuScaleShiftedGenMet = []
							for weightIndex in range(0,8):
								MuMuScaleShiftedGenMet.append(produceScaleShiftUncertainty(tree,cuts,weightIndex)/scaleNormalizations[weightIndex])
							
							
							if ISRNormalization > 0:
								MuMuISRUpGenMet = MuMuISRUpGenMet*ISRNormalization/ISRNormalizationUp
								MuMuISRDownGenMet = MuMuISRDownGenMet*ISRNormalization/ISRNormalizationDown
								
								MuMuPileupHighGenMet = MuMuPileupHighGenMet * denominator/denominatorHighPU
								MuMuPileupLowGenMet = MuMuPileupLowGenMet * denominator/denominatorLowPU
							
							cuts = cuts.replace("genMet","met")
							
							MuMusignalYield = 0.5* (MuMusignalYieldMet + MuMusignalYieldGenMet)
							MuMusignalEfficiency = 0.5* (MuMusignalEfficiencyMet + MuMusignalEfficiencyGenMet)
							MuMuMCEvents = 0.5* (MuMuMCEventsMet + MuMuMCEventsGenMet)
							MuMuISRUp = 0.5* (MuMuISRUpMet + MuMuISRUpGenMet) 
							MuMuISRDown = 0.5* (MuMuISRDownMet + MuMuISRDownGenMet) 
							MuMuJESUp = 0.5* (MuMuJESUpMet + MuMuJESUpGenMet) 
							MuMuJESDown = 0.5* (MuMuJESDownMet + MuMuJESDownGenMet) 
							MuMuLeptonFastSimScaleFactorShifted = 0.5* (MuMuLeptonFastSimScaleFactorShiftedMet + MuMuLeptonFastSimScaleFactorShiftedGenMet) 
							MuMuLeptonFullSimScaleFactorShifted = 0.5* (MuMuLeptonFullSimScaleFactorShiftedMet + MuMuLeptonFullSimScaleFactorShiftedGenMet) 
							MuMuPileupHigh = 0.5* (MuMuPileupHighMet + MuMuPileupHighGenMet) 
							MuMuPileupLow = 0.5* (MuMuPileupLowMet + MuMuPileupLowGenMet) 
							MuMuBTagHeavy = 0.5* (MuMuBTagHeavyMet + MuMuBTagHeavyGenMet) 
							MuMuBTagLight = 0.5* (MuMuBTagLightMet + MuMuBTagLightGenMet)
							
							MuMuScaleShifted = []
							for weightIndex in range(0,8):
								MuMuScaleShifted.append(0.5*(MuMuScaleShiftedMet[weightIndex]+MuMuScaleShiftedGenMet[weightIndex]))
							
							UnscaledYield = MuMusignalYield / (MuMuTriggerEff * scalingLumi)
							UnscaledYieldMuMu = UnscaledYield
							
							if UnscaledYield > 0:
								MuMuISRUncertainty = max(abs(MuMuISRUp-UnscaledYield)/UnscaledYield,abs(MuMuISRDown-UnscaledYield)/UnscaledYield) 
								MuMuJESUncertainty = max(abs(MuMuJESUp-UnscaledYield)/UnscaledYield,abs(MuMuJESDown-UnscaledYield)/UnscaledYield) 
								MuMuPileupUncertainty = max(abs(MuMuPileupHigh-UnscaledYield)/UnscaledYield,abs(MuMuPileupLow-UnscaledYield)/UnscaledYield) 
								MuMuLeptonFullSimUncertainty = abs(MuMuLeptonFullSimScaleFactorShifted-UnscaledYield)/UnscaledYield
								MuMuLeptonFastSimUncertainty = abs(MuMuLeptonFastSimScaleFactorShifted-UnscaledYield)/UnscaledYield
								MuMuBTagHeavyUncertainty = abs(MuMuBTagHeavy-UnscaledYield)/UnscaledYield
								MuMuBTagLightUncertainty = abs(MuMuBTagLight-UnscaledYield)/UnscaledYield
								MuMuMetUncertainty = abs(MuMusignalYieldMet-MuMusignalYield)/MuMusignalYield
								MuMuStatUncertainty = sqrt(MuMuMCEvents)
								
								MuMuScaleUncertainty = 0
								for weightIndex in range(0,8):
									MuMuScaleUncertainty = max (MuMuScaleUncertainty, abs(MuMuScaleShifted[weightIndex]-UnscaledYield)/UnscaledYield)
									
							else:
								MuMuISRUncertainty = 0.
								MuMuJESUncertainty = 0. 
								MuMuPileupUncertainty = 0.
								MuMuLeptonFullSimUncertainty = 0.
								MuMuLeptonFastSimUncertainty = 0.
								MuMuBTagHeavyUncertainty = 0.
								MuMuBTagLightUncertainty = 0.
								MuMuMetUncertainty = 0.
								MuMuStatUncertainty = 0.
								MuMuScaleUncertainty = 0.
							
							
							
					counts = {}
					SFYield = EEsignalYield + MuMusignalYield
					if SFYield > 0:
						if EEMCEvents+MuMuMCEvents > 0:
							StatUncertainty = sqrt(EEMCEvents+MuMuMCEvents)/(EEMCEvents+MuMuMCEvents)
						else:
							StatUncertainty = 0
						
						Mean = UnscaledYieldEE + UnscaledYieldMuMu
						### JES Uncertainty
						JESUp = EEJESUp + MuMuJESUp
						JESDown = EEJESDown + MuMuJESDown
						
						if Mean > 0:
							JESUncertainty = max(abs(JESUp-Mean)/Mean,abs(JESDown-Mean)/Mean)
						else:
							JESUncertainty = 0
						
						METYield = 	EEsignalYieldMet + MuMusignalYieldMet
						GenMETYield = 	EEsignalYieldGenMet + MuMusignalYieldGenMet
						METUncertainty = 0.5*abs(METYield-GenMETYield)/SFYield
												
						
						### Lepton FastSim Uncertainty						
						LeptonFastSimScaleFactorShifted = EELeptonFastSimScaleFactorShifted + MuMuLeptonFastSimScaleFactorShifted
						
						if Mean > 0:
							LeptonFastSimUncertainty = abs(LeptonFastSimScaleFactorShifted-Mean)/Mean
						else:
							LeptonFastSimUncertainty = 0
						
						### Lepton FullSim Uncertainty						
						LeptonFastSimScaleFactorShifted = EELeptonFastSimScaleFactorShifted + MuMuLeptonFastSimScaleFactorShifted
						
						if Mean > 0:
							LeptonFullSimUncertainty = abs(LeptonFastSimScaleFactorShifted-Mean)/Mean
						else:
							LeptonFullSimUncertainty = 0
						
						###  Pileup Uncertainty
						
						PileupHigh = EEPileupHigh + MuMuPileupHigh
						PileupLow = EEPileupLow + MuMuPileupLow
						
						if Mean > 0:
							PileupUncertainty = max(abs(PileupHigh-Mean)/Mean,abs(PileupLow-Mean)/Mean)
						else:
							PileupUncertainty = 0
						
						### ISR Uncertainty
						ISRUp = EEISRUp + MuMuISRUp
						ISRDown = EEISRDown + MuMuISRDown
						
						if Mean > 0:
							ISRUncertainty = max(abs(ISRUp-Mean)/Mean,abs(ISRDown-Mean)/Mean)
						else:
							ISRUncertainty = 0
							
						### BTag Uncertainty
						BTagHeavy = EEBTagHeavy + MuMuBTagHeavy
						
						BTagLight = EEBTagLight + MuMuBTagLight
						
						if Mean > 0:
							BTagUncertaintyHeavy = abs(BTagHeavy-Mean)/Mean
							BTagUncertaintyLight = abs(BTagLight-Mean)/Mean
						else:
							BTagUncertaintyHeavy = 0
							BTagUncertaintyLight = 0
							
						ScaleUncertainty = 0
						if Mean > 0:
							for weightIndex in range(0,8):
								ScaleUncertainty = max (ScaleUncertainty, abs(MuMuScaleShifted[weightIndex]+EEScaleShifted[weightIndex]-Mean)/Mean)
							
						SystUncertainty = sqrt(ScaleUncertainty**2 + JESUncertainty**2 + METUncertainty**2 + LeptonFastSimUncertainty**2 + LeptonFullSimUncertainty**2 + PileupUncertainty**2 + ISRUncertainty**2 + BTagUncertaintyLight**2 + BTagUncertaintyHeavy**2  + triggerEffUncertainty**2 + lumiUncertainty**2)
						
					else:
						SFYield = 0 
						StatUncertainty = 0
						SystUncertainty = 0
				
					counts["Signal"] = {"totalEventYield":totalEventYield,"totalEventYieldUncert":totalEventYieldUncert,"SFYield":SFYield,"StatUncertainty":StatUncertainty,"SystUncertainty":SystUncertainty}		
					outFilePkl = open("shelvesCutFlowTables/%s_%s_%s_%s.pkl"%(fileName,region,mllCutRegion,NllCutRegion),"w")
					pickle.dump(counts, outFilePkl)
					outFilePkl.close()
										
		
