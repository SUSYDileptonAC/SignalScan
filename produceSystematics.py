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


ROOT.gStyle.SetOptStat(0)

massCuts = {
			"20To60":"mll < 60 && mll > 20",
			"60To86":"mll < 86 && mll > 60",
			"86To96":"mll < 96 && mll > 86",
			"96To150":"mll < 150 && mll > 96",
			"150To200":"mll < 200 && mll > 150",
			"200To300":"mll > 200 && mll < 300",
			"300To400":"mll > 300 && mll < 400",
			"Above400":"mll > 400 ",
			}
			
nLLCuts = {
			"default":"nLL > 0",
			"lowNll":"nLL < 21",
			"highNll":"nLL > 21",
			}
			
HTCuts = {
			"lowHT":"ht < 400",
			"mediumHT":"ht > 400 && ht < 800",
			"highHT":"ht > 800",
			}
			
MT2Cuts = {
			"lowMT2":"MT2 < 80",
			"highMT2":"MT2 > 80",
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
	
			
			   
def plot():
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
	

	path = locations.dataSetPathSignalNLL
	
	print path
	
	lumi = 35867.
		
	#~ m_b_min = 700
	m_b_min = 1000
	#~ m_b_min = 1150
	m_b_max = 1100
	#~ m_b_max = 1600
	m_n_min = 150

	
	EETrees = readTrees(path, "EE")
	EMTrees = readTrees(path, "EMu")
	MMTrees = readTrees(path, "MuMu")
		
	massRegions = ["20To60","60To86","86To96","96To150","150To200","200To300","300To400","Above400"]
	nLLRegions = ["lowNll","highNll"]
	MT2Regions = ["highMT2"]
	
	signalBins = []
	signalCuts = {}
	
	EETriggerEff = triggerEffs.inclusive.effEE.val
	EMuTriggerEff = triggerEffs.inclusive.effEM.val
	MuMuTriggerEff = triggerEffs.inclusive.effMM.val
	RSFOF = rSFOFDirect.inclusive.val
	
	triggerEffUncertainty = 0.03
	lumiUncertainty = 0.026

	for massRegion in massRegions:
		for nLLRegion in nLLRegions:
			for MT2Region in MT2Regions:
				signalBins.append("%s_%s_%s"%(massRegion,nLLRegion,MT2Region))
				signalCuts["%s_%s_%s"%(massRegion,nLLRegion,MT2Region)] = "%s && %s && %s"%(massCuts[massRegion],nLLCuts[nLLRegion],MT2Cuts[MT2Region])
	
	
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
			
			#~ print m_n
			m_neutralino_2 = str(m_n)
			counts = {}
			
			sampleName = "T6bbllslepton_msbottom_%s_mneutralino_%s"%(m_sbottom,m_neutralino_2)
			
			
			denominator = denominatorHisto.GetBinContent(denominatorHisto.GetXaxis().FindBin(int(sampleName.split("_")[2])),denominatorHisto.GetYaxis().FindBin(int(sampleName.split("_")[4])))
			denominatorLowPU = denominatorHistoLowPU.GetBinContent(denominatorHistoLowPU.GetXaxis().FindBin(int(sampleName.split("_")[2])),denominatorHistoLowPU.GetYaxis().FindBin(int(sampleName.split("_")[4])))
			denominatorHighPU = denominatorHistoHighPU.GetBinContent(denominatorHistoHighPU.GetXaxis().FindBin(int(sampleName.split("_")[2])),denominatorHistoHighPU.GetYaxis().FindBin(int(sampleName.split("_")[4])))
						
			ISRNormalization = ISRNormalizationHisto.GetBinContent(ISRNormalizationHisto.GetXaxis().FindBin(int(sampleName.split("_")[2])),denominatorHisto.GetYaxis().FindBin(int(sampleName.split("_")[4])))
			ISRNormalizationUp = ISRNormalizationHistoUp.GetBinContent(ISRNormalizationHisto.GetXaxis().FindBin(int(sampleName.split("_")[2])),denominatorHisto.GetYaxis().FindBin(int(sampleName.split("_")[4])))
			ISRNormalizationDown = ISRNormalizationHistoDown.GetBinContent(ISRNormalizationHisto.GetXaxis().FindBin(int(sampleName.split("_")[2])),denominatorHisto.GetYaxis().FindBin(int(sampleName.split("_")[4])))
			
			scaleNormalizations = []
			
			for weightIndex in range(1,9):
				scaleNormalizations.append(scaleNormalizationHistos["ScaleWeight"+str(weightIndex)].GetBinContent(scaleNormalizationHistos["ScaleWeight"+str(weightIndex)].GetXaxis().FindBin(int(sampleName.split("_")[2])),scaleNormalizationHistos["ScaleWeight"+str(weightIndex)].GetYaxis().FindBin(int(sampleName.split("_")[4]))))
			
			xsection = getattr(sbottom_masses, M_SBOTTOM).cross_section13TeV
			scalingLumi = lumi*xsection/(denominator*ISRNormalization)
			weight = denominator*ISRNormalization
			
			for signalBin in signalBins:
				
				#~ cuts = "weight*bTagWeight*ISRCorrection*leptonFastSimScaleFactor1*leptonFastSimScaleFactor2*leptonFullSimScaleFactor1*leptonFullSimScaleFactor2*(abs(deltaPhiJetMet1) > 0.4  && abs(deltaPhiJetMet2) > 0.4 && nJets > 1 && met > 150 && nUnmatchedJets == 0 && deltaR > 0.1 && chargeProduct < 0 && ((pt1 > 25 && pt2 > 20) || (pt1 > 20 && pt2 > 25)) && abs(eta1) < 2.4 && abs(eta2) < 2.4 && %s)"%(signalCuts[signalBin])
				#~ cuts = "bTagWeight*ISRCorrection*leptonFullSimScaleFactor1*leptonFullSimScaleFactor2*(abs(deltaPhiJetMet1) > 0.4  && abs(deltaPhiJetMet2) > 0.4 && nJets > 1 && met > 150 && nUnmatchedJets == 0 && deltaR > 0.1 && chargeProduct < 0 && ((pt1 > 25 && pt2 > 20) || (pt1 > 20 && pt2 > 25)) && abs(eta1) < 2.4 && abs(eta2) < 2.4 && %s)"%(signalCuts[signalBin])
				cuts = "bTagWeight*ISRCorrection*leptonFullSimScaleFactor1*leptonFullSimScaleFactor2*leptonFastSimScaleFactor1*leptonFastSimScaleFactor2*(met / caloMet < 5 && nBadMuonJets == 0 && pt > 25 && abs(deltaPhiJetMet1) > 0.4  && abs(deltaPhiJetMet2) > 0.4 && nJets > 1 && met > 150 && nUnmatchedJets == 0 && deltaR > 0.1 && chargeProduct < 0 && ((pt1 > 25 && pt2 > 20) || (pt1 > 20 && pt2 > 25)) && abs(eta1) < 2.4 && abs(eta2) < 2.4 && %s)"%(signalCuts[signalBin])
				unweightedCuts = "(met / caloMet < 5 && nBadMuonJets == 0 && pt > 25 && abs(deltaPhiJetMet1) > 0.4  && abs(deltaPhiJetMet2) > 0.4 && nJets > 1 && met > 150 && nUnmatchedJets == 0 && chargeProduct < 0 && deltaR > 0.1 && ((pt1 > 25 && pt2 > 20) || (pt1 > 20 && pt2 > 25)) && abs(eta1) < 2.4 && abs(eta2) < 2.4 && %s)"%(signalCuts[signalBin])
				
				histos = {}
				
				
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
							EEScaleUncertainty = 0
						
						
						#~ EEsystUncertainty = sqrt(EEMetUncertainty**2+EEBTagHeavyUncertainty**2+EEBTagLightUncertainty**2+EEJESUncertainty**2+EELeptonFullSimUncertainty**2+EEPileupUncertainty**2+EEISRUncertainty**2+triggerEffUncertainty**2+lumiUncertainty**2)
						EEsystUncertainty = sqrt(EEScaleUncertainty**2+EEMetUncertainty**2+EEBTagHeavyUncertainty**2+EEBTagLightUncertainty**2+EEJESUncertainty**2+EELeptonFastSimUncertainty**2+EELeptonFullSimUncertainty**2+EEPileupUncertainty**2+EEISRUncertainty**2+triggerEffUncertainty**2+lumiUncertainty**2)
						#~ EEsystUncertainty = sqrt(EEMetUncertainty**2+EEJESUncertainty**2+EELeptonFastSimUncertainty**2+EELeptonFullSimUncertainty**2+EEPileupUncertainty**2+EEISRUncertainty**2+triggerEffUncertainty**2+lumiUncertainty**2)
						
						counts["%s_EE"%signalBin] = {"Val":EEsignalYield,"MCEvents":EEMCEvents,"SignalEfficiency":EEsignalEfficiency,"StatUncertainty":EEStatUncertainty,"TotSystUncertainty":EEsystUncertainty,
						"JESUncertainty":EEJESUncertainty,"JESMean":UnscaledYield*EETriggerEff,"JESUp":EEJESUp*EETriggerEff,"JESDown":EEJESDown*EETriggerEff,
						"LeptonFastSimUncertainty":EELeptonFastSimUncertainty,"LeptonFastSimMean":UnscaledYield*EETriggerEff,"LeptonFastSimScaleFactorShifted":EELeptonFastSimScaleFactorShifted*EETriggerEff,
						"LeptonFullSimUncertainty":EELeptonFullSimUncertainty,"LeptonFullSimMean":UnscaledYield*EETriggerEff,"LeptonFullSimScaleFactorShifted":EELeptonFullSimScaleFactorShifted*EETriggerEff,
						"PileupUncertainty":EEPileupUncertainty,"PileupMean":UnscaledYield*EETriggerEff,"PileupHigh":EEPileupHigh*EETriggerEff,"PileupLow":EEPileupLow*EETriggerEff,
						"ISRUncertainty":EEISRUncertainty,"ISRMean":UnscaledYield*EETriggerEff,"ISRUp":EEISRUp*EETriggerEff,"ISRDown":EEISRDown*EETriggerEff,
						"MetUncertainty":EEMetUncertainty,"Met":EEsignalYieldMet,"GenMet":EEsignalYieldGenMet,
						"ScaleUncertainty":EEScaleUncertainty,"ScaleShifted1":EEScaleShifted[0]*EETriggerEff,"ScaleShifted2":EEScaleShifted[1]*EETriggerEff,"ScaleShifted3":EEScaleShifted[2]*EETriggerEff,"ScaleShifted4":EEScaleShifted[3]*EETriggerEff,"ScaleShifted5":EEScaleShifted[4]*EETriggerEff,"ScaleShifted6":EEScaleShifted[5]*EETriggerEff,"ScaleShifted7":EEScaleShifted[6]*EETriggerEff,"ScaleShifted8":EEScaleShifted[7]*EETriggerEff,
						"BTagHeavyUncertainty":EEBTagHeavyUncertainty,"BTagMean":UnscaledYield*EETriggerEff,"BTagHeavy":EEBTagHeavy*EETriggerEff,"BTagLightUncertainty":EEBTagLightUncertainty,"BTagLight":EEBTagLight*EETriggerEff}

				
				for sample, tree in EMTrees.iteritems():
					if sample == sampleName:
						EMusignalYieldMet = signalYields(tree,cuts, EMuTriggerEff * scalingLumi)
						EMusignalEfficiencyMet = produceSignalEfficiency(tree, cuts, weight)
						EMuMCEventsMet = MCEventYield(tree,unweightedCuts)
						EMuISRUpMet,EMuISRDownMet = produceISRUncertainty(tree,cuts)
						EMuJESUpMet,EMuJESDownMet= produceJESUncertainty(tree,cuts)
						EMuLeptonFastSimScaleFactorShiftedMet = produceLeptonFastSimUncertainty(tree,cuts)
						EMuLeptonFullSimScaleFactorShiftedMet = produceLeptonFullSimUncertainty(tree,cuts)
						EMuPileupHighMet,EMuPileupLowMet = producePileupUncertainty(tree,cuts)
						EMuBTagHeavyMet,EMuBTagLightMet = produceBTagUncertainty(tree,cuts)
						
						EMuScaleShiftedMet = []
						for weightIndex in range(0,8):
							EMuScaleShiftedMet.append(produceScaleShiftUncertainty(tree,cuts,weightIndex)/scaleNormalizations[weightIndex])
						
						if ISRNormalization > 0:
							EMuISRUpMet = EMuISRUpMet*ISRNormalization/ISRNormalizationUp
							EMuISRDownMet = EMuISRDownMet*ISRNormalization/ISRNormalizationDown
							
							EMuPileupHighMet = EMuPileupHighMet * denominator/denominatorHighPU
							EMuPileupLowMet = EMuPileupLowMet * denominator/denominatorLowPU
						
						
						cuts = cuts.replace("met","genMet")
						
						EMusignalYieldGenMet = signalYields(tree,cuts, EMuTriggerEff * scalingLumi)
						EMusignalEfficiencyGenMet = produceSignalEfficiency(tree, cuts, weight)
						EMuMCEventsGenMet = MCEventYield(tree,unweightedCuts)
						EMuISRUpGenMet,EMuISRDownGenMet = produceISRUncertainty(tree,cuts)
						EMuJESUpGenMet,EMuJESDownGenMet= produceJESUncertainty(tree,cuts)
						EMuLeptonFastSimScaleFactorShiftedGenMet = produceLeptonFastSimUncertainty(tree,cuts)
						EMuLeptonFullSimScaleFactorShiftedGenMet = produceLeptonFullSimUncertainty(tree,cuts)
						EMuPileupHighGenMet,EMuPileupLowGenMet = producePileupUncertainty(tree,cuts)
						EMuBTagHeavyGenMet,EMuBTagLightGenMet = produceBTagUncertainty(tree,cuts)
						
						EMuScaleShiftedGenMet = []
						for weightIndex in range(0,8):
							EMuScaleShiftedGenMet.append(produceScaleShiftUncertainty(tree,cuts,weightIndex)/scaleNormalizations[weightIndex])
						
						if ISRNormalization > 0:
							EMuISRUpGenMet = EMuISRUpGenMet*ISRNormalization/ISRNormalizationUp
							EMuISRDownGenMet = EMuISRDownGenMet*ISRNormalization/ISRNormalizationDown
							
							EMuPileupHighGenMet = EMuPileupHighGenMet * denominator/denominatorHighPU
							EMuPileupLowGenMet = EMuPileupLowGenMet * denominator/denominatorLowPU
						
						cuts = cuts.replace("genMet","met")
						
						EMusignalYield = 0.5* (EMusignalYieldMet + EMusignalYieldGenMet)
						EMusignalEfficiency = 0.5* (EMusignalEfficiencyMet + EMusignalEfficiencyGenMet)
						EMuMCEvents = 0.5* (EMuMCEventsMet + EMuMCEventsGenMet)
						EMuISRUp = 0.5* (EMuISRUpMet + EMuISRUpGenMet) 
						EMuISRDown = 0.5* (EMuISRDownMet + EMuISRDownGenMet) 
						EMuJESUp = 0.5* (EMuJESUpMet + EMuJESUpGenMet) 
						EMuJESDown = 0.5* (EMuJESDownMet + EMuJESDownGenMet) 
						EMuLeptonFastSimScaleFactorShifted = 0.5* (EMuLeptonFastSimScaleFactorShiftedMet + EMuLeptonFastSimScaleFactorShiftedGenMet) 
						EMuLeptonFullSimScaleFactorShifted = 0.5* (EMuLeptonFullSimScaleFactorShiftedMet + EMuLeptonFullSimScaleFactorShiftedGenMet) 
						EMuPileupHigh = 0.5* (EMuPileupHighMet + EMuPileupHighGenMet) 
						EMuPileupLow = 0.5* (EMuPileupLowMet + EMuPileupLowGenMet) 
						EMuBTagHeavy = 0.5* (EMuBTagHeavyMet + EMuBTagHeavyGenMet) 
						EMuBTagLight = 0.5* (EMuBTagLightMet + EMuBTagLightGenMet)
						
						EMuScaleShifted = []
						for weightIndex in range(0,8):
							EMuScaleShifted.append(0.5*(EMuScaleShiftedMet[weightIndex]+EMuScaleShiftedGenMet[weightIndex]))
						
						UnscaledYield = EMusignalYield / (EMuTriggerEff * scalingLumi)
						
						if UnscaledYield > 0:
							EMuISRUncertainty = max(abs(EMuISRUp-UnscaledYield)/UnscaledYield,abs(EMuISRDown-UnscaledYield)/UnscaledYield) 
							EMuJESUncertainty = max(abs(EMuJESUp-UnscaledYield)/UnscaledYield,abs(EMuJESDown-UnscaledYield)/UnscaledYield) 
							EMuPileupUncertainty = max(abs(EMuPileupHigh-UnscaledYield)/UnscaledYield,abs(EMuPileupLow-UnscaledYield)/UnscaledYield) 
							EMuLeptonFullSimUncertainty = abs(EMuLeptonFullSimScaleFactorShifted-UnscaledYield)/UnscaledYield
							EMuLeptonFastSimUncertainty = abs(EMuLeptonFastSimScaleFactorShifted-UnscaledYield)/UnscaledYield
							EMuBTagHeavyUncertainty = abs(EMuBTagHeavy-UnscaledYield)/UnscaledYield
							EMuBTagLightUncertainty = abs(EMuBTagLight-UnscaledYield)/UnscaledYield
							EMuMetUncertainty = abs(EMusignalYieldMet-EMusignalYield)/EMusignalYield
							EMuStatUncertainty = sqrt(EMuMCEvents)
							
							EMuScaleUncertainty = 0
							for weightIndex in range(0,8):
								EMuScaleUncertainty = max (EMuScaleUncertainty, abs(EMuScaleShifted[weightIndex]-UnscaledYield)/UnscaledYield)
								
						else:
							EMuISRUncertainty = 0.
							EMuJESUncertainty = 0. 
							EMuPileupUncertainty = 0.
							EMuLeptonFullSimUncertainty = 0.
							EMuLeptonFastSimUncertainty = 0.
							EMuBTagHeavyUncertainty = 0.
							EMuBTagLightUncertainty = 0.
							EMuMetUncertainty = 0.
							EMuStatUncertainty = 0.
							EMuScaleUncertainty = 0
						
						
						#~ EMusystUncertainty = sqrt(EMuMetUncertainty**2+EMuBTagHeavyUncertainty**2+EMuBTagLightUncertainty**2+EMuJESUncertainty**2+EMuLeptonFullSimUncertainty**2+EMuPileupUncertainty**2+EMuISRUncertainty**2+triggerEffUncertainty**2+lumiUncertainty**2)
						EMusystUncertainty = sqrt(EMuScaleUncertainty**2+EMuMetUncertainty**2+EMuBTagHeavyUncertainty**2+EMuBTagLightUncertainty**2+EMuJESUncertainty**2+EMuLeptonFastSimUncertainty**2+EMuLeptonFullSimUncertainty**2+EMuPileupUncertainty**2+EMuISRUncertainty**2+triggerEffUncertainty**2+lumiUncertainty**2)
						#~ EMusystUncertainty = sqrt(EMuMetUncertainty**2+EMuJESUncertainty**2+EMuLeptonFastSimUncertainty**2+EMuLeptonFullSimUncertainty**2+EMuPileupUncertainty**2+EMuISRUncertainty**2+triggerEffUncertainty**2+lumiUncertainty**2)
						
						counts["%s_EMu"%signalBin] = {"Val":EMusignalYield,"MCEvents":EMuMCEvents,"SignalEfficiency":EMusignalEfficiency,"StatUncertainty":EMuStatUncertainty,"TotSystUncertainty":EMusystUncertainty,
						"JESUncertainty":EMuJESUncertainty,"JESMean":UnscaledYield*EMuTriggerEff,"JESUp":EMuJESUp*EMuTriggerEff,"JESDown":EMuJESDown*EMuTriggerEff,
						"LeptonFastSimUncertainty":EMuLeptonFastSimUncertainty,"LeptonFastSimMean":UnscaledYield*EMuTriggerEff,"LeptonFastSimScaleFactorShifted":EMuLeptonFastSimScaleFactorShifted*EMuTriggerEff,
						"LeptonFullSimUncertainty":EMuLeptonFullSimUncertainty,"LeptonFullSimMean":UnscaledYield*EMuTriggerEff,"LeptonFullSimScaleFactorShifted":EMuLeptonFullSimScaleFactorShifted*EMuTriggerEff,
						"PileupUncertainty":EMuPileupUncertainty,"PileupMean":UnscaledYield*EMuTriggerEff,"PileupHigh":EMuPileupHigh*EMuTriggerEff,"PileupLow":EMuPileupLow*EMuTriggerEff,
						"ISRUncertainty":EMuISRUncertainty,"ISRMean":UnscaledYield*EMuTriggerEff,"ISRUp":EMuISRUp*EMuTriggerEff,"ISRDown":EMuISRDown*EMuTriggerEff,
						"MetUncertainty":EMuMetUncertainty,"Met":EMusignalYieldMet,"GenMet":EMusignalYieldGenMet,
						"ScaleUncertainty":EMuScaleUncertainty,"ScaleShifted1":EMuScaleShifted[0]*EMuTriggerEff,"ScaleShifted2":EMuScaleShifted[1]*EMuTriggerEff,"ScaleShifted3":EMuScaleShifted[2]*EMuTriggerEff,"ScaleShifted4":EMuScaleShifted[3]*EMuTriggerEff,"ScaleShifted5":EMuScaleShifted[4]*EMuTriggerEff,"ScaleShifted6":EMuScaleShifted[5]*EMuTriggerEff,"ScaleShifted7":EMuScaleShifted[6]*EMuTriggerEff,"ScaleShifted8":EMuScaleShifted[7]*EMuTriggerEff,
						"BTagHeavyUncertainty":EMuBTagHeavyUncertainty,"BTagMean":UnscaledYield*EMuTriggerEff,"BTagHeavy":EMuBTagHeavy*EMuTriggerEff,"BTagLightUncertainty":EMuBTagLightUncertainty,"BTagLight":EMuBTagLight*EMuTriggerEff}
						
				
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
						
						
						#~ MuMusystUncertainty = sqrt(MuMuMetUncertainty**2+MuMuBTagHeavyUncertainty**2+MuMuBTagLightUncertainty**2+MuMuJESUncertainty**2+MuMuLeptonFullSimUncertainty**2+MuMuPileupUncertainty**2+MuMuISRUncertainty**2+triggerEffUncertainty**2+lumiUncertainty**2)
						MuMusystUncertainty = sqrt(MuMuScaleUncertainty**2+MuMuMetUncertainty**2+MuMuBTagHeavyUncertainty**2+MuMuBTagLightUncertainty**2+MuMuJESUncertainty**2+MuMuLeptonFastSimUncertainty**2+MuMuLeptonFullSimUncertainty**2+MuMuPileupUncertainty**2+MuMuISRUncertainty**2+triggerEffUncertainty**2+lumiUncertainty**2)
						#~ MuMusystUncertainty = sqrt(MuMuMetUncertainty**2+MuMuJESUncertainty**2+MuMuLeptonFastSimUncertainty**2+MuMuLeptonFullSimUncertainty**2+MuMuPileupUncertainty**2+MuMuISRUncertainty**2+triggerEffUncertainty**2+lumiUncertainty**2)
						
						counts["%s_MuMu"%signalBin] = {"Val":MuMusignalYield,"MCEvents":MuMuMCEvents,"SignalEfficiency":MuMusignalEfficiency,"StatUncertainty":MuMuStatUncertainty,"TotSystUncertainty":MuMusystUncertainty,
						"JESUncertainty":MuMuJESUncertainty,"JESMean":UnscaledYield*MuMuTriggerEff,"JESUp":MuMuJESUp*MuMuTriggerEff,"JESDown":MuMuJESDown*MuMuTriggerEff,
						"LeptonFastSimUncertainty":MuMuLeptonFastSimUncertainty,"LeptonFastSimMean":UnscaledYield*MuMuTriggerEff,"LeptonFastSimScaleFactorShifted":MuMuLeptonFastSimScaleFactorShifted*MuMuTriggerEff,
						"LeptonFullSimUncertainty":MuMuLeptonFullSimUncertainty,"LeptonFullSimMean":UnscaledYield*MuMuTriggerEff,"LeptonFullSimScaleFactorShifted":MuMuLeptonFullSimScaleFactorShifted*MuMuTriggerEff,
						"PileupUncertainty":MuMuPileupUncertainty,"PileupMean":UnscaledYield*MuMuTriggerEff,"PileupHigh":MuMuPileupHigh*MuMuTriggerEff,"PileupLow":MuMuPileupLow*MuMuTriggerEff,
						"ISRUncertainty":MuMuISRUncertainty,"ISRMean":UnscaledYield*MuMuTriggerEff,"ISRUp":MuMuISRUp*MuMuTriggerEff,"ISRDown":MuMuISRDown*MuMuTriggerEff,
						"ISRUncertainty":MuMuISRUncertainty,"ISRMean":UnscaledYield*MuMuTriggerEff,"ISRUp":MuMuISRUp*MuMuTriggerEff,"ISRDown":MuMuISRDown*MuMuTriggerEff,
						"MetUncertainty":MuMuMetUncertainty,"Met":MuMusignalYieldMet,"GenMet":MuMusignalYieldGenMet,
						"ScaleUncertainty":MuMuScaleUncertainty,"ScaleShifted1":MuMuScaleShifted[0]*MuMuTriggerEff,"ScaleShifted2":MuMuScaleShifted[1]*MuMuTriggerEff,"ScaleShifted3":MuMuScaleShifted[2]*MuMuTriggerEff,"ScaleShifted4":MuMuScaleShifted[3]*MuMuTriggerEff,"ScaleShifted5":MuMuScaleShifted[4]*MuMuTriggerEff,"ScaleShifted6":MuMuScaleShifted[5]*MuMuTriggerEff,"ScaleShifted7":MuMuScaleShifted[6]*MuMuTriggerEff,"ScaleShifted8":MuMuScaleShifted[7]*MuMuTriggerEff,
						"BTagHeavyUncertainty":MuMuBTagHeavyUncertainty,"BTagMean":UnscaledYield*MuMuTriggerEff,"BTagHeavy":MuMuBTagHeavy*MuMuTriggerEff,"BTagLightUncertainty":MuMuBTagLightUncertainty,"BTagLight":MuMuBTagLight*MuMuTriggerEff}
						
			
			outFilePkl = open("shelvesSystematics/%s.pkl"%(sampleName),"w")
			#~ outFilePkl = open("shelvesSystematics17fb/%s.pkl"%(sampleName),"w")
			pickle.dump(counts, outFilePkl)
			outFilePkl.close()
			
			m_n += stepsize
		m_b += stepsize
			
			

				
	

# This method just waits for a button to be pressed
def waitForInput():
    raw_input("Press any key to continue!")
    return

# entry point
#-------------
if (__name__ == "__main__"):
    # use option parser to allow verbose mode
    
    # start
    plot()
 
