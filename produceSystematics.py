####################################################################################################
# FastSim lepton scaling factors were not in included in the trees                                 #
# -> One has to run over each event to scale it properly                                           #
# For the slepton model (fixed neutralino 1): Events have to be scaled to Z/slepton BR either via  #
# the number of generated sleptons (new) or the number of generated leptons 					   #
# (with a neutralino 2 as mother) and the mother PDG ID (old)                                      #
####################################################################################################

import sys
sys.path.append('cfg/')
from frameworkStructure import pathes
sys.path.append(pathes.basePath)

from corrections import triggerEffs, rSFOF


import ROOT
from ROOT import TH1F
from math import sqrt
attic = []


ROOT.gStyle.SetOptStat(0)

signalRegionCuts = {
			"EdgeLegacy":"((nJets >= 2 && met > 150) || (nJets >= 3 && met > 100)) && mll < 70 && mll > 20 && abs(eta1) < 1.4 && abs(eta2) < 1.4 && deltaR > 0.3",
			"LowMassLowNll":"nJets >= 2 && met > 150 && mll < 81 && mll > 20 && nLL < 21 && deltaR > 0.1",
			"LowMassHighNll":"nJets >= 2 && met > 150 && mll < 81 && mll > 20 && nLL >= 21 && deltaR > 0.1",
			"HighMassLowNll":"nJets >= 2 && met > 150 && mll > 101 && nLL < 21 && deltaR > 0.1",
			"HighMassHighNll":"nJets >= 2 && met > 150 && mll > 101 && nLL >= 21 && deltaR > 0.1",
			}
			

def readTreeFromFile(path, dileptonCombination):
	"""
	helper functionfrom argparse import ArgumentParser
	path: path to .root file containing simulated events
	dileptonCombination: EMu, EMu, or EMu for electron-electron, electron-muon, or muon-muon events

	returns: tree containing events for on sample and dileptonCombination
	"""
	from ROOT import TChain
	result = TChain()
	result.Add("%s/cutsV33DileptonFinalTrees/%sDileptonTree"%(path, dileptonCombination))
	return result
	
def getFilePathsAndSampleNames(path):
	"""
	helper function
	path: path to directory containing all sample files

	returns: dict of smaple names -> path of .root file (for all samples in path)
	"""
	result = []
	from glob import glob
	from re import match
	result = {}
	for filePath in glob("%s/sw8014*.root"%path):
		sampleName = match(".*sw8014v1009.processed\.(.*).root", filePath).groups()[0]
		#for the python enthusiats: yield sampleName, filePath is more efficient here :)
		result[sampleName] = filePath
	return result
	
def totalNumberOfGeneratedEvents(path):
	"""
	path: path to directory containing all sample files

	returns dict samples names -> number of simulated events in source sample
			(note these include events without EMu EMu EMu signature, too )
	"""
	from ROOT import TFile
	result = {}
	for sampleName, filePath in getFilePathsAndSampleNames(path).iteritems():
		rootFile = TFile(filePath, "read")
		result[sampleName] = rootFile.FindObjectAny("analysis paths").GetBinContent(1)
	return result
	
def readTrees(path, dileptonCombination):
	"""
	path: path to directory containing all sample files
	dileptonCombination: "EMu", "EMu", or pyroot"EMu" for electron-electron, electron-muon, or muon-muon events

	returns: dict of sample names ->  trees containing events (for all samples for one dileptonCombination)
	"""
	result = {}
	for sampleName, filePath in getFilePathsAndSampleNames(path).iteritems():		
		result[sampleName] = readTreeFromFile(filePath, dileptonCombination)
		
	return result
	
	
def createHistoFromTree(tree, variable, weight, nBins, firstBin, lastBin, nEvents = -1):
	"""
	tree: tree to create histo from)
	variable: variable to plot (must be a branch of the tree)
	weight: weights to apply (e.g. "var1*(var2 > 15)" will use weights from var1 and cut on var2 > 15
	nBins, firstBin, lastBin: number of bins, first bin and last bin (same as in TH1F constructor)
	nEvents: number of events to process (-1 = all)
	"""
	from ROOT import TH1F
	from random import randint
	from sys import maxint
	if nEvents < 0:
		nEvents = maxint
	#make a random name you could give something meaningfull here,
	#but that would make this less readable
	name = "%x"%(randint(0, maxint))
	result = TH1F(name, "", nBins, firstBin, lastBin)
	result.Sumw2()
	tree.Draw("%s>>%s"%(variable, name), weight, "goff", nEvents)
	for i in range(0,nBins+1):
		if TMath.IsNaN(result.GetBinContent(i)) or result.GetBinContent(i)<0:
			result.SetBinContent(i,0)
	return result
	
def signalYields(tree,cuts,scalingLumi):
	histo = createHistoFromTree(tree, "bTagWeight", cuts, 400, 0, 40)
	histo.Scale(scalingLumi)
	yields = float(histo.Integral())
	return yields
			
def produceSignalEfficiency(tree,cuts,weight):
	yields = [float(createHistoFromTree(tree, "bTagWeight", cuts, 400, 0, 4).Integral())]
	result = yields[0]/weight
	return result
			
def statisticalUncertainty(tree,cuts):
	histo = createHistoFromTree(tree, "bTagWeight", cuts, 400, 0, 4)
	Int = histo.Integral(0,histo.GetNbinsX()+1)
	return Int	

def producePileupUncertainty(tree,cuts):
	
	uncertaintySources = ["weightUp*","weightDown*"]
	
	yields = []
	for index, source in enumerate(uncertaintySources):
		if index == 0:
			cuts = cuts.replace("weight*",source)
		else:
			cuts = cuts.replace(uncertaintySources[index-1],source)
		yields.append(float(createHistoFromTree(tree, "bTagWeight", cuts, 400, 0, 4).Integral()))

	return yields[0],yields[1]	
	
def produceJESUncertainty(tree,cuts):

	cuts = cuts.replace("nJets","nShiftedJetsJESUp")	
	cuts = cuts.replace("met","metJESUp")	
	cuts = cuts.replace("nLL","nLLJESUp")	
	yields = [float(createHistoFromTree(tree, "bTagWeight", cuts, 400, 0, 4).Integral())]
	
	cuts = cuts.replace("nShiftedJetsJESUp","nShiftedJetsJESDown")	
	cuts = cuts.replace("metJESUp","metJESDown")	
	cuts = cuts.replace("nLLJESUp","nLLJESDown")	
	yields.append(float(createHistoFromTree(tree, "bTagWeight", cuts, 400, 0, 4).Integral()))


	return yields[0],yields[1]

			
def produceISRUncertainty(tree,cuts):
	
	cuts = cuts.replace("ISRCorrection","(ISRCorrection+ISRUncertainty)")
	yields = [float(createHistoFromTree(tree, "bTagWeight", cuts, 400, 0, 4).Integral())]

	cuts = cuts.replace("(ISRCorrection+ISRUncertainty)","(ISRCorrection-ISRUncertainty)")
	yields.append(float(createHistoFromTree(tree, "bTagWeight", cuts, 400, 0, 4).Integral()))
	
	return yields[0],yields[1]
			
def produceLeptonFastSimUncertainty(tree,cuts):
	
	result = 0.
	cuts = cuts.replace("leptonFastSimScaleFactor1*leptonFastSimScaleFactor2*","")
	result = float(createHistoFromTree(tree, "bTagWeight", cuts, 400, 0, 4).Integral())
	
	return result
	
def produceLeptonFullSimUncertainty(tree,cuts):
	
	result = 0.
	cuts = cuts.replace("leptonFullSimScaleFactor1*leptonFullSimScaleFactor2*","")
	result = float(createHistoFromTree(tree, "bTagWeight", cuts, 400, 0, 4).Integral())
	
	return result
	
def produceBTagUncertainty(tree,cuts):
	
	cuts = cuts.replace("bTagWeight","bTagWeight * (1 + bTagWeightErrHeavy)")
	yields = [float(createHistoFromTree(tree, "bTagWeight", cuts, 400, 0, 4).Integral())]
	
	cuts = cuts.replace("bTagWeight * (1 + bTagWeightErrHeavy)","bTagWeight * (1 + bTagWeightErrLight)")
	yields.append(float(createHistoFromTree(tree, "bTagWeight", cuts, 400, 0, 4).Integral()))
	
	
	return yields[0],yields[1]
	

		
if (__name__ == "__main__"):
	import sys
	sys.path.append('cfg/')
	from frameworkStructure import pathes
	sys.path.append(pathes.basePath)
	from sys import argv
	import pickle	
	from ROOT import TCanvas, TPad, TH1F, TH1I, THStack, TLegend, TF1, TH2F, TH2D, TFile, TMath
	import ratios
	from defs import sbottom_masses
	from math import sqrt
	
	denominatorFile = TFile("T6bbllsleptonDenominatorHisto.root")
	denominatorHisto = TH2F(denominatorFile.Get("massScan"))
	ISRNormalizationHisto = TH2F(denominatorFile.Get("ISRNormalization"))
	ISRNormalizationHistoUp = TH2F(denominatorFile.Get("ISRNormalizationUp"))
	ISRNormalizationHistoDown = TH2F(denominatorFile.Get("ISRNormalizationDown"))

	lumi = 12900
		
	m_b = int(argv[1])
	m_n_max = 900
	
	if m_b < 800:
		step_size = 25
	else:
		step_size = 50
		
	path = "/disk1/user/schomakers/trees/sw8014v1009_NLL"
	m_n_min = 150
	
	signalRegions = ["EdgeLegacy","LowMassLowNll","LowMassHighNll","HighMassLowNll","HighMassHighNll"]

	for region in signalRegions:
		signalRegionCut = signalRegionCuts[region]
		
		cuts = "weight*ISRCorrection*leptonFastSimScaleFactor1*leptonFastSimScaleFactor2*leptonFullSimScaleFactor1*leptonFullSimScaleFactor2*bTagWeight*( nUnmatchedJets == 0 && chargeProduct < 0 && ((pt1 > 25 && pt2 > 20) || (pt1 > 20 && pt2 > 25)) && abs(eta1) < 2.4 && abs(eta2) < 2.4 && %s)"%(signalRegionCut)
		unweightedCuts = "( nUnmatchedJets == 0 && chargeProduct < 0 && ((pt1 > 25 && pt2 > 20) || (pt1 > 20 && pt2 > 25)) && abs(eta1) < 2.4 && abs(eta2) < 2.4 && %s)"%(signalRegionCut)
		
		if region == "edgeLegacy":
			EETriggerEff = triggerEffs.central.effEE.val
			EMuTriggerEff = triggerEffs.central.effEM.val
			MuMuTriggerEff = triggerEffs.central.effMM.val
			RSFOF = rSFOF.central.val
		else:
			EETriggerEff = triggerEffs.inclusive.effEE.val
			EMuTriggerEff = triggerEffs.inclusive.effEM.val
			MuMuTriggerEff = triggerEffs.inclusive.effMM.val
			RSFOF = rSFOF.inclusive.val

	
		triggerEffUncertainty = 0.05
		lumiUncertainty = 0.062
		electronUncertainty = 0.05
		muonUncertainty = 0.03
					
	
	
		M_SBOTTOM = "m_b_"+str(m_b)
		m_sbottom = str(m_b)
		j = 0
		m_n_2 = m_n_min
		while m_n_min + j*step_size <= m_n_max:
			m_n_2 = m_n_min + j*step_size
			m_neutralino_2 = str(m_n_min + j*step_size)
			
			if m_b > m_n_2:
				print "m_n: "+m_neutralino_2
				print cuts
				
				if not ((m_b == 775 and m_n_2 == 750) or (m_b == 800 and m_n_2 == 150) or (m_b == 950 and m_n_2 == 900) or (m_b == 950 and m_n_2 == 850) or (m_b == 950 and m_n_2 == 550) or (m_b == 950 and m_n_2 == 500) or (m_b == 950 and m_n_2 == 300) or (m_b == 950 and m_n_2 == 250)):
			
					denominator = denominatorHisto.GetBinContent(denominatorHisto.GetXaxis().FindBin(m_b),denominatorHisto.GetYaxis().FindBin(m_n_2))
					ISRNormalization = ISRNormalizationHisto.GetBinContent(ISRNormalizationHisto.GetXaxis().FindBin(m_b),ISRNormalizationHisto.GetYaxis().FindBin(m_n_2))
					ISRNormalizationUp = ISRNormalizationHistoUp.GetBinContent(ISRNormalizationHisto.GetXaxis().FindBin(m_b),ISRNormalizationHisto.GetYaxis().FindBin(m_n_2))
					ISRNormalizationDown = ISRNormalizationHistoDown.GetBinContent(ISRNormalizationHisto.GetXaxis().FindBin(m_b),ISRNormalizationHisto.GetYaxis().FindBin(m_n_2))
					
					EMutrees = readTrees(path, "EMu")	
					EEtrees = readTrees(path, "EE")	
					MuMutrees = readTrees(path, "MuMu")	
				
					sampleName = "T6bbllslepton_msbottom_%s_mneutralino_%s"%(m_sbottom,m_neutralino_2)
					fileName = "T6bbllslepton_msbottom_%s_mneutralino_%s"%(m_sbottom,m_neutralino_2)
		
					xsection = getattr(sbottom_masses, M_SBOTTOM).cross_section13TeV
					
					scalingLumi = lumi*xsection*ISRNormalization/denominator
					weight = denominator*ISRNormalization
					
					for sample, tree in EEtrees.iteritems():
						if sample == sampleName:
							EEsignalYieldMet = signalYields(tree,cuts, EETriggerEff * scalingLumi)
							EEsignalEfficiencyMet = produceSignalEfficiency(tree, cuts, weight)
							EEMCEventsMet = statisticalUncertainty(tree,unweightedCuts)
							EEISRUpMet,EEISRDownMet = produceISRUncertainty(tree,cuts)
							EEJESUpMet,EEJESDownMet= produceJESUncertainty(tree,cuts)
							EELeptonNoFastSimScaleFactorMet = produceLeptonFastSimUncertainty(tree,cuts)
							EELeptonNoFullSimScaleFactorMet = produceLeptonFullSimUncertainty(tree,cuts)
							EEPileupUpMet,EEPileupDownMet = producePileupUncertainty(tree,cuts)
							EEBTagHeavyMet,EEBTagLightMet = produceBTagUncertainty(tree,cuts)
							
							if ISRNormalization > 0:
								EEISRUpMet = EEISRUpMet*ISRNormalizationUp/ISRNormalization
								EEISRDownMet = EEISRDownMet*ISRNormalizationDown/ISRNormalization
							
							
							cuts = cuts.replace("met","genMet")
							
							EEsignalYieldGenMet = signalYields(tree,cuts, EETriggerEff * scalingLumi)
							EEsignalEfficiencyGenMet = produceSignalEfficiency(tree, cuts, weight)
							EEMCEventsGenMet = statisticalUncertainty(tree,unweightedCuts)
							EEISRUpGenMet,EEISRDownGenMet = produceISRUncertainty(tree,cuts)
							EEJESUpGenMet,EEJESDownGenMet= produceJESUncertainty(tree,cuts)
							EELeptonNoFastSimScaleFactorGenMet = produceLeptonFastSimUncertainty(tree,cuts)
							EELeptonNoFullSimScaleFactorGenMet = produceLeptonFullSimUncertainty(tree,cuts)
							EEPileupUpGenMet,EEPileupDownGenMet = producePileupUncertainty(tree,cuts)
							EEBTagHeavyGenMet,EEBTagLightGenMet = produceBTagUncertainty(tree,cuts)
							
							if ISRNormalization > 0:
								EEISRUpGenMet = EEISRUpGenMet*ISRNormalizationUp/ISRNormalization
								EEISRDownGenMet = EEISRDownGenMet*ISRNormalizationDown/ISRNormalization
							
							cuts = cuts.replace("genMet","met")
							
							EEsignalYield = 0.5* (EEsignalYieldMet + EEsignalYieldGenMet)
							EEsignalEfficiency = 0.5* (EEsignalEfficiencyMet + EEsignalEfficiencyGenMet)
							EEMCEvents = 0.5* (EEMCEventsMet + EEMCEventsGenMet)
							EEISRUp = 0.5* (EEISRUpMet + EEISRUpGenMet) 
							EEISRDown = 0.5* (EEISRDownMet + EEISRDownGenMet) 
							EEJESUp = 0.5* (EEJESUpMet + EEJESUpGenMet) 
							EEJESDown = 0.5* (EEJESDownMet + EEJESDownGenMet) 
							EELeptonNoFastSimScaleFactor = 0.5* (EELeptonNoFastSimScaleFactorMet + EELeptonNoFastSimScaleFactorGenMet) 
							EELeptonNoFullSimScaleFactor = 0.5* (EELeptonNoFullSimScaleFactorMet + EELeptonNoFullSimScaleFactorGenMet) 
							EEPileupUp = 0.5* (EEPileupUpMet + EEPileupUpGenMet) 
							EEPileupDown = 0.5* (EEPileupDownMet + EEPileupUpGenMet) 
							EEBTagHeavy = 0.5* (EEBTagHeavyMet + EEBTagHeavyGenMet) 
							EEBTagLight = 0.5* (EEBTagLightMet + EEBTagLightGenMet)
							
							UnscaledYield = EEsignalYield / (EETriggerEff * scalingLumi)
							
							if UnscaledYield > 0:
								EEISRUncertainty = max(abs(EEISRUp-UnscaledYield)/UnscaledYield,abs(EEISRDown-UnscaledYield)/UnscaledYield) 
								EEJESUncertainty = max(abs(EEJESUp-UnscaledYield)/UnscaledYield,abs(EEJESDown-UnscaledYield)/UnscaledYield) 
								EEPileupUncertainty = max(abs(EEPileupUp-UnscaledYield)/UnscaledYield,abs(EEPileupDown-UnscaledYield)/UnscaledYield) 
								EELeptonFullSimUncertainty = abs(EELeptonNoFullSimScaleFactor-UnscaledYield)/UnscaledYield
								EELeptonFastSimUncertainty = abs(EELeptonNoFastSimScaleFactor-UnscaledYield)/UnscaledYield
								EEBTagHeavyUncertainty = abs(EEBTagHeavy-UnscaledYield)/UnscaledYield
								EEBTagLightUncertainty = abs(EEBTagLight-UnscaledYield)/UnscaledYield
								EEMetUncertainty = abs(EEsignalYieldMet-EEsignalYield)/EEsignalYield
								EEStatUncertainty = sqrt(EEMCEvents)
							else:
								EEISRUncertainty = 0.
								EEJESUncertainty = 0. 
								EEPileupUncertainty = 0.
								EELeptonNoFullSimUncertainty = 0.
								EELeptonNoFastSimUncertainty = 0.
								EEBTagHeavyUncertainty = 0.
								EEBTagLightUncertainty = 0.
								EEMetUncertainty = 0.
								EEStatUncertainty = 0.
							
							
							EEsystUncertainty = sqrt(EEMetUncertainty**2+EEBTagHeavyUncertainty**2+EEBTagLightUncertainty**2+EEJESUncertainty**2+EELeptonFastSimUncertainty**2+EELeptonFullSimUncertainty**2+EEPileupUncertainty**2+EEISRUncertainty**2+triggerEffUncertainty**2+lumiUncertainty**2)
							
							counts = {}
							counts["EE"] = {"EEval":EEsignalYield,"EEMCEvents":EEMCEvents,"EEsignalEfficiency":EEsignalEfficiency,"EEstatUncertainty":EEStatUncertainty,"EETotSystUncertainty":EEsystUncertainty,
							"EEJESUncertainty":EEJESUncertainty,"EEJESMean":UnscaledYield*EETriggerEff,"EEJESUp":EEJESUp*EETriggerEff,"EEJESDown":EEJESDown*EETriggerEff,
							"EELeptonFastSimUncertainty":EELeptonFastSimUncertainty,"EELeptonFastSimMean":UnscaledYield*EETriggerEff,"EELeptonNoFastSimScaleFactor":EELeptonNoFastSimScaleFactor*EETriggerEff,
							"EELeptonFullSimUncertainty":EELeptonFullSimUncertainty,"EELeptonFullSimMean":UnscaledYield*EETriggerEff,"EELeptonNoFullSimScaleFactor":EELeptonNoFullSimScaleFactor*EETriggerEff,
							"EEpileupUncertainty":EEPileupUncertainty,"EEPileupMean":UnscaledYield*EETriggerEff,"EEPileupUp":EEPileupUp*EETriggerEff,"EEPileupDown":EEPileupDown*EETriggerEff,
							"EEISRUncertainty":EEISRUncertainty,"EEISRMean":UnscaledYield*EETriggerEff,"EEISRUp":EEISRUp*EETriggerEff,"EEISRDown":EEISRDown*EETriggerEff,
							"EEMetUncertainty":EEMetUncertainty,"EEMet":EEsignalYieldMet,"EEGenMet":EEsignalYieldGenMet,
							"EEBTagHeavyUncertainty":EEBTagHeavyUncertainty,"EEbTagMean":UnscaledYield*EETriggerEff,"EEbTagHeavy":EEBTagHeavy*EETriggerEff,"EEbTagLightUncertainty":EEBTagLightUncertainty,"EEbTagLight":EEBTagLight*EETriggerEff}
							
							outFilePkl = open("shelves/%s_%s_EE.pkl"%(fileName,region),"w")
							pickle.dump(counts, outFilePkl)
							outFilePkl.close()
					
					for sample, tree in EMutrees.iteritems():
						if sample == sampleName:
							EMusignalYieldMet = signalYields(tree,cuts, EMuTriggerEff * scalingLumi)
							EMusignalEfficiencyMet = produceSignalEfficiency(tree, cuts, weight)
							EMuMCEventsMet = statisticalUncertainty(tree,unweightedCuts)
							EMuISRUpMet,EMuISRDownMet = produceISRUncertainty(tree,cuts)
							EMuJESUpMet,EMuJESDownMet= produceJESUncertainty(tree,cuts)
							EMuLeptonNoFastSimScaleFactorMet = produceLeptonFastSimUncertainty(tree,cuts)
							EMuLeptonNoFullSimScaleFactorMet = produceLeptonFullSimUncertainty(tree,cuts)
							EMuPileupUpMet,EMuPileupDownMet = producePileupUncertainty(tree,cuts)
							EMuBTagHeavyMet,EMuBTagLightMet = produceBTagUncertainty(tree,cuts)
							
							if ISRNormalization > 0:
								EMuISRUpMet = EMuISRUpMet*ISRNormalizationUp/ISRNormalization
								EMuISRDownMet = EMuISRDownMet*ISRNormalizationDown/ISRNormalization
							
							
							cuts = cuts.replace("met","genMet")
							
							EMusignalYieldGenMet = signalYields(tree,cuts, EMuTriggerEff * scalingLumi)
							EMusignalEfficiencyGenMet = produceSignalEfficiency(tree, cuts, weight)
							EMuMCEventsGenMet = statisticalUncertainty(tree,unweightedCuts)
							EMuISRUpGenMet,EMuISRDownGenMet = produceISRUncertainty(tree,cuts)
							EMuJESUpGenMet,EMuJESDownGenMet= produceJESUncertainty(tree,cuts)
							EMuLeptonNoFastSimScaleFactorGenMet = produceLeptonFastSimUncertainty(tree,cuts)
							EMuLeptonNoFullSimScaleFactorGenMet = produceLeptonFullSimUncertainty(tree,cuts)
							EMuPileupUpGenMet,EMuPileupDownGenMet = producePileupUncertainty(tree,cuts)
							EMuBTagHeavyGenMet,EMuBTagLightGenMet = produceBTagUncertainty(tree,cuts)
							
							if ISRNormalization > 0:
								EMuISRUpGenMet = EMuISRUpGenMet*ISRNormalizationUp/ISRNormalization
								EMuISRDownGenMet = EMuISRDownGenMet*ISRNormalizationDown/ISRNormalization
							
							cuts = cuts.replace("genMet","met")
							
							EMusignalYield = 0.5* (EMusignalYieldMet + EMusignalYieldGenMet)
							EMusignalEfficiency = 0.5* (EMusignalEfficiencyMet + EMusignalEfficiencyGenMet)
							EMuMCEvents = 0.5* (EMuMCEventsMet + EMuMCEventsGenMet)
							EMuISRUp = 0.5* (EMuISRUpMet + EMuISRUpGenMet) 
							EMuISRDown = 0.5* (EMuISRDownMet + EMuISRDownGenMet) 
							EMuJESUp = 0.5* (EMuJESUpMet + EMuJESUpGenMet) 
							EMuJESDown = 0.5* (EMuJESDownMet + EMuJESDownGenMet) 
							EMuLeptonNoFastSimScaleFactor = 0.5* (EMuLeptonNoFastSimScaleFactorMet + EMuLeptonNoFastSimScaleFactorGenMet) 
							EMuLeptonNoFullSimScaleFactor = 0.5* (EMuLeptonNoFullSimScaleFactorMet + EMuLeptonNoFullSimScaleFactorGenMet) 
							EMuPileupUp = 0.5* (EMuPileupUpMet + EMuPileupUpGenMet) 
							EMuPileupDown = 0.5* (EMuPileupDownMet + EMuPileupUpGenMet) 
							EMuBTagHeavy = 0.5* (EMuBTagHeavyMet + EMuBTagHeavyGenMet) 
							EMuBTagLight = 0.5* (EMuBTagLightMet + EMuBTagLightGenMet)
							
							UnscaledYield = EMusignalYield / (EMuTriggerEff * scalingLumi)
							
							if UnscaledYield > 0:
								EMuISRUncertainty = max(abs(EMuISRUp-UnscaledYield)/UnscaledYield,abs(EMuISRDown-UnscaledYield)/UnscaledYield) 
								EMuJESUncertainty = max(abs(EMuJESUp-UnscaledYield)/UnscaledYield,abs(EMuJESDown-UnscaledYield)/UnscaledYield) 
								EMuPileupUncertainty = max(abs(EMuPileupUp-UnscaledYield)/UnscaledYield,abs(EMuPileupDown-UnscaledYield)/UnscaledYield) 
								EMuLeptonFullSimUncertainty = abs(EMuLeptonNoFullSimScaleFactor-UnscaledYield)/UnscaledYield
								EMuLeptonFastSimUncertainty = abs(EMuLeptonNoFastSimScaleFactor-UnscaledYield)/UnscaledYield
								EMuBTagHeavyUncertainty = abs(EMuBTagHeavy-UnscaledYield)/UnscaledYield
								EMuBTagLightUncertainty = abs(EMuBTagLight-UnscaledYield)/UnscaledYield
								EMuMetUncertainty = abs(EMusignalYieldMet-EMusignalYield)/EMusignalYield
								EMuStatUncertainty = sqrt(EMuMCEvents)
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
							
							
							EMusystUncertainty = sqrt(EMuMetUncertainty**2+EMuBTagHeavyUncertainty**2+EMuBTagLightUncertainty**2+EMuJESUncertainty**2+EMuLeptonFastSimUncertainty**2+EMuLeptonFullSimUncertainty**2+EMuPileupUncertainty**2+EMuISRUncertainty**2+triggerEffUncertainty**2+lumiUncertainty**2)
							
							counts = {}
							counts["EMu"] = {"EMuval":EMusignalYield,"EMuMCEvents":EMuMCEvents,"EMusignalEfficiency":EMusignalEfficiency,"EMustatUncertainty":EMuStatUncertainty,"EMuTotSystUncertainty":EMusystUncertainty,
							"EMuJESUncertainty":EMuJESUncertainty,"EMuJESMean":UnscaledYield*EMuTriggerEff,"EMuJESUp":EMuJESUp*EMuTriggerEff,"EMuJESDown":EMuJESDown*EMuTriggerEff,
							"EMuLeptonFastSimUncertainty":EMuLeptonFastSimUncertainty,"EMuLeptonFastSimMean":UnscaledYield*EMuTriggerEff,"EMuLeptonNoFastSimScaleFactor":EMuLeptonNoFastSimScaleFactor*EMuTriggerEff,
							"EMuLeptonFullSimUncertainty":EMuLeptonFullSimUncertainty,"EMuLeptonFullSimMean":UnscaledYield*EMuTriggerEff,"EMuLeptonNoFullSimScaleFactor":EMuLeptonNoFullSimScaleFactor*EMuTriggerEff,
							"EMupileupUncertainty":EMuPileupUncertainty,"EMuPileupMean":UnscaledYield*EMuTriggerEff,"EMuPileupUp":EMuPileupUp*EMuTriggerEff,"EMuPileupDown":EMuPileupDown*EMuTriggerEff,
							"EMuISRUncertainty":EMuISRUncertainty,"EMuISRMean":UnscaledYield*EMuTriggerEff,"EMuISRUp":EMuISRUp*EMuTriggerEff,"EMuISRDown":EMuISRDown*EMuTriggerEff,
							"EMuMetUncertainty":EMuMetUncertainty,"EMuMet":EMusignalYieldMet,"EMuGenMet":EMusignalYieldGenMet,
							"EMuBTagHeavyUncertainty":EMuBTagHeavyUncertainty,"EMubTagMean":UnscaledYield*EMuTriggerEff,"EMubTagHeavy":EMuBTagHeavy*EMuTriggerEff,"EMubTagLightUncertainty":EMuBTagLightUncertainty,"EMubTagLight":EMuBTagLight*EMuTriggerEff}
							
							outFilePkl = open("shelves/%s_%s_EMu.pkl"%(fileName,region),"w")
							pickle.dump(counts, outFilePkl)
							outFilePkl.close()
					
					for sample, tree in MuMutrees.iteritems():
						if sample == sampleName:
							MuMusignalYieldMet = signalYields(tree,cuts, MuMuTriggerEff * scalingLumi)
							MuMusignalEfficiencyMet = produceSignalEfficiency(tree, cuts, weight)
							MuMuMCEventsMet = statisticalUncertainty(tree,unweightedCuts)
							MuMuISRUpMet,MuMuISRDownMet = produceISRUncertainty(tree,cuts)
							MuMuJESUpMet,MuMuJESDownMet= produceJESUncertainty(tree,cuts)
							MuMuLeptonNoFastSimScaleFactorMet = produceLeptonFastSimUncertainty(tree,cuts)
							MuMuLeptonNoFullSimScaleFactorMet = produceLeptonFullSimUncertainty(tree,cuts)
							MuMuPileupUpMet,MuMuPileupDownMet = producePileupUncertainty(tree,cuts)
							MuMuBTagHeavyMet,MuMuBTagLightMet = produceBTagUncertainty(tree,cuts)
							
							if ISRNormalization > 0:
								MuMuISRUpMet = MuMuISRUpMet*ISRNormalizationUp/ISRNormalization
								MuMuISRDownMet = MuMuISRDownMet*ISRNormalizationDown/ISRNormalization
							
							
							cuts = cuts.replace("met","genMet")
							
							MuMusignalYieldGenMet = signalYields(tree,cuts, MuMuTriggerEff * scalingLumi)
							MuMusignalEfficiencyGenMet = produceSignalEfficiency(tree, cuts, weight)
							MuMuMCEventsGenMet = statisticalUncertainty(tree,unweightedCuts)
							MuMuISRUpGenMet,MuMuISRDownGenMet = produceISRUncertainty(tree,cuts)
							MuMuJESUpGenMet,MuMuJESDownGenMet= produceJESUncertainty(tree,cuts)
							MuMuLeptonNoFastSimScaleFactorGenMet = produceLeptonFastSimUncertainty(tree,cuts)
							MuMuLeptonNoFullSimScaleFactorGenMet = produceLeptonFullSimUncertainty(tree,cuts)
							MuMuPileupUpGenMet,MuMuPileupDownGenMet = producePileupUncertainty(tree,cuts)
							MuMuBTagHeavyGenMet,MuMuBTagLightGenMet = produceBTagUncertainty(tree,cuts)
							
							if ISRNormalization > 0:
								MuMuISRUpGenMet = MuMuISRUpGenMet*ISRNormalizationUp/ISRNormalization
								MuMuISRDownGenMet = MuMuISRDownGenMet*ISRNormalizationDown/ISRNormalization
							
							cuts = cuts.replace("genMet","met")
							
							MuMusignalYield = 0.5* (MuMusignalYieldMet + MuMusignalYieldGenMet)
							MuMusignalEfficiency = 0.5* (MuMusignalEfficiencyMet + MuMusignalEfficiencyGenMet)
							MuMuMCEvents = 0.5* (MuMuMCEventsMet + MuMuMCEventsGenMet)
							MuMuISRUp = 0.5* (MuMuISRUpMet + MuMuISRUpGenMet) 
							MuMuISRDown = 0.5* (MuMuISRDownMet + MuMuISRDownGenMet) 
							MuMuJESUp = 0.5* (MuMuJESUpMet + MuMuJESUpGenMet) 
							MuMuJESDown = 0.5* (MuMuJESDownMet + MuMuJESDownGenMet) 
							MuMuLeptonNoFastSimScaleFactor = 0.5* (MuMuLeptonNoFastSimScaleFactorMet + MuMuLeptonNoFastSimScaleFactorGenMet) 
							MuMuLeptonNoFullSimScaleFactor = 0.5* (MuMuLeptonNoFullSimScaleFactorMet + MuMuLeptonNoFullSimScaleFactorGenMet) 
							MuMuPileupUp = 0.5* (MuMuPileupUpMet + MuMuPileupUpGenMet) 
							MuMuPileupDown = 0.5* (MuMuPileupDownMet + MuMuPileupUpGenMet) 
							MuMuBTagHeavy = 0.5* (MuMuBTagHeavyMet + MuMuBTagHeavyGenMet) 
							MuMuBTagLight = 0.5* (MuMuBTagLightMet + MuMuBTagLightGenMet)
							
							UnscaledYield = MuMusignalYield / (MuMuTriggerEff * scalingLumi)
							
							if UnscaledYield > 0:
								MuMuISRUncertainty = max(abs(MuMuISRUp-UnscaledYield)/UnscaledYield,abs(MuMuISRDown-UnscaledYield)/UnscaledYield) 
								MuMuJESUncertainty = max(abs(MuMuJESUp-UnscaledYield)/UnscaledYield,abs(MuMuJESDown-UnscaledYield)/UnscaledYield) 
								MuMuPileupUncertainty = max(abs(MuMuPileupUp-UnscaledYield)/UnscaledYield,abs(MuMuPileupDown-UnscaledYield)/UnscaledYield) 
								MuMuLeptonFullSimUncertainty = abs(MuMuLeptonNoFullSimScaleFactor-UnscaledYield)/UnscaledYield
								MuMuLeptonFastSimUncertainty = abs(MuMuLeptonNoFastSimScaleFactor-UnscaledYield)/UnscaledYield
								MuMuBTagHeavyUncertainty = abs(MuMuBTagHeavy-UnscaledYield)/UnscaledYield
								MuMuBTagLightUncertainty = abs(MuMuBTagLight-UnscaledYield)/UnscaledYield
								MuMuMetUncertainty = abs(MuMusignalYieldMet-MuMusignalYield)/MuMusignalYield
								MuMuStatUncertainty = sqrt(MuMuMCEvents)
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
							
							
							MuMusystUncertainty = sqrt(MuMuMetUncertainty**2+MuMuBTagHeavyUncertainty**2+MuMuBTagLightUncertainty**2+MuMuJESUncertainty**2+MuMuLeptonFastSimUncertainty**2+MuMuLeptonFullSimUncertainty**2+MuMuPileupUncertainty**2+MuMuISRUncertainty**2+triggerEffUncertainty**2+lumiUncertainty**2)
							
							counts = {}
							counts["MuMu"] = {"MuMuval":MuMusignalYield,"MuMuMCEvents":MuMuMCEvents,"MuMusignalEfficiency":MuMusignalEfficiency,"MuMustatUncertainty":MuMuStatUncertainty,"MuMuTotSystUncertainty":MuMusystUncertainty,
							"MuMuJESUncertainty":MuMuJESUncertainty,"MuMuJESMean":UnscaledYield*MuMuTriggerEff,"MuMuJESUp":MuMuJESUp*MuMuTriggerEff,"MuMuJESDown":MuMuJESDown*MuMuTriggerEff,
							"MuMuLeptonFastSimUncertainty":MuMuLeptonFastSimUncertainty,"MuMuLeptonFastSimMean":UnscaledYield*MuMuTriggerEff,"MuMuLeptonNoFastSimScaleFactor":MuMuLeptonNoFastSimScaleFactor*MuMuTriggerEff,
							"MuMuLeptonFullSimUncertainty":MuMuLeptonFullSimUncertainty,"MuMuLeptonFullSimMean":UnscaledYield*MuMuTriggerEff,"MuMuLeptonNoFullSimScaleFactor":MuMuLeptonNoFullSimScaleFactor*MuMuTriggerEff,
							"MuMupileupUncertainty":MuMuPileupUncertainty,"MuMuPileupMean":UnscaledYield*MuMuTriggerEff,"MuMuPileupUp":MuMuPileupUp*MuMuTriggerEff,"MuMuPileupDown":MuMuPileupDown*MuMuTriggerEff,
							"MuMuISRUncertainty":MuMuISRUncertainty,"MuMuISRMean":UnscaledYield*MuMuTriggerEff,"MuMuISRUp":MuMuISRUp*MuMuTriggerEff,"MuMuISRDown":MuMuISRDown*MuMuTriggerEff,
							"MuMuISRUncertainty":MuMuISRUncertainty,"MuMuISRMean":UnscaledYield*MuMuTriggerEff,"MuMuISRUp":MuMuISRUp*MuMuTriggerEff,"MuMuISRDown":MuMuISRDown*MuMuTriggerEff,
							"MuMuMetUncertainty":MuMuMetUncertainty,"MuMuMet":MuMusignalYieldMet,"MuMuGenMet":MuMusignalYieldGenMet,
							"MuMuBTagHeavyUncertainty":MuMuBTagHeavyUncertainty,"MuMubTagMean":UnscaledYield*MuMuTriggerEff,"MuMubTagHeavy":MuMuBTagHeavy*MuMuTriggerEff,"MuMubTagLightUncertainty":MuMuBTagLightUncertainty,"MuMubTagLight":MuMuBTagLight*MuMuTriggerEff}
							
							outFilePkl = open("shelves/%s_%s_MuMu.pkl"%(fileName,region),"w")
							pickle.dump(counts, outFilePkl)
							outFilePkl.close()
							
					
							
								
				
				
			
			j += 1
		
