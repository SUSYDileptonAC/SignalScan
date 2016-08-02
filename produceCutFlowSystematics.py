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

sleptonPoints = {
				"1":{"m_b":700,"m_n":175},
				"2":{"m_b":700,"m_n":500},
				"3":{"m_b":750,"m_n":300},
				"4":{"m_b":800,"m_n":200},
				"5":{"m_b":800,"m_n":600},
				}

#~ mllCuts = {
			#~ "noMllCut":"1>0",
			#~ "lowMll":"mll < 70 && mll > 20",
			#~ "highMll":"mll > 120",
			#~ }
#~ 
mllCuts = {
			"noMllCut":"1>0",
			"lowMll":"p4.M() < 70 && p4.M() > 20",
			"highMll":"p4.M() > 120",
			}

NllCuts = {
			"inclusiveNll":"1 > 0",
			"lowNll":"nLL < 21",
			"highNll":"nLL >= 21",
			}

signalRegionCuts = {
			"basicCut":"p4.M() > 20",
			"signalRegionCut":"nJets >= 2 && met > 150",
			"met150Cut":"met > 150",
			"nJets2Cut":"nJets >= 2",
			}			
			

def readPDFHistsFromFile(filePath, dileptonCombination, pdfSet):
	from ROOT import TH1F, TFile
	from random import randint
	from sys import maxint
	name1 = "%x"%(randint(0, maxint))	
	name2 = "%x"%(randint(0, maxint))	
	name3 = "%x"%(randint(0, maxint))	
	result = {}
	f1 = TFile(filePath,"READ")
	
	firstBin = 20
	lastBin = 300
	nBins = 280/5
	histoMean = TH1F(name1,name1, nBins, firstBin, lastBin)
	histoMean = (f1.Get("%s_%sDileptonTree_%s"%(pdfSet,dileptonCombination,"mean"))).Clone()
	result["mean"] = histoMean.Clone()
	result["up"] = f1.Get("%s_%sDileptonTree_%s"%(pdfSet,dileptonCombination,"up")).Clone(name2)
	result["down"] = f1.Get("%s_%sDileptonTree_%s"%(pdfSet,dileptonCombination,"down")).Clone(name3)
	f1.Close()
	print result
	return result

def readPDFHists(path,dileptonCombination):
	
	result = {}
	for sampleName, filePath in getFilePathsAndSampleNames(path).iteritems():
		result[sampleName] = {}
		print readPDFHistsFromFile(filePath, dileptonCombination,"CT10")
		result[sampleName]["CT10"] = readPDFHistsFromFile(filePath, dileptonCombination,"CT10")
		result[sampleName]["MSTW"] = readPDFHistsFromFile(filePath, dileptonCombination,"MSTW")
		result[sampleName]["NNPDF"] = readPDFHistsFromFile(filePath, dileptonCombination,"NNPDF")
		
	return result	

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
		#~ sampleName = match(".*sw538v.*\.processed.*\.(.*).root", filePath).groups()[0]
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
	#~ histo = createHistoFromTree(tree, "mll", cuts, 300, 0, 300)
	histo = createHistoFromTree(tree, "bTagWeight", cuts, 400, 0, 40)
	histo.Scale(scalingLumi)
	yields = float(histo.Integral())
	return yields
			
def produceSignalEfficiency(tree,cuts,weight):
	#~ yields = [float(createHistoFromTree(tree, "mll", cuts, 300, 0, 300).Integral())]
	yields = [float(createHistoFromTree(tree, "bTagWeight", cuts, 400, 0, 4).Integral())]
	result = yields[0]/weight
	return result
			
def statisticalUncertainty(tree,cuts):
	#~ histo = createHistoFromTree(tree, "mll", cuts, 300, 0, 300)
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
		#~ yields.append(float(createHistoFromTree(tree, "mll", cuts, 300, 0, 300).Integral()))
		yields.append(float(createHistoFromTree(tree, "bTagWeight", cuts, 400, 0, 4).Integral()))

	return yields[0],yields[1]	
	
def produceJESUncertainty(tree,cuts):

	cuts = cuts.replace("nJets","nShiftedJetsJESUp")	
	cuts = cuts.replace("met","metJESUp")	
	cuts = cuts.replace("nLL","nLLJESUp")	
	#~ yields.append(float(createHistoFromTree(tree, "mll", cuts, 300, 0, 300).Integral()))
	yields = [float(createHistoFromTree(tree, "bTagWeight", cuts, 400, 0, 4).Integral())]
	
	cuts = cuts.replace("nShiftedJetsJESUp","nShiftedJetsJESDown")	
	cuts = cuts.replace("metJESUp","metJESDown")	
	cuts = cuts.replace("nLLJESUp","nLLJESDown")	
	#~ yields.append(float(createHistoFromTree(tree, "mll", cuts, 300, 0, 300).Integral()))
	yields.append(float(createHistoFromTree(tree, "bTagWeight", cuts, 400, 0, 4).Integral()))


	return yields[0],yields[1]

			
def produceISRUncertainty(tree,cuts):
	
	cuts = cuts.replace("ISRCorrection","(ISRCorrection+ISRUncertainty)")
	#~ yields.append(float(createHistoFromTree(tree, "mll", cuts, 300, 0, 300).Integral()))
	yields = [float(createHistoFromTree(tree, "bTagWeight", cuts, 400, 0, 4).Integral())]

	cuts = cuts.replace("(ISRCorrection+ISRUncertainty)","(ISRCorrection-ISRUncertainty)")
	#~ yields.append(float(createHistoFromTree(tree, "mll", cuts, 300, 0, 300).Integral()))
	yields.append(float(createHistoFromTree(tree, "bTagWeight", cuts, 400, 0, 4).Integral()))
	
	return yields[0],yields[1]
			
def produceLeptonFastSimUncertainty(tree,cuts):
	
	result = 0.
	cuts = cuts.replace("leptonFastSimScaleFactor1*leptonFastSimScaleFactor2*","")
	#~ yields.append(float(createHistoFromTree(tree, "mll", cuts, 300, 0, 300).Integral()))
	result = float(createHistoFromTree(tree, "bTagWeight", cuts, 400, 0, 4).Integral())
	
	return result
	
def produceLeptonFullSimUncertainty(tree,cuts):
	
	result = 0.
	cuts = cuts.replace("leptonFullSimScaleFactor1*leptonFullSimScaleFactor2*","")
	#~ yields.append(float(createHistoFromTree(tree, "mll", cuts, 300, 0, 300).Integral()))
	result = float(createHistoFromTree(tree, "bTagWeight", cuts, 400, 0, 4).Integral())
	
	return result
	
def produceBTagUncertainty(tree,cuts):
	
	cuts = cuts.replace("bTagWeight","bTagWeight * (1 + bTagWeightErrHeavy)")
	#~ yields.append(float(createHistoFromTree(tree, "mll", cuts, 300, 0, 300).Integral()))
	yields = [float(createHistoFromTree(tree, "bTagWeight", cuts, 400, 0, 4).Integral())]
	
	cuts = cuts.replace("bTagWeight * (1 + bTagWeightErrHeavy)","bTagWeight * (1 + bTagWeightErrLight)")
	#~ yields.append(float(createHistoFromTree(tree, "mll", cuts, 300, 0, 300).Integral()))
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
		
	#~ path = "/disk1/user/schomakers/trees/sw8014v1009_NLL"
	path = "/disk1/user/schomakers/trees/sw8014v1009"
	
	#~ CutRegions = ["basicCut","signalRegionCut","met150Cut","nJets2Cut"]
	CutRegions = ["basicCut","met150Cut","nJets2Cut"]
	#~ CutRegions = ["signalRegionCut"]
	mllCutRegions = ["noMllCut","lowMll","highMll"]
	#~ NllCutRegions = ["inclusiveNll","lowNll","highNll"]
	NllCutRegions = ["inclusiveNll"]
	
	points = ["1","2","3","4","5"]
	
	BR_factor =  0.25      + 0.5   + 0.25*(0.13+0.06*0.124+0.0011*(0.015+0.113+0.312))

	for point in points:
		m_b = sleptonPoints[point]["m_b"]
		m_sbottom = str(m_b)
		M_SBOTTOM = "m_b_"+str(m_b)
		m_n_2 =sleptonPoints[point]["m_n"]
		m_neutralino_2 = str(m_n_2)
		
		xsection = getattr(sbottom_masses, M_SBOTTOM).cross_section13TeV
		
		sampleName = "T6bbllslepton_msbottom_%s_mneutralino_%s"%(m_sbottom,m_neutralino_2)
		fileName = "T6bbllslepton_msbottom_%s_mneutralino_%s"%(m_sbottom,m_neutralino_2)
		
		triggerEffUncertainty = 0.05
		lumiUncertainty = 0.062
		electronUncertainty = 0.05
		muonUncertainty = 0.03
		
		denominator = denominatorHisto.GetBinContent(denominatorHisto.GetXaxis().FindBin(m_b),denominatorHisto.GetYaxis().FindBin(m_n_2))
		ISRNormalization = ISRNormalizationHisto.GetBinContent(ISRNormalizationHisto.GetXaxis().FindBin(m_b),ISRNormalizationHisto.GetYaxis().FindBin(m_n_2))
		ISRNormalizationUp = ISRNormalizationHistoUp.GetBinContent(ISRNormalizationHisto.GetXaxis().FindBin(m_b),ISRNormalizationHisto.GetYaxis().FindBin(m_n_2))
		ISRNormalizationDown = ISRNormalizationHistoDown.GetBinContent(ISRNormalizationHisto.GetXaxis().FindBin(m_b),ISRNormalizationHisto.GetYaxis().FindBin(m_n_2))
		
		totalEventYield = lumi*xsection*BR_factor
		totalEventYieldUncert = lumi*xsection*BR_factor/sqrt(denominator)
		print totalEventYield
		print totalEventYieldUncert
	
		for region in CutRegions:
			signalRegionCut = signalRegionCuts[region]
			 
			for mllCutRegion in mllCutRegions:
				mllCut = mllCuts[mllCutRegion]
				
				for NllCutRegion in NllCutRegions:
					NllCut = NllCuts[NllCutRegion]
			
			
					cuts = "weight*ISRCorrection*leptonFastSimScaleFactor1*leptonFastSimScaleFactor2*leptonFullSimScaleFactor1*leptonFullSimScaleFactor2*bTagWeight*( nUnmatchedJets == 0 && chargeProduct < 0 && ((pt1 > 25 && pt2 > 20) || (pt1 > 20 && pt2 > 25)) && abs(eta1) < 2.4 && abs(eta2) < 2.4 && deltaR > 0.1 && %s && %s && %s)"%(signalRegionCut,mllCut,NllCut)
					unweightedCuts = "( nUnmatchedJets == 0 && chargeProduct < 0 && ((pt1 > 25 && pt2 > 20) || (pt1 > 20 && pt2 > 25)) && abs(eta1) < 2.4 && abs(eta2) < 2.4 && deltaR > 0.1 && %s && %s && %s)"%(signalRegionCut,mllCut,NllCut)
					
					
					
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
			
	
					EEtrees = readTrees(path, "EE")	
					MuMutrees = readTrees(path, "MuMu")	
					
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
							
							UnscaledYieldEE = EEsignalYield / (EETriggerEff * scalingLumi)
							
							if UnscaledYieldEE > 0:
								EEISRUncertainty = max(abs(EEISRUp-UnscaledYieldEE)/UnscaledYieldEE,abs(EEISRDown-UnscaledYieldEE)/UnscaledYieldEE) 
								EEJESUncertainty = max(abs(EEJESUp-UnscaledYieldEE)/UnscaledYieldEE,abs(EEJESDown-UnscaledYieldEE)/UnscaledYieldEE) 
								EEPileupUncertainty = max(abs(EEPileupUp-UnscaledYieldEE)/UnscaledYieldEE,abs(EEPileupDown-UnscaledYieldEE)/UnscaledYieldEE) 
								EELeptonFullSimUncertainty = abs(EELeptonNoFullSimScaleFactor-UnscaledYieldEE)/UnscaledYieldEE
								EELeptonFastSimUncertainty = abs(EELeptonNoFastSimScaleFactor-UnscaledYieldEE)/UnscaledYieldEE
								EEBTagHeavyUncertainty = abs(EEBTagHeavy-UnscaledYieldEE)/UnscaledYieldEE
								EEBTagLightUncertainty = abs(EEBTagLight-UnscaledYieldEE)/UnscaledYieldEE
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
							
							UnscaledYieldMuMu = MuMusignalYield / (MuMuTriggerEff * scalingLumi)
							
							if UnscaledYieldMuMu > 0:
								MuMuISRUncertainty = max(abs(MuMuISRUp-UnscaledYieldMuMu)/UnscaledYieldMuMu,abs(MuMuISRDown-UnscaledYieldMuMu)/UnscaledYieldMuMu) 
								MuMuJESUncertainty = max(abs(MuMuJESUp-UnscaledYieldMuMu)/UnscaledYieldMuMu,abs(MuMuJESDown-UnscaledYieldMuMu)/UnscaledYieldMuMu) 
								MuMuPileupUncertainty = max(abs(MuMuPileupUp-UnscaledYieldMuMu)/UnscaledYieldMuMu,abs(MuMuPileupDown-UnscaledYieldMuMu)/UnscaledYieldMuMu) 
								MuMuLeptonFullSimUncertainty = abs(MuMuLeptonNoFullSimScaleFactor-UnscaledYieldMuMu)/UnscaledYieldMuMu
								MuMuLeptonFastSimUncertainty = abs(MuMuLeptonNoFastSimScaleFactor-UnscaledYieldMuMu)/UnscaledYieldMuMu
								MuMuBTagHeavyUncertainty = abs(MuMuBTagHeavy-UnscaledYieldMuMu)/UnscaledYieldMuMu
								MuMuBTagLightUncertainty = abs(MuMuBTagLight-UnscaledYieldMuMu)/UnscaledYieldMuMu
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
							
							
							
					counts = {}
					SFYield = EEsignalYield + MuMusignalYield
					if SFYield > 0:
						if EEMCEvents+MuMuMCEvents > 0:
							StatUncertainty = sqrt(EEMCEvents+MuMuMCEvents)/(EEMCEvents+MuMuMCEvents)
						else:
							StatUncertainty = 0
						
						Mean = UnscaledYieldEE*EETriggerEff + UnscaledYieldMuMu*MuMuTriggerEff
						### JES Uncertainty
						JESUp = EEJESUp*EETriggerEff + MuMuJESUp*MuMuTriggerEff
						JESDown = EEJESDown*EETriggerEff + MuMuJESDown*MuMuTriggerEff
						
						if Mean > 0:
							JESUncertainty = max(abs(JESUp-Mean)/Mean,abs(JESDown-Mean)/Mean)
						else:
							JESUncertainty = 0
						
						METYield = 	EEsignalYieldMet + MuMusignalYieldMet
						GenMETYield = 	EEsignalYieldGenMet + MuMusignalYieldGenMet
						METUncertainty = 0.5*abs(METYield-GenMETYield)/SFYield
						
						
						
						### Lepton FastSim Uncertainty						
						LeptonNoFastSimScaleFactor = EELeptonNoFastSimScaleFactor*EETriggerEff + MuMuLeptonNoFastSimScaleFactor*MuMuTriggerEff
						
						if Mean > 0:
							LeptonFastSimUncertainty = abs(LeptonNoFastSimScaleFactor-Mean)/Mean
						else:
							LeptonFastSimUncertainty = 0
						
						### Lepton FullSim Uncertainty						
						LeptonNoFullSimScaleFactor = EELeptonNoFullSimScaleFactor*EETriggerEff + MuMuLeptonNoFullSimScaleFactor*MuMuTriggerEff
						
						if Mean > 0:
							LeptonFullSimUncertainty = abs(LeptonNoFullSimScaleFactor-Mean)/Mean
						else:
							LeptonFullSimUncertainty = 0
						
						###  Pileup Uncertainty
						
						PileupUp = EEPileupUp*EETriggerEff + MuMuPileupUp*MuMuTriggerEff
						PileupDown = EEPileupDown*EETriggerEff + MuMuPileupDown*MuMuTriggerEff
						
						if Mean > 0:
							PileupUncertainty = max(abs(PileupUp-Mean)/Mean,abs(PileupDown-Mean)/Mean)
						else:
							PileupUncertainty = 0
						
						### ISR Uncertainty
						ISRUp = EEISRUp*EETriggerEff + MuMuISRUp*MuMuTriggerEff
						ISRDown = EEISRDown*EETriggerEff + MuMuISRDown*MuMuTriggerEff
						
						if Mean > 0:
							ISRUncertainty = max(abs(ISRUp-Mean)/Mean,abs(ISRDown-Mean)/Mean)
						else:
							ISRUncertainty = 0
							
						### BTag Uncertainty
						BTagHeavy = EEBTagHeavy*EETriggerEff + MuMuBTagHeavy*MuMuTriggerEff
						
						BTagLight = EEBTagLight*EETriggerEff + MuMuBTagLight*MuMuTriggerEff
						
						if Mean > 0:
							BTagUncertaintyHeavy = abs(BTagHeavy-Mean)/Mean
							BTagUncertaintyLight = abs(BTagLight-Mean)/Mean
						else:
							BTagUncertaintyHeavy = 0
							BTagUncertaintyLight = 0
							
						SystUncertainty = sqrt(JESUncertainty**2 + METUncertainty**2 + LeptonFastSimUncertainty**2 + LeptonFullSimUncertainty**2 + PileupUncertainty**2 + ISRUncertainty**2 + BTagUncertaintyLight**2 + BTagUncertaintyHeavy**2  + triggerEffUncertainty**2 + lumiUncertainty**2)

					else:
						SFYield = 0 
						StatUncertainty = 0
						SystUncertainty = 0
				
					counts["Signal"] = {"totalEventYield":totalEventYield,"totalEventYieldUncert":totalEventYieldUncert,"SFYield":SFYield,"StatUncertainty":StatUncertainty,"SystUncertainty":SystUncertainty}		
					outFilePkl = open("shelvesCutFlowTables/%s_%s_%s_%s.pkl"%(fileName,region,mllCutRegion,NllCutRegion),"w")
					pickle.dump(counts, outFilePkl)
					outFilePkl.close()
										
		
