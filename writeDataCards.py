#! /usr/bin/env python

### Tool to write datacards from the yields and systematic uncertainties
### determined before

import ROOT
import os
import pickle

import sys
sys.path.append("/.automount/home/home__home4/institut_1b/schomakers/FrameWork/frameWorkBase")
from messageLogger import messageLogger as log
from ROOT import TCanvas, TPad, TH1F, TH2F, TH1I, THStack, TLegend, TMath, TF1
import pickle
from defs import sbottom_masses, getRunRange,theCuts
from math import sqrt

from corrections import rSFOFDirect,rSFOFTrig, rOutIn
from centralConfig import zPredictions, regionsToUse, runRanges,OnlyZPredictions


from optparse import OptionParser


def readPickle(name,regionName,runName,MC=False):
	
	if MC:
		if os.path.isfile("/.automount/home/home__home4/institut_1b/schomakers/FrameWork/frameWorkBase/shelves/%s_%s_%s_MC.pkl"%(name,regionName,runName)):
			result = pickle.load(open("/.automount/home/home__home4/institut_1b/schomakers/FrameWork/frameWorkBase/%s_%s_%s_MC.pkl"%(name,regionName,runName),"rb"))
		else:
			print "/.automount/home/home__home4/institut_1b/schomakers/FrameWork/frameWorkBase/shelves/%s_%s_%s.pkl not found, exiting"%(name,regionName,runName) 		
			sys.exit()		
	else:
		if os.path.isfile("/.automount/home/home__home4/institut_1b/schomakers/FrameWork/frameWorkBase/shelves/%s_%s_%s.pkl"%(name,regionName,runName)):
			result = pickle.load(open("/.automount/home/home__home4/institut_1b/schomakers/FrameWork/frameWorkBase/shelves/%s_%s_%s.pkl"%(name,regionName,runName),"rb"))
		else:
			print "/.automount/home/home__home4/institut_1b/schomakers/FrameWork/frameWorkBase/shelves/%s_%s_%s.pkl not found, exiting"%(name,regionName,runName) 		
			sys.exit()

	return result


### load pickles for the systematics
def loadPickles(path):
	from glob import glob
	result = {}
	for pklPath in glob(path):
		pklFile = open(pklPath, "r")
		result.update(pickle.load(pklFile))
	return result
	
def getWeightedAverage(val1,err1,val2,err2):
	
	weightedAverage = (val1/(err1**2) +val2/(err2**2))/(1./err1**2+1./err2**2)
	weightedAverageErr = 1./(1./err1**2+1./err2**2)**0.5
	
	return weightedAverage, weightedAverageErr
	

### get the results and predictions from the counting experiment in all mass regions
def getResults(shelvesNll,shelvesOnZ):
	
	
	NLLRegions = ["lowNLL","highNLL"]
	massRegions = ["20To60","60To86","96To150","150To200","200To300","300To400","Above400"]	
	MT2Regions = ["highMT2"]
	result = {}
	
	region = "inclusive"
	
	#~ result["rSFOF"] = getattr(rSFOF,region).val
	#~ result["rSFOFErr"] = getattr(rSFOF,region).err
	result["rSFOFDirect"] = getattr(rSFOFDirect,region).val
	result["rSFOFDirectErr"] = getattr(rSFOFDirect,region).err
	result["rSFOFTrig"] = getattr(rSFOFTrig,region).val
	result["rSFOFTrigErr"] = getattr(rSFOFTrig,region).err
	#~ result["rEEOF"] = getattr(rEEOF,region).val
	#~ result["rEEOFErr"] = getattr(rEEOF,region).err
	#~ result["rMMOF"] = getattr(rMMOF,region).val
	#~ result["rMMOFErr"] = getattr(rMMOF,region).err
	
	#~ result["onZPrediction_highNLL"] = shelvesOnZ["86To96_highNll_highMT2_SF"] - shelvesOnZ["86To96_highNll_highMT2_OF"]
	#~ result["onZPrediction_lowNLL"] = shelvesOnZ["86To96_lowNll_highMT2_SF"]	- shelvesOnZ["86To96_lowNll_highMT2_OF"]
	
	result["onZPrediction_highNLL"] = OnlyZPredictions.MT2.SF.highNLL.val
	result["onZPrediction_lowNLL"] = OnlyZPredictions.MT2.SF.lowNLL.val
	result["onZPrediction_highNLL_Err"] = OnlyZPredictions.MT2.SF.highNLL.err
	result["onZPrediction_lowNLL_Err"] = OnlyZPredictions.MT2.SF.lowNLL.err				
	
	
	for selection in NLLRegions:
		result[selection] = {}
		for MT2Region in MT2Regions:
			for massRegion in massRegions:
				if massRegion == "Above400":
					massRegionLabel = "mass400"
				else:
					massRegionLabel = "mass"+massRegion
			
				result[selection]["%s_%s_EE"%(MT2Region,massRegion)] = shelvesNll[selection][getattr(theCuts.mt2Cuts,MT2Region).name+"_"+getattr(theCuts.massCuts,massRegionLabel).name]["EE"]
				result[selection]["%s_%s_MM"%(MT2Region,massRegion)] = shelvesNll[selection][getattr(theCuts.mt2Cuts,MT2Region).name+"_"+getattr(theCuts.massCuts,massRegionLabel).name]["MM"]
				result[selection]["%s_%s_SF"%(MT2Region,massRegion)] = shelvesNll[selection][getattr(theCuts.mt2Cuts,MT2Region).name+"_"+getattr(theCuts.massCuts,massRegionLabel).name]["EE"] + shelvesNll[selection][getattr(theCuts.mt2Cuts,MT2Region).name+"_"+getattr(theCuts.massCuts,massRegionLabel).name]["MM"]
				result[selection]["%s_%s_OF"%(MT2Region,massRegion)] = shelvesNll[selection][getattr(theCuts.mt2Cuts,MT2Region).name+"_"+getattr(theCuts.massCuts,massRegionLabel).name]["EM"]
				result[selection]["%s_%s_OFRMuEScaled"%(MT2Region,massRegion)] = shelvesNll[selection][getattr(theCuts.mt2Cuts,MT2Region).name+"_"+getattr(theCuts.massCuts,massRegionLabel).name]["EMRMuEScaled"]
				result[selection]["%s_%s_OFRMuEScaledUp"%(MT2Region,massRegion)] = shelvesNll[selection][getattr(theCuts.mt2Cuts,MT2Region).name+"_"+getattr(theCuts.massCuts,massRegionLabel).name]["EMRMuEScaledUp"]
				result[selection]["%s_%s_OFRMuEScaledDown"%(MT2Region,massRegion)] = shelvesNll[selection][getattr(theCuts.mt2Cuts,MT2Region).name+"_"+getattr(theCuts.massCuts,massRegionLabel).name]["EMRMuEScaledDown"]
				
				result[selection]["%s_%s_PredFactSF"%(MT2Region,massRegion)] = result[selection]["%s_%s_OFRMuEScaled"%(MT2Region,massRegion)]*getattr(rSFOFTrig,region).val
				if result[selection]["%s_%s_OF"%(MT2Region,massRegion)] > 0:
					result[selection]["%s_%s_PredFactStatErrSF"%(MT2Region,massRegion)] = result[selection]["%s_%s_OF"%(MT2Region,massRegion)]**0.5*result[selection]["%s_%s_PredFactSF"%(MT2Region,massRegion)]/result[selection]["%s_%s_OF"%(MT2Region,massRegion)]
					result[selection]["%s_%s_PredFactSystErrSF"%(MT2Region,massRegion)] = result[selection]["%s_%s_OFRMuEScaled"%(MT2Region,massRegion)]*(getattr(rSFOFTrig,region).err**2 + max(abs(result[selection]["%s_%s_OFRMuEScaled"%(MT2Region,massRegion)] - result[selection]["%s_%s_OFRMuEScaledUp"%(MT2Region,massRegion)])/result[selection]["%s_%s_OFRMuEScaled"%(MT2Region,massRegion)],abs(result[selection]["%s_%s_OFRMuEScaled"%(MT2Region,massRegion)] - result[selection]["%s_%s_OFRMuEScaledDown"%(MT2Region,massRegion)])/result[selection]["%s_%s_OFRMuEScaled"%(MT2Region,massRegion)])**2)**0.5		
				else:
					result[selection]["%s_%s_PredFactStatErrSF"%(MT2Region,massRegion)] = 1.8
					result[selection]["%s_%s_PredFactSystErrSF"%(MT2Region,massRegion)] = 0
				
				if result[selection]["%s_%s_OF"%(MT2Region,massRegion)] > 0:
					result[selection]["%s_%s_RSFOF_Fact"%(MT2Region,massRegion)] = result[selection]["%s_%s_PredFactSF"%(MT2Region,massRegion)] / result[selection]["%s_%s_OF"%(MT2Region,massRegion)]
					result[selection]["%s_%s_RSFOF_Fact_Err"%(MT2Region,massRegion)] = result[selection]["%s_%s_PredFactSystErrSF"%(MT2Region,massRegion)] / result[selection]["%s_%s_OF"%(MT2Region,massRegion)]
				else:
					result[selection]["%s_%s_RSFOF_Fact"%(MT2Region,massRegion)] = 0.
					result[selection]["%s_%s_RSFOF_Fact_Err"%(MT2Region,massRegion)] = 0.
				
				if result[selection]["%s_%s_OF"%(MT2Region,massRegion)] > 0:
					result[selection]["%s_%s_RSFOF_Combined"%(MT2Region,massRegion)],result[selection]["%s_%s_RSFOF_Combined_Err"%(MT2Region,massRegion)] = getWeightedAverage(result[selection]["%s_%s_RSFOF_Fact"%(MT2Region,massRegion)],result[selection]["%s_%s_RSFOF_Fact_Err"%(MT2Region,massRegion)],getattr(rSFOFDirect,region).val,getattr(rSFOFDirect,region).err)
				else:
					result[selection]["%s_%s_RSFOF_Combined"%(MT2Region,massRegion)] = getattr(rSFOFDirect,region).val
					result[selection]["%s_%s_RSFOF_Combined_Err"%(MT2Region,massRegion)] = getattr(rSFOFDirect,region).err
				
				result[selection]["%s_%s_PredSF"%(MT2Region,massRegion)] = result[selection]["%s_%s_OF"%(MT2Region,massRegion)]*result[selection]["%s_%s_RSFOF_Combined"%(MT2Region,massRegion)]
				if result[selection]["%s_%s_OF"%(MT2Region,massRegion)] > 0:
					result[selection]["%s_%s_PredStatErrSF"%(MT2Region,massRegion)] = result[selection]["%s_%s_OF"%(MT2Region,massRegion)]**0.5*result[selection]["%s_%s_RSFOF_Combined"%(MT2Region,massRegion)]
				else:
					result[selection]["%s_%s_PredStatErrSF"%(MT2Region,massRegion)] = 1.8
				result[selection]["%s_%s_PredSystErrSF"%(MT2Region,massRegion)] = result[selection]["%s_%s_OF"%(MT2Region,massRegion)]*result[selection]["%s_%s_RSFOF_Combined_Err"%(MT2Region,massRegion)]
				
				#~ 
				result[selection]["%s_%s_ZPredSF"%(MT2Region,massRegion)] = result["onZPrediction_%s"%selection]*getattr(getattr(rOutIn,massRegionLabel),region).val
				result[selection]["%s_%s_ZPredErrSF"%(MT2Region,massRegion)] = ((result["onZPrediction_%s"%selection]*getattr(getattr(rOutIn,massRegionLabel),region).err)**2 + result["onZPrediction_%s"%selection] * getattr(getattr(rOutIn,massRegionLabel),region).val**2 )**0.5
				
				result[selection]["%s_%s_RarePredSF"%(MT2Region,massRegion)] = shelvesOnZ["%s_%s_SF"%(massRegionLabel,selection)] - shelvesOnZ["%s_%s_OF"%(massRegionLabel,selection)]
				result[selection]["%s_%s_RarePredSF_Up"%(MT2Region,massRegion)] = shelvesOnZ["%s_%s_SF_Up"%(massRegionLabel,selection)] - shelvesOnZ["%s_%s_OF_Up"%(massRegionLabel,selection)]
				result[selection]["%s_%s_RarePredSF_Down"%(MT2Region,massRegion)] = shelvesOnZ["%s_%s_SF_Down"%(massRegionLabel,selection)] - shelvesOnZ["%s_%s_OF_Down"%(massRegionLabel,selection)]
				result[selection]["%s_%s_RarePredErrSF"%(MT2Region,massRegion)] = max(abs(result[selection]["%s_%s_RarePredSF_Up"%(MT2Region,massRegion)]-result[selection]["%s_%s_RarePredSF"%(MT2Region,massRegion)]),abs(result[selection]["%s_%s_RarePredSF_Down"%(MT2Region,massRegion)]-result[selection]["%s_%s_RarePredSF"%(MT2Region,massRegion)]))
		
	
	
	return result
	

			   
def writeDataCards(systematics):
	from messageLogger import messageLogger as log
	from ROOT import TCanvas, TPad, TH1F, TH2F, TH1I, THStack, TLegend, TMath, TF1
	import pickle
	from defs import sbottom_masses
	from math import sqrt
		
	runRange = getRunRange(runRanges.name)
	lumi = runRange.lumi

	path = "shelvesSystematics"	
	
	generalSignalLabel = "T6bbllslepton"
	
	OnZPickle = loadPickles("shelvesMT2/RareOnZ_Powheg.pkl")
	
	countingShelves= {"NLL":readPickle("cutAndCountNLL",regionsToUse.signal.inclusive.name , runRanges.name),"onZ":OnZPickle}
	
	### get the results from the pkls
	results = getResults(countingShelves["NLL"],countingShelves["onZ"])
	#~ results = getResults(countingShelves["NLL"])
	
	
	#~ massRegions = ["LowMass","ZMass","HighMass100To200","HighMass200To400","HighMassHighAbove400"]
	#~ nLLRegions = ["lowNll","highNll"]
	#~ MT2Regions = ["lowMT2","highMT2"]
	massRegions = ["20To60","60To86","96To150","150To200","200To300","300To400","Above400"]
	nLLRegions = ["lowNll","highNll"]
	MT2Regions = ["highMT2"]
	
	signalBins = []

	
	for massRegion in massRegions:
		for nLLRegion in nLLRegions:
			for MT2Region in MT2Regions:
				signalBins.append("%s_%s_%s"%(massRegion,nLLRegion,MT2Region))
	
	

	### mass range
	m_n_min = 150
	m_b_min = 700
	m_b_max = 1600
		
	TriggerEffUncertainty = 0.03
	LumiUncertainty = 0.026
	FastSimUncertainty = 0.04
						
	### loop over mass points
	m_b = m_b_min
	stepsize = 25	
	while m_b <= m_b_max:
		print m_b
		
		M_SBOTTOM = "m_b_"+str(m_b)
		m_sbottom = str(m_b)
		xsection = getattr(sbottom_masses, M_SBOTTOM).cross_section13TeV
		
		m_n = m_n_min
		
		while m_n < m_b:
			
			if m_b < 800:
				stepsize = 25
			else:
				stepsize = 50
	
			m_neutralino_2 = str(m_n)
		

			Pickles = {}
			Yields = {}
			MCEvents = {}
			statUncertainties = {}
			systUncertainties = {}
			JES = {}
			LeptonFastSim = {}
			LeptonFullSim = {}
			Pileup = {}
			ISR= {}				
			Met= {}	
			BTag= {}	
			ScaleShift= {}	
						
			
			### get the yields and uncertainties for each region
			for region in signalBins:
				### get the pickles	
				Pickles["%s_%s"%(m_sbottom,m_neutralino_2)] = loadPickles("%s/%s_msbottom_%s_mneutralino_%s.pkl"%(path,generalSignalLabel,m_sbottom,m_neutralino_2))
				
				### get the yields				
				Yields["EE_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EE"%region]["Val"]
				Yields["EMu_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EMu"%region]["Val"]
				Yields["MuMu_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_MuMu"%region]["Val"]
				Yields["SFOF_%s"%region] = Yields["EE_%s"%region] + Yields["MuMu_%s"%region] - Yields["EMu_%s"%region]
				Yields["SFOF_%s"%region] = max(Yields["SFOF_%s"%region],0)
				
				MCEvents["EE_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EE"%region]["MCEvents"]
				MCEvents["EMu_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EMu"%region]["MCEvents"]
				MCEvents["MuMu_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_MuMu"%region]["MCEvents"]
				MCEvents["SFOF_%s"%region] = MCEvents["EE_%s"%region] + MCEvents["MuMu_%s"%region] - MCEvents["EMu_%s"%region]
				MCEvents["SFOF_%s"%region] = max(MCEvents["SFOF_%s"%region],0)
				
				if MCEvents["SFOF_%s"%region] > 0:
					statUncertainties["SFOF_%s"%region] = sqrt(MCEvents["EE_%s"%region]+MCEvents["EMu_%s"%region]+MCEvents["MuMu_%s"%region])/MCEvents["SFOF_%s"%region]
				else:
					statUncertainties["SFOF_%s"%region] = 0
					
				### JES Uncertainty
				JES["Mean_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EE"%region]["JESMean"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_MuMu"%region]["JESMean"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EMu"%region]["JESMean"] 
				
				JES["JESUp_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EE"%region]["JESUp"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_MuMu"%region]["JESUp"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EMu"%region]["JESUp"] 
				JES["JESDown_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EE"%region]["JESDown"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_MuMu"%region]["JESDown"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EMu"%region]["JESDown"] 
				
				if JES["Mean_%s"%region] > 0:
					systUncertainties["JESUncertainty_%s"%region] = max(abs(JES["JESUp_%s"%region]-JES["Mean_%s"%region])/JES["Mean_%s"%region],abs(JES["JESDown_%s"%region]-JES["Mean_%s"%region])/JES["Mean_%s"%region])
				else:
					systUncertainties["JESUncertainty_%s"%region] = 0
				
				
				
				### Lepton FastSim Uncertainty
				
				LeptonFastSim["Mean_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EE"%region]["LeptonFastSimMean"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_MuMu"%region]["LeptonFastSimMean"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EMu"%region]["LeptonFastSimMean"] 
				LeptonFastSim["MeanShifted_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EE"%region]["LeptonFastSimMean"] * (1+FastSimUncertainty) + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_MuMu"%region]["LeptonFastSimMean"]*(1+FastSimUncertainty)  - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EMu"%region]["LeptonFastSimMean"]*(1+FastSimUncertainty)  
				
				if LeptonFastSim["Mean_%s"%region] > 0:
					systUncertainties["LeptonFastSimUncertainty_%s"%region] = abs(LeptonFastSim["MeanShifted_%s"%region]-LeptonFastSim["Mean_%s"%region])/LeptonFastSim["Mean_%s"%region]
				else:
					systUncertainties["LeptonFastSimUncertainty_%s"%region] = 0
					
				### Lepton FullSim Uncertainty
				
				LeptonFullSim["Mean_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EE"%region]["LeptonFullSimMean"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_MuMu"%region]["LeptonFullSimMean"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EMu"%region]["LeptonFullSimMean"] 
				LeptonFullSim["MeanShifted_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EE"%region]["LeptonFullSimScaleFactorShifted"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_MuMu"%region]["LeptonFullSimScaleFactorShifted"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EMu"%region]["LeptonFullSimScaleFactorShifted"] 
				
				if LeptonFullSim["Mean_%s"%region] > 0:
					systUncertainties["LeptonFullSimUncertainty_%s"%region] = abs(LeptonFullSim["MeanShifted_%s"%region]-LeptonFullSim["Mean_%s"%region])/LeptonFullSim["Mean_%s"%region]
				else:
					systUncertainties["LeptonFullSimUncertainty_%s"%region] = 0
				
				#~ systUncertainties["LeptonFullSimUncertainty_%s"%region] = 0.07
					
				
				###  Pileup Uncertainty
				Pileup["Mean_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EE"%region]["PileupMean"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_MuMu"%region]["PileupMean"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EMu"%region]["PileupMean"] 
				
				Pileup["PileupHigh_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EE"%region]["PileupHigh"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_MuMu"%region]["PileupHigh"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EMu"%region]["PileupHigh"] 
				Pileup["PileupLow_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EE"%region]["PileupLow"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_MuMu"%region]["PileupLow"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EMu"%region]["PileupLow"] 
				
				if Pileup["Mean_%s"%region] > 0:
					systUncertainties["PileupUncertainty_%s"%region] = max(abs(Pileup["PileupHigh_%s"%region]-Pileup["Mean_%s"%region])/Pileup["Mean_%s"%region],abs(Pileup["PileupLow_%s"%region] -Pileup["Mean_%s"%region])/Pileup["Mean_%s"%region] )
				else:
					systUncertainties["PileupUncertainty_%s"%region] = 0
					
				systUncertainties["PileupUncertainty_%s"%region] = 0.02
				
				### ISR Uncertainty
				ISR["Mean_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EE"%region]["ISRMean"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_MuMu"%region]["ISRMean"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EMu"%region]["ISRMean"] 
				
				ISR["ISRUp_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EE"%region]["ISRUp"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_MuMu"%region]["ISRUp"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EMu"%region]["ISRUp"]  
				ISR["ISRDown_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EE"%region]["ISRDown"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_MuMu"%region]["ISRDown"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EMu"%region]["ISRDown"]  
				
				if ISR["Mean_%s"%region] > 0:
					systUncertainties["ISRUncertainty_%s"%region] = max(abs(ISR["ISRUp_%s"%region]-ISR["Mean_%s"%region])/ISR["Mean_%s"%region],abs(ISR["ISRDown_%s"%region]-ISR["Mean_%s"%region])/ISR["Mean_%s"%region])
				else:
					systUncertainties["ISRUncertainty_%s"%region] = 0
					
				### BTag Uncertainty
				BTag["Mean_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EE"%region]["BTagMean"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_MuMu"%region]["BTagMean"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EMu"%region]["BTagMean"] 
				
				if BTag["Mean_%s"%region] > 0:
					BTag["BTagHeavy_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EE"%region]["BTagHeavy"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_MuMu"%region]["BTagHeavy"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EMu"%region]["BTagHeavy"]  
					BTag["BTagLight_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EE"%region]["BTagLight"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_MuMu"%region]["BTagLight"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EMu"%region]["BTagLight"]  
				
					systUncertainties["BTagHeavyUncertainty_%s"%region] = abs(BTag["BTagHeavy_%s"%region]-BTag["Mean_%s"%region])/BTag["Mean_%s"%region]
					systUncertainties["BTagLightUncertainty_%s"%region] = abs(BTag["BTagLight_%s"%region]-BTag["Mean_%s"%region])/BTag["Mean_%s"%region]	
				else:
					systUncertainties["BTagHeavyUncertainty_%s"%region] = 0
					systUncertainties["BTagLightUncertainty_%s"%region] = 0
					
				### FastSim Met Uncertainty					
				Met["Met_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EE"%region]["Met"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_MuMu"%region]["Met"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EMu"%region]["Met"]  
				Met["GenMet_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EE"%region]["GenMet"] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_MuMu"%region]["GenMet"] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EMu"%region]["GenMet"]  
				
				if Yields["SFOF_%s"%region] > 0:
					systUncertainties["MetUncertainty_%s"%region] = 0.5*abs(Met["Met_%s"%region]-Met["GenMet_%s"%region])/Yields["SFOF_%s"%region]
				else:
					systUncertainties["MetUncertainty_%s"%region] = 0
					
				### Trigger eff uncertainty
				if Yields["SFOF_%s"%region] > 0:
					systUncertainties["TriggerEffUncertainty_%s"%region] = sqrt((Yields["EE_%s"%region]*TriggerEffUncertainty)**2 + (Yields["MuMu_%s"%region]*TriggerEffUncertainty)**2 + (Yields["EMu_%s"%region]*TriggerEffUncertainty)**2)/Yields["SFOF_%s"%region]
				else:
					systUncertainties["TriggerEffUncertainty_%s"%region] = 0
					
				### Scale uncertainty
				systUncertainties["ScaleUncertainty_%s"%region] = 0
				if BTag["Mean_%s"%region] > 0:
					
					for scaleIndex in range(1,9):
						ScaleShift["ScaleShift%s"%str(scaleIndex)] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EE"%region]["ScaleShifted%s"%str(scaleIndex)] + Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_MuMu"%region]["ScaleShifted%s"%str(scaleIndex)] - Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["%s_EMu"%region]["ScaleShifted%s"%str(scaleIndex)]
						systUncertainties["ScaleUncertainty_%s"%region] = max(systUncertainties["ScaleUncertainty_%s"%region], abs( (ScaleShift["ScaleShift%s"%str(scaleIndex)] - BTag["Mean_%s"%region])/BTag["Mean_%s"%region]))
				
							
				
				
			### create the datacards for each mass bin
			for massRegion in massRegions:									
				
				### number of signal bins = 2 (ttbar-like, non-ttbar like)
				n_bins = 2 
				
				### number of background processes (FS, Z, and Rare background)
				n_processes = 3
				
				n_nuicance_parameters = 18 
				#~ n_nuicance_parameters = 14 
				
				Name= "T6bbllslepton_%s_%s_%s"%(m_sbottom,m_neutralino_2,massRegion)
				DataCard = open("DataCards/%s.txt"%Name,'w')
				
				### information for us
				DataCard.write("# sbottom = %s \n"%m_sbottom)
				DataCard.write("# neutralino 2 = %s \n"%m_neutralino_2)
				DataCard.write("# mass region = %s \n"%massRegion)
				DataCard.write("# Xsection = %s \n"%xsection)
				
				### bins, processes and nuisance parameters for the tool
				DataCard.write("imax %s number of bins \n"%n_bins)
				DataCard.write("jmax %s number of processes minus 1 \n"%n_processes)
				DataCard.write("kmax %s number of nuisance parameters \n"%n_nuicance_parameters)
				
				### observed event numbers in each bin
				DataCard.write("------------------------------------------------------------------------------------------- \n")
				DataCard.write("bin          LowNll   HighNll  \n")
				DataCard.write("observation  %s       %s       \n"%(results["lowNLL"]["highMT2_%s_SF"%massRegion],results["highNLL"]["highMT2_%s_SF"%massRegion]))
					
				### Background estimation and signal events for each mass point
				### For the flavor symmetric background the OF events * rSFOF have to be put in here
				### without rounding the results. This has to be exactly the same as the numbers put
				### into the lines on the statistic uncertainty for OF background
				DataCard.write("------------------------------------------------------------------------------------------------------------------------------ \n")
				DataCard.write("bin                                         LowNll   	  LowNll        LowNll        LowNll        HighNll         HighNll         HighNll         HighNll \n")
				DataCard.write("process                                     SUSY          FS            onZ           Rare          SUSY            FS              onZ             Rare   \n")
				DataCard.write("process                                     0             1             2             3             0               1               2               3      \n")
				DataCard.write("rate                                        %s            %s            %s            %s            %s              %s              %s              %s      \n"%(Yields["SFOF_%s_lowNll_highMT2"%massRegion],results["lowNLL"]["highMT2_%s_OF"%massRegion]*results["lowNLL"]["highMT2_%s_RSFOF_Combined"%massRegion],results["lowNLL"]["highMT2_%s_ZPredSF"%massRegion],results["lowNLL"]["highMT2_%s_RarePredSF"%massRegion],Yields["SFOF_%s_highNll_highMT2"%massRegion],results["highNLL"]["highMT2_%s_OF"%massRegion]*results["highNLL"]["highMT2_%s_RSFOF_Combined"%massRegion],results["highNLL"]["highMT2_%s_ZPredSF"%massRegion],results["highNLL"]["highMT2_%s_RarePredSF"%massRegion]))
				DataCard.write("------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n")
				
				
				DataCard.write("SigLumiUncertainty               lnN        %s            -            -             -             %s               -               -               -             \n"%(str(1+LumiUncertainty ),str(1+LumiUncertainty )))
				DataCard.write("SigTriggerEffUncertainty         lnN        %s            -            -             -             %s               -               -               -             \n"%(str(1+systUncertainties["TriggerEffUncertainty_%s_lowNll_highMT2"%massRegion] ),str(1+systUncertainties["TriggerEffUncertainty_%s_highNll_highMT2"%massRegion] )))
				DataCard.write("SigJESUncertainty                lnN        %s            -            -             -             %s               -               -               -             \n"%(str(1+systUncertainties["JESUncertainty_%s_lowNll_highMT2"%massRegion] ),str(1+systUncertainties["JESUncertainty_%s_highNll_highMT2"%massRegion] )))
				DataCard.write("SigLeptonFullSimUncertainty      lnN        %s            -            -             -             %s               -               -               -             \n"%(str(1+systUncertainties["LeptonFullSimUncertainty_%s_lowNll_highMT2"%massRegion] ),str(1+systUncertainties["LeptonFullSimUncertainty_%s_highNll_highMT2"%massRegion] )))
				DataCard.write("SigLeptonFastSimUncertainty      lnN        %s            -            -             -             %s               -               -               -             \n"%(str(1+systUncertainties["LeptonFastSimUncertainty_%s_lowNll_highMT2"%massRegion] ),str(1+systUncertainties["LeptonFastSimUncertainty_%s_highNll_highMT2"%massRegion] )))
				DataCard.write("SigPileupUncertainty             lnN        %s            -            -             -             %s               -               -               -             \n"%(str(1+systUncertainties["PileupUncertainty_%s_lowNll_highMT2"%massRegion] ),str(1+systUncertainties["PileupUncertainty_%s_highNll_highMT2"%massRegion] )))
				DataCard.write("SigISRUncertainty                lnN        %s            -            -             -             %s               -               -               -             \n"%(str(1+systUncertainties["ISRUncertainty_%s_lowNll_highMT2"%massRegion] ),str(1+systUncertainties["ISRUncertainty_%s_highNll_highMT2"%massRegion] )))
				DataCard.write("SigMetUncertainty                lnN        %s            -            -             -             %s               -               -               -             \n"%(str(1+systUncertainties["MetUncertainty_%s_lowNll_highMT2"%massRegion] ),str(1+systUncertainties["MetUncertainty_%s_highNll_highMT2"%massRegion] )))
				DataCard.write("SigBTagHeavyUncertainty          lnN        %s            -            -             -             %s               -               -               -             \n"%(str(1+systUncertainties["BTagHeavyUncertainty_%s_lowNll_highMT2"%massRegion] ),str(1+systUncertainties["BTagHeavyUncertainty_%s_highNll_highMT2"%massRegion] )))
				DataCard.write("SigBTagLightUncertainty          lnN        %s            -            -             -             %s               -               -               -             \n"%(str(1+systUncertainties["BTagLightUncertainty_%s_lowNll_highMT2"%massRegion] ),str(1+systUncertainties["BTagLightUncertainty_%s_highNll_highMT2"%massRegion] )))
				DataCard.write("SigScaleUncertainty              lnN        %s            -            -             -             %s               -               -               -             \n"%(str(1+systUncertainties["ScaleUncertainty_%s_lowNll_highMT2"%massRegion] ),str(1+systUncertainties["ScaleUncertainty_%s_highNll_highMT2"%massRegion] )))

				DataCard.write("SigStatUncertLowNll_%s           lnN        %s            -            -             -              -               -               -               -             \n"%(massRegion,str(1+statUncertainties["SFOF_%s_lowNll_highMT2"%massRegion] )))
				DataCard.write("SigStatUncertHighNll_%s          lnN         -            -            -             -             %s               -               -               -             \n"%(massRegion,str(1+statUncertainties["SFOF_%s_highNll_highMT2"%massRegion])))

				### OF yield * RSFOF as prediction -> stat. uncertainty on FS background is gamma uncertainty and depends on the OF statistics
				### not correlated between bins
				DataCard.write("OFStatUncertLowNll_%s             gmN %s     -             %s            -            -            -               -               -               -             \n"%(massRegion,str(results["lowNLL"]["highMT2_%s_OF"%massRegion]),str(results["lowNLL"]["highMT2_%s_RSFOF_Combined"%massRegion])))
				DataCard.write("OFStatUncertHighNll_%s            gmN %s     -             -             -            -            -               %s              -               -             \n"%(massRegion,str(results["highNLL"]["highMT2_%s_OF"%massRegion]),str(results["highNLL"]["highMT2_%s_RSFOF_Combined"%massRegion])))
				
				### Uncertainty on RSFOF -> syst uncertainty on flavor symmetric background. Correlated between bins
				DataCard.write("RSFOFUncert                       lnN        -            %s             -             -             -              %s              -              -             \n"%(str(1+results["lowNLL"]["highMT2_%s_RSFOF_Combined_Err"%massRegion]),str(1+results["highNLL"]["highMT2_%s_RSFOF_Combined_Err"%massRegion])))
				
				DataCard.write("OnZUncert                        lnN        -             -            %s             -              -               -              %s             -             \n"%(str(1+results["lowNLL"]["highMT2_%s_ZPredErrSF"%massRegion]/results["lowNLL"]["highMT2_%s_ZPredSF"%massRegion]),str(1+results["highNLL"]["highMT2_%s_ZPredErrSF"%massRegion]/results["highNLL"]["highMT2_%s_ZPredSF"%massRegion])))
				
				DataCard.write("RareUncert                       lnN        -             -            -             %s             -               -              -              %s            \n"%(str(1+results["lowNLL"]["highMT2_%s_RarePredErrSF"%massRegion]/results["lowNLL"]["highMT2_%s_RarePredSF"%massRegion]),str(1+results["highNLL"]["highMT2_%s_RarePredErrSF"%massRegion]/results["highNLL"]["highMT2_%s_RarePredSF"%massRegion])))
				
						
			m_n += stepsize		
		m_b += stepsize
			
	

		
				
		
	
									

# entry point
#-------------
if (__name__ == "__main__"):
    parser = OptionParser()
    parser.add_option("-s", "--systematics", dest="systematics", action="store_true", default=False,
                                  help="plot maps for systematic uncertainties")
    (opts, args) = parser.parse_args()

   
    writeDataCards(opts.systematics)
 
