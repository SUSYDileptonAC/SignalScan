#!/usr/bin/env python

import sys
sys.path.append('cfg/')
from frameworkStructure import pathes
sys.path.append(pathes.basePath)

import pickle	
from ROOT import TCanvas, TPad, TH1F, TH1I, THStack, TLegend, TF1, TH2F, TH2D, TFile, TMath
import ratios
from defs import sbottom_masses
from math import sqrt
from locations import locations
from corrections import triggerEffs, rSFOF
from helpers import readTrees, totalNumberOfGeneratedEvents,createHistoFromTree


etaCuts = {
			"central":"abs(eta1) < 1.4 && abs(eta2) < 1.4",
			"forward":"(((abs(eta1) < 1.4 || abs(eta1) > 1.6) && (abs(eta2) < 1.4 || abs(eta2) > 1.6)) && 1.6 <= TMath::Max(abs(eta1),abs(eta2)))",
			"BothEndcap":"abs(eta1) > 1.6 && abs(eta2) > 1.6",
			"Inclusive":"abs(eta1) < 2.4 && abs(eta2) < 2.4 && ((abs(eta1) < 1.4 || abs(eta1) > 1.6) && (abs(eta2) < 1.4 || abs(eta2) > 1.6))"
			}
			
mllCuts = {
			"lowMass":"p4.M() < 70 && p4.M() > 20",
			"belowZ":"p4.M() < 81 && p4.M() > 70",
			"onZ":"p4.M() < 101 && p4.M() > 81",
			"aboveZ":"p4.M() < 120 && p4.M() > 101",
			"highMass":"p4.M() > 120",
			}
			
bTagCuts = {
			"inclusiveBTags":"nBJets >= 0",
			"noBTag":"nBJets == 0",
			"geOneBTag":"nBJets > 0",
			}

	
def signalYields(tree,cuts):
	histo = createHistoFromTree(tree, "p4.M()", cuts, 50, 0, 500,verbose=False)
	yields = float(histo.Integral(0,histo.GetNbinsX()+1))
	histo.Delete()
	return yields
	

def producePileupUncertainty(tree,cuts):
	
	cuts = cuts.replace("weight*","weightUp*")
	histoUp = createHistoFromTree(tree, "p4.M()", cuts, 50, 0, 500,verbose=False)
	yields = [float(histoUp.Integral())]
	histoUp.Delete()
	
	cuts = cuts.replace("weightUp*","weightDown*")
	histoDown = createHistoFromTree(tree, "p4.M()", cuts, 50, 0, 500,verbose=False)
	yields.append(float(histoDown.Integral()))
	histoDown.Delete()

	return yields[0],yields[1]	
	
def produceJESUncertainty(tree,cuts):

	cuts = cuts.replace("nJets","nShiftedJetsJESUp")	
	cuts = cuts.replace("met","metJESUp")	
	histoUp = createHistoFromTree(tree, "p4.M()", cuts, 50, 0, 500,verbose=False)
	yields = [float(histoUp.Integral(0,histoUp.GetNbinsX()+1))]
	histoUp.Delete()
	
	cuts = cuts.replace("nShiftedJetsJESUp","nShiftedJetsJESDown")	
	cuts = cuts.replace("metJESUp","metJESDown")		
	histoDown = createHistoFromTree(tree, "p4.M()", cuts, 50, 0, 500,verbose=False)
	yields.append(float(histoDown.Integral(0,histoDown.GetNbinsX()+1)))
	histoDown.Delete()
	
	return yields[0],yields[1]

			
def produceISRUncertainty(tree,cuts):
		
	cuts = "(1 + ISRUncertainty)*%s"%cuts
	histoUp = createHistoFromTree(tree, "p4.M()", cuts, 50, 0, 500,verbose=False)
	yields = [float(histoUp.Integral(0,histoUp.GetNbinsX()+1))]
	histoUp.Delete()

	cuts.replace("(1 + ISRUncertainty)","(1 - ISRUncertainty)")
	histoDown = createHistoFromTree(tree, "p4.M()", cuts, 50, 0, 500,verbose=False)
	yields.append(float(histoDown.Integral(0,histoDown.GetNbinsX()+1)))
	histoDown.Delete()
	
	return yields[0],yields[1]
			
def produceLeptonFastSimUncertainty(tree,cuts):
	
	result = 0.
	cuts = cuts.replace("leptonFastSimScaleFactor1*leptonFastSimScaleFactor2*","")
	histo = createHistoFromTree(tree, "p4.M()", cuts, 50, 0, 500,verbose=False)
	result = histo.Integral(0,histo.GetNbinsX()+1)
	histo.Delete()
	
	return result
	
def produceLeptonFullSimUncertainty(tree,cuts):
	
	result = 0.
	cuts = cuts.replace("leptonFullSimScaleFactor1*leptonFullSimScaleFactor2*","")
	histo = createHistoFromTree(tree, "p4.M()", cuts, 50, 0, 500,verbose=False)
	result = histo.Integral(0,histo.GetNbinsX()+1)
	histo.Delete()
	
	return result
	
def produceBTagUncertainty(tree,cuts):
	
	cuts = cuts.replace("bTagWeight","bTagWeight * (1 + bTagWeightErrHeavy)")
	histoUp = createHistoFromTree(tree, "p4.M()", cuts, 50, 0, 500,verbose=False)
	yields = [float(histoUp.Integral(0,histoUp.GetNbinsX()+1))]
	histoUp.Delete()
	
	cuts = cuts.replace("bTagWeight * (1 + bTagWeightErrHeavy)","bTagWeight * (1 + bTagWeightErrLight)")
	histoDown = createHistoFromTree(tree, "p4.M()", cuts, 50, 0, 500,verbose=False)
	yields.append(float(histoDown.Integral(0,histoDown.GetNbinsX()+1)))
	histoDown.Delete()
	
	return yields[0],yields[1]
	

		
if (__name__ == "__main__"):
	
	denominatorFile = TFile("T6bbllsleptonDenominatorHisto.root")
	denominatorHisto = TH2F(denominatorFile.Get("massScan"))

	lumi = 2260
	triggerEffUncertainty = 0.05
	lumiUncertainty = 0.027
		
	m_b_max = 700
	m_b_min = 450
	m_n_max = 650
	m_n_min = 150

		
	path = locations.dataSetPath
	
	etaRegions = ["central","forward"]
	mllRegions = ["lowMass","belowZ","onZ","aboveZ","highMass"]
	bTagBins = ["inclusiveBTags","noBTag","geOneBTag"]
	
	EMutrees = readTrees(path, "EMu")	
	EEtrees = readTrees(path, "EE")	
	MuMutrees = readTrees(path, "MuMu")	
	
	m_b = m_b_min
	m_b_stepsize = 25	
	while m_b <= m_b_max:
		print m_b
			
		M_SBOTTOM = "m_b_"+str(m_b)
		m_sbottom = str(m_b)
		m_n_2 = m_n_min
		while m_n_2 <= m_n_max:
			m_neutralino_2 = str(m_n_2)
			
			if m_n_2 < 300:
				m_n_2_stepsize = 25
			else:
				m_n_2_stepsize = 50
			
			if m_b > m_n_2:
				print "m_n: "+m_neutralino_2
				
				denominator = denominatorHisto.GetBinContent(denominatorHisto.GetXaxis().FindBin(m_b),denominatorHisto.GetYaxis().FindBin(m_n_2))
				
				sampleName = "T6bbllslepton_msbottom_%s_mneutralino_%s"%(m_sbottom,m_neutralino_2)
				fileName = "T6bbllslepton_msbottom_%s_mneutralino_%s"%(m_sbottom,m_neutralino_2)
	
				xsection = getattr(sbottom_masses, M_SBOTTOM).cross_section13TeV
				
				scalingLumi = lumi*xsection/denominator

				for etaRegion in etaRegions:
					etaCut = etaCuts[etaRegion]
					
					if etaRegion == "Barrel":
						EETriggerEff = triggerEffs.central.effEE.val
						EMuTriggerEff = triggerEffs.central.effEM.val
						MuMuTriggerEff = triggerEffs.central.effMM.val
						RSFOF = rSFOF.central.val
					else:
						EETriggerEff = triggerEffs.forward.effEE.val
						EMuTriggerEff = triggerEffs.forward.effEM.val
						MuMuTriggerEff = triggerEffs.forward.effMM.val
						RSFOF = rSFOF.forward.val
					
					for mllRegion in mllRegions:
						mllCut = mllCuts[mllRegion]
						
						for bTagBin in bTagBins:
							bTagCut = bTagCuts[bTagBin]	
				
				
							suffix = etaRegion+"_"+mllRegion+"_"+bTagBin
							cuts = "weight*leptonFastSimScaleFactor1*leptonFastSimScaleFactor2*leptonFullSimScaleFactor1*leptonFullSimScaleFactor2*bTagWeight*(chargeProduct < 0 && pt1 > 20 && pt2 > 20 && %s && %s && %s && ((met > 100 && nJets >= 3) ||  (met > 150 && nJets >=2)) && abs(eta1) < 2.4 && abs(eta2) < 2.4 && deltaR > 0.3)"%(etaCut,mllCut,bTagCut)
							unweightedCuts = "(chargeProduct < 0 && pt1 > 20 && pt2 > 20 && %s && %s && %s && ((met > 100 && nJets >= 3) ||  (met > 150 && nJets >=2)) && abs(eta1) < 2.4 && abs(eta2) < 2.4 && deltaR > 0.3)"%(etaCut,mllCut,bTagCut)

							counts = {}
							
							for sample, tree in EEtrees.iteritems():
								if sample == sampleName:
									EEsignalYield = signalYields(tree,cuts) * EETriggerEff * scalingLumi
									EEsignalEfficiency = EEsignalYield/ (lumi*xsection)
									EEMCEvents = signalYields(tree,unweightedCuts)										
									EEISRUp,EEISRDown = produceISRUncertainty(tree,cuts)									
									EEJESUp,EEJESDown= produceJESUncertainty(tree,cuts)
									EENoLeptonFastSimScaleFactor = produceLeptonFastSimUncertainty(tree,cuts)
									EENoLeptonFullSimScaleFactor = produceLeptonFullSimUncertainty(tree,cuts)
									EEPileupUp,EEPileupDown = producePileupUncertainty(tree,cuts)
									EEBTagHeavy,EEBTagLight = produceBTagUncertainty(tree,cuts)
									
									counts = {}
									counts["EE"] = {"EEval":EEsignalYield,"EEMCEvents":EEMCEvents,"EEsignalEfficiency":EEsignalEfficiency,
									"EEJESUp":EEJESUp * EETriggerEff * scalingLumi,"EEJESDown":EEJESDown * EETriggerEff * scalingLumi,
									"EENoLeptonFastSimScaleFactor":EENoLeptonFastSimScaleFactor * EETriggerEff * scalingLumi,
									"EENoLeptonFullSimScaleFactor":EENoLeptonFullSimScaleFactor * EETriggerEff * scalingLumi,
									"EEPileupUp":EEPileupUp * EETriggerEff * scalingLumi,"EEPileupDown":EEPileupDown * EETriggerEff * scalingLumi,
									"EEISRUp":EEISRUp * EETriggerEff * scalingLumi,"EEISRDown":EEISRDown * EETriggerEff * scalingLumi,
									"EEbTagHeavy":EEBTagHeavy * EETriggerEff * scalingLumi,"EEbTagLight":EEBTagLight * EETriggerEff * scalingLumi}
									
							
							for sample, tree in EMutrees.iteritems():
								if sample == sampleName:
									EMusignalYield = signalYields(tree,cuts) * EMuTriggerEff * scalingLumi * RSFOF
									EMusignalEfficiency = EMusignalYield/ (lumi*xsection)
									EMuMCEvents = signalYields(tree,unweightedCuts)									
									EMuISRUp,EMuISRDown = produceISRUncertainty(tree,cuts)									
									EMuJESUp,EMuJESDown= produceJESUncertainty(tree,cuts)
									EMuNoLeptonFastSimScaleFactor = produceLeptonFastSimUncertainty(tree,cuts)
									EMuNoLeptonFullSimScaleFactor = produceLeptonFullSimUncertainty(tree,cuts)
									EMuPileupUp,EMuPileupDown = producePileupUncertainty(tree,cuts)
									EMuBTagHeavy,EMuBTagLight = produceBTagUncertainty(tree,cuts)
						
									counts["EMu"] = {"EMuval":EMusignalYield,"EMuMCEvents":EMuMCEvents,"EMusignalEfficiency":EMusignalEfficiency,
									"EMuJESUp":EMuJESUp * EMuTriggerEff * scalingLumi * RSFOF,"EMuJESDown":EMuJESDown * EMuTriggerEff * scalingLumi * RSFOF,
									"EMuNoLeptonFastSimScaleFactor":EMuNoLeptonFastSimScaleFactor * EMuTriggerEff * scalingLumi * RSFOF,
									"EMuNoLeptonFullSimScaleFactor":EMuNoLeptonFullSimScaleFactor * EMuTriggerEff * scalingLumi * RSFOF,
									"EMuPileupUp":EMuPileupUp * EMuTriggerEff * scalingLumi * RSFOF,"EMuPileupDown":EMuPileupDown * EMuTriggerEff * scalingLumi * RSFOF,
									"EMuISRUp":EMuISRUp * EMuTriggerEff * scalingLumi * RSFOF,"EMuISRDown":EMuISRDown * EMuTriggerEff * scalingLumi * RSFOF,
									"EMubTagHeavy":EMuBTagHeavy * EMuTriggerEff * scalingLumi * RSFOF,"EMubTagLight":EMuBTagLight * EMuTriggerEff * scalingLumi * RSFOF}
				
							
							for sample, tree in MuMutrees.iteritems():
								if sample == sampleName:
									MuMusignalYield = signalYields(tree,cuts) * MuMuTriggerEff * scalingLumi
									MuMusignalEfficiency = MuMusignalYield/ (lumi*xsection)
									MuMuMCEvents = signalYields(tree,unweightedCuts)										
									MuMuISRUp,MuMuISRDown = produceISRUncertainty(tree,cuts)									
									MuMuJESUp,MuMuJESDown= produceJESUncertainty(tree,cuts)
									MuMuNoLeptonFastSimScaleFactor = produceLeptonFastSimUncertainty(tree,cuts)
									MuMuNoLeptonFullSimScaleFactor = produceLeptonFullSimUncertainty(tree,cuts)
									MuMuPileupUp,MuMuPileupDown = producePileupUncertainty(tree,cuts)
									MuMuBTagHeavy,MuMuBTagLight = produceBTagUncertainty(tree,cuts)
									
									counts["MuMu"] = {"MuMuval":MuMusignalYield,"MuMuMCEvents":MuMuMCEvents,"MuMusignalEfficiency":MuMusignalEfficiency,
									"MuMuJESUp":MuMuJESUp * MuMuTriggerEff * scalingLumi,"MuMuJESDown":MuMuJESDown * MuMuTriggerEff * scalingLumi,
									"MuMuNoLeptonFastSimScaleFactor":MuMuNoLeptonFastSimScaleFactor * MuMuTriggerEff * scalingLumi,
									"MuMuNoLeptonFullSimScaleFactor":MuMuNoLeptonFullSimScaleFactor * MuMuTriggerEff * scalingLumi,
									"MuMuPileupUp":MuMuPileupUp * MuMuTriggerEff * scalingLumi,"MuMuPileupDown":MuMuPileupDown * MuMuTriggerEff * scalingLumi,
									"MuMuISRUp":MuMuISRUp * MuMuTriggerEff * scalingLumi,"MuMuISRDown":MuMuISRDown * MuMuTriggerEff * scalingLumi,
									"MuMubTagHeavy":MuMuBTagHeavy * MuMuTriggerEff * scalingLumi,"MuMubTagLight":MuMuBTagLight * MuMuTriggerEff * scalingLumi}
									
							outFilePkl = open("shelves/%s_%s.pkl"%(fileName,suffix),"w")
							pickle.dump(counts, outFilePkl)
							outFilePkl.close()
									
										
								
			
		
			m_n_2 += m_n_2_stepsize		
		m_b += m_b_stepsize
		
