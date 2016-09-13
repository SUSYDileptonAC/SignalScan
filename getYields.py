#!/usr/bin/env python

### Routine to get the signal yields for each mass point and signal region

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


### dictionaries for the different cuts that are combined for the signal regions
etaCuts = {
			"central":"abs(eta1) < 1.4 && abs(eta2) < 1.4",
			"forward":"(((abs(eta1) < 1.4 || abs(eta1) > 1.6) && (abs(eta2) < 1.4 || abs(eta2) > 1.6)) && 1.6 <= TMath::Max(abs(eta1),abs(eta2)))",
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

### get the yields for a tree with certain cuts	
def signalYields(tree,cuts):
	histo = createHistoFromTree(tree, "p4.M()", cuts, 50, 0, 500,verbose=False)
	yields = float(histo.Integral(1,histo.GetNbinsX()))
	histo.Delete()
	return yields

		
if (__name__ == "__main__"):
	
	### File to normalize the signal MC
	denominatorFile = TFile("T6bbllsleptonDenominatorHisto.root")
	denominatorHisto = TH2F(denominatorFile.Get("massScan"))

	lumi = 2260
	
	### range of the scan
	m_b_max = 700
	m_b_min = 450
	m_n_max = 650
	m_n_min = 150
	
	m_b_stepsize = 25
			
	path = locations.dataSetPath
	
	### regions that are combined for the signal regions
	etaRegions = ["central","forward"]
	mllRegions = ["lowMass","belowZ","onZ","aboveZ","highMass"]
	bTagBins = ["inclusiveBTags","noBTag","geOneBTag"]
	
	EMutrees = readTrees(path, "EMu")	
	EEtrees = readTrees(path, "EE")	
	MuMutrees = readTrees(path, "MuMu")	
	
	### loop over the mass combinations
	m_b = m_b_min	
	while m_b <= m_b_max:
		print m_b
			
		M_SBOTTOM = "m_b_"+str(m_b)
		m_sbottom = str(m_b)
		m_n_2 = m_n_min
		while m_n_2 <= m_n_max:
			m_neutralino_2 = str(m_n_2)
			
			### step size for the neutralino 2 varies depending on its mass
			if m_n_2 < 300:
				m_n_2_stepsize = 25
			else:
				m_n_2_stepsize = 50
			
			if m_b > m_n_2:
				print "m_n: "+m_neutralino_2
				
				### normalization for this mass point
				denominator = denominatorHisto.GetBinContent(denominatorHisto.GetXaxis().FindBin(m_b),denominatorHisto.GetYaxis().FindBin(m_n_2))
				
				sampleName = "T6bbllslepton_msbottom_%s_mneutralino_%s"%(m_sbottom,m_neutralino_2)
				fileName = "T6bbllslepton_msbottom_%s_mneutralino_%s"%(m_sbottom,m_neutralino_2)
	
				xsection = getattr(sbottom_masses, M_SBOTTOM).cross_section13TeV
				
				scalingLumi = lumi*xsection/denominator

				### loop over signal regions
				for etaRegion in etaRegions:
					etaCut = etaCuts[etaRegion]
					
					### fetch correction factors
					EETriggerEff = getattr(triggerEffs,etaRegion).effEE.val
					EMuTriggerEff = getattr(triggerEffs,etaRegion).effEM.val
					MuMuTriggerEff = getattr(triggerEffs,etaRegion).effMM.val
					RSFOF = getattr(rSFOF,etaRegion).val
					
					
					for mllRegion in mllRegions:
						mllCut = mllCuts[mllRegion]
						
						for bTagBin in bTagBins:
							bTagCut = bTagCuts[bTagBin]	
				
				
							suffix = etaRegion+"_"+mllRegion+"_"+bTagBin
							
							### use cuts without scale factors to get the pure MC event numbers for the statistical uncertainty
							cuts = "weight*leptonFastSimScaleFactor1*leptonFastSimScaleFactor2*leptonFullSimScaleFactor1*leptonFullSimScaleFactor2*bTagWeight*(chargeProduct < 0 && pt1 > 20 && pt2 > 20 && %s && %s && %s && ((met > 100 && nJets >= 3) ||  (met > 150 && nJets >=2)) && abs(eta1) < 2.4 && abs(eta2) < 2.4 && deltaR > 0.3)"%(etaCut,mllCut,bTagCut)
							unweightedCuts = "(chargeProduct < 0 && pt1 > 20 && pt2 > 20 && %s && %s && %s && ((met > 100 && nJets >= 3) ||  (met > 150 && nJets >=2)) && abs(eta1) < 2.4 && abs(eta2) < 2.4 && deltaR > 0.3)"%(etaCut,mllCut,bTagCut)

							counts = {}
							
							### fetch the tree with the right mass point and dilepton combination
							for sample, tree in EEtrees.iteritems():
								if sample == sampleName:
									EEsignalYield = signalYields(tree,cuts) * EETriggerEff * scalingLumi
									EEsignalEfficiency = EEsignalYield/ (lumi*xsection)
									EEMCEvents = signalYields(tree,unweightedCuts)										
									
									counts["EE"] = {"EEval":EEsignalYield,"EEMCEvents":EEMCEvents,"EEsignalEfficiency":EEsignalEfficiency}
									
							
							for sample, tree in EMutrees.iteritems():
								if sample == sampleName:
									EMusignalYield = signalYields(tree,cuts) * EMuTriggerEff * scalingLumi * RSFOF
									EMusignalEfficiency = EMusignalYield/ (lumi*xsection)
									EMuMCEvents = signalYields(tree,unweightedCuts)									
									
									counts["EMu"] = {"EMuval":EMusignalYield,"EMuMCEvents":EMuMCEvents,"EMusignalEfficiency":EMusignalEfficiency}
									
							
							for sample, tree in MuMutrees.iteritems():
								if sample == sampleName:
									MuMusignalYield = signalYields(tree,cuts) * MuMuTriggerEff * scalingLumi
									MuMusignalEfficiency = MuMusignalYield/ (lumi*xsection)
									MuMuMCEvents = signalYields(tree,unweightedCuts)						
									
									counts["MuMu"] = {"MuMuval":MuMusignalYield,"MuMuMCEvents":MuMuMCEvents,"MuMusignalEfficiency":MuMusignalEfficiency}
							
							### store the yields 		
							outFilePkl = open("shelvesYields/%s_%s.pkl"%(fileName,suffix),"w")
							pickle.dump(counts, outFilePkl)
							outFilePkl.close()		
								
			
		
			m_n_2 += m_n_2_stepsize		
		m_b += m_b_stepsize
		
