
import ROOT
from array import *

def readTreeFromFile(path):
	"""
	helper functionfrom argparse import ArgumentParser
	path: path to .root file containing simulated events
	dileptonCombination: EMu, EMu, or EMu for electron-electron, electron-muon, or muon-muon events

	returns: tree containing events for on sample and dileptonCombination
	"""
	from ROOT import TChain
	result = TChain()
	result.Add("%s/cutsV33DileptonMiniAODSignalNominatorFinalTrees/Tree"%(path))
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
	for filePath in glob("%s/T6bbllslepton*.root"%path):
		sampleName = match(".*T6bbllslepton(.*).root", filePath).groups()[0]
		#for the python enthusiats: yield sampleName, filePath is more efficient here :)
		result[sampleName] = filePath
	return result
	
	
def readTrees(path):
	"""
	path: path to directory containing all sample files
	dileptonCombination: "EMu", "EMu", or pyroot"EMu" for electron-electron, electron-muon, or muon-muon events

	returns: dict of sample names ->  trees containing events (for all samples for one dileptonCombination)
	"""
	result = {}
	for sampleName, filePath in getFilePathsAndSampleNames(path).iteritems():		
		result[sampleName] = readTreeFromFile(filePath)
		
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
	return result




if (__name__ == "__main__"):
	from ROOT import TCanvas, TPad, TH1F, TH1I, THStack, TLegend, TF1, TH2F, TH2D, TFile
	
	#~ firstBinX = 387.5
	#~ lastBinX = 962.5
	#~ nBinsX = 23
	#~ 
	#~ firstBinY = 115.
	#~ lastBinY = 912.5
	#~ nBinsY = 34
	
	nBinsX = 20
	x_bins =  array('d',[387.5, 412.5, 437.5, 462.5, 487.5, 512.5, 537.5, 562.5, 587.5, 612.5, 637.5, 662.5, 687.5, 712.5, 737.5, 762.5, 787.5, 825., 875., 925., 975.])
	nBinsY = 29
	y_bins =  array('d',[137.5,162.5,187.5, 212.5, 237.5, 262.5, 287.5, 312.5, 337.5, 362.5, 387.5, 412.5, 437.5, 462.5, 487.5, 512.5, 537.5, 562.5, 587.5, 612.5, 637.5, 662.5, 687.5, 712.5, 737.5, 775., 825., 875., 925., 975.])

	
	sampleName = "Masses"
	path = "/disk1/user/schomakers/trees/SignalNominator/"
	
	histo2D = TH2F("massScan", "masses", nBinsX, x_bins, nBinsY, y_bins)
	histoISRNormalization = TH2F("ISRNormalization", "ISRNormalization", nBinsX, x_bins, nBinsY, y_bins)
	histoISRNormalizationUp = TH2F("ISRNormalizationUp", "ISRNormalizationUp", nBinsX, x_bins, nBinsY, y_bins)
	histoISRNormalizationDown = TH2F("ISRNormalizationDown", "ISRNormalizationDown", nBinsX, x_bins, nBinsY, y_bins)
	
	f1 = TFile("T6bbllsleptonDenominatorHisto.root","RECREATE")
	
	trees = readTrees(path)
	for sample, tree in trees.iteritems():
		if sample == sampleName:
			m_b_min = 400.
			i = 0
			while m_b_min + i *25. <= 950:
				m_b = m_b_min + i *25.
				m_n_min = 120.
				j = 0
				#~ while m_n_min + j *10. <= 140:
					#~ m_n =  m_n_min + j *10.
					#~ if m_n < m_b:
						#~ cuts = "mSbottom == %s && mNeutralino2 > %s - 5 && mNeutralino2 < %s + 5"%(str(m_b),str(m_n),str(m_n))
						#~ 
						#~ histo = createHistoFromTree(tree, "mSbottom", cuts, 300, 0, 1000)
						#~ events = histo.Integral()
						#~ print "msbottom: "+str(m_b)+", mneutralino: "+str(m_n)+", Events: "+str(events)
						#~ histo2D.SetBinContent(histo2D.GetXaxis().FindBin(m_b),histo2D.GetYaxis().FindBin(m_n),events)
					#~ j += 1
				k = 0
				m_n_min_2 = 150.
				while m_n_min_2 + k *25. <= 900:
					m_n =  m_n_min_2 + k *25.
					if m_n < m_b:
						cuts = "mSbottom == %s && mNeutralino2 > %s - 5 && mNeutralino2 < %s + 5"%(str(m_b),str(m_n),str(m_n))
						cutsISR = "(mSbottom == %s && mNeutralino2 > %s - 5 && mNeutralino2 < %s + 5)*ISRCorrection"%(str(m_b),str(m_n),str(m_n))
						cutsISRUp = "(mSbottom == %s && mNeutralino2 > %s - 5 && mNeutralino2 < %s + 5)*(ISRCorrection+ISRUncertainty)"%(str(m_b),str(m_n),str(m_n))
						cutsISRDown = "(mSbottom == %s && mNeutralino2 > %s - 5 && mNeutralino2 < %s + 5)*(ISRCorrection-ISRUncertainty)"%(str(m_b),str(m_n),str(m_n))
						
						histo = createHistoFromTree(tree, "mSbottom", cuts, 300, 0, 1000)
						histoISR = createHistoFromTree(tree, "mSbottom", cutsISR, 300, 0, 1000)
						histoISRUp = createHistoFromTree(tree, "mSbottom", cutsISRUp, 300, 0, 1000)
						histoISRDown = createHistoFromTree(tree, "mSbottom", cutsISRDown, 300, 0, 1000)
						events = histo.Integral()
						eventsISRWeighted = histoISR.Integral()
						eventsISRWeightedUp = histoISRUp.Integral()
						eventsISRWeightedDown = histoISRDown.Integral()
						print "msbottom: "+str(m_b)+", mneutralino: "+str(m_n)+", Events: "+str(events)
						#~ print "msbottom: "+str(m_b)+", mneutralino: "+str(m_n)+", ISRNormalization: "+str(events/eventsISRWeighted)
						histo2D.SetBinContent(histo2D.GetXaxis().FindBin(m_b),histo2D.GetYaxis().FindBin(m_n),events)
						if events > 0:
							histoISRNormalization.SetBinContent(histo2D.GetXaxis().FindBin(m_b),histo2D.GetYaxis().FindBin(m_n),events/eventsISRWeighted)
							histoISRNormalizationUp.SetBinContent(histo2D.GetXaxis().FindBin(m_b),histo2D.GetYaxis().FindBin(m_n),events/eventsISRWeightedUp)
							histoISRNormalizationDown.SetBinContent(histo2D.GetXaxis().FindBin(m_b),histo2D.GetYaxis().FindBin(m_n),events/eventsISRWeightedDown)
						else:
							histoISRNormalization.SetBinContent(histo2D.GetXaxis().FindBin(m_b),histo2D.GetYaxis().FindBin(m_n),0.)
							histoISRNormalizationUp.SetBinContent(histo2D.GetXaxis().FindBin(m_b),histo2D.GetYaxis().FindBin(m_n),0.)
							histoISRNormalizationDown.SetBinContent(histo2D.GetXaxis().FindBin(m_b),histo2D.GetYaxis().FindBin(m_n),0.)
					
					k += 1
				i += 1
	
	histo2D.Write()
	histoISRNormalization.Write()
	histoISRNormalizationUp.Write()
	histoISRNormalizationDown.Write()
	f1.Close()
	
