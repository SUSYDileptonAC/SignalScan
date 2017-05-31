
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
	
	nBinsX = 21
	x_bins =  array('d',[ 687.5, 712.5, 737.5, 762.5, 787.5, 825., 875., 925., 975., 1025., 1075., 1125., 1175., 1225., 1275., 1325., 1375., 1425., 1475., 1525., 1575., 1625.])
	#~ x_bins =  array('d',[ 687.5, 712.5, 737.5, 762.5, 787.5, 825., 875., 925., 975., 1025., 1075., 1125., 1175., 1225., 1275., 1325., 1375., 1425., 1475., 1525.])
	nBinsY = 41
	y_bins =  array('d',[137.5,162.5,187.5, 212.5, 237.5, 262.5, 287.5, 312.5, 337.5, 362.5, 387.5, 412.5, 437.5, 462.5, 487.5, 512.5, 537.5, 562.5, 587.5, 612.5, 637.5, 662.5, 687.5, 712.5, 737.5, 775., 825., 875., 925., 975., 1025., 1075., 1125., 1175., 1225., 1275., 1325., 1375., 1425., 1475., 1525., 1575.])
	#~ y_bins =  array('d',[137.5,162.5,187.5, 212.5, 237.5, 262.5, 287.5, 312.5, 337.5, 362.5, 387.5, 412.5, 437.5, 462.5, 487.5, 512.5, 537.5, 562.5, 587.5, 612.5, 637.5, 662.5, 687.5, 712.5, 737.5, 775., 825., 875., 925., 975., 1025., 1075., 1125., 1175., 1225., 1275., 1325., 1375., 1425., 1475.])

	
	sampleName1 = "Nominator1"
	sampleName2 = "Nominator2"
	path = "/disk1/user/schomakers/trees/SignalNominator/"
	
	histo2D = TH2F("massScan", "masses", nBinsX, x_bins, nBinsY, y_bins)
	histo2DHighPU = TH2F("massScanHighPU", "massesHighPU", nBinsX, x_bins, nBinsY, y_bins)
	histo2DLowPU = TH2F("massScanLowPU", "massesLowPU", nBinsX, x_bins, nBinsY, y_bins)
	histoISRNormalization = TH2F("ISRNormalization", "ISRNormalization", nBinsX, x_bins, nBinsY, y_bins)
	histoISRNormalizationUp = TH2F("ISRNormalizationUp", "ISRNormalizationUp", nBinsX, x_bins, nBinsY, y_bins)
	histoISRNormalizationDown = TH2F("ISRNormalizationDown", "ISRNormalizationDown", nBinsX, x_bins, nBinsY, y_bins)
	
	histosScaleWeights = {}
	for i in range (1,9):
		histosScaleWeights["ScaleWeight"+str(i)] = TH2F("ScaleWeight"+str(i), "ScaleWeight"+str(i), nBinsX, x_bins, nBinsY, y_bins)
	
	f1 = TFile("T6bbllsleptonDenominatorHisto7.root","RECREATE")
	
	m_b_min = 700
	m_b_max = 1500
	m_n_min = 150
	
	trees = readTrees(path)
	for sample, tree in trees.iteritems():
		if sample == sampleName1:
			m_b = m_b_min
			
			while m_b <= m_b_max:
				print m_b
				if m_b < 800:
					stepsize = 25
				else:
					stepsize = 50
			
				m_n = m_n_min
				
				while m_n < m_b:
					
					if not ( m_b > 1500 or (m_b == 775 and m_n == 750) or (m_b == 800 and m_n == 150) or (m_b == 950 and m_n == 250) or (m_b == 950 and m_n == 300) or (m_b == 950 and m_n == 500) or (m_b == 950 and m_n == 550) or (m_b == 950 and m_n == 850) or (m_b == 950 and m_n == 900)):	
					
						cuts = "mSbottom == %s && mNeutralino2 > %s - 5 && mNeutralino2 < %s + 5"%(str(m_b),str(m_n),str(m_n))
						cutsISR = "(mSbottom == %s && mNeutralino2 > %s - 5 && mNeutralino2 < %s + 5)*ISRCorrection"%(str(m_b),str(m_n),str(m_n))
						cutsISRUp = "(mSbottom == %s && mNeutralino2 > %s - 5 && mNeutralino2 < %s + 5)*(ISRCorrection+ISRUncertainty)"%(str(m_b),str(m_n),str(m_n))
						cutsISRDown = "(mSbottom == %s && mNeutralino2 > %s - 5 && mNeutralino2 < %s + 5)*(ISRCorrection-ISRUncertainty)"%(str(m_b),str(m_n),str(m_n))
						
						cutsHighPU = "mSbottom == %s && mNeutralino2 > %s - 5 && mNeutralino2 < %s + 5 && nVertices > 16"%(str(m_b),str(m_n),str(m_n))
						cutsLowPU = "mSbottom == %s && mNeutralino2 > %s - 5 && mNeutralino2 < %s + 5 && nVertices <= 16"%(str(m_b),str(m_n),str(m_n))
						
						
						
						histo = createHistoFromTree(tree, "mSbottom", cuts, 50, 0, 1700)
						histoISR = createHistoFromTree(tree, "mSbottom", cutsISR, 50, 0, 1700)
						histoISRUp = createHistoFromTree(tree, "mSbottom", cutsISRUp, 50, 0, 1700)
						histoISRDown = createHistoFromTree(tree, "mSbottom", cutsISRDown, 50, 0, 1700)
						histoHighPU = createHistoFromTree(tree, "mSbottom", cutsHighPU, 50, 0, 1700)
						histoLowPU = createHistoFromTree(tree, "mSbottom", cutsLowPU, 50, 0, 1700)
						
												
						events = histo.Integral()
						eventsISRWeighted = histoISR.Integral()
						eventsISRWeightedUp = histoISRUp.Integral()
						eventsISRWeightedDown = histoISRDown.Integral()
						eventsISRHighPU = histoHighPU.Integral()
						eventsISRLowPU = histoLowPU.Integral()
							
						print "msbottom: "+str(m_b)+", mneutralino: "+str(m_n)+", Events: "+str(events)
						#~ print "msbottom: "+str(m_b)+", mneutralino: "+str(m_n)+", ISRNormalization: "+str(events/eventsISRWeighted)
						histo2D.SetBinContent(histo2D.GetXaxis().FindBin(m_b),histo2D.GetYaxis().FindBin(m_n),events)
						histo2DHighPU.SetBinContent(histo2DHighPU.GetXaxis().FindBin(m_b),histo2DHighPU.GetYaxis().FindBin(m_n),eventsISRHighPU)
						histo2DLowPU.SetBinContent(histo2DHighPU.GetXaxis().FindBin(m_b),histo2DLowPU.GetYaxis().FindBin(m_n),eventsISRLowPU)
						
						if events > 0:
							#~ histoISRNormalization.SetBinContent(histo2D.GetXaxis().FindBin(m_b),histo2D.GetYaxis().FindBin(m_n),events/eventsISRWeighted)
							#~ histoISRNormalizationUp.SetBinContent(histo2D.GetXaxis().FindBin(m_b),histo2D.GetYaxis().FindBin(m_n),events/eventsISRWeightedUp)
							#~ histoISRNormalizationDown.SetBinContent(histo2D.GetXaxis().FindBin(m_b),histo2D.GetYaxis().FindBin(m_n),events/eventsISRWeightedDown)
							histoISRNormalization.SetBinContent(histo2D.GetXaxis().FindBin(m_b),histo2D.GetYaxis().FindBin(m_n),eventsISRWeighted/events)
							histoISRNormalizationUp.SetBinContent(histo2D.GetXaxis().FindBin(m_b),histo2D.GetYaxis().FindBin(m_n),eventsISRWeightedUp/events)
							histoISRNormalizationDown.SetBinContent(histo2D.GetXaxis().FindBin(m_b),histo2D.GetYaxis().FindBin(m_n),eventsISRWeightedDown/events)
							
							for weightIndex in range (1,9):
								cutsScaleWeight= "(%s)*scaleWeight%s"%(cuts,str(weightIndex))
								histoScaleWeight = createHistoFromTree(tree, "mSbottom", cutsScaleWeight, 50, 0, 1700)
								eventsScaleWeight = histoScaleWeight.Integral()
								
								histosScaleWeights["ScaleWeight"+str(weightIndex)].SetBinContent(histo2D.GetXaxis().FindBin(m_b),histo2D.GetYaxis().FindBin(m_n),eventsScaleWeight/events)
							
						else:
							histoISRNormalization.SetBinContent(histo2D.GetXaxis().FindBin(m_b),histo2D.GetYaxis().FindBin(m_n),0.)
							histoISRNormalizationUp.SetBinContent(histo2D.GetXaxis().FindBin(m_b),histo2D.GetYaxis().FindBin(m_n),0.)
							histoISRNormalizationDown.SetBinContent(histo2D.GetXaxis().FindBin(m_b),histo2D.GetYaxis().FindBin(m_n),0.)
							
							for weightIndex in range (1,9):
								
								histosScaleWeights["ScaleWeight"+str(weightIndex)].SetBinContent(histo2D.GetXaxis().FindBin(m_b),histo2D.GetYaxis().FindBin(m_n),0.)
						
						
					m_n += stepsize
				m_b += stepsize
	
	m_b_min = 775
	m_b_max = 1600
	m_n_min = 150
	for sample, tree in trees.iteritems():
		if sample == sampleName2:
			m_b = m_b_min
			
			while m_b <= m_b_max:
				print m_b
				if m_b < 800:
					stepsize = 25
				else:
					stepsize = 50
			
				m_n = m_n_min
				
				while m_n < m_b:
					
					if m_b > 1500 or (m_b == 775 and m_n == 750) or (m_b == 800 and m_n == 150) or (m_b == 950 and m_n == 250) or (m_b == 950 and m_n == 300) or (m_b == 950 and m_n == 500) or (m_b == 950 and m_n == 550) or (m_b == 950 and m_n == 850) or (m_b == 950 and m_n == 900):	
					
						cuts = "mSbottom == %s && mNeutralino2 > %s - 5 && mNeutralino2 < %s + 5"%(str(m_b),str(m_n),str(m_n))
						cutsISR = "(mSbottom == %s && mNeutralino2 > %s - 5 && mNeutralino2 < %s + 5)*ISRCorrection"%(str(m_b),str(m_n),str(m_n))
						cutsISRUp = "(mSbottom == %s && mNeutralino2 > %s - 5 && mNeutralino2 < %s + 5)*(ISRCorrection+ISRUncertainty)"%(str(m_b),str(m_n),str(m_n))
						cutsISRDown = "(mSbottom == %s && mNeutralino2 > %s - 5 && mNeutralino2 < %s + 5)*(ISRCorrection-ISRUncertainty)"%(str(m_b),str(m_n),str(m_n))
						
						cutsHighPU = "mSbottom == %s && mNeutralino2 > %s - 5 && mNeutralino2 < %s + 5 && nVertices > 16"%(str(m_b),str(m_n),str(m_n))
						cutsLowPU = "mSbottom == %s && mNeutralino2 > %s - 5 && mNeutralino2 < %s + 5 && nVertices <= 16"%(str(m_b),str(m_n),str(m_n))
						
						
						
						histo = createHistoFromTree(tree, "mSbottom", cuts, 50, 0, 1700)
						histoISR = createHistoFromTree(tree, "mSbottom", cutsISR, 50, 0, 1700)
						histoISRUp = createHistoFromTree(tree, "mSbottom", cutsISRUp, 50, 0, 1700)
						histoISRDown = createHistoFromTree(tree, "mSbottom", cutsISRDown, 50, 0, 1700)
						histoHighPU = createHistoFromTree(tree, "mSbottom", cutsHighPU, 50, 0, 1700)
						histoLowPU = createHistoFromTree(tree, "mSbottom", cutsLowPU, 50, 0, 1700)
						
												
						events = histo.Integral()
						eventsISRWeighted = histoISR.Integral()
						eventsISRWeightedUp = histoISRUp.Integral()
						eventsISRWeightedDown = histoISRDown.Integral()
						eventsISRHighPU = histoHighPU.Integral()
						eventsISRLowPU = histoLowPU.Integral()
							
						print "msbottom: "+str(m_b)+", mneutralino: "+str(m_n)+", Events: "+str(events)
						#~ print "msbottom: "+str(m_b)+", mneutralino: "+str(m_n)+", ISRNormalization: "+str(events/eventsISRWeighted)
						histo2D.SetBinContent(histo2D.GetXaxis().FindBin(m_b),histo2D.GetYaxis().FindBin(m_n),events)
						histo2DHighPU.SetBinContent(histo2DHighPU.GetXaxis().FindBin(m_b),histo2DHighPU.GetYaxis().FindBin(m_n),eventsISRHighPU)
						histo2DLowPU.SetBinContent(histo2DHighPU.GetXaxis().FindBin(m_b),histo2DLowPU.GetYaxis().FindBin(m_n),eventsISRLowPU)
						
						if events > 0:
							#~ histoISRNormalization.SetBinContent(histo2D.GetXaxis().FindBin(m_b),histo2D.GetYaxis().FindBin(m_n),events/eventsISRWeighted)
							#~ histoISRNormalizationUp.SetBinContent(histo2D.GetXaxis().FindBin(m_b),histo2D.GetYaxis().FindBin(m_n),events/eventsISRWeightedUp)
							#~ histoISRNormalizationDown.SetBinContent(histo2D.GetXaxis().FindBin(m_b),histo2D.GetYaxis().FindBin(m_n),events/eventsISRWeightedDown)
							histoISRNormalization.SetBinContent(histo2D.GetXaxis().FindBin(m_b),histo2D.GetYaxis().FindBin(m_n),eventsISRWeighted/events)
							histoISRNormalizationUp.SetBinContent(histo2D.GetXaxis().FindBin(m_b),histo2D.GetYaxis().FindBin(m_n),eventsISRWeightedUp/events)
							histoISRNormalizationDown.SetBinContent(histo2D.GetXaxis().FindBin(m_b),histo2D.GetYaxis().FindBin(m_n),eventsISRWeightedDown/events)
							
							for weightIndex in range (1,9):
								cutsScaleWeight= "(%s)*scaleWeight%s"%(cuts,str(weightIndex))
								histoScaleWeight = createHistoFromTree(tree, "mSbottom", cutsScaleWeight, 50, 0, 1700)
								eventsScaleWeight = histoScaleWeight.Integral()
								
								histosScaleWeights["ScaleWeight"+str(weightIndex)].SetBinContent(histo2D.GetXaxis().FindBin(m_b),histo2D.GetYaxis().FindBin(m_n),eventsScaleWeight/events)
							
						else:
							histoISRNormalization.SetBinContent(histo2D.GetXaxis().FindBin(m_b),histo2D.GetYaxis().FindBin(m_n),0.)
							histoISRNormalizationUp.SetBinContent(histo2D.GetXaxis().FindBin(m_b),histo2D.GetYaxis().FindBin(m_n),0.)
							histoISRNormalizationDown.SetBinContent(histo2D.GetXaxis().FindBin(m_b),histo2D.GetYaxis().FindBin(m_n),0.)
							
							for weightIndex in range (1,9):
								
								histosScaleWeights["ScaleWeight"+str(weightIndex)].SetBinContent(histo2D.GetXaxis().FindBin(m_b),histo2D.GetYaxis().FindBin(m_n),0.)
						
						
					m_n += stepsize
				m_b += stepsize
	
	histo2D.Write()
	histo2DHighPU.Write()
	histo2DLowPU.Write()
	histoISRNormalization.Write()
	histoISRNormalizationUp.Write()
	histoISRNormalizationDown.Write()
	for weightIndex in range (1,9):
		histosScaleWeights["ScaleWeight"+str(weightIndex)].Write()
	f1.Close()
	
