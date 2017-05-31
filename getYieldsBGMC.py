#! /usr/bin/env python

import sys
sys.path.append("/.automount/home/home__home4/institut_1b/schomakers/FrameWork/frameWorkBase")

from helpers import readTrees, totalNumberOfGeneratedEvents,createHistoFromTree, Process, TheStack
from locations import locations

from math import sqrt
from array import array

import ROOT

from sys import argv
import pickle	
from ROOT import TCanvas, TPad, TH1F, TH1I, THStack, TLegend, TF1, TH2F, TH2D, TFile, TMath
import ratios
from defs import sbottom_masses,Backgrounds,Signals,Region,Regions,Plot,getRunRange,getRegion, getPlot
from corrections import triggerEffs, rSFOF, rSFOFDirect

from math import sqrt
from setTDRStyle import setTDRStyle

from ConfigParser import ConfigParser
config = ConfigParser()
config.read("/.automount/home/home__home4/institut_1b/schomakers/FrameWork/SubmitScripts/Input/Master80X_MC.ini")


ROOT.gStyle.SetOptStat(0)

massCuts = {
			"mass20To50":"mll < 50 && mll > 20",
			"mass20To60":"mll < 60 && mll > 20",
			"mass20To70":"mll < 70 && mll > 20",
			"mass20To81":"mll < 81 && mll > 20",
			"mass20To86":"mll < 86 && mll > 20",
			"mass50To81":"mll < 81 && mll > 50",
			"mass50To86":"mll < 86 && mll > 50",
			"mass60To81":"mll < 81 && mll > 60",
			"mass60To86":"mll < 86 && mll > 60",
			"mass70To81":"mll < 81 && mll > 70",
			"mass70To86":"mll < 86 && mll > 70",
			"mass81To101":"mll < 101 && mll > 81",
			"mass86To96":"mll < 96 && mll > 86",
			"mass96To150":"mll < 150 && mll > 96",
			"mass96To200":"mll < 200 && mll > 96",
			"mass101To150":"mll < 150 && mll > 101",
			"mass101To200":"mll < 200 && mll > 101",
			"mass150To200":"mll < 200 && mll > 150",
			"mass200To300":"mll > 200 && mll < 300",
			"mass200To400":"mll > 200 && mll < 400",
			"mass300To400":"mll > 300 && mll < 400",
			"massAbove300":"mll > 300 ",
			"mass400":"mll > 400 ",
			"edgeMass":"mll < 70  && mll > 20 ",
			"lowMass":"mll > 20 && mll < 86 ",
			"highMass":"mll > 96 ",
			"highMassOld":"mll > 101 ",
			}
#~ 
#~ massCuts = {
			#~ "LowMass":"mll < 81 && mll > 20",
			#~ "ZMass":"mll > 81 && mll < 101",
			#~ "HighMass":"mll > 101 ",
			#~ }			
			
nLLCuts = {
			"default":"nLL > 0",
			"lowNLL":"nLL < 21",
			"highNLL":"nLL > 21",
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
			
			   
def plot():
	import ratios
	
	import pickle
	
	signalDenominatorFile = TFile("/.automount/home/home__home4/institut_1b/schomakers/FrameWork/SignalScan/T6bbllsleptonDenominatorHisto3.root")
	denominatorHisto = TH2F(signalDenominatorFile.Get("massScan"))
	

	path = locations.dataSetPathNLL
	#~ path = locations.dataSetPath
	
	
	runRange = getRunRange("Run2016_36fb")
	#~ runRange = getRunRange("Run2016_17fb_Unblinded")
	
	#~ totalEvents = totalNumberOfGeneratedEvents(path)

	
	EETrees = readTrees(path, "EE")
	EMTrees = readTrees(path, "EMu")
	MMTrees = readTrees(path, "MuMu")
	
		
	massRegions = ["mass20To60","mass60To86","mass60To81","mass81To101","mass86To96","mass96To150","mass101To150","mass150To200","mass200To300","mass300To400","mass400","lowMass","highMass","highMassOld","edgeMass"]
	#~ massRegions = ["HighMass"]
	#~ massRegions = ["86To96"]
	#~ nLLRegions = ["default","lowNll","highNll"]
	#~ nLLRegions = ["highNll"]
	#~ massRegions = ["edgeMass"]
	nLLRegions = ["lowNLL","highNLL"]
	HTRegions = ["lowHT","mediumHT","highHT"]
	#~ MT2Regions = ["lowMT2","highMT2"]
	MT2Regions = ["highMT2"]
	
	signalBins = []
	signalCuts = {}

	
	#~ for massRegion in massRegions:
		#~ for nLLRegion in nLLRegions:
			#~ for HTRegion in HTRegions:
				#~ signalBins.append("%s_%s_%s"%(massRegion,nLLRegion,HTRegion))
				#~ signalCuts["%s_%s_%s"%(massRegion,nLLRegion,HTRegion)] = "%s && %s && %s"%(massCuts[massRegion],nLLCuts[nLLRegion],HTCuts[HTRegion])
	#~ for massRegion in massRegions:
		#~ for nLLRegion in nLLRegions:
			#~ for MT2Region in MT2Regions:
				#~ signalBins.append("%s_%s_%s"%(massRegion,nLLRegion,MT2Region))
				#~ signalCuts["%s_%s_%s"%(massRegion,nLLRegion,MT2Region)] = "%s && %s && %s"%(massCuts[massRegion],nLLCuts[nLLRegion],MT2Cuts[MT2Region])
	for massRegion in massRegions:
		for nLLRegion in nLLRegions:
			signalBins.append("%s_%s"%(massRegion,nLLRegion))
			signalCuts["%s_%s"%(massRegion,nLLRegion)] = "%s && %s"%(massCuts[massRegion],nLLCuts[nLLRegion])
	#~ for massRegion in massRegions:
		#~ signalBins.append("%s"%massRegion)
		#~ signalCuts["%s"%massRegion] = "%s"%massCuts[massRegion]
	
	#~ plot = getPlot("mllPlot")
	plot = getPlot("mllPlotROutIn")
	#~ selection = getRegion("SignalInclusive")
	#~ selection = getRegion("SignalCentralOld")
	selection = getRegion("SignalHighMT2DeltaPhiJetMet")
	plot.addRegion(selection)
	plot.cuts = plot.cuts % runRange.runCut
	plot.cuts = plot.cuts.replace("p4.M()","mll")	
	plot.variable = plot.variable.replace("p4.M()","mll")	
		
	counts = {}
				
	#~ backgrounds = ["Rare","SingleTop","TT_Powheg","Diboson","DrellYanTauTau","DrellYan"]
	backgrounds = ["TTZNonFS","RareNonFS","ZZNonFS","WZNonFS","DrellYan"]
	#~ backgrounds = ["RareWZOnZ","RareZZOnZ","RareTTZOnZ","RareRestOnZ"]
	#~ backgrounds = ["TT_Powheg"]
	eventCounts = totalNumberOfGeneratedEvents(path)
	
	defaultCut = plot.cuts
	
	for signalBin in signalBins:
		
		plot.cuts = defaultCut.replace("chargeProduct < 0","chargeProduct < 0 && %s"%(signalCuts[signalBin]))	
		plot.cuts = plot.cuts.replace("metFilterSummary > 0 &&","")	
		plot.cuts = plot.cuts.replace("triggerSummary > 0 &&","")	
		plot.cuts = plot.cuts.replace("p4.Pt()","pt")
		#~ plot.cuts = plot.cuts.replace("genWeight*weight*","")	
		#~ plot.cuts = plot.cuts.replace("genWeight*weight*","fullSimScaleFactor1*fullSimScaleFactor2*")	
		#~ plot.cuts = plot.cuts.replace("leptonFullSimScaleFactor1*leptonFullSimScaleFactor2*","")	
		print plot.cuts
		#~ plot.cuts = plot.cuts.replace("met","metJESUp")	
		#~ plot.cuts = plot.cuts.replace("nJets","nShiftedJetsJESUp")	
		#~ plot.cuts = plot.cuts.replace("nLL","nLLJESUp")	
		#~ plot.cuts = plot.cuts.replace("MT2","MT2JESUp")	
				
		#~ plot.cuts = plot.cuts.replace("met","metJESDown")	
		#~ plot.cuts = plot.cuts.replace("nJets","nShiftedJetsJESDown")	
		#~ plot.cuts = plot.cuts.replace("nLL","nLLJESDown")	
		#~ plot.cuts = plot.cuts.replace("MT2","MT2JESDown")	
				
		processes = []
		for background in backgrounds:
			processes.append(Process(getattr(Backgrounds,background),eventCounts))
		
		
		stackEE = TheStack(processes,runRange.lumi,plot,EETrees,"None",1.0,1.0,1.0,doTopReweighting=True)
		stackMM = TheStack(processes,runRange.lumi,plot,MMTrees,"None",1.0,1.0,1.0,doTopReweighting=True)
		
		histoEE = stackEE.theHistogram		
		histoMM = stackMM.theHistogram
		
		histoEEUp = stackEE.theHistogramXsecUp		
		histoMMUp = stackMM.theHistogramXsecUp
		
		histoEEDown = stackEE.theHistogramXsecDown		
		histoMMDown = stackMM.theHistogramXsecDown
		
		eventCountSF = histoEE.GetEntries()+histoMM.GetEntries()
		histoEE.Scale(triggerEffs.inclusive.effEE.val)
		histoMM.Scale(triggerEffs.inclusive.effMM.val)
		
		histoEEUp.Scale(triggerEffs.inclusive.effEE.val)
		histoMMUp.Scale(triggerEffs.inclusive.effMM.val)
		
		histoEEDown.Scale(triggerEffs.inclusive.effEE.val)
		histoMMDown.Scale(triggerEffs.inclusive.effMM.val)
		
		#~ histoEE.Scale(triggerEffs.central.effEE.val)
		#~ histoMM.Scale(triggerEffs.central.effMM.val)

		#~ histoEEUp.Scale(triggerEffs.central.effEE.val)
		#~ histoMMUp.Scale(triggerEffs.central.effMM.val)

		#~ histoEEDown.Scale(triggerEffs.central.effEE.val)
		#~ histoMMDown.Scale(triggerEffs.central.effMM.val)

		stackEM = TheStack(processes,runRange.lumi,plot,EMTrees,"None",1.0,1.0,1.0,doTopReweighting=True)
		histoEM = stackEM.theHistogram		
		histoEMUp = stackEM.theHistogramXsecUp		
		histoEMDown = stackEM.theHistogramXsecDown		
		histoEM.Scale(triggerEffs.inclusive.effEM.val)
		histoEMUp.Scale(triggerEffs.inclusive.effEM.val)
		histoEMDown.Scale(triggerEffs.inclusive.effEM.val)
		#~ histoEM.Scale(triggerEffs.central.effEM.val)
		#~ histoEMUp.Scale(triggerEffs.central.effEM.val)
		#~ histoEMDown.Scale(triggerEffs.central.effEM.val)
		
		eventYieldSF = histoEE.Integral(0,-1)+histoMM.Integral(0,-1)
		eventYieldSFUp = histoEEUp.Integral(0,-1)+histoMMUp.Integral(0,-1)
		eventYieldSFDown = histoEEDown.Integral(0,-1)+histoMMDown.Integral(0,-1)
		eventYieldOF = histoEM.Integral(0,-1)
		eventYieldOFUp = histoEMUp.Integral(0,-1)
		eventYieldOFDown = histoEMDown.Integral(0,-1)
		
		#~ eventCountOF = histoEM.GetEntries()
		
		counts["%s_SF"%signalBin] = eventYieldSF
		counts["%s_SF_Up"%signalBin] = eventYieldSFUp
		counts["%s_SF_Down"%signalBin] = eventYieldSFDown
		counts["MCEvents_%s_SF"%signalBin] = eventCountSF
		counts["%s_OF"%signalBin] = eventYieldOF*rSFOFDirect.inclusive.val
		counts["%s_OF_Up"%signalBin] = eventYieldOFUp*rSFOFDirect.inclusive.val
		counts["%s_OF_Down"%signalBin] = eventYieldOFDown*rSFOFDirect.inclusive.val
		#~ counts["%s_OF"%signalBin] = eventYieldOF*rSFOFDirect.central.val
		#~ counts["%s_OF_Up"%signalBin] = eventYieldOFUp*rSFOFDirect.central.val
		#~ counts["%s_OF_Down"%signalBin] = eventYieldOFDown*rSFOFDirect.central.val

	#~ outFilePkl = open("shelvesMT2/TTbar_40fb.pkl","w")
	#~ outFilePkl = open("shelves/FullMC_40fb_OldOnZ.pkl","w")
	#~ outFilePkl = open("shelvesMT2/FullMC_40fb_MassStudy_DeltaPhi.pkl","w")
	#~ outFilePkl = open("shelvesMT2/FullMC_legacy_17fb.pkl","w")
	outFilePkl = open("shelvesMT2/RareOnZ_Powheg.pkl","w")
	#~ outFilePkl = open("shelvesMT2/RareOnZ_3.pkl","w")
	#~ outFilePkl = open("shelves/OnZBG_ICHEP_36fb.pkl","w")
	#~ outFilePkl = open("shelves/OnZBG_legacy_36fb.pkl","w")
	#~ outFilePkl = open("shelves/TTZ.pkl","w")
	pickle.dump(counts, outFilePkl)
	outFilePkl.close()
	
	
			
			

				
	

# This method just waits for a button to be pressed
def waitForInput():
    raw_input("Press any key to continue!")
    return

# entry point
#-------------
if (__name__ == "__main__"):
    # use option parser to allow verbose mode
    
    # start
    #~ plot(opts.base_path, int(opts.nEvents), opts.observable)
    plot()
 
