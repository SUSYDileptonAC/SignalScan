#! /usr/bin/env python

import ROOT
import os
import pickle

import sys
sys.path.append('cfg/')
from frameworkStructure import pathes
sys.path.append(pathes.basePath)
from messageLogger import messageLogger as log
from ROOT import TCanvas, TPad, TH1F, TH2F, TH1I, THStack, TLegend, TMath, TF1
import pickle
from defs import sbottom_masses
from math import sqrt

from corrections import rSFOF, rOutIn
from centralConfig import zPredictions, regionsToUse, runRanges


from optparse import OptionParser

def readPickle(name,regionName,runName,MC=False):
	
	if MC:
		if os.path.isfile("../frameWorkBase/shelves/%s_%s_%s_MC.pkl"%(name,regionName,runName)):
			result = pickle.load(open("shelves/%s_%s_%s_MC.pkl"%(name,regionName,runName),"rb"))
		else:
			print "../frameWorkBase/shelves/%s_%s_%s.pkl not found, exiting"%(name,regionName,runName) 		
			sys.exit()		
	else:
		if os.path.isfile("../frameWorkBase/shelves/%s_%s_%s.pkl"%(name,regionName,runName)):
			result = pickle.load(open("../frameWorkBase/shelves/%s_%s_%s.pkl"%(name,regionName,runName),"rb"))
		else:
			print "../frameWorkBase/shelves/%s_%s_%s.pkl not found, exiting"%(name,regionName,runName) 		
			sys.exit()

	return result


def loadPickles(path):
	from glob import glob
	result = {}
	for pklPath in glob(path):
		pklFile = open(pklPath, "r")
		result.update(pickle.load(pklFile))
	return result
	

def getResults(shelve,region,selection):
	
	result = {}
	
	result["rSFOF"] = getattr(rSFOF,region).val
	result["rSFOFErr"] = getattr(rSFOF,region).err
	
	result["lowMassSF"] = shelve[region][selection]["edgeMass"]["EE"] + shelve[region][selection]["edgeMass"]["MM"]
	result["lowMassOF"] = shelve[region][selection]["edgeMass"]["EM"]
	result["lowMassPredSF"] = result["lowMassOF"]*getattr(rSFOF,region).val
	result["lowMassPredStatErrSF"] = result["lowMassOF"]**0.5*getattr(rSFOF,region).val
	result["lowMassPredSystErrSF"] = result["lowMassOF"]*getattr(rSFOF,region).err
	
	result["lowMassZPredSF"] = getattr(getattr(zPredictions,selection).SF,region).val*getattr(rOutIn.lowMass,region).val
	result["lowMassZPredErrSF"] = ((getattr(getattr(zPredictions,selection).SF,region).val*getattr(rOutIn.lowMass,region).err)**2 + (getattr(getattr(zPredictions,selection).SF,region).err*getattr(rOutIn.lowMass,region).val)**2 )**0.5
	

	
	result["belowZSF"] = shelve[region][selection]["belowZ"]["EE"] + shelve[region][selection]["belowZ"]["MM"]
	result["belowZOF"] = shelve[region][selection]["belowZ"]["EM"]
	
	result["belowZPredSF"] = result["belowZOF"]*getattr(rSFOF,region).val
	result["belowZPredStatErrSF"] = result["belowZOF"]**0.5*getattr(rSFOF,region).val
	result["belowZPredSystErrSF"] = result["belowZOF"]*getattr(rSFOF,region).err
	
	result["belowZZPredSF"] = getattr(getattr(zPredictions,selection).SF,region).val*getattr(rOutIn.belowZ,region).val
	result["belowZZPredErrSF"] = ((getattr(getattr(zPredictions,selection).SF,region).val*getattr(rOutIn.belowZ,region).err)**2 + (getattr(getattr(zPredictions,selection).SF,region).err*getattr(rOutIn.belowZ,region).val)**2 )**0.5
	
	
	result["aboveZSF"] = shelve[region][selection]["aboveZ"]["EE"] + shelve[region][selection]["aboveZ"]["MM"]
	result["aboveZOF"] = shelve[region][selection]["aboveZ"]["EM"]
	
	result["aboveZPredSF"] = result["aboveZOF"]*getattr(rSFOF,region).val
	result["aboveZPredStatErrSF"] = result["aboveZOF"]**0.5*getattr(rSFOF,region).val
	result["aboveZPredSystErrSF"] = result["aboveZOF"]*getattr(rSFOF,region).err
	
	result["aboveZZPredSF"] = getattr(getattr(zPredictions,selection).SF,region).val*getattr(rOutIn.aboveZ,region).val
	result["aboveZZPredErrSF"] = ((getattr(getattr(zPredictions,selection).SF,region).val*getattr(rOutIn.aboveZ,region).err)**2 + (getattr(getattr(zPredictions,selection).SF,region).err*getattr(rOutIn.aboveZ,region).val)**2 )**0.5
	
	
	result["highMassSF"] = shelve[region][selection]["highMass"]["EE"] + shelve[region][selection]["highMass"]["MM"]
	result["highMassOF"] = shelve[region][selection]["highMass"]["EM"]
	
	result["highMassPredSF"] = result["highMassOF"]*getattr(rSFOF,region).val
	result["highMassPredStatErrSF"] = result["highMassOF"]**0.5*getattr(rSFOF,region).val
	result["highMassPredSystErrSF"] = result["highMassOF"]*getattr(rSFOF,region).err
	
	result["highMassZPredSF"] = getattr(getattr(zPredictions,selection).SF,region).val*getattr(rOutIn.highMass,region).val
	result["highMassZPredErrSF"] = ((getattr(getattr(zPredictions,selection).SF,region).val*getattr(rOutIn.highMass,region).err)**2 + (getattr(getattr(zPredictions,selection).SF,region).err*getattr(rOutIn.highMass,region).val)**2 )**0.5
	
	
	result["onZSF"] = shelve[region][selection]["zMass"]["EE"] + shelve[region][selection]["zMass"]["MM"]
	result["onZOF"] = shelve[region][selection]["zMass"]["EM"]	
	
	result["onZPredSF"] = result["onZOF"]*getattr(rSFOF,region).val
	result["onZPredStatErrSF"] = result["onZOF"]**0.5*getattr(rSFOF,region).val
	result["onZPredSystErrSF"] = result["onZOF"]*getattr(rSFOF,region).err
	
	result["onZZPredSF"] = getattr(getattr(zPredictions,selection).SF,region).val
	result["onZZPredErrSF"] = getattr(getattr(zPredictions,selection).SF,region).err
	
	return result

	
			   
def writeDataCards():
	
	lumi = 2260.
	printLumi = "2.26"
	
	systUncertainty = 0.15
	
	path = "shelvesYields"
	
	generalSignalLabel = "T6bbllslepton"
	
	CentrNoBTag = {}
	CentrGeOneBTag = {}
	ForwNoBTag = {}
	ForwGeOneBTag = {}
	
	rSFOFCentral = rSFOF.central.val
	rSFOFForward = rSFOF.forward.val
	rSFOFCentralUnc = rSFOF.central.err
	rSFOFForwardUnc = rSFOF.forward.err
	
	countingShelves = {"central": readPickle("cutAndCount",regionsToUse.signal.central.name,runRanges.name), "forward":readPickle("cutAndCount",regionsToUse.signal.forward.name,runRanges.name)}	
		
	resultsCentralNoBTags = getResults(countingShelves,"central","noBTags")
	resultsForwardNoBTags = getResults(countingShelves,"forward","noBTags")
	resultsCentralGeOneBTags = getResults(countingShelves,"central","geOneBTags")
	resultsForwardGeOneBTags = getResults(countingShelves,"forward","geOneBTags")	

	
	etaRegions = ["central","forward"]
	massRegions = ["lowMass","belowZ","onZ","aboveZ","highMass"]
	bTagRegions = ["noBTag","geOneBTag"]
	
	regions = []
	
	for etaRegion in etaRegions:
		for massRegion in massRegions:
			for bTagRegion in bTagRegions:
				regions.append("%s_%s_%s"%(etaRegion,massRegion,bTagRegion))
	

	m_n_min = 150
	m_b_min = 400
	m_b_max = 950
				
	
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
			
			if not ((m_b == 775 and m_n == 750) or (m_b == 800 and m_n == 150)):
		
			

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
				BTag= {}			
				
				
				for region in regions:
	
					Pickles["%s_%s"%(m_sbottom,m_neutralino_2)] = loadPickles("%s/%s_msbottom_%s_mneutralino_%s_%s.pkl"%(path,generalSignalLabel,m_sbottom,m_neutralino_2,region))
					
					Yields["EE_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["EE"]["EEval"]
					Yields["EMu_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuval"]
					Yields["MuMu_%s"%region] =  Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuval"]
					Yields["SFOF_%s"%region] = Yields["EE_%s"%region] + Yields["MuMu_%s"%region] - Yields["EMu_%s"%region]
					Yields["SFOF_%s"%region] = max(Yields["SFOF_%s"%region],0)
					
					MCEvents["EE_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["EE"]["EEMCEvents"]
					MCEvents["EMu_%s"%region] = Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuMCEvents"]
					MCEvents["MuMu_%s"%region] =  Pickles["%s_%s"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuMCEvents"]
					MCEvents["SFOF_%s"%region] = MCEvents["EE_%s"%region] + MCEvents["MuMu_%s"%region] - MCEvents["EMu_%s"%region]
					MCEvents["SFOF_%s"%region] = max(MCEvents["SFOF_%s"%region],0)
					
					if MCEvents["SFOF_%s"%region] > 0:
						statUncertainties["SFOF_%s"%region] = sqrt(MCEvents["EE_%s"%region]+MCEvents["EMu_%s"%region]+MCEvents["MuMu_%s"%region])/MCEvents["SFOF_%s"%region]
					else:
						statUncertainties["SFOF_%s"%region] = 0						
					
				
				ValueNoBTagCentral = {}
				ValueGeOneBTagCentral = {}
				ValueNoBTagForward = {}
				ValueGeOneBTagForward = {}
				
				StatUncertaintyNoBTagCentral = {}
				StatUncertaintyGeOneBTagCentral = {}
				StatUncertaintyNoBTagForward = {}
				StatUncertaintyGeOneBTagForward = {}
				
						
				
				
				
				for massRegion in massRegions:
					ValueNoBTagCentral[massRegion] = Yields["SFOF_central_%s_noBTag"%massRegion]
					ValueGeOneBTagCentral[massRegion] = Yields["SFOF_central_%s_geOneBTag"%massRegion]
					ValueNoBTagForward[massRegion] = Yields["SFOF_forward_%s_noBTag"%massRegion]
					ValueGeOneBTagForward[massRegion] = Yields["SFOF_forward_%s_geOneBTag"%massRegion]
					
					
					StatUncertaintyNoBTagCentral[massRegion]  = statUncertainties["SFOF_central_%s_noBTag"%massRegion]
					StatUncertaintyGeOneBTagCentral[massRegion] = statUncertainties["SFOF_central_%s_geOneBTag"%massRegion]
					StatUncertaintyNoBTagForward[massRegion] = statUncertainties["SFOF_forward_%s_noBTag"%massRegion]
					StatUncertaintyGeOneBTagForward[massRegion] = statUncertainties["SFOF_forward_%s_geOneBTag"%massRegion]
					
					
					n_bins = 4
					n_processes = 2
					n_nuicance_parameters = 11
					
					Name= "T6bbllslepton_%s_%s_%s"%(m_sbottom,m_neutralino_2,massRegion)
					DataCard = open("DataCardsFlatSystematics/%s.txt"%Name,'w')
					
					DataCard.write("# sbottom = %s \n"%m_sbottom)
					DataCard.write("# neutralino 2 = %s \n"%m_neutralino_2)
					DataCard.write("# mass region = %s \n"%massRegion)
					DataCard.write("# Xsection = %s \n"%xsection)
					DataCard.write("imax %s number of bins \n"%n_bins)
					DataCard.write("jmax %s number of processes minus 1 \n"%n_processes)
					DataCard.write("kmax %s number of nuisance parameters \n"%n_nuicance_parameters)
					DataCard.write("------------------------------------------------------------------------------------------- \n")
					DataCard.write("bin          CentrNoBTag   CentrGeOneBTag  ForwNoBTag   ForwGeOneBTag  \n")
					DataCard.write("observation  %s            %s              %s           %s  \n"%(resultsCentralNoBTags["%sSF"%massRegion],resultsCentralGeOneBTags["%sSF"%massRegion],resultsForwardNoBTags["%sSF"%massRegion],resultsForwardGeOneBTags["%sSF"%massRegion]))
					DataCard.write("-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n")
					DataCard.write("bin                                     CentrNoBTag   CentrNoBTag   CentrNoBTag  CentrGeOneBTag   CentrGeOneBTag   CentrGeOneBTag  ForwNoBTag   ForwNoBTag   ForwNoBTag  ForwGeOneBTag   ForwGeOneBTag   ForwGeOneBTag   \n")
					DataCard.write("process                                 SUSY          ZJets         OF           SUSY             ZJets            OF              SUSY         ZJets        OF          SUSY            ZJets           OF     \n")
					DataCard.write("process                                 0             1             2            0                1                2               0            1            2           0               1               2      \n")
					DataCard.write("rate                                    %s            %s            %s           %s               %s               %s              %s           %s           %s          %s             %s              %s      \n"%(ValueNoBTagCentral[massRegion],resultsCentralNoBTags["%sZPredSF"%massRegion],resultsCentralNoBTags["%sOF"%massRegion]*rSFOFCentral,ValueGeOneBTagCentral[massRegion],resultsCentralGeOneBTags["%sZPredSF"%massRegion],resultsCentralGeOneBTags["%sOF"%massRegion]*rSFOFCentral,ValueNoBTagForward[massRegion],resultsForwardNoBTags["%sZPredSF"%massRegion],resultsForwardNoBTags["%sOF"%massRegion]*rSFOFForward,ValueGeOneBTagForward[massRegion],resultsForwardGeOneBTags["%sZPredSF"%massRegion],resultsForwardGeOneBTags["%sOF"%massRegion]*rSFOFForward))
					DataCard.write("-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n")
					
					DataCard.write("SigSystUncertainty              lnN     %s            -             -            %s               -                -               %s           -            -           %s              -               -  \n"%(str(1+systUncertainty),str(1+systUncertainty ),str(1+systUncertainty ),str(1+systUncertainty )))
						
					DataCard.write("SigStatUncertCentrNoBTag_%s     lnN     %s            -             -            -                -                -               -            -            -           -               -               -  \n"%(massRegion,str(1+StatUncertaintyNoBTagCentral[massRegion])))
					DataCard.write("SigStatUncertCentrGeOneBTag_%s  lnN     -             -             -            %s               -                -               -            -            -           -               -               -  \n"%(massRegion,str(1+StatUncertaintyGeOneBTagCentral[massRegion])))
					DataCard.write("SigStatUncertForwNoBTag_%s      lnN     -             -             -            -                -                -               %s           -            -           -               -               -  \n"%(massRegion,str(1+StatUncertaintyNoBTagForward[massRegion])))
					DataCard.write("SigStatUncertForwGeOneBTag_%s   lnN     -             -             -            -                -                -               -            -            -           %s              -               -  \n"%(massRegion,str(1+StatUncertaintyGeOneBTagForward[massRegion])))
											
					DataCard.write("OFStatUncertCentrNoBTag_%s      gmN %s  -             -             %s           -                -                -               -            -            -           -               -               -  \n"%(massRegion,str(resultsCentralNoBTags["%sOF"%massRegion]),str(rSFOFCentral)))
					DataCard.write("OFStatUncertCentrGeOneBTag_%s   gmN %s  -             -             -            -                -               %s               -            -            -           -               -               -  \n"%(massRegion,str(resultsCentralGeOneBTags["%sOF"%massRegion]),str(rSFOFCentral)))
					DataCard.write("OFStatUncertForwNoBTag_%s       gmN %s  -             -             -            -                -                -               -            -            %s          -               -               -  \n"%(massRegion,str(resultsForwardNoBTags["%sOF"%massRegion]),str(rSFOFForward)))
					DataCard.write("OFStatUncertForwGeOneBTag_%s    gmN %s  -             -             -            -                -                -               -            -            -           -               -              %s  \n"%(massRegion,str(resultsForwardGeOneBTags["%sOF"%massRegion]),str(rSFOFForward)))
					
					DataCard.write("RSFOFUncert                     lnN     -             -            %s            -                -               %s               -            -           %s           -               -              %s  \n"%(str(1+rSFOFCentralUnc),str(1+rSFOFCentralUnc),str(1+rSFOFForwardUnc),str(1+rSFOFForwardUnc)))
					
					DataCard.write("ZJetsUncert                     lnN     -             %s            -            -                %s               -               -            %s           -           -               %s              -  \n"%(str(1+resultsCentralNoBTags["%sZPredErrSF"%massRegion]),str(1+resultsCentralGeOneBTags["%sZPredErrSF"%massRegion]),str(1+resultsForwardNoBTags["%sZPredErrSF"%massRegion]),str(1+resultsForwardGeOneBTags["%sZPredErrSF"%massRegion])))
					
							
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
    parser = OptionParser()
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=False,
                                  help="talk about everything")
   
    writeDataCards()
 
