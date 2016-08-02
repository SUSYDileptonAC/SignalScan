#! /usr/bin/env python

import ROOT


from optparse import OptionParser


def loadPickles(path):
	from glob import glob
	import pickle
	result = {}
	for pklPath in glob(path):
		pklFile = open(pklPath, "r")
		result.update(pickle.load(pklFile))
	return result
	
			   
def writeDataCards():
	import sys
	sys.path.append('cfg/')
	from frameworkStructure import pathes
	sys.path.append(pathes.basePath)
	from messageLogger import messageLogger as log
	from ROOT import TCanvas, TPad, TH1F, TH2F, TH1I, THStack, TLegend, TMath, TF1
	import pickle
	from defs import sbottom_masses
	from math import sqrt

	
	path = "shelves"
	
	generalSignalLabel = "T6bbllslepton"
	
	LowMassLowNll = {}
	LowMassHighNll = {}
	HighMassLowNll = {}
	HighMassHighNll = {}
	
	rSFOFInclusive = 1.091
	rSFOFInclusiveUnc = 0.023
	
	LowMassLowNll["observation"] = 1417
	LowMassLowNll["TTBar"] = 1260
	LowMassLowNll["TTBarUncertainty"] = 0.028
	LowMassLowNll["DY"] = 8.0
	LowMassLowNll["DYUncertainty"] = 0.4
	
	LowMassHighNll["observation"] = 135
	LowMassHighNll["TTBar"] = 97
	LowMassHighNll["TTBarUncertainty"] = 0.095
	LowMassHighNll["DY"] = 4.3
	LowMassHighNll["DYUncertainty"] = 0.40
	
	HighMassLowNll["observation"] = 2347
	HighMassLowNll["TTBar"] = 2233
	HighMassLowNll["TTBarUncertainty"] = 0.021
	HighMassLowNll["DY"] = 4.5
	HighMassLowNll["DYUncertainty"] = 0.4
	
	HighMassHighNll["observation"] = 285
	HighMassHighNll["TTBar"] = 191
	HighMassHighNll["TTBarUncertainty"] = 0.072
	HighMassHighNll["DY"] = 2.4
	HighMassHighNll["DYUncertainty"] = 0.42
	
	
	
	regions = ["LowMassLowNll","LowMassHighNll","HighMassLowNll","HighMassHighNll"]
	
	

	m_n_min = 150
	m_b_min = 500
	m_b_max = 950
	
	
	TriggerEffUncertainty = 0.05
	PDFUncertainty = 0.
	LumiUncertainty = 0.062
	ElectronUncertainty = 0.05
	MuonUncertainty = 0.01
				
	
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
			
			if not ((m_b == 775 and m_n == 750) or (m_b == 800 and m_n == 150) or (m_b == 950 and m_n == 900) or (m_b == 950 and m_n == 850) or (m_b == 950 and m_n == 550) or (m_b == 950 and m_n == 500) or (m_b == 950 and m_n == 300) or (m_b == 950 and m_n == 250)):
		
			

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
				Met= {}			
				
				
				for region in regions:
					
					Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)] = loadPickles("%s/%s_msbottom_%s_mneutralino_%s_%s_EE.pkl"%(path,generalSignalLabel,m_sbottom,m_neutralino_2,region))
					Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)] = loadPickles("%s/%s_msbottom_%s_mneutralino_%s_%s_EMu.pkl"%(path,generalSignalLabel,m_sbottom,m_neutralino_2,region))
					Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)] = loadPickles("%s/%s_msbottom_%s_mneutralino_%s_%s_MuMu.pkl"%(path,generalSignalLabel,m_sbottom,m_neutralino_2,region))
					
					Yields["EE_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EEval"]
					Yields["EMu_%s"%region] = Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuval"]
					Yields["MuMu_%s"%region] =  Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuval"]
					Yields["SFOF_%s"%region] = Yields["EE_%s"%region] + Yields["MuMu_%s"%region] - Yields["EMu_%s"%region]
					Yields["SFOF_%s"%region] = max(Yields["SFOF_%s"%region],0)
					
					MCEvents["EE_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EEMCEvents"]
					MCEvents["EMu_%s"%region] = Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuMCEvents"]
					MCEvents["MuMu_%s"%region] =  Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuMCEvents"]
					MCEvents["SFOF_%s"%region] = MCEvents["EE_%s"%region] + MCEvents["MuMu_%s"%region] - MCEvents["EMu_%s"%region]
					MCEvents["SFOF_%s"%region] = max(MCEvents["SFOF_%s"%region],0)
					
					if MCEvents["SFOF_%s"%region] > 0:
						statUncertainties["SFOF_%s"%region] = sqrt(MCEvents["EE_%s"%region]+MCEvents["EMu_%s"%region]+MCEvents["MuMu_%s"%region])/MCEvents["SFOF_%s"%region]
					else:
						statUncertainties["SFOF_%s"%region] = 0
					
					
					### JES Uncertainty
					JES["Mean_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EEJESMean"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuJESMean"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuJESMean"] 
					
					JES["JESUp_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EEJESUp"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuJESUp"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuJESUp"] 
					JES["JESDown_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EEJESDown"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuJESDown"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuJESDown"] 
					
					if JES["Mean_%s"%region] > 0:
						systUncertainties["JESUncertainty_%s"%region] = max(abs(JES["JESUp_%s"%region]-JES["Mean_%s"%region])/JES["Mean_%s"%region],abs(JES["JESDown_%s"%region]-JES["Mean_%s"%region])/JES["Mean_%s"%region])
					else:
						systUncertainties["JESUncertainty_%s"%region] = 0
					
					
					
					### Lepton FastSim Uncertainty
					
					LeptonFastSim["Mean_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EELeptonFastSimMean"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuLeptonFastSimMean"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuLeptonFastSimMean"] 
					LeptonFastSim["MeanUnscaled_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EELeptonNoFastSimScaleFactor"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuLeptonNoFastSimScaleFactor"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuLeptonNoFastSimScaleFactor"] 
					
					if LeptonFastSim["Mean_%s"%region] > 0:
						systUncertainties["LeptonFastSimUncertainty_%s"%region] = abs(LeptonFastSim["MeanUnscaled_%s"%region]-LeptonFastSim["Mean_%s"%region])/LeptonFastSim["Mean_%s"%region]
					else:
						systUncertainties["LeptonFastSimUncertainty_%s"%region] = 0
						
					### Lepton FullSim Uncertainty
					
					LeptonFullSim["Mean_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EELeptonFullSimMean"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuLeptonFullSimMean"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuLeptonFullSimMean"] 
					LeptonFullSim["MeanUnscaled_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EELeptonNoFullSimScaleFactor"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuLeptonNoFullSimScaleFactor"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuLeptonNoFullSimScaleFactor"] 
					
					if LeptonFullSim["Mean_%s"%region] > 0:
						systUncertainties["LeptonFullSimUncertainty_%s"%region] = abs(LeptonFullSim["MeanUnscaled_%s"%region]-LeptonFullSim["Mean_%s"%region])/LeptonFullSim["Mean_%s"%region]
					else:
						systUncertainties["LeptonFullSimUncertainty_%s"%region] = 0
					
					systUncertainties["LeptonFullSimUncertainty_%s"%region] = 0.07
						
					
					###  Pileup Uncertainty
					Pileup["Mean_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EEPileupMean"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuPileupMean"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuPileupMean"] 
					
					Pileup["PileupUp_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EEPileupUp"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuPileupUp"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuPileupUp"] 
					Pileup["PileupDown_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EEPileupDown"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuPileupDown"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuPileupDown"] 
					
					if Pileup["Mean_%s"%region] > 0:
						systUncertainties["pileupUncertainty_%s"%region] = max(abs(Pileup["PileupUp_%s"%region]-Pileup["Mean_%s"%region])/Pileup["Mean_%s"%region],abs(Pileup["PileupDown_%s"%region] -Pileup["Mean_%s"%region])/Pileup["Mean_%s"%region] )
					else:
						systUncertainties["pileupUncertainty_%s"%region] = 0
					
					### ISR Uncertainty
					ISR["Mean_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EEISRMean"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuISRMean"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuISRMean"] 
					
					ISR["ISRUp_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EEISRUp"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuISRUp"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuISRUp"]  
					ISR["ISRDown_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EEISRDown"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuISRDown"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuISRDown"]  
					
					if ISR["Mean_%s"%region] > 0:
						systUncertainties["ISRUncertainty_%s"%region] = max(abs(ISR["ISRUp_%s"%region]-ISR["Mean_%s"%region])/ISR["Mean_%s"%region],abs(ISR["ISRDown_%s"%region]-ISR["Mean_%s"%region])/ISR["Mean_%s"%region])
					else:
						systUncertainties["ISRUncertainty_%s"%region] = 0
						
					### BTag Uncertainty
					BTag["Mean_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EEISRMean"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuISRMean"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuISRMean"] 
					
					BTag["BTagHeavy_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EEbTagHeavy"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMubTagHeavy"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMubTagHeavy"]  
					BTag["BTagLight_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EEbTagLight"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMubTagLight"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMubTagLight"]  
					
					if BTag["Mean_%s"%region] > 0:
						systUncertainties["BTagHeavyUncertainty_%s"%region] = abs(BTag["BTagHeavy_%s"%region]-BTag["Mean_%s"%region])/BTag["Mean_%s"%region]
						systUncertainties["BTagLightUncertainty_%s"%region] = abs(BTag["BTagLight_%s"%region]-BTag["Mean_%s"%region])/BTag["Mean_%s"%region]	
					else:
						systUncertainties["BTagHeavyUncertainty_%s"%region] = 0
						systUncertainties["BTagLightUncertainty_%s"%region] = 0
						
					### FastSim Met Uncertainty					
					Met["Met_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EEMet"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuMet"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuMet"]  
					Met["GenMet_%s"%region] = Pickles["%s_%s_EE"%(m_sbottom,m_neutralino_2)]["EE"]["EEGenMet"] + Pickles["%s_%s_MuMu"%(m_sbottom,m_neutralino_2)]["MuMu"]["MuMuGenMet"] - Pickles["%s_%s_EMu"%(m_sbottom,m_neutralino_2)]["EMu"]["EMuGenMet"]  
					
					if Yields["SFOF_%s"%region] > 0:
						systUncertainties["MetUncertainty_%s"%region] = 0.5*abs(Met["Met_%s"%region]-Met["GenMet_%s"%region])/Yields["SFOF_%s"%region]
					else:
						systUncertainties["MetUncertainty_%s"%region] = 0
						
					### Trigger eff uncertainty
					if Yields["SFOF_%s"%region] > 0:
						systUncertainties["TriggerEffUncertainty_%s"%region] = sqrt((Yields["EE_%s"%region]*TriggerEffUncertainty)**2 + (Yields["MuMu_%s"%region]*TriggerEffUncertainty)**2 + (Yields["EMu_%s"%region]*TriggerEffUncertainty)**2)/Yields["SFOF_%s"%region]
					else:
						systUncertainties["TriggerEffUncertainty_%s"%region] = 0
							
	
	
				
				n_bins = 4
				n_processes = 2
				n_nuicance_parameters = 19
				
				Name= "T6bbllslepton_%s_%s"%(m_sbottom,m_neutralino_2)
				DataCard = open("DataCards/%s.txt"%Name,'w')
				
				DataCard.write("# sbottom = %s \n"%m_sbottom)
				DataCard.write("# neutralino 2 = %s \n"%m_neutralino_2)
				DataCard.write("# Xsection = %s \n"%xsection)
				DataCard.write("imax %s number of bins \n"%n_bins)
				DataCard.write("jmax %s number of processes minus 1 \n"%n_processes)
				DataCard.write("kmax %s number of nuisance parameters \n"%n_nuicance_parameters)
				DataCard.write("------------------------------------------------------------------------------------------- \n")
				DataCard.write("bin          LowMassLowNll   LowMassHighNll   HighMassLowNll    HighMassHighNll \n")
				DataCard.write("observation   %s              %s               %s                %s \n"%(LowMassLowNll["observation"],LowMassHighNll["observation"],HighMassLowNll["observation"],HighMassHighNll["observation"]))
				DataCard.write("-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n")
				DataCard.write("bin                                     LowMassLowNll    LowMassLowNll    LowMassLowNll   LowMassHighNll   LowMassHighNll   LowMassHighNll  HighMassLowNll   HighMassLowNll   HighMassLowNll  HighMassHighNll   HighMassHighNll   HighMassHighNll   \n")
				DataCard.write("process                                 SUSY             ZJets            OF              SUSY             ZJets            OF              SUSY             ZJets            OF              SUSY             ZJets            OF     \n")
				DataCard.write("process                                 0                1                2               0                1                2               0                1                2               0                1                2      \n")
				DataCard.write("rate                                    %s               %s               %s              %s               %s               %s              %s               %s               %s              %s               %s               %s     \n"%(Yields["SFOF_LowMassLowNll"],LowMassLowNll["DY"],LowMassLowNll["TTBar"]*rSFOFInclusive,Yields["SFOF_LowMassHighNll"],LowMassHighNll["DY"],LowMassHighNll["TTBar"]*rSFOFInclusive,Yields["SFOF_HighMassLowNll"],HighMassLowNll["DY"],HighMassLowNll["TTBar"]*rSFOFInclusive,Yields["SFOF_HighMassHighNll"],HighMassHighNll["DY"],HighMassHighNll["TTBar"]*rSFOFInclusive))
				DataCard.write("-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n")
				
				DataCard.write("SigTriggerEffUncertainty        lnN     %s               -                -               %s               -                -               %s               -                -               %s               -                -  \n"%(str(1+TriggerEffUncertainty ),str(1+TriggerEffUncertainty ),str(1+TriggerEffUncertainty ),str(1+TriggerEffUncertainty )))
				DataCard.write("SigLumiUncertainty              lnN     %s               -                -               %s               -                -               %s               -                -               %s               -                -  \n"%(str(1+LumiUncertainty ),str(1+LumiUncertainty ),str(1+LumiUncertainty ),str(1+LumiUncertainty )))
				DataCard.write("SigJESUncertainty               lnN     %s               -                -               %s               -                -               %s               -                -               %s               -                -  \n"%(str(1+systUncertainties["JESUncertainty_LowMassLowNll"]) ,str(1+systUncertainties["JESUncertainty_LowMassHighNll"]) ,str(1+systUncertainties["JESUncertainty_HighMassLowNll"]) ,str(1+systUncertainties["JESUncertainty_HighMassHighNll"] )))
				DataCard.write("SigLeptonFullSimUncertainty     lnN     %s               -                -               %s               -                -               %s               -                -               %s               -                -  \n"%(str(1+systUncertainties["LeptonFullSimUncertainty_LowMassLowNll"]) ,str(1+systUncertainties["LeptonFullSimUncertainty_LowMassHighNll"]) ,str(1+systUncertainties["LeptonFullSimUncertainty_HighMassLowNll"]) ,str(1+systUncertainties["LeptonFullSimUncertainty_HighMassHighNll"] )))
				DataCard.write("SigPileupUncertainty            lnN     %s               -                -               %s               -                -               %s               -                -               %s               -                -  \n"%(str(1+systUncertainties["pileupUncertainty_LowMassLowNll"]) ,str(1+systUncertainties["pileupUncertainty_LowMassHighNll"]) ,str(1+systUncertainties["pileupUncertainty_HighMassLowNll"]) ,str(1+systUncertainties["pileupUncertainty_HighMassHighNll"] )))
				DataCard.write("SigISRUncertainty               lnN     %s               -                -               %s               -                -               %s               -                -               %s               -                -  \n"%(str(1+systUncertainties["ISRUncertainty_LowMassLowNll"]) ,str(1+systUncertainties["ISRUncertainty_LowMassHighNll"]) ,str(1+systUncertainties["ISRUncertainty_HighMassLowNll"]) ,str(1+systUncertainties["ISRUncertainty_HighMassHighNll"] )))
				DataCard.write("SigBTagHeavyUncertainty         lnN     %s               -                -               %s               -                -               %s               -                -               %s               -                -  \n"%(str(1+systUncertainties["BTagHeavyUncertainty_LowMassLowNll"]) ,str(1+systUncertainties["BTagHeavyUncertainty_LowMassHighNll"]) ,str(1+systUncertainties["BTagHeavyUncertainty_HighMassLowNll"]) ,str(1+systUncertainties["BTagHeavyUncertainty_HighMassHighNll"] )))
				DataCard.write("SigBTagLightUncertainty         lnN     %s               -                -               %s               -                -               %s               -                -               %s               -                -  \n"%(str(1+systUncertainties["BTagLightUncertainty_LowMassLowNll"]) ,str(1+systUncertainties["BTagLightUncertainty_LowMassHighNll"]) ,str(1+systUncertainties["BTagLightUncertainty_HighMassLowNll"]) ,str(1+systUncertainties["BTagLightUncertainty_HighMassHighNll"] )))
				DataCard.write("SigMetUncertainty               lnN     %s               -                -               %s               -                -               %s               -                -               %s               -                -  \n"%(str(1+systUncertainties["MetUncertainty_LowMassLowNll"]) ,str(1+systUncertainties["MetUncertainty_LowMassHighNll"]) ,str(1+systUncertainties["MetUncertainty_HighMassLowNll"]) ,str(1+systUncertainties["MetUncertainty_HighMassHighNll"] )))
				
				
				DataCard.write("SigStatUncertLowMassLowNll      lnN     %s               -                -               -                -                -               -                -                -               -                -                -  \n"%(str(1+statUncertainties["SFOF_LowMassLowNll"])))
				DataCard.write("SigStatUncertLowMassHighNll     lnN     -                -                -               %s               -                -               -                -                -               -                -                -  \n"%(str(1+statUncertainties["SFOF_LowMassHighNll"])))
				DataCard.write("SigStatUncertHighMassLowNll     lnN     -                -                -               -                -                -               %s               -                -               -                -                -  \n"%(str(1+statUncertainties["SFOF_HighMassLowNll"])))
				DataCard.write("SigStatUncertHighMassHighNll    lnN     -                -                -               -                -                -               -                -                -               %s               -                -  \n"%(str(1+statUncertainties["SFOF_HighMassHighNll"])))
										
				DataCard.write("OFSystUncertLowMassLowNll       gmN %s  -                -               %s               -                -                -               -                -                -               -                -                -  \n"%(str(LowMassLowNll["TTBar"]),str(rSFOFInclusive)))
				DataCard.write("OFSystUncertLowMassHighNll      gmN %s  -                -                -               -                -                %s              -                -                -               -                -                -  \n"%(str(LowMassHighNll["TTBar"]),str(rSFOFInclusive)))
				DataCard.write("OFSystUncertHighMassLowNll      gmN %s  -                -                -               -                -                -               -                -               %s               -                -                -  \n"%(str(HighMassLowNll["TTBar"]),str(rSFOFInclusive)))
				DataCard.write("OFSystUncertHighMassHighNll     gmN %s  -                -                -               -                -                -               -                -                -               -                -                %s \n"%(str(HighMassHighNll["TTBar"]),str(rSFOFInclusive)))
				
				DataCard.write("RSFOFUncert                     lnN     -                -               %s               -                -                %s              -                -               %s               -                -                %s  \n"%(str(1+rSFOFInclusiveUnc),str(1+rSFOFInclusiveUnc),str(1+rSFOFInclusiveUnc),str(1+rSFOFInclusiveUnc)))
				
				DataCard.write("ZJetsUncert                     lnN     -                %s               -               -                %s               -               -                %s              -                -                %s               - \n"%(str(1+LowMassLowNll["DYUncertainty"]),str(1+LowMassHighNll["DYUncertainty"]),str(1+HighMassLowNll["DYUncertainty"]),str(1+HighMassHighNll["DYUncertainty"])))
					
					
							
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
 
