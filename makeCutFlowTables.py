#!/usr/bin/env python



def saveTable(table, name):
	tabFile = open("tab/table_%s.tex"%name, "w")
	tabFile.write(table)
	tabFile.close()

	#~ print table
	
def loadPickles(path):
	from glob import glob
	import pickle
	result = {}
	for pklPath in glob(path):
		pklFile = open(pklPath, "r")
		result.update(pickle.load(pklFile))
	return result
	


				  


def main():
	from math import sqrt
	
	#~ Slepton: 900-150, 900-500, 1000-300, 1200-200, 1200-1000
	m_sbottom = 1200
	M_SBOTTOM = "m_b_"+str(m_sbottom)
	m_neutralino_2 = 1000
	
	fileName = "T6bbllslepton_msbottom_%s_mneutralino_%s"%(m_sbottom,m_neutralino_2)
	tableName = "CutFlow_T6bbllslepton_msbottom_%s_mneutralino_%s"%(m_sbottom,m_neutralino_2)
	tableName2 = "CutFlow_T6bbllslepton_msbottom_%s_mneutralino_%s_2"%(m_sbottom,m_neutralino_2)
		
	CutString = ["basicCut","mll20Cut","pt25Cut","signalRegionCut","met150Cut","nJets2Cut","deltaPhiCut"]
	
	
	Pickles = {}
	
	#~ lineTemplate = r"%s &  %.1f $\pm$ %.1f  $\pm$ %.1f & %.1f $\pm$ %.1f  $\pm$ %.1f \\"+"\n"
	#~ lineTemplate = r"%s &%.1f$\pm$%.1f$\pm$%.1f&%.1f$\pm$%.1f$\pm$%.1f\\"+"\n"
	#~ lineTemplate2 = r"%s &%.1f$\pm$%.1f&%.1f$\pm$%.1f\\"+"\n"
	lineTemplate = r"%s &%.1f$\pm$%.1f$\pm$%.1f\\"+"\n"
	lineTemplate2 = r"%s &%.1f$\pm$%.1f\\"+"\n"
	#~ lineTemplateAllEvents = r"%s & \multicolumn{2}{c}{%.0f} \\"+"\n"
	lineTemplateAllEvents = r"%s & %.0f \\"+"\n"
	#~ lineTemplate2 = r"%s & %.1f$\pm$%.1f & %.1f$\pm$%.1f & %.1f$\pm$%.1f \\"+"\n"
	horizontalLine = r"\hline"+"\n"
	
	for Cut in CutString:
		if Cut in ["basicCut","mll20Cut","pt25Cut"]:
			MllCuts = ["noMllCut"]
			NllCuts = ["inclusiveNll"]
		else:
			MllCuts = ["noMllCut","20To60","60To86","86To96","96To150","150To200","200To300","300To400","Above400"]
			NllCuts = ["inclusiveNll", "lowNll", "highNll"]
		for MllCut in MllCuts:
			for NllCut in NllCuts:
				Pickles["%s_%s_%s_%s"%(fileName,Cut,MllCut,NllCut)] = loadPickles("shelvesCutFlowTables/%s_%s_%s_%s.pkl"%(fileName,Cut,MllCut,NllCut))
	
	tableTemplate =r"""
\begin{tabular}{l|c}
\hline
slepton model, m$_{\tilde{b}}=1200$, m$_{\tilde{\chi}_{2}^{0}}=1000$ GeV & SF events \\
\hline
%s
\hline
\end{tabular}
"""

	table =""
				
	table += lineTemplateAllEvents%("Total dilepton events", Pickles["%s_basicCut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["totalEventYield"])
	table += lineTemplate%("$\geq 2$ leptons ($\ell^{\pm}\ell^{\mp}$), p$_{T}>25(20)$ GeV", 
							Pickles["%s_basicCut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_basicCut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_basicCut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_basicCut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SystUncertainty"]*Pickles["%s_basicCut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += lineTemplate%("m$_{\ell\ell}>20$ GeV", 
							Pickles["%s_mll20Cut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_mll20Cut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_mll20Cut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_mll20Cut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SystUncertainty"]*Pickles["%s_mll20Cut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							)	
							
	table += lineTemplate%("p$_{T}^{\ell\ell}>25$ GeV", 
							Pickles["%s_pt25Cut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_pt25Cut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_pt25Cut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_pt25Cut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SystUncertainty"]*Pickles["%s_pt25Cut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							)
	
	table += lineTemplate%("jet-multiplicity $\geq 2$", 
							Pickles["%s_nJets2Cut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_nJets2Cut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_nJets2Cut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_nJets2Cut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SystUncertainty"]*Pickles["%s_nJets2Cut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += lineTemplate%("$\Delta\phi$(jet$_{1,2}$,MET)$>0.4$", 
							Pickles["%s_deltaPhiCut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_deltaPhiCut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_deltaPhiCut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_deltaPhiCut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SystUncertainty"]*Pickles["%s_deltaPhiCut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += lineTemplate%("MET $> 150$ GeV", 
							Pickles["%s_met150Cut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_met150Cut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_met150Cut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_met150Cut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SystUncertainty"]*Pickles["%s_met150Cut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							)
	
	table += lineTemplate%("MT2 $> 80$ GeV", 
							Pickles["%s_signalRegionCut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SystUncertainty"]*Pickles["%s_signalRegionCut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += horizontalLine
	table += "\multicolumn{2}{c}{ \textbf{ttbar like}} \\"+"\n"
	table += horizontalLine
	table += lineTemplate%("20 $<$ m$_{\ell\ell} <$ 60 GeV", 
							Pickles["%s_signalRegionCut_20To60_lowNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_20To60_lowNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_20To60_lowNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_20To60_lowNll"%(fileName)]["Signal"]["SystUncertainty"]*Pickles["%s_signalRegionCut_20To60_lowNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += lineTemplate%("60 $<$ m$_{\ell\ell} <$ 86 GeV", 
							Pickles["%s_signalRegionCut_60To86_lowNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_60To86_lowNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_60To86_lowNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_60To86_lowNll"%(fileName)]["Signal"]["SystUncertainty"]*Pickles["%s_signalRegionCut_60To86_lowNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += lineTemplate%("96 $<$ m$_{\ell\ell} <$ 150 GeV", 
							Pickles["%s_signalRegionCut_96To150_lowNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_96To150_lowNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_96To150_lowNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_96To150_lowNll"%(fileName)]["Signal"]["SystUncertainty"]*Pickles["%s_signalRegionCut_96To150_lowNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += lineTemplate%("150 $<$ m$_{\ell\ell} <$ 200 GeV", 
							Pickles["%s_signalRegionCut_150To200_lowNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_150To200_lowNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_150To200_lowNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_150To200_lowNll"%(fileName)]["Signal"]["SystUncertainty"]*Pickles["%s_signalRegionCut_150To200_lowNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += lineTemplate%("200 $<$ m$_{\ell\ell} <$ 300 GeV", 
							Pickles["%s_signalRegionCut_200To300_lowNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_200To300_lowNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_200To300_lowNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_200To300_lowNll"%(fileName)]["Signal"]["SystUncertainty"]*Pickles["%s_signalRegionCut_200To300_lowNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += lineTemplate%("300 $<$ m$_{\ell\ell} <$ 400 GeV", 
							Pickles["%s_signalRegionCut_300To400_lowNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_300To400_lowNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_300To400_lowNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_300To400_lowNll"%(fileName)]["Signal"]["SystUncertainty"]*Pickles["%s_signalRegionCut_300To400_lowNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += lineTemplate%("m$_{\ell\ell} >$ 400 GeV", 
							Pickles["%s_signalRegionCut_Above400_lowNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_Above400_lowNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_Above400_lowNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_Above400_lowNll"%(fileName)]["Signal"]["SystUncertainty"]*Pickles["%s_signalRegionCut_Above400_lowNll"%(fileName)]["Signal"]["SFYield"],
							)
	table += horizontalLine						
	table += "\multicolumn{2}{c}{ \textbf{non ttbar like}} \\"+"\n"
	table += horizontalLine
	table += lineTemplate%("20 $<$ m$_{\ell\ell} <$ 60 GeV", 
							Pickles["%s_signalRegionCut_20To60_highNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_20To60_highNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_20To60_highNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_20To60_highNll"%(fileName)]["Signal"]["SystUncertainty"]*Pickles["%s_signalRegionCut_20To60_highNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += lineTemplate%("60 $<$ m$_{\ell\ell} <$ 86 GeV", 
							Pickles["%s_signalRegionCut_60To86_highNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_60To86_highNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_60To86_highNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_60To86_highNll"%(fileName)]["Signal"]["SystUncertainty"]*Pickles["%s_signalRegionCut_60To86_highNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += lineTemplate%("96 $<$ m$_{\ell\ell} <$ 150 GeV", 
							Pickles["%s_signalRegionCut_96To150_highNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_96To150_highNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_96To150_highNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_96To150_highNll"%(fileName)]["Signal"]["SystUncertainty"]*Pickles["%s_signalRegionCut_96To150_highNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += lineTemplate%("150 $<$ m$_{\ell\ell} <$ 200 GeV", 
							Pickles["%s_signalRegionCut_150To200_highNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_150To200_highNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_150To200_highNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_150To200_highNll"%(fileName)]["Signal"]["SystUncertainty"]*Pickles["%s_signalRegionCut_150To200_highNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += lineTemplate%("200 $<$ m$_{\ell\ell} <$ 300 GeV", 
							Pickles["%s_signalRegionCut_200To300_highNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_200To300_highNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_200To300_highNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_200To300_highNll"%(fileName)]["Signal"]["SystUncertainty"]*Pickles["%s_signalRegionCut_200To300_highNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += lineTemplate%("300 $<$ m$_{\ell\ell} <$ 400 GeV", 
							Pickles["%s_signalRegionCut_300To400_highNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_300To400_highNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_300To400_highNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_300To400_highNll"%(fileName)]["Signal"]["SystUncertainty"]*Pickles["%s_signalRegionCut_300To400_highNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += lineTemplate%("m$_{\ell\ell} >$ 400 GeV", 
							Pickles["%s_signalRegionCut_Above400_highNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_Above400_highNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_Above400_highNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_Above400_highNll"%(fileName)]["Signal"]["SystUncertainty"]*Pickles["%s_signalRegionCut_Above400_highNll"%(fileName)]["Signal"]["SFYield"],
							)								
	
							

	saveTable(tableTemplate%(table), tableName)
	
	
	tableTemplate2 =r"""
\begin{tabular}{l|c}
\hline
slepton model, m$_{\tilde{b}}=1200$, m$_{\tilde{\chi}_{2}^{0}}=1000$ GeV & SF events \\
\hline
%s
\hline
\end{tabular}
"""

	table =""
				
				
	table += lineTemplateAllEvents%("Total dilepton events", Pickles["%s_basicCut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["totalEventYield"])
	table += lineTemplate2%("$\geq 2$ leptons ($\ell^{\pm}\ell^{\mp}$), p$_{T}>25(20)$ GeV", 
							Pickles["%s_basicCut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_basicCut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_basicCut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += lineTemplate2%("m$_{\ell\ell}>20$ GeV", 
							Pickles["%s_mll20Cut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_mll20Cut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_mll20Cut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							)	
							
	table += lineTemplate2%("p$_{T}^{\ell\ell}>25$ GeV", 
							Pickles["%s_pt25Cut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_pt25Cut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_pt25Cut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							)
	
	table += lineTemplate2%("jet-multiplicity $\geq 2$", 
							Pickles["%s_nJets2Cut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_nJets2Cut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_nJets2Cut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += lineTemplate2%("$\Delta\phi$(jet$_{1,2}$,MET)$>0.4$", 
							Pickles["%s_deltaPhiCut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_deltaPhiCut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_deltaPhiCut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += lineTemplate2%("MET $> 150$ GeV", 
							Pickles["%s_met150Cut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_met150Cut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_met150Cut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							)
	
	table += lineTemplate2%("MT2 $> 80$ GeV", 
							Pickles["%s_signalRegionCut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_noMllCut_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += horizontalLine
	table += "\multicolumn{2}{c}{ \textbf{ttbar like}} \\"+"\n"
	table += horizontalLine
	table += lineTemplate2%("20 $<$ m$_{\ell\ell} <$ 60 GeV", 
							Pickles["%s_signalRegionCut_20To60_lowNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_20To60_lowNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_20To60_lowNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += lineTemplate2%("60 $<$ m$_{\ell\ell} <$ 86 GeV", 
							Pickles["%s_signalRegionCut_60To86_lowNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_60To86_lowNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_60To86_lowNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += lineTemplate2%("96 $<$ m$_{\ell\ell} <$ 150 GeV", 
							Pickles["%s_signalRegionCut_96To150_lowNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_96To150_lowNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_96To150_lowNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += lineTemplate2%("150 $<$ m$_{\ell\ell} <$ 200 GeV", 
							Pickles["%s_signalRegionCut_150To200_lowNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_150To200_lowNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_150To200_lowNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += lineTemplate2%("200 $<$ m$_{\ell\ell} <$ 300 GeV", 
							Pickles["%s_signalRegionCut_200To300_lowNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_200To300_lowNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_200To300_lowNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += lineTemplate2%("300 $<$ m$_{\ell\ell} <$ 400 GeV", 
							Pickles["%s_signalRegionCut_300To400_lowNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_300To400_lowNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_300To400_lowNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += lineTemplate2%("m$_{\ell\ell} >$ 400 GeV", 
							Pickles["%s_signalRegionCut_Above400_lowNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_Above400_lowNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_Above400_lowNll"%(fileName)]["Signal"]["SFYield"],
							)
	table += horizontalLine						
	table += "\multicolumn{2}{c}{ \textbf{non ttbar like}} \\"+"\n"
	table += horizontalLine
	table += lineTemplate2%("20 $<$ m$_{\ell\ell} <$ 60 GeV", 
							Pickles["%s_signalRegionCut_20To60_highNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_20To60_highNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_20To60_highNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += lineTemplate2%("60 $<$ m$_{\ell\ell} <$ 86 GeV", 
							Pickles["%s_signalRegionCut_60To86_highNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_60To86_highNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_60To86_highNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += lineTemplate2%("96 $<$ m$_{\ell\ell} <$ 150 GeV", 
							Pickles["%s_signalRegionCut_96To150_highNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_96To150_highNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_96To150_highNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += lineTemplate2%("150 $<$ m$_{\ell\ell} <$ 200 GeV", 
							Pickles["%s_signalRegionCut_150To200_highNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_150To200_highNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_150To200_highNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += lineTemplate2%("200 $<$ m$_{\ell\ell} <$ 300 GeV", 
							Pickles["%s_signalRegionCut_200To300_highNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_200To300_highNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_200To300_highNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += lineTemplate2%("300 $<$ m$_{\ell\ell} <$ 400 GeV", 
							Pickles["%s_signalRegionCut_300To400_highNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_300To400_highNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_300To400_highNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += lineTemplate2%("m$_{\ell\ell} >$ 400 GeV", 
							Pickles["%s_signalRegionCut_Above400_highNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_Above400_highNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_Above400_highNll"%(fileName)]["Signal"]["SFYield"],
							)	
								
	

	saveTable(tableTemplate2%(table), tableName2)
	
	



main()
