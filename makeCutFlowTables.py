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
	
	#~ Slepton: 700-200, 700-500, 750-300, 800-175, 800-600
	m_sbottom = 800
	M_SBOTTOM = "m_b_"+str(m_sbottom)
	m_neutralino_2 = 200
	
	fileName = "T6bbllslepton_msbottom_%s_mneutralino_%s"%(m_sbottom,m_neutralino_2)
	tableName = "CutFlow_T6bbllslepton_msbottom_%s_mneutralino_%s"%(m_sbottom,m_neutralino_2)
	tableName2 = "CutFlow_T6bbllslepton_msbottom_%s_mneutralino_%s_2"%(m_sbottom,m_neutralino_2)
		
	CutString = ["basicCut","signalRegionCut","met150Cut","nJets2Cut"]
	MllCuts = ["noMllCut","lowMll","highMll"]
	NllCuts = ["inclusiveNll", "lowNll", "highNll"]
	
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
		for MllCut in MllCuts:
			for NllCut in NllCuts:
				Pickles["%s_%s_%s_%s"%(fileName,Cut,MllCut,NllCut)] = loadPickles("shelvesCutFlowTables/%s_%s_%s_%s.pkl"%(fileName,Cut,MllCut,NllCut))
	
	tableTemplate =r"""
\begin{tabular}{l|c}
\hline
slepton model, m$_{\tilde{b}}=800$, m$_{\tilde{\chi}_{2}^{0}}=200$ GeV & SF events \\
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
	
						
	table += horizontalLine
	table += horizontalLine
	
	table += lineTemplate%(" \textbf{low mass}", 
							Pickles["%s_basicCut_lowMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_basicCut_lowMll_inclusiveNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_basicCut_lowMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_basicCut_lowMll_inclusiveNll"%(fileName)]["Signal"]["SystUncertainty"]*Pickles["%s_basicCut_lowMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += horizontalLine
	table += horizontalLine
	
	table += lineTemplate%("jet-multiplicity $\geq 2$", 
							Pickles["%s_nJets2Cut_lowMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_nJets2Cut_lowMll_inclusiveNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_nJets2Cut_lowMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_nJets2Cut_lowMll_inclusiveNll"%(fileName)]["Signal"]["SystUncertainty"]*Pickles["%s_nJets2Cut_lowMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							)
	
	table += horizontalLine
							
	table += lineTemplate%("MET $> 150$ GeV", 
							Pickles["%s_met150Cut_lowMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_met150Cut_lowMll_inclusiveNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_met150Cut_lowMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_met150Cut_lowMll_inclusiveNll"%(fileName)]["Signal"]["SystUncertainty"]*Pickles["%s_met150Cut_lowMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += horizontalLine
	
	table += lineTemplate%("Signal region", 
							Pickles["%s_signalRegionCut_lowMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_lowMll_inclusiveNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_lowMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_lowMll_inclusiveNll"%(fileName)]["Signal"]["SystUncertainty"]*Pickles["%s_signalRegionCut_lowMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += horizontalLine
	table += horizontalLine
	
	table += lineTemplate%("ttbar like", 
							Pickles["%s_signalRegionCut_lowMll_lowNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_lowMll_lowNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_lowMll_lowNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_lowMll_lowNll"%(fileName)]["Signal"]["SystUncertainty"]*Pickles["%s_signalRegionCut_lowMll_lowNll"%(fileName)]["Signal"]["SFYield"],
							)
	table += horizontalLine
	
	table += lineTemplate%("non-ttbar like", 
							Pickles["%s_signalRegionCut_lowMll_highNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_lowMll_highNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_lowMll_highNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_lowMll_highNll"%(fileName)]["Signal"]["SystUncertainty"]*Pickles["%s_signalRegionCut_lowMll_highNll"%(fileName)]["Signal"]["SFYield"],
							)
	
	
	table += horizontalLine
	table += horizontalLine
	
	table += lineTemplate%(" \textbf{high mass}", 
							Pickles["%s_basicCut_highMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_basicCut_highMll_inclusiveNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_basicCut_highMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_basicCut_highMll_inclusiveNll"%(fileName)]["Signal"]["SystUncertainty"]*Pickles["%s_basicCut_highMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += horizontalLine
	table += horizontalLine
	
	table += lineTemplate%("jet-multiplicity $\geq 2$", 
							Pickles["%s_nJets2Cut_highMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_nJets2Cut_highMll_inclusiveNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_nJets2Cut_highMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_nJets2Cut_highMll_inclusiveNll"%(fileName)]["Signal"]["SystUncertainty"]*Pickles["%s_nJets2Cut_highMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							)
	
	table += horizontalLine
							
	table += lineTemplate%("MET $> 150$ GeV", 
							Pickles["%s_met150Cut_highMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_met150Cut_highMll_inclusiveNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_met150Cut_highMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_met150Cut_highMll_inclusiveNll"%(fileName)]["Signal"]["SystUncertainty"]*Pickles["%s_met150Cut_highMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += horizontalLine
	
	table += lineTemplate%("Signal region", 
							Pickles["%s_signalRegionCut_highMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_highMll_inclusiveNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_highMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_highMll_inclusiveNll"%(fileName)]["Signal"]["SystUncertainty"]*Pickles["%s_signalRegionCut_highMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += horizontalLine
	table += horizontalLine
	
	table += lineTemplate%("ttbar like", 
							Pickles["%s_signalRegionCut_highMll_lowNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_highMll_lowNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_highMll_lowNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_highMll_lowNll"%(fileName)]["Signal"]["SystUncertainty"]*Pickles["%s_signalRegionCut_highMll_lowNll"%(fileName)]["Signal"]["SFYield"],
							)
	table += horizontalLine
	
	table += lineTemplate%("non-ttbar like", 
							Pickles["%s_signalRegionCut_highMll_highNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_highMll_highNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_highMll_highNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_highMll_highNll"%(fileName)]["Signal"]["SystUncertainty"]*Pickles["%s_signalRegionCut_highMll_highNll"%(fileName)]["Signal"]["SFYield"],
							)
	
							

	saveTable(tableTemplate%(table), tableName)
	
	
	tableTemplate2 =r"""
\begin{tabular}{l|c}
\hline
slepton model, m$_{\tilde{b}}=800$, m$_{\tilde{\chi}_{2}^{0}}=200$ GeV & SF events \\
 \\
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
	
						
	table += horizontalLine
	table += horizontalLine
	
	table += lineTemplate2%(" \textbf{low mass}", 
							Pickles["%s_basicCut_lowMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_basicCut_lowMll_inclusiveNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_basicCut_lowMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += horizontalLine
	table += horizontalLine
	
	table += lineTemplate2%("jet-multiplicity $\geq 2$", 
							Pickles["%s_nJets2Cut_lowMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_nJets2Cut_lowMll_inclusiveNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_nJets2Cut_lowMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							)
	
	table += horizontalLine
							
	table += lineTemplate2%("MET $> 150$ GeV", 
							Pickles["%s_met150Cut_lowMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_met150Cut_lowMll_inclusiveNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_met150Cut_lowMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += horizontalLine
	
	table += lineTemplate2%("Signal region", 
							Pickles["%s_signalRegionCut_lowMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_lowMll_inclusiveNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_lowMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += horizontalLine
	table += horizontalLine
	
	table += lineTemplate2%("ttbar like", 
							Pickles["%s_signalRegionCut_lowMll_lowNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_lowMll_lowNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_lowMll_lowNll"%(fileName)]["Signal"]["SFYield"],
							)
	table += horizontalLine
	
	table += lineTemplate2%("non-ttbar like", 
							Pickles["%s_signalRegionCut_lowMll_highNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_lowMll_highNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_lowMll_highNll"%(fileName)]["Signal"]["SFYield"],
							)
	
	
	table += horizontalLine
	table += horizontalLine
	
	table += lineTemplate2%(" \textbf{high mass}", 
							Pickles["%s_basicCut_highMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_basicCut_highMll_inclusiveNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_basicCut_highMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += horizontalLine
	table += horizontalLine
	
	table += lineTemplate2%("jet-multiplicity $\geq 2$", 
							Pickles["%s_nJets2Cut_highMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_nJets2Cut_highMll_inclusiveNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_nJets2Cut_highMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							)
	
	table += horizontalLine
							
	table += lineTemplate2%("MET $> 150$ GeV", 
							Pickles["%s_met150Cut_highMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_met150Cut_highMll_inclusiveNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_met150Cut_highMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += horizontalLine
	
	table += lineTemplate2%("Signal region", 
							Pickles["%s_signalRegionCut_highMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_highMll_inclusiveNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_highMll_inclusiveNll"%(fileName)]["Signal"]["SFYield"],
							)
							
	table += horizontalLine
	table += horizontalLine
	
	table += lineTemplate2%("ttbar like", 
							Pickles["%s_signalRegionCut_highMll_lowNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_highMll_lowNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_highMll_lowNll"%(fileName)]["Signal"]["SFYield"],
							)
	table += horizontalLine
	
	table += lineTemplate2%("non-ttbar like", 
							Pickles["%s_signalRegionCut_highMll_highNll"%(fileName)]["Signal"]["SFYield"],
							Pickles["%s_signalRegionCut_highMll_highNll"%(fileName)]["Signal"]["StatUncertainty"]*Pickles["%s_signalRegionCut_highMll_highNll"%(fileName)]["Signal"]["SFYield"],
							)
	
							

	saveTable(tableTemplate2%(table), tableName2)
	
	



main()
