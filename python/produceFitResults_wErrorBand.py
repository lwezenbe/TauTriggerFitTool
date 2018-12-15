from ROOT import *
from array import array
gROOT.SetBatch(True)
from math import sqrt


files = ("../data/NTuple_Data_Run2017BCDEF_31Mar2018_SSsubtraction_VVLooseWP2017v2.root", "../data/NTuple_DYJetsToLL_12Apr2018_v1Andext1v1_12062018_puWeightsANDtauEScorrectionIncluded_OStauGenMatched_VVLooseWP2017v2.root") 

gStyle.SetFrameLineWidth(1)
gStyle.SetPadBottomMargin(0.13)
gStyle.SetPadLeftMargin(0.15)
gStyle.SetPadTopMargin(0.09)
gStyle.SetPadRightMargin(0.05)

triggers = ["ditau", "mutau", "etau"]
types = ["DATA", "MC"]
WPs = ["vvlooseTauMVA", "vlooseTauMVA", "looseTauMVA", "mediumTauMVA", "tightTauMVA", "vtightTauMVA", "vvtightTauMVA"]
tauDMs = ["dm0", "dm1", "dm10"]


hPtDen = [[],[],[]]
hPtNum = [[],[],[]] 

hPtDenDM0 = [[],[],[]]
hPtNumDM0 = [[],[],[]] 
hPtDenDM1 = [[],[],[]]
hPtNumDM1 = [[],[],[]] 
hPtDenDM10 = [[],[],[]]
hPtNumDM10 = [[],[],[]] 

	
for ipath, trigger in enumerate(triggers):
	
	if (ipath == 0 ): #ditau
		edges =[]
		for i in range( 20, 50, 2) :
			edges.append( float(i) )
		for i in range( 50, 75, 5 ) :
			edges.append( float(i) )
		for i in range( 75, 100, 25) :
   			edges.append( float(i) )
		for i in range(100, 150, 50 ) :
			edges.append( float(i) )
		for i in range(150, 500, 300 ) :
 	  	 	edges.append( float(i) )
	if (ipath == 1 ): #mutau
		edges =[]
		for i in range( 20, 50, 2) :
			edges.append( float(i) )
		for i in range( 50, 75, 5 ) :
			edges.append( float(i) )
		for i in range( 75, 100, 25) :
   			edges.append( float(i) )
		for i in range(100, 150, 50 ) :
			edges.append( float(i) )
		for i in range(150, 500, 300 ) :
 	  	 	edges.append( float(i) )
	if (ipath == 2): #etau
		edges =[]
		for i in range( 20, 30, 5) :
			edges.append( float(i) )
		for i in range( 30, 40, 10) :
			edges.append( float(i) )
		for i in range( 40, 100, 20) :
   			edges.append( float(i) )
   		for i in range(100, 140, 40 ) :
			edges.append( float(i) )
		for i in range(140, 500, 310 ) :
			edges.append( float(i) )
	

	print trigger,"trigger", "binning: N=",  len(edges), edges
	
	hPtNum.append([])
	hPtDen.append([])
	# per DM
	hPtNumDM0.append([])
	hPtDenDM0.append([])
	hPtNumDM1.append([])
	hPtDenDM1.append([])
	hPtNumDM10.append([])
	hPtDenDM10.append([])
	
	for index, typ in enumerate(types):
		hPtNum[ipath].append([])
		hPtDen[ipath].append([])
		# per DM	
		hPtNumDM0[ipath].append([])
		hPtDenDM0[ipath].append([])
		hPtNumDM1[ipath].append([])
		hPtDenDM1[ipath].append([])
		hPtNumDM10[ipath].append([])
		hPtDenDM10[ipath].append([])
		
		for ind, wp in enumerate(WPs):
			histoname = "histo_" + typ + "_" + wp + "_" + trigger
			
			hPtDen[ipath][index].append(TH1F (histoname + "_Den", "", len(edges)-1, array('f',edges)))
			hPtNum[ipath][index].append(TH1F (histoname + "_Num", "", len(edges)-1, array('f',edges)))
			# per DM
			hPtDenDM0[ipath][index].append(TH1F (histoname + "_DenDM0", "", len(edges)-1, array('f',edges)))
			hPtNumDM0[ipath][index].append(TH1F (histoname + "_NumDM0", "", len(edges)-1, array('f',edges)))
			hPtDenDM1[ipath][index].append(TH1F (histoname + "_DenDM1", "", len(edges)-1, array('f',edges)))
			hPtNumDM1[ipath][index].append(TH1F (histoname + "_NumDM1", "", len(edges)-1, array('f',edges)))
			hPtDenDM10[ipath][index].append(TH1F (histoname + "_DenDM10", "", len(edges)-1, array('f',edges)))
			hPtNumDM10[ipath][index].append(TH1F (histoname + "_NumDM10", "", len(edges)-1, array('f',edges)))

# This function is taken from Tyler Ruggles SF tool
# https://github.com/truggles/TauTriggerSFs2017/blob/master/python/helpers.py#L9-L40
# Function to create TH1Fs from TGraphAsymmErrors
# This does not preserve the asymmetric errors, only
# bin width and value and does a rough approximation
# on symmetric errors.
def getTH1FfromTGraphAsymmErrors( asym, name ) :

    # Holding vals for TH1F binning and y-vals
    xSpacing = array( 'd', [] )
    yVals = array( 'd', [] )
    yErrors = array( 'd', [] )

    nVals = asym.GetN()
    x, y = Double(0.), Double(0.)
    xEPlus, xEMin = 0., 0.
    yEPlus, yEMin = 0., 0.

    for n in range( nVals ) :
        asym.GetPoint( n, x, y )
        xEPlus = asym.GetErrorXhigh( n )
        xEMin = asym.GetErrorXlow( n )
        yEPlus = asym.GetErrorYhigh( n )
        yEMin = asym.GetErrorYlow( n )
        xSpacing.append( x-xEMin )
        yVals.append( y )
        # To simplify, take asymm errors and go to approximation
        # of symmetric for TH1
        yErrors.append( sqrt(yEPlus**2 + yEMin**2) )

    # Don't forget to add the high end of last bin
    xSpacing.append( x+xEPlus )
    
    outH = TH1F( name, name, len(xSpacing)-1, xSpacing )
    for bin in range( 1, outH.GetNbinsX()+1 ) :
        outH.SetBinContent( bin, yVals[bin-1] )
        outH.SetBinError( bin, yErrors[bin-1] )
    return outH
    
def getScaleFactor(func):

	SF = TGraphAsymmErrors()
	for i in range(20, 450):
		SF.SetPoint(i, i, (func[0].Eval(i)/func[1].Eval(i)))
		
	SF.Draw("A*")
	SF.GetXaxis().SetLimits(18,600)
	SF.GetXaxis().SetMoreLogLabels()
	SF.SetMarkerStyle(20)
	SF.GetXaxis().SetTitle("Offline p_{T}^{#tau} [GeV]")
	SF.GetYaxis().SetTitle("SF: Data/MC")
	return SF

 
for index, filename in enumerate(files):
	print  "filename", filename
	
	file = TFile.Open(filename)
	tree = file.Get('TagAndProbe')
	triggerNamesTree = file.Get("triggerNames")	
	
	print "Populating histograms"

	Nevts = 0
	for iEv in range (0, tree.GetEntries()):
		tree.GetEntry(iEv)
			
		tauPt = tree.tauPt
		HLTPt = tree.hltPt
		tauEta = tree.tauEta
		tauDM = tree.tauDM
		isOS = tree.isOS
		Nvtx = tree.Nvtx
		
		if ("MC" in filename):
			puweight = tree.puweight
		else:
			puweight=1
		
		vvlooseWP = tree.byVVLooseIsolationMVArun2017v2DBoldDMwLT2017	
		vlooseWP = tree.byVLooseIsolationMVArun2017v2DBoldDMwLT2017	
		looseWP = tree.byLooseIsolationMVArun2017v2DBoldDMwLT2017	
		mediumWP = tree.byMediumIsolationMVArun2017v2DBoldDMwLT2017	
		tightWP = tree.byTightIsolationMVArun2017v2DBoldDMwLT2017
		vtightWP = tree.byVTightIsolationMVArun2017v2DBoldDMwLT2017	
		vvtightWP = tree.byVVTightIsolationMVArun2017v2DBoldDMwLT2017	
		
		hasHLTmutauPath_13 = tree.hasHLTmutauPath_13
		hasHLTetauPath_13 = tree.hasHLTetau_Path_13
		hasHLTditauPath_9or10or11 = tree.hasHLTditauPath_9or10or11
		
		Nevents = tree.EventNumber
		Nevts =Nevts + 1
		
		#bkgSubW = 1. if tree.isOS else -1.
		weight = tree.bkgSubW*puweight
		
		HLTpaths = [hasHLTditauPath_9or10or11, hasHLTmutauPath_13 , hasHLTetauPath_13 ]
		WPoints = [vvlooseWP, vlooseWP, looseWP, mediumWP, tightWP, vtightWP, vvtightWP]
		
		for WPind, WP in enumerate(WPoints):
			if(WP > 0):
				for ipath, trigger in enumerate(HLTpaths):
					hPtDen[ipath][index][WPind].Fill(tauPt, weight)
					if(tauDM ==0):
						hPtDenDM0[ipath][index][WPind].Fill(tauPt, weight)
					elif(tauDM ==1):
						hPtDenDM1[ipath][index][WPind].Fill(tauPt, weight)
					elif(tauDM ==10):
						hPtDenDM10[ipath][index][WPind].Fill(tauPt, weight)
				if ( HLTpaths[0] >0 ):
					hPtNum[0][index][WPind].Fill(tauPt, weight)
					if(tauDM ==0):
						hPtNumDM0[0][index][WPind].Fill(tauPt, weight)
					elif(tauDM ==1):
						hPtNumDM1[0][index][WPind].Fill(tauPt, weight)
					elif(tauDM ==10):
						hPtNumDM10[0][index][WPind].Fill(tauPt, weight)
				if( HLTpaths[1] >0 ): 
					hPtNum[1][index][WPind].Fill(tauPt, weight)
					if(tauDM ==0):
						hPtNumDM0[1][index][WPind].Fill(tauPt, weight)
					elif(tauDM ==1):
						hPtNumDM1[1][index][WPind].Fill(tauPt, weight)
					elif(tauDM ==10):
						hPtNumDM10[1][index][WPind].Fill(tauPt, weight)
				if( HLTpaths[2] >0 ):
					hPtNum[2][index][WPind].Fill(tauPt, weight)
					if(tauDM ==0):
						hPtNumDM0[2][index][WPind].Fill(tauPt, weight)
					elif(tauDM ==1):
						hPtNumDM1[2][index][WPind].Fill(tauPt, weight)
					elif(tauDM ==10):
						hPtNumDM10[2][index][WPind].Fill(tauPt, weight)
				


# definition of fit function and its initial parameters
f1 = []
f2 = [[],[]]
for index, typ in enumerate(types):
	f1.append(TF1( 'f1'+typ, '[5] - ROOT::Math::crystalball_cdf(-x, [0], [1], [2], [3])*([4])' ))
	if(index ==0):
		f1[index].SetLineColor( kBlue)
	else:
		f1[index].SetLineColor( kRed)
	f1[index].SetParName( 0, "alpha" )
	f1[index].SetParName( 1, "n" )
	f1[index].SetParName( 2, "simga" )
	f1[index].SetParName( 3, "x0" )
	f1[index].SetParName( 4, "scale" )
	f1[index].SetParName( 5, "y-rise" )

for idm, DM in enumerate(tauDMs):
	f2.append([])
	for index, typ in enumerate(types):
		f2[idm].append(TF1( 'f2_'+ DM  +"_" + typ, '[5] - ROOT::Math::crystalball_cdf(-x, [0], [1], [2], [3])*([4])' ))
		if(idm ==0): f2[idm][index].SetLineColor( kBlue )
		elif(idm ==1): f2[idm][index].SetLineColor( kRed )
		elif(idm ==2): f2[idm][index].SetLineColor( kGreen+3 )		
		f2[idm][index].SetParName( 0, "alpha" )
		f2[idm][index].SetParName( 1, "n" )
		f2[idm][index].SetParName( 2, "simga" )
		f2[idm][index].SetParName( 3, "x0" )
		f2[idm][index].SetParName( 4, "scale" )
		f2[idm][index].SetParName( 5, "y-rise" )

#creates the ratio pad in the same canvas
def createCanvasPads():
    c = TCanvas("c", "canvas", 800, 800)
    # Upper histogram plot is pad1
    pad1 = TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0.1)  # joins upper and lower plot
    pad1.SetLeftMargin( gPad.GetLeftMargin() * 1.5 )
    pad1.SetRightMargin( gPad.GetRightMargin() * 1.5 )
    pad1.SetGridx()
    pad1.Draw()
    # Lower ratio plot is pad2
    c.cd()  # returns to main canvas before defining pad2
    pad2 = TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
    pad2.SetTopMargin(0.1)  # joins upper and lower plot
    pad2.SetBottomMargin(0.2)
    pad2.SetGridx()
    pad2.SetLeftMargin( gPad.GetLeftMargin() * 1.5 )
    pad2.SetRightMargin( gPad.GetRightMargin() * 1.5 )
    pad2.Draw()
    return c, pad1, pad2


def createRelativeErrors(fitresult, g_efficiency, f1):

	# confidence interval#
	values = fitresult.GetConfidenceIntervals(0.68, False)
	interval = TGraphErrors(g_efficiency.GetN())
	ratio = TGraphAsymmErrors(g_efficiency.GetN())
	ratio2 = TGraphAsymmErrors(g_efficiency.GetN())
	for i in range(0, g_efficiency.GetN()):
   		interval.SetPoint(i, g_efficiency.GetX()[i], f1[index].Eval(g_efficiency.GetX()[i] ))
		interval.SetPointError(i, 0, values[i] )
		ratio.SetPoint(i, g_efficiency.GetX()[i], (f1[index].Eval(g_efficiency.GetX()[i]) - values[i])/f1[index].Eval(g_efficiency.GetX()[i]))
		ratio2.SetPoint(i, g_efficiency.GetX()[i], (f1[index].Eval(g_efficiency.GetX()[i]) + values[i])/f1[index].Eval(g_efficiency.GetX()[i]))
		
	ratio.GetXaxis().SetTitleSize(0.05)
	ratio.GetYaxis().SetTitleSize(0.05)
	ratio.GetXaxis().SetTitleOffset(1.1)
	ratio.GetYaxis().SetTitleOffset(1.1)
	ratio.SetLineWidth(2)
	ratio2.SetLineWidth(2)
	ratio.Draw("A*")
	ratio2.Draw("*same")
	ratio.SetMarkerStyle(20)
	ratio2.SetMarkerStyle(20)
	ratio.GetYaxis().SetRangeUser(0.8,1.2)
	ratio.SetTitle("")
	ratio.SetTitle(trigger +"Path_" + wp +"_"+ typ)
	ratio.GetYaxis().SetTitle("#frac{fit #pm error}{fit}")
	ratio.GetXaxis().SetTitle("Offline p_{T}^{#tau} [GeV]")
	ratio2.SetMarkerColor(2)
	ratio.GetXaxis().SetMoreLogLabels()
	
	return ratio , ratio2

outputname = "tauTriggerEfficiencies2017_final.root"
file = TFile( "../data/"+outputname, 'recreate')

# efficiency calculation after filling the histograms for 3 different triggers for each WPs of DATA and MC
for ipath, trigger in enumerate(triggers):
	for WPind, wp in enumerate(WPs):
		for index, typ in enumerate(types):
			nbins = hPtNum[ipath][index][WPind].GetNbinsX()
			# checks the content of the bins such that Den > Num
			for binid in range(0, nbins + 1):
				if(hPtNum[ipath][index][WPind].GetBinContent(binid) > hPtDen[ipath][index][WPind].GetBinContent(binid)):
					hPtNum[ipath][index][WPind].SetBinContent(binid, hPtDen[ipath][index][WPind].GetBinContent(binid))
				if(hPtNumDM0[ipath][index][WPind].GetBinContent(binid) > hPtDenDM0[ipath][index][WPind].GetBinContent(binid)):
					hPtNumDM0[ipath][index][WPind].SetBinContent(binid, hPtDenDM0[ipath][index][WPind].GetBinContent(binid))
				if(hPtNumDM1[ipath][index][WPind].GetBinContent(binid) > hPtDenDM1[ipath][index][WPind].GetBinContent(binid)):
					hPtNumDM1[ipath][index][WPind].SetBinContent(binid, hPtDenDM1[ipath][index][WPind].GetBinContent(binid))
				if(hPtNumDM10[ipath][index][WPind].GetBinContent(binid) > hPtDenDM10[ipath][index][WPind].GetBinContent(binid)):
					hPtNumDM10[ipath][index][WPind].SetBinContent(binid, hPtDenDM10[ipath][index][WPind].GetBinContent(binid))
			# setting the fit parameters for various cases
			if(ipath == 0): # ditau
				f1[index].SetParameter( 0, 0.2)
				f1[index].SetParameter( 1, 5.0 )
				f1[index].SetParameter( 2, 7.0 )
				f1[index].SetParameter( 3, -30.)
				f1[index].SetParameter( 4, 1.0 )
				f1[index].SetParameter( 5, 1.0)
				for idm, DM in enumerate(tauDMs):
					if(idm ==0 or idm == 1):
						f2[idm][index].SetParameter( 0, 0.2)
						f2[idm][index].SetParameter( 1, 5.0 )
						f2[idm][index].SetParameter( 2, 7.0 )
						f2[idm][index].SetParameter( 3, -30.)
						f2[idm][index].SetParameter( 4, 1.0 )
						f2[idm][index].SetParameter( 5, 1.0)
					else:
						f2[idm][index].SetParameter( 0, 0.8)
						f2[idm][index].SetParameter( 1, 5.0)
						f2[idm][index].SetParameter( 2, 7.0)
						f2[idm][index].SetParameter( 3, -30.)
						f2[idm][index].SetParameter( 4, 1.0 )
						f2[idm][index].SetParameter( 5, 1.0)
			if(ipath == 1): # mutau
				f1[index].SetParameter( 0, 0.2)
				f1[index].SetParameter( 1, 5.0 )
				f1[index].SetParameter( 2, 7.0 )
				f1[index].SetParameter( 3, -30.)
				f1[index].SetParameter( 4, 1.0 )
				f1[index].SetParameter( 5, 1.0)
				for idm, DM in enumerate(tauDMs):
					if(idm ==0 or idm == 1):
						f2[idm][index].SetParameter( 0, 0.5)
						f2[idm][index].SetParameter( 1, 5.0 )
						f2[idm][index].SetParameter( 2, 7.0 )
						f2[idm][index].SetParameter( 3, -20.)
						f2[idm][index].SetParameter( 4, 1.0 )
						f2[idm][index].SetParameter( 5, 1.0)
					else:
						f2[idm][index].SetParameter( 0, 0.8)
						f2[idm][index].SetParameter( 1, 10.0 )
						f2[idm][index].SetParameter( 2, 7.0 )
						f2[idm][index].SetParameter( 3, -25.)
						f2[idm][index].SetParameter( 4, 1.0 )
						f2[idm][index].SetParameter( 5, 1.0)				
					
			if(ipath == 2): # etau
				f1[index].SetParameter( 0, 0.8)
				f1[index].SetParameter( 1, 8.0 )
				f1[index].SetParameter( 2, 7.0 )
				f1[index].SetParameter( 3, -20.)
				f1[index].SetParameter( 4, 1.0)
				f1[index].SetParameter( 5, 1.0)
				for idm, DM in enumerate(tauDMs):
					f2[idm][index].SetParameter( 0, 1.0)
					f2[idm][index].SetParameter( 1, 10.0 )
					f2[idm][index].SetParameter( 2, 8.0 )
					f2[idm][index].SetParameter( 3, -25.)
					f2[idm][index].SetParameter( 4, 1.0 )
					f2[idm][index].SetParameter( 5, 1.0)
			
			g_efficiency =TGraphAsymmErrors()
			h_errBand68 = TH1F(histoname+"_CL68","histo of 0.68 confidence band", 480, 20, 500)
			g_errBand68 = TGraphErrors()
			
			g_efficiency.BayesDivide(hPtNum[ipath][index][WPind],hPtDen[ipath][index][WPind]) #,"cl=0.683 b(1,1) mode")		
			h_efficiency = getTH1FfromTGraphAsymmErrors(g_efficiency,"histo_" + trigger + "ErrorBand_" + wp +"_"+ typ )
					
			# write the histograms/graphs into the output ROOT file before the fit
			h_efficiency.Write("histo_"+ trigger +"Efficiency_" + wp +"_"+ typ)
			g_efficiency.Write("graph_"+ trigger +"Efficiency_" + wp +"_"+ typ)
			
			print "Fit is performed for", trigger, "trigger in", wp ,"WP for", typ 
			print "Fit parameters:", f1[index].GetParameter(0), f1[index].GetParameter(1), f1[index].GetParameter(2), f1[index].GetParameter(3), f1[index].GetParameter(4), f1[index].GetParameter(5)    
			fit_result = g_efficiency.Fit('f1'+ typ, 'S')
			
			TVirtualFitter.GetFitter().GetConfidenceIntervals(h_errBand68, 0.68)
			
			for i in range(0, g_efficiency.GetN()):
				g_errBand68.SetPoint(i, g_efficiency.GetX()[i], 0)
			TVirtualFitter.GetFitter().GetConfidenceIntervals(g_errBand68, 0.68)

			# Set the title of the histograms/graphs and their axes
			g_efficiency.SetTitle(trigger +"Path_" + wp +"_"+ typ)
			g_efficiency.GetYaxis().SetTitle("Efficiency")
			g_efficiency.GetXaxis().SetTitle("Offline p_{T}^{#tau} [GeV]")
			g_efficiency.SetTitle(trigger +"Path_" + wp +"_"+ typ)
			g_efficiency.GetYaxis().SetTitle("Efficiency")
			g_efficiency.GetXaxis().SetTitle("Offline p_{T}^{#tau} [GeV]")
			h_efficiency.SetTitle(trigger +"Path_" + wp +"_"+ typ)
			h_efficiency.GetYaxis().SetTitle("Efficiency")
			h_efficiency.GetXaxis().SetTitle("Offline p_{T}^{#tau} [GeV]")
			h_errBand68.SetTitle(trigger +"Path_" + wp +"_"+ typ)
			h_errBand68.GetYaxis().SetTitle("Efficiency")
			h_errBand68.GetXaxis().SetTitle("Offline p_{T}^{#tau} [GeV]")
			g_errBand68.SetTitle(trigger +"Path_" + wp +"_"+ typ)
			g_errBand68.GetYaxis().SetTitle("Efficiency")
			g_errBand68.GetXaxis().SetTitle("Offline p_{T}^{#tau} [GeV]")
			
			# write the histograms/graphs into the output ROOT file after the fit
			g_efficiency.Write("grFit_"+ trigger +"Efficiency_" + wp +"_"+ typ)
			h_errBand68.Write("histo_"+ trigger +"ErrorBand_" + wp +"_"+ typ)
			g_errBand68.Write("graph_"+ trigger +"ErrorBand_" + wp +"_"+ typ)
			fit_result.Write("fitResult_"+ trigger +"Path_" + wp +"_"+ typ)
			

			#======== Relative error of the fit: "fit +/- error/ fit " ===================
			relativeError = TGraphAsymmErrors()
			relativeErrorUP = TGraphAsymmErrors()
			relativeErrorDown = TGraphAsymmErrors()
			relativeErrorUP, relativeErrorDown = createRelativeErrors(fit_result, g_efficiency, f1)
			relativeErrorUP.Draw("P")
			relativeErrorDown.Draw("P")
			relativeErrorUP.Write("relativeErrorUp_"+ trigger  + wp +"_"+ typ)
			relativeErrorDown.Write("relativeErrorDown_"+ trigger  + wp +"_"+ typ)
			
			# per DM efficiencies
			g_efficiencyDM0 =TGraphAsymmErrors()
			g_efficiencyDM1 =TGraphAsymmErrors()
			g_efficiencyDM10 =TGraphAsymmErrors()
			g_efficiencyDM0.BayesDivide(hPtNumDM0[ipath][index][WPind],hPtDenDM0[ipath][index][WPind]) 		
			g_efficiencyDM1.BayesDivide(hPtNumDM1[ipath][index][WPind],hPtDenDM1[ipath][index][WPind])		
			g_efficiencyDM10.BayesDivide(hPtNumDM10[ipath][index][WPind],hPtDenDM10[ipath][index][WPind]) 		

			
			#print "Fit is performed for", trigger, "trigger in", wp ,"WP for", typ , " per DM"
			fit_result2 = g_efficiencyDM0.Fit('f2_dm0'+"_" + typ, 'S')
			fit_result3 = g_efficiencyDM1.Fit('f2_dm1'+ "_" +typ, 'S')	
			fit_result4 = g_efficiencyDM10.Fit('f2_dm10' +"_" + typ, 'S')
			
			g_efficiencyDM0.Write("graph_"+ trigger +"Efficiency_dm0_" + wp +"_"+ typ)
			g_efficiencyDM1.Write("graph_"+ trigger +"Efficiency_dm1_" + wp +"_"+ typ)
			g_efficiencyDM10.Write("graph_"+ trigger +"Efficiency_dm10_" + wp +"_"+ typ)
			
		# Getting Scale Factors
		SF = TGraphAsymmErrors()
		SF = getScaleFactor(f1)
		SF.Write( 'ScaleFactor_'+ trigger + '_' + wp + '_DataMC')
		
		# Getting Scale Factors per decay mode
		for idm, DM in enumerate(tauDMs):
			SF_dm = TGraphErrors()
			SF_dm = getScaleFactor(f2[idm])
			SF_dm.Write( 'ScaleFactor_'+ trigger + "_" + DM + '_' + wp + '_DataMC')
	  
file.Close()
print "The output ROOT file has been created: ../data/" + outputname

