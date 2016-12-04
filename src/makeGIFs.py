import os
import sys
import optparse
import ROOT
import commands
import getpass
import pickle
import numpy

from PIL import Image, ImageSequence
from images2gif import writeGif

ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
#ROOT.gROOT.SetBatch()

def replaceBadCharacters(inputStr):
    newStr = inputStr
    for token in ['+','-','*',' ','#','{','(',')','}','@']:
        newStr=newStr.replace(token,'')
    return newStr

"""
save distributions to file
"""
def saveToShapesFile(outFile,shapeColl,directory=''):
    fOut=ROOT.TFile.Open(outFile,'UPDATE')
    if len(directory)==0:
        fOut.cd()
    else:
        if not fOut.Get(directory):
            fOut.mkdir(directory)
        outDir=fOut.Get(directory)
        outDir.cd()
    for key in shapeColl:
        #remove bin labels
        shapeColl[key].GetXaxis().Clear()

        #convert to TH1D (projections are TH1D)
        if not shapeColl[key].InheritsFrom('TH1D') :
            h=ROOT.TH1D()
            shapeColl[key].Copy(h)
            shapeColl[key]=h

        shapeColl[key].Write(key,ROOT.TObject.kOverwrite)
    fOut.Close()

"""
steer the script
"""

#configuration
usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-i', '--input',     dest='input',     help='input plotter',                            default=None,            type='string')
parser.add_option('-o', '--output',    dest='output',    help='output directory',                         default='datacards',     type='string')
parser.add_option('-n', '--outname',   dest='outname',   help='output file name',                         default='shapes',     type='string')
parser.add_option('--makesens', dest='makeSens',         help='make local sensitivity validation (very slow)', default=False, action='store_true')
parser.add_option('--novalidation', dest='noValidation', help='do not make any validation plots',              default=False, action='store_true')
(opt, args) = parser.parse_args()

def make3DPlotGIF(plot,n3DScan=120,text="") :
    canvas=ROOT.TCanvas("","",600,800)
    canvas.cd();
    ROOT.gStyle.SetPalette(ROOT.kBird)
    rotScanImgs=[]
    status=ROOT.TLatex()

    for theta in xrange(0,359,int(360/n3DScan)) :
        # mlb vs. mass plot for 3D scan
        plot.Draw("SURF1")
        status.DrawLatexNDC(0.75,0.02,text)
        ROOT.gPad.SetPhi(theta)

        canvas.SaveAs("%s/%s_%i.png"%(opt.output,opt.theta))
        canvas.Clear()

        # save images to proper arrays
        img=Image.open("%s/%s_%i.png"%(opt.output,plot.getName(),theta))
        rotScanImgs+=[img]

    writeGif("%s/%s.gif"%(opt.output,plot.getName()),rotScanImgs,duration=.05)



#def main():
#
#    # prepare output directory
#    os.system('mkdir -p %s'%opt.output)
#
#    # get data and nominal expectations
#    fIn=ROOT.TFile.Open(opt.input)
#    fIn.Get._creates=True
#    systfIn=None
#    if opt.systInput:
#        systfIn=ROOT.TFile.Open(opt.systInput)
#        systfIn.Get._creates=True
#
#    # prepare output ROOT file
#    if not opt.skipMakingShapes :
#        outFile='%s/%s.root'%(opt.output,opt.outname)
#        fOut=ROOT.TFile.Open(outFile,'UPDATE')
#        fOut.Close()
#
#    ##################
#    # BEGIN MORPHING #
#    ##################
#    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
#
#    print "\n Creating morphed dists\n"
#
#    kOrangeList=[-4,-2,0,1,-3,7,8,10,9,2]
#    widths=map(float,[wid.replace('p','.').replace('w','') for wid in rawWidList])
#    massInfo=[(169.5,opt.systInput,'t#bar{t} m=169.5'),
#              (172.5,opt.input    ,'t#bar{t}'),
#              (175.5,opt.systInput,'t#bar{t} m=175.5'),
#              (169.5,opt.systInput,'tW m=169.5'),
#              (172.5,opt.input    ,'tW'),
#              (175.5,opt.systInput,'tW m=175.5')]
#    minMT=169.5
#    nomMT=172.5
#    maxMT=175.5
#    minGammaT=0.5*1.324
#    nomGammaT=1.324
#    maxGammaT=5.0*1.324
#    distToTitle={
#            "mlb":    "Inclusive M_{lb}",
#            "minmlb": "Minimum M_{lb}",
#            "mdrmlb": "#DeltaR-Filtered M_{lb}",
#            "incmlb": "Inclusive M_{lb}",
#            "sncmlb": "Semi-Inclusive M_{lb}",
#            "mt2mlb": "M_{T2}^{lb} Strategy"
#           }
#
#    canvas=ROOT.TCanvas()
#    procHasUnc={}
#
#    # make workspace for each signal process
#    first=True
#    for sig,dist,lbCat,ch,cat in [(a,b,c,d,e)
#            for a in (['tbart'] if not opt.allMorphs else signalList)
#            for b in (['incmlb'] if not opt.allMorphs else distList)
#            for c in (['highpt'] if not opt.allMorphs else lbCatList)
#            for d in (['EM'] if not opt.allMorphs else lfsList)
#            for e in (['2b'] if not opt.allMorphs else catList)]:
#        modSig=replaceBadCharacters(sig)
#        distTitle="M_{lb}" if dist not in distToTitle else distToTitle[dist]
#
#        if sig == signalList[0] :
#            procHasUnc["%s%s%s_%s"%(lbCat,ch,cat,dist)]={}
#
#        # get the proper mass information for this signal
#        masses=[]
#        for mass,fname,distname in massInfo :
#            if modSig in replaceBadCharacters(distname) :
#                masses += [(mass,fname,distname)]
#
#        # create workspace
#        ws=ROOT.RooWorkspace('sigws_%s%s%s_%s_%s'%(lbCat,ch,cat,modSig,dist))
#        ws.factory('x[0,300]')
#
#        ##########################################
#        # produce systematic uncertainties in ws #
#        ##########################################
#        downShapes, upShapes = {}, {}
#
#        #nominal expectations
#        exp={}
#        if opt.addSigs :
#            _,exp=getMergedDists(fIn,('%s%s%s_%s_'%(lbCat,ch,cat,dist)),rawSignalList,'top',widList,nomWid,signalList)
#        else :
#            _,exp=getMergedDists(fIn,('%s%s%s_%s_'%(lbCat,ch,cat,dist)),None,'',widList,nomWid,signalList)
#
#        # sample systematics
#        _,genVarShapes = getMergedDists(fIn,
#                ('%s%s%s_%s_'%(lbCat,ch,cat,dist)),
#                (rawSignalList if opt.addSigs else None),
#                ('top' if opt.addSigs else ''),
#                widList,nomWid,systSignalList)
#        _,altExp       = getMergedDists(systfIn,
#                ('%s%s%s_%s_'%(lbCat,ch,cat,dist)),
#                (rawSignalList if opt.addSigs else None),
#                ('top' if opt.addSigs else ''),
#                widList,nomWid,systSignalList)
#
#        for systVar, procsToApply, normalize, useAltShape, projectRelToNom in sampleSysts:
#            downShapes[systVar]={}
#            upShapes[systVar]={}
#
#            if sig == signalList[0] :
#                procHasUncName="%s%s%s_%s"%(lbCat,ch,cat,dist)
#                procHasUnc[procHasUncName][systVar] = []
#
#            for iproc in procsToApply:
#                if sig not in iproc and sig != signalList[0] : continue
#                nomH=exp[iproc]
#
#                #check which shape to use
#                if useAltShape:
#
#                    #get directly from another file
#                    downH  = altExp[ procsToApply[iproc][0] ]
#                    if len( procsToApply[iproc] ) > 1 :
#                        upH = altExp[ procsToApply[iproc][1] ]
#                    else:
#                        #if only one variation is available, mirror it
#                        upH = downH.Clone( '%s%sUp'%(iproc,systVar) )
#                        for xbin in xrange(1,upH.GetNbinsX()+1):
#                            diff=upH.GetBinContent(xbin)-nomH.GetBinContent(xbin)
#                            upH.SetBinContent(xbin,nomH.GetBinContent(xbin)-diff)
#                else:
#
#                    #project from 2D histo (re-weighted from nominal sample)
#                    ybinUp, ybinDown = -1, -1
#                    for ybin in xrange(1,genVarShapes[ iproc ].GetNbinsY()+1):
#                        label = genVarShapes[ iproc ].GetYaxis().GetBinLabel(ybin)
#                        if procsToApply[iproc][0] in label : ybinDown=ybin
#                        if procsToApply[iproc][1] in label : ybinUp=ybin
#
#                    downH = genVarShapes[ iproc ].ProjectionX('%s%sDown'%(iproc,systVar), ybinDown, ybinDown)
#                    upH   = genVarShapes[ iproc ].ProjectionX('%s%sUp'%(iproc,systVar),   ybinUp,   ybinUp)
#
#                # use do down/up x nom to generate the variation, then mirror it
#                if projectRelToNom:
#                    ratioH=downH.Clone()
#                    ratioH.Divide(upH)
#                    for xbin in xrange(1,nomH.GetNbinsX()+1):
#                        nomVal=nomH.GetBinContent(xbin)
#                        varVal = ratioH.GetBinContent(xbin) * nomVal
#                        upH.SetBinContent(xbin, varVal)
#                        if upH.GetBinContent(xbin) <= 0 :
#                            upH.SetBinContent(xbin,1e-010)
#                        varVal = varVal- nomVal
#                        downH.SetBinContent(xbin, nomVal-varVal)
#                        if downH.GetBinContent(xbin) <= 0 :
#                            downH.SetBinContent(xbin,1e-010)
#
#                #normalize (shape only variation is considered)
#                if normalize : downH.Scale( nomH.Integral()/downH.Integral() )
#                if normalize : upH.Scale( nomH.Integral()/upH.Integral() )
#
#                #check if variation is meaningful
#                accept = acceptVariationForDataCard(nomH=nomH, upH=upH, downH=downH)
#                if not accept : continue
#
#                #save
#                if sig == signalList[0] :
#                    procHasUncName="%s%s%s_%s"%(lbCat,ch,cat,dist)
#                    procHasUnc[procHasUncName][systVar]+=[iproc]
#
#                if sig not in iproc : continue
#                downShapes[systVar][iproc]=downH
#                upShapes[systVar][iproc]=upH
#
#
#
#
#        for systVar, normalize, useAltShape, projectRelToNom in genSysts:
#            dirForSyst=('%sup%s%s%s_%s_'%(systVar,lbCat,lfs,cat,dist))
#            if "MEmuF" in systVar : dirForSyst=('gen3%s%s%s_%s_'%(lbCat,lfs,cat,dist))
#            if "MEmuR" in systVar : dirForSyst=('gen5%s%s%s_%s_'%(lbCat,lfs,cat,dist))
#            if "MEtot" in systVar : dirForSyst=('gen6%s%s%s_%s_'%(lbCat,lfs,cat,dist))
#
#            _,genVarShapesUp = getMergedDists((systfIn if useAltShape else fIn),
#                    dirForSyst,
#                    (rawSignalList if opt.addSigs else None),
#                    ('top' if opt.addSigs else ''),
#                    widList,nomWid,signalList)
#
#            dirForSyst=('%sup%s%s%s_%s_'%(systVar,lbCat,lfs,cat,dist))
#            if "MEmuF" in systVar : dirForSyst=('gen4%s%s%s_%s_'%(lbCat,lfs,cat,dist))
#            if "MEmuR" in systVar : dirForSyst=('gen8%s%s%s_%s_'%(lbCat,lfs,cat,dist))
#            if "MEtot" in systVar : dirForSyst=('gen10%s%s%s_%s_'%(lbCat,lfs,cat,dist))
#
#            _,genVarShapesDn = getMergedDists((systfIn if useAltShape else fIn),
#                    dirForSyst,
#                    (rawSignalList if opt.addSigs else None),
#                    ('top' if opt.addSigs else ''),
#                    widList,nomWid,signalList)
#
#            downShapes[systVar]={}
#            upShapes[systVar]={}
#
#            if sig == signalList[0] :
#                procHasUncName="%s%s%s_%s"%(lbCat,ch,cat,dist)
#                procHasUnc[procHasUncName][systVar] = []
#
#            for iproc in exp:
#                if sig not in iproc and sig != signalList[0] : continue
#                nomH=exp[iproc]
#
#                if not iproc in genVarShapesDn :
#                    continue
#                #get directly from another file
#                downH  = genVarShapesDn[iproc]
#                if iproc in genVarShapesUp:
#                    upH = genVarShapesUp[iproc]
#                else:
#                    #if only one variation is available, mirror it
#                    upH = downH.Clone( '%s%sUp'%(iproc,systVar) )
#                    for xbin in xrange(1,upH.GetNbinsX()+1):
#                        diff=upH.GetBinContent(xbin)-nomH.GetBinContent(xbin)
#                        upH.SetBinContent(xbin,nomH.GetBinContent(xbin)-diff)
#
#                # use do down/up x nom to generate the variation, then mirror it
#                if projectRelToNom:
#                    ratioH=downH.Clone()
#                    ratioH.Divide(upH)
#                    for xbin in xrange(1,nomH.GetNbinsX()+1):
#                        nomVal=nomH.GetBinContent(xbin)
#                        varVal = ratioH.GetBinContent(xbin) * nomVal
#                        upH.SetBinContent(xbin, varVal)
#                        varVal = varVal- nomVal
#                        downH.SetBinContent(xbin, nomVal-varVal)
#
#                #normalize (shape only variation is considered)
#                if normalize : downH.Scale( nomH.Integral()/downH.Integral() )
#                if normalize : upH.Scale( nomH.Integral()/upH.Integral() )
#
#                #check if variation is meaningful
#                accept = acceptVariationForDataCard(nomH=nomH, upH=upH, downH=downH)
#                if not accept : continue
#
#                #save
#                if sig in iproc and 'V' not in iproc:
#                    downShapes[systVar][iproc]=downH
#                    upShapes[systVar][iproc]=upH
#
#                if sig == signalList[0] :
#                    procHasUncName="%s%s%s_%s"%(lbCat,ch,cat,dist)
#                    print procHasUncName,',',systVar,',',iproc,',',upH.Integral()
#                    procHasUnc[procHasUncName][systVar]+=[iproc]
#
#
#
#        ######################
#        # Produce morph grid #
#        ######################
#
#        # initialize grid with width and mass
#        widthDim=ROOT.RooBinning(len(widths),0,len(widths)-1)
#        massDim =ROOT.RooBinning(len(masses),0,len(masses)-1)
#        refGrid =ROOT.RooMomentMorphND.Grid(widthDim,massDim)
#
#        massAxes=[ws.var('x').frame()]*len(widths)
#        widAxes =[ws.var('x').frame()]*len(masses)
#        for imass in xrange(0,len(masses)):
#            mass,url,proc=masses[imass]
#
#            if modSig not in replaceBadCharacters(proc) : continue
#
#            tfIn=ROOT.TFile.Open(url)
#            for iwid in xrange(0,len(widths)):
#
#                #get histogram and convert to a PDF
#                dirname='%s%s%s_%s_%3.1fw'%(lbCat,ch,cat,dist,widths[iwid])
#                h=tfIn.Get(dirname+"/"+dirname+'_'+proc)
#                if not isinstance(h, ROOT.TH1):
#                    print "\t\tProblem with pdf ",dirName,"/",dirName,"_",proc,"  | Continuing..."
#                    continue
#                name='%s_m%d'%(dirname,int(10*mass))
#                data=ROOT.RooDataHist(name,name,ROOT.RooArgList(ws.var("x")),h)
#                pdf=ROOT.RooHistPdf(name+"_pdf",name+"_pdf",ROOT.RooArgSet(ws.var("x")),data)
#                pdf.plotOn(massAxes[iwid],
#                        ROOT.RooFit.Name("%3.1f"%mass),
#                        ROOT.RooFit.LineColor(ROOT.kOrange+kOrangeList[imass*3]))
#                pdf.plotOn(widAxes[imass],
#                        ROOT.RooFit.Name("%3.1f"%widths[iwid]),
#                        ROOT.RooFit.LineColor(ROOT.kOrange+kOrangeList[iwid]))
#                getattr(ws,'import')(pdf,ROOT.RooCmdArg())
#
#                #add pdf to the grid
#                print 'Adding',pdf.GetName(),'@ (',iwid,',',imass,')'
#                insertVector=ROOT.vector('int')(2)
#                insertVector[0]=iwid
#                insertVector[1]=imass
#                refGrid.addPdf(ws.pdf(pdf.GetName()),insertVector)
#
#            tfIn.Close()
#
#        #################################
#        # create pdf and save workspace #
#        #################################
#        ws.factory('alpha[%i,%i]'%(0,len(widths)-1))
#        ws.factory('beta[%i,%i]'%(0,len(masses)-1))
#        ws.var('alpha').setVal(1)
#        ws.var('beta' ).setVal(1)
#        wsArgList=ROOT.RooArgList( ws.var('alpha'), ws.var('beta') )
#        for systVar in upShapes :
#            if len(upShapes[systVar])==0 : continue
#            ws.factory('%s[-1,1]'%(systVar))
#            ws.var('%s'%systVar).setVal(0)
#            wsArgList.add( ws.var('%s'%systVar) , True )
#        pdf=ROOT.RooMomentMorphND('widmorphpdf','widmorphpdf',
#                                  wsArgList,
#                                  ROOT.RooArgList( ws.var('x') ),
#                                  refGrid,
#                                  ROOT.RooMomentMorphND.Linear)
#        pdf.useHorizontalMorphing(False)
#        getattr(ws,'import')(pdf,ROOT.RooCmdArg())
#
#        # save workspace to shapes
#        outFile='%s/shapes.root'%(opt.output)
#        fOut=ROOT.TFile.Open(outFile,'UPDATE')
#        fOut.cd()
#        ws.Write()
#        fOut.Close()
#        print "Closing outfile"
#
#
#        ###############################
#        # Produce morphing animations #
#        ###############################
#        if not opt.allMorphs and not first : continue
#        if opt.noValidation : continue
#        ROOT.gStyle.SetOptStat(0)
#        status=ROOT.TLatex()
#
#        zLimList={ "mlb": 3.5e-04,
#                "asens" : 3.5e-04,
#                "bsens" : 3.5e-04 }
#        n2DScan=60
#        n3DScan=120
#
#        # setup morph validation arrays
#        rotScanAlphaImgs  = []
#        rotScanBetaImgs   = []
#        alphaScanImgs     = []
#        betaScanImgs      = []
#        sensAlphaScanImgs = []
#        sensBetaScanImgs  = []
#
#        print "here1"
#        os.mkdir('%s/gifplots_%s%s%s_%s_%s'%(opt.output,lbCat,ch,cat,sig,dist))
#
#        # NOTE: we want to reopen the outfile so that we see what Combine sees.
#        #       This is inefficient but it is important validation!
#        fOut=ROOT.TFile.Open(outFile)
#        workspace=fOut.Get("sigws_%s%s%s_%s_%s"%(lbCat,ch,cat,modSig,dist))
#        morphPDF= workspace.pdf("widmorphpdf")
#
#        alphaVar= workspace.var("alpha")
#        betaVar = workspace.var("beta")
#        mlbVar  = workspace.var("x")
#
#        # this is used separately, so we can use pdf
#        dMorDAlpha=pdf.derivative(alphaVar)
#        dMorDBeta =pdf.derivative(betaVar)
#        dMorDMlb  =pdf.derivative(mlbVar)
#
#        #import pdb
#        #pdb.set_trace()
#
#        print "Begin morph validation!"
#
#        #################################
#        # create morphed PDF validation #
#        #################################
#        histos=[mlbVar.frame()]*len(widths)
#        for i in range(0,len(widths)) :
#            print "Setting value of alpha!"
#            alphaVar.setVal(widths[i])
#
#            leg=ROOT.TLegend(0.65,0.17,0.90,0.37)
#            leg.SetHeader("#Gamma_{t} = %s#times#Gamma_{SM}"%rawWidList[i])
#            print "Validating morphed pdfs!"
#
#            for imass in range(0,len(masses)) :
#                print "Setting value of beta!"
#                betaVar.setVal(imass)
#                morphPDF.plotOn(histos[i],
#                        ROOT.RooFit.Name("%3.1f"%masses[imass][0]),
#                        ROOT.RooFit.LineColor(ROOT.kOrange+kOrangeList[imass*3]))
#                leg.AddEntry("%3.1f"%(masses[imass][0]),"m_{t} = %3.1f GeV"%(masses[imass][0]),"L")
#
#            histos[i].GetXaxis().SetTitle("%s [GeV]"%distTitle)
#            histos[i].GetYaxis().SetTitle("Morphed density: lepton-jet pairs")
#            histos[i].SetTitle("")
#            histos[i].SetLineWidth(1)
#            histos[i].Draw()
#            leg.Draw()
#
#
#            canvas.SetLogy()
#            canvas.Update()
#            canvas.SaveAs("%s/morphValidation_%s%s%s_%s_%s.pdf"%(opt.output,lbCat,ch,cat,dist,modWidList[i]))
#            canvas.SetGrayscale(True)
#            canvas.Update()
#            canvas.SaveAs("%s/gScale__morphValidation_%s%s%s_%s_%s.pdf"%(opt.output,lbCat,ch,cat,dist,modWidList[i]))
#            canvas.SetGrayscale(False)
#
#        ##########################################
#        # create 3D plots that rotate the camera #
#        ##########################################
#        canvas.SetLogx(False)
#        canvas.SetLogy(False)
#        canvas.SetLogz(True)
#        ROOT.gStyle.SetPalette(52)
#
#        histAlpha = mlbVar.createHistogram("Morphing against width variations",alphaVar)
#        morphPDF.fillHistogram(histAlpha,ROOT.RooArgList(mlbVar,alphaVar))
#        histAlpha.GetYaxis().SetTitleOffset(1.1)
#        histAlpha.GetYaxis().SetTitle("Grid width index")
#        histAlpha.GetXaxis().SetTitleOffset(1.35)
#        histAlpha.GetXaxis().SetTitle("%s [GeV]"%distTitle)
#        histAlpha.SetTitle("")
#
#        histBeta  = mlbVar.createHistogram("Morphing against mass variations",betaVar)
#        morphPDF.fillHistogram(histBeta ,ROOT.RooArgList(mlbVar,betaVar))
#        histBeta.GetYaxis().SetTitleOffset(1.1)
#        histBeta.GetYaxis().SetTitle("Grid mass index")
#        histBeta.GetXaxis().SetTitleOffset(1.35)
#        histBeta.GetXaxis().SetTitle("%s [GeV]"%distTitle)
#        histBeta.SetTitle("")
#
#        alphaVar.setVal(1)
#        betaVar.setVal(1)
#
#        for theta in xrange(0,359,int(360/n3DScan)) :
#            # mlb vs. mass plot for 3D scan
#            if theta == 0 : zLimList["mlb"]=histAlpha.GetMaximum()*1.1
#            histAlpha.SetMaximum(zLimList["mlb"])
#            histAlpha.Draw("SURF1")
#            status.DrawLatexNDC(0.75,0.02,"m_{t}=%3.2f GeV"%(betaVar.getVal()*(maxMT-minMT)/betaVar.getMax()+minMT))
#            ROOT.gPad.SetPhi(theta)
#
#            if first:
#                canvas.SetLeftMargin(canvas.GetRightMargin()*1.1)
#
#            canvas.SaveAs("%s/gifplots_%s_%s/rotscanalpha_%s_%i.png"%(opt.output,sig,dist,sig,theta))
#            canvas.Clear()
#
#            # mlb vs. width plot for 3D scan
#            histBeta.SetMaximum(zLimList["mlb"])
#            histBeta.Draw("SURF1")
#            status.DrawLatexNDC(0.75,0.02,"#Gamma_{t}=%3.2f GeV"%(alphaVar.getVal()*(maxGammaT-minGammaT)/alphaVar.getMax()+minGammaT))
#            ROOT.gPad.SetPhi(theta)
#
#            if first:
#                canvas.SetLeftMargin(canvas.GetRightMargin()*1.1)
#
#            canvas.SaveAs("%s/gifplots_%s_%s/rotscanbeta_%s_%i.png"%(opt.output,sig,dist,sig,theta))
#            canvas.Clear()
#
#            # save images to proper arrays
#            imgAlpha=Image.open("%s/gifplots_%s_%s/rotscanalpha_%s_%i.png"%(opt.output,sig,dist,sig,theta))
#            rotScanAlphaImgs+=[imgAlpha]
#            imgBeta=Image.open("%s/gifplots_%s_%s/rotscanbeta_%s_%i.png"%(opt.output,sig,dist,sig,theta))
#            rotScanBetaImgs+=[imgBeta]
#
#        ##########################################
#        # create 2D plots that scan in variables #
#        ##########################################
#        canvas.SetLogx(False)
#        canvas.SetLogy(False)
#        canvas.SetLogz(False)
#
#        if first:
#            canvas.SetRightMargin(canvas.GetRightMargin()*2)
#            canvas.SetLeftMargin(canvas.GetLeftMargin()*2)
#
#        for step in xrange(0,n2DScan) :
#            alphaStep=alphaVar.getMin()+(alphaVar.getMax()-alphaVar.getMin())/n2DScan*step
#            betaStep = betaVar.getMin()+( betaVar.getMax()- betaVar.getMin())/n2DScan*step
#
#            alphaVar.setVal(alphaStep);
#            betaVar.setVal( betaStep);
#
#            # book histos
#            histAlpha = mlbVar.createHistogram("Morphing against mass variations %i"%step,alphaVar)
#            histBeta  = mlbVar.createHistogram("Morphing against width variations %i"%step,betaVar)
#            sensHistAlpha = ROOT.TH2D("sensAlpha %i"%step,"Sensitivity against mass variations",
#                    30,0,300,
#                    n2DScan, alphaVar.getMin(), alphaVar.getMax())
#            sensHistBeta  = ROOT.TH2D("sensBeta %i"%step,"Sensitivity against width variations",
#                    30,0,300,
#                    n2DScan, betaVar.getMin(), betaVar.getMax())
#
#            morphPDF.fillHistogram(histAlpha,ROOT.RooArgList(mlbVar,alphaVar))
#            morphPDF.fillHistogram(histBeta ,ROOT.RooArgList(mlbVar,betaVar))
#
#            # mlb histo against mass
#            histAlpha.SetTitle("")
#            histAlpha.GetXaxis().SetTitle("%s [GeV]"%distTitle)
#            histAlpha.GetYaxis().SetTitle("Grid width index")
#            histAlpha.GetZaxis().SetRangeUser(0,zLimList["mlb"])
#            histAlpha.Draw("COLZ")
#            status.DrawLatexNDC(0.25,0.02,"m_{t}=%3.2f GeV"%(betaVar.getVal()*(maxMT-minMT)/betaVar.getMax()+minMT))
#
#            canvas.SaveAs("%s/gifplots_%s_%s/betascan_%s_%i.png"%(opt.output,sig,dist,sig,step))
#            canvas.Clear()
#
#            # mlb histo against width
#            histBeta.SetTitle("")
#            histBeta.GetYaxis().SetTitle("Grid mass index")
#            histBeta.GetXaxis().SetTitle("%s [GeV]"%distTitle)
#            histBeta.GetZaxis().SetRangeUser(0,zLimList["mlb"])
#            histBeta.Draw("COLZ")
#            status.DrawLatexNDC(0.25,0.02,"#Gamma_{t}=%3.2f GeV"%(alphaVar.getVal()*(maxGammaT-minGammaT)/alphaVar.getMax()+minGammaT))
#
#            canvas.SaveAs("%s/gifplots_%s_%s/alphascan_%s_%i.png"%(opt.output,sig,dist,sig,step))
#            canvas.Clear()
#
#
#            #####################################
#            # create 2D local sensitivity scans #
#            #####################################
#            for ia,ib,ix in [(ia,ib,ix) for ia in range(1,n2DScan+1)
#                    for ib in range(1,n2DScan+1)
#                    for ix in range(1,31)] :
#                if not opt.makeSens : break
#
#                # get all the bin numbers
#                gifBinA=sensHistAlpha.FindBin(ix*3-1.5,alphaStep)
#                gifBinB=sensHistAlpha.FindBin(ix*3-1.5, betaStep)
#                binNumA=sensHistAlpha.GetBin(ix,ia)
#                binNumB= sensHistBeta.GetBin(ix,ib)
#
#                # prepare alpha sensitivity plot
#                alphaVar.setVal(sensHistAlpha.GetXaxis().GetBinCenter(binNumA));
#                betaVar.setVal(  sensHistBeta.GetXaxis().GetBinCenter(gifBinB));
#                mlbVar.setVal(  sensHistAlpha.GetYaxis().GetBinCenter(binNumA));
#                f=pdf.getVal()
#                if f!=0 :
#                    dfdx=pdf.derivative(ws.var('x')).getVal()
#                    sensHistAlpha.Fill(binNumA,(1/f)*(dfdx**2))
#
#                # prepare beta sensitivity plots
#                alphaVar.setVal(sensHistAlpha.GetXaxis().GetBinCenter(gifBinA));
#                betaVar.setVal(  sensHistBeta.GetXaxis().GetBinCenter(binNumB));
#                mlbVar.setVal(   sensHistBeta.GetYaxis().GetBinCenter(binNumB));
#                f=pdf.getVal()
#                if f!=0 :
#                    dfdx=pdf.derivative(ws.var('x')).getVal()
#                    sensHistBeta.Fill(binNumB,(1/f)*(dfdx**2))
#
#            if opt.makeSens :
#                ###################################
#                # sensitivity of mlb against mass #
#                ###################################
#                if step == 0 : zLimList["asens"]=sensHistAlpha.GetMaximum()*2
#                sensHistAlpha.SetTitle("")
#                sensHistAlpha.GetXaxis().SetTitle("%s [GeV]"%distTitle)
#                sensHistAlpha.GetYaxis().SetTitle("Generator-level width [GeV]")
#                sensHistAlpha.GetZaxis().SetRangeUser(0,zLimList["asens"])
#                sensHistAlpha.Draw("COLZ")
#                status.DrawLatexNDC(0.25,0.02,"m_{t}=%3.2f GeV"%(betaStep))
#
#                print sensHistAlpha.GetMaximum()
#
#                canvas.SaveAs("%s/gifplots_%s_%s/sensbetascan_%s_%i.png"%(opt.output,sig,dist,sig,step))
#                canvas.Clear()
#
#                ####################################
#                # sensitivity of mlb against width #
#                ####################################
#                if step == 0 : zLimList["bsens"]=sensHistBeta.GetMaximum()*2
#                sensHistBeta.SetTitle("")
#                sensHistBeta.GetXaxis().SetTitle("%s [GeV]"%distTitle)
#                sensHistBeta.GetYaxis().SetTitle("Generator-level mass [GeV]")
#                sensHistBeta.GetZaxis().SetRangeUser(0,zLimList["bsens"])
#                sensHistBeta.Draw("COLZ")
#                status.DrawLatexNDC(0.25,0.02,"#Gamma_{t}=%3.2f GeV"%(alphaStep))
#
#                canvas.SaveAs("%s/gifplots_%s_%s/sensalphascan_%s_%i.png"%(opt.output,sig,dist,sig,step))
#                canvas.Clear()
#
#                # save images
#                sensImgAlpha=Image.open("%s/gifplots_%s_%s/sensalphascan_%s_%i.png"%(opt.output,sig,dist,sig,step))
#                sensAlphaScanImgs+=[sensImgAlpha]
#                sensImgBeta=Image.open("%s/gifplots_%s_%s/sensbetascan_%s_%i.png"%(opt.output,sig,dist,sig,step))
#                sensBetaScanImgs+=[sensImgBeta]
#
#
#            # save images
#            imgAlpha=Image.open("%s/gifplots_%s_%s/alphascan_%s_%i.png"%(opt.output,sig,dist,sig,step))
#            alphaScanImgs+=[imgAlpha]
#            imgBeta=Image.open("%s/gifplots_%s_%s/betascan_%s_%i.png"%(opt.output,sig,dist,sig,step))
#            betaScanImgs+=[imgBeta]
#
#
#        # make gif files from arrays
#        writeGif("%s/gifplots_%s_%s/alphascan_%s.gif"%(opt.output,sig,dist,sig),alphaScanImgs,duration=.05)
#        writeGif("%s/gifplots_%s_%s/betascan_%s.gif"%( opt.output,sig,dist,sig), betaScanImgs,duration=.05)
#        writeGif("%s/gifplots_%s_%s/rotscanalpha_%s.gif"%(opt.output,sig,dist,sig),rotScanAlphaImgs,duration=.1)
#        writeGif("%s/gifplots_%s_%s/rotscanbeta_%s.gif"%( opt.output,sig,dist,sig), rotScanBetaImgs,duration=.1)
#
#        if opt.makeSens :
#            writeGif("%s/gifplots_%s_%s/sensalphascan_%s.gif"%(opt.output,sig,dist,sig),sensAlphaScanImgs,duration=.05)
#            writeGif("%s/gifplots_%s_%s/sensbetascan_%s.gif"%(opt.output,sig,dist,sig),sensBetaScanImgs,duration=.05)
#
#        first=False
#
#"""
#for execution from another script
#"""
#if __name__ == "__main__":
#    sys.exit(main())
