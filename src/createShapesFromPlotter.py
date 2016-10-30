import os
import sys
import optparse
import ROOT
import commands
import getpass
import pickle
import numpy


ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
#ROOT.gROOT.SetBatch()

def replaceBadCharacters(inputStr):
    newStr = inputStr
    for token in ['+','-','*',' ','#','{','(',')','}','@']:
        newStr=newStr.replace(token,'')
    return newStr

"""
test if variation is significant enough i.e. if sum_{bins} |var-nom| > tolerance
"""
def acceptVariationForDataCard(nomH,upH,downH,tol=1e-2):
    diffUp,diffDown=0,0
    for xbin in xrange(1,nomH.GetNbinsX()):
        val,valUp,valDown=nomH.GetBinContent(xbin),upH.GetBinContent(xbin),downH.GetBinContent(xbin)
        diffUp+=ROOT.TMath.Abs(valUp-val)
        diffDown+=ROOT.TMath.Abs(valDown-val)
    accept = True if (diffUp>tol or diffDown>tol) else False
    return accept

"""
in case of multiple signals, remove others
"""
def filterShapeList(exp,signalList,rawSignalList):
    newExp={}
    for key in exp:

        matchFound=False
        for rawKey in rawSignalList:
            if rawKey in key:
                matchFound=True
        if matchFound and not key in signalList : continue

        newExp[key]=exp[key]
    return newExp


"""
get distributions from file
"""
def getDistsFrom(directory,addList=None,addTitle='',keyFilter=''):
    obs=None
    exp={}
    dirName=directory.GetName()
    addedHistProduced=False
    for key in directory.GetListOfKeys():
        if len(keyFilter)>0 and key.GetName()!='%s_%s'%(dirName,keyFilter) : continue
        obj=directory.Get(key.GetName())
        if not obj.InheritsFrom('TH1') : continue
        inAddList = False;
        if obj.GetName()==dirName :
            obs=obj.Clone('data_obs')
            obs.SetDirectory(0)
        else :
            newName=obj.GetName().split(dirName+'_')[-1]
            for token in ['+','-','*',' ','#','{','(',')','}','@']:
                newName=newName.replace(token,'')
            if addList and newName in addList :
                newName='added' if addTitle == '' else addTitle
                if addedHistProduced :
                    exp[newName].Add(obj)
                    addedHistProduced = True;
                    inAddList=True
                else :
                    exp[newName] = obj.Clone(newName)
            else :
                exp[newName]=obj.Clone(newName)
                exp[newName].SetDirectory(0)

            #newName = ('added' if addTitle == '' else addTitle) if inAddList and addedHistProduced else newName
            if obj.InheritsFrom('TH2'):
                for xbin in xrange(1,exp[newName].GetNbinsX()+1):
                    for ybin in xrange(1,exp[newName].GetNbinsY()+1):
                        binContent=exp[newName].GetBinContent(xbin,ybin)
                        if binContent>0: continue
                        newBinContent=ROOT.TMath.Max(ROOT.TMath.Abs(binContent),1e-3)
                        exp[newName].SetBinContent(xbin,ybin,newBinContent)
                        exp[newName].SetBinError(xbin,ybin,newBinContent)
            else:
                for xbin in xrange(1,exp[newName].GetNbinsX()+1):
                    binContent=exp[newName].GetBinContent(xbin)
                    if binContent>0: continue
                    newBinContent=ROOT.TMath.Max(ROOT.TMath.Abs(binContent),1e-3)
                    exp[newName].SetBinContent(xbin,newBinContent)
                    exp[newName].SetBinError(xbin,newBinContent)
    return obs,exp

"""
merge many directories with the same directory base-name
  '.' --> 'p'
"""
def getMergedDists(fIn,basedir='',addList=None,addTitle='',
                   widList=None,nomWid='',sigList=None,keyFilter=''):
    obs=None
    exp={}
    fIn.Get._creates=True
    for wid in widList :
        dirName='%s%s'%(basedir,wid)
        #print dirName
        directory = fIn.Get(dirName)
        if wid == nomWid :
            obs,mergeExp=getDistsFrom(directory=directory,addList=addList,addTitle=addTitle,keyFilter=keyFilter)
            for key in mergeExp :
                newName=key
                if key in sigList :
                    newName=('%s%s'%(key,wid.replace('.','p')))
                exp[newName]=mergeExp[key].Clone(newName)
        else :
            _,mergeExp=getDistsFrom(directory=directory,addList=addList,addTitle=addTitle,keyFilter=keyFilter)
            for key in mergeExp :
                if key in sigList :
                    newName=('%s%s'%(key,wid.replace('.','p')))
                    exp[newName]=mergeExp[key].Clone(newName)
                else : continue
    return obs,exp


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
make an MC truth dataset (sigs+bkgs)
"""
def makeMCTruthHist(hypothesis,sigList,dists,infile=None,dirhandle="",filtername="",exthypo=""):
    outputHist=None
    firstLoop=True
    print "Producing MC Truth data for hypothesis %s"%hypothesis

    tdists=dists.copy()

    #replace with whatever we find in the provided dir
    if dirhandle != "" and (infile is not None) and filtername != "" and exthypo!="":
        _,pdists=getMergedDists(infile,dirhandle,None,"",[exthypo.replace('p','.')],"1.0w",sigList)
        for key in pdists :
            if filtername not in key : continue
            tdists[key.replace(filtername,'')+hypothesis] = pdists[key].Clone()
            tdists[key.replace(filtername,'')+hypothesis].SetDirectory(0)

    for sig in sigList :
        sigHist=tdists["%s%s"%(sig,hypothesis)].Clone()
        if firstLoop :
            outputHist=sigHist.Clone()
            firstLoop=False
        else :
            outputHist.Add(sigHist)
    for dist in tdists :
        isSig=False
        for sig in sigList :
            if sig in dist and 'V' not in dist: isSig=True
        if isSig : continue
        bkgHist=tdists[dist].Clone()
        outputHist.Add(bkgHist)
        outputHist.SetName("data_obs")

    outputHist.SetDirectory(0)

    return outputHist

"""
steer the script
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--input',     dest='input',     help='input plotter',                            default=None,            type='string')
    parser.add_option(      '--systInput', dest='systInput', help='input plotter for systs from alt samples', default=None,            type='string')
    parser.add_option('-d', '--dists',     dest='distList',  help='distribution',        default='minmlb,mdrmlb,incmlb,sncmlb,mt2mlb', type='string')
    parser.add_option('-s', '--signal',    dest='signal',    help='signal (csv)',                             default='tbart,tW',      type='string')
    parser.add_option('-c', '--cat',       dest='cat',       help='categories (csv)',                         default='1b,2b',         type='string')
    parser.add_option('-o', '--output',    dest='output',    help='output directory',                         default='datacards',     type='string')
    parser.add_option('-n', '--outname',   dest='outname',   help='output file name',                         default='shapes',     type='string')
    parser.add_option(      '--addSigs',   dest='addSigs',   help='signal processes to add',                  default=False,     action='store_true')
    parser.add_option(      '--lfs',       dest='lfsInput',  help='lepton final states to consider',          default='EE,EM,MM',      type='string')
    parser.add_option(      '--lbCat',     dest='lbCat',     help='pt categories to consider',                default='highpt,lowpt',  type='string')
    parser.add_option(      '--truth', dest='truthDataset',  help='make data out of MC truth',                default='',              type='string')
    parser.add_option(      '--trx4', dest='truthExtDataset',help='make data out of 4xMC truth',              default='',              type='string')
    parser.add_option(      '--trx4Width', dest='truthExtWidth',help='4xMC truth has this width',             default='',              type='string')
    parser.add_option(      '--trx4Input', dest='truthExtfIn'  ,help='4xMC truth has this infile',            default='',              type='string')
    parser.add_option(      '--min', dest='minIndex',help='start at this card',              default=-5,              type=int)
    parser.add_option(      '--max', dest='maxIndex',help='end at this card',                default=10000,           type=int)
    parser.add_option('--noshapes', dest='skipMakingShapes', help='jump straight to morphing',                     default=False, action='store_true')
    parser.add_option('--nomorph',  dest='skipMorphing',     help='do not morph signal dists',                     default=False, action='store_true')
    parser.add_option('--allmorph', dest='allMorphs',        help='make morph validation plots for all cats',      default=False, action='store_true')
    parser.add_option('--makesens', dest='makeSens',         help='make local sensitivity validation (very slow)', default=False, action='store_true')
    parser.add_option('--novalidation', dest='noValidation', help='do not make any validation plots',              default=False, action='store_true')
    parser.add_option('-w', '--wids',  dest='widList',       help='signal widths',
            default='0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.5,4.0',
            type='string')
    (opt, args) = parser.parse_args()

    # parse the dists to consider
    distList = opt.distList.split(',')

    # parse the channels, pt(l,b), b categories to consider
    lfsList=opt.lfsInput.split(',')
    lbCatList=opt.lbCat.split(',')
    catList=opt.cat.split(',')

    # what are our signal processes?
    rawSignalList=opt.signal.split(',')
    signalList = { 'top' } if opt.addSigs else rawSignalList[:]

    # what are our widths?
    rawWidList = opt.widList.split(',')
    widList = rawWidList[:]
    for i in xrange(0,len(rawWidList)) :
        widList[i] = ('%sw'%(rawWidList[i]))
    nomWid="1.0w"
    modNomWid="1p0w"

    modWidList=widList[:]
    for i in xrange(0,len(widList)) :
        modWidList[i]=widList[i].replace('.','p')

    # HARDCODE systematics
    tWRateSystList=['tW']
    ttRateSystList=['tbart']
    tWIsSig=('tW'    in rawSignalList)
    ttIsSig=('tbart' in rawSignalList)
    for wid in modWidList :
        if tWIsSig : tWRateSystList += ["tW%s"%wid]
        if ttIsSig : ttRateSystList += ["tbart%s"%wid]

    rateSysts=[
          ('lumi_13TeV',       1.062,    'lnN',    []                   ,['Multijetsdata']),
          ('DYnorm_th',        1.30,     'lnN',    ['DY']  ,[]),
          ('Wnorm_th',         1.30,     'lnN',    ['W']   ,[]),
          ('tWnorm_th',        1.054,    'lnN',    tWRateSystList,[]),
          ('tnorm_th',         1.044,    'lnN',    ['tch']              ,[]),
          ('VVnorm_th',        1.20,     'lnN',    ['Multiboson']       ,[]),
          ('tbartVnorm_th',    1.30,     'lnN',    ['tbartV']           ,[]),
          ('sel',              1.02,     'lnN',    tWRateSystList+ttRateSystList+['DY','W','tch','Multiboson','tbartV'], []),
    ]

    Mtop       =['tbartm=169.5','tbartm=175.5','tWm=169.5','tWm=175.5']
    ttParton_tt=['tbartscaledown','tbartscaleup']
    ttParton_tW=['tWscaledown','tWscaleup']
    #ME={    "muF": ['%sgen%i'%(sig,ind) for sig in ['tbart'] for ind in [3,4]],
    #        "muR": ['%sgen%i'%(sig,ind) for sig in ['tbart'] for ind in [5,8]],
    #        "tot": ['%sgen%i'%(sig,ind) for sig in ['tbart'] for ind in [6,10]]
    #        }
    tWinterf=['tWDS']
    amcNLO  =['tbartamcnloFxFx']
    Herwig  =['tbartHerwig']

    systSignalList=Mtop+ttParton_tt+amcNLO+Herwig+tWinterf+ttParton_tW

    MtopMap={}
    ttPartonMap={}
    twPartonMap={}
    tWinterfMap={}
    amcNLOMap={}
    HerwigMap={}

    for wid in modWidList :
        MtopMap[    'tbart%s'%wid]=['tbartm=169.5%s'  %(wid),'tbartm=175.5%s'%(wid)]
        MtopMap[    'tW%s'   %wid]=['tWm=169.5%s'     %(wid),'tWm=175.5%s'   %(wid)]

        ttPartonMap['tbart%s'%wid]=['tbartscaledown%s'%(wid),'tbartscaleup%s'%(wid)]
        twPartonMap['tW%s'   %wid]=['tWscaledown%s'   %(wid),'tWscaleup%s'   %(wid)]

        amcNLOMap['tbart%s'%wid]=['tbartamcnloFxFx%s'%(wid)]

        HerwigMap['tbart%s'%wid]=['tbartHerwig%s'%(wid)]

        tWinterfMap['tW%s'   %wid]=['tWDS%s'   %(wid)]


    sampleSysts=[
          #ttbar modelling
          ('Mtop'          , MtopMap    , True , True, False),
          ('ttPartonShower', ttPartonMap, False, True, False),
          ('tWQCDScale'    , twPartonMap, False, True, False),
          ('Herwig'        , HerwigMap  , False, True, True),
          ('amcnloFxFx'    , amcNLOMap  , False, True, True),
          #tWinterference
          ('tWttinterf'    , tWinterfMap, False, True, True),
    ]


    genSysts=[
        ('jes',   True,False,False),
        ('les',   False,False,False),
        ('ltag',  False,False,False),
        ('trig',  False,False,False),
        #('sel',   False,False,False),
        ('toppt',  True,False,False),
        ('jer',   False,False,False),
        ('btag',  False,False,False),
        ('pu',     True,False,False),
        ('MEqcdscale',  True,False,False),
        ('PDF',   False,False,False)
    ]

    # prepare output directory
    os.system('mkdir -p %s'%opt.output)

    anCat=''
    for subDir in opt.input.split('/'):
        if 'analysis_' not in subDir: continue
        anCat=subDir.replace('analysis_','')

    # get data and nominal expectations
    fIn=ROOT.TFile.Open(opt.input)
    fIn.Get._creates=True
    systfIn=None
    if opt.systInput:
        systfIn=ROOT.TFile.Open(opt.systInput)
        systfIn.Get._creates=True

    # prepare output ROOT file
    if not opt.skipMakingShapes :
        outFile='%s/%s.root'%(opt.output,opt.outname)
        fOut=ROOT.TFile.Open(outFile,'UPDATE')
        fOut.Close()

    # keep track of progress
    cardIndex=-1
    numCards =(len(rawWidList))*len(lfsList)*len(lbCatList)*len(catList)*len(distList)

    #loop over lepton final states
    for (lbCat,lfs,cat,dist) in [(a,b,c,d) for a in lbCatList
            for b in lfsList
            for c in catList
            for d in distList] :
        if opt.skipMakingShapes : break

        if cardIndex >= opt.maxIndex : break
        if cardIndex+len(rawWidList) < opt.minIndex:
            cardIndex+=len(rawWidList)
            continue

        obs=ROOT.TH1F('','',100,0,100)
        exp={}
        #nominal expectations
        if opt.addSigs :
            obs,exp=getMergedDists(fIn,('%s%s%s_%s_'%(lbCat,lfs,cat,dist)),rawSignalList,'top',widList,nomWid,signalList)
        else :
            obs,exp=getMergedDists(fIn,('%s%s%s_%s_'%(lbCat,lfs,cat,dist)),None,'',widList,nomWid,signalList)

        if not opt.truthDataset=="" :
            obs=makeMCTruthHist(opt.truthDataset,signalList,exp,
                    (ROOT.TFile(opt.truthExtfIn) if opt.truthExtfIn != "" else None),
                    ('%s%s%s_%s_'%(lbCat,lfs,cat,dist)),
                    opt.truthExtDataset,opt.truthExtWidth)
        #exp=filterShapeList(exp,signalList,rawSignalList)

        nomShapes=exp.copy()
        nomShapes['data_obs']=obs

        for tSig in signalList :
            nomShapes['%sNOM'%(tSig)] = exp['%s%s'%(tSig,modNomWid)].Clone()

        saveToShapesFile(outFile,nomShapes,('%s%s%s_%s'%(lbCat,lfs,cat,dist)))

        #loop over categories, widths
        for wid in modWidList:

            cardIndex+=1
            if cardIndex < opt.minIndex :
                print 'Skipping %s datacards for %s%s%s_%s \t\t [%i/%i]'%(dist,lbCat,lfs,cat,wid,cardIndex,numCards)
                continue
            if cardIndex > opt.maxIndex :
                break
            print 'Initiating %s datacard for %s%s%s_%s \t\t [%i/%i]'%(dist,lbCat,lfs,cat,wid,cardIndex,numCards)


            #start the datacard
            datacard=open('%s/datacard__%s_%s%s%s_%s.dat'%(opt.output,wid,lbCat,lfs,cat,dist),'w')
            datacard.write('#\n')
            datacard.write('# Generated by %s with git hash %s for analysis category %s%s%s_%s\n' % (getpass.getuser(),
                commands.getstatusoutput('git log --pretty=format:\'%h\' -n 1')[1],
                lbCat, lfs, cat, wid) )

            datacard.write('#\n')
            datacard.write('imax *\n')
            datacard.write('jmax *\n')
            datacard.write('kmax *\n')
            datacard.write('-'*50+'\n')
            datacard.write('shapes *        * shapes.root %s%s%s_%s/$PROCESS %s%s%s_%s_$SYSTEMATIC/$PROCESS\n'%(lbCat,lfs,cat,dist,lbCat,lfs,cat,dist))
            datacard.write('-'*50+'\n')
            datacard.write('bin 1\n')
            datacard.write('observation %3.1f\n' % obs.Integral())
            datacard.write('-'*50+'\n')

            #expectations
            datacard.write('\t\t\t %32s'%'bin')
            for i in xrange(0,len(exp)+(2 if opt.addSigs else 0)-(len(modWidList)-1)*len(signalList)+len(signalList)):
                datacard.write('%15s'%'1')
            datacard.write('\n')
            datacard.write('\t\t\t %32s'%'process')
            for sig in signalList: datacard.write('%15s'%('%s%s'%(sig,modNomWid)))
            if modNomWid == wid :
                for sig in signalList: datacard.write('%15s'%('%sNOM'%(sig)))
            else :
                for sig in signalList: datacard.write('%15s'%('%s%s'%(sig,wid)))
            for proc in exp:
                isSig=False
                for sig in modWidList :
                    if sig in proc and 'V' not in proc: isSig=True
                if isSig : continue
                datacard.write('%15s'%proc)
            datacard.write('\n')
            datacard.write('\t\t\t %32s'%'process')
            for i in xrange(0,2*len(signalList)) : datacard.write('%15s'%str(i+1-2*len(signalList)))
            i=0
            for proc in exp:
                isSig=False
                for sig in modWidList :
                    if sig in proc and 'V' not in proc: isSig=True
                if isSig : continue
                i=i+1
                datacard.write('%15s'%str(i))
            datacard.write('\n')
            datacard.write('\t\t\t %32s'%'rate')
            for sig in signalList :
                datacard.write('%15s'%('%3.2f'%(exp['%s%s'%(sig,modNomWid)].Integral())))
            for sig in signalList :
                datacard.write('%15s'%('%3.2f'%(exp['%s%s'%(sig,wid)].Integral())))
            for proc in exp:
                isSig=False
                for sig in modWidList :
                    if sig in proc and 'V' not in proc: isSig=True
                if isSig : continue
                datacard.write('%15s'%('%3.2f'%exp[proc].Integral()))
            datacard.write('\n')
            datacard.write('-'*50+'\n')


            #rate systematics
            try:
                jetCat=cat[:-2] if cat.endswith('t') else cat
                rateSysts.append( ('MultiJetsNorm%s%s'%(jetCat,anCat),  1+qcdNorm[jetCat][1],                       'lnN',    ['Multijetsdata']    ,[]) )
            except:
                pass

            for syst,val,pdf,whiteList,blackList in rateSysts:
                if syst in ['sel','DYnorm_th'] : syst="%s%s"%(syst,lfs)


                datacard.write('%32s %8s'%(syst,pdf))
                entryTxt=''
                try:
                    entryTxt='%15s'%('%3.3f/%3.3f'%(ROOT.TMath.Max(val[0],0.01),val[1]))
                except:
                    entryTxt='%15s'%('%3.3f'%val)

                for sig in signalList:
                    newSig='%s%s'%(sig,modNomWid)
                    if (len(whiteList)==0 and not newSig in blackList) or newSig in whiteList:
                        datacard.write(entryTxt)
                    else:
                        datacard.write('%15s'%'-')
                for sig in signalList:
                    newSig='%s%s'%(sig,wid)
                    if (len(whiteList)==0 and not newSig in blackList) or newSig in whiteList:
                        datacard.write(entryTxt)
                    else:
                        datacard.write('%15s'%'-')
                for proc in exp:
                    isSig=False
                    for sig in modWidList :
                        if sig in proc and 'V' not in proc: isSig=True
                    if isSig : continue
                    if (len(whiteList)==0 and not proc in blackList) or proc in whiteList:
                        datacard.write(entryTxt)
                    else:
                        datacard.write('%15s'%'-')
                datacard.write('\n')

            #generator level systematics
            if systfIn is None :
                datacard.close()
                continue


            # sample systematics
            _,genVarShapes = getMergedDists(fIn,
                    ('%s%s%s_%s_'%(lbCat,lfs,cat,dist)),
                    (rawSignalList if opt.addSigs else None),
                    ('top' if opt.addSigs else ''),
                    widList,nomWid,systSignalList)
            _,altExp       = getMergedDists(systfIn,
                    ('%s%s%s_%s_'%(lbCat,lfs,cat,dist)),
                    (rawSignalList if opt.addSigs else None),
                    ('top' if opt.addSigs else ''),
                    widList,nomWid,systSignalList)

            for systVar, procsToApply, normalize, useAltShape, projectRelToNom in sampleSysts:

                #prepare shapes and check if variation is significant
                downShapes, upShapes = {}, {}

                for iproc in procsToApply:
                    nomH=exp[iproc]

                    #check which shape to use
                    if useAltShape:

                        #get directly from another file
                        downH  = altExp[ procsToApply[iproc][0] ]
                        if len( procsToApply[iproc] ) > 1 :
                            upH    = altExp[ procsToApply[iproc][1] ]
                        else:
                            #if only one variation is available, mirror it
                            upH = downH.Clone( '%s%sUp'%(iproc,systVar) )
                            for xbin in xrange(1,upH.GetNbinsX()+1):
                                diff=upH.GetBinContent(xbin)-nomH.GetBinContent(xbin)
                                if nomH.GetBinContent(xbin)-diff >= 0 :
                                    upH.SetBinContent(xbin,nomH.GetBinContent(xbin)-diff)
                                else :
                                    upH.SetBinContent(xbin,upH.GetBinContent(xbin-1))
                    else:

                        #project from 2D histo (re-weighted from nominal sample)
                        ybinUp, ybinDown = -1, -1
                        for ybin in xrange(1,genVarShapes[ iproc ].GetNbinsY()+1):
                            label = genVarShapes[ iproc ].GetYaxis().GetBinLabel(ybin)
                            if procsToApply[iproc][0] in label : ybinDown=ybin
                            if procsToApply[iproc][1] in label : ybinUp=ybin

                        downH = genVarShapes[ iproc ].ProjectionX('%s%sDown'%(iproc,systVar), ybinDown, ybinDown)
                        upH   = genVarShapes[ iproc ].ProjectionX('%s%sUp'%(iproc,systVar),   ybinUp,   ybinUp)

                    # use do down/up x nom to generate the variation, then mirror it
                    if projectRelToNom:
                        ratioH=downH.Clone()
                        ratioH.Divide(upH)
                        for xbin in xrange(1,nomH.GetNbinsX()+1):
                            nomVal=nomH.GetBinContent(xbin)
                            varVal = ratioH.GetBinContent(xbin) * nomVal
                            upH.SetBinContent(xbin, varVal)
                            if upH.GetBinContent(xbin) <= 0 :
                                upH.SetBinContent(xbin,1e-010)
                            varVal = varVal- nomVal
                            downH.SetBinContent(xbin, nomVal-varVal)
                            if downH.GetBinContent(xbin) <= 0 :
                                downH.SetBinContent(xbin,1e-010)

                    #normalize (shape only variation is considered)
                    if normalize : downH.Scale( nomH.Integral()/downH.Integral() )
                    if normalize : upH.Scale( nomH.Integral()/upH.Integral() )

                    #check if variation is meaningful
                    accept = acceptVariationForDataCard(nomH=nomH, upH=upH, downH=downH)
                    if not accept : continue

                    #save
                    downShapes[iproc]=downH
                    upShapes[iproc]=upH

                    if modNomWid in iproc:
                        for tSig in signalList :
                            downShapes['%sNOM'%tSig]=downH
                            upShapes['%sNOM'%tSig]=upH

                #check if something has been accepted
                if len(upShapes)==0 : continue

                #export to shapes file
                saveToShapesFile(outFile,downShapes,lbCat+lfs+cat+"_"+dist+"_"+systVar+'Down')
                saveToShapesFile(outFile,upShapes,lbCat+lfs+cat+"_"+dist+"_"+systVar+'Up')

                #write to datacard
                datacard.write('%26s shape'%systVar)
                for sig in signalList:
                    if ("%s%s"%(sig,wid)) in procsToApply and ("%s%s"%(sig,wid)) in upShapes:
                        if "Mtop" in systVar:
                            datacard.write('%15s'%'0.167')
                        else:
                            datacard.write('%15s'%'1')
                    else:
                        datacard.write('%15s'%'-')
                for sig in signalList:
                    if ("%s%s"%(sig,wid)) in procsToApply and ("%s%s"%(sig,wid)) in upShapes:
                        if "Mtop" in systVar:
                            datacard.write('%15s'%'0.167')
                        else:
                            datacard.write('%15s'%'1')
                    else:
                        datacard.write('%15s'%'-')
                for proc in exp:
                    isSig=False
                    for sig in modWidList :
                        if sig in proc and 'V' not in proc: isSig=True
                    if isSig : continue
                    if proc in procsToApply and proc in upShapes:
                        if "Mtop" in systVar:
                            datacard.write('%15s'%'0.167')
                        else:
                            datacard.write('%15s'%'1')
                    else:
                        datacard.write('%15s'%'-')
                datacard.write('\n')


            #
            for systVar, normalize, useAltShape, projectRelToNom in genSysts:

                genVarShapesUp,genVarShapesDn=None,None

                if 'PDF'==systVar:

                    #fill in 2D histos y-axis is the replica number
                    genVarShapesDn2D={}
                    nPdfReplicas=100
                    startIdx=11
                    for igen in xrange(startIdx,nPdfReplicas+startIdx+1):
                        dirForSyst=('gen%d%s%s%s_%s_'%(igen,lbCat,lfs,cat,dist))
                        _,iGenVarShapesDn = getMergedDists((systfIn if useAltShape else fIn),
                                                           dirForSyst,
                                                           (rawSignalList if opt.addSigs else None),
                                                           ('top' if opt.addSigs else ''),
                                                           widList,nomWid,signalList)
                        for proc in iGenVarShapesDn:
                            nbins=iGenVarShapesDn[proc].GetXaxis().GetNbins()
                            if not proc in genVarShapesDn2D:
                                name=iGenVarShapesDn[proc].GetName()+'2D'
                                title=iGenVarShapesDn[proc].GetTitle()
                                xmin=iGenVarShapesDn[proc].GetXaxis().GetXmin()
                                xmax=iGenVarShapesDn[proc].GetXaxis().GetXmax()
                                genVarShapesDn2D[proc]=ROOT.TH2D(name,title,nbins,xmin,xmax,nPdfReplicas,0,nPdfReplicas)
                                genVarShapesDn2D[proc].SetDirectory(0)
                            for xbin in xrange(1,nbins+1):
                                genVarShapesDn2D[proc].SetBinContent(xbin,igen-startIdx+1,iGenVarShapesDn[proc].GetBinContent(xbin))

                    #fill the down shapes with the bin contents being the mean-RMS of the replicas
                    genVarShapesDn={}
                    for proc in genVarShapesDn2D:
                        #down shape
                        name=genVarShapesDn2D[proc].GetName()[:-2]+"PDF"
                        genVarShapesDn[proc]=genVarShapesDn2D[proc].ProjectionX(name,1,1)
                        genVarShapesDn[proc].SetDirectory(0)
                        genVarShapesDn[proc].Reset('ICE')

                        #project to compute the mean and RMS
                        for xbin in xrange(1,genVarShapesDn2D[proc].GetXaxis().GetNbins()+1):
                            tmp=genVarShapesDn2D[proc].ProjectionY("tmp",xbin,xbin)

                            tmean=numpy.zeros(tmp.GetXaxis().GetNbins())
                            for txbin in xrange(1,tmp.GetXaxis().GetNbins()+1) :
                                tmean[txbin-1]=tmp.GetBinContent(txbin)

                            mean=numpy.mean(tmean)
                            rms=numpy.std(tmean)

                            tmp.Delete()
                            genVarShapesDn[proc].SetBinContent(xbin,mean-rms)

                        #all done, can remove this histo from memory
                        genVarShapesDn2D[proc].Delete()

                else:

                    genVarShapesUp,genVarShapesDn={},{}
                    #construct an envelope from the 6 QCD scale variations, then mirror it
                    if 'MEqcdscale' in systVar:
                        for igen in [3,5,6,4,8,10]:
                            dirForSyst=('gen%d%s%s%s_%s_'%(igen,lbCat,lfs,cat,dist))
                            _,igenVarShapes = getMergedDists((systfIn if useAltShape else fIn),
                                                             dirForSyst,
                                                             (rawSignalList if opt.addSigs else None),
                                                             ('top' if opt.addSigs else ''),
                                                             widList,nomWid,signalList)
                            for proc in igenVarShapes:
                                if proc not in genVarShapesDn:
                                    name=igenVarShapes[proc].GetName()+'envdn'
                                    genVarShapesDn[proc]=exp[proc].Clone(name)
                                    genVarShapesDn[proc].SetDirectory(0)
                                    name=igenVarShapes[proc].GetName()+'envup'
                                    genVarShapesUp[proc]=exp[proc].Clone(name)
                                    genVarShapesUp[proc].SetDirectory(0)
                                for xbin in xrange(1,igenVarShapes[proc].GetNbinsX()+1):
                                    yvar = igenVarShapes[proc].GetBinContent(xbin)
                                    ydn  = genVarShapesDn[proc].GetBinContent(xbin)
                                    yup  = genVarShapesUp[proc].GetBinContent(xbin)
                                    if yvar>yup:
                                        genVarShapesUp[proc].SetBinContent(xbin,yvar)
                                    elif yvar<ydn:
                                        genVarShapesDn[proc].SetBinContent(xbin,yvar)

                    #proceed as usual
                    else:
                        dirForSyst=('%sup%s%s%s_%s_'%(systVar,lbCat,lfs,cat,dist))
                        _,genVarShapesUp = getMergedDists((systfIn if useAltShape else fIn),
                                                        dirForSyst,
                                                        (rawSignalList if opt.addSigs else None),
                                                        ('top' if opt.addSigs else ''),
                                                        widList,nomWid,signalList)

                        dirForSyst=('%sup%s%s%s_%s_'%(systVar,lbCat,lfs,cat,dist))
                        _,genVarShapesDn = getMergedDists((systfIn if useAltShape else fIn),
                                                         dirForSyst,
                                                        (rawSignalList if opt.addSigs else None),
                                                        ('top' if opt.addSigs else ''),
                                                        widList,nomWid,signalList)

                #prepare shapes and check if variation is significant
                downShapes, upShapes = {}, {}

                for iproc in exp:
                    nomH=exp[iproc]

                    if not iproc in genVarShapesDn :
                        continue
                    #get directly from another file
                    downH  = genVarShapesDn[iproc]
                    if genVarShapesUp is not None and iproc in genVarShapesUp :
                        upH = genVarShapesUp[iproc]
                    else:
                        #if only one variation is available, mirror it
                        upH = downH.Clone( '%s%sUp'%(iproc,systVar) )
                        for xbin in xrange(1,upH.GetNbinsX()+1):
                            diff=upH.GetBinContent(xbin)-nomH.GetBinContent(xbin)
                            upH.SetBinContent(xbin,nomH.GetBinContent(xbin)-diff)

                    # use do down/up x nom to generate the variation, then mirror it
                    if projectRelToNom:
                        ratioH=downH.Clone()
                        ratioH.Divide(upH)
                        for xbin in xrange(1,nomH.GetNbinsX()+1):
                            nomVal=nomH.GetBinContent(xbin)
                            varVal = ratioH.GetBinContent(xbin) * nomVal
                            upH.SetBinContent(xbin, varVal)
                            varVal = varVal- nomVal
                            downH.SetBinContent(xbin, nomVal-varVal)

                    #normalize (shape only variation is considered)
                    if normalize : downH.Scale( nomH.Integral()/downH.Integral() )
                    if normalize : upH.Scale( nomH.Integral()/upH.Integral() )

                    #check if variation is meaningful
                    accept = acceptVariationForDataCard(nomH=nomH, upH=upH, downH=downH)
                    if not accept : continue

                    #save
                    downShapes[iproc]=downH
                    upShapes[iproc]=upH

                    if modNomWid in iproc:
                        for tSig in signalList :
                            downShapes['%sNOM'%tSig]=downH
                            upShapes['%sNOM'%tSig]=upH

                #check if something has been accepted
                if len(upShapes)==0 : continue

                #export to shapes file
                if systVar == "trig" :
                    systVar=systVar+lfs
                extraRateSyst={}
                if systVar in ['jes']:
                    for proc in downShapes:
                        nevtsup,nevtsdn=upShapes[proc].Integral(),downShapes[proc].Integral()
                        if nevtsdn==0 : continue
                        rateVarInProc=1.0+0.5*ROOT.TMath.Abs(1.0-nevtsup/nevtsdn)
                        extraRateSyst[proc]=rateVarInProc
                saveToShapesFile(outFile,downShapes,lbCat+lfs+cat+"_"+dist+"_"+systVar+'Down')
                saveToShapesFile(outFile,upShapes,lbCat+lfs+cat+"_"+dist+"_"+systVar+'Up')

                #write to datacard
                datacard.write('%26s shape'%systVar)
                for sig in signalList:
                    if ("%s%s"%(sig,modNomWid)) in exp and ("%s%s"%(sig,modNomWid)) in upShapes:
                        datacard.write('%15s'%'1')
                    else:
                        datacard.write('%15s'%'-')
                for sig in signalList:
                    if ("%s%s"%(sig,wid)) in exp and ("%s%s"%(sig,wid)) in upShapes:
                        datacard.write('%15s'%'1')
                    else:
                        datacard.write('%15s'%'-')
                for proc in exp:
                    isSig=False
                    for sig in modWidList :
                        if sig in proc and 'V' not in proc: isSig=True
                    if isSig : continue
                    if proc in upShapes:
                        datacard.write('%15s'%'1')
                    else:
                        datacard.write('%15s'%'-')
                datacard.write('\n')

                #add extra rate systs
                if len(extraRateSyst)>0:
                    datacard.write('%26s lnN'%(systVar+'rate'))
                    for sig in signalList:
                        if ("%s%s"%(sig,modNomWid)) in exp and ("%s%s"%(sig,modNomWid)) in extraRateSyst:
                            datacard.write('%15s'%('%3.3f'%extraRateSyst[sig+modNomWid]))
                        else:
                            datacard.write('%15s'%'-')
                    for sig in signalList:
                        if ("%s%s"%(sig,wid)) in exp and ("%s%s"%(sig,wid)) in extraRateSyst:
                            datacard.write('%15s'%('%3.3f'%extraRateSyst[sig+wid]))
                        else:
                            datacard.write('%15s'%'-')
                    for proc in exp:
                        isSig=False
                        for sig in modWidList :
                            if sig in proc and 'V' not in proc: isSig=True
                        if isSig : continue
                        if proc in extraRateSyst:
                            datacard.write('%15s'%('%3.3f'%extraRateSyst[proc]))
                        else:
                            datacard.write('%15s'%'-')
                    datacard.write('\n')


            #all done
            datacard.close()

    ##################
    # BEGIN MORPHING #
    ##################

    if systfIn is None : return
    if opt.skipMorphing : return

    import CMS_lumi
    import tdrStyle
    tdrStyle.setTDRStyle()

    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

    print "\n Creating morphed dists\n"

    kOrangeList=[-4,-2,0,1,-3,7,8,10,9,2]
    widths=map(float,[wid.replace('p','.').replace('w','') for wid in rawWidList])
    massInfo=[(169.5,opt.systInput,'t#bar{t} m=169.5'),
              (172.5,opt.input    ,'t#bar{t}'),
              (175.5,opt.systInput,'t#bar{t} m=175.5'),
              (169.5,opt.systInput,'tW m=169.5'),
              (172.5,opt.input    ,'tW'),
              (175.5,opt.systInput,'tW m=175.5')]
    minMT=169.5
    nomMT=172.5
    maxMT=175.5
    minGammaT=0.5*1.324
    nomGammaT=1.324
    maxGammaT=5.0*1.324
    distToTitle={
            "mlb":    "Inclusive M_{lb}",
            "minmlb": "Minimum M_{lb}",
            "mdrmlb": "#DeltaR-Filtered M_{lb}",
            "incmlb": "Inclusive M_{lb}",
            "sncmlb": "Semi-Inclusive M_{lb}",
            "mt2mlb": "M_{T2}^{lb} Strategy"
           }

    canvas=ROOT.TCanvas()
    procHasUnc={}

    # make workspace for each signal process
    first=True
    for sig,dist,lbCat,ch,cat in [(a,b,c,d,e)
            for a in (['tbart'] if not opt.allMorphs else signalList)
            for b in (['incmlb'] if not opt.allMorphs else distList)
            for c in (['highpt'] if not opt.allMorphs else lbCatList)
            for d in (['EM'] if not opt.allMorphs else lfsList)
            for e in (['2b'] if not opt.allMorphs else catList)]:
        modSig=replaceBadCharacters(sig)
        distTitle="M_{lb}" if dist not in distToTitle else distToTitle[dist]

        if sig == signalList[0] :
            procHasUnc["%s%s%s_%s"%(lbCat,ch,cat,dist)]={}

        # get the proper mass information for this signal
        masses=[]
        for mass,fname,distname in massInfo :
            if modSig in replaceBadCharacters(distname) :
                masses += [(mass,fname,distname)]

        # create workspace
        ws=ROOT.RooWorkspace('sigws_%s%s%s_%s_%s'%(lbCat,ch,cat,modSig,dist))
        ws.factory('x[0,300]')

        ##########################################
        # produce systematic uncertainties in ws #
        ##########################################
        downShapes, upShapes = {}, {}

        #nominal expectations
        exp={}
        if opt.addSigs :
            _,exp=getMergedDists(fIn,('%s%s%s_%s_'%(lbCat,ch,cat,dist)),rawSignalList,'top',widList,nomWid,signalList)
        else :
            _,exp=getMergedDists(fIn,('%s%s%s_%s_'%(lbCat,ch,cat,dist)),None,'',widList,nomWid,signalList)

        # sample systematics
        _,genVarShapes = getMergedDists(fIn,
                ('%s%s%s_%s_'%(lbCat,ch,cat,dist)),
                (rawSignalList if opt.addSigs else None),
                ('top' if opt.addSigs else ''),
                widList,nomWid,systSignalList)
        _,altExp       = getMergedDists(systfIn,
                ('%s%s%s_%s_'%(lbCat,ch,cat,dist)),
                (rawSignalList if opt.addSigs else None),
                ('top' if opt.addSigs else ''),
                widList,nomWid,systSignalList)

        for systVar, procsToApply, normalize, useAltShape, projectRelToNom in sampleSysts:
            downShapes[systVar]={}
            upShapes[systVar]={}

            if sig == signalList[0] :
                procHasUncName="%s%s%s_%s"%(lbCat,ch,cat,dist)
                procHasUnc[procHasUncName][systVar] = []

            for iproc in procsToApply:
                if sig not in iproc and sig != signalList[0] : continue
                nomH=exp[iproc]

                #check which shape to use
                if useAltShape:

                    #get directly from another file
                    downH  = altExp[ procsToApply[iproc][0] ]
                    if len( procsToApply[iproc] ) > 1 :
                        upH = altExp[ procsToApply[iproc][1] ]
                    else:
                        #if only one variation is available, mirror it
                        upH = downH.Clone( '%s%sUp'%(iproc,systVar) )
                        for xbin in xrange(1,upH.GetNbinsX()+1):
                            diff=upH.GetBinContent(xbin)-nomH.GetBinContent(xbin)
                            upH.SetBinContent(xbin,nomH.GetBinContent(xbin)-diff)
                else:

                    #project from 2D histo (re-weighted from nominal sample)
                    ybinUp, ybinDown = -1, -1
                    for ybin in xrange(1,genVarShapes[ iproc ].GetNbinsY()+1):
                        label = genVarShapes[ iproc ].GetYaxis().GetBinLabel(ybin)
                        if procsToApply[iproc][0] in label : ybinDown=ybin
                        if procsToApply[iproc][1] in label : ybinUp=ybin

                    downH = genVarShapes[ iproc ].ProjectionX('%s%sDown'%(iproc,systVar), ybinDown, ybinDown)
                    upH   = genVarShapes[ iproc ].ProjectionX('%s%sUp'%(iproc,systVar),   ybinUp,   ybinUp)

                # use do down/up x nom to generate the variation, then mirror it
                if projectRelToNom:
                    ratioH=downH.Clone()
                    ratioH.Divide(upH)
                    for xbin in xrange(1,nomH.GetNbinsX()+1):
                        nomVal=nomH.GetBinContent(xbin)
                        varVal = ratioH.GetBinContent(xbin) * nomVal
                        upH.SetBinContent(xbin, varVal)
                        if upH.GetBinContent(xbin) <= 0 :
                            upH.SetBinContent(xbin,1e-010)
                        varVal = varVal- nomVal
                        downH.SetBinContent(xbin, nomVal-varVal)
                        if downH.GetBinContent(xbin) <= 0 :
                            downH.SetBinContent(xbin,1e-010)

                #normalize (shape only variation is considered)
                if normalize : downH.Scale( nomH.Integral()/downH.Integral() )
                if normalize : upH.Scale( nomH.Integral()/upH.Integral() )

                #check if variation is meaningful
                accept = acceptVariationForDataCard(nomH=nomH, upH=upH, downH=downH)
                if not accept : continue

                #save
                if sig == signalList[0] :
                    procHasUncName="%s%s%s_%s"%(lbCat,ch,cat,dist)
                    procHasUnc[procHasUncName][systVar]+=[iproc]

                if sig not in iproc : continue
                downShapes[systVar][iproc]=downH
                upShapes[systVar][iproc]=upH




        for systVar, normalize, useAltShape, projectRelToNom in genSysts:
            dirForSyst=('%sup%s%s%s_%s_'%(systVar,lbCat,lfs,cat,dist))
            if "MEmuF" in systVar : dirForSyst=('gen3%s%s%s_%s_'%(lbCat,lfs,cat,dist))
            if "MEmuR" in systVar : dirForSyst=('gen5%s%s%s_%s_'%(lbCat,lfs,cat,dist))
            if "MEtot" in systVar : dirForSyst=('gen6%s%s%s_%s_'%(lbCat,lfs,cat,dist))

            _,genVarShapesUp = getMergedDists((systfIn if useAltShape else fIn),
                    dirForSyst,
                    (rawSignalList if opt.addSigs else None),
                    ('top' if opt.addSigs else ''),
                    widList,nomWid,signalList)

            dirForSyst=('%sup%s%s%s_%s_'%(systVar,lbCat,lfs,cat,dist))
            if "MEmuF" in systVar : dirForSyst=('gen4%s%s%s_%s_'%(lbCat,lfs,cat,dist))
            if "MEmuR" in systVar : dirForSyst=('gen8%s%s%s_%s_'%(lbCat,lfs,cat,dist))
            if "MEtot" in systVar : dirForSyst=('gen10%s%s%s_%s_'%(lbCat,lfs,cat,dist))

            _,genVarShapesDn = getMergedDists((systfIn if useAltShape else fIn),
                    dirForSyst,
                    (rawSignalList if opt.addSigs else None),
                    ('top' if opt.addSigs else ''),
                    widList,nomWid,signalList)

            downShapes[systVar]={}
            upShapes[systVar]={}

            if sig == signalList[0] :
                procHasUncName="%s%s%s_%s"%(lbCat,ch,cat,dist)
                procHasUnc[procHasUncName][systVar] = []

            for iproc in exp:
                if sig not in iproc and sig != signalList[0] : continue
                nomH=exp[iproc]

                if not iproc in genVarShapesDn :
                    continue
                #get directly from another file
                downH  = genVarShapesDn[iproc]
                if iproc in genVarShapesUp:
                    upH = genVarShapesUp[iproc]
                else:
                    #if only one variation is available, mirror it
                    upH = downH.Clone( '%s%sUp'%(iproc,systVar) )
                    for xbin in xrange(1,upH.GetNbinsX()+1):
                        diff=upH.GetBinContent(xbin)-nomH.GetBinContent(xbin)
                        upH.SetBinContent(xbin,nomH.GetBinContent(xbin)-diff)

                # use do down/up x nom to generate the variation, then mirror it
                if projectRelToNom:
                    ratioH=downH.Clone()
                    ratioH.Divide(upH)
                    for xbin in xrange(1,nomH.GetNbinsX()+1):
                        nomVal=nomH.GetBinContent(xbin)
                        varVal = ratioH.GetBinContent(xbin) * nomVal
                        upH.SetBinContent(xbin, varVal)
                        varVal = varVal- nomVal
                        downH.SetBinContent(xbin, nomVal-varVal)

                #normalize (shape only variation is considered)
                if normalize : downH.Scale( nomH.Integral()/downH.Integral() )
                if normalize : upH.Scale( nomH.Integral()/upH.Integral() )

                #check if variation is meaningful
                accept = acceptVariationForDataCard(nomH=nomH, upH=upH, downH=downH)
                if not accept : continue

                #save
                if sig in iproc and 'V' not in iproc:
                    downShapes[systVar][iproc]=downH
                    upShapes[systVar][iproc]=upH

                if sig == signalList[0] :
                    procHasUncName="%s%s%s_%s"%(lbCat,ch,cat,dist)
                    print procHasUncName,',',systVar,',',iproc,',',upH.Integral()
                    procHasUnc[procHasUncName][systVar]+=[iproc]



        ######################
        # Produce morph grid #
        ######################

        # initialize grid with width and mass
        widthDim=ROOT.RooBinning(len(widths),0,len(widths)-1)
        massDim =ROOT.RooBinning(len(masses),0,len(masses)-1)
        refGrid =ROOT.RooMomentMorphND.Grid(widthDim,massDim)

        # add systematics binnings to grid
        nShapes=0
        for systVar in upShapes :
            if len(upShapes[systVar])==0 : continue
            binning=ROOT.RooBinning(3,-1,1,"%s"%systVar)
            refGrid.addBinning(binning)
            nShapes += 1

        massAxes=[ws.var('x').frame()]*len(widths)
        widAxes =[ws.var('x').frame()]*len(masses)
        for imass in xrange(0,len(masses)):
            mass,url,proc=masses[imass]

            if modSig not in replaceBadCharacters(proc) : continue

            tfIn=ROOT.TFile.Open(url)
            for iwid in xrange(0,len(widths)):

                #get histogram and convert to a PDF
                dirname='%s%s%s_%s_%3.1fw'%(lbCat,ch,cat,dist,widths[iwid])
                h=tfIn.Get(dirname+"/"+dirname+'_'+proc)
                if not isinstance(h, ROOT.TH1):
                    print "\t\tProblem with pdf ",dirName,"/",dirName,"_",proc,"  | Continuing..."
                    continue
                name='%s_m%d'%(dirname,int(10*mass))
                data=ROOT.RooDataHist(name,name,ROOT.RooArgList(ws.var("x")),h)
                pdf=ROOT.RooHistPdf(name+"_pdf",name+"_pdf",ROOT.RooArgSet(ws.var("x")),data)
                pdf.plotOn(massAxes[iwid],
                        ROOT.RooFit.Name("%3.1f"%mass),
                        ROOT.RooFit.LineColor(ROOT.kOrange+kOrangeList[imass*3]))
                pdf.plotOn(widAxes[imass],
                        ROOT.RooFit.Name("%3.1f"%widths[iwid]),
                        ROOT.RooFit.LineColor(ROOT.kOrange+kOrangeList[iwid]))
                getattr(ws,'import')(pdf,ROOT.RooCmdArg())

                #add pdf to the grid
                print 'Adding',pdf.GetName(),'@ (',iwid,',',imass,')'
                insertVector=ROOT.vector('int')(2+nShapes)
                insertVector[0]=iwid
                insertVector[1]=imass
                refGrid.addPdf(ws.pdf(pdf.GetName()),insertVector)

                isyst=0
                for systVar in upShapes :
                    if mass != 172.5 : break;
                    # get systematic histos
                    systHUp=None
                    systHDn=None
                    for sigwid in upShapes[systVar] :
                        if ("%3.1f"%(widths[iwid])).replace('.','p') in sigwid :
                            systHUp=  upShapes[systVar][sigwid]
                            systHDn=downShapes[systVar][sigwid]
                            break
                    if systHUp is None or systHDn is None: continue

                    # make up RooAbsPdf
                    upName='%s_m%d_%sUp'%(dirname,int(10*mass),systVar)
                    upData=ROOT.RooDataHist(upName,upName,ROOT.RooArgList(ws.var("x")),systHUp)
                    upPdf=ROOT.RooHistPdf(upName+"_pdf",upName+"_pdf",ROOT.RooArgSet(ws.var("x")),upData)
                    getattr(ws,'import')(upPdf,ROOT.RooCmdArg())

                    # make down RooAbsPdf
                    dnName='%s_m%d_%sDown'%(dirname,int(10*mass),systVar)
                    dnData=ROOT.RooDataHist(dnName,dnName,ROOT.RooArgList(ws.var("x")),systHDn)
                    dnPdf=ROOT.RooHistPdf(dnName+"_pdf",dnName+"_pdf",ROOT.RooArgSet(ws.var("x")),dnData)
                    getattr(ws,'import')(dnPdf,ROOT.RooCmdArg())

                    insertVectorShape=ROOT.vector('int')(2+nShapes)
                    insertVectorShape[0]=iwid
                    insertVectorShape[1]=imass

                    # add pdfs to grid, setting all other systs to nominal
                    insertVectorShape[2+isyst]=1
                    refGrid.addPdf(ws.pdf(upPdf.GetName()),insertVectorShape)
                    insertVectorShape[2+isyst]=-1
                    refGrid.addPdf(ws.pdf(dnPdf.GetName()),insertVectorShape)
                    isyst+=1

            tfIn.Close()

        ###############################
        # create input PDF validation #
        ###############################
        # all masses for one width
        for i in range(0,len(massAxes)) :
            if opt.noValidation : break
            canvas.cd()
            canvas.SetLogy()

            massAxes[i].GetXaxis().SetTitle("%s [GeV]"%distTitle)
            massAxes[i].GetYaxis().SetTitle("Probability density: lepton-jet pairs")
            massAxes[i].SetTitle("")
            massAxes[i].SetLineWidth(1)
            massAxes[i].Draw()

            leg=ROOT.TLegend(0.65,0.17,0.90,0.37)
            leg.SetHeader("#Gamma_{t} = %s#times#Gamma_{SM}"%rawWidList[i])
            for mass,_,_ in masses :
                leg.AddEntry("%3.1f"%(mass),"m_{t} = %3.1f GeV"%(mass),"L")
            leg.Draw()

            CMS_lumi.relPosX = 0.180
            CMS_lumi.extraText = "Simulation Preliminary"
            CMS_lumi.extraOverCmsTextSize=0.50
            CMS_lumi.lumiTextSize = 0.55
            CMS_lumi.CMS_lumi(canvas,4,0)
            canvas.Update()
            canvas.SaveAs("%s/massInputValidation_%s%s%s_%s_%s.pdf"%(opt.output,lbCat,ch,cat,dist,modWidList[i]))
            canvas.SetGrayscale(True)
            canvas.Update()
            canvas.SaveAs("%s/gScale__massInputValidation_%s%s%s_%s_%s.pdf"%(opt.output,lbCat,ch,cat,dist,modWidList[i]))
            canvas.SetGrayscale(False)

        # all widths for one mass
        canvas.Clear()
        for i in range(0,len(widAxes)) :
            if opt.noValidation : break
            canvas.cd()
            canvas.SetLogy()

            widAxes[i].GetXaxis().SetTitle("%s [GeV]"%distTitle)
            widAxes[i].GetYaxis().SetTitle("Probability density: lepton-jet pairs")
            widAxes[i].SetTitle("")
            widAxes[i].SetLineWidth(1)
            widAxes[i].Draw()

            leg=ROOT.TLegend(0.65,0.17,0.90,0.37)
            leg.SetHeader("m_{t} = %3.1f GeV"%masses[i][0])
            for wid in widths :
                leg.AddEntry("%3.1f"%(wid),"#Gamma_{t} = %3.1f#times#Gamma_{SM}"%(wid),"L")
            leg.Draw()

            CMS_lumi.relPosX = 0.180
            CMS_lumi.extraText = "Simulation Preliminary"
            CMS_lumi.extraOverCmsTextSize=0.50
            CMS_lumi.lumiTextSize = 0.55
            CMS_lumi.CMS_lumi(canvas,4,0)
            canvas.Update()
            canvas.SaveAs("%s/widthInputValidation_%s%s%s_%s_%s.pdf"%(opt.output,lbCat,ch,cat,dist,("%3.1f"%masses[i][0]).replace('.','p')))
            canvas.SetGrayscale(True)
            canvas.Update()
            canvas.SaveAs("%s/gScale__widthInputValidation_%s%s%s_%s_%s.pdf"%(opt.output,lbCat,ch,cat,dist,("%3.1f"%masses[i][0]).replace('.','p')))
            canvas.SetGrayscale(False)

        #################################
        # create pdf and save workspace #
        #################################
        ws.factory('alpha[%i,%i]'%(0,len(widths)-1))
        ws.factory('beta[%i,%i]'%(0,len(masses)-1))
        ws.var('alpha').setVal(1)
        ws.var('beta' ).setVal(1)
        wsArgList=ROOT.RooArgList( ws.var('alpha'), ws.var('beta') )
        for systVar in upShapes :
            if len(upShapes[systVar])==0 : continue
            ws.factory('%s[-1,1]'%(systVar))
            ws.var('%s'%systVar).setVal(0)
            wsArgList.add( ws.var('%s'%systVar) , True )
        pdf=ROOT.RooMomentMorphND('widmorphpdf','widmorphpdf',
                                  wsArgList,
                                  ROOT.RooArgList( ws.var('x') ),
                                  refGrid,
                                  ROOT.RooMomentMorphND.Linear)
        pdf.useHorizontalMorphing(False)
        getattr(ws,'import')(pdf,ROOT.RooCmdArg())

        # save workspace to shapes
        outFile='%s/shapes.root'%(opt.output)
        fOut=ROOT.TFile.Open(outFile,'UPDATE')
        fOut.cd()
        ws.Write()
        fOut.Close()
        print "Closing outfile"


        ###############################
        # Produce morphing animations #
        ###############################
        if not opt.allMorphs and not first : continue
        if opt.noValidation : continue
        ROOT.gStyle.SetOptStat(0)
        status=ROOT.TLatex()

        zLimList={ "mlb": 3.5e-04,
                "asens" : 3.5e-04,
                "bsens" : 3.5e-04 }
        n2DScan=60
        n3DScan=120

        from PIL import Image, ImageSequence
        from images2gif import writeGif

        # setup morph validation arrays
        rotScanAlphaImgs  = []
        rotScanBetaImgs   = []
        alphaScanImgs     = []
        betaScanImgs      = []
        sensAlphaScanImgs = []
        sensBetaScanImgs  = []

        print "here1"
        os.mkdir('%s/gifplots_%s%s%s_%s_%s'%(opt.output,lbCat,ch,cat,sig,dist))

        # NOTE: we want to reopen the outfile so that we see what Combine sees.
        #       This is inefficient but it is important validation!
        fOut=ROOT.TFile.Open(outFile)
        workspace=fOut.Get("sigws_%s%s%s_%s_%s"%(lbCat,ch,cat,modSig,dist))
        morphPDF= workspace.pdf("widmorphpdf")

        alphaVar= workspace.var("alpha")
        betaVar = workspace.var("beta")
        mlbVar  = workspace.var("x")

        # this is used separately, so we can use pdf
        dMorDAlpha=pdf.derivative(alphaVar)
        dMorDBeta =pdf.derivative(betaVar)
        dMorDMlb  =pdf.derivative(mlbVar)

        #import pdb
        #pdb.set_trace()

        print "Begin morph validation!"

        #################################
        # create morphed PDF validation #
        #################################
        histos=[mlbVar.frame()]*len(widths)
        for i in range(0,len(widths)) :
            print "Setting value of alpha!"
            alphaVar.setVal(widths[i])

            leg=ROOT.TLegend(0.65,0.17,0.90,0.37)
            leg.SetHeader("#Gamma_{t} = %s#times#Gamma_{SM}"%rawWidList[i])
            print "Validating morphed pdfs!"

            for imass in range(0,len(masses)) :
                print "Setting value of beta!"
                betaVar.setVal(imass)
                morphPDF.plotOn(histos[i],
                        ROOT.RooFit.Name("%3.1f"%masses[imass][0]),
                        ROOT.RooFit.LineColor(ROOT.kOrange+kOrangeList[imass*3]))
                leg.AddEntry("%3.1f"%(masses[imass][0]),"m_{t} = %3.1f GeV"%(masses[imass][0]),"L")

            histos[i].GetXaxis().SetTitle("%s [GeV]"%distTitle)
            histos[i].GetYaxis().SetTitle("Morphed density: lepton-jet pairs")
            histos[i].SetTitle("")
            histos[i].SetLineWidth(1)
            histos[i].Draw()
            leg.Draw()

            CMS_lumi.relPosX = 0.180
            CMS_lumi.extraText = "Simulation Preliminary"
            CMS_lumi.extraOverCmsTextSize=0.50
            CMS_lumi.lumiTextSize = 0.55
            CMS_lumi.CMS_lumi(canvas,4,0)

            canvas.SetLogy()
            canvas.Update()
            canvas.SaveAs("%s/morphValidation_%s%s%s_%s_%s.pdf"%(opt.output,lbCat,ch,cat,dist,modWidList[i]))
            canvas.SetGrayscale(True)
            canvas.Update()
            canvas.SaveAs("%s/gScale__morphValidation_%s%s%s_%s_%s.pdf"%(opt.output,lbCat,ch,cat,dist,modWidList[i]))
            canvas.SetGrayscale(False)

        ##########################################
        # create 3D plots that rotate the camera #
        ##########################################
        canvas.SetLogx(False)
        canvas.SetLogy(False)
        canvas.SetLogz(True)
        ROOT.gStyle.SetPalette(52)

        histAlpha = mlbVar.createHistogram("Morphing against width variations",alphaVar)
        morphPDF.fillHistogram(histAlpha,ROOT.RooArgList(mlbVar,alphaVar))
        histAlpha.GetYaxis().SetTitleOffset(1.1)
        histAlpha.GetYaxis().SetTitle("Grid width index")
        histAlpha.GetXaxis().SetTitleOffset(1.35)
        histAlpha.GetXaxis().SetTitle("%s [GeV]"%distTitle)
        histAlpha.SetTitle("")

        histBeta  = mlbVar.createHistogram("Morphing against mass variations",betaVar)
        morphPDF.fillHistogram(histBeta ,ROOT.RooArgList(mlbVar,betaVar))
        histBeta.GetYaxis().SetTitleOffset(1.1)
        histBeta.GetYaxis().SetTitle("Grid mass index")
        histBeta.GetXaxis().SetTitleOffset(1.35)
        histBeta.GetXaxis().SetTitle("%s [GeV]"%distTitle)
        histBeta.SetTitle("")

        CMS_lumi.relPosX = 0.180
        CMS_lumi.cmsTextOffset=0
        CMS_lumi.extraText = "Simulation Preliminary"
        CMS_lumi.extraOverCmsTextSize=0.50
        CMS_lumi.lumiTextSize = 0.55

        alphaVar.setVal(1)
        betaVar.setVal(1)

        for theta in xrange(0,359,int(360/n3DScan)) :
            # mlb vs. mass plot for 3D scan
            if theta == 0 : zLimList["mlb"]=histAlpha.GetMaximum()*1.1
            histAlpha.SetMaximum(zLimList["mlb"])
            histAlpha.Draw("SURF1")
            status.DrawLatexNDC(0.75,0.02,"m_{t}=%3.2f GeV"%(betaVar.getVal()*(maxMT-minMT)/betaVar.getMax()+minMT))
            ROOT.gPad.SetPhi(theta)

            CMS_lumi.CMS_lumi(canvas,4,0)
            if first:
                canvas.SetLeftMargin(canvas.GetRightMargin()*1.1)

            canvas.SaveAs("%s/gifplots_%s_%s/rotscanalpha_%s_%i.png"%(opt.output,sig,dist,sig,theta))
            canvas.Clear()

            # mlb vs. width plot for 3D scan
            histBeta.SetMaximum(zLimList["mlb"])
            histBeta.Draw("SURF1")
            status.DrawLatexNDC(0.75,0.02,"#Gamma_{t}=%3.2f GeV"%(alphaVar.getVal()*(maxGammaT-minGammaT)/alphaVar.getMax()+minGammaT))
            ROOT.gPad.SetPhi(theta)

            CMS_lumi.CMS_lumi(canvas,4,0)
            if first:
                canvas.SetLeftMargin(canvas.GetRightMargin()*1.1)

            canvas.SaveAs("%s/gifplots_%s_%s/rotscanbeta_%s_%i.png"%(opt.output,sig,dist,sig,theta))
            canvas.Clear()

            # save images to proper arrays
            imgAlpha=Image.open("%s/gifplots_%s_%s/rotscanalpha_%s_%i.png"%(opt.output,sig,dist,sig,theta))
            rotScanAlphaImgs+=[imgAlpha]
            imgBeta=Image.open("%s/gifplots_%s_%s/rotscanbeta_%s_%i.png"%(opt.output,sig,dist,sig,theta))
            rotScanBetaImgs+=[imgBeta]

        ##########################################
        # create 2D plots that scan in variables #
        ##########################################
        canvas.SetLogx(False)
        canvas.SetLogy(False)
        canvas.SetLogz(False)

        CMS_lumi.cmsTextOffset=0.1
        if first:
            canvas.SetRightMargin(canvas.GetRightMargin()*2)
            canvas.SetLeftMargin(canvas.GetLeftMargin()*2)

        for step in xrange(0,n2DScan) :
            alphaStep=alphaVar.getMin()+(alphaVar.getMax()-alphaVar.getMin())/n2DScan*step
            betaStep = betaVar.getMin()+( betaVar.getMax()- betaVar.getMin())/n2DScan*step

            alphaVar.setVal(alphaStep);
            betaVar.setVal( betaStep);

            # book histos
            histAlpha = mlbVar.createHistogram("Morphing against mass variations %i"%step,alphaVar)
            histBeta  = mlbVar.createHistogram("Morphing against width variations %i"%step,betaVar)
            sensHistAlpha = ROOT.TH2D("sensAlpha %i"%step,"Sensitivity against mass variations",
                    30,0,300,
                    n2DScan, alphaVar.getMin(), alphaVar.getMax())
            sensHistBeta  = ROOT.TH2D("sensBeta %i"%step,"Sensitivity against width variations",
                    30,0,300,
                    n2DScan, betaVar.getMin(), betaVar.getMax())

            morphPDF.fillHistogram(histAlpha,ROOT.RooArgList(mlbVar,alphaVar))
            morphPDF.fillHistogram(histBeta ,ROOT.RooArgList(mlbVar,betaVar))

            # mlb histo against mass
            histAlpha.SetTitle("")
            histAlpha.GetXaxis().SetTitle("%s [GeV]"%distTitle)
            histAlpha.GetYaxis().SetTitle("Grid width index")
            histAlpha.GetZaxis().SetRangeUser(0,zLimList["mlb"])
            histAlpha.Draw("COLZ")
            status.DrawLatexNDC(0.25,0.02,"m_{t}=%3.2f GeV"%(betaVar.getVal()*(maxMT-minMT)/betaVar.getMax()+minMT))

            CMS_lumi.CMS_lumi(canvas,4,0)
            canvas.SaveAs("%s/gifplots_%s_%s/betascan_%s_%i.png"%(opt.output,sig,dist,sig,step))
            canvas.Clear()

            # mlb histo against width
            histBeta.SetTitle("")
            histBeta.GetYaxis().SetTitle("Grid mass index")
            histBeta.GetXaxis().SetTitle("%s [GeV]"%distTitle)
            histBeta.GetZaxis().SetRangeUser(0,zLimList["mlb"])
            histBeta.Draw("COLZ")
            status.DrawLatexNDC(0.25,0.02,"#Gamma_{t}=%3.2f GeV"%(alphaVar.getVal()*(maxGammaT-minGammaT)/alphaVar.getMax()+minGammaT))

            CMS_lumi.CMS_lumi(canvas,4,0)
            canvas.SaveAs("%s/gifplots_%s_%s/alphascan_%s_%i.png"%(opt.output,sig,dist,sig,step))
            canvas.Clear()


            #####################################
            # create 2D local sensitivity scans #
            #####################################
            for ia,ib,ix in [(ia,ib,ix) for ia in range(1,n2DScan+1)
                    for ib in range(1,n2DScan+1)
                    for ix in range(1,31)] :
                if not opt.makeSens : break

                # get all the bin numbers
                gifBinA=sensHistAlpha.FindBin(ix*3-1.5,alphaStep)
                gifBinB=sensHistAlpha.FindBin(ix*3-1.5, betaStep)
                binNumA=sensHistAlpha.GetBin(ix,ia)
                binNumB= sensHistBeta.GetBin(ix,ib)

                # prepare alpha sensitivity plot
                alphaVar.setVal(sensHistAlpha.GetXaxis().GetBinCenter(binNumA));
                betaVar.setVal(  sensHistBeta.GetXaxis().GetBinCenter(gifBinB));
                mlbVar.setVal(  sensHistAlpha.GetYaxis().GetBinCenter(binNumA));
                f=pdf.getVal()
                if f!=0 :
                    dfdx=pdf.derivative(ws.var('x')).getVal()
                    sensHistAlpha.Fill(binNumA,(1/f)*(dfdx**2))

                # prepare beta sensitivity plots
                alphaVar.setVal(sensHistAlpha.GetXaxis().GetBinCenter(gifBinA));
                betaVar.setVal(  sensHistBeta.GetXaxis().GetBinCenter(binNumB));
                mlbVar.setVal(   sensHistBeta.GetYaxis().GetBinCenter(binNumB));
                f=pdf.getVal()
                if f!=0 :
                    dfdx=pdf.derivative(ws.var('x')).getVal()
                    sensHistBeta.Fill(binNumB,(1/f)*(dfdx**2))

            if opt.makeSens :
                ###################################
                # sensitivity of mlb against mass #
                ###################################
                if step == 0 : zLimList["asens"]=sensHistAlpha.GetMaximum()*2
                sensHistAlpha.SetTitle("")
                sensHistAlpha.GetXaxis().SetTitle("%s [GeV]"%distTitle)
                sensHistAlpha.GetYaxis().SetTitle("Generator-level width [GeV]")
                sensHistAlpha.GetZaxis().SetRangeUser(0,zLimList["asens"])
                sensHistAlpha.Draw("COLZ")
                status.DrawLatexNDC(0.25,0.02,"m_{t}=%3.2f GeV"%(betaStep))

                print sensHistAlpha.GetMaximum()

                CMS_lumi.CMS_lumi(canvas,4,0)
                canvas.SaveAs("%s/gifplots_%s_%s/sensbetascan_%s_%i.png"%(opt.output,sig,dist,sig,step))
                canvas.Clear()

                ####################################
                # sensitivity of mlb against width #
                ####################################
                if step == 0 : zLimList["bsens"]=sensHistBeta.GetMaximum()*2
                sensHistBeta.SetTitle("")
                sensHistBeta.GetXaxis().SetTitle("%s [GeV]"%distTitle)
                sensHistBeta.GetYaxis().SetTitle("Generator-level mass [GeV]")
                sensHistBeta.GetZaxis().SetRangeUser(0,zLimList["bsens"])
                sensHistBeta.Draw("COLZ")
                status.DrawLatexNDC(0.25,0.02,"#Gamma_{t}=%3.2f GeV"%(alphaStep))

                CMS_lumi.CMS_lumi(canvas,4,0)
                canvas.SaveAs("%s/gifplots_%s_%s/sensalphascan_%s_%i.png"%(opt.output,sig,dist,sig,step))
                canvas.Clear()

                # save images
                sensImgAlpha=Image.open("%s/gifplots_%s_%s/sensalphascan_%s_%i.png"%(opt.output,sig,dist,sig,step))
                sensAlphaScanImgs+=[sensImgAlpha]
                sensImgBeta=Image.open("%s/gifplots_%s_%s/sensbetascan_%s_%i.png"%(opt.output,sig,dist,sig,step))
                sensBetaScanImgs+=[sensImgBeta]


            # save images
            imgAlpha=Image.open("%s/gifplots_%s_%s/alphascan_%s_%i.png"%(opt.output,sig,dist,sig,step))
            alphaScanImgs+=[imgAlpha]
            imgBeta=Image.open("%s/gifplots_%s_%s/betascan_%s_%i.png"%(opt.output,sig,dist,sig,step))
            betaScanImgs+=[imgBeta]


        # make gif files from arrays
        writeGif("%s/gifplots_%s_%s/alphascan_%s.gif"%(opt.output,sig,dist,sig),alphaScanImgs,duration=.05)
        writeGif("%s/gifplots_%s_%s/betascan_%s.gif"%( opt.output,sig,dist,sig), betaScanImgs,duration=.05)
        writeGif("%s/gifplots_%s_%s/rotscanalpha_%s.gif"%(opt.output,sig,dist,sig),rotScanAlphaImgs,duration=.1)
        writeGif("%s/gifplots_%s_%s/rotscanbeta_%s.gif"%( opt.output,sig,dist,sig), rotScanBetaImgs,duration=.1)

        if opt.makeSens :
            writeGif("%s/gifplots_%s_%s/sensalphascan_%s.gif"%(opt.output,sig,dist,sig),sensAlphaScanImgs,duration=.05)
            writeGif("%s/gifplots_%s_%s/sensbetascan_%s.gif"%(opt.output,sig,dist,sig),sensBetaScanImgs,duration=.05)

        first=False

    print "\n Morphed dists created, workspace saved. Producing datacards: \n"

    ####################################
    # PRODUCE WIDTH/MASS SCAN DATACARD #
    ####################################
    fIn=ROOT.TFile.Open(opt.input)
    for lbCat,ch,cat,dist in [(lbCat,ch,cat,dist)
            for lbCat in (['highpt'] if not opt.allMorphs else lbCatList)
            for ch    in (['EM'] if not opt.allMorphs else lfsList)
            for cat   in (['2b'] if not opt.allMorphs else catList)
            for dist  in (['incmlb'] if not opt.allMorphs else distList)] :
        # get signal/backgrounds
        obs,exp=getMergedDists(fIn,('%s%s%s_%s_'%(lbCat,ch,cat,dist)),None,'',widList,nomWid,signalList)
        if not opt.truthDataset=="" :
            obs=makeMCTruthHist(opt.truthDataset,signalList,exp)

        # prepare datacard
        datacard=open('%s/datacard_widfit__%s%s%s_%s.dat'%(opt.output,lbCat,ch,cat,dist),'w')

        datacard.write('#\n')
        datacard.write('# Generated by %s with git hash %s for analysis category %s%s\n' % (getpass.getuser(),
            commands.getstatusoutput('git log --pretty=format:\'%h\' -n 1')[1],
            cat, ch) )
        datacard.write('#\n')

        # header
        datacard.write('imax *\n')
        datacard.write('jmax *\n')
        datacard.write('kmax *\n')

        # shapes locations
        datacard.write('-'*50+'\n')
        for proc in exp:
            isSig=False
            for sig in modWidList :
                if sig in proc and 'V' not in proc: isSig=True
            if isSig : continue
            datacard.write('shapes %10s * shapes.root %s%s%s_%s/$PROCESS %s%s%s_%s_$SYSTEMATIC/$PROCESS\n'%(proc,lbCat,ch,cat,dist,lbCat,ch,cat,dist))

        datacard.write('shapes data_obs   * shapes.root %s%s%s_%s/$PROCESS\n'%(lbCat,ch,cat,dist))

        for sig in signalList:
            datacard.write('shapes %10s * shapes.root sigws_%s%s%s_$PROCESS_%s:widmorphpdf\n'%(sig,lbCat,ch,cat,dist))
        datacard.write('-'*50+'\n')

        # bin names
        datacard.write('bin 1\n')
        datacard.write('observation %3.1f\n' % obs.Integral())

        # observations
        datacard.write('-'*50+'\n')
        datacard.write('\t\t\t %32s'%'bin')
        for i in xrange(0,len(exp)-(len(widths)-1)*len(signalList)): datacard.write('%15s'%'1')
        datacard.write('\n')
        datacard.write('\t\t\t %32s'%'process')
        for sig in signalList:
            datacard.write('%15s'%(sig))
        for proc in exp:
            isSig=False
            for sig in modWidList :
                if sig in proc and 'V' not in proc: isSig=True
            if isSig : continue
            datacard.write('%15s'%proc)
        datacard.write('\n')
        datacard.write('\t\t\t %32s'%'process')
        for i in xrange(0,len(signalList)) : datacard.write('%15s'%str(-1*(len(signalList)-1-i)))
        i=0
        for proc in exp:
            isSig=False
            for sig in modWidList :
                if sig in proc and 'V' not in proc: isSig=True
            if isSig : continue
            i=i+1
            datacard.write('%15s'%str(i))
        datacard.write('\n')
        datacard.write('\t\t\t %32s'%'rate')
        for sig in signalList :
            datacard.write('%15s'%('%3.2f'%(exp['%s%s'%(sig,modNomWid)].Integral())))
        for proc in exp:
            isSig=False
            for sig in modWidList :
                if sig in proc and 'V' not in proc: isSig=True
            if isSig : continue
            datacard.write('%15s'%('%3.2f'%exp[proc].Integral()))
        datacard.write('\n')
        datacard.write('-'*50+'\n')


        # rate uncertainties
        for syst,val,pdf,whiteList,blackList in rateSysts:

            datacard.write('%32s %8s'%(syst,pdf))
            entryTxt=''
            try:
                entryTxt='%15s'%('%3.3f/%3.3f'%(ROOT.TMath.Max(val[0],0.01),val[1]))
            except:
                entryTxt='%15s'%('%3.3f'%val)

            for sig in signalList:
                newSig='%s1p0w'%(sig)
                if (len(whiteList)==0 and not newSig in blackList) or newSig in whiteList:
                    datacard.write(entryTxt)
                else:
                    datacard.write('%15s'%'-')
            for proc in exp:
                isSig=False
                for sig in modWidList :
                    if sig in proc and 'V' not in proc: isSig=True
                if isSig : continue
                if (len(whiteList)==0 and not proc in blackList) or proc in whiteList:
                    datacard.write(entryTxt)
                else:
                    datacard.write('%15s'%'-')
            datacard.write('\n')


        # write sample systematics to datacard
        for systVar, procsToApply, normalize, useAltShape, projectRelToNom in sampleSysts:
            datacard.write('%26s shape'%systVar)
            for sigproc in signalList:
                if ("%s1p0w"%(sigproc)) in procsToApply and ("%s1p0w"%(sigproc)) in procHasUnc["%s%s%s_%s"%(lbCat,ch,cat,dist)][systVar] :
                    datacard.write('%15s'%'1')
                else:
                    datacard.write('%15s'%'-')
            for proc in exp:
                isSig=False
                for modwid in modWidList :
                    if modwid in proc : isSig=True
                if isSig : continue
                if proc in procsToApply and proc in procHasUnc["%s%s%s_%s"%(lbCat,ch,cat,dist)][systVar] :
                    datacard.write('%15s'%'1')
                else:
                    datacard.write('%15s'%'-')
            datacard.write('\n')


        # write gen-level systematics to datacard
        for systVar, normalize, useAltShape, projectRelToNom in genSysts:
            datacard.write('%26s shape'%systVar)
            for sigproc in signalList:
                if ("%s1p0w"%(sigproc)) in procHasUnc["%s%s%s_%s"%(lbCat,ch,cat,dist)][systVar]:
                    datacard.write('%15s'%'1')
                else:
                    datacard.write('%15s'%'-')
            for proc in exp:
                isSig=False
                for modwid in modWidList :
                    if modwid in proc : isSig=True
                if isSig : continue
                if proc in procHasUnc["%s%s%s_%s"%(lbCat,ch,cat,dist)][systVar] :
                    datacard.write('%15s'%'1')
                else:
                    datacard.write('%15s'%'-')
            datacard.write('\n')

        datacard.close()


"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
