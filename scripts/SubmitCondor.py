import os,sys
from optparse import OptionParser
import math

pwd=os.environ['PWD']

# get options
parser = OptionParser(description='Submit condor jobs for anaSubstructure calls.')
parser.add_option('--indir',       action='store', dest='indir',    default=pwd, help='Input LHE sample directory')
parser.add_option('--outdir',      action='store', dest='outdir',   default=pwd+"/output/",
        help='Location of output directory')
parser.add_option('--min',     action='store', dest='minJob',      default=-1,
        help='Job number to start at')
parser.add_option('--max',     action='store', dest='maxJob',      default=-1,
        help='Maximum number of jobs to perform')
parser.add_option('--hList',       action='store', dest='h',        default='',  help='List of config h fields')
parser.add_option('--jList',       action='store', dest='j',        default='',  help='List of config J fields')
parser.add_option('--tList',       action='store', dest='t',        default='',  help='List of config temps')
parser.add_option('--sigList',     action='store', dest='sig',      default='',  help='List of config sigmas')
parser.add_option('--mcStepsList', action='store', dest='mcsteps',  default='',  help='List of config # MC steps')
parser.add_option('--dimList',     action='store', dest='dim',      default='',  help='List of config dimensions')
parser.add_option('--depthList',   action='store', dest='depth',    default='',  help='List of config depths')
parser.add_option('--exe',         action='store', dest='ex',       default="src/runIsingModel.cpp",
        help='Location of ROOT macro')
parser.add_option('--outeos',      action='store', dest='outeos',
        default="root://cmseos.fnal.gov:///store/user/ecoleman/HausdorffIsingModel/output/",
        help='Location of output eos directory')

(options, args) = parser.parse_args()
cmssw_base = os.environ['CMSSW_BASE']

# check that input directory is specified
if options.indir == "":
    print "ERROR: no --indir specified. Exiting..."
    quit()

# check that output directory exists
options.outdir = os.path.abspath(options.outdir)
if not os.path.exists(options.outdir):
    os.system('mkdir -p ' + options.outdir + '/stdout/')

from math import ceil
#nFilesPerLHE=int(ceil(float(options.maxEvents)/float(options.evPerJob)))

hList=options.h.split(',')
jList=options.j.split(',')
tList=options.t.split(',')
sigList=options.sig.split(',')
dimList=options.dim.split(',')
depList=options.depth.split(',')
mcSteps=options.mcsteps.split(',')


# prepare all configs
nJob=0
for h,j,t,sig,mcsteps,dim,depth in [(a,b,c,d,e,f,g)
        for a in hList
        for b in jList
        for c in tList
        for d in sigList
        for e in mcSteps
        for f in dimList
        for g in depList] :
    print nJob
    if h == "" or j == "" or t=="" or sig=="" or mcsteps=="" or dim=="" or depth=="" : continue
    if nJob>float(options.maxJob) : break
    if nJob<float(options.minJob) :
        nJob+=1
        continue

    current_name = "dim"+dim+"_h"+h+"_j"+j+"_t"+t+"_s"+sig+"_m"+mcsteps+"_dep"+depth
    current_conf = open(current_name.replace('.','p') + '.condor','w')
    current_shel = open(current_name.replace('.','p') + '.sh','w')

    # write condor config
    conf_tmpl = open('./condor/CondorConf.tmpl.condor')
    for line in conf_tmpl:
        if 'CMSSWBASE'     in line: line = line.replace('CMSSWBASE',     cmssw_base)
        if 'OUTDIR'        in line: line = line.replace('OUTDIR',        options.outeos)
        if 'OUTPUT_PATH'   in line: line = line.replace('OUTPUT_PATH',   "%s/stdout/"%options.outdir)
        if 'INDIR'         in line: line = line.replace('INDIR',         "%s/"%options.indir)
        if 'PARAM_H'       in line: line = line.replace('EXEC',          options.ex)
        if 'PARAM_H'       in line: line = line.replace('PARAM_H',       h          )
        if 'PARAM_J'       in line: line = line.replace('PARAM_J',       j          )
        if 'PARAM_T'       in line: line = line.replace('PARAM_T',       t          )
        if 'PARAM_SIG'     in line: line = line.replace('PARAM_SIG',     sig        )
        if 'PARAM_MCSTEPS' in line: line = line.replace('PARAM_MCSTEPS', mcsteps    )
        if 'PARAM_DIM'     in line: line = line.replace('PARAM_DIM',     dim        )
        if 'PARAM_DEPTH'   in line: line = line.replace('PARAM_DEPTH',   depth      )
        if 'NAME'          in line: line = line.replace('NAME',
                "dim"+dim+"_h"+h+"_j"+j+"_t"+t+"_s"+sig+"_m"+mcsteps+"_dep"+depth)

        current_conf.write(line)

    conf_tmpl.close()

    # write shell script
    shel_tmpl = open('./condor/CondorShel.tmpl.sh')
    for line in shel_tmpl:
        if 'CMSSWBASE'     in line: line = line.replace('CMSSWBASE',     cmssw_base)
        if 'OUTDIR'        in line: line = line.replace('OUTDIR',        options.outeos)
        if 'OUTPUT_PATH'   in line: line = line.replace('OUTPUT_PATH',   "%s/stdout/"%options.outdir)
        if 'INDIR'         in line: line = line.replace('INDIR',         "%s/"%options.indir)
        if 'PARAM_H'       in line: line = line.replace('EXEC',          options.ex)
        if 'PARAM_H'       in line: line = line.replace('PARAM_H',       h          )
        if 'PARAM_J'       in line: line = line.replace('PARAM_J',       j          )
        if 'PARAM_T'       in line: line = line.replace('PARAM_T',       t          )
        if 'PARAM_SIG'     in line: line = line.replace('PARAM_SIG',     sig        )
        if 'PARAM_MCSTEPS' in line: line = line.replace('PARAM_MCSTEPS', mcsteps    )
        if 'PARAM_DIM'     in line: line = line.replace('PARAM_DIM',     dim        )
        if 'PARAM_DEPTH'   in line: line = line.replace('PARAM_DEPTH',   depth      )
        if 'NAME'          in line: line = line.replace('NAME',
                "dim"+dim+"_h"+h+"_j"+j+"_t"+t+"_s"+sig+"_m"+mcsteps+"_dep"+depth)

        current_shel.write(line)

    shel_tmpl.close()

    current_conf.close()
    current_shel.close()
    nJob+=1

# end loop so that Python writes files
# run jobs
nJob=0
for h,j,t,sig,mcsteps,dim,depth in [(a,b,c,d,e,f,g)
        for a in options.h.split(',')
        for b in options.j.split(',')
        for c in options.t.split(',')
        for d in options.sig.split(',')
        for f in options.mcsteps.split(',')
        for f in options.dim.split(',')
        for g in options.depth.split(',')] :
    if h == "" or j == "" or t=="" or sig=="" or mcsteps=="" or dim=="" or depth=="" : continue
    if nJob>options.max : break
    if nJob<options.min :
        nJob+=1
        continue

    current_name = "dim"+dim+"_h"+h+"_j"+j+"_t"+t+"_s"+sig+"_m"+mcsteps+"_dep"+depth
    current_conf = current_name.replace('.','p') + '.condor'

    os.chdir(options.outdir + '/')
    os.system('condor_submit ' + current_conf)
    nJob+=1
