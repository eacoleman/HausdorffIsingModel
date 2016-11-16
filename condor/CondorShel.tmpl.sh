#!/bin/sh -f

# Efficient PATH/DISPLAY settings 
export DISPLAY=0
export PATH=/bin:/usr/bin:/usr/local/bin:/usr/krb5/bin:/usr/afsws/bin:/usr/krb5/bin/aklog

# Establish ROOT environment 
cd CMSSWBASE/src/ 
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
rehash

# Get resources and scripts
files=(interface output scripts src steerScript.sh)
for file in ${files[*]} ; do
    cp -r INDIR/* ${_CONDOR_SCRATCH_DIR}/
done


cd ${_CONDOR_SCRATCH_DIR}

# Set permissions
chmod 777 * 

# Run model, plot
sh steerScript RUN \
    PARAM_H \
    PARAM_J \
    PARAM_T \
    PARAM_SIG \
    PARAM_MCSTEPS \
    PARAM_DIM \
    PARAM_DEPTH
sh steerScript BASE_GIFS

# Copy results to output directory
xrdcp */*.{png,jpg,root} OUTDIR

# Clear output to prevent unwanted transfers
rm -rf condor 
rm -rf interface 
rm -rf output 
rm -rf scripts 
rm -rf src 
rm *.{png,jpg,gif,root} 
rm */*.{png,jpg,gif,root} 
