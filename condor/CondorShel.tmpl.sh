#!/bin/tcsh -f

# Efficient PATH/DISPLAY settings 
setenv DISPLAY 0
setenv PATH /bin:/usr/bin:/usr/local/bin:/usr/krb5/bin:/usr/afsws/bin:/usr/krb5/bin/aklog

# Establish ROOT environment 
cd CMSSWBASE/src/ 
source /cvmfs/cms.cern.ch/cmsset_default.csh
eval `scramv1 runtime -csh`
rehash

# Get resources and scripts
cp -r INDIR/* ${_CONDOR_SCRATCH_DIR}/

cd ${_CONDOR_SCRATCH_DIR}

# Set permissions
chmod 777 * 

# Run model, plot
sh steerScript RUN_ISING_MODEL \
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
