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
files=(src src/interface) 
for dir in ${files[*]} ; do
    mkdir -p $_CONDOR_SCRATCH_DIR/$dir
for file in $(xrdfs root://cmseos.fnal.gov ls -u /store/user/ecoleman/HausdorffIsingModel/$dir/) ; do
    xrdcp $file ${_CONDOR_SCRATCH_DIR}/$dir/
done
done


cd ${_CONDOR_SCRATCH_DIR}

# Set permissions
chmod 777 * 

# Run model, plot
root -l -b -q "EXEC(PARAM_DIM, PARAM_DEPTH, PARAM_T, PARAM_SIG, PARAM_H, PARAM_J, PARAM_MCSTEPS, 40)"

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
