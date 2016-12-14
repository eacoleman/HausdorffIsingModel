#!/bin/sh -f

# Efficient PATH/DISPLAY settings 
export DISPLAY=0
export PATH=/bin:/usr/bin:/usr/local/bin:/usr/krb5/bin:/usr/afsws/bin:/usr/krb5/bin/aklog

# Establish ROOT environment 
cd CMSSWBASE/src/ 
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`

# Get resources and scripts
echo ""
echo ""
echo " --- COPYING FILES ---"
files=(src src/interface) 
for dir in ${files[*]} ; do
    mkdir -p $_CONDOR_SCRATCH_DIR/$dir
for file in $(xrdfs root://cmseos.fnal.gov ls -u /store/user/ecoleman/HausdorffIsingModel/$dir/) ; do
    xrdcp $file ${_CONDOR_SCRATCH_DIR}/$dir/
done
done


cd ${_CONDOR_SCRATCH_DIR}

# Run model, plot
echo ""
echo ""
echo " --- RUNNING EXE ---"
root -l -b -q "EXEC(PARAM_DIM, PARAM_DEPTH, PARAM_T, PARAM_SIG, PARAM_H, PARAM_J, PARAM_MCSTEPS, 40)"

# Copy results to output directory
echo ""
echo ""
echo " --- MOVING FILES ---"
ls -u $_CONDOR_SCRATCH_DIR
for file in $(ls -u $_CONDOR_SCRATCH_DIR/) ; do 
    if [[ "condor" =~ "$file" ]] ; then 
        continue
    fi
    echo "\nMoving $file"
    xrdcp $file OUTDIR
done

# Clear output to prevent unwanted transfers
rm -rf condor 
rm -rf interface 
rm -rf output 
rm -rf scripts 
rm -rf src 
rm *.{png,jpg,gif,root} 
rm */*.{png,jpg,gif,root} 
