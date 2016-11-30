###########################################################
# HausdorffIsingModel Steer Script                        #
# Author: Evan Coleman, 2016                              #
#                                                         #
# Simulations of the Ising model as a function of the     #
# lattice Hausdorff dimension                             #
###########################################################

if [[ "$1" == "" ]] ; then
    echo "*****************************************"
    echo "* HausdorffIsingModel Steering          *"
    echo "* - MAKE                                *"
    echo "* - TEST                                *"
    echo "* - JOBS                                *"
    echo "* - GIF                                 *"
    echo "* - RUN                                 *"
    echo "* - WWW                                 *"
    echo "*****************************************"
    exit 0 
fi

# Sample settings (to be modified)
hList="1,-1"
jList="1,-1"
tList="0.01,0.10,0.50,1.00"
sigList="-1,-0.5,0,0.5,1,2"
mcStepsList="10000"
dimList="0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0"


######################### Nothing below should change #########################


case $1 in

#################################### MAKE #####################################
MAKE )
echo "Compiling simulation code:"

make
make test

;;

#################################### TEST #####################################
TEST )
echo "Testing the simulation code:"

root -l -b -q src/testIsingModel.cpp 

;;

#################################### JOBS #####################################
JOBS )
echo "Starting grid jobs:"

python/SubmitCondor.py \
    --indir ${PWD} \
    --outdir ${PWD}/output/ \
    --hList $hList \
    --jList $jList \
    --tList $tList \
    --sigList $sigList \
    --mcSList $mcSList \
    --dimList $dimList \
    --depthList $depthList

;;

#################################### GIF  #####################################
GIF )
echo "Making GIF animations:"

sh scripts/plotsToGIF.sh

;;

#################################### RUN  #####################################
RUN )
echo "Running simulation for configuration $2 $3 $4 $5 $6 $7 $8:"

runIsingModel $2 $3 $4 $5 $6 $7 $8

;;

#################################### WWW  #####################################
WWW )
echo "Sending plots to CERNWEB"

scp -r output/* ecoleman@lxplus.cern.ch:~/www/HausdorffIsingModel

;;
esac
