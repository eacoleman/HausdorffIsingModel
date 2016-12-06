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

# Sample settings 
t=("0.5" "10" "0.5")
h=("0" "2" "0.5")
j=("-5" "5" "0.5")
s=("0" "3" "0.5")
H=("0.05" "4" "0.05")

hList=$(awk "BEGIN{ for (i=${h[0]}; i <= ${h[1]}; i+= ${h[2]}) printf(\"%.2f,\",i); }")
jList=$(awk "BEGIN{ for (i=${j[0]}; i <= ${j[1]}; i+= ${j[2]}) printf(\"%.2f,\",i); }") 
tList=$(awk "BEGIN{ for (i=${t[0]}; i <= ${t[1]}; i+= ${t[2]}) printf(\"%.2f,\",i); }") 
sigList=$(awk "BEGIN{ for (i=${s[0]}; i <= ${s[1]}; i+= ${s[2]}) printf(\"%.2f,\",i); }")
dimList=$(awk "BEGIN{ for (i=${H[0]}; i <= ${H[1]}; i+= ${H[2]}) printf(\"%.2f,\",i); }")

depList="1,2,3,4"
mcStepsList="10,100,1000,10000"

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

echo h: ${hList}
echo j: ${jList}
echo t: ${tList}
echo s: ${sigList}
echo H: ${dimList}
echo d: ${depList}
echo m: ${mcStepsList}

python scripts/SubmitCondor.py \
    --indir ${PWD} \
    --outdir ${PWD}/output/ \
    --hList ${hList} \
    --jList ${jList} \
    --tList ${tList} \
    --sigList ${sigList} \
    --mcStepsList $mcStepsList \
    --dimList $dimList \
    --depthList ${depList} \
    --min 0 --max 100

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

mkdir -p output/LatticeTests/
cd output/

tar -czvf output.tar.gz *
scp output.tar.gz ecoleman@lxplus.cern.ch:~/www/HausdorffIsingModel/

;;
esac
