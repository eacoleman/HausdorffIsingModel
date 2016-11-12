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
    echo "* - JOBS                                *"
    echo "* - PLOT                                *"
    echo "*****************************************"
    exit 0 
fi


case $1 in

#################################### JOBS #####################################
JOBS )

;;

#################################### PLOT #####################################
PLOT )

;;
esac
