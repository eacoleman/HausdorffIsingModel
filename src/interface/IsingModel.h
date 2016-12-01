/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * IsingModel.cpp                                                              *
 * Author: Evan Coleman and Brad Marston, 2016                                 *
 *                                                                             *
 * Class structure for Ising model simulation. Key characteristics:            *
 *  - Generates lattice of a given Hausdorff dimension and depth               *
 *  - Multithreaded computation of partition function                          *
 *  - Multithreaded Monte Carlo steps                                          *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
//#include <boost/thread/thread.hpp>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <iostream>
#include <cmath>
#include "TRandom3.h"
#include "TGraph.h"

class IsingModel {
    public :
        // Constructors, destructor
        IsingModel();
        virtual ~IsingModel();
        
        // Settings
        void setDebug             (const bool dbg   ) {debug = dbg;}
        void setNumThreads        (const int num    );
        void setNumMCSteps        (const int num    );
        void setLatticeDepth      (const int num    );
        void setHausdorffDimension(const double dim );
        void setHausdorffMethod   (char* const  hmtd);
        void setMCMethod          (char* const  mcmd);
        void setInteractionSigma  (const double sig );   
        void setTemperature       (const double tkbT);
        void setCouplingConsts    (const double H,
                                   const double J); 

        const std::vector<int> getSpinArray();
        const std::vector<int> getLatticeDimensions();
        const std::string      getHausdorffMethod() 
                                    {return hausdorffMethod ;}
        const std::string      getMCMethod() 
                                    {return mcMethod;}
        const int    getNumThreads()         {return nThreads        ;}
        const int    getNumSpins()           {return nSpins          ;}
        const int    getLatticeDepth()       {return latticeDepth    ;}
        const double getHausdorffDimension() {return hausdorffDim    ;}
        const double getHausdorffSlices()    {return hausdorffSlices ;}
        const double getHausdorffScale()     {return hausdorffScale  ;}
        const double getInteractionSigma()   {return interactionSigma;}   
        const double getNumMCSteps()         {return nMCSteps        ;}
        const std::vector<double> getMCInfo(){return mcInfo          ;}
        const std::vector<double> getHybridInfo(){return hybridInfo  ;}
        
        // Observables
        const int    getMagnetization();
        const double getFreeEnergy(
                const std::vector<int>& flips=std::vector<int>());
        const double getFreeEnergy(const int flip);
        const double computePartitionFunction(
                const int start=0,
                const std::vector<int>& flips=std::vector<int>());

        // Shorthand definitions
        const double K() {return J/kbT                     ;}
        const double h() {return H/kbT                     ;}
        const int    m() {return getMagnetization()        ;}
        const double Z() {return computePartitionFunction();}

        // Simulation
        void setup();
        void reset();
        void status();
        void runMonteCarlo();
        void randomizeSpins();
        void setAllSpins(const int direction=1);

        // Plots
        TGraph* getConvergenceGr();

    private :
        // Spins
        struct spin {
           int S;
           bool active;
           std::vector<double> coords;
           bool operator < (const spin& ts) const {
                for(int iS=0; iS < coords.size(); iS++) {
                    if(coords.at(iS) < ts.coords.at(iS)) return true;
                    else if(coords.at(iS) > ts.coords.at(iS)) return false;
                    else continue;
                }
                return false;
            }
           bool operator <= (const spin& ts) const {
                return !(*this > ts);
            }
           bool operator > (const spin& ts) const {
                for(int iS=0; iS < coords.size(); iS++) {
                    if(coords.at(iS) > ts.coords.at(iS)) return true;
                    else if(coords.at(iS) < ts.coords.at(iS)) return false;
                    else continue;
                }
                return false;
           }
           bool operator >= (const spin& ts) const {
                return !(*this < ts);
           }
        };

        typedef struct spin spin;
        
        // Settings
        bool   debug=false;
        std::vector<spin> spinArray;
        std::vector<int > latticeDimensions;
        int    latticeDepth=1;
        int    nThreads=1;
        int    nSpins=0;
        double interactionSigma=1;
        double hausdorffDim=1;
        double hausdorffSlices=2;
        double hausdorffScale=1/3;
        int    nMCSteps=10000;
        std::string hausdorffMethod="SCALING";
        std::string mcMethod="HEATBATH";

        // Thermodynamic variables
        double kbT=1;
        double H=1;
        double J=1;
        double freeEnergy=0;

        // Observables
        int magnetization = 0;
        
        // Simulation
        bool   hasBeenSetup=false;
        double metropolisStep(TRandom3* rNG);
        double heatBathStep(TRandom3* rNG);
        void   hybridStep(const double rng, const std::vector<int>& spinFlips);
        double getDistanceSq(const spin i1, const spin i2);
        void   nextPermutation(std::vector<int>& tvN, const int max);
        void   addSpins(const int depth, 
                        const std::vector<double>& x0, 
                        const std::vector<double>& x1);
        std::vector<double> mcInfo;
        std::vector<double> hybridInfo;
         
        // C++ utils
        void swap(spin *a, spin *b);
        void QuickSort(std::vector<spin>& vec, int left, int right);
        int QSPartition(std::vector<spin>& vec, int left, int right);
};
