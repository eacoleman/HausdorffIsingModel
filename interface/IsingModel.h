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
#include <pthread.h>
#include <cstdlib>
#include <vector>
#include <iostream>

class IsingModel {
    public :
        // Constructors, destructor
        IsingModel();
        virtual ~IsingModel() = {};
        
        // Settings
        void setNumThreads        (const int num=1   );
        void setNumSpins          (const int num=-1  );
        void setHausdorffDimension(const double dim=1);
        void setInteractionSigma  (const double sig=1);    
        void setCouplingConsts    (const double H, const double J);    
        void setTemperature       (const double tkbT);

        const std::vector<int> getSpinArray();
        const std::vector<int> getLatticeDimensions();
        const int    getNumThreads()         = {return nThreads        };
        const int    getNumSpins()           = {return nSpins          };
        const double getHausdorffDimension() = {return hausdorffDim    };
        const double getInteractionSigma()   = {return interactionSigma};    
        
        // Observables
        const int    getMagnetization();
        const double getFreeEnergy(
                const std::vector<int>& flips=std::vector<int>());
        const double computePartitionFunction(
                const int start=0,
                const std::vector<int>& flips=std::vector<int>());

        // Shorthand definitions
        const double K()={return J/kbT                     };
        const double h()={return H/kbT                     };
        const int    m()={return getMagnetization()        };
        const double Z()={return computePartitionFunction()};

        // Simulation
        void setup();
        void runMonteCarlo();

    private :
        // Spins
        typedef struct {
           int S;
           bool active;
           std::vector<int> coords;
        } spin;

        // Settings
        std::vector<spin> spinArray;
        std::vector<int > latticeDimensions;
        int    nThreads=1;
        int    nSpins=0;
        double hausdorffDim=1;
        double interactionSigma=1;

        // Thermodynamic variables
        double kbT=1;
        double H=1;
        double J=1;
        double latticeSpacing=1;

        // Observables
        int magnetization = 0;
        
        // Simulation
        bool   hasBeenSetup=false;
        double computePartitionFunction();
         
};
