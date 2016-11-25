/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * IsingModel.cpp                                                              *
 * Author: Evan Coleman and Brad Marston, 2016                                 *
 *                                                                             *
 * Class definitions for Ising model simulation. Key characteristics:          *
 *  - Generates lattice of a given Hausdorff dimension and depth               *
 *  - Multithreaded computation of partition function                          *
 *  - Multithreaded Monte Carlo steps                                          *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "interface/IsingModel.h"
#include "TRandom3.h"

// Constructors/destructors implemented simply
// (because of number of options)
IsingModel::~IsingModel() {};


/* (void) setNumThreads
 *    | How many threads to use at a time.
 *  I | (int) number of threads
 */
void IsingModel::setNumThreads(const int num) {
    if(nThreads < 1) return;
    nThreads = num;
    hasBeenSetup=false;
}


/* (void) setNumMCSteps 
 *    | How many MC steps to perform 
 *  I | (int) number of steps 
 */
void IsingModel::setNumMCSteps(const int num) {
    if(num < 1) return;
    nMCSteps = num;
    hasBeenSetup=false;
}


/* (void) setLatticeDepth
 *    | How many steps to simulate into the fractal lattice 
 *  I | (int) depth 
 */
void IsingModel::setLatticeDepth(const int num) {
    latticeDepth=num;
    hasBeenSetup=false;
}


/* (void) setInteractionSigma
 *    | The distance coupling exponent on the J_ij 
 *  I | (int) number of spins 
 */
void IsingModel::setInteractionSigma(const double sig) {
    interactionSigma=sig;
    hasBeenSetup=false;
}


/* (void) setHausdorffDimension
 *    | Set the lattice Hausdorff dimension.
 *  I | (double) dimension to use 
 */
void IsingModel::setHausdorffDimension(const double dim) {
    if(dim <= 0) return;
    hausdorffDim=dim;
    hasBeenSetup=false;
}


/* (void) setHausdorffMethod
 *    | Set the lattice Hausdorff scaling method.
 *  I | (double) method to use:
 *    |         - SCALING   = modify separation
 *    |         - SPLITTING = modify # divisions
 */
void IsingModel::setHausdorffMethod(char* const hmtd) {
    hausdorffMethod=hmtd;
    hasBeenSetup=false;
}


/* (void) setMCMethod 
 *    | Set the MC method.
 *  I | (double) method to use:
 *    |         - METROPOLIS (no multithread) 
 *    |         - HEATBATH   (multithread)
 */
void IsingModel::setMCMethod(char* const mcmd) {
    mcMethod=mcmd;
    hasBeenSetup=false;
}


/* (void) setCouplingConsts
 *    | Set the values of H,J in the hamiltonian 
 *  I | (double) value of H, magnetic field coupling 
 *    | (double) value of J, neighbor couplings
 */
void IsingModel::setCouplingConsts(const double tH, const double tJ) {
    H=tH;
    J=tJ;
    hasBeenSetup=false;
}


/* (void) setTemperature
 *    | Set the temperature of the system
 *  I | (double) value of k_B * T (>0) to use 
 */
void IsingModel::setTemperature(const double tkbT) {
    if (tkbT < 0) return;
    kbT=tkbT;
    hasBeenSetup=false;
}


/* (vector<int>) getSpinArray 
 *    | Returns an array of the spins (+1,-1, or 0)
 */
const std::vector<int> IsingModel::getSpinArray() {
    std::vector<int> spins(0);
    for(const auto &it : spinArray) {
        spins.push_back(it.active ? it.S : 0);
    }
    return spins;
}


/* (vector<int>) getLatticeDimensions 
 *    | Returns an array of the number of spins along
 *    | each lattice edge
 */
const std::vector<int> IsingModel::getLatticeDimensions() {
    return latticeDimensions;
}


/* (int) getMagnetization() 
 *    | Returns the magnetization of the lattice 
 */
const int IsingModel::getMagnetization() {
    int mag=0;
    for(const auto &it : spinArray) {
       mag += it.S*it.active;
    }
    return mag;
}


/* (double) getFreeEnergy 
 *    | Get the free energy of the system (no multithread)
 *  I | (int) single spin to flip 
 */
const double IsingModel::getFreeEnergy(const int flip) {
    std::vector<int> tflips;
    tflips.push_back(flip);
    return getFreeEnergy(tflips);
}


/* (double) getFreeEnergy 
 *    | Get the free energy of the system (no multithread)
 *  I | (vector<int> (default: empty)) array of spin indices to flip 
 */
const double IsingModel::getFreeEnergy(const std::vector<int>& flips) {
    double energy=0;
    for(int i=0; i<nSpins; i++) {
        spin s=spinArray.at(i);
        if (!s.active) continue;
        bool spinFlip=
            (std::find(flips.begin(),flips.end(),i)!=flips.end() ? -1 : 1); 

        energy += h()*s.S*spinFlip;

        for(int j=0; j<nSpins; j++) {
            spin ts=spinArray.at(j);
            if(!ts.active) continue;
            if(i==j)       continue;
            bool tspinFlip=
                (std::find(flips.begin(),flips.end(),j)!=flips.end() ? -1 : 1); 

            energy += K()*pow(getDistanceSq(s,ts),interactionSigma/2)
                *s.S*ts.S*spinFlip*tspinFlip;
        }
    }
    return energy;
}


/* (double) computePartitionFunction 
 *    | Get the partition function of the system (no multithread)
 *  I | (int (default: 0)) index to start the trace at 
 *    | (vector<int> (default: empty)) array of spin indices to flip 
 */
const double IsingModel::computePartitionFunction(const int start, const std::vector<int>& flips) {
    double Z=0;
    
    std::vector<int> newflips = flips; 
    newflips.push_back(start);

    // add e^-BH with S_i = +1 and S_i = -1
    Z+=exp(getFreeEnergy(flips));
    Z+=exp(getFreeEnergy(newflips));

    // branch into computations with history 
    // S_i = +1 and S_i = -1
    Z+=computePartitionFunction(start+1,newflips);
    Z+=computePartitionFunction(start+1,flips);

    return Z; 
}


/* (void) nextPermutation
 *    | For a vector of indices and some uniform maximum index,
 *    | get a new tuple of indices from the current state by 
 *    | adding one to the lowest index and "carrying" over.
 *    | Return tvN[0]=-1 when finished, for the loop end condition
 *  I | (std::vector<int>) the array of indices to permute 
 *    | (int) the maximum array index
 */
void IsingModel::nextPermutation(std::vector<int> tvN, const int max) {
    for(size_t i=0; i<tvN.size(); i++) {
        if(tvN.at(i) < max) {
            tvN.at(i)++;
            return;
        }
        else if(i< tvN.size()-1) {
            tvN.at(i)=0;
        }
        else if(i==tvN.size()-1) {
            tvN.at(0)=-1;
            return;
        }
    }
}

/* (void) addSpins 
 *    | Adds spins to the spinArray by isolating each smallest hypercube
 *    | making up the lattice at a given depth 
 *  I | (int) depth to build to
 *    | (double) vector of coordinates to start fractal at
 *    | (double) vector of coordinates to end fractal at
 */
void IsingModel::addSpins(const int depth,
        const std::vector<double>& x0, 
        const std::vector<double>& x1) {
        
    double delta    = abs(x1.at(0)-x0.at(0));
    //double sliceLen = hausdorffScale * delta;
    //double spaceLen = (delta - sliceLen*hausdorffSlices)/(hausdorffSlices-1);

    // Create an array of length dimension, p
    // Containing arrays of length depth, d
    // Which will contain our p*d indices 
    std::vector<int> vN(latticeDimensions.size()*depth);

    // Loop over all valid positions for a spin hypercube
    for(std::fill(vN.begin(),vN.end(),0); 
        vN.at(0) != -1; 
        nextPermutation(vN,hausdorffSlices-1)) {
        // Current position is lower corner of hypercube we produce 
        // (think in terms of the origin for the unit hypercube, [0,1]^p) 
        std::vector<double> cPos(latticeDimensions.size());
        for(size_t iDim=0; iDim < latticeDimensions.size(); iDim++) {
            cPos.at(iDim) += x0.at(iDim);
            
            for(int iDepth=0; iDepth < depth; iDepth++) {
                int depthVal = vN.at(iDim*depth+iDepth);
                double depthScale = pow(hausdorffScale,depth-iDepth)*delta;
                cPos.at(iDim) += (1 + (1/hausdorffScale-hausdorffSlices)/(hausdorffSlices-1))
                                 *depthScale
                                 *depthVal;
            }
        }

        // - place spins at each corner of a hypercube
        // - cube will have side length s^d and bottom corner at cPos
        // - loop will go through e.g. for p=2 (0,0) , (0,1) , (1,0) , (1,1)
        std::vector<int> cubePoints(latticeDimensions.size());
        for(std::fill(cubePoints.begin(),cubePoints.end(),0);
            cubePoints.at(0) != -1;
            nextPermutation(cubePoints,1)) {
           
            spin ts;
            ts.active=1;
            ts.S=1;
            ts.coords=std::vector<double>(latticeDimensions.size());
            for(size_t index; index < cubePoints.size(); index++) {
                ts.coords.at(index) = cubePoints.at(index) * pow(hausdorffScale,depth)*delta
                                        + cPos.at(index);
            }

            spinArray.push_back(ts);
            nSpins++;
        } 
    }
}


/* (void) setup 
 *    | Prepare the class object for simulation 
 */
void IsingModel::setup() {
    
    // Calculate the lattice dimensions from the input
    // Hausdorff dimension
    // (See my notes for derivation)
    if(hausdorffMethod=="SCALING") {
        hausdorffScale=pow(hausdorffSlices,-1/hausdorffDim);
    }

    // Check that the dimensions are reasonable
    if(hausdorffScale > 1/hausdorffSlices) {
        std::cout<<"ERROR: Invalid Hausdorff scaling"<<std::endl;
        exit(EXIT_FAILURE); 
    }
    
    // Keep track of the number of site coordinates along each axis
    for(int i=0; i<ceil(hausdorffDim); i++) {
        latticeDimensions.push_back(0);
    }

    // Generate the lattice array
    std::vector<double> x0;
    std::vector<double> x1;
    for(size_t i=0; i<latticeDimensions.size(); i++) {
        x0.push_back(0);
        x1.push_back(1);
    }

    addSpins(latticeDepth,x0,x1);

    hasBeenSetup=true;
}


/* (void) runMonteCarlo 
 *    | Run the Monte Carlo simulation (spin-flipping) 
 */
void IsingModel::runMonteCarlo() {
    if(!hasBeenSetup) {
        std::cout<<"ERROR: Object has not been setup!"<<std::endl;
        exit(EXIT_FAILURE); 
    }
    
    // Various utils 
    TRandom3* rNG = new TRandom3(); 

    freeEnergy=getFreeEnergy();
    std::vector<boost::thread*> threads;

    // Start performing MC steps
    for(int i=0; i < nMCSteps; i++) {
             if(mcMethod=="METROPOLIS")   
                 freeEnergy = metropolisStep(rNG->Uniform());
        else if(mcMethod=="NNMETROPOLIS") 
                 nnMetropolisStep();
        else if(mcMethod=="HEATBATH") {
            
            // Prepare threads
            int currentNThreads=0;
            int nSpinsPerThread = floor(nSpins/nThreads);
            
            // Create a vector of spin indices
            int popVectorSize=nSpins;
            std::vector<int> popVector(nSpins);
            for(int j=0; j < nSpins; j++) popVector.push_back(j);
            
            // Loop through threads
            for(int iThread=0; iThread < nThreads; iThread++) {

                // Generate array of spins to flip for thread
                // by ripping apart vector of spin indices 
                std::vector<int> spinFlips(nSpinsPerThread);
                for(int j=0; j < nSpinsPerThread; ){ 
                    int index=rNG->Integer(popVectorSize);
                    spinFlips.push_back(popVector.at(index));
                    popVector.erase(popVector.begin()+index);
                    popVectorSize--;
                }

                // Save thread info outside of scope
                // and submit the thread
                boost::thread * tThread
                    = new boost::thread(&IsingModel::heatBathStep,this,rNG->Uniform(),spinFlips);
                threads.push_back(tThread);

                currentNThreads++;
            }

            // Let them run in sequence
            for(int iThread; iThread<currentNThreads; iThread++) {
                threads.at(iThread)->join();
            }

            freeEnergy = getFreeEnergy();
        }

    }

    delete rNG;
}


/* (void) metropolisStep 
 *    | Perform one run over the lattice, attempting spin-flips 
 */
double IsingModel::metropolisStep(const double rng) {

    // loop over spins
    double newE=0; 
    for(int i=0; i < nSpins; i++) {
        double tE = getFreeEnergy(i);
        bool spinFlip=false;

        if(tE-freeEnergy<0) spinFlip = true;
        else spinFlip = (rng < exp((freeEnergy-tE)/kbT));

        if(spinFlip) {
            spinArray.at(i).S=-spinArray.at(i).S;
            mcInfo.push_back(tE-freeEnergy);
        }
    }

    return newE;
}


/* (void) nnMetropolisStep 
 *    | Perform one run over the lattice, attempting spin-flips 
 *    | SIMPLE: Nearest-neighbors only!
 */
void IsingModel::nnMetropolisStep() {
    // TODO: Implement or remove
}


/* (void) simpleMonteCarloStep 
 *    | Perform one run over the lattice, given a group of spins 
 */
void IsingModel::heatBathStep(double rng, std::vector<int> spinFlips) {

    double tE = getFreeEnergy(spinFlips);
    bool spinFlip=false;

    if(tE-freeEnergy<0) spinFlip=true;
    else spinFlip = (rng < exp((freeEnergy-tE)/kbT));
        
    if(spinFlip) {
        for(size_t i=0; i<spinFlips.size(); i++) {
            spinArray.at(i).S=-spinArray.at(i).S;
        }
    }

    mcInfo.push_back(tE-freeEnergy);
}
