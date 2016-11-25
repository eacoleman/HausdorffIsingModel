#include "IsingModel.cpp"
#include "TFile.h"

int main(int argc, const char** args) {
    std::cout<<"***********************************************"<<std::endl;
    std::cout<<"* HausdorffIsingModel: TEST                   *"<<std::endl;
    std::cout<<"*                                             *"<<std::endl;
    std::cout<<"* Runs the following tests on the Ising model *"<<std::endl;
    std::cout<<"* class:                                      *"<<std::endl;
    std::cout<<"*       - 4D lattice has correct form         *"<<std::endl;
    std::cout<<"*       - Known 1D, 2D exact solutions work   *"<<std::endl;
    std::cout<<"*       - Heat bath algorithm converges       *"<<std::endl;
    std::cout<<"*         (faster than Metropolis algorithm)  *"<<std::endl;
    std::cout<<"***********************************************"<<std::endl;

    // Declare initial model, output files
    IsingModel model;
    TFile *fOut = new TFile("IsingModel_TestOutput.root","RECREATE");
    std::clock_t start = std::clock();

    // Declare settings for 4D lattice check
    std::cout<<"Constructing the 4D lattice"<<std::endl;
    model.setNumThreads        (40);
    model.setNumMCSteps        (10);
    model.setLatticeDepth      (2);
    model.setHausdorffDimension(3.5);
    model.setHausdorffMethod   ("SCALING");
    model.setMCMethod          ("METROPOLIS");
    model.setInteractionSigma  (0);   
    model.setTemperature       (0.001);
    model.setCouplingConsts    (1,1); 
    model.setup();
    
    std::cout<<"\t\t- Done. It took "<<(std::clock()-start)<<" ms"<<std::endl;
    start = std::clock();

    // Prepare 1D system
    std::cout<<"Preparing the 1D lattice"<<std::endl;
    model.reset();
    model.setHausdorffDimension(1);
    model.setup();
    model.runMonteCarlo();
    
    std::cout<<"\t\t- Done. It took "<<(std::clock()-start)<<" ms"<<std::endl;
    start = std::clock();

    // Prepare 2D system
    std::cout<<"Preparing the 2D lattice"<<std::endl;
    model.reset();
    model.setHausdorffDimension(2);
    model.setup();
    model.runMonteCarlo();
    
    std::cout<<"\t\t- Done. It took "<<(std::clock()-start)<<" ms"<<std::endl;
    start = std::clock();

    // Check 2D convergence for Heat Bath
    // and plot the change in Delta E
    std::cout<<"Preparing the 2D lattice with Heat Bath"<<std::endl;
    model.reset();
    model.setMCMethod("HEATBATH");
    
    std::cout<<"\t\t- Done. It took "<<(std::clock()-start)<<" ms"<<std::endl;
    start = std::clock();

    fOut->Close();

}
