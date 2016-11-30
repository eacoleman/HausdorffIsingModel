#include "IsingModel.cpp"
#include "TFile.h"
#include "TGraph.h"
#include "TMultiGraph.h"
    

std::clock_t start = std::clock();
double getTimeDelta() {
    double value=((float) std::clock()-start)/1000000;
    std::cout<<"\t\t- Done. It took "<<value<<" s"<<std::endl;
    start = std::clock();
    return value; 
}

bool niceAssert(TString statement, bool isTrue) {
    std::cout<<statement.Data()<<": "
             <<(isTrue ? "SUCCESS" : "FAILED")
             <<std::endl;
    return isTrue;
}

void testIsingModel() {
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

    // Declare settings for 4D lattice check
    std::cout<<"***********************************************"<<std::endl;
    std::cout<<"* Constructing the 4D lattice                 *"<<std::endl;
    std::cout<<"***********************************************"<<std::endl;
    model.setDebug             (true);
    model.setNumThreads        (40);
    model.setNumMCSteps        (10);
    model.setLatticeDepth      (4);
    model.setHausdorffDimension(3.5);
    model.setHausdorffMethod   ("SCALING");
    model.setMCMethod          ("METROPOLIS");
    model.setInteractionSigma  (0);   
    model.setTemperature       (0.001);
    model.setCouplingConsts    (1,1); 
    model.setup();
        getTimeDelta();
    //model.getFreeEnergy();
        getTimeDelta();
    model.status();
        getTimeDelta();


    std::cout<<"\nTESTS:"<<std::endl;
    niceAssert("Number of spins is 2^p*n^(pd)",pow(2,ceil(model.getHausdorffDimension()))
                                               *pow(model.getHausdorffSlices(),
                                                    model.getLatticeDepth()
                                                    *ceil(model.getHausdorffDimension()))
                                                == model.getNumSpins());
    niceAssert("Model magnetization is +1 * nSpins", model.getMagnetization()
                                                      == model.getSpinArray().size());

    

    // Prepare 1D system
    std::cout<<"\n\n***********************************************"<<std::endl;
    std::cout<<"* Preparing the 1D lattice                    *"<<std::endl;
    std::cout<<"***********************************************"<<std::endl;
    model.reset();
    model.setHausdorffDimension(1);
    model.setup();
        getTimeDelta();
    model.runMonteCarlo();
        getTimeDelta();
    model.status();
        getTimeDelta();





    // Prepare 2D system
    std::cout<<"\n\n***********************************************"<<std::endl;
    std::cout<<"* Preparing the 2D lattice                    *"<<std::endl;
    std::cout<<"***********************************************"<<std::endl;
    model.reset();
    model.setHausdorffDimension(1.5);
    model.setNumMCSteps(1);
    model.setup();
        getTimeDelta();
    model.runMonteCarlo();
        getTimeDelta();
    model.status();
        getTimeDelta();

    TGraph *metGr = (TGraph*) model.getConvergenceGr()->Clone("Metropolis");
    metGr->SetTitle("Metropolis");
    metGr->SetMarkerStyle(20);
    metGr->SetMarkerColor(2);



    // Check 2D convergence for Heat Bath
    // and plot the change in Delta E
    std::cout<<"\n\n***********************************************"<<std::endl;
    std::cout<<"* Preparing the 2D lattice with Heat Bath     *"<<std::endl;
    std::cout<<"***********************************************"<<std::endl;
    model.reset();
    model.setMCMethod("HEATBATH");
    model.setup();
        getTimeDelta();
    model.runMonteCarlo();
        getTimeDelta();
    model.status();
        getTimeDelta();


    TGraph *hbtGr = (TGraph*) model.getConvergenceGr()->Clone("HeatBath");
    hbtGr->SetTitle("Heat Bath");
    hbtGr->SetMarkerStyle(20);
    hbtGr->SetMarkerColor(4);




    // Check 2D convergence for Hybrid
    // and plot the change in Delta E
    std::cout<<"\n\n***********************************************"<<std::endl;
    std::cout<<"* Preparing the 2D lattice with Hybrid Mode   *"<<std::endl;
    std::cout<<"***********************************************"<<std::endl;
    model.reset();
    model.setMCMethod("HYBRID");
    model.setup();
        getTimeDelta();
    model.runMonteCarlo();
        getTimeDelta();
    model.status();
        getTimeDelta();

    TGraph *hrbGr = (TGraph*) model.getConvergenceGr()->Clone("Hybrid");
    hrbGr->SetTitle("Hybrid");
    hrbGr->SetMarkerStyle(20);
    hrbGr->SetMarkerColor(6);


    // Prepare convergence graph
    TMultiGraph* convergenceGr = new TMultiGraph();
    convergenceGr->SetName("ConvergenceGraph");
    convergenceGr->Add(metGr);
    convergenceGr->Add(hbtGr);
    convergenceGr->Add(hrbGr);


    //convergenceGr->GetXaxis()->SetTitle("Monte Carlo step");
    //convergenceGr->GetYaxis()->SetTitle("#Delta F");


    fOut->cd();
    convergenceGr->Write();


    fOut->Close();

}
