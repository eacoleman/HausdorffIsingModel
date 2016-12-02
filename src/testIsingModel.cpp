#include "IsingModel.cpp"
#include "TFile.h"
#include "TCanvas.h"
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
    model.setHausdorffDimension(2.5);
    model.setHausdorffMethod   ("SCALING");
    model.setMCMethod          ("METROPOLIS");
    model.setInteractionSigma  (0);   
    model.setTemperature       (0.001);
    model.setCouplingConsts    (1,1); /*
    model.setup();
        getTimeDelta();
    //model.getEffHamiltonian();
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
                                                      == model.getSpinArray().size());*/

/*
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
*/




    // Prepare 2D system
    std::cout<<"\n\n***********************************************"<<std::endl;
    std::cout<<"* Preparing the 2D lattice                    *"<<std::endl;
    std::cout<<"***********************************************"<<std::endl;
    model.reset();
    model.setHausdorffDimension(1.5);
    model.setNumMCSteps(100);
    model.setCouplingConsts(1,-100);
    model.setTemperature(0.001);
    model.setup();
    model.randomizeSpins();
    std::cout<<"\t\t- Magnetization: "<<model.getMagnetization()<<std::endl;
        getTimeDelta();
    model.runMonteCarlo();
        getTimeDelta();
    model.status();
        getTimeDelta();

    TGraph *metConGr = (TGraph*) model.getConvergenceGr()->Clone("Metropolis");
    metConGr->SetTitle("Metropolis");
    metConGr->SetMarkerStyle(20);
    metConGr->SetLineColor(2);
    metConGr->SetFillStyle(0);
    metConGr->SetMarkerColor(2);



    // Check 2D convergence for Heat Bath
    // and plot the change in Delta E
    std::cout<<"\n\n***********************************************"<<std::endl;
    std::cout<<"* Preparing the 2D lattice with Heat Bath     *"<<std::endl;
    std::cout<<"***********************************************"<<std::endl;
    model.reset();
    model.setMCMethod("HEATBATH");
    model.setup();
    model.randomizeSpins();
    std::cout<<"\t\t- Magnetization: "<<model.getMagnetization()<<std::endl;
        getTimeDelta();
    model.runMonteCarlo();
        getTimeDelta();
    model.status();
        getTimeDelta();


    TGraph *hbtConGr = (TGraph*) model.getConvergenceGr()->Clone("HeatBath");
    hbtConGr->SetTitle("Heat Bath");
    hbtConGr->SetMarkerStyle(20);
    hbtConGr->SetLineColor(4);
    hbtConGr->SetFillStyle(0);
    hbtConGr->SetMarkerColor(4);




    // Check 2D convergence for Hybrid
    // and plot the change in Delta E
    std::cout<<"\n\n***********************************************"<<std::endl;
    std::cout<<"* Preparing the 2D lattice with Hybrid Mode   *"<<std::endl;
    std::cout<<"***********************************************"<<std::endl;
    model.reset();
    model.setMCMethod("HYBRID");
    model.setNumMCSteps(100);
    model.setup();
    model.randomizeSpins();
    std::cout<<"\t\t- Magnetization: "<<model.getMagnetization()<<std::endl;
        getTimeDelta();
    model.runMonteCarlo();
        getTimeDelta();
    model.status();
        getTimeDelta();

    TGraph *hrbConGr = (TGraph*) model.getConvergenceGr()->Clone("Hybrid");
    hrbConGr->SetTitle("Hybrid");
    hrbConGr->SetMarkerStyle(20);
    hrbConGr->SetLineColor(6);
    hrbConGr->SetFillStyle(0);
    hrbConGr->SetMarkerColor(6);



    std::cout<<"\n\n***********************************************"<<std::endl;
    std::cout<<"* Preparing validation plots                  *"<<std::endl;
    std::cout<<"***********************************************"<<std::endl;

    // Prepare drawing utils
    TCanvas *C = new TCanvas("cnv","cnv",800,600);
    C->cd();

    // Prepare convergence graph
    TMultiGraph* convergenceGr = new TMultiGraph();
    convergenceGr->SetName("ConvergenceGraph");
    convergenceGr->SetTitle("Convergence of MC Minimization");
    convergenceGr->Add(metConGr);
    convergenceGr->Add(hbtConGr);
    convergenceGr->Add(hrbConGr);
    convergenceGr->Draw("APL");
    gPad->Update();
    convergenceGr->GetXaxis()->SetTitle("Monte Carlo step");
    convergenceGr->GetYaxis()->SetTitle("#Sigma|#Delta(#betaH)|");
    convergenceGr->GetYaxis()->SetTitleOffset(1.5);
    gPad->Modified();

    //C->SetLogy();
    C->BuildLegend();
    C->SaveAs("ConvergenceGr.pdf");
    C->SaveAs("ConvergenceGr.png");


    

    fOut->cd();
    convergenceGr->Write();


    fOut->Close();

}
