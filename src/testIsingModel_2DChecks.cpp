#include "IsingModel.cpp"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph2D.h"
    

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

void testIsingModel_2DChecks() {
    std::cout<<"***********************************************"<<std::endl;
    std::cout<<"* HausdorffIsingModel: TEST                   *"<<std::endl;
    std::cout<<"*                                             *"<<std::endl;
    std::cout<<"* Runs the following tests on the Ising model *"<<std::endl;
    std::cout<<"* class:                                      *"<<std::endl;
    std::cout<<"*       - Known 2D exact solutions work       *"<<std::endl;
    std::cout<<"***********************************************"<<std::endl;

    // Declare initial model, output files
    IsingModel model;
    TFile *fOut = new TFile("IsingModel_TestOutput_2DChecks_F2.root","RECREATE");

    // Prepare 2D system
    std::cout<<"\n\n***********************************************"<<std::endl;
    std::cout<<"* Preparing the 2D lattice                    *"<<std::endl;
    std::cout<<"***********************************************"<<std::endl;
    model.setDebug             (true);
    model.setNumThreads        (40);
    model.setNumMCSteps        (10);
    model.setLatticeDepth      (4);
    model.setHausdorffMethod   ("SCALING");
    model.setMCMethod          ("METROPOLIS");
    model.setInteractionSigma  (0);   
    model.setHausdorffDimension(1.5);
    model.setNumMCSteps(60);
    model.setCouplingConsts(0,1);
    model.setTemperature(0.001);
    model.setup();
    model.randomizeSpins();
    std::cout<<"\t\t- Magnetization: "<<model.getMagnetization()<<std::endl;
        getTimeDelta();
    model.runMonteCarlo();
        getTimeDelta();
    model.status();
        getTimeDelta();


    std::vector<double> hausdorffDims;
    std::vector<double> temps;
    std::vector<double> magnetizations;
    std::vector<double> energies;

    for(double i=1.25; i < 1.75; i += 0.2) {
        model.setHausdorffDimension(i);

        for(double j=0.1; j < 5; j += 0.25) {
            model.setTemperature(j);
            
            double mag=0;
            double en=0;
            for(int k=0; k < 3; k++) {
                model.reset();
                model.setup();
                model.randomizeSpins();
                model.runMonteCarlo();

                mag+=model.getMagnetization()/3;
                en+=model.getEffHamiltonian()/3;
            }

            hausdorffDims.push_back(i);
            temps.push_back(j);
            magnetizations.push_back(mag);
            energies.push_back(en);
        }
    }


    std::cout<<"\n\n***********************************************"<<std::endl;
    std::cout<<"* Preparing validation plots                  *"<<std::endl;
    std::cout<<"***********************************************"<<std::endl;

    // Prepare magnetization graph
    TGraph2D *magGraph = new TGraph2D();
    for(int i=0; i < hausdorffDims.size(); i++) {
        magGraph->SetPoint(i,hausdorffDims.at(i),temps.at(i),magnetizations.at(i));
        //std::cout<<"("<<hausdorffDims.at(i)<<", "<<temps.at(i)<<", "
        //              <<magnetizations.at(i)<<")"<<std::endl;
    }

    gStyle->SetPalette(1);
    magGraph->SetTitle("Magnetization: #sigma = 0, J = 1");
    magGraph->Draw("SURF1");
    gPad->Update();
    magGraph->GetXaxis()->SetTitle("Hausdorff dimension");
    magGraph->GetYaxis()->SetTitle("Temperature (k_{B}T)");
    magGraph->GetZaxis()->SetTitle("Magnetization");
    magGraph->GetXaxis()->SetTitleOffset(1.5);
    magGraph->GetYaxis()->SetTitleOffset(2.2);
    magGraph->GetZaxis()->SetTitleOffset(1.5);
    gPad->Modified();
    gPad->SaveAs("2DCheck_MagGraph_F2.pdf");

    // Prepare energy graph
    TGraph2D *energyGraph = new TGraph2D();
    for(int i=0; i < hausdorffDims.size(); i++) {
        energyGraph->SetPoint(i,hausdorffDims.at(i),temps.at(i),energies.at(i));
    }

    gStyle->SetPalette(1);
    energyGraph->SetTitle("#beta H: #sigma = 0, J = 1");
    energyGraph->Draw("SURF1");
    gPad->Update();
    energyGraph->GetXaxis()->SetTitle("Hausdorff dimension");
    energyGraph->GetYaxis()->SetTitle("Temperature (k_{B}T)");
    energyGraph->GetZaxis()->SetTitle("#beta H");
    energyGraph->GetXaxis()->SetTitleOffset(1.5);
    energyGraph->GetYaxis()->SetTitleOffset(2.2);
    energyGraph->GetZaxis()->SetTitleOffset(1.5);
    gPad->Modified();
    gPad->SaveAs("2DCheck_energyGraph_F2.pdf");

  


    fOut->cd();
    magGraph->Write();
    energyGraph->Write();

    fOut->Close();

}