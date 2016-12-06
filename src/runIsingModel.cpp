#include "IsingModel.cpp"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TTree.h"

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


void runIsingModel(Double_t HDIM, 
                   Int_t DEPTH, 
                   Double_t KBT, 
                   Double_t SIGMA, 
                   Double_t COUPLING_H, 
                   Double_t COUPLING_J, 
                   Int_t NMCSTEPS, 
                   Int_t NTHREADS) {
    /*
     *  Make the ntuple 
     */
    TString outName = TString(HDIM)+"_d"+TString(DEPTH)
                                   +"_t"+TString(KBT)
                                   +"_s"+TString(SIGMA)
                                   +"_h"+TString(COUPLING_H)
                                   +"_j"+TString(COUPLING_J)
                                   +"_m"+TString(NMCSTEPS)
                                   +"_"+TString(NTHREADS);
    TFile *outFile = new TFile(TString(outName)+".root","RECREATE");
    TTree *outTree = new TTree("HausdorffIsingModel","Simulated data for HausdorffIsingModel");

    Int_t    mag              =0;
    Int_t    magInit          =0;
    Int_t    numSpins         =0;
    Int_t    latticeDepth     =0;
    Int_t    numMCSteps       =0; 
    Int_t    hausdorffScale   =0;
    Double_t hausdorffSpacing =0;
    Double_t hausdorffDim     =0;
    Double_t hausdorffSpacing =0; 
    Double_t effH             =0;
    Double_t effHInit         =0;
    Double_t Z                =0; 
    Double_t h                =0;
    Double_t J                =0;
    Double_t sig              =0;
    Double_t kbT              =0;
    TString  MCMethod         ="METROPOLIS";

    outTree->Branch("m",        &mag);
    outTree->Branch("m_o",      &magInit);
    outTree->Branch("Ham",      &effH);
    outTree->Branch("Ham_o",    &effHInit);
    outTree->Branch("Z",        &Z);

    outTree->Branch("h",        &h);
    outTree->Branch("J",        &J);
    outTree->Branch("sigma",    &sig);
    outTree->Branch("kbT",      &kbT);

    outTree->Branch("hScale",   &hausdorffScale);
    outTree->Branch("hSpacing", &hausdorffSpacing);
    outTree->Branch("hDim",     &hausdorffDim);
    outTree->Branch("numSpins", &numSpins);
    outTree->Branch("depth",    &latticeDepth);

    outTree->Branch("numSteps", &numMCSteps);
    outTree->Branch("MCMethod", &hausdorffMethod);

    /*
     *  Make the model
     */
    IsingModel model;
    model.setDebug             (true);
    model.setNumThreads        (NTHREADS);
    model.setNumMCSteps        (NMCSTEPS);
    model.setLatticeDepth      (DEPTH);
    model.setHausdorffDimension(HDIM);
    model.setHausdorffMethod   ("SCALING");
    model.setMCMethod          ("METROPOLIS");
    model.setInteractionSigma  (SIGMA);   
    model.setTemperature       (KBT);
    model.setCouplingConsts    (COUPLING_H,COUPLING_J); 

    /*
     *  Run the model
     */
    model.setup();
        getTimeDelta();
    model.randomizeSpins();
        magInit =model.getMagnetization();
        effHInit=model.getEffHamiltonian();
        getTimeDelta();
    model.runMonteCarlo();
        getTimeDelta();
    model.status();
        getTimeDelta();

    /*
     *  Store the results
     */
    h               = model.H();
    J               = model.J();
    mag             = model.m();
    effH            = model.getEffHamiltonian();
    latticeDepth    = model.getLatticeDepth();
    hausdorffDim    = model.getHausdorffDimension();
    hausdorffSlices = model.getHausdorffSlices();
    hausdorffSpacing= model.getHausdorffScale();
    sig             = model.getInteractionSigma();
    kbT             = model.kbT();
    numMCSteps      = model.getNumMCSteps();
    numSpins        = model.getNumSpins();
    MCMethod        = TString(model.getMCMethod().data());

    outTree->Fill();

    /*
     *  Write the output 
     */
    outFile->cd();
    TGraph *convGr = (TGraph*) model.getConvergenceGr()->Clone();
    convGr->Write();
    outTree->Write();
    outFile->Close();




}
