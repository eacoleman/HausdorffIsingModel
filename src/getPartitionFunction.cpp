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


void getPartitionFunction(Double_t HDIM, 
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
    std::cout<<"\t - Making model"<<std::endl;
    char name[256];
    sprintf(name, "%1.2f_d%i_%3.2ft_%1.2fs_%2.2fh_%2.2fj_%im_%i_PARTITION",
                HDIM,DEPTH,KBT,SIGMA,
                COUPLING_H,COUPLING_J,NMCSTEPS,NTHREADS);
    TFile *outFile = new TFile(TString(name).ReplaceAll(".","-")+".root","RECREATE");
    TTree *outTree = new TTree("HausdorffIsingModel","Partition function data for HausdorffIsingModel");

    Int_t    tmagInit          =0;
    Int_t    tnumSpins         =0;
    Int_t    tlatticeDepth     =0;
    Int_t    thausdorffSlices  =0;
    Double_t thausdorffSpacing =0;
    Double_t thausdorffDim     =0;
    Double_t teffHInit         =0;
    Double_t tZ                =0; 
    Double_t th                =0;
    Double_t tJ                =0;
    Double_t tsig              =0;
    Double_t tkbT              =0;

    outTree->Branch("m_o",      &tmagInit);
    outTree->Branch("Ham_o",    &teffHInit);
    outTree->Branch("Z",        &tZ);

    outTree->Branch("h",        &th);
    outTree->Branch("J",        &tJ);
    outTree->Branch("sigma",    &tsig);
    outTree->Branch("kbT",      &tkbT);

    outTree->Branch("hSlices",   &thausdorffSlices);
    outTree->Branch("hSpacing", &thausdorffSpacing);
    outTree->Branch("hDim",     &thausdorffDim);
    outTree->Branch("numSpins", &tnumSpins);
    outTree->Branch("depth",    &tlatticeDepth);

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
    std::cout<<"\t - Running model"<<std::endl;
    model.setup();
        getTimeDelta();
    model.setAllSpins(1);
        tmagInit =model.getMagnetization();
        teffHInit=model.getEffHamiltonian();
        getTimeDelta();
    tZ=model.getZ();
    std::cout<<"\t\t- Partition function is "<<tZ<<std::endl;
        getTimeDelta();
    model.status();
        getTimeDelta();

    /*
     *  Store the results
     */
    th               = model.getH();
    tJ               = model.getJ();
    tlatticeDepth    = model.getLatticeDepth();
    thausdorffDim    = model.getHausdorffDimension();
    thausdorffSlices = model.getHausdorffSlices();
    thausdorffSpacing= model.getHausdorffScale();
    tsig             = model.getInteractionSigma();
    tkbT             = model.getkbT();
    tnumSpins        = model.getNumSpins();

    outTree->Fill();

    /*
     *  Write the output 
     */
    outFile->cd();
    outTree->Write();
    outFile->Close();

}
