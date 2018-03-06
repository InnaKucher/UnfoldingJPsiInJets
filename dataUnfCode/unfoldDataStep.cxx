#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;

#include "THnSparse.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"
#include "TSystem.h"

#include "TRandom.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#endif

using namespace std;

void unfold(bool doPrompt = true, bool doMid = true, Int_t iterMax = 8, Int_t stepNumber = 1) {
  
  Int_t iterMin =1;
  Int_t iterDef = iterMax - 2;
  
  string testInputName = "";
  string trainInputName = "";
  string outputName = "";

  if(doPrompt && doMid){
    testInputName = "optDataUnf/unfInput/data/meas_data_prompt_mid_allBins.root";
    trainInputName = Form("optDataUnf/unfInput/step%i/response_4D_prompt_midRapidity.root",stepNumber);
    outputName = Form("optDataUnf/unfOutput/step%i/UnfoldedDistributions_Prompt_Mid_8iter.root",stepNumber);
  }
  if(doPrompt && !doMid){
    testInputName = "optDataUnf/unfInput/data/meas_data_prompt_fow_allBins.root";
    trainInputName = Form("optDataUnf/unfInput/step%i/response_4D_prompt_fwdRapidity.root",stepNumber);
    outputName = Form("optDataUnf/unfOutput/step%i/UnfoldedDistributions_Prompt_Fwd_8iter.root",stepNumber);
  }
    
  if(!doPrompt && doMid){
    testInputName = "optDataUnf/unfInput/data/meas_data_nonprompt_mid_allBins.root";
    trainInputName = Form("optDataUnf/unfInput/step%i/response_4D_nonprompt_midRapidity.root",stepNumber);
    outputName = Form("optDataUnf/unfOutput/step%i/UnfoldedDistributions_NonPrompt_Mid_8iter.root",stepNumber);
  }
  if(!doPrompt && !doMid){
    testInputName = "optDataUnf/unfInput/data/meas_data_nonprompt_fow_allBins.root";
    trainInputName = Form("optDataUnf/unfInput/step%i/response_4D_nonprompt_fwdRapidity.root",stepNumber);
    outputName = Form("optDataUnf/unfOutput/step%i/UnfoldedDistributions_NonPrompt_Fwd_8iter.root",stepNumber);
  }
  

  // data
  TFile *f_measured = new TFile(testInputName.c_str());
  f_measured->ls();
  TH2D *hMeasured = (TH2D*)f_measured->Get("h_Meas;1");

  //mc
  
  TFile *f = new TFile(trainInputName.c_str());
  RooUnfoldResponse *resp = (RooUnfoldResponse*)f->Get("resp;1");
  TH2F *hSmear = (TH2F*)f->Get("fh2RespDimM;1");
  TH2F *hTrue = (TH2F*)f->Get("fh2RespDimT;1");
    
  hSmear->Sumw2();
  hTrue->Sumw2();

  RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;//kNoError;//
  //RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kNoError;
  
  Int_t nIterTmp = iterMax - iterMin+1;
  const Int_t nIter = nIterTmp;
  TH2D *hReco[nIter];
  TH2D *hFolded[nIter];
  TMatrixD covmat[nIter];
  RooUnfoldBayes unfold[nIter];

  TH1D *hJetPtTrue;
  TH1D *hJetPtMeas;
  TH1D *hJetPtUnf[nIter];
  TH1D *hJetPtFol[nIter];
  
  for(Int_t iter = iterMin; iter<=iterMax; iter++) {

    // RooUnfoldBayes    unfold (resp, hMeas, iter);
    unfold[iter-iterMin] = RooUnfoldBayes(resp, hMeasured, iter);

    // unfold using response matrix
    
    hReco[iter-iterMin] = (TH2D*)unfold[iter-iterMin].Hreco(errorTreatment);
    hReco[iter-iterMin]->SetName(Form("hReco_Iter%d",iter));
    unfold[iter-iterMin].Print();

    // fold using just unfolded result
    
    hFolded[iter-iterMin] = (TH2D*)resp->ApplyToTruth(hReco[iter-iterMin]);
    hFolded[iter-iterMin]->Sumw2();
    hFolded[iter-iterMin]->SetName(Form("hFolded_Iter%d",iter));

    covmat[iter-iterMin].ResizeTo(unfold[iter-iterMin].Ereco(RooUnfold::kCovariance));
    covmat[iter-iterMin] = unfold[iter-iterMin].Ereco(RooUnfold::kCovariance);
    
    // jet pt distributions save here
    
    hJetPtUnf[iter-iterMin] = dynamic_cast<TH1D*>(hReco[iter-iterMin]->ProjectionY(Form("hJetPtUnf_Iter%d",iter)));
    hJetPtFol[iter-iterMin] = dynamic_cast<TH1D*>(hFolded[iter-iterMin]->ProjectionY(Form("hJetPtFol_Iter%d",iter)));

    hJetPtUnf[iter-iterMin]->Sumw2();
    hJetPtFol[iter-iterMin]->Sumw2();
    
  }

  // response matrix folded with true distribution (test sample)
  
  TH2D *hPriorFolded = (TH2D*)resp->ApplyToTruth(hTrue);
  hPriorFolded->SetName("hPriorFolded");
  
  
  Int_t min_PriorFolded = hPriorFolded->GetYaxis()->FindBin(25.0+0.00001);
  Int_t max_PriorFolded = hPriorFolded->GetYaxis()->FindBin(35.0-0.00001);
  TH1D *hPriorFolded_x = dynamic_cast<TH1D*>(hPriorFolded->ProjectionX("hPriorFolded_x",min_PriorFolded,max_PriorFolded));
   
  
  // 1d projections from unfolded 2d distribution
  TH1D *hRecoP[2];
  hRecoP[0] = dynamic_cast<TH1D*>(hReco[iterDef-iterMin+2]->ProjectionX("hRecoP_x"));
  hRecoP[1] = dynamic_cast<TH1D*>(hReco[iterDef-iterMin+2]->ProjectionY("hRecoP_y"));

  // 1d projections from folded 2d distribution (obtained using new unfolded distribution)
  TH1D *hFoldedP[2];
  hFoldedP[0] = dynamic_cast<TH1D*>(hFolded[iterDef-iterMin+2]->ProjectionX("hFoldedP_x"));
  hFoldedP[1] = dynamic_cast<TH1D*>(hFolded[iterDef-iterMin+2]->ProjectionY("hFoldedP_y"));

  // 1d projections of measured z and jet pt (test sample)
  TH1D *hSP[2];
  hSP[0] = dynamic_cast<TH1D*>(hMeasured->ProjectionX("hSP_x"));
  hSP[1] = dynamic_cast<TH1D*>(hMeasured->ProjectionY("hSP_y"));

  Int_t min_Measured = hMeasured->GetYaxis()->FindBin(25.0+0.00001);
  Int_t max_Measured = hMeasured->GetYaxis()->FindBin(35.0-0.00001);
  TH1D *hMeasured_x = dynamic_cast<TH1D*>(hMeasured->ProjectionX("hMeasured_x",min_Measured,max_Measured));
    
  // 1d projections of true z and jet pt (test sampel)
  TH1D *hTP[2];
  hTP[0] = dynamic_cast<TH1D*>(hTrue->ProjectionX("hTP_x"));
  hTP[1] = dynamic_cast<TH1D*>(hTrue->ProjectionY("hTP_y"));

  hSP[0]->Sumw2();
  hFoldedP[0]->Sumw2();

  hSP[1]->Sumw2();
  hFoldedP[1]->Sumw2();

  hTP[0]->Sumw2();
  hRecoP[0]->Sumw2();

  hTP[1]->Sumw2();
  hRecoP[1]->Sumw2();

  hJetPtTrue = (TH1D*)hTP[1]->Clone();
  hJetPtTrue->SetName("hJetPtTrue");
  hJetPtMeas = (TH1D*)hSP[1]->Clone();
  hJetPtMeas->SetName("hJetPtMeas");
  
  const Int_t nPtBins = 3;
  
  Double_t ptmin[nPtBins] = {15.,25.,35.};
  Double_t ptmax[nPtBins] = {25.,35.,45.};
    
  TH1D *hMUnf[nPtBins][nIter];
  TH1D *hMFol[nPtBins][nIter];
  TH1D *hMTru[nPtBins];
  TH1D *hMSme[nPtBins];
  TH1D *hMPri[nPtBins];

  for(Int_t i = 0; i<nPtBins; i++) {
    Int_t min = hTrue->GetYaxis()->FindBin(ptmin[i]+0.00001);
    Int_t max = hTrue->GetYaxis()->FindBin(ptmax[i]-0.00001);
    hMTru[i] = dynamic_cast<TH1D*>(hTrue->ProjectionX(Form("hMTru_%d",i),min,max));

    min = hMeasured->GetYaxis()->FindBin(ptmin[i]+0.00001);
    max = hMeasured->GetYaxis()->FindBin(ptmax[i]-0.00001);
    hMSme[i] = dynamic_cast<TH1D*>(hMeasured->ProjectionX(Form("hMSme_%d",i),min,max));
    
    for(Int_t iter = iterMin; iter<=iterMax; iter++) {
      Int_t min2 = hReco[iter-iterMin]->GetYaxis()->FindBin(ptmin[i]+0.00001);
      Int_t max2 = hReco[iter-iterMin]->GetYaxis()->FindBin(ptmax[i]-0.00001);
      hMUnf[i][iter-iterMin] = dynamic_cast<TH1D*>(hReco[iter-iterMin]->ProjectionX(Form("hMUnf_%d_Iter%d",i,iter),min2,max2));

      min2 = hFolded[iter-iterMin]->GetYaxis()->FindBin(ptmin[i]+0.00001);
      max2 = hFolded[iter-iterMin]->GetYaxis()->FindBin(ptmax[i]-0.00001);
      hMFol[i][iter-iterMin] = dynamic_cast<TH1D*>(hFolded[iter-iterMin]->ProjectionX(Form("hMFol_%d_Iter%d",i,iter),min2,max2));
    }
  }


  TFile *fout = new TFile(outputName.c_str(),"RECREATE");
  hTrue->Write("fh2TrueMC");
  hMeasured->Write("fh2MeasData");
  resp->Write();
  hPriorFolded->Write();

  hJetPtTrue->Write();
  hJetPtMeas->Write();
  
  for(Int_t iter = iterMin; iter<=iterMax; iter++) {
    hReco[iter-iterMin]->Write();
    hFolded[iter-iterMin]->Write();
    covmat[iter-iterMin].Write(Form("covmat%d",iter));
    unfold[iter-iterMin].Write(Form("unfold%d",iter));
    hJetPtUnf[iter-iterMin]->Write();
    hJetPtFol[iter-iterMin]->Write();
  }


  for(Int_t i = 0; i<nPtBins; i++) {
    hMSme[i]->Write();
    hMTru[i]->Write();
    for(Int_t iter = iterMin; iter<=iterMax; iter++) {
      if(hMUnf[i][iter-iterMin]) hMUnf[i][iter-iterMin]->Write();
      if(hMFol[i][iter-iterMin]) hMFol[i][iter-iterMin]->Write();
    }
  }

  
  fout->Write();
  fout->Close();
}

void unfoldDataStep(Int_t step = 1){

  //prompt
  unfold(true,true,8,step);
  unfold(true,false,8,step);

  //nonprompt
  unfold(false,true,8,step);
  unfold(false,false,8,step);
  
}
