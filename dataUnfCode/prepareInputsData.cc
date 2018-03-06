#include "TROOT.h"
#include "TKey.h"
#include "TFile.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TObjArray.h"
#include "THStack.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TMath.h"
#include "TCut.h"
#include "TPaletteAxis.h"
#include "TBranchElement.h"
#include "TPaveStats.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TProfile.h"
#include "THnSparse.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <algorithm>

using namespace std;

void prepare(bool doPrompt = false, bool doMid = true, Int_t stepNumber = 1){

  string filename = "";
  string outputfile = "";
  string filenamePrevStep = "";
  
  if(doPrompt){
    filename = "../../TreesForUnfolding/tree_MCJPSIPR_PP_NoBkg_AccEff_JEC.root";
    
    if(doMid) outputfile = Form("../unfInput/step%i/unfolding_4D_prompt_midRapidity.root",stepNumber);
    if(!doMid) outputfile = Form("../unfInput/step%i/unfolding_4D_prompt_fwdRapidity.root",stepNumber);
    
    if(stepNumber > 1 && doMid) filenamePrevStep = Form("../unfOutput/step%i/UnfoldedDistributions_Prompt_Mid_8iter.root",stepNumber-1);
    if(stepNumber > 1 && !doMid) filenamePrevStep = Form("../unfOutput/step%i/UnfoldedDistributions_Prompt_Fwd_8iter.root",stepNumber-1);
  }
  
    
  if(!doPrompt){
    filename = "../../TreesForUnfolding/tree_MCJPSINOPR_PP_NoBkg_AccEff_JEC.root";
    
    if(doMid) outputfile = Form("../unfInput/step%i/unfolding_4D_nonprompt_midRapidity.root",stepNumber);
    if(!doMid) outputfile = Form("../unfInput/step%i/unfolding_4D_nonprompt_fwdRapidity.root",stepNumber);
    
    if(stepNumber > 1 && doMid) filenamePrevStep = Form("../unfOutput/step%i/UnfoldedDistributions_NonPrompt_Mid_8iter.root",stepNumber-1);
    if(stepNumber > 1 && !doMid) filenamePrevStep = Form("../unfOutput/step%i/UnfoldedDistributions_NonPrompt_Fwd_8iter.root",stepNumber-1);
  }
  
  TFile *file = new TFile(filename.c_str());
  TTree *t_unf = (TTree*)file->Get("treeForUnfolding");

  int n_entries = t_unf->GetEntries();

  int nbins_z = 5;
  float min_z = 0.;
  float max_z = 1.0;
    
  Double_t *z_frac_lowJetPtBin = new Double_t[nbins_z];
  Double_t *z_frac_centJetPtBin = new Double_t[nbins_z];
  Double_t *z_frac_highJetPtBin = new Double_t[nbins_z];
  
  // if this is step 2, 3, 4 ... take the z unfolded and use it for the prior, else put z prior to flat
  
  if(stepNumber > 1) {
    TFile *filePrevStep = new TFile(filenamePrevStep.c_str());
    TH1D *h_z_unf_lowJetPtBin = (TH1D*)filePrevStep->Get("hMUnf_0_Iter3;1");
    TH1D *h_z_unf_centJetPtBin = (TH1D*)filePrevStep->Get("hMUnf_1_Iter3;1");
    TH1D *h_z_unf_highJetPtBin = (TH1D*)filePrevStep->Get("hMUnf_2_Iter3;1");

    TH1D *h_z_allBins = (TH1D*)h_z_unf_lowJetPtBin->Clone();
    h_z_allBins->Add(h_z_unf_centJetPtBin);
    h_z_allBins->Add(h_z_unf_highJetPtBin);

    h_z_allBins->Scale(1/h_z_allBins->Integral());
    
    h_z_unf_lowJetPtBin->Scale(1/h_z_unf_lowJetPtBin->Integral());
    h_z_unf_centJetPtBin->Scale(1/h_z_unf_centJetPtBin->Integral());
    h_z_unf_highJetPtBin->Scale(1/h_z_unf_highJetPtBin->Integral());
        
    for(int i = 1; i <= nbins_z; i++){
      z_frac_lowJetPtBin[i-1] = h_z_unf_lowJetPtBin->GetBinContent(i);
      z_frac_centJetPtBin[i-1] = h_z_unf_centJetPtBin->GetBinContent(i);
      z_frac_highJetPtBin[i-1] = h_z_unf_highJetPtBin->GetBinContent(i);
    }

  }
  else{
    for(int i = 1; i <= nbins_z; i++){
      z_frac_lowJetPtBin[i-1] = 1.0;
      z_frac_centJetPtBin[i-1] = 1.0;
      z_frac_highJetPtBin[i-1] = 1.0;
    }
    
  }
  
  //J/Psi variables:
  
  float jp_pt = 0.;
  float jp_eta = 0.;
  float jp_mass = 0.;
     
  t_unf->SetBranchAddress("jp_pt", &jp_pt);
  t_unf->SetBranchAddress("jp_eta", &jp_eta);
  t_unf->SetBranchAddress("jp_mass", &jp_mass);

  //jet reco variables:

  float jt_pt = 0.;
  float jt_eta = 0.;
  
  t_unf->SetBranchAddress("jt_pt", &jt_pt);
  t_unf->SetBranchAddress("jt_eta", &jt_eta);

  float corr_AccEff = 0.;
  float corr_AccEff_weight =0.;
  t_unf->SetBranchAddress("corr_AccEff", &corr_AccEff);
    
  //z corrected

  float z = 0.;
  
  t_unf->SetBranchAddress("z", &z);

  //jet gen variables

  float jt_ref_pt = 0.;
  
  t_unf->SetBranchAddress("jt_ref_pt", &jt_ref_pt);

  //gen z:

  float gen_z = 0.;
  
  t_unf->SetBranchAddress("gen_z", &gen_z);

  // histograms for renormalizing response matrix:
   
  TH2F * hRespZJetPtCentBin = new TH2F ("hRespZJetPtCentBin","z_gen vs z_raw for 25 < refpt < 40",nbins_z,min_z,max_z,nbins_z,min_z,max_z);
  TH2F * hRespZJetPtLowBin = new TH2F ("hRespZJetPtLowBin","z_gen vs z_raw for 10 < refpt < 25",nbins_z,min_z,max_z,nbins_z,min_z,max_z);
  TH2F * hRespZJetPtHighBin = new TH2F ("hRespZJetPtHighBin","z_gen vs z_raw for 40 < refpt < 55",nbins_z,min_z,max_z,nbins_z,min_z,max_z);

  TH1D * h_z_gen = new TH1D ("h_z_gen","z_gen for 25 < refpt < 55",nbins_z,min_z,max_z);
  TH1D * hGenZJetPtCentBin = new TH1D ("hGenZJetPtCentBin","z_gen for 25 < refpt < 40",nbins_z,min_z,max_z);
  TH1D * hGenZJetPtLowBin = new TH1D ("hGenZJetPtLowBin","z_gen for 10 < refpt < 25",nbins_z,min_z,max_z);
  TH1D * hGenZJetPtHighBin = new TH1D ("hGenZJetPtHighBin","z_gen for 40 < refpt < 55",nbins_z,min_z,max_z);
    
  int nbins_jetpt = 3;
  float min_jetpt = 15.;
  float max_jetpt = 45.;

  TH2F * h_jetpPtGen_jetPtReco = new TH2F ("h_jetpPtGen_jetPtReco","gen vs reco corr jet pt",nbins_jetpt,min_jetpt,max_jetpt,nbins_jetpt,min_jetpt,max_jetpt);

  hRespZJetPtCentBin->Sumw2();
  hRespZJetPtLowBin->Sumw2();
  hRespZJetPtHighBin->Sumw2();
  h_jetpPtGen_jetPtReco->Sumw2();

  // 4D histogram creation
  
  int fDim = 4.;
  Double_t fValue[fDim];
  
  Int_t* bins_sparce = new Int_t[fDim];
  Double_t *xmin_sparce = new Double_t[fDim];
  Double_t *xmax_sparce = new Double_t[fDim];

  bins_sparce[0] = 3;
  xmin_sparce[0] = 15.0;
  xmax_sparce[0] = 45.0;

  bins_sparce[1] = 5;
  xmin_sparce[1] = 0.;
  xmax_sparce[1] = 1.0;

  bins_sparce[2] = 3;
  xmin_sparce[2] = 15.0;
  xmax_sparce[2] = 45.0;

  bins_sparce[3] = 5;
  xmin_sparce[3] = 0.;
  xmax_sparce[3] = 1.0;

  //initial not normalized 4D  
  THnSparseF * fSparse = new THnSparseF("hs", "hs", fDim, bins_sparce, xmin_sparce, xmax_sparce);
  fSparse->Sumw2();
  fSparse->CalculateErrors();

  //4D with z gen flat
  THnSparseF * fSparse_norm = new THnSparseF("hs_norm_z", "hs_norm_z", fDim, bins_sparce, xmin_sparce, xmax_sparce);
  fSparse_norm->Sumw2();
  fSparse_norm->CalculateErrors();

  //4D with z gen flat and gen pt flat
  THnSparseF * fSparse_norm_genJetPt = new THnSparseF("hs_norm_z_genJetPt", "hs_norm_z_genJetPt", fDim, bins_sparce, xmin_sparce, xmax_sparce);
  fSparse_norm_genJetPt->Sumw2();
  fSparse_norm_genJetPt->CalculateErrors();

  float jp_pt_cut = 0.;
  
  for(int nEv = 0; nEv < n_entries; nEv++){
     
    t_unf->GetEntry(nEv);
    
    // check event content

    if(doMid) jp_pt_cut = 6.5;
    if(!doMid) jp_pt_cut = 3.0;
    
    if(jp_pt < jp_pt_cut || jp_pt > 35.0) continue;
    
    if(doMid && TMath::Abs(jp_eta) > 1.6) continue;
    
    if(!doMid){
      if(TMath::Abs(jp_eta) < 1.6 || TMath::Abs(jp_eta) > 2.4) continue;
    }
            
    if(jp_mass < 2.6 || jp_mass > 3.5) continue;

    if(jt_pt < 15.0 || jt_pt > 45.0) continue;
    if(jt_ref_pt < 15.0 || jt_ref_pt > 45.0) continue;
    if(TMath::Abs(jt_eta) > 2.4) continue;
    
    double gen_z2 = gen_z;
    if(gen_z == 1.) gen_z2 = 0.99;

    double z2 = z;
    if(z == 1.) z2 = 0.99;

    if(gen_z2 >=1. ) continue;
    if(z2 >= 1. ) continue;
        
    if(jt_ref_pt >= 25.0 && jt_ref_pt < 35.){
      hRespZJetPtCentBin->Fill(gen_z2,z2);
      hGenZJetPtCentBin->Fill(gen_z2);
    }
    
    if(jt_ref_pt >= 15.0 && jt_ref_pt < 25.0){
      hRespZJetPtLowBin->Fill(gen_z2,z2);
      hGenZJetPtLowBin->Fill(gen_z2);
    }
    
    if(jt_ref_pt >= 35.0 && jt_ref_pt < 45.0){
      hRespZJetPtHighBin->Fill(gen_z2,z2);
      hGenZJetPtHighBin->Fill(gen_z2);
    }
        	  
    h_jetpPtGen_jetPtReco->Fill(jt_ref_pt,jt_pt);
    h_z_gen->Fill(gen_z2);
    
    //response matrix fill
    fValue[0] = jt_ref_pt;
    fValue[1] = gen_z2;
    fValue[2] = jt_pt;
    fValue[3] = z2;

    corr_AccEff_weight = 1/corr_AccEff;
    fSparse->Fill(fValue,corr_AccEff_weight);
    
  }

  // check 4D content

  cout << "number of bins :" << fSparse->GetNbins() <<endl;

  Double_t* x = new Double_t[fDim + 1];
  Double_t* x0 = new Double_t[fDim + 1];
  Int_t *bins = new Int_t[fDim];

  Double_t* vals_in_zGen_lowJetPtBin = new Double_t[5];
  Double_t* vals_in_zGen_centJetPtBin = new Double_t[5];
  Double_t* vals_in_zGen_highJetPtBin = new Double_t[5];
  Double_t* vals_in_getJetPtBin_init = new Double_t[3];

  vals_in_getJetPtBin_init[0] = 0.;
  vals_in_getJetPtBin_init[1] =0.;
  vals_in_getJetPtBin_init[2] =0.;
  
  vals_in_zGen_lowJetPtBin[0] = 0.;
  vals_in_zGen_lowJetPtBin[1] = 0.;
  vals_in_zGen_lowJetPtBin[2] = 0.;
  vals_in_zGen_lowJetPtBin[3] = 0.;
  vals_in_zGen_lowJetPtBin[4] = 0.;

  vals_in_zGen_centJetPtBin[0] = 0.;
  vals_in_zGen_centJetPtBin[1] = 0.;
  vals_in_zGen_centJetPtBin[2] = 0.;
  vals_in_zGen_centJetPtBin[3] = 0.;
  vals_in_zGen_centJetPtBin[4] = 0.;

  vals_in_zGen_highJetPtBin[0] = 0.;
  vals_in_zGen_highJetPtBin[1] = 0.;
  vals_in_zGen_highJetPtBin[2] = 0.;
  vals_in_zGen_highJetPtBin[3] = 0.;
  vals_in_zGen_highJetPtBin[4] = 0.;

  cout << "fSparse->GetNbins() = " << fSparse->GetNbins() << endl;

  // put z gen to be flat in each gen jet pt bin
  
  for (Long64_t count = 0; count < fSparse->GetNbins(); ++count) {
    x[fDim] = fSparse->GetBinContent(count, bins);

    x[1] = fSparse->GetAxis(1)->GetBinCenter(bins[1]);
    x[0] = fSparse->GetAxis(0)->GetBinCenter(bins[0]);

    if(x[0] > 19 && x[0] < 21){
      vals_in_getJetPtBin_init[0] += x[fDim];
    }
    if(x[0] > 29 && x[0] < 31){
      vals_in_getJetPtBin_init[1] += x[fDim];
    }
    if(x[0] > 39 && x[0] < 41){
      vals_in_getJetPtBin_init[2] += x[fDim];
    }
    
    
    if(x[1] > 0.09 && x[1] < 0.11){

      if(x[0] > 19 && x[0] < 21){
	vals_in_zGen_lowJetPtBin[0] += x[fDim]; 
      }
      if(x[0] > 29 && x[0] < 31){
	vals_in_zGen_centJetPtBin[0] += x[fDim];
      }
      if(x[0] > 39 && x[0] < 41){
	vals_in_zGen_highJetPtBin[0] += x[fDim];
      }
    }
    
    if(x[1] > 0.29 && x[1] < 0.31){
      if(x[0] > 19 && x[0] < 21){
	vals_in_zGen_lowJetPtBin[1] += x[fDim];
      }
      if(x[0] > 29 && x[0] < 31){
	vals_in_zGen_centJetPtBin[1] += x[fDim];
      }
      if(x[0] > 39 && x[0] < 41){
	vals_in_zGen_highJetPtBin[1] += x[fDim];
      }
    }
    
    if(x[1] > 0.49 && x[1] < 0.51){

      if(x[0] > 19 && x[0] < 21){
	vals_in_zGen_lowJetPtBin[2] += x[fDim];
      }
      if(x[0] > 29 && x[0] < 31){
	vals_in_zGen_centJetPtBin[2] += x[fDim];
      }
      if(x[0] > 39 && x[0] < 41){
	vals_in_zGen_highJetPtBin[2] += x[fDim];
      }
    }
    
    if(x[1] > 0.69 && x[1] < 0.71){
      if(x[0] > 19 && x[0] < 21){
	vals_in_zGen_lowJetPtBin[3] += x[fDim];
      }
      if(x[0] > 29 && x[0] < 31){
	vals_in_zGen_centJetPtBin[3] += x[fDim];
      }
      if(x[0] > 39 && x[0] < 41){
	vals_in_zGen_highJetPtBin[3] += x[fDim];
      }    
    }
    
    if(x[1] > 0.89 && x[1] < 0.91){
      if(x[0] > 19 && x[0] < 21){
	vals_in_zGen_lowJetPtBin[4] += x[fDim];
      }
      if(x[0] > 29 && x[0] < 31){
	vals_in_zGen_centJetPtBin[4] += x[fDim];
      }
      if(x[0] > 39 && x[0] < 41){
	vals_in_zGen_highJetPtBin[4] += x[fDim];
      }
    }
    
  }

  
  // normalize so that z gen is flat (if this is step 1) or z gen = z unfolded from the previous step

  Double_t* x2 = new Double_t[fDim + 1];
  Double_t* err2 = new Double_t[fDim + 1];
  
  Int_t *bins2 = new Int_t[fDim];
  
  double this4DBinVal = 0.;
  double this4DBinErr = 0.;

  Double_t* new_vals_in_zGen_lowJetPtBin = new Double_t[5];
  Double_t* new_vals_in_zGen_centJetPtBin = new Double_t[5];
  Double_t* new_vals_in_zGen_highJetPtBin = new Double_t[5];
  
  new_vals_in_zGen_lowJetPtBin[0] = 0.;
  new_vals_in_zGen_lowJetPtBin[1] = 0.;
  new_vals_in_zGen_lowJetPtBin[2] = 0.;
  new_vals_in_zGen_lowJetPtBin[3] = 0.;
  new_vals_in_zGen_lowJetPtBin[4] = 0.;

  new_vals_in_zGen_centJetPtBin[0] = 0.;
  new_vals_in_zGen_centJetPtBin[1] = 0.;
  new_vals_in_zGen_centJetPtBin[2] = 0.;
  new_vals_in_zGen_centJetPtBin[3] = 0.;
  new_vals_in_zGen_centJetPtBin[4] = 0.;

  new_vals_in_zGen_highJetPtBin[0] = 0.;
  new_vals_in_zGen_highJetPtBin[1] = 0.;
  new_vals_in_zGen_highJetPtBin[2] = 0.;
  new_vals_in_zGen_highJetPtBin[3] = 0.;
  new_vals_in_zGen_highJetPtBin[4] = 0.;
  
  
  for (Long64_t count2 = 0; count2 < fSparse->GetNbins(); ++count2) {
    x2[fDim] = fSparse->GetBinContent(count2, bins2);
    err2[fDim] = fSparse->GetBinError(count2);
    
    x2[1] = fSparse->GetAxis(1)->GetBinCenter(bins2[1]);
    x2[0] = fSparse->GetAxis(0)->GetBinCenter(bins2[0]);

    if(x2[1] > 0.09 && x2[1] < 0.11){

      if(x2[0] > 19 && x2[0] < 21){
	if(vals_in_zGen_lowJetPtBin[0] > 0.){
	  this4DBinVal = x2[fDim]*z_frac_lowJetPtBin[0]/(vals_in_zGen_lowJetPtBin[0]);
	  this4DBinErr = err2[fDim]*z_frac_lowJetPtBin[0]/(vals_in_zGen_lowJetPtBin[0]);
	  new_vals_in_zGen_lowJetPtBin[0] += this4DBinVal;
	}
	else{
	  this4DBinVal = 0.;
	  this4DBinErr = 0.;
	}
      }
      if(x2[0] > 29 && x2[0] < 31){
	if(vals_in_zGen_centJetPtBin[0] > 0){
	  this4DBinVal = x2[fDim]*z_frac_centJetPtBin[0]/(vals_in_zGen_centJetPtBin[0]);
	  this4DBinErr = err2[fDim]*z_frac_centJetPtBin[0]/(vals_in_zGen_centJetPtBin[0]);
	  new_vals_in_zGen_centJetPtBin[0] += this4DBinVal;
	}
	else{
	  this4DBinVal = 0.;
	  this4DBinErr = 0.;
	}
      }
      if(x2[0] > 39 && x2[0] < 41){
	if(vals_in_zGen_highJetPtBin[0] > 0.){
	  this4DBinVal = x2[fDim]*z_frac_highJetPtBin[0]/(vals_in_zGen_highJetPtBin[0]);
	  this4DBinErr = err2[fDim]*z_frac_highJetPtBin[0]/(vals_in_zGen_highJetPtBin[0]);
	  new_vals_in_zGen_highJetPtBin[0] += this4DBinVal;
	}
	else{
	  this4DBinVal = 0.;
	  this4DBinErr = 0.;
	}
      }
    }

    if(x2[1] > 0.29 && x2[1] < 0.31){

      if(x2[0] > 19 && x2[0] < 21){
	if(vals_in_zGen_lowJetPtBin[1] > 0){
	  this4DBinVal = x2[fDim]*z_frac_lowJetPtBin[1]/(vals_in_zGen_lowJetPtBin[1]);
	  this4DBinErr = err2[fDim]*z_frac_lowJetPtBin[1]/(vals_in_zGen_lowJetPtBin[1]);
	  new_vals_in_zGen_lowJetPtBin[1] += this4DBinVal;
	}
	else{
	  this4DBinVal = 0.;
	  this4DBinErr = 0.;
	}
      }
      if(x2[0] > 29 && x2[0] < 31){
	if(vals_in_zGen_centJetPtBin[1] > 0.){
	  this4DBinVal = x2[fDim]*z_frac_centJetPtBin[1]/(vals_in_zGen_centJetPtBin[1]);
	  this4DBinErr = err2[fDim]*z_frac_centJetPtBin[1]/(vals_in_zGen_centJetPtBin[1]);
	  new_vals_in_zGen_centJetPtBin[1] += this4DBinVal;
	}
	else{
	  this4DBinVal = 0.;
	  this4DBinErr = 0.;
	}
      }
      if(x2[0] > 39 && x2[0] < 41){
	if(vals_in_zGen_highJetPtBin[1] > 0.){
	  this4DBinVal = x2[fDim]*z_frac_highJetPtBin[1]/(vals_in_zGen_highJetPtBin[1]);
	  this4DBinErr = err2[fDim]*z_frac_highJetPtBin[1]/(vals_in_zGen_highJetPtBin[1]);
	  new_vals_in_zGen_highJetPtBin[1] += this4DBinVal;
	}
	else{
	  this4DBinVal = 0.;
	  this4DBinErr = 0.;
	}
      }
      
      
    }
    
    if(x2[1] > 0.49 && x2[1] < 0.51){

      if(x2[0] > 19 && x2[0] < 21){
	if(vals_in_zGen_lowJetPtBin[2] > 0.){
	  this4DBinVal = x2[fDim]*z_frac_lowJetPtBin[2]/(vals_in_zGen_lowJetPtBin[2]);
	  this4DBinErr = err2[fDim]*z_frac_lowJetPtBin[2]/(vals_in_zGen_lowJetPtBin[2]);
	  new_vals_in_zGen_lowJetPtBin[2] += this4DBinVal;
	}
	else{
	  this4DBinVal = 0.;
	  this4DBinErr = 0.;
	}
      }
      if(x2[0] > 29 && x2[0] < 31){
	if(vals_in_zGen_centJetPtBin[2] > 0.){
	  this4DBinVal = x2[fDim]*z_frac_centJetPtBin[2]/(vals_in_zGen_centJetPtBin[2]);
	  this4DBinErr = err2[fDim]*z_frac_centJetPtBin[2]/(vals_in_zGen_centJetPtBin[2]);
	  new_vals_in_zGen_centJetPtBin[2] += this4DBinVal;
	}
	else{
	  this4DBinVal = 0.;
	  this4DBinErr = 0.;
	}
      }
      if(x2[0] > 39 && x2[0] < 41){
	if(vals_in_zGen_highJetPtBin[2] > 0.){
	  this4DBinVal = x2[fDim]*z_frac_highJetPtBin[2]/(vals_in_zGen_highJetPtBin[2]);
	  this4DBinErr = err2[fDim]*z_frac_highJetPtBin[2]/(vals_in_zGen_highJetPtBin[2]);
	  new_vals_in_zGen_highJetPtBin[2] += this4DBinVal;
	}
	else{
	  this4DBinVal = 0.;
	  this4DBinErr = 0.;
	}
      }
      
      
    }
    
    if(x2[1] > 0.69 && x2[1] < 0.71){

      if(x2[0] > 19 && x2[0] < 21){
	if(vals_in_zGen_lowJetPtBin[3] > 0.){
	  this4DBinVal = x2[fDim]*z_frac_lowJetPtBin[3]/(vals_in_zGen_lowJetPtBin[3]);
	  this4DBinErr = err2[fDim]*z_frac_lowJetPtBin[3]/(vals_in_zGen_lowJetPtBin[3]);
	  new_vals_in_zGen_lowJetPtBin[3] += this4DBinVal;
	}
	else{
	  this4DBinVal = 0.;
	  this4DBinErr = 0.;
	}
      }
      if(x2[0] > 29 && x2[0] < 31){
	if(vals_in_zGen_centJetPtBin[3] > 0.){
	  this4DBinVal = x2[fDim]*z_frac_centJetPtBin[3]/(vals_in_zGen_centJetPtBin[3]);
	  this4DBinErr = err2[fDim]*z_frac_centJetPtBin[3]/(vals_in_zGen_centJetPtBin[3]);
	  new_vals_in_zGen_centJetPtBin[3] += this4DBinVal;
	}
	else{
	  this4DBinVal = 0.;
	  this4DBinErr = 0.;
	}
      }
      if(x2[0] > 39 && x2[0] < 41){
	if(vals_in_zGen_highJetPtBin[3] > 0.){
	  this4DBinVal = x2[fDim]*z_frac_highJetPtBin[3]/(vals_in_zGen_highJetPtBin[3]);
	  this4DBinErr = err2[fDim]*z_frac_highJetPtBin[3]/(vals_in_zGen_highJetPtBin[3]);
	  new_vals_in_zGen_highJetPtBin[3] += this4DBinVal;
	}
	else{
	  this4DBinVal = 0.;
	  this4DBinErr = 0.;
	}
      }
      
    }
    
    if(x2[1] > 0.89 && x2[1] < 0.91){

      if(x2[0] > 19 && x2[0] < 21){
	if(vals_in_zGen_lowJetPtBin[4] > 0.){
	  this4DBinVal = x2[fDim]*z_frac_lowJetPtBin[4]/(vals_in_zGen_lowJetPtBin[4]);
	  this4DBinErr = err2[fDim]*z_frac_lowJetPtBin[4]/(vals_in_zGen_lowJetPtBin[4]);
	  new_vals_in_zGen_lowJetPtBin[4] += this4DBinVal;
	}
	else{
	  this4DBinVal = 0.;
	  this4DBinErr = 0.;
	}
      }
      if(x2[0] > 29 && x2[0] < 31){
	if(vals_in_zGen_centJetPtBin[4] > 0.){
	  this4DBinVal = x2[fDim]*z_frac_centJetPtBin[4]/(vals_in_zGen_centJetPtBin[4]);
	  this4DBinErr = err2[fDim]*z_frac_centJetPtBin[4]/(vals_in_zGen_centJetPtBin[4]);
	  new_vals_in_zGen_centJetPtBin[4] += this4DBinVal;
	}
	else{
	  this4DBinVal = 0.;
	  this4DBinErr = 0.;
	}
      }
      if(x2[0] > 39 && x2[0] < 41){
	if(vals_in_zGen_highJetPtBin[4] > 0){
	  this4DBinVal = x2[fDim]*z_frac_highJetPtBin[4]/(vals_in_zGen_highJetPtBin[4]);
	  this4DBinErr = err2[fDim]*z_frac_highJetPtBin[4]/(vals_in_zGen_highJetPtBin[4]);
	  new_vals_in_zGen_highJetPtBin[4] += this4DBinVal;
	}
	else{
	  this4DBinVal = 0.;
	  this4DBinErr = 0.;
	}
      }
            
    }
     
    //set bin content for normalized 4D histo
    fSparse_norm->GetBin(bins2);
    fSparse_norm->SetBinContent(count2, this4DBinVal);
    fSparse_norm->SetBinError(count2, this4DBinErr);
  }


  //get the sum of gen jet pt in each bin of gen jet pt 
  
  Double_t* x3 = new Double_t[fDim + 1];
  Int_t *bins3 = new Int_t[fDim];
   
  Double_t* vals_in_jetPtGen_all = new Double_t[3];
  vals_in_jetPtGen_all[0] = 0.;
  vals_in_jetPtGen_all[1] = 0.;
  vals_in_jetPtGen_all[2] = 0.;

  for (Long64_t count3 = 0; count3 < fSparse_norm->GetNbins(); ++count3) {
    x3[fDim] = fSparse_norm->GetBinContent(count3, bins3);

    x3[0] = fSparse_norm->GetAxis(0)->GetBinCenter(bins3[0]);

    if(x3[0] > 19 && x3[0] < 21){
      vals_in_jetPtGen_all[0] += x3[fDim];
    }
    if(x3[0] > 29 && x3[0] < 31){
      vals_in_jetPtGen_all[1] += x3[fDim];
    }
    if(x3[0] > 39 && x3[0] < 41){
      vals_in_jetPtGen_all[2] += x3[fDim];
    }
  }

  // normalize gen jet pt, so it corresponds to gen jet pt

  Double_t* x4 = new Double_t[fDim + 1];
  Double_t* err4 = new Double_t[fDim + 1];

  Int_t *bins4 = new Int_t[fDim];

  double this4DBinVal2 = 0.;
  double this4DBinErr2 = 0.;

  Double_t* new_vals_in_zGen_jetPt = new Double_t[5];
  new_vals_in_zGen_jetPt[0] = 0.;
  new_vals_in_zGen_jetPt[1] = 0.;
  new_vals_in_zGen_jetPt[2] = 0.;

  for (Long64_t count4 = 0; count4 < fSparse_norm->GetNbins(); ++count4) {

    x4[fDim] = fSparse_norm->GetBinContent(count4, bins4);
    err4[fDim] = fSparse_norm->GetBinError(count4);

    x4[0] = fSparse_norm->GetAxis(0)->GetBinCenter(bins4[0]);

    if(x4[0] > 19 && x4[0] < 21){
      this4DBinVal2 = x4[fDim]*vals_in_getJetPtBin_init[0]/vals_in_jetPtGen_all[0];
      this4DBinErr2 = err4[fDim]*vals_in_getJetPtBin_init[0]/vals_in_jetPtGen_all[0];
      new_vals_in_zGen_jetPt[0] += this4DBinVal2;
    }

    if(x4[0] > 29 && x4[0] < 31){
      this4DBinVal2 = x4[fDim]*vals_in_getJetPtBin_init[1]/vals_in_jetPtGen_all[1];
      this4DBinErr2 = err4[fDim]*vals_in_getJetPtBin_init[1]/vals_in_jetPtGen_all[1];
      new_vals_in_zGen_jetPt[1] += this4DBinVal2;
    }

    if(x4[0] > 39 && x4[0] < 41){
      this4DBinVal2 = x4[fDim]*vals_in_getJetPtBin_init[2]/vals_in_jetPtGen_all[2];
      this4DBinErr2 = err4[fDim]*vals_in_getJetPtBin_init[2]/vals_in_jetPtGen_all[2];
      new_vals_in_zGen_jetPt[2] += this4DBinVal2;
    }

    fSparse_norm_genJetPt->GetBin(bins4);
    fSparse_norm_genJetPt->SetBinContent(count4, this4DBinVal2);
    fSparse_norm_genJetPt->SetBinError(count4, this4DBinErr2);
  }

  TFile *file_forProf_varBin = new TFile(outputfile.c_str(),"RECREATE");

  h_jetpPtGen_jetPtReco->Write();
  
  fSparse->Write();
  fSparse_norm->Write();
  fSparse_norm_genJetPt->Write();
  
  h_z_gen->Write();

  hGenZJetPtCentBin->Write();
  hGenZJetPtLowBin->Write();
  hGenZJetPtHighBin->Write();
  
  hRespZJetPtCentBin->Write();
  hRespZJetPtLowBin->Write();
  hRespZJetPtHighBin->Write();

  file_forProf_varBin->Write();
  file_forProf_varBin->Close();
  
}

void prepareInputsData(Int_t step = 1){

  prepare(false,true,step);
  prepare(false,false,step);

  prepare(true,true,step);
  prepare(true,false,step);
  

}
