#include "ConstManager.h"
#include "TMath.h"
#include "StarClassLibrary/SystemOfUnits.h"

ClassImp(ConstManager)

// event cut
Int_t ConstManager::Ncentralities = 16;
Int_t ConstManager::CENT_BINS = 16;
Int_t ConstManager::FIRST_CENT =  16 - ConstManager::CENT_BINS;
std::map<float,float> ConstManager::mVzMaxMap = ConstManager::createVzMaxMap();
Float_t ConstManager::mVrMax = 2.0;
Float_t ConstManager::mVzMax = 70.0;
Int_t ConstManager::mMatchedToFMin = 2;

// track cut
std::map<float,float> ConstManager::mSigScaleMap = ConstManager::createSigScaleMap();
Float_t ConstManager::mDcaEPMax[5] = {3.0,1.0,1.0,1.0,1.0}; // for event plane reconstruction: 3.0 for 200GeV, 1.0 for BES
Float_t ConstManager::mDcaTrMax = 1.0; // for pion, kaon, proton mDcaTrMax = 1.0 for flow;
Float_t ConstManager::mDcaTrMax_phi = 2.0; // for phi meson mDcaTrMax = 2.0 to fill a tree and apply an additional cut
Int_t ConstManager::mHitsDedxMin = 5;
Int_t ConstManager::mHitsFitTPCMin = 15;
Int_t ConstManager::mHitsMaxTPCMin = 0;
Float_t ConstManager::mHitsRatioTPCMin = 0.52;
Float_t ConstManager::mEtaMax = 1.0;
Float_t ConstManager::mPrimPtMin[5] = {0.15,0.2,0.2,0.2,0.2}; // for event plane reconstruction and for pion, kaon, proton: 0.15 for 200 GeV, 0.2 for BES
Float_t ConstManager::mGlobPtMin = 0.1; // for phi, Lambda, K0s
Float_t ConstManager::mPrimPtMax = 2.0;
Float_t ConstManager::mPrimPtWeight = 2.0;
Float_t ConstManager::mPrimMomMax = 10.0; // also use for gMom
Float_t ConstManager::mMass2Min = -10.0;
Double_t ConstManager::MAGFIELDFACTOR = kilogauss;
Int_t ConstManager::mTrackMin = 2;
Int_t ConstManager::mTrackMin_Full = 4;
Float_t ConstManager::mToFYLocalMax = 1.8;
Float_t ConstManager::mToFZLocalMax = 1.8;
Float_t ConstManager::mNSigmaElectronMax = 2.5;
Float_t ConstManager::mNSigmaPionMax = 2.5;
Float_t ConstManager::mNSigmaKaonMax = 2.5;
Float_t ConstManager::mNSigmaProtonMax = 2.5;
Float_t ConstManager::mMassPion = 0.13957;
Float_t ConstManager::mMassKaon = 0.49368;
Float_t ConstManager::mMassProton = 0.93827;
Float_t ConstManager::D_M0_DE = 1.875613;   // Deuteron
Float_t ConstManager::D_M0_TR = 2.808921;   // Triton
Float_t ConstManager::D_M0_HE3 = 2.809414;   // Helium-3
Float_t ConstManager::D_M0_AL = 3.727379;   // Alpha
Float_t ConstManager::mVzVpdDiffMax = 3.0;
Float_t ConstManager::mEta_Gap[4] = {0.05,0.10,0.20,0.50};

// used constant
Float_t ConstManager::mShiftOrder2[5] = {2.0, 4.0, 6.0, 8.0, 10.0};
Float_t ConstManager::mShiftOrder3[5] = {3.0, 6.0, 9.0, 12.0, 15.0};
std::map<float,float> ConstManager::mExScaleMap = ConstManager::createExScaleMap();

//   change pT bin from 11 to 14  line 52,53,58,59,61,62,81 changed
// pt bin
//                              0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 10 , 11 , 12 , 13 , 14
Float_t ConstManager::pt_low[14] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8};
Float_t ConstManager::pt_up[14]  = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0};
Float_t ConstManager::pt_low_spectra[32] = {0.15,0.25,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.2,3.4,3.8};
Float_t ConstManager::pt_up_spectra[32]  = {0.25,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.2,3.4,3.8,4.2};
// x and y range
//                              0 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 ,  8 ,  9 ,  10 
Float_t ConstManager::x_low[14] = {-0.4,-0.6,-0.6,-0.6,-0.6,-0.6,-0.6,-0.8,-0.8,-1.0,-1.0,-1.0,-1.0,-1.0};
Float_t ConstManager::x_up[14]  = { 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.6, 1.6, 2.0, 2.0, 2.0, 2.0, 2.4};

Float_t ConstManager::y_low[14] = {-0.15,-0.2,-0.3,-0.45,-0.6,-0.8,-1.0,-1.3,-1.5,-1.7,-1.8,-1.8,-1.8,-1.8};
Float_t ConstManager::y_up[14]  = { 0.15, 0.2, 0.2, 0.20, 0.3, 0.4, 0.6, 1.0, 1.0, 1.0, 1.4, 1.4, 1.4, 1.4};

// Centrality bin
Int_t ConstManager::cent_low[4] = {0,7,4,0}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
Int_t ConstManager::cent_up[4]  = {8,8,6,3}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
TString ConstManager::Centrality_01[4] = {"0080","0010","1040","4080"};
TString ConstManager::Centrality_23[4] = {"0070","0010","1040","4070"};

// phi-Psi bin
Double_t ConstManager::phi_Psi2_low[7] = {0.0,TMath::Pi()/14.0,2.0*TMath::Pi()/14.0,3.0*TMath::Pi()/14.0,4.0*TMath::Pi()/14.0,5.0*TMath::Pi()/14.0,6.0*TMath::Pi()/14.0};
Double_t ConstManager::phi_Psi2_up[7]  = {TMath::Pi()/14.0,2.0*TMath::Pi()/14.0,3.0*TMath::Pi()/14.0,4.0*TMath::Pi()/14.0,5.0*TMath::Pi()/14.0,6.0*TMath::Pi()/14.0,7.0*TMath::Pi()/14.0};
Double_t ConstManager::phi_Psi3_low[7] = {0.0,TMath::Pi()/21.0,2.0*TMath::Pi()/21.0,3.0*TMath::Pi()/21.0,4.0*TMath::Pi()/21.0,5.0*TMath::Pi()/21.0,6.0*TMath::Pi()/21.0};
Double_t ConstManager::phi_Psi3_up[7]  = {TMath::Pi()/21.0,2.0*TMath::Pi()/21.0,3.0*TMath::Pi()/21.0,4.0*TMath::Pi()/21.0,5.0*TMath::Pi()/21.0,6.0*TMath::Pi()/21.0,7.0*TMath::Pi()/21.0};

Double_t ConstManager::Psi2_low[3] = {-3.0*TMath::Pi()/2.0,-1.0*TMath::Pi()/2.0,1.0*TMath::Pi()/2.0};
Double_t ConstManager::Psi2_up[3]  = {-1.0*TMath::Pi()/2.0, 1.0*TMath::Pi()/2.0,3.0*TMath::Pi()/2.0};
Double_t ConstManager::Psi3_low[5] = {-4.0*TMath::Pi()/3.0,-3.0*TMath::Pi()/3.0,-1.0*TMath::Pi()/3.0,1.0*TMath::Pi()/3.0,3.0*TMath::Pi()/3.0};
Double_t ConstManager::Psi3_up[5]  = {-3.0*TMath::Pi()/3.0,-1.0*TMath::Pi()/3.0, 1.0*TMath::Pi()/3.0,3.0*TMath::Pi()/3.0,4.0*TMath::Pi()/3.0};

Int_t ConstManager::pt_total = 14;
Int_t ConstManager::pt_total_spectra = 32;
Int_t ConstManager::Centrality_total = 4;
Int_t ConstManager::Centrality_start = 0;
Int_t ConstManager::Centrality_stop  = 4;

Int_t ConstManager::EtaGap_total = 4;
Int_t ConstManager::EtaGap_start = 0;
Int_t ConstManager::EtaGap_stop  = 1;

Int_t ConstManager::Charge_total = 2;  // shaowei
Int_t ConstManager::Charge_start = 0;
Int_t ConstManager::Charge_stop  = 2;  // shaowei


TString ConstManager::Energy[5] = {"200GeV","54GeV","39GeV","7GeV","15GeV"};

// mix event
Int_t ConstManager::Bin_Centrality = 9;
Int_t ConstManager::Bin_VertexZ = 10;
Int_t ConstManager::Bin_Phi_Psi = 5;
Int_t ConstManager::Buffer_depth = 5;
Int_t ConstManager::Flag_ME = 0; // 0 for same event, 1 for mixed event
TString ConstManager::MixEvent[2] = {"SE","ME"};
