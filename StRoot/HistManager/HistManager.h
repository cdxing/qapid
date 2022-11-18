#ifndef HistManager_h
#define HistManager_h

#include "TVector2.h"
#include "TVector3.h"
#include "TString.h"
#include "StMessMgr.h"

class TProfile2D;
class TProfile;
class TH1F;
class TH1D;
class TH2D;
class StPicoTrack;
class StPicoDst;

class HistManager
{
  public:
    HistManager();
    ~HistManager();


    // QAPID
    void InitQAPID();
    void WriteQAPID();
    // Event level QA
    void FillEventQA(TVector3 PrimaryVertex, Double_t RefMult, Double_t TofMult, Double_t TrackMult); 
    void FillEventQaCut(TVector3 PrimaryVertex, Double_t RefMult, Double_t TofMult, Double_t TrackMult); 
    void FillEventCent(Int_t centrality); // 
    void FillEventCut(Int_t CutID); // 
    // Track level QA
    void FillTrackQA(StPicoTrack *PicoTrack, TVector3 pVtx); // 
    void FillTrackPhysics(StPicoTrack *track);
    void FillTrackTof(StPicoDst* pico, StPicoTrack *track);
    void FillTrackCut(Int_t CutID); // 
    // PID 
    void FillTriton(StPicoDst *pico, StPicoTrack *PicoTrack, Double_t y_mid);  
    void FillDeuteron(StPicoDst *pico, StPicoTrack *PicoTrack, Double_t y_mid);  
    void FillProtonETof(StPicoTrack *track, Double_t y_mid); // ETof test
    void FillProton(StPicoDst *pico, StPicoTrack *PicoTrack, Double_t y_mid);  
    void FillKaon(StPicoDst *pico, StPicoTrack *PicoTrack, Double_t y_mid); 
    void FillPion(StPicoDst *pico, StPicoTrack *PicoTrack, Double_t y_mid); 
    void FillPIDMult(Int_t, Int_t, Int_t, Int_t, Int_t); 

    // Scalor Product method

  private:

    //Event level QA Plots
    TH1D*                h_eventCheck;
    TH1D*                h_centralities;
    TH1D*                h_zvtx;
    TH1D*                h_trackmult;
    TH1D*                h_refmult;
    TH1D*                h_tofmult;
    TH2D*                h2_trans_vtx;
    TH2D*                h2_trans_vtx_cut;
    TH2D*                h2_refmult_vs_trackmult;
    TH2D*                h2_tofmult_vs_trackmult;
    TH2D*                h2_tofmult_vs_refmult;

    //Track level QA Plots
    TH1D*                h_trackCheck;
    TH1D*                h_nhits_dEdx;
    TH1D*                h_nhits;
    TH1D*                h_nhits_ratio;
    TH1D*                h_DCA;
    TH1D*                h_pT;
    TH1D*                h_eta;
    TH1D*                h_phi;
    TH2D*                h2_dEdx_vs_qp;
    TH2D*                h2_dEdx_vs_qp_half;
    TH2D*                h2_beta_vs_qp;
    TH2D*                h2_m2_vs_qp;

    //PID Plots
    TH1D*                h_mult_pp;
    TH1D*                h_mult_pm;
    TH1D*                h_mult_kp;
    TH1D*                h_mult_km;
    TH1D*                h_mult_pr;
    TH1D*                h_mult_pbar;
    //TH1D*                h_mult_de;
    //TH1D*                h_mult_tr;

    TH1D*                h_pT_pp;
    TH1D*                h_pT_pm;
    TH1D*                h_pT_kp;
    TH1D*                h_pT_km;
    TH1D*                h_pT_pr;
    TH1D*                h_pT_pbar;
    TH1D*                h_pT_de;
    TH1D*                h_pT_tr;

    TH1D*                h_dndy_pp;
    TH1D*                h_dndy_pm;
    TH1D*                h_dndy_kp;
    TH1D*                h_dndy_km;
    TH1D*                h_dndy_pr;
    TH1D*                h_dndy_pbar;
    TH1D*                h_dndy_de;
    TH1D*                h_dndy_tr;

    TH1D*                h_eta_pp;
    TH1D*                h_eta_pm;
    TH1D*                h_eta_kp;
    TH1D*                h_eta_km;
    TH1D*                h_eta_pr;
    TH1D*                h_eta_pbar;
    TH1D*                h_eta_de;
    TH1D*                h_eta_tr;

    TH1D*                h_phi_pp;
    TH1D*                h_phi_pm;
    TH1D*                h_phi_kp;
    TH1D*                h_phi_km;
    TH1D*                h_phi_pr;
    TH1D*                h_phi_pbar;
    TH1D*                h_phi_de;
    TH1D*                h_phi_tr;

    TH2D*                h2_pT_vs_yCM_pp;
    TH2D*                h2_pT_vs_yCM_pm;
    TH2D*                h2_pT_vs_yCM_kp;
    TH2D*                h2_pT_vs_yCM_km;
    TH2D*                h2_pT_vs_yCM_pr;
    TH2D*                h2_pT_vs_yCM_pr_ETof;
    TH2D*                h2_pT_vs_yCM_pbar;
    TH2D*                h2_pT_vs_yCM_pbar_ETof;
    TH2D*                h2_pT_vs_yCM_de;
    TH2D*                h2_pT_vs_yCM_tr;

    TH2D*                h2_dEdx_vs_qp_pp;
    TH2D*                h2_dEdx_vs_qp_pm;
    TH2D*                h2_dEdx_vs_qp_kp;
    TH2D*                h2_dEdx_vs_qp_km;
    TH2D*                h2_dEdx_vs_qp_pr;
    TH2D*                h2_dEdx_vs_qp_pbar;
    TH2D*                h2_dEdx_vs_qp_de;
    TH2D*                h2_dEdx_vs_qp_tr;

    TH2D*                h2_beta_vs_qp_pp;
    TH2D*                h2_beta_vs_qp_pm;
    TH2D*                h2_beta_vs_qp_kp;
    TH2D*                h2_beta_vs_qp_km;
    TH2D*                h2_beta_vs_qp_pr;
    TH2D*                h2_beta_vs_qp_pbar;
    TH2D*                h2_beta_vs_qp_de;
    TH2D*                h2_beta_vs_qp_tr;

    TH2D*                h2_m2_vs_qp_pp;
    TH2D*                h2_m2_vs_qp_pm;
    TH2D*                h2_m2_vs_qp_kp;
    TH2D*                h2_m2_vs_qp_km;
    TH2D*                h2_m2_vs_qp_pr;
    TH2D*                h2_m2_vs_qp_pbar;
    TH2D*                h2_m2_vs_qp_de;
    TH2D*                h2_m2_vs_qp_tr;


    static TString mVStr[2];

    ClassDef(HistManager,1)
};

#endif
