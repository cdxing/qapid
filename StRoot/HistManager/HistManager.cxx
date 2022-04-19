#include "HistManager.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TString.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "../StPicoEvent/StPicoTrack.h"
#include "../StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "../ConstManager/ConstManager.h"

ClassImp(HistManager)

    TString HistManager::mVStr[2] = {"pos","neg"};

//---------------------------------------------------------------------------------

HistManager::HistManager()
{
}

//---------------------------------------------------------------------------------

HistManager::~HistManager()
{
}


//----------------------------------------------------------------------------

void HistManager::InitQAPID()
{
  // QA plots
  h_eventCheck = new TH1D("h_eventCheck","Event number after each cut;;Events", 6, 0, 6);
  h_eventCheck->GetXaxis()->SetBinLabel(1,"No Cuts");
  h_eventCheck->GetXaxis()->SetBinLabel(2,"Event Cut");
  h_eventCheck->GetXaxis()->SetBinLabel(3,"Centrality Cut");
  /*h_eventCheck->GetXaxis()->SetBinLabel(2,"Bad Run Cut");
  h_eventCheck->GetXaxis()->SetBinLabel(3,"Minbias Cut");
  h_eventCheck->GetXaxis()->SetBinLabel(4,"V_{Z} Cut");
  h_eventCheck->GetXaxis()->SetBinLabel(5,"V_{r} Cut");
  h_eventCheck->GetXaxis()->SetBinLabel(6,"Centrality Cut");*/

  h_trackCheck = new TH1D("h_trackCheck","Track number after each cut;;Tracks", 3, 0, 3);
  h_trackCheck->GetXaxis()->SetBinLabel(1,"Event Cuts Only");
  h_trackCheck->GetXaxis()->SetBinLabel(2,"QA Cuts");
  h_trackCheck->GetXaxis()->SetBinLabel(3,"PID Cuts");

  h_centralities = new TH1D("h_centralities", "Centralities;Centrality ID;Events", ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS);
  h_centralities->GetXaxis()->SetTitle("Centrality bin");
  h_centralities->GetYaxis()->SetTitle("# of events");

  //h_zvtx = new TH1D("h_zvtx","Primary Vertex Position in z;Distance (cm);Events", 100, 190, 210);
  h_zvtx = new TH1D("h_zvtx","Primary Vertex Position in z;Distance (cm);Events", 840, -210, 210);

  h2_trans_vtx = new TH2D("h2_trans_vtx","Primary Vertex after V_{z} Cut;x (cm);y (cm)", 500, -5, 5, 500, -5, 5);
  h2_trans_vtx_cut = new TH2D("h2_trans_vtx_cut","Final Primary Vertices;x (cm);y (cm)", 500, -5, 5, 500, -5, 5);

  h_refmult = new TH1D("h_refmult","Reference multiplicity",1001,-0.5,1000.5);
  h_refmult->GetXaxis()->SetTitle("RefMult");
  h_refmult->GetXaxis()->SetTitle("# of events");

  h_tofmult = new TH1D("h_tofmult","TOF multiplicity",1001,-0.5,1000.5);
  h_tofmult->GetXaxis()->SetTitle("TofMult");
  h_tofmult->GetXaxis()->SetTitle("# of events");

  h_trackmult = new TH1D("h_trackmult","Actual track multiplicity",1501,-0.5,1500.5);
  h_trackmult->GetXaxis()->SetTitle("TrackMult");
  h_trackmult->GetXaxis()->SetTitle("# of events");

  h2_refmult_vs_trackmult = new TH2D("h2_refmult_vs_trackmult","RefMult vs. Actual track multiplicity;TrackMult;RefMult",1501,-0.5,1500.5,1001,-0.5,1000.5);

  h2_tofmult_vs_trackmult = new TH2D("h2_tofmult_vs_trackmult","TofMult vs. Actual track multiplicity;TrackMult;TofMult",1501,-0.5,1500.5,1001,-0.5,1000.5);

  h2_tofmult_vs_refmult = new TH2D("h2_tofmult_vs_refmult","TofMult vs. RefMult;RefMult;TofMult",1001,-0.5,1000.5,1001,-0.5,1000.5);

  h_nhits       = new TH1D("h_nhits", "nHits;Number of hits;Tracks", 50, 0, 50);
  h_nhits_dEdx  = new TH1D("h_nhits_dEdx","nHitsdEdx;Number of hits;Tracks", 50, 0, 50);
  h_nhits_ratio = new TH1D("h_nhits_ratio","nhitsFit/nhitsPoss;nhitsFit/nhitsPoss;Tracks",200,0.0,2.0);
  h_DCA         = new TH1D("h_DCA","DCA (cm);DCA (cm);Tracks",100,0.0,10.0);
  
  h_pT = new TH1D("h_pT","p_{T};p_{T};Tracks",1000,0.0,10.0);
  h_eta = new TH1D("h_eta","#eta;#eta;Tracks",600,-6.0,6.0);
  h_phi = new TH1D("h_phi","#phi (Radian);#phi;Tracks",1000,-1.5*TMath::Pi(),1.5*TMath::Pi());
  h2_dEdx_vs_qp = new TH2D("h2_dEdx_vs_qp", "dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
  h2_dEdx_vs_qp_half = new TH2D("h2_dEdx_vs_qp_half", "dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, 0, 4, 500, 0, 12);
  h2_beta_vs_qp = new TH2D("h2_beta_vs_qp","1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  h2_m2_vs_qp = new TH2D("h2_m2_vs_qp", "m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)", 400, -4, 4, 400, -0.1, 1.5);

  h_pT_pp = new TH1D("h_pT_pp","#pi^{+} p_{T};p_{T} (GeV/c);Tracks",1000,0.0,5.0);
  h_pT_pm = new TH1D("h_pT_pm","#pi^{-} p_{T}; p_{T} (GeV/c);Tracks",1000,0.0,5.0);
  h_pT_kp = new TH1D("h_pT_kp","K^{+} p_{T}; p_{T} (GeV/c);Tracks",1000,0.0,5.0);
  h_pT_km = new TH1D("h_pT_km","K^{-} p_{T}; p_{T} (GeV/c);Tracks",1000,0.0,5.0);
  h_pT_pr = new TH1D("h_pT_pr","Proton p_{T};p_{T} (GeV/c);Tracks",1000,0.0,5.0);
  //h_pT_de = new TH1D("h_pT_de","Deuteron p_{T};p_{T} (GeV/c);Tracks",1000,0.0,5.0);
  //h_pT_tr = new TH1D("h_pT_tr","Triton p_{T};p_{T} (GeV/c);Tracks",1000,0.0,5.0);

  h_eta_pp = new TH1D("h_eta_pp","#pi^{+} #eta;#eta;Tracks",500,-5.0,5.0);
  h_eta_pm = new TH1D("h_eta_pm","#pi^{-} #eta;#eta;Tracks",500,-5.0,5.0);
  h_eta_kp = new TH1D("h_eta_kp","K^{+} #eta;#eta;Tracks",500,-5.0,5.0);
  h_eta_km = new TH1D("h_eta_km","K^{-} #eta;#eta;Tracks",500,-5.0,5.0);
  h_eta_pr = new TH1D("h_eta_pr","Proton #eta;#eta;Tracks",500,-5.0,5.0);
  //h_eta_de = new TH1D("h_eta_de","Deuteron #eta;#eta;Tracks",500,-5.0,5.0);
  //h_eta_tr = new TH1D("h_eta_tr","Triton #eta;#eta;Tracks",500,-5.0,5.0);

  h_dndy_pp = new TH1D("h_dndy_pp", "#pi^{+} Raw Rapidity Spectrum;y;dN/dy", 80, -2, 2);
  h_dndy_pm = new TH1D("h_dndy_pm", "#pi^{-} Raw Rapidity Spectrum;y;dN/dy", 80, -2, 2);
  h_dndy_kp = new TH1D("h_dndy_kp", "K^{+} Raw Rapidity Spectrum;y;dN/dy",   80, -2, 2);
  h_dndy_km = new TH1D("h_dndy_km", "K^{-} Raw Rapidity Spectrum;y;dN/dy",   80, -2, 2);
  h_dndy_pr = new TH1D("h_dndy_pr", "Proton Raw Rapidity Spectrum;y;dN/dy",  80, -2, 2);
  //h_dndy_de = new TH1D("h_dndy_de", "Deuteron Raw Rapidity Spectrum;y;dN/dy",  80, -2, 2);
  //h_dndy_tr = new TH1D("h_dndy_tr", "Triton Raw Rapidity Spectrum;y;dN/dy",  80, -2, 2);

  h_phi_pp = new TH1D("h_phi_pp","#pi^{+} #phi (Radian);#phi;Tracks",1000,-1.5*TMath::Pi(),1.5*TMath::Pi());
  h_phi_pm = new TH1D("h_phi_pm","#pi^{-} #phi (Radian);#phi;Tracks",1000,-1.5*TMath::Pi(),1.5*TMath::Pi());
  h_phi_kp = new TH1D("h_phi_kp","K^{+} #phi (Radian);#phi;Tracks",1000,-1.5*TMath::Pi(),1.5*TMath::Pi());
  h_phi_km = new TH1D("h_phi_km","K^{-} #phi (Radian);#phi;Tracks",1000,-1.5*TMath::Pi(),1.5*TMath::Pi());
  h_phi_pr = new TH1D("h_phi_pr","Proton #phi (Radian);#phi;Tracks",1000,-1.5*TMath::Pi(),1.5*TMath::Pi());
  //h_phi_de = new TH1D("h_phi_de","Deuteron #phi (Radian);#phi;Tracks",1000,-1.5*TMath::Pi(),1.5*TMath::Pi());
  //h_phi_tr = new TH1D("h_phi_tr","Triton #phi (Radian);#phi;Tracks",1000,-1.5*TMath::Pi(),1.5*TMath::Pi());

  h_mult_pp = new TH1D("h_mult_pp","#pi^{+} track multiplicity;#pi^{+} Mult;Events",1001,-0.5,1000.5);
  h_mult_pm = new TH1D("h_mult_pm","#pi^{-} track multiplicity;#pi^{-} Mult;Events",1001,-0.5,1000.5);
  h_mult_kp = new TH1D("h_mult_kp","K^{#plus} track multiplicity;K^{+} Mult;Events",1001,-0.5,1000.5);
  h_mult_km = new TH1D("h_mult_km","K^{-} track multiplicity;K^{-} Mult;Events",1001,-0.5,1000.5);
  h_mult_pr = new TH1D("h_mult_pr","Proton track multiplicity;Proton Mult;Events",1001,-0.5,1000.5);
  //h_mult_de = new TH1D("h_mult_de","Deuteron track multiplicity;Deuteron Mult;Events",1001,-0.5,1000.5);
  //h_mult_tr = new TH1D("h_mult_tr","Triton track multiplicity;Triton Mult;Events",1001,-0.5,1000.5);
  h2_dEdx_vs_qp_pp = new TH2D("h2_dEdx_vs_qp_pp", "#pi^{+} dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
  h2_dEdx_vs_qp_pm = new TH2D("h2_dEdx_vs_qp_pm", "#pi^{-} dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
  h2_dEdx_vs_qp_kp = new TH2D("h2_dEdx_vs_qp_kp", "K^{+} dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
  h2_dEdx_vs_qp_km = new TH2D("h2_dEdx_vs_qp_km", "K^{-} dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
  h2_dEdx_vs_qp_pr = new TH2D("h2_dEdx_vs_qp_pr", "Proton dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
  //h2_dEdx_vs_qp_de = new TH2D("h2_dEdx_vs_qp_de", "Deuteron dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
  //h2_dEdx_vs_qp_tr = new TH2D("h2_dEdx_vs_qp_tr", "Triton dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);

  h2_beta_vs_qp_pp = new TH2D("h2_beta_vs_qp_pp","#pi^{+} 1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  h2_beta_vs_qp_pm = new TH2D("h2_beta_vs_qp_pm","#pi^{-} 1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  h2_beta_vs_qp_kp = new TH2D("h2_beta_vs_qp_kp","K^{+} 1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  h2_beta_vs_qp_km = new TH2D("h2_beta_vs_qp_km","K^{-} 1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  h2_beta_vs_qp_pr = new TH2D("h2_beta_vs_qp_pr","Proton 1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  //h2_beta_vs_qp_de = new TH2D("h2_beta_vs_qp_de","Deuteron 1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
  //h2_beta_vs_qp_tr = new TH2D("h2_beta_vs_qp_tr","Triton 1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);

  h2_m2_vs_qp_pp = new TH2D("h2_m2_vs_qp_pp", "#pi^{+} m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)", 400, -4, 4, 400, -0.1, 1.5);
  h2_m2_vs_qp_pm = new TH2D("h2_m2_vs_qp_pm", "#pi^{-} m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)", 400, -4, 4, 400, -0.1, 1.5);
  h2_m2_vs_qp_kp = new TH2D("h2_m2_vs_qp_kp", "K^{+} m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)",   400, -4, 4, 400, -0.1, 1.5);
  h2_m2_vs_qp_km = new TH2D("h2_m2_vs_qp_km", "K^{-} m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)",   400, -4, 4, 400, -0.1, 1.5);
  h2_m2_vs_qp_pr = new TH2D("h2_m2_vs_qp_pr", "Proton m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)",  400, -4, 4, 400, -0.1, 1.5);
  //h2_m2_vs_qp_de = new TH2D("h2_m2_vs_qp_de", "Deuteron m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)",  400, -4, 4, 400, -0.1, 1.5);
  //h2_m2_vs_qp_tr = new TH2D("h2_m2_vs_qp_tr", "Triton m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)",  400, -4, 4, 400, -0.1, 1.5);

  h2_pT_vs_yCM_pp = new TH2D("h2_pT_vs_yCM_pp", "#pi^{+};y-y_{mid};p_{T} (GeV/c)", 300, -1.2, 1.2, 300, 0, 3.0);
  h2_pT_vs_yCM_pm = new TH2D("h2_pT_vs_yCM_pm", "#pi^{-};y-y_{mid};p_{T} (GeV/c)", 300, -1.2, 1.2, 300, 0, 3.0);
  h2_pT_vs_yCM_kp = new TH2D("h2_pT_vs_yCM_kp", "K^{+};y-y_{mid};p_{T} (GeV/c)",   300, -1.2, 1.2, 300, 0, 3.0);
  h2_pT_vs_yCM_km = new TH2D("h2_pT_vs_yCM_km", "K^{-};y-y_{mid};p_{T} (GeV/c)",   300, -1.2, 1.2, 300, 0, 3.0);
  h2_pT_vs_yCM_pr = new TH2D("h2_pT_vs_yCM_pr", "Proton;y-y_{mid};p_{T} (GeV/c)",  300, -1.2, 1.2, 300, 0, 3.0);
  //h2_pT_vs_yCM_de = new TH2D("h2_pT_vs_yCM_de", "Deuteron;y-y_{mid};p_{T} (GeV/c)",  300, -1.2, 1.2, 300, 0, 3.0);
  //h2_pT_vs_yCM_tr = new TH2D("h2_pT_vs_yCM_tr", "Triton;y-y_{mid};p_{T} (GeV/c)",  300, -1.2, 1.2, 300, 0, 3.0);

}


//----------------------------------------------------------------------------
// Event QA plots
void HistManager::FillEventQA(TVector3 PrimaryVertex, Double_t RefMult, Double_t TofMult, Double_t TrackMult)
{
    Float_t x_vtx = PrimaryVertex.X();
    Float_t y_vtx = PrimaryVertex.Y();
    Float_t z_vtx = PrimaryVertex.Z();
    
    h_zvtx->Fill(z_vtx);
    h2_trans_vtx->Fill(x_vtx,y_vtx);
    h_trackmult->Fill(TrackMult);
    h_refmult->Fill(RefMult);
    h_tofmult->Fill(TofMult);
    h2_refmult_vs_trackmult->Fill(TrackMult,RefMult);
    h2_tofmult_vs_trackmult->Fill(TrackMult,TofMult);
    h2_tofmult_vs_refmult->Fill(RefMult,TofMult);
}

void HistManager::FillEventQaCut(TVector3 PrimaryVertex, Double_t RefMult, Double_t TofMult, Double_t TrackMult)
{
    Float_t x_vtx = PrimaryVertex.X();
    Float_t y_vtx = PrimaryVertex.Y();
    Float_t z_vtx = PrimaryVertex.Z();
    
    h2_trans_vtx_cut->Fill(x_vtx,y_vtx);
}

void HistManager::FillEventCent(Int_t centrality) 
{
	h_centralities->Fill(centrality);
}

void HistManager::FillEventCut(Int_t CutID)
{
	h_eventCheck->Fill(CutID);
}
// Track QA plots
void HistManager::FillTrackQA(StPicoTrack *track, TVector3 pVtx)
{
  h_nhits->Fill(track->nHits());
  h_nhits_dEdx->Fill(track->nHitsDedx());
  h_nhits_ratio->Fill((Double_t)track->nHitsFit()/track->nHitsPoss());
  h_DCA->Fill(track->gDCA(pVtx).Mag());
}
// Track physics plots
void HistManager::FillTrackPhysics(StPicoTrack *track)
{
  Double_t d_pT = track->pPt();
  //Double_t d_px = track->pMom().X();
  //Double_t d_py = track->pMom().Y();
  //Double_t d_pz = track->pMom().Z();
  Double_t d_mom = track->pMom().Mag();
  Double_t d_eta = track->pMom().Eta();
  Double_t d_phi = track->pMom().Phi();
  Short_t  s_charge = track->charge();
  h_phi->Fill(d_phi);
  h_eta->Fill(d_eta);
  h_pT->Fill(d_pT);
  h2_dEdx_vs_qp->Fill(s_charge * d_mom, track->dEdx());
}

// Track physics plots
void HistManager::FillTrackTof(StPicoDst *pico, StPicoTrack *track)
{
  Double_t d_tofBeta = -999.0;
  Double_t d_m2 = 0.0;
  Short_t  s_charge = track->charge();
  Double_t d_mom = track->pMom().Mag();
  Int_t trackTofIndex = track->bTofPidTraitsIndex();
  if(trackTofIndex >= 0)
    d_tofBeta = pico->btofPidTraits(trackTofIndex)->btofBeta();

  if(d_tofBeta != -999.0)
    {
      d_m2 = d_mom*d_mom*( (1.0 / (d_tofBeta*d_tofBeta)) - 1.0 );
      h2_beta_vs_qp->Fill(s_charge*d_mom, 1.0/d_tofBeta);
      h2_m2_vs_qp->Fill(s_charge*d_mom, d_m2);
    }
  else 
    {
    	std::cout << "no tof info" << std::endl;
    }
}
void HistManager::FillTrackCut(Int_t CutID)
{
	h_trackCheck->Fill(CutID);
}
//----------------------------------------------------------------------------
//PID
void HistManager::FillProton(StPicoDst *pico, StPicoTrack *track, Double_t y_mid)
{
  Double_t d_pT = track->pPt();
  Double_t d_px = track->pMom().X();
  Double_t d_py = track->pMom().Y();
  Double_t d_pz = track->pMom().Z();
  Double_t d_eta = track->pMom().Eta();
  Double_t d_phi = track->pMom().Phi();
  Double_t d_mom = track->pMom().Mag();
  Double_t d_dEdx = track->dEdx();
  Short_t  s_charge = track->charge();
  TLorentzVector ltrack;
  ltrack.SetXYZM(d_px,d_py,d_pz,ConstManager::mMassProton);
  Double_t mRapidity = ltrack.Rapidity();

  Double_t d_m2 = -999.0;
  Double_t d_tofBeta = -999.0;
  Int_t trackTofIndex = track->bTofPidTraitsIndex();
  if(trackTofIndex >= 0)
    d_tofBeta = pico->btofPidTraits(trackTofIndex)->btofBeta();
  if(d_tofBeta != -999.0)
    {
      d_m2 = d_mom*d_mom*( (1.0 / (d_tofBeta*d_tofBeta)) - 1.0 );
    }
  
  h_eta_pr->Fill(d_eta);
  h_phi_pr->Fill(d_phi);
  h_pT_pr->Fill(d_pT);
  h_dndy_pr->Fill(mRapidity);
  h2_pT_vs_yCM_pr->Fill(mRapidity - y_mid, d_pT);
  h2_dEdx_vs_qp_pr->Fill(s_charge*d_mom, track->dEdx());
  h2_beta_vs_qp_pr->Fill(s_charge*d_mom, 1.0/d_tofBeta);
  h2_m2_vs_qp_pr->Fill(s_charge*d_mom, d_m2);
}

void HistManager::FillKaon(StPicoDst *pico, StPicoTrack *track, Double_t y_mid)
{
  Double_t d_pT = track->pPt();
  Double_t d_px = track->pMom().X();
  Double_t d_py = track->pMom().Y();
  Double_t d_pz = track->pMom().Z();
  Double_t d_eta = track->pMom().Eta();
  Double_t d_phi = track->pMom().Phi();
  Double_t d_mom = track->pMom().Mag();
  Double_t d_dEdx = track->dEdx();
  Short_t  s_charge = track->charge();
  TLorentzVector ltrack;
  ltrack.SetXYZM(d_px,d_py,d_pz,ConstManager::mMassKaon);
  Double_t mRapidity = ltrack.Rapidity();

  Double_t d_m2 = -999.0;
  Double_t d_tofBeta = -999.0;
  Int_t trackTofIndex = track->bTofPidTraitsIndex();
  if(trackTofIndex >= 0)
    d_tofBeta = pico->btofPidTraits(trackTofIndex)->btofBeta();
  if(d_tofBeta != -999.0)
    {
      d_m2 = d_mom*d_mom*( (1.0 / (d_tofBeta*d_tofBeta)) - 1.0 );
    }
  if(s_charge > 0)
  {
    h_eta_kp->Fill(d_eta);
    h_phi_kp->Fill(d_phi);
    h_pT_kp->Fill(d_pT);
    h_dndy_kp->Fill(mRapidity);
    h2_pT_vs_yCM_kp->Fill(mRapidity - y_mid, d_pT);
    h2_dEdx_vs_qp_kp->Fill(s_charge*d_mom, track->dEdx());
    h2_beta_vs_qp_kp->Fill(s_charge*d_mom, 1.0/d_tofBeta);
    h2_m2_vs_qp_kp->Fill(s_charge*d_mom, d_m2);
  }
  if(s_charge < 0)
  {
    h_eta_km->Fill(d_eta);
    h_phi_km->Fill(d_phi);
    h_pT_km->Fill(d_pT);
    h_dndy_km->Fill(mRapidity);
    h2_pT_vs_yCM_km->Fill(mRapidity - y_mid, d_pT);
    h2_dEdx_vs_qp_km->Fill(s_charge*d_mom, track->dEdx());
    h2_beta_vs_qp_km->Fill(s_charge*d_mom, 1.0/d_tofBeta);
    h2_m2_vs_qp_km->Fill(s_charge*d_mom, d_m2);
  }
}
void HistManager::FillPion(StPicoDst *pico, StPicoTrack *track, Double_t y_mid)
{
  Double_t d_pT = track->pPt();
  Double_t d_px = track->pMom().X();
  Double_t d_py = track->pMom().Y();
  Double_t d_pz = track->pMom().Z();
  Double_t d_eta = track->pMom().Eta();
  Double_t d_phi = track->pMom().Phi();
  Double_t d_mom = track->pMom().Mag();
  Double_t d_dEdx = track->dEdx();
  Short_t  s_charge = track->charge();
  TLorentzVector ltrack;
  ltrack.SetXYZM(d_px,d_py,d_pz,ConstManager::mMassPion);
  Double_t mRapidity = ltrack.Rapidity();

  Double_t d_m2 = -999.0;
  Double_t d_tofBeta = -999.0;
  Int_t trackTofIndex = track->bTofPidTraitsIndex();
  if(trackTofIndex >= 0)
    d_tofBeta = pico->btofPidTraits(trackTofIndex)->btofBeta();
  if(d_tofBeta != -999.0)
    {
      d_m2 = d_mom*d_mom*( (1.0 / (d_tofBeta*d_tofBeta)) - 1.0 );
    }
  if(s_charge > 0)
  {
    h_eta_pp->Fill(d_eta);
    h_phi_pp->Fill(d_phi);
    h_pT_pp->Fill(d_pT);
    h_dndy_pp->Fill(mRapidity);
    h2_pT_vs_yCM_pp->Fill(mRapidity - y_mid, d_pT);
    h2_dEdx_vs_qp_pp->Fill(s_charge*d_mom, track->dEdx());
    h2_beta_vs_qp_pp->Fill(s_charge*d_mom, 1.0/d_tofBeta);
    h2_m2_vs_qp_pp->Fill(s_charge*d_mom, d_m2);
  }
  if(s_charge < 0)
  {
    h_eta_pm->Fill(d_eta);
    h_phi_pm->Fill(d_phi);
    h_pT_pm->Fill(d_pT);
    h_dndy_pm->Fill(mRapidity);
    h2_pT_vs_yCM_pm->Fill(mRapidity - y_mid, d_pT);
    h2_dEdx_vs_qp_pm->Fill(s_charge*d_mom, track->dEdx());
    h2_beta_vs_qp_pm->Fill(s_charge*d_mom, 1.0/d_tofBeta);
    h2_m2_vs_qp_pm->Fill(s_charge*d_mom, d_m2);
  }
}
void HistManager::FillPIDMult(Int_t N_pp, Int_t N_pm, Int_t N_kp, Int_t N_km, Int_t N_pr) 
{
  h_mult_pp->Fill(N_pp);
  h_mult_pm->Fill(N_pm);
  h_mult_kp->Fill(N_kp);
  h_mult_km->Fill(N_km);
  h_mult_pr->Fill(N_pr);
}
//----------------------------------------------------------------------------

void HistManager::WriteQAPID()
{
    h_eventCheck->Write();
    h_trackCheck->Write();
    h_centralities->Write();
    h_zvtx->Write();
    h_trackmult->Write();
    h_refmult->Write();
    h_tofmult->Write();
    h2_trans_vtx->Write();
    h2_trans_vtx_cut->Write();
    h2_refmult_vs_trackmult->Write();
    h2_tofmult_vs_trackmult->Write();
    h2_tofmult_vs_refmult->Write();

    h_nhits->Write();
    h_nhits_dEdx->Write();
    h_nhits_ratio->Write();
    h_DCA->Write();
    
    h_phi->Write();
    h_eta->Write();
    h_pT->Write();
    h2_dEdx_vs_qp->Write();
    
    h2_beta_vs_qp->Write();                                               
    h2_m2_vs_qp->Write();

    h_mult_pp->Write();
    h_mult_pm->Write();
    h_mult_kp->Write();
    h_mult_km->Write();
    h_mult_pr->Write();
    h_eta_pr->Write();
    h_phi_pr->Write();
    h_pT_pr->Write();    
    h_dndy_pr->Write();
    h2_pT_vs_yCM_pr->Write();
    h2_dEdx_vs_qp_pr->Write();
    h2_beta_vs_qp_pr->Write();
    h2_m2_vs_qp_pr->Write();

    h_eta_kp->Write();
    h_phi_kp->Write();
    h_pT_kp->Write();
    h_dndy_kp->Write();
    h2_pT_vs_yCM_kp->Write();
    h2_dEdx_vs_qp_kp->Write();
    h2_beta_vs_qp_kp->Write();
    h2_m2_vs_qp_kp->Write();

    h_eta_km->Write();
    h_phi_km->Write();
    h_pT_km->Write();
    h_dndy_km->Write();
    h2_pT_vs_yCM_km->Write();
    h2_dEdx_vs_qp_km->Write();
    h2_beta_vs_qp_km->Write();
    h2_m2_vs_qp_km->Write();
    
    h_eta_pp->Write();
    h_phi_pp->Write();
    h_pT_pp->Write();
    h_dndy_pp->Write();
    h2_pT_vs_yCM_pp->Write();
    h2_dEdx_vs_qp_pp->Write();
    h2_beta_vs_qp_pp->Write();
    h2_m2_vs_qp_pp->Write();

    h_eta_pm->Write();
    h_phi_pm->Write();
    h_pT_pm->Write();
    h_dndy_pm->Write();
    h2_pT_vs_yCM_pm->Write();
    h2_dEdx_vs_qp_pm->Write();
    h2_beta_vs_qp_pm->Write();
    h2_m2_vs_qp_pm->Write();
}

//----------------------------------------------------------------------------

