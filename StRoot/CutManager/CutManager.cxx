#include "CutManager.h"
#include "StRoot/ConstManager/ConstManager.h"
#include "StRoot/StPicoEvent/StPicoDst.h"  // shaowei  18c
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"

// Bichsel header
#include "StRoot/StBichsel/Bichsel.h"

#include "StMessMgr.h"

ClassImp(CutManager)

// Bichsel Function
Double_t bichselZ(Double_t *x,Double_t *par) 
{
  Double_t pove   = x[0];
  Double_t poverm = pove/par[0];
  return TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(poverm),par[1]));
}
    //StRefMultCorr* CutManager::mRefMultCorr = NULL; // shaowei
    //---------------------------------------------------------------------------------

CutManager::CutManager(ConfigReader configs)
{
    mConfigs = configs;
}

//---------------------------------------------------------------------------------

CutManager::~CutManager()
{
}

//---------------------------------------------------------------------------------

bool CutManager::isGoodTrigger(StPicoDst *pico)
{
    //Bool_t b_good_trig = false;
    Bool_t b_good_trig = true; // disable trigger cut
    StPicoEvent *event = pico->event();
    std::vector<UInt_t> triggerIDs = event->triggerIds();
    for (UInt_t i = 0; i < triggerIDs.size(); i++)
	   { if ( mConfigs.triggersMatch(triggerIDs[i]) ) {b_good_trig = true;} }

    /*for(auto trg: triggers)
    {
        if(event->isTrigger(trg)) return kTRUE;
    }*/
    if(b_good_trig) return kTRUE;
}

bool CutManager::passEventCut(StPicoDst *pico)
{
    //if(!CutManager::isGoodTrigger(pico)) return kFALSE;
    //No Tirgger Cut for now.

    //std::cout<< "event cut test 0 "<<std::endl;
    StPicoEvent *event = pico->event();
    if(!event)
    {
        return kFALSE;
    }
    //std::cout<< "event cut test 1 "<<std::endl;
    //std::cout<< "config vertex z low "<< mConfigs.z_vtx_low <<std::endl;

    // initialize StRefMultCorr
    const Float_t vx = event->primaryVertex().X();
    const Float_t vy = event->primaryVertex().Y();
    const Float_t vz = event->primaryVertex().Z();
    // event vertex cut
    // vz cut
    if(vz < mConfigs.z_vtx_low || vz > mConfigs.z_vtx_high)
    {
        return kFALSE;
    }
    //std::cout<< "event cut test 2 "<<std::endl;
    //if(event->btofTrayMultiplicity()<2)return kFALSE;
    // vr cut
    if(mConfigs.fixed_target==1){ // FXT
    	if(sqrt(vx*vx+(vy+2)*(vy+2)) > mConfigs.r_vtx)
    	{
    	    return kFALSE;
    	}
    } else 
    if(mConfigs.fixed_target==0){ // COL
    	if(sqrt(vx*vx+vy*vy) > mConfigs.r_vtx)
    	{
    	    return kFALSE;
    	}
	    
    }
    //std::cout<< "event cut test 3 "<<std::endl;
    // vz-vzVpd cut for 200 GeV

    return kTRUE;
}

//---------------------------------------------------------------------------------
bool CutManager::isTofTrack(StPicoDst *pico, StPicoTrack *track)
{
  Bool_t b_tofTrack = false;
  Double_t d_tofBeta = -999.0;
  //Double_t d_m2 = 0.0;
  Int_t trackTofIndex = track->bTofPidTraitsIndex();
  if(trackTofIndex >= 0)
    d_tofBeta = pico->btofPidTraits(trackTofIndex)->btofBeta();

  if(d_tofBeta != -999.0)
    {
      b_tofTrack = true;
      /*d_m2 = d_mom*d_mom*( (1.0 / (d_tofBeta*d_tofBeta)) - 1.0 );
      h2_beta_vs_qp->Fill(s_charge*d_mom, 1.0/d_tofBeta);
      h2_m2_vs_qp->Fill(s_charge*d_mom, d_m2);*/
    }
  return b_tofTrack;
}

bool CutManager::isETofTrack(StPicoDst *pico, StPicoTrack *track)
{
  Bool_t b_etofTrack = false;
  //Double_t d_etofBeta = -999.0;
  //Double_t d_m2 = 0.0;
  Int_t trackTofIndex = track->eTofPidTraitsIndex();
  if(trackTofIndex >= 0) b_etofTrack = true;
    //d_etofBeta = pico->betofPidTraits(trackTofIndex)->betofBeta();

  //if(d_etofBeta != -999.0)
    //{
      //b_etofTrack = true;
      ///*d_m2 = d_mom*d_mom*( (1.0 / (d_etofBeta*d_etofBeta)) - 1.0 );
      //h2_beta_vs_qp->Fill(s_charge*d_mom, 1.0/d_etofBeta);
      //h2_m2_vs_qp->Fill(s_charge*d_mom, d_m2);*/
    //}
  return b_etofTrack;
}

bool CutManager::passTrackBasic(StPicoTrack *track)
{
    // nHitsFit cut
    if(track->nHitsFit() < mConfigs.nHits)     //  
    {
        return kFALSE;
    }

    // nHitsRatio cut
    if(track->nHitsMax() <=  mConfigs.nHits )   // 
    {
        return kFALSE;
    }
    if((Float_t)track->nHitsFit()/(Float_t)track->nHitsMax() < mConfigs.nHits_ratio)  // 0.52
    {
        return kFALSE;
    }

    // eta cut
    /*Float_t eta = track->pMom().PseudoRapidity();
    if(fabs(eta) > ConstManager::mEtaMax)  // 1
    {
        return kFALSE;
    }*/

    return kTRUE;
}


bool CutManager::passTrackEP(StPicoTrack *track, float dca)
{
    if(!track) return kFALSE;
    if(!passTrackBasic(track)) return kFALSE;

    // dca cut for event plane reconstruction: 200GeV = 3.0, BES = 1.0
    //if(dca > ConstManager::mDcaEPMax[mEnergy])   // change by shaowei  1
    if(dca > 1.0)
    {
        return kFALSE;
    }

    // pt cut 0.2 - 2.0 GeV/c
    Float_t pt = track->pMom().Perp();
    Float_t p  = track->pMom().Mag();
    if(!(pt > ConstManager::mPrimPtMin[mEnergy] && pt < ConstManager::mPrimPtMax && p < ConstManager::mPrimMomMax)) // 0.2<pt<2  p<10
    {
        return kFALSE;
    }

    return kTRUE;
}
//---------------------------------------------------------------------------------
// PID
bool CutManager::isTriton(StPicoDst *pico, StPicoTrack *track)
{
  //=========================================================
  //          Bichsel Function Setup
  //=========================================================
  Double_t log2dx = 1.0;
  Double_t xStart = 0.01;
  Double_t xStop  = 3.0;
  Int_t npx = 10000;
  //                      Mass  log2(dx)
  Double_t params[2] = {  1.0,   log2dx  };

  params[0] = ConstManager::D_M0_TR;
  TF1 *bichselZ_tr = new TF1(Form("BichselZ_tr_log2dx_%i",(int)log2dx),bichselZ,xStart,xStop,2);
  if (!bichselZ_tr) { std::cout << "De function error" << std::endl; return 1; }
  bichselZ_tr->SetParameters(params); 
  bichselZ_tr->SetNpx(npx);

  Bool_t triton   = false;
  Double_t d_tofBeta = -999.0;
  Double_t d_m2 = -999.0;
  Double_t d_mom = track->pMom().Mag();
  Short_t  s_charge = track->charge();
  Float_t d_dEdx = track->dEdx();
  Double_t d_zTriton = (s_charge == 1) ? TMath::Log(d_dEdx / bichselZ_tr->Eval(d_mom)) : -999.0;
  triton = (d_zTriton > mConfigs.z_tr_low) && (d_zTriton < mConfigs.z_tr_high);
  Int_t trackTofIndex = track->bTofPidTraitsIndex();
  if(trackTofIndex >= 0)
    d_tofBeta = pico->btofPidTraits(trackTofIndex)->btofBeta();
  if(d_tofBeta != -999.0)
    {
      d_m2 = d_mom*d_mom*( (1.0 / (d_tofBeta*d_tofBeta)) - 1.0 );
    }

//  triton = ((d_zTriton > mConfigs.z_tr_low) && (d_zTriton < mConfigs.z_tr_high) ) ||
//  	 ( (d_tofBeta != -999.0)&&
//		   (d_zTriton > mConfigs.z_tr_low) &&
//		    (d_zTriton < mConfigs.z_tr_high) &&
//		    (d_m2 > mConfigs.m2_tr_low) &&
//		    (d_m2 < mConfigs.m2_tr_high)
//	 );
		  // TRITON
		  if (d_mom >= 1.0 && d_mom < 4.0)
		    {
		      if (d_mom >= 1.0 && d_mom < 1.1 && d_zTriton > -0.332011 && d_zTriton < 0.251103) triton = true;
		      else if (d_mom >= 1.1 && d_mom < 1.2 && d_zTriton > -0.310412 && d_zTriton < 0.296090) triton = true;
		      else if (d_mom >= 1.2 && d_mom < 1.3 && d_zTriton > -0.293322 && d_zTriton < 0.334467) triton = true;
		      else if (d_mom >= 1.3 && d_mom < 1.4 && d_zTriton > -0.270550 && d_zTriton < 0.373857) triton = true;
		      else if (d_mom >= 1.4 && d_mom < 1.5 && d_zTriton > -0.248412 && d_zTriton < 0.406237) triton = true;
		      else if (d_mom >= 1.5 && d_mom < 1.6 && d_zTriton > -0.228044 && d_zTriton < 0.333261) triton = true;
		      else if (d_mom >= 1.6 && d_mom < 1.7 && d_zTriton > -0.210093 && d_zTriton < 0.343588) triton = true;
		      else if (d_mom >= 1.7 && d_mom < 1.8 && d_zTriton > -0.190900 && d_zTriton < 0.332586) triton = true;
		      else if (d_mom >= 1.8 && d_mom < 1.9 && d_zTriton > -0.183153 && d_zTriton < 0.334197) triton = true;
		      else if (d_mom >= 1.9 && d_mom < 2.0 && d_zTriton > -0.166020 && d_zTriton < 0.323303) triton = true;
		      else if (d_mom >= 2.0 && d_mom < 2.1 && d_zTriton > -0.102334 && d_zTriton < 0.307724) triton = true;
		      else if (d_mom >= 2.1 && d_mom < 2.2 && d_zTriton > -0.091053 && d_zTriton < 0.294345) triton = true;
		      else if (d_mom >= 2.2 && d_mom < 2.3 && d_zTriton > -0.076457 && d_zTriton < 0.285978) triton = true;
		      else if (d_mom >= 2.3 && d_mom < 2.4 && d_zTriton > -0.055669 && d_zTriton < 0.253769) triton = true;
		      else if (d_mom >= 2.4 && d_mom < 2.5 && d_zTriton > -0.035848 && d_zTriton < 0.254487) triton = true;
		      else if (d_mom >= 2.5 && d_mom < 2.6 && d_zTriton > -0.027266 && d_zTriton < 0.249350) triton = true;
		      else if (d_mom >= 2.6 && d_mom < 2.7 && d_zTriton > -0.028152 && d_zTriton < 0.236713) triton = true;
		      else if (d_mom >= 2.7 && d_mom < 2.8 && d_zTriton > -0.027867 && d_zTriton < 0.227672) triton = true;
		      else if (d_mom >= 2.8 && d_mom < 2.9 && d_zTriton > -0.024675 && d_zTriton < 0.222215) triton = true;
		      else if (d_mom >= 2.9 && d_mom < 3.0 && d_zTriton > -0.019179 && d_zTriton < 0.227362) triton = true;
		      else if (d_mom >= 3.0 && d_mom < 3.1 && d_zTriton > -0.013267 && d_zTriton < 0.236052) triton = true;
		      else if (d_mom >= 3.1 && d_mom < 3.2 && d_zTriton > -0.007851 && d_zTriton < 0.246071) triton = true;
		      else if (d_mom >= 3.2 && d_mom < 3.3 && d_zTriton > -0.006311 && d_zTriton < 0.254907) triton = true;
		      else if (d_mom >= 3.3 && d_mom < 3.4 && d_zTriton > 0.019834 && d_zTriton < 0.244291) triton = true;
		      else if (d_mom >= 3.4 && d_mom < 3.5 && d_zTriton > 0.031221 && d_zTriton < 0.273652) triton = true;
		      else if (d_mom >= 3.5 && d_mom < 3.6 && d_zTriton > 0.068248 && d_zTriton < 0.257484) triton = true;
		      else if (d_mom >= 3.6 && d_mom < 3.7 && d_zTriton > 0.088804 && d_zTriton < 0.260799) triton = true;
		      else if (d_mom >= 3.7 && d_mom < 3.8 && d_zTriton > 0.091490 && d_zTriton < 0.271776) triton = true;
		      else if (d_mom >= 3.8 && d_mom < 3.9 && d_zTriton > 0.106161 && d_zTriton < 0.285652) triton = true;
		      else if (d_mom >= 3.9 && d_mom < 4.0 && d_zTriton > 0.103653 && d_zTriton < 0.299234) triton = true;
		    }
		  //else if (tofTrack)
		  else if (d_tofBeta != -999.0)
		    {
		      if (d_zTriton > mConfigs.z_tr_low &&
			  d_zTriton < mConfigs.z_tr_high &&
			  d_m2 > mConfigs.m2_tr_low &&
			  d_m2 < mConfigs.m2_tr_high)
			triton = true;
		    }

  return triton;
}

bool CutManager::isDeuteron(StPicoDst *pico, StPicoTrack *track)
{
  //=========================================================
  //          Bichsel Function Setup
  //=========================================================
  Double_t log2dx = 1.0;
  Double_t xStart = 0.01;
  Double_t xStop  = 3.0;
  Int_t npx = 10000;
  //                      Mass  log2(dx)
  Double_t params[2] = {  1.0,   log2dx  };

  params[0] = ConstManager::D_M0_DE;
  TF1 *bichselZ_de = new TF1(Form("BichselZ_de_log2dx_%i",(int)log2dx),bichselZ,xStart,xStop,2);
  if (!bichselZ_de) { std::cout << "De function error" << std::endl; return 1; }
  bichselZ_de->SetParameters(params); 
  bichselZ_de->SetNpx(npx);

  Bool_t deuteron   = false;
  Double_t d_tofBeta = -999.0;
  Double_t d_m2 = -999.0;
  Double_t d_mom = track->pMom().Mag();
  Short_t  s_charge = track->charge();
  Float_t d_dEdx = track->dEdx();
  Double_t d_zDeuteron = (s_charge == 1) ? TMath::Log(d_dEdx / bichselZ_de->Eval(d_mom)) : -999.0;
  deuteron = (d_zDeuteron > mConfigs.z_de_low) && (d_zDeuteron < mConfigs.z_de_high);
  Int_t trackTofIndex = track->bTofPidTraitsIndex();
  if(trackTofIndex >= 0)
    d_tofBeta = pico->btofPidTraits(trackTofIndex)->btofBeta();
  if(d_tofBeta != -999.0)
    {
      d_m2 = d_mom*d_mom*( (1.0 / (d_tofBeta*d_tofBeta)) - 1.0 );
    }

//  deuteron = ((d_zDeuteron > mConfigs.z_de_low) && (d_zDeuteron < mConfigs.z_de_high) ) ||
//  	 ( (d_tofBeta != -999.0)&&
//		   (d_zDeuteron > mConfigs.z_de_low) &&
//		    (d_zDeuteron < mConfigs.z_de_high) &&
//		    (d_m2 > mConfigs.m2_de_low) &&
//		    (d_m2 < mConfigs.m2_de_high)
//	 );

		  if (d_mom >= 0.4 && d_mom < 3.0)
		    {
		      if (d_mom >= 0.4 && d_mom < 0.5 && d_zDeuteron > -0.476112 && d_zDeuteron < 0.248539) deuteron = true;
		      else if (d_mom >= 0.5 && d_mom < 0.6 && d_zDeuteron > -0.445644 && d_zDeuteron < 0.311067) deuteron = true;
		      else if (d_mom >= 0.6 && d_mom < 0.7 && d_zDeuteron > -0.43008 && d_zDeuteron < 0.331624) deuteron = true;
		      else if (d_mom >= 0.7 && d_mom < 0.8 && d_zDeuteron > -0.416061 && d_zDeuteron < 0.341399) deuteron = true;
		      else if (d_mom >= 0.8 && d_mom < 0.9 && d_zDeuteron > -0.404842 && d_zDeuteron < 0.338091) deuteron = true;
		      else if (d_mom >= 0.9 && d_mom < 1.0 && d_zDeuteron > -0.37419 && d_zDeuteron < 0.337724) deuteron = true;
		      else if (d_mom >= 1.0 && d_mom < 1.1 && d_zDeuteron > -0.32986 && d_zDeuteron < 0.332241) deuteron = true;
		      else if (d_mom >= 1.1 && d_mom < 1.2 && d_zDeuteron > -0.332995 && d_zDeuteron < 0.325582) deuteron = true;
		      else if (d_mom >= 1.2 && d_mom < 1.3 && d_zDeuteron > -0.306145 && d_zDeuteron < 0.319532) deuteron = true;
		      else if (d_mom >= 1.3 && d_mom < 1.4 && d_zDeuteron > -0.275987 && d_zDeuteron < 0.313227) deuteron = true;
		      else if (d_mom >= 1.4 && d_mom < 1.5 && d_zDeuteron > -0.250464 && d_zDeuteron < 0.301911) deuteron = true;
		      else if (d_mom >= 1.5 && d_mom < 1.6 && d_zDeuteron > -0.215135 && d_zDeuteron < 0.302149) deuteron = true;
		      else if (d_mom >= 1.6 && d_mom < 1.7 && d_zDeuteron > -0.176733 && d_zDeuteron < 0.308644) deuteron = true;
		      else if (d_mom >= 1.7 && d_mom < 1.8 && d_zDeuteron > -0.160866 && d_zDeuteron < 0.29673) deuteron = true;
		      else if (d_mom >= 1.8 && d_mom < 1.9 && d_zDeuteron > -0.149249 && d_zDeuteron < 0.281362) deuteron = true;
		      else if (d_mom >= 1.9 && d_mom < 2.0 && d_zDeuteron > -0.0830817 && d_zDeuteron < 0.273483) deuteron = true;
		      else if (d_mom >= 2.0 && d_mom < 2.1 && d_zDeuteron > -0.065219 && d_zDeuteron < 0.269654) deuteron = true;
		      else if (d_mom >= 2.1 && d_mom < 2.2 && d_zDeuteron > -0.04952 && d_zDeuteron < 0.265074) deuteron = true;
		      else if (d_mom >= 2.2 && d_mom < 2.3 && d_zDeuteron > -0.0358834 && d_zDeuteron < 0.258749) deuteron = true;
		      else if (d_mom >= 2.3 && d_mom < 2.4 && d_zDeuteron > -0.0218641 && d_zDeuteron < 0.25294) deuteron = true;
		      else if (d_mom >= 2.4 && d_mom < 2.5 && d_zDeuteron > -0.0114193 && d_zDeuteron < 0.244108) deuteron = true;
		      else if (d_mom >= 2.5 && d_mom < 2.6 && d_zDeuteron > -0.000659632 && d_zDeuteron < 0.205416) deuteron = true;
		      else if (d_mom >= 2.6 && d_mom < 2.7 && d_zDeuteron > 0.010662  && d_zDeuteron < 0.198006) deuteron = true;
		      else if (d_mom >= 2.7 && d_mom < 2.8 && d_zDeuteron > 0.0203815 && d_zDeuteron < 0.189092) deuteron = true;
		      else if (d_mom >= 2.8 && d_mom < 2.9 && d_zDeuteron > 0.0313737 && d_zDeuteron < 0.181285) deuteron = true;
		      else if (d_mom >= 2.9 && d_mom < 3.0 && d_zDeuteron > 0.0446902 && d_zDeuteron < 0.174561) deuteron = true;
		    }
		  //else if (tofTrack)
		  else if (d_tofBeta != -999.0)
		    {
		      if (d_zDeuteron > mConfigs.z_de_low &&
			  d_zDeuteron < mConfigs.z_de_high &&
			  d_m2 > mConfigs.m2_de_low &&
			  d_m2 < mConfigs.m2_de_high)
			deuteron = true;
		    }
  return deuteron;
}

bool CutManager::isProton(StPicoDst *pico, StPicoTrack *track)
{
  Double_t d_TPCnSigmaProton = track->nSigmaProton();
  Bool_t proton   = false;
  Double_t d_tofBeta = -999.0;
  Double_t d_m2 = -999.0;
  Double_t d_mom = track->pMom().Mag();
  Short_t  s_charge = track->charge();
  Int_t trackTofIndex = track->bTofPidTraitsIndex();
  if(trackTofIndex >= 0)
    d_tofBeta = pico->btofPidTraits(trackTofIndex)->btofBeta();
  if(d_tofBeta != -999.0)
    {
      d_m2 = d_mom*d_mom*( (1.0 / (d_tofBeta*d_tofBeta)) - 1.0 );
    }

  //proton = ((d_TPCnSigmaProton > mConfigs.nSig_pr_low) && (d_TPCnSigmaProton < mConfigs.nSig_pr_high) && (s_charge > 0) && d_mom < 1.3) ||
  proton = ((d_TPCnSigmaProton > mConfigs.nSig_pr_low) && (d_TPCnSigmaProton < mConfigs.nSig_pr_high) && d_mom < 1.3) ||
  	 ((d_TPCnSigmaProton > mConfigs.nSig_pr_low) &&
         (d_TPCnSigmaProton < mConfigs.nSig_pr_high) &&
         (d_m2 > 0.8) &&
         (d_m2 < 1.2))
	  ;
  return proton;
}

bool CutManager::isKaon(StPicoDst *pico, StPicoTrack *track)
{
  Double_t d_TPCnSigmaKaon   = track->nSigmaKaon();
  //Short_t  s_charge = track->charge();
  Bool_t kaon   = false;
  Double_t d_tofBeta = -999.0;
  Double_t d_m2 = -999.0;
  Double_t d_mom = track->pMom().Mag();
  Int_t trackTofIndex = track->bTofPidTraitsIndex();
  if(trackTofIndex >= 0)
    d_tofBeta = pico->btofPidTraits(trackTofIndex)->btofBeta();
  if(d_tofBeta != -999.0)
    {
      d_m2 = d_mom*d_mom*( (1.0 / (d_tofBeta*d_tofBeta)) - 1.0 );
    }
  kaon = (d_TPCnSigmaKaon > mConfigs.nSig_ka_low) &&
         (d_TPCnSigmaKaon < mConfigs.nSig_ka_high) &&
         (d_m2 > mConfigs.m2_ka_low) &&
         (d_m2 < mConfigs.m2_ka_high);

  return kaon;
}

bool CutManager::isPion(StPicoDst *pico, StPicoTrack *track)
{
  Double_t d_TPCnSigmaPion   = track->nSigmaPion();
  //Short_t  s_charge = track->charge();
  Bool_t pion   = false;
  Double_t d_tofBeta = -999.0;
  Double_t d_m2 = -999.0;
  Double_t d_mom = track->pMom().Mag();
  Int_t trackTofIndex = track->bTofPidTraitsIndex();
  if(trackTofIndex >= 0)
    d_tofBeta = pico->btofPidTraits(trackTofIndex)->btofBeta();
  if(d_tofBeta != -999.0)
    {
      d_m2 = d_mom*d_mom*( (1.0 / (d_tofBeta*d_tofBeta)) - 1.0 );
    }
  pion = (d_TPCnSigmaPion > mConfigs.nSig_pi_low) &&
         (d_TPCnSigmaPion < mConfigs.nSig_pi_high) &&
         (d_m2 > mConfigs.m2_pi_low) &&
         (d_m2 < mConfigs.m2_pi_high);

  return pion;
}
//---------------------------------------------------------------------------------

Int_t CutManager::getCentrality(int gRefMult)
{
    int centrality=-99;
    	int centHigh[16];
    	int centLow[16] ;
    if (mConfigs.sqrt_s_NN == 3.0)
    {
	//std::cout << "3GeV centrality" << std::endl;
    	int centHigh_3p0GeV[16]={6,8,11,15,20,25,32,40,49,59,71,85,100,118,141,195};
    	int centLow_3p0GeV[16] ={5,7,9, 12,16,21,26,33,41,50,60,72,86 ,101,119,142};
        for(int i = 0; i<16 ; i++)
        {
		centHigh[i] = centHigh_3p0GeV[i];
		centLow[i]  = centLow_3p0GeV[i];
        }
    }
    if      (gRefMult>=centLow[15] && gRefMult<=centHigh[15]) centrality=15;
    else if (gRefMult>=centLow[14] && gRefMult<=centHigh[14]) centrality=14;
    else if (gRefMult>=centLow[13] && gRefMult<=centHigh[13]) centrality=13;
    else if (gRefMult>=centLow[12] && gRefMult<=centHigh[12]) centrality=12;
    else if (gRefMult>=centLow[11] && gRefMult<=centHigh[11]) centrality=11;
    else if (gRefMult>=centLow[10] && gRefMult<=centHigh[10]) centrality=10;
    else if (gRefMult>=centLow[9] && gRefMult<=centHigh[9]) centrality=9;
    else if (gRefMult>=centLow[8] && gRefMult<=centHigh[8]) centrality=8;
    else if (gRefMult>=centLow[7] && gRefMult<=centHigh[7]) centrality=7;
    else if (gRefMult>=centLow[6] && gRefMult<=centHigh[6]) centrality=6;
    else if (gRefMult>=centLow[5] && gRefMult<=centHigh[5]) centrality=5;
    else if (gRefMult>=centLow[4] && gRefMult<=centHigh[4]) centrality=4;
    else if (gRefMult>=centLow[3] && gRefMult<=centHigh[3]) centrality=3;
    else if (gRefMult>=centLow[2] && gRefMult<=centHigh[2]) centrality=2;
    else if (gRefMult>=centLow[1] && gRefMult<=centHigh[1]) centrality=1;
    else if (gRefMult>=centLow[0] && gRefMult<=centHigh[0]) centrality=0;
    else centrality = 16;

    return centrality;

}

Int_t CutManager::getMatchedToF()
{
    return mMatchedToF;
}

Int_t CutManager::getNpirm()
{
    return mN_prim;
}

Int_t CutManager::getNnonprim()
{
    return mN_non_prim;
}

