#include "CutManager.h"
#include "StRoot/ConstManager/ConstManager.h"
#include "StRoot/StPicoEvent/StPicoDst.h"  // shaowei  18c
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StMessMgr.h"

ClassImp(CutManager)

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
    Bool_t b_good_trig = false;
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
    if(sqrt(vx*vx+(vy+2)*(vy+2)) > mConfigs.r_vtx)
    {
        return kFALSE;
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
bool CutManager::isProton( StPicoTrack *track)
{
  Double_t d_TPCnSigmaProton = track->nSigmaProton();
  Short_t  s_charge = track->charge();
  Bool_t proton = (d_TPCnSigmaProton > mConfigs.nSig_pr_low) && (d_TPCnSigmaProton < mConfigs.nSig_pr_high) && (s_charge > 0);
  return proton;
}

bool CutManager::isKaon(StPicoDst *pico, StPicoTrack *track)
{
  Double_t d_TPCnSigmaKaon   = track->nSigmaKaon();
  Short_t  s_charge = track->charge();
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
  Short_t  s_charge = track->charge();
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

