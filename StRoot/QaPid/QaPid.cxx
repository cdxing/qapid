#include "QaPid.h"
#include "StRoot/CutManager/CutManager.h"
#include "StRoot/ConstManager/ConstManager.h"
//#include "../HistManager/HistManager.h"
#include "StRoot/HistManager/HistManager.h"
// Load STARLibrary header files
#include "StMaker.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoHelix.h"
#include "StRoot/StPicoEvent/StPicoBbcHit.h"
#include "StRoot/StPicoEvent/StPicoEpdHit.h"
#include "StRoot/StPicoEvent/StPicoBTofHit.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StRoot/StPicoEvent/StPicoTrackCovMatrix.h"
//#include "StRoot/run/run.h"
#include "StThreeVectorF.hh"
#include "TH1F.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TFile.h"
#include "TVector3.h"
#include "TMath.h"
#include "StMessMgr.h"
#include <algorithm>
#include <array>
//#include "StRoot/StRefMultCorr/StRefMultCorr.h"
//#include "StRoot/StRefMultCorr/CentralityMaker.h"
ClassImp(QaPid)

    //StRefMultCorr* QaPid::mRefMultCorr = NULL;
    //-----------------------------------------------------------------------------
    QaPid::QaPid(const char* name, StPicoDstMaker *picoMaker, char* jobid, std::string configFileName)
: StMaker(name)
{
    configs.read(configFileName);
    mPicoDstMaker = picoMaker;
    mPicoDst = 0;
    mEnergy = 3;

    mOutPut_QAPID=Form("%s_QAPID.root",jobid);
}

//----------------------------------------------------------------------------- 
QaPid::~QaPid()
{ /*  */ }

//----------------------------------------------------------------------------- 
Int_t QaPid::Init() 
{
    //if(!mRefMultCorr)
    //{
    //    mRefMultCorr = CentralityMaker::instance()->getRefMultCorr();
    //}
    
    mCutManager = new CutManager(configs);
    mHistManager = new HistManager();
    mHistManager->InitQAPID();
    mFile_QAPID = new TFile(mOutPut_QAPID.Data(),"RECREATE");
    mFile_QAPID->cd();
    return 0;

}

//----------------------------------------------------------------------------- 
Int_t QaPid::Finish() 
{


    if(mOutPut_QAPID != "")
    {
        mFile_QAPID->cd();
        mHistManager->WriteQAPID();
        mFile_QAPID->Close();
    }
    return kStOK;
}

//----------------------------------------------------------------------------- 
void QaPid::Clear(Option_t *opt) {
}

//----------------------------------------------------------------------------- 
Int_t QaPid::Make() 
{
    if(!mPicoDstMaker) 
    {
        LOG_WARN << " No PicoDstMaker! Skip! " << endm;
        return kStWarn;
    }

    mPicoDst = mPicoDstMaker->picoDst();
    if(!mPicoDst) 
    {
        LOG_WARN << " No PicoDst! Skip! " << endm;
        return kStWarn;
    }

    mPicoEvent = (StPicoEvent*)mPicoDst->event();
    if(!mPicoEvent)
    {
        LOG_WARN << " No PicoEvent! Skip! " << endm;
        return kStWarn;
    }
    const Int_t nTracks = mPicoDst->numberOfTracks();
    Int_t TrkMult = 0;
    for(Int_t i = 0; i < nTracks; i++) // track loop for TrkMult
    {
        StPicoTrack *track = (StPicoTrack*)mPicoDst->track(i);
        if(!track->isPrimary()) continue; // Only Primary Tracks
        mHistManager->FillTrackQA(track,(TVector3)mPicoEvent->primaryVertex());
        //StPicoPhysicalHelix helix = track->helix(mField);
        //Float_t dca = helix.geometricSignedDistance(mVertexPos);
	TrkMult ++;
    }

    // RefMult
    Int_t runId = mPicoEvent->runId();

    //cout << "runID = " << runId << endl;
    Int_t refMult = mPicoEvent->refMult();
    Float_t vz = mPicoEvent->primaryVertex().Z();
    Float_t vx = mPicoEvent->primaryVertex().X();
    Float_t vy = mPicoEvent->primaryVertex().Y();

    Float_t vzvpd = mPicoEvent->vzVpd();
    Int_t TOF_Mul = mPicoEvent->btofTrayMultiplicity();
    Int_t nMatchedToF = mPicoEvent->nBTOFMatch();
    Float_t zdcX = mPicoEvent->ZDCx();
    
    mHistManager->FillEventQA(mPicoEvent->primaryVertex(),refMult,TOF_Mul,TrkMult);
    mHistManager->FillEventCut(0);
    // vz sign
    Int_t vz_sign;
    if(vz > 0.0)
    {
        vz_sign = 0;
    }
    else
    {
        vz_sign = 1;
    }

    // runIndex
    //const int runIndex = GetRunIndex(runId);

    // Event Cut
    if(mCutManager->passEventCut(mPicoDst)) // event cut
    {
	
        mHistManager->FillEventQaCut(mPicoEvent->primaryVertex(),refMult,TOF_Mul,TrkMult);
        mHistManager->FillEventCut(1);
        //mRefMultCorr->init(runId);
        //mRefMultCorr->initEvent(refMult, vz, zdcX);
        //const Int_t cent9 = mRefMultCorr->getCentralityBin9();
        //const Double_t reweight = mRefMultCorr->getWeight();
	//std::cout << "refMult: " << refMult << std::endl;
	//std::cout << "TrkMult: " << TrkMult << std::endl;
        const int cent16 = mCutManager->getCentrality(TrkMult);
	
	//std::cout << "cent16: " << cent16 << std::endl;
        const double reweight = 1.0;
        if(cent16 >  15 || cent16 < 0) return 0;
        mHistManager->FillEventCent(cent16);
        mHistManager->FillEventCut(2);

        //const Int_t nToFMatched = mCutManager->getMatchedToF();
        TVector3 mVertexPos = mPicoDst->event()->primaryVertex();
        float mField = mPicoEvent->bField();
        Int_t N_pp = 0, N_pm = 0, N_kp = 0, N_km = 0, N_pr = 0;

        //cout << "nTracks = " << nTracks << endl;
        for(Int_t i = 0; i < nTracks; i++) // track loop
        {
            StPicoTrack *track = (StPicoTrack*)mPicoDst->track(i);
            if(!track->isPrimary()) continue; // Only Primary Tracks
            mHistManager->FillTrackCut(0);
            if(!mCutManager->passTrackBasic(track)) continue;
            mHistManager->FillTrackCut(1);
            mHistManager->FillTrackPhysics(track );
		
            //StPicoPhysicalHelix helix = track->helix(mField);
            //Float_t dca = helix.geometricSignedDistance(mVertexPos);
            Short_t  s_charge = track->charge();
            Float_t dca=track->gDCA(mVertexPos).Mag();
	    if(mCutManager->isTofTrack(mPicoDst,track)) mHistManager->FillTrackTof(mPicoDst,track);
	    if(mCutManager->isProton(mPicoDst,track))
	    {
	   	mHistManager->FillProton(mPicoDst,track,configs.y_mid); 
		N_pr++;
	    }
	    if(mCutManager->isKaon(mPicoDst,track))
	    {
	   	mHistManager->FillKaon(mPicoDst,track,configs.y_mid); 
		if(s_charge >0) N_kp++;
		else N_km++;
	    }
	    if(mCutManager->isPion(mPicoDst,track))
	    {
	   	mHistManager->FillPion(mPicoDst,track,configs.y_mid); 
		if(s_charge >0) N_pp++;
		else N_pm++;
	    }
        } // track loop
	mHistManager->FillPIDMult(N_pp , N_pm , N_kp , N_km , N_pr ); 
    } // event cut

    return kStOK;
}
/*
int QaPid::GetRunIndex(int runID)
{
    int runIndex=-999;
    for(int i=0; i<2704; i++)
    {
        if(runID==numbers[i])
        {
            runIndex=i;
        }
    }
    if(runIndex == -999) cout << "Run numbers are not found!!!" << endl;
    return runIndex;
}

int QaPid::Centrality(int gRefMult )
{
    int centrality;
    int centFull[9]={4, 9,17,30,50,78, 116,170,205};
    if      (gRefMult>=centFull[8]) centrality=8;
    else if (gRefMult>=centFull[7]) centrality=7;
    else if (gRefMult>=centFull[6]) centrality=6;
    else if (gRefMult>=centFull[5]) centrality=5;
    else if (gRefMult>=centFull[4]) centrality=4;
    else if (gRefMult>=centFull[3]) centrality=3;
    else if (gRefMult>=centFull[2]) centrality=2;
    else if (gRefMult>=centFull[1]) centrality=1;
    else if (gRefMult>=centFull[0]) centrality=0;
    else centrality = 9;

    return centrality;
}*/
