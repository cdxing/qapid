#ifndef QaPid_h
#define QaPid_h

#include "StMaker.h"
#include "TString.h"
#include "TFile.h"
// Configuration file reader
#include "../ConfigReader/ConfigReader.h"

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class CutManager;
class HistManager;
class StRefMultCorr;
class StPicoBTofPidTraits;

// run QA
//class TProfile;
//class TH1F;
//class TH2F;

class QaPid : public StMaker {
    public:
        QaPid(const char *name, StPicoDstMaker *picoMaker, char* jobid, std::string configFileName);
        virtual ~QaPid();

        virtual Int_t Init();
        virtual Int_t Make();
        virtual void  Clear(Option_t *opt="");
        virtual Int_t Finish();

        //int Centrality(int gRefMult);
        //int GetRunIndex(int runId);

        //bool isGoodTrigger(StPicoEvent const*) const;

    private:
	ConfigReader configs;
        static StRefMultCorr *mRefMultCorr;
        StPicoDstMaker *mPicoDstMaker;
        StPicoDst      *mPicoDst;
        StPicoEvent *mPicoEvent;
        CutManager  *mCutManager;
        HistManager *mHistManager;


        Int_t mEnergy;

        TString mOutPut_QAPID;

        TFile *mFile_QAPID;


        ClassDef(QaPid, 1)
};

#endif
