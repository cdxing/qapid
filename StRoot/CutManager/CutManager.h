#ifndef CutManager_h
#define CutManager_h

#include "TObject.h"
#include "TString.h"
#include "../ConfigReader/ConfigReader.h"
//#include <array>

class StPicoDst;
class StPicoTrack;
class StPicoBTofPidTraits;
class StRefMultCorr;
//class ConfigReader;

using namespace std;


class CutManager : public TObject
{
  public:
    CutManager(ConfigReader configs);
    ~CutManager();

    bool isGoodTrigger(StPicoDst* );
    bool passEventCut(StPicoDst* );
    bool isTofTrack(StPicoDst*, StPicoTrack* );
    bool passTrackEP(StPicoTrack*, float);
    bool passTrackBasic(StPicoTrack* );
    bool isProton( StPicoTrack* );
    bool isKaon(StPicoDst*, StPicoTrack* );
    bool isPion(StPicoDst*, StPicoTrack* );

    Int_t getCentrality(int gRefMult);
    Int_t getMatchedToF();
    Int_t getNpirm();
    Int_t getNnonprim();

  private:
    ConfigReader mConfigs;
    Int_t mEnergy;
    Int_t mMatchedToF;
    Int_t mN_prim;
    Int_t mN_non_prim;

    ClassDef(CutManager,1)
};
#endif
