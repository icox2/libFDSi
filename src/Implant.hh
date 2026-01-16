#ifndef FDSI_IMPLANT_HH
#define FDSI_IMPLANT_HH

#include <string>

#include "libpixie/measurement.hh"

#include "Definitions.hh"

namespace FDSi {
  const int MAX_DSSDS = 2;
  const int NUM_GAINS = 3; // HG, MG, LG
  const int MAX_DSSD_STRIPS = 128;

  class ImplantChannel {
    public:
      virtual void SetMeas(PIXIE::Measurement &meas, int indx) {}
  };      

  class pspmtCal {
    public:
      bool loaded = false;
      float *ens[2304];
      float *posxs[2304];
      float *posys[2304];
      
      float *cal_xposx;
      float *cal_xposy;
  
      float *cal_yposx;
      float *cal_yposy;

      int lens[2304];
    public:
      pspmtCal() : loaded(false) { cal_xposx=new float[128*48*48]; cal_xposy=new float[128*48*48]; cal_yposx=new float[128*48*48]; cal_yposy=new float[128*48*48];} ;
      void read_calgraphs(std::string file);
      void make_map();
      std::pair<float,float> interp_en(float val, int xind, int yind);
      int find_indx(float x, float y, int ien, int axis);

      float interp(float val, float *x, float *y);
      std::pair<float,float> interp2(float x, float y, float x1, float y1, float x2, float y2, float x3, float y3, float x4, float y4);
      std::pair<float,float> cal(float en, float x, float y); 
  };

  class Implant : ImplantChannel {
    public:
      double TOF;
      double dE;
      double dE2; 
      unsigned long long time;

      double energy[3][5];

      double xpos[3];
      double ypos[3];
      double xcal[3];
      double ycal[3];

      double thresh[5];
      int upperthresh;

      bool present = false;
      bool valid = false;
      bool validTOF = false;
      bool validCut = false;
      int nHits = 0;
      int cutID = 0;

      unsigned long long start;
      unsigned long long stop;

      int fired[5];

      uint32_t tracelen;
      uint16_t *trace; 
    public:
      Implant();
      void SetMeas(PIXIE::Measurement &meas, int indx);
      int firedAnodes();
      int firedDynode();
      int EnergySum(int type);
      int SetPos(int type, pspmtCal &cal);
      void Calibrate(int type, pspmtCal &cal);
      void Reset();
  };

  class SiPMImplant : ImplantChannel {
  public:
    double TOF;
    double dE;
    unsigned long long time;

    int energy[8][8];
    int dynodeEn;
    int dynodeFired;
    int nFired;

    double xpos;
    double ypos;
    double xcal;
    double ycal;

    double thresh;
    int upperthresh;

    bool present = false;
    bool valid = false;
    bool validCut = false;
    int nHits = 0;
    int cutID = 0;

    int fired[8][8];

    uint32_t tracelen;
    uint16_t *trace; 
  public:
    SiPMImplant();
    void SetMeas(PIXIE::Measurement &meas, int indx);
    void SetDynodeMeas(PIXIE::Measurement &meas);
    void SetAnodeMeas(PIXIE::Measurement &meas, int indx);
    int firedAnodes();
    int firedDynode();
    int EnergySum();
    int SetPos();
    void Calibrate(int type, pspmtCal &cal);
    void Reset();
  };

  class DSSDConf {
  public:
    int numStrips; // Total number of DSSD strips
    std::string conf_name;
  };

  struct DSSDhit {
    int indx; // strip index combines DSSD number and gain 
    double energy;
    unsigned long long time;
  };

  class DSSDBeta : ImplantChannel {
  public:
    unsigned long long time;
    
    // Get the largest energy from front (HG), back (MG) and back (LG)
    // Ideally for implants, best TKE resolution will be the LG
    // front (HG) should saturate.. energy not needed but strip # is..
    // Arranged by [DSSD][Gain][Strip]
    double energy; // decay energy from HG

    // xpos determined by strip with highest energy (or saturated) from HG 
    int xpos;
    // ypos determined by strip with highest energy (or saturated) from LG 
    int ypos;
    // tells which dssd the imp was detected in (last along beam to give signal)
    int zpos; 
    // position calibration should not be needed..

    double upperthresh;

    bool present = false;
    bool valid = false;
    int nHits = 0;

    double tdiff = -999; //correlation time
    int nCuts;
    int firstCut;
    int cutIDs[MAX_BETACUTS];
    double tdiffs[MAX_BETACUTS];
  
    bool promptGamma;
    int gammaIndex;
    double gammaEnergy;

  public:
    DSSDBeta();

    void Reset();
  };

  class DSSDImplant : ImplantChannel {
  public:
    double TOF;
    double dE;
    double dE2; 
    unsigned long long time;

    unsigned long long start;
    unsigned long long stop;
    
    // Get the largest energy from front (HG), back (MG) and back (LG)
    // Ideally for implants, best TKE resolution will be the LG
    // front (HG) should saturate.. energy not needed but strip # is..
    double energy; // total implant energy from LG
    double tke;

    // xpos determined by strip with highest energy (or saturated) from HG 
    int xpos;
    // ypos determined by strip with highest energy (or saturated) from LG 
    int ypos;
    // tells which dssd the imp was detected in (last along beam to give signal)
    int zpos; 
    // also store the first DSSD x-y to track beam path
    int xpos0;
    int ypos0;

    double lowerthresh = 1;

    bool present = false;
    bool validTOF = false;
    bool valid = false;
    bool validCut = false;
    int nHits = 0;
    int cutID = 0;
  
    bool promptGamma;
    int gammaIndex;
    double gammaEnergy;

  public:
    DSSDImplant();

    void Reset();
  };

  class DSSD : ImplantChannel {
    public:
      std::vector<DSSDhit> hits;
      double thresh[MAX_DSSDS][NUM_GAINS][MAX_DSSD_STRIPS]; // thresholds for individual strips
      double upperthresh[NUM_GAINS]; // upper threshold for decay
      
      DSSD();
      void SetMeas(PIXIE::Measurement &meas, int indx);
      
      // Take in both beta and implant then fill corresponding info
      void DSSDHitAnalysis(DSSDBeta *beta, DSSDImplant *implant);

      void Reset(){ hits.clear(); }
  };

  class Beta : ImplantChannel {
    public:
      unsigned long long time = 0;
      bool present = false;
      bool valid = false;
      int nHits = 0;
      double tdiff = -999; //correlation time
      int nCuts;
      int firstCut;
      int cutIDs[MAX_BETACUTS];
      double tdiffs[MAX_BETACUTS];

      double energy[3][5];
      int fired[5];

      double xpos[3];
      double ypos[3];
      double xcal[3];
      double ycal[3];

      double thresh[5];
      int upperthresh;

      uint16_t *trace[5];
      int tracelen;
    
      bool promptGamma;
      int gammaIndex;
      double gammaEnergy;

    public:
      Beta();
      void SetMeas(PIXIE::Measurement &meas, int indx);
      int firedDynode();
      int firedAnodes();
      int EnergySum(int type);
      void Reset();
      int SetPos(int type, pspmtCal &cal);
      void Calibrate(int type, pspmtCal &cal);
  };

  class SiPMBeta : ImplantChannel {
  public:
    unsigned long long time;

    int energy[8][8];
    int dynodeEn;
    int dynodeFired;
    int nFired;

    double xpos;
    double ypos;
    double xcal;
    double ycal;

    double thresh;
    int upperthresh;

    bool present = false;
    bool valid = false;
    bool validCut = false;
    int nHits;
    double tdiff = -999; //correlation time
    int nCuts;
    int firstCut;
    int cutIDs[MAX_BETACUTS];
    double tdiffs[MAX_BETACUTS];

    int fired[8][8];

    uint32_t tracelen;
    uint16_t *trace;

    bool promptGamma;
    int gammaIndex;
    double gammaEnergy;
    
  public:
    SiPMBeta();
    void SetMeas(PIXIE::Measurement &meas, int indx);
    void SetDynodeMeas(PIXIE::Measurement &meas);
    void SetAnodeMeas(PIXIE::Measurement &meas, int indx);
    int firedAnodes();
    int firedDynode();
    int EnergySum();
    int SetPos();
    void Calibrate(int type, pspmtCal &cal);
    void Reset();
  };
  
  class IonTrigger : ImplantChannel { 
    public:
      double thresh = 0;
      double energy = 0;
      unsigned long long time = 0;
      bool valid = 0;
      int fired = 0;
    public:
      IonTrigger() : energy(0), time(0), valid(false), fired(0) {}
      void SetMeas(PIXIE::Measurement &meas, int indx);
      void Reset();
  };

  class Pin : ImplantChannel {
    public:
      double energy[2];
      unsigned long long time[2];
      int fired[2];
      bool valid = false;
      double ecal;

      double gain;
      double offset;

      Pin();
      void SetMeas(PIXIE::Measurement &meas, int indx);
      void Reset();
      void SetCal(double g, double o);
      void Calibrate();

  };

  class Scint : ImplantChannel {
    public:
      double energy[4];
      unsigned long long time[4];
      int fired[4];
      double thresh[4];
    int nchans = 4;

      unsigned long long avtime = 0;
      bool valid = false;
      Scint();
    Scint(int nch);
      void SetMeas(PIXIE::Measurement &meas, int indx);
      void Reset();
  };

  class PPAC  : ImplantChannel {
    public:
      unsigned long long time[5];
      double energy[5];
    double thresh[5];
      int fired[5];
    int present;

      unsigned long long avtime;
      bool valid = false;
      PPAC();
    PPAC(int nch); 
      void SetMeas(PIXIE::Measurement &meas, int indx);
      void Reset();
      void validate();
  };

  class ImplantEvent {
    public:
      //maps (chan,slot,crate)-> index
      int indexMap[FD_MAX_CRATES][FD_MAX_SLOTS_PER_CRATE][FD_MAX_CHANNELS_PER_BOARD];
    int subIndexMap[FD_MAX_CRATES][FD_MAX_SLOTS_PER_CRATE][FD_MAX_CHANNELS_PER_BOARD];
      //maps (ID)->(object,index)
      ImplantChannel **implantMap[FD_MAX_CRATES][FD_MAX_SLOTS_PER_CRATE][FD_MAX_CHANNELS_PER_BOARD];

      int nStored;
      int impCtr = 0;
    int sipmImpCtr = 0;
      int betaCtr = 0;
    int sipmBetaCtr = 0;
    Implant *imps;
    SiPMImplant *sipmImps;
    SiPMBeta *sipmBetas;
    DSSDImplant *dssdImps;
    DSSDBeta *dssdBetas;
    Beta *betas;

    bool sipm_present = false;
    bool dssd_present = false;
    Beta *beta = NULL;
    Implant *imp = NULL;
    SiPMImplant *sipmImp = NULL;
    SiPMBeta *sipmBeta = NULL;
    DSSD *dssd = NULL; // this dssd only stores the hits which are processed later 
    DSSDBeta *dssdBeta = NULL;
    DSSDImplant *dssdImp = NULL;

      IonTrigger *fit;     
      IonTrigger *rit;     
      Pin *pin[4];
      Scint *scint[4];
      PPAC *ppac[6];

      int nPins;
      int nScints;
      int nPPACs;

    public:
      ImplantEvent();
    ImplantEvent(int nst, bool sipm = false, bool dssd_pres = false);
      ImplantEvent(int nPins, int nScints, int nPPACs);
    void Init(int nst, bool sipm = false, bool dssd_pres = false);
      void ReadConf(std::string conffile);
      void Set(PIXIE::Measurement &meas) {
        ImplantChannel **ch = implantMap[meas.crateID][meas.slotID][meas.channelNumber];
        int indx = indexMap[meas.crateID][meas.slotID][meas.channelNumber];
        int subindx = subIndexMap[meas.crateID][meas.slotID][meas.channelNumber];
        if (ch != NULL && subindx >= 0) {
          (*ch)->SetMeas(meas, subindx);
        }
      }
      void SetBetaThresh(int indx, float val);
      void SetImplantThresh(int indx, float val);
      void SetDSSDThresh(int indx, float val);
    void SetPPACThresh(int indx, int subindx, float val);

      void Reset();
  };

}

#endif
