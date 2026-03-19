#include "Event.hh"

#include "libpixie/event.hh"
#include "libpixie/reader.hh"

namespace FDSi {
const PIXIE::Event &Event::Set(const PIXIE::Reader &reader, int eventInd) {
  const PIXIE::Event &e = reader.events[eventInd];
  degai.nHits = 0;
  labr.nHits = 0;

  time = reader.measurements[e.fMeasurements[0]].eventTime;

  for (int i = 0; i < e.nMeas; ++i) {
    const PIXIE::Measurement &meas = reader.measurements[e.fMeasurements[i]];
    ////////////////////////
    // DEGAi PROCESSING //
    ////////////////////////
    // for each measurement, check if it's part of Clover
    int cloverID =
        degai.conf.CloverIDMap[meas.crateID][meas.slotID][meas.channelNumber];
    int crystal =
        degai.conf.CrystalMap[meas.crateID][meas.slotID][meas.channelNumber];

    if (!(cloverID == -1 || crystal == -1)) {
      degai.AddHit(meas, cloverID, crystal, e.fMeasurements[i]);
    } else {
      ////////////////////////
      // LaBr PROCESSING //
      ////////////////////////
      cloverID =
          labr.conf.CloverIDMap[meas.crateID][meas.slotID][meas.channelNumber];
      crystal =
          labr.conf.CrystalMap[meas.crateID][meas.slotID][meas.channelNumber];
      if (!(cloverID == -1 || crystal == -1)) {
        labr.AddHit(meas, cloverID, crystal);
      }
    }
  }

  // Process Hits -> Gammas (addback) here
  degai.MakeGammas();
  if (degai.nHits == MAX_HITS) {
    std::cout << "Warning! More than " << MAX_HITS << " hits" << std::endl;
    for (int j = 0; j < degai.nHits; ++j) {
      degai.hits[j].Print();
    }
  }
  labr.MakeGammas();
  if (labr.nHits == MAX_HITS) {
    std::cout << "Warning! More than " << MAX_HITS << " hits" << std::endl;
    for (int j = 0; j < labr.nHits; ++j) {
      labr.hits[j].Print();
    }
  }

  return e;
}
} // namespace FDSi
