// singularity exec ~/Images/CCMadgraph.img make
// singularity exec ~/Images/CCMadgraph.img ./pythia
// zcat tsW.hepmc | singularity exec ~/Images/CCMadgraph.img /code/MG5_aMC_v2_6_0/Delphes/DelphesHepMC ./Cards/delphes_card_CMS.dat tsW.root -

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"

using Pythia8::Pythia;

int main()
{
  Pythia pythia("./pythia-xml");

  pythia.readString("StandardModel:Vtb = 0.0");
  pythia.readString("StandardModel:Vts = 1.0");

  pythia.readString("Top:gg2ttbar = on");
  pythia.readString("Top:qqbar2ttbar = on");

  pythia.readString("Main:numberOfEvents = 10000");

  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:frameType = 1");
  pythia.readString("Beams:eCM = 13000");

  pythia.settings.listChanged();
  HepMC::Pythia8ToHepMC ToHepMC;
  HepMC::IO_GenEvent ascii_io("tsW.hepmc", std::ios::out);

  pythia.init();
  for (int iEvent = 0; iEvent < 10000; ++iEvent) {
    pythia.next();
    HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
    ToHepMC.fill_next_event(pythia, hepmcevt);
    ascii_io << hepmcevt;
    delete hepmcevt;
  }

  pythia.stat();
  
  return 0;
}
