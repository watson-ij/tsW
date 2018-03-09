// singularity exec ~/Images/Madgraph.img make
// singularity exec ~/Images/Madgraph.img ./pythia
// cat tsW.hepmc | singularity exec ~/Images/Madgraph.img /code/MG5_aMC_v2_6_0/Delphes/DelphesHepMC ./Cards/delphes_card_CMS.dat tsW.root -

#include "unistd.h"
#include <string>

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"

using std::string;

using Pythia8::Pythia;

static bool doS = false;
static bool doH = false;
static bool doLam = false;
static int nEvents = 1000;

int main(int argc, char *argv[])
{
  int c;
  while ((c = getopt(argc, argv, "sbhln:")) != -1) {
    switch (c) {
    case 's':
      doS = true;
      break;
    case 'b':
      doS = false;
      break;
    case 'h':
      doH = true;
      break;
    case 'l':
      doLam = true;
      break;
    case 'n':
      nEvents = atoi(optarg);
      break;
    case '?':
      std::cout << "Bad Option: " << optopt << std::endl;
      exit(12);
    }
  }

  if ((argc - optind) != 1) {
    std::cout << "Requires output hepmc file argument. Optional -s flag for tt->sWsW generation, -b flag for tt->bWbW generation, -h flag for tt->bWsW, -l flag for Lambda_b requirement" << std::endl;
    return 1;
  }

  string out = string(argv[optind]);
  
  Pythia pythia("./pythia-xml");

  if (doS) {
    std::cout << "t->sW generation" << std::endl;
    pythia.readString("StandardModel:Vtb = 0.0");
    pythia.readString("StandardModel:Vts = 1.0");
    pythia.readString("StandardModel:Vtd = 0.0");
    doH = false;
  } else if (doH)  { //doH
    std::cout << "tt->sWbW generation" << std::endl;
    pythia.readString("StandardModel:Vtb = 0.7071");
    pythia.readString("StandardModel:Vts = 0.7071");
    pythia.readString("StandardModel:Vtd = 0.0");
  } else { //doB
    std::cout << "t->bW generation" << std::endl;
    pythia.readString("StandardModel:Vtb = 1.0");
    pythia.readString("StandardModel:Vts = 0.0");
    pythia.readString("StandardModel:Vtd = 0.0");
  }
  std::cout << nEvents << " events" << std::endl;
  std::cout << std::endl;
  std::cout << "-------------------------------------------------------------------------" << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;

  pythia.readString("Top:gg2ttbar = on");
  pythia.readString("Top:qqbar2ttbar = on");

  pythia.readString("Main:numberOfEvents = 20000");

  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:frameType = 1");
  pythia.readString("Beams:eCM = 13000");

  pythia.readString("24:onMode = 0"); // turn all W decays off
  pythia.readString("24:onIfAny = 13"); // turn on mu decay mode
  pythia.readString("24:onIfAny = 11"); // turn on electron decay mode

  pythia.settings.listChanged();
  HepMC::Pythia8ToHepMC ToHepMC;
  HepMC::IO_GenEvent ascii_io(out.c_str(), std::ios::out);

  pythia.init();
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
    while (true) {
      pythia.next();
      HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
      ToHepMC.fill_next_event(pythia, hepmcevt);
      if (doH || doL) { // check bWsW topology
	bool seenS = false, seenB = false, seenL = false;
	for (auto p = hepmcevt->particles_begin(); p != hepmcevt->particles_end(); ++p) {
	  if (abs((*p)->status()) == 23 && abs((*p)->pdg_id()) == 3)
	    seenS = true;
	  if (abs((*p)->status()) == 23 && abs((*p)->pdg_id()) == 5)
	    seenB = true;
	  if (abs((*p)->pdg_id()) == 5122)
	    seenB = true;
	}

	// discard events without an s and b
	if (doH) {
	  if (!seenS || !seenB) {
	    delete hepmcevt;
	    continue;
	  }
	} else if (doL) {
	  if (!seenL) {
	    delete hepmcevt;
	    continue;
	  }
	} else {
	  cout << "Unreachable" << endl;
	}
      }
      ascii_io << hepmcevt;
      delete hepmcevt;
      break;
    }
  }

  pythia.stat();
  
  return 0;
}
