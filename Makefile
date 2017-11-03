pythia_tsW : pythia_tsW.cc
	g++ pythia_tsW.cc -o pythia_tsW `/code/MG5_aMC_v2_6_0/HEPTools/pythia8/bin/pythia8-config --cflags --libs` -I/code/MG5_aMC_v2_6_0/HEPTools/hepmc/include -L/code/MG5_aMC_v2_6_0/HEPTools/hepmc/lib -lHepMC -ldl

pythia_tbW : pythia_tbW.cc
	g++ pythia_tbW.cc -o pythia_tbW `/code/MG5_aMC_v2_6_0/HEPTools/pythia8/bin/pythia8-config --cflags --libs` -I/code/MG5_aMC_v2_6_0/HEPTools/hepmc/include -L/code/MG5_aMC_v2_6_0/HEPTools/hepmc/lib -lHepMC -ldl

run_tsW : pythia_tsW
	LD_LIBRARY_PATH=/code/MG5_aMC_v2_6_0/HEPTools/hepmc/lib:$LD_LIBRARY_PATH ./pythia_tsW

run_tbW : pythia_tbW
	LD_LIBRARY_PATH=/code/MG5_aMC_v2_6_0/HEPTools/hepmc/lib:$LD_LIBRARY_PATH ./pythia_tbW
