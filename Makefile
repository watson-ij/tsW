pythia : pythia.cc
	g++ pythia.cc -o pythia `/code/MG5_aMC_v2_6_0/HEPTools/pythia8/bin/pythia8-config --cflags --libs` -I/code/MG5_aMC_v2_6_0/HEPTools/hepmc/include -L/code/MG5_aMC_v2_6_0/HEPTools/hepmc/lib -lHepMC -ldl

pythia_tbW : pythia_tbW.cc
	g++ pythia_tbW.cc -o pythia_tbW `/code/MG5_aMC_v2_6_0/HEPTools/pythia8/bin/pythia8-config --cflags --libs` -I/code/MG5_aMC_v2_6_0/HEPTools/hepmc/include -L/code/MG5_aMC_v2_6_0/HEPTools/hepmc/lib -lHepMC -ldl

run : pythia
	LD_LIBRARY_PATH=/code/MG5_aMC_v2_6_0/HEPTools/hepmc/lib:$LD_LIBRARY_PATH ./pythia

run_tbW : pythia_tbW
	LD_LIBRARY_PATH=/code/MG5_aMC_v2_6_0/HEPTools/hepmc/lib:$LD_LIBRARY_PATH ./pythia_tbW
