# Use singularity image to run make, i.e.
# singularity ~/Images/Madgraph.img make run_tsW_dilep

pythia_tsW : pythia_tsW.cc
	g++ pythia_tsW.cc -o pythia_tsW `/code/MG5_aMC_v2_6_0/HEPTools/pythia8/bin/pythia8-config --cflags --libs` -I/code/MG5_aMC_v2_6_0/HEPTools/hepmc/include -L/code/MG5_aMC_v2_6_0/HEPTools/hepmc/lib -lHepMC -ldl

pythia_tbW : pythia_tbW.cc
	g++ pythia_tbW.cc -o pythia_tbW `/code/MG5_aMC_v2_6_0/HEPTools/pythia8/bin/pythia8-config --cflags --libs` -I/code/MG5_aMC_v2_6_0/HEPTools/hepmc/include -L/code/MG5_aMC_v2_6_0/HEPTools/hepmc/lib -lHepMC -ldl

run_tsW : pythia_tsW
	LD_LIBRARY_PATH=/code/MG5_aMC_v2_6_0/HEPTools/hepmc/lib:$LD_LIBRARY_PATH ./pythia_tsW

run_tbW : pythia_tbW
	LD_LIBRARY_PATH=/code/MG5_aMC_v2_6_0/HEPTools/hepmc/lib:$LD_LIBRARY_PATH ./pythia_tbW

pythia_tsW_dilep : pythia_tsW_dilep.cc
	g++ pythia_tsW_dilep.cc -o pythia_tsW_dilep `/code/MG5_aMC_v2_6_0/HEPTools/pythia8/bin/pythia8-config --cflags --libs` -I/code/MG5_aMC_v2_6_0/HEPTools/hepmc/include -L/code/MG5_aMC_v2_6_0/HEPTools/hepmc/lib -lHepMC -ldl

pythia_tbW_dilep : pythia_tbW_dilep.cc
	g++ pythia_tbW_dilep.cc -o pythia_tbW_dilep `/code/MG5_aMC_v2_6_0/HEPTools/pythia8/bin/pythia8-config --cflags --libs` -I/code/MG5_aMC_v2_6_0/HEPTools/hepmc/include -L/code/MG5_aMC_v2_6_0/HEPTools/hepmc/lib -lHepMC -ldl

run_tsW_dilep : pythia_tsW_dilep
	LD_LIBRARY_PATH=/code/MG5_aMC_v2_6_0/HEPTools/hepmc/lib:$LD_LIBRARY_PATH ./pythia_tsW_dilep
	/code/MG5_aMC_v2_6_0/Delphes/DelphesHepMC delphes_card_CMS.tcl tsW_dilep.root tsW_dilep.hepmc

run_tbW_dilep : pythia_tbW_dilep
	LD_LIBRARY_PATH=/code/MG5_aMC_v2_6_0/HEPTools/hepmc/lib:$LD_LIBRARY_PATH ./pythia_tbW_dilep
	/code/MG5_aMC_v2_6_0/Delphes/DelphesHepMC delphes_card_CMS.tcl tbW_dilep.root tbW_dilep.hepmc
