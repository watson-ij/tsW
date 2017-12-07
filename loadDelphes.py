import ROOT

ROOT.gInterpreter.AddIncludePath('/code/MG5_aMC_v2_6_0/Delphes/external/')
ROOT.gInterpreter.AddIncludePath('/code/MG5_aMC_v2_6_0/Delphes/external/ExRootAnalysis')
ROOT.gInterpreter.AddIncludePath('/code/MG5_aMC_v2_6_0/Delphes/')
ROOT.gInterpreter.AddIncludePath('/code/MG5_aMC_v2_6_0/Delphes/classes/')
ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
ROOT.gSystem.Load("/code/MG5_aMC_v2_6_0/Delphes/libDelphes.so")
