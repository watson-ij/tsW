# Crab generation

Use CMSSW_8_0_30 for generation.

Use:

```
crab submit --config=pset_generation.py
```

to submit and

```
crab status
```

to check the status.

When its done, add files to and run:

```
cmsRun RECO_RAW2DIGI_L1Reco_RECO_EI_PAT.py
```

To get the RECO that can be run with hadAOD