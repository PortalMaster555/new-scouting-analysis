import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
import mplhep
import hist
import glob
import pickle
from tqdm import tqdm

# MUON = "Vtx"
# TRIGGER = ["DST_PFScouting_DoubleMuonVtx", "DST_PFScouting_DoubleMuonNoVtx"]
# pickleString = "initial_studies_2025/histos_2025_TrgOR_Vtx.pkl"

MUON = "NoVtx"
TRIGGER = ["DST_PFScouting_DoubleMuonNoVtx"]
pickleString = "initial_studies_2025/histos_2025_TrgNoVtx_NoVtx.pkl"


## Open the file
selected_branches = []
selected_branches += TRIGGER
selected_branches.append("nScoutingMuon%sDisplacedVertex"%(MUON))
selected_branches.append("nScoutingMuon%s"%(MUON))
selected_branches.append("ScoutingMuon%s_pt"%(MUON))
selected_branches.append("ScoutingMuon%s_eta"%(MUON))
selected_branches.append("ScoutingMuon%s_phi"%(MUON))
selected_branches.append("ScoutingMuon%s_charge"%(MUON))


indir = "/eos/user/f/fernance/DST-Muons/samples"
outdir = "/afs/cern.ch/user/d/dspringb/private/CMSSW_15_0_5/src/new-scouting-analysis"
isRereadingFullFile = False
if isRereadingFullFile:
    files = glob.glob(f"{indir}/*.root")
    pickleIndex = 0
    for f in tqdm(files, desc="Progress"):
        events = uproot.open(f)["Events"].arrays(selected_branches, library="ak")
        print(f"> Opened sample with {len(events)} events") # 1_763_889 events

        ## Filter by trigger
        #TrgOR
        mask = ak.zeros_like(events[TRIGGER[0]], dtype=bool)
        for trigger in TRIGGER:
            mask = np.array(mask) | np.array(events[trigger])
        events = events[mask]
        
        print(f"> {len(events)} events survive the filter") # 500_206 events survive from Vtx TrgOR
        # 356_664 events survive NoVtx TrgNoVtx
        print("~\nFirst event: events[0]")
        print(events[0])
        with open (outdir+"/large_pickles/events[%s]Pickle.pkl"%(MUON), "wb") as pickleOut:
            pickle.dump(events, pickleOut)
        pickleIndex += 1
else: # if not rereading the original unfiltered file
    with open (outdir+"/large_pickles/events%sPickle.pkl"%(MUON), "rb") as pickleIn:
        events = pickle.load(pickleIn)
    print(f"> Opened filtered sample with {len(events)} events") # 500_206 events Vtx TrgOR, # 356_664 events NoVtx TrgNoVtx

## Filter by number of muons
events = events[ events["nScoutingMuon%s"%(MUON)] > 1]

print(ak.fields(events[0])) # get list of fields again since i forgot
print(events[0]["nScoutingMuon%sDisplacedVertex"%(MUON)])