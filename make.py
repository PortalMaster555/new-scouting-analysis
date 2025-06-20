import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
import mplhep
import hist
import glob
import pickle
from tqdm import tqdm

MUON = "Vtx"
TRIGGER = ["DST_PFScouting_DoubleMuonVtx", "DST_PFScouting_DoubleMuonNoVtx"]
pickleString = "initial_studies_2025/histos_2025_TrgOR_Vtx.pkl"

# MUON = "NoVtx"
# TRIGGER = ["DST_PFScouting_DoubleMuonVtx", "DST_PFScouting_DoubleMuonNoVtx"]
# pickleString = "initial_studies_2025/histos_2025_TrgOR_NoVtx.pkl"


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
files = glob.glob(f"{indir}/*.root")
for f in tqdm(files, desc="Progress"):
    events = uproot.open(f)["Events"].arrays(selected_branches, library="ak")
    print(f"> Opened sample with {len(events)} events")
    for event in events:
        print(event)