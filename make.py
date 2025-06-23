import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
import mplhep
import hist
import glob
import pickle
import sys
from tqdm import tqdm

try:
    if sys.argv[1].lower() == "--vtx":
        MUON = "Vtx"
        TRIGGER = ["DST_PFScouting_DoubleMuonVtx"]
        pickleString = "zoom_2025/histos_2025_TrgVtx_Vtx.pkl"
    elif sys.argv[1].lower() == "--novtx":
        MUON = "NoVtx"
        TRIGGER = ["DST_PFScouting_DoubleMuonNoVtx"]
        pickleString = "zoom_2025/histos_2025_TrgNoVtx_NoVtx.pkl"
    else:
        print("Please supply only --vtx or --novtx.")
        sys.exit(1)
except IndexError:
    print("Please supply a flag (--vtx, --novtx).")
    sys.exit(1)

isRereadingFullFile = False
indir = "/eos/user/f/fernance/DST-Muons/samples"
outdir = "/afs/cern.ch/user/d/dspringb/private/CMSSW_15_0_5/src/new-scouting-analysis"
try:
    if sys.argv[2].lower() == "--reread":
        isRereadingFullFile = True
except IndexError: # if --reread not specified, assume using same file.
    print("Using stored data from "+ outdir + "/large_pickles/events[%s]Pickle.pkl"%(MUON))
finally:
    if isRereadingFullFile:
        ## Open the file
        selected_branches = []
        selected_branches += TRIGGER
        selected_branches.append("nScoutingMuon%sDisplacedVertex"%(MUON))
        selected_branches.append("nScoutingMuon%s"%(MUON))
        selected_branches.append("ScoutingMuon%s_pt"%(MUON))
        selected_branches.append("ScoutingMuon%s_eta"%(MUON))
        selected_branches.append("ScoutingMuon%s_phi"%(MUON))
        selected_branches.append("ScoutingMuon%s_charge"%(MUON))
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
    print(f"> Events with 2+ muons: {len(events)}") # 500_206 events Vtx TrgOR, # 356_664 events NoVtx TrgNoVtx
    # print(ak.fields(events[0])) # get list of fields again since i forgot
    # print(events[0]["nScoutingMuon%sDisplacedVertex"%(MUON)])

    ## Filter by having a displaced vertex reconstruction
    events = events[events["nScoutingMuon%sDisplacedVertex"%(MUON)] > 0]
    print(f"> Events with a displaced vertex reco: {len(events)}") 
    
    ## Print all charges as a test
    for i in tqdm(range(len(events))):
        print(events["ScoutingMuon%s_charge"%(MUON)][i])
        print(events["nScoutingMuon%s"%(MUON)])
        print("\n")