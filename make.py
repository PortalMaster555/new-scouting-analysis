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

np.set_printoptions(precision=8, suppress=False)

try:
    if sys.argv[1].lower() == "--vtx":
        MUON = "Vtx"
        TRIGGER = ["DST_PFScouting_DoubleMuonVtx"]
        # pickleString = "zoom_2025/histos_2025_TrgVtx_Vtx.pkl"
    elif sys.argv[1].lower() == "--novtx":
        MUON = "NoVtx"
        TRIGGER = ["DST_PFScouting_DoubleMuonNoVtx"]
        # pickleString = "zoom_2025/histos_2025_TrgNoVtx_NoVtx.pkl"
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
    print("Using stored data from "+ outdir + "/large_pickles/events%sPickle.pkl"%(MUON))
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

        selected_branches.append("nScoutingMuon%sVtxIndx"%(MUON))
        selected_branches.append("ScoutingMuon%sVtxIndx_vtxIndx"%(MUON))
        selected_branches.append("ScoutingMuon%sDisplacedVertex_isValidVtx"%(MUON))

        # Access
        selected_branches.append("ScoutingMuon%s_nScoutingMuon%sVtxIndx" % (MUON, MUON))
        selected_branches.append("ScoutingMuon%s_oScoutingMuon%sVtxIndx" % (MUON, MUON))

        selected_branches.append("ScoutingMuonNoVtxDisplacedVertex_x")
        selected_branches.append("ScoutingMuonNoVtxDisplacedVertex_y")
        selected_branches.append("ScoutingMuonNoVtxDisplacedVertex_z")

        selected_branches.append("ScoutingMuonNoVtx_trk_vx")
        selected_branches.append("ScoutingMuonNoVtx_trk_vy")
        selected_branches.append("ScoutingMuonNoVtx_trk_vz")

        selected_branches.append("nScoutingPrimaryVertex")
        selected_branches.append("ScoutingPrimaryVertex_x")
        selected_branches.append("ScoutingPrimaryVertex_y")
        selected_branches.append("ScoutingPrimaryVertex_z")

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
            with open (outdir+"/large_pickles/events%sPickle.pkl"%(MUON), "wb") as pickleOut:
                pickle.dump(events, pickleOut)
            pickleIndex += 1
    else: # if not rereading the original unfiltered file
        with open (outdir+"/large_pickles/events%sPickle.pkl"%(MUON), "rb") as pickleIn:
            events = pickle.load(pickleIn)
        print(f"> Opened filtered sample with {len(events)} events") # 500_206 events Vtx TrgOR, # 356_664 events NoVtx TrgNoVtx
###################################################


# Remove muons where |eta| is greater than/eq to 2.4
print("*Allow |eta| < 2.4*")
# mask_pt = events["ScoutingMuonNoVtx_pt"] >= 5
mask_eta = (abs(events["ScoutingMuonNoVtx_eta"]) < 2.4)
# combined_mask = mask_pt & mask_eta
combined_mask = mask_eta
# print(combined_mask)
muon_fields = [
"_pt", "_eta", "_phi", "_charge",
"_trk_vx", "_trk_vy", "_trk_vz"]

for field in muon_fields:
    key = f"ScoutingMuonNoVtx{field}"
    # print(len(events[key]))
    events[key] = events[key][combined_mask]
events["nScoutingMuonNoVtx"] = ak.num(events["ScoutingMuonNoVtx_pt"])
print("> Filter by eta successful.")


## Filter by number of muons
print("*Allow nMuons >= 2*")
events = events[ events["nScoutingMuon%s"%(MUON)] > 1]
print(f"> Significant events with 2+ muons: {len(events)}") 
# print(ak.fields(events[0])) # get list of fields again since i forgot
# print(events[0]["nScoutingMuon%sDisplacedVertex"%(MUON)])


## Filter by having a displaced vertex reconstruction
print("*Allow nDisplacedVertex >= 1*")
events = events[events["nScoutingMuon%sDisplacedVertex"%(MUON)] >= 1]
print(f"> Events with a displaced vertex reco: {len(events)}") 
# events = events[events["ScoutingMuon%sDisplacedVertex_isValidVtx"%(MUON)] == True]
# print(f"> Events with a valid displaced vertex reco: {len(events)}") 

## Filter by requiring at least one pair of opposite charges
print("*Allow opposite-charge-containing events*")
charges = events[f"ScoutingMuon{MUON}_charge"]
has_pos = ak.any(charges > 0, axis=1)
has_neg = ak.any(charges < 0, axis=1)
keep_mask = has_pos & has_neg
events = events[keep_mask]
print(f"> Events with at least one pair of opposite charges: {len(events)}") 

# Addtl. pre-filtering goes here.

### BEGIN POST-FILTER ITERATION ###
lxy_range = (0, 7.5)
h_lxy = hist.new.Reg(100, lxy_range[0], lxy_range[1], name="lxy", label="lxy").Double()

# nVtxIndxString = "ScoutingMuon%s_nScoutingMuon%sVtxIndx" % (MUON, MUON)
oVtxIndxString = "ScoutingMuon%s_oScoutingMuon%sVtxIndx" % (MUON, MUON)

# for i in tqdm(range(10)):
for i in tqdm(range(len(events))):  
    print("~~~~~~~~~")
    nMuons = events["nScoutingMuon%s"%(MUON)][i]
    print("Num Muons:", nMuons)
    print("Num Displaced Vertices:", events["nScoutingMuon%sDisplacedVertex"%(MUON)][i])
    # print("nScoutingMuon%s_VtxIndx:"%(MUON), events["nScoutingMuon%sVtxIndx"%(MUON)][i])
    print("ScoutingMuon%sVtxIndx_vtxIndx:"%(MUON), events["ScoutingMuon%sVtxIndx_vtxIndx"%(MUON)][i])
    print("ScoutingMuon%sDisplacedVertex_isValidVtx:"%(MUON), events["ScoutingMuon%sDisplacedVertex_isValidVtx"%(MUON)][i])

    print("Charges:", events["ScoutingMuon%s_charge"%(MUON)][i])

    # print(nVtxIndxString, i, events[nVtxIndxString][i])
    # print(oVtxIndxString, i, events[oVtxIndxString][i])
    # print("ScoutingMuon%sVtxIndx_vtxIndx:"%(MUON), events["ScoutingMuon%sVtxIndx_vtxIndx"%(MUON)][i])

    # print("* Begin *")
    vertexListByMuonIndex = [] # corresponds directly to, for example, charge entries
    oVtxIndxArray = events[oVtxIndxString][i]
    for n in range(nMuons):
        # print("The %dth entry in oVtxIndxArray is:"%(n), oVtxIndxArray[n])
        current_offset = oVtxIndxArray[n]
        if n != (nMuons-1): # if NOT the final entry, use explicit list slicing
            # print("%d is not the final entry"%(n))
            next_offset = oVtxIndxArray[n+1]
            vertexSlice = events["ScoutingMuon%sVtxIndx_vtxIndx"%(MUON)][i][current_offset:next_offset] # exclusive
        else: # if it is the final entry, go to end of list
            vertexSlice = events["ScoutingMuon%sVtxIndx_vtxIndx"%(MUON)][i][current_offset::]
        # print("Vertex Slice for %d is:"%(n), vertexSlice)
        vertexListByMuonIndex.append(vertexSlice)
    # print("*  End  *")
    vertexArrayByMuonIndex = ak.Array(vertexListByMuonIndex)
    print("Vertices:", vertexArrayByMuonIndex)


##################################################
'''

            lowestEtas = events[f"ScoutingMuon{MUON}_eta"][i][indices]
            lowestPhis = events[f"ScoutingMuon{MUON}_phi"][i][indices]
            eta1, eta2 = lowestEtas[0], lowestEtas[1]
            phi1, phi2 = lowestPhis[0], lowestPhis[1]
            deltaEta = np.abs(eta1 - eta2)
            deltaPhi = np.abs(phi1 - phi2)
            deltaPhi = np.where(deltaPhi > np.pi, 2 * np.pi - deltaPhi, deltaPhi)
            deltaR = np.sqrt(deltaEta**2 + deltaPhi**2)
            if deltaR > 0.2: # remove duplicate muons
                print("deltaR > 0.2.")
                print("TO-DO: FILTER THE CORRECT PV FROM THE LIST OF SIZE:", events["nScoutingPrimaryVertex"][i])
                pv_x = events["ScoutingPrimaryVertex_x"][i][0]
                pv_y = events["ScoutingPrimaryVertex_y"][i][0]
                print("TO-DO: GET CORRECT DISPLACED VERTEX FROM LIST")
                sv_x = events["ScoutingMuon%sDisplacedVertex_x"%(MUON)][i][0]
                sv_y = events["ScoutingMuon%sDisplacedVertex_y"%(MUON)][i][0]
                
                dx = sv_x - pv_x
                dy = sv_y - pv_y
                lxy = np.sqrt(dx**2 + dy**2)
                h_lxy.fill(lxy=lxy)

with open (outdir+"/large_pickles/events%sLxyPickle.pkl"%(MUON), "wb") as pickleOut:
    pickle.dump(h_lxy, pickleOut)
    pickle.dump(lxy_range, pickleOut)
'''