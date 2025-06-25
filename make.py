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

        # Selection
        selected_branches.append("ScoutingMuon%sDisplacedVertex_chi2"%(MUON))
        selected_branches.append("ScoutingMuon%sDisplacedVertex_ndof"%(MUON))

        selected_branches.append("ScoutingMuon%sDisplacedVertex_x"%(MUON))
        selected_branches.append("ScoutingMuon%sDisplacedVertex_y"%(MUON))
        selected_branches.append("ScoutingMuon%sDisplacedVertex_z"%(MUON))

        selected_branches.append("ScoutingMuon%s_trk_vx"%(MUON))
        selected_branches.append("ScoutingMuon%s_trk_vy"%(MUON))
        selected_branches.append("ScoutingMuon%s_trk_vz"%(MUON))

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
h_lxy_sidebands = hist.new.Reg(100, lxy_range[0], lxy_range[1], name="lxy_sidebands", label="lxy_sidebands").Double()
h_lxy_peak = hist.new.Reg(100, lxy_range[0], lxy_range[1], name="lxy_peak", label="lxy_peak").Double()

# nVtxIndxString = "ScoutingMuon%s_nScoutingMuon%sVtxIndx" % (MUON, MUON)
oVtxIndxString = "ScoutingMuon%s_oScoutingMuon%sVtxIndx" % (MUON, MUON)

# for i in tqdm(range(10)):
indexErrRejected = 0
rejected = 0
for i in tqdm(range(len(events))):  
# for i in tqdm(range(4080, 4101)):  
    # 4100 first instance of [[0, 1], [0, 1, 2], [0, 1, 2]]
    # 4097 is interesting because it has -1,1,1,1 and [0,0] -> matches first two (nice test of the code!)
    nMuons = events["nScoutingMuon%s"%(MUON)][i]

    print("~~~~~~~~~")
    print("Num Muons:", nMuons)
    print("Num Displaced Vertices:", events["nScoutingMuon%sDisplacedVertex"%(MUON)][i])
    print("ScoutingMuon%sVtxIndx_vtxIndx:"%(MUON), events["ScoutingMuon%sVtxIndx_vtxIndx"%(MUON)][i])
    print("ScoutingMuon%sDisplacedVertex_isValidVtx:"%(MUON), events["ScoutingMuon%sDisplacedVertex_isValidVtx"%(MUON)][i])
    print("Charges:", events["ScoutingMuon%s_charge"%(MUON)][i])

    vertexListByMuonIndex = [] # corresponds directly to, for example, charge entries
    oVtxIndxArray = events[oVtxIndxString][i]
    for n in range(nMuons):
        print("The %dth entry in oVtxIndxArray is:"%(n), oVtxIndxArray[n])
        current_offset = oVtxIndxArray[n]
        if n != (nMuons-1): # if NOT the final entry, use explicit list slicing
            # print("%d is not the final entry"%(n))
            next_offset = oVtxIndxArray[n+1]
            vertexSlice = events["ScoutingMuon%sVtxIndx_vtxIndx"%(MUON)][i][current_offset:next_offset] # exclusive
        else: # if it is the final entry, go to end of list
            vertexSlice = events["ScoutingMuon%sVtxIndx_vtxIndx"%(MUON)][i][current_offset::]
        print("Vertex Slice for %d is:"%(n), vertexSlice)
        vertexListByMuonIndex.append(vertexSlice)
    vertexArrayByMuonIndex = ak.Array(vertexListByMuonIndex)
    # flatten a bit using really weird syntax
    vertexArrayByMuonIndex = ak.Array([
        -999 if len(subarray) == 0 else subarray[0] if len(subarray) == 1 else subarray
        for subarray in vertexArrayByMuonIndex
    ])
    print("%d Vertices:"%(i), vertexArrayByMuonIndex)

    ### choice code for multiple vertices
    chi2s = events["ScoutingMuon%sDisplacedVertex_chi2"%(MUON)][i]
    ndofs = events["ScoutingMuon%sDisplacedVertex_ndof"%(MUON)][i]
    scores = chi2s/ndofs
    print("chi^2:", chi2s)
    print("ndof:", ndofs)
    print("scores", scores)

    '''This is about to get really ugly'''
    validVertices = np.sort(np.unique(ak.to_numpy(events["ScoutingMuon%sVtxIndx_vtxIndx"%(MUON)][i])))
    if len(validVertices) == 0:
        continue
    print("Valid Vertices:", validVertices) # [0, 1, 2] for example
    isGoodVertex = [score <= 10 for score in scores]

    print("Is a good vertex?", isGoodVertex)
    try: 
        validVertices = validVertices[isGoodVertex] # Returns list of vertices that pass the score threshold.
    except IndexError: #Event 58005, any more?
        print("Event %d threw an Index Error; Size of scores != size of vertex list"%(i))
        continue
        indexErrRejected +=1
    print("Now valid vertices:", validVertices)

    ###

    # stopgap
    try:
        maxVertexIdx = ak.max(vertexArrayByMuonIndex)
        print("maxVertexIdx value is:", maxVertexIdx)
        muonIndexArray = []
        for vertexIdx in range(maxVertexIdx + 1):
            # ak.where returns a tuple -> unpack
            muonIndexArray.append(ak.where(vertexArrayByMuonIndex == vertexIdx)[0]) 
        muonIndexArray = ak.Array(muonIndexArray)
            # print("Index Array:", indexArray)
            # indexArray is of the form [[all indices for vtx 0], [all indices for vtx 1], ...]
            # so if vtx 0 is the vertex for two muons (for instance) then it is just [[0, 1]]
        print("muonIndexArray:", muonIndexArray)

        # From https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideOfflinePrimaryVertexProduction:
        '''
            The primary vertex collection is sorted according to the sum of the Pt squared of the tracks associated to each vertex,
                such that the vertex with largest sum, likely to be the "signal" vertex, appears first.
                Some justification of this can be found in the Higgs note CMS-AN-11-129. 
        '''

        # get best PV according to sum(pt**2):
        pv_x = events["ScoutingPrimaryVertex_x"][i][0]
        pv_y = events["ScoutingPrimaryVertex_y"][i][0]
        # print("PV_x", pv_x)
        # print("PV_y", pv_y)

        for j, indices in enumerate(muonIndexArray): # j is the VERTEX INDEX
            # indices is the LIST OF MUON INDICES BELONGING TO EACH VERTEX
            sv_x = events["ScoutingMuon%sDisplacedVertex_x"%(MUON)][i][j]
            sv_y = events["ScoutingMuon%sDisplacedVertex_y"%(MUON)][i][j]
            print("SV_x %d is"%(j), sv_x)
            print("SV_y %d is"%(j), sv_y)

            dx = sv_x - pv_x
            dy = sv_y - pv_y
            lxy = np.sqrt(dx**2 + dy**2)
            h_lxy.fill(lxy=lxy)

            # Compute invariant mass stuff
            pt1, pt2 = events["ScoutingMuon%s_pt"%(MUON)][i][indices]
            eta1, eta2 = events["ScoutingMuon%s_eta"%(MUON)][i][indices]
            phi1, phi2 = events["ScoutingMuon%s_phi"%(MUON)][i][indices]

            mu_mass = 0.105
            pz1 = pt1 * np.sinh(eta1)
            pz2 = pt2 * np.sinh(eta2)
            px1 = pt1 * np.cos(phi1)
            px2 = pt2 * np.cos(phi2)
            py1 = pt1 * np.sin(phi1)
            py2 = pt2 * np.sin(phi2)
            E1  = np.sqrt(px1**2 + py1**2 + pz1**2 + mu_mass**2)
            E2  = np.sqrt(px2**2 + py2**2 + pz2**2 + mu_mass**2)
            
            invariant_mass = np.sqrt((E1 + E2)**2 - (px1 + px2)**2 - (py1 + py2)**2 - (pz1 + pz2)**2)
            print("Invariant mass in GeV is ", invariant_mass)
            if (invariant_mass >= 2.4 and invariant_mass < 2.9) or (invariant_mass > 3.3 and invariant_mass <= 3.8):
                h_lxy_sidebands.fill(lxy_sidebands=lxy)
            elif (invariant_mass >= 2.9 and invariant_mass <= 3.3):
                h_lxy_peak.fill(lxy_peak=lxy)
    except ValueError:
        rejected += 1
    print("Loop execution finished!")
with open (outdir+"/large_pickles/events%sLxyPickle.pkl"%(MUON), "wb") as pickleOut:
    pickle.dump(h_lxy, pickleOut)
    pickle.dump(lxy_range, pickleOut)
    pickle.dump(h_lxy_peak, pickleOut)
    pickle.dump(h_lxy_sidebands, pickleOut)

print("score mismatches:", indexErrRejected) 
print("many choice rejects: ", rejected)

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
'''

'''  # isGoodVertex = ak.Array([
    #     # False if score > 10 else True
    #     False if score > 1 else True # extremely aggressive cut for demonstration purposes
    #     for score in scores
    # ])
    # print("isGoodVertex:", isGoodVertex)
    # bestScore = ak.min(scores)
    # bestScoreIdx = ak.argmin(scores)
    # print("Minimum of scores", bestScore, bestScoreIdx)

    # print("TO-DO: FILTER VERTICES BY SCORE FOR SCALAR ARRAYS")
    # if ak.num(vertexArrayByMuonIndex, axis=1, highlevel=False) is None: # if a scalar array (good)
    #     for muon_idx, vertexNumber in enumerate(vertexArrayByMuonIndex):
    #         if not isGoodVertex[vertexNumber]:
    #             vertexArrayByMuonIndex[muon_idx] = -999
    # else:
    #     for muon_idx, subelement in enumerate(vertexArrayByMuonIndex):
    #         if isGoodVertex[muon_idx] and ak.any(subelement == bestScoreIdx):
    #             vertexArrayByMuonIndex[muon_idx] = bestScoreIdx
    #         else: # if it is not a good vertex or if it does not contain the best vertex
    #             vertexArrayByMuonIndex[muon_idx] = -999
'''