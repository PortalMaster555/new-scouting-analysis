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

outdir = "/afs/cern.ch/user/d/dspringb/private/CMSSW_15_0_5/src/new-scouting-analysis"
with open (outdir+"/large_pickles/events%sLxyPickle.pkl"%(MUON), "rb") as pickleIn:
    h_lxy = pickle.load(pickleIn)
    lxy_range = pickle.load(pickleIn)

plt.style.use(mplhep.style.CMS)
fig, ax = plt.subplots(figsize=(10,8))
h_lxy.plot(ax=ax)
ax.set_xlim(lxy_range)
bin_values = h_lxy.values()

# nonzero bin values only because many bins are 0
nonzero_vals = bin_values[bin_values > 0]
yMin = 10**np.floor(np.log10(nonzero_vals.min()))
yMax = 10**np.ceil(np.log10(nonzero_vals.max()))

ax.set_ylim(yMin, yMax)
ax.set_xlabel("lxy")
ax.set_ylabel("Number of events")
# ax.set_xscale("log")
ax.set_yscale("log")
fig.savefig("img/hist_lxy_temp_%s.png"%(MUON), dpi=300)

