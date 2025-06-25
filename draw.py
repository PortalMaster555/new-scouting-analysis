import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
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
    h_lxy_peak = pickle.load(pickleIn)
    h_lxy_sidebands = pickle.load(pickleIn)



# Fit an exponential to the two curves
''' Definition from other file:
h_lxy_sidebands = hist.new.Reg(100, lxy_range[0], lxy_range[1], name="lxy_sidebands", label="lxy_sidebands").Double()
h_lxy_peak = hist.new.Reg(100, lxy_range[0], lxy_range[1], name="lxy_peak", label="lxy_peak").Double()
'''
def func(x, a, b):
    return a * x**(-b)

peak_bin_values = h_lxy_peak.values()
peak_bin_centers = h_lxy_peak.axes[0].centers
sidebands_bin_values = h_lxy_sidebands.values()
sidebands_bin_centers = h_lxy_sidebands.axes[0].centers

print(peak_bin_values)
print(peak_bin_centers)
print(sidebands_bin_values)


# fit in log space using np
center_min = 1
center_max = 20
logx = np.log(peak_bin_centers[center_min:center_max])
logy = np.log(peak_bin_values[center_min:center_max])
slope, intercept = np.polyfit(logx, logy, 1)
pk_a = np.exp(intercept)
pk_b = -slope
x = np.linspace(peak_bin_centers[center_min], peak_bin_centers[center_max-1], num = center_max-center_min) # where power law is strongest fit
pk_y = func(x, pk_a, pk_b)

logx = np.log(sidebands_bin_centers[center_min:center_max])
logy = np.log(sidebands_bin_values[center_min:center_max])
slope, intercept = np.polyfit(logx, logy, 1)
sb_a = np.exp(intercept)
sb_b = -slope
x = np.linspace(sidebands_bin_centers[center_min], sidebands_bin_centers[center_max-1], num = center_max-center_min) # where power law is strongest fit
sb_y = func(x, sb_a, sb_b)

plt.style.use(hep.style.CMS)
fig, ax = plt.subplots(figsize=(10,8))
hep.cms.label("Preliminary", data=True, year='2025', com='13.6', ax=ax, loc=2)
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


##### 

print("Plotting dual plot...")
# Dual peak/sidebands plot
plt.style.use(hep.style.CMS)
fig, ax = plt.subplots(figsize=(10,8))
hep.cms.label("Preliminary", data=True, year='2025', com='13.6', ax=ax, loc=2)
# h_lxy.plot(ax=ax, label="Full lxy")

h_lxy_peak.plot(ax=ax, label="Peak mass lxy")
h_lxy_sidebands.plot(ax=ax, label="Sidebands mass lxy")

plt.plot(x, pk_y, label="Peak curvefit", color="blue")
plt.plot(x, sb_y, label="Sidebands curvefit", color="red")

ax.legend(loc='center right', fontsize = 16, frameon = False, ncol=1)
ax.set_xlim(1e-2, 7)
# ax.set_xlim(lxy_range)
bin_values = h_lxy_peak.values()
# nonzero bin values only because many bins are 0
nonzero_vals = bin_values[bin_values > 0]
yMin = 1
yMax = 10**np.ceil(np.log10(nonzero_vals.max()))

ax.set_ylim(yMin, yMax)
ax.set_xlabel("lxy")
ax.set_ylabel("Number of events")

ax.set_xscale("linear")
ax.set_yscale("log")

textstr1 = (
    "Peak " + r"$ax^{-b}$" + "\n"
    r"$a=$" + f"{pk_a}\n"
    r"$b=$" + f"{pk_b}\n"
    "\nSidebands " + r"$ax^{-b}$" + "\n"
    r"$a=$" + f"{sb_a}\n"
    r"$b=$" + f"{sb_b}\n"
)
txt1 = ax.text(
    0.70, 0.98, textstr1,
    ha='right', va='top',
    fontsize=12,
    transform=ax.transAxes,
    bbox=dict(boxstyle='square,pad=0.3', facecolor='white', edgecolor='black', linewidth=1.5),
)

fig.savefig("img/hist_lxy_pk_side_%s.png"%(MUON), dpi=300)
