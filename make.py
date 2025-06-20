import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
import mplhep
import hist
import glob
import pickle
from tqdm import tqdm

indir = "/eos/user/f/fernance/DST-Muons/samples"
# files = glob.glob(f"{indir}/*.root")
file = uproot.open(indir+"/e0e9c734-dead-4179-b699-f6b75d0e4502.root")
print(file)
for event in file["Events"]:
    print(event)