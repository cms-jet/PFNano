#!/usr/bin/env python

import os
import sys
import argparse

import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
hep.set_style("CMS")

check_files = [
    'nano_data2016.root',
    'nano_data2017.root',
    'nano_data2018abc.root',
    'nano_data2018d.root',
    'nano_mc2016.root',
    'nano_mc2017.root',
    'nano_mc2018.root',
]

check_hists = [
    ("FatJet_btagDDBvLV2", 0, 1, 30),
    ("FatJet_btagDDCvLV2", 0, 1, 30),
    ("FatJet_btagDDCvBV2", 0, 1, 30),
    ("FatJet_btagDDBvL", 0, 1, 30),
    ("FatJet_btagDDCvL", 0, 1, 30),
    ("FatJet_btagDDCvB", 0, 1, 30),
    ("FatJet_particleNetMD_Xbb", 0, 1, 30),
    ("FatJet_particleNetMD_Xcc", 0, 1, 30),
    ("FatJet_particleNetMD_QCD", 0, 1, 30),
    ("FatJet_pt", 300, 800, 30),
    ("FatJet_eta", -2.5, 2.5, 30),
    ("FatJet_mass", 300, 800, 30),
    ("FatJet_msoftdrop", 300, 800, 30),   
]

up_files = [uproot.open(f) for f in check_files]

try:
    os.mkdir('valplots')
except:
    pass

print("Plotting:")
for hist_spec in check_hists:
    fig, ax = plt.subplots()
    print(f"  {hist_spec[0]}")
    for fname, rfile in zip(check_files, up_files):
        pts = rfile['Events'].array(hist_spec[0]).flatten()
        ax.hist(pts, bins=np.linspace(hist_spec[1], hist_spec[2], hist_spec[3]), label=fname.split(".")[0].split("_")[1], density=True, histtype='step')
    hep.cms.text("Simulation Preliminary")
    ax.set_xlabel(hist_spec[0])
    ax.legend()
    fig.savefig(f'valplots/{hist_spec[0]}.png')

print("Checking Coffea/NanoEvents compatibility:")

import coffea
import uproot
import awkward1 as ak
import numpy as np
from coffea.nanoevents import NanoEventsFactory

for fn in check_files: 
    print(f"  {fn}")
    def nano_evts(fname):
        factory = NanoEventsFactory.from_file(
            fname,
            entry_start=0, entry_stop=10000,
            metadata={"dataset": ""},
        )
        return factory.events()

    evts = nano_evts(fn)
    print("  Subjets pt ordered (should be false):", ak.all(evts.SubJet.pt[:, :-1] - evts.SubJet.pt[:, 1:] >= 0))
    print("  SDmass from subjets hist:")
    print(np.histogram(ak.to_numpy(ak.flatten(evts.FatJet.subjets.sum().mass, None)), bins=np.linspace(40, 200, 20))[0])



    
