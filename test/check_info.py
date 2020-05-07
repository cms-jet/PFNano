import uproot

f = uproot.open("nano106X_on_mini106X_2017_mc_NANO.root")
print('File contents:')
print('    ', f.keys())

print('File/Events contents:')
for key in f['Events'].keys():
    #if key.startswith("Jet_") or key.startswith("FatJet_") or key.startswith("SubJet_"):
    pass
    if "Jet" in key:
        print("   ", key)

interest = [
    "Jet_Proba",
    "FatJet_Proba",
    "SubJet_Proba",
    "Jet_nBHadrons",
    "FatJet_nBHadrons",
    "SubJet_nBHadrons",
    "nFatJet",
    "FatJet_Proba",
    "FatJetSV_jetIdx",
    "FatJetSV_mass",
    #"FatJetSV_phirel",
    "nJetSV",
    "JetSV_jetIdx",
    "JetSV_mass",
    #"JetSV_phirel",
    "SV_mass",
    'FatJetSV_mass',
    'FatJetSV_pt',
    'FatJetSV_ntracks',
    'FatJetSV_chi2',
    'FatJetSV_normchi2',
    'FatJetSV_dxy',
    'FatJetSV_dxysig',
    'FatJetSV_d3d',
    'FatJetSV_d3dsig',
    'FatJetSV_costhetasvpv',
    'FatJetSV_phirel',
    'FatJetSV_ptrel',
    'FatJetSV_deltaR',
    'FatJetSV_enration',
    'FatJetSV_jetIdx',
]    

for var in interest:
    print(var)
    print(f['Events'].arrays([var])[var])

# for ev in f['Events'].arrays(['Jet_Proba'])['Jet_Proba']:
#     print(ev)
# print("Jet")
# 
