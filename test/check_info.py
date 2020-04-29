import uproot

f = uproot.open("nano106X_on_mini106X_2017_mc_NANO.root")
print('File contents:')
print('    ', f.keys())

print('File/Events contents:')
for key in f['Events'].keys():
    #if key.startswith("Jet_") or key.startswith("FatJet_") or key.startswith("SubJet_"):
     if "Proba" in key:
        print("   ", key)

interest = [
    "Jet_Proba",
    "FatJet_Proba",
    "SubJet_Proba",
    "Jet_nBHadrons",
    "FatJet_nBHadrons",
    "SubJet_nBHadrons",
]    

for var in interest:
    print(var)
    print(f['Events'].arrays([var])[var])

# for ev in f['Events'].arrays(['Jet_Proba'])['Jet_Proba']:
#     print(ev)
# print("Jet")
# 
