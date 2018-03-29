
# import ROOT in batch mode
import sys
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

# load FWLite C++ libraries
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.FWLiteEnabler.enable()

print 'Hello!'

# load FWlite python libraries
from DataFormats.FWLite import Handle, Events

muons, muonLabel = Handle("std::vector<pat::Muon>"), "slimmedMuons"
electrons, electronLabel = Handle("std::vector<pat::Electron>"), "slimmedElectrons"
photons, photonLabel = Handle("std::vector<pat::Photon>"), "slimmedPhotons"
taus, tauLabel = Handle("std::vector<pat::Tau>"), "slimmedTaus"
tauLabelb = "slimmedTausBoosted"
jets = Handle("std::vector<pat::Jet>")
fatjets, fatjetLabel = Handle("std::vector<pat::Jet>"), "slimmedJetsAK8"
fatgenjets, fatgenjetLabel = Handle("std::vector<reco::GenJet>"), "slimmedGenJetsAK8"
mets, metLabel = Handle("std::vector<pat::MET>"), "slimmedMETs"
vertices, vertexLabel = Handle("std::vector<reco::Vertex>"), "offlineSlimmedPrimaryVertices"
verticesScore = Handle("edm::ValueMap<float>")
seenIt = {} # list of things we've seen (so that we dump them in full only once)

# open file (you can use 'edmFileUtil -d /store/whatever.root' to get the physical file name)
events = Events("root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/20000/001C74A0-B4D6-E711-BD4B-FA163EB4F61D.root")

print 'Events read...'

for iev,event in enumerate(events):

    #print 'iev = ', iev
    #if iev >= 10: break
    event.getByLabel(muonLabel, muons)
    event.getByLabel(electronLabel, electrons)
    event.getByLabel(photonLabel, photons)
    event.getByLabel(tauLabel, taus)
    event.getByLabel(fatjetLabel, fatjets)
    event.getByLabel(fatgenjetLabel, fatgenjets)
    event.getByLabel(metLabel, mets)
    event.getByLabel(vertexLabel, vertices)
    event.getByLabel(vertexLabel, verticesScore)

    #

    # Jets (AK4, CHS and Puppi)
    if False : 
        for jetLabel, algo in ("slimmedJets", "CHS"), ("slimmedJetsPuppi", "PUPPI"):
            event.getByLabel(jetLabel, jets)
            for i,j in enumerate(jets.product()):
                if j.pt() < 20: continue
                print "jet %s %3d: pt %5.1f (raw pt %5.1f, matched-calojet pt %5.1f), eta %+4.2f, btag CSVIVFv2 %.3f, CMVAv2 %.3f, pileup mva disc %+.2f" % (
                    algo, i, j.pt(), j.pt()*j.jecFactor('Uncorrected'), j.userFloat("caloJetMap:pt") if algo == "CHS" else -99.0, j.eta(), max(0,j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")), max(0,j.bDiscriminator("pfCombinedMVAV2BJetTags")), j.userFloat("pileupJetId:fullDiscriminant") if algo == "CHS" else -99)
                if 'jetAk4'+algo not in seenIt:
                    constituents = [ j.daughter(i2) for i2 in xrange(j.numberOfDaughters()) ]
                    constituents.sort(key = lambda c:c.pt(), reverse=True)
                    for i2, cand in enumerate(constituents):
                        if i2 > 12:
                                print "         ....."
                                break
                        print "         constituent %3d: pt %6.2f, dz(pv) %+.3f, pdgId %+3d, hcal energy fraction %.2f, puppi weight %.3f " % (i2,cand.pt(),cand.dz(PV.position()),cand.pdgId(),cand.hcalFraction(),cand.puppiWeight())
                    print "   btag discriminators:"
                    for btag in j.getPairDiscri():
                        print  "\t%s %s" % (btag.first, btag.second)
                    print "   userFloats:"
                    for ufl in j.userFloatNames():
                        print  "\t%s %s" % (ufl, j.userFloat(ufl))
                    seenIt['jetAk4'+algo] = True

    # Fat AK8 Jets
    if len(fatjets.product()) < 1 or fatjets.product()[0].pt() < 200:
        continue

    print "\nEvent %d: run %6d, lumi %4d, event %12d" % (iev,event.eventAuxiliary().run(), event.eventAuxiliary().luminosityBlock(),event.eventAuxiliary().event())
    
    for i,j in enumerate(fatgenjets.product()):
        print "genAK8 %3d: pt %5.1f " % (
            i, j.pt(), )
        
    for i,j in enumerate(fatjets.product()):
        genpt = 0.0
        if j.genJet() != 0 and j.genJet() != None :
            genpt = j.genJet().pt()
        print "jetAK8 %3d: pt %5.1f (raw pt %5.1f), eta %+4.2f, mass %5.1f ungroomed, %5.4f softdrop, N2 = %5.2f, N3 = %5.2f., genpt = %5.2f " % (
            i, j.pt(), j.pt()*j.jecFactor('Uncorrected'), j.eta(), j.mass(), j.userFloat('ak8PFJetsPuppiSoftDropMass'), j.userFloat('ak8PFJetsPuppiSoftDropValueMap:nb1AK8PuppiSoftDropN2'), j.userFloat('ak8PFJetsPuppiSoftDropValueMap:nb1AK8PuppiSoftDropN3'), genpt)
        # To get the constituents of the AK8 jets, you have to loop over all of the
        # daughters recursively. To save space, the first two constituents are actually
        # the Soft Drop SUBJETS, which will then point to their daughters.
        # The remaining constituents are those constituents removed by soft drop but
        # still in the AK8 jet.

        if True : 
            if True:
                constituents = []
                constituentsGroomed = []
                for ida in xrange( j.numberOfDaughters() ) :
                    cand = j.daughter(ida)
                    if cand.numberOfDaughters() == 0 :
                        constituents.append( cand )
                    else :
                        for jda in xrange( cand.numberOfDaughters() ) :
                            cand2 = cand.daughter(jda)
                            constituents.append( cand2 )
                            val = ROOT.TLorentzVector( cand2.px(), cand2.py(), cand2.pz(), cand2.energy())
                            constituentsGroomed.append( val * cand2.puppiWeight() )
                constituents.sort(key = lambda c:c.pt(), reverse=True)
                #groomedJet = copy.copy( constituentsGroomed[0] )
                if len(constituentsGroomed) > 0:
                    groomedJet = sum( constituentsGroomed , ROOT.TLorentzVector())
                for i2, cand in enumerate(constituents):
                    if i2 >4:
                                print "         ....."
                                break
                    print "         constituent %3d: pt %6.2f, pdgId %+3d, #dau %+3d" % (i2,cand.pt(),cand.pdgId(), cand.numberOfDaughters())
                if groomedJet : 
                    print "         reclustered    : pt %6.2f, mass %6.4f" % (groomedJet.Perp(),groomedJet.M())
                print "   btag discriminators:"
                for btag in j.getPairDiscri():
                    print  "\t%s %s" % (btag.first, btag.second)
                print "   userFloats:"
                for ufl in j.userFloatNames():
                    print  "\t%s %s" % (ufl, j.userFloat(ufl))
                seenIt['jetAk8'] = True
            # Print Subjets
            if  j.hasSubjets('SoftDropPuppi'):
                wSubjets = j.subjets('SoftDropPuppi')
                for iw,wsub in enumerate( wSubjets ) :
                    print "   w subjet %3d: pt %5.1f (raw pt %5.1f), eta %+4.2f, mass %5.1f" % (
                        iw, wsub.pt(), wsub.pt()*wsub.jecFactor('Uncorrected'), wsub.eta(), wsub.mass()
                        )
                    print "   \tbtag discriminators:"
                    for btag in wsub.getPairDiscri():
                        print  "\t\t%s %s" % (btag.first, btag.second)
                    print "   \tuserFloats:"
                    for ufl in wsub.userFloatNames():
                        print  "\t\t%s %s" % (ufl, wsub.userFloat(ufl))
                    seenIt['jetAk8SD'] = True
