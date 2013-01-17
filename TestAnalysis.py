from DataFormats.FWLite import Events, Handle
from ROOT import TH1D, TCanvas, std

inputfiles = 'testdimuon.root'

events = Events(inputfiles)

chicand_h = Handle('vector<pat::CompositeCandidate>')
chicand_l = ('chiCandProducer','chicand')

piZero_h = Handle('std::vector<std::vector<float> >')
piZero_l = ('chiCandProducer','piZeroRejectCand')

histo = TH1D("chib","ChiB invariant mass", 100, 9.6,11.1)
#invmcVect = std.vector('float')(1000)

for event in events:

	event.getByLabel(chicand_l, chicand_h)
	event.getByLabel(piZero_l,piZero_h)

	if (not chicand_h.isValid() or not piZero_h.isValid()):
        	continue

	chiColl = chicand_h.product()
	piZero = piZero_h.product()
	select = 1
	chi_index = 0
	for cand in chiColl:
#		if( cand.hasUserData("invmc") ):
#			invmcVect = cand.invmc()
		photon_index = 0
		dimuon_index = 1
		if( cand.daughter(1).name() == "photon" ):
			photon_index = 1
			dimuon_index = 0
		if (cand.daughter(photon_index).pt() < 1.0 ):
			select = 0
		if( cand.daughter(dimuon_index).mass() < 8.6 ):		
			select = 0
		if( cand.daughter(dimuon_index).mass() > 11.4 ):		
			select = 0
		if( abs(cand.daughter(dimuon_index).rapidity() ) > 1.0):
			select = 0
		if( abs(cand.daughter(photon_index).eta() ) > 1.0):
			select = 0
		# dz selection is applied at producer stage with a value of 0.5

		# deltamass selection is applied at producer stage with a value of 0.5

		for pi in piZero[chi_index]:
			if( pi > 0.115 and pi < 0.155):
				select = 0
		if(select):
			histo.Fill( cand.mass() - cand.daughter(dimuon_index).mass() + 9.46 )
		chi_index += 1

canvas = TCanvas()
histo.Draw("E")
canvas.SaveAs("testinvm.root")
canvas.SaveAs("testinvm.png")

