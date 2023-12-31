from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from ROOT import TH1D, TFile, TLorentzVector, TMath, TTree
from math import *
from optparse import OptionParser
from array import array
import os
import fnmatch

#########################

inFile ='/collab/project/snowmass21/data/muonc/fmeloni/LegacyProductions/before29Jul23/DataMuC_MuColl_v1/electronGun/reco/electronGun_reco_0.slcio'

parser = OptionParser()
parser.add_option('-i','--inputFile', help='--inFile electronGun_reco_0.slcio',
                  type=str, default='electronGun_reco_0.slcio')
parser.add_option('-o', '--outDir', help='--outDir ./',
                  type=str, default='./')
(options, args) = parser.parse_args()

# declare histograms
arrBins_R = array('d', (0., 10., 20., 31., 51., 74., 102.,
                        127., 150., 200., 250., 340., 450., 554.))
arrBins_pT = array('d', (0., 0.5, 1., 1.5, 2., 2.5, 3.,
                         3.5, 4., 5., 6., 7., 8., 10., 20., 30., 50.))
# arrBins_pT = array('d', (0., 0.5, 1., 1.5, 2., 2.5, 3.,
#                         3.5, 4., 5., 6., 7., 8., 10., 20., 30., 50., 75., 100., 250., 500., 1000., 1500.))
arrBins_theta = array('d', (30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180.,
                            90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180.))

h_truth_Rprod = TH1D('truth_Rprod', 'truth_Rprod',
                     len(arrBins_R)-1, arrBins_R)  # mm
h_truth_pT = TH1D('truth_pT', 'truth_pT', len(arrBins_pT)-1, arrBins_pT)
h_truth_theta = TH1D('truth_theta', 'truth_theta',
                     len(arrBins_theta)-1, arrBins_theta)

h_trk_pT = TH1D('trk_pT', 'trk_pT', len(arrBins_pT)-1, arrBins_pT)
h_trk_theta = TH1D('trk_theta', 'trk_theta',
                   len(arrBins_theta)-1, arrBins_theta)
h_trk_Rprod = TH1D('trk_Rprod',
                   'trk_Rprod', len(arrBins_R)-1, arrBins_R)  # mm

h_ele_pT = TH1D('ele_pT', 'ele_pT', len(arrBins_pT)-1, arrBins_pT)
h_ele_theta = TH1D('ele_theta', 'ele_theta',
                   len(arrBins_theta)-1, arrBins_theta)
h_ele_Rprod = TH1D('ele_Rprod',
                   'ele_Rprod', len(arrBins_R)-1, arrBins_R)  # mm

histos_list = [h_truth_Rprod, h_truth_pT, h_truth_theta,
               h_trk_pT, h_trk_theta, h_trk_Rprod,
               h_ele_pT, h_ele_theta, h_ele_Rprod]

for histo in histos_list:
    histo.SetDirectory(0)


# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(inFile)

Bfield = 3.56  # T
# loop over all events in the file
for ievt, event in enumerate(reader):

    pfoCollection = event.getCollection('PandoraPFOs')
    trkCollection = event.getCollection('SiTracks_Refitted') 
    mcpCollection = event.getCollection('MCParticle')
    
   # relationCollection = event.getCollection('MCParticle,SiTracks_Refitted')
   # relation = UTIL.LCRelationNavigator(relationCollection)


    for mcp in mcpCollection:

        charge = mcp.getCharge()
        status = mcp.getGeneratorStatus()

        if fabs(charge) > 0:
            if fabs(mcp.getPDG()) == 11:
                vx = mcp.getVertex()
                rprod = sqrt(vx[0]*vx[0]+vx[1]*vx[1])
                dp3 = mcp.getMomentum()
                tlv = TLorentzVector()
                tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], mcp.getEnergy())
                
                goodtheta = False
                
                if tlv.Theta() > 30.*TMath.Pi()/180. and tlv.Theta() < 150.*TMath.Pi()/180.:
                    goodtheta = True

                if tlv.Perp() > 1 and not mcp.isDecayedInTracker() and goodtheta:
                    h_truth_Rprod.Fill(rprod)
                    h_truth_pT.Fill(tlv.Perp())
                    h_truth_theta.Fill(tlv.Theta())
                    tracks = relation.getRelatedToObjects(mcp)

                    if len(tracks) > 0:
                        track = tracks[0]
                        h_trk_Rprod.Fill(rprod)
                        h_trk_pT.Fill(tlv.Perp())
                        h_trk_theta.Fill(tlv.Theta())

                    for pfo in pfoCollection:
                        if fabs(pfo.getType()) == 11:
                            dp3 = pfo.getMomentum()
                            tlv_pfo = TLorentzVector()
                            tlv_pfo.SetPxPyPzE(
                                dp3[0], dp3[1], dp3[2], pfo.getEnergy())

                            pfotracks = pfo.getTracks()
                            if len(pfotracks) > 0:
                                trk = pfotracks[0]
                                print(trk)
              
                            if track == trk:
                                h_ele_Rprod.Fill(rprod)
                                h_ele_pT.Fill(tlv.Perp())
                                h_ele_theta.Fill(tlv.Theta())
reader.close()

# write histograms
output_file = TFile(options.outDir + "ntup_elePFO0.root", 'RECREATE')
for histo in histos_list:
    histo.Write()
output_file.Close()
