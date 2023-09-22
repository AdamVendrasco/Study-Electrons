from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from ROOT import *
from math import *
from array import array
import os
import glob
#########################

#declare array
arrBins_R = array('d', (0., 10., 20., 31., 51., 74., 102.,
                        127., 150., 200., 250., 340., 450., 554.))
arrBins_pT = array('d', (0., 0.5, 1., 1.5, 2., 2.5, 3.,
                         3.5, 4., 5., 6., 7., 8., 10., 20., 30., 50.))
arrBins_pT = array('d', (0., 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 6., 7., 8., 10., 20., 30., 50., 75., 100., 250., 500., 1000., 1500.))
arrBins_theta = array('d', (30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180.,
                            90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180.))
arr_all_p= array('d', (0,5,10,15,20,25,30,35,100,200,300,1000,1100,1500,2000,2100,2200))


# declare histograms
h_truth_Rprod = TH1D('truth_Rprod', 'truth_Rprod', len(arrBins_R)-1, arrBins_R)  # mm
h_truth_pT = TH1D('truth_pT', 'truth_pT', len(arrBins_pT)-1, arrBins_pT)
h_truth_theta = TH1D('truth_theta', 'truth_theta',len(arrBins_theta)-1, arrBins_theta)

h_trk_pT = TH1D('trk_pT', 'trk_pT', len(arrBins_pT)-1, arrBins_pT)
h_trk_theta = TH1D('trk_theta','trk_theta',len(arrBins_theta)-1, arrBins_theta)
h_trk_Rprod = TH1D('trk_Rprod','trk_Rprod', len(arrBins_R)-1, arrBins_R)  # mm

h_ele_pT = TH1D('ele_pT', 'ele_pT', len(arrBins_pT)-1, arrBins_pT)
h_ele_theta = TH1D('ele_theta', 'ele_theta',len(arrBins_theta)-1, arrBins_theta)
h_ele_Rprod = TH1D('ele_Rprod','ele_Rprod', len(arrBins_R)-1, arrBins_R)  # mm



#No BIB
h_all_p = TH1D('h_all_p','All Particles PDGID, No BIB',300,0,2200)

#Neutrons
h_pfo_neut_pT= TH1D("h_neut_pT","Neutron_pT",100,0,0.05)
h_pfo_neut_eta =TH1D("h_neut_eta","Neutron_eta",200,0,0.05)
h_pfo_neut_phi=TH1D("h_neut_phi","Neutron_phi",200,0,0.05)
h_pfo_neut_eta_v_phi= TH2D("h_pfo_neut_eta_v_phi","h_pfo_neut_eta_v_phi",20,-4,4,20,-3.5,3.5)
 
#Electrons
h_pfo_ele_theta = TH1D("h_pfo_ele_theta","h_pfo_ele_theta",100,0,10)
h_pfo_ele_eta = TH1D("h_pfo_ele_eta","h_pfo_ele_eta",100,-4,4)
h_pfo_ele_pT = TH1D("h_pfo_ele_pT","h_pfo_ele_pT",100,-10,1200)
h_pfo_ele_phi=TH1D("h_pfo_ele_phi","h_pfo_ele_phi",100,-4,4)
h_pfo_ele_eta_v_phi= TH2D("h_pfo_ele_eta_v_phi","h_pfo_ele_eta_v_phi",15,-4,4,15,-4,4)
       
#Resolution Hists
h_ele_etaRes=TH1D("h_ele_etaRes","h_ele_etaRes",50,-4,4)
h_ele_phiRes=TH1D("h_ele_phiRes","h_ele_phiRes",50,-4,4)
h_ele_thetaRes=TH1D("h_ele_thetaRes","h_ele_thetaRes",50,-4,4)
h_ele_pTRes=TH1D("h_ele_pTRes","h_ele_pTRes",50,-100,1200)



#Create a reader and open an LCIO file
#Find all files matching the directory pattern
histos_list = [h_truth_Rprod, h_truth_pT, h_truth_theta,
               h_trk_pT, h_trk_theta, h_trk_Rprod,
               h_ele_pT, h_ele_theta, h_ele_Rprod, h_all_p]

reader = IOIMPL.LCFactory.getInstance().createLCReader()
directory_pattern ='/collab/project/snowmass21/data/muonc/fmeloni/LegacyProductions/before29Jul23/DataMuC_MuColl_v1/electronGun/reco/electronGun_reco*'
#directory_pattern = '/collab/project/snowmass21/data/muonc/fmeloni/DataMuC_MuColl_v1/electronGun/reco/electronGun_reco_0.slcio' 
file_paths = glob.glob(directory_pattern)
reader.open(file_paths)

#Initalizing entries
pfo_array=[]
truth_array=[]
pdgid11=0
pdgid2112=0
Bfield = 3.56  # T
count=0
elecount=0
elecount2=0
count3=0

# looping over all events in the file

canvas=TCanvas("canvas","canvas",600,600)
for ievt, event in enumerate(reader):

    print("Processing event " + str(ievt))
    pfoCollection = event.getCollection('PandoraPFOs')
    trkCollection = event.getCollection('SiTracks_Refitted')    
#    relationCollection = event.getCollection('MCParticle_SiTracks_Refitted')
#    relation = UTIL.LCRelationNavigator(relationCollection)

    mcpCollection = event.getCollection('MCParticle')

#within file for only MC particles   
    for mcp in mcpCollection:
        charge = mcp.getCharge()
        status = mcp.getGeneratorStatus()
        count +=1
        if mcp.getPDG() not in truth_array:
            truth_array.append(mcp.getPDG())
        

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
                    #tracks = relation.getRelatedToObjects(mcp)

                    #if len(tracks) > 0:
                        #track = tracks[0]
                        #h_trk_Rprod.Fill(rprod)
                        #h_trk_pT.Fill(tlv.Perp())
                        #h_trk_theta.Fill(tlv.Theta())

                    for pfo in pfoCollection:
                        if fabs(pfo.getType()) == 11:
                            dp3 = pfo.getMomentum()
                            tlv_pfo = TLorentzVector()
                            tlv_pfo.SetPxPyPzE(
                                dp3[0], dp3[1], dp3[2], pfo.getEnergy())

                            pfotracks = pfo.getTracks()
                            if len(pfotracks) > 0:
                                trk = pfotracks[0]
        
                            #if track == trk:
                                #h_ele_Rprod.Fill(rprod)
                                #h_ele_pT.Fill(tlv.Perp())
                               #h_ele_theta.Fill(tlv.Theta()) 

#For ALL RECONSTRUCTED particles, not just MC/electrons
    for pfo in pfoCollection:
        h_all_p.Fill(pfo.getType())
        count3+=1
        if pfo.getType() not in pfo_array:      
            pfo_array.append(pfo.getType())

#electrons
        if fabs(pfo.getType()) == 11:
            pdgid11 = pdgid11+1
            dp3_ele=pfo.getMomentum() 
            tlv_pfo_ele = TLorentzVector()
            tlv_pfo_ele.SetPxPyPzE(dp3_ele[0], dp3_ele[1], dp3_ele[2], pfo.getEnergy())
            
            h_pfo_ele_pT.Fill(tlv_pfo_ele.Perp())
            h_pfo_ele_eta.Fill(tlv_pfo_ele.Eta())
            h_pfo_ele_phi.Fill(tlv_pfo_ele.Phi())
            h_pfo_ele_eta_v_phi.Fill(tlv_pfo_ele.Eta(),tlv_pfo_ele.Phi())
            h_pfo_ele_theta.Fill(tlv_pfo_ele.Theta())        

            h_ele_etaRes.Fill(tlv_pfo_ele.Eta() - tlv.Eta())
            h_ele_phiRes.Fill(tlv_pfo_ele.Phi() - tlv.Phi())
            h_ele_thetaRes.Fill(tlv_pfo_ele.Theta() - tlv.Theta())
            h_ele_pTRes.Fill(tlv_pfo_ele.Perp() - tlv.Perp())


#neutrons
        if fabs(pfo.getType()) == 2112:
            pdgid2112 = pdgid2112+1
            dp3_neut = pfo.getMomentum()
            tlv_pfo_neut = TLorentzVector()
            tlv_pfo_neut.SetPxPyPzE(dp3_neut[0], dp3_neut[1], dp3_neut[2], pfo.getEnergy())
           
            h_pfo_neut_pT.Fill(tlv.Perp())
            h_pfo_neut_eta.Fill(tlv.Eta())
            h_pfo_neut_phi.Fill(tlv.Phi())        
            h_pfo_neut_eta_v_phi.Fill(tlv_pfo_neut.Eta(),tlv_pfo_neut.Phi())

#Draw Histogram 

#Res Plots

#h_ele_etaRes.Draw()
#h_ele_phiRes.Draw()
#h_ele_thetaRes.Draw()
h_ele_pTRes.Draw()


h_all_p.GetXaxis().SetTitle("PFO Particle pdgID")
h_all_p.GetYaxis().SetTitle("Count")
h_all_p.GetXaxis().SetNdivisions(505)
#h_all_p.Draw()

  #Neutrons
h_pfo_neut_pT.SetTitle("Neutron pT")
h_pfo_neut_pT.GetXaxis().SetTitle("PFO Neutron pT [GeV]")
h_pfo_neut_pT.GetYaxis().SetTitle("")
h_pfo_neut_pT.GetXaxis().SetNdivisions(5)
#h_pfo_neut_pT.Draw()

h_pfo_neut_eta_v_phi.SetTitle("h_pfo_neut_eta_v_phi")
h_pfo_neut_eta_v_phi.GetYaxis().SetTitle("Phi")
h_pfo_neut_eta_v_phi.GetXaxis().SetTitle("Eta")
#h_pfo_neut_eta_v_phi.Draw()

#Electrons
h_pfo_ele_theta.SetTitle("PFO electron theta")
h_pfo_ele_theta.GetYaxis().SetTitle("")
h_pfo_ele_theta.GetXaxis().SetTitle("PFO Particle theta")
h_pfo_ele_theta.GetXaxis().SetNdivisions(505)
#h_pfo_ele_theta.Draw()

h_pfo_ele_eta.SetTitle("PFO electron eta")
h_pfo_ele_eta.GetYaxis().SetTitle("")
h_pfo_ele_eta.GetXaxis().SetTitle("PFO Electron Eta")
h_pfo_ele_eta.GetXaxis().SetNdivisions(10)
#h_pfo_ele_eta.Draw()

h_pfo_ele_pT.SetTitle("PFO electron pT")
h_pfo_ele_pT.GetYaxis().SetTitle("")
h_pfo_ele_pT.GetXaxis().SetTitle("PFO Electron pT [GeV]")
h_pfo_ele_eta.GetXaxis().SetNdivisions(10)
#h_pfo_ele_pT.Draw()


h_pfo_ele_phi.SetTitle("PFO electron phi")
h_pfo_ele_phi.GetYaxis().SetTitle("")
h_pfo_ele_phi.GetXaxis().SetTitle("PFO Electron Phi")
#h_pfo_ele_phi.Draw()

h_pfo_ele_eta_v_phi.SetTitle("h_pfo_ele_eta_v_phi")
h_pfo_ele_eta_v_phi.GetYaxis().SetTitle("PFO Electron Phi")
h_pfo_ele_eta_v_phi.GetXaxis().SetTitle("PFO Electron Eta")
#circle=TEllipse(0.0, 0.0, 0.3, 0.3) 
#h_pfo_ele_eta_v_phi.Draw("COLZ")


reader.close()
canvas.SaveAs("h_ele_pTRes.png")
# write histograms
#output_file = TFile(options.outDir + "ntup_elePFO_reco.root", 'RECREATE')
#for histo in histos_list:
#    histo.Write()
#output_file.Close()
