import ROOT
import numpy as np
import random
import math
ROOT.gROOT.SetBatch(1)

res = 0.05     #(%)
p_thresh = 0.05
thresh = 0.03 #(GeV)
proton_thresh = 0.05
mean = 1
sig_per_yr = 4500
filename = "/icarus/data/users/marshalc/ntuples/numi_wTarget.root"
beamdump_nu = ROOT.TVector3(0.39984193, 0.70939696, 0.58041570)

#kdar signal momentum from higher stats file
tf1 = ROOT.TFile("/icarus/data/users/marshalc/ntuples/numi_KDAR_signal.root")
tree1 = tf1.Get("tree")
for entry in tree1:
    p_nu = ROOT.TVector3(entry.nuPx, entry.nuPy, entry.nuPz)
    

#start of function
def TTree_loop(tree, hdict_sig, hdict_bkg):

    N = tree.GetEntries()
    n = 0                           

    n_signal = 0
    for entry in tree:
        
        if not n % 10000:
            print( "Event %d of %d..." % (n,N) )
        n += 1


        is_cc = False
        muon_E = None
        Emu_rec = None
        mu_p = None                
        true_nup = ROOT.TVector3(entry.nuPx, entry.nuPy, entry.nuPz)
        reco1_sum = ROOT.TVector3(0., 0., 0.)
        reco2_sum = ROOT.TVector3(0., 0., 0.)
        reco3_sum = ROOT.TVector3(0., 0., 0.)
        for i in range (entry.nfsp):
            if entry.pdg[i] > 9999: continue # skip over nuclear fragments, bindinos, etc.

            if entry.pdg[i] == 13:
                is_cc = True
                muon_E = entry.E[i]
                muon_ke = (muon_E - 0.105)
                mu_p = ROOT.TVector3( entry.px[i], entry.py[i], entry.pz[i] )
                x = random.gauss(mean, res)
                Emu_rec = (muon_ke) * x
                

            # Reco 1 = add every final state particle
            reco1_sum[0] += entry.px[i]
            reco1_sum[1] += entry.py[i]
            reco1_sum[2] += entry.pz[i]

            # Reco 2 = everything but neutrons
            if entry.pdg[i] != 2112:
                reco2_sum[0] += entry.px[i]
                reco2_sum[1] += entry.py[i]
                reco2_sum[2] += entry.pz[i]
    
                # Reco 3 = no neutrons, proton w/ threshold                                 
                if entry.pdg[i] != 2212 or (entry.E[i] - 0.93826 > proton_thresh):
                    reco3_sum[0] += (entry.px[i])
                    reco3_sum[1] += (entry.py[i])
                    reco3_sum[2] += (entry.pz[i])


        # Make a cut on CC events
        if is_cc:
              
            nu_p = ROOT.TVector3(entry.nuPx, entry.nuPy, entry.nuPz)
            theta = mu_p.Angle(nu_p)
            # True q and omega
            omega = (entry.Enu - muon_E)
            q = (nu_p - mu_p).Mag()
              
            # signal fills
            if entry.Enu >= 0.2355317 and entry.Enu <= 0.2355319:
                n_signal += 1
                hdict_sig["omega"].Fill(omega*1000.)
                hdict_sig["Emu"].Fill( (muon_E)*1000. ) # convert to MeV
                hdict_sig["omega and q"].Fill(omega*1000., q*1000)
                hdict_sig["omega and muon angle"].Fill(omega*1000., (theta))
                hdict_sig["omega vs q (truth panel)"].Fill((omega)*1000, (q)*1000)
                hdict_sig["true neutrino angles"].Fill(math.atan2(nu_p.x(), nu_p.z()), math.atan2(nu_p.y(), nu_p.z()))
                hdict_sig["reconstructed neutrino angle"].Fill(math.atan2(reco1_sum.x(),reco1_sum.z()), math.atan2(reco1_sum.y(),reco1_sum.z()))
                hdict_sig["reco.neutrino angle (no neutron)"].Fill(math.atan2(reco2_sum.x(),reco2_sum.z()), math.atan2(reco2_sum.y(),reco2_sum.z()))
                hdict_sig["reco. neutrino angle (signal proton threshold)"].Fill(math.atan2(reco3_sum.x(),reco3_sum.z()), math.atan2(reco3_sum.y(),reco3_sum.z()))
                hdict_sig["reco. assumption 1 signal wrst. beam dump"].Fill(math.cos(reco1_sum.Angle(beamdump_nu)))
                hdict_sig["reco. assumption 2 wrst. beam dump"].Fill(math.cos(reco2_sum.Angle(beamdump_nu)))
                hdict_sig["reco. assumption 3 wrst. beam dump"].Fill(math.cos(reco3_sum.Angle(beamdump_nu)))
                hdict_sig["muon energy vs. beam dump angle (sig 1) "].Fill(muon_E*1000, math.cos(reco1_sum.Angle(beamdump_nu)))
                hdict_sig["muon energy vs. beam dump angle (sig 2) "].Fill(muon_E*1000, math.cos(reco2_sum.Angle(beamdump_nu)))
                hdict_sig["muon energy vs. beam dump angle (sig 3) "].Fill(muon_E*1000, math.cos(reco3_sum.Angle(beamdump_nu)))
                hdict_sig["omega q stack (truth)"].Fill(omega * 1000, q * 1000)
                if muon_E < 0.236:
                    hdict_sig["muon energy (with threshold) vs. beamdump angle (sig 1)"].Fill(muon_E*1000, math.cos(reco1_sum.Angle(beamdump_nu)))
                    hdict_sig["muon energy (with threshold) vs. beamdump angle (sig 2)"].Fill(muon_E*1000, math.cos(reco2_sum.Angle(beamdump_nu)))
                    hdict_sig["muon energy (with threshold) vs. beamdump angle (sig 3)"].Fill(muon_E*1000, math.cos(reco3_sum.Angle(beamdump_nu)))
                if muon_ke > thresh:
                    hdict_sig["omega(rec) and q"].Fill((entry.Enu - (Emu_rec + 0.105))*1000,  q*1000)
                    hdict_sig["omega vs. q (reco. panel)"].Fill((entry.Enu - (Emu_rec + 0.105))*1000,  q*1000)
                    hdict_sig["fractional residual"].Fill((Emu_rec - (muon_ke))/(muon_ke), (muon_ke))
                    hdict_sig["rec-muon energy and true muon E"].Fill((muon_E - 0.105)*1000, (Emu_rec)*1000)
                    hdict_sig["muon energy and angle"].Fill((Emu_rec)*1000, (theta))
                    hdict_sig["omega q stack (reco)"].Fill((entry.Enu - (Emu_rec + 0.105))*1000,  q*1000)
            else:
                #background fills
                hdict_bkg["neutrino angle wrst. Z"].Fill(math.atan2(nu_p.x(), nu_p.z()), math.atan2(nu_p.y(), nu_p.z()))
                hdict_bkg["reco neutrino angle"].Fill(math.atan2(reco1_sum.x(),reco1_sum.z()), math.atan2(reco1_sum.y(),reco1_sum.z()))
                hdict_bkg["reco. neutrino angle (no neutron)"].Fill(math.atan2(reco2_sum.x(),reco2_sum.z()), math.atan2(reco2_sum.y(),reco2_sum.z()))
                hdict_bkg["reco. neutrino angle (proton threshold)"].Fill(math.atan2(reco3_sum.x(),reco3_sum.z()), math.atan2(reco3_sum.y(),reco3_sum.z()))
                hdict_bkg["beam dump angle (reco. 1)"].Fill(math.cos(reco1_sum.Angle(beamdump_nu)))
                hdict_bkg["beam dump angle (reco. 2)"].Fill(math.cos(reco2_sum.Angle(beamdump_nu)))
                hdict_bkg["beam dump angle (reco. 3)"].Fill(math.cos(reco3_sum.Angle(beamdump_nu)))    
                if muon_E < 0.236:
                    hdict_bkg["beam dump angle (muon thresh reco. 1)"].Fill(math.cos(reco1_sum.Angle(beamdump_nu)))
                    hdict_bkg["beam dump angle (muon thresh reco. 2)"].Fill(math.cos(reco2_sum.Angle(beamdump_nu)))
                    hdict_bkg["beam dump angle (muon thresh reco. 3)"].Fill(math.cos(reco3_sum.Angle(beamdump_nu)))
                    hdict_bkg["muon energy (with threshold) vs. beamdump angle (1)"].Fill(muon_E*1000, math.cos(reco1_sum.Angle(beamdump_nu)))
                    hdict_bkg["muon energy (with threshold) vs. beamdump angle (2)"].Fill(muon_E*1000, math.cos(reco2_sum.Angle(beamdump_nu)))
                    hdict_bkg["muon energy (with threshold) vs. beamdump angle (3)"].Fill(muon_E*1000, math.cos(reco3_sum.Angle(beamdump_nu)))
                hdict_bkg["muon energy vs. beam dump angle (1) "].Fill(muon_E*1000, math.cos(reco1_sum.Angle(beamdump_nu)))
                hdict_bkg["muon energy vs. beam dump angle (2) "].Fill(muon_E*1000, math.cos(reco2_sum.Angle(beamdump_nu)))
                hdict_bkg["muon energy vs. beam dump angle (3)"].Fill(muon_E*1000, math.cos(reco3_sum.Angle(beamdump_nu)))
                if math.cos(p_nu.Angle(beamdump_nu)) > .99:         
                        hdict_bkg["omega vs. q (bkg)"].Fill((0.236 - muon_E) * 1000, ((p_nu - mu_p).Mag())*1000)
                        hdict_bkg["omega vs. q (stack)"].Fill((0.236 - muon_E) * 1000, ((p_nu - mu_p).Mag())*1000)
                        hdict_bkg["omega v. q (truth panel)"].Fill((0.236 - muon_E) * 1000, ((p_nu - mu_p).Mag())*1000)
                        if muon_ke > thresh:
                            hdict_bkg["omega and q (reco panel)"].Fill((0.236 - (Emu_rec + 0.105))*1000, ((p_nu - mu_p).Mag())*1000)
                            hdict_bkg["omega q (reco.)"].Fill((0.236 - (Emu_rec + 0.105))*1000, ((p_nu - mu_p).Mag())*1000)
                            hdict_bkg["omega q (reco stack)"].Fill((0.236 - (Emu_rec + 0.105))*1000, ((p_nu - mu_p).Mag())*1000)
    return n_signal

if __name__ == "__main__":
               
    #background plots (cosmetic changes here)
    hdict_bkg = {}
    #bkg plot (1) (projection of true angle b/w nu and ICARUS in XZ plane vs. YZ plane)
    hdict_bkg["neutrino angle wrst. Z"] = ROOT.TH2D("neutrino angle wrst. Z", ";#theta_{#nu (XZ)} ;#theta_{#nu (YZ)}", 136, 0., 1., 136, 0., 1.)
    #bkg plot (2)  (direction reco. assumption 1: angle b/w nu (bkg) and ICARUS in XZ vs. YZ plane)
    hdict_bkg["reco neutrino angle"] =  ROOT.TH2D("reco neutrino angle", ";#theta_{#nu (XZ)} ;#theta_{#nu (YZ)}", 136, -1.6, 1.6, 136, -1.6, 1.6)
    #bkg plot (3) (direction reco. assumption 2: angle b/w nu (bkg) and ICARUS in XZ vs. YZ plane)
    hdict_bkg["reco. neutrino angle (no neutron)"] = ROOT.TH2D("reco. neutrino angle (no neutron)", ";#theta_{#nu (XZ)} ;#theta_{#nu (YZ)}", 136, -1.6, 1.6, 136, -1.6, 1.6)
    #bkg plot (4) (direction reco. assumption 3: angle b/w nu (bkg) and ICARUS in XZ vs. YZ plane)
    hdict_bkg["reco. neutrino angle (proton threshold)"] = ROOT.TH2D("reco. neutrino angle (proton threshold)", ";#theta_{#nu (XZ)} ;#theta_{#nu (YZ)}", 136, -1.6, 1.6, 136, -1.6, 1.6)
    #bkg plot (5); for next 3 plots: cos of angle between neutrinos for each direction reco. assumption and NuMI beam dump
    hdict_bkg["beam dump angle (reco. 1)"] = ROOT.TH1D("beam dump angle (reco. 1)", ";#theta_{#nu}", 136, -1., 1.)
    hdict_bkg["beam dump angle (reco. 2)"] = ROOT.TH1D("beam dump angle (reco. 2)", ";#theta_{#nu}", 136, -1., 1.)
    hdict_bkg["beam dump angle (reco. 3)"] = ROOT.TH1D("beam dump angle (reco. 3)", ";#theta_{#nu}", 136, -1., 1.)
    #bkg plot (6); for next 3 plots: cos of angle b/w neutrinos for each direction reco. assumption and beam dump w/236 MeV muon thresh
    hdict_bkg["beam dump angle (muon thresh reco. 1)"] = ROOT.TH1D("beam dump angle (muon thresh reco. 1)", ";#theta_{#nu}", 136, -1., 1.)
    hdict_bkg["beam dump angle (muon thresh reco. 2)"] = ROOT.TH1D("beam dump angle (muon thresh reco. 2)", ";#theta_{#nu}", 136, -1., 1.)
    hdict_bkg["beam dump angle (muon thresh reco. 3)"] = ROOT.TH1D("beam dump angle (muon thresh reco. 3)", ";#theta_{#nu}", 136, -1., 1.)
    #bkg plot (7); muon E vs. cos of angle b/w nu (direction reco. 1) and beam dump
    hdict_bkg["muon energy vs. beam dump angle (1) "] = ROOT.TH2D("muon energy vs. beam dump angle (1) ", ";E_{#mu} MeV ;cos(#theta_{#nu})", 136, 100., 250., 136, -1., 1.)
    #bkg plot (8); muon E vs. cos of angle b/w nu (direction reco. 2) and beam dump
    hdict_bkg["muon energy vs. beam dump angle (2) "] = ROOT.TH2D("muon energy vs. beam dump angle (2) ", ";E_{#mu} MeV ;cos(#theta_{#nu})", 136, 100., 250., 136, -1., 1.)
    #bkg plot (9); muon E vs. cos of angle b/w nu (direction reco. 3) and beam dump
    hdict_bkg["muon energy vs. beam dump angle (3)"] = ROOT.TH2D("muon energy vs. beam dump angle (3)", ";E_{#mu} MeV ;cos(#theta_{#nu})", 136, 100., 250., 136, -1., 1.)
    #bkg plot (10); muon E (w/236 MeV Threshold) vs. cos of angle b/w nu (direction reco. 1) and beam dump
    hdict_bkg["muon energy (with threshold) vs. beamdump angle (1)"] = ROOT.TH2D("muon energy (with threshold) vs. beamdump angle (1)", ";E_{#mu} MeV ;cos(#theta_{#nu})", 136, 100, 240., 136, -1., 1.)
    #bkg plot (11); muon E (w/236 MeV Threshold) vs. cos of angle b/w nu (direction reco. 2) and beam dump
    hdict_bkg["muon energy (with threshold) vs. beamdump angle (2)"] = ROOT.TH2D("muon energy (with threshold) vs. beamdump angle (2)", ";E_{#mu} MeV ;cos(#theta_{#nu})", 136, 100, 240., 136, -1., 1.)
    #bkg plot (12); muon E (w/236 MeV Threshold) vs. cos of angle b/w nu (direction reco. 3) and beam dump
    hdict_bkg["muon energy (with threshold) vs. beamdump angle (3)"] = ROOT.TH2D("muon energy (with threshold) vs. beamdump angle (3)", ";E_{#mu} MeV ;cos(#theta_{#nu})", 136, 100, 240., 136, -1., 1.)
    #bkg plot (13); omega vs. q (2D)
    hdict_bkg["omega vs. q (bkg)"] = ROOT.TH2D("omega vs. q (bkg)",";#omega = E_{#nu} - E_{#mu} (MeV), ;q=|p_{#nu}-p_{#mu}| (MeV/c)", 136, 0., 140., 136, 0., 450. )
    #bkg plot (14); reco. omega vs. q (2D)
    hdict_bkg["omega q (reco.)"] = ROOT.TH2D( "omega and q (reco.)", ";#omega = E_{#nu} - E_{#mu} (MeV), ;q=|p_{#nu}-p_{#mu}| (MeV/c)", 136, 0., 140., 136, 0., 450. )
    #bkg plot (15); omega vs. q paneled plot (truth level)
    hdict_bkg["omega v. q (truth panel)"] = ROOT.TH2D("omega vs. q (stack)", ";#omega = E_{#nu}-E_{#mu} (MeV) ;q =|p_{#nu}-p_{#mu}| (MeV/c)", 80, 0., 140., 10, 0., 500.)
    #bkg plot (16); omega vs. q paneled plot with reco. muon E
    hdict_bkg["omega and q (reco panel)"]= ROOT.TH2D("omega and q (reco panel)", ";#omega = E_{#nu} - E_{#mu (rec)} (MeV) ;q=|p_{#nu}-p_{#mu}| (MeV/c)", 80, 0., 140., 10, 0., 500.)
    #stacked plot 1; omega vs. q (truth level for stacked plot w/signal)
    hdict_bkg["omega vs. q (stack)"] = ROOT.TH2D("omega vs. q (stack)", ";#omega = E_{#nu}-E_{#mu} (MeV) ;q =|p_{#nu}-p_{#mu}| (MeV/c)", 80, 0., 140., 10, 0., 500.)
    #stacked plot 2; omega vs. q paneled plot with reco. muon E for stacked plot w/reco. signal
    hdict_bkg["omega q (reco stack)"] = ROOT.TH2D("omega q (reco stack)", ";#omega = E_{#nu} - E_{#mu (rec)} (MeV) ;q=|p_{#nu}-p_{#mu}| (MeV/c)", 80, 0., 140., 10, 0., 500.)
    
    #signal plots (cosmetic changes here:)
    hdict_sig = {}
    #sig plot (1); omega
    hdict_sig["omega"] = ROOT.TH1D("omega",";#omega = E_{#nu} - E_{#mu} (MeV)", 136, 0., 136. )
    #sig plot (2); muon energy
    hdict_sig["Emu"] = ROOT.TH1D( "E_{#mu}", ";E_{#mu} (MeV)", 136, 0., 236. )
    #sig plot (3); omega vs. q (truth)
    hdict_sig["omega and q"] = ROOT.TH2D( "omega and q", ";#omega = E_{#nu} - E_{#mu} (MeV), ;q=|p_{#nu}-p_{#mu}| (MeV/c)", 136, 0., 140., 136, 0., 450. )
    #sig plot (4); omega vs. muon deflection angle
    hdict_sig["omega and muon angle"] = ROOT.TH2D( "omega and muon angle", ";#omega = E_{#nu} - E_{#mu} (MeV) ;muon angle", 136, 0., 140., 136, 0., 5. )
    #sig plot (5); muon energy vs. muon deflection angle
    hdict_sig["muon energy and angle"] = ROOT.TH2D( "Muon energy and angle ", ";E_{#mu} (MeV) ;muon angle (rad)", 136, 0., 300., 136, 0., 5.)
    #sig plot (6); reco. muon energy vs. muon deflection angle
    hdict_sig["reconstructed muon energy and muon angle"] = ROOT.TH2D( "reconstructed muon energy and muon angle", ";E_{#mu (rec)} = x * E_{#mu} ;muon angle (rad)", 136, 0 , 300., 136, 0. , 5.)
    #sig plot (7); reco. muon energy vs. true muon energy
    hdict_sig["rec-muon energy and true muon E"] = ROOT.TH2D( "rec-muon energy and true muon E", "; E_{#mu} (MeV) ;E_{#mu (rec)} = x * E_{#mu} (MeV) ", 136, 0. , 350., 136, 0. , 350.)
    #sig plot (8); fractional residual of reco. muon energy and muon energy
    hdict_sig["fractional residual"] = ROOT.TH2D("frational residual", ";(E_{#mu (rec)}-E_{#mu})/E_{#mu}, ;E_{#mu}(GeV) ", 136, -1., 1., 136, 0., 0.5)
    #sig plot (9); reco. omega vs. q (5%fractional res. and 30 MeV thresh)
    hdict_sig["omega(rec) and q"] = ROOT.TH2D("omega (rec) and q", ";#omega = E_{#nu} - E_{#mu (rec)} (MeV) ;q=|p_{#nu}-p_{#mu}| (MeV/c)", 136, 0., 140., 136, 0., 450.)
    #sig plot (10); omega vs. q paneled plot
    hdict_sig["omega vs q (truth panel)"] = ROOT.TH2D("omega vs q (truth panel)", ";#omega = E_{#nu}-E_{#mu} (MeV) ;q =|p_{#nu}-p_{#mu}| (MeV/c)", 136, 0., 140., 10, 0., 500.)
    #sig plot (11); reco. omega vs. q paneled plot
    hdict_sig["omega vs. q (reco. panel)"] = ROOT.TH2D("omega vs. q (reco. panel)", ";#omega = E_{#nu}-E_{#mu} (MeV) ;q =|p_{#nu}-p_{#mu}| (MeV/c)", 50, 0., 140., 10, 0., 500.)
    #sig plot (12); true neutrino direction wrst. ICARUS in XZ vs. YZ plane
    hdict_sig["true neutrino angles"] = ROOT.TH2D("true neutrino angles", ";#theta_{#nu (XZ)} ;#theta_{#nu (YZ)}", 136, 0., 1., 136, 0., 1.)
    #sig plot (13); angle of nu direction reco. 1 in XZ vs. YZ plane wrst. ICARUS
    hdict_sig["reconstructed neutrino angle"] = ROOT.TH2D("reconstructed neutrino angle", ";#theta_{#nu (XZ)} ;#theta_{#nu (YZ)}", 136, -1.6, 1.6, 136, -1.6, 1.6)
    #sig plot (14); angle of reco. nu direction 2 in XZ vs. YZ plane wrst. ICARUS
    hdict_sig["reco.neutrino angle (no neutron)"] = ROOT.TH2D("reco.neutrino angle (no neutron)", ";#theta_{#nu (XZ)} ;#theta_{#nu (YZ)}", 136, -1.6, 1.6, 136, -1.6, 1.6)
    #sig plot (15); angle of reco. nu direction 3 in XZ vs. YZ plane wrst. ICARUS
    hdict_sig["reco. neutrino angle (signal proton threshold)"] = ROOT.TH2D("reco. neutrino angle (signal proton threshold)", ";#theta_{#nu (XZ)} ;#theta_{#nu (YZ)}", 136, -1.6, 1.6, 136, -1.6, 1.6)
    #sig plot (16); for next 3 plots: cos of angle between neutrinos for each direction reco. assumption and NuMI beam dump
    hdict_sig["reco. assumption 1 signal wrst. beam dump"] = ROOT.TH1D("reco. assumption 1 signal wrst. beam dump" , ";#theta_{#nu}", 136, -1., 1.)
    hdict_sig["reco. assumption 2 wrst. beam dump"] = ROOT.TH1D("reco. assumption 2 wrst. beam dump", ";#theta_{#nu}", 136, -1., 1.)
    hdict_sig["reco. assumption 3 wrst. beam dump"] = ROOT.TH1D("reco. assumption 3 wrst. beam dump", ";#theta_{#nu}", 136, -1., 1.)
    #sig plot (17); muon energy vs. cos of angle b/w reco. (1) nu direction and beam dump
    hdict_sig["muon energy vs. beam dump angle (sig 1) "] = ROOT.TH2D("muon energy vs. beam dump angle (sig 1) ", ";E_{#mu} MeV ;cos(#theta_{#nu})", 136, 80., 250., 136, -1., 1.)
    #sig plot (18); muon energy vs. cos of angle b/w reco. (2) nu direction and beam dump
    hdict_sig["muon energy vs. beam dump angle (sig 2) "] = ROOT.TH2D("muon energy vs. beam dump angle (sig 2) ", ";E_{#mu} MeV ;cos(#theta_{#nu})", 136, 80., 250., 136, -1., 1.)
    #sig plot (19); muon energy vs. cos of angle b/w reco. (3) nu direction and beam dump
    hdict_sig["muon energy vs. beam dump angle (sig 3) "] = ROOT.TH2D("muon energy vs. beam dump angle (sig 3) ", ";E_{#mu} MeV ;cos(#theta_{#nu})", 136, 80., 250., 136, -1., 1.)
    #sig plot (20); muon energy vs. cos of angle b/w reco. (1) nu direction and beam dump (with 236 MeV muon thresh)
    hdict_sig["muon energy (with threshold) vs. beamdump angle (sig 1)"] = ROOT.TH2D("muon energy (with threshold) vs. beamdump angle (sig 1)", ";E_{#mu} MeV ;cos(#theta_{#nu})", 136, 100, 240., 136, -1., 1.)
    #sig plot (21); muon energy vs. cos of angle b/w reco. (2) nu direction and beam dump (with 236 MeV muon thresh)
    hdict_sig["muon energy (with threshold) vs. beamdump angle (sig 2)"] = ROOT.TH2D("muon energy (with threshold) vs. beamdump angle (sig 2)", ";E_{#mu} MeV ;cos(#theta_{#nu})", 136, 100, 240., 136, -1., 1.)
    #sig plot (22); muon energy vs. cos of angle b/w reco. (3) nu direction and beam dump (with 236 MeV muon thresh)
    hdict_sig["muon energy (with threshold) vs. beamdump angle (sig 3)"] = ROOT.TH2D("muon energy (with threshold) vs. beamdump angle (sig 3)", ";E_{#mu} MeV ;cos(#theta_{#nu})", 136, 100, 240., 136, -1., 1.)
    #stack plot 1;
    hdict_sig["omega q stack (truth)"] = ROOT.TH2D("omega q stack (truth)", ";#omega = E_{#nu}-E_{#mu} (MeV) ;q =|p_{#nu}-p_{#mu}| (MeV/c)", 80, 0., 140., 10, 0., 500.)
    #stack plot 2;
    hdict_sig["omega q stack (reco)"] = ROOT.TH2D("omega q stack (reco)", ";#omega = E_{#nu}-E_{#mu} (MeV) ;q =|p_{#nu}-p_{#mu}| (MeV/c)", 80, 0., 140., 10, 0., 500.)

    # run it
    tf = ROOT.TFile( filename )
    tree = tf.Get( "tree" )
    n_signal = TTree_loop( tree, hdict_sig, hdict_bkg )

    norm = 1.*sig_per_yr / n_signal
          
    # Make an output file for the histograms
    tfout = ROOT.TFile( "plots.root", "RECREATE" ) # overwrite the file if it already exists
    hdict_sig["omega"].Write()
    
    #signal drawings
    #sig plot (1)
    c = ROOT.TCanvas()
    hdict_sig["omega"].Scale(norm)
    for b in range (1, hdict_sig["omega"].GetNbinsX()+1):
        content = hdict_sig["omega"].GetBinContent(b)
        error = content**0.5
        hdict_sig["omega"].SetBinError(b, error)
    hdict_sig["omega"].SetLineWidth(1)
    hdict_sig["omega"].SetMarkerSize(0.5)
    hdict_sig["omega"].Draw("hist E")
    c.Print( "omega.png" )
    
    #sig plot (2)
    c1 = ROOT.TCanvas()
    hdict_sig["Emu"].Scale(norm)
    for b in range (1, hdict_sig["Emu"].GetNbinsX()+1):
        content = hdict_sig["Emu"].GetBinContent(b)
        error = content**0.5
        hdict_sig["Emu"].SetBinError(b, error)
    hdict_sig["Emu"].SetLineWidth(1)
    hdict_sig["Emu"].SetMarkerSize(0.5)
    hdict_sig["Emu"].Draw("hist E")
    c1.Print("emu.png")
    
    #sig plot (3)
    c2 = ROOT.TCanvas()
    hdict_sig["omega and q"].Scale(norm)
    hdict_sig["omega and q"].Draw("colz")
    c2.Print("omega-q.png")
    
    #sig plot (4)
    c3 = ROOT.TCanvas()
    hdict_sig["omega and muon angle"].Scale(norm)
    hdict_sig["omega and muon angle"].Draw("colz")
    c3.Print("omega-angle.png")
    
    #sig plot (5)
    c4 = ROOT.TCanvas()
    hdict_sig["muon energy and angle"].Scale(norm)
    hdict_sig["muon energy and angle"].Draw("colz")
    c4.Print("emu-angle.png")
    
    #sig plot (6)
    c5 = ROOT.TCanvas()
    hdict_sig["reconstructed muon energy and muon angle"].Scale(norm)
    hdict_sig["reconstructed muon energy and muon angle"].Draw("colz")
    c5.Print("emu.rec-angle.png")
    
    #sig plot (7)
    c6 = ROOT.TCanvas()
    hdict_sig["rec-muon energy and true muon E"].Scale(norm)
    hdict_sig["rec-muon energy and true muon E"].Draw("colz")
    c6.Print("emur-emu.png")
    
    #sig plot (8)
    c7 = ROOT.TCanvas()
    hdict_sig["fractional residual"].Scale(norm)
    hdict_sig["fractional residual"].Draw("colz")
    c7.Print("erec.png")
    
    #sig plot (9)
    c8 = ROOT.TCanvas()
    hdict_sig["omega(rec) and q"].Scale(norm)
    hdict_sig["omega(rec) and q"].Draw("colz")
    c8.Print("omegaq.rec.png")
    
    #sig plot (10)
    c11 = ROOT.TCanvas("omega and q", "multipads", 1000, 500)
    c11.Divide(5, 2, 0, 0)
    hdict_sig["omega vs q (truth panel)"].Scale(norm)
    maximum = hdict_sig["omega vs q (truth panel)"].GetMaximum()
    for i in range (0, 10):
        c11.cd(i+1)
        hq = hdict_sig["omega vs q (truth panel)"].ProjectionX('name_%d' % i, i, i)
        hq.GetYaxis().SetRangeUser(0., maximum*1.2)
        for b in range (1, hq.GetNbinsX()+1):
            content = hq.GetBinContent(b)
            error = ((content)**0.5)
            hq.SetBinError(b, error)
        hq.SetMarkerSize(0.5)
        hq.SetLineWidth(1)
        hq.Draw("hist E")
    c11.Print("omegaq-multiloop.png")
    
    #sig plot (11)
    c12 = ROOT.TCanvas("omega vs. q (rec)", "multipads", 1000, 500)
    c12.Divide(5,2,0,0)
    hdict_sig["omega vs. q (reco. panel)"].Scale(norm)
    maximum = hdict_sig["omega vs. q (reco. panel)"].GetMaximum()
    for i in range (0, 10):
        c12.cd(i+1)
        hqrec = hdict_sig["omega vs. q (reco. panel)"].ProjectionX('name_%d' % i, i, i)
        hqrec.GetYaxis().SetRangeUser(0., maximum*1.2)
        for b in range (1, hqrec.GetNbinsX()+1):
            content = hqrec.GetBinContent(b)
            error = (content)**0.5
            hqrec.SetBinError(b, error)
        hqrec.SetLineWidth(1)
        hqrec.SetMarkerSize(0.5)
        hqrec.Draw("hist E")  
    c12.Print("omegaq-rec-multi.png")
     
    #sig plot (12)
    c13 = ROOT.TCanvas()      
    hdict_sig["true neutrino angles"].Scale(norm)
    hdict_sig["true neutrino angles"].Draw("colz")
    c13.Print("true_signal_angles.png")
    
    #sig plot (13)
    c14 = ROOT.TCanvas()     
    hdict_sig["reconstructed neutrino angle"] .Scale(norm)
    hdict_sig["reconstructed neutrino angle"] .Draw("colz")
    c14.Print("recoFS_signal_angle.png")
    
    #sig plot (14)
    c15 = ROOT.TCanvas()     
    hdict_sig["reco.neutrino angle (no neutron)"].Scale(norm)
    hdict_sig["reco.neutrino angle (no neutron)"].Draw("colz")
    c15.Print("recoFS_signal_neu.png")
    
    #sig plot (15)
    c16 = ROOT.TCanvas()       
    hdict_sig["reco. neutrino angle (signal proton threshold)"].Scale(norm)
    hdict_sig["reco. neutrino angle (signal proton threshold)"].Draw("colz")
    c16.Print("signal_proton_thresh.png")
       
    #sig plot (16)    
    c17 = ROOT.TCanvas()
    leg = ROOT.TLegend()
    leg.AddEntry(hdict_sig["reco. assumption 1 signal wrst. beam dump"],"Sum of all FS particles","f")
    hdict_sig["reco. assumption 1 signal wrst. beam dump"].Scale(norm)
    hdict_sig["reco. assumption 1 signal wrst. beam dump"].SetLineColor(ROOT.kRed)
    hdict_sig["reco. assumption 1 signal wrst. beam dump"].Draw("hist")
           
    hdict_sig["reco. assumption 2 wrst. beam dump"] .Scale(norm)
    leg.AddEntry(hdict_sig["reco. assumption 2 wrst. beam dump"] ,"Sum of all FS particles excluding neutrons","f")
    hdict_sig["reco. assumption 2 wrst. beam dump"] .SetLineColor(ROOT.kBlue)
    hdict_sig["reco. assumption 2 wrst. beam dump"] .Draw("hist same")
           
    hdict_sig["reco. assumption 3 wrst. beam dump"].Scale(norm)
    leg.AddEntry(hdict_sig["reco. assumption 3 wrst. beam dump"],"Sum of all FS particles with 50 MeV proton threshold","f")
    hdict_sig["reco. assumption 3 wrst. beam dump"].SetLineColor(ROOT.kGreen+2)
    hdict_sig["reco. assumption 3 wrst. beam dump"].Draw("hist same")
    leg.Draw()
    c17.Print("signal_beamdump_angles.png")
    
    #sig plot (17)
    c18 = ROOT.TCanvas()      
    hdict_sig["muon energy vs. beam dump angle (sig 1) "].Scale(norm)
    hdict_sig["muon energy vs. beam dump angle (sig 1) "].Draw("colz")
    c18.Print("sig_muonE_angle1.png")
    
    #sig plot (18)
    c19 = ROOT.TCanvas()      
    hdict_sig["muon energy vs. beam dump angle (sig 2) "].Scale(norm)
    hdict_sig["muon energy vs. beam dump angle (sig 2) "].Draw("colz")
    c19.Print("sig_muonE_angle2.png")
    
    #sig plot (19)
    c20 = ROOT.TCanvas()    
    hdict_sig["muon energy vs. beam dump angle (sig 3) "].Scale(norm)
    hdict_sig["muon energy vs. beam dump angle (sig 3) "].Draw("colz")
    c20.Print("sig_muonE_angle3.png")
    
    #sig plot (20)
    c21 = ROOT.TCanvas()      
    hdict_sig["muon energy (with threshold) vs. beamdump angle (sig 1)"].Scale(1/ hdict_sig["muon energy (with threshold) vs. beamdump angle (sig 1)"].Integral())
    hdict_sig["muon energy (with threshold) vs. beamdump angle (sig 1)"].Draw("colz")
    c21.Print("sig_muonE_thresh_angle1.png")
    
    #sig plot (21)
    c22 = ROOT.TCanvas()       
    hdict_sig["muon energy (with threshold) vs. beamdump angle (sig 2)"].Scale(1/hdict_sig["muon energy (with threshold) vs. beamdump angle (sig 2)"].Integral())
    hdict_sig["muon energy (with threshold) vs. beamdump angle (sig 2)"].Draw("colz")
    c22.Print("sig_muonE_thresh_angle2.png")
    
    #sig plot (22)
    c23 = ROOT.TCanvas()      
    hdict_sig["muon energy (with threshold) vs. beamdump angle (sig 3)"].Scale(1/hdict_sig["muon energy (with threshold) vs. beamdump angle (sig 3)"].Integral())
    hdict_sig["muon energy (with threshold) vs. beamdump angle (sig 3)"].Draw("colz")
    c23.Print("sig_muonE_thresh_angle3.png")
    
    #background drawings  
    #bkg plot (1)
    u1 = ROOT.TCanvas()
    hdict_bkg["neutrino angle wrst. Z"].Scale(norm)
    hdict_bkg["neutrino angle wrst. Z"].Draw("colz")
    u1.Print("true_background_angle.png")
    
    #bkg plot (2)
    u2 = ROOT.TCanvas()
    hdict_bkg["reco neutrino angle"].Scale(norm)
    hdict_bkg["reco neutrino angle"].Draw("colz")
    u2.Print("recoFS_background_angle.png")
    
    #bkg plot (3)
    u3 = ROOT.TCanvas()           
    hdict_bkg["reco. neutrino angle (no neutron)"].Scale(norm)
    hdict_bkg["reco. neutrino angle (no neutron)"].Draw("colz")
    u3.Print("reco_background_neu.png")
     
    #bkg plot (4)
    u4 = ROOT.TCanvas()         
    hdict_bkg["reco. neutrino angle (proton threshold)"].Scale(norm)
    hdict_bkg["reco. neutrino angle (proton threshold)"].Draw("colz")
    u4.Print("background_proton_thresh.png")
    
    #bkg plot (5)       
    u5 = ROOT.TCanvas()
    legend = ROOT.TLegend()
    legend.AddEntry(hdict_bkg["beam dump angle (reco. 1)"],"Sum of all FS particles","f")
    hdict_bkg["beam dump angle (reco. 1)"].Scale(norm)
    hdict_bkg["beam dump angle (reco. 1)"].SetLineColor(ROOT.kRed)
    hdict_bkg["beam dump angle (reco. 1)"].Draw("hist")
           
    hdict_bkg["beam dump angle (reco. 2)"].Scale(norm)
    legend.AddEntry(hdict_bkg["beam dump angle (reco. 2)"],"Sum of all FS particles excluding neutrons","f")
    hdict_bkg["beam dump angle (reco. 2)"].SetLineColor(ROOT.kBlue)
    hdict_bkg["beam dump angle (reco. 2)"].Draw("hist same")
           
    hdict_bkg["beam dump angle (reco. 3)"].Scale(norm)
    legend.AddEntry(hdict_bkg["beam dump angle (reco. 3)"],"Sum of all FS particles with 50 MeV proton threshold","f")
    hdict_bkg["beam dump angle (reco. 3)"].SetLineColor(ROOT.kGreen+2)
    hdict_bkg["beam dump angle (reco. 3)"].Draw("hist same")
    legend.Draw()
    u5.Print("background_beamdump_angles.png")
     
    #bkg plot (6)         
    u6 = ROOT.TCanvas()
    legend = ROOT.TLegend()
    legend.AddEntry(hdict_bkg["beam dump angle (muon thresh reco. 1)"],"Sum of all FS particles","f")
    hdict_bkg["beam dump angle (muon thresh reco. 1)"].Scale(norm)
    hdict_bkg["beam dump angle (muon thresh reco. 1)"].SetLineColor(ROOT.kRed)
    hdict_bkg["beam dump angle (muon thresh reco. 1)"].Draw("hist")
           
    hdict_bkg["beam dump angle (muon thresh reco. 2)"].Scale(norm)
    legend.AddEntry(hdict_bkg["beam dump angle (muon thresh reco. 2)"],"Sum of all FS particles excluding neutrons","f")
    hdict_bkg["beam dump angle (muon thresh reco. 2)"].SetLineColor(ROOT.kBlue)
    hdict_bkg["beam dump angle (muon thresh reco. 2)"].Draw("hist same")
           
    hdict_bkg["beam dump angle (muon thresh reco. 3)"].Scale(norm)
    legend.AddEntry(hdict_bkg["beam dump angle (muon thresh reco. 3)"],"Sum of all FS particles with 50 MeV proton threshold","f")
    hdict_bkg["beam dump angle (muon thresh reco. 3)"].SetLineColor(ROOT.kGreen+2)
    hdict_bkg["beam dump angle (muon thresh reco. 3)"].Draw("hist same")
    legend.Draw()
    u6.Print("background_beamdump_thresh.png")
    
    #bkg plot (7)
    u7 = ROOT.TCanvas()         
    hdict_bkg["muon energy vs. beam dump angle (1) "].Scale(norm)
    hdict_bkg["muon energy vs. beam dump angle (1) "].Draw("colz")
    u7.Print("bg_muonE_angle1.png")
    
    #bkg plot (8)
    u8 = ROOT.TCanvas()       
    hdict_bkg["muon energy vs. beam dump angle (2) "].Scale(norm)
    hdict_bkg["muon energy vs. beam dump angle (2) "].Draw("colz")
    u8.Print("bg_muonE_angle2.png")
    
    #bkg plot (9)
    u9 = ROOT.TCanvas()         
    hdict_bkg["muon energy vs. beam dump angle (3)"].Scale(norm)
    hdict_bkg["muon energy vs. beam dump angle (3)"].Draw("colz")
    u9.Print("bg_muonE_angle3.png")
    
    #bkg plot (10)
    u10 = ROOT.TCanvas()         
    hdict_bkg["muon energy (with threshold) vs. beamdump angle (1)"].Scale(1/hdict_bkg["muon energy (with threshold) vs. beamdump angle (1)"].Integral())
    hdict_bkg["muon energy (with threshold) vs. beamdump angle (1)"].Draw("colz")
    u10.Print("bg_muonE_thresh_angle1.png")
    
    #bkg plot (11)
    u11 = ROOT.TCanvas()          
    hdict_bkg["muon energy (with threshold) vs. beamdump angle (2)"].Scale(1/hdict_bkg["muon energy (with threshold) vs. beamdump angle (2)"].Integral())
    hdict_bkg["muon energy (with threshold) vs. beamdump angle (2)"].Draw("colz")
    u11.Print("bg_muonE_thresh_angle2.png")
    
    #bkg plot (12)
    u12 = ROOT.TCanvas()          
    hdict_bkg["muon energy (with threshold) vs. beamdump angle (3)"].Scale(1/hdict_bkg["muon energy (with threshold) vs. beamdump angle (3)"].Integral())
    hdict_bkg["muon energy (with threshold) vs. beamdump angle (3)"].Draw("colz")
    u12.Print("bg_muonE_thresh_angle3.png")
    
    #bkg plot (13)
    u13 = ROOT.TCanvas()
    hdict_bkg["omega vs. q (bkg)"].Scale(norm)
    hdict_bkg["omega vs. q (bkg)"].Draw("colz")
    u13.Print("omegaq_bkg.png")
    
    #bkg plot (14)
    u16 = ROOT.TCanvas()
    hdict_bkg["omega q (reco.)"].Scale(norm)
    hdict_bkg["omega q (reco.)"].Draw("colz")
    u16.Print("omegaq_bkg_reco.png")
     
    #bkg plot (15)
    u14 = ROOT.TCanvas("omega vs. q (true bkg)", "multipads", 1000, 500)
    u14.Divide(5,2,0,0)
    hdict_bkg["omega v. q (truth panel)"].Scale(norm)
    maxi1 = hdict_bkg["omega v. q (truth panel)"].GetMaximum()
    for i in range (0, 10):
        u14.cd(i+1)
        hq_bkg1 = hdict_bkg["omega v. q (truth panel)"].ProjectionX('name_%d' % i, i, i)
        hq_bkg1.GetYaxis().SetRangeUser(0., maxi1*1.2)
        for b in range (1, hq_bkg1.GetNbinsX()+1):
            content = hq_bkg1.GetBinContent(b)
            error = (content)**0.5
            hq_bkg1.SetBinError(b, error)
        hq_bkg1.SetLineWidth(1)
        hq_bkg1.SetMarkerSize(0.5)   
        hq_bkg1.Draw("hist e")                  
    u14.Print("omegaq_bkg_true_panel.png")
    
    #bkg plot (16)
    u15 = ROOT.TCanvas("omega vs. q (reco bkg)", "multipads", 1000, 500)
    u15.Divide(5,2,0,0)
    hdict_bkg["omega and q (reco panel)"].Scale(norm)
    max2 = hdict_bkg["omega and q (reco panel)"].GetMaximum()
    for i in range (0, 10):
        u15.cd(i+1)
        omegaq_bkg2 = hdict_bkg["omega and q (reco panel)"].ProjectionX('name_%d' % i, i, i)
        omegaq_bkg2.GetYaxis().SetRangeUser(0., max2*1.2)
        for b in range (1, omegaq_bkg2.GetNbinsX()+1):
            content = omegaq_bkg2.GetBinContent(b)
            error = (content)**0.5
            omegaq_bkg2.SetBinError(b, error)
        omegaq_bkg2.SetLineWidth(1)
        omegaq_bkg2.SetMarkerSize(0.5)   
        omegaq_bkg2.Draw("hist e")            
    u15.Print("omegaq_reco_bkg_panel.png")
    
    #stacked plot 1
    u17 = ROOT.TCanvas("omega vs. q (true signal and background stack)", "multipads", 1000, 500)
    u17.Divide(5,2,0,0)
    hdict_sig["omega q stack (truth)"].Scale(norm)
    hdict_bkg["omega vs. q (stack)"].Scale(norm)
    maxi1 = hdict_sig["omega q stack (truth)"].GetMaximum()
    maxi2 = hdict_bkg["omega vs. q (stack)"].GetMaximum()
    true_stack = {}
    for i in range (0, 10):
        u17.cd(i+1)
        true_stack[i] = ROOT.THStack("true_stack", " " )
        hq_bkg1 = hdict_bkg["omega vs. q (stack)"].ProjectionX('name_%d' % i, i, i)
        hq_bkg1.SetFillColor(ROOT.kBlue)
        hq_bkg1.SetLineWidth(1)
        hq_bkg1.SetMarkerSize(0.5)  
        for b in range (1, hq_bkg1.GetNbinsX()+1):
            content = hq_bkg1.GetBinContent(b)
            error = (content)**0.5
            hq_bkg1.SetBinError(b, error)
        true_stack[i].Add(hq_bkg1)
        hq_sig1 = hdict_sig["omega q stack (truth)"].ProjectionX('name_%c' % i, i, i)
        hq_sig1.SetFillColor(ROOT.kRed)
        hq_sig1.SetLineWidth(1)
        hq_sig1.SetMarkerSize(0.5)  
        for b in range (1, hq_sig1.GetNbinsX()+1):
            content = hq_sig1.GetBinContent(b)
            error = (content)**0.5
            hq_sig1.SetBinError(b, error)
        true_stack[i].Add(hq_sig1)  
        true_stack[i].Draw("hist e")
        true_stack[i].SetMaximum(maxi2*1.2)       
        true_stack[i].Draw("hist e")
    u17.Print("omegaq_truth_stack.png")
    
    #stacked plot 2
    y1 = ROOT.TCanvas("omega vs. q (reco. signal and background stack)", "multipads", 1000, 500)
    y1.Divide(5,2,0,0)
    hdict_sig["omega q stack (reco)"].Scale(norm)
    hdict_bkg["omega q (reco stack)"].Scale(norm)
    maximum_1 = hdict_bkg["omega q (reco stack)"].GetMaximum()
    reco_stack = {}
    for i in range (0, 10):
        y1.cd(i+1)
        reco_stack[i] = ROOT.THStack("reco stack", " " )
        hq_bkg_reco1 = hdict_bkg["omega q (reco stack)"].ProjectionX('name_%d' % i, i, i)
        hq_bkg_reco1.SetFillColor(ROOT.kBlue)
        hq_bkg_reco1.SetLineWidth(1)
        hq_bkg_reco1.SetMarkerSize(0.5)  
        for b in range (1, hq_bkg_reco1.GetNbinsX()+1):
            content = hq_bkg_reco1.GetBinContent(b)
            error = (content)**0.5
            hq_bkg_reco1.SetBinError(b, error)
        reco_stack[i].Add(hq_bkg_reco1)
        hq_sig_reco1 = hdict_sig["omega q stack (reco)"].ProjectionX('name_%c' % i, i, i)
        for b in range (1, hq_sig_reco1.GetNbinsX()+1):
            content = hq_sig_reco1.GetBinContent(b)
            error = (content)**0.5
            hq_sig_reco1.SetBinError(b, error)
        hq_sig_reco1.SetFillColor(ROOT.kRed)
        hq_sig_reco1.SetLineWidth(1)
        hq_sig_reco1.SetMarkerSize(0.5)  
        reco_stack[i].Add(hq_sig_reco1)
        reco_stack[i].Draw("hist e")
        reco_stack[i].SetMaximum(maximum_1*1.2)
        reco_stack[i].Draw("hist e")   
    y1.Print("omegaq_reco_stack.png")
    


