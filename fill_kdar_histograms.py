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
                hdict_sig["muon energy and angle"].Fill((Emu_rec)*1000, (theta))
                hdict_sig["rec-muon energy and true muon E"].Fill((muon_ke)*1000, (Emu_rec)*1000)
                hdict_sig["fractional residual"].Fill((Emu_rec - (muon_ke))/(muon_ke), (muon_ke))
                hdict_sig["omega(rec) and q"].Fill((entry.Enu - (Emu_rec + 0.105))*1000,  q*1000)
                hdict_sig["omega (MeV) vs. q (MeV/C)"].Fill(omega*1000, q*1000)
                hdict_sig["omega (MeV) vs. q (MeV/c) (rec)"].Fill((entry.Enu - (Emu_rec + 0.105))*1000,  q*1000)
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
                if muon_E < 0.236:
                    hdict_sig["muon energy (with threshold) vs. beamdump angle (sig 1)"].Fill(muon_E*1000, math.cos(reco1_sum.Angle(beamdump_nu)))
                    hdict_sig["muon energy (with threshold) vs. beamdump angle (sig 2)"].Fill(muon_E*1000, math.cos(reco2_sum.Angle(beamdump_nu)))
                    hdict_sig["muon energy (with threshold) vs. beamdump angle (sig 3)"].Fill(muon_E*1000, math.cos(reco3_sum.Angle(beamdump_nu)))
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
    return n_signal

if __name__ == "__main__":
               
    #background plots
    hdict_bkg = {}
    hdict_bkg["neutrino angle wrst. Z"] = ROOT.TH2D("neutrino angle wrst. Z", ";#theta_{#nu (XZ)} ;#theta_{#nu (YZ)}", 136, 0., 1., 136, 0., 1.)   
    hdict_bkg["reco neutrino angle"] =  ROOT.TH2D("reco neutrino angle", ";#theta_{#nu (XZ)} ;#theta_{#nu (YZ)}", 136, -1.6, 1.6, 136, -1.6, 1.6)
    hdict_bkg["reco. neutrino angle (no neutron)"] = ROOT.TH2D("reco. neutrino angle (no neutron)", ";#theta_{#nu (XZ)} ;#theta_{#nu (YZ)}", 136, -1.6, 1.6, 136, -1.6, 1.6)
    hdict_bkg["reco. neutrino angle (proton threshold)"] = ROOT.TH2D("reco. neutrino angle (proton threshold)", ";#theta_{#nu (XZ)} ;#theta_{#nu (YZ)}", 136, -1.6, 1.6, 136, -1.6, 1.6)
    hdict_bkg["beam dump angle (reco. 1)"] = ROOT.TH1D("beam dump angle (reco. 1)", ";#theta_{#nu}", 136, -1., 1.)
    hdict_bkg["beam dump angle (reco. 2)"] = ROOT.TH1D("beam dump angle (reco. 2)", ";#theta_{#nu}", 136, -1., 1.)
    hdict_bkg["beam dump angle (reco. 3)"] = ROOT.TH1D("beam dump angle (reco. 3)", ";#theta_{#nu}", 136, -1., 1.)
    hdict_bkg["beam dump angle (muon thresh reco. 1)"] = ROOT.TH1D("beam dump angle (muon thresh reco. 1)", ";#theta_{#nu}", 136, -1., 1.)
    hdict_bkg["beam dump angle (muon thresh reco. 2)"] = ROOT.TH1D("beam dump angle (muon thresh reco. 2)", ";#theta_{#nu}", 136, -1., 1.)
    hdict_bkg["beam dump angle (muon thresh reco. 3)"] = ROOT.TH1D("beam dump angle (muon thresh reco. 3)", ";#theta_{#nu}", 136, -1., 1.)
    hdict_bkg["muon energy vs. beam dump angle (1) "] = ROOT.TH2D("muon energy vs. beam dump angle (1) ", ";E_{#mu} MeV ;cos(#theta_{#nu})", 136, 80., 250., 136, -1., 1.)
    hdict_bkg["muon energy vs. beam dump angle (2) "] = ROOT.TH2D("muon energy vs. beam dump angle (2) ", ";E_{#mu} MeV ;cos(#theta_{#nu})", 136, 80., 250., 136, -1., 1.)
    hdict_bkg["muon energy vs. beam dump angle (3)"] = ROOT.TH2D("muon energy vs. beam dump angle (3)", ";E_{#mu} MeV ;cos(#theta_{#nu})", 136, 80., 250., 136, -1., 1.)
    hdict_bkg["muon energy (with threshold) vs. beamdump angle (1)"] = ROOT.TH2D("muon energy (with threshold) vs. beamdump angle (1)", ";E_{#mu} MeV ;cos(#theta_{#nu})", 136, 80, 250., 136, -1., -1.)
    hdict_bkg["muon energy (with threshold) vs. beamdump angle (2)"] = ROOT.TH2D("muon energy (with threshold) vs. beamdump angle (2)", ";E_{#mu} MeV ;cos(#theta_{#nu})", 136, 80, 250., 136, -1., -1.)
    hdict_bkg["muon energy (with threshold) vs. beamdump angle (3)"] = ROOT.TH2D("muon energy (with threshold) vs. beamdump angle (3)", ";E_{#mu} MeV ;cos(#theta_{#nu})", 136, 80, 250., 136, -1., -1.)
    
    #signal plots   
    hdict_sig = {}
    hdict_sig["omega"] = ROOT.TH1D("omega",";#omega = E_{#nu} - E_{#mu} (MeV)", 136, 0., 136. )
    hdict_sig["Emu"] = ROOT.TH1D( "E_{#mu}", ";E_{#mu} (MeV)", 136, 0., 236. )
    hdict_sig["omega and q"] = ROOT.TH2D( "omega and q", ";#omega = E_{#nu} - E_{#mu} (MeV), ;q=|p_{#nu}-p_{#mu}| (MeV/c)", 136, 0., 140., 136, 0., 450. )
    hdict_sig["omega and muon angle"] = ROOT.TH2D( "omega and muon angle", ";#omega = E_{#nu} - E_{#mu} (MeV) ;muon angle", 136, 0., 140., 136, 0., 5. )
    hdict_sig["muon energy and angle"] = ROOT.TH2D( "Muon energy and angle ", ";E_{#mu} (MeV) ;muon angle (rad)", 136, 0., 300., 136, 0., 5.)
    hdict_sig["reconstructed muon energy and muon angle"] = ROOT.TH2D( "reconstructed muon energy and muon angle", ";E_{#mu (rec)} = x * E_{#mu} ;muon angle (rad)", 136, 0 , 300., 136, 0. , 5.)
    hdict_sig["rec-muon energy and true muon E"] = ROOT.TH2D( "rec-muon energy and true muon E", "; E_{#mu} (MeV) ;E_{#mu (rec)} = x * E_{#mu} (MeV) ", 136, 0. , 350., 136, 0. , 350.)
    hdict_sig["fractional residual"] = ROOT.TH2D("frational residual", ";(E_{#mu (rec)}-E_{#mu})/E_{#mu}, ;E_{#mu}(GeV) ", 136, -1., 1., 136, 0., 0.5)
    hdict_sig["omega(rec) and q"] = ROOT.TH2D("omega (rec) and q", ";#omega = E_{#nu} - E_{#mu (rec)} (MeV) ;q=|p_{#nu}-p_{#mu}| (MeV/c)", 136, 0., 140., 136, 0., 450.)
    hdict_sig["omega (MeV) vs. q (MeV/C)"] = ROOT.TH2D("omega (MeV) vs. q (MeV/c)", ";#omega = E_{#nu}-E_{#mu} (MeV) ;q =|p_{#nu}-p_{#mu}| (MeV/c)", 136, 0., 140., 10, 0., 500.)
    hdict_sig["omega (MeV) vs. q (MeV/c) (rec)"] = ROOT.TH2D("omega (MeV) vs. q (MeV/c) (rec)", ";#omega = E_{#nu}-E_{#mu} (MeV) ;q =|p_{#nu}-p_{#mu}| (MeV/c)", 50, 0., 140., 10, 0., 500.)
    hdict_sig["true neutrino angles"] = ROOT.TH2D("true neutrino angles", ";#theta_{#nu (XZ)} ;#theta_{#nu (YZ)}", 136, 0., 1., 136, 0., 1.)
    hdict_sig["reconstructed neutrino angle"] = ROOT.TH2D("reconstructed neutrino angle", ";#theta_{#nu (XZ)} ;#theta_{#nu (YZ)}", 136, -1.6, 1.6, 136, -1.6, 1.6)
    hdict_sig["reco.neutrino angle (no neutron)"] = ROOT.TH2D("reco.neutrino angle (no neutron)", ";#theta_{#nu (XZ)} ;#theta_{#nu (YZ)}", 136, -1.6, 1.6, 136, -1.6, 1.6)
    hdict_sig["reco. neutrino angle (signal proton threshold)"] = ROOT.TH2D("reco. neutrino angle (signal proton threshold)", ";#theta_{#nu (XZ)} ;#theta_{#nu (YZ)}", 136, -1.6, 1.6, 136, -1.6, 1.6)
    hdict_sig["reco. assumption 1 signal wrst. beam dump"] = ROOT.TH1D("reco. assumption 1 signal wrst. beam dump" , ";#theta_{#nu}", 136, -1., 1.)
    hdict_sig["reco. assumption 2 wrst. beam dump"] = ROOT.TH1D("reco. assumption 2 wrst. beam dump", ";#theta_{#nu}", 136, -1., 1.)
    hdict_sig["reco. assumption 3 wrst. beam dump"] = ROOT.TH1D("reco. assumption 3 wrst. beam dump", ";#theta_{#nu}", 136, -1., 1.)
    hdict_sig["muon energy vs. beam dump angle (sig 1) "] = ROOT.TH2D("muon energy vs. beam dump angle (sig 1) ", ";E_{#mu} MeV ;cos(#theta_{#nu})", 136, 80., 250., 136, -1., 1.)
    hdict_sig["muon energy vs. beam dump angle (sig 2) "] = ROOT.TH2D("muon energy vs. beam dump angle (sig 2) ", ";E_{#mu} MeV ;cos(#theta_{#nu})", 136, 80., 250., 136, -1., 1.)
    hdict_sig["muon energy vs. beam dump angle (sig 3) "] = ROOT.TH2D("muon energy vs. beam dump angle (sig 3) ", ";E_{#mu} MeV ;cos(#theta_{#nu})", 136, 80., 250., 136, -1., 1.)
    hdict_sig["muon energy (with threshold) vs. beamdump angle (sig 1)"] = ROOT.TH2D("muon energy (with threshold) vs. beamdump angle (sig 1)", ";E_{#mu} MeV ;cos(#theta_{#nu})", 136, 80, 250., 136, -1., -1.)
    hdict_sig["muon energy (with threshold) vs. beamdump angle (sig 2)"] = ROOT.TH2D("muon energy (with threshold) vs. beamdump angle (sig 2)", ";E_{#mu} MeV ;cos(#theta_{#nu})", 136, 80, 250., 136, -1., -1.)
    hdict_sig["muon energy (with threshold) vs. beamdump angle (sig 3)"] = ROOT.TH2D("muon energy (with threshold) vs. beamdump angle (sig 3)", ";E_{#mu} MeV ;cos(#theta_{#nu})", 136, 80, 250., 136, -1., -1.)

    # run it
    tf = ROOT.TFile( filename )
    tree = tf.Get( "tree" )
    n_signal = TTree_loop( tree, hdict_sig, hdict_bkg )

    norm = 1.*sig_per_yr / n_signal
          
    # Make an output file for the histograms
    tfout = ROOT.TFile( "plots.root", "RECREATE" ) # overwrite the file if it already exists
    hdict_sig["omega"].Write()
    #signal drawings
    c = ROOT.TCanvas()
    hdict_sig["omega"].Scale(norm)
    for b in range (1, hdict_sig["omega"].GetNbinsX()+1):
        content = hdict_sig["omega"].GetBinContent(b)
        error = content**0.5
        hdict_sig["omega"].SetBinError(b, error)
    hdict_sig["omega"].Draw("hist E") # draw it as a histogram
    c.Print( "omega.png" ) # write it out to a file
    
    hdict_sig["Emu"].Scale(norm)
    for b in range (1, hdict_sig["Emu"].GetNbinsX()+1):
        content = hdict_sig["Emu"].GetBinContent(b)
        error = content**0.5
        hdict_sig["Emu"].SetBinError(b, error)
    hdict_sig["Emu"].Draw("hist E")
    c.Print("emu.png")
    
    hdict_sig["omega and q"].Scale(norm)
    hdict_sig["omega and q"].Draw("colz")
    c.Print("omega-q.png")
    
    hdict_sig["omega and muon angle"].Scale(norm)
    hdict_sig["omega and muon angle"].Draw("colz")
    c.Print("omega-angle.png")
    
    hdict_sig["muon energy and angle"].Scale(norm)
    hdict_sig["muon energy and angle"].Draw("colz")
    c.Print("emu-angle.png")
       
    hdict_sig["reconstructed muon energy and muon angle"].Scale(norm)
    hdict_sig["reconstructed muon energy and muon angle"].Draw("colz")
    c.Print("emu.rec-angle.png")
    
    hdict_sig["rec-muon energy and true muon E"].Scale(norm)
    hdict_sig["rec-muon energy and true muon E"].Draw("colz")
    c.Print("emur-emu.png")
    
    hdict_sig["fractional residual"].Scale(norm)
    hdict_sig["fractional residual"].Draw("colz")
    c.Print("erec.png")
    
    hdict_sig["omega(rec) and q"].Scale(norm)
    hdict_sig["omega(rec) and q"].Draw("colz")
    c.Print("omegaq.rec.png")
    
    c11 = ROOT.TCanvas("omega and q", "multipads", 1000, 500)
    c11.Divide(5, 2, 0, 0)
    hdict_sig["omega (MeV) vs. q (MeV/C)"].Scale(norm)
    maximum = hdict_sig["omega (MeV) vs. q (MeV/C)"].GetMaximum()
    for i in range (0, 10):
        c11.cd(i+1)
        hq = hdict_sig["omega (MeV) vs. q (MeV/C)"].ProjectionX('name_%d' % i, i, i)
        hq.GetYaxis().SetRangeUser(0., maximum*1.2)
        for b in range (1, hq.GetNbinsX()+1):
            content = hq.GetBinContent(b)
            error = ((content)**0.5)
            hq.SetBinError(b, error)
        hq.SetMarkerSize(0.5)
        hq.SetLineWidth(1)
        hq.Draw("hist E")
    c11.Print("omegaq-multiloop.png")
    
    
    c2 = ROOT.TCanvas("omega vs. q (rec)", "multipads", 1000, 500)
    c2.Divide(5,2,0,0)
    hdict_sig["omega (MeV) vs. q (MeV/c) (rec)"].Scale(norm)
    maximum = hdict_sig["omega (MeV) vs. q (MeV/c) (rec)"].GetMaximum()
    for i in range (0, 10):
        c2.cd(i+1)
        hqrec = hdict_sig["omega (MeV) vs. q (MeV/c) (rec)"].ProjectionX('name_%c' % i, i, i)
        hqrec.GetYaxis().SetRangeUser(0., maximum*1.2)
        for b in range (1, hqrec.GetNbinsX()+1):
            content = hqrec.GetBinContent(b)
            error = (content)**0.5
            hqrec.SetBinError(b, error)
        hqrec.SetLineWidth(1)
        hqrec.SetMarkerSize(0.5)
        hqrec.Draw("hist E")  
    c2.Print("omegaq-rec-multi.png")
           
    hdict_sig["true neutrino angles"].Scale(norm)
    hdict_sig["true neutrino angles"].Draw("colz")
    c.Print("true_signal_angles.png")
           
    hdict_sig["reconstructed neutrino angle"] .Scale(norm)
    hdict_sig["reconstructed neutrino angle"] .Draw("colz")
    c.Print("recoFS_signal_angle.png")
           
    hdict_sig["reco.neutrino angle (no neutron)"].Scale(norm)
    hdict_sig["reco.neutrino angle (no neutron)"].Draw("colz")
    c.Print("recoFS_signal_neu.png")
           
    hdict_sig["reco. neutrino angle (signal proton threshold)"].Scale(norm)
    hdict_sig["reco. neutrino angle (signal proton threshold)"].Draw("colz")
    c.Print("signal_proton_thresh.png")
           
    c12 = ROOT.TCanvas()
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
    c12.Print("signal_beamdump_angles.png")
           
    hdict_sig["muon energy vs. beam dump angle (sig 1) "].Scale(norm)
    hdict_sig["muon energy vs. beam dump angle (sig 1) "].Draw("colz")
    c.Print("sig_muonE_angle1.png")
           
    hdict_sig["muon energy vs. beam dump angle (sig 2) "].Scale(norm)
    hdict_sig["muon energy vs. beam dump angle (sig 2) "].Draw("colz")
    c.Print("sig_muonE_angle2.png")
           
    hdict_sig["muon energy vs. beam dump angle (sig 3) "].Scale(norm)
    hdict_sig["muon energy vs. beam dump angle (sig 3) "].Draw("colz")
    c.Print("sig_muonE_angle3.png")
           
    hdict_sig["muon energy (with threshold) vs. beamdump angle (sig 1)"].Scale(1/ hdict_sig["muon energy (with threshold) vs. beamdump angle (sig 1)"].Integral())
    hdict_sig["muon energy (with threshold) vs. beamdump angle (sig 1)"].Draw("colz")
    c.Print("sig_muonE_thresh_angle1.png")
           
    hdict_sig["muon energy (with threshold) vs. beamdump angle (sig 2)"].Scale(1/hdict_sig["muon energy (with threshold) vs. beamdump angle (sig 2)"].Integral())
    hdict_sig["muon energy (with threshold) vs. beamdump angle (sig 2)"].Draw("colz")
    c.Print("sig_muonE_thresh_angle2.png")
           
    hdict_sig["muon energy (with threshold) vs. beamdump angle (sig 3)"].Scale(1/hdict_sig["muon energy (with threshold) vs. beamdump angle (sig 3)"].Integral())
    hdict_sig["muon energy (with threshold) vs. beamdump angle (sig 3)"].Draw("colz")
    c.Print("sig_muonE_thresh_angle3.png")
    
    #background drawings   
    u = ROOT.TCanvas()
    hdict_bkg["neutrino angle wrst. Z"].Scale(norm)
    hdict_bkg["neutrino angle wrst. Z"].Draw("colz")
    u.Print("true_background_angle.png")
               
    hdict_bkg["reco neutrino angle"].Scale(norm)
    hdict_bkg["reco neutrino angle"].Draw("colz")
    u.Print("recoFS_background_angle.png")
               
    hdict_bkg["reco. neutrino angle (no neutron)"].Scale(norm)
    hdict_bkg["reco. neutrino angle (no neutron)"].Draw("colz")
    u.Print("reco_background_neu.png")
               
    hdict_bkg["reco. neutrino angle (proton threshold)"].Scale(norm)
    hdict_bkg["reco. neutrino angle (proton threshold)"].Draw("colz")
    u.Print("background_proton_thresh.png")
           
    c13 = ROOT.TCanvas()
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
    c13.Print("background_beamdump_angles.png")
               
    c16 = ROOT.TCanvas()
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
    c16.Print("background_beamdump_thresh.png")
               
    hdict_bkg["muon energy vs. beam dump angle (1) "].Scale(norm)
    hdict_bkg["muon energy vs. beam dump angle (1) "].Draw("colz")
    u.Print("bg_muonE_angle1.png")
               
    hdict_bkg["muon energy vs. beam dump angle (2) "].Scale(norm)
    hdict_bkg["muon energy vs. beam dump angle (2) "].Draw("colz")
    u.Print("bg_muonE_angle2.png")
               
    hdict_bkg["muon energy vs. beam dump angle (3)"].Scale(norm)
    hdict_bkg["muon energy vs. beam dump angle (3)"].Draw("colz")
    u.Print("bg_muonE_angle3.png")
               
    hdict_bkg["muon energy (with threshold) vs. beamdump angle (1)"].Scale(1/hdict_bkg["muon energy (with threshold) vs. beamdump angle (1)"].Integral())
    hdict_bkg["muon energy (with threshold) vs. beamdump angle (1)"].Draw("colz")
    u.Print("bg_muonE_thresh_angle1.png")
               
    hdict_bkg["muon energy (with threshold) vs. beamdump angle (2)"].Scale(1/hdict_bkg["muon energy (with threshold) vs. beamdump angle (2)"].Integral())
    hdict_bkg["muon energy (with threshold) vs. beamdump angle (2)"].Draw("colz")
    u.Print("bg_muonE_thresh_angle2.png")
               
    hdict_bkg["muon energy (with threshold) vs. beamdump angle (3)"].Scale(1/hdict_bkg["muon energy (with threshold) vs. beamdump angle (3)"].Integral())
    hdict_bkg["muon energy (with threshold) vs. beamdump angle (3)"].Draw("colz")
    u.Print("bg_muonE_thresh_angle3.png")
               
               
          












