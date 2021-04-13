import ROOT
# "Batch mode" means do not display graphics on the screen
ROOT.gROOT.SetBatch(1) 

if __name__ == "__main__":

    # Open the TFile that contains the TTree of events
    tf = ROOT.TFile( "out.root" )
    # Get the TTree object that contains the neutrino events
    tree = tf.Get( "tree" )

    # Declare some histograms
    hOmega = ROOT.TH1D( "omega", ";#omega = E_{#nu} - E_{#mu} (MeV)", 136, 0., 136. )

    # Loop through the neutrino interactions
    for entry in tree:
        # evt, globalEvt, spill = book-keeping indices
        # Enu = the neutrino energy
        # nuPdg = the PDG code of the neutrino (14 = numu, 12 = nue, -14 = anti-numu)
        # x, y, Q2, W = the event kinematics: Bjorken x, inelasticity, four-momentum transfer, hadronic invariant mass
        # tgt = the PDG code of the struck nucleus
        # struck = PDG code of the struck nucleon within the nucleus
        # scat = a code for what type of interaction (1 = QE, 3 = DIS, 4 = RES, 5 = COH, 7 = nu+e, 10 = MEC)
        # nfsp = number of final state particles
        # Branches for each FS particle:
        #    pdg = PDG code of final-state particle
        #    px, py, pz = Momentum of final-state particles (z = neutrino)
        #    E = total energy of final-state particle (including rest mass)

        # is the event charged-current?
        is_cc = False
        muon_E = None

        # Loop over final-state particles
        for i in range(entry.nfsp):
            if entry.pdg[i] == 13:
                is_cc = True
                muon_E = entry.E[i]

        # Make a cut on CC events
        if is_cc:

            # Fill the histogram
            hOmega.Fill( (entry.Enu - muon_E)*1000. ) # convert to MeV

    # Make an output file for the histograms
    tfout = ROOT.TFile( "plots.root", "RECREATE" ) # overwrite the file if it already exists
    hOmega.Write()

    # Make a canvas to draw stuff on
    c = ROOT.TCanvas()
    hOmega.Draw("hist") # draw it as a histogram
    c.Print( "omega.png" ) # write it out to a file
