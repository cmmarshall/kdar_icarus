void flatten_reco()
{

  const char * recf = "/icarus/data/users/mmr/CAFMaker/v09_20_00/NuMIBeamWindow/ICARUS_prod_2020A_00_numioffaxis_v09_10_01_reco2/ICARUS_prod_2020A_00_numioffaxis_v09_10_01_reco2_CAFMaker_out_flat.root";
  TFile * tf = new TFile( recf );
  TTree * tree = (TTree*) tf->Get( "recTree" );

  TFile * outFile = new TFile( "background.root", "RECREATE" );
  TTree * outTree = new TTree( "tree", "output tree" );

  static const int NNU = 20;
  static const int NPRIM = 100;

  // CAF variables
  int nnu;
  float nuE[NNU], nuPxi[NNU], nuPyi[NNU], nuPzi[NNU];
  int nuPdgCode[NNU], nuScat[NNU], nPrim[NNU];
  float nuQ2[NNU], bjX[NNU], intY[NNU], nuW[NNU]; 
  int primPdg[NPRIM];
  float primPx[NPRIM], primPy[NPRIM], primPz[NPRIM], primE[NPRIM];

  tree->SetBranchAddress( "rec.mc.nu..length", &nnu );
  tree->SetBranchAddress( "rec.mc.nu.pdg", &nuPdgCode );
  tree->SetBranchAddress( "rec.mc.nu.E", &nuE );
  tree->SetBranchAddress( "rec.mc.nu.momentum.x", &nuPxi );
  tree->SetBranchAddress( "rec.mc.nu.momentum.y", &nuPyi );
  tree->SetBranchAddress( "rec.mc.nu.momentum.z", &nuPzi );
  tree->SetBranchAddress( "rec.mc.nu.Q2", &nuQ2 );
  tree->SetBranchAddress( "rec.mc.nu.w", &nuW );
  tree->SetBranchAddress( "rec.mc.nu.bjorkenX", &bjX );
  tree->SetBranchAddress( "rec.mc.nu.inelasticityY", &intY );
  tree->SetBranchAddress( "rec.mc.nu.genie_intcode", &nuScat );
  //tree->SetBranchAddress( "rec.mc.nu.hitnuc", hitNuc );
  tree->SetBranchAddress( "rec.mc.nu.prim..length", &nPrim );
  tree->SetBranchAddress( "rec.mc.nu.prim.pdg", &primPdg );
  tree->SetBranchAddress( "rec.mc.nu.prim.genE", &primE );
  tree->SetBranchAddress( "rec.mc.nu.prim.genp.x", &primPx );
  tree->SetBranchAddress( "rec.mc.nu.prim.genp.y", &primPy );
  tree->SetBranchAddress( "rec.mc.nu.prim.genp.z", &primPz );


  // Output tree
  int evtNo, spillNo, globalEvtNo, nuPdg, nfsp, scat;
  float Enu, x, y, Q2, W, nuPx, nuPy, nuPz;

  int pdg[100];
  float px[100], py[100], pz[100], E[100];

  outTree->Branch( "evt", &evtNo, "evt/I" );
  outTree->Branch( "globalEvt", &globalEvtNo, "globalEvt/I" );
  outTree->Branch( "spill", &spillNo, "spill/I" );
  outTree->Branch( "Enu", &Enu, "Enu/F" );
  outTree->Branch( "nuPdg", &nuPdg, "nuPdg/I" );
  outTree->Branch( "nuPx", &nuPx, "nuPx/F" );
  outTree->Branch( "nuPy", &nuPy, "nuPy/F" );
  outTree->Branch( "nuPz", &nuPz, "nuPz/F" );
  outTree->Branch( "x", &x, "x/F" );
  outTree->Branch( "y", &y, "y/F" );
  outTree->Branch( "Q2", &Q2, "Q2/F" );
  outTree->Branch( "W", &W, "W/F" );
  //outTree->Branch( "tgt", &tgt, "tgt/I" );
  //outTree->Branch( "struck", &struck, "struck/I" );
  outTree->Branch( "scat", &scat, "scat/I" );
  outTree->Branch( "nfsp", &nfsp, "nfsp/I" );
  outTree->Branch( "pdg", pdg, "pdg[nfsp]/I" );
  outTree->Branch( "px", px, "px[nfsp]/F" );
  outTree->Branch( "py", py, "py[nfsp]/F" );
  outTree->Branch( "pz", pz, "pz[nfsp]/F" );
  outTree->Branch( "E", E, "E[nfsp]/F" );

  globalEvtNo = 0;

  const int N = tree->GetEntries();
  for( int ii = 0; ii < N; ++ii ) {
    if( ii % 1000 == 0 ) printf( "Spill %d of %d...\n", ii, N );
    tree->GetEntry(ii);

    spillNo = ii;
    evtNo = 0;
    int primIdx = 0;

    // loop over neutrinos
    for( int i = 0; i < nnu; ++ i ) {
      Enu = nuE[i];
      nuPdg = nuPdgCode[i];
      nuPx = nuPxi[i];
      nuPy = nuPyi[i];
      nuPz = nuPzi[i];
      x = bjX[i];
      y = intY[i];
      W = nuW[i];
      Q2 = nuQ2[i];
      scat = nuScat[i];
      nfsp = nPrim[i];
      for( int p = 0; p < nfsp; ++p ) {
        pdg[p] = primPdg[primIdx+p];
        E[p] = primE[primIdx+p];
        px[p] = primPx[primIdx+p];
        py[p] = primPy[primIdx+p];
        pz[p] = primPz[primIdx+p];
        //printf( "     primary %d p = (%1.1f, %1.1f, %1.1f) E = %1.2f\n", pdg[p], px[p], py[p], pz[p], E[p] );
      }
      primIdx += nfsp;    
      outTree->Fill();
      evtNo++;
      globalEvtNo++;
    }
  }

  outTree->Write();

}
