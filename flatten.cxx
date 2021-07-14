#include <string>
#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TIterator.h>
#include <TLorentzVector.h>
#include <TH2.h>
#include <TH3.h>

#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Ntuple/NtpMCFormat.h"
#include "Framework/Ntuple/NtpMCTreeHeader.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Utils/CmdLnArgParser.h"
#include "Framework/Conventions/Units.h"

using namespace genie;

std::string fOutFileName;
std::string fInFileDir;
std::string fInFilePrefix;
std::string fAngleFile;
int first_run, last_run;
bool dump_dir = false;
bool rotate = false;

//___________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{

  CmdLnArgParser parser(argc,argv);

  // ghep file name
  if( parser.OptionExists('i') ) {
    fInFileDir = parser.ArgAsString('i');
  } else {
    exit(1);
  }

  if( parser.OptionExists('p') ) {
    fInFilePrefix = parser.ArgAsString('p');
  } else {
    exit(1);
  }

  if( parser.OptionExists('f') ) {
    first_run = parser.ArgAsInt('f');
  } else {
    exit(1);
  }

  if( parser.OptionExists('l') ) {
    last_run = parser.ArgAsInt('l');
  } else {
    exit(1);
  }

  if( parser.OptionExists('o') ) {
    fOutFileName = parser.ArgAsString('o');
  } else {
    exit(1);
  }

  if( parser.OptionExists('a') ) {
    fAngleFile = parser.ArgAsString('a');
    rotate = true;
  } else {
    rotate = false;
  }

  if( parser.OptionExists('d') ) {
    rotate = true;
    dump_dir = true;
  }


}


//___________________________________________________________________
int main(int argc, char ** argv)
{

  GetCommandLineArgs (argc, argv);

  TFile * outFile = new TFile( fOutFileName.c_str(), "RECREATE" );
  TTree * tree = new TTree( "tree", "output ttree" );

  int evtNo, spillNo, globalEvtNo, nuPdg, nfsp, tgt, scat;
  float Enu, nuPx, nuPy, nuPz, x, y, Q2, W;

  int pdg[100];
  float px[100], py[100], pz[100], E[100];

  tree->Branch( "evt", &evtNo, "evt/I" );
  tree->Branch( "globalEvt", &globalEvtNo, "globalEvt/I" );
  tree->Branch( "spill", &spillNo, "spill/I" );
  tree->Branch( "Enu", &Enu, "Enu/F" );
  tree->Branch( "nuPdg", &nuPdg, "nuPdg/I" );
  tree->Branch( "nuPx", &nuPx, "nuPx/F" );
  tree->Branch( "nuPy", &nuPy, "nuPy/F" );
  tree->Branch( "nuPz", &nuPz, "nuPz/F" );
  tree->Branch( "x", &x, "x/F" );
  tree->Branch( "y", &y, "y/F" );
  tree->Branch( "Q2", &Q2, "Q2/F" );
  tree->Branch( "W", &W, "W/F" );
  tree->Branch( "tgt", &tgt, "tgt/I" );
  tree->Branch( "scat", &scat, "scat/I" );
  tree->Branch( "nfsp", &nfsp, "nfsp/I" );
  tree->Branch( "pdg", pdg, "pdg[nfsp]/I" );
  tree->Branch( "px", px, "px[nfsp]/F" );
  tree->Branch( "py", py, "py[nfsp]/F" );
  tree->Branch( "pz", pz, "pz[nfsp]/F" );
  tree->Branch( "E", E, "E[nfsp]/F" );

  //-- open the ROOT file and get the TTree & its header
  TChain * gTree = new TChain( "gtree", "gtree" );

  for( int i = first_run; i <= last_run; ++i ) {
    gTree->Add( Form("%s/%s.%d.ghep.root",fInFileDir.c_str(),fInFilePrefix.c_str(),i) );
  }

  if(!gTree) return 1;

  TVector3 nudir(0.,0.,1.);
  TH3D * angles[4];
  if( rotate ) {
    if( dump_dir ) {
      // NuMI beam is 58 mrad down and ~405 mrad off of the ICARUS z axis in the +x direction
      // If you make a unit vector out of that you get this, the NuMI beam axis in the ICARUS coordinate system
      TVector3 NuMI_beam_dir( 0.39446818847294496, -0.05320733704983909, 0.9173678801696394 );
      TVector3 icarus_position_numi(448.2, 7962.6, 79519.1); // cm, NuMI beam coordinates i.e. z axis is beam, origin is target (DocDB 19294) 
      TVector3 dump_position_numi(0., 0., 71500.0); // cm, dump is 715m downstream of target
      TVector3 kdar_dir = (icarus_position_numi-dump_position_numi); // in the NuMI beam coordinates
      kdar_dir.RotateUz(NuMI_beam_dir);
      nudir = kdar_dir.Unit();
    } else { // use an angle file
      TFile * tf = new TFile( fAngleFile.c_str() );
      angles[0] = (TH3D*) tf->Get("angles_numu");
      angles[1] = (TH3D*) tf->Get("angles_nue");
      angles[2] = (TH3D*) tf->Get("angles_anumu");
      angles[3] = (TH3D*) tf->Get("angles_anue");
    }
  }

  NtpMCEventRecord * mcrec = 0;
  gTree->SetBranchAddress("gmcrec", &mcrec);

  int nev = gTree->GetEntries();

  // Loop over all events
  for( int i = 0; i < nev; i++ ) {

    if( i % 10000 == 0 ) std::cout << "Event " << i << " of " << nev << "..." << std::endl;

    gTree->GetEntry(i);

    // Fill event accouting
    evtNo = i;
    globalEvtNo = i;
    spillNo = i;

    // get the GENIE event
    EventRecord &  event = *(mcrec->event);

    double xsec = event.XSec() / (1.0 * units::cm2);
    //std::cout << "unit conv " << (1.*units::cm2) << " xsec " << xsec << std::endl;

    Interaction *in = event.Summary();

    TLorentzVector lep = in->Kine().FSLeptonP4();
    TLorentzVector nu = *(in->InitState().GetProbeP4(kRfLab));
    TLorentzVector q = nu-lep;

    // interaction-level variables
    x = -q.Mag2()/(2*0.939*q.E());
    y = q.E() / nu.E();
    Q2 = -q.Mag2();
    W = sqrt(0.939*0.939 + 2.*q.E()*0.939 + q.Mag2());
    scat = in->ProcInfo().ScatteringTypeId();
    tgt = in->InitState().Tgt().Pdg();
    //struck = in->InitState().Tgt().HitNucPdg();
    Enu = in->InitState().ProbeE(kRfLab);
    nuPdg = in->InitState().ProbePdg();

    if( rotate && !dump_dir ) { // event-by-event rotation
      int nuidx = -1;
      if     ( nuPdg ==  14 ) nuidx = 0;
      else if( nuPdg ==  12 ) nuidx = 1;
      else if( nuPdg == -14 ) nuidx = 2;
      else if( nuPdg == -12 ) nuidx = 3;
      else continue; // this should never happen

      int ebin = angles[nuidx]->GetXaxis()->FindBin(Enu);
      double cosx, cosy;
      angles[nuidx]->GetXaxis()->SetRange(ebin,ebin);
      TH2D * proj = (TH2D*)angles[nuidx]->Project3D("zy");
      proj->GetRandom2(cosx,cosy);
      nuPx = Enu * cosx;
      nuPy = Enu * cosy;
      nuPz = sqrt(Enu*Enu - nuPx*nuPx - nuPy*nuPy);
    } else {
      nuPx = (Enu*nudir).x();
      nuPy = (Enu*nudir).y();
      nuPz = (Enu*nudir).z();
    }

    // Loop over all particles in this event, fill particle variables
    GHepParticle * p = 0;
    TIter event_iter(&event);

    //char * evtstr[100];
    //int ntotfsp = 0;

    nfsp = 0;
    while((p=dynamic_cast<GHepParticle *>(event_iter.Next()))) {

      //evtstr[ntotfsp] = Form("%d %d E = %1.3f", p->Status(), p->Pdg(), p->P4()->E() );
      //ntotfsp += 1;

      if (p->Status() == kIStStableFinalState ) {
        pdg[nfsp] = p->Pdg();
        TVector3 mom = p->P4()->Vect();
        mom.RotateUz(nudir);
        px[nfsp] = mom.x();
        py[nfsp] = mom.y();
        pz[nfsp] = mom.z();
        E[nfsp] = p->P4()->E();
        ++nfsp;
      }
    }// end loop over particles	
    //if( nfsp == 1 ) {
      //for( int a = 0; a < ntotfsp; ++a ) std::cout << evtstr[a] << std::endl;
      //std::cout << std::endl;
    //}

    // clear current mc event record
    mcrec->Clear();

    tree->Fill();

  }//end loop over events

  // close input GHEP event file
  //inF->Close();

  outFile->cd();
  tree->Write();

  outFile->Close();

  return 0;
}


//___________________________________________________________________
