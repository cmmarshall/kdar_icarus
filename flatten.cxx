#include <string>
#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TIterator.h>
#include <TLorentzVector.h>

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
int first_run, last_run;

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
  //tree->Branch( "struck", &struck, "struck/I" );
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

  // NuMI beam is 58 mrad down and ~405 mrad off of the ICARUS z axis in the +x direction
  // If you make a unit vector out of that you get this, the NuMI beam axis in the ICARUS coordinate system
  TVector3 NuMI_beam_dir( 0.39446818847294496, -0.05320733704983909, 0.9173678801696394 );
  TVector3 icarus_position_numi(448.2, 7962.6, 79519.1); // cm, NuMI beam coordinates i.e. z axis is beam, origin is target (DocDB 19294) 
  TVector3 dump_position_numi(0., 0., 71500.0); // cm, dump is 715m downstream of target
  TVector3 kdar_dir = (icarus_position_numi-dump_position_numi); // in the NuMI beam coordinates
  kdar_dir.RotateUz(NuMI_beam_dir);
  TVector3 nudir = kdar_dir.Unit();

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

    nuPx = (Enu*nudir).x();
    nuPy = (Enu*nudir).y();
    nuPz = (Enu*nudir).z();

    // Loop over all particles in this event, fill particle variables
    GHepParticle * p = 0;
    TIter event_iter(&event);

    nfsp = 0;
    while((p=dynamic_cast<GHepParticle *>(event_iter.Next()))) {

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
