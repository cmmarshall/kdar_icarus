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

  int evtNo, spillNo, globalEvtNo, nuPdg, nfsp, tgt, scat, struck;
  double Enu, x, y, Q2, W;

  int pdg[100];
  double px[100], py[100], pz[100], E[100];

  tree->Branch( "evt", &evtNo, "evt/I" );
  tree->Branch( "globalEvt", &globalEvtNo, "globalEvt/I" );
  tree->Branch( "spill", &spillNo, "spill/I" );
  tree->Branch( "Enu", &Enu, "Enu/D" );
  tree->Branch( "nuPdg", &nuPdg, "nuPdg/I" );
  tree->Branch( "x", &x, "x/D" );
  tree->Branch( "y", &y, "y/D" );
  tree->Branch( "Q2", &Q2, "Q2/D" );
  tree->Branch( "W", &W, "W/D" );
  tree->Branch( "tgt", &tgt, "tgt/I" );
  tree->Branch( "struck", &struck, "struck/I" );
  tree->Branch( "scat", &scat, "scat/I" );
  tree->Branch( "nfsp", &nfsp, "nfsp/I" );
  tree->Branch( "pdg", pdg, "pdg[nfsp]/I" );
  tree->Branch( "px", px, "px[nfsp]/D" );
  tree->Branch( "py", py, "py[nfsp]/D" );
  tree->Branch( "pz", pz, "pz[nfsp]/D" );
  tree->Branch( "E", E, "E[nfsp]/D" );

  //-- open the ROOT file and get the TTree & its header
  TChain * gTree = new TChain( "gtree", "gtree" );

  for( int i = first_run; i <= last_run; ++i ) {
    gTree->Add( Form("%s/%s.%d.ghep.root",fInFileDir.c_str(),fInFilePrefix.c_str(),i) );
  }

  if(!gTree) return 1;

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
    struck = in->InitState().Tgt().HitNucPdg();
    Enu = in->InitState().ProbeE(kRfLab);
    nuPdg = in->InitState().ProbePdg();

    // Loop over all particles in this event, fill particle variables
    GHepParticle * p = 0;
    TIter event_iter(&event);

    nfsp = 0;
    while((p=dynamic_cast<GHepParticle *>(event_iter.Next()))) {

      if (p->Status() == kIStStableFinalState ) {
        pdg[nfsp] = p->Pdg();
        TLorentzVector mom = *(p->P4());
        px[nfsp] = mom.X();
        py[nfsp] = mom.Y();
        pz[nfsp] = mom.Z();
        E[nfsp] = mom.E();
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
