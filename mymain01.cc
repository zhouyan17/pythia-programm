// Headers and Namespaces.
#include <iostream>
#include "Pythia8/Pythia.h" // Include Pythia headers.
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "TH1.h"
#include "TVirtualPad.h"
#include "TApplication.h"
#include "TFile.h"
using namespace Pythia8;
int main(int argc, char* argv[]) {
    
// Set up generation.
Pythia pythia;
pythia.readFile(argv[1]);
pythia.init();
    
// Create the ROOT application environment.
TApplication theApp("hist", &argc, argv);
TFile* outFile = new TFile(argv[2], "RECREATE");

fastjet::JetDefinition jetDef( fastjet::antikt_algorithm, 0.4, fastjet::E_scheme, fastjet::Best);
std::vector <fastjet::PseudoJet> fjInputs;
    
// Histograms of parton level
TH1F *HMassR = new TH1F("HMassR","Mass of new resonances", 100, 0., 1000);
TH1F *HptR = new TH1F("HptR","pT of new resonances", 100, 0., 1000);
TH1F *HptD1 = new TH1F("HptD1","pT of daughter 1", 100, 0., 1000);
TH1F *HptD2 = new TH1F("HptD2","pT of daughter 2", 100, 0., 1000);
TH1F *HetaR = new TH1F("HetaR","eta of new resonances", 100, 0., 2.);
TH1F *HetaD1 = new TH1F("HetaD1","eta of daughter 1", 100, 0., 2.);
TH1F *HetaD2 = new TH1F("HetaD2","eta of daughter 2", 100, 0., 2.);

// Histograms of particle level
TH1F *HMassC = new TH1F("HMassC","Mass of dijet candidates", 100, 0., 1000);
TH1F *HptC = new TH1F("HptC","pT of dijet candidates", 100, 0., 1000);
TH1F *HptJ1 = new TH1F("HptJ1","pT of leading jets", 100, 0, 1000);
TH1F *HptJ2 = new TH1F("HptJ2","pT of sub leading jets", 100, 0, 1000);
TH1F *HetaC = new TH1F("HetaC","eta of dijet candidates", 100, 0., 2.);
TH1F *HetaJ1 = new TH1F("HetaJ1","eta of leading jets", 100, 0., 2.);
TH1F *HetaJ2 = new TH1F("HetaJ2","eta of sub leading jets", 100, 0., 2.);
TH1F *HYstar = new TH1F("HYstar","Ystar", 100, 0., 0.6);

for (int iEvent = 0; iEvent < 10000; ++iEvent)
{
    pythia.next();
    int iex = 0;
    for (int i = 0; i < pythia.event.size(); ++i)
    {
        if (pythia.event[i].id() == 4000001 || pythia.event[i].id() == 4000002)
        {
            iex = i;
        }
    }
    int iD1 = pythia.event[iex].daughter1();
    int iD2 = pythia.event[iex].daughter2();
    
    HMassR->Fill(pythia.event[iex].m());
    HptR->Fill(pythia.event[iex].p().pT());
    HptD1->Fill(pythia.event[iD1].p().pT());
    HptD2->Fill(pythia.event[iD2].p().pT());
    HetaR->Fill(pythia.event[iex].p().eta());
    HetaD1->Fill(pythia.event[iD1].p().eta());
    HetaD2->Fill(pythia.event[iD2].p().eta());
    
    fjInputs.clear();
    for (int i = 0; i < pythia.event.size(); ++i)
    {
        if (!pythia.event[i].isFinal()) continue;

          // No neutrinos or DM.
        if ( pythia.event[i].idAbs() == 12 || pythia.event[i].idAbs() == 14 || pythia.event[i].idAbs() == 16 || pythia.event[i].idAbs() == 52)
        continue;

          // Only |eta| < 3.6.
          //      if (abs(pythia.event[i].eta()) > 3.6) continue;
          // still don't know how to cut.

          // Store as input to Fastjet.
        fjInputs.push_back( fastjet::PseudoJet( pythia.event[i].px(),
        pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e() ) );
    }
        // Run Fastjet algorithm.
        vector <fastjet::PseudoJet> inclusiveJets, sortedJets;
        fastjet::ClusterSequence clustSeq(fjInputs, jetDef);

        // Extract inclusive jets sorted by pT (note minimum pT of 20.0 GeV).
        inclusiveJets = clustSeq.inclusive_jets(20.0);
        sortedJets    = sorted_by_pt(inclusiveJets);

        // cut the events where jets have too little pT
        //    if (sortedJets[0].pt() <= 420 || sortedJets[1].pt() <= 150)
        //      continue;
    
        if(sortedJets.size() < 2) continue;

        Vec4 pJ1(sortedJets[0].px(), sortedJets[0].py(), sortedJets[0].pz(), sortedJets[0].e());
        Vec4 pJ2(sortedJets[1].px(), sortedJets[1].py(), sortedJets[1].pz(), sortedJets[1].e());
        Vec4 pC = pJ1 + pJ2;

        double Ystar = (sortedJets[0].rap() - sortedJets[1].rap()) / 2;

        // cut the events where jet eta or ystar is too large
        //    if (abs(etaC) >= 2 || abs(Ystar) >= 0.6) continue;

        HMassC->Fill(pC.mCalc());
        HptC->Fill(pC.pT());
        HptJ1->Fill(sortedJets[0].pt());
        HptJ2->Fill(sortedJets[1].pt());
        HetaC->Fill(pC.eta());
        HetaJ1->Fill(sortedJets[0].eta());
        HetaJ2->Fill(sortedJets[1].eta());
        HYstar->Fill(Ystar);
        //cout << "i = " << i << ", id = " << pythia.event[i].id() << endl;
    //cout << pythia.event[iTop].charge();
}
HMassC->Draw();
gPad -> WaitPrimitive();
HMassR->Write();
HptR->Write();
HptD1->Write();
HptD2->Write();
HetaR->Write();
HetaD1->Write();
HetaD2->Write();
HMassC->Write();
HptC->Write();
HptJ1->Write();
HptJ2->Write();
HetaC->Write();
HetaJ1->Write();
HetaJ2->Write();
HYstar->Write();
delete outFile;

pythia.stat();
return 0;
}
// End main program with error-free return.
