// mymain00.cpp
// Date: 11 Aug 2025
// Author: Mete AYDIN

// Keywords:
//    Pythia8
//    FastJet
//    Jet Clustering

// Description:
// This program simulates proton-proton collisions at 13 TeV center-of-mass energy
// using the Pythia8 event generator. It performs jet clustering with the FastJet
// package using the anti-kt algorithm and outputs kinematic properties of the two
// highest transverse momentum jets for each event. The simulation includes hard QCD
// processes with a minimum transverse momentum threshold.
//
// The program is intended as a basic demonstration of combining Pythia8 and FastJet
// for jet analysis in high-energy physics.


// --- Required Libraries ---
#include "Pythia8/Pythia.h" // Main Pythia8 library for event generation and particle simulation
#include "fastjet/ClusterSequence.hh" // FastJet library for jet clustering algorithms
#include <iostream>
#include <vector>
#include <cmath>

using namespace Pythia8;
using namespace fastjet;
using namespace std;

// --- Main function --- 
// Initialize Pythia with proton-proton collisions at 13 TeV,
// set QCD processes and event generation parameters,
// then start event loop.
int main() {

    Pythia pythia; // Create and configure Pythia event generator

    // Settings
    pythia.readString("Beams:idA  = 2212"); // The PDG id code for the first incoming particle (id: 2212 'proton')
    pythia.readString("Beams:idB  = 2212"); // The PDG id code for the second incoming particle (id: 2212 'proton')
    pythia.readString("Beams:eCM = 13000"); // Collision CM (13 TeV)
    pythia.readString("HardQCD:all = on"); // Common switch for the group of all hard QCD 2→2 processes, as listed separately in the following
    pythia.readString("PhaseSpace:pTHatMin = 30"); // The minimum invariant pT (transverse momentum --> min 30 GeV)
    pythia.readString("Next:numberShowInfo = 0"); // The number of events to list the Info information 
    pythia.readString("Next:numberShowProcess = 1"); // The number of events to list the process record
    pythia.readString("Next:numberShowEvent = 1"); // The number of events to list the event record
    pythia.readString("PartonLevel:MPI = on"); // Master switch for multiparton interactions
    pythia.readString("HadronLevel:all = on"); // If off then stop the generation after the hard process and parton-level activity has been generated, but before the hadron-level steps
    
    pythia.init(); // Initialize the generator with the above settings

     // --- Number of events to generate ---
    int nEvents = 20;  

    // --- Define jet clustering algorithm with radius parameter R ---
    double R = 0.8;
    JetDefinition jetDef(antikt_algorithm, R);

    // --- The event loop ---
    for (int iEvent = 0; iEvent < nEvents; ++iEvent){
    
        if(!pythia.next()) continue;

        vector<PseudoJet> particles; // Vector to store final-state particles for jet clustering
        
        // --- Loop over all particles in the current event ---
        for (int i = 0; i < pythia.event.size(); ++i){ 

            if (pythia.event[i].isFinal() && pythia.event[i].isVisible()){ // Select only final, visible particles (stable and detectable)

                particles.push_back(PseudoJet( // Add particle four-momentum to the vector as a PseudoJet object
 
                    pythia.event[i].px(),
                    pythia.event[i].py(),
                    pythia.event[i].pz(),
                    pythia.event[i].e()
                
                ));

            }

        }
        
        ClusterSequence cs(particles, jetDef); // Perform jet clustering on the collected particles using the specified jet definition
        vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets()); // Extract inclusive jets sorted by transverse momentum (pT) in descending order

        cout << "Event " << iEvent << ": " << jets.size() << " jets found" << endl;
        
        // --- If at least two jets are found, print properties of the leading two jets ---
        if (jets.size() >= 2){

            for (int j = 0; j < 2; ++j){

                double px = jets[j].px();
                double py = jets[j].py();
                double pz = jets[j].pz();
                double E  = jets[j].e();
                double pt  = jets[j].pt();
                double phi = jets[j].phi();
                double eta = jets[j].eta();
                double mass = jets[j].m();
                double p   = sqrt(px*px + py*py + pz*pz);

                cout << " Jet " << j+1 << ":" << endl;
                cout << "  mass = " << mass << endl;
                cout << "  px = " << px << ", py = " << py << ", pz = " << pz << ", E = " << E << endl;
                cout << "  pT = " << pt << ", eta = " << eta << ", phi = " << phi << endl;

            }

        }

    } // End of event loop

    // --- Print final statistics from Pythia run ---
    pythia.stat();

    return 0;

}
