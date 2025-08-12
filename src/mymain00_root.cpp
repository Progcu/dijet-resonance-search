// mymain00_root.cpp
// Date: 11 Aug 2025
// Author: Mete AYDIN

// Keywords:
//    ROOT
//    Invariant Mass
//    Data Analysis
//    Graphing

// Description:
// This ROOT macro reads jet data from proton-proton collision simulation output,
// compares invariant masses calculated using FastJet and formula-based methods,
// and visualizes the results through graphs. The macro compares invariant masses
// for Jet-1, Jet-2, and combined dijet systems to assess consistency between
// different calculation approaches.

// --- Required Libraries ---

// --- Standard C++ libraries for I/O, string manipulation, math, and data structures ---
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <cctype>
#include <iomanip>
#include <vector>

// --- ROOT libraries for plotting graphs, canvases, legends, and axis customization ---
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>

using namespace std;

// --- Returns a small horizontal offset (jitter) to separate overlapping points on the X-axis ---
double jitter(int i) {
    static double offsets[] = {-0.05, 0.05};
    return offsets[i % 2];
}

template<typename T> // Template allows function to work with any data type

// --- Generic function to find min and max in any numeric vector ---
// --- Works in-place by modifying minVal and maxVal through reference parameters ---
void findMinMax(const vector<T> &v, double &minVal, double &maxVal) {

    if (v.empty()) return;
    minVal = v[0];
    maxVal = v[0];

    for (auto &val : v) {
        if (val < minVal) minVal = val;
        if (val > maxVal) maxVal = val;
    }
}

// --- Returns the square root of x, but safely handles invalid or non-positive inputs ---
static double safe_sqrt(double x) { // 'static' limits this function's visibility to the current source file (internal linkage)
    if (!isfinite(x) || x <= 0.0) return 0.0;
    return sqrt(x);
}

// --- Try to open the input file and verify it can be read ---
void compareInvariantMass() {
    const string filename = "/home/calcifer/Desktop/dijet-resonance-search/data/mymain00.txt";
    ifstream in(filename);
    if (!in.is_open()) {
        cerr << "Failed to open file: " << filename << endl;
        return; 
    }

    cout << fixed << setprecision(3); // Set output format, initialize variables for jet/event data, and prepare containers for mass calculations

    string line;
    int eventNumber = -1;
    int currentJet = -1; // Index of the jet being processed (0 or 1)
    bool haveJet[2] = {false, false};
    bool printed = false;
    
    // Arrays to store kinematic properties for up to two jets
    double fastMass[2] = {0.,0.};
    double px[2] = {0.,0.}, py[2] = {0.,0.}, pz[2] = {0.,0.}, E[2] = {0.,0.};
    double eta[2] = {0.,0.}, phi[2] = {0.,0.};
    
    // Vectors to store calculated masses for later comparison and analysis
    vector<double> fastJet1Masses, fourVectorJet1Masses; 
    vector<double> fastJet2Masses, fourVectorJet2Masses; 
    vector<double> fourVectorDijetMasses, colliderExpDijetMasses; 
    
    // --- Lambda function to print event information and calculate invariant masses ---
    auto printEvent = [&](){

        if (eventNumber < 0) return;

        cout << "Event-" << eventNumber << endl;
        cout << "FastJet Jet-1 invariant mass: " << fastMass[0] << endl;
        cout << "FastJet Jet-2 invariant mass: " << fastMass[1] << endl;

        // --- Calculate invariant mass using energy and momentum components for each jet ---
        double m1sq = haveJet[0] ? (E[0]*E[0] - (px[0]*px[0] + py[0]*py[0] + pz[0]*pz[0])) : -1.0;
        double m2sq = haveJet[1] ? (E[1]*E[1] - (px[1]*px[1] + py[1]*py[1] + pz[1]*pz[1])) : -1.0;
        double fourvectorJet1 = safe_sqrt(m1sq); 
        double fourvectorJet2 = safe_sqrt(m2sq); 

        cout << "Four vector Jet-1 invariant mass: " << fourvectorJet1 << endl;
        cout << "Four vector Jet-2 invariant mass: " << fourvectorJet2 << endl;

        // --- Calculate dijet invariant mass from combined energy and momentum vectors ---
        double dijet = 0.0;
        if (haveJet[0] && haveJet[1]) { 
            
            double Et = E[0] + E[1];
            double pxt = px[0] + px[1];
            double pyt = py[0] + py[1];
            double pzt = pz[0] + pz[1];
            
            dijet = safe_sqrt(Et*Et - (pxt*pxt + pyt*pyt + pzt*pzt));
        }

        cout << "Four vector Dijet invariant mass: " << dijet << endl;

        // --- Calculate dijet invariant mass using collider experiment formula based on pt, eta, and phi ---
        double collExpMass = 0.0;
        if (haveJet[0] && haveJet[1]) {

            double pt1 = sqrt(px[0]*px[0] + py[0]*py[0]);
            double pt2 = sqrt(px[1]*px[1] + py[1]*py[1]);
            double deltaEta = eta[0] - eta[1];
            double deltaPhi = phi[0] - phi[1];

            // --- Normalize deltaPhi to range [-pi, pi] ---    
            while (deltaPhi > M_PI) deltaPhi -= 2*M_PI;
            while (deltaPhi < -M_PI) deltaPhi += 2*M_PI;

            double collExpmassSq = 2 * pt1 * pt2 * ( cosh(deltaEta) - cos(deltaPhi) ); 
            
            collExpMass = safe_sqrt(collExpmassSq);
        }
        cout << "Collider experiment Dijet invariant mass: " << collExpMass << endl;

        cout << "---------------------------------------------------------------" << endl;

        // --- Store values for plotting ---
        fastJet1Masses.push_back(fastMass[0]);
        fourVectorJet1Masses.push_back(fourvectorJet1);
        fastJet2Masses.push_back(fastMass[1]);
        fourVectorJet2Masses.push_back(fourvectorJet2);
        fourVectorDijetMasses.push_back(dijet);
        colliderExpDijetMasses.push_back(collExpMass);

        printed = true;
    };

    // --- Lambda function to reset all event-related variables before processing a new event ---    
    auto resetEvent = [&](){
        currentJet = -1;
        haveJet[0] = haveJet[1] = false;
        printed = false;
        fastMass[0] = fastMass[1] = 0.0;
        px[0]=px[1]=py[0]=py[1]=pz[0]=pz[1]=E[0]=E[1]=0.0;
        eta[0]=eta[1]=phi[0]=phi[1]=0.0;
    };
    
    // --- Read the file line by line until EOF ---
    while (getline(in, line)) {
        size_t posEv = line.find("Event ");
        
        // --- Check if the current line contains the keyword "Event " ---
        if (posEv != string::npos) {
            size_t i = posEv + 6; // Find the position after "Event " and skip any whitespace
            
            while (i < line.size() && isspace((unsigned char)line[i])) ++i; 
            
            // --- Check if the next character is a minus sign or digit (start of event number) ---
            if (i < line.size() && (line[i]=='-' || isdigit((unsigned char)line[i]))) {
                bool neg = false;
                
                // --- If negative sign, mark event number as negative and move forward ---
                if (line[i] == '-') { neg = true; ++i; }
                int val = 0;
                bool any = false;
                
                // --- Parse consecutive digits to extract the integer event number ---
                while (i < line.size() && isdigit((unsigned char)line[i])) {
                    any = true; val = val*10 + (line[i]-'0'); ++i;
                }
                
                // --- If at least one digit was found, handle previous event printing and reset ---
                if (any) {
                    
                    if (eventNumber >= 0 && haveJet[0] && haveJet[1] && !printed) printEvent();
                    eventNumber = neg ? -val : val;
                    resetEvent();
                    continue;
                }
            }
        }

        // --- Check if the current line contains the keyword "Jet" ---
        size_t posJet = line.find("Jet");
        if (posJet != string::npos) {
            size_t j = posJet + 3; // Move index to the position just after "Jet"
            
            // --- Skip any characters until a digit is found (jet number) ---
            while (j < line.size() && !isdigit((unsigned char)line[j])) ++j;
            
            // --- If digit found, extract it as jet index ---
            if (j < line.size() && isdigit((unsigned char)line[j])) {
                int idx = (line[j] - '0') - 1;
                
                // --- If index is 0 or 1 (valid for two jets), set currentJet accordingly and continue ---
                if (idx >= 0 && idx <= 1) { currentJet = idx; continue; }
            }
        }
        
        // --- Search for the substring "mass = " in the current line ---
        size_t posMass = line.find("mass = ");
        if (posMass != string::npos && currentJet >= 0) {
            double m = 0.0;
            sscanf(line.c_str() + posMass, "mass = %lf", &m); // Parse the mass value from the line starting at posMass
            fastMass[currentJet] = m; // Assign the extracted mass value to the corresponding jet's fastMass
            continue;
        }

        // --- Search for "px = " in the line and check if a valid jet index is set ---
        size_t posPx = line.find("px = ");
        if (posPx != string::npos && currentJet >= 0) {
            double a,b,c,d; // Variables to hold px, py, pz and E values
            int read = sscanf(line.c_str() + posPx, "px = %lf, py = %lf, pz = %lf, E = %lf", &a, &b, &c, &d); // Parse momentum components and energy from the line
            
            // --- If all four values were successfully read ---
            if (read == 4) {
                px[currentJet] = a; py[currentJet] = b; pz[currentJet] = c; E[currentJet] = d; // Store the momentum components and energy for the current jet
                haveJet[currentJet] = true;
            }
            continue;
        }
        
        // --- Search for "mass = " in the line and check if a valid jet index is set ---
        size_t posEta = line.find("eta = ");
        if (posEta != string::npos && currentJet >= 0) {
            double e = 0.0;
            sscanf(line.c_str() + posEta, "eta = %lf", &e); // Parse eta value from the line
            eta[currentJet] = e;
            continue;
        }

        // --- Search for "phi = " and verify current jet index ---
        size_t posPhi = line.find("phi = ");
        if (posPhi != string::npos && currentJet >= 0) {
            double p = 0.0;
            sscanf(line.c_str() + posPhi, "phi = %lf", &p); // Parse phi value from the line
            phi[currentJet] = p;
            
            if (haveJet[0] && haveJet[1] && !printed) printEvent(); // If both jets have data and event not printed yet, output event info
            continue;
        }
    }

    // --- Ensure last event data is printed if not already done ---
    if (eventNumber >= 0 && haveJet[0] && haveJet[1] && !printed) printEvent();

    in.close(); // Close the input file after reading all lines

    int n = fastJet1Masses.size(); // Number of events for plotting

    // --- Find min and max of four vector dijet masses and collider experiment dijet masses for y-axis limits ---
    double yMin, yMax;
    findMinMax(fourVectorDijetMasses, yMin, yMax);
    double yMin2, yMax2;
    findMinMax(colliderExpDijetMasses, yMin2, yMax2);
    
    // --- Update yMin and yMax to cover full data range ---
    if (yMax2 > yMax) yMax = yMax2;
    
    // --- Add 10% margin for better plot appearance ---
    double margin = (yMax - yMin)*0.1;
    yMin -= margin;
    yMax += margin;

    TCanvas *c1 = new TCanvas("c1","Jet 1 Mass Comparison",800,600); // Create a new canvas for plotting Jet 1 mass comparison
    TGraph *gr_fastJet1 = new TGraph(n); // Create graph for FastJet Jet 1 masses with 'n' points
    TGraph *gr_fourvectorJet1 = new TGraph(n); // Create graph for four vector-calculated Jet 1 masses with 'n' points 

    // --- Loop over all events to set graph points with jitter to avoid overlap ---
    for (int i=0; i<n; ++i) {
        gr_fastJet1->SetPoint(i, i + jitter(0), fastJet1Masses[i]);
        gr_fourvectorJet1->SetPoint(i, i + jitter(1), fourVectorJet1Masses[i]);
    }

    // --- Customize FastJet (Jet-1) graph markers and color ---
    gr_fastJet1->SetMarkerStyle(20);
    gr_fastJet1->SetMarkerColor(kRed);
    gr_fastJet1->SetTitle("Jet 1 Mass Comparison");
    gr_fastJet1->GetXaxis()->SetTitle("Event Index");
    gr_fastJet1->GetYaxis()->SetTitle("Invariant Mass (GeV)");
    gr_fastJet1->Draw("AP");
    
    // --- Customize four vector (Jet-1) graph markers and color ---
    gr_fourvectorJet1->SetMarkerStyle(21);
    gr_fourvectorJet1->SetMarkerColor(kBlue);
    gr_fourvectorJet1->Draw("P SAME");

    // --- Create and draw legend explaining the graph points ---
    TLegend *leg1 = new TLegend(0.6,0.7,0.9,0.9);
    leg1->AddEntry(gr_fastJet1,"FastJet Jet 1","p");
    leg1->AddEntry(gr_fourvectorJet1,"Four vector Jet 1","p");
    leg1->Draw();

    c1->Update();  // Update the canvas to display the plot

    TCanvas *c2 = new TCanvas("c2","Jet 2 Mass Comparison",800,600); // Create a new canvas for plotting Jet 2 mass comparison
    TGraph *gr_fastJet2 = new TGraph(n); // Create graph for FastJet Jet 2 masses with 'n' points
    TGraph *gr_fourvectorJet2 = new TGraph(n); // Create graph for four vector-calculated Jet 2 masses with 'n' points 

    // --- Loop over all events to set graph points with jitter to avoid overlap --- 
    for (int i=0; i<n; ++i) {
        gr_fastJet2->SetPoint(i, i + jitter(0), fastJet2Masses[i]);
        gr_fourvectorJet2->SetPoint(i, i + jitter(1), fourVectorJet2Masses[i]);
    }

    // --- Customize FastJet (Jet-2) graph markers and color ---
    gr_fastJet2->SetMarkerStyle(20);
    gr_fastJet2->SetMarkerColor(kRed);  
    gr_fastJet2->SetTitle("Jet 2 Mass Comparison");
    gr_fastJet2->GetXaxis()->SetTitle("Event Index");
    gr_fastJet2->GetYaxis()->SetTitle("Invariant Mass (GeV)");
    gr_fastJet2->Draw("AP");

    // --- Customize four vector (Jet-2) graph markers and color ---
    gr_fourvectorJet2->SetMarkerStyle(21);
    gr_fourvectorJet2->SetMarkerColor(kBlue);
    gr_fourvectorJet2->Draw("P SAME");

    // --- Create and draw legend explaining the graph points ---
    TLegend *leg2 = new TLegend(0.6,0.7,0.9,0.9);
    leg2->AddEntry(gr_fastJet2,"FastJet Jet 2","p");
    leg2->AddEntry(gr_fourvectorJet2,"Formula Jet 2","p");
    leg2->Draw();

    c2->Update(); // Update the canvas to display the plot

    TCanvas *c3 = new TCanvas("c3","Dijet Mass Comparison",800,600); // Create canvas for dijet mass comparison
    TGraph *gr_fourvectorDijet = new TGraph(n); // Create graph for four vector-calculated dijet masses with 'n' points 
    TGraph *gr_collExpDijet = new TGraph(n); // Create graph for collider experiment-calculated dijet masses with 'n' points

    // --- Fill graphs with dijet mass data for each event ---
    for (int i=0; i<n; ++i) {
        gr_fourvectorDijet->SetPoint(i, i, fourVectorDijetMasses[i]);
        gr_collExpDijet->SetPoint(i, i, colliderExpDijetMasses[i]);
    }

    // --- Customize four vector-calculated dijet graph markers and color ---
    gr_fourvectorDijet->SetMarkerStyle(20);
    gr_fourvectorDijet->SetMarkerColor(kRed);
    gr_fourvectorDijet->SetTitle("Dijet Mass Comparison");
    gr_fourvectorDijet->GetXaxis()->SetTitle("Event Index");
    gr_fourvectorDijet->GetYaxis()->SetTitle("Invariant Mass (GeV)");

    // --- Set y-axis limits based on previously computed min and max values ---
    gr_fourvectorDijet->GetYaxis()->SetLimits(yMin, yMax);
    gr_fourvectorDijet->GetYaxis()->SetRangeUser(yMin, yMax);
    
    // --- Draw the formula dijet graph with axes and points ---
    gr_fourvectorDijet->Draw("AP"); 

    // --- Customize collider experiment-calculated dijet graph markers and color ---
    gr_collExpDijet->SetMarkerStyle(21);
    gr_collExpDijet->SetMarkerColor(kBlue);
    gr_collExpDijet->Draw("P SAME");

    // --- Create and draw legend explaining the graph points ---
    TLegend *leg3 = new TLegend(0.6,0.7,0.9,0.9);
    leg3->AddEntry(gr_fourvectorDijet,"Four Vector Dijet Mass","p");
    leg3->AddEntry(gr_collExpDijet,"Coll. Exp.   Dijet Mass","p");
    leg3->Draw();

    c3->Update(); // Update the canvas to display the plot
}
