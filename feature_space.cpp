#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <numeric>
#include <set>
#include <tuple>
#include <string>
#include <vector>
#include <algorithm>
#include <TFile.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <THStack.h>
#include <TVector.h>
#include <TVectorT.h>
#include <TVectorF.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TBranch.h>
#include <TMath.h>
#include <TF1.h>
#include <TMultiDimFit.h>
#include <TPaveText.h>
#include <TLegend.h>
using namespace std;

// A function to create the x% truncated mean from a vector of ADC counts
float TM(std::vector<float> adcs, float truncation_percent) {
    float size = adcs.size();
    float truncated_size = size * (1 - truncation_percent);

    // Check if the truncated size is already an integer
    if (truncated_size - std::floor(truncated_size) == 0) {
        truncated_size = ceil(size * (1-truncation_percent));
    }
    // Sort the ADC counts in ascending order
    std::sort(adcs.begin(), adcs.end());

    // Cut off the truncate amount hits
    float track_sum = 0;
    for (int i = 0; i < truncated_size; i++) {
        track_sum += adcs.at(i);
    }
    // Calculate truncated mean
    float truncated_mean = track_sum / truncated_size;
    // Then return the truncated mean and its error
    return truncated_mean;
}

// Return a Truncated mean from each xBin in a TH2 object
std::vector<float> fit2DHistogramTM(TH2 *h2) {
    int nbinsx = h2->GetNbinsX();

    //Record the truncate dmean in each bin
    std::vector<float> truncated_means;
    //Loop through the bins
    for (int i = 1; i <= nbinsx; i++) {
        // Turn each xbin into a 1d histogram
        TH1D* projectionY = h2->ProjectionY("", i, i);
        int nbinsy = projectionY->GetNbinsX();
        std::vector<float> temp_bin;
        // Loop through the ybins of the 1d histogram
        for (int j = 1; j <= nbinsy; j++) {
            // Dump the hits in each bin into a temporary vector
            int content = projectionY->GetBinContent(j);
            // Get the value at the center of the current bin
            float bin_val = projectionY->GetXaxis()->GetBinCenter(j);
            // Content returns how many entries are in this bin. So we need to push back the bin value
            // as many times as there are entries in the bin
            for (int k=0; k<content; k++) {
                temp_bin.push_back(bin_val);
            }
        }
        // Special case if the xbin has no hits
        if (temp_bin.size() == 0) {
            truncated_means.push_back(0);
            delete projectionY;
        }
        else {
            // Calculate the truncated mean and uncertainty for that bin
            float truncated_mean = TM(temp_bin, 0.4);
            truncated_means.push_back(truncated_mean);
            delete projectionY;
        }
    }
    return truncated_means;
}

// Return a Truncated mean from each xBin in a TH2 object
std::vector<std::vector<float>> fit2DHistogramLandau(TH2 *h2, const std::string &option="") {
    int nbinsx = h2->GetNbinsX();
    const char* h2_title = h2->GetTitle();

    // Make a figure where all the bin porjections are plotted with the fit on top along with labels.

    // Record the parameters of the Landau fits <MPV, sigma, MPV_err, sigma_err>
    std::vector<std::vector<float>> landau_params;
    //Loop through the bins
    for (int i = 1; i <= nbinsx; i++) {
        // Turn each xbin into a 1d histogram
        TH1D* projectionY = h2->ProjectionY("", i, i);

        // Fit a landau to the xbin histogram
        projectionY->Fit("landau");

        // Retrieve the fit result
        TF1* fitResult = projectionY->GetFunction("landau");\
        float mpv = fitResult->GetParameter(1); // Most probable value (MPV)
        float sigma = fitResult->GetParameter(2); // Sigma (width of the distribution)
        float mpvError = fitResult->GetParError(1); // Error on MPV
        float sigmaError = fitResult->GetParError(2); // Error on Sigma
        std::vector<float> results = {mpv, sigma, mpvError, sigmaError};
        landau_params.push_back(results);

        if (option == "DRAW") { // Create a canvas on which to draw the fit
            TCanvas *canvas = new TCanvas(Form("Momentum bin%d from %s", i, h2_title), Form("Momentum Bin %d Projection and Landau Fit from %s", i, h2_title), 800, 600);
            TPaveText *label = new TPaveText(0.6, 0.7, 0.9, 0.9, "NDC");
            label->AddText(Form("Chi^2: %f", fitResult->GetChisquare()));
            label->AddText(Form("NDF: %d", fitResult->GetNDF()));
            label->AddText(Form("MPV: %f +/- %.2f", mpv, mpvError));
            label->AddText(Form("Scale: %f +/- %.2f", sigma, sigmaError));
            projectionY->SetTitle(Form("Landau fit for Bin %d from %s", i, h2_title));
            projectionY->Draw();
            fitResult->Draw("same");
            label->Draw("same");
            canvas->Write();
            delete canvas;
            delete label;
        }
        delete fitResult;
    }
    delete h2_title;
    return landau_params;
}

// A function to go from trkMdt_stationPhi to sector
int stationPhiToSector(int stationIndex, int stationPhi) {
    if (stationIndex==0||stationIndex==2||stationIndex==4) {
        return (2 * stationPhi - 1);
    }
    else if (stationIndex==1||stationIndex==3||stationIndex==5) {
        return (2 * stationPhi);
    }
    return 0;
}

// A function to convert transverse momentum to total |p|
float ptToP(float pt, float eta) {
    float p = exp(-eta);
    p = 2 * atan(p);
    p = pt / sin(p);
    return p;
}


void featureSpace(std::string indir, const::std::string& outfile) {

    // Set a random seed
    std::random_device rd;

    // Initialize the random number generator
    std::mt19937 gen(rd());
    std::mt19937 gen2(rd());

    // Define a distribution for numbers between 1 and 3200
    std::uniform_int_distribution<int> dis(1, 3200);
    std::uniform_int_distribution<int> dis2(1, 5700);

    // Initalize the directory that was input to this function and retrieve 
        // the file within
    TSystemDirectory dir(indir.c_str(), indir.c_str());
    TList *files = dir.GetListOfFiles();

    // Intiialize the TChain
    TChain chain("BasicTesterTree");

    // Iterate through the files and add them to the chain if they are valid root files
    if (files) {
        TSystemFile *file;
        TString fname;
        TIter next(files);
        while ((file = (TSystemFile*)next())) {
            fname = file->GetName();
            std::cout << "Opened: " << fname << endl;
            if (!file->IsDirectory() && fname.EndsWith(".root")) {
                std::string result = indir+std::string(fname.Data());
                std:: cout << "File path: " << result << endl;
                chain.Add(result.c_str());
                std::cout << "Added to chain" << endl;
            }
        }
    }

    // Disable all branches from being read
    chain.SetBranchStatus("*", false);

    // Re-enable relevant branches
    chain.SetBranchStatus("trkMdt_muonsLink", true);
    chain.SetBranchStatus("trkMdt_Adc", true);
    chain.SetBranchStatus("trkMdt_driftR", true);
    chain.SetBranchStatus("muons_eta", true);
    chain.SetBranchStatus("trkMdt_stationIndex", true);
    chain.SetBranchStatus("trkMdt_stationEta", true);
    chain.SetBranchStatus("trkMdt_stationPhi", true);
    chain.SetBranchStatus("trkMdt_multiLayer", true);
    chain.SetBranchStatus("trkMdt_tubeLayer", true);
    chain.SetBranchStatus("trkMdt_tube", true);
    chain.SetBranchStatus("trkMdt_HitType", true);
    chain.SetBranchStatus("muons_author", true);
    chain.SetBranchStatus("muons_pt", true);
    chain.SetBranchStatus("muons_eta", true);
    chain.SetBranchStatus("muons_phi", true);
    chain.SetBranchStatus("muons_d0", true);
    chain.SetBranchStatus("muons_z0", true);
    chain.SetBranchStatus("muons_q", true);
    chain.SetBranchStatus("MSTracks_phi", true);
    chain.SetBranchStatus("MStrkMdt_MSTracksLink", true);



    // Define variables to hold branch values
    std::vector<float> *eta = nullptr;
    std::vector<int> *ADC_counts = nullptr;
    std::vector<unsigned short> *trk_link = nullptr;
    std::vector<float> *drift_radius = nullptr;
    std::vector<unsigned char> *station_index = nullptr;
    std::vector<signed char> *station_eta = nullptr;
    std::vector<unsigned char> *station_phi = nullptr;
    std::vector<unsigned char> *multi_layer = nullptr;
    std::vector<unsigned char> *tube_layer = nullptr;
    std::vector<unsigned char> *tubes = nullptr;
    std::vector<short> *hit_type = nullptr;
    std::vector<unsigned short> *authors = nullptr;
    std::vector<float> *pts = nullptr;
    std::vector<float> *etas = nullptr;
    std::vector<float> *phis = nullptr;
    std::vector<float> *d0s = nullptr;
    std::vector<float> *z0s = nullptr;
    std::vector<int> *charges = nullptr;
    std::vector<float> *MS_phis = nullptr;
    std::vector<unsigned char> *MS_trk_link = nullptr;

    // Set the address to each branch in the chain
    chain.SetBranchAddress("muons_eta", &eta);
    chain.SetBranchAddress("trkMdt_Adc", &ADC_counts);
    chain.SetBranchAddress("trkMdt_muonsLink", &trk_link);
    chain.SetBranchAddress("trkMdt_driftR", &drift_radius);
    chain.SetBranchAddress("trkMdt_stationIndex", &station_index);
    chain.SetBranchAddress("trkMdt_stationEta", &station_eta);
    chain.SetBranchAddress("trkMdt_stationPhi", &station_phi);
    chain.SetBranchAddress("trkMdt_multiLayer", &multi_layer);
    chain.SetBranchAddress("trkMdt_tubeLayer", &tube_layer);
    chain.SetBranchAddress("trkMdt_tube", &tubes);
    chain.SetBranchAddress("trkMdt_HitType", &hit_type);
    chain.SetBranchAddress("muons_author", &authors);
    chain.SetBranchAddress("muons_pt", &pts);
    chain.SetBranchAddress("muons_eta", &etas);
    chain.SetBranchAddress("muons_phi", &phis);
    chain.SetBranchAddress("muons_d0", &d0s);
    chain.SetBranchAddress("muons_z0", &z0s);
    chain.SetBranchAddress("muons_q", &charges);
    chain.SetBranchAddress("MSTracks_phi", &MS_phis);
    chain.SetBranchAddress("MStrkMdt_MSTracksLink", &MS_trk_link);


    // A map to get the length of every mdt tube and its position
    // [stationIndex: [sector: [Multilayer: Layer]]]
    // Within a chamber, multiLayer1 is closer to IP and within each multilayer, layers are numbered in ascending order with layer1 closest to IP Fig 3-7
    // Actual MDT propagation length is the active tube length (additional 60mm for electronics)
    
    // An array for the radial distance from IP to chamber center in the XY plane from Fig 4-4 in mm
    float radial_extents[6] = {4926, 4525, 7116, 8070, 9477, 10530}; // Same ordering as trkMdt_stationIndex (BIL, BIS, BML, BMS, BOL, BOS)
    
    // A map to get the radial distance from IP to the midpoint of the MDT tube wire in the XY plane in mm
    std::map<int, std::map<int, std::map<int, float>>> tube_positions; // Structure is [trkMdt_stationIndex: [MultiLayer: [Layer: radial distance in XY plane]]]
    tube_positions[0] = { // BIL Type 1,2,3,4
        {1, {{1, radial_extents[0] - (0.5*170 + 15)}, {2, radial_extents[0] - (0.5*170 + 15 * (1 + pow(3,0.5)))}, {3, radial_extents[0] - (0.5*170 + 15 * (1 + 2*pow(3,0.5)))}, {4, radial_extents[0] - (0.5*170 + 15 * (1 + 3*pow(3,0.5)))}}},
        {2, {{1, radial_extents[0] + (0.5*170 + 15)}, {2, radial_extents[0] + (0.5*170 + 15 * (1 + pow(3,0.5)))}, {3, radial_extents[0] + (0.5*170 + 15 * (1 + 2*pow(3,0.5)))}, {4, radial_extents[0] + (0.5*170 + 15 * (1 + 3*pow(3,0.5)))}}}
    };
    tube_positions[1] = { // BIS Type 2,3 *** BIS Type 1 chambers have a different structure ***
        {1, {{1, radial_extents[1] - (0.5*8 + 15)}, {2, radial_extents[1] - (0.5*8 + 15 * (1 + pow(3,0.5)))}, {3, radial_extents[1] - (0.5*8 + 15 * (1 + 2*pow(3,0.5)))}, {4, radial_extents[1] - (0.5*8 + 15 * (1 + 3*pow(3,0.5)))}}},
        {2, {{1, radial_extents[1] + (0.5*8 + 15)}, {2, radial_extents[1] + (0.5*8 + 15 * (1 + pow(3,0.5)))}, {3, radial_extents[1] + (0.5*8 + 15 * (1 + 2*pow(3,0.5)))}, {4, radial_extents[1] + (0.5*8 + 15 * (1 + 3*pow(3,0.5)))}}}
    };
    tube_positions[2] = { // BML Type 1,2
        {1, {{1, radial_extents[2] - (0.5*317 + 15)}, {2, radial_extents[2] - (0.5*317 + 15*(1 + pow(3,0.5)))}, {3, radial_extents[2] - (0.5*317 + 15*(1 + 2*pow(3,0.5)))}}},
        {2, {{1, radial_extents[2] + (0.5*317 + 15)}, {2, radial_extents[2] + (0.5*317 + 15*(1 + pow(3,0.5)))}, {3, radial_extents[2] + (0.5*317 + 15*(1 + 2*pow(3,0.5)))}}}
    };
    tube_positions[3] = { // BMS Type 1,2,3
        {1, {{1, radial_extents[3] - (0.5*170 + 15)}, {2, radial_extents[3] - (0.5*170 + 15*(1 + pow(3,0.5)))}, {3, radial_extents[3] - (0.5*170 + 15*(1 + 2*pow(3,0.5)))}}},
        {2, {{1, radial_extents[3] + (0.5*170 + 15)}, {2, radial_extents[3] + (0.5*170 + 15*(1 + pow(3,0.5)))}, {3, radial_extents[3] + (0.5*170 + 15*(1 + 2*pow(3,0.5)))}}}
    };
    tube_positions[4] = { // BOL Type 1,2
        {1, {{1, radial_extents[4] - (0.5*317 + 15)}, {2, radial_extents[4] - (0.5*317 + 15*(1 + pow(3,0.5)))}, {3, radial_extents[4] - (0.5*317 + 15*(1 + 2*pow(3,0.5)))}}},
        {2, {{1, radial_extents[4] + (0.5*317 + 15)}, {2, radial_extents[4] + (0.5*317 + 15*(1 + pow(3,0.5)))}, {3, radial_extents[4] + (0.5*317 + 15*(1 + 2*pow(3,0.5)))}}}
    };
    tube_positions[5] = { // BOS Type 1,2,3
        {1, {{1, radial_extents[5] - (0.5*317 + 15)}, {2, radial_extents[5] - (0.5*317 + 15*(1 + pow(3,0.5)))}, {3, radial_extents[5] - (0.5*317 + 15*(1 + 2*pow(3,0.5)))}}},
        {2, {{1, radial_extents[5] + (0.5*317 + 15)}, {2, radial_extents[5] + (0.5*317 + 15*(1 + pow(3,0.5)))}, {3, radial_extents[5] + (0.5*317 + 15*(1 + 2*pow(3,0.5)))}}}
    };

    // NOTE: For large chambers, sector = 2*trkMdt_station - 1. For small chambers, sector = 2*trkMdt_station.

    // A map for the chamber type based on trkMdt_stationIndex and trkMdt_stationEta
    std::map<int, std::map<int, std::map<int, int>>> chamber_types; // Structure is [trkMdt_stationIndex: [trkMdt_stationEta: [sector: chamber type]]]
    chamber_types[0] = { // BIL has less wide chambers in sectors 11 and 15 (trkMdt_stationPhi = 6, 8)
        {1, {{1, 1}, {2, 1}, {3, 1}, {4, 1}, {5, 1}, {6, 2}, {7, 1}, {8, 2}}}, 
        {2, {{1, 3}, {2, 3}, {3, 3}, {4, 3}, {5, 3}, {6, 4}, {7, 3}, {8, 4}}}, 
        {3, {{1, 3}, {2, 3}, {3, 3}, {4, 3}, {5, 3}, {6, 4}, {7, 3}, {8, 4}}}, 
        {4, {{1, 3}, {2, 3}, {3, 3}, {4, 3}, {5, 3}, {6, 4}, {7, 3}, {8, 4}}}, 
        {5, {{1, 3}, {2, 3}, {3, 3}, {4, 3}, {5, 3}, {6, 4}, {7, 3}, {8, 4}}}, 
        {6, {{1, 3}, {2, 3}, {3, 3}, {4, 3}, {5, 3}, {6, 4}, {7, 3}, {8, 4}}}
    };
    chamber_types[1] = {
        {1, {{1, 2}, {2, 2}, {3, 2}, {4, 2}, {5, 2}, {6, 2}, {7, 2}, {8, 2}}},
        {2, {{1, 2}, {2, 2}, {3, 2}, {4, 2}, {5, 2}, {6, 2}, {7, 2}, {8, 2}}},
        {3, {{1, 2}, {2, 2}, {3, 2}, {4, 2}, {5, 2}, {6, 2}, {7, 2}, {8, 2}}},
        {4, {{1, 2}, {2, 2}, {3, 2}, {4, 2}, {5, 2}, {6, 2}, {7, 2}, {8, 2}}},
        {5, {{1, 2}, {2, 2}, {3, 2}, {4, 2}, {5, 2}, {6, 2}, {7, 2}, {8, 2}}},
        {6, {{1, 2}, {2, 2}, {3, 2}, {4, 2}, {5, 2}, {6, 2}, {7, 2}, {8, 2}}},
        {7, {{1, 3}, {2, 3}, {3, 3}, {4, 3}, {5, 3}, {6, 3}, {7, 3}, {8, 3}}},
        {8, {{1, 1}, {2, 1}, {3, 1}, {4, 1}, {5, 1}, {6, 1}, {7, 1}, {8, 1}}}

    };
    chamber_types[2] = {
        {1, {{1, 1}, {2, 1}, {3, 1}, {4, 1}, {5, 1}, {6, 1}, {7, 1}, {8, 1}}},
        {2, {{1, 1}, {2, 1}, {3, 1}, {4, 1}, {5, 1}, {6, 1}, {7, 1}, {8, 1}}},
        {3, {{1, 1}, {2, 1}, {3, 1}, {4, 1}, {5, 1}, {6, 1}, {7, 1}, {8, 1}}},
        {4, {{1, 2}, {2, 2}, {3, 2}, {4, 2}, {5, 2}, {6, 2}, {7, 2}, {8, 2}}},
        {5, {{1, 1}, {2, 1}, {3, 1}, {4, 1}, {5, 1}, {6, 1}, {7, 1}, {8, 1}}},
        {6, {{1, 1}, {2, 1}, {3, 1}, {4, 1}, {5, 1}, {6, 1}, {7, 1}, {8, 1}}}
    };
    chamber_types[3] = { // No BMS in sectors 12 and 14 (trkMdt_stationIndex = 6, 7), they are replaced with BMF
        {1, {{1, 3}, {2, 3}, {3, 3}, {4, 3}, {5, 3}, {8, 3}}},
        {2, {{1, 2}, {2, 2}, {3, 2}, {4, 2}, {5, 2}, {8, 2}}},
        {3, {{1, 3}, {2, 3}, {3, 3}, {4, 3}, {5, 3}, {8, 3}}},
        {4, {{1, 2}, {2, 2}, {3, 2}, {4, 2}, {5, 2}, {8, 2}}},
        {5, {{1, 1}, {2, 1}, {3, 1}, {4, 1}, {5, 1}, {8, 1}}},
        {6, {{1, 2}, {2, 2}, {3, 2}, {4, 2}, {5, 2}, {8, 2}}}
    };
    chamber_types[4] = {
        {1, {{1, 1}, {2, 1}, {3, 1}, {4, 1}, {5, 1}, {6, 1}, {7, 1}, {8, 1}}},
        {2, {{1, 2}, {2, 2}, {3, 2}, {4, 2}, {5, 2}, {6, 2}, {7, 2}, {8, 2}}},
        {3, {{1, 2}, {2, 2}, {3, 2}, {4, 2}, {5, 2}, {6, 2}, {7, 2}, {8, 2}}},
        {4, {{1, 2}, {2, 2}, {3, 2}, {4, 2}, {5, 2}, {6, 2}, {7, 2}, {8, 2}}},
        {5, {{1, 2}, {2, 2}, {3, 2}, {4, 2}, {5, 2}, {6, 2}, {7, 2}, {8, 2}}},
        {6, {{1, 1}, {2, 1}, {3, 1}, {4, 1}, {5, 1}, {6, 1}, {7, 1}, {8, 1}}}
    };
    chamber_types[5] = { // NO BOS in sectors 12 and 14 (trkMdt_stationIndex = 6, 7), they are replaced with BOF
        {1, {{1, 3}, {2, 3}, {3, 3}, {4, 3}, {5, 3}, {8, 3}}},
        {2, {{1, 3}, {2, 3}, {3, 3}, {4, 3}, {5, 3}, {8, 3}}},
        {3, {{1, 3}, {2, 3}, {3, 3}, {4, 3}, {5, 3}, {8, 3}}},
        {4, {{1, 3}, {2, 3}, {3, 3}, {4, 3}, {5, 3}, {8, 3}}},
        {5, {{1, 3}, {2, 3}, {3, 3}, {4, 3}, {5, 3}, {8, 3}}},
        {6, {{1, 2}, {2, 2}, {3, 2}, {4, 2}, {5, 2}, {8, 2}}}
    };

    // A map for the active length of MDT tubes 
    std::map<int, std::map<int, int>> tube_lengths; // Strcuture is [trkMdt_stationIndex: [chamber type: MDT tube length]]
    tube_lengths[0] = {{1, 2700}, {2, 1580}, {3, 2700}, {4, 1580}};
    tube_lengths[1] = {{1, 1000}, {2, 1700}, {3, 1700}};
    tube_lengths[2] = {{1, 3580}, {2, 3580}};
    tube_lengths[3] = {{1, 3100}, {2, 3100}, {3, 3100}};
    tube_lengths[4] = {{1, 4990}, {2, 4990}};
    tube_lengths[5] = {{1, 3800}, {2, 3800}, {3, 3800}};

    // A nested vector to record the number of hits in each chamber as <station index <station Phi <staion etas>>>
    std::vector<std::vector<std::vector<int>>> chambers = vector<vector<vector<int>>>(6, vector<vector<int>>(8, vector<int>(16, 0)));


    // All histograms for characterization
    TH1F *h_dprop = new TH1F("h_dprop", "Propagation Distance Along MDT tubes;Propagation Distance(mm)", 200, 0, 10000);
    TH2F *h2_tubelength_dprop = new TH2F("tubelength_dprop", "Propagation Distance vs. Tube Length; Tube Length(mm); Propagation Distance(mm)",50,0,5000,100,-5000,5000);
    TH1F *h_phi = new TH1F("phi", "Phi Distribution of Muons; Phi; Entries", 252, -3.15, 3.15);
    TH1F *h_delta_phi = new TH1F("Phi bending", "Change in Phi from Solenoid Magnet; Delta Phi; Entries", 50, -0.2, 0.2);
    TH1F *h_exit_phi = new TH1F("Exit Phi", "Phi After Solenoid; Phi; Entries", 252, -3.15, 3.15);
    TH1F *h_xangle = new TH1F("xangle", "X angle; x; Entries", 40, -1, 1);
    TH1F *h_temp_phi = new TH1F("Temp Phi 1D", "Rotated Phi Distribution of Muons; Phi; Entries", 252, -3.15, 3.15);
    TH1F *h_a = new TH1F("a", "Displacement from tube center; a(mm); Entries", 800,-4000,4000);
    TH2F *h2_dr = new TH2F("Global Drift Radius ADC Distribution", "ADC Distribution vs. Drift Radius; Drift Radius (mm); ADC Counts", 15, 0, 15, 400, 0, 400);
    TH2F *h2_p = new TH2F("Global Momentum ADC Distribution", "ADC Distribution vs. Momentum All Barrel Hits; Momentum(GeV) ; ADC Counts", 16, 20, 100, 400, 0, 400);
    TH2F *h2_track_p = new TH2F("Global Momentum Track MPV", "ADC Distribution vs. Momentum All Tracks; Momentum(GeV); Track MPV", 16, 20, 100, 400, 0, 400);
    TH2F *h2_adc = new TH2F("adc_dprop", "ADC Distribution vs. Propagation Distance; Propagation Distance; ADC Counts", 25, 0, 5000, 400, 0, 400);
    
    // ADC vs. propagation distance distributions for all different chamber types
    TH2F *h2_adc1 = new TH2F("adc_dprop1", "ADC Distribution vs. Propagation Distance; Propagation Distance; ADC Counts", 50, 0, 5000, 400, 0, 400);
    TH2F *h2_adc2 = new TH2F("adc_dprop2", "ADC Distribution vs. Propagation Distance; Propagation Distance; ADC Counts", 50, 0, 5000, 400, 0, 400);
    TH2F *h2_adc3 = new TH2F("adc_dprop3", "ADC Distribution vs. Propagation Distance; Propagation Distance; ADC Counts", 50, 0, 5000, 400, 0, 400);
    TH2F *h2_adc4 = new TH2F("adc_dprop4", "ADC Distribution vs. Propagation Distance; Propagation Distance; ADC Counts", 50, 0, 5000, 400, 0, 400);
    TH2F *h2_adc5 = new TH2F("adc_dprop5", "ADC Distribution vs. Propagation Distance; Propagation Distance; ADC Counts", 50, 0, 5000, 400, 0, 400);
    TH2F *h2_adc6 = new TH2F("adc_dprop6", "ADC Distribution vs. Propagation Distance; Propagation Distance; ADC Counts", 50, 0, 5000, 400, 0, 400);
    TH2F *h2_adc7 = new TH2F("adc_dprop7", "ADC Distribution vs. Propagation Distance; Propagation Distance; ADC Counts", 50, 0, 5000, 400, 0, 400);
    TH2F *h2_adc8 = new TH2F("adc_dprop8", "ADC Distribution vs. Propagation Distance; Propagation Distance; ADC Counts", 50, 0, 5000, 400, 0, 400);
    TH2F *h2_adc9 = new TH2F("adc_dprop9", "ADC Distribution vs. Propagation Distance; Propagation Distance; ADC Counts", 50, 0, 5000, 400, 0, 400);
    TH2F *h2_adc10 = new TH2F("adc_dprop10", "ADC Distribution vs. Propagation Distance; Propagation Distance; ADC Counts", 50, 0, 5000, 400, 0, 400);
    TH2F *h2_adc11 = new TH2F("adc_dprop11", "ADC Distribution vs. Propagation Distance; Propagation Distance; ADC Counts", 50, 0, 5000, 400, 0, 400);
    TH2F *h2_adc12 = new TH2F("adc_dprop12", "ADC Distribution vs. Propagation Distance; Propagation Distance; ADC Counts", 50, 0, 5000, 400, 0, 400);
    TH2F *h2_adc13 = new TH2F("adc_dprop13", "ADC Distribution vs. Propagation Distance; Propagation Distance; ADC Counts", 50, 0, 5000, 400, 0, 400);
    TH2F *h2_adc14 = new TH2F("adc_dprop14", "ADC Distribution vs. Propagation Distance; Propagation Distance; ADC Counts", 50, 0, 5000, 400, 0, 400);
    TH2F *h2_adc15 = new TH2F("adc_dprop15", "ADC Distribution vs. Propagation Distance; Propagation Distance; ADC Counts", 50, 0, 5000, 400, 0, 400);
    TH2F *h2_adc16 = new TH2F("adc_dprop16", "ADC Distribution vs. Propagation Distance; Propagation Distance; ADC Counts", 50, 0, 5000, 400, 0, 400);
    TH2F *h2_adc17 = new TH2F("adc_dprop17", "ADC Distribution vs. Propagation Distance; Propagation Distance; ADC Counts", 50, 0, 5000, 400, 0, 400);
    TH2F *h2_exit_coordinates = new TH2F("Exit Coordinates", "(x, y) coordinates on edge of solenoid; x(mm); y(mm)", 260, -1300, 1300, 260, -1300, 1300);
    TH2F *h2_temp_phi = new TH2F("Temp Phi 2D", "Mock Phi vs. Sector; Sector; Phi", 16,0,16,40,-3.15,3.15);

    TH1F *h_tube1 = new TH1F("tube entries", "Number of Entries in Random Tubes; Tube; Entries", 12, 0, 12);
    TH3F *h3_rotated_phi = new TH3F("Mock Phi", "Mock Phi vs. Sector; Muon Phi; Rotated Phi; Sector",40,-3.15,3.15, 40,-3.15,3.15, 16,0,16);
    TH2F *h2_long_propd_eta = new TH2F("Unphysical Propagation Distance Eta", "Eta for Unphysical Propagation Distances; Propagation Distance (mm); Muon Eta", 200, -10000, 10000, 20,-2,2);
    TH2F *h2_long_propd_phi = new TH2F("Unphysical Propagation Distance Phi", "Phi for Unphysical Propagation Distances; Propagation Distance (mm); Muon Phi0", 200, -10000, 10000, 20,-3.15,3.15);
    TH2F *h2_long_propd_pt = new TH2F("Unphysical Propagation Distance pt", "pt for Unphysical Propagation Distances; Propagation Distance (mm); pt (GeV)", 200, -10000, 10000, 20,0,100);
    TH2F *h2_long_propd_delta_phi = new TH2F("Unphysical Propagation Distance Delta Phi", "Change in Phi for Unphysical Propagation Distances; Propagation Distance (mm); Delta Phi", 200, -10000, 10000, 400,-2,2);
    TH2F *h2_long_propd_d0 = new TH2F("Unphysical Propagation Distance d0", "d0 for Unphysical Propagation Distances; Propagation Distance (mm); d0 (mm)", 200, -10000, 10000, 1000,-100,100);
    TH2F *h2_long_propd_z0 = new TH2F("Unphysical Propagation Distance z0", "z0 for Unphysical Propagation Distances; Propagation Distance (mm); z0 (mm)", 200, -10000, 10000, 1000,-100,100);
    TH3F *h3_long_propd_chamber = new TH3F("Unphysical Propagation Distance by Chamber", "CHamber Info; Chamber Eta; Chamber Phi; Station Index", 16,-8,8, 8,0,8, 6,0,6);

    // Make histograms for the chamber of interest in 2 p bins (20-26 GeV) and (55-85 GeV)
    TH2F *h_dprop_low = new TH2F("dprop_low", "ADC Distribution vs. Propagation Distance in Chamber |BIL:Sector 5:Eta 2| (20-26GeV Muons); Propagation Distance; ADC Counts", 25, 0, 5000, 400, 0, 400);
    TH2F *h_dprop_high = new TH2F("dprop_high", "ADC Distribution vs. Propagation Distance in Chamber |BIL:Sector 5:Eta 2| (30-45GeV Muons); Propagation Distance; ADC Counts", 25, 0, 5000, 400, 0, 400);
    TH2F *h_dr_low = new TH2F("dr_low", "ADC Distribution vs. Drift radius in Chamber |BIL:Sector 5:Eta 2| (20-26GeV Muons); Drift radius; ADC Counts", 15, 0, 15, 400, 0, 400);
    TH2F *h_dr_high = new TH2F("dr_high", "ADC Distribution vs. Drift radius in Chamber |BIL:Sector 5:Eta 2| (30-45 GeV Muons); Drift radius; ADC Counts", 15, 0, 15, 400, 0, 400);
    TH2F *h_small_dprop = new TH2F("small_dprop", "ADC Distribution vs. Drift radius in Chamber |BIL:Sector 5:Eta 2| (30-45GeV Muons) (Propagation Distance 500-1000mm)); Drift radius; ADC Counts", 15, 0, 15, 400, 0, 400);
    TH2F *h_small_dprop_small_dr = new TH2F("small_dprop_small_dr", "ADC Distribution vs. Momentum in Chamber |BIL:Sector 5:Eta 2| (Drift Radius 3-7mm) (Propagation Distance 500-1000mm)); Momentum(GeV) ; ADC Counts", 16, 20, 100, 400, 0, 400);

    // A Histogram of global adc distribution down sampled for similar statistics to isolated bin
    TH2F *h2_downsampled_adc_dr = new TH2F("downsampled ADC vs. Drift Radius", "Global ADC Distribution vs. Drift Radius (Downsampled); Drift Radius(mm); ADC Counts", 15, 0, 15, 400, 0, 400);
    TH2F *h2_downsampled_adc_p = new TH2F("downsampled ADC vs. Momentum", "Global ADC Distribution vs. Momentum (Downsampled); Momentum(GeV); ADC Counts", 20, 0, 100, 400, 0, 400);

    std::vector<float> temp_track;
    std::vector<float> temp_track_stations;
    int small_dprop_dr_count = 0;

    // Loop trhough teh TChain
    Long64_t nEntries = chain.GetEntries();
    for (Int_t i=0; i < nEntries; i++) {

        // A progress indicator
        if (i % 10000 == 0) {
            std::cout << i << "/" << nEntries << " processed" << endl;
        }

        // Get the event data
        chain.GetEntry(i);

        // Make sure the current muon has some ADC counts recorded
        if (ADC_counts->size() == 0) {
            continue;
        }

        // Mark the current muon link
        int muon_link = 0;
        float muon_pt = pts->at(muon_link);
        float muon_eta = etas->at(muon_link);
        float muon_phi = phis->at(muon_link);
        float muon_d0 = d0s->at(muon_link);
        float muon_z0 = z0s->at(muon_link);
        int muon_charge = charges->at(muon_link);
        float muon_p = ptToP(muon_pt, muon_eta);
        

        // Loop through the muon hnits
        for (unsigned int n=0; n<ADC_counts->size(); n++) {
            
            // SKip empty muons
            if (ADC_counts->size()==0) { continue; }

            // Fill the hit into the current track
            temp_track.push_back(ADC_counts->at(n));
            temp_track_stations.push_back(station_index->at(n));

            // Skip hits not in the barrel
            if (station_index->at(n) > 5) {
                continue;
            }

            // Fill the general adc vs drift radius distribution
            h2_dr->Fill(drift_radius->at(n), ADC_counts->at(n));
            h2_p->Fill(muon_p, ADC_counts->at(n));

            // Downsample
            int random_number = dis(gen);
            if (random_number == 1411) {
                h2_downsampled_adc_dr->Fill(drift_radius->at(n), ADC_counts->at(n));
            }
            int random_number2 = dis2(gen2);
            if (random_number2 == 1411) {
                h2_downsampled_adc_p->Fill(muon_p, ADC_counts->at(n));
            }

            // Check of we have a new muon
            if (trk_link->at(n) != muon_link) {
                std::sort(temp_track_stations.begin(), temp_track_stations.end());
                if ( temp_track_stations.at(temp_track_stations.size()-1) < 6) { // The track is contained in the barrel if this is met
                    h2_track_p->Fill(muon_p, TM(temp_track, 0.6));
                    temp_track.clear();
                    temp_track_stations.clear();
                }
                
                // Update muon
                muon_link = trk_link->at(n);
                muon_pt = pts->at(muon_link);
                muon_eta = etas->at(muon_link);
                muon_phi = phis->at(muon_link);
                muon_p = ptToP(muon_pt, muon_eta);
            }

            // Skip outlier hits
            if (hit_type->at(n) > 60) { continue; }

            // Get the length of the current mdt tube
            int chamber_type = chamber_types[station_index->at(n)][abs(station_eta->at(n))][station_phi->at(n)];
            int tube_length = tube_lengths[station_index->at(n)][chamber_type];
            
            // Make sure that the chamber is instantiated
            if (chamber_type == 0) { continue; }

            // Find the location of the tube where the hit occured
            float radial_distance = tube_positions[station_index->at(n)][multi_layer->at(n)][tube_layer->at(n)];
            int sector = stationPhiToSector(station_index->at(n), station_phi->at(n));

            // use this temp phi to rotate all of the tracks as if they are moving vertically for geometry purposes 
            float temp_phi;
            if (muon_phi > 0) {
                temp_phi = muon_phi - (sector-1) *(2*M_PI / 16) + M_PI/2;
            } 
            else if (muon_phi < 0) {
                temp_phi = muon_phi + (2*M_PI - (sector-1) *(2*M_PI / 16)) + M_PI/2;
            }

            // Keep track of hits per chamber
            if (station_eta->at(n) < 0) {
                chambers.at(static_cast<int>(station_index->at(n))).at(static_cast<int>(station_phi->at(n)-1)).at(static_cast<int>(station_eta->at(n)+8)) += 1;
            } else if (station_eta->at(n) > 0) {
                chambers.at(static_cast<int>(station_index->at(n))).at(static_cast<int>(station_phi->at(n)-1)).at(static_cast<int>(station_eta->at(n)+7)) += 1;
            }

            // Compute the bending of the particle in phi before the MS and get the propagation distance
            float larmor_radius = muon_pt * (5.344 * pow(10, -19)) / ((1.602 * pow(10, -19)) * (2)) * 1000; // In mm
            float ID_radius = 1150; //mm
            float length_through_solenoid = ID_radius / sin(2 * atan(exp(-muon_eta)));
            float delta_phi = muon_charge * 4 * asin( sqrt( (larmor_radius - sqrt(pow(larmor_radius, 2) - 0.25 * pow(length_through_solenoid, 2))) / (2 * larmor_radius)) );
            float exit_phi = temp_phi + delta_phi; 
            float exit_x = ID_radius * cos(temp_phi + (0.5 * delta_phi));
            float exit_y = ID_radius * sin(temp_phi + (0.5 * delta_phi));
            float slope = tan(exit_phi);
            float y_intercept = exit_y - slope * exit_x;
            float theta = M_PI / 2;
            float x = temp_phi - theta;
            float a = (radial_distance - y_intercept) / tan(exit_phi); // This is the displacement from the the middle of the MDT along its central axis
            float distance_to_propagate;

            // The readout is on opposite sides for large sectors and small sectors
            if (station_index->at(n) % 2 == 0) {
                distance_to_propagate = (tube_length/2) + a;
            }
            else if (station_index->at(n) % 2 == 1) {
                distance_to_propagate = (tube_length/2) - a;
            }

            // Skip unphysical hits
            if (distance_to_propagate < 0 || distance_to_propagate > tube_length) {
                continue;
            }

            // Fill binned histograms for specific chamber
            if (station_index->at(n)==0 && station_phi->at(n)==3 && station_eta->at(n)==2) {
                if (muon_p > 20 && muon_p < 26) { 
                    h_dprop_low->Fill(distance_to_propagate, ADC_counts->at(n));
                    h_dr_low->Fill(drift_radius->at(n), ADC_counts->at(n));
                } else if (muon_p > 30 && muon_p < 45) {
                    h_dprop_high->Fill(distance_to_propagate, ADC_counts->at(n));
                    h_dr_high->Fill(drift_radius->at(n), ADC_counts->at(n));
                    if (distance_to_propagate > 500 && distance_to_propagate < 1000) {
                        h_small_dprop->Fill(drift_radius->at(n), ADC_counts->at(n));
                    }
                }
            }

            // Include adjacent stations as well
            if ((station_index->at(n)==0 && station_phi->at(n)==3) && station_eta->at(n)==2) {// (station_eta->at(n)==1 || station_eta->at(n)==2 || station_eta->at(n)==3)) {
                if (distance_to_propagate > 500 && distance_to_propagate < 1000) {
                    if (drift_radius->at(n) > 3.0 && drift_radius->at(n) < 7.0) {
                        small_dprop_dr_count++;
                        h_small_dprop_small_dr->Fill(muon_p, ADC_counts->at(n));
                    }
                }
            }

            // Fill histograms
            h_temp_phi->Fill(temp_phi);
            h2_temp_phi->Fill(sector, temp_phi);
            h3_rotated_phi->Fill(muon_phi, temp_phi, sector);
            h_dprop->Fill(distance_to_propagate);
            h2_tubelength_dprop->Fill(tube_length, distance_to_propagate);
            h_xangle->Fill(x);
            h_phi->Fill(temp_phi);
            h_delta_phi->Fill(delta_phi);
            h_exit_phi->Fill(exit_phi);
            h_a->Fill(a);
            h2_adc->Fill(distance_to_propagate, ADC_counts->at(n));
            h2_exit_coordinates->Fill(exit_x, exit_y);

            // Fill Propagation distance histograms for specific chamber types
            if (station_index->at(n)==0){
                if (chamber_type==1) {
                    h2_adc1->Fill(distance_to_propagate, ADC_counts->at(n));
                }
                else if (chamber_type==2) {
                    h2_adc2->Fill(distance_to_propagate, ADC_counts->at(n));
                }
                if (chamber_type==3) {
                    h2_adc3->Fill(distance_to_propagate, ADC_counts->at(n));
                }
                else if (chamber_type==4) {
                    h2_adc4->Fill(distance_to_propagate, ADC_counts->at(n));
                }
            }
            else if (station_index->at(n)==1) {
                if (chamber_type==1) {
                    h2_adc5->Fill(distance_to_propagate, ADC_counts->at(n));
                }
                else if (chamber_type==2) {
                    h2_adc6->Fill(distance_to_propagate, ADC_counts->at(n));
                }
                if (chamber_type==3) {
                    h2_adc7->Fill(distance_to_propagate, ADC_counts->at(n));
                }
            }
            else if (station_index->at(n)==2) {
                if (chamber_type==1) {
                    h2_adc8->Fill(distance_to_propagate, ADC_counts->at(n));
                }
                else if (chamber_type==2) {
                    h2_adc9->Fill(distance_to_propagate, ADC_counts->at(n));
                }
            }
            else if (station_index->at(n)==3) {
                if (chamber_type==1) {
                    h2_adc10->Fill(distance_to_propagate, ADC_counts->at(n));
                }
                else if (chamber_type==2) {
                    h2_adc11->Fill(distance_to_propagate, ADC_counts->at(n));
                }
                if (chamber_type==3) {
                    h2_adc12->Fill(distance_to_propagate, ADC_counts->at(n));
                }
            }
            else if (station_index->at(n)==4) {
                if (chamber_type==1) {
                    h2_adc13->Fill(distance_to_propagate, ADC_counts->at(n));
                }
                else if (chamber_type==2) {
                    h2_adc14->Fill(distance_to_propagate, ADC_counts->at(n));
                }
            }
            else if (station_index->at(n)==5) {
                if (chamber_type==1) {
                    h2_adc15->Fill(distance_to_propagate, ADC_counts->at(n));
                }
                else if (chamber_type==2) {
                    h2_adc16->Fill(distance_to_propagate, ADC_counts->at(n));
                }
                if (chamber_type==3) {
                    h2_adc17->Fill(distance_to_propagate, ADC_counts->at(n));
                }
            }

            // Fill some random tubes for tube entry
            if (station_index->at(n) == 0) {
                if (station_eta->at(n) == 1 && station_phi->at(n) == 6 && tubes->at(n) == 10) {
                    h_tube1->Fill(0);
                }
                if (station_eta->at(n) == -4 && station_phi->at(n) == 2 && tubes->at(n) == 10) {
                    h_tube1->Fill(1);
                }
            }
            else if (station_index->at(n) == 1) {
                if (station_eta->at(n) == 2 && station_phi->at(n) == 7 && tubes->at(n) == 10) {
                    h_tube1->Fill(2);
                }
                if (station_eta->at(n) == -6 && station_phi->at(n) == 4 && tubes->at(n) == 10) {
                    h_tube1->Fill(3);
                }
            }
            else if (station_index->at(n) == 2) {
                if (station_eta->at(n) == -3 && station_phi->at(n) == 1 && tubes->at(n) == 10) {
                    h_tube1->Fill(4);
                }
                if (station_eta->at(n) == 5 && station_phi->at(n) == 8 && tubes->at(n) == 10) {
                    h_tube1->Fill(5);
                }
            }
            else if (station_index->at(n) == 3) {
                if (station_eta->at(n) == -2 && station_phi->at(n) == 2 && tubes->at(n) == 10) {
                    h_tube1->Fill(6);
                }
                if (station_eta->at(n) == 4 && station_phi->at(n) == 7 && tubes->at(n) == 10) {
                    h_tube1->Fill(7);
                }
            }
            else if (station_index->at(n) == 4) {
                if (station_eta->at(n) == 1 && station_phi->at(n) == 4 && tubes->at(n) == 10) {
                    h_tube1->Fill(8);
                }
                if (station_eta->at(n) == 5 && station_phi->at(n) == 1 && tubes->at(n) == 10) {
                    h_tube1->Fill(9);
                }
            }
            else if (station_index->at(n) == 5) {
                if (station_eta->at(n) == -3 && station_phi->at(n) == 3 && tubes->at(n) == 10) {
                    h_tube1->Fill(10);
                }
                if (station_eta->at(n) == 2 && station_phi->at(n) == 5 && tubes->at(n) == 10) {
                    h_tube1->Fill(11);
                }
            }


            /***
            if (distance_to_propagate > tube_length || distance_to_propagate < 0) {
                h2_long_propd_eta->Fill(distance_to_propagate, muon_eta);
                h2_long_propd_phi->Fill(distance_to_propagate, temp_phi);
                h2_long_propd_pt->Fill(distance_to_propagate, muon_pt);
                h2_long_propd_delta_phi->Fill(distance_to_propagate, delta_phi);
                h2_long_propd_d0->Fill(distance_to_propagate, muon_d0);
                h2_long_propd_z0->Fill(distance_to_propagate, muon_z0);
                h3_long_propd_chamber->Fill(station_eta->at(n), station_phi->at(n), station_index->at(n));
            }
            ***/
        }

        // FIll the track here if it is the end of the event but there is no change in muon link(last muon)
        
        if ( temp_track_stations.at(temp_track_stations.size()-1) < 6) { // The track is contained in the barrel if this is met
            h2_track_p->Fill(muon_p, TM(temp_track, 0.6));
            temp_track.clear();
            temp_track_stations.clear();
        }
    }

    // Print out the number of entries ine very chamber
    int maxe = 0;
    for (unsigned int i=0; i<chambers.size(); i++) {
        for (unsigned int j=0; j<chambers.at(i).size(); j++) {
            for (unsigned int k=0; k<chambers.at(i).at(j).size(); k++) {
                cout << "There are " << chambers.at(i).at(j).at(k) << " entries in stationIndex: " << i << " stationPhi: " << j << "stationEta: " << k << endl;
                maxe = max(maxe, chambers.at(i).at(j).at(k));
            }
        }
    }
    cout << "Max hits in a chamber is " << maxe << endl;

    // Create the output file
    TFile *output_file = new TFile(outfile.c_str(), "RECREATE");

    std::vector<float> adc_dprop_TMs = fit2DHistogramTM(h2_adc);
    TGraph *g_adc = new TGraph();
    for (unsigned int i=0; i<adc_dprop_TMs.size(); i++) {
        g_adc->AddPoint((200*i) + 100, adc_dprop_TMs.at(i));
    }
    std::vector<float> dprop_low_TMs = fit2DHistogramTM(h_dprop_low);
    TGraph *g_dprop_low = new TGraph();
    for (unsigned int i=0; i<dprop_low_TMs.size(); i++) {
        g_dprop_low->AddPoint((200*i) + 100, dprop_low_TMs.at(i));
    }
    std::vector<float> dprop_high_TMs = fit2DHistogramTM(h_dprop_high);
    TGraph *g_dprop_high = new TGraph();
    for (unsigned int i=0; i<dprop_high_TMs.size(); i++) {
        g_dprop_high->AddPoint((200*i) + 100, dprop_high_TMs.at(i));
    }
    std::vector<float> dr_low_TMs = fit2DHistogramTM(h_dr_low);
    TGraph *g_dr_low = new TGraph();
    for (unsigned int i=0; i<dr_low_TMs.size(); i++) {
        g_dr_low->AddPoint((200*i) + 100, dr_low_TMs.at(i));
    }
    std::vector<float> dr_high_TMs = fit2DHistogramTM(h_dr_high);
    TGraph *g_dr_high = new TGraph();
    for (unsigned int i=0; i<dr_high_TMs.size(); i++) {
        g_dr_high->AddPoint((200*i) + 100, dr_high_TMs.at(i));
    }

    // Find the scale factor behavior for ADC in isolated momentum, propagation distance, and chamber
    std::vector<std::vector<float>>small_dprop_params = fit2DHistogramLandau(h_small_dprop);
    TGraphErrors *g_small_dprop_scale_factors = new TGraphErrors();
    for (unsigned int i=0; i<small_dprop_params.size(); i++) {
        g_small_dprop_scale_factors->AddPointError((i * 1) + 0.5, small_dprop_params.at(i).at(1), 0, small_dprop_params.at(i).at(3));
    }
    g_small_dprop_scale_factors->SetTitle("Landau Scale Factor for ADC Distribution in |BIL:Sector 5:Eta 1-3| (30-45GeV Muons) (Propagation Distance 500-1000mm); Drift Radius(mm); ADC Distribution Scale Factor");
    g_small_dprop_scale_factors->SetMarkerSize(0.5);
    g_small_dprop_scale_factors->SetMarkerStyle(8);
    g_small_dprop_scale_factors->Write();

    // FInd Scale Factors and MPVs for adc vs. Momentum in isolated drift radius, propagation distance and chamber
    std::vector<std::vector<float>>small_dprop_small_dr_params = fit2DHistogramLandau(h_small_dprop_small_dr, "DRAW");
    TGraphErrors *g_small_dprop_small_dr_scale_factors = new TGraphErrors();
    TGraphErrors *g_small_dprop_small_dr_mpvs = new TGraphErrors();
    for (unsigned int i=0; i<small_dprop_small_dr_params.size(); i++) {
        g_small_dprop_small_dr_scale_factors->AddPointError((i * 5) + 22.5, small_dprop_small_dr_params.at(i).at(1), 0, small_dprop_small_dr_params.at(i).at(3));
        g_small_dprop_small_dr_mpvs->AddPointError((i * 5) + 22.5, small_dprop_small_dr_params.at(i).at(0), 0, small_dprop_small_dr_params.at(i).at(2));
    }
    g_small_dprop_small_dr_scale_factors->SetTitle("Landau Scale Factor for ADC Distribution in |BIL:Sector 5:Eta 1-3| (Drift Radius 3-7mm) (Propagation Distance 500-1000mm); Momentum (GeV); ADC Distribution Scale Factor");
    g_small_dprop_small_dr_scale_factors->SetMarkerSize(0.5);
    g_small_dprop_small_dr_scale_factors->SetMarkerStyle(8);
    g_small_dprop_small_dr_scale_factors->Write();
    g_small_dprop_small_dr_mpvs->SetTitle("Landau MPV for ADC/Momentum Distribution in |BIL:Sector 5:Eta 1-3| (Drift Radius 3-7mm) (Propagation Distance 500-1000mm); Momentum(GeV); ADC Distribution MPV");
    g_small_dprop_small_dr_mpvs->SetName("MPV in |BIL:Sector 5:Eta 1-3| (Drift Radius 3-7mm) (Propagation Distance 500-1000mm)");
    g_small_dprop_small_dr_mpvs->SetLineStyle(0);
    g_small_dprop_small_dr_mpvs->SetMarkerSize(0.5);
    g_small_dprop_small_dr_mpvs->SetMarkerStyle(8);
    g_small_dprop_small_dr_mpvs->Write();
    
    // Find the scale factor behavior for ADC in isolated momentum, propagation distance, and chamber
    std::vector<std::vector<float>>all_hits_params = fit2DHistogramLandau(h2_dr);
    TGraphErrors *g_all_hits_scale_factors = new TGraphErrors();
    for (unsigned int i=0; i<all_hits_params.size(); i++) {
        g_all_hits_scale_factors->AddPointError((i * 1) + 0.5, all_hits_params.at(i).at(1), 0, all_hits_params.at(i).at(3));
    }
    g_all_hits_scale_factors->SetTitle("Landau Scale Factor for ADC Distribution in Barrel; Drift Radius(mm); ADC Distribution Scale Factor");
    g_all_hits_scale_factors->SetMarkerSize(0.5);
    g_all_hits_scale_factors->SetMarkerStyle(8);
    g_all_hits_scale_factors->Write();
    
    // Find the MPV behavior for ADC in all hits
    std::vector<std::vector<float>>all_hits_params_p_dist = fit2DHistogramLandau(h2_p);
    TGraphErrors *g_all_hits_mpvs = new TGraphErrors();
    for (unsigned int i=0; i<all_hits_params_p_dist.size(); i++) {
        g_all_hits_mpvs->AddPointError((i * 5) + 22.5, all_hits_params_p_dist.at(i).at(0), 0, all_hits_params_p_dist.at(i).at(2));
    }
    g_all_hits_mpvs->SetTitle("Landau MPV for ADC/Momentum Distribution in Barrel; Momentum (GeV); ADC Distribution MPV");
    g_all_hits_mpvs->SetName("MPV from all barrel hits");
    g_all_hits_mpvs->SetMarkerSize(0.5);
    g_all_hits_mpvs->SetMarkerStyle(8);
    g_all_hits_mpvs->SetLineStyle(0);
    g_all_hits_mpvs->Write();

    // Find the scale factors for global adc distribution at varius drift radius downsampled to match isolated bin
    std::vector<std::vector<float>>all_downsampled_hits_params_p_dist = fit2DHistogramLandau(h2_downsampled_adc_p);
    TGraphErrors *g_all_downsampled_hits_mpvs = new TGraphErrors();
    for (unsigned int i=0; i<all_downsampled_hits_params_p_dist.size(); i++) {
        g_all_downsampled_hits_mpvs->AddPointError((i * 5) + 22.5, all_downsampled_hits_params_p_dist.at(i).at(0), 0, all_downsampled_hits_params_p_dist.at(i).at(2));
    }
    g_all_downsampled_hits_mpvs->SetTitle("Landau MPV for Downsampled ADC/Momentum Distribution in Barrel; Momentum (GeV); ADC Distribution MPV");
    g_all_downsampled_hits_mpvs->SetName("MPV from downsampled barrel hits");
    g_all_downsampled_hits_mpvs->SetMarkerSize(0.5);
    g_all_downsampled_hits_mpvs->SetMarkerStyle(8);
    g_all_downsampled_hits_mpvs->SetLineStyle(0);
    g_all_downsampled_hits_mpvs->Write();

    // Find the scale factors for global adc distribution at varius drift radius downsampled to match isolated bin
    std::vector<std::vector<float>>all_downsampled_hits_params = fit2DHistogramLandau(h2_downsampled_adc_dr);
    TGraphErrors *g_all_downsampled_hits_scale_factors = new TGraphErrors();
    for (unsigned int i=0; i<all_downsampled_hits_params.size(); i++) {
        g_all_downsampled_hits_scale_factors->AddPointError((i * 1) + 0.5, all_downsampled_hits_params.at(i).at(1), 0, all_downsampled_hits_params.at(i).at(3));
    }
    g_all_downsampled_hits_scale_factors->SetTitle("Landau MPV for Downsampled ADC/Drift Radius Distribution in Barrel; Drift Radius(mm); ADC Distribution MPV");
    g_all_downsampled_hits_scale_factors->SetMarkerSize(0.5);
    g_all_downsampled_hits_scale_factors->SetMarkerStyle(8);
    g_all_downsampled_hits_scale_factors->Write();

    // Create overlayed MPVs
    TCanvas *c_mpvs = new TCanvas("MPV Behaviors", "MPV vs. Momentum; Momentum(GeV); ADC MPV", 800, 600);
    g_small_dprop_small_dr_mpvs->SetLineColor(kRed);
    g_small_dprop_small_dr_mpvs->SetMarkerColor(kRed);
    g_all_hits_mpvs->SetLineColor(kBlue);
    g_all_hits_mpvs->SetMarkerColor(kBlue);
    g_small_dprop_small_dr_mpvs->Draw("AP");
    g_all_hits_mpvs->Draw("Psame");
    c_mpvs->Write();

    // Fit a line to the adc vs. propagation distance distribution
    g_adc->Fit("pol1", "S");
    TF1* fitResult = g_adc->GetFunction("pol1");
    fitResult->SetLineColor(kRed);
    fitResult->SetLineStyle(2);
    gStyle->SetOptFit(1111);
    g_adc->SetTitle("ADC 40% Truncated Mean vs. Propagation Distance");
    g_adc->GetXaxis()->SetTitle("Propagation Distance(mm)");
    g_adc->GetYaxis()->SetTitle("ADC 40% Truncated Mean");
    TCanvas* c1 = new TCanvas("ADC vs Propagation Distance", "ADC vs Distance;  Propagation Distance (mm); ADC Count", 800, 600);
    c1->SetTitle("ADC Truncated Mean vs. Propagation Distance; Propagation Distance; ADC Truncated Mean");
    g_adc->SetMarkerSize(8);
    g_adc->Draw("AP");
    fitResult->Draw("same");
    g_adc->SetMarkerSize(8);
    c1->Write();

    // Write the chamber of interest histograms
    h_dprop_low->Write();
    h_dprop_high->Write();
    h_dr_low->Write();
    h_dr_high->Write();
    h_small_dprop->Write();
    h_small_dprop_small_dr->Write();
    h2_downsampled_adc_dr->Write();
    h2_downsampled_adc_p->Write();
    h2_p->Write();
    h2_track_p->Write();
    h2_dr->Write();

    // Write the raw histograms
    h_dprop->Write();
    h2_tubelength_dprop->SetOption("COLZ");
    h2_tubelength_dprop->Write();
    h_phi->Write();
    h2_temp_phi->Write();
    h_delta_phi->Write();
    h_exit_phi->Write();
    h_xangle->Write();
    h_a->Write();
    h_tube1->Write();
    h2_adc->SetOption("COLZ");
    h2_adc1->Write();
    h2_adc2->Write();
    h2_adc3->Write();
    h2_adc4->Write();
    h2_adc5->Write();
    h2_adc6->Write();
    h2_adc7->Write();
    h2_adc8->Write();
    h2_adc9->Write();
    h2_adc10->Write();
    h2_adc11->Write();
    h2_adc12->Write();
    h2_adc13->Write();
    h2_adc14->Write();
    h2_adc15->Write();
    h2_adc16->Write();
    h2_adc17->Write();
    h2_exit_coordinates->Write();
    g_adc->Write();
    h3_rotated_phi->Write();
    h2_long_propd_eta->Write();
    h2_long_propd_phi->Write();
    h2_long_propd_pt->Write();
    h2_long_propd_delta_phi->Write();
    h2_long_propd_d0->Write();
    h2_long_propd_z0->Write();
    h3_long_propd_chamber->Write();
    output_file->Close();
    return;
}

int main(int argc, char* argv[]) {
    // Check for correct number of command line arguments
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_directory> <output_file>" << std::endl;
        return 1;
    }
    // Parse the command line arguments
    std::string indir(argv[1]);
    std::string outfile(argv[2]);

    // Execute the primary function
    featureSpace(indir, outfile);
    return 0;
}