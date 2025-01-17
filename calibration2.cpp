#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
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
#include <TGraph2DErrors.h>
#include <TLine.h>
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
#include <TF2.h>
#include <TF3.h>
#include <TMultiDimFit.h>
#include <TPaveText.h>
#include <TLegend.h>
using namespace std;

// A function to evaluate the 2D correction function for a given ADC
float correctionFunction(float adc, float drift_radius, float propagation_distance, std::vector<float> params) {
    return (params.at(0) * TMath::Exp(params.at(1)*propagation_distance)) * (params.at(2) + drift_radius*params.at(3)
     + TMath::Power(drift_radius, 2)*params.at(4) + TMath::Power(drift_radius, 3)*params.at(5) + 
     TMath::Power(drift_radius, 4)*params.at(6));
}

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

// Return the truncated mean of a 1D histogram
float fit1DHistogramTM(TH1 *h1) {
    int nbins = h1->GetNbinsX();
    double hist_entries = h1->GetEntries();
    if (hist_entries==0) { return 0; }
    int cutoff = std::floor(hist_entries * 0.6);
    int entry_count = 0;
    float tm = 0;
    
    for (int i=1; i<=nbins; i++) {
        float binContent = h1->GetBinContent(i);
        float binLowEdge = h1->GetBinLowEdge(i);
        float binWidth = h1->GetBinWidth(i);
        float binCenter = binLowEdge + 0.5 * binWidth;
        if (entry_count + binContent < cutoff) {
            tm += binCenter * binContent;
            entry_count += binContent;
        }
        else {
            float remaining = cutoff - entry_count;
            tm += binCenter * remaining;
            entry_count = cutoff;
            break;
        }
    }
    cout << "Cumulative sum is " << tm << " and cutoff is " << cutoff << " out of " << hist_entries << " entries " << endl;
    tm /= cutoff;
    cout << "Returning " << tm << endl;
    return tm; 
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
            // delete canvas;
            // delete label;
        }
    }
    return landau_params;
}

// Return a Truncated mean from each xBin in a TH2 object
std::vector<std::vector<float>> fit3DHistogramLandau(TH3 *h3, const std::string &option="") {
    int nbinsx = h3->GetNbinsX();
    int nbinsy = h3->GetNbinsY();
    const char* h_title = h3->GetTitle();

    // Record the parameters of the Landau fits <MPV, sigma, MPV_err, sigma_err>
    std::vector<std::vector<float>> landau_params;

    //Loop through the bins
    for (int i = 1; i<=nbinsx; i++) { // Loop through x bins
        for (int j = 1; j<=nbinsy; j++) { // Loop through y bins
            // Turn each xbin into a 1d histogram
            TH1D* projectionZ = h3->ProjectionZ("", i, i, j, j);
            // cout << "before fit1D" << endl;
            float tm = fit1DHistogramTM(projectionZ);

            // cout << "After fit1D" << endl;
            float mpv;
            float sigma;
            float mpvError;
            float sigmaError;

            if (tm == 0) {
                mpv = -1000;
                mpvError = 100000000;
                sigma = 1000;
                sigmaError = 100000000;
                std::vector<float> empty_results = {mpv, sigma, mpvError, sigmaError};
                landau_params.push_back(empty_results);
                continue;
            }
            // Fit a landau to the xbin histogram
            projectionZ->Fit("landau");

            // Retrieve the fit result
            TF1* fitResult = projectionZ->GetFunction("landau");
            mpv = fitResult->GetParameter(1); // Most probable value (MPV)
            sigma = fitResult->GetParameter(2); // Sigma (width of the distribution)
            mpvError = fitResult->GetParError(1); // Error on MPV
            sigmaError = fitResult->GetParError(2); // Error on Sigma
            std::vector<float> results = {mpv, sigma, mpvError, sigmaError};
            landau_params.push_back(results);

            if (option == "DRAW") { // Create a canvas on which to draw the fit
                TCanvas *canvas = new TCanvas(Form("Momentum bin x%d|y%d from %s", i, j, h_title), Form("Momentum Bin x%d|y%d Projection and Landau Fit from %s", i, j, h_title), 800, 600);
                TPaveText *label = new TPaveText(0.6, 0.7, 0.9, 0.9, "NDC");
                label->AddText(Form("Chi^2: %f", fitResult->GetChisquare()));
                label->AddText(Form("NDF: %d", fitResult->GetNDF()));
                label->AddText(Form("MPV: %f +/- %.2f", mpv, mpvError));
                label->AddText(Form("Scale: %f +/- %.2f", sigma, sigmaError));
                label->AddText(Form("40% Truncated Mean: %f", tm));
                projectionZ->SetTitle(Form("Landau fit for Bin x%d|y%d from %s", i, j, h_title));
                projectionZ->Draw();
                fitResult->Draw("same");
                label->Draw("same");
                TLine *line = new TLine(tm, 0, tm, projectionZ->GetMaximum());
                line->SetLineColor(kBlue);
                line->SetLineStyle(2); // Dashed line
                line->Draw("same");
                canvas->Write();
                delete canvas;
                delete label;
            }
        }
    }
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

// The function which we calibrate in every bin
Double_t fitFunc(Double_t *vals, Double_t *pars) {
    Float_t x = vals[0];
    Float_t y = vals[1]; 
    Double_t f = (pars[0] * TMath::Exp(pars[1]*x) + pars[2]) * (pars[3] + y*pars[4] + TMath::Power(y,2)*pars[5] + TMath::Power(y,3)*pars[6] + TMath::Power(y,4)*pars[7]);
    return f;
}


void calibration2(std::string indir, const::std::string& outfile) {

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

    /***
    All histograms for characterization.
    Propagation Distance bins will be 500mm wide from 0-5k
    Drift radius bins from { (0,3.75),(3.75,7.5),(7.5,11.25),(11.25,15) }
    ***/
    int num_dr_bins = 5; // Set the number of bins each spatial bin will have in drift radius
    int num_dprop_bins = 5; // Set the number of bins each spatial bin will have in propagation distance
    std::map<int, std::map<int, std::map<int, TH3F*>>> binHits;
    for (unsigned int h=0; h<3; h++) { // Loop through BIL, BML, BOL
        for (unsigned int i=0; i<4; i++) { // Loop through drift radius bins
            for (unsigned int j=0; j<10; j++) { // Loop through propagation distance bins
                std::stringstream ss;
                ss << "adc_R" << h << "dr" << i << "dprop" << j; 
                std::string hist_name = ss.str();

                ss.str("");
                float dr_start = i*3.75;
                float dr_end = (i+1)*3.75;
                float dprop_start = j*500;
                float dprop_end = (j+1)*500;
                ss << "ADC Distribution for Drift Radius (" << dr_start << "-" << dr_end << ")mm and Propagation Distance (" << dprop_start << "-" << dprop_end << ")mm; Drift Radius(mm); Propagation Distance(mm); ADC Counts";
                std::string hist_title = ss.str();
                binHits[h][i][j] = new TH3F(hist_name.c_str(), hist_title.c_str(), num_dr_bins, dr_start, dr_end, num_dprop_bins, dprop_start, dprop_end, 400, 0, 400);
            }
        }
    }

    // Loop through the TChain
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
        int muon_author = authors->at(muon_link);
        int muon_charge = charges->at(muon_link);
        float muon_p = ptToP(muon_pt, muon_eta);
        int station_bin;
        int dr_bin;
        int dprop_bin;
        

        // Loop through the muon hnits
        for (unsigned int n=0; n<ADC_counts->size(); n++) {

            // SKip empty muons
            if (ADC_counts->size()==0) { continue; }

            // Skip hits not in sector 5 and outlier hits
            if ((station_index->at(n) % 2 )!=0 || station_phi->at(n)==3 || hit_type->at(n) > 60 || std::abs(muon_eta) > 1 || muon_author != 1) {
                continue;
            }

            // Check of we have a new muon
            if (trk_link->at(n) != muon_link) {
                // Update muon
                muon_link = trk_link->at(n);
                muon_pt = pts->at(muon_link);
                muon_eta = etas->at(muon_link);
                muon_phi = phis->at(muon_link);
                muon_p = ptToP(muon_pt, muon_eta);
                muon_author = authors->at(muon_link);
            }

            // Only take muons in the 30-40GeV bin
            if (muon_p < 30 || muon_p > 40) { continue; }
            

            // Get the length of the current mdt tube
            int chamber_type = chamber_types[station_index->at(n)][abs(station_eta->at(n))][station_phi->at(n)];
            int tube_length = tube_lengths[station_index->at(n)][chamber_type];
            // Make sure that the chamber is instantiated
            if (chamber_type == 0) { continue; }

            // Find the location of the tube where the hit occured
            float radial_distance = tube_positions[station_index->at(n)][multi_layer->at(n)][tube_layer->at(n)];
            int sector = stationPhiToSector(station_index->at(n), station_phi->at(n));

            float temp_phi; // use this temp phi to rotate all of the tracks as if they are moving vertically for geometry purposes 
            if (muon_phi > 0) {
                temp_phi = muon_phi - (sector-1) *(2*M_PI / 16) + M_PI/2;
            } 
            else if (muon_phi < 0) {
                temp_phi = muon_phi + (2*M_PI - (sector-1) *(2*M_PI / 16)) + M_PI/2;
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
            float a = (radial_distance - y_intercept) / tan(exit_phi); // This is the displacement from the the middle of the MDT along its central axis
            float distance_to_propagate;

            // The readout is on opposite sides for large sectors and small sectors
            if (station_index->at(n) % 2 == 0) {
                distance_to_propagate = (tube_length/2) + a;
            }
            else if (station_index->at(n) % 2 == 1) {
                distance_to_propagate = (tube_length/2) - a;
            }
            station_bin = station_index->at(n) / 2;
            dr_bin = drift_radius->at(n) / 3.75;
            dprop_bin = distance_to_propagate / 500;
            if (dprop_bin > 9 || dprop_bin < 0 || dr_bin > 3 || dr_bin < 0) { continue; }
            binHits[station_bin][dr_bin][dprop_bin]->Fill(drift_radius->at(n), distance_to_propagate, ADC_counts->at(n));
        }
    }


    // Create the output file
    TFile *output_file = new TFile(outfile.c_str(), "RECREATE");

    // Fit each bin histogram to a 2D function now
    std::map<int, std::map<int, std::map<int, vector<float>>>> binFuncs; // Store the function parameters here
    for (unsigned int h = 0; h < 3; h++) { // Loop through BIL, BML, BOL
        for (unsigned int i = 0; i < 4; i++) { // Loop through drift radius bins
            for (unsigned int j = 0; j < 10; j++) { // Loop through propagation distance bins
                cout << "bin R" << h << "dr" << i << "dprop" << j << " has " << binHits[h][i][j]->GetEntries() << " entries" << endl;
                float xmin = i * 3.75;
                float xmax = (i + 1) * 3.75;
                float ymin = j * 500;
                float ymax = (j + 1) * 500;

                // Get the MPVs in every bin.
                std::vector<std::vector<float>> bin_params = fit3DHistogramLandau(binHits[h][i][j], "DRAW");
                
                // SKip this bin if it has no entries
                if (binHits[h][i][j]->GetEntries()==0) { continue; }

                // Creaye the TGraph2D for MPv's to fit correction function to
                std::stringstream ss;
                ss << "adc_R" << h << "dr" << i << "dprop" << j; 
                std::string graph_name = ss.str();
                ss.str("");
                ss << "ADC Distribution for Drift Radius (" << xmin << "-" << xmax << ")mm and Propagation Distance (" << ymin << "-" << ymax << ")mm; Drift Radius(mm); Propagation Distance(mm); ADC Counts";
                std::string graph_title = ss.str();
                TGraph2DErrors *g= new TGraph2DErrors();
                g->SetNameTitle(graph_name.c_str(), graph_title.c_str());
                for (unsigned int k=0; k<bin_params.size(); k++) { //Loop though bin fits
                    int xbin = ceil(k / num_dprop_bins);
                    int ybin = ceil(k % num_dprop_bins);
                    float x = xmin + ((3.75 / num_dr_bins) * (xbin + 0.5)); 
                    float y = ymin + ((500 / num_dr_bins) * (ybin + 0.5));
                    g->AddPointError(x, y, bin_params.at(k).at(0), 0, 0, bin_params.at(k).at(2));
                }

                TF2 *hfunc = new TF2("histogram bin correction", 
                                    "([0] * TMath::Exp([1]*y)) * ([2] + x*[3] + TMath::Power(x,2)*[4] + TMath::Power(x,3)*[5] + TMath::Power(x,4)*[6])", xmin, xmax, ymin, ymax); 
                TF2 *gfunc = new TF2("graph bin correction", 
                                    "([0] * TMath::Exp([1]*y)) * ([2] + x*[3] + TMath::Power(x,2)*[4] + TMath::Power(x,3)*[5] + TMath::Power(x,4)*[6])", xmin, xmax, ymin, ymax); 
                hfunc->SetParameters(30, -0.00002, 2.2, 1, -0.5, 0.012, -0.0003);
                gfunc->SetParameters(30, -0.00002, 2.2, 1, -0.5, 0.012, -0.0003);

                TCanvas *canvas = new TCanvas(Form("histogram_R%d_dr%d_dprop%d", h, i, j), "", 800, 600);
                binHits[h][i][j]->Fit(hfunc);
                g->Fit(gfunc);
                g->GetXaxis()->SetRangeUser(xmin, xmax);
                g->GetYaxis()->SetRangeUser(ymin, ymax);
                binHits[h][i][j]->Draw("COLZ");
                hfunc->Draw("SAME SURF");
                gfunc->Draw("SAME SURF");



                // Create TPaveText to display the function parameters and their uncertainties
                TPaveText *pt = new TPaveText(0.1, 0.7, 0.3, 0.9, "NDC");
                pt->SetFillColor(0);
                pt->SetTextAlign(12);
                pt->AddText("f(x, y) = (p0*exp(p1*y)) * (p2 + p3*x + p4*x^2 + p5*x^3 + p6*x^4)");
                pt->AddText("Histogram Fit");
                pt->AddText(Form("Chi^2: %f", gfunc->GetChisquare()));
                pt->AddText(Form("NDF: %d", gfunc->GetNDF()));

                for (int p = 0; p < hfunc->GetNpar(); p++) {
                    double param = hfunc->GetParameter(p);
                    double error = hfunc->GetParError(p);
                    pt->AddText(Form("p%d = %.6f ± %.6f", p, param, error));
                }

                pt->AddText("MPV Fit");
                pt->AddText(Form("Chi^2: %f", gfunc->GetChisquare()));
                pt->AddText(Form("NDF: %d", gfunc->GetNDF()));

                for (int p = 0; p < gfunc->GetNpar(); p++) {
                    double param = gfunc->GetParameter(p);
                    double error = gfunc->GetParError(p);
                    pt->AddText(Form("p%d = %.6f ± %.6f", p, param, error));
                }
                pt->Draw("same");

                canvas->Write(); // Write the canvas with the fit drawn and parameters displayed
                binHits[h][i][j]->Write(); // Optionally write the histogram itself

                // Get the fit parameters
                vector<float> params;
                for (int p = 0; p < gfunc->GetNpar(); p++) {
                    params.push_back(gfunc->GetParameter(p));
                }
                binFuncs[h][i][j] = params;

                TCanvas *canvas2 = new TCanvas(Form("GraphFit_R%d_dr%d_dprop%d", h, i, j), "", 800, 600);
                g->Draw();
                gfunc->Draw("same");
                canvas2->Write();

                
                delete g;
                delete hfunc;
                delete gfunc;
                delete pt;
                delete canvas;
                delete canvas2;
                
            }
        }
    }

    /* Now we have all the binned functions. Reiterate through all events keeping track of all hits on track.
    Apply appropriate correction function to each hit and compute truncated mean of each track from corrected
    hits. 
    */ 

    std::vector<float> temp_track;
    std::vector<float> temp_track_stations;
    float track_tm = 0;

    // A histogram for global adc distribution
    TH3F *h3_global_adc_corrected = new TH3F("global adc disytribution", "ADC Distribution after 2D correction; Drift radius(mm); Propagation Distance(mm)", 15,0,15,100,0,5000,160,0,4);

    // Keep the track estimates here
    TH1F *h_Track_TMs = new TH1F("Track Estimator", "40% Truncated mean Track Estimator from 2D Correction (30-40GeV muons); Truncated Mean; Tracks", 200, 0, 4);

    // Now loop through a second time with the corrected ADCs and make tracks/estimators
    for (unsigned int i=0; i<nEntries; i++) {

        // A progress indicator
        if (i % 10000 == 0) {
            std::cout << i << "/" << nEntries << " processed" << endl;
        }
        // cout << "Event " << i << endl;

        // Get the event data
        chain.GetEntry(i);

        // Make sure the current muon has some ADC counts recorded
        if (ADC_counts->size() == 0) {
            continue;
        }

        // cout << "b1" << endl;
        // Mark the current muon link
        int muon_link = 0;
        float muon_pt = pts->at(muon_link);
        float muon_eta = etas->at(muon_link);
        float muon_phi = phis->at(muon_link);
        float muon_d0 = d0s->at(muon_link);
        float muon_z0 = z0s->at(muon_link);
        int muon_author = authors->at(muon_link);
        int muon_charge = charges->at(muon_link);
        float muon_p = ptToP(muon_pt, muon_eta);
        int station_bin;
        int dr_bin;
        int dprop_bin;

        // cout << "b2" << endl;
        for (unsigned int n=0; n<ADC_counts->size(); n++) {

            // cout << "b3" << endl;
            // Check of we have a new muon
            if (trk_link->at(n) != muon_link) {
                // cout << "b4" << endl;
                track_tm = TM(temp_track, 0.4);
                if (!temp_track_stations.empty()) {
                    auto maxIt = std::max_element(temp_track_stations.begin(), temp_track_stations.end());
                    int maxStation = *maxIt;
                    if (maxStation <= 5 && muon_p > 30 && muon_p < 40) {
                        h_Track_TMs->Fill(track_tm);
                    }
                }
                temp_track.clear();
                temp_track_stations.clear();
                // cout << "b5" << endl;

                // Update muon
                muon_link = trk_link->at(n);
                muon_pt = pts->at(muon_link);
                muon_eta = etas->at(muon_link);
                muon_phi = phis->at(muon_link);
                muon_p = ptToP(muon_pt, muon_eta);
                muon_author = authors->at(muon_link);
            }

            // Skip hits not in sector 5 and outlier hits
            if ((station_index->at(n) % 2 )!=0 || station_phi->at(n)==3 || hit_type->at(n) > 60 || std::abs(muon_eta) > 1 || muon_author != 1) {
                continue;
            }
            // CUt on momentum here
            if (muon_p < 30 || muon_p > 40) { continue; }

            // cout << "b6" << endl;
            // Get the length of the current mdt tube
            int chamber_type = chamber_types[station_index->at(n)][abs(station_eta->at(n))][station_phi->at(n)];
            int tube_length = tube_lengths[station_index->at(n)][chamber_type];
            
            // Make sure that the chamber is instantiated
            if (chamber_type == 0) { continue; }

            // Find the location of the tube where the hit occured
            float radial_distance = tube_positions[station_index->at(n)][multi_layer->at(n)][tube_layer->at(n)];
            int sector = stationPhiToSector(station_index->at(n), station_phi->at(n));

            // cout << "b7" << endl;
            float temp_phi; // use this temp phi to rotate all of the tracks as if they are moving vertically for geometry purposes 
            if (muon_phi > 0) {
                temp_phi = muon_phi - (sector-1) *(2*M_PI / 16) + M_PI/2;
            } 
            else if (muon_phi < 0) {
                temp_phi = muon_phi + (2*M_PI - (sector-1) *(2*M_PI / 16)) + M_PI/2;
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
            float a = (radial_distance - y_intercept) / tan(exit_phi); // This is the displacement from the the middle of the MDT along its central axis
            float distance_to_propagate;

            // cout << "b8" << endl;
            // The readout is on opposite sides for large sectors and small sectors
            if (station_index->at(n) % 2 == 0) {
                distance_to_propagate = (tube_length/2) + a;
            }
            else if (station_index->at(n) % 2 == 1) {
                distance_to_propagate = (tube_length/2) - a;
            }

            // cout << "b9" << endl;
            // Retrieve the correction function here
            station_bin = station_index->at(n) / 2;
            dr_bin = drift_radius->at(n) / 3.75;
            dprop_bin = distance_to_propagate / 500;
            std::vector<float> correction_params = binFuncs[station_bin][dr_bin][dprop_bin];
            if (correction_params.size() == 0) { continue; }
            float correction_factor = correctionFunction(ADC_counts->at(n), drift_radius->at(n), distance_to_propagate, correction_params);
            float corrected_adc = ADC_counts->at(n) / correction_factor;
            h3_global_adc_corrected->Fill(drift_radius->at(n), distance_to_propagate, corrected_adc);
            temp_track.push_back(corrected_adc);
            temp_track_stations.push_back(station_index->at(n));
            // cout << "b10" << endl;
        }
        // cout << "b11" << endl;
        // Add the last muon in the event
        auto maxIt = std::max_element(temp_track_stations.begin(), temp_track_stations.end());
        // If the track has no valid hits, skip
        if (temp_track.size()==0) { continue; }
        track_tm = TM(temp_track, 0.6);
        int maxStation = *maxIt;
        if (maxStation <= 5) {
            h_Track_TMs->Fill(track_tm);
        }
        temp_track.clear();
        temp_track_stations.clear();
        // cout << "b12" << endl;
    }
    
    h3_global_adc_corrected->Write();
    h_Track_TMs->Write();

    // Fit a Gaussian function to the histogram
    h_Track_TMs->Fit("gaus");

    // Retrieve the fit function
    TF1 *fitResult = h_Track_TMs->GetFunction("gaus");

    // Extract fit parameters
    double mean = fitResult->GetParameter(1); // Mean of the Gaussian
    double sigma = fitResult->GetParameter(2); // Sigma (standard deviation) of the Gaussian
    double meanError = fitResult->GetParError(1); // Error on the mean
    double sigmaError = fitResult->GetParError(2); // Error on sigma

    // Optionally draw the histogram with the fit
    TCanvas *c1 = new TCanvas("Track Estimators Fit", "Fit Canvas", 800, 600);
    h_Track_TMs->Draw();
    fitResult->Draw("same");
    TPaveText *label = new TPaveText(0.6, 0.7, 0.9, 0.9, "NDC");
    label->AddText(Form("Sigma: %f +/-", sigma, sigmaError));
    label->AddText(Form("Mean: %f +/-", mean, meanError));
    label->AddText(Form("Chisquared: %f", fitResult->GetChisquare()));
    label->AddText(Form("NDF: %d", fitResult->GetNDF()));
    label->Draw("same");
    c1->Write();


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
    calibration2(indir, outfile);
    return 0;
}