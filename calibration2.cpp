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
float correctionFunction(float drift_radius, float propagation_distance, std::vector<float> params) {
    return (params.at(0) * std::exp(params.at(1)*propagation_distance)) * 
    (params.at(2) + drift_radius*params.at(3) +
     std::pow(drift_radius, 2)*params.at(4) + 
     std::pow(drift_radius, 3)*params.at(5) + 
     std::pow(drift_radius, 4)*params.at(6));
}

// A function to create the x% truncated mean from a vector of ADC counts
float TM(std::vector<float> adcs, float truncation_percent) {
    float size = adcs.size();
    float truncated_size = ceil(size * (1-truncation_percent));
    
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

// Return the truncated mean of a 1D histogram and the error
std::pair<float, float> fit1DHistogramTM(TH1 *h1, const std::string &option="") {
    int nbins = h1->GetNbinsX();
    double hist_entries = h1->GetEntries();
    if (hist_entries==0) { return std::make_pair(0.0, 0.0); }
    int cutoff = std::floor(hist_entries * 0.6);
    int entry_count = 0;
    float tm;
    float last_entry;
    float remaining = 0;
    
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
            remaining = cutoff - entry_count;
            tm += binCenter * remaining;
            entry_count = cutoff;
            last_entry = binCenter * remaining;
            break;
        }
    }

    if (option=="W") { // return also the standard error from a winsorized mean
        float winsorized_mean = (tm + (hist_entries - cutoff) * (last_entry)) / hist_entries;
        float winsorized_square_deviation = 0;
        entry_count = 0;
        for (unsigned int i=1; i<=nbins; i++) { // Loop over histogram entries
            float binContent = h1->GetBinContent(i);
            float binLowEdge = h1->GetBinLowEdge(i);
            float binWidth = h1->GetBinWidth(i);
            float binCenter = binLowEdge + 0.5 * binWidth;
            if (entry_count + binContent < cutoff) {
                winsorized_square_deviation += pow((binCenter * binContent) - winsorized_mean, 2);
                entry_count += binContent;
            }
            else {
                winsorized_square_deviation += remaining * pow((last_entry - winsorized_mean), 2);
                break;
            }
        }
        float winsorized_mean_standard_error = winsorized_square_deviation / TMath::Sqrt((hist_entries - cutoff) * (hist_entries - cutoff - 1));
        tm /= cutoff;
        return std::make_pair(tm, winsorized_mean_standard_error);
    } else {
        tm /= cutoff;
        return std::make_pair(tm, 0.0); 
    }
}

/*
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
*/

// Return a Truncated mean from each xBin in a TH2 object
std::vector<std::vector<float>> fit2DHistogramTM(TH2 *h2, const std::string &option="") {
    int nbinsx = h2->GetNbinsX();
    const char* h2_title = h2->GetTitle();

    // Make a figure where all the bin porjections are plotted with the fit on top along with labels.

    // Record the parameters of the Landau fits <MPV, sigma, MPV_err, sigma_err>
    std::vector<std::vector<float>> landau_params;
    //Loop through the bins
    for (int i = 1; i <= nbinsx; i++) {
        // Turn each xbin into a 1d histogram
        TH1D* projectionY = h2->ProjectionY("", i, i);
        std::pair<float, float> tm = fit1DHistogramTM(projectionY, "w");
        int hist_entries = projectionY->GetEntries();
        if (hist_entries==0) {
            std::vector<float> results = {tm.first, 0.0, tm.second, 0.0};
            landau_params.push_back(results);
            continue;
        }


        // Fit a landau to the xbin histogram
        projectionY->Fit("landau");

        // Retrieve the fit result
        TF1* fitResult = projectionY->GetFunction("landau");
        float mpv = fitResult->GetParameter(1); // Most probable value (MPV)
        float sigma = fitResult->GetParameter(2); // Sigma (width of the distribution)
        float mpvError = fitResult->GetParError(1); // Error on MPV
        float sigmaError = fitResult->GetParError(2); // Error on Sigma
        std::vector<float> results = {tm.first, 0.0, tm.second, 0.0};
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
std::vector<std::vector<float>> fit3DHistogramTM(TH3 *h3, const std::string &option="") {
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
            int hist_entries = projectionZ->GetEntries();
            if (hist_entries == 0) { // Make sure the TH1 is not empty
                landau_params.push_back({0.0, 0.0, -10000.0, -10000.0});
                continue;
            }
            std::pair<float, float> binTM = fit1DHistogramTM(projectionZ);
            float tm = binTM.first;
            float tm_error = binTM.second;
            float mpv;
            float sigma;
            float mpvError;
            float sigmaError;

            // Fit a landau to the xbin histogram
            projectionZ->Fit("landau");

            // Retrieve the fit result
            TF1* fitResult = projectionZ->GetFunction("landau");
            mpv = fitResult->GetParameter(1); // Most probable value (MPV)
            sigma = fitResult->GetParameter(2); // Sigma (width of the distribution)
            mpvError = fitResult->GetParError(1); // Error on MPV
            sigmaError = fitResult->GetParError(2); // Error on Sigma
            std::vector<float> results = {tm, sigma, tm_error, sigmaError};
            landau_params.push_back(results);

            if (option == "DRAW") { // Create a canvas on which to draw the fit
                TCanvas *canvas = new TCanvas(Form("Momentum bin x%d|y%d from %s", i, j, h_title), Form("Momentum Bin x%d|y%d Projection and Landau Fit from %s", i, j, h_title), 800, 600);
                TPaveText *label = new TPaveText(0.6, 0.7, 0.9, 0.9, "NDC");
                label->AddText(Form("Chi^2: %f", fitResult->GetChisquare()));
                label->AddText(Form("NDF: %d", fitResult->GetNDF()));
                label->AddText(Form("MPV: %f +/- %.2f", mpv, mpvError));
                label->AddText(Form("Scale: %f +/- %.2f", sigma, sigmaError));
                label->AddText(Form("40% Truncated Mean: %f +/- %f", tm, tm_error));
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

// A function to convert a TH3 to a vector of params for the correction function fit to the TH3
std::vector<float> processTMs(TH3* histogram, const std::string& name) {
    std::vector<std::vector<float>> tms = fit3DHistogramTM(histogram, "DRAW");
    if (tms.empty()) {
        std::cout << name << "_tms is empty!" << std::endl;
        return {0.0,0.0,0.0,0.0,0.0,0.0};
    }
    
    TGraph2DErrors* g_tms = new TGraph2DErrors();

    int pointsAdded = 0;
    for (unsigned int i = 0; i < tms.size(); i++) {
        int xbin = (i < 10) ? 0 : i / 10;
        int ybin = i % 10;
        float x = 0.75 + xbin * 1.5;
        float y = 250.0 + ybin * 500.0;
        if (tms[i][0] == 0.0) { continue; } // Skip the point if the TM is identically zero
        g_tms->AddPointError(x, y, tms[i][0], 0.0, 0.0, tms[i][0]);
        std::cout << "Added TM " << tms[i][0] << " at point " << x << "," << y << std::endl;
        pointsAdded++;
    }

    std::cout << "Total points added: " << pointsAdded << std::endl;
    std::cout << "Number of points in graph: " << g_tms->GetN() << std::endl;

    TF2* func = new TF2(("graph " + name + " correction").c_str(), 
                        "([0] * TMath::Exp([1]*y)) * ([2] + x*[3] + TMath::Power(x,2)*[4] + TMath::Power(x,3)*[5] + TMath::Power(x,4)*[6])", 0, 15, 0, 5000);
    func->SetParameters(140, -0.00002, 0.58, 0.29, -0.05, 0.0035, -0.0001);
    
    TCanvas* canvas = new TCanvas();
    g_tms->Draw();
    g_tms->Fit(func);
    func->Draw("same SURF");

    TPaveText* pt = new TPaveText(0.1, 0.7, 0.3, 0.9, "NDC");
    pt->SetFillColor(0);
    pt->SetTextAlign(12);
    pt->AddText("f(x, y) = (p0*exp(p1*y)) * (p2 + p3*x + p4*x^2 + p5*x^3 + p6*x^4)");
    pt->AddText(Form("Chi^2: %f", func->GetChisquare()));
    pt->AddText(Form("NDF: %d", func->GetNDF()));

    std::vector<float> params;
    for (int p = 0; p < func->GetNpar(); p++) {
        double param = func->GetParameter(p);
        double error = func->GetParError(p);
        params.push_back(param);
        pt->AddText(Form("p%d = %.6f +/- %.6f", p, param, error));
    }
    pt->Draw("same");

    canvas->Write();
    g_tms->Write();
    
    delete canvas;
    delete g_tms;
    delete func;
    delete pt;
    return params;
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

// Check that the track only has hits in preferred regions
int trackRegionCheck(int sector, int minEta, int maxEta, std::vector<int> stations, std::vector<int> etas) {
    if (sector % 2 == 1) { // Make sure the track is in a large sector
        for (unsigned int i=0; i<stations.size(); i++) {
            if (stations[i] % 2 == 1) {
                return 0;
            }
        }
    } else if (sector % 2 == 0) { // Make sure the track is in a small sector
        for (unsigned int i=0; i<stations.size(); i++) {
            if (stations[i] % 2 == 0) {
                return 0;
            }
        }
    }
    for (unsigned int i=0; i<etas.size(); i++) { // Make sure that all the hits are in the allowed eta range
        if (etas[i] < minEta || etas[i] > maxEta) {
            return 0;
        }
    }
    return 1;
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
    chain.SetBranchStatus("runNumber", true);



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
    Int_t runNumber;

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
    chain.SetBranchAddress("runNumber", &runNumber);



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

    // Create a single hsitogram to which we will extract MPVs from bins to fit a single global ADC correction
    TH3F *h3_global_adc_raw = new TH3F("global adc distribution before correcction", "ADC Distribution BEFORE 2D correction; Drift radius(mm); Propagation Distance(mm)", 10,0,15,10,0,5000,100,0,400);

    // Create a histogram fro each of the three stations:
    TH3F *h3_inner = new TH3F("ADC distribution in BIL", "ADC Distribution BEFORE 2D correction in BIL|Sector 5|Eta 1-3; Drift Radius(mm); Propagation Distance(mm); ADC Counts", 10,0,15,10,0,5000,100,0,400);
    TH3F *h3_middle = new TH3F("ADC distribution in BML", "ADC Distribution BEFORE 2D correction in BML|Sector 5|Eta 1-3; Drift Radius(mm); Propagation Distance(mm); ADC Counts",10,0,15,10,0,5000,100,0,400);
    TH3F *h3_outer = new TH3F("ADC distribution in BOL", "ADC Distribution BEFORE 2D correction in BOL|Sector 5|Eta 1-3; Drift Radius(mm); Propagation Distance(mm); ADC Counts", 10,0,15,10,0,5000,100,0,400);

    // Histograms for tube adc in each chamber
    //BIL has 240 or 288 (36*2*4=240 for BIL sector 5 Eta=1)
    TH2F *h2_BIL_tubes = new TH2F("ADC Distribution By Tube - BIL", "ADC Distribution by tube in BIL Sector 5, Eta=1 Chamber; Tube; ADC Counts", 288,1,289, 400,0,400);
    //BML has 288 or 336 (56*2*3=288 for BML sector 5 Eta=1)
    TH2F *h2_BML_tubes = new TH2F("ADC Distribution By Tube - BML", "ADC Distribution by tube in BML Sector 5, Eta=1 Chamber; Tube; ADC Counts", 336,1,337, 400,0,400);
    //BOL has 336 or 432 (72*2*3=384 for BOL sector 5 Eta=1)
    TH2F *h2_BOL_tubes = new TH2F("ADC Distribution By Tube - BOL", "ADC Distribution by tube in BOL Sector 5, Eta=1 Chamber; Tube; ADC Counts", 432,1,433, 400,0,400);

    // Make Histograms for adc vs. drift radius in every propagation distance bin in BOL Type 1 chamber 500mm at a time
    // and adc vs dprop in every drift radius in BOL Type 1 chamber 1.5mm at a time
    const int NUM_DR_BINS = 10;  // 0-15mm in 1.5mm bins
    const double DR_BIN_WIDTH = 1.5;
    const int NUM_DPROP_BINS = 10;  // 0-5000mm in 500mm bins
    const double DPROP_BIN_WIDTH = 500.0;

    std::vector<TH2F*> h2_dprop_dr_vec_raw;
    std::vector<TH2F*> h2_dr_dprop_vec_raw;
    std::vector<TH2F*> h2_dprop_dr_vec_corrected;
    std::vector<TH2F*> h2_dr_dprop_vec_corrected;

    // Create histograms for ADC vs. propagation distance for each drift radius bin
    for (int i = 0; i < NUM_DR_BINS; i++) {
        double dr_min = i * DR_BIN_WIDTH;
        double dr_max = (i + 1) * DR_BIN_WIDTH;
        
        // Raw histogram
        std::string hist_name_raw = "h2_dprop_dr_" + std::to_string(i) + "_raw";
        std::string hist_title_raw = "ADC vs. dprop BOL Type 1 Raw, Drift Radius (" + 
                                    std::to_string(dr_min) + "-" + std::to_string(dr_max) + "mm)";
        
        TH2F *h2_raw = new TH2F(hist_name_raw.c_str(), 
                                (hist_title_raw + "; Propagation Distance(mm); ADC Counts").c_str(), 
                                20, 0, 5000, 200, 0, 400);
        
        h2_dprop_dr_vec_raw.push_back(h2_raw);

        // Corrected histogram
        std::string hist_name_corrected = "h2_dprop_dr_" + std::to_string(i) + "_corrected";
        std::string hist_title_corrected = "ADC vs. dprop BOL Type 1 Corrected, Drift Radius (" + 
                                        std::to_string(dr_min) + "-" + std::to_string(dr_max) + "mm)";
        
        TH2F *h2_corrected = new TH2F(hist_name_corrected.c_str(), 
                                    (hist_title_corrected + "; Propagation Distance(mm); ADC Counts").c_str(), 
                                    20, 0, 5000, 200, 0, 4);
        
        h2_dprop_dr_vec_corrected.push_back(h2_corrected);
    }

    // Create histograms for ADC vs. drift radius for each propagation distance bin
    for (int i = 0; i < NUM_DPROP_BINS; i++) {
        double dprop_min = i * DPROP_BIN_WIDTH;
        double dprop_max = (i + 1) * DPROP_BIN_WIDTH;
        
        // Raw histogram
        std::string hist_name_raw = "h2_dr_dprop_" + std::to_string(i) + "_raw";
        std::string hist_title_raw = "ADC vs. dr BOL Type 1 Raw, Propagation Distance (" + 
                                    std::to_string(dprop_min) + "-" + std::to_string(dprop_max) + "mm)";
        
        TH2F *h2_raw = new TH2F(hist_name_raw.c_str(), 
                                (hist_title_raw + "; Drift Radius(mm); ADC Counts").c_str(), 
                                10, 0, 15, 200, 0, 400);
        
        h2_dr_dprop_vec_raw.push_back(h2_raw);

        // Corrected histogram
        std::string hist_name_corrected = "h2_dr_dprop_" + std::to_string(i) + "_corrected";
        std::string hist_title_corrected = "ADC vs. dr BOL Type 1 Corrected, Propagation Distance (" + 
                                        std::to_string(dprop_min) + "-" + std::to_string(dprop_max) + "mm)";
        
        TH2F *h2_corrected = new TH2F(hist_name_corrected.c_str(), 
                                    (hist_title_corrected + "; Drift Radius(mm); ADC Counts").c_str(), 
                                    10, 0, 15, 200, 0, 4);
        
        h2_dr_dprop_vec_corrected.push_back(h2_corrected);
    }



    // Vectors here of hits in each station for single chamber TM
    std::vector<float> BIL_hits;
    std::vector<float> BML_hits;
    std::vector<float> BOL_hits;

    /***
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

        // Restrict to a single run
        // if (runNumber != 473747) { continue; }

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

        // SKip empty muons
        if (ADC_counts->size()==0) { continue; }

        // Loop through the muon hits
        for (unsigned int n=0; n<ADC_counts->size(); n++) {
            

            // Skip hits not in sector 5, |eta|>1 and outlier hits
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
            } else if (muon_phi < 0) {
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
            
            // Place the hits into the right histograms
            station_bin = station_index->at(n) / 2;
            dr_bin = drift_radius->at(n) / 3.75;
            dprop_bin = distance_to_propagate / 500;
            if (dprop_bin > 9 || dprop_bin < 0 || dr_bin > 3 || dr_bin < 0) { continue; }
            binHits[station_bin][dr_bin][dprop_bin]->Fill(drift_radius->at(n), distance_to_propagate, ADC_counts->at(n));
            h3_global_adc_raw->Fill(drift_radius->at(n), distance_to_propagate, ADC_counts->at(n));

            int hist_dr_bin = drift_radius->at(n) / 1.5;
            int hist_dprop_bin = distance_to_propagate / 500.0;
            if (station_index->at(n)==4 && distance_to_propagate < tube_length && distance_to_propagate > 0) {
                h2_dprop_dr_vec_raw.at(hist_dr_bin)->Fill(distance_to_propagate, ADC_counts->at(n));
                h2_dr_dprop_vec_raw.at(hist_dprop_bin)->Fill(drift_radius->at(n), ADC_counts->at(n));
            }

            // For stationed corrections restrict the chambers which are accepted
            if (station_eta->at(n) == 1) {
                float tube_num;
                if (station_index->at(n) == 0) {
                    h3_inner->Fill(drift_radius->at(n), distance_to_propagate, ADC_counts->at(n));
                    tube_num = (multi_layer->at(n) - 1) * (4*36) + (tube_layer->at(n) - 1) * 30 + tubes->at(n);
                    h2_BIL_tubes->Fill(tube_num, ADC_counts->at(n));
                    BIL_hits.push_back(ADC_counts->at(n));
                }
                else if (station_index->at(n) == 2) {
                    h3_middle->Fill(drift_radius->at(n), distance_to_propagate, ADC_counts->at(n));
                    tube_num = (multi_layer->at(n) - 1) * (3*56) + (tube_layer->at(n) - 1) * 48 + tubes->at(n);
                    h2_BML_tubes->Fill(tube_num, ADC_counts->at(n));
                    BML_hits.push_back(ADC_counts->at(n));
                }
                else if (station_index->at(n) == 4) {
                    h3_outer->Fill(drift_radius->at(n), distance_to_propagate, ADC_counts->at(n));
                    tube_num = (multi_layer->at(n) - 1) * (3*72) + (tube_layer->at(n) - 1) * 72 + tubes->at(n);
                    h2_BOL_tubes->Fill(tube_num, ADC_counts->at(n));
                    BOL_hits.push_back(ADC_counts->at(n));
                }
            }
        }
    }


    // Create the output file
    TFile *output_file = new TFile(outfile.c_str(), "RECREATE");

    /*
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
                std::vector<std::vector<float>> bin_params = fit3DHistogramTM(binHits[h][i][j], "DRAW");
                
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
    */

    // Draw the drift radius distribution in poropagtaion distance bins and vice versa
    for (unsigned int i=0; i<h2_dprop_dr_vec_raw.size(); i++) {
        h2_dprop_dr_vec_raw.at(i)->Write();
    }
    for (unsigned int i=0; i<h2_dr_dprop_vec_raw.size(); i++) {
        h2_dr_dprop_vec_raw.at(i)->Write();
    }

    // One single global correction here
    std::vector<std::vector<float>> global_tms = fit3DHistogramTM(h3_global_adc_raw, "DRAW");
    if (global_tms.empty()) {
        std::cout << "global_tms is empty!" << std::endl;
    }
    TGraph2DErrors *g_global_tms = new TGraph2DErrors();

    int pointsAdded = 0;
    for (unsigned int i=0; i<global_tms.size(); i++) {
        int xbin;
        if (i < 10) {
            xbin = 0;
        } else {
            xbin = i / 10;
        }
    
        int ybin = i % 10;
        float x = 0.75 + xbin * (1.5);
        float y = 250.0 + ybin * (500.0);
        g_global_tms->AddPointError(x, y, global_tms.at(i).at(0), 0.0, 0.0, global_tms.at(i).at(0));
        cout << "Added TM " << global_tms.at(i).at(0) << " at point " << x << "," << y << endl;
        pointsAdded++;
    }

    cout << "Total points added: " << pointsAdded << endl;
    cout << "Number of points in graph: " << g_global_tms->GetN() << endl;

    TF2 *gfunc = new TF2("graph global correction", 
                        "([0] * TMath::Exp([1]*y)) * ([2] + x*[3] + TMath::Power(x,2)*[4] + TMath::Power(x,3)*[5] + TMath::Power(x,4)*[6])", 0, 15, 0, 5000); 
    gfunc->SetParameters(140, -0.00002, 0.58, 0.29, -0.05, 0.0035, -0.0001);
    TCanvas *canvas = new TCanvas();
    g_global_tms->Draw();
    g_global_tms->Fit(gfunc);
    gfunc->Draw("same SURF");

    // Create TPaveText to display the function parameters and their uncertainties
    TPaveText *pt = new TPaveText(0.1, 0.7, 0.3, 0.9, "NDC");
    pt->SetFillColor(0);
    pt->SetTextAlign(12);
    pt->AddText("f(x, y) = (p0*exp(p1*y)) * (p2 + p3*x + p4*x^2 + p5*x^3 + p6*x^4)");
    pt->AddText(Form("Chi^2: %f", gfunc->GetChisquare()));
    pt->AddText(Form("NDF: %d", gfunc->GetNDF()));

    for (int p = 0; p < gfunc->GetNpar(); p++) {
        double param = gfunc->GetParameter(p);
        double error = gfunc->GetParError(p);
        pt->AddText(Form("p%d = %.6f ± %.6f", p, param, error));
    }
    pt->Draw("same");

    canvas->Write();
    g_global_tms->Write();
    delete canvas;
    delete g_global_tms;
    delete gfunc;
    delete pt;

    // Now fit a function for every station
    // BIL first
    // Main function or where you want to call the processing
    std::vector<float> BIL_correction = processTMs(h3_inner, "inner");
    std::vector<float> BML_correction = processTMs(h3_middle, "middle");
    std::vector<float> BOL_correction = processTMs(h3_outer, "outer");

    cout << "\n\nStation corrections obtained\n\n " << endl;

    cout << "b1" << endl;
    // Get the tube level MPVs here process
    std::vector<std::vector<float>> BIL_tube_landaus = fit2DHistogramTM(h2_BIL_tubes);
    std::vector<float> BIL_tube_tms;
    for (unsigned int i=0; i<BIL_tube_landaus.size(); i++) {
        cout << "BIL pushing back tube " << i << endl;
        BIL_tube_tms.push_back(BIL_tube_landaus.at(i).at(0));
    }
    std::vector<std::vector<float>> BML_tube_landaus = fit2DHistogramTM(h2_BML_tubes);
    std::vector<float> BML_tube_tms;
    for (unsigned int i=0; i<BML_tube_landaus.size(); i++) {
        cout << "BML pushing back tube " << i << endl;
        BML_tube_tms.push_back(BML_tube_landaus.at(i).at(0));
    }
    std::vector<std::vector<float>> BOL_tube_landaus = fit2DHistogramTM(h2_BOL_tubes);
    std::vector<float> BOL_tube_tms;
    for (unsigned int i=0; i<BOL_tube_landaus.size(); i++) {
        cout << "BOL pushing back tube " << i << endl;
        BOL_tube_tms.push_back(BOL_tube_landaus.at(i).at(0));
    }
    float BIL_total_tm = TM(BIL_hits, 0.4);
    float BML_total_tm = TM(BML_hits, 0.4);
    float BOL_total_tm = TM(BOL_hits, 0.4);

    cout << "b2" << endl;

    // Get the global correction parameters
    // std::vector<float> global_correction_params;
    /*
    for (int p = 0; p < gfunc->GetNpar(); p++) {
        global_correction_params.push_back(gfunc->GetParameter(p));
        cout << "p" << p << " is " << gfunc->GetParameter(p) << endl;
    }
    */

    /*
    Now we have all the binned functions. Reiterate through all events keeping track of all hits on track.
    Apply appropriate correction function to each hit and compute truncated mean of each track from corrected
    hits. 
    */ 
    std::vector<float> temp_track;
    std::vector<int> temp_track_stations;
    std::vector<int> temp_track_etas;
    float track_tm = 0;

    // A histogram for global adc distribution
    TH3F *h3_global_adc_corrected = new TH3F("global adc distribution after correction", "ADC Distribution AFTER 2D correction; Drift radius(mm); Propagation Distance(mm)", 15,-15,15,50,-5000,5000,160,-100,100);

    // Keep the track estimates here
    TH1F *h_Track_TMs = new TH1F("Track Estimator", "40% Truncated mean Track Estimator from 2D Correction (60-80GeV muons); Truncated Mean; Tracks", 120, 0.6, 1.4);

    // Now loop through a second time with the corrected ADCs and make tracks/estimators
    for (unsigned int i=0; i<nEntries; i++) {

        // A progress indicator
        if (i % 10000 == 0) {
            std::cout << i << "/" << nEntries << " processed" << endl;
        }

        // Get the event data
        chain.GetEntry(i);

        // Restrict to a single run
        // if (runNumber != 473747) { continue; }

        // Make sure the current muon has some ADC counts recorded
        if (ADC_counts->size() == 0) { continue; }

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

        for (unsigned int n=0; n<ADC_counts->size(); n++) {

            // Check if we have a new muon
            if (trk_link->at(n) != muon_link) {
                
                // Make sure the track isn't empty and fill the old muon
                if (!temp_track_stations.empty()) {
                    track_tm = TM(temp_track, 0.4);
                    if (trackRegionCheck(5, 1, 1, temp_track_stations, temp_track_etas) == 1) {
                        h_Track_TMs->Fill(track_tm);
                    }
                }
                temp_track.clear();
                temp_track_stations.clear();
                temp_track_etas.clear();
                
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

            // Cut on momentum here
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

            // Retrieve the correction function here
            /*
            station_bin = station_index->at(n) / 2;
            dr_bin = drift_radius->at(n) / 3.75;
            dprop_bin = distance_to_propagate / 500;
            std::vector<float> correction_params = binFuncs[station_bin][dr_bin][dprop_bin];
            if (correction_params.size() == 0) { continue; }
            float correction_factor = correctionFunction(ADC_counts->at(n), drift_radius->at(n), distance_to_propagate, correction_params);
            float corrected_adc = ADC_counts->at(n) / correction_factor;
            h3_global_adc_corrected->Fill(drift_radius->at(n), distance_to_propagate, corrected_adc);
            */

            // Below is the global correction factor
            // float correction_factor = correctionFunction(drift_radius->at(n), distance_to_propagate, global_correction_params);

            // Here is the correction factor for stationed correction, along with the corrective tube factor
            float correction_factor = 0;
            float tube_num;
            float tube_factor;
            if (station_index->at(n) == 0) {
                correction_factor = correctionFunction(drift_radius->at(n), distance_to_propagate, BIL_correction);
                tube_num = (multi_layer->at(n) - 1) * (4*36) + (tube_layer->at(n) - 1) * 36 + tubes->at(n);
                tube_factor = BIL_tube_tms.at(tube_num - 1) / BIL_total_tm;
            } else if (station_index->at(n) == 2) {
                correction_factor = correctionFunction(drift_radius->at(n), distance_to_propagate, BML_correction);
                tube_num = (multi_layer->at(n) - 1) * (3*56) + (tube_layer->at(n) - 1) * 56 + tubes->at(n);
                tube_factor = BML_tube_tms.at(tube_num - 1) / BML_total_tm;
            } else if (station_index->at(n) == 4) {
                correction_factor = correctionFunction(drift_radius->at(n), distance_to_propagate, BOL_correction);
                tube_num = (multi_layer->at(n) - 1) * (3*72) + (tube_layer->at(n) - 1) * 72 + tubes->at(n);
                tube_factor = BOL_tube_tms.at(tube_num - 1) / BOL_total_tm;
            }

            float corrected_adc = ADC_counts->at(n) * (tube_factor / correction_factor);
            // float corrected_adc = ADC_counts->at(n) / correction_factor;

            // cout << "ADC is " << ADC_counts->at(n) << " correction factor is " << correction_factor << " for drift radius " << drift_radius->at(n) << " and propagation distance " << distance_to_propagate << " and corrected ADC is " << corrected_adc << endl;
            h3_global_adc_corrected->Fill(drift_radius->at(n), distance_to_propagate, corrected_adc);
            temp_track.push_back(corrected_adc);
            temp_track_stations.push_back(station_index->at(n));
            temp_track_etas.push_back(station_eta->at(n));

            // FIll the corrected adc histograms here
            int hist_dr_bin = drift_radius->at(n) / 1.5;
            int hist_dprop_bin = distance_to_propagate / 500.0;
            if (station_index->at(n)==4 && distance_to_propagate < tube_length && distance_to_propagate > 0) {
                h2_dprop_dr_vec_corrected.at(hist_dr_bin)->Fill(distance_to_propagate, corrected_adc);
                h2_dr_dprop_vec_corrected.at(hist_dprop_bin)->Fill(drift_radius->at(n), corrected_adc);
            }
        }

        // Add the last muon in the event
        // Make sure the track isn't empty
        if (!temp_track_stations.empty()) {
            track_tm = TM(temp_track, 0.4);
            // Ensure the track is in the barrel
            if (trackRegionCheck(5, 1, 1, temp_track_stations, temp_track_etas) == 1) {
                h_Track_TMs->Fill(track_tm);
            }
        }
        temp_track.clear();
        temp_track_stations.clear();
        temp_track_etas.clear();
    }
    cout << "2nd pass finished" << endl;

    // Draw the drift radius Distribution in poropagtaion distance bins and vice versa after corrections
    for (unsigned int i=0; i<h2_dprop_dr_vec_corrected.size(); i++) {
        h2_dprop_dr_vec_corrected.at(i)->Write();
    }
    for (unsigned int i=0; i<h2_dr_dprop_vec_corrected.size(); i++) {
        h2_dr_dprop_vec_corrected.at(i)->Write();
    }
    
    h3_inner->Write();
    h3_middle->Write();
    h3_outer->Write();
    h3_global_adc_corrected->Write();
    h2_BIL_tubes->Write();
    h2_BML_tubes->Write();
    h2_BOL_tubes->Write();
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
    label->AddText(Form("Sigma: %f +/- %f", sigma, sigmaError));
    label->AddText(Form("Mean: %f +/- %f", mean, meanError));
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