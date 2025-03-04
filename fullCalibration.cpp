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
#include <TLatex.h>
#include <TLegend.h>
using namespace std;


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
    if (truncated_mean==0) {
        cout << "Track sum is " << track_sum << " and track size is " << truncated_size << endl;
    }
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

// Return the Landau MPV and error of a 1D histogram
std::pair<float, float> fit1DHistogramLandau(TH1 *h1, const std::string &option="") {
    
    // Record the parameters of the Landau fits <MPV, sigma, MPV_err, sigma_err>
    std::pair<float, float> landau_params;
    //Loop through the bins
    h1->Fit("landau");

    // Retrieve the fit result
    TF1* fitResult = h1->GetFunction("landau");
    if (!fitResult) {
        std::cerr << "Error: Landau fit failed or function not found!" << std::endl;
        return std::make_pair(-9999.0, -9999.0);  // Return an error value
    }
    float mpv = fitResult->GetParameter(1); // Most probable value (MPV)
    float mpvError = fitResult->GetParError(1); // Error on MPV
    std::pair<float, float> results = std::make_pair(mpv, mpvError);
    return results;
}

// Return vector with <mpv, scake, mpvErr, scaleErr> from each xBin in a TH2 object
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

// Return a Truncated mean from each xyBin in a TH3 object
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
            float sigma = 0;
            float sigmaError = 0;

            std::vector<float> results = {tm, sigma, tm_error, sigmaError};
            landau_params.push_back(results);
        }
    }
    return landau_params;
}

// A function to convert a TH3 to a vector of params for the correction function fit to the TH3
std::vector<float> processTMs(TH3* histogram, const std::string& name) {
    std::vector<std::vector<float>> tms = fit3DHistogramTM(histogram);
    double xMax = histogram->GetXaxis()->GetXmax();
    double yMax = histogram->GetYaxis()->GetXmax();
    cout << "xmax is " << xMax << " and ymax is " << yMax << endl;
    if (tms.empty()) {
        std::cout << name << "_tms is empty!" << std::endl;
        return {0.0,0.0,0.0,0.0,0.0,0.0};
    }
    
    TGraph2DErrors* g_tms = new TGraph2DErrors();
    for (unsigned int i = 0; i < tms.size(); i++) {
        int xbin = (i < 10) ? 0 : i / 10;
        int ybin = i % 10;
        float xbin_size = xMax / 10.0;
        float ybin_size = yMax / 10.0;
        float x = (xbin_size * 0.5) + (xbin * xbin_size);
        float y = (ybin_size * 0.5) + (ybin * ybin_size);
        if (tms[i][0] == 0.0) { continue; } // Skip the point if the TM is identically zero
        g_tms->AddPointError(x, y, tms[i][0], 0.0, 0.0, tms[i][0]);
    }

    TF2* func = new TF2(("graph " + name + " correction").c_str(), 
                        "([0] * TMath::Exp([1]*y)) * ([2] + x*[3] + TMath::Power(x,2)*[4] + TMath::Power(x,3)*[5] + TMath::Power(x,4)*[6])", 0, 15, 0, 5000);
    func->SetParameters(140, -0.00002, 0.58, 0.29, -0.05, 0.0035, -0.0001);
    
    TCanvas* canvas = new TCanvas();
    g_tms->Draw();
    g_tms->Fit(func, "Q");
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

// A function to evaluate the 2D correction function for a given ADC
float correctionFunction(float drift_radius, float propagation_distance, std::vector<float> params) {
    return (params.at(0) * std::exp(params.at(1)*propagation_distance)) * 
    (params.at(2) + drift_radius*params.at(3) +
     std::pow(drift_radius, 2)*params.at(4) + 
     std::pow(drift_radius, 3)*params.at(5) + 
     std::pow(drift_radius, 4)*params.at(6));
}

// A function to calculate the invariant mass of two particles created in a collision
float invariantMass(float pt1, float pt2, float eta1, float eta2, float phi1, float phi2) {
    return sqrt(2 * pt1 * pt2 * (cosh(eta1 - eta2) - cos(phi1 - phi2)));
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

// A function to find which momentum bin the track goes in
int getMomentumBin(float p) {
    if (p < 20) {
    return -1;
    } else if (p < 260) {
        return ((p - 20) / 10);
    } else {
        return -1;
    }
}

// The function which we calibrate in every bin
Double_t fitFunc(Double_t *vals, Double_t *pars) {
    Float_t x = vals[0];
    Float_t y = vals[1]; 
    Double_t f = (pars[0] * TMath::Exp(pars[1]*x) + pars[2]) * (pars[3] + y*pars[4] + TMath::Power(y,2)*pars[5] + TMath::Power(y,3)*pars[6] + TMath::Power(y,4)*pars[7]);
    return f;
}

// DEFINE CONSTANTS HERE
const double zmass = 91.0;
const double zmass_window = 10.0;

void tubeCalibration(std::string indir, const::std::string& outfile) {

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

    
    // Create the output file
    TFile *output_file = new TFile(outfile.c_str(), "RECREATE");

    // Loop trhough events and create track estimates without any corrections
    std::vector<float> temp_track;
    std::vector<int> temp_track_stations;
    std::vector<int> temp_track_etas;
    float track_tm = 0;

    // Make histograms to show the adc distribution vs. drift radius, propagation distance and momentum
    TH2F *h2_dr_all = new TH2F("h2_dr_all", "ADC Distribution vs. Drift Radius; Drift Radius(mm); ADC Counts", 15, 0, 15, 400, 0, 400);
    TH2F *h2_dprop_all = new TH2F("h2_dprop_all", "ADC Distribution vs. Propagation Distance; Propagation Distance(mm); ADC Counts", 50, 0, 5000, 400, 0, 400);
    TH2F *h2_momentum_all = new TH2F("h2_momentum_all", "ADC Distribution vs. Momentum; Momentum (GeV); ADC Counts", 100, 0, 200, 400, 0, 400);
    TH1F *h_tube_tms = new TH1F("h_tube_tms", "ADC Count Truncated Means Across All Tubes; 40% Truncated Mean; Tubes / bin", 200, 0, 200);
    TH1F *h_chamber_tms = new TH1F("h_chamber_tms", "ADC Count Truncated Means Across All Chambers; 40% Truncated Mean; Chambers / bin", 100, 0, 2);
    
    // A plot of uncorrected ADC distribution
    TH1F *h_uncorrected_ADC = new TH1F("Uncorrected ADC", "Global ADC Distribution Before Any Correction", 400, 0, 400);

    // Keep the track estimates here
    TH1F *h_Track_TMs = new TH1F("Track Estimator", "40% Truncated mean Track Estimator TUBE+2D CORRECTION (All Momenta); Truncated Mean; Tracks(AU)", 400, 0, 4);
    TH1F *h_Track_TMs_1 = new TH1F("0-14 Hits", "40% Truncated mean Track Estimator TUBE+2D CORRECTION (All Momenta) 0-14 Hits; Truncated Mean; Tracks(AU)", 400, 0, 4);
    TH1F *h_Track_TMs_2 = new TH1F("15-29 Hits", "40% Truncated mean Track Estimator TUBE+2D CORRECTION (All Momenta) 15-29 Hits; Truncated Mean; Tracks(AU)", 400, 0, 4);
    TH1F *h_Track_TMs_3 = new TH1F("30-39 Hits", "40% Truncated mean Track Estimator TUBE+2D CORRECTION (All Momenta) 30-39 Hits; Truncated Mean; Tracks(AU)", 400, 0, 4);
    TH1F *h_Track_TMs_4 = new TH1F("40+ Hits", "40% Truncated mean Track Estimator TUBE+2D CORRECTION (All Momenta) 40+ Hits; Truncated Mean; Tracks(AU)", 400, 0, 4);
    TH1F *h_Track_Sizes = new TH1F("Track Sizes", "Track Size; Hits Used; Tracks", 50, 0, 50);
    TH1F *h_momentum = new TH1F("Muon Momentum", "Muon Momentum; Momentum(GeV); Muons / bin", 250, 0, 500);

    // Smaller bins
    const int n_mini_bins = 11;  // Number of histograms for 0-4, 5-9, ..., 55-59
    TH1F* h_Track_TMs_mini[n_mini_bins];
    // Loop to create histograms for each range
    for (int i = 0; i < n_mini_bins; ++i) {
        int lower = 5 + (i * 5);
        int upper = lower + 4;
        TString name = TString::Format("%d-%d Hits", lower, upper);
        TString title = TString::Format("40%% Truncated mean Track Estimator TUBE+2D CORRECTION (All Momenta) %d-%d Hits; Truncated Mean; Tracks(AU)", lower, upper);
        
        h_Track_TMs_mini[i] = new TH1F(name, title, 400, 0, 4);
    }

    // A Tgraph for binned Bethe Bloch
    TGraphErrors *g_BetheBloch= new TGraphErrors();
    TGraphErrors *g_BetheBloch_bg= new TGraphErrors();
    TGraphErrors *g_BetheBlochTM= new TGraphErrors();
    TGraphErrors *g_BetheBlochLandau= new TGraphErrors();

    g_BetheBloch->SetNameTitle("Gaussian Fit to Track de/dx distributions", "Bethe Bloch Curve 1M Z->MuMu; Momentum(GeV); Track Estimator dE/dx");
    g_BetheBloch_bg->SetNameTitle("Gaussian Fit to Track de/dx distributions BetaGamma Axis", "Bethe Bloch Curve 1M Z->MuMu; Beta Gamma; Track Estimator dE/dx (AU)");
    g_BetheBlochTM->SetTitle("Bethe Bloch Curve Truncated Means; Momentum(GeV); Track Estimator Truncated Mean dE/dx");
    g_BetheBlochLandau->SetTitle("Bethe Bloch Curve Landau Fits; Momentum(GeV); Track Estimator Landau MPV dE/dx");
    std::vector<TH1F*> BetheBlochPoints(24);

    // A TGraph for binned Bethe Bloch of hits
    TGraphErrors *g_BetheBlochHitLevelTM = new TGraphErrors();
    TGraphErrors *g_BetheBlochHitLevelLandau = new TGraphErrors();
    TGraphErrors *g_BetheBlochHitLevelLandauBC = new TGraphErrors();
    g_BetheBlochHitLevelTM->SetTitle("Bethe Bloch Curve Hit Level; Momentum(GeV); All hits Truncated Mean");
    g_BetheBlochHitLevelLandau->SetTitle("Bethe Bloch Curve Hit Level; Momentum(GeV); All hits Landau MPV");
    g_BetheBlochHitLevelLandauBC->SetTitle("Bethe Bloch Curve Hit Level BEFORE CORRECTION; Momentum(GeV); All hits Landau MPV");
    std::vector<TH1F*> BetheBlochPointsHitLevel(24);

    std::vector<TH1F*> BetheBlochPointsHitLevelBC(24);

    // Define the initial momentum range and bin count
    int numBins = 150;
    double binWidth = 10.0; // 10 GeV range per histogram
    double startMomentum = 20.0;

    // Loop over the first dimension
    for (size_t i = 0; i < BetheBlochPoints.size(); ++i) {
        // Create histograms for this set
        double minRange = startMomentum + i * binWidth;
        double maxRange = minRange + binWidth;
            
        // Create histogram and assign it to the vector
        TString histName = Form("BetheBlochHist_%lu", i);
        BetheBlochPoints[i] = new TH1F(histName, Form("Track de/dx Distribution %lu-%lu GeV", (unsigned long)minRange, (unsigned long)maxRange),
                                        numBins, 0.5, 1.5);
        TString histName_hl = Form("HitLevelBetheBlochHist_%lu", i);
        BetheBlochPointsHitLevel[i] = new TH1F(histName_hl, Form("Hit Level ADC Count Distribution %lu-%lu GeV", (unsigned long)minRange,
                                        (unsigned long)maxRange), numBins, 0.5, 1.5);
        TString histName_hl_BC = Form("HitLevelBetheBlochHist_%lu BEFORE CORRECTION", i);
        BetheBlochPointsHitLevelBC[i] = new TH1F(histName_hl_BC, Form("Hit Level ADC Count Distribution BEFORE CORRECTION %lu-%lu GeV", (unsigned long)minRange,
                                        (unsigned long)maxRange), numBins, 50, 350);
    }

    // Create a map from stationIndex->stationPhi->stationEta->multilayer->tubelayer->tube->Hits
    std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, std::vector<float>>>>>>> tubeHitMap;

    // A map from chamber to (drift radius, propagation distance) coordinates
    std::map<int, std::map<int, std::map<int, TH3F*>>> tubeCorrectedChamberDistributions;

    // Also create a map with the tube truncated means
    std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, float>>>>>> tubeTMMap;

    // Populate the map with keys
    for (int stationIndex = 0; stationIndex <= 5; ++stationIndex) {
        for (int stationPhi = 1; stationPhi <= 16; ++stationPhi) {
            for (int stationEta = -8; stationEta <= 8; ++stationEta) {
                // Set histogram name
                TString histName = Form("hist_stationIndex%d_stationPhi%d_stationEta%d",
                                                stationIndex, stationPhi, stationEta);
                // Get Chamber type and tube length
                int ctype = chamber_types[stationIndex][stationEta][stationPhiToSector(stationIndex, stationPhi)];
                int tlength = tube_lengths[stationIndex][ctype];
                // cReate histogram and place into map
                TH3F *hist = new TH3F(histName, "Drift radius vs Propagation distance",
                                      10, 0, 15, 10, 0, tlength, 400, 0, 4);
                tubeCorrectedChamberDistributions[stationIndex][stationPhi][stationEta] = hist;
                if (stationEta == 0) continue;  // Skip stationEta == 0
                for (int multilayer = 1; multilayer <= 2; ++multilayer) { 
                    for (int tubelayer = 1; tubelayer <= 4; ++tubelayer) {
                        for (int tube = 1; tube <= 72; ++tube) {
                            // Initialize the vector to store hits (as an example, empty initially)
                            tubeHitMap[stationIndex][stationPhi][stationEta][multilayer][tubelayer][tube] = {};
                        }
                    }
                }
            }
        }
    }

    // Get entries in the total tchain
    unsigned int nEntries = chain.GetEntries();
    
    // Now loop through a second time with the corrected ADCs and make tracks/estimators
    for (unsigned int i=0; i<nEntries; i++) {

        // A progress indicator
        if (i % 10000 == 0) {
            std::cout << i << "/" << nEntries << " processed" << endl;
        }

        // Get the event data
        chain.GetEntry(i);

        // Make sure the current muon has some ADC counts recorded
        if (ADC_counts->size() == 0) { continue; }

        // Loop through Momentum and extract the two best z's
        std::pair<int, int> zpeak_muons (-1,-1);
        int current_closest = 0;
        for (unsigned long n=0; n<pts->size() - 1; n++) { // Loop trhough muons
            // Make sure the first muon has p > 20GeV
            if (ptToP(pts->at(n), etas->at(n)) < 20) { continue; }
            for (unsigned long j=n+1; j<pts->size(); j++) { // Get other muon for a pair
                // Make sure the second in the pair has p > 20GeV
                if (ptToP(pts->at(j), etas->at(j)) < 20) { continue; }

                float mumu_mass = invariantMass(pts->at(n), pts->at(j), etas->at(n), // invariant mass
                etas->at(j), phis->at(n), phis->at(j));
                // Skip pairs not within zmass peak
                if (abs(zmass - mumu_mass) > zmass_window) { 
                    continue; }
                else if (abs(zmass - mumu_mass) < abs(zmass - current_closest)) {
                    current_closest = mumu_mass;
                    zpeak_muons.first = n;
                    zpeak_muons.second = j;
                }
            }
        }

        // Fill the good track momenta
        if (zpeak_muons.first!=-1 && zpeak_muons.second!=-1) {
            h_momentum->Fill(ptToP(pts->at(zpeak_muons.first), etas->at(zpeak_muons.first)));
            h_momentum->Fill(ptToP(pts->at(zpeak_muons.second), etas->at(zpeak_muons.second)));
        }

        // If there isn't a valid pair skip
        if (zpeak_muons.first==-1 || zpeak_muons.second==-1) { continue; }

        // Mark the current muon link
        int muon_link = 0;
        float muon_pt = pts->at(muon_link);
        float muon_eta = etas->at(muon_link);
        float muon_phi = phis->at(muon_link);
        int muon_author = authors->at(muon_link);
        int muon_charge = charges->at(muon_link);
        float muon_p = ptToP(muon_pt, muon_eta);
        int momentum_bin = getMomentumBin(muon_p);
        

        for (unsigned int n=0; n<ADC_counts->size(); n++) {

            // Check if we have a new muon
            if (trk_link->at(n) != muon_link) {
                
                // Update muon
                muon_link = trk_link->at(n);
                muon_pt = pts->at(muon_link);
                muon_eta = etas->at(muon_link);
                muon_phi = phis->at(muon_link);
                muon_p = ptToP(muon_pt, muon_eta);
                muon_author = authors->at(muon_link);
                momentum_bin = getMomentumBin(muon_p);
            }

            // Skip hits not in sector 5 and outlier hits
            if (hit_type->at(n) > 60 || std::abs(muon_eta) > 1 || muon_author != 1 || station_index->at(n) > 5) {
                continue;
            }

            // Only consider the muons from zmass
            if (muon_link != zpeak_muons.first && muon_link != zpeak_muons.second) { continue; }

            // FIll the raw adc histogram
            h_uncorrected_ADC->Fill(ADC_counts->at(n));

            // Fill the tubeHitMap
            tubeHitMap[station_index->at(n)]
            [station_phi->at(n)][station_eta->at(n)]
            [multi_layer->at(n)][tube_layer->at(n)]
            [tubes->at(n)].push_back(ADC_counts->at(n));
            h2_dr_all->Fill(drift_radius->at(n), ADC_counts->at(n));
            h2_momentum_all->Fill(muon_p, ADC_counts->at(n));
            if (momentum_bin == -1) {
                continue;
            } else {
                BetheBlochPointsHitLevelBC[momentum_bin]->Fill(ADC_counts->at(n));
            }
        }
    }

    cout << "Pass 1/3 done" << endl;

    // Now loop through the tubeHitMap and get the TM in every hit
    for (int stationIndex = 0; stationIndex <= 5; ++stationIndex) {
        for (int stationPhi = 1; stationPhi <= 16; ++stationPhi) {
            for (int stationEta = -8; stationEta <= 8; ++stationEta) {
                if (stationEta == 0) continue;  // Skip stationEta == 0
                for (int multilayer = 1; multilayer <= 2; ++multilayer) {
                    for (int tubelayer = 1; tubelayer <= 4; ++tubelayer) {
                        for (int tube = 1; tube <= 72; ++tube) {
                            // Initialize the vector to store hits (as an example, empty initially)
                            float temp_tube_tm = TM(tubeHitMap[stationIndex][stationPhi][stationEta][multilayer][tubelayer][tube], 0.4);
                            tubeTMMap[stationIndex][stationPhi][stationEta][multilayer][tubelayer][tube] = temp_tube_tm;
                            h_tube_tms->Fill(temp_tube_tm);
                        }
                    }
                }
            }
        }
    }

    // A plot of adc distribution after tube correction
    TH1F *h_tube_corrected_ADC = new TH1F("tube corrected ADC", "Global ADC Distribution After Tube Level Correction", 400, 0, 4);

    // A 3D histogram of the the adc after tube correction vs. drift radius and propagation distance
    TH3F *h3_global_tube_corrected_adc = new TH3F("global ADC vs. Drift Radius & Propagation Distance", "Global ADC Distribution vs. Drift Radius and Propagation Distance; Drift Radius(mm); Propagation Distance(mm); ADC Counts / bin", 30, 0, 15, 30, 0, 5000, 400, 0, 4);

    // Now loop through events again dividing every ADC count by the truncated mean of its tube
    for (unsigned int i=0; i<nEntries; i++) {

        // A progress indicator
        if (i % 10000 == 0) {
            std::cout << i << "/" << nEntries << " processed" << endl;
        }

        // Get the event data
        chain.GetEntry(i);

        // Make sure the current muon has some ADC counts recorded
        if (ADC_counts->size() == 0) { continue; }

        // Loop through Momentum and extract the two best z's
        std::pair<int, int> zpeak_muons (-1,-1);
        int current_closest = 0;
        for (unsigned long n=0; n<pts->size() - 1; n++) { // Loop trhough muons
            // Make sure the first muon has p > 20GeV
            if (ptToP(pts->at(n), etas->at(n)) < 20) { continue; }
            for (unsigned long j=n+1; j<pts->size(); j++) { // Get other muon for a pair
                // Make sure the second in the pair has p > 20GeV
                if (ptToP(pts->at(j), etas->at(j)) < 20) { continue; }

                float mumu_mass = invariantMass(pts->at(n), pts->at(j), etas->at(n), // invariant mass
                etas->at(j), phis->at(n), phis->at(j));
                // Skip pairs not within zmass peak
                if (abs(zmass - mumu_mass) > zmass_window) { 
                    continue; }
                else if (abs(zmass - mumu_mass) < abs(zmass - current_closest)) {
                    current_closest = mumu_mass;
                    zpeak_muons.first = n;
                    zpeak_muons.second = j;
                }
            }
        }

        // If there isn't a valid pair skip
        if (zpeak_muons.first==-1 || zpeak_muons.second==-1) { continue; }

        // Mark the current muon link
        int muon_link = 0;
        float muon_pt = pts->at(muon_link);
        float muon_eta = etas->at(muon_link);
        float muon_phi = phis->at(muon_link);
        int muon_author = authors->at(muon_link);
        int muon_charge = charges->at(muon_link);
        float muon_p = ptToP(muon_pt, muon_eta);

        for (unsigned int n=0; n<ADC_counts->size(); n++) {

            // Check if we have a new muon
            if (trk_link->at(n) != muon_link) {
                
                // Update muon
                muon_link = trk_link->at(n);
                muon_pt = pts->at(muon_link);
                muon_eta = etas->at(muon_link);
                muon_phi = phis->at(muon_link);
                muon_p = ptToP(muon_pt, muon_eta);
                muon_author = authors->at(muon_link);
            }

            // Skip hits not in sector 5 and outlier hits
            if (hit_type->at(n) > 60 || std::abs(muon_eta) > 1 || muon_author != 1 || station_index->at(n) > 5) {
                continue;
            }

            // Only consider the muons from zmass
            if (muon_link != zpeak_muons.first && muon_link != zpeak_muons.second) { continue; }

            // Find propagation distance
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
            h2_dprop_all->Fill(distance_to_propagate, ADC_counts->at(n));

            // The readout is on opposite sides for large sectors and small sectors
            if (station_index->at(n) % 2 == 0) {
                distance_to_propagate = (tube_length/2) + a;
            }
            else if (station_index->at(n) % 2 == 1) {
                distance_to_propagate = (tube_length/2) - a;
            }

            // Fill the tubeHitMap
            if (tubeHitMap[station_index->at(n)][station_phi->at(n)][station_eta->at(n)][multi_layer->at(n)][tube_layer->at(n)]
            [tubes->at(n)].size() < 50) {
                continue;
            }
            float tube_tm = tubeTMMap[station_index->at(n)]
            [station_phi->at(n)][station_eta->at(n)]
            [multi_layer->at(n)][tube_layer->at(n)]
            [tubes->at(n)];
            float corrected_ADC = ADC_counts->at(n) / tube_tm;
            h_tube_corrected_ADC->Fill(corrected_ADC);
            h3_global_tube_corrected_adc->Fill(drift_radius->at(n), distance_to_propagate, corrected_ADC);

            // Fill the corrected adc coordinate map
            tubeCorrectedChamberDistributions[station_index->at(n)][station_phi->at(n)][station_eta->at(n)]->Fill(drift_radius->at(n), distance_to_propagate, corrected_ADC);
        }
    }

    tubeHitMap.clear();

    cout << "2nd pass finished" << endl;

    // Get a single 2D function in the entire tube corrected barrel
    std::vector<float> global_func = processTMs(h3_global_tube_corrected_adc, h3_global_tube_corrected_adc->GetName());

    
    // Now Get the 2D functions in every tube corrected chamber
    std::map<int, std::map<int, std::map<int, std::vector<float>>>> chamber_funcs;
    for (int stationIndex = 0; stationIndex <= 5; ++stationIndex) {
        for (int stationPhi = 1; stationPhi <= 16; ++stationPhi) {
            for (int stationEta = -8; stationEta <= 8; ++stationEta) {
                if (stationEta==0) { continue; }
                TH3F *temp_hist = tubeCorrectedChamberDistributions[stationIndex][stationPhi][stationEta];
                // Make sure there are hits in the current chamber
                if (temp_hist->GetEntries()==0) {
                    continue;
                }
                std::vector<float> zValues;

                // Loop through all bins and collect z-values
                for (int x = 1; x <= temp_hist->GetNbinsX(); ++x) {
                    for (int y = 1; y <= temp_hist->GetNbinsY(); ++y) {
                        for (int z = 1; z <= temp_hist->GetNbinsZ(); ++z) {
                            float binContent = temp_hist->GetBinContent(x, y, z);
                            if (binContent != 0) {  // Exclude empty bins
                                zValues.push_back(temp_hist->GetZaxis()->GetBinCenter(z) * binContent);
                            }
                        }
                    }
                }

                // Sort in descending order
                std::sort(zValues.begin(), zValues.end());
                float cutoff = floor((1-0.4) * zValues.size());
                int hit_count = 0;
                float sum = 0.0;
                if (hit_count < cutoff) {
                    sum += zValues[hit_count];
                    hit_count++;
                }
                sum /= cutoff;
                h_chamber_tms->Fill(sum);

                std::vector<float> temp_chamber_params = processTMs(temp_hist, temp_hist->GetName());
                chamber_funcs[stationIndex][stationPhi][stationEta] = temp_chamber_params;
            }
        }
    }

    // A plot of the ADC distribution after tube + 2D correction
    TH1F *h_full_corrected_ADC = new TH1F("Fully corrected ADC", "Global ADC Distribution After Tube Level + 2D Correction", 400, 0, 4);
    TH2F *h2_full_corrected_ADC_dr = new TH2F("Fully corrected ADC vs. Drift Radius", "Global ADC Distribution vs. Drift Radius After Tube Level + 2D Correction; Drift Radius(mm); ADC Counts", 15, 0, 15, 400, 0, 4);
    TH2F *h2_full_corrected_ADC_dprop = new TH2F("Fully corrected ADC vs. Propagation Distance", "Global ADC Distribution vs. Propagation Distance After Tube Level + 2D Correction; Propagation Distance(mm); ADC Counts", 50, 0, 5000, 400, 0, 4);

    cout << "3rd pass" << endl;
    // Now Loop through the events again applying the 2D corrections to tube corrected hits and find tracks
    for (unsigned int i=0; i<nEntries; i++) {

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


        // Loop through Momentum and extract the two best z's
        std::pair<int, int> zpeak_muons (-1,-1);
        int current_closest = 0;
        for (unsigned long n=0; n<pts->size() - 1; n++) { // Loop trhough muons
            // Make sure the first muon has p > 20GeV
            if (ptToP(pts->at(n), etas->at(n)) < 20) { continue; }
            for (unsigned long j=n+1; j<pts->size(); j++) { // Get other muon for a pair
                // Make sure the second in the pair has p > 20GeV
                if (ptToP(pts->at(j), etas->at(j)) < 20) { continue; }

                float mumu_mass = invariantMass(pts->at(n), pts->at(j), etas->at(n), // invariant mass
                etas->at(j), phis->at(n), phis->at(j));
                // Skip pairs not within zmass peak
                if (abs(zmass - mumu_mass) > zmass_window) { 
                    continue; }
                else if (abs(zmass - mumu_mass) < abs(zmass - current_closest)) {
                    current_closest = mumu_mass;
                    zpeak_muons.first = n;
                    zpeak_muons.second = j;
                }
            }
        }

        // If there isn't a valid pair skip
        if (zpeak_muons.first==-1 || zpeak_muons.second==-1) { continue; }

        // Reset current track
        temp_track.clear();
        temp_track_stations.clear();
        temp_track_etas.clear();

        // Mark the current muon link
        int muon_link = 0;
        float muon_pt = pts->at(muon_link);
        float muon_eta = etas->at(muon_link);
        float muon_phi = phis->at(muon_link);
        int muon_author = authors->at(muon_link);
        int muon_charge = charges->at(muon_link);
        float muon_p = ptToP(muon_pt, muon_eta);
        int momentum_bin = getMomentumBin(muon_p);

        // Cut on momentum here
        if (momentum_bin == -1) {
            continue;
        }

        for (unsigned int n=0; n<ADC_counts->size(); n++) {

            // Check if we have a new muon
            if (trk_link->at(n) != muon_link) {
                
                // Make sure the track isn't empty and fill the old muon
                if (!temp_track_stations.empty()) {
                    track_tm = TM(temp_track, 0.4);

                    // Make sure the track is contained in the barrel
                    auto max_station = std::max_element(temp_track_stations.begin(),temp_track_stations.end());
                    if (*max_station < 6) {
                        h_Track_TMs->Fill(track_tm);
                        BetheBlochPoints.at(momentum_bin)->Fill(track_tm);
                        if (temp_track.size() < 15) {
                            h_Track_TMs_1->Fill(track_tm);
                        } else if (temp_track.size() >= 15 && temp_track.size() < 30) {
                            h_Track_TMs_2->Fill(track_tm);
                        } else if (temp_track.size() >= 30 && temp_track.size() < 39) {
                            h_Track_TMs_3->Fill(track_tm);
                        } else if (temp_track.size() >= 40) {
                            h_Track_TMs_4->Fill(track_tm);
                        }
                        h_Track_Sizes->Fill(temp_track.size());

                        // Fill the smaller track size binned histograms
                        if (temp_track.size()==0) {
                            continue;
                        }
                        int mini_bin = (temp_track.size() - 5) / 5;
                        if (mini_bin > 10 || mini_bin < 0) { continue; }
                        else {
                            h_Track_TMs_mini[mini_bin]->Fill(track_tm);
                        }
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
                momentum_bin = getMomentumBin(muon_p);
            };

            // Skip hits not in sector 5 and outlier hits
            if (hit_type->at(n) > 60 || std::abs(muon_eta) > 1 || muon_author != 1 || station_index->at(n) > 5) {
                continue;
            }

            // Only consider the muons from zmass
            if (muon_link != zpeak_muons.first && muon_link != zpeak_muons.second) { continue; }

            // if (muon_p < 30 || muon_p > 40) { continue; }

            // Find propagation distance
            // Get the length of the current mdt tube
            int chamber_type = chamber_types[station_index->at(n)][abs(station_eta->at(n))][station_phi->at(n)];
            int tube_length = tube_lengths[station_index->at(n)][chamber_type];

            // Make sure that the chamber is instantiated
            if (chamber_type == 0) {
                continue;
            }

            // Find the location of the tube where the hit occured
            float radial_distance = tube_positions[station_index->at(n)][multi_layer->at(n)][tube_layer->at(n)];
            int sector = stationPhiToSector(station_index->at(n), station_phi->at(n));
            float temp_phi; // use this temp phi to rotate all of the tracks as if they are in sector 5 for geometry purposes 
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

            // Fill the tubeTMMap
            float tube_tm = tubeTMMap[station_index->at(n)]
            [station_phi->at(n)][station_eta->at(n)]
            [multi_layer->at(n)][tube_layer->at(n)]
            [tubes->at(n)];

            float once_corrected_ADC = ADC_counts->at(n) / tube_tm;

            if (tubeCorrectedChamberDistributions[station_index->at(n)][station_phi->at(n)][station_eta->at(n)]->GetEntries() == 0) {
                continue;
            }
            
            // Below is the corrected ADC if you use a unique function in every chamber
            float twice_corrected_ADC = once_corrected_ADC / correctionFunction(drift_radius->at(n),
            distance_to_propagate, chamber_funcs[station_index->at(n)][station_phi->at(n)][station_eta->at(n)]);
            h_full_corrected_ADC->Fill(twice_corrected_ADC);
            h2_full_corrected_ADC_dr->Fill(drift_radius->at(n), twice_corrected_ADC);
            h2_full_corrected_ADC_dprop->Fill(distance_to_propagate, twice_corrected_ADC);

            // Below is the corrected ADC if you use a single fucntion for every chamber
            // float twice_corrected_ADC = once_corrected_ADC / correctionFunction(drift_radius->at(n),
            // distance_to_propagate, global_func);
            
            // Cut on momentum here
            if (momentum_bin == -1) {
                continue;
            }

            BetheBlochPointsHitLevel[momentum_bin]->Fill(twice_corrected_ADC);

            temp_track.push_back(twice_corrected_ADC);
            temp_track_stations.push_back(station_index->at(n));
            temp_track_etas.push_back(station_eta->at(n));
        }

        // Add the last muon in the event
        // Make sure the track isn't empty
        if (!temp_track_stations.empty()) {
            
            // Get track TM
            float track_tm = TM(temp_track, 0.4);

            // Make sure the track is contained in the barrel
            auto max_station = std::max_element(temp_track_stations.begin(),temp_track_stations.end());
            if (*max_station < 6) {
                h_Track_TMs->Fill(track_tm);
                BetheBlochPoints.at(momentum_bin)->Fill(track_tm);
                if (temp_track.size() < 15) {
                    h_Track_TMs_1->Fill(track_tm);
                } else if (temp_track.size() >= 15 && temp_track.size() < 30) {
                    h_Track_TMs_2->Fill(track_tm);
                    // BetheBlochPoints.at(momentum_bin)->Fill(track_tm);
                } else if (temp_track.size() >= 30 && temp_track.size() < 39) {
                    h_Track_TMs_3->Fill(track_tm);
                    // BetheBlochPoints.at(momentum_bin)->Fill(track_tm);
                } else if (temp_track.size() >= 40) {
                    h_Track_TMs_4->Fill(track_tm);
                    // BetheBlochPoints.at(momentum_bin)->Fill(track_tm);
                }
                h_Track_Sizes->Fill(temp_track.size());

                // Fill the smaller track size binned histograms
                int mini_bin = (temp_track.size() - 5) / 5;
                cout  << "mini bin " << mini_bin << " for " << temp_track.size() << " hits on track" << endl;
                if (mini_bin > 10 || mini_bin < 0) { continue; }
                else {
                    h_Track_TMs_mini[mini_bin]->Fill(track_tm);
                }
            }
        }
        temp_track.clear();
        temp_track_stations.clear();
        temp_track_etas.clear();
    }

    cout << "Starting Bethe Bloch" << endl;
    gStyle->SetOptFit(1);
    //Loop through Momentum Bins to Get average track_tm in each bin
    for (unsigned int i=0; i<BetheBlochPoints.size(); i++) {

        // ADd points to all graphs for this momentum bin
        float momentum_bin_center = 25 + i*10;

        // Get the histogram from each momentum Bin and fit a gaussian
        TH1F *temp_hist2 = BetheBlochPoints[i];
        if (temp_hist2->GetEntries() > 0) {
            TF1 *BB_bin_fit = new TF1(Form("fit_%d", i), "gaus", temp_hist2->GetXaxis()->GetXmin(), temp_hist2->GetXaxis()->GetXmax());
            temp_hist2->Fit(BB_bin_fit, "RQ");
            temp_hist2->Write();
            

            // Extract fit parameters
            double mean = BB_bin_fit->GetParameter(1);  // Mean of Gaussian
            double sigma = BB_bin_fit->GetParameter(2); // Standard deviation of Gaussian
            std::pair<float, float> tms = fit1DHistogramTM(temp_hist2);
            std::pair<float, float> mpvs = fit1DHistogramLandau(temp_hist2);

            // Add point to Bethe-Bloch graph
            g_BetheBloch->AddPointError(momentum_bin_center, mean, 0, BB_bin_fit->GetParError(1));
            g_BetheBloch_bg->AddPointError(momentum_bin_center / 0.1057, mean, 0, BB_bin_fit->GetParError(1));
            g_BetheBlochTM->AddPointError(momentum_bin_center, tms.first, 0, tms.second);
            g_BetheBlochLandau->AddPointError(momentum_bin_center, mpvs.first, 0, mpvs.second);
        }

        // Get the hit level histograms
        TH1F *hist_corrected_hit_level = BetheBlochPointsHitLevel[i];
        TH1F *hist_raw_hit_level = BetheBlochPointsHitLevelBC[i];
        if (hist_corrected_hit_level->GetEntries() > 0) {
            // Extract fit parameters
            std::pair<float, float> hit_level_tms = fit1DHistogramTM(hist_corrected_hit_level);
            std::pair<float, float> hit_level_mpvs = fit1DHistogramLandau(hist_corrected_hit_level);
            std::pair<float, float> hit_level_mpvs_BC = fit1DHistogramLandau(hist_raw_hit_level);
            g_BetheBlochHitLevelTM->AddPointError(momentum_bin_center, hit_level_tms.first, 0, hit_level_tms.second);
            g_BetheBlochHitLevelLandau->AddPointError(momentum_bin_center, hit_level_mpvs.first, 0, hit_level_mpvs.second);
            g_BetheBlochHitLevelLandauBC->AddPointError(momentum_bin_center, hit_level_mpvs_BC.first, 0, hit_level_mpvs_BC.second);
        }
        
        hist_corrected_hit_level->Write();
        hist_raw_hit_level->Write();
    }

    g_BetheBloch->Write();
    g_BetheBloch_bg->Write();
    g_BetheBlochTM->Write();
    g_BetheBlochLandau->Write();
    g_BetheBlochHitLevelTM->Write();
    g_BetheBlochHitLevelLandau->Write();
    g_BetheBlochHitLevelLandauBC->Write();

    // Create a canvas to draw the beta gamma bounds on the Bethe Bloch curve
    TCanvas *c_bg = new TCanvas("c_bg", "Bethe-Bloch Curve", 800, 600);
    g_BetheBloch_bg->SetNameTitle("Gaussian Fit to Track de/dx distributions BetaGamma Axis", "Bethe Bloch Curve 1M Z->MuMu; Beta Gamma; Track Estimator dE/dx (AU)");
    g_BetheBloch_bg->SetMarkerStyle(7);
    g_BetheBloch_bg->SetMarkerColor(kBlue);
    g_BetheBloch_bg->Draw("AP");

    // Draw vertical lines at x = 200 and x = 1000
    TLine *line1 = new TLine(200, gPad->GetUymin(), 200, gPad->GetUymax());
    TLine *line2 = new TLine(1000, gPad->GetUymin(), 1000, gPad->GetUymax());

    // Set line styles and colors
    line1->SetLineColor(kTeal+3);
    line1->SetLineStyle(2);
    line1->SetLineWidth(2);
    line1->Draw();

    line2->SetLineColor(kGreen+3);
    line2->SetLineStyle(2);
    line2->SetLineWidth(2);
    line2->Draw();

    // Update the canvas
    c_bg->Update();
    c_bg->Write();

    // Characterization
    h2_dr_all->Write();
    h2_dprop_all->Write();
    h2_momentum_all->Write();
    h_momentum->Write();

    TF1 *tube_fitFunc = new TF1("tube_fitFunc", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[1])/[4])^2)", 0, 4);
    h_tube_tms->Fit(tube_fitFunc, "R");
    h_tube_tms->Write();
    TF1 *chamber_fitFunc = new TF1("chamber_fitFunc", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[1])/[4])^2)", 0, 4);
    h_chamber_tms->Fit(chamber_fitFunc, "R");
    h_chamber_tms->Write();

    // Write the uncorrected track estimators
    h_Track_TMs->Write();
    h_Track_Sizes->Write();
    h_uncorrected_ADC->Write();
    h_tube_corrected_ADC->Write();
    h_full_corrected_ADC->Write();
    h2_full_corrected_ADC_dr->Write();
    h2_full_corrected_ADC_dprop->Write();
    h3_global_tube_corrected_adc->Write();

    // Fit the fully corrected ADC distribution to a landau
    TCanvas *c_landau = new TCanvas();
    h_full_corrected_ADC->Draw();
    TF1 *corrected_landau = new TF1("corrected landau", "landau");
    h_full_corrected_ADC->Fit(corrected_landau);
    corrected_landau->Draw("Same");
    TPaveText *landau_label = new TPaveText(0.3, 0.6, 0.9, 0.9, "NDC");
    landau_label->AddText(Form("MPV: %f±%f", corrected_landau->GetParameter(1), corrected_landau->GetParError(1)));
    landau_label->AddText(Form("Width: %f±%f", corrected_landau->GetParameter(2), corrected_landau->GetParError(2)));
    landau_label->Draw("same");
    c_landau->Write();
    delete c_landau;
    delete corrected_landau;
    delete landau_label;

    // Fit a double Gaussian function to the histogram
    TF1 *tm_fitFunc = new TF1("tm_fitFunc", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[1])/[4])^2)", 0, 4);
    tm_fitFunc->SetParameters(0.05, 1, 0.08, 0.005, 0.08);
    h_Track_TMs->Fit(tm_fitFunc, "R");

    // Extract fit parameters
    double mean = tm_fitFunc->GetParameter(1);
    double meanErr = tm_fitFunc->GetParError(1);
    double sigma1 = tm_fitFunc->GetParameter(2);
    double sigma1Err = tm_fitFunc->GetParError(2);
    double amp1 = tm_fitFunc->GetParameter(0);
    double amp1Err = tm_fitFunc->GetParError(0);
    double amp2 = tm_fitFunc->GetParameter(3);
    double amp2Err = tm_fitFunc->GetParError(3);
    double sigma2 = tm_fitFunc->GetParameter(4);
    double sigma2Err = tm_fitFunc->GetParError(4);
    
    // Choose the Gaussian with the larger amplitude
    double chosenSigma, chosenSigmaErr, resolution, resErr;
    if (amp1 > amp2) {
        chosenSigma = sigma1;
        chosenSigmaErr = sigma1Err;
    } else {
        chosenSigma = sigma2;
        chosenSigmaErr = sigma2Err;
    }

    // Calculate relative resolution and its uncertainty
    resolution = chosenSigma / mean;
    resErr = resolution * sqrt(pow(chosenSigmaErr / chosenSigma, 2) + pow(meanErr / mean, 2));

    // Optionally draw the histogram with the fit
    TCanvas *c1 = new TCanvas("Track Estimators Fit", "Fit Canvas", 800, 600);
    h_Track_TMs->Draw();
    tm_fitFunc->Draw("same");
    TPaveText *label = new TPaveText(0.3, 0.6, 0.9, 0.9, "NDC");
    label->AddText(Form("amp1: %f +/- %f", amp1, amp1Err));
    label->AddText(Form("amp2: %f +/- %f", amp2, amp2Err));
    label->AddText(Form("Sigma 1: %f +/- %f", sigma1, sigma1Err));
    label->AddText(Form("Sigma 2: %f +/- %f", sigma2, sigma2Err));
    label->AddText(Form("Mean: %f +/- %f", mean, meanErr));
    label->AddText(Form("resolution: %f ± %f", resolution, resErr));
    label->AddText(Form("Chisquared: %f", tm_fitFunc->GetChisquare()));
    label->AddText(Form("NDF: %d", tm_fitFunc->GetNDF()));
    label->Draw("same");
    c1->Write();

    h_Track_TMs_1->Write();
    h_Track_TMs_2->Write();
    h_Track_TMs_3->Write();
    h_Track_TMs_4->Write();
    
    /*
    // Define histograms array
    TH1 *histograms[4] = {h_Track_TMs_4, h_Track_TMs_3, h_Track_TMs_2, h_Track_TMs_1};
    TCanvas *c2 = new TCanvas("Track Estimators Fit Different Size", "Different track Sizes", 800, 600);
    TLegend *leg = new TLegend(0.1, 0.65, 0.43, 0.9);
    leg->Draw();
    histograms[0]->SetLineColor(1);
    histograms[0]->Scale(1.0 / histograms[0]->Integral());
    histograms[0]->GetYaxis()->SetRangeUser(0,1);
    histograms[0]->Draw();


    // Loop through histograms for estimators from different number of hits on track
    for (int i = 0; i < n_mini_bins; i++) {
        histograms[i]->Scale(1.0 / histograms[i]->Integral());
        if (i > 0) {
            histograms[i]->GetYaxis()->SetRangeUser(0,1);
            histograms[i]->Draw("same");
        }
        TCanvas *c_hist = new TCanvas(Form("Track Estimators for %s", histograms[i]->GetName()), "Different Track Sizes", 800, 600);
        c_hist->cd();
        histograms[i]->Draw();

        // Define a double Gaussian fit function with the same mean
        TF1 *fitFunc = new TF1("fitFunc", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[1])/[4])^2)", 0, 4);
        fitFunc->SetLineColor(i+1);
        fitFunc->SetParameters(0.05, 1, 0.08, 0.005, 0.08); // Initial guesses: amplitude1, mean, sigma1, amplitude2, sigma2
        // Constrain amplitudes to be positive
        fitFunc->SetParLimits(0, 0, 1); // Constrain amplitude1 to be positive
        fitFunc->SetParLimits(1, 0, 4);
        fitFunc->SetParLimits(2, 0, 100);
        fitFunc->SetParLimits(3, 0, 1); // Constrain amplitude2 to be positive
        fitFunc->SetParLimits(4, 0, 100);

        // Fit the histogram
        histograms[i]->Fit(fitFunc, "R");
        histograms[i]->SetLineColor(i+1);
        fitFunc->Draw("same");
        c2->cd();
        fitFunc->Draw("same");

        // Extract fit parameters
        double mean = fitFunc->GetParameter(1);
        double sigma1 = fitFunc->GetParameter(2);
        double sigma2 = fitFunc->GetParameter(4);
        double amp1 = fitFunc->GetParameter(0);
        double amp2 = fitFunc->GetParameter(3);
        double amp1Err = fitFunc->GetParError(0);
        double amp2Err = fitFunc->GetParError(3);

        // Extract fit parameter uncertainties
        double meanErr = fitFunc->GetParError(1);
        double sigma1Err = fitFunc->GetParError(2);
        double sigma2Err = fitFunc->GetParError(4);

        // Choose the Gaussian with the larger amplitude
        double chosenSigma, chosenSigmaErr, resolution, resError;
        if (amp1 > amp2) {
            chosenSigma = sigma1;
            chosenSigmaErr = sigma1Err;
        } else {
            chosenSigma = sigma2;
            chosenSigmaErr = sigma2Err;
        }

        // Calculate relative resolution and its uncertainty
        resolution = chosenSigma / mean;
        resError = resolution * sqrt(pow(chosenSigmaErr / chosenSigma, 2) + pow(meanErr / mean, 2));

        // Create a label for the current track size plot
        TPaveText *temp_label = new TPaveText(0.3, 0.6, 0.9, 0.9, "NDC");
        temp_label->AddText(Form("Chi^2: %f", fitFunc->GetChisquare()));
        temp_label->AddText(Form("NDF: %d", fitFunc->GetNDF()));
        temp_label->AddText(Form("resolution: %f ± %f", resolution, resError));
        temp_label->AddText(Form("Mean: %f ± %f", mean, meanErr));
        temp_label->AddText(Form("Amp1: %f ± %f | Sigma1: %f ± %f", amp1, amp1Err, sigma1, sigma1Err));
        temp_label->AddText(Form("Amp2: %f ± %f | Sigma2: %f ± %f", amp2, amp2Err, sigma2, sigma2Err));
        c_hist->cd();
        temp_label->Draw("same");
        c_hist->Write();
        delete temp_label;
        delete c_hist;
        
        // Add TPaveText to display stats
        TString label1 = Form("%s: #sigma/#mu - %.3f +/- %.3f", histograms[i]->GetName(), resolution, resError);
        leg->AddEntry(fitFunc, label1.Data(), "l");
    }
    c2->cd();
    leg->Draw();

    // Save the canvas
    c2->Write();
    */

    // Loop through track sized binned histograms
    TCanvas *c2 = new TCanvas("Track Estimators Fit Different Size", "Different track Sizes", 800, 600);
    TLegend *leg = new TLegend(0.1, 0.65, 0.43, 0.9);
    leg->Draw();
    h_Track_TMs_mini[0]->SetLineColor(1);
    h_Track_TMs_mini[0]->Scale(1.0 / h_Track_TMs_mini[0]->Integral());
    h_Track_TMs_mini[0]->GetYaxis()->SetRangeUser(0,1);
    h_Track_TMs_mini[0]->Draw();

    TGraphErrors *g_Resolution_sizes = new TGraphErrors();
    g_Resolution_sizes->SetNameTitle("Estimator Resolution vs. Hits on Track", "Relative Resolution vs. Number of Hits on Track; Number of Hits on Track; Relative Resolution");
    g_Resolution_sizes->SetMarkerStyle(7);

    // Loop through histograms for estimators from different number of hits on track
    for (int i = 0; i < n_mini_bins; i++) {
        h_Track_TMs_mini[i]->Scale(1.0 / h_Track_TMs_mini[i]->Integral());
        if (i > 0) {
            h_Track_TMs_mini[i]->GetYaxis()->SetRangeUser(0,1);
            h_Track_TMs_mini[i]->Draw("same");
        }
        TCanvas *c_hist = new TCanvas(Form("Track Estimators for %s", h_Track_TMs_mini[i]->GetName()), "Different Track Sizes", 800, 600);
        c_hist->cd();
        h_Track_TMs_mini[i]->Draw();

        // Define a double Gaussian fit function with the same mean
        TF1 *fitFunc = new TF1("fitFunc", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[1])/[4])^2)", 0, 4);
        fitFunc->SetLineColor(i+1);
        fitFunc->SetParameters(0.05, 1, 0.08, 0.005, 0.08); // Initial guesses: amplitude1, mean, sigma1, amplitude2, sigma2
        
        // Constrain amplitudes to be positive
        fitFunc->SetParLimits(0, 0, 1); // Constrain amplitude1 to be positive
        fitFunc->SetParLimits(1, 0, 4);
        fitFunc->SetParLimits(2, 0, 100);
        fitFunc->SetParLimits(3, 0, 1); // Constrain amplitude2 to be positive
        fitFunc->SetParLimits(4, 0, 100);

        // Fit the histogram
        h_Track_TMs_mini[i]->Fit(fitFunc, "R");
        h_Track_TMs_mini[i]->SetLineColor(i+1);
        fitFunc->Draw("same");
        c2->cd();
        fitFunc->Draw("same");

        // Extract fit parameters
        double mean = fitFunc->GetParameter(1);
        double sigma1 = fitFunc->GetParameter(2);
        double sigma2 = fitFunc->GetParameter(4);
        double amp1 = fitFunc->GetParameter(0);
        double amp2 = fitFunc->GetParameter(3);
        double amp1Err = fitFunc->GetParError(0);
        double amp2Err = fitFunc->GetParError(3);

        // Extract fit parameter uncertainties
        double meanErr = fitFunc->GetParError(1);
        double sigma1Err = fitFunc->GetParError(2);
        double sigma2Err = fitFunc->GetParError(4);

        // Choose the Gaussian with the larger amplitude
        double chosenSigma, chosenSigmaErr, resolution, resError;
        if (amp1 > amp2) {
            chosenSigma = sigma1;
            chosenSigmaErr = sigma1Err;
        } else {
            chosenSigma = sigma2;
            chosenSigmaErr = sigma2Err;
        }

        // Calculate relative resolution and its uncertainty
        resolution = chosenSigma / mean;
        resError = resolution * sqrt(pow(chosenSigmaErr / chosenSigma, 2) + pow(meanErr / mean, 2));

        cout << "Adding point " << resolution << "+\\-" << resError << " at " << 2.5 + (5 * i) << " hits on track" << endl;
        g_Resolution_sizes->AddPointError(7.5 + (5 * i), resolution, 0, resError);

        // Create a label for the current track size plot
        TPaveText *temp_label = new TPaveText(0.3, 0.6, 0.9, 0.9, "NDC");
        temp_label->AddText(Form("Chi^2: %f", fitFunc->GetChisquare()));
        temp_label->AddText(Form("NDF: %d", fitFunc->GetNDF()));
        temp_label->AddText(Form("resolution: %f ± %f", resolution, resError));
        temp_label->AddText(Form("Mean: %f ± %f", mean, meanErr));
        temp_label->AddText(Form("Amp1: %f ± %f | Sigma1: %f ± %f", amp1, amp1Err, sigma1, sigma1Err));
        temp_label->AddText(Form("Amp2: %f ± %f | Sigma2: %f ± %f", amp2, amp2Err, sigma2, sigma2Err));
        c_hist->cd();
        temp_label->Draw("same");
        c_hist->Write();
        delete temp_label;
        delete c_hist;
        
        // Add TPaveText to display stats
        TString label1 = Form("%s: #sigma/#mu - %.3f +/- %.3f", h_Track_TMs_mini[i]->GetName(), resolution, resError);
        leg->AddEntry(fitFunc, label1.Data(), "l");
    }
    c2->cd();
    leg->Draw();

    // Save the canvas
    c2->Write();

    g_Resolution_sizes->Write();

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
    tubeCalibration(indir, outfile);
    return 0;
}