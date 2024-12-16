#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
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
#include <TBranch.h>
#include <TMath.h>
#include <TF1.h>
#include <TPaveText.h>
#include <TLegend.h>
using namespace std;



// A function to check if a hit is in the barrel
bool isInBarrel(std::set<int> stations) {
    // Station indices are saved as hexidecimal strings. These are the strings for the layers listed above
    std::set<int> barrel_layers = {0, 1, 2, 3, 4, 5};
    // See if any of the stations in stations are not in the barrel
    for (const auto& station : stations) {
        if (barrel_layers.find(station) == barrel_layers.end()) {
            return false; // If any element is not found in barrel_layers, return false
        }
    }
    // All stations are in the barrel if the above loop is psased without a return
    return true;
}

// A function which cuts on a specific pseudorapidity
bool isInValidEta(std::vector<float> etas, float cutoff) {
    for (const auto& eta : etas) {
        if (abs(eta) > cutoff) {
            return false;
        }
    }
    return true;
}

// A function to convert transverse momentum to total |p|
float ptToP(float pt, float eta) {
    float p = exp(-eta);
    p = 2 * atan(p);
    p = pt / sin(p);
    return p;
}

// A function to convert momentum to Beta*gammae
float pToBG(float p) {
    float muon_rest_mass = 0.105658; //GeV/c
    float bg = p / muon_rest_mass;
    return bg;
}

// A function to calculate the invariant mass of two particles created in a collision
float invariantMass(float pt1, float pt2, float eta1, float eta2, float phi1, float phi2) {
    return sqrt(2 * pt1 * pt2 * (cosh(eta1 - eta2) - cos(phi1 - phi2)));
}

// A function to check if a track has hits in overlapping regions
std::vector<int> overlaps(std::set<int> stations) {
    // Set the 1st, 2nd, and 3rd element to 1 if those layers have overlapping hits respectively
    std::vector<int> status = {0,0,0};

    // Check if there is overlap in barrel
    std::set<int> inner_layer = {0, 1};
    std::set<int> middle_layer = {2, 3};
    std::set<int> outer_layer = {4, 5};
    if (std::includes(stations.begin(), stations.end(), inner_layer.begin(), inner_layer.end()) == true) {
        status.at(0) = 1; 
    }
    if (std::includes(stations.begin(), stations.end(), middle_layer.begin(), middle_layer.end()) == true) {
        status.at(1) = 1; 
    }
    if (std::includes(stations.begin(), stations.end(), outer_layer.begin(), outer_layer.end()) == true) {
        status.at(2) = 1; 
    }
    return status;
}

// A fucntion to count how many layers a track passes through
int layerCounter(std::map<int, std::vector<float>> hitMap) {
    int layer_count = 0;
    if (!hitMap[0].empty() || !hitMap[1].empty()) { layer_count++; }
    if (!hitMap[2].empty() || !hitMap[3].empty()) { layer_count++; }
    if (!hitMap[4].empty() || !hitMap[5].empty()) { layer_count++; }
    return layer_count;
}

// A function to test if a muon event is a full track extending through all barrel layers
bool isFullTrack(std::map<int, std::vector<float>> hitMap) {
    int layer_count = layerCounter(hitMap);
    return (layer_count == 3);
}

// A function to build a track hit vector without removing overlaps
std::vector<float> buildTracks(std::map<int, std::vector<float>> hitMap) {
    std::vector<float> track;
    // Retrieve the hits from every barrel station
    std::vector<float> BIL_hits = hitMap[0];
    std::vector<float> BIS_hits = hitMap[1];
    std::vector<float> BML_hits = hitMap[2];
    std::vector<float> BMS_hits = hitMap[3];
    std::vector<float> BOL_hits = hitMap[4];
    std::vector<float> BOS_hits = hitMap[5];

    // Glue all the hits together
    track.insert(track.end(), BIL_hits.begin(), BIL_hits.end());
    track.insert(track.end(), BIS_hits.begin(), BIS_hits.end());
    track.insert(track.end(), BML_hits.begin(), BML_hits.end());
    track.insert(track.end(), BMS_hits.begin(), BMS_hits.end());
    track.insert(track.end(), BOL_hits.begin(), BOL_hits.end());
    track.insert(track.end(), BOS_hits.begin(), BOS_hits.end());

    return track;
}

// An overloaded function to build a track hit vector accounting for overlap regions
std::vector<float> buildTracks(std::map<int, std::vector<float>> hitMap, std::set<int> stations) {
    std::vector<float> track;
    std::vector<float> BIL_hits = hitMap[0];
    std::vector<float> BIS_hits = hitMap[1];
    std::vector<float> BML_hits = hitMap[2];
    std::vector<float> BMS_hits = hitMap[3];
    std::vector<float> BOL_hits = hitMap[4];
    std::vector<float> BOS_hits = hitMap[5];

    // Glue all the hits together
    track.insert(track.end(), BIL_hits.begin(), BIL_hits.end());
    track.insert(track.end(), BIS_hits.begin(), BIS_hits.end());
    track.insert(track.end(), BML_hits.begin(), BML_hits.end());
    track.insert(track.end(), BMS_hits.begin(), BMS_hits.end());
    track.insert(track.end(), BOL_hits.begin(), BOL_hits.end());
    track.insert(track.end(), BOS_hits.begin(), BOS_hits.end());

    // Get the overlap statuses in the inner and middle layers
    std::vector<int> status = overlaps(stations);
    //If there is no overlap in the inner region. Add both chamber types hits (one will be empty)
    if (status.at(0) == 0) {
        track.insert(track.end(), BIL_hits.begin(), BIL_hits.end());
        track.insert(track.end(), BIS_hits.begin(), BIS_hits.end());
    }
    // If there is overlap only use hits from the BIL
    else if (status.at(0) == 1) {
        track.insert(track.end(), BIL_hits.begin(), BIL_hits.end());
    }
    //If there is no overlap in the middle region. Add both chamber types hits (one will be empty)
    if (status.at(1) == 0) {
        track.insert(track.end(), BML_hits.begin(), BML_hits.end());
        track.insert(track.end(), BMS_hits.begin(), BMS_hits.end());
    } 
    // If there is overlap only use hits from the BML
    else if (status.at(1) == 1) {
        track.insert(track.end(), BML_hits.begin(), BML_hits.end());
    }
    if (status.at(2) == 0) {
        track.insert(track.end(), BOL_hits.begin(), BOL_hits.end());
        track.insert(track.end(), BOS_hits.begin(), BOS_hits.end());
    } 
    // If there is overlap only use hits from the BOL
    else if (status.at(2) == 1) {
        track.insert(track.end(), BOL_hits.begin(), BOL_hits.end());
    }
    return track;
}

// A function to build a tracklet hit vector accounting for overlap regions
std::vector<float> buildTracklets(std::map<int, std::vector<float>> hitMap, std::set<int> stations) {
    std::vector<float> tracklet;
    std::vector<float> BIL_hits = hitMap[0];
    std::vector<float> BIS_hits = hitMap[1];
    std::vector<float> BML_hits = hitMap[2];
    std::vector<float> BMS_hits = hitMap[3];
    // Get the overlap statuses in the inner and middle layers
    std::vector<int> status = overlaps(stations);
    //If there is no overlap in the inner region. Add both chamber types hits (one will be empty)
    if (status.at(0) == 0) {
        tracklet.insert(tracklet.end(), BIL_hits.begin(), BIL_hits.end());
        tracklet.insert(tracklet.end(), BIS_hits.begin(), BIS_hits.end());
    }
    // If there is overlap only use hits from the BIL
    else if (status.at(0) == 1) {
        tracklet.insert(tracklet.end(), BIL_hits.begin(), BIL_hits.end());
    }
    //If there is no overlap in the middle region. Add both chamber types hits (one will be empty)
    if (status.at(1) == 0) {
        tracklet.insert(tracklet.end(), BML_hits.begin(), BML_hits.end());
        tracklet.insert(tracklet.end(), BMS_hits.begin(), BMS_hits.end());
    } 
    // If there is overlap only use hits from the BML
    else if (status.at(1) == 1) {
        tracklet.insert(tracklet.end(), BML_hits.begin(), BML_hits.end());
    }
    return tracklet;
}

// Build a vector of bin centers
std::vector<float> buildBinCenters(double start, double stop, double binSize) {
    double binRange = stop - start;
    int numBins = static_cast<int>(binRange / binSize);
    // See if the binSize divides the bin range
    if (binRange - (numBins * binSize) > 1e-6) { // Check for floating-point precision
        std::cout << "bin size must divide bin range" << std::endl;
        return {};
    }
    // Build the bin center vector
    std::vector<float> binCenters;
    for (int i = 1; i <= numBins; i++) {
        binCenters.push_back(static_cast<float>(start + ((i - 0.5) * binSize)));
    }
    return binCenters;
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

// An overload function to create the x% truncated mean and its error from a vector of ADC counts
std::tuple<float, float> TM(std::vector<float> adcs, float truncation_percent, const char* option) {
    std::set<std::string> options = {"e","E"};
    if (options.find(option) == options.end()) {
        std::cout << "Invalid option" << std::endl;
        exit(EXIT_FAILURE); 
    }
    std::vector<int> new_adcs;
    float size = adcs.size();
    float truncated_size = size * (1 - truncation_percent);
    float points_truncated = size - truncated_size;

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
        new_adcs.push_back(adcs.at(i));
    }

    // Add replace the truncated points with the greatest untruncated value
    // This will be used to compute winsorized mean stabdard error
    new_adcs.insert(new_adcs.end(), points_truncated, adcs.at(truncated_size));
    float winsorized_mean_stderror = TMath::StdDev(new_adcs.begin(), new_adcs.end());
    float truncated_mean_stderror = winsorized_mean_stderror / ((1 - truncation_percent) * sqrt(size));
    // Calculate truncated mean
    float truncated_mean = track_sum / truncated_size;
    // Then return the truncated mean and its error
    return (std::make_tuple(truncated_mean, truncated_mean_stderror));
}

// Check if a muon is a CB or standalone muon
int checkAuthor(unsigned short author) {
    //if (author == 1 || author == 5) {
    if (author == 1) {
        return 1;
    }
    else {
        return 0;
    }
}

// Check if the tube hit is in the valid radius
bool isInValidRadius(float radius) {
    if (radius >= 1 && radius <= 14) {
        return true;
    }
    return false;
}

// Return the Landau fit parameters and errors of each xBin in a TH2 object
std::vector<std::tuple<float, float>> fit2DHistogramLandau(TH2 *h2) {
    int nbinsx = h2->GetNbinsX();
    std::vector<std::tuple<float, float>> fitParameters;

    for (int i = 1; i <= nbinsx; i++) {
        TH1D* projectionY = h2->ProjectionY("", i, i);
        TF1* landauFit = new TF1("landauFit", "landau", projectionY->GetXaxis()->GetXmin(), projectionY->GetXaxis()->GetXmax());

        // Fit the projection to the Landau function
        projectionY->Fit(landauFit, "RQ");

        // Retrieve the fit parameters
        float mpv = landauFit->GetParameter(1);
        float mpvErr = landauFit->GetParError(1);
        std::tuple<float, float> binFit = std::make_tuple(mpv, mpvErr);

        fitParameters.push_back(binFit);

        delete projectionY;
        delete landauFit;
    }
    return fitParameters;
}

// Return a Truncated mean from each xBin in a TH2 object
std::vector<std::tuple<float, float>> fit2DHistogramTM(TH2 *h2) {
    int nbinsx = h2->GetNbinsX();

    //Record the truncate dmean in each bin
    std::vector<std::tuple<float, float>> truncated_means;
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
            truncated_means.push_back(std::make_tuple(0, 10000));
            delete projectionY;
        }
        else {
            // Calculate the truncated mean and uncertainty for that bin
            std::tuple<float, float> truncated_mean = TM(temp_bin, 0.4, "E");
            truncated_means.push_back(truncated_mean);
            delete projectionY;
        }
    }
    return truncated_means;
}

// Fit a track(let) to a Landau distribution for a dE/dx estimate
float fitTrackLandau(std::vector<float> track) {
    // Create a histogram and fit for the track and check if it exists first
    TH1* existing_h = (TH1*)gDirectory->Get("hlan");
    if (existing_h) {
        delete existing_h;
    }

    TH1F *hlan = new TH1F("hlan", "hlan", 100,-5,5);
    TF1 *f = new TF1("f", "landau", -5,5);
    // Set initial parameter values for the Landau function
    f->SetParameters(2, 1, 0.1);
    for (auto entry: track) {
        hlan->Fill(entry);
    }
    hlan->Fit(f, "QR");

    // Retrieve the fit parameters
    float mpv = f->GetParameter(1);
    return mpv;
}

// Give the median of a TH1 object
Double_t median1D(TH1 *h1) {
    Int_t num_bins = h1->GetXaxis()->GetNbins();
    Double_t *x = new Double_t[num_bins];
    Double_t *y = new Double_t[num_bins];
    for (Int_t i=0; i<num_bins; i++) {
        x[i] = h1->GetXaxis()->GetBinCenter(i+1);
        y[i] = h1->GetBinContent(i+1);
    }
    Double_t median = TMath::Median(num_bins, x, y);
    delete [] x;
    delete [] y;
    return median;
}

// Find the median for each slice along x axis of a TH2 object h2
std::vector<float> median2D(TH2 *h2) {
    std::vector<float> medians;
    Int_t num_bins = h2->GetXaxis()->GetNbins();
    for (Int_t i=1; i<=num_bins; i++) {
        TH1 *h1 = h2->ProjectionY("",i,i);
        float median = median1D(h1);
        std::vector<float>::iterator it = medians.end();
        medians.insert(it, median);
    }
    return medians;
}

// Define a general nth degree polynomial
float polynomial(int degree, float drift_radius, const Double_t *params) {
    float value = 0;
    for (int i = 0; i <= degree; i++) {
        value += params[i] * pow(drift_radius, i);
    }
    return value;
}

// An overloaded polynomial function which takes a vector of parameters rather than an array
float polynomial(int degree, float drift_radius, std::vector<Double_t> params) {
    float value = 0;
    for (int i = 0; i <= degree; i++) {
        value += params.at(i) * pow(drift_radius, i);
    }
    return value;
}

// A function to create a TPaveText label for a polynomial fit
TPaveText* pol_fit_label(int degree, TF1 *f) {
    // Retrieve the fit results
    const Double_t *parameters = f->GetParameters();
    const Double_t *parErrors = f->GetParErrors();
    
    // Create the label and Add the relevant info (fit parameters, errors, chi^2)
    TPaveText* fit_label = new TPaveText(0.7,0.1,0.9,0.3,"NDC");
    fit_label->SetTextAlign(12);
    fit_label->AddText(Form("Fit Parameters"));
    for (int i = 0; i <= degree; i++) {
        fit_label->AddText(Form("p%d: %f +/- %f", i, parameters[i], parErrors[i]));
    }
    fit_label->AddText(Form("Chi-Squared: %.4f", f->GetChisquare()));
    fit_label->AddText(Form("NDF: %d", f->GetNDF()));

    return fit_label;
}

// A function to create a TPaveText label for a gaussian fit
TPaveText* gaus_fit_label(TF1 *f) {
    // Retrieve the fit results
    TPaveText *fit_label = new TPaveText(0.777,0.753,0.977,0.953,"NDC");
    fit_label->SetTextSize(0.05);
    fit_label->AddText(Form("Chi^2 / NDF: %f", f->GetChisquare() / f->GetNDF()));
    fit_label->AddText(Form("Scale: %f +/- %f", f->GetParameter(0), f->GetParError(0)));
    fit_label->AddText(Form("Mean: %f +/- %f", f->GetParameter(1), f->GetParError(1)));
    fit_label->AddText(Form("Std dev: %f +/- %f", f->GetParameter(2), f->GetParError(2)));
    fit_label->AddText(Form("Relative Res: %f", f->GetParameter(2) / f->GetParameter(1)));
    return fit_label;
}

// A function to create a TPaveText label for a fit of two gaussians
TPaveText* double_gaus_fit_label(TF1 *f) {
    // Retrieve the fit results
    TPaveText *fit_label = new TPaveText(0.777,0.753,0.977,0.953,"NDC");
    fit_label->SetTextSize(0.05);
    fit_label->AddText(Form("Chi^2 / NDF: %f", f->GetChisquare() / f->GetNDF()));
    fit_label->AddText(Form("Mean: %f +/- %f", f->GetParameter(1), f->GetParError(1)));
    fit_label->AddText(Form("1st Peak Scale: %f +/- %f", f->GetParameter(0), f->GetParError(0)));
    fit_label->AddText(Form("1st Peak Std dev: %f +/- %f", f->GetParameter(2), f->GetParError(2)));
    fit_label->AddText(Form("1st Peak Relative Res: %f", f->GetParameter(2) / f->GetParameter(1)));
    fit_label->AddText(Form("2nd Peak Scale: %f +/- %f", f->GetParameter(3), f->GetParError(3)));
    fit_label->AddText(Form("2nd Peak Std dev: %f +/- %f", f->GetParameter(4), f->GetParError(4)));
    fit_label->AddText(Form("2nd Peak Relative Res: %f", f->GetParameter(4) / f->GetParameter(1)));
    return fit_label;
}

// Builds a set which holds each MDT chamber in the inner, middle and outer barrel as
// a tuple in the format (stationIndex, stationEta, stationPhi)
// Define possible chamber indices
std::map<std::tuple<std::string, int, int>, std::string> buildChambers() { 
    std::map<std::tuple<std::string, int, int>, std::string> chamberDict;

    // Define possible values for chamberIndex, chamberEta, and chamberPhi
    std::vector<std::string> chamberIndices = {"BIL", "BIS", "BML", "BMS", "BOL", "BOS"};
    std::vector<int> etaRange, phiRange;

    // Populate eta and phi ranges
    for (int i = -8; i <= 8; ++i) {
        etaRange.emplace_back(i);
    }
    for (int i = 1; i <= 8; ++i) {
        phiRange.emplace_back(i);
    }
    
    // Generate all possible chambers
    for (const auto& index : chamberIndices) {
        for (const auto& eta : etaRange) {
            for (const auto& phi : phiRange) {
                std::tuple<std::string, int, int> chamberTuple = std::make_tuple(index, eta, phi);
                std::string chamberInfo = index + " / Eta:" + std::to_string(eta) + " / Phi:" + std::to_string(phi);
                chamberDict[chamberTuple] = chamberInfo;
            }
        }
    }
    return chamberDict;
}

// Converts the stationIndex to a text string from whatever hexadecmial value it has
std::map<std::string, std:: string> stationTranslator() {
    std::map<std::string, std:: string> stations;
    stations["0"] = "BIL";
    stations["1"] = "BIS";
    stations["2"] = "BML";
    stations["3"] = "BMS";
    stations["4"] = "BOL";
    stations["5"] = "BOS";
    return stations;
}

// A function for standard deviatoin of a vector
float standardDeviation(std::vector<float> input) {
    float mean = std::accumulate(input.begin(), input.end(), 0.0) / input.size();
    float stdev = 0;
    for (unsigned int i=0; i<input.size(); i++) {
        stdev += pow(input.at(i) - mean, 2);
    }
    stdev = stdev / input.size();
    stdev = pow(stdev, 0.5);
    return stdev;
}

// A function to create the MPV estimator vs. Beta*gamma plot with log scale points
TGraphErrors *BetheBlochEstimatorCurve(TH1F* Momentum_distribution, vector<float> TrackMomenta, vector<float> TrackTMs, float start_p, float end_p, int BB_num_bins) {
    
    // Make evenly populated momentum bins from the Momentum distribution
    float entries_per_bin;
    for (int i=floor(start_p); i<floor(end_p); i++) {
        entries_per_bin += Momentum_distribution->GetBinContent(i); 
    }
    entries_per_bin /= BB_num_bins;

    // Find the edges of the bins
    vector<float> BB_bin_centers;
    vector<float> BB_bin_edges;
    float current_entries = 0;
    float last_edge = start_p;
    for (int i=floor(start_p); i<floor(end_p); i++) {
        current_entries+=Momentum_distribution->GetBinContent(i);
        if (current_entries > entries_per_bin) {
            BB_bin_edges.push_back(i);
            BB_bin_centers.push_back((i + last_edge) / 2);
            current_entries = 0;
            last_edge = i;
        }
    }
    if (BB_bin_centers.size() < BB_num_bins) {
        BB_bin_edges.push_back(end_p);
        BB_bin_centers.push_back((last_edge + end_p ) / 2);
    }

    // Vector structure is:   <p bin<Momenta<Track Instance>, TMs<Track Instance>>>
    std::vector<std::pair<std::vector<float>, std::vector<float>>> bin_TMs(10, {std::vector<float>(), std::vector<float>()});
    //Loop through the tracks
    for (unsigned int i=0; i<TrackMomenta.size(); i++) {
        if (TrackMomenta.at(i) > end_p || TrackMomenta.at(i) < start_p) { continue; }

        // Place the track into the right p bin
        for (unsigned int n=0; n<BB_bin_edges.size(); n++) {
            if (TrackMomenta.at(i) < BB_bin_edges.at(n)) {
                bin_TMs.at(n).first.push_back(TrackMomenta.at(i));
                bin_TMs.at(n).second.push_back(TrackTMs.at(i));
                break;
            }
        }
    }

    // Populate the TGraph
    TGraphErrors *BBCurve = new TGraphErrors();
    Double_t bin_mean_p;
    Double_t bin_mean_p_uncertainty;
    Double_t bin_mean_TM;
    Double_t bin_mean_TM_uncertainty;
    for (unsigned long i=0; i<BB_bin_centers.size(); i++) {
        bin_mean_p = std::accumulate(bin_TMs.at(i).first.begin(), bin_TMs.at(i).first.end(), 0.0) / bin_TMs.at(i).first.size();
        bin_mean_p_uncertainty = standardDeviation(bin_TMs.at(i).first);
        bin_mean_TM = std::accumulate(bin_TMs.at(i).second.begin(), bin_TMs.at(i).second.end(), 0.0) / bin_TMs.at(i).second.size();
        bin_mean_TM_uncertainty = standardDeviation(bin_TMs.at(i).second);
        //BBCurve->AddPointError(bin_mean_p, bin_mean_TM, bin_mean_p_uncertainty, bin_mean_TM_uncertainty);
        BBCurve->AddPoint(bin_mean_p, bin_mean_TM);
        BBCurve->SetPointError(i, bin_mean_p_uncertainty, bin_mean_TM_uncertainty);
    }
    return BBCurve;
}



/*
* Primary Operations: Perform calibrations to account for ADC count dependency on drift radius
* Plot these fits, their paramters and the data from which they were constructed for each chamber
* Then apply these corrections. Show the global effect of the corrections
* And finally produce 20, 40% truncated means of the corrected ADC counts for muon tracklets
*/
void adcVsRadius(std::string indir, const::std::string& outfile, const:: std::string& chamber_level_outfile) {

    // The number of xbins to be used in TH2 objects later
    const float numbins = 30;
    const float zmass = 91;

    // The bins of momentum we want to consider muons for
    float low_momentum_edge = 23;
    float high_momentum_edge = 70;
    float low_interval = 3;
    float high_interval = 15;
    std::vector<std::pair<Double_t, Double_t>> momentum_bins = {
        {low_momentum_edge - low_interval, low_momentum_edge + low_interval},
        {high_momentum_edge - high_interval, high_momentum_edge + high_interval}
    };
    
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
    chain.SetBranchStatus("trkMdt_HitType", true);
    chain.SetBranchStatus("muons_author", true);
    chain.SetBranchStatus("muons_pt", true);
    chain.SetBranchStatus("muons_eta", true);
    chain.SetBranchStatus("muons_phi", true);

    // Define variables to hold branch values
    std::vector<float> *eta = nullptr;
    std::vector<int> *ADC_counts = nullptr;
    std::vector<unsigned short> *trk_link = nullptr;
    std::vector<float> *drift_radius = nullptr;
    std::vector<unsigned char> *station_index = nullptr;
    std::vector<char> *station_eta = nullptr;
    std::vector<unsigned char> *station_phi = nullptr;
    std::vector<short> *hit_type = nullptr;
    std::vector<unsigned short> *authors = nullptr;
    std::vector<float> *pts = nullptr;
    std::vector<float> *etas = nullptr;
    std::vector<float> *phis = nullptr;

    // Set the address to each branch in the chain
    chain.SetBranchAddress("muons_eta", &eta);
    chain.SetBranchAddress("trkMdt_Adc", &ADC_counts);
    chain.SetBranchAddress("trkMdt_muonsLink", &trk_link);
    chain.SetBranchAddress("trkMdt_driftR", &drift_radius);
    chain.SetBranchAddress("trkMdt_stationIndex", &station_index);
    chain.SetBranchAddress("trkMdt_stationEta", &station_eta);
    chain.SetBranchAddress("trkMdt_stationPhi", &station_phi);
    chain.SetBranchAddress("trkMdt_HitType", &hit_type);
    chain.SetBranchAddress("muons_author", &authors);
    chain.SetBranchAddress("muons_pt", &pts);
    chain.SetBranchAddress("muons_eta", &etas);
    chain.SetBranchAddress("muons_phi", &phis);

    // Create the map from hexadecimal labels to decimal label strings for the stationIndex
    std::map<std::string, std::string> stationStrings = stationTranslator();

    // Create the map from chamber info to the chamber label
    std::map<std::tuple<std::string, int, int>, std::string> chamber_info = buildChambers();

    // Create a map which links chamber labels to an array containing their ADC count
    // and drift radii. Ex. map[BOS_1_3] = (<10.6, 4.2,...>, <124, 142,...>)
    std::map<std::string, std::array<std::tuple<std::vector<float>, std::vector<float>>, 2>> chamber_data;
    for (const auto& entry: chamber_info) {
        const std::string& key = entry.second;
        chamber_data[key] = {std::make_tuple(std::vector<float>(), std::vector<float>()),
                             std::make_tuple(std::vector<float>(), std::vector<float>())};
    }


    // Initialize the 2D histogram and the TGraph scatter plot
    std::string lowbin_title = "Global ADC vs. Drift Radius (" 
    + to_string(low_momentum_edge - low_interval) + "-" + to_string(low_momentum_edge + low_interval) + ")GeV; Drift Radius (mm); ADC Count";
    std::string highbin_title = "Global ADC vs. Drift Radius (" 
    + to_string(high_momentum_edge - high_interval) + "-" + to_string(high_momentum_edge + high_interval) + ")GeV; Drift Radius (mm); ADC Count";

    TH2F* h = new TH2F("h", "ADC vs. Drift Radius (MuonTesterRun456714); Drift Radius (mm); ADC Count", numbins, 0, 15, 400, 0, 400);
    TH2F* hglobal_lowbin = new TH2F("hglobal_lowp", lowbin_title.c_str(), numbins, 0, 15, 400, 0, 400);
    TH2F* hglobal_highbin = new TH2F("hglobal_highp", highbin_title.c_str(), numbins, 0, 15, 400, 0, 400);
    TGraph *g1 = new TGraph();
    g1->SetMarkerStyle(5);

    // Create a histogram to visualize the MDT hit distribution with drift radius
    TH1F* hradius = new TH1F("hradius", "Tube Hits vs. Drift Radius (MuonTesterRun456714); Drift Radius(mm); Tube Hits", numbins, 0, 15);

    // Read the relevant branches for all entries
    Long64_t nEntries = chain.GetEntries();
    for (Int_t i=0; i < nEntries; i++) {
        // A progress indicator
        
        if (i % 10000 == 0) {
            std::cout << i << "/" << nEntries << " processed" << endl;
        }
        
        // Get the entry data
        chain.GetEntry(i);

        // Make sure the current muon has some ADC counts recorded
        if (ADC_counts->size() == 0) {
            continue;
        }

        // Initialize the muon link
        int muon_link = 0;
        int author = authors->at(muon_link);
        float muon_pt = pts->at(muon_link);
        float muon_eta = etas->at(muon_link);
        
        // Loop through all hits in the entry
        for (unsigned long n = 0; n < ADC_counts->size(); n++) {
            // Get the index station as a string and translate it to the text label
            // If it isn't in the barrel skip it
            std::string parsed_station_index = std::to_string(station_index->at(n));
            std::string translated_station_index;
            if (stationStrings.find(parsed_station_index) != stationStrings.end()) {
                translated_station_index = stationStrings[parsed_station_index];
            } else {
                continue;
            }

            // Parse the other parts of the chmaber ID and form it into a chamber string
            // which will identify which histogram to put the hit into
            int stationEta = station_eta->at(n);
            int stationPhi = station_phi->at(n);
            std::string chamber_label = chamber_info[std::make_tuple(translated_station_index, stationEta, stationPhi)];
            
            /*
            // Skip over sMDT data (Don't for simulated data)
            std::string sMDT = "BIS / Eta:7";
            if (chamber_label.rfind(sMDT, 0) == 0) {
                continue;
            }
            */
            
            // Place the hit into the correct momentum range
            int momentum_range = -1;
            if (ptToP(muon_pt, muon_eta) >= momentum_bins.at(0).first && ptToP(muon_pt, muon_eta) < momentum_bins.at(0).second) {
                momentum_range = 0;
            } else if (ptToP(muon_pt, muon_eta) >= momentum_bins.at(1).first && ptToP(muon_pt, muon_eta) < momentum_bins.at(1).second) {
                momentum_range = 1;
            }

            // Check if we have a new muon
            if (trk_link->at(n) != muon_link) {
                // Update the muon link
                muon_link = trk_link->at(n);
                author = authors->at(muon_link);
                muon_pt = pts->at(muon_link);
                muon_eta = etas->at(muon_link);

                // Cut on hits where ADC counts < 50, author, eta and pt
                if (ADC_counts->at(n) < 50 || hit_type->at(n) > 60  || checkAuthor(author) == 0 || muon_pt < 5 || abs(muon_eta) > 1 || momentum_range == -1) {
                    continue;
                }
                else {
                    // Fill the histogram and scatter plot if the above conditions are satisfied
                    h->Fill(drift_radius->at(n), ADC_counts->at(n));
                    // If in the low momentum bin, fill the lowbin histogram and the high one if in the high range
                    (momentum_range & 1 ? hglobal_highbin : hglobal_lowbin)->Fill(drift_radius->at(n), ADC_counts->at(n));
                    g1->AddPoint(drift_radius->at(n), ADC_counts->at(n));
                    // Add the data to the end of the data vectors in the tuple for the appropriate chamber
                    //std::get<0>(chamber_data[chamber_label]).emplace_back(drift_radius->at(n));
                    //std::get<1>(chamber_data[chamber_label]).emplace_back(ADC_counts->at(n));
                    auto& [drift_radii, adc_counts] = chamber_data[chamber_label][momentum_range];
                    drift_radii.emplace_back(drift_radius->at(n));
                    adc_counts.emplace_back(ADC_counts->at(n));

                    // Also fill the tube hits vs. Drift Radius histogram
                    hradius->Fill(drift_radius->at(n));
                }
            // Now check the same conditions for when it is not a new muon
            } else { 
                if (ADC_counts->at(n) < 50 || checkAuthor(author) == 0 || muon_pt < 5 || abs(muon_eta) > 1 || momentum_range == -1) { continue; } 
                else {
                h->Fill(drift_radius->at(n), ADC_counts->at(n));
                // If in the low momentum bin, fill the lowbin histogram and the high one if in the high range
                (momentum_range & 1 ? hglobal_highbin : hglobal_lowbin)->Fill(drift_radius->at(n), ADC_counts->at(n));
                g1->AddPoint(drift_radius->at(n), ADC_counts->at(n));
                // Add the data to the end of the data vectors in the tuple for the appropriate chamber
                //std::get<0>(chamber_data[chamber_label]).emplace_back(drift_radius->at(n));
                //std::get<1>(chamber_data[chamber_label]).emplace_back(ADC_counts->at(n));
                auto& [drift_radii, adc_counts] = chamber_data[chamber_label][momentum_range];
                drift_radii.emplace_back(drift_radius->at(n));
                adc_counts.emplace_back(ADC_counts->at(n));

                // Also fill the tube hits vs. Drift Radius histogram
                hradius->Fill(drift_radius->at(n));
                }
            }
        }
    }

    // Fit a Landau distribution to each drift radius bin and return the result here
    std::vector<std::tuple<float, float>> globalFits = fit2DHistogramTM(h);
    std::vector<std::tuple<float, float>> lowp_globalFits = fit2DHistogramTM(hglobal_lowbin);
    std::vector<std::tuple<float, float>> highp_globalFits = fit2DHistogramTM(hglobal_highbin);

    // Create the bin centers vector(mm / bins)
    float bin_size = 15 / numbins;
    std::vector<float> v_bin_centers = buildBinCenters(0, 15, bin_size);

    // We want to exclude the first and last mm of drift radius bins
    int excluded_bins = round(1 / bin_size);

    // Initialize the arrays for the bin MPVs, their errors and their drift radius
    int good_bins = v_bin_centers.size() - (2 * excluded_bins);
    Double_t globalMPVs[good_bins];
    Double_t globalMPVErrs[good_bins];
    Double_t lowp_globalMPVs[good_bins];
    Double_t lowp_globalMPVErrs[good_bins];
    Double_t highp_globalMPVs[good_bins];
    Double_t highp_globalMPVErrs[good_bins];
    Double_t bin_centers[good_bins];
    Double_t bin_centers_errs[good_bins];

    // Loop through the drift radius bins we're concerned with and add to bin_centers and MPVs
    for (unsigned long j = excluded_bins; j < v_bin_centers.size() - excluded_bins; j++) {
        globalMPVs[j - excluded_bins] = std::get<0>(globalFits.at(j));
        globalMPVErrs[j - excluded_bins] = std::get<1>(globalFits.at(j));
        lowp_globalMPVs[j - excluded_bins] = std::get<0>(lowp_globalFits.at(j));
        lowp_globalMPVErrs[j - excluded_bins] = std::get<1>(lowp_globalFits.at(j));
        highp_globalMPVs[j - excluded_bins] = std::get<0>(highp_globalFits.at(j));
        highp_globalMPVErrs[j - excluded_bins] = std::get<1>(highp_globalFits.at(j));
        bin_centers[j - excluded_bins] = v_bin_centers.at(j);
        bin_centers_errs[j - excluded_bins] = 0;
    }

    // Construct a second TGraph to which we will fit a polynomial from the relevant MPVs
    TGraphErrors *g2 = new TGraphErrors(good_bins, bin_centers, globalMPVs, bin_centers_errs, globalMPVErrs);
    TGraphErrors *g2_low = new TGraphErrors(good_bins, bin_centers, lowp_globalMPVs, bin_centers_errs, lowp_globalMPVErrs);
    TGraphErrors *g2_high = new TGraphErrors(good_bins, bin_centers, highp_globalMPVs, bin_centers_errs, highp_globalMPVErrs);


    // Set limits on the x and y ranges for the TGraph and make the marker style
    g2->GetXaxis()->SetLimits(0,15);
    g2->SetMinimum(0);
    g2->SetMaximum(160);
    g2->SetTitle("ADC MPV vs. Drift Radius (MuonTesterRun456714) | (pol4); Drift Radius (mm); ADC Count MPV");
    g2->SetMarkerStyle(7);
    
    g2_low->GetXaxis()->SetLimits(0,15);
    g2_low->SetMinimum(0);
    g2_low->SetMaximum(160);
    g2_low->SetTitle(lowbin_title.c_str());
    g2_low->SetMarkerStyle(7);

    g2_high->GetXaxis()->SetLimits(0,15);
    g2_high->SetMinimum(0);
    g2_high->SetMaximum(160);
    g2_high->SetTitle(highbin_title.c_str());
    g2_high->SetMarkerStyle(7);

    
    // Fit a 4th degree polynomial to the global MPV data
    TF1* fglobalpol4 = new TF1("fglobalpol4","pol4", 1,14);
    g2->Fit(fglobalpol4, "RQ");
    TF1* fglobalpol4low = new TF1("fglobalpol4low","pol4", 1,14);
    g2_low->Fit(fglobalpol4low, "RQ");
    TF1* fglobalpol4high = new TF1("fglobalpol4high","pol4", 1,14);
    g2_high->Fit(fglobalpol4high, "RQ");
    g2_low->Write();
    g2_high->Write();

    // retrieve the global polynomial fit parameters for calibrating ADC counts later
    const Double_t *globalPol4Params = fglobalpol4->GetParameters();
    
    // Retrieve the fits for the different momentum bins
    const Double_t *lowGlobalPol4Params = fglobalpol4low->GetParameters();
    const Double_t *highGlobalPol4Params = fglobalpol4high->GetParameters();

    // Initialize vectors which will holds the fit parameters for each chamber
    std::tuple<std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>> chamber_fits;

    // Open the chamber level TFile
    TFile *chamber_level_file = new TFile(chamber_level_outfile.c_str(), "RECREATE");

    // Write the tube hits vs. Drift radius histogram
    hradius->Write();

    // Create a map from a chamber label to momentum bin to polynomial fit
    std::map<std::string, std::map<int, std::vector<Double_t>>> chamber_polynomials;

    // Create a vector to track the number of hits per chamber used for calibration
    std::vector<int> hits_per_chamber;
    // Create a Histogram which will hold all of the Chi-Squared/NDF statistics of each chamber
    TH1F *hchamber_chis = new TH1F("hchamber_chis", "Chi-Squared by Chamber (MuonTesterRun456714); Chi-squared; Chamber", 20, 0, 10);

    // CReate Graphs for the distribution of MPV vs drift radius in each momnetum bins
    TH2F *h_low_chamberlevel = new TH2F("h_low_chamberlevel", "MPV Distribution by Drift Radius (20-26GeV Muons); Drift Radius(mm); MPV", 30, 0, 15, 400, 0, 400);
    TH2F *h_high_chamberlevel = new TH2F("h_high_chamberlevel", "MPV Distribution by Drift Radius (55-85GeV Muons); Drift Radius(mm); MPV", 30, 0, 15, 400, 0, 400);
    TH3F *h3_low_chamberlevel = new TH3F("h3_low_chamberlevel", "MPV Distribution by Drift Radius (20-26GeV Muons); Drift Radius(mm); MPV; Bin Hits", 30, 0, 15, 400, 0, 400, 20, 0, 400);
    TH3F *h3_high_chamberlevel = new TH3F("h3_high_chamberlevel", "MPV Distribution by Drift Radius (55-85GeV Muons); Drift Radius(mm); MPV; Bin Hits", 30, 0, 15, 400, 0, 400, 20, 0, 400);
    

    // Now Create a TGraph and TH2F object for every chamber and populate them with the corresponding chamber data
    // These TGraphs will fit a landau distribution to each drift radius bin in the TH2F 
    // This collection of points will then be fit to a polynomial

    for (const auto& entry: chamber_data) {
        // Get the chamber label
        const std::string& key = entry.first;

        // Go through the momentum bins
        for (int p_bin =0; p_bin < 2; p_bin++) {
            // Retrieve the vectors for drift radius and adc count in the curretn momentum bin and chamber
            const auto&[drift_radius_data, ADC_count_data] = chamber_data[key][p_bin];
        
            // retrieve the Vectors of data for the current chamber and cast them as TVectorF objects
            // vector<float> drift_radius_data = std::get<0>(chamber_data[key]);
            // vector<float> ADC_count_data = std::get<1>(chamber_data[key]);

            // Check if there are any ADC counts in the chamber
            if (ADC_count_data.size() == 0) { continue; }

            // Record the number of hits to be used for this calibration
            hits_per_chamber.emplace_back(drift_radius_data.size());

            // Title the plot for this chamber and initialize the TH2F
            std::string graph_title_with_axes = key + " (Bin: " + to_string(static_cast<int>(momentum_bins.at(p_bin).first)) + "-" +\
             to_string(static_cast<int>(momentum_bins.at(p_bin).second)) + "GeV/c);Drift Radius(mm); ADC counts MPV";
            TH2F* hist = new TH2F((key + "_bin" + to_string(p_bin)).c_str(), graph_title_with_axes.c_str(), numbins, 0, 15, 400, 0, 400);

            // Add the vectors to the histogram
            for (long unsigned hit_i = 0; hit_i < drift_radius_data.size(); ++hit_i) {
                hist->Fill(drift_radius_data.at(hit_i), ADC_count_data.at(hit_i));
            }

            // Fit a Landau distribution to each bin and return the result here
            std::vector<std::tuple<float, float>> chamberFitParams = fit2DHistogramTM(hist);
            
            // Get the number of htis in each drift bin
            std::vector<float> drift_bin_sizes (v_bin_centers.size(), 0);
            int n_ybins = hist->GetNbinsY();
            for (unsigned long xbin = 0; xbin < v_bin_centers.size(); xbin++) {
                for (int ybin=1; ybin <=n_ybins; ybin++) {
                    drift_bin_sizes.at(xbin) += hist->GetBinContent(xbin, ybin);
                }
            }

            // Loop through the MPVs in the interesting drift radius domain and add to bin_centers 
            // with their corresponding MPV
            // We want to exclude the first two bins and the last two bins
            Double_t chamberMPVs[good_bins];
            Double_t chamberMPVErrs[good_bins];
            for (unsigned long j = excluded_bins; j < v_bin_centers.size() - excluded_bins; j++) {
                chamberMPVs[j - excluded_bins] = std::get<0>(chamberFitParams.at(j));
                chamberMPVErrs[j - excluded_bins] = std::get<1>(chamberFitParams.at(j));
                // Fill the MPV vs Drift radius distribution histogram for the correct momentum bin
                p_bin == 0 ? h_low_chamberlevel->Fill(v_bin_centers.at(j), chamberMPVs[j - excluded_bins]) :
                h_high_chamberlevel->Fill(v_bin_centers.at(j), chamberMPVs[j - excluded_bins]);
                p_bin == 0 ? h3_low_chamberlevel->Fill(v_bin_centers.at(j), chamberMPVs[j - excluded_bins], drift_bin_sizes.at(j)) :
                h3_high_chamberlevel->Fill(v_bin_centers.at(j), chamberMPVs[j - excluded_bins], drift_bin_sizes.at(j));
            }

            // Loop through the drift radius bins we're concerned with and add to bin_centers and MPVs
            for (unsigned long j = excluded_bins; j < v_bin_centers.size() - excluded_bins; j++) {
                globalMPVs[j - excluded_bins] = std::get<0>(globalFits.at(j));
                globalMPVErrs[j - excluded_bins] = std::get<1>(globalFits.at(j));
                lowp_globalMPVs[j - excluded_bins] = std::get<0>(lowp_globalFits.at(j));
                lowp_globalMPVErrs[j - excluded_bins] = std::get<1>(lowp_globalFits.at(j));
                highp_globalMPVs[j - excluded_bins] = std::get<0>(highp_globalFits.at(j));
                highp_globalMPVErrs[j - excluded_bins] = std::get<1>(highp_globalFits.at(j));
                bin_centers[j - excluded_bins] = v_bin_centers.at(j);
                bin_centers_errs[j - excluded_bins] = 0;
                
            }

            // Construct the TGraph from the arrays created above
            TGraphErrors *g = new TGraphErrors(good_bins, bin_centers, chamberMPVs, bin_centers_errs, chamberMPVErrs);
            g->SetNameTitle(key.c_str(), graph_title_with_axes.c_str());
            g->SetMarkerStyle(7);
            g->SetMinimum(0);

            // Now perform the fit to these MPV plots for each chamber and overlay it on the TH2D
            TCanvas *c2chamber = new TCanvas((key + "_bin" + std::to_string(p_bin)).c_str(), graph_title_with_axes.c_str(), 800, 600);
            TF1* fchamber = new TF1("fchamber","pol4", 1,14);
            g->Fit(fchamber, "RQ");
            TF1 *fittedChamberFunc = g->GetFunction("fchamber");
            TPaveText* chamber_fit_label = pol_fit_label(4, fchamber);
            hist->Draw("colz");
            g->Draw("P same");
            fittedChamberFunc->Draw("same");
            chamber_fit_label->Draw("same");
            c2chamber->Write();

            std::get<0>(chamber_fits).push_back(fchamber->GetParameter(0));
            std::get<1>(chamber_fits).push_back(fchamber->GetParameter(1));
            std::get<2>(chamber_fits).push_back(fchamber->GetParameter(2));
            std::get<3>(chamber_fits).push_back(fchamber->GetParameter(3));
            std::get<4>(chamber_fits).push_back(fchamber->GetParameter(4));
            
            // Create an std::vector of polynomial parameters to map to the chamber label
            std::vector<Double_t> chamber_params = {fchamber->GetParameter(0), fchamber->GetParameter(1), fchamber->GetParameter(2), \
            fchamber->GetParameter(3), fchamber->GetParameter(4)};
            chamber_polynomials[key][p_bin] = chamber_params;

            // Create a canvas so we can draw the TGraph with the preferred options and write
            // the canvas to the chamber level fileß
            TCanvas *cchamber = new TCanvas((key + "_bin" + std::to_string(p_bin)).c_str(), graph_title_with_axes.c_str(), 800, 600);
            cchamber->cd();
            g->Draw("AP");
            chamber_fit_label->Draw();
            // cchamber->Write();

            // Extract the goodness of fit to display it for all chambers
            hchamber_chis->Fill(fchamber->GetChisquare());
        }
    }

    // KS Test between high and low bin
    Double_t KS_result = h_low_chamberlevel->KolmogorovTest(h_high_chamberlevel);
    cout << "\n\n" <<  KS_result << "\n\n" << endl;


    // Write the Chi-Squared File to the chamber file
    hchamber_chis->Write();
    h_low_chamberlevel->Write();
    h_high_chamberlevel->Write();
    h3_low_chamberlevel->Write();
    h3_high_chamberlevel->Write();
    

    std::cout << "\n\n\nChamber Fits Performed!" << std::endl;


    /*--------------------------------------------------------------------------------------------------------------*\
    |                                                                                                                |
    |                                 At this point we have performed a fit to every chamber                         |
    |                                                                                                                |
    \*--------------------------------------------------------------------------------------------------------------*/
    
    // Extract each parameter as a vector from the fit tuple and convert to an array
    std::vector<float> v_param0 = std::get<0>(chamber_fits);
    std::vector<float> v_param1 = std::get<1>(chamber_fits);
    std::vector<float> v_param2 = std::get<2>(chamber_fits);
    std::vector<float> v_param3 = std::get<3>(chamber_fits);
    std::vector<float> v_param4 = std::get<4>(chamber_fits);

    // Create the TH1F histograms
    TH1F* h_param0 = new TH1F("fp0", "Fit parameter 0 Across all chambers; p0; chambers", 50, 50, 100);
    TH1F* h_param1 = new TH1F("fp1", "Fit parameter 1 Across all chambers; p1; chambers", 40, 20, 60);
    TH1F* h_param2 = new TH1F("fp2", "Fit parameter 2 Across all chambers; p2; chambers", 60, -15, 0);
    TH1F* h_param3 = new TH1F("fp3", "Fit parameter 3 Across all chambers; p3; chambers", 50, 0, 1);
    TH1F* h_param4 = new TH1F("fp4", "Fit parameter 4 Across all chambers; p4; chambers", 30, -0.03, 0);
    
    // Get the maximum number of hits used ina  claibrastion for binning purposes
    auto maxElementIter = std::max_element(hits_per_chamber.begin(), hits_per_chamber.end());
    int max_hits = *maxElementIter;

    int max_bin = max_hits - (max_hits % 1000) + 1000;
    int bins = max_bin / 1000;
    
    // Create the Th2F's for hits vs calibration parameters
    TH2F* hits_p0 = new TH2F("fp0_hits", "Hit's used for Calibration vs Calibrated p0; p0; Hits", 10, 50, 100, bins, 0, max_bin);
    TH2F* hits_p1 = new TH2F("fp1_hits", "Hit's used for Calibration vs Calibrated p1; p1; Hits", 10, 20, 60, bins, 0, max_bin);
    TH2F* hits_p2 = new TH2F("fp2_hits", "Hit's used for Calibration vs Calibrated p2; p2; Hits", 15, -15, 0, bins, 0, max_bin);
    TH2F* hits_p3 = new TH2F("fp3_hits", "Hit's used for Calibration vs Calibrated p3; p3; Hits", 10, 0, 1, bins, 0, max_bin);
    TH2F* hits_p4 = new TH2F("fp4_hits", "Hit's used for Calibration vs Calibrated p4; p4; Hits", 10, -0.03, 0, bins, 0, max_bin);

    // Set the bin contents from the std::vector<float>'s
    for (long unsigned int i=0; i<v_param0.size(); i++) {
        h_param0->Fill(v_param0.at(i));
        h_param1->Fill(v_param1.at(i));
        h_param2->Fill(v_param2.at(i));
        h_param3->Fill(v_param3.at(i));
        h_param4->Fill(v_param4.at(i));

        // Fill the TH2Fs for hits vs parameters
        hits_p0->Fill(v_param0.at(i), hits_per_chamber.at(i));
        hits_p1->Fill(v_param1.at(i), hits_per_chamber.at(i));
        hits_p2->Fill(v_param2.at(i), hits_per_chamber.at(i));
        hits_p3->Fill(v_param3.at(i), hits_per_chamber.at(i));
        hits_p4->Fill(v_param4.at(i), hits_per_chamber.at(i));
    }

    // Wite these histograms to the outfile
    h_param0->Write();
    h_param1->Write();
    h_param2->Write();
    h_param3->Write();
    h_param4->Write();
    hits_p0->Write();
    hits_p1->Write();
    hits_p2->Write();
    hits_p3->Write();
    hits_p4->Write();

    //Close the chamber level file
    chamber_level_file->Close();

    // Create the new TGraph and TH2 for which the calibrated values will be input using the polynomials
    TH2F* hcalpol4 = new TH2F("hcalpol4", "ADC vs. Drift Radius Calibrated (MuonTesterRun456714) | \
    (pol4); Drift Radius (mm); Calibrated ADC counts", numbins, 0, 15, 150, 0, 5);

    // A histogram just for the momnentum and pT of every muon eregardless of zmass requirements
    TH2F *h_p_pt = new TH2F("h_p_pt", "Momentum and Transverse Momentum of Muons; Momentum (GeV/c); pT (GeV/c)", 100,0,200, 100,0,200);
    // A histogram of the number of muons per event
    TH1F *h_muon_count = new TH1F("h_muon_count", "Muons per event; Muons; Events", 10,0,10);

    // Characterize the high and low momentum bins
    TH1F *h_low_eta = new TH1F("h_low_eta", "Eta for Muons from Zmumu in Low Momentum Bin; Eta; Tracks", 40,-1,1);
    TH1F *h_low_phi = new TH1F("h_low_phi", "Phi for Muons from Zmumu in Low Momentum Bin; Phi; Tracks", 20, -3.1418, 3.1418);
    TH1F *h_high_eta = new TH1F("h_high_eta", "Eta for Muons from Zmumu in High Momentum Bin; Eta; Tracks", 40,-1,1);
    TH1F *h_high_phi = new TH1F("h_high_phi", "Phi for Muons from Zmumu in High Momentum Bin; Phi; TRacks", 20, -3.1418, 3.1418);

    // Keep track of the track(let) dE/dx estimators
    std::vector<float> tracklet_truncated_means_20;
    std::vector<float> tracklet_truncated_means_40;
    std::vector<float> track_truncated_means_20;
    std::vector<float> track_truncated_means_40;

    // A subset of tracks which come from muons belonging to a muon pair with invariant mass
    // compatible with a Z
    std::vector<float> zmuon_tracklet_truncated_means_20;
    std::vector<float> zmuon_tracklet_truncated_means_40;
    std::vector<float> zmuon_track_truncated_means_20;
    std::vector<float> zmuon_track_truncated_means_40;
    std::vector<float> zmuon_tracklet_pts;
    std::vector<float> zmuon_track_pts;
    std::vector<float> zmuon_tracklet_totalp;
    std::vector<float> zmuon_track_totalp;
    std::vector<int> zmuon_tracklet_hits;
    std::vector<int> zmuon_track_hits;
    std::vector<float> zmuon_eta;
    std::vector<float> zmuon_phi;

    // Track muons from compatible invariant masses within near a high and low momentums chosen
    std::vector<float> low_momentum_edge_zmuon_tracklet_truncated_means_40;
    std::vector<float> high_momentum_edge_zmuon_tracklet_truncated_means_40;
    std::vector<float> low_momentum_edge_zmuon_track_truncated_means_40;
    std::vector<float> high_momentum_edge_zmuon_track_truncated_means_40;

    // Keep track of track(let) pTs
    std::vector<float> tracklet_pts;
    std::vector<float> track_pts;

    // Keep track of how many hits each instance of the estimator uses
    std::vector<int> tracklet_hits;
    std::vector<int> track_hits;

    // Keep track of track phi and eta
    std::vector<float> track_phis;
    std::vector<float> track_etas;

    /*
    Now loop through the TChain again (which is hot in cache so it should be faster)
    and apply the new fit to the adc counts based on their drift radii
    
    Simultaneously, select all the tracklets and apply appropriate processing
    */

    // Now do calibrations and record track(lets)
    std::cout << "\n\n\nCalibrating\n\n\n" << endl;
    for (Int_t i=0; i < nEntries; i++) {

        // A progress indicator
        if (i % 10000 == 0) {
            std::cout << i << "/" << nEntries << " processed" << endl;
        }

        // Get the event data
        chain.GetEntry(i);
        std::vector<int> current_event_Z_status(pts->size(), 0);

        // Make sure the current muon has some ADC counts recorded
        if (ADC_counts->size() == 0) {
            continue;
        }
        h_muon_count->Fill(pts->size());

        // Loop through all the muons in the event, keeping track of which pair if any is
        // closest to the invariant mass of the Z peak
        std::pair<int, int> zpeak_muons (-1,-1);
        int current_closest = 0;
        for (unsigned long n=0; n<pts->size() - 1; n++) {
            h_p_pt->Fill(ptToP(pts->at(n), etas->at(n)), pts->at(n));
            for (unsigned long j=n+1; j<pts->size(); j++) {
                float mumu_mass = invariantMass(pts->at(n), pts->at(j), etas->at(n),
                etas->at(j), phis->at(n), phis->at(j));
                // Skip pairs not within 20GeV of the zmass
                if (mumu_mass < zmass - 20 || mumu_mass > zmass + 20) {
                    continue;
                }
                else if (abs(zmass - mumu_mass) < abs(zmass - current_closest)) {
                    current_closest = mumu_mass;
                    zpeak_muons.first = n;
                    zpeak_muons.second = j;
                }
            }
            h_p_pt->Fill(ptToP(pts->at(pts->size()-1), etas->at(etas->size()-1)), pts->at(pts->size()-1));
        }

        // Mark the current muon link
        int muon_link = 0;
        int current_author = authors->at(muon_link);
        float muon_pt = pts->at(muon_link);
        float muon_eta = etas->at(muon_link);
        float muon_phi = phis->at(muon_link);

        // Keep track of what stations the current track passes through
        std::set<int> stations;
        std::vector<float> inner_large_hits;
        std::vector<float> inner_small_hits;
        std::vector<float> middle_large_hits;
        std::vector<float> middle_small_hits;
        std::vector<float> outer_large_hits;
        std::vector<float> outer_small_hits;
        std::map<int, std::vector<float>> station_hit_map;

        // Keep track of the current track's calibrated ADC counts
        std::vector<float> calibrated_tracklet;
        std::vector<float> calibrated_track;

        // Loop through all hits in the entry
        for (unsigned long n = 0; n < ADC_counts->size(); n++) {

            // Get the exact chamber and retrieve it's fit so we can compare chamber level calibrations
            std::string translated_station_index;
            std::string parsed_station_index = std::to_string(station_index->at(n));
            if (stationStrings.find(parsed_station_index) != stationStrings.end()) {
                translated_station_index = stationStrings[parsed_station_index];
            } else {
                continue;
            }



            // Place the hit into the correct momentum range
            int momentum_range = -1;
            if (ptToP(muon_pt, muon_eta) >= momentum_bins.at(0).first && ptToP(muon_pt, muon_eta) < momentum_bins.at(0).second) {
                momentum_range = 0;
            } else if (ptToP(muon_pt, muon_eta) >= momentum_bins.at(1).first && ptToP(muon_pt, muon_eta) < momentum_bins.at(1).second) {
                momentum_range = 1;
            }

            // Construct the chamber label
            std::string chamber_label = chamber_info[std::make_tuple(translated_station_index, station_eta->at(n), station_phi->at(n))];
            
            // Skip chambers / momneutm bins for which no polynomial was formed
            // if (chamber_polynomials.find(chamber_label) == chamber_polynomials.end()) { continue; } // SKip chambers without a polynomial
            // if (chamber_polynomials[chamber_label].find(momentum_range) == chamber_polynomials[chamber_label].end()) { continue;} // SKip bins without a polynomial

            /* Don't check for sMDTs in simulation data
            std::string sMDT1 = "BIS / Eta:7";
            std::string sMDT2 = "BIS / Eta:-7";
            if (chamber_label.rfind(sMDT1, 0) == 0 || chamber_label.rfind(sMDT2, 0) == 0) {
                cout << "drift radius is " << drift_radius->at(n) << endl;
                continue;
            }
            */
            

            // Initialize the calibrated ADC at the chamber level binned in momentum
            // float pol4_calibrated_ADC = ADC_counts->at(n) / polynomial(4, drift_radius->at(n), chamber_polynomials[chamber_label][momentum_range]);

            // Initialize the calibrated ADC at the global level unbinned in momentum
            //float pol4_calibrated_ADC = ADC_counts->at(n) / polynomial(4, drift_radius->at(n), globalPol4Params);
            
            float pol4_calibrated_ADC = (momentum_range & 1) ? ADC_counts->at(n) / polynomial(4, drift_radius->at(n), lowGlobalPol4Params) :
             ADC_counts->at(n) / polynomial(4, drift_radius->at(n), highGlobalPol4Params);

            // Check if we have a new muon
            if (trk_link->at(n) != muon_link) {

                // Update the muon link and author
                muon_link = trk_link->at(n);
                current_author = authors->at(muon_link);
                muon_pt = pts->at(muon_link);
                muon_eta = etas->at(muon_link);
                muon_phi = phis->at(muon_link);

                // Cut on empty/noise hits, outlier hits, and author drift radius)
                if (ADC_counts->at(n) < 50 || hit_type->at(n) > 60 || checkAuthor(current_author) == 0 \
                || muon_pt < 5 || abs(muon_eta) > 1 || !isInValidRadius(drift_radius->at(n))) {
                    continue;
                } else { // Fill the histogram and scatter plot if the above conditions are satisfied and calibrate the hit
                    hcalpol4->Fill(drift_radius->at(n), pol4_calibrated_ADC);
                }

                // Build the tracklet then compute the estimators. Also record tracklet pT
                calibrated_tracklet = buildTracklets(station_hit_map, stations);
                if (calibrated_tracklet.size() > 0) {
                    tracklet_truncated_means_20.emplace_back(TM(calibrated_tracklet, 0.2));
                    tracklet_truncated_means_40.emplace_back(TM(calibrated_tracklet, 0.4));
                    tracklet_hits.emplace_back(calibrated_tracklet.size());
                    tracklet_pts.emplace_back(muon_pt);
                    // CHeck if the muon is part of a muon pair compatible with an invariant mass at the z peak
                    if (zpeak_muons.first == muon_link || zpeak_muons.second == muon_link) {
                        zmuon_tracklet_truncated_means_20.emplace_back(TM(calibrated_tracklet, 0.2));
                        zmuon_tracklet_truncated_means_40.emplace_back(TM(calibrated_tracklet, 0.4));
                        zmuon_tracklet_hits.emplace_back(calibrated_tracklet.size());
                        zmuon_tracklet_pts.emplace_back(muon_pt);
                        zmuon_tracklet_totalp.emplace_back(ptToP(muon_pt, muon_eta));
                        zmuon_eta.emplace_back(muon_eta);
                        zmuon_phi.emplace_back(muon_phi);
                        // Check if this specific muon is within some range of the extreme momenta
                        if (abs(ptToP(muon_pt, muon_eta) - low_momentum_edge) < low_interval) {
                            low_momentum_edge_zmuon_tracklet_truncated_means_40.emplace_back(TM(calibrated_tracklet, 0.4));
                            h_low_eta->Fill(muon_eta);
                            h_low_phi->Fill(muon_phi);
                        } else if (abs(ptToP(muon_pt, muon_eta) - high_momentum_edge) < high_interval) {
                            high_momentum_edge_zmuon_tracklet_truncated_means_40.emplace_back(TM(calibrated_tracklet, 0.4));
                            h_high_eta->Fill(muon_eta);
                            h_high_phi->Fill(muon_phi);
                        }
                        
                    }
                }

                // If we have a full track, also build the track and its estimators. Also record track pT
                if (isFullTrack(station_hit_map)) {
                    calibrated_track = buildTracks(station_hit_map);
                    track_truncated_means_20.emplace_back(TM(calibrated_track, 0.2));
                    track_truncated_means_40.emplace_back(TM(calibrated_track, 0.4));
                    track_hits.emplace_back(calibrated_track.size());
                    track_pts.emplace_back(muon_pt);
                    track_phis.emplace_back(muon_phi);
                    track_etas.emplace_back(muon_eta);
                    if (zpeak_muons.first == muon_link || zpeak_muons.second == muon_link) {
                        zmuon_track_truncated_means_20.emplace_back(TM(calibrated_track, 0.2));
                        zmuon_track_truncated_means_40.emplace_back(TM(calibrated_track, 0.4));
                        zmuon_track_hits.emplace_back(calibrated_track.size());
                        zmuon_track_pts.emplace_back(muon_pt);
                        zmuon_track_totalp.emplace_back(ptToP(muon_pt, muon_eta));
                        zmuon_eta.emplace_back(muon_eta);
                        zmuon_phi.emplace_back(muon_phi);
                        if (abs(muon_pt - low_momentum_edge) < low_interval) {
                            low_momentum_edge_zmuon_track_truncated_means_40.emplace_back(TM(calibrated_track, 0.4));
                            h_low_eta->Fill(muon_eta);
                            h_low_phi->Fill(muon_phi);
                        } else if (abs(muon_pt - high_momentum_edge) < high_interval) {
                            high_momentum_edge_zmuon_track_truncated_means_40.emplace_back(TM(calibrated_track, 0.4));
                            h_high_eta->Fill(muon_eta);
                            h_high_phi->Fill(muon_phi);
                        }
                    }
                }

                // After all operations, start the new track(let) off with its first hit and reset all the track variables
                stations.clear();
                station_hit_map.clear();
                calibrated_tracklet.clear();
                calibrated_track.clear();
                stations.insert(static_cast<int>(station_index->at(n)));
                station_hit_map[static_cast<int>(station_index->at(n))].push_back(pol4_calibrated_ADC);
            
            // Now check the same cuts for hits when the muon hasn't changed
            } else if (checkAuthor(current_author) == 0 || ADC_counts->at(n) < 50 || hit_type->at(n) > 60 || \
            !isInValidRadius(drift_radius->at(n)) || muon_pt < 5 || abs(muon_eta) > 1) {
                continue;
            } else { // Fill the histogram and station hit map if the above conditions are satisfied
                hcalpol4->Fill(drift_radius->at(n), pol4_calibrated_ADC);
                stations.insert(static_cast<int>(station_index->at(n)));
                station_hit_map[static_cast<int>(station_index->at(n))].push_back(pol4_calibrated_ADC);
            }
        }
        // At the end of the event we want to add the last muon's tracklet and build the estimators
        // Check if the track is contained in the barrel
        if (isInBarrel(stations) == true) {
            calibrated_tracklet = buildTracklets(station_hit_map, stations);
            if (calibrated_tracklet.size() > 0) {
                tracklet_truncated_means_20.emplace_back(TM(calibrated_tracklet, 0.2));
                tracklet_truncated_means_40.emplace_back(TM(calibrated_tracklet, 0.4));
                tracklet_hits.emplace_back(calibrated_tracklet.size());
                tracklet_pts.emplace_back(muon_pt);
                if (zpeak_muons.first == muon_link || zpeak_muons.second == muon_link) {
                    zmuon_tracklet_truncated_means_20.emplace_back(TM(calibrated_tracklet, 0.2));
                    zmuon_tracklet_truncated_means_40.emplace_back(TM(calibrated_tracklet, 0.4));
                    zmuon_tracklet_hits.emplace_back(calibrated_tracklet.size());
                    zmuon_tracklet_pts.emplace_back(muon_pt);
                    zmuon_tracklet_totalp.emplace_back(ptToP(muon_pt, muon_eta));
                    zmuon_eta.emplace_back(muon_eta);
                    zmuon_phi.emplace_back(muon_phi);
                    if (abs(muon_pt - low_momentum_edge) < low_interval) {
                        low_momentum_edge_zmuon_tracklet_truncated_means_40.emplace_back(TM(calibrated_tracklet, 0.4));
                        h_low_eta->Fill(muon_eta);
                        h_low_phi->Fill(muon_phi);
                    } else if (abs(muon_pt - high_momentum_edge) < high_interval) {
                        high_momentum_edge_zmuon_tracklet_truncated_means_40.emplace_back(TM(calibrated_tracklet, 0.4));
                        h_high_eta->Fill(muon_eta);
                        h_high_phi->Fill(muon_phi);
                    }
                }
            }

            // If we have a full track, also build the track and its estimators. Also record track pT
            if (isFullTrack(station_hit_map)) {
                calibrated_track = buildTracks(station_hit_map);
                track_truncated_means_20.emplace_back(TM(calibrated_track, 0.2));
                track_truncated_means_40.emplace_back(TM(calibrated_track, 0.4));
                track_hits.emplace_back(calibrated_track.size());
                track_pts.emplace_back(muon_pt);
                track_phis.emplace_back(muon_phi);
                track_etas.emplace_back(muon_eta);
                if (zpeak_muons.first == muon_link || zpeak_muons.second == muon_link) {
                    zmuon_track_truncated_means_20.emplace_back(TM(calibrated_track, 0.2));
                    zmuon_track_truncated_means_40.emplace_back(TM(calibrated_track, 0.4));
                    zmuon_track_hits.emplace_back(calibrated_track.size());
                    zmuon_track_pts.emplace_back(muon_pt);
                    zmuon_track_totalp.emplace_back(ptToP(muon_pt, muon_eta));
                    zmuon_eta.emplace_back(muon_eta);
                    zmuon_phi.emplace_back(muon_phi);
                    if (abs(muon_pt - low_momentum_edge) < low_interval) {
                        low_momentum_edge_zmuon_track_truncated_means_40.emplace_back(TM(calibrated_track, 0.4));
                        h_low_eta->Fill(muon_eta);
                        h_low_phi->Fill(muon_phi);
                    } else if (abs(muon_pt - high_momentum_edge) < high_interval) {
                        high_momentum_edge_zmuon_track_truncated_means_40.emplace_back(TM(calibrated_track, 0.4));
                        h_high_eta->Fill(muon_eta);
                        h_high_phi->Fill(muon_phi);
                    }
                }
            }
        }
    }
    
    
    //------------------------------------------------------------------------------------------------------//
    // At this point we have created all of the tracks and all of the chambers are populated with their hits//
    //------------------------------------------------------------------------------------------------------//


    // We'll create a csv file of track de/dx's along with track pt, phi, eta, hits on track for NN simulation.
    std::ofstream csvFile("MuonTraining.csv");
    
    // The data to be written to the csv file
    csvFile << "dedx,nhits,pt,eta,phi\n";
    for (size_t i = 0; i < track_hits.size(); ++i) {
        csvFile << track_truncated_means_40[i] << "," << track_hits[i] << "," << track_pts[i] \
        << "," << track_etas[i] << "," << track_phis[i] << "\n";
    }
    csvFile.close();

    // Instead we need to get the landau fit to hcalpol4
    std::vector<std::tuple<float, float>> calibratedGlobalFitParams = fit2DHistogramTM(hcalpol4);
    Double_t calibratedGlobalMPVs[good_bins];
    Double_t calibratedGlobalMPVErrs[good_bins];

    // Loop through the MPVs we're concerned with and add to bin_centers and their corresponding mpv
    // We want to exclude the first and last mm of drift radius bins
    for (unsigned long j = excluded_bins; j < v_bin_centers.size() - excluded_bins; j++) {
        calibratedGlobalMPVs[j - excluded_bins] = std::get<0>(calibratedGlobalFitParams.at(j));
        calibratedGlobalMPVErrs[j - excluded_bins] = std::get<1>(calibratedGlobalFitParams.at(j));
    }

    // Construct a second TGraph to which we will fit a polynomial from the relevant Landau's
    TGraphErrors *g2calpol4 = new TGraphErrors(good_bins, bin_centers, calibratedGlobalMPVs, bin_centers_errs, calibratedGlobalMPVErrs);

    // Set the graph limits and axes/title
    g2calpol4->GetXaxis()->SetLimits(0,15);
    g2calpol4->SetMinimum(0);
    g2calpol4->SetMaximum(2);
    g2calpol4->SetTitle("ADC MPV vs. Drift Radius Calibrated (MuonTesterRun456714) | (pol4); Drift Radius (mm); ADC Count MPV");
    g2calpol4->SetMarkerStyle(7);

    // Construct the histograms for tracklet de/dx estimators
    TH1F *h20_tracklets = new TH1F("h20_tracklets", "Track (1/2 stations) dE/dx Estimator; MPV Estimator (20% Truncated Mean); Tracks", 100, 0, 2);
    TH1F *h40_tracklets = new TH1F("h40_tracklets", "Track (1/2 stations) dE/dx Estimator; MPV Estimator (40% Truncated Mean); Tracks", 100, 0, 2);
    TH1F *h20_tracks = new TH1F("h20_tracks", "Track dE/dx Estimator; MPV Estimator (20% Truncated Mean); Tracks", 100, 0, 2);
    TH1F *h40_tracks = new TH1F("h40_tracks", "Track dE/dx Estimator; MPV Estimator (40% Truncated Mean); Tracks", 100, 0, 2);
    for (unsigned long i = 0; i < tracklet_truncated_means_20.size(); i++) {
        h20_tracklets->Fill(tracklet_truncated_means_20.at(i));
        h40_tracklets->Fill(tracklet_truncated_means_40.at(i));
    }
    for (unsigned long i = 0; i < track_truncated_means_20.size(); i++) {
        h20_tracks->Fill(track_truncated_means_20.at(i));
        h40_tracks->Fill(track_truncated_means_40.at(i));
    }

    // Construct the histograms for tracklet de/dx estimators
    TH1F *h20_zmuon_tracklets = new TH1F("h20_zmuon_tracklets", "Track (1/2 stations) dE/dx Estimator; MPV Estimator (20% Truncated Mean); Tracks", 100, 0, 2);
    TH1F *h40_zmuon_tracklets = new TH1F("h40_zmuon_tracklets", "Track (1/2 stations) dE/dx Estimator; MPV Estimator (40% Truncated Mean); Tracks", 100, 0, 2);
    TH1F *h20_zmuon_tracks = new TH1F("h20_zmuon_tracks", "Track dE/dx Estimator; MPV Estimator (20% Truncated Mean); Tracks", 100, 0, 2);
    TH1F *h40_zmuon_tracks = new TH1F("h40_zmuon_tracks", "Track dE/dx Estimator; MPV Estimator (40% Truncated Mean); Tracks", 100, 0, 2);
    TH1F *h40_low_momentum_edge_zmuon_tracks = new TH1F("h40_low_momentum_edge_zmuon_tracks", ("Track dE/dx Estimator for Muons ~"+to_string(static_cast<int>(low_momentum_edge))+"GeV/c; MPV Estimator (40% Truncated Mean); Tracks").c_str(), 100, 0, 2);
    TH1F *h40_high_momentum_edge_zmuon_tracks = new TH1F("h40_high_momentum_edge_zmuon_tracks", ("Track dE/dx Estimator for Muons ~"+to_string(static_cast<int>(high_momentum_edge))+"GeV/c; MPV Estimator (40% Truncated Mean); Tracks").c_str(), 100, 0, 2);
    TH1F *h_zmuon_track_pts = new TH1F("h_zmuon_track_pts", "Transverse Momentum for Muons from Zmumu; pT (GeV/c); Tracks", 100,0,200);
    TH1F *h_zmuon_track_totalp = new TH1F("h_zmuon_track_totalp", "Total Momentum for Muons from Zmumu; p (GeV/c); Tracks", 200,0,200);
    TH1F *h_zmuon_eta = new TH1F("h_zmuon_eta", "Eta for Muons from Zmumu; Eta; Tracks", 40,-1,1);
    TH1F *h_zmuon_phi = new TH1F("h_zmuon_phi", "Phi for Muons from Zmumu; Phi; TRacks", 20, -3.1418, 3.1418);

    for (unsigned long i = 0; i < zmuon_tracklet_truncated_means_40.size(); i++) {
        h20_zmuon_tracklets->Fill(zmuon_tracklet_truncated_means_20.at(i));
        h40_zmuon_tracklets->Fill(zmuon_tracklet_truncated_means_40.at(i));
    }
    for (unsigned long i = 0; i < zmuon_track_truncated_means_40.size(); i++) {
        h20_zmuon_tracks->Fill(zmuon_track_truncated_means_20.at(i));
        h40_zmuon_tracks->Fill(zmuon_track_truncated_means_40.at(i));
        h_zmuon_track_totalp->Fill(zmuon_tracklet_totalp.at(i));
        h_zmuon_track_pts->Fill(zmuon_track_pts.at(i));
    }
    for (unsigned long i=0; i < low_momentum_edge_zmuon_track_truncated_means_40.size(); i++) {
        h40_low_momentum_edge_zmuon_tracks->Fill(low_momentum_edge_zmuon_track_truncated_means_40.at(i));
    }
    for (unsigned long i=0; i < high_momentum_edge_zmuon_track_truncated_means_40.size(); i++) {
        h40_high_momentum_edge_zmuon_tracks->Fill(high_momentum_edge_zmuon_track_truncated_means_40.at(i));
    }
    for (unsigned long i=0; i<zmuon_eta.size(); i++) {
        h_zmuon_eta->Fill(zmuon_eta.at(i));
        h_zmuon_phi->Fill(zmuon_phi.at(i));
    }

    // Create the output file and write all TH1F, TH2F and TGraphs to it
    TFile *output_file = new TFile(outfile.c_str(), "RECREATE");

    // Write the histograms for the Z->MuMu muons 
    h40_low_momentum_edge_zmuon_tracks->Write();
    h40_high_momentum_edge_zmuon_tracks->Write();
    h_zmuon_track_pts->Write();
    h_zmuon_track_totalp->Write();
    h_zmuon_eta->Write();
    h_zmuon_phi->Write();
    h_p_pt->Write();
    h_muon_count->Write();
    h_low_eta->Write();
    h_low_phi->Write();
    h_high_eta->Write();
    h_high_phi->Write();

    cout << "b1" << endl;
    // Call the function to make the Bethe-Bloch Curve for 
    TGraph *BBCurve = BetheBlochEstimatorCurve(h_zmuon_track_totalp, zmuon_track_totalp, zmuon_track_truncated_means_40, 20, 100, 10);
    BBCurve->SetTitle("Bethe Bloch Curve for Z Muons; Momentum (GeV); Average dE/dX Estimator");
    BBCurve->Write();
    cout << "b2" << endl;
    
    // Fit a single gaussian to the 40% TM estimator for Z peak muons
    // Draw the high edge estimator
    TCanvas *c40_Zmumu_single_gaus = new TCanvas("Single gaussian estimator");
    TF1 *tm40_Zmumu_single_gaus_high_edge = new TF1("tm_40_Zmumu_single_gaus_high_edge", "gaus", h40_high_momentum_edge_zmuon_tracks->GetXaxis()->GetXmin()\
    , h40_high_momentum_edge_zmuon_tracks->GetXaxis()->GetXmax());
    tm40_Zmumu_single_gaus_high_edge->SetParameters(0.2, 1, 0.1);
    h40_high_momentum_edge_zmuon_tracks->SetLineColor(kRed);
    h40_high_momentum_edge_zmuon_tracks->Fit(tm40_Zmumu_single_gaus_high_edge, "R+");
    h40_high_momentum_edge_zmuon_tracks->Draw("SAME");
    tm40_Zmumu_single_gaus_high_edge->Draw("SAME");
    TPaveText *tm40_Zmumu_label_single_gaus_high_edge = gaus_fit_label(tm40_Zmumu_single_gaus_high_edge);
    tm40_Zmumu_label_single_gaus_high_edge->Draw("SAME");
    
    // Draw the low edge estimator
    TF1 *tm40_Zmumu_single_gaus_low_edge = new TF1("tm_40_Zmumu_single_gaus_low_edge", "gaus", h40_low_momentum_edge_zmuon_tracks->GetXaxis()->GetXmin()\
    , h40_low_momentum_edge_zmuon_tracks->GetXaxis()->GetXmax());
    tm40_Zmumu_single_gaus_low_edge->SetParameters(0.2, 1, 0.1);
    h40_low_momentum_edge_zmuon_tracks->SetLineColor(kBlue);
    h40_low_momentum_edge_zmuon_tracks->Fit(tm40_Zmumu_single_gaus_low_edge, "R+");
    h40_low_momentum_edge_zmuon_tracks->Draw("SAME");
    TPaveText *tm40_Zmumu_label_single_gaus_low_edge = gaus_fit_label(tm40_Zmumu_single_gaus_low_edge);
    tm40_Zmumu_single_gaus_low_edge->Draw("SAME");
    tm40_Zmumu_label_single_gaus_low_edge->Draw("SAME");
    c40_Zmumu_single_gaus->Write();

    // Fit a double gaussian to the 40% TM estimator for Z peak muons
    TCanvas *c40_Zmumu_double_gaus_low_edge = new TCanvas("low bin double gaussian estimator");
    TF1 *tm40_Zmumu_double_gaus_low_edge = new TF1("tm_40_Zmumu_double_gaus_low_edge", "[0]*exp(-0.5*(x-[1])*(x-[1])/([2]*[2])) +\
     [3]*exp(-0.5*(x-[1])*(x-[1])/([4]*[4]))", h40_low_momentum_edge_zmuon_tracks->GetXaxis()->GetXmin()\
    , h40_low_momentum_edge_zmuon_tracks->GetXaxis()->GetXmax());
    tm40_Zmumu_double_gaus_low_edge->SetParameters(0.2, 1, 0.1, 0.2, 1.2, 0.1);
    h40_low_momentum_edge_zmuon_tracks->Fit(tm40_Zmumu_double_gaus_low_edge, "R");
    h40_low_momentum_edge_zmuon_tracks->Draw();
    TPaveText *tm40_Zmumu_label_double_gaus_low_edge = double_gaus_fit_label(tm40_Zmumu_double_gaus_low_edge);
    tm40_Zmumu_label_double_gaus_low_edge->Draw();
    c40_Zmumu_double_gaus_low_edge->Write();

    
    // Fit a double gaussian to the 40% TM estimator for Z peak muons
    TCanvas *c40_Zmumu_double_gaus_high_edge = new TCanvas("high bin double gaussian estimator");
    TF1 *tm40_Zmumu_double_gaus_high_edge = new TF1("tm_40_Zmumu_double_gaus_high_edge", "[0]*exp(-0.5*(x-[1])*(x-[1])/([2]*[2])) +\
     [3]*exp(-0.5*(x-[1])*(x-[1])/([4]*[4]))", h40_high_momentum_edge_zmuon_tracks->GetXaxis()->GetXmin()\
    , h40_high_momentum_edge_zmuon_tracks->GetXaxis()->GetXmax());
    tm40_Zmumu_double_gaus_high_edge->SetParameters(0.2, 1, 0.1, 0.2, 1.2, 0.1);
    h40_high_momentum_edge_zmuon_tracks->Fit(tm40_Zmumu_double_gaus_high_edge, "R");
    h40_high_momentum_edge_zmuon_tracks->Draw();
    TPaveText *tm40_Zmumu_label_double_gaus_high_edge = double_gaus_fit_label(tm40_Zmumu_double_gaus_high_edge);
    tm40_Zmumu_label_double_gaus_high_edge->Draw();
    c40_Zmumu_double_gaus_high_edge->Write();
    
    // Fit a gaussian to the 20% truncated mean histograms for TRACKLETS
    TCanvas *c20_tracklets = new TCanvas("20% Tracklets");
    TF1 *tm20_tracklets = new TF1("tm20_tracklets", "gaus", 0.5, 1.5);
    h20_tracklets->Fit(tm20_tracklets, "RQ");
    h20_tracklets->Draw();
    // Construct a fit label giving the relative resolution
    TPaveText *tm20_tracklet_label = gaus_fit_label(tm20_tracklets);
    tm20_tracklet_label->Draw();
    c20_tracklets->Write();
    h20_tracklets->Write();


    // Repeat for the 20% truncated mean TRACKS
    TCanvas *c20_tracks = new TCanvas("20% Tracks");
    TF1 *tm20_tracks = new TF1("tm20_tracks", "gaus", 0.5, 1.5);
    h20_tracks->Fit(tm20_tracks, "RQ");
    h20_tracks->Draw();
    TPaveText *tm20_track_label = gaus_fit_label(tm20_tracks);
    tm20_track_label->Draw();
    c20_tracks->Write();
    h20_tracks->Write();


    // Repeat for the 40% truncated mean TRACKLETS
    TCanvas *c40_tracklets = new TCanvas("40%TM Tracklets");
    TF1 *tm40_tracklets = new TF1("tm40_tracklets", "gaus", 0.5, 1.5);
    h40_tracklets->Fit(tm40_tracklets, "RQ");
    h40_tracklets->Draw();
    TPaveText *tm40_tracklet_label = gaus_fit_label(tm40_tracklets);
    tm40_tracklet_label->Draw();
    c40_tracklets->Write();
    h40_tracklets->Write();

    // Repeat for the 40% truncated mean TRACKS
    TCanvas *c40_tracks = new TCanvas("40%TM Tracks");
    TF1 *tm40_tracks = new TF1("tm40_tracks", "gaus", 0.5, 1.5);
    h40_tracks->Fit(tm40_tracks, "RQ");
    h40_tracks->Draw();
    TPaveText *tm40_track_label = gaus_fit_label(tm40_tracks);
    tm40_track_label->Draw();
    c40_tracks->Write();
    h40_tracks->Write();

    // Write the 2D histograms of ADC counts vs drift Radius
    h->Write();
    TF1* fpol4 = new TF1("fpol4","pol4", 1,14);
    g2->Fit(fpol4, "RQ");
    TF1 *fittedpol4 = g2->GetFunction("fpol4");
    TPaveText *pt_pol4 = pol_fit_label(4, fpol4);
    h->Draw("colz");
    g2->Draw("P same");
    fittedpol4->Draw("same");
    pt_pol4->Draw("same");
    TCanvas *c2d = gPad->GetCanvas();
    c2d->Write();

    // Write the momentum binned histograms of ADC counts vs drift radius
    hglobal_lowbin->Write();
    hglobal_highbin->Write();

    TCanvas* c_global_lowbin = new TCanvas("global lowbin correction");
    c_global_lowbin->SetName("global_lowbin");
    TF1* lowbin_pol4 = new TF1("lowbin_pol4","pol4", 1,14);
    g2_low->Fit(lowbin_pol4, "RQ");
    TPaveText *pt_lowbin = pol_fit_label(4, lowbin_pol4);
    hglobal_lowbin->Draw();
    lowbin_pol4->Draw("Same");
    pt_lowbin->Draw("same");
    c_global_lowbin->Write();

    TCanvas* c_global_highbin = new TCanvas("global highbin correction");
    TF1* highbin_pol4 = new TF1("highbin_pol4","pol4", 1,14);
    g2_high->Fit(highbin_pol4, "RQ");
    TPaveText *pt_highbin = pol_fit_label(4, highbin_pol4);
    hglobal_highbin->Draw();
    g2_high->Draw("Same");
    highbin_pol4->Draw("Same");
    pt_highbin->Draw("same");
    c_global_highbin->Write();

    TCanvas* c_global_comparison = new TCanvas("global_comparison", "Global Comparison", 800, 600);

    g2_low->SetMarkerStyle(20);  // Full circle
    g2_low->SetMarkerColor(kRed);
    g2_low->SetMarkerSize(0.4);
    g2_low->SetTitle("Global Drift Radius Correction in Momentum Bins");
    g2_low->Draw("AP"); 

    g2_high->SetMarkerStyle(21);  // Full square
    g2_high->SetMarkerColor(kBlue);
    g2_high->SetMarkerSize(0.4);
    g2_high->Draw("P same"); 

    lowbin_pol4->SetLineColor(kRed);
    lowbin_pol4->SetLineWidth(2);
    lowbin_pol4->Draw("same");
    highbin_pol4->SetLineColor(kBlue);
    highbin_pol4->SetLineWidth(2);
    highbin_pol4->Draw("same");

    // Create a legend
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);  // (x1, y1, x2, y2) in NDC coordinates
    legend->SetBorderSize(0);  // No border
    legend->SetFillStyle(0);   // Transparent background

    // Add entries to the legend
    legend->AddEntry(g2_low, "20-26GeV", "p");
    legend->AddEntry(g2_high, "55-85GeV", "p");
    legend->AddEntry(lowbin_pol4, "Low Bin Fit", "l");
    legend->AddEntry(highbin_pol4, "High Bin Fit", "l");

    // Draw the legend
    legend->Draw();

    c_global_comparison->Write();

    

    // Draw the calibrated TH2 ADC count graph with fit to constant
    hcalpol4->Write();
    TF1* fcalpol4 = new TF1("fcalpol4","[0]", 1,14);
    g2calpol4->Fit(fcalpol4, "RQ");
    TF1 *fittedcalpol4 = g2calpol4->GetFunction("fcalpol4");
    TPaveText *pt_pol4_cal = pol_fit_label(0, fcalpol4);
    hcalpol4->Draw("colz");
    g2calpol4->Draw("P same");
    fittedcalpol4->Draw("same");
    pt_pol4_cal->Draw("same");
    TCanvas *c2dcal = gPad->GetCanvas();
    c2dcal->Write();

    // The ADC curve to be calibrated with error bars (pol4)
    TCanvas *cpol4 = new TCanvas("ADC Polynomial Correction");
    g2->Draw("AP");
    TPaveText* global_pol4_fit_label = pol_fit_label(4, fglobalpol4);
    global_pol4_fit_label->Draw();
    cpol4->Write();

    // The ADC curve after calibration with 4th degree polynomial
    TCanvas *cpol4cal = new TCanvas("Calibrated MPVs");
    g2calpol4->Draw("AP");
    pt_pol4_cal->Draw();
    cpol4cal->Write();

    // Create Histograms for the size of tracklets and tracks
    TH1D *h_tracklet_sizes = new TH1D("h_tracklet_size", "Hits on Track (1/2 station tracks); Hits; Tracks", 50, 0, 50);
    TH1D *h_track_sizes = new TH1D("h_track_size", "Hits on track; Hits; Tracks", 50, 0, 50);

    // Create 2D histograms for track(let) sizes and their momentums
    TH2F *h_tracklet_pts = new TH2F("h_tracklet_pts", "Hits on Track (1/2 stations) vs. pT; Hits on Track;  Tracklet pT (GeV)", 50, 0, 50, 50, 0, 50);
    TH2F *h_track_pts = new TH2F("h_track_pts", "Hits on Track vs. pT; Hits on track; Track pT (GeV)", 50, 0, 50, 50, 0, 50);

    // Fill the hits on track histograms
    for (unsigned int i = 0; i < tracklet_hits.size(); i++) {
        h_tracklet_sizes->Fill(tracklet_hits.at(i));
        h_tracklet_pts->Fill(tracklet_hits.at(i), tracklet_pts.at(i));
    }
    for (unsigned int i = 0; i < track_hits.size(); i++) {
        h_track_sizes->Fill(track_hits.at(i));
        h_track_pts->Fill(track_hits.at(i), track_pts.at(i));
    }

    // Write the track(let) sizes histograms
    h_tracklet_sizes->Write();
    h_track_sizes->Write();
    h_tracklet_pts->Write();
    h_track_pts->Write();
    
    // Close the TFile
    output_file->Close();
}




int main(int argc, char* argv[]) {
    // Check for correct number of command line arguments
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <input_directory> <outfile> <chamber outfile>" << std::endl;
        return 1;
    }
    // Parse the command line arguments
    std::string indir(argv[1]);
    std::string outfile(argv[2]);
    std::string chamberoutfile(argv[3]);
    // Execute the primary function
    adcVsRadius(indir, outfile, chamberoutfile);
    return 0;
}
