#include <iostream>
#include <cmath>
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

// A function to test if a muon event is a full track extending through all barrel layers
bool isFullTrack(std::map<int, std::vector<float>> hitMap) {
    if ((!hitMap[0].empty() || !hitMap[1].empty()) && (!hitMap[2].empty() || !hitMap[3].empty()) \
    && (!hitMap[4].empty() || !hitMap[4].empty())) {
        return true;
    }
    else {return false;}
}

// A function to build a track hit vector
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
    for (int i = 0; i < numBins; i++) {
        binCenters.push_back(static_cast<float>(start + (i + 0.5) * binSize));
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
    float truncated_mean = track_sum / truncated_size;
    return (truncated_mean);
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

// Return the Landau fit parameters and errors of each xBin in a TH2 object
std::vector<std::tuple<float, float, float, float, float, float>> fit2DHistogramLandau(TH2 *h2) {
    int nbinsx = h2->GetNbinsX();
    std::vector<std::tuple<float, float, float, float, float, float>> fitParameters;

    for (int i = 1; i <= nbinsx; i++) {
        TH1D* projectionY = h2->ProjectionY("", i, i);
        TF1* landauFit = new TF1("landauFit", "landau", projectionY->GetXaxis()->GetXmin(), projectionY->GetXaxis()->GetXmax());

        // Set initial parameter values for the Landau function
        landauFit->SetParameters(100, 100, 30);

        // Fit the projection to the Landau function
        projectionY->Fit(landauFit, "Q");

        // Retrieve the fit parameters
        float norm = landauFit->GetParameter(0);
        float normErr = landauFit->GetParError(0);
        float mpv = landauFit->GetParameter(1);
        float mpvErr = landauFit->GetParError(1);
        float width = landauFit->GetParameter(2);
        float widthErr = landauFit->GetParError(2);
        std::tuple<float, float, float, float, float, float> binFit = std::make_tuple(norm, normErr, mpv, mpvErr, width, widthErr);

        fitParameters.push_back(binFit);

        delete projectionY;
        delete landauFit;

    }
    return fitParameters;
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

// Define a general 4th degree polynomial
float pol4_fit(float drift_radius, float p0, float p1, float p2, float p3, float p4) {
    return (p0 + p1*drift_radius + (p2*pow(drift_radius,2)) + (p3*pow(drift_radius,3)) \
     + (p4*pow(drift_radius,4)));
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




/*
* Primary Operations: Perform calibrations to account for ADC count dependency on drift radius
* Plot these fits, their paramters and the data from which they were constructed for each chamber
Then apply these corrections. Show the global effect of the corrections
* And finally produce 20, 40% truncated means of the corrected ADC counts for muon tracklets
*/
void adcVsRadius(std::string indir, const::std::string& outfile, const:: std::string& chamber_level_outfile) {
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

    // Create the map from hexadecimal labels to decimal label strings for the stationIndex
    std::map<std::string, std::string> stationStrings = stationTranslator();

    // Create the map from chamber info to the chamber label
    std::map<std::tuple<std::string, int, int>, std::string> chamber_info = buildChambers();

    // Create a map which links chamber labels to an array containing their ADC count
    // and drift radii. Ex. map[BOS_1_3] = (<10.6, 4.2,...>, <124, 142,...>)
    std::map<std::string, std::tuple<std::vector<float>, std::vector<float>>> chamber_data;
    for (const auto& entry: chamber_info) {
        const std::string& key = entry.second;
        chamber_data[key] = std::tuple<std::vector<float>, std::vector<float>>();
    }

    // Initialize the 2D histogram and the TGraph scatter plot
    TH2F* h = new TH2F("h", "ADC vs. Drift Radius; Drift Radius (mm); ADC Count MPV", 30, 0, 15, 400, 0, 400);
    TGraph *g1 = new TGraph();
    g1->SetMarkerStyle(5);

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
        
        // Loop through all hits in the entry
        for (unsigned long n = 0; n < ADC_counts->size(); n++) {
            // Get the index station as a string and translated it to the text label
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

            // Check if we have a new muon
            if (trk_link->at(n) != muon_link) {
                // Update the muon link
                muon_link = trk_link->at(n);
                // Cut on empty hits
                if (!ADC_counts->at(n)) {
                    continue;
                }
                else {
                    // Fill the histogram and scatter plot if the above conditions are satisfied
                    h->Fill(drift_radius->at(n), ADC_counts->at(n));
                    g1->AddPoint(drift_radius->at(n), ADC_counts->at(n));
                    // Add the data to the end of the data vectors in the tuple for the appropriate chamber
                    std::get<0>(chamber_data[chamber_label]).emplace_back(drift_radius->at(n));
                    // Convert the ADC counts to floats so we can use the TVectorF object to construct tgraphs later
                    float temp_counts = (float)ADC_counts->at(n);
                    std::get<1>(chamber_data[chamber_label]).emplace_back(temp_counts);
                }
            // Now check the same conditions for when it is not a new muon
            } else if (!ADC_counts->at(n)) {
                continue;
            } else {
                h->Fill(drift_radius->at(n), ADC_counts->at(n));
                g1->AddPoint(drift_radius->at(n), ADC_counts->at(n));
                // Add the data to the end of the data vectors in the tuple for the appropriate chamber
                std::get<0>(chamber_data[chamber_label]).emplace_back(drift_radius->at(n));
                // Convert the ADC counts to floats so we can use the TVectorF object to construct tgraphs later
                float temp_counts = (float)ADC_counts->at(n);
                std::get<1>(chamber_data[chamber_label]).emplace_back(temp_counts);
            }
        }
    }

    // Fit a Landau distribution to each bin and return the result here
    std::vector<std::tuple<float, float, float, float, float, float>> globalFitParams = fit2DHistogramLandau(h);

    // We're not using medians as our y-values for global calibration
    /*
    // medians has the medians in all drift radius bins
    std::vector<float>medians = median2D(h);
    Double_t valid_medians[medians.size()-4];
    Double_t bin_centers[medians.size()-4];
    */

    // create the bin centers vector and TVectorF
    std::vector<float> v_bin_centers = buildBinCenters(1, 14, 0.5);
    TVectorF f_bin_centers(v_bin_centers.size(), v_bin_centers.data());

    Double_t globalMPVs[v_bin_centers.size()];
    Double_t globalMPVErrs[v_bin_centers.size()];
    Double_t bin_centers[v_bin_centers.size()];
    Double_t bin_centers_errs[v_bin_centers.size()];
    // Loop through the medians we're concerned with and add to bin_centers and their corresponding median
    // We want to exclude the first two bins and the last two bins
    for (unsigned long j = 2; j < v_bin_centers.size() + 2; j++) {
        std::cout << "The MPV at " << v_bin_centers.at(j-2) << " mm is " << std::get<2>(globalFitParams.at(j)) \
        << " and the error is " << std::get<3>(globalFitParams.at(j)) << std::endl;
        globalMPVs[j-2] = std::get<2>(globalFitParams.at(j));
        globalMPVErrs[j-2] = std::get<3>(globalFitParams.at(j));
        bin_centers[j-2] = v_bin_centers.at(j-2);
        bin_centers_errs[j-2] = 0;
    }
    
    
    // Construct a second TGraph to which we will fit a polynomial from the relevant MPVs
    TGraphErrors *g2 = new TGraphErrors(v_bin_centers.size(), bin_centers, globalMPVs, bin_centers_errs, globalMPVErrs);
    // Set limits on the x and y ranges for the TGraph and make the marker style
    TAxis *xaxis = g2->GetXaxis();
    xaxis->SetLimits(0,14);
    g2->SetMinimum(0);
    g2->SetMaximum(160);
    g2->SetTitle("ADC medians vs. Drift Radius; Drift Radius (mm); ADC Count MPV");
    g2->SetMarkerStyle(7);
    
    // Fit a 4th degree polynomial to the median data
    TF1* f1 = new TF1("f1","pol4", 1,14);
    g2->Fit(f1, "RQ");

    // retrieve the global fit parameters
    float p0 = f1->GetParameter(0);
    float p1 = f1->GetParameter(1);
    float p2 = f1->GetParameter(2);
    float p3 = f1->GetParameter(3);
    float p4 = f1->GetParameter(4);
    float global_fit_quality = f1->GetChisquare() / f1->GetNDF();
    const Double_t *global_parErrors = f1->GetParErrors();

    // Create the fit label to be displayed
    TPaveText* global_fit_label = new TPaveText(0.7,0.1,0.9,0.3,"NDC");
    global_fit_label->SetTextAlign(12);
    global_fit_label->AddText(Form("Fit Parameters"));
    global_fit_label->AddText(Form("p0: %f | Err: %f", p0, global_parErrors[0]));
    global_fit_label->AddText(Form("p1: %f | Err: %f", p1, global_parErrors[1]));
    global_fit_label->AddText(Form("p2: %f | Err: %f", p2, global_parErrors[2]));     
    global_fit_label->AddText(Form("p3: %f | Err: %f", p3, global_parErrors[3]));
    global_fit_label->AddText(Form("p4: %f | Err: %f", p4, global_parErrors[4]));
    global_fit_label->AddText(Form("ChiSquare / NDF %f", global_fit_quality));


    // Initialize vectors which will holds the fit parameters for each chamber
    std::tuple<std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>> chamber_fits;

    // Open the chamber level TFile
    TFile *chamber_level_file = new TFile(chamber_level_outfile.c_str(), "RECREATE");

    // Create a vector to track the number of hits per chamber used for calibration
    std::vector<int> hits_per_chamber;



    // Now Create a TGraph object for every chamber and populate it with the correesponding TVectors
    for (const auto& entry: chamber_data) {
        // Get the chamber label
        const std::string& key = entry.first;
        // retrieve the Vectors of data for the current chamber and cast them as TVectorF objects
        vector<float> drift_radius_data = std::get<0>(chamber_data[key]);
        vector<float> ADC_count_data = std::get<1>(chamber_data[key]);
        if (ADC_count_data.size() == 0) {
            continue;
        }

        // Record the number of hits to be used for this calibration
        hits_per_chamber.emplace_back(drift_radius_data.size());

        // Title the plot for this chamber and initialize the TH2F
        std::string graph_title_with_axes = key + ";Drift Radius(mm); ADC Counts";
        TH2F* hist = new TH2F(key.c_str(), graph_title_with_axes.c_str(), 30, 0, 15, 400, 0, 400);

        // Add the vectors to the histogram
        for (long unsigned hit_i = 0; hit_i < drift_radius_data.size(); ++hit_i) {
            hist->Fill(drift_radius_data.at(hit_i), ADC_count_data.at(hit_i));
        }
        hist->Write();

        // Fit a Landau distribution to each bin and return the result here
        std::vector<std::tuple<float, float, float, float, float, float>> chamberFitParams = fit2DHistogramLandau(hist);

        Double_t chamberMPVs[v_bin_centers.size()];
        Double_t chamberMPVErrs[v_bin_centers.size()];
        // Loop through the medians we're concerned with and add to bin_centers and their corresponding median
        // We want to exclude the first two bins and the last two bins
        for (unsigned long j = 2; j < v_bin_centers.size() + 2; j++) {
            chamberMPVs[j-2] = std::get<2>(chamberFitParams.at(j));
            chamberMPVErrs[j-2] = std::get<3>(chamberFitParams.at(j));
        }

        // Construct the TGraph from the TVectors
        TGraphErrors *g = new TGraphErrors(v_bin_centers.size(), bin_centers, chamberMPVs, bin_centers_errs, chamberMPVErrs);
        g->SetNameTitle(key.c_str(), graph_title_with_axes.c_str());
        g->SetMarkerStyle(7);
        g->SetMinimum(0);

        // Now perform the fit to these median plots for each chamber
        TF1* fchamber = new TF1("fchamber","pol4", 1,14);
        g->Fit(fchamber, "RQ");

        // Extract the fit parameters
        const Double_t *pErrors = fchamber->GetParErrors();
        float p0_temp = fchamber->GetParameter(0);
        float p1_temp = fchamber->GetParameter(1);
        float p2_temp = fchamber->GetParameter(2);
        float p3_temp = fchamber->GetParameter(3);
        float p4_temp = fchamber->GetParameter(4);
        float chamber_fit_quality = fchamber->GetChisquare() / fchamber->GetNDF();

        // Add the fit parameters to the paramater vectors for all chambers
        std::get<0>(chamber_fits).emplace_back(p0_temp);
        std::get<1>(chamber_fits).emplace_back(p1_temp);
        std::get<2>(chamber_fits).emplace_back(p2_temp);
        std::get<3>(chamber_fits).emplace_back(p3_temp);
        std::get<4>(chamber_fits).emplace_back(p4_temp);
        
        // Create a canvas so we can draw the TGraph with the preferred options and write
        // the canvas to the chamber level fileß
        TCanvas *cchamber = new TCanvas(key.c_str(), graph_title_with_axes.c_str(), 800, 600);
        g->Draw("AP");

        // Construct the fit label and draw it
        TPaveText* chamber_fit_label = new TPaveText(0.7,0.1,0.9,0.3,"NDC");
        chamber_fit_label->SetTextAlign(12);
        const char* formatted_label_0 = Form("Fit Parameters");
        const char* formatted_label_1 = Form("p0: %f | Err: %f", p0_temp, pErrors[0]);
        const char* formatted_label_2 = Form("p1: %f | Err: %f", p1_temp, pErrors[1]);
        const char* formatted_label_3 = Form("p2: %f | Err: %f", p2_temp, pErrors[2]);
        const char* formatted_label_4 = Form("p3: %f | Err: %f", p3_temp, pErrors[3]);
        const char* formatted_label_5 = Form("p4: %f | Err: %f", p4_temp, pErrors[4]);
        const char* formatted_label_6 = Form("ChiSquare / NDF %f", chamber_fit_quality);
        
        chamber_fit_label->AddText(formatted_label_0);
        chamber_fit_label->AddText(formatted_label_1);
        chamber_fit_label->AddText(formatted_label_2);
        chamber_fit_label->AddText(formatted_label_3);     
        chamber_fit_label->AddText(formatted_label_4);
        chamber_fit_label->AddText(formatted_label_5);
        chamber_fit_label->AddText(formatted_label_6);
        chamber_fit_label->Draw();
        cchamber->cd();

        // Write the TCanvas to the file
        cchamber->Write();
    }
    
    // Extract each parameter as a vector from the fit tuple and convert to an array
    std::vector<float> v_param0 = std::get<0>(chamber_fits);
    std::vector<float> v_param1 = std::get<1>(chamber_fits);
    std::vector<float> v_param2 = std::get<2>(chamber_fits);
    std::vector<float> v_param3 = std::get<3>(chamber_fits);
    std::vector<float> v_param4 = std::get<4>(chamber_fits);

    // Create the TH1F histograms
    TH1F* h_param0 = new TH1F("fp0", "Fit parameter 0 Across all chambers", 50, 20, 120);
    TH1F* h_param1 = new TH1F("fp1", "Fit parameter 1 Across all chambers", 50, 0, 100);
    TH1F* h_param2 = new TH1F("fp2", "Fit parameter 2 Across all chambers", 50, -30, 20);
    TH1F* h_param3 = new TH1F("fp3", "Fit parameter 3 Across all chambers", 50, -1, 3);
    TH1F* h_param4 = new TH1F("fp4", "Fit parameter 4 Across all chambers", 40, -0.1, 0.1);
    
    // Get the maximum number of hits used ina  claibrastion for binning purposes
    auto maxElementIter = std::max_element(hits_per_chamber.begin(), hits_per_chamber.end());
    int max_hits = *maxElementIter;

    int max_bin = max_hits - (max_hits % 1000) + 1000;
    int bins = max_bin / 1000;
    
    // Create the Th2F's for hits vs calibration parameters
    TH2F* hits_p0 = new TH2F("fp0_hits", "Hit's used for Calibration vs Calibrated p0; p0; Hits", 10, 20, 120, bins, 0, max_bin);
    TH2F* hits_p1 = new TH2F("fp1_hits", "Hit's used for Calibration vs Calibrated p1; p1; Hits", 10, 0, 100, bins, 0, max_bin);
    TH2F* hits_p2 = new TH2F("fp2_hits", "Hit's used for Calibration vs Calibrated p2; p2; Hits", 10, -30, 20, bins, 0, max_bin);
    TH2F* hits_p3 = new TH2F("fp3_hits", "Hit's used for Calibration vs Calibrated p3; p3; Hits", 10, -1, 3, bins, 0, max_bin);
    TH2F* hits_p4 = new TH2F("fp4_hits", "Hit's used for Calibration vs Calibrated p4; p4; Hits", 10, -0.1, 0.1, bins, 0, max_bin);

    // Set the bin contents from the std::vector<float>'s
    for (long unsigned int i=0; i<v_param0.size(); i++) {
        h_param0->Fill(v_param0.at(i));
        h_param1->Fill(v_param1.at(i));
        h_param2->Fill(v_param2.at(i));
        h_param3->Fill(v_param3.at(i));
        h_param4->Fill(v_param4.at(i));

        // Fill the TH2Fs for hits vs paramters
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

    // Create the new TGraph and TH2 for which the calibrated values will be input
    TH2F* hcal = new TH2F("hcal", "ADC vs. Drift Radius Calibrated; Drift Radius (mm); Calibrated ADC counts", 30, 0, 15, 150, 0, 5);
    TGraph* g1cal = new TGraph();
    g1cal->SetMarkerStyle(7);







    /*
    Now loop through the TChain again (which is hot in cache so it should be faster)
    and apply the new fit to the adc counts based on their drift radii
    
    Simultaneously, select all the tracklets and apply appropriate processing
    */

    // Keep track of the estimators when truncating means by 20 and 40 percent
    std::vector<float> tracklet_truncated_means_20;
    std::vector<float> tracklet_truncated_means_40;
    std::vector<float> track_truncated_means_20;
    std::vector<float> track_truncated_means_40;

    // Keep track of how many hits each instance of the estimator uses
    std::vector<int> tracklet_hits;
    std::vector<int> track_hits;

    // Now do calibrations and record tracklets
    std::cout << "\n\n\nCalibrating\n\n\n" << endl;
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
        int current_author = authors->at(muon_link);
        

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

            // Apply the fit function to the adc counts
            float calibrated_ADC = ADC_counts->at(n) / pol4_fit(drift_radius->at(n), p0, p1, p2, p3, p4);

            // Check if we have a new muon
            if (trk_link->at(n) != muon_link) {
                
                // Update the muon link and author
                muon_link = trk_link->at(n);
                current_author = authors->at(muon_link);

                // Cut on empty hits, outlier hits, and author (only CM or standalone)
                if (!(*ADC_counts)[n] || (*hit_type)[n] > 60 || checkAuthor(current_author) == 0) {
                    continue;
                }
                else {
                    // Fill the histogram and scatter plot if the above conditions are satisfied
                    hcal->Fill(drift_radius->at(n), calibrated_ADC);
                    g1cal->AddPoint(drift_radius->at(n), calibrated_ADC);
                }

                // Check if the track is contained in the barrel
                if (isInBarrel(stations) == true) {
                    // Build the tracklet then compute
                    calibrated_tracklet = buildTracklets(station_hit_map, stations);
                    if (calibrated_tracklet.size() > 0) {
                        tracklet_truncated_means_20.emplace_back(TM(calibrated_tracklet, 0.2));
                        tracklet_truncated_means_40.emplace_back(TM(calibrated_tracklet, 0.4));
                        tracklet_hits.emplace_back(calibrated_tracklet.size());
                    }
                    // If we have a full track, also build the track and its estimators
                    if (isFullTrack(station_hit_map)) {
                        calibrated_track = buildTracks(station_hit_map);
                        track_truncated_means_20.emplace_back(TM(calibrated_track, 0.2));
                        track_truncated_means_40.emplace_back(TM(calibrated_track, 0.4));
                        track_hits.emplace_back(calibrated_track.size());
                    }
                }

                // After all operations, start the new track(let) off with its first hit and reset all the track variables
                stations = {};
                station_hit_map.clear();
                calibrated_tracklet.clear();
                calibrated_track.clear();
                stations.insert(static_cast<int>(station_index->at(n)));
                station_hit_map[static_cast<int>(station_index->at(n))].push_back(calibrated_ADC);
                
            // Now check the same cuts for hits wehn the muon hasn't changed
            } else if (!(*ADC_counts)[n] || (*hit_type)[n] > 60 || checkAuthor(current_author) == 0) {
                continue;
            } else {
                // Fill the histogram and scatter plot if the above conditions are satisfied
                hcal->Fill(drift_radius->at(n), calibrated_ADC);
                g1cal->AddPoint(drift_radius->at(n), calibrated_ADC);

                // Make sure the stations set contains this hit's station
                stations.insert(static_cast<int>(station_index->at(n)));
                
                // Add the calibrated ADC count to the track vector and the station track specifically
                calibrated_tracklet.push_back(calibrated_ADC);
                station_hit_map[static_cast<int>(station_index->at(n))].push_back(calibrated_ADC);
            }
        }
        // At the end of the event we want to add the last muon's tracklet and build the estimators
        //// Check if the track is contained in the barrel
        if (isInBarrel(stations) == true) {
            calibrated_tracklet = buildTracklets(station_hit_map, stations);
            if (calibrated_tracklet.size() > 0) {
                tracklet_truncated_means_20.emplace_back(TM(calibrated_tracklet, 0.2));
                tracklet_truncated_means_40.emplace_back(TM(calibrated_tracklet, 0.4));
                tracklet_hits.emplace_back(calibrated_tracklet.size());
            }
            // If we have a full track, also build the track and its estimators
            if (isFullTrack(station_hit_map)) {
                calibrated_track = buildTracks(station_hit_map);
                track_truncated_means_20.emplace_back(TM(calibrated_track, 0.2));
                track_truncated_means_40.emplace_back(TM(calibrated_track, 0.4));
                track_hits.emplace_back(calibrated_track.size());
            }
        }
    }
    
    //------------------------------------------------------------------------------------------------------//
    // At this point we have created all of the tracks and all of the chambers are populated with their hits//
    //------------------------------------------------------------------------------------------------------//


    // medians has the medians in all drift radius bins
    std::vector<float>medians_cal = median2D(hcal);
    Double_t valid_medians_cal[medians_cal.size()-4];

    // Loop through the medians we're concerned with and build the calibrated median array
    for (unsigned long j=2; j<medians_cal.size()-2; j++) {
        valid_medians_cal[j-2] = medians_cal.at(j);
    }

    // Construct a second TGraph to which we will fit a polynomial from the relevant medians
    TGraph *g2cal = new TGraph(v_bin_centers.size(), bin_centers, valid_medians_cal);

    // Set limits on the x and y ranges for the TGraph
    TAxis *xaxiscal = g2cal->GetXaxis();
    xaxiscal->SetLimits(0,15);

    // Note these extrema are different because we are expecting a fit around the constant 1
    g2cal->SetMinimum(0);
    g2cal->SetMaximum(2);
    g2cal->SetTitle("ADC medians vs. Drift Radius Calibrated; Drift Radius (mm); ADC Count MPV(Arbitrary Units)");
    g2cal->SetMarkerStyle(7);

    // Fit a 4th degree polynomial to the median data
    TF1* fcal = new TF1("fcal","[0]", 1,14);
    g2cal->Fit(fcal, "RQ");
    float constant_fit = fcal->GetParameter(0);

    // Construct the histograms for tracklet de/dx estimators
    TH1F *h20_tracklets = new TH1F("h20_tracklets", "Tracklet dE/dx Estimator; 20% Truncated Mean; Tracklets", 100, 0, 2);
    TH1F *h40_tracklets = new TH1F("h40_tracklets", "Tracklet dE/dx Estimator; 40% Truncated Mean; Tracklets", 100, 0, 2);
    TH1F *h20_tracks = new TH1F("h20_tracks", "Track dE/dx Estimator; 20% Truncated Mean; Tracks", 100, 0, 2);
    TH1F *h40_tracks = new TH1F("h40_tracks", "Track dE/dx Estimator; 40% Truncated Mean; Tracks", 100, 0, 2);
    for (unsigned long i = 0; i < tracklet_truncated_means_20.size(); i++) {
        h20_tracklets->Fill(tracklet_truncated_means_20.at(i));
        h40_tracklets->Fill(tracklet_truncated_means_40.at(i));
    }
    for (unsigned long i = 0; i < track_truncated_means_20.size(); i++) {
        h20_tracks->Fill(track_truncated_means_20.at(i));
        h40_tracks->Fill(track_truncated_means_40.at(i));
    }


    
    // Create the output file and write all TH1F, TH2F and TGraphs to it
    TFile *output_file = new TFile(outfile.c_str(), "RECREATE");
    TCanvas *c20_tracklets = new TCanvas();

    // Fit a gaussian to the 20% truncated mean histograms for TRACKLETS
    TF1 *tm20_tracklets = new TF1("tm20_tracklets", "gaus", 0.5, 1.5);
    h20_tracklets->Fit(tm20_tracklets, "RQ");
    h20_tracklets->Draw();
    
    // Construct a fit label giving the relative resolution
    TPaveText *tm20_tracklet_label = new TPaveText(0.777,0.753,0.977,0.953,"NDC");
    tm20_tracklet_label->AddText(Form("Mean: %f", tm20_tracklets->GetParameter(1)));
    tm20_tracklet_label->AddText(Form("Std dev: %f", tm20_tracklets->GetParameter(2)));
    tm20_tracklet_label->AddText(Form("Relative Res: %f", tm20_tracklets->GetParameter(2) / tm20_tracklets->GetParameter(1)));
    tm20_tracklet_label->Draw();
    c20_tracklets->Write();


    // Repeat for the 20% truncated mean TRACKS
    TCanvas *c20_tracks = new TCanvas();
    TF1 *tm20_tracks = new TF1("tm20_tracks", "gaus", 0.5, 1.5);
    h20_tracks->Fit(tm20_tracks, "RQ");
    h20_tracks->Draw();
    TPaveText *tm20_track_label = new TPaveText(0.777,0.753,0.977,0.953,"NDC");
    tm20_track_label->AddText(Form("Mean: %f", tm20_tracks->GetParameter(1)));
    tm20_track_label->AddText(Form("Std dev: %f", tm20_tracks->GetParameter(2)));
    tm20_track_label->AddText(Form("Relative Res: %f", tm20_tracks->GetParameter(2) / tm20_tracks->GetParameter(1)));
    tm20_track_label->Draw();
    c20_tracks->Write();

    // Repeat for the 40% truncated mean TRACKLETS
    TCanvas *c40_tracklets = new TCanvas();
    TF1 *tm40_tracklets = new TF1("tm40_tracklets", "gaus", 0.5, 1.5);
    h40_tracklets->Fit(tm40_tracklets, "RQ");
    h40_tracklets->Draw();
    TPaveText *tm40_tracklet_label = new TPaveText(0.777,0.753,0.977,0.953,"NDC");
    tm40_tracklet_label->AddText(Form("Mean: %f", tm40_tracklets->GetParameter(1)));
    tm40_tracklet_label->AddText(Form("Std dev: %f", tm40_tracklets->GetParameter(2)));
    tm40_tracklet_label->AddText(Form("Relative Res: %f", tm40_tracklets->GetParameter(2) / tm40_tracklets->GetParameter(1)));
    tm40_tracklet_label->Draw();
    c40_tracklets->Write();

    // Repeat for the 40% truncated mean TRACKS
    TCanvas *c40_tracks = new TCanvas();
    TF1 *tm40_tracks = new TF1("tm40_tracks", "gaus", 0.5, 1.5);
    h40_tracks->Fit(tm40_tracks, "RQ");
    h40_tracks->Draw();
    TPaveText *tm40_track_label = new TPaveText(0.777,0.753,0.977,0.953,"NDC");
    tm40_track_label->AddText(Form("Mean: %f", tm40_tracks->GetParameter(1)));
    tm40_track_label->AddText(Form("Std dev: %f", tm40_tracks->GetParameter(2)));
    tm40_track_label->AddText(Form("Relative Res: %f", tm40_tracks->GetParameter(2) / tm40_tracks->GetParameter(1)));
    tm40_track_label->Draw();
    c40_tracks->Write();


    // Write the 2D histograms of ADC counts vs drift Radius
    h->Write();
    hcal->Write();

    // The ADC curve to be calibrated with error barsm
    TCanvas *c = new TCanvas();
    g2->Draw("AP");
    global_fit_label->Draw();
    c->Write();

    // The ADC curve after calibration
    TCanvas *c1 = new TCanvas();
    g2cal->Draw("AP");

    // Add the fit label
    TPaveText *pt_cal = new TPaveText(0.7,0.1,0.9,0.3,"NDC");
    char text_cal[100];
    sprintf(text_cal, "p0: %f", constant_fit);
    pt_cal->AddText(text_cal);
    pt_cal->Draw();
    c1->Write();

    // Create Histograms for the size of tracklets and tracks
    TH1D *h_tracklet_sizes = new TH1D("h_tracklet_size", "Tracklet Sizes; Hits", 50, 0, 50);
    TH1D *h_track_sizes = new TH1D("h_track_size", "Track Sizes; Hits", 50, 0, 50);
    for (unsigned int i = 0; i < tracklet_hits.size(); i++) {
        h_tracklet_sizes->Fill(tracklet_hits.at(i));
    }
    for (unsigned int i = 0; i < track_hits.size(); i++) {
        h_track_sizes->Fill(track_hits.at(i));
    }

    // Draw and save the number of hits in the track(let)s
    TCanvas *c_trackletInstances = new TCanvas();
    h_tracklet_sizes->Draw();
    c_trackletInstances->Write();

    TCanvas *c_trackInstances = new TCanvas();
    h_track_sizes->Draw();
    c_trackInstances->Write();
    
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
