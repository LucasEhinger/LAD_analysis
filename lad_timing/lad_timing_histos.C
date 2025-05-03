#include <TFile.h>
#include <TH1.h> // Include ROOT histogram header
#include <TTree.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

const static int NPLANES                      = 6;    // Number of planes
const static int NBARS                        = 11;   // Number of bars
const static int MAXDATA                      = 1000; // Maximum number of data points
const static std::string PLANE_NAMES[NPLANES] = {"000", "001", "100", "101", "200", "REFBAR"}; // Names of the planes
const static string PLANE_NAMES[NPLANES]      = {"000", "001", "100", "101", "200", "REFBAR"}; // Names of the planes

struct histRange {
  int nxbins;       // Number of bins along x-axis
  double xmin;      // Minimum value along x-axis
  double xmax;      // Maximum value along x-axis
  int nybins;       // Number of bins along y-axis (optional for 2D histograms)
  double ymin;      // Minimum value along y-axis (optional for 2D histograms)
  double ymax;      // Maximum value along y-axis (optional for 2D histograms)
  std::string opt;  // Options for histogram creation (e.g., "COLZ", "E")
};

const static histRange hr_TDC_raw = {100, 0, 40000, 0, 0, 0, ""}; // Histogram range for TDC raw data
const static histRange hr_ADC_raw = {100, 0, 500, 0, 0, 0, ""}; // Histogram range for ADC raw data
const static histRange hr_TDC_time = {100, 0, 40000, 0, 0, 0, ""}; // Histogram range for TDC time
const static histRange hr_ADC_Difftime = {100, -40, 40, 0, 0, 0, ""}; // Histogram range for ADC time
const static histRange hr_Edep = {100, 0, 300, 0, 0, 0, ""}; // Histogram range for ADC energy deposition

void lad_timing_histos(const char *filename) {
  // Open the ROOT file
  TFile *file = TFile::Open(filename, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Cannot open file " << filename << std::endl;
    return;
  }

  // Get the tree T
  TTree *tree = (TTree *)file->Get("T");
  if (!tree) {
    std::cerr << "Error: Tree T not found in file " << filename << std::endl;
    file->Close();
    return;
  }

  // Define variables to hold branch data
  Double_t TopTdcCounter[NPLANES][MAXDATA];      // 2D array for planes and data points
  Double_t TopTdcTimeRaw[NPLANES][MAXDATA];         // 2D array for planes and data points
  Int_t NdataTopTdcTime;                         // Number of data points for each plane
  Double_t TopAdcCounter[NPLANES][MAXDATA];      // 2D array for planes and data points
  Double_t TopAdcPulseTimeRaw[NPLANES][MAXDATA]; // 2D array for planes and data points
  Int_t NdataTopAdcPulseTimeRaw;                 // Number of data points for each plane
  Double_t BtmTdcCounter[NPLANES][MAXDATA];      // 2D array for planes and data points
  Double_t BtmTdcTimeRaw[NPLANES][MAXDATA];         // 2D array for planes and data points
  Int_t NdataBtmTdcTime;                         // Number of data points for each plane
  Double_t BtmAdcCounter[NPLANES][MAXDATA];      // 2D array for planes and data points
  Double_t BtmAdcPulseTimeRaw[NPLANES][MAXDATA]; // 2D array for planes and data points
  Int_t NdataBtmAdcPulseTimeRaw;                 // Number of data points for each plane
  Double_t GoodHitTimeAvg[NPLANES][NBARS];       // 2D array for planes and data points
  Double_t GoodHitTimeDiff[NPLANES][NBARS];
  Double_t HodoHitEdep[NPLANES][NBARS];

  // Set branch addresses
  for (int i = 0; i < NPLANES; ++i) {
    std::string branchNameTopTdcCounter      = "H.ladhod." + PLANE_NAMES[i] + ".TopTdcCounter";
    std::string branchNameTopTdcTimeRaw         = "H.ladhod." + PLANE_NAMES[i] + ".TopTdcTimeRaw";
    std::string branchNameTopAdcCounter      = "H.ladhod." + PLANE_NAMES[i] + ".TopAdcCounter";
    std::string branchNameTopAdcPulseTimeRaw = "H.ladhod." + PLANE_NAMES[i] + ".TopAdcPulseTimeRaw";
    std::string branchNameBtmTdcCounter      = "H.ladhod." + PLANE_NAMES[i] + ".BtmTdcCounter";
    std::string branchNameBtmTdcTimeRaw         = "H.ladhod." + PLANE_NAMES[i] + ".BtmTdcTimeRaw";
    std::string branchNameBtmAdcCounter      = "H.ladhod." + PLANE_NAMES[i] + ".BtmAdcCounter";
    std::string branchNameBtmAdcPulseTimeRaw = "H.ladhod." + PLANE_NAMES[i] + ".BtmAdcPulseTimeRaw";
    std::string branchNameGoodHitTimeAvg     = "H.ladhod." + PLANE_NAMES[i] + ".GoodHitTimeAvg";
    std::string branchNameGoodHitTimeDiff    = "H.ladhod." + PLANE_NAMES[i] + ".GoodHitTimeDiff";
    std::string branchNameHodoHitEdep        = "H.ladhod." + PLANE_NAMES[i] + ".HodoHitEdep";

    tree->SetBranchAddress(branchNameTopTdcTimeRaw.c_str(), TopTdcTimeRaw[i]);
    tree->SetBranchAddress(branchNameTopTdcCounter.c_str(), TopTdcCounter[i]);
    tree->SetBranchAddress(branchNameTopAdcPulseTimeRaw.c_str(), TopAdcPulseTimeRaw[i]);
    tree->SetBranchAddress(branchNameTopAdcCounter.c_str(), TopAdcCounter[i]);
    tree->SetBranchAddress(branchNameBtmTdcTimeRaw.c_str(), BtmTdcTimeRaw[i]);
    tree->SetBranchAddress(branchNameBtmTdcCounter.c_str(), BtmTdcCounter[i]);
    tree->SetBranchAddress(branchNameBtmAdcPulseTimeRaw.c_str(), BtmAdcPulseTimeRaw[i]);
    tree->SetBranchAddress(branchNameBtmAdcCounter.c_str(), BtmAdcCounter[i]);
    tree->SetBranchAddress(branchNameGoodHitTimeAvg.c_str(), &GoodHitTimeAvg[i]);
    tree->SetBranchAddress(branchNameGoodHitTimeDiff.c_str(), &GoodHitTimeDiff[i]);
    tree->SetBranchAddress(branchNameHodoHitEdep.c_str(), &HodoHitEdep[i]);
  }
  // Set branch addresses for Ndata
  tree->SetBranchAddress("Ndata.H.ladhod.TopTdcTimeRaw", &NdataTopTdcTime);
  tree->SetBranchAddress("Ndata.H.ladhod.TopAdcPulseTimeRaw", &NdataTopAdcPulseTimeRaw);
  tree->SetBranchAddress("Ndata.H.ladhod.BtmTdcTimeRaw", &NdataBtmTdcTime);
  tree->SetBranchAddress("Ndata.H.ladhod.BtmAdcPulseTimeRaw", &NdataBtmAdcPulseTimeRaw);

  // Define Histograms
  TH1D *histTopTdcTimeRaw[NPLANES][NBARS];
  TH1D *histTopAdcPulseTimeRaw[NPLANES][NBARS];
  TH1D *histBtmTdcTimeRaw[NPLANES][NBARS];
  TH1D *histBtmAdcPulseTimeRaw[NPLANES][NBARS];
  TH1D *histGoodHitTimeAvg[NPLANES][NBARS];
  TH1D *histGoodHitTimeDiff[NPLANES][NBARS];
  TH1D *histHodoHitEdep[NPLANES][NBARS];

  for (int i = 0; i < NPLANES; ++i) {
    for (int j = 0; j < NBARS; ++j) {
      histTopTdcTimeRaw[i][j] =
        new TH1D(("histTopTdcTimeRaw_" + PLANE_NAMES[i] + "_Bar" + to_string(j)).c_str(),
             ("Top Raw TDC Time for Plane " + PLANE_NAMES[i] + " Bar " + to_string(j)).c_str(),
             hr_TDC_raw.nxbins, hr_TDC_raw.xmin, hr_TDC_raw.xmax);
      histTopAdcPulseTimeRaw[i][j] =
        new TH1D(("histTopAdcPulseTimeRaw_" + PLANE_NAMES[i] + "_Bar" + to_string(j)).c_str(),
             ("Top ADC Pulse Time Raw for Plane " + PLANE_NAMES[i] + " Bar " + to_string(j)).c_str(),
             hr_ADC_raw.nxbins, hr_ADC_raw.xmin, hr_ADC_raw.xmax);
      histBtmTdcTimeRaw[i][j] =
        new TH1D(("histBtmTdcTimeRaw_" + PLANE_NAMES[i] + "_Bar" + to_string(j)).c_str(),
             ("Bottom Raw TDC Time for Plane " + PLANE_NAMES[i] + " Bar " + to_string(j)).c_str(),
             hr_TDC_raw.nxbins, hr_TDC_raw.xmin, hr_TDC_raw.xmax);
      histBtmAdcPulseTimeRaw[i][j] =
        new TH1D(("histBtmAdcPulseTimeRaw_" + PLANE_NAMES[i] + "_Bar" + to_string(j)).c_str(),
             ("Bottom ADC Pulse Time Raw for Plane " + PLANE_NAMES[i] + " Bar " + to_string(j)).c_str(),
             hr_ADC_raw.nxbins, hr_ADC_raw.xmin, hr_ADC_raw.xmax);
      histGoodHitTimeAvg[i][j] =
        new TH1D(("histGoodHitTimeAvg_" + PLANE_NAMES[i] + "_Bar" + to_string(j)).c_str(),
             ("Good Hit Time Avg for Plane " + PLANE_NAMES[i] + " Bar " + to_string(j)).c_str(),
             hr_TDC_time.nxbins, hr_TDC_time.xmin, hr_TDC_time.xmax);
      histGoodHitTimeDiff[i][j] =
        new TH1D(("histGoodHitTimeDiff_" + PLANE_NAMES[i] + "_Bar" + to_string(j)).c_str(),
             ("Good Hit Time Diff for Plane " + PLANE_NAMES[i] + " Bar " + to_string(j)).c_str(),
             hr_ADC_Difftime.nxbins, hr_ADC_Difftime.xmin, hr_ADC_Difftime.xmax);
      histHodoHitEdep[i][j] =
        new TH1D(("histHodoHitEdep_" + PLANE_NAMES[i] + "_Bar" + to_string(j)).c_str(),
             ("Hodo Hit Energy Deposition for Plane " + PLANE_NAMES[i] + " Bar " + to_string(j)).c_str(),
             hr_Edep.nxbins, hr_Edep.xmin, hr_Edep.xmax);
    }
  }

  // Loop over all entries in the tree
  Long64_t nEntries = tree->GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);
    for (int plane = 0; plane < NPLANES; ++plane) {
      // Fill raw data histograms
      for (int j = 0; j < NdataTopTdcTime; ++j) {
        histTopTdcTimeRaw[plane][static_cast<int>(TopTdcCounter[plane][j])]->Fill(TopTdcTimeRaw[plane][j]);
      }
      for (int j = 0; j < NdataTopAdcPulseTimeRaw; ++j) {
        histTopAdcPulseTimeRaw[plane][static_cast<int>(TopAdcCounter[plane][j])]->Fill(TopAdcPulseTimeRaw[plane][j]);
      }
      for (int j = 0; j < NdataBtmTdcTime; ++j) {
        histBtmTdcTimeRaw[plane][static_cast<int>(BtmTdcCounter[plane][j])]->Fill(BtmTdcTimeRaw[plane][j]);
      }
      for (int j = 0; j < NdataBtmAdcPulseTimeRaw; ++j) {
        histBtmAdcPulseTimeRaw[plane][static_cast<int>(BtmAdcCounter[plane][j])]->Fill(BtmAdcPulseTimeRaw[plane][j]);
      }
      for (int bar = 0; bar < NBARS; ++bar) {
        histGoodHitTimeAvg[plane][bar]->Fill(GoodHitTimeAvg[plane][bar]);
        histGoodHitTimeDiff[plane][bar]->Fill(GoodHitTimeDiff[plane][bar]);
        histHodoHitEdep[plane][bar]->Fill(HodoHitEdep[plane][bar]);
      }
    }
  }

  //Make histograms for each plane
  // Create histograms for each plane by summing over bars
  TH1D *histPlaneTopTdcTimeRaw[NPLANES];
  TH1D *histPlaneTopAdcPulseTimeRaw[NPLANES];
  TH1D *histPlaneBtmTdcTimeRaw[NPLANES];
  TH1D *histPlaneBtmAdcPulseTimeRaw[NPLANES];
  TH1D *histPlaneGoodHitTimeAvg[NPLANES];
  TH1D *histPlaneGoodHitTimeDiff[NPLANES];
  TH1D *histPlaneHodoHitEdep[NPLANES];

  for (int plane = 0; plane < NPLANES; ++plane) {
    histPlaneTopTdcTimeRaw[plane] =
      new TH1D(("histPlaneTopTdcTimeRaw_" + PLANE_NAMES[plane]).c_str(),
               ("Top Raw TDC Time for Plane " + PLANE_NAMES[plane]).c_str(),
               hr_TDC_raw.nxbins, hr_TDC_raw.xmin, hr_TDC_raw.xmax);
    histPlaneTopAdcPulseTimeRaw[plane] =
      new TH1D(("histPlaneTopAdcPulseTimeRaw_" + PLANE_NAMES[plane]).c_str(),
               ("Top ADC Pulse Time Raw for Plane " + PLANE_NAMES[plane]).c_str(),
               hr_ADC_raw.nxbins, hr_ADC_raw.xmin, hr_ADC_raw.xmax);
    histPlaneBtmTdcTimeRaw[plane] =
      new TH1D(("histPlaneBtmTdcTimeRaw_" + PLANE_NAMES[plane]).c_str(),
               ("Bottom Raw TDC Time for Plane " + PLANE_NAMES[plane]).c_str(),
               hr_TDC_raw.nxbins, hr_TDC_raw.xmin, hr_TDC_raw.xmax);
    histPlaneBtmAdcPulseTimeRaw[plane] =
      new TH1D(("histPlaneBtmAdcPulseTimeRaw_" + PLANE_NAMES[plane]).c_str(),
               ("Bottom ADC Pulse Time Raw for Plane " + PLANE_NAMES[plane]).c_str(),
               hr_ADC_raw.nxbins, hr_ADC_raw.xmin, hr_ADC_raw.xmax);
    histPlaneGoodHitTimeAvg[plane] =
      new TH1D(("histPlaneGoodHitTimeAvg_" + PLANE_NAMES[plane]).c_str(),
               ("Good Hit Time Avg for Plane " + PLANE_NAMES[plane]).c_str(),
               hr_TDC_time.nxbins, hr_TDC_time.xmin, hr_TDC_time.xmax);
    histPlaneGoodHitTimeDiff[plane] =
      new TH1D(("histPlaneGoodHitTimeDiff_" + PLANE_NAMES[plane]).c_str(),
               ("Good Hit Time Diff for Plane " + PLANE_NAMES[plane]).c_str(),
               hr_ADC_Difftime.nxbins, hr_ADC_Difftime.xmin, hr_ADC_Difftime.xmax);
    histPlaneHodoHitEdep[plane] =
      new TH1D(("histPlaneHodoHitEdep_" + PLANE_NAMES[plane]).c_str(),
               ("Hodo Hit Energy Deposition for Plane " + PLANE_NAMES[plane]).c_str(),
               hr_Edep.nxbins, hr_Edep.xmin, hr_Edep.xmax);

    // Sum over bars for each plane
    for (int bar = 0; bar < NBARS; ++bar) {
      histPlaneTopTdcTimeRaw[plane]->Add(histTopTdcTimeRaw[plane][bar]);
      histPlaneTopAdcPulseTimeRaw[plane]->Add(histTopAdcPulseTimeRaw[plane][bar]);
      histPlaneBtmTdcTimeRaw[plane]->Add(histBtmTdcTimeRaw[plane][bar]);
      histPlaneBtmAdcPulseTimeRaw[plane]->Add(histBtmAdcPulseTimeRaw[plane][bar]);
      histPlaneGoodHitTimeAvg[plane]->Add(histGoodHitTimeAvg[plane][bar]);
      histPlaneGoodHitTimeDiff[plane]->Add(histGoodHitTimeDiff[plane][bar]);
      histPlaneHodoHitEdep[plane]->Add(histHodoHitEdep[plane][bar]);
    }
  }

  // Save histograms to a new ROOT file
  TFile *outputFile = new TFile("lad_timing_histos.root", "RECREATE");
  if (!outputFile || outputFile->IsZombie()) {
    std::cerr << "Error: Cannot create output file lad_timing_histos.root" << std::endl;
    return;
  }
  // Write histograms to the output file
  // Create directories for each histogram type and subdirectories for each plane
  outputFile->cd();
  for (int plane = 0; plane < NPLANES; ++plane) {
    // Create subdirectories for each plane
    std::string planeDirName = "Plane_" + PLANE_NAMES[plane];
    outputFile->mkdir(planeDirName.c_str());
    outputFile->cd(planeDirName.c_str());

    // Write histograms for the current plane
    histPlaneTopTdcTimeRaw[plane]->Write();
    histPlaneTopAdcPulseTimeRaw[plane]->Write();
    histPlaneBtmTdcTimeRaw[plane]->Write();
    histPlaneBtmAdcPulseTimeRaw[plane]->Write();
    histPlaneGoodHitTimeAvg[plane]->Write();
    histPlaneGoodHitTimeDiff[plane]->Write();
    histPlaneHodoHitEdep[plane]->Write();

    // Create subdirectories for each bar within the plane
    for (int bar = 0; bar < NBARS; ++bar) {
      std::string barDirName = "Bar_" + std::to_string(bar);
      outputFile->mkdir(barDirName.c_str());
      outputFile->cd(barDirName.c_str());

      // Write histograms for the current bar
      histTopTdcTimeRaw[plane][bar]->Write();
      histTopAdcPulseTimeRaw[plane][bar]->Write();
      histBtmTdcTimeRaw[plane][bar]->Write();
      histBtmAdcPulseTimeRaw[plane][bar]->Write();
      histGoodHitTimeAvg[plane][bar]->Write();
      histGoodHitTimeDiff[plane][bar]->Write();
      histHodoHitEdep[plane][bar]->Write();

      // Return to the plane directory
      outputFile->cd("..");
    }

    // Return to the root directory
    outputFile->cd("..");
  }

  // Write and close the output file
  outputFile->Write();
  outputFile->Close();
  // Close the file
  file->Close();
}