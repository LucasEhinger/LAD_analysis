#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <iostream>
#include <string>
#include <map>
#include <tuple>
#include <cmath>
#include <TMath.h>
#include <string>

using namespace std;

void check_lad_positions(string energy= "400") {
  // Open the ROOT file
  string filepath     = "/volatile/hallc/c-lad/ehingerl/G4_LAD/carlos_proton/";
  string infile_name  = filepath + "raw/" + "ScanLAD_proton_" + energy + "MeV_10k_20240205.root";

  TFile *file = TFile::Open(infile_name.c_str(), "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Error opening file" << std::endl;
    return;
  }

  // Get the tree from the file
  TTree *tree;
  file->GetObject("hodoposition", tree);
  if (!tree) {
    std::cerr << "Error getting tree" << std::endl;
    file->Close();
    return;
  }

  // Define variables to hold branch data
  std::vector<double> *vXbar = nullptr;
  std::vector<double> *vYbar = nullptr;
  std::vector<double> *vZbar = nullptr;
  std::vector<double> *vTbar = nullptr;
  std::vector<int> *vPaddle = nullptr;
  std::vector<int> *vTrackID = nullptr;

  // Set branch addresses
  tree->SetBranchAddress("vXbar", &vXbar);
  tree->SetBranchAddress("vYbar", &vYbar);
  tree->SetBranchAddress("vZbar", &vZbar);
  tree->SetBranchAddress("vTbar", &vTbar);
  tree->SetBranchAddress("vPaddle", &vPaddle);
  tree->SetBranchAddress("vTrackID", &vTrackID);

  // Create a map to store sums and counts for each paddle
  std::map<int, std::tuple<double, double, double, double, double, double, double, double, int>> paddleData;

  // Loop through the tree entries
  Long64_t nEntries = tree->GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);

    // Loop through the vectors
    for (size_t j = 0; j < vXbar->size(); ++j) {
      if (vTrackID->at(j) == 1) {
        int paddle = vPaddle->at(j);
        double x = vXbar->at(j);
        double y = vYbar->at(j);
        double z = vZbar->at(j);
        double t = vTbar->at(j);
        double radius = sqrt(x * x + y * y);
        double xz = sqrt(x * x + z * z);
        double theta = TMath::ATan2(TMath::Sqrt(x * x + y * y), z);
        double phi = TMath::ATan2(y, x);

        auto &data = paddleData[paddle];
        std::get<0>(data) += x;
        std::get<1>(data) += y;
        std::get<2>(data) += z;
        std::get<3>(data) += radius;
        std::get<4>(data) += xz;
        std::get<5>(data) += t;
        std::get<6>(data) += theta;
        std::get<7>(data) += phi;
        std::get<8>(data) += 1;
      }
    }
  }

  // Calculate and print the averages for each paddle
  for (const auto &entry : paddleData) {
    int paddle = entry.first;
    const auto &data = entry.second;
    double avgX = std::get<0>(data) / std::get<8>(data);
    double avgY = std::get<1>(data) / std::get<8>(data);
    double avgZ = std::get<2>(data) / std::get<8>(data);
    double avgRadius = std::get<3>(data) / std::get<8>(data);
    double avgXZ = std::get<4>(data) / std::get<8>(data);
    double avgT = std::get<5>(data) / std::get<8>(data);
    double avgTheta = std::get<6>(data) / std::get<8>(data);
    double avgPhi = std::get<7>(data) / std::get<8>(data);

    std::cout << "Paddle: " << paddle
              << ", Avg X: " << avgX
              << ", Avg Y: " << avgY
              << ", Avg Z: " << avgZ
              << ", Avg Radius: " << avgRadius
              << ", Avg sqrt(X*X + Z*Z): " << avgXZ
              << ", Avg T: " << avgT
              << ", Avg Theta: " << avgTheta
              << ", Avg Phi: " << avgPhi
              << std::endl;
  }

  // Close the file
  file->Close();
}