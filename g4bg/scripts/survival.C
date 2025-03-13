#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>

using namespace std;
using namespace ROOT;

void survival() {
  // Open the first ROOT file
  TFile *file1 = TFile::Open("file1.root");
  if (!file1 || file1->IsZombie()) {
    printf("Error opening file1.root\n");
    return;
  }

  // Open the second ROOT file
  TFile *file2 = TFile::Open("file2.root");
  if (!file2 || file2->IsZombie()) {
    printf("Error opening file2.root\n");
    return;
  }

  // Get the histogram from the first file
  TH1 *hist1 = (TH1*)file1->Get("histogram_name");
  if (!hist1) {
    printf("Error getting histogram from file1.root\n");
    return;
  }

  // Get the histogram from the second file
  TH1 *hist2 = (TH1*)file2->Get("histogram_name");
  if (!hist2) {
    printf("Error getting histogram from file2.root\n");
    return;
  }

  // Create a canvas to draw the histograms
  TCanvas *canvas = new TCanvas("canvas", "Overlay Histograms", 800, 600);

  // Draw the first histogram
  hist1->SetLineColor(kRed);
  hist1->Draw();

  // Draw the second histogram on top of the first
  hist2->SetLineColor(kBlue);
  hist2->Draw("SAME");

  // Update the canvas to show the histograms
  canvas->Update();

  // Save the canvas as an image file
  canvas->SaveAs("overlay_histograms.png");

  // Close the ROOT files
  file1->Close();
  file2->Close();
}