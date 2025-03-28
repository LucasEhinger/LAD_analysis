#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TString.h>
#include <cmath>
#include <vector>

void saveHistogramsToPDF(const char *inputFileName, const char *outputFileName,
                         const std::vector<TString> &histogramNames, int nCols, int nRows) {
  TFile *inputFile = TFile::Open(inputFileName, "READ");
  if (!inputFile || inputFile->IsZombie()) {
    printf("Error: Cannot open input file %s\n", inputFileName);
    return;
  }

  TCanvas canvas("canvas", "Histograms", 800, 600);
  canvas.Divide(nCols, nRows); // Dynamic layout based on plotsPerPage

  int padIndex = 1;
  canvas.Print((TString(outputFileName) + "[").Data()); // Open the PDF

  for (const auto &histName : histogramNames) {
    TH1 *hist = dynamic_cast<TH1 *>(inputFile->Get(histName));
    if (!hist) {
      printf("Warning: Histogram %s not found in file.\n", histName.Data());
      continue;
    }

    canvas.cd(padIndex);
    hist->Draw();

    if (padIndex == nCols * nRows) {
      canvas.Print(outputFileName); // Save the current page
      canvas.Clear();
      canvas.Divide(nCols, nRows);
      padIndex = 1;
    } else {
      padIndex++;
    }
  }

  if (padIndex != 1) {
    canvas.Print(outputFileName); // Save the last page if not full
  }

  canvas.Print((TString(outputFileName) + "]").Data()); // Close the PDF
  inputFile->Close();
}

void histos_to_pdf() {
  const char *inputFileName =
      "/u/home/ehingerl/hallc/software/lad_replay/ROOTfiles/COSMICS/LAD_wGEM_cosmic_hall_296_-1.root";
  const char *outputFileName          = "files/LAD_wGEM_cosmic_hall_296_-1_histos.pdf";
  std::vector<TString> histogramNames = {
      "h1_gem_Nlayers_hit",      "h1_gem_Nlayers_hitu",     "h1_gem_Nlayers_hitv",     "h1_gem_Nlayers_hituv",
      "h2_gem_NstripsU_layer",   "h2_gem_NstripsV_layer",   "h2_gem_NclustU_layer",    "h2_gem_NclustV_layer",
      "h1_gem_clustWidthU_0",    "h1_gem_clustWidthV_0",    "h1_gem_clustWidthU_1",    "h1_gem_clustWidthV_1",
      "h1_gem_clustSampMaxU_0",  "h1_gem_clustSampMaxV_0",  "h1_gem_clustSampMaxU_1",  "h1_gem_clustSampMaxV_1",
      "h1_gem_clustADCMaxU_0",   "h1_gem_clustADCMaxV_0",   "h1_gem_clustADCMaxU_1",   "h1_gem_clustADCMaxV_1",
      "h1_gem_clustADCSumU_0",   "h1_gem_clustADCSumV_0",   "h1_gem_clustADCSumU_1",   "h1_gem_clustADCSumV_1",
      "h1_gem_clustTimeMeanU_0", "h1_gem_clustTimeMeanV_0", "h1_gem_clustTimeMeanU_1", "h1_gem_clustTimeMeanV_1",
      "h2_gem_2dhit_0",          "h2_gem_2dhit_1",          "h1_gem_nhits_layer0",     "h1_gem_nhits_layer1",
      "h1_gem_time_0",           "h1_gem_time_1",           "h1_gem_ADCMean_0",        "h1_gem_ADCMean_1",
      "h1_gem_ADCAsym_0",        "h1_gem_ADCAsym_1",        "h1_gem_TimeDiff_0",       "h1_gem_TimeDiff_1",
      "h1_gem_TimeCorr_0",       "h1_gem_TimeCorr_1",       "h1_gem_ntracks",          "h1_gem_track_t",
      "h1_gem_track_dt",         "h1_gem_track_d0",         "h1_gem_stripsfiredU_m0",  "h1_gem_stripsfiredV_m0",
      "h1_gem_stripsfiredU_m1",  "h1_gem_stripsfiredV_m1",  "h2_gem_stripU_adc_0",       "h2_gem_stripU_adc_1",
      "h2_gem_stripV_adc_0",       "h2_gem_stripV_adc_1"};

  saveHistogramsToPDF(inputFileName, outputFileName, histogramNames, 2, 2);
}
