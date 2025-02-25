// Description: Reads in replayed root file, and creates histograms of the waveforms for each. Currently works for lad hodoscope data.
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <iostream>

using namespace ROOT;
using namespace std;

const int MAX_WAVEFORM_LENGTH = 5000;
const int NUM_EVENTS_TO_USE   = 50; // Number of events to show the events for. -1 for all of them.

static const Double_t adcDynamicRange = 1000.0;                     // Units of mV
static const Double_t nAdcChan        = 4096.0;                     // Units of ADC channels
static const Double_t adcChanTomV     = adcDynamicRange / nAdcChan; // Units of mV/ADC Chan
static const Double_t adcImpedence    = 50.0;                       // FADC input impedence in units of Ohms
static const Double_t adcTimeSample   = 4.0;                        // Length of FADC time sample in units of ns
static const Double_t adcTimeRes      = 0.0625;                     // FADC time resolution in units of ns
static const Double_t adcChanTopC     = (adcDynamicRange / 1000 / nAdcChan) * (adcTimeSample / adcImpedence);
// (1000 mV / 4096 adc channels) * (4 ns time sample / 50 ohms input resistance) = ~0.020 pc/channel

void cosmic_waveform_histos(int run_number) {
  // Open the input file
  string input_string =
      "/volatile/hallc/c-lad/ehingerl/ROOTfiles/COSMICS/LAD_cosmic_hall_" + to_string(run_number) + "_-1.root";
  TFile *inputFile = TFile::Open(input_string.c_str(), "READ");

  // Get the TTree from the input file
  TTree *tree = dynamic_cast<TTree *>(inputFile->Get("T"));

  // Create a new output file
  string output_string = "../histos/waveform_histos_" + to_string(run_number) + "_output.root";
  TFile *outputFile    = TFile::Open(output_string.c_str(), "RECREATE");
  // Create directories for each type of histogram
  TDirectory *indivEventWaveform = outputFile->mkdir("Indiv_Event_Waveform");

  Double_t Btm_ADC_Waveform[MAX_WAVEFORM_LENGTH] = {0};
  Double_t Top_ADC_Waveform[MAX_WAVEFORM_LENGTH] = {0};

  // Create directories for each type of histogram
  tree->SetBranchAddress("L.hod.200.adcBtmSampWaveform", &Btm_ADC_Waveform);
  tree->SetBranchAddress("L.hod.200.adcTopSampWaveform", &Top_ADC_Waveform);

  // Loop over the TTree entries
  for (int i_evt = 0; i_evt < tree->GetEntries(); i_evt++) {
    tree->GetEntry(i_evt);

    // if (num_fired_pmts == 16 && !plottedWaveform) {
    if (i_evt < NUM_EVENTS_TO_USE ||
        NUM_EVENTS_TO_USE == -1) { // this is dumb inside a for loop, but different conditions were used to look for
                                   // specific "bad events"
      TDirectory *eventWaveformDir = indivEventWaveform->mkdir(Form("Event_%d", i_evt));
      eventWaveformDir->cd();

      int indx = 0;
      while (Btm_ADC_Waveform[indx] != 0) {
        int btm_paddle_indx = (int)Btm_ADC_Waveform[indx++] - 1;
        int n_samples       = (int)Btm_ADC_Waveform[indx++];
        TH1D h1_btm_waveform(Form("Btm_Waveform_%d", btm_paddle_indx), "Btm_Waveform; Time [ns]; ADC Value [mV]",
                             n_samples, 0, n_samples * adcTimeSample);
        for (int i = 0; i < n_samples; i++) {
          h1_btm_waveform.Fill(i * adcTimeSample, Btm_ADC_Waveform[indx + i] * adcChanTomV);
        }
        h1_btm_waveform.Write();
        indx += n_samples;
      }
      indx = 0;
      while (Top_ADC_Waveform[indx] != 0) {
        int top_paddle_indx = (int)Top_ADC_Waveform[indx++] - 1;
        int n_samples       = (int)Top_ADC_Waveform[indx++];
        TH1D h1_top_waveform(Form("Top_Waveform_%d", top_paddle_indx), "Top_Waveform; Time [ns]; ADC Value [mV]",
                             n_samples, 0, n_samples * adcTimeSample);
        for (int i = 0; i < n_samples; i++) {
          h1_top_waveform.Fill(i * adcTimeSample, Top_ADC_Waveform[indx + i] * adcChanTomV);
        }
        h1_top_waveform.Write();
        indx += n_samples;
      }
    }
    if (i_evt % 10000 == 0) {
      cout << "Processed " << i_evt << " events" << endl;
    }
  }

  // Close the output file
  outputFile->Close();
  // Close the input file
  inputFile->Close();
}