#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <iostream>

using namespace ROOT;
using namespace std;

const int NUM_BARS   = 8;
const int N_DATA_MAX = 100;

const int MAX_WAVEFORM_LENGTH = 5000;

static const int PED_N_BINS       = 50;
static const int PED_MIN          = 0;   // mV
static const int PED_MAX          = 100; // mV
static const int AMP_N_BINS       = 50;
static const int AMP_MIN          = 0;    // mV
static const int AMP_MAX          = 2000; // mV
static const int INTEGRAL_N_BINS  = 50;
static const int INTEGRAL_MIN     = 0;    // pC
static const int INTEGRAL_MAX     = 1000; // pC
static const int TIME_N_BINS      = 500;
static const int TIME_MIN         = 0;    // ns
static const int TIME_MAX         = 1000; // ns
static const int TIME_DIFF_N_BINS = 50;
static const int TIME_DIFF_MIN    = -100; // ns
static const int TIME_DIFF_MAX    = 100;  // ns
static const int TIME_WINDOW_MIN  = 000;  // ns
static const int TIME_WINDOW_MAX  = 2000; // ns

static const Double_t adcDynamicRange = 1000.0;                     // Units of mV
static const Double_t nAdcChan        = 4096.0;                     // Units of ADC channels
static const Double_t adcChanTomV     = adcDynamicRange / nAdcChan; // Units of mV/ADC Chan
static const Double_t adcImpedence    = 50.0;                       // FADC input impedence in units of Ohms
static const Double_t adcTimeSample   = 4.0;                        // Length of FADC time sample in units of ns
static const Double_t adcTimeRes      = 0.0625;                     // FADC time resolution in units of ns
static const Double_t adcChanTopC     = (adcDynamicRange / 1000 / nAdcChan) * (adcTimeSample / adcImpedence);
// (1000 mV / 4096 adc channels) * (4 ns time sample / 50 ohms input resistance) = ~0.020 pc/channel

static const int MIN_NUM_PMTS_FIRED = 10;
static const int MAX_NOISY_PMTS     = 22;

void cosmic_histos(int run_number) {
  // int run_number = 70;
  // int run_number = 22;

  // Open the input file
  string input_string =
      "/volatile/hallc/c-lad/ehingerl/ROOTfiles/COSMICS/LAD_cosmic_" + to_string(run_number) + "_-1.root";
  TFile *inputFile = TFile::Open(input_string.c_str(), "READ");

  // Get the TTree from the input file
  TTree *tree = dynamic_cast<TTree *>(inputFile->Get("T"));

  // Create a new output file
  string output_string = "cosmic_histos_" + to_string(run_number) + "_output.root";
  TFile *outputFile    = TFile::Open(output_string.c_str(), "RECREATE");
  // Create directories for each type of histogram
  TDirectory *pedDir             = outputFile->mkdir("Ped");
  TDirectory *pulseIntDir        = outputFile->mkdir("Pulse_Int");
  TDirectory *eventIntDir        = pulseIntDir->mkdir("Event_Int");
  TDirectory *pulseAmpDir        = outputFile->mkdir("Pulse_Amp");
  TDirectory *eventAmpDir        = pulseAmpDir->mkdir("Event_Amp");
  TDirectory *pulseTimeDir       = outputFile->mkdir("Pulse_Time");
  TDirectory *topBtmCompDir      = outputFile->mkdir("2D_Top_Btm_Comparison");
  TDirectory *multDir            = outputFile->mkdir("Mult");
  TDirectory *indivEventDir      = outputFile->mkdir("Indiv_Event");
  TDirectory *indivEventWaveform = outputFile->mkdir("Indiv_Event_Waveform");
  TDirectory *pulseTimeDiffDir   = outputFile->mkdir("Pulse_Time_Diff");
  TDirectory *intAmpCompDir      = outputFile->mkdir("2D_Int_Amp_Comparison");
  TDirectory *noisyPMTDir        = outputFile->mkdir("Noisy_PMTs");

  // Create histograms
  TH1F *h1_Btm_Ped[NUM_BARS];
  TH1F *h1_Btm_Pulse_Int[NUM_BARS];
  TH1F *h1_Btm_Pulse_Amp[NUM_BARS];
  TH1F *h1_Btm_Pulse_Time[NUM_BARS];

  TH1F *h1_Top_Ped[NUM_BARS];
  TH1F *h1_Top_Pulse_Int[NUM_BARS];
  TH1F *h1_Top_Pulse_Amp[NUM_BARS];
  TH1F *h1_Top_Pulse_Time[NUM_BARS];

  TH1F *h1_Btm_Evt_Int[NUM_BARS];
  TH1F *h1_Top_Evt_Int[NUM_BARS];
  TH1F *h1_Btm_Evt_Amp[NUM_BARS];
  TH1F *h1_Top_Evt_Amp[NUM_BARS];

  TH1F *h1_Top_Mult[NUM_BARS];
  TH1F *h1_Btm_Mult[NUM_BARS];
  TH1F *h1_Tot_PMTsFired;

  TH2D *h2_Pulse_Int[NUM_BARS];
  TH2D *h2_Pulse_Amp[NUM_BARS];
  TH2D *h2_Pulse_Time[NUM_BARS];
  TH2D *h2_Btm_Pulse_Int_Amp[NUM_BARS];
  TH2D *h2_Top_Pulse_Int_Amp[NUM_BARS];
  // TH2D *h2_Btm_Evt_Int_Amp[NUM_BARS];
  // TH2D *h2_Top_Evt_Int_Amp[NUM_BARS];

  TH2D *h2_Max_Adc_Time_Diff;
  TH2D *h2_Max_Adc_Time_Diff_Normalized;
  TH1F *h1_Max_Adc_Time_Diff[NUM_BARS];

  TH1F *h1_Btm_Noisy_PMTs;
  TH1F *h1_Top_Noisy_PMTs;

  for (int i = 0; i < NUM_BARS; i++) {
    h1_Btm_Ped[i] =
        new TH1F(Form("Btm_Ped_%d", i), Form("Btm_Ped_%d; Pulse Ped; Counts", i), PED_N_BINS, PED_MIN, PED_MAX);
    h1_Btm_Pulse_Int[i]  = new TH1F(Form("Btm_Pulse_Int_%d", i), Form("Btm_Pulse_Int_%d; Pulse Int [pC]; Counts", i),
                                    INTEGRAL_N_BINS, INTEGRAL_MIN, INTEGRAL_MAX);
    h1_Btm_Pulse_Amp[i]  = new TH1F(Form("Btm_Pulse_Amp_%d", i), Form("Btm_Pulse_Amp_%d; Pulse Amp [mV]; Counts", i),
                                    AMP_N_BINS, AMP_MIN, AMP_MAX);
    h1_Btm_Pulse_Time[i] = new TH1F(Form("Btm_Pulse_Time_%d", i), Form("Btm_Pulse_Time_%d; Pulse Time [ns]; Counts", i),
                                    TIME_N_BINS, TIME_MIN, TIME_MAX);

    h1_Top_Ped[i] =
        new TH1F(Form("Top_Ped_%d", i), Form("Top_Ped_%d; Pulse Ped; Counts", i), PED_N_BINS, PED_MIN, PED_MAX);
    h1_Top_Pulse_Int[i]  = new TH1F(Form("Top_Pulse_Int_%d", i), Form("Top_Pulse_Int_%d; Pulse Int [pC]; Counts", i),
                                    INTEGRAL_N_BINS, INTEGRAL_MIN, INTEGRAL_MAX);
    h1_Top_Pulse_Amp[i]  = new TH1F(Form("Top_Pulse_Amp_%d", i), Form("Top_Pulse_Amp_%d; Pulse Amp [mV]; Counts", i),
                                    AMP_N_BINS, AMP_MIN, AMP_MAX);
    h1_Top_Pulse_Time[i] = new TH1F(Form("Top_Pulse_Time_%d", i), Form("Top_Pulse_Time_%d; Pulse Time [ns]; Counts", i),
                                    TIME_N_BINS, TIME_MIN, TIME_MAX);

    h1_Btm_Evt_Int[i] = new TH1F(Form("Btm_Evt_Int_%d", i), Form("Btm_Evt_Int_%d; Pulse Int [pC]; Counts", i),
                                 INTEGRAL_N_BINS, INTEGRAL_MIN, INTEGRAL_MAX);
    h1_Top_Evt_Int[i] = new TH1F(Form("Top_Evt_Int_%d", i), Form("Top_Evt_Int_%d; Pulse Int [pC]; Counts", i),
                                 INTEGRAL_N_BINS, INTEGRAL_MIN, INTEGRAL_MAX);
    h1_Btm_Evt_Amp[i] = new TH1F(Form("Btm_Evt_Amp_%d", i), Form("Btm_Evt_Amp_%d; Pulse Amp [mV]; Counts", i),
                                 AMP_N_BINS, AMP_MIN, AMP_MAX);
    h1_Top_Evt_Amp[i] = new TH1F(Form("Top_Evt_Amp_%d", i), Form("Top_Evt_Amp_%d; Pulse Amp [mV]; Counts", i),
                                 AMP_N_BINS, AMP_MIN, AMP_MAX);

    h1_Top_Mult[i] = new TH1F(Form("Top_Mult_%d", i), Form("Top_Mult_%d; Paddle Multiplicity; Counts", i), NUM_BARS + 1,
                              -0.5, NUM_BARS + 0.5);
    h1_Btm_Mult[i] = new TH1F(Form("Btm_Mult_%d", i), Form("Btm_Mult_%d; Paddle Multiplicity; Counts", i), NUM_BARS + 1,
                              -0.5, NUM_BARS + 0.5);

    h2_Pulse_Int[i]  = new TH2D(Form("Pulse_Int_%d", i), Form("Pulse_Int_%d; Btm [pC]; Top [pC]", i), INTEGRAL_N_BINS,
                                INTEGRAL_MIN, INTEGRAL_MAX, INTEGRAL_N_BINS, INTEGRAL_MIN, INTEGRAL_MAX);
    h2_Pulse_Amp[i]  = new TH2D(Form("Pulse_Amp_%d", i), Form("Pulse_Amp_%d; Btm [pC]; Top [pC]", i), AMP_N_BINS,
                                AMP_MIN, AMP_MAX, AMP_N_BINS, AMP_MIN, AMP_MAX);
    h2_Pulse_Time[i] = new TH2D(Form("Pulse_Time_%d", i), Form("Pulse_Time_%d; Btm [ns]; Top [ns]", i), TIME_N_BINS,
                                TIME_MIN, TIME_MAX, TIME_N_BINS, TIME_MIN, TIME_MAX);

    h2_Btm_Pulse_Int_Amp[i] =
        new TH2D(Form("Btm_Pulse_Int_Amp_%d", i), Form("Btm_Pulse_Int_Amp_%d; Int  [pC]; Amp [mV]", i), INTEGRAL_N_BINS,
                 INTEGRAL_MIN, INTEGRAL_MAX, AMP_N_BINS, AMP_MIN, AMP_MAX);
    h2_Top_Pulse_Int_Amp[i] =
        new TH2D(Form("Top_Pulse_Int_Amp_%d", i), Form("Top_Pulse_Int_Amp_%d; Int [pC]; Amp [mV]", i), INTEGRAL_N_BINS,
                 INTEGRAL_MIN, INTEGRAL_MAX, AMP_N_BINS, AMP_MIN, AMP_MAX);
    // h2_Btm_Evt_Int_Amp[i]   = new TH2D(Form("Btm_Evt_Int_Amp_%d", i), Form("Btm_Evt_Int_Amp_%d; Int; Amp", i),
    //                                    INTEGRAL_N_BINS, INTEGRAL_MIN, INTEGRAL_MAX, AMP_N_BINS, AMP_MIN, AMP_MAX);
    // h2_Top_Evt_Int_Amp[i]   = new TH2D(Form("Top_Evt_Int_Amp_%d", i), Form("Top_Evt_Int_Amp_%d; Int; Amp", i),
    //                                    INTEGRAL_N_BINS, INTEGRAL_MIN, INTEGRAL_MAX, AMP_N_BINS, AMP_MIN, AMP_MAX);

    h1_Max_Adc_Time_Diff[i] =
        new TH1F(Form("Max_Adc_Time_Diff_%d", i), Form("Adc_Time_Diff_%d; t_top - t_btm [ns]; Counts", i),
                 TIME_DIFF_N_BINS, TIME_DIFF_MIN, TIME_DIFF_MAX);
  }
  h1_Tot_PMTsFired =
      new TH1F("Tot_PMTs_Fired", "Tot_PMTs_Fired; Num PMTs Fired; Counts", NUM_BARS * 2 + 1, -0.5, NUM_BARS * 2 + 0.5);
  h2_Max_Adc_Time_Diff = new TH2D("Max_Adc_Time_Diff", "Adc_Time_Diff; Bar Number; t_top - t_btm [ns]", NUM_BARS, -0.5,
                                  NUM_BARS - 0.5, TIME_DIFF_N_BINS, TIME_DIFF_MIN, TIME_DIFF_MAX);
  h2_Max_Adc_Time_Diff_Normalized = new TH2D(
      "Max_Adc_Time_Diff_Normalized", "Adc_Time_Diff_Normalized; Bar Number; t_top - t_btm - #Delta t_avg [ns]",
      NUM_BARS, -0.5, NUM_BARS - 0.5, TIME_DIFF_N_BINS, TIME_DIFF_MIN, TIME_DIFF_MAX);

  h1_Btm_Noisy_PMTs = new TH1F("Btm_Noisy_PMTs", "Btm_Noisy_PMTs; Bar Number; Counts", NUM_BARS, -0.5, NUM_BARS - 0.5);
  h1_Top_Noisy_PMTs = new TH1F("Top_Noisy_PMTs", "Top_Noisy_PMTs; Bar Number; Counts", NUM_BARS, -0.5, NUM_BARS - 0.5);

  // Loop over the TTree entries
  Double_t Btm_ADC_Counter[N_DATA_MAX];
  Double_t Btm_ADC_Ped[N_DATA_MAX];
  Double_t Btm_ADC_Pulse_Int[N_DATA_MAX];
  Double_t Btm_ADC_Pulse_Amp[N_DATA_MAX];
  Double_t Btm_ADC_Pulse_Time[N_DATA_MAX];

  Double_t Top_ADC_Counter[N_DATA_MAX];
  Double_t Top_ADC_Ped[N_DATA_MAX];
  Double_t Top_ADC_Pulse_Int[N_DATA_MAX];
  Double_t Top_ADC_Pulse_Amp[N_DATA_MAX];
  Double_t Top_ADC_Pulse_Time[N_DATA_MAX];

  Double_t Btm_ADC_Waveform[MAX_WAVEFORM_LENGTH] = {0};
  Double_t Top_ADC_Waveform[MAX_WAVEFORM_LENGTH] = {0};

  Double_t Btm_ADC_Int_Max[NUM_BARS]  = {0};
  Double_t Top_ADC_Int_Max[NUM_BARS]  = {0};
  Double_t Btm_ADC_Amp_Max[NUM_BARS]  = {0};
  Double_t Top_ADC_Amp_Max[NUM_BARS]  = {0};
  Double_t Btm_ADC_Time_Max[NUM_BARS] = {0};
  Double_t Top_ADC_Time_Max[NUM_BARS] = {0};

  Double_t Btm_ADC_Time_Avg[NUM_BARS] = {0};
  Double_t Top_ADC_Time_Avg[NUM_BARS] = {0};

  int top_mult[NUM_BARS] = {0};
  int btm_mult[NUM_BARS] = {0};

  int n_data_top;
  int n_data_btm;
  bool plottedWaveform = false;

  tree->SetBranchAddress("L.hod.000.BtmAdcCounter", Btm_ADC_Counter);
  tree->SetBranchAddress("L.hod.000.BtmAdcPed", Btm_ADC_Ped);
  tree->SetBranchAddress("L.hod.000.BtmAdcPulseInt", Btm_ADC_Pulse_Int);
  tree->SetBranchAddress("L.hod.000.BtmAdcPulseAmp", Btm_ADC_Pulse_Amp);
  tree->SetBranchAddress("L.hod.000.BtmAdcPulseTime", Btm_ADC_Pulse_Time);

  tree->SetBranchAddress("L.hod.000.TopAdcCounter", Top_ADC_Counter);
  tree->SetBranchAddress("L.hod.000.TopAdcPed", Top_ADC_Ped);
  tree->SetBranchAddress("L.hod.000.TopAdcPulseInt", Top_ADC_Pulse_Int);
  tree->SetBranchAddress("L.hod.000.TopAdcPulseAmp", Top_ADC_Pulse_Amp);
  tree->SetBranchAddress("L.hod.000.TopAdcPulseTime", Top_ADC_Pulse_Time);

  tree->SetBranchAddress("Ndata.L.hod.000.BtmAdcCounter", &n_data_btm);
  tree->SetBranchAddress("Ndata.L.hod.000.TopAdcCounter", &n_data_top);

  tree->SetBranchAddress("L.hod.000.adcBtmSampWaveform", &Btm_ADC_Waveform);
  tree->SetBranchAddress("L.hod.000.adcTopSampWaveform", &Top_ADC_Waveform);

  // Get average times of top and bottom
  int n_ev_used_top[NUM_BARS] = {0};
  int n_ev_used_btm[NUM_BARS] = {0};
  for (int i_evt = 0; i_evt < tree->GetEntries(); i_evt++) {
    tree->GetEntry(i_evt);
    for (int i = 0; i < NUM_BARS; i++) {
      Btm_ADC_Int_Max[i]  = 0;
      Top_ADC_Int_Max[i]  = 0;
      Btm_ADC_Amp_Max[i]  = 0;
      Top_ADC_Amp_Max[i]  = 0;
      Btm_ADC_Time_Max[i] = 0;
      Top_ADC_Time_Max[i] = 0;
    }

    for (int i = 0; i < n_data_btm; i++) {
      int paddle_indx = Btm_ADC_Counter[i] - 1;
      if (Btm_ADC_Pulse_Int[i] > Btm_ADC_Int_Max[paddle_indx]) {
        Btm_ADC_Int_Max[paddle_indx]  = Btm_ADC_Pulse_Int[i];
        Btm_ADC_Amp_Max[paddle_indx]  = Btm_ADC_Pulse_Amp[i];
        Btm_ADC_Time_Max[paddle_indx] = Btm_ADC_Pulse_Time[i];
      }
    }
    for (int i = 0; i < n_data_top; i++) {
      int paddle_indx = Top_ADC_Counter[i] - 1;
      if (Top_ADC_Pulse_Int[i] > Top_ADC_Int_Max[paddle_indx]) {
        Top_ADC_Int_Max[paddle_indx]  = Top_ADC_Pulse_Int[i];
        Top_ADC_Amp_Max[paddle_indx]  = Top_ADC_Pulse_Amp[i];
        Top_ADC_Time_Max[paddle_indx] = Top_ADC_Pulse_Time[i];
      }
    }

    for (int i = 0; i < NUM_BARS; i++) {
      if (Btm_ADC_Time_Max[i] != 0) {
        n_ev_used_btm[i]++;
        Btm_ADC_Time_Avg[i] += Btm_ADC_Time_Max[i];
      }
      if (Top_ADC_Time_Max[i] != 0) {
        n_ev_used_top[i]++;
        Top_ADC_Time_Avg[i] += Top_ADC_Time_Max[i];
      }
    }
  }
  for (int i = 0; i < NUM_BARS; i++) {
    Btm_ADC_Time_Avg[i] /= n_ev_used_btm[i];
    Top_ADC_Time_Avg[i] /= n_ev_used_top[i];
  }

  // Loop over the TTree entries
  for (int i_evt = 0; i_evt < tree->GetEntries(); i_evt++) {
    tree->GetEntry(i_evt);

    // Reset arrays
    for (int i = 0; i < NUM_BARS; i++) {
      Btm_ADC_Int_Max[i]  = 0;
      Top_ADC_Int_Max[i]  = 0;
      Btm_ADC_Amp_Max[i]  = 0;
      Top_ADC_Amp_Max[i]  = 0;
      Btm_ADC_Time_Max[i] = 0;
      Top_ADC_Time_Max[i] = 0;
      top_mult[i]         = 0;
      btm_mult[i]         = 0;
    }

    for (int i = 0; i < n_data_btm; i++) {
      // Apply unit conversions
      // I think they're actually already applied in the replay.
      // Btm_ADC_Ped[i]        = Btm_ADC_Ped[i] * adcChanTomV;
      // Btm_ADC_Pulse_Int[i]  = Btm_ADC_Pulse_Int[i] * adcChanTopC;
      // Btm_ADC_Pulse_Amp[i]  = Btm_ADC_Pulse_Amp[i] * adcChanTomV;
      // Btm_ADC_Pulse_Time[i] = Btm_ADC_Pulse_Time[i] * adcTimeRes;

      // if (i_evt == 0) {
      //   cout << "Pulse Int: " << Btm_ADC_Pulse_Int[i] << " Pulse Amp: " << Btm_ADC_Pulse_Amp[i]
      //        << " Pulse Time: " << Btm_ADC_Pulse_Time[i] << endl;
      // }
      // Fill histograms
      int paddle_indx = Btm_ADC_Counter[i] - 1;
      if (Btm_ADC_Pulse_Time[i] < TIME_WINDOW_MIN || Btm_ADC_Pulse_Time[i] > TIME_WINDOW_MAX)
        continue;
      btm_mult[paddle_indx]++;
      h1_Btm_Ped[paddle_indx]->Fill(Btm_ADC_Ped[i]);
      h1_Btm_Pulse_Int[paddle_indx]->Fill(Btm_ADC_Pulse_Int[i]);
      h1_Btm_Pulse_Amp[paddle_indx]->Fill(Btm_ADC_Pulse_Amp[i]);
      h1_Btm_Pulse_Time[paddle_indx]->Fill(Btm_ADC_Pulse_Time[i]);
      h2_Btm_Pulse_Int_Amp[paddle_indx]->Fill(Btm_ADC_Pulse_Int[i], Btm_ADC_Pulse_Amp[i]);
      if (Btm_ADC_Pulse_Int[i] > Btm_ADC_Int_Max[paddle_indx]) {
        Btm_ADC_Int_Max[paddle_indx]  = Btm_ADC_Pulse_Int[i];
        Btm_ADC_Amp_Max[paddle_indx]  = Btm_ADC_Pulse_Amp[i];
        Btm_ADC_Time_Max[paddle_indx] = Btm_ADC_Pulse_Time[i];
      }
    }
    for (int i = 0; i < n_data_top; i++) {
      // Apply unit conversions
      // I think they're actually already applied in the replay.
      // Top_ADC_Ped[i]        = Top_ADC_Ped[i] * adcChanTomV;
      // Top_ADC_Pulse_Int[i]  = Top_ADC_Pulse_Int[i] * adcChanTopC;
      // Top_ADC_Pulse_Amp[i]  = Top_ADC_Pulse_Amp[i] * adcChanTomV;
      // Top_ADC_Pulse_Time[i] = Top_ADC_Pulse_Time[i] * adcTimeRes;
      // Fill histograms
      int paddle_indx = Top_ADC_Counter[i] - 1;
      if (Top_ADC_Pulse_Time[i] < TIME_WINDOW_MIN || Top_ADC_Pulse_Time[i] > TIME_WINDOW_MAX)
        continue;
      top_mult[paddle_indx]++;
      h1_Top_Ped[paddle_indx]->Fill(Top_ADC_Ped[i]);
      h1_Top_Pulse_Int[paddle_indx]->Fill(Top_ADC_Pulse_Int[i]);
      h1_Top_Pulse_Amp[paddle_indx]->Fill(Top_ADC_Pulse_Amp[i]);
      h1_Top_Pulse_Time[paddle_indx]->Fill(Top_ADC_Pulse_Time[i]);
      h2_Top_Pulse_Int_Amp[paddle_indx]->Fill(Top_ADC_Pulse_Int[i], Top_ADC_Pulse_Amp[i]);
      if (Top_ADC_Pulse_Int[i] > Top_ADC_Int_Max[paddle_indx]) {
        Top_ADC_Int_Max[paddle_indx]  = Top_ADC_Pulse_Int[i];
        Top_ADC_Amp_Max[paddle_indx]  = Top_ADC_Pulse_Amp[i];
        Top_ADC_Time_Max[paddle_indx] = Top_ADC_Pulse_Time[i];
      }
    }

    // Count number of fired PMTs per event
    int num_fired_pmts = 0;
    for (int i = 0; i < NUM_BARS; i++) {
      h1_Top_Mult[i]->Fill(top_mult[i]);
      h1_Btm_Mult[i]->Fill(btm_mult[i]);
      if (top_mult[i] > 0)
        num_fired_pmts++;
      if (btm_mult[i] > 0)
        num_fired_pmts++;

      if (Top_ADC_Int_Max[i] == 0 || Btm_ADC_Int_Max[i] == 0)
        continue;
      h2_Max_Adc_Time_Diff->Fill(i, Top_ADC_Time_Max[i] - Btm_ADC_Time_Max[i]);
      h2_Max_Adc_Time_Diff_Normalized->Fill(i, (Top_ADC_Time_Max[i] - Btm_ADC_Time_Max[i]) -
                                                   (Top_ADC_Time_Avg[i] - Btm_ADC_Time_Avg[i]));
      h1_Max_Adc_Time_Diff[i]->Fill(Top_ADC_Time_Max[i] - Btm_ADC_Time_Max[i]);
    }
    h1_Tot_PMTsFired->Fill(num_fired_pmts);

    if (num_fired_pmts == 2 * NUM_BARS || (num_fired_pmts == 2 * NUM_BARS - 1) && run_number == 67) { //run 67 had a faulty PMT
      for (int i = 0; i < NUM_BARS; i++) {
        h1_Btm_Evt_Int[i]->Fill(Btm_ADC_Int_Max[i]);
        h1_Top_Evt_Int[i]->Fill(Top_ADC_Int_Max[i]);
        h1_Btm_Evt_Amp[i]->Fill(Btm_ADC_Amp_Max[i]);
        h1_Top_Evt_Amp[i]->Fill(Top_ADC_Amp_Max[i]);
        // h2_Btm_Evt_Int_Amp[i]->Fill(Btm_ADC_Int_Max[i], Btm_ADC_Pulse_Amp[i]);
        // h2_Top_Evt_Int_Amp[i]->Fill(Top_ADC_Int_Max[i], Top_ADC_Pulse_Amp[i]);
      }
    }

    for (int i = 0; i < n_data_btm; i++) {
      int btm_paddle_indx = Btm_ADC_Counter[i] - 1;
      for (int j = 0; j < n_data_top; j++) {
        int top_paddle_indx = Top_ADC_Counter[j] - 1;
        if (btm_paddle_indx != top_paddle_indx)
          continue;
        if (Btm_ADC_Pulse_Time[i] < TIME_WINDOW_MIN || Btm_ADC_Pulse_Time[i] > TIME_WINDOW_MAX ||
            Top_ADC_Pulse_Time[j] < TIME_WINDOW_MIN || Top_ADC_Pulse_Time[j] > TIME_WINDOW_MAX)
          continue;
        h2_Pulse_Int[btm_paddle_indx]->Fill(Btm_ADC_Pulse_Int[i], Top_ADC_Pulse_Int[j]);
        h2_Pulse_Amp[btm_paddle_indx]->Fill(Btm_ADC_Pulse_Amp[i], Top_ADC_Pulse_Amp[j]);
        h2_Pulse_Time[btm_paddle_indx]->Fill(Btm_ADC_Pulse_Time[i], Top_ADC_Pulse_Time[j]);
      }
    }

    if (num_fired_pmts >= MIN_NUM_PMTS_FIRED && i_evt < 10) {
      TH2D h2_pulse_times(Form("Pulse_Times_%d", i_evt), "Pulse_Times; Pulse Time; Bar Number", TIME_N_BINS, TIME_MIN,
                          TIME_MAX, NUM_BARS, 0, NUM_BARS);
      for (int i = 0; i < n_data_btm; i++) {
        int btm_paddle_indx = Btm_ADC_Counter[i] - 1;
        h2_pulse_times.Fill(Btm_ADC_Pulse_Time[i], btm_paddle_indx);
      }
      for (int i = 0; i < n_data_top; i++) {
        int top_paddle_indx = Top_ADC_Counter[i] - 1;
        h2_pulse_times.Fill(Top_ADC_Pulse_Time[i], top_paddle_indx, 2);
      }
      indivEventDir->cd();
      h2_pulse_times.Write();
    }

    // if (num_fired_pmts == 16 && !plottedWaveform) {
    if (i_evt == 1) {
      plottedWaveform = true;
      indivEventWaveform->cd();
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

    if (num_fired_pmts <= MAX_NOISY_PMTS) {
      for (int i = 0; i < NUM_BARS; i++) {
        if (top_mult[i] != 0)
          h1_Top_Noisy_PMTs->Fill(i);
        if (btm_mult[i] != 0)
          h1_Btm_Noisy_PMTs->Fill(i);
      }
    }
  }

  // Write histograms to the corresponding directories
  for (int i = 0; i < NUM_BARS; i++) {
    pedDir->cd();
    h1_Btm_Ped[i]->Write();

    pulseIntDir->cd();
    h1_Btm_Pulse_Int[i]->Write();

    eventIntDir->cd();
    h1_Btm_Evt_Int[i]->Write();
    h1_Top_Evt_Int[i]->Write();

    pulseAmpDir->cd();
    h1_Btm_Pulse_Amp[i]->Write();

    eventAmpDir->cd();
    h1_Btm_Evt_Amp[i]->Write();
    h1_Top_Evt_Amp[i]->Write();

    pulseTimeDir->cd();
    h1_Btm_Pulse_Time[i]->Write();

    pedDir->cd();
    h1_Top_Ped[i]->Write();

    pulseIntDir->cd();
    h1_Top_Pulse_Int[i]->Write();

    pulseAmpDir->cd();
    h1_Top_Pulse_Amp[i]->Write();

    pulseTimeDir->cd();
    h1_Top_Pulse_Time[i]->Write();

    multDir->cd();
    h1_Top_Mult[i]->Write();
    h1_Btm_Mult[i]->Write();

    topBtmCompDir->cd();
    h2_Pulse_Int[i]->Write();
    h2_Pulse_Amp[i]->Write();
    h2_Pulse_Time[i]->Write();

    pulseTimeDiffDir->cd();
    h1_Max_Adc_Time_Diff[i]->Write();

    intAmpCompDir->cd();
    h2_Btm_Pulse_Int_Amp[i]->Write();
    h2_Top_Pulse_Int_Amp[i]->Write();
    // h2_Btm_Evt_Int_Amp[i]->Write();
    // h2_Top_Evt_Int_Amp[i]->Write();
  }
  multDir->cd();
  h1_Tot_PMTsFired->Write();

  pulseTimeDiffDir->cd();
  h2_Max_Adc_Time_Diff->Write();
  h2_Max_Adc_Time_Diff_Normalized->Write();

  noisyPMTDir->cd();
  h1_Btm_Noisy_PMTs->Write();
  h1_Top_Noisy_PMTs->Write();

  // Close the output file
  outputFile->Close();
  // Close the input file
  inputFile->Close();
}