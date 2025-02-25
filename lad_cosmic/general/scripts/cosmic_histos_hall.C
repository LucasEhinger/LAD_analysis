// Plotting scripts used for cosmic runs in Hall C.

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <iostream>

using namespace ROOT;
using namespace std;

const int NUM_BARS         = 11;
const int N_DATA_MAX       = 100;
// const int NUM_PLANES       = 5;
// const string plane_names[] = {"000", "001", "100", "101", "200"};
const int NUM_PLANES       = 1;
const string plane_names[] = {"200"};

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

static const int T_DIFF_MIN = -50;
static const int T_DIFF_MAX = 50;

static const Double_t adcDynamicRange = 1000.0;                     // Units of mV
static const Double_t nAdcChan        = 4096.0;                     // Units of ADC channels
static const Double_t adcChanTomV     = adcDynamicRange / nAdcChan; // Units of mV/ADC Chan
static const Double_t adcImpedence    = 50.0;                       // FADC input impedence in units of Ohms
static const Double_t adcTimeSample   = 4.0;                        // Length of FADC time sample in units of ns
static const Double_t adcTimeRes      = 0.0625;                     // FADC time resolution in units of ns
static const Double_t adcChanTopC     = (adcDynamicRange / 1000 / nAdcChan) * (adcTimeSample / adcImpedence);
// (1000 mV / 4096 adc channels) * (4 ns time sample / 50 ohms input resistance) = ~0.020 pc/channel

void cosmic_histos_hall(int run_number) {
  // int run_number = 70;
  // int run_number = 22;

  // Open the input file
  string input_string =
      "/volatile/hallc/c-lad/ehingerl/ROOTfiles/COSMICS/LAD_cosmic_hall_" + to_string(run_number) + "_-1.root";
  TFile *inputFile = TFile::Open(input_string.c_str(), "READ");

  // Get the TTree from the input file
  TTree *tree = dynamic_cast<TTree *>(inputFile->Get("T"));

  // Create a new output file
  string output_string = "../histos/cosmic_histos_" + to_string(run_number) + "_output.root";
  TFile *outputFile    = TFile::Open(output_string.c_str(), "RECREATE");
  // Create directories for each type of histogram
  TDirectory *pedDir             = outputFile->mkdir("Ped");
  TDirectory *pulseIntDir        = outputFile->mkdir("Pulse_Int");
  TDirectory *eventIntDir        = pulseIntDir->mkdir("Event_Int");
  TDirectory *pulseAmpDir        = outputFile->mkdir("Pulse_Amp");
  TDirectory *eventAmpDir        = pulseAmpDir->mkdir("Event_Amp");
  TDirectory *pulseTimeDir       = outputFile->mkdir("Pulse_Time");
  TDirectory *eventTimeDir       = pulseTimeDir->mkdir("Event_Time");
  TDirectory *topBtmCompDir      = outputFile->mkdir("2D_Top_Btm_Comparison");
  TDirectory *multDir            = outputFile->mkdir("Mult");
  TDirectory *indivEventDir      = outputFile->mkdir("Indiv_Event");
  TDirectory *indivEventWaveform = outputFile->mkdir("Indiv_Event_Waveform");
  TDirectory *pulseTimeDiffDir   = outputFile->mkdir("Pulse_Time_Diff");
  TDirectory *intAmpCompDir      = outputFile->mkdir("2D_Int_Amp_Comparison");

  // Create histograms
  TH1F *h1_Btm_Ped[NUM_PLANES][NUM_BARS];
  TH1F *h1_Btm_Pulse_Int[NUM_PLANES][NUM_BARS];
  TH1F *h1_Btm_Pulse_Amp[NUM_PLANES][NUM_BARS];
  TH1F *h1_Btm_Pulse_Time[NUM_PLANES][NUM_BARS];

  TH1F *h1_Top_Ped[NUM_PLANES][NUM_BARS];
  TH1F *h1_Top_Pulse_Int[NUM_PLANES][NUM_BARS];
  TH1F *h1_Top_Pulse_Amp[NUM_PLANES][NUM_BARS];
  TH1F *h1_Top_Pulse_Time[NUM_PLANES][NUM_BARS];

  TH1F *h1_Btm_Evt_Int[NUM_PLANES][NUM_BARS];
  TH1F *h1_Top_Evt_Int[NUM_PLANES][NUM_BARS];
  TH1F *h1_Btm_Evt_Amp[NUM_PLANES][NUM_BARS];
  TH1F *h1_Top_Evt_Amp[NUM_PLANES][NUM_BARS];
  TH1F *h1_Btm_Evt_Time[NUM_PLANES][NUM_BARS];
  TH1F *h1_Top_Evt_Time[NUM_PLANES][NUM_BARS];

  TH1F *h1_Top_Mult[NUM_PLANES][NUM_BARS];
  TH1F *h1_Btm_Mult[NUM_PLANES][NUM_BARS];
  TH1F *h1_Tot_PMTsFired[NUM_PLANES];

  TH2D *h2_Pulse_Int[NUM_PLANES][NUM_BARS];
  TH2D *h2_Pulse_Amp[NUM_PLANES][NUM_BARS];
  TH2D *h2_Pulse_Time[NUM_PLANES][NUM_BARS];
  TH2D *h2_Btm_Pulse_Int_Amp[NUM_PLANES][NUM_BARS];
  TH2D *h2_Top_Pulse_Int_Amp[NUM_PLANES][NUM_BARS];

  TH2D *h2_Max_Adc_Time_Diff[NUM_PLANES];
  TH2D *h2_Max_Adc_Time_Diff_Normalized[NUM_PLANES];
  TH1F *h1_Max_Adc_Time_Diff[NUM_PLANES][NUM_BARS];


  for (int plane_idx = 0; plane_idx < NUM_PLANES; plane_idx++) {
    for (int i = 0; i < NUM_BARS; i++) {
      h1_Btm_Ped[plane_idx][i] =
          new TH1F(Form("Btm_Ped_%s_%d", plane_names[plane_idx].c_str(), i), 
                   Form("Btm_Ped_%s_%d; Pulse Ped; Counts", plane_names[plane_idx].c_str(), i),
                   PED_N_BINS, PED_MIN, PED_MAX);
      h1_Btm_Pulse_Int[plane_idx][i] = new TH1F(Form("Btm_Pulse_Int_%s_%d", plane_names[plane_idx].c_str(), i),
                                                Form("Btm_Pulse_Int_%s_%d; Pulse Int [pC]; Counts", plane_names[plane_idx].c_str(), i),
                                                INTEGRAL_N_BINS, INTEGRAL_MIN, INTEGRAL_MAX);
      h1_Btm_Pulse_Amp[plane_idx][i] =
          new TH1F(Form("Btm_Pulse_Amp_%s_%d", plane_names[plane_idx].c_str(), i),
                   Form("Btm_Pulse_Amp_%s_%d; Pulse Amp [mV]; Counts", plane_names[plane_idx].c_str(), i), AMP_N_BINS, AMP_MIN, AMP_MAX);
      h1_Btm_Pulse_Time[plane_idx][i] = new TH1F(Form("Btm_Pulse_Time_%s_%d", plane_names[plane_idx].c_str(), i),
                                                 Form("Btm_Pulse_Time_%s_%d; Pulse Time [ns]; Counts", plane_names[plane_idx].c_str(), i),
                                                 TIME_N_BINS, TIME_MIN, TIME_MAX);

      h1_Top_Ped[plane_idx][i] =
          new TH1F(Form("Top_Ped_%s_%d", plane_names[plane_idx].c_str(), i), 
                   Form("Top_Ped_%s_%d; Pulse Ped; Counts", plane_names[plane_idx].c_str(), i),
                   PED_N_BINS, PED_MIN, PED_MAX);
      h1_Top_Pulse_Int[plane_idx][i] = new TH1F(Form("Top_Pulse_Int_%s_%d", plane_names[plane_idx].c_str(), i),
                                                Form("Top_Pulse_Int_%s_%d; Pulse Int [pC]; Counts", plane_names[plane_idx].c_str(), i),
                                                INTEGRAL_N_BINS, INTEGRAL_MIN, INTEGRAL_MAX);
      h1_Top_Pulse_Amp[plane_idx][i] =
          new TH1F(Form("Top_Pulse_Amp_%s_%d", plane_names[plane_idx].c_str(), i),
                   Form("Top_Pulse_Amp_%s_%d; Pulse Amp [mV]; Counts", plane_names[plane_idx].c_str(), i), AMP_N_BINS, AMP_MIN, AMP_MAX);
      h1_Top_Pulse_Time[plane_idx][i] = new TH1F(Form("Top_Pulse_Time_%s_%d", plane_names[plane_idx].c_str(), i),
                                                 Form("Top_Pulse_Time_%s_%d; Pulse Time [ns]; Counts", plane_names[plane_idx].c_str(), i),
                                                 TIME_N_BINS, TIME_MIN, TIME_MAX);

      h1_Btm_Evt_Int[plane_idx][i] = new TH1F(Form("Btm_Evt_Int_%s_%d", plane_names[plane_idx].c_str(), i),
                                              Form("Btm_Evt_Int_%s_%d; Pulse Int [pC]; Counts", plane_names[plane_idx].c_str(), i),
                                              INTEGRAL_N_BINS, INTEGRAL_MIN, INTEGRAL_MAX);
      h1_Top_Evt_Int[plane_idx][i] = new TH1F(Form("Top_Evt_Int_%s_%d", plane_names[plane_idx].c_str(), i),
                                              Form("Top_Evt_Int_%s_%d; Pulse Int [pC]; Counts", plane_names[plane_idx].c_str(), i),
                                              INTEGRAL_N_BINS, INTEGRAL_MIN, INTEGRAL_MAX);
      h1_Btm_Evt_Amp[plane_idx][i] =
          new TH1F(Form("Btm_Evt_Amp_%s_%d", plane_names[plane_idx].c_str(), i),
                   Form("Btm_Evt_Amp_%s_%d; Pulse Amp [mV]; Counts", plane_names[plane_idx].c_str(), i), AMP_N_BINS, AMP_MIN, AMP_MAX);
      h1_Top_Evt_Amp[plane_idx][i] =
          new TH1F(Form("Top_Evt_Amp_%s_%d", plane_names[plane_idx].c_str(), i),
                   Form("Top_Evt_Amp_%s_%d; Pulse Amp [mV]; Counts", plane_names[plane_idx].c_str(), i), AMP_N_BINS, AMP_MIN, AMP_MAX);

      h1_Btm_Evt_Time[plane_idx][i] =
          new TH1F(Form("Btm_Evt_Time_%s_%d", plane_names[plane_idx].c_str(), i),
                   Form("Btm_Evt_Time_%s_%d; Pulse Time [ns]; Counts", plane_names[plane_idx].c_str(), i), TIME_N_BINS, TIME_MIN, TIME_MAX);
      h1_Top_Evt_Time[plane_idx][i] =
          new TH1F(Form("Top_Evt_Time_%s_%d", plane_names[plane_idx].c_str(), i),
                   Form("Top_Evt_Time_%s_%d; Pulse Time [ns]; Counts", plane_names[plane_idx].c_str(), i), TIME_N_BINS, TIME_MIN, TIME_MAX);

      h1_Top_Mult[plane_idx][i] = new TH1F(Form("Top_Mult_%s_%d", plane_names[plane_idx].c_str(), i),
                                           Form("Top_Mult_%s_%d; Paddle Multiplicity; Counts", plane_names[plane_idx].c_str(), i),
                                           NUM_BARS + 1, -0.5, NUM_BARS + 0.5);
      h1_Btm_Mult[plane_idx][i] = new TH1F(Form("Btm_Mult_%s_%d", plane_names[plane_idx].c_str(), i),
                                           Form("Btm_Mult_%s_%d; Paddle Multiplicity; Counts", plane_names[plane_idx].c_str(), i),
                                           NUM_BARS + 1, -0.5, NUM_BARS + 0.5);

      h2_Pulse_Int[plane_idx][i] =
          new TH2D(Form("Pulse_Int_%s_%d", plane_names[plane_idx].c_str(), i), 
                   Form("Pulse_Int_%s_%d; Btm [pC]; Top [pC]", plane_names[plane_idx].c_str(), i),
                   INTEGRAL_N_BINS, INTEGRAL_MIN, INTEGRAL_MAX, INTEGRAL_N_BINS, INTEGRAL_MIN, INTEGRAL_MAX);
      h2_Pulse_Amp[plane_idx][i] =
          new TH2D(Form("Pulse_Amp_%s_%d", plane_names[plane_idx].c_str(), i), 
                   Form("Pulse_Amp_%s_%d; Btm [pC]; Top [pC]", plane_names[plane_idx].c_str(), i),
                   AMP_N_BINS, AMP_MIN, AMP_MAX, AMP_N_BINS, AMP_MIN, AMP_MAX);
      h2_Pulse_Time[plane_idx][i] =
          new TH2D(Form("Pulse_Time_%s_%d", plane_names[plane_idx].c_str(), i), 
                   Form("Pulse_Time_%s_%d; Btm [ns]; Top [ns]", plane_names[plane_idx].c_str(), i),
                   TIME_N_BINS, TIME_MIN, TIME_MAX, TIME_N_BINS, TIME_MIN, TIME_MAX);

      h2_Btm_Pulse_Int_Amp[plane_idx][i] =
          new TH2D(Form("Btm_Pulse_Int_Amp_%s_%d", plane_names[plane_idx].c_str(), i),
                   Form("Btm_Pulse_Int_Amp_%s_%d; Int  [pC]; Amp [mV]", plane_names[plane_idx].c_str(), i), INTEGRAL_N_BINS, INTEGRAL_MIN,
                   INTEGRAL_MAX, AMP_N_BINS, AMP_MIN, AMP_MAX);
      h2_Top_Pulse_Int_Amp[plane_idx][i] =
          new TH2D(Form("Top_Pulse_Int_Amp_%s_%d", plane_names[plane_idx].c_str(), i),
                   Form("Top_Pulse_Int_Amp_%s_%d; Int [pC]; Amp [mV]", plane_names[plane_idx].c_str(), i), INTEGRAL_N_BINS, INTEGRAL_MIN,
                   INTEGRAL_MAX, AMP_N_BINS, AMP_MIN, AMP_MAX);

      h1_Max_Adc_Time_Diff[plane_idx][i] =
          new TH1F(Form("Max_Adc_Time_Diff_%s_%d", plane_names[plane_idx].c_str(), i),
                   Form("Adc_Time_Diff_%s_%d; t_top - t_btm [ns]; Counts", plane_names[plane_idx].c_str(), i), TIME_DIFF_N_BINS,
                   TIME_DIFF_MIN, TIME_DIFF_MAX);
    }
    h1_Tot_PMTsFired[plane_idx] =
        new TH1F(Form("Tot_PMTs_Fired_%s", plane_names[plane_idx].c_str()), 
                 Form("Tot_PMTs_Fired_%s; Num PMTs Fired; Counts", plane_names[plane_idx].c_str()),
                 NUM_BARS * 2 + 1, -0.5, NUM_BARS * 2 + 0.5);
    h2_Max_Adc_Time_Diff[plane_idx] = new TH2D(
        Form("Max_Adc_Time_Diff_%s", plane_names[plane_idx].c_str()), 
        Form("Adc_Time_Diff_%s; Bar Number; t_top - t_btm [ns]", plane_names[plane_idx].c_str()),
        NUM_BARS, -0.5, NUM_BARS - 0.5, TIME_DIFF_N_BINS, TIME_DIFF_MIN, TIME_DIFF_MAX);
  }

  // Create variables to hold the data
  Double_t Btm_ADC_Counter[NUM_PLANES][N_DATA_MAX]    = {0};
  Double_t Btm_ADC_Ped[NUM_PLANES][N_DATA_MAX]        = {0};
  Double_t Btm_ADC_Pulse_Int[NUM_PLANES][N_DATA_MAX]  = {0};
  Double_t Btm_ADC_Pulse_Amp[NUM_PLANES][N_DATA_MAX]  = {0};
  Double_t Btm_ADC_Pulse_Time[NUM_PLANES][N_DATA_MAX] = {0};

  Double_t Top_ADC_Counter[NUM_PLANES][N_DATA_MAX]    = {0};
  Double_t Top_ADC_Ped[NUM_PLANES][N_DATA_MAX]        = {0};
  Double_t Top_ADC_Pulse_Int[NUM_PLANES][N_DATA_MAX]  = {0};
  Double_t Top_ADC_Pulse_Amp[NUM_PLANES][N_DATA_MAX]  = {0};
  Double_t Top_ADC_Pulse_Time[NUM_PLANES][N_DATA_MAX] = {0};

  Double_t Btm_ADC_Waveform[NUM_PLANES][MAX_WAVEFORM_LENGTH] = {0};
  Double_t Top_ADC_Waveform[NUM_PLANES][MAX_WAVEFORM_LENGTH] = {0};

  Double_t Btm_ADC_Int_Max[NUM_PLANES][NUM_BARS]  = {0};
  Double_t Top_ADC_Int_Max[NUM_PLANES][NUM_BARS]  = {0};
  Double_t Btm_ADC_Amp_Max[NUM_PLANES][NUM_BARS]  = {0};
  Double_t Top_ADC_Amp_Max[NUM_PLANES][NUM_BARS]  = {0};
  Double_t Btm_ADC_Time_Max[NUM_PLANES][NUM_BARS] = {0};
  Double_t Top_ADC_Time_Max[NUM_PLANES][NUM_BARS] = {0};
  Double_t Btm_ADC_Ped_Max[NUM_PLANES][NUM_BARS]  = {0};
  Double_t Top_ADC_Ped_Max[NUM_PLANES][NUM_BARS]  = {0};

  int top_mult[NUM_PLANES][NUM_BARS] = {0};
  int btm_mult[NUM_PLANES][NUM_BARS] = {0};
  double t_diff[NUM_PLANES][NUM_BARS] = {0};

  int n_data_top[NUM_PLANES];
  int n_data_btm[NUM_PLANES];
  bool plottedWaveform[NUM_PLANES] = {false};

  for (int plane_idx = 0; plane_idx < NUM_PLANES; plane_idx++) {
    tree->SetBranchAddress(Form("L.hod.%s.BtmAdcCounter", plane_names[plane_idx].c_str()), Btm_ADC_Counter[plane_idx]);
    tree->SetBranchAddress(Form("L.hod.%s.BtmAdcPed", plane_names[plane_idx].c_str()), Btm_ADC_Ped[plane_idx]);
    tree->SetBranchAddress(Form("L.hod.%s.BtmAdcPulseInt", plane_names[plane_idx].c_str()),
                           Btm_ADC_Pulse_Int[plane_idx]);
    tree->SetBranchAddress(Form("L.hod.%s.BtmAdcPulseAmp", plane_names[plane_idx].c_str()),
                           Btm_ADC_Pulse_Amp[plane_idx]);
    tree->SetBranchAddress(Form("L.hod.%s.BtmAdcPulseTime", plane_names[plane_idx].c_str()),
                           Btm_ADC_Pulse_Time[plane_idx]);

    tree->SetBranchAddress(Form("L.hod.%s.TopAdcCounter", plane_names[plane_idx].c_str()), Top_ADC_Counter[plane_idx]);
    tree->SetBranchAddress(Form("L.hod.%s.TopAdcPed", plane_names[plane_idx].c_str()), Top_ADC_Ped[plane_idx]);
    tree->SetBranchAddress(Form("L.hod.%s.TopAdcPulseInt", plane_names[plane_idx].c_str()),
                           Top_ADC_Pulse_Int[plane_idx]);
    tree->SetBranchAddress(Form("L.hod.%s.TopAdcPulseAmp", plane_names[plane_idx].c_str()),
                           Top_ADC_Pulse_Amp[plane_idx]);
    tree->SetBranchAddress(Form("L.hod.%s.TopAdcPulseTime", plane_names[plane_idx].c_str()),
                           Top_ADC_Pulse_Time[plane_idx]);

    tree->SetBranchAddress(Form("Ndata.L.hod.%s.BtmAdcCounter", plane_names[plane_idx].c_str()),
                           &n_data_btm[plane_idx]);
    tree->SetBranchAddress(Form("Ndata.L.hod.%s.TopAdcCounter", plane_names[plane_idx].c_str()),
                           &n_data_top[plane_idx]);

    tree->SetBranchAddress(Form("L.hod.%s.adcBtmSampWaveform", plane_names[plane_idx].c_str()),
                           Btm_ADC_Waveform[plane_idx]);
    tree->SetBranchAddress(Form("L.hod.%s.adcTopSampWaveform", plane_names[plane_idx].c_str()),
                           Top_ADC_Waveform[plane_idx]);
  }

  // Loop over the TTree entries
  for (int i_evt = 0; i_evt < tree->GetEntries(); i_evt++) {
    tree->GetEntry(i_evt);
    for (int plane_idx = 0; plane_idx < NUM_PLANES; plane_idx++) {
      // Reset arrays
      for (int i = 0; i < NUM_BARS; i++) {
        Btm_ADC_Int_Max[plane_idx][i]  = 0;
        Top_ADC_Int_Max[plane_idx][i]  = 0;
        Btm_ADC_Amp_Max[plane_idx][i]  = 0;
        Top_ADC_Amp_Max[plane_idx][i]  = 0;
        Btm_ADC_Time_Max[plane_idx][i] = 0;
        Top_ADC_Time_Max[plane_idx][i] = 0;
        top_mult[plane_idx][i]         = 0;
        btm_mult[plane_idx][i]         = 0;
      }

      //Populate arrays with only 1 hit per paddle
      for (int i = 0; i < n_data_btm[plane_idx]; i++) {
        // Fill histograms
        int paddle_indx = Btm_ADC_Counter[plane_idx][i] - 1;
        if (Btm_ADC_Pulse_Time[plane_idx][i] < TIME_WINDOW_MIN || Btm_ADC_Pulse_Time[plane_idx][i] > TIME_WINDOW_MAX)
          continue;
        btm_mult[plane_idx][paddle_indx]++;
        if (Btm_ADC_Pulse_Int[plane_idx][i] > Btm_ADC_Int_Max[plane_idx][paddle_indx]) {
          Btm_ADC_Int_Max[plane_idx][paddle_indx]  = Btm_ADC_Pulse_Int[plane_idx][i];
          Btm_ADC_Amp_Max[plane_idx][paddle_indx]  = Btm_ADC_Pulse_Amp[plane_idx][i];
          Btm_ADC_Time_Max[plane_idx][paddle_indx] = Btm_ADC_Pulse_Time[plane_idx][i];
        }
      }
      for (int i = 0; i < n_data_top[plane_idx]; i++) {
        int paddle_indx = Top_ADC_Counter[plane_idx][i] - 1;
        if (Top_ADC_Pulse_Time[plane_idx][i] < TIME_WINDOW_MIN || Top_ADC_Pulse_Time[plane_idx][i] > TIME_WINDOW_MAX)
          continue;
        top_mult[plane_idx][paddle_indx]++;
        if (Top_ADC_Pulse_Int[plane_idx][i] > Top_ADC_Int_Max[plane_idx][paddle_indx]) {
          Top_ADC_Int_Max[plane_idx][paddle_indx]  = Top_ADC_Pulse_Int[plane_idx][i];
          Top_ADC_Amp_Max[plane_idx][paddle_indx]  = Top_ADC_Pulse_Amp[plane_idx][i];
          Top_ADC_Time_Max[plane_idx][paddle_indx] = Top_ADC_Pulse_Time[plane_idx][i];
        }
      }

      // Fill histograms (all hits)
      for (int i = 0; i < n_data_btm[plane_idx]; i++) {
        int paddle_indx = Btm_ADC_Counter[plane_idx][i] - 1;
        if (btm_mult[plane_idx][paddle_indx] == 0 || top_mult[plane_idx][paddle_indx] == 0) {
          continue;
        }

        h1_Btm_Ped[plane_idx][paddle_indx]->Fill(Btm_ADC_Ped[plane_idx][i]);
        h1_Btm_Pulse_Int[plane_idx][paddle_indx]->Fill(Btm_ADC_Pulse_Int[plane_idx][i]);
        h1_Btm_Pulse_Amp[plane_idx][paddle_indx]->Fill(Btm_ADC_Pulse_Amp[plane_idx][i]);
        h1_Btm_Pulse_Time[plane_idx][paddle_indx]->Fill(Btm_ADC_Pulse_Time[plane_idx][i]);
        h2_Btm_Pulse_Int_Amp[plane_idx][paddle_indx]->Fill(Btm_ADC_Pulse_Int[plane_idx][i],
                                                           Btm_ADC_Pulse_Amp[plane_idx][i]);
      }
      for (int i = 0; i < n_data_top[plane_idx]; i++) {
        int paddle_indx = Top_ADC_Counter[plane_idx][i] - 1;
        if (btm_mult[plane_idx][paddle_indx] == 0 || top_mult[plane_idx][paddle_indx] == 0) {
          continue;
        }

        h1_Top_Ped[plane_idx][paddle_indx]->Fill(Top_ADC_Ped[plane_idx][i]);
        h1_Top_Pulse_Int[plane_idx][paddle_indx]->Fill(Top_ADC_Pulse_Int[plane_idx][i]);
        h1_Top_Pulse_Amp[plane_idx][paddle_indx]->Fill(Top_ADC_Pulse_Amp[plane_idx][i]);
        h1_Top_Pulse_Time[plane_idx][paddle_indx]->Fill(Top_ADC_Pulse_Time[plane_idx][i]);
        h2_Top_Pulse_Int_Amp[plane_idx][paddle_indx]->Fill(Top_ADC_Pulse_Int[plane_idx][i],
                                                           Top_ADC_Pulse_Amp[plane_idx][i]);
      }
      
      // Fill histograms (all paddles. Max hit per paddle)
      for (int i = 0; i < NUM_BARS; i++) {
        if (Top_ADC_Int_Max[plane_idx][i] > 0 && Btm_ADC_Int_Max[plane_idx][i] > 0) {
          t_diff[plane_idx][i] = Top_ADC_Time_Max[plane_idx][i] - Btm_ADC_Time_Max[plane_idx][i];
          if (t_diff[plane_idx][i] < T_DIFF_MIN || t_diff[plane_idx][i] > T_DIFF_MAX) {
            continue;
          }
          h1_Btm_Evt_Int[plane_idx][i]->Fill(Btm_ADC_Int_Max[plane_idx][i]);
          h1_Top_Evt_Int[plane_idx][i]->Fill(Top_ADC_Int_Max[plane_idx][i]);
          h1_Btm_Evt_Amp[plane_idx][i]->Fill(Btm_ADC_Amp_Max[plane_idx][i]);
          h1_Top_Evt_Amp[plane_idx][i]->Fill(Top_ADC_Amp_Max[plane_idx][i]);
          h1_Btm_Evt_Time[plane_idx][i]->Fill(Btm_ADC_Time_Max[plane_idx][i]);
          h1_Top_Evt_Time[plane_idx][i]->Fill(Top_ADC_Time_Max[plane_idx][i]);
          h1_Max_Adc_Time_Diff[plane_idx][i]->Fill(Top_ADC_Time_Max[plane_idx][i] - Btm_ADC_Time_Max[plane_idx][i]);
          h2_Max_Adc_Time_Diff[plane_idx]->Fill(i, Top_ADC_Time_Max[plane_idx][i] - Btm_ADC_Time_Max[plane_idx][i]);
        }
      }

      // Count number of fired PMTs per event
      int num_fired_pmts = 0;
      for (int i = 0; i < NUM_BARS; i++) {
        h1_Top_Mult[plane_idx][i]->Fill(top_mult[plane_idx][i]);
        h1_Btm_Mult[plane_idx][i]->Fill(btm_mult[plane_idx][i]);
        if (top_mult[plane_idx][i] > 0)
          num_fired_pmts++;
        if (btm_mult[plane_idx][i] > 0)
          num_fired_pmts++;

        if (Top_ADC_Int_Max[plane_idx][i] == 0 || Btm_ADC_Int_Max[plane_idx][i] == 0)
          continue;
        h2_Max_Adc_Time_Diff[plane_idx]->Fill(i, Top_ADC_Time_Max[plane_idx][i] - Btm_ADC_Time_Max[plane_idx][i]);
        h1_Max_Adc_Time_Diff[plane_idx][i]->Fill(Top_ADC_Time_Max[plane_idx][i] - Btm_ADC_Time_Max[plane_idx][i]);
      }
      h1_Tot_PMTsFired[plane_idx]->Fill(num_fired_pmts);

     
      // Fill 2D histograms
      for (int i = 0; i < n_data_btm[plane_idx]; i++) {
        int btm_paddle_indx = Btm_ADC_Counter[plane_idx][i] - 1;
        for (int j = 0; j < n_data_top[plane_idx]; j++) {
          int top_paddle_indx = Top_ADC_Counter[plane_idx][j] - 1;
          if (btm_paddle_indx != top_paddle_indx)
            continue;
          if (Btm_ADC_Pulse_Time[plane_idx][i] < TIME_WINDOW_MIN ||
              Btm_ADC_Pulse_Time[plane_idx][i] > TIME_WINDOW_MAX ||
              Top_ADC_Pulse_Time[plane_idx][j] < TIME_WINDOW_MIN || Top_ADC_Pulse_Time[plane_idx][j] > TIME_WINDOW_MAX)
            continue;
          h2_Pulse_Int[plane_idx][btm_paddle_indx]->Fill(Btm_ADC_Pulse_Int[plane_idx][i],
                                                         Top_ADC_Pulse_Int[plane_idx][j]);
          h2_Pulse_Amp[plane_idx][btm_paddle_indx]->Fill(Btm_ADC_Pulse_Amp[plane_idx][i],
                                                         Top_ADC_Pulse_Amp[plane_idx][j]);
          h2_Pulse_Time[plane_idx][btm_paddle_indx]->Fill(Btm_ADC_Pulse_Time[plane_idx][i],
                                                          Top_ADC_Pulse_Time[plane_idx][j]);
        }
      }

      
      // if (num_fired_pmts == 16 && !plottedWaveform) {
      if (i_evt == 1) {
        plottedWaveform[plane_idx] = true;
        indivEventWaveform->cd();
        int indx = 0;
        while (Btm_ADC_Waveform[plane_idx][indx] != 0) {
          int btm_paddle_indx = (int)Btm_ADC_Waveform[plane_idx][indx++] - 1;
          int n_samples       = (int)Btm_ADC_Waveform[plane_idx][indx++];
          TH1D h1_btm_waveform(Form("Btm_Waveform_%d", btm_paddle_indx), "Btm_Waveform; Time [ns]; ADC Value [mV]",
                               n_samples, 0, n_samples * adcTimeSample);
          for (int i = 0; i < n_samples; i++) {
            h1_btm_waveform.Fill(i * adcTimeSample, Btm_ADC_Waveform[plane_idx][indx + i] * adcChanTomV);
          }
          h1_btm_waveform.Write();
          indx += n_samples;
        }
        indx = 0;
        while (Top_ADC_Waveform[plane_idx][indx] != 0) {
          int top_paddle_indx = (int)Top_ADC_Waveform[plane_idx][indx++] - 1;
          int n_samples       = (int)Top_ADC_Waveform[plane_idx][indx++];
          TH1D h1_top_waveform(Form("Top_Waveform_%d", top_paddle_indx), "Top_Waveform; Time [ns]; ADC Value [mV]",
                               n_samples, 0, n_samples * adcTimeSample);
          for (int i = 0; i < n_samples; i++) {
            h1_top_waveform.Fill(i * adcTimeSample, Top_ADC_Waveform[plane_idx][indx + i] * adcChanTomV);
          }
          h1_top_waveform.Write();
          indx += n_samples;
        }
      }
      if (i_evt % 10000 == 0) {
        cout << "Processed " << i_evt << " events" << endl;
      }
    }
  }

  // Write histograms to the corresponding directories
  for (int plane_idx = 0; plane_idx < NUM_PLANES; plane_idx++) {
    for (int i = 0; i < NUM_BARS; i++) {
      pedDir->cd();
      h1_Btm_Ped[plane_idx][i]->Write();
      h1_Top_Ped[plane_idx][i]->Write();

      pulseIntDir->cd();
      h1_Btm_Pulse_Int[plane_idx][i]->Write();
      h1_Top_Pulse_Int[plane_idx][i]->Write();

      eventIntDir->cd();
      h1_Btm_Evt_Int[plane_idx][i]->Write();
      h1_Top_Evt_Int[plane_idx][i]->Write();

      pulseAmpDir->cd();
      h1_Btm_Pulse_Amp[plane_idx][i]->Write();
      h1_Top_Pulse_Amp[plane_idx][i]->Write();

      eventAmpDir->cd();
      h1_Btm_Evt_Amp[plane_idx][i]->Write();
      h1_Top_Evt_Amp[plane_idx][i]->Write();

      pulseTimeDir->cd();
      h1_Btm_Pulse_Time[plane_idx][i]->Write();
      h1_Top_Pulse_Time[plane_idx][i]->Write();

      eventTimeDir->cd();
      h1_Btm_Evt_Time[plane_idx][i]->Write();
      h1_Top_Evt_Time[plane_idx][i]->Write();

      multDir->cd();
      h1_Top_Mult[plane_idx][i]->Write();
      h1_Btm_Mult[plane_idx][i]->Write();

      topBtmCompDir->cd();
      h2_Pulse_Int[plane_idx][i]->Write();
      h2_Pulse_Amp[plane_idx][i]->Write();
      h2_Pulse_Time[plane_idx][i]->Write();

      pulseTimeDiffDir->cd();
      h1_Max_Adc_Time_Diff[plane_idx][i]->Write();

      intAmpCompDir->cd();
      h2_Btm_Pulse_Int_Amp[plane_idx][i]->Write();
      h2_Top_Pulse_Int_Amp[plane_idx][i]->Write();
    }

    multDir->cd();
    h1_Tot_PMTsFired[plane_idx]->Write();

    pulseTimeDiffDir->cd();
    h2_Max_Adc_Time_Diff[plane_idx]->Write();


  }

  // Close the output file
  outputFile->Close();
  // Close the input file
  inputFile->Close();
}