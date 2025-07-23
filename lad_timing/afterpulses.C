#include <Rtypes.h> // For Double_t
#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TTree.h>
#include <iostream>
#include <string>

const int MAX_DATA   = 500;
const int CACHE_SIZE = 100 * 1024 * 1024; // 100 MB cache size

const int N_PLANES                      = 5;
const int N_PADDLES                     = 11;
const std::string plane_names[N_PLANES] = {"000", "001", "100", "101", "200"};
const int N_SIDES                       = 2;
const std::string side_names[N_SIDES]   = {"Top", "Btm"};

const char spec_prefix = 'P'; // Spectrometer to replay

const double kBIG   = std::numeric_limits<double>::max();
const double TDC2NS = 0.09766; // TDC to ns conversion factor

template <typename T> void add_branch(TTree *tree, const char *branch_name, T *branch_data) {
  // Add a branch to the tree
  tree->SetBranchStatus(branch_name, 1); // Enable the branch
  tree->SetBranchAddress(branch_name, branch_data);
  tree->AddBranchToCache(branch_name, kTRUE);
};

void afterpulses() {

  gROOT->SetBatch(kTRUE);
  TChain *T = new TChain("T");
  T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23471_0_1_-1.root");
  T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23471_10_11_-1.root");
  T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23471_12_13_-1.root");
  T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23471_14_15_-1.root");
  T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23471_16_17_-1.root");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23471_18_19_-1.root");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23471_20_21_-1.root");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23471_22_23_-1.root");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23471_2_3_-1.root");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23471_24_25_-1.root");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23471_4_5_-1.root");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23471_6_7_-1.root");
  // T->Add("/cache/hallc/c-lad/analysis/ehingerl/online_v2/LAD_COIN_23471_8_9_-1.root");
  // TTree *T             = dynamic_cast<TTree *>(file->Get("T"));

  Double_t tdc_time[N_SIDES][N_PLANES][MAX_DATA];
  Double_t tdc_counter[N_SIDES][N_PLANES][MAX_DATA];
  Double_t adc_time[N_SIDES][N_PLANES][MAX_DATA];
  Double_t adc_counter[N_SIDES][N_PLANES][MAX_DATA];
  Int_t nData_tdc[N_SIDES][N_PLANES]              = {0};
  Int_t nData_adc[N_SIDES][N_PLANES]              = {0};
  int nData_tot_adc[N_SIDES][N_PLANES][N_PADDLES] = {0};
  int nData_tot_tdc[N_SIDES][N_PLANES][N_PADDLES] = {0};

  for (int side = 0; side < N_SIDES; ++side) {
    for (int plane = 0; plane < N_PLANES; ++plane) {
      add_branch(T,
                 Form("%c.ladhod.%s.%sTdcTimeRaw", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
                 &tdc_time[side][plane]);
      add_branch(T,
                 Form("%c.ladhod.%s.%sAdcPulseTime", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
                 &adc_time[side][plane]);
      add_branch(T,
                 Form("%c.ladhod.%s.%sTdcCounter", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
                 &tdc_counter[side][plane]);
      add_branch(T,
                 Form("%c.ladhod.%s.%sAdcCounter", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
                 &adc_counter[side][plane]);
      add_branch(
          T, Form("Ndata.%c.ladhod.%s.%sAdcCounter", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
          &nData_adc[side][plane]);
      add_branch(
          T, Form("Ndata.%c.ladhod.%s.%sTdcCounter", spec_prefix, plane_names[plane].c_str(), side_names[side].c_str()),
          &nData_tdc[side][plane]);
    }
  }

  Double_t last_tdc_time[N_SIDES][N_PLANES][N_PADDLES];
  Double_t last_adc_time[N_SIDES][N_PLANES][N_PADDLES];
  std::fill_n(&last_tdc_time[0][0][0], N_SIDES * N_PLANES * N_PADDLES, kBIG);
  std::fill_n(&last_adc_time[0][0][0], N_SIDES * N_PLANES * N_PADDLES, kBIG);

  TH1F *hit_time_diff_adc[N_SIDES][N_PLANES][N_PADDLES];
  TH1F *hit_time_diff_tdc[N_SIDES][N_PLANES][N_PADDLES];
  for (int side = 0; side < N_SIDES; ++side) {
    for (int plane = 0; plane < N_PLANES; ++plane) {
      for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
        hit_time_diff_adc[side][plane][paddle] =
            new TH1F(Form("hit_time_diff_%s_%s_%d", side_names[side].c_str(), plane_names[plane].c_str(), paddle),
                     Form("Hit time difference for %s %s paddle %d", side_names[side].c_str(),
                          plane_names[plane].c_str(), paddle),
                     100, 0, 500);
        hit_time_diff_adc[side][plane][paddle]->GetXaxis()->SetTitle("Time difference (ns)");
        hit_time_diff_adc[side][plane][paddle]->GetYaxis()->SetTitle("Counts");

        hit_time_diff_tdc[side][plane][paddle] =
            new TH1F(Form("hit_time_diff_tdc_%s_%s_%d", side_names[side].c_str(), plane_names[plane].c_str(), paddle),
                     Form("Hit time difference for TDC %s %s paddle %d", side_names[side].c_str(),
                          plane_names[plane].c_str(), paddle),
                     100, 0, 1000);
        hit_time_diff_tdc[side][plane][paddle]->GetXaxis()->SetTitle("Time difference (ns)");
        hit_time_diff_tdc[side][plane][paddle]->GetYaxis()->SetTitle("Counts");
      }
    }
  }

  Long64_t nentries = T->GetEntries();
  // nentries = 10000; // For testing purposes, limit to 10,000 entries
  printf("Processing %lld entries...\n", nentries);
  for (Long64_t i = 0; i < nentries; ++i) {
    T->GetEntry(i);

    std::fill_n(&last_tdc_time[0][0][0], N_SIDES * N_PLANES * N_PADDLES, kBIG);
    std::fill_n(&last_adc_time[0][0][0], N_SIDES * N_PLANES * N_PADDLES, kBIG);

    for (int side = 0; side < N_SIDES; ++side) {
      for (int plane = 0; plane < N_PLANES; ++plane) {
        bool has_hit_tdc[N_PADDLES] = {false};
        bool has_hit_adc[N_PADDLES] = {false};
        for (int hit = 0; hit < nData_tdc[side][plane]; ++hit) {
          Double_t curr_time = tdc_time[side][plane][hit] * TDC2NS; // Convert TDC time to ns
          Int_t curr_counter = tdc_counter[side][plane][hit] - 1;
          if (curr_counter < 0 || curr_counter >= N_PADDLES) {
            printf("Warning: Invalid TDC counter %d for side %d, plane %d, hit "
                   "%d, time %f\n",
                   curr_counter, side, plane, hit, curr_time);
            continue;
          }
          if (last_tdc_time[side][plane][curr_counter] == kBIG) {
            last_tdc_time[side][plane][curr_counter] = curr_time;
            has_hit_tdc[curr_counter]                = true;
          } else {
            Double_t time_diff = curr_time - last_tdc_time[side][plane][curr_counter];
            // Fill the histogram with the time difference
            hit_time_diff_tdc[side][plane][curr_counter]->Fill(time_diff);
            last_tdc_time[side][plane][curr_counter] = curr_time;
          }
        }
        // Loop over all ADC hits for this side/plane
        for (int hit = 0; hit < nData_adc[side][plane]; ++hit) {
          Double_t curr_time = adc_time[side][plane][hit];
          Int_t curr_counter = adc_counter[side][plane][hit] - 1;
          if (curr_counter < 0 || curr_counter >= N_PADDLES) {
            printf("Warning: Invalid ADC counter %d for side %d, plane %d\n", curr_counter, side, plane);
            continue;
          }
          if (last_adc_time[side][plane][curr_counter] == kBIG) {
            last_adc_time[side][plane][curr_counter] = curr_time;
            has_hit_adc[curr_counter]                = true;
          } else {
            Double_t time_diff = curr_time - last_adc_time[side][plane][curr_counter];
            // Fill the histogram with the time difference
            hit_time_diff_adc[side][plane][curr_counter]->Fill(time_diff);
            last_adc_time[side][plane][curr_counter] = curr_time;
          }
        }
        // Count the total number of hits for each paddle
        for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
          if (has_hit_tdc[paddle]) {
            nData_tot_tdc[side][plane][paddle] = nData_tot_tdc[side][plane][paddle] + 1;
          }
          if (has_hit_adc[paddle]) {
            nData_tot_adc[side][plane][paddle] = nData_tot_adc[side][plane][paddle] + 1;
          }
        }
      }
    }

    if (i % 1000 == 0 || i == nentries - 1) {
      float progress = 100.0f * (i + 1) / nentries;
      printf("\rProgress: %.2f%%", progress);
      fflush(stdout);
    }
  }
  printf("\nProcessing complete.\n");

  TFile *outfile = TFile::Open("/home/ehingerl/hallc/analysis/lad_timing/files/output_histograms.root", "RECREATE");

  // Create separate canvases for ADC and TDC histograms for each plane and side
  for (int plane = 0; plane < N_PLANES; ++plane) {
    for (int side = 0; side < N_SIDES; ++side) {
      // ADC canvas
      std::string adc_canvas_name = Form("adc_canvas_%s_%s", plane_names[plane].c_str(), side_names[side].c_str());
      TCanvas *adc_canvas         = new TCanvas(adc_canvas_name.c_str(), adc_canvas_name.c_str(), 1800, 1200);
      adc_canvas->Divide(3, 4); // 3 columns x 4 rows = 12 pads

      for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
        int pad_num = paddle + 1;
        adc_canvas->cd(pad_num);
        hit_time_diff_adc[side][plane][paddle]->SetLineColor(kBlue);
        hit_time_diff_adc[side][plane][paddle]->Draw();

        // Fit with linear PDF: P(t) = a + b*t
        TF1 *lin_fit_adc = new TF1(Form("lin_fit_adc_%d_%d_%d", side, plane, paddle), "[0] + [1]*x", 21, 500);
        lin_fit_adc->SetParameters(1.0, 0.0); // Initial guesses
        lin_fit_adc->SetParNames("Intercept", "Slope");
        hit_time_diff_adc[side][plane][paddle]->Fit(lin_fit_adc, "RQ");
        lin_fit_adc->SetLineColor(kGreen + 2);
        lin_fit_adc->Draw("same");
        // std::cout << "ADC Linear Fit (side=" << side << ", plane=" << plane
        //       << ", paddle=" << paddle
        //       << "): Intercept=" << lin_fit_adc->GetParameter(0)
        //       << ", Slope=" << lin_fit_adc->GetParameter(1) << std::endl;
        // // Fit with exponential PDF: P(t) = lambda * exp(-lambda * t)
        // TF1 *exp_fit_adc = new TF1(Form("exp_fit_adc_%d_%d_%d", side, plane, paddle), "[0]*exp(-[1]*x)", 21, 500);
        // exp_fit_adc->SetParameters(1.0, 0.01); // Initial guesses
        // exp_fit_adc->SetParNames("Norm", "lambda");
        // hit_time_diff_adc[side][plane][paddle]->Fit(exp_fit_adc, "RQ");
        // exp_fit_adc->SetLineColor(kGreen + 2);
        // exp_fit_adc->Draw("same");
        // std::cout << "ADC Fit lambda (side=" << side << ", plane=" << plane
        //           << ", paddle=" << paddle
        //           << "): " << exp_fit_adc->GetParameter(1) << std::endl;

        // Draw normalized exponential function with lambda = hist->Integral() /
        // (x axis range) / nData_tot
        double x_min         = hit_time_diff_adc[side][plane][paddle]->GetXaxis()->GetXmin();
        double x_max         = hit_time_diff_adc[side][plane][paddle]->GetXaxis()->GetXmax();
        int first_bin        = 1;
        int last_bin         = hit_time_diff_adc[side][plane][paddle]->GetNbinsX();
        double hist_integral = hit_time_diff_adc[side][plane][paddle]->Integral(first_bin, last_bin);
        double bin_width     = hit_time_diff_adc[side][plane][paddle]->GetXaxis()->GetBinWidth(1);
        double lambda_custom = hist_integral / (x_max - x_min) / nData_tot_adc[side][plane][paddle];

        // std::cout << "Custom lambda (side=" << side << ", plane=" << plane
        //           << ", paddle=" << paddle << "): " << lambda_custom
        //           << std::endl;
        // Normalization so that the area under the exponential matches the
        // histogram integral

        // Calculate normalization so that the custom exponential has the same
        // integral as the fitted one
        // double fit_integral = exp_fit_adc->Integral(x_min, x_max);
        double fit_integral = lin_fit_adc->Integral(x_min, x_max);
        double custom_integral =
            (lambda_custom > 0) ? (1.0 / lambda_custom) * (exp(-lambda_custom * x_min) - exp(-lambda_custom * x_max))
                                : (x_max - x_min); // fallback for lambda=0

        double norm_custom = (custom_integral != 0) ? fit_integral / custom_integral : 1.0;

        // TF1 *custom_exp =
        //     new TF1(Form("custom_exp_adc_%d_%d_%d", side, plane, paddle),
        //             "[0]*exp(-[1]*x)", x_min, x_max);
        // custom_exp->SetParameters(norm_custom, lambda_custom);
        // custom_exp->SetLineColor(kOrange + 2);
        // custom_exp->SetLineStyle(2);
        // custom_exp->Draw("same");
        // // Draw scaler rate label
        // double scaler_rate = hist_integral / (x_max - x_min) /
        //                      nData_tot_adc[side][plane][paddle] *
        //                      pow(10, 6); // Convert to kHz
        // TLatex *rate_label = new TLatex(
        //     0.55 * x_max,
        //     0.95 * hit_time_diff_adc[side][plane][paddle]->GetMaximum(),
        //     Form("Scaler Rate: %.3g kHz", scaler_rate));
        // rate_label->SetTextSize(0.04);
        // rate_label->SetTextColor(kBlack);
        // rate_label->Draw();
      }
      // if (adc_canvas->GetPad(12)) {
      //   adc_canvas->cd(12);
      //   TLegend *legend = new TLegend(0.15, 0.7, 0.85, 0.9);
      //   legend->AddEntry(hit_time_diff_adc[side][plane][0], "Data", "l");
      //   legend->AddEntry("exp_fit_adc_0_0_0", "Fit", "l");
      //   legend->AddEntry("custom_exp_adc_0_0_0", "Expectation", "l");
      //   legend->Draw();
      // }
      adc_canvas->Write();
      delete adc_canvas;

      // TDC canvas
      std::string tdc_canvas_name = Form("tdc_canvas_%s_%s", plane_names[plane].c_str(), side_names[side].c_str());
      TCanvas *tdc_canvas         = new TCanvas(tdc_canvas_name.c_str(), tdc_canvas_name.c_str(), 1800, 1200);
      tdc_canvas->Divide(3, 4);

      for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
        int pad_num = paddle + 1;
        tdc_canvas->cd(pad_num);
        hit_time_diff_tdc[side][plane][paddle]->SetLineColor(kRed);
        hit_time_diff_tdc[side][plane][paddle]->Draw();

        // Fit with exponential PDF: P(t) = lambda * exp(-lambda * t)
        TF1 *exp_fit_tdc = new TF1(Form("exp_fit_tdc_%d_%d_%d", side, plane, paddle), "[0]*exp(-[1]*x)", 21, 1000);
        exp_fit_tdc->SetParameters(1.0, 0.01); // Initial guesses
        exp_fit_tdc->SetParNames("Norm", "lambda");
        hit_time_diff_tdc[side][plane][paddle]->Fit(exp_fit_tdc, "RQ");
        exp_fit_tdc->SetLineColor(kMagenta + 2);
        exp_fit_tdc->Draw("same");

        TF1 *custom_exp_tdc =
            new TF1(Form("custom_exp_tdc_%d_%d_%d", side, plane, paddle), "[0]*exp(-[1]*x)", 21, 1000);
        double x_min_tdc         = hit_time_diff_tdc[side][plane][paddle]->GetXaxis()->GetXmin();
        double x_max_tdc         = hit_time_diff_tdc[side][plane][paddle]->GetXaxis()->GetXmax();
        int first_bin_tdc        = 1;
        int last_bin_tdc         = hit_time_diff_tdc[side][plane][paddle]->GetNbinsX();
        double hist_integral_tdc = hit_time_diff_tdc[side][plane][paddle]->Integral(first_bin_tdc, last_bin_tdc);
        double bin_width_tdc     = hit_time_diff_tdc[side][plane][paddle]->GetXaxis()->GetBinWidth(1);
        double lambda_custom_tdc = hist_integral_tdc / (x_max_tdc - x_min_tdc) / nData_tot_tdc[side][plane][paddle] *
                                   bin_width_tdc * bin_width_tdc;

        double fit_integral_tdc = exp_fit_tdc->Integral(x_min_tdc, x_max_tdc);
        double custom_integral_tdc =
            (lambda_custom_tdc > 0) ? (1.0 / lambda_custom_tdc) *
                                          (exp(-lambda_custom_tdc * x_min_tdc) - exp(-lambda_custom_tdc * x_max_tdc))
                                    : (x_max_tdc - x_min_tdc);

        double norm_custom_tdc = (custom_integral_tdc != 0) ? fit_integral_tdc / custom_integral_tdc : 1.0;

        // custom_exp_tdc->SetParameters(norm_custom_tdc, lambda_custom_tdc);
        // custom_exp_tdc->SetLineColor(kOrange + 2);
        // custom_exp_tdc->SetLineStyle(2);
        // custom_exp_tdc->Draw("same");

        // // Draw scaler rate label
        // double scaler_rate = hist_integral_tdc / (x_max_tdc - x_min_tdc) /
        //                      nData_tot_tdc[side][plane][paddle] *
        //                      pow(10, 6); // Convert to kHz
        // TLatex *rate_label = new TLatex(
        //     0.55 * x_max_tdc,
        //     0.95 * hit_time_diff_tdc[side][plane][paddle]->GetMaximum(),
        //     Form("Scaler Rate: %.3g kHz", scaler_rate));
        // rate_label->SetTextSize(0.04);
        // rate_label->SetTextColor(kBlack);
        // rate_label->Draw();
      }
      // if (tdc_canvas->GetPad(12)) {
      //   tdc_canvas->cd(12);
      //   TLegend *legend = new TLegend(0.15, 0.7, 0.85, 0.9);
      //   legend->AddEntry(hit_time_diff_tdc[side][plane][0], "Data", "l");
      //   legend->AddEntry("exp_fit_tdc_0_0_0", "Fit", "l");
      //   legend->AddEntry("custom_exp_tdc_0_0_0", "Expectation", "l");
      //   legend->Draw();
      // }
      tdc_canvas->Write();
      delete tdc_canvas;
    }
  }

  // // Optionally, write the histograms individually as well
  // for (int plane = 0; plane < N_PLANES; ++plane) {
  //   for (int side = 0; side < N_SIDES; ++side) {
  //     for (int paddle = 0; paddle < N_PADDLES; ++paddle) {
  //       hit_time_diff_adc[side][plane][paddle]->Write();
  //       hit_time_diff_tdc[side][plane][paddle]->Write();
  //     }
  //   }
  // }
  outfile->Close();
  return;
}