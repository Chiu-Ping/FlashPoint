#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TLegend.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

void plot_dndt(const char *inputName = "sipm_timing.root",
               const char *allOutput = "dNdt_all_SiPMs.pdf",
               const char *overlayOutput = "dNdt_overlay.pdf") {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  std::unique_ptr<TFile> input(TFile::Open(inputName, "READ"));
  if (!input || input->IsZombie()) {
    std::cerr << "Cannot open " << inputName << std::endl;
    return;
  }

  auto *eventCounter = dynamic_cast<TH1 *>(input->Get("event_counter"));
  const double numberOfEvents = eventCounter ? eventCounter->Integral() : 0.0;
  if (numberOfEvents <= 0.0) {
    std::cerr << "event_counter is missing or empty" << std::endl;
    return;
  }

  std::vector<std::unique_ptr<TH1>> histograms;
  for (int sipm = 0;; ++sipm) {
    const std::string name = "sipm" + std::to_string(sipm) + "_time";
    auto *source = dynamic_cast<TH1 *>(input->Get(name.c_str()));
    if (!source)
      break;

    auto histogram = std::unique_ptr<TH1>(
        dynamic_cast<TH1 *>(source->Clone((name + "_dndt").c_str())));
    histogram->SetDirectory(nullptr);
    histogram->Scale(1.0 / (numberOfEvents * histogram->GetBinWidth(1)));
    histogram->SetTitle(("SiPM " + std::to_string(sipm) +
                         ";time [ns];dN_{PE}/dt [PE event^{-1} ns^{-1}]")
                            .c_str());
    histograms.push_back(std::move(histogram));
  }

  if (histograms.empty()) {
    std::cerr << "No sipmN_time histograms found" << std::endl;
    return;
  }

  const int columns =
      static_cast<int>(std::ceil(std::sqrt(histograms.size())));
  const int rows = static_cast<int>(
      std::ceil(static_cast<double>(histograms.size()) / columns));

  TCanvas allCanvas("allCanvas", "Per-SiPM dN/dt", 420 * columns,
                    330 * rows);
  allCanvas.Divide(columns, rows, 0.002, 0.002);
  for (std::size_t i = 0; i < histograms.size(); ++i) {
    allCanvas.cd(static_cast<int>(i) + 1);
    gPad->SetGrid();
    histograms[i]->SetLineWidth(2);
    histograms[i]->GetXaxis()->SetRangeUser(0, 40);
    histograms[i]->Draw("hist");
  }
  allCanvas.Print(allOutput);

  TCanvas overlayCanvas("overlayCanvas", "SiPM dN/dt overlay", 1100, 750);
  overlayCanvas.SetGrid();
  TLegend legend(0.78, 0.55, 0.94, 0.90);
  legend.SetNColumns(histograms.size() > 12 ? 2 : 1);
    histograms[0]->GetXaxis()->SetRangeUser(0, 40);
  for (std::size_t i = 0; i < histograms.size(); ++i) {
    histograms[i]->SetLineColor(TColor::GetColorPalette(
        static_cast<int>(i * 254 / std::max<std::size_t>(1, histograms.size() - 1))));
    histograms[i]->SetLineWidth(2);
    histograms[i]->Draw(i == 0 ? "hist" : "hist same");
    legend.AddEntry(histograms[i].get(),
                    ("SiPM " + std::to_string(i)).c_str(), "l");
  }
  legend.Draw();
  overlayCanvas.Print(overlayOutput);

  std::cout << "Plotted " << histograms.size() << " SiPMs from "
            << numberOfEvents << " events" << std::endl;
}
