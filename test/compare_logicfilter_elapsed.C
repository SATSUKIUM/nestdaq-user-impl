#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TNamed.h>
#include <TStyle.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

void compare_logicfilter_elapsed() {
    struct Sample {
        std::string key, label;
        int color;
        std::vector<double> values;
        double mean, median;
    };
    std::vector<Sample> samples = {
        {"cpu", "CPU", kBlack, {}, 0, 0},
        {"hbf1", "GPU: 1 HBF per batch", kAzure + 2, {}, 0, 0},
        {"hbf4", "GPU: 4 HBFs per batch", kOrange + 7, {}, 0, 0}};
    const std::string prefix = "test/logicfilter-elapsed-comparison";
    std::ofstream summary(prefix + "-summary.csv");
    summary << "mode,samples,mean_us,median_us,p95_us,p99_us,min_us,max_us\n" << std::setprecision(15);
    for (auto& s : samples) {
        const std::string path = "test/logicfilter-elapsed-0-" + s.key + ".txt";
        std::ifstream input(path);
        if (!input) throw std::runtime_error("Cannot open " + path);
        for (std::string line; std::getline(input, line);) {
            if (line.empty() || line[0] == '#') continue;
            double x;
            std::istringstream row(line);
            if (!(row >> x) || !std::isfinite(x) || x < 0) throw std::runtime_error("Invalid sample in " + path);
            s.values.push_back(x);
        }
        if (s.values.empty()) throw std::runtime_error("Empty sample: " + path);
        std::sort(s.values.begin(), s.values.end());
        const size_t n = s.values.size();
        s.mean = std::accumulate(s.values.begin(), s.values.end(), 0.) / n;
        s.median = n % 2 ? s.values[n/2] : (s.values[n/2-1] + s.values[n/2]) / 2.;
        summary << s.key << ',' << n << ',' << s.mean << ',' << s.median << ','
                << s.values[size_t(.95*(n-1))] << ',' << s.values[size_t(.99*(n-1))] << ','
                << s.values.front() << ',' << s.values.back() << '\n';
        std::cout << s.label << ": N=" << n << ", mean=" << s.mean << " us, median=" << s.median << " us\n";
    }
    gStyle->SetOptStat(0);
    TFile output((prefix + ".root").c_str(), "RECREATE");
    TNamed("measurement", "All samples; ConditionalRun; Receive excluded; Send included; no warmup or outlier removal").Write();
    TCanvas canvas("comparison", "ConditionalRun elapsed time comparison", 800, 600);
    canvas.SetLeftMargin(.13); canvas.SetBottomMargin(.13); canvas.SetTopMargin(.09);
    auto draw = [&](const std::string& suffix, double low, double high, bool logy) {
        canvas.Clear(); canvas.SetLogy(logy);
        std::vector<std::unique_ptr<TH1D>> histograms;
        double peak = 0;
        for (const auto& s : samples) {
            auto h = std::make_unique<TH1D>((s.key + suffix).c_str(),
                "ConditionalRun elapsed time;Elapsed time excluding Receive [#mus];TFs [% / bin]", 200, low, high);
            for (double x : s.values) h->Fill(x);
            h->Scale(100. / s.values.size());
            h->SetLineWidth(3); h->SetLineColor(s.color);
            peak = std::max(peak, h->GetMaximum());
            histograms.push_back(std::move(h));
        }
        TLegend legend(.40, .69, .89, .89);
        legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextSize(.027);
        for (size_t i = 0; i < samples.size(); ++i) {
            auto& h = histograms[i];
            h->SetMinimum(logy ? .002 : 0.);
            h->SetMaximum(peak * (logy ? 6. : 1.4));
            h->Draw(i ? "HIST SAME" : "HIST");
            legend.AddEntry(h.get(), TString::Format("%s (mean %.2f #mus)", samples[i].label.c_str(), samples[i].mean), "l");
            h->Write();
        }
        legend.Draw();
        TLatex annotation; annotation.SetNDC(); annotation.SetTextSize(.025);
        annotation.DrawLatex(.14, .93, "Receive excluded; Send included; all samples retained");
        for (size_t i = 0; i < samples.size(); ++i) {
            size_t outside = 0;
            for (double x : samples[i].values) if (x < low || x >= high) ++outside;
            annotation.SetTextColor(samples[i].color);
            annotation.DrawLatex(.43, .65 - i*.036, TString::Format("%s: N=%zu, outside range=%zu",
                samples[i].label.c_str(), samples[i].values.size(), outside));
        }
        canvas.RedrawAxis();
        canvas.SaveAs((prefix + suffix + ".png").c_str());
        canvas.SaveAs((prefix + suffix + ".pdf").c_str());
        canvas.Write(("canvas" + suffix).c_str());
    };
    draw("", 0, 8000, false);
    draw("-zoom", 400, 1200, false);
    draw("-full", 0, 8000, false);
    output.Close();
}
