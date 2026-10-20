#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TNamed.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

void plot_logicfilter_elapsed(const char* input = "test/logicfilter-elapsed-0.txt") {
    std::ifstream stream(input);
    if (!stream) throw std::runtime_error("Cannot open elapsed time text");
    std::vector<double> values;
    for (std::string line; std::getline(stream, line);) {
        if (line.empty() || line[0] == '#') continue;
        double value;
        std::istringstream row(line);
        if (!(row >> value) || !std::isfinite(value) || value < 0)
            throw std::runtime_error("Invalid elapsed time");
        values.push_back(value);
    }
    if (values.empty()) throw std::runtime_error("No elapsed times");
    std::sort(values.begin(), values.end());
    const double mean = std::accumulate(values.begin(), values.end(), 0.) / values.size();
    const size_t middle = values.size() / 2;
    const double median = values.size() % 2 ? values[middle] : (values[middle - 1] + values[middle]) / 2;
    const double p95 = values[size_t(.95 * (values.size() - 1))];
    const double p99 = values[size_t(.99 * (values.size() - 1))];
    const auto path = std::filesystem::path(input);
    const auto prefix = (path.parent_path() / (path.stem().string() + "-distribution")).string();
    const double upper = std::ceil(std::max(1., values.back() * 1.02) / 100.) * 100.;
    gStyle->SetOptStat(0);
    TFile output((prefix + ".root").c_str(), "RECREATE");
    TNamed("source", input).Write();
    TNamed("measurement", "Microseconds; Receive excluded; Send included; all samples, no warmup exclusion").Write();
    TH1D full("elapsed_all", "ConditionalRun elapsed time;Elapsed time excluding Receive [#mus];TFs / bin", 200, 0, upper);
    TH1D zoom("elapsed_zoom", "ConditionalRun elapsed time (main range);Elapsed time excluding Receive [#mus];TFs / bin", 200, 500, 1000);
    for (double x : values) { full.Fill(x); zoom.Fill(x); }
    for (auto* h : {&full, &zoom}) { h->SetLineWidth(3); h->SetLineColor(kAzure + 2); }
    TCanvas canvas("elapsed_canvas", "ConditionalRun elapsed time", 800, 600);
    canvas.SetLeftMargin(.13); canvas.SetBottomMargin(.13);
    auto draw = [&](TH1D& h, bool logarithmic, const std::string& suffix) {
        canvas.Clear(); canvas.SetLogy(logarithmic);
        h.SetMinimum(logarithmic ? .5 : 0.);
        h.SetMaximum(h.GetMaximum() * (logarithmic ? 5. : 1.2));
        h.Draw("HIST");
        TPaveText stats(.52, .62, .88, .87, "NDC");
        stats.SetFillColor(kWhite); stats.SetBorderSize(0); stats.SetTextAlign(12); stats.SetTextSize(.029);
        stats.AddText(TString::Format("TFs: %zu (all samples)", values.size()));
        stats.AddText(TString::Format("Mean: %.2f #mus", mean));
        stats.AddText(TString::Format("Median: %.2f #mus", median));
        stats.AddText(TString::Format("95%%: %.2f #mus, 99%%: %.2f #mus", p95, p99));
        stats.AddText(TString::Format("Maximum: %.2f #mus", values.back()));
        if (!logarithmic)
            stats.AddText(TString::Format("Outside displayed range: %.0f TFs", h.GetBinContent(0) + h.GetBinContent(201)));
        stats.Draw();
        canvas.SaveAs((prefix + suffix + ".png").c_str());
        canvas.SaveAs((prefix + suffix + ".pdf").c_str());
        canvas.Write((std::string("canvas") + suffix).c_str());
    };
    draw(full, true, "");
    draw(zoom, false, "-zoom");
    full.Write(); zoom.Write(); output.Close();
    std::ofstream summary(prefix + "-summary.txt");
    summary << std::setprecision(12) << "samples " << values.size() << '\n'
            << "mean_us " << mean << '\n' << "median_us " << median << '\n'
            << "p95_us " << p95 << '\n' << "p99_us " << p99 << '\n'
            << "min_us " << values.front() << '\n' << "max_us " << values.back() << '\n';
    std::cout << "Plotted all " << values.size() << " samples; mean=" << mean
              << " us, median=" << median << " us; prefix=" << prefix << '\n';
}
