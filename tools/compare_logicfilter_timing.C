// root -l -b -q 'tools/compare_logicfilter_timing.C("cpu.root","gpu.root","/tmp/lf-timing",10)'
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLeaf.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace lf_timing {
const std::array<const char*, 14> counts = {{"input_parts", "input_bytes", "hbf_count", "segment_count",
    "trigger_input_bytes", "trigger_input_words", "payload_words", "leading_tdc_words",
    "trailing_tdc_words", "tdc_words", "triggers", "output_parts", "output_bytes", "tf_id"}};
const std::array<const char*, 10> times = {{"conditional_us", "receive_us", "trigger_us",
    "cleanup_host_us", "exec_us", "send_us", "accounting_us", "processing_us", "other_us", "elapsed_no_receive_us"}};
struct Row {
    std::array<ULong64_t, 14> data;
    std::array<double, 10> time;
};
struct Run {
    std::vector<Row> rows;
    Long64_t calls{};
    std::string signals, expression;
    int width{};
};
TLeaf* leaf(TTree* tree, const char* name) {
    auto* value = tree->GetLeaf(name);
    if (!value) throw std::runtime_error(std::string("Missing branch: ") + name);
    return value;
}
Run read(const char* path, bool gpu) {
    std::unique_ptr<TFile> file(TFile::Open(path, "READ"));
    if (!file || file->IsZombie()) throw std::runtime_error(std::string("Cannot open ") + path);
    auto* completed = file->Get<TParameter<int>>("completed");
    auto* tree = file->Get<TTree>("logicfilter");
    auto* version = file->Get<TParameter<int>>("schema_version");
    if (!tree || !version || version->GetVal() != 1 || !completed || completed->GetVal() != 1)
        throw std::runtime_error(std::string("Incomplete/invalid profile: ") + path);
    Run run;
    auto* signals = file->Get<TNamed>("trigger_signals");
    auto* expression = file->Get<TNamed>("trigger_expression");
    auto* width = file->Get<TParameter<int>>("trigger_width");
    if (!signals || !expression || !width) throw std::runtime_error("Missing trigger configuration");
    run.signals = signals->GetTitle(); run.expression = expression->GetTitle(); run.width = width->GetVal();
    run.calls = tree->GetEntries();
    for (Long64_t i = 0; i < run.calls; ++i) {
        if (tree->GetEntry(i) <= 0) throw std::runtime_error("Cannot read profile entry");
        if (leaf(tree, "iteration")->GetValueLong64() != i) throw std::runtime_error("Nonconsecutive iteration");
        if (bool(leaf(tree, "cuda_enabled")->GetValue()) != gpu) throw std::runtime_error("Wrong CPU/GPU implementation");
        if (!leaf(tree, "received")->GetValue()) continue; // Exclude timeout calls from statistics.
        if (!leaf(tree, "sent")->GetValue()) throw std::runtime_error("Received TF was not sent");
        Row row;
        for (size_t j = 0; j < counts.size(); ++j) row.data[j] = leaf(tree, counts[j])->GetValueLong64();
        for (size_t j = 0; j < times.size() - 1; ++j) {
            row.time[j] = leaf(tree, times[j])->GetValue();
            if (!std::isfinite(row.time[j]) || row.time[j] < -1e-5) throw std::runtime_error("Invalid elapsed time");
        }
        if (row.data[9] != row.data[7] + row.data[8] || row.data[9] > row.data[6] || row.data[6] > row.data[5])
            throw std::runtime_error("Inconsistent word counts");
        if (std::abs(row.time[2] - row.time[3] - row.time[4]) > 1e-5 ||
            std::abs(row.time[7] - row.time[2] - row.time[8]) > 1e-5 ||
            std::abs(row.time[0] - row.time[1] - row.time[5] - row.time[6] - row.time[7]) > 1e-5)
            throw std::runtime_error("Inconsistent timing decomposition");
        row.time[9] = row.time[0] - row.time[1];
        run.rows.push_back(row);
    }
    if (run.rows.empty()) throw std::runtime_error("Profile has no received TFs");
    return run;
}
double quantile(std::vector<double> values, double fraction) {
    std::sort(values.begin(), values.end());
    return values[size_t(fraction * (values.size() - 1))];
}
}

void compare_logicfilter_timing(const char* cpuPath, const char* gpuPath,
                               const char* outputDir = "/tmp/logicfilter-timing", int warmup = 10, int triggerSelection = -1) {
    using namespace lf_timing;
    try {
        auto cpu = read(cpuPath, false), gpu = read(gpuPath, true);
        if (cpu.signals != gpu.signals || cpu.expression != gpu.expression || cpu.width != gpu.width)
            throw std::runtime_error("CPU/GPU trigger configurations differ");
        if (cpu.rows.size() != gpu.rows.size()) throw std::runtime_error("CPU/GPU TF counts differ");
        for (size_t i = 0; i < cpu.rows.size(); ++i)
            if (cpu.rows[i].data != gpu.rows[i].data)
                throw std::runtime_error("CPU/GPU data counts/TF ID differ at received row " + std::to_string(i));
        if (warmup < 0 || size_t(warmup) >= cpu.rows.size()) throw std::runtime_error("Invalid warmup count");
        if (triggerSelection < -1 || triggerSelection > 1) throw std::runtime_error("Trigger selection must be -1 (all), 0 (zero), or 1 (at least one)");
        std::vector<size_t> selected;
        for (size_t i = warmup; i < cpu.rows.size(); ++i)
            if (triggerSelection == -1 || (cpu.rows[i].data[10] > 0) == bool(triggerSelection)) selected.push_back(i);
        if (selected.empty()) throw std::runtime_error("No TFs satisfy trigger selection");
        const std::string selection = triggerSelection == -1 ? "all" : triggerSelection ? "triggers >= 1" : "triggers = 0";
        std::filesystem::create_directories(outputDir);
        const auto prefix = std::string(outputDir) + "/";
        std::ofstream summary(prefix + "summary.csv");
        summary << "metric,cpu_mean_us,gpu_mean_us,cpu_median_us,gpu_median_us,cpu_p95_us,gpu_p95_us,cpu_over_gpu,sample_count\n";
        if (!summary) throw std::runtime_error("Cannot create summary.csv");

        for (size_t j = 0; j < times.size(); ++j) {
            std::array<std::vector<double>, 2> values;
            for (size_t i : selected) {
                values[0].push_back(cpu.rows[i].time[j]); values[1].push_back(gpu.rows[i].time[j]);
            }
            auto mean = [](const std::vector<double>& v) { return std::accumulate(v.begin(), v.end(), 0.0) / v.size(); };
            const auto a = mean(values[0]), b = mean(values[1]);
            summary << times[j] << ',' << a << ',' << b << ',' << quantile(values[0], .5) << ',' << quantile(values[1], .5)
                    << ',' << quantile(values[0], .95) << ',' << quantile(values[1], .95) << ',' << (b > 0 ? a/b : 0) << ',' << selected.size() << '\n';
            std::cout << times[j] << ": CPU=" << a << " us, GPU=" << b << " us, CPU/GPU=" << (b > 0 ? a/b : 0) << '\n';

        }
        TFile graphs((prefix + "comparison.root").c_str(), "RECREATE");
        TNamed("selection", selection.c_str()).Write();
        TParameter<Long64_t>("selected_tfs", selected.size()).Write();
        TParameter<int>("warmup_tfs", warmup).Write();
        for (bool triggerOnly : {false, true}) {
            const size_t metric = triggerOnly ? 2 : 9;
            const std::string stem = triggerOnly ? "trigger" : "elapsed_no_receive";
            const std::string title = (triggerOnly ? std::string("LogicFilter trigger processing") : std::string("LogicFilter elapsed time excluding Receive")) + " (" + selection + ")";
            const std::string yLabel = triggerOnly ? "Trigger processing elapsed time [us/TF]" : "Elapsed time excluding Receive [us/TF]";
            TGraph cpuGraph, gpuGraph;
            cpuGraph.SetName(("cpu_" + stem + "_vs_tdc").c_str());
            gpuGraph.SetName(("gpu_" + stem + "_vs_tdc").c_str());
            cpuGraph.SetMarkerColor(kBlue); gpuGraph.SetMarkerColor(kRed);
            cpuGraph.SetLineWidth(3); gpuGraph.SetLineWidth(3);
            cpuGraph.SetMarkerStyle(7); gpuGraph.SetMarkerStyle(7);
            double maximum = 0;
            for (size_t i : selected) {
                // Only type 0x0b words are eligible for Mark's channel/time lookup.
                cpuGraph.SetPoint(cpuGraph.GetN(), cpu.rows[i].data[7], cpu.rows[i].time[metric]);
                gpuGraph.SetPoint(gpuGraph.GetN(), gpu.rows[i].data[7], gpu.rows[i].time[metric]);
                maximum = std::max({maximum, cpu.rows[i].time[metric], gpu.rows[i].time[metric]});
            }
            TCanvas scatter((stem + "_scatter").c_str(), title.c_str(), 800, 600);
            scatter.SetLeftMargin(.14);
            TMultiGraph multi;
            multi.Add(static_cast<TGraph*>(cpuGraph.Clone()), "P");
            multi.Add(static_cast<TGraph*>(gpuGraph.Clone()), "P");
            multi.SetTitle((title + ";TDC+TOT words per TF (type 0x0b);" + yLabel).c_str());
            multi.Draw("A");
            TLegend legend(.18, .75, .38, .88);
            legend.AddEntry(&cpuGraph, "CPU", "p"); legend.AddEntry(&gpuGraph, "GPU", "p"); legend.Draw();
            scatter.SaveAs((prefix + stem + "_vs_tdc.png").c_str());
            scatter.SaveAs((prefix + stem + "_vs_tdc.pdf").c_str());
            cpuGraph.Write(); gpuGraph.Write(); scatter.Write();
            TH1D cpuHist(("cpu_" + stem + "_us").c_str(), (title + ";" + yLabel + ";TFs").c_str(), 200, 0, maximum * 1.01 + 1);
            TH1D gpuHist(("gpu_" + stem + "_us").c_str(), (title + ";" + yLabel + ";TFs").c_str(), 200, 0, maximum * 1.01 + 1);
            for (size_t i : selected) {
                cpuHist.Fill(cpu.rows[i].time[metric]); gpuHist.Fill(gpu.rows[i].time[metric]);
            }
            cpuHist.SetLineColor(kBlue); gpuHist.SetLineColor(kRed);
            cpuHist.SetLineWidth(3); gpuHist.SetLineWidth(3);
            TCanvas distribution((stem + "_distribution").c_str(), title.c_str(), 800, 600);
            cpuHist.SetMaximum(std::max(cpuHist.GetMaximum(), gpuHist.GetMaximum()) * 1.15);
            cpuHist.Draw("HIST"); gpuHist.Draw("HIST SAME"); legend.Draw();
            distribution.SaveAs((prefix + stem + "_distribution.png").c_str());
            distribution.SaveAs((prefix + stem + "_distribution.pdf").c_str());
            cpuHist.Write(); gpuHist.Write(); distribution.Write();
        }
        std::cout << "MATCHED DATA: " << cpu.rows.size() << " TFs; excluded " << warmup
                  << " warmup TFs; selected " << selected.size() << " TFs (" << selection << "); ConditionalRun calls CPU=" << cpu.calls << ", GPU=" << gpu.calls << '\n';
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        gSystem->Exit(1);
    }
}
