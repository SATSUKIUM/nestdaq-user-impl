// Analyze GPU Exec internals from an instrumented run and a control run.
#include "compare_logicfilter_timing.C"
#include "../LogicFilter/ExecProfileData.h"
#include <THStack.h>

void analyze_gpu_exec(const char* detailPath, const char* controlPath,
    const char* outputDir = "/tmp/gpu-exec-breakdown", int warmup = 10) {
    using namespace lf_timing;
    try {
        auto detailed = read(detailPath, true), control = read(controlPath, true);
        if (detailed.rows.size() != control.rows.size() || detailed.signals != control.signals ||
            detailed.expression != control.expression || detailed.width != control.width)
            throw std::runtime_error("Detailed/control inputs or configurations differ");
        for (size_t i = 0; i < detailed.rows.size(); ++i)
            if (detailed.rows[i].data != control.rows[i].data) throw std::runtime_error("Detailed/control TF data differ");
        if (warmup < 0 || size_t(warmup) >= detailed.rows.size()) throw std::runtime_error("Invalid warmup count");
        std::vector<std::string> hostNames, deviceNames, countNames;
#define PUSH_HOST(field) hostNames.push_back(#field);
        LF_EXEC_HOST_FIELDS(PUSH_HOST)
#undef PUSH_HOST
        hostNames.pop_back(); // total_host_us is the denominator, not another stage.
        hostNames.push_back("outer_wrapper_host_us");
#define PUSH_DEVICE(field) deviceNames.push_back(#field);
        LF_EXEC_DEVICE_FIELDS(PUSH_DEVICE)
#undef PUSH_DEVICE
#define PUSH_COUNT(field) countNames.push_back(#field);
        LF_EXEC_COUNT_FIELDS(PUSH_COUNT)
#undef PUSH_COUNT
        struct DetailRow { std::vector<double> host, device, count; };
        std::vector<DetailRow> rows;
        std::unique_ptr<TFile> input(TFile::Open(detailPath));
        auto* enabled = input->Get<TParameter<int>>("exec_detail_enabled");
        auto* tree = input->Get<TTree>("logicfilter");
        auto* eventMode = input->Get<TParameter<int>>("exec_cuda_events_enabled");
        const bool cudaEvents = !eventMode || eventMode->GetVal() != 0;
        if (!enabled || enabled->GetVal() != 1) throw std::runtime_error("Exec detail profiling was not enabled");
        for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
            tree->GetEntry(entry);
            if (!leaf(tree, "received")->GetValue()) continue;
            if (!leaf(tree, "gpu_exec_enabled")->GetValue()) throw std::runtime_error("Missing Exec details");
            DetailRow row;
            auto value = [&](const std::string& name) { return leaf(tree, ("gpu_exec_" + name).c_str())->GetValue(); };
            const double internalTotal = value("total_host_us"), outerTotal = leaf(tree, "exec_us")->GetValue();
            for (size_t j = 0; j < hostNames.size() - 1; ++j) row.host.push_back(value(hostNames[j]));
            row.host.push_back(outerTotal - internalTotal);
            for (const auto& name : deviceNames) row.device.push_back(value(name));
            for (const auto& name : countNames) row.count.push_back(value(name));
            for (const auto& values : {row.host, row.device, row.count})
                for (double x : values)
                    if (!std::isfinite(x) || x < -1e-4) throw std::runtime_error("Invalid Exec detail value");
            if (std::abs(std::accumulate(row.host.begin(), row.host.end(), 0.0) - outerTotal) > 1e-5)
                throw std::runtime_error("Host stage sum does not reproduce outer Exec time");
            if (row.count[0] != leaf(tree, "segment_count")->GetValue() ||
                row.count[1] != 8 * leaf(tree, "trigger_input_words")->GetValue() ||
                row.count[2] != 4 * leaf(tree, "hbf_count")->GetValue() ||
                row.count[3] != 4 * leaf(tree, "triggers")->GetValue())
                throw std::runtime_error("Exec transfer counters disagree with processed input/output");
            rows.push_back(std::move(row));
        }
        if (rows.size() != detailed.rows.size()) throw std::runtime_error("Detail row count mismatch");
        input.reset();
        std::array<std::vector<size_t>, 3> groups;
        for (size_t i = warmup; i < rows.size(); ++i) {
            groups[detailed.rows[i].data[10] > 0 ? 0 : 1].push_back(i); groups[2].push_back(i);
        }
        for (const auto& group : groups) if (group.empty()) throw std::runtime_error("Empty trigger group");
        const std::array<std::string, 3> groupNames = {{"triggers_ge_1", "triggers_eq_0", "all"}};
        std::filesystem::create_directories(outputDir);
        const std::string prefix = std::string(outputDir) + "/";
        std::ofstream summary(prefix + "exec_stages.csv"), countsCsv(prefix + "exec_counts.csv"), overhead(prefix + "profiling_overhead.csv");
        summary << "group,tf_count,domain,stage,mean_us,median_us,p95_us,share_of_exec_percent\n";
        countsCsv << "group,tf_count,counter,mean_per_tf,total\n";
        overhead << "group,tf_count,control_exec_mean_us,detailed_exec_mean_us,relative_change_percent,instrumentation_host_mean_us\n";
        if (!summary || !countsCsv || !overhead) throw std::runtime_error("Cannot create Exec summaries");
        auto mean = [](const std::vector<double>& v) { return std::accumulate(v.begin(), v.end(), 0.0) / v.size(); };
        std::array<std::vector<double>, 3> hostMeans, deviceMeans;
        std::array<double, 3> execMeans{};
        for (size_t group = 0; group < groups.size(); ++group) {
            std::vector<double> exec, reference;
            for (auto i : groups[group]) { exec.push_back(detailed.rows[i].time[4]); reference.push_back(control.rows[i].time[4]); }
            execMeans[group] = mean(exec);
            for (int domain = 0; domain < (cudaEvents ? 2 : 1); ++domain) {
                const auto& names = domain ? deviceNames : hostNames;
                for (size_t stage = 0; stage < names.size(); ++stage) {
                    std::vector<double> values;
                    for (auto i : groups[group]) values.push_back(domain ? rows[i].device[stage] : rows[i].host[stage]);
                    const double average = mean(values);
                    (domain ? deviceMeans[group] : hostMeans[group]).push_back(average);
                    summary << groupNames[group] << ',' << groups[group].size() << ',' << (domain ? "cuda_event" : "host") << ','
                            << names[stage] << ',' << average << ',' << quantile(values, .5) << ',' << quantile(values, .95) << ',';
                    if (!domain) summary << average / execMeans[group] * 100;
                    summary << '\n';
                }
            }
            for (size_t stage = 0; stage < countNames.size(); ++stage) {
                double sum = 0; for (auto i : groups[group]) sum += rows[i].count[stage];
                countsCsv << groupNames[group] << ',' << groups[group].size() << ',' << countNames[stage] << ',' << sum / groups[group].size() << ',' << sum << '\n';
            }
            overhead << groupNames[group] << ',' << groups[group].size() << ',' << mean(reference) << ',' << mean(exec)
                     << ',' << (mean(exec) / mean(reference) - 1) * 100 << ',' << hostMeans[group][11] << '\n';
            std::cout << groupNames[group] << ": Exec control=" << mean(reference) << " us, detailed=" << mean(exec)
                      << " us; instrumentation=" << hostMeans[group][11] << " us\n";
        }
        TFile output((prefix + "gpu_exec.root").c_str(), "RECREATE");
        TParameter<int>("cuda_events_enabled", cudaEvents).Write();
        TNamed("host_note", "Host intervals sum to Exec wall time. instrumentation_host_us covers event/setup calls; timer overhead and scheduling perturbation are not fully removed.").Write();
        TNamed("device_note", "CUDA-event stream intervals overlap host time; include scheduling gaps; do not add to host intervals").Write();
        // Group the disjoint host stages into readable stacks; full detail remains in CSV/ROOT.
        const std::array<const char*, 6> labels = {{"Input copy + buffer", "Mark launch + synchronization", "Scan pipeline submission",
            "Result copies + resize", "Event/setup measurement calls", "Other + outer wrapper"}};
        const std::array<int, 6> colors = {{kAzure+1, kOrange+1, kGreen+1, kMagenta-3, kGray+1, kYellow+1}};
        std::array<std::array<double, 6>, 3> grouped{};
        for (size_t g = 0; g < 3; ++g) {
            const auto& h = hostMeans[g];
            grouped[g] = {{h[0]+h[1], h[2]+h[3], h[4]+h[5]+h[6]+h[7], h[8]+h[9]+h[10], h[11], h[12]+h[13]}};
        }
        for (bool fraction : {false, true}) {
            const std::string stem = fraction ? "host_fraction" : "host_time";
            std::vector<std::unique_ptr<TH1D>> histograms;
            TCanvas canvas(stem.c_str(), "GPU Exec host breakdown", 800, 600);
            canvas.SetLeftMargin(.13); canvas.SetBottomMargin(.14); canvas.SetTopMargin(.28);
            THStack stack((stem + "_stack").c_str(), ((std::string("GPU Exec host intervals (") + (cudaEvents ? "CUDA events enabled" : "host timers only") + ");;") +
                std::string(fraction ? "Share of Exec elapsed time [%]" : "Mean Exec elapsed time [us/TF]")).c_str());
            TLegend legend(.12, .75, .89, .9); legend.SetNColumns(2); legend.SetTextSize(.025);
            for (size_t stage = 0; stage < labels.size(); ++stage) {
                auto h = std::make_unique<TH1D>((stem + "_" + std::to_string(stage)).c_str(), labels[stage], 3, .5, 3.5);
                h->SetDirectory(nullptr); h->SetFillColor(colors[stage]); h->SetLineColor(colors[stage]); h->SetLineWidth(3);
                h->SetBarWidth(.65); h->SetBarOffset(.175);
                for (size_t g = 0; g < 3; ++g) {
                    h->GetXaxis()->SetBinLabel(g+1, (g == 0 ? "triggers #geq 1" : g == 1 ? "triggers = 0" : "All TFs"));
                    h->SetBinContent(g+1, fraction ? grouped[g][stage] / execMeans[g] * 100 : grouped[g][stage]);
                }
                stack.Add(h.get()); legend.AddEntry(h.get(), labels[stage], "f"); h->Write(); histograms.push_back(std::move(h));
            }
            stack.SetMaximum(fraction ? 105 : *std::max_element(execMeans.begin(), execMeans.end()) * 1.08);
            stack.Draw("BAR"); legend.Draw();
            canvas.SaveAs((prefix + stem + ".png").c_str()); canvas.SaveAs((prefix + stem + ".pdf").c_str());
            stack.Write(); canvas.Write();
        }
        const std::array<const char*, 8> deviceLabels = {{"Input H2D", "Mark", "SubTR", "LUT", "Excl. scan", "Compaction", "Count D2H", "Hits D2H"}};
        for (size_t g = 0; cudaEvents && g < 3; ++g) {
            const std::string name = "cuda_intervals_" + groupNames[g];
            TH1D h(name.c_str(), ("CUDA-event intervals (" + groupNames[g] + ");;Mean stream interval [us/TF]").c_str(), 8, .5, 8.5);
            h.SetDirectory(nullptr); h.SetFillColor(kAzure+1); h.SetLineColor(kAzure+1); h.SetLineWidth(3);
            h.SetBarWidth(.7); h.SetBarOffset(.15); h.SetStats(false);
            for (size_t stage = 0; stage < deviceLabels.size(); ++stage) {
                h.GetXaxis()->SetBinLabel(stage+1, deviceLabels[stage]); h.SetBinContent(stage+1, deviceMeans[g][stage]);
            }
            TCanvas canvas((name + "_canvas").c_str(), "CUDA event intervals", 800, 600);
            canvas.SetLeftMargin(.14); canvas.SetBottomMargin(.15);
            h.GetXaxis()->SetLabelSize(.032); h.Draw("BAR");
            canvas.SaveAs((prefix + name + ".png").c_str()); canvas.SaveAs((prefix + name + ".pdf").c_str());
            h.Write(); canvas.Write();
        }
        // Detailed phase distributions retain the requested 200-bin convention.
        for (size_t stage = 0; stage < hostNames.size(); ++stage) {
            double maximum = 0; for (auto i : groups[0]) maximum = std::max(maximum, rows[i].host[stage]);
            TH1D h(("host_distribution_" + hostNames[stage]).c_str(), hostNames[stage].c_str(), 200, 0, maximum * 1.01 + 1);
            h.SetLineWidth(3); for (auto i : groups[0]) h.Fill(rows[i].host[stage]); h.Write();
        }
        std::cout << "PASS: stage sums, transfer counts, matched TF data; " << rows.size() << " TFs, warmup=" << warmup << '\n';
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n'; gSystem->Exit(1);
    }
}
