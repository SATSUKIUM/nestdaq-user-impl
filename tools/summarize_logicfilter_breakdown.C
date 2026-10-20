// Summarize each implementation's own ConditionalRun time; no new DAQ run.
#include "compare_logicfilter_timing.C"
#include <THStack.h>

void summarize_logicfilter_breakdown(const char* cpuPath, const char* gpuPath,
    const char* outputDir = "/tmp/logicfilter-breakdown", int warmup = 10) {
    using namespace lf_timing;
    try {
        auto cpu = read(cpuPath, false), gpu = read(gpuPath, true);
        if (cpu.rows.size() != gpu.rows.size() || cpu.signals != gpu.signals ||
            cpu.expression != gpu.expression || cpu.width != gpu.width)
            throw std::runtime_error("CPU/GPU inputs/configuration differ");
        if (warmup < 0 || size_t(warmup) >= cpu.rows.size()) throw std::runtime_error("Invalid warmup count");
        for (size_t i = 0; i < cpu.rows.size(); ++i)
            if (cpu.rows[i].data != gpu.rows[i].data) throw std::runtime_error("CPU/GPU TF data differ");
        std::filesystem::create_directories(outputDir);
        const std::string prefix = std::string(outputDir) + "/";
        std::ofstream csv(prefix + "breakdown.csv");
        csv << "implementation,group,tf_count,denominator,stage,mean_us,percent\n";
        if (!csv) throw std::runtime_error("Cannot write breakdown.csv");
        TFile output((prefix + "breakdown.root").c_str(), "RECREATE");
        TNamed("note", "Receive-excluded times include accounting; ROOT Fill itself is excluded. GPU cleanup is host-call time.").Write();
        TParameter<int>("warmup_tfs", warmup).Write();
        const std::array<const char*, 6> stageNames = {{"exec", "cleanup_host", "other", "send", "accounting", "receive"}};
        const std::array<const char*, 6> stageLabels = {{"Exec (incl. transfer/sync on GPU)", "Cleanup (host call)",
            "Other processing", "Send incl. wait", "Accounting (instrumentation)", "Receive incl. wait"}};
        const std::array<size_t, 6> metric = {{4, 3, 8, 5, 6, 1}};
        const std::array<int, 6> colors = {{kAzure+1, kOrange+1, kGreen+1, kMagenta-3, kGray+1, kYellow+1}};
        const std::array<const char*, 3> groupNames = {{"triggers_ge_1", "triggers_eq_0", "all"}};
        for (int variant : {0, 1}) {
            const auto& run = variant ? gpu : cpu;
            const std::string implementation = variant ? "GPU" : "CPU";
            std::array<std::vector<size_t>, 3> groups;
            for (size_t i = warmup; i < run.rows.size(); ++i) {
                groups[run.rows[i].data[10] > 0 ? 0 : 1].push_back(i);
                groups[2].push_back(i);
            }
            std::array<std::array<double, 6>, 3> means{};
            std::array<double, 3> elapsed{}, conditional{};
            for (size_t group = 0; group < groups.size(); ++group) {
                if (groups[group].empty()) throw std::runtime_error("Empty trigger group");
                for (auto i : groups[group]) {
                    for (size_t stage = 0; stage < metric.size(); ++stage) means[group][stage] += run.rows[i].time[metric[stage]];
                    elapsed[group] += run.rows[i].time[9];
                    conditional[group] += run.rows[i].time[0];
                }
                for (auto& value : means[group]) value /= groups[group].size();
                elapsed[group] /= groups[group].size(); conditional[group] /= groups[group].size();
                const auto total = std::accumulate(means[group].begin(), means[group].end(), 0.0);
                if (std::abs(total - conditional[group]) > 1e-5 ||
                    std::abs(total - means[group][5] - elapsed[group]) > 1e-5)
                    throw std::runtime_error("Stage totals do not reproduce ConditionalRun time");
                for (bool includeReceive : {false, true}) {
                    const double denominator = includeReceive ? conditional[group] : elapsed[group];
                    for (size_t stage = 0; stage < (includeReceive ? 6u : 5u); ++stage)
                        csv << implementation << ',' << groupNames[group] << ',' << groups[group].size() << ','
                            << (includeReceive ? "conditional_us" : "elapsed_no_receive_us") << ',' << stageNames[stage]
                            << ',' << means[group][stage] << ',' << 100 * means[group][stage] / denominator << '\n';
                }
                std::cout << implementation << ' ' << groupNames[group] << ": n=" << groups[group].size()
                    << ", excluding Receive=" << elapsed[group] << " us, Exec=" << means[group][0]
                    << " us (" << 100 * means[group][0] / elapsed[group] << "%), Cleanup host=" << means[group][1]
                    << " us (" << 100 * means[group][1] / elapsed[group] << "%)\n";
            }
            // Categorical bins represent the three groups, not elapsed-time histogram bins.
            for (bool fraction : {false, true}) {
                const std::string name = implementation + (fraction ? "_fraction" : "_time");
                std::vector<std::unique_ptr<TH1D>> histograms;
                TCanvas canvas(name.c_str(), "ConditionalRun stage breakdown", 800, 600);
                canvas.SetLeftMargin(.12); canvas.SetBottomMargin(.14); canvas.SetTopMargin(.27);
                THStack stack((name + "_stack").c_str(), (implementation + " ConditionalRun breakdown (Receive excluded);;"
                    + (fraction ? std::string("Share of elapsed time [%]") : std::string("Mean elapsed time [ms/TF]"))).c_str());
                TLegend legend(.12, .75, .89, .9);
                legend.SetNColumns(2);
                legend.SetTextSize(.027);
                for (size_t stage = 0; stage < 5; ++stage) {
                    auto h = std::make_unique<TH1D>((name + "_" + stageNames[stage]).c_str(), stageLabels[stage], 3, .5, 3.5);
                    h->SetDirectory(nullptr); h->SetLineWidth(3);
                    h->SetBarWidth(.65); h->SetBarOffset(.175);
                    h->SetLineColor(colors[stage]); h->SetFillColor(colors[stage]);
                    for (size_t group = 0; group < groups.size(); ++group) {
                        const std::string label = (group == 0 ? "triggers #geq 1" : group == 1 ? "triggers = 0" : "All TFs")
                            + std::string(" (n=") + std::to_string(groups[group].size()) + ")";
                        h->GetXaxis()->SetBinLabel(group+1, label.c_str());
                        h->SetBinContent(group+1, fraction ? 100 * means[group][stage] / elapsed[group] : means[group][stage] / 1000);
                    }
                    stack.Add(h.get()); legend.AddEntry(h.get(), stageLabels[stage], "f");
                    h->Write(); histograms.push_back(std::move(h));
                }
                stack.SetMaximum(fraction ? 105 : *std::max_element(elapsed.begin(), elapsed.end()) / 1000 * 1.08);
                stack.Draw("BAR");
                stack.GetXaxis()->SetLabelSize(.038);
                stack.GetYaxis()->SetTitleSize(.04);
                legend.Draw();
                canvas.SaveAs((prefix + name + "_breakdown.png").c_str());
                canvas.SaveAs((prefix + name + "_breakdown.pdf").c_str());
                stack.Write(); canvas.Write();
            }
        }
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        gSystem->Exit(1);
    }
}
