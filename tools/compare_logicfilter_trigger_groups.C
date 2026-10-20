// No new DAQ run: split the recorded per-ConditionalRun data by triggers >= 1 / = 0.
#include "compare_logicfilter_timing.C"

void compare_logicfilter_trigger_groups(const char* cpuPath, const char* gpuPath,
    const char* outputDir = "/tmp/logicfilter-trigger-groups", int warmup = 10) {
    using namespace lf_timing;
    try {
        const std::string prefix = std::string(outputDir) + "/";
        compare_logicfilter_timing(cpuPath, gpuPath, (prefix + "trigger_present").c_str(), warmup, 1);
        compare_logicfilter_timing(cpuPath, gpuPath, (prefix + "trigger_absent").c_str(), warmup, 0);
        auto cpu = read(cpuPath, false), gpu = read(gpuPath, true);
        std::array<std::vector<size_t>, 2> selected;
        for (size_t i = warmup; i < cpu.rows.size(); ++i) selected[cpu.rows[i].data[10] > 0].push_back(i);
        if (selected[0].size() + selected[1].size() != cpu.rows.size() - warmup)
            throw std::runtime_error("Trigger partition does not cover all TFs");
        std::ofstream summary(prefix + "groups_summary.csv");
        summary << "group,tf_count,metric,cpu_mean_us,gpu_mean_us,cpu_median_us,gpu_median_us,cpu_p95_us,gpu_p95_us,cpu_over_gpu\n";
        std::ofstream dataSummary(prefix + "groups_data.csv");
        dataSummary << "group,tf_count,mean_input_bytes,mean_tdc_tot_words_type_0x0b,mean_triggers,total_triggers\n";
        if (!summary || !dataSummary) throw std::runtime_error("Cannot write group summaries");
        auto mean = [](const std::vector<double>& v) { return std::accumulate(v.begin(), v.end(), 0.0) / v.size(); };
        for (int group : {1, 0}) {
            const char* label = group ? "triggers_ge_1" : "triggers_eq_0";
            const auto& indices = selected[group];
            double bytes = 0, words = 0;
            ULong64_t totalTriggers = 0;
            for (auto i : indices) {
                bytes += cpu.rows[i].data[1]; words += cpu.rows[i].data[7]; totalTriggers += cpu.rows[i].data[10];
            }
            dataSummary << label << ',' << indices.size() << ',' << bytes / indices.size() << ','
                        << words / indices.size() << ',' << double(totalTriggers) / indices.size() << ',' << totalTriggers << '\n';
            for (size_t metric = 0; metric < times.size(); ++metric) {
                std::array<std::vector<double>, 2> values;
                for (auto i : indices) {
                    values[0].push_back(cpu.rows[i].time[metric]); values[1].push_back(gpu.rows[i].time[metric]);
                }
                const double a = mean(values[0]), b = mean(values[1]);
                summary << label << ',' << indices.size() << ',' << times[metric] << ',' << a << ',' << b << ','
                        << quantile(values[0], .5) << ',' << quantile(values[1], .5) << ','
                        << quantile(values[0], .95) << ',' << quantile(values[1], .95) << ',' << (b > 0 ? a / b : 0) << '\n';
            }
        }
        TFile output((prefix + "trigger_groups.root").c_str(), "RECREATE");
        TNamed("selection_note", "Trigger presence is a proxy, not a measured beam spill gate").Write();
        TParameter<int>("warmup_tfs", warmup).Write();
        TParameter<Long64_t>("trigger_present_tfs", selected[1].size()).Write();
        TParameter<Long64_t>("trigger_absent_tfs", selected[0].size()).Write();
        for (bool triggerOnly : {false, true}) {
            const size_t metric = triggerOnly ? 2 : 9;
            const std::string stem = triggerOnly ? "trigger" : "elapsed_no_receive";
            const std::string yLabel = triggerOnly ? "Trigger processing elapsed time [us/TF]" : "Elapsed time excluding Receive [us/TF]";
            double maximum = 0;
            for (size_t i = warmup; i < cpu.rows.size(); ++i)
                maximum = std::max({maximum, cpu.rows[i].time[metric], gpu.rows[i].time[metric]});
            TCanvas canvas((stem + "_trigger_groups").c_str(), "Comparison by trigger presence", 800, 600);
            canvas.SetLeftMargin(.14);
            std::vector<std::unique_ptr<TH1D>> histograms;
            double ymax = 0;
            for (int variant : {0, 1}) {
                for (int group : {1, 0}) {
                    const std::string label = std::string(variant ? "GPU" : "CPU") + (group ? ", triggers >= 1" : ", triggers = 0");
                    auto h = std::make_unique<TH1D>((stem + "_" + std::to_string(variant) + "_" + std::to_string(group)).c_str(),
                        ("Elapsed-time comparison by trigger presence;" + yLabel + ";Fraction of TFs / bin").c_str(),
                        200, 0, maximum * 1.01 + 1);
                    h->SetDirectory(nullptr);
                    for (auto i : selected[group]) h->Fill(variant ? gpu.rows[i].time[metric] : cpu.rows[i].time[metric]);
                    h->Scale(1.0 / h->Integral());
                    h->SetLineColor(variant ? kRed : kBlue);
                    h->SetLineStyle(group ? 1 : 2);
                    h->SetLineWidth(3);
                    h->SetStats(false);
                    h->SetTitle(("Elapsed-time comparison by trigger presence;" + yLabel + ";Fraction of TFs / bin").c_str());
                    ymax = std::max(ymax, h->GetMaximum());
                    histograms.push_back(std::move(h));
                }
            }
            TLegend legend(.51, .68, .88, .88);
            size_t index = 0;
            for (int variant : {0, 1}) {
                for (int group : {1, 0}) {
                    auto* h = histograms[index].get();
                    h->SetMaximum(ymax * 1.25);
                    h->Draw(index ? "HIST SAME" : "HIST");
                    const auto label = std::string(variant ? "GPU" : "CPU") + (group ? ", triggers >= 1" : ", triggers = 0")
                        + " (n=" + std::to_string(selected[group].size()) + ")";
                    legend.AddEntry(h, label.c_str(), "l");
                    h->Write();
                    ++index;
                }
            }
            legend.Draw();
            canvas.SaveAs((prefix + stem + "_trigger_groups.png").c_str());
            canvas.SaveAs((prefix + stem + "_trigger_groups.pdf").c_str());
            canvas.Write();
        }
        // Acquisition-order view: do not infer physical beam time from replay time.
        TGraph nTriggers, tdcWords;
        nTriggers.SetName("triggers_by_received_tf"); tdcWords.SetName("tdc_words_by_received_tf");
        for (size_t i = warmup; i < cpu.rows.size(); ++i) {
            nTriggers.SetPoint(nTriggers.GetN(), i, cpu.rows[i].data[10]);
            tdcWords.SetPoint(tdcWords.GetN(), i, cpu.rows[i].data[7]);
        }
        nTriggers.SetLineWidth(3); tdcWords.SetLineWidth(3);
        nTriggers.SetLineColor(kBlue); tdcWords.SetLineColor(kBlue);
        TCanvas sequence("acquisition_order", "Trigger and TDC activity in acquisition order", 800, 600);
        sequence.Divide(1, 2);
        sequence.cd(1);
        nTriggers.SetTitle("Trigger activity (not a measured spill gate);Received TF index;Triggers / TF");
        nTriggers.Draw("AL");
        sequence.cd(2);
        tdcWords.SetTitle("TDC activity;Received TF index;TDC+TOT words (type 0x0b) / TF");
        tdcWords.Draw("AL");
        sequence.SaveAs((prefix + "activity_by_tf.png").c_str());
        sequence.SaveAs((prefix + "activity_by_tf.pdf").c_str());
        nTriggers.Write(); tdcWords.Write(); sequence.Write();
        std::cout << "PARTITION: triggers >= 1: " << selected[1].size() << ", triggers = 0: " << selected[0].size()
                  << ", total after warmup: " << (cpu.rows.size() - warmup) << '\n';
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        gSystem->Exit(1);
    }
}
