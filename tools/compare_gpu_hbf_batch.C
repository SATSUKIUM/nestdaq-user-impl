#include "compare_logicfilter_timing.C"

void compare_gpu_hbf_batch(const char* sequentialPath, const char* batchedPath,
    const char* outputDir = "/tmp/hbf-batch-comparison", int warmup = 10) {
    using namespace lf_timing;
    auto sequential = read(sequentialPath, true), batched = read(batchedPath, true);
    if (sequential.rows.size() != batched.rows.size() || sequential.signals != batched.signals ||
        sequential.expression != batched.expression || sequential.width != batched.width)
        throw std::runtime_error("Sequential/batched run configurations differ");
    for (size_t i = 0; i < sequential.rows.size(); ++i)
        if (sequential.rows[i].data != batched.rows[i].data)
            throw std::runtime_error("Sequential/batched TF data differ");
    if (warmup < 0 || size_t(warmup) >= sequential.rows.size()) throw std::runtime_error("Invalid warmup");
    const std::array<const char*, 2> paths = {{sequentialPath, batchedPath}};
    for (size_t index = 0; index < paths.size(); ++index) {
        TFile file(paths[index]);
        auto* mode = file.Get<TParameter<int>>("hbf_batch_size");
        if (!mode || mode->GetVal() != (index == 0 ? 1 : 4))
            throw std::runtime_error("Expected batch-size 1/4 metadata");
    }
    std::array<std::vector<size_t>, 3> groups;
    for (size_t i = warmup; i < sequential.rows.size(); ++i) {
        groups[sequential.rows[i].data[10] > 0 ? 0 : 1].push_back(i);
        groups[2].push_back(i);
    }
    std::filesystem::create_directories(outputDir);
    const std::string prefix = std::string(outputDir) + "/";
    std::ofstream csv(prefix + "summary.csv");
    csv << "group,tf_count,metric,sequential_mean_us,batched_mean_us,reduction_percent,sequential_median_us,batched_median_us,sequential_p95_us,batched_p95_us\n";
    const std::array<const char*, 3> names = {{"triggers_ge_1", "triggers_eq_0", "all"}};
    const std::array<const char*, 3> labels = {{"triggers #geq 1", "triggers = 0", "All TFs"}};
    const std::array<size_t, 4> metrics = {{4, 3, 2, 9}};
    TFile output((prefix + "comparison.root").c_str(), "RECREATE");
    TNamed("note", "Same GPU binary and input; batch-size 1 versus 4; Receive excluded from elapsed time; Send included.").Write();
    TParameter<int>("warmup_tfs", warmup).Write();
    for (size_t metric : metrics) {
        const std::string name = times[metric];
        TH1D before((name + "_before").c_str(), "", 3, 0, 3);
        TH1D after((name + "_after").c_str(), "", 3, 0, 3);
        for (size_t group = 0; group < groups.size(); ++group) {
            if (groups[group].empty()) throw std::runtime_error("Empty trigger group");
            std::vector<double> a, b;
            for (size_t i : groups[group]) { a.push_back(sequential.rows[i].time[metric]); b.push_back(batched.rows[i].time[metric]); }
            const double ma = std::accumulate(a.begin(), a.end(), 0.) / a.size();
            const double mb = std::accumulate(b.begin(), b.end(), 0.) / b.size();
            csv << names[group] << ',' << a.size() << ',' << name << ',' << ma << ',' << mb << ','
                << 100 * (ma - mb) / ma << ',' << quantile(a, .5) << ',' << quantile(b, .5) << ','
                << quantile(a, .95) << ',' << quantile(b, .95) << '\n';
            before.SetBinContent(group + 1, ma); after.SetBinContent(group + 1, mb);
            before.GetXaxis()->SetBinLabel(group + 1, labels[group]);
            after.GetXaxis()->SetBinLabel(group + 1, labels[group]);
            std::cout << names[group] << ' ' << name << ": " << ma << " -> " << mb
                      << " us (" << 100 * (ma - mb) / ma << "% reduction)\n";
        }
        const std::string title = metric == 9 ? "ConditionalRun elapsed time excluding Receive" :
            metric == 4 ? "GPU Exec time" : metric == 3 ? "GPU cleanup host time" : "GPU trigger time";
        before.SetTitle((title + ";;Mean elapsed time [#mus/TF]").c_str());
        before.SetStats(false); after.SetStats(false);
        before.SetLineWidth(3); after.SetLineWidth(3);
        before.SetLineColor(kAzure + 2); after.SetLineColor(kOrange + 7);
        before.SetFillColor(kAzure - 9); after.SetFillColor(kOrange - 2);
        before.SetBarWidth(.35); after.SetBarWidth(.35);
        before.SetBarOffset(.1); after.SetBarOffset(.5);
        before.SetMinimum(0); before.SetMaximum(1.4 * std::max(before.GetMaximum(), after.GetMaximum()));
        TCanvas canvas((name + "_canvas").c_str(), title.c_str(), 800, 600);
        canvas.SetLeftMargin(.13); canvas.SetBottomMargin(.13);
        before.Draw("BAR"); after.Draw("BAR SAME");
        TLegend legend(.40, .76, .88, .88); legend.SetBorderSize(0);
        legend.AddEntry(&before, "1 HBF per batch", "f"); legend.AddEntry(&after, "4 HBFs per batch", "f"); legend.Draw();
        canvas.SaveAs((prefix + name + ".pdf").c_str()); canvas.SaveAs((prefix + name + ".png").c_str());
        before.Write(); after.Write(); canvas.Write();
        if (metric == 4 || metric == 9) {
            double high = 0;
            for (size_t i : groups[2]) high = std::max({high, sequential.rows[i].time[metric], batched.rows[i].time[metric]});
            TH1D a((name + "_distribution_before").c_str(), (title + ";Elapsed time [#mus/TF];TFs / bin").c_str(), 200, 0, high * 1.02 + 1);
            TH1D b((name + "_distribution_after").c_str(), "", 200, 0, high * 1.02 + 1);
            for (size_t i : groups[2]) { a.Fill(sequential.rows[i].time[metric]); b.Fill(batched.rows[i].time[metric]); }
            a.SetStats(false); b.SetStats(false); a.SetLineWidth(3); b.SetLineWidth(3);
            a.SetLineColor(kAzure + 2); b.SetLineColor(kOrange + 7);
            TCanvas distribution((name + "_distribution").c_str(), title.c_str(), 800, 600);
            distribution.SetLeftMargin(.13); distribution.SetBottomMargin(.13); distribution.SetLogy();
            a.SetMinimum(.5); a.SetMaximum(5 * std::max(a.GetMaximum(), b.GetMaximum()));
            a.Draw("HIST"); b.Draw("HIST SAME");
            TLegend l(.50,.76,.88,.88); l.SetBorderSize(0);
            l.AddEntry(&a,"1 HBF per batch","l"); l.AddEntry(&b,"4 HBFs per batch","l"); l.Draw();
            distribution.SaveAs((prefix + name + "_distribution.pdf").c_str());
            distribution.SaveAs((prefix + name + "_distribution.png").c_str());
            a.Write(); b.Write(); distribution.Write();
        }
    }
}
