#include <fstream>
#include <iostream>
#include <vector>
#include <memory>
#include "FilterTimeFrameSliceBySomething.h"
#include "FilterTimeFrameSliceABC.icxx"
#include "fairmq/runDevice.h"
#include "utility/MessageUtil.h"
#include "UnpackTdc.h"
#include "SubTimeFrameHeader.h"
#include "TimeFrameHeader.h"

// 書き出しファイル名用
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>

#define DEBUG 0

using nestdaq::FilterTimeFrameSliceBySomething;
namespace bpo = boost::program_options;

FilterTimeFrameSliceBySomething::FilterTimeFrameSliceBySomething()
{
    // 書き出しファイル名の設定
    /*
    std::string outputFilename_;
    public:
        FilterTimeFrameSliceBySomething(){
            auto now = std::chrono::system_clock::now();
            std::time_t t = std::chrono::system_clock::to_time_t(now);
            std::tm tm = *std::localtime(&t);

            std::ostringstream filename;
            filename << "/home/nestdaq/kashima/run/output_FilterTimeFrameSliceBySomething_"
                    << std::put_time(&tm, "%Y%m%d_%H%M%S")
                    << ".txt";
            outputFilename_ = filename.str();
        }
    */
}

bool FilterTimeFrameSliceBySomething::ProcessSlice(TTF& tf)
{
    auto t_start = std::chrono::high_resolution_clock::now();
    double elapsed_ms = 0; // このsliceの処理にかかった時間[ms]
    // elapsed_ms += std::chrono::duration<double, std::milli>(t_end - t_start).count();
    #if 0

    int numSTF = tf.size();
    std::cout << "Number of SubTimeFrames: " << numSTF << std::endl;
    for(int istf = 0; istf < numSTF; ++istf){
        auto& stf = *(tf[istf]);
        auto header = stf.GetHeader();
        std::cout << " SubTimeFrame " << istf
                  << " FEM ID: " << std::hex << header->femId
                  << " Type: " << std::dec << static_cast<int>(header->femType)
                  << " Num HB: " << stf.size()
                  << std::endl;
    }
    #endif
    //std::cout << "関数内部" << std::endl;

    // std::ofstream outFile(outputFilename_, std::ios::app);


    std::vector<std::unique_ptr<int>> t1tof_r_hits;
    std::vector<std::unique_ptr<int>> t1tof_l_hits;
    std::vector<std::unique_ptr<int>> utof_r_hits;
    std::vector<std::unique_ptr<int>> utof_l_hits;

    bool tofCondition = false;

    auto tfHeader = tf.GetHeader();
    //std::cout << "TimeFrameヘッダーサイズ: " << sizeof(*tfHeader) << " bytes" << std::endl;

    for (auto& SubTimeFrame : tf) {
        //std::cout << "SubTimeFrame内" << std::endl;
        auto header = SubTimeFrame->GetHeader();
        auto& hbf = SubTimeFrame->at(0); // SubTimeFrameの上から１つ目のHB

        //std::cout << "SubTimeFrameヘッダーサイズ: " << sizeof(*header) << " bytes" << std::endl;
        uint64_t nData = hbf->GetNumData();

        for (int i = 0; i < nData; ++i) {
            if (header->femType == SubTimeFrame::TDC64H) {
                TDC64H::tdc64 tdc;
                TDC64H::Unpack(hbf->UncheckedAt(i), &tdc);

                if (header->femId == 0xc0a802aa) {
                    //std::cout << "64H_utofv1: " << tdc.tdc << std::endl;
                    //std::cout << "64H_t1v1: " << tdc.tdc << std::endl;
                }
            } else if (header->femType == SubTimeFrame::TDC64L) {
                TDC64L::tdc64 tdc;
                TDC64L::Unpack(hbf->UncheckedAt(i), &tdc);
                if (header->femId == 0xc0802aa9 && (tdc.ch == 10 || tdc.ch == 12)) {
                    //std::cout << "64L_utof: " << tdc.tdc << std::endl;
                    //std::cout << "64L_t1: " << tdc.tdc << std::endl;
                }
            } else if (header->femType == SubTimeFrame::TDC64H_V3) {
                TDC64H_V3::tdc64 tdc;
                TDC64H_V3::Unpack(hbf->UncheckedAt(i), &tdc);
                if (header->femId == 0xc0a802aa) {
                    if (tdc.ch == 12) {
                        // TOF
                        t1tof_r_hits.push_back(std::make_unique<int>(tdc.tdc));
                    } else if (tdc.ch == 10) {
                        // TOF
                        t1tof_l_hits.push_back(std::make_unique<int>(tdc.tdc));
                    }
                } else if (header->femId == 0xc0a802a9) {
                    if (tdc.ch == 10) {
                        // TOF
                        utof_r_hits.push_back(std::make_unique<int>(tdc.tdc));
                    } else if (tdc.ch == 8) {
                        // TOF
                        utof_l_hits.push_back(std::make_unique<int>(tdc.tdc));
                    }
                }
            } else if (header->femType == SubTimeFrame::TDC64L_V3) {
                TDC64L_V3::tdc64 tdc;
                TDC64L_V3::Unpack(hbf->UncheckedAt(i), &tdc);
                // Add processing for TDC64L_V3 if needed
            }

            /* マルチプリシティ前段処理
            if (wireMap.find(header->femId) != wireMap.end()) {
                Flt::Wire_map& wireMapEntry = wireMap[header->femId];
                GeoIDs[header->femId].push_back(wireMapEntry.geo);
                std::cout << "Wire ID: " << wireMapEntry.id << std::endl;
            } */
        }
    }

    // デバッグ用
/*     flt->PrintList(t1tof_r_hits, "t1tof_r_hits");
    flt->PrintList(t1tof_l_hits, "t1tof_l_hits");
    flt->PrintList(utof_r_hits, "utof_r_hits");
    flt->PrintList(utof_l_hits, "utof_l_hits"); */

    // 面ごとの各ペアの平均値の計算と表示
    if ((!t1tof_r_hits.empty() || !t1tof_l_hits.empty()) && (!utof_r_hits.empty() || !utof_l_hits.empty())) {


        //auto t1_averages = flt->CalculateAveragePairs(t1tof_r_hits, t1tof_l_hits);
        //auto utof_averages = flt->CalculateAveragePairs(utof_r_hits, utof_l_hits);

        // アルゴリズム改善
        auto t1_averages = flt->CalculateAveragePairs_limited(t1tof_r_hits, t1tof_l_hits);
        auto utof_averages = flt->CalculateAveragePairs_limited(utof_r_hits, utof_l_hits);


        // std::ofstream outFile("/home/nestdaq/kashima/run/output.txt", std::ios::app);
        



        /*
        outFile << "T1_hits L: " << std::dec <<t1tof_l_hits.size() << " R: " << std::dec<< t1tof_r_hits.size() << " avg_pair: " << std::dec << t1_averages.size() << std::endl;
        outFile << "UTOF_hits L: " << std::dec <<utof_l_hits.size() << " R: " << std::dec<< utof_r_hits.size() << "avg_pair: " << std::dec << utof_averages.size()  << std::endl;
        outFile.close();
        */





        /*
        auto t1_L_max = *std::max_element(t1tof_l_hits.begin(), t1tof_l_hits.end());
        auto t1_L_min = *std::min_element(t1tof_l_hits.begin(), t1tof_l_hits.end());

        outFile << "T1_hits L: " << t1tof_l_hits.size()
                << " R: " << t1tof_r_hits.size()
                << " avg_pair: " << t1_averages.size()
                << " maxtdc-mintdc(means_range): " << (t1_L_max - t1_L_min)
                << std::endl;

        auto utof_L_max = *std::max_element(utof_l_hits.begin(), utof_l_hits.end());
        auto utof_L_min = *std::min_element(utof_l_hits.begin(), utof_l_hits.end());

        outFile << "UTOF_hits L: " << utof_l_hits.size()
                << " R: " << utof_r_hits.size()
                << " avg_pair: " << utof_averages.size()
                << " maxtdc-mintdc(means_range): " << (utof_L_max - utof_L_min)
                << std::endl;
        */


        // TDC幅の出力用
        /*
        if (!t1tof_l_hits.empty()) {
            auto it_max = std::max_element(t1tof_l_hits.begin(), t1tof_l_hits.end(),
                [](const std::unique_ptr<int>& a, const std::unique_ptr<int>& b){
                    return *a < *b;
                });
            auto it_min = std::min_element(t1tof_l_hits.begin(), t1tof_l_hits.end(),
                [](const std::unique_ptr<int>& a, const std::unique_ptr<int>& b){
                    return *a < *b;
                });

            int t1_L_max = **it_max;
            int t1_L_min = **it_min;

            outFile << "T1_hits L: " << t1tof_l_hits.size()
                    << " R: " << t1tof_r_hits.size()
                    << " avg_pair: " << t1_averages.size()
                    << " maxtdc-mintdc(means_range): " << (t1_L_max - t1_L_min)
                    << std::endl;
        }

        if (!utof_l_hits.empty()) {
            auto it_max = std::max_element(utof_l_hits.begin(), utof_l_hits.end(),
                [](const std::unique_ptr<int>& a, const std::unique_ptr<int>& b){
                    return *a < *b;
                });
            auto it_min = std::min_element(utof_l_hits.begin(), utof_l_hits.end(),
                [](const std::unique_ptr<int>& a, const std::unique_ptr<int>& b){
                    return *a < *b;
                });

            int utof_L_max = **it_max;
            int utof_L_min = **it_min;

            outFile << "UTOF_hits L: " << utof_l_hits.size()
                    << " R: " << utof_r_hits.size()
                    << " avg_pair: " << utof_averages.size()
                    << " maxtdc-mintdc(means_range): " << (utof_L_max - utof_L_min)
                    << std::endl;
        }
        outFile.close();
        */
        

        flt->CalculateAndPrintTOF(utof_averages, t1_averages);

        for (const auto& t1_avg : t1_averages) {
            for (const auto& utof_avg : utof_averages) {
                //最小値 最大値の順番
                if (flt->CheckAllTOFConditions(*t1_avg, *utof_avg, t1tof_r_hits, t1tof_l_hits, utof_r_hits, utof_l_hits, -130000, -125000)) {
                    tofCondition = true;
                    //std::cout << "合格tof: " << *utof_avg - *t1_avg << std::endl;
                    //std::cout << "min: -130000" << " " << "max: -125000" << std::endl;
                    break;
                }
           }
            if (tofCondition) break;
        }
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    elapsed_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();
    #if 0
    std::ofstream outFile("/home/nestdaq/kashima/run/output.txt", std::ios::app);
    outFile << "this_slice_size: " << tf.GetRealLength()
            << " elapsed_ms_per_slice: " << elapsed_ms
            << " hits_T1: " << (t1tof_r_hits.size() + t1tof_l_hits.size())/2.0
            << " hits_UTOF: " << (utof_r_hits.size() + utof_l_hits.size())/2.0
            << std::endl;
    outFile.close();
    #endif

    if (tofCondition) {
        /*// TOF フィルターデバッグ用 TDC値分布を見る用
        std::ofstream outFile("tof_output.txt", std::ios::app);
        double tof = *utof_ave - *t1tof_ave;
        outFile << *utof_ave << " " << *t1tof_ave << " " << tof << "\n";
        outFile.close();*/

        //std::cout << "関数内tofbool値: true" << std::endl;
        return true;
    }

    return false;
}

void addCustomOptions(bpo::options_description& options)
{
    using opt = FilterTimeFrameSliceBySomething::OptionKey;

    options.add_options()
        (opt::InputChannelName.data(),
         bpo::value<std::string>()->default_value("in"),
         "Name of the input channel")
        (opt::OutputChannelName.data(),
         bpo::value<std::string>()->default_value("out"),
         "Name of the output channel")
        (opt::DQMChannelName.data(),
         bpo::value<std::string>()->default_value("dqm"),
         "Name of the data quality monitoring channel")
        (opt::PollTimeout.data(),
         bpo::value<std::string>()->default_value("1"),
         "Timeout of polling (in msec)")
        (opt::SplitMethod.data(),
         bpo::value<std::string>()->default_value("1"),
         "STF split method");
}

std::unique_ptr<fair::mq::Device> getDevice(fair::mq::ProgOptions& /*config*/)
{
    return std::make_unique<FilterTimeFrameSliceBySomething>();
}
