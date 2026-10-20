# run000400.dat全件テスト

入力は`~/GitHub/nestdaq-run/data/run000400.dat`の全43362 TF。
TFBFilePlayerへ`--max-iterations`を指定せず、`--wait 1`でEOFまで再生する。
テストスクリプトの既定も全件・wait=1に変更した。
入力レコードの境界と長さを検証してTF数を数え、FileSinkにだけ全件数の
`--max-iteration 43362`を渡す。PlayerはEOFでも受信完了ackまでソケットを維持し、
ByTrackの全受信記録とFileSinkの正常終了でキューの消化を確認する。

ByTrackの本来の`ProcessSlice()`は`#if 0`で囲み、`false`だけ返す実装を有効にした。
トラック再構成は行わず全SLICEを落とす。TrgTime/SLICE受信記録はこの処理より前で
行うため、LogicFilterの結果とスライスの対応の照合は可能。
FileSinkにはMETA等が届く。これはトラック処理やフィルター通過率の検証ではない。

## トリガー設定

`~/GitHub/nestdaq-run/LogicFilter/triggerlogic.sh`は実行/sourceせず、
同じ値をLogicFilterの起動引数で渡す。Redisは更新しない。
使用した完全な引数は各出力ディレクトリの`LogicFilter.command.json`へ保存する。

```sh
--trigger-signals '(0 0xc0a802a9 8 0) (1 0xc0a802a9 10 0) (2-0 0xc0a802aa 16 -12) (2-0 0xc0a802aa 17 -12) (2-1 0xc0a802aa 18 -12) (2-1 0xc0a802aa 19 -12) (2-2 0xc0a802aa 20 -12) (2-2 0xc0a802aa 21 -12) (2-3 0xc0a802aa 22 -12) (2-3 0xc0a802aa 23 -12) (2-4 0xc0a802aa 24 -12) (2-4 0xc0a802aa 25 -12) (2-5 0xc0a802aa 27 -12) (2-5 0xc0a802aa 28 -12)'
--trigger-expression '(0 & 1) & 2'
--trigger-width 30
```

`0xc0a802a9`は192.168.2.169、`0xc0a802aa`は192.168.2.170。
グループ2は6組のch対のどれかが成立したときに立つ。
各対は同じサブグループ（2-0～2-5）なので対内はAND、対同士はOR。
グループ0、1、2のすべてが成立するbinを判定し、そのrising edgeをTrgTimeにする。

## 再現コマンド

既存のビルドディレクトリを使う場合も、ByTrackとPlayerを再ビルドすること。

```sh
cmake --build build --target TFBFilePlayer FilterTimeFrameSliceByTrack LogicFilter -j4
cmake --build /tmp/nestdaq-trgtime-cpu --target TFBFilePlayer FilterTimeFrameSliceByTrack LogicFilter -j4
python3 tools/run_trgtime_pipeline.py --build-dir ./build \
  --run-dir "$HOME/GitHub/nestdaq-run" --output-dir /tmp/nestdaq-full-safe-gpu-control \
  --wait 1 --performance --timeout 1800
python3 tools/run_trgtime_pipeline.py --build-dir /tmp/nestdaq-trgtime-cpu \
  --run-dir "$HOME/GitHub/nestdaq-run" --output-dir /tmp/nestdaq-full-safe-cpu \
  --wait 1 --performance --timeout 1800
python3 tools/run_trgtime_pipeline.py --build-dir ./build \
  --run-dir "$HOME/GitHub/nestdaq-run" --output-dir /tmp/nestdaq-full-safe-gpu-details \
  --wait 1 --performance --gpu-exec-details --timeout 1800
python3 tools/compare_trgtime.py /tmp/nestdaq-full-safe-cpu/trgtime.csv /tmp/nestdaq-full-safe-gpu-control/trgtime.csv
python3 tools/compare_trgtime.py /tmp/nestdaq-full-safe-gpu-control/trgtime.csv /tmp/nestdaq-full-safe-gpu-details/trgtime.csv
```

性能集計では全件受信記録から最初の10 TFを除く。グラフは800×600、線幅3、
時間分布は200 bins、PDF/PNG/ROOTを出力する。
`Elapsed Time`はConditionalRunのReceive時間を除いた値。
今回のwait=1・トラック処理無効の条件は、以前のwait=0・トラック処理有効の
2000 TFテストと異なるので、両者の差をコードの改善率とは扱わない。

## 結果

CPU版・GPU版・GPU詳細計測版の3実行とも、全5デバイスが正常終了した。
受信記録は43362メッセージ、173448 META、338701 TrgTime、338701受信SLICE。
CPU対GPU、GPU対詳細計測GPUの両比較で、TF・時刻・type・順序を含め全件一致。
各runのPlayerコマンドで`--max-iterations`がなく`--wait 1`であることも検査した。

先頭10 TFを除く43352 TFを、トリガーあり21032 TF／なし22320 TFへ分割した。
詳細計測なしのCPU/GPUについて、TFあたり平均は以下のとおり。

| 条件 | 指標 | CPU [us/TF] | GPU [us/TF] | CPU/GPU |
|---|---|---:|---:|---:|
| トリガーあり | Exec | 2784.2 | 478.2 | 5.82 |
| トリガーなし | Exec | 2066.0 | 432.3 | 4.78 |
| トリガーあり | ConditionalRun−Receive | 4430.1 | 851.1 | 5.21 |
| トリガーなし | ConditionalRun−Receive | 3670.5 | 701.5 | 5.23 |

GPU詳細計測版のトリガーあり群はExec平均478.7 us/TF。
主な内訳は入力・記述子H2Dが274.5 us（57.3%）、Mark起動が8.7 us（1.8%）、
件数D2Hが121.4 us（25.4%）、ヒット時刻D2Hが17.6 us（3.7%）。
件数D2HはMarkとScanの完了待ちも含む。Mark起動回数4回/TF、
生データコピー156回/TF、記述子4992 bytes/TFを確認した。
区間和と入力・出力カウンターの整合検査も全受信TFで通過した。
詳細計測なしとの別run間の平均差は、トリガーあり+0.11%、なし+0.46%。

テスト中の一時点では、GPU LogicFilterのRSS約615 MiB、ByTrack約258 MiBだった。
これはピーク値ではない。ByTrackは元のトラック処理を無効化したままにしている。

### 出力

- CPUの生記録: `/tmp/nestdaq-full-safe-cpu/`
- GPUの生記録: `/tmp/nestdaq-full-safe-gpu-control/`
- GPU詳細計測の生記録: `/tmp/nestdaq-full-safe-gpu-details/`
- GPU Exec内訳: `/tmp/nestdaq-full-safe-exec-analysis/host_time.pdf`、`host_fraction.pdf`、`exec_stages.csv`
- CPU/GPU比較: `/tmp/nestdaq-full-safe-comparison/groups_summary.csv`
- トリガーありのElapsed Time対TDC数: `/tmp/nestdaq-full-safe-comparison/trigger_present/elapsed_no_receive_vs_tdc.pdf`
- トリガー有無のElapsed Time比較: `/tmp/nestdaq-full-safe-comparison/elapsed_no_receive_trigger_groups.pdf`
- ConditionalRun内訳: `/tmp/nestdaq-full-safe-breakdown/CPU_time_breakdown.pdf`、`GPU_time_breakdown.pdf`

グラフの再生成コマンド（再測定不要）:

```sh
root -l -b -q 'tools/analyze_gpu_exec.C("/tmp/nestdaq-full-safe-gpu-details/logicfilter.root","/tmp/nestdaq-full-safe-gpu-control/logicfilter.root","/tmp/nestdaq-full-safe-exec-analysis",10)'
root -l -b -q 'tools/compare_logicfilter_trigger_groups.C("/tmp/nestdaq-full-safe-cpu/logicfilter.root","/tmp/nestdaq-full-safe-gpu-control/logicfilter.root","/tmp/nestdaq-full-safe-comparison",10)'
root -l -b -q 'tools/summarize_logicfilter_breakdown.C("/tmp/nestdaq-full-safe-cpu/logicfilter.root","/tmp/nestdaq-full-safe-gpu-control/logicfilter.root","/tmp/nestdaq-full-safe-breakdown",10)'
```

TF数事前検査のテストを追加し、途中切れ・未知magic・不正ヘッダー長・Filterヘッダー
付き入力を検証した。Pythonテスト計14件、CPU/GPUビルド、`git diff --check`は通過。
