# 4 HBFをまとめたGPUトリガー処理

## ブランチとコミット

作業前に`LogicFilterMore32GPU`の変更を次の2件へ保存した。

- `b5a73fd`: quietなConditionalRun計測、STOP時のテキスト保存、分布作成マクロ。
- `d873419`: ByTrackのトラック処理無効化とPlayerログ削減。

その後`verify-logicfilter-cpu-gpu-trgtime`へcheckoutし、`feceedf`で4 HBF一括処理を実装した。
このブランチはROOT詳細計測・TrgTime照合用の検証ブランチである。
通常運用のインストール済みバイナリと`bin/LogicFilter`は更新していない。
検証時は`build/LogicFilter/GPU/LogicFilter`を直接使う。

## 領域の配置

各HBFの論理TimeRegion長は131072 bins（4 ns/bin）を維持し、GPUの確保領域を4 HBF分にする。
合計524288 binsなのでCUBのint32要素数にも収まる。
`PrepareHbfBatch()`は総bin数がINT_MAXを超える設定を拒否する。

```text
MainTimeRegion: [HBF0のN bins][HBF1のN bins][HBF2のN bins][HBF3のN bins]
SubTimeRegion : [HBF0: S領域][HBF1: S領域][HBF2: S領域][HBF3: S領域]
各HBFのS領域: [sub0のN bins][sub1のN bins] ... [sub(S-1)のN bins]
```

N=131072、S=6の場合、主なGPU配列は以下のとおり。単位はMiB。

| 配列 | 1 HBF | 4 HBF |
|---|---:|---:|
| MainTimeRegion | 0.5 | 2 |
| SubTimeRegion合計 | 3 | 12 |
| edge flags | 0.5 | 2 |
| exclusive-scan offsets | 0.5 | 2 |
| 最大ヒット時刻配列 | 0.25 | 1 |
| 合計 | 4.75 | 19 |

上記に入力・記述子・CUB作業領域・LUT等が加わる。
容量は必要なときだけ拡張し、使用HBF数が減っても確保済み領域を再利用する。
CUBの作業領域は使用中のbin数で問い合わせ、最大の確保量を保持する。

## Markと境界

`MarkSegment`に`hbfIndex`を追加した。生データは従来どおり各segからGPUへ直接コピーし、
ホスト側の詰め直しは行わない。全HBFの記述子をまとめてMarkを1回起動する。
各blockは担当segを二分探索し、`hbfIndex`からそのHBFのMain/SubTimeRegionを選ぶ。

`markWord()`へ渡す領域ポインタはHBFの先頭、サイズは1 HBF分のNのままにしている。
このため遅延や時間窓がHBF境界を越えて隣のHBFへ書き込むことはない。
サブグループの初期化・集約もHBF番号とサブグループ番号を別々に求めるため、
別HBFに入った信号同士でANDが成立することはない。

## exclusive scanとcompaction

HBFごとのsegmented scanは使わず、全HBFのflag配列に通常のexclusive scanを1回かける。

1. LUT判定のrising-edge flagは各HBF先頭binで必ず0とする（CPUのbin 0除外規則）。
2. 残りのbinは、そのHBF内の直前binだけと比較する。
3. 全flagにCUB ExclusiveSumをかけ、全体の出力位置を求める。
4. compactionは全体内のbin番号をuint32_tで書き込む。
5. 全体の件数4 bytesを1回D2Hコピーし、その件数の時刻だけを1回D2Hコピーする。
6. ホストで`HBF番号 = globalBin / N`、`HBF内時刻 = globalBin % N`へ戻す。

prefix sumをHBF境界でリセットしなくても、出力配列全体の位置として使うので問題ない。
結果はglobalBinの昇順であり、HBF順と各HBF内の時刻順が保たれる。
ヒットなしHBFの出力ブロックも空のまま残し、META/TrgTimeとの対応を維持する。
各HBFのrising edge数上限はfloor(N/2)で、確保量はHBF数×floor(N/2)。
奇数Nの場合も単純なfloor(全bin数/2)ではなく、HBFごとの上限を合計している。

全カーネルとコピーは同じデフォルトストリームへ投入する。
MarkごとのcudaDeviceSynchronizeは不要で、最後の同期版D2Hコピーで完成を待つ。

## 比較用オプション

```text
--hbf-batch-size 1   HBFを1つずつ処理（比較基準）
--hbf-batch-size 4   最大4 HBFを一括処理（既定）
```

1～4を指定できる。HBFが不足する最後のバッチも実個数で処理する。
検証runner側は`--gpu-hbf-batch-size 1` / `4`で渡す。
ROOTには`hbf_batch_size`を保存した。
`gpu_exec_scan_calls`を追加し、件数D2H量を`4 bytes × scan_calls`で検査する。
従来の`hbf_count`は処理したHBF数、`mark_calls`は入力seg数で、意味は変更しない。

## GPU回帰テスト

`LogicFilter/tests/GpuHbfBatchTest.cu`で以下を確認した。

- HBF末尾の窓・次HBF先頭の成立・bin 0除外・空HBF。
- 4 TDC形式、IP/遅延/データ長の混在、39 seg/HBF、複数blockと端数block。
- 1 HBFずつ実行した結果と、4 HBF一括結果の一致。
- 4→2→1→3→4 HBFと変化したときの領域・CUB作業領域の再利用。
- HBFをまたいだサブグループANDの禁止と、同じHBF内のAND成立。
- N=1、5、64で最大rising-edge密度と出力容量。
- 空バッチ領域の古いヒットの消去と、不正サイズの拒否。
- 詳細なし／ホスト計測／CUDAイベント計測で同じ結果。

```sh
/usr/local/cuda/bin/nvcc -std=c++17 -arch=sm_86 -rdc=true -I. -ILogicFilter/GPU \
  -I"$HOME/local/nestdaq/include" LogicFilter/tests/GpuHbfBatchTest.cu \
  "$HOME/local/nestdaq/lib/libcudamap.a" -o /tmp/gpu-hbf-batch-test
/tmp/gpu-hbf-batch-test
compute-sanitizer --tool memcheck --error-exitcode 1 /tmp/gpu-hbf-batch-test
```

GPUテストは通過、compute-sanitizerはERROR SUMMARY: 0 errors。
既存GpuExecProfilerTest、Python 14 tests、CPU/GPU LogicFilterビルドも通過した。

## run000400.dat全件の結果

同じGPUバイナリで1 HBFずつ／4 HBF一括／4 HBF一括＋ホスト詳細計測を順番に実行した。
全件43362 TF、173448 META、338701 TrgTime、338701受信SLICE。
3 runとも全5デバイスが正常終了し、TrgTimeのTF・時刻・type・順序がすべて一致。
1 HBFずつの結果は、前回保存したCPUの全件結果とも一致した。

既存のRedis/tmuxへ接続しない専用IPCの検証パイプラインを使用。
ByTrackのProcessSliceはfalseを返す構成、TFBFilePlayerはwait=1、max-iterations指定なし。
最初の10 TFを除く43352 TFを集計（トリガーあり21032、なし22320）。

| 条件 | 指標 | 1 HBFずつ [us/TF] | 4 HBF一括 [us/TF] | 短縮率 |
|---|---|---:|---:|---:|
| トリガーあり | Exec | 537.2 | 436.9 | 18.7% |
| トリガーなし | Exec | 496.3 | 399.7 | 19.5% |
| 全体 | Exec | 516.1 | 417.8 | 19.1% |
| トリガーあり | Cleanup＋Exec | 665.0 | 501.9 | 24.5% |
| トリガーなし | Cleanup＋Exec | 617.3 | 465.5 | 24.6% |
| トリガーあり | ConditionalRun−Receive | 977.8 | 813.8 | 16.8% |
| トリガーなし | ConditionalRun−Receive | 817.2 | 663.6 | 18.8% |
| 全体 | ConditionalRun−Receive | 895.1 | 736.5 | 17.7% |

各条件1回の別run間比較なので実行時のばらつきを含む。
以前の別バイナリ・別runの平均値を今回の改善率の基準にはしていない。

詳細記録から、4 HBF一括では次の値を確認した。

| 指標 | TFあたり |
|---|---:|
| 入力seg数 | 156 |
| Markカーネル起動 | 1 |
| Scan一式 | 1 |
| 件数D2H | 4 bytes |
| 記述子H2D | 4992 bytes |

トリガーあり群の詳細計測Exec平均439.9 usのうち、入力＋記述子H2Dは284.8 us、
件数D2Hは123.9 us（Mark/scanの完了待ちも含む）。生データコピー156回/TFは残る。
ホスト詳細計測あり／なしのExec平均差は、トリガーあり+0.67%、全体+0.21%。

## 出力と再現

- 1 HBFずつ: `/tmp/nestdaq-hbf1-full-control/`
- 4 HBF一括: `/tmp/nestdaq-hbf4-full-control/`
- 4 HBF一括＋詳細: `/tmp/nestdaq-hbf4-full-details/`
- 比較図・CSV・ROOT: `test/hbf-batch-comparison/`
- 4 HBF一括の内部内訳: `test/hbf-batch-exec/`

```sh
cmake --build build --target LogicFilter TFBFilePlayer FilterTimeFrameSliceByTrack -j4
python3 tools/run_trgtime_pipeline.py --build-dir ./build --run-dir "$HOME/GitHub/nestdaq-run" \
  --output-dir /tmp/nestdaq-hbf1-full-control --wait 1 --performance --gpu-hbf-batch-size 1 --timeout 1800
python3 tools/run_trgtime_pipeline.py --build-dir ./build --run-dir "$HOME/GitHub/nestdaq-run" \
  --output-dir /tmp/nestdaq-hbf4-full-control --wait 1 --performance --gpu-hbf-batch-size 4 --timeout 1800
python3 tools/run_trgtime_pipeline.py --build-dir ./build --run-dir "$HOME/GitHub/nestdaq-run" \
  --output-dir /tmp/nestdaq-hbf4-full-details --wait 1 --performance --gpu-hbf-batch-size 4 --gpu-exec-details --timeout 1800
python3 tools/compare_trgtime.py /tmp/nestdaq-hbf1-full-control/trgtime.csv /tmp/nestdaq-hbf4-full-control/trgtime.csv
root -l -b -q 'tools/compare_gpu_hbf_batch.C("/tmp/nestdaq-hbf1-full-control/logicfilter.root","/tmp/nestdaq-hbf4-full-control/logicfilter.root","test/hbf-batch-comparison",10)'
root -l -b -q 'tools/analyze_gpu_exec.C("/tmp/nestdaq-hbf4-full-details/logicfilter.root","/tmp/nestdaq-hbf4-full-control/logicfilter.root","test/hbf-batch-exec",10)'
```

比較はトリガーあり／なし／全体の平均値と、Exec/Elapsed Timeの分布を出力する。
キャンバス800×600、線幅3、時間分布200 bins、PNG/PDF/ROOT。
平均値の棒グラフは時間分布ではなく3カテゴリのbinを使う。
