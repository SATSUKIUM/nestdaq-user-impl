# GPU Exec内部の時間内訳

4 HBFをまとめた新しい実装と全件比較は[HBF一括処理](logicfilter-hbf-batch.md)を参照。

`--exec-details true`を追加し、GPU版のExec内部をホスト側の経過時間で分割した。
既定は無効。`--performance-output`で記録を有効にした場合だけ使用できる。
ROOTの`logicfilter` TTreeへ`gpu_exec_*`ブランチを追加し、各ConditionalRun内の
全Execの値を合算して保存する。既存の粗い時間項目とCPU版との互換性は維持する。

## segから直接転送し、MarkをHBF時刻ごとに起動

`Exec()`は各segのホストポインタからGPU入力バッファの独立領域へ
`cudaMemcpyAsync`で直接コピーする。ホスト上で生データを詰め直す処理はない。
別途`offset/count/femId/femType/firstBlock`の記述子を転送し、全segを扱う
Markカーネルを1回だけ起動する。同じデフォルトストリームでコピー→Mark→Scanを
順序付ける。ページ可能ホストメモリのため、コピーAPIが完全非同期とは限らない。
各blockは記述子を二分探索して担当segを特定し、大きいsegには複数blockを割り当てる。
入力・記述子のGPUバッファは容量不足時だけ拡張して再利用する。

39基板×4 HBFでは、データコピーは156回/TFのまま、Mark起動が156回から4回/TFへ減る。
記述子コピーは4回/TF（この環境で合計4992 bytes/TF）追加される。
ROOTの`gpu_exec_mark_calls`は従来どおり入力seg数で、新しい
`gpu_exec_mark_launches`が実際のカーネル起動回数、`gpu_exec_metadata_h2d_bytes`が
記述子の転送量を表す。従来の`gpu_exec_h2d_bytes`は生データだけを数える。
時間の`h2d_host_us`には生データと記述子の両方を含む。

同じ2000 TF・wait=0・先頭10 TF除外、詳細計測なしで比較した。
変更前は直前のMark同期削除版。別run間の比較なので実行時のばらつきも含む。

| 条件 | segごとにMark起動 [us/TF] | HBF時刻ごとにMark起動 [us/TF] | 短縮率 |
|---|---:|---:|---:|
| トリガー1件以上 | 1286.6 | 472.0 | 63.3% |
| トリガー0件 | 1219.8 | 422.3 | 65.4% |
| 全体 | 1252.8 | 446.8 | 64.3% |

詳細計測ありのトリガーあり群ではExec 490.1 us/TF、入力・記述子H2Dが301.5 us、
Mark起動が8.8 us、件数D2Hが105.1 us（MarkとScanの完了待ちも含む）。
詳細計測の対照runとの差は同群+3.8%。
この変更はコピー方式・GPUバッファ配置・カーネル起動単位を同時に変えるため、
改善率をカーネル起動費用だけに帰属させることはできない。

4 TDC形式・39入力×4 HBF・計測3モードの回帰テストに加え、異なるFEM/IPと
TDC形式の混在、不均等なデータ長、複数block、空seg、未登録FEMをテストし、
従来のseg単位Markと同じ結果になることを確認した。CPU/GPU版ともビルド成功。
実データは変更前／変更後／詳細計測で2000 TF・8000 META・16126 TrgTime・
16126 slicesが順序込みで一致した。詳細記録でもMark起動4回/TFを確認した。

出力:

- 対照: `/tmp/nestdaq-gpu-batch-mark-control/logicfilter.root`
- 詳細: `/tmp/nestdaq-gpu-batch-mark-details/logicfilter.root`
- 解析（PDF/PNG/CSV/ROOT）: `/tmp/nestdaq-gpu-batch-mark-analysis/`
- 前後の詳細計測なし比較: 同解析ディレクトリの`before_after.csv`

## Markごとの同期削除後（seg単位起動）

39 STF（読み出し基板）×各4 HBFという構造に合わせ、同じHBF時刻の39基板分を
MarkしてからScanする。Mark末尾の`cudaDeviceSynchronize()`を削除した。
同一デフォルトストリームの転送→Mark→次の転送→Scanという順序で、
再利用する`devMarkWords`への上書きとTimeRegionへのアクセス順序を保証する。
Scanの同期版D2Hで結果を取得してから次のHBFへ進むため、Markに明示的な同期は不要。
Mark内の明示的デバイス同期は156回/TFから0回/TFへ減った。
初期化等の同期は残し、入力コピー156回/TFと結果取得4回/TFは変更していない。

同じ2000 TF・wait=0・先頭10 TF除外で、詳細計測なしの平均Exec時間を比較した。
変更前後は別runなので、差には実行時のばらつきも含む。

| 条件 | 変更前 [us/TF] | 変更後 [us/TF] | 短縮率 |
|---|---:|---:|---:|
| トリガー1件以上 | 1374.0 | 1286.6 | 6.4% |
| トリガー0件 | 1336.4 | 1219.8 | 8.7% |
| 全体 | 1355.0 | 1252.8 | 7.5% |

詳細計測ありのトリガーあり群では、Exec 1321.5 us/TFのうち入力H2Dが864.6 us
（65.4%）、Mark起動が298.6 us（22.6%）。変更前の同期待ちの多くが次のH2Dへ
移っており、同期区間666.7 usが丸ごと短縮されるわけではなかった。
詳細計測による対照runとの差は同群で+2.7%。次の改善候補は入力転送の一括化。

39入力で異なるヒットを交互に書き込み、4 HBFでTimeRegionを再使用する回帰テストを、
4 TDC形式・計測3モードで実行して通過した。実データでも変更前／変更後／詳細計測の
2000 TF・8000 META・16126 TrgTime・16126 slicesが順序を含めて一致した。

出力:

- 対照: `/tmp/nestdaq-gpu-no-mark-sync-control/logicfilter.root`
- 詳細: `/tmp/nestdaq-gpu-no-mark-sync-details/logicfilter.root`
- 解析（PDF/PNG/CSV/ROOT）: `/tmp/nestdaq-gpu-no-mark-sync-analysis/`

`gpu_exec_mark_sync_host_us`は旧データとの互換性のため残し、変更後は0となる。

## 変更前の主な結果（Markごとに同期）

同じrun000400.datの先頭2000 TF、wait=0。最初の10 TFを除く1990 TFを集計し、
トリガーあり982 TF・なし1008 TFに分けた。入力・TF順序・転送バイト数の整合を検査し、
計測なし／ホスト計測／CUDAイベント併用の各実行で、TrgTime 16126件が全件一致した。

トリガーありの、ホスト区間計測だけを有効にした実行の平均:

| Exec内の区間 | 平均 [us/TF] | Exec内の割合 |
|---|---:|---:|
| Mark後のcudaDeviceSynchronize | 666.7 | 45.6% |
| 入力データのH2D cudaMemcpy | 343.4 | 23.5% |
| Markカーネル起動＋cudaGetLastError | 288.5 | 19.8% |
| 件数D2Hコピー（scan完了待ちも含む） | 59.9 | 4.1% |
| CUB exclusive scanのホスト呼び出し | 22.4 | 1.5% |
| ヒット時刻配列D2Hコピー | 16.9 | 1.2% |
| SubTimeRegion集約カーネル起動 | 7.7 | 0.5% |
| LUT/rising-edge判定カーネル起動 | 7.8 | 0.5% |
| compactionカーネル起動 | 7.7 | 0.5% |
| その他（バッファ管理・resize・ループ・計測処理等） | 約40 | 約2.7% |

詳細計測時のExec全体は1460.6 us/TF。
入力転送＋Mark起動＋Mark後同期に約88.9%が集中している。
計測対象では1 TFあたりExecを4回呼び、合計156回の入力コピー・Mark起動・同期を行う。
トリガーあり群の入力H2Dは平均401466 bytes/TFで、1コピーあたり約2.6 kB。
トリガーなし群でも呼び出し数は156回/TFで、コピーするバイト数は平均52473 bytes/TF。
この繰り返しの回数を減らす処理構成が、改善箇所の優先候補となる。
今回、転送のまとめ方や同期の回数など、処理アルゴリズム自体は変更していない。

## ホスト計測とCUDA計測の意味

主計測は`std::chrono::high_resolution_clock`で、API呼び出しをホストから見た経過時間。
カーネルは非同期なので、SubTR・LUT・scan・compactionの各ホスト値は
カーネル計算時間そのものではない。これらの完了待ちは後続の件数D2Hに含まれる。
H2D APIはホスト側のステージングや先行処理待ちを含み得る。
Mark後の同期は、先行する転送・Markカーネルの完了待ちとCUDAランタイム処理を含む。
したがって666.7 usを純粋なGPUカーネル計算時間とは解釈しない。

ホスト計測では追加のCUDA同期を行わず、既存のセグメントごとの同期を維持する。
`buffer_host_us`には容量確認・必要な場合のcudaFree/cudaMallocを含む。
`hits_resize_host_us`はホストヒット配列のresize。
`other_host_us`は内部タイマーの対象外のループ・呼び出し等の残差。
`total_host_us`はExec内部の入口～末尾で、外側の`exec_us`との差は解析時に
`outer_wrapper_host_us`として補う。ホスト各項目の合計は外側のExec時間と一致する。

`--exec-cuda-events true`も指定すると、同じストリームにCUDAイベントを挿入し、
H2D、Mark、SubTR、LUT、scan、compaction、件数D2H、時刻配列D2Hの区間を記録する。
区間間に同期は追加せず、最後の既存コピーが完了した後に最終イベントを待って回収する。
イベントを使う場合だけ、末尾のcudaEventSynchronizeが追加される。
イベント区間にはホスト投入の隙間も含まれ得るため、純粋な演算時間ではない。
CUDAイベント値とホスト時間は重複しており、足し合わせないこと。

## 計測による影響

ホスト計測の対照実行と計測実行を同じバイナリで順番に行った結果:

| 条件 | 詳細なしExec [us/TF] | ホスト詳細ありExec [us/TF] | 変化 |
|---|---:|---:|---:|
| トリガー1件以上 | 1374.0 | 1460.6 | +6.3% |
| トリガー0件 | 1336.4 | 1376.9 | +3.0% |
| 全体 | 1355.0 | 1418.2 | +4.7% |

これは各条件1回ずつの観測差で、タイマーの負荷と実行時のばらつきを含む。
イベントを全セグメントへ挿入した別実行では、トリガーありのExecが約1400 usから
約2816 usへ増加した。約2倍の変化があるため、このモードの値を通常運転の
内訳として使わず、補助計測として区別した。

`instrumentation_host_us`に直接記録するのはイベント準備・record・回収等の呼び出し時間。
全タイマーの読み出しや呼び出しの費用、実行順序への影響が全部この項目に入るわけではない。
この項目だけを差し引いて「計測なしの正確な内訳」を復元することはできない。

## 再現手順

```sh
cmake --build build --target LogicFilter -j4
python3 tools/run_trgtime_pipeline.py \
  --build-dir ./build --run-dir "$HOME/GitHub/nestdaq-run" \
  --output-dir /tmp/gpu-exec-control-run --iterations 2000 --wait 0 --performance --gpu-hbf-batch-size 1
python3 tools/run_trgtime_pipeline.py \
  --build-dir ./build --run-dir "$HOME/GitHub/nestdaq-run" \
  --output-dir /tmp/gpu-exec-detail-run --iterations 2000 --wait 0 --performance --gpu-hbf-batch-size 1 --gpu-exec-details
root -l -b -q 'tools/analyze_gpu_exec.C("/tmp/gpu-exec-detail-run/logicfilter.root","/tmp/gpu-exec-control-run/logicfilter.root","/tmp/gpu-exec-analysis",10)'
```

CUDAイベントを併用する補助計測は`--gpu-exec-cuda-events`を指定する。
これは`--gpu-exec-details`も有効にする。CUDAイベントはExecをまたいで再利用し、
最初の確保等の影響を除くため、今回も全体の先頭10 TFを統計から除外した。

今回の出力:

- 主計測: `/tmp/nestdaq-gpu-exec-host-details/logicfilter.root`
- 主計測の対照: `/tmp/nestdaq-gpu-exec-host-control/logicfilter.root`
- 主計測の解析: `/tmp/nestdaq-gpu-exec-host-analysis/`
- イベント併用の記録: `/tmp/nestdaq-gpu-exec-details/logicfilter.root`
- イベント併用の解析: `/tmp/nestdaq-gpu-exec-analysis/`

解析ディレクトリには以下を保存する。

- `exec_stages.csv`: トリガーあり／なし／全体の各区間の平均・中央値・95%点・割合。
- `exec_counts.csv`: Mark回数、H2D/D2Hバイト数、バッファ拡張回数。
- `profiling_overhead.csv`: 対照実行と詳細実行のExec時間の比較。
- `host_time.pdf/png`、`host_fraction.pdf/png`: ホスト時間の構成図。
- `gpu_exec.root`: 構成図、区間別200 binsの分布ヒストグラム、キャンバス。
- イベント併用時のみ`cuda_intervals_*.pdf/png`: CUDAイベント区間の平均。

キャンバス800×600、線幅3。構成図の分類binと、時間分布の200 binsは別の設定。
GPU回帰テストでは、詳細なし・ホストのみ・イベント併用で4種のTDC形式の出力が
一致すること、空のHBF、転送件数・バイト数を確認した。CPU版もビルド確認済み。
実データの全受信TFについて、区間和、入力セグメント数とMark回数、H2D量、
件数のD2H量（4 bytes×Exec回数）、ヒットD2H量（4 bytes×トリガー数）の整合を検査済み。
