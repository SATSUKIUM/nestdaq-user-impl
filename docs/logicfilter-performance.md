# LogicFilterのCPU/GPU処理時間・データ量

CPU版・GPU版とも`--performance-output /absolute/path/logicfilter.root`で記録を有効にする。
空なら環境変数`NESTDAQ_LOGICFILTER_PERFORMANCE_OUTPUT`を参照し、それも空なら記録しない。
PreRunでROOTファイルを新規作成し、PostRunで保存・closeする。同じパスは上書きされる。

## 計測区間

`LogicFilter`クラスに`std::chrono::high_resolution_clock::time_point`を持たせ、
ConditionalRun内の受信・cleanup・Exec・送信の前後で`now()`を取得する。
CPU/GPUの保存形式は共通の`LogicFilter/PerformanceRecorder.h`で定義している。
TTree名は`logicfilter`、1行がConditionalRunの1回の呼び出し。時間の単位はすべてus。

| branch | 内容 |
|---|---|
| `conditional_us` | 関数入口からローカルメッセージの破棄後まで。ROOT Fill自身は除外 |
| `receive_us` | Receive呼び出し。入力待ちを含む |
| `cleanup_host_us` | 各HBFでCleanUpTimeRegion／CleanUpSubTimeRegionを呼ぶホスト経過時間の合計 |
| `exec_us` | 各HBFのExec呼び出し時間の合計 |
| `trigger_us` | cleanup開始～Exec終了の合計。`cleanup_host_us + exec_us` |
| `send_us` | 主出力のPoller作成・poll・Send・再試行待ちの合計 |
| `accounting_us` | データ量・TDCワードを数える計測用処理の時間 |
| `processing_us` | `conditional_us - receive_us - send_us - accounting_us` |
| `elapsed_no_receive_us`（比較時に算出） | `conditional_us - receive_us`。送信待ち・集計時間も含む |
| `other_us` | `processing_us - trigger_us`。HBF構築・メッセージ生成・コピー・ログなど |

GPUのcleanupは非同期カーネル起動なので、`cleanup_host_us`はGPUのcleanup実行時間ではない。
現在のGPU版ExecはMarkの同期、および圧縮ヒット配列の同期コピーを行ってから戻るため、
`trigger_us`にはcleanupのGPU実行・転送・完了待ちも含まれる。
計測のための追加cudaDeviceSynchronizeは挿入していない。
DQMが有効な場合、そのコピー・転送待ちは`other_us`に含まれる。

`received=false`の行は受信タイムアウト等。全呼び出しの記録を残すが、比較ツールは
`received=true`だけを集計する。計測用ROOT書き込みの時間は処理時間へ含めない。
TDC集計はExec後に実行し、計測前に入力データをスキャンしてキャッシュを暖めることを避ける。
計測自身や下流の追跡・CSV・ROOT出力による影響は残るので、純粋なカーネル性能ではなく
このトポロジー内でのホスト側経過時間を比較する。

## 処理データ量と識別子

| branch | 内容 |
|---|---|
| `iteration`, `tf_id`, `cuda_enabled` | 呼び出し連番、TimeFrame ID、GPU版かどうか |
| `received`, `receive_result`, `sent` | 受信成功・Receive戻り値・主出力送信成功 |
| `input_parts`, `input_bytes` | inPartsの要素数、全partのGetSize()合計 |
| `hbf_count` | HBF処理ループの回数（bsize_min） |
| `segment_count` | Execに渡したHBFIndexセグメント数の合計 |
| `trigger_input_bytes`, `trigger_input_words` | Execへ渡したセグメントのバイト数・完全な64bitワード数。HBFヘッダー等も含む |
| `payload_words` | HBFヘッダーを除いた64bitワード数。HB等の制御ワードも含む |
| `leading_tdc_words` | TDC64H/Lと各V3形式のleading TDCワード数（type=0x0b） |
| `trailing_tdc_words` | 同じ形式のtrailing TDCワード数（type=0x0d） |
| `tdc_words` | leading＋trailing。設定チャンネルへの一致や時間窓による選別前の個数 |
| `triggers` | LogicFilterが検出したトリガー数の合計 |
| `output_parts`, `output_bytes` | Send直前の主出力メッセージのpart数とバイト数 |

ROOTファイルには`implementation`、`trigger_signals`、`trigger_expression`、`trigger_width`、
`clock`、`clock_is_steady`、`schema_version`、`completed`も保存する。
高分解能クロックの単調性は環境依存なので、`clock_is_steady`も記録する。
`completed=1`は正常にPostRunで保存されたもの。異常終了・Resetによる保存は比較しない。

## 再実行と比較

[TrgTime比較手順](logicfilter-trgtime-comparison.md)のCPU/GPUビルドを更新した後、
既存の実行スクリプトに`--performance`を追加する。

```sh
cmake --build /tmp/nestdaq-trgtime-cpu --target LogicFilter TFBFilePlayer -j4
cmake --build build --target LogicFilter TFBFilePlayer -j4  # このbuildはCUDA_ENABLE=1で構成済み

python3 tools/run_trgtime_pipeline.py \
  --build-dir /tmp/nestdaq-trgtime-cpu --run-dir "$HOME/GitHub/nestdaq-run" \
  --output-dir /tmp/lf-timing-cpu-run --performance
python3 tools/run_trgtime_pipeline.py \
  --build-dir ./build --run-dir "$HOME/GitHub/nestdaq-run" \
  --output-dir /tmp/lf-timing-gpu-run --performance

python3 tools/compare_trgtime.py \
  /tmp/lf-timing-cpu-run/trgtime.csv /tmp/lf-timing-gpu-run/trgtime.csv
root -l -b -q 'tools/compare_logicfilter_timing.C("/tmp/lf-timing-cpu-run/logicfilter.root","/tmp/lf-timing-gpu-run/logicfilter.root","/tmp/lf-timing-comparison",10)'
```

最後の引数10は、最初の受信10 TFを統計から除くwarmup件数。生データは全件保存する。
CPUとGPUは並列に走らせず、同じ入力範囲・設定・送信間隔で順番に実行する。
既定は`wait=0`。全件処理後のack通知までPlayerを維持し、終了時の未配送を防ぐ。
`conditional_us`にはReceiveの呼び出し時間が含まれる。
トリガー比較には主に`trigger_us`、デバイス内処理比較には`processing_us`を見る。

ROOT比較マクロは、正常保存、TF数・順序・各データ量・トリガー数・設定の一致、
時間内訳の整合を検査した上で以下を出力する。

- `summary.csv`: 各時間の平均・中央値・95%点とCPU/GPU平均比。
- `elapsed_no_receive_vs_tdc.png`: type 0x0bのTDC+TOTワード数に対する、`conditional_us - receive_us`。
- `trigger_vs_tdc.png`: 同じ横軸に対する、cleanup＋ExecのElapsed Time。
- `elapsed_no_receive_distribution.png`: Receiveを除いたElapsed Timeの分布。
- `trigger_distribution.png`: CPU/GPUトリガー処理時間の分布。
- `comparison.root`: 上記グラフ・ヒストグラム・キャンバス。

## 以前の測定結果（wait=10 ms）

RTX 3070 Ti、CUDA architecture 86。前回と同じ入力の先頭2000 TFをCPU/GPUで順に処理した。
両方ともConditionalRunは2006回（受信成功2000回、残りは待機）で、各TFの入力・出力サイズ、
HBF数、TDCワード数、トリガー数が一致した。下流のTrgTimeも16126件すべて一致。
先頭10 TFを除いた1990 TFの平均は次のとおり。

| 項目 | CPU [ms/TF] | GPU [ms/TF] |
|---|---:|---:|
| ConditionalRun全体 | 10.178 | 10.175 |
| 受信待ちを含むReceive | 5.970 | 8.310 |
| トリガー処理 | 4.126 | 1.613 |
| その他の処理 | 0.049 | 0.161 |
| 処理合計（受信・送信待ち・集計を除く） | 4.176 | 1.774 |

トリガー処理のCPU/GPU平均比は2.56、処理合計は2.35。
10ms間隔で入力しているため、ConditionalRun全体の差をスループット向上とは解釈しない。
これらは現状のビルド設定での値。CPU CXX_FLAGSには`-O3`がある一方、
LogicFilter.cuのCUDA_FLAGSにはホスト向け`-O3`が指定されていない。
ホスト処理の最適化条件も揃えた比較は、今回の変更範囲では行っていない。

記録ファイルと集計結果:

- CPU: `/tmp/nestdaq-lf-timing-cpu/logicfilter.root`
- GPU: `/tmp/nestdaq-lf-timing-gpu/logicfilter.root`
- 比較: `/tmp/nestdaq-lf-timing-comparison/summary.csv`
- グラフ: `/tmp/nestdaq-lf-timing-comparison/trigger_vs_tdc.png`
- ROOTグラフ: `/tmp/nestdaq-lf-timing-comparison/comparison.root`

ワード集計の回帰テスト（ROOT環境で実行）:

```sh
c++ -I. LogicFilter/tests/PerformanceRecorderTest.cxx \
  $(root-config --cflags --libs) -o /tmp/lf-performance-counter-test
/tmp/lf-performance-counter-test
```

HBFヘッダー除外、leading/trailing識別、非対応FEM形式の除外、HBF間の合算を検査する。

## グラフの軸の定義（wait=0で再測定）

横軸は`leading_tdc_words`を使い、`TDC+TOT words per TF (type 0x0b)`と表記する。
TDCとTOTを別々に数えた値ではなく、1つの64bitワードを1件として数える。
ROOTの既存`tdc_words`（0x0b＋0x0d）は互換性のため残すが、横軸には使用しない。

主グラフの縦軸は`Elapsed time excluding Receive [us/TF]`で、
`elapsed_no_receive_us = conditional_us - receive_us`として比較マクロ内で算出する。
今回はReceiveだけを引くので、送信待ちやデータ量集計時間は含まれる。
ROOTへのFill自体は従来どおり`conditional_us`に含まれない。
別の`trigger_vs_tdc.png`は`Trigger processing elapsed time [us/TF]`と表記し、
cleanup開始～Exec終了の合計を表示する。両者を混同しないこと。

### wait=0の再測定結果

同じ入力2000 TFを使い、CPU/GPUとも全件をByTrackで処理してFileSinkに書き込んだ。
`--completion-ack-file`による完了通知でPlayerを維持したため、送信終了時の欠落はない。
データ量・TF ID・トリガー数が全行一致し、下流TrgTime 16126件も一致した。
最初の10 TFを除いた1990 TFの平均:

| 項目 | CPU [ms/TF] | GPU [ms/TF] | CPU/GPU |
|---|---:|---:|---:|
| Receiveを除くElapsed Time | 4.137 | 1.672 | 2.47 |
| トリガー処理（cleanup＋Exec） | 4.055 | 1.418 | 2.86 |
| ConditionalRun全体 | 4.159 | 1.717 | 2.42 |
| Receive呼び出し | 0.022 | 0.045 | — |

今回の出力先:

- CPU生データ: `/tmp/nestdaq-lf-wait0-cpu/logicfilter.root`
- GPU生データ: `/tmp/nestdaq-lf-wait0-gpu/logicfilter.root`
- 集計: `/tmp/nestdaq-lf-wait0-comparison/summary.csv`
- Receive除外グラフ: `/tmp/nestdaq-lf-wait0-comparison/elapsed_no_receive_vs_tdc.png`
- トリガー処理グラフ: `/tmp/nestdaq-lf-wait0-comparison/trigger_vs_tdc.png`
- ROOTグラフ: `/tmp/nestdaq-lf-wait0-comparison/comparison.root`

Playerの`--completion-ack-file`は既定では空で、既存の運用は変更しない。
指定時はEOFまたはmax-iterations到達後、指定ファイルが作成されるかSTOPを要求されるまで
RUNNING状態を維持する。既存のackファイルを指定すると誤完了を避けるため初期化エラーになる。
比較スクリプト側には全体タイムアウトがあり、異常時には起動したプロセスを終了させる。

## ConditionalRunのトリガー有無で分けた比較

再計測せず、wait=0の保存済みROOTから`triggers >= 1`と`triggers == 0`に分ける。
受信成功行だけを使い、最初の10 TFを全体から除外してから分類する。
CPU/GPUでTFごとのトリガー数とデータ量が一致することも従来どおり検査する。

```sh
root -l -b -q 'tools/compare_logicfilter_trigger_groups.C("/tmp/nestdaq-lf-wait0-cpu/logicfilter.root","/tmp/nestdaq-lf-wait0-gpu/logicfilter.root","/tmp/nestdaq-lf-trigger-groups",10)'
```

Receiveを除くElapsed Time（ConditionalRun − Receive）の結果:

| 条件 | TF数 | CPU平均 [ms] | GPU平均 [ms] | CPU 95%点 [ms] | GPU 95%点 [ms] | CPU/GPU平均比 |
|---|---:|---:|---:|---:|---:|---:|
| トリガー1件以上 | 982 | 4.518 | 1.749 | 5.035 | 2.380 | 2.58 |
| トリガー0件 | 1008 | 3.766 | 1.597 | 3.860 | 2.434 | 2.36 |

cleanup＋Execだけの平均は、トリガーありでCPU 4.416 ms・GPU 1.441 ms（平均比3.06）、
トリガーなしでCPU 3.703 ms・GPU 1.396 ms（平均比2.65）。
平均type 0x0b TDC+TOTワード数は、ありで48196件/TF、なしで5931件/TF。
トリガーあり群の平均入力サイズは約403 kB、なし群は約54 kB。

分類はConditionalRunで見つかったトリガー数に基づく。slow extractionのビーム状態を
直接記録したspillゲートによる分類とは区別する。ビーム供給中でも条件不成立のTFは
トリガー0件群に入る可能性がある。今回の収録順表示には高活動・低活動の繰り返しが見えるが、
横軸は受信TFのインデックスとし、再生時の経過時間を収録時のビーム時間に置き換えない。

出力（`/tmp/nestdaq-lf-trigger-groups/`）:

- `groups_summary.csv`: 各群の全時間項目の平均・中央値・95%点・CPU/GPU平均比。
- `groups_data.csv`: 群ごとのTF数、平均入力バイト数・TDCワード数・トリガー数。
- `elapsed_no_receive_trigger_groups.pdf/png`: CPU/GPU×トリガー有無の4分布。
- `trigger_trigger_groups.pdf/png`: cleanup＋Execについての4分布。
- `activity_by_tf.pdf/png`: 取得順のトリガー数とtype 0x0bワード数。
- `trigger_groups.root`: 上記のヒストグラム・グラフ・キャンバスと選択件数。
- `trigger_present/`と`trigger_absent/`: 各群だけのCPU/GPU比較、TDCワード数との散布図、ROOT。

全キャンバス800×600、線幅3、ヒストグラム200 bins。
4分布を重ねるグラフは、群ごとの件数の違いを除いて形を比較できるよう、
各ヒストグラムの面積を1へ正規化する。個別群のヒストグラムはTF件数のまま保存する。
982＋1008＝1990、全群の合計トリガー数16126、群の件数で重み付けした平均が
未分割の平均を再現することを検査済み。

既存マクロを直接使う場合は、第6引数で`1`（1件以上）、`0`（0件）、`-1`（全件、既定）を指定できる。

## 各実装のConditionalRun内で時間を使っている処理

保存済みwait=0のROOTデータを使い、再計測せず時間内訳を集計した。
受信成功2000 TFの先頭10 TFを除く1990 TFが対象。
グラフの割合はReceiveを除いた時間を分母とする。CSVにはReceiveを含む分母の割合も併記する。

トリガーあり982 TFの平均:

| 区間 | CPU [ms/TF] | CPU内の割合 | GPU [ms/TF] | GPU内の割合 |
|---|---:|---:|---:|---:|
| Exec | 3.452 | 76.4% | 1.380 | 78.9% |
| Cleanup（ホスト呼び出し） | 0.965 | 21.4% | 0.061 | 3.5% |
| その他の処理 | 0.051 | 1.1% | 0.174 | 9.9% |
| 送信（待ちを含む） | 0.012 | 0.3% | 0.019 | 1.1% |
| 計測用データ集計 | 0.039 | 0.9% | 0.115 | 6.6% |

CPUはExecとCleanupで約97.7%を占める。トリガーなしでもCleanupは0.962 ms、
Execは2.742 msで、トリガーが見つからないTFにも相応の固定処理コストがある。
優先して詳しく調べる区間はExecと、毎HBFの領域初期化。

GPUも主な区間はExec。ただしExecにはMark、サブ領域集約、LUT判定、exclusive scan、
compaction、入力・ヒット配列転送、GPU同期が含まれ、非同期cleanupの実行時間の一部もここに入る。
この粗い計測では内部を分離していない。追加したExec詳細計測は末尾のリンクを参照。cleanup_hostはホスト側で観測した呼び出し時間。

その他の処理にはMakeBlockMap、BuildHBF、ヒット配列コピー、送信フラグ更新、
出力ヘッダー・メッセージ生成、DQM経路、ログ、ローカルメッセージの破棄などが含まれる。
これら個別の関数の時間も、現行データでは合算されている。
計測用データ集計は計測負荷であり、performance-output無効時は実行されない。
ROOT Fill自体はConditionalRunの計測区間から外している。

受信を含むConditionalRun全体に対するReceiveの割合は、トリガーありでCPU約0.54%、
GPU約2.72%。今回のwait=0条件ではReceive・主出力Sendの割合は小さい。

```sh
root -l -b -q 'tools/summarize_logicfilter_breakdown.C("/tmp/nestdaq-lf-wait0-cpu/logicfilter.root","/tmp/nestdaq-lf-wait0-gpu/logicfilter.root","/tmp/nestdaq-lf-breakdown",10)'
```

出力先`/tmp/nestdaq-lf-breakdown/`:

- `CPU_time_breakdown.pdf/png`、`GPU_time_breakdown.pdf/png`: 実装ごとの平均時間の積み上げ。
- `CPU_fraction_breakdown.pdf/png`、`GPU_fraction_breakdown.pdf/png`: 実装ごとの時間割合。
- `breakdown.csv`: トリガーあり・なし・全体について、各区間の平均時間と割合。
- `breakdown.root`: 以上のグラフ・構成要素・キャンバス。

各図は800×600、線幅3。3本の棒は「トリガーあり」「トリガーなし」「全体」の分類を表す。
これは分類別の積み上げ図なので、時間分布ヒストグラムの200 bins設定とは独立。
CPU/GPU×3群×2種類の分母について、内訳が計測時間および100%を再現することを確認済み。

GPU Exec内部の詳細計測を追加した。主計測はホスト区間だけを使い、CUDAイベント併用は
負荷の大きい補助モードとして分けている。
[Exec内部の内訳・計測負荷・再現手順](logicfilter-gpu-exec-profile.md)を参照。
