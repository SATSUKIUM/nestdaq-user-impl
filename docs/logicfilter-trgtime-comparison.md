# CPU/GPU LogicFilterの受信TrgTime比較

ブランチ: `verify-logicfilter-cpu-gpu-trgtime`

`FilterTimeFrameSliceABC::ConditionalRun()`で、受信した全TrgTimeと各スライスに
対応付けたTrgTimeを、`ProcessSlice()`の採否判定前にCSVへ記録する。
時刻は4 ns単位の整数のまま保存し、トリガー種別・TF ID・順序も比較する。
ログの時刻、処理時間、run番号、追跡処理の採否は比較に含めない。

## ビルド

リポジトリのルートで実行する。CPU/GPUのビルド先を分け、バイナリの取り違えを防ぐ。
インストールせず、次のスクリプトからビルド先の実行ファイルを直接使用できる。

```sh
cmake -S . -B /tmp/nestdaq-trgtime-cpu \
  -DCMAKE_INSTALL_PREFIX="$HOME/local/nestdaq" \
  -DCMAKE_PREFIX_PATH="$HOME/local/nestdaq;$HOME/local/nestdaq/src/uhbook" \
  -DCMAKE_LOGICFILTER_CUDA_ENABLE=0 -DCMAKE_CUDA_ARCHITECTURES=86
cmake --build /tmp/nestdaq-trgtime-cpu -j4 --target \
  LogicFilter FilterTimeFrameSliceByTrack FileSink TFBFilePlayer TimeFrameSlicerByLogicTiming

cmake -S . -B /tmp/nestdaq-trgtime-gpu \
  -DCMAKE_INSTALL_PREFIX="$HOME/local/nestdaq" \
  -DCMAKE_PREFIX_PATH="$HOME/local/nestdaq;$HOME/local/nestdaq/src/uhbook" \
  -DCMAKE_LOGICFILTER_CUDA_ENABLE=1 -DCMAKE_CUDA_ARCHITECTURES=86
cmake --build /tmp/nestdaq-trgtime-gpu -j4 --target \
  LogicFilter FilterTimeFrameSliceByTrack FileSink TFBFilePlayer TimeFrameSlicerByLogicTiming
```

ROOTを自動検出できない環境では`-DROOT_DIR="$HOME/local/root/cmake"`も指定する。
GPU実行には対応するCUDAデバイスとドライバが必要。

## 隔離した同一トポロジーで比較

```sh
python3 tools/run_trgtime_pipeline.py \
  --build-dir /tmp/nestdaq-trgtime-cpu \
  --run-dir "$HOME/GitHub/nestdaq-run" \
  --output-dir /tmp/trgtime-cpu-run

python3 tools/run_trgtime_pipeline.py \
  --build-dir /tmp/nestdaq-trgtime-gpu \
  --run-dir "$HOME/GitHub/nestdaq-run" \
  --output-dir /tmp/trgtime-gpu-run

python3 tools/compare_trgtime.py \
  /tmp/trgtime-cpu-run/trgtime.csv /tmp/trgtime-gpu-run/trgtime.csv
```

出力ディレクトリは新しいパスを指定する。両方の実行が`COMPLETE`で終了したことを確認する。
スクリプトは、TFBFilePlayer → LogicFilter → TimeFrameSlicerByLogicTiming →
FilterTimeFrameSliceByTrack → FileSinkを、`/tmp`配下の専用IPCソケットで接続する。
既存のtmuxセッション、Redis設定、インストール済みバイナリは変更しない。
FileSinkの出力、ROOT出力、デバイスログ、実行コマンドJSONも出力ディレクトリへ保存する。

設定は今回確認したTracking構成の固定値を使用する。元のshellをsourceしないため、
今後設定を変更した場合はスクリプト内の指定も揃える。

- 入力: `RUN_DIR/data/run000400.dat`の全件（`--input`で変更、`--iterations N`指定時のみ件数制限）。
  全件モードではレコード境界を検査してTF数を数え、Playerに`--max-iterations`を渡さない。
- UTOF 8/10 chと、AAモジュール16/17、18/19、20/21、22/23、24/25、27/28 chのサブグループ。
- 式: `(0 & 1) & 2`、width: 30、AAモジュールのdelay: -12。
- スライサーの範囲: -250～250（4 ns単位）。
- channel map: `mapdata_T103_20260909.csv`。
- 検出器設定: `param_DetGeom001`、`param_DCTdcCalib_run0602`、`param_DCDrift_run0602`。

既定は`--wait 1`（1 ms間隔）。`--wait 0`で送信間隔なしにできる。Playerへ`--completion-ack-file`を渡し、
EOFまたは件数上限に到達してもソケットを維持する。実行スクリプトは、ByTrackが全メッセージ・META・
スライス・TrgTimeを処理したことをCSVから確認してackファイルを作成する。
FileSinkも`--max-iteration`で全件の書き込み後に終了させるため、静止時間だけで
キューが空になったとは判断しない。受信件数と正常終了も最後に再検査する。
TFBFilePlayerの新しいオプションを使うため、Playerも再ビルドすること。

## 既存のlauncherで記録する場合

ByTrackの起動オプション`--trgtime-output /絶対パス/cpu.csv`、またはmq-paramの
`trgtime-output`で出力先を指定する。CPU/GPUでは異なるファイル名にする。
空のままなら環境変数`NESTDAQ_TRGTIME_OUTPUT`を参照し、それも空なら記録しない。
既存tmuxサーバーは任意の環境変数を自動継承しない場合があるため、launcherでは
mq-paramまたは明示的な起動オプションでの指定を推奨する。

mq-param設定例（デバイスのInitTask前に適用）:

```sh
param FilterTimeFrameSliceByTrack trgtime-output /tmp/trgtime-cpu.csv
```

同一入力・トリガー設定・スライサー設定・TF数で実行し、下流まで処理が完了してから
正常にSTOPする。ファイルはPreRunで上書きされ、PostRunでENDレコードを出力する。
`kill -9`や異常終了した記録は比較しない。

## CSVと判定

```text
# nestdaq-trgtime-v1,time_unit=4ns
MESSAGE,message_index,number_of_TFs,number_of_LF_blocks
META,message_index,meta_index,tf_id,lf_tf_id,number_of_triggers,number_of_time_blocks
TRG,message_index,meta_index,tf_id,time_block_index,trigger_index,time_4ns,type
SLICE,message_index,meta_index,tf_id,slice_index,time_4ns,type
END,total_messages,total_metas,total_triggers,total_slices
```

TRGは受信した全トリガー、SLICEは実際にProcessSliceへ渡す対応値を記録する。
トリガー0件のMETAも保存するため、ゼロ件TFの欠落も検知できる。
比較は全レコードを順番どおり照合し、最初の差分を表示する。
終了レコードがない場合、件数不整合、トリガーまたはスライスが全くない場合は失敗する。
同じ途中範囲同士の一致は全入力範囲の一致を保証しないため、実行スクリプト側の
受信件数検査も併用する。

## この環境での検証

CPU版は実データの先頭2000 TFを5デバイス経由で処理し、全プロセス正常終了。
受信記録は2000メッセージ、8000 META、16126 TrgTime、16126スライス。
全スライスの時刻・種別が対応する受信TRGと一致した。
RTX 3070 Ti（sm_86、driver 615.71.09）でGPU版も実行し、CPU版と全レコード一致。
CPU版の再実行を含む3ファイルのSHA-256も同一だった。
GPUとIPCソケットはサンドボックス内から利用できなかったため、実行は権限承認後に
サンドボックス外で行った。既存DAQのRedis・tmuxは使用していない。

```text
MATCH: messages=2000, metas=8000, TrgTime=16126, slices=16126
SHA-256: 017590a2aa77584884fd5ae77d9bf77afb69e887aae75fc529c0d76d81797412
```

今回の記録・ログ・実行コマンド:

- CPU: `/tmp/nestdaq-trgtime-cpu-reference/`
- CPU再実行: `/tmp/nestdaq-trgtime-cpu-repeat/`
- GPU: `/tmp/nestdaq-trgtime-gpu-reference/`

各ディレクトリの`trgtime.csv`が比較対象。今回の入力範囲と設定での一致確認であり、
送信待機とCSV出力を含むため性能測定には使用しない。

比較ツールのテスト:

```sh
python3 -m unittest discover -s tools/tests -v
```

全件確認用にByTrackの元の`ProcessSlice()`を`#if 0`で囲み、`false`だけ返す実装を
有効にしている。TrgTimeと受信スライスの記録はその呼び出し前なので照合は継続できる。
トラック再構成とスライスの出力は無効。詳細は[全件テスト](logicfilter-full-file-test.md)参照。

GPU版は既定で最大4 HBFを一括処理する。`--gpu-hbf-batch-size 1`でHBF単位処理と
比較できる。境界・出力形式・全件の照合は[HBF一括処理](logicfilter-hbf-batch.md)参照。
