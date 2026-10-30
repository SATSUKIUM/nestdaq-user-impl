# LogicFilter: cuda-visible-devices

GPU版の文字列オプション `--cuda-visible-devices` を追加しました。

```bash
LogicFilter --cuda-visible-devices 0
LogicFilter --cuda-visible-devices 1
LogicFilter --cuda-visible-devices GPU-xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
LogicFilter --cuda-visible-devices 1,0
```

空文字（デフォルト）は `setenv` を呼ばず、継承した `CUDA_VISIBLE_DEVICES` を維持します。
未設定と空の環境変数は区別します。非空オプションは継承環境変数を上書きします。
CUDAの数値・UUID・可視リストを独自の数値検証なしで渡します。
`CUDA_VISIBLE_DEVICES=1` の物理GPU 1は通常プロセス内では論理GPU 0です。
カンマ区切りは可視リストであり、複数GPUへの並列分散機能ではありません。
可視リスト先頭を使用する既存挙動を保ち、`cudaSetDevice` は追加していません。

InitTaskでは、指定値と適用後の環境変数を最初のCUDA APIより前にログ出力します。
その後、チェック付きの `cudaGetDevice` / `cudaGetDeviceProperties` /
`cudaDeviceGetPCIBusId` で論理番号・GPU名・PCI bus IDを取得してログ出力します。
setenvまたはCUDA APIが失敗した場合はerrorログと例外で初期化を失敗させます。
別GPUへのフォールバックは追加していません。CUDA自身がリストの一部だけを
受理する場合の挙動はCUDAに従います。

同一プロセスでINIT TASKを繰り返す場合、実効可視リストが同じなら受け付けます。
最初のCUDA API呼び出し直前にプロセス寿命の設定を固定するため、CUDA初期化に
失敗した後も別リストへの変更は拒否します。変更時にはプロセスを再起動してください。
同じ割当を表す別表記でも文字列が異なれば拒否します。
再初期化時にはホスト固定メモリとdecoderの既存領域を解放して再確保します。

## Redis / mq-param

`nestdaq-run/mq-param/mq-param_TFBFP-LF-TFS-Nothing.sh` に無効のコメント例を追加しました。
共通パラメータ・起動数・トポロジーは変更していません。

```bash
param LogicFilter hbf-batch-size 4
# 次回計測準備時に必要な行を有効化
# param LogicFilter-0 cuda-visible-devices 0
# param LogicFilter-1 cuda-visible-devices 1
# param LogicFilter-0 cuda-visible-devices GPU-xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
```

共通キー `parameters:LogicFilter` → 個別キー `parameters:LogicFilter-0` の順で
既存ParameterConfigPluginが設定するため、個別値が優先されます。
CLIとの優先順位も既存プラグインに任せます（Redis読取りが既存プロパティを更新）。
専用Redisクライアントやシェルでの環境変数展開は追加していません。

手元のNestDAQ `eda728718362cc0ee169bea9b111ffc3fda7ca6b` の
`plugins/ParameterConfigPlugin.cxx` を確認しました。カンマなしは文字列ですが、
カンマ入りは `Parse` → `ToArray` で `vector<string>` に変換されます。
LogicFilterは通常文字列で取得し、その型の場合だけ配列をカンマで連結して受け取ります。
この既存プラグインは空白も区切りにするため、Redisでは `1,0` のような通常の
区切り表記を使用してください。CUDAへ渡す元の任意文字列を完全に保存する仕組みではありません。

`launcher/Tracking.sh` → mq-param → `start_device.sh`、および
`-P parameter_config --parameter-config-uri tcp://127.0.0.1:5921/2` も確認済みです。
スクリプト冒頭にはflushdbがあるため、検証ではsourceしていません。
UUIDは実行ホストで調べた値を使用してください。他PCのUUIDを流用しないでください。

## CUDA初期化の順序とスレッド

`LogicFilter` コンストラクタは `new Trigger()` とKTimer生成だけです。
`Trigger() = default` のメンバ初期化はCPU側のポインタ・数値・vectorのみで、
KTimerもchronoのみです。Initオーバーライドはありません。
従来の最初のCUDA APIはInitTask内の `SetTimeRegion` にある `cudaFree` です。
今回の設定処理はInitTask先頭で、情報確認API・SetTimeRegion・cudaMallocHostのすべてより前です。
Triggerとcudamapの確認対象ソースには先行するCUDA初期化処理は見つかりませんでした。

手元のFairMQソースではDeviceRunnerがプラグインを初期化してから
`RunStateMachine()` を呼びます。Deviceの状態ハンドラからInitTaskWrapperと
RunWrapperが呼ばれ、StateMachine::ProcessWorkの同じ処理スレッドで実行されます。
ConditionalRunからTrigger::ExecBatchを呼び、GPU用の独自ワーカースレッドはありません。
デストラクタでCUDA領域を解放します。解放が別スレッドでも、そのスレッドの既定論理
デバイス0はプロセスの同じ可視リスト先頭を指します。
確認したFairMQとNestDAQプラグインソースには先行CUDA APIは見つかりませんでした。
将来CUDAを使うプラグインを追加する場合、この順序条件を再確認する必要があります。

## 検証

作業用コピーの `nestdaq-user-impl/` で次を実行しました（installなし）。

```bash
NESTDAQ=/home/kashima/local/nestdaq cmake -S . -B build \
  -DCMAKE_INSTALL_PREFIX=/home/kashima/local/nestdaq \
  '-DCMAKE_PREFIX_PATH=/home/kashima/local/nestdaq;/home/kashima/local/nestdaq/src/uhbook;/home/kashima/local/root' \
  -DCMAKE_LOGICFILTER_CUDA_ENABLE=1 -DCMAKE_CUDA_ARCHITECTURES=86
NESTDAQ=/home/kashima/local/nestdaq cmake --build build --target LogicFilter --parallel 2
build/LogicFilter/GPU/LogicFilter --help
```

CUDA有効ビルド成功、helpにCUDA_VISIBLE_DEVICESとの対応と空文字の意味が表示されました。
ビルドには既存のFairMQ非推奨API等の警告があります。
`configure.log`, `build.log`, `build-final.log`, `help.log` に結果を保存しています。
CPU版の変更はありません。

作業資料ディレクトリで実行した追加確認:

```bash
python3 tests/check_configure_cuda.py
bash -n nestdaq-run/mq-param/mq-param_TFBFP-LF-TFS-Nothing.sh
nvidia-smi -L
```

シェル構文確認成功。ConfigureCudaの実装本体を抽出した模擬APIテスト11ケース成功。
未指定・継承・空の継承値・明示上書き・UUID・CLI/Redisのカンマリスト・setenv失敗・
3つのCUDA情報取得APIそれぞれの失敗を確認しました。
成功ケースでは同一プロセスで同じ設定の再実行と設定変更拒否も確認しています。
模擬CUDA API側で環境変数の適用順序も検証しています。結果は `unit-tests.log`。

`nvidia-smi -L` はNVIDIAドライバと通信できず失敗しました。
実GPUでの数値/UUID/無効指定、実FairMQでのINIT TASK繰返し、隔離Redisとの統合、
2GPU・2プロセスの実割当は未確認です。模擬テストは実GPU/FairMQ統合試験の代替ではありません。
DAQの起動・RUN、Redisへの接続・変更、install、デプロイ、pushは実施していません。

## 作業場所と変更の持ち運び

元リポジトリは両方とも作業開始時・完了時に変更なしです。
すべての成果物は `workspace/20261022/logicfilter-cuda-visible-devices/` にあります。
作業用ディレクトリは元の追跡ファイルをコピーしたものです（gitリポジトリではありません）。

- nestdaq-user-impl: `LogicFilterMore32GPU`, `17bc9438967feda5ff0e659021d6576c0e53d582`
- nestdaq-run: `main`, `0edba7c077561b5526938150c80bf7ac1f656f5d`

変更ファイルはGPU版 `LogicFilter.cu`, `Trigger.cxx`, この説明のコピー
`LogicFilter/GPU/CudaVisibleDevices.md` と、mq-paramのコメント例です。
`nestdaq-user-impl.patch` / `nestdaq-run.patch` はそれぞれの基準からの未コミット差分で、
各リポジトリに `git apply --check <patch>` で確認して運搬できます。
反映先にローカル変更がある場合はその内容を保持してレビューしてください。
