# GPU版LogicFilterの実装

## 読む順序

`LogicFilter.cu`の`ConditionalRun()` → `BuildHBF()` → `Trigger.cxx`の`ExecBatch()` →
`devMarkBatch()` / `markWord()` → `Scan()`の順で追うと、入力からTrgTimeまでを確認できる。
CPU版のソースは変更していない。

## 入力の単位

TFを受信したホスト側で`MakeBlockMap()`が基板ごとのメッセージ位置を整理する。
`BuildHBF(..., i)`は各基板の同じ位置iのHBFを集め、`HBFIndex`の配列を作る。
`HBFIndex`は`femType`、`femId`、元メッセージ内のデータポインタ、バイト数だけを持つ。
生データはホスト上で詰め直さない。

今回の39 STF × 各4 HBFなら、既定の`hbf-batch-size 4`では156入力を持つ
`ExecBatch()`を1回呼ぶ。`hbf-batch-size 1`なら39入力ずつ4回呼ぶ。
HBFごとに独立したTimeRegion/SubTimeRegionを使い、バッチ全体をまとめて初期化・処理する。
`SetTimeRegion(size, nHbf)`のsizeは1 HBF分のbin数で、nHbf個分の容量を確保する。
`PrepareHbfBatch(nHbf)`は実際に処理する個数を切り替えるだけで、メモリ確保や転送はしない。
最後が端数の場合も、確保済み容量内で処理する。
`bsize_min`は各STFのblock_map要素数の最小値であり、今回の入力形式では4。
区切りを別要素として持つ入力形式では、この値をそのまま物理HBF数と解釈できない。

## 初期化時の転送

`InitializeFEEAddrDecoder()`でIPアドレス・chをトリガー設定へ対応付ける。
設定はbit、サブグループ番号、遅延、幅を持つ。
設定配列をGPUへコピーし、デコーダー内の設定配列ポインタをGPUポインタへ差し替えて
デコーダー自身もコピーする。カーネルはCPUのstd::mapを参照しない。

`MakeTable()`はCPUで論理式のLUTを作り、`GetLutBytes()`でバイト配列へ変換して
`devLut`へコピーする。std::vectorオブジェクトをカーネルへ渡す処理はない。
CPUのLUT生成オブジェクトはこの関数内だけで使う。

## ExecBatch：直接コピーと1回のMark起動

まずホストで各segのワード数と必要block数を計算し、`MarkSegment`を作る。

| フィールド | 内容 |
|---|---|
| offset | GPU入力バッファ内の先頭位置。uint64_tワード単位 |
| count | uint64_tワード数 |
| fem | 基板IPアドレス |
| type | TDCモジュール形式 |
| firstBlock | このsegに割り当てる先頭CUDA block番号 |
| hbfIndex | バッチ内のHBF番号。独立した領域を選ぶ |

1 blockは256 threads。各segへ`ceil(count/256)`個のblockを割り当てる。
例えばAが513 words、Bが2 wordsなら、次の配置になる。

| seg | offset | count | firstBlock | 担当block |
|---|---:|---:|---:|---|
| A | 0 | 513 | 0 | 0, 1, 2 |
| B | 513 | 2 | 3 | 3 |

GPUの入力バッファと記述子バッファは、容量が不足した場合だけ拡張する。
ワード用バッファは必要ワード数の2倍を確保し、小さなデータ量増加での再確保を減らす。
各segの元ポインタから`devMarkWords + offset`へ`cudaMemcpyAsync`で直接コピーする。
記述子配列も転送し、全seg分の合計block数で`devMarkBatch`を1回起動する。
39基板×4 HBFをまとめると、生データコピー156回/TF、記述子コピー1回/TF、Mark起動1回/TFとなる。
`cudaMemcpyAsync`の入力がページ可能メモリなら、ホスト側で待つ場合がある。

## Mark：segの識別と時間窓の更新

`devMarkBatch()`は`blockIdx.x`以下の最大の`firstBlock`を二分探索し、担当segを特定する。
同じblock内のthreadsは同じIP・TDC形式・HBF番号を使う。
`hbfIndex`で領域先頭をずらすので、別HBFの時間窓やサブグループは混ざらない。
seg内のワード位置は次の式になる。

```cpp
(blockIdx.x - seg.firstBlock) * blockDim.x + threadIdx.x
```

`markWord()`は担当位置が`count`以上なら何もしない。
TDC64H/Lと各V3形式からch・時刻を取り出し、type 0x0bのTDC+TOTワードだけを使う。
TOT自体はトリガー判定には使わない。
IP・chからデコーダーで設定を引き、遅延と幅を適用する。
窓の両端は包含し、幅は整数除算したhalfWidthで扱う。例えばwidth=30なら
`[hit-15, hit+15]`の31 binsであり、従来CPU版の仕様を維持している。

通常グループは対象binのMainTimeRegionへ`atomicOr`でbitを立てる。
サブグループは、初期状態で立てておいた各信号のbitを`atomicAnd`で落とす。
同じbinへ複数threadが到達するため、atomic操作が必要になる。

SubTimeRegionは初期化時、n信号なら下位n bitsを1にする。n=32はUINT32_MAXを使う。
すべての必要信号のbitが落ちて値0になったbinを、`ScanSubTRandMarkMainTR()`で
対応するMainTimeRegionのbitへ反映する。これでサブグループ内AND、同じMain bitへ
対応するサブグループ間ORを表現する。

## Scan：LUT・exclusive scan・compaction

1. サブグループの成立をMainTimeRegionへ反映する。
2. 各binで`devLut[timeRegion[bin]]`を調べ、直前binが不成立かつ現在binが成立ならflag=1。
3. CUBの`DeviceScan::ExclusiveSum`でflagの出力位置を計算する。
4. flag=1のthreadだけが`hits[offsets[bin]] = bin`を書き込む。
5. 件数だけをD2Hコピーし、ホスト配列をresizeして必要件数の時刻をD2Hコピーする。

例えばflagが`[0,1,0,1]`ならoffsetは`[0,0,1,1]`で、結果は`[1,3]`になる。
出力はバッチ内のglobal bin順。CPU版に合わせ、各HBFのbin 0は出力しない。
HBF境界をまたいで直前binと比較しない。
rising edgeは隣接できないため、出力容量は`nHbf * floor(TimeRegionSize/2)`で足りる。
ホストで`bin / TimeRegionSize`と`bin % TimeRegionSize`に分けて、
HBF番号とHBF内のTrgTimeへ戻す。ヒットのないHBFの空配列も維持する。

## 順序とバッファの寿命

Clear・各コピー・Mark・サブグループ集約・LUT・CUB scan・compactionは
すべて同じデフォルトストリームへ投入する。投入順序でGPU内の依存関係を守るので、
Markごとの`cudaDeviceSynchronize()`はない。
件数の同期版D2Hコピーで先行処理の完了を待つ。ExecBatchが戻る時には入力のコピーも
終わっているので、次のHBFバッチでGPUバッファやTimeRegionを再利用できる。
元メッセージとホスト記述子も、その完了までは維持される。
初期化時のSetTimeRegionには従来どおり明示的な同期がある。

## 不要コードの削除

GPU版から、CPU用のチャンネルstd::map群、Entry/EntryTo/ClearEntry/CheckEntryFEM、
Mark長のsetter/getter、InitParam、ホストTimeRegionと読み戻しgetter、未使用ホスト配列、
旧seg単位Mark経路、古いTEST_MAIN_TRIGを削除した。Scanと設定内部用helperはprivateにした。
HBFIndexの未使用メッセージ番号・TF番号、DataBlockの読み出されないフィールド、
未使用のMakeBlockMap2/CheckData/RemoveData、無効化済みコードと専用timerも削除した。
初期化中だけ必要なデコーダー・LUT生成器・登録用vectorはローカル変数にした。
既存のヘッダー生成・データ抑制・DQM出力・運用ログに使う値は残している。

## 確認

```sh
cmake --build build --target LogicFilter -j4
/usr/local/cuda/bin/nvcc -std=c++17 -arch=sm_86 -rdc=true -I. -ILogicFilter/GPU \
  -I"$HOME/local/nestdaq/include" LogicFilter/tests/GpuBatchedMarkTest.cu \
  "$HOME/local/nestdaq/lib/libcudamap.a" -o /tmp/gpu-batched-mark-test
/tmp/gpu-batched-mark-test
```

4 TDC形式、39入力×4 HBF、異なるIP・形式・長さ、複数block、空入力、未登録IP、
32信号サブグループの全成立／1信号欠損を検査する。不要コード削除後にビルドとテストは通過。
全run000400.datの338701 TrgTime一致は削除前の処理で確認済みであり、今回再実行はしていない。

### 初期化処理の分担

`SetTimeRegion(size, nHbf)`はnHbf個分のTimeRegion本体の確保とゼロ初期化を行う。
続けて`InitializeScanBuffers()`を呼び、立ち上がりフラグ、exclusive-scanの
オフセットと作業領域、圧縮したヒットと個数の格納領域を準備する。
TimeRegionのサイズを変えたときも`InitializeScanBuffers()`を呼び直す。
初期化されていないサイズでのScanは例外にする。

`CleanUpSubTimeRegion()`は引数を取らず、初期化時に転送した各SubTimeRegionの
信号数を使ってビットを復元する。cleanupごとの設定のH2Dコピーは行わない。

### 4 HBFバッチの確認

`LogicFilter/tests/GpuHbfBatchTest.cu`で1 HBFずつの結果と比較する。
4 TDC形式、39基板×4 HBF、HBF境界、空HBF、SubTimeRegionの独立性、
処理HBF数4→2→1→3→4での領域再利用、奇数サイズと最大ヒット密度を確認する。
初期化時に1〜nHbf個のscanに必要な作業領域サイズを問い合わせ、最大値で一度だけ確保する。
実行中は確保済み領域を使い回す。

```sh
/usr/local/cuda/bin/nvcc -std=c++17 -arch=sm_86 -rdc=true -I. -ILogicFilter/GPU \
  -I"$HOME/local/nestdaq/include" LogicFilter/tests/GpuHbfBatchTest.cu \
  "$HOME/local/nestdaq/lib/libcudamap.a" -o /tmp/gpu-hbf-batch-test
/tmp/gpu-hbf-batch-test
```

`hbf-batch-size`はNestDAQのparameter_configに合わせ文字列オプションとし、1〜4を許可する。
既存のelapsedファイル名は維持しているため、1/4を別runで比較するときは各STOP後に退避する。

取り込み後の確認：run000400.datの先頭1024 TFをwait=1、独立IPCで再生した。
1 HBF/4 HBFで1024 MESSAGE、4096 META、7879 TrgTime、7879 SLICEが完全一致し、
全デバイスが正常終了した。結果は`/tmp/logicfilter-hbf-port-{1,4}`に保存している。
