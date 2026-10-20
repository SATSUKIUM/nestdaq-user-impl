# Tracking.shでのelapsed計測

`verify-logicfilter-cpu-gpu-trgtime`のGPU LogicFilterは、受信成功後から
Send終了までの時間をマイクロ秒でvectorに蓄積し、PostRunでテキストへ保存する。
ROOTへの実行中の書き込み、Exec詳細計測、CUDAイベント計測は使用しない。
メモリ検査は独立したテストコマンドであり、通常runには含まれない。
1サンプルはConditionalRun 1回（通常4 HBFを含むTF）であり、1 HBFの時間ではない。

## 起動設定

`~/GitHub/nestdaq-run/mq-param/mq-param_TFBFP-LF-TFS-Nothing.sh`の有効な設定を以下にする。
同じ項目が複数あれば後の設定が優先されるため重複させない。

```sh
param TFBFilePlayer-0 in-file ./data/run000400.dat wait 1
param LogicFilter hbf-batch-size 4
```

`max-iterations`は指定しない。比較対象は`hbf-batch-size 1`とする。
Tracking.shはmq-paramをsourceするときRedisのパラメータDBを初期化するため、
起動前にredis-cliで手動設定しても消される。
`performance-output`、`exec-details`、`exec-cuda-events`はGPUの起動設定から削除する。

start_device.shは`$NESTDAQUSERIMPL/bin/LogicFilter`があれば優先し、
なければ`$NESTDAQ/bin/LogicFilter`を使う。
更新したバイナリが選択されることを確認してから、通常通り起動・RUNする。

```sh
cd ~/GitHub/nestdaq-run
bash launcher/Tracking.sh
```

ファイル全体の処理と下流の処理が完了したら、通常のSTOP操作を行う。
STOP前に次のlauncherを起動するとtmuxセッションが終了されるため、保存完了ログを確認する。

## 出力

起動ディレクトリによらずビルド時のソースツリー直下のtestに出力する。
LogicFilter-0の場合は以下。

- `test/logicfilter-elapsed-0-hbf1.txt`
- `test/logicfilter-elapsed-0-hbf4.txt`

同じID・batch数で再測定すると上書きされる。既存の`logicfilter-elapsed-0.txt`は保持する。
Receiveは除外するが、Sendの下流待ちは含むため、両runのトポロジーと入力waitを揃える。
以前の約665 usとの比較もConditionalRun単位で行う。

```sh
cd ~/GitHub/nestdaq-user-impl
root -l -b -q 'test/plot_logicfilter_elapsed.C("test/logicfilter-elapsed-0-hbf1.txt")'
root -l -b -q 'test/plot_logicfilter_elapsed.C("test/logicfilter-elapsed-0-hbf4.txt")'
```

800x600、200 bins、line width 3でPNG/PDF/ROOTを出力する。
以前のROOTによる詳細計測・箱庭比較は`6379b07`時点の実装と結果を参照。
