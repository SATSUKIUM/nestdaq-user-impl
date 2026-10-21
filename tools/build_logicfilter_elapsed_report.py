#!/usr/bin/env python3
"""Build a self-contained HTML report from the simple elapsed comparison plots."""
import base64
import csv
import html
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TEST = ROOT / "test"
PREFIX = "logicfilter-elapsed-comparison"


def embedded(path):
    mime = "image/png" if path.suffix == ".png" else "application/pdf"
    return f"data:{mime};base64," + base64.b64encode(path.read_bytes()).decode("ascii")


def main():
    with (TEST / f"{PREFIX}-summary.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    names = {"cpu": "CPU", "hbf1": "GPU・1 HBFずつ", "hbf4": "GPU・4 HBFまとめて"}
    cpu_mean = float(next(row["mean_us"] for row in rows if row["mode"] == "cpu"))
    table = []
    tails = []
    for row in rows:
        mode = row["mode"]
        values = [float(line) for line in (TEST / f"logicfilter-elapsed-0-{mode}.txt").read_text().splitlines()
                  if line.strip() and not line.lstrip().startswith("#")]
        assert len(values) == int(row["samples"]), "Summary and input sample counts differ"
        table.append(f'<tr><th scope="row">{names[mode]}</th><td>{int(row["samples"]):,}</td>'
                     + "".join(f'<td>{float(row[key]):,.2f}</td>' for key in ("mean_us", "median_us", "p95_us", "p99_us"))
                     + f'<td>{cpu_mean / float(row["mean_us"]):.2f}倍</td></tr>')
        outside = sum(x < 0 or x >= 8000 for x in values)
        tails.append(f'{names[mode]}：範囲外{outside:,}件、最大{max(values):,.2f} µs')
    figures = []
    for suffix, title, note in [
        ("-full", "CPU・GPUの分布比較", "横軸0〜8,000 µs、縦軸は全サンプルに対する割合（% / bin）。両軸とも線形です。"),
        ("-zoom", "GPU付近の拡大図", "横軸400〜1,200 µs。CPUの分布は表示範囲外です。平均値は全サンプルから計算しています。"),
    ]:
        stem = PREFIX + suffix
        figures.append(f'''<section><h2>{title}</h2><p>{note}</p>
        <img src="{embedded(TEST / (stem + '.png'))}" alt="{title}" width="800" height="600">
        <p class="downloads"><a download="{stem}.pdf" href="{embedded(TEST / (stem + '.pdf'))}">PDFを保存</a>
        <a download="{stem}.png" href="{embedded(TEST / (stem + '.png'))}">PNGを保存</a></p></section>''')
    source_files = "、".join(f'<code>logicfilter-elapsed-0-{mode}.txt</code>' for mode in names)
    document = '''<!doctype html><html lang="ja"><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>LogicFilter 処理時間の比較</title>
<style>
*{box-sizing:border-box}body{margin:0;background:#f4f6f8;color:#202934;font:16px/1.75 system-ui,sans-serif}
main{max-width:1040px;margin:auto;padding:32px 20px 60px}h1{font-size:28px;line-height:1.4}h2{font-size:21px;margin:0 0 12px}
section{background:white;border:1px solid #dbe1e8;border-radius:10px;padding:24px;margin:24px 0}
p{margin:10px 0}img{display:block;width:100%;max-width:800px;height:auto;margin:16px auto}
.scroll{overflow-x:auto}table{width:100%;border-collapse:collapse;white-space:nowrap;font-variant-numeric:tabular-nums}
th,td{padding:10px 12px;border-bottom:1px solid #dbe1e8;text-align:right}th:first-child{text-align:left}thead{background:#f4f6f8}
a{color:#155ca2}.downloads{display:flex;gap:24px}code{overflow-wrap:anywhere;font-size:13px}.small{font-size:14px;color:#4c5968}
@media print{body{background:white}main{padding:0}section{break-inside:avoid}.downloads{display:none}}
</style><main><h1>LogicFilter 処理時間の比較</h1>
<p>Tracking.shで測定したCPU版、GPU版（1 HBFずつ／4 HBFまとめて）のConditionalRunの分布です。</p>
<section><h2>計測条件と平均値</h2>
<p>受信成功後からSend完了までを計測し、実行中はvectorへ蓄積、STOP時にテキストへ保存しました。
Receiveの待ち時間を除き、Sendの下流待ち時間を含みます。1サンプルはConditionalRun 1回（通常4 HBFを含むTF）の時間です。</p>
<div class="scroll"><table><thead><tr><th>処理方式</th><th>サンプル数</th><th>平均 [µs]</th><th>中央値 [µs]</th><th>95%点 [µs]</th><th>99%点 [µs]</th><th>CPU比の速度</th></tr></thead><tbody>'''
    document += "".join(table) + '''</tbody></table></div>
<p class="small">ウォームアップ・外れ値の除外は行っていません。ヒストグラムは200 bins、線幅3、TCanvasは800×600です。</p></section>'''
    document += "".join(figures)
    document += '<section><h2>表示範囲外のサンプル</h2><p>0〜8,000 µsの図では以下のサンプルが範囲外になります。平均値・分位点にはすべて含めています。</p><ul>'
    document += "".join(f"<li>{html.escape(t)}</li>" for t in tails)
    document += f'</ul><p class="small">入力：{source_files}</p></section><p class="small">画像とダウンロード用PDFはこのHTMLに埋め込み済みです。HTML単体をコピーして閲覧できます。</p></main></html>'
    output = TEST / f"{PREFIX}.html"
    output.write_text(document)
    print(output)
    print(f"{output.stat().st_size:,} bytes; 2 embedded plots with PDF and PNG downloads")


if __name__ == "__main__":
    main()
