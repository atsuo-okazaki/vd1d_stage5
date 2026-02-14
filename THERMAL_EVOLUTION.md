# 熱的不安定性の時間進化解析

## 概要

S カーブ解析を拡張し、シミュレーション時間発展に沿った熱的不安定性の診断を行います。

## 1. Fortran 側の変更（軽量）

### thermal_stability_*.dat へのサマリー追加

`output_thermal_stability_profile` 出力の末尾に、不安定セル数と半径範囲を追加：

```
# Summary: n_unstable = 42, r_unstable [cm] = [ 1.23e12, 5.67e13]
```

- 既存の per-cell 出力（dQ+/dT, dQ-/dT, stability）は変更なし
- 計算コストは thermal_stability_freq（既定 nt/10）のタイミングのみ

## 2. Python 事後解析スクリプト（ゼロ計算コスト）

`scripts/analyze_thermal_evolution.py`

**必要環境**: numpy（matplotlib はプロット用、オプション）

### サブコマンド

| コマンド | 機能 |
|----------|------|
| `limit_cycle` | mdot_inner.dat からアウトバーストを検出 |
| `track r_cgs [scurve.dat]` | 指定半径での (Σ, T) 進化経路を抽出。S カーブ重ね書き可 |
| `stability` | thermal_stability_*.dat から n_unstable の時間変化を集約 |

### 使用例

```bash
# アウトバースト検出（閾値 = 中央値の2倍）
python scripts/analyze_thermal_evolution.py limit_cycle -o mdot_bursts.pdf

# r=1e12 cm での (Σ, T) 経路、S カーブ上にプロット
# 事前に: ./scurve で scurve.dat を生成
python scripts/analyze_thermal_evolution.py track 1e12 scurve.dat -o track.dat

# 熱的不安定セル数の時間発展
python scripts/analyze_thermal_evolution.py stability -o stability_summary.dat
```

### limit_cycle オプション

- `-t, --threshold`: アウトバースト判定閾値（中央値の倍数、既定 2.0）
- `-s, --min-sep`: 連続アウトバースト間の最小ステップ数（既定 10）

## 3. ワークフロー

1. 通常どおり `./vd1d_stage5 ad1d.in` を実行
2. `do_thermal_stability_output = .true.` で thermal_stability_*.dat を出力
3. 計算終了後、`analyze_thermal_evolution.py` で事後解析
4. S カーブ上への経路プロット: 希望する半径で `./scurve` を実行し、`track` コマンドに `scurve.dat` を指定
