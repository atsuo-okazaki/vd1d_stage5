# use_irradiation_delay 機能

## 概要

到着時刻 (arrival time) モデル: 各時間ステップで
1. 降着量 mdot(t) → L_acc(t) を計算
2. その照射がディスクに届く時刻 t_arrive = t + τ を計算
3. (t_arrive, L_acc) をバッファに記録
4. 現在時刻 t での L_irr は「t に到着する照射」を補間して取得

## 物理

- 時刻 t で発生した L_acc は τ 後に届く
- L_irr(t_target) = 時刻 t_target に到着する照射（補間）
- t_target < 最早到着時刻 → 0（まだ何も届いていない）
- τ の与え方: tau_irr_lag_mode で選択

## 遅延時間の与え方

| tau_irr_lag_mode | 遅延 τ の計算 |
|------------------|---------------|
| **explicit** | τ = tau_irr_lag_nd × t0 （固定値） |
| **viscous**  | τ = tau_irr_lag_nd × (r_in² / ν_in) （内縁の粘性時間） |

- explicit: 光伝播遅延などの簡易モデル
- viscous: 降着流の粘性時間スケール（時間とともに変化）

## 設定（ad1d.in / disk_mode）

```namelist
use_irradiation_delay = .true.
tau_irr_lag_mode     = 'explicit'   ! または 'viscous'
tau_irr_lag_nd       = 100.0        ! explicit: 100×t0; viscous: f_delay
```

## 動作

1. 各ステップで (t_arrive, L_acc) をリングバッファに記録（将来の到着のみ）
2. L_irr(t) = 時刻 t に到着する照射を sample_Lacc_at(t) で補間取得
3. チェックポイント (CKPT_VER=4) でバッファを保存・復元

## 設計上の利点

- **バッファ長**: 将来の到着時刻だけを覚えておけばよいため、nt を大きくする必要はない
- **到着の前後**: 粘性遅延などで delay が変動すると、発生順と到着順が逆転し得る。`sample_Lacc_at` 内で到着時刻順にソートしてから補間し、前後する到着にも対応

## 注意

- use_be_decretion = .true. のときは L_star をそのまま使い、遅延モデルは適用しない
- 現状のバッファ長は nt（init_irradiation_buffer(nt)）。設計上は max_delay をカバーする長さで十分

---

## use_finite_irradiation_source（有限サイズ照射源）

LOH24 Eq.16 は点光源を仮定（Appendix A）。恒星の有限半径 R_star を考慮するには、
grazing angle を加法的に扱う（Chiang & Goldreich; LOH24 の表面形状と整合）。

- **use_finite_irradiation_source = .true.** のとき:
  - α_total = α_star + α_flare
  - α_star = 0.4 × R_star / r（有限恒星円盤）
  - α_flare = ξ × dY/dξ = r × d(H/r)/dr（flaring、LOH24 に既存）
  - α_flare ≥ α_star: Qirr = Qirr_LOH24 × (α_star+α_flare)/α_flare
  - α_flare < α_star: Qirr = (α_star+α_flare)/α_star × Qirr_flat_star（滑らか、発散なし）
- Be disk のように恒星照射が支配的な場合は .true. 推奨
- 影（shadow）の扱いは別途検討
