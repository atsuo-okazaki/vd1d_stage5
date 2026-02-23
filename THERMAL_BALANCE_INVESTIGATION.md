# 熱収支 (Qvis+Qirr vs Qrad) 不整合の調査結果

## 問題の症状
最近の変更により、Qvis+Qirr と Qrad が大きく食い違ったまま計算が進行するようになった。

## 調査で特定した原因

### 1. Qvis_and_Qrad_KFM2008 の Hubeny 型ブリッジ公式の係数誤り (主要因)

**disk_thermal_mod.f90** で、`Qvis_and_Qrad_KFM2008` の放射冷却 Qrad のブリッジ公式が変更された:

```fortran
! 現在の実装 (誤り)
Qrad_loc = (16.0_dp * sbc * T_loc**4) / &
           ( 3.0_dp * tauR_eff + 2.0_dp / max(tauP_eff, 1.0e-30_dp) )
```

**極限の検証:**
- **光学厚 (tauR >> 1)**: Qrad → 16/(3*tauR) × sbc*T^4 = (32/3)×sbc*T^4/(κR*Σ)
  - 正しい極限: **(64/3)×sbc*T^4/(κR*Σ)** → 係数が **1/2** になっている
- **光学薄 (tauP << 1)**: Qrad → 16×tauP/2 × sbc*T^4 = 4×κP×Σ×sbc*T^4
  - 正しい極限: **2×κP×Σ×sbc*T^4** → 係数が **2倍** になっている

両極限で系統的な誤差があるため、根探索で得た T が実際の熱平衡とずれ、見かけ上 Qvis+Qirr ≠ Qrad のように見える。

**正しい係数:**
```
Qrad = 32×sbc×T^4 / (3×tauR + 8/tauP)
```
- 光学厚: 32/(3*tauR) = 64/(3×κR×Σ) ✓
- 光学薄: 32×tauP/8 = 2×κP×Σ ✓

### 2. output_mod での Qvis/Qrad の再計算式の不一致

**output_mod.f90** では `Qvis_and_Qrad`（別ルーチン）で Qvis, Qrad を再計算している。一方、熱平衡求解と `build_Qvis_Qrad_from_solution` は `heating_cooling_cell` → `Qvis_and_Qrad_KFM2008` を使用している。

- `Qvis_and_Qrad`: 調和ブレンド w=τR/(1+τR)
- `Qvis_and_Qrad_KFM2008`: Hubeny 型 (上記式)

式が異なるため、出力時に表示する Qvis/Qrad が、求解に使った値と一致しない可能性がある。表示の一貫性のためには、いずれかに統一することが望ましい。

### 3. enforce_equilibrium_by_bisection 失敗時の整合性

`thermal_solve_stageA` で、`enforce_equilibrium_by_bisection` が ierr≠0 で終了した場合、T が更新されても H, rho などは再計算されず不整合が残っていた。（→ 解決済み：bisection 後に常に `heating_cooling_cell` を呼ぶよう変更）

### 4. irradiation_mod の有限恒星補正（加法型）

`use_finite_irradiation_source` のとき、Qirr の算出式は検討の結果、加法型 `Qirr = pref × (α_flare + α_star) / ξ^2` に変更されている。旧式のスケーリング型補正から加法型への変更は意図的な設計判断である。

## 推奨修正

### 最優先: Qrad の係数修正 (disk_thermal_mod.f90)

```fortran
! 修正後
Qrad_loc = (32.0_dp * sbc * T_loc**4) / &
           ( 3.0_dp * tauR_eff + 8.0_dp / max(tauP_eff, 1.0e-30_dp) )
```

### オプション: output_mod の一貫性

出力も `heating_cooling_cell` 経由で計算するか、あるいは保存済みの `Qvis(it,i)`, `Qrad(it,i)` をそのまま使うと、求解と表示が一致する。（→ 実施済み：保存済み値を使用するように変更済み）
