# チェックポイント機能の分析レポート

## 現状：保存されている変数 (checkpoint_mod.f90)

| 変数 | it-1 | it | 備考 |
|------|------|-----|------|
| sigmat | ✓ | ✓ | |
| nu_conv | ✓ | ✓ | |
| Tmid | ✓ | ✓ | |
| H | ✓ | ✓ | |
| rho | ✓ | ✓ | |
| kappaR | ✓ | ✓ | |
| tauR | ✓ | ✓ | |
| Qirr | ✓ | ✓ | |
| dYdXi | ✓ | ✓ | |
| Qvis | ✓ | ✓ | V3 |
| Qrad | ✓ | ✓ | V3 |
| is_shadow | ✓ | ✓ | V3 |

## 実装済み（V3）

- **Qvis, Qrad, is_shadow**: checkpoint_mod.f90 CKPT_VER=3 で保存・復元済み
- **state_io_mod.f90**: shadow_cur(:) = is_shadow(it,:) で復元

## 不足していた変数（V2 以前）

### 1. **Qvis, Qrad**（重要度：中）→ V3で対応済み

- **旧状態**: 保存されていなかった
- **影響**: `load_state_from_history` が `Qvis_cur`, `Qrad_cur` を `Qvis(it,:)`, `Qrad(it,:)` から読み込むが、restart 時は 0 または未初期化のまま
- **evolve の挙動**: `Qvis_prev`, `Qrad_prev` として反復の初期値に使われる。0 だと最初のステップで収束がやや遅れる可能性あり
- **output**: `output_disk_structure_dtl` は `Qvis_and_Qrad` で再計算しているため、出力自体は可能

### 2. **is_shadow**（重要度：高）→ V3で対応済み

- **旧状態**: 保存されていなかった
- **影響**: `load_state_from_history` で `shadow_cur(:) = .false.` に固定されている（`is_shadow(it,:)` の読み込みがコメントアウト）
- **問題**: `use_irradiation = .true.` の場合、shadow が disk の幾何 (H/r) から計算される。restart 時は shadow が常に false となり、**照射が過大評価される**
- **evolve**: `shadow_cur` は thermal 反復で重要な役割を持つ

### 3. **k_iter, m_iter**（重要度：低）

- **現状**: 保存されていない
- **影響**: 反復回数の診断用。restart 後は 0 になるが、物理計算には影響しない

### 4. **r_edge, i_edge**（重要度：低、use_be_decretion 時のみ）

- **現状**: 保存されていない
- **影響**: `use_be_decretion = .true.` の場合、外縁追跡に使われる。restart 直後は未設定になる可能性

### 5. **irradiation バッファ**（重要度：中）→ V4で対応済み

- **旧状態**: t_hist, mdot_hist はチェックポイントに含まれていなかった
- **V4**: use_irradiation_delay 時に save_irradiation_buffer/load_irradiation_buffer で保存・復元

## 推奨する修正

### 必須（Qvis, Qrad, is_shadow の追加）

1. **checkpoint_mod.f90** の write/read に以下を追加：
   - `Qvis(it_chk-1,:)`, `Qvis(it_chk,:)`
   - `Qrad(it_chk-1,:)`, `Qrad(it_chk,:)`
   - `is_shadow(it_chk-1,:)`, `is_shadow(it_chk,:)`（logical の I/O に注意）

2. **checkpoint バージョン**: 既存チェックポイントとの互換性のため、`CKPT_VER` を 3 に上げる

3. **state_io_mod.f90**: `shadow_cur(:) = is_shadow(it,:)` のコメントアウトを解除（is_shadow を checkpoint に保存すれば有効化可能）

### オプション（r_edge, i_edge）

- `use_be_decretion` 時に `r_edge(it_restart)`, `i_edge(it_restart)` を保存・復元

### オプション（irradiation バッファ）

- `use_irradiation_delay = .true.` を使用する場合、`irradiation_mod` のバッファ保存・復元用サブルーチンを追加
