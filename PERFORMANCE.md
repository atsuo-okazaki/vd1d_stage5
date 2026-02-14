# パフォーマンス最適化ガイド

## 実施済みの変更

### 1. チェックポイント完成（CKPT_VER=3）
- Qvis, Qrad, is_shadow を保存・復元
- 照射あり再開時の一貫性を確保

### 2. Makefile リリース設定
- **-fcheck=all を削除**: 配列境界チェックで 2〜5 倍程度遅延していたため削除
- **-march=native**: CPU 最適化（M1 では NEON 等を有効化）
- デバッグビルド: `make debug` で `-fcheck=all` 付きビルド

### 3. 10 倍動径範囲
- `ad1d_wide.in`: rout_value=1000, rin_value=1, nr=401
- ロググリッドにより 1–100 と同等の解像度を維持

## 出力頻度の切り分け（run_control.nml）

| パラメータ | 用途 | 目安 | 備考 |
|------------|------|------|------|
| **outfreq** | disk_t*.dat（円盤構造・アニメ用） | 頻繁 | 動画用のため現状維持 |
| **chkfreq** | チェックポイント | 疎らでOK（最低10回） | 既定: nt/10 |
| **thermal_stability_freq** | 熱安定性解析出力 | 疎らでOK（最低10回） | 既定: nt/10 |

例（t_sim_end = 3e7, nt ≈ 20000 のとき）:
- outfreq = 0 → 自動で約 20 ステップごと（disk 出力）
- chkfreq = 0 → 自動で約 2000 ステップごと（checkpoint、10回）
- thermal_stability_freq = 0 → 自動で約 2000 ステップごと（10回）

## さらなる高速化のヒント

### OpenMP スレッド数
```bash
export OMP_NUM_THREADS=10   # M1 Max など（実コア数に合わせる）
./vd1d_stage5 ad1d.in
```

### 出力頻度の明示指定
- `outfreq = 20`: disk 構造の出力間隔（アニメ用、nt=20000 なら約20で約1000枚）
- `chkfreq = 2000`: checkpoint（nt=20000 なら10回）
- `thermal_stability_freq = 2000`: 熱安定性出力（10回）

### -march=native が失敗する場合
古い gfortran では `-march=native` が使えない場合があります。
Makefile の FFLAGS から `-march=native` を削除してください。
