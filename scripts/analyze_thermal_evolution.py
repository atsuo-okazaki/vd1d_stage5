#!/usr/bin/env python3
"""
Thermal evolution analysis for 1D accretion disk simulation.

Post-processing script (zero runtime cost during simulation):
  - Limit cycle detection from mdot_inner.dat
  - (Sigma, T) evolutionary track on S-curve
  - Thermal stability evolution from thermal_stability_*.dat

Requirements: numpy (matplotlib optional for plotting)

Usage:
  python analyze_thermal_evolution.py limit_cycle              # detect outbursts
  python analyze_thermal_evolution.py track r_cgs [scurve.dat]  # Sigma-T track at r
  python analyze_thermal_evolution.py stability                 # n_unstable vs t
"""

import argparse
import glob
import re
import sys
from pathlib import Path

import numpy as np

# Optional matplotlib for plotting
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_PLT = True
except ImportError:
    HAS_PLT = False


def load_mdot_inner(fname='mdot_inner.dat'):
    """Load mdot_inner.dat. Returns t_nd, t_s, mdot_nd, mdot_phys."""
    data = []
    with open(fname) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 4:
                data.append([float(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])])
    if not data:
        return None, None, None, None
    arr = np.array(data)
    return arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]


def load_disk_file(fname, cols=None):
    """Load disk_t*.dat. cols: dict of name->index (1-based). Returns dict of arrays."""
    data = []
    with open(fname) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 12:
                data.append([float(x) for x in parts[:15]])
    if not data:
        return {}
    arr = np.array(data)
    # Default: xi=0, r=1, sigma=2, ..., Tmid=11 (0-based: 10)
    default_cols = {'xi': 0, 'r': 1, 'sigma': 2, 'H': 3, 'Tmid': 11}  # disk_t*.dat: 15 cols
    c = cols or default_cols
    return {k: arr[:, v] for k, v in c.items()}


def load_thermal_stability(fname):
    """Load thermal_stability_*.dat. Returns r, Sigma, T, stability, and summary if present."""
    data = []
    summary = {}
    with open(fname) as f:
        for line in f:
            if line.startswith('#'):
                if 'Summary' in line or 'n_unstable' in line:
                    m = re.search(r'n_unstable\s*=\s*(\d+)', line)
                    if m:
                        summary['n_unstable'] = int(m.group(1))
                    m = re.search(r'r_unstable\s*\[cm\]\s*=\s*\[\s*([\d.eE+-]+)\s*,\s*([\d.eE+-]+)\s*\]', line)
                    if m:
                        summary['r_min'] = float(m.group(1))
                        summary['r_max'] = float(m.group(2))
                continue
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 10:
                data.append([float(parts[0]), float(parts[1]), float(parts[2]),
                            float(parts[3]), float(parts[4]), float(parts[5]),
                            float(parts[6]), float(parts[7]), float(parts[8]),
                            int(parts[9])])
    if not data:
        return None, None, None, None, summary
    arr = np.array(data)
    return arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 9], summary


def detect_limit_cycle(mdot_phys, t_s, threshold_ratio=2.0, min_sep=10):
    """
    Detect outbursts (limit cycle) from mdot time series.
    threshold_ratio: outburst = mdot > median * threshold_ratio
    min_sep: minimum steps between outbursts
    Returns list of (t_start, t_peak, t_end, mdot_peak) for each outburst.
    """
    if mdot_phys is None or len(mdot_phys) < 3:
        return []
    med = np.median(mdot_phys)
    if med <= 0:
        return []
    thresh = med * threshold_ratio
    in_outburst = mdot_phys > thresh
    # Find runs of True
    bursts = []
    i = 0
    while i < len(in_outburst):
        if in_outburst[i]:
            start = i
            peak = i
            peak_val = mdot_phys[i]
            while i < len(in_outburst) and in_outburst[i]:
                if mdot_phys[i] > peak_val:
                    peak = i
                    peak_val = mdot_phys[i]
                i += 1
            end = i
            if end - start >= 1:  # at least 2 points
                bursts.append((t_s[start], t_s[peak], t_s[end-1], peak_val))
            # Require min_sep before next burst
            i = min(i + min_sep, len(in_outburst))
        else:
            i += 1
    return bursts


def cmd_limit_cycle(args):
    """Detect and report limit cycle (outbursts)."""
    t_nd, t_s, mdot_nd, mdot_phys = load_mdot_inner(args.mdot_file)
    if mdot_phys is None:
        print(f"ERROR: No data in {args.mdot_file}")
        return 1
    bursts = detect_limit_cycle(mdot_phys, t_s,
                              threshold_ratio=args.threshold,
                              min_sep=args.min_sep)
    print(f"Limit cycle analysis: {args.mdot_file}")
    print(f"  Total steps: {len(mdot_phys)}, time range: {t_s[0]:.3e} - {t_s[-1]:.3e} [s]")
    print(f"  mdot range: {mdot_phys.min():.3e} - {mdot_phys.max():.3e} [g/s]")
    print(f"  Detected {len(bursts)} outburst(s) (threshold = median x {args.threshold})")
    for i, (t0, tpeak, t1, mpeak) in enumerate(bursts, 1):
        print(f"    Outburst {i}: t=[{t0:.3e}, {t1:.3e}] s, peak at t={tpeak:.3e} s, mdot_peak={mpeak:.3e} g/s")
    if args.output and HAS_PLT:
        plt.figure(figsize=(8, 5))
        plt.semilogy(t_s, mdot_phys, 'b-', lw=0.8, label='Mdot')
        for t0, tpeak, t1, mpeak in bursts:
            plt.axvspan(t0, t1, alpha=0.2, color='red')
        plt.xlabel('t [s]')
        plt.ylabel('Mdot [g/s]')
        plt.title('Limit cycle detection')
        plt.legend()
        plt.tight_layout()
        plt.savefig(args.output, dpi=150)
        print(f"  Plot saved: {args.output}")
    return 0


def cmd_track(args):
    """Extract (Sigma, T) evolutionary track at r_cgs from disk_t*.dat."""
    files = sorted(glob.glob('disk_t*.dat'))
    if not files:
        print("ERROR: No disk_t*.dat files found")
        return 1
    r_target = float(args.r_cgs)
    # Extract it from filename: disk_t00000100.dat -> 100
    pattern = re.compile(r'disk_t(\d+)\.dat')
    tracks = []
    for f in files:
        m = pattern.match(Path(f).name)
        if not m:
            continue
        it = int(m.group(1))
        d = load_disk_file(f, {'r': 1, 'sigma': 2, 'Tmid': 10})
        if not d:
            continue
        r, sigma, T = d['r'], d['sigma'], d['Tmid']
        # Interpolate at r_target
        if r_target <= r[0]:
            s, t = sigma[0], T[0]
        elif r_target >= r[-1]:
            s, t = sigma[-1], T[-1]
        else:
            idx = np.searchsorted(r, r_target)
            w = (r_target - r[idx-1]) / (r[idx] - r[idx-1])
            s = sigma[idx-1] + w * (sigma[idx] - sigma[idx-1])
            t = T[idx-1] + w * (T[idx] - T[idx-1])
        tracks.append((it, s, t))
    if not tracks:
        print("ERROR: No valid data in disk files")
        return 1
    tracks = np.array(tracks)
    np.savetxt(args.output or 'track_r{:.0e}.dat'.format(r_target),
               tracks, header='it  Sigma[g/cm2]  T[K]  (at r={:.2e} cm)'.format(r_target),
               fmt='%d %.6e %.6e')
    print(f"Track at r={r_target:.2e} cm: {len(tracks)} points -> {args.output or 'track_r*.dat'}")
    if args.scurve and Path(args.scurve).exists() and HAS_PLT:
        scurve = np.loadtxt(args.scurve)
        if len(scurve.shape) == 2 and scurve.shape[1] >= 2:
            plt.figure(figsize=(8, 6))
            plt.plot(scurve[:, 1], scurve[:, 0], 'k-', lw=1.5, label='S-curve')
            sc = plt.scatter(tracks[:, 2], tracks[:, 1], c=range(len(tracks)), cmap='viridis', s=10)
            plt.colorbar(sc, label='time index')
            plt.xlabel('T [K]')
            plt.ylabel('Sigma [g/cm^2]')
            plt.xscale('log')
            plt.yscale('log')
            plt.title(f'Evolutionary track on S-curve at r={r_target:.2e} cm')
            plt.legend()
            plt.tight_layout()
            outfig = (args.output or 'track_r{:.0e}.dat'.format(r_target)).replace('.dat', '.pdf')
            plt.savefig(outfig, dpi=150)
            print(f"  Plot saved: {outfig}")
    return 0


def cmd_stability(args):
    """Summarize thermal stability evolution from thermal_stability_*.dat."""
    files = sorted(glob.glob('thermal_stability_*.dat'))
    # Exclude all_roots
    files = [f for f in files if 'all_roots' not in f]
    if not files:
        print("ERROR: No thermal_stability_*.dat files found")
        return 1
    pattern = re.compile(r'thermal_stability_(\d+)\.dat')
    results = []
    for f in files:
        m = pattern.search(f)
        if not m:
            continue
        it = int(m.group(1))
        r, sigma, T, stab, summary = load_thermal_stability(f)
        n_unst = summary.get('n_unstable', np.sum(stab == -1) if stab is not None else 0)
        # Read t from first # line if possible
        t_nd = np.nan
        with open(f) as fp:
            for line in fp:
                if 't_nd' in line:
                    mt = re.search(r't_nd\s*=\s*([\d.e+-]+)', line)
                    if mt:
                        t_nd = float(mt.group(1))
                    break
        results.append((it, t_nd, n_unst))
    if not results:
        print("ERROR: No valid thermal stability data")
        return 1
    res = np.array(results)
    out = args.output or 'thermal_stability_summary.dat'
    np.savetxt(out, res, header='it  t_nd  n_unstable', fmt='%d %.6e %d')
    print(f"Stability summary: {len(results)} snapshots -> {out}")
    if HAS_PLT:
        plt.figure(figsize=(8, 4))
        valid = ~np.isnan(res[:, 1])
        if np.any(valid):
            plt.plot(res[valid, 1], res[valid, 2], 'b-o', ms=3)
            plt.xlabel('t (dimensionless)')
        else:
            plt.plot(res[:, 0], res[:, 2], 'b-o', ms=3)
            plt.xlabel('it')
        plt.ylabel('n_unstable')
        plt.title('Thermal instability evolution')
        plt.tight_layout()
        plt.savefig(out.replace('.dat', '.pdf'), dpi=150)
        print(f"  Plot saved: {out.replace('.dat', '.pdf')}")
    return 0


def main():
    parser = argparse.ArgumentParser(description='Thermal evolution analysis')
    sub = parser.add_subparsers(dest='cmd', required=True)
    # limit_cycle
    p1 = sub.add_parser('limit_cycle', help='Detect outbursts from mdot_inner.dat')
    p1.add_argument('-f', '--mdot-file', default='mdot_inner.dat')
    p1.add_argument('-t', '--threshold', type=float, default=2.0, help='Outburst threshold (x median)')
    p1.add_argument('-s', '--min-sep', type=int, default=10, help='Min steps between bursts')
    p1.add_argument('-o', '--output', help='Output plot file')
    # track
    p2 = sub.add_parser('track', help='Sigma-T track at r from disk_t*.dat')
    p2.add_argument('r_cgs', help='Radius [cm]')
    p2.add_argument('scurve', nargs='?', help='Optional: scurve_*.dat for overlay')
    p2.add_argument('-o', '--output', help='Output data file')
    # stability
    p3 = sub.add_parser('stability', help='n_unstable vs t from thermal_stability_*.dat')
    p3.add_argument('-o', '--output', help='Output summary file')
    args = parser.parse_args()
    if args.cmd == 'limit_cycle':
        return cmd_limit_cycle(args)
    if args.cmd == 'track':
        return cmd_track(args)
    if args.cmd == 'stability':
        return cmd_stability(args)
    return 1


if __name__ == '__main__':
    sys.exit(main())
