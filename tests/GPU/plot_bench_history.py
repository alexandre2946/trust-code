#!/usr/bin/env python3
"""
Plot evolution of 'Average time per time step' metric from BENCH TU files
across git history. One plot per test case, x-axis = time (commit dates),
one curve per hostname (only hostnames with >= MIN_POINTS data points).
Scientific paper style.
"""
import subprocess, re, os
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime

GPU_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = subprocess.check_output(
    ["git", "rev-parse", "--show-toplevel"], text=True, cwd=GPU_DIR
).strip()

MIN_POINTS = 3

# ---------------------------------------------------------------------------
# Scientific paper style
# ---------------------------------------------------------------------------
plt.rcParams.update({
    # Font: Latin Modern Roman (LaTeX-like)
    'font.family':       'serif',
    'font.serif':        ['Latin Modern Roman', 'DejaVu Serif', 'Times New Roman'],
    'mathtext.fontset':  'cm',

    # Font sizes (suitable for single-column figures)
    'font.size':         11,
    'axes.titlesize':    13,
    'axes.labelsize':    12,
    'xtick.labelsize':   10,
    'ytick.labelsize':   10,
    'legend.fontsize':   9,

    # Lines & markers
    'lines.linewidth':   1.8,
    'lines.markersize':  5,

    # Axes
    'axes.linewidth':    0.8,
    'axes.grid':         True,
    'grid.alpha':        0.25,
    'grid.linewidth':    0.5,
    'grid.linestyle':    '--',

    # Ticks
    'xtick.direction':   'in',
    'ytick.direction':   'in',
    'xtick.major.size':  5,
    'ytick.major.size':  5,
    'xtick.minor.size':  3,
    'ytick.minor.size':  3,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'xtick.top':         True,
    'ytick.right':       True,

    # Legend
    'legend.frameon':       True,
    'legend.framealpha':    0.9,
    'legend.edgecolor':     '0.8',
    'legend.fancybox':      False,
    'legend.borderpad':     0.4,
    'legend.handlelength':  2.0,

    # Figure
    'figure.dpi':        150,
    'savefig.dpi':       300,
    'savefig.bbox':      'tight',
    'savefig.pad_inches': 0.05,
})

# Distinct markers for up to ~12 curves
MARKERS = ['o', 's', '^', 'D', 'v', 'P', 'X', 'p', 'h', '*', '<', '>']

# Colorblind-friendly palette (Tol's bright)
COLORS = [
    '#4477AA',  # blue
    '#EE6677',  # red
    '#228833',  # green
    '#CCBB44',  # yellow
    '#66CCEE',  # cyan
    '#AA3377',  # purple
    '#BBBBBB',  # grey
    '#EE8866',  # orange
    '#44BB99',  # teal
    '#332288',  # indigo
    '#882255',  # wine
    '#117733',  # dark green
]

# ---------------------------------------------------------------------------
# Explicit rename mapping: new_path -> old_path (from commit d40abe2f)
# ---------------------------------------------------------------------------
RENAME_MAP = {
    "tests/GPU/ColdLegMixing/ColdLegMixing_BENCH.TU.is157091_cc86":    "tests/GPU/ColdLegMixing/ColdLegMixing_BENCH.TU.ref_a6000",
    "tests/GPU/ColdLegMixing/ColdLegMixing_BENCH.TU.topaze_cc80":      "tests/GPU/ColdLegMixing/ColdLegMixing_BENCH.TU.ref_a100",
    "tests/GPU/DomainFlowLES/DomainFlowLES_BENCH.TU.is157091_cc86":    "tests/GPU/DomainFlowLES/DomainFlowLES_BENCH.TU.ref_a6000",
    "tests/GPU/DomainFlowLES/DomainFlowLES_BENCH.TU.topaze_cc80":      "tests/GPU/DomainFlowLES/DomainFlowLES_BENCH.TU.ref_a100",
    "tests/GPU/GPU4/GPU4_BENCH.TU.topaze_cc80":                        "tests/GPU/GPU4/GPU4_BENCH.TU.ref_a100",
    "tests/GPU/JEL_bous/JEL_bous_BENCH.TU.is157091_cc86":              "tests/GPU/JEL_bous/JEL_bous_BENCH.TU.ref_a6000",
    "tests/GPU/JEL_bous/JEL_bous_BENCH.TU.topaze_cc80":                "tests/GPU/JEL_bous/JEL_bous_BENCH.TU.ref_a100",
    "tests/GPU/OpenMP_Iterateur/OpenMP_Iterateur_BENCH_AmgX.TU.is157091_cc86":  "tests/GPU/OpenMP_Iterateur/OpenMP_Iterateur_BENCH_AmgX.TU.ref_a6000",
    "tests/GPU/OpenMP_Iterateur/OpenMP_Iterateur_BENCH_AmgX.TU.is246827_cc86":  "tests/GPU/OpenMP_Iterateur/OpenMP_Iterateur_BENCH_AmgX.TU.ref_a3000",
    "tests/GPU/OpenMP_Iterateur/OpenMP_Iterateur_BENCH_AmgX.TU.jean-zay_cc70":  "tests/GPU/OpenMP_Iterateur/OpenMP_Iterateur_BENCH_AmgX.TU.ref_v100-32g",
    "tests/GPU/OpenMP_Iterateur/OpenMP_Iterateur_BENCH_AmgX.TU.petra_cc86":     "tests/GPU/OpenMP_Iterateur/OpenMP_Iterateur_BENCH_AmgX.TU.ref_a5000",
    "tests/GPU/OpenMP_Iterateur/OpenMP_Iterateur_BENCH_AmgX.TU.topaze_cc80":    "tests/GPU/OpenMP_Iterateur/OpenMP_Iterateur_BENCH_AmgX.TU.ref_a100",
    "tests/GPU/OpenMP_Iterateur/OpenMP_Iterateur_BENCH_rocALUTION.TU.adastra_gfx90a": "tests/GPU/OpenMP_Iterateur/OpenMP_Iterateur_BENCH_rocALUTION.TU.ref_MI250X",
    "tests/GPU/OpenMP_QC/OpenMP_QC_BENCH.TU.is157091_cc86":            "tests/GPU/OpenMP_QC/OpenMP_QC_BENCH.TU.ref_a6000",
    "tests/GPU/OpenMP_QC/OpenMP_QC_BENCH.TU.topaze_cc80":              "tests/GPU/OpenMP_QC/OpenMP_QC_BENCH.TU.ref_a100",
    "tests/GPU/thermohydraulique_VEF_DNS/thermohydraulique_VEF_DNS_BENCH.TU.is157091_cc86": "tests/GPU/thermohydraulique_VEF_DNS/thermohydraulique_VEF_DNS_BENCH.TU.ref_a100",
    "tests/GPU/thermohydraulique_VEF_DNS/thermohydraulique_VEF_DNS_BENCH.TU.topaze_cc80":   "tests/GPU/thermohydraulique_VEF_DNS/thermohydraulique_VEF_DNS_BENCH.TU.ref_a6000",
}

# ---------------------------------------------------------------------------
# Hostname -> human-readable GPU label for legends
# ---------------------------------------------------------------------------
HOST_LABELS = {
    'adastra_gfx90a':      'Adastra MI250X',
    'adastra_gfx942':      'Adastra MI300A',
    'b200_cc100':           'B200',
    'b300_cc100':           'B300',
    'calypso-grace_cc90':   'Grace Hopper',
    'h200_cc90':            'H200',
    'irene-amd-ccrt_cc70':  'Irene V100',
    'is157091_cc86':        'A6000',
    'is246827_cc86':        'A3000',
    'is247793_gfx1100':     'RX 7600',
    'jean-zay_cc70':        'Jean-Zay V100',
    'jean-zay_cc90':        'Jean-Zay H100',
    'lumi_gfx90a':          'LUMI MI250X',
    'marenostrum_cc90':     'MareNostrum H100',
    'orcus-amd_cc70':       'Orcus V100',
    'orcus-amd_cc90':       'Orcus H100',
    'petra_cc86':           'Petra A5000',
    'topaze_cc80':          'Topaze A100',
}

def host_label(hostname):
    return HOST_LABELS.get(hostname, hostname)

# ---------------------------------------------------------------------------
# Git helpers
# ---------------------------------------------------------------------------
def git(*args):
    return subprocess.check_output(
        ["git", "-C", REPO_ROOT] + list(args),
        text=True, stderr=subprocess.DEVNULL
    )

def get_tracked_bench_files():
    out = git("ls-files", "--", "tests/GPU/*_BENCH.TU.*",
              "tests/GPU/OpenMP_Iterateur/OpenMP_Iterateur.TU.*").strip().split('\n')
    tests = defaultdict(list)
    for f in out:
        if '/PAR_' in f or '.swp' in f or f.endswith('.old'):
            continue
        rel = f.replace('tests/GPU/', '', 1)
        d = os.path.dirname(rel)
        host = f.split('.TU.')[-1]
        tests[d].append((f, host))
    return tests

def extract_avg_time_step(content):
    m = re.search(r'Average time per time step:\s+([\d.]+)', content)
    if m:
        return float(m.group(1))
    m = re.search(r'Secondes / pas de temps\s+([\d.]+)', content)
    if m:
        return float(m.group(1))
    return None

def extract_dof(content):
    m = re.search(r'Total number of elements used for the calculation:\s+(\d+)', content)
    if m:
        return int(m.group(1))
    return None

def get_history_for_path(repo_path):
    try:
        raw = git("log", "--format=%H|%ai|%s", "--", repo_path)
    except subprocess.CalledProcessError:
        return []
    entries = []
    for line in raw.strip().split('\n'):
        if not line:
            continue
        parts = line.split('|', 2)
        if len(parts) == 3:
            entries.append((parts[0], parts[1][:10], parts[2], repo_path))
    return entries

def get_full_history(repo_path):
    entries = get_history_for_path(repo_path)
    old_path = RENAME_MAP.get(repo_path)
    if old_path:
        entries = entries + get_history_for_path(old_path)
    seen = set()
    unique = []
    for e in reversed(entries):
        if e[0] not in seen:
            seen.add(e[0])
            unique.append(e)
    return unique

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    tests = get_tracked_bench_files()
    print(f"Found {len(tests)} test cases:")
    for t in sorted(tests):
        print(f"  {t}: {len(tests[t])} hostnames")

    os.makedirs(os.path.join(GPU_DIR, "bench_plots"), exist_ok=True)

    # Collect all valid test cases data first
    all_plot_data = []  # list of (test_case, host_series)

    for test_case in sorted(tests):
        files = tests[test_case]
        print(f"\nProcessing {test_case}...")

        # DOF is shared across all TU files in the same directory: find it once
        shared_dof = None
        for repo_path, _ in files:
            disk_path = os.path.join(REPO_ROOT, repo_path)
            if os.path.isfile(disk_path):
                with open(disk_path) as f:
                    shared_dof = extract_dof(f.read())
                if shared_dof is not None:
                    break
        if shared_dof is None:
            print(f"  No DOF found in any file for {test_case}, skipping")
            continue

        host_series = {}
        for repo_path, hostname in files:
            entries = get_full_history(repo_path)
            if not entries:
                continue

            dof = shared_dof
            series = []
            prev_val = None
            for chash, date, subj, fpath in entries:
                try:
                    content = git("show", f"{chash}:{fpath}")
                except subprocess.CalledProcessError:
                    continue
                val = extract_avg_time_step(content)
                if val is not None and val != prev_val and val > 0:
                    mdof_s = dof / (val * 1e6)
                    series.append((datetime.strptime(date, "%Y-%m-%d"), mdof_s))
                    prev_val = val

            if len(series) >= MIN_POINTS:
                host_series[hostname] = series
                print(f"  {host_label(hostname):20s} ({hostname}): {len(series)} points "
                      f"[{dof:,} DOF -> {series[-1][1]:.2f} MDOF/s] "
                      f"({series[0][0].strftime('%Y-%m-%d')} -> "
                      f"{series[-1][0].strftime('%Y-%m-%d')})")

        if not host_series:
            print(f"  No hostname with >= {MIN_POINTS} points, skipping")
            continue

        all_plot_data.append((test_case, host_series))

    if not all_plot_data:
        print("No data to plot.")
        return

    # ---- Build a global hostname -> (color, marker) mapping --------------
    all_hosts = sorted({h for _, hs in all_plot_data for h in hs})
    host_style = {
        host: (COLORS[i % len(COLORS)], MARKERS[i % len(MARKERS)])
        for i, host in enumerate(all_hosts)
    }

    # ---- Individual plots ------------------------------------------------
    for test_case, host_series in all_plot_data:
        fig, ax = plt.subplots(figsize=(7, 4.5))
        _plot_single(ax, test_case, host_series, host_style)
        outfile = os.path.join(GPU_DIR, "bench_plots",
                               f"{test_case.replace('/', '_')}.png")
        fig.savefig(outfile)
        plt.close(fig)
        print(f"  Saved {outfile}")

    # ---- Combined figure -------------------------------------------------
    n = len(all_plot_data)
    ncols = min(3, n)
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(7 * ncols, 4.5 * nrows),
                             squeeze=False)

    for idx, (test_case, host_series) in enumerate(all_plot_data):
        row, col = divmod(idx, ncols)
        ax = axes[row][col]
        _plot_single(ax, test_case, host_series, host_style)

    # Hide unused subplots
    for idx in range(n, nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row][col].set_visible(False)

    fig.suptitle('GPU Benchmark — Performance Evolution (MDOF/s)',
                 fontsize=15, fontweight='bold', y=1.01)
    plt.tight_layout()

    combined = os.path.join(GPU_DIR, "bench_plots", "all_benchmarks.png")
    fig.savefig(combined)
    plt.close(fig)
    print(f"\n  Saved combined figure: {combined}")
    print("\nDone! All plots saved in bench_plots/")


def _plot_single(ax, test_case, host_series, host_style):
    """Draw one test case on a given axes."""
    for host in sorted(host_series):
        series = host_series[host]
        dates  = [s[0] for s in series]
        values = [s[1] for s in series]
        color, marker = host_style[host]
        ax.plot(dates, values,
                marker=marker,
                color=color,
                label=host_label(host),
                markeredgecolor='white',
                markeredgewidth=0.4)

    ax.set_xlabel('Date')
    ax.set_ylabel('Performance (MDOF/s)')
    ax.set_title(test_case.replace('_', ' '))
    ax.set_ylim(bottom=0)

    # Date formatting
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
    ax.xaxis.set_minor_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b\n%Y'))
    for lbl in ax.get_xticklabels():
        lbl.set_rotation(30)
        lbl.set_ha('right')

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.20),
              ncol=min(4, len(host_series)),
              columnspacing=1.0, handletextpad=0.5)

if __name__ == '__main__':
    main()
