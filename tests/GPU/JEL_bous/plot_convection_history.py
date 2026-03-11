#!/usr/bin/env python3
"""
Plot evolution of 'Convection operator' time/step metric from
JEL_bous_BENCH.TU.is157091_cc86 (A6000) across git history.

New format (>=99e7eac):
  Convection operator                      | 0.04150604      |  7.3        | 4

Old format (<99e7eac):
  Dont operateurs convection        0.244672 18% (4 appels/pas de temps)
  -> 4th whitespace-delimited token
"""
import subprocess, re, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime

GPU_DIR   = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = subprocess.check_output(
    ["git", "rev-parse", "--show-toplevel"], text=True, cwd=GPU_DIR
).strip()

# Current and historical path for A6000
NEW_PATH = "tests/GPU/JEL_bous/JEL_bous_BENCH.TU.is157091_cc86"
OLD_PATH = "tests/GPU/JEL_bous/JEL_bous_BENCH.TU.ref_a6000"

# ---------------------------------------------------------------------------
# Scientific paper style (mirrored from plot_bench_history.py)
# ---------------------------------------------------------------------------
plt.rcParams.update({
    'font.family':       'serif',
    'font.serif':        ['Latin Modern Roman', 'DejaVu Serif', 'Times New Roman'],
    'mathtext.fontset':  'cm',
    'font.size':         11,
    'axes.titlesize':    13,
    'axes.labelsize':    12,
    'xtick.labelsize':   10,
    'ytick.labelsize':   10,
    'legend.fontsize':   9,
    'lines.linewidth':   1.8,
    'lines.markersize':  5,
    'axes.linewidth':    0.8,
    'axes.grid':         True,
    'grid.alpha':        0.25,
    'grid.linewidth':    0.5,
    'grid.linestyle':    '--',
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
    'legend.frameon':       True,
    'legend.framealpha':    0.9,
    'legend.edgecolor':     '0.8',
    'legend.fancybox':      False,
    'figure.dpi':        150,
    'savefig.dpi':       300,
    'savefig.bbox':      'tight',
    'savefig.pad_inches': 0.05,
})

COLOR  = '#4477AA'   # Tol bright blue
MARKER = 'o'

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def git(*args):
    return subprocess.check_output(
        ["git", "-C", REPO_ROOT] + list(args),
        text=True, stderr=subprocess.DEVNULL
    )

def extract_convection(content):
    """Return convection time/step from either new or old TU format."""
    # New format: pipe-separated table
    m = re.search(r'Convection operator\s*\|\s*([\d.eE+\-]+)', content)
    if m:
        return float(m.group(1))
    # Old format: "Dont operateurs convection   VALUE   %  ..."
    m = re.search(r'Dont operateurs convection\s+([\d.eE+\-]+)', content)
    if m:
        return float(m.group(1))
    return None

def get_history(path):
    try:
        raw = git("log", "--format=%H|%ai|%s", "--", path)
    except subprocess.CalledProcessError:
        return []
    entries = []
    for line in raw.strip().split('\n'):
        if not line:
            continue
        parts = line.split('|', 2)
        if len(parts) == 3:
            entries.append((parts[0], parts[1][:10], parts[2], path))
    return entries

# ---------------------------------------------------------------------------
# Collect data across git history (new path + old path merged)
# ---------------------------------------------------------------------------
entries = get_history(NEW_PATH) + get_history(OLD_PATH)

# Deduplicate by commit hash, keep chronological order
seen = set()
unique = []
for e in reversed(entries):
    if e[0] not in seen:
        seen.add(e[0])
        unique.append(e)

print(f"Total unique commits: {len(unique)}")

series = []    # list of (datetime, value_s)
prev_val = None
for chash, date_str, subj, fpath in unique:
    try:
        content = git("show", f"{chash}:{fpath}")
    except subprocess.CalledProcessError:
        continue
    val = extract_convection(content)
    if val is not None and val > 0 and val != prev_val:
        dt = datetime.strptime(date_str, "%Y-%m-%d")
        series.append((dt, val, subj))
        print(f"  {date_str}  {val:.6f} s/step   {subj[:70]}")
        prev_val = val

if not series:
    print("No data found.")
    exit(1)

dates  = [s[0] for s in series]
values = [s[1] for s in series]

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 4.5))

ax.plot(dates, values,
        marker=MARKER,
        color=COLOR,
        label='A6000 (is157091)',
        markeredgecolor='white',
        markeredgewidth=0.5,
        zorder=3)

# Annotate notable changes (commits where value differs >5% from previous)
prev = values[0]
for dt, val, subj in series[1:]:
    pct = (val - prev) / prev * 100
    if abs(pct) > 5:
        ax.annotate(
            f'{pct:+.0f}%',
            xy=(dt, val),
            xytext=(0, 10 if pct > 0 else -15),
            textcoords='offset points',
            fontsize=7,
            color='#228833' if pct < 0 else '#EE6677',
            ha='center',
        )
    prev = val

ax.set_xlabel('Date')
ax.set_ylabel('Convection operator — time per step (s)')
ax.set_title('JEL_bous \u2014 Convection operator performance history\n'
             '(A6000, NVIDIA RTX A6000, cc86)')
ax.set_ylim(bottom=0)

# Date axis formatting
ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
ax.xaxis.set_minor_locator(mdates.MonthLocator())
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b\n%Y'))
for lbl in ax.get_xticklabels():
    lbl.set_rotation(30)
    lbl.set_ha('right')

ax.legend(loc='upper right')

plt.tight_layout()
outfile = os.path.join(GPU_DIR, "bench_plots", "JEL_bous_convection_history.png")
os.makedirs(os.path.dirname(outfile), exist_ok=True)
fig.savefig(outfile)
plt.close(fig)
print(f"\nSaved: {outfile}")
