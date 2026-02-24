#!/usr/bin/env python3
"""
Plot a stacked bar chart of 'Average time per time step' breakdown
from *_BENCH.TU.* files in the current directory.

X-axis  : GPU/CPU architectures (sorted fastest → slowest)
Y-axis  : Average time per time step (s)
Segments: individual operation timings listed between
          'Linear solver resolutions Ax=B' and 'Other operations'

Handles both the new English table format and the old French format.
Scientific paper style consistent with ../plot_bench_history.py.
"""

import re
import os
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# Scientific paper style (mirrors ../plot_bench_history.py)
# ---------------------------------------------------------------------------
plt.rcParams.update({
    'font.family':        'serif',
    'font.serif':         ['Latin Modern Roman', 'DejaVu Serif', 'Times New Roman'],
    'mathtext.fontset':   'cm',
    'font.size':          11,
    'axes.titlesize':     13,
    'axes.labelsize':     12,
    'xtick.labelsize':    10,
    'ytick.labelsize':    10,
    'legend.fontsize':    9,
    'lines.linewidth':    1.8,
    'axes.linewidth':     0.8,
    'axes.grid':          False,   # managed per-axis below
    'xtick.direction':    'in',
    'ytick.direction':    'in',
    'xtick.major.size':   5,
    'ytick.major.size':   5,
    'xtick.minor.size':   3,
    'ytick.minor.size':   3,
    'xtick.major.width':  0.8,
    'ytick.major.width':  0.8,
    'xtick.top':          True,
    'ytick.right':        True,
    'legend.frameon':     True,
    'legend.framealpha':  0.9,
    'legend.edgecolor':   '0.8',
    'legend.fancybox':    False,
    'legend.borderpad':   0.4,
    'figure.dpi':         150,
    'savefig.dpi':        300,
    'savefig.bbox':       'tight',
    'savefig.pad_inches': 0.05,
})

# ---------------------------------------------------------------------------
# Hostname → human-readable GPU label  (same mapping as plot_bench_history.py)
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

# Base machine name → short label (for parallel CPU runs, hostname = <base>x<N>)
BASE_LABELS = {
    'is157091':   'A6000 host',
    'is246827':   'A3000 host',
    'topaze':     'Topaze',
    'jean-zay':   'Jean-Zay',
    'adastra':    'Adastra',
    'lumi':       'LUMI',
}


def host_label(hostname):
    """Return a human-readable label for a hostname."""
    # Detect parallel CPU suffix: hostname = <base>x<N>
    m = re.match(r'^(.+?)x(\d+)$', hostname)
    if m:
        base, n = m.group(1), m.group(2)
        base_lbl = BASE_LABELS.get(base, HOST_LABELS.get(base, base))
        return f"{base_lbl} \u00d7{n} MPI"
    return HOST_LABELS.get(hostname, hostname)


# ---------------------------------------------------------------------------
# Canonical operation names – in table order from new-format TU files.
# These are the segments of the stacked bar.
# ---------------------------------------------------------------------------
OPERATIONS = [
    "Linear solver resolutions Ax=B",
    "Matrix assembly for implicit scheme",
    "Convection operator",
    "Diffusion operator",
    "Gradient operator",
    "Divergence operator",
    "Source terms",
    "Update ::mettre_a_jour",
    "Solver for implicit diffusion",
    "Computation of the time step dt",
    "Turbulence model::update",
    "Post-treatment operations",
    "Other operations",
]

# Mapping from old French keyword (after "Dont ") to canonical name
OLD_TO_NEW = {
    'solveurs Ax=B':               'Linear solver resolutions Ax=B',
    'assemblage matrice_implicite': 'Matrix assembly for implicit scheme',
    'operateurs convection':        'Convection operator',
    'operateurs diffusion':         'Diffusion operator',
    'operateurs gradient':          'Gradient operator',
    'operateurs divergence':        'Divergence operator',
    'operateurs source':            'Source terms',
    'mettre_a_jour':                'Update ::mettre_a_jour',
    'solveur diffusion_implicite':  'Solver for implicit diffusion',
    'calcul dt':                    'Computation of the time step dt',
    'modele turbulence':            'Turbulence model::update',
    'operations postraitement':     'Post-treatment operations',
    'calcul divers':                'Other operations',
}

# ---------------------------------------------------------------------------
# Colors for each operation segment (colorblind-friendly)
# ---------------------------------------------------------------------------
OP_COLORS = {
    "Linear solver resolutions Ax=B":    '#4477AA',
    "Matrix assembly for implicit scheme": '#332288',
    "Convection operator":               '#EE6677',
    "Diffusion operator":                '#CC3311',
    "Gradient operator":                 '#228833',
    "Divergence operator":               '#44BB99',
    "Source terms":                      '#99DDFF',
    "Update ::mettre_a_jour":            '#CCBB44',
    "Solver for implicit diffusion":     '#AA3377',
    "Computation of the time step dt":   '#BBBBBB',
    "Turbulence model::update":          '#EE8866',
    "Post-treatment operations":         '#66CCEE',
    "Other operations":                  '#882255',
}

# ---------------------------------------------------------------------------
# TU file parsers
# ---------------------------------------------------------------------------

def parse_new_format(content):
    """Parse new-format TU files (English, table with '|' separators)."""
    avg_time = None
    m = re.search(r'Average time per time step:\s+([\d.eE+-]+)', content)
    if m:
        avg_time = float(m.group(1))

    operations = {}
    header_found = False
    in_table = False

    for line in content.splitlines():
        if 'Standard counter description' in line and '|' in line:
            header_found = True
            continue
        if header_found and re.match(r'\s*-{40,}', line):
            in_table = True
            header_found = False
            continue
        if in_table:
            if '|' not in line:
                continue          # skip blank or non-table lines
            parts = line.split('|')
            if len(parts) >= 2:
                name = parts[0].strip()
                try:
                    val = float(parts[1].strip().split()[0])
                    if name:
                        operations[name] = val
                        if name == 'Other operations':
                            break   # stop after last entry
                except (ValueError, IndexError):
                    pass

    return avg_time, operations


def parse_old_format(content):
    """Parse old French-format TU files."""
    avg_time = None
    m = re.search(r'Secondes / pas de temps\s+([\d.eE+-]+)', content)
    if m:
        avg_time = float(m.group(1))

    operations = {}
    for old_key, canonical in OLD_TO_NEW.items():
        m = re.search(rf'Dont\s+{re.escape(old_key)}\s+([\d.eE+-]+)', content)
        if m:
            operations[canonical] = float(m.group(1))

    return avg_time, operations


def parse_tu_file(filepath):
    """Detect format and return (avg_time_s, ops_dict)."""
    try:
        with open(filepath, encoding='utf-8', errors='replace') as fh:
            content = fh.read()
    except OSError as e:
        print(f"  Warning: cannot read {filepath}: {e}")
        return None, {}

    if 'Average time per time step:' in content:
        return parse_new_format(content)
    if 'Secondes / pas de temps' in content:
        return parse_old_format(content)
    return None, {}


# ---------------------------------------------------------------------------
# Architecture key helpers
# ---------------------------------------------------------------------------

def make_arch_key(basename):
    """
    Extract (arch_key, hostname) from a BENCH TU filename.

    Examples
    --------
    DomainFlowLES_BENCH.TU.topaze_cc80       -> ('topaze_cc80',     'topaze_cc80')
    PAR_DomainFlowLES_BENCH.TU.is157091x32   -> ('PAR/is157091x32', 'is157091x32')
    DEC_DomainFlowLES_BENCH.TU.topaze_cc80   -> ('DEC/topaze_cc80', 'topaze_cc80')
    """
    if '.TU.' not in basename:
        return None, None
    hostname = basename.split('.TU.', 1)[1]
    file_prefix = basename.split('_BENCH.TU.')[0]  # e.g. 'PAR_DomainFlowLES'

    for prefix in ('PAR_', 'DEC_', 'SEQ_'):
        if file_prefix.startswith(prefix):
            tag = prefix.rstrip('_')
            return f"{tag}/{hostname}", hostname

    return hostname, hostname


def arch_label(key):
    """Human-readable label for an arch_key."""
    if '/' in key:
        tag, host = key.split('/', 1)
        lbl = host_label(host)
        if tag == 'PAR':
            return f"{lbl}\n(CPU parallel)"
        return f"{lbl}\n({tag})"
    return host_label(key)


# ---------------------------------------------------------------------------
# Per-directory plot
# ---------------------------------------------------------------------------

def plot_directory(directory):
    """Collect all *_BENCH.TU.<hostname> files in *directory* and produce
    one stacked-bar PNG saved as bench_breakdown.png inside that directory."""

    test_name = os.path.basename(directory)
    print(f"\n=== {test_name} ===")

    candidates = glob.glob(os.path.join(directory, '*_BENCH.TU.*'))
    candidates = [
        f for f in candidates
        if not f.endswith('.swp') and not f.endswith('.old')
    ]

    data = {}   # arch_key -> (avg_time, ops_dict)

    for filepath in sorted(candidates):
        basename = os.path.basename(filepath)
        arch_key, hostname = make_arch_key(basename)
        if arch_key is None:
            continue

        avg_time, ops = parse_tu_file(filepath)
        if avg_time is None:
            print(f"  Skip {basename}: no timing data found")
            continue
        if not ops:
            print(f"  Skip {basename}: no per-operation breakdown found")
            continue

        data[arch_key] = (avg_time, ops)
        print(f"  {basename}: avg={avg_time:.4f} s, {len(ops)} operations")

    if not data:
        print(f"  No data to plot in {directory}.")
        return

    # Sort architectures by total avg time (fastest → slowest)
    sorted_archs = sorted(data.keys(), key=lambda k: data[k][0])

    # Keep only operations present in at least one file (preserve table order)
    present_ops = [op for op in OPERATIONS
                   if any(op in data[arch][1] for arch in sorted_archs)]

    # --- Build figure -------------------------------------------------------
    n_archs = len(sorted_archs)
    fig_w = max(9, n_archs * 1.4 + 5)   # extra room for legend on the right
    fig, ax = plt.subplots(figsize=(fig_w, 6))

    x = np.arange(n_archs)
    bar_width = 0.6
    bottoms = np.zeros(n_archs)

    for op_name in present_ops:
        values = np.array([data[arch][1].get(op_name, 0.0)
                           for arch in sorted_archs])
        color = OP_COLORS.get(op_name, '#888888')
        ax.bar(x, values, bar_width,
               bottom=bottoms,
               label=op_name,
               color=color,
               edgecolor='white',
               linewidth=0.4)
        bottoms += values

    # Mark the total avg_time with a horizontal tick so any gap is visible
    totals = np.array([data[arch][0] for arch in sorted_archs])
    ax.scatter(x, totals,
               marker='_', s=300, color='black', zorder=5,
               linewidths=2.0, label='Total avg. time/step')

    # --- Axes decoration ----------------------------------------------------
    x_labels = [arch_label(arch) for arch in sorted_archs]
    ax.set_xticks(x)
    ax.set_xticklabels(x_labels, ha='center')

    ax.set_ylabel('Avg. time / time step (s)')
    ax.set_title(f'{test_name} — Time per Step Breakdown by Architecture')
    ax.set_ylim(bottom=0)

    # Horizontal grid only
    ax.yaxis.grid(True, alpha=0.3, linewidth=0.5, linestyle='--')
    ax.set_axisbelow(True)

    # Legend outside the plot on the right
    ax.legend(loc='upper left',
              bbox_to_anchor=(1.01, 1.0),
              borderaxespad=0.,
              ncol=1,
              handlelength=1.5)

    plt.tight_layout()

    outfile = os.path.join(directory, 'bench_breakdown.png')
    fig.savefig(outfile)
    plt.close(fig)
    print(f"  Saved → {outfile}")


# ---------------------------------------------------------------------------
# Main – iterate over every subdirectory that contains BENCH.TU files
# ---------------------------------------------------------------------------

def main():
    # Collect unique subdirectories that own at least one *_BENCH.TU.* file
    all_files = glob.glob(os.path.join(SCRIPT_DIR, '**', '*_BENCH.TU.*'),
                          recursive=True)
    dirs = sorted({os.path.dirname(f) for f in all_files
                   if not f.endswith('.swp') and not f.endswith('.old')})

    if not dirs:
        print("No *_BENCH.TU.* files found anywhere under", SCRIPT_DIR)
        return

    for directory in dirs:
        plot_directory(directory)


if __name__ == '__main__':
    main()
