"""
Dead-zone / overlap map builder + scalar dead-zone-fraction extractor.

Requires a per-hit table (one row per hit, not the aggregated summary JSON)
with at minimum:
    structure    -- geometry name (e.g. "kresling", "barrel", "plate")
    run          -- run identifier, as a string (e.g. "3", "0_2")
    x, y, z      -- global hit position (mm), OR
    u_local, v_local, facet_id  -- pre-computed local facet coords, if your
                                    SD tags them
    path_X0      -- path length in radiation lengths for that hit

This is exactly the schema extract_hits_for_deadzone.py writes (it reads
"structure" from the PionEvents tree's structureTag branch, and "run" from
the source .root filename). If you're loading hits from somewhere else,
adjust load_hits() and make sure "structure"/"run" columns exist one way
or another.

Two entry points:
  - build_maps(...)          -> 2D ratio/overlap maps + saves a plot (diagnostic)
  - dead_zone_fraction(...)  -> single scalar per run, for feeding into the
                                same summary table aggregate_and_predict.py
                                builds for eta and B(N)

NOTE ON EXCLUDED GEOMETRIES
----------------------------
The "plate" geometry is intentionally excluded from both entry points below
(alongside "barrel", which is always excluded as the baseline itself, not a
geometry to map). This mirrors the exclusion already applied upstream in
deadzone_analysis.ps1, which only resolves the .root files it's pointed at
-- this is a second line of defense in case a stray plate hits CSV from an
older run is still sitting in hits_dir.
"""

import glob
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Structure-name PREFIXES skipped outright when scanning hits_dir: "barrel"
# is the baseline itself (never mapped against itself), "plate" is excluded
# from this analysis per the same exclusion applied in deadzone_analysis.ps1.
# Matched as a prefix (via is_structure()), not an exact match, because
# structureTag values in the tree carry suffixes like "_reference" /
# "_deployed" (e.g. "barrel_reference", "kresling_deployed") rather than
# the bare names.
EXCLUDED_STRUCTURE_PREFIXES = ("barrel", "plate")

REQUIRED_COLUMNS = {"structure", "run", "path_X0"}


def is_structure(value: str, prefix: str) -> bool:
    """True if `value` (a structureTag, e.g. 'barrel_reference') identifies
    the geometry named by `prefix` (e.g. 'barrel'). Matches the prefix
    followed by either nothing or an underscore, case-insensitively, so
    'barrel' matches 'barrel_reference' and 'barrel' but not 'barrelish'."""
    value = value.lower()
    prefix = prefix.lower()
    return value == prefix or value.startswith(prefix + "_")


def is_excluded_structure(value: str) -> bool:
    return any(is_structure(value, p) for p in EXCLUDED_STRUCTURE_PREFIXES)


def load_hits(path):
    # Swap for pd.read_hdf(...) or uproot.open(...).arrays(library="pd")
    # depending on what your pipeline actually writes. Whatever the source,
    # the result needs "structure" and "run" columns -- see module docstring.
    df = pd.read_csv(path, dtype={"run": str})
    missing = REQUIRED_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(
            f"'{path}' is missing required column(s) {sorted(missing)}. "
            "Was this file produced by the current extract_hits_for_deadzone.py? "
            "Older hits CSVs (pre structure/run columns) need re-extracting."
        )
    return df


def unroll(df):
    theta = np.degrees(np.arctan2(df["y"], df["x"])) % 360
    z = df["z"]
    return theta, z
    # If your sensitive detector tags local facet (u, v, facet_id), prefer
    # those directly -- they resolve hinge lines much more cleanly than
    # global (theta, z) once facet curvature is involved.


def build_maps(fold_hits: pd.DataFrame, barrel_hits: pd.DataFrame,
               n_theta_bins=72, n_z_bins=40):
    """Returns (ratio_map, mean_pathX0_map, theta_edges, z_edges)."""
    theta_f, z_f = unroll(fold_hits)
    theta_b, z_b = unroll(barrel_hits)

    theta_edges = np.linspace(0, 360, n_theta_bins + 1)
    z_edges = np.linspace(
        min(z_f.min(), z_b.min()), max(z_f.max(), z_b.max()), n_z_bins + 1
    )

    count_f, _, _ = np.histogram2d(theta_f, z_f, bins=[theta_edges, z_edges])
    count_b, _, _ = np.histogram2d(theta_b, z_b, bins=[theta_edges, z_edges])

    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(count_b > 0, count_f / count_b, np.nan)

    sum_pathX0, _, _ = np.histogram2d(
        theta_f, z_f, bins=[theta_edges, z_edges], weights=fold_hits["path_X0"]
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        mean_pathX0 = np.where(count_f > 0, sum_pathX0 / count_f, np.nan)

    return ratio, mean_pathX0, theta_edges, z_edges


def dead_zone_fraction(ratio_map: np.ndarray, threshold: float = 0.5) -> float:
    """Fraction of geometrically-covered bins (barrel has hits) where the
    fold geometry's hit density falls below `threshold` x the barrel's.
    This is the scalar to log per (structure, N) run."""
    covered = ~np.isnan(ratio_map)
    if covered.sum() == 0:
        return float("nan")
    dead = covered & (ratio_map < threshold)
    return dead.sum() / covered.sum()


def plot_maps(ratio, mean_pathX0, theta_edges, z_edges, out_path):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    im0 = axes[0].imshow(
        ratio.T, origin="lower", aspect="auto",
        extent=[theta_edges[0], theta_edges[-1], z_edges[0], z_edges[-1]],
        cmap="RdBu_r", vmin=0, vmax=2
    )
    axes[0].set_title("Hit-density ratio vs. barrel\n(< 1 = Relative Under-coverage, > 1 = overlap-ish)")
    axes[0].set_xlabel("theta (deg)"); axes[0].set_ylabel("z (mm)")
    fig.colorbar(im0, ax=axes[0])

    im1 = axes[1].imshow(
        mean_pathX0.T, origin="lower", aspect="auto",
        extent=[theta_edges[0], theta_edges[-1], z_edges[0], z_edges[-1]],
        cmap="viridis"
    )
    axes[1].set_title("Mean path length (X0) per bin\n(elevated = overlap seam)")
    axes[1].set_xlabel("theta (deg)"); axes[1].set_ylabel("z (mm)")
    fig.colorbar(im1, ax=axes[1])

    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close(fig)


def _load_all_hits(hits_dir: str) -> pd.DataFrame:
    """Loads every *_hits.csv in hits_dir, tags each row's source file, and
    concatenates. Structure/run come from the columns written by
    extract_hits_for_deadzone.py, not from filenames."""
    frames = []
    for path in sorted(glob.glob(f"{hits_dir}/*.csv")):
        name = Path(path).name
        if name.endswith(("_dead_zone_summary.csv", "_dead_zone_summary_aggregate.csv",
                           "dead_zone_N_sweep.csv")):
            continue  # skip this script's own output tables if re-run in the same dir
        df = load_hits(path)
        df["__source_file"] = name
        frames.append(df)
    if not frames:
        return pd.DataFrame(columns=["structure", "run", "path_X0", "__source_file"])
    return pd.concat(frames, ignore_index=True)


def sweep_dead_zone_fraction(hits_dir: str, barrel_hits_path: str,
                              threshold: float = 0.5,
                              out_csv: str = "dead_zone_N_sweep.csv",
                              save_maps: bool = False):
    """Runs dead_zone_fraction over every per-hit CSV in hits_dir (grouped by
    its "structure" and "run" columns), against a single shared
    barrel_hits_path baseline. Writes a tidy CSV keyed the same way as the
    other sweep scripts, so it can be merged into aggregate_and_predict.py's
    output on (structure, N).

    NOTE: "N" in the output table is the hits CSV's "run" column, taken
    as-is (as a string). If your sweep uses "run" to mean something else
    than the N-value (e.g. a repeated-seed index within a fixed N), rename/
    remap it before calling this, or post-process the output CSV's "N" column.

    Structures matching EXCLUDED_STRUCTURE_PREFIXES (currently "barrel" and
    "plate", matched as prefixes -- see is_structure()) are skipped, same
    as build_original_four_maps()."""
    barrel_hits = load_hits(barrel_hits_path)

    all_hits = _load_all_hits(hits_dir)
    if all_hits.empty:
        print(f"\nNo hits CSVs found in '{hits_dir}'. Nothing to summarize.")
        return pd.DataFrame(columns=["structure", "N", "dead_zone_fraction"])

    rows = []
    for (structure, run), fold_hits in all_hits.groupby(["structure", "run"]):
        if is_excluded_structure(str(structure)):
            print(f"skipping excluded structure '{structure}' (run {run})")
            continue

        ratio, mean_pathX0, theta_edges, z_edges = build_maps(fold_hits, barrel_hits)
        dz = dead_zone_fraction(ratio, threshold=threshold)

        if save_maps:
            plot_maps(ratio, mean_pathX0, theta_edges, z_edges,
                       f"{hits_dir}/map_{structure}_N{run}.png")

        rows.append({
            "structure": structure,
            "N": run,
            "dead_zone_fraction": dz,
        })

    if not rows:
        print(f"\nNo structures left in '{hits_dir}' after excluding "
              f"{sorted(EXCLUDED_STRUCTURE_PREFIXES)}. Nothing to summarize.")
        return pd.DataFrame(columns=["structure", "N", "dead_zone_fraction"])

    df = pd.DataFrame(rows).sort_values(["structure", "N"])
    df.to_csv(out_csv, index=False)
    print(df.to_string(index=False))
    return df


def build_original_four_maps(hits_dir: str, barrel_hits_dir: str = None,
                              threshold: float = 0.5):
    """
    For the ORIGINAL geometries (Yoshimura, Miura, Kresling -- single
    fold1p00 state) run as repeated MC trials -- not the N-sweep table
    sweep_dead_zone_fraction() builds. Loads every *_hits.csv in hits_dir
    and groups rows by their "structure" and "run" columns (written by
    extract_hits_for_deadzone.py from the PionEvents tree's structureTag
    branch and the source .root filename, respectively) -- no filename
    convention required.

    IMPORTANT -- per-run baseline matching:
    Each fold-geometry run is compared against the barrel run with the SAME
    run value (e.g. kresling run "5" is mapped against barrel run "5", NOT
    run "3" or any other run). This matters because each run is an
    independent MC trial with its own random seed -- the barrel file itself
    carries its own statistical fluctuations, so comparing fold-run-N
    against a DIFFERENT barrel run injects noise from two uncorrelated MC
    samples into the ratio, which has nothing to do with the fold
    geometry's actual dead-zone behavior. It would inflate the apparent
    run-to-run variance and can bias the mean if any one barrel run happens
    to be a statistical outlier. Matching run-for-run keeps that source of
    noise out of the comparison.

    barrel_hits_dir defaults to hits_dir (i.e. the barrel hits CSVs live
    alongside the fold-geometry files) -- pass a different directory only
    if your barrel baselines are stored elsewhere.

    Saves one diagnostic plot PER RUN (each run gets its own PNG so
    repeated MC trials don't overwrite each other's diagnostic map), and
    returns two tables:

      - a per-run table (one row per structure/run, written to
        original_four_dead_zone_summary.csv) -- the raw values, each row
        computed against its own matching-run barrel baseline
      - a per-structure aggregate table (mean/std/sem/count across runs,
        written to original_four_dead_zone_summary_aggregate.csv) -- the
        statistically-validated number you'd actually quote, given multiple
        MC repeats per structure

    "barrel" is always skipped as a fold-geometry candidate (it's the
    baseline, not something to map against itself), and "plate" is skipped
    as well -- see EXCLUDED_STRUCTURE_PREFIXES at the top of this file.
    Matching is by prefix (via is_structure()), so a structureTag like
    "barrel_reference" is recognized as "barrel", and "kresling_deployed"
    is recognized as the "kresling" structure (not excluded, just labeled
    with its full tag in the output).
    """
    if barrel_hits_dir is None:
        barrel_hits_dir = hits_dir

    all_hits = _load_all_hits(hits_dir)
    barrel_pool = all_hits if barrel_hits_dir == hits_dir else _load_all_hits(barrel_hits_dir)

    # Build a {run: barrel DataFrame} lookup up front, so each fold group can
    # find its matching-run baseline (see docstring above for why this must
    # be per-run rather than a single shared barrel file). Matched by prefix
    # since structureTag may be "barrel_reference" rather than bare "barrel".
    is_barrel_mask = barrel_pool["structure"].apply(lambda s: is_structure(str(s), "barrel"))
    barrel_rows = barrel_pool[is_barrel_mask]
    if barrel_rows.empty:
        seen = sorted(barrel_pool["structure"].unique()) if not barrel_pool.empty else []
        raise FileNotFoundError(
            f"No rows identifying as 'barrel' (exact or '<name>_...' prefix) found in "
            f"'{barrel_hits_dir}'. structure values actually present: {seen}. "
            "Can't build per-run baselines."
        )
    barrel_by_run = {run: df for run, df in barrel_rows.groupby("run")}
    print(f"Found barrel baselines (structureTag e.g. "
          f"'{barrel_rows['structure'].iloc[0]}') for runs: {sorted(barrel_by_run)}")

    rows = []
    fold_hits_pool = all_hits[~is_barrel_mask] if barrel_hits_dir == hits_dir else all_hits[
        ~all_hits["structure"].apply(lambda s: is_structure(str(s), "barrel"))
    ]
    for (structure, run), fold_hits in fold_hits_pool.groupby(["structure", "run"]):
        if is_excluded_structure(structure):
            print(f"skipping excluded structure '{structure}' (run {run})")
            continue

        if run not in barrel_by_run:
            print(f"skipping {structure} run {run}: no matching barrel run found in "
                  f"'{barrel_hits_dir}' (available barrel runs: "
                  f"{sorted(barrel_by_run)}) -- refusing to fall back "
                  "to a different run's barrel baseline.")
            continue

        barrel_hits = barrel_by_run[run]

        ratio, mean_pathX0, theta_edges, z_edges = build_maps(fold_hits, barrel_hits)
        dz = dead_zone_fraction(ratio, threshold=threshold)

        out_png = f"{hits_dir}/map_{structure}_run{run}_vs_barrel.png"
        plot_maps(ratio, mean_pathX0, theta_edges, z_edges, out_png)

        rows.append({"structure": structure, "run": run, "dead_zone_fraction": dz, "map": out_png})
        print(f"{structure} run {run}: dead_zone_fraction={dz:.4f}  (threshold={threshold}, "
              f"barrel run {run})  -> {out_png}")

    if not rows:
        print(f"\nNo (structure, run) groups in '{hits_dir}' survived exclusion/matching "
              f"(excluded structure prefixes: {sorted(EXCLUDED_STRUCTURE_PREFIXES)}). "
              "Nothing to summarize -- check that hits CSVs have 'structure'/'run' columns "
              "and matching barrel runs exist.")
        empty = pd.DataFrame(columns=["structure", "run", "dead_zone_fraction", "map"])
        return empty, empty

    df = pd.DataFrame(rows).sort_values(["structure", "run"])
    df.to_csv(f"{hits_dir}/original_four_dead_zone_summary.csv", index=False)

    # Aggregate across the repeated MC runs per structure -- this is the
    # number to actually quote/plot, with an uncertainty that reflects
    # run-to-run MC variance rather than a single (possibly unlucky) run.
    agg = (
        df.groupby("structure")["dead_zone_fraction"]
        .agg(mean="mean", std="std", n="count")
        .reset_index()
    )
    # std is NaN for n=1 (no variance to compute); sem is NaN there too.
    agg["sem"] = agg["std"] / np.sqrt(agg["n"])
    agg = agg.sort_values("structure")
    agg.to_csv(f"{hits_dir}/original_four_dead_zone_summary_aggregate.csv", index=False)

    print("\nPer-structure aggregate over MC runs:")
    print(agg.to_string(index=False))

    return df, agg


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Dead-zone / overlap mapping from per-hit CSVs."
    )
    parser.add_argument(
        "mode", choices=["original", "sweep"],
        help="'original' = the single-state geometries (yoshimura/miura/"
             "kresling vs. barrel), grouped by the 'structure'/'run' columns "
             "extract_hits_for_deadzone.py writes; 'plate' is excluded. "
             "'sweep' = the N-sweep table (many structure/run rows treated "
             "as structure/N; 'plate' is excluded here too)."
    )
    parser.add_argument(
        "hits_dir",
        help="Directory CONTAINING the *_hits.csv files (not a path to one CSV itself)."
    )
    parser.add_argument(
        "barrel_hits_path_or_dir",
        help="For mode='original': directory containing hits CSVs with "
             "structure=='barrel' rows, one run each, matched run-for-run "
             "against each fold geometry's run (defaults to hits_dir itself "
             "if you pass the same directory here). "
             "For mode='sweep': path to the single barrel reference "
             "*_hits.csv file (the sweep table isn't per-run repeated MC "
             "trials, so a single shared baseline applies there)."
    )
    parser.add_argument("--threshold", type=float, default=0.5)
    args = parser.parse_args()

    if args.mode == "original":
        build_original_four_maps(args.hits_dir, args.barrel_hits_path_or_dir, threshold=args.threshold)
    else:
        sweep_dead_zone_fraction(args.hits_dir, args.barrel_hits_path_or_dir,
                                  threshold=args.threshold, save_maps=True)