#!/usr/bin/env python3
"""
Join trajectories and merge i-pi per-iteration outputs into a single JSON.

- Trajectories expected: simulation{iter}.pos_{bead:02d}.xyz
- Outputs expected:     simulation{iter}.out
- Uses ipi.read_trajectory and ipi.read_output (your i-pi fork)

Behavior:
 - For each bead (0..nbeads-1) concatenates all existing simulation{iter}.pos_{bead:02d}.xyz
   in increasing iteration order into complete_simulation.pos_{bead:02d}.xyz and removes sources.
 - Reads simulation{i}.out for i in 0..Niterations-1. Expects ipi.read_output(...) to return
   either a dict or a tuple whose first element is the dict. Each dict maps keys -> 1D arrays.
 - The merged output_data dict is built by appending arrays. For 'time' and 'step'/'steps'
   keys the appended values are shifted by the last value currently present so sequences continue.
 - Writes root/complete_simulation.json and removes original simulation{i}.out files on success.

Usage:
  python join_and_merge.py --root N20/0/ --niterations 10 --nbeads 16
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np
import sys

# import ipi from your local fork
import ipi  # assumes ipi.read_trajectory and ipi.read_output are available


def concat_trajectories(root: Path, niterations: int, nbeads: int):
    """Concatenate per-iteration per-bead trajectory files into per-bead complete files."""
    for bead in range(nbeads):
        out_name = f"complete_simulation.pos_{bead:02d}.xyz"
        out_path = root / out_name
        any_appended = False
        with open(out_path, "wb") as fout:
            for it in range(niterations):
                src = root / f"simulation{it}.pos_{bead:02d}.xyz"
                if not src.exists():
                    # skip missing per-iteration bead file
                    continue
                any_appended = True
                sys.stdout.write(f"Appending {src.name} -> {out_name}\n")
                with open(src, "rb") as fin:
                    while True:
                        chunk = fin.read(1 << 20)
                        if not chunk:
                            break
                        fout.write(chunk)
        if any_appended:
            # remove sources for this bead
            for it in range(nbeads):  # nothing; handled below iteratively
                pass
            # delete source files for this bead across iterations
            for it in range(niterations):
                src = root / f"simulation{it}.pos_{bead:02d}.xyz"
                try:
                    if src.exists():
                        src.unlink()
                except Exception as e:
                    print(f"Warning: could not remove {src}: {e}", file=sys.stderr)
        else:
            # if no per-iteration files were found, remove the empty output file
            try:
                out_path.unlink()
            except Exception:
                pass


def normalize_read_output(ret):
    """Normalize ipi.read_output return to a dict {key: np.ndarray}.

    ipi.read_output may return a dict or (dict, desc) or similar. We take first
    element if tuple, then ensure values are numpy 1D float arrays.
    """
    if isinstance(ret, tuple) and len(ret) >= 1:
        candidate = ret[0]
    else:
        candidate = ret

    if not isinstance(candidate, dict):
        raise TypeError(f"Unexpected ipi.read_output return type: {type(ret)}")

    out = {}
    for k, v in candidate.items():
        arr = np.asarray(v)
        if arr.ndim == 0:
            arr = arr.reshape(1)
        out[str(k)] = arr.astype(float)
    return out


def merge_outputs_to_json(root: Path, niterations: int, out_json_name: str = "complete_simulation.json"):
    """
    Merge simulation{i}.out files into a single dict, shifting 'time' and 'step(s)' columns.
    Writes JSON file with lists for each key and deletes original .out files on success.
    """
    merged = {}  # key -> numpy 1D array
    time_keys = ("time", "t")
    step_keys = ("step", "steps")

    # iterate iterations in order
    for i in range(niterations):
        fname = root / f"simulation{i}.out"
        if not fname.exists():
            print(f"Warning: {fname} not found; skipping iteration {i}")
            continue
        print(f"Reading {fname.name}")
        try:
            ret = ipi.read_output(str(fname))
        except Exception as e:
            print(f"Error: ipi.read_output failed for {fname}: {e}", file=sys.stderr)
            continue

        try:
            data_dict = normalize_read_output(ret)
        except Exception as e:
            print(f"Error: cannot normalize output from {fname}: {e}", file=sys.stderr)
            continue

        if not merged:
            # first valid file: copy arrays
            for k, v in data_dict.items():
                merged[k] = v.copy()
        else:
            # append per key; if key missing in merged, initialize empty array
            for k, v in data_dict.items():
                key = k
                if key not in merged:
                    merged[key] = np.array([], dtype=float)
                if key.lower() in time_keys:
                    # shift by last time in merged
                    last = merged[key][-1] if merged[key].size > 0 else 0.0
                    merged[key] = np.append(merged[key], v + last)
                elif key.lower() in step_keys:
                    last = merged[key][-1] if merged[key].size > 0 else 0.0
                    merged[key] = np.append(merged[key], v + last)
                else:
                    merged[key] = np.append(merged[key], v)

            # If merged has keys that are not present in current data_dict, pad missing columns with NaN
            for existing_key in list(merged.keys()):
                if existing_key not in data_dict:
                    # determine how many rows current iteration had (take any present column)
                    if len(data_dict) > 0:
                        rep_len = next(iter(data_dict.values())).size
                    else:
                        rep_len = 0
                    if rep_len > 0:
                        merged[existing_key] = np.append(merged[existing_key], np.full(rep_len, np.nan))

    if not merged:
        print("No output data collected; nothing written.")
        return

    # Convert to lists and dump JSON
    json_data = {k: v.tolist() for k, v in merged.items()}
    out_path = root / out_json_name
    with open(out_path, "w") as fh:
        json.dump(json_data, fh, indent=2)
    print(f"Wrote merged output JSON to {out_path}")

    # delete original .out files
    for i in range(niterations):
        fname = root / f"simulation{i}.out"
        try:
            if fname.exists():
                fname.unlink()
        except Exception as e:
            print(f"Warning: could not remove {fname}: {e}", file=sys.stderr)


def parse_args():
    p = argparse.ArgumentParser(description="Join per-iteration bead trajs and merge per-iteration outputs.")
    p.add_argument("--root", "-r", default=".", help="root directory containing simulation*.pos_*.xyz and simulation*.out")
    p.add_argument("--niterations", "-n", type=int, required=True, help="number of iterations (Niterations)")
    p.add_argument("--nbeads", "-b", type=int, required=True, help="number of beads (nbeads)")
    p.add_argument("--out-json", default="complete_simulation.json", help="output merged JSON filename")
    return p.parse_args()


def main():
    args = parse_args()
    root = Path(args.root)
    niterations = args.niterations
    nbeads = args.nbeads

    if not root.exists():
        print(f"Error: root directory {root} does not exist", file=sys.stderr)
        sys.exit(1)

    # 1) join trajectory files per bead
    concat_trajectories(root, niterations, nbeads)

    # 2) merge outputs into JSON and delete per-iteration .out files
    merge_outputs_to_json(root, niterations, out_json_name=args.out_json)

    print("Done.")


if __name__ == "__main__":
    main()