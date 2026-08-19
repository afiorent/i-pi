#!/usr/bin/env python3
"""
Update RPQA without local minimization — choose best bead from min{iteration}.out.

Changes in this version:
 - When assembling the init file, the script now reads the previous simulation snapshot:
     simulation{iteration-1}.pos_{:02d}.xyz
   for all beads (except the best bead which still prefers min{iteration}.pos_{best_bead:02d}.xyz).
 - For iteration == 0 the previous-snapshot source is delocalization.pos_{:02d}.xyz as requested.
 - Prints the list of per-bead final potentials and the chosen best bead + energy.
 - No minimization is performed.

Usage:
  python update_RPQA.py --iteration 0 --nbeads 32 --root . --final-root . --template RPQA_template.xml
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
from jinja2 import Template
from ase.io import write

import ipi  # assumes your i-pi fork exposes read_output and read_trajectory


def parse_args():
    p = argparse.ArgumentParser(description="Update RPQA using min{iteration}.out to pick best bead.")
    p.add_argument("--iteration", "-i", type=int, required=True, help="current iteration index")
    p.add_argument("--nbeads", "-b", type=int, default=32, help="number of beads")
    p.add_argument("--hbar2-0", type=float, default=5.0, help="start hbar2")
    p.add_argument("--hbar2-1", type=float, default=1.0, help="end hbar2")
    p.add_argument("--N-iteration", type=int, default=10, help="number of hbar2 points")
    p.add_argument("--root", "-r", default=".", help="root directory where files live")
    p.add_argument("--template", "-t", default="RPQA_template.xml", help="jinja2 template path")
    p.add_argument("--final-root", "-o", default=".", help="where to write init and input files")
    p.add_argument("--simu-steps", type=int, default=10000, help="number of simulation steps (also used as QA_steps)")
    p.add_argument("--temperature", type=float, default=10.0, help="simulation temperature")
    return p.parse_args()


def read_min_output(min_out_path: str):
    """Read min{iteration}.out and return the first element if tuple, or dict directly."""
    ret = ipi.read_output(min_out_path)
    if isinstance(ret, tuple) and len(ret) >= 1:
        out = ret[0]
    else:
        out = ret
    if not isinstance(out, dict):
        raise RuntimeError(f"ipi.read_output returned unexpected type for {min_out_path}: {type(out)}")
    return out


def get_final_bead_potentials_from_output(output_dict: dict, expected_nbeads: int | None = None):
    """
    Extract final bead potentials using the explicit key 'bead_potentials' when available.
    Returns a 1D numpy array of final potentials per bead.
    """
    key = "bead_potentials"
    if key not in output_dict:
        raise KeyError(f"'{key}' not found in min output dictionary")
    arr = np.asarray(output_dict[key], dtype=float)
    if arr.ndim == 2:
        last = arr[-1, :]
    elif arr.ndim == 1:
        last = arr
    else:
        last = arr.reshape(-1)
    if expected_nbeads is not None and last.size != expected_nbeads:
        print(f"Warning: expected {expected_nbeads} bead potentials but found {last.size}", file=sys.stderr)
    return np.asarray(last)


def assemble_init_file(iteration: int, best_bead: int, nbeads: int, root: str, final_root: str):
    """
    Build init_{iteration}.xyz containing one frame per bead:
     - for bead != best_bead:
         * if iteration==0 use delocalization.pos_{:02d}.xyz (previous snapshot for iteration 0)
         * else use simulation{iteration-1}.pos_{:02d}.xyz (previous iteration snapshot)
     - for best_bead use last frame of min{iteration}.pos_{best_bead:02d}.xyz (preferred),
       falling back to the corresponding previous-simulation/delocalization file if min file missing.
    Returns (N_atoms, init_path, frames_list).
    """
    frames = []
    for bead in range(nbeads):
        if bead == best_bead:
            # prefer min{iteration} pos for best bead
            min_pos = os.path.join(root, f"min{iteration}.pos_{bead:02d}.xyz")
            if os.path.exists(min_pos):
                frames_min = ipi.read_trajectory(min_pos)
                atoms = frames_min[-1]
            else:
                # fallback to previous simulation snapshot logic
                if iteration == 0:
                    sim_pos = os.path.join(root, f"delocalization.pos_{bead:02d}.xyz")
                else:
                    sim_pos = os.path.join(root, f"simulation{iteration-1}.pos_{bead:02d}.xyz")
                if not os.path.exists(sim_pos):
                    raise FileNotFoundError(f"Neither {min_pos} nor {sim_pos} exist")
                frames_sim = ipi.read_trajectory(sim_pos)
                atoms = frames_sim[-1]
        else:
            # non-best beads: use previous simulation snapshot (iteration-1) or delocalization for iteration 0
            if iteration == 0:
                sim_pos = os.path.join(root, f"delocalization.pos_{bead:02d}.xyz")
            else:
                sim_pos = os.path.join(root, f"simulation{iteration-1}.pos_{bead:02d}.xyz")
            if not os.path.exists(sim_pos):
                raise FileNotFoundError(f"{sim_pos} not found")
            frames_sim = ipi.read_trajectory(sim_pos)
            atoms = frames_sim[-1]
        frames.append(atoms.copy())

    sample = frames[0]
    if hasattr(sample, "get_global_number_of_atoms"):
        N = sample.get_global_number_of_atoms()
    else:
        N = len(sample)

    init_name = os.path.join(final_root, f"init_{iteration}.xyz")
    write(init_name, frames)
    return N, init_name, frames


def build_fixbeads(N: int, best_bead: int):
    return (N * best_bead + np.arange(N)).tolist()


def build_hbar2_pair(i: int, N_iteration: int, hbar2_start: float, hbar2_end: float):
    sqrt_vals = np.linspace(np.sqrt(hbar2_start), np.sqrt(hbar2_end), N_iteration)
    hbar2_vals = sqrt_vals ** 2
    if i < 0:
        i = 0
    if i >= N_iteration:
        i = N_iteration - 1
    j = i + 1
    if j >= N_iteration:
        j = N_iteration - 1
    return float(hbar2_vals[i]), float(hbar2_vals[j])


def main():
    args = parse_args()
    iteration = args.iteration
    nbeads = args.nbeads
    final_root = args.final_root
    root_dir = args.root
    template_path = args.template
    simu_steps = args.simu_steps
    temperature = args.temperature

    if simu_steps <= 0:
        print("Error: --simu-steps must be > 0", file=sys.stderr)
        sys.exit(1)

    # 1) read min{iteration}.out and extract bead potentials
    min_out = os.path.join(root_dir, f"min{iteration}.out")
    if not os.path.exists(min_out):
        print(f"Error: {min_out} not found", file=sys.stderr)
        sys.exit(1)

    out_dict = read_min_output(min_out)

    try:
        final_pots = get_final_bead_potentials_from_output(out_dict, nbeads)
    except KeyError as e:
        print(f"Error extracting bead potentials: {e}", file=sys.stderr)
        sys.exit(1)

    # Print all bead energies and choose the best bead
    energies_list = [float(x) for x in final_pots]
    best_bead = int(np.argmin(final_pots))
    best_pot = float(final_pots[best_bead])

    print("Bead final potentials (index: energy):")
    for idx, val in enumerate(energies_list):
        mark = "*" if idx == best_bead else " "
        print(f"{mark} [{idx:2d}] {val:.12g}")

    print(f"Chosen best bead: {best_bead} with final potential = {best_pot:.12g}")

    # 2) assemble init file (using min pos for best bead, previous simulation snapshots otherwise)
    N, init_path, frames = assemble_init_file(iteration, best_bead, nbeads, root_dir, final_root)
    print(f"Wrote init file: {init_path} (N_atoms={N})")

    # 3) compute fixbeads
    fixbeads = build_fixbeads(N, best_bead)

    # 4) pick hbar2 pair
    hbar2_0, hbar2_1 = build_hbar2_pair(iteration, args.N_iteration, args.hbar2_0, args.hbar2_1)

    # 5) render template
    tpl_text = open(template_path).read()
    tpl = Template(tpl_text)
    rendered = tpl.render(
        simu_steps = simu_steps,
        QA_steps = simu_steps,
        temperature = temperature,
        init_file=os.path.basename(init_path),
        nbeads=nbeads,
        hbar2_0=hbar2_0,
        hbar2_1=hbar2_1,
        sqrt_scale=True,
        stride=1000,
        fixbeads=fixbeads,
        best_bead_index=best_bead,
        best_energy_eV=best_pot,
        best_positions=frames[best_bead].get_positions().tolist(),
        n_atoms=N,
        fsimulation=f"simulation{iteration}",
    )

    os.makedirs(final_root, exist_ok=True)
    out_path = os.path.join(final_root, f"input_{iteration}.xml")
    with open(out_path, "w") as f:
        f.write(rendered)

    print("fixbeads:", fixbeads)
    print(f"hbar2 pair used: {hbar2_0:.6e}, {hbar2_1:.6e}")
    print(f"Wrote rendered template to {out_path}")


if __name__ == "__main__":
    main()