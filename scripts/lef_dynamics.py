"""
1D Loop Extrusion Factor (LEF) dynamics engine.

Simulates cohesin-like loop extruding factors on a 1D lattice
with CTCF barriers (orientation-dependent stalling).

Based on:
  Fudenberg et al. 2016, Cell Reports 15:2038-2049
  Banigan & Mirny 2020, eLife 9:e53558

This module provides a pure-Python implementation that can be used
standalone or integrated with polychrom's 3D polymer simulation.
"""

import numpy as np
from typing import Optional, Tuple, List


class LEFSimulator:
    """
    Simulate loop-extruding factors (LEFs) on a 1D chromatin lattice.

    Each LEF has two arms that extrude in opposite directions.
    LEFs stall at CTCF sites (orientation-dependent) and at other LEFs.
    LEFs have a finite lifetime and are replaced upon unbinding.

    Parameters
    ----------
    N : int
        Number of lattice sites (monomers).
    n_lefs : int
        Number of LEFs to simulate.
    lifetime : float
        Mean LEF lifetime in simulation steps.
    ctcf_positions : array-like
        Positions of CTCF binding sites on the lattice.
    ctcf_orientations : array-like
        Orientation of each CTCF site: +1 (forward) or -1 (reverse).
        A +1 site blocks the left arm of a LEF (moving leftward).
        A -1 site blocks the right arm of a LEF (moving rightward).
    ctcf_capture : float
        Probability that a LEF stalls when encountering an occupied CTCF site.
    ctcf_release : float
        Probability per step that a stalled LEF releases from CTCF.
    rng_seed : int, optional
        Random seed for reproducibility.
    """

    def __init__(
        self,
        N: int,
        n_lefs: int,
        lifetime: float,
        ctcf_positions: np.ndarray,
        ctcf_orientations: np.ndarray,
        ctcf_capture: float = 0.9,
        ctcf_release: float = 0.003,
        rng_seed: Optional[int] = None,
    ):
        self.N = N
        self.n_lefs = n_lefs
        self.lifetime = lifetime
        self.ctcf_capture = ctcf_capture
        self.ctcf_release = ctcf_release
        self.rng = np.random.RandomState(rng_seed)

        # Build CTCF stalling maps
        # left_stall[i] = True if a left-moving arm should stall at position i
        # right_stall[i] = True if a right-moving arm should stall at position i
        self.left_stall = np.zeros(N, dtype=bool)
        self.right_stall = np.zeros(N, dtype=bool)

        for pos, ori in zip(ctcf_positions, ctcf_orientations):
            if 0 <= pos < N:
                if ori == +1:
                    self.left_stall[pos] = True   # blocks leftward-moving arm
                elif ori == -1:
                    self.right_stall[pos] = True  # blocks rightward-moving arm

        # LEF state arrays
        # Each LEF has: left_pos, right_pos, age, left_stalled, right_stalled
        self.left_pos = np.zeros(n_lefs, dtype=int)
        self.right_pos = np.zeros(n_lefs, dtype=int)
        self.ages = np.zeros(n_lefs, dtype=int)
        self.left_stalled = np.zeros(n_lefs, dtype=bool)
        self.right_stalled = np.zeros(n_lefs, dtype=bool)

        # Occupancy grid for mutual exclusion
        self.occupied = np.zeros(N, dtype=int)  # 0 = free, >0 = occupied by LEF arm

        # Initialize LEFs at random positions
        self._initialize_lefs()

    def _initialize_lefs(self):
        """Place all LEFs at random positions along the lattice."""
        for i in range(self.n_lefs):
            self._place_lef(i)

    def _place_lef(self, lef_idx: int):
        """Place a single LEF at a random position."""
        # Remove old occupancy
        if self.left_pos[lef_idx] != self.right_pos[lef_idx]:
            self.occupied[self.left_pos[lef_idx]] = 0
            self.occupied[self.right_pos[lef_idx]] = 0

        # Find a free position
        for _ in range(100):
            pos = self.rng.randint(0, self.N)
            if self.occupied[pos] == 0:
                adj = min(pos + 1, self.N - 1)
                if self.occupied[adj] == 0 or adj == pos:
                    self.left_pos[lef_idx] = pos
                    self.right_pos[lef_idx] = adj
                    self.ages[lef_idx] = 0
                    self.left_stalled[lef_idx] = False
                    self.right_stalled[lef_idx] = False
                    self.occupied[pos] = lef_idx + 1
                    if adj != pos:
                        self.occupied[adj] = lef_idx + 1
                    return

        # Fallback: place anyway (may overlap)
        pos = self.rng.randint(0, self.N - 1)
        self.left_pos[lef_idx] = pos
        self.right_pos[lef_idx] = pos + 1
        self.ages[lef_idx] = 0
        self.left_stalled[lef_idx] = False
        self.right_stalled[lef_idx] = False

    def step(self):
        """
        Advance the simulation by one step.

        Each LEF:
        1. Ages by 1 step
        2. Unbinds with probability 1/lifetime (or if past lifetime)
        3. Each arm attempts to extrude by 1 position
        4. Arms stall at CTCF sites or upon collision with other LEFs
        """
        unbind_prob = 1.0 / self.lifetime

        # Rebuild occupancy
        self.occupied[:] = 0
        for i in range(self.n_lefs):
            lp, rp = self.left_pos[i], self.right_pos[i]
            if 0 <= lp < self.N:
                self.occupied[lp] = i + 1
            if 0 <= rp < self.N:
                self.occupied[rp] = i + 1

        # Process each LEF
        perm = self.rng.permutation(self.n_lefs)
        for i in perm:
            self.ages[i] += 1

            # --- Unbinding ---
            if self.rng.random() < unbind_prob:
                # Remove and re-place
                self.occupied[self.left_pos[i]] = 0
                if self.right_pos[i] != self.left_pos[i]:
                    self.occupied[self.right_pos[i]] = 0
                self._place_lef(i)
                continue

            # --- Left arm extrusion (moves leftward) ---
            if not self.left_stalled[i]:
                new_left = self.left_pos[i] - 1
                if new_left >= 0 and self.occupied[new_left] == 0:
                    # Check CTCF stalling
                    if self.left_stall[new_left] and self.rng.random() < self.ctcf_capture:
                        self.left_stalled[i] = True
                    else:
                        self.occupied[self.left_pos[i]] = 0
                        self.left_pos[i] = new_left
                        self.occupied[new_left] = i + 1
                elif new_left >= 0 and self.occupied[new_left] != 0:
                    pass  # blocked by another LEF
                else:
                    pass  # at boundary
            else:
                # Stalled at CTCF — release with small probability
                if self.rng.random() < self.ctcf_release:
                    self.left_stalled[i] = False

            # --- Right arm extrusion (moves rightward) ---
            if not self.right_stalled[i]:
                new_right = self.right_pos[i] + 1
                if new_right < self.N and self.occupied[new_right] == 0:
                    if self.right_stall[new_right] and self.rng.random() < self.ctcf_capture:
                        self.right_stalled[i] = True
                    else:
                        self.occupied[self.right_pos[i]] = 0
                        self.right_pos[i] = new_right
                        self.occupied[new_right] = i + 1
                elif new_right < self.N and self.occupied[new_right] != 0:
                    pass  # blocked
                else:
                    pass  # at boundary
            else:
                if self.rng.random() < self.ctcf_release:
                    self.right_stalled[i] = False

    def get_bonds(self) -> List[Tuple[int, int]]:
        """Return current LEF-bridged monomer pairs as a list of (i, j) tuples."""
        bonds = []
        for i in range(self.n_lefs):
            lp, rp = self.left_pos[i], self.right_pos[i]
            if 0 <= lp < self.N and 0 <= rp < self.N and lp != rp:
                bonds.append((lp, rp))
        return bonds

    def get_bond_arrays(self) -> Tuple[np.ndarray, np.ndarray]:
        """Return bonds as two arrays: left_positions, right_positions."""
        bonds = self.get_bonds()
        if not bonds:
            return np.array([], dtype=int), np.array([], dtype=int)
        bonds = np.array(bonds)
        return bonds[:, 0], bonds[:, 1]

    def get_loop_sizes(self) -> np.ndarray:
        """Return array of current loop sizes (in monomers)."""
        sizes = self.right_pos - self.left_pos
        return sizes[sizes > 0]

    def run(self, n_steps: int) -> List[List[Tuple[int, int]]]:
        """
        Run simulation for n_steps and return bond trajectories.

        Returns
        -------
        trajectories : list of list of (int, int)
            Bond pairs at each step.
        """
        trajectories = []
        for _ in range(n_steps):
            self.step()
            trajectories.append(self.get_bonds())
        return trajectories
