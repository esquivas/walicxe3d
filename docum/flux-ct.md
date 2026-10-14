# Constrained transport (flux-CT) for MHD with AMR in Walicxe3D

Developer notes on the face-centred constrained-transport implementation of
the MHD solver in Walicxe3D. They cover the algorithm, how it is mapped onto
the block-based AMR, the checks built into the code, the verification tests,
and the open issues.

*Status (14 Oct 2026):* implemented and tested in the project directory
`flux-ct/` (branch `dev`, commit `2ecbc9b`). The routines live in local
copies of `amr.f90`, `boundary.f90`, `flux_ct.f90`, `godunov.f90`,
`hlld.f90`, `ghosts.f90`, `diagnostics.f90`, `main.f90`, ..., which
override `source/` through `VPATH`. They have not been moved to `source/`
yet; see [Pending work](#11-known-issues-and-pending-work).

---

## Contents

1. [Why constrained transport](#1-why-constrained-transport)
2. [Data layout](#2-data-layout)
3. [Time integration (VL+CT)](#3-time-integration-vlct)
4. [Ghost cells and ghost faces](#4-ghost-cells-and-ghost-faces)
5. [EMF synchronisation between blocks](#5-emf-synchronisation-between-blocks)
6. [Refinement and coarsening of the face field](#6-refinement-and-coarsening-of-the-face-field)
7. [Domain boundary conditions for the faces](#7-domain-boundary-conditions-for-the-faces)
8. [Parameters and build](#8-parameters-and-build)
9. [Built-in checks and diagnostics](#9-built-in-checks-and-diagnostics)
10. [Verification tests](#10-verification-tests)
11. [Known issues and pending work](#11-known-issues-and-pending-work)
12. [Map of routines](#12-map-of-routines)
13. [Other changes to the MHD/AMR code made during this work](#13-other-changes-to-the-mhdamr-code-made-during-this-work)
14. [References](#14-references)

---

## 1. Why constrained transport

The MHD solver (HLLD Riemann solver, second-order MUSCL–Hancock-type
predictor–corrector) originally kept $\nabla\cdot\mathbf{B}$ under control
with Powell's 8-wave source terms. The next attempt was flux-CD
(field-interpolated constrained transport on cell-centred B, as in the
uniform-grid code Guacho). On a uniform mesh flux-CD keeps
$\nabla\cdot\mathbf{B}$ at roundoff. With AMR, however, every refinement
or coarsening of cell-centred B creates divergence that flux-CD then
*preserves* (it conserves the discrete divergence, whatever its value).
That produced persistent stripes and seams in the pressure at former
level jumps.

The standard cure is to store the **face-normal** components of
$\mathbf{B}$ and evolve them with **edge-centred EMFs** (Evans & Hawley
1988). The discrete divergence of each cell is then conserved exactly by
construction. With AMR it also needs a divergence-free prolongation and
restriction (Tóth & Roe 2002, Balsara 2001) and consistent EMFs at level
jumps. The implementation follows Athena (Gardiner & Stone 2005, 2008;
Stone et al. 2008; Stone & Gardiner 2009), adapted to Walicxe's cubic
blocks and its data structures.

---

## 2. Data layout

| Array | Shape | Meaning |
|---|---|---|
| `BF(nb,c,i,j,k)` | `(nbMaxProc, 3, nxmin:nxmax, nymin:nymax, nzmin:nzmax)` | Face-normal B, component `c`, at face `i+1/2` (for `c=1`), `j+1/2` (`c=2`), `k+1/2` (`c=3`) |
| `BFP` | same | Half-step (predictor) copy |
| `U(:,6:8,...)`, `UP(:,6:8,...)` | as usual | Cell-centred B = **average of the two faces** (Athena eqs. 26–28) |
| `EMF(e,i,j,k,nb)` | `(3, 0:n, 0:n, 0:n, nbMaxProc)` | Edge EMFs, "right" indexing like the fluxes (below) |
| `EW` | same as `EMF` | Weights used while synchronising the EMFs (§5) |
| `ER(e,i,j,k)` | `(3, 0:n+1, ...)` | Cell-centred reference field $-\mathbf{v}\times\mathbf{B}$ of the block being processed |

* `BF` uses the **same index convention as the fluxes** `FC/GC/HC`: index
  `i` is the face between cells `i` and `i+1`. Index `0` is the block's own
  low face and index `n` its high face (`n = ncells_x`).
* EMF positions:
  * `EMF(1,i,j,k)` = $E_x$ at $(i,\ j+\tfrac12,\ k+\tfrac12)$, `i=1..n`, `j,k=0..n`
  * `EMF(2,i,j,k)` = $E_y$ at $(i+\tfrac12,\ j,\ k+\tfrac12)$
  * `EMF(3,i,j,k)` = $E_z$ at $(i+\tfrac12,\ j+\tfrac12,\ k)$
* **Invariant:** after every step, $U_{6} = \tfrac12(B_{x,i-1/2}+B_{x,i+1/2})$,
  and likewise for $U_7, U_8$, in every physical cell. The thermal energy
  absorbs any change of $B^2/2$ whenever the cell-centred B is recomputed
  from the faces (refinement, coarsening), so total energy is conserved.
* **Initial conditions** are set from a vector potential (face B from
  differences of $A$ along the face edges), so the initial face divergence
  is zero to roundoff. Examples: `orzag-tang.f90`, `field_loop.f90`.
* **Cells must be cubic** ($\Delta x=\Delta y=\Delta z$), and so must the
  blocks (`ncells_x = ncells_y = ncells_z`). The prolongation and several
  helper routines rely on this.

Memory: `BF`, `BFP`, `EMF` and `EW` are included in the per-block memory
estimate `block_arrays` in `parameters.f90`, which sets `nbMaxProc`.

---

## 3. Time integration (VL+CT)

Predictor–corrector of Stone & Gardiner (2009), as in `godunov.f90`:

```
call boundary_bf (1, BF, U)                    ! ghost faces at time n
stage 1 (dt/2, first-order fluxes):
   for each block: fluxes (HLLD); upwindStep; flux_ct_emf
   sync_emf                                    ! all blocks agree on edge EMFs
   for each block: flux_ct_update(dt/2)        ! BFP = BF + CT update, UP(6:8)
call boundary_bf (1, BFP, UP)                  ! ghost faces at n+1/2
stage 2 (dt, second-order fluxes from the half-step state):
   for each block: store_face_fluxes; fluxes; upwindStep; flux_ct_emf
   sync_emf
   for each block: flux_ct_update(dt)          ! BFP = BF(n) + CT update with stage-2 EMFs
   flux_ct_advance                             ! BF = BFP
apply_flux_correction                          ! hydro refluxing, B excluded
```

### 3.1 Normal B in the Riemann problem
In HLLD the normal field at an interface is not reconstructed. The face
value is used for both states: `pl(6) = pr(6) = BF(...)` in stage 1 and
`BFP` in stage 2, for the x, y and z sweeps (`hlld.f90`).

The GS05 EMFs need fluxes on **one extra transverse layer**, so the HLLD
loops run over `1-e .. ncells+e` in the transverse directions (`e = 1` when
flux-CT is on). That in turn needs edge and corner ghost cells (§4).

### 3.2 Edge EMFs (Gardiner & Stone 2005)
`flux_ct_emf(locIndx)`:

* Face EMFs come from the induction fluxes:
  * x-faces: $E_z=-F(B_y)$, $E_y=+F(B_z)$
  * y-faces: $E_z=+G(B_x)$, $E_x=-G(B_z)$
  * z-faces: $E_y=-H(B_x)$, $E_x=+H(B_y)$
* The cell reference field is $\mathbf{E}^r=-\mathbf{v}\times\mathbf{B}$ from `PRIM`.
* Each edge EMF (`emf_edge`) is the average of the four surrounding face
  EMFs plus four upwinded "face minus cell" corrections, chosen by the
  sign of the mass flux on the corresponding face (`upwind`; the average
  is used when the mass flux is 0). These are Athena eqs. 79–81.

### 3.3 CT update
`flux_ct_update(locIndx, dtp)` updates all faces of the block (Athena
eqs. 16–18), always starting from `BF` at time $n$:

$$B_{x,i+1/2}^{*} = B_{x,i+1/2}^{n} - \frac{\Delta t}{\Delta y}\left(E_{z,j+1/2}-E_{z,j-1/2}\right) + \frac{\Delta t}{\Delta z}\left(E_{y,k+1/2}-E_{y,k-1/2}\right)$$

and cyclically for $B_y$ and $B_z$. It then sets `UP(6:8)` to the face
averages. Each EMF enters the two faces that share it with opposite signs,
so the divergence of every cell is preserved exactly, provided every copy
of a shared edge EMF is identical in all blocks. §5 makes sure of that.

### 3.4 Flux correction
Hydro refluxing at fine–coarse faces (`flux_correction.f90`) skips the
magnetic components (`dF(6:8) = 0` when flux-CT is on). For B, the
equivalent of refluxing is the EMF restriction of §5.

---

## 4. Ghost cells and ghost faces

### 4.1 Edge and corner ghost cells (U)
The GS05 EMFs and the extra transverse flux layer need edge and corner
ghost cells. `normalBoundary` already exchanges directions in order
(x, then y, then z). With the optional argument `edges=.true.` of
`layerLimits`, the transverse range of the y- and z-layers is extended over
the ghost cells filled earlier. Edge and corner ghosts between same-level
blocks and at non-periodic boundaries are then filled without diagonal
messages.

At level jumps, `interpGhost` (limited linear coarse→fine) fills the face
ghosts. `extrapolate_edges` (zero gradient) then fills the adjacent edge
and corner ghosts.

### 4.2 Ghost faces (`boundary_bf(depth, bfield, uvars)`)
It uses the same x→y→z sweep and the same index sets as the cells:

* **Same level (including periodic):** the neighbour's layer is copied.
  The low ghost layer includes index `0`, the block's own low face, which
  is overwritten by the neighbour's copy. The difference found is recorded
  in `bf_shared_err` (a check that should be exactly 0).
* **Level jump:** the tangential ghost faces get the average of the two
  adjacent cell-centred values (`faces_from_cells`). The normal face at
  index 0 is the block's own face and is not touched.
* **Non-periodic domain boundary:** `bc_faces` (§7).

Ghost faces are used by the HLLD normal field on the extra transverse
layer and by the transverse slopes of the prolongation (§6).

### 4.3 Ghost-refresh hook
`refineZone` (`amr.f90`) refines one level per pass. The blocks it creates
have no ghost cells or ghost faces when the next pass refines them. `amr`
cannot call `boundary` or `boundary_bf` directly, because that would be a
circular module dependency:
`amr` ← `boundaries` → `userconds` (`user.f90`) → user modules → `amr`.

The solution is a procedure pointer in `amr`:

```fortran
procedure(ghost_refresher), pointer :: refresh_ghosts => null()
```

`main.f90` registers it right after `initmain`:
`refresh_ghosts => refresh_ghosts_all` (module `ghosts`, file `ghosts.f90`,
compiled after `boundary.o` and `flux_ct.o`). `refresh_ghosts_all` does
`boundary(1,U)`, `boundary_bf(1,BF,U)` and `calcPrimsAll(U,PRIM,CELLS_GHOST)`.
`refineZone` calls it at the start of each pass when it is associated. The
main loop calls `refresh_ghosts_all` directly before the diagnostics and
`admesh`. Projects that do not register the hook keep the old behaviour.

---

## 5. EMF synchronisation between blocks

Every edge is shared by up to four blocks, and each computes its own copy
of the EMF from its own fluxes and ghost data. All copies must be
**identical**; otherwise the faces shared by two blocks drift apart, and
the overwrite in `boundary_bf` turns the difference into divergence. At
level jumps the coarse EMFs must also equal the average of the fine ones,
so that each coarse face stays equal to the average of its four fine faces
(the B analogue of refluxing; Tóth & Roe 2002, sect. 6).

`sync_emf` (called by all ranks after `flux_ct_emf` on every block, in
both stages) does:

```
EW = 1 on all local edges
for lev = maxlev, ..., 1:
    sync_same_level (SYNC_AVERAGE, lev)
    sync_same_level (SYNC_MAXW,    lev)
    if lev > 1: restrict_fine_emf (lev-1)
```

* **`sync_same_level(mode, lev)`:** pairs every same-level face neighbour
  (high face of *a* with low face of *b*) in an x, then y, then z sweep.
  It exchanges the tangential EMFs of that face together with their weights
  (`pack_face_emf`, MPI_SENDRECV), and combines them with
  `combine_face_emf`, which is symmetric: both owners call it in the same
  (a, b) order and get bit-identical results.
  * `SYNC_AVERAGE`: weighted average, weights added. The pairs of one pass
    always merge disjoint sets of copies, so with four blocks around an
    edge all copies end with the full average (weight 4).
  * `SYNC_MAXW`: the copy with the larger weight wins.
* **Why both passes.** At an inner corner of a finer region an edge is
  shared by **three** same-level blocks and one coarser block. The two
  fine blocks that are only diagonal neighbours never exchange directly.
  After the averaging pass one copy has weight 2 and the other two weight
  3; the max-weight pass completes it. (This configuration first appeared
  in the dynamic-AMR test r5 and produced $\nabla\cdot\mathbf{B}$ at the
  block corners.)
* **`restrict_fine_emf(clev)`:** for every coarse face of a level-`clev`
  block whose neighbours are four finer blocks, each fine block averages
  its two collinear fine edge EMFs onto the coarse edges of its quadrant.
  The coarse block takes them with weight `WFINE = 1000`, which wins over
  any same-level value in later passes.
* **Level by level, finest first.** A level is final before it is
  restricted, so restriction cascades correctly with three or more levels.
* With equal weights the arithmetic equals a plain average ((a+b)/2, then
  (2a'+2b')/4), so uniform-mesh runs are bit-identical to the
  unsynchronised scheme.

---

## 6. Refinement and coarsening of the face field

### 6.1 Prolongation (Tóth & Roe 2002, sect. 3)
`refineBlock` calls `prolong_cell_faces` for every parent cell. Indices
follow T&R: along the normal $-2, 0, +2$ are the low parent face, the fine
mid-plane face and the high parent face; $\pm1$ are the transverse fine
cells. For cubic cells all face areas are equal, so the formulas, written
in fluxes in the paper, hold directly for field values.

* **Exterior fine faces (eq. 4).** Each parent face gets limited
  transverse slopes (`average(..., limiter_type)`, same limiter as the
  hydro), using the father's neighbouring faces including its ghost faces:
  $$u^{\pm2,j,k} = U^{\pm2} + \tfrac14\left(j\,\partial_y U + k\,\partial_z U\right),\quad j,k=\pm1,$$
  so the four fine faces average exactly to the parent face (eq. 3).
* **Interior fine faces (eqs. 8–12)** (`tr_interior`):
  $$u^{0,j,k} = \tfrac12\left(u^{2,j,k}+u^{-2,j,k}\right) + U_{xx} + k\,\Delta z^2 V_{xyz} + j\,\Delta y^2 W_{xyz}$$
  with $U_{xx}=\tfrac18\sum_{ijk}\left(ij\,v^{i,2j,k}+ik\,w^{i,j,2k}\right)$
  and $\Delta y^2 U_{xyz}=\tfrac1{16}\sum_{ijk} ijk\,u^{2i,j,k}$ (cubic
  cells), and cyclically for $v$ and $w$. All eight fine cells get the
  divergence of the parent cell (zero). The third-order terms vanish for
  eq.-4 exterior faces but not for copied ones (6.3).
* The children's cell-centred B is recomputed from their faces
  (`faces_to_cells`), followed by `calcPrimsBlock`. Mass, momentum and
  total energy are interpolated conservatively as before (limited linear),
  so the thermal energy absorbs the change in $B^2/2$.
* The father's ghost faces must be current. `refresh_ghosts_all` is called
  before `admesh` in the main loop, and through the hook in `refineZone`.

### 6.2 Restriction on coarsening
`coarseBlocks`: each child averages each group of four fine faces onto the
father's faces of its octant (`restrict_child_faces`, eq. 3). Remote
children send them by MPI right after their U data (same tag, same rank
pair, so the messages arrive in order) and the father stores them with
`put_child_faces`. The father's `U(6:8)` is then recomputed from its faces
before `calcPrimsBlock`. No fix-up is needed: the father's faces equal the
neighbour's faces (same level) or their restriction (coarse–fine), because
the children's faces already satisfied that.

### 6.3 Fix-up after refinement (`fix_new_faces`)
A new block takes its exterior faces from the father (eq. 4). T&R
require instead that, where finer data already exists, the fine faces be
copied from it. `admesh` (and each pass of `refineZone`) records the list
of fathers about to be refined (global, from `refineFlags`). After the
refinement and `syncBlockLists`, and before coarsening and load balancing,
`fix_new_faces(fathers, nf)` handles three cases:

1. **New block next to an old same-level block.** The shared face plane
   is replaced by the old block's (evolved) faces.
2. **Two new blocks with different fathers.** Their copies are averaged.
   Each group of four fine faces still averages to the father's face, so
   the result is consistent.
3. **New block next to finer blocks.** With three or more levels, refining
   a level-$L{+}1$ block can force (2:1 balance) the refinement of its
   level-$L$ neighbour in the same pass. The new level-$(L{+}1)$ child then
   faces new level-$(L{+}2)$ blocks. Its face must take the average of each
   2×2 group of their faces, which is exactly the old fine block's faces,
   not the values prolonged from its own father. (Missing this case gave
   a coarse–fine mismatch, and later divergence, in the 3-level test r7.)

Afterwards the 12 interior faces of every parent cell of every new block
are recomputed with eqs. 8–12 (`reset_interior_faces`, now with nonzero
third-order terms), followed by `faces_to_cells` and `calcPrimsBlock`.
For cells whose exterior faces did not change, this reproduces the
prolonged values bit for bit. The routine logs the largest face change
found (`fix_new_faces: it, max change of shared faces`).

A new block's finer neighbours are always new themselves. A father next
to an *old* block two levels finer would violate 2:1 balance.

---

## 7. Domain boundary conditions for the faces

`bc_faces` (called from `boundary_bf` for `NEIGH_BOUNDARY`) fills the same
index set as the same-level copy (`layerLimits` with edges):

| BC | Tangential faces | Normal ghost faces |
|---|---|---|
| `BC_FREE_FLOW` | copy of the boundary cell | zero gradient (copy of the wall face) |
| `BC_REFLECTIVE` | copy of the mirror cell | antisymmetric about the wall: $B_n(l)=-B_n(2l_w-l)$ (perfectly conducting wall, as in Athena) |
| `BC_PERIODIC` | handled as same-level neighbours | |

The wall face itself (index 0 on the low side, `ncells` on the high side)
is the block's own face, evolved by CT, and is never modified. Other BC
types abort with a message.

Two consistent changes were made to the cell boundaries in
`normalBoundary` (`boundary.f90`):
* **Reflective:** the normal component of the cell-centred B is flipped
  as well (previously only the normal velocity was).
* **Free flow:** the high sides (RIGHT/BACK/TOP) used a mirror index
  (`2n-i+1`), so with two ghost layers the outer ghost received cell
  `n-1`. They now copy cell `n`, as the low sides already did (true zero
  gradient). This also affects pure-hydro runs.

User boundary conditions (`userBoundary`) act on cells only. Any user
condition that sets B must also set the faces consistently. This is not
handled yet.

---

## 8. Parameters and build

`parameters.f90` (project copy in `flux-ct/`):

| Parameter | Meaning |
|---|---|
| `enable_flux_ct = .true.` | Face-centred CT (requires `BFIELD = Y` in the Makefile). Mutually exclusive with `eight_wave` (checked at start-up). |
| `enable_flux_correction` | Hydro refluxing at fine–coarse faces (B excluded under flux-CT) |
| `enable_cons_diagnostics` | Global conservation and divergence diagnostics (§9) |
| `limiter_type` | Limiter for the reconstruction **and** for the slopes of the face prolongation (e.g. `LIMITER_WOODWARD` = MC) |
| `refineThres`, `coarseThres` | Refinement thresholds (relative jumps of $\rho$ and $p$) |
| `refine_passive` | Also refine on **absolute** jumps of the first passive scalar (used by the field-loop test, where $\rho$ and $p$ are uniform) |
| `ncells_block` | Cells per block side; must give cubic blocks and cells |

Makefile: `flux_ct.o` after `hlld.o` and before `godunov.o`; `ghosts.o`
after `flux_ct.o`; `diagnostics.o` before `main.o`; `flux_correction.o`
after `boundary.o`. The test problems are user modules
(`MODULES_USER = orzag-tang.o field_loop.o`).

---

## 9. Built-in checks and diagnostics

`conservation_diagnostics` (every `diag_every` steps) writes
`diagdir/diagnostics.dat` with these columns: `it, time, mass, passive,
px, py, pz, E, Bx, By, Bz, Emag, <|divB|dx>/<|B|>, max|divB|dx`, and the
number of blocks per level. With flux-CT the divergence is the **face**
divergence $\sum (B_{f+}-B_{f-})/\Delta$. It should stay at roundoff
($\sim10^{-18}$ for $B\sim10^{-3}$, i.e. $\sim10^{-15}$ relative).

The log has one line with three checks:

| Check | Meaning | Expected |
|---|---|---|
| \|B_cell − ⟨B_face⟩\| | `U(6:8)` vs average of the two faces, all cells | exactly 0 |
| shared-face diff | largest difference between the two copies of a same-level shared face, found by `boundary_bf` since the previous report | 0 (≲1e-19 with 3 levels) |
| cf | coarse face vs average of its four fine faces (`check_cf_faces`) | roundoff |

`fix_new_faces` logs the largest face change at each refinement. A
temporary `z-symmetry` line (max \|vz\|, max \|Bz\|) and `zvar_check` are
in place for the investigation described in §11.

Note for post-processing: a "corner" divergence computed from cell-centred
B within a block equals the average of the face divergence of the
surrounding cells only if **all three** terms are included. Leaving out
$\partial_z B_z$ in a 2D slice measures $-\partial_z B_z$ instead.

---

## 10. Verification tests

Field loop (Gardiner & Stone 2005): $\rho=p=1$, $\mathbf{v}=(2,1,0)$,
$A_z=A_0\max(R_0-r,0)$, $A_0=10^{-3}$, $R_0=0.3$, box $2\times1\times0.5$,
periodic, finest grid 128×64×32, MC limiter, $t_{\rm fin}=2$ (two
crossings). Orszag–Tang for the uniform-mesh solver.

| Run | Setup | Result |
|---|---|---|
| OT uniform | 1 level | face divB at roundoff, mismatch 0, shared-face 0, energy conserved to ~1e-14 |
| loop, uniform | 1 level, minmod / MC | $E_{\rm mag}(2)/E_{\rm mag}(0)$ = 0.55 (minmod), 0.85 (MC); the small core asymmetry is limiter/upwind behaviour (reversed-flow run identical to 1e-15) |
| r2 | static mesh, half the box refined, IC re-imposed | divB roundoff on both levels, cf ≤ 2e-18, mass and energy ~1e-15; first test of the EMF sync and restriction |
| r3 | as r2, refined half **prolonged** (T&R) | at $t=0$: cf exactly 0, divB 2e-19; $E_{\rm mag}$ 1% lower than r2 (smoother field) |
| r4 | r3 + fix-up | bit-identical to r3 (all neighbours new, no change) |
| r5 | dynamic AMR (passive-scalar criterion) | **failed:** divB at block corners, from the 3-blocks-around-an-edge configuration → weighted EMF sync (§5) |
| r6 | r5 with the weighted sync | divB ~5e-15 relative through 16 refine/coarsen events, flat; shared ≤ 8e-21, cf ≤ 3e-20; $E_{\rm mag}$ ratio 0.847 |
| r7 | 3 levels (8³ blocks), dynamic | **failed:** cf mismatch at forced refinements → divB later; case 3 of the fix-up (§6.3) |
| r7b | r7 with the fix | max divB·dx 5.2e-18 over the run, cf ≤ 1.3e-18, mass/energy ~1e-16/1e-15, $E_{\rm mag}$ ratio 0.79 |
| r7c | r7b with frozen mesh | max \|vz\| and \|Bz\| exactly 0 (see §11) |
| r8 | 2 levels, free flow in x, reflective in z, loop leaves the box | $E_{\rm mag}$ → 1e-19 with no reflection; divB ≤ 5e-18 during the crossing; shared and cf ≤ 7e-21; mesh fully coarsens afterwards |

Lessons from the failures. Each one was caught by a built-in check
**before** it showed up in the solution:
* A missing call to `flux_ct_update` and a missing `+e` loop bound in HLLD
  were caught by the B_cell − ⟨B_face⟩ check.
* The 3-block edge (r5) was caught by the shared-face check, which was
  then turned into divergence by the overwrite in `boundary_bf`.
* The forced refinement next to finer blocks (r7) was caught by the cf
  check, at the same step as the event.

---

## 11. Known issues and pending work

* **z-asymmetry in dynamic 3-level runs (under investigation).** In a
  problem that is uniform in z, $v_z$ and $B_z$ should stay exactly 0. With
  a frozen 3-level mesh they do (r7c). With dynamic refinement they appear
  right after the first refinement event and grow to $|B_z|/B_0\sim5\times10^{-6}$
  and $|v_z|\sim10^{-7}$ by $t=2$ (r7b). With 2 levels the effect is much
  smaller ($\sim3\times10^{-11}$ relative). The face divergence is not
  affected. `zvar_check` (temporary) logs the largest z-variation after
  each stage of `admesh` and after load balancing, to identify the
  operation responsible (run r7d).
* **Move to `source/`:** rename `enable_flux_cd` → `enable_flux_ct` in all
  projects, add `ghosts.o` and `flux_ct.o` to their Makefiles, register the
  hook in `main.f90`, wrap the BFIELD-only lines in `#ifdef BFIELD`.
* **Orszag–Tang with dynamic AMR**, the first nonlinear MHD test with the
  full machinery.
* **User boundary conditions and sources that modify B** must update the
  faces as well.
* **Performance:** the synchronisation passes loop over all global blocks
  with `neighbors()` calls (O(N) each); `getOwner`/`find` searches; and
  blocking MPI in the boundary exchange. Candidates: a bID → owner/index
  lookup table, non-blocking exchanges, and optimised compiler flags
  (development uses `-O3 -g -traceback -check bounds -warn all`).
* `faces_from_cells` at level jumps is first order. It only affects ghost
  faces (HLLD normal field on the extra layer, prolongation slopes), not
  the divergence.

---

## 12. Map of routines

**`flux_ct.f90`** (module `flux_ct`, only under `BFIELD`)

| Routine | Purpose |
|---|---|
| `flux_ct_emf` | Reference field and GS05 edge EMFs of one block |
| `emf_edge`, `upwind` | Edge EMF from four face EMFs with upwinded corrections |
| `flux_ct_update` | CT update `BFP = BF + ...`, `UP(6:8)` = face averages |
| `flux_ct_advance` | `BF = BFP` |
| `sync_emf` | Level-by-level synchronisation (§5) |
| `sync_same_level` | Same-level pairwise exchange, average or max-weight |
| `combine_face_emf` | Symmetric combination of two packed face copies |
| `restrict_fine_emf` | Fine → coarse edge EMFs, weight `WFINE` |
| `pack_face_emf`, `unpack_face_emf`, `get_face`, `put_face` | Tangential EMFs (+ weights) of one face |
| `boundary_bf` | Ghost faces: same level, level jumps, domain boundaries |
| `faces_from_cells` | Tangential ghost faces from cell values (level jumps) |
| `bc_faces` | Free-flow / reflective ghost faces |
| `own_low_face`, `normal_face` | Helpers for the checks |
| `check_cf_faces` | Coarse face vs average of fine faces |

**`amr.f90`**

| Routine | Purpose |
|---|---|
| `refresh_ghosts` (pointer) | Ghost-refresh hook for `refineZone` |
| `prolong_cell_faces` | Exterior faces (eq. 4) + interior (`tr_interior`) of one parent cell |
| `tr_interior` | T&R eqs. 8–12 |
| `reset_interior_faces` | Interior faces of one parent cell from the block's own exterior faces |
| `fix_new_faces` | Post-refinement fix-up (§6.3) |
| `get_plane`, `put_plane` | Normal-B plane of a block face |
| `faces_to_cells` | `U(6:8)` from the faces |
| `restrict_child_faces`, `put_child_faces` | Restriction on coarsening |
| `zvar_check` | Temporary z-symmetry probe |

**Others:** `ghosts.f90` (`refresh_ghosts_all`); `boundary.f90`
(`layerLimits(..., edges)`, `extrapolate_edges`, reflective/free-flow cell
BCs); `hlld.f90` (face normal B, extra transverse layer); `godunov.f90`
(call sequence of §3); `doBalance` in `amr.f90` (transfers `BF` with `U`);
`diagnostics.f90`.

---

## 13. Other changes to the MHD/AMR code made during this work

These were found while validating MHD with Orszag–Tang (Oct 2026) and are
committed in `source/` or the project directories:

* **`hydro_core.f90`:**
  * Fixed parentheses in `cfast`, and added `max(0, ...)` guards in
    `cfast`/`cfastX`.
  * `calcPrims` uses `intent(inout)` and correct ghost bounds.
  * `average` (limiter) is now a module function. The Van Albada limiter
    is fixed (0 when the slopes have opposite signs).
  * New `limited_slope`.
* **`hlld.f90`:**
  * The degenerate $U_R^*$ branch used left-state variables.
  * The degeneracy test is now relative (`abs(den) < 1e-8*pst`).
  * Removed obsolete `PASSIVES` directives.
  * Split the flux loops per sweep direction.
* **`amr.f90`:**
  * `neighborLevel` now wraps periodic boundaries. Periodic BCs were
    previously never applied across the box.
  * `getOwner` returns `MPI_PROC_NULL` when a block is not found.
  * The refinement criterion uses relative jumps.
  * Refined children use limited linear interpolation instead of
    injection, with a positivity fallback.
* **`boundary.f90`:**
  * Coarse → fine ghost cells by limited linear interpolation
    (`coarseBox`, `interpGhost`).
  * Guards for `MPI_PROC_NULL`.
* **`flux_correction.f90` (new):** fine–coarse refluxing (mass conserved to
  roundoff with AMR).
* **`diagnostics.f90` (new):** conservation and divergence diagnostics.
* **`main.f90`:** checkpoint on demand (create `datadir/checkpoint` to
  write output 999); ghost refresh before `admesh`.

---

## 14. References

* Balsara, D. S. 2001, *J. Comput. Phys.* 174, 614 — divergence-free adaptive mesh refinement for MHD.
* Evans, C. R. & Hawley, J. F. 1988, *ApJ* 332, 659 — constrained transport.
* Gardiner, T. A. & Stone, J. M. 2005, *J. Comput. Phys.* 205, 509 — unsplit Godunov + CT, upwind EMFs (GS05).
* Gardiner, T. A. & Stone, J. M. 2008, *J. Comput. Phys.* 227, 4123 — 3D extension.
* Miyoshi, T. & Kusano, K. 2005, *J. Comput. Phys.* 208, 315 — HLLD.
* Stone, J. M. et al. 2008, *ApJS* 178, 137 — Athena.
* Stone, J. M. & Gardiner, T. A. 2009, *New Astron.* 14, 139 — VL+CT integrator.
* Tóth, G. & Roe, P. L. 2002, *J. Comput. Phys.* 180, 736 — divergence- and curl-preserving prolongation and restriction.
