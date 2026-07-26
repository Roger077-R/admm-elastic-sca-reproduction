# ADMM Elastic — SCA 2016 Code Guide

This repository contains the sample code for the 2016 SCA paper
**[ADMM ⊇ Projective Dynamics: Fast Simulation of General Constitutive Models](https://diglib.eg.org/server/api/core/bitstreams/f588facd-3719-44cd-8fbd-2d756af18031/content)**
by Rahul Narain, Matthew Overby, and George E. Brown.

The project applies the alternating direction method of multipliers (ADMM) to
implicit integration of deformable bodies. It retains the parallel local/global
structure of projective dynamics while supporting nonlinear elastic energies and
hard constraints. The included programs demonstrate cloth, volumetric elasticity,
external forces, and collisions.

> [!IMPORTANT]
> This is the original SCA 2016 implementation. It is useful for studying the
> paper and reproducing its examples, but it is legacy research code. The later
> [TVCG implementation](https://github.com/mattoverby/admm-elastic) includes
> self-collision and other improvements that are not present here.

## Contents

- [Algorithm overview](#algorithm-overview)
- [Paper-to-code map](#paper-to-code-map)
- [How to read the code](#how-to-read-the-code)
- [Forces and material models](#forces-and-material-models)
- [Repository layout](#repository-layout)
- [Build on Windows](#build-on-windows)
- [Build on Ubuntu or Debian](#build-on-ubuntu-or-debian)
- [Run the examples](#run-the-examples)
- [Controls](#controls)
- [Scene and solver configuration](#scene-and-solver-configuration)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)

## Algorithm overview

### 1. Backward Euler as minimization

Let \(x^n\) and \(v^n\) be the current positions and velocities, \(M\) the
diagonal mass matrix, and \(\Delta t\) the timestep. After explicit forces such
as gravity and wind update the velocity, the unconstrained prediction is

$$
\tilde{x} = x^n + \Delta t\,v^n.
$$

Implicit Euler can then be written as the minimization

$$
x^{n+1} = \underset{x}{\operatorname{argmin}}\;
\frac{1}{2\Delta t^2}\lVert M^{1/2}(x-\tilde{x})\rVert^2
+ \sum_i U_i(D_i x).
$$

Each \(D_i\) extracts a small element-local quantity from the global position
vector. For example, it can produce a spring edge, a triangle deformation
gradient, or a tetrahedron deformation gradient. \(U_i\) is the corresponding
elastic energy or constraint.

In the code, `System::step()` forms `x_bar` (the paper's \(\tilde{x}\)) and
`M_xbar`, while `m_D` stores all \(D_i\) matrices stacked by rows.

### 2. ADMM variable splitting

The method introduces local variables

$$
z_i = D_i x
$$

and scaled dual variables \(u_i\). Stacking all local quantities gives
\(z=Dx\) and \(u\). A diagonal weight matrix \(W\), assembled from per-force
weights \(w_i\), controls the augmented-Lagrangian penalty and convergence.

The main state has a direct representation in `System`:

| Mathematical quantity | Code |
| --- | --- |
| \(x\), global positions | `System::m_x` and local `curr_x` |
| \(v\), global velocities | `System::m_v` |
| \(M\), diagonal mass matrix | `System::m_masses` |
| \(D\), stacked reduction matrix | `System::m_D` |
| \(W\), diagonal ADMM weights | `System::m_W_diag` |
| \(z\), local primal variables | `System::curr_z` |
| \(u\), scaled dual variables | `System::curr_u` |

### 3. Local step

For fixed \(x\), every energy term independently solves a small proximal
problem of the form

$$
z_i \leftarrow \underset{z}{\operatorname{argmin}}\;
\Delta t^2 U_i(z) + \frac{w_i^2}{2}
\lVert D_i x-z+u_i\rVert^2,
$$

followed by its dual update. This work is parallel over forces:

```cpp
Dx = m_D * curr_x;

#pragma omp parallel for
for (int i = 0; i < n_forces; ++i) {
    forces[i]->project(dt, Dx, curr_u, curr_z);
}
```

`Force::project()` is therefore the key extension point. Simple projective
dynamics energies use an analytic projection. Neo-Hookean, StVK, and Fung
materials reduce the deformation gradient with an SVD and solve a small
singular-value proximal problem using L-BFGS.

### 4. Global step

For fixed \(z\) and \(u\), the position update is the sparse linear solve

$$
\left(M+\Delta t^2 D^T W^T W D\right)x =
M\tilde{x}+\Delta t^2D^TW^TW(z-u).
$$

The left-hand matrix stays constant while the topology, timestep, masses, and
weights remain unchanged. `System::initialize()` factorizes it once with
Eigen's `SimplicialLDLT`; each ADMM iteration in `System::step()` only rebuilds
the right-hand side and calls `solver.solve()`.

This prefactored global solve, together with parallel local projections, is the
main reason the method retains the speed and structure of projective dynamics.
For projective-dynamics energies, choosing \(w_i=\sqrt{k_i}\) gives nearly the
same iteration; for affine constraint manifolds the two methods are identical.

### 5. Finish the timestep

After the configured number of ADMM iterations, the code updates

$$
v^{n+1}=\frac{x^{n+1}-x^n}{\Delta t}, \qquad x^n\leftarrow x^{n+1}.
$$

`SimContext::update()` then copies the solver positions back into the render
meshes and refreshes the scene objects.

## Paper-to-code map

| Paper concept | Implementation |
| --- | --- |
| Implicit objective and Algorithm 1 | `deps/admm-elastic-sca/src/system/System.cpp` |
| Global variables, matrices, and factorization | `deps/admm-elastic-sca/src/system/System.hpp` |
| \(D_i\), \(w_i\), \(U_i\), and local proximal interface | `Force::get_selector()` and `Force::project()` |
| Parallel local step | OpenMP loop in `System::step()` |
| Prefactored global step | `System::initialize()` and `solver.solve()` |
| Triangle deformation energies | `TriangleForce.cpp` |
| Tetrahedral deformation energies | `TetForce.cpp` |
| Neo-Hookean and StVK proximal solves | `HyperElasticTet`, `NHProx`, `StVKProx` |
| Scene/XML to solver assembly | `src/SimContext.cpp` and `src/ForceBuilder.cpp` |
| Visualization and input | `deps/mclscene/src/Application.cpp` |

## How to read the code

The shortest useful route through the project is the windy-flag example:

1. Start at `samples/windyflag/windyflag.cpp`. `main()` loads `cloth.xml`,
   initializes the simulation, creates the viewer, and registers a callback that
   changes the wind strength.
2. Read `samples/windyflag/cloth.xml`. The scene section defines the cloth mesh,
   pole, mass, material, and force references. The `admmelastic` section defines
   gravity, triangle strain, bending, timestep, and iteration count.
3. Follow `SimContext::load()`. It parses solver/force descriptions and asks
   `mclscene::SceneManager` to create the scene objects.
4. Follow `ForceBuilder::admm_build_object()`. It copies dynamic mesh vertices
   into `System::m_x`, assigns masses, and creates one or more `Force` objects
   per triangle, edge, hinge, or tetrahedron.
5. Read `System::initialize()`. It initializes forces, stacks their selector
   matrices into `m_D`, assembles `W` and `M`, and factorizes the global matrix.
6. Read one call to `System::step()`: explicit-force prediction, local ADMM
   projections, global sparse solve, and position/velocity update.
7. Return through `SimContext::update()` to `Application::display()`, where the
   updated meshes are rendered and keyboard/mouse events are handled.

After this path, inspect one analytic local model such as
`LimitedTriangleStrain::project()`, followed by the nonlinear
`HyperElasticTet::project()` path.

## Forces and material models

All implicit forces derive from `admm::Force`. Each force supplies rows of the
global selector matrix, its ADMM weight, and its local projection/proximal solve.

| Model | Main class | Purpose |
| --- | --- | --- |
| Edge spring | `Spring` | Preserves an edge's rest length with a quadratic energy. |
| Triangle strain | `LimitedTriangleStrain` | Limits the singular values of a triangle deformation gradient. |
| Triangle area | `TriArea` | Projects triangle deformation toward an area range. |
| Cloth bending | `BendForce` | Penalizes bending over adjacent triangle hinges. |
| Linear tetrahedral strain | `LinearTetStrain` | Projective-dynamics-style tetrahedral deformation. |
| Volume preservation | `TetVolume` | Restricts tetrahedron volume to a configured range. |
| Neo-Hookean tetrahedron | `HyperElasticTet` + `NHProx` | Nonlinear volumetric elasticity solved in singular-value space. |
| StVK tetrahedron | `HyperElasticTet` + `StVKProx` | Saint Venant–Kirchhoff elasticity solved with a local L-BFGS proximal step. |
| Fung triangle | `FungTriangle` + `FungProx` | Nonlinear triangle constitutive energy. |
| Static/moving anchor | `StaticAnchor`, `MovingAnchor` | Hard or controllable positional constraints. |
| Static obstacle collision | `CollisionForce` | Projects vertices out of configured collision shapes. |
| Gravity/constant acceleration | `ExplicitForce` | Updates velocity before the implicit solve. |
| Wind | `WindForce` | Applies face-normal-dependent forces to triangle meshes. |

`ExplicitForce` and `WindForce` are not ADMM energy terms. They modify velocity
before the local/global iterations.

## Repository layout

```text
.
├── CMakeLists.txt                 # Top-level libraries and four sample targets
├── src/
│   ├── SimContext.*              # XML/scene/solver bridge
│   └── ForceBuilder.*            # Builds forces and masses from meshes
├── samples/
│   ├── windyflag/                # Cloth strain, bending, gravity, wind
│   ├── bunnyexpand/              # StVK recovery from extreme deformation
│   ├── plinkopony/               # Linear tet model and obstacle collisions
│   └── poordillo/                # Neo-Hookean model and moving anchors
└── deps/
    ├── admm-elastic-sca/         # Core ADMM solver and force models
    └── mclscene/                 # Scene loading, meshes, OpenGL viewer
```

Eigen, cppoptlib, trimesh2, TetGen, SOIL2, and pugixml are vendored under
`deps/`. GLFW, GLEW, OpenGL, and a C/C++ compiler are system dependencies.

## Build on Windows

The recommended setup is an **MSYS2 UCRT64** shell so that the compiler and all
binary libraries use the same runtime.

### 1. Install dependencies

From an MSYS2 UCRT64 shell:

```bash
pacman -Syu
pacman -S --needed \
  mingw-w64-ucrt-x86_64-gcc \
  mingw-w64-ucrt-x86_64-cmake \
  mingw-w64-ucrt-x86_64-ninja \
  mingw-w64-ucrt-x86_64-glfw \
  mingw-w64-ucrt-x86_64-glew
```

Close and reopen the UCRT64 shell if `pacman -Syu` requests it, then rerun the
installation command.

### 2. Configure and build

Run these commands from the repository root:

```bash
cmake --fresh -S . -B build-ucrt -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build-ucrt --parallel 4
```

The current code has been successfully configured with CMake 4.4 and built with
GCC 15 and Ninja on 64-bit Windows.

### 3. Run

From the same UCRT64 shell:

```bash
./build-ucrt/samples/windyflag.exe
```

Using the same shell is important because its `PATH` contains the UCRT64 GLFW,
GLEW, GCC runtime, and OpenMP DLLs. From PowerShell, add the same runtime first:

```powershell
$env:PATH = "C:\tools\msys64\ucrt64\bin;$env:PATH"
.\build-ucrt\samples\windyflag.exe
```

Adjust `C:\tools\msys64` if MSYS2 is installed elsewhere.

## Build on Ubuntu or Debian

Install the compiler, CMake, Ninja, and rendering dependencies:

```bash
sudo apt update
sudo apt install build-essential cmake ninja-build \
  libglfw3-dev libglew-dev libgl1-mesa-dev libglu1-mesa-dev
```

Configure and build from the repository root:

```bash
cmake --fresh -S . -B build-ucrt -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build-ucrt --parallel "$(nproc)"
```

Run a sample:

```bash
./build-ucrt/samples/windyflag
```

The Linux instructions reflect the project's declared dependencies but have not
been verified in this repository update.

### Incremental rebuild

After editing source code, keep the existing configuration and run:

```bash
cmake --build build-ucrt --parallel 4
```

Use `cmake --fresh` again after changing compilers, generators, dependency
locations, or important CMake options.

## Run the examples

Windows binaries are generated under `build-ucrt/samples/`:

| Executable | Demonstration |
| --- | --- |
| `windyflag.exe` | Strain-limited cloth with bending, gravity, anchors, and adjustable wind. |
| `bunnyexpand.exe` | StVK tetrahedral bunny recovering from randomized or collapsed positions. |
| `plinkopony.exe` | Linear tetrahedral horse falling through a field of cylindrical obstacles. |
| `poordillo.exe` | Neo-Hookean armadillo with interactive moving hand/foot constraints. |

PowerShell example:

```powershell
.\build-ucrt\samples\windyflag.exe
```

Linux example:

```bash
./build-ucrt/samples/windyflag
```

## Controls

Click the viewer window first so that it receives keyboard input.

### Common controls

| Input | Action |
| --- | --- |
| `Space` | Start or pause continuous simulation. |
| `P` | Advance one simulation step while paused. |
| Hold left mouse button + move | Rotate the camera. |
| Mouse wheel | Zoom in or out. |
| `S` | Toggle saving every rendered frame to PNG. |
| `Esc` | Close the viewer. |

Screenshots are written to the top-level `build-ucrt/` directory as timestamped
`screenshot_*.png` files.

### Example-specific controls

| Example | Input | Action |
| --- | --- | --- |
| Windy flag | `W` | Toggle normal and high wind strength. |
| Poor dillo | `H` | Toggle the hand control-point anchors. |
| Poor dillo | `F` | Toggle the foot control-point anchors. |

The bunny and plinko examples only use the common viewer controls.

## Scene and solver configuration

Each example combines an `mclScene` section with an `admmelastic` section.
Objects refer to force definitions by name:

```xml
<Object name="cloth1" type="plane">
    <width value="30" />
    <length value="20" />
    <Mass value=".5" />
    <Force value="admmstyle" />
    <Force value="bend" />
</Object>

<admmelastic>
    <Force name="admmstyle" type="TriangleStrain">
        <limit value=".95 1.05" />
        <stiffness value="100" />
    </Force>

    <solver>
        <iterations value="30" />
        <timestep value="0.04" />
        <realtime value="0" />
        <verbose value="1" />
    </solver>
</admmelastic>
```

### Important parameters

| XML parameter | Meaning |
| --- | --- |
| `Mass` | Total object mass in kilograms; distributed over mesh vertices. |
| `density_weighted_mass` | Use area/volume-weighted masses (`1`) or uniform masses (`0`). |
| Object `Force` | Name of a force definition applied to every relevant mesh element. |
| `iterations` | Number of local/global ADMM iterations per simulation timestep. |
| `timestep` | Fixed simulation timestep \(\Delta t\), in seconds. |
| `realtime` | If true, take enough fixed steps to cover the current rendered-frame time. |
| `verbose` | Solver output level. |
| `stiffness` | Energy stiffness for springs, strain, bending, or volume forces. |
| `limit` | Minimum and maximum allowed singular values for triangle strain. |
| `range_min`, `range_max` | Allowed volume range for `volpres`. |
| `mu`, `lambda` | Lamé parameters for Neo-Hookean and StVK tetrahedra. |
| `max_iterations` | Maximum L-BFGS iterations for a nonlinear local proximal solve. |
| `weight_scale` | Multiplier for the default linear-tet ADMM weight. |
| `direction` | Acceleration vector for gravity or direction/magnitude for wind. |

Changing masses, the timestep, topology, or force weights requires rebuilding
the global matrix with `System::initialize()` or `System::recompute_weights()`.
Changing only an explicit force direction does not.

## Extending the solver

To add a new implicit energy:

1. Derive a class from `admm::Force`.
2. In `initialize()`, compute rest-state data and a suitable ADMM weight.
3. In `get_selector()`, append the rows of \(D_i\) and the corresponding
   diagonal entries of \(W\).
4. In `project()`, solve the local proximal problem and update the relevant
   slices of `z` and `u`.
5. Teach `ForceBuilder` how to construct the new force from XML.

The local problem should remain low-dimensional and element-local to preserve
the algorithm's parallelism. For isotropic tetrahedral energies, the paper
recommends reducing the proximal solve to the three singular values of the
deformation gradient.

## Troubleshooting

### CMake reports a version or policy error

Use CMake 3.10 or newer. The project declares compatibility through CMake 4.0
and has been tested with CMake 4.4:

```bash
cmake --version
cmake --fresh -S . -B build-ucrt -G Ninja
```

### GLFW, GLEW, or OpenGL is not found

Install the development packages for the active compiler environment. On
Windows, launch the UCRT64 shell and confirm that `cmake`, `g++`, GLFW, and GLEW
all come from the UCRT64 prefix. On Linux, install the packages listed above and
rerun the fresh configure command.

### A Windows executable reports a missing DLL

Run it from the MSYS2 UCRT64 shell, or put the UCRT64 `bin` directory on
`PATH`. Do not copy arbitrary DLLs from another MinGW distribution into the
build directory.

### CMake selects the wrong compiler or architecture

Do not mix MSVC, 32-bit MinGW, MinGW64 MSVCRT, and MinGW64 UCRT libraries in one
build. Open the desired shell, remove or fresh-configure the build cache, and
confirm the compiler recorded in `build-ucrt/CMakeCache.txt`.

### Build files are absent from `git status`

This is intentional. `/build/`, `/build-ucrt/`, and `/.cmake4-check/` are
ignored because they contain machine-specific caches, objects, DLL metadata,
and executables.

## Known limitations

- This is the SCA 2016 code, not the extended TVCG solver.
- Self-collision is not implemented here.
- The XML/scene layer was designed for the paper examples rather than as a
  stable public file format.
- The code uses fixed ADMM iteration counts by default; residual-based early
  stopping is mentioned in `System::step()` but not enabled.
- Some vendored libraries are old and may emit warnings with modern compilers.

## Citation

If this code or method contributes to academic work, cite the SCA 2016 paper:

```bibtex
@inproceedings{Narain2016,
  author    = {Narain, Rahul and Overby, Matthew and Brown, George E.},
  title     = {{ADMM} $\supseteq$ Projective Dynamics: Fast Simulation of General Constitutive Models},
  booktitle = {Proceedings of the ACM SIGGRAPH/Eurographics Symposium on Computer Animation},
  series    = {SCA '16},
  year      = {2016},
  isbn      = {978-3-905674-61-3},
  location  = {Zurich, Switzerland},
  pages     = {21--28},
  numpages  = {8},
  doi       = {10.2312/sca.20161219},
  publisher = {Eurographics Association}
}
```

Project page: <https://www.cse.iitd.ac.in/~narain/admm-pd/>

## License

Copyright (c) 2017 University of Minnesota.

This project is distributed under the BSD 2-Clause License. See
[`LICENSE.txt`](LICENSE.txt) for the complete license text.
