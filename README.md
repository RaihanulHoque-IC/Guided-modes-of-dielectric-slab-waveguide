# Guided Modes of Dielectric Slab Waveguide with Anisotropic Core

**Master's Project** | Universität Paderborn — Fachgebiet Theoretische Elektrotechnik

**Submitted to:** Prof. Dr. Jens Förstner  
**Supervisors:** Dr. Manfred Hammer, Henna Farheen  
**Presented by:** Shaiera Washima, Rajdeep Roy, Muhammed Raihanul Hoque  
**Date:** 3rd July 2024

---

## Table of Contents
- [Overview](#overview)
- [Waveguide Structure](#waveguide-structure)
- [Theory](#theory)
  - [Isotropic Core](#isotropic-core)
  - [Anisotropic Core (X-cut)](#anisotropic-core-x-cut)
  - [Anisotropic Core (Z-cut)](#anisotropic-core-z-cut)
- [Numerical Method](#numerical-method)
- [Results](#results)
- [Concluding Remarks](#concluding-remarks)
- [References](#references)

---

## Overview

This project implements and analyzes **guided modes** of a three-layer dielectric slab waveguide with both **isotropic** and **anisotropic** (Lithium Niobate, LiNbO₃) core media. The analytic procedures for modal analysis are outlined and implemented in **MATLAB**, with results validated against the OMS and TFLN online solvers.

Key analyses performed:
- Dependence of effective indices on waveguide thickness
- Dependence of effective indices on oblique propagation angle (θ)
- Characterization of mode polarization
- Mode hybridization in anisotropic media

---

## Waveguide Structure

```
        ↑ x
        |          n₁ (cladding)
  a ────|──────────────────────────
        |    a      ε̂ (core)
  0 ────|──────────────────────────────→ z
        y          n₃ (substrate)
```

- **Guiding principle:** Total internal reflection, requiring `nf > ns, nc`
- **Non-zero evanescent field** exists in the cladding and substrate
- The structure is **homogeneous along y and z**: `∂y ε = 0`, `∂z ε = 0`
- Fields are constant along y and vary harmonically with z:

```
∂y E = 0,  ∂y H = 0,  ∂z = −jβ

(E, H)(x, z) = (Ē, H̄)(x) · e^(−jβz)
```

Where:
- `β` — Propagation constant
- `Neff = β / k₀` — Effective index
- `(Ē, H̄)` — Mode profile

---

## Theory

### Core Permittivity Tensors

| Type | Tensor |
|------|--------|
| **Isotropic** | `diag(n², n², n²)` |
| **Anisotropic X-cut** | `diag(no², ne², no²)` |
| **Anisotropic Z-cut** | `diag(ne², no², no²)` |

For **oblique propagation** with rotation angle θ around the x-axis, the X-cut permittivity tensor becomes:

```
ε̂_f^x = R(θ) · ε̂_f,0^x · R^T(θ)

      ⎡ εx   0    0  ⎤
    = ⎢  0   εy   δ  ⎥
      ⎣  0   δ    εz ⎦
```

where `δ = (ne² − no²) sin θ cos θ` is the off-diagonal coupling term responsible for mode hybridization.

For **Z-cut**, the rotation leaves the tensor unchanged: `ε̂_f^z = ε̂_f,0^z`.

---

### Isotropic Core

The wave equations decouple into pure **TE** and **TM** modes:

**TE mode** (principal component Ēy):
```
∂²Ēy/∂x² + (ω²ε − β²)Ēy = 0
H̄z = (−1/jωμ₀)(∂Ēy/∂x),   H̄x = −(β/ωμ₀)Ēy
```

**TM mode** (principal component H̄y):
```
∂²H̄y/∂x² + (ω²ε − β²)H̄y = 0
Ēz = (1/jωε)(∂H̄y/∂x),     Ēx = (β/ωε)H̄y
```

**Field solutions (piecewise):**
```
Ēy(x) = { E₁ e^(−k₁x),              x ≥ a
         { E₂cos(k₂x) + E₃sin(k₂x), 0 ≤ x ≤ a
         { E₄ e^(k₃x),               x ≤ 0
```

Guided modes exist only when:
```
n₁ ≤ n₃ ≤ β/k₀ ≤ n₂
```

---

### Anisotropic Core (X-cut)

The off-diagonal term `δ ≠ 0` **couples** the TE and TM equations:

```
∂²Ey/∂x² − jωμ(δ/εz)∂x Hy + [ω²με₀(εy − δ²/εz) − β²]Ey = 0

∂²Hy/∂x² − jωε₀δ ∂x Ey + [ω²με₀εz − β²(εz/εx)]Hy = 0
```

**Matrix form:**
```
A ∂x² ψ + B ∂x ψ + C ψ = 0

where ψ = [Ey, Hy]^T
```

**Ansatz for the anisotropic core:**
```
ψ(x) = Σ (m=1 to 4) Am · ψm · e^(−jkm x)
```

The `km` and `ψm` are obtained by solving the eigenvalue problem `MV = λV`.

**Polarization mixing parameter Π:**
```
Π = (−Re ∫ Ey* Hx dx) / (−Re ∫ (Ex* Hy − Ey* Hx) dx)
```
- `Π = 1` → Pure TE
- `Π = 0` → Pure TM
- `0 < Π < 1` → Hybrid mode

---

### Anisotropic Core (Z-cut)

With `δ = 0`, the wave equations **decouple** (no hybridization):

```
∂²Ey/∂x² + (ω²με₀εy − β²)Ey = 0

∂²Hy/∂x² + (ω²με₀εz − β²(εz/εx))Hy = 0
```

---

## Numerical Method

### Isotropic Case
1. Apply boundary conditions (continuity of `Ēy, H̄z` for TE; `H̄y, Ēz` for TM) at `x = 0` and `x = a`
2. Assemble 4×4 matrix `M(β)`
3. Find roots of `det(M(β)) = 0` using the **Bisection method**
4. Reconstruct field profiles from the corresponding null eigenvectors

### Anisotropic Case
1. Compute core wavenumbers `km` and eigenvectors `ψm` from the 4×4 eigenvalue problem
2. Apply boundary conditions → 8 linear equations in 8 unknowns
3. Find `β` where the **product of eigenvalues** of `M(β)` equals zero (Minima-finding method)
4. Extract field amplitudes from the null eigenvector of `M(β)`

---

## Results

All simulations use vacuum wavelength **λ = 1.55 μm** with Lithium Niobate parameters:
`no = 2.1837`, `ne = 2.1227`

### Isotropic — Neff vs Thickness

| Config | n₁ | n₂ | n₃ |
|--------|-----|-----|-----|
| Asymmetric | 1.0 | 3.45 | 1.45 |
| Symmetric | 1.45 | 3.45 | 1.45 |

- TE modes (blue) have higher Neff than TM modes (black)
- Neff increases monotonically with waveguide half-thickness `a`
- Example at `a = 0.22 μm`: TE₀ Neff = 2.8051, TM₀ Neff = 1.8747
- Example at `a = 0.44 μm`: TE₀ Neff = 3.1973, TM₀ Neff = 3.0180

---

### Anisotropic X-cut — Neff vs Thickness (θ = 45°)

**Asymmetric (n₁ = 1, n₃ = 1.4483), a = 0.45 μm:**
| Mode | Π | Neff |
|------|---|------|
| (a) | > 0.99 (near pure TE) | 1.8919 |
| (b) | < 0.01 (near pure TM) | 1.7337 |

**Symmetric (n₁ = n₃ = 1.4483), a = 0.45 μm:**
| Mode | Π | Neff |
|------|---|------|
| (a) | > 0.99 | 1.9122 |
| (b) | < 0.01 | 1.8073 |

**Mode hybridization region** (asymmetric, a ≈ 1.03 μm):
| Mode | Π | Neff |
|------|---|------|
| (a) | 0.36 | 2.0713 |
| (b) | 0.64 | 2.0705 |

> In the symmetric case (n₁ = n₃), modes remain near-purely polarized (Π ≈ 0 or 1) across all thicknesses — no mode crossing/hybridization occurs.

---

### Anisotropic X-cut — Neff vs θ

- For small `a` (0.45 μm): Neff varies smoothly with θ; modes remain near-purely polarized
- For large `a` near the hybridization thickness (1.03 μm): strong polarization exchange occurs near θ = ±45°
- In the symmetric waveguide, hybridization is absent even at larger thicknesses

---

### Anisotropic Z-cut — Neff vs Thickness

- No mode hybridization (δ = 0 → decoupled equations)
- Clear TE/TM separation maintained across all thicknesses
- Example at `a = 0.45 μm` (asymmetric):
  - TE₀: Neff = 1.9219
  - TM₀: Neff = 1.7029
- Example at `a = 0.45 μm` (symmetric):
  - TE₀: Neff = 1.9412
  - TM₀: Neff = 1.773

---

## Concluding Remarks

- Analytic modal analysis procedures for a **three-layer dielectric slab waveguide** with isotropic and anisotropic (Lithium Niobate) core have been outlined and implemented in **MATLAB**
- **X-cut LiNbO₃** introduces off-diagonal permittivity terms that couple TE and TM equations, leading to **hybrid modes** and **mode hybridization** near thickness crossover points
- **Z-cut LiNbO₃** preserves a diagonal permittivity tensor under oblique propagation, resulting in **no hybridization**
- Results have been validated against the **OMS** and **TFLN** online solvers:
  - https://www.siio.eu/oms.html
  - https://www.siio.eu/tflns.html

---

## References

1. Lena Ebers et al. "Flexible source of correlated photons based on LNOI rib waveguides". *Journal of Physics: Photonics* 4.2 (2022), p. 025001.
2. M. Hammer. "1-D mode solver for dielectric multilayer slab waveguides". *Software online from SiIo.Eu* (2022). https://www.siio.eu/oms.html
3. Press WH, Vetterling WT, Teukolsky SA, Flannery BP. *Numerical Recipes*. Cambridge University Press, London, England; 1988.
4. Lecture "Fields and Waves" — Prof. Dr. Jens Förstner, TET, University of Paderborn
5. Lecture "Optical Waveguide Theory" — Dr. Manfred Hammer, University of Paderborn
