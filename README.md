# sqSTFT.jl — Short-Time Fourier Transform (STFT) & Synchrosqueezing (SST)

## Overview
This Julia implementation computes a **Short-Time Fourier Transform (STFT)** and its **synchrosqueezed** variant (SST), returning both the spectrogram and its reassigned, sharpened time–frequency representation. It is adapted from Hau-tieng Wu’s synchrosqueezing routine (originally in MATLAB) and translated/extended in Julia by Naoki Saito, Haotian Bao.

## Features
- STFT with arbitrary, **odd-length** analysis window `h` and its derivative `Dh`
- Frequency reassignment via **synchrosqueezing** to concentrate TF energy
- Configurable frequency resolution (`N`) and hop size (via `t`)
- Optional progress timing/tracing
- Example plotting snippet provided (commented)

## Function Signature
```julia
tfr, rtfr = sqSTFT(x, t, N, h, Dh; trace=false)
```

### Arguments
- `x::AbstractVector`  
  Real or complex signal (one column only).
- `t::AbstractVector{Int}`  
  Time indices (e.g., `1:dt:len`). Must be **regularly sampled**; `dt` controls hop size.
- `N::Int`  
  Number of frequency bins (ideally a power of two for speed; the code warns otherwise).
- `h::AbstractMatrix`  
  **Column** window with **odd** length `hrow`, `hcol == 1`. `h[0]` (center) conceptually forced to 1.
- `Dh::AbstractMatrix`  
  Derivative of `h` with the **same shape** as `h`.
- `trace::Bool` (default `false`)  
  Print progress and timing info.

### Returns
- `tfr::Matrix{ComplexF64}`  
  STFT (size `N × length(t)`), computed via FFT over frequency.
- `rtfr::Matrix{ComplexF64}`  
  Synchrosqueezed TF representation (size `floor(Int,N/2) × length(t)`).

## Requirements & Installation
The script installs dependencies itself, but you can pre-add them:

```julia
using Pkg
Pkg.add.(["TickTock", "JLD", "SpecialFunctions", "Plots", "FFTW", "Statistics", "LinearAlgebra"])
```

**Julia:** 1.6+ recommended  
**Packages:** `FFTW`, `LinearAlgebra`, `Statistics`, `SpecialFunctions`, `Plots`, `TickTock`, `JLD`

Also ensure the helper files are available in the same directory:
- `hermf.jl`  — Hermite functions window generator
- `nround.jl` — Custom rounding utility

Load them via:
```julia
include("hermf.jl")
include("nround.jl")
```

## How It Works (High Level)
1. **Windowed Framing**  
   For each time index `tᵢ`, the routine extracts the local segment `x[tᵢ+τ]` within window half-length `Lh`, forms two columns:
   - STFT numerator: `x .* conj(h)`
   - Derivative track: `x .* conj(Dh)`  
   Both are **normalized by ‖h‖** over the active support to keep amplitudes consistent.

2. **Frequency Transform**  
   Apply `fft(., 1)` across rows to obtain `tfr` (STFT) and `tf3` (needed for instantaneous frequency estimation).

3. **Instantaneous Frequency Estimation**  
   Where `tfr` is non-zero and above a **threshold** `1e-8 * mean(abs2, segment)`, compute the frequency shift index
   \[
   \hat{\jmath} = \mathrm{round}\Big(\Im\big(N\; tf3 / (2\pi\, tfr)\big)\Big)
   \]
   (ties rounded up) and wrap into valid bin range.

4. **Synchrosqueezing (Reassignment)**  
   Accumulate STFT energy from bin `j` into reassigned bin `\(\hat{\jmath}\)` to build `rtfr` (only up to `N/2` retained for positive frequencies).

This sharpened `rtfr` yields ridges with improved readability for AM–FM signals compared to a standard spectrogram.

## Parameter Tips
- **N (frequency bins):** Use a power of two (e.g., 1024, 2048) for FFT speed.
- **Window length (h):** Longer windows improve frequency resolution but smear time; ensure **odd length**.
- **Hop size (t):** Smaller step (denser `t`) improves time resolution at higher compute cost.
- **Threshold:** The code auto-sets a tiny fraction of mean energy; adjust if needed for very low or high SNR.

## Common Pitfalls & Troubleshooting
- **“X must have only one column!”**  
  Supply a **vector** signal (not a 2-D array).
- **“T must only have one row!” / irregular `t`**  
  Provide a 1-D index vector with **constant spacing** (e.g., `1:dt:len`).
- **“H must be … odd length”**  
  Ensure `size(h) == (odd, 1)` and likewise for `Dh`.
- **Artifacts / empty `rtfr`:**  
  Check normalization, verify `H`/`DH` pair matches (both centered, same length), and consider larger `N` or different window lengths.

## Performance Notes
- The routine warns if `N` isn’t a power of two; prefer powers of two to leverage FFTW performance.
- Avoid unnecessary allocations by reusing arrays if you adapt the code for batch processing.
- `rtfr` returns only **positive frequencies** (`N/2` rows). If you need full spectrum, you can adapt the accumulation loop.

## Acknowledgments
- Original synchrosqueezing concept and MATLAB routine: Hau-tieng Wu / F. Auger et al.  
- Julia translation and modifications: Naoki Saito, Haotian Bao.

## License
The header notes legacy confidentiality from the original code comment block. If you plan to distribute or publish derivatives, verify your rights and the original authors’ licensing terms. For academic use, please credit the authors accordingly.
