# ForwardRCS — forward real Compton scattering off the proton

Empirical evaluation of the forward Compton scattering amplitudes of the proton from photoabsorption cross
sections by dispersion relations, with the forward sum rules (Baldin, GDH and higher order) and the forward
observables.

**Code:** Oleksii Gryniuk (main author), Vladimir Pascalutsa.
The code was written for the two papers below; please cite them when you use it (see also `CITATION.cff`):

1. O. Gryniuk, F. Hagelstein, V. Pascalutsa, *Evaluation of the forward Compton scattering off protons:
   Spin-independent amplitude*, Phys. Rev. D **92**, 074031 (2015),
   [doi:10.1103/PhysRevD.92.074031](https://doi.org/10.1103/PhysRevD.92.074031),
   [arXiv:1508.07952](https://arxiv.org/abs/1508.07952).
2. O. Gryniuk, F. Hagelstein, V. Pascalutsa, *Evaluation of the forward Compton scattering off protons:
   II. Spin-dependent amplitude and observables*, Phys. Rev. D **94**, 034043 (2016),
   [doi:10.1103/PhysRevD.94.034043](https://doi.org/10.1103/PhysRevD.94.034043),
   [arXiv:1604.00789](https://arxiv.org/abs/1604.00789).

## What is computed

The two forward proton Compton-scattering amplitudes, $f(\nu)$ and $g(\nu)$,
with $\nu$ the photon lab energy, are determined from exact once-subtracted dispersion relations:

$$f(\nu) = -\frac{\alpha}{M} + \frac{\nu^2}{2\pi^2} \int_{\nu_0}^\infty \frac{\sigma(\nu') d\nu'}{\nu'^2-\nu^2},\qquad
g(\nu) = -\frac{\alpha\kappa^2}{2M^2}\nu - \frac{\nu^3}{4\pi^2} \int_{\nu_0}^\infty \frac{\Delta\sigma(\nu') d\nu'}{\nu'(\nu'^2-\nu^2)},$$

where $\alpha$ is the fine-structure constant, $M$ and $\kappa$ are the proton mass and anomalous magnetic moment, $\sigma = (\sigma_{3/2}+\sigma_{1/2})/2$ is the total
photoabsorption cross section and $\Delta\sigma = \sigma_{3/2}-\sigma_{1/2}$ the helicity difference.  The following sum rules are evaluated too:

| | |
|---|---|
| Baldin | $\alpha_{E1}+\beta_{M1} = \frac{1}{2\pi^2}\int \sigma/\nu^2 d\nu$ |
| 4th, 6th order | $\frac{1}{2\pi^2}\int \sigma/\nu^4 d\nu$, $\frac{1}{2\pi^2}\int \sigma/\nu^6 d\nu$ |
| GDH | $\int \Delta\sigma/\nu d\nu$ compared with $2\pi^2\alpha\kappa^2/M^2 = 204.78\ldots$ μb |
| forward spin polarizabilities | $\gamma_0 = -\frac{1}{4\pi^2}\int \Delta\sigma/\nu^3 d\nu$, $\overline\gamma_0 = -\frac{1}{4\pi^2}\int \Delta\sigma/\nu^5 d\nu$ |

Forward observables: $d\sigma/d\Omega(\theta=0) = |f|^2+|g|^2$ and
$\Sigma_{2z} = -2 \mathrm{Re}(f g^*)/(|f|^2+|g|^2)$.

The cross section is represented in three regions: a partial-wave table (MAID, SAID) from the pion threshold
$\nu_0$ up to the $2\pi$ threshold $E_1 = 0.309$ GeV, a fit (Breit–Wigner resonances plus background) up to
2 GeV, and a Regge form above. Errors are propagated linearly through the covariance of the fit parameters.

**Units and conventions of the code.** Photon lab energy $\nu$ in GeV; cross sections in μb;
amplitudes in μb·GeV (divide $|f|^2$ by $10(\hbar c)^2$ for nb/sr); sum rules in $10^{-4} \mathrm{fm}^{3,5,7}$
(Baldin, 4th, 6th), μb (GDH), $10^{-6} \mathrm{fm}^4$ ($\gamma_0$) and $10^{-6} \mathrm{fm}^6$ ($\overline\gamma_0$).
The polarized data column is $\Delta\sigma/2 = (\sigma_{3/2}-\sigma_{1/2})/2$.

## Quick start

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python tot.py                              # fit I: sum rules by region (print_all)
.venv/bin/python pol.py                              # polarized: GDH, gamma0, gamma0bar
.venv/bin/python -m unittest discover -s tests       # regression test, 2-3 min
.venv/bin/python unpol2026.py; .venv/bin/python pol2026.py; .venv/bin/python amp2026.py   # 2026 update, ~5 min each
```

Run everything from the repository root (the data are read with relative paths).

## Layout

| | |
|---|---|
| `compton.py` | constants, kinematics, thresholds $E_0$ (π⁰p), $E_1$ (2π⁰), $E_2$ = 2 GeV |
| `val.py`, `common.py` | value-with-error class, I/O, interpolation |
| `fit.py` | fit models: Breit–Wigner sum + background, Regge form, smooth junction; tree-level pion-production cross sections (`BornT`, `BornP`) |
| `integrate.py` | regions, principal-value dispersion integrals, sum rules, error propagation |
| `tot.py`, `pol.py`, `amp.py` | the 2015/2016 evaluations: unpolarized fits I and II and the Ahrens fit; polarized fit; forward observables |
| `suppl.py`, `plot*.py`, `gen_figs.sh` | tables of the APS Supplemental Material, figures (from `out/`) |
| `data/tot`, `data/pol` | the 2015/2016 databases (fit I: `mdata.dat`, fit II: `sdata.dat`, polarized: `pol/data.dat`) |
| `out/` | outputs of the 2016-era code, kept for reference |
| `tests/test_regression.py` | regression test (pinned numbers, published tables, audit fixes) |
| `unpol2026.py`, `pol2026.py`, `amp2026.py`, `data/2026`, `results/2026` | the 2026 update (below) |

`integral_c.py`, `pion.py`, `gdh_data.py` and `perr.py` are auxiliary scripts of the original code
(`perr.py` needs a module `custom.py` that is not part of the frozen code).

## History

* 2015–2016: Python 2 code by O. Gryniuk for papers 1 and 2. The first commit of this repository is the code as
  frozen on 7 May 2020, unchanged.
* 2026: Python 3 port (behaviour-preserving), an audit with fixes, and an update of the evaluation with new data and
  global fits. Every change is a separate commit with its checks in the message.

### Fixes to the 2015/2016 evaluation

| commit | fix | effect on published numbers |
|---|---|---|
| a3cb856 | precomputed threshold integrals: no node above $E_1$; linear (not cubic) interpolation, as in the runs behind the published tables | none on the published grid; the frozen code was wrong by up to 7 μb·GeV within 10 MeV below $E_1$ |
| d1b58b3 | crash when $\nu$ equals a table energy | none |
| cf4ad76 | $\Sigma_{2z}$ errors (relative error stored as absolute) | none (the papers use $-2\,\mathrm{Re}\,fg^*$) |
| 6f84e53, eafc153, f49f08b | error matrices from the fitted points; MAID/SAID anchor points of fits I, II with their 2 % error; $E_1$ belongs to the PWA table | α+β: 14.29 → 14.293 (I), 13.85 → 13.848 (II) |
| aa0da58 | converged quadrature above $E_1$; quadrature error no longer added to the errors | ≤ 3×10⁻³ in Re f (near the 1.70 GeV resonance) |

Current numbers of the 2015/2016 evaluation (all within the published errors; the published Baldin errors 0.27 and
0.22 were a χ²-weighted misfit measure summed over regions, the code propagates the covariance):

| | fit I | fit II | published (I / II) |
|---|---|---|---|
| α+β [10⁻⁴ fm³] | 14.293(127) | 13.848(126) | 14.29(27) / 13.85(22) |
| 4th order [10⁻⁴ fm⁵] | 6.082(91) | 6.008(91) | 6.08(12) / 6.01(11) |
| 6th order [10⁻⁴ fm⁷] | 4.359(78) | 4.421(80) | 4.36(9) / 4.42(8) |

| | polarized | published |
|---|---|---|
| $I_\mathrm{GDH}$ [μb] | 204.46(21.37) | 204.5(21.4) |
| $\gamma_0$ [10⁻⁶ fm⁴] | −92.88(10.39) | −92.9(10.5) |
| $\overline\gamma_0$ [10⁻⁶ fm⁶] | 48.43(8.23) | 48.4(8.2) |

## 2026 update

Global fits replace the two unpolarized fits of 2015 and the polarized fit of 2016; the code above stays for
reproducing the papers. The update is described in, and used by, V. Pascalutsa, *Chiral-EFT prediction for proton
polarizabilities tested by the new Compton-scattering data* (2026, in preparation; its code and data:
[HiWARPs/ComptonBChPT](https://github.com/HiWARPs/ComptonBChPT)); please cite it together with papers 1 and 2 when
you use the 2026 results. The update follows an audit of the 2016 evaluation (normalization errors ignored,
inconsistent error matrices, the fit I/fit II difference) and uses the data and analyses that appeared since.

**Unpolarized (`unpol2026.py`).** One fit to all 322 data points above the 2π threshold (Armstrong 72, LEGS, MacCormick 96,
Bartalini 08, Meyer 70, Michalowski 77, Caldwell 73/78, H1 95, ZEUS 02) with one normalization nuisance per data set
(multiplying the model), the GHP model (six Breit–Wigner terms with background up to 2 GeV, Regge form above), all
parameters fitted together. Below $E_1$: the average of MAID2007 and SAID SM22, with error 2 % ⊕ half their
difference, as one correlated mode shared with the continuity point at $E_1$ (the data pull it down by
0.9σ = 1.9 %). χ²/dof = 448.8/295; the covariance is scaled by $S^2$, $S = 1.23$. The fit I/fit II difference of 2015
(0.44 in α+β) came mostly from Armstrong vs MacCormick on the high-energy side of the Δ; the global fit uses both,
with fitted normalizations.

**Polarized (`pol2026.py`).** Data: A2 2026 (0.2–1.4 GeV, 222 points), MAMI 2001 (all points from 0.2 GeV), ELSA
2003/2004, each with a scale nuisance. Below 0.2 GeV: the MAID2007/SM22 average as one error mode. From 0.2 GeV to
$W$ = 2.5 GeV: relativistic Breit–Wigner terms with a smooth background; above: the Regge fit of Bass, Skurzok and
Moskal (2018), chosen for its smaller uncertainty (the Regge parametrization of paper 2, Helbing (2006), gives an
$I_\mathrm{GDH}$ lower by 2.0 μb with an error of ±21.4 μb; variant in `results/2026/pol_variants.json`). χ²/dof = 209.9/266. The second error is the spread of the partial-wave analyses below 0.2 GeV
(MAID2007, SM22, CM12), quoted separately.

Results: value ± (exp+fit) ± (systematic). The systematic is the data selection (below) for the unpolarized and the
partial-wave model below 0.2 GeV for the polarized quantities. The unpolarized sum rules have in addition a one-sided,
upward systematic (in parentheses): the elastic channel γp → γp, which no hadronic cross section contains (below).

| | 2026 | GHP 2015/2016 |
|---|---|---|
| α+β [10⁻⁴ fm³] | **14.01 ± 0.18 ± 0.10 (+0.07)** | 14.0(2) (average of fits I and II) |
| 4th order [10⁻⁴ fm⁵] | 5.95 ± 0.09 ± 0.05 (+0.03) | 6.08(12) / 6.01(11) |
| 6th order [10⁻⁴ fm⁷] | 4.29 ± 0.08 ± 0.04 (+0.02) | 4.36(9) / 4.42(8) |
| $I_\mathrm{GDH}$ [μb] | **209.8 ± 9.2 ± 1.9** | 204.5(21.4) |
| $\gamma_0$ [10⁻⁶ fm⁴] | **−102.5 ± 8.0 ± 6.2** | −92.9(10.5) |
| $\overline\gamma_0$ [10⁻⁶ fm⁶] | 44.8 ± 5.8 ± 8.3 | 48.4(8.2) |
| $\gamma_0$, GDH sum rule imposed | −98.6 ± 3.6 ± 5.0 | |
| $\overline\gamma_0$, GDH sum rule imposed | 46.6 ± 4.8 ± 7.9 | |

**Data selection.** The 2015 fits I and II differ in the data used on the high-energy side of the Δ
(0.31–0.47 GeV: Armstrong 72 and LEGS in fit I, MacCormick 96 in fit II). Refitting the 2026 model to the two
selections gives α+β = 14.30(21) and 14.10(21): with fitted normalizations the 0.44 gap shrinks to 0.20. Locally the
shape difference remains: at 0.40 GeV the two selections give σ = 294 and 272 μb, and Re f differs by up to
0.6 μb·GeV at 0.45–0.5 GeV. The global fit follows MacCormick there. Half the difference of the two refits is quoted
as a separate data-selection systematic: 0.10 in α+β, 0.05 and 0.04 in the 4th- and 6th-order sum rules, up to
0.32 μb·GeV in Re f (0.45–0.5 GeV; below 0.3 GeV at most 0.03), and correspondingly for the forward observables
(`results/2026/unpol_selection.txt`, `results/2026/obs_selection.txt`). An analysis with independent data
(exclusive channels on the Δ tail) is under way.

**Channels absent from the database (one-sided).** The sum rules hold for the total photoabsorption cross section,
and a channel left out can only raise them, since every cross section is positive. At the order in α at which the
polarizabilities are defined, no hadronic channel is missing: below the 2π threshold the PWA table contains the only
ones, π⁰p and π⁺n, and above it the data are total photoabsorption measurements. Absent from all of them are the
electromagnetic final states, dominated by elastic scattering γp → γp (pair production is pure QED and is excluded by
the definition of the polarizabilities). Its cross section, from the parameter-free NNLO BχPT prediction of the 2026
paper (`data/2026/unpol_sigma_el.dat`; 3.2 μb at the Δ peak, 0.6 % of σ, and within 15 % of the estimate
4π · 0.8 · (M²/s) dσ/dΩ(0) from the forward amplitude of this code), adds 0.047 to α+β from the pion threshold to
0.45 GeV (0.034 of it beyond the Born and π⁰-pole terms) and at most 0.022 above, where σ_el ≤ r σ with r = 0.54 %,
its value at 0.45 GeV (it falls to about 0.1 % at 10 GeV). The sum, +0.07 (4th and 6th order: +0.03, +0.02), is quoted
as a one-sided systematic; the central values stay hadronic. It is not added to them because it is of order α, like
the radiative corrections to the hadronic cross sections, which have either sign. Moreover, below the pion threshold
the Thomson cross section makes its contribution to ∫σ/ν² infrared divergent: at this order the polarizabilities can
only be defined together with the radiative corrections. The shift of f (dispersive, σ_el from 0, where the integral
for Re f converges) is in `results/2026/unpol_missing.txt`, that of the observables in `results/2026/obs_missing.txt`
(+5.8 nb/sr, or 0.9 %, in dσ/dΩ(0) at the Δ peak). It is also in the forward points, as `miss_f` in
`results/2026/forward_points.json`: +0.0007 to +0.009 μb·GeV at 50–150 MeV, about a quarter of the exp+fit error.

Error budget of α+β: PWA mode 0.140, normalizations 0.043, data 0.011, scale factor 0.106 (in quadrature). Model
choices not included in the error move α+β by less than 0.1 (normalization widths 2 %/5 %: −0.10/+0.06; junction at
1.6 GeV: +0.08; data below $E_1$ added: −0.06; see `results/2026/unpol_variants.txt`). For comparison, A2 2026
quote $I_\mathrm{GDH}$ = 210 ± 2 ± 18 μb.

Amplitudes [μb·GeV] (exp+fit errors in parentheses; the data-selection systematic and the one-sided elastic-channel
shift of Re f in the next rows; the PWA-model systematic of Re g is in `results/2026/pol_g_pwa_model.txt`):

| ν [GeV] | 0.05 | 0.10 | 0.15 | 0.20 | 0.30 | 0.50 | 1.0 | 2.2 |
|---|---|---|---|---|---|---|---|---|
| Re f | −2.846(2) | −2.222(11) | −0.444(40) | 1.025(45) | 2.033(67) | −4.143(54) | −8.97(21) | −10.66(29) |
| ± data selection | 0.001 | 0.006 | 0.022 | 0.027 | 0.034 | 0.31 | 0.10 | 0.09 |
| elastic channel (one-sided) | +0.001 | +0.003 | +0.009 | +0.019 | +0.036 | −0.005 | −0.03 | −0.04 |
| Re g | −0.2625(3) | −0.538(3) | −0.525(24) | −2.414(40) | −2.978(35) | 0.32(13) | 1.46(31) | 3.16(60) |

Tables in `results/2026`: `unpol_f.txt` (Re f, Im f), `unpol_selection.txt` (its data-selection systematic),
`unpol_missing.txt` (its elastic-channel shift),
`pol_g.txt` (Re g, Im g), `pol_g_pwa_model.txt`, `obs.txt` ($d\sigma/d\Omega(0)$, $-2\,\mathrm{Re}\,fg^*$, $\Sigma_{2z}$, by
`amp2026.py`), `obs_selection.txt`, `obs_missing.txt`, the sum rules with their error budgets
(`*_sumrules.json`) and the model variants (`*_variants.*`). The databases with the source of every point and the
normalizations are `data/2026/db_unpolarized.csv` and `data/2026/pol_db.csv`.

**Inputs.**
- Unpolarized data: T. A. Armstrong *et al.*, Phys. Rev. D 5, 1640 (1972); M. MacCormick *et al.*, Phys. Rev. C 53, 41 (1996);
  O. Bartalini *et al.* (GRAAL), Phys. At. Nucl. 71, 75 (2008); D. O. Caldwell *et al.*, Phys. Rev. Lett. 40, 1222 (1978);
  S. Aid *et al.* (H1), Z. Phys. C 69, 27 (1995); S. Chekanov *et al.* (ZEUS), Nucl. Phys. B 627, 3 (2002); Meyer 70,
  Michalowski 77 and Caldwell 73 as in the PDG compilation; LEGS as cited in paper 1 (report BNL-64382, not found on INSPIRE).
- Polarized data: P. Pedroni *et al.* (A2), *Measurement of the Gerasimov-Drell-Hearn integrand for the proton and deuteron from
  200 to 1400 MeV*, Phys. Rev. C 114, 045203 (2026), arXiv:2604.14385; J. Ahrens *et al.* (GDH), Phys. Rev. Lett. 87, 022003 (2001);
  H. Dutz *et al.* (GDH), Phys. Rev. Lett. 91, 192001 (2003) and 93, 032003 (2004).
- Partial waves: D. Drechsel, S. S. Kamalov, L. Tiator, *Unitary Isobar Model – MAID2007*, Eur. Phys. J. A 34, 69 (2007);
  W. J. Briscoe *et al.*, *Extended SAID partial-wave analysis of pion photoproduction* (SM22), Phys. Rev. C 108, 065205 (2023);
  R. L. Workman *et al.*, *Unified Chew-Mandelstam SAID analysis of pion photoproduction data* (CM12), Phys. Rev. C 86, 015202 (2012).
  The SM22 cross sections were obtained from the public SAID interface (`data/2026/*said_sm22.dat`).
- Elastic channel: σ(γp → γp) of the NNLO amplitude of V. Lensky and V. Pascalutsa, *Predictive powers of chiral
  perturbation theory in Compton scattering off protons*, Eur. Phys. J. C 65, 195 (2010), as evaluated in the 2026 paper.
- Regge (polarized): S. D. Bass, M. Skurzok, P. Moskal, *Updating spin-dependent Regge intercepts*, Phys. Rev. C 98, 025209 (2018);
  variant: K. Helbing, *The Gerasimov-Drell-Hearn Sum Rule*, Prog. Part. Nucl. Phys. 57, 405 (2006), as in paper 2.

**Unverified inputs** (marked in the databases): the source and normalization (3 %) of the LEGS points; the normalization of
Bartalini 08 (5 %); the correlation of the Bass parameters (0.37, chosen to reproduce their quoted error); the error
composition of ELSA 2004; the independence of the MAMI 2001 and A2 2026 normalizations.

## License

Code: MIT (`LICENSE`). Data and result tables: CC BY 4.0 (`LICENSE-DATA`); the experimental data are
transcriptions of published measurements, cite their original sources.
