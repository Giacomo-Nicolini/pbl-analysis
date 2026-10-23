# PBLH literature: annotated notes

Search of 2026-10-15 (Cowork, web search + Crossref/publisher verification); second-wave papers checked against the
PDFs on 2026-10-22. BibTeX: `pblh_literature_2026-10.bib`
(50 entries, all verified; author lists marked "and others" to be completed). Kljun App. B references:
`pblh_collection.bib`. Footprint models: `footprint_collection.bib`.
Priorities: **K** = key (read in full, cite), U = useful, B = background. OA = open access.
Search was in English; Chinese-language literature and conference-only work were not covered.

## 1. Positioning: is our study new?

- The Kljun et al. (2015) App. B PBLH is already used at network scale **without evaluation**:
  - ICOS: `Rebmann2018` states ICOS applies Kljun 2015, which needs h. The ETC tool `icos-etc/FFP-utils` (GitHub,
    grey literature) uses an internally modelled PBLH by default, with an optional ERA5 hourly file interpolated to
    half-hours.
  - FLUXNET2015 footprint climatologies (`Fang2024`): h "calculated by methods from Kljun"; records with zm ≥ 0.8h
    removed. No sensitivity test.
  - AmeriFlux footprint representativeness (`Chu2021`, in footprint bib): h following Kljun 2015 from BASE data. Tests
    of the 80 % contour and σv, not of h. Night footprints ≈ 45 % longer than day.
  - Campbell EasyFlux (`Zhou2026`): Nieuwstadt + B&G (Newton–Raphson), no validation of h (declined when a referee
    asked). Describes γ as the "dry adiabatic lapse rate (commonly 9.8 × 10⁻³ K m⁻¹)" (checked in the PDF):
    numerically ≈ our default, conceptually B&G's γ is the free-atmosphere θ gradient.
  - ERA5 used directly as h input in recent studies (`Moutahir2026`).
- Closest prior work:
  - `Beamesderfer2023` (K): ceilometer MLH vs EC H, LE, VPD at 5 AmeriFlux sites; surface variables explain 38–76 % of
    daily-max MLH; NARR PBLH ≈ 2× ceilometer at semi-arid sites. Statistical, daytime maximum only, no slab model.
  - `Chu2026` (K): ML on cumulative H, LE since sunrise + lower-tropospheric stability → Doppler-lidar CBL height at ARM
    SGP (R² ≈ 0.84). Convective only, no physical benchmark, no ERA5.
  - `Sinclair2022` (K): FI-Hyy (ICOS). ERA5 vs radiosondes and MWR, stratified by EC z/L. No PBLH computed from EC.
  - `Helbig2021` (K): community review; PBLH observations "not routinely shared through FLUXNET or the regional
    observation networks"; footprints depend on ABLH, most at night, early morning, tall towers and zm > 0.1h. Names
    the gap we fill. Lists Swedish ICOS sites (SE-Htm, SE-Nor, SE-Svb) with ceilometers next to EC towers.
- **Conclusion:** no published study evaluates EC-derived PBLH against ERA5 and measurements across a network and over
  all stabilities, nor propagates realistic, stability-dependent h errors through FFP climatologies.

## 2. ERA5 BLH: definition and known errors

- Definition (checked in the PDF, 2026-10-17): `ECMWF2016` Sect. 3.10.1, Eq. 3.90. A post-processing diagnostic,
  independent of the turbulence scheme: the lowest level where the bulk Richardson number between the lowest model
  level and that height reaches 0.25, computed from virtual dry static energy; surface friction ignored and 2 m wind
  set to zero (Vogelezang & Holtslag 1996 algorithm as simplified by `Seidel2012`); linear interpolation between
  levels. **No surface parcel excess** (an earlier version of these notes said otherwise; corrected). The mixed-layer
  height used inside the turbulence scheme comes from an entraining parcel model and is not the output BLH. ERA5 =
  Cy41r2, 0.25°, 137 levels (`Hersbach2020`). Same formulation as `Seidel2012` (Beljaars co-author). A single Ri_c
  gives a stability-dependent bias; optimal values 0.24 strongly stable, 0.31 weakly stable, 0.39 unstable, from four
  campaigns with PBLH from turbulence, θ-gradient, LLJ and parcel methods (`Zhang2014`; values depend on the lower
  reference height z_s: 0.24 for z_s = 40 m, 0.31 for 80 m, 0.39 for the surface-layer top); ≈ 0.12 over Arctic sea
  ice (`Akansu2023`).
- Convective daytime: near-global radiosondes, ERA5 best of four reanalyses, ≈ −130 m mean bias (`Guo2021`, K); land
  2011–2021, 370 radiosonde stations, all launch times: radiosonde − ERA5 = +95.7 m on average (ERA5 too shallow),
  150 m in MAM to 64 m in SON, 309 m at 17 LST to 1.8 m at 00 LST, s.d. of the difference 807 m at 17 LST (`Guo2024`,
  U; radiosonde PBLH by bulk Ri = 0.25). Guo2024 also provides an ERA5 product bias-corrected by random forest
  (3-hourly, 0.25°, land, 2011–2021; Zenodo 10.5281/zenodo.6498004). In Europe the sign varies: ERA5 too deep at 12 UTC Europe-wide
  (`Arganaraz2026`, K: RMSE 235–470 m), slightly too shallow at Potenza (`SalcedoBosch2026`, K: bias −109 m by day),
  good agreement at Hyytiälä (`Sinclair2022`), within 10–35 % of ALICENET ceilometer MAL (`Bellini2025`, U).
- Stable / transitions: ERA5 too deep in very stable conditions and poor in the evening transition (18 UTC)
  (`Sinclair2022`); too deep over sea ice (`Xi2024`, U: MOSAiC soundings with the same Ri = 0.25 criterion; centred
  RMSE 201 m and r = 0.61 for sea-ice concentration > 90 %; differences largest in neutral layers and the cold season;
  ERA5 BLH errors follow its temperature and wind-profile errors more than humidity); night RMSE small in absolute terms (≈ 60–100 m) but ERA5
  explains only ≈ 10 % of radiosonde night variance (`Arganaraz2026`); reanalysis deeper than sondes mainly at night
  (`Seidel2012`, K).
- Method of the reference matters as much as ERA5 error: IGRA gradient-method BLH ≈ 1150 m (night) / 1750 m (day) vs
  ERA5 ≈ 350 / 1150 m in Europe, mostly definitional (`Madonna2021`, K); radiosonde methods differ by hundreds of
  metres (`Seidel2010`, K).
- Complex terrain and coast: largest errors, often daytime overestimation by hundreds of m to > 1 km (`Wei2025`,
  `CollaudCoen2014` for COSMO-2, `Chen2023`) — all B.
- Implication for us: expect regime-dependent ERA5 errors, relatively largest in shallow stable PBLs, i.e. where FFP is
  most sensitive to h. Neither ERA5 nor EC-derived h can be treated as truth in that regime.

## 3. Measured PBLH: what each method class sees

- Reviews: `Kotthaus2023` (K, the anchor for our "different tops" argument), `Emeis2008` (B), `Dang2019` (B, aerosol
  lidar algorithms).
- Turbulence-based (Doppler lidar σw/ε, sodar): tracks the currently mixing layer → natural reference for EC-derived
  PBLH. Doppler-lidar ABL classification `Manninen2018` (K, basis of the Cloudnet product) and precursor `Harvey2013`.
- Aerosol-based (ceilometer, lidar): agrees in the mature CBL; late in the morning, high in the afternoon and at night
  (residual layer). Ceilometer vs Doppler lidar at JOYCE: +300 m (morning) to +600 m (late afternoon) (`Schween2014`, K);
  European ALC evaluation (`Haeffelin2012`, K). Algorithms: STRAT (`Morille2007`), Pathfinder (`deBruine2017`),
  PathfinderTURB (`Poltera2017`: MeteoSwiss, Lufft CHM15k; CBL height and continuous aerosol layer; mean bias 27 m vs
  expert manual detection and 53 m vs radiosonde bulk Ri at Payerne), CABAM (`Kotthaus2018`: first open algorithm for
  the Vaisala CL31, handles precipitation and classifies clouds; good against AMDAR inversions in the afternoon, while
  at night and in the morning transition the inversion usually marks the residual-layer top — relevant for the
  Vaisala CL51 at SE-Nor).
- Thermodynamic (radiosonde bulk Ri/parcel, MWR): closest in definition to ERA5. Regime-aware sounding algorithm
  `Liu2010` (U): regime from θ(level 5) − θ(level 2) against δs (1.0 K over land); CBL and neutral layers by a parcel
  scan from above 150 m agl (δu = 0.5 K, then first θ gradient ≥ 4 K km⁻¹), SBL by the top of the surface inversion
  or the low-level-jet nose, whichever is lower; uses θ, not θv; checked against lidar and sodar; 58,286 soundings
  from 14 campaigns. IGRA derived parameters (`Durre2008`, U; PDF not available, entry checked against Crossref only). MWR best vs sondes at Payerne (`CollaudCoen2014`).
- Multi-instrument sites collocated with ICOS: Hyytiälä (`Sinclair2022`, `Franck2021`), Swedish sites (`Helbig2021`).
- Recommendation (to discuss): primary validation in daytime convective conditions (EC z/L < 0); turbulence-based
  height (or ceilometer MLH only late morning–mid afternoon) as reference for EC-derived PBLH; radiosonde bulk Ri for
  ERA5; stable and transition hours reported separately with a regime classification (`Liu2010`, `Manninen2018`).
  This matches the existing plan ("Comparing different method classes").

## 4. Parameterisations

- Stable BL height from surface fluxes: large scatter, site-dependent coefficients; surface bulk Ri best when profiles
  exist; well-defined h only ≈ 20 % of the time on CASES-99 (`Vickers2004`, K). Zilitinkevich–Mironov offset for
  shallow layers; dimensional-analysis formula without f performs best (`Steeneveld2007`, K; detail from the conference
  version, check in the journal paper). `Zilitinkevich2002` (U), on Cabauw data: the critical bulk Ri depends on the
  free-flow stability N (fit Ri_Bc = 0.1371 + 0.0024 N/|f|); among diagnostic SBL equations their multi-limit h_CQE
  (u*, L, N, f) fits best (r = 0.60), the Zilitinkevich-scale formula h = 0.4 u*² |f B_s|^(−1/2) — the same scaling as
  our B1 Nieuwstadt — the worst (r = 0.27) (Table IV; quote correlations only: some bias/rms entries in that table
  are internally inconsistent); a relaxation (prognostic) equation towards h_CQE is recommended. Small-sample confirmations: `Sun2021`, `Akansu2023` (B).
  → Relevant to the open point "bound B1 by B3 on windy stable nights".
- Convective slab: order of magnitude right but strongly γ-dependent; without subsidence/advection a slab model
  overestimates by ≈ 70 % (BLLAST case, `Pietersen2015`, U). Global evaluation of a slab model (CLASS4GL) against
  18,385 morning–afternoon IGRA sounding pairs from 134 stations: bias of the ML growth rate 10.1 m h⁻¹, r = 0.53
  (`Wouters2019`, U); initial γθ from the sounding levels above h, initial jump from the capping-inversion value minus
  the mixed-layer value → recipe for our γ option 2/3 from ERA5 profiles. Surface heterogeneity
  and γ can cancel H effects (`Brugger2018`, U).
- Pre-processors (AERMET, OML): no recent peer-reviewed evaluation found (EPA reports, HARMO proceedings only).
  `Helmis2012` (B): model PBLH generally higher than sodar-RASS/ceilometer in Athens.

## 5. Footprint sensitivity to h

- The only formal test: `Kljun2015` Table 5 — h ± 20 % → peak location ∓ ≈ 3.7 % and peak value ± ≈ 3.6 % in stable
  conditions, < 1 % neutral/convective; single footprints only, not climatologies. Validity limits 20 z0 < zm < 0.8h,
  zm/L ≥ −15.5.
- Network applications compute h with Kljun formulas and never test it (`Chu2021`, `Fang2024`; ICOS per `Rebmann2018`,
  whose target-area criteria — ≥ 90 % in stable conditions for ≥ 50 % of the time, ≥ 80 % cumulative — are exactly the
  metrics our phase 5 perturbs).
- Tracer / LES validations (all in `footprint_collection.bib`: `Markkanen2009`, `VandeBoer2013`, `Heidbach2017`,
  `Arriga2017`, `Kumari2020`, `Coates2021`) quantify structural model error (tens of %), mostly daytime/unstable, never
  with h as an uncertain input. `VandeBoer2013` is a methodological precedent for ensemble footprint uncertainty.
- Grey literature: McCombs & Hiscox, AMS 2017 poster — FFP insensitive to h for zm < ≈ 10 m. No journal version found.
- Gap: realistic, stability-dependent h errors (from our phases 3–4) propagated through FFP climatologies, with
  effects on fetch, target-area fraction and validity flags. Tests whether the 2015 "minor sensitivity" holds at
  realistic error ranges.

## 6. Reference datasets for measured PBLH

| Dataset | What | Coverage | Notes |
|---|---|---|---|
| IGRA v2.2 (NOAA NCEI) | soundings + derived `MIXHGT` (parcel), `INVHGT` | > 2800 stations, mostly 00/12 UTC | public; low vertical resolution limits bulk Ri (`Seidel2010`) |
| E-PROFILE (EUMETNET) | ALC attenuated backscatter, cloud base | ≈ 280 ceilometers, 22 countries, since ≈ 2013 | no operational MLH product: run STRAT/Pathfinder/CABAM ourselves; archive at CEDA |
| ACTRIS Cloudnet | Doppler lidar, ceilometer, MWR, cloud radar | ACTRIS sites (Hyytiälä, Jülich, SIRTA, Cabauw…) | CC BY 4.0; check the portal for the Doppler-lidar ABL product status |
| ALICENET (Italy) | ceilometer backscatter; MAL/CAL heights | Italian network | layer products not yet public; relevant for Italian ICOS sites |
| Guo2024 product | merged radiosonde/ERA5/GLDAS PBLH | global land, 3-hourly, 0.25° | extra comparison dataset (ESSD) |

## 7. Papers to download into `C:\WORK\bibliography_new_win11\`

Open-access papers I can read online; PDFs are still welcome. **Priority downloads (key, not open access or
uncertain):**

1. `Seidel2012` — JGR, 10.1029/2012JD018143
2. `Seidel2010` — JGR, 10.1029/2009JD013680
3. `Vickers2004` — JAM, 10.1175/JAM2160.1
4. `Steeneveld2007` — JAMC, 10.1175/JAM2454.1
5. `Haeffelin2012` — BLM, 10.1007/s10546-011-9643-z
6. `Helbig2021` — AFM, 10.1016/j.agrformet.2021.108509
7. `Beamesderfer2023` — AFM, 10.1016/j.agrformet.2023.109687
8. `Manninen2018` — JGR, 10.1029/2017JD028169
9. `ECMWF2016` — IFS Cy41r2 Part IV (downloaded 2026-10-17; BLH: Sect. 3.10.1, Eq. 3.90)

Key open-access papers (PDF optional): `Kotthaus2023`, `Sinclair2022`, `Guo2021`, `Madonna2021`, `Arganaraz2026`,
`SalcedoBosch2026`, `Schween2014`, `Chu2026`, `Fang2024`, `Rebmann2018`, `Hersbach2020`.

Second wave (useful): `Zilitinkevich2002`, `Wouters2019`, `Pietersen2015`, `Liu2010`, `Zhang2014`, `Guo2024`,
`Xi2024`, `Poltera2017`, `Kotthaus2018` — downloaded and checked against the PDFs on 2026-10-22 (notes above
updated). `Durre2008`: no access to the PDF; bibliographic entry from Crossref.

Naming convention in the folder: `Author_Year_Journal_Title.pdf`.

## 8. Left out / to check

- Not verified, not in the bib: Weidinger et al. 2008 (Időjárás; surface-flux mixing height vs AERMOD vs sondes in
  Hungary; convective growth sensitive to γ). Zhang, Seidel & Zhang 2013 (J. Climate, PBLH trends over Europe).
- Excluded as not final: Mulet-Benzo et al. (WES preprint wes-2025-155, lidar vs ERA5 offshore/coastal).
- `footprint_collection.bib`: entry `McCombs2019` (BLM, point-to-grid conversion) carries the DOI of `Arriga2017`
  (10.1016/j.agrformet.2017.04.006). To fix in the reference manager.
