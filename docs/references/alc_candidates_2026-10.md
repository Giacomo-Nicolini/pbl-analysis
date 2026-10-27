# Ceilometer / lidar candidates near ICOS ecosystem stations (search of 2026-10-26)

Scope: the 40 class 1–2 ICOS ecosystem stations with target areas in `config/sites.csv`. Sources: CEDA catalogue
records of E-PROFILE stations (one per instrument, coordinates from the record), ALICENET (Bellini et al. 2024, AMT,
Table 1), FMI network paper (Hirsikko et al. 2014, AMT 7:1351), national service pages. Distances are great-circle
from those coordinates. Distance classes as in `docs/project_plan.md` (< 5 km all regimes; 5–25 km daytime
convective only). Searches by web (e-profile.eu and CEDA directory listings not reachable), so the list may be
incomplete where noted.

**Main finding: no operator publishes an operational mixing-layer-height product.** E-PROFILE distributes cloud base
and attenuated backscatter only (EUMETNET; ICOS/ACTRIS/E-PROFILE testbed still evaluating MLH algorithms, EGU22-7378;
KNMI Deep-Pathfinder in implementation, EMS2025-700). RMI has an internal research MLH algorithm; ALICENET MAL/CAL
layers are "planned to be released". Apart from Cloudnet ε (turbulence) and the Lanzhot BL-View L3 product, every
candidate below means running an aerosol-MLH retrieval ourselves on backscatter.

## Within 5 km

| ICOS | Instrument site | Operator | Type | km | Notes |
|---|---|---|---|---|---|
| FI-Hyy | Hyytiälä SMEAR II | FMI / UHel | HALO Doppler lidar (Cloudnet) | 0.6 | done (spec 07) |
| SE-Nor | Norunda | ICOS/ACTRIS | CL51 (Cloudnet) | 0.1 | backscatter only |
| SE-Svb | Svartberget | SLU (E-PROFILE since 2024-12; SLU archive 2018–24) | CL51 | 0.0 | backscatter only |
| CH-Dav | Davos | MeteoSwiss (E-PROFILE, since 2021-07) | CHM15k | 1.0 | no MLH in record; MeteoSwiss has PathfinderTURB (Poltera2017) — ask; alpine valley |
| FI-Sod | Sodankylä | FMI | HALO Doppler lidar + CT25K | 1.2 | research data, confirmed to 2019 only; check Cloudnet |
| FI-Lom | Kenttärova | FMI (E-PROFILE) + Cloudnet | CL31; Cloudnet ε (spec 07 inventory) | 1.8 | Cloudnet ε 2022-08 → 2023-12, from 2026-04 |
| FI-Sii | Hyytiälä SMEAR II | FMI / UHel | HALO Doppler lidar (Cloudnet) | 5.2 | just over 5 km |
| DE-RuS | Jülich JOYCE | Univ. Cologne (E-PROFILE + Cloudnet) | CHM15k; Cloudnet ε | 5.3 | daytime convective only (rule) |
| CZ-Lnz | Lanzhot grassland | CzechGlobe (E-PROFILE) | CL61, BL-View L3 MLH | < 5 (to confirm) | data offered by PI (2026-10) |

## 5–25 km (daytime convective only)

| ICOS | Instrument site | Operator | Type | km | Notes |
|---|---|---|---|---|---|
| CZ-BK1 | Lysá hora | CHMI | CL31 | 8.1 | mountain summit (≈ 1320 m) vs tower at 881 m: probably not representative |
| FR-Gri | Trappes | Météo-France | Mini-MPL | 9.1 | closer than SIRTA |
| SE-Deg | Svartberget | SLU | CL51 | 13.4 | backscatter only |
| BE-Maa | Diepenbeek | RMI | CL51 | ≈ 16 if in Diepenbeek | CEDA coordinates (50.55, 5.27) look wrong: verify (OSCAR) |
| FR-Fon | Melun | Météo-France | CL31 | 16.6 | |
| IT-Cp2 | Castel di Guido | CNR-ISAC (ALICENET) | CHM15k | 21.9 | coastal; ALICENET layers not released |
| FR-Gri | SIRTA Palaiseau | IPSL (Cloudnet) | CHM15k, CL31, Doppler lidar | 23.4 | Cloudnet ε from 2025-01 |
| DE-HoH | Magdeburg | DWD | CHM15k | ≈ 24.5 (coords ±1 km) | borderline |
| BE-Dor | Humain | RMI | CL51 | ≈ 24–25 | borderline; RMI internal MLH |
| FR-Lus | Poitiers | Météo-France | CL31 | 24.9 | borderline |

## Nothing within 25 km (nearest)

BE-Bra (Uccle 58), BE-Lon (Uccle 39), BE-Vie (Humain 54), DE-Geb (Harzgerode 63; DWD network larger than
E-PROFILE: ask DWD), DE-Tha (Chemnitz 52; idem), DK-Sor (none found; listing incomplete), NL-Loo (Cabauw 60),
FR-Bil (Bordeaux 43), FR-FBn (Hyères 41), FR-Hes (Strasbourg 45), FR-Lqu (Clermont-Ferrand 36), FR-Pue (Montpellier
35), IT-BCi (Potenza 66), IT-MBo and IT-Ren (> 150; Innsbruck-area ceilometers not checked), IT-SR2 (Monte Cimone
61), UK-AMo (Eskdalemuir 53), SE-Htm (Hörby 30), SE-Myc (Kroppefjäll 27), SE-Sto (Naimakka 107), NO-Hur (Oslo 52),
GL-ZaF (none).

## Low-effort follow-ups

1. Cloudnet inventory: add Sodankylä (FMI HALO) if it is a Cloudnet site.
2. MeteoSwiss: is PathfinderTURB (or another MLH) run on the Davos CHM15k, and is it available?
3. RMI: MLH from Humain and Diepenbeek on request; verify the Diepenbeek coordinates.
4. DWD: any ceilometer (in or outside E-PROFILE) near DE-Geb, DE-Tha; MLH product?
