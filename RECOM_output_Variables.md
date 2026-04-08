# REcoM Model Output Variables

This document provides a comprehensive list of output variables available in the REcoM (Regulated Ecosystem Model) biogeochemical ocean model.

---

## Table of Contents

1. [Air-Sea Gas Exchange](#air-sea-gas-exchange)
2. [Carbon Chemistry](#carbon-chemistry)
3. [Atmospheric Inputs](#atmospheric-inputs)
4. [Benthic Processes](#benthic-processes)
5. [Phytoplankton Primary Production](#phytoplankton-primary-production)
   - [Nanophytoplankton](#nanophytoplankton)
   - [Diatoms](#diatoms)
   - [Coccolithophores](#coccolithophores)
   - [Phaeocystis](#phaeocystis)
6. [Zooplankton Grazing](#zooplankton-grazing)
   - [Microzooplankton](#microzooplankton)
   - [Mesozooplankton](#mesozooplankton)
   - [Macrozooplankton](#macrozooplankton)

---

## Air-Sea Gas Exchange

### CO2 Exchange

| Variable | Description | Units |
|----------|-------------|-------|
| `alphaCO2` | CO2 solubility coefficient in seawater | mol/kg/atm |
| `Kw` | Air-sea piston velocity (gas transfer velocity) | m/s |
| `dpCO2s` | Difference between oceanic and atmospheric pCO2 (ocean - atmosphere) | μatm |
| `pCO2s` | Partial pressure of CO2 in surface ocean | μatm |
| `CO2f` | CO2 flux into surface water (positive = ocean uptake) | mmolC/m²/d |

### O2 Exchange

| Variable | Description | Units |
|----------|-------------|-------|
| `O2f` | O2 flux into surface water (positive = ocean uptake) | mmolO/m²/d |

---

## Carbon Chemistry

| Variable | Description | Units |
|----------|-------------|-------|
| `Hp` | Mean concentration of H⁺ ions in surface water (ocean acidity) | mol/kg |

---

## Atmospheric Inputs

| Variable | Description | Units |
|----------|-------------|-------|
| `aFe` | Atmospheric iron deposition flux | μmolFe/m²/s |
| `aN` | Atmospheric dissolved inorganic nitrogen (DIN) deposition flux | mmolN/m²/s |

---

## Benthic Processes

Benthic compartments represent material accumulated at the seafloor from sinking particulate matter.

| Variable | Description | Units |
|----------|-------------|-------|
| `benN` | Nitrogen stored in benthic pool | mmol |
| `benC` | Carbon stored in benthic pool | mmol |
| `benSi` | Silicon stored in benthic pool | mmol |
| `benCalc` | Calcite (CaCO₃) stored in benthic pool | mmol |

---

## Phytoplankton Primary Production

### Nanophytoplankton

Small phytoplankton cells (typically 2-20 μm), representing the general small phytoplankton community.

| Variable | Description | Units |
|----------|-------------|-------|
| `NPPn` | Net primary production (carbon fixation minus respiration) | mmolC/m²/d |
| `GPPn` | Gross primary production (total carbon fixation) | mmolC/m²/d |
| `NNAn` | Net nitrogen assimilation rate | mmolN/m²/d |
| `Chldegn` | Chlorophyll degradation rate | 1/d |

### Diatoms

Large phytoplankton with silica frustules, important for export production.

| Variable | Description | Units |
|----------|-------------|-------|
| `NPPd` | Net primary production | mmolC/m²/d |
| `GPPd` | Gross primary production | mmolC/m²/d |
| `NNAd` | Net nitrogen assimilation rate | mmolN/m²/d |
| `Chldegd` | Chlorophyll degradation rate | 1/d |

### Coccolithophores

Phytoplankton that produce calcite plates (coccoliths), important for calcium carbonate cycle.

| Variable | Description | Units |
|----------|-------------|-------|
| `NPPc` | Net primary production | mmolC/m²/d |
| `GPPc` | Gross primary production | mmolC/m²/d |
| `NNAc` | Net nitrogen assimilation rate | mmolN/m²/d |
| `Chldegc` | Chlorophyll degradation rate | 1/d |

### Phaeocystis

Colonial phytoplankton forming gelatinous colonies, important in polar and coastal regions.

| Variable | Description | Units |
|----------|-------------|-------|
| `NPPp` | Net primary production | mmolC/m²/d |
| `GPPp` | Gross primary production | mmolC/m²/d |
| `NNAp` | Net nitrogen assimilation rate | mmolN/m²/d |
| `Chldegp` | Chlorophyll degradation rate | 1/d |

---

## Zooplankton Grazing

Grazing fluxes represent carbon transfer through the food web. Fluxes are given both as total (accounting for grazing efficiency) and by prey type (representing actual prey loss).

### Microzooplankton

Small heterotrophic protists (typically 20-200 μm) that graze on small phytoplankton.

| Variable | Description | Units |
|----------|-------------|-------|
| `grazmicro_tot` | Total grazing flux (accounting for grazing efficiency) | mmolC/m²/d |
| `grazmicro_n` | Grazing on nanophytoplankton (equals nanophyto loss to microzoo) | mmolC/m²/d |
| `grazmicro_d` | Grazing on diatoms (equals diatom loss to microzoo) | mmolC/m²/d |
| `grazmicro_c` | Grazing on coccolithophores (equals coccolith loss to microzoo) | mmolC/m²/d |
| `grazmicro_p` | Grazing on Phaeocystis (equals Phaeocystis loss to microzoo) | mmolC/m²/d |

### Mesozooplankton

Medium-sized zooplankton (typically 0.2-20 mm) including copepods and other crustaceans.

| Variable | Description | Units |
|----------|-------------|-------|
| `grazmeso_tot` | Total grazing flux (accounting for grazing efficiency) | mmolC/m²/d |
| `grazmeso_n` | Grazing on nanophytoplankton (equals nanophyto loss to mesozoo) | mmolC/m²/d |
| `grazmeso_d` | Grazing on diatoms (equals diatom loss to mesozoo) | mmolC/m²/d |
| `grazmeso_c` | Grazing on coccolithophores (equals coccolith loss to mesozoo) | mmolC/m²/d |
| `grazmeso_p` | Grazing on Phaeocystis (equals Phaeocystis loss to mesozoo) | mmolC/m²/d |
| `grazmeso_det` | Grazing on detritus group 1 (equals detritus 1 loss to mesozoo) | mmolC/m²/d |
| `grazmeso_mic` | Grazing on microzooplankton (equals microzoo loss to mesozoo) | mmolC/m²/d |
| `grazmeso_det2` | Grazing on detritus group 2 (equals detritus 2 loss to mesozoo) | mmolC/m²/d |

### Macrozooplankton

Large zooplankton (typically >20 mm) including krill and larger crustaceans.

| Variable | Description | Units |
|----------|-------------|-------|
| `grazmacro_tot` | Total grazing flux (accounting for grazing efficiency) | mmolC/m²/d |
| `grazmacro_n` | Grazing on nanophytoplankton (equals nanophyto loss to macrozoo) | mmolC/m²/d |
| `grazmacro_d` | Grazing on diatoms (equals diatom loss to macrozoo) | mmolC/m²/d |
| `grazmacro_c` | Grazing on coccolithophores (equals coccolith loss to macrozoo) | mmolC/m²/d |
| `grazmacro_p` | Grazing on Phaeocystis (equals Phaeocystis loss to macrozoo) | mmolC/m²/d |
| `grazmacro_mes` | Grazing on mesozooplankton (equals mesozoo loss to macrozoo) | mmolC/m²/d |
| `grazmacro_det` | Grazing on detritus group 1 (equals detritus 1 loss to macrozoo) | mmolC/m²/d |
| `grazmacro_mic` | Grazing on microzooplankton (equals microzoo loss to macrozoo) | mmolC/m²/d |
| `grazmacro_det2` | Grazing on detritus group 2 (equals detritus 2 loss to macrozoo) | mmolC/m²/d |

---

## Notes

### Variable Naming Convention
- Suffix `n` = nanophytoplankton
- Suffix `d` = diatoms
- Suffix `c` = coccolithophores
- Suffix `p` = Phaeocystis
- Suffix `mic` = microzooplankton
- Suffix `mes` = mesozooplankton
- Suffix `det` = detritus group 1
- Suffix `det2` = detritus group 2

### Unit Abbreviations
- `μatm` = micro-atmospheres (10⁻⁶ atm)
- `mmol` = millimoles (10⁻³ mol)
- `μmol` = micromoles (10⁻⁶ mol)
- `/d` = per day
- `/s` = per second
- `m²` = square meter

### Grazing Efficiency
Grazing fluxes without efficiency (`grazmeso_n`, `grazmicro_d`, etc.) represent the actual amount of prey consumed and therefore equal the loss rate of that prey type. The total grazing flux (`grazmeso_tot`, etc.) accounts for grazing efficiency, representing the carbon actually assimilated by the grazer (the rest becomes detritus or dissolved organic matter).

---

## Carbonate Chemistry System (3D Variables)

Detailed carbonate chemistry variables throughout the water column.

| Variable | Description | Units | Notes |
|----------|-------------|-------|-------|
| `PAR` | Photosynthetically Active Radiation | W/m² | Light energy available for photosynthesis (400-700nm wavelength), primary driver of phytoplankton growth |
| `CO2` | Aqueous CO2 concentration | mol/m³ | Dissolved CO2 in seawater, primary inorganic carbon form and substrate for photosynthesis |
| `pH` | Seawater pH | total scale | Measure of hydrogen ion concentration, affects carbonate chemistry and organism physiology (typically 7.5-8.5 in ocean) |
| `pCO2` | Partial pressure of CO2 | µatm | CO2 fugacity in equilibrium with atmosphere, determines air-sea CO2 flux direction and magnitude |
| `HCO3` | Bicarbonate ion concentration | mol/m³ | Major dissolved inorganic carbon species (~90% of DIC), pH buffer |
| `CO3` | Carbonate ion concentration | mol/m³ | Secondary DIC species, building block for calcium carbonate (CaCO₃) formation |
| `OmegaC` | Calcite saturation state (Ω) | dimensionless | Ratio of [Ca²⁺][CO3²⁻] to calcite solubility product; Ω>1 favors precipitation, reduced by ocean acidification |
| `kspc` | Calcite solubility product | mol²/kg² | Equilibrium constant for CaCO₃ dissolution, temperature/pressure dependent |
| `rhoSW` | In-situ seawater density | kg/m³ | Density at ambient temperature, salinity, and pressure; affects stratification and mixing |

---

## Particle Dynamics

Physical properties and sinking behavior of detrital particles.

| Variable | Description | Units | Notes |
|----------|-------------|-------|-------|
| `rho_det1` | Density of particle class 1 | kg/m³ | Mass per unit volume of small/slow-sinking particles; determines sinking speed via Stokes' law |
| `rho_det2` | Density of particle class 2 | kg/m³ | Mass per unit volume of large/fast-sinking particles; higher density = faster sinking = efficient carbon export |
| `scaling_visc` | Viscosity-based scaling factor | dimensionless | Adjusts particle sinking based on water viscosity (warmer water = lower viscosity = faster sinking) |
| `wsink_det1` | Sinking velocity of particle class 1 | m/s | Vertical descent rate of small detrital particles; typically 1-50 m/day for small aggregates |
| `wsink_det2` | Sinking velocity of particle class 2 | m/s | Vertical descent rate of large detrital particles; typically 50-200+ m/day for large aggregates/fecal pellets |

---

## Benthic Nutrient Fluxes (Extended)

Sediment-water interface nutrient release from remineralization processes.

| Variable | Description | Units | Notes |
|----------|-------------|-------|-------|
| `DIC_bf` | Dissolved Inorganic Carbon bottom flux | mmolC/m³ | Release of CO2/HCO3/CO3 from sediment remineralization; important in shallow seas |
| `DIN_bf` | Dissolved Inorganic Nitrogen bottom flux | mmolN/m³ | Release of NH4⁺/NO3⁻ from sediment remineralization; major N source in coastal/shelf systems |
| `Alk_bf` | Alkalinity bottom flux | meq/m³ | Release of acid-neutralizing capacity; dominated by NH4⁺ production and CaCO₃ dissolution |
| `DSi_bf` | Dissolved Silicate bottom flux | mmolSi/m³ | Release of Si(OH)4 from opal dissolution; critical nutrient for diatom growth |

---

## Zooplankton Respiration (Extended)

Metabolic CO2 production by different zooplankton size classes.

| Variable | Description | Units | Notes |
|----------|-------------|-------|-------|
| `respmicro` | Microzooplankton respiration rate | mmolC/m²/d | Small grazers (<0.2 mm): ciliates, heterotrophic dinoflagellates; rapid nutrient recycling |
| `respmeso` | Mesozooplankton respiration rate | mmolC/m²/d | Medium-sized grazers (0.2-2 mm): copepods, euphausiids, chaetognaths |
| `respmacro` | Macrozooplankton respiration rate | mmolC/m²/d | Large grazers (>2 mm): large euphausiids, amphipods, jellyfish; produces fast-sinking fecal pellets |

---

## Calcium Carbonate Dynamics

Calcification and dissolution processes.

| Variable | Description | Units | Notes |
|----------|-------------|-------|-------|
| `calcdiss` | Calcite (CaCO₃) dissolution rate | mmolC/m²/d | Conversion of solid CaCO₃ to dissolved Ca²⁺ and CO3²⁻; increases alkalinity, accelerated by ocean acidification |
| `calcif` | Calcification rate | mmolC/m²/d | Precipitation of calcium carbonate shells/tests by coccolithophores, foraminifera, pteropods; decreases alkalinity |

---

## Phytoplankton Aggregation

Formation of larger particles from small cells, converting slow-sinking cells to faster-sinking aggregates.

| Variable | Description | Units | Notes |
|----------|-------------|-------|-------|
| `aggn` | Nanophytoplankton aggregation | mmolC/m²/d | Collision-driven clumping of 2-20 µm cells; enhanced by stickiness (TEP), turbulence, high cell density |
| `aggd` | Diatom aggregation | mmolC/m²/d | Formation of diatom-dominated marine snow; major pathway for carbon export during bloom collapse; Si:C ballasting enhances sinking (100-200 m/day) |
| `aggc` | Coccolithophore aggregation | mmolC/m²/d | Formation of coccolith-bearing aggregates; CaCO₃ ballasting enhances export efficiency; creates "white waters" during blooms |
| `aggp` | Phaeocystis aggregation | mmolC/(m²*d) | Colony formation of colonial haptophyte; forms mucilaginous colonies up to 3 cm; dominant in polar/subpolar regions |

---

## Dissolved Organic Carbon (DOC) Excretion

Release of dissolved organic compounds by phytoplankton (passive/active release, typically 5-20% of gross production).

| Variable | Description | Units | Notes |
|----------|-------------|-------|-------|
| `docexn` | Nanophytoplankton DOC excretion | mmolC/m²/d | Fuels microbial loop; increases under nutrient stress or viral lysis |
| `docexd` | Diatom DOC excretion | mmolC/m²/d | Can be substantial during stationary phase or Si limitation; includes transparent exopolymer particles (TEP) |
| `docexc` | Coccolithophore DOC excretion | mmolC/m²/d | Released during calcification/growth; links organic and inorganic carbon cycles; may increase under high pCO₂ |
| `docexp` | Phaeocystis DOC excretion | mmolC/(m²*d) | Release of polysaccharides forming colonial matrix; major mucus producer; high C:N ratio, relatively refractory |

---

## Phytoplankton Respiration (Extended)

Metabolic CO2 release by autotrophs (maintenance metabolism and biosynthesis costs).

| Variable | Description | Units | Notes |
|----------|-------------|-------|-------|
| `respn` | Nanophytoplankton respiration | mmolC/m²/d | Determines net vs gross primary production; typically 10-30% of GPP, increases with temperature |
| `respd` | Diatom respiration | mmolC/m²/d | Metabolic costs of rapid growth and Si metabolism; higher in fast-growing bloom-forming species |
| `respc` | Coccolithophore respiration | mmolC/(m²*d) | Includes calcification energy costs (~20% additional energy vs non-calcifiers); may increase under acidification stress |
| `respp` | Phaeocystis respiration | mmolC/(m²*d) | Metabolic costs of colonial vs solitary forms; colonial forms may have lower specific respiration rates |

---

## Net Primary Production by Functional Group (3D)

Net carbon fixation after respiration losses (gross photosynthesis minus respiration and exudation).

| Variable | Description | Units | Ecological Role |
|----------|-------------|-------|-----------------|
| `NPPn3D` | Nanophytoplankton NPP | mmolC/m³/d | Sustained baseline production in oligotrophic waters; dominant in stratified, nutrient-poor conditions |
| `NPPd3D` | Diatom NPP | mmolC/m³/d | Bloom-forming, episodic production pulses; dominant in upwelling, mixed, nutrient-rich waters; r-strategists |
| `NPPc3D` | Coccolithophore NPP | mmolC/m³/d | Links organic and inorganic carbon pumps; bloom in stratified, post-bloom nutrient-depleted waters |
| `NPPp3D` | Phaeocystis NPP | mmolC/(m³*d) | High-latitude spring bloom specialist; tolerates low light, cold temperatures, ice-edge blooms |

---

## Temperature Effects on Photosynthesis

Temperature-dependent growth rate modifiers (Q10-based or Arrhenius growth rate multipliers).

| Variable | Description | Units | Notes |
|----------|-------------|-------|-------|
| `TTemp_phyto` | Nanophytoplankton temperature factor | per day | Growth rate scaling with temperature; often have broader temperature tolerance than specialists |
| `TTemp_diatoms` | Diatom temperature factor | per day | Accelerates growth in warmer waters; often have lower temperature optimum than flagellates |
| `TTemp_cocco` | Coccolithophore temperature factor | per day | Often have higher temperature optima than diatoms; may expand poleward with warming |
| `TTemp_phaeo` | Phaeocystis temperature factor | per day | Adapted to low temperatures (0-10°C optimum); climate warming may reduce competitive advantage |

---

## CO2 Effects on Photosynthesis

pCO2-dependent growth modifiers representing ocean acidification impacts on carbon concentration mechanisms (CCM).

| Variable | Description | Units | Notes |
|----------|-------------|-------|-------|
| `TPhyCO2` | Nanophytoplankton CO2 effect | per day | Some species benefit from higher CO2 (reduced CCM energy cost); generally small positive or neutral effect |
| `TDiaCO2` | Diatom CO2 effect | per day | Mixed responses; some species C3-like, others have efficient CCMs; high CO2 may reduce competitive advantage |
| `TCoccoCO2` | Coccolithophore CO2 effect | per day | High pCO2 may benefit photosynthesis but impair calcification; net effect species-specific, generally negative at >1000 µatm |
| `TPhaeoCO2` | Phaeocystis CO2 effect | per day | May benefit from CO2 enrichment in cold waters; limited experimental data compared to other groups |

---

## Nutrient Limitation Factors

Multi-nutrient limitation of photosynthesis following Liebig's Law of the Minimum (0-1 dimensionless multipliers).

| Variable | Description | Units | Notes |
|----------|-------------|-------|-------|
| `TqLF_phyto` | Nanophytoplankton nutrient limitation | per day | Minimum of N, P, Fe limitation; often N-limited in stratified waters, Fe-limited in HNLC regions |
| `TqLF_diatoms` | Diatom nutrient limitation | per day | Minimum of N, P, Si, Fe limitation; Si limitation unique to diatoms, controls bloom termination |
| `TqLF_cocco` | Coccolithophore nutrient limitation | per day | Minimum of N, P limitation (no Si requirement); can exploit low-nutrient conditions after diatom bloom collapse |
| `TqLF_phaeo` | Phaeocystis nutrient limitation | per day | Minimum of N, P, Fe limitation; can form luxury nutrient reserves in colonial form; often blooms early when nutrients replete |

---

## Light Limitation Factors

Photosynthesis-irradiance (P-I) curve response factors (0-1 dimensionless multipliers).

| Variable | Description | Units | Notes |
|----------|-------------|-------|-------|
| `TCphotLL_phyto` | Nanophytoplankton light limitation | per day | Reduces growth in low light, photoinhibition at high light; small cells acclimate faster to changing light |
| `TCphotLL_dia` | Diatom light limitation | per day | P-I curve with potentially higher light saturation point; often shade-adapted in coastal turbid waters; larger cells = slower photoacclimation |
| `TCphotLL_cocco` | Coccolithophore light limitation | per day | P-I curve for high-light adapted species; often prefer high irradiance, stratified waters; bloom in clear waters with deep light penetration |
| `TCphotLL_phaeo` | Phaeocystis light limitation | per day | P-I curve for low-light adapted polar species; low light saturation point; can grow under sea ice with <1% surface irradiance |

---

## Carbon-Specific Photosynthesis Rates

Actual photosynthetic carbon fixation rates (product of all limiting factors × maximum rate).

| Variable | Description | Units | Notes |
|----------|-------------|-------|-------|
| `TCphot_phyto` | Nanophytoplankton C-specific photosynthesis | per day (d⁻¹) | Gross carbon fixation per unit biomass; typically 0.5-3.0 d⁻¹ in natural populations |
| `TCphot_diatoms` | Diatom C-specific photosynthesis | per day (d⁻¹) | Can reach very high rates (>5 d⁻¹) during blooms; r-strategists with high maximum growth rates |
| `TCphot_cocco` | Coccolithophore C-specific photosynthesis | per day (d⁻¹) | Lower than diatoms but sustained in stratified conditions; additional energy for CaCO₃ formation reduces net efficiency |
| `TCphot_phaeo` | Phaeocystis C-specific photosynthesis | per day (d⁻¹) | Efficient at low temperature and light; colonial forms may have lower per-cell rates but higher bloom biomass |

---

## Silicate Assimilation

Diatom-specific nutrient uptake for frustule (shell) construction.

| Variable | Description | Units | Notes |
|----------|-------------|-------|-------|
| `TSi_assimDia` | Diatom silicate assimilation rate | per day (d⁻¹) | Uptake of Si(OH)4; Si:C ratio varies 0.1-1.0 depending on species and conditions; Si depletion terminates diatom blooms and shifts community to flagellates |

---

**Model Reference**: REcoM (Regulated Ecosystem Model) - Biogeochemical ocean model component  
**Last Updated**: 2025

