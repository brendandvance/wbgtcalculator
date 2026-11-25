## Intro

Calculator used to find the WBGT for various ADs in the Robert's Street Park assignment

## Method

- WBGT = 0.7 Tnw + 0.2 Tg + 0.1 Ta (ISO 7243)
- Natural wet-bulb: Hunter & Minyard (1999)
- Globe temperature: heat balance on a 150 mm black globe at 1.1 m, with a canopy in the sky view
- Sunlight: Haurwitz clear-sky model; wind: log wind profile from the ECCC July normal

## Sources

| Input | Source |
|---|---|
| Heat Warning threshold, 31 °C / humidex 40 | City of Toronto, Heat Relief Strategy (ECCC criteria) |
| Psychrometric wet-bulb | R. Stull (2011), J. Appl. Meteor. Climatol. 50, 2267–2269 |
| Clear-sky sunlight | B. Haurwitz (1945), J. Meteorol. 2, 154–166 |
| Direct/diffuse split, globe and ground properties, natural wet-bulb | C. H. Hunter and C. O. Minyard (1999), WSRC-MS-99-00757 |
| July wind at 10 m (3.6 m/s) | Environment Canada, Canadian Climate Normals 1981–2010 |
| Log wind profile, roughness z0 = 0.5 m | WMO Guide to Instruments and Methods of Observation (WMO-No. 8) |
| Globe height 1.1 m, WBGT formula | ISO 7243:2017; AIHA |
| Membrane transmission 4–10%, absorption 13–22% | Birdair, “Tensile membrane structures: energy use and lighting” |
| Light panel reflectance ≥ 0.65 | ENERGY STAR cool roof criteria |
| Canopy wind coefficient h = 5.7 + 3.8u | J. A. Duffie and W. A. Beckman (2013), Solar Engineering of Thermal Processes |
| Mist cooling 2–4 °C at 40 °C / 15% RH | Arizona State University (2023), outdoor misting study |
| Terracotta tube cooling 2 °C | Ibrahim et al., porous ceramic evaporative cooling |
| Ground temperatures in sun and shade | Team 207 site readings, Oct 12, 2026 (CDS Table A2) |

Hunter and Minyard print the natural wet-bulb solar coefficient as 0.021. The script uses 0.0021, because 0.021 contradicts the paper's own statement that Tnw stays between Tw and Ta.
