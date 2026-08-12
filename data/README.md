# Input data — provenance, licence and preparation

Everything in this directory is an **input**. Files the code computes are written to
`../results/`, figures and tables to `../output/`.

| File | Content | Source | Licence |
|---|---|---|---|
| `NHIS.xlsx` | daily influenza care episodes by date, province, sex and age band, National Health Insurance Service (NHIS) of Korea | [Open Government Data Portal, dataset 15089429](https://www.data.go.kr/data/15089429/fileData.do) | "이용허락범위 제한 없음" — no restriction on use. **Not included in this repository**, and neither is the tidy table derived from it (`../results/NHIS_Flu.csv`): download the workbook and run step 0 (below) |
| `Pops_Dec2022.csv` | resident registration population by single year of age, December 2022 | Ministry of the Interior and Safety (MOIS), resident registration statistics | public statistics, attribution requested. **EUC-KR encoded**, not UTF-8 |
| `contact_prem2017.csv` | 6 × 6 contact matrix for Korea, aggregated from Prem et al. (2017) | Prem K, Cook AR, Jit M. *PLoS Comput Biol* 2017;13:e1005697 (supplementary `MUestimates_all_locations_*.xlsx`) | derived product; the 16 × 16 original is **not** redistributed here |
| `contact_prem2021.csv` | 6 × 6 contact matrix for Korea, aggregated from Prem et al. (2021) | Prem K, et al. *PLoS Comput Biol* 2021;17:e1009098 (`output/syntheticmatrices/contact_all.rdata`, entry `KOR`) | as above |

The survey-based Korean contact matrix (Son et al. 2025) used for the main analysis is a
6 × 6 matrix written directly in the source files rather than stored here.

---

## NHIS workbook: download and rename

1. Open <https://www.data.go.kr/data/15089429/fileData.do> and download the XLSX file.
   Its published name is Korean —
   `국민건강보험공단_감염성질환(인플루엔자) 의료이용정보_<YYYYMMDD>.xlsx`
   (National Health Insurance Service — infectious disease (influenza) healthcare
   utilisation, followed by the cut-off date).
2. **Rename it to `NHIS.xlsx`** and place it in this directory. The Korean file name is not
   used anywhere in the code, so the analysis runs on systems that cannot display Hangul.
3. Convert it to the tidy table the rest of the pipeline reads:

   ```bash
   julia --project=. code/00_prepare_nhis.jl     # data/NHIS.xlsx -> results/NHIS_Flu.csv
   ```

   The step takes about 20 seconds and prints the row count, the date range and the case
   totals per age band, which is the quickest check that the right workbook was picked up.

### Which release of the workbook

The analysis in the paper used the release with a cut-off of **2023-12-31**. Neither the workbook
nor the table derived from it is redistributed here, so re-running the NHIS analysis starts with
downloading the file — and the portal is updated periodically, so what you download today is a
later release with more rows.

A newer release runs through the same code without modification, and the seasons analysed here
(2009-10 to 2018-19 and 2022-23) are unaffected by rows added after 2023. It will not, however,
reproduce the published numbers to the last digit. `code/00_prepare_nhis.jl` prints the
fingerprint below and warns when it does not match, so you always know which release you are on:

| | release used in the paper |
|---|---|
| rows | 820,906 |
| date range | 2006-01-01 to 2023-12-31 |
| total care episodes | 23,656,347 |
| cases by age band 1-6 | 3,386,152 / 5,113,473 / 3,556,167 / 6,560,712 / 3,308,547 / 1,731,296 |
| MD5 of the `results/NHIS_Flu.csv` it produces | `9aa2a7204046a77b4ff768a0e683f8eb` |

### Workbook layout

Sheet 1 (outpatient and inpatient, primary and secondary diagnosis) is the one used.
`code/00_prepare_nhis.jl` reads it by column position rather than by the Korean headers:

| Column | Content |
|---|---|
| 1 | treatment start date (yyyy-mm-dd) |
| 2 | province code (11 = Seoul, …) |
| 3 | sex |
| 4 | age band, written as `<digit>. <label>`, digits 1–6 = 0-5, 6-11, 12-17, 18-44, 45-64, 65+ |
| 5 | number of care episodes |

The analysis pools provinces and sexes, so only date, age band and count are used
downstream; the other two columns are carried into `results/NHIS_Flu.csv` so that the CSV
remains a faithful long-format copy of the source.

---

## Rebuilding the Prem contact matrices

```bash
python3 code/make_prem_contact.py --prem2017 <dir with MUestimates_*.xlsx> \
                                  --prem2021 <dir with contact_all.rdata>
```

Two transformations are applied, and both matter:

1. **Aggregation.** The originals use sixteen 5-year bands; they are aggregated to the six
   bands of this study weighted by the single-year population in `Pops_Dec2022.csv`.
2. **Transposition.** The convention here is `C[i, j]` with **column = infector** (the force
   of infection is `C * I`), whereas Prem reports `row = ego (participant)`. This was checked
   against the reciprocity relation `C_ij N_j = C_ji N_i`: the violation is 0.011 in the
   column convention against 1.02 in the row convention.

## Data availability

The claims data are public at the portal link above and the population data are published by
MOIS. No individual-level data are used: every input is an aggregated daily count.
