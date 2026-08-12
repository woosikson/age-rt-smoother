# age-rt-smoother / code/00_prepare_nhis.jl — STEP 0: NHIS workbook -> tidy CSV
# ─────────────────────────────────────────────────────────────────────────────
# Reads the National Health Insurance Service (NHIS) influenza workbook
#   data/NHIS.xlsx      (see data/README.md for the download link and the rename step)
# and writes the tidy long table every later step consumes:
#   results/NHIS_Flu.csv    columns: date, adm1, sex, agg, N
#
# Workbook layout (as published; column order is used rather than the Korean headers)
#   sheet 1  outpatient + inpatient, primary + secondary diagnosis   <- the sheet we use
#   sheet 2  inpatient only, primary + secondary diagnosis
#   sheet 3  province code table
#   columns  1 treatment start date (yyyy-mm-dd)
#            2 province code (11 = Seoul, ...)
#            3 sex, written in Korean
#            4 age band, written as "<digit>. <label>" with digits 1-6
#              1 = 0-5, 2 = 6-11, 3 = 12-17, 4 = 18-44, 5 = 45-64, 6 = 65+
#            5 number of care episodes
#
# The analysis pools provinces and sexes, so only date, age band and count are used
# downstream; adm1 and sex are kept in the CSV so that the file stays a faithful
# long-format copy of the source.
#
# Run: julia --project=. code/00_prepare_nhis.jl
# ─────────────────────────────────────────────────────────────────────────────
using XLSX, CSV, DataFrames, Dates, Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const DDIR = joinpath(ROOT, "data")
const RDIR = joinpath(ROOT, "results")
isdir(RDIR) || mkpath(RDIR)

const XLSX_PATH = joinpath(DDIR, "NHIS.xlsx")
const OUT_PATH  = joinpath(RDIR, "NHIS_Flu.csv")

isfile(XLSX_PATH) || error("""
    NHIS.xlsx not found at $XLSX_PATH
    Download the workbook from https://www.data.go.kr/data/15089429/fileData.do ,
    rename it to NHIS.xlsx and put it in data/ (see data/README.md).""")

# The male label starts with U+B0A8 and the female label with U+C5EC. Comparing the
# first character keeps this source file pure ASCII while still parsing the data.
const MALE_FIRST_CHAR = '\ub0a8'   # first character of the male label in the source file
sex_code(s) = first(strip(String(s))) == MALE_FIRST_CHAR ? "M" : "F"

# Age band cell looks like "<digit>. <label>"; the leading digit is the band index.
band_index(s) = parse(Int, first(strip(String(s))))

println("reading $(basename(XLSX_PATH)) ...")
raw = XLSX.readtable(XLSX_PATH, 1) |> DataFrame     # sheet 1, first row is the header
@printf("  %d rows x %d columns\n", nrow(raw), ncol(raw))
ncol(raw) == 5 || error("expected 5 columns in sheet 1, found $(ncol(raw))")

out = DataFrame(
    date = Date.(string.(raw[!, 1])),
    adm1 = Int.(raw[!, 2]),
    sex  = sex_code.(raw[!, 3]),
    agg  = band_index.(raw[!, 4]),
    N    = Int.(raw[!, 5]))

sort!(out, [:date, :adm1, :sex, :agg])
all(1 .<= out.agg .<= 6) || error("age band index outside 1-6")

CSV.write(OUT_PATH, out)
@printf("wrote %s  (%d rows, %s to %s, %d cases)\n", "results/" * basename(OUT_PATH),
        nrow(out), string(minimum(out.date)), string(maximum(out.date)), sum(out.N))
@printf("  cases by age band 1-6: %s\n",
        join([string(b, "=", sum(out.N[out.agg .== b])) for b in 1:6], " "))

# The portal re-publishes this dataset periodically. The paper used the 2023-12-31 release; a
# later one runs through the same code but will not reproduce results/NHIS_Flu.csv byte for byte.
const PAPER_ROWS, PAPER_LAST, PAPER_CASES = 820906, Date(2023, 12, 31), 23656347
if nrow(out) != PAPER_ROWS || maximum(out.date) != PAPER_LAST || sum(out.N) != PAPER_CASES
    @warn """This is not the release used in the paper (2023-12-31 cut-off:
             820,906 rows, last date 2023-12-31, 23,656,347 care episodes).
             The analysed seasons are unaffected, but the published numbers are reproduced
             exactly only from the shipped results/NHIS_Flu.csv. See data/README.md."""
end
