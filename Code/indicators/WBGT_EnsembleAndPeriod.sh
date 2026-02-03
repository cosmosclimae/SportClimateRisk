#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

ROOT="/mnt/c/Data/WBGT"
OUT="${ROOT}/ENSEMBLE_Foot"
TMP="${OUT}/_tmp2"
mkdir -p "${OUT}" "${TMP}"

MODELS=(
  "GFDL-ESM4"
  "IPSL-CM6A-LR"
  "MPI-ESM1-2-HR"
  "MRI-ESM2-0"
  "UKESM1-0-LL"
)

SSPS=( "ssp126" "ssp245" "ssp585" )

log(){ echo "[INFO] $*"; }
die(){ echo "[ERROR] $*" >&2; exit 1; }
need(){ [[ -f "$1" ]] || die "Missing file: $1"; }

# helper to build path
fpath() {
  local model="$1" scen="$2" file="$3"
  echo "${ROOT}/${model}/${scen}/${file}"
}

# =========================
# 1) BASELINE 1991–2020
#    historical (1991–2014) + ssp126 (2015–2020)
# =========================
log "Building baseline 1991–2020 = historical(1991–2014) + ssp126(2015–2020)"

BASE_FILES=()

for m in "${MODELS[@]}"; do
  H="$(fpath "${m}" "historical" "WBGT_mon_${m}_historical_1991-2014.nc")"
  S="$(fpath "${m}" "ssp126"     "WBGT_mon_${m}_ssp126_2015-2100.nc")"
  need "${H}"
  need "${S}"

  T1="$(mktemp --tmpdir="${TMP}" --suffix=.nc)"
  T2="$(mktemp --tmpdir="${TMP}" --suffix=.nc)"
  OUT_M="${TMP}/WBGT_mon_${m}_baseline_1991-2020.nc"

  cdo -L -O -selyear,1991/2014 "${H}" "${T1}"
  cdo -L -O -selyear,2015/2020 "${S}" "${T2}"
  cdo -L -O mergetime "${T1}" "${T2}" "${OUT_M}"

  rm -f "${T1}" "${T2}"
  BASE_FILES+=( "${OUT_M}" )
done

cdo -L -O ensmean "${BASE_FILES[@]}" "${OUT}/WBGT_mon_baseline_1991-2020_ensmean.nc"
cdo -L -O ensmin  "${BASE_FILES[@]}" "${OUT}/WBGT_mon_baseline_1991-2020_ensmin.nc"
cdo -L -O ensmax  "${BASE_FILES[@]}" "${OUT}/WBGT_mon_baseline_1991-2020_ensmax.nc"

log "Baseline ensembles done"

# =========================
# 2) FUTURE 2071–2100 (per SSP)
# =========================
for ssp in "${SSPS[@]}"; do
  log "Processing future ${ssp} (2071–2100)"

  FUT_FILES=()

  for m in "${MODELS[@]}"; do
    if [[ "${ssp}" == "ssp126" ]]; then
      F="$(fpath "${m}" "${ssp}" "WBGT_mon_${m}_${ssp}_2015-2100.nc")"
    else
      F="$(fpath "${m}" "${ssp}" "WBGT_mon_${m}_${ssp}_2020-2100.nc")"
    fi
    need "${F}"

    OUT_M="${TMP}/WBGT_mon_${m}_${ssp}_2071-2100.nc"
    cdo -L -O -selyear,2071/2100 "${F}" "${OUT_M}"
    FUT_FILES+=( "${OUT_M}" )
  done

  cdo -L -O ensmean "${FUT_FILES[@]}" "${OUT}/WBGT_mon_${ssp}_2071-2100_ensmean.nc"
  cdo -L -O ensmin  "${FUT_FILES[@]}" "${OUT}/WBGT_mon_${ssp}_2071-2100_ensmin.nc"
  cdo -L -O ensmax  "${FUT_FILES[@]}" "${OUT}/WBGT_mon_${ssp}_2071-2100_ensmax.nc"

done

log "All ensemble products written to ${OUT}"
