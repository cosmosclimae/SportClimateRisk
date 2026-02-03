#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

############################################
# WBGT pipeline (CDO) - NEX-GDDP-CMIP6
#
# Inputs (daily): tasmin, tasmax, tas (K), hurs (%), rsds (W m-2)
# - RH diurne via e constant from RHmean & Tmean
# - Tw via Stull (2011)
# - Tg proxy: Tg = 0.01498*RSDS_Jm2 + 1.184*T - 0.0789*RH - 2.739
# - WBGT = 0.7*Tw + 0.2*Tg + 0.1*T
# - Outputs (monthly):
#   * ndays_wbgt_{morn|aft|eve}_<bin> : days/month in each WBGT class
#   * wbgt_{morn|aft|eve}_mean : monthly mean WBGT (degC)
############################################

IN_ROOT="/mnt/e/CMIP6/NEX-GDDP"
OUT_ROOT="/mnt/c/Data/WBGT"

if [ "$#" -lt 2 ]; then
  echo "Usage: $0 MODEL SCEN1 [SCEN2 ...]" >&2
  exit 1
fi

MODEL="$1"; shift
SCENS=("$@")

echo "==== WBGT pipeline for model: ${MODEL} ===="
echo "Scenarios: ${SCENS[*]}"
echo

for SCEN in "${SCENS[@]}"; do
  echo "---- Scenario: ${SCEN} ----"

  IN_DIR="${IN_ROOT}/${MODEL}/${SCEN}"
  [ -d "${IN_DIR}" ] || { echo "  [WARN] Missing ${IN_DIR} -> skip" >&2; continue; }

  OUT_DIR="${OUT_ROOT}/${MODEL}/${SCEN}"
  mkdir -p "${OUT_DIR}"

  YEARLY_MON_FILES=()

  for TASMAX_FILE in "${IN_DIR}/tasmax/"*.nc; do
    BASENAME="$(basename "${TASMAX_FILE}")"

    BASENAME_TASMIN="${BASENAME/tasmax_day/tasmin_day}"
    BASENAME_TAS="${BASENAME/tasmax_day/tas_day}"
    BASENAME_HURS="${BASENAME/tasmax_day/hurs_day}"
    BASENAME_RSDS="${BASENAME/tasmax_day/rsds_day}"

    TASMIN_FILE="${IN_DIR}/tasmin/${BASENAME_TASMIN}"
    TAS_FILE="${IN_DIR}/tas/${BASENAME_TAS}"
    HURS_FILE="${IN_DIR}/hurs/${BASENAME_HURS}"
    RSDS_FILE="${IN_DIR}/rsds/${BASENAME_RSDS}"

    if [ ! -f "${TASMIN_FILE}" ] || [ ! -f "${TAS_FILE}" ] || [ ! -f "${HURS_FILE}" ] || [ ! -f "${RSDS_FILE}" ]; then
      echo "  [WARN] Missing companion file(s) for ${BASENAME} -> skip year" >&2
      continue
    fi

    YEAR="$(echo "${BASENAME}" | grep -oE '[0-9]{4}' | tail -n1)"
    [ -n "${YEAR}" ] || { echo "  [WARN] No year in ${BASENAME} -> skip" >&2; continue; }

    OUT_MON_YEAR="${OUT_DIR}/wbgt_mon_${MODEL}_${SCEN}_${YEAR}.nc"
    YEARLY_MON_FILES+=("${OUT_MON_YEAR}")

    if [ -f "${OUT_MON_YEAR}" ]; then
      echo "  [INFO] Exists ${OUT_MON_YEAR} -> skip"
      continue
    fi

    echo "  * Processing year ${YEAR}..."

    TMP_DAY="$(mktemp --suffix=.nc)"
    TMP_MON="$(mktemp --suffix=.nc)"

    # 1) construire WBGT journaliers + flags journaliers (sortie daily)
    cdo -L -z zip_5 \
      -expr,"\
        # Kelvin -> Celsius
        tmin = tasmin - 273.15; \
        tmax = tasmax - 273.15; \
        tmean = tas - 273.15; \
        rhm = max(min(hurs,99),1); \
        rsds_j = rsds ; \
        # Saturation vapor pressure (hPa), Tetens
        es_tmean = 6.112*exp((17.67*tmean)/(tmean+243.5)); \
        e = (rhm/100.) * es_tmean; \
        es_tmin = 6.112*exp((17.67*tmin)/(tmin+243.5)); \
        es_tmax = 6.112*exp((17.67*tmax)/(tmax+243.5)); \
        rh_morn = max(min(100.*e/es_tmin,99),1); \
        rh_aft  = max(min(100.*e/es_tmax,99),1); \
        rh_eve  = max(min(100.*e/es_tmean,99),1); \

        # Tw Stull (2011)
        tw_morn = tmin*atan(0.151977*sqrt(rh_morn+8.313659)) + atan(tmin+rh_morn) - atan(rh_morn-1.676331) + 0.00391838*pow(rh_morn,1.5)*atan(0.023101*rh_morn) - 4.686035; \
        tw_aft  = tmax*atan(0.151977*sqrt(rh_aft +8.313659)) + atan(tmax+rh_aft ) - atan(rh_aft -1.676331) + 0.00391838*pow(rh_aft ,1.5)*atan(0.023101*rh_aft ) - 4.686035; \
        tw_eve  = tmean*atan(0.151977*sqrt(rh_eve +8.313659)) + atan(tmean+rh_eve) - atan(rh_eve-1.676331) + 0.00391838*pow(rh_eve,1.5)*atan(0.023101*rh_eve) - 4.686035; \

        # Tg proxy (RSDS=0 morn/eve)
        tg_morn = 1.184*tmin  - 0.0789*rh_morn - 2.739; \
        tg_aft  = 0.01498*rsds_j + 1.184*tmax - 0.0789*rh_aft  - 2.739; \
        tg_eve  = 1.184*tmean - 0.0789*rh_eve  - 2.739; \

        # WBGT
        wbgt_morn = 0.7*tw_morn + 0.2*tg_morn + 0.1*tmin; \
        wbgt_aft  = 0.7*tw_aft  + 0.2*tg_aft  + 0.1*tmax; \
        wbgt_eve  = 0.7*tw_eve  + 0.2*tg_eve  + 0.1*tmean; \

        # Flags journaliers par classes
        nd_morn_c0 = wbgt_morn < 18; \
        nd_morn_c1 = (wbgt_morn >= 18) * (wbgt_morn < 24); \
        nd_morn_c2 = (wbgt_morn >= 24) * (wbgt_morn < 26); \
        nd_morn_c3 = (wbgt_morn >= 26) * (wbgt_morn < 28); \
        nd_morn_c4 = (wbgt_morn >= 28) * (wbgt_morn < 30); \
        nd_morn_c5 = (wbgt_morn >= 30) * (wbgt_morn < 32); \
        nd_morn_c6 = wbgt_morn >= 32; \

        nd_aft_c0 = wbgt_aft < 18; \
        nd_aft_c1 = (wbgt_aft >= 18) * (wbgt_aft < 24); \
        nd_aft_c2 = (wbgt_aft >= 24) * (wbgt_aft < 26); \
        nd_aft_c3 = (wbgt_aft >= 26) * (wbgt_aft < 28); \
        nd_aft_c4 = (wbgt_aft >= 28) * (wbgt_aft < 30); \
        nd_aft_c5 = (wbgt_aft >= 30) * (wbgt_aft < 32); \
        nd_aft_c6 = wbgt_aft >= 32; \

        nd_eve_c0 = wbgt_eve < 18; \
        nd_eve_c1 = (wbgt_eve >= 18) * (wbgt_eve < 24); \
        nd_eve_c2 = (wbgt_eve >= 24) * (wbgt_eve < 26); \
        nd_eve_c3 = (wbgt_eve >= 26) * (wbgt_eve < 28); \
        nd_eve_c4 = (wbgt_eve >= 28) * (wbgt_eve < 30); \
        nd_eve_c5 = (wbgt_eve >= 30) * (wbgt_eve < 32); \
        nd_eve_c6 = wbgt_eve >= 32; \
      " \
      -merge "${TASMAX_FILE}" "${TASMIN_FILE}" "${TAS_FILE}" "${HURS_FILE}" "${RSDS_FILE}" \
      "${TMP_DAY}"

    # 2) agréger en mensuel :
    # - ndays = monsum(flags)
    # - mean WBGT = monmean(wbgt_*)
    cdo -L -z zip_5 \
      -merge \
        -monsum -selname,nd_morn_c0,nd_morn_c1,nd_morn_c2,nd_morn_c3,nd_morn_c4,nd_morn_c5,nd_morn_c6,\
nd_aft_c0,nd_aft_c1,nd_aft_c2,nd_aft_c3,nd_aft_c4,nd_aft_c5,nd_aft_c6,\
nd_eve_c0,nd_eve_c1,nd_eve_c2,nd_eve_c3,nd_eve_c4,nd_eve_c5,nd_eve_c6 "${TMP_DAY}" \
        -monmean -selname,wbgt_morn,wbgt_aft,wbgt_eve "${TMP_DAY}" \
      "${TMP_MON}"

    # 3) renommer variables (lisibles)
    cdo -L -z zip_5 \
      -chname,nd_morn_c0,ndays_wbgt_morn_lt18,\
nd_morn_c1,ndays_wbgt_morn_18_24,\
nd_morn_c2,ndays_wbgt_morn_24_26,\
nd_morn_c3,ndays_wbgt_morn_26_28,\
nd_morn_c4,ndays_wbgt_morn_28_30,\
nd_morn_c5,ndays_wbgt_morn_30_32,\
nd_morn_c6,ndays_wbgt_morn_ge32,\
nd_aft_c0,ndays_wbgt_aft_lt18,\
nd_aft_c1,ndays_wbgt_aft_18_24,\
nd_aft_c2,ndays_wbgt_aft_24_26,\
nd_aft_c3,ndays_wbgt_aft_26_28,\
nd_aft_c4,ndays_wbgt_aft_28_30,\
nd_aft_c5,ndays_wbgt_aft_30_32,\
nd_aft_c6,ndays_wbgt_aft_ge32,\
nd_eve_c0,ndays_wbgt_eve_lt18,\
nd_eve_c1,ndays_wbgt_eve_18_24,\
nd_eve_c2,ndays_wbgt_eve_24_26,\
nd_eve_c3,ndays_wbgt_eve_26_28,\
nd_eve_c4,ndays_wbgt_eve_28_30,\
nd_eve_c5,ndays_wbgt_eve_30_32,\
nd_eve_c6,ndays_wbgt_eve_ge32,\
wbgt_morn,wbgt_morn_mean,\
wbgt_aft,wbgt_aft_mean,\
wbgt_eve,wbgt_eve_mean \
      "${TMP_MON}" "${OUT_MON_YEAR}"

    # 4) attributs (units + long_name)
    cdo -L -z zip_5 \
      -setattribute,wbgt_morn_mean@units="degC" \
      -setattribute,wbgt_aft_mean@units="degC" \
      -setattribute,wbgt_eve_mean@units="degC" \
      -setattribute,wbgt_morn_mean@long_name="Monthly mean Wet Bulb Globe Temperature (morning proxy)" \
      -setattribute,wbgt_aft_mean@long_name="Monthly mean Wet Bulb Globe Temperature (afternoon proxy)" \
      -setattribute,wbgt_eve_mean@long_name="Monthly mean Wet Bulb Globe Temperature (evening proxy)" \
      \
      -setattribute,ndays_wbgt_morn_lt18@units="days" \
      -setattribute,ndays_wbgt_morn_18_24@units="days" \
      -setattribute,ndays_wbgt_morn_24_26@units="days" \
      -setattribute,ndays_wbgt_morn_26_28@units="days" \
      -setattribute,ndays_wbgt_morn_28_30@units="days" \
      -setattribute,ndays_wbgt_morn_30_32@units="days" \
      -setattribute,ndays_wbgt_morn_ge32@units="days" \
      -setattribute,ndays_wbgt_aft_lt18@units="days" \
      -setattribute,ndays_wbgt_aft_18_24@units="days" \
      -setattribute,ndays_wbgt_aft_24_26@units="days" \
      -setattribute,ndays_wbgt_aft_26_28@units="days" \
      -setattribute,ndays_wbgt_aft_28_30@units="days" \
      -setattribute,ndays_wbgt_aft_30_32@units="days" \
      -setattribute,ndays_wbgt_aft_ge32@units="days" \
      -setattribute,ndays_wbgt_eve_lt18@units="days" \
      -setattribute,ndays_wbgt_eve_18_24@units="days" \
      -setattribute,ndays_wbgt_eve_24_26@units="days" \
      -setattribute,ndays_wbgt_eve_26_28@units="days" \
      -setattribute,ndays_wbgt_eve_28_30@units="days" \
      -setattribute,ndays_wbgt_eve_30_32@units="days" \
      -setattribute,ndays_wbgt_eve_ge32@units="days" \
      \
      "${OUT_MON_YEAR}" "${OUT_MON_YEAR}.tmp" && mv "${OUT_MON_YEAR}.tmp" "${OUT_MON_YEAR}"

    rm -f "${TMP_DAY}" "${TMP_MON}"
  done

  if [ "${#YEARLY_MON_FILES[@]}" -eq 0 ]; then
    echo "  [WARN] No yearly files for ${MODEL}/${SCEN} -> skip merge" >&2
    continue
  fi

  start_year=$(printf '%s\n' "${YEARLY_MON_FILES[@]}" | grep -oE '[0-9]{4}' | sort | head -n1)
  end_year=$(printf '%s\n' "${YEARLY_MON_FILES[@]}" | grep -oE '[0-9]{4}' | sort | tail -n1)

  FINAL_OUT="${OUT_DIR}/WBGT_mon_${MODEL}_${SCEN}_${start_year}-${end_year}.nc"
  echo "  -> Merging into: ${FINAL_OUT}"
  cdo -L -z zip_5 mergetime "${YEARLY_MON_FILES[@]}" "${FINAL_OUT}"

  echo "  -> Removing yearly intermediates..."
  rm -f "${YEARLY_MON_FILES[@]}"

  echo "  [OK] Scenario ${SCEN} done."
  echo
done

echo "==== All done for model: ${MODEL} ===="
