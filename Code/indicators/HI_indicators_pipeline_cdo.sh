#!/usr/bin/env bash
set -euo pipefail

############################################
# SPORT indicators pipeline (CDO) - NEX-GDDP-CMIP6
#
# Lit annuellement :
#   - tasmax (K)
#   - hurs   (%)
#   - pr     (kg m-2 s-1 ou mm/j selon NEX-GDDP)
#
# Produit 3 indicateurs, chacun dans SON NetCDF :
#
# 1) Heat Index (HI) → fichier HI_mon_MODEL_SCEN.nc
#    Variables :
#       ndHI32 : nb jours/mois HI > 32°C
#       ndHI38 : nb jours/mois HI > 38°C
#       ndHI40 : nb jours/mois HI > 40°C  (Danger)
#       ndHI46 : nb jours/mois HI > 46°C  (Extreme danger)
#
# 2) Pluie > 20 mm/j → fichier R20_mon_MODEL_SCEN.nc
#    Variable :
#       ndR20  : nb jours/mois pr > 20 mm/j
#
# 3) Vagues de chaleur → fichier HW_mon_MODEL_SCEN.nc
#    Définition :
#       vague = séquence ≥3 jours consécutifs avec tasmax_C > 35°C
#    Approx via runmean( hw_raw, 3j )
#    Variable :
#       ndHW  : nb jours/mois "en vague de chaleur"
#
# Entrées : fichiers journaliers par année, par scénario :
#   /.../MODEL/SCEN/tasmax/*.nc
#   /.../MODEL/SCEN/hurs/*.nc
#   /.../MODEL/SCEN/pr/*.nc
#
# Sorties intermédiaires : 1 fichier mensuel/an/indicateur :
#   OUT/SPORT/MODEL/SCEN/HI_mon_MODEL_SCEN_YYYY.nc
#   OUT/SPORT/MODEL/SCEN/R20_mon_MODEL_SCEN_YYYY.nc
#   OUT/SPORT/MODEL/SCEN/HW_mon_MODEL_SCEN_YYYY.nc
#
# Sorties finales par scénario (toute la période dispo) :
#   HI_mon_MODEL_SCEN.nc
#   R20_mon_MODEL_SCEN.nc
#   HW_mon_MODEL_SCEN.nc
#
# Post-traitement :
#   - Histo 1991–2020 (historical + ssp126) :
#       <OUT_ROOT>/MODEL/HI_mon_MODEL_HISTO_1991-2020.nc
#       <OUT_ROOT>/MODEL/R20_mon_MODEL_HISTO_1991-2020.nc
#       <OUT_ROOT>/MODEL/HW_mon_MODEL_HISTO_1991-2020.nc
#
#   - Futur 2020–2100 pour ssp126/245/585 :
#       <OUT_ROOT>/MODEL/HI_mon_MODEL_sspXXX_2020-2100.nc
#       idem pour R20 / HW
#
# USAGE :
#   ./sport_indicators_pipeline_cdo.sh MODEL SCEN1 [SCEN2 ...]
#
# EXEMPLE :
#   ./sport_indicators_HI_pipeline_cdo.sh GFDL-ESM4 historical ssp126 ssp245 ssp585
############################################

IN_ROOT="/mnt/e/CMIP6/NEX-GDDP"
OUT_ROOT="/mnt/c/Data/SPORT"

if [ "$#" -lt 2 ]; then
  echo "Usage: $0 MODEL SCEN1 [SCEN2 ...]" >&2
  exit 1
fi

MODEL="$1"
shift
SCENS=("$@")

shopt -s nullglob

echo "==== SPORT indicators pipeline for model: ${MODEL} ===="
echo "Scenarios: ${SCENS[*]}"
echo

for SCEN in "${SCENS[@]}"; do
  echo "---- Scenario: ${SCEN} ----"

  IN_DIR="${IN_ROOT}/${MODEL}/${SCEN}"
  if [ ! -d "${IN_DIR}" ]; then
    echo "  [WARN] Input dir not found: ${IN_DIR}  -> skipping scenario" >&2
    continue
  fi

  OUT_DIR="${OUT_ROOT}/${MODEL}/${SCEN}"
  mkdir -p "${OUT_DIR}"

  YEARLY_HI_FILES=()
  
  # Boucle sur les fichiers tasmax annuels
  for TAS_FILE in "${IN_DIR}/tasmax/"*.nc; do
    BASENAME="$(basename "${TAS_FILE}")"

    BASENAME_HURS="${BASENAME/tasmax_day/hurs_day}"
    HURS_FILE="${IN_DIR}/hurs/${BASENAME_HURS}"

  
    if [ ! -f "${HURS_FILE}" ]; then
      echo "  [WARN] Missing hurs file for: ${TAS_FILE} -> expected ${HURS_FILE}, skipping year" >&2
      continue
    fi

    YEAR="$(echo "${BASENAME}" | grep -oE '[0-9]{4}' | tail -n1)"
    if [ -z "${YEAR}" ]; then
      echo "  [WARN] Could not extract year from filename: ${BASENAME}, skipping" >&2
      continue
    fi

    OUT_HI_YEAR="${OUT_DIR}/HI_mon_${MODEL}_${SCEN}_${YEAR}.nc"
   
    YEARLY_HI_FILES+=("${OUT_HI_YEAR}")
    
    # Si tout existe déjà, on saute
    if [ -f "${OUT_HI_YEAR}" ]; then
      echo "  [INFO] Monthly SPORT files already exist for ${YEAR}, skipping computation."
      continue
    fi

 echo "  * Processing year ${YEAR}..."

    TMP_DIR="$(mktemp -d)"
    DAILY_MERGE="${TMP_DIR}/daily_merge_${YEAR}.nc"
    DAILY_FLAGS="${TMP_DIR}/daily_flags_${YEAR}.nc"
    MON_HI="${TMP_DIR}/mon_hi_${YEAR}.nc"

    ############################################
    # 1) Merge quotidien tasmax + hurs
    ############################################
    cdo -L -z zip_5 merge "${TAS_FILE}" "${HURS_FILE}" "${DAILY_MERGE}"

    ############################################
    # 2) Expr quotidien : HI (corrigé NOAA)
    ############################################
cdo -L -z zip_5 \
  -expr,"\
    tas_c = tasmax-273.15; \
    rh    = max(min(hurs,99),1); \
    T_F   = tas_c*9.0/5.0 + 32.0; \
    HI_F_raw  = -42.379 \
                + 2.04901523*T_F \
                + 10.14333127*rh \
                - 0.22475541*T_F*rh \
                - 0.00683783*T_F*T_F \
                - 0.05481717*rh*rh \
                + 0.00122874*T_F*T_F*rh \
                + 0.00085282*T_F*rh*rh \
                - 0.00000199*T_F*T_F*rh*rh; \
    HI_F = (T_F >= 80.0 && rh >= 40.0) ? HI_F_raw : T_F; \
    HI_C = (HI_F - 32.0)*5.0/9.0; \
    ndHI32 = HI_C > 32.0; \
    ndHI38 = HI_C > 38.0; \
    ndHI40 = HI_C > 40.0; \
    ndHI46 = HI_C > 46.0; \
  " \
  "${DAILY_MERGE}" \
  "${DAILY_FLAGS}"


    ############################################
    # 3) Agrégation mensuelle HI (monsum)
    ############################################
    cdo -L -z zip_5 \
      -monsum \
      -select,name=ndHI32,ndHI38,ndHI40,ndHI46 \
      "${DAILY_FLAGS}" \
      "${MON_HI}"

    ############################################
    # 4) Écriture du NetCDF annuel HI
    ############################################
    cdo -L -z zip_5 -select,name=ndHI32,ndHI38,ndHI40,ndHI46 \
      "${MON_HI}" "${OUT_HI_YEAR}"

    rm -rf "${TMP_DIR}"
    

  done

  # Si aucun fichier, on saute le merge
  if [ "${#YEARLY_HI_FILES[@]}" -eq 0 ]; then
    echo "  [WARN] No yearly SPORT files for ${MODEL}/${SCEN}, skipping scenario merge." >&2
    continue
  fi

  ############################################
  # 6) Merge temporel par indicateur (scénario)
  ############################################
  HI_SCEN_OUT="${OUT_DIR}/HI_mon_${MODEL}_${SCEN}.nc"
 
  echo "  -> Merging yearly HI files into: ${HI_SCEN_OUT}"
  cdo -L -z zip_5 mergetime "${YEARLY_HI_FILES[@]}" "${HI_SCEN_OUT}"

  echo "  -> Removing yearly monthly intermediates..."
  rm -f "${YEARLY_HI_FILES[@]}" 

  echo "  [OK] Scenario ${SCEN} done."
  echo

done

############################################
# 7) Post-traitement : HISTO 1991–2020 et FUTUR 2020–2100
############################################

echo "---- Post-processing HISTO and FUTURE files ----"

INDICATORS=("HI")

# Histo : historical + ssp126 -> 1991-2020
for IND in "${INDICATORS[@]}"; do
  HIST_FILE="${OUT_ROOT}/${MODEL}/historical/${IND}_mon_${MODEL}_historical.nc"
  SSP126_FILE="${OUT_ROOT}/${MODEL}/ssp126/${IND}_mon_${MODEL}_ssp126.nc"

  if [ -f "${HIST_FILE}" ] && [ -f "${SSP126_FILE}" ]; then
    TMP_MERGE="${OUT_ROOT}/${MODEL}/${IND}_mon_${MODEL}_HISTO_tmp.nc"
    HISTO_OUT="${OUT_ROOT}/${MODEL}/${IND}_mon_${MODEL}_HISTO_1991-2020.nc"

    echo "  -> Building HISTO 1991-2020 for ${IND} (historical + ssp126)"
    cdo -L -z zip_5 mergetime "${HIST_FILE}" "${SSP126_FILE}" "${TMP_MERGE}"
    cdo -L -z zip_5 selyear,1991/2020 "${TMP_MERGE}" "${HISTO_OUT}"
    rm -f "${TMP_MERGE}"
  else
    echo "  [WARN] Missing HIST or ssp126 for ${IND}, skipping HISTO build." >&2
  fi
done

# Futur 2020–2100 pour chaque scénario
for IND in "${INDICATORS[@]}"; do
  for SCEN in ssp126 ssp245 ssp585; do
    SCEN_FILE="${OUT_ROOT}/${MODEL}/${SCEN}/${IND}_mon_${MODEL}_${SCEN}.nc"
    if [ -f "${SCEN_FILE}" ]; then
      FUT_OUT="${OUT_ROOT}/${MODEL}/${IND}_mon_${MODEL}_${SCEN}_2020-2100.nc"
      echo "  -> Building FUTURE 2020-2100 for ${IND} / ${SCEN}"
      cdo -L -z zip_5 selyear,2020/2100 "${SCEN_FILE}" "${FUT_OUT}"
    fi
  done
done

echo "==== All done for model: ${MODEL} ===="

