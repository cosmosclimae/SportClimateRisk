## =========================
## 0. LIBRAIRIES & PATHS
## =========================
library(terra)
library(sf)
library(ggplot2)
library(viridis)
library(patchwork)
library(grid)      # pour unit()

## Dossiers à adapter
wbgt_dir  <- "C:/Data/WBGT/ENSEMBLE"   # WBGT (multi-subds)
sport_dir <- "C:/Data/SPORT/ensemble"  # R20 / HW / HI
shape_dir <- "C:/Data/SHAPE/NE110"

## =========================
## 1. SHAPEFILES & THÈME COMMUN
## =========================
coast   <- st_read(file.path(shape_dir, "ne_110m_coastline.shp"))
land    <- st_read(file.path(shape_dir, "ne_110m_land.shp"))

# Projection Robinson
crs_rob <- "ESRI:54030"
coast_r <- st_transform(coast, crs_rob)
land_r  <- st_transform(land,  crs_rob)

theme_map <- theme_void() +
  theme(
    plot.title       = element_text(size = 10, face = "bold", hjust = 0.5),
    legend.position  = "bottom",
    legend.key.width = unit(1.4, "cm")
  )

## =========================
## 2. HELPERS GÉNÉRIQUES
## =========================

## 2.1 Raster -> saison (fichiers multi-variables via subds)
r_season <- function(file, var, season = c("JJA","DJF")) {
  season <- match.arg(season)
  r <- rast(file, subds = var)
  n <- nlyr(r)
  nyear <- n / 12
  if (n %% 12 != 0) stop("n not multiple de 12 in ", file, " (var=", var, ")")

  months <- rep(1:12, length.out = n)
  idx <- if (season == "JJA") {
    which(months %in% c(6,7,8))
  } else {
    which(months %in% c(12,1,2))
  }

  r_mean <- sum(r[[idx]]) / nyear   # nb moyen de jours / saison
  r_mean[r_mean == 0] <- NA         # masque les zéros
  project(r_mean, crs_rob)
}

## 2.2 Raster -> saison (une seule variable par fichier)
r_season_1var <- function(file, season = c("JJA","DJF")) {
  season <- match.arg(season)
  r <- rast(file)     # 1 seule variable dans le NetCDF
  n <- nlyr(r)
  nyear <- n / 12
  if (n %% 12 != 0) stop("n not multiple de 12 in ", file)

  months <- rep(1:12, length.out = n)
  idx <- if (season == "JJA") {
    which(months %in% c(6,7,8))
  } else {
    which(months %in% c(12,1,2))
  }

  r_mean <- sum(r[[idx]]) / nyear
  r_mean[r_mean == 0] <- NA
  project(r_mean, crs_rob)
}

## 2.3 Raster -> data.frame pour ggplot
to_df <- function(r, nm = "value") {
  as.data.frame(r, xy = TRUE, na.rm = TRUE) |>
    setNames(c("x", "y", nm))
}

## 2.4 Panel "jours par saison" (viridis)
make_panel_days <- function(df, title, vmax, legend_title,
                            option = "plasma") {
  breaks_vec <- c(0, 5, 10, 20, 30, 40, 50, 60, 80)
  breaks_vec <- breaks_vec[breaks_vec <= vmax]

  ggplot() +
    geom_sf(data = land_r, fill = "#242222", color = NA) +
    geom_raster(data = df, aes(x = x, y = y, fill = days)) +
    geom_sf(data = coast_r, color = "black", size = 0.1) +
    scale_fill_viridis(
      option = option,
      limits = c(0, vmax),
      breaks = breaks_vec,
      name   = legend_title
    ) +
    coord_sf(crs = crs_rob, expand = FALSE) +
    theme_map +
    labs(title = title)
}

## 2.5 Panel delta (gradient biaxe bleu-blanc-rouge)
make_panel_delta <- function(df, title, limits_days) {
  ggplot() +
    geom_sf(data = land_r, fill = "#242222", color = NA) +
    geom_raster(data = df, aes(x = x, y = y, fill = ddays)) +
    geom_sf(data = coast_r, color = "black", size = 0.1) +
    scale_fill_gradient2(
      low      = "#2166ac",
      mid      = "white",
      high     = "#b2182b",
      midpoint = 0,
      limits   = c(-limits_days, limits_days),
      name     = "Δ days (PM – AM)"
    ) +
    coord_sf(crs = crs_rob, expand = FALSE) +
    theme_map +
    labs(title = title)
}

## 2.6 Normalisation min-max globale (toutes couches)
normalize_indicator <- function(r) {
  if (!inherits(r, "SpatRaster")) {
    stop("normalize_indicator attend un terra::SpatRaster.")
  }
  
  rg <- terra::global(r, fun = "range", na.rm = TRUE)
  xmin <- min(rg[, 1], na.rm = TRUE)
  xmax <- max(rg[, 2], na.rm = TRUE)
  
  if (!is.finite(xmin) || !is.finite(xmax) || xmax <= xmin) {
    warning("Impossible de normaliser : xmax <= xmin ou valeurs non finies. Retourne NA.")
    return(r * NA)
  }
  
  (r - xmin) / (xmax - xmin)
}

## 2.7 Vérification de géométrie
check_geom_same <- function(...) {
  rast_list <- list(...)
  ref <- rast_list[[1]]
  for (i in seq_along(rast_list)[-1]) {
    if (!terra::compareGeom(ref, rast_list[[i]], stopOnError = FALSE)) {
      stop("Les rasters n'ont pas la même géométrie (extent / résolution / CRS).")
    }
  }
}

## 2.8 Construction SCRI
build_scri <- function(r_wbgt_aft_n, r_hi_n, r_rain_n, r_hw_n,
                       w_wbgt = 0.25, w_hi = 0.25,
                       w_rain = 0.25, w_hw = 0.25) {
  
  check_geom_same(r_wbgt_aft_n, r_hi_n, r_rain_n, r_hw_n)
  w_sum <- w_wbgt + w_hi + w_rain + w_hw
  
  scri <- (w_wbgt * r_wbgt_aft_n +
           w_hi   * r_hi_n +
           w_rain * r_rain_n +
           w_hw   * r_hw_n)
  
  # si tu veux forcer 0-1 même si w_sum != 1 :
  # scri <- scri / w_sum
  
  scri
}

## =========================
## 3. WBGT : FIGURE 2 & 3 (3×3 JJA / DJF)
## =========================

## fichiers WBGT (à adapter)
hist_f <- file.path(wbgt_dir, "WBGT_mon_baseline_1991-2020_ensmean.nc")
mid_f  <- file.path(wbgt_dir, "WBGT_mon_ssp585_2031-2060_ensmean.nc")
late_f <- file.path(wbgt_dir, "WBGT_mon_ssp585_2071-2100_ensmean.nc")

## mapping des subds WBGT
vars <- c(
  Morning   = "ndays_wbgt_morn_ge32",
  Afternoon = "ndays_wbgt_aft_ge32",
  Evening   = "ndays_wbgt_eve_ge32"
)

make_panel_wbgt <- function(df, title, limits_days) {
  make_panel_days(
    df           = df,
    title        = title,
    vmax         = limits_days,
    legend_title = "Days per season (WBGT ≥ 32°C)",
    option       = "plasma"
  )
}

make_fig_wbgt <- function(season = "JJA") {
  files3      <- c(hist_f, mid_f, late_f)
  labels_rows <- c("Morning", "Afternoon", "Evening")
  labels_cols <- c("1991–2020", "2050 (2031–2060)", "2100 (2071–2100)")

  ## 1) échelle commune
  layers <- list()
  for (v in vars) {
    for (f in files3) {
      layers[[length(layers) + 1]] <- r_season(f, v, season)
    }
  }
  all_vals <- unlist(lapply(layers, values))
  all_vals <- all_vals[is.finite(all_vals)]
  vmax <- min(80, ceiling(max(all_vals, na.rm = TRUE)))

  ## 2) construire les 9 panels
  panels <- list()
  i <- 1
  row_id <- 1

  for (vname in names(vars)) {
    var <- vars[[vname]]
    col_id <- 1

    for (f in files3) {
      r  <- r_season(f, var, season)
      df <- to_df(r, "days")

      title <- paste0(labels_rows[row_id], " – ", labels_cols[col_id])
      panels[[i]] <- make_panel_wbgt(df, title, vmax)

      i      <- i + 1
      col_id <- col_id + 1
    }
    row_id <- row_id + 1
  }

  ## 3) assembler 3×3 + légende unique
  fig <- (panels[[1]] | panels[[2]] | panels[[3]]) /
         (panels[[4]] | panels[[5]] | panels[[6]]) /
         (panels[[7]] | panels[[8]] | panels[[9]])

  fig + plot_layout(guides = "collect") & theme(legend.position = "bottom")
}

## Figure 2 : JJA
fig2 <- make_fig_wbgt("JJA")
ggsave("C:/Data/WBGT/Figure2_WBGT_JJA_3x3.png",
       fig2, width = 12, height = 11, dpi = 300)

## Figure 3 : DJF
fig3 <- make_fig_wbgt("DJF")
ggsave("C:/Data/WBGT/Figure3_WBGT_DJF_3x3.png",
       fig3, width = 12, height = 11, dpi = 300)

## =========================
## 4. WBGT Δ (PM – AM) : FIGURE 4 (2×3)
## =========================

## JJA
r_hist_jja_aft  <- r_season(hist_f, vars[["Afternoon"]], "JJA")
r_hist_jja_morn <- r_season(hist_f, vars[["Morning"]],   "JJA")

r_mid_jja_aft   <- r_season(mid_f,  vars[["Afternoon"]], "JJA")
r_mid_jja_morn  <- r_season(mid_f,  vars[["Morning"]],   "JJA")

r_late_jja_aft  <- r_season(late_f, vars[["Afternoon"]], "JJA")
r_late_jja_morn <- r_season(late_f, vars[["Morning"]],   "JJA")

r_hist_jja_delta <- r_hist_jja_aft - r_hist_jja_morn
r_mid_jja_delta  <- r_mid_jja_aft  - r_mid_jja_morn
r_late_jja_delta <- r_late_jja_aft - r_late_jja_morn

## DJF
r_hist_djf_aft  <- r_season(hist_f, vars[["Afternoon"]], "DJF")
r_hist_djf_morn <- r_season(hist_f, vars[["Morning"]],   "DJF")

r_mid_djf_aft   <- r_season(mid_f,  vars[["Afternoon"]], "DJF")
r_mid_djf_morn  <- r_season(mid_f,  vars[["Morning"]],   "DJF")

r_late_djf_aft  <- r_season(late_f, vars[["Afternoon"]], "DJF")
r_late_djf_morn <- r_season(late_f, vars[["Morning"]],   "DJF")

r_hist_djf_delta <- r_hist_djf_aft - r_hist_djf_morn
r_mid_djf_delta  <- r_mid_djf_aft  - r_mid_djf_morn
r_late_djf_delta <- r_late_djf_aft - r_late_djf_morn

## échelle commune symétrique
vals_all_delta <- c(
  values(r_hist_jja_delta), values(r_mid_jja_delta), values(r_late_jja_delta),
  values(r_hist_djf_delta), values(r_mid_djf_delta), values(r_late_djf_delta)
)
vals_all_delta <- vals_all_delta[is.finite(vals_all_delta)]
vmax_delta <- ceiling(max(abs(vals_all_delta), na.rm = TRUE))

## data.frames
df_hist_jja <- to_df(r_hist_jja_delta, "ddays")
df_mid_jja  <- to_df(r_mid_jja_delta,  "ddays")
df_late_jja <- to_df(r_late_jja_delta, "ddays")

df_hist_djf <- to_df(r_hist_djf_delta, "ddays")
df_mid_djf  <- to_df(r_mid_djf_delta,  "ddays")
df_late_djf <- to_df(r_late_djf_delta, "ddays")

## panels
p_hist_jja <- make_panel_delta(df_hist_jja, "(a) JJA, 1991–2020",            vmax_delta)
p_mid_jja  <- make_panel_delta(df_mid_jja,  "(b) JJA, 2031–2060 (SSP5–8.5)", vmax_delta)
p_late_jja <- make_panel_delta(df_late_jja, "(c) JJA, 2071–2100 (SSP5–8.5)", vmax_delta)

p_hist_djf <- make_panel_delta(df_hist_djf, "(d) DJF, 1991–2020",            vmax_delta)
p_mid_djf  <- make_panel_delta(df_mid_djf,  "(e) DJF, 2031–2060 (SSP5–8.5)", vmax_delta)
p_late_djf <- make_panel_delta(df_late_djf, "(f) DJF, 2071–2100 (SSP5–8.5)", vmax_delta)

fig4 <- (p_hist_jja | p_mid_jja | p_late_jja) /
        (p_hist_djf | p_mid_djf | p_late_djf)

fig4 <- fig4 + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("C:/Data/WBGT/Figure4_WBGT_delta_PM_AM.png",
       fig4, width = 12, height = 7.5, dpi = 300)

## =========================
## 5. R20 & HEATWAVES : FIGURES 6 & 7
## =========================

make_fig_indicator <- function(prefix, legend_title, outfile_png) {
  f_hist <- file.path(sport_dir, paste0(prefix, "_HISTO_1991-2020.nc"))
  f_mid  <- file.path(sport_dir, paste0(prefix, "_ssp585_2031-2060.nc"))
  f_late <- file.path(sport_dir, paste0(prefix, "_ssp585_2071-2100.nc"))

  ## extraire toutes les saisons pour échelle commune
  r_hist_jja <- r_season_1var(f_hist, "JJA")
  r_mid_jja  <- r_season_1var(f_mid,  "JJA")
  r_late_jja <- r_season_1var(f_late, "JJA")

  r_hist_djf <- r_season_1var(f_hist, "DJF")
  r_mid_djf  <- r_season_1var(f_mid,  "DJF")
  r_late_djf <- r_season_1var(f_late, "DJF")

  vals_all <- c(
    values(r_hist_jja), values(r_mid_jja), values(r_late_jja),
    values(r_hist_djf), values(r_mid_djf), values(r_late_djf)
  )
  vals_all <- vals_all[is.finite(vals_all)]
  vmax <- min(80, ceiling(max(vals_all, na.rm = TRUE)))

  ## data.frames
  df_hist_jja <- to_df(r_hist_jja, "days")
  df_mid_jja  <- to_df(r_mid_jja,  "days")
  df_late_jja <- to_df(r_late_jja, "days")

  df_hist_djf <- to_df(r_hist_djf, "days")
  df_mid_djf  <- to_df(r_mid_djf,  "days")
  df_late_djf <- to_df(r_late_djf, "days")

  ## panels
  p_hist_jja <- make_panel_days(df_hist_jja, "(a) JJA, 1991–2020",            vmax, legend_title)
  p_mid_jja  <- make_panel_days(df_mid_jja,  "(b) JJA, 2031–2060 (SSP5–8.5)", vmax, legend_title)
  p_late_jja <- make_panel_days(df_late_jja, "(c) JJA, 2071–2100 (SSP5–8.5)", vmax, legend_title)

  p_hist_djf <- make_panel_days(df_hist_djf, "(d) DJF, 1991–2020",            vmax, legend_title)
  p_mid_djf  <- make_panel_days(df_mid_djf,  "(e) DJF, 2031–2060 (SSP5–8.5)", vmax, legend_title)
  p_late_djf <- make_panel_days(df_late_djf, "(f) DJF, 2071–2100 (SSP5–8.5)", vmax, legend_title)

  fig <- (p_hist_jja | p_mid_jja | p_late_jja) /
         (p_hist_djf | p_mid_djf | p_late_djf)

  fig <- fig + plot_layout(guides = "collect") & theme(legend.position = "bottom")

  ggsave(outfile_png, fig, width = 12, height = 7.5, dpi = 300)
}

## Figure 6 : R20 (P > 20 mm/jour)
make_fig_indicator(
  prefix       = "R20_mon_ENSEMBLE",
  legend_title = "Days per season (P > 20 mm/day)",
  outfile_png  = "C:/Data/WBGT/Figure6_R20_gt20.png"
)

## Figure 7 : Heatwave days
make_fig_indicator(
  prefix       = "HW_mon_ENSEMBLE",
  legend_title = "Days per season in heatwaves",
  outfile_png  = "C:/Data/WBGT/Figure7_HW_days.png"
)

## =========================
## 6. HEAT INDEX > 40°C : FIGURE 5
## =========================

## fichiers HI (à adapter si besoin)
f_hist_hi <- file.path(sport_dir, "HI_mon_ENSEMBLE_HISTO_1991-2020.nc")
f_mid_hi  <- file.path(sport_dir, "HI_mon_ENSEMBLE_ssp585_2031-2060.nc")
f_late_hi <- file.path(sport_dir, "HI_mon_ENSEMBLE_ssp585_2071-2100.nc")

## nom de la variable HI>40 dans les NetCDF
var_hi <- "ndHI40"   # <-- change si différent

r_season_hi <- function(file, var, season = c("JJA","DJF")) {
  season <- match.arg(season)
  r <- rast(file, subds = var)
  n <- nlyr(r)
  nyear <- n / 12
  if (n %% 12 != 0) stop("n not multiple de 12 in ", file)

  months <- rep(1:12, length.out = n)
  idx <- if (season == "JJA") which(months %in% c(6,7,8)) else which(months %in% c(12,1,2))

  r_mean <- sum(r[[idx]]) / nyear
  r_mean[r_mean == 0] <- NA
  project(r_mean, crs_rob)
}

## extraire toutes les saisons pour échelle commune
layers_hi <- list(
  r_season_hi(f_hist_hi, var_hi, "JJA"),
  r_season_hi(f_mid_hi,  var_hi, "JJA"),
  r_season_hi(f_late_hi, var_hi, "JJA"),
  r_season_hi(f_hist_hi, var_hi, "DJF"),
  r_season_hi(f_mid_hi,  var_hi, "DJF"),
  r_season_hi(f_late_hi, var_hi, "DJF")
)
vals_all_hi <- unlist(lapply(layers_hi, values))
vals_all_hi <- vals_all_hi[is.finite(vals_all_hi)]
vmax_hi <- min(80, ceiling(max(vals_all_hi, na.rm = TRUE)))

make_panel_hi <- function(df, title, vmax) {
  make_panel_days(
    df           = df,
    title        = title,
    vmax         = vmax,
    legend_title = "Days per season (HI ≥ 40°C)",
    option       = "magma"
  )
}

## JJA
r_hist_jja_hi <- layers_hi[[1]]
r_mid_jja_hi  <- layers_hi[[2]]
r_late_jja_hi <- layers_hi[[3]]

df_hist_jja_hi <- to_df(r_hist_jja_hi, "days")
df_mid_jja_hi  <- to_df(r_mid_jja_hi,  "days")
df_late_jja_hi <- to_df(r_late_jja_hi, "days")

p_hist_jja_hi <- make_panel_hi(df_hist_jja_hi, "(a) JJA, 1991–2020",            vmax_hi)
p_mid_jja_hi  <- make_panel_hi(df_mid_jja_hi,  "(b) JJA, 2031–2060 (SSP5–8.5)", vmax_hi)
p_late_jja_hi <- make_panel_hi(df_late_jja_hi, "(c) JJA, 2071–2100 (SSP5–8.5)", vmax_hi)

## DJF
r_hist_djf_hi <- layers_hi[[4]]
r_mid_djf_hi  <- layers_hi[[5]]
r_late_djf_hi <- layers_hi[[6]]

df_hist_djf_hi <- to_df(r_hist_djf_hi, "days")
df_mid_djf_hi  <- to_df(r_mid_djf_hi,  "days")
df_late_djf_hi <- to_df(r_late_djf_hi, "days")

p_hist_djf_hi <- make_panel_hi(df_hist_djf_hi, "(d) DJF, 1991–2020",            vmax_hi)
p_mid_djf_hi  <- make_panel_hi(df_mid_djf_hi,  "(e) DJF, 2031–2060 (SSP5–8.5)", vmax_hi)
p_late_djf_hi <- make_panel_hi(df_late_djf_hi, "(f) DJF, 2071–2100 (SSP5–8.5)", vmax_hi)

fig5 <- (p_hist_jja_hi | p_mid_jja_hi | p_late_jja_hi) /
        (p_hist_djf_hi | p_mid_djf_hi | p_late_djf_hi)

fig5 <- fig5 + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("C:/Data/WBGT/Figure5_HI_gt40.png",
       fig5, width = 12, height = 7.5, dpi = 300)

## =========================
## 0. LIBRAIRIES & PATHS
## =========================
library(terra)
library(sf)
library(ggplot2)
library(viridis)
library(patchwork)
library(grid)      # pour unit()

## Dossiers à adapter
wbgt_dir  <- "C:/Data/WBGT/ENSEMBLE"
sport_dir <- "C:/Data/SPORT/ensemble"
shape_dir <- "C:/Data/SHAPE/NE110"

## =========================
## 1. SHAPEFILES & THÈME COMMUN
## =========================
coast   <- st_read(file.path(shape_dir, "ne_110m_coastline.shp"))
land    <- st_read(file.path(shape_dir, "ne_110m_land.shp"))

crs_rob <- "ESRI:54030"
coast_r <- st_transform(coast, crs_rob)
land_r  <- st_transform(land,  crs_rob)

theme_map <- theme_void() +
  theme(
    plot.title        = element_text(size = 10, face = "bold", hjust = 0.5),
    legend.position   = "bottom",
    legend.title      = element_text(size = 9),
    legend.text       = element_text(size = 8),
    legend.key.width  = unit(1.6, "cm"),
    legend.key.height = unit(0.4, "cm")
  )

## =========================
## 2. HELPERS GÉNÉRIQUES
## =========================

## 2.1 Raster -> saison (multi-subds, ex : WBGT)
r_season <- function(file, var, season = c("JJA","DJF")) {
  season <- match.arg(season)
  r <- rast(file, subds = var)
  n <- nlyr(r)
  nyear <- n / 12
  if (n %% 12 != 0) stop("n not multiple of 12 in ", file, " (var = ", var, ")")

  months <- rep(1:12, length.out = n)
  idx <- if (season == "JJA") which(months %in% c(6,7,8)) else which(months %in% c(12,1,2))

  r_mean <- sum(r[[idx]]) / nyear
  r_mean[r_mean == 0] <- NA
  project(r_mean, crs_rob)
}

## 2.2 Raster -> saison (1 variable par fichier)
r_season_1var <- function(file, season = c("JJA","DJF")) {
  season <- match.arg(season)
  r <- rast(file)
  n <- nlyr(r)
  nyear <- n / 12
  if (n %% 12 != 0) stop("n not multiple of 12 in ", file)

  months <- rep(1:12, length.out = n)
  idx <- if (season == "JJA") which(months %in% c(6,7,8)) else which(months %in% c(12,1,2))

  r_mean <- sum(r[[idx]]) / nyear
  r_mean[r_mean == 0] <- NA
  project(r_mean, crs_rob)
}

## 2.3 Raster -> saison (HI, subds = var_hi)
var_hi <- "ndHI40"   # adapte si différent dans tes NetCDF HI

r_season_hi <- function(file, var = var_hi, season = c("JJA","DJF")) {
  season <- match.arg(season)
  r <- rast(file, subds = var)
  n <- nlyr(r)
  nyear <- n / 12
  if (n %% 12 != 0) stop("n not multiple of 12 in ", file)

  months <- rep(1:12, length.out = n)
  idx <- if (season == "JJA") which(months %in% c(6,7,8)) else which(months %in% c(12,1,2))

  r_mean <- sum(r[[idx]]) / nyear
  r_mean[r_mean == 0] <- NA
  project(r_mean, crs_rob)
}

## 2.4 Raster -> data.frame
to_df <- function(r, nm = "value") {
  as.data.frame(r, xy = TRUE, na.rm = TRUE) |>
    setNames(c("x", "y", nm))
}

## 2.5 Panel "jours par saison" avec saturation à vmax
make_panel_days <- function(df, title, vmax, legend_title, option = "plasma") {
  # saturer > vmax
  df$days[df$days > vmax] <- vmax

  breaks_vec <- c(0, 5, 10, 20, 30, 40, 60, 80)
  breaks_vec <- breaks_vec[breaks_vec <= vmax]

  ggplot() +
    geom_sf(data = land_r, fill = "#242222", color = NA) +
    geom_raster(data = df, aes(x = x, y = y, fill = days)) +
    geom_sf(data = coast_r, color = "black", size = 0.1) +
    scale_fill_viridis(
      option = option,
      limits = c(0, vmax),
      breaks = breaks_vec,
      name   = legend_title
    ) +
    coord_sf(crs = crs_rob, expand = FALSE) +
    theme_map +
    labs(title = title)
}

## 2.6 Panel delta (PM–AM) avec clipping symétrique
make_panel_delta <- function(df, title, limits_days) {
  df$ddays[df$ddays >  limits_days] <-  limits_days
  df$ddays[df$ddays < -limits_days] <- -limits_days

  br <- c(-limits_days, -round(limits_days/2),
          0,
          round(limits_days/2), limits_days)

  ggplot() +
    geom_sf(data = land_r, fill = "#242222", color = NA) +
    geom_raster(data = df, aes(x = x, y = y, fill = ddays)) +
    geom_sf(data = coast_r, color = "black", size = 0.1) +
    scale_fill_gradient2(
      low      = "#2166ac",
      mid      = "white",
      high     = "#b2182b",
      midpoint = 0,
      limits   = c(-limits_days, limits_days),
      breaks   = br,
      name     = "Δ days (PM – AM)"
    ) +
    coord_sf(crs = crs_rob, expand = FALSE) +
    theme_map +
    labs(title = title)
}

## 2.7 Normalisation min-max globale (sur toutes couches)
normalize_indicator <- function(r) {
  if (!inherits(r, "SpatRaster")) {
    stop("normalize_indicator attend un terra::SpatRaster.")
  }
  rg <- terra::global(r, fun = "range", na.rm = TRUE)
  xmin <- min(rg[, 1], na.rm = TRUE)
  xmax <- max(rg[, 2], na.rm = TRUE)

  if (!is.finite(xmin) || !is.finite(xmax) || xmax <= xmin) {
    warning("Impossible de normaliser : xmax <= xmin ou valeurs non finies. Retourne NA.")
    return(r * NA)
  }
  (r - xmin) / (xmax - xmin)
}

## 2.8 Vérifier géométrie
check_geom_same <- function(...) {
  rast_list <- list(...)
  ref <- rast_list[[1]]
  for (i in seq_along(rast_list)[-1]) {
    if (!terra::compareGeom(ref, rast_list[[i]], stopOnError = FALSE)) {
      stop("Les rasters n'ont pas la même géométrie.")
    }
  }
}

## 2.9 SCRI (combinaison pondérée)
build_scri <- function(r_wbgt_aft_n, r_hi_n, r_rain_n, r_hw_n,
                       w_wbgt = 0.25, w_hi = 0.25,
                       w_rain = 0.25, w_hw = 0.25) {

  check_geom_same(r_wbgt_aft_n, r_hi_n, r_rain_n, r_hw_n)
  scri <- (w_wbgt * r_wbgt_aft_n +
           w_hi   * r_hi_n +
           w_rain * r_rain_n +
           w_hw   * r_hw_n)
  scri
}

## =========================
## 3. WBGT : FIGURES 2–3 (3×3 JJA / DJF)
## =========================

hist_f_585 <- file.path(wbgt_dir, "WBGT_mon_baseline_1991-2020_ensmean.nc")
mid_f_585  <- file.path(wbgt_dir, "WBGT_mon_ssp585_2031-2060_ensmean.nc")
late_f_585 <- file.path(wbgt_dir, "WBGT_mon_ssp585_2071-2100_ensmean.nc")

wbgt_vars <- c(
  Morning   = "ndays_wbgt_morn_ge32",
  Afternoon = "ndays_wbgt_aft_ge32",
  Evening   = "ndays_wbgt_eve_ge32"
)

make_panel_wbgt <- function(df, title, vmax) {
  make_panel_days(
    df           = df,
    title        = title,
    vmax         = vmax,
    legend_title = "Days per season (WBGT ≥ 32°C)",
    option       = "plasma"
  )
}

make_fig_wbgt <- function(season = "JJA") {
  files3      <- c(hist_f_585, mid_f_585, late_f_585)
  labels_rows <- c("Morning", "Afternoon", "Evening")
  labels_cols <- c("1991–2020", "2050 (2031–2060)", "2100 (2071–2100)")

  # échelle commune
  layers <- list()
  for (v in wbgt_vars) {
    for (f in files3) {
      layers[[length(layers) + 1]] <- r_season(f, v, season)
    }
  }
  all_vals <- unlist(lapply(layers, values))
  all_vals <- all_vals[is.finite(all_vals)]
  vmax <- min(80, ceiling(max(all_vals, na.rm = TRUE)))

  panels <- list()
  i <- 1
  row_id <- 1
  for (vname in names(wbgt_vars)) {
    var <- wbgt_vars[[vname]]
    col_id <- 1
    for (f in files3) {
      r  <- r_season(f, var, season)
      df <- to_df(r, "days")
      title <- paste0(labels_rows[row_id], " – ", labels_cols[col_id])
      panels[[i]] <- make_panel_wbgt(df, title, vmax)
      i      <- i + 1
      col_id <- col_id + 1
    }
    row_id <- row_id + 1
  }

  fig <- (panels[[1]] | panels[[2]] | panels[[3]]) /
         (panels[[4]] | panels[[5]] | panels[[6]]) /
         (panels[[7]] | panels[[8]] | panels[[9]])

  fig + plot_layout(guides = "collect") & theme(legend.position = "bottom")
}

fig2 <- make_fig_wbgt("JJA")
ggsave("C:/Data/WBGT/Figure2_WBGT_JJA_3x3.png",
       fig2, width = 12, height = 11, dpi = 300)

fig3 <- make_fig_wbgt("DJF")
ggsave("C:/Data/WBGT/Figure3_WBGT_DJF_3x3.png",
       fig3, width = 12, height = 11, dpi = 300)

## =========================
## 4. Δ WBGT (PM – AM) : FIGURE 4
## =========================

# JJA
r_hist_jja_aft  <- r_season(hist_f_585, wbgt_vars[["Afternoon"]], "JJA")
r_hist_jja_morn <- r_season(hist_f_585, wbgt_vars[["Morning"]],   "JJA")
r_mid_jja_aft   <- r_season(mid_f_585,  wbgt_vars[["Afternoon"]], "JJA")
r_mid_jja_morn  <- r_season(mid_f_585,  wbgt_vars[["Morning"]],   "JJA")
r_late_jja_aft  <- r_season(late_f_585, wbgt_vars[["Afternoon"]], "JJA")
r_late_jja_morn <- r_season(late_f_585, wbgt_vars[["Morning"]],   "JJA")

r_hist_jja_delta <- r_hist_jja_aft - r_hist_jja_morn
r_mid_jja_delta  <- r_mid_jja_aft  - r_mid_jja_morn
r_late_jja_delta <- r_late_jja_aft - r_late_jja_morn

# DJF
r_hist_djf_aft  <- r_season(hist_f_585, wbgt_vars[["Afternoon"]], "DJF")
r_hist_djf_morn <- r_season(hist_f_585, wbgt_vars[["Morning"]],   "DJF")
r_mid_djf_aft   <- r_season(mid_f_585,  wbgt_vars[["Afternoon"]], "DJF")
r_mid_djf_morn  <- r_season(mid_f_585,  wbgt_vars[["Morning"]],   "DJF")
r_late_djf_aft  <- r_season(late_f_585, wbgt_vars[["Afternoon"]], "DJF")
r_late_djf_morn <- r_season(late_f_585, wbgt_vars[["Morning"]],   "DJF")

r_hist_djf_delta <- r_hist_djf_aft - r_hist_djf_morn
r_mid_djf_delta  <- r_mid_djf_aft  - r_mid_djf_morn
r_late_djf_delta <- r_late_djf_aft - r_late_djf_morn

vals_all_delta <- c(
  values(r_hist_jja_delta), values(r_mid_jja_delta), values(r_late_jja_delta),
  values(r_hist_djf_delta), values(r_mid_djf_delta), values(r_late_djf_delta)
)
vals_all_delta <- vals_all_delta[is.finite(vals_all_delta)]
vmax_delta <- ceiling(max(abs(vals_all_delta), na.rm = TRUE))

df_hist_jja <- to_df(r_hist_jja_delta, "ddays")
df_mid_jja  <- to_df(r_mid_jja_delta,  "ddays")
df_late_jja <- to_df(r_late_jja_delta, "ddays")
df_hist_djf <- to_df(r_hist_djf_delta, "ddays")
df_mid_djf  <- to_df(r_mid_djf_delta,  "ddays")
df_late_djf <- to_df(r_late_djf_delta, "ddays")

p_hist_jja <- make_panel_delta(df_hist_jja, "(a) JJA, 1991–2020",            vmax_delta)
p_mid_jja  <- make_panel_delta(df_mid_jja,  "(b) JJA, 2031–2060 (SSP5–8.5)", vmax_delta)
p_late_jja <- make_panel_delta(df_late_jja, "(c) JJA, 2071–2100 (SSP5–8.5)", vmax_delta)

p_hist_djf <- make_panel_delta(df_hist_djf, "(d) DJF, 1991–2020",            vmax_delta)
p_mid_djf  <- make_panel_delta(df_mid_djf,  "(e) DJF, 2031–2060 (SSP5–8.5)", vmax_delta)
p_late_djf <- make_panel_delta(df_late_djf, "(f) DJF, 2071–2100 (SSP5–8.5)", vmax_delta)

fig4 <- (p_hist_jja | p_mid_jja | p_late_jja) /
        (p_hist_djf | p_mid_djf | p_late_djf)

fig4 <- fig4 + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("C:/Data/WBGT/Figure4_WBGT_delta_PM_AM.png",
       fig4, width = 12, height = 7.5, dpi = 300)

## =========================
## 5. R20 & HW : FIGURES 6–7
## =========================

make_fig_indicator <- function(prefix, legend_title, outfile_png) {
  f_hist <- file.path(sport_dir, paste0(prefix, "_HISTO_1991-2020.nc"))
  f_mid  <- file.path(sport_dir, paste0(prefix, "_ssp585_2031-2060.nc"))
  f_late <- file.path(sport_dir, paste0(prefix, "_ssp585_2071-2100.nc"))

  r_hist_jja <- r_season_1var(f_hist, "JJA")
  r_mid_jja  <- r_season_1var(f_mid,  "JJA")
  r_late_jja <- r_season_1var(f_late, "JJA")
  r_hist_djf <- r_season_1var(f_hist, "DJF")
  r_mid_djf  <- r_season_1var(f_mid,  "DJF")
  r_late_djf <- r_season_1var(f_late, "DJF")

  vals_all <- c(
    values(r_hist_jja), values(r_mid_jja), values(r_late_jja),
    values(r_hist_djf), values(r_mid_djf), values(r_late_djf)
  )
  vals_all <- vals_all[is.finite(vals_all)]
  vmax <- min(80, ceiling(max(vals_all, na.rm = TRUE)))

  df_hist_jja <- to_df(r_hist_jja, "days")
  df_mid_jja  <- to_df(r_mid_jja,  "days")
  df_late_jja <- to_df(r_late_jja, "days")
  df_hist_djf <- to_df(r_hist_djf, "days")
  df_mid_djf  <- to_df(r_mid_djf,  "days")
  df_late_djf <- to_df(r_late_djf, "days")

  p_hist_jja <- make_panel_days(df_hist_jja, "(a) JJA, 1991–2020",            vmax, legend_title)
  p_mid_jja  <- make_panel_days(df_mid_jja,  "(b) JJA, 2031–2060 (SSP5–8.5)", vmax, legend_title)
  p_late_jja <- make_panel_days(df_late_jja, "(c) JJA, 2071–2100 (SSP5–8.5)", vmax, legend_title)
  p_hist_djf <- make_panel_days(df_hist_djf, "(d) DJF, 1991–2020",            vmax, legend_title)
  p_mid_djf  <- make_panel_days(df_mid_djf,  "(e) DJF, 2031–2060 (SSP5–8.5)", vmax, legend_title)
  p_late_djf <- make_panel_days(df_late_djf, "(f) DJF, 2071–2100 (SSP5–8.5)", vmax, legend_title)

  fig <- (p_hist_jja | p_mid_jja | p_late_jja) /
         (p_hist_djf | p_mid_djf | p_late_djf)

  fig <- fig + plot_layout(guides = "collect") & theme(legend.position = "bottom")

  ggsave(outfile_png, fig, width = 12, height = 7.5, dpi = 300)
}

make_fig_indicator(
  prefix       = "R20_mon_ENSEMBLE",
  legend_title = "Days per season (P > 20 mm/day)",
  outfile_png  = "C:/Data/SPORT/Figure6_R20_gt20.png"
)

make_fig_indicator(
  prefix       = "HW_mon_ENSEMBLE",
  legend_title = "Days per season in heatwaves",
  outfile_png  = "C:/Data/SPORT/Figure7_HW_days.png"
)

## =========================
## 6. HI > 40 °C : FIGURE 5
## =========================

f_hist_hi <- file.path(sport_dir, "HI_mon_ENSEMBLE_HISTO_1991-2020.nc")
f_mid_hi  <- file.path(sport_dir, "HI_mon_ENSEMBLE_ssp585_2031-2060.nc")
f_late_hi <- file.path(sport_dir, "HI_mon_ENSEMBLE_ssp585_2071-2100.nc")

layers_hi <- list(
  r_season_hi(f_hist_hi, var_hi, "JJA"),
  r_season_hi(f_mid_hi,  var_hi, "JJA"),
  r_season_hi(f_late_hi, var_hi, "JJA"),
  r_season_hi(f_hist_hi, var_hi, "DJF"),
  r_season_hi(f_mid_hi,  var_hi, "DJF"),
  r_season_hi(f_late_hi, var_hi, "DJF")
)
vals_all_hi <- unlist(lapply(layers_hi, values))
vals_all_hi <- vals_all_hi[is.finite(vals_all_hi)]
vmax_hi <- min(80, ceiling(max(vals_all_hi, na.rm = TRUE)))

make_panel_hi <- function(df, title, vmax) {
  make_panel_days(
    df           = df,
    title        = title,
    vmax         = vmax,
    legend_title = "Days per season (HI ≥ 40°C)",
    option       = "magma"
  )
}

r_hist_jja_hi <- layers_hi[[1]]
r_mid_jja_hi  <- layers_hi[[2]]
r_late_jja_hi <- layers_hi[[3]]
r_hist_djf_hi <- layers_hi[[4]]
r_mid_djf_hi  <- layers_hi[[5]]
r_late_djf_hi <- layers_hi[[6]]

df_hist_jja_hi <- to_df(r_hist_jja_hi, "days")
df_mid_jja_hi  <- to_df(r_mid_jja_hi,  "days")
df_late_jja_hi <- to_df(r_late_jja_hi, "days")
df_hist_djf_hi <- to_df(r_hist_djf_hi, "days")
df_mid_djf_hi  <- to_df(r_mid_djf_hi,  "days")
df_late_djf_hi <- to_df(r_late_djf_hi, "days")

p_hist_jja_hi <- make_panel_hi(df_hist_jja_hi, "(a) JJA, 1991–2020",            vmax_hi)
p_mid_jja_hi  <- make_panel_hi(df_mid_jja_hi,  "(b) JJA, 2031–2060 (SSP5–8.5)", vmax_hi)
p_late_jja_hi <- make_panel_hi(df_late_jja_hi, "(c) JJA, 2071–2100 (SSP5–8.5)", vmax_hi)
p_hist_djf_hi <- make_panel_hi(df_hist_djf_hi, "(d) DJF, 1991–2020",            vmax_hi)
p_mid_djf_hi  <- make_panel_hi(df_mid_djf_hi,  "(e) DJF, 2031–2060 (SSP5–8.5)", vmax_hi)
p_late_djf_hi <- make_panel_hi(df_late_djf_hi, "(f) DJF, 2071–2100 (SSP5–8.5)", vmax_hi)

fig5 <- (p_hist_jja_hi | p_mid_jja_hi | p_late_jja_hi) /
        (p_hist_djf_hi | p_mid_djf_hi | p_late_djf_hi)

fig5 <- fig5 + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("C:/Data/WBGT/Figure5_HI_gt40.png",
       fig5, width = 12, height = 7.5, dpi = 300)

## =========================
## 7. SCRI multi-scenarios + historique : FIGURES 8–9
## =========================

scens_fut       <- c("ssp126", "ssp245", "ssp585")
scen_labels_fut <- c("SSP1–2.6", "SSP2–4.5", "SSP5–8.5")
seasons         <- c("JJA", "DJF")
period_years    <- list(
  mid  = "2031-2060",
  late = "2071-2100"
)

# fichiers HISTO
f_hist_wbgt <- file.path(wbgt_dir, "WBGT_mon_baseline_1991-2020_ensmean.nc")
f_hist_r20  <- file.path(sport_dir, "R20_mon_ENSEMBLE_HISTO_1991-2020.nc")
f_hist_hw   <- file.path(sport_dir, "HW_mon_ENSEMBLE_HISTO_1991-2020.nc")

## -- helpers lecture SCRI --

get_wbgt_aft_hist <- function(season) {
  r_season(f_hist_wbgt, wbgt_vars[["Afternoon"]], season)
}
get_wbgt_aft_fut <- function(scen, period, season) {
  years <- period_years[[period]]
  f <- file.path(wbgt_dir,
                 sprintf("WBGT_mon_%s_%s_ensmean.nc", scen, years))
  r_season(f, wbgt_vars[["Afternoon"]], season)
}

f_hist_hi_scri <- file.path(sport_dir, "HI_mon_ENSEMBLE_HISTO_1991-2020.nc")
get_hi_hist <- function(season) {
  r_season_hi(f_hist_hi_scri, var_hi, season)
}
get_hi_fut <- function(scen, period, season) {
  years <- period_years[[period]]
  f <- file.path(sport_dir,
                 sprintf("HI_mon_ENSEMBLE_%s_%s.nc", scen, years))
  r_season_hi(f, var_hi, season)
}

get_rain_hist <- function(season) {
  r_season_1var(f_hist_r20, season)
}
get_rain_fut <- function(scen, period, season) {
  years <- period_years[[period]]
  f <- file.path(sport_dir,
                 sprintf("R20_mon_ENSEMBLE_%s_%s.nc", scen, years))
  r_season_1var(f, season)
}

get_hw_hist <- function(season) {
  r_season_1var(f_hist_hw, season)
}
get_hw_fut <- function(scen, period, season) {
  years <- period_years[[period]]
  f <- file.path(sport_dir,
                 sprintf("HW_mon_ENSEMBLE_%s_%s.nc", scen, years))
  r_season_1var(f, season)
}

meta <- data.frame(
  scen   = character(),
  period = character(),
  season = character(),
  stringsAsFactors = FALSE
)
wbgt_list <- list()
hi_list   <- list()
rain_list <- list()
hw_list   <- list()

add_case <- function(scen, period, season, wbgt, hi, rain, hw) {
  idx <- length(wbgt_list) + 1L
  wbgt_list[[idx]] <<- wbgt
  hi_list[[idx]]   <<- hi
  rain_list[[idx]] <<- rain
  hw_list[[idx]]   <<- hw
  meta[idx, "scen"]   <<- scen
  meta[idx, "period"] <<- period
  meta[idx, "season"] <<- season
}

# baseline
for (s in seasons) {
  wbgt_s <- get_wbgt_aft_hist(s)
  hi_s   <- get_hi_hist(s)
  rain_s <- get_rain_hist(s)
  hw_s   <- get_hw_hist(s)
  add_case("hist", "hist", s, wbgt_s, hi_s, rain_s, hw_s)
}

# futurs
for (p in names(period_years)) {
  for (sc in scens_fut) {
    for (s in seasons) {
      wbgt_s <- get_wbgt_aft_fut(sc, p, s)
      hi_s   <- get_hi_fut(sc, p, s)
      rain_s <- get_rain_fut(sc, p, s)
      hw_s   <- get_hw_fut(sc, p, s)
      add_case(sc, p, s, wbgt_s, hi_s, rain_s, hw_s)
    }
  }
}

wbgt_stack <- rast(wbgt_list)
hi_stack   <- rast(hi_list)
rain_stack <- rast(rain_list)
hw_stack   <- rast(hw_list)

wbgt_stack_n <- normalize_indicator(wbgt_stack)
hi_stack_n   <- normalize_indicator(hi_stack)
rain_stack_n <- normalize_indicator(rain_stack)
hw_stack_n   <- normalize_indicator(hw_stack)

scri_list <- vector("list", nrow(meta))
for (i in seq_len(nrow(meta))) {
  scri_list[[i]] <- build_scri(
    wbgt_stack_n[[i]],
    hi_stack_n[[i]],
    rain_stack_n[[i]],
    hw_stack_n[[i]]
  )
}

scri_panel <- function(r, title) {
  df <- to_df(r, "scri")
  ggplot() +
    geom_sf(data = land_r, fill = "#242222", color = NA) +
    geom_raster(data = df, aes(x = x, y = y, fill = scri)) +
    geom_sf(data = coast_r, color = "black", size = 0.1) +
    scale_fill_viridis(
      option = "magma",
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.1),
      name   = "SCRI (0–1)"
    ) +
    coord_sf(crs = crs_rob, expand = FALSE) +
    theme_map +
    theme(
      legend.title = element_text(size = 9),
      legend.text  = element_text(size = 8)
    ) +
    labs(title = title)
}

get_idx <- function(scen, period, season) {
  which(meta$scen == scen &
        meta$period == period &
        meta$season == season)
}

## Figure 8 : 2031–2060 (mid) + baseline
idx_hist_jja  <- get_idx("hist",   "hist", "JJA")
idx_126_jja_m <- get_idx("ssp126", "mid",  "JJA")
idx_245_jja_m <- get_idx("ssp245", "mid",  "JJA")
idx_585_jja_m <- get_idx("ssp585", "mid",  "JJA")

idx_hist_djf  <- get_idx("hist",   "hist", "DJF")
idx_126_djf_m <- get_idx("ssp126", "mid",  "DJF")
idx_245_djf_m <- get_idx("ssp245", "mid",  "DJF")
idx_585_djf_m <- get_idx("ssp585", "mid",  "DJF")

p8a <- scri_panel(scri_list[[idx_hist_jja]],  "(a) JJA, 1991–2020")
p8b <- scri_panel(scri_list[[idx_126_jja_m]], "(b) JJA, SSP1–2.6, 2031–2060")
p8c <- scri_panel(scri_list[[idx_245_jja_m]], "(c) JJA, SSP2–4.5, 2031–2060")
p8d <- scri_panel(scri_list[[idx_585_jja_m]], "(d) JJA, SSP5–8.5, 2031–2060")

p8e <- scri_panel(scri_list[[idx_hist_djf]],  "(e) DJF, 1991–2020")
p8f <- scri_panel(scri_list[[idx_126_djf_m]], "(f) DJF, SSP1–2.6, 2031–2060")
p8g <- scri_panel(scri_list[[idx_245_djf_m]], "(g) DJF, SSP2–4.5, 2031–2060")
p8h <- scri_panel(scri_list[[idx_585_djf_m]], "(h) DJF, SSP5–8.5, 2031–2060")

fig8 <- (p8a | p8b | p8c | p8d) /
        (p8e | p8f | p8g | p8h)

fig8 <- fig8 + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("C:/Data/SPORT/Figure8_SCRI_2031-2060.png",
       fig8, width = 14, height = 7.5, dpi = 300)

## Figure 9 : 2071–2100 (late) + baseline
idx_126_jja_l <- get_idx("ssp126", "late", "JJA")
idx_245_jja_l <- get_idx("ssp245", "late", "JJA")
idx_585_jja_l <- get_idx("ssp585", "late", "JJA")

idx_126_djf_l <- get_idx("ssp126", "late", "DJF")
idx_245_djf_l <- get_idx("ssp245", "late", "DJF")
idx_585_djf_l <- get_idx("ssp585", "late", "DJF")

p9a <- scri_panel(scri_list[[idx_hist_jja]],  "(a) JJA, 1991–2020")
p9b <- scri_panel(scri_list[[idx_126_jja_l]], "(b) JJA, SSP1–2.6, 2071–2100")
p9c <- scri_panel(scri_list[[idx_245_jja_l]], "(c) JJA, SSP2–4.5, 2071–2100")
p9d <- scri_panel(scri_list[[idx_585_jja_l]], "(d) JJA, SSP5–8.5, 2071–2100")

p9e <- scri_panel(scri_list[[idx_hist_djf]],  "(e) DJF, 1991–2020")
p9f <- scri_panel(scri_list[[idx_126_djf_l]], "(f) DJF, SSP1–2.6, 2071–2100")
p9g <- scri_panel(scri_list[[idx_245_djf_l]], "(g) DJF, SSP2–4.5, 2071–2100")
p9h <- scri_panel(scri_list[[idx_585_djf_l]], "(h) DJF, SSP5–8.5, 2071–2100")

fig9 <- (p9a | p9b | p9c | p9d) /
        (p9e | p9f | p9g | p9h)

fig9 <- fig9 + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("C:/Data/SPORT/Figure9_SCRI_2071-2100.png",
       fig9, width = 14, height = 7.5, dpi = 300)
