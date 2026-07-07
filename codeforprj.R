# ============================================================
# Wildfire dNBR analysis
# Burn + Recovery, with Recovery masked by Burn pixels
# ============================================================

library(terra)
library(ggplot2)

# ============================================================
# 0. הגדרות כלליות
# ============================================================

# CRS: UTM Zone 36N
my_crs <- "EPSG:32636"

# סף להגדרת פיקסלים שנפגעו מהשריפה
# אפשר לשנות ל-0.20 או 0.27 אם רוצים סינון מחמיר יותר
burn_threshold <- 0.10


# ============================================================
# 1. פונקציות עזר
# ============================================================

# הגדרת CRS לרסטר
set_my_crs <- function(r) {
  terra::crs(r) <- my_crs
  return(r)
}

# חישוב NBR
calc_nbr <- function(b08, b12) {
  nbr <- (b08 - b12) / (b08 + b12)
  return(nbr)
}

# סינון פיקסלים שנפגעו מהשריפה
# פיקסלים מעל הסף נשארים, כל היתר הופכים ל-NA
filter_burn_pixels <- function(r, threshold) {
  r_fire <- terra::ifel(r > threshold, r, NA)
  names(r_fire) <- "burn_dNBR"
  return(r_fire)
}


# ============================================================
# 2. טעינת הקבצים
# ============================================================

# ----------------------------
# Eshtaol 2025
# ----------------------------

pre_eshtaol_b08 <- terra::rast("eshtaol2025/2025-04-21_B08_(Raw).tiff")
pre_eshtaol_b12 <- terra::rast("eshtaol2025/2025-04-21_B12_(Raw).tiff")

post_eshtaol_b08 <- terra::rast("eshtaol2025/2025-05-09_B08_(Raw).tiff")
post_eshtaol_b12 <- terra::rast("eshtaol2025/2025-05-09_B12_(Raw).tiff")

rec_eshtaol_b08 <- terra::rast("eshtaol2025/2026-05-14_B08_(Raw).tiff")
rec_eshtaol_b12 <- terra::rast("eshtaol2025/2026-05-14_B12_(Raw).tiff")


# ----------------------------
# Jerusalem 2021
# ----------------------------

pre_jerus_b08 <- terra::rast("jerus2021/2021-08-13_B08_(Raw).tiff")
pre_jerus_b12 <- terra::rast("jerus2021/2021-08-13_B12_(Raw).tiff")

post_jerus_b08 <- terra::rast("jerus2021/2021-08-28_B08_(Raw).tiff")
post_jerus_b12 <- terra::rast("jerus2021/2021-08-28_B12_(Raw).tiff")

rec_jerus_b08 <- terra::rast("jerus2021/2022-08-23_B08_(Raw).tiff")
rec_jerus_b12 <- terra::rast("jerus2021/2022-08-23_B12_(Raw).tiff")


# ============================================================
# 3. הגדרת CRS לכל הקבצים
# ============================================================

pre_eshtaol_b08  <- set_my_crs(pre_eshtaol_b08)
pre_eshtaol_b12  <- set_my_crs(pre_eshtaol_b12)
post_eshtaol_b08 <- set_my_crs(post_eshtaol_b08)
post_eshtaol_b12 <- set_my_crs(post_eshtaol_b12)
rec_eshtaol_b08  <- set_my_crs(rec_eshtaol_b08)
rec_eshtaol_b12  <- set_my_crs(rec_eshtaol_b12)

pre_jerus_b08  <- set_my_crs(pre_jerus_b08)
pre_jerus_b12  <- set_my_crs(pre_jerus_b12)
post_jerus_b08 <- set_my_crs(post_jerus_b08)
post_jerus_b12 <- set_my_crs(post_jerus_b12)
rec_jerus_b08  <- set_my_crs(rec_jerus_b08)
rec_jerus_b12  <- set_my_crs(rec_jerus_b12)


# ============================================================
# 4. קביעת grid ייחוס לכל אתר
# ============================================================

# לכל אתר נשתמש ב-B08 של pre-fire כ-grid ייחוס
ref_eshtaol <- pre_eshtaol_b08
ref_jerus   <- pre_jerus_b08


# ============================================================
# 5. התאמת כל השכבות לאותו grid
# ============================================================

# ----------------------------
# Eshtaol
# ----------------------------

pre_eshtaol_b08_aligned <- pre_eshtaol_b08

pre_eshtaol_b12_aligned <- terra::resample(
  pre_eshtaol_b12,
  ref_eshtaol,
  method = "bilinear"
)

post_eshtaol_b08_aligned <- terra::resample(
  post_eshtaol_b08,
  ref_eshtaol,
  method = "bilinear"
)

post_eshtaol_b12_aligned <- terra::resample(
  post_eshtaol_b12,
  ref_eshtaol,
  method = "bilinear"
)

rec_eshtaol_b08_aligned <- terra::resample(
  rec_eshtaol_b08,
  ref_eshtaol,
  method = "bilinear"
)

rec_eshtaol_b12_aligned <- terra::resample(
  rec_eshtaol_b12,
  ref_eshtaol,
  method = "bilinear"
)


# ----------------------------
# Jerusalem
# ----------------------------

pre_jerus_b08_aligned <- pre_jerus_b08

pre_jerus_b12_aligned <- terra::resample(
  pre_jerus_b12,
  ref_jerus,
  method = "bilinear"
)

post_jerus_b08_aligned <- terra::resample(
  post_jerus_b08,
  ref_jerus,
  method = "bilinear"
)

post_jerus_b12_aligned <- terra::resample(
  post_jerus_b12,
  ref_jerus,
  method = "bilinear"
)

rec_jerus_b08_aligned <- terra::resample(
  rec_jerus_b08,
  ref_jerus,
  method = "bilinear"
)

rec_jerus_b12_aligned <- terra::resample(
  rec_jerus_b12,
  ref_jerus,
  method = "bilinear"
)


# ============================================================
# 6. חישוב NBR
# NBR = (NIR - SWIR) / (NIR + SWIR)
# B08 = NIR
# B12 = SWIR
# ============================================================

# ----------------------------
# Eshtaol
# ----------------------------

nbr_pre_eshtaol <- calc_nbr(
  pre_eshtaol_b08_aligned,
  pre_eshtaol_b12_aligned
)

nbr_post_eshtaol <- calc_nbr(
  post_eshtaol_b08_aligned,
  post_eshtaol_b12_aligned
)

nbr_rec_eshtaol <- calc_nbr(
  rec_eshtaol_b08_aligned,
  rec_eshtaol_b12_aligned
)


# ----------------------------
# Jerusalem
# ----------------------------

nbr_pre_jerus <- calc_nbr(
  pre_jerus_b08_aligned,
  pre_jerus_b12_aligned
)

nbr_post_jerus <- calc_nbr(
  post_jerus_b08_aligned,
  post_jerus_b12_aligned
)

nbr_rec_jerus <- calc_nbr(
  rec_jerus_b08_aligned,
  rec_jerus_b12_aligned
)


# ============================================================
# 7. חישוב dNBR
# ============================================================

# ----------------------------
# Burn severity
# ----------------------------

dNBR_burn_eshtaol <- nbr_pre_eshtaol - nbr_post_eshtaol
dNBR_burn_jerus   <- nbr_pre_jerus   - nbr_post_jerus

names(dNBR_burn_eshtaol) <- "burn_dNBR"
names(dNBR_burn_jerus)   <- "burn_dNBR"


# ----------------------------
# Recovery
# ----------------------------

dNBR_recovery_eshtaol <- nbr_pre_eshtaol - nbr_rec_eshtaol
dNBR_recovery_jerus   <- nbr_pre_jerus   - nbr_rec_jerus

names(dNBR_recovery_eshtaol) <- "recovery_dNBR"
names(dNBR_recovery_jerus)   <- "recovery_dNBR"


# ============================================================
# 8. סינון פיקסלים לפי BURN בלבד
# ============================================================

# כאן נבחרים הפיקסלים שהושפעו מהשריפה לפי dNBR burn.
# הפיקסלים האלו יהיו ה-mask גם עבור recovery.

dNBR_burn_eshtaol_fire_only <- filter_burn_pixels(
  dNBR_burn_eshtaol,
  burn_threshold
)

dNBR_burn_jerus_fire_only <- filter_burn_pixels(
  dNBR_burn_jerus,
  burn_threshold
)


# ============================================================
# 9. החלת אותו mask של BURN על RECOVERY
# ============================================================

# כלומר: Recovery נשאר רק באותם פיקסלים שעברו את סינון ה-BURN.

dNBR_recovery_eshtaol_fire_only <- terra::mask(
  dNBR_recovery_eshtaol,
  dNBR_burn_eshtaol_fire_only
)

names(dNBR_recovery_eshtaol_fire_only) <- "recovery_dNBR"


dNBR_recovery_jerus_fire_only <- terra::mask(
  dNBR_recovery_jerus,
  dNBR_burn_jerus_fire_only
)

names(dNBR_recovery_jerus_fire_only) <- "recovery_dNBR"


# ============================================================
# 10. המרה ל-data.frame עבור ggplot
# ============================================================

df_burn_eshtaol <- as.data.frame(
  dNBR_burn_eshtaol_fire_only,
  xy = TRUE,
  na.rm = TRUE
)

df_recovery_eshtaol <- as.data.frame(
  dNBR_recovery_eshtaol_fire_only,
  xy = TRUE,
  na.rm = TRUE
)

df_burn_jerus <- as.data.frame(
  dNBR_burn_jerus_fire_only,
  xy = TRUE,
  na.rm = TRUE
)

df_recovery_jerus <- as.data.frame(
  dNBR_recovery_jerus_fire_only,
  xy = TRUE,
  na.rm = TRUE
)


# ============================================================
# 11. מפות ggplot - כל מפה בנפרד
# ============================================================

# ----------------------------
# Eshtaol - Burn
# ----------------------------

gg_burn_eshtaol <- ggplot(
  df_burn_eshtaol,
  aes(x = x, y = y, fill = burn_dNBR)
) +
  geom_raster() +
  coord_equal() +
  scale_fill_gradientn(
    colours = rev(heat.colors(20)),
    name = "Burn dNBR"
  ) +
  labs(
    title = "Burn severity - Eshtaol",
    x = "X",
    y = "Y"
  ) +
  theme_minimal()


# ----------------------------
# Eshtaol - Recovery
# ----------------------------

gg_recovery_eshtaol <- ggplot(
  df_recovery_eshtaol,
  aes(x = x, y = y, fill = recovery_dNBR)
) +
  geom_raster() +
  coord_equal() +
  scale_fill_gradientn(
    colours = terrain.colors(20),
    name = "Recovery dNBR"
  ) +
  labs(
    title = "Recovery in burned pixels - Eshtaol",
    x = "X",
    y = "Y"
  ) +
  theme_minimal()


# ----------------------------
# Jerusalem - Burn
# ----------------------------

gg_burn_jerus <- ggplot(
  df_burn_jerus,
  aes(x = x, y = y, fill = burn_dNBR)
) +
  geom_raster() +
  coord_equal() +
  scale_fill_gradientn(
    colours = rev(heat.colors(20)),
    name = "Burn dNBR"
  ) +
  labs(
    title = "Burn severity - Jerusalem",
    x = "X",
    y = "Y"
  ) +
  theme_minimal()


# ----------------------------
# Jerusalem - Recovery
# ----------------------------

gg_recovery_jerus <- ggplot(
  df_recovery_jerus,
  aes(x = x, y = y, fill = recovery_dNBR)
) +
  geom_raster() +
  coord_equal() +
  scale_fill_gradientn(
    colours = terrain.colors(20),
    name = "Recovery dNBR"
  ) +
  labs(
    title = "Recovery in burned pixels - Jerusalem",
    x = "X",
    y = "Y"
  ) +
  theme_minimal()


# ============================================================
# 12. הצגת מפות ggplot
# ============================================================

gg_burn_eshtaol
gg_recovery_eshtaol

gg_burn_jerus
gg_recovery_jerus


# ============================================================
# 13. בדיקה נוספת עם terra::plot
# ============================================================

par(mfrow = c(2, 2))

terra::plot(
  dNBR_burn_eshtaol_fire_only,
  asp = 1,
  main = "Burn severity - Eshtaol",
  col = rev(heat.colors(20))
)

terra::plot(
  dNBR_recovery_eshtaol_fire_only,
  asp = 1,
  main = "Recovery in burned pixels - Eshtaol",
  col = terrain.colors(20)
)

terra::plot(
  dNBR_burn_jerus_fire_only,
  asp = 1,
  main = "Burn severity - Jerusalem",
  col = rev(heat.colors(20))
)

terra::plot(
  dNBR_recovery_jerus_fire_only,
  asp = 1,
  main = "Recovery in burned pixels - Jerusalem",
  col = terrain.colors(20)
)

par(mfrow = c(1, 1))


# ============================================================
# 14. בדיקת מספר פיקסלים לא-NA בכל שכבה
# ============================================================

cat("\nNon-NA pixels after filtering:\n")

cat(
  "Burn Eshtaol:",
  terra::global(!is.na(dNBR_burn_eshtaol_fire_only), "sum", na.rm = TRUE)[1, 1],
  "\n"
)

cat(
  "Recovery Eshtaol after burn mask:",
  terra::global(!is.na(dNBR_recovery_eshtaol_fire_only), "sum", na.rm = TRUE)[1, 1],
  "\n"
)

cat(
  "Burn Jerusalem:",
  terra::global(!is.na(dNBR_burn_jerus_fire_only), "sum", na.rm = TRUE)[1, 1],
  "\n"
)

cat(
  "Recovery Jerusalem after burn mask:",
  terra::global(!is.na(dNBR_recovery_jerus_fire_only), "sum", na.rm = TRUE)[1, 1],
  "\n"
)

plot(dNBR_burn_jerus)

class(dNBR_burn_jerus)

terra::plot(
  dNBR_burn_jerus,
  asp = 1,
  main = "Burn Severity - Jerusalem",
  col = rev(heat.colors(20))
)

dNBR_burn_jerus
install.packages("ggplot2")
library(terra)
library(ggplot2)



# לבחור איזו שריפה לבדוק

names(dNBR_burn_jerus) <- "dNBR"

df_dnbr <- as.data.frame(
  dNBR_burn_jerus,
  xy = TRUE,
  na.rm = TRUE
)

# גרף התפלגות ערכי dNBR
ggplot(df_dnbr, aes(x = dNBR)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 0.10, linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = 0.20, linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = 0.27, linetype = "dashed", linewidth = 1) +
  labs(
    title = "Distribution of dNBR values - Jerusalem",
    subtitle = "Dashed lines: possible fire thresholds",
    x = "dNBR",
    y = "Number of pixels"
  ) +
  theme_minimal()

burn_threshold <- 0.1

# כל מה שמתחת לסף הופך ל-NA
dNBR_burn_jerus_fire_only <- terra::ifel(
  dNBR_burn_jerus > burn_threshold,
  dNBR_burn_jerus,
  NA
)

