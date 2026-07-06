library(terra)

# --- הגדרת ה-CRS ---
# אנחנו עובדים ב-UTM 36N שזה 32636
my_crs <- "EPSG:32636"

# הערה: בדקנו שהקבצים המקוריים כבר ב-32636, אז השורות הבאות בפועל לא
# משנות כלום - נשארות רק לתיעוד/ביטחון.

# טעינה של כל הקבצים עם הגדרה מפורשת של CRS
pre_eshtaol_b08 <- terra::rast("eshtaol2025/2025-04-21_B08_(Raw).tiff"); crs(pre_eshtaol_b08) <- my_crs
pre_eshtaol_b12 <- terra::rast("eshtaol2025/2025-04-21_B12_(Raw).tiff"); crs(pre_eshtaol_b12) <- my_crs
post_eshtaol_b08 <- terra::rast("eshtaol2025/2025-05-09_B08_(Raw).tiff"); crs(post_eshtaol_b08) <- my_crs
post_eshtaol_b12 <- terra::rast("eshtaol2025/2025-05-09_B12_(Raw).tiff"); crs(post_eshtaol_b12) <- my_crs
rec_eshtaol_b08 <- terra::rast("eshtaol2025/2026-05-14_B08_(Raw).tiff"); crs(rec_eshtaol_b08) <- my_crs
rec_eshtaol_b12 <- terra::rast("eshtaol2025/2026-05-14_B12_(Raw).tiff"); crs(rec_eshtaol_b12) <- my_crs

pre_jerus_b08 <- terra::rast("jerus2021/2021-08-13_B08_(Raw).tiff"); crs(pre_jerus_b08) <- my_crs
pre_jerus_b12 <- terra::rast("jerus2021/2021-08-13_B12_(Raw).tiff"); crs(pre_jerus_b12) <- my_crs
post_jerus_b08 <- terra::rast("jerus2021/2021-08-28_B08_(Raw).tiff"); crs(post_jerus_b08) <- my_crs
post_jerus_b12 <- terra::rast("jerus2021/2021-08-28_B12_(Raw).tiff"); crs(post_jerus_b12) <- my_crs
rec_jerus_b08 <- terra::rast("jerus2021/2022-08-23_B08_(Raw).tiff"); crs(rec_jerus_b08) <- my_crs
rec_jerus_b12 <- terra::rast("jerus2021/2022-08-23_B12_(Raw).tiff"); crs(rec_jerus_b12) <- my_crs

# --- רסמפל (Bilinear) ---
# אשתאול
pre_eshtaol_b12_res  <- terra::resample(pre_eshtaol_b12,  pre_eshtaol_b08,  method = "bilinear")
post_eshtaol_b12_res <- terra::resample(post_eshtaol_b12, post_eshtaol_b08, method = "bilinear")
rec_eshtaol_b12_res  <- terra::resample(rec_eshtaol_b12,  rec_eshtaol_b08,  method = "bilinear")

# ירושלים
pre_jerus_b12_res  <- terra::resample(pre_jerus_b12,  pre_jerus_b08,  method = "bilinear")
post_jerus_b12_res <- terra::resample(post_jerus_b12, post_jerus_b08, method = "bilinear")
rec_jerus_b12_res  <- terra::resample(rec_jerus_b12,  rec_jerus_b08,  method = "bilinear")

# --- חישוב NBR ---
# אשתאול
nbr_pre_eshtaol  <- (pre_eshtaol_b08  - pre_eshtaol_b12_res)  / (pre_eshtaol_b08  + pre_eshtaol_b12_res)
nbr_post_eshtaol <- (post_eshtaol_b08 - post_eshtaol_b12_res) / (post_eshtaol_b08 + post_eshtaol_b12_res)
nbr_rec_eshtaol  <- (rec_eshtaol_b08  - rec_eshtaol_b12_res)  / (rec_eshtaol_b08  + rec_eshtaol_b12_res)

# ירושלים
nbr_pre_jerus  <- (pre_jerus_b08  - pre_jerus_b12_res)  / (pre_jerus_b08  + pre_jerus_b12_res)
nbr_post_jerus <- (post_jerus_b08 - post_jerus_b12_res) / (post_jerus_b08 + post_jerus_b12_res)
nbr_rec_jerus  <- (rec_jerus_b08  - rec_jerus_b12_res)  / (rec_jerus_b08  + rec_jerus_b12_res)

# --- חישוב dNBR ---
# חומרה (Burn Severity)
dNBR_burn_eshtaol <- nbr_pre_eshtaol - nbr_post_eshtaol
dNBR_burn_jerus   <- nbr_pre_jerus   - nbr_post_jerus

# שיקום (Recovery - שנה אחרי)
dNBR_recovery_eshtaol <- nbr_pre_eshtaol - nbr_rec_eshtaol
dNBR_recovery_jerus   <- nbr_pre_jerus   - nbr_rec_jerus

# --- שלב 1: trim - הסרת שורות/עמודות חיצוניות שכולן NA ---
dNBR_burn_eshtaol     <- terra::trim(dNBR_burn_eshtaol)
dNBR_burn_jerus       <- terra::trim(dNBR_burn_jerus)
dNBR_recovery_eshtaol <- terra::trim(dNBR_recovery_eshtaol)
dNBR_recovery_jerus   <- terra::trim(dNBR_recovery_jerus)

# --- שלב 2: crop לתוך מלבן פנימי "נקי" ---
# אם ה-AOI המקורי לא היה מלבן, trim לבד לא מספיק כי יישארו NA
# באלכסון בתוך הקצוות. הפונקציה הזו מוצאת את המלבן הגדול ביותר
# שכולו valid data (ללא NA בכלל) ע"י כיווץ הדרגתי מהקצוות.
crop_to_valid_rect <- function(r, step = 0.01) {
  shrink <- 0
  repeat {
    e <- ext(r) 
    dx <- (e[2] - e[1]) * shrink
    dy <- (e[4] - e[3]) * shrink
    e2 <- ext(e[1] + dx, e[2] - dx, e[3] + dy, e[4] - dy)
    r_crop <- crop(r, e2)
    if (!any(is.na(values(r_crop)))) return(r_crop)
    shrink <- shrink + step
    if (shrink >= 0.5) {
      warning("לא נמצא מלבן נקי לגמרי - מחזיר את הגרסה המכווצת ביותר שנוסתה")
      return(r_crop)
    }
  }
}

dNBR_burn_eshtaol_clean     <- crop_to_valid_rect(dNBR_burn_eshtaol)
dNBR_burn_jerus_clean       <- crop_to_valid_rect(dNBR_burn_jerus)
dNBR_recovery_eshtaol_clean <- crop_to_valid_rect(dNBR_recovery_eshtaol)
dNBR_recovery_jerus_clean   <- crop_to_valid_rect(dNBR_recovery_jerus)

# נגדיר layout של 2 על 2
par(mfrow = c(2, 2))

# פלוט לכל 4 השכבות - עם asp=1 כדי שלא יתעוות
plot(dNBR_burn_eshtaol_clean, asp=1, main="Burn Severity - Eshtaol", col=rev(heat.colors(20)))
plot(dNBR_burn_jerus_clean, asp=1, main="Burn Severity - Jerusalem", col=rev(heat.colors(20)))
plot(dNBR_recovery_eshtaol_clean, asp=1, main="Recovery (1yr) - Eshtaol", col=terrain.colors(20))
plot(dNBR_recovery_jerus_clean, asp=1, main="Recovery (1yr) - Jerusalem", col=terrain.colors(20))

# החזרת ה-par למצב רגיל
par(mfrow = c(1, 1))

# רשימה של כל השכבות שחישבנו
layers_list <- list(
  Burn_Eshtaol = dNBR_burn_eshtaol_clean,
  Burn_Jerus = dNBR_burn_jerus_clean,
  Recovery_Eshtaol = dNBR_recovery_eshtaol_clean,
  Recovery_Jerus = dNBR_recovery_jerus_clean
)

# לופ פשוט לבדיקה
for (layer_name in names(layers_list)) {
  na_count <- global(layers_list[[layer_name]], fun=function(x) sum(is.na(x)))
  cat("שכבה:", layer_name, "-> מספר פיקסלים NA:", na_count[1,1], "\n")
}
