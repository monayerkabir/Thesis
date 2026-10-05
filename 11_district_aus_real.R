# ============================================================
# 11_district_aus_real.R
# Real BBS district Aus data, harvest years 2006-2025 (no simulation)
#
# Stage 1: heat -> Aus yield, district FE + year FE, SEs clustered
#          by district. Full period + split (2007-14 vs 2015-25)
# Stage 2: heat (June-July, C_late) -> yield -> stunted, per MICS
#          round and pooled. Indirect effect bootstrapped by district.
#
# Needs in D:/Monayer:
#   district_day_full.rds, ch_aus_v1.rds, aus_district_real_2006_2025.csv
# Output: results_district_real/
# ============================================================
library(dplyr); library(tidyr); library(purrr); library(readr)
library(lubridate); library(sandwich); library(lmtest)

setwd("D:/Monayer")
OUT <- "results_district_real"; dir.create(OUT, showWarnings = FALSE)
set.seed(2026)
B <- 200                                   # bootstrap draws

stages <- tibble::tribble(
  ~stage, ~start_md, ~end_md, ~t_hot,
  "s1", "03-15", "04-15", 35,
  "s2", "04-16", "05-31", 33,
  "s3", "06-01", "06-30", 35,
  "s4", "07-01", "07-31", 30) %>%
  mutate(start_i = as.integer(sub("-", "", start_md)),
         end_i   = as.integer(sub("-", "", end_md)))

# ---------- 1. District x year heat counts ----------
#   A = days MaxT > Yoshida limit;  C = days MaxT > district's own 90th pct
dd <- readRDS("district_day_full.rds") %>%
  mutate(year = year(date), mmdd = month(date) * 100L + mday(date))

heat <- map_dfr(seq_len(nrow(stages)), function(i) {
  s <- stages[i, ]
  x <- dd %>% filter(mmdd >= s$start_i, mmdd <= s$end_i)
  p <- x %>% group_by(district) %>%
    summarise(p90 = quantile(maxt_pred, 0.90, na.rm = TRUE), .groups = "drop")
  x %>% left_join(p, by = "district") %>% group_by(district, year) %>%
    summarise(A = sum(maxt_pred > s$t_hot, na.rm = TRUE),
              C = sum(maxt_pred > p90,     na.rm = TRUE), .groups = "drop") %>%
    mutate(stage = s$stage)
}) %>%
  pivot_wider(names_from = stage, values_from = c(A, C), names_sep = "_") %>%
  mutate(A_season = A_s1 + A_s2 + A_s3 + A_s4,
         C_season = C_s1 + C_s2 + C_s3 + C_s4,
         C_late   = C_s3 + C_s4)             # June + July

# ---------- 2. District production (all real) ----------
prod <- read_csv("aus_district_real_2006_2025.csv", show_col_types = FALSE) %>%
  mutate(log_yield = ifelse(area_acres > 0 & prod_tons > 0,
                            log(prod_tons / area_acres), NA_real_)) %>%
  group_by(district) %>%
  mutate(yield_anom = log_yield - mean(log_yield, na.rm = TRUE)) %>%
  ungroup()

cat("Name check (should be empty):",
    setdiff(unique(prod$district), unique(heat$district)), "\n")

panel <- prod %>% inner_join(heat, by = c("district", "harvest_year" = "year"))
real  <- panel %>% filter(!is.na(log_yield))
cat("Stage 1 sample:", nrow(real), "district-years,", n_distinct(real$district),
    "districts,", paste(range(real$harvest_year), collapse = "-"), "\n")

# ---------- 3. STAGE 1: heat -> yield ----------
fe_fit <- function(rhs, df = real, label = "full") {
  m  <- lm(as.formula(paste("log_yield ~", rhs,
                            "+ factor(district) + factor(harvest_year)")), data = df)
  ct <- coeftest(m, vcov = vcovCL(m, cluster = ~district))
  keep <- !grepl("factor|Intercept", rownames(ct))
  tibble(period = label, model = rhs, term = rownames(ct)[keep],
         est = ct[keep, 1], se = ct[keep, 2], p = ct[keep, 4],
         pct_per_day = 100 * (exp(ct[keep, 1]) - 1), n = nobs(m))
}
specs <- c("C_season", "A_season", "C_late",
           "C_s1 + C_s2 + C_s3 + C_s4", "A_s1 + A_s2 + A_s3 + A_s4")
early <- filter(real, harvest_year <= 2014)
late  <- filter(real, harvest_year >= 2015)

s1 <- bind_rows(
  map_dfr(specs, ~ fe_fit(.x, real,  "full")),
  map_dfr(specs, ~ fe_fit(.x, early, "2007-2014")),
  map_dfr(specs, ~ fe_fit(.x, late,  "2015-2025")))
write_csv(s1, file.path(OUT, "stage1_heat_to_yield.csv"))

show <- function(x) print(as.data.frame(mutate(x, across(where(is.numeric), ~ signif(.x, 3)))),
                          row.names = FALSE)
cat("\n=== STAGE 1: % change in Aus yield per extra hot day ===\n")
for (pp in c("full", "2007-2014", "2015-2025")) {
  cat("\n--", pp, "--\n"); show(filter(s1, period == pp) %>% select(-period))
}

# ---------- 4. Children: attach heat + yield for each window ----------
ch <- readRDS("ch_aus_v1.rds") %>% filter(age_24plus == 1, !is.na(stunted))

attach_win <- function(j) {
  sy <- paste0("aus_season_year_y", j)
  ch %>% select(child_id, district, all_of(sy)) %>%
    rename(season_year = all_of(sy)) %>%
    left_join(panel %>% select(district, harvest_year, C_late, yield_anom),
              by = c("district", "season_year" = "harvest_year")) %>%
    rename_with(~ paste0(.x, "_y", j), c(season_year, C_late, yield_anom))
}
chw <- reduce(lapply(1:3, attach_win), left_join, by = c("child_id", "district"),
              .init = ch %>% select(child_id, district, survey_year, stunted,
                                    cage_num, chweight))

cat("\nChildren by round:\n"); print(table(chw$survey_year))
cat("Share with yield attached, by window:\n")
print(chw %>% group_by(survey_year) %>%
        summarise(across(starts_with("yield_anom_y"), ~ round(mean(!is.na(.x)), 2))))

# ---------- 5. STAGE 2: mediation ----------
# a: yield_anom ~ heat + district FE + season-year FE        (weighted)
# b: stunted ~ yield_anom + heat + age + district FE + season-year FE (LPM, weighted)
mediate <- function(df, j) {
  h  <- paste0("C_late_y", j); y <- paste0("yield_anom_y", j)
  sy <- paste0("season_year_y", j)
  d  <- df %>% filter(!is.na(.data[[h]]), !is.na(.data[[y]]), !is.na(chweight))
  if (nrow(d) < 100 || n_distinct(d[[sy]]) < 2) return(NULL)
  fa <- as.formula(paste(y, "~", h, "+ factor(district) + factor(", sy, ")"))
  fb <- as.formula(paste("stunted ~", y, "+", h,
                         "+ cage_num + factor(district) + factor(", sy, ")"))
  est <- function(dd) {
    a  <- coef(lm(fa, data = dd, weights = chweight))[h]
    mb <- lm(fb, data = dd, weights = chweight)
    c(a = unname(a), b = unname(coef(mb)[y]), direct = unname(coef(mb)[h]))
  }
  pt    <- est(d)
  dists <- unique(d$district)
  by_d  <- split(d, d$district)
  bs <- replicate(B, {
    pick <- sample(dists, replace = TRUE)
    bd <- bind_rows(Map(function(x, k) { x$district <- paste0(x$district[1], "_", k); x },
                        by_d[pick], seq_along(pick)))
    e <- tryCatch(est(bd), error = function(e) c(a = NA, b = NA, direct = NA))
    c(e, indirect = unname(e["a"] * e["b"]))
  })
  ci <- function(v) quantile(v, c(.025, .975), na.rm = TRUE)
  tibble(window = paste0("y", j), n = nrow(d), n_districts = length(dists),
         a = pt["a"],      a_lo = ci(bs["a", ])[1],      a_hi = ci(bs["a", ])[2],
         b = pt["b"],      b_lo = ci(bs["b", ])[1],      b_hi = ci(bs["b", ])[2],
         direct = pt["direct"], d_lo = ci(bs["direct", ])[1], d_hi = ci(bs["direct", ])[2],
         indirect = pt["a"] * pt["b"],
         ind_lo = ci(bs["indirect", ])[1], ind_hi = ci(bs["indirect", ])[2])
}

s2_all <- list()
for (yr in c(2012, 2019, 2025, NA)) {
  lab  <- ifelse(is.na(yr), "pooled", as.character(yr))
  dsub <- if (is.na(yr)) chw else filter(chw, survey_year == yr)
  cat("\nRunning mediation:", lab, "...\n")
  s2 <- map_dfr(1:3, ~ mediate(dsub, .x))
  if (nrow(s2) == 0) { cat("  not enough data\n"); next }
  s2 <- mutate(s2, round = lab, .before = 1)
  s2_all[[lab]] <- s2
  cat("=== STAGE 2:", lab, "===\n"); show(s2)
}
write_csv(bind_rows(s2_all), file.path(OUT, "stage2_mediation.csv"))

cat("\nReading guide:\n",
    " a  = change in log yield per extra June-July hot day (expect < 0)\n",
    " b  = change in P(stunted) per 1.0 log-yield anomaly (expect < 0;\n",
    "      b * 0.1 = effect of a 10% higher yield)\n",
    " indirect = a * b, per extra hot day; CI from district bootstrap\n")
cat("\nDone. Results in", OUT, "\n")
