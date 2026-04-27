dir0 <- ifelse(Sys.info()[["sysname"]] == "Windows", "D:/", "/work/sph-huangj")
# source(file.path(dir0, "scripts", "f", "0conf_ML.R"))
setwd(file.path(dir0, "analysis", "ems120"))
pacman::p_load(readxl, writexl, data.table, tidyverse, scales, RColorBrewer, reticulate, lubridate, patchwork, zoo, broom, forcats, circlize) 
invisible(lapply(c("phe.f.R", "plot.f.R"), \(f) source(file.path(dir0, "scripts", "f", f))))
dir.dat <- "D:/data/ems120"
years <- 2013:2024
vars.basic <- c("电话", "地址", "地址类型", "开始受理时刻", "派车时间", "去程时间", "现场时间", "返程时间", "急救时间", "疾病类型", "接车地址经度", "接车地址纬度")
vars.basic.alias <- c("^病人电话号码|^联系电话.1|^联系电话", "^接车地址$|^接车地点$|^现场地址$", "地址类型",
	"^开始受理时刻|^开始时刻|^摘机时刻|^收到指令时刻", "^派车时间|^受理调度时间", "^去程时间|^去程在途时间",
	"^现场时间|^现场救援时间|^现场治疗时间|^现场急救时间", "^返程时间|^返程在途时间", "^急救时间|^急救反应时间", "疾病类型", "接车地址经度", "接车地址纬度")
vars.dxs <- c("性别", "年龄", "呼救原因", "病种判断", "病因", "伤病程度", "症状", "主诉", "病史", "初步诊断", "补充诊断")
vars.dxs.alias <- c("^性别$|病人性别|患者性别", "^年龄$|病人年龄|患者年龄", "^呼救原因|^呼叫原因", "病种判断", "^病因$|辅助诊断",
	"伤病程度|病情分级", "^症状$|患者症状", "^主诉$|病情\\(主诉\\)", "^病史$|现病史", "^初步诊断$", "初步诊断2|补充诊断")
vars <- c(vars.basic, vars.dxs); vars.alias <- c(vars.basic.alias, vars.dxs.alias)
vars.time <- c("开始受理时刻", "驶向现场时刻", "到达现场时刻", "病人上车时刻", "到达医院时刻")
dxs <- list(
	"Traffic" = "创伤-交通事故", "Intoxication" = "理化中毒",
	"Trauma" = c("创伤-暴力事件", "创伤-跌倒", "创伤-高处坠落", "创伤-其他原因", "其他-昏迷"),
	"CVD" = c("其他-胸闷", "神经系统疾病-脑卒中", "神经系统疾病-其他疾病", "心血管系统疾病-其他疾病", "心血管系统疾病-胸痛"),
	"Respiratory" = "呼吸系统疾病", "Mental" = "精神病", "NCD-Other" = c("泌尿系统疾病", "消化系统疾病", "内分泌系统疾病"),
	"Other" = c("妇产科", "儿科", "其他-其他症状"), "Death" = "其他-死亡"
)
dxs.raw <- unlist(dxs, use.names = FALSE)
dxs.grp <- setdiff(names(dxs), "Other") 
dxs.grp.color <- setNames(rainbow(length(dxs.grp), s = 0.8, v = 0.85), dxs.grp)
map_grp <- stack(dxs) %>% setNames(c("dx_raw", "dx_grp"))
grp_use <- c("low", "high")
roll3 <- function(x) zoo::rollmean(x, 3, fill = NA, align = "center")
sig_star <- function(p) dplyr::case_when(p < .0001 ~ "***", p < .001 ~ "**", p < .01 ~ "*", TRUE ~ "")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 读入数据
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (file.exists("dat0.list.rds")) dat0.list <- readRDS("dat0.list.rds")

dat0.list <- list()
for (year in years) {
	cat(sprintf("========== %d ==========\n", year))
	dat <- read_excel(file.path(dir.dat, "120数据", "清洗后数据", paste0(year, ".xlsx"))); names(dat) <- trimws(names(dat))
	col_names <- names(dat); new_names <- col_names; missing_vars <- c()
	for (i in seq_along(vars)) {
		all_m <- grep(vars.alias[i], col_names, value=TRUE)
		if (!length(all_m)) missing_vars <- c(missing_vars, vars[i]) else {
			best_m <- all_m[1]
			if (length(all_m) > 1) for (p in strsplit(vars.alias[i], "\\|")[[1]]) { m <- grep(p, all_m, value=TRUE); if (length(m)) { best_m <- m[1]; break } }
			new_names[match(best_m, col_names)] <- vars[i]
		}
	}
	if (length(missing_vars)) cat(sprintf("-> 警告：%d年缺 '%s'\n", year, paste(missing_vars, collapse="' 和 '")))
	names(dat) <- new_names
	dat <- dat %>% filter(!is.na(电话), !is.na(性别), !if_all(c(主诉, 病史, 初步诊断, 补充诊断), is.na)) %>% mutate( # 🏮
		across(any_of(vars.time), \(x) ymd_hms(trimws(x))),
		派车时间 = if ("派车时间" %in% names(.)) as.numeric(派车时间) else if (all(c("驶向现场时刻","开始受理时刻") %in% names(.))) as.numeric(驶向现场时刻 - 开始受理时刻, units="secs") else NA_real_,
		去程时间 = if ("去程时间" %in% names(.)) as.numeric(去程时间) else if (all(c("到达现场时刻","驶向现场时刻") %in% names(.))) as.numeric(到达现场时刻 - 驶向现场时刻, units="secs") else NA_real_,
		现场时间 = if ("现场时间" %in% names(.)) as.numeric(现场时间) else if (all(c("病人上车时刻","到达现场时刻") %in% names(.))) as.numeric(病人上车时刻 - 到达现场时刻, units="secs") else NA_real_,
		返程时间 = if ("返程时间" %in% names(.)) as.numeric(返程时间) else if (all(c("到达医院时刻","病人上车时刻") %in% names(.))) as.numeric(到达医院时刻 - 病人上车时刻, units="secs") else NA_real_,
		急救时间 = if ("急救时间" %in% names(.)) as.numeric(急救时间) else if (all(c("派车时间","去程时间","现场时间","返程时间") %in% names(.))) 派车时间 + 去程时间 + 现场时间 + 返程时间 else NA_real_
	) %>% select(any_of(vars)) %>% mutate(
		电话 = as.character(sub("^0+", "", 电话)), 电话 = ifelse(nchar(电话) == 11, 电话, NA_character_),
		phone = ifelse(is.na(电话), NA_character_, substring(电话, 4, 11)),
		年龄 = as.numeric(gsub("岁$", "", 年龄)), 时刻 = 开始受理时刻, 日期 = as.Date(时刻), hour = hour(时刻)
	) %>% select(-时刻)
	dat0.list[[as.character(year)]] <- dat
}
saveRDS(dat0.list, "dat0.list.rds")
write_xlsx("2019.train_dx.xlsx", dat0.list[["2019"]][1:10000, intersect(c(vars.dxs, "疾病类型"), names(dat0.list[["2019"]])), drop=FALSE])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Dx和Phone机器学习🩺
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (file.exists("dat.list.rds")) dat.list <- readRDS("dat.list.rds")
if (file.exists("dat1.list.rds")) dat1.list <- readRDS("dat1.list.rds")

reticulate::source_python("D:/scripts/main/ems120.py")
dat.list <- list()
for (year in years) {
	cat("Processing year:", year, "\n")
	key <- as.character(year); outfile <- paste0(key, ".xlsx")
	if (file.exists(outfile)) { dat.list[[key]] <- read_excel(outfile); next }
	res <- tryCatch({
		dat <- dat0.list[[key]]
		py_phone <- eval_phone_batch(dat$phone) # 🐂🐎
		dat <- dat %>% mutate(phone.sco = as.numeric(unlist(py_phone$sco)), phone.sco.reason = as.character(unlist(py_phone$reason)))
		dat_py <- dat %>% select(any_of(vars.dxs)) # 不能把datetime列传给python
		py_dx <- eval_dx_batch(dat_py, data_name = key) # 🐂🐎
		dat <- dat %>% mutate(
			疾病分类.关键词 = as.character(unlist(py_dx$kw)), 疾病分类.关键词.理由 = as.character(unlist(py_dx$kw_reason)),
			疾病分类.ML = as.character(unlist(py_dx$ml)), 疾病分类.ML.理由 = as.character(unlist(py_dx$ml_reason))
		)
		write_xlsx(dat, outfile); dat
	}, error = \(e) { cat("ERROR at year", key, "\n"); print(e); NULL })
	dat.list[[key]] <- res; rm(res); gc()
}
saveRDS(dat.list, "dat.list.rds")
lapply(dat.list, names)
for (year in 2014:2021) {
	tab <- table(dat.list[[as.character(year)]]$疾病类型, dat.list[[as.character(year)]]$疾病分类.ML, useNA = "no") # print(tab)
	concordance <- sum(diag(tab)) / sum(tab)
	cat("Year:", year, " | Concordance:", sprintf("%.1f%%", concordance * 100), "\n\n")
}

dat1.list <- list()
for (year in years) {
	year <- as.character(year); dat1 <- dat.list[[year]]
	if (is.null(dat1) || nrow(dat1)==0) { dat1.list[[year]] <- NULL; next }
	dat1.list[[year]] <- dat1 %>% add_count(电话, name="n_tel") %>% add_count(疾病分类.关键词, name="n_kw") %>% add_count(疾病分类.ML, name="n_ml") %>%
		mutate(电话 = ifelse(!is.na(电话) & n_tel > 5, NA_character_, 电话),
			疾病分类.关键词 = ifelse(!is.na(疾病分类.关键词) & n_kw < 50, NA_character_, trimws(疾病分类.关键词)),
			疾病分类.ML = ifelse(!is.na(疾病分类.ML) & n_ml < 50, NA_character_, trimws(疾病分类.ML)),
			phone.luck = case_when(is.na(phone.sco) ~ NA_character_, phone.sco <= 2 ~ "low", phone.sco <= 7 ~ "middle", phone.sco <= 10 ~ "high", TRUE ~ NA_character_),
			phone.luck = factor(phone.luck, levels=c("low", "middle", "high")),
			dx_raw = trimws(疾病分类.ML)
		) %>% left_join(map_grp, by = "dx_raw") %>% select(-n_tel, -n_kw, -n_ml)
}
saveRDS(dat1.list, "dat1.list.rds")
lapply(dat1.list, function(dat) table(dat$dx_grp, useNA = "ifany"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 表1 🦋
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
the_table1 <- purrr::imap_dfr(dat1.list, \(d, y) if (is.null(d) || !nrow(d)) tibble() else tibble(
  Year = as.integer(y),
  `EMS calls, N` = format(nrow(d), big.mark = ","),
  `Unique caller numbers, N` = format(n_distinct(d$电话, na.rm = TRUE), big.mark = ","),
  `Unique / total, %` = sprintf("%.2f%%", 100 * n_distinct(d$电话, na.rm = TRUE) / nrow(d)),
  `Age, mean (SD)` = sprintf("%.1f (%.1f)", mean(d$年龄, na.rm = TRUE), sd(d$年龄, na.rm = TRUE)),
  `Female, n (%)` = sprintf("%s (%.1f%%)", format(sum(d$性别 == "女", na.rm = TRUE), big.mark = ","), 100 * mean(d$性别 == "女", na.rm = TRUE)),
  `Low-luck, n (%)` = sprintf("%s (%.1f%%)", format(sum(d$phone.luck == "low", na.rm = TRUE), big.mark = ","), 100 * mean(d$phone.luck == "low", na.rm = TRUE)),
  `High-luck, n (%)` = sprintf("%s (%.1f%%)", format(sum(d$phone.luck == "high", na.rm = TRUE), big.mark = ","), 100 * mean(d$phone.luck == "high", na.rm = TRUE))
))
write_xlsx(the_table1, "Table1.xlsx"); data.table::fwrite(the_table1, "Tab1.txt", sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 图1. call-volume adjusted phenotype trends
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
yrs1 <- 2017:2024

dat.fig1.week <- bind_rows(lapply(yrs1, \(y){
	dat1.list[[as.character(y)]] %>%
		filter(!is.na(dx_grp), !is.na(日期), dx_grp %in% dxs.grp) %>%
		transmute(year = y, 日期 = as.Date(日期), dx = factor(as.character(dx_grp), levels = dxs.grp)) %>%
		mutate(week_start = floor_date(日期, "week", week_start = 1)) %>%
		group_by(year, week_start, dx) %>% summarise(call_count = n(), days = n_distinct(日期), .groups = "drop") %>%
		filter(days == 7)
}))

dat.fig1.week2 <- dat.fig1.week %>%
	group_by(year, week_start) %>%
	mutate(total_calls_week = sum(call_count), pct_week = call_count / total_calls_week) %>%
	ungroup()

dat.fig1.month <- bind_rows(lapply(yrs1, \(y){
	dat1.list[[as.character(y)]] %>%
		filter(!is.na(dx_grp), !is.na(日期), dx_grp %in% dxs.grp) %>%
		transmute(year = y, 日期 = as.Date(日期), dx = factor(as.character(dx_grp), levels = dxs.grp)) %>%
		mutate(month_start = floor_date(日期, "month")) %>%
		count(year, month_start, dx, name = "n") %>%
		group_by(year, month_start) %>%
		mutate(total_calls_month = sum(n), pct = n / total_calls_month) %>%
		ungroup()
}))

p1A <- ggplot(dat.fig1.week2, aes(week_start, pct_week, color = dx, group = dx)) +
	geom_line(linewidth = .9) +
	facet_wrap(~dx, scales = "free_y", ncol = 2) +
	scale_color_manual(values = dxs.grp.color[dxs.grp], guide = "none") +
	scale_y_continuous(labels = percent_format(accuracy = 1)) +
	scale_x_date(breaks = as.Date(paste0(yrs1, "-01-01")), labels = yrs1) +
	labs(title = "A. Weekly proportion of EMS calls by phenotype", x = NULL, y = "Percentage of weekly EMS calls") +
	fig_theme(base_size = 11) +
	theme(strip.text = element_text(face = "bold", size = 12), axis.text = element_text(face = "bold"),
		axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1),
		plot.title = element_text(face = "bold", size = 14))

p1B <- ggplot(dat.fig1.month, aes(month_start, fct_rev(dx), fill = pct)) +
	geom_tile() +
	scale_fill_viridis_c(name = NULL, labels = percent_format(accuracy = 1),
		guide = guide_colorbar(title = NULL, barwidth = unit(14, "cm"), barheight = unit(.5, "cm"))) +
	scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = c(0, 0)) +
	labs(title = "B. Monthly composition of the EMS spectrum", x = NULL, y = NULL) +
	fig_theme(base_size = 11) +
	theme(axis.text = element_text(face = "bold"), axis.text.y = element_text(size = 11),
		plot.title = element_text(face = "bold", size = 14),
		legend.position = "bottom", legend.direction = "horizontal")

Fig1 <- p1A / p1B + plot_layout(heights = c(2.2, 1))
Fig1
ggsave("Fig1.png", Fig1, width = 12, height = 10, dpi = 600)

save_xlsx("Fig1.data.xlsx", list(
	weekly_adjusted_data = dat.fig1.week2,
	weekly_total = dat.fig1.week2 %>% distinct(year, week_start, total_calls_week) %>% arrange(year, week_start),
	monthly_data = dat.fig1.month,
	monthly_total = dat.fig1.month %>% distinct(year, month_start, total_calls_month) %>% arrange(year, month_start),
	weekly_adjusted_summary = dat.fig1.week2 %>% group_by(dx) %>%
		summarise(weeks = n(), mean_weekly_pct = round(mean(pct_week), 4),
			sd_weekly_pct = round(sd(pct_week), 4),
			min_weekly_pct = round(min(pct_week), 4), max_weekly_pct = round(max(pct_week), 4), .groups = "drop"),
	yearly_adjusted_summary = dat.fig1.week2 %>% group_by(year, dx) %>%
		summarise(year_total_calls = sum(call_count), mean_weekly_pct = round(mean(pct_week), 4), .groups = "drop")
))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 图2. 🎭📱两组人的12年发病频率
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cap2 <- 1.2; xN <- 1.15
s2 <- bind_rows(lapply(years, \(y){
	d0 <- dat1.list[[as.character(y)]]
	if (is.null(d0) || !nrow(d0)) return(tibble())
	d0 <- d0 %>% filter(phone.luck %in% grp_use, !is.na(疾病分类.ML)) %>%
		transmute(group = as.character(phone.luck), dx_raw = trimws(疾病分类.ML)) %>%
		left_join(map_grp, by = "dx_raw") %>% filter(dx_grp %in% dxs.grp)
	tg <- table(d0$dx_grp, d0$group); ta <- table(d0$dx_grp); Nl <- sum(tg[, "low"]); Nh <- sum(tg[, "high"]); Tall <- sum(ta)
	tibble(year = y, disease = dxs.grp, n_low = as.numeric(tg[dxs.grp, "low"]), n_high = as.numeric(tg[dxs.grp, "high"]),
		N_low = Nl, N_high = Nh, pct_all = as.numeric(ta[dxs.grp]) / Tall,
		pct_low = as.numeric(tg[dxs.grp, "low"]) / Nl, pct_high = as.numeric(tg[dxs.grp, "high"]) / Nh) %>%
		mutate(enrich_low = pct_low / pct_all, enrich_high = pct_high / pct_all, RR = pct_high / pct_low,
			p = purrr::pmap_dbl(list(n_high, N_high, n_low, N_low), \(a, A, b, B) suppressWarnings(chisq.test(matrix(c(a, A-a, b, B-b), 2))$p.value)),
			sig = !is.na(p) & p < .01, lo = pmax(enrich_low, 1 / cap2), hi = pmin(enrich_high, cap2), lf = enrich_low < 1 / cap2, hf = enrich_high > cap2)
}))
# 🏮
Fig2 <- wrap_plots(lapply(dxs.grp, \(dx){
	col <- dxs.grp.color[dx]; d <- s2 %>% filter(disease == dx)
	ggplot(d, aes(y = year)) + geom_vline(xintercept = 1) +
		geom_vline(xintercept = xN, linetype = "dashed", color = "grey40", linewidth = .6) +
		geom_segment(aes(x = 1, xend = lo, yend = year), linetype = "dashed", color = "grey70", linewidth = .9) +
		geom_segment(aes(x = 1, xend = hi, yend = year), linetype = "dashed", color = col, linewidth = .9) +
		geom_point(aes(x = lo), color = "grey50", size = 3) + geom_point(aes(x = hi), color = col, size = 3) +
		geom_text(data = d %>% filter(sig), aes(x = hi, label = "*"), hjust = -.2, vjust = .3, size = 5, fontface = "bold") +
		geom_text(data = d %>% filter(lf), aes(x = lo, label = "<"), hjust = 1.2) +
		geom_text(data = d %>% filter(hf), aes(x = hi, label = ">"), hjust = 0) +
		geom_text(aes(x = xN - .005, label = n_low), hjust = 1, size = 3.2, fontface = "bold", color = "grey60") +
		geom_text(aes(x = xN + .005, label = n_high), hjust = 0, size = 3.2, fontface = "bold", color = col) +
		scale_x_continuous(limits = c(.9, 1.2)) + scale_y_continuous(breaks = years, labels = years) +
		labs(title = dx, x = NULL, y = NULL) + theme_minimal() +
		theme(axis.text = element_text(face = "bold"), axis.title = element_text(face = "bold"), axis.line = element_line(), legend.position = "none")
}), nrow = 4, ncol = 2)
Fig2; ggsave("Fig2.png", Fig2, width = 12, height = 10, dpi = 600)

save_xlsx("Fig2.data.xlsx", list(
	panel_data = s2 %>% select(year, disease, n_low, n_high, pct_all, pct_low, pct_high, enrich_low, enrich_high, RR, p, sig, lo, hi),
	summary = s2 %>% group_by(disease) %>% summarise(mean_RR = round(mean(RR, na.rm = TRUE), 3), min_RR = round(min(RR, na.rm = TRUE), 3), max_RR = round(max(RR, na.rm = TRUE), 3), n_sig_years = sum(sig, na.rm = TRUE), .groups = "drop")
))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 图3. 急救🚑时间 (Onsite)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dxs1 <- dxs.grp
dat3 <- bind_rows(lapply(years, \(y){
	d0 <- dat1.list[[as.character(y)]]
	if (is.null(d0) || !nrow(d0)) return(tibble())
	d0 %>% filter(dx_grp %in% dxs1) %>%
		transmute(Year = y, dx_grp = factor(as.character(dx_grp), levels = dxs1), phone.luck = as.character(phone.luck), onsite = suppressWarnings(as.numeric(现场时间) / 60))
}))

dat3a <- dat3 %>% group_by(dx_grp, Year) %>% summarise(mean = mean(onsite, na.rm = TRUE), .groups = "drop")
dat3b <- dat3 %>% filter(phone.luck %in% c("low", "high")) %>%
	group_by(Year, split = phone.luck) %>% summarise(mean = mean(onsite, na.rm = TRUE), .groups = "drop") %>%
	pivot_wider(names_from = split, values_from = mean) %>%
	mutate(diff = high - low,
		p = pmap_dbl(list(Year), \(yy){
			d0 <- dat3 %>% filter(Year == yy, phone.luck %in% c("low", "high"))
			suppressWarnings(tryCatch(wilcox.test(d0$onsite[d0$phone.luck == "high"], d0$onsite[d0$phone.luck == "low"])$p.value, error = \(e) NA_real_))
		}),
		p_adj = p.adjust(p, "BH"), sig = !is.na(p_adj) & p_adj < .01)

gm <- mean(dat3b$high, na.rm = TRUE); rng <- range(c(dat3b$low, dat3b$high), na.rm = TRUE); eps <- .05 * diff(rng)
dat3b <- dat3b %>% mutate(star_x = ifelse(high >= low, high + eps, low - eps))

# A图：Y轴改成 Year
Fig3a <- ggplot(dat3a, aes(x = mean, y = Year, color = dx_grp)) +
	geom_point(size = 3) +
	scale_color_manual(values = dxs.grp.color[dxs1], name = NULL, drop = FALSE) +
	scale_y_continuous(breaks = years, labels = years) +
	labs(title = "A. Mean onsite time by phenotype and year", x = "Time (mins)", y = "Year") +
	theme_minimal(base_size = 12) +
	theme(axis.title = element_text(face = "bold"), axis.text = element_text(face = "bold"), axis.line = element_line(), plot.title = element_text(face = "bold"))

Fig3b <- ggplot(dat3b, aes(y = Year)) +
	geom_segment(aes(x = low, xend = high, yend = Year), color = "black", linewidth = .5) +
	geom_point(aes(x = low), color = "grey50", size = 3) +
	geom_point(aes(x = high), color = "blue", size = 3, shape = 17) +
	geom_vline(xintercept = gm, color = "blue", linetype = "dashed", linewidth = .9) +
	geom_text(data = dplyr::filter(dat3b, sig), aes(x = star_x, label = "*"), fontface = "bold", size = 5, vjust = .35) +
	labs(title = "B. High vs low onsite time", x = "Time (mins)", y = "Year") +
	scale_y_continuous(breaks = years, labels = years) +
	theme_minimal(base_size = 12) +
	theme(axis.title = element_text(face = "bold"), axis.text = element_text(face = "bold"), axis.line = element_line(), plot.title = element_text(face = "bold"))

Fig3 <- Fig3a | Fig3b
Fig3; ggsave("Fig3.png", Fig3, width = 12, height = 10, dpi = 600)

save_xlsx("Fig3.data.xlsx", list(
	panelA_data = dat3a,
	panelB_data = dat3b,
	by_dx_high_low = dat3 %>% filter(phone.luck %in% c("low", "high")) %>% group_by(dx_grp, phone.luck) %>% summarise(mean = round(mean(onsite, na.rm = TRUE), 2), sd = round(sd(onsite, na.rm = TRUE), 2), .groups = "drop")
))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 图4. 疫情影响 🚫
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
date_begin <- as.Date("2022-03-14"); date_end <- as.Date("2022-03-20"); flank_days <- 7
dat1 <- dat1.list[[as.character(year(date_begin))]] %>%
	transmute(日期 = as.Date(日期), dx_grp = factor(as.character(dx_grp), levels = dxs.grp), phone.luck = as.character(phone.luck)) %>%
	filter(dx_grp %in% dxs.grp, phone.luck %in% c("low", "high"), between(日期, date_begin - flank_days, date_end + flank_days)) %>%
	count(dx_grp, phone.luck, 日期, name = "count") %>%
	complete(dx_grp, phone.luck, 日期, fill = list(count = 0)) %>%
	group_by(dx_grp, phone.luck) %>% arrange(日期) %>% mutate(count3 = roll3(count)) %>% ungroup() %>%
	mutate(post = as.integer(between(日期, date_begin, date_end)), high = as.integer(phone.luck == "high"))
dat4A <- dat1

p4A <- ggplot(dat1, aes(日期, count3, color = phone.luck, group = interaction(dx_grp, phone.luck))) +
	geom_vline(xintercept = c(date_begin, date_end), linetype = "dashed", color = "orange", linewidth = 1) +
	geom_line(linewidth = 1, na.rm = TRUE) +
	scale_color_manual(values = c(low = "grey60", high = "#D55E00"), guide = "none") +
	facet_wrap(~dx_grp, scales = "free_y", ncol = 1) +
	scale_x_date(labels = date_format("%b %d", locale = "en")) +
	scale_y_continuous(breaks = pretty_breaks(n = 2)) +
	labs(title = "A. EMS calls around 03/2022 PHSM", x = NULL, y = "3-day rolling mean") +
	fig_theme(base_size = 12) +
	theme( 
		axis.title = element_text(face = "bold", size = 15),
		axis.text = element_text(face = "bold", size = 13),
		strip.text = element_text(face = "bold", size = 14),
		plot.title = element_text(face = "bold", size = 16)
	)

dat1 <- dat1 %>% filter(日期 <= date_end) %>% group_by(dx_grp) %>% group_modify(\(.x, .y){
	fit <- tryCatch(glm(count ~ high + post + high:post, family = poisson(), data = .x), error = \(e) NULL)
	if (is.null(fit)) return(tibble(did_rr = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_))
	x <- broom::tidy(fit) %>% filter(term == "high:post")
	tibble(did_rr = exp(x$estimate), lo = exp(x$estimate - 1.96 * x$std.error), hi = exp(x$estimate + 1.96 * x$std.error), p = x$p.value)
}) %>% ungroup() %>% mutate(sig = sig_star(p), dx_grp = factor(dx_grp, levels = dxs.grp)) 
dat4B <- dat1

p4B <- ggplot(dat1, aes(x = did_rr, y = fct_rev(dx_grp), color = dx_grp)) +
	geom_vline(xintercept = 1, linetype = "dashed") +
	geom_errorbar(aes(xmin = lo, xmax = hi), width = .2, orientation = "y", linewidth = .8) +
	geom_point(size = 3) +
	geom_text(aes(label = sig), hjust = -.4, fontface = "bold") +
	scale_color_manual(values = dxs.grp.color[dxs.grp], guide = "none") +
	labs(title = "B. DID effect during PHSM", x = "DID rate ratio", y = NULL) +
	fig_theme(base_size = 12) +
	theme(
		axis.title = element_text(face = "bold", size = 15),
		axis.text = element_text(face = "bold", size = 13),
		plot.title = element_text(face = "bold", size = 16)
	)

date_half <- as.Date("2022-11-11"); date_full <- as.Date("2022-12-07")

dat1 <- bind_rows(
	dat1.list[[as.character(year(date_half))]] %>%
		transmute(日期 = as.Date(日期), dx_grp = factor(as.character(dx_grp), levels = dxs.grp), phone.luck = as.character(phone.luck)) %>%
		filter(dx_grp %in% dxs.grp, phone.luck %in% c("low", "high"), between(日期, date_half - 10, as.Date("2022-12-31"))),
	dat1.list[[as.character(year(date_half) + 1)]] %>%
		transmute(日期 = as.Date(日期), dx_grp = factor(as.character(dx_grp), levels = dxs.grp), phone.luck = as.character(phone.luck)) %>%
		filter(dx_grp %in% dxs.grp, phone.luck %in% c("low", "high"), between(日期, as.Date("2023-01-01"), date_full + 35))
) %>%
	count(dx_grp, phone.luck, 日期, name = "count") %>%
	complete(dx_grp, phone.luck, 日期, fill = list(count = 0)) %>%
	group_by(dx_grp, phone.luck) %>%
	arrange(日期) %>%
	mutate(count3 = roll3(count)) %>%
	ungroup() %>%
	mutate(period = factor(case_when(
		日期 < date_half ~ "Pre-open",
		日期 < date_full ~ "Initial relaxation",
		TRUE ~ "Full opening"
	), levels = c("Pre-open", "Initial relaxation", "Full opening")))

dat4C <- dat1

p4C <- ggplot(dat4C, aes(日期, count3, color = phone.luck, group = interaction(dx_grp, phone.luck))) +
	geom_vline(xintercept = c(date_half, date_full), linetype = "dashed", color = "orange", linewidth = 1) +
	geom_line(linewidth = 1, na.rm = TRUE) +
	scale_color_manual(values = c(low = "grey60", high = "#D55E00"), guide = "none") +
	facet_wrap(~dx_grp, scales = "free_y", ncol = 1) +
	scale_x_date(labels = date_format("%b %d", locale = "en")) +
	scale_y_continuous(breaks = pretty_breaks(n = 2)) +
	labs(title = "C. Reopening period and EMS rebound pattern", x = NULL, y = "3-day rolling mean") +
	fig_theme(base_size = 12) +
	theme(
		axis.title = element_text(face = "bold", size = 15),
		axis.text = element_text(face = "bold", size = 13),
		strip.text = element_text(face = "bold", size = 14),
		plot.title = element_text(face = "bold", size = 16)
	)

dat4D <- dat4C %>%
	group_by(dx_grp, phone.luck, period) %>%
	summarise(mean_n = mean(count, na.rm = TRUE), .groups = "drop") %>%
	group_by(dx_grp, phone.luck) %>%
	mutate(rr_vs_pre = mean_n / mean_n[period == "Pre-open"]) %>%
	ungroup() %>%
	filter(period == "Full opening") %>%
	select(dx_grp, phone.luck, period, rr_vs_pre) %>%
	pivot_wider(names_from = phone.luck, values_from = rr_vs_pre) %>%
	mutate(dx_grp = factor(as.character(dx_grp), levels = dxs.grp))

p4D <- ggplot(dat4D, aes(y = fct_rev(dx_grp))) +
	geom_vline(xintercept = 1, linetype = "dashed") +
	geom_segment(aes(x = low, xend = high, yend = fct_rev(dx_grp)), color = "black", linewidth = .5) +
	geom_point(aes(x = low), color = "grey60", shape = 17, size = 3.2) +
	geom_point(aes(x = high), color = "#D55E00", shape = 17, size = 3.2) +
	scale_y_discrete(limits = rev(dxs.grp)) +
	labs(title = "D. Relative rebound compared with pre-open baseline", x = "Rate ratio vs pre-open", y = NULL) +
	fig_theme(base_size = 12) +
	theme(
		axis.title = element_text(face = "bold", size = 15),
		axis.text = element_text(face = "bold", size = 13),
		plot.title = element_text(face = "bold", size = 16),
		legend.position = "none"
	)

Fig4 <- (p4A | p4C) / plot_spacer() / (p4B | p4D) + plot_layout(heights = c(2.2, 0.12, 1))
Fig4
ggsave("Fig4.png", Fig4, width = 12, height = 12, dpi = 600)

save_xlsx("Fig4.data.xlsx", list(
	Fig4_panelA = dat4A %>% select(dx_grp, phone.luck, 日期, count, count3, post, high),
	Fig4_panelB = dat4B,
	Fig4_panelC = dat4C %>% select(dx_grp, phone.luck, 日期, count, count3, period),
	Fig4_panelD = dat4D
))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 图5 房价 sensitivity analyses
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pacman::p_load(sf)
year_use <- "2021"; dist_max_m <- 1000
house_sf <- read_excel(file.path(dir.dat, "深圳房价.xlsx")) %>% transmute(house.id = row_number(), house.price = as.numeric(房价), Lon = as.numeric(Lon), Lat = as.numeric(Lat)) %>% filter(is.finite(Lon), is.finite(Lat), is.finite(house.price)) %>% st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE) %>% st_transform(3857)
X_sf <- dat1.list[[year_use]] %>% transmute(X.id = row_number(), 地址类型, phone.sco = as.numeric(phone.sco), phone.luck = as.character(phone.luck), lon = as.numeric(接车地址经度), lat = as.numeric(接车地址纬度)) %>% filter(地址类型 == "住宅区", phone.luck %in% c("low", "high"), is.finite(lon), is.finite(lat)) %>% st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>% st_transform(3857)
idx <- st_nearest_feature(X_sf, house_sf)
dat.figS5 <- X_sf %>% mutate(house.price = house_sf$house.price[idx], dist_m = as.numeric(st_distance(X_sf, house_sf[idx, ], by_element = TRUE))) %>% st_drop_geometry() %>% mutate(house.price = if_else(dist_m <= dist_max_m, house.price, NA_real_))
dat_bin <- dat.figS5 %>% filter(is.finite(house.price), is.finite(phone.sco)) %>% mutate(logp = log10(house.price)); h <- hist(dat_bin$logp, breaks = 10, plot = FALSE)
s6 <- dat_bin %>% mutate(bin = cut(logp, breaks = h$breaks, include.lowest = TRUE)) %>% group_by(bin) %>% summarise(x = round(mean(range(logp)), 2), n = n(), y = round(mean(phone.sco), 2), sd = round(sd(phone.sco), 2), .groups = "drop") %>% arrange(x)
#🏮
png("Fig5.png", width = 10, height = 5, units = "in", res = 300); par(mar = c(5,4,3,5)+.2, font.lab = 2, font.axis = 2)
hh <- hist(dat_bin$logp, breaks = h$breaks, freq = TRUE, col = "grey85", border = "grey40", main = "", xlab = "log10(house price)", ylab = "")
par(new = TRUE); plot(s6$x, s6$y, ylim = range(c(s6$y-s6$sd, s6$y+s6$sd), na.rm = TRUE), xlim = range(hh$breaks), axes = FALSE, xlab = NA, ylab = NA, pch = 16, cex = 1.2, col = "blue")
arrows(s6$x, s6$y-s6$sd, s6$x, s6$y+s6$sd, angle = 90, code = 3, length = .05, col = "grey60")
axis(side = 4, font.axis = 2); mtext(side = 4, line = 3, "Phone luck score", col = "blue", font = 2)
dev.off()
save_xlsx("Fig5.data.xlsx", list(raw = dat.fig5, binned = s6))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 图S1. raw phenotype trends
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pS1A <- ggplot(dat.fig1.week, aes(week_start, call_count, color = dx, group = dx)) +
	geom_line(linewidth = .9) +
	facet_wrap(~dx, scales = "free_y", ncol = 2) +
	scale_color_manual(values = dxs.grp.color[dxs.grp], guide = "none") +
	scale_x_date(breaks = as.Date(paste0(yrs1, "-01-01")), labels = yrs1) +
	labs(title = "A. Weekly EMS demand by phenotype", x = NULL, y = "Calls per complete week") +
	fig_theme(base_size = 11) +
	theme(strip.text = element_text(face = "bold", size = 12), axis.text = element_text(face = "bold"),
		axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1),
		plot.title = element_text(face = "bold", size = 14))

pS1B <- ggplot(dat.fig1.month, aes(month_start, fct_rev(dx), fill = pct)) +
	geom_tile() +
	scale_fill_viridis_c(name = NULL, labels = percent_format(accuracy = 1),
		guide = guide_colorbar(title = NULL, barwidth = unit(14, "cm"), barheight = unit(.5, "cm"))) +
	scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = c(0, 0)) +
	labs(title = "B. Monthly composition of the EMS spectrum", x = NULL, y = NULL) +
	fig_theme(base_size = 11) +
	theme(axis.text = element_text(face = "bold"), axis.text.y = element_text(size = 11),
		plot.title = element_text(face = "bold", size = 14),
		legend.position = "bottom", legend.direction = "horizontal")

FigS1 <- pS1A / pS1B + plot_layout(heights = c(2.2, 1))
FigS1
ggsave("FigS1.png", FigS1, width = 12, height = 10, dpi = 600)

save_xlsx("FigS1.data.xlsx", list(
	weekly_count_data = dat.fig1.week2,
	monthly_data = dat.fig1.month,
	weekly_summary = dat.fig1.week2 %>% group_by(dx) %>%
		summarise(weeks = n(), mean_weekly_calls = round(mean(call_count), 1),
			sd_weekly_calls = round(sd(call_count), 1),
			mean_weekly_pct = round(mean(pct_week), 4), sd_weekly_pct = round(sd(pct_week), 4), .groups = "drop"),
	yearly_summary = dat.fig1.week2 %>% group_by(year, dx) %>%
		summarise(year_total_calls = sum(call_count), mean_weekly_calls = round(mean(call_count), 1),
			mean_weekly_pct = round(mean(pct_week), 4), .groups = "drop")
))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 图S2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_dx_trend <- function(dat_list, yrs = 2017:2024, dx_var = "dx_grp", dxs = dxs.grp, time.unit = "weekly", y.unit = "auto", var.group = "dxs", cap = 1000, out_png = NA, width = 12, height = 9, dpi = 300) {
  time.unit <- match.arg(time.unit, c("weekly", "hourly")); y.unit <- if (y.unit == "auto") ifelse(time.unit == "weekly", "count", "pct") else match.arg(y.unit, c("count", "pct"))
  dat1 <- map_dfr(yrs, \(y){ d <- dat_list[[as.character(y)]]; if (is.null(d) || !nrow(d)) return(tibble()); d %>% filter(!is.na(.data[[dx_var]]), !is.na(hour), .data[[dx_var]] %in% dxs, between(hour, 0, 23)) %>% transmute(year = y, hour = as.integer(hour), dx = factor(as.character(.data[[dx_var]]), levels = dxs)) %>% count(year, dx, hour, name = "call_count") %>% complete(year, dx, hour = 0:23, fill = list(call_count = 0)) %>% group_by(year, dx) %>% mutate(pct = call_count / sum(call_count)) %>% ungroup() })
  plots <- lapply(seq_along(dxs), \(i){ dat2 <- dat1 %>% filter(dx == dxs[i]); ggplot(dat2, aes(hour, pct, color = factor(year), group = year)) + geom_line(linewidth = .9) + geom_point(size = 1.6) + scale_color_manual(values = rainbow(length(yrs), s = .8, v = .85), breaks = yrs, name = NULL, drop = FALSE) + scale_x_continuous(breaks = seq(0, 22, 2), labels = sprintf("%02d", seq(0, 22, 2))) + scale_y_continuous(labels = percent_format(accuracy = 1)) + labs(title = dxs[i], x = NULL, y = if (i %% 2 == 1) "Percentage" else NULL) + theme_minimal(base_size = 11) + theme(axis.title = element_text(face = "bold"), axis.text = element_text(face = "bold"), axis.line = element_line(), axis.text.x = element_text(angle = 45, hjust = 1)) })
  p <- wrap_plots(plots, nrow = ceiling(length(dxs)/2), ncol = 2, guides = "collect") & theme(legend.position = "bottom", legend.text = element_text(face = "bold", size = 12)) & guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(linewidth = 2)))
  if (!is.na(out_png)) ggsave(out_png, p, width = width, height = height, dpi = dpi); invisible(p)
}
FigS2 <- plot_dx_trend(dat1.list, yrs = 2017:2024, dx_var = "dx_grp", dxs = dxs.grp, time.unit = "hourly", y.unit = "pct", var.group = "years", out_png = "FigS2.png"); FigS2


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 图S3. Other EMS workflow intervals by phone-score group
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vars_time2 <- c("派车时间", "去程时间", "返程时间")
labs_time2 <- c("Dispatch", "Driving", "Return")

dat1 <- bind_rows(lapply(years, \(y){
	d0 <- dat1.list[[as.character(y)]]
	if (is.null(d0) || !nrow(d0)) return(tibble())
	d0 %>% filter(phone.luck %in% c("low", "high")) %>%
		transmute(
			Year = y,
			phone.luck = as.character(phone.luck),
			Dispatch = suppressWarnings(as.numeric(派车时间) / 60),
			Driving = suppressWarnings(as.numeric(去程时间) / 60),
			Return = suppressWarnings(as.numeric(返程时间) / 60)
		)
}))

dat1 <- dat1 %>%
	pivot_longer(cols = c(Dispatch, Driving, Return), names_to = "interval", values_to = "time_min") %>%
	filter(is.finite(time_min))

dat2 <- dat1 %>%
	group_by(interval, Year, phone.luck) %>%
	summarise(mean = mean(time_min, na.rm = TRUE), .groups = "drop") %>%
	pivot_wider(names_from = phone.luck, values_from = mean) %>%
	mutate(diff = high - low)

FigS3 <- wrap_plots(lapply(c("Dispatch", "Driving", "Return"), \(vv){
	d <- dat2 %>% filter(interval == vv)
	ggplot(d, aes(y = Year)) +
		geom_segment(aes(x = low, xend = high, yend = Year), color = "black", linewidth = .5) +
		geom_point(aes(x = low), color = "grey50", size = 3) +
		geom_point(aes(x = high), color = "blue", size = 3, shape = 17) +
		scale_y_continuous(breaks = years, labels = years) +
		labs(title = vv, x = "Time (mins)", y = "Year") +
		theme_minimal(base_size = 12) +
		theme(
			axis.title = element_text(face = "bold"),
			axis.text = element_text(face = "bold"),
			axis.line = element_line(),
			plot.title = element_text(face = "bold")
		)
}), nrow = 1)

FigS3; ggsave("FigS3.png", FigS3, width = 11, height = 4.2, dpi = 600)

save_xlsx("FigS3.data.xlsx", list(
	panel_data = dat2,
	summary = dat2 %>% group_by(interval) %>% summarise(
		mean_low = round(mean(low, na.rm = TRUE), 2),
		mean_high = round(mean(high, na.rm = TRUE), 2),
		mean_diff = round(mean(diff, na.rm = TRUE), 2),
		.groups = "drop"
	)
))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 图S4. Raw disease-category trends
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_dx <- function(dat_list, years, dx_var, dxs.cn, level = c("raw", "group"), show = c("percent", "count")){
	level <- match.arg(level); show <- match.arg(show)
	dxs_raw <- trimws(as.character(unlist(dxs.cn, use.names = FALSE)))
	map_grp2 <- stack(dxs.cn) %>%
		setNames(c("dx_raw", "group")) %>%
		mutate(across(everything(), \(x) trimws(as.character(x))))

	dat <- purrr::map_dfr(years, \(y){
		d <- dat_list[[as.character(y)]]
		if (is.null(d) || !nrow(d)) return(tibble())

		d0 <- d %>%
			transmute(dx_raw = trimws(as.character(.data[[dx_var]]))) %>%
			filter(!is.na(dx_raw), dx_raw != "", dx_raw != "NA", dx_raw != "NaN")

		out <- if (level == "raw") {
			d0 %>%
				filter(dx_raw %in% dxs_raw) %>%
				count(group = dx_raw, name = "count")
		} else {
			d0 %>%
				left_join(map_grp2, by = "dx_raw") %>%
				count(group, name = "count")
		}

		out %>%
			mutate(group = trimws(as.character(group))) %>%
			filter(!is.na(group), group != "", group != "NA", group != "NaN") %>%
			mutate(year = y, pct = count / sum(count))
	})

	lev <- dat %>%
		filter(year == max(years), !is.na(group), group != "", group != "NA") %>%
		arrange(desc(pct)) %>%
		pull(group) %>%
		unique()

	dat <- dat %>%
		filter(!is.na(group), group != "", group != "NA") %>%
		mutate(group = factor(group, levels = lev), y = if (show == "percent") pct else count) %>%
		filter(!is.na(group))

	p <- ggplot(dat, aes(year, y, color = group, group = group)) +
		geom_line(linewidth = 1) +
		geom_point(size = 2) +
		scale_x_continuous(breaks = years) +
		scale_color_discrete(drop = TRUE, na.translate = FALSE) +
		{if (show == "percent") scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) else scale_y_continuous(expand = c(0, 0))} +
		labs(x = "Year", y = if (show == "percent") "Percentage" else "Count", color = NULL) +
		theme_minimal(base_size = 12)

	list(plot = p, data = dat)
}

tmpS4 <- plot_dx(dat1.list, years, "疾病分类.ML", dxs, level = "raw", show = "percent")
FigS4 <- tmpS4$plot
FigS4; ggsave("FigS4.png", FigS4, width = 11, height = 5.2, dpi = 300)
save_xlsx("FigS4.data.xlsx", tmpS4$data)
