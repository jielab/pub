`%||%` <- function(a, b) if (!is.null(a)) a else b

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 helpers used in main LE8 analysis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
read_proxy_txt <- function(file) {
	if (!file.exists(file)) return(list())
	txt <- readLines(file, warn = FALSE)
	txt <- txt[grepl("^\\s*-\\s*", txt)]
	if (!length(txt)) return(list())
	out <- lapply(txt, function(z) {
		z2 <- sub("^\\s*-\\s*", "", z)
		group <- sub("\\s*\\(.*$", "", z2)
		items <- sub("^.*\\):\\s*", "", z2)
		items <- if (nzchar(items)) strsplit(items, ",")[[1]] else character(0)
		items <- trimws(items)
		items <- items[items != ""]
		list(group = group, items = items)
	})
	names(out) <- sapply(out, `[[`, "group")
	out
}

get_proxy_items <- function(x, group) {
	if (group %in% names(x)) unique(x[[group]]$items) else character(0)
}

get_biom_type <- function(x) {
	dplyr::case_when(
		grepl("^bb_", x) ~ "clinical_chem",
		grepl("^bc_", x) ~ "blood_count",
		TRUE ~ "protein"
	)
}

fcidx2 <- function(y, lp) {
	out <- tryCatch({
		if (exists("fcidx", mode = "function")) fcidx(y, lp) else survival::concordance(y ~ lp)$concordance
	}, error = function(e) NA_real_)
	as.numeric(out)
}

mk_x <- function(dat, vars) {
	vars <- intersect(vars, names(dat))
	if (!length(vars)) return(NULL)
	x <- dat[, vars, drop = FALSE]
	for (j in seq_along(x)) x[[j]] <- suppressWarnings(as.numeric(x[[j]]))
	x <- as.data.frame(x)
	ok <- sapply(x, function(v) sum(is.finite(v)) >= 20 && is.finite(stats::var(v, na.rm = TRUE)) && stats::var(v, na.rm = TRUE) > 0)
	x <- x[, ok, drop = FALSE]
	if (!ncol(x)) return(NULL)
	as.matrix(x)
}

fit_ridge_cox_pred <- function(datTR, datTE, vars_x, t2e.var, event.var, vars.basic = NULL, nfolds.inner = 3) {
	vars0 <- unique(c(vars_x, vars.basic)); vars0 <- intersect(vars0, intersect(names(datTR), names(datTE)))
	if (!length(vars0)) return(rep(NA_real_, nrow(datTE)))
	need <- c(t2e.var, event.var, vars0)
	dtr <- datTR[stats::complete.cases(datTR[, need, drop = FALSE]), need, drop = FALSE]
	dte <- datTE[stats::complete.cases(datTE[, need, drop = FALSE]), need, drop = FALSE]
	if (nrow(dtr) < 50 || nrow(dte) < 5) return(rep(NA_real_, nrow(datTE)))
	xtr <- stats::model.matrix(stats::as.formula(paste0("~", paste(paste0("`", vars0, "`"), collapse = "+"))), dtr)[, -1, drop = FALSE]
	xte <- stats::model.matrix(stats::as.formula(paste0("~", paste(paste0("`", vars0, "`"), collapse = "+"))), dte)[, -1, drop = FALSE]
	ytr <- survival::Surv(dtr[[t2e.var]], dtr[[event.var]])
	fit <- tryCatch(glmnet::cv.glmnet(xtr, ytr, family = "cox", alpha = 0, nfolds = nfolds.inner), error = function(e) NULL)
	if (is.null(fit)) return(rep(NA_real_, nrow(datTE)))
	p <- rep(NA_real_, nrow(datTE)); idx <- seq_len(nrow(dte)); p[idx] <- as.numeric(stats::predict(fit, newx = xte, s = "lambda.min")); p
}

pick_cluster_reps <- function(datTR, vars, k = 23) {
	x <- mk_x(datTR, vars); if (is.null(x)) return(character(0))
	vars <- colnames(x); if (length(vars) <= k) return(vars)
	cm <- suppressWarnings(cor(x, use = "pairwise.complete.obs", method = "spearman")); cm[!is.finite(cm)] <- 0
	d <- as.dist(1 - abs(cm)); hc <- stats::hclust(d, method = "average"); cl <- stats::cutree(hc, k = k)
	out <- lapply(sort(unique(cl)), function(g) { vv <- names(cl)[cl == g]; if (length(vv) == 1) return(vv); subm <- abs(cm[vv, vv, drop = FALSE]); vv[which.max(rowMeans(subm, na.rm = TRUE))] })
	unique(unlist(out))
}

pick_top_univ <- function(datTR, vars, t2e.var, event.var, k = 23) {
	vars <- intersect(vars, names(datTR)); if (!length(vars)) return(character(0))
	score1 <- sapply(vars, function(v) {
		d0 <- datTR[, c(t2e.var, event.var, v), drop = FALSE]; d0 <- d0[stats::complete.cases(d0), , drop = FALSE]
		if (nrow(d0) < 50) return(NA_real_)
		d0[[v]] <- suppressWarnings(as.numeric(d0[[v]]))
		if (!is.finite(var(d0[[v]], na.rm = TRUE)) || var(d0[[v]], na.rm = TRUE) == 0) return(NA_real_)
		fit <- tryCatch(survival::coxph(stats::as.formula(paste0("survival::Surv(", t2e.var, ",", event.var, ") ~ `", v, "`")), data = d0), error = function(e) NULL)
		if (is.null(fit)) return(NA_real_)
		sm <- summary(fit)$coefficients; if (!nrow(sm)) return(NA_real_); abs(sm[1, "z"])
	})
	score1 <- sort(score1, decreasing = TRUE, na.last = NA); head(names(score1), min(k, length(score1)))
}

eval_domain_proxy_cv <- function(dat2, fold_var, biom_vars, yvar, nfolds.inner = 3) {
	biom_vars <- intersect(biom_vars, names(dat2)); if (!length(biom_vars) || !yvar %in% names(dat2)) return(data.frame(cor = NA_real_, r2 = NA_real_))
	res <- lapply(sort(unique(dat2[[fold_var]])), function(fd) {
		dtr <- dat2[dat2[[fold_var]] != fd, , drop = FALSE]; dte <- dat2[dat2[[fold_var]] == fd, , drop = FALSE]
		need <- c(yvar, biom_vars)
		dtr <- dtr[stats::complete.cases(dtr[, need, drop = FALSE]), need, drop = FALSE]; dte <- dte[stats::complete.cases(dte[, need, drop = FALSE]), need, drop = FALSE]
		if (nrow(dtr) < 50 || nrow(dte) < 10) return(data.frame(cor = NA_real_, r2 = NA_real_))
		xtr <- mk_x(dtr, biom_vars); xte <- mk_x(dte, biom_vars); if (is.null(xtr) || is.null(xte)) return(data.frame(cor = NA_real_, r2 = NA_real_))
		common <- intersect(colnames(xtr), colnames(xte)); if (length(common) < 2) return(data.frame(cor = NA_real_, r2 = NA_real_))
		xtr <- xtr[, common, drop = FALSE]; xte <- xte[, common, drop = FALSE]; ytr <- as.numeric(dtr[[yvar]]); yte <- as.numeric(dte[[yvar]])
		fit <- tryCatch(glmnet::cv.glmnet(xtr, ytr, family = "gaussian", alpha = 0, nfolds = nfolds.inner), error = function(e) NULL)
		if (is.null(fit)) return(data.frame(cor = NA_real_, r2 = NA_real_))
		pred <- as.numeric(stats::predict(fit, newx = xte, s = "lambda.min"))
		cor0 <- suppressWarnings(stats::cor(pred, yte, method = "spearman", use = "complete.obs"))
		r20 <- 1 - sum((yte - pred)^2, na.rm = TRUE) / sum((yte - mean(ytr, na.rm = TRUE))^2, na.rm = TRUE)
		data.frame(cor = cor0, r2 = r20)
	})
	do.call(rbind, res)
}

eval_proxy_suite <- function(dat2, Y, Ys.lst, t2e.var, event.var, vars.basic.use0 = NULL, cut_cor = 0.35, nfolds.inner = 3, out_prefix = NULL) {
	if (is.null(out_prefix)) out_prefix <- Y
	file_ns <- paste0(Y, ".proxy.NS.txt"); file_ys <- paste0(Y, ".proxy.YS.txt"); file_ysp <- paste0(Y, ".proxy.YSP.txt")
	if (!all(file.exists(c(file_ns, file_ys, file_ysp)))) { message("skip eval for ", Y, ": proxy txt missing"); return(NULL) }
	if (!"fold_id" %in% names(dat2)) { message("skip eval for ", Y, ": fold_id missing"); return(NULL) }
	px.ns <- read_proxy_txt(file_ns); px.ys <- read_proxy_txt(file_ys); px.ysp <- read_proxy_txt(file_ysp)
	ns <- get_proxy_items(px.ns, "NS.sco"); ys <- unique(unlist(lapply(px.ys, `[[`, "items"))); ysp <- get_proxy_items(px.ysp, "YSP")
	if (!length(ns) || !length(ysp)) { message("skip eval for ", Y, ": empty NS or YSP"); return(NULL) }
	overlap <- intersect(ns, ysp); ns_only <- setdiff(ns, ysp); ysp_only <- setdiff(ysp, ns); k <- length(ysp)

	# membership + anchoring
	domain_names <- intersect(c("diet", "pa", "smoke", "bmi", "nonhdl", "hba1c", "bp", "sleep"), names(px.ysp))
	ysp_annot <- data.frame(biomarker = ysp, in_plus = ysp %in% get_proxy_items(px.ysp, "plus"), stringsAsFactors = FALSE)
	ysp_annot$source_domains <- sapply(ysp, function(b) { hit <- domain_names[sapply(domain_names, function(g) b %in% px.ysp[[g]]$items)]; if (!length(hit)) "" else paste(hit, collapse = ";") })
	ysp_annot$n_source_domains <- ifelse(ysp_annot$source_domains == "", 0, stringr::str_count(ysp_annot$source_domains, ";") + 1)
	save_xlsx(paste0(out_prefix, ".Fig5.out.xlsx"), list(ysp_annotation = ysp_annot))
	write.csv(ysp_annot, paste0(out_prefix, ".ysp_annotation.csv"), row.names = FALSE)

	tab_membership <- data.frame(biomarker = unique(c(ns, ysp)), in_NS = unique(c(ns, ysp)) %in% ns, in_YSP = unique(c(ns, ysp)) %in% ysp, type = get_biom_type(unique(c(ns, ysp))), stringsAsFactors = FALSE)
	tab_membership$class <- ifelse(tab_membership$in_NS & tab_membership$in_YSP, "NS_YSP_overlap", ifelse(tab_membership$in_NS & !tab_membership$in_YSP, "NS_only", "YSP_only"))
	write.csv(tab_membership, paste0(out_prefix, ".proxy_membership.csv"), row.names = FALSE)

	# size + composition
	tab_size <- data.frame(set = c("NS", "NS-only", "Overlap", "YSP"), n = c(length(ns), length(ns_only), length(overlap), length(ysp)))
	p_size <- ggplot2::ggplot(tab_size, ggplot2::aes(set, n)) + ggplot2::geom_col(width = .7, fill = "grey40") + ggplot2::theme_minimal(base_size = 16) + ggplot2::labs(title = paste0(unname(Ys.lst[Y]), ": proxy size reduction"), x = NULL, y = "Number of biomarkers")
	save_plot(p_size, paste0(out_prefix, ".proxy_size_bar.png"), 7, 5)
	tab_comp <- tab_membership %>% dplyr::filter(class %in% c("NS_YSP_overlap", "NS_only")) %>% dplyr::count(class, type, name = "n")
	p_comp <- ggplot2::ggplot(tab_comp, ggplot2::aes(class, n, fill = type)) + ggplot2::geom_col(width = .7) + ggplot2::theme_minimal(base_size = 16) + ggplot2::labs(title = paste0(unname(Ys.lst[Y]), ": composition of selected biomarkers"), x = NULL, y = "Count")
	save_plot(p_comp, paste0(out_prefix, ".proxy_composition.png"), 7, 5)

	# benchmark against simple compact baselines
	folds <- sort(unique(dat2$fold_id))
	bench <- lapply(folds, function(fd) {
		dtr <- dat2[dat2$fold_id != fd, , drop = FALSE]; dte <- dat2[dat2$fold_id == fd, , drop = FALSE]
		ns_cluster23 <- pick_cluster_reps(dtr, ns, k = k); ns_univ23 <- pick_top_univ(dtr, ns, t2e.var, event.var, k = k)
		p_ns67 <- fit_ridge_cox_pred(dtr, dte, ns, t2e.var, event.var, vars.basic = vars.basic.use0, nfolds.inner = nfolds.inner)
		p_ysp <- fit_ridge_cox_pred(dtr, dte, ysp, t2e.var, event.var, vars.basic = vars.basic.use0, nfolds.inner = nfolds.inner)
		p_cl23 <- fit_ridge_cox_pred(dtr, dte, ns_cluster23, t2e.var, event.var, vars.basic = vars.basic.use0, nfolds.inner = nfolds.inner)
		p_uv23 <- fit_ridge_cox_pred(dtr, dte, ns_univ23, t2e.var, event.var, vars.basic = vars.basic.use0, nfolds.inner = nfolds.inner)
		get_c <- function(pred) { idx <- which(is.finite(pred) & !is.na(dte[[t2e.var]]) & !is.na(dte[[event.var]])); if (length(idx) < 10) return(NA_real_); fcidx2(survival::Surv(dte[[t2e.var]][idx], dte[[event.var]][idx]), pred[idx]) }
		data.frame(fold = fd, model = c("NS.full", "YSP", "NS.compact.cluster", "NS.compact.univ"), cidx = c(get_c(p_ns67), get_c(p_ysp), get_c(p_cl23), get_c(p_uv23)), n_marker = c(length(ns), length(ysp), length(ns_cluster23), length(ns_univ23)))
	}) %>% bind_rows()
	bench_sum <- bench %>% dplyr::group_by(model) %>% dplyr::summarise(mean_c = mean(cidx, na.rm = TRUE), sd_c = sd(cidx, na.rm = TRUE), mean_n = mean(n_marker, na.rm = TRUE), .groups = "drop")
	save_xlsx(paste0(out_prefix, ".Fig4.out.xlsx"), list(benchmark_fold = bench, benchmark_summary = bench_sum))
	p_bench <- ggplot2::ggplot(bench, ggplot2::aes(x = model, y = cidx, fill = model)) + ggplot2::geom_violin(alpha = .3, width = .9) + ggplot2::geom_boxplot(width = .12, outlier.shape = NA, alpha = .8) + ggplot2::geom_jitter(width = .08, size = 1.5, alpha = .7) + ggplot2::theme_minimal(base_size = 15) + ggplot2::labs(title = paste0(unname(Ys.lst[Y]), ": benchmark against simple ", k, "-marker baselines"), x = NULL, y = "Test-set C-index") + ggplot2::theme(legend.position = "none")
	save_plot(p_bench, paste0(out_prefix, ".proxy_benchmark.png"), 8, 5.5)

	# domain proxy quality
	tab_domain <- lapply(domain_names, function(g) {
		yvar <- paste0(g, ".pts"); bv <- get_proxy_items(px.ysp, g); tmp <- eval_domain_proxy_cv(dat2, "fold_id", bv, yvar, nfolds.inner = nfolds.inner)
		data.frame(domain = g, n_biom = length(bv), mean_cor = mean(tmp$cor, na.rm = TRUE), mean_r2 = mean(tmp$r2, na.rm = TRUE), median_cor = median(tmp$cor, na.rm = TRUE), median_r2 = median(tmp$r2, na.rm = TRUE))
	}) %>% bind_rows()
	save_xlsx(paste0(out_prefix, ".Fig6.out.xlsx"), list(domain_proxy_quality = tab_domain))
	p_domain <- ggplot2::ggplot(tab_domain, ggplot2::aes(x = reorder(domain, mean_cor), y = mean_cor)) + ggplot2::geom_col(width = .7, fill = "grey40") + ggplot2::coord_flip() + ggplot2::theme_minimal(base_size = 15) + ggplot2::labs(title = paste0(unname(Ys.lst[Y]), ": held-out proxy quality for LE8 domains"), x = NULL, y = "Mean held-out Spearman correlation")
	save_plot(p_domain, paste0(out_prefix, ".domain_proxy_quality.png"), 7, 5.5)

	# network, white background
	vars_net <- intersect(unique(c(ns, ysp)), names(dat2)); sum_node <- data.frame(node_class = c("NS_only", "YSP_kept"), median_degree = NA_real_, median_betweenness = NA_real_)
	if (length(vars_net) >= 3) {
		dnet <- dat2 %>% dplyr::select(dplyr::any_of(vars_net)) %>% dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric))
		cor_mat <- suppressWarnings(stats::cor(dnet, use = "pairwise.complete.obs", method = "spearman")); diag(cor_mat) <- 0
		idx <- which(abs(cor_mat) >= cut_cor, arr.ind = TRUE); idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]
		if (nrow(idx)) {
			edges <- data.frame(from = rownames(cor_mat)[idx[, 1]], to = colnames(cor_mat)[idx[, 2]], rho = cor_mat[idx], stringsAsFactors = FALSE)
			verts <- data.frame(name = unique(c(edges$from, edges$to)), stringsAsFactors = FALSE)
			verts$node_class <- ifelse(verts$name %in% ysp, "YSP_kept", ifelse(verts$name %in% ns_only, "NS_only", "other"))
			g <- igraph::graph_from_data_frame(edges, vertices = verts, directed = FALSE)
			igraph::V(g)$degree <- igraph::degree(g); igraph::V(g)$betweenness <- igraph::betweenness(g, normalized = TRUE)
			tab_node <- data.frame(biomarker = igraph::V(g)$name, degree = igraph::V(g)$degree, betweenness = igraph::V(g)$betweenness, node_class = igraph::V(g)$node_class, stringsAsFactors = FALSE)
			write.csv(tab_node, paste0(out_prefix, ".network_node_metrics.csv"), row.names = FALSE)
			sum_node <- tab_node %>% dplyr::filter(node_class %in% c("YSP_kept", "NS_only")) %>% dplyr::group_by(node_class) %>% dplyr::summarise(median_degree = median(degree, na.rm = TRUE), median_betweenness = median(betweenness, na.rm = TRUE), .groups = "drop")
			p_net <- ggraph::ggraph(g, layout = "fr") +
				ggraph::geom_edge_link(ggplot2::aes(alpha = abs(rho)), colour = "grey75", show.legend = FALSE) +
				ggraph::geom_node_point(ggplot2::aes(size = degree, color = node_class), alpha = .9) +
				ggraph::geom_node_text(ggplot2::aes(label = ifelse(node_class == "YSP_kept", name, "")), repel = TRUE, size = 3) +
				ggplot2::scale_color_manual(values = c("YSP_kept" = "#F8766D", "NS_only" = "#00BFC4", "other" = "grey75")) +
				ggplot2::scale_size_continuous(range = c(2.5, 8)) +
				ggplot2::theme_void(base_size = 14) +
				ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", color = NA), panel.background = ggplot2::element_rect(fill = "white", color = NA), legend.background = ggplot2::element_rect(fill = "white", color = NA), legend.key = ggplot2::element_rect(fill = "white", color = NA), plot.title = ggplot2::element_text(face = "bold")) +
				ggplot2::labs(title = paste0(unname(Ys.lst[Y]), ": NS vs YSP biomarker network"), color = NULL, size = "Degree")
			save_plot(p_net, paste0(out_prefix, ".network_ysp_vs_ns.png"), 10, 8)
		}
	}

	# summary text
	get_bench_mean <- function(md) { x <- bench_sum$mean_c[bench_sum$model == md]; if (!length(x)) NA_real_ else x }
	get_node_val <- function(cls, vn) { x <- sum_node[sum_node$node_class == cls, vn]; if (!length(x)) NA_real_ else x }
	sum_txt <- c(
		paste0(unname(Ys.lst[Y]), " proxy summary"),
		paste0("NS selected ", length(ns), " biomarkers."),
		paste0("YSP selected ", length(ysp), " biomarkers."),
		paste0("Overlap between NS and YSP = ", length(overlap), "."),
		paste0("Removed from NS by YSP = ", length(ns_only), "."),
		paste0("YSP-only biomarkers newly introduced = ", length(ysp_only), "."),
		paste0(sum(ysp_annot$n_source_domains > 0), "/", nrow(ysp_annot), " YSP biomarkers were anchored to at least one LE8 domain proxy set."),
		paste0(sum(ysp_annot$in_plus), "/", nrow(ysp_annot), " YSP biomarkers were also retained in the plus set."),
		paste0("Mean test C-index: NS.full = ", signif(get_bench_mean("NS.full"), 4), ", YSP = ", signif(get_bench_mean("YSP"), 4), ", NS.compact.cluster = ", signif(get_bench_mean("NS.compact.cluster"), 4), ", NS.compact.univ = ", signif(get_bench_mean("NS.compact.univ"), 4), "."),
		paste0("Median network degree in NS-only biomarkers = ", signif(get_node_val("NS_only", "median_degree"), 3), "."),
		paste0("Median network degree in YSP-kept biomarkers = ", signif(get_node_val("YSP_kept", "median_degree"), 3), "."),
		paste0("Median held-out domain proxy correlation across LE8 domains = ", signif(median(tab_domain$mean_cor, na.rm = TRUE), 3), "."),
		"Note: similar proxy counts across LE8 domains are partly driven by the selection cap (pmax_each), so proxy validity should be judged by held-out domain prediction rather than raw counts alone."
	)
	writeLines(sum_txt, paste0(out_prefix, ".summary_for_results.txt")); cat(paste(sum_txt, collapse = "\n"), "\n")
	invisible(list(ns = ns, ys = ys, ysp = ysp, overlap = overlap, ns_only = ns_only, ysp_only = ysp_only, bench = bench, bench_test = bench, bench_sum = bench_sum, sensitivity = data.frame(), tab_domain = tab_domain, ysp_annot = ysp_annot, tab_membership = tab_membership, net = if (exists("tab_node")) tab_node else data.frame(), edges = if (exists("edges")) edges else data.frame(), sum_node = sum_node))
}

proxy_seek <- function(datTR, datTE, vars.biom, mode = c("YS","YSP","NS"), item_names = names.le8, t2e.var = NULL, event.var = NULL, pmax_each = 20, pmax_plus = 30, pf_ys = 0.6, ns_top = NULL, ysp_top_ys = NULL, screen_top_each = 120, allow_overlap_ys = TRUE, SHOW = FALSE, PLOT = FALSE, SAVE = NULL, title = NULL, dat.plot = NULL, ys0 = NULL, nfolds.inner = 3) {
	mode <- match.arg(mode)
	if (is.null(dat.plot)) dat.plot <- datTR
	show_res <- function(obj, title = NULL) {
		txt <- c(if (!is.null(title)) paste0("================ ", title, " ================"), unlist(lapply(names(obj$biom_list), \(nm) sprintf(" - %s (%d): %s", nm, length(obj$biom_list[[nm]]), if (length(obj$biom_list[[nm]])) paste(obj$biom_list[[nm]], collapse = ", ") else "None"))))
		if (SHOW) cat(paste(txt, collapse = "\n"), "\n")
		if (!is.null(SAVE)) writeLines(txt, SAVE)
		if (PLOT) {
			df <- data.frame(Biomarker = unlist(obj$biom_list), Group = rep(names(obj$biom_list), lengths(obj$biom_list)))
			if (nrow(df) >= 3) {
				mat <- scale(as.matrix(dat.plot[, unique(df$Biomarker), drop = FALSE])); mat[is.na(mat)] <- 0
				pc <- prcomp(cor(mat, use = "pairwise.complete.obs"), center = TRUE, scale. = TRUE)
				pdat <- merge(data.frame(Biomarker = rownames(pc$x), PC1 = pc$x[,1], PC2 = pc$x[,2]), df, by = "Biomarker")
				p <- ggplot(pdat, aes(PC1, PC2, color = Group, label = Biomarker)) + geom_point(size = 4, alpha = .8) + ggrepel::geom_text_repel(size = 3, max.overlaps = 20, show.legend = FALSE) + theme_minimal(base_size = 14) + scale_color_brewer(palette = "Set1")
				if (any(table(pdat$Group) >= 3)) p <- p + stat_ellipse(data = subset(pdat, Group %in% names(which(table(pdat$Group) >= 3))), aes(fill = Group), geom = "polygon", alpha = .1, type = "norm", show.legend = FALSE) + scale_fill_brewer(palette = "Set1")
				print(p)
			}
		}
		invisible(obj)
	}
	valid_biom <- intersect(vars.biom, intersect(names(datTR), names(datTE)))
	valid_biom <- valid_biom[sapply(datTR[, valid_biom, drop = FALSE], \(x) { x <- suppressWarnings(as.numeric(x)); sum(is.finite(x)) >= 20 && is.finite(var(x, na.rm = TRUE)) && var(x, na.rm = TRUE) > 0 })]
	if (length(valid_biom) < 2) return(NULL)
	xtr <- scale(as.matrix(datTR[, valid_biom, drop = FALSE])); cen <- attr(xtr, "scaled:center"); scl <- attr(xtr, "scaled:scale"); scl[!is.finite(scl) | scl == 0] <- 1
	xte <- scale(as.matrix(datTE[, valid_biom, drop = FALSE]), center = cen, scale = scl); xtr[is.na(xtr)] <- 0; xte[is.na(xte)] <- 0
	mk_sco <- function(raw_tr, raw_te, nm) {
		m <- mean(raw_tr, na.rm = TRUE); s <- sd(raw_tr, na.rm = TRUE); if (!is.finite(s) || s == 0) s <- 1
		tr <- setNames(data.frame((raw_tr - m) / s), nm); te <- setNames(data.frame((raw_te - m) / s), nm)
		tr[[1]][!is.finite(tr[[1]])] <- 0; te[[1]][!is.finite(te[[1]])] <- 0
		list(tr = tr, te = te)
	}
	.screen_keep <- function(keep, y, n_top = screen_top_each) {
		keep <- intersect(keep, valid_biom)
		if (!length(keep)) return(character(0))
		y <- suppressWarnings(as.numeric(y)); idx <- which(is.finite(y))
		if (length(idx) < 50) return(keep)
		sc <- sapply(keep, function(v) {
			x <- xtr[idx, v]
			if (all(!is.finite(x)) || stats::var(x, na.rm = TRUE) == 0) return(NA_real_)
			abs(suppressWarnings(stats::cor(x, y[idx], use = "pairwise.complete.obs", method = "spearman")))
		})
		sc <- sort(sc, decreasing = TRUE, na.last = NA)
		if (!length(sc)) return(keep)
		head(names(sc), min(length(sc), n_top))
	}
	fit_gaus <- function(keep, y, nm, pmax0 = pmax_each, screen_top = screen_top_each) {
		keep <- .screen_keep(keep, y, screen_top)
		if (length(keep) < 2 || sum(!is.na(y)) < 100) return(list(tr = setNames(data.frame(rep(0, nrow(datTR))), nm), te = setNames(data.frame(rep(0, nrow(datTE))), nm), sel = character(0)))
		idx <- which(!is.na(y))
		fit <- tryCatch(glmnet::cv.glmnet(xtr[idx, keep, drop = FALSE], y[idx], family = "gaussian", alpha = 1, nfolds = nfolds.inner), error = \(e) NULL)
		if (is.null(fit)) return(list(tr = setNames(data.frame(rep(0, nrow(datTR))), nm), te = setNames(data.frame(rep(0, nrow(datTE))), nm), sel = character(0)))
		b <- coef(fit, s = "lambda.min")[-1,,drop = FALSE]; b <- setNames(as.numeric(b), rownames(b)); sel <- intersect(names(b)[b != 0], keep)
		if (!length(sel)) {
			sc <- sapply(keep, function(v) abs(suppressWarnings(stats::cor(xtr[idx, v], y[idx], use = "pairwise.complete.obs", method = "spearman"))))
			sc <- sort(sc, decreasing = TRUE, na.last = NA)
			sel <- head(names(sc), min(length(sc), pmax0))
		}
		if (length(sel) > pmax0) sel <- sel[order(abs(b[sel] %||% 0), decreasing = TRUE)][1:pmax0]
		if (!length(sel)) return(list(tr = setNames(data.frame(rep(0, nrow(datTR))), nm), te = setNames(data.frame(rep(0, nrow(datTE))), nm), sel = character(0)))
		bb <- b[sel]; bb[!is.finite(bb)] <- 0
		if (all(bb == 0)) bb <- rep(1, length(sel))
		tmp <- mk_sco(as.numeric(xtr[, sel, drop = FALSE] %*% bb), as.numeric(xte[, sel, drop = FALSE] %*% bb), nm)
		list(tr = tmp$tr, te = tmp$te, sel = sel)
	}
	fit_cox <- function(keep, nm, pf = NULL, pmax0 = NULL, screen_top = NULL) {
		if (is.null(t2e.var) || is.null(event.var) || !all(c(t2e.var, event.var) %in% names(datTR))) return(NULL)
		if (!is.null(screen_top) && is.finite(screen_top)) {
			idx0 <- which(!is.na(datTR[[t2e.var]]) & !is.na(datTR[[event.var]]))
			if (length(idx0) >= 50) {
				sc <- sapply(intersect(keep, valid_biom), function(v) {
					d0 <- data.frame(t = datTR[[t2e.var]][idx0], y = datTR[[event.var]][idx0], x = xtr[idx0, v])
					if (stats::var(d0$x, na.rm = TRUE) == 0) return(NA_real_)
					fit1 <- tryCatch(survival::coxph(survival::Surv(t, y) ~ x, data = d0), error = function(e) NULL)
					if (is.null(fit1)) return(NA_real_)
					sm <- summary(fit1)$coefficients; if (!nrow(sm)) return(NA_real_); abs(sm[1, "z"])
				})
				sc <- sort(sc, decreasing = TRUE, na.last = NA)
				if (length(sc)) keep <- head(names(sc), min(length(sc), screen_top))
			}
		}
		if (length(keep) < 2) return(NULL)
		idx <- which(!is.na(datTR[[t2e.var]]) & !is.na(datTR[[event.var]])); if (length(idx) < 50) return(NULL)
		fit <- tryCatch(glmnet::cv.glmnet(xtr[idx, keep, drop = FALSE], survival::Surv(datTR[[t2e.var]][idx], datTR[[event.var]][idx]), family = "cox", alpha = 1, nfolds = nfolds.inner, penalty.factor = if (is.null(pf)) rep(1, length(keep)) else pf, pmax = if (is.null(pmax0)) length(keep) else min(pmax0, length(keep))), error = \(e) NULL)
		if (is.null(fit)) return(NULL)
		b <- coef(fit, s = "lambda.min"); b <- setNames(as.numeric(b), rownames(b)); sel <- intersect(names(b)[b != 0], keep)
		if (!is.null(pmax0) && length(sel) > pmax0) sel <- sel[order(abs(b[sel]), decreasing = TRUE)[1:pmax0]]
		tr0 <- setNames(data.frame(rep(0, nrow(datTR))), nm); te0 <- setNames(data.frame(rep(0, nrow(datTE))), nm)
		if (!length(sel)) return(list(tr = tr0, te = te0, sel = character(0), biom_list = setNames(list(character(0)), sub("^biom\\.|\\.sco$", "", nm))))
		tmp <- mk_sco(as.numeric(xtr[, sel, drop = FALSE] %*% b[sel]), as.numeric(xte[, sel, drop = FALSE] %*% b[sel]), nm)
		list(tr = tmp$tr, te = tmp$te, sel = sel, biom_list = setNames(list(sel), sub("^biom\\.|\\.sco$", "", nm)))
	}
	get_ys <- function() {
		outn <- paste0("biom.sco.", item_names, ".biom")
		tr <- as.data.frame(matrix(0, nrow(datTR), length(item_names))); te <- as.data.frame(matrix(0, nrow(datTE), length(item_names)))
		colnames(tr) <- colnames(te) <- outn; biom_list <- list(); avail <- valid_biom
		for (i in seq_along(item_names)) {
			pts <- paste0(item_names[i], ".pts")
			if (!pts %in% names(datTR)) { biom_list[[item_names[i]]] <- character(0); next }
			keep0 <- if (allow_overlap_ys) valid_biom else avail
			obj <- fit_gaus(keep0, datTR[[pts]], outn[i], pmax_each, screen_top_each)
			biom_list[[item_names[i]]] <- obj$sel
			if (length(obj$sel)) {
				if (!allow_overlap_ys) avail <- setdiff(avail, obj$sel)
				tr[, i] <- obj$tr[[1]]; te[, i] <- obj$te[[1]]
			}
		}
		tmp <- mk_sco(rowMeans(tr, na.rm = TRUE), rowMeans(te, na.rm = TRUE), "biom.YS.sco")
		tr$biom.YS.sco <- tmp$tr[[1]]; te$biom.YS.sco <- tmp$te[[1]]
		list(tr = tr[,"biom.YS.sco", drop = FALSE], te = te[,"biom.YS.sco", drop = FALSE], biom_list = biom_list, sel = unique(unlist(biom_list)))
	}
	if (mode == "YS") return(show_res(with(get_ys(), list(tr = tr, te = te, biom_list = biom_list)), if (is.null(title)) "YS proxies" else title))
	if (mode == "NS") {
		keep_ns <- if (is.null(ns_top)) valid_biom else valid_biom
		res0 <- fit_cox(keep_ns, "biom.NS.sco", pmax0 = if (is.null(ns_top)) NULL else min(ns_top, length(keep_ns)), screen_top = ns_top)
		if (is.null(res0)) return(NULL)
		return(show_res(list(tr = res0$tr, te = res0$te, biom_list = res0$biom_list), if (is.null(title)) "NS proxies" else title))
	}
	ys <- if (is.null(ys0)) get_ys() else ys0
	ys_sel <- ys$sel
	if (!is.null(ysp_top_ys) && length(ys_sel) > ysp_top_ys) ys_sel <- ys_sel[seq_len(ysp_top_ys)]
	plus <- fit_cox(setdiff(valid_biom, ys$sel), "biom.plus.sco", pmax0 = pmax_plus, screen_top = max(c(ns_top, pmax_plus), na.rm = TRUE))
	if (is.null(plus)) plus <- list(tr = data.frame(biom.plus.sco = rep(0, nrow(datTR))), te = data.frame(biom.plus.sco = rep(0, nrow(datTE))), sel = character(0))
	union_sel <- unique(c(ys_sel, plus$sel))
	if (!length(union_sel)) return(show_res(list(tr = ys$tr, te = ys$te, biom_list = c(ys$biom_list, list(plus = character(0), YSP = character(0)))), if (is.null(title)) "YSP proxies" else title))
	pf <- ifelse(union_sel %in% ys_sel, pf_ys, 1)
	ysp <- fit_cox(union_sel, "biom.YSP.sco", pf = pf, pmax0 = length(union_sel))
	if (is.null(ysp) || !length(ysp$sel)) ysp <- list(tr = data.frame(biom.YSP.sco = ys$tr$biom.YS.sco), te = data.frame(biom.YSP.sco = ys$te$biom.YS.sco), sel = ys_sel)
	show_res(list(tr = cbind(ys$tr, plus$tr, ysp$tr), te = cbind(ys$te, plus$te, ysp$te), biom_list = c(ys$biom_list, list(plus = plus$sel, YSP = ysp$sel))), if (is.null(title)) "YSP proxies" else title)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 plotting helpers for LE8 manuscript
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.group_cols <- c("LE8-anchored" = "#56B1F7", "Plus-only" = "#E69F00", "NS-only" = "#00BFC4")

fig_theme <- function(base_size = 16) {
	ggplot2::theme_minimal(base_size = base_size) +
		ggplot2::theme(
			plot.title = ggplot2::element_text(face = "bold", size = base_size * 1.05, hjust = 0),
			axis.title = ggplot2::element_text(face = "bold"),
			axis.text = ggplot2::element_text(face = "bold", color = "black"),
			legend.title = ggplot2::element_text(face = "bold"),
			legend.text = ggplot2::element_text(face = "bold"),
			panel.grid.minor = ggplot2::element_blank(),
			plot.margin = ggplot2::margin(12, 16, 12, 16)
		)
}

panel_title <- function(letter, title) ifelse(is.null(letter) || !nzchar(letter), title, paste0(letter, ". ", title))

note_plot <- function(title = NULL, note = "") {
	ggplot2::ggplot() +
		ggplot2::annotate("text", x = 0, y = 0, label = note, size = 5, fontface = "bold") +
		ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1) +
		ggplot2::labs(title = title %||% "") +
		ggplot2::theme_void(base_size = 14) +
		ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0))
}

.proxy_num <- function(x) suppressWarnings(as.numeric(x))
.scale01 <- function(x, to = c(0.35, 1)) {
	x <- as.numeric(x)
	if (!length(x) || all(!is.finite(x))) return(rep(mean(to), length(x)))
	rg <- range(x, na.rm = TRUE)
	if (!is.finite(diff(rg)) || diff(rg) == 0) return(rep(mean(to), length(x)))
	to[1] + (x - rg[1]) / diff(rg) * diff(to)
}
prepare_eval_groups <- function(eval.res, dat2 = NULL, Y = NULL) {
	if (is.null(eval.res)) return(NULL)
	ysp_annot <- eval.res$ysp_annot %||% data.frame()
	tab_membership <- eval.res$tab_membership %||% data.frame()
	net <- eval.res$net %||% data.frame()
	edges <- eval.res$edges %||% data.frame()
	ysp <- if (nrow(ysp_annot)) unique(ysp_annot$biomarker) else character(0)
	plus_only <- if (nrow(ysp_annot) && "in_plus" %in% names(ysp_annot)) unique(ysp_annot$biomarker[ysp_annot$in_plus]) else character(0)
	ns_only <- if (nrow(tab_membership)) unique(tab_membership$biomarker[tab_membership$class == "NS_only"]) else character(0)
	anchored <- setdiff(ysp, plus_only)
	all_selected <- unique(c(anchored, plus_only, ns_only))
	group_df <- data.frame(
		biomarker = all_selected,
		group = dplyr::case_when(all_selected %in% ns_only ~ "NS-only", all_selected %in% plus_only ~ "Plus-only", TRUE ~ "LE8-anchored"),
		stringsAsFactors = FALSE
	)
	group_df$group <- factor(group_df$group, levels = c("LE8-anchored", "Plus-only", "NS-only"))
	domain_names <- c("diet", "pa", "smoke", "bmi", "nonhdl", "hba1c", "bp", "sleep")
	group_df$primary_domain <- NA_character_
	if (nrow(ysp_annot) && "source_domains" %in% names(ysp_annot)) {
		i <- match(group_df$biomarker, ysp_annot$biomarker)
		src <- ysp_annot$source_domains[i]
		group_df$primary_domain <- ifelse(is.na(src) | src == "", NA_character_, vapply(strsplit(src, ";", fixed = TRUE), `[`, character(1), 1))
	}
	group_df$primary_domain[group_df$group == "NS-only"] <- "NS-only"
	domain_cor_all <- NULL
	if (!is.null(dat2) && nrow(group_df)) {
		dom_vars <- paste0(domain_names, ".pts")
		dom_vars <- dom_vars[dom_vars %in% names(dat2)]
		if (length(dom_vars) >= 2) {
			domain_cor_all <- sapply(dom_vars, function(v) {
				sapply(group_df$biomarker, function(b) {
					if (!all(c(b, v) %in% names(dat2))) return(NA_real_)
					d0 <- dat2[, c(b, v), drop = FALSE]
					d0 <- d0[stats::complete.cases(d0), , drop = FALSE]
					if (nrow(d0) < 30) return(NA_real_)
					suppressWarnings(stats::cor(.proxy_num(d0[[b]]), .proxy_num(d0[[v]]), method = "spearman"))
				})
			})
			if (is.null(dim(domain_cor_all))) domain_cor_all <- matrix(domain_cor_all, ncol = length(dom_vars), dimnames = list(group_df$biomarker, gsub("\\.pts$", "", dom_vars)))
			else {
				colnames(domain_cor_all) <- gsub("\\.pts$", "", dom_vars)
				rownames(domain_cor_all) <- group_df$biomarker
			}
		}
	}
	net2 <- if (nrow(net)) dplyr::left_join(net, group_df, by = "biomarker") else NULL
	if (!is.null(net2) && nrow(net2)) net2$group <- factor(as.character(net2$group), levels = names(.group_cols))
	list(group_df = group_df, domain_cor_all = domain_cor_all, net2 = net2, edges = edges)
}

plot_perf <- function(ml.res, Ylab = NULL) {
	target <- intersect(c("M.1sco", "M.8each", "M.raw.only.bp", "M.biom.YS", "M.biom.NS", "M.biom.YSP"), unique(as.character(ml.res$model)))
	d <- as.data.frame(ml.res)[, intersect(c("fold", "model", "cidx"), names(as.data.frame(ml.res))), drop = FALSE]
	d <- d[is.finite(d$cidx) & as.character(d$model) %in% target, , drop = FALSE]
	if (!nrow(d)) return(note_plot(paste0("Model performance", if (!is.null(Ylab)) paste0(": ", Ylab) else ""), "No model results"))
	ord <- d %>% dplyr::group_by(model) %>% dplyr::summarise(med = stats::median(cidx, na.rm = TRUE), .groups = "drop") %>% dplyr::arrange(med)
	d$model <- factor(as.character(d$model), levels = ord$model)
	pal <- c("M.1sco" = "#F8766D", "M.8each" = "#00BF7D", "M.raw.only.bp" = "#C77CFF", "M.biom.YS" = "#A3A500", "M.biom.YSP" = "#00B0F6", "M.biom.NS" = "#E76BF3")
	ggplot2::ggplot(d, ggplot2::aes(cidx, model, fill = model)) +
		ggplot2::geom_violin(alpha = 0.35, linewidth = 0.8, color = "black", trim = FALSE) +
		ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.95, linewidth = 0.8) +
		ggplot2::geom_point(position = ggplot2::position_jitter(height = 0.08, width = 0), size = 2.1, alpha = 0.7, color = "black") +
		ggplot2::scale_fill_manual(values = pal, guide = "none") +
		ggplot2::labs(title = paste0("Predictive performance", if (!is.null(Ylab)) paste0(": ", Ylab) else ""), x = "C-index", y = NULL) +
		fig_theme(16)
}


strip_plot_title <- function(p) p + ggplot2::labs(title = NULL)

collect_group_proteins <- function(eval.res, dat2 = NULL) {
	obj <- prepare_eval_groups(eval.res, dat2)
	if (is.null(obj) || !nrow(obj$group_df)) return(list())
	setNames(lapply(names(.group_cols), function(g) obj$group_df$biomarker[obj$group_df$group == g]), names(.group_cols))
}

run_trend_groups <- function(dat1, eval.res, t2e.var, event.var, ylab = NULL, top_n = 18, year_max = 14.4) {
	gps <- collect_group_proteins(eval.res, dat1)
	plot_lst <- list(); out <- list()
	xlab <- paste0("Years before ", ifelse(is.null(ylab), "event", tolower(ylab)), " diagnosis")
	for (g in names(gps)) {
		prots <- intersect(gps[[g]], names(dat1))
		if (length(prots) < 2) {
			plot_lst[[g]] <- note_plot(g, "<2 proteins")
			out[[paste0(g, ".loess")]] <- data.frame()
			out[[paste0(g, ".cluster")]] <- data.frame()
			out[[paste0(g, ".zdat")]] <- data.frame()
			next
		}
		prots <- head(prots, top_n)
		tr <- plot_trend(dat = dat1, proteins = prots, t2e.var = t2e.var, event.var = event.var, year_max = year_max, top_order = prots)
		plot_lst[[g]] <- strip_plot_title(tr$p_bubble + ggplot2::labs(x = xlab, y = NULL, subtitle = g))
		out[[paste0(g, ".loess")]] <- tr$loess_df %||% data.frame()
		out[[paste0(g, ".cluster")]] <- tr$cluster_df %||% data.frame()
		out[[paste0(g, ".zdat")]] <- tr$zdat %||% data.frame()
	}
	list(plot = patchwork::wrap_plots(plot_lst, nrow = 1), out = out)
}

run_pca_distance <- function(eval.res, dat1) {
	obj <- prepare_eval_groups(eval.res, dat1)
	if (is.null(obj) || is.null(obj$domain_cor_all)) return(NULL)
	mat <- obj$domain_cor_all
	mat <- mat[stats::complete.cases(mat), , drop = FALSE]
	if (nrow(mat) < 3 || ncol(mat) < 2) return(NULL)
	pc <- stats::prcomp(mat, center = TRUE, scale. = TRUE)
	dp <- as.data.frame(pc$x[, 1:2, drop = FALSE]); names(dp) <- c("PC1", "PC2")
	dp$biomarker <- rownames(dp)
	dp <- dplyr::left_join(dp, obj$group_df, by = "biomarker")
	anc <- dp[dp$group == "LE8-anchored", c("PC1", "PC2"), drop = FALSE]
	if (nrow(anc)) {
		ctr <- c(mean(anc$PC1, na.rm = TRUE), mean(anc$PC2, na.rm = TRUE))
		dp$dist_to_le8 <- sqrt((dp$PC1 - ctr[1])^2 + (dp$PC2 - ctr[2])^2)
	} else dp$dist_to_le8 <- NA_real_
	list(dp = dp, obj = obj, pc = pc)
}

plot_pca_space <- function(eval.res, dat1, panel = NULL, ylab = NULL) {
	x <- run_pca_distance(eval.res, dat1)
	if (is.null(x)) return(note_plot(if (!is.null(panel)) paste0(panel, ". ", ylab %||% "") else ylab, "No PCA data"))
	dp <- x$dp; pc <- x$pc
	ggplot2::ggplot(dp, ggplot2::aes(PC1, PC2, color = group)) +
		ggplot2::geom_point(size = 3.2, alpha = 0.9) +
		ggplot2::stat_ellipse(ggplot2::aes(group = group), linewidth = 0.9, alpha = 0.22, show.legend = FALSE, na.rm = TRUE) +
		ggrepel::geom_text_repel(ggplot2::aes(label = biomarker), size = 3, max.overlaps = 25, show.legend = FALSE) +
		ggplot2::scale_color_manual(values = .group_cols, drop = FALSE) +
		ggplot2::labs(title = if (!is.null(panel)) panel_title(panel, ylab %||% "") else NULL, x = paste0("PC1 (", sprintf("%.1f", 100 * summary(pc)$importance[2, 1]), "%)"), y = paste0("PC2 (", sprintf("%.1f", 100 * summary(pc)$importance[2, 2]), "%)"), color = NULL) +
		fig_theme(16)
}

plot_pca_distance <- function(eval.res, dat1, panel = NULL, ylab = NULL) {
	x <- run_pca_distance(eval.res, dat1)
	if (is.null(x)) return(note_plot(if (!is.null(panel)) paste0(panel, ". ", ylab %||% "") else ylab, "No PCA data"))
	dp <- x$dp
	ggplot2::ggplot(dp, ggplot2::aes(group, dist_to_le8, fill = group)) +
		ggplot2::geom_violin(alpha = 0.28, color = NA, trim = FALSE) +
		ggplot2::geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.9, linewidth = 0.8) +
		ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.08, height = 0), size = 2, alpha = 0.7) +
		ggplot2::scale_fill_manual(values = .group_cols, drop = FALSE, guide = "none") +
		ggplot2::labs(title = if (!is.null(panel)) panel_title(panel, ylab %||% "") else NULL, x = NULL, y = "Euclidean distance in PCA space") +
		fig_theme(16)
}

plot_network_graph <- function(eval.res, panel = NULL, ylab = NULL, top_label = 25) {
	obj <- prepare_eval_groups(eval.res)
	net2 <- obj$net2; edges <- obj$edges
	if (is.null(net2) || !nrow(net2) || is.null(edges) || !nrow(edges)) return(note_plot(if (!is.null(panel)) paste0(panel, ". ", ylab %||% "") else ylab, "No network"))
	verts <- net2 %>% dplyr::rename(name = biomarker)
	g <- igraph::graph_from_data_frame(edges, vertices = verts, directed = FALSE)
	lab <- names(sort(igraph::degree(g), decreasing = TRUE))[seq_len(min(top_label, igraph::gorder(g)))]
	ggraph::ggraph(g, layout = "fr") +
		ggraph::geom_edge_link(ggplot2::aes(alpha = abs(rho)), colour = "grey75", show.legend = FALSE, linewidth = 0.5) +
		ggraph::geom_node_point(ggplot2::aes(size = degree, color = group), alpha = 0.92) +
		ggraph::geom_node_text(ggplot2::aes(label = ifelse(name %in% lab, name, "")), repel = TRUE, size = 3, show.legend = FALSE) +
		ggplot2::scale_color_manual(values = .group_cols, drop = FALSE) +
		ggplot2::scale_size_continuous(range = c(2.5, 8)) +
		ggplot2::labs(title = if (!is.null(panel)) panel_title(panel, ylab %||% "") else NULL, color = NULL, size = "Degree") +
		ggplot2::theme_void(base_size = 14) +
		ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0), legend.position = "bottom")
}

run_reactome_groups <- function(eval.res, orgdb = org.Hs.eg.db, p_cutoff = 0.2, q_cutoff = 0.2, min_g = 3) {
	obj <- prepare_eval_groups(eval.res)
	if (is.null(obj) || !nrow(obj$group_df)) return(list(tab = data.frame(), map = data.frame()))
	conv <- function(xs) {
		xs <- unique(xs)
		xs <- xs[!grepl("^(bb_|bc_)", xs)]
		xs <- xs[grepl("^[A-Za-z0-9._-]+$", xs)]
		xs
	}
	maps <- list(); tabs <- list()
	for (g in names(.group_cols)) {
		genes <- conv(obj$group_df$biomarker[obj$group_df$group == g])
		if (length(genes) < min_g) next
		map <- tryCatch(clusterProfiler::bitr(genes, fromType = "SYMBOL", toType = c("ENTREZID", "SYMBOL"), OrgDb = orgdb), error = function(e) NULL)
		if (is.null(map) || !nrow(map)) next
		map <- unique(map[, c("SYMBOL", "ENTREZID"), drop = FALSE])
		enr <- tryCatch(ReactomePA::enrichPathway(gene = unique(map$ENTREZID), organism = "human", pvalueCutoff = p_cutoff, qvalueCutoff = q_cutoff, readable = TRUE), error = function(e) NULL)
		if (is.null(enr)) next
		tab <- as.data.frame(enr)
		if (!nrow(tab)) next
		tab$group <- g
		tab$gene_n <- length(unique(map$SYMBOL))
		tabs[[g]] <- tab
		maps[[g]] <- transform(map, group = g)
	}
	list(tab = dplyr::bind_rows(tabs), map = dplyr::bind_rows(maps))
}

plot_reactome_integrated <- function(eval.res, panel = NULL, ylab = NULL, top_n = 8) {
	rx <- run_reactome_groups(eval.res)
	d <- rx$tab
	if (!nrow(d)) return(note_plot(if (!is.null(panel)) paste0(panel, ". ", ylab %||% "") else ylab, "No Reactome terms"))
	d <- d %>% dplyr::filter(is.finite(p.adjust)) %>% dplyr::group_by(group) %>% dplyr::arrange(p.adjust, .by_group = TRUE) %>% dplyr::slice_head(n = top_n) %>% dplyr::ungroup() %>% dplyr::mutate(score = -log10(p.adjust), Description = stringr::str_trunc(Description, 55))
	d$Description <- factor(d$Description, levels = rev(unique(d$Description)))
	ggplot2::ggplot(d, ggplot2::aes(score, Description, fill = group)) +
		ggplot2::geom_col() +
		ggplot2::facet_wrap(~group, scales = "free_y") +
		ggplot2::scale_fill_manual(values = .group_cols, drop = FALSE, guide = "none") +
		ggplot2::labs(title = if (!is.null(panel)) panel_title(panel, ylab %||% "") else NULL, x = expression(-log[10](adjusted~p)), y = NULL) +
		fig_theme(15)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Figure 5
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.reactome_name_col <- function(d) intersect(c("Description", "description", "term", "pathway", "ID"), names(d))[1]
.reactome_p_col <- function(d) intersect(c("p.adjust", "padj", "FDR", "qvalue", "p_adj"), names(d))[1]

.collect_reactome_long <- function(eval.res) {
	if (is.null(eval.res) || !length(eval.res)) return(data.frame())

	parse_set_from_name <- function(nm) {
		nm2 <- tolower(nm)
		if (grepl("anch", nm2)) return("Anchored")
		if (grepl("plus", nm2)) return("Plus-only")
		if (grepl("ns", nm2)) return("NS-only")
		for (x in c("diet","pa","smoke","bmi","nonhdl","hba1c","bp","sleep")) if (grepl(x, nm2)) return(x)
		NA_character_
	}

	out <- lapply(names(eval.res), function(nm) {
		x <- eval.res[[nm]]
		if (!is.data.frame(x) || !nrow(x)) return(NULL)
		nc <- .reactome_name_col(x); pc <- .reactome_p_col(x)
		if (is.na(nc) || is.na(pc)) return(NULL)

		set_col <- intersect(c("set","group","component","source_set","source_domain"), names(x))[1]
		set0 <- if (!is.na(set_col)) as.character(x[[set_col]]) else parse_set_from_name(nm)

		d <- data.frame(
			pathway = as.character(x[[nc]]),
			padj = suppressWarnings(as.numeric(x[[pc]])),
			set = set0,
			source_object = nm,
			stringsAsFactors = FALSE
		)
		d <- d[is.finite(d$padj) & !is.na(d$set) & d$pathway != "", , drop = FALSE]
		if (!nrow(d)) return(NULL)
		d
	})
	d <- dplyr::bind_rows(out)
	if (!nrow(d)) return(data.frame())
	d$score <- -log10(pmax(d$padj, 1e-300))
	d
}

plot_reactome_integrated_heatmap <- function(eval.res, Ylab = NULL, top_n = 26) {
	d <- .collect_reactome_long(eval.res)
	if (!nrow(d)) return(note_plot("Integrated Reactome enrichment", "No Reactome tables found in eval.res"))

	col_order <- c("Anchored","Plus-only","NS-only","diet","pa","smoke","bmi","nonhdl","hba1c","bp","sleep")
	d <- d[d$set %in% col_order, , drop = FALSE]
	if (!nrow(d)) return(note_plot("Integrated Reactome enrichment", "No usable Reactome sets"))

	mat0 <- d %>%
		dplyr::group_by(pathway, set) %>%
		dplyr::summarise(score = max(score, na.rm = TRUE), .groups = "drop") %>%
		tidyr::pivot_wider(names_from = set, values_from = score)

	mat0 <- as.data.frame(mat0)
	rownames(mat0) <- mat0$pathway
	mat0$pathway <- NULL

	for (cc in setdiff(col_order, colnames(mat0))) mat0[[cc]] <- NA_real_
	mat0 <- mat0[, col_order[col_order %in% colnames(mat0)], drop = FALSE]

	if (!nrow(mat0) || !ncol(mat0)) return(note_plot("Integrated Reactome enrichment", "Empty heatmap matrix"))

	biom_cols <- intersect(c("Anchored","Plus-only","NS-only"), colnames(mat0))
	top_path <- rownames(mat0)[order(apply(mat0[, biom_cols, drop = FALSE], 1, function(z) max(z, na.rm = TRUE)), decreasing = TRUE)]
	top_path <- top_path[seq_len(min(top_n, length(top_path)))]
	mat1 <- mat0[top_path, , drop = FALSE]

	dom_group <- apply(mat1[, biom_cols, drop = FALSE], 1, function(z) {
		if (all(!is.finite(z))) return(NA_character_)
		names(which.max(z))[1]
	})
	row_anno <- data.frame(dominant_group = factor(dom_group, levels = c("Anchored","Plus-only","NS-only")))
	rownames(row_anno) <- rownames(mat1)

	col_type <- ifelse(colnames(mat1) %in% c("Anchored","Plus-only","NS-only"), "Biomarker group", "LE8 component")
	col_anno <- data.frame(type = factor(col_type, levels = c("Biomarker group","LE8 component")))
	rownames(col_anno) <- colnames(mat1)

	ann_colors <- list(
		dominant_group = c("Anchored" = "#F4A3A3", "Plus-only" = "#E78AC3", "NS-only" = "#00BFC4"),
		type = c("Biomarker group" = "#B8B000", "LE8 component" = "#00BFFF")
	)

	pheatmap::pheatmap(
		mat = as.matrix(mat1),
		color = colorRampPalette(c("#D9ECF2", "#F1C27D", "#D7191C"))(100),
		cluster_rows = TRUE,
		cluster_cols = FALSE,
		scale = "none",
		border_color = NA,
		annotation_row = row_anno,
		annotation_col = col_anno,
		annotation_colors = ann_colors,
		main = paste0("Integrated Reactome enrichment", if (!is.null(Ylab)) paste0(": ", Ylab) else ""),
		angle_col = 45,
		cellwidth = 26,
		cellheight = 18,
		fontsize = 10,
		fontsize_row = 9,
		fontsize_col = 10,
		silent = TRUE
	)
}

get_reactome_integrated_out <- function(eval.res, top_n = 26) {
	d <- .collect_reactome_long(eval.res)
	if (!nrow(d)) return(list(reactome_long = data.frame(), reactome_matrix = data.frame(), row_annotation = data.frame(), col_annotation = data.frame()))

	col_order <- c("Anchored","Plus-only","NS-only","diet","pa","smoke","bmi","nonhdl","hba1c","bp","sleep")
	d <- d[d$set %in% col_order, , drop = FALSE]

	mat0 <- d %>%
		dplyr::group_by(pathway, set) %>%
		dplyr::summarise(score = max(score, na.rm = TRUE), .groups = "drop") %>%
		tidyr::pivot_wider(names_from = set, values_from = score)

	mat0 <- as.data.frame(mat0)
	rownames(mat0) <- mat0$pathway
	mat0$pathway <- NULL
	for (cc in setdiff(col_order, colnames(mat0))) mat0[[cc]] <- NA_real_
	mat0 <- mat0[, col_order[col_order %in% colnames(mat0)], drop = FALSE]

	biom_cols <- intersect(c("Anchored","Plus-only","NS-only"), colnames(mat0))
	top_path <- rownames(mat0)[order(apply(mat0[, biom_cols, drop = FALSE], 1, function(z) max(z, na.rm = TRUE)), decreasing = TRUE)]
	top_path <- top_path[seq_len(min(top_n, length(top_path)))]
	mat1 <- mat0[top_path, , drop = FALSE]

	row_anno <- data.frame(
		pathway = rownames(mat1),
		dominant_group = apply(mat1[, biom_cols, drop = FALSE], 1, function(z) if (all(!is.finite(z))) NA_character_ else names(which.max(z))[1]),
		row.names = NULL
	)
	col_anno <- data.frame(
		column = colnames(mat1),
		type = ifelse(colnames(mat1) %in% c("Anchored","Plus-only","NS-only"), "Biomarker group", "LE8 component"),
		row.names = NULL
	)

	list(
		reactome_long = d,
		reactome_matrix = cbind(pathway = rownames(mat1), as.data.frame(mat1)),
		row_annotation = row_anno,
		col_annotation = col_anno
	)
}











# backward-compatible alias
# plot_reactome_groups <- plot_reactome_integrated
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 likely-unused functions kept at the end for manual deletion
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plot_delta <- function(eval.res, Ylab = NULL) {
	if (is.null(eval.res$bench) || !nrow(eval.res$bench)) return(note_plot("Fold-wise delta C-index", "No paired comparison"))
	dd <- eval.res$bench %>%
		dplyr::filter(model %in% c("YSP", "NS.compact.cluster", "NS.compact.univ", "NS.full")) %>%
		dplyr::select(fold, model, cidx) %>%
		tidyr::pivot_wider(names_from = model, values_from = cidx)
	if (!all(c("YSP", "NS.compact.cluster", "NS.compact.univ", "NS.full") %in% names(dd))) return(note_plot("Fold-wise delta C-index", "Comparison models unavailable"))
	dd <- dd %>%
		dplyr::mutate(
			`YSP - NS.compact.cluster` = YSP - `NS.compact.cluster`,
			`YSP - NS.compact.univ` = YSP - `NS.compact.univ`,
			`YSP - NS.full` = YSP - `NS.full`
		) %>%
		dplyr::select(fold, starts_with("YSP - ")) %>%
		tidyr::pivot_longer(-fold, names_to = "compare", values_to = "delta")

	ggplot2::ggplot(dd, ggplot2::aes(compare, delta, fill = compare)) +
		ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "grey45", linewidth = 0.9) +
		ggplot2::geom_boxplot(alpha = 0.8, width = 0.72, outlier.shape = NA, linewidth = 0.9) +
		ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.08, height = 0), size = 2.4, alpha = 0.75) +
		ggplot2::stat_summary(fun = mean, geom = "point", shape = 23, size = 4, fill = "yellow", color = "black") +
		ggplot2::scale_fill_manual(values = c("#F8766D", "#00BA38", "#619CFF"), guide = "none") +
		ggplot2::labs(title = panel_title("B", paste0("YSP versus agnostic baselines", if (!is.null(Ylab)) paste0(": ", Ylab) else "")), x = NULL, y = expression(Delta*"C-index")) +
		fig_theme(17) +
		ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 18, hjust = 1))
}

plot_sensitivity <- function(eval.res, Ylab = NULL) {
	d <- eval.res$sensitivity %||% data.frame()
	if (!nrow(d)) return(note_plot("Sensitivity analyses", "No sensitivity analysis"))
	d$model <- factor(d$model, levels = unique(d$model))
	ggplot2::ggplot(d, ggplot2::aes(model, cidx, fill = model)) +
		ggplot2::geom_boxplot(alpha = 0.8, width = 0.6, outlier.shape = NA, linewidth = 0.9) +
		ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.07, height = 0), size = 2.6, alpha = 0.8) +
		ggplot2::labs(title = panel_title("C", paste0("Sensitivity analyses", if (!is.null(Ylab)) paste0(": ", Ylab) else "")), x = NULL, y = "Test-set C-index") +
		fig_theme(17) +
		ggplot2::theme(legend.position = "none")
}

plot_component_cor <- function(dat, Ylab = NULL) {
	if (!all(vars.le8 %in% names(dat))) return(note_plot("LE8 component correlation", "Required LE8 variables are missing"))
	cm <- stats::cor(dplyr::select(dat, tidyselect::all_of(vars.le8)), use = "pairwise.complete.obs")
	as.data.frame(as.table(cm)) %>%
		dplyr::rename(var1 = Var1, var2 = Var2, cor = Freq) %>%
		ggplot2::ggplot(ggplot2::aes(var2, var1, fill = cor)) +
		ggplot2::geom_tile(color = "white") +
		ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", cor)), fontface = "bold", size = 4.2) +
		ggplot2::scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0, limits = c(-1, 1)) +
		ggplot2::coord_equal() +
		ggplot2::labs(title = paste0("LE8 component correlation", if (!is.null(Ylab)) paste0(": ", Ylab) else ""), x = NULL, y = NULL, fill = "r") +
		fig_theme(16) +
		ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 0))
}

plot_domain_quality <- function(eval.res, Ylab = NULL) {
	d <- eval.res$tab_domain %||% data.frame()
	if (!nrow(d)) return(note_plot("Held-out quality of LE8 domain proxies", "No data"))
	d <- d[is.finite(d$mean_cor) & is.finite(d$n_biom), , drop = FALSE]
	if (!nrow(d)) return(note_plot("Held-out quality of LE8 domain proxies", "No valid data"))
	ord <- c("sleep", "pa", "diet", "smoke", "hba1c", "bp", "bmi", "nonhdl")
	d$domain <- factor(d$domain, levels = ord[ord %in% d$domain])
	ggplot2::ggplot(d, ggplot2::aes(mean_cor, domain, color = mean_r2, size = n_biom)) +
		ggplot2::geom_segment(ggplot2::aes(x = 0, xend = mean_cor, yend = domain), linewidth = 1.2, alpha = 0.7, color = "grey70") +
		ggplot2::geom_point() +
		ggplot2::scale_color_gradient(low = "#ABD9E9", high = "#D7191C") +
		ggplot2::labs(title = panel_title("A", paste0("Quality of LE8 domain proxies", if (!is.null(Ylab)) paste0(": ", Ylab) else "")), x = "Held-out Spearman correlation", y = NULL, color = expression(R^2), size = "N proteins") +
		fig_theme(17)
}

plot_ysp_pca <- function(eval.res, dat2, Ylab = NULL) {
	obj <- prepare_eval_groups(eval.res, dat2)
	if (is.null(obj) || is.null(obj$domain_cor_all)) return(note_plot(panel_title("B", "PCA of selected biomarkers"), "No domain-correlation matrix"))
	mat <- obj$domain_cor_all
	mat <- mat[stats::complete.cases(mat), , drop = FALSE]
	if (nrow(mat) < 3 || ncol(mat) < 2) return(note_plot(panel_title("B", "PCA of selected biomarkers"), "No PCA data"))
	pc <- stats::prcomp(mat, center = TRUE, scale. = TRUE)
	dp <- as.data.frame(pc$x[, 1:2, drop = FALSE]); names(dp) <- c("PC1", "PC2")
	dp$biomarker <- rownames(dp)
	dp <- dplyr::left_join(dp, obj$group_df, by = "biomarker")

	ggplot2::ggplot(dp, ggplot2::aes(PC1, PC2, color = group)) +
		ggplot2::geom_point(size = 3, alpha = 0.9) +
		ggplot2::stat_ellipse(ggplot2::aes(group = group), linewidth = 0.9, alpha = 0.25, show.legend = FALSE) +
		ggplot2::scale_color_manual(values = .group_cols, drop = FALSE) +
		ggrepel::geom_text_repel(ggplot2::aes(label = biomarker), size = 3, max.overlaps = 20, show.legend = FALSE) +
		ggplot2::labs(
			title = panel_title("B", paste0("PCA of selected biomarkers", if (!is.null(Ylab)) paste0(": ", Ylab) else "")),
			x = paste0("PC1 (", sprintf("%.1f", 100 * summary(pc)$importance[2, 1]), "%)"),
			y = paste0("PC2 (", sprintf("%.1f", 100 * summary(pc)$importance[2, 2]), "%)")
		) +
		fig_theme(17)
}

plot_fig2c_corr_heatmap_signed <- function(eval.res, dat2, Ylab = NULL) {
	obj <- prepare_eval_groups(eval.res, dat2)
	mat <- obj$domain_cor_all
	if (is.null(mat) || !nrow(mat)) return(note_plot(panel_title("C", "Signed domain-correlation heatmap"), "No matrix"))
	mat <- mat[stats::complete.cases(mat), , drop = FALSE]
	if (!nrow(mat)) return(note_plot(panel_title("C", "Signed domain-correlation heatmap"), "No complete cases"))
	dd <- as.data.frame(as.table(mat))
	names(dd) <- c("biomarker", "domain", "cor")
	ann <- obj$group_df[, c("biomarker", "group"), drop = FALSE]
	dd <- dplyr::left_join(dd, ann, by = "biomarker")
	ord <- ann %>% dplyr::arrange(group, biomarker) %>% pull(biomarker)
	dd$biomarker <- factor(dd$biomarker, levels = ord)
	ggplot2::ggplot(dd, ggplot2::aes(domain, biomarker, fill = cor)) +
		ggplot2::geom_tile(color = "white") +
		ggplot2::scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0, limits = c(-1, 1)) +
		ggplot2::labs(title = panel_title("C", paste0("Signed LE8-domain correlations", if (!is.null(Ylab)) paste0(": ", Ylab) else "")), x = NULL, y = NULL, fill = "rho") +
		fig_theme(16) +
		ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
}

plot_top_ysp <- function(eval.res, Ylab = NULL, top_n = 15) {
	obj <- prepare_eval_groups(eval.res)
	net2 <- obj$net2
	if (is.null(net2) || !nrow(net2)) return(note_plot(panel_title("A", "Top biomarkers by network centrality"), "No network"))
	if (!"degree" %in% names(net2)) return(note_plot(panel_title("A", "Top biomarkers by network centrality"), "No centrality metrics"))
	d <- net2 %>%
		dplyr::mutate(degree = .proxy_num(degree), betweenness = .proxy_num(betweenness)) %>%
		dplyr::arrange(dplyr::desc(degree), dplyr::desc(betweenness), biomarker) %>%
		dplyr::slice_head(n = top_n)
	d$biomarker <- factor(d$biomarker, levels = rev(d$biomarker))
	ggplot2::ggplot(d, ggplot2::aes(degree, biomarker, fill = group)) +
		ggplot2::geom_col() +
		ggplot2::scale_fill_manual(values = .group_cols, drop = FALSE) +
		ggplot2::labs(title = paste0("Top biomarkers by network centrality", if (!is.null(Ylab)) paste0(": ", Ylab) else ""), x = "Degree", y = NULL, fill = NULL) +
		fig_theme(17)
}

plot_network_metrics <- function(eval.res, Ylab = NULL) {
	obj <- prepare_eval_groups(eval.res)
	net2 <- obj$net2
	edges <- obj$edges
	if (is.null(net2) || !nrow(net2) || is.null(edges) || !nrow(edges)) return(list(NULL, note_plot(panel_title("B", "Connectivity strength of selected biomarkers"), "No network")))
	edge_long <- edges %>%
		dplyr::transmute(biomarker = from, abs_rho = abs(rho)) %>%
		dplyr::bind_rows(edges %>% dplyr::transmute(biomarker = to, abs_rho = abs(rho)))
	d <- net2 %>%
		dplyr::left_join(edge_long %>% dplyr::group_by(biomarker) %>% dplyr::summarise(mean_abs_rho = mean(abs_rho, na.rm = TRUE), .groups = "drop"), by = "biomarker") %>%
		dplyr::mutate(mean_abs_rho = tidyr::replace_na(mean_abs_rho, 0), degree = .proxy_num(degree))
	d$biomarker <- factor(d$biomarker, levels = d %>% dplyr::arrange(degree, mean_abs_rho) %>% pull(biomarker))
	p <- ggplot2::ggplot(d, ggplot2::aes(mean_abs_rho, biomarker, color = group, size = degree)) +
		ggplot2::geom_point(alpha = 0.9) +
		ggplot2::scale_color_manual(values = .group_cols, drop = FALSE) +
		ggplot2::labs(title = paste0("Connectivity strength of selected biomarkers", if (!is.null(Ylab)) paste0(": ", Ylab) else ""), x = "Mean absolute correlation in network", y = NULL, color = NULL, size = "Degree") +
		fig_theme(16)
	list(NULL, p)
}

plot_fig3c_dist <- function(eval.res, Ylab = NULL) {
	obj <- prepare_eval_groups(eval.res)
	d <- obj$group_df
	if (is.null(d) || !nrow(d)) return(note_plot(panel_title("C", "Distribution of selected biomarkers"), "No biomarker groups"))
	d <- d %>% dplyr::count(group, primary_domain, name = "n")
	ggplot2::ggplot(d, ggplot2::aes(primary_domain, n, fill = group)) +
		ggplot2::geom_col(position = "stack") +
		ggplot2::scale_fill_manual(values = .group_cols, drop = FALSE) +
		ggplot2::labs(title = panel_title("C", paste0("Distribution of selected biomarkers", if (!is.null(Ylab)) paste0(": ", Ylab) else "")), x = NULL, y = "N biomarkers", fill = NULL) +
		fig_theme(17) +
		ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
}
