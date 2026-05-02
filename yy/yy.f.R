
# ☯
plot_yy <- function(dat, b2e, b2e.breaks, Xlab, Xs, Xs.fix = FALSE, Ylab.1, Ylab.2,
		Yadj, Xs.curve, curve.df = 3, Y.sd = 1, skip = c(2, 1), v_line = 0, pt = 1e-04, sd.bar = FALSE) {
	if (is.logical(Xs.curve)) Xs.curve <- NULL
	dat$b2e <- dat[[b2e]]
	yes_breaks <- !(isFALSE(b2e.breaks) || is.null(b2e.breaks))
	if (yes_breaks) dat$b2e[!is.na(dat$b2e) & (dat$b2e < min(b2e.breaks) | dat$b2e > max(b2e.breaks))] <- NA
	colrs <- grDevices::hcl(h = seq(15, 375, length.out = length(Xs) + 1)[1:length(Xs)], c = 95, l = 50)
	op <- par(mar = c(5, 5, 3.5, 5)); on.exit(par(op), add = TRUE)

	if (inherits(dat$b2e, "Date")) {
		if (yes_breaks) {
			myhist <- hist(dat$b2e, breaks = b2e.breaks, plot = FALSE)
			breaks <- b2e.breaks
		} else {
			myhist <- hist(dat$b2e, breaks = "years", plot = FALSE)
			breaks <- as.Date(myhist$breaks)
		}
		mids <- as.Date(myhist$mids); axes <- FALSE
	} else {
		if (yes_breaks) {
			myhist <- hist(dat$b2e, breaks = b2e.breaks, plot = FALSE)
			breaks <- b2e.breaks
		} else {
			myhist <- hist(dat$b2e, plot = FALSE)
			breaks <- myhist$breaks
		}
		mids <- myhist$mids; axes <- TRUE
	}

	bin.colrs <- ifelse(mids < v_line, adjustcolor("grey70", alpha.f = 0.5), adjustcolor("#8EC5FF", alpha.f = 0.5))
	hist(dat$b2e, breaks = breaks, border = "white", col = bin.colrs, main = "", axes = axes, xlab = Xlab, ylab = Ylab.1, font.lab = 2)
	grp <- cut(dat$b2e, breaks = breaks, include_lowest = TRUE, right = TRUE)
	avg_list <- sd_list <- yhat_list <- vector("list", length(Xs)); names(avg_list) <- names(sd_list) <- names(yhat_list) <- Xs
	p_vec <- b_vec <- rep(NA_real_, length(Xs))
	draw_fit <- rep(FALSE, length(Xs))
	fit_labs <- fit_cols <- character(0)

	for (i in seq_along(Xs)) {
		if (!(Xs[i] %in% names(dat)) || is.null(dat[[Xs[i]]])) next
		dat$Y <- dat[[Xs[i]]]
		if (all(is.na(dat$Y))) next
		if (Yadj) dat$Y <- residuals(lm(Y ~ sex + splines::ns(age, df = 3) + splines::ns(bmi, df = 2) + tdi + PC1 + PC2, data = dat, na.action = na.exclude))
		dat$Y <- as.vector(scale(dat$Y))
		avg <- tapply(dat$Y, grp, function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE))
		sdv <- tapply(dat$Y, grp, function(x) if (sum(is.finite(x)) < 2) NA_real_ else sd(x, na.rm = TRUE))
		avg[c(head(seq_along(avg), skip[1]), tail(seq_along(avg), skip[2]))] <- NA
		sdv[c(head(seq_along(sdv), skip[1]), tail(seq_along(sdv), skip[2]))] <- NA
		if (!is.null(Xs.fix) && as.character(i) %in% names(Xs.fix)) avg[as.numeric(names(Xs.fix[[as.character(i)]]))] <- unlist(Xs.fix[[as.character(i)]])

		ok <- is.finite(mids) & is.finite(avg)
		fit1 <- lm(Y ~ splines::ns(b2e, df = curve.df), data = dat)
		fit2 <- if (sum(ok) >= 4) stats::smooth.spline(x = mids[ok], y = avg[ok], df = curve.df) else NULL
		yhat2 <- if (!is.null(fit2)) predict(fit2, x = mids)$y else rep(NA_real_, length(mids))
		yhat2[c(head(seq_along(yhat2), skip[1]), tail(seq_along(yhat2), skip[2]))] <- NA
		summ1 <- summary(fit1); b1 <- summ1$coef[2, 1]; p1 <- summ1$coef[2, 4]
		has_spec <- is.numeric(Xs.curve) && length(Xs.curve) > 0
		yes_spec <- has_spec && i %in% Xs.curve
		yes_sig <- isTRUE(is.finite(p1) && p1 < pt)
		draw_fit[i] <- if (has_spec) yes_spec else yes_sig
		if (draw_fit[i]) {
			fit_labs <- c(fit_labs, sprintf("%s\n \u03B2 = %.3f\n p = %.2E", Xs[i], b1, p1))
			fit_cols <- c(fit_cols, colrs[i])
		}
		avg_list[[i]] <- avg; sd_list[[i]] <- sdv; yhat_list[[i]] <- yhat2; b_vec[i] <- b1; p_vec[i] <- p1
	}

	if (sd.bar) {
		ylo <- unlist(Map(function(m, s) m - Y.sd * s, avg_list, sd_list))
		yhi <- unlist(Map(function(m, s) m + Y.sd * s, avg_list, sd_list))
		yr <- range(c(ylo, yhi), na.rm = TRUE)
	} else {
		yr <- range(unlist(avg_list), na.rm = TRUE)
	}
	pad <- max(0.05, 0.12 * diff(yr))
	ylim2 <- range(pretty(yr + c(-pad, pad), n = 5))

	par(new = TRUE)
	plot(mids, rep(NA_real_, length(mids)), type = "n", xlim = range(myhist$breaks), ylim = ylim2, axes = FALSE, xlab = NA, ylab = NA)
	for (i in seq_along(Xs)) {
		avg <- avg_list[[i]]
		if (is.null(avg)) next
		lines(mids, as.numeric(avg), type = "b", col = colrs[i], pch = 16, cex = 2, lwd = 2, lty = 1)
		if (sd.bar) {
			sdv <- sd_list[[i]]
			arrows(mids, avg - Y.sd * sdv, mids, avg + Y.sd * sdv, angle = 90, code = 3, length = 0.04, col = adjustcolor(colrs[i], alpha.f = 0.6), lwd = 1.2)
		}
		if (draw_fit[i]) lines(mids, yhat_list[[i]], col = colrs[i], lwd = 2, lty = 3)
	}

	axis(4)
	mtext(Ylab.2, side = 4, line = 3, cex = 1.2, font = 2)
	abline(h = 0, col = adjustcolor("black", 0.7), lwd = 1.5, lty = 2)
	if (!is.na(v_line)) segments(x0 = v_line, y0 = ylim2[1], x1 = v_line, y1 = ylim2[2] * 0.85, col = adjustcolor("black", 0.7), lwd = 2, lty = 3)
	legend("topright", inset = c(0.002, 0.01), xpd = FALSE, legend = Xs, ncol = 1, col = colrs, lty = 2, pch = 10, cex = 0.86, text.col = colrs, text.font = 2,
		x.intersp = 0.05, y.intersp = 0.90, pt.cex = 0.75, seg.len = 0.35, box.lty = 1, box.lwd = 1.2, box.col = "grey80")	
	if (length(fit_labs) > 0)
		legend("left", bty = "n", inset = c(-0.06, 0), legend = paste0("\n", c(fit_labs, "\n")), ncol = 1, y.intersp = 1.02, cex = 0.95, text.col = fit_cols, xpd = NA)
	
	return(lapply(avg_list, round, 3))
}
