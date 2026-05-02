library(data.table)
library(ape)
library(eoffice)
args <- commandArgs(TRUE)
root <- normalizePath(if (length(args)) args[1] else Sys.getenv("BALD_RES", "D:/bald/res"), winslash = "/", mustWork = FALSE)
phy0 <- file.path(root, "phy"); out0 <- file.path(root, "plot")
dir.create(out0, recursive = TRUE, showWarnings = FALSE)

ppt_full <- file.path(out0, "s8_tree_full.pptx")
ppt_main <- file.path(out0, "s8_tree_main.pptx")
if (file.exists(ppt_full)) file.remove(ppt_full)
if (file.exists(ppt_main)) file.remove(ppt_main)

min_boot <- 70
scale_bar_len <- 0.005
arch_col <- "#1f78b4"
anc_col <- "#33a02c"
hi_col <- "#d73027"

style <- function(tag) {
  if (tag == "main") {
    list(dot_r = 0.0035, cex_tip = 0.36, cex_boot = 0.36,
         lwd_tree = 0.45, lwd_guide = 0.22,
         r_lab = 1.10, r_text = 1.18, r_lim = 1.30,
         ppt_w = 3.6, ppt_h = 3.6)
  } else {
    list(dot_r = 0.0028, cex_tip = 0.27, cex_boot = 0.28,
         lwd_tree = 0.34, lwd_guide = 0.16,
         r_lab = 1.08, r_text = 1.15, r_lim = 1.32,
         ppt_w = 4.0, ppt_h = 3.8)
  }
}

topptx_safe <- function(filename, width, height, append = FALSE) {
  if ("append" %in% names(formals(topptx))) {
    topptx(filename = filename, width = width, height = height, append = append)
  } else {
    if (append) warning("topptx() has no append argument; output may be overwritten")
    topptx(filename = filename, width = width, height = height)
  }
}

read_phyml_tree <- function(f) {
  x <- readLines(f, warn = FALSE)
  x <- x[nzchar(trimws(x))]
  i <- grep("^\\(", x)
  if (!length(i)) stop("no tree found: ", f)
  read.tree(text = x[i[1]])
}

desc_tips <- function(tr, node) {
  nt <- Ntip(tr)
  kid <- tr$edge[tr$edge[, 1] == node, 2]
  unlist(lapply(kid, function(k) if (k <= nt) k else desc_tips(tr, k)), use.names = FALSE)
}

parent_node <- function(tr, node) {
  x <- tr$edge[tr$edge[, 2] == node, 1]
  if (length(x)) x[1] else NA_integer_
}

sector_range <- function(a) {
  a <- sort((a + 2 * pi) %% (2 * pi))
  if (length(a) <= 1) return(c(a, a))
  g <- c(diff(a), a[1] + 2 * pi - a[length(a)])
  j <- which.max(g)
  s <- a[(j %% length(a)) + 1]
  e <- a[j]
  if (e < s) e <- e + 2 * pi
  c(s, e)
}

tip_type <- function(meta, labs) {
  out <- rep(NA_character_, length(labs))
  names(out) <- labs
  if (!is.null(meta) && all(c("label", "type") %in% names(meta))) {
    idx <- match(labs, meta$label)
    out[!is.na(idx)] <- as.character(meta$type[idx[!is.na(idx)]])
  }
  out[labs %in% c("Vindija", "Altai", "Chagyr", "Denisova")] <- "archaic"
  out[labs == "Ancestral"] <- "ancestral"
  out[is.na(out)] <- "modern"
  out
}

highlight_node <- function(tr, meta) {
  nt <- Ntip(tr)
  labs <- tr$tip.label
  tp <- tip_type(meta, labs)
  arch <- labs[tp == "archaic"]
  modern <- labs[tp == "modern"]
  if (!length(arch) || !length(modern)) return(NA_integer_)

  candidates <- function(all_arch = TRUE) {
    rbindlist(lapply((nt + 1):(nt + tr$Nnode), function(node) {
      tips <- labs[desc_tips(tr, node)]
      ok_arch <- if (all_arch) all(arch %in% tips) else any(tips %in% arch)
      support <- suppressWarnings(as.numeric(tr$node.label[node - nt]))
      if (!ok_arch || !any(tips %in% modern) || any(tips == "Ancestral") ||
          is.na(support) || support < min_boot) return(NULL)
      data.table(node = node, support = support, n_tips = length(tips),
                 n_mod = sum(tips %in% modern))
    }), fill = TRUE)
  }

  z <- candidates(TRUE)
  if (!nrow(z)) z <- candidates(FALSE)
  if (!nrow(z)) return(NA_integer_)
  setorder(z, n_tips, -support, -n_mod)
  z$node[1]
}

draw_tree <- function(tr, meta = NULL, tag = c("full", "main")) {
  tag <- match.arg(tag)
  st <- style(tag)

  tr$tip.label <- trimws(tr$tip.label)
  if ("Ancestral" %in% tr$tip.label) {
    tr <- tryCatch(root(tr, outgroup = "Ancestral", resolve.root = TRUE), error = function(e) tr)
  }
  tr <- ladderize(tr)

  nt <- Ntip(tr)
  ia <- which(tr$tip.label == "Ancestral")
  dep <- node.depth.edgelength(tr)[1:nt]
  r_tree <- max(dep[setdiff(seq_len(nt), ia)], na.rm = TRUE)
  if (!is.finite(r_tree)) r_tree <- max(dep, na.rm = TRUE)

  r_lab <- r_tree * st$r_lab
  r_text <- r_tree * st$r_text
  r_lim <- r_tree * st$r_lim

  par(mar = c(0.2, 0.2, 0.2, 0.2), xpd = NA, pty = "s")
  plot.phylo(tr, type = "fan", use.edge.length = TRUE, show.tip.label = FALSE,
             no.margin = TRUE, edge.width = st$lwd_tree,
             x.lim = c(-r_lim, r_lim), y.lim = c(-r_lim, r_lim))

  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  xx <- pp$xx
  yy <- pp$yy
  xt <- xx[1:nt]
  yt <- yy[1:nt]
  ang <- atan2(yt, xt)
  labs <- tr$tip.label
  tp <- tip_type(meta, labs)

  hi <- highlight_node(tr, meta)
  if (!is.na(hi)) {
    idx <- desc_tips(tr, hi)
    aa <- sector_range(atan2(yy[idx], xx[idx]))
    th <- seq(aa[1], aa[2], length.out = 500)
    r0 <- sqrt(xx[hi]^2 + yy[hi]^2)
    polygon(c(r0 * cos(th), rev(r_lab * cos(th))),
            c(r0 * sin(th), rev(r_lab * sin(th))),
            col = adjustcolor(hi_col, alpha.f = 0.22), border = NA)
  }

  segments(xt, yt, r_lab * cos(ang), r_lab * sin(ang),
           lty = 3, col = "grey45", lwd = st$lwd_guide)
  symbols(xt, yt, circles = rep(r_tree * st$dot_r, nt), inches = FALSE,
          add = TRUE, bg = "black", fg = "black", lwd = 0.15)

  lab_col <- rep("black", nt)
  lab_col[tp == "archaic"] <- arch_col
  lab_col[tp == "ancestral"] <- anc_col

  deg <- ang * 180 / pi
  flip <- deg < -90 | deg > 90
  srt <- ifelse(flip, deg + 180, deg)
  for (i in seq_len(nt)) {
    text(r_text * cos(ang[i]), r_text * sin(ang[i]), labs[i],
         srt = srt[i], adj = if (flip[i]) c(1, 0.5) else c(0, 0.5),
         cex = st$cex_tip, col = lab_col[i])
  }

  if (!is.na(hi)) {
    boot <- suppressWarnings(as.numeric(tr$node.label))
    nd <- nt + seq_len(tr$Nnode)
    show <- unique(c(parent_node(tr, hi), hi))
    idx <- match(show[show > nt], nd)
    idx <- idx[!is.na(idx) & !is.na(boot[idx])]
    if (length(idx)) nodelabels(boot[idx], node = nd[idx], frame = "n",
                                cex = st$cex_boot, col = "black")
  }

  usr <- par("usr")
  add.scale.bar(x = usr[2] - 0.14 * diff(usr[1:2]),
                y = usr[3] + 0.06 * diff(usr[3:4]),
                length = scale_bar_len, lwd = 0.8, cex = 0.65)
}

files <- list.files(phy0, pattern = "\\.phy_phyml_tree\\.txt$", recursive = TRUE, full.names = TRUE)
if (!length(files)) quit(save = "no", status = 0)
n_ppt <- c(full = 0L, main = 0L)

for (treef in files) {
  base <- sub("\\.phy_phyml_tree\\.txt$", "", basename(treef))
  tag <- if (grepl("\\.main$", base)) "main" else "full"
  st <- style(tag)
  metaf <- file.path(dirname(treef), paste0(base, ".meta.tsv"))
  meta <- if (file.exists(metaf)) fread(metaf) else NULL
  tr <- tryCatch(read_phyml_tree(treef), error = function(e) NULL)
  if (is.null(tr)) next

  tf <- tempfile(fileext = ".png")
  png(tf, width = st$ppt_w, height = st$ppt_h, units = "in", res = 300)
  draw_tree(tr, meta, tag)
  topptx_safe(if (tag == "main") ppt_main else ppt_full,
              width = st$ppt_w, height = st$ppt_h,
              append = n_ppt[tag] > 0L)
  dev.off(); unlink(tf)

  n_ppt[tag] <- n_ppt[tag] + 1L
}
