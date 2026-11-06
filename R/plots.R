
plot_results <- function(input, output) {
  input <- normalizePath(input, mustWork = TRUE)
  if (!dir.exists(output)) dir.create(output, recursive = TRUE)
  read_table <- function(name) read.csv(file.path(input, paste0(name, ".csv")))
  colors <- c("#0072B2", "#D55E00", "#009E73")
  render <- function(name, width, height, draw) {
    for (format in c("pdf", "png")) {
      path <- file.path(output, paste0(name, ".", format))
      if (format == "pdf") {
        grDevices::pdf(path, width = width, height = height,
                       family = "Helvetica", useDingbats = FALSE)
      } else {
        grDevices::png(path, width = width, height = height, units = "in", res = 300)
      }
      tryCatch({
        graphics::par(mar = c(4.5, 4.5, 2.5, 1), mgp = c(2.7, .8, 0),
                      las = 1, bty = "l", cex = 1)
        draw()
      }, finally = grDevices::dev.off())
    }
  }
  positions <- read_table("positions")
  selected <- c("central_sparse90", "high_sparse90")
  position_gains <- vapply(selected, function(rule) {
    d <- positions[positions$rule == rule, ]
    d$gain[order(d$label)]
  }, numeric(6))
  render("gain_by_MTD", 7.5, 4.8, function() {
    graphics::matplot(1:6, position_gains, type = "b", lty = 1, pch = c(16, 17),
      col = colors[1:2], ylim = range(c(0, position_gains)) * 1.15,
      xlab = "True MTD dose level", ylab = "PCS change (percentage points)",
      main = "Benefit by MTD location", xaxt = "n")
    graphics::axis(1, at = 1:6)
    graphics::abline(h = 0, col = "grey65", lty = 3)
    graphics::legend("topright", c("Central prior, 90% table", "High prior, 90% table"),
      col = colors[1:2], pch = c(16, 17), lty = 1, bty = "n", cex = .85)
  })
  summary <- read_table("summary")
  summary <- summary[summary$category == "bracketing", ]
  render("gain_and_changes", 8, 4.8, function() {
    graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.2, 2.5, 1))
    for (i in 1:2) {
      prior <- c("central", "high")[i]
      rules <- paste0(prior, c("_sparse80", "_sparse90", "_original"))
      d <- summary[match(paste(prior, rules), paste(summary$prior, summary$rule)), ]
      graphics::plot(d$override, d$PCS_gain, pch = 16, col = colors[i],
        xlim = c(0, max(d$override) * 1.2), ylim = c(0, max(d$PCS_gain) * 1.2),
        xlab = "Changed recommendations (%)", ylab = "PCS gain (percentage points)",
        main = paste(c("Central", "High")[i], "prior"))
      graphics::lines(d$override, d$PCS_gain, col = colors[i])
      graphics::text(d$override, d$PCS_gain, labels = c("80%", "90%", "Full"),
                     pos = c(2, 1, 3), cex = .8)
    }
  })
  render("prior_sensitivity", 7.5, 5, function() {
    rules <- c("central_sparse80", "central_sparse90", "high_sparse80", "high_sparse90")
    priors <- c("central", "high", "equal")
    gain <- vapply(priors, function(prior) {
      d <- summary[summary$prior == prior, ]
      d$PCS_gain[match(rules, d$rule)]
    }, numeric(length(rules)))
    graphics::par(mar = c(4.5, 8.5, 2.5, 1))
    graphics::plot(range(c(0, gain)) + c(-.05, .1), c(.5, 4.5), type = "n",
      yaxt = "n", ylab = "", xlab = "PCS change (percentage points)",
      main = "Performance under different location priors")
    graphics::axis(2, at = 4:1,
      labels = c("Central, 80% table", "Central, 90% table", "High, 80% table", "High, 90% table"),
      tick = FALSE)
    graphics::abline(v = 0, col = "grey65", lty = 3)
    for (j in seq_along(priors)) {
      graphics::points(gain[, j], 4:1 + (j - 2) * .14, pch = c(16, 17, 15)[j], col = colors[j])
    }
    graphics::legend("bottomright", c("Central prior", "High prior", "Equal prior"),
      col = colors, pch = c(16, 17, 15), bty = "n", cex = .85)
  })
  paired <- read_table("paired_changes")
  paired <- paired[paired$rule %in% selected, ]
  paired <- paired[match(selected, paired$rule), ]
  outcomes <- 10 * t(as.matrix(paired[c("corrected", "lost", "other")]))
  render("changed_recommendations", 7.5, 5, function() {
    graphics::barplot(outcomes, beside = TRUE, col = colors,
      names.arg = c("Central prior", "High prior"),
      ylim = c(0, max(outcomes) * 1.35),
      ylab = "Expected number per 1000 trials", main = "Consequences of the 90% tables")
    graphics::legend("topright", c("Corrected", "Correct selection lost", "Other change"),
      fill = colors, bty = "n", cex = .85)
  })
  clinical <- read_table("protocol_consequences")
  clinical <- clinical[clinical$prior %in% c("MN_0.5", "DD_0.5"), ]
  rules <- c("primary_MN", "expanded_MN", "uncertain_MN", "uncertain_sources")
  outcomes <- vapply(c("MN_0.5", "DD_0.5"), function(prior) {
    d <- clinical[clinical$prior == prior, ]
    -1000 * d$DLT_change[match(rules, d$rule)]
  }, numeric(length(rules)))
  render("expansion_consequences", 9, 5, function() {
    graphics::par(mar = c(4.5, 10, 2.5, 1))
    graphics::barplot(t(outcomes[4:1, ]), beside = TRUE, horiz = TRUE, col = colors[1:2],
      names.arg = rev(c("Fixed endpoints", "Flexible endpoints", "Uncertain preferences", "Uncertain sources")),
      xlim = c(min(0, outcomes) - .1, max(outcomes) * 1.4),
      xlab = "Expected reduction in DLTs per 1000 future trials",
      main = "Modeled consequences during dose expansion")
    graphics::abline(v = 0, col = "grey65", lty = 3)
    graphics::legend("bottomright", c("MN source, power 0.5", "DD source, power 0.5"),
      fill = colors[1:2], bty = "n", cex = .8)
  })
  files <- list.files(output, pattern = "\\.(pdf|png)$", full.names = TRUE)
  invisible(files)
}
