#install.packages("MASS")
library(MASS)

conc <- c(800, 400, 200, 100, 50)
n <- rep(20, 5)

pos_A <- c(20, 14, 7, 3, 3)
pos_B <- c(20, 17, 2, 1, 0)
pos_C <- c(20, 14, 8, 0, 0)
pos_D <- c(20, 17, 2, 1, 0)

calc_lod95 <- function(pos, n, conc) {
  model <- glm(cbind(pos, n - pos) ~ log2(conc),
               family = binomial(link = "probit"))
  b0 <- coef(model)[1]
  b1 <- coef(model)[2]

  lod95_log2 <- (qnorm(0.95) - b0) / b1
  lod95 <- 2^lod95_log2

  vcov_mat <- vcov(model)
  var_b0 <- vcov_mat[1, 1]
  var_b1 <- vcov_mat[2, 2]
  cov_b0b1 <- vcov_mat[1, 2]

  d_b0 <- -1 / b1
  d_b1 <- -(qnorm(0.95) - b0) / b1^2
  var_lod95 <- d_b0^2 * var_b0 + d_b1^2 * var_b1 + 2 * d_b0 * d_b1 * cov_b0b1
  se_lod95 <- sqrt(var_lod95)

  ci_lower <- 2^(lod95_log2 - 1.96 * se_lod95)
  ci_upper <- 2^(lod95_log2 + 1.96 * se_lod95)

  return(list(model = model, lod95 = lod95, lod95_log2 = lod95_log2,
              ci_lower = ci_lower, ci_upper = ci_upper))
}

plot_probit <- function(pos, n, conc, name) {

  rate <- pos / n
  res <- calc_lod95(pos, n, conc)
  model <- res$model
  lod95 <- res$lod95
  lod95_log2 <- res$lod95_log2
  ci_lower <- res$ci_lower
  ci_upper <- res$ci_upper

  x_min <- log2(50) - 0.5
  x_max <- max(log2(800), log2(ci_upper)) + 1.5

  x_seq_conc <- 2^seq(x_min, x_max, length.out = 200)

  pred <- predict(model, newdata = data.frame(conc = x_seq_conc),
                  type = "link", se.fit = TRUE)
  fit_prob <- pnorm(pred$fit)
  lower_prob <- pnorm(pred$fit - 1.96 * pred$se.fit)
  upper_prob <- pnorm(pred$fit + 1.96 * pred$se.fit)

  plot(log2(conc), rate,
       xlab = "Concentration (copies/reaction)",
       ylab = "Detection rate",
       xaxt = "n", xlim = c(x_min, x_max), ylim = c(0, 1),
       pch = 19, cex = 1.2, cex.lab = 1.1,
       main = paste("Genotype", name))

  axis(1, at = log2(c(50, 100, 200, 400, 800)),
       labels = c(50, 100, 200, 400, 800))

  polygon(c(log2(x_seq_conc), rev(log2(x_seq_conc))),
          c(lower_prob, rev(upper_prob)),
          col = rgb(0, 0, 1, 0.15), border = NA)

  lines(log2(x_seq_conc), fit_prob, col = "blue", lwd = 2)

  abline(h = 0.95, lty = 2, col = "red")

  abline(v = lod95_log2, lty = 3, col = "red")

  text(x = x_max - 0.2, y = 0.85,
       labels = paste0("LoD95 = ", round(lod95, 0),
                       "\n(95% CI: ", round(ci_lower, 0), "-", round(ci_upper, 0), ")"),
       adj = 1, cex = 0.8, col = "red")
}

pdf("qPCR.20rep.lod95.pdf")
par(mfrow = c(2, 2))
plot_probit(pos_A, n, conc, "A")
plot_probit(pos_B, n, conc, "B")
plot_probit(pos_C, n, conc, "C")
plot_probit(pos_D, n, conc, "D")
par(mfrow = c(1, 1))
dev.off()
