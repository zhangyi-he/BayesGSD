#' @title Shifting from overly aggressive to more conservative early stopping in Bayesian group sequential designs
#' @author Zhangyi He, Feng Yu, Suzie Cro, Laurent Billot

#' This R script implements the operating characteristics evaluation for the manuscript.

setwd("~/Dropbox/Jeffery He/iResearch/Publications/2024/ZH2026-ClinTrial-StatAdv-BayesGSD-AlphaSpending/Code/Code v1.0")

if (!require("jsonlite")) {
  install.packages("jsonlite")
  library("jsonlite")
}
# if (!require("CalvinBayes")) {
#   devtools::install_github("rpruim/CalvinBayes")
#   library("CalvinBayes")
# }
if (!require("corrplot")) {
  install.packages("corrplot")
  library("corrplot")
}
if (!require("ggplot2")) {
  install.packages("ggplot2")
  library("ggplot2")
}
if (!require("ggsci")) {
  install.packages("ggsci")
  library("ggsci")
}
if (!require("plot3D")) {
  install.packages("plot3D")
  library("plot3D")
}
if (!require("rpact")) {
  install.packages("rpact")
  library("rpact")
}

#' label panels
fig_label <- function(text, region = "figure", pos = "topleft", cex = NULL, ...) {

  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", "left", "center", "right", "bottomleft", "bottom", "bottomright"))

  if (region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from = "in", to = "user")
    y <- grconvertY(c(0, ds[2]), from = "in", to = "user")

    # fragment of the device we use to plot
    if (region == "figure") {
      # account for the fragment of the device that
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    }
  }

  # much simpler if in plotting region
  if (region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }

  sw <- strwidth(text, cex = cex) * 60/100
  sh <- strheight(text, cex = cex) * 60/100

  x1 <- switch(pos,
               topleft     = x[1] + sw,
               left        = x[1] + sw,
               bottomleft  = x[1] + sw,
               top         = (x[1] + x[2])/2,
               center      = (x[1] + x[2])/2,
               bottom      = (x[1] + x[2])/2,
               topright    = x[2] - sw,
               right       = x[2] - sw,
               bottomright = x[2] - sw)

  y1 <- switch(pos,
               topleft     = y[2] - sh,
               top         = y[2] - sh,
               topright    = y[2] - sh,
               left        = (y[1] + y[2])/2,
               center      = (y[1] + y[2])/2,
               right       = (y[1] + y[2])/2,
               bottomleft  = y[1] + sh,
               bottom      = y[1] + sh,
               bottomright = y[1] + sh)

  old.par <- par(xpd = NA)
  on.exit(par(old.par))

  text(x1, y1, text, cex = cex, ...)

  return(invisible(c(x,y)))
}

################################################################################

# Figure 1
jpeg(file = "./fig11.jpeg", width = 12, height = 9, units = "in", res = 500)
par(mar = c(5.5, 5, 5.5, 5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
HP <- getDesignGroupSequential(kMax = 5,
                               alpha = 0.025,
                               beta = 0.2,
                               sided = 1,
                               typeOfDesign = "HP")
P <- getDesignGroupSequential(kMax = 5,
                              alpha = 0.025,
                              beta = 0.2,
                              sided = 1,
                              typeOfDesign = "P")
OF <- getDesignGroupSequential(kMax = 5,
                               alpha = 0.025,
                               beta = 0.2,
                               sided = 1,
                               typeOfDesign = "OF")
asP <- getDesignGroupSequential(kMax = 5,
                                alpha = 0.025,
                                beta = 0.2,
                                sided = 1,
                                typeOfDesign = "asP")
asOF <- getDesignGroupSequential(kMax = 5,
                                 alpha = 0.025,
                                 beta = 0.2,
                                 sided = 1,
                                 typeOfDesign = "asOF")

plot(x = (1:5) / 5, y = HP$criticalValues,
     type = "b", pch = 19, lwd = 3, col = pal_nejm("default", alpha = 1)(8)[1],
     bty = "n", xaxt = "n", xlim = c(1/5, 1), ylim = c(1.5, 5),
     main = "Stopping boundaries", xlab = "Information rate", ylab = "Critical value")
axis(1, at = (1:5) / 5, labels = c("1/5", "2/5", "3/5", "4/5", "1"))
lines(x = (1:5) / 5, y = P$criticalValues,
      type = "b", pch = 19, lwd = 3, col = pal_nejm("default", alpha = 1)(8)[2])
lines(x = (1:5) / 5, y = OF$criticalValues,
      type = "b", pch = 19, lwd = 3, col = pal_nejm("default", alpha = 1)(8)[3])
lines(x = (1:5) / 5, y = asP$criticalValues,
      type = "b", pch = 19, lwd = 3, col = pal_nejm("default", alpha = 1)(8)[4])
lines(x = (1:5) / 5, y = asOF$criticalValues,
      type = "b", pch = 19, lwd = 3, col = pal_nejm("default", alpha = 1)(8)[5])
abline(h = qnorm(1 - 0.025), lty = 2, lwd = 2)
legend("topright", legend = c("Haybittle-Peto", "Pocock", "O'Brien-Fleming", "Pocock-type spending", "O'Brien-Fleming-type spending"),
       col = pal_nejm("default", alpha = 1)(8)[1:5], lty = 1, lwd = 3, bty = "n", cex = 1.5)
dev.off()

############################################################

# conventional Bayesian GSD with 2 interim analyses
sampleSize <- (163:187) * 2
treatmentEffect <- c(0.40, 0.25)
allocationRatio <- c(1, 1)
infoFraction <- (1:3) / 3
probThreshold <- seq(from = 0.980, to = 0.999, by = 0.0005)
type1Error <- matrix(data = NA, nrow = length(sampleSize), ncol = length(probThreshold))
type2Error <- matrix(data = NA, nrow = length(sampleSize), ncol = length(probThreshold))
stopProbH0 <- array(data = NA, dim = c(length(sampleSize), length(probThreshold), length(infoFraction)))
stopProbH1 <- array(data = NA, dim = c(length(sampleSize), length(probThreshold), length(infoFraction)))
for (i in 1:length(sampleSize)) {
  for (j in 1:length(probThreshold)) {
    config <- list(
      target = "op_char",
      n_threads = 8,
      effect_grid_size = 1001,
      arm_group_type = 0,
      prior = list(
        list(
          list(1, 1, 1)
        ),
        list(
          list(1, 1, 1)
        )
      ),
      allocation_ratio = as.list(allocationRatio),
      information_fracs = as.list(infoFraction),
      total_final_sample_size = sampleSize[i],
      theta = as.list(treatmentEffect),
      bin_size = list(1, 1, 1),
      prob_threshs = list(1e-08, 1e-08, 1e-08),
      criteria = list(
        list(0, "diff", 0, 1, list(-1.01, 0.00), ">",
             list(probThreshold[j], probThreshold[j], probThreshold[j]))
      ),
      decisions = list(
        list(1, "eff1", "post", list(0))
      ),
      print_info = FALSE
    )
    write_json(config, "config_custom.json", auto_unbox = TRUE, pretty = TRUE)

    system2(command = "./main",
            args = c("-c", "config_custom.json", "-o", "output_custom.json"))

    result <- fromJSON("output_custom.json")
    type1Error[i, j] <- result$p_eff_any_stage0
    type2Error[i, j] <- 1 - result$p_eff_any_stage1
    stopProbH0[i, j, ] <- unlist(lapply(X = result$p_efc0, FUN = function(x) {x[1, 2]}))
    stopProbH1[i, j, ] <- unlist(lapply(X = result$p_efc1, FUN = function(x) {x[1, 2]}))
  }
}
file.remove("config_custom.json", "output_custom.json")

save(sampleSize, treatmentEffect, allocationRatio, infoFraction, probThreshold,
     type1Error, type2Error, stopProbH0, stopProbH1,
     file = "~/Dropbox/Jeffery He/iResearch/Publications/2024/ZH2026-ClinTrial-StatAdv-BayesGSD-AlphaSpending/Output/Output v1.0/HYPRESS_3Stage_convBayesGSD.RData")

####################

# Figure A1
jpeg(file = "./fig21.jpeg", width = 12, height = 9, units = "in", res = 500)
par(mar = c(3.5, 7.5, 5.5, 5.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
load("~/Dropbox/Jeffery He/iResearch/Publications/2024/ZH2026-ClinTrial-StatAdv-BayesGSD-AlphaSpending/Output/Output v1.0/HYPRESS_3Stage_convBayesGSD.RData")
type1Error <- type1Error * 100
rownames(type1Error) <- sampleSize
colnames(type1Error) <- formatC(x = probThreshold, digits = 4, format = "f")

zlim <- range(type1Error, na.rm = TRUE)
n_col <- 200
prop <- min(max((0.025 * 100 - zlim[1]) / diff(zlim), 0), 1)
n1 <- max(1, round(n_col * prop))
n2 <- max(1, n_col - n1)
cols <- c(colorRampPalette(c(pal_nejm("default", alpha = 1)(8)[2], "white"))(n1),
          colorRampPalette(c("white", pal_nejm("default", alpha = 1)(8)[1]))(n2))

corrplot(corr = type1Error,
         method = "square",
         type = "full",
         col = cols,
         is.corr = FALSE,
         title = "",
         diag = TRUE,
         mar = c(1, 1, 2, 0),
         tl.col  = "black",
         tl.srt  = 90,
         cl.pos  = "r",
         cl.lim  = zlim)
title("Type I Error Rate", line = 1.0, cex.main = 2.0)
mtext("Probability threshold", side = 1, line = -0.5, cex = 1.1)
mtext("Total sample size", side = 2, line = 6.5, cex = 1.1)
dev.off()

####################

# Figure A2
jpeg(file = "./fig22.jpeg", width = 12, height = 9, units = "in", res = 500)
par(mar = c(3.5, 7.5, 5.5, 5.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
load("~/Dropbox/Jeffery He/iResearch/Publications/2024/ZH2026-ClinTrial-StatAdv-BayesGSD-AlphaSpending/Output/Output v1.0/HYPRESS_3Stage_convBayesGSD.RData")
power <- (1 - type2Error) * 100
rownames(power) <- sampleSize
colnames(power) <- formatC(x = probThreshold, digits = 4, format = "f")

zlim <- range(power, na.rm = TRUE)
n_col <- 200
prop <- min(max((0.800 * 100 - zlim[1]) / diff(zlim), 0), 1)
n1 <- max(1, round(n_col * prop))
n2 <- max(1, n_col - n1)
cols <- c(colorRampPalette(c(pal_nejm("default", alpha = 1)(8)[2], "white"))(n1),
          colorRampPalette(c("white", pal_nejm("default", alpha = 1)(8)[1]))(n2))

corrplot(corr = power,
         method = "square",
         type = "full",
         col = cols,
         is.corr = FALSE,
         title = "",
         diag = TRUE,
         mar = c(1, 1, 2, 0),
         tl.col  = "black",
         tl.srt  = 90,
         cl.pos  = "r",
         cl.lim  = zlim)
title("Power", line = 1.0, cex.main = 2.0)
mtext("Probability threshold", side = 1, line = -0.5, cex = 1.1)
mtext("Total sample size", side = 2, line = 6.5, cex = 1.1)
dev.off()

####################

# Figure 2
jpeg(file = "./fig23.jpeg", width = 12, height = 9, units = "in", res = 500)
par(mar = c(5.5, 5, 5.5, 5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
HP <- getDesignGroupSequential(kMax = 3,
                               alpha = 0.025,
                               beta = 0.2,
                               sided = 1,
                               typeOfDesign = "HP")
getSampleSizeRates(design = HP,
                   groups = 2,
                   pi1 = 0.25,
                   pi2 = 0.40,
                   allocationRatioPlanned = 1)
asP <- getDesignGroupSequential(kMax = 3,
                                alpha = 0.025,
                                beta = 0.2,
                                sided = 1,
                                typeOfDesign = "asP")
getSampleSizeRates(design = asP,
                   groups = 2,
                   pi1 = 0.25,
                   pi2 = 0.40,
                   allocationRatioPlanned = 1)
asOF <- getDesignGroupSequential(kMax = 3,
                                 alpha = 0.025,
                                 beta = 0.2,
                                 sided = 1,
                                 typeOfDesign = "asOF")
getSampleSizeRates(design = asOF,
                   groups = 2,
                   pi1 = 0.25,
                   pi2 = 0.40,
                   allocationRatioPlanned = 1)
load("~/Dropbox/Jeffery He/iResearch/Publications/2024/ZH2026-ClinTrial-StatAdv-BayesGSD-AlphaSpending/Output/Output v1.0/HYPRESS_3Stage_convBayesGSD.RData")
for (i in 1:length(sampleSize)) {
  if (any(type2Error[i, which(type1Error[i, ] <= 0.025)] <= 0.200)) {
    index <- c(i, min(intersect(which(type1Error[i, ] <= 0.025), which(type2Error[i, ] <= 0.200))))
    break
  }
}
print(sampleSize[index[1]])
print(probThreshold[index[2]])
print(type1Error[index[1], index[2]])
print(1 - type2Error[index[1], index[2]])
cumSampleSize <- ceiling(infoFraction * sampleSize[index[1]] / 2) * 2
print(sum(stopProbH0[index[1], index[2], 1:2] * cumSampleSize[1:2]) +
        (1 - sum(stopProbH0[index[1], index[2], 1:2])) * cumSampleSize[3])
print(sum(stopProbH1[index[1], index[2], 1:2] * cumSampleSize[1:2]) +
        (1 - sum(stopProbH1[index[1], index[2], 1:2])) * cumSampleSize[3])

plot(x = (0:3) / 3, y = cumsum(c(0, stopProbH0[index[1], index[2], ])),
     type = "b", pch = 19, lwd = 3, col = pal_nejm("default", alpha = 1)(8)[1],
     bty = "n", xaxt = "n", xlim = c(0, 1), ylim = c(0, 0.025),
     main = "Type I error spending", xlab = "Information rate", ylab = "Cumulative alpha spending")
axis(1, at = (0:3) / 3, labels = c("0", "1/3", "2/3", "1"))
lines(x = (0:3) / 3, y = c(0, HP$alphaSpent),
      type = "b", pch = 19, lwd = 3, col = pal_nejm("default", alpha = 1)(8)[2])
lines(x = (0:3) / 3, y = c(0, asP$alphaSpent),
      type = "b", pch = 19, lwd = 3, col = pal_nejm("default", alpha = 1)(8)[3])
lines(x = (0:3) / 3, y = c(0, asOF$alphaSpent),
      type = "b", pch = 19, lwd = 3, col = pal_nejm("default", alpha = 1)(8)[4])
legend("topleft", legend = c("BayesGSD (conventional)", "Haybittle-Peto", "Pocock-type spending", "O'Brien-Fleming-type spending"),
       col = pal_nejm("default", alpha = 1)(8)[1:4], lty = 1, lwd = 3, bty = "n", cex = 1.5)
dev.off()

############################################################

# conventional Bayesian GSD with 4 interim analyses
sampleSize <- (178:202) * 2
treatmentEffect <- c(0.40, 0.25)
allocationRatio <- c(1, 1)
infoFraction <- (1:5) / 5
probThreshold <- seq(from = 0.980, to = 0.999, by = 0.0005)
type1Error <- matrix(data = NA, nrow = length(sampleSize), ncol = length(probThreshold))
type2Error <- matrix(data = NA, nrow = length(sampleSize), ncol = length(probThreshold))
stopProbH0 <- array(data = NA, dim = c(length(sampleSize), length(probThreshold), length(infoFraction)))
stopProbH1 <- array(data = NA, dim = c(length(sampleSize), length(probThreshold), length(infoFraction)))
for (i in 1:length(sampleSize)) {
  for (j in 1:length(probThreshold)) {
    config <- list(
      target = "op_char",
      n_threads = 8,
      effect_grid_size = 1001,
      arm_group_type = 0,
      prior = list(
        list(
          list(1, 1, 1)
        ),
        list(
          list(1, 1, 1)
        )
      ),
      allocation_ratio = as.list(allocationRatio),
      information_fracs = as.list(infoFraction),
      total_final_sample_size = sampleSize[i],
      theta = as.list(treatmentEffect),
      bin_size = list(1, 1, 1, 1, 1),
      prob_threshs = list(1e-08, 1e-08, 1e-08, 1e-08, 1e-08),
      criteria = list(
        list(0, "diff", 0, 1, list(-1.01, 0.00), ">",
             list(probThreshold[j], probThreshold[j],
                  probThreshold[j], probThreshold[j], probThreshold[j]))
      ),
      decisions = list(
        list(1, "eff1", "post", list(0))
      ),
      print_info = FALSE
    )
    write_json(config, "config_custom.json", auto_unbox = TRUE, pretty = TRUE)

    system2(command = "./main",
            args = c("-c", "config_custom.json", "-o", "output_custom.json"))

    result <- fromJSON("output_custom.json")
    type1Error[i, j] <- result$p_eff_any_stage0
    type2Error[i, j] <- 1 - result$p_eff_any_stage1
    stopProbH0[i, j, ] <- unlist(lapply(X = result$p_efc0, FUN = function(x) {x[1, 2]}))
    stopProbH1[i, j, ] <- unlist(lapply(X = result$p_efc1, FUN = function(x) {x[1, 2]}))
  }
}
file.remove("config_custom.json", "output_custom.json")

save(sampleSize, treatmentEffect, allocationRatio, infoFraction, probThreshold,
     type1Error, type2Error, stopProbH0, stopProbH1,
     file = "~/Dropbox/Jeffery He/iResearch/Publications/2024/ZH2026-ClinTrial-StatAdv-BayesGSD-AlphaSpending/Output/Output v1.0/HYPRESS_5Stage_convBayesGSD.RData")

####################

# Figure A3
jpeg(file = "./fig301.jpeg", width = 12, height = 9, units = "in", res = 500)
par(mar = c(3.5, 7.5, 5.5, 5.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
load("~/Dropbox/Jeffery He/iResearch/Publications/2024/ZH2026-ClinTrial-StatAdv-BayesGSD-AlphaSpending/Output/Output v1.0/HYPRESS_5Stage_convBayesGSD.RData")
type1Error <- type1Error * 100
rownames(type1Error) <- sampleSize
colnames(type1Error) <- formatC(x = probThreshold, digits = 4, format = "f")

zlim <- range(type1Error, na.rm = TRUE)
n_col <- 200
prop <- min(max((0.025 * 100 - zlim[1]) / diff(zlim), 0), 1)
n1 <- max(1, round(n_col * prop))
n2 <- max(1, n_col - n1)
cols <- c(colorRampPalette(c(pal_nejm("default", alpha = 1)(8)[2], "white"))(n1),
          colorRampPalette(c("white", pal_nejm("default", alpha = 1)(8)[1]))(n2))

corrplot(corr = type1Error,
         method = "square",
         type = "full",
         col = cols,
         is.corr = FALSE,
         title = "",
         diag = TRUE,
         mar = c(1, 1, 2, 0),
         tl.col  = "black",
         tl.srt  = 90,
         cl.pos  = "r",
         cl.lim  = zlim)
title("Type I Error Rate", line = 1.0, cex.main = 2.0)
mtext("Probability threshold", side = 1, line = -0.5, cex = 1.1)
mtext("Total sample size", side = 2, line = 6.5, cex = 1.1)
dev.off()

####################

# Figure A4
jpeg(file = "./fig302.jpeg", width = 12, height = 9, units = "in", res = 500)
par(mar = c(3.5, 7.5, 5.5, 5.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
load("~/Dropbox/Jeffery He/iResearch/Publications/2024/ZH2026-ClinTrial-StatAdv-BayesGSD-AlphaSpending/Output/Output v1.0/HYPRESS_5Stage_convBayesGSD.RData")
power <- (1 - type2Error) * 100
rownames(power) <- sampleSize
colnames(power) <- formatC(x = probThreshold, digits = 4, format = "f")

zlim <- range(power, na.rm = TRUE)
n_col <- 200
prop <- min(max((0.800 * 100 - zlim[1]) / diff(zlim), 0), 1)
n1 <- max(1, round(n_col * prop))
n2 <- max(1, n_col - n1)
cols <- c(colorRampPalette(c(pal_nejm("default", alpha = 1)(8)[2], "white"))(n1),
          colorRampPalette(c("white", pal_nejm("default", alpha = 1)(8)[1]))(n2))

corrplot(corr = power,
         method = "square",
         type = "full",
         col = cols,
         is.corr = FALSE,
         title = "",
         diag = TRUE,
         mar = c(1, 1, 2, 0),
         tl.col  = "black",
         tl.srt  = 90,
         cl.pos  = "r",
         cl.lim  = zlim)
title("Power", line = 1.0, cex.main = 2.0)
mtext("Probability threshold", side = 1, line = -0.5, cex = 1.1)
mtext("Total sample size", side = 2, line = 6.5, cex = 1.1)
dev.off()

####################

#
jpeg(file = "./fig303.jpeg", width = 12, height = 9, units = "in", res = 500)
par(mar = c(5.5, 5, 5.5, 5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
HP <- getDesignGroupSequential(kMax = 5,
                               alpha = 0.025,
                               beta = 0.2,
                               sided = 1,
                               typeOfDesign = "HP")
asP <- getDesignGroupSequential(kMax = 5,
                                alpha = 0.025,
                                beta = 0.2,
                                sided = 1,
                                typeOfDesign = "asP")
asOF <- getDesignGroupSequential(kMax = 5,
                                 alpha = 0.025,
                                 beta = 0.2,
                                 sided = 1,
                                 typeOfDesign = "asOF")
load("~/Dropbox/Jeffery He/iResearch/Publications/2024/ZH2026-ClinTrial-StatAdv-BayesGSD-AlphaSpending/Output/Output v1.0/HYPRESS_5Stage_convBayesGSD.RData")
for (i in 1:length(sampleSize)) {
  if (any(type2Error[i, which(type1Error[i, ] <= 0.025)] <= 0.200)) {
    index <- c(i, min(intersect(which(type1Error[i, ] <= 0.025), which(type2Error[i, ] <= 0.200))))
    break
  }
}
print(sampleSize[index[1]])
print(probThreshold[index[2]])

plot(x = (0:5) / 5, y = cumsum(c(0, stopProbH0[index[1], index[2], ])),
     type = "b", pch = 19, lwd = 3, col = pal_nejm("default", alpha = 1)(8)[1],
     bty = "n", xaxt = "n", xlim = c(0, 1), ylim = c(0, 0.025),
     main = "Type I error spending", xlab = "Information rate", ylab = "Cumulative alpha spending")
axis(1, at = (0:5) / 5, labels = c("0", "1/5", "2/5", "3/5", "4/5", "1"))
lines(x = (0:5) / 5, y = c(0, HP$alphaSpent),
      type = "b", pch = 19, lwd = 3, col = pal_nejm("default", alpha = 1)(8)[2])
lines(x = (0:5) / 5, y = c(0, asP$alphaSpent),
      type = "b", pch = 19, lwd = 3, col = pal_nejm("default", alpha = 1)(8)[3])
lines(x = (0:5) / 5, y = c(0, asOF$alphaSpent),
      type = "b", pch = 19, lwd = 3, col = pal_nejm("default", alpha = 1)(8)[4])
legend("topleft", legend = c("BayesGSD (conventional)", "Haybittle-Peto", "Pocock-type spending", "O'Brien-Fleming-type spending"),
       col = pal_nejm("default", alpha = 1)(8)[1:4], lty = 1, lwd = 3, bty = "n", cex = 1.5)
dev.off()

############################################################

# conventional Bayesian GSD with 4 interim analyses (refinement 1: early vs late phase)
sampleSize <- 368
treatmentEffect <- c(0.40, 0.25)
allocationRatio <- c(1, 1)
infoFraction <- (1:5) / 5
earlyThreshold <- seq(from = 0.980, to = 0.999, by = 0.0005)
lateThreshold <- seq(from = 0.980, to = 0.999, by = 0.0005)
type1Error <- matrix(data = NA, nrow = length(earlyThreshold), ncol = length(lateThreshold))
type2Error <- matrix(data = NA, nrow = length(earlyThreshold), ncol = length(lateThreshold))
stopProbH0 <- array(data = NA, dim = c(length(earlyThreshold), length(lateThreshold), length(infoFraction)))
stopProbH1 <- array(data = NA, dim = c(length(earlyThreshold), length(lateThreshold), length(infoFraction)))
for (i in 1:length(earlyThreshold)) {
  for (j in 1:length(lateThreshold)) {
    config <- list(
      target = "op_char",
      n_threads = 8,
      effect_grid_size = 1001,
      arm_group_type = 0,
      prior = list(
        list(
          list(1, 1, 1)
        ),
        list(
          list(1, 1, 1)
        )
      ),
      allocation_ratio = as.list(allocationRatio),
      information_fracs = as.list(infoFraction),
      total_final_sample_size = sampleSize,
      theta = as.list(treatmentEffect),
      bin_size = list(1, 1, 1, 1, 1),
      prob_threshs = list(1e-08, 1e-08, 1e-08, 1e-08, 1e-08),
      criteria = list(
        list(0, "diff", 0, 1, list(-1.01, 0.00), ">",
             list(earlyThreshold[i], earlyThreshold[i],
                  lateThreshold[j], lateThreshold[j], lateThreshold[j]))
      ),
      decisions = list(
        list(1, "eff1", "post", list(0))
      ),
      print_info = FALSE
    )
    write_json(config, "config_custom.json", auto_unbox = TRUE, pretty = TRUE)

    system2(command = "./main",
            args = c("-c", "config_custom.json", "-o", "output_custom.json"))

    result <- fromJSON("output_custom.json")
    type1Error[i, j] <- result$p_eff_any_stage0
    type2Error[i, j] <- 1 - result$p_eff_any_stage1
    stopProbH0[i, j, ] <- unlist(lapply(X = result$p_efc0, FUN = function(x) {x[1, 2]}))
    stopProbH1[i, j, ] <- unlist(lapply(X = result$p_efc1, FUN = function(x) {x[1, 2]}))
  }
}
file.remove("config_custom.json", "output_custom.json")

save(sampleSize, treatmentEffect, allocationRatio,
     infoFraction, earlyThreshold, lateThreshold,
     type1Error, type2Error, stopProbH0, stopProbH1,
     file = "~/Dropbox/Jeffery He/iResearch/Publications/2024/ZH2026-ClinTrial-StatAdv-BayesGSD-AlphaSpending/Output/Output v1.0/HYPRESS_5Stage_convBayesGSDwR1.RData")

####################

# Figure A5
jpeg(file = "./fig321.jpeg", width = 12, height = 9, units = "in", res = 500)
par(mar = c(5.5, 5, 5.5, 5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
load("~/Dropbox/Jeffery He/iResearch/Publications/2024/ZH2026-ClinTrial-StatAdv-BayesGSD-AlphaSpending/Output/Output v1.0/HYPRESS_5Stage_convBayesGSDwR1.RData")
type1Error <- type1Error * 100
rownames(type1Error) <- formatC(x = earlyThreshold, digits = 4, format = "f")
colnames(type1Error) <- formatC(x = lateThreshold, digits = 4, format = "f")

zlim <- range(type1Error, na.rm = TRUE)
n_col <- 200
prop <- min(max((0.025 * 100 - zlim[1]) / diff(zlim), 0), 1)
n1 <- max(1, round(n_col * prop))
n2 <- max(1, n_col - n1)
cols <- c(colorRampPalette(c(pal_nejm("default", alpha = 1)(8)[2], "white"))(n1),
          colorRampPalette(c("white", pal_nejm("default", alpha = 1)(8)[1]))(n2))

image2D(z = type1Error,
        x = earlyThreshold,
        y = lateThreshold,
        col = cols,
        xlab = "Probability threshold for the early phase",
        ylab = "Probability threshold for the late phase",
        main = "Typer I Error Rate")
contour2D(z = type1Error,
          x = earlyThreshold,
          y = lateThreshold,
          col = "black",
          lwd = 1,
          add = TRUE)

# corrplot(corr = type1Error,
#          method = "square",
#          type = "full",
#          col = cols,
#          is.corr = FALSE,
#          title = "Type I Error Rate",
#          diag = TRUE,
#          mar = c(1, 1, 2, 0),
#          tl.col  = "black",
#          tl.srt  = 90,
#          cl.pos  = "r",
#          cl.lim  = zlim)
# mtext("Probability threshold for the late phase", side = 1, line = 4.0, cex = 1.1)
# mtext("Probability threshold for the early phase", side = 2, line = -1.5, cex = 1.1)
dev.off()

####################

# Figure A6
jpeg(file = "./fig322.jpeg", width = 12, height = 9, units = "in", res = 500)
par(mar = c(5.5, 5, 5.5, 5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
load("~/Dropbox/Jeffery He/iResearch/Publications/2024/ZH2026-ClinTrial-StatAdv-BayesGSD-AlphaSpending/Output/Output v1.0/HYPRESS_5Stage_convBayesGSDwR1.RData")
power <- (1 - type2Error) * 100
rownames(power) <- formatC(x = earlyThreshold, digits = 4, format = "f")
colnames(power) <- formatC(x = lateThreshold, digits = 4, format = "f")

zlim <- range(power, na.rm = TRUE)
n_col <- 200
prop <- min(max((0.800 * 100 - zlim[1]) / diff(zlim), 0), 1)
n1 <- max(1, round(n_col * prop))
n2 <- max(1, n_col - n1)
cols <- c(colorRampPalette(c(pal_nejm("default", alpha = 1)(8)[2], "white"))(n1),
          colorRampPalette(c("white", pal_nejm("default", alpha = 1)(8)[1]))(n2))

image2D(z = power,
        x = earlyThreshold,
        y = lateThreshold,
        col = cols,
        xlab = "Probability threshold for the early phase",
        ylab = "Probability threshold for the late phase",
        main = "Power")
contour2D(z = power,
          x = earlyThreshold,
          y = lateThreshold,
          col = "black",
          lwd = 1,
          add = TRUE)

# corrplot(corr = power,
#          method = "square",
#          type = "full",
#          col = cols,
#          is.corr = FALSE,
#          title = "Power",
#          diag = TRUE,
#          mar = c(1, 1, 2, 0),
#          tl.col  = "black",
#          tl.srt  = 90,
#          cl.pos  = "r",
#          cl.lim  = zlim)
# mtext("Probability threshold for the late phase", side = 1, line = 4.0, cex = 1.1)
# mtext("Probability threshold for the early phase", side = 2, line = -1.5, cex = 1.1)
dev.off()

####################

# Figure A7
jpeg(file = "./fig323.jpeg", width = 12, height = 9, units = "in", res = 500)
par(mar = c(5.5, 5, 5.5, 5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
load("~/Dropbox/Jeffery He/iResearch/Publications/2024/ZH2026-ClinTrial-StatAdv-BayesGSD-AlphaSpending/Output/Output v1.0/HYPRESS_5Stage_convBayesGSDwR1.RData")
expectedSampleSizeH0 <- matrix(data = NA, nrow = length(earlyThreshold), ncol = length(lateThreshold))
sampleSize <- ceiling(infoFraction * sampleSize / 2) * 2
for (i in 1:length(earlyThreshold)) {
  for (j in 1:length(lateThreshold)) {
    prob <- stopProbH0[i, j, 1:(length(infoFraction) - 1)]
    prob <- c(prob, 1 - sum(prob))
    expectedSampleSizeH0[i, j] <- sum(sampleSize * prob)
  }
}
rownames(expectedSampleSizeH0) <- formatC(x = earlyThreshold, digits = 4, format = "f")
colnames(expectedSampleSizeH0) <- formatC(x = lateThreshold, digits = 4, format = "f")

# zlim <- range(expectedSampleSizeH0, na.rm = TRUE)
# n_col <- 200
# prop <- min(max((0.025 * 100 - zlim[1]) / diff(zlim), 0), 1)
# n1 <- max(1, round(n_col * prop))
# n2 <- max(1, n_col - n1)
# cols <- c(colorRampPalette(c(pal_nejm("default", alpha = 1)(8)[2], "white"))(n1),
#           colorRampPalette(c("white", pal_nejm("default", alpha = 1)(8)[1]))(n2))

image2D(z = expectedSampleSizeH0,
        x = earlyThreshold,
        y = lateThreshold,
        col = colorRampPalette(pal_nejm("default")(8)[2:1])(200),
        xlab = "Probability threshold for the early phase",
        ylab = "Probability threshold for the late phase",
        main = "Expected sample size under the null")
contour2D(z = expectedSampleSizeH0,
          x = earlyThreshold,
          y = lateThreshold,
          col = "black",
          lwd = 1,
          add = TRUE)

# corrplot(corr = expectedSampleSizeH0,
#          method = "square",
#          type = "full",
#          col = cols,
#          is.corr = FALSE,
#          title = "Expected sample size under the null",
#          diag = TRUE,
#          mar = c(1, 1, 2, 0),
#          tl.col  = "black",
#          tl.srt  = 90,
#          cl.pos  = "r",
#          cl.lim  = zlim)
# mtext("Probability threshold for the late phase", side = 1, line = 4.0, cex = 1.1)
# mtext("Probability threshold for the early phase", side = 2, line = -1.5, cex = 1.1)
dev.off()

####################

# Figure A8
jpeg(file = "./fig324.jpeg", width = 12, height = 9, units = "in", res = 500)
par(mar = c(5.5, 5, 5.5, 5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
load("~/Dropbox/Jeffery He/iResearch/Publications/2024/ZH2026-ClinTrial-StatAdv-BayesGSD-AlphaSpending/Output/Output v1.0/HYPRESS_5Stage_convBayesGSDwR1.RData")
expectedSampleSizeH1 <- matrix(data = NA, nrow = length(earlyThreshold), ncol = length(lateThreshold))
sampleSize <- ceiling(infoFraction * sampleSize / 2) * 2
for (i in 1:length(earlyThreshold)) {
  for (j in 1:length(lateThreshold)) {
    prob <- stopProbH1[i, j, 1:(length(infoFraction) - 1)]
    prob <- c(prob, 1 - sum(prob))
    expectedSampleSizeH1[i, j] <- sum(sampleSize * prob)
  }
}
rownames(expectedSampleSizeH1) <- formatC(x = earlyThreshold, digits = 4, format = "f")
colnames(expectedSampleSizeH1) <- formatC(x = lateThreshold, digits = 4, format = "f")

# zlim <- range(expectedSampleSizeH1, na.rm = TRUE)
# n_col <- 200
# prop <- min(max((0.025 * 100 - zlim[1]) / diff(zlim), 0), 1)
# n1 <- max(1, round(n_col * prop))
# n2 <- max(1, n_col - n1)
# cols <- c(colorRampPalette(c(pal_nejm("default", alpha = 1)(8)[2], "white"))(n1),
#           colorRampPalette(c("white", pal_nejm("default", alpha = 1)(8)[1]))(n2))

image2D(z = expectedSampleSizeH1,
        x = earlyThreshold,
        y = lateThreshold,
        col = colorRampPalette(pal_nejm("default")(8)[2:1])(200),
        xlab = "Probability threshold for the early phase",
        ylab = "Probability threshold for the late phase",
        main = "Expected sample size under the alternative")
contour2D(z = expectedSampleSizeH1,
          x = earlyThreshold,
          y = lateThreshold,
          col = "black",
          lwd = 1,
          add = TRUE)

# corrplot(corr = expectedSampleSizeH1,
#          method = "square",
#          type = "full",
#          col = cols,
#          is.corr = FALSE,
#          title = "Expected sample size under the alternative",
#          diag = TRUE,
#          mar = c(1, 1, 2, 0),
#          tl.col  = "black",
#          tl.srt  = 90,
#          cl.pos  = "r",
#          cl.lim  = zlim)
# mtext("Probability threshold for the late phase", side = 1, line = 4.0, cex = 1.1)
# mtext("Probability threshold for the early phase", side = 2, line = -1.5, cex = 1.1)
dev.off()

########################################

load("~/Dropbox/Jeffery He/iResearch/Publications/2024/ZH2026-ClinTrial-StatAdv-BayesGSD-AlphaSpending/Output/Output v1.0/HYPRESS_5Stage_convBayesGSDwR1.RData")
index <- NULL
for (j in 1:length(lateThreshold)) {
  if (any(type1Error[, j] <= 0.025))
    index <- rbind(index, c(min(which(type1Error[, j] <= 0.025)), j))
}
# index <- index[which(type2Error[index] <= 0.200), ]

OC <- data.frame(earlyThreshold = earlyThreshold[index[, 1]],
                 lateThreshold = lateThreshold[index[, 2]],
                 type1Error = type1Error[index],
                 power = (1 - type2Error[index]))
print(OC)

stopProbH0 <- matrix(data = NA, nrow = nrow(index), ncol = length(infoFraction))
stopProbH1 <- matrix(data = NA, nrow = nrow(index), ncol = length(infoFraction))
for (i in 1:nrow(index)) {
  probThreshold <-
    c(earlyThreshold[index[i, 1]], earlyThreshold[index[i, 1]],
      lateThreshold[index[i, 2]], lateThreshold[index[i, 2]], lateThreshold[index[i, 2]])
  config <- list(
    target = "op_char",
    n_threads = 8,
    effect_grid_size = 1001,
    arm_group_type = 0,
    prior = list(
      list(
        list(1, 1, 1)
      ),
      list(
        list(1, 1, 1)
      )
    ),
    allocation_ratio = as.list(allocationRatio),
    information_fracs = as.list(infoFraction),
    total_final_sample_size = sampleSize,
    theta = as.list(treatmentEffect),
    bin_size = list(1, 1, 1, 1, 1),
    prob_threshs = list(1e-08, 1e-08, 1e-08, 1e-08, 1e-08),
    criteria = list(
      list(0, "diff", 0, 1, list(-1.01, 0.00), ">", as.list(probThreshold))
    ),
    decisions = list(
      list(1, "eff1", "post", list(0))
    ),
    print_info = FALSE
  )
  write_json(config, "config_custom.json", auto_unbox = TRUE, pretty = TRUE)

  system2(command = "./main",
          args = c("-c", "config_custom.json", "-o", "output_custom.json"))

  result <- fromJSON("output_custom.json")
  alpha <- result$p_eff_any_stage0

  while (alpha <= 0.025) {
    OC$earlyThreshold[i] <- probThreshold[1]
    OC$lateThreshold[i] <- probThreshold[3]
    OC$type1Error[i] <- result$p_eff_any_stage0
    OC$power[i] <- result$p_eff_any_stage1
    stopProbH0[i, ] <- unlist(lapply(X = result$p_efc0, FUN = function(x) {x[1, 2]}))
    stopProbH1[i, ] <- unlist(lapply(X = result$p_efc1, FUN = function(x) {x[1, 2]}))

    probThreshold <- probThreshold - c(0.0001, 0.0001, 0, 0, 0)
    config <- list(
      target = "op_char",
      n_threads = 8,
      effect_grid_size = 1001,
      arm_group_type = 0,
      prior = list(
        list(
          list(1, 1, 1)
        ),
        list(
          list(1, 1, 1)
        )
      ),
      allocation_ratio = as.list(allocationRatio),
      information_fracs = as.list(infoFraction),
      total_final_sample_size = sampleSize,
      theta = as.list(treatmentEffect),
      bin_size = list(1, 1, 1, 1, 1),
      prob_threshs = list(1e-08, 1e-08, 1e-08, 1e-08, 1e-08),
      criteria = list(
        list(0, "diff", 0, 1, list(-1.01, 0.00), ">", as.list(probThreshold))
      ),
      decisions = list(
        list(1, "eff1", "post", list(0))
      ),
      print_info = FALSE
    )
    write_json(config, "config_custom.json", auto_unbox = TRUE, pretty = TRUE)

    system2(command = "./main",
            args = c("-c", "config_custom.json", "-o", "output_custom.json"))

    result <- fromJSON("output_custom.json")
    alpha <- result$p_eff_any_stage0
  }

  if (OC$power[i] >= 0.80) {
    next
  } else {
    OC <- OC[1:(i - 1), ]
    stopProbH0 <- stopProbH0[1:(i - 1), ]
    stopProbH1 <- stopProbH1[1:(i - 1), ]
    break
  }
}
file.remove("config_custom.json", "output_custom.json")

save(sampleSize, OC, stopProbH0, stopProbH1,
     file = "~/Dropbox/Jeffery He/iResearch/Publications/2024/ZH2026-ClinTrial-StatAdv-BayesGSD-AlphaSpending/Output/Output v1.0/HYPRESS_5Stage_convBayesGSDwR1_OC.RData")

####################

# Figure 3
jpeg(file = "./fig325.jpeg", width = 12, height = 9, units = "in", res = 500)
par(mar = c(5.5, 5, 5.5, 5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
load("~/Dropbox/Jeffery He/iResearch/Publications/2024/ZH2026-ClinTrial-StatAdv-BayesGSD-AlphaSpending/Output/Output v1.0/HYPRESS_5Stage_convBayesGSDwR1_OC.RData")
cumSampleSize <- ceiling(infoFraction * sampleSize / 2) * 2
cumAlpha <- NULL
for (i in 1:nrow(stopProbH0)) {
  print(sum(stopProbH0[i, 1:4]))
  print(sum(stopProbH0[i, 1:4] * cumSampleSize[1:4]) + (1 - sum(stopProbH0[i, 1:4])) * cumSampleSize[5])
  cumAlpha <- rbind(cumAlpha, cumsum(c(0, stopProbH0[i, ])))
}
print(cumAlpha)
asP <- getDesignGroupSequential(kMax = 5,
                                alpha = 0.025,
                                beta = 0.2,
                                sided = 1,
                                typeOfDesign = "asP")
asOF <- getDesignGroupSequential(kMax = 5,
                                 alpha = 0.025,
                                 beta = 0.2,
                                 sided = 1,
                                 typeOfDesign = "asOF")

plot(range(c(0, infoFraction)), range(cumAlpha), type = "n",
     bty = "n", xaxt = "n", xlim = c(0, 1), ylim = c(0, 0.025),
     main = "Type I error spending", xlab = "Information rate", ylab = "Cumulative alpha spending")
axis(1, at = (0:5) / 5, labels = c("0", "1/5", "2/5", "3/5", "4/5", "1"))
for (i in 1:nrow(cumAlpha)) {
  lines(x = c(0, infoFraction), y = cumAlpha[i, ],
        type = "b", pch = 19, lwd = 3, col = adjustcolor(colorRampPalette(pal_nejm("default")(8)[1:2])(nrow(cumAlpha))[i], alpha.f = 0.25))
}
lines(x = c(0, infoFraction), y = cumAlpha[1, ],
      type = "b", pch = 19, lwd = 3, col = adjustcolor(colorRampPalette(pal_nejm("default")(8)[1:2])(nrow(cumAlpha))[1], alpha.f = 1))
lines(x = c(0, infoFraction), y = cumAlpha[nrow(cumAlpha), ],
      type = "b", pch = 19, lwd = 3, col = adjustcolor(colorRampPalette(pal_nejm("default")(8)[1:2])(nrow(cumAlpha))[nrow(cumAlpha)], alpha.f = 1))
lines(x = (0:5) / 5, y = c(0, asP$alphaSpent),
      type = "b", pch = 19, lwd = 3, col = pal_nejm("default", alpha = 1)(8)[3])
lines(x = (0:5) / 5, y = c(0, asOF$alphaSpent),
      type = "b", pch = 19, lwd = 3, col = pal_nejm("default", alpha = 1)(8)[4])
legend("topleft", legend = c("BayesGSD (strategy 1 with highest power)", "BayesGSD (strategy 1 with lowest power)", "Pocock-type spending", "O'Brien-Fleming-type spending"),
       col = pal_nejm("default", alpha = 1)(8)[1:4], lty = 1, lwd = 3, bty = "n", cex = 1.5)
dev.off()

####################

# Figure A9
jpeg(file = "./fig326.jpeg", width = 12, height = 9, units = "in", res = 500)
par(mar = c(5.5, 5, 5.5, 5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
load("~/Dropbox/Jeffery He/iResearch/Publications/2024/ZH2026-ClinTrial-StatAdv-BayesGSD-AlphaSpending/Output/Output v1.0/HYPRESS_5Stage_convBayesGSDwR1_OC.RData")
cumSampleSize <- ceiling(infoFraction * sampleSize / 2) * 2
cumPower <- NULL
for (i in 1:nrow(stopProbH1)) {
  print(sum(stopProbH1[i, 1:4]))
  print(sum(stopProbH1[i, 1:4] * cumSampleSize[1:4]) + (1 - sum(stopProbH1[i, 1:4])) * cumSampleSize[5])
  cumPower <- rbind(cumPower, cumsum(c(0, stopProbH1[i, ])))
}
print(cumPower)

plot(range(c(0, infoFraction)), range(cumPower), type = "n",
     bty = "n", xaxt = "n", xlim = c(0, 1), ylim = c(0, max(cumPower)),
     main = "Cumulative power", xlab = "Information rate", ylab = "Cumulative power")
axis(1, at = (0:5) / 5, labels = c("0", "1/5", "2/5", "3/5", "4/5", "1"))
for (i in 1:nrow(cumPower)) {
  lines(x = c(0, infoFraction), y = cumPower[i, ],
        type = "b", pch = 19, lwd = 3, col = adjustcolor(colorRampPalette(pal_nejm("default")(8)[1:2])(nrow(cumPower))[i], alpha.f = 0.25))
}
lines(x = c(0, infoFraction), y = cumPower[1, ],
      type = "b", pch = 19, lwd = 3, col = adjustcolor(colorRampPalette(pal_nejm("default")(8)[1:2])(nrow(cumPower))[1], alpha.f = 1))
lines(x = c(0, infoFraction), y = cumPower[nrow(cumPower), ],
      type = "b", pch = 19, lwd = 3, col = adjustcolor(colorRampPalette(pal_nejm("default")(8)[1:2])(nrow(cumPower))[nrow(cumPower)], alpha.f = 1))
legend("topleft", legend = c("BayesGSD (strategy 1 with highest power)", "BayesGSD (strategy 1 with lowest power)"),
       col = pal_nejm("default", alpha = 1)(8)[1:4], lty = 1, lwd = 3, bty = "n", cex = 1.5)
dev.off()

############################################################

# conventional Bayesian GSD with 4 interim analyses (refinement 2: posterior vs predictive)
sampleSize <- 368
treatmentEffect <- c(0.40, 0.25)
allocationRatio <- c(1, 1)
infoFraction <- (1:5) / 5
postThreshold <- seq(from = 0.975, to = 0.992, by = 0.0005)
predThreshold <- seq(from = 0.825, to = 0.992, by = 0.0005)
type1Error <- matrix(data = NA, nrow = length(postThreshold), ncol = length(predThreshold))
type2Error <- matrix(data = NA, nrow = length(postThreshold), ncol = length(predThreshold))
stopProbH0 <- array(data = NA, dim = c(length(postThreshold), length(predThreshold), length(infoFraction)))
stopProbH1 <- array(data = NA, dim = c(length(postThreshold), length(predThreshold), length(infoFraction)))
for (i in 1:length(postThreshold)) {
  for (j in 1:length(predThreshold)) {
    config <- list(
      target = "op_char",
      n_threads = 8,
      effect_grid_size = 1001,
      arm_group_type = 0,
      prior = list(
        list(
          list(1, 1, 1)
        ),
        list(
          list(1, 1, 1)
        )
      ),
      allocation_ratio = as.list(allocationRatio),
      information_fracs = as.list(infoFraction),
      total_final_sample_size = sampleSize,
      theta = as.list(treatmentEffect),
      bin_size = list(1, 1, 1, 1, 1),
      prob_threshs = list(1e-16, 1e-16, 1e-16, 1e-16, 1e-16),
      criteria = list(
        list(0, "diff", 0, 1, list(-1.01, 0.00), ">",
             list(predThreshold[j], predThreshold[j],
                  predThreshold[j], predThreshold[j], postThreshold[i]))
      ),
      decisions = list(
        list(1, "eff1", "pred", list(0))
      ),
      print_info = FALSE
    )
    write_json(config, "config_custom.json", auto_unbox = TRUE, pretty = TRUE)

    system2(command = "./main",
            args = c("-c", "config_custom.json", "-o", "output_custom.json"))

    result <- fromJSON("output_custom.json")
    type1Error[i, j] <- result$p_eff_any_stage0
    type2Error[i, j] <- 1 - result$p_eff_any_stage1
    stopProbH0[i, j, ] <- unlist(lapply(X = result$p_efc0, FUN = function(x) {x[1, 2]}))
    stopProbH1[i, j, ] <- unlist(lapply(X = result$p_efc1, FUN = function(x) {x[1, 2]}))
  }
}
file.remove("config_custom.json", "output_custom.json")

save(sampleSize, treatmentEffect, allocationRatio,
     infoFraction, postThreshold, predThreshold,
     type1Error, type2Error, stopProbH0, stopProbH1,
     file = "~/Dropbox/Jeffery He/iResearch/Publications/2024/ZH2026-ClinTrial-StatAdv-BayesGSD-AlphaSpending/Output/Output v1.0/HYPRESS_5Stage_convBayesGSDwR2.RData")

####################

# Figure A10
jpeg(file = "./fig341.jpeg", width = 12, height = 9, units = "in", res = 500)
par(mar = c(5.5, 5, 5.5, 5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
load("~/Dropbox/Jeffery He/iResearch/Publications/2024/ZH2026-ClinTrial-StatAdv-BayesGSD-AlphaSpending/Output/Output v1.0/HYPRESS_5Stage_convBayesGSDwR2.RData")
type1Error <- type1Error * 100
rownames(type1Error) <- formatC(x = postThreshold, digits = 4, format = "f")
colnames(type1Error) <- formatC(x = predThreshold, digits = 4, format = "f")

zlim <- range(type1Error, na.rm = TRUE)
n_col <- 200
prop <- min(max((0.025 * 100 - zlim[1]) / diff(zlim), 0), 1)
n1 <- max(1, round(n_col * prop))
n2 <- max(1, n_col - n1)
cols <- c(colorRampPalette(c(pal_nejm("default", alpha = 1)(8)[2], "white"))(n1),
          colorRampPalette(c("white", pal_nejm("default", alpha = 1)(8)[1]))(n2))

image2D(z = type1Error,
        x = postThreshold,
        y = predThreshold,
        col = cols,
        xlab = "Probability threshold for the posterior rule",
        ylab = "Probability threshold for the predictive rule",
        main = "Typer I Error Rate")
contour2D(z = type1Error,
          x = postThreshold,
          y = predThreshold,
          col = "black",
          lwd = 1,
          add = TRUE)

# corrplot(corr = type1Error,
#          method = "square",
#          type = "full",
#          col = cols,
#          is.corr = FALSE,
#          title = "Type I Error Rate",
#          diag = TRUE,
#          mar = c(1, 1, 2, 0),
#          tl.col  = "black",
#          tl.srt  = 90,
#          cl.pos  = "r",
#          cl.lim  = zlim)
# mtext("Probability threshold for the predictive rule", side = 1, line = 4.0, cex = 1.1)
# mtext("Probability threshold for the posterior rule", side = 2, line = -1.5, cex = 1.1)
dev.off()

####################

# Figure A11
jpeg(file = "./fig342.jpeg", width = 12, height = 9, units = "in", res = 500)
par(mar = c(5.5, 5, 5.5, 5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
load("~/Dropbox/Jeffery He/iResearch/Publications/2024/ZH2026-ClinTrial-StatAdv-BayesGSD-AlphaSpending/Output/Output v1.0/HYPRESS_5Stage_convBayesGSDwR2.RData")
power <- (1 - type2Error) * 100
rownames(power) <- formatC(x = postThreshold, digits = 4, format = "f")
colnames(power) <- formatC(x = predThreshold, digits = 4, format = "f")

zlim <- range(power, na.rm = TRUE)
n_col <- 200
prop <- min(max((0.800 * 100 - zlim[1]) / diff(zlim), 0), 1)
n1 <- max(1, round(n_col * prop))
n2 <- max(1, n_col - n1)
cols <- c(colorRampPalette(c(pal_nejm("default", alpha = 1)(8)[2], "white"))(n1),
          colorRampPalette(c("white", pal_nejm("default", alpha = 1)(8)[1]))(n2))

image2D(z = power,
        x = postThreshold,
        y = predThreshold,
        col = cols,
        xlab = "Probability threshold for the posterior rule",
        ylab = "Probability threshold for the predictive rule",
        main = "Power")
contour2D(z = power,
          x = postThreshold,
          y = predThreshold,
          col = "black",
          lwd = 1,
          add = TRUE)

# corrplot(corr = power,
#          method = "square",
#          type = "full",
#          col = cols,
#          is.corr = FALSE,
#          title = "Power",
#          diag = TRUE,
#          mar = c(1, 1, 2, 0),
#          tl.col  = "black",
#          tl.srt  = 90,
#          cl.pos  = "r",
#          cl.lim  = zlim)
# mtext("Probability threshold for the predictive rule", side = 1, line = 4.0, cex = 1.1)
# mtext("Probability threshold for the posterior rule", side = 2, line = -1.5, cex = 1.1)
dev.off()

####################

# Figure A12
jpeg(file = "./fig343.jpeg", width = 12, height = 9, units = "in", res = 500)
par(mar = c(5.5, 5, 5.5, 5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
load("~/Dropbox/Jeffery He/iResearch/Publications/2024/ZH2026-ClinTrial-StatAdv-BayesGSD-AlphaSpending/Output/Output v1.0/HYPRESS_5Stage_convBayesGSDwR2.RData")
expectedSampleSizeH0 <- matrix(data = NA, nrow = length(postThreshold), ncol = length(predThreshold))
sampleSize <- ceiling(infoFraction * sampleSize / 2) * 2
for (i in 1:length(postThreshold)) {
  for (j in 1:length(predThreshold)) {
    prob <- stopProbH0[i, j, 1:(length(infoFraction) - 1)]
    prob <- c(prob, 1 - sum(prob))
    expectedSampleSizeH0[i, j] <- sum(sampleSize * prob)
  }
}
rownames(expectedSampleSizeH0) <- formatC(x = postThreshold, digits = 4, format = "f")
colnames(expectedSampleSizeH0) <- formatC(x = predThreshold, digits = 4, format = "f")

# zlim <- range(expectedSampleSizeH0, na.rm = TRUE)
# n_col <- 200
# prop <- min(max((0.025 * 100 - zlim[1]) / diff(zlim), 0), 1)
# n1 <- max(1, round(n_col * prop))
# n2 <- max(1, n_col - n1)
# cols <- c(colorRampPalette(c(pal_nejm("default", alpha = 1)(8)[2], "white"))(n1),
#           colorRampPalette(c("white", pal_nejm("default", alpha = 1)(8)[1]))(n2))

image2D(z = expectedSampleSizeH0,
        x = postThreshold,
        y = predThreshold,
        col = colorRampPalette(pal_nejm("default")(8)[2:1])(200),
        xlab = "Probability threshold for the posterior rule",
        ylab = "Probability threshold for the predictive rule",
        main = "Expected sample size under the null")
contour2D(z = expectedSampleSizeH0,
          x = postThreshold,
          y = predThreshold,
          col = "black",
          lwd = 1,
          add = TRUE)

# corrplot(corr = expectedSampleSizeH0,
#          method = "square",
#          type = "full",
#          col = cols,
#          is.corr = FALSE,
#          title = "Expected sample size under the null",
#          diag = TRUE,
#          mar = c(1, 1, 2, 0),
#          tl.col  = "black",
#          tl.srt  = 90,
#          cl.pos  = "r",
#          cl.lim  = zlim)
# mtext("Probability threshold for the predictive rule", side = 1, line = 4.0, cex = 1.1)
# mtext("Probability threshold for the posterior rule", side = 2, line = -1.5, cex = 1.1)
dev.off()

####################

# Figure A13
jpeg(file = "./fig344.jpeg", width = 12, height = 9, units = "in", res = 500)
par(mar = c(5.5, 5, 5.5, 5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
load("~/Dropbox/Jeffery He/iResearch/Publications/2024/ZH2026-ClinTrial-StatAdv-BayesGSD-AlphaSpending/Output/Output v1.0/HYPRESS_5Stage_convBayesGSDwR2.RData")
expectedSampleSizeH1 <- matrix(data = NA, nrow = length(postThreshold), ncol = length(predThreshold))
sampleSize <- ceiling(infoFraction * sampleSize / 2) * 2
for (i in 1:length(postThreshold)) {
  for (j in 1:length(predThreshold)) {
    prob <- stopProbH1[i, j, 1:(length(infoFraction) - 1)]
    prob <- c(prob, 1 - sum(prob))
    expectedSampleSizeH1[i, j] <- sum(sampleSize * prob)
  }
}
rownames(expectedSampleSizeH1) <- formatC(x = postThreshold, digits = 4, format = "f")
colnames(expectedSampleSizeH1) <- formatC(x = predThreshold, digits = 4, format = "f")

# zlim <- range(expectedSampleSizeH1, na.rm = TRUE)
# n_col <- 200
# prop <- min(max((0.025 * 100 - zlim[1]) / diff(zlim), 0), 1)
# n1 <- max(1, round(n_col * prop))
# n2 <- max(1, n_col - n1)
# cols <- c(colorRampPalette(c(pal_nejm("default", alpha = 1)(8)[2], "white"))(n1),
#           colorRampPalette(c("white", pal_nejm("default", alpha = 1)(8)[1]))(n2))

image2D(z = expectedSampleSizeH1,
        x = postThreshold,
        y = predThreshold,
        col = colorRampPalette(pal_nejm("default")(8)[2:1])(200),
        xlab = "Probability threshold for the posterior rule",
        ylab = "Probability threshold for the predictive rule",
        main = "Expected sample size under the alternative")
contour2D(z = expectedSampleSizeH1,
          x = postThreshold,
          y = predThreshold,
          col = "black",
          lwd = 1,
          add = TRUE)

# corrplot(corr = expectedSampleSizeH1,
#          method = "square",
#          type = "full",
#          col = cols,
#          is.corr = FALSE,
#          title = "Expected sample size under the alternative",
#          diag = TRUE,
#          mar = c(1, 1, 2, 0),
#          tl.col  = "black",
#          tl.srt  = 90,
#          cl.pos  = "r",
#          cl.lim  = zlim)
# mtext("Probability threshold for the predictive rule", side = 1, line = 4.0, cex = 1.1)
# mtext("Probability threshold for the posterior rule", side = 2, line = -1.5, cex = 1.1)
dev.off()

########################################

load("~/Dropbox/Jeffery He/iResearch/Publications/2024/ZH2026-ClinTrial-StatAdv-BayesGSD-AlphaSpending/Output/Output v1.0/HYPRESS_5Stage_convBayesGSDwR2.RData")
index <- NULL
for (i in 1:length(postThreshold)) {
  if (any(type1Error[i, ] <= 0.025))
    index <- rbind(index, c(i, min(which(type1Error[i, ] <= 0.025))))
}
# index <- index[which(type2Error[index] <= 0.200), ]

OC <- data.frame(postThreshold = postThreshold[index[, 1]],
                 predThreshold = predThreshold[index[, 2]],
                 type1Error = type1Error[index],
                 power = (1 - type2Error[index]))
print(OC)

stopProbH0 <- matrix(data = NA, nrow = nrow(index), ncol = length(infoFraction))
stopProbH1 <- matrix(data = NA, nrow = nrow(index), ncol = length(infoFraction))
for (i in 1:nrow(index)) {
  probThreshold <-
    c(predThreshold[index[i, 2]], predThreshold[index[i, 2]],
      predThreshold[index[i, 2]], predThreshold[index[i, 2]], postThreshold[index[i, 1]])
  config <- list(
    target = "op_char",
    n_threads = 8,
    effect_grid_size = 1001,
    arm_group_type = 0,
    prior = list(
      list(
        list(1, 1, 1)
      ),
      list(
        list(1, 1, 1)
      )
    ),
    allocation_ratio = as.list(allocationRatio),
    information_fracs = as.list(infoFraction),
    total_final_sample_size = sampleSize,
    theta = as.list(treatmentEffect),
    bin_size = list(1, 1, 1, 1, 1),
    prob_threshs = list(1e-16, 1e-16, 1e-16, 1e-16, 1e-16),
    criteria = list(
      list(0, "diff", 0, 1, list(-1.01, 0.00), ">", as.list(probThreshold))
    ),
    decisions = list(
      list(1, "eff1", "pred", list(0))
    ),
    print_info = FALSE
  )
  write_json(config, "config_custom.json", auto_unbox = TRUE, pretty = TRUE)

  system2(command = "./main",
          args = c("-c", "config_custom.json", "-o", "output_custom.json"))

  result <- fromJSON("output_custom.json")
  alpha <- result$p_eff_any_stage0

  while (alpha <= 0.025) {
    OC$postThreshold[i] <- probThreshold[5]
    OC$predThreshold[i] <- probThreshold[1]
    OC$type1Error[i] <- result$p_eff_any_stage0
    OC$power[i] <- result$p_eff_any_stage1
    stopProbH0[i, ] <- unlist(lapply(X = result$p_efc0, FUN = function(x) {x[1, 2]}))
    stopProbH1[i, ] <- unlist(lapply(X = result$p_efc1, FUN = function(x) {x[1, 2]}))

    probThreshold <- probThreshold - c(0.0001, 0.0001, 0.0001, 0.0001, 0)
    config <- list(
      target = "op_char",
      n_threads = 8,
      effect_grid_size = 1001,
      arm_group_type = 0,
      prior = list(
        list(
          list(1, 1, 1)
        ),
        list(
          list(1, 1, 1)
        )
      ),
      allocation_ratio = as.list(allocationRatio),
      information_fracs = as.list(infoFraction),
      total_final_sample_size = sampleSize,
      theta = as.list(treatmentEffect),
      bin_size = list(1, 1, 1, 1, 1),
      prob_threshs = list(1e-16, 1e-16, 1e-16, 1e-16, 1e-16),
      criteria = list(
        list(0, "diff", 0, 1, list(-1.01, 0.00), ">", as.list(probThreshold))
      ),
      decisions = list(
        list(1, "eff1", "pred", list(0))
      ),
      print_info = FALSE
    )
    write_json(config, "config_custom.json", auto_unbox = TRUE, pretty = TRUE)

    system2(command = "./main",
            args = c("-c", "config_custom.json", "-o", "output_custom.json"))

    result <- fromJSON("output_custom.json")
    alpha <- result$p_eff_any_stage0
  }

  if (OC$power[i] >= 0.80) {
    next
  } else {
    OC <- OC[1:(i - 1), ]
    stopProbH0 <- stopProbH0[1:(i - 1), ]
    stopProbH1 <- stopProbH1[1:(i - 1), ]
    break
  }
}
file.remove("config_custom.json", "output_custom.json")

save(sampleSize, OC, stopProbH0, stopProbH1,
     file = "~/Dropbox/Jeffery He/iResearch/Publications/2024/ZH2026-ClinTrial-StatAdv-BayesGSD-AlphaSpending/Output/Output v1.0/HYPRESS_5Stage_convBayesGSDwR2_OC.RData")

####################

# Figure 4
jpeg(file = "./fig345.jpeg", width = 12, height = 9, units = "in", res = 500)
par(mar = c(5.5, 5, 5.5, 5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
load("~/Dropbox/Jeffery He/iResearch/Publications/2024/ZH2026-ClinTrial-StatAdv-BayesGSD-AlphaSpending/Output/Output v1.0/HYPRESS_5Stage_convBayesGSDwR2_OC.RData")
cumSampleSize <- ceiling(infoFraction * sampleSize / 2) * 2
cumAlpha <- NULL
for (i in 1:nrow(stopProbH0)) {
  print(sum(stopProbH0[i, 1:4]))
  print(sum(stopProbH0[i, 1:4] * cumSampleSize[1:4]) + (1 - sum(stopProbH0[i, 1:4])) * cumSampleSize[5])
  cumAlpha <- rbind(cumAlpha, cumsum(c(0, stopProbH0[i, ])))
}
print(cumAlpha)
asP <- getDesignGroupSequential(kMax = 5,
                                alpha = 0.025,
                                beta = 0.2,
                                sided = 1,
                                typeOfDesign = "asP")
asOF <- getDesignGroupSequential(kMax = 5,
                                 alpha = 0.025,
                                 beta = 0.2,
                                 sided = 1,
                                 typeOfDesign = "asOF")

plot(range(c(0, infoFraction)), range(cumAlpha), type = "n",
     bty = "n", xaxt = "n", xlim = c(0, 1), ylim = c(0, 0.025),
     main = "Type I error spending", xlab = "Information rate", ylab = "Cumulative alpha spending")
axis(1, at = (0:5) / 5, labels = c("0", "1/5", "2/5", "3/5", "4/5", "1"))
for (i in 1:nrow(cumAlpha)) {
  lines(x = c(0, infoFraction), y = cumAlpha[i, ],
        type = "b", pch = 19, lwd = 3, col = adjustcolor(colorRampPalette(pal_nejm("default")(8)[1:2])(nrow(cumAlpha))[i], alpha.f = 0.25))
}
lines(x = c(0, infoFraction), y = cumAlpha[1, ],
      type = "b", pch = 19, lwd = 3, col = adjustcolor(colorRampPalette(pal_nejm("default")(8)[1:2])(nrow(cumAlpha))[1], alpha.f = 1))
lines(x = c(0, infoFraction), y = cumAlpha[nrow(cumAlpha), ],
      type = "b", pch = 19, lwd = 3, col = adjustcolor(colorRampPalette(pal_nejm("default")(8)[1:2])(nrow(cumAlpha))[nrow(cumAlpha)], alpha.f = 1))
lines(x = (0:5) / 5, y = c(0, asP$alphaSpent),
      type = "b", pch = 19, lwd = 3, col = pal_nejm("default", alpha = 1)(8)[3])
lines(x = (0:5) / 5, y = c(0, asOF$alphaSpent),
      type = "b", pch = 19, lwd = 3, col = pal_nejm("default", alpha = 1)(8)[4])
legend("topleft", legend = c("BayesGSD (strategy 2 with highest power)", "BayesGSD (strategy 2 with lowest power)", "Pocock-type spending", "O'Brien-Fleming-type spending"),
       col = pal_nejm("default", alpha = 1)(8)[1:4], lty = 1, lwd = 3, bty = "n", cex = 1.5)
dev.off()

####################

# Figure A15
jpeg(file = "./fig346.jpeg", width = 12, height = 9, units = "in", res = 500)
par(mar = c(5.5, 5, 5.5, 5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
load("~/Dropbox/Jeffery He/iResearch/Publications/2024/ZH2026-ClinTrial-StatAdv-BayesGSD-AlphaSpending/Output/Output v1.0/HYPRESS_5Stage_convBayesGSDwR2_OC.RData")
cumSampleSize <- ceiling(infoFraction * sampleSize / 2) * 2
cumPower <- NULL
for (i in 1:nrow(stopProbH1)) {
  print(sum(stopProbH1[i, 1:4]))
  print(sum(stopProbH1[i, 1:4] * cumSampleSize[1:4]) + (1 - sum(stopProbH1[i, 1:4])) * cumSampleSize[5])
  cumPower <- rbind(cumPower, cumsum(c(0, stopProbH1[i, ])))
}
print(cumPower)

plot(range(c(0, infoFraction)), range(cumPower), type = "n",
     bty = "n", xaxt = "n", xlim = c(0, 1), ylim = c(0, max(cumPower)),
     main = "Cumulative power", xlab = "Information rate", ylab = "Cumulative power")
axis(1, at = (0:5) / 5, labels = c("0", "1/5", "2/5", "3/5", "4/5", "1"))
for (i in 1:nrow(cumPower)) {
  lines(x = c(0, infoFraction), y = cumPower[i, ],
        type = "b", pch = 19, lwd = 3, col = adjustcolor(colorRampPalette(pal_nejm("default")(8)[1:2])(nrow(cumPower))[i], alpha.f = 0.25))
}
lines(x = c(0, infoFraction), y = cumPower[1, ],
      type = "b", pch = 19, lwd = 3, col = adjustcolor(colorRampPalette(pal_nejm("default")(8)[1:2])(nrow(cumPower))[1], alpha.f = 1))
lines(x = c(0, infoFraction), y = cumPower[nrow(cumPower), ],
      type = "b", pch = 19, lwd = 3, col = adjustcolor(colorRampPalette(pal_nejm("default")(8)[1:2])(nrow(cumPower))[nrow(cumPower)], alpha.f = 1))
legend("topleft", legend = c("BayesGSD (strategy 2 with highest power)", "BayesGSD (strategy 2 with lowest power)"),
       col = pal_nejm("default", alpha = 1)(8)[1:4], lty = 1, lwd = 3, bty = "n", cex = 1.5)
dev.off()

################################################################################
