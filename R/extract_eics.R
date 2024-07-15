# This script include some needed functions to extract EIC (they are from DecoMetDIA)

#' plotFeatureEIC function
#'
#' @param xset xcms_object
#' @param is.check.eic Boolean
#' @param extra integer
#'
#' @importFrom xcms groupval
#'
#' @return EIC plot
#'
plotFeatureEIC <- function(xset, is.check.eic = TRUE, extra = 30) {
  rt.range.data <- range(xset@rt)
  
  rtmin.origin <- unname(apply(groupval(xset, value = 'rtmin'), 1, median))
  rtmax.origin <- unname(apply(groupval(xset, value = 'rtmax'), 1, median))
  
  rtmin <- rtmin.origin - extra
  rtmax <- rtmax.origin + extra
  
  rtmin[which(rtmin < rt.range.data[1])] <- rt.range.data[1]
  rtmax[which(rtmax > rt.range.data[2])] <- rt.range.data[2]
  
  rt.range.eic <- cbind( rtmin, rtmax)
  eics <- getEIC(xset, rtrange = rt.range.eic, sampleidx = seq(length(xset@filepaths)), groupidx = seq(nrow(rt.range.eic)))
  
  cat('Ploting EICs ... ')
  smpnames <- names(eics@eic)
  clr <- rainbow(length(smpnames), end = 0.85)
  rgbvec <- pmin(col2rgb(clr) + 153, 255)
  clr_gray <- rgb(rgbvec[1, ], rgbvec[2, ], rgbvec[3, ], max = 255)
  
  
  dirEics <- file.path('featureEICs', c('all', smpnames))
  names(dirEics) <- c('all', smpnames)
  sapply(dirEics, function(fp) {
    if (!file.exists(fp)) {
      dir.create(fp, recursive = TRUE)
    }
  })
  
  fn.all <- file.path(dirEics['all'], paste0(eics@groupnames, '.png'))
  
  if (length(smpnames) > 1) {
    sapply(seq_along(eics@groupnames), function(idx.ft) {
      max.int <- max(unname(sapply(smpnames, function(smpname) {
        max(eics@eic[[smpname]][[idx.ft]][, 2])
      })))
      png(file = fn.all[idx.ft], width = 960, height = 480)
      plot(0, 0, type = "n", xlim = rt.range.eic[idx.ft,], ylim = c(0, max.int),
           xlab = "Retention Time (seconds)", ylab = "Intensity"
      )
      sapply(seq_along(smpnames), function(idx.smp) {
        ft <- eics@eic[[smpnames[idx.smp]]][[idx.ft]]
        lines(ft, col = clr_gray[idx.smp])
        idx.ft.peak <- which(ft[, 1] >= rtmin.origin[idx.ft] & ft[, 1] <= rtmax.origin[idx.ft])
        lines(ft[idx.ft.peak, ], col = clr[idx.smp])
        if (is.check.eic) {
          abline(v = c(rtmin.origin[idx.ft], rtmax.origin[idx.ft]), col = clr[idx.smp])
        }
      })
      legend('topright', smpnames, col = clr, lty = 1)
      title(eics@groupnames[idx.ft])
      par(new = FALSE)
      dev.off()
    })
  }
  
  sapply(seq_along(smpnames), function(idx.smp) {
    sapply(seq_along(eics@groupnames), function(idx.ft) {
      png(file = file.path(dirEics[smpnames[idx.smp]], paste0(eics@groupnames[idx.ft], '.png')), width = 960, height = 480)
      ft <- eics@eic[[smpnames[idx.smp]]][[idx.ft]]
      plot(ft, col = 'black')
      lines(ft, col = clr_gray[1])
      idx.ft.peak <- which(ft[, 1] >= rtmin.origin[idx.ft] & ft[, 1] <= rtmax.origin[idx.ft])
      lines(ft[idx.ft.peak, ], col = clr[1])
      if (is.check.eic) {
        abline(v = c(rtmin.origin[idx.ft], rtmax.origin[idx.ft]), col = clr[1])
      }
      title(eics@groupnames[idx.ft])
      dev.off()
    })
  })
  return(eics)
}

#-------------------------------------------------------------------------------

#' get TICs
#'
#' @param xsetc I should check
#' @param fn.tic character pdf name
#' @param rt character c("raw", "corrected")
#'
#' @return TICs
#'
GetTICs <- function(xsetc = NULL, fn.tic = "TICs.pdf", rt = c("raw", "corrected")) {
  fn.smps <- xsetc@filepaths
  num.smp <- length(fn.smps)
  TIC <- vector('list', num.smp)
  for (i in 1:num.smp) {
    cat(fn.smps[i], '\n')
    if (!is.null(xsetc) && rt == 'corrected') {
      rtcor <- xsetc@rt$corrected[[i]]
    }
    else rtcor <- NULL
    TIC[[i]] <- GetTIC(fn.smps[i], rtcor = rtcor)
  }
  pdf(fn.tic, width = 16, height = 10)
  cols <- rainbow(num.smp)
  lty = 1:num.smp
  pch = 1:num.smp
  xlim = range(sapply(TIC, function(x) range(x[, 1], na.rm = T)))
  ylim = range(sapply(TIC, function(x) range(x[, 2], na.rm = T)))
  plot(0, 0,
       type = "n",
       xlim = xlim,
       ylim = ylim,
       main = "Total Ion Chromatograms",
       xlab = "Retention Time",
       ylab = "TIC")
  for (i in 1:num.smp) {
    tic <- TIC[[i]]
    points(tic[, 1], tic[, 2], col = cols[i], pch = pch[i],
           type = "l")
  }
  legend("topright", paste(basename(fn.smps)),
         col = cols, lty = lty, pch = pch)
  dev.off()
  invisible(TIC)
}

#-------------------------------------------------------------------------------

#' get TIC
#'
#' @param fn character files paths list
#' @param rtcor character c("raw", "corrected")
#'
#' @return TIC
#'
GetTIC <- function (fn, rtcor = NULL) {
  xraw <- xcmsRaw(fn)
  if (is.null(rtcor)) {
    rt <- xraw@scantime
  }
  else {
    rt <- rtcor
  }
  int <- rawEIC(xraw, mzrange = range(xraw@env$mz))$intensity
  cbind(rt, int)
}


