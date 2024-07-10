# this script to get the dia features as done by decometdia

# modified or extracted functions from Decometdia

#' get dia features from a mzML file
#'
#' @param obj new generated dia.Feature object
#' @param files list_of_strings mzML files paths
#' @param d.out path results path
#' @param polarity string 'POS' or 'NEG'
#' @param lc string Liquid chromatoghraphy column used for chromatograms separation
#' @param ppm.pd integer max mz devaition in consecutive scans
#' @param sn integer singl to noise ration
#' @param peakwidth integer range of the chromatogram peak width in seconds
#' @param mzdiff integer
#' @param minfrac integer
#' @param method.rector character
#' @param is.ms1.filter Boolean
#' @param thr.ms1.filter Boolean or NULL
#' @param filed.ms1.filter character
#' @param is.plot.eic.feature Boolean
#' @param is.check.eic Boolean
#' @param nSlaves integer
#' @param bw integer
#' @param mzwid integer
#'
#' @importFrom xcms xcmsSet groups groupval
#' @importFrom CAMERA xsAnnotate groupFWHM findIsotopes findAdducts getPeaklist
#' @importFrom BiocParallel SnowParam
#' @importFrom utils write.csv
#'
#' @return dia.feature.RData and many other xcms objects
#'
GetFeatures <- function(obj, files, d.out,
           polarity = "positive", lc = "RPLC",
           ppm.pd = 6, sn = 0, peakwidth = c(5,30),
           mzdiff = 0.0003, minfrac = 0.5, method.rector = "obiwarp",
           is.ms1.filter = FALSE, thr.ms1.filter = NULL,
           filed.ms1.filter = "into",
           is.plot.eic.feature = TRUE, is.check.eic = FALSE,
           nSlaves = 6, bw = 5, mzwid = 0.015) {
    filed.ms1.filter <- match.arg(filed.ms1.filter)
    method.rector <- match.arg(method.rector)
    is.xcms3 <- FALSE
    if (packageVersion("xcms") >= '1.50.1') {
      is.xcms3 <- TRUE
      # require(BiocParallel)
      param <- SnowParam(workers = nSlaves, type = 'SOCK')
    }

    fn.lc.list <- paste0(lc, '_POS.csv')
    rules.camera <- read.csv(file.path(system.file(package = "DecoMetDIA", 'rules'), fn.lc.list))

    cat('Detect peaks with xcms ...\n')

    if (is.xcms3) {
      xset <- xcmsSet(files, method = 'centWave', ppm = ppm.pd, snthr = sn, peakwidth = peakwidth, mzdiff = mzdiff, BPPARAM = param)
    } else {
      xset <- xcmsSet(files,  method = 'centWave', ppm = ppm.pd, snthr = sn, peakwidth = peakwidth, mzdiff = mzdiff, nSlaves = nSlaves)
    }

    save(xset, file = paste0(d.out, '/xset-centWave.Rda'))
    sclassv <- as.character(sampclass(xset))
    # output the detected peak num
    for (i in 1:length( xset@filepaths)) {
      cat( i, xset@filepaths[ i], "\t [",  sclassv[ i], "]", sep = "")
      cat( "  --> ", length( which( xset@peaks[, "sample"] == i)), " Features. \n")
    }
    # output sample group information
    cat("\n Sample classes -> ", levels(sampclass(xset)), "\n")


    # require(CAMERA)
    if (length(files) == 1) {
      ftname <- apply(xset@peaks, 1, function(dr) {
        paste("M", sprintf('%.0f', dr['mz']),
              "T", sprintf('%.0f', dr['rt']),
              sep = "")
      })
      is.duplicated <- TRUE
      while (is.duplicated) {
        i <- 2
        idx.duplicated <- which(duplicated(ftname))
        if (length(idx.duplicated) > 0) {
          ftname[idx.duplicated] <- paste0(ftname[idx.duplicated], '_', i)
          i <- i + 1
        } else {
          is.duplicated <- FALSE
        }
      }

      pt <- cbind('name' = ftname, xset@peaks)
      col.nm <- colnames(pt)
      col.nm[c(2,5)] <- c('mzmed', 'rtmed')
      colnames(pt) <- (col.nm)

      write.csv(pt, file = paste0( d.out, "/PeakTable-raw.csv"), row.names = T)

      # for one sample only, no peak alignment
      xa <- xsAnnotate(xset, polarity = polarity)
      xa <- groupFWHM(xa)
      xa <- findIsotopes(xa)
      xa <- findAdducts(xa, rules = rules.camera, polarity = polarity)
      peaklist.filtered <- peaklist.anno <- cbind('name' = ftname,
                                                  getPeaklist(xa, intval = filed.ms1.filter))
      col.nm <- colnames(peaklist.anno)
      print(col.nm)
      col.nm[c(2,5)] <- c('mzmed', 'rtmed')
      colnames(peaklist.anno) <- col.nm
      colnames(peaklist.filtered) <- col.nm

      write.csv(peaklist.anno, paste0( d.out, '/PeakTable-annotated.csv'),
                row.names = FALSE)

      nr.remove <- NULL
      nr.isotope <- which(grepl('\\]\\[M\\+\\d+', peaklist.anno[, 'isotopes']))
      nr.ms1.filter <- which(xset@peaks[, filed.ms1.filter] < thr.ms1.filter)
      if (is.ms1.filter & (length(nr.ms1.filter) > 0 | length(nr.isotope) > 0)) {
        nr.remove <- unique(c(nr.isotope, nr.ms1.filter))
        xset@peaks <- xset@peaks[-nr.remove, , drop = FALSE]
        peaklist.filtered <- peaklist.anno[-nr.remove, , drop = FALSE]
      }
      peaklist.filtered <- cbind('ft.idx' = seq(nrow(peaklist.filtered)), peaklist.filtered)
      write.csv(peaklist.filtered, paste0( d.out, '/PeakTable-annotated_filtered.csv'),
                row.names = FALSE)

      tmp <- lapply(seq(nrow(xset@peaks)), function(idx.r) {
        r <- xset@peaks[idx.r, ]
        r.gp <- c(r[1:6], 1, 1)
        names(r.gp) <- c("mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax", "npeaks", "10")
        r.gp
      } )

      xset@groups <- do.call(rbind,tmp)
      xset@groupidx <- lapply(seq(nrow(xset@groups)), function(i) {
        i
      })
      if (is.null(nr.remove)) {
        obj@peaks <- data.frame(xset@peaks)
      } else {
        obj@peaks <- data.frame(xset@peaks[-nr.isotope, , drop = FALSE])
      }

      obj@peakgroup$raw <- lapply(seq(nrow(obj@peaks)),
                                  function(nr) {
                                    obj@peaks[nr, , drop = FALSE]
                                  })
      obj@peakgroup$corrected <- obj@peakgroup$raw
      obj@samplenames <- sampnames(xset)
      obj@filepaths <- xset@filepaths
      obj@smpclass <- sclassv
      if(is.plot.eic.feature) {
        eics <- plotFeatureEIC(xset, is.check.eic = is.check.eic, extra = 30)
      }
      xset3 <- xset
      save(xset3, file = paste0( d.out, '/xset-centWave-final.Rda'))
    } else {
      # for multiple samples, peak alignment

      xset <- group(xset, minfrac = minfrac)
      pdf(paste0('retcor-', method.rector, '.pdf'))
      xset2 <- switch(method.rector,
                      'obiwarp' = retcor(xset, method = method.rector, plottype = 'deviation', profStep = 0.1),
                      'peakgroups' = retcor(xset, method = method.rector, plottype = 'deviation')
      )
      dev.off()

      CatLineSeparator("Creating retention time corrected TIC's ...")
      TICS.rtcor <- GetTICs(xset2, fn.tic = 'TICs-rtcor.pdf',
                            rt = 'corrected')  # total ion current
      save(xset2, TICS.rtcor, file = paste0( d.out, '/xset-centWave-retcor.Rda'))

      xset2 <- group(xset2, bw = bw, mzwid = mzwid, minfrac = minfrac)

      groupmat <- groups(xset2)
      pt <- data.frame(cbind('name' = groupnames(xset2), groupmat,
                             groupval(xset2, "medret", 'into')),
                       row.names = NULL)
      write.csv(pt, file = paste0( d.out, "/PeakTable-raw.csv"), row.names = T)
      gc()

      cat(nrow(groupmat), "aligned features. \n")

      xset3 <- fillPeaks(xset2)

      xa <- xsAnnotate(xset3, polarity= polarity, nSlaves = nSlaves)
      xa <- groupFWHM(xa)
      xa <- findIsotopes(xa)
      xa <- findAdducts(xa, rules = rules.camera, polarity = polarity)
      pt <- data.frame(cbind('name' = groupnames(xset2), groupmat,
                             groupval(xset3, "medret", 'into')),
                       row.names = NULL)
      pt.filter <- groupval(xset3, value = filed.ms1.filter)
      colnames(pt.filter) <- paste0(filed.ms1.filter, '_',
                                    make.names(colnames(pt.filter)))
      rownames(pt.filter) <- NULL
      pt.anno <- getPeaklist(xa, intval = filed.ms1.filter)[, c('isotopes', 'adduct', 'pcgroup')]
      peaklist.filtered <- peaklist.anno <- cbind(pt, pt.filter, pt.anno)
      nr.isotope <- which(grepl('\\]\\[M\\+\\d+', peaklist.anno[, 'isotopes']))
      write.csv(peaklist.anno, paste0( d.out, '/PeakTable-annotated.csv'), row.names = FALSE)

      value.ms1.filter <- unname(apply(pt.filter, 1, max))
      nr.ms1.filter <- which(value.ms1.filter < thr.ms1.filter)
      if (is.ms1.filter & (length(nr.ms1.filter) > 0 | length(nr.isotope) > 0)) {
        nr.remove <- unique(c(nr.isotope, nr.ms1.filter))
        xset3@groups <- xset3@groups[-nr.remove, , drop = FALSE]
        xset3@groupidx <- xset3@groupidx[-nr.remove]
        peaklist.filtered <- peaklist.anno[-nr.remove, , drop = FALSE]
      }
      peaklist.filtered <- cbind('ft.idx' = seq(nrow(peaklist.filtered)), peaklist.filtered)
      write.csv(peaklist.filtered, paste0( d.out, '/PeakTable-annotated_filtered.csv'),
                row.names = FALSE)

      obj@peaks <- data.frame(xset3@peaks)

      grpIdx.pre <- xset3@groupidx
      grpIdx <- lapply(seq_along(grpIdx.pre), function(idx) {
        Index <- grpIdx.pre[[idx]][order(abs(obj@peaks[grpIdx.pre[[idx]], 'rt'] - median(obj@peaks[grpIdx.pre[[idx]], 'rt'])))]
        Index[match(seq_along(sampnames(xset)), obj@peaks[Index,'sample'])]
      })

      obj@peakgroup$corrected <- lapply(grpIdx, function(idx) {
        obj@peaks[idx, , drop = FALSE]
      })

      obj@rt$corrected <- xset3@rt$corrected
      obj@rt$raw <- xset3@rt$raw

      obj@peakgroup$raw <- lapply(obj@peakgroup$corrected, function(peak) {
        peak.raw.list <- lapply(1:nrow(peak), function(y) {
          rt.idx <- which.min(abs(peak[y, "rt"] - obj@rt$corrected[[peak[y, "sample"]]]))
          peak[y, "rt"] <- obj@rt$raw[[peak[y, "sample"]]][rt.idx]
          rtmin.idx <- which.min(abs(peak[y, "rtmin"] - obj@rt$corrected[[peak[y, "sample"]]]))
          peak[y, "rtmin"] <- obj@rt$raw[[peak[y, "sample"]]][rtmin.idx]
          rtmax.idx <- which.min(abs(peak[y, "rtmax"] - obj@rt$corrected[[peak[y, "sample"]]]))
          peak[y, "rtmax"] <- obj@rt$raw[[peak[y, "sample"]]][rtmax.idx]
          peak[y, ]
        })
        peak.raw <- do.call(rbind, peak.raw.list)
      })

      obj@samplenames <- sampnames(xset3)
      obj@filepaths <- xset3@filepaths
      obj@smpclass <- sclassv
      if(is.plot.eic.feature) {
        eics <- plotFeatureEIC(xset3, is.check.eic = is.check.eic, extra = 30)
      }

      save(xset3, file = paste0( d.out, '/xset-centWave-final.Rda'))
    }
    return(obj)
  }

#-------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------

#' separate info by line
#'
#' @param info character info to be printed
#'
#' @return print info separeted by line
#'
CatLineSeparator <- function(info) {
  info.sep <- paste(rep('-', 38), collapse  = '')
  cat('', info.sep, info, info.sep, '', sep = '\n')
}

#-------------------------------------------------------------------------------

# diaSpectra
# diaMethods

# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

#' get dia_features and process the data using xcms
#'
#' @param d.in character
#' @param d.out path results path
#' MS1 peaks detection parameters
#' @param nSlaves integer
#' @param peakwidth integer range of the chromatogram peak width in seconds
#' @param sn integer singl to noise ration
#' @param ppm.pd integer max mz devaition in consecutive scans
#' @param polarity string 'POS' or 'NEG'
#' @param lc string Liquid chromatoghraphy column used for chromatograms separation
#' @param minfrac integer
#' @param bw integer
#' @param mzwid integer
#' @param mzdiff integer
#' @param is.ms1.filter Boolean
#' @param thr.ms1.filter Boolean or NULL
#' @param filed.ms1.filter character
#' @param is.plot.eic.feature Boolean
#' @param is.check.eic Boolean
#' MS2 peaks detections and processed using the SWATH (dictionary) file
#' @param is.PlotEIC.MS2 Boolean
#' @param snthr.ms2 integer
#' @param fp.swath.setup SWATH file path
#' @param ppm.ms2.mtp integer
#' @param int.filter integer
#' @param minfrac.vote integer
#' @param is.plot.eic.spec Boolean
#' @param is.output.spec.smp Boolean
#' @param is.span.constant Boolean
#' @param rerun Boolean
#' @param rerun.ms2 Boolean
#' @param is.ms1.only Boolean
#' @param is.dec.all Boolean
#' @param is.dec.smoothed Boolean
#' @param is.norm.consensus Boolean
#' Parameter I added it, if you have several files, you can choose to process them in correlation or every file separately
#' @param correlation Boolean
#'
#' @importFrom methods new
#' @importFrom parallel detectCores
#' @importFrom BiocParallel SnowParam
#' @importFrom tools file_path_sans_ext
#'
#' @return dia.features and xcms objects
#'
DecoMetDIA_modified_funct <- function(d.in,
                       d.out, nSlaves = 6,
                       # peak detection setup
                       peakwidth = c(5, 30), mzdiff = 0.001, sn = 0, ppm.pd = 6,
                       polarity = 'positive', lc = 'RPLC',  minfrac = 0.5,
                       bw = 5, mzwid = 0.015, is.ms1.filter = FALSE, thr.ms1.filter = NULL,
                       filed.ms1.filter = c('maxo', 'into', 'intb'), is.plot.eic.feature = F, is.check.eic = FALSE,
                       # SWATH MS2 setup
                       is.PlotEIC.MS2 = FALSE, snthr.ms2 = 5, fp.swath.setup,
                       ppm.ms2.mtp = 15, int.filter = 50, minfrac.vote = 0.5,
                       is.plot.eic.spec = F, is.output.spec.smp = TRUE, is.span.constant = FALSE,
                       rerun = FALSE, rerun.ms2 = FALSE, is.ms1.only = FALSE,
                       is.dec.all = FALSE, is.dec.smoothed = FALSE, is.norm.consensus = TRUE,
                       # I added a parameter to process the files at separately or in correlation
                       correlation = TRUE ) {
  span <- 0.3
  mz.diff <- 0.005
  wd0 <- getwd();

  files <- list.files(d.in, recursive = TRUE, full.names = TRUE, pattern = '(?i)mzML$')

  if(correlation == TRUE ){
    d.out2 <- paste0(wd0, '/', d.out);
    if (!dir.exists(d.out2)) {
      dir.create(d.out2, recursive = TRUE)
    };

    nSlaves <- min(detectCores() - 1, nSlaves, length(file))
    bpparam <- SnowParam(workers = nSlaves, type = 'SOCK')
    CatLineSeparator('Detecting and aligning features ...')
    fn.skip <- paste0(d.out2 ,'/dia.feature.RData')

    if ((!rerun) & file.exists(fn.skip)) {
      cat('using existing results:', fn.skip, '...\n')
      load(fn.skip)
    } else {
      dia.feature <- GetFeatures(new('diaFeature'), files, d.out2,
                                 polarity = polarity, lc = lc,
                                 ppm.pd = ppm.pd, sn = sn, peakwidth = peakwidth,
                                 mzdiff = mz.diff, minfrac = minfrac,
                                 method.rector = 'obiwarp',
                                 is.ms1.filter = is.ms1.filter,
                                 thr.ms1.filter = thr.ms1.filter,
                                 filed.ms1.filter = filed.ms1.filter,
                                 is.plot.eic.feature = is.plot.eic.feature,
                                 is.check.eic = is.check.eic,
                                 nSlaves = nSlaves, bw = bw, mzwid = mzwid)

      save(dia.feature, file = fn.skip)
    }

  } else {
    for (i in 1:length(files)) {
      one_file <- files[[i]]
      d.out2 <- paste0(wd0, '/', d.out, '/', file_path_sans_ext(basename(one_file)))
      if (!dir.exists(d.out2)) {
        dir.create(d.out2, recursive = TRUE)
      };

      nSlaves <- min(detectCores() - 1, nSlaves, length(file))
      bpparam <- SnowParam(workers = nSlaves, type = 'SOCK')
      CatLineSeparator('Detecting and aligning features ...')
      fn.skip <- paste0( d.out2 ,'/dia.feature.RData')

      if ((!rerun) & file.exists(fn.skip)) {
        cat('using existing results:', fn.skip, '...\n')
        load(fn.skip)
      } else {
        dia.feature <- GetFeatures(new('diaFeature'), one_file, d.out2,
                                   polarity = polarity, lc = lc,
                                   ppm.pd = ppm.pd, sn = sn, peakwidth = peakwidth,
                                   mzdiff = mz.diff, minfrac = minfrac,
                                   method.rector = 'obiwarp',
                                   is.ms1.filter = is.ms1.filter,
                                   thr.ms1.filter = thr.ms1.filter,
                                   filed.ms1.filter = filed.ms1.filter,
                                   is.plot.eic.feature = is.plot.eic.feature,
                                   is.check.eic = is.check.eic,
                                   nSlaves = nSlaves, bw = bw, mzwid = mzwid)

        save(dia.feature, file = fn.skip)
      }
    }

  }

  # if(is.ms1.only) {
  #   return(NULL)
  # }
  # dt.peaktable <- read.csv('PeakTable-annotated_filtered.csv', stringsAsFactors = FALSE)
  #
  # cat('ALL peak detection work are finished!!\n')
  #
  # t.start <- Sys.time()
  # CatLineSeparator('Start processing SWATH data...')
  #
  #
  # spectra.all <- bplapply(files, function(fp.smp) {
  #   # fp.smp <- files
  #   print(fp.smp)
  #   suppressMessages(require(DecoMetDIA))
  #   spec.dia <- DecoMetDIA::GetDIASpectra(new('diaSpectra'), fp.smp, dia.feature,
  #                                         dt.peaktable=dt.peaktable,
  #                                         d.spec = file.path(d.out, 'SpecDIA'),
  #                                         peakwidth = peakwidth, snthr = snthr.ms2,
  #                                         int.filter = int.filter, ppm.ms2.mtp = ppm.ms2.mtp,
  #                                         is.PlotEIC.MS2 = is.PlotEIC.MS2,
  #                                         span = span,
  #                                         rerun = rerun.ms2,
  #                                         is.span.constant = is.span.constant,
  #                                         fp.swath = fp.swath.setup,
  #                                         is.dec.all = is.dec.all,
  #                                         is.dec.smoothed = is.dec.smoothed)
  #   return(spec.dia)
  # },BPPARAM = bpparam)
  # save(spectra.all, file = 'dataSpectraDIA.Rda')
  #
  # spec.consensus <- ConsensusSpec(spectra.all, dia.feature@smpclass,
  #                                 nrow(dt.peaktable), minfrac.vote,
  #                                 is.norm.consensus = is.norm.consensus)
  # names(spec.consensus) <- dt.peaktable$name
  # spec.consensus <- spec.consensus[!sapply(spec.consensus, is.null)]
  # sink(file.path(d.out, 'spec_consensus.msp'))
  # for (nm.ft in names(spec.consensus)) {
  #   for (nm.grp in unique(dia.feature@smpclass)) {
  #     spec <- spec.consensus[[nm.ft]][[nm.grp]]
  #     if (!is.null(spec)) {
  #       OutputSpec(spec, which(dt.peaktable$name == nm.ft), dt.peaktable, grp = nm.grp)
  #     }
  #   }
  # }
  # sink()
  #
  # if (is.output.spec.smp) {
  #   names(spectra.all) <- files
  #   for (fp.smp in files) {
  #     idx.smp <- which(dia.feature@filepaths %in% fp.smp)
  #     print(idx.smp)
  #     if (length(idx.smp) == 0) {
  #       if (file.exists(fp.smp)) {
  #         idx.smp <- which(basename(dia.feature@filepaths) %in% basename(fp.smp))
  #       } else {
  #         stop(paste0("file not exist: \n  ", fp.smp))
  #       }
  #     }
  #
  #     d.spec <- file.path(file.path(d.out, 'SpecDIA'), dia.feature@smpclass[idx.smp],  gsub('(?i).mzxml', '', basename(fp.smp)))
  #     if (!dir.exists(d.spec)) {
  #       dir.create(d.spec, recursive = TRUE)
  #     }
  #     info.all <- spectra.all[[fp.smp]]
  #     for(nm.spec in do.call(rbind, lapply(info.all, names))[1,]) {
  #       fn.spec <- file.path(d.spec, paste0(nm.spec ,'.msp'))
  #       sink(fn.spec)
  #       for (idx.pg in seq_along(info.all)) {
  #         if (!is.null(info.all[[idx.pg]])) {
  #           spec.output <- info.all[[idx.pg]][[nm.spec]]
  #           if (!is.null(spec.output)) {
  #           }
  #           OutputSpec(spec.output, idx.pg, dt.peaktable)
  #         }
  #       }
  #       sink()
  #     }
  #   }
  # }
  #
  # t.end <- Sys.time()
  # cat('ALL WORK DONE!!')
  # cat(t.start, '\n', t.end, '\n')
  # setwd(wd0)
  # return(file.path(d.out, 'SpecDIA'))

}
