# Scott Pletcher's Life Table Function

LT <-function (data, window = 2) 
{
  ## Here data is three columns from FormatDataFromTreatment 
  ## Age, Events, EventType
  ## This will only accept a matrix with three
  ## columns.  
  ## Three column format is (time,number of events, event).
  ## event = 0 is alive and event = 1 is dead.
  ## Assumes censored are right censored.
  ## Note that the first observation point (setup) must be included
  ## in the data file.  Note also that by definition there
  ## can be no events at this time (the time number of events must
  ## be zero for this interval.
  
  ## First we need to determine at what points the observations
  ## were made.  Thus, we need to keep the zeros, and pick out
  ## unique entries
  
  tmp <- data[data[, 3] == 1, ] ##Parse out lines that show deaths
  times <- sort(unique(tmp[, 1])) ##Get time points when deaths happened
  tmp <- data[data[, 3] == 0, ] ##Parse out lines that show censored
  censored.times <- sort(unique(tmp[, 1])) ##Get time points when censor happened
  
  if (length(censored.times) == 0) ##See if there is censor. 
    any.censored <- FALSE
  else any.censored <- TRUE
  
  ## Again, the first observation is when the cohorts were established
  ## There is one less interval than age entry.
  ## Deaths observed at time t_i are assumed to occur over the
  ## interval t_(i-1) - > t_i, while censored observations
  ## entered at time t_i are assumed to have occured at time
  ## t_i, meaning these individuals were at risk over the interval
  ## t_(i-1) -> t_i.
  ##Note that this previous definition of censoring applies to flies, which
  ##are lost at the time of measuring mortality, versus worms, which are
  ##lost before they are measured.
  
  interval.starts <- times[-(length(times))] ##All time points except for the last one
  interval.ends <- times[-1] ##All time points except for the first one
  interval.starts <- c(interval.starts, interval.ends[length(interval.ends)]) ##This is for calculating interval widths
  interval.ends <- c(interval.ends, Inf) ##This is for calculating interval widths
  interval.widths <- interval.ends - interval.starts
  interval.midpoints <- (interval.ends + interval.starts)/2
  
  ## Now make the dx in terms of intervals.
  ## Remember according to our data structure the first
  ## entry is the set up date.  Thus, the second entry
  ## in this column reflects the number of deaths in the first
  ## interval, and so on...  The second entry in the data is
  ## the first element of interval.ends, so we will use that
  ## as the marker.
  
  ## Group the deaths together in case there are
  ## multiple lines per age.
  
  dx <- rep(0, length(interval.ends))
  for (i in 1:length(interval.ends)) {
    tmp <- data[data[, 1] == interval.ends[i] & data[, 3] == 
                  1, 2]
    dx[i] <- sum(tmp)
  }
  
  ## Same for censored.  We don't require an entry for
  ## observations of zero censored observations (althoug
  ## it is okay if they are included), but we
  ## do require that the censored observations occur on a
  ## day when the cohorts were observed for deaths (although
  ## there need not be any deaths on this day, it must be
  ## included in the data set with zero as the entry).
  
  cx <- rep(0, length(interval.ends))
  if (any.censored) {
    for (i in censored.times) {
      if (i != (min(data[, 1]))) {
        if (sum(i == interval.ends) != 1) {
          print("ERROR!!  There is not an entered interval for an observed censored event!!")
          print("Is this because there is an entry for censored data on the setup day??")
          return(i)
        }
        tmp <- data[data[, 1] == i & data[, 3] == 0, 
                    2]
        cx[interval.ends == i] <- sum(tmp)
      }
    }
  }
  
  ## ###################################################
  ## Now remove all trailing zeros
  
  max.t1 <- 1:length(dx)
  max.t2 <- 1:length(cx)
  chooser <- dx > 0
  if(sum(chooser)>0) {
    max.t1 <- max.t1[chooser]
    max.t1 <- max(max.t1)
  }
  else {
    max.t1<- -99
  }
  chooser <- cx > 0
  if(sum(chooser)>0) {
    max.t2 <- max.t2[chooser]
    max.t2 <- max(max.t2)
  }
  else
    max.t2<- -99
  max.t3 <- max(c(max.t1, max.t2))
  
  ## we will still have an error here if there are no deaths and 
  ## no censored, i.e., if there are chambers with no animals.
  
  ## Remeber to leave one last zero to define the end
  ## of the final interval!
  
  if (max.t3 < max(length(dx), length(cx))) {
    max.t3 <- max.t3 + 1
  }
  
  interval.starts <- interval.starts[1:max.t3]
  interval.ends <- interval.ends[1:max.t3]
  interval.midpoints <- interval.midpoints[1:max.t3]
  interval.widths <- interval.widths[1:max.t3]
  
  dx <- dx[1:max.t3]
  cx <- cx[1:max.t3]
  
  ## ne_i is the number entering the ith interval.
  ## ne_i=ne_(i-1)-d_(i-1)-c_(i-1).
  ## This is due to the way we record our data.  The number of deaths
  ## recorded on day 'i' imply these individuals died in the
  ## previous interval (t_(i-1)->t_i).  However, the number of
  ## censored indiviudals on t_i reflect individuals that are lost
  ## at that exact time, and therefore were at risk in the previous
  ## interval, but not the subsequent one.  This is accounted for
  ## in the cx vector, which is aligned with interval.ends, as opposed
  ## to interval.starts like dx.
  
  ne <- rep(0, length(interval.starts))
  ne[1] <- sum(c(dx, cx))
  for (i in 2:length(ne)) {
    ne[i] <- ne[i - 1] - dx[i - 1] - cx[i - 1]
  }
  nrisk <- ne - (0.5) * dx
  nrisk <- ne
  
  ## nrisk is the number at risk in the ith interval.
  ## individuals dying in the interval are assumed to have done so
  ## half-way through.  Note that we do not adjust for censoring here
  ## because we know exactly when the censoring happend, thus censored
  ## individuals are exposed throughout the previous interval, but not
  ## at all in the current one. 
  
  qx <- dx/nrisk
  px <- 1 - qx
  
  ## Now get Sx from px.  Note that this is an estimate
  ## of the survival function at time t_i, thus the point times
  ## for this is interval.starts, not interval.ends.
  
  Sx <- rep(0, length(interval.starts))
  Sx[1] <- 1
  
  ## The rest of the measures deal with intervals and do not extend beyond 
  ## the observed interval ends.
  for (i in 2:length(interval.starts)) Sx[i] <- Sx[i - 1] * 
    px[i - 1]
  Lx <- (ne - (0.5 * dx)) * interval.widths
  Lx[is.nan(Lx)] <- 0
  Tx <- rev(cumsum(rev(Lx)))
  ex <- Tx/ne
  fx <- (Sx * qx)/interval.widths
  hx2 <- (2 * qx)/(interval.widths * (1 + px))
  hx1 <- ((-1) * log(px))/interval.widths
  
  tmp <- (qx/(nrisk * px))
  tmp <- cumsum(tmp)
  tmp2 <- Sx * Sx
  seSx <- sqrt(tmp2 * tmp)
  tmp <- (qx/(nrisk * px))
  tmp <- cumsum(tmp)
  tmp2 <- ((Sx * qx) * (Sx * qx))/interval.widths
  tmp3 <- px/(nrisk * qx)
  sefx <- rep(0, length(fx))
  sefx[1] <- tmp2[1] * tmp3[1]
  for (i in 2:(length(fx))) sefx[i] <- sqrt(tmp2[i] * (tmp[i - 
                                                             1] + tmp3[i]))
  tmp <- 0.5 * hx1 * interval.widths
  tmp <- (1 - tmp * tmp)
  tmp2 <- (hx1 * hx1)/(nrisk * qx)
  sehx1 <- tmp * tmp2
  tmp <- 0.5 * hx2 * interval.widths
  tmp <- (1 - tmp * tmp)
  tmp2 <- (hx2 * hx2)/(nrisk * qx)
  sehx2 <- tmp * tmp2
  lnhx1 <- log(hx1)
  lnhx2 <- log(hx2)
  upper.95ci.hx1 <- log(hx1 + 1.6 * sehx1)
  lower.95ci.hx1 <- log(hx1 - 1.6 * sehx1)
  upper.95ci.hx2 <- log(hx2 + 1.6 * sehx2)
  lower.95ci.hx2 <- log(hx2 - 1.6 * sehx2)
  indx <- (-1) * length(hx1)
  indx <- c(indx, indx + 1)
  x <- interval.midpoints[is.finite(hx1)]
  y <- hx1[is.finite(hx1)]
  smoothlnux <- lowess(x, y, f = 1/6, delta = 0.1)
  tmp.int <- rep(NA, length(interval.midpoints))
  names(tmp.int) <- interval.midpoints
  tmp.int[as.character(smoothlnux$x)] <- smoothlnux$y
  smoothlnux <- log(c(smoothlnux$y, NA, NA))
  smoothlnux <- log(tmp.int)
  results <- data.frame(interval.starts, interval.ends, interval.widths, 
                        dx, cx, ne, nrisk, qx, px, ex, Sx, seSx, interval.midpoints, 
                        fx, sefx, hx1, sehx1, hx2, sehx2, lnhx1, upper.95ci.hx1, 
                        lower.95ci.hx1, lnhx2, upper.95ci.hx2, lower.95ci.hx2, 
                        smoothlnux)
  names(results) <- c("IntStart", "IntEnd", "IntWidth", "dx", 
                      "cx", "ne", "nrisk", "qx", "px", "ex", "Sx", "seSx", 
                      "IntMid", "fx", "sefx", "hx", "sehx", "hx2", "sehx2", 
                      "ln", "seUCI", "seLCI", "lnux2", "U95CIlnux2", "L95CIlnux2", 
                      "Sm")
  ## For now I will remove the last row because
  ## it is generally an infinite interval width.
  ##results[-nrow(results),]
  
  ## Actually, i will keep the last row because it has the last of the 
  ## survivorship numbers.
  results
}
