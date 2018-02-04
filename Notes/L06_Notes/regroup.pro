FUNCTION regroup, Y, X, DX=dx, DY=dy, SORT=sort, BINWIDTH=binwidth, $
                MINBIN=minbin, VARERR=varerr, BIN_DY=bin_dy, BIN_N=bin_n, $
                BIN_X=bin_x, BIN_L=bin_l, BIN_U=bin_u, bin_map=bin_map, $
                GEOM=geom, SNR=snr

; ----------------------------------------------------------
;+
; NAME:
;       REGROUP
;
; PURPOSE:
;       Rebin (regroup) binned/spaced one dimensional data
;
; AUTHOR:
;       Simon Vaughan (U.Leicester) 
;
; CALLING SEQUENCE:
;       result = REGROUP(y, x, binwidth = -1.5)
;
; INPUTS:
;      x - vector containing abcissae.
;      y - vector containing ordinates.
;       
;
; OPTIONAL INPUTS:  
;      DX        - (float) spacing of x positions
;      DY        - (vector) errors on ordinates
;      BINWIDTH  - (float) bin width (default=5.0)
;      MINBIN    - (integer) minimum no. data points per bin
;                    (default=1) 
;      VARERR    - (logical) whether to calculate error based
;                    on sqrt(variance) within bin (default=no)
;      SORT      - (logical) sort into accending order of x?
;      GEOM      - (logical) use geometric centre of bin? (detault=no)
;      SNR       - (float) specific a minimum signal/noise ratio
;
; OUTPUTS:
;      BIN_y     - (vector)binned ordinates
;
; OPTIONAL OUTPUTS:  
;      BIN_x   - binned abcissae
;      BIN_l   - lower edge of abcissae bin
;      BIN_u   - upper edge of abcissae bin
;      BIN_dy  - error on ordinate bin
;      BIN_n   - no. points contained in bin
;      BIN_map - mapping of output bins for each input bin
;
; DETAILS:
;      'regroups' evenly spaced or binned data into larger bins. If the
;      data comprise y[i] (with i=0,1,2,...,N-1) taken at evenly spaced
;      x[i] then GROUP calculates the average y over wider intervals
;      of x. If the x are sampled every dx = 5 we might want to
;      'group' these into bins of width 15 using the BINWIDTH=15
;      keyword to set the width of the output bins. 
;
;      Unlike the rebinning performed by LIN_REBIN and LOG_REBIN,
;      'grouping' respects the input bin widths. We assume the data are
;      already binned into contiguous bins of width dx, and the output
;      bins are simply 'groupings' of the existing bins. Another way
;      to think about this is that the output bins contain
;      integer numbers of complete input bins, the start and end
;      points of output bins always fall at the (existing) boundaries
;      of input bins. 
;
;      There are three choices for grouping: linear, logarithmic and
;      adaptive. For linear grouping use BINWIDTH > 0. For example,
;      BINWIDTH = 10 will regroup the data into bins of width 10 (or similar if
;      BINWIDTH is not an exactly multiple of DX). For logarithmic
;      grouping using BINWIDTH < -1. For example, BINWIDTH = -1.2 will
;      regroup the data every factor of 1.2, i.e. bins run roughly X
;      --> 1.2*X. As output bins are forced to completely contain in
;      integer number (>=MINBIN) of input bins, the bins may become
;      much larger than this where 1.2*X is smaller than DX.
;
;      'Adaptive' regrouping means the data are group together such
;      that the signal/noise in each output bin is roughly constant,
;      so the output bin widths may vary in width - being narrow where
;      the signal-to-noise is high in the input, and broad where it is
;      low. To use adaptive regrouping set the SNR keyword to the
;      desired Signal/Noise Ratio. For example, usinh SNR=5.0 will
;      regroup the data to value[j]/error[j] >= 5.0 in each bin j,
;      where value[j] is the mean of the input data points
;      contributing to output bin j, and error[j] is the standard
;      error on the mean.
;
;      If the VARERR keyword is set the error is calculated from the
;      standard error on the mean of the data within the bin
;      (i.e. error = sqrt[variance]/sqrt[n]).  
;      The bin locations (along X) are defined as follows.
;      The lower and upper bin boundaries - BIN_L and BIN_U - should
;      fall on esiting boundaries of the input bins, and these output
;      bins should be contiguous, i.e. BIN_U[i] = BIN_L[i+1]. For
;      linear binning the bin centre BIN_X is the half-way point of
;      the bin interval; for logarithmic binning it is the geometrical
;      half-way point of the bin interval. 
;
; EXAMPLE USAGE:
;
;      n = 400
;      x = INDGEN(n) + 1 
;      y = 100*exp(-x/80.0) 
;      y = POIDEV(y, SEED=seed)
;      dy = SQRT(y)
;      plot, x, y, /ylog, yrange=[0.1,200], psym=10, /xstyle, /ystyle
;      
;      bin_y = REGROUP(y, x, dy=dy, dx=1.0, bin_x=bin_x, bin_l=bin_l, $
;      		     bin_u=bin_u, snr=10, bin_dy=bin_dy, $
;      		     bin_map=bin_map)
;      bin_y = REGROUP(y, x, dy=dy, dx=1.0, bin_x=bin_x, bin_l=bin_l, $
;      		     binwidth=-1.2, bin_u=bin_u, bin_dy=bin_dy, $
;      		     bin_map=bin_map)      
;      bin_y = REGROUP(y, x, /varerr, binwidth=-1.2, dx=1.0, bin_x=bin_x, $
;      		     bin_l=bin_l, bin_u=bin_u, bin_dy=bin_dy, $
;      		     bin_map=bin_map)      
;
;      for i = 0,N_ELEMENTS(bin_x)-1 do plots,[bin_l[i],bin_u[i]],[bin_y[i],bin_y[i]], color=150, thick=2.0
;      for i = 0,N_ELEMENTS(bin_x)-2 do plots,[bin_u[i],bin_u[i]],[bin_y[i],bin_y[i+1]], color=150, thick=2.0     
;      plot_err, bin_x, bin_y, bin_dy, color=150, thick=2.0

;
; HISTORY:
;      24/09/2010 - v1.0 - first working version
;      01/10/2010 - v1.1 - added GEOM keyword
;      10/03/2011 - v1.2 - renamed REGROUP. Added 'adaptive' binning
;                           using new SNR keyword.
;      16/07/2012 - v1.3 - Added check on LAST_BIN to watch-out for a 
;                           glitch in the lowest bin.
;
; NOTES:
;      Although this is slightly bad practice, the input arguments X
;      and Y are sorted into accending X order during running of REGROUP,
;      if the SORT keyword is set. 
;
;      We currently assume that the data are evenly sampled/binned
;      with a bin spacing of DX, with no ties in the X values. 
;
;      WARNING: still not degugged main loop for linear
;      binning... check 30/10/2010
;
;-
; ----------------------------------------------------------

; options for compilation (recommended by RSI)
  COMPILE_OPT idl2, HIDDEN  
  
; watch out for errors
  on_error, 2

; ---------------------------------------------------------
; Check arguments

; is the data array well-defined?

  n = N_ELEMENTS(y)
  if (n lt 2) then BEGIN
      PRINT,'** Not enough data Y in REGROUP.'
      RETURN, -1
  endif

; is the X array empty, if so generate a dummy

  if NOT KEYWORD_SET(x) then x = INDGEN(n)+1

; check X and Y are the same size

  if (N_ELEMENTS(x) ne n) then BEGIN
      PRINT,'** X and Y of different size in REGROUP.'
      RETURN, -1
  endif

; is the minimum number of points per bin supplied?

  if (N_ELEMENTS(minbin) eq 0) then minbin=1

; compute a typical input bin spacing, if not supplied

  if NOT KEYWORD_SET(dx) then dx = FLOAT(ABS(x[1] - x[0]))
  if (N_ELEMENTS(dx) ne 1) then BEGIN
      PRINT,'** DX incorrectly specified in REGROUP.'
      RETURN, -1
  endif

; fill out the DY error array if needed

  if (N_ELEMENTS(dy) eq 1) then dy = MAKE_ARRAY(n, VALUE=dy)

; check that MINBIN is integer type

  type = SIZE(minbin, /TYPE)
  if (type lt 2 or type gt 3) then BEGIN
      PRINT,'** MINBIN is not integer in REGROUP.'
      RETURN, -1
  endif

; if using 'adaptive' regrouping (fixed SNR) set BINW = DX 

  if KEYWORD_SET(snr) then binwidth = dx
 
; determine what type of binning to apply: lin/log?

  if NOT KEYWORD_SET(binwidth) then BEGIN
       PRINT,'** BINWIDTH or SNR must be set in REGROUP.'
       RETURN, -1
  endif

  if (binwidth gt 0) then BEGIN
      qlog = 0
  endif
  if (binwidth lt -1) then BEGIN
      qlog = 1
  endif
  if (binwidth le 0 and binwidth ge -1) then BEGIN
      PRINT,'** BINWIDTH incorrectly specified in REGROUP.'
      RETURN, -1
  endif


; ---------------------------------------------------------
; Main routine

; sort into accending order if requested

  if KEYWORD_SET(sort) then BEGIN
      indx = SORT(x)
      x = x[indx]
      y = y[indx]
      if (N_ELEMENTS(dy) eq n) then dy = dy[indx]
  endif

; re-scale the x axis - differently for LIN and LOG scaling - to make
; the binning in x easy.

  if NOT QLOG then BEGIN
      binw = CEIL(binwidth / dx)
      xx = (x - x[0]) / dx
  endif else BEGIN
      binw = ALOG10(-binwidth)
      xx = ALOG10(x/x[0])
  endelse

; determine the bin locations by defining the lower/upper boundary of each
; bin. 

; First establish arrays to contain the boundaries.
; The BIN array lists the partition points for the binning, i.e. the
; index of the highest X point in each new (output) bin.
; (These arrays are made larger than needed and are trimmed down after
; being populated.)

  bin = MAKE_ARRAY(n, /LONG)
  bin_u = MAKE_ARRAY(n, /FLOAT)
  bin_l = MAKE_ARRAY(n, /FLOAT)

; Define the lower edge of the lowest (output) bin

  j = 0
  bin[j] = 0
  bin_l[j] = x[0] - 0.5*dx

; Loop over all input data points. Assign each input point to a bin,
; j. Each time a bin is 'full' we move to the next bin. A bin is full if
; it contains at least MINBIN data points *and* is separated from the
; next lowest X value by at least one (new) binwidth *or* has a
; signal/noise >= SNR. The SNR case is handled by a seperate loop.

  no_in_bin = 0

  case KEYWORD_SET(snr) of

; grouping done on the basis of 'minbin' and 'binwidth'

      0: begin
          for i = 0, n-1 do BEGIN
              no_in_bin = no_in_bin + 1
              if (no_in_bin lt minbin) then CONTINUE               
              last_bin = xx[bin[j]]
              if (i eq 0) then last_bin = last_bin - binw - 1.0
              if (xx[i] - last_bin lt binw) then CONTINUE      
              no_in_bin = 0
              bin_u[j] = x[i] + 0.5*dx
              bin_l[j+1] = bin_u[j]
              j = j + 1
              bin[j] = i
          endfor
      end

; grouping done on the basis of SNR

      1: begin
          cur_m = 0             ; running (current) mean
          cur_s = 0             ; running variance

          for i = 0, n-1 do BEGIN
              no_in_bin = no_in_bin + 1
              old_m = cur_m
              cur_m = old_m + (y[i] - old_m)/FLOAT(no_in_bin)
              sigma = cur_m * cur_m
              if (N_ELEMENTS(dy) gt 0) then BEGIN
                  cur_s = cur_s + dy[i] * dy[i]
                  sigma = SQRT(cur_s)/FLOAT(no_in_bin)
              endif else BEGIN
                  if (no_in_bin le 1) then CONTINUE
                  cur_s = cur_s + (y[i] - old_m)*(y[i] - cur_m)
                  sigma = SQRT(cur_s/FLOAT(no_in_bin-1))/SQRT(no_in_bin)
              endelse

; check if the SNR criterion has been met 

              if (sigma le 0.0)  then CONTINUE
              if (ABS(cur_m)/sigma lt snr) then CONTINUE
              no_in_bin = 0
              cur_m = 0.0
              cur_s = 0.0
              bin_u[j] = x[i] + 0.5*dx
              bin_l[j+1] = bin_u[j]
              j = j + 1
              bin[j] = i
          endfor
      end

  endcase

; Define the upper edge of the highest (output) bin

  bin_u[j] = x[n-1] + 0.5*dx

; Remove unused elements of the arrays

  bin = bin[0:(j-1)]
  bin_u = bin_u[0:(j-1)]
  bin_l = bin_l[0:(j-1)]

  n_bins = N_ELEMENTS(bin)

; partition the data into bins

  bin_map = VALUE_LOCATE(bin_l, x)

 ; set up arrays for binned data

  bin_x = MAKE_ARRAY(n_bins, /float, value=0.0)
  bin_y = MAKE_ARRAY(n_bins, /float, value=0.0)
  bin_dy = MAKE_ARRAY(n_bins, /float, value=0.0)
  bin_std = MAKE_ARRAY(n_bins, /float, value=0.0)
  bin_n = MAKE_ARRAY(n_bins, /long, value=0)

; within each bin compute the mean, variance and
; combined error (if requested)

  for i = 0, n_bins-1 do BEGIN
      mask = WHERE(bin_map eq i, count)
      if (count eq 0) then CONTINUE
      bin_n[i] = count
      bin_y[i] = MEAN(y[mask])
      bin_std[i] = STDDEV(y[mask])
      if (N_ELEMENTS(dy) eq n) then bin_dy[i] = SQRT(TOTAL(dy[mask]^2))
  endfor

; compute the combined 'error' 

   if (N_ELEMENTS(dy) eq n) then bin_dy = bin_dy / (bin_n > 1)

; set the bin 'error' to the standard error within the bin
; (if requested)

  if KEYWORD_SET(varerr) then bin_dy = bin_std / SQRT(bin_n > 1)

; define the bin 'centre'

  if KEYWORD_SET(geom) then BEGIN
      bin_x = SQRT(bin_u * bin_l)
  endif else BEGIN
      bin_x = (bin_u + bin_l) * 0.5
  endelse

; ---------------------------------------------------------
; Return to user

  RETURN, bin_y

END
