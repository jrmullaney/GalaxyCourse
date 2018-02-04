PRO photz

  spec = mrdfits('spec-0285-51930-0183.fits', 1, head)
  pson, page_s=[21,20], filen='Spec1.ps'
  plot, 10d^spec.loglam, spec.flux, $
        xra=[min(10d^spec.loglam), max(10d^spec.loglam)], /xst, $
        yra=[min(spec.flux), max(spec.flux)], /yst, $
        xtit='Wavelength (A)', ytit='Flux'
  psoff
  psconvert, 'Spec1'

  pson, page_s=[21,20], filen='Spec_NoAx.ps'
  plot, 10d^spec.loglam, spec.flux, $
        xra=[min(10d^spec.loglam), max(10d^spec.loglam)], xst=5, $
        yra=[min(spec.flux), max(spec.flux)], yst=5, $
        xtit='Wavelength (A)', ytit='Flux'
  psoff
  psconvert, 'Spec_NoAx'

  
  ;Bin by 10:
  y10 = regroup(spec.flux, 10d^spec.loglam, binwidth=70, bin_x=x10) 
  pson, page_s=[21,20], filen='Spec10.ps'
  plot, x10, y10, $
        xra=[min(10d^spec.loglam), max(10d^spec.loglam)], /xst, $
        yra=[min(spec.flux), max(spec.flux)], /yst, $
        xtit='Wavelength (A)', ytit='Flux', psym=1
  psoff
  psconvert, 'Spec10'


  ;Bin by 100:
  y10 = regroup(spec.flux, 10d^spec.loglam, binwidth=500, bin_x=x10) 
  pson, page_s=[21,20], filen='Spec100.ps'
  plot, x10, y10, $
        xra=[min(10d^spec.loglam), max(10d^spec.loglam)], /xst, $
        yra=[min(spec.flux), max(spec.flux)], /yst, $
        xtit='Wavelength (A)', ytit='Flux', psym=1
  psoff
  psconvert, 'Spec100'

END
