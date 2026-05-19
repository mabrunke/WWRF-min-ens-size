
; *******************************************************************
; W. Scheftic
; Computes the continuous ranked probability score (CRPS)
; for an ensemble of forecasts.
; returns scalar rps_sum of the same type as the input forecast.
; Input:
; fcstens[member] - ordered [lowest to highest] forecast members
;                 -   1d array (float or double).
; ob              - observed value (float or double).
function ensRPS1(fcstens[*]:snumeric, ob:snumeric)
local nens, rps_sum, fcstcur, Fo_y, F_y, fcstnxt, Dy1, Dy2, rps, n
begin
  nens = num(.not.ismissing(fcstens))  ; ensemble size
  fcst_typ = typeof(fcstens)           ; type of input array
  rps_sum = totype(0., fcst_typ)       ; sum of all rps values to be returned
  if(nens .lt. 1) then
    rps_sum = fcstens@_FillValue
  else
    fcstnxt = fcstens(0)               ; value of next forecast
    if(ob .lt. fcstnxt) then
    ; Case for observation being less than all members
      Dy1 = fcstnxt - ob               ; width b/w consecutive forecasts (or b/w forecast and ob)
                                       ;   in units of target variable
      rps = Dy1                        ; current rps value
      rps_sum = rps_sum + rps
    end if

    do n = 1, nens - 1
      fcstcur = fcstnxt                ; value of next forecast
      fcstnxt = fcstens(n)
      F_y = n / totype(nens, fcst_typ) ; percentile value of next forecast
      if(fcstnxt .lt. ob) then
      ; Case for current and next forecast being less than observation
	    Fo_y = totype(0., fcst_typ)    ; percentile value of current forecast
	    Dy1 = fcstnxt - fcstcur
        rps = (F_y - Fo_y) * (F_y - Fo_y) * Dy1
      else
        if(fcstcur .ge. ob) then
        ; Case for current and next forecast being greater than observation
	     Fo_y = 1.
	     Dy1 = fcstnxt - fcstcur
	     s = (F_y - Fo_y) * (F_y - Fo_y) * Dy1
        else
        ; Case for observation being between current and next forecast
	     Dy1 = ob - fcstcur
	     Dy2 = fcstnxt - ob            ; 2nd width b/w next forecast and ob for split case
	     rps = F_y * F_y * Dy1 + (F_y - 1.) * (F_y - 1.) * Dy2
	    end if
      end if

      rps_sum = rps_sum + rps
    end do

    if(fcstnxt .lt. ob) then
    ; Case for observation being greater than all members
      Dy1 = ob - fcstnxt
      rps = Dy1
      rps_sum = rps_sum + rps
    end if
  end if

  return(rps_sum)
end
