load "$NCARG_NCARG/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_NCARG/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_NCARG/nclscripts/csm/contributed.ncl"
begin

  ; Read in MSSSunc calculated from the spread in forecast MSE using Eqs. 5-7 and reference MSE derived from PRISM climatology for
  ; 1981-2000
  fid = addfile("wrfcf_fullrand_MSSSunc_HUC2-18.Tmax.nc", "r")
  ss = fid->MSSS_unc
  delete(fid)

  ; Calculate skill score spread as the difference between the 90th and 10th percentiles
  ss_sprd = new((/ 6, 200 /), "float")
  do n = 0, 199
    do t = 0, 5
      stat = stat_dispersion(ss(:, n, t), False)
      hiss = stat(13)
      loss = stat(3)
      ss_sprd(t, n) = abs(hiss - loss)
    end do
  end do

  ; For each lead day, find the minimum ensemble size as the first instance where the skill score drops below a certain fraction of
  ; the skill score for an ensemble size of 1
  do r = 0, 5

    ; For MSSSunc, the computationally inexpensive threshold is 15% of ss_sprd(0).
    thresh = 0.15 * ss_sprd(r, 0)
    iabv = ind(ss_sprd(r, :) .le. thresh)
    print(r + " " + nens(iabv(0)))
    delete(iabv)

    ; For MSSSunc, the stringent threshold is 1% of ss_sprd(0).
    thresh = 0.10 * ss_sprd(r, 0)
    iabv = ind(ss_sprd(r, :) .le. thresh)
    print(r + " " + nens(iabv(0)))
    delete(iabv)

    ; More stringent thresholds produce minimum ensemble size estimates > 100 which are not accurate.
    thresh = 0.05 * ss_sprd(r, 0)
    iabv = ind(ss_sprd(r, :) .le. thresh)
    print(r + " " + nens(iabv(0)))
    delete(iabv)

    thresh = 0.01 * ss_sprd(r, 0)
    iabv = ind(ss_sprd(r, :) .le. thresh)
    print(r + " " + nens(iabv(0)))
    delete(iabv)

   end do

end

