load "$NCARG_NCARG/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_NCARG/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_NCARG/nclscripts/csm/contributed.ncl"
begin

  ; Read in MSSSensm calculated from forecast MSE using Eqs. 3-4 and reference MSE derived from PRISM climatology for 1981-2000
  fid = addfile("wrfcf_fullrand_MSSSensm_HUC2-18.Tmax.nc", "r")
  ss = fid->MSSS_ensm
  delete(fid)

  ; Calculate skill score differences of every ensemble size from that of 200
  dss200 = new((/ 6, 200 /), "float")
  do n = 0, 199
    dss200(:, n) = abs(ss(n, :) - ss(199, :))
  end do

  ; For each lead day, find the minimum ensemble size as the first instance where the skill score drops below a certain fraction of
  ; the skill score for an ensemble size of 1
  do r = 0, 5

    ; For MSSSensm, the computationally inexpensive threshold is 5% of dss200(0).
    thresh = 0.05 * dss200(r, 0)
    iabv = ind(dss200(r, :) .le. thresh)
    print(r + " " + nens(iabv(0)))
    delete(iabv)

    ; For MSSSensm, the stringent threshold is 1% of dss200(0).
    thresh = 0.01 * dss200(r, 0)
    iabv = ind(dss200(r, :) .le. thresh)
    print(r + " " + nens(iabv(0)))
    delete(iabv)

  end do

end

