load "$NCARG_NCARG/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_NCARG/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_NCARG/nclscripts/csm/contributed.ncl"
load "shapefile_utils.ncl"
load "ensRPS1.ncl"
begin

  Tmax_obs = new((/ 121, 6, 10 /), "float") ; (# rand set, # ens, # days, # subregions)
  ; Read in observations into Tmax_obs

  Tmax_climo = new((/ 121, 200, 6, 10 /), "float") ; (# days, # ens, # days, # subregions)
  ; Read in obs. climatology into Tmax_climo

  dimNames = (/ "set", "ens", "time" /)
  dimSizes = (/ 400, 200, 6 /)
  dimUnlim = (/ False, False, False /)

  rps = new((/ 400, 200, 121, 6, 10 /), "float")
  rpr = new((/ 121, 6, 10 /), "float")
  avg_rpss = new((/ 200, 6 /), "float")

  ; Read in calculated CRPSs for each randomization set
  do e = 0, 399
    fid = addfile("wrfcf_fullrand" + e + "_rps_HUC2-18.Tmax.nc", "r")
    rps(e, :, :, :, :) = fid->RPS
    delete(fid)
  end do
  rps!0 = "set"
  rps!1 = "ens"
  rps!2 = "day"
  rps!3 = "lead"
  rps!4 = "basin"

    ; Calculate climatological CRPS
    do t = 0, 5
      do s = 0, 9

        do dy = 1, ndays

         Tmax_ref = Tmax_climo(tb + dy - 1, :, t, s)
         ireal = ind(.not. ismissing(Tmax_ref))
         nreal = dimsizes(ireal)
         if(.not. all(ismissing(ireal))) then
           Tmax_sort = Tmax_ref(ireal)
           qsort(Tmax_sort)            ; Sort climatological values in order
           rpr(tb + dy - 1, t, s) = ensRPS1(Tmax_sort, Tmax_obs(tb + dy - 1, t, s)) ; ensRPS1 calculates CRPS
           delete(Tmax_sort)
         end if
         delete(ireal)
         delete(Tmax_ref)

        end do

      end do
    end do

  end do

  do t = 0, 5 ; over lead day

    rps_t = rps(ens|:, set|:, day|:, lead|t, basin|:)
    avg_rpss(:, t) = 1. - (dim_avg_n(rps_t, (/ 1, 2, 3 /)) / avg(rpr(:, t, :)))
    delete(rps_t)

  end do

  ; Write out CRPSS to file
  systemfunc("rm -f wrfcf_fullrand_RPSS_HUC2-18.Tmax.nc")
  fout = addfile("wrfcf_fullrand_RPSS_HUC2-18.Tmax.nc", "c")
  setfileoption(fout, "DefineMode", True)

  filedimdef(fout, dimNames, dimSizes, dimUnlim)
  filevardef(fout, "RPSS_avg", typeof(avg_rpss), dimNames(1:2))
  setfileoption(fout, "DefineMode", False)

  fout->RPSS_avg = (/ avg_rpss /)
  delete(fout)

end

