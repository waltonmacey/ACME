;================================================;
;  test11-lat-lon.ncl
;  produces the lat-lon plots for test 11
;================================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"   
; ================================================;

begin

  ; open file and read in data(time,lev,lat,lon) from 0 to n-1
  f    = addfile("./movies/dcmip2016_test21.nc","r")

	lat  = f->lat
	lon  = f->lon
	lev  = f->lev
  time = f->time

	nlat = getfilevardimsizes(f, "lat" )
	nlon = getfilevardimsizes(f, "lon" )
	nlev = getfilevardimsizes(f, "lev" )
  nt   = dimsizes(time)

  ; We want the output at day 6 or 12.
  t_index = nt-1
  if(isvar("t")) then
    t_index = t         ; use t from command line if possible
  end if
  print ("t_index="+t_index)

  ; We require plots of q1, q2, q3 and q4 at 4900 m.
  ; For 60 vertical levels, this corresponds to model level 35

  modheight = nlev-4

  ; Input tracers

	u = f->u (t_index,modheight,:,:)
	v = f->v (t_index,modheight,:,:)
	Q = f->Q (t_index,modheight,:,:)
	T = f->T (t_index,modheight,:,:)
  ps= f->ps(t_index,nlev-1,:,:)
  
  
  ; find lat,lon of minimum pressure
  dims = dimsizes(ps)
  ps_1d = ndtooned(ps)      ; convert 2D array to 1D for use in maxind
  inds = ind_resolve(minind (ps_1d), dims)    ; convert 1D array back to 2D
  ilat = inds(0,0)        ; select the latitude index where the X array is at its' maximum  
  ilon = inds(0,1)        ; select the longitude index where the X array is at its' maximum
  lat_min = ps&lat(ilat)   ; insert the latitude index into the lat coordinate variable
  lon_min = ps&lon(ilon)   ; insert the longitude index into the lon coordinate variable
  print("Minimum pressure located at "+lat_min+", "+lon_min)

  ; Delete tracer units and long names

	;Q@units = " "
	;Q@long_name = " "

  ; Setup the plot
	plot  = new (4, graphic)               			; define plot - need 4 panels

	res1                      = True
	res1@gsnDraw              = True        		; panel plot
	res1@gsnFrame             = True           ; don't draw yet
	res1@cnFillOn             = True
	res1@cnLinesOn            = False
	res1@gsnSpreadColors      = True
	res1@lbLabelAutoStride    = True
	res1@gsnCenterString      = ""
	res1@tiMainString         = "horizontal wind speed   t="+sprintf("%.1g", time(t_index)) +" days"
  res1@gsnMaximize              = True
	res1@vpWidthF             = 0.75
	res1@vpHeightF            = 0.50
  
	res1@cnInfoLabelOn        = False           ; don't give each
	res1@cnLineLabelsOn       = False           ; panel its own
	res1@lbLabelBarOn         = False           ; label bar
	res1@sfXArray             = lon           	; uses lon as plot x-axis
	res1@sfYArray             = lat         		; uses lat for y axis
	res1@trYReverse           = False      			; reverses y-axis, false
	res1@tiYAxisString        = "lat"              ; y-axis title

	res1@tiXAxisString        = "lon"        			; x-axis title
  res1@txFont               = 22 ; helvetica-bold
  res1@tiMainFontHeightF    = 0.02
  res1@tmXBLabelFontHeightF = 0.008
  res1@tiXAxisFontHeightF   = 0.02
  res1@tiYAxisFontHeightF   = 0.02

  res1@gsnPaperOrientation = "portrait"

  ; zoom in on a region
  w = 30.0 ;dg
  res1@trXMinF =  lon_min - w*1.5
  res1@trXMaxF =  lon_min + w*1.5
  res1@trYMinF =  lat_min - w
  res1@trYMaxF =  lat_min + w

  ; set output type and name

  wks_type                  = "pdf"
  wks_type@wkPaperHeightF   = 8
  wks_type@wkPaperWidthF    = 11
  ;wks_type@wkOrientation = "landscape"
	wks  = gsn_open_wks(wks_type,"winds_t"+t_index)

  ; use rainbow colormap with a few of color substitutions

  gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")      ; choose default colormap
  colors       = gsn_retrieve_colormap(wks)   ; retrieve color map for editing
  colors(25 ,:)= (/ 1, 1, 1/)                 ; white     for tracer-free region
  colors(2  ,:)= (/ 0.95, 0.95, 0.95 /)          ; gray for undershoots
  colors(255,:)= (/ 0.4, 0.0, 0.0 /)          ; brick-red for overshoots
  gsn_define_colormap(wks,colors)

  speed = sqrt(u*u+v*v)
  ;copy_VarMeta(u,speed)

	res1@gsnCenterString = ""
	plot = gsn_csm_contour(wks,speed(:,:),res1)		; plot Q1

  ;frame(wks)

end


  
