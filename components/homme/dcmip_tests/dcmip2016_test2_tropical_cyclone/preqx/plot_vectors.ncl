; These files are loaded by default in NCL V6.2.0 and newer
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


begin

  f    = addfile("./movies/dcmip2016_test21.nc","r")
  time = f->time
  nt   = dimsizes(time)
	nlev = getfilevardimsizes(f, "lev" )

  t_index = nt-1
  if(isvar("t")) then
    t_index = t         ; use t from command line if possible
  end if
  print ("t_index="+t_index)

  modheight = nlev-4
	u = f->u (t_index,modheight,:,:)
	v = f->v (t_index,modheight,:,:)
  speed = sqrt(u^2+v^2)
  copy_VarCoords(u,speed)                         ; copy coord vars to speed

  wks = gsn_open_wks("png","panel")               ; send graphics to PNG file
  plot = new(2,graphic)                           ; create a plot array

  res                     = True
  res@gsnDraw             = False           ; don't draw
  res@gsnFrame            = False           ; don't advance frame
  res@cnFillOn            = True            ; turn on color
  res@cnFillPalette       = "gui_default"   ; set color map
  res@cnLinesOn           = False           ; no contour lines
  res@gsnLeftString       = "Speed"         ; change left string
  res@gsnRightString      = u@units         ; assign right string
  res@mpFillOn            = False           ; no map fill

  vecres                  = True            ; vector only resources
  vecres@gsnDraw          = False           ; don't draw
  vecres@gsnFrame         = False           ; don't advance frame
  vecres@vcGlyphStyle     = "CurlyVector"   ; curly vectors
  vecres@vcRefMagnitudeF  = 20              ; define vector ref mag
  vecres@vcRefLengthF     = 0.045           ; define length of vec ref
  vecres@gsnRightString   = " "             ; turn off right string
  vecres@gsnLeftString    = " "             ; turn off left string
  vecres@tiXAxisString    = " "             ; turn off axis label 
  vecres@vcRefAnnoOrthogonalPosF = -.535    ; move ref vector into plot

  plotA   = gsn_csm_contour_map_ce(wks,speed,res)
  plotB   = gsn_csm_vector(wks,u,v,vecres)
  overlay(plotA,plotB)                      ; result will be plotA
  plot(0) = plotA                           ; now assign plotA to array

  res@gsnLeftString       = "Wind"          ; change left string
  plotC   = gsn_csm_contour_map_ce(wks,u,res)
  plotD   = gsn_csm_vector(wks,u,v,vecres)
  overlay(plotC,plotD)                      ; result is plotC
  plot(1) = plotC                           ; now assign plotC to array

  resP                = True                ; panel only resources
  resP@gsnMaximize    = True                ; maximize plots
  gsn_panel(wks,plot,(/2,1/),resP)          ; now draw as one plot

end
