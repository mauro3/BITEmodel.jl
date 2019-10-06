; read .grid files
; fn -- filename
; All others are output:
; da -- data array
; header -- header as string
pro READ_AGR, fn, da, header=header, xx=xx, yy=yy, ncols=ncols, nrows=nrows, $
              xllcorner=xllcorner, yllcorner=yllcorner, $
              cellsize=cellsize, nodata_value=nodata_value
header=strarr(6) & openr,1, fn & readf,1, header
ncols=long(strmid(header(0),6,40)) & nrows=long(strmid(header(1),6,40))
da=intarr(ncols,nrows) & readf,1, da & close, 1
xllcorner=double(strmid(header(2),10,40))
yllcorner=double(strmid(header(3),10,40))
cellsize=double(strmid(header(4),9,40))
nodata_value=double(strmid(header(5),13,40))
a=cellsize/2d & xx=lindgen(ncols)*cellsize+xllcorner+a
yy=lindgen(nrows)*cellsize+yllcorner+a
da=rotate(da,7)		; lower left corner is da(0,0)

END  ; {READ_AGR}

; .......................

; write .grid files
pro WRITE_AGR, fn, da, header=header, xx=xx, yy=yy, ncols=ncols, nrows=nrows, $
              xllcorner=xllcorner, yllcorner=yllcorner, $
              cellsize=cellsize, nodata_value=nodata_value, format=format
n=size(da) & dar=rotate(da,7)
if KEYWORD_SET(header) then begin
    openw, 4, fn
    for i=0,5 do printf, 4, header(i) & for i=0l,n(2)-1 do printf, 4, dar(*,i)
    close, 4
endif else begin
    if not KEYWORD_SET(xllcorner) then if KEYWORD_SET(xx) then xllcorner=xx(0)-(xx(1)-xx(0))/2
    if not KEYWORD_SET(yllcorner) then if KEYWORD_SET(yy) then yllcorner=yy(0)-(yy(1)-yy(0))/2
    if not KEYWORD_SET(ncols) then ncols=n(1)
    if not KEYWORD_SET(nrows) then nrows=n(2)
    openw, 4, fn
    printf, 4, 'NCOLS       ', strtrim(ncols,2)
    printf, 4, 'NROWS       ', strtrim(nrows,2)
    printf, 4, 'XLLCORNER   ', xllcorner
    printf, 4, 'YLLCORNER   ', yllcorner
    printf, 4, 'CELLSIZE    ', cellsize
    printf, 4, 'NODATA_VALUE', nodata_value
      for i=0l,n(2)-1 do printf, 4, dar(*,i) & close, 4
endelse

END  ; {WRITE_AGR}


; .......................

; read .bin files
pro READ_BINGR,fn,grid,header=header,xx=xx,yy=yy,nc=nc,nr=nr, $
              xll=xll,yll=yll,cs=cs,noval=noval

 openr,u,fn,/get_lun & x=fltarr(12) & readu,u,x
 nc=x(0) &  nr=x(1) &  xll=x(2) & yll=x(3) & cs=x(4) & noval=x(5)
 ;; HEADER-string erstellen
 header=strarr(6)
 header(0)='ncols     '+string(nc   ,fo='(i12)')
 header(1)='nrows     '+string(nr   ,fo='(i12)')
 header(2)='xllcorner '+string(xll  ,fo='(d14.6)')
 header(3)='yllcorner '+string(yll  ,fo='(d14.6)')
 header(4)='cellsize  '+string(cs   ,fo='(d16.9)')
 header(5)='noval     '+string(noval,fo='(d12.2)')

 ;; DATEN einlesen
 grid=fltarr(x(0),x(1)) & readu,u,grid
 free_lun,u
 grid=rotate(grid,7)

 ;; xx und yy erstellen
 a=cs/2d & xx=lindgen(nc)*cs+xll+a
 yy=lindgen(nr)*cs+yll+a & xxb=xx & yyb=yy

end

; .......................
; write .bin files
pro WRITE_BINGR,fn,grid,header=header,xx=xx,yy=yy,nc=nc,nr=nr, $
              xll=xll,yll=yll,cs=cs,noval=noval

n=size(grid)
if n(0) ne 2 then stop, '%% WRITE_BINGR :  GRID has to be a 2-dim array'
x=fltarr(12)
if KEYWORD_SET(header) then begin
  if n_elements(header) ne 6 then stop,'%% WRITE_BINGR :  wrong dimension of HEADER'
  for i=0,5 do begin
    t1=strsplit(header(i),' ')
    n1=n_elements(t1)
    x(i)=t1(n1-1)
  endfor
endif else begin
  ;; Falls "header" nicht vorhanden: angegebene Werte übernehmen
  if not KEYWORD_SET(nc) then x(0)=n(1) else x(0)=nc             ;; nc
  if not KEYWORD_SET(nr) then x(1)=n(2) else x(1)=nr             ;; nr
  if not KEYWORD_SET(xll) then stop, '%% WRITE_AGR :  XLL is not defined' $
     else  x(2)=xll
  if not KEYWORD_SET(yll) then stop, '%% WRITE_AGR :  YLL is not defined' $
    else x(3)=yll
  if not KEYWORD_SET(cs) then $
    stop, '%% WRITE_AGR :  CS is not defined' else x(4)=cs       ;; cs
  if not KEYWORD_SET(noval) then $
    stop, '%% WRITE_AGR :  NOVAL is not defined'else x(5)=noval  ;; noval
endelse
;;--- .bin-File schreiben
OPENW,unit,fn,/Get_Lun
WriteU,unit,x,float(rotate(grid,7))
FREE_LUN,unit

end


;;;;;;;;;;;;;;

;; .agr read-write round trip works:
READ_AGR,'wiki.agr',tt,header=head;,xx=xx,yy=yy,ncols=nc,nrows=nr,xllcorner=xll,yllcorner=yll,cellsize=cells,nodata_value=noval
WRITE_AGR,'wiki2.agr',tt,header=head ;,xx=xx,yy=yy,nc=nc,nr=nr,xll=xll,yll=yll,cs=cells,noval=noval

;; .bin read-write round trip errors:
WRITE_BINGR,'wiki.bin',tt,header=head ;,xx=xx,yy=yy,nc=nc,nr=nr,xll=xll,yll=yll,cs=cells,noval=noval
READ_BINGR,'wiki.bin',tt,header=head ;,xx=xx,yy=yy,nc=nc,nr=nr,xll=xll,yll=yll,cs=cells,noval=noval
;WRITE_AGR,'wiki3.agr',tt,header=head ;,xx=xx,yy=yy,nc=nc,nr=nr,xll=xll,yll=yll,cs=cells,noval=noval

end
