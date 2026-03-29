'reinit'
'set display color white'
'c'
*'run rgbset.gs'
'set rgb 21 255 240 245'
'set rgb 22 238 210 238'
'set rgb 23 125 38 205'
'set rgb 24 122 103 238'

'set rgb 31 0 0 205'
'set rgb 32 0 0 238'
'set rgb 33 0 0 255'
'set rgb 34 30 144 255'

'set rgb 41 118 238 198'
'set rgb 42 0 205 102'
'set rgb 43 0 238 118'
'set rgb 44 192 255 62'

'set rgb 51 255 255 0'
'set rgb 52 255 235 0'
'set rgb 53 255 215 0'
'set rgb 54 238 180 34'

path='D:\'


'open 'path'fm_grav2d_50m.ctl'
'set parea 1 10 2 6'
'enable print 'path'cb.gmf'

*'set z 1 32'
'set xlopts 1 1 0.16'
'set ylopts 1 1 0.16'
'set mproj off'
*'set x 256 513'

*'set x 129 230'
'set x 513 1024'

'q file'
rec=sublin(result,5)
_endtime=subwrd(rec,12)

runscript = 1
dis_t = 1


*'set z 1 60'
*'set lev 0 6000'

while(runscript)

'set t ' dis_t
'q dims'
say result
rec=sublin(result,5)
_analysis=subwrd(rec,6)

tim = (subwrd(rec,9)-1)*60

say 'Time is ' _analysis

'set grads off'
'set grid off'
if(1) 
say ok
'set gxout shaded'
*'set cint 1'
'set clevs -15 -14 -13 -12 -11 -10  -9  -8  -7  -6  -5  -4  -3  -2 -1'
'set rbcols 21 22 23 24 31 32 33 34 41 42 43 44 51 52 53 54 0'
*'set rbcols 9 14 4 11 5 13 3 10 7 12 8 2 6 1 0'
*'set xlabs 25.6|32.0|38.4|44.8|51.2'

*'set xlabs 0|||||5|||||10|||||15|||||20'


*'set xlabs 0|5|10|15|20'
'set ylabs 0|2|4|6'
*'d smth9(th-300)'

*'d th-th0'
'd pow(pfm,5./7.)*pow(100000.,2./7.)/rhofm/287.-pow(PRESSURE_B(t=1),5./7.)*pow(100000.,2./7.)/RHO_BASE(t=1)/287.'

*'d pow(pfm,5./7.)*pow(100000.,2./7.)/rhofm/287'


'cbarn 0.8 0 5.5 0.3'
endif

'set gxout contour'
'set clevs -13 -12 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 '
*'set cint 1'
'set rbcols 9 14 4 11 5 13 3 10 7 12 8 2 6 1 '
'set grid off'
if(dis_t = 1 )
'set clab off'
endif

*'d th-th0'
'd pow(pfm,5./7.)*pow(100000.,2./7.)/rhofm/287.-pow(PRESSURE_B(t=1),5./7.)*pow(100000.,2./7.)/RHO_BASE(t=1)/287.'

*'d pow(pfm,5./7.)*pow(100000.,2./7.)/rhofm/287'

'set strsiz .2'
'set string 1 l 1'
'set font 0'
'draw title 50m grav2d   t='tim'[sec]'
'set strsiz .2'
'draw ylab  z[km]' 
'draw xlab  x[km]' 
'print'

'printim 'path'cb50m'tim's-f5.jpg'

pull dummy
'c'
if ( dis_t=_endtime )
 runscript=0
endif 
dis_t = dis_t + 1
endwhile
'disable print'
*'reinit'
*'cb-energy.gs'



