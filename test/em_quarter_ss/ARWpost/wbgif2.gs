'reinit'
'set rgb 61 255 230 230'
'set rgb 62 255 200 200'
'set rgb 63 248 160 160'
'set rgb 64 230 140 140'
'set rgb 65 230 112 112'
'set rgb 66 230  80  80'   
'set rgb 67 200  60  60'   
'set rgb 68 180  40  40'
'set rgb 69 164  32  32'

'c'

'open D:\fm\wb100m.ctl'
'set mproj off'
'set display color white'
'c'
'set parea  1.8 9 0.8 8'
'set xlopts 1 1 0.16'
'set ylopts 1 1 0.16'

*'set x 1 202'
'set z 1 135'

it=1
while(it<=16)
'set t 'it
tim=(it-1)*60
'set gxout shaded'
'set grads off'

*'set rbcols 5 0 61 62 63 64 65 66 67 68 69'
*'set clevs 326.9 327 327.1 327.2 327.3 327.4 327.5 327.6'
'set rbcols 0 61 62 63 64 65 66 67 68 69'
*'set clevs 327.1 327.2 327.3 327.4 327.5 327.6'
'set grid off'
'set xlabs 0||200||400||600||800||1000'
*'set ylabs 0||200||400||600||800||1000'
 
*'d th'
*'d pow(pfm,5./7.)*pow(100000.,2./7.)/rhofm/287.-pow(PRESSURE_B(t=1),5./7.)*pow(100000.,2./7.)/RHO_BASE(t=1)/287.'
'set cmin 300.01'
'd pow(pfm,5./7.)*pow(100000.,2./7.)/rhofm/287'

'cbarn 1.2 1 9.3 4.25'
'set gxout contour'
'set grads off'
'set clab off'
*'set rbcols 5 0 61 62 63 64 65 66 67 68 69'
*'set clevs 326.9 327 327.1 327.2 327.3 327.4 327.5 327.6'
'set rbcols 0 61 62 63 64 65 66 67 68 69'
*'set clevs 327.1 327.2 327.3 327.4 327.5 327.6'
'set grid off'
*'d pow(pfm,5./7.)*pow(100000.,2./7.)/rhofm/287.-pow(PRESSURE_B(t=1),5./7.)*pow(100000.,2./7.)/RHO_BASE(t=1)/287.'
'set cmin 300.01'
'd pow(pfm,5./7.)*pow(100000.,2./7.)/rhofm/287'

'set gxout vector'
*'set arrscl 1 2'
'set ccolor 1'
'set grid off'
'd skip(u,5,5);skip(w,5,5)'

'set strsiz 0.25'
'draw title    Fmodel warm bubble  dx=dz=100m   'tim's  '
'set string 1 c 3 90'
'draw string 1. 4.25 z[m]'
'set string 1 c 3 0'
'draw string 5.5 0.3  x[m]'
*'draw xlab x[m]'
'printim d:\fm\wb_'tim's.gif'
pull dummy
'c'
it=it+1
endwhile
'disable print'
'reinit'




