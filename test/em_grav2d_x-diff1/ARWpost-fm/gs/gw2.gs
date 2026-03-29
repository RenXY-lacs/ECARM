'open d:\fm\gw1000m.ctl'
'set display color white'
'c'
'enable print D:\fm\gw_fm_smth9.gmf'
'set z 1 10'
'set lev 0 10000'
'set parea 1 10.5 2 6'

'set xlopts 1 1 0.16'
'set ylopts 1 1 0.16'
'set mproj off'


'set t  1 11'
*'define thp=th-th0'
say result
pull dummy

it=1
while(it<=13)
'set t 'it
*'set xlabs 0|||||100|||||200|||||300'
'set ylabs 0|||||5|||||10'

'set xlopts 1 1 0.16'
'set ylopts 1 1 0.16'
'set mproj off'
'set grads off'
'set grid off'

if(it=1)
'set cint 0.001'
else
'set cint 0.0005'
endif

tim=(it-1)*300
'set clab off'
'set grads off'
'set ccolor 1'
'set cstyle 5'
'set clevs 0'
'd const(pow(pfm,5./7.)*pow(100000.,2./7.)/rhofm/287.-pow(PRESSURE_B(t=1),5./7.)*pow(100000.,2./7.)/RHO_BASE(t=1)/287.,0,-u)'
say result
pull dummy
'set ccolor 1'
'set cstyle 1'

if(it=1)
'set cint 0.001'
else
'set cint 0.0005'
endif


'set cmin 0.0000001'
'd const(pow(pfm,5./7.)*pow(100000.,2./7.)/rhofm/287.-pow(PRESSURE_B(t=1),5./7.)*pow(100000.,2./7.)/RHO_BASE(t=1)/287.,0,-u)'
say result
pull dummy
'set ccolor 1'
'set cstyle 3'

if(it=1)
'set cint 0.001'
else
'set cint 0.0005'
endif

'set cmax -0.00001 '
'd const(pow(pfm,5./7.)*pow(100000.,2./7.)/rhofm/287.-pow(PRESSURE_B(t=1),5./7.)*pow(100000.,2./7.)/RHO_BASE(t=1)/287.,0,-u)'
say result
'set strsiz .2'
'set string 1 l 6'
'draw string 5 6.25    t='tim'[sec]'
'set strsiz .2'
'set string 1 c 3 90'
'draw string 0.5  4. z[km]'
'set string 1 c 3 0'
'draw string 5.5 1.5 x[km]'


'print'
'printim D:\fm\gw_fm_'tim's.gif'
pull dummy

'c'
it=it+1
endwhile
'reinit'

