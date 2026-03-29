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

path='D:\fm\'


'open 'path'fm_ss2.ctl'

*'set parea 1 10 2 6'

'enable print 'path'supercell.gmf'

var='pow(pfm,5./7.)*pow(100000.,2./7.)/rhofm/287.04-pow(Pfm(t=1),5./7.)*pow(100000.,2./7.)/RHO_BASE(t=1)/287.04'
*var='pow(pfm,5./7.)*pow(100000.,2./7.)/rhofm/287.04-pow(PRESSURE_B(t=1),5./7.)*pow(100000.,2./7.)/RHO_BASE(t=1)/287.04'
*var='pow(pfm,5./7.)*pow(100000.,2./7.)/rhofm/287.04'
'set xlopts 1 1 0.16'
'set ylopts 1 1 0.16'
'set mproj off'

'set z 1'



'q file'
rec=sublin(result,5)
_endtime=subwrd(rec,12)

runscript = 1
dis_t = 1
while(runscript)
'set t ' dis_t


'q dims'
say result
rec=sublin(result,5)
_analysis=subwrd(rec,6)

tim = (subwrd(rec,9)-1)*5

say 'Time is ' _analysis

'set xlabs 0|||||20|||||40|||||60|||||80km'
'set ylabs 0|||||20|||||40|||||60|||||80km'

'set gxout shaded'
'set grads off'

'set grid off'

'set cint 0.1'
*'d 'var
*'cbarn 0.8 0 5.5 0.3'

'set gxout contour'
'set cthick 6'
'set cint 0.05'
'd 'var
say result

'set ccolor 1'
'set cthick 2'
'set gxout stream'
'd u-u(t=1);v-v(t=1)'

'set font 0'
'draw title    supercell   'tim'm  '

'printim 'path'supercell_'tim's.gif x1000 y800' 


pull dummy
'c'
if ( dis_t=_endtime )
 runscript=0
endif 
dis_t = dis_t + 1
endwhile
'disable print'
*'reinit'




