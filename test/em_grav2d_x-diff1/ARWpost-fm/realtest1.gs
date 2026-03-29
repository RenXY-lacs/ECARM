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




 
*test='fm_real_full'
test='fm_real_full'
'open 'test'.ctl'



*'set x 1 61' 
*'set y 1 61' 
*'set x 10 51' 
*'set y 10 51'  


it=25
'set t 'it
 
timf=(it-1)*60
  
'c'

*'set xlab off'
*'set ylab off'
'set gxout shaded'
'set cmin 0.1'
 

*'run  rgbset.gs'
'run d:\rgbset.gs'



'set grads off' 
'set xlopts 1 1 0.16'
'set ylopts 1 1 0.16'
'set gxout shaded'
 
'set grads off'
'set gxout stream'
'set z 10'
'd u;v'
'cbar'
 
'draw title FM uv z=10' 
'printim d:\uv-fm-'timf'mn.jpg'
*pull dummy
'c'
'set grads off'
'set z 1'
'd hgt'
'cbar'
'draw title terrain'
'printim d:\terrain.jpg'
*pull dummy
'c'
'set grads off'
'set z 10'
'd w'
'cbar'
'draw title w z=20'
'printim d:\w-z20-'timf'mn.jpg'

 'c'
 
'set grads off'
'set gxout stream'
'set z 34'
'd v'
'cbar'
'draw title FM v z=34' 
'printim d:\v-fm-z34-'timf'mn.jpg'
*pull dummy
'c'


'set grads off'
'set gxout stream'
'set z 10'
'd u'
'cbar'
 
'draw title FM u z=10' 
'printim d:\u-fm-z10-'timf'mn.jpg'
*pull dummy
'c'

'set grads off'
'set gxout shaded'
'set z 1'
'd pfm'
'cbar'
 
'draw title FM pfm z=1 t='timf'mn' 
'printim d:\p-fm-z1-'timf'mn.jpg'
*pull dummy
'c'



'set lon 160'
'set lat 34'
'set z 1 50'
'set t 1 13'

'd tloop(pfm/rhofm-pfm(t=1)/rhofm(t=1)))/284.07'
'cbar'
'printim d:\tdifference-2hr.jpg'


'c'
'd tloop(pfm-pfm(t=1))'
'draw title pfm-pfm(t=1)' 
'cbar'
'printim d:\pdifference-2hr.jpg'
*'reinit'




