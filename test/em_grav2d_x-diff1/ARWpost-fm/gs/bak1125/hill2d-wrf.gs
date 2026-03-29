'reinit'
'set display color white'
'c'

'open d:\hillwrf200m.ctl'
*'open d:\hillwrf200m-1w6nodamp.ctl'
'set t 7'
'set cint 0.3'

'set grads off'
'set xlopts 1 4 0.18'
'set ylopts 1 4 0.18'


'd w'
'draw title WRF w dx=200  2160s ztop=16000 nodamp '
'printim d:\w-200m-ztop16000nodamp.jpg'
