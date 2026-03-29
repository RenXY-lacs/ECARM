'reinit'
'set display color white'
'c'


'open d:\fm_hill2km_openbdy.ctl'

path='d:\hill-openbdy\'
it=1
while(it<=61)


it2=(it-1)*600


'set t 'it' '


*'set cint 0.07'

'set cint 0.02'



'set mproj off'

'set grads off'
'set xlopts 1 4 0.18'
'set ylopts 1 4 0.18'
'set clab off'

'd w'
*'draw title FM w dx=2000  21600s ztop=16000 nodamp diff2 '

'draw title   W [m/s]     t='it2's  ztop=16000 nodamp diff2  '

*'printim d:\fm-w-2000m-ztop16000nodamp-kh300v75-diff2.jpg'



'printim 'path'hill_'it2's.jpg' 
'print'
pull dummy
'c'
it=it+6
endwhile

'disable print'

*'close 1'
*'reinit'