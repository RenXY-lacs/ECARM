'reinit'
'set display color white'
'c'
'open fm_real_full.ctl' 
it=1
while(it<=25)
ith=it-1 
'set t 'it

'set z 1'

'set grads off'
'set grid off'
'set cthick 6'
*'set cint 50'
'set gxout contour'
'd pfm'
'draw title pressure 'ith'hr '
'printim  pic/pfmt2-'ith'hr.jpg x1024y768'
'c'


'set grads off'
'set grid off'
'set cthick 6'
'set gxout contour'
'd rhofm'
'draw title density rho 'ith'hr '
'printim  pic/rhofm-'ith'hr.jpg x1024y768'
'c'

'set grads off'
'set grid off'
'set cthick 6'
'set gxout contour'
'd u'
'draw title u 'ith'hr '
'printim  pic/u-'ith'hr.jpg x1024y768'
'c'

'set grads off'
'set grid off'
'set cthick 6'
'set gxout contour'
'd v'
'draw title v 'ith'hr '
'printim  pic/v-'ith'hr.jpg x1024y768'
'c'


'set grads off'
'set grid off'
'set cthick 6'
'set gxout contour'
'set z 2' 
'd w'
'draw title w 'ith'hr '
'printim  pic/w-'ith'hr.jpg x1024y768'
'c'

'set z 1' 
'set grads off'
'set grid off'
'set cthick 6'
'set gxout contour'
'd qvapor'
'draw title qvapor 'ith'hr '
'printim  pic/qvapor-'ith'hr.jpg x1024y768'
'c'


'set grads off'
'set grid off'
'set cthick 2'
'set gxout stream'
'd u;v'
'draw title u;v 'ith'hr '
'printim  pic/stream-'ith'hr.jpg x1024y768'
'c'


'set grads off'
'set grid off'
'set gxout shaded'
'd rainnc+rainc'
'draw title rain 'ith'hr ' 
'cbarn'
'printim pic/rain-'ith'hr.jpg x1024y768'
'c'







it=it+1
endwhile
'quit'

