'reinit'
'set display color white'
'c'
test='test'
'! mkdir pic'test''
'open fm_real_'test'.ctl' 
'q file'
subdim=sublin(result,5)
nt=subwrd(subdim,12)
say nt 'nt' 
it=1
while(it<=nt)
ith=it-1 
'set t 'it

'set z 3'

'set grads off'
'set grid off'
'set cthick 6'
*'set cint 50'
'set gxout contour'
'd pfm'
'draw title pressure 'ith'hr '
'printim  pic'test'/pfmt2-'ith'hr.jpg x1024y768'
'c'


'set grads off'
'set grid off'
'set cthick 6'
'set gxout contour'
'd rhofm'
'draw title density rho 'ith'hr '
'printim  pic'test'/rhofm-'ith'hr.jpg x1024y768'
'c'

'set grads off'
'set grid off'
'set cthick 6'
'set gxout contour'
'd pfm/rhofm/287.04'
'draw title  Tv 'ith'hr '
'printim  pic'test'/tv-'ith'hr.jpg x1024y768'
'c'




'set grads off'
'set grid off'
'set cthick 6'
'set gxout contour'
'd u'
'draw title u 'ith'hr '
'printim  pic'test'/u-'ith'hr.jpg x1024y768'
'c'

'set grads off'
'set grid off'
'set cthick 6'
'set gxout contour'
'd v'
'draw title v 'ith'hr '
'printim  pic'test'/v-'ith'hr.jpg x1024y768'
'c'


'set grads off'
'set grid off'
'set cthick 6'
'set gxout contour'
'set z 2' 
'd w'
'draw title w 'ith'hr '
'printim  pic'test'/w-'ith'hr.jpg x1024y768'
'c'


'set grads off'
'set grid off'
'set cthick 6'
'set gxout contour'
'd qvapor'
'draw title qvapor 'ith'hr '
'printim  pic'test'/qvapor-'ith'hr.jpg x1024y768'
'c'


'set grads off'
'set grid off'
'set cthick 2'
'set gxout stream'
'd u;v'
'draw title u;v 'ith'hr '
'printim  pic'test'/stream-'ith'hr.jpg x1024y768'
'c'

'set z 1' 

'set grads off'
'set grid off'
'set gxout shaded'
'd rainnc+rainc'
'draw title rain 'ith'hr ' 
'cbarn'
'printim pic'test'/rain-'ith'hr.jpg x1024y768'
'c'



'set grads off'
'set grid off'
'set gxout shaded'
'd rainnc'
'draw title rainnc 'ith'hr ' 
'cbarn'
'printim pic'test'/rainnc-'ith'hr.jpg x1024y768'
'c'

'set grads off'
'set grid off'
'set gxout shaded'
'd rainc'
'draw title rainc 'ith'hr ' 
'cbarn'
'printim pic'test'/rainc-'ith'hr.jpg x1024y768'
'c'


it=it+1
endwhile
'quit'

