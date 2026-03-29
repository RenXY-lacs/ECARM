'reinit'
'set display color white'
'c'
test='ter'
'! mkdir pic'test''
'open fm_real_'test'.ctl' 
'q file'
subdim=sublin(result,5)
nt=subwrd(subdim,12)
say nt 'nt' 

varx='u v w pfm rhofm qvapor u;v pfm/rhofm/287.04'
varxc='u v w pfm rhofm qvapor stream TV'


iv=1
while(iv<=8)
var=subwrd(varx,iv) 
varc=subwrd(varxc,iv) 

'! mkdir pic'test'/'varc''

say 'varc' varc 

iz=1
while(iz<=5)
'set z 'iz 

it=1
while(it<=nt)
ith=it-1
ith2=ith*5 
'set t 'it




'set grads off'
'set grid off'
'set cthick 6'
*'set cint 50'
'set gxout contour'
'd 'var 
'draw title 'varc'  'iz'   'ith2'mn '
'printim  pic'test'/'varc'_z'iz'-'ith2'mn.jpg x1024y768'
'c'




it=it+1
endwhile

iz=iz+1
endwhile


iv=iv+1
endwhile 


'set z 1' 



varrx='rainnc rainc rainnc+rainc'
varrxc='rainnc rainc rain'

iv=1
while(iv<=3)
var=subwrd(varrx,iv) 
varc=subwrd(varrxc,iv) 

'! mkdir pic'test'/'varc''
say 'varc' varc  

it=1
while(it<=nt)
ith=it-1
ith2=ith*5 
'set t 'it


'set grads off'
'set grid off'
'set gxout shaded'
'd 'var' ' 
'draw title 'varc'  'ith2'mn ' 
'cbarn'
'printim pic'test'/'varc'-'ith2'mn.jpg x1024y768'
'c'


 

it=it+1
endwhile


iv=iv+1
endwhile 

'quit'

