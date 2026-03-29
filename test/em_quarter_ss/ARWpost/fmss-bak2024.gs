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

 
test='fm_ss-dx500' 
*test='fm_ss_openbdy-1hr'
'open 'test'.ctl'


it=13
'set t 'it
stime= (it-1)*600 



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


 
'c'

'set xlab off'
'set ylab off'
'set gxout shaded'
'set cmin 0.1'
 
*'set clevs  0.1 2 3 5 10 15 20 30 50 70 100 120 150'
*'set ccols  0 9 14 4 11 5 13 3 10 7 12 8 2 6'

'run rgbset.gs'
'set clevs .1 .2 .4 .8 1.6 3.2 6.4 12.8 25.6 51.6 '
'set ccols 0  33 34 35 36 37 38 39 7 8 2'


'set grads off' 
'set xlopts 1 1 0.16'
'set ylopts 1 1 0.16'
'set gxout shaded'
'set cmin 0.01'

'set grads off'
'd rainnc'
say result
'cbarn'

'draw title rainnc FM supercell 'stime's khdif=1000 kvdif=500  ' 
'printim  pic/rainnc-fm-ss-'stime's.jpg'
*pull dummy
'c'
'set gxout stat'
'd rainnc'
say result
  


'set mproj off'

'set xlab on'
'set ylab on'

'set gxout contour'


xvar='w v u pfm rhofm'
xvardes='W V U pressure rho'

iv=1

while(iv<=5) 

var=subwrd(xvar,iv) 
vardes=subwrd(xvardes,iv) 

'set y 1 61' 
iz=1
while(iz<41) 
'set z 'iz 

'set xlopts 1 1 0.16'
'set ylopts 1 1 0.16'



'set grads off'
'set grid off'
 

*'set cint 1'
'd 'var''
'draw title z= 'iz'  'stime's fm   'vardes'   ' 
'printim  pic/fm_'var'_z'iz'_'stime's.jpg'
'c'


iz=iz+8
endwhile 




 

'set z 1 40'

iy=22
while(iy<=35) 
*iy=80
*while(iy<=87) 

'set y 'iy 


'set xlopts 1 1 0.16'
'set ylopts 1 1 0.16'
'set grads off'
'set grid off'

*'set cint 1'
'd 'var' '

'draw title  'var' y='iy' 'stime's fm  '

say result
'printim pic/'var'-y'iy'-'stime's.jpg'
*pull dummy 
iy=iy+3
'c'
endwhile 


iv=iv+1
endwhile
 
'c'
'quit'



