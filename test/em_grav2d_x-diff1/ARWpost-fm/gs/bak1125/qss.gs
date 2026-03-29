path='D:\fm\'
'open 'path'fm_ss_diff2.ctl'

** script needs height level input data
** will look for .75 , 1.5. 4 and 8 km data to plot

'set display color white'
'set xlab off'
'set ylab off'
'run rgbset.gs'

var='pow(pfm,5./7.)*pow(100000.,2./7.)/rhofm/287.04-pow(Pfm(t=1),5./7.)*pow(100000.,2./7.)/RHO_BASE(t=1)/287.04'
*var='pow(pfm,5./7.)*pow(100000.,2./7.)/rhofm/287.04-pow(PRESSURE_B(t=1),5./7.)*pow(100000.,2./7.)/RHO_BASE(t=1)/287.04'
*var='pow(pfm,5./7.)*pow(100000.,2./7.)/rhofm/287.04'


say 'Create gif images as well (1=yes ; 0=no)'
pull ans
frame = 1

'q file'
rec=sublin(result,5)
_endtime=subwrd(rec,12)
_endlev=subwrd(rec,9)

runscript = 1
dis_t = 2

while(runscript)
dis_z = 1
'set t ' dis_t
'q dims'
rec=sublin(result,5)
_analysis=subwrd(rec,6)
say 'Time is ' _analysis

runlevel = 1
while(runlevel)
'set z ' dis_z
'q dims'
rec=sublin(result,4)
newz=subwrd(rec,6)
say 'Level is ' newz

found_level = 0
** Looking for right level to plot
if (newz=2)
 found_level=1
endif
if (newz=4)
 found_level=1
endif
if (newz=9)
 found_level=1
endif
if (newz=17)
 found_level=1
endif

newzkm= (newz*500-250)/1000 

if (found_level)
say 'Plotting'
'c'
'set grads off'
'set gxout shaded'
'd 'var' '
'set gxout contour'
'd p'
'set strsiz .2'
'set string 1 l 6'
'draw string 2.7 8.35   pressure (mb,lines), Theta (color)'
'set strsiz .15'
'set string 1 l 3'
'draw string 6.5 8.0 ' _analysis


'draw string 2.5 8.0 z=' newzkm ' km'
'run cbar.gs'
if(ans)
*'printim d:\fm\supcell_pic_diff2-test\qss'frame'-z='newzkm'km.gif gif '
frame=frame+1
endif

pull dummy
'c'
'set grads off'
'set gxout shaded'
'd w'
'set gxout stream'
*'d u-u(t=1);v-v(t=1)'\
'd u;v'
'set strsiz .2'
'set string 1 l 6'
'draw string 4.7 8.35   W (color)'
'set strsiz .15'
'set string 1 l 3'
'draw string 6.5 8.0 ' _analysis
'draw string 2.5 8.0 z=' newz ' km'
'run cbar.gs'
if(ans)
'printim d:\fm\supcell_pic_diff2-test\qss'frame'-z='newzkm'km.gif gif '
frame=frame+1
endif
pull dummy
'c'




endif

if ( dis_z=_endlev )
 runlevel=0
endif

dis_z = dis_z + 1
endwhile


dis_t = dis_t + 2
if ( dis_t>_endtime )
 runscript=0
endif 

endwhile

