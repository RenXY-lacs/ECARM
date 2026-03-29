'reinit'
'set display color white'
'c'

*hx='37   39   41   43   45   47   49   51   53   55   57   59   61   63   65   67   69   71   73   75   77   79   81  83   85   87   89   91   93   95   97   99  101  103  105  107'
*hy='   10.810811   12.360939   14.265335   16.638935   19.646366   23.529411   28.653296   35.587189   45.248867   59.171597   80.000000  112.359543  163.934433  243.902451  344.827576  400.000000  344.827576  243.902451  163.934433  112.359543   80.000000   59.171597   45.248867   35.587189   28.653296   23.529411   19.646366   16.638935   14.265335   12.360939   10.810811    9.532888    8.467401    7.570023    6.807352    6.153846'
hx='37   39   41   43   45   47   49   51   53   55   57   59   61   63   65   67   69   71   73   75   77   79   81  83   85   87   89   91   93   95   97   99  101  103  105  107'
hy='   10.810811   12.360939   14.265335   16.638935   19.646366   23.529411   28.653296   35.587189   45.248867   59.171597   80.000000  112.359543  163.934433  243.902451  344.827576  400.000000  344.827576  243.902451  163.934433  112.359543   80.000000   59.171597   45.248867   35.587189   28.653296   23.529411   19.646366   16.638935   14.265335   12.360939   10.810811    9.532888    8.467401    7.570023    6.807352    6.153846'

dx=2km

var='w'
path='D:\fm\'
'enable print 'path'hill_'var'_'dx'.gmf'

'open 'path'fm_hill_2km.ctl'

*'set lev 0 10000'
*'set x 17 116'
*'set z 1 5'
*'set x 90 110'

it=1
while(it<=4)


if(dx=20m)
it2=(it-1)*60
'set cint 6'
'set xlabs 0|1|2'
endif

if(dx=200m)
it2=(it-1)*60
'set cint 0.3'
*'set xlabs 0|10|20'
'set xlabs -20|0|20'
endif

if(dx=2km)
it2=(it-1)*600
*'set cint 0.07'
'set xlabs -200|-100|0|100|200'
endif

if(dx=20km)
it2=(it-1)*6000
'set cint 0.004'
'set xlabs -2000|-1000|0|1000|2000'
endif

if(dx=200km)
it2=(it-1)*60000
'set cint 0.0003'
'set xlabs  -20000|-10000|0|10000|20000'
endif

'set t 'it' '

*'set lev 0 10000'
*'set z 1 5'
*'set lev 0 15000'


'set mproj off'
'set grads off'
'set xlopts 1 4 0.18'
'set ylopts 1 4 0.18'
*'set ylabs 0 |||||5|||||10'
'set ylabs 0 |||||8|||||16'
'set ccolor 1'
'set grid off'
*'set cmin 0.1'
'd w'
say result

'set font 0'
'set annot  1 3'
'draw title   W [m/s]     t='it2's  '
'set strsiz 0.5'
'draw xlab X [km]' 
'set font 0'
'draw ylab Z [km] '


*  画山脉廓线
if(0=0)
pnt=1

while(pnt<=31)




*wx1=subwrd(hx,pnt) 
*wx2=subwrd(hx,pnt+1)

*for 199 point
wx1=subwrd(hx,pnt) +32
wx2=subwrd(hx,pnt+1) +32

wy1=subwrd(hy,pnt)
wy2=subwrd(hy,pnt+1)




'q w2xy 'wx1' 'wy1' '
x1=subwrd(result,3)
y1=subwrd(result,6)

'q w2xy 'wx2' 'wy2' '
x2=subwrd(result,3)
y2=subwrd(result,6)
'draw line 'x1' 'y1' 'x2' 'y2' '
pnt=pnt+1
endwhile

endif


'printim 'path'hill_'it2's.jpg' 
'print'
pull dummy
'c'
it=it+1
endwhile

'disable print'

'close 1'
'reinit'


