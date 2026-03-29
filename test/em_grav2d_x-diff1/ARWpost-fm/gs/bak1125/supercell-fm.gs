'reinit'
'set rgb 61 255 230 230'
'set rgb 62 255 200 200'
'set rgb 63 248 160 160'
'set rgb 64 230 140 140'
'set rgb 65 230 112 112'
'set rgb 66 230  80  80'   
'set rgb 67 200  60  60'   
'set rgb 68 180  40  40'
'set rgb 69 164  32  32'

'c'

'open D:\fm\fm_ss_diff2.ctl'
'set mproj off'
'set display color white'
'c'
'set parea  1.8 9 0.8 8'
'set xlopts 1 1 0.16'
'set ylopts 1 1 0.16'



it=1
while(it<=31)
'set t 'it

'set gxout shaded'
'set grads off'

'set grid off'

'd pow(pfm,5./7.)*pow(100000.,2./7.)/rhofm/287.-pow(PRESSURE_B(t=1),5./7.)*pow(100000.,2./7.)/RHO_BASE(t=1)/287.'


'draw title    supercell   'tim's  '

'printim d:\fm\supcell_pic_diff2\supercell_'tim's.gif'
pull dummy
'c'
it=it+1
endwhile
'disable print'
'reinit'




