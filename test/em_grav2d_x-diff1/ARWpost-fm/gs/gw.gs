'open d:\fm\gw1000m.ctl'

'enable print d:\fm\gw_s_e7_smth9.gmf'

'set z 1 10'
*'set lev 0 10000'

*'set x 901 1201'
'set parea 0.5 10.5 2 6'


it=1
while(it<=11)
'set t 'it

tim=(it-1)*300


'set clab off'
'set grads off'
'set ccolor 1'
*'set cstyle 5'

'set black -0.0001  0.0001'
'd pow(pfm,5./7.)*pow(100000.,2./7.)/rhofm/287.-pow(PRESSURE_B(t=1),5./7.)*pow(100000.,2./7.)/RHO_BASE(t=1)/287.'


say result
'draw title 'tim's '
'print'


'printim D:\ideal-test\gravity_wave\Console1\Console1\gw_'tim's.gif'
pull dummy



'c'
it=it+1
endwhile
'reinit'

