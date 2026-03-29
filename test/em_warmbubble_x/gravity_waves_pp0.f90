program warmbubble
implicit none
   REAL    , PARAMETER :: g = 9.81  ! acceleration due to gravity (m {s}^-2)

   REAL    , PARAMETER :: r_d          = 287.
   REAL    , PARAMETER :: r_v          = 461.6
   REAL    , PARAMETER :: cp           = 7.*r_d/2.
   REAL    , PARAMETER :: cv           = cp-r_d
   REAL    , PARAMETER :: cpv          = 4.*r_v
   REAL    , PARAMETER :: cvv          = cpv-r_v
   REAL    , PARAMETER :: cvpm         = -cv/cp
   REAL    , PARAMETER :: rcv          = r_d/cv
   REAL    , PARAMETER :: rcp          = r_d/cp
   REAL    , PARAMETER :: rovg         = r_d/g
   REAL    , PARAMETER :: c2           = cp * rcv
   REAL    , PARAMETER :: cpovcv       = cp/(cp-r_d)
   REAL    , PARAMETER :: cvovcp       = 1./cpovcv
   REAL    , PARAMETER :: rvovrd       = r_v/r_d

   REAL    , PARAMETER :: reradius     = 1./6370.e03 

   REAL    , PARAMETER :: p1000mb      = 100000.   ,   &
                  th0 = 300.0     ,    &  !300.0
                  p00 = 100000.  ,    &
				     zt = 10000.0,     &
				     dx = 1000.0 ,    &
					 dz = 1000.0  ,   &
					 dth0 = 0.01  ,   &
                       xm = 100000. ,  &
					   xr =   5000. ,  &
					   N2 =   0.0001,   &
					   epssm = 0.1 

integer, parameter :: im =  301  ,   &
                      km =   11 ,    &
					  nz =    201,   &
                      nz2= 301
					  
    REAL    , PARAMETER ::   dzs = 50.    


real,dimension(nz) :: th,p,rho,t
real,dimension(nz) :: zz

real,dimension(im,nz) :: thz,pz,rhoz,tz
real,dimension(im,nz) :: xz,rho0z,p0z,dthz,Th0z

real,dimension(im,km-1) :: ps,rhos,rho0s,p0s,zs,Ths,dth,Th0s
real,dimension(im,km) ::psgm,z8w,p8w,rho8w,th8w,p08w,rho08w,pp8w,rhop8w
real,dimension(im) :: pst,x
real,dimension(km-1) :: sgmh
real,dimension(km) ::sgm

real  :: pi,DLnTh,ptop

integer :: i,k,it



! 20220720 
! ztop,ht,znw,znu,z8w,z,   &
real :: ztop
real,dimension(km) ::znw,znu
real,dimension(im,km) :: z2
real,dimension(im) :: ht



ht = 0.0 



pi = 3.1415926  !2*asin(1.)

 


DLnTh = N2*dzs/g
print*, 'DLnTh = ', DLnTh
do k=1,nz
     TH(k)=log(Th0)+(k-1)*DlnTh
enddo

Th = exp(Th)

! 生产input_sounding文件用于WRF
open(30,file='input_sounding.txt')

write(30,*) '  1000.000     300.00        0.00'
do k=1,nz2

    write(30,'(2x,5f12.2)')   (k-1)*dzs, exp(log(Th0)+(k-1)*DlnTh),0.0, 20.0,0.0 
enddo

close(30) 

stop

do k = 1,nz
zz(k)   = (k-1)*dzs
enddo


do i=1,im
  x(i)=dx*(i-1)
enddo


do k=1,nz-1

  do i=1,im 
    ! if(abs(x(i)-xm) <= 10*xr ) then        ! important for  peoriod bdy
         Thz(i,k) = Th(k) ! +  dth0*sin(pi*zz(k)/zt)/(1.+(x(i)-xm)*(x(i)-xm)/(xr*xr)) 
    !  else
	!     Thz(i,k) = Th(k)  
    !endif
  enddo
enddo

Thz(:,nz)=th(nz)


pz(:,1) = p00
tz(:,1) = thz(:,1)*(pz(:,1)/p00)**(r_d/cp)
rhoz(:,1) = pz(:,1)/r_d/tz(:,1)

do i=1,im
 do k=2,nz
    rhoz(i,k) = rhoz(i,k-1)
    do it=1,20
      pz(i,k) = pz(i,k-1) - 0.5*dzs*(rhoz(i,k)+rhoz(i,k-1))*g
      rhoz(i,k) = 1./((r_d/p1000mb)*th(k)*((pz(i,k)/p1000mb)**cvpm))
    enddo
 enddo

enddo

!print*, pz(:,nz)

ptop=pz(1,nz)

do k = 1,km
p08w(:,k)   = pz(1,(k-1)*20+1)
rho08w(:,k)   = rhoz(1,(k-1)*20+1)
z8w(:,k)   = ((k-1)*20)*dzs
enddo




do k=1,nz-1

  do i=1,im 
    ! if(abs(x(i)-xm) <= 10*xr ) then        ! important for  peoriod bdy
         Thz(i,k) = Th(k) +  dth0*sin(pi*zz(k)/zt)/(1.+(x(i)-xm)*(x(i)-xm)/(xr*xr)) 
    !  else
	!     Thz(i,k) = Th(k)  
    !endif
  enddo
enddo

Thz(:,nz)=th(nz)


do k = 1,km
th8w(:,k)   = thz(:,(k-1)*20+1)
enddo

!print*, maxval(abs(th8w))
!stop



do k = 1,km
p8w(:,k)   =  p08w(:,k) 
!rho8w(:,k)   = rho08w(:,k)*(1.-(Th8w(:,k)-Th8w(1,k))/Th8w(1,k)) 
rho8w(:,k)=(p8w(:,k)**(5./7.))*(100000.**(2./7.))/th8w(:,k)/r_d

enddo


do k = 1,km-1
p0s(:,k)   = (p08w(:,k) + p08w(:,k+1) )/2.
rho0s(:,k)   = (rho08w(:,k) + rho08w(:,k+1) )/2.
enddo


do k = 1,km-1
ps(:,k)   = (p8w(:,k) + p8w(:,k+1) )/2.
rhos(:,k)   = (rho8w(:,k) + rho8w(:,k+1) )/2.
enddo



pst(:) = pz(:,1)-ptop != p00-ptop


do k=1,km
  sgm(k)  =  (p8w(1,k)-p8w(1,km))/pst(1)
enddo
print*, 'sgm= '
print*, sgm





 open(   unit    = 10,                  &
          file    = 'fmodel_input.grd',     &
	      convert = 'big_endian',           &
	      form    = 'unformatted'   )

!write(10) sgm,pst,z8w,rho0s,p0s,rhos,ps

write(10) sgm,pst,z8w,rho0s,p0s,rhos,ps

close(10)
print*, rhos(50,:)-rho0s(50,:)

print*, dth(50,1:km-1)





!   -------------------------------


ztop = z8w(1,km)

do k=1,km
  znw(k)  =  z8w(1,k)
enddo

do k=1,km-1
 znu(k) = 0.5*(znw(k+1)+znw(k))
 
enddo

 ! ht这里是0，不能真有地形，z8w也是无地形的。
do k=1,km-1
       do i=1,im   
	       z2(i,k)  = znu(k) *((ztop - ht(i))/ztop) + ht(i) 
       enddo 
   ENDDO
   
   

open(17, file = 'fmodel_input_gw.grd',     &
	   convert = 'big_endian',      &
	      form = 'unformatted'   )


write(17) ztop,ht,znw,znu,z8w,z2,   &
          rho0s,p0s,    &
		  rhos,ps

close(17)



open(17,file = 'fmodel_input_gw.txt'  )


write(17,*) ztop,ht,znw,znu,z8w,z2,   &
          rho0s,p0s,    &
		  rhos,ps

close(17)


print*, znw


print*, 'stop 9999'    
end



      
      