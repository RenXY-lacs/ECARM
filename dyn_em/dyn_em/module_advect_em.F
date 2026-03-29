
!WRF:MODEL_LAYER:DYNAMICS
!
#if ( defined(ADVECT_KERNEL) )
! cpp -traditional -C -P -DADVECT_KERNEL module_advect_em.F > advection_kernel.f90
! gfortran -ffree-form -ffree-line-length-none advection_kernel.f90
! ./a.out
MODULE advection_kernel
   TYPE grid_config_rec_type
      INTEGER :: scalar_adv_opt = 0
      INTEGER :: h_sca_adv_order = 5
      INTEGER :: v_sca_adv_order = 3
      LOGICAL :: periodic_x = .false.
      LOGICAL :: periodic_y = .false.
      LOGICAL :: symmetric_xs = .false.
      LOGICAL :: symmetric_xe = .false.
      LOGICAL :: symmetric_ys = .false.
      LOGICAL :: symmetric_ye = .false.
      LOGICAL :: open_xs = .false.
      LOGICAL :: open_xe = .false.
      LOGICAL :: open_ys = .false.
      LOGICAL :: open_ye = .false.
      LOGICAL :: specified = .true.
      LOGICAL :: nested = .false.
      LOGICAL :: polar = .false.
   END TYPE grid_config_rec_type
   CHARACTER (LEN=256) :: wrf_err_message
CONTAINS
!----------------------------------------------------------------
SUBROUTINE wrf_error_fatal ( message )
   IMPLICIT NONE
   CHARACTER(LEN=*) , INTENT(IN) :: message
   PRINT *,'advect_scalar_pd: FATAL MESSAGE = ',TRIM(message)
   STOP 12345
END SUBROUTINE wrf_error_fatal
!----------------------------------------------------------------
SUBROUTINE init ( config_flags )
   IMPLICIT NONE
   TYPE (grid_config_rec_type) :: config_flags
   config_flags%h_sca_adv_order = 5
   config_flags%v_sca_adv_order = 3
   config_flags%periodic_x = .true.
   config_flags%periodic_y = .true.
   config_flags%symmetric_xs = .false.
   config_flags%symmetric_xe = .false.
   config_flags%symmetric_ys = .false.
   config_flags%symmetric_ye = .false.
   config_flags%open_xs = .false.
   config_flags%open_xe = .false.
   config_flags%open_ys = .false.
   config_flags%open_ye = .false.
   config_flags%specified = .false.
   config_flags%nested = .false.
END SUBROUTINE init
!----------------------------------------------------------------
SUBROUTINE tophat ( field, num_scalars , &
                                ids, ide, jds, jde, kds, kde, &
                                ims, ime, jms, jme, kms, kme, &
                                its, ite, jts, jte, kts, kte )
   IMPLICIT NONE
   INTEGER , INTENT(IN ) :: num_scalars ,     ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte
   REAL , DIMENSION( ims:ime , kms:kme , jms:jme , num_scalars ) , INTENT(OUT) :: field
   INTEGER :: i, j, k , n

   field = 0

   DO n = 1 , num_scalars
   DO j = jts , jte
   DO k = kts , kte
   DO i = its , ite
      IF ( i .gt. 35 .and. i.lt. 55 ) THEN
         field (i,k,j,n) = 1.
      END IF
   END DO
   END DO
   END DO
   END DO
END SUBROUTINE tophat
!----------------------------------------------------------------
SUBROUTINE column (loop , data_list, its,ite)
   IMPLICIT NONE
   INTEGER , INTENT(IN) :: loop, its, ite
   REAL , INTENT(IN) , DIMENSION(its:ite) :: data_list
   INTEGER , DIMENSION(its:ite) :: data_int
   INTEGER :: i
   CHARACTER (len = 10 ) :: filename

   IF ( loop.EQ.0 ) THEN
      OPEN (unit=7,file = "x_locations.txt" , &
      form = "formatted" , &
      access = "sequential" )
      
      DO i = its,ite
         write (7,*) i
      END DO
      close (7)
   END IF

   WRITE(filename,fmt='(i6.6,".txt")') loop
   OPEN (unit=7,file = filename , &
   form = "formatted" , &
   access = "sequential" )
   
   data_int = NINT(data_list * 100 )
   DO i = its,ite
      write (7,*) data_int(i)
   END DO
   close (7)
   
END SUBROUTINE column
!----------------------------------------------------------------
#elif ( ! defined(ADVECT_KERNEL) )

MODULE module_advect_em

  USE module_bc
  USE module_model_constants
  USE module_wrf_error

CONTAINS

!-------------------------------------------------------------------------------

     
!!+==============================================================================

!  计算动能engk

      SUBROUTINE calcu_engk(engk,u,v,w,config_flags,                       &
                                  ids, ide, jds, jde, kds, kde,  &
                                  ims, ime, jms, jme, kms, kme,  &
                                  its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE

      INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
      INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
      INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte

      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( OUT )   ::  engk
      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( IN )    ::  u
      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( IN )    ::  v
      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( IN )    ::  w
	  TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      INTEGER      :: i,j,k,k_start,k_end,i_start,i_end,j_start,j_end
	 
	  REAL :: ub,vb 

    LOGICAL :: specified
	
	 integer :: ktf 
	 
	  ktf=MIN(kte,kde-1)

    specified = .false.

    if(config_flags%specified .or. config_flags%nested) specified = .true.


    

      k_start = kts
      k_end   = MIN(kte,kde-1)
	  
	  
  !  for mpi 扩1      
    i_start = MAX(its-1,ids)
    i_end   = MIN(ite+1,ide-1)
    j_start = MAX(jts-1,jds)
    j_end   = MIN(jte+1,jde-1) 
	  

	   
	   
	   
	!  IF ( config_flags%open_xs ) i_start = MAX(ids,its-1)
    !  IF ( config_flags%open_xe ) i_end   = MIN(ide-1,ite+1)

    !  IF ( config_flags%periodic_x ) i_start = its-1
    !  IF ( config_flags%periodic_x ) i_end = ite+1 
	  
	  
    !  IF ( specified ) i_start = MAX(ids+2,its-1)
    !  IF ( specified ) i_end   = MIN(ide-2,ite+1)	  
    !  IF ( specified ) j_start = MAX(jds+2,jts-1)
    !  IF ( specified ) j_end   = MIN(jte+1,jde-2)  
	  
       engk  = 0.0 
	   
	   
      DO      j = j_start, j_end
      DO  k = k_start, k_end
      DO  i = i_start, i_end
          !  engk(i,k,j) = ( u(i,k,j)*u(i,k,j) + u(i+1,k,j)*u(i+1,k,j) + v(i,k,j)*v(i,k,j) + v(i,k,j+1)*v(i,k,j+1) + w(i,k,j)*w(i,k,j) + w(i,k+1,j)*w(i,k+1,j) ) * 0.5 
			  engk(i,k,j) = ( u(i,k,j)**2 + u(i+1,k,j)**2 + v(i,k,j)**2 + v(i,k,j+1)**2 + w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25 
      ENDDO
      ENDDO
      ENDDO
	  

if(.false.) then ! 20240706	  
! 边角处理未扩1 ，暂时	  

	  
      IF ( specified .AND. its .LE. ids+1 .AND. .NOT. config_flags%periodic_x ) THEN
	  j_start = MAX(jds+2,jts)
	  j_end   = MIN(jte,jde-2)
        DO j = j_start, j_end
        DO k=kts,ktf
           i = ids 
         
           ub = u(i ,k,j)
           IF (u(i+1,k,j) .LT. 0.) ub = u(i+1,k,j)
            engk(i,k,j) = ( ub**2 + u(i+1,k,j)**2 + v(i,k,j)**2 + v(i,k,j+1)**2 + w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25 
        ENDDO
        ENDDO
		
		
		
		
		
      ENDIF	  
	  
      IF ( specified .AND. ite .GE. ide-1 .AND. .NOT. config_flags%periodic_x ) THEN
	  j_start = MAX(jds+2,jts)
	  j_end   = MIN(jte,jde-2)
	  
        DO j = j_start, j_end
        DO k=kts,ktf
           i =  ide-1  
         
           ub = u(i+1 ,k,j)
          IF (u(i,k,j) .GT. 0.) ub = u(i,k,j) 
            engk(i,k,j) = ( ub**2 + u(i,k,j)**2 + v(i,k,j)**2 + v(i,k,j+1)**2 + w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25 
        ENDDO
        ENDDO
		
      ENDIF		  

     IF ( specified .AND. jts .LE. jds+1 ) THEN
     i_start = MAX(ids+2,its)
     i_end   = MIN(ide-2,ite)	 
	 
        j = jds 
        DO k=kts,ktf
        DO i = i_start, i_end
 
           vb = v(i,k,j)
           IF (v(i,k,j+1) .LT. 0.) vb = v(i,k,j+1)

               engk(i,k,j) = ( u(i,k,j)**2 + u(i+1,k,j)**2 + vb**2 + v(i,k,j+1)**2 + w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25 

        ENDDO
        ENDDO
      ENDIF

      IF ( specified .AND. jte .GE. jde-1 ) THEN
	  
	       i_start = MAX(ids+2,its)
           i_end   = MIN(ide-2,ite)	
	 
	 
        j = jde-1
        DO k=kts,ktf
        DO i = i_start, i_end

   
           vb = v(i,k,j+1)
           IF (v(i,k,j) .GT. 0.) vb = v(i,k,j)
	   
              engk(i,k,j) = ( u(i,k,j)**2 + u(i+1,k,j)**2 + v(i,k,j)**2 + vb**2 + w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25  

        ENDDO
        ENDDO
      ENDIF
	  
!  西南角 	
	
      IF ( specified .AND. its .LE. ids+1 .AND. .NOT. config_flags%periodic_x .AND. jts .LE. jds+1  ) THEN
	 	 
	    i = ids 
        j = jds 
		
        DO k=kts,ktf
 

           ub = u(i ,k,j)
           vb = v(i,k,j)
		   IF (u(i+1,k,j) .LT. 0.) ub = u(i+1,k,j)
           IF (v(i,k,j+1) .LT. 0.) vb = v(i,k,j+1)   

              engk(i,k,j) = ( ub**2 + u(i+1,k,j)**2 + vb**2 + v(i,k,j+1)**2 + w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25  

 
        ENDDO
      ENDIF	  
	  
!  东南角 
    IF ( specified .AND. ite .GE. ide-1 .AND. .NOT. config_flags%periodic_x .AND. jts .LE. jds+1  ) THEN
	
	
       j = jds 
        DO k=kts,ktf
           i =  ide-1  
         
           ub = u(i+1 ,k,j)
		     vb = v(i,k,j)
			 
          IF (u(i,k,j) .GT. 0.) ub = u(i,k,j) 
		   IF (v(i,k,j+1) .LT. 0.) vb = v(i,k,j+1)    
		  
		  
            engk(i,k,j) = (  u(i,k,j)**2 + ub**2 +  vb**2 + v(i,k,j+1)**2 + w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25 
     
        ENDDO
		
      ENDIF		  
	  
	  
	  
!  西北角 	
	
      IF ( specified .AND. its .LE. ids+1 .AND. .NOT. config_flags%periodic_x .AND. jte .GE. jde-1  ) THEN
	 	 
	    i = ids 
        j = jde-1 
		
        DO k=kts,ktf
 

           ub = u(i ,k,j)
       
		   IF (u(i+1,k,j) .LT. 0.) ub = u(i+1,k,j)
     
		   
		   
		   vb = v(i,k,j+1)
           IF (v(i,k,j) .GT. 0.) vb = v(i,k,j)
		   
		   

              engk(i,k,j) = ( ub**2 + u(i+1,k,j)**2 + v(i,k,j)**2 + vb**2 +  w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25  

 
        ENDDO
      ENDIF	  	  


!  东北角 
    IF ( specified .AND. ite .GE. ide-1 .AND. .NOT. config_flags%periodic_x .AND. jte .GE. jde-1  ) THEN
	
	
       j = jde-1 
        DO k=kts,ktf
           i =  ide-1  
         
           ub = u(i+1 ,k,j)
 	 
          IF (u(i,k,j) .GT. 0.) ub = u(i,k,j) 
 


		   vb = v(i,k,j+1)
           IF (v(i,k,j) .GT. 0.) vb = v(i,k,j) 		   
		  
		  
            engk(i,k,j) = (  u(i,k,j)**2 + ub**2 +  vb**2 + v(i,k,j+1)**2 + w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25 
        ENDDO
   
		
      ENDIF		  	 
	  
	
endif ! 20240706	  
	
    END SUBROUTINE calcu_engk




	  
	  
!   calculate_sgmdot 子程序用于：计算sigma坐标的垂直速度sgmdot	  	  
!_______________________________________
	  
     SUBROUTINE calculate_sgmdot( sgmdot,                 &
                        u, v, w,                       &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   ! Input data
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) ::     u,     &
                                                                      v, &
                                                                      w
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      !   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: mut
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) :: sgmdot

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   ! Local data
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
		 
	     INTEGER :: i_endu, j_endv
				 
				 
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         
         
         REAL , DIMENSION( kms:kme ) :: sgmh,dsgm,dsgmh
      !   REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: zzh,dzzh
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: sx,sy !varz,varzy,varzys
         REAL :: c1x,c1y,c2x,c2y,mrdx,sxxs,syys
    LOGICAL :: specified

	specified = .false.
	if(config_flags%specified .or. config_flags%nested) specified = .true.


		  
   
     ktf=MIN(kte,kde-1)
 

 	  DO k=kts,ktf
          sgmh(k) = 0.5*(zh(k)+zh(k+1))
          dsgm(k) = zh(k+1) - zh(k) 
      ENDDO


      i_start = its
      i_end   = MIN(ide-1,ite)
	  
	  j_start = jts
      j_end   = MIN(jte,jde-1)
	  
	  
	  i_endu =  MIN(ide,ite+1)
	  j_endv =  MIN(jde,jte+1)
	  
	  ! 20231222 为了specified 边界等于0
	   sgmdot = 0.0 




! zs 对于周期和对称边界有值。

 
		  
		  sx = 0.0
		  sy = 0.0 
	

      i_start = its  
      i_end   = MIN(ide-1,ite)
	  
	  j_start = jts
      j_end   = MIN(jte,jde-1)
	  
      DO  j = j_start, j_end
      DO  k = kts, ktf
      DO  i = i_start, i_endu 
	  
	  ! 计算边界值，如果是物理边界，zs(i-1,j) 有值
	      c1x = (msft(i,j)*(ztop-sgmh(k))/(ztop-zs(i,j)) + msft(i-1,j)*(ztop-sgmh(k))/(ztop-zs(i-1,j)))*0.5
	      
	      c2x = u(i,k,j)*rdx*(zs(i,j)-zs(i-1,j))
	   
	      sx(i,k,j) =   c1x*c2x 

	
	  
	  ENDDO
      ENDDO
      ENDDO 



      DO  j = j_start, j_endv
      DO  k = kts, ktf
      DO  i = i_start, i_end 
	  
	  
	  
	   
	      c1y = (msft(i,j)*(ztop-sgmh(k))/(ztop-zs(i,j)) + msft(i,j-1)*(ztop-sgmh(k))/(ztop-zs(i,j-1)))*0.5
	  
	
	      c2y = v(i,k,j)*rdy*(zs(i,j)-zs(i,j-1))
	
	      sy(i,k,j) =   c1y*c2y 
		  
	
	  
	  ENDDO
      ENDDO
      ENDDO 
	  

	  sxxs = 0.0 
	  syys = 0.0 



   !   IF ( config_flags%open_xs .or. specified ) i_start = MAX(ids+1,its)
   !   IF ( config_flags%open_xe .or. specified ) i_end   = MIN(ide-1,ite)


      IF ( config_flags%periodic_x ) i_start = its
      IF ( config_flags%periodic_x ) i_end = ite


    !  IF ( config_flags%open_ys .or. specified ) j_start = MAX(jds+1,jts)
    !  IF ( config_flags%open_ye .or. specified ) j_end   = MIN(jde-1,jte)

	  
!20230204  ! 增加计算边界，暂时简单处理。sxxs,syys暂时为0，不适用实际地形。

      i_start = its
      i_end   = MIN(ide-1,ite)
	  
	  j_start = jts
      j_end   = MIN(jte,jde-1)	  
	  
	  
     DO  j = j_start, j_end
      DO  k = kts+1, ktf
      DO  i = i_start, i_end 
	  
    	!  sxxs = 0.5*((sx(i-1,k-1,j) + sx(i,k-1,j))*dsgm(k) +(sx(i-1,k,j) + sx(i,k,j))*dsgm(k-1) ) /(dsgm(k-1) + dsgm(k) ) 
	  
	     ! syys = 0.5*((sy(i,k-1,j-1) + sy(i,k-1,j))*dsgm(k) +(sy(i,k,j-1) + sy(i,k,j))*dsgm(k-1) ) /(dsgm(k-1) + dsgm(k) ) 
	  !    20220722
          sxxs = 0.5*((sx(i+1,k-1,j) + sx(i,k-1,j))*dsgm(k) +(sx(i+1,k,j) + sx(i,k,j))*dsgm(k-1) ) /(dsgm(k-1) + dsgm(k) ) 
	  
	      syys = 0.5*((sy(i,k-1,j+1) + sy(i,k-1,j))*dsgm(k) +(sy(i,k,j+1) + sy(i,k,j))*dsgm(k-1) ) /(dsgm(k-1) + dsgm(k) ) 
	  
	  
	      sgmdot(i,k,j) = - sxxs - syys + ztop/(ztop-zs(i,j))*w(i,k,j)
	  

		
	  ENDDO
      ENDDO
      ENDDO

  
 
  !   边界条件。
	  

      DO  j = j_start, j_end
      DO  i = i_start, i_end 	  
	  
	        sgmdot(i,kts,j) =  0. !
	       sgmdot(i,kte,j) =  0. !
		   
		   
      ENDDO
      ENDDO 

    !write(6,*)   ' sgmdot(1,2,1),w(1,2,1),sx(1,2,1),sy(1,2,1),msft(1,0),msft(0,1),zs(1,0),zs(0,1)='	  
	!write(6,*)    sgmdot(1,2,1),w(1,2,1),sx(1,2,1),sy(1,2,1),msft(1,0),msft(0,1),zs(1,0),zs(0,1)

      END SUBROUTINE calculate_sgmdot


SUBROUTINE calculate_sgmdot2( sgmdot,                 &
                        u, v, w,                       &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   ! Input data
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) ::     u,     &
                                                                      v, &
                                                                      w
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      !   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: mut
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) :: sgmdot

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   ! Local data
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         
         
         REAL , DIMENSION( kms:kme ) :: sgmh,dsgm,dsgmh
      !   REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: zzh,dzzh
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: sx,sy !varz,varzy,varzys
         REAL :: c1x,c1y,c2x,c2y,mrdx,sxxs,syys
    LOGICAL :: specified

	specified = .false.
	if(config_flags%specified .or. config_flags%nested) specified = .true.


		  
   
     ktf=MIN(kte,kde-1)
 

 	  DO k=kts,ktf
          sgmh(k) = 0.5*(zh(k)+zh(k+1))
          dsgm(k) = zh(k+1) - zh(k) 
      ENDDO


  
      i_start = its
      i_end   = ite
	  
	  j_start = jts
      j_end   = MIN(jte,jde-1)
	  

     ! IF ( config_flags%open_xs .or. config_flags%specified ) i_start = MAX(ids+1,its)
     ! IF ( config_flags%open_xe .or. config_flags%specified ) i_end   = MIN(ide-1,ite)

      IF ( config_flags%open_xs .or. specified ) i_start = MAX(ids+1,its)
      IF ( config_flags%open_xe .or. specified ) i_end   = MIN(ide-1,ite)


      IF ( config_flags%periodic_x ) i_start = its
      IF ( config_flags%periodic_x ) i_end = ite


      IF ( config_flags%open_ys .or. specified ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified ) j_end   = MIN(jde-1,jte)

 
		  
		  sx = 0.0
		  sy = 0.0 
		  
      DO  j = j_start, j_end
      DO  k = kts, ktf
      DO  i = i_start, i_end 
	  
	  
	  
	      c1x = (msft(i,j)*(ztop-sgmh(k))/(ztop-zs(i,j)) + msft(i-1,j)*(ztop-sgmh(k))/(ztop-zs(i-1,j)))*0.5
	      c1y = (msft(i,j)*(ztop-sgmh(k))/(ztop-zs(i,j)) + msft(i,j-1)*(ztop-sgmh(k))/(ztop-zs(i,j-1)))*0.5
	  
	      c2x = u(i,k,j)*rdx*(zs(i,j)-zs(i-1,j))
	      c2y = v(i,k,j)*rdy*(zs(i,j)-zs(i,j-1))
	  
	      sx(i,k,j) =   c1x*c2x 
	      sy(i,k,j) =   c1y*c2y 
		  
	
	  
	  ENDDO
      ENDDO
      ENDDO 
	  
  
	  sxxs = 0.0 
	  syys = 0.0 
	  
      DO  j = j_start, j_end
      DO  k = kts, kte
      DO  i = i_start, i_end 
	  
    	  sxxs = 0.5*((sx(i-1,k-1,j) + sx(i,k-1,j))*dsgm(k) +(sx(i-1,k,j) + sx(i,k,j))*dsgm(k-1) ) /(dsgm(k-1) + dsgm(k) ) 
	  
	      syys = 0.5*((sy(i,k-1,j-1) + sy(i,k-1,j))*dsgm(k) +(sy(i,k,j-1) + sy(i,k,j))*dsgm(k-1) ) /(dsgm(k-1) + dsgm(k) ) 
	  

	  
	      sgmdot(i,k,j) = - sxxs - syys + ztop/(ztop-zs(i,j))*w(i,k,j)
	  

		
	  ENDDO
      ENDDO
      ENDDO    	 

 		   
		   

      END SUBROUTINE calculate_sgmdot2
	  



SUBROUTINE calcu_divg_bak( var, divg,                 &
                        u, v, sgmdot,                       &
                     config_flags,                  &
                        msft,zh,zs,fzm,fzp,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   ! Input data
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var,u,     &
                                                                      v, &
                                                                      sgmdot
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      !   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: mut
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: divg

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh,fzm,fzp


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   ! Local data
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx
         
        
       !  REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: zzh,dzzh
	   
	   
	   
	   
	   
		 
		    REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) ::  varw
		 
		 
		 
		 


      REAL , DIMENSION( kms:kme ) :: sgmh,dsgm
       REAL , DIMENSION( ims:ime , jms:jme )  ::  rmsft,rmm
	   
	   REAL :: ub,vb 
				 
				 

    LOGICAL :: specified

   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

   
     ktf=MIN(kte,kde-1)
 

	  DO k=kts,kde-1
          sgmh(k) = 0.5*(zh(k)+zh(k+1))
          dsgm(k) = zh(k+1) - zh(k) 
      ENDDO


 



	  

   ! 20230228  for rho couple scalar 
      CALL set_physical_bc3d( var, 'p', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              its, ite, jts, jte, kts, kte,       &
                              its, ite, jts, jte, kts, kte )
	  
	  

     divg = 0.0 
	 
!  ---- x direction -----------
 

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)


      IF ( .NOT. config_flags%periodic_x ) THEN
        IF ( config_flags%open_xs .or. specified ) i_start = MAX(ids+1,its)
        IF ( config_flags%open_xe .or. specified ) i_end   = MIN(ide-2,ite)
      ENDIF

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
  
    !  	divg(i,k,j) =  -  (  0.5*((ztop-zs(i,j))*var(i,k,j)*rmsft(i,j)  + (ztop-zs(i+1,j))*var(i+1,k,j)*rmsft(i+1,j) ) *u(i+1,k,j) &
    !                -   0.5*((ztop-zs(i-1,j))*var(i-1,k,j)*rmsft(i-1,j)  + (ztop-zs(i,j))*var(i,k,j)*rmsft(i,j) ) *u(i,k,j) )*rdx

      	divg(i,k,j) =  -  (  0.5*((ztop-zs(i,j))*var(i,k,j)/msft(i,j)  + (ztop-zs(i+1,j))*var(i+1,k,j)/msft(i+1,j) ) *u(i+1,k,j) &
                     -   0.5*((ztop-zs(i-1,j))*var(i-1,k,j)/msft(i-1,j)  + (ztop-zs(i,j))*var(i,k,j)/msft(i,j) ) *u(i,k,j) )*rdx
      ENDDO
      ENDDO
      ENDDO
	  
	  
	IF ( config_flags%open_xs  .and.  its == ids ) THEN 
		
	!	i = ids  !  用于方法1 
		
      DO  j = j_start, j_end
      DO k=kts,ktf
  
  
 !  方法1 ，直接改造原来的格式。不用i-1点的var
     ! 	divg(i,k,j) =  -  (  0.5*((ztop-zs(i,j))*var(i,k,j)/msft(i,j)  + (ztop-zs(i+1,j))*var(i+1,k,j)/msft(i+1,j) ) *u(i+1,k,j) &
     !                -    (  (ztop-zs(i,j))*var(i,k,j)/msft(i,j) ) *u(i,k,j) )*rdx


					 
	! 参考wrf	,这一方法的使用，改进非常大。避免了hill2km_x试验的西边界问题。			 
	         ub = MIN( 0.5*(u(its,k,j)+u(its+1,k,j)), 0. )
		 
         divg(its,k,j) =                      &
               - (ztop-zs(its,j))/msft(its,j)*rdx*(                                             &
                       ub*(   var(its+1,k,j)                 &
                            - var(its  ,k,j)   ) +           &
                       var(its,k,j)*(u(its+1,k,j)-u(its,k,j))  &       ! modified 
                                                                )				 
					 
      ENDDO
      ENDDO
 
	  
	ENDIF	  
	

	  
	    IF ( config_flags%open_xe  .and.  ite == ide ) THEN 
		
	 !	i = ite-1
		
		 i_end   = MIN(ite,ide-1)
		 
		
      DO  j = j_start, j_end
      DO k=kts,ktf
 
  
      ! 	divg(i,k,j) =  -  (  ((ztop-zs(i,j))*var(i,k,j)/msft(i,j)    ) *u(i+1,k,j) &
      !                -   0.5*((ztop-zs(i-1,j))*var(i-1,k,j)/msft(i-1,j)  + (ztop-zs(i,j))*var(i,k,j)/msft(i,j) ) *u(i,k,j) )*rdx
	
	     ub = MAX( 0.5*(u(ite-1,k,j)+u(ite,k,j)), 0. )	 
		 
         divg(i_end,k,j) =                  &
               - (ztop-zs(i_end,j))/msft(i_end,j)*rdx*(                                               &
                       ub*(  var(i_end  ,k,j)                  &
                           - var(i_end-1,k,j) ) +              &
                       var(i_end,k,j)*(u(ite,k,j)-u(ite-1,k,j))  &   ! modified  
                                                                    )				 
					 
					 
      ENDDO
      ENDDO
  

	ENDIF	
	  
	  
	  
	  

!  ---- y direction -----------

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
      ! Polar boundary conditions are like open or specified
      IF ( config_flags%open_ys .or. specified .or. config_flags%polar ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified .or. config_flags%polar ) j_end   = MIN(jde-2,jte) 	  

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end

   !   	divg(i,k,j) = divg(i,k,j) -   ( 0.5*((ztop-zs(i,j))*var(i,k,j)*rmsft(i,j)  + (ztop-zs(i,j+1))*var(i,k,j+1)*rmsft(i,j+1) ) *v(i,k,j+1)    &
   !                  -   0.5*((ztop-zs(i,j-1))*var(i,k,j-1)*rmsft(i,j-1)  + (ztop-zs(i,j))*var(i,k,j)*rmsft(i,j) ) *v(i,k,j) )*rdy
					 
 		 
    	divg(i,k,j) = divg(i,k,j) -   ( 0.5*((ztop-zs(i,j))*var(i,k,j)/msft(i,j)  + (ztop-zs(i,j+1))*var(i,k,j+1)/msft(i,j+1) ) *v(i,k,j+1)    &
                 -   0.5*((ztop-zs(i,j-1))*var(i,k,j-1)/msft(i,j-1)  + (ztop-zs(i,j))*var(i,k,j)/msft(i,j) ) *v(i,k,j) )*rdy					 

   

       ENDDO
      ENDDO
      ENDDO
	  
	  

	  
	    IF ( config_flags%open_ys  .and.  jts == jds ) THEN 
		
	!	j = jds 
		
      DO k=kts,ktf
      DO  i = i_start, i_end		
		
    !  	divg(i,k,j) = divg(i,k,j) -   ( 0.5*((ztop-zs(i,j))*var(i,k,j)/msft(i,j)  + (ztop-zs(i,j+1))*var(i,k,j+1)/msft(i,j+1) ) *v(i,k,j+1)    &
    !                 -    ( (ztop-zs(i,j))*var(i,k,j)/msft(i,j) ) *v(i,k,j) )*rdy
	
	
	
	      vb = MIN( 0.5*(v(i,k,jts)+v(i,k,jts+1)), 0. )
		  
        divg(i,k,jts) =    divg(i,k,jts)                &
               - (ztop-zs(i,jts))/msft(i,jts)*rdy*(                                             &
                       vb*(  var(i,k,jts+1)                  &
                           - var(i,k,jts  ) ) +              &
                       var(i,k,jts)*(v(i,k,jts+1)-v(i,k,jts))  &     ! modified  
                                                               )
					 
      ENDDO
      ENDDO
	    ENDIF	
		
 
	  
	    IF ( config_flags%open_ye  .and.  jte == jde ) THEN 
		
	 ! 	j = jte-1
		
		j_end   = MIN(jte,jde-1)
		
      DO k=kts,ktf
      DO  i = i_start, i_end		
		
	  !      divg(i,k,j) = divg(i,k,j) -   (  ((ztop-zs(i,j))*var(i,k,j)/msft(i,j) )   *v(i,k,j+1)    &
       !              -   0.5*((ztop-zs(i,j-1))*var(i,k,j-1)/msft(i,j-1)  + (ztop-zs(i,j))*var(i,k,j)/msft(i,j) ) *v(i,k,j) )*rdy					 
	
		  vb = MAX( 0.5*(v(i,k,jte-1)+v(i,k,jte)), 0. ) 
		 
         divg(i,k,j_end) =    divg(i,k,j_end)               &
               -  (ztop-zs(i,j_end))/msft(i,j_end)*rdy*(                                               &
                       vb*(   var(i,k,j_end  )                 &
                            - var(i,k,j_end-1) ) +             &
                       var(i,k,j_end)*(v(i,k,jte)-v(i,k,jte-1))  &    ! modified   
                                                                    )
					 
      ENDDO
      ENDDO
	    ENDIF

		
	  
	  
	  	call average_vertical_w( varw, var,                 &
                     config_flags,                  &
                         fzm,fzp,                    &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )	


!  ---- z  direction -----------

  !  计划拆分为简单公式，处理变量从半层插值到整层，及边界处理问题
  
  
        i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
	  
      DO  j = j_start, j_end
      DO  i = i_start, i_end
        rmsft(i,j)  = 1./msft(i,j) 
        rmm(i,j) = 1./(msft(i,j)*msft(i,j) )
      ENDDO
      ENDDO             

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end



            divg(i,k,j) = divg(i,k,j) -    ( varw(i,k+1,j)* rmm(i,j) * sgmdot(i,k+1,j) - varw(i,k,j)* rmm(i,j) * sgmdot(i,k,j)   ) /dsgm(k) * (ztop-zs(i,j))


      ENDDO
      ENDDO
      ENDDO



      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
         	divg(i,k,j) = divg(i,k,j)  * msft(i,j)*msft(i,j)/(ztop-zs(i,j)) 
      ENDDO
      ENDDO
      ENDDO

      END SUBROUTINE calcu_divg_bak






!  vortx,vorty,vortz 为求取 涡度相关项的子程序 

			
   
      SUBROUTINE calcu_vorty( vorty,                 &
                        u, v, w,                       &
                     config_flags,                  &
                        msft,msfty,xlat,zh,zs,fzm,fzp,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   ! Input data

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: u,     &
                                                                      v,w
																	  
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      !   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: mut
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: vorty

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,msfty,xlat,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh,fzm,fzp


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   ! Local data
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx
  
       !  REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: zzh,dzzh
	   ! 20240706 多余
       !  REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: varz,varzx,varzxs

       REAL , DIMENSION( kms:kme ) :: sgmh,dsgm,dsgmh
       REAL , DIMENSION( ims:ime , jms:jme )  ::  rmsft,rmm
	   
	   REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  ::  uw,uwy,vds,vss,udxzs,tmp,wvdsh !vs, 

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme )   :: vorty1,vorty2,vorty21, vorty22
		 
		 
       LOGICAL :: specified
	   
	   integer,dimension(3) ::  fdim(3)
	  
 
	   vorty = 0.0 
	   vorty1 = 0.0
	   vorty2 = 0.0
	   vorty21 = 0.0 
	   vorty22 = 0.0 

! 1) 判定 是否 给定边界
   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

 !2) 设定 do循环的 垂直 层次   
	

	 
       ktf=MIN(kte,kde-1)

	  DO k=kts,kde-1
          sgmh(k) = 0.5*(zh(k)+zh(k+1))
          dsgm(k) = zh(k+1) - zh(k) 
      ENDDO

 	  DO k=2,kde-1
             dsgmh(k) = 0.5*(dsgm(k) + dsgm(k-1))   ! new 20220711
			 
      ENDDO
	  

       
      i_start = its
      i_end   =  MIN(ite,ide-1)

      j_start = jts
      j_end   = MIN(jte,jde-1)
	  



if(.false.) then 
      DO  j = j_start, j_end
      DO  k = kts, ktf
      DO  i = i_start, i_end


   
       uw(i,k,j) = ( u(i,k,j)**2 + u(i+1,k,j)**2 + w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25 
	  
      ENDDO
      ENDDO
      ENDDO


endif 

call    calcu_engk_uw(uw,u,v,w,config_flags,                       &
                                  ids, ide, jds, jde, kds, kde,  &
                                  ims, ime, jms, jme, kms, kme,  &
                                  its, ite, jts, jte, kts, kte  )	  

  

! 1) 第一项 


	uwy = 0. 
	
		!write(6,*)  'uw grady' 
	
      call  calcu_grady( uw, uwy,                 &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )



      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
     
      	vorty1(i,k,j)=  msfty(i,j)*(-1.0*uwy(i,k,j))   

      ENDDO
      ENDDO
      ENDDO
		  
! 2) 第二项 

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = jte

    !  IF ( config_flags%open_ys ) j_start = MAX(jds+1,jts)
     ! IF ( config_flags%open_ye ) j_end   = MIN(jde-1,jte)
     ! IF ( specified ) j_start = MAX(jds+1,jts)
     ! IF ( specified ) j_end   = MIN(jde-1,jte)
	 
	 

	
	! 周期边界
   !   IF ( config_flags%periodic_x ) i_start = its 
    !  IF ( config_flags%periodic_x ) i_end = ite
	
       IF ( config_flags%periodic_y ) j_start = jts           ! 周期边界，因为设置了外围，可以计算全部its至ite范围。
      IF ( config_flags%periodic_y ) j_end = jte	  
	  

	  
	  
!      南北边界另外单算。 
      IF ( config_flags%open_ys .or. specified ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified ) j_end   = MIN(jde-1,jte)
	  
	!  IF ( config_flags%polar ) j_start = MAX(jds+1,jts)
    !  IF ( config_flags%polar ) j_end   = MIN(jde-1,jte)  
	  
 !   开边界、给定边界
!     东西边界只有vorty21, vorty22另外算。
      IF ( config_flags%open_xs  .or. specified ) i_start = MAX(ids+1,its)   ! 开边界的ids,ide, 需要特别处理
      IF ( config_flags%open_xe  .or. specified ) i_end   = MIN(ide-1,ite)
	  

	  
		  
      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
    
         !   vorty(i,k,j) = vorty(i,k,j) + 0.5*( (u(i,k,j)+u(i,k,j+1))*0.5 * (v(i,k,j) - v(i-1,k,j))*rdx   &
         !                                          +  (u(i+1,k,j)+u(i+1,k,j+1))*0.5 * (v(i+1,k,j) - v(i,k,j))*rdx  )
		 
		 ! 20221204
             vorty21(i,k,j) =   0.5*( (u(i,k,j)+u(i,k,j-1))*0.5 * (v(i,k,j) - v(i-1,k,j))*rdx   &
                                                   +  (u(i+1,k,j)+u(i+1,k,j-1))*0.5 * (v(i+1,k,j) - v(i,k,j))*rdx  )       

      ENDDO
      ENDDO
      ENDDO
	  
	  
    !  IF ( specified .AND. jts .LE. jds+1 ) THEN
     !   j = jds+1
	!	DO k=kts,ktf
	!	DO  i = i_start, i_end
     !         vb = v(i,k,j-1) 
       
		!ENDDO
		!ENDDO
	  
      !ENDIF 
	  
	  
	  


	   
    !  DO  j = j_start, j_end
     ! DO k=kts+1,ktf
     ! DO  i = i_start, i_end
     
     ! 	vs(i,k,j) = fzm(k)*v(i,k,j) + fzp(k)*v(i,k-1,j)

     ! ENDDO
    !  ENDDO
	  
     ! DO  i = i_start, i_end
    !    vs(i,kts,j) = vs(i,kts+1,j)    ! low  boundary
    !    vs(i,kte,j) = vs(i,ktf,j) 
    !  ENDDO
		
    !  ENDDO


! 20221204
    !  i_start = its
    !  i_end   =  MIN(ite,ide-1)

    !  j_start = jts
    !  j_end   = jte
	  


      DO  j = j_start, j_end
      DO k=2,ktf 
      DO  i = i_start, i_end  
    ! dv/dsgm
      	vds(i,k,j) = ( v(i,k,j) - v(i,k-1,j) )/dsgmh(k) 
	
		
      ENDDO
      ENDDO
      ENDDO	 
	  
	  
	 ! 上下边界 
	  DO  j = j_start, j_end
      DO  i = i_start, i_end

    ! dv/dsgm
      	vds(i,kte,j) = vds(i,kte-1 ,j) 
        vds(i,1,j) =  2.0*( v(i,1,j) - 0.0 )/dsgmh(2)   ! surface is 0       !vds(i,kts+1,j)
	
		
      ENDDO
      ENDDO
	  
	       

	   
	   vss = 0.0 

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
     
      	!vss(i,k,j) = 0.5*( vs(i,k,j) + vs(i,k+1,j) ) 
          vss(i,k,j) = 0.5*( vds(i,k,j) + vds(i,k+1,j) ) 
      ENDDO
      ENDDO
      ENDDO
	  

	  
 !   开边界、给定边界
!     东西边界 vorty21, vorty22另外算。	  
!	  IF ( config_flags%open_xs  .or. specified ) i_start = MAX(ids+1,its)   
 !     IF ( config_flags%open_xe  .or. specified ) i_end   = MIN(ide-1,ite)  
	  
	udxzs = 0.0 
	! zs 的水平梯度处理需要优化，为不变量。
	
			j_end =  jte+1   ! 20240706    
	  
      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
     
      	udxzs(i,k,j) = 0.5*( u(i,k,j)*(zs(i,j)-zs(i-1,j))*rdx + u(i+1,k,j)*(zs(i+1,j)-zs(i,j))*rdx )

      ENDDO
      ENDDO
      ENDDO


	  
      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
     
      	udxzs(i,k,j) = (ztop-sgmh(k))/(ztop-zs(i,j))*udxzs(i,k,j)

      ENDDO
      ENDDO
      ENDDO
	  
	  
	        CALL set_physical_bc3d( udxzs, 'v', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              its, ite, jts, jte, kts, kte,       &
                              its, ite, jts, jte, kts, kte ) 
							  
							  


				j_end =  jte   ! 20240706   

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
     
      	vorty22(i,k,j)=    - 0.5*( udxzs(i,k,j)  + udxzs(i,k,j+1)  ) * vss(i,k,j) 
!                                 
      ENDDO
      ENDDO
      ENDDO


      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
     
      	vorty2(i,k,j)= - msfty(i,j)*(vorty21(i,k,j) + vorty22(i,k,j))   

      ENDDO
      ENDDO
      ENDDO
	  
	  
	  
	  
	  
	  



   !   j_start = jts
   !   j_end   = jte
	  

	  
	  
! vort3 
     i_start = its
     i_end   =  MIN(ite,ide-1)
	 
	 
   !  j_start = jts
   !  j_end   = MIN(jte,jde-1)

	  
    !  IF ( config_flags%open_ys .or. specified ) j_start = MAX(jds+1,jts)
    !  IF ( config_flags%open_ye .or. specified ) j_end   = MIN(jde-1,jte)

	  
	  DO  j = j_start, j_end
      DO k=kts,ktf+1
      DO  i = i_start, i_end
         	tmp(i,k,j) = 0.5*( w(i,k,j-1) + w(i,k,j) )  * vds(i,k,j) 
			

		 
      ENDDO
      ENDDO
      ENDDO  
		  

	    
 	  DO  j = j_start, j_end
      DO k=kts,ktf 
      DO  i = i_start, i_end 
          	wvdsh(i,k,j) =  (tmp(i,k,j) + tmp(i,k+1,j) ) *0.5 
			
      ENDDO
      ENDDO
      ENDDO   

!垂直方向

	  !  可以设置可分配动态数组释放一些数据。

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
 
      	vorty(i,k,j) = vorty1(i,k,j) + vorty2(i,k,j)  -   (ztop/(ztop-zs(i,j)) + ztop/(ztop-zs(i,j-1)))*0.5*wvdsh(i,k,j)

      ENDDO
      ENDDO
      ENDDO	  
	  

! 注意实例模拟时，取值范围。		
		  
      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end

      ! u(i,j)           u(i+1，j)
      !          v(i,j)
      !
      ! u(i,j-1)         u(i+1,j-1） 


       
                !     m(i,j)
                !     v(i,j) 
                !      m(i,j-1) 

  
        vorty(i,k,j) = vorty(i,k,j)  - 0.25* ( u(i+1,k,j)*u(i+1,k,j)   + u(i,k,j)*u(i,k,j)  + u(i,k,j-1)*u(i,k,j-1)  + u(i,k,j-1)*u(i,k,j-1)  )  & 
                      *(msft(i,j)-msft(i,j-1))*rdy    &
                                        + v(i,k,j) * ( u(i,k,j)*(msft(i,j)-msft(i-1,j))*rdx   &
                                                   +  u(i,k,j-1)*(msft(i,j-1)-msft(i-1,j-1))*rdx   &
                                                   +  u(i+1,k,j-1)*(msft(i+1,j-1)-msft(i,j-1))*rdx    &
                                                   +  u(i+1,k,j)*(msft(i+1,j)-msft(i,j))*rdx  ) * 0.25 
     



      ENDDO
      ENDDO
      ENDDO	  

  
      if(.true.)then   ! false for ideal test  now 20220830, 202307,20240111
	  
      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
  
          vorty(i,k,j) = vorty(i,k,j)  -  reradius*(  v(i,k,j)*0.25*(w(i,k,j)+w(i,k,j+1)+w(i,k+1,j)+w(i,k+1,j+1) )          &
          
           + 0.5*( (u(i,k,j-1)*u(i,k,j-1)+u(i+1,k,j-1)*u(i+1,k,j-1))*0.5*tan(xlat(i,j)*degrad) +  (u(i,k,j)*u(i,k,j)+u(i+1,k,j)*u(i+1,k,j))*0.5*tan(xlat(i,j)*degrad) ) ) 
                                       
       

      ENDDO
      ENDDO
      ENDDO	
	  
     endif
 
      END SUBROUTINE calcu_vorty






      
!  vortx,vorty,vortz 为求取 涡度相关项的子程序 
   
      SUBROUTINE calcu_vortz(  vortz,                 &
                        u, v, w,                      &
                     config_flags,                  &
                        msft,xlat,zh,zs,fzm,fzp,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   ! Input data

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: u,v,w
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      !   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: mut
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: vortz

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,xlat,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh,fzm,fzp


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop
		 
		 
		    REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: ufull

   ! Local data
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx
  


       REAL , DIMENSION( kms:kme ) :: sgmh,dsgm,dsgmh
       REAL , DIMENSION( ims:ime , jms:jme )  ::  rmsft,rmm
	   
       ! 下列数组，可以设置些动态数组，减少内存占用。
	   REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  ::  uv,vortz1,vortz21,udxzs,temp,vortz22,vdyzs,vortz31,vortz32,vortz4


       LOGICAL :: specified

! 1) 判定 是否 给定边界
   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

 !2) 设定 do循环的 垂直 层次   
      ktf=MIN(kte,kde-1)
	  
 

	  DO k=kts,kde-1
          sgmh(k) = 0.5*(zh(k)+zh(k+1))
          dsgm(k) = zh(k+1) - zh(k) 
      ENDDO

 	  DO k=2,kde-1
             dsgmh(k) = 0.5*(dsgm(k) + dsgm(k-1))   ! new 20220711
			 
      ENDDO

	  
       
      i_start = its
      i_end   =  MIN(ite,ide-1)

      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
	  
	  
if (.false.) then 

        uv = 0.0 
      DO      j = j_start, j_end
      DO  k = kts, ktf
      DO  i = i_start, i_end

      

              uv(i,k,j) = ( u(i,k,j)**2 + u(i+1,k,j)**2 + v(i,k,j)**2 + v(i,k,j+1)**2 ) * 0.25 
     
      ENDDO
      ENDDO
      ENDDO

endif 

call    calcu_engk_uv(uv,u,v,w,config_flags,                       &
                                  ids, ide, jds, jde, kds, kde,  &
                                  ims, ime, jms, jme, kms, kme,  &
                                  its, ite, jts, jte, kts, kte  ) 
	  

! 1) 第一项 
	 

		vortz1 = 0.0 
		! last will couple -1.0 
        call  calcu_gradz( uv, vortz1,                    &
                      config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )



! 2) 第二项 

		 	call average_vertical_w_x( ufull, u,                 &
                     config_flags,                  &
                         fzm,fzp,                    &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )	

  !IF ( config_flags%periodic_x )  then 
						
      CALL set_physical_bc3d( ufull, 'u', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              its, ite, jts, jte, kts, kte,       &
                              its, ite, jts, jte, kts, kte )	

!endif							  
		
	  vortz21 = 0.0 
		
      DO  j = j_start, j_end
	  
	  
      DO k=1,kte 
      DO  i = i_start, i_end
    
            vortz21(i,k,j) =  0.5*(ufull(i,k,j) * (w(i,k,j)-w(i-1,k,j))*rdx                               &
                                                   + ufull(i+1,k,j) * (w(i+1,k,j)-w(i,k,j))*rdx ) 

        
           

      ENDDO
      ENDDO
	  


      ENDDO


        udxzs = 0.0 
		

		
      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
       ! ud(zs)/dx
      	udxzs(i,k,j) = 0.5*( u(i,k,j)*(zs(i,j)-zs(i-1,j))*rdx + u(i+1,k,j)*(zs(i+1,j)-zs(i,j))*rdx )

      ENDDO
      ENDDO
      ENDDO
	  

		
      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
     
      	temp(i,k,j) = (ztop-sgmh(k))/(ztop-zs(i,j))*udxzs(i,k,j) * (w(i,k+1,j)-w(i,k,j))/dsgm(k)

      ENDDO
      ENDDO
      ENDDO



      vortz22 = 0.0
      DO  j = j_start, j_end
      DO k=kts+1,ktf
      DO  i = i_start, i_end
     
      	vortz22(i,k,j) = fzm(k)*temp(i,k,j) + fzp(k)*temp(i,k-1,j)

      ENDDO
      ENDDO
      DO  i = i_start, i_end
        vortz22(i,kts,j) = vortz22(i,kts+1,j)      ! ??
        vortz22(i,kte,j) = vortz22(i,ktf,j) 
		
      ENDDO
		
      ENDDO


      
! 3) vortz第3项 
       !1) part1 
	   
     vortz31 = 0.0 
      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
    
            vortz31(i,k,j) =  0.5*((v(i,k,j)+v(i,k+1,j))*0.5 * (w(i,k,j)-w(i,k,j-1))*rdy                           &
                                                   + (v(i,k,j+1)+v(i,k+1,j+1))*0.5 * (w(i,k,j+1)-w(i,k,j))*rdy )
        

      ENDDO
      ENDDO
      ENDDO
	  


      !   v*d(zs)/dy * dw 从半层垂直插值到整层
	  vdyzs = 0.0 
	  
      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
      ! v*d(zs)/dy

      ! y方向跳点： v(i,j-1)-----w(i,j-1)----v(i,j)---w(i,j)---v(i,j+1)----w(i,j+1)
    
            vdyzs(i,k,j) =  (ztop-sgmh(k))/(ztop-zs(i,j))* ( v(i,k,j)* ( zs(i,j)- zs(i,j-1) )*rdy +  v(i,k,j+1)* ( zs(i,j+1)- zs(i,j) )*rdy)*0.5   &
                              * (w(i,k+1,j)-w(i,k,j))/dsgm(k)
							  
							  


      ENDDO
      ENDDO
      ENDDO
	  

			

			

      !从半层垂直插值到整层
	  
	  vortz32 = 0.0  
      DO  j = j_start, j_end
      DO k=kts+1,ktf
      DO  i = i_start, i_end
            vortz32(i,k,j) =  fzm(k)*vdyzs(i,k,j)  + fzp(k)*vdyzs(i,k-1,j)   
      ENDDO
      ENDDO
	  
      DO  i = i_start, i_end
        vortz32(i,kts,j) =  vortz32(i,kts+1,j)
        vortz32(i,kte,j) =  vortz32(i,ktf,j) 
      ENDDO
	  
      ENDDO
	
      ! vortz第四项，曲率项

			
			
	  vortz4 = 0.0 		
			
      DO  j = j_start, j_end
      DO k=kts+1,ktf
      DO  i = i_start, i_end
            vortz4(i,k,j) =  fzm(k)*uv(i,k,j)*2.0  + fzp(k)*uv(i,k-1,j)*2.0   ! 动能的2倍。 20220722      
      ENDDO
      ENDDO
	  
      DO  i = i_start, i_end
        vortz4(i,kts,j) =  vortz4(i,kts+1,j)
        vortz4(i,kte,j) =  vortz4(i,ktf,j) 
      ENDDO
	  
      ENDDO


	  !  可以设置可分配动态数组释放一些数据。
	  
	  
 	
       vortz = 0.0 	
	  
      DO  j = j_start, j_end
      DO k=kts,ktf+1 
      DO  i = i_start, i_end
 
      	vortz(i,k,j) =  -1.0*vortz1(i,k,j) -msft(i,j)*(vortz21(i,k,j)-vortz22(i,k,j) + vortz31(i,k,j) - vortz32(i,k,j)  )   &    
			  + reradius*vortz4(i,k,j)         ! 20240226 

      ENDDO
      ENDDO
      ENDDO	  
	  



      END SUBROUTINE calcu_vortz







	
	

	  	  
!  vortx,vorty,vortz 为求取 涡度相关项的子程序  
   
      SUBROUTINE calcu_vortx( vortx,                 &
                        u, v, w,                    &
                     config_flags,                  &
                        msft,xlat,zh,zs,fzm,fzp,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   ! Input data
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) ::u,v,w
                                                                     
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      !   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: mut
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT( OUT) :: vortx

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,xlat,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh,fzm,fzp


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   ! Local data
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx

         
        
       !  REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: zzh,dzzh
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: vortx1,vortx21,vortx22 


      REAL , DIMENSION( kms:kme ) :: sgmh,dsgm,dsgmh
       REAL , DIMENSION( ims:ime , jms:jme )  ::  rmsft,rmm
	   
	   REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  ::  vw,us,tmp,vwx,uds,wudsh

    LOGICAL :: specified



! 1) 判定 是否 给定边界
   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

 !2) 设定 do循环的 垂直 层次   
      ktf=MIN(kte,kde-1)
	  
 
	 
	  DO k=kts,kde-1
          sgmh(k) = 0.5*(zh(k)+zh(k+1))
          dsgm(k) = zh(k+1) - zh(k) 
      ENDDO

 	  DO k=2,kde-1
             dsgmh(k) = 0.5*(dsgm(k) + dsgm(k-1))   ! new 20220711
			 
      ENDDO
	  
	 ! 3) 设定  x、 y do循环的 范围 
       
      i_start = its
      i_end   =  MIN(ite,ide-1)

      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
	  
if(.false.) then   

      DO      j = j_start, j_end
      DO  k = kts, ktf
      DO  i = i_start, i_end

      
     !   vw(i,k,j) = ( u(i,k,j)*u(i,k,j) + u(i+1,k,j)*u(i+1,k,j) + v(i,k,j)*v(i,k,j) + v(i,k,j+1)*v(i,k,j+1) + w(i,k,j)*w(i,k,j) + w(i,k+1,j)*w(i,k+1,j) ) * 0.5 
      !   vw(i,k,j) = (v(i,k,j)*v(i,k,j) + v(i,k,j+1)*v(i,k,j+1) + w(i,k,j)*w(i,k,j) + w(i,k+1,j)*w(i,k+1,j) ) * 0.5 
        vw(i,k,j) = ( v(i,k,j)**2 + v(i,k,j+1)**2 + w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25 
     
      ENDDO
      ENDDO
      ENDDO

endif 
	  
	  
      call    calcu_engk_vw(vw,u,v,w,config_flags,                       &
                                  ids, ide, jds, jde, kds, kde,  &
                                  ims, ime, jms, jme, kms, kme,  &
                                  its, ite, jts, jte, kts, kte  )


	  
 
        vortx1  = 0.0 
	! 注意梯度有负号，用是有时需要变正，p方程和旋转项。	
      call  calcu_gradx( vw, vortx1,                 &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )


	  
! 4)  修订 开边界、给定边界、循环边界的x方向的 do循环 的范围 


      i_start = its
      i_end   = ite  


      IF ( config_flags%open_xs ) i_start = MAX(ids+1,its)   ! 开边界的ids,ide, 需要特别处理
      IF ( config_flags%open_xe ) i_end   = MIN(ide-1,ite)
	  
	  
    !  IF ( specified ) i_start = MAX(ids+1,its)            ! 给定边界的，ids，ide，不做处理，ids+1, ide-1  二阶算法。
    !  IF ( specified ) i_end   = MIN(ide-1,ite)
	! 20240303 
       IF ( specified ) i_start = MAX(ids+2,its)            ! 给定边界的，次外层也需要单独处理，进行迎风计算。
       IF ( specified ) i_end   = MIN(ide-2,ite)
	
	  
      IF ( config_flags%periodic_x ) i_start = its           ! 周期边界，因为设置了外围，可以计算全部its至ite范围。
      IF ( config_flags%periodic_x ) i_end = ite

! 南北边界，单独处理，vort1，vort3不单独处理
      IF ( config_flags%open_ys .or. specified ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified ) j_end   = MIN(jde-1,jte)
	  
	  



     vortx21  = 0.0 
         DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end

    
    
      	vortx21(i,k,j) =  0.5*( (v(i-1,k,j+1)+ v(i,k,j+1))*0.5 * (u(i,k,j+1) - u(i,k,j) )*rdy   &
                             +  (v(i-1,k,j)  + v(i,k,j) )*0.5  * (u(i,k,j)   - u(i,k,j-1) )*rdy )

        
        ! j+1
        ! j   v_xave(i,j) :  (v(i-1,k,j)+v(i,k,j))*0.5

        ! u(i,j) is u(i-1/2,j)
        ! v(i,j) is v(i,j-1/2)
        ! T(i,j) is T(i,j)


        !  下： v_xave(i,j)*dyU(i,j) :   (v(i-1,k,j)+v(i,k,j))*0.5 * (u(i,k,j) - u(i,k,j-1))*rdy 
         !上：  v_xave(i,j+1)*dyU(i,j+1) :   (v(i-1,k,j+1)+v(i,k,j+1))*0.5 * (u(i,k,j+1) - u(i,k,j))*rdy 

         !  0.5*( (v(i-1,k,j+1)+v(i,k,j+1))*0.5 * (u(i,k,j+1) - u(i,k,j))*rdy  +  (v(i-1,k,j)+v(i,k,j))*0.5 * (u(i,k,j) - u(i,k,j-1))*rdy )
      ENDDO
      ENDDO
      ENDDO







	   vortx22  = 0.0 
	  
      DO  j = j_start, j_end
	  
	  
      DO k=2,ktf-1
      DO  i = i_start, i_end

        
        ! 右侧    (ztop-sgmh(k))/(ztop-zs(i,j))* 0.5*(v(i,k,j+1)*(zs(i,j+1) - zs(i,j))*rdy + v(i,k,j)*(zs(i,j) - zs(i,j-1))*rdy)

           ! 左侧 (ztop-sgmh(k))/(ztop-zs(i-1,j)) *0.5*(v(i-1,k,j+1)*(zs(i-1,j+1) - zs(i-1,j))*rdy + v(i-1,k,j)*(zs(i-1,j) - zs(i-1,j-1))*rdy) 


           vortx22(i,k,j) =  - 0.5*((ztop-sgmh(k))/(ztop-zs(i,j) ) *0.5*(v(i,k,j+1)  *(zs(i,j+1)   - zs(i,j)  )*rdy + v(i,k,j)  *(zs(i,j)   - zs(i,j-1)   )*rdy)   &
                              +    (ztop-sgmh(k))/(ztop-zs(i-1,j)) *0.5*(v(i-1,k,j+1)*(zs(i-1,j+1) - zs(i-1,j))*rdy + v(i-1,k,j)*(zs(i-1,j) - zs(i-1,j-1) )*rdy)   )   &
                   *  0.5*  ( (u(i,k,j)   - u(i,k-1,j) )/dsgmh(k)  & 
				            + (u(i,k+1,j) - u(i,k,j) )/dsgmh(k+1) )    ! delat_sgm_u(i,k,j)_ave_sgm   



      ENDDO
      ENDDO
	!  下边界  
	   k = kts
	   DO  i = i_start, i_end  
	              vortx22(i,k,j) =  - 0.5*(   &
                  (ztop-sgmh(k))/(ztop-zs(i,j))   *0.5*(v(i,k,j+1)*(zs(i,j+1)     - zs(i,j))  *rdy + v(i,k,j)  *(zs(i,j)   - zs(i,j-1))*rdy)   &
             +    (ztop-sgmh(k))/(ztop-zs(i-1,j)) *0.5*(v(i-1,k,j+1)*(zs(i-1,j+1) - zs(i-1,j))*rdy + v(i-1,k,j)*(zs(i-1,j) - zs(i-1,j-1))*rdy) &
               )   &
                   *  0.5*  (2.0*(u(i,k,j) - 0.  )/dsgm(kts)    &
                               + (u(i,k+1,j) - u(i,k,j) )/dsgmh(k+1) )    ! delat_sgm_u(i,k,j)_ave_sgm   

	    enddo 

	!  上边界  du(kte) = du(kte-1)
	   k = ktf
	   DO  i = i_start, i_end  
           vortx22(i,k,j) =  - 0.5*((ztop-sgmh(k))/(ztop-zs(i,j) ) *0.5*(v(i,k,j+1)  *(zs(i,j+1)   - zs(i,j)  )*rdy + v(i,k,j)  *(zs(i,j)   - zs(i,j-1)   )*rdy)   &
                              +    (ztop-sgmh(k))/(ztop-zs(i-1,j)) *0.5*(v(i-1,k,j+1)*(zs(i-1,j+1) - zs(i-1,j))*rdy + v(i-1,k,j)*(zs(i-1,j) - zs(i-1,j-1) )*rdy)   )   &
                   *(u(i,k,j)-u(i,k-1,j))/dsgmh(k)  
				               

	    enddo 
      ENDDO


		! 20220727 
		vortx = 0.0 
		
      j_start = jts
      j_end   = MIN(jte,jde-1)

	  

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end

    
    
      	vortx(i,k,j) =  (-1.0*vortx1(i,k,j) -  vortx21(i,k,j) - vortx22(i,k,j) )  * (msft(i,j)+msft(i+1,j))*0.5   ! m_x 

      ENDDO
      ENDDO
      ENDDO


	  
	  ! uds 对于开边界条件，不要计算边界值。202230117
	  
	  
 !  i_start = its
 !  i_end   = ite
   j_start = jts
   j_end   = min(jte,jde-1)

   IF ( config_flags%open_xs ) i_start = MAX(ids+1,its)
   IF ( config_flags%open_xe ) i_end   = MIN(ide-1,ite)
!   这里差分与y无关，open_y, j也是1从开始。
   IF ( config_flags%open_ys .or. specified ) i_start = MAX(ids+1,its)
   IF ( config_flags%open_ye .or. specified ) i_end   = MIN(ide-1,ite)
      IF ( config_flags%periodic_x ) i_start = its
      IF ( config_flags%periodic_x ) i_end = ite
	  
	  
	  
      DO  j = j_start, j_end
      DO k=2,ktf      ! 2021 12 21 
      DO  i = i_start, i_end

    
    ! du/dsgm
      	uds(i,k,j) = ( u(i,k,j) - u(i,k-1,j) )/dsgmh(k) 

      ENDDO
      ENDDO
      ENDDO	 
	  
	  
	 ! 上下边界 
	  DO  j = j_start, j_end
      DO  i = i_start, i_end

    ! du/dsgm
      	uds(i,kte,j) = uds(i,kte-1,j) 
        uds(i,kts,j) = uds(i,kts+1,j)
      ENDDO
      ENDDO
	  
	  
	  DO  j = j_start, j_end
      DO k=kts,kte  ! 202207
      DO  i = i_start, i_end

    
    
      	tmp(i,k,j) = 0.5*( w(i-1,k,j) + w(i,k,j) )  * uds(i,k,j) 

      ENDDO
      ENDDO
	  
	  
      ENDDO  
	  
 	  DO  j = j_start, j_end
      DO k=kts,ktf 
      DO  i = i_start, i_end

    
    
      	wudsh(i,k,j) =  (tmp(i,k,j) + tmp(i,k+1,j) ) *0.5 

      ENDDO
      ENDDO
      ENDDO    

	  
      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
 
      	vortx(i,k,j) = vortx(i,k,j)  -   (ztop/(ztop-zs(i,j))+ztop/(ztop-zs(i-1,j)))*0.5*wudsh(i,k,j)

      ENDDO
      ENDDO
      ENDDO	  




	  
      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
  
      !v(i-1,j+1)            v(i,j+1)
      !          u(i,j)
      !  
      !v(i-1,j)            v(i,j)

           vortx(i,k,j) = vortx(i,k,j)  - 0.25* ( v(i-1,k,j)*v(i-1,k,j)   + v(i,k,j)*v(i,k,j)  + v(i-1,k,j+1)*v(i-1,k,j+1)  + v(i,k,j+1)*v(i,k,j+1)  )  & 
                       *(msft(i,j)-msft(i-1,j))*rdx    &
                             + u(i,k,j) * ( v(i,k,j)*(msft(i,j)-msft(i,j-1))*rdy    &
                                                   +  v(i-1,k,j)*(msft(i-1,j)-msft(i-1,j-1))*rdy    &
                                                   +  v(i,k,j+1)*(msft(i,j+1)-msft(i,j))*rdy    &
                                                   +  v(i-1,k,j+1)*(msft(i-1,j+1)-msft(i-1,j))*rdy  ) * 0.25 
       
        
        !(ztop/(ztop-zs(i,j))+ztop/(ztop-zs(i-1,j)))*0.5*wudsh(i,k,j)

      ENDDO
      ENDDO
      ENDDO	  


if (.true.) then  ! template not add for supercell 20231229 ,20240111
      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
  

   !  vortx(i,k,j) = vortx(i,k,j)  +  u(i,k,j)*reradius*(0.25*( v(i-1,k,j)   + v(i,k,j)  + v(i-1,k,j+1)  + v(i,k,j+1)  )*tan(xlat(i,j)*degrad) &
    !                                    - 0.25*(w(i-1,k+1,j)+w(i-1,k,j)+w(i,k+1,j)+w(i,k,j)) )
       
        
    vortx(i,k,j) = vortx(i,k,j)  +  u(i,k,j)*reradius*(0.25*( v(i-1,k,j)   + v(i,k,j)  + v(i-1,k,j+1)  + v(i,k,j+1)  )*tan(xlat(i,j)*degrad) &
                                 - 0.25*(w(i-1,k+1,j)+w(i-1,k,j)+w(i,k+1,j)+w(i,k,j)) )

	
		

      ENDDO
      ENDDO
      ENDDO	   
 endif 
 


      END SUBROUTINE calcu_vortx	
















	  
	  
	!   calcu_gradx 为求x方向梯度的子程序，变量var可以是动能engk或气压p，或扰动气压p'。 grady,gradz,为y和z方向。
     
      SUBROUTINE calcu_gradx( var, gradx,                 &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   ! Input data 

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var    
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      !   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: mut
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: gradx

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   ! Local data 
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx            
       !  REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: zzh,dzzh
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: varz 

        REAL , DIMENSION( kms:kme ) :: sgmh,dsgm,dsgmh
		
		 REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: vardx 
		
          REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  ::  varbc
		
        LOGICAL :: specified

! 1) 判定 是否 给定边界
   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

 !2) 设定 do循环的 垂直 层次   
      ktf=MIN(kte,kde-1)

 	  DO k=kts,kde-1
          sgmh(k) = 0.5*(zh(k)+zh(k+1))
          dsgm(k) = zh(k+1) - zh(k) 
      ENDDO

 	  DO k=2,kde-1
             dsgmh(k) = 0.5*(dsgm(k) + dsgm(k-1))   ! new 20220711
			 
      ENDDO
	  
	  

         varbc  = var
	  
      CALL set_physical_bc3d( varbc, 'p', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              its, ite, jts, jte, kts, kte,       &
                              its, ite, jts, jte, kts, kte )



	 ! 3) 设定  x、 y do循环的 范围 
       
      i_start = its-1
      i_end   =  MIN(ite+1,ide-1)

      j_start = jts-1
      j_end   = MIN(jte+1,jde-1) 
	  

	  varz = 0.0 
	  
 DO  j = j_start, j_end	  
	 DO k=2,ktf
      DO  i = i_start, i_end
    
      	varz(i,k,j) =  (varbc(i,k,j) - varbc(i,k-1,j) ) / dsgmh(k)
      ENDDO
      ENDDO
 ! 上、下边界赋值,外插
      DO  i = i_start, i_end
       varz(i,kts,j) =  varz(i,kts+1,j)  !varz(i,kts+3,j)  -3.0*varz(i,kts+2,j) + 3.0*varz(i,kts+1,j) 
      varz(i,kte,j) =  varz(i,ktf,j)  ! varz(i,kte-3,j)  -3.0*varz(i,kte-2,j) + 3.0*varz(i,kte-1,j)  
      ENDDO
	  
	  
      ENDDO	  
  




      CALL set_physical_bc3d( varz, 'w', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              its, ite, jts, jte, kts, kte,       &
                              its, ite, jts, jte, kts, kte )
							  
							  

      CALL set_physical_bc2d( msft, 't', config_flags,  &
                              ids, ide, jds, jde,      &
                              ims, ime, jms, jme,      &
                              its, ite, jts, jte,      &
                              its, ite, jts, jte  )

      CALL set_physical_bc2d( zs, 't', config_flags,  &
                              ids, ide, jds, jde,      &
                              ims, ime, jms, jme,      &
                              its, ite, jts, jte,      &
                              its, ite, jts, jte  )	





	  
	gradx = 0.0 

	
	   ! for 气压方程的第二部分的梯度平均，用到外围值。
	!   j_start = MAX(jts-1,jds)
   !    j_end   = MIN(jte+1,jde-1)
      ! x-point
      i_start = MAX(its-1,ids)
      i_end   = MIN(ite+1, ide)
	  
	  
      j_start = jts
      j_end   = MIN(jte,jde-1)	

	! 取消测边界外插，因为，给定边界条件不需要外插，周期边条件也不需要。开边界，倾向计算计划用迎风方案或最终倾向量直接外插等方法，也不需要中间量的外插。  
	  
	  
! 4)  修订 开边界、给定边界、循环边界的x方向的 do循环 的范围 

      IF ( config_flags%open_xs ) i_start = MAX(ids+1,its)   ! 开边界的ids,ide, 需要特别处理
      IF ( config_flags%open_xe ) i_end   = MIN(ide-1,ite)
      IF ( specified ) i_start = MAX(ids+1,its-1)            ! 给定边界的，ids，ide，不做处理，ids+1, ide-1  二阶算法。
      IF ( specified ) i_end   = MIN(ide-1,ite+1)
	  
      IF ( config_flags%periodic_x ) i_start = its           ! 周期边界，因为设置了外围，可以计算全部its至ite范围。
      IF ( config_flags%periodic_x ) i_end = ite


!！  下面两行多余，南北开边界，不影响x方向上的梯度计算。
 !     IF ( config_flags%open_ys .or. specified ) j_start = MAX(jds+1,jts)
 !     IF ( config_flags%open_ye .or. specified ) j_end   = MIN(jde-1,jte)
!WRF
!  now the east-west (x) pressure gradient

 !  itf=ite
 !  jtf=MIN(jte,jde-1)
 !  ktf=MIN(kte,kde-1)
 !  i_start = its
 !  j_start = jts
 !  IF ( (config_flags%open_xs .or. specified .or. &
!           config_flags%nested ) .and. its == ids ) i_start = its+1
 !  IF ( (config_flags%open_xe .or. specified .or. &
!           config_flags%nested ) .and. ite == ide ) itf = itf-1
 !  IF ( config_flags%periodic_x ) i_start = its
 !  IF ( config_flags%periodic_x ) itf=ite
!WRF	  

		
  loop_j :     DO  j = j_start, j_end
	
		 
	 
      DO  k = kts, ktf
      DO  i = i_start, i_end 
	  
	    mrdx= 0.5* ( msft(i-1,j) + msft(i,j) ) *rdx
		
	    c1 = (msft(i,j)*(ztop-sgmh(k))/(ztop-zs(i,j)) + msft(i-1,j)*(ztop-sgmh(k))/(ztop-zs(i-1,j)))*0.5
	    c2 = rdx*(zs(i,j)-zs(i-1,j))*( (varz(i-1,k,j) + varz(i,k,j)) + (varz(i-1,k+1,j) + varz(i,k+1,j)) ) *0.25	  
	


    !   gradx(i,k,j) =  - mrdx*(varbc(i,k,j)-varbc(i-1,k,j))  - c1*c2
	! 2024 
	   gradx(i,k,j) =  - mrdx*(varbc(i,k,j)-varbc(i-1,k,j))  + c1*c2
	  

	  ENDDO
      ENDDO	 

      ENDDO    loop_j
	  


      END SUBROUTINE calcu_gradx



 
      SUBROUTINE calcu_grady( var0, grady,                 &
 !                       u, v, w,                       &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   ! Input data
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var0

        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags


         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT( OUT) :: grady
		 
		 
         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   ! Local data
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdy
         
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: varz  !,varzy,varzys
		   REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: var 

       REAL , DIMENSION( kms:kme ) :: sgmh,dsgm,dsgmh
    LOGICAL :: specified

   ! 1) 判定 是否 给定边界
   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

 !2) 设定 do循环的 垂直 层次   
      ktf=MIN(kte,kde-1)


 

	  DO k=kts,kde-1
          sgmh(k) = 0.5*(zh(k)+zh(k+1))
          dsgm(k) = zh(k+1) - zh(k) 
      ENDDO

 	  DO k=2,kde-1
             dsgmh(k) = 0.5*(dsgm(k) + dsgm(k-1))   ! new 20220711
			 
      ENDDO



	 ! 3) 设定  x、 y do循环的 范围 
       
      i_start = its-1
      i_end   =  MIN(ite+1,ide-1)

      j_start = jts-1
      j_end   = MIN(jte+1,jde-1)
	  

       var = var0
	  
      CALL set_physical_bc3d( var, 'p', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              its, ite, jts, jte, kts, kte,       &
                              its, ite, jts, jte, kts, kte )
							  
							  

      IF ( config_flags%periodic_y ) j_start = jts            ! 20220708
      IF ( config_flags%periodic_y ) j_end = jte 


       varz = 0.0 

      DO  j = j_start, j_end
      DO k=2,ktf
      DO  i = i_start, i_end
    
      	varz(i,k,j) =  (var(i,k,j) - var(i,k-1,j) ) / dsgmh(k)
      ENDDO
      ENDDO
      ENDDO
	  
	  

! 上、下边界赋值
      DO      j = j_start, j_end
  
      DO  i = i_start, i_end
		varz(i,kts,j) =  varz(i,kts+1,j)  !varz(i,kts+3,j)  -3.0*varz(i,kts+2,j) + 3.0*varz(i,kts+1,j) 
	   
		varz(i,kte,j) =  varz(i,kte-1,j)  ! varz(i,kte-3,j)  -3.0*varz(i,kte-2,j) + 3.0*varz(i,kte-1,j)  
      ENDDO
      ENDDO
	  
	

	  

   !   DO j = j_start, j_end
   !   DO k=kts ,ktf
    !  DO  i = i_start, i_end
 
	!	varzys(i,k,j) =  ( varz(i,k,j-1) + varz(i,k,j) + varz(i,k+1,j-1) + varz(i,k+1,j)  ) *0.25    !   (varz(i+1,k,j) + varz(i,k,j))*0.5  + (varz(i+1,k+1,j) + varz(i,k+1,j))*0.5
    !  ENDDO
     ! ENDDO
     ! ENDDO
	 
	 ! hanshenyou 20220801
	 CALL set_physical_bc3d( varz, 'w', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              its, ite, jts, jte, kts, kte,       &
                              its, ite, jts, jte, kts, kte )
							  
							  

      CALL set_physical_bc2d( msft, 't', config_flags,  &
                              ids, ide, jds, jde,      &
                              ims, ime, jms, jme,      &
                              its, ite, jts, jte,      &
                              its, ite, jts, jte  )

      CALL set_physical_bc2d( zs, 't', config_flags,  &
                              ids, ide, jds, jde,      &
                              ims, ime, jms, jme,      &
                              its, ite, jts, jte,      &
                              its, ite, jts, jte  )	
	  



 ! 4)  修订 开边界、给定边界、循环边界的x方向的 do循环 的范围 
     ! y梯度与x无关。
    !  IF ( config_flags%open_xs ) i_start = MAX(ids+1,its)   ! 开边界的ids,ide, 需要特别处理
    !  IF ( config_flags%open_xe ) i_end   = MIN(ide-1,ite) 
	
      i_start = its
      i_end   =  MIN(ite,ide-1)
   
 
	   
	   ! for 气压方程的第二部分的梯度平均，用到外围值。
	   j_start = MAX(jts-1,jds)
       j_end   = MIN(jte+1,jde-1)
	  
      IF ( config_flags%periodic_y ) j_start = jts          ! 周期边界，因为设置了外围，该子程序外再设置周期边界条件？。
      IF ( config_flags%periodic_y ) j_end = jte



    IF ( config_flags%open_ys .or. specified ) j_start = MAX(jds+1,jts-1)
    IF ( config_flags%open_ye .or. specified ) j_end   = MIN(jde-1,jte+1)
	
	
	  

      grady = 0.0 
  
      DO  j = j_start, j_end
      DO  k = kts, ktf
      DO  i = i_start, i_end 
	  
	  
	     mrdy= 0.5* ( msft(i,j-1) + msft(i,j) ) *rdy
	     c1 = (msft(i,j)*(ztop-sgmh(k))/(ztop-zs(i,j)) + msft(i,j-1)*(ztop-sgmh(k))/(ztop-zs(i,j-1)))*0.5 
	     c2 = rdy*(zs(i,j)-zs(i,j-1))*( varz(i,k,j-1) + varz(i,k,j) + varz(i,k+1,j-1) + varz(i,k+1,j)  ) *0.25
	   
	   
	!  grady(i,k,j) =       - mrdy*(var(i,k,j)-var(i,k,j-1)) - c1*c2
	! 2024   
	  grady(i,k,j) =       - mrdy*(var(i,k,j)-var(i,k,j-1)) + c1*c2
		 
	  ENDDO
      ENDDO
      ENDDO    
	  



      END SUBROUTINE calcu_grady


  
     SUBROUTINE calcu_gradz( var, gradz,                 &
 !                       u, v, w,                       &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   ! Input data
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var
	!	 ,u,     &
    !                                                                  v, &
    !                                                                  w
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      !   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: mut
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: gradz

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   ! Local data
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx
         
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: varz,varzy,varzys




     REAL , DIMENSION( kms:kme ) :: sgmh,dsgm,dsgmh
    LOGICAL :: specified


	specified = .false.
	if(config_flags%specified .or. config_flags%nested) specified = .true.
   
     ktf=MIN(kte,kde-1)
 
 ! sgmh,dsgm,dsgmh待改为参数传入到子程序。
	  DO k=kts,kde-1
          sgmh(k) = 0.5*(zh(k)+zh(k+1))
          dsgm(k) = zh(k+1) - zh(k) 
      ENDDO

 	  DO k=2,kde-1
             dsgmh(k) = 0.5*(dsgm(k) + dsgm(k-1))   ! new 20220711
			 
      ENDDO

    


      i_start = its
      i_end   = MIN(ide-1,ite)



      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
	  

      DO  j = j_start, j_end
      DO k=2,ktf
      DO  i = i_start, i_end
    
      	varz(i,k,j) =  (var(i,k,j) - var(i,k-1,j) ) / dsgmh(k)
      ENDDO
      ENDDO
      ENDDO

! 上、下边界赋值
! 20240121 ，梯度值的1阶外插， 适合气压全量，密度全量吗?
      DO      j = j_start, j_end
  
      DO  i = i_start, i_end
	  
       varz(i,kts,j) =  varz(i,kts+1,j)  !varz(i,kts+3,j)  -3.0*varz(i,kts+2,j) + 3.0*varz(i,kts+1,j) 
	   
      varz(i,kte,j) =  varz(i,kte-1,j)  ! varz(i,kte-3,j)  -3.0*varz(i,kte-2,j) + 3.0*varz(i,kte-1,j)  
       	  
      ENDDO
      ENDDO
	  
	  
	  
	  
! 因为，w、scalar的垂直平流 ，没有开边界和给定边界的影响，边界值必须计算垂直平流。

	  

	  gradz = 0.0
	  
      DO  j = j_start, j_end
      DO  k = kts, ktf
      DO  i = i_start, i_end 
			gradz(i,k,j) =  - ztop/(ztop-zs(i,j))*varz(i,k,j)
	  ENDDO
      ENDDO
      ENDDO    
	  
	   
	  

      END SUBROUTINE calcu_gradz
	  
	  



! 20220710
 SUBROUTINE average_vertical_w( varw, var,                 &
                     config_flags,                  &
                         fzm,fzp,                    &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   ! Input data
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var

        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: varw
		 
		 	  REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) ::  fzm,fzp
			  
			  


   ! Local data
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx
         
     !    REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: varz,varzy,varzys




  !   REAL , DIMENSION( kms:kme ) :: sgmh,dsgm,dsgmh
    LOGICAL :: specified


	specified = .false.
	if(config_flags%specified .or. config_flags%nested) specified = .true.
   
     ktf=MIN(kte,kde-1)
 
 


      i_start = its
      i_end   = MIN(ide-1,ite)



      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
	  

      varw  = 0. 
	  
	  
! 4)  修订 开边界、给定边界、循环边界的x方向的 do循环 的范围 

    !  IF ( config_flags%open_xs ) i_start = MAX(ids+1,its)   ! 开边界的ids,ide, 需要特别处理
    !  IF ( config_flags%open_xe ) i_end   = MIN(ide-1,ite)

    !  IF ( specified ) i_start = MAX(ids+1,its)            ! 给定边界的，ids，ide，不做处理，ids+1, ide-1  二阶算法。
    !  IF ( specified ) i_end   = MIN(ide-2,ite)
   

    !  IF ( config_flags%periodic_x ) i_start = its           ! 周期边界，因为设置了外围，可以计算全部its至ite范围。
    !  IF ( config_flags%periodic_x ) i_end = ite


    !  IF ( config_flags%open_ys .or. specified ) j_start = MAX(jds+1,jts)
    !  IF ( config_flags%open_ye .or. specified ) j_end   = MIN(jde-2,jte)	  
	  
	  

      DO  j = j_start, j_end
      DO k=kts+1,ktf
      DO  i = i_start, i_end
       !   fzm(k)*var(i,k,j)  + fzp(k)*var(i,k-1,j)   
      	varw(i,k,j) =  fzm(k)*var(i,k,j)  + fzp(k)*var(i,k-1,j)     
      ENDDO
      ENDDO
      ENDDO

! 上、下边界赋值
      DO      j = j_start, j_end
  
      DO  i = i_start, i_end
       varw(i,kts,j) =   2.0* var(i,kts,j)  -  varw(i,kts+1,j)  
	   
      varw(i,kte,j) =  2.0* var(i,kte-1,j)  -  varw(i,kte-1,j)   
      ENDDO
      ENDDO
	   
	  

      END SUBROUTINE average_vertical_w 



 SUBROUTINE average_vertical_pw( varw, var,pfm_base,rhofm_base,                 &
                     config_flags,                  &
                         fzm,fzp,                    &
						          dnw,  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   ! Input data
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var,pfm_base,rhofm_base

        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: varw
		 
		 	  REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) ::  fzm,fzp,dnw
			  
		real :: pi,temp,theta	  


   ! Local data
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx
         
     !    REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: varz,varzy,varzys




  !   REAL , DIMENSION( kms:kme ) :: sgmh,dsgm,dsgmh
    LOGICAL :: specified


	specified = .false.
	if(config_flags%specified .or. config_flags%nested) specified = .true.
   
     ktf=MIN(kte,kde-1)
 
 


      i_start = its
      i_end   = MIN(ide-1,ite)



      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
	  

      varw  = 0. 
	  
	  
! 4)  修订 开边界、给定边界、循环边界的x方向的 do循环 的范围 

    !  IF ( config_flags%open_xs ) i_start = MAX(ids+1,its)   ! 开边界的ids,ide, 需要特别处理
    !  IF ( config_flags%open_xe ) i_end   = MIN(ide-1,ite)

    !  IF ( specified ) i_start = MAX(ids+1,its)            ! 给定边界的，ids，ide，不做处理，ids+1, ide-1  二阶算法。
    !  IF ( specified ) i_end   = MIN(ide-2,ite)
   

    !  IF ( config_flags%periodic_x ) i_start = its           ! 周期边界，因为设置了外围，可以计算全部its至ite范围。
    !  IF ( config_flags%periodic_x ) i_end = ite


    !  IF ( config_flags%open_ys .or. specified ) j_start = MAX(jds+1,jts)
    !  IF ( config_flags%open_ye .or. specified ) j_end   = MIN(jde-2,jte)	  
	  
	  

      DO  j = j_start, j_end
      DO k=kts+1,ktf
      DO  i = i_start, i_end
       !   fzm(k)*var(i,k,j)  + fzp(k)*var(i,k-1,j)   
      	varw(i,k,j) =  exp(fzm(k)*log(var(i,k,j))  + fzp(k)*log(var(i,k-1,j)))     
      ENDDO
      ENDDO
      ENDDO

! 上、下边界赋值
      DO      j = j_start, j_end
  
      DO  i = i_start, i_end
     !  varw(i,kts,j) =   2.0* var(i,kts,j)  -  varw(i,kts+1,j)  
	   
      varw(i,kte,j) =  2.0* var(i,kte-1,j)  -  varw(i,kte-1,j)   
      ENDDO
      ENDDO
	   
      DO      j = j_start, j_end
  
      DO  i = i_start, i_end
	    temp =  pfm_base(i,kts,j)/rhofm_base(i,kts,j)/r_d
		pi    = (var(i,kts,j)/p0)**(r_d/Cp)
		theta = temp/pi-6.5/1000.0*0.5*dnw(kts) 
		varw(i,kts,j) = p0*(temp/theta)**(Cp/r_d) + var(i,kts,j) - pfm_base(i,kts,j)
		
   !    varw(i,kts,j) =   2.0* var(i,kts,j)  -  varw(i,kts+1,j)  
	   
     ! varw(i,kte,j) =  2.0* var(i,kte-1,j)  -  varw(i,kte-1,j)   
      ENDDO
      ENDDO	  
	 
 

      END SUBROUTINE average_vertical_pw
	  
	  
! 20220710
 SUBROUTINE average_vertical_pws( varw, var, rhofm_base,                 &
                     config_flags,                  &
                         fzm,fzp,                    &
						          dnw,  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   ! Input data
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var, rhofm_base

        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: varw
		 
		 	  REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) ::  fzm,fzp,dnw
			  
		real :: pi,temp,theta	  


   ! Local data
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx
         
     !    REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: varz,varzy,varzys




  !   REAL , DIMENSION( kms:kme ) :: sgmh,dsgm,dsgmh
    LOGICAL :: specified


	specified = .false.
	if(config_flags%specified .or. config_flags%nested) specified = .true.
   
     ktf=MIN(kte,kde-1)
 
 


      i_start = its
      i_end   = MIN(ide-1,ite)



      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
	  

      varw  = 0. 
	  
	  
! 4)  修订 开边界、给定边界、循环边界的x方向的 do循环 的范围 

    !  IF ( config_flags%open_xs ) i_start = MAX(ids+1,its)   ! 开边界的ids,ide, 需要特别处理
    !  IF ( config_flags%open_xe ) i_end   = MIN(ide-1,ite)

    !  IF ( specified ) i_start = MAX(ids+1,its)            ! 给定边界的，ids，ide，不做处理，ids+1, ide-1  二阶算法。
    !  IF ( specified ) i_end   = MIN(ide-2,ite)
   

    !  IF ( config_flags%periodic_x ) i_start = its           ! 周期边界，因为设置了外围，可以计算全部its至ite范围。
    !  IF ( config_flags%periodic_x ) i_end = ite


    !  IF ( config_flags%open_ys .or. specified ) j_start = MAX(jds+1,jts)
    !  IF ( config_flags%open_ye .or. specified ) j_end   = MIN(jde-2,jte)	  
	  
	  

      DO  j = j_start, j_end
      DO k=kts+1,ktf
      DO  i = i_start, i_end
       !   fzm(k)*var(i,k,j)  + fzp(k)*var(i,k-1,j)   
     ! 	varw(i,k,j) =  fzm(k)*var(i,k,j)  + fzp(k)*var(i,k-1,j)   
        varw(i,k,j) =  exp(fzm(k)*log(var(i,k,j))  + fzp(k)*log(var(i,k-1,j)))   		
      ENDDO
      ENDDO
      ENDDO

! 上、下边界赋值
      DO      j = j_start, j_end
  
      DO  i = i_start, i_end
    !   varw(i,kts,j) =   2.0* var(i,kts,j)  -  varw(i,kts+1,j)  
	   
      varw(i,kte,j) =  2.0* var(i,kte-1,j)  -  varw(i,kte-1,j)   
      ENDDO
      ENDDO
	   
      DO      j = j_start, j_end
  
      DO  i = i_start, i_end
	    temp =  var(i,kts,j)/rhofm_base(i,kts,j)/r_d
		pi    = (var(i,kts,j)/p0)**(r_d/Cp)
		theta = temp/pi-6.5/1000.0*0.5*dnw(kts) 
		varw(i,kts,j) = p0*(temp/theta)**(Cp/r_d)   
		
   !    varw(i,kts,j) =   2.0* var(i,kts,j)  -  varw(i,kts+1,j)  
	   
     ! varw(i,kte,j) =  2.0* var(i,kte-1,j)  -  varw(i,kte-1,j)   
      ENDDO
      ENDDO	  

      END SUBROUTINE average_vertical_pws







	  
SUBROUTINE average_vertical_w_x( varw, var,                 &
                     config_flags,                  &
                         fzm,fzp,                    &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   ! Input data
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var

        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: varw
		 
		 	  REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) ::  fzm,fzp
			  
			  


   ! Local data
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx
         
     !    REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: varz,varzy,varzys




  !   REAL , DIMENSION( kms:kme ) :: sgmh,dsgm,dsgmh
    LOGICAL :: specified


	specified = .false.
	if(config_flags%specified .or. config_flags%nested) specified = .true.
   
     ktf=MIN(kte,kde-1)
 
 


      i_start = its
      i_end   = MIN(ide,ite)



      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
	  

      varw  = 0. 
	  
	  
! 4)  修订 开边界、给定边界、循环边界的x方向的 do循环 的范围 
! 20240706 
!      IF ( config_flags%open_xs ) i_start = MAX(ids+1,its)   ! 开边界的ids,ide, 需要特别处理
 !     IF ( config_flags%open_xe ) i_end   = MIN(ide-1,ite)


! 20240706 
    !  IF ( specified ) i_start = MAX(ids+1,its)            ! 给定边界的，ids，ide，不做处理，ids+1, ide-1  二阶算法。
    !  IF ( specified ) i_end   = MIN(ide-2,ite)
   

      IF ( config_flags%periodic_x ) i_start = its           ! 周期边界，因为设置了外围，可以计算全部its至ite范围。
      IF ( config_flags%periodic_x ) i_end = ite

! 20240706 
 !     IF ( config_flags%open_ys .or. specified ) j_start = MAX(jds+1,jts)
 !     IF ( config_flags%open_ye .or. specified ) j_end   = MIN(jde-2,jte)	  
   	  
	  

      DO  j = j_start, j_end
      DO k=kts+1,ktf
      DO  i = i_start, i_end
       !   fzm(k)*var(i,k,j)  + fzp(k)*var(i,k-1,j)   
      	varw(i,k,j) =  fzm(k)*var(i,k,j)  + fzp(k)*var(i,k-1,j)     
      ENDDO
      ENDDO
      ENDDO

! 上、下边界赋值
      DO      j = j_start, j_end
  
      DO  i = i_start, i_end
       varw(i,kts,j) =   2.0* var(i,kts,j)  -  varw(i,kts+1,j)  
	   
      varw(i,kte,j) =  2.0* var(i,kte-1,j)  -  varw(i,kte-1,j)   
      ENDDO
      ENDDO
	   
	  

      END SUBROUTINE average_vertical_w_x 
	  
	  
	  
	  
	  
	  
SUBROUTINE advect_u_fm   ( tendency,             &
                        u, v, w,                   &
                        time_step, config_flags,  &
	                    msfux, msfuy, msfvx, msfvy,    &
                        msft,xlat,zs,    &
                        zh,fzm, fzp,                      &
                        rdx, rdy, ztop,                &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

   IMPLICIT NONE
  		

  
   ! Input data
   
   TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

   INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: u, &
                                                                      v, &
                                                                      w																	  


   
   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) :: tendency  


   REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msfux,  &
                                                                    msfuy,  &
                                                                    msfvx,  &
                                                                    msfvy
																	
																	
   REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft, &
                                                                    xlat, &
																	zs

   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: fzm,  &     !hanshenyou
                                                                  fzp,  &     ! hanshenyou  
                                                                  zh

   REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy,  &
																  ztop
																  
   INTEGER ,                                     INTENT(IN   ) :: time_step

   
  
   
   
   ! Local data
   
   INTEGER :: i, j, k, itf, jtf, ktf
   INTEGER :: i_start, i_end, j_start, j_end
   INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f

   
   LOGICAL :: degrade_xs, degrade_ys
   LOGICAL :: degrade_xe, degrade_ye

   LOGICAL :: specified
   
    REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) :: engk, engk_gradx,vortx    
	
	   REAL    :: mrdx, mrdy, ub, vb, uw, vw, dvm, dvp
	   INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip   
	real :: cbx 
	
	 LOGICAL :: fix_bdy 
	 

    IF ( config_flags%open_xs .and. config_flags%open_ys) THEN
              fix_bdy = .true.
			! 	fix_bdy = .false.
	ENDIf 
	
	
	
	cbx = 25.0 
!	zh,zs,fzm,fzp,

   specified = .false.

   if(config_flags%specified .or. config_flags%nested) specified = .true.





   ktf=MIN(kte,kde-1)


 


call    calcu_engk(engk,u,v,w,config_flags,                       &
                                  ids, ide, jds, jde, kds, kde,  &
                                  ims, ime, jms, jme, kms, kme,  &
                                  its, ite, jts, jte, kts, kte  )




   
call    calcu_gradx( engk, engk_gradx,                 &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )
						
   

call   calcu_vortx( vortx,                 &
                        u, v, w,                    &
                     config_flags,                  &
                        msft,xlat,zh,zs,fzm,fzp,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )
						
						

! 这里，边界条件在gradx和vortx中处理，开边界赋值为0，需要另外处理，次外层给定边界条件暂时忽略。
! engk_gradx，vortx ,open-y时，j=1时不为0. 
     tendency = tendency + 	 engk_gradx  +   vortx  

!  radiative lateral boundary condition in x for normal velocity (u)  

!-------
		
	 
		

      IF ( (config_flags%open_xs) .and. its == ids ) THEN

        j_start = jts
        j_end   = MIN(jte,jde-1)

        DO j = j_start, j_end
        DO k = kts, ktf
		
         ! openbdy  wrf  		
		
	     ub = MIN(u(its,k,j)-cbx, 0.)   ! bak 

		
    	   tendency(its,k,j) =     - rdx*ub*(u(its+1,k,j) - u(its,k,j))			  
 
		  
        ENDDO
        ENDDO

      ENDIF
	  
	  
	  

      IF ( (config_flags%open_xe) .and. ite == ide ) THEN

        j_start = jts
        j_end   = MIN(jte,jde-1)

        DO j = j_start, j_end
        DO k = kts, ktf
       	!wrf  
		  ub = MAX(u(ite,k,j)+cbx, 0.)
	 
          tendency(ite,k,j) =                                  &
                      - rdx*ub*(u(ite,k,j) - u(ite-1,k,j))
		
		  
        ENDDO
        ENDDO

      ENDIF
! WRF
!  pick up the rest of the horizontal radiation boundary conditions.
!  (these are the computations that don't require 'cb')
!  first, set to index ranges

      i_start = its
      i_end   = MIN(ite,ide)
      imin    = ids
      imax    = ide-1

      IF (config_flags%open_xs) THEN
        i_start = MAX(ids+1, its)
        imin = ids
      ENDIF
      IF (config_flags%open_xe) THEN
        i_end = MIN(ite,ide-1)
        imax = ide-1
      ENDIF

   IF( (config_flags%open_ys) .and. (jts == jds)) THEN

      DO i = i_start, i_end

         mrdy=msfux(i,jts)*rdy       ! ADT eqn 44, 2nd term on RHS
         ip = MIN( imax, i   )
         im = MAX( imin, i-1 )

         DO k=kts,ktf

    
		!wrf  
		  vw = 0.5*(v(ip,k,jts)+v(im,k,jts)) 
          vb = MIN( vw, 0. )
	  
		
          tendency(i,k,jts) =    tendency(i,k,jts) - mrdy*vb*(u(i,k,jts+1)-u(i,k,jts))   
		
	
						   
         ENDDO
      ENDDO

   ENDIF

   IF( (config_flags%open_ye) .and. (jte == jde)) THEN

      DO i = i_start, i_end

         mrdy=msfux(i,jte-1)*rdy     ! ADT eqn 44, 2nd term on RHS
         ip = MIN( imax, i   )
         im = MAX( imin, i-1 )

         DO k=kts,ktf

	 
	 ! wrf 
          vw = 0.5*(v(ip,k,jte)+v(im,k,jte))
          vb = MAX( vw, 0. )
	  
		  
         tendency(i,k,jte-1) =   tendency(i,k,jte-1)   -mrdy*vb*(u(i,k,jte-1)-u(i,k,jte-2))   
			 
         ENDDO
      ENDDO

   ENDIF
   
 


END SUBROUTINE advect_u_fm


SUBROUTINE advect_v_fm   ( tendency,            &
                        u, v, w,                   &
                        time_step, config_flags,  &
					   msfux, msfuy, msfvx, msfvy,    &
                        msft,msfty,xlat,zs,    &
                        zh,fzm, fzp,                      &
                        rdx, rdy, ztop,                &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

   IMPLICIT NONE
   
   ! Input data
   
    ! Input data
   
   TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

   INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: u, &
                                                                      v,   &
                                                                      w	


   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) :: tendency


   REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msfux,  &
                                                                    msfuy,  &
                                                                    msfvx,  &
                                                                    msfvy
																	
																	
   REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft, msfty, &
                                                                    xlat, &
																	zs

   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: fzm,  &     !hanshenyou
                                                                  fzp,  &     ! hanshenyou  
                                                                  zh

   REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy,  &
																  ztop
																  
   INTEGER ,                                     INTENT(IN   ) :: time_step

   

   ! Local data
       REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) :: engk, engk_grady, vorty   
	   
      REAL    :: mrdx, mrdy, ub, vb, uw, vw, dup, dum
	     INTEGER :: jmin, jmax, jp, jm, imin, imax
	  
   INTEGER :: i, j, k, itf, jtf, ktf
   INTEGER :: i_start, i_end, j_start, j_end
   INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f

   
   LOGICAL :: degrade_xs, degrade_ys
   LOGICAL :: degrade_xe, degrade_ye

   LOGICAL :: specified
   
   	real :: cbx 
	
	LOGICAL :: fix_bdy 
	 

    IF ( config_flags%open_xs .and. config_flags%open_ys) THEN
		fix_bdy = .true.
	  ! 	fix_bdy = .false.
	ENDIf 
	
	
	cbx = 25.0  

   specified = .false.

   if(config_flags%specified .or. config_flags%nested) specified = .true.

       ktf=MIN(kte,kde-1) 
	  




!  higher order flux has a 5 or 7 point stencil, so compute
!  bounds so we can switch to second order flux close to the boundary




call    calcu_engk(engk,u,v,w,config_flags,                       &
                                  ids, ide, jds, jde, kds, kde,  &
                                  ims, ime, jms, jme, kms, kme,  &
                                  its, ite, jts, jte, kts, kte  )



	
		!write(6,*)  'engk grady' 
	
call    calcu_grady( engk, engk_grady,                 &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )



!? 需要吗？ 20230117
	        CALL set_physical_bc3d( engk_grady, 'v', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              its, ite, jts, jte, kts, kte,       &
                              its, ite, jts, jte, kts, kte ) 
							  
							  
							  
						  

call   calcu_vorty( vorty,                 &
                        u, v, w,                    &
                     config_flags,                  &
                        msft,msfty,xlat,zh,zs,fzm,fzp,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )





! for 'v' 
      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = jte

      IF ( config_flags%open_ys ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye ) j_end   = MIN(jde-1,jte)
   ! IF ( specified ) j_start = MAX(jds+2,jts)
   !   IF ( specified ) j_end   = MIN(jde-2,jte)
      IF ( config_flags%polar ) j_start = MAX(jds+1,jts)
      IF ( config_flags%polar ) j_end   = MIN(jde-1,jte)

	  
	  
		
        DO j = j_start, j_end
        DO k = kts, ktf
	    DO i = i_start, i_end
		
			tendency(i,k,j)  = tendency(i,k,j)  + engk_grady(i,k,j) + vorty(i,k,j) 
	    ENDDO
        ENDDO
		ENDDO
		
	


		  
	


		
	if(fix_bdy) then 
	! fix bdy 
 
	else
	
	
   !  radiative lateral boundary condition in y for normal velocity (v)

      IF ( (config_flags%open_ys) .and. jts == jds ) THEN

        i_start = its
        i_end   = MIN(ite,ide-1)

        DO i = i_start, i_end
        DO k = kts, ktf
		
      !wrf
          vb = MIN(v(i,k,jts)-cbx, 0.)
		  
         tendency(i,k,jts) =                    &
                     - rdy*vb*(v(i,k,jts+1) - v(i,k,jts))
		
		  
					  
        ENDDO
        ENDDO

      ENDIF

      IF ( (config_flags%open_ye) .and. jte == jde ) THEN

        i_start = its
        i_end   = MIN(ite,ide-1)

        DO i = i_start, i_end
        DO k = kts, ktf
		
      !wrf
          vb = MAX(v(i,k,jte)+cbx, 0.)
		  
          tendency(i,k,jte) =                   &
                    - rdy*vb*(v(i,k,jte) - v(i,k,jte-1))
			  
        ENDDO
        ENDDO

      ENDIF

!  pick up the rest of the horizontal radiation boundary conditions.
!  (these are the computations that don't require 'cb'.
!  first, set to index ranges

      j_start = jts
      j_end   = MIN(jte,jde)

      jmin    = jds
      jmax    = jde-1

      IF (config_flags%open_ys) THEN
          j_start = MAX(jds+1, jts)
          jmin = jds
      ENDIF
      IF (config_flags%open_ye) THEN
          j_end = MIN(jte,jde-1)
          jmax = jde-1
      ENDIF

!  compute x (u) conditions for v, w, or scalar

   IF( (config_flags%open_xs) .and. (its == ids)) THEN

      DO j = j_start, j_end

         mrdx=msfvy(its,j)*rdx       ! ADT eqn 45, 1st term on RHS
         jp = MIN( jmax, j   )
         jm = MAX( jmin, j-1 )

         DO k=kts,ktf

  !  wrf 
          uw = 0.5*(u(its,k,jp)+u(its,k,jm))
          ub = MIN( uw, 0. )
       

	! 2023 ,不能仿照wrf，因为wrf为delt（uv），这里是udelta（v),其他类似。
	
		       !         tendency(its,k,j) =    -mrdx*ub*(v(its+1,k,j)-v(its,k,j))  
              tendency(its,k,j) =  tendency(its,k,j) -mrdx*ub*(v(its+1,k,j)-v(its,k,j))   

		   
         ENDDO
      ENDDO

   ENDIF

   IF( (config_flags%open_xe) .and. (ite == ide) ) THEN
      DO j = j_start, j_end

         mrdx=msfvy(ite-1,j)*rdx     ! ADT eqn 45, 1st term on RHS
         jp = MIN( jmax, j   )
         jm = MAX( jmin, j-1 )

         DO k=kts,ktf

   !wrf 
          uw = 0.5*(u(ite,k,jp)+u(ite,k,jm))
          ub = MAX( uw, 0. )
	   
		  

     !   tendency(ite-1,k,j) =    -mrdx* ub*(v(ite-1,k,j)-v(ite-2,k,j))   
        tendency(ite-1,k,j) =   tendency(ite-1,k,j)  -mrdx* ub*(v(ite-1,k,j)-v(ite-2,k,j))   

         ENDDO
      ENDDO

   ENDIF
   
  						
	endif  



	
		  
END SUBROUTINE advect_v_fm



SUBROUTINE advect_w_fm   ( tendency,   w_old,         &
                        u, v, w,                   &
                        time_step, config_flags,  &
                        msft,xlat,zs,    &
                        zh,fzm, fzp,                      &
                        rdx, rdy, ztop,                &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

   IMPLICIT NONE
   
   ! Input data
   
    ! Input data
   
   TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

   INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: u, &
                                                                      v, &
                                                                      w,w_old	


   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) :: tendency

   REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft, &
                                                                    xlat, &
																	zs

   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: fzm,  &     !hanshenyou
                                                                  fzp,  &     ! hanshenyou  
                                                                  zh

   REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy,  &
																  ztop
																  
   INTEGER ,                                     INTENT(IN   ) :: time_step



   ! Local data
   
   INTEGER :: i, j, k, itf, jtf, ktf
   INTEGER :: i_start, i_end, j_start, j_end, k_start, k_end 
!   INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
   
   
      INTEGER :: jmin, jmax, jp, jm, imin, imax
	  REAL    :: mrdx, mrdy, ub, vb, uw, vw
	  

       REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) :: engk, engk_gradz, vortz   
   
   !LOGICAL :: degrade_xs, degrade_ys
   !LOGICAL :: degrade_xe, degrade_ye

   LOGICAL :: specified

   specified = .false.

   if(config_flags%specified .or. config_flags%nested) specified = .true.









call    calcu_engk(engk,u,v,w,config_flags,                       &
                                  ids, ide, jds, jde, kds, kde,  &
                                  ims, ime, jms, jme, kms, kme,  &
                                  its, ite, jts, jte, kts, kte  )



call    calcu_gradz( engk, engk_gradz,                 &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )


call   calcu_vortz( vortz,                 &
                        u, v, w,                    &
                     config_flags,                  &
                        msft,xlat,zh,zs,fzm,fzp,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

	! tendency = tendency + engk_gradz + vortz


     ktf=MIN(kte,kde-1)


   
      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)


	  
	! 垂直平流，engk,engk_gradz,vortz已经计算了边界值，计算这些量时使用了边界的u，v，值。 20230608
	!  如果是周期边界条件，已经完成了计算
	 
	     k_start = kts
         k_end   = kte

	DO j = j_start, j_end
	
	
	DO k = k_start, k_end
     DO i = i_start, i_end
	     tendency(i,k,j)  = tendency(i,k,j)  + engk_gradz(i,k,j)  + vortz(i,k,j)  
 
 

     ENDDO
     ENDDO
	 
     ENDDO 
	 


! 重新计算开边界值。20230608
   
 !  pick up the the horizontal radiation boundary conditions.
!  (these are the computations that don't require 'cb'.
!  first, set to index ranges

if(.true.) then 

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)

   IF( (config_flags%open_xs) .and. (its == ids)) THEN

       DO j = j_start, j_end
       DO k = kts+1, ktf

         uw = 0.5*(fzm(k)*(u(its,k  ,j)+u(its+1,k  ,j)) +  &
                   fzp(k)*(u(its,k-1,j)+u(its+1,k-1,j))   )
         ub = MIN( uw, 0. )

         tendency(its,k,j) =                    &
               - rdx*(                                             &
                       ub*(w_old(its+1,k,j) - w_old(its,k,j)) +    &
                       w(its,k,j)*(                                &
                       fzm(k)*(u(its+1,k  ,j)-u(its,k  ,j))+     &
                       fzp(k)*(u(its+1,k-1,j)-u(its,k-1,j)))     &
                                                                  )
       ENDDO
       ENDDO

       k = ktf+1
       DO j = j_start, j_end

         uw = 0.5*( (2.-fzm(k-1))*(u(its,k-1,j)+u(its+1,k-1,j))   &
                   -fzp(k-1)*(u(its,k-2,j)+u(its+1,k-2,j))   )
         ub = MIN( uw, 0. )

         tendency(its,k,j) =                     &
               - rdx*(                                             &
                       ub*(w_old(its+1,k,j) - w_old(its,k,j)) +    &
                       w(its,k,j)*(                                &
                             (2.-fzm(k-1))*(u(its+1,k-1,j)-u(its,k-1,j))-  &
                             fzp(k-1)*(u(its+1,k-2,j)-u(its,k-2,j)))  &
                                                                  )
       ENDDO

   ENDIF

   IF( (config_flags%open_xe) .and. (ite == ide)) THEN

       DO j = j_start, j_end
       DO k = kts+1, ktf

         uw = 0.5*(fzm(k)*(u(ite-1,k  ,j)+u(ite,k  ,j)) +  &
                   fzp(k)*(u(ite-1,k-1,j)+u(ite,k-1,j))   )
         ub = MAX( uw, 0. )

         tendency(i_end,k,j) =                    &
               - rdx*(                                                 &
                       ub*(w_old(i_end,k,j) - w_old(i_end-1,k,j)) +    &
                       w(i_end,k,j)*(                                  &
                            fzm(k)*(u(ite,k  ,j)-u(ite-1,k  ,j)) +   &
                            fzp(k)*(u(ite,k-1,j)-u(ite-1,k-1,j)))    &
                                                                    )
       ENDDO
       ENDDO

       k = ktf+1
       DO j = j_start, j_end

         uw = 0.5*( (2.-fzm(k-1))*(u(ite-1,k-1,j)+u(ite,k-1,j))    &
                   -fzp(k-1)*(u(ite-1,k-2,j)+u(ite,k-2,j))   )
         ub = MAX( uw, 0. )

         tendency(i_end,k,j) =                      &
               - rdx*(                                                 &
                       ub*(w_old(i_end,k,j) - w_old(i_end-1,k,j)) +    &
                       w(i_end,k,j)*(                                  &
                               (2.-fzm(k-1))*(u(ite,k-1,j)-u(ite-1,k-1,j)) -   &
                               fzp(k-1)*(u(ite,k-2,j)-u(ite-1,k-2,j)))    &
                                                                    )
       ENDDO

   ENDIF


   IF( (config_flags%open_ys) .and. (jts == jds)) THEN

       DO i = i_start, i_end
       DO k = kts+1, ktf

         vw = 0.5*( fzm(k)*(v(i,k  ,jts)+v(i,k  ,jts+1)) +  &
                    fzp(k)*(v(i,k-1,jts)+v(i,k-1,jts+1))   )
         vb = MIN( vw, 0. )

         tendency(i,k,jts) =                     &
               - rdy*(                                             &
                       vb*(w_old(i,k,jts+1) - w_old(i,k,jts)) +    &
                       w(i,k,jts)*(                                &
                       fzm(k)*(v(i,k  ,jts+1)-v(i,k  ,jts))+     &
                       fzp(k)*(v(i,k-1,jts+1)-v(i,k-1,jts)))     &
                                                                )
       ENDDO
       ENDDO

       k = ktf+1
       DO i = i_start, i_end
         vw = 0.5*( (2.-fzm(k-1))*(v(i,k-1,jts)+v(i,k-1,jts+1))    &
                   -fzp(k-1)*(v(i,k-2,jts)+v(i,k-2,jts+1))   )
         vb = MIN( vw, 0. )

         tendency(i,k,jts) =                      &
               - rdy*(                                             &
                       vb*(w_old(i,k,jts+1) - w_old(i,k,jts)) +    &
                       w(i,k,jts)*(                                &
                          (2.-fzm(k-1))*(v(i,k-1,jts+1)-v(i,k-1,jts))-     &
                          fzp(k-1)*(v(i,k-2,jts+1)-v(i,k-2,jts)))     &
                                                                )
       ENDDO

   ENDIF

   IF( (config_flags%open_ye) .and. (jte == jde) ) THEN

       DO i = i_start, i_end
       DO k = kts+1, ktf

         vw = 0.5*( fzm(k)*(v(i,k  ,jte-1)+v(i,k  ,jte)) +  &
                    fzp(k)*(v(i,k-1,jte-1)+v(i,k-1,jte))   )
         vb = MAX( vw, 0. )

         tendency(i,k,j_end) =                    &
               - rdy*(                                                 &
                       vb*(w_old(i,k,j_end) - w_old(i,k,j_end-1)) +    &
                       w(i,k,j_end)*(                                  &
                            fzm(k)*(v(i,k  ,jte)-v(i,k  ,jte-1))+    &
                            fzp(k)*(v(i,k-1,jte)-v(i,k-1,jte-1)))    &
                                                                      )
       ENDDO
       ENDDO

       k = ktf+1
       DO i = i_start, i_end

         vw = 0.5*( (2.-fzm(k-1))*(v(i,k-1,jte-1)+v(i,k-1,jte))    &
                   -fzp(k-1)*(v(i,k-2,jte-1)+v(i,k-2,jte))   )
         vb = MAX( vw, 0. )

         tendency(i,k,j_end) =                     &
               - rdy*(                                                 &
                       vb*(w_old(i,k,j_end) - w_old(i,k,j_end-1)) +    &
                       w(i,k,j_end)*(                                  &
                               (2.-fzm(k-1))*(v(i,k-1,jte)-v(i,k-1,jte-1))-    &
                               fzp(k-1)*(v(i,k-2,jte)-v(i,k-2,jte-1)))    &
                                                                      )
       ENDDO

   ENDIF	

 endif 



END SUBROUTINE advect_w_fm	  
	  
	  
SUBROUTINE advect_scalar_pd_fm_old   ( field, field_old, tendency,    &
                                h_tendency, z_tendency,        &
                                u, v, sigmdot,                   &
                                c1, c2,                        &
                                mut, mub, mu_old,              &
                                time_step, config_flags,       &
                                tenddec,                       &
                                msfux, msfuy, msfvx, msfvy,    &
                                msftx, msfty,                  &
								msft,   &
                                fzm, fzp, ztop,zs,                       &
                                rdx, rdy, rdzw, dt,            &
                                ids, ide, jds, jde, kds, kde,  &
                                ims, ime, jms, jme, kms, kme,  &
                                its, ite, jts, jte, kts, kte  )

!  this is a first cut at a positive definite advection option
!  for scalars in WRF.  This version is memory intensive ->
!  we save 3d arrays of x, y and z both high and low order fluxes
!  (six in all).  Alternatively, we could sweep in a direction
!  and lower the cost considerably.

!  uses the Smolarkiewicz MWR 1989 approach, with addition of first-order
!  fluxes initially

!  WCS, 3 December 2002, 24 February 2003

   IMPLICIT NONE
   
   ! Input data
   
   TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

   LOGICAL ,                 INTENT(IN   ) :: tenddec  ! tendency flag

   INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: field,     &
                                                                      field_old, &
                                                                      u,    &
                                                                      v,    &
                                                                      sigmdot

   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: mut, mub, mu_old
   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) :: tendency
   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(  OUT) :: h_tendency, z_tendency

   REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msfux,  &
                                                                    msfuy,  &
                                                                    msfvx,  &
                                                                    msfvy,  &
                                                                    msftx,  &
                                                                    msfty,msft,zs  

   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: fzm,  &
                                                                  fzp,  &
                                                                  rdzw, &
                                                                  c1,   &
                                                                  c2

   REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy,  &
                                                                  dt,ztop 
   INTEGER ,                                     INTENT(IN   ) :: time_step

   ! Local data
   
   INTEGER :: i, j, k, itf, jtf, ktf
   INTEGER :: i_start, i_end, j_start, j_end
   INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
   INTEGER :: jmin, jmax, jp, jm, imin, imax

   REAL    :: mrdx, mrdy, ub, vb, uw, vw, mu

!  storage for high and low order fluxes

   REAL,  DIMENSION( its-1:ite+2, kts:kte, jts-1:jte+2  ) :: fqx, fqy, fqz

   REAL,  DIMENSION( its-1:ite+2, kts:kte, jts-1:jte+2  ) :: fqxl, fqyl, fqzl

   INTEGER :: horz_order, vert_order
   
   LOGICAL :: degrade_xs, degrade_ys
   LOGICAL :: degrade_xe, degrade_ye

   INTEGER :: jp1, jp0, jtmp

   REAL,DIMENSION( its-1:ite+2, kts:kte, jts-1:jte+2  ) :: flux_out, ph_low
   REAL :: scale
   !REAL :: flux_out, ph_low, scale
   REAL, PARAMETER :: eps=1.e-20


! definition of flux operators, 3rd, 4th, 5th or 6th order

   REAL    :: flux3, flux4, flux5, flux6, flux_upwind
   REAL    :: q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua, vel, cr

      flux4(q_im2, q_im1, q_i, q_ip1, ua) =                     &
            (7./12.)*(q_i + q_im1) - (1./12.)*(q_ip1 + q_im2)

      flux3(q_im2, q_im1, q_i, q_ip1, ua) =                     &
           flux4(q_im2, q_im1, q_i, q_ip1, ua) +                &
           sign(1,time_step)*sign(1.,ua)*(1./12.)*((q_ip1 - q_im2)-3.*(q_i-q_im1))

      flux6(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua) =       &
            (37./60.)*(q_i+q_im1) - (2./15.)*(q_ip1+q_im2)      &
            +(1./60.)*(q_ip2+q_im3)

      flux5(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua) =       &
           flux6(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua)    &
            -sign(1,time_step)*sign(1.,ua)*(1./60.)*(           &
              (q_ip2-q_im3)-5.*(q_ip1-q_im2)+10.*(q_i-q_im1) )

      flux_upwind(q_im1, q_i, cr ) = 0.5*min( 1.0,(cr+abs(cr)))*q_im1 &
                                    +0.5*max(-1.0,(cr-abs(cr)))*q_i

!      flux_upwind(q_im1, q_i, cr ) = 0.5*(1.+sign(1.,cr))*q_im1 &
!                                    +0.5*(1.-sign(1.,cr))*q_i
!      flux_upwind(q_im1, q_i, cr ) = 0.

    REAL     :: dx,dy,dz

    LOGICAL, PARAMETER :: pd_limit = .true.

! set order for the advection schemes

!  write(6,*) ' in pd advection routine '

    ! Empty arrays just in case:
    IF (config_flags%polar) THEN
       fqx(:,:,:)  = 0.
       fqy(:,:,:)  = 0.
       fqz(:,:,:)  = 0.
       fqxl(:,:,:) = 0.
       fqyl(:,:,:) = 0.
       fqzl(:,:,:) = 0.
    END IF

  ktf=MIN(kte,kde-1)
  horz_order = config_flags%h_sca_adv_order
  vert_order = config_flags%v_sca_adv_order

!  determine boundary mods for flux operators
!  We degrade the flux operators from 3rd/4th order
!   to second order one gridpoint in from the boundaries for
!   all boundary conditions except periodic and symmetry - these
!   conditions have boundary zone data fill for correct application
!   of the higher order flux stencils

   degrade_xs = .true.
   degrade_xe = .true.
   degrade_ys = .true.
   degrade_ye = .true.

!  begin with horizontal flux divergence
!  here is the choice of flux operators


  horizontal_order_test : IF( horz_order == 6 ) THEN

 

   ELSE IF( horz_order == 2 ) THEN

   IF( config_flags%periodic_x   .or. &
       config_flags%symmetric_xs .or. &
       (its > ids+1)                ) degrade_xs = .false.
   IF( config_flags%periodic_x   .or. &
       config_flags%symmetric_xe .or. &
       (ite < ide-2)                ) degrade_xe = .false.
   IF( config_flags%periodic_y   .or. &
       config_flags%symmetric_ys .or. &
       (jts > jds+1)                ) degrade_ys = .false.
   IF( config_flags%periodic_y   .or. &
       config_flags%symmetric_ye .or. &
       (jte < jde-2)                ) degrade_ye = .false.

!--  y flux compute; these bounds are for periodic and sym b.c.

      ktf=MIN(kte,kde-1)
      i_start = its-1
      i_end   = MIN(ite,ide-1)+1
      j_start = jts-1
      j_end   = MIN(jte,jde-1)+1

!--  modify loop bounds if open or specified

      IF(degrade_xs) i_start = its
      IF(degrade_xe) i_end   = MIN(ite,ide-1)
      IF(degrade_ys) j_start = MAX(jts,jds+1)
      IF(degrade_ye) j_end = MIN(jte,jde-2)

!  compute fluxes, 2nd order, y flux

      DO j = j_start, j_end+1
        DO k=kts,ktf
        DO i = i_start, i_end
 
		   
		   fqy(i,k,j) =  0.5*v(i,k,j)*          &
                  ((ztop-zs(i,j))*field(i,k,j)/msft(i,j) + (ztop-zs(i,j-1))*field(i,k,j-1)/msft(i,j-1))				  
				  
		   
        ENDDO
        ENDDO
      ENDDO

!  next, x flux



!   20230818 
	!  j_start = jts 
    ! j_end   = MIN(jte,jde-1) 
      IF(degrade_xs) i_start = its
      IF(degrade_xe) i_end   = MIN(ite,ide-1)
      IF(degrade_ys) j_start = jts  
      IF(degrade_ye) j_end = MIN(jte,jde-1)  
	  


      DO j = j_start, j_end
        DO k=kts,ktf
        DO i = i_start, i_end+1
          !  dx = 2./(msfty(i,j)+msfty(i-1,j))/rdx ! ADT eqn 48 d/dx
          !  mu = 0.5*((c1(k)*mut(i,j)+c2(k))+(c1(k)*mut(i-1,j)+c2(k)))
           vel =  u(i,k,j)
            cr = vel*dt/dx 
          !  fqxl(i,k,j) = mu*(dx/dt)*flux_upwind(field_old(i-1,k,j), field_old(i,k,j  ), cr)
		      fqxl(i,k,j) =  (dx/dt)*flux_upwind((ztop-zs(i-1,j))*field(i-1,k,j)/msft(i-1,j), (ztop-zs(i,j))*field(i,k,j)/msft(i,j), cr)
		  
          !  fqx( i,k,j ) = 0.5*ru(i,k,j)*          &
           !         (field(i,k,j)+field(i-1,k,j))
				  
            fqx( i,k,j ) = 0.5*u(i,k,j)*          &
                  ((ztop-zs(i,j))*field(i,k,j)/msft(i,j) + (ztop-zs(i-1,j))*field(i-1,k,j)/msft(i-1,j))				  
				  
				  

           fqx(i,k,j) = fqx(i,k,j) - fqxl(i,k,j)
        ENDDO
        ENDDO
      ENDDO

!--- end of 2nd order horizontal flux calculation

   ELSE

      WRITE ( wrf_err_message , * ) 'module_advect: advect_scalar_pd, h_order not known ',horz_order
      CALL wrf_error_fatal ( TRIM( wrf_err_message ) )

   ENDIF horizontal_order_test






















!  pick up the rest of the horizontal radiation boundary conditions.
!  (these are the computations that don't require 'cb'.
!  first, set to index ranges

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)

!  compute x (u) conditions for v, w, or scalar

   IF( (config_flags%open_xs) .and. (its == ids) ) THEN

       DO j = j_start, j_end
       DO k = kts, ktf
         ub = MIN( 0.5*(u(its,k,j)+u(its+1,k,j)), 0. )
         tendency(its,k,j) = tendency(its,k,j)                     &
             !  - rdx*(                                             &
             !          ub*(   field_old(its+1,k,j)                 &
             !               - field_old(its  ,k,j)   ) +           &
             !          field(its,k,j)*(ru(its+1,k,j)-ru(its,k,j))  &
             !                                                   )
                 - (ztop-zs(its,j))/msft(its,j)*rdx*(                                             &
                       ub*(   field_old(its+1,k,j)                 &
                            - field_old(its  ,k,j)   ) +           &
                       field_old(its,k,j)*(u(its+1,k,j)-u(its,k,j))  &       ! modified 
                                                                ) * msft(its,j)*msft(its,j)/(ztop-zs(its,j))  																	
																
																
       ENDDO
       ENDDO

   ENDIF

   IF( (config_flags%open_xe) .and. (ite == ide) ) THEN

       DO j = j_start, j_end
       DO k = kts, ktf
         ub = MAX( 0.5*(u(ite-1,k,j)+u(ite,k,j)), 0. )
         tendency(i_end,k,j) = tendency(i_end,k,j)                   &
        !       - rdx*(                                               &
         !              ub*(  field_old(i_end  ,k,j)                  &
         !                  - field_old(i_end-1,k,j) ) +              &
         !              field(i_end,k,j)*(ru(ite,k,j)-ru(ite-1,k,j))  &
         !                                                           )
																	
                - (ztop-zs(i_end,j))/msft(i_end,j)*rdx*(                                               &
                       ub*(  field_old(i_end  ,k,j)                  &
                           - field_old(i_end-1,k,j) ) +              &
                       field_old(i_end,k,j)*(u(ite,k,j)-u(ite-1,k,j))  &   ! modified  
                                                                    )* msft(i_end,j)*msft(i_end,j)/(ztop-zs(i_end,j)) 																		
																	
       ENDDO
       ENDDO

   ENDIF

   IF( (config_flags%open_ys) .and. (jts == jds) ) THEN

       DO i = i_start, i_end
       DO k = kts, ktf
         vb = MIN( 0.5*(v(i,k,jts)+v(i,k,jts+1)), 0. )
         tendency(i,k,jts) = tendency(i,k,jts)                     &
               - (ztop-zs(i,jts))/msft(i,jts)*rdy*(                                             &
                       vb*(  field_old(i,k,jts+1)                  &
                           - field_old(i,k,jts  ) ) +              &
                       field_old(i,k,jts)*(v(i,k,jts+1)-v(i,k,jts))  &     ! modified  
                                                               )* msft(i,jts)*msft(i,jts)/(ztop-zs(i,jts)) 


       ENDDO
       ENDDO

   ENDIF

   IF( (config_flags%open_ye) .and. (jte == jde)) THEN

       DO i = i_start, i_end
       DO k = kts, ktf
         vb = MAX( 0.5*(v(i,k,jte-1)+v(i,k,jte)), 0. )
         tendency(i,k,j_end) = tendency(i,k,j_end)                   &
               -  (ztop-zs(i,j_end))/msft(i,j_end)*rdy*(                                               &
                       vb*(   field_old(i,k,j_end  )                 &
                            - field_old(i,k,j_end-1) ) +             &
                       field_old(i,k,j_end)*(v(i,k,jte)-v(i,k,jte-1))  &    ! modified   
                                                                    )* msft(i,j_end)*msft(i,j_end)/(ztop-zs(i,j_end)) 
       ENDDO
       ENDDO

   ENDIF

 !  IF( (config_flags%polar) .and. (jts == jds) ) THEN

   

  ! ENDIF

  ! IF( (config_flags%polar) .and. (jte == jde)) THEN
 

  ! ENDIF




  ! write(6,*) 'advect_scalar_pd_fm maxval(abs(tendency) ) after openbdy  :', maxval(abs(tendency) ) 














!-------------------- vertical advection

!-- loop bounds for periodic or sym conditions

      i_start = its-1
      i_end   = MIN(ite,ide-1)+1
      j_start = jts-1
      j_end   = MIN(jte,jde-1)+1

!-- loop bounds for open or specified conditions

    IF(degrade_xs) i_start = MAX(its-1,ids)
    IF(degrade_xe) i_end   = MIN(ite+1,ide-1)
    IF(degrade_ys) j_start = MAX(jts-1,jds)
    IF(degrade_ye) j_end   = MIN(jte+1,jde-1)

    vert_order_test : IF (vert_order == 6) THEN    

 
   ELSE IF (vert_order == 2) THEN    

      DO j = j_start, j_end

         DO i = i_start, i_end
           fqz(i,1,j)  = 0.
           fqzl(i,1,j) = 0.
           fqz(i,kde,j)  = 0.
           fqzl(i,kde,j) = 0.
         ENDDO

         DO k=kts+1,ktf
         DO i = i_start, i_end

         !  dz = 2./(rdzw(k)+rdzw(k-1))
         !  mu = 0.5*((c1(k)*mut(i,j)+c2(k))+(c1(k)*mut(i,j)+c2(k)))
          ! vel = rom(i,k,j)
         !  cr = vel*dt/dz/mu
         !  fqzl(i,k,j) = mu*(dz/dt)*flux_upwind(field_old(i,k-1,j), field_old(i,k,j  ), cr)
        !   fqz(i,k,j)=rom(i,k,j)*(fzm(k)*field(i,k,j)+fzp(k)*field(i,k-1,j))
         !  fqz(i,k,j) = fqz(i,k,j) - fqzl(i,k,j)
		   
	   fqz(i,k,j) = sigmdot(i,k,j) *(fzm(k)*field(i,k,j)+fzp(k)*field(i,k-1,j))*(ztop-zs(i,j))/msft(i,j)    !这里不是m**2
		   

        ENDDO
        ENDDO

      ENDDO

   ELSE

      WRITE (wrf_err_message,*) ' advect_scalar_pd, v_order not known ',vert_order
      CALL wrf_error_fatal ( wrf_err_message )

   ENDIF vert_order_test
   
   
   
   
   
   
   
   
   
   
   
   
   
   

   IF (pd_limit) THEN

! positive definite filter

   i_start = its-1
   i_end   = MIN(ite,ide-1)+1
   j_start = jts-1
   j_end   = MIN(jte,jde-1)+1

!-- loop bounds for open or specified conditions

   IF(degrade_xs) i_start = MAX(its-1,ids)
   IF(degrade_xe) i_end   = MIN(ite+1,ide-1)
   IF(degrade_ys) j_start = MAX(jts-1,jds)
   IF(degrade_ye) j_end   = MIN(jte+1,jde-1)

   IF(config_flags%specified .or. config_flags%nested) THEN
     IF (degrade_xs) i_start = MAX(its-1,ids+1)
     IF (degrade_xe) i_end   = MIN(ite+1,ide-2)
     IF (degrade_ys) j_start = MAX(jts-1,jds+1)
     IF (degrade_ye) j_end   = MIN(jte+1,jde-2)
   END IF

   IF(config_flags%open_xs) THEN
     IF (degrade_xs) i_start = MAX(its-1,ids+1)
   END IF
   IF(config_flags%open_xe) THEN
     IF (degrade_xe) i_end   = MIN(ite+1,ide-2)
   END IF
   IF(config_flags%open_ys) THEN
     IF (degrade_ys) j_start = MAX(jts-1,jds+1)
   END IF
   IF(config_flags%open_ye) THEN
     IF (degrade_ye) j_end   = MIN(jte+1,jde-2)
   END IF
   ! ADT note:
   ! We don't want to change j_start and j_end
   ! for polar BC's since we want to calculate
   ! fluxes for directions other than y at the
   ! edge

!-- here is the limiter...

   DO j=j_start, j_end
   DO k=kts, ktf
#ifdef XEON_SIMD
!DIR$ simd
#else
!DIR$ vector always
#endif
   DO i=i_start, i_end

     ph_low(i,k,j) = field_old(i,k,j) &
                - dt*msft(i,j)*msft(i,j)/(ztop-zs(i,j)) *(               &
                       rdx*(fqx(i+1,k,j)-fqx(i,k,j)) +     &
                       rdy*(fqy(i,k,j+1)-fqy(i,k,j))     &
                      +rdzw(k)*(fqz(i,k+1,j)-fqz(i,k,j)) )  
                 
   ENDDO
   ENDDO
   ENDDO

   DO j=j_start, j_end
   DO k=kts, ktf
!DIR$ vector always
   DO i=i_start, i_end

     flux_out(i,k,j) = dt* msft(i,j)*msft(i,j)/(ztop-zs(i,j)) *( &
                                rdx*(  max(0.,fqx (i+1,k,j))      &
                                      -min(0.,fqx (i  ,k,j)) )    &
                               +rdy*(  max(0.,fqy (i,k,j+1))      &
                                      -min(0.,fqy (i,k,j  )) )   &
 !                              +rdzw(k)*(  min(0.,fqz (i,k+1,j))      &
 !                                     -max(0.,fqz (i,k  ,j)) )   )
                               +rdzw(k)*(  max(0.,fqz (i,k+1,j))      &
                                      -min(0.,fqz (i,k  ,j)) )   )									  
									  
									  
									  

   ENDDO
   ENDDO
   ENDDO 
   
      	

   DO j=j_start, j_end
   DO k=kts, ktf
!DIR$ vector always
   DO i=i_start, i_end
     IF( flux_out(i,k,j) .gt. ph_low(i,k,j) ) THEN
       scale = max(0.,ph_low(i,k,j)/(flux_out(i,k,j)+eps))
       IF( fqx (i+1,k,j) .gt. 0.) fqx(i+1,k,j) = scale*fqx(i+1,k,j)
       IF( fqx (i  ,k,j) .lt. 0.) fqx(i  ,k,j) = scale*fqx(i  ,k,j)
       IF( fqy (i,k,j+1) .gt. 0.) fqy(i,k,j+1) = scale*fqy(i,k,j+1)
       IF( fqy (i,k,j  ) .lt. 0.) fqy(i,k,j  ) = scale*fqy(i,k,j  )
!  note: z flux is opposite sign in mass coordinate because
!  vertical coordinate decreases with increasing k
       IF( fqz (i,k+1,j) .lt. 0.) fqz(i,k+1,j) = scale*fqz(i,k+1,j)
       IF( fqz (i,k  ,j) .gt. 0.) fqz(i,k  ,j) = scale*fqz(i,k  ,j)

     END IF

   ENDDO
   ENDDO
   ENDDO

   END IF
   
   


! add in the pd-limited flux divergence

  i_start = its
  i_end   = MIN(ite,ide-1)
  j_start = jts
  j_end   = MIN(jte,jde-1)

  DO j = j_start, j_end
  DO k = kts, ktf
!DEC$ vector always
  DO i = i_start, i_end

     tendency (i,k,j) = tendency(i,k,j)                           &
                            -rdzw(k)* msft(i,j)*msft(i,j)/(ztop-zs(i,j))*( fqz (i,k+1,j)-fqz (i,k,j)  &
                                   !   +fqzl(i,k+1,j)-fqzl(i,k,j)  &
									  )

  ENDDO
  ENDDO
  ENDDO
  
  !write(6,*) 'advect_scalar_pd_fm maxval(abs(tendency) ) fqz :', maxval(abs(tendency(ids:ide-1,kts:kte-1,jds:jde-1)) ) 
  
  IF(tenddec) THEN
  DO j = j_start, j_end
  DO k = kts, ktf
  DO i = i_start, i_end

     z_tendency (i,k,j) = 0. -rdzw(k)* msft(i,j)*msft(i,j)/(ztop-zs(i,j))*( fqz (i,k+1,j)-fqz (i,k,j)  &
                                    !  +fqzl(i,k+1,j)-fqzl(i,k,j)    &
									  )

  ENDDO
  ENDDO
  ENDDO
  END IF

! x flux divergence
!
  IF(degrade_xs) i_start = MAX(its,ids+1)
  IF(degrade_xe) i_end   = MIN(ite,ide-2)

  DO j = j_start, j_end
  DO k = kts, ktf
!DEC$ vector always  
  DO i = i_start, i_end

     ! Un-"canceled" map scale factor, ADT Eq. 48
     tendency (i,k,j) = tendency(i,k,j)                           &
               - msft(i,j)*msft(i,j)/(ztop-zs(i,j))*( rdx*( fqx (i+1,k,j)-fqx (i,k,j)     &
                               !    +fqxl(i+1,k,j)-fqxl(i,k,j)   &
								   ) )

  ENDDO
  ENDDO
  ENDDO
  
   ! write(6,*) 'advect_scalar_pd_fm maxval(abs(tendency) ) fqx :', maxval(abs(tendency(ids:ide-1,kts:kte-1,jds:jde-1)) ) 

  IF(tenddec) THEN
  DO j = j_start, j_end
  DO k = kts, ktf
  DO i = i_start, i_end

     h_tendency (i,k,j) = 0.                                      &
               - msft(i,j)*msft(i,j)/(ztop-zs(i,j))*( rdx*( fqx (i+1,k,j)-fqx (i,k,j)     &
                              !     +fqxl(i+1,k,j)-fqxl(i,k,j)     &
								   ))

  ENDDO
  ENDDO
  ENDDO
  END IF

! y flux divergence
!
  i_start = its
  i_end   = MIN(ite,ide-1)
  IF(degrade_ys) j_start = MAX(jts,jds+1)
  IF(degrade_ye) j_end   = MIN(jte,jde-2)

  DO j = j_start, j_end
  DO k = kts, ktf
!DEC$ vector always
  DO i = i_start, i_end

     ! Un-"canceled" map scale factor, ADT Eq. 48
     ! It is correct to use msftx (and not msfty), per W. Skamarock, 20080606
     tendency (i,k,j) = tendency(i,k,j)                           &
               -  msft(i,j)*msft(i,j)/(ztop-zs(i,j))*( rdy*( fqy (i,k,j+1)-fqy (i,k,j)     &
             !                      +fqyl(i,k,j+1)-fqyl(i,k,j)       &
								   ))

  ENDDO
  ENDDO
  ENDDO
  
  !write(6,*) 'advect_scalar_pd_fm maxval(abs(tendency) ) :', maxval(abs(tendency(ids:ide-1,kts:kte-1,jds:jde-1)) ) 
  

  IF(tenddec) THEN
  DO j = j_start, j_end
  DO k = kts, ktf
  DO i = i_start, i_end

     h_tendency (i,k,j) = h_tendency (i,k,j)                      &
               -  msft(i,j)*msft(i,j)/(ztop-zs(i,j))*( rdy*( fqy (i,k,j+1)-fqy (i,k,j)     &
                             !      +fqyl(i,k,j+1)-fqyl(i,k,j)    &
								  ) )

  ENDDO
  ENDDO
  ENDDO
  END IF

END SUBROUTINE advect_scalar_pd_fm_old



SUBROUTINE advect_scalar_fm2   ( field, field_old, tendency,    &
                             !   h_tendency, z_tendency,        &
                                u, v, sigmdot,                   &
                           !     c1, c2,                        &
                            !    mut, mub, mu_old,              &
                                time_step, config_flags,       &
								msft,   &
                                fzm, fzp, ztop,zs,                       &
                                rdx, rdy, rdzw, dt,            &
                                ids, ide, jds, jde, kds, kde,  &
                                ims, ime, jms, jme, kms, kme,  &
                                its, ite, jts, jte, kts, kte  )

!  this is a first cut at a positive definite advection option
!  for scalars in WRF.  This version is memory intensive ->
!  we save 3d arrays of x, y and z both high and low order fluxes
!  (six in all).  Alternatively, we could sweep in a direction
!  and lower the cost considerably.

!  uses the Smolarkiewicz MWR 1989 approach, with addition of first-order
!  fluxes initially

!  WCS, 3 December 2002, 24 February 2003

   IMPLICIT NONE
   
   ! Input data
   
   TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

 

   INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: field,     &
                                                                      field_old, &
                                                                      u,    &
                                                                      v,    &
                                                                      sigmdot

 
   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) :: tendency
 

   REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) ::    msft,zs  

   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: fzm,  &
                                                                  fzp,  &
                                                                  rdzw 

   REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy,  &
                                                                  dt,ztop 
   INTEGER ,                                     INTENT(IN   ) :: time_step

   ! Local data
   
   INTEGER :: i, j, k, itf, jtf, ktf
   INTEGER :: i_start, i_end, j_start, j_end
   INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
   INTEGER :: jmin, jmax, jp, jm, imin, imax

   REAL    :: mrdx, mrdy, ub, vb, uw, vw, mu

!  storage for high and low order fluxes

   REAL,  DIMENSION( its-1:ite+2, kts:kte, jts-1:jte+2  ) :: fqx, fqy, fqz

   REAL,  DIMENSION( its-1:ite+2, kts:kte, jts-1:jte+2  ) :: fqxl, fqyl, fqzl

   INTEGER :: horz_order, vert_order
   
   LOGICAL :: degrade_xs, degrade_ys
   LOGICAL :: degrade_xe, degrade_ye

   INTEGER :: jp1, jp0, jtmp

   REAL,DIMENSION( its-1:ite+2, kts:kte, jts-1:jte+2  ) :: flux_out, ph_low
   REAL :: scale
   !REAL :: flux_out, ph_low, scale
   REAL, PARAMETER :: eps=1.e-20


! definition of flux operators, 3rd, 4th, 5th or 6th order

   REAL    :: flux3, flux4, flux5, flux6, flux_upwind
   REAL    :: q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua, vel, cr

      flux4(q_im2, q_im1, q_i, q_ip1, ua) =                     &
            (7./12.)*(q_i + q_im1) - (1./12.)*(q_ip1 + q_im2)

      flux3(q_im2, q_im1, q_i, q_ip1, ua) =                     &
           flux4(q_im2, q_im1, q_i, q_ip1, ua) +                &
           sign(1,time_step)*sign(1.,ua)*(1./12.)*((q_ip1 - q_im2)-3.*(q_i-q_im1))

      flux6(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua) =       &
            (37./60.)*(q_i+q_im1) - (2./15.)*(q_ip1+q_im2)      &
            +(1./60.)*(q_ip2+q_im3)

      flux5(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua) =       &
           flux6(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua)    &
            -sign(1,time_step)*sign(1.,ua)*(1./60.)*(           &
              (q_ip2-q_im3)-5.*(q_ip1-q_im2)+10.*(q_i-q_im1) )

      flux_upwind(q_im1, q_i, cr ) = 0.5*min( 1.0,(cr+abs(cr)))*q_im1 &
                                    +0.5*max(-1.0,(cr-abs(cr)))*q_i

!      flux_upwind(q_im1, q_i, cr ) = 0.5*(1.+sign(1.,cr))*q_im1 &
!                                    +0.5*(1.-sign(1.,cr))*q_i
!      flux_upwind(q_im1, q_i, cr ) = 0.

    REAL     :: dx,dy,dz

    LOGICAL, PARAMETER :: pd_limit = .false.

! set order for the advection schemes

!  write(6,*) ' in pd advection routine '

    ! Empty arrays just in case:
    IF (config_flags%polar) THEN
       fqx(:,:,:)  = 0.
       fqy(:,:,:)  = 0.
       fqz(:,:,:)  = 0.
       fqxl(:,:,:) = 0.
       fqyl(:,:,:) = 0.
       fqzl(:,:,:) = 0.
    END IF

  ktf=MIN(kte,kde-1)
  horz_order = config_flags%h_sca_adv_order
  vert_order = config_flags%v_sca_adv_order

!  determine boundary mods for flux operators
!  We degrade the flux operators from 3rd/4th order
!   to second order one gridpoint in from the boundaries for
!   all boundary conditions except periodic and symmetry - these
!   conditions have boundary zone data fill for correct application
!   of the higher order flux stencils

   degrade_xs = .true.
   degrade_xe = .true.
   degrade_ys = .true.
   degrade_ye = .true.

!  begin with horizontal flux divergence
!  here is the choice of flux operators


  horizontal_order_test : IF( horz_order == 6 ) THEN

 

   ELSE IF( horz_order == 2 ) THEN

   IF( config_flags%periodic_x   .or. &
       config_flags%symmetric_xs .or. &
       (its > ids+1)                ) degrade_xs = .false.
   IF( config_flags%periodic_x   .or. &
       config_flags%symmetric_xe .or. &
       (ite < ide-2)                ) degrade_xe = .false.
   IF( config_flags%periodic_y   .or. &
       config_flags%symmetric_ys .or. &
       (jts > jds+1)                ) degrade_ys = .false.
   IF( config_flags%periodic_y   .or. &
       config_flags%symmetric_ye .or. &
       (jte < jde-2)                ) degrade_ye = .false.

!--  y flux compute; these bounds are for periodic and sym b.c.

      ktf=MIN(kte,kde-1)
      i_start = its-1
      i_end   = MIN(ite,ide-1)+1
      j_start = jts-1
      j_end   = MIN(jte,jde-1)+1

!--  modify loop bounds if open or specified

      IF(degrade_xs) i_start = its
      IF(degrade_xe) i_end   = MIN(ite,ide-1)
      IF(degrade_ys) j_start = MAX(jts,jds+1)
      IF(degrade_ye) j_end = MIN(jte,jde-2)

!  compute fluxes, 2nd order, y flux

      DO j = j_start, j_end+1
        DO k=kts,ktf
        DO i = i_start, i_end
 
		   
		   fqy(i,k,j) =  0.5*v(i,k,j)*          &
                  ((ztop-zs(i,j))*field(i,k,j)/msft(i,j) + (ztop-zs(i,j-1))*field(i,k,j-1)/msft(i,j-1))				  
				  
		   
        ENDDO
        ENDDO
      ENDDO

!  next, x flux

      DO j = j_start, j_end
        DO k=kts,ktf
        DO i = i_start, i_end+1
 
				  
            fqx( i,k,j ) = 0.5*u(i,k,j)*          &
                  ((ztop-zs(i,j))*field(i,k,j)/msft(i,j) + (ztop-zs(i-1,j))*field(i-1,k,j)/msft(i-1,j))				  
				  
				  
 
        ENDDO
        ENDDO
      ENDDO

!--- end of 2nd order horizontal flux calculation

   ELSE

      WRITE ( wrf_err_message , * ) 'module_advect: advect_scalar_pd, h_order not known ',horz_order
      CALL wrf_error_fatal ( TRIM( wrf_err_message ) )

   ENDIF horizontal_order_test






















!  pick up the rest of the horizontal radiation boundary conditions.
!  (these are the computations that don't require 'cb'.
!  first, set to index ranges

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)

!  compute x (u) conditions for v, w, or scalar

   IF( (config_flags%open_xs) .and. (its == ids) ) THEN

       DO j = j_start, j_end
       DO k = kts, ktf
         ub = MIN( 0.5*(u(its,k,j)+u(its+1,k,j)), 0. )
         tendency(its,k,j) = tendency(its,k,j)                     &
             !  - rdx*(                                             &
             !          ub*(   field_old(its+1,k,j)                 &
             !               - field_old(its  ,k,j)   ) +           &
             !          field(its,k,j)*(ru(its+1,k,j)-ru(its,k,j))  &
             !                                                   )
                 - (ztop-zs(its,j))/msft(its,j)*rdx*(                                             &
                       ub*(   field_old(its+1,k,j)                 &
                            - field_old(its  ,k,j)   ) +           &
                       field_old(its,k,j)*(u(its+1,k,j)-u(its,k,j))  &       ! modified 
                                                                ) * msft(its,j)*msft(its,j)/(ztop-zs(its,j))  																	
																
																
       ENDDO
       ENDDO

   ENDIF

   IF( (config_flags%open_xe) .and. (ite == ide) ) THEN

       DO j = j_start, j_end
       DO k = kts, ktf
         ub = MAX( 0.5*(u(ite-1,k,j)+u(ite,k,j)), 0. )
         tendency(i_end,k,j) = tendency(i_end,k,j)                   &
        !       - rdx*(                                               &
         !              ub*(  field_old(i_end  ,k,j)                  &
         !                  - field_old(i_end-1,k,j) ) +              &
         !              field(i_end,k,j)*(ru(ite,k,j)-ru(ite-1,k,j))  &
         !                                                           )
																	
                - (ztop-zs(i_end,j))/msft(i_end,j)*rdx*(                                               &
                       ub*(  field_old(i_end  ,k,j)                  &
                           - field_old(i_end-1,k,j) ) +              &
                       field_old(i_end,k,j)*(u(ite,k,j)-u(ite-1,k,j))  &   ! modified  
                                                                    )* msft(i_end,j)*msft(i_end,j)/(ztop-zs(i_end,j)) 																		
																	
       ENDDO
       ENDDO

   ENDIF

   IF( (config_flags%open_ys) .and. (jts == jds) ) THEN

       DO i = i_start, i_end
       DO k = kts, ktf
         vb = MIN( 0.5*(v(i,k,jts)+v(i,k,jts+1)), 0. )
         tendency(i,k,jts) = tendency(i,k,jts)                     &
               - (ztop-zs(i,jts))/msft(i,jts)*rdy*(                                             &
                       vb*(  field_old(i,k,jts+1)                  &
                           - field_old(i,k,jts  ) ) +              &
                       field_old(i,k,jts)*(v(i,k,jts+1)-v(i,k,jts))  &     ! modified  
                                                               )* msft(i,jts)*msft(i,jts)/(ztop-zs(i,jts)) 


       ENDDO
       ENDDO

   ENDIF

   IF( (config_flags%open_ye) .and. (jte == jde)) THEN

       DO i = i_start, i_end
       DO k = kts, ktf
         vb = MAX( 0.5*(v(i,k,jte-1)+v(i,k,jte)), 0. )
         tendency(i,k,j_end) = tendency(i,k,j_end)                   &
               -  (ztop-zs(i,j_end))/msft(i,j_end)*rdy*(                                               &
                       vb*(   field_old(i,k,j_end  )                 &
                            - field_old(i,k,j_end-1) ) +             &
                       field_old(i,k,j_end)*(v(i,k,jte)-v(i,k,jte-1))  &    ! modified   
                                                                    )* msft(i,j_end)*msft(i,j_end)/(ztop-zs(i,j_end)) 
       ENDDO
       ENDDO

   ENDIF

 !  IF( (config_flags%polar) .and. (jts == jds) ) THEN

   

  ! ENDIF

  ! IF( (config_flags%polar) .and. (jte == jde)) THEN
 

  ! ENDIF




  ! write(6,*) 'advect_scalar_pd_fm maxval(abs(tendency) ) after openbdy  :', maxval(abs(tendency) ) 














!-------------------- vertical advection

!-- loop bounds for periodic or sym conditions

      i_start = its-1
      i_end   = MIN(ite,ide-1)+1
      j_start = jts-1
      j_end   = MIN(jte,jde-1)+1

!-- loop bounds for open or specified conditions

    IF(degrade_xs) i_start = MAX(its-1,ids)
    IF(degrade_xe) i_end   = MIN(ite+1,ide-1)
    IF(degrade_ys) j_start = MAX(jts-1,jds)
    IF(degrade_ye) j_end   = MIN(jte+1,jde-1)

    vert_order_test : IF (vert_order == 6) THEN    

 
   ELSE IF (vert_order == 2) THEN    

      DO j = j_start, j_end

         DO i = i_start, i_end
           fqz(i,1,j)  = 0.
           fqzl(i,1,j) = 0.
           fqz(i,kde,j)  = 0.
           fqzl(i,kde,j) = 0.
         ENDDO

         DO k=kts+1,ktf
         DO i = i_start, i_end

         !  dz = 2./(rdzw(k)+rdzw(k-1))
         !  mu = 0.5*((c1(k)*mut(i,j)+c2(k))+(c1(k)*mut(i,j)+c2(k)))
          ! vel = rom(i,k,j)
         !  cr = vel*dt/dz/mu
         !  fqzl(i,k,j) = mu*(dz/dt)*flux_upwind(field_old(i,k-1,j), field_old(i,k,j  ), cr)
        !   fqz(i,k,j)=rom(i,k,j)*(fzm(k)*field(i,k,j)+fzp(k)*field(i,k-1,j))
         !  fqz(i,k,j) = fqz(i,k,j) - fqzl(i,k,j)
		   
	   fqz(i,k,j) = sigmdot(i,k,j) *(fzm(k)*field(i,k,j)+fzp(k)*field(i,k-1,j))*(ztop-zs(i,j))/msft(i,j)    !这里不是m**2
		   

        ENDDO
        ENDDO

      ENDDO

   ELSE

      WRITE (wrf_err_message,*) ' advect_scalar_pd, v_order not known ',vert_order
      CALL wrf_error_fatal ( wrf_err_message )

   ENDIF vert_order_test
   
   
   


  i_start = its
  i_end   = MIN(ite,ide-1)
  j_start = jts
  j_end   = MIN(jte,jde-1)

  DO j = j_start, j_end
  DO k = kts, ktf
!DEC$ vector always
  DO i = i_start, i_end

     tendency (i,k,j) = tendency(i,k,j)                           &
                            -rdzw(k)* msft(i,j)*msft(i,j)/(ztop-zs(i,j))*( fqz (i,k+1,j)-fqz (i,k,j)  &
                                   !   +fqzl(i,k+1,j)-fqzl(i,k,j)  &
									  )

  ENDDO
  ENDDO
  ENDDO
  
  !write(6,*) 'advect_scalar_pd_fm maxval(abs(tendency) ) fqz :', maxval(abs(tendency(ids:ide-1,kts:kte-1,jds:jde-1)) ) 
  
 

! x flux divergence
!
  IF(degrade_xs) i_start = MAX(its,ids+1)
  IF(degrade_xe) i_end   = MIN(ite,ide-2)

  DO j = j_start, j_end
  DO k = kts, ktf
!DEC$ vector always  
  DO i = i_start, i_end

     ! Un-"canceled" map scale factor, ADT Eq. 48
     tendency (i,k,j) = tendency(i,k,j)                           &
               - msft(i,j)*msft(i,j)/(ztop-zs(i,j))*( rdx*( fqx (i+1,k,j)-fqx (i,k,j)     &
                               !    +fqxl(i+1,k,j)-fqxl(i,k,j)   &
								   ) )

  ENDDO
  ENDDO
  ENDDO
  
 

! y flux divergence
!
  i_start = its
  i_end   = MIN(ite,ide-1)
  IF(degrade_ys) j_start = MAX(jts,jds+1)
  IF(degrade_ye) j_end   = MIN(jte,jde-2)

  DO j = j_start, j_end
  DO k = kts, ktf
!DEC$ vector always
  DO i = i_start, i_end

     ! Un-"canceled" map scale factor, ADT Eq. 48
     ! It is correct to use msftx (and not msfty), per W. Skamarock, 20080606
     tendency (i,k,j) = tendency(i,k,j)                           &
               -  msft(i,j)*msft(i,j)/(ztop-zs(i,j))*( rdy*( fqy (i,k,j+1)-fqy (i,k,j)     &
             !                      +fqyl(i,k,j+1)-fqyl(i,k,j)       &
								   ))

  ENDDO
  ENDDO
  ENDDO
  
  !write(6,*) 'advect_scalar_pd_fm maxval(abs(tendency) ) :', maxval(abs(tendency(ids:ide-1,kts:kte-1,jds:jde-1)) ) 
  

 

END SUBROUTINE advect_scalar_fm2
!----------------------------------------------------------------

! 改变流程，先边界，后。

SUBROUTINE calcu_divg2( var, divg,                 &
                        u, v, sgmdot,                       &
                     config_flags,                  &
                        msft,zh,zs,fzm,fzp,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   ! Input data
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var,u,     &
                                                                      v, &
                                                                      sgmdot
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      !   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: mut
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: divg

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh,fzm,fzp


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   ! Local data
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx
         
        
       !  REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: zzh,dzzh
	   
	   
	   
	   
	   
		 
		    REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) ::  varw
		 
		 
		 
		 


      REAL , DIMENSION( kms:kme ) :: sgmh,dsgm
       REAL , DIMENSION( ims:ime , jms:jme )  ::  rmsft,rmm
	   
	   REAL :: ub,vb 
				 
				 

    LOGICAL :: specified

   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

   
     ktf=MIN(kte,kde-1)
 

	  DO k=kts,kde-1
          sgmh(k) = 0.5*(zh(k)+zh(k+1))
          dsgm(k) = zh(k+1) - zh(k) 
      ENDDO


 



	  

   ! 20230228  for rho couple scalar 
      CALL set_physical_bc3d( var, 'p', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              its, ite, jts, jte, kts, kte,       &
                              its, ite, jts, jte, kts, kte )
	  
	  

     divg = 0.0 
	 
	 
! boundary	 
	 
!  ---- x direction -----------
 

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)


      IF ( .NOT. config_flags%periodic_x ) THEN
        IF ( config_flags%open_xs .or. specified ) i_start = MAX(ids+1,its)
        IF ( config_flags%open_xe .or. specified ) i_end   = MIN(ide-2,ite)
      ENDIF

      
	  
	  
	IF ( config_flags%open_xs  .and.  its == ids ) THEN 
		
	!	i = ids  !  用于方法1 
		
      DO  j = j_start, j_end
      DO k=kts,ktf
   
 					 
	! 参考wrf	,这一方法的使用，改进非常大。避免了hill2km_x试验的西边界问题。			 
	         ub = MIN( 0.5*(u(its,k,j)+u(its+1,k,j)), 0. )
		 
         divg(its,k,j) =                      &
               - (ztop-zs(its,j))/msft(its,j)*rdx*(                                             &
                       ub*(   var(its+1,k,j)                 &
                            - var(its  ,k,j)   ) +           &
                       var(its,k,j)*(u(its+1,k,j)-u(its,k,j))  &       ! modified 
                                                                )				 
					 
      ENDDO
      ENDDO
 
	  
	ENDIF	  
	

	  
	    IF ( config_flags%open_xe  .and.  ite == ide ) THEN 
		
	 !	i = ite-1
		
		 i_end   = MIN(ite,ide-1)
		 
		
      DO  j = j_start, j_end
      DO k=kts,ktf
 
  
     
	     ub = MAX( 0.5*(u(ite-1,k,j)+u(ite,k,j)), 0. )	 
		 
         divg(i_end,k,j) =                  &
               - (ztop-zs(i_end,j))/msft(i_end,j)*rdx*(                                               &
                       ub*(  var(i_end  ,k,j)                  &
                           - var(i_end-1,k,j) ) +              &
                       var(i_end,k,j)*(u(ite,k,j)-u(ite-1,k,j))  &   ! modified  
                                                                    )				 
					 
					 
      ENDDO
      ENDDO
  

	ENDIF	
	  
	  
	  
	  

!  ---- y direction -----------

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
      ! Polar boundary conditions are like open or specified
      IF ( config_flags%open_ys .or. specified .or. config_flags%polar ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified .or. config_flags%polar ) j_end   = MIN(jde-2,jte) 	  

       
	  
	    IF ( config_flags%open_ys  .and.  jts == jds ) THEN 
		
	!	j = jds 
		
      DO k=kts,ktf
      DO  i = i_start, i_end		

	
	      vb = MIN( 0.5*(v(i,k,jts)+v(i,k,jts+1)), 0. )
		  
        divg(i,k,jts) =    divg(i,k,jts)                &
               - (ztop-zs(i,jts))/msft(i,jts)*rdy*(                                             &
                       vb*(  var(i,k,jts+1)                  &
                           - var(i,k,jts  ) ) +              &
                       var(i,k,jts)*(v(i,k,jts+1)-v(i,k,jts))  &     ! modified  
                                                               )
					 
      ENDDO
      ENDDO
	    ENDIF	
		
 
	  
	    IF ( config_flags%open_ye  .and.  jte == jde ) THEN 
		
	 ! 	j = jte-1
		
		j_end   = MIN(jte,jde-1)
		
      DO k=kts,ktf
      DO  i = i_start, i_end		
		
	 
	
		  vb = MAX( 0.5*(v(i,k,jte-1)+v(i,k,jte)), 0. ) 
		 
         divg(i,k,j_end) =    divg(i,k,j_end)               &
               -  (ztop-zs(i,j_end))/msft(i,j_end)*rdy*(                                               &
                       vb*(   var(i,k,j_end  )                 &
                            - var(i,k,j_end-1) ) +             &
                       var(i,k,j_end)*(v(i,k,jte)-v(i,k,jte-1))  &    ! modified   
                                                                    )
					 
      ENDDO
      ENDDO
	    ENDIF
		
		








!!!!!!!!!!!!!!!!!!!!!!!!!!!	 
	 
	 
	 
	 
!  ---- z  direction -----------	  
	  
	  	call average_vertical_w( varw, var,                 &
                     config_flags,                  &
                         fzm,fzp,                    &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )	



  !  计划拆分为简单公式，处理变量从半层插值到整层，及边界处理问题
  
  
        i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
	  
      DO  j = j_start, j_end
      DO  i = i_start, i_end
        rmsft(i,j)  = 1./msft(i,j) 
        rmm(i,j) = 1./(msft(i,j)*msft(i,j) )
      ENDDO
      ENDDO             

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end



            divg(i,k,j) = divg(i,k,j) -    ( varw(i,k+1,j)* rmm(i,j) * sgmdot(i,k+1,j) - varw(i,k,j)* rmm(i,j) * sgmdot(i,k,j)   ) /dsgm(k) * (ztop-zs(i,j))


      ENDDO
      ENDDO
      ENDDO
	 
	 
	 
	 
	 
	 
	 
!  ---- x direction -----------
 

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)


      IF ( .NOT. config_flags%periodic_x ) THEN
        IF ( config_flags%open_xs .or. specified ) i_start = MAX(ids+1,its)
        IF ( config_flags%open_xe .or. specified ) i_end   = MIN(ide-2,ite)
      ENDIF

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
  
    

      	divg(i,k,j) =  divg(i,k,j) -  (  0.5*((ztop-zs(i,j))*var(i,k,j)/msft(i,j)  + (ztop-zs(i+1,j))*var(i+1,k,j)/msft(i+1,j) ) *u(i+1,k,j) &
                     -   0.5*((ztop-zs(i-1,j))*var(i-1,k,j)/msft(i-1,j)  + (ztop-zs(i,j))*var(i,k,j)/msft(i,j) ) *u(i,k,j) )*rdx
      ENDDO
      ENDDO
      ENDDO
	  
	

!  ---- y direction -----------

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
      ! Polar boundary conditions are like open or specified
      IF ( config_flags%open_ys .or. specified .or. config_flags%polar ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified .or. config_flags%polar ) j_end   = MIN(jde-2,jte) 	  

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end

   !   	divg(i,k,j) = divg(i,k,j) -   ( 0.5*((ztop-zs(i,j))*var(i,k,j)*rmsft(i,j)  + (ztop-zs(i,j+1))*var(i,k,j+1)*rmsft(i,j+1) ) *v(i,k,j+1)    &
   !                  -   0.5*((ztop-zs(i,j-1))*var(i,k,j-1)*rmsft(i,j-1)  + (ztop-zs(i,j))*var(i,k,j)*rmsft(i,j) ) *v(i,k,j) )*rdy
					 
 		 
    	divg(i,k,j) = divg(i,k,j) -   ( 0.5*((ztop-zs(i,j))*var(i,k,j)/msft(i,j)  + (ztop-zs(i,j+1))*var(i,k,j+1)/msft(i,j+1) ) *v(i,k,j+1)    &
                 -   0.5*((ztop-zs(i,j-1))*var(i,k,j-1)/msft(i,j-1)  + (ztop-zs(i,j))*var(i,k,j)/msft(i,j) ) *v(i,k,j) )*rdy					 

   

       ENDDO
      ENDDO
      ENDDO
	  
	  
 
	
		


  
        i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
         	divg(i,k,j) = divg(i,k,j)  * msft(i,j)*msft(i,j)/(ztop-zs(i,j)) 
      ENDDO
      ENDDO
      ENDDO

      END SUBROUTINE calcu_divg2


SUBROUTINE calcu_divg_test( var, divg,                 &
                        u, v, sgmdot,                       &
                     config_flags,                  &
                        msft,zh,zs,fzm,fzp,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   ! Input data
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var,u,     &
                                                                      v, &
                                                                      sgmdot
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      !   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: mut
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: divg

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh,fzm,fzp


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   ! Local data
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx
         
        
       !  REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: zzh,dzzh
	   
	   
	   
	   
	   
		 
		    REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) ::  varw
		 
		 
		 
		 


      REAL , DIMENSION( kms:kme ) :: sgmh,dsgm
       REAL , DIMENSION( ims:ime , jms:jme )  ::  rmsft,rmm
	   
	   REAL :: ub,vb 
				 
				 

    LOGICAL :: specified

   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

   
     ktf=MIN(kte,kde-1)
 

	  DO k=kts,kde-1
          sgmh(k) = 0.5*(zh(k)+zh(k+1))
          dsgm(k) = zh(k+1) - zh(k) 
      ENDDO


 



	  

   ! 20230228  for rho couple scalar 
      CALL set_physical_bc3d( var, 'p', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              its, ite, jts, jte, kts, kte,       &
                              its, ite, jts, jte, kts, kte )
	  
	  

     divg = 0.0 
	 
	 
	 
	 
	  	call average_vertical_w( varw, var,                 &
                     config_flags,                  &
                         fzm,fzp,                    &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )	


!  ---- z  direction -----------

  !  计划拆分为简单公式，处理变量从半层插值到整层，及边界处理问题
  
  
        i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
	  
      DO  j = j_start, j_end
      DO  i = i_start, i_end
        rmsft(i,j)  = 1./msft(i,j) 
        rmm(i,j) = 1./(msft(i,j)*msft(i,j) )
      ENDDO
      ENDDO             

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end



            divg(i,k,j) = divg(i,k,j) -    ( varw(i,k+1,j)* rmm(i,j) * sgmdot(i,k+1,j) - varw(i,k,j)* rmm(i,j) * sgmdot(i,k,j)   ) /dsgm(k) * (ztop-zs(i,j))


      ENDDO
      ENDDO
      ENDDO	 
	 
	 
	 
	 
	 
!  ---- x direction -----------
 

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)


      IF ( .NOT. config_flags%periodic_x ) THEN
        IF ( config_flags%open_xs .or. specified ) i_start = MAX(ids+1,its)
        IF ( config_flags%open_xe .or. specified ) i_end   = MIN(ide-2,ite)
      ENDIF

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
  
    !  	divg(i,k,j) =  -  (  0.5*((ztop-zs(i,j))*var(i,k,j)*rmsft(i,j)  + (ztop-zs(i+1,j))*var(i+1,k,j)*rmsft(i+1,j) ) *u(i+1,k,j) &
    !                -   0.5*((ztop-zs(i-1,j))*var(i-1,k,j)*rmsft(i-1,j)  + (ztop-zs(i,j))*var(i,k,j)*rmsft(i,j) ) *u(i,k,j) )*rdx

      	divg(i,k,j) =  divg(i,k,j) -  (  0.5*((ztop-zs(i,j))*var(i,k,j)/msft(i,j)  + (ztop-zs(i+1,j))*var(i+1,k,j)/msft(i+1,j) ) *u(i+1,k,j) &
                     -   0.5*((ztop-zs(i-1,j))*var(i-1,k,j)/msft(i-1,j)  + (ztop-zs(i,j))*var(i,k,j)/msft(i,j) ) *u(i,k,j) )*rdx
      ENDDO
      ENDDO
      ENDDO
	  
	  
	IF ( config_flags%open_xs  .and.  its == ids ) THEN 
		
	!	i = ids  !  用于方法1 
		
      DO  j = j_start, j_end
      DO k=kts,ktf
  
  
 !  方法1 ，直接改造原来的格式。不用i-1点的var
     ! 	divg(i,k,j) =  -  (  0.5*((ztop-zs(i,j))*var(i,k,j)/msft(i,j)  + (ztop-zs(i+1,j))*var(i+1,k,j)/msft(i+1,j) ) *u(i+1,k,j) &
     !                -    (  (ztop-zs(i,j))*var(i,k,j)/msft(i,j) ) *u(i,k,j) )*rdx


					 
	! 参考wrf	,这一方法的使用，改进非常大。避免了hill2km_x试验的西边界问题。			 
	         ub = MIN( 0.5*(u(its,k,j)+u(its+1,k,j)), 0. )
		 
         divg(its,k,j) =   divg(its,k,j)                   &
               - (ztop-zs(its,j))/msft(its,j)*rdx*(                                             &
                       ub*(   var(its+1,k,j)                 &
                            - var(its  ,k,j)   ) +           &
                       var(its,k,j)*(u(its+1,k,j)-u(its,k,j))  &       ! modified 
                                                                )				 
					 
      ENDDO
      ENDDO
 
	  
	ENDIF	  
	

	  
	    IF ( config_flags%open_xe  .and.  ite == ide ) THEN 
		
	 !	i = ite-1
		
		 i_end   = MIN(ite,ide-1)
		 
		
      DO  j = j_start, j_end
      DO k=kts,ktf
 
  
      ! 	divg(i,k,j) =  -  (  ((ztop-zs(i,j))*var(i,k,j)/msft(i,j)    ) *u(i+1,k,j) &
      !                -   0.5*((ztop-zs(i-1,j))*var(i-1,k,j)/msft(i-1,j)  + (ztop-zs(i,j))*var(i,k,j)/msft(i,j) ) *u(i,k,j) )*rdx
	
	     ub = MAX( 0.5*(u(ite-1,k,j)+u(ite,k,j)), 0. )	 
		 
         divg(i_end,k,j) =   divg(i_end,k,j)                &
               - (ztop-zs(i_end,j))/msft(i_end,j)*rdx*(                                               &
                       ub*(  var(i_end  ,k,j)                  &
                           - var(i_end-1,k,j) ) +              &
                       var(i_end,k,j)*(u(ite,k,j)-u(ite-1,k,j))  &   ! modified  
                                                                    )				 
					 
					 
      ENDDO
      ENDDO
  

	ENDIF	
	  
	  
	  
	  

!  ---- y direction -----------

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
      ! Polar boundary conditions are like open or specified
      IF ( config_flags%open_ys .or. specified .or. config_flags%polar ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified .or. config_flags%polar ) j_end   = MIN(jde-2,jte) 	  

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end

   !   	divg(i,k,j) = divg(i,k,j) -   ( 0.5*((ztop-zs(i,j))*var(i,k,j)*rmsft(i,j)  + (ztop-zs(i,j+1))*var(i,k,j+1)*rmsft(i,j+1) ) *v(i,k,j+1)    &
   !                  -   0.5*((ztop-zs(i,j-1))*var(i,k,j-1)*rmsft(i,j-1)  + (ztop-zs(i,j))*var(i,k,j)*rmsft(i,j) ) *v(i,k,j) )*rdy
					 
 		 
    	divg(i,k,j) = divg(i,k,j) -   ( 0.5*((ztop-zs(i,j))*var(i,k,j)/msft(i,j)  + (ztop-zs(i,j+1))*var(i,k,j+1)/msft(i,j+1) ) *v(i,k,j+1)    &
                 -   0.5*((ztop-zs(i,j-1))*var(i,k,j-1)/msft(i,j-1)  + (ztop-zs(i,j))*var(i,k,j)/msft(i,j) ) *v(i,k,j) )*rdy					 

   

       ENDDO
      ENDDO
      ENDDO
	  
	  

	  
	    IF ( config_flags%open_ys  .and.  jts == jds ) THEN 
		
	!	j = jds 
		
      DO k=kts,ktf
      DO  i = i_start, i_end		
		
    !  	divg(i,k,j) = divg(i,k,j) -   ( 0.5*((ztop-zs(i,j))*var(i,k,j)/msft(i,j)  + (ztop-zs(i,j+1))*var(i,k,j+1)/msft(i,j+1) ) *v(i,k,j+1)    &
    !                 -    ( (ztop-zs(i,j))*var(i,k,j)/msft(i,j) ) *v(i,k,j) )*rdy
	
	
	
	      vb = MIN( 0.5*(v(i,k,jts)+v(i,k,jts+1)), 0. )
		  
        divg(i,k,jts) =    divg(i,k,jts)                &
               - (ztop-zs(i,jts))/msft(i,jts)*rdy*(                                             &
                       vb*(  var(i,k,jts+1)                  &
                           - var(i,k,jts  ) ) +              &
                       var(i,k,jts)*(v(i,k,jts+1)-v(i,k,jts))  &     ! modified  
                                                               )
					 
      ENDDO
      ENDDO
	    ENDIF	
		
 
	  
	    IF ( config_flags%open_ye  .and.  jte == jde ) THEN 
		
	 ! 	j = jte-1
		
		j_end   = MIN(jte,jde-1)
		
      DO k=kts,ktf
      DO  i = i_start, i_end		
		
	  !      divg(i,k,j) = divg(i,k,j) -   (  ((ztop-zs(i,j))*var(i,k,j)/msft(i,j) )   *v(i,k,j+1)    &
       !              -   0.5*((ztop-zs(i,j-1))*var(i,k,j-1)/msft(i,j-1)  + (ztop-zs(i,j))*var(i,k,j)/msft(i,j) ) *v(i,k,j) )*rdy					 
	
		  vb = MAX( 0.5*(v(i,k,jte-1)+v(i,k,jte)), 0. ) 
		 
         divg(i,k,j_end) =    divg(i,k,j_end)               &
               -  (ztop-zs(i,j_end))/msft(i,j_end)*rdy*(                                               &
                       vb*(   var(i,k,j_end  )                 &
                            - var(i,k,j_end-1) ) +             &
                       var(i,k,j_end)*(v(i,k,jte)-v(i,k,jte-1))  &    ! modified   
                                                                    )
					 
      ENDDO
      ENDDO
	    ENDIF

		
	  
	  
        i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1) 



      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
         	divg(i,k,j) = divg(i,k,j)  * msft(i,j)*msft(i,j)/(ztop-zs(i,j)) 
      ENDDO
      ENDDO
      ENDDO

      END SUBROUTINE calcu_divg_test
	  
	  
	  
! 改变流程，先边界，后。

SUBROUTINE calcu_divg_bak2024( var, divg,                 &
                        u, v, sgmdot,                       &
                     config_flags,                  &
                        msft,zh,zs,fzm,fzp,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   ! Input data
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var,u,     &
                                                                      v, &
                                                                      sgmdot
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      !   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: mut
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: divg

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh,fzm,fzp


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   ! Local data
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx
         
        
       !  REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: zzh,dzzh
	   
	   
	   
	   
	   
		 
		    REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) ::  varw
		 
		 
		 
		 


      REAL , DIMENSION( kms:kme ) :: sgmh,dsgm
       REAL , DIMENSION( ims:ime , jms:jme )  ::  rmm    ! rmsft 
	   
	   REAL :: ub,vb 
				 
				 
		!20230816		 
				 
	   REAL,  DIMENSION( its-1:ite+2, kts:kte, jts-1:jte+2  ) :: fqx, fqy, fqz
	   
	       REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) ::  field

   LOGICAL :: degrade_xs, degrade_ys
   LOGICAL :: degrade_xe, degrade_ye
   
   
		   

    LOGICAL :: specified




   degrade_xs = .true.
   degrade_xe = .true.
   degrade_ys = .true.
   degrade_ye = .true.


   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

   
     ktf=MIN(kte,kde-1)
 

	  DO k=kts,kde-1
          sgmh(k) = 0.5*(zh(k)+zh(k+1))
          dsgm(k) = zh(k+1) - zh(k) 
      ENDDO


 



	  

   ! 20230228  for rho couple scalar 
      CALL set_physical_bc3d( var, 'p', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              its, ite, jts, jte, kts, kte,       &
                              its, ite, jts, jte, kts, kte )
	  
	  

     divg = 0.0 
	 
	 
! boundary	 
	 
!  ---- x direction -----------
 

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)


      IF ( .NOT. config_flags%periodic_x ) THEN
        IF ( config_flags%open_xs .or. specified ) i_start = MAX(ids+1,its)
        IF ( config_flags%open_xe .or. specified ) i_end   = MIN(ide-2,ite)
      ENDIF

      
	  
	  
	IF ( config_flags%open_xs  .and.  its == ids ) THEN 
		
	!	i = ids  !  用于方法1 
		
      DO  j = j_start, j_end
      DO k=kts,ktf
   
 					 
	! 参考wrf	,这一方法的使用，改进非常大。避免了hill2km_x试验的西边界问题。			 
	         ub = MIN( 0.5*(u(its,k,j)+u(its+1,k,j)), 0. )
		 
         divg(its,k,j) =                      &
               - (ztop-zs(its,j))/msft(its,j)*rdx*(                                             &
                       ub*(   var(its+1,k,j)                 &
                            - var(its  ,k,j)   ) +           &
                       var(its,k,j)*(u(its+1,k,j)-u(its,k,j))  &       ! modified 
                                                                )				 
					 
      ENDDO
      ENDDO
 
	  
	ENDIF	  
!	write(6,*) '1 divg(1,1,1) = ', divg(1,1,1)  

	  
	    IF ( config_flags%open_xe  .and.  ite == ide ) THEN 
		
	 !	i = ite-1
		
		 i_end   = MIN(ite,ide-1)
		 
		
      DO  j = j_start, j_end
      DO k=kts,ktf
 
  
     
	     ub = MAX( 0.5*(u(ite-1,k,j)+u(ite,k,j)), 0. )	 
		 
         divg(i_end,k,j) =                  &
               - (ztop-zs(i_end,j))/msft(i_end,j)*rdx*(                                               &
                       ub*(  var(i_end  ,k,j)                  &
                           - var(i_end-1,k,j) ) +              &
                       var(i_end,k,j)*(u(ite,k,j)-u(ite-1,k,j))  &   ! modified  
                                                                    )				 
					 
					 
      ENDDO
      ENDDO
  

	ENDIF	
	  
	  
	!  write(6,*) '2 divg(1,1,1) = ', divg(1,1,1)  
	  

!  ---- y direction -----------

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
      ! Polar boundary conditions are like open or specified
      IF ( config_flags%open_ys .or. specified .or. config_flags%polar ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified .or. config_flags%polar ) j_end   = MIN(jde-2,jte) 	  

       
	  
	    IF ( config_flags%open_ys  .and.  jts == jds ) THEN 
		
	!	j = jds 
		
      DO k=kts,ktf
      DO  i = i_start, i_end		

	
	      vb = MIN( 0.5*(v(i,k,jts)+v(i,k,jts+1)), 0. )
		  
        divg(i,k,jts) =    divg(i,k,jts)                &
               - (ztop-zs(i,jts))/msft(i,jts)*rdy*(                                             &
                       vb*(  var(i,k,jts+1)                  &
                           - var(i,k,jts  ) ) +              &
                       var(i,k,jts)*(v(i,k,jts+1)-v(i,k,jts))  &     ! modified  
                                                               )
					 
      ENDDO
      ENDDO
	    ENDIF	
		
  !      write(6,*) '3 divg(1,1,1) = ', divg(1,1,1) 
	  
	    IF ( config_flags%open_ye  .and.  jte == jde ) THEN 
		
	 ! 	j = jte-1
		
		j_end   = MIN(jte,jde-1)
		
      DO k=kts,ktf
      DO  i = i_start, i_end		
		
	 
	
		  vb = MAX( 0.5*(v(i,k,jte-1)+v(i,k,jte)), 0. ) 
		 
         divg(i,k,j_end) =    divg(i,k,j_end)               &
               -  (ztop-zs(i,j_end))/msft(i,j_end)*rdy*(                                               &
                       vb*(   var(i,k,j_end  )                 &
                            - var(i,k,j_end-1) ) +             &
                       var(i,k,j_end)*(v(i,k,jte)-v(i,k,jte-1))  &    ! modified   
                                                                    )
					 
      ENDDO
      ENDDO
	    ENDIF
		
		
     !     write(6,*) '4 divg(1,1,1) = ', divg(1,1,1) 







!!!!!!!!!!!!!!!!!!!!!!!!!!!	 
	 
      field  = var  	 
	 
	 
!  ---- z  direction -----------	




! fqz

!-------------------- vertical advection




	  
!-- loop bounds for periodic or sym conditions

      i_start = its-1
      i_end   = MIN(ite,ide-1)+1
      j_start = jts-1
      j_end   = MIN(jte,jde-1)+1

!-- loop bounds for open or specified conditions

    IF(degrade_xs) i_start = MAX(its-1,ids)
    IF(degrade_xe) i_end   = MIN(ite+1,ide-1)
    IF(degrade_ys) j_start = MAX(jts-1,jds)
    IF(degrade_ye) j_end   = MIN(jte+1,jde-1)

    

      DO j = j_start, j_end

         DO i = i_start, i_end
           fqz(i,1,j)  = 0.
        !   fqzl(i,1,j) = 0.
           fqz(i,kde,j)  = 0.
          ! fqzl(i,kde,j) = 0.
         ENDDO

         DO k=kts+1,ktf
         DO i = i_start, i_end


		   
	   fqz(i,k,j) = sgmdot(i,k,j) *(fzm(k)*field(i,k,j)+fzp(k)*field(i,k-1,j))*(ztop-zs(i,j))/msft(i,j)    !这里不是m**2
		   

        ENDDO
        ENDDO

      ENDDO

    
   



! fqz end 
  
	  
	  	call average_vertical_w( varw, var,                 &
                     config_flags,                  &
                         fzm,fzp,                    &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )	



  !  计划拆分为简单公式，处理变量从半层插值到整层，及边界处理问题
  
  
        i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
	  
      DO  j = j_start, j_end
      DO  i = i_start, i_end
      !  rmsft(i,j)  = 1./msft(i,j) 
        rmm(i,j) = 1./(msft(i,j)*msft(i,j) )
      ENDDO
      ENDDO             

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end


           ! 20230817测试不用fqz，用原来方案。
             divg(i,k,j) = divg(i,k,j) -    ( varw(i,k+1,j)* rmm(i,j) * sgmdot(i,k+1,j) - varw(i,k,j)* rmm(i,j) * sgmdot(i,k,j)   ) /dsgm(k) * (ztop-zs(i,j))

            !      divg(i,k,j) = divg(i,k,j) -rdzw(k)*( fqz (i,k+1,j)-fqz (i,k,j) )
									  
			 !       divg(i,k,j) = divg(i,k,j) - ( fqz (i,k+1,j)-fqz (i,k,j) ) /dsgm(k)		

 

			 
									  
      ENDDO
      ENDDO
      ENDDO
	 
  !write(6,*) 'max divgw =  ', maxval(divg) , maxloc(divg) 
  !write(6,*) 'min divgw =  ', minval(divg) , minloc(divg) 
!!!!!!!!!!!!!!!!!!!!!
! fqx,fqy

  
		 
        
   IF( config_flags%periodic_x   .or. &
       config_flags%symmetric_xs .or. &
       (its > ids+1)                ) degrade_xs = .false.
   IF( config_flags%periodic_x   .or. &
       config_flags%symmetric_xe .or. &
       (ite < ide-2)                ) degrade_xe = .false.
   IF( config_flags%periodic_y   .or. &
       config_flags%symmetric_ys .or. &
       (jts > jds+1)                ) degrade_ys = .false.
   IF( config_flags%periodic_y   .or. &
       config_flags%symmetric_ye .or. &
       (jte < jde-2)                ) degrade_ye = .false.

!--  y flux compute; these bounds are for periodic and sym b.c.

      ktf=MIN(kte,kde-1)
      i_start = its-1
      i_end   = MIN(ite,ide-1)+1
      j_start = jts-1
      j_end   = MIN(jte,jde-1)+1

!--  modify loop bounds if open or specified

      IF(degrade_xs) i_start = its
      IF(degrade_xe) i_end   = MIN(ite,ide-1)
      IF(degrade_ys) j_start = MAX(jts,jds+1)
      IF(degrade_ye) j_end = MIN(jte,jde-2)


      DO j = j_start, j_end+1
        DO k=kts,ktf
        DO i = i_start, i_end
 
		   
	!	   fqy(i,k,j) =  0.5*v(i,k,j)*          &
    !              ((ztop-zs(i,j))*field(i,k,j)/msft(i,j) + (ztop-zs(i,j-1))*field(i,k,j-1)/msft(i,j-1))	

				  
	 ! 20240228	    fqy(i,k,j) =    0.5*((ztop-zs(i,j-1))*var(i,k,j-1)*rmsft(i,j-1)  + (ztop-zs(i,j))*var(i,k,j)*rmsft(i,j) ) *v(i,k,j) 	  
		   
		    fqy(i,k,j) =    0.5*((ztop-zs(i,j-1))*var(i,k,j-1)/msft(i,j-1)  + (ztop-zs(i,j))*var(i,k,j)/msft(i,j) ) *v(i,k,j) 	
        ENDDO
        ENDDO
      ENDDO
	  
	  
	  ! 20230818 14:12 
	  j_start = jts 
      j_end   = MIN(jte,jde-1) 
	  
	  
      DO j = j_start, j_end
        DO k=kts,ktf
        DO i = i_start, i_end+1
 
				  
      !      fqx( i,k,j ) = 0.5*u(i,k,j)*          &
      !            ((ztop-zs(i,j))*field(i,k,j)/msft(i,j) + (ztop-zs(i-1,j))*field(i-1,k,j)/msft(i-1,j))			

		 	 fqx( i,k,j ) =  	  0.5*((ztop-zs(i-1,j))*var(i-1,k,j)/msft(i-1,j)  + (ztop-zs(i,j))*var(i,k,j)/msft(i,j) ) *u(i,k,j)   
 
 
 
 
        ENDDO
        ENDDO
      ENDDO
	  
	  


!!!!!!!!!!!!!!!!!!!!!!!	 
	 
	 
	 
!  ---- x direction -----------
 

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)


      IF ( .NOT. config_flags%periodic_x ) THEN
        IF ( config_flags%open_xs .or. specified ) i_start = MAX(ids+1,its)
        IF ( config_flags%open_xe .or. specified ) i_end   = MIN(ide-2,ite)
      ENDIF

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
  
    

       !  	divg(i,k,j) =  divg(i,k,j) -  (  0.5*((ztop-zs(i,j))*var(i,k,j)/msft(i,j)  + (ztop-zs(i+1,j))*var(i+1,k,j)/msft(i+1,j) ) *u(i+1,k,j) &
         !             -   0.5*((ztop-zs(i-1,j))*var(i-1,k,j)/msft(i-1,j)  + (ztop-zs(i,j))*var(i,k,j)/msft(i,j) ) *u(i,k,j) )*rdx
					 
     	divg(i,k,j) =  divg(i,k,j) -   rdx*( fqx (i+1,k,j)-fqx (i,k,j)   )		
 					 
					 
      ENDDO
      ENDDO
      ENDDO
	  
 !write(6,*) 'max divgw x =  ', maxval(divg) , maxloc(divg) 
  !write(6,*) 'min divgw x =  ', minval(divg) , minloc(divg) 

!  ---- y direction -----------

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
      ! Polar boundary conditions are like open or specified
      IF ( config_flags%open_ys .or. specified .or. config_flags%polar ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified .or. config_flags%polar ) j_end   = MIN(jde-2,jte) 	  

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end

    !   	divg(i,k,j) = divg(i,k,j) -   ( 0.5*((ztop-zs(i,j))*var(i,k,j)*rmsft(i,j)  + (ztop-zs(i,j+1))*var(i,k,j+1)*rmsft(i,j+1) ) *v(i,k,j+1)    &
    !                                -   0.5*((ztop-zs(i,j-1))*var(i,k,j-1)*rmsft(i,j-1)  + (ztop-zs(i,j))*var(i,k,j)*rmsft(i,j) ) *v(i,k,j) )*rdy
					 
 		 
     	divg(i,k,j) = divg(i,k,j) -    rdy*( fqy (i,k,j+1)-fqy (i,k,j) )  					 

    

       ENDDO
      ENDDO
      ENDDO
	  
	  
  !write(6,*) 'max divgw x y=  ', maxval(divg) , maxloc(divg) 
  !write(6,*) 'min divgw x y=  ', minval(divg) , minloc(divg)     
	
		


  
        i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
         	divg(i,k,j) = divg(i,k,j)  * msft(i,j)*msft(i,j)/(ztop-zs(i,j)) 
      ENDDO
      ENDDO
      ENDDO
	  
  !write(6,*) 'max divgw x y/z=  ', maxval(divg) , maxloc(divg) 
  !write(6,*) 'min divgw x y/z=  ', minval(divg) , minloc(divg)  	 

      END SUBROUTINE calcu_divg_bak2024

	  
SUBROUTINE advect_scalar_pd_fm   ( field, field_old, tendency,    &
                                h_tendency, z_tendency,        &
                                u, v, sigmdot,                   &
                                c1, c2,                        &
                                mut, mub, mu_old,              &
                                time_step, config_flags,       &
                                tenddec,                       &
                                msfux, msfuy, msfvx, msfvy,    &
                                msftx, msfty,                  &
								msft,   &
                                fzm, fzp, ztop,zs,                       &
                                rdx, rdy, rdzw, dt,            &
                                ids, ide, jds, jde, kds, kde,  &
                                ims, ime, jms, jme, kms, kme,  &
                                its, ite, jts, jte, kts, kte  )

!  this is a first cut at a positive definite advection option
!  for scalars in WRF.  This version is memory intensive ->
!  we save 3d arrays of x, y and z both high and low order fluxes
!  (six in all).  Alternatively, we could sweep in a direction
!  and lower the cost considerably.

!  uses the Smolarkiewicz MWR 1989 approach, with addition of first-order
!  fluxes initially

!  WCS, 3 December 2002, 24 February 2003

   IMPLICIT NONE
   
   ! Input data
   
   TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

   LOGICAL ,                 INTENT(IN   ) :: tenddec  ! tendency flag

   INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: field,     &
                                                                      field_old, &
                                                                      u,    &
                                                                      v,    &
                                                                      sigmdot

   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: mut, mub, mu_old
   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) :: tendency
   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(  OUT) :: h_tendency, z_tendency

   REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msfux,  &
                                                                    msfuy,  &
                                                                    msfvx,  &
                                                                    msfvy,  &
                                                                    msftx,  &
                                                                    msfty,msft,zs  

   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: fzm,  &
                                                                  fzp,  &
                                                                  rdzw, &
                                                                  c1,   &
                                                                  c2

   REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy,  &
                                                                  dt,ztop 
   INTEGER ,                                     INTENT(IN   ) :: time_step

   ! Local data
   
   INTEGER :: i, j, k, itf, jtf, ktf
   INTEGER :: i_start, i_end, j_start, j_end
   INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
   INTEGER :: jmin, jmax, jp, jm, imin, imax

   REAL    :: mrdx, mrdy, ub, vb, uw, vw, mu

!  storage for high and low order fluxes

   REAL,  DIMENSION( its-1:ite+2, kts:kte, jts-1:jte+2  ) :: fqx, fqy, fqz

   REAL,  DIMENSION( its-1:ite+2, kts:kte, jts-1:jte+2  ) :: fqxl, fqyl, fqzl

   INTEGER :: horz_order, vert_order
   
   LOGICAL :: degrade_xs, degrade_ys
   LOGICAL :: degrade_xe, degrade_ye

   INTEGER :: jp1, jp0, jtmp

   REAL,DIMENSION( its-1:ite+2, kts:kte, jts-1:jte+2  ) :: flux_out, ph_low
   REAL :: scale
   !REAL :: flux_out, ph_low, scale
   REAL, PARAMETER :: eps=1.e-20


! definition of flux operators, 3rd, 4th, 5th or 6th order

   REAL    :: flux3, flux4, flux5, flux6, flux_upwind
   REAL    :: q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua, vel, cr

      flux4(q_im2, q_im1, q_i, q_ip1, ua) =                     &
            (7./12.)*(q_i + q_im1) - (1./12.)*(q_ip1 + q_im2)

      flux3(q_im2, q_im1, q_i, q_ip1, ua) =                     &
           flux4(q_im2, q_im1, q_i, q_ip1, ua) +                &
           sign(1,time_step)*sign(1.,ua)*(1./12.)*((q_ip1 - q_im2)-3.*(q_i-q_im1))

      flux6(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua) =       &
            (37./60.)*(q_i+q_im1) - (2./15.)*(q_ip1+q_im2)      &
            +(1./60.)*(q_ip2+q_im3)

      flux5(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua) =       &
           flux6(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua)    &
            -sign(1,time_step)*sign(1.,ua)*(1./60.)*(           &
              (q_ip2-q_im3)-5.*(q_ip1-q_im2)+10.*(q_i-q_im1) )

      flux_upwind(q_im1, q_i, cr ) = 0.5*min( 1.0,(cr+abs(cr)))*q_im1 &
                                    +0.5*max(-1.0,(cr-abs(cr)))*q_i

!      flux_upwind(q_im1, q_i, cr ) = 0.5*(1.+sign(1.,cr))*q_im1 &
!                                    +0.5*(1.-sign(1.,cr))*q_i
!      flux_upwind(q_im1, q_i, cr ) = 0.

    REAL     :: dx,dy,dz

    LOGICAL, PARAMETER :: pd_limit = .true.

! set order for the advection schemes

!  write(6,*) ' in pd advection routine '

    ! Empty arrays just in case:
    IF (config_flags%polar) THEN
       fqx(:,:,:)  = 0.
       fqy(:,:,:)  = 0.
       fqz(:,:,:)  = 0.
       fqxl(:,:,:) = 0.
       fqyl(:,:,:) = 0.
       fqzl(:,:,:) = 0.
    END IF

  ktf=MIN(kte,kde-1)
  horz_order = config_flags%h_sca_adv_order
  vert_order = config_flags%v_sca_adv_order

!  determine boundary mods for flux operators
!  We degrade the flux operators from 3rd/4th order
!   to second order one gridpoint in from the boundaries for
!   all boundary conditions except periodic and symmetry - these
!   conditions have boundary zone data fill for correct application
!   of the higher order flux stencils

   degrade_xs = .true.
   degrade_xe = .true.
   degrade_ys = .true.
   degrade_ye = .true.

!  begin with horizontal flux divergence
!  here is the choice of flux operators


  horizontal_order_test : IF( horz_order == 6 ) THEN

 

   ELSE IF( horz_order == 2 ) THEN

   IF( config_flags%periodic_x   .or. &
       config_flags%symmetric_xs .or. &
       (its > ids+1)                ) degrade_xs = .false.
   IF( config_flags%periodic_x   .or. &
       config_flags%symmetric_xe .or. &
       (ite < ide-2)                ) degrade_xe = .false.
   IF( config_flags%periodic_y   .or. &
       config_flags%symmetric_ys .or. &
       (jts > jds+1)                ) degrade_ys = .false.
   IF( config_flags%periodic_y   .or. &
       config_flags%symmetric_ye .or. &
       (jte < jde-2)                ) degrade_ye = .false.

!--  y flux compute; these bounds are for periodic and sym b.c.

      ktf=MIN(kte,kde-1)
      i_start = its-1
      i_end   = MIN(ite,ide-1)+1
      j_start = jts-1
      j_end   = MIN(jte,jde-1)+1

!--  modify loop bounds if open or specified

      IF(degrade_xs) i_start = its
      IF(degrade_xe) i_end   = MIN(ite,ide-1)
      IF(degrade_ys) j_start = MAX(jts,jds+1)
      IF(degrade_ye) j_end = MIN(jte,jde-2)

!  compute fluxes, 2nd order, y flux

      DO j = j_start, j_end+1
        DO k=kts,ktf
        DO i = i_start, i_end
 
 
           dy = 2./(msftx(i,j)+msftx(i,j-1))/rdy  ! ADT eqn 48 d/dy  ??? 
        !  mu = 0.5*((c1(k)*mut(i,j)+c2(k))+(c1(k)*mut(i,j-1)+c2(k)))
          vel =  v(i,k,j)
          cr = vel*dt/dy       !/mu
        !  fqyl(i,k,j) = mu*(dy/dt)*flux_upwind(field_old(i,k,j-1), field_old(i,k,j  ), cr)

 
		            fqyl(i,k,j) =  (dy/dt)*flux_upwind((ztop-zs(i,j-1))*field(i,k,j-1)/msft(i,j-1), (ztop-zs(i,j))*field(i,k,j)/msft(i,j), cr)
		   
		   fqy(i,k,j) =  0.5*v(i,k,j)*          &
                  ((ztop-zs(i,j))*field(i,k,j)/msft(i,j) + (ztop-zs(i,j-1))*field(i,k,j-1)/msft(i,j-1))				  
				  
		       fqy(i,k,j) = fqy(i,k,j) - fqyl(i,k,j)
			   
			   
			   
        ENDDO
        ENDDO
      ENDDO

!  next, x flux



!   20230818 
	!  j_start = jts 
    ! j_end   = MIN(jte,jde-1) 
      IF(degrade_xs) i_start = its
      IF(degrade_xe) i_end   = MIN(ite,ide-1)
      IF(degrade_ys) j_start = jts  
      IF(degrade_ye) j_end = MIN(jte,jde-1)  
	  


      DO j = j_start, j_end
        DO k=kts,ktf
        DO i = i_start, i_end+1
            dx = 2./(msfty(i,j)+msfty(i-1,j))/rdx ! ADT eqn 48 d/dx   ???
          !  mu = 0.5*((c1(k)*mut(i,j)+c2(k))+(c1(k)*mut(i-1,j)+c2(k)))
           vel =  u(i,k,j)
         cr = vel*dt/dx !/mu
          fqxl(i,k,j) =  (dx/dt)*flux_upwind(field_old(i-1,k,j), field_old(i,k,j  ), cr)
          !  fqx( i,k,j ) = 0.5*ru(i,k,j)*          &
           !         (field(i,k,j)+field(i-1,k,j))
				  
            fqx( i,k,j ) = 0.5*u(i,k,j)*          &
                  ((ztop-zs(i,j))*field(i,k,j)/msft(i,j) + (ztop-zs(i-1,j))*field(i-1,k,j)/msft(i-1,j))				  
				  
				  

          fqx(i,k,j) = fqx(i,k,j) - fqxl(i,k,j)
        ENDDO
        ENDDO
      ENDDO

!--- end of 2nd order horizontal flux calculation

   ELSE

      WRITE ( wrf_err_message , * ) 'module_advect: advect_scalar_pd, h_order not known ',horz_order
      CALL wrf_error_fatal ( TRIM( wrf_err_message ) )

   ENDIF horizontal_order_test






















!  pick up the rest of the horizontal radiation boundary conditions.
!  (these are the computations that don't require 'cb'.
!  first, set to index ranges

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)

!  compute x (u) conditions for v, w, or scalar

   IF( (config_flags%open_xs) .and. (its == ids) ) THEN

       DO j = j_start, j_end
       DO k = kts, ktf
         ub = MIN( 0.5*(u(its,k,j)+u(its+1,k,j)), 0. )
         tendency(its,k,j) = tendency(its,k,j)                     &
             !  - rdx*(                                             &
             !          ub*(   field_old(its+1,k,j)                 &
             !               - field_old(its  ,k,j)   ) +           &
             !          field(its,k,j)*(ru(its+1,k,j)-ru(its,k,j))  &
             !                                                   )
                 - (ztop-zs(its,j))/msft(its,j)*rdx*(                                             &
                       ub*(   field_old(its+1,k,j)                 &
                            - field_old(its  ,k,j)   ) +           &
                       field_old(its,k,j)*(u(its+1,k,j)-u(its,k,j))  &       ! modified 
                                                                ) * msft(its,j)*msft(its,j)/(ztop-zs(its,j))  																	
																
																
       ENDDO
       ENDDO

   ENDIF

   IF( (config_flags%open_xe) .and. (ite == ide) ) THEN

       DO j = j_start, j_end
       DO k = kts, ktf
         ub = MAX( 0.5*(u(ite-1,k,j)+u(ite,k,j)), 0. )
         tendency(i_end,k,j) = tendency(i_end,k,j)                   &
        !       - rdx*(                                               &
         !              ub*(  field_old(i_end  ,k,j)                  &
         !                  - field_old(i_end-1,k,j) ) +              &
         !              field(i_end,k,j)*(ru(ite,k,j)-ru(ite-1,k,j))  &
         !                                                           )
																	
                - (ztop-zs(i_end,j))/msft(i_end,j)*rdx*(                                               &
                       ub*(  field_old(i_end  ,k,j)                  &
                           - field_old(i_end-1,k,j) ) +              &
                       field_old(i_end,k,j)*(u(ite,k,j)-u(ite-1,k,j))  &   ! modified  
                                                                    )* msft(i_end,j)*msft(i_end,j)/(ztop-zs(i_end,j)) 																		
																	
       ENDDO
       ENDDO

   ENDIF

   IF( (config_flags%open_ys) .and. (jts == jds) ) THEN

       DO i = i_start, i_end
       DO k = kts, ktf
         vb = MIN( 0.5*(v(i,k,jts)+v(i,k,jts+1)), 0. )
         tendency(i,k,jts) = tendency(i,k,jts)                     &
               - (ztop-zs(i,jts))/msft(i,jts)*rdy*(                                             &
                       vb*(  field_old(i,k,jts+1)                  &
                           - field_old(i,k,jts  ) ) +              &
                       field_old(i,k,jts)*(v(i,k,jts+1)-v(i,k,jts))  &     ! modified  
                                                               )* msft(i,jts)*msft(i,jts)/(ztop-zs(i,jts)) 


       ENDDO
       ENDDO

   ENDIF

   IF( (config_flags%open_ye) .and. (jte == jde)) THEN

       DO i = i_start, i_end
       DO k = kts, ktf
         vb = MAX( 0.5*(v(i,k,jte-1)+v(i,k,jte)), 0. )
         tendency(i,k,j_end) = tendency(i,k,j_end)                   &
               -  (ztop-zs(i,j_end))/msft(i,j_end)*rdy*(                                               &
                       vb*(   field_old(i,k,j_end  )                 &
                            - field_old(i,k,j_end-1) ) +             &
                       field_old(i,k,j_end)*(v(i,k,jte)-v(i,k,jte-1))  &    ! modified   
                                                                    )* msft(i,j_end)*msft(i,j_end)/(ztop-zs(i,j_end)) 
       ENDDO
       ENDDO

   ENDIF

 !  IF( (config_flags%polar) .and. (jts == jds) ) THEN

   

  ! ENDIF

  ! IF( (config_flags%polar) .and. (jte == jde)) THEN
 

  ! ENDIF




  ! write(6,*) 'advect_scalar_pd_fm maxval(abs(tendency) ) after openbdy  :', maxval(abs(tendency) ) 














!-------------------- vertical advection

!-- loop bounds for periodic or sym conditions

      i_start = its-1
      i_end   = MIN(ite,ide-1)+1
      j_start = jts-1
      j_end   = MIN(jte,jde-1)+1

!-- loop bounds for open or specified conditions

    IF(degrade_xs) i_start = MAX(its-1,ids)
    IF(degrade_xe) i_end   = MIN(ite+1,ide-1)
    IF(degrade_ys) j_start = MAX(jts-1,jds)
    IF(degrade_ye) j_end   = MIN(jte+1,jde-1)

    vert_order_test : IF (vert_order == 6) THEN    

 
   ELSE IF (vert_order == 2) THEN    

      DO j = j_start, j_end

         DO i = i_start, i_end
           fqz(i,1,j)  = 0.
           fqzl(i,1,j) = 0.
           fqz(i,kde,j)  = 0.
           fqzl(i,kde,j) = 0.
         ENDDO

         DO k=kts+1,ktf
         DO i = i_start, i_end

            dz = 2./(rdzw(k)+rdzw(k-1))
         !  mu = 0.5*((c1(k)*mut(i,j)+c2(k))+(c1(k)*mut(i,j)+c2(k)))
            vel = sigmdot(i,k,j)  !rom(i,k,j)
           cr = vel*dt/dz  !/mu
         !  fqzl(i,k,j) = mu*(dz/dt)*flux_upwind(field_old(i,k-1,j), field_old(i,k,j  ), cr)
        !   fqz(i,k,j)=rom(i,k,j)*(fzm(k)*field(i,k,j)+fzp(k)*field(i,k-1,j))
         
		  fqzl(i,k,j) = (dz/dt)*flux_upwind(field_old(i,k-1,j), field_old(i,k,j  ), cr)  *(ztop-zs(i,j))/msft(i,j)   

		
	   fqz(i,k,j) = sigmdot(i,k,j) *(fzm(k)*field(i,k,j)+fzp(k)*field(i,k-1,j))*(ztop-zs(i,j))/msft(i,j)    !这里不是m**2
		  fqz(i,k,j) = fqz(i,k,j) - fqzl(i,k,j)  

        ENDDO
        ENDDO

      ENDDO

   ELSE

      WRITE (wrf_err_message,*) ' advect_scalar_pd, v_order not known ',vert_order
      CALL wrf_error_fatal ( wrf_err_message )

   ENDIF vert_order_test
   
   
   
   
   
   
   
   
   
   
   
   
   
   

   IF (pd_limit) THEN

! positive definite filter

   i_start = its-1
   i_end   = MIN(ite,ide-1)+1
   j_start = jts-1
   j_end   = MIN(jte,jde-1)+1

!-- loop bounds for open or specified conditions

   IF(degrade_xs) i_start = MAX(its-1,ids)
   IF(degrade_xe) i_end   = MIN(ite+1,ide-1)
   IF(degrade_ys) j_start = MAX(jts-1,jds)
   IF(degrade_ye) j_end   = MIN(jte+1,jde-1)

   IF(config_flags%specified .or. config_flags%nested) THEN
     IF (degrade_xs) i_start = MAX(its-1,ids+1)
     IF (degrade_xe) i_end   = MIN(ite+1,ide-2)
     IF (degrade_ys) j_start = MAX(jts-1,jds+1)
     IF (degrade_ye) j_end   = MIN(jte+1,jde-2)
   END IF

   IF(config_flags%open_xs) THEN
     IF (degrade_xs) i_start = MAX(its-1,ids+1)
   END IF
   IF(config_flags%open_xe) THEN
     IF (degrade_xe) i_end   = MIN(ite+1,ide-2)
   END IF
   IF(config_flags%open_ys) THEN
     IF (degrade_ys) j_start = MAX(jts-1,jds+1)
   END IF
   IF(config_flags%open_ye) THEN
     IF (degrade_ye) j_end   = MIN(jte+1,jde-2)
   END IF
   ! ADT note:
   ! We don't want to change j_start and j_end
   ! for polar BC's since we want to calculate
   ! fluxes for directions other than y at the
   ! edge

!-- here is the limiter...

   DO j=j_start, j_end
   DO k=kts, ktf
#ifdef XEON_SIMD
!DIR$ simd
#else
!DIR$ vector always
#endif
   DO i=i_start, i_end
!   old before 20230823
  !   ph_low(i,k,j) = field_old(i,k,j) &
  !              - dt*msft(i,j)*msft(i,j)/(ztop-zs(i,j)) *(               &
   !                    rdx*(fqx(i+1,k,j)-fqx(i,k,j)) +     &
   !                    rdy*(fqy(i,k,j+1)-fqy(i,k,j))     &
   !                   +rdzw(k)*(fqz(i,k+1,j)-fqz(i,k,j)) )  
                 

   !  ph_low(i,k,j) = ((c1(k)*mub(i,j)+c2(k))+(c1(k)*mu_old(i,j)))*field_old(i,k,j) &
   !             - dt*( msftx(i,j)*msfty(i,j)*(               &
   !                    rdx*(fqxl(i+1,k,j)-fqxl(i,k,j)) +     &
  !                     rdy*(fqyl(i,k,j+1)-fqyl(i,k,j))  )    &
   !                   +msfty(i,j)*rdzw(k)*(fqzl(i,k+1,j)-fqzl(i,k,j)) )

! 20230823 
     ph_low(i,k,j) = field_old(i,k,j) &
                - dt*msft(i,j)*msft(i,j)/(ztop-zs(i,j)) *(               &
                       rdx*(fqxl(i+1,k,j)-fqxl(i,k,j)) +     &
                       rdy*(fqyl(i,k,j+1)-fqyl(i,k,j))     &
                      +rdzw(k)*(fqzl(i,k+1,j)-fqzl(i,k,j)) )  

   ENDDO
   ENDDO
   ENDDO

   DO j=j_start, j_end
   DO k=kts, ktf
!DIR$ vector always
   DO i=i_start, i_end

     flux_out(i,k,j) = dt* msft(i,j)*msft(i,j)/(ztop-zs(i,j)) *( &
                                rdx*(  max(0.,fqx (i+1,k,j))      &
                                      -min(0.,fqx (i  ,k,j)) )    &
                               +rdy*(  max(0.,fqy (i,k,j+1))      &
                                      -min(0.,fqy (i,k,j  )) )   &
 !                              +rdzw(k)*(  min(0.,fqz (i,k+1,j))      &
 !                                     -max(0.,fqz (i,k  ,j)) )   )
                               +rdzw(k)*(  max(0.,fqz (i,k+1,j))      &
                                      -min(0.,fqz (i,k  ,j)) )   )									  
									  
									  
									  

   ENDDO
   ENDDO
   ENDDO 
   
 
   DO j=j_start, j_end
   DO k=kts, ktf
!DIR$ vector always
   DO i=i_start, i_end
     IF( flux_out(i,k,j) .gt. ph_low(i,k,j) ) THEN
       scale = max(0.,ph_low(i,k,j)/(flux_out(i,k,j)+eps))
       IF( fqx (i+1,k,j) .gt. 0.) fqx(i+1,k,j) = scale*fqx(i+1,k,j)
       IF( fqx (i  ,k,j) .lt. 0.) fqx(i  ,k,j) = scale*fqx(i  ,k,j)
       IF( fqy (i,k,j+1) .gt. 0.) fqy(i,k,j+1) = scale*fqy(i,k,j+1)
       IF( fqy (i,k,j  ) .lt. 0.) fqy(i,k,j  ) = scale*fqy(i,k,j  )
!  note: z flux is opposite sign in mass coordinate because
!  vertical coordinate decreases with increasing k
       IF( fqz (i,k+1,j) .lt. 0.) fqz(i,k+1,j) = scale*fqz(i,k+1,j)
       IF( fqz (i,k  ,j) .gt. 0.) fqz(i,k  ,j) = scale*fqz(i,k  ,j)

     END IF

   ENDDO
   ENDDO
   ENDDO

   END IF
   
   


! add in the pd-limited flux divergence

  i_start = its
  i_end   = MIN(ite,ide-1)
  j_start = jts
  j_end   = MIN(jte,jde-1)

  DO j = j_start, j_end
  DO k = kts, ktf
!DEC$ vector always
  DO i = i_start, i_end

     tendency (i,k,j) = tendency(i,k,j)                           &
                            -rdzw(k)* msft(i,j)*msft(i,j)/(ztop-zs(i,j))*( fqz (i,k+1,j)-fqz (i,k,j)  &
                                      +fqzl(i,k+1,j)-fqzl(i,k,j)  &
									  )

  ENDDO
  ENDDO
  ENDDO
  
  !write(6,*) 'advect_scalar_pd_fm maxval(abs(tendency) ) fqz :', maxval(abs(tendency(ids:ide-1,kts:kte-1,jds:jde-1)) ) 
  
  IF(tenddec) THEN
  DO j = j_start, j_end
  DO k = kts, ktf
  DO i = i_start, i_end

     z_tendency (i,k,j) = 0. -rdzw(k)* msft(i,j)*msft(i,j)/(ztop-zs(i,j))*( fqz (i,k+1,j)-fqz (i,k,j)  &
                                     +fqzl(i,k+1,j)-fqzl(i,k,j)    &
									  )

  ENDDO
  ENDDO
  ENDDO
  END IF

! x flux divergence
!
  IF(degrade_xs) i_start = MAX(its,ids+1)
  IF(degrade_xe) i_end   = MIN(ite,ide-2)

  DO j = j_start, j_end
  DO k = kts, ktf
!DEC$ vector always  
  DO i = i_start, i_end

     ! Un-"canceled" map scale factor, ADT Eq. 48
     tendency (i,k,j) = tendency(i,k,j)                           &
               - msft(i,j)*msft(i,j)/(ztop-zs(i,j))*( rdx*( fqx (i+1,k,j)-fqx (i,k,j)     &
                                   +fqxl(i+1,k,j)-fqxl(i,k,j)   &
								   ) )

  ENDDO
  ENDDO
  ENDDO
  
   ! write(6,*) 'advect_scalar_pd_fm maxval(abs(tendency) ) fqx :', maxval(abs(tendency(ids:ide-1,kts:kte-1,jds:jde-1)) ) 

  IF(tenddec) THEN
  DO j = j_start, j_end
  DO k = kts, ktf
  DO i = i_start, i_end

     h_tendency (i,k,j) = 0.                                      &
               - msft(i,j)*msft(i,j)/(ztop-zs(i,j))*( rdx*( fqx (i+1,k,j)-fqx (i,k,j)     &
                                   +fqxl(i+1,k,j)-fqxl(i,k,j)     &
								   ))

  ENDDO
  ENDDO
  ENDDO
  END IF

! y flux divergence
!
  i_start = its
  i_end   = MIN(ite,ide-1)
  IF(degrade_ys) j_start = MAX(jts,jds+1)
  IF(degrade_ye) j_end   = MIN(jte,jde-2)

  DO j = j_start, j_end
  DO k = kts, ktf
!DEC$ vector always
  DO i = i_start, i_end

     ! Un-"canceled" map scale factor, ADT Eq. 48
     ! It is correct to use msftx (and not msfty), per W. Skamarock, 20080606
     tendency (i,k,j) = tendency(i,k,j)                           &
               -  msft(i,j)*msft(i,j)/(ztop-zs(i,j))*( rdy*( fqy (i,k,j+1)-fqy (i,k,j)     &
                                 +fqyl(i,k,j+1)-fqyl(i,k,j)       &
								   ))

  ENDDO
  ENDDO
  ENDDO
  
  !write(6,*) 'advect_scalar_pd_fm maxval(abs(tendency) ) :', maxval(abs(tendency(ids:ide-1,kts:kte-1,jds:jde-1)) ) 
  

  IF(tenddec) THEN
  DO j = j_start, j_end
  DO k = kts, ktf
  DO i = i_start, i_end

     h_tendency (i,k,j) = h_tendency (i,k,j)                      &
               -  msft(i,j)*msft(i,j)/(ztop-zs(i,j))*( rdy*( fqy (i,k,j+1)-fqy (i,k,j)     &
                                +fqyl(i,k,j+1)-fqyl(i,k,j)    &
								  ) )

  ENDDO
  ENDDO
  ENDDO
  END IF

END SUBROUTINE advect_scalar_pd_fm



SUBROUTINE advect_scalar_mono_fm   ( field, field_old, tendency,    &
                                h_tendency, z_tendency,        &
                                u, v, sigmdot,                   &
                                c1, c2,                        &
                                mut, mub, mu_old,              &
                                time_step, config_flags,       &
                                tenddec,                       &
                                msfux, msfuy, msfvx, msfvy,    &
                                msftx, msfty,                  &
								msft,   &
                                fzm, fzp, ztop,zs,                       &
                                rdx, rdy, rdzw, dt,            &
                                ids, ide, jds, jde, kds, kde,  &
                                ims, ime, jms, jme, kms, kme,  &
                                its, ite, jts, jte, kts, kte  )

!  this is a first cut at a positive definite advection option
!  for scalars in WRF.  This version is memory intensive ->
!  we save 3d arrays of x, y and z both high and low order fluxes
!  (six in all).  Alternatively, we could sweep in a direction
!  and lower the cost considerably.

!  uses the Smolarkiewicz MWR 1989 approach, with addition of first-order
!  fluxes initially

!  WCS, 3 December 2002, 24 February 2003

   IMPLICIT NONE
   
   ! Input data
   
   TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

   LOGICAL ,                 INTENT(IN   ) :: tenddec  ! tendency flag

   INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: field,     &
                                                                      field_old, &
                                                                      u,    &
                                                                      v,    &
                                                                      sigmdot

   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: mut, mub, mu_old
   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) :: tendency
   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(  OUT) :: h_tendency, z_tendency

   REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msfux,  &
                                                                    msfuy,  &
                                                                    msfvx,  &
                                                                    msfvy,  &
                                                                    msftx,  &
                                                                    msfty,msft,zs  

   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: fzm,  &
                                                                  fzp,  &
                                                                  rdzw, &
                                                                  c1,   &
                                                                  c2

   REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy,  &
                                                                  dt,ztop 
   INTEGER ,                                     INTENT(IN   ) :: time_step

   ! Local data
   
   INTEGER :: i, j, k, itf, jtf, ktf
   INTEGER :: i_start, i_end, j_start, j_end
   INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
   INTEGER :: jmin, jmax, jp, jm, imin, imax

   REAL    :: mrdx, mrdy, ub, vb, uw, vw, mu

!  storage for high and low order fluxes

 
   REAL,  DIMENSION( its-2:ite+2, kts:kte, jts-2:jte+2  ) :: fqx, fqy, fqz
   REAL,  DIMENSION( its-2:ite+2, kts:kte, jts-2:jte+2  ) :: fqxl, fqyl, fqzl
   REAL,  DIMENSION( its-2:ite+2, kts:kte, jts-2:jte+2  ) :: qmin, qmax
   REAL,  DIMENSION( its-2:ite+2, kts:kte, jts-2:jte+2  ) :: scale_in, scale_out
   REAL :: ph_upwind


   INTEGER :: horz_order, vert_order
   
   LOGICAL :: degrade_xs, degrade_ys
   LOGICAL :: degrade_xe, degrade_ye

   INTEGER :: jp1, jp0, jtmp

   
   
   
   
   REAL :: flux_out, ph_low, flux_in, ph_hi, scale  
 
   REAL, PARAMETER :: eps=1.e-20


! definition of flux operators, 3rd, 4th, 5th or 6th order

   REAL    :: flux3, flux4, flux5, flux6, flux_upwind
   REAL    :: q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua, vel, cr

      flux4(q_im2, q_im1, q_i, q_ip1, ua) =                     &
            (7./12.)*(q_i + q_im1) - (1./12.)*(q_ip1 + q_im2)

      flux3(q_im2, q_im1, q_i, q_ip1, ua) =                     &
           flux4(q_im2, q_im1, q_i, q_ip1, ua) +                &
           sign(1,time_step)*sign(1.,ua)*(1./12.)*((q_ip1 - q_im2)-3.*(q_i-q_im1))

      flux6(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua) =       &
            (37./60.)*(q_i+q_im1) - (2./15.)*(q_ip1+q_im2)      &
            +(1./60.)*(q_ip2+q_im3)

      flux5(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua) =       &
           flux6(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua)    &
            -sign(1,time_step)*sign(1.,ua)*(1./60.)*(           &
              (q_ip2-q_im3)-5.*(q_ip1-q_im2)+10.*(q_i-q_im1) )

      flux_upwind(q_im1, q_i, cr ) = 0.5*min( 1.0,(cr+abs(cr)))*q_im1 &
                                    +0.5*max(-1.0,(cr-abs(cr)))*q_i

!      flux_upwind(q_im1, q_i, cr ) = 0.5*(1.+sign(1.,cr))*q_im1 &
!                                    +0.5*(1.-sign(1.,cr))*q_i
!      flux_upwind(q_im1, q_i, cr ) = 0.

    REAL     :: dx,dy,dz
	
	
!  LOGICAL, PARAMETER :: pd_limit = .true.
    LOGICAL, PARAMETER :: mono_limit = .true.
! set order for the advection schemes

!  write(6,*) ' in pd advection routine '

    ! Empty arrays just in case:


  ktf=MIN(kte,kde-1)
  horz_order = config_flags%h_sca_adv_order
  vert_order = config_flags%v_sca_adv_order


  do j=jts-2,jte+2
  do k=kts,kte
  do i=its-2,ite+2
    qmin(i,k,j) = field_old(i,k,j)
    qmax(i,k,j) = field_old(i,k,j)
    scale_in(i,k,j) = 1.
    scale_out(i,k,j) = 1.
    fqx(i,k,j) = 0.
    fqy(i,k,j) = 0.
    fqz(i,k,j) = 0.
    fqxl(i,k,j) = 0.
    fqyl(i,k,j) = 0.
    fqzl(i,k,j) = 0.
  enddo
  enddo
  enddo









!  determine boundary mods for flux operators
!  We degrade the flux operators from 3rd/4th order
!   to second order one gridpoint in from the boundaries for
!   all boundary conditions except periodic and symmetry - these
!   conditions have boundary zone data fill for correct application
!   of the higher order flux stencils

   degrade_xs = .true.
   degrade_xe = .true.
   degrade_ys = .true.
   degrade_ye = .true.

!  begin with horizontal flux divergence
!  here is the choice of flux operators


  horizontal_order_test : IF( horz_order == 6 ) THEN

 

   ELSE IF( horz_order == 2 ) THEN

   IF( config_flags%periodic_x   .or. &
       config_flags%symmetric_xs .or. &
       (its > ids+1)                ) degrade_xs = .false.
   IF( config_flags%periodic_x   .or. &
       config_flags%symmetric_xe .or. &
       (ite < ide-2)                ) degrade_xe = .false.
   IF( config_flags%periodic_y   .or. &
       config_flags%symmetric_ys .or. &
       (jts > jds+1)                ) degrade_ys = .false.
   IF( config_flags%periodic_y   .or. &
       config_flags%symmetric_ye .or. &
       (jte < jde-2)                ) degrade_ye = .false.

!--  y flux compute; these bounds are for periodic and sym b.c.

      ktf=MIN(kte,kde-1)
      i_start = its-1
      i_end   = MIN(ite,ide-1)+1
      j_start = jts-1
      j_end   = MIN(jte,jde-1)+1

!--  modify loop bounds if open or specified

      IF(degrade_xs) i_start = its
      IF(degrade_xe) i_end   = MIN(ite,ide-1)
      IF(degrade_ys) j_start = MAX(jts,jds+1)
      IF(degrade_ye) j_end = MIN(jte,jde-2)

!  compute fluxes, 2nd order, y flux

      DO j = j_start, j_end+1
        DO k=kts,ktf
        DO i = i_start, i_end
 
 
           dy = 2./(msftx(i,j)+msftx(i,j-1))/rdy  ! ADT eqn 48 d/dy  ??? 
        !  mu = 0.5*((c1(k)*mut(i,j)+c2(k))+(c1(k)*mut(i,j-1)+c2(k)))
          vel =  v(i,k,j)
          cr = vel*dt/dy       !/mu
        !  fqyl(i,k,j) = mu*(dy/dt)*flux_upwind(field_old(i,k,j-1), field_old(i,k,j  ), cr)

 
		            fqyl(i,k,j) =  (dy/dt)*flux_upwind((ztop-zs(i,j-1))*field(i,k,j-1)/msft(i,j-1), (ztop-zs(i,j))*field(i,k,j)/msft(i,j), cr)
		   
		   fqy(i,k,j) =  0.5*v(i,k,j)*          &
                  ((ztop-zs(i,j))*field(i,k,j)/msft(i,j) + (ztop-zs(i,j-1))*field(i,k,j-1)/msft(i,j-1))				  
				  
		       fqy(i,k,j) = fqy(i,k,j) - fqyl(i,k,j)
			   
          if(cr.gt. 0) then
             qmax(i,k,j)  = amax1(qmax(i,k,j),field_old(i,k,j-1))
             qmin(i,k,j)  = amin1(qmin(i,k,j),field_old(i,k,j-1))
          else
             qmax(i,k,j-1)  = amax1(qmax(i,k,j-1),field_old(i,k,j))
             qmin(i,k,j-1)  = amin1(qmin(i,k,j-1),field_old(i,k,j))
          end if

		  
			   
        ENDDO
        ENDDO
      ENDDO

!  next, x flux



!   20230818 
	!  j_start = jts 
    ! j_end   = MIN(jte,jde-1) 
      IF(degrade_xs) i_start = its
      IF(degrade_xe) i_end   = MIN(ite,ide-1)
      IF(degrade_ys) j_start = jts  
      IF(degrade_ye) j_end = MIN(jte,jde-1)  
	  


      DO j = j_start, j_end
        DO k=kts,ktf
        DO i = i_start, i_end+1
            dx = 2./(msfty(i,j)+msfty(i-1,j))/rdx ! ADT eqn 48 d/dx   ???
          !  mu = 0.5*((c1(k)*mut(i,j)+c2(k))+(c1(k)*mut(i-1,j)+c2(k)))
           vel =  u(i,k,j)
         cr = vel*dt/dx !/mu
          fqxl(i,k,j) =  (dx/dt)*flux_upwind((ztop-zs(i-1,j))*field_old(i-1,k,j)/msft(i-1,j), (ztop-zs(i,j))*field_old(i,k,j  )/msft(i,j), cr)
		  
		! fqyl(i,k,j) =  (dy/dt)*flux_upwind((ztop-zs(i,j-1))*field(i,k,j-1)/msft(i,j-1), (ztop-zs(i,j))*field(i,k,j)/msft(i,j), cr)  
		  
		  
          !  fqx( i,k,j ) = 0.5*ru(i,k,j)*          &
           !         (field(i,k,j)+field(i-1,k,j))
				  
            fqx( i,k,j ) = 0.5*u(i,k,j)*          &
                  ((ztop-zs(i,j))*field(i,k,j)/msft(i,j) + (ztop-zs(i-1,j))*field(i-1,k,j)/msft(i-1,j))				  
				  
				  

          fqx(i,k,j) = fqx(i,k,j) - fqxl(i,k,j)
		  
		  
               if(cr.gt. 0) then
                  qmax(i,k,j)  = amax1(qmax(i,k,j),field_old(i-1,k,j))
                  qmin(i,k,j)  = amin1(qmin(i,k,j),field_old(i-1,k,j))
                else
                  qmax(i-1,k,j)  = amax1(qmax(i-1,k,j),field_old(i,k,j))
                  qmin(i-1,k,j)  = amin1(qmin(i-1,k,j),field_old(i,k,j))
                end if

				
		  
		  
        ENDDO
        ENDDO
      ENDDO

!--- end of 2nd order horizontal flux calculation

   ELSE

      WRITE ( wrf_err_message , * ) 'module_advect: advect_scalar_pd, h_order not known ',horz_order
      CALL wrf_error_fatal ( TRIM( wrf_err_message ) )

   ENDIF horizontal_order_test






















!  pick up the rest of the horizontal radiation boundary conditions.
!  (these are the computations that don't require 'cb'.
!  first, set to index ranges

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)

!  compute x (u) conditions for v, w, or scalar

   IF( (config_flags%open_xs) .and. (its == ids) ) THEN

       DO j = j_start, j_end
       DO k = kts, ktf
         ub = MIN( 0.5*(u(its,k,j)+u(its+1,k,j)), 0. )
         tendency(its,k,j) = tendency(its,k,j)                     &
             !  - rdx*(                                             &
             !          ub*(   field_old(its+1,k,j)                 &
             !               - field_old(its  ,k,j)   ) +           &
             !          field(its,k,j)*(ru(its+1,k,j)-ru(its,k,j))  &
             !                                                   )
                 - (ztop-zs(its,j))/msft(its,j)*rdx*(                                             &
                       ub*(   field_old(its+1,k,j)                 &
                            - field_old(its  ,k,j)   ) +           &
                       field_old(its,k,j)*(u(its+1,k,j)-u(its,k,j))  &       ! modified 
                                                                ) * msft(its,j)*msft(its,j)/(ztop-zs(its,j))  																	
																
																
       ENDDO
       ENDDO

   ENDIF

   IF( (config_flags%open_xe) .and. (ite == ide) ) THEN

       DO j = j_start, j_end
       DO k = kts, ktf
         ub = MAX( 0.5*(u(ite-1,k,j)+u(ite,k,j)), 0. )
         tendency(i_end,k,j) = tendency(i_end,k,j)                   &
        !       - rdx*(                                               &
         !              ub*(  field_old(i_end  ,k,j)                  &
         !                  - field_old(i_end-1,k,j) ) +              &
         !              field(i_end,k,j)*(ru(ite,k,j)-ru(ite-1,k,j))  &
         !                                                           )
																	
                - (ztop-zs(i_end,j))/msft(i_end,j)*rdx*(                                               &
                       ub*(  field_old(i_end  ,k,j)                  &
                           - field_old(i_end-1,k,j) ) +              &
                       field_old(i_end,k,j)*(u(ite,k,j)-u(ite-1,k,j))  &   ! modified  
                                                                    )* msft(i_end,j)*msft(i_end,j)/(ztop-zs(i_end,j)) 																		
																	
       ENDDO
       ENDDO

   ENDIF

   IF( (config_flags%open_ys) .and. (jts == jds) ) THEN

       DO i = i_start, i_end
       DO k = kts, ktf
         vb = MIN( 0.5*(v(i,k,jts)+v(i,k,jts+1)), 0. )
         tendency(i,k,jts) = tendency(i,k,jts)                     &
               - (ztop-zs(i,jts))/msft(i,jts)*rdy*(                                             &
                       vb*(  field_old(i,k,jts+1)                  &
                           - field_old(i,k,jts  ) ) +              &
                       field_old(i,k,jts)*(v(i,k,jts+1)-v(i,k,jts))  &     ! modified  
                                                               )* msft(i,jts)*msft(i,jts)/(ztop-zs(i,jts)) 


       ENDDO
       ENDDO

   ENDIF

   IF( (config_flags%open_ye) .and. (jte == jde)) THEN

       DO i = i_start, i_end
       DO k = kts, ktf
         vb = MAX( 0.5*(v(i,k,jte-1)+v(i,k,jte)), 0. )
         tendency(i,k,j_end) = tendency(i,k,j_end)                   &
               -  (ztop-zs(i,j_end))/msft(i,j_end)*rdy*(                                               &
                       vb*(   field_old(i,k,j_end  )                 &
                            - field_old(i,k,j_end-1) ) +             &
                       field_old(i,k,j_end)*(v(i,k,jte)-v(i,k,jte-1))  &    ! modified   
                                                                    )* msft(i,j_end)*msft(i,j_end)/(ztop-zs(i,j_end)) 
       ENDDO
       ENDDO

   ENDIF

 !  IF( (config_flags%polar) .and. (jts == jds) ) THEN

   

  ! ENDIF

  ! IF( (config_flags%polar) .and. (jte == jde)) THEN
 

  ! ENDIF




  ! write(6,*) 'advect_scalar_pd_fm maxval(abs(tendency) ) after openbdy  :', maxval(abs(tendency) ) 














!-------------------- vertical advection

!-- loop bounds for periodic or sym conditions

      i_start = its-1
      i_end   = MIN(ite,ide-1)+1
      j_start = jts-1
      j_end   = MIN(jte,jde-1)+1

!-- loop bounds for open or specified conditions

    IF(degrade_xs) i_start = MAX(its-1,ids)
    IF(degrade_xe) i_end   = MIN(ite+1,ide-1)
    IF(degrade_ys) j_start = MAX(jts-1,jds)
    IF(degrade_ye) j_end   = MIN(jte+1,jde-1)

    vert_order_test : IF (vert_order == 6) THEN    

 
   ELSE IF (vert_order == 2) THEN    

      DO j = j_start, j_end

         DO i = i_start, i_end
           fqz(i,1,j)  = 0.
           fqzl(i,1,j) = 0.
           fqz(i,kde,j)  = 0.
           fqzl(i,kde,j) = 0.
         ENDDO

         DO k=kts+1,ktf
         DO i = i_start, i_end

            dz = 2./(rdzw(k)+rdzw(k-1))
         !  mu = 0.5*((c1(k)*mut(i,j)+c2(k))+(c1(k)*mut(i,j)+c2(k)))
            vel = sigmdot(i,k,j)  !rom(i,k,j)
           cr = vel*dt/dz  !/mu
         !  fqzl(i,k,j) = mu*(dz/dt)*flux_upwind(field_old(i,k-1,j), field_old(i,k,j  ), cr)
        !   fqz(i,k,j)=rom(i,k,j)*(fzm(k)*field(i,k,j)+fzp(k)*field(i,k-1,j))
         
		  fqzl(i,k,j) = (dz/dt)*flux_upwind(field_old(i,k-1,j), field_old(i,k,j  ), cr)  *(ztop-zs(i,j))/msft(i,j)   

		
	   fqz(i,k,j) = sigmdot(i,k,j) *(fzm(k)*field(i,k,j)+fzp(k)*field(i,k-1,j))*(ztop-zs(i,j))/msft(i,j)    !这里不是m**2
		  fqz(i,k,j) = fqz(i,k,j) - fqzl(i,k,j)  
		  
		  
		 if(cr.gt. 0) then
             qmax(i,k,j)  = amax1(qmax(i,k,j),field_old(i,k-1,j))
             qmin(i,k,j)  = amin1(qmin(i,k,j),field_old(i,k-1,j))
          else
             qmax(i,k-1,j)  = amax1(qmax(i,k-1,j),field_old(i,k,j))
             qmin(i,k-1,j)  = amin1(qmin(i,k-1,j),field_old(i,k,j))
          end if

		  
		  

        ENDDO
        ENDDO

      ENDDO

   ELSE

      WRITE (wrf_err_message,*) ' advect_scalar_pd, v_order not known ',vert_order
      CALL wrf_error_fatal ( wrf_err_message )

   ENDIF vert_order_test
   
   
   
   
   
   
   
   
   
   
   
   
   
   

! IF (pd_limit) THEN

! positive definite filter


   IF (mono_limit) THEN

! montonic filter

   i_start = its-1
   i_end   = MIN(ite,ide-1)+1
   j_start = jts-1
   j_end   = MIN(jte,jde-1)+1

!-- loop bounds for open or specified conditions

   IF(degrade_xs) i_start = MAX(its-1,ids)
   IF(degrade_xe) i_end   = MIN(ite+1,ide-1)
   IF(degrade_ys) j_start = MAX(jts-1,jds)
   IF(degrade_ye) j_end   = MIN(jte+1,jde-1)

   IF(config_flags%specified .or. config_flags%nested) THEN
     IF (degrade_xs) i_start = MAX(its-1,ids+1)
     IF (degrade_xe) i_end   = MIN(ite+1,ide-2)
     IF (degrade_ys) j_start = MAX(jts-1,jds+1)
     IF (degrade_ye) j_end   = MIN(jte+1,jde-2)
   END IF

   IF(config_flags%open_xs) THEN
     IF (degrade_xs) i_start = MAX(its-1,ids+1)
   END IF
   IF(config_flags%open_xe) THEN
     IF (degrade_xe) i_end   = MIN(ite+1,ide-2)
   END IF
   IF(config_flags%open_ys) THEN
     IF (degrade_ys) j_start = MAX(jts-1,jds+1)
   END IF
   IF(config_flags%open_ye) THEN
     IF (degrade_ye) j_end   = MIN(jte+1,jde-2)
   END IF
   ! ADT note:
   ! We don't want to change j_start and j_end
   ! for polar BC's since we want to calculate
   ! fluxes for directions other than y at the
   ! edge

!-- here is the limiter...




   DO j=j_start, j_end
   DO k=kts, ktf
   DO i=i_start, i_end

     ph_upwind = field_old(i,k,j) &
                - dt*msft(i,j)*msft(i,j)/(ztop-zs(i,j)) *( (                &
                          rdx*(fqxl(i+1,k,j)-fqxl(i,k,j)) +     &
                          rdy*(fqyl(i,k,j+1)-fqyl(i,k,j))  )    &
                         +rdzw(k)*(fqzl(i,k+1,j)-fqzl(i,k,j)) )

     flux_in = -dt*msft(i,j)*msft(i,j)/(ztop-zs(i,j)) *( (                   &
                               rdx*(  min(0.,fqx (i+1,k,j))      &
                                     -max(0.,fqx (i  ,k,j)) )    &
                              +rdy*(  min(0.,fqy (i,k,j+1))      &
                                     -max(0.,fqy (i,k,j  )) ) )  &
               + rdzw(k)*(  max(0.,fqz (i,k+1,j))      &
                                     -min(0.,fqz (i,k  ,j)) )   )
									 
									 
									 
									 

     ph_hi =  qmax(i,k,j) - ph_upwind
	 
	 
     IF( flux_in .gt. ph_hi ) scale_in(i,k,j) = max(0.,ph_hi/(flux_in+eps))


     flux_out = dt*msft(i,j)*msft(i,j)/(ztop-zs(i,j)) * (  (                    &
                                rdx*(  max(0.,fqx (i+1,k,j))      &
                                      -min(0.,fqx (i  ,k,j)) )    &
                               +rdy*(  max(0.,fqy (i,k,j+1))      &
                                      -min(0.,fqy (i,k,j  )) ) )  &
                + rdzw(k)*(  min(0.,fqz (i,k+1,j))      &
                                      -max(0.,fqz (i,k  ,j)) )   )

     ph_low = ph_upwind - qmin(i,k,j)
     IF( flux_out .gt. ph_low ) scale_out(i,k,j) = max(0.,ph_low/(flux_out+eps))

   ENDDO
   ENDDO
   ENDDO

   DO j=j_start, j_end
   DO k=kts, ktf
   DO i=i_start, i_end+1
       IF( fqx (i,k,j) .gt. 0.) then
         fqx(i,k,j) = min(scale_in(i,k,j),scale_out(i-1,k,j))*fqx(i,k,j)
       ELSE
         fqx(i,k,j) = min(scale_out(i,k,j),scale_in(i-1,k,j))*fqx(i,k,j)
       ENDIF
   ENDDO
   ENDDO
   ENDDO

   DO j=j_start, j_end+1
   DO k=kts, ktf
   DO i=i_start, i_end
       IF( fqy (i,k,j) .gt. 0.) then
         fqy(i,k,j) = min(scale_in(i,k,j),scale_out(i,k,j-1))*fqy(i,k,j)
       ELSE
         fqy(i,k,j) = min(scale_out(i,k,j),scale_in(i,k,j-1))*fqy(i,k,j)
       ENDIF
   ENDDO
   ENDDO
   ENDDO

   DO j=j_start, j_end
   DO k=kts+1, ktf
   DO i=i_start, i_end
       IF( fqz (i,k,j) .lt. 0.) then
         fqz(i,k,j) = min(scale_in(i,k,j),scale_out(i,k-1,j))*fqz(i,k,j)
       ELSE
         fqz(i,k,j) = min(scale_out(i,k,j),scale_in(i,k-1,j))*fqz(i,k,j)
       ENDIF
   ENDDO
   ENDDO
   ENDDO

   END IF








   
! add in the mono-limited flux divergence
! we need to fix this for open b.c set ***********

  i_start = its
  i_end   = MIN(ite,ide-1)
  j_start = jts
  j_end   = MIN(jte,jde-1)

  DO j = j_start, j_end
  DO k = kts, ktf
!DEC$ vector always
  DO i = i_start, i_end

     tendency (i,k,j) = tendency(i,k,j)                           &
                            -rdzw(k)* msft(i,j)*msft(i,j)/(ztop-zs(i,j))*( fqz (i,k+1,j)-fqz (i,k,j)  &
                                      +fqzl(i,k+1,j)-fqzl(i,k,j)  &
									  )

  ENDDO
  ENDDO
  ENDDO
  
   
  IF(tenddec) THEN
  DO j = j_start, j_end
  DO k = kts, ktf
  DO i = i_start, i_end

     z_tendency (i,k,j) = 0. -rdzw(k)* msft(i,j)*msft(i,j)/(ztop-zs(i,j))*( fqz (i,k+1,j)-fqz (i,k,j)  &
                                     +fqzl(i,k+1,j)-fqzl(i,k,j)    &
									  )

  ENDDO
  ENDDO
  ENDDO
  END IF

! x flux divergence
!
  IF(degrade_xs) i_start = MAX(its,ids+1)
  IF(degrade_xe) i_end   = MIN(ite,ide-2)

  DO j = j_start, j_end
  DO k = kts, ktf
!DEC$ vector always  
  DO i = i_start, i_end

     ! Un-"canceled" map scale factor, ADT Eq. 48
     tendency (i,k,j) = tendency(i,k,j)                           &
               - msft(i,j)*msft(i,j)/(ztop-zs(i,j))*( rdx*( fqx (i+1,k,j)-fqx (i,k,j)     &
                                   +fqxl(i+1,k,j)-fqxl(i,k,j)   &
								   ) )

  ENDDO
  ENDDO
  ENDDO
  
 
  IF(tenddec) THEN
  DO j = j_start, j_end
  DO k = kts, ktf
  DO i = i_start, i_end

     h_tendency (i,k,j) = 0.                                      &
               - msft(i,j)*msft(i,j)/(ztop-zs(i,j))*( rdx*( fqx (i+1,k,j)-fqx (i,k,j)     &
                                   +fqxl(i+1,k,j)-fqxl(i,k,j)     &
								   ))

  ENDDO
  ENDDO
  ENDDO
  END IF

! y flux divergence
!
  i_start = its
  i_end   = MIN(ite,ide-1)
  IF(degrade_ys) j_start = MAX(jts,jds+1)
  IF(degrade_ye) j_end   = MIN(jte,jde-2)

  DO j = j_start, j_end
  DO k = kts, ktf
!DEC$ vector always
  DO i = i_start, i_end

     ! Un-"canceled" map scale factor, ADT Eq. 48
     ! It is correct to use msftx (and not msfty), per W. Skamarock, 20080606
     tendency (i,k,j) = tendency(i,k,j)                           &
               -  msft(i,j)*msft(i,j)/(ztop-zs(i,j))*( rdy*( fqy (i,k,j+1)-fqy (i,k,j)     &
                                 +fqyl(i,k,j+1)-fqyl(i,k,j)       &
								   ))

  ENDDO
  ENDDO
  ENDDO
  
  !write(6,*) 'advect_scalar_pd_fm maxval(abs(tendency) ) :', maxval(abs(tendency(ids:ide-1,kts:kte-1,jds:jde-1)) ) 
  

  IF(tenddec) THEN
  DO j = j_start, j_end
  DO k = kts, ktf
  DO i = i_start, i_end

     h_tendency (i,k,j) = h_tendency (i,k,j)                      &
               -  msft(i,j)*msft(i,j)/(ztop-zs(i,j))*( rdy*( fqy (i,k,j+1)-fqy (i,k,j)     &
                                +fqyl(i,k,j+1)-fqyl(i,k,j)    &
								  ) )

  ENDDO
  ENDDO
  ENDDO
  END IF

END SUBROUTINE advect_scalar_mono_fm

	  
SUBROUTINE calcu_divg_p( var, divg,pfm_base,rhofm_base,                 &
                        u, v, sgmdot,                       &
                     config_flags,                  &
                        msft,zh,zs,fzm,fzp,dnw,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   ! Input data
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var,u,v,    &
                                                                      pfm_base,rhofm_base, &
                                                                      sgmdot
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      !   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: mut
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: divg

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh,fzm,fzp,dnw


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   ! Local data
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx
         
        
       !  REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: zzh,dzzh
	   
	   
	   
	   
	   
		 
		    REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) ::  varw
		 
		 
		 
		 


      REAL , DIMENSION( kms:kme ) :: sgmh,dsgm
       REAL , DIMENSION( ims:ime , jms:jme )  ::  rmsft,rmm
	   
	   REAL :: ub,vb 
				 
				 
		!20230816		 
				 
	   REAL,  DIMENSION( its-1:ite+2, kts:kte, jts-1:jte+2  ) :: fqx, fqy, fqz
	   
	       REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) ::  field

   LOGICAL :: degrade_xs, degrade_ys
   LOGICAL :: degrade_xe, degrade_ye
   
   
		   

    LOGICAL :: specified




   degrade_xs = .true.
   degrade_xe = .true.
   degrade_ys = .true.
   degrade_ye = .true.


   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

   
     ktf=MIN(kte,kde-1)
 

	  DO k=kts,kde-1
          sgmh(k) = 0.5*(zh(k)+zh(k+1))
          dsgm(k) = zh(k+1) - zh(k) 
      ENDDO


 



	  

   ! 20230228  for rho couple scalar 
      CALL set_physical_bc3d( var, 'p', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              its, ite, jts, jte, kts, kte,       &
                              its, ite, jts, jte, kts, kte )
	  
	  

     divg = 0.0 
	 
	 
! boundary	 
	 
!  ---- x direction -----------
 

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)


      IF ( .NOT. config_flags%periodic_x ) THEN
        IF ( config_flags%open_xs .or. specified ) i_start = MAX(ids+1,its)
        IF ( config_flags%open_xe .or. specified ) i_end   = MIN(ide-2,ite)
      ENDIF

      
	  
	  
	IF ( config_flags%open_xs  .and.  its == ids ) THEN 
		
	!	i = ids  !  用于方法1 
		
      DO  j = j_start, j_end
      DO k=kts,ktf
   
 					 
	! 参考wrf	,这一方法的使用，改进非常大。避免了hill2km_x试验的西边界问题。			 
	         ub = MIN( 0.5*(u(its,k,j)+u(its+1,k,j)), 0. )
		 
         divg(its,k,j) =                      &
               - (ztop-zs(its,j))/msft(its,j)*rdx*(                                             &
                       ub*(   var(its+1,k,j)                 &
                            - var(its  ,k,j)   ) +           &
                       var(its,k,j)*(u(its+1,k,j)-u(its,k,j))  &       ! modified 
                                                                )				 
					 
      ENDDO
      ENDDO
 
	  
	ENDIF	  
!	write(6,*) '1 divg(1,1,1) = ', divg(1,1,1)  

	  
	    IF ( config_flags%open_xe  .and.  ite == ide ) THEN 
		
	 !	i = ite-1
		
		 i_end   = MIN(ite,ide-1)
		 
		
      DO  j = j_start, j_end
      DO k=kts,ktf
 
  
     
	     ub = MAX( 0.5*(u(ite-1,k,j)+u(ite,k,j)), 0. )	 
		 
         divg(i_end,k,j) =                  &
               - (ztop-zs(i_end,j))/msft(i_end,j)*rdx*(                                               &
                       ub*(  var(i_end  ,k,j)                  &
                           - var(i_end-1,k,j) ) +              &
                       var(i_end,k,j)*(u(ite,k,j)-u(ite-1,k,j))  &   ! modified  
                                                                    )				 
					 
					 
      ENDDO
      ENDDO
  

	ENDIF	
	  
	  
	!  write(6,*) '2 divg(1,1,1) = ', divg(1,1,1)  
	  

!  ---- y direction -----------

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
      ! Polar boundary conditions are like open or specified
      IF ( config_flags%open_ys .or. specified .or. config_flags%polar ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified .or. config_flags%polar ) j_end   = MIN(jde-2,jte) 	  

       
	  
	    IF ( config_flags%open_ys  .and.  jts == jds ) THEN 
		
	!	j = jds 
		
      DO k=kts,ktf
      DO  i = i_start, i_end		

	
	      vb = MIN( 0.5*(v(i,k,jts)+v(i,k,jts+1)), 0. )
		  
        divg(i,k,jts) =    divg(i,k,jts)                &
               - (ztop-zs(i,jts))/msft(i,jts)*rdy*(                                             &
                       vb*(  var(i,k,jts+1)                  &
                           - var(i,k,jts  ) ) +              &
                       var(i,k,jts)*(v(i,k,jts+1)-v(i,k,jts))  &     ! modified  
                                                               )
					 
      ENDDO
      ENDDO
	    ENDIF	
		
 
	  
	    IF ( config_flags%open_ye  .and.  jte == jde ) THEN 
		
	 ! 	j = jte-1
		
		j_end   = MIN(jte,jde-1)
		
      DO k=kts,ktf
      DO  i = i_start, i_end		
		
	 
	
		  vb = MAX( 0.5*(v(i,k,jte-1)+v(i,k,jte)), 0. ) 
		 
         divg(i,k,j_end) =    divg(i,k,j_end)               &
               -  (ztop-zs(i,j_end))/msft(i,j_end)*rdy*(                                               &
                       vb*(   var(i,k,j_end  )                 &
                            - var(i,k,j_end-1) ) +             &
                       var(i,k,j_end)*(v(i,k,jte)-v(i,k,jte-1))  &    ! modified   
                                                                    )
					 
      ENDDO
      ENDDO
	    ENDIF
		
		
     !     write(6,*) '4 divg(1,1,1) = ', divg(1,1,1) 







!!!!!!!!!!!!!!!!!!!!!!!!!!!	 
	 
      field  = var  	 
	 
	 
!  ---- z  direction -----------	




! fqz

!-------------------- vertical advection




	  
!-- loop bounds for periodic or sym conditions

      i_start = its-1
      i_end   = MIN(ite,ide-1)+1
      j_start = jts-1
      j_end   = MIN(jte,jde-1)+1

!-- loop bounds for open or specified conditions

    IF(degrade_xs) i_start = MAX(its-1,ids)
    IF(degrade_xe) i_end   = MIN(ite+1,ide-1)
    IF(degrade_ys) j_start = MAX(jts-1,jds)
    IF(degrade_ye) j_end   = MIN(jte+1,jde-1)

    

      DO j = j_start, j_end

         DO i = i_start, i_end
           fqz(i,1,j)  = 0.
        !   fqzl(i,1,j) = 0.
           fqz(i,kde,j)  = 0.
          ! fqzl(i,kde,j) = 0.
         ENDDO

         DO k=kts+1,ktf
         DO i = i_start, i_end


		   
	   fqz(i,k,j) = sgmdot(i,k,j) *(fzm(k)*field(i,k,j)+fzp(k)*field(i,k-1,j))*(ztop-zs(i,j))/msft(i,j)    !这里不是m**2
		   

        ENDDO
        ENDDO

      ENDDO

    
   



! fqz end 
  
	  
	  	call average_vertical_pw( varw, var,pfm_base,rhofm_base,                 &
                     config_flags,                  &
                         fzm,fzp,dnw,                    &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )	



  !  计划拆分为简单公式，处理变量从半层插值到整层，及边界处理问题
  
  
        i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
	  
      DO  j = j_start, j_end
      DO  i = i_start, i_end
        rmsft(i,j)  = 1./msft(i,j) 
        rmm(i,j) = 1./(msft(i,j)*msft(i,j) )
      ENDDO
      ENDDO             

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end


           ! 20230817测试不用fqz，用原来方案。
             divg(i,k,j) = divg(i,k,j) -    ( varw(i,k+1,j)* rmm(i,j) * sgmdot(i,k+1,j) - varw(i,k,j)* rmm(i,j) * sgmdot(i,k,j)   ) /dsgm(k) * (ztop-zs(i,j))

            !      divg(i,k,j) = divg(i,k,j) -rdzw(k)*( fqz (i,k+1,j)-fqz (i,k,j) )
									  
			 !       divg(i,k,j) = divg(i,k,j) - ( fqz (i,k+1,j)-fqz (i,k,j) ) /dsgm(k)					  
									  
      ENDDO
      ENDDO
      ENDDO
	 
 
  
		 
        
   IF( config_flags%periodic_x   .or. &
       config_flags%symmetric_xs .or. &
       (its > ids+1)                ) degrade_xs = .false.
   IF( config_flags%periodic_x   .or. &
       config_flags%symmetric_xe .or. &
       (ite < ide-2)                ) degrade_xe = .false.
   IF( config_flags%periodic_y   .or. &
       config_flags%symmetric_ys .or. &
       (jts > jds+1)                ) degrade_ys = .false.
   IF( config_flags%periodic_y   .or. &
       config_flags%symmetric_ye .or. &
       (jte < jde-2)                ) degrade_ye = .false.

!--  y flux compute; these bounds are for periodic and sym b.c.

      ktf=MIN(kte,kde-1)
      i_start = its-1
      i_end   = MIN(ite,ide-1)+1
      j_start = jts-1
      j_end   = MIN(jte,jde-1)+1

!--  modify loop bounds if open or specified

      IF(degrade_xs) i_start = its
      IF(degrade_xe) i_end   = MIN(ite,ide-1)
      IF(degrade_ys) j_start = MAX(jts,jds+1)
      IF(degrade_ye) j_end = MIN(jte,jde-2)


      DO j = j_start, j_end+1
        DO k=kts,ktf
        DO i = i_start, i_end
 
		   
	!	   fqy(i,k,j) =  0.5*v(i,k,j)*          &
    !              ((ztop-zs(i,j))*field(i,k,j)/msft(i,j) + (ztop-zs(i,j-1))*field(i,k,j-1)/msft(i,j-1))	

				  
	 	    fqy(i,k,j) =    0.5*((ztop-zs(i,j-1))*var(i,k,j-1)*rmsft(i,j-1)  + (ztop-zs(i,j))*var(i,k,j)*rmsft(i,j) ) *v(i,k,j) 	  
		   
        ENDDO
        ENDDO
      ENDDO
	  
	  
	  ! 20230818 14:12 
	  j_start = jts 
      j_end   = MIN(jte,jde-1) 
	  
	  
      DO j = j_start, j_end
        DO k=kts,ktf
        DO i = i_start, i_end+1
 
				  
      !      fqx( i,k,j ) = 0.5*u(i,k,j)*          &
      !            ((ztop-zs(i,j))*field(i,k,j)/msft(i,j) + (ztop-zs(i-1,j))*field(i-1,k,j)/msft(i-1,j))			

		 	 fqx( i,k,j ) =  	  0.5*((ztop-zs(i-1,j))*var(i-1,k,j)/msft(i-1,j)  + (ztop-zs(i,j))*var(i,k,j)/msft(i,j) ) *u(i,k,j)   
 
        ENDDO
        ENDDO
      ENDDO
	  
	  


!!!!!!!!!!!!!!!!!!!!!!!	 
	 
	 
	 
!  ---- x direction -----------
 

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)


      IF ( .NOT. config_flags%periodic_x ) THEN
        IF ( config_flags%open_xs .or. specified ) i_start = MAX(ids+1,its)
        IF ( config_flags%open_xe .or. specified ) i_end   = MIN(ide-2,ite)
      ENDIF

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
  
    

       !  	divg(i,k,j) =  divg(i,k,j) -  (  0.5*((ztop-zs(i,j))*var(i,k,j)/msft(i,j)  + (ztop-zs(i+1,j))*var(i+1,k,j)/msft(i+1,j) ) *u(i+1,k,j) &
         !             -   0.5*((ztop-zs(i-1,j))*var(i-1,k,j)/msft(i-1,j)  + (ztop-zs(i,j))*var(i,k,j)/msft(i,j) ) *u(i,k,j) )*rdx
					 
     	divg(i,k,j) =  divg(i,k,j) -   rdx*( fqx (i+1,k,j)-fqx (i,k,j)   )					 
					 
					 
      ENDDO
      ENDDO
      ENDDO
	  
	!   write(6,*) '6 divg(1,1,1) = ', divg(1,1,1) 

!  ---- y direction -----------

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
      ! Polar boundary conditions are like open or specified
      IF ( config_flags%open_ys .or. specified .or. config_flags%polar ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified .or. config_flags%polar ) j_end   = MIN(jde-2,jte) 	  

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end

    !   	divg(i,k,j) = divg(i,k,j) -   ( 0.5*((ztop-zs(i,j))*var(i,k,j)*rmsft(i,j)  + (ztop-zs(i,j+1))*var(i,k,j+1)*rmsft(i,j+1) ) *v(i,k,j+1)    &
    !                                -   0.5*((ztop-zs(i,j-1))*var(i,k,j-1)*rmsft(i,j-1)  + (ztop-zs(i,j))*var(i,k,j)*rmsft(i,j) ) *v(i,k,j) )*rdy
					 
 		 
     	divg(i,k,j) = divg(i,k,j) -    rdy*( fqy (i,k,j+1)-fqy (i,k,j) )  					 

   

       ENDDO
      ENDDO
      ENDDO
	  
	  
     
	
		


  
        i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
         	divg(i,k,j) = divg(i,k,j)  * msft(i,j)*msft(i,j)/(ztop-zs(i,j)) 
      ENDDO
      ENDDO
      ENDDO
	  
 

      END SUBROUTINE calcu_divg_p
	  
	  
	  
	  
	  
	  
	  
	  
	  
SUBROUTINE calcu_divg_ps( var, divg,rhofm_save,                 &
                        u, v, sgmdot,                       &
                     config_flags,                  &
                        msft,zh,zs,fzm,fzp,dnw,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   ! Input data
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var,u,v,    &
                                                                       rhofm_save, &
                                                                      sgmdot
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      !   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: mut
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: divg

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh,fzm,fzp,dnw


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   ! Local data
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx
         
        
       !  REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: zzh,dzzh
	   
	   
	   
	   
	   
		 
		    REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) ::  varw
		 
		 
		 
		 


      REAL , DIMENSION( kms:kme ) :: sgmh,dsgm
       REAL , DIMENSION( ims:ime , jms:jme )  ::  rmsft,rmm
	   
	   REAL :: ub,vb 
				 
				 
		!20230816		 
				 
	   REAL,  DIMENSION( its-1:ite+2, kts:kte, jts-1:jte+2  ) :: fqx, fqy, fqz
	   
	       REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) ::  field

   LOGICAL :: degrade_xs, degrade_ys
   LOGICAL :: degrade_xe, degrade_ye
   
   
		   

    LOGICAL :: specified




   degrade_xs = .true.
   degrade_xe = .true.
   degrade_ys = .true.
   degrade_ye = .true.


   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

   
     ktf=MIN(kte,kde-1)
 

	  DO k=kts,kde-1
          sgmh(k) = 0.5*(zh(k)+zh(k+1))
          dsgm(k) = zh(k+1) - zh(k) 
      ENDDO


 



	  

   ! 20230228  for rho couple scalar 
      CALL set_physical_bc3d( var, 'p', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              its, ite, jts, jte, kts, kte,       &
                              its, ite, jts, jte, kts, kte )
	  
	  

     divg = 0.0 
	 
	 
! boundary	 
	 
!  ---- x direction -----------
 

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)


      IF ( .NOT. config_flags%periodic_x ) THEN
        IF ( config_flags%open_xs .or. specified ) i_start = MAX(ids+1,its)
        IF ( config_flags%open_xe .or. specified ) i_end   = MIN(ide-2,ite)
      ENDIF

      
	  
	  
	IF ( config_flags%open_xs  .and.  its == ids ) THEN 
		
	!	i = ids  !  用于方法1 
		
      DO  j = j_start, j_end
      DO k=kts,ktf
   
 					 
	! 参考wrf	,这一方法的使用，改进非常大。避免了hill2km_x试验的西边界问题。			 
	         ub = MIN( 0.5*(u(its,k,j)+u(its+1,k,j)), 0. )
		 
         divg(its,k,j) =                      &
               - (ztop-zs(its,j))/msft(its,j)*rdx*(                                             &
                       ub*(   var(its+1,k,j)                 &
                            - var(its  ,k,j)   ) +           &
                       var(its,k,j)*(u(its+1,k,j)-u(its,k,j))  &       ! modified 
                                                                )				 
					 
      ENDDO
      ENDDO
 
	  
	ENDIF	  
!	write(6,*) '1 divg(1,1,1) = ', divg(1,1,1)  

	  
	    IF ( config_flags%open_xe  .and.  ite == ide ) THEN 
		
	 !	i = ite-1
		
		 i_end   = MIN(ite,ide-1)
		 
		
      DO  j = j_start, j_end
      DO k=kts,ktf
 
  
     
	     ub = MAX( 0.5*(u(ite-1,k,j)+u(ite,k,j)), 0. )	 
		 
         divg(i_end,k,j) =                  &
               - (ztop-zs(i_end,j))/msft(i_end,j)*rdx*(                                               &
                       ub*(  var(i_end  ,k,j)                  &
                           - var(i_end-1,k,j) ) +              &
                       var(i_end,k,j)*(u(ite,k,j)-u(ite-1,k,j))  &   ! modified  
                                                                    )				 
					 
					 
      ENDDO
      ENDDO
  

	ENDIF	
	  
	  
 
!  ---- y direction -----------

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
      ! Polar boundary conditions are like open or specified
      IF ( config_flags%open_ys .or. specified .or. config_flags%polar ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified .or. config_flags%polar ) j_end   = MIN(jde-2,jte) 	  

       
	  
	    IF ( config_flags%open_ys  .and.  jts == jds ) THEN 
		
	!	j = jds 
		
      DO k=kts,ktf
      DO  i = i_start, i_end		

	
	      vb = MIN( 0.5*(v(i,k,jts)+v(i,k,jts+1)), 0. )
		  
        divg(i,k,jts) =    divg(i,k,jts)                &
               - (ztop-zs(i,jts))/msft(i,jts)*rdy*(                                             &
                       vb*(  var(i,k,jts+1)                  &
                           - var(i,k,jts  ) ) +              &
                       var(i,k,jts)*(v(i,k,jts+1)-v(i,k,jts))  &     ! modified  
                                                               )
					 
      ENDDO
      ENDDO
	    ENDIF	
		
  !      write(6,*) '3 divg(1,1,1) = ', divg(1,1,1) 
	  
	    IF ( config_flags%open_ye  .and.  jte == jde ) THEN 
		
	 ! 	j = jte-1
		
		j_end   = MIN(jte,jde-1)
		
      DO k=kts,ktf
      DO  i = i_start, i_end		
		
	 
	
		  vb = MAX( 0.5*(v(i,k,jte-1)+v(i,k,jte)), 0. ) 
		 
         divg(i,k,j_end) =    divg(i,k,j_end)               &
               -  (ztop-zs(i,j_end))/msft(i,j_end)*rdy*(                                               &
                       vb*(   var(i,k,j_end  )                 &
                            - var(i,k,j_end-1) ) +             &
                       var(i,k,j_end)*(v(i,k,jte)-v(i,k,jte-1))  &    ! modified   
                                                                    )
					 
      ENDDO
      ENDDO
	    ENDIF
		
		
 





!!!!!!!!!!!!!!!!!!!!!!!!!!!	 
	 
      field  = var  	 
	 
	 
!  ---- z  direction -----------	




! fqz

!-------------------- vertical advection




	  
!-- loop bounds for periodic or sym conditions

      i_start = its-1
      i_end   = MIN(ite,ide-1)+1
      j_start = jts-1
      j_end   = MIN(jte,jde-1)+1

!-- loop bounds for open or specified conditions

    IF(degrade_xs) i_start = MAX(its-1,ids)
    IF(degrade_xe) i_end   = MIN(ite+1,ide-1)
    IF(degrade_ys) j_start = MAX(jts-1,jds)
    IF(degrade_ye) j_end   = MIN(jte+1,jde-1)

    

      DO j = j_start, j_end

         DO i = i_start, i_end
           fqz(i,1,j)  = 0.
        !   fqzl(i,1,j) = 0.
           fqz(i,kde,j)  = 0.
          ! fqzl(i,kde,j) = 0.
         ENDDO

         DO k=kts+1,ktf
         DO i = i_start, i_end


		   
	   fqz(i,k,j) = sgmdot(i,k,j) *(fzm(k)*field(i,k,j)+fzp(k)*field(i,k-1,j))*(ztop-zs(i,j))/msft(i,j)    !这里不是m**2
		   

        ENDDO
        ENDDO

      ENDDO

    
   



! fqz end 
  
	  
	  	call average_vertical_pws( varw, var, rhofm_save,                 &
                     config_flags,                  &
                         fzm,fzp,dnw,                    &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )	



  !  计划拆分为简单公式，处理变量从半层插值到整层，及边界处理问题
  
  
        i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
	  
      DO  j = j_start, j_end
      DO  i = i_start, i_end
        rmsft(i,j)  = 1./msft(i,j) 
        rmm(i,j) = 1./(msft(i,j)*msft(i,j) )
      ENDDO
      ENDDO             

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end


           ! 20230817测试不用fqz，用原来方案。
             divg(i,k,j) = divg(i,k,j) -    ( varw(i,k+1,j)* rmm(i,j) * sgmdot(i,k+1,j) - varw(i,k,j)* rmm(i,j) * sgmdot(i,k,j)   ) /dsgm(k) * (ztop-zs(i,j))

            !      divg(i,k,j) = divg(i,k,j) -rdzw(k)*( fqz (i,k+1,j)-fqz (i,k,j) )
									  
			 !       divg(i,k,j) = divg(i,k,j) - ( fqz (i,k+1,j)-fqz (i,k,j) ) /dsgm(k)					  
									  
      ENDDO
      ENDDO
      ENDDO
	 
 

  
		 
        
   IF( config_flags%periodic_x   .or. &
       config_flags%symmetric_xs .or. &
       (its > ids+1)                ) degrade_xs = .false.
   IF( config_flags%periodic_x   .or. &
       config_flags%symmetric_xe .or. &
       (ite < ide-2)                ) degrade_xe = .false.
   IF( config_flags%periodic_y   .or. &
       config_flags%symmetric_ys .or. &
       (jts > jds+1)                ) degrade_ys = .false.
   IF( config_flags%periodic_y   .or. &
       config_flags%symmetric_ye .or. &
       (jte < jde-2)                ) degrade_ye = .false.

!--  y flux compute; these bounds are for periodic and sym b.c.

      ktf=MIN(kte,kde-1)
      i_start = its-1
      i_end   = MIN(ite,ide-1)+1
      j_start = jts-1
      j_end   = MIN(jte,jde-1)+1

!--  modify loop bounds if open or specified

      IF(degrade_xs) i_start = its
      IF(degrade_xe) i_end   = MIN(ite,ide-1)
      IF(degrade_ys) j_start = MAX(jts,jds+1)
      IF(degrade_ye) j_end = MIN(jte,jde-2)


      DO j = j_start, j_end+1
        DO k=kts,ktf
        DO i = i_start, i_end
 
		   
	!	   fqy(i,k,j) =  0.5*v(i,k,j)*          &
    !              ((ztop-zs(i,j))*field(i,k,j)/msft(i,j) + (ztop-zs(i,j-1))*field(i,k,j-1)/msft(i,j-1))	

				  
	 	    fqy(i,k,j) =    0.5*((ztop-zs(i,j-1))*var(i,k,j-1)*rmsft(i,j-1)  + (ztop-zs(i,j))*var(i,k,j)*rmsft(i,j) ) *v(i,k,j) 	  
		   
        ENDDO
        ENDDO
      ENDDO
	  
	  
	  ! 20230818 14:12 
	  j_start = jts 
      j_end   = MIN(jte,jde-1) 
	  
	  
      DO j = j_start, j_end
        DO k=kts,ktf
        DO i = i_start, i_end+1
 
				  
      !      fqx( i,k,j ) = 0.5*u(i,k,j)*          &
      !            ((ztop-zs(i,j))*field(i,k,j)/msft(i,j) + (ztop-zs(i-1,j))*field(i-1,k,j)/msft(i-1,j))			

		 	 fqx( i,k,j ) =  	  0.5*((ztop-zs(i-1,j))*var(i-1,k,j)/msft(i-1,j)  + (ztop-zs(i,j))*var(i,k,j)/msft(i,j) ) *u(i,k,j)   
 
        ENDDO
        ENDDO
      ENDDO
	  
	  


!!!!!!!!!!!!!!!!!!!!!!!	 
	 
	 
	 
!  ---- x direction -----------
 

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)


      IF ( .NOT. config_flags%periodic_x ) THEN
        IF ( config_flags%open_xs .or. specified ) i_start = MAX(ids+1,its)
        IF ( config_flags%open_xe .or. specified ) i_end   = MIN(ide-2,ite)
      ENDIF

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
  
    

       !  	divg(i,k,j) =  divg(i,k,j) -  (  0.5*((ztop-zs(i,j))*var(i,k,j)/msft(i,j)  + (ztop-zs(i+1,j))*var(i+1,k,j)/msft(i+1,j) ) *u(i+1,k,j) &
         !             -   0.5*((ztop-zs(i-1,j))*var(i-1,k,j)/msft(i-1,j)  + (ztop-zs(i,j))*var(i,k,j)/msft(i,j) ) *u(i,k,j) )*rdx
					 
     	divg(i,k,j) =  divg(i,k,j) -   rdx*( fqx (i+1,k,j)-fqx (i,k,j)   )					 
					 
					 
      ENDDO
      ENDDO
      ENDDO
	  
	!   write(6,*) '6 divg(1,1,1) = ', divg(1,1,1) 

!  ---- y direction -----------

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
      ! Polar boundary conditions are like open or specified
      IF ( config_flags%open_ys .or. specified .or. config_flags%polar ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified .or. config_flags%polar ) j_end   = MIN(jde-2,jte) 	  

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end

    !   	divg(i,k,j) = divg(i,k,j) -   ( 0.5*((ztop-zs(i,j))*var(i,k,j)*rmsft(i,j)  + (ztop-zs(i,j+1))*var(i,k,j+1)*rmsft(i,j+1) ) *v(i,k,j+1)    &
    !                                -   0.5*((ztop-zs(i,j-1))*var(i,k,j-1)*rmsft(i,j-1)  + (ztop-zs(i,j))*var(i,k,j)*rmsft(i,j) ) *v(i,k,j) )*rdy
					 
 		 
     	divg(i,k,j) = divg(i,k,j) -    rdy*( fqy (i,k,j+1)-fqy (i,k,j) )  					 

   

       ENDDO
      ENDDO
      ENDDO
	  
	  
     
	
		


  
        i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
         	divg(i,k,j) = divg(i,k,j)  * msft(i,j)*msft(i,j)/(ztop-zs(i,j)) 
      ENDDO
      ENDDO
      ENDDO
	  
 

      END SUBROUTINE calcu_divg_ps	  
	  
	  
	  
SUBROUTINE calcu_divg_pres( var, divg,                 &
                        u, v, sgmdot,                       &
                     config_flags,                  &
                        msft,zh,zs,fzm,fzp,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   ! Input data
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var,u,     &
                                                                      v, &
                                                                      sgmdot
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      !   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: mut
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: divg

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh,fzm,fzp


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   ! Local data
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx
         
        
       !  REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: zzh,dzzh
	   
	   
	   
	   
	   
		 
		    REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) ::  varw
		 
		 
		 
		 


      REAL , DIMENSION( kms:kme ) :: sgmh,dsgm
       REAL , DIMENSION( ims:ime , jms:jme )  ::  rmsft,rmm
	   
	   REAL :: ub,vb 
				 
				 
		!20230816		 
				 
	   REAL,  DIMENSION( its-1:ite+2, kts:kte, jts-1:jte+2  ) :: fqx, fqy, fqz
	   
	       REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) ::  field

   LOGICAL :: degrade_xs, degrade_ys
   LOGICAL :: degrade_xe, degrade_ye
   
   
		   

    LOGICAL :: specified




   degrade_xs = .true.
   degrade_xe = .true.
   degrade_ys = .true.
   degrade_ye = .true.


   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

   
     ktf=MIN(kte,kde-1)
 

	  DO k=kts,kde-1
          sgmh(k) = 0.5*(zh(k)+zh(k+1))
          dsgm(k) = zh(k+1) - zh(k) 
      ENDDO


 ! hanshenyou 20240128 
   fqx = 0.0 
   fqy = 0.0 
   fqz = 0.0 
	  

   ! 20230228  for rho couple scalar 
      CALL set_physical_bc3d( var, 'p', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              its, ite, jts, jte, kts, kte,       &
                              its, ite, jts, jte, kts, kte )
	  
	  

     divg = 0.0 
	 
	 
! boundary	 
	 
!  ---- x direction -----------
 

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)


      IF ( .NOT. config_flags%periodic_x ) THEN
        IF ( config_flags%open_xs .or. specified ) i_start = MAX(ids+1,its)
        IF ( config_flags%open_xe .or. specified ) i_end   = MIN(ide-2,ite)
      ENDIF

      
	  
	  
	IF ( config_flags%open_xs  .and.  its == ids ) THEN 
		
	!	i = ids  !  用于方法1 
		
      DO  j = j_start, j_end
      DO k=kts,ktf
   
 					 
	! 参考wrf	,这一方法的使用，改进非常大。避免了hill2km_x试验的西边界问题。			 
	         ub = MIN( 0.5*(u(its,k,j)+u(its+1,k,j)), 0. )
		 
         divg(its,k,j) =                      &
               - (1.0-zs(its,j)/ztop)/msft(its,j)*rdx*(                                             &
                       ub*(   var(its+1,k,j)                 &
                            - var(its  ,k,j)   ) +           &
                       var(its,k,j)*(u(its+1,k,j)-u(its,k,j))  &       ! modified 
                                                                )				 
					 
      ENDDO
      ENDDO
 
	  
	ENDIF	  
!	write(6,*) '1 divg(1,1,1) = ', divg(1,1,1)  

	  
	    IF ( config_flags%open_xe  .and.  ite == ide ) THEN 
		
	 !	i = ite-1
		
		 i_end   = MIN(ite,ide-1)
		 
		
      DO  j = j_start, j_end
      DO k=kts,ktf
 
  
     
	     ub = MAX( 0.5*(u(ite-1,k,j)+u(ite,k,j)), 0. )	 
		 
         divg(i_end,k,j) =                  &
               - (1.0-zs(i_end,j)/ztop)/msft(i_end,j)*rdx*(                                               &
                       ub*(  var(i_end  ,k,j)                  &
                           - var(i_end-1,k,j) ) +              &
                       var(i_end,k,j)*(u(ite,k,j)-u(ite-1,k,j))  &   ! modified  
                                                                    )				 
					 
					 
      ENDDO
      ENDDO
  

	ENDIF	
	  
	  
 
	  

!  ---- y direction -----------

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
      ! Polar boundary conditions are like open or specified
      IF ( config_flags%open_ys .or. specified .or. config_flags%polar ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified .or. config_flags%polar ) j_end   = MIN(jde-2,jte) 	  

       
	  
	    IF ( config_flags%open_ys  .and.  jts == jds ) THEN 
		
	!	j = jds 
		
      DO k=kts,ktf
      DO  i = i_start, i_end		

	
	      vb = MIN( 0.5*(v(i,k,jts)+v(i,k,jts+1)), 0. )
		  
        divg(i,k,jts) =    divg(i,k,jts)                &
               - (1.0-zs(i,jts)/ztop)/msft(i,jts)*rdy*(                                             &
                       vb*(  var(i,k,jts+1)                  &
                           - var(i,k,jts  ) ) +              &
                       var(i,k,jts)*(v(i,k,jts+1)-v(i,k,jts))  &     ! modified  
                                                               )
					 
      ENDDO
      ENDDO
	    ENDIF	
		
  !      write(6,*) '3 divg(1,1,1) = ', divg(1,1,1) 
	  
	    IF ( config_flags%open_ye  .and.  jte == jde ) THEN 
		
	 ! 	j = jte-1
		
		j_end   = MIN(jte,jde-1)
		
      DO k=kts,ktf
      DO  i = i_start, i_end		
		
	 
	
		  vb = MAX( 0.5*(v(i,k,jte-1)+v(i,k,jte)), 0. ) 
		 
         divg(i,k,j_end) =    divg(i,k,j_end)               &
               -  (1.0-zs(i,j_end)/ztop)/msft(i,j_end)*rdy*(                                               &
                       vb*(   var(i,k,j_end  )                 &
                            - var(i,k,j_end-1) ) +             &
                       var(i,k,j_end)*(v(i,k,jte)-v(i,k,jte-1))  &    ! modified   
                                                                    )
					 
      ENDDO
      ENDDO
	    ENDIF
		
 



!!!!!!!!!!!!!!!!!!!!!!!!!!!	 
	 
      field  = var  	 
	 
	 
!  ---- z  direction -----------	




! fqz

!-------------------- vertical advection




	  
!-- loop bounds for periodic or sym conditions

      i_start = its-1
      i_end   = MIN(ite,ide-1)+1
      j_start = jts-1
      j_end   = MIN(jte,jde-1)+1

!-- loop bounds for open or specified conditions

    IF(degrade_xs) i_start = MAX(its-1,ids)
    IF(degrade_xe) i_end   = MIN(ite+1,ide-1)
    IF(degrade_ys) j_start = MAX(jts-1,jds)
    IF(degrade_ye) j_end   = MIN(jte+1,jde-1)

    

      DO j = j_start, j_end
      DO i = i_start, i_end
		 
           fqz(i,1,j)  = 0.
        !   fqzl(i,1,j) = 0.
           fqz(i,kde,j)  = 0.
          ! fqzl(i,kde,j) = 0.
         ENDDO

         DO k=kts+1,ktf
         DO i = i_start, i_end


		   
	   fqz(i,k,j) = sgmdot(i,k,j) *(fzm(k)*field(i,k,j)+fzp(k)*field(i,k-1,j))*(1.0-zs(i,j)/ztop)/msft(i,j)    !这里不是m**2
		   

        ENDDO
        ENDDO

      ENDDO

    
   



! fqz end 
  
	  
	   	call average_vertical_w( varw, var,                 &
                      config_flags,                  &
                         fzm,fzp,                    &
                         ids, ide, jds, jde, kds, kde,  &
                         ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )	

	  !	call average_vertical_w_log( varw, var,                 &
      !               config_flags,                  &
      !                   fzm,fzp,                    &
       !                 ids, ide, jds, jde, kds, kde,  &
      !                  ims, ime, jms, jme, kms, kme,  &
       !                 its, ite, jts, jte, kts, kte  )	

  !  计划拆分为简单公式，处理变量从半层插值到整层，及边界处理问题
  
  
      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
	  
      DO  j = j_start, j_end
      DO  i = i_start, i_end
        rmsft(i,j)  = 1./msft(i,j)    ! 周期边界条件出错。
        rmm(i,j) = 1./(msft(i,j)*msft(i,j) )
      ENDDO
      ENDDO             

      DO  j = j_start, j_end
 
     DO k=kts,ktf
	 
      DO  i = i_start, i_end


           ! 20230817测试不用fqz，用原来方案。
             divg(i,k,j) = divg(i,k,j) -    ( varw(i,k+1,j)* rmm(i,j) * sgmdot(i,k+1,j) - varw(i,k,j)* rmm(i,j) * sgmdot(i,k,j)   ) /dsgm(k) * (1.0-zs(i,j)/ztop)

            !      divg(i,k,j) = divg(i,k,j) -rdzw(k)*( fqz (i,k+1,j)-fqz (i,k,j) )
									  
			 !       divg(i,k,j) = divg(i,k,j) - ( fqz (i,k+1,j)-fqz (i,k,j) ) /dsgm(k)		

  

			 
									  
      ENDDO
      ENDDO
      ENDDO
	  
 !write(6,*) 'max divg w    =  ', maxval(divg) , maxloc(divg) 
 !write(6,*) 'min divg w    =  ', minval(divg) , minloc(divg) 	 
  
		 
        
   IF( config_flags%periodic_x   .or. &
       config_flags%symmetric_xs .or. &
       (its > ids+1)                ) degrade_xs = .false.
   IF( config_flags%periodic_x   .or. &
       config_flags%symmetric_xe .or. &
       (ite < ide-2)                ) degrade_xe = .false.
   IF( config_flags%periodic_y   .or. &
       config_flags%symmetric_ys .or. &
       (jts > jds+1)                ) degrade_ys = .false.
   IF( config_flags%periodic_y   .or. &
       config_flags%symmetric_ye .or. &
       (jte < jde-2)                ) degrade_ye = .false.

!--  y flux compute; these bounds are for periodic and sym b.c.

      ktf=MIN(kte,kde-1)
      i_start = its-1
      i_end   = MIN(ite,ide-1)+1
      j_start = jts-1
      j_end   = MIN(jte,jde-1)+1

!--  modify loop bounds if open or specified

      IF(degrade_xs) i_start = its
      IF(degrade_xe) i_end   = MIN(ite,ide-1)
      IF(degrade_ys) j_start = MAX(jts,jds+1)
      IF(degrade_ye) j_end = MIN(jte,jde-1) ! hanshenyou 


      DO j = j_start, j_end 
        DO k=kts,ktf
        DO i = i_start, i_end
 
		   
	 

				  
	 	    fqy(i,k,j) =    0.5*((1.0-zs(i,j-1)/ztop)*var(i,k,j-1)*rmsft(i,j-1)  + (1.0-zs(i,j)/ztop)*var(i,k,j)*rmsft(i,j) ) *v(i,k,j) 	  
		   
        ENDDO
        ENDDO
      ENDDO
	  

 
  

  
	  ! 20230818 14:12 
	  j_start = jts 
      j_end   = MIN(jte,jde-1) 
	  
	  
      DO j = j_start, j_end
  
 
      DO k=kts,ktf
        DO i = i_start, i_end
 
	 	

		 	 fqx( i,k,j ) =  	  0.5*((1.0-zs(i-1,j)/ztop)*var(i-1,k,j)/msft(i-1,j)  + (1.0-zs(i,j)/ztop)*var(i,k,j)/msft(i,j) ) *u(i,k,j)   
 
  
 
        ENDDO
        ENDDO
      ENDDO
	  
 
 
	 
	 
	 
!  ---- x direction -----------
 

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)


      IF ( .NOT. config_flags%periodic_x ) THEN
        IF ( config_flags%open_xs .or. specified ) i_start = MAX(ids+1,its)
        IF ( config_flags%open_xe .or. specified ) i_end   = MIN(ide-2,ite)
      ENDIF

      DO  j = j_start, j_end
 
      DO k=kts,ktf
      DO  i = i_start, i_end
  
    

 
					 
     	divg(i,k,j) =  divg(i,k,j) -   rdx*( fqx (i+1,k,j)-fqx (i,k,j)   )		
 					 
					 
      ENDDO
      ENDDO
      ENDDO
	  
  !write(6,*) 'max divg w + x =  ', maxval(divg(ids:ide-1,kds:kde-1,jds:jde-1)) , maxloc(divg(ids:ide-1,kds:kde-1,jds:jde-1)) 
  !write(6,*) 'min divg w + x =  ', minval(divg(ids:ide-1,kds:kde-1,jds:jde-1)) , minloc(divg(ids:ide-1,kds:kde-1,jds:jde-1))  
 
 
	

!  ---- y direction -----------

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
      ! Polar boundary conditions are like open or specified
      IF ( config_flags%open_ys .or. specified .or. config_flags%polar ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified .or. config_flags%polar ) j_end   = MIN(jde-2,jte) 	  

      DO  j = j_start, j_end
 
      DO k=kts,ktf
      DO  i = i_start, i_end

 
 		 
     	divg(i,k,j) = divg(i,k,j) -    rdy*( fqy (i,k,j+1)-fqy (i,k,j) )  					 

    

       ENDDO
      ENDDO
      ENDDO
	  
	  
  !write(6,*) 'max divg w x y=  ', maxval(divg) , maxloc(divg) 
  !write(6,*) 'min divg w x y=  ', minval(divg) , minloc(divg)     

  
      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)

      DO  j = j_start, j_end
 
      DO k=kts,ktf
      DO  i = i_start, i_end
         	divg(i,k,j) = divg(i,k,j)  * msft(i,j)*msft(i,j)/(1.0-zs(i,j)/ztop) 
      ENDDO
      ENDDO
      ENDDO
	  
 
 
	  
 !write(6,*) 'max divgw x y/z=  ', maxval(divg) , maxloc(divg) 
 !write(6,*) 'min divgw x y/z=  ', minval(divg) , minloc(divg)  	 
 
 
 

      END SUBROUTINE calcu_divg_pres  
	  
! 20220710
 SUBROUTINE average_vertical_w_log( varw, var,                 &
                     config_flags,                  &
                         fzm,fzp,                    &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   ! Input data
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var

        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: varw
		 
		 	  REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) ::  fzm,fzp
			  
			  


   ! Local data
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx
         
     !    REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: varz,varzy,varzys




  !   REAL , DIMENSION( kms:kme ) :: sgmh,dsgm,dsgmh
    LOGICAL :: specified


	specified = .false.
	if(config_flags%specified .or. config_flags%nested) specified = .true.
   
     ktf=MIN(kte,kde-1)
 
 


      i_start = its
      i_end   = MIN(ide-1,ite)



      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
	  

      varw  = 0. 
	  
	  
! 4)  修订 开边界、给定边界、循环边界的x方向的 do循环 的范围 

    !  IF ( config_flags%open_xs ) i_start = MAX(ids+1,its)   ! 开边界的ids,ide, 需要特别处理
    !  IF ( config_flags%open_xe ) i_end   = MIN(ide-1,ite)

    !  IF ( specified ) i_start = MAX(ids+1,its)            ! 给定边界的，ids，ide，不做处理，ids+1, ide-1  二阶算法。
    !  IF ( specified ) i_end   = MIN(ide-2,ite)
   

    !  IF ( config_flags%periodic_x ) i_start = its           ! 周期边界，因为设置了外围，可以计算全部its至ite范围。
    !  IF ( config_flags%periodic_x ) i_end = ite


    !  IF ( config_flags%open_ys .or. specified ) j_start = MAX(jds+1,jts)
    !  IF ( config_flags%open_ye .or. specified ) j_end   = MIN(jde-2,jte)	  
	  
	  

      DO  j = j_start, j_end
      DO k=kts+1,ktf
      DO  i = i_start, i_end
       !   fzm(k)*var(i,k,j)  + fzp(k)*var(i,k-1,j)   
      	varw(i,k,j) = exp(  fzm(k)*log(var(i,k,j) ) + fzp(k)*log(var(i,k-1,j) ) )     
      ENDDO
      ENDDO
      ENDDO

! 上、下边界赋值
      DO      j = j_start, j_end
  
      DO  i = i_start, i_end
       varw(i,kts,j) =   2.0* var(i,kts,j)  -  varw(i,kts+1,j)  
	   
      varw(i,kte,j) =  2.0* var(i,kte-1,j)  -  varw(i,kte-1,j)   
      ENDDO
      ENDDO
	   
	  

      END SUBROUTINE average_vertical_w_log 	






! 20240303 
     SUBROUTINE calcu_engk_vw(engk,u,v,w,config_flags,                       &
                                  ids, ide, jds, jde, kds, kde,  &
                                  ims, ime, jms, jme, kms, kme,  &
                                  its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE

      INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
      INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
      INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte

      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( OUT )   ::  engk
      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( IN )    ::  u
      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( IN )    ::  v
      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( IN )    ::  w
	  TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      INTEGER      :: i,j,k,k_start,k_end,i_start,i_end,j_start,j_end

      REAL :: ub,vb 
	  
	  
     LOGICAL :: specified
	 
	 	 
	 integer :: ktf 
	 
	  ktf=MIN(kte,kde-1)
	  

    specified = .false.

    if(config_flags%specified .or. config_flags%nested) specified = .true.




      k_start = kts
      k_end   = MIN(kte,kde-1)
       
  !  for mpi 扩1      
    i_start = MAX(its-1,ids)
    i_end   = MIN(ite+1,ide-1)
    j_start = MAX(jts-1,jds)
    j_end   = MIN(jte+1,jde-1) 
	  

	   
	   
	   
!	  IF ( config_flags%open_xs ) i_start = MAX(ids+1,its-1)
!      IF ( config_flags%open_xe ) i_end   = MIN(ide-1,ite+1)
!      IF ( specified ) i_start = MAX(ids+2,its-1)
!      IF ( specified ) i_end   = MIN(ide-2,ite+1)
!      IF ( config_flags%periodic_x ) i_start = its-1
!      IF ( config_flags%periodic_x ) i_end = ite+1 
!	  
!      IF ( specified ) j_start = MAX(jds+2,jts-1)
!      IF ( specified ) j_end   = MIN(jte+1,jde-2)  
	  
       engk  = 0.0  
	  

      DO      j = j_start, j_end
      DO  k = k_start, k_end
      DO  i = i_start, i_end
         
			  engk(i,k,j) = ( v(i,k,j)**2 + v(i,k,j+1)**2 + w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25 
      ENDDO
      ENDDO
      ENDDO
	  
	  
if(.false.) then ! 20240706	
	  
      IF ( specified .AND. its .LE. ids+1 .AND. .NOT. config_flags%periodic_x ) THEN
	  j_start = MAX(jds+2,jts)
	  j_end   = MIN(jte,jde-2)
        DO j = j_start, j_end
        DO k=kts,ktf
           i = ids 
         
     
            engk(i,k,j) = (   v(i,k,j)**2 + v(i,k,j+1)**2 + w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25 
        ENDDO
        ENDDO
		
		
		
		
		
      ENDIF	  
	  
      IF ( specified .AND. ite .GE. ide-1 .AND. .NOT. config_flags%periodic_x ) THEN
	  j_start = MAX(jds+2,jts)
	  j_end   = MIN(jte,jde-2)
	  
        DO j = j_start, j_end
        DO k=kts,ktf
           i =  ide-1  
 
            engk(i,k,j) = (   v(i,k,j)**2 + v(i,k,j+1)**2 + w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25 
        ENDDO
        ENDDO
		
      ENDIF		  	  

	

     IF ( specified .AND. jts .LE. jds+1 ) THEN
     i_start = MAX(ids+2,its)
     i_end   = MIN(ide-2,ite)	 
	 
        j = jds 
        DO k=kts,ktf
        DO i = i_start, i_end
 
           vb = v(i,k,j)
           IF (v(i,k,j+1) .LT. 0.) vb = v(i,k,j+1)

               engk(i,k,j) = (  vb**2 + v(i,k,j+1)**2 + w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25 

        ENDDO
        ENDDO
      ENDIF

      IF ( specified .AND. jte .GE. jde-1 ) THEN
	  
	       i_start = MAX(ids+2,its)
           i_end   = MIN(ide-2,ite)	
	 
	 
        j = jde-1
        DO k=kts,ktf
        DO i = i_start, i_end

   
           vb = v(i,k,j+1)
           IF (v(i,k,j) .GT. 0.) vb = v(i,k,j)
	   
              engk(i,k,j) = ( v(i,k,j)**2 + vb**2 + w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25  

        ENDDO
        ENDDO
      ENDIF
	  
!  西南角 	
	
      IF ( specified .AND. its .LE. ids+1 .AND. .NOT. config_flags%periodic_x .AND. jts .LE. jds+1  ) THEN
	 	 
	    i = ids 
        j = jds 
		
        DO k=kts,ktf
 

   
           vb = v(i,k,j)
	
           IF (v(i,k,j+1) .LT. 0.) vb = v(i,k,j+1)   

              engk(i,k,j) = ( vb**2 + v(i,k,j+1)**2 + w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25  

 
        ENDDO
      ENDIF	  
	  
!  东南角 
    IF ( specified .AND. ite .GE. ide-1 .AND. .NOT. config_flags%periodic_x .AND. jts .LE. jds+1  ) THEN
	
	
       j = jds 
        DO k=kts,ktf
           i =  ide-1  
         
     
		     vb = v(i,k,j)
			 
   
		   IF (v(i,k,j+1) .LT. 0.) vb = v(i,k,j+1)    
		  
		  
            engk(i,k,j) = (    vb**2 + v(i,k,j+1)**2 + w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25 
        ENDDO
      
		
      ENDIF		  
	  
	  
	  
!  西北角 	
	
      IF ( specified .AND. its .LE. ids+1 .AND. .NOT. config_flags%periodic_x .AND. jte .GE. jde-1  ) THEN
	 	 
	    i = ids 
        j = jde-1 
		
        DO k=kts,ktf
 

        
		   
		   
		   vb = v(i,k,j+1)
           IF (v(i,k,j) .GT. 0.) vb = v(i,k,j)
		   
		   

              engk(i,k,j) = ( v(i,k,j)**2 + vb**2 +  w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25  

 
        ENDDO
      ENDIF	  	  


!  东北角 
    IF ( specified .AND. ite .GE. ide-1 .AND. .NOT. config_flags%periodic_x .AND. jte .GE. jde-1  ) THEN
	
	
       j = jde-1 
        DO k=kts,ktf
           i =  ide-1  
         
          
 


		   vb = v(i,k,j+1)
           IF (v(i,k,j) .GT. 0.) vb = v(i,k,j) 		   
		  
		  
            engk(i,k,j) = (   v(i,k,j+1)**2 + w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25 
        ENDDO
    
		
      ENDIF		  	 
	  
endif ! 20240706
	  
    END SUBROUTINE calcu_engk_vw
	
	

     SUBROUTINE calcu_engk_uw(engk,u,v,w,config_flags,                       &
                                  ids, ide, jds, jde, kds, kde,  &
                                  ims, ime, jms, jme, kms, kme,  &
                                  its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE

      INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
      INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
      INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte

      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( OUT )   ::  engk
      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( IN )    ::  u
      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( IN )    ::  v
      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( IN )    ::  w
	  TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      INTEGER      :: i,j,k,k_start,k_end,i_start,i_end,j_start,j_end
	  
	  REAL :: ub,vb 
	  
	  

     LOGICAL :: specified
	 
	 	 
	 integer :: ktf 
	 
	  ktf=MIN(kte,kde-1)
	  

    specified = .false.

    if(config_flags%specified .or. config_flags%nested) specified = .true.




      k_start = kts
      k_end   = MIN(kte,kde-1)
       
  !  for mpi 扩1      
    i_start = MAX(its-1,ids)
    i_end   = MIN(ite+1,ide-1)
    j_start = MAX(jts-1,jds)
    j_end   = MIN(jte+1,jde-1) 
	  

	   
	   
	   
!	  IF ( config_flags%open_xs ) i_start = MAX(ids+1,its-1)
!      IF ( config_flags%open_xe ) i_end   = MIN(ide-1,ite+1)
!      IF ( specified ) i_start = MAX(ids+2,its-1)
!      IF ( specified ) i_end   = MIN(ide-2,ite+1)
!      IF ( config_flags%periodic_x ) i_start = its-1
!      IF ( config_flags%periodic_x ) i_end = ite+1 
	  
!      IF ( specified ) j_start = MAX(jds+2,jts-1)
!      IF ( specified ) j_end   = MIN(jte+1,jde-2)  
	  
       engk  = 0.0  
	  

      DO      j = j_start, j_end
      DO  k = k_start, k_end
      DO  i = i_start, i_end
  
			  engk(i,k,j) = ( u(i,k,j)**2 + u(i+1,k,j)**2   + w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25 
      ENDDO
      ENDDO
      ENDDO
	  
	  
	  
if(.false.) then ! 20240706	
	  
      IF ( specified .AND. its .LE. ids+1 .AND. .NOT. config_flags%periodic_x ) THEN
	  j_start = MAX(jds+2,jts)
	  j_end   = MIN(jte,jde-2)
        DO j = j_start, j_end
        DO k=kts,ktf
           i = ids 
         
           ub = u(i ,k,j)
           IF (u(i+1,k,j) .LT. 0.) ub = u(i+1,k,j)
            engk(i,k,j) = ( ub**2 + u(i+1,k,j)**2  + w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25 
        ENDDO
        ENDDO
		
		
		
		
		
      ENDIF	  
	  
      IF ( specified .AND. ite .GE. ide-1 .AND. .NOT. config_flags%periodic_x ) THEN
	  j_start = MAX(jds+2,jts)
	  j_end   = MIN(jte,jde-2)
	  
        DO j = j_start, j_end
        DO k=kts,ktf
           i =  ide-1  
         
           ub = u(i+1 ,k,j)
          IF (u(i,k,j) .GT. 0.) ub = u(i,k,j) 
            engk(i,k,j) = ( ub**2 + u(i,k,j)**2  + w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25 
        ENDDO
        ENDDO
		
      ENDIF		  

     IF ( specified .AND. jts .LE. jds+1 ) THEN
     i_start = MAX(ids+2,its)
     i_end   = MIN(ide-2,ite)	 
	 
        j = jds 
        DO k=kts,ktf
        DO i = i_start, i_end
 
 

               engk(i,k,j) = ( u(i,k,j)**2 + u(i+1,k,j)**2 +  w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25 

        ENDDO
        ENDDO
      ENDIF

      IF ( specified .AND. jte .GE. jde-1 ) THEN
	  
	       i_start = MAX(ids+2,its)
           i_end   = MIN(ide-2,ite)	
	 
	 
        j = jde-1
        DO k=kts,ktf
        DO i = i_start, i_end

   
          
	   
              engk(i,k,j) = ( u(i,k,j)**2 + u(i+1,k,j)**2  + w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25  

        ENDDO
        ENDDO
      ENDIF
	  
!  西南角 	
	
      IF ( specified .AND. its .LE. ids+1 .AND. .NOT. config_flags%periodic_x .AND. jts .LE. jds+1  ) THEN
	 	 
	    i = ids 
        j = jds 
		
        DO k=kts,ktf
 

           ub = u(i ,k,j)
          
		   IF (u(i+1,k,j) .LT. 0.) ub = u(i+1,k,j)
         

              engk(i,k,j) = ( ub**2 + u(i+1,k,j)**2 +  w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25  

 
        ENDDO
      ENDIF	  
	  
!  东南角 
    IF ( specified .AND. ite .GE. ide-1 .AND. .NOT. config_flags%periodic_x .AND. jts .LE. jds+1  ) THEN
	
	
       j = jds 
        DO k=kts,ktf
           i =  ide-1  
         
           ub = u(i+1 ,k,j)
	 
			 
          IF (u(i,k,j) .GT. 0.) ub = u(i,k,j) 
  
		  
		  
            engk(i,k,j) = (  u(i,k,j)**2 + ub**2 +    w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25 
        ENDDO
     
		
      ENDIF		  
	  
	  
	  
!  西北角 	
	
      IF ( specified .AND. its .LE. ids+1 .AND. .NOT. config_flags%periodic_x .AND. jte .GE. jde-1  ) THEN
	 	 
	    i = ids 
        j = jde-1 
		
        DO k=kts,ktf
 
           ub = u(i ,k,j)
       
		   IF (u(i+1,k,j) .LT. 0.) ub = u(i+1,k,j)
     
              engk(i,k,j) = ( ub**2 + u(i+1,k,j)**2 +   w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25  
 
        ENDDO
      ENDIF	  	  


!  东北角 
    IF ( specified .AND. ite .GE. ide-1 .AND. .NOT. config_flags%periodic_x .AND. jte .GE. jde-1  ) THEN
	
	
       j = jde-1 
        DO k=kts,ktf
           i =  ide-1  
         
           ub = u(i+1 ,k,j)
 	 
          IF (u(i,k,j) .GT. 0.) ub = u(i,k,j) 
 
            engk(i,k,j) = (  u(i,k,j)**2 + ub**2 +    w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25 
        ENDDO
    
		
      ENDIF		  	 
	  

endif ! 20240706
	  
    END SUBROUTINE calcu_engk_uw



     SUBROUTINE calcu_engk_uv(engk,u,v,w,config_flags,                       &
                                  ids, ide, jds, jde, kds, kde,  &
                                  ims, ime, jms, jme, kms, kme,  &
                                  its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE

      INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
      INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
      INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte

      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( OUT )   ::  engk
      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( IN )    ::  u
      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( IN )    ::  v
      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( IN )    ::  w
	  TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      INTEGER      :: i,j,k,k_start,k_end,i_start,i_end,j_start,j_end


    REAL :: ub,vb 
	
     LOGICAL :: specified
	 
	 integer :: ktf 
	 
	  ktf=MIN(kte,kde-1)
	  
    specified = .false.

    if(config_flags%specified .or. config_flags%nested) specified = .true.




      k_start = kts
      k_end   = MIN(kte,kde-1)
       
  !  for mpi 扩1      
    i_start = MAX(its-1,ids)
    i_end   = MIN(ite+1,ide-1)
    j_start = MAX(jts-1,jds)
    j_end   = MIN(jte+1,jde-1) 
	  

	   
	   
	   
!	  IF ( config_flags%open_xs ) i_start = MAX(ids+1,its-1)
!      IF ( config_flags%open_xe ) i_end   = MIN(ide-1,ite+1)
!      IF ( specified ) i_start = MAX(ids+2,its-1)
!      IF ( specified ) i_end   = MIN(ide-2,ite+1)
!     IF ( config_flags%periodic_x ) i_start = its-1
!      IF ( config_flags%periodic_x ) i_end = ite+1 
	  
!      IF ( specified ) j_start = MAX(jds+2,jts-1)
!      IF ( specified ) j_end   = MIN(jte+1,jde-2)  
	  
       engk  = 0.0 

      DO      j = j_start, j_end
      DO  k = k_start, k_end
      DO  i = i_start, i_end
    
			  engk(i,k,j) = ( u(i,k,j)**2 + u(i+1,k,j)**2 + v(i,k,j)**2 + v(i,k,j+1)**2   ) * 0.25 
      ENDDO
      ENDDO
      ENDDO
	  
	  
	  
if(.false.) then ! 20240706	 
	  
      IF ( specified .AND. its .LE. ids+1 .AND. .NOT. config_flags%periodic_x ) THEN
	  j_start = MAX(jds+2,jts)
	  j_end   = MIN(jte,jde-2)
        DO j = j_start, j_end
        DO k=kts,ktf
           i = ids 
         
           ub = u(i ,k,j)
           IF (u(i+1,k,j) .LT. 0.) ub = u(i+1,k,j)
            engk(i,k,j) = ( ub**2 + u(i+1,k,j)**2 + v(i,k,j)**2 + v(i,k,j+1)**2  ) * 0.25 
        ENDDO
        ENDDO
		
		
		
		
		
      ENDIF	  
	  
      IF ( specified .AND. ite .GE. ide-1 .AND. .NOT. config_flags%periodic_x ) THEN
	  j_start = MAX(jds+2,jts)
	  j_end   = MIN(jte,jde-2)
	  
        DO j = j_start, j_end
        DO k=kts,ktf
           i =  ide-1  
         
           ub = u(i+1 ,k,j)
          IF (u(i,k,j) .GT. 0.) ub = u(i,k,j) 
            engk(i,k,j) = ( ub**2 + u(i,k,j)**2 + v(i,k,j)**2 + v(i,k,j+1)**2   ) * 0.25 
        ENDDO
        ENDDO
		
      ENDIF		  

     IF ( specified .AND. jts .LE. jds+1 ) THEN
     i_start = MAX(ids+2,its)
     i_end   = MIN(ide-2,ite)	 
	 
        j = jds 
        DO k=kts,ktf
        DO i = i_start, i_end
 
           vb = v(i,k,j)
           IF (v(i,k,j+1) .LT. 0.) vb = v(i,k,j+1)

               engk(i,k,j) = ( u(i,k,j)**2 + u(i+1,k,j)**2 + vb**2 + v(i,k,j+1)**2   ) * 0.25 

        ENDDO
        ENDDO
      ENDIF

      IF ( specified .AND. jte .GE. jde-1 ) THEN
	  
	       i_start = MAX(ids+2,its)
           i_end   = MIN(ide-2,ite)	
	 
	 
        j = jde-1
        DO k=kts,ktf
        DO i = i_start, i_end

   
           vb = v(i,k,j+1)
           IF (v(i,k,j) .GT. 0.) vb = v(i,k,j)
	   
              engk(i,k,j) = ( u(i,k,j)**2 + u(i+1,k,j)**2 + v(i,k,j)**2 + vb**2  ) * 0.25  

        ENDDO
        ENDDO
      ENDIF
	  
!  西南角 	
	
      IF ( specified .AND. its .LE. ids+1 .AND. .NOT. config_flags%periodic_x .AND. jts .LE. jds+1  ) THEN
	 	 
	    i = ids 
        j = jds 
		
        DO k=kts,ktf
 

           ub = u(i ,k,j)
           vb = v(i,k,j)
		   IF (u(i+1,k,j) .LT. 0.) ub = u(i+1,k,j)
           IF (v(i,k,j+1) .LT. 0.) vb = v(i,k,j+1)   

              engk(i,k,j) = ( ub**2 + u(i+1,k,j)**2 + vb**2 + v(i,k,j+1)**2   ) * 0.25  

 
        ENDDO
      ENDIF	  
	  
!  东南角 
    IF ( specified .AND. ite .GE. ide-1 .AND. .NOT. config_flags%periodic_x .AND. jts .LE. jds+1  ) THEN
	
	
       j = jds 
        DO k=kts,ktf
           i =  ide-1  
         
           ub = u(i+1 ,k,j)
		     vb = v(i,k,j)
			 
          IF (u(i,k,j) .GT. 0.) ub = u(i,k,j) 
		   IF (v(i,k,j+1) .LT. 0.) vb = v(i,k,j+1)    
		  
		  
            engk(i,k,j) = (  u(i,k,j)**2 + ub**2 +  vb**2 + v(i,k,j+1)**2  ) * 0.25 
        ENDDO
     
		
      ENDIF		  
	  
	  
	  
!  西北角 	
	
      IF ( specified .AND. its .LE. ids+1 .AND. .NOT. config_flags%periodic_x .AND. jte .GE. jde-1  ) THEN
	 	 
	    i = ids 
        j = jde-1 
		
        DO k=kts,ktf
 

           ub = u(i ,k,j)
       
		   IF (u(i+1,k,j) .LT. 0.) ub = u(i+1,k,j)
     
		   
		   
		   vb = v(i,k,j+1)
           IF (v(i,k,j) .GT. 0.) vb = v(i,k,j)
		   
		   

              engk(i,k,j) = ( ub**2 + u(i+1,k,j)**2 + v(i,k,j)**2 + vb**2  ) * 0.25  

 
        ENDDO
      ENDIF	  	  


!  东北角 
    IF ( specified .AND. ite .GE. ide-1 .AND. .NOT. config_flags%periodic_x .AND. jte .GE. jde-1  ) THEN
	
	
       j = jde-1 
        DO k=kts,ktf
           i =  ide-1  
         
           ub = u(i+1 ,k,j)
 	 
          IF (u(i,k,j) .GT. 0.) ub = u(i,k,j) 
 


		   vb = v(i,k,j+1)
           IF (v(i,k,j) .GT. 0.) vb = v(i,k,j) 		   
		  
		  
            engk(i,k,j) = (  u(i,k,j)**2 + ub**2 +  vb**2 + v(i,k,j+1)**2  ) * 0.25 
        ENDDO
      
		
      ENDIF		

	  
endif ! 20240706
	  
	  
    END SUBROUTINE calcu_engk_uv




     SUBROUTINE calcu_grady_new( var0, grady,                 &
 !                       u, v, w,                       &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   ! Input data
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var0

        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags


         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT( OUT) :: grady
		 
		 
         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   ! Local data
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdy
         
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: varz  !,varzy,varzys
		   REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: var 

       REAL , DIMENSION( kms:kme ) :: sgmh,dsgm,dsgmh
    LOGICAL :: specified

   ! 1) 判定 是否 给定边界
   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

 !2) 设定 do循环的 垂直 层次   
      ktf=MIN(kte,kde-1)


 

	  DO k=kts,kde-1
          sgmh(k) = 0.5*(zh(k)+zh(k+1))
          dsgm(k) = zh(k+1) - zh(k) 
      ENDDO

 	  DO k=2,kde-1
             dsgmh(k) = 0.5*(dsgm(k) + dsgm(k-1))   ! new 20220711
			 
      ENDDO



	 ! 3) 设定  x、 y do循环的 范围 
       
      i_start = its
      i_end   =  MIN(ite,ide-1)

      j_start = jts
      j_end   = MIN(jte,jde-1)
	  

       var = var0
	  
      CALL set_physical_bc3d( var, 'p', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              its, ite, jts, jte, kts, kte,       &
                              its, ite, jts, jte, kts, kte )
							  
							  

      IF ( config_flags%periodic_y ) j_start = jts            ! 20220708
      IF ( config_flags%periodic_y ) j_end = jte 


       varz = 0.0 

      DO  j = j_start, j_end
      DO k=2,ktf
      DO  i = i_start, i_end
    
      	varz(i,k,j) =  (var(i,k,j) - var(i,k-1,j) ) / dsgmh(k)
      ENDDO
      ENDDO
      ENDDO
	  
	  

! 上、下边界赋值
      DO      j = j_start, j_end
  
      DO  i = i_start, i_end
		varz(i,kts,j) =  varz(i,kts+1,j)  !varz(i,kts+3,j)  -3.0*varz(i,kts+2,j) + 3.0*varz(i,kts+1,j) 
	   
		varz(i,kte,j) =  varz(i,kte-1,j)  ! varz(i,kte-3,j)  -3.0*varz(i,kte-2,j) + 3.0*varz(i,kte-1,j)  
      ENDDO
      ENDDO
	  

 

	  

   !   DO j = j_start, j_end
   !   DO k=kts ,ktf
    !  DO  i = i_start, i_end
 
	!	varzys(i,k,j) =  ( varz(i,k,j-1) + varz(i,k,j) + varz(i,k+1,j-1) + varz(i,k+1,j)  ) *0.25    !   (varz(i+1,k,j) + varz(i,k,j))*0.5  + (varz(i+1,k+1,j) + varz(i,k+1,j))*0.5
    !  ENDDO
     ! ENDDO
     ! ENDDO
	 
	 ! hanshenyou 20220801
	 CALL set_physical_bc3d( varz, 'w', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              its, ite, jts, jte, kts, kte,       &
                              its, ite, jts, jte, kts, kte )
							  
							  

      CALL set_physical_bc2d( msft, 't', config_flags,  &
                              ids, ide, jds, jde,      &
                              ims, ime, jms, jme,      &
                              its, ite, jts, jte,      &
                              its, ite, jts, jte  )

      CALL set_physical_bc2d( zs, 't', config_flags,  &
                              ids, ide, jds, jde,      &
                              ims, ime, jms, jme,      &
                              its, ite, jts, jte,      &
                              its, ite, jts, jte  )	
	  



 ! 4)  修订 开边界、给定边界、循环边界的x方向的 do循环 的范围 
     ! y梯度与x无关。
    !  IF ( config_flags%open_xs ) i_start = MAX(ids+1,its)   ! 开边界的ids,ide, 需要特别处理
    !  IF ( config_flags%open_xe ) i_end   = MIN(ide-1,ite) 
	
      i_start = its
      i_end   =  MIN(ite,ide-1)
      j_start = jts
      j_end   = jte   
	  
      IF ( config_flags%periodic_y ) j_start = jts          ! 周期边界，因为设置了外围，
      IF ( config_flags%periodic_y ) j_end = jte



    IF ( config_flags%open_ys .or. specified ) j_start = MAX(jds+1,jts)
    IF ( config_flags%open_ye .or. specified ) j_end   = MIN(jde-1,jte)
	  

      grady = 0.0 
  
      DO  j = j_start, j_end
      DO  k = kts, ktf
      DO  i = i_start, i_end 
	  
	  
	     mrdy= 0.5* ( msft(i,j-1) + msft(i,j) ) *rdy
	   !  c1 = (msft(i,j)*(ztop-sgmh(k))/(ztop-zs(i,j)) + msft(i,j-1)*(ztop-sgmh(k))/(ztop-zs(i,j-1)))*0.5 
	   !  c2 = rdy*(zs(i,j)-zs(i,j-1))*( varz(i,k,j-1) + varz(i,k,j) + varz(i,k+1,j-1) + varz(i,k+1,j)  ) *0.25
	   !  grady(i,k,j) =       - mrdy*(var(i,k,j)-var(i,k,j-1)) - c1*c2
		 
	    
		 
		    c2  =   (   msft(i,j)  *(ztop-sgmh(k))/(ztop-zs(i,j))   )  *   (   varz(i,k,j)  + varz(i,k+1,j)  )               &  
		         +  (   msft(i,j-1)*(ztop-sgmh(k))/(ztop-zs(i,j-1 ))  )*   (   varz(i,k,j-1)  + varz(i,k+1,j-1)  )  
			   
			    c1 = rdy*(zs(i,j)-zs(i,j-1)) 
				
				
		 
	  !  grady(i,k,j) =       - mrdy*(var(i,k,j)-var(i,k,j-1)) - c1*c2*0.25
	  
	  ! 2024 
	    grady(i,k,j) =       - mrdy*(var(i,k,j)-var(i,k,j-1)) +  c1*c2*0.25
		


		 
	  ENDDO
      ENDDO
      ENDDO    
	  


  
	  

      END SUBROUTINE calcu_grady_new






	!   calcu_gradx 为求x方向梯度的子程序，变量var可以是动能engk或气压p，或扰动气压p'。 grady,gradz,为y和z方向。
     
      SUBROUTINE calcu_gradx_new( var, gradx,                 &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   ! Input data 

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var    
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      !   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: mut
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: gradx

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   ! Local data 
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx            
       !  REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: zzh,dzzh
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: varz 

        REAL , DIMENSION( kms:kme ) :: sgmh,dsgm,dsgmh
		
		 REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: vardx 
		
          REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  ::  varbc
		
        LOGICAL :: specified

! 1) 判定 是否 给定边界
   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

 !2) 设定 do循环的 垂直 层次   
      ktf=MIN(kte,kde-1)

 	  DO k=kts,kde-1
          sgmh(k) = 0.5*(zh(k)+zh(k+1))
          dsgm(k) = zh(k+1) - zh(k) 
      ENDDO

 	  DO k=2,kde-1
             dsgmh(k) = 0.5*(dsgm(k) + dsgm(k-1))   ! new 20220711
			 
      ENDDO
	  
	  

         varbc  = var
	  
      CALL set_physical_bc3d( varbc, 'p', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              its, ite, jts, jte, kts, kte,       &
                              its, ite, jts, jte, kts, kte )



! 3) 设定  x、 y do循环的 范围 	 
 
      ! x-point
      i_start = its
      i_end   = MIN(ite,ide-1)
	  
	  
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  

	  varz = 0.0 
	  
 DO  j = j_start, j_end	  
	 DO k=2,ktf
      DO  i = i_start, i_end
    
      	varz(i,k,j) =  (varbc(i,k,j) - varbc(i,k-1,j) ) / dsgmh(k)
      ENDDO
      ENDDO
 ! 上、下边界赋值,外插
      DO  i = i_start, i_end
       varz(i,kts,j) =  varz(i,kts+1,j)  !varz(i,kts+3,j)  -3.0*varz(i,kts+2,j) + 3.0*varz(i,kts+1,j) 
      varz(i,kte,j) =  varz(i,ktf,j)  ! varz(i,kte-3,j)  -3.0*varz(i,kte-2,j) + 3.0*varz(i,kte-1,j)  
      ENDDO
	  
	  
      ENDDO	  
  




      CALL set_physical_bc3d( varz, 'w', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              its, ite, jts, jte, kts, kte,       &
                              its, ite, jts, jte, kts, kte )
							  
							  

      CALL set_physical_bc2d( msft, 't', config_flags,  &
                              ids, ide, jds, jde,      &
                              ims, ime, jms, jme,      &
                              its, ite, jts, jte,      &
                              its, ite, jts, jte  )

      CALL set_physical_bc2d( zs, 't', config_flags,  &
                              ids, ide, jds, jde,      &
                              ims, ime, jms, jme,      &
                              its, ite, jts, jte,      &
                              its, ite, jts, jte  )	





	  
	gradx = 0.0 

	

      ! x-point
      i_start = its
      i_end   = ite 
	  
	  
      j_start = jts
      j_end   = MIN(jte,jde-1)	

	! 取消测边界外插，因为，给定边界条件不需要外插，周期边条件也不需要。开边界，倾向计算计划用迎风方案或最终倾向量直接外插等方法，也不需要中间量的外插。  
	  
	  
! 4)  修订 开边界、给定边界、循环边界的x方向的 do循环 的范围 

      IF ( config_flags%open_xs ) i_start = MAX(ids+1,its)   ! 开边界的ids,ide, 需要特别处理
      IF ( config_flags%open_xe ) i_end   = MIN(ide-1,ite)
      IF ( specified ) i_start = MAX(ids+1,its)            ! 给定边界的，ids，ide，不做处理，ids+1, ide-1  二阶算法。
      IF ( specified ) i_end   = MIN(ide-1,ite)
	  
      IF ( config_flags%periodic_x ) i_start = its           ! 周期边界，因为设置了外围，可以计算全部its至ite范围。
      IF ( config_flags%periodic_x ) i_end = ite


!！  下面两行多余，南北开边界，不影响x方向上的梯度计算。
 !     IF ( config_flags%open_ys .or. specified ) j_start = MAX(jds+1,jts)
 !     IF ( config_flags%open_ye .or. specified ) j_end   = MIN(jde-1,jte)
!WRF
!  now the east-west (x) pressure gradient

 !  itf=ite
 !  jtf=MIN(jte,jde-1)
 !  ktf=MIN(kte,kde-1)
 !  i_start = its
 !  j_start = jts
 !  IF ( (config_flags%open_xs .or. specified .or. &
!           config_flags%nested ) .and. its == ids ) i_start = its+1
 !  IF ( (config_flags%open_xe .or. specified .or. &
!           config_flags%nested ) .and. ite == ide ) itf = itf-1
 !  IF ( config_flags%periodic_x ) i_start = its
 !  IF ( config_flags%periodic_x ) itf=ite
!WRF	  

		
  loop_j :     DO  j = j_start, j_end
	
		 
	 
      DO  k = kts, ktf
      DO  i = i_start, i_end 
	  
	    mrdx= 0.5* ( msft(i-1,j) + msft(i,j) ) *rdx
		
	 !   c1 = (msft(i,j)*(ztop-sgmh(k))/(ztop-zs(i,j)) + msft(i-1,j)*(ztop-sgmh(k))/(ztop-zs(i-1,j)))*0.5
	 !   c2 = rdx*(zs(i,j)-zs(i-1,j))*( (varz(i-1,k,j) + varz(i,k,j)) + (varz(i-1,k+1,j) + varz(i,k+1,j)) ) *0.25	  
	


		    c2  =   (   msft(i,j)  *(ztop-sgmh(k))/(ztop-zs(i,j))   )  *   (   varz(i,k,j)  + varz(i,k+1,j)  )               &  
		         +  (   msft(i-1,j)*(ztop-sgmh(k))/(ztop-zs(i-1,j ))  )*   (   varz(i-1,k,j)  + varz(i-1,k+1,j)  )  
			   
			    c1 = rdx*(zs(i,j)-zs(i-1,j)) 
				

    !   gradx(i,k,j) =  - mrdx*(varbc(i,k,j)-varbc(i-1,k,j))  - c1*c2
	! 2024 
	   gradx(i,k,j) =  - mrdx*(varbc(i,k,j)-varbc(i-1,k,j))  + c1*c2
	  

	  ENDDO
      ENDDO	 

      ENDDO    loop_j
	  


      END SUBROUTINE calcu_gradx_new



! 改变流程，先边界，后。

SUBROUTINE calcu_divg( var, divg,                 &
                        u, v, sgmdot,                       &
                     config_flags,                  &
                        msft,zh,zs,fzm,fzp,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   ! Input data
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var,u,     &
                                                                      v, &
                                                                      sgmdot
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      !   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: mut
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: divg

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh,fzm,fzp


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   ! Local data
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx
         
        
       !  REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: zzh,dzzh
	   
	   
	   
	   
	   
		 
		    REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) ::  varw
		 
		 
		 
		 


      REAL , DIMENSION( kms:kme ) :: sgmh,dsgm
       REAL , DIMENSION( ims:ime , jms:jme )  ::  rmm    ! rmsft 
	   
	   REAL :: ub,vb 
				 
				 
		!20230816		 
				 
	   REAL,  DIMENSION( its-1:ite+2, kts:kte, jts-1:jte+2  ) :: fqx, fqy, fqz
	   
	       REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) ::  field

   LOGICAL :: degrade_xs, degrade_ys
   LOGICAL :: degrade_xe, degrade_ye
   
   
		   

    LOGICAL :: specified




   degrade_xs = .true.
   degrade_xe = .true.
   degrade_ys = .true.
   degrade_ye = .true.


   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

   
     ktf=MIN(kte,kde-1)
 

	  DO k=kts,kde-1
          sgmh(k) = 0.5*(zh(k)+zh(k+1))
          dsgm(k) = zh(k+1) - zh(k) 
      ENDDO


 



	  

   ! 20230228  for rho couple scalar 
      CALL set_physical_bc3d( var, 'p', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              its, ite, jts, jte, kts, kte,       &
                              its, ite, jts, jte, kts, kte )
	  
	  

     divg = 0.0 
	 
	 
! boundary	 
	 
!  ---- x direction -----------
 

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)


      IF ( .NOT. config_flags%periodic_x ) THEN
        IF ( config_flags%open_xs .or. specified ) i_start = MAX(ids+1,its)
        IF ( config_flags%open_xe .or. specified ) i_end   = MIN(ide-2,ite)
      ENDIF

      
	  
	  
	IF ( config_flags%open_xs  .and.  its == ids ) THEN 
		
	!	i = ids  !  用于方法1 
		
      DO  j = j_start, j_end
      DO k=kts,ktf
   
 					 
	! 参考wrf	,这一方法的使用，改进非常大。避免了hill2km_x试验的西边界问题。			 
	         ub = MIN( 0.5*(u(its,k,j)+u(its+1,k,j)), 0. )
		 
         divg(its,k,j) =                      &
               - (ztop-zs(its,j))/msft(its,j)*rdx*(                                             &
                       ub*(   var(its+1,k,j)                 &
                            - var(its  ,k,j)   ) +           &
                       var(its,k,j)*(u(its+1,k,j)-u(its,k,j))  &       ! modified 
                                                                )				 
					 
      ENDDO
      ENDDO
 
	  
	ENDIF	  
!	write(6,*) '1 divg(1,1,1) = ', divg(1,1,1)  

	  
	    IF ( config_flags%open_xe  .and.  ite == ide ) THEN 
		
	 !	i = ite-1
		
		 i_end   = MIN(ite,ide-1)
		 
		
      DO  j = j_start, j_end
      DO k=kts,ktf
 
  
     
	     ub = MAX( 0.5*(u(ite-1,k,j)+u(ite,k,j)), 0. )	 
		 
         divg(i_end,k,j) =                  &
               - (ztop-zs(i_end,j))/msft(i_end,j)*rdx*(                                               &
                       ub*(  var(i_end  ,k,j)                  &
                           - var(i_end-1,k,j) ) +              &
                       var(i_end,k,j)*(u(ite,k,j)-u(ite-1,k,j))  &   ! modified  
                                                                    )				 
					 
					 
      ENDDO
      ENDDO
  

	ENDIF	
	  
	  
	!  write(6,*) '2 divg(1,1,1) = ', divg(1,1,1)  
	  

!  ---- y direction -----------

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
      ! Polar boundary conditions are like open or specified
      IF ( config_flags%open_ys .or. specified .or. config_flags%polar ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified .or. config_flags%polar ) j_end   = MIN(jde-2,jte) 	  

       
	  
	    IF ( config_flags%open_ys  .and.  jts == jds ) THEN 
		
	!	j = jds 
		
      DO k=kts,ktf
      DO  i = i_start, i_end		

	
	      vb = MIN( 0.5*(v(i,k,jts)+v(i,k,jts+1)), 0. )
		  
        divg(i,k,jts) =    divg(i,k,jts)                &
               - (ztop-zs(i,jts))/msft(i,jts)*rdy*(                                             &
                       vb*(  var(i,k,jts+1)                  &
                           - var(i,k,jts  ) ) +              &
                       var(i,k,jts)*(v(i,k,jts+1)-v(i,k,jts))  &     ! modified  
                                                               )
					 
      ENDDO
      ENDDO
	    ENDIF	
		
  !      write(6,*) '3 divg(1,1,1) = ', divg(1,1,1) 
	  
	    IF ( config_flags%open_ye  .and.  jte == jde ) THEN 
		
	 ! 	j = jte-1
		
		j_end   = MIN(jte,jde-1)
		
      DO k=kts,ktf
      DO  i = i_start, i_end		
		
	 
	
		  vb = MAX( 0.5*(v(i,k,jte-1)+v(i,k,jte)), 0. ) 
		 
         divg(i,k,j_end) =    divg(i,k,j_end)               &
               -  (ztop-zs(i,j_end))/msft(i,j_end)*rdy*(                                               &
                       vb*(   var(i,k,j_end  )                 &
                            - var(i,k,j_end-1) ) +             &
                       var(i,k,j_end)*(v(i,k,jte)-v(i,k,jte-1))  &    ! modified   
                                                                    )
					 
      ENDDO
      ENDDO
	    ENDIF
		
		
     !     write(6,*) '4 divg(1,1,1) = ', divg(1,1,1) 







!!!!!!!!!!!!!!!!!!!!!!!!!!!	 
	 
      field  = var  	 
	 
	 
!  ---- z  direction -----------	




! fqz

!-------------------- vertical advection




	  
!-- loop bounds for periodic or sym conditions

      i_start = its 
      i_end   = MIN(ite,ide-1) 
      j_start = jts 
      j_end   = MIN(jte,jde-1) 

!-- loop bounds for open or specified conditions

    IF(degrade_xs) i_start = MAX(its-1,ids)
    IF(degrade_xe) i_end   = MIN(ite+1,ide-1)
    IF(degrade_ys) j_start = MAX(jts-1,jds)
    IF(degrade_ye) j_end   = MIN(jte+1,jde-1)

    

      DO j = j_start, j_end

         DO i = i_start, i_end
           fqz(i,1,j)  = 0.
        !   fqzl(i,1,j) = 0.
           fqz(i,kde,j)  = 0.
          ! fqzl(i,kde,j) = 0.
         ENDDO

         DO k=kts+1,ktf
         DO i = i_start, i_end


		   
	   fqz(i,k,j) = sgmdot(i,k,j) *(fzm(k)*field(i,k,j)+fzp(k)*field(i,k-1,j))*(ztop-zs(i,j))/msft(i,j)    !这里不是m**2
		   

        ENDDO
        ENDDO

      ENDDO

    
   



! fqz end 
  
	  
	  	call average_vertical_w( varw, var,                 &
                     config_flags,                  &
                         fzm,fzp,                    &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )	



  !  计划拆分为简单公式，处理变量从半层插值到整层，及边界处理问题
  
  
        i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
	  
      DO  j = j_start, j_end
      DO  i = i_start, i_end
      !  rmsft(i,j)  = 1./msft(i,j) 
        rmm(i,j) = 1./(msft(i,j)*msft(i,j) )
      ENDDO
      ENDDO             

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end


           ! 20230817测试不用fqz，用原来方案。
             divg(i,k,j) = divg(i,k,j) -    ( varw(i,k+1,j)* rmm(i,j) * sgmdot(i,k+1,j) - varw(i,k,j)* rmm(i,j) * sgmdot(i,k,j)   ) /dsgm(k) * (ztop-zs(i,j))

            !      divg(i,k,j) = divg(i,k,j) -rdzw(k)*( fqz (i,k+1,j)-fqz (i,k,j) )
									  
			 !       divg(i,k,j) = divg(i,k,j) - ( fqz (i,k+1,j)-fqz (i,k,j) ) /dsgm(k)		

 

			 
									  
      ENDDO
      ENDDO
      ENDDO
	 
 
		 
        
   IF( config_flags%periodic_x   .or. &
       config_flags%symmetric_xs .or. &
       (its > ids+1)                ) degrade_xs = .false.
   IF( config_flags%periodic_x   .or. &
       config_flags%symmetric_xe .or. &
       (ite < ide-2)                ) degrade_xe = .false.
   IF( config_flags%periodic_y   .or. &
       config_flags%symmetric_ys .or. &
       (jts > jds+1)                ) degrade_ys = .false.
   IF( config_flags%periodic_y   .or. &
       config_flags%symmetric_ye .or. &
       (jte < jde-2)                ) degrade_ye = .false.

!--  y flux compute; these bounds are for periodic and sym b.c.

      ktf=MIN(kte,kde-1)
	  
	  
	    ! 20240706
		
      i_start = its
      i_end   = MIN(ite,ide-1) 
      j_start = jts
      j_end   = MIN(jte+1,jde) 

!--  modify loop bounds if open or specified

  !    IF(degrade_xs) i_start = its
  !    IF(degrade_xe) i_end   = MIN(ite,ide-1)
  !    IF(degrade_ys) j_start = MAX(jts,jds+1)
  !    IF(degrade_ye) j_end = MIN(jte,jde-2)


      DO j = j_start, j_end 
        DO k=kts,ktf
        DO i = i_start, i_end
 
		   
	!	   fqy(i,k,j) =  0.5*v(i,k,j)*          &
    !              ((ztop-zs(i,j))*field(i,k,j)/msft(i,j) + (ztop-zs(i,j-1))*field(i,k,j-1)/msft(i,j-1))	

				  
	 ! 20240228	    fqy(i,k,j) =    0.5*((ztop-zs(i,j-1))*var(i,k,j-1)*rmsft(i,j-1)  + (ztop-zs(i,j))*var(i,k,j)*rmsft(i,j) ) *v(i,k,j) 	  
		   
		    fqy(i,k,j) =    0.5*((ztop-zs(i,j-1))*var(i,k,j-1)/msft(i,j-1)  + (ztop-zs(i,j))*var(i,k,j)/msft(i,j) ) *v(i,k,j) 	
        ENDDO
        ENDDO
      ENDDO
	  
	  
	  ! 20240706
      i_start = its
      i_end   = MIN(ite+1,ide) 
      j_start = jts
      j_end   = MIN(jte,jde-1) 
	  
	  
      DO j = j_start, j_end
        DO k=kts,ktf
        DO i = i_start, i_end 
 
				  
      !      fqx( i,k,j ) = 0.5*u(i,k,j)*          &
      !            ((ztop-zs(i,j))*field(i,k,j)/msft(i,j) + (ztop-zs(i-1,j))*field(i-1,k,j)/msft(i-1,j))			

		 	 fqx( i,k,j ) =  	  0.5*((ztop-zs(i-1,j))*var(i-1,k,j)/msft(i-1,j)  + (ztop-zs(i,j))*var(i,k,j)/msft(i,j) ) *u(i,k,j)   
 
 
 
 
        ENDDO
        ENDDO
      ENDDO
	  
	  


!!!!!!!!!!!!!!!!!!!!!!!	 
	 
	 
	 
!  ---- x direction -----------
 

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)


      IF ( .NOT. config_flags%periodic_x ) THEN
        IF ( config_flags%open_xs .or. specified ) i_start = MAX(ids+1,its)
        IF ( config_flags%open_xe .or. specified ) i_end   = MIN(ide-2,ite)
      ENDIF

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
  
    

       !  	divg(i,k,j) =  divg(i,k,j) -  (  0.5*((ztop-zs(i,j))*var(i,k,j)/msft(i,j)  + (ztop-zs(i+1,j))*var(i+1,k,j)/msft(i+1,j) ) *u(i+1,k,j) &
         !             -   0.5*((ztop-zs(i-1,j))*var(i-1,k,j)/msft(i-1,j)  + (ztop-zs(i,j))*var(i,k,j)/msft(i,j) ) *u(i,k,j) )*rdx
					 
     	divg(i,k,j) =  divg(i,k,j) -   rdx*( fqx (i+1,k,j)-fqx (i,k,j)   )		
 					 
					 
      ENDDO
      ENDDO
      ENDDO
	  
 !write(6,*) 'max divgw x =  ', maxval(divg) , maxloc(divg) 
  !write(6,*) 'min divgw x =  ', minval(divg) , minloc(divg) 

 
	 
	 
!  ---- y direction -----------

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
      ! Polar boundary conditions are like open or specified
      IF ( config_flags%open_ys .or. specified .or. config_flags%polar ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified .or. config_flags%polar ) j_end   = MIN(jde-2,jte) 	  

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end

    !   	divg(i,k,j) = divg(i,k,j) -   ( 0.5*((ztop-zs(i,j))*var(i,k,j)*rmsft(i,j)  + (ztop-zs(i,j+1))*var(i,k,j+1)*rmsft(i,j+1) ) *v(i,k,j+1)    &
    !                                -   0.5*((ztop-zs(i,j-1))*var(i,k,j-1)*rmsft(i,j-1)  + (ztop-zs(i,j))*var(i,k,j)*rmsft(i,j) ) *v(i,k,j) )*rdy
					 
 		 
     	divg(i,k,j) = divg(i,k,j) -    rdy*( fqy (i,k,j+1)-fqy (i,k,j) )  					 

    

       ENDDO
      ENDDO
      ENDDO
	  
	  
  !write(6,*) 'max divgw x y=  ', maxval(divg) , maxloc(divg) 
  !write(6,*) 'min divgw x y=  ', minval(divg) , minloc(divg)     
	
 		


  
        i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
         	divg(i,k,j) = divg(i,k,j)  * msft(i,j)*msft(i,j)/(ztop-zs(i,j)) 
      ENDDO
      ENDDO
      ENDDO
	  
  !write(6,*) 'max divgw x y/z=  ', maxval(divg) , maxloc(divg) 
  !write(6,*) 'min divgw x y/z=  ', minval(divg) , minloc(divg)  	 

      END SUBROUTINE calcu_divg




SUBROUTINE calcu_engmass(tmass,teng,tengk,tie,tpe,       &
                    config_flags,                  &
                    u,v,w,p,rho,                   &
                    zz,rdx,rdy,                    &
                    ids, ide, jds, jde, kds, kde,  &
                    ims, ime, jms, jme, kms, kme,  &
                    its, ite, jts, jte, kts, kte)

 IMPLICIT NONE


TYPE(grid_config_rec_type)    ,           INTENT(IN   ) :: config_flags                                                                        

  INTEGER,INTENT(IN) :: ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte

  REAL,INTENT(INOUT) :: tmass,teng,tengk,tie,tpe
  REAL,DIMENSION(ims:ime , kms:kme , jms:jme), INTENT(IN) :: u,v,w,p,rho,zz
  REAL,INTENT(IN) :: rdx ,rdy
  INTEGER :: i,k,j,i_end,j_end
  REAL,DIMENSION(ims:ime , kms:kme , jms:jme) :: rhoengk,internal_energy,potential_energy,vtj,z_at_mass,engk
  REAL    :: sum_vtj 


 call    calcu_engk(engk,u,v,w,config_flags,                       &                    
                                   ids, ide, jds, jde, kds, kde,  &                     
                                   ims, ime, jms, jme, kms, kme,  &                     
                                   its, ite, jts, jte, kts, kte  ) 
i_end = MIN(ite,ide-1)
j_end = MIN(jte,jde-1)

DO j=jts,j_end
  DO k=kts,kte-1 
    DO i=its,i_end
          z_at_mass(i,k,j) = (zz(i,k+1,j) + zz(i,k,j))/2
    ENDDO
  ENDDO
ENDDO


  vtj = 0.0
  
DO j=jts,j_end
  DO k=kts,kte-1 
    DO i=its,i_end
          vtj(i,k,j) = (1.0/rdx)*(1.0/rdy)* (zz(i,k+1,j) - zz(i,k,j))
!                     vtj(i,k,j) = (zz(i,k+1,j) - zz(i,k,j))  
   ENDDO
  ENDDO
ENDDO

  sum_vtj = 0.0

DO j=jts,j_end
  DO k=kts,kte-1 
    DO i=its,i_end
      sum_vtj = sum_vtj + vtj(i,k,j)
    ENDDO
  ENDDO
ENDDO

    write(6,*) 'sum_vtj = ', sum_vtj 
!  call cint(tmass,rho,vtj,sum_vtj,                 &
!                    ids, ide, jds, jde, kds, kde,  &
!                    ims, ime, jms, jme, kms, kme,  &
!                    its, ite, jts, jte, kts, kte)

  rhoengk = rho*engk
  potential_energy  = rho*z_at_mass*g
  internal_energy = p*cv/r_d

!  call cint(tengk,rhoengk,vtj,sum_vtj,             &
!                    ids, ide, jds, jde, kds, kde,  &
!                    ims, ime, jms, jme, kms, kme,  &
!                    its, ite, jts, jte, kts, kte)
!
!  call cint(tie,internal_energy,vtj,sum_vtj,       &
!                    ids, ide, jds, jde, kds, kde,  &
!                    ims, ime, jms, jme, kms, kme,  &
!                    its, ite, jts, jte, kts, kte)
!
!  call cint(tpe,potential_energy,vtj,sum_vtj,      &
!                    ids, ide, jds, jde, kds, kde,  &
!                    ims, ime, jms, jme, kms, kme,  &
!                    its, ite, jts, jte, kts, kte)
   tmass = 0.0
   tengk = 0.0
   tie  = 0.0
   tpe = 0.0


DO j=jts,j_end
  DO k=kts,kte-1 
    DO i=its,i_end
       tmass =  tmass + rho(i,k,j)*vtj(i,k,j)
       tengk =  tengk + rhoengk(i,k,j)*vtj(i,k,j)
       tie =    tie + internal_energy(i,k,j)*vtj(i,k,j)
       tpe =  tpe + potential_energy(i,k,j)*vtj(i,k,j)
    ENDDO
  ENDDO
ENDDO
   tmass = tmass/sum_vtj
   tengk = tengk/sum_vtj
   tie  = tie/sum_vtj
   tpe = tpe/sum_vtj 
  
  teng = tengk + tie + tpe
 
END SUBROUTINE calcu_engmass


!  计算动能engk

      SUBROUTINE calcu_engk_test(engk,rho,rhop,u,v,w,config_flags,                       &
                                  ids, ide, jds, jde, kds, kde,  &
                                  ims, ime, jms, jme, kms, kme,  &
                                  its, ite, jts, jte, kts, kte  )

! hanshenyou 2025 test ideal conservation ，delete bdy specified. 
      IMPLICIT NONE

      INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
      INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
      INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte

      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( OUT )   ::  engk
      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( IN )    ::  u
      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( IN )    ::  v
      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( IN )    ::  w
      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( IN )    :: rho,rhop
	  TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      INTEGER      :: i,j,k,k_start,k_end,i_start,i_end,j_start,j_end
	 
	  REAL :: ub,vb 

    LOGICAL :: specified
	
	 integer :: ktf 
	 
	  ktf=MIN(kte,kde-1)

    specified = .false.

    if(config_flags%specified .or. config_flags%nested) specified = .true.


    

      k_start = kts
      k_end   = MIN(kte,kde-1)
	  
	  
  !  for mpi 扩1      
    i_start = MAX(its-1,ids)
    i_end   = MIN(ite+1,ide-1)
    j_start = MAX(jts-1,jds)
    j_end   = MIN(jte+1,jde-1) 
	  
  
	   	   

	  
       engk  = 0.0 
	   
	   
      DO      j = j_start, j_end
      DO  k = k_start, k_end
      DO  i = i_start, i_end
        
			  engk(i,k,j) =  ( u(i,k,j)**2 + u(i+1,k,j)**2 + v(i,k,j)**2 + v(i,k,j+1)**2 + w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25 
              engk(i,k,j) = (rho(i,k,j)+rhop(i,k,j))*engk(i,k,j) 
      ENDDO
      ENDDO
      ENDDO
	  
	  
	
    END SUBROUTINE calcu_engk_test
	  
END MODULE module_advect_em

#endif
