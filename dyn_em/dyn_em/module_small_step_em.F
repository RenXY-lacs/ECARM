
!WRF:MODEL_LAYER:DYNAMICS
!

!  SMALL_STEP code for the geometric height coordinate model
!
!---------------------------------------------------------------------------

MODULE module_small_step_em

   USE module_configure
   USE module_model_constants
   USE module_advect_em, only: calcu_gradx,calcu_grady,calcu_gradz,calcu_divg,calcu_divg_pres,calcu_divg_p,  &
                                  calculate_sgmdot,calculate_sgmdot2,calcu_engk,     &
								   average_vertical_w,average_vertical_pw        !,testy
   

CONTAINS




SUBROUTINE small_step_prep_fm( u_1, u_2, v_1, v_2, w_1, w_2, &          ! hanshenyou
                            rho_1, rho_2, p_1, p_2,         &  ! hanshenyou
                            u_save, v_save, w_save,       &
                            rho_save, p_save,          &    ! hanshenyou
                            rk_step,                      &
                            ids,ide, jds,jde, kds,kde,    &
                            ims,ime, jms,jme, kms,kme,    &
                            its,ite, jts,jte, kts,kte    )


!20220704
  IMPLICIT NONE  ! religion first

! declarations for the stuff coming in

  INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
  INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
  INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte

  INTEGER,      INTENT(IN   )    :: rk_step

  REAL, DIMENSION(ims:ime, kms:kme, jms:jme),INTENT(INOUT) :: u_1,   &
                                                              v_1,   &
                                                              w_1,   &
                                                              rho_1,   &
                                                              p_1

  REAL, DIMENSION(ims:ime, kms:kme, jms:jme),INTENT(  OUT) :: u_save,   &
                                                              v_save,   &
                                                              w_save,   &
                                                              rho_save,   &
                                                              p_save

  REAL, DIMENSION(ims:ime, kms:kme, jms:jme),INTENT(INOUT) :: u_2,   &
                                                              v_2,   &
                                                              w_2,   &
                                                              rho_2,   &
                                                              p_2



! local variables

  INTEGER :: i, j, k
  INTEGER :: i_start, i_end, j_start, j_end, k_start, k_end
  INTEGER :: i_endu, j_endv


!<DESCRIPTION>
!
!  small_step_prep prepares the prognostic variables for the small timestep.
!  This includes switching time-levels in the arrays and computing coupled
!  perturbation variables for the small timestep
!  (i.e. mu*u" = mu(t)*u(t)-mu(*)*u(*); mu*u" is advanced during the small
!  timesteps
!
!</DESCRIPTION>

    i_start = its
    i_end   = min(ite,ide-1)
    j_start = jts
    j_end   = min(jte,jde-1)
    k_start = kts
    k_end = min(kte,kde-1)

    i_endu = ite
    j_endv = jte

    !  if this is the first RK step, reset *_1 to *_2
    !  (we are replacing the t-dt fields with the time t fields)

    IF ((rk_step == 1) ) THEN


      DO j=j_start, j_end
      DO k=k_start, k_end
      DO i=i_start, i_endu
        u_1(i,k,j) = u_2(i,k,j)
      ENDDO
      ENDDO
      ENDDO

      DO j=j_start, j_endv
      DO k=k_start, k_end
      DO i=i_start, i_end
        v_1(i,k,j) = v_2(i,k,j)
      ENDDO
      ENDDO
      ENDDO

      DO j=j_start, j_end
      DO k=k_start, min(kde,kte)
      DO i=i_start, i_end
        w_1(i,k,j)  = w_2(i,k,j)
      ENDDO
      ENDDO
      ENDDO
	  
	  
      DO j=j_start, j_end
      DO k=k_start, k_end
      DO i=i_start, i_end
        rho_1(i,k,j) = rho_2(i,k,j)
         p_1(i,k,j) = p_2(i,k,j)
      ENDDO
      ENDDO
      ENDDO	  
	  

    ELSE  ! (rk_step /= 1)




    END IF

    ! set up the small timestep variables

  !    DO j=j_start, j_end
  !    DO i=i_start, i_end
  !      ww_save(i,kde,j) = 0.
  !      ww_save(i,1,j) = 0.
  !    ENDDO
  !    ENDDO



    DO j=j_start, j_end
    DO k=k_start, k_end
    DO i=i_start, i_endu
      u_save(i,k,j) = u_2(i,k,j)
	  u_2(i,k,j) = u_1(i,k,j) -u_2(i,k,j)
    ENDDO
    ENDDO
    ENDDO

 

    DO j=j_start, j_endv
    DO k=k_start, k_end
    DO i=i_start, i_end
      v_save(i,k,j) = v_2(i,k,j)
	  v_2(i,k,j) = v_1(i,k,j) -v_2(i,k,j)

    ENDDO
    ENDDO
    ENDDO

    DO j=j_start, j_end
    DO k=k_start, k_end
    DO i=i_start, i_end
      rho_save(i,k,j) = rho_2(i,k,j)

    ENDDO
    ENDDO
    ENDDO
	
	DO j=j_start, j_end
    DO k=k_start, k_end
    DO i=i_start, i_end
	  p_save(i,k,j) = p_2(i,k,j)

    ENDDO
    ENDDO
    ENDDO
	

    DO j=j_start, j_end
    DO k=k_start, min(kde,kte)
 !   DO k=k_start, kde
    DO i=i_start, i_end
      w_save(i,k,j) = w_2(i,k,j)
	  w_2(i,k,j) = w_1(i,k,j) - w_2(i,k,j)
 
    ENDDO
    ENDDO
    ENDDO
	
	

	  
	  
	DO j=j_start, j_end
    DO k=k_start, k_end
    DO i=i_start, i_end
	  rho_2(i,k,j) = rho_1(i,k,j) -  rho_2(i,k,j) 
      p_2(i,k,j) = p_1(i,k,j) -  p_2(i,k,j)  
	  
	  
	  
    ENDDO
    ENDDO
    ENDDO

 
END SUBROUTINE small_step_prep_fm

!-------------------------------------------------------------------------


SUBROUTINE small_step_finish( u_2, u_1, v_2, v_1, w_2, w_1,    &
                              rho_2, rho_1, p_2, p_1, ww, ww1,   &
                              mu_2, mu_1,                      &
                              mut, muts, muu, muus, muv, muvs, &
                              c1h, c2h, c1f, c2f,              &
                              c3h, c4h, c3f, c4f,              &
                              u_save, v_save, w_save,          &
                              rho_save, p_save, mu_save,        &
                              msfux, msfuy, msfvx, msfvy,      &
                              msftx, msfty,                    &
                              h_diabatic,                      &
                              number_of_small_timesteps,dts,   &
                              rk_step, rk_order,               &
                              ids,ide, jds,jde, kds,kde,       &
                              ims,ime, jms,jme, kms,kme,       &
                              its,ite, jts,jte, kts,kte       )


  IMPLICIT NONE  ! religion first

!  stuff passed in

  INTEGER,                  INTENT(IN   ) :: ids,ide, jds,jde, kds,kde
  INTEGER,                  INTENT(IN   ) :: ims,ime, jms,jme, kms,kme
  INTEGER,                  INTENT(IN   ) :: its,ite, jts,jte, kts,kte
  INTEGER,                  INTENT(IN   ) :: number_of_small_timesteps
  INTEGER,                  INTENT(IN   ) :: rk_step, rk_order
  REAL,                     INTENT(IN   ) :: dts

  REAL,   DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN   ) :: u_1, &
                                                                 v_1, &
                                                                 w_1, &
                                                                 rho_1, &
                                                                 ww1, &
                                                                 p_1

  REAL,   DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT) :: u_2, &
                                                                 v_2, &
                                                                 w_2, &
                                                                 rho_2, &
                                                                 ww,  &
                                                                 p_2

  REAL, DIMENSION(ims:ime, kms:kme, jms:jme),INTENT(IN   ) :: u_save,   &
                                                              v_save,   &
                                                              w_save,   &
                                                              rho_save,   &
                                                              p_save,  &
                                                             h_diabatic

  REAL,   DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: muus, muvs
  REAL,   DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: mu_2, mu_1
  REAL,   DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: mut, muts, &
                                                        muu, muv, mu_save
  REAL,   DIMENSION(ims:ime, jms:jme), INTENT(IN   ) :: msfux, msfuy, &
                                                        msfvx, msfvy, &
                                                        msftx, msfty

  REAL, DIMENSION( kms:kme ),INTENT(IN   ) :: c1h, c2h, c1f, c2f, &
                                              c3h, c4h, c3f, c4f


! local stuff

  INTEGER         :: i,j,k
  INTEGER :: i_start, i_end, j_start, j_end, i_endu, j_endv

!<DESCRIPTION>
!
!  small_step_finish reconstructs the full uncoupled prognostic variables
!  from the coupled perturbation variables used in the small timesteps.
!
!</DESCRIPTION>

  i_start = its
  i_end   = min(ite,ide-1)
  j_start = jts
  j_end   = min(jte,jde-1)

  i_endu = ite
  j_endv = jte




!    addition of time level t back into variables

  DO j = j_start, j_endv
  DO k = kds, kde-1
  DO i = i_start, i_end
		v_2(i,k,j) = v_2(i,k,j) + v_save(i,k,j)  

  ENDDO
  ENDDO
  ENDDO

  DO j = j_start, j_end
  DO k = kds, kde-1
  DO i = i_start, i_endu
  	u_2(i,k,j) = u_2(i,k,j) + u_save(i,k,j) 

  ENDDO
  ENDDO
  ENDDO

     
	 
  DO j = j_start, j_end
  DO k = kds, kde
  DO i = i_start, i_end

  
   w_2(i,k,j) = w_2(i,k,j) + w_save(i,k,j) 
	
  ENDDO
  ENDDO
  ENDDO
  
  
  DO j = j_start, j_end
  DO k = kds, kde-1
  DO i = i_start, i_end
  
  		p_2(i,k,j) = p_2(i,k,j) + p_save(i,k,j)
		
  ENDDO
  ENDDO
  ENDDO
  
  
  DO j = j_start, j_end
  DO k = kds, kde-1
  DO i = i_start, i_end
  
		rho_2(i,k,j) = rho_2(i,k,j) + rho_save(i,k,j)

  ENDDO
  ENDDO
  ENDDO
  

END SUBROUTINE small_step_finish



SUBROUTINE advance_uv_fm   ( u, u_tend, v, v_tend,            &            
                        p,alt,cqu,cqv,                                  &
                        config_flags,spec_zone,               &
                        msft,zs,                              &
                        zh,fzm, fzp,                          &
                        rdx, rdy, dts, ztop,                  &
                        ids, ide, jds, jde, kds, kde,         &
                        ims, ime, jms, jme, kms, kme,         &
                        its, ite, jts, jte, kts, kte  )

   IMPLICIT NONE

   ! Input data
   
   TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags
   
   INTEGER,      INTENT(IN   )    :: spec_zone
   
   INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: u_tend, &
                                                                      v_tend, &
                                                                      p,   &
                                                                      alt,  &
																	  cqu,cqv																  
   
   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) :: u,v

   REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft, &
																	zs

   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: fzm,  &     !hanshenyou
                                                                  fzp,  &     ! hanshenyou  
                                                                  zh

   REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy,  &
																  dts,  &
																  ztop
																  
  
   
   ! Local data
   
 
   
    REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) :: tendency 
	
!	zh,zs,fzm,fzp,

 !  specified = .false.

      INTEGER :: i,j,k, i_start, i_end, j_start, j_end, k_start, k_end
      INTEGER :: i_endu, j_endv, k_endw
      INTEGER :: i_start_up, i_end_up, j_start_up, j_end_up
      INTEGER :: i_start_vp, i_end_vp, j_start_vp, j_end_vp
      INTEGER :: i_start_u_tend, i_end_u_tend, j_start_v_tend, j_end_v_tend
!<DESCRIPTION>
!
!  advance_uv advances the explicit perturbation horizontal momentum
!  equations (u,v) by adding in the large-timestep tendency along with
!  the small timestep pressure gradient tendency.
!
!</DESCRIPTION>
!  now, the real work.
!  set the loop bounds taking into account boundary conditions.
    IF( config_flags%nested .or. config_flags%specified ) THEN
      i_start = max( its,ids+spec_zone )
      i_end   = min( ite,ide-spec_zone-1 )
      j_start = max( jts,jds+spec_zone )
      j_end   = min( jte,jde-spec_zone-1 )
      k_start = kts
      k_end   = min( kte, kde-1 )
      i_endu = min( ite,ide-spec_zone )
      j_endv = min( jte,jde-spec_zone )
      k_endw = k_end
      IF( config_flags%periodic_x) THEN
        i_start = its
        i_end   = min(ite,ide-1)
        i_endu = ite
      ENDIF
    ELSE
      i_start = its
      i_end   = min(ite,ide-1)
      j_start = jts
      j_end   = min(jte,jde-1)
      k_start = kts
      k_end   = kte-1
      i_endu = ite
      j_endv = jte
      k_endw = k_end
    ENDIF
      i_start_up = i_start
      i_end_up   = i_endu
      j_start_up = j_start
      j_end_up   = j_end
      i_start_vp = i_start
      i_end_vp   = i_end
      j_start_vp = j_start
      j_end_vp   = j_endv
      IF ( (config_flags%open_xs   .or.     &
            config_flags%symmetric_xs   )   &
            .and. (its == ids) )            &
                 i_start_up = i_start_up + 1
      IF ( (config_flags%open_xe    .or.  &
            config_flags%symmetric_xe   ) &
             .and. (ite == ide) )         &
                 i_end_up   = i_end_up - 1
      IF ( (config_flags%open_ys    .or.   &
            config_flags%symmetric_ys  .or.   &
            config_flags%polar   )  &
                     .and. (jts == jds) )  &
                 j_start_vp = j_start_vp + 1
      IF ( (config_flags%open_ye     .or. &
            config_flags%symmetric_ye  .or.   &
            config_flags%polar   )  &
            .and. (jte == jde) )          &
                 j_end_vp   = j_end_vp - 1
      i_start_u_tend = i_start
      i_end_u_tend   = i_endu
      j_start_v_tend = j_start
      j_end_v_tend   = j_endv
      IF ( config_flags%symmetric_xs .and. (its == ids) ) &
                     i_start_u_tend = i_start_u_tend+1
      IF ( config_flags%symmetric_xe .and. (ite == ide) ) &
                     i_end_u_tend = i_end_u_tend-1
      IF ( config_flags%symmetric_ys .and. (jts == jds) ) &
                     j_start_v_tend = j_start_v_tend+1
      IF ( config_flags%symmetric_ye .and. (jte == jde) ) &
                     j_end_v_tend = j_end_v_tend-1



 	 
    tendency = 0.0 
	call    calcu_gradx( p, tendency,                 &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )


 

  u_outer_j_loop: DO j = j_start, j_end
   DO k = k_start, k_end
   DO i = i_start_u_tend, i_end_u_tend
     u(i,k,j) = u(i,k,j) + dts*u_tend(i,k,j)
	 
	 
   ENDDO
   ENDDO
    ENDDO u_outer_j_loop
	
	
	
	!write(6,*) j_start, j_end
!write(6,*) i_start_u_tend, i_end_u_tend
!write(6,*)  'u_tend(i_end_u_tend,k_start,j_end-3:j_end)'
    !write(6,*) u_tend(i_end_u_tend,k_start,j_end-3:j_end)
	
	
	
	
	
  u_outer_j_loop2: DO j = j_start, j_end
   DO k = k_start, k_end
     DO i = i_start_up, i_end_up
	 
	 
       u(i,k,j)=u(i,k,j) + dts*(alt(i-1,k,j)+alt(i,k,j))*cqu(i,k,j)*0.5*tendency(i,k,j)
	   
	   
	   
     ENDDO
   ENDDO
   ENDDO u_outer_j_loop2


		
    tendency = 0.0 
	
	!write(6,*) 'small p grady ' 
	
	call    calcu_grady( p, tendency,                 &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )   
   

 
	   

! now v
  v_outer_j_loop: DO j = j_start_v_tend, j_end_v_tend
  
  
   DO k = k_start, k_end
   DO i = i_start, i_end
	
   
     v(i,k,j) = v(i,k,j) + dts*v_tend(i,k,j)
   ENDDO
   ENDDO

	

     IF (     ( j >= j_start_vp)  &
       .and.( j <= j_end_vp  ) )  THEN
	   

     DO k = k_start, k_end
       DO i = i_start, i_end
	   
         v(i,k,j)=v(i,k,j)  + dts*(alt(i,k,j-1)+alt(i,k,j))*cqv(i,k,j)*0.5*tendency(i,k,j)
		 
		 
		 
       ENDDO
     ENDDO
   END IF
  ENDDO  v_outer_j_loop
  
 

END SUBROUTINE advance_uv_fm






!-----------------------------------------------------------------------
						
SUBROUTINE calc_p_rho_fm( al, p, pfm,  pb,                  &
                       alt, rho , rhob, pm1,      &
                       smdiv,             &
                       non_hydrostatic, step,        &
                       ids, ide, jds, jde, kds, kde, &
                       ims, ime, jms, jme, kms, kme, &
                       its,ite, jts,jte, kts,kte    )

  IMPLICIT NONE  ! religion first

! declarations for the stuff coming in

  INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
  INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
  INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte

  INTEGER,      INTENT(IN   )    :: step

  REAL, DIMENSION(ims:ime, kms:kme, jms:jme),INTENT(  OUT) :: al, alt,  &
                                                               p

  REAL, DIMENSION(ims:ime, kms:kme, jms:jme),INTENT(IN   ) ::  pfm, rho,  &
                                                               pb, rhob

  REAL, DIMENSION(ims:ime, kms:kme, jms:jme),INTENT(INOUT) ::  pm1


  REAL,                                       INTENT(IN   ) ::  smdiv


  LOGICAL, INTENT(IN   )  :: non_hydrostatic

! local variables

  INTEGER :: i, j, k
  INTEGER :: i_start, i_end, j_start, j_end, k_start, k_end
  REAL    :: ptmp

!<DESCRIPTION>
!
!  For the nonhydrostatic option,
!  calc_p_rho computes the perturbation inverse density and
!  perturbation pressure from the hydrostatic relation and the
!  linearized equation of state, respectively.
!
!  For the hydrostatic option,
!  calc_p_rho computes the perturbation pressure, perturbation density,
!  and perturbation geopotential
!  from the vertical coordinate definition, linearized equation of state
!  and the hydrostatic relation, respectively.
!
!  forward weighting of the pressure (divergence damping) is also
!  computed here.
!
!</DESCRIPTION>

   i_start = its
   i_end   = min(ite,ide-1)
   j_start = jts
   j_end   = min(jte,jde-1)
   k_start = kts
   k_end = min(kte,kde-1)


     DO j=j_start, j_end
     DO k=k_start, k_end
     DO i=i_start, i_end

!  al computation is all dry, so ok with moisture

      al(i,k,j) = rho(i,k,j) - rhob(i,k,j)

!  this is temporally linearized p, no moisture correction needed

      p(i,k,j) = pfm(i,k,j) - pb(i,k,j)
	  alt(i,k,j) = 1.0/rho(i,k,j)  
	  

     ENDDO
     ENDDO
     ENDDO


END SUBROUTINE calc_p_rho_fm



!----------------------------------------------------------------------

SUBROUTINE calc_p_rho_fm_new ( moist, n_moist,      &
                            al, alt,p,rhofm,rhofm_base,pfm,pfm_base,    &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            its, ite, jts, jte, kts, kte          )

  IMPLICIT NONE
   
   ! Input data



  INTEGER ,          INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                      ims, ime, jms, jme, kms, kme, &
                                      its, ite, jts, jte, kts, kte

  INTEGER ,          INTENT(IN   ) :: n_moist
 

  REAL, DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(IN   ) :: rhofm,rhofm_base,pfm,pfm_base 


  REAL, DIMENSION( ims:ime , kms:kme , jms:jme, n_moist ), INTENT(IN   ) :: moist

  REAL, DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(  OUT) :: al, alt,p


  ! Local stuff

  INTEGER :: i, j, k, kk, itf, jtf, ktf, ispe
  REAL    :: qtot, qf1, qf2, qvf
  
  real :: rhot
  
  REAL, DIMENSION( its:ite) :: temp,cpovcv_v
 ! REAL    :: pfu, phm, pfd

!<DESCRIPTION>
!
! For the nonhydrostatic option, calc_p_rho_phi calculates the
! diagnostic quantities pressure and (inverse) density from the
! prognostic variables using the equation of state.
!
! For the hydrostatic option, calc_p_rho_phi calculates the
! diagnostic quantities (inverse) density and geopotential from the
! prognostic variables using the equation of state and the hydrostatic
! equation.
!
!</DESCRIPTION>

  itf=MIN(ite,ide-1)
  jtf=MIN(jte,jde-1)
  ktf=MIN(kte,kde-1)

#ifndef INTELMKL
  cpovcv_v = cpovcv
#endif


          al = 0.0 
		 


  !   write(6,*) 'PARAM_FIRST_SCALAR,n_moist=',PARAM_FIRST_SCALAR,n_moist


      moist_nonhydro : IF (n_moist >= PARAM_FIRST_SCALAR ) THEN  

 
        DO j=jts,jtf

          DO k=kts,ktf  
            DO i=its,itf

            qtot = 0.
            DO ispe=PARAM_FIRST_SCALAR,n_moist
              qtot = qtot + moist(i,k,j,ispe)
            ENDDO
            qf2 = 1.
            qf1 = qtot*qf2
			
			! rhofm为干密度
			rhot   = rhofm(i,k,j)*(1.0+qtot) 
		 
			! alt为干密度倒数，水汽后面用cqu，cqv，cqw加。
            alt(i,k,j) = 1.0/rhofm(i,k,j)  
			
			
			! 	 rhofm_base  初始化可能为干密度，湿密度，对于的p基本态要考虑。
          	
            ! 考虑温度和湿度变化，基本态应该以干密度为合适，气压用干气压。	
   !	 扰动密度和扰动气压的初始态，考虑应该保持静力平衡。			
		!	al(i,k,j) = rhot  - rhofm_base(i,k,j)  
		
		!	29240422 测试初始场全密度作为基本态，这时rhofm_base应为初始全密度，计算形式不变
		     al(i,k,j) = rhot  - rhofm_base(i,k,j)  
			
			
	
     ! 	20240116 积分16小时13分30秒
	 !	 if( abs(al(i,k,j) ) > 0.5 ) then  !supercell bubble not zero 
	 !	 	write(6,*) i,k,j, rhofm(i,k,j)  ,rhofm_base(i,k,j)
	 !	 	write(6,*)   rhot  
	 !	 	write(6,*)  moist(i,k  ,j,2:4)
	 !	 	write(6,*)  pfm(i,k,j) - pfm_base(i,k,j) ,pfm(i,k,j) ,  pfm_base(i,k,j) 
	 !	 	!write(6,*) 10,3,10, rhofm(10,3,10)  ,rhofm_base(10,3,10), moist(10,3,10,2:4)
	 !	  	stop  'stop at calc_p_rho_fm_new' 
	 !	 endif
			
			

            p(i,k,j) = pfm(i,k,j) - pfm_base(i,k,j) 

            ENDDO
          ENDDO
        ENDDO
 

      ELSE moist_nonhydro      

        DO j=jts,jtf

          DO k=kts,ktf-1   
            DO i=its,itf

			al(i,k,j) = rhofm(i,k,j)  - rhofm_base(i,k,j) 
            alt(i,k,j) = 1.0/rhofm(i,k,j)  
      !      p(i,k,j) = pfm(i,k+1,j) - pfm_base(i,k,j)  
               p(i,k,j) = pfm(i,k,j) - pfm_base(i,k,j)   ! 20231215 
            ENDDO
          ENDDO

        ENDDO
	  




      END IF moist_nonhydro
	  


END SUBROUTINE calc_p_rho_fm_new 








SUBROUTINE calcu_abc_z( u,v,w, p,rho,w_tend,dp,drho,  &            
                        alt,cqw,cqx,rho_save,p_save, rhob,pfmb,                                 &
						a,b,c,wa,wb,ra,rb,pa,pb,pc,pd,pe,                &
                        config_flags,spec_zone,               &
                        msft,zs,                              &
                        zh,fzm, fzp,dnw,                          &
                        rdx, rdy, dts, ztop,                  &
                        ids, ide, jds, jde, kds, kde,         &
                        ims, ime, jms, jme, kms, kme,         &
                        its, ite, jts, jte, kts, kte  )

   IMPLICIT NONE
   
   
   ! rho,rho_save,rhob(rho0)

   ! Input data
   
   TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags
   
   INTEGER,      INTENT(IN   )    :: spec_zone
   
   INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte


  
   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN) :: w
    REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN  ) :: rho,p
	
	
   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: u,v

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: w_tend, drho,dp

    !20220630
   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) ::  alt,cqw,cqx,rho_save,p_save, rhob,pfmb
   
  
   REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft, &
																	zs

   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: fzm,  &     !hanshenyou
                                                                  fzp, dnw,  &     ! hanshenyou  
                                                                  zh

   REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy,  &
																  dts,  &
																  ztop
																  
  
     REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  ::  sigmdot

     REAL , DIMENSION( ims:ime , jms:jme ) :: rmm,rmsft,zs1 

     REAL , DIMENSION( kms:kme ) :: sgmh,dsgm,ax,bx 

 

   ! Local data
   
 
 !  specified = .false.

      INTEGER :: i,j,k, i_start, i_end, j_start, j_end, k_start, k_end
	  
	
	! implicit calculate w 
	
	 REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: dwf,drhof,dpf,dpf1
	 
	 
	  REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) :: w2
	  
	  
	  
      REAL, DIMENSION( ims:ime, kms:kme, jms:jme )  ::    dwsc,dpsc,drhosc,    &
                                                           ws,ps,rs
													       
      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(OUT) ::  a,b,c,wa,wb,ra,rb,pa,pb,pc,pd,pe 
														    
														
	   REAL, DIMENSION( ims:ime, kms:kme, jms:jme )  ::   pz 


       REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) ::  altw,pw,rhopw,rhox
       
       REAL ::  epssm 
       REAL ::  alfa,beta,gama
       REAL ::  c1,c2
	   REAL ::  cq1
	 real ::  tb,tu

       integer :: ktf



      ! gama = r_d/cv+1
       gama = rcv+1

       !!write(6,*)  'gamma  =  ' ,gama 

       ! epssm ????? 
       epssm = 1.0   ! 0.2 										
	 
       alfa = 0.5*(1.-epssm)
       beta = 0.5*(1.+epssm)
!<DESCRIPTION>
!
!  advance_w advances the implicit w and geopotential equations.
!
!</DESCRIPTION>
!  set loop limits.
!  Currently set for periodic boundary conditions
      i_start = its
      i_end   = min(ite,ide-1)
      j_start = jts
      j_end   = min(jte,jde-1)
	  
      k_start = kts
      k_end   = kte-1
	  
! hanshenyou 20231219 true
! 20230204
if (.true.) then 	  
      IF ( .NOT. config_flags%periodic_x )THEN
         IF ( config_flags%specified .or. config_flags%nested ) then
           i_start = max(its,ids+1)
           i_end   = min(ite,ide-2)
         ENDIF
      ENDIF
      IF ( config_flags%specified .or. config_flags%nested ) then
         j_start = max(jts,jds+1)
         j_end   = min(jte,jde-2)
      ENDIF

endif



     ktf = MIN(kte,kde-1)
   


 	  DO k=kts,ktf
          sgmh(k) = 0.5*(zh(k)+zh(k+1))
          dsgm(k) = zh(k+1) - zh(k) 
      ENDDO

 

      DO  j = j_start, j_end
      DO  i = i_start, i_end
        rmsft(i,j)  = 1./msft(i,j) 
        rmm(i,j) = 1./(msft(i,j)*msft(i,j) )
		zs1(i,j)  = zs(i,j)/ztop 
      ENDDO
      ENDDO


 
	 rhox = rho_save*(1.+cqx )  ! 20240705
	 
   if(.false.) then 
     !write(6,*) 'rho_save(ide/2,1:3,jde/2)i=', rho_save(ide/2,1:3,jde/2) 
     !write(6,*) 'cqx(ide/2,1:3,jde/2)i=', cqx(ide/2,1:3,jde/2)
     

      !write(6,*) 'rhox_(ide/2,1:3,jde/2)i=', rhox(ide/2,1:3,jde/2)
      !write(6,*) 'alt(ide/2,1:3,jde/2)i=', alt(ide/2,1:3,jde/2)
      !write(6,*) '1./rhox(ide/2,1:3,jde/2)i=', 1./rhox(ide/2,1:3,jde/2)
 
	 	call average_vertical_w( altw, 1./rhox,                 &
                     config_flags,                  &
                         fzm,fzp,                    &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )
						
						
    else
		 	call average_vertical_w( altw, alt,                 &
                     config_flags,                  &
                         fzm,fzp,                    &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )
						
						
       DO j = j_start, j_end
       DO k = kts, kte
       DO i = i_start, i_end
          cq1 = 1./(1.+cqw(i,k,j))
          altw(i,k,j) =    altw(i,k,j)*cq1 	   
       ENDDO
       ENDDO
       ENDDO					
	
	
	endif 
	
	
	
		

    pz = 0.0
	call    calcu_gradz( p, pz,                 &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )





  !   wa = dts*beta*dsgmzrw

    wa = 0.0 
	DO j = j_start, j_end
   DO k = kts, kte
     DO i = i_start, i_end
	 
 
        wa(i,k,j) =   dts*beta*altw(i,k,j)*(-1.)*(ztop/(ztop-zs(i,j)))

		
		   
     ENDDO
   ENDDO
   
   ENDDO

 
    DO  j = j_start, j_end
    DO  k = kts+1, kte-1
    DO  i = i_start, i_end     
      wa(i,k,j) = wa(i,k,j)*2.0/(dsgm(k)+dsgm(k-1))
    enddo
    enddo
    enddo




      wb = 0.0 
    !    j_start = MAX(jts,jds+1)
    !  j_end   = MIN(jte,jde-1)
      k_start = kts
      k_end   = MIN(kte,kde)
   !   i_start = MAX(its,ids+1)
   !   i_end   = MIN(ite,ide-1)

    DO  j = j_start, j_end
    DO  k = kts+1,kde-1  ! k_start, k_end
    DO  i = i_start, i_end
       wb(i,k,j) =   -1.*beta*dts*altw(i,k,j)*g
    enddo
    enddo
    enddo

 


!   初始化pc,pd
     pc = 0.0
	 pd = 0.0 
	 
   pc = 0.5*beta*dts*(1.0-gama)*rhob*g
   

   
   
   
  !   pd = 0.5*beta*dts*(gama-1.0)*rho_save
  
  ! cqw(1)= 0 ,cqw(kte) = 0
  
  
    DO  j = j_start, j_end
    DO  k = kts+1,kde-1  ! k_start, k_end
    DO  i = i_start, i_end
	
       pd(i,k,j)= 0.5*beta*dts*(gama-1.0)*rho_save(i,k,j)* ( 1.+ cqx(i,k,j)      ) 

	 enddo
    enddo
    enddo
   
  
	 
       pe = 0.0 



						
						
   !   IF ( config_flags%periodic_y ) j_start = jts
    !  IF ( config_flags%periodic_y ) j_end = jte 
	  
		
  
		
      DO  j = j_start, j_end
      DO  k = kts   , kte
      DO  i = i_start, i_end
            pe(i,k,j) =  -1.0 * altw(i,k,j)*pz(i,k,j)    ! pz have had  vertical metre  in gradz so *-1.0 
      ENDDO
      ENDDO 
      ENDDO


 
	   
!   specified = .false.

!! 1) 判定 是否 给定边界
  ! specified = .false.
  ! if(config_flags%specified .or. config_flags%nested) specified = .true.

 !2) 设定 do循环的 垂直 层次   
  



      k_end = MIN(kte,kde-1)

   
      

	  ra = 0.0 
	  pa = 0.0
	  
    DO  j = j_start, j_end
    DO  k = k_start, k_end
    DO  i = i_start, i_end

     ra(i,k,j) = dts*msft(i,j)*msft(i,j)* beta/dsgm(k)
  !   pa(i,k,j) = -1.0*ra(i,k,j)*gama
	  pa(i,k,j) = -1.0*gama*dts*msft(i,j)*msft(i,j)* beta/dsgm(k)
	enddo
    enddo
    enddo
	
	if(.false.)then 
	
	call average_vertical_pw( pb, p_save,pfmb,rhob,               &
                     config_flags,                  &
                         fzm,fzp,dnw,                    &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )
						
	endif 
	
	
	  call average_vertical_w( pb, p_save,                 &
                     config_flags,                  &
                         fzm,fzp,                    &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )	
						
						
						
						
	if(.false.)then 					
						
	call average_vertical_w( rb, rho_save,                 &
                     config_flags,                  &
                         fzm,fzp,                    &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )	
						
						
	endif 
	
	
  call average_vertical_w( rb, rhox,                 &
                     config_flags,                  &
                         fzm,fzp,                    &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )	
						
						
! cqw(1)= 0 ,cqw(kte) = 0
    DO  j = j_start, j_end
    DO  k = kts, kte
    DO  i = i_start, i_end
	 
	! rb(i,k,j) = (1.+cqw(i,k,j))*rb(i,k,j) *rmm(i,j)  
	  rb(i,k,j) =  rb(i,k,j) *rmm(i,j)  
	 
     pb(i,k,j) = (ztop/(ztop-zs(i,j)))*pb(i,k,j) *rmm(i,j)     
    enddo
    enddo
    enddo

	
	


      k_start = kts+1
      k_end   = MIN(kte,kde-1)

 


    DO  j = j_start, j_end
    DO  k = k_start, k_end
    DO  i = i_start, i_end     
    
       tb=dsgm(k)/(dsgm(k-1)+dsgm(k))
       tu=1.-tb

      a(i,k,j) = wa(i,k,j)*(pa(i,k-1,j)*pb(i,k-1,j)+pc(i,k-1,j)+pd(i,k-1,j)*pe(i,k-1,j)) - wb(i,k,j)*tb*(ra(i,k-1,j)*rb(i,k-1,j))

      c(i,k,j) = wa(i,k,j)*(pa(i,k,j)*pb(i,k+1,j)+pc(i,k,j)+pd(i,k,j)*pe(i,k+1,j))+wb(i,k,j)*tu*(ra(i,k,j)*rb(i,k+1,j))
      b(i,k,j) = wa(i,k,j)*(pc(i,k,j)+pd(i,k,j)*pe(i,k,j)-pa(i,k,j)*pb(i,k,j)-pa(i,k-1,j)*pb(i,k,j)-pc(i,k-1,j)-pd(i,k-1,j)*pe(i,k,j))    &
              + wb(i,k,j)*(tb*(ra(i,k-1,j)*rb(i,k,j))-tu*(ra(i,k,j)*rb(i,k,j))) -1.
			  
			  

		
    enddo
    enddo
    enddo


     

 !write(6,*) 'small a(ide/2,1:3,jde/2)', a(ide/2,1:3,jde/2) 
 !write(6,*) 'small b(ide/2,1:3,jde/2)', b(ide/2,1:3,jde/2)  
  !write(6,*) 'small c(ide/2,1:3,jde/2)', c(ide/2,1:3,jde/2)  

		  
 !write(6,*) 'small wa(ide/2,1:3,jde/2)', wa(ide/2,1:3,jde/2) 
 !write(6,*) 'small wb(ide/2,1:3,jde/2)', wb(ide/2,1:3,jde/2)  
   
 !write(6,*) 'small ra(ide/2,1:3,jde/2)', ra(ide/2,1:3,jde/2) 
 !write(6,*) 'small rb(ide/2,1:3,jde/2)', rb(ide/2,1:3,jde/2)  

  !write(6,*) 'small rho_save(ide/2,1:3,jde/2)', rho_save(ide/2,1:3,jde/2)  
  !write(6,*) 'small rmm(ide/2,jde/2)', rmm(ide/2,jde/2) 	
  
END SUBROUTINE calcu_abc_z


SUBROUTINE sumflux ( ru, rv, ww, rho,                            &
                     u_lin, v_lin, ww_lin,rho_lin,                   &
                   !  muu, muv,                               &
                   !  c1h, c2h, c1f, c2f,                     &
                   !  c3h, c4h, c3f, c4f,                     &
                     ru_m, rv_m, ww_m, rho_m,epssm,                &
                  !   msfux, msfuy, msfvx, msfvx_inv, msfvy,  &
                     iteration , number_of_small_timesteps,  &
                     ids,ide, jds,jde, kds,kde,              &
                     ims,ime, jms,jme, kms,kme,              &
                     its,ite, jts,jte, kts,kte              )


  IMPLICIT NONE  ! religion first

! declarations for the stuff coming in

  INTEGER,      INTENT(IN   )    :: number_of_small_timesteps
  INTEGER,      INTENT(IN   )    :: iteration
  INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
  INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
  INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte

  REAL, DIMENSION(ims:ime, kms:kme, jms:jme),  INTENT(IN   ) :: ru, &
                                                                rv, &
                                                                ww, &
                                                                u_lin,  &
                                                                v_lin,  rho,rho_lin, &
                                                                ww_lin


  REAL, DIMENSION(ims:ime, kms:kme, jms:jme) , INTENT(INOUT) :: ru_m, &
                                                                rv_m, rho_m, &
                                                                ww_m
!  REAL, DIMENSION(ims:ime, jms:jme) , INTENT(IN   ) :: muu, muv,      &
!                                                       msfux, msfuy,  &
!                                                       msfvx, msfvy, msfvx_inv

 ! REAL, DIMENSION( kms:kme ),INTENT(IN   ) :: c1h, c2h, c1f, c2f, &
 !                                             c3h, c4h, c3f, c4f

  INTEGER :: mini, minj, mink


  REAL, INTENT(IN   )  ::  epssm
  INTEGER   :: i,j,k


!<DESCRIPTION>
!
!  update the small-timestep time-averaged mass fluxes;  these
!  are needed for consistent mass-conserving scalar advection.
!
!</DESCRIPTION>

    IF (iteration == 1 )THEN
      DO  j = jts, jte
      DO  k = kts, kte
      DO  i = its, ite
        ru_m(i,k,j)  = 0.
        rv_m(i,k,j)  = 0.
        ww_m(i,k,j)  = 0.
		
      ENDDO
      ENDDO
      ENDDO
	  
	   rho_m = 0.0 
    ENDIF
	
      DO  j = jts, jte
      DO  k = kts, kte
      DO  i = its, ite
       rho_m(i,k,j)  = rho_m(i,k,j)  + rho(i,k,j) 
		
      ENDDO
      ENDDO
      ENDDO
	  
	! rho_m   = rho_m  + rho 

  mini = min(ide-1,ite)
  minj = min(jde-1,jte)
  mink = min(kde-1,kte)


    DO  j = jts, minj
    DO  k = kts, mink
    DO  i = its, mini
      ru_m(i,k,j)  = ru_m(i,k,j) + ru(i,k,j)
      rv_m(i,k,j)  = rv_m(i,k,j) + rv(i,k,j)
      ww_m(i,k,j)  = ww_m(i,k,j) + ww(i,k,j)
	  
    ENDDO
    ENDDO
    ENDDO
	
 

    IF (ite .GT. mini) THEN
      DO  j = jts, minj
      DO  k = kts, mink
      DO  i = mini+1, ite
        ru_m(i,k,j)  = ru_m(i,k,j) + ru(i,k,j)
      ENDDO
      ENDDO
      ENDDO
    END IF
    IF (jte .GT. minj) THEN
      DO  j = minj+1, jte
      DO  k = kts, mink
      DO  i = its, mini
 rv_m(i,k,j)  = rv_m(i,k,j) + rv(i,k,j)
      ENDDO
      ENDDO
      ENDDO
    END IF
    IF ( kte .GT. mink) THEN
      DO  j = jts, minj
      DO  k = mink+1, kte
      DO  i = its, mini
        ww_m(i,k,j)  = ww_m(i,k,j) + ww(i,k,j)
      ENDDO
      ENDDO
      ENDDO
    END IF

  IF (iteration == number_of_small_timesteps) THEN


      DO  j = jts, jte
      DO  k = kts, kte
      DO  i = its, ite
       rho_m(i,k,j)  = rho_m(i,k,j) / number_of_small_timesteps + rho_lin(i,k,j) 
		
      ENDDO
      ENDDO
      ENDDO



    DO  j = jts, minj
    DO  k = kts, mink
    DO  i = its, mini
      ru_m(i,k,j)  = ru_m(i,k,j) / number_of_small_timesteps   &
  !                   + (c1h(k)*muu(i,j)+c2h(k))*u_lin(i,k,j)/msfuy(i,j)
				      + u_lin(i,k,j)
					  
      rv_m(i,k,j)  = rv_m(i,k,j) / number_of_small_timesteps   &
                     + v_lin(i,k,j)
!                     + (c1h(k)*muv(i,j)+c2h(k))*v_lin(i,k,j)*msfvx_inv(i,j)					 
      ww_m(i,k,j)  = ww_m(i,k,j) / number_of_small_timesteps   &
                     + ww_lin(i,k,j)
    ENDDO
    ENDDO
    ENDDO


    IF (ite .GT. mini) THEN
      DO  j = jts, minj
      DO  k = kts, mink
      DO  i = mini+1, ite
        ru_m(i,k,j)  = ru_m(i,k,j) / number_of_small_timesteps   &
                     + u_lin(i,k,j) 
				!  + (c1h(k)*muu(i,j)+c2h(k))*u_lin(i,k,j)/msfuy(i,j)	 
					 
      ENDDO
      ENDDO
      ENDDO
    END IF
    IF (jte .GT. minj) THEN
      DO  j = minj+1, jte
      DO  k = kts, mink
      DO  i = its, mini
 rv_m(i,k,j)  = rv_m(i,k,j) / number_of_small_timesteps   &
                     + v_lin(i,k,j)
					 
	! + (c1h(k)*muv(i,j)+c2h(k))*v_lin(i,k,j)*msfvx_inv(i,j)				 
					 
      ENDDO
      ENDDO
      ENDDO
    END IF
    IF ( kte .GT. mink) THEN
      DO  j = jts, minj
      DO  k = mink+1, kte
      DO  i = its, mini
        ww_m(i,k,j)  = ww_m(i,k,j) / number_of_small_timesteps   &
                     + ww_lin(i,k,j)
      ENDDO
      ENDDO
      ENDDO
    END IF

  ENDIF


END SUBROUTINE sumflux



SUBROUTINE advance_w_p_rho_fm_abc( u,v,w, p,rho,w_tend,dp,drho,  &            
                        alt,cqw,cqx,rho_save,p_save, rhob,pfmb,                                 &
						a,b,c,wa,wb,ra,rb,pa,pb,pc,pd,pe,                         &
                        u_save,v_save,w_save,z,engk1,                                             &   !  2025
                        config_flags,spec_zone,               &
                        msft,zs,                              &
                        zh,fzm, fzp,dnw,                          &
                        rdx, rdy, dts, ztop, z_at_w,                 &
                        ids, ide, jds, jde, kds, kde,         &
                        ims, ime, jms, jme, kms, kme,         &
                        its, ite, jts, jte, kts, kte  )

   IMPLICIT NONE
   
   
   ! rho,rho_save,rhob(rho0)

   ! Input data
   
   TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags
   
   INTEGER,      INTENT(IN   )    :: spec_zone
   
   INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   )  :: a,b,c,wa,wb,ra,rb,pa,pb,pc,pd,pe
  
   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) :: w
    REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT   ) :: rho,p
	
	
   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: u,v

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: w_tend, drho,dp

    !20220630
   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) ::  alt,cqw,cqx,rho_save,p_save, rhob,pfmb
   
  
   REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft, &
																	zs

   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: fzm,  &     !hanshenyou
                                                                  fzp,dnw,  &     ! hanshenyou  
                                                                  zh

   REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy,  &
																  dts,  &
																  ztop
																  
																  
	 REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) ::  z_at_w															  
																  
  
     REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  ::  sigmdot

     REAL , DIMENSION( ims:ime , jms:jme ) :: rmm,rmsft,zs1 

     REAL , DIMENSION( kms:kme ) :: sgmh,dsgm,ax,bx 



 

   ! Local data
   
 
 !  specified = .false.

      INTEGER :: i,j,k, i_start, i_end, j_start, j_end, k_start, k_end
	  
	
	! implicit calculate w 
	
	 REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: dwf,drhof,dpf,dpf1
	 
	 
	  REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) :: w2
	  
	  
	  
      REAL, DIMENSION( ims:ime, kms:kme, jms:jme )  ::    dwsc,dpsc,drhosc,    &
                                                           ws,ps,rs
														
	   REAL, DIMENSION( ims:ime, kms:kme, jms:jme )  ::   pz,d 


       REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) ::  altw,pw,rhopw
       
       REAL ::  epssm 
       REAL ::  alfa,beta,gama
       REAL ::  c1,c2
	   REAL :: cq1
	   
	   
	 real ::  tb,tu

       integer :: ktf
	   
	   
	  ! 20230510  for w damping  
	  REAL, DIMENSION (kts:kte) :: dampwt
      real :: htop,hbot,hdepth,hk
      real    :: pi,dampmag
      REAL :: muthk, muthkm1

! 2025 hanshenyou 
      REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: u_save,v_save,w_save,z,engk1
       REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  ::  dengk,engk,engk2,letot
	   REAL :: co1,co2,co3





      ! gama = r_d/cv+1
       gama = rcv+1

       !!write(6,*)  'gamma  =  ' ,gama 

       ! epssm ????? 
       epssm = 1.0   ! 0.2 										
	 
       alfa = 0.5*(1.-epssm)
       beta = 0.5*(1.+epssm)
!<DESCRIPTION>
!
!  advance_w advances the implicit w and geopotential equations.
!
!</DESCRIPTION>
!  set loop limits.
!  Currently set for periodic boundary conditions
      i_start = its
      i_end   = min(ite,ide-1)
      j_start = jts
      j_end   = min(jte,jde-1)
	  
      k_start = kts
      k_end   = kte-1
	  
	  ! 20231219 这里有可能有问题 true 
	  ! 20230204
	  ! 20240114 true to false 
	  if (.true.) then   ! spec_bdy_final_fm 才赋值边界！
	  
      IF ( .NOT. config_flags%periodic_x )THEN
         IF ( config_flags%specified .or. config_flags%nested ) then
           i_start = max(its,ids+1)
           i_end   = min(ite,ide-2)
         ENDIF
      ENDIF
      IF ( config_flags%specified .or. config_flags%nested ) then
         j_start = max(jts,jds+1)
         j_end   = min(jte,jde-2)
      ENDIF

    endif ! 20230204 
	
	

     ktf = MIN(kte,kde-1)
   
   
   
   
  
   


 	  DO k=kts,ktf
          sgmh(k) = 0.5*(zh(k)+zh(k+1))
          dsgm(k) = zh(k+1) - zh(k) 
      ENDDO

 

      DO  j = j_start, j_end
      DO  i = i_start, i_end
        rmsft(i,j)  = 1./msft(i,j) 
        rmm(i,j) = 1./(msft(i,j)*msft(i,j) )
		zs1(i,j)  = zs(i,j)/ztop 
      ENDDO
      ENDDO




! 2025 hanshenyou 
	!	call calcu_engk( engk1,u+u_save,v+v_save,w+w_save,config_flags,     &
    !                     ids, ide, jds, jde, kds, kde,  &
    !                     ims, ime, jms, jme, kms, kme,  &
    !                     its, ite, jts, jte, kts, kte  )
     !   engk1 =(rho_save+rho)*engk1 



	  
	 	call average_vertical_w( altw, alt,                 &
                     config_flags,                  &
                         fzm,fzp,                    &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )	
						
	   DO j = j_start, j_end
       DO k = kts, kte
       DO i = i_start, i_end
          cq1 = 1./(1.+cqw(i,k,j))
          altw(i,k,j) =    altw(i,k,j)*cq1 	   
       ENDDO
       ENDDO
       ENDDO					
						
						
						

	 	call average_vertical_w( rhopw, rho,                 &
                     config_flags,                  &
                         fzm,fzp,                    &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )					

    pz = 0.0
	call    calcu_gradz( p, pz,                 &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )



  
 

	 dwf = 0.0 
	DO j = j_start, j_end
    DO k = kts, kte
    DO i = i_start, i_end
	
 
	   
	  dwf(i,k,j) = pz(i,k,j)*altw(i,k,j)    - rhopw(i,k,j)*altw(i,k,j) *g 	     	  
	   
    ENDDO
    ENDDO
    ENDDO
 
 
     w2=alfa*w 

     call    calculate_sgmdot( sigmdot,                 &
                        u, v, w2,                       &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  ) 

 			

    drhof = 0.0
    dpf   = 0.0 
     call calcu_divg( rho_save, drhof,                 &
                        u, v, sigmdot,                       &
                     config_flags,                  &
                        msft,zh,zs,fzm,fzp,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )
						
						
 
!20240131 测试山脉波稳定性
 ! call calcu_divg_pres      
    call calcu_divg(p_save,dpf,           & 
                       u, v, sigmdot,                       &
                   config_flags,                  &
                      msft,zh,zs,fzm,fzp,                    &
                       rdx,rdy,ztop,                  &
                       ids, ide, jds, jde, kds, kde,  &
                     ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )
						
 !    call calcu_divg_p( p_save, dpf,pfmb, rhob,                 &
 !                       u, v, sigmdot,                       &
 !                    config_flags,                  &
 !                       msft,zh,zs,fzm,fzp,dnw,                    &
  !                      rdx,rdy,ztop,                  &
  !                      ids, ide, jds, jde, kds, kde,  &
   !                     ims, ime, jms, jme, kms, kme,  &
   !                     its, ite, jts, jte, kts, kte  )



	dpf1 = 0.0 
	
 	  call     calculate_ptend_part2_fm_small_test( p, dpf1,                 &
!  call     calculate_ptend_part2_fm_small_tao2( p, dpf1,                 &
                        u, v, w2,sigmdot,                       &
						rho_save,rhob,         &
						alt,     &
                     config_flags,                  &
                        msft,zs,zh,fzm,fzp,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )



 

  ! dpf= gama*dpf + alfa*(1.-gama)*rho0*g*wc  + (gama-1)*rho_save*dpf1
   
   
	 dpf= cpovcv*dpf + dpf1					
						
 	   !write(6,*) ' maxval(abs(w_tend)) ,maxloc(abs(w_tend)):'
       !write(6,*) maxval(abs(w_tend)) ,maxloc(abs(w_tend)) 


     dwsc =  w_tend  + alfa*dwf 
	 
	 
	 
	      !write(6,*) 'maxval(abs(w)),maxval(abs(w_tend)) ,maxval(abs(dwf)):'
          !write(6,*) maxval(abs(w)),maxval(abs(w_tend)) ,maxval(abs(dwf))
		 

		 
 
	    !write(6,*) '==dp(z,1:4)= :', dp(ide/2,1:4,jde/2) 
	    !write(6,*) '==dpf(z,1:4)= :', dpf(ide/2,1:4,jde/2) 
		!write(6,*)  'max dp(pfm_tend) =  ', maxval(dp),maxloc(dp) 
     	!write(6,*)  'min dp(pfm_tend) =  ', minval(dp),minloc(dp) 
	   
     dpsc = dp + dpf
	 
     drhosc = drho + drhof


	  

     ws = dts*dwsc + w   
     ps = dts*dpsc + p
     rs = dts*drhosc+rho
	 
   ! write(6,*) 'maxval(abs(ws)),maxval(abs(rho)),maxval(abs(ps)):'
   !      write(6,*) maxval(abs(ws)),maxval(abs(rho)),maxval(abs(ps))

   ! write(6,*) 'maxval(abs(a)),maxval(abs(b)),maxval(abs(c)),maxval(abs(d)):'
   !      write(6,*) maxval(abs(a)),maxval(abs(b)),maxval(abs(c)),maxval(abs(d)) 
						
	!	stop

		
      IF ( config_flags%periodic_y ) j_start = jts
      IF ( config_flags%periodic_y ) j_end = jte 
	
       

!------------ for bdy w  -----------------------
	
!	     call    calculate_sgmdot( sigmdot,                 &
  call    calculate_sgmdot2( sigmdot,                 &     ! 20250714 
                     u, v, w2,                       &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  ) 
						
if(.false.) then !2025014 						
       CALL set_w_surface_fm_small( config_flags,              &
                           w, zs,  u, v,          &
                            rdx,  rdy,                  &
                            msft,                          &
                              ids, ide, jds, jde, kds, kde,                    &
                              ims, ime, jms, jme, kms, kme,                    &
                              its, ite, jts, jte, kts, kte                     )
    !  j_start = MAX(jts,jds+1)
    !  j_end   = MIN(jte,jde-2)
    !  i_start = MAX(its,ids+1)
    !  i_end   = MIN(ite,ide-2)
endif
	  
    DO  j = j_start, j_end
 
    DO  i = i_start, i_end
  
	   
    !    w(i,1,j) = 0.
	   w(i,1,j) = -sigmdot(i,1,j)/(ztop/(ztop-zs(i,j)))/beta  
      w(i,kte,j) = 0.
    ENDDO
    ENDDO


 

   !   j_start = MAX(jts,jds+1)
    !  j_end   = MIN(jte,jde-2)
      k_start = kts+1
      k_end   = MIN(kte,kde-1)
    !  i_start = MAX(its,ids+1)
    !  i_end   = MIN(ite,ide-2)
 

    
    DO  j = j_start, j_end
    DO  k = k_start, k_end
    DO  i = i_start, i_end     
    
       tb=dsgm(k)/(dsgm(k-1)+dsgm(k))
       tu=1.-tb


			
			  

        d(i,k,j) = -1.*(wa(i,k,j)*(ps(i,k,j)-ps(i,k-1,j))+wb(i,k,j)*(tb*rs(i,k-1,j)+tu*rs(i,k,j)) + ws(i,k,j))
		
  	
		
		
    enddo
    enddo
    enddo
	
 


     k=2
    DO  j = j_start, j_end
    DO  i = i_start, i_end     
      d(i,k,j) = d(i,k,j)-a(i,k,j)*w(i,k-1,j)
    enddo
    enddo
	
 
     
     k=kte-1
    DO  j = j_start, j_end
    DO  i = i_start, i_end     
      d(i,k,j) = d(i,k,j)-c(i,k,j)*w(i,k+1,j)
    enddo
    enddo

 
		  
		  
    DO  j = j_start, j_end
    DO  i = i_start, i_end  
	
	
     k = 2
     ax(k)=c(i,k,j)/b(i,k,j)
     bx(k)=d(i,k,j)/b(i,k,j)
	 
	 
	

     do k = 3,kte-2
       ax(k) = c(i,k,j)/(b(i,k,j) - a(i,k,j)*ax(k-1))
	   
	
	   
	   
     enddo

     do k = 3,kte-2
       bx(k) = (d(i,k,j)-a(i,k,j)*bx(k-1))/(b(i,k,j) - a(i,k,j)*ax(k-1))
	
    
	   
     enddo

     k=kte-1
     w(i,k,j) = (d(i,k,j) - a(i,k,j)*bx(k-1))/(b(i,k,j) - a(i,k,j)*ax(k-1))
	 
 
	
	   
     do k = kte-2,2,-1
	 
	 
       w(i,k,j) = bx(k) - ax(k)*w(i,k+1,j)
	
 
     enddo
 

    enddo
    enddo
	
	
	
	
	pi = 4.*atan(1.)
   dampmag = dts*config_flags%dampcoef
   hdepth=config_flags%zdamp
   
	
	
	
	

	
	
   IF (config_flags%damp_opt .eq. 3) THEN
   
      DO  j = j_start, j_end	
   
      DO k=k_end+1,2,-1
      DO i=i_start, i_end
          htop=ztop
          hk=z_at_w(i,k,j)
          hbot=htop-hdepth
          dampwt(k) = 0.
          if(hk .ge. hbot)then
            dampwt(k) = dampmag*sin(0.5*pi*(hk-hbot)/hdepth)*sin(0.5*pi*(hk-hbot)/hdepth)
          endif
          w(i,k,j) = (w(i,k,j) - dampwt(k)*w(i,k,j))/(1.+dampwt(k))
      ENDDO
      ENDDO


    

  ENDDO 
	
    ENDIF	
	
	
	
	
	
	
	


 
      k_start = kts 


   if(.false.) then    ! 密度计算的守恒算法选择
   
   ! overfloat 
    do  j = j_start, j_end
    do  k = k_start, k_end
    do  i = i_start, i_end     
         rho(i,k,j) = rs(i,k,j)  + ra(i,k,j)*(rb(i,k+1,j)*w(i,k+1,j)-rb(i,k,j)*w(i,k,j))
    enddo
    enddo
    enddo
	
	
	else
	
		     call    calculate_sgmdot( sigmdot,                 &
                        u, v, w,                       &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  ) 
						
 
    drhof = 0.0  ! add 20250711

		call calcu_divg( rho_save, drhof,                 &
                        u, v, sigmdot,                       &
                     config_flags,                  &
                        msft,zh,zs,fzm,fzp,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )
						
						
	
	  drhosc= drho+drhof 

   !  rho = dts*drhosc+rho 
	  
    do  j = j_start, j_end
    do  k = k_start, k_end
    do  i = i_start, i_end     
         rho(i,k,j) = dts*drhosc(i,k,j)+rho(i,k,j) 
    enddo
    enddo
	enddo
	
	endif
	

 if(.false.)  then
	 
    do  j = j_start, j_end
    do  k = k_start, k_end
    do  i = i_start, i_end     
        p(i,k,j) = ps(i,k,j)    &
 	         + pa(i,k,j)*(pb(i,k+1,j)*w(i,k+1,j)-pb(i,k,j)*w(i,k,j))  &  
                       +pc(i,k,j)*(w(i,k,j)+w(i,k+1,j))          &
	 				  +pd(i,k,j)*(pe(i,k,j)*w(i,k,j)+pe(i,k+1,j)*w(i,k+1,j))

  


    enddo
    enddo
    enddo

 else  

		call calcu_engk( engk, u_save,v_save,w_save,config_flags,     &
                         ids, ide, jds, jde, kds, kde,  &
                         ims, ime, jms, jme, kms, kme,  &
                         its, ite, jts, jte, kts, kte  )
		call calcu_engk( engk2,u+u_save,v+v_save,w+w_save,config_flags,     &
                         ids, ide, jds, jde, kds, kde,  &
                         ims, ime, jms, jme, kms, kme,  &
                         its, ite, jts, jte, kts, kte  )

    	 engk2 =(rho_save+rho)*engk2 ! rho is new 
     	 dengk = (engk2-engk1)/dts 

    	call calculate_sgmdot( sigmdot,                     &
                              u+u_save, v+v_save, w+w_save, &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                              ids, ide, jds, jde, kds, kde,  &
                              ims, ime, jms, jme, kms, kme,  &
                              its, ite, jts, jte, kts, kte  )

		co1 = cv/r_d
		co2 = 1./co1
		co3 = 1.+co1
        CALL calcu_divg( rho_save*engk+co3*p_save+rho_save*g*z,   &
						 letot,                 &
                         u+u_save, v+v_save, sigmdot,                       &
                      config_flags,                  &
                        msft,zh,zs,fzm,fzp,                    &
                        rdx,rdy,ztop,                  &
                         ids, ide, jds, jde, kds, kde,  &
                         ims, ime, jms, jme, kms, kme,  &
                         its, ite, jts, jte, kts, kte  )

 

      p = dts*(-dengk-drhosc*g*z+letot)*co2+p  
endif

	 
	 



	  

END SUBROUTINE advance_w_p_rho_fm_abc



!!!!!

SUBROUTINE   calculate_ptend_part2_fm_small_test( pp, tendency,                 &
                        u, v, w,sgmdot,                       &
						rho,rhob,         &
						alt,             &
                     config_flags,                  &
                        msft,zs,zh,fzm,fzp,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   ! Input data
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte
											  
											  
		REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: pp

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: u,     &
                                                                      v, w,&
                                                                      sgmdot
																	  
																	  
	    REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: rho,rhob,alt  ! ,pb
		

		 
		 
																		 
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      !   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: mut
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT( INOUT) :: tendency

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,zs


         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh,fzm,fzp

         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   ! Local data
   
   !		 REAL , DIMENSION( ims:ime , kms:kme , jms:jme )   ::  alt 
   
         INTEGER :: i, j, k, ktf
         INTEGER :: i_start, i_end, j_start, j_end
      !   INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
      !   INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
       !  INTEGER :: jp1, jp0, jtmp
       !  REAL :: c1,c2,mrdx
         

	   
	   REAL , DIMENSION( ims:ime , kms:kme , jms:jme )   :: px,py,pz,altw

      !  REAL    :: ur, ul, ub, vr, vl, vb
			

    LOGICAL :: specified
	
     
	 


 
   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.



 
   
 	   px = 0.0
	   py = 0.0 
	   pz = 0.0
	   
!write(6,*)   '====  in calculate_ptend_part2_fm_small_test  ===='

	  call    calcu_gradx( pp, px,                 &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )



call    calcu_grady(pp, py,                 &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )


call    calcu_gradz( pp, pz,                 &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )


   call average_vertical_w( altw, alt,                 &
                     config_flags,                  &
                         fzm,fzp,                    &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

! 梯度转化为差分
 	   px = -1.0 * px
	   py = -1.0 * py
	   pz = -1.0 * pz

   i_start = its
   j_start = jts
   i_end=MIN(ite,ide-1)
   j_end=MIN(jte,jde-1)
   ktf=MIN(kte,kde-1)


     DO j = j_start, j_end
     DO k = 1, kte-1
     DO i = i_start, i_end
       
  ! 2023 rho(i,k,j)*

		
       tendency(i,k,j) =  tendency(i,k,j)  +  rcv* (                                  &  
	   	                   rho(i,k,j)*0.5* ( w(i,k,j)*altw(i,k,j)*pz(i,k,j) + w(i,k+1,j)*altw(i,k+1,j)*pz(i,k+1,j) )     &
                      	 - rhob(i,k,j) * g * 0.5* (  w(i,k,j) + w(i,k+1,j)   )           &  				   
	                                          )

      ENDDO
      ENDDO

      ENDDO 

  


	
   
!  x (u) advection

   i_start = its
   j_start = jts
   i_end=MIN(ite,ide-1)
   j_end=MIN(jte,jde-1)

   IF ( (config_flags%open_xs .or. specified) .and. its == ids ) i_start = its+1
   IF ( (config_flags%open_xe .or. specified) .and. ite == ide )  i_end = MIN(ite,ide-2)  



	   DO j = j_start, j_end
	   
     DO k = 1, kte-1
     DO i = i_start, i_end 
	 
	  tendency(i,k,j) =  tendency(i,k,j)  +      &
	 
	       rcv*rho(i,k,j)* 0.25*(                                                 &
		                       u(i+1,k,j) * ( alt(i+1,k,j) + alt(i,k,j) )   * px(i+1,k,j)      &
							 + u(i,k,j)   *  ( alt(i,k,j) + alt(i-1,k,j) )  * px(i,k,j)     )        
	 
	 
	 ENDDO
     ENDDO
						
   ENDDO


   
!  y (v) advection

   i_start = its
   i_end=MIN(ite,ide-1)  
   
   j_start = jts
   j_end=MIN(jte,jde-1)


 
   IF ( (config_flags%open_ys .or. specified) .and. jts == jds ) j_start = jts+1
   IF ( (config_flags%open_ye .or. specified) .and. jte == jde ) j_end = MIN(jte,jde-2) 

 


	 DO j = j_start, j_end
	   
     DO k = 1, kte-1
     DO i = i_start, i_end




	  tendency(i,k,j) =  tendency(i,k,j)  +      &
	 	       rcv*rho(i,k,j)* 0.25*(                                                 &
		                       v(i,k,j+1) * ( alt(i,k,j+1) + alt(i,k,j) )   * py(i,k,j+1)          & 
							 + v(i,k,j)   *  ( alt(i,k,j)  + alt(i,k,j-1) ) * py(i,k,j)      )        
	 
	 
	 ENDDO
     ENDDO
						
   ENDDO
   
   
   
 ! 20230209 测试扰动边界必要性。  
if(.false.) then 
   
   
  

   !  lateral open boundary conditions,
!  start with north and south (y) boundaries

   i_start = its
   i_end=MIN(ite,ide-1)

   !  south

   IF ( (config_flags%open_ys) .and. jts == jds ) THEN

     j=jts

       DO k = kts, kte-1
     DO i = i_start, i_end
	   
	   ! 当vb大于0时，该倾向为0 
	   
     !    vb =.5*( v(i,k,j+1)  + v(i,k,j) )   
      !   vl=amin1(vb,0.)
      !   ph_tend(i,k,j)=ph_tend(i,k,j)-rdy*(c1f(k)*mut(i,j)+c2f(k))*(      &
      !                        +vl*(ph_old(i,k,j+1)-ph_old(i,k,j)))
							  
							  
							  
        tendency(i,k,j) =  tendency(i,k,j)  +      &
	 
	       rcv*rho(i,k,j)* 0.5* v(i,k,j+1) * ( alt(i,k,j+1) + alt(i,k,j) )   * py(i,k,j+1)    

 



!	  tendency(i,k,j) =  tendency(i,k,j)  +      &
!	 	       rcv*rho(i,k,j)* 0.25*(                                                 &
!		                       v(i,k,j+1) * ( alt(i,k,j+1) + alt(i,k,j) )   * py(i,k,j+1)          & 
!							 + v(i,k,j)   *  ( alt(i,k,j)  + alt(i,k,j-1) ) * py(i,k,j)      )        
	 		   
							  
       ENDDO
     ENDDO

   END IF

   


   ! north

   IF ( (config_flags%open_ye) .and. jte == jde ) THEN

     j=jte-1

       DO k = kts, kte-1 
     DO i = i_start, i_end
	   
	   
       ! vb=.5*( fnm(kz)*(v(i,kz  ,j+1)+v(i,kz  ,j))   &
        !       +fnp(kz)*(v(i,kz-1,j+1)+v(i,kz-1,j)) )
			   
		!    vb =.5*( v(i,k,j+1)  + v(i,k,j) )  	   
			   
       ! vr=amax1(vb,0.)
		
        tendency(i,k,j) =  tendency(i,k,j)  +      &
	 	       rcv*rho(i,k,j)* 0.5*   v(i,k,j) *  ( alt(i,k,j)  + alt(i,k,j-1) ) * py(i,k,j) 
       ENDDO
     ENDDO

   END IF
   
   
   
   
   
      !  now the east and west (y) boundaries

   j_start = jts
   j_end=MIN(jte,jde-1)

   !  west

   IF ( (config_flags%open_xs) .and. its == ids ) THEN

     i=its

	   DO j = j_start, j_end
     DO k = kts, kte-1 
	 
       !  ub =.5*( u(i+1,k,j)   + u(i,k,j) )    
       !  ul=amin1(ub,0.)

	  tendency(i,k,j) =  tendency(i,k,j)  +      &
	 	       rcv*rho(i,k,j)* 0.5* u(i+1,k,j) * ( alt(i+1,k,j) + alt(i,k,j) )   * px(i+1,k,j)   

       ENDDO

        
     ENDDO

   END IF

   ! east

   IF ( (config_flags%open_xe) .and. ite == ide ) THEN

     i = ite-1

     DO j = j_start, j_end
     DO k = kts, kte-1 
	 
       !   ub =.5*( u(i+1,k,j)   + u(i,k,j) )  
			   
      !  ur=amax1(ub,0.)


	  tendency(i,k,j) =  tendency(i,k,j)  +      &
	 	       rcv*rho(i,k,j)* 0.5* u(i,k,j)  *  ( alt(i,k,j) + alt(i-1,k,j) )  * px(i,k,j)    

       ENDDO



     ENDDO

   END IF




 

 endif

 
 

END SUBROUTINE calculate_ptend_part2_fm_small_test


!!!!!

SUBROUTINE   calculate_ptend_part2_fm_small_tao2( pp, tendency,                 &
                        u, v, w,sgmdot,                       &
						rho,rhob,         &
						alt,             &
                     config_flags,                  &
                        msft,zs,zh,fzm,fzp,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   ! Input data
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte
											  
											  
		REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: pp

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: u,     &
                                                                      v, w,&
                                                                      sgmdot
																	  
																	  
	    REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: rho,rhob,alt  ! ,pb
		

		 
		 
																		 
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      !   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: mut
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT( INOUT) :: tendency

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,zs


         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh,fzm,fzp

         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   ! Local data
   
   !		 REAL , DIMENSION( ims:ime , kms:kme , jms:jme )   ::  alt 
   
         INTEGER :: i, j, k, ktf
         INTEGER :: i_start, i_end, j_start, j_end
      !   INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
      !   INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
       !  INTEGER :: jp1, jp0, jtmp
       !  REAL :: c1,c2,mrdx
         

	   
	   REAL , DIMENSION( ims:ime , kms:kme , jms:jme )   :: px,py,pz,altw

      !  REAL    :: ur, ul, ub, vr, vl, vb
			

    LOGICAL :: specified
	
     
	 


 
   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.



 
   
 	   px = 0.0
	   py = 0.0 
	   pz = 0.0
	   
	  !!!write(6,*)   '====in calculate_ptend_part2_fm===='
	  call    calcu_gradx( pp, px,                 &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )



call    calcu_grady(pp, py,                 &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )


call    calcu_gradz( pp, pz,                 &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )


   call average_vertical_w( altw, alt,                 &
                     config_flags,                  &
                         fzm,fzp,                    &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

! 梯度转化为差分
 	   px = -1.0 * px
	   py = -1.0 * py
	   pz = -1.0 * pz

   i_start = its
   j_start = jts
   i_end=MIN(ite,ide-1)
   j_end=MIN(jte,jde-1)
   ktf=MIN(kte,kde-1)


     DO j = j_start, j_end
     DO k = 1, kte-1
     DO i = i_start, i_end
       

		
       tendency(i,k,j) =  tendency(i,k,j)  +  rcv* (                                  &  
	   	                   rho(i,k,j)*0.5* ( w(i,k,j)*altw(i,k,j)*pz(i,k,j) + w(i,k+1,j)*altw(i,k+1,j)*pz(i,k+1,j) )     &
        !              	 - rhob(i,k,j) * g * 0.5* (  w(i,k,j) + w(i,k+1,j)   )           &  				   
	                                          )

      ENDDO
      ENDDO

      ENDDO 

  


	
   
!  x (u) advection

   i_start = its
   j_start = jts
   i_end=MIN(ite,ide-1)
   j_end=MIN(jte,jde-1)

   IF ( (config_flags%open_xs .or. specified) .and. its == ids ) i_start = its+1
   IF ( (config_flags%open_xe .or. specified) .and. ite == ide )  i_end = MIN(ite,ide-2)  



	   DO j = j_start, j_end
	   
     DO k = 1, kte-1
     DO i = i_start, i_end 
	 
	  tendency(i,k,j) =  tendency(i,k,j)  +      &
	 
	       rcv*rho(i,k,j)* 0.25*(                                                 &
		                       u(i+1,k,j) * ( alt(i+1,k,j) + alt(i,k,j) )   * px(i+1,k,j)      &
							 + u(i,k,j)   *  ( alt(i,k,j) + alt(i-1,k,j) )  * px(i,k,j)     )        
	 
	 
	 ENDDO
     ENDDO
						
   ENDDO


   
!  y (v) advection

   i_start = its
   i_end=MIN(ite,ide-1)  
   
   j_start = jts
   j_end=MIN(jte,jde-1)


 
   IF ( (config_flags%open_ys .or. specified) .and. jts == jds ) j_start = jts+1
   IF ( (config_flags%open_ye .or. specified) .and. jte == jde ) j_end = MIN(jte,jde-2) 

 


	 DO j = j_start, j_end
	   
     DO k = 1, kte-1
     DO i = i_start, i_end




	  tendency(i,k,j) =  tendency(i,k,j)  +      &
	 	       rcv*rho(i,k,j)* 0.25*(                                                 &
		                       v(i,k,j+1) * ( alt(i,k,j+1) + alt(i,k,j) )   * py(i,k,j+1)          & 
							 + v(i,k,j)   *  ( alt(i,k,j)  + alt(i,k,j-1) ) * py(i,k,j)      )        
	 
	 
	 ENDDO
     ENDDO
						
   ENDDO
   
 

 
 

END SUBROUTINE calculate_ptend_part2_fm_small_tao2
!-----------------------------------------------------------------------

SUBROUTINE calc_p_damp( p, pm1,      &
                       smdiv,             &
                       step,        &
                       ids, ide, jds, jde, kds, kde, &
                       ims, ime, jms, jme, kms, kme, &
                       its,ite, jts,jte, kts,kte    )

  IMPLICIT NONE  ! religion first

! declarations for the stuff coming in

  INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
  INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
  INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte

  INTEGER,      INTENT(IN   )    :: step



  REAL, DIMENSION(ims:ime, kms:kme, jms:jme),INTENT(INOUT) :: p,pm1


  REAL,                                       INTENT(IN   ) :: smdiv



! local variables

  INTEGER :: i, j, k
  INTEGER :: i_start, i_end, j_start, j_end, k_start, k_end
  REAL    :: ptmp

!<DESCRIPTION>
!
!  For the nonhydrostatic option,
!  calc_p_rho computes the perturbation inverse density and
!  perturbation pressure from the hydrostatic relation and the
!  linearized equation of state, respectively.
!
!  For the hydrostatic option,
!  calc_p_rho computes the perturbation pressure, perturbation density,
!  and perturbation geopotential
!  from the vertical coordinate definition, linearized equation of state
!  and the hydrostatic relation, respectively.
!
!  forward weighting of the pressure (divergence damping) is also
!  computed here.
!
!</DESCRIPTION>

   i_start = its
   i_end   = min(ite,ide-1)
   j_start = jts
   j_end   = min(jte,jde-1)
   k_start = kts
   k_end = min(kte,kde-1)

     
!  divergence damping setup
     IF (step == 0) then   ! we're initializing small timesteps
       DO j=j_start, j_end
       DO k=k_start, k_end
       DO i=i_start, i_end
         pm1(i,k,j)=p(i,k,j)
       ENDDO
       ENDDO
       ENDDO
     ELSE                     ! we're in the small timesteps 
       DO j=j_start, j_end    ! and adding div damping component
       DO k=k_start, k_end
       DO i=i_start, i_end
         ptmp = p(i,k,j)
         p(i,k,j) = p(i,k,j) + smdiv*(p(i,k,j)-pm1(i,k,j))
         pm1(i,k,j) = ptmp
       ENDDO
       ENDDO
       ENDDO
     END IF
END SUBROUTINE calc_p_damp

   SUBROUTINE set_w_surface_fm_small( config_flags,              &
                            w, ht, u, v,    &
							rdx, rdy,        &
                            msft,                               &
                            ids, ide, jds, jde, kds, kde,                &
                            ims, ime, jms, jme, kms, kme,                &
                            its, ite, jts, jte, kts, kte                 )
     implicit none

     TYPE( grid_config_rec_type ) config_flags

     INTEGER ,       INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                      ims, ime, jms, jme, kms, kme, &
                                      its, ite, jts, jte, kts, kte

     REAL ,       INTENT(IN   )  :: rdx, rdy   !, cf1, cf2, cf3

     REAL , DIMENSION(  ims:ime , kms:kme, jms:jme ) ,                      &
                                                 INTENT(IN   ) ::  u,       &
                                                                   v

     REAL , DIMENSION(  ims:ime , kms:kme, jms:jme ) ,                      &
                                              INTENT(INOUT) ::  w

     REAL , DIMENSION(  ims:ime , jms:jme ) , INTENT(IN   ) ::  ht,         &
                                                                msft
 
 


     INTEGER :: i,j,k
     INTEGER :: ip1,im1,jp1,jm1
     INTEGER :: ip1_limit,im1_limit,jp1_limit,jm1_limit

!  set kinematic lower boundary condition on W

!  Comments on directional map scale factors:
!  Chain rule: if Z=Z(X,Y) [true at the surface] then
!  dZ/dt = dZ/dX * dX/dt + dZ/dY * dY/dt, U=dX/dt, V=dY/dt
!  using capitals to denote actual values
!  in mapped values, u=U, v=V, z=Z, 1/dX=mx/dx, 1/dY=my/dy
!     => w = dz/dt = mx u dz/dx + my v dz/dy
!  [where dz/dx is just the surface height change between x
!   gridpoints, and dz/dy is the change between y gridpoints]
!  [NB - cf1, cf2 and cf3 do vertical weighting of u or v values
!   nearest the surface]

! get indices for points next to edge of domain

     jm1_limit = jds        ! No periodic BC's
     jp1_limit = jde-1
     im1_limit = ids
     ip1_limit = ide-1

     IF ( config_flags%periodic_x ) THEN
       im1_limit = ids-1
       ip1_limit = ide
     ENDIF

     IF ( config_flags%periodic_y ) THEN
       jm1_limit = jds-1
       jp1_limit = jde
     ENDIF

     DO j = jts,min(jte,jde-1)
       jm1 = max(j-1, jm1_limit)
       jp1 = min(j+1, jp1_limit)
     DO i = its,min(ite,ide-1)
       im1 = max(i-1, im1_limit)
       ip1 = min(i+1, ip1_limit)

         w(i,1,j)=  msft(i,j)*                              &
                  .5*rdy*(                                   &
                           (ht(i,jp1)-ht(i,j  ))             &
          *(v(i,1,j+1))    &
                          +(ht(i,j  )-ht(i,jm1))             &
          *(v(i,1,j  ))  ) &
                    +msft(i,j)*                             &
                  .5*rdx*(                                   &
                           (ht(ip1,j)-ht(i,j  ))             &
          *(u(i+1,1,j))    &
                          +(ht(i  ,j)-ht(im1,j))             &
          *(u(i  ,1,j))  )

      ENDDO
      ENDDO

  

  END SUBROUTINE set_w_surface_fm_small 
  

END MODULE module_small_step_em
