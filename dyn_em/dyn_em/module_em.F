
!WRF:MODEL_LAYER:DYNAMICS
!

MODULE module_em

   USE module_model_constants
   
   USE module_advect_em, only: advect_u_fm, advect_v_fm, advect_w_fm, calcu_divg,calcu_divg_pres,calculate_sgmdot,   &
              advect_scalar_fm2,advect_scalar_pd_fm,advect_scalar_mono_fm,calcu_divg_p  		! hanshenyou  2023 template sub  
   
   USE module_big_step_utilities_em, only: grid_config_rec_type, calculate_full, couple_momentum, calc_mu_uv, calc_ww_cp, &
        calc_cq, calc_alt, calc_php, set_tend, rhs_ph, &
    horizontal_pressure_gradient_fm, pg_buoy_w_fm,    &                ! hanshenyou 
	 w_damp, perturbation_coriolis, coriolis, curvature, horizontal_diffusion, &
        horizontal_diffusion_3dmp, vertical_diffusion_u, &
 vertical_diffusion_v, vertical_diffusion, vertical_diffusion_3dmp, sixth_order_diffusion, rk_rayleigh_damp, &
        theta_relaxation, vertical_diffusion_mp, zero_tend, zero_tend2d    &   !hanshenyou
		,horizontal_pressure_gradient_fm_base   &   !hanshenyou 
	     , coriolis_fm,calculate_ptend_part2_fm_fullp, calculate_ptend_part2_fm_test,   &  !  , vortx,vorty,vortz        ! hanshenyou  
           vertical_diffusion_rho ! hanshenyou


   USE module_state_description, only: param_first_scalar, p_qr, p_qv, p_qc, p_qg, p_qi, p_qs, tiedtkescheme,ntiedtkescheme, heldsuarez, &
        positivedef, gdscheme, g3scheme, gfscheme, kfetascheme, mskfscheme, monotonic, wenopd_scalar, weno_scalar, weno_mom

   USE module_damping_em, only: held_suarez_damp

   USE module_dm

   USE module_llxy

   USE module_domain, ONLY : domain, get_ijk_from_grid

   USE module_configure, ONLY: grid_config_rec_type, model_config_rec, model_to_grid_config_rec

  USE module_bc   ! add 20230603 for bc3d
  
CONTAINS

!------------------------------------------------------------------------

SUBROUTINE rk_step_prep  ( config_flags, rk_step,           &
                           u, v, w, t, ph, mu,              &
                           c1h, c2h, c1f, c2f, moist,       &
                           ru, rv, rw, ww, php, alt,        &
                           muu, muv,                        &
                           mub, mut, phb, pb, p, al, alb,   &
                           cqu, cqv, cqw,cqx,                   &
                           msfux, msfuy,                    &
                           msfvx, msfvx_inv, msfvy,         &
                           msftx, msfty,                    &
                           fnm, fnp, dnw, rdx, rdy,         &
                           n_moist,                         &
                           ids, ide, jds, jde, kds, kde,    &
                           ims, ime, jms, jme, kms, kme,    &
                           its, ite, jts, jte, kts, kte    )

   IMPLICIT NONE


   !  Input data.

   TYPE(grid_config_rec_type   ) ,   INTENT(IN   ) :: config_flags

   INTEGER ,       INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                    ims, ime, jms, jme, kms, kme, &
                                    its, ite, jts, jte, kts, kte

   INTEGER ,       INTENT(IN   ) :: n_moist, rk_step

   REAL ,          INTENT(IN   ) :: rdx, rdy

   REAL , DIMENSION(  ims:ime , kms:kme, jms:jme ) ,                      &
                                               INTENT(IN   ) ::  u,       &
                                                                 v,       &
                                                                 w,       &
                                                                 t,       &
                                                                 ph,      &
                                                                 phb,     &
                                                                 pb,      &
                                                                 al,      &
                                                                 alb

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme  ) ,                     &
                                               INTENT(  OUT) ::  ru,      &
                                                                 rv,      &
                                                                 rw,      &
                                                                 ww,      &
                                                                 php,     &
                                                                 cqu,     &
                                                                 cqv,     &
                                                                 cqw,cqx,     &
                                                                 alt

   REAL , DIMENSION(  ims:ime , kms:kme, jms:jme ) ,                      &
                                               INTENT(IN   ) ::  p
                                                                 



   REAL , DIMENSION( ims:ime, kms:kme, jms:jme, n_moist ), INTENT(   IN) :: &
                                                           moist

   REAL , DIMENSION( ims:ime , jms:jme ) ,    INTENT(IN   ) :: msftx,     &
                                                               msfty,     &
                                                               msfux,     &
                                                               msfuy,     &
                                                               msfvx,     &
                                                               msfvx_inv, &
                                                               msfvy,     &
                                                               mu,        &
                                                               mub

   REAL , DIMENSION( ims:ime , jms:jme ) ,    INTENT(  OUT) :: muu,    &
                                                               muv,    &
                                                               mut

   REAL , DIMENSION( kms:kme ) ,    INTENT(IN   ) :: fnm, fnp, dnw
   REAL , DIMENSION( kms:kme ) ,    INTENT(IN   ) :: c1h, c2h, c1f, c2f

   integer :: k


!<DESCRIPTION>
!
!  rk_step_prep prepares a number of diagnostic quantities
!  in preperation for a Runge-Kutta timestep.  subroutines called
!  by rk_step_prep calculate
!
!  (1) total column dry air mass (mut, call to calculate_full)
!
!  (2) total column dry air mass at u and v points
!      (muu, muv, call to calculate_mu_uv)
!
!  (3) mass-coupled velocities for advection
!      (ru, rv, and rw, call to couple_momentum)
!
!  (4) omega (call to calc_ww_cp)
!
!  (5) moisture coefficients (cqu, cqv, cqw, call to calc_cq)
!
!  (6) inverse density (alt, call to calc_alt)
!
!  (7) geopotential at pressure points (php, call to calc_php)
!
!</DESCRIPTION>

   CALL calculate_full( mut, mub, mu,             &
                        ids, ide, jds, jde, 1, 2, &
                        ims, ime, jms, jme, 1, 1, &
                        its, ite, jts, jte, 1, 1 )

   CALL calc_mu_uv ( config_flags,                  &
                     mu, mub, muu, muv,             &
                     ids, ide, jds, jde, kds, kde,  &
                     ims, ime, jms, jme, kms, kme,  &
                     its, ite, jts, jte, kts, kte  )

   CALL couple_momentum( muu, ru, u, msfuy,             &
                         muv, rv, v, msfvx, msfvx_inv,  &
                         mut, rw, w, msfty,             &
                         c1h, c2h, c1f, c2f,            &
                         ids, ide, jds, jde, kds, kde,  &
                         ims, ime, jms, jme, kms, kme,  &
                         its, ite, jts, jte, kts, kte  )

!  new call, couples V with mu, also has correct map factors.  WCS, 3 june 2001
   CALL calc_ww_cp ( u, v, mu, mub, c1h, c2h, ww,     &
                     rdx, rdy, msftx, msfty,          &
                     msfux, msfuy, msfvx, msfvx_inv,  &
                     msfvy, dnw,                      &
                     ids, ide, jds, jde, kds, kde,    &
                     ims, ime, jms, jme, kms, kme,    &
                     its, ite, jts, jte, kts, kte    )

   CALL calc_cq ( moist, cqu, cqv, cqw, cqx, n_moist, &
                  ids, ide, jds, jde, kds, kde,  &
                  ims, ime, jms, jme, kms, kme,  &
                  its, ite, jts, jte, kts, kte  )

   CALL calc_alt ( alt, al, alb,                 &
                   ids, ide, jds, jde, kds, kde, &
                   ims, ime, jms, jme, kms, kme, &
                   its, ite, jts, jte, kts, kte )

   CALL calc_php ( php, ph, phb,                 &
                   ids, ide, jds, jde, kds, kde, &
                   ims, ime, jms, jme, kms, kme, &
                   its, ite, jts, jte, kts, kte )

END SUBROUTINE rk_step_prep

!------------------------
!------------------------------------------------------------------------

! 2023 现在只计算cqu，cqv，cqw，其他用了calc_p_rho_fm_new 
SUBROUTINE rk_step_prep_fm( config_flags, alt,        &
                           al, alb,   &
	                       cqu, cqv, cqw, cqx,                  & 
						   moist,       &
				          n_moist,                         &
                           ids, ide, jds, jde, kds, kde,    &
                           ims, ime, jms, jme, kms, kme,    &
                           its, ite, jts, jte, kts, kte    )

   IMPLICIT NONE


   !  Input data.

   TYPE(grid_config_rec_type   ) ,   INTENT(IN   ) :: config_flags

   INTEGER ,       INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                    ims, ime, jms, jme, kms, kme, &
                                    its, ite, jts, jte, kts, kte



   REAL , DIMENSION(  ims:ime , kms:kme, jms:jme ) ,                      &
                                               INTENT(IN   ) ::  al,      &
                                                                 alb

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme  ) ,                     &
                                               INTENT(  OUT) ::    alt,   &
                                                                 cqu,     &
                                                                 cqv,     &
                                                                 cqw,cqx  
																 
   INTEGER ,          INTENT(IN   ) :: n_moist
   

   REAL, DIMENSION( ims:ime, kms:kme , jms:jme , n_moist ), INTENT(IN   ) :: moist 


  INTEGER :: i, j, k
  INTEGER :: i_start, i_end, j_start, j_end, k_start, k_end
 

 !  CALL calc_alt ( alt, al, alb,                 &
 !                  ids, ide, jds, jde, kds, kde, &
 !                  ims, ime, jms, jme, kms, kme, &
 !                  its, ite, jts, jte, kts, kte )
				   

if(.true.) then  ! 测试cqu等的影响。 20230416
				   
   CALL calc_cq ( moist, cqu, cqv, cqw, cqx,n_moist, &
                  ids, ide, jds, jde, kds, kde,  &
                  ims, ime, jms, jme, kms, kme,  &
                  its, ite, jts, jte, kts, kte  )	
				  
				  


 
else
	cqu=1.0
	cqv=1.0
	cqw=0.0
	cqx = 0.0 
	
endif






 

 			  


END SUBROUTINE rk_step_prep_fm
!-------------------------	
 			 
SUBROUTINE rk_tendency ( config_flags, rk_step,                           &
                         u_tend, v_tend, w_tend, pfm_tend, rhofm_tend,      &
                         u_tendf, v_tendf, w_tendf, pfm_tendf, rhofm_tendf, &
						 u_save, v_save, w_save, pfm_save, rhofm_save,        &    ! 20231223 zero that 
                         RTHFTEN,                                           &
                   !      ru, rv, rw, ww,                                    &
                         u, v, w, rhofm, pfm,                                  &
						 pfm_base,rhofm_base,                          &                          
   						 sigmdot,                                         & 
                         u_old, v_old, w_old, rhofm_old, pfm_old,              &
                         h_diabatic, phb,t_init,                          &
                         mu, mut, muu, muv, mub, c1h, c2h, c1f, c2f,      &
                        al, alt, p,                    &                          ! hanshenyou here p is disturb p 
						 pb, php, cqu, cqv, cqw,cqx,              &
                         u_base, v_base, t_base, qv_base, z_base,         &
						 msft, zs,zh,                                          &     
						 fzm,fzp, ztop ,                                     &    
                        msfux, msfuy, msfvx, msfvx_inv,                  &
                        msfvy, msftx, msfty,                             &
                         clat, f, e, sina, cosa,                          &
                         fnm, fnp, rdn, rdnw,dnw,                             &
                         dt, rdx, rdy, khdif, kvdif, xkmhd, xkhh,         &
                         diff_6th_opt, diff_6th_factor,                   &
                         adv_opt,                                         &
                         dampcoef,zdamp,damp_opt,rad_nudge,               &
                         cf1, cf2, cf3, cfn, cfn1, n_moist,               &
                         non_hydrostatic, top_lid,                        &
           !              u_frame, v_frame,                                &
                         ids, ide, jds, jde, kds, kde,                    &
                         ims, ime, jms, jme, kms, kme,                    &
                         its, ite, jts, jte, kts, kte,                    &
                         max_vert_cfl, max_horiz_cfl)

   IMPLICIT NONE

   !  Input data.

   TYPE(grid_config_rec_type)    ,           INTENT(IN   ) :: config_flags

   INTEGER ,               INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                            ims, ime, jms, jme, kms, kme, &
                                            its, ite, jts, jte, kts, kte

   LOGICAL ,               INTENT(IN   ) :: non_hydrostatic, top_lid

   INTEGER ,               INTENT(IN   ) :: n_moist, rk_step


   REAL , DIMENSION( ims:ime , kms:kme, jms:jme  ) ,              &
                                        INTENT(INOUT  ) ::  sigmdot  

   REAL , DIMENSION( ims:ime , kms:kme, jms:jme  ) ,              &
                                        INTENT(IN   ) :: u,       &
                                                         v,       &
                                                         w,       &
                                                         rhofm,pfm,pfm_base,rhofm_base,      &    
                                                         u_old,   &
                                                         v_old,   &
                                                         w_old,   &
                                                         rhofm_old,   &
                                                         pfm_old,  &
                                                         phb,     &
                                                         al,      &
                                                         alt,     &
                                                         p,       &
                                                         pb,      &
                                                         php,     &
                                                         cqu,     &
                                                         cqv,     &
														 cqw,cqx,    &
                                                         t_init,  &
                                                         xkmhd,   &
                                                         xkhh,    &
                                                         h_diabatic

   REAL , DIMENSION( ims:ime , kms:kme, jms:jme  ) ,              &
                                        INTENT(OUT  ) :: u_tend, &
                                                         v_tend, &
                                                         w_tend, &
                                                         rhofm_tend,  &
                                                         pfm_tend, &
                                                         RTHFTEN 

   REAL , DIMENSION( ims:ime , kms:kme, jms:jme  ) ,               &
                                        INTENT(INOUT) :: u_tendf, &
                                                         v_tendf, &
                                                         w_tendf, &
                                                         rhofm_tendf,  &
                                                         pfm_tendf 
														 
														 
   REAL , DIMENSION( ims:ime , kms:kme, jms:jme  ) ,               &
                                        INTENT(INOUT) ::  u_save, v_save, w_save,  &
										                  pfm_save, rhofm_save     ! 20231223 zero that 

 !  REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(  OUT) :: mu_tend, &
 !                                                                   mu_save



   REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) ::  msft,zs     ! hanshenyou 
    
   REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msfux,   &
                                                                    msfuy,   &
                                                                    msfvx,   &
                                                                    msfvx_inv,   &
                                                                    msfvy,   &
                                                                    msftx,   &
                                                                    msfty,   &
                                                                    clat,    &
                                                                    f,       &
                                                                    e,       &
                                                                    sina,    &
                                                                    cosa,    &
                                                                    mu,      &
                                                                    mut,     &
                                                                    mub,     &
                                                                    muu,     &
                                                                    muv


   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: fzm,     &
                                                                  fzp, &
                                                                  zh
																  
																  
   REAL ,                                      INTENT(IN   ) ::   ztop
   
	   
   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: fnm,     &
                                                                  fnp,     &
                                                                  rdn,     &
                                                                  rdnw,dnw,   &
                                                                  u_base,  &
                                                                  v_base,  &
                                                                  t_base,  &
                                                                  qv_base, &
                                                                  z_base,  &
                                                                  c1h,     &
                                                                  c2h,     &
                                                                  c1f,     &
                                                                  c2f

   REAL ,                                      INTENT(IN   ) :: rdx,     &
                                                                rdy,     &
                                                                dt,      &
                                                      !          u_frame, &
                                                       !         v_frame, &
                                                                khdif,   &
                                                                kvdif
   INTEGER, INTENT( IN ) :: diff_6th_opt
   REAL,    INTENT( IN ) :: diff_6th_factor
   INTEGER, INTENT( IN ) :: adv_opt

   INTEGER, INTENT( IN ) :: damp_opt, rad_nudge

   REAL, INTENT( IN ) :: zdamp, dampcoef

   REAL, INTENT( OUT ) :: max_horiz_cfl
   REAL, INTENT( OUT ) :: max_vert_cfl

   REAL    :: kdift, khdq, kvdq, cfn, cfn1, cf1, cf2, cf3
   INTEGER :: i,j,k
   INTEGER :: time_step
   
   
   REAL , DIMENSION( kms:kme ) :: sgmh,dsgm,dsgmh

   INTEGER,dimension(3)  :: dimmax 	
		
   !20220703 
    !  REAL , DIMENSION( ims:ime , kms:kme, jms:jme  ) ::  pp,rhop

!<DESCRIPTION>
!
!  rk_tendency computes the large-timestep tendency terms in the
!  momentum, thermodynamic (theta), and geopotential equations.  
!  These terms include:
!
!  (1) advection (for u, v, w, theta - calls to advect_u, advect_v,
!                 advect_w, and advact_scalar).
!
!  (2) geopotential equation terms (advection and "gw" - call to rhs_ph).
!
!  (3) buoyancy term in vertical momentum equation (call to pg_buoy_w).
!
!  (4) Coriolis and curvature terms in u,v,w momentum equations
!      (calls to subroutines coriolis, curvature)
!
!  (5) 3D diffusion on coordinate surfaces.
!
!</DESCRIPTION>

   CALL zero_tend ( u_tend,                      &
                    ids, ide, jds, jde, kds, kde, &
                    ims, ime, jms, jme, kms, kme, &
                    its, ite, jts, jte, kts, kte )

   CALL zero_tend ( v_tend,                      &
                    ids, ide, jds, jde, kds, kde, &
                    ims, ime, jms, jme, kms, kme, &
                    its, ite, jts, jte, kts, kte )

   CALL zero_tend ( w_tend,                      &
                    ids, ide, jds, jde, kds, kde, &
                    ims, ime, jms, jme, kms, kme, &
                    its, ite, jts, jte, kts, kte )

   CALL zero_tend ( rhofm_tend,                       &
                    ids, ide, jds, jde, kds, kde, &
                    ims, ime, jms, jme, kms, kme, &
                    its, ite, jts, jte, kts, kte )

   CALL zero_tend ( pfm_tend,                      &
                    ids, ide, jds, jde, kds, kde, &
                    ims, ime, jms, jme, kms, kme, &
                    its, ite, jts, jte, kts, kte )
			

 
					
!   debug  template  zero_tend only on its-ite and so on  not ims to ime  20230107 

  ! 20231224 hanshenyou  				
   CALL zero_tend ( u_save,                      &
                    ids, ide, jds, jde, kds, kde, &
                    ims, ime, jms, jme, kms, kme, &
                    its, ite, jts, jte, kts, kte )

   CALL zero_tend ( v_save,                      &
                    ids, ide, jds, jde, kds, kde, &
                    ims, ime, jms, jme, kms, kme, &
                    its, ite, jts, jte, kts, kte )

   CALL zero_tend ( w_save,                      &
                    ids, ide, jds, jde, kds, kde, &
                    ims, ime, jms, jme, kms, kme, &
                    its, ite, jts, jte, kts, kte )

   CALL zero_tend ( rhofm_save,                       &
                    ids, ide, jds, jde, kds, kde, &
                    ims, ime, jms, jme, kms, kme, &
                    its, ite, jts, jte, kts, kte )

   CALL zero_tend ( pfm_save,                      &
                    ids, ide, jds, jde, kds, kde, &
                    ims, ime, jms, jme, kms, kme, &
                    its, ite, jts, jte, kts, kte )


   !  advection tendencies
   CALL nl_get_time_step ( 1, time_step )

 
 
	 
	!  hanshenyou  energy form  
     CALL advect_u_fm (u_tend,  &
                       u, v, w,                   &
                       time_step, config_flags,  &
			           msfux, msfuy, msfvx, msfvy,    &
                       msft,clat,zs,                &                       ! clat -> xlat 
                       zh,fzm, fzp,                      &
                      rdx, rdy, ztop,                &
                   ids, ide, jds, jde, kds, kde, &
                   ims, ime, jms, jme, kms, kme, &
                   its, ite, jts, jte, kts, kte )
				   
				   

	   

     CALL advect_v_fm ( v_tend, &
                       u, v, w,                   &
                       time_step, config_flags,  &
		               msfux, msfuy, msfvx, msfvy,    &
                       msft,msfty,clat,zs,             &                   ! clat -> xlat 
                       zh,fzm, fzp,                      &
                      rdx, rdy, ztop,                &
                   ids, ide, jds, jde, kds, kde, &
                   ims, ime, jms, jme, kms, kme, &
                   its, ite, jts, jte, kts, kte )
				   


	 
	 
	!   calculate_sgmdot 子程序用于：计算sigma坐标的垂直速度sgmdot	  	  
!_______________________________________
	
! hanshenyou 20211201  sigmdot	
 
	 
    call    calculate_sgmdot( sigmdot,                 &
                        u, v, w,                       &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  ) 
						
						

	 
   ! 这里pfm是全气压

     CALL calcu_divg( pfm, pfm_tend,                 &
    ! CALL calcu_divg_pres( pfm, pfm_tend,                 &
                        u, v, sigmdot,                       &
                     config_flags,                  &
                        msft,zh,zs,fzm,fzp,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

						



	 
	 
	
     CALL calcu_divg( rhofm, rhofm_tend,                 &
                        u, v, sigmdot,                       &
                     config_flags,                  &
                        msft,zh,zs,fzm,fzp,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

       


						
						


	 
         

	 
              pfm_tend  =   cpovcv*pfm_tend   !  20221201 
	 
 
 
	 
 	 
 
 ! 这里p是扰动气压
 ! 扰动气压p改为 pfm全气压，20230724测试全气压，斜压波试验有问题
  !hanshenyou 20231209 ,又改回扰动气压p，测试实例模拟 
  
 if(.false.)then ! 全量，目前适合无地形
 

	CALL horizontal_pressure_gradient_fm( u, u_tend, v, v_tend,            &            
                        pfm,alt,cqu,cqv,                                      &    ! p ? pfm 
                        config_flags,               &
                        msft,zs,                              &
                        zh,fzm, fzp,                          &
                        rdx, rdy,            &
				!		dts,               &
						ztop,                  &
                        ids, ide, jds, jde, kds, kde,         &
                        ims, ime, jms, jme, kms, kme,         &
                        its, ite, jts, jte, kts, kte  ) 
 
	CALL  calculate_ptend_part2_fm_fullp( pfm, p,pfm_tend,                 &
                        u, v, w,sigmdot,                       &
						rhofm,rhofm_base,           &          
						alt,     &
                     config_flags,                  &
                        msft,zs,zh,fzm,fzp,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )
else
     if(.false.) then !若无地形，1D参考态，这一部分应该为0。若有地形则可能不为0。密度没有用参考态，没有线性化？
	 ! 测试有地形，比如山脉波？
!	CALL horizontal_pressure_gradient_fm_base( u, u_tend, v, v_tend,            &            
!                        pfm_base,alt,cqu,cqv,                                      &    ! p ? pfm 
!                        config_flags,               &

						
						
						
	endif 					
	! 适合 1D 参考态， 扰动气压				
	CALL horizontal_pressure_gradient_fm( u, u_tend, v, v_tend,            &            
                        p,alt,cqu,cqv,                                      &    ! p ? pfm 
                        config_flags,               &
                        msft,zs,                              &
                        zh,fzm, fzp,                          &
                        rdx, rdy,            &
				!		dts,               &
						ztop,                  &
                        ids, ide, jds, jde, kds, kde,         &
                        ims, ime, jms, jme, kms, kme,         &
                        its, ite, jts, jte, kts, kte  )


	 
   
	 
	 

	CALL  calculate_ptend_part2_fm_test( p, pfm_tend,                 &
                        u, v, w,sigmdot,                       &
						rhofm,rhofm_base,            &          
						alt,     &
                     config_flags,                  &
                        msft,zs,zh,fzm,fzp,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )				
						
	   !write(6,*) ' after add calculate_ptend_part2_fm_test  max divg pfm_tend(i,k,j)= ', pfm_tend(i,k,j),i,k,j 
	   
	   !  	  dimmax = maxloc(pfm_tend(its:ite,kts:kte, jts:jte)) 
       !  i = dimmax(1) 
       !  k = dimmax(2)
       ! j = dimmax(3) 	  
	
 

	 
 !write(6,*) 'pfm_tend(dimmax(1)-1:dimmax(1)+1 , dimmax(2):dimmax(2)+1,dimmax(3)-2:dimmax(3)+2)  = :' 
 !write(6,*)           pfm_tend(dimmax(1)-1:dimmax(1)+1 , dimmax(2):dimmax(2)+1,dimmax(3)-2:dimmax(3)+2)  
						
 ! 扰动气压p改为 pfm全气压，20230724测试全气压，斜压波试验有问题
    
 !hanshenyou 20231209 ,又改回扰动气压p，测试实例模拟。若为pfm，则为全气压。



endif 


	CALL pg_buoy_w_fm( w, w_tend,            &            
                        p,alt,al,cqw,                                    &
                        config_flags,              &
                        msft,zs,                              &
                        zh,fzm, fzp,                          &
                        rdx, rdy,       &
					!	dts,     &
						ztop,                  &
                        ids, ide, jds, jde, kds, kde,         &
                        ims, ime, jms, jme, kms, kme,         &
                        its, ite, jts, jte, kts, kte  )
 
  !  write(6,*) ' w_tend(ide/3,:,jde/3) : '
  !  write(6,*)   w_tend(ide/3,:,jde/3) 
  !  stop 
	 
   

   ! hanshenyou this function is not finish  w_damp_fm 
   !  CALL w_damp   ( w_tend, max_vert_cfl,             &
   !                   max_horiz_cfl,                    &
   !                   u, v, ww, w, mut, c1f, c2f, rdnw, &
   !                   rdx, rdy, msfux, msfuy, msfvx,    &
   !                   msfvy, dt, config_flags,          &
   !                   ids, ide, jds, jde, kds, kde,     &
   !                   ims, ime, jms, jme, kms, kme,     &
   !                   its, ite, jts, jte, kts, kte     )

    CALL advect_w_fm ( w_tend,  w, &
                       u, v, w,                   &
                       time_step, config_flags,  &
                       msft,clat,zs,    &                   ! clat -> xlat 
                       zh,fzm, fzp,                      &
                      rdx, rdy, ztop,                &
                   ids, ide, jds, jde, kds, kde, &
                   ims, ime, jms, jme, kms, kme, &
                   its, ite, jts, jte, kts, kte )


  	!  dimmax = maxloc(w_tend(its:ite,kts:kte, jts:jte)) 
    !     i = dimmax(1) 
    !     k = dimmax(2)
     !   j = dimmax(3) 	  
	
	 !  write(6,*)  
	  ! write(6,*) ' add advect_w_fm   max w_tend(i,k,j)', w_tend(i,k,j),i,k,j 
	 ! write(6,*) u(i,k,j),v(i,k,j)





     
 
      call  coriolis_fm ( u, v, w, u_tend, v_tend, w_tend, &
                      config_flags,                          &
                      msftx, msfty, msfux, msfuy,            &
                      msfvx, msfvy,                          &
                      f, e, sina, cosa, fzm, fzp,            &
                      ids, ide, jds, jde, kds, kde,          &
                      ims, ime, jms, jme, kms, kme,          &
                      its, ite, jts, jte, kts, kte          )
					  
					  



	 
 
	  

  !   CALL curvature ( ru, rv, rw, u, v, w,            &
   !                    ru_tend,  rv_tend,  rw_tend,    &
   !                    config_flags,                   &
   !                    msfux, msfuy, msfvx, msfvy,     &
    !                   msftx, msfty,                   &
    !                   clat, fnm, fnp, rdx, rdy,       &
     !                  ids, ide, jds, jde, kds, kde,   &
     !                  ims, ime, jms, jme, kms, kme,   &
     !                  its, ite, jts, jte, kts, kte   )


! Damping option added for Held-Suarez test (also uses lw option HELDSUAREZ)
    !  IF (config_flags%ra_lw_physics == HELDSUAREZ) THEN
    !     CALL held_suarez_damp ( u_tend, v_tend,               &   
    !                             u,v,pfm,pb,                     &
     !!                            ids, ide, jds, jde, kds, kde,   &
    !                             ims, ime, jms, jme, kms, kme,   &
    !                             its, ite, jts, jte, kts, kte   )
    !  END IF


!**************************************************************
!
!  Next, the terms that we integrate only with forward-in-time
!  (evaluate with time t variables).
!
!**************************************************************

  forward_step: IF( rk_step == 1 ) THEN
  
     diff_opt9 : IF (config_flags%diff_opt .eq. 9) THEN
	 
	 
if(.false.) then  ! old diffusion 
	 
	 
  call diffusion_2ord_fm( u_tendf, config_flags,              &
                                     u, rdx,                           &
                                     ids, ide, jds, jde, kds, kde,        &
                                     ims, ime, jms, jme, kms, kme,        &
                                     its, ite, jts, jte, kts, kte        )
  
   call diffusion_2ord_fm( w_tendf, config_flags,              &
                                     w, rdx,                           &
                                     ids, ide, jds, jde, kds, kde,        &
                                     ims, ime, jms, jme, kms, kme,        &
                                     its, ite, jts, jte, kts, kte        )
									 
									 
  call diffusion_2ord_fm( rhofm_tendf, config_flags,              &
                                     al, rdx,                           &   !rho 
                                     ids, ide, jds, jde, kds, kde,        &
                                     ims, ime, jms, jme, kms, kme,        &
                                     its, ite, jts, jte, kts, kte        )	

      !   pp = pfm - pfm_base
		 

     ! CALL set_physical_bc3d( pp, 'p', config_flags,               &
     !                         ids, ide, jds, jde, kds, kde,       &
     !                         ims, ime, jms, jme, kms, kme,       &
     !                          its, ite, jts, jte, kts, kte,       &
     !                         its, ite, jts, jte, kts, kte )
							  
  call diffusion_2ord_fm( pfm_tendf, config_flags,              &
                                     p, rdx,                           &
                                     ids, ide, jds, jde, kds, kde,        &
                                     ims, ime, jms, jme, kms, kme,        &
                                     its, ite, jts, jte, kts, kte        )				


else

call diffusion_2ord_fmh( 'm',u_tendf, config_flags,              &
                                     u, rdx,                           &
                                     ids, ide, jds, jde, kds, kde,        &
                                     ims, ime, jms, jme, kms, kme,        &
                                     its, ite, jts, jte, kts, kte        )
  
   call diffusion_2ord_fmh( 'w',w_tendf, config_flags,              &
                                     w, rdx,                           &
                                     ids, ide, jds, jde, kds, kde,        &
                                     ims, ime, jms, jme, kms, kme,        &
                                     its, ite, jts, jte, kts, kte        )
									 
									 
  call diffusion_2ord_fmh( 'm',rhofm_tendf, config_flags,              &
                                     al, rdx,                           &   !rho 
                                     ids, ide, jds, jde, kds, kde,        &
                                     ims, ime, jms, jme, kms, kme,        &
                                     its, ite, jts, jte, kts, kte        )	


							  
  call diffusion_2ord_fmh( 'm',pfm_tendf, config_flags,              &
                                     p, rdx,                           &
                                     ids, ide, jds, jde, kds, kde,        &
                                     ims, ime, jms, jme, kms, kme,        &
                                     its, ite, jts, jte, kts, kte        )
									 
	! vertical_diffusion : 								 

									 
	call diffusion_2ord_fmv( 'u',u_tendf, config_flags,              &
                                     u, rdx,                           &
                                     ids, ide, jds, jde, kds, kde,        &
                                     ims, ime, jms, jme, kms, kme,        &
                                     its, ite, jts, jte, kts, kte        )
  
   call diffusion_2ord_fmv( 'w',w_tendf, config_flags,              &
                                     w, rdx,                           &
                                     ids, ide, jds, jde, kds, kde,        &
                                     ims, ime, jms, jme, kms, kme,        &
                                     its, ite, jts, jte, kts, kte        )
									 
			! name 'p'  vflux(i,0) = 0  20220828 						 
  call diffusion_2ord_fmv( 'm',rhofm_tendf, config_flags,              &
                                     al, rdx,                           &   !rho 
                                     ids, ide, jds, jde, kds, kde,        &
                                     ims, ime, jms, jme, kms, kme,        &
                                     its, ite, jts, jte, kts, kte        )	


							  
  call diffusion_2ord_fmv( 'p',pfm_tendf, config_flags,              &
                                     p, rdx,                           &
                                     ids, ide, jds, jde, kds, kde,        &
                                     ims, ime, jms, jme, kms, kme,        &
                                     its, ite, jts, jte, kts, kte        )		





endif

									 
      END IF diff_opt9
	  
	  
  if (.true.) then 
    diff_opt8 : IF (config_flags%diff_opt .eq. 8) THEN
	
	
      CALL horizontal_diffusion ('u', u-u_old, u_tendf, mut,               &
                                        c1h, c2h, config_flags,         &
                                        msfux, msfuy, msfvx, msfvx_inv, &
                                        msfvy,msftx, msfty,             &
                                        khdif, xkmhd, rdx, rdy,         &
                                        ids, ide, jds, jde, kds, kde,   &
                                        ims, ime, jms, jme, kms, kme,   &
                                        its, ite, jts, jte, kts, kte   )
										
 

        CALL horizontal_diffusion ('v', v-v_old, v_tendf, mut,               &
                                        c1h, c2h, config_flags,         &
                                        msfux, msfuy, msfvx, msfvx_inv, &
                                        msfvy,msftx, msfty,             &
                                        khdif, xkmhd, rdx, rdy,         &
                                        ids, ide, jds, jde, kds, kde,   &
                                        ims, ime, jms, jme, kms, kme,   &
                                        its, ite, jts, jte, kts, kte   )

        CALL horizontal_diffusion ('w', w-w_old, w_tendf, mut,               &
                                        c1f, c2f, config_flags,         &
                                        msfux, msfuy, msfvx, msfvx_inv, &
                                        msfvy,msftx, msfty,             &
                                        khdif, xkmhd, rdx, rdy,         &
                                        ids, ide, jds, jde, kds, kde,   &
                                        ims, ime, jms, jme, kms, kme,   &
                                        its, ite, jts, jte, kts, kte   )
		
	 
			   
        CALL horizontal_diffusion ('p', rhofm-rhofm_old, rhofm_tendf, mut,               &
                                        c1f, c2f, config_flags,         &
                                        msfux, msfuy, msfvx, msfvx_inv, &
                                        msfvy,msftx, msfty,             &
                                        khdif, xkhh, rdx, rdy,         &              ! xkmhd,改为xkhh，与水汽同样乘以pr数。
                                        ids, ide, jds, jde, kds, kde,   &
                                        ims, ime, jms, jme, kms, kme,   &
                                        its, ite, jts, jte, kts, kte   )
										
 

         CALL horizontal_diffusion ('p', pfm-pfm_old, pfm_tendf, mut,               &
                                         c1f, c2f, config_flags,         &
                                         msfux, msfuy, msfvx, msfvx_inv, &
                                         msfvy,msftx, msfty,             &
                                         khdif, xkhh, rdx, rdy,         &
                                         ids, ide, jds, jde, kds, kde,   &
                                         ims, ime, jms, jme, kms, kme,   &
                                         its, ite, jts, jte, kts, kte   )	
	
	
	
	
	
	endif   diff_opt8 
  
  endif
  
  

    diff_opt1 : IF (config_flags%diff_opt .eq. 1) THEN
	
 
	  
   
        CALL horizontal_diffusion ('u', u, u_tendf, mut,               &
                                        c1h, c2h, config_flags,         &
                                        msfux, msfuy, msfvx, msfvx_inv, &
                                        msfvy,msftx, msfty,             &
                                        khdif, xkmhd, rdx, rdy,         &
                                        ids, ide, jds, jde, kds, kde,   &
                                        ims, ime, jms, jme, kms, kme,   &
                                        its, ite, jts, jte, kts, kte   )
										
 

        CALL horizontal_diffusion ('v', v, v_tendf, mut,               &
                                        c1h, c2h, config_flags,         &
                                        msfux, msfuy, msfvx, msfvx_inv, &
                                        msfvy,msftx, msfty,             &
                                        khdif, xkmhd, rdx, rdy,         &
                                        ids, ide, jds, jde, kds, kde,   &
                                        ims, ime, jms, jme, kms, kme,   &
                                        its, ite, jts, jte, kts, kte   )

        CALL horizontal_diffusion ('w', w, w_tendf, mut,               &
                                        c1f, c2f, config_flags,         &
                                        msfux, msfuy, msfvx, msfvx_inv, &
                                        msfvy,msftx, msfty,             &
                                        khdif, xkmhd, rdx, rdy,         &
                                        ids, ide, jds, jde, kds, kde,   &
                                        ims, ime, jms, jme, kms, kme,   &
                                        its, ite, jts, jte, kts, kte   )
		
		       if(.true.) then 
			   
        CALL horizontal_diffusion ('p', al, rhofm_tendf, mut,               &
                                        c1f, c2f, config_flags,         &
                                        msfux, msfuy, msfvx, msfvx_inv, &
                                        msfvy,msftx, msfty,             &
                                        khdif, xkmhd, rdx, rdy,         &              ! xkmhd,改为xkhh，与水汽同样乘以pr数。 2025,khdif no used,xkhh is xkmhd/pr
                                        ids, ide, jds, jde, kds, kde,   &
                                        ims, ime, jms, jme, kms, kme,   &
                                        its, ite, jts, jte, kts, kte   )
										
          endif

        if(.true.) then  ! 2025 hanshenyou for grav2d test 

         CALL horizontal_diffusion ('p', p, pfm_tendf, mut,               &
                                         c1f, c2f, config_flags,         &
                                         msfux, msfuy, msfvx, msfvx_inv, &
                                         msfvy,msftx, msfty,             &
                                         khdif, xkmhd, rdx, rdy,         &
                                         ids, ide, jds, jde, kds, kde,   &
                                         ims, ime, jms, jme, kms, kme,   &
                                         its, ite, jts, jte, kts, kte   )
										
         
			 
		endif 
		
 
   
     !   khdq = 3.*khdif
     !   CALL horizontal_diffusion_3dmp ( 'm', al, rhofm_tendf, mut,            &
     !                                    c1h, c2h,                        &
      !                                   config_flags, rhofm_base,            &
     !                                    msfux, msfuy, msfvx, msfvx_inv,  &
      !                                   msfvy, msftx, msfty,             &
     !                                    khdq , xkhh, rdx, rdy,           &
     !                                    ids, ide, jds, jde, kds, kde,    &
     !                                    ims, ime, jms, jme, kms, kme,    &
     !                                    its, ite, jts, jte, kts, kte    )
	!									 
     !   CALL horizontal_diffusion_3dmp ( 'm', p, pfm_tendf, mut,            &
     !                                    c1h, c2h,                        &
      !                                   config_flags, pfm_base,            &
     !                                    msfux, msfuy, msfvx, msfvx_inv,  &
     !                                    msfvy, msftx, msfty,             &
     !                                    khdq , xkhh, rdx, rdy,           &
     !                                    ids, ide, jds, jde, kds, kde,    &
     !                                    ims, ime, jms, jme, kms, kme,    &
     !                                    its, ite, jts, jte, kts, kte    )										 
										 
										 

       pbl_test : IF (config_flags%bl_pbl_physics .eq. 0) THEN
	 
	 
	 
     ! pbl_test : IF ( .false. ) THEN             ! real test 20231108  
	  
	  
	  
	  
          CALL vertical_diffusion_u ( u, u_tendf, config_flags,      &
                                      u_base, c1h, c2h,               &
                                      alt, muu, rdn, rdnw, kvdif,     &
                                      ids, ide, jds, jde, kds, kde,   &
                                      ims, ime, jms, jme, kms, kme,   &
                                      its, ite, jts, jte, kts, kte   )

          CALL vertical_diffusion_v ( v, v_tendf, config_flags,      &
                                      v_base, c1h, c2h,               &
                                      alt, muv, rdn, rdnw, kvdif,     &
                                      ids, ide, jds, jde, kds, kde,   &
                                      ims, ime, jms, jme, kms, kme,   &
                                      its, ite, jts, jte, kts, kte   )
	
	
	
	

			 
			 



		

      
          CALL vertical_diffusion ( 'w', w, w_tendf, config_flags,      &
                                    c1f, c2f,                            &
                                    alt, mut, rdn, rdnw, kvdif,          &
                                    ids, ide, jds, jde, kds, kde,        &
                                    ims, ime, jms, jme, kms, kme,        &
                                    its, ite, jts, jte, kts, kte        )
									
									


                   






									
 							
   if(.true.) then   ! template for grav2d 2025 
	

									  
          CALL vertical_diffusion_rho ( 'r', al, rhofm_tendf, config_flags,      &
                                    c1f, c2f,                            &
                                    alt, mut, rdn, rdnw, kvdif,          &
                                    ids, ide, jds, jde, kds, kde,        &
                                    ims, ime, jms, jme, kms, kme,        &
                                    its, ite, jts, jte, kts, kte        )



                  endif 


              if(.true.)  then 
          CALL vertical_diffusion_rho ( 'r', p, pfm_tendf, config_flags,      &
                                     c1f, c2f,                            &
                                    alt, mut, rdn, rdnw, kvdif,          &
                                    ids, ide, jds, jde, kds, kde,        &
                                    ims, ime, jms, jme, kms, kme,        &
                                    its, ite, jts, jte, kts, kte        )

	
             endif 									

kvdq = 3.*kvdif 

  
        !  CALL vertical_diffusion_3dmp ( rhofm, rhofm_tendf, config_flags, rhofm_base,     &
        !                                 c1h, c2h,                             &
        !                                 alt, mut, rdn, rdnw, kvdq ,           &
        !                                 ids, ide, jds, jde, kds, kde,         &
        !                                 ims, ime, jms, jme, kms, kme,         &
        !                                 its, ite, jts, jte, kts, kte         )
		!								 
										 
										 
										 
										 

        ENDIF pbl_test



    END IF diff_opt1

 

   

  END IF forward_step
  
  
 
END SUBROUTINE rk_tendency

!-------------------------------------------------------------------------------

SUBROUTINE rk_addtend_dry ( ru_tend, rv_tend, rw_tend, ph_tend, t_tend,      &
                            ru_tendf, rv_tendf, rw_tendf, ph_tendf, t_tendf, &
                            u_save, v_save, w_save, ph_save, t_save,         &
                            mu_tend, mu_tendf, rk_step, c1, c2,              &
                            h_diabatic, mut, msftx, msfty, msfux, msfuy,     &
                            msfvx, msfvx_inv, msfvy,                         &
                            ids,ide, jds,jde, kds,kde,                       &
                            ims,ime, jms,jme, kms,kme,                       &
                            ips,ipe, jps,jpe, kps,kpe,                       &
                            its,ite, jts,jte, kts,kte                       )

   IMPLICIT NONE

   !  Input data.

   INTEGER ,               INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                            ims, ime, jms, jme, kms, kme, &
                                            ips, ipe, jps, jpe, kps, kpe, &
                                            its, ite, jts, jte, kts, kte
   INTEGER ,               INTENT(IN   ) :: rk_step

   REAL , DIMENSION( ims:ime , kms:kme, jms:jme  ) , INTENT(INOUT) :: ru_tend, &
                                                                      rv_tend, &
                                                                      rw_tend, &
                                                                      ph_tend, &
                                                                      t_tend,  &
                                                                      ru_tendf, &
                                                                      rv_tendf, &
                                                                      rw_tendf, &
                                                                      ph_tendf, &
                                                                      t_tendf

   REAL , DIMENSION( ims:ime , jms:jme  ) , INTENT(INOUT) :: mu_tend, &
                                                             mu_tendf

   REAL , DIMENSION( ims:ime , kms:kme, jms:jme  ) , INTENT(IN   ) ::  u_save,  &
                                                                       v_save,  &
                                                                       w_save,  &
                                                                      ph_save,  &
                                                                       t_save,  &
                                                                      h_diabatic

   REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: mut,       &
                                                                    msftx,     &
                                                                    msfty,     &
                                                                    msfux,     &
                                                                    msfuy,     &
                                                                    msfvx,     &
                                                                    msfvx_inv, &
                                                                    msfvy

   REAL , DIMENSION( kms:kme ) ,         INTENT(IN   ) :: c1, c2

! Local
   INTEGER :: i, j, k

!<DESCRIPTION>
!
! rk_addtend_dry constructs the full large-timestep tendency terms for
! momentum (u,v,w), theta and geopotential equations.   This is accomplished
! by combining the physics tendencies (in *tendf; these are computed
! the first RK substep, held fixed thereafter) with the RK tendencies
! (in *tend, these include advection, pressure gradient, etc;
! these change each rk substep).  Output is in *tend.
!
!</DESCRIPTION>

!  Finally, add the forward-step tendency to the rk_tendency

! u/v/w/save contain bc tendency that needs to be multiplied by msf
! (u by msfuy, v by msfvx)
!  before adding it to physics tendency (*tendf)
! For momentum we need the final tendency to include an inverse msf
! physics/bc tendency needs to be divided, advection tendency already has it

! For scalars we need the final tendency to include an inverse msf (msfty)
! advection tendency is OK, physics/bc tendency needs to be divided by msf

   DO j = jts,MIN(jte,jde-1)
   DO k = kts,kte-1
   DO i = its,ite
     ! multiply by my to uncouple u
     IF(rk_step == 1)ru_tendf(i,k,j) = ru_tendf(i,k,j) +  u_save(i,k,j)*msfuy(i,j)
     ! divide by my to couple u
     ru_tend(i,k,j) = ru_tend(i,k,j) + ru_tendf(i,k,j)/msfuy(i,j)
   ENDDO
   ENDDO
   ENDDO

   DO j = jts,jte
   DO k = kts,kte-1
   DO i = its,MIN(ite,ide-1)
     ! multiply by mx to uncouple v
     IF(rk_step == 1)rv_tendf(i,k,j) = rv_tendf(i,k,j) +  v_save(i,k,j)*msfvx(i,j)
     ! divide by mx to couple v
     rv_tend(i,k,j) = rv_tend(i,k,j) + rv_tendf(i,k,j)*msfvx_inv(i,j)
   ENDDO
   ENDDO
   ENDDO

   DO j = jts,MIN(jte,jde-1)
   DO k = kts,kte
   DO i = its,MIN(ite,ide-1)
     ! multiply by my to uncouple w
     IF(rk_step == 1)rw_tendf(i,k,j) = rw_tendf(i,k,j) +  w_save(i,k,j)*msfty(i,j)
     ! divide by my to couple w
     rw_tend(i,k,j) = rw_tend(i,k,j) + rw_tendf(i,k,j)/msfty(i,j)
     IF(rk_step == 1)ph_tendf(i,k,j) = ph_tendf(i,k,j) +  ph_save(i,k,j)
     ! divide by my to couple scalar
     ph_tend(i,k,j) = ph_tend(i,k,j) + ph_tendf(i,k,j)/msfty(i,j)
   ENDDO
   ENDDO
   ENDDO

   DO j = jts,MIN(jte,jde-1)
   DO k = kts,kte-1
   DO i = its,MIN(ite,ide-1)
     IF(rk_step == 1)t_tendf(i,k,j) = t_tendf(i,k,j) +  t_save(i,k,j)
     ! divide by my to couple theta
      t_tend(i,k,j) =  t_tend(i,k,j) +  t_tendf(i,k,j)/msfty(i,j)  &
                                     +  (c1(k)*mut(i,j)+c2(k))*h_diabatic(i,k,j)/msfty(i,j)
     ! divide by my to couple heating
   ENDDO
   ENDDO
   ENDDO

   DO j = jts,MIN(jte,jde-1)
   DO i = its,MIN(ite,ide-1)
! mu tendencies not coupled with 1/msf
      MU_TEND(i,j) =  MU_TEND(i,j) +  MU_TENDF(i,j)
   ENDDO
   ENDDO

END SUBROUTINE rk_addtend_dry



!-------------------------------------------------------------------------------


SUBROUTINE rk_addtend_dry_fm ( u_tend, v_tend, w_tend, pfm_tend, rhofm_tend,      &
                               u_tendf, v_tendf, w_tendf, pfm_tendf, rhofm_tendf, &
                            u_save, v_save, w_save, pfm_save, rhofm_save,         &
							     rk_step,                      &
                      !      mu_tend, mu_tendf, rk_step, c1, c2,              &
                      !      h_diabatic, mut, msftx, msfty, msfux, msfuy,     &
                      !      msfvx, msfvx_inv, msfvy,                         &
                            ids,ide, jds,jde, kds,kde,                       &
                            ims,ime, jms,jme, kms,kme,                       &
                            ips,ipe, jps,jpe, kps,kpe,                       &
                            its,ite, jts,jte, kts,kte                       )

   IMPLICIT NONE

   !  Input data.

   INTEGER ,               INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                            ims, ime, jms, jme, kms, kme, &
                                            ips, ipe, jps, jpe, kps, kpe, &
                                            its, ite, jts, jte, kts, kte
   INTEGER ,               INTENT(IN   ) :: rk_step

   REAL , DIMENSION( ims:ime , kms:kme, jms:jme  ) , INTENT(INOUT) :: u_tend, &
                                                                      v_tend, &
                                                                      w_tend, &
                                                                      pfm_tend, &
                                                                      rhofm_tend, &
                                                                      u_tendf, &
                                                                      v_tendf, &
                                                                      w_tendf, &
                                                                      pfm_tendf, &
                                                                      rhofm_tendf 



   REAL , DIMENSION( ims:ime , kms:kme, jms:jme  ) , INTENT(IN) ::  u_save, v_save, w_save, pfm_save, rhofm_save
  
! Local
   INTEGER :: i, j, k

!<DESCRIPTION>
!
! rk_addtend_dry constructs the full large-timestep tendency terms for
! momentum (u,v,w), theta and geopotential equations.   This is accomplished
! by combining the physics tendencies (in *tendf; these are computed
! the first RK substep, held fixed thereafter) with the RK tendencies
! (in *tend, these include advection, pressure gradient, etc;
! these change each rk substep).  Output is in *tend.
!
!</DESCRIPTION>

!  Finally, add the forward-step tendency to the rk_tendency

! u/v/w/save contain bc tendency that needs to be multiplied by msf
! (u by msfuy, v by msfvx)
!  before adding it to physics tendency (*tendf)
! For momentum we need the final tendency to include an inverse msf
! physics/bc tendency needs to be divided, advection tendency already has it

! For scalars we need the final tendency to include an inverse msf (msfty)
! advection tendency is OK, physics/bc tendency needs to be divided by msf




   DO j = jts,MIN(jte,jde-1)
   DO k = kts,kte-1
   DO i = its,ite

  !   IF(rk_step == 1)ru_tendf(i,k,j) = ru_tendf(i,k,j) +  u_save(i,k,j)*msfuy(i,j)
         IF(rk_step == 1) u_tendf(i,k,j) =  u_tendf(i,k,j) +  u_save(i,k,j)  
     u_tend(i,k,j) = u_tend(i,k,j) + u_tendf(i,k,j)  
   ENDDO
   ENDDO
   ENDDO

   DO j = jts,jte
   DO k = kts,kte-1
   DO i = its,MIN(ite,ide-1)
     ! multiply by mx to uncouple v
   !  IF(rk_step == 1)rv_tendf(i,k,j) = rv_tendf(i,k,j) +  v_save(i,k,j)*msfvx(i,j)
      IF(rk_step == 1)    v_tendf(i,k,j) = v_tendf(i,k,j) +  v_save(i,k,j) 
     v_tend(i,k,j) = v_tend(i,k,j) + v_tendf(i,k,j)  
   ENDDO
   ENDDO
   ENDDO

! w
   DO j = jts,MIN(jte,jde-1)
   DO k = kts,kte
   DO i = its,MIN(ite,ide-1)
   !  ph_tend(i,k,j) = ph_tend(i,k,j) + ph_tendf(i,k,j) 
      IF(rk_step == 1)    w_tendf(i,k,j) = w_tendf(i,k,j) +  w_save(i,k,j)
	   w_tend(i,k,j) = w_tend(i,k,j) + w_tendf(i,k,j)  
   ENDDO
   ENDDO
   ENDDO


! pfm 
   DO j = jts,MIN(jte,jde-1)
   DO k = kts,kte-1
   DO i = its,MIN(ite,ide-1)
    ! IF(rk_step == 1)t_tendf(i,k,j) = t_tendf(i,k,j) +  t_save(i,k,j)
     ! divide by my to couple theta
	 IF(rk_step == 1)    pfm_tendf(i,k,j) = pfm_tendf(i,k,j) +  pfm_save(i,k,j) 
	 
      pfm_tend(i,k,j) =  pfm_tend(i,k,j) +  pfm_tendf(i,k,j)  

   ENDDO
   ENDDO
   ENDDO
! rhofm 

   DO j = jts,MIN(jte,jde-1)
   DO k = kts,kte-1
   DO i = its,MIN(ite,ide-1)
    ! IF(rk_step == 1)t_tendf(i,k,j) = t_tendf(i,k,j) +  t_save(i,k,j)
     ! divide by my to couple theta
	IF(rk_step == 1)  rhofm_tendf(i,k,j) =   rhofm_tendf(i,k,j) +  rhofm_save(i,k,j) 
	 
   rhofm_tend(i,k,j) =  rhofm_tend(i,k,j) +  rhofm_tendf(i,k,j)   

   ENDDO
   ENDDO
   ENDDO

END SUBROUTINE rk_addtend_dry_fm



!-------------------------------------------------------------------------------

SUBROUTINE rk_scalar_tend_fm ( scs, sce, config_flags,          &
                            tenddec,                        &
                            rk_step, dt,                     &
							rho,                              &       !is grid%rhofm_1
						    msft,zs,zh,fzm,fzp,                    &  ! hanshenyou
                                ztop,                  &     ! hanshenyou     
                            u, v, w,                        &
							mut, mub, mu_old,    &
                            c1h, c2h, alt,                   &
                            scalar_old, scalar,              &
                            scalar_tends, advect_tend,       &
                            h_tendency, z_tendency,          &
                            RQVFTEN,                         &
                            base, moist_step, fnm, fnp,      &
                            msfux, msfuy, msfvx, msfvx_inv,  &
                            msfvy, msftx, msfty,             &
                            rdx, rdy, rdn, rdnw,             &
                            khdif, kvdif, xkmhd,             &
                            diff_6th_opt, diff_6th_factor,   &
                            adv_opt,                         &
                            phb, ph,                         &
                            mix2_off, mix6_off,              &
                            ids, ide, jds, jde, kds, kde,    &
                            ims, ime, jms, jme, kms, kme,    &
                            its, ite, jts, jte, kts, kte    )

   IMPLICIT NONE

   !  Input data.

   TYPE(grid_config_rec_type   ) ,   INTENT(IN   ) :: config_flags

   LOGICAL ,                INTENT(IN   ) :: tenddec ! tendency term

   INTEGER ,                INTENT(IN   ) :: rk_step, scs, sce
   INTEGER ,                INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                             ims, ime, jms, jme, kms, kme, &
                                             its, ite, jts, jte, kts, kte

   LOGICAL , INTENT(IN   ) :: moist_step

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme , scs:sce ),                &
                                         INTENT(IN   )  :: scalar, scalar_old

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme , scs:sce ),                      &
                                         INTENT(INOUT)  :: scalar_tends
                                                    
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme  ), INTENT(INOUT) :: advect_tend
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme  ), INTENT(  OUT) :: h_tendency, z_tendency

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme  ), INTENT(OUT  ) :: RQVFTEN

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme  ), INTENT(IN   ) ::     u,  &
                                                                      v,  &
                                                                      w,  &
																	  rho,  & 
                                                                      xkmhd,  &
                                                                      alt,    &
                                                                      phb,    &
                                                                      ph

   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: fnm,  &
                                                                  fnp,  &
                                                                  rdn,  &
                                                                  rdnw, &
                                                                  base

   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: c1h, c2h

   REAL , DIMENSION( ims:ime , jms:jme ) ,       INTENT(IN   ) :: msfux,    &
                                                                  msfuy,    &
                                                                  msfvx,    &
                                                                  msfvx_inv,    &
                                                                  msfvy,    &
                                                                  msftx,    &
                                                                  msfty,    &
                                                                  mub,     &
                                                                  mut,     &
                                                                  mu_old

   REAL ,                                        INTENT(IN   ) :: rdx,     &
                                                                  rdy,     &
                                                                  khdif,   &
                                                                  kvdif

   INTEGER, INTENT( IN ) :: diff_6th_opt
   REAL,    INTENT( IN ) :: diff_6th_factor

   REAL ,                                        INTENT(IN   ) :: dt

   INTEGER, INTENT(IN   ) :: adv_opt          
   LOGICAL, INTENT(IN   ) :: mix2_off, mix6_off

   ! Local data
  
   INTEGER :: im, i,j,k
   INTEGER :: time_step

   REAL    :: khdq, kvdq, tendency
   
   
   !2023
   



         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,zs
         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh,fzm,fzp
         REAL ,                                        INTENT(IN   ) :: ztop



   REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: sigmdot,rhoscalar 

  INTEGER :: i_start, i_end, j_start, j_end, k_start, k_end
  
  

!<DESCRIPTION>
!
! rk_scalar_tend calls routines that computes scalar tendency from advection
! and 3D mixing (TKE or fixed eddy viscosities).
!
!</DESCRIPTION>


   khdq = khdif/prandtl
   kvdq = kvdif/prandtl
   
   
   
       call    calculate_sgmdot( sigmdot,                 &
                        u, v, w,                       &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  ) 
						
						
						

   scalar_loop : DO im = scs, sce

     CALL zero_tend ( advect_tend(ims,kms,jms),     &
                      ids, ide, jds, jde, kds, kde, &
                      ims, ime, jms, jme, kms, kme, &
                      its, ite, jts, jte, kts, kte )

     CALL zero_tend ( h_tendency(ims,kms,jms),      &
                      ids, ide, jds, jde, kds, kde, &
                      ims, ime, jms, jme, kms, kme, &
                      its, ite, jts, jte, kts, kte )

     CALL zero_tend ( z_tendency(ims,kms,jms),      &
                      ids, ide, jds, jde, kds, kde, &
                      ims, ime, jms, jme, kms, kme, &
                      its, ite, jts, jte, kts, kte )

     CALL nl_get_time_step ( 1, time_step )

  
   rhoscalar = 0.0 
   
     !  for mpi 扩1   

    k_start = kts
    k_end   = MIN(kte,kde-1)	 
    i_start = MAX(its-1,ids)
    i_end   = MIN(ite+1,ide-1)
    j_start = MAX(jts-1,jds)
    j_end   = MIN(jte+1,jde-1) 
	
	
   
      DO      j = j_start, j_end
      DO  k = k_start, k_end
      DO  i = i_start, i_end
   
		rhoscalar(i,k,j) =  rho(i,k,j)*scalar(i,k,j,im)

   ENDDO
   ENDDO
   ENDDO
   
   ! 20230603 
          CALL set_physical_bc3d( rhoscalar, 'p', config_flags,               &
                               ids, ide, jds, jde, kds, kde,       &
                               ims, ime, jms, jme, kms, kme,       &
                               its, ite, jts, jte, kts, kte,       &
                               its, ite, jts, jte, kts, kte )
	 
	 
   
     IF( (rk_step == 3) .and. (adv_opt == POSITIVEDEF) ) THEN  !202304 ! 20230524
    !   IF( (rk_step == 3)  ) THEN    ! 2023 cut hanshenyou 
   
   
	 
		 
        CALL advect_scalar_pd_fm       ( rhoscalar ,             &
                                      rhoscalar ,         &
                                      advect_tend(ims,kms,jms),           &
                                      h_tendency(ims,kms,jms),            &
                                      z_tendency(ims,kms,jms),            &
                                      u, v, sigmdot, c1h, c2h,               &
                                      mut, mub, mu_old,                   &
                                      time_step, config_flags, tenddec,   &
                                      msfux, msfuy, msfvx, msfvy,         &
                                      msftx, msfty,               &
									  msft,   &
									  fnm, fnp,ztop,zs,             &
                                      rdx, rdy, rdnw,dt,                  &
                                      ids, ide, jds, jde, kds, kde,       &
                                      ims, ime, jms, jme, kms, kme,       &
                                      its, ite, jts, jte, kts, kte     )
									  
									  
			




             ELSE IF( (rk_step == 3) .and. (adv_opt == MONOTONIC) ) THEN ! 20230824 
   

                    CALL advect_scalar_mono_fm( rhoscalar ,             &
                                      rhoscalar ,         &
                                      advect_tend(ims,kms,jms),           &
                                      h_tendency(ims,kms,jms),            &
                                      z_tendency(ims,kms,jms),            &
                                      u, v, sigmdot, c1h, c2h,               &
                                      mut, mub, mu_old,                   &
                                      time_step, config_flags, tenddec,   &
                                      msfux, msfuy, msfvx, msfvy,         &
                                      msftx, msfty,               &
									  msft,   &
									  fnm, fnp,ztop,zs,             &
                                      rdx, rdy, rdnw,dt,                  &
                                      ids, ide, jds, jde, kds, kde,       &
                                      ims, ime, jms, jme, kms, kme,       &
                                      its, ite, jts, jte, kts, kte     )
                      

			

 
	 ELSE
	 
	  if(.true.) then 
		     call calcu_divg( rhoscalar,        &
    !  call advect_scalar_fm ( scalar(ims,kms,jms,im),        &
                                 advect_tend(ims,kms,jms),      &
                        u, v, sigmdot,                       &
                     config_flags,                  &
                        msft,zh,zs,fzm,fzp,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )		 
	 

         else
		 
        CALL advect_scalar_fm2     ( rhoscalar ,             &
                                      rhoscalar ,         &
                                      advect_tend(ims,kms,jms),           &
                                    !  h_tendency(ims,kms,jms),            &
                                    !  z_tendency(ims,kms,jms),            &
                                      u, v, sigmdot,      &
                                      time_step, config_flags,    &
                                      msft,   &
									  fnm, fnp,ztop,zs,             &
                                      rdx, rdy, rdnw,dt,                  &
                                      ids, ide, jds, jde, kds, kde,       &
                                      ims, ime, jms, jme, kms, kme,       &
                                      its, ite, jts, jte, kts, kte     )		 
		 
		 
		 endif 
	 


	 
								 

      END IF

     IF((config_flags%cu_physics == GDSCHEME .OR. config_flags%cu_physics == G3SCHEME .OR. &
         config_flags%cu_physics == GFSCHEME .OR.    &
         config_flags%cu_physics == KFETASCHEME .OR. config_flags%cu_physics == MSKFSCHEME .OR. &
         config_flags%cu_physics == TIEDTKESCHEME .OR. config_flags%cu_physics == NTIEDTKESCHEME  ) &      ! Tiedtke
                     .and. moist_step .and. ( im == P_QV) ) THEN

        CALL set_tend( RQVFTEN, advect_tend, msfty,    &
                       ids, ide, jds, jde, kds, kde,   &
                       ims, ime, jms, jme, kms, kme,   &
                       its, ite, jts, jte, kts, kte      )
     ENDIF

     rk_step_1: IF( rk_step == 1 ) THEN

       diff_opt1 : IF (config_flags%diff_opt .eq. 1) THEN

       IF (.not. mix2_off)                                                   &
   !      CALL horizontal_diffusion ( 'm', scalar(ims,kms,jms,im),            &
         CALL horizontal_diffusion ( 'm', rhoscalar(ims,kms,jms ),            &	    ! test 20230603 	 需要获得耦合密度的倾向
                                          scalar_tends(ims,kms,jms,im), mut, &
                                          c1h, c2h, config_flags,            &
                                          msfux, msfuy, msfvx, msfvx_inv,    &
                                          msfvy, msftx, msfty,               &
                                          khdq , xkmhd, rdx, rdy,            &
                                          ids, ide, jds, jde, kds, kde,      &
                                          ims, ime, jms, jme, kms, kme,      &
                                          its, ite, jts, jte, kts, kte      )

       pbl_test : IF (config_flags%bl_pbl_physics .eq. 0) THEN

         IF( (moist_step) .and. ( im == P_QV)) THEN
		 
                ! vertical_diffusion_mp 需要减去基本态，所以用q，不耦合密度。
				! 倾向需要耦合密度。alt = 1.0/rhofm  
				! scalar_tends 是 耦合密度的比湿
				
            CALL vertical_diffusion_mp ( scalar(ims,kms,jms,im),       &
                                         scalar_tends(ims,kms,jms,im), &
                                         config_flags, base, c1h, c2h, &
                                         alt, mut, rdn, rdnw, kvdq ,   &
                                         ids, ide, jds, jde, kds, kde, &
                                         ims, ime, jms, jme, kms, kme, &
                                         its, ite, jts, jte, kts, kte )


         ELSE

            CALL vertical_diffusion (  'm', scalar(ims,kms,jms,im),       &
                                            scalar_tends(ims,kms,jms,im), &
                                            config_flags, c1h, c2h,       &
                                            alt, mut, rdn, rdnw, kvdq,    &
                                            ids, ide, jds, jde, kds, kde, &
                                            ims, ime, jms, jme, kms, kme, &
                                            its, ite, jts, jte, kts, kte )

         END IF

      ENDIF pbl_test
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  

    ENDIF diff_opt1
	
	if(.false.) then  ! for wrf

    IF ( (diff_6th_opt .NE. 0 ) .and. (.not. mix6_off) )              &
      CALL sixth_order_diffusion( 'm', scalar(ims,kms,jms,im),        &
                                       scalar_tends(ims,kms,jms,im),  &
                                       mut, dt, config_flags, c1h,c2h,&
                                       diff_6th_opt, diff_6th_factor, &
                                       phb, ph,                       &
                                       rdx, rdy,                      &
                                       msftx, msfty,                  &
                                       msfux, msfuy,                  &
                                       msfvx, msfvy,                  &
                                       ids, ide, jds, jde, kds, kde,  &
                                       ims, ime, jms, jme, kms, kme,  &
                                       its, ite, jts, jte, kts, kte )

  
   endif 
   
  ENDIF rk_step_1

 END DO scalar_loop

END SUBROUTINE rk_scalar_tend_fm

!-------------------------------------------------------------------------------

SUBROUTINE q_diabatic_add_fm ( scs, sce,                        &
                            dt,  rho,   &
						! 	mut, c1, c2,                 &
                            qv_diabatic, qc_diabatic,        &
                            scalar_tends,                    &
                            ids, ide, jds, jde, kds, kde,    &
                            ims, ime, jms, jme, kms, kme,    &
                            its, ite, jts, jte, kts, kte    )

   IMPLICIT NONE

   !  Input data.

   INTEGER ,                INTENT(IN   ) :: scs, sce
   INTEGER ,                INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                             ims, ime, jms, jme, kms, kme, &
                                             its, ite, jts, jte, kts, kte



   REAL, DIMENSION(ims:ime, kms:kme, jms:jme ),                &
                                         INTENT(IN   )  :: rho 
										 
										 
!   REAL,     DIMENSION( ims:ime, jms:jme )                        , &
!             INTENT(IN   )   ::                                 mut

 !  REAL , DIMENSION( kms:kme ) , INTENT(IN   ) :: c1, c2

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme ),                &
                                         INTENT(IN   )  :: qv_diabatic, qc_diabatic

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme , scs:sce ),                &
                                         INTENT(INOUT)  :: scalar_tends

   REAL ,                                        INTENT(IN   ) :: dt

   ! Local data
  
   INTEGER :: im, i,j,k

!<DESCRIPTION>
!
!</DESCRIPTION>

!!!       print *,'  q_diab add: '
!!!       print *,its,MIN(ite,ide-1)
!!!       print *,jts,MIN(jte,jde-1)
!!!       print *,kts,kte-1

   scalar_loop : DO im = scs, sce

     IF( im.eq.p_qv )THEN
!!!       print *,'  qv '
       DO j = jts,MIN(jte,jde-1)
       DO k = kts,kte-1
       DO i = its,MIN(ite,ide-1)
         scalar_tends(i,k,j,im) = scalar_tends(i,k,j,im) + qv_diabatic(i,k,j)*rho(i,k,j)
       ENDDO
       ENDDO
       ENDDO
     ENDIF
     IF( im.eq.p_qc )THEN
!!!       print *,'  qc '
       DO j = jts,MIN(jte,jde-1)
       DO k = kts,kte-1
       DO i = its,MIN(ite,ide-1)
         scalar_tends(i,k,j,im) = scalar_tends(i,k,j,im) + qc_diabatic(i,k,j)*rho(i,k,j) 
       ENDDO
       ENDDO
       ENDDO
     ENDIF

   END DO scalar_loop

END SUBROUTINE q_diabatic_add_fm

!-------------------------------------------------------------------------------

SUBROUTINE q_diabatic_subtr( scs, sce,                       &
                            dt,                     &
                            qv_diabatic, qc_diabatic,        &
                            scalar,              &
                            ids, ide, jds, jde, kds, kde,    &
                            ims, ime, jms, jme, kms, kme,    &
                            its, ite, jts, jte, kts, kte    )

   IMPLICIT NONE

   !  Input data.

   INTEGER ,                INTENT(IN   ) :: scs, sce
   INTEGER ,                INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                             ims, ime, jms, jme, kms, kme, &
                                             its, ite, jts, jte, kts, kte

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme ),                &
                                         INTENT(IN   )  :: qv_diabatic, qc_diabatic

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme , scs:sce ),                &
                                         INTENT(INOUT)  :: scalar

   REAL ,                                        INTENT(IN   ) :: dt

   ! Local data
  
   INTEGER :: im, i,j,k

!<DESCRIPTION>
!
!</DESCRIPTION>


!!!       print *,'  q_diab subtr, dt = : ',dt
!!!       print *,its,MIN(ite,ide-1)
!!!       print *,jts,MIN(jte,jde-1)
!!!       print *,kts,kte-1

   scalar_loop : DO im = scs, sce

     IF( im.eq.p_qv )THEN
!!!       print *,'  qv '
       DO j = jts,MIN(jte,jde-1)
       DO k = kts,kte-1
       DO i = its,MIN(ite,ide-1)
         scalar(i,k,j,im) = scalar(i,k,j,im) - dt*qv_diabatic(i,k,j)
       ENDDO
       ENDDO
       ENDDO
     ENDIF
     IF( im.eq.p_qc )THEN
!!!       print *,'  qc '
       DO j = jts,MIN(jte,jde-1)
       DO k = kts,kte-1
       DO i = its,MIN(ite,ide-1)
         scalar(i,k,j,im) = scalar(i,k,j,im) - dt*qc_diabatic(i,k,j)
       ENDDO
       ENDDO
       ENDDO
     ENDIF

   END DO scalar_loop

END SUBROUTINE q_diabatic_subtr


!-------------------------------------------------------------------------------

SUBROUTINE rk_update_scalar( scs, sce,                      &
                             scalar_1, scalar_2, sc_tend,   &
                             advh_t, advz_t,                &
                             advect_tend,                   &
                             h_tendency, z_tendency,        &
                          !   msftx, msfty, c1, c2,          &
                          !   mu_old, mu_new, mu_base,       &
						  rho_old,rho_new,    &
                             rk_step, dt, spec_zone,        &
                             config_flags,                  &
                             tenddec,                      &
                             ids, ide, jds, jde, kds, kde,  &
                             ims, ime, jms, jme, kms, kme,  &
                             its, ite, jts, jte, kts, kte  )

   IMPLICIT NONE

   !  Input data.

   TYPE(grid_config_rec_type   ) ,   INTENT(IN   ) :: config_flags

   LOGICAL ,                INTENT(IN   ) :: tenddec

   INTEGER ,                INTENT(IN   ) :: scs, sce, rk_step, spec_zone
   INTEGER ,                INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                             ims, ime, jms, jme, kms, kme, &
                                             its, ite, jts, jte, kts, kte

   REAL,                    INTENT(IN   ) :: dt

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme , scs:sce),                &
         INTENT(INOUT)                                  :: scalar_1,    &
                                                           scalar_2

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme , scs:sce),                &
         INTENT(IN)                                     :: sc_tend

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme ),                &
         INTENT(IN)                                  :: advect_tend

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme ), INTENT(INOUT) , OPTIONAL :: advh_t,  advz_t ! accumulating for output
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: h_tendency, z_tendency ! from rk_scalar_tend


   REAL, DIMENSION(ims:ime, kms:kme, jms:jme ),                &
         INTENT(IN)                                     :: rho_old,rho_new
		 
 !  REAL, DIMENSION(ims:ime, jms:jme  ), INTENT(IN   ) ::  mu_old,  &
 !                                                         mu_new,  &
 !                                                         mu_base, &
 !                                                         msftx,   &
 !                                                         msfty

 !  REAL, DIMENSION(kms:kme ), INTENT(IN   ) ::  c1, c2

   INTEGER :: i,j,k,im
   REAL    :: sc_middle, msfsq
   REAL, DIMENSION(its:ite) :: muold, munew

   REAL, DIMENSION(its:ite, kts:kte, jts:jte  ) :: tendency

   INTEGER :: i_start,i_end,j_start,j_end,k_start,k_end
   INTEGER :: i_start_spc,i_end_spc,j_start_spc,j_end_spc,k_start_spc,k_end_spc

!<DESCRIPTION>
!
!  rk_scalar_update advances the scalar equation given the time t value
!  of the scalar and the scalar tendency.  
!
!</DESCRIPTION>


!
!  set loop limits.

      i_start = its
      i_end   = min(ite,ide-1)
      j_start = jts
      j_end   = min(jte,jde-1)
      k_start = kts
      k_end   = kte-1

      i_start_spc = i_start
      i_end_spc   = i_end
      j_start_spc = j_start
      j_end_spc   = j_end
      k_start_spc = k_start
      k_end_spc   = k_end

    IF( config_flags%nested .or. config_flags%specified ) THEN
      IF( .NOT. config_flags%periodic_x)i_start = max( its,ids+spec_zone )
      IF( .NOT. config_flags%periodic_x)i_end   = min( ite,ide-spec_zone-1 )
      j_start = max( jts,jds+spec_zone )
      j_end   = min( jte,jde-spec_zone-1 )
      k_start = kts
      k_end   = min( kte, kde-1 )
    ENDIF

    IF ( rk_step == 1 ) THEN

      !  replace t-dt values (in scalar_1) with t values scalar_2,
      !  then compute new values by adding tendency to values at t

      DO  im = scs,sce

       DO  j = jts, min(jte,jde-1)
       DO  k = kts, min(kte,kde-1)
       DO  i = its, min(ite,ide-1)
           tendency(i,k,j) = 0.
       ENDDO
       ENDDO
       ENDDO
   
    !   DO  j = j_start,j_end
    !   DO  k = k_start,k_end
    !   DO  i = i_start,i_end
    !      ! scalar was coupled with my
     !      tendency(i,k,j) = advect_tend(i,k,j) * msfty(i,j)
    !   ENDDO
    !   ENDDO
    !   ENDDO
   
       DO  j = j_start_spc,j_end_spc
       DO  k = k_start_spc,k_end_spc
       DO  i = i_start_spc,i_end_spc
           tendency(i,k,j) = advect_tend(i,k,j) + sc_tend(i,k,j,im)
       ENDDO
       ENDDO
       ENDDO
   
      DO  j = jts, min(jte,jde-1)

    !  DO  i = its, min(ite,ide-1)
    !    MUOLD(i) = MU_OLD(i,j) + MU_BASE(i,j)
   !     MUNEW(i) = MU_NEW(i,j) + MU_BASE(i,j)
   !   ENDDO

      DO  k = kts, min(kte,kde-1)
      DO  i = its, min(ite,ide-1)

        scalar_1(i,k,j,im) = scalar_2(i,k,j,im)
    !    scalar_2(i,k,j,im) = ((c1(k)*muold(i)+c2(k))*scalar_1(i,k,j,im)   &
    !                         + dt*tendency(i,k,j))/(c1(k)*munew(i)+c2(k))
							 
     !   scalar_2(i,k,j,im) =  scalar_1(i,k,j,im)   + dt*tendency(i,k,j)  					 

           scalar_2(i,k,j,im) =  (rho_old(i,k,j)*scalar_1(i,k,j,im)   &
                             + dt*tendency(i,k,j))/rho_new(i,k,j)	
							 
 
			
!             20230806 正定限制
		    scalar_2(i,k,j,im)  =  max(0.0, scalar_2(i,k,j,im)    ) 							 

      ENDDO
      ENDDO
      ENDDO

      ENDDO

    ELSE
 
      !  just compute new values, scalar_1 already at time t.

      DO  im = scs, sce

       DO  j = jts, min(jte,jde-1)
       DO  k = kts, min(kte,kde-1)
       DO  i = its, min(ite,ide-1)
           tendency(i,k,j) = 0.
       ENDDO
       ENDDO
       ENDDO
   
   !    DO  j = j_start,j_end
    !   DO  k = k_start,k_end
    !   DO  i = i_start,i_end
    !       ! scalar was coupled with my
    !       tendency(i,k,j) = advect_tend(i,k,j) * msfty(i,j)
    !   ENDDO
    !   ENDDO
   !    ENDDO
   
       DO  j = j_start_spc,j_end_spc
       DO  k = k_start_spc,k_end_spc
       DO  i = i_start_spc,i_end_spc
         !  tendency(i,k,j) = tendency(i,k,j) + sc_tend(i,k,j,im)
		   tendency(i,k,j) = advect_tend(i,k,j) + sc_tend(i,k,j,im)
       ENDDO
       ENDDO
       ENDDO
	   
	   
	 

      DO  j = jts, min(jte,jde-1)

    !  DO  i = its, min(ite,ide-1)
    !    MUOLD(i) = MU_OLD(i,j) + MU_BASE(i,j)
    !    MUNEW(i) = MU_NEW(i,j) + MU_BASE(i,j)
     ! ENDDO

      DO  k = kts, min(kte,kde-1)
      DO  i = its, min(ite,ide-1)

      !  scalar_2(i,k,j,im) = ((c1(k)*muold(i)+c2(k))*scalar_1(i,k,j,im)   &
     !                        + dt*tendency(i,k,j))/(c1(k)*munew(i)+c2(k))
							 
		!  scalar_2(i,k,j,im) =  	scalar_1(i,k,j,im) + dt*tendency(i,k,j)	 
 
		
           scalar_2(i,k,j,im) =  (rho_old(i,k,j)*scalar_1(i,k,j,im)   &
                             + dt*tendency(i,k,j))/rho_new(i,k,j)	
							 

!             20230806 正定限制
		    scalar_2(i,k,j,im)  =  max(0.0, scalar_2(i,k,j,im)    ) 	


		   
		   

      ENDDO
      ENDDO


     !  call  positive_definite_filter ( a,                          &
      !                                ids,ide, jds,jde, kds,kde,  &
      !                                ims,ime, jms,jme, kms,kme,  &
      !                                its,ite, jts,jte, kts,kte  )
	  

      ! This is separated from the k/i-loop above for better performance
      IF ( PRESENT(advh_t) .AND. PRESENT(advz_t) ) THEN
        IF(tenddec.and.rk_step.eq.config_flags%rk_ord) THEN
          DO  k = kts, min(kte,kde-1)
          DO  i = its, min(ite,ide-1)

         !   advh_t(i,k,j) = advh_t(i,k,j) + (dt*h_tendency(i,k,j)* msfty(i,j))/(c1(k)*munew(i)+c2(k))
         !   advz_t(i,k,j) = advz_t(i,k,j) + (dt*z_tendency(i,k,j)* msfty(i,j))/(c1(k)*munew(i)+c2(k))
			
			
	 !           advh_t(i,k,j) = advh_t(i,k,j) +  dt*h_tendency(i,k,j)  
     !       advz_t(i,k,j) = advz_t(i,k,j) +  dt*z_tendency(i,k,j) 


          advh_t(i,k,j) =  (rho_old(i,k,j)*advh_t(i,k,j)   &
                             + dt*h_tendency(i,k,j) )/rho_new(i,k,j)

           advz_t(i,k,j) =  (rho_old(i,k,j)*advz_t(i,k,j)  &
                             + dt*z_tendency(i,k,j) )/rho_new(i,k,j)							 

          ENDDO
          ENDDO
        END IF
      END IF
	  
	  
	  
      
      ENDDO
	  
	  
  
	  

      ENDDO

    END IF
	

		 

END SUBROUTINE rk_update_scalar

!-------------------------------------------------------------------------------

SUBROUTINE rk_update_scalar_pd( scs, sce,                      &
                                scalar, sc_tend,               &
                                c1, c2,                        &
                                mu_old, mu_new, mu_base,       &
                                rk_step, dt, spec_zone,        &
                                config_flags,                  &
                                ids, ide, jds, jde, kds, kde,  &
                                ims, ime, jms, jme, kms, kme,  &
                                its, ite, jts, jte, kts, kte  )

   IMPLICIT NONE

   !  Input data.

   TYPE(grid_config_rec_type   ) ,   INTENT(IN   ) :: config_flags

   INTEGER ,                INTENT(IN   ) :: scs, sce, rk_step, spec_zone
   INTEGER ,                INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                             ims, ime, jms, jme, kms, kme, &
                                             its, ite, jts, jte, kts, kte

   REAL,                    INTENT(IN   ) :: dt

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme , scs:sce),                &
         INTENT(INOUT)                                  :: scalar,      &
                                                           sc_tend

   REAL, DIMENSION(ims:ime, jms:jme  ), INTENT(IN   ) ::  mu_old,  &
                                                          mu_new,  &
                                                          mu_base

   REAL, DIMENSION(kms:kme ), INTENT(IN   ) ::  c1, c2

   INTEGER :: i,j,k,im
   REAL    :: sc_middle, msfsq
   REAL, DIMENSION(its:ite) :: muold, munew

   REAL, DIMENSION(its:ite, kts:kte, jts:jte  ) :: tendency

   INTEGER :: i_start,i_end,j_start,j_end,k_start,k_end
   INTEGER :: i_start_spc,i_end_spc,j_start_spc,j_end_spc,k_start_spc,k_end_spc

!<DESCRIPTION>
!
!  rk_scalar_update advances the scalar equation given the time t value
!  of the scalar and the scalar tendency.  
!
!</DESCRIPTION>


!
!  set loop limits.

      i_start = its
      i_end   = min(ite,ide-1)
      j_start = jts
      j_end   = min(jte,jde-1)
      k_start = kts
      k_end   = kte-1

      i_start_spc = i_start
      i_end_spc   = i_end
      j_start_spc = j_start
      j_end_spc   = j_end
      k_start_spc = k_start
      k_end_spc   = k_end

    IF( config_flags%nested .or. config_flags%specified ) THEN
      IF( .NOT. config_flags%periodic_x)i_start = max( its,ids+spec_zone )
      IF( .NOT. config_flags%periodic_x)i_end   = min( ite,ide-spec_zone-1 )
      j_start = max( jts,jds+spec_zone )
      j_end   = min( jte,jde-spec_zone-1 )
      k_start = kts
      k_end   = min( kte, kde-1 )
    ENDIF

      DO  im = scs, sce

       DO  j = jts, min(jte,jde-1)
       DO  k = kts, min(kte,kde-1)
       DO  i = its, min(ite,ide-1)
           tendency(i,k,j) = 0.
       ENDDO
       ENDDO
       ENDDO
   
       DO  j = j_start_spc,j_end_spc
       DO  k = k_start_spc,k_end_spc
       DO  i = i_start_spc,i_end_spc
           tendency(i,k,j) = tendency(i,k,j) + sc_tend(i,k,j,im)
           sc_tend(i,k,j,im) = 0.
       ENDDO
       ENDDO
       ENDDO

      DO  j = jts, min(jte,jde-1)

      DO  i = its, min(ite,ide-1)
        MUOLD(i) = MU_OLD(i,j) + MU_BASE(i,j)
        MUNEW(i) = MU_NEW(i,j) + MU_BASE(i,j)
      ENDDO

      DO  k = kts, min(kte,kde-1)
      DO  i = its, min(ite,ide-1)

        scalar(i,k,j,im) = ((c1(k)*muold(i)+c2(k))*scalar(i,k,j,im)   &
                             + dt*tendency(i,k,j))/(c1(k)*munew(i)+c2(k))
      ENDDO
      ENDDO
      ENDDO

      ENDDO

END SUBROUTINE rk_update_scalar_pd

!------------------------------------------------------------
!-------------------------------------------------------------------------------

SUBROUTINE rk_update_scalar_pd_fm( scs, sce,                      &
                                scalar, sc_tend,               &
								rho_old,rho_new,             & 
                 !               c1, c2,                        &
                 !              mu_old, mu_new, mu_base,       &
                                rk_step, dt, spec_zone,        &
                                config_flags,                  &
                                ids, ide, jds, jde, kds, kde,  &
                                ims, ime, jms, jme, kms, kme,  &
                                its, ite, jts, jte, kts, kte  )

   IMPLICIT NONE

   !  Input data.

   TYPE(grid_config_rec_type   ) ,   INTENT(IN   ) :: config_flags

   INTEGER ,                INTENT(IN   ) :: scs, sce, rk_step, spec_zone
   INTEGER ,                INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                             ims, ime, jms, jme, kms, kme, &
                                             its, ite, jts, jte, kts, kte

   REAL,                    INTENT(IN   ) :: dt

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme , scs:sce),                &
         INTENT(INOUT)                                  :: scalar,      &
                                                           sc_tend


REAL, DIMENSION(ims:ime, kms:kme, jms:jme ),                &
         INTENT(IN)                                  ::  rho_old,rho_new 

   INTEGER :: i,j,k,im
   REAL    :: sc_middle, msfsq
   REAL, DIMENSION(its:ite) :: muold, munew

   REAL, DIMENSION(its:ite, kts:kte, jts:jte  ) :: tendency

   INTEGER :: i_start,i_end,j_start,j_end,k_start,k_end
   INTEGER :: i_start_spc,i_end_spc,j_start_spc,j_end_spc,k_start_spc,k_end_spc

!<DESCRIPTION>
!
!  rk_scalar_update advances the scalar equation given the time t value
!  of the scalar and the scalar tendency.  
!
!</DESCRIPTION>


!
!  set loop limits.

      i_start = its
      i_end   = min(ite,ide-1)
      j_start = jts
      j_end   = min(jte,jde-1)
      k_start = kts
      k_end   = kte-1

      i_start_spc = i_start
      i_end_spc   = i_end
      j_start_spc = j_start
      j_end_spc   = j_end
      k_start_spc = k_start
      k_end_spc   = k_end

    IF( config_flags%nested .or. config_flags%specified ) THEN
      IF( .NOT. config_flags%periodic_x)i_start = max( its,ids+spec_zone )
      IF( .NOT. config_flags%periodic_x)i_end   = min( ite,ide-spec_zone-1 )
      j_start = max( jts,jds+spec_zone )
      j_end   = min( jte,jde-spec_zone-1 )
      k_start = kts
      k_end   = min( kte, kde-1 )
    ENDIF
	
	
	
	     
	  

      DO  im = scs, sce

       DO  j = jts, min(jte,jde-1)
       DO  k = kts, min(kte,kde-1)
       DO  i = its, min(ite,ide-1)
           tendency(i,k,j) = 0.
       ENDDO
       ENDDO
       ENDDO
   
       DO  j = j_start_spc,j_end_spc
       DO  k = k_start_spc,k_end_spc
       DO  i = i_start_spc,i_end_spc
           tendency(i,k,j) = tendency(i,k,j) + sc_tend(i,k,j,im)
           sc_tend(i,k,j,im) = 0.
       ENDDO
       ENDDO
       ENDDO

      DO  j = jts, min(jte,jde-1)

 

      DO  k = kts, min(kte,kde-1)
      DO  i = its, min(ite,ide-1)

        scalar(i,k,j,im) =  (rho_old(i,k,j)*scalar(i,k,j,im)   &
                             + dt*tendency(i,k,j) ) /rho_new(i,k,j) 
      ENDDO
      ENDDO
      ENDDO

      ENDDO
	  
   

END SUBROUTINE rk_update_scalar_pd_fm

!------------------------------------------------------------




SUBROUTINE init_zero_tendency_fm(u_tendf, v_tendf, w_tendf, pfm_tendf, rhofm_tendf,  &
                              t_tendf,  tke_tendf,             &
                         !      mu_tendf,           &
                              moist_tendf,chem_tendf,scalar_tendf,     &
                              tracer_tendf,n_tracer,                   &
                              n_moist,n_chem,n_scalar,rk_step,         &
                              ids, ide, jds, jde, kds, kde,            &
                              ims, ime, jms, jme, kms, kme,            &
                              its, ite, jts, jte, kts, kte             )
!-----------------------------------------------------------------------
   IMPLICIT NONE
!-----------------------------------------------------------------------

   INTEGER ,       INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                    ims, ime, jms, jme, kms, kme, &
                                    its, ite, jts, jte, kts, kte

   INTEGER ,       INTENT(IN   ) :: n_moist,n_chem,n_scalar,n_tracer,rk_step

   REAL , DIMENSION( ims:ime , kms:kme, jms:jme  ) , INTENT(INOUT) ::  &
                                                             u_tendf, &
                                                             v_tendf, &
                                                             w_tendf, &
                                                             pfm_tendf, rhofm_tendf,  &
                                                              t_tendf, &
                                                            tke_tendf

 !  REAL , DIMENSION( ims:ime , jms:jme  ) , INTENT(INOUT) ::  mu_tendf

   REAL , DIMENSION(ims:ime, kms:kme, jms:jme, n_moist),INTENT(INOUT)::&
                                                          moist_tendf

   REAL , DIMENSION(ims:ime, kms:kme, jms:jme, n_chem ),INTENT(INOUT)::&
                                                          chem_tendf
   REAL , DIMENSION(ims:ime, kms:kme, jms:jme, n_tracer ),INTENT(INOUT)::&
                                                          tracer_tendf

   REAL , DIMENSION(ims:ime, kms:kme, jms:jme, n_scalar ),INTENT(INOUT)::&
                                                          scalar_tendf

! LOCAL VARS

   INTEGER :: im, ic, is

!<DESCRIPTION>
!
! init_zero_tendency
! sets tendency arrays to zero for all prognostic variables.
!
!</DESCRIPTION>


   CALL zero_tend ( u_tendf,                        &
                    ids, ide, jds, jde, kds, kde,    &
                    ims, ime, jms, jme, kms, kme,    &
                    its, ite, jts, jte, kts, kte     )

   CALL zero_tend ( v_tendf,                        &
                    ids, ide, jds, jde, kds, kde,    &
                    ims, ime, jms, jme, kms, kme,    &
                    its, ite, jts, jte, kts, kte     )

   CALL zero_tend ( w_tendf,                        &
                    ids, ide, jds, jde, kds, kde,    &
                    ims, ime, jms, jme, kms, kme,    &
                    its, ite, jts, jte, kts, kte     )

   CALL zero_tend ( pfm_tendf,                        &
                    ids, ide, jds, jde, kds, kde,    &
                    ims, ime, jms, jme, kms, kme,    &
                    its, ite, jts, jte, kts, kte     )
					

   CALL zero_tend ( rhofm_tendf,                        &
                    ids, ide, jds, jde, kds, kde,    &
                    ims, ime, jms, jme, kms, kme,    &
                    its, ite, jts, jte, kts, kte     )
					

   CALL zero_tend ( t_tendf,                         &
                    ids, ide, jds, jde, kds, kde,    &
                    ims, ime, jms, jme, kms, kme,    &
                    its, ite, jts, jte, kts, kte     )

   CALL zero_tend ( tke_tendf,                       &
                    ids, ide, jds, jde, kds, kde,    &
                    ims, ime, jms, jme, kms, kme,    &
                    its, ite, jts, jte, kts, kte     )

 !  CALL zero_tend2d( mu_tendf,                        &
  !                  ids, ide, jds, jde, kds, kds,    &
  !                  ims, ime, jms, jme, kms, kms,    &
 !                   its, ite, jts, jte, kts, kts     )

!   DO im=PARAM_FIRST_SCALAR,n_moist
   DO im=1,n_moist                      ! make sure first one is zero too
      CALL zero_tend ( moist_tendf(ims,kms,jms,im),  &
                       ids, ide, jds, jde, kds, kde, &
                       ims, ime, jms, jme, kms, kme, &
                       its, ite, jts, jte, kts, kte  )
   ENDDO

!   DO ic=PARAM_FIRST_SCALAR,n_chem
   DO ic=1,n_chem                       ! make sure first one is zero too
      CALL zero_tend ( chem_tendf(ims,kms,jms,ic),   &
                       ids, ide, jds, jde, kds, kde, &
                       ims, ime, jms, jme, kms, kme, &
                       its, ite, jts, jte, kts, kte  )
   ENDDO

!   DO ic=PARAM_FIRST_SCALAR,n_tracer
   DO ic=1,n_tracer                     ! make sure first one is zero too
      CALL zero_tend ( tracer_tendf(ims,kms,jms,ic), &
                       ids, ide, jds, jde, kds, kde, &
                       ims, ime, jms, jme, kms, kme, &
                       its, ite, jts, jte, kts, kte  )
   ENDDO

!   DO ic=PARAM_FIRST_SCALAR,n_scalar
   DO ic=1,n_scalar                       ! make sure first one is zero too
      CALL zero_tend ( scalar_tendf(ims,kms,jms,ic),   &
                       ids, ide, jds, jde, kds, kde, &
                       ims, ime, jms, jme, kms, kme, &
                       its, ite, jts, jte, kts, kte  )
   ENDDO

END SUBROUTINE init_zero_tendency_fm

!===================================================================


SUBROUTINE dump_data( a, field, io_unit,            &
                      ims, ime, jms, jme, kms, kme, &
                      ids, ide, jds, jde, kds, kde )
implicit none
integer ::  ims, ime, jms, jme, kms, kme, &
            ids, ide, jds, jde, kds, kde
real, dimension(ims:ime, kms:kme, jds:jde) :: a
character :: field
integer :: io_unit

integer :: is,ie,js,je,ks,ke

!<DESCRIPTION
!
! quick and dirty debug io utility
!
!</DESCRIPTION

is = ids
ie = ide-1
js = jds
je = jde-1
ks = kds
ke = kde-1

if(field == 'u') ie = ide
if(field == 'v') je = jde
if(field == 'w') ke = kde

write(io_unit) is,ie,ks,ke,js,je
write(io_unit) a(is:ie, ks:ke, js:je)

end subroutine dump_data

!-----------------------------------------------------------------------

SUBROUTINE calculate_phy_tend (config_flags,c1,c2,                     &
                     mut,muu,muv,pi3d,                                 &
                     RTHRATEN,                                         &
                     RUBLTEN,RVBLTEN,RTHBLTEN,                         &
                     RQVBLTEN,RQCBLTEN,RQIBLTEN,                       &
                     RUCUTEN,RVCUTEN,RTHCUTEN,                         &
                     RQVCUTEN,RQCCUTEN,RQRCUTEN,                       &
                     RQICUTEN,RQSCUTEN,                                &
                     RUSHTEN,RVSHTEN,RTHSHTEN,                         &
                     RQVSHTEN,RQCSHTEN,RQRSHTEN,                       &
                     RQISHTEN,RQSSHTEN,RQGSHTEN,                       &
                     RUNDGDTEN,RVNDGDTEN,RTHNDGDTEN,RQVNDGDTEN,        &
                     RMUNDGDTEN,                                       &
                     scalar, scalar_tend, num_scalar,                  &
                     tracer, tracer_tend, num_tracer,                  &
                     ids,ide, jds,jde, kds,kde,                        &
                     ims,ime, jms,jme, kms,kme,                        &
                     its,ite, jts,jte, kts,kte                         )
!-----------------------------------------------------------------------
      IMPLICIT NONE

      TYPE(grid_config_rec_type), INTENT(IN)     ::      config_flags

      INTEGER,  INTENT(IN   )   ::          num_scalar, num_tracer
      INTEGER,  INTENT(IN   )   ::          ids,ide, jds,jde, kds,kde, &
                                            ims,ime, jms,jme, kms,kme, &
                                            its,ite, jts,jte, kts,kte

      REAL,     DIMENSION( ims:ime, kms:kme, jms:jme )               , &
                INTENT(IN   )   ::                               pi3d

      REAL,     DIMENSION( kms:kme )                                 , &
                INTENT(IN   )   ::                               c1, c2
                                                                 
      REAL,     DIMENSION( ims:ime, jms:jme )                        , &
                INTENT(IN   )   ::                                mut, &
                                                                  muu, &
                                                                  muv
      
                                                           
! radiation

      REAL,     DIMENSION( ims:ime, kms:kme, jms:jme ),                &
                INTENT(INOUT)   ::                           RTHRATEN

! cumulus

      REAL,     DIMENSION( ims:ime , kms:kme , jms:jme ),              &
                INTENT(INOUT)   ::                                     &
                                                              RUCUTEN, &
                                                              RVCUTEN, &
                                                             RTHCUTEN, &
                                                             RQVCUTEN, &
                                                             RQCCUTEN, &
                                                             RQRCUTEN, &
                                                             RQICUTEN, &
                                                             RQSCUTEN, &
                                                              RUSHTEN, &
                                                              RVSHTEN, &
                                                             RTHSHTEN, &
                                                             RQVSHTEN, &
                                                             RQCSHTEN, &
                                                             RQRSHTEN, &
                                                             RQISHTEN, &
                                                             RQSSHTEN, &
                                                             RQGSHTEN

! pbl

      REAL,     DIMENSION( ims:ime, kms:kme, jms:jme )               , &
                INTENT(INOUT)   ::                            RUBLTEN, &
                                                              RVBLTEN, &
                                                             RTHBLTEN, &
                                                             RQVBLTEN, &
                                                             RQCBLTEN, &
                                                             RQIBLTEN

! fdda

      REAL,     DIMENSION( ims:ime, kms:kme, jms:jme )               , &
                INTENT(INOUT)   ::                            RUNDGDTEN, &
                                                              RVNDGDTEN, &
                                                             RTHNDGDTEN, &
                                                             RQVNDGDTEN
      REAL,     DIMENSION( ims:ime, jms:jme )               , &
                INTENT(INOUT)   ::                           RMUNDGDTEN

! 4d arrays

    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_tracer),INTENT(INOUT)   :: tracer
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_tracer),INTENT(INOUT)   :: tracer_tend
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_scalar),INTENT(INOUT)   :: scalar
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_scalar),INTENT(INOUT)   :: scalar_tend

      INTEGER :: i,k,j, im
      INTEGER :: itf,ktf,jtf,itsu,jtsv

!-----------------------------------------------------------------------

!<DESCRIPTION>
!
!  calculate_phy_tend couples the physics tendencies to the column mass (mu),
!  because prognostic equations are in flux form, but physics tendencies are
!  computed for uncoupled variables.
!
!</DESCRIPTION>

      itf=MIN(ite,ide-1)
      jtf=MIN(jte,jde-1)
      ktf=MIN(kte,kde-1)
      itsu=MAX(its,ids+1)
      jtsv=MAX(jts,jds+1)

! radiation

   IF (config_flags%ra_lw_physics .gt. 0 .or. config_flags%ra_sw_physics .gt. 0) THEN

      DO J=jts,jtf
      DO K=kts,ktf
      DO I=its,itf
         RTHRATEN(I,K,J)=(c1(k)*mut(I,J)+c2(k))*RTHRATEN(I,K,J)
      ENDDO
      ENDDO
      ENDDO

   ENDIF

! cumulus

   IF (config_flags%cu_physics .gt. 0) THEN

      DO J=jts,jtf
      DO I=its,itf
      DO K=kts,ktf
         RUCUTEN(I,K,J) =(c1(k)*mut(I,J)+c2(k))*RUCUTEN(I,K,J)
         RVCUTEN(I,K,J) =(c1(k)*mut(I,J)+c2(k))*RVCUTEN(I,K,J)
         RTHCUTEN(I,K,J)=(c1(k)*mut(I,J)+c2(k))*RTHCUTEN(I,K,J)
         RQVCUTEN(I,K,J)=(c1(k)*mut(I,J)+c2(k))*RQVCUTEN(I,K,J)
      ENDDO
      ENDDO
      ENDDO

      IF (P_QC .ge. PARAM_FIRST_SCALAR)THEN
         DO J=jts,jtf
         DO I=its,itf
         DO K=kts,ktf
            RQCCUTEN(I,K,J)=(c1(k)*mut(I,J)+c2(k))*RQCCUTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF (P_QR .ge. PARAM_FIRST_SCALAR)THEN
         DO J=jts,jtf
         DO I=its,itf
         DO K=kts,ktf
            RQRCUTEN(I,K,J)=(c1(k)*mut(I,J)+c2(k))*RQRCUTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF (P_QI .ge. PARAM_FIRST_SCALAR)THEN
         DO J=jts,jtf
         DO I=its,itf
         DO K=kts,ktf
            RQICUTEN(I,K,J)=(c1(k)*mut(I,J)+c2(k))*RQICUTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF(P_QS .ge. PARAM_FIRST_SCALAR)THEN
         DO J=jts,jtf
         DO I=its,itf
         DO K=kts,ktf
            RQSCUTEN(I,K,J)=(c1(k)*mut(I,J)+c2(k))*RQSCUTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

   ENDIF

! shallow cumulus

   IF (config_flags%shcu_physics .gt. 0) THEN

      DO J=jts,jtf
      DO I=its,itf
      DO K=kts,ktf
         RUSHTEN(I,K,J) =(c1(k)*mut(I,J)+c2(k))*RUSHTEN(I,K,J)
         RVSHTEN(I,K,J) =(c1(k)*mut(I,J)+c2(k))*RVSHTEN(I,K,J)
         RTHSHTEN(I,K,J)=(c1(k)*mut(I,J)+c2(k))*RTHSHTEN(I,K,J)
         RQVSHTEN(I,K,J)=(c1(k)*mut(I,J)+c2(k))*RQVSHTEN(I,K,J)
      ENDDO
      ENDDO
      ENDDO

      IF (P_QC .ge. PARAM_FIRST_SCALAR)THEN
         DO J=jts,jtf
         DO I=its,itf
         DO K=kts,ktf
            RQCSHTEN(I,K,J)=(c1(k)*mut(I,J)+c2(k))*RQCSHTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF (P_QR .ge. PARAM_FIRST_SCALAR)THEN
         DO J=jts,jtf
         DO I=its,itf
         DO K=kts,ktf
            RQRSHTEN(I,K,J)=(c1(k)*mut(I,J)+c2(k))*RQRSHTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF (P_QI .ge. PARAM_FIRST_SCALAR)THEN
         DO J=jts,jtf
         DO I=its,itf
         DO K=kts,ktf
            RQISHTEN(I,K,J)=(c1(k)*mut(I,J)+c2(k))*RQISHTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF(P_QS .ge. PARAM_FIRST_SCALAR)THEN
         DO J=jts,jtf
         DO I=its,itf
         DO K=kts,ktf
            RQSSHTEN(I,K,J)=(c1(k)*mut(I,J)+c2(k))*RQSSHTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF(P_QG .ge. PARAM_FIRST_SCALAR)THEN
         DO J=jts,jtf
         DO I=its,itf
         DO K=kts,ktf
            RQGSHTEN(I,K,J)=(c1(k)*mut(I,J)+c2(k))*RQGSHTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

   ENDIF

! pbl

   IF (config_flags%bl_pbl_physics .gt. 0) THEN

      DO J=jts,jtf
      DO K=kts,ktf
      DO I=its,itf
         RUBLTEN(I,K,J) =(c1(k)*mut(I,J)+c2(k))*RUBLTEN(I,K,J)
         RVBLTEN(I,K,J) =(c1(k)*mut(I,J)+c2(k))*RVBLTEN(I,K,J)
         RTHBLTEN(I,K,J)=(c1(k)*mut(I,J)+c2(k))*RTHBLTEN(I,K,J)
      ENDDO
      ENDDO
      ENDDO

      IF (P_QV .ge. PARAM_FIRST_SCALAR) THEN
         DO J=jts,jtf
         DO K=kts,ktf
         DO I=its,itf
            RQVBLTEN(I,K,J)=(c1(k)*mut(I,J)+c2(k))*RQVBLTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF (P_QC .ge. PARAM_FIRST_SCALAR) THEN
         DO J=jts,jtf
         DO K=kts,ktf
         DO I=its,itf
           RQCBLTEN(I,K,J)=(c1(k)*mut(I,J)+c2(k))*RQCBLTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF (P_QI .ge. PARAM_FIRST_SCALAR) THEN
         DO J=jts,jtf
         DO K=kts,ktf
         DO I=its,itf
            RQIBLTEN(I,K,J)=(c1(k)*mut(I,J)+c2(k))*RQIBLTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

    ENDIF

! fdda
! note fdda u and v tendencies are staggered, also only interior points have muu/muv,
!   so only couple those

   IF (config_flags%grid_fdda .gt. 0) THEN

      DO J=jts,jtf
      DO K=kts,ktf
      DO I=itsu,itf
!     if( i == itf/2 .AND. j == jtf/2 .AND. k == ktf/2 ) &
!     write(*,'(a,3i6,e15.5)') 'u_ten before=',i,k,j, RUNDGDTEN(i,k,j)
         RUNDGDTEN(I,K,J) =(c1(k)*muu(I,J)+c2(k))*RUNDGDTEN(I,K,J)
!        if( i == itf/2 .AND. j == jtf/2 .AND. k==ktf/2 ) &
!          write(*,'(a,2f15.5)') 'mu, muu=',(c1(k)*mut(I,J)+c2(k)), (c1(k)*muu(i,j)+c2(k))
!     if( i == itf/2 .AND. j == jtf/2 .AND. k == ktf/2 ) &
!     write(*,'(a,3i6,e15.5)') 'u_ten after=',i,k,j, RUNDGDTEN(i,k,j)
!     if( RUNDGDTEN(i,k,j) > 30.0 ) write(*,*) 'IKJ=',i,k,j
      ENDDO
      ENDDO
      ENDDO
!     write(*,'(a,e15.5)') 'u_ten MAXIMUM after=', maxval(RUNDGDTEN)
      DO J=jtsv,jtf
      DO K=kts,ktf
      DO I=its,itf
         RVNDGDTEN(I,K,J) =(c1(k)*muv(I,J)+c2(k))*RVNDGDTEN(I,K,J)
      ENDDO
      ENDDO
      ENDDO
      DO J=jts,jtf
      DO K=kts,ktf
      DO I=its,itf
!     if( i == itf/2 .AND. j == jtf/2 .AND. k == ktf/2 ) &
!     write(*,'(a,3i6,e15.5)') 'th before=',i,k,j, RTHNDGDTEN(I,K,J)
         RTHNDGDTEN(I,K,J)=(c1(k)*mut(I,J)+c2(k))*RTHNDGDTEN(I,K,J)
!        RMUNDGDTEN(I,J) - no coupling
!     if( i == itf/2 .AND. j == jtf/2 .AND. k == ktf/2 ) &
!     write(*,'(a,3i6,e15.5)') 'th after=',i,k,j, RTHNDGDTEN(I,K,J)
      ENDDO
      ENDDO
      ENDDO

      IF (P_QV .ge. PARAM_FIRST_SCALAR) THEN
         DO J=jts,jtf
         DO K=kts,ktf
         DO I=its,itf
            RQVNDGDTEN(I,K,J)=(c1(k)*mut(I,J)+c2(k))*RQVNDGDTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

    ENDIF

! 4d couple scalar tendencies that have only uncoupled physics tendencies at this point

   DO im = PARAM_FIRST_SCALAR,num_scalar
         DO J=jts,jtf
         DO K=kts,ktf
         DO I=its,itf
            scalar_tend(I,K,J,im)=(c1(k)*mut(I,J)+c2(k))*scalar_tend(I,K,J,im)
         ENDDO
         ENDDO
         ENDDO
   ENDDO

   DO im = PARAM_FIRST_SCALAR,num_tracer
         DO J=jts,jtf
         DO K=kts,ktf
         DO I=its,itf
            tracer_tend(I,K,J,im)=(c1(k)*mut(I,J)+c2(k))*tracer_tend(I,K,J,im)
         ENDDO
         ENDDO
         ENDDO
   ENDDO


END SUBROUTINE calculate_phy_tend


!-----------------------------------------------------------------------

SUBROUTINE calculate_phy_tend_fm (config_flags,        & 
                        rho,             &                  ! hanshenyou
             !        c1,c2,                     &
             !        mut,muu,muv,pi3d,                                 &
                     RTHRATEN,                                         &
                     RUBLTEN,RVBLTEN,RTHBLTEN,                         &
                     RQVBLTEN,RQCBLTEN,RQIBLTEN,                       &
                     RUCUTEN,RVCUTEN,RTHCUTEN,                         &
                     RQVCUTEN,RQCCUTEN,RQRCUTEN,                       &
                     RQICUTEN,RQSCUTEN,                                &
                     RUSHTEN,RVSHTEN,RTHSHTEN,                         &
                     RQVSHTEN,RQCSHTEN,RQRSHTEN,                       &
                     RQISHTEN,RQSSHTEN,RQGSHTEN,                       &
                     RUNDGDTEN,RVNDGDTEN,RTHNDGDTEN,RQVNDGDTEN,        &
                     RMUNDGDTEN,                                       &
                     scalar, scalar_tend, num_scalar,                  &
                     tracer, tracer_tend, num_tracer,                  &
                     ids,ide, jds,jde, kds,kde,                        &
                     ims,ime, jms,jme, kms,kme,                        &
                     its,ite, jts,jte, kts,kte                         )
!-----------------------------------------------------------------------
      IMPLICIT NONE

      TYPE(grid_config_rec_type), INTENT(IN)     ::      config_flags

      INTEGER,  INTENT(IN   )   ::          num_scalar, num_tracer
      INTEGER,  INTENT(IN   )   ::          ids,ide, jds,jde, kds,kde, &
                                            ims,ime, jms,jme, kms,kme, &
                                            its,ite, jts,jte, kts,kte

 !     REAL,     DIMENSION( ims:ime, kms:kme, jms:jme )               , &
 !               INTENT(IN   )   ::                               pi3d

  !    REAL,     DIMENSION( kms:kme )                                 , &
  !              INTENT(IN   )   ::                               c1, c2
                                                                 
  !    REAL,     DIMENSION( ims:ime, jms:jme )                        , &
  !              INTENT(IN   )   ::                                mut, &
  !                                                                muu, &
  !                                                                muv
        REAL,     DIMENSION( ims:ime, kms:kme, jms:jme ),                &
                INTENT(IN)   ::                           rho  
                                                           
! radiation

      REAL,     DIMENSION( ims:ime, kms:kme, jms:jme ),                &
                INTENT(INOUT)   ::                           RTHRATEN

! cumulus

      REAL,     DIMENSION( ims:ime , kms:kme , jms:jme ),              &
                INTENT(INOUT)   ::                                     &
                                                              RUCUTEN, &
                                                              RVCUTEN, &
                                                             RTHCUTEN, &
                                                             RQVCUTEN, &
                                                             RQCCUTEN, &
                                                             RQRCUTEN, &
                                                             RQICUTEN, &
                                                             RQSCUTEN, &
                                                              RUSHTEN, &
                                                              RVSHTEN, &
                                                             RTHSHTEN, &
                                                             RQVSHTEN, &
                                                             RQCSHTEN, &
                                                             RQRSHTEN, &
                                                             RQISHTEN, &
                                                             RQSSHTEN, &
                                                             RQGSHTEN

! pbl

      REAL,     DIMENSION( ims:ime, kms:kme, jms:jme )               , &
                INTENT(INOUT)   ::                            RUBLTEN, &
                                                              RVBLTEN, &
                                                             RTHBLTEN, &
                                                             RQVBLTEN, &
                                                             RQCBLTEN, &
                                                             RQIBLTEN

! fdda

      REAL,     DIMENSION( ims:ime, kms:kme, jms:jme )               , &
                INTENT(INOUT)   ::                            RUNDGDTEN, &
                                                              RVNDGDTEN, &
                                                             RTHNDGDTEN, &
                                                             RQVNDGDTEN
      REAL,     DIMENSION( ims:ime, jms:jme )               , &
                INTENT(INOUT)   ::                           RMUNDGDTEN

! 4d arrays

    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_tracer),INTENT(INOUT)   :: tracer
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_tracer),INTENT(INOUT)   :: tracer_tend
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_scalar),INTENT(INOUT)   :: scalar
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_scalar),INTENT(INOUT)   :: scalar_tend

      INTEGER :: i,k,j, im
      INTEGER :: itf,ktf,jtf,itsu,jtsv

!-----------------------------------------------------------------------

!<DESCRIPTION>
!
!  calculate_phy_tend couples the physics tendencies to the column mass (mu),
!  because prognostic equations are in flux form, but physics tendencies are
!  computed for uncoupled variables.
!
!</DESCRIPTION>

      itf=MIN(ite,ide-1)
      jtf=MIN(jte,jde-1)
      ktf=MIN(kte,kde-1)
      itsu=MAX(its,ids+1)
      jtsv=MAX(jts,jds+1)

! radiation

   IF (config_flags%ra_lw_physics .gt. 0 .or. config_flags%ra_sw_physics .gt. 0) THEN

      DO J=jts,jtf
      DO K=kts,ktf
      DO I=its,itf
         RTHRATEN(I,K,J)=rho(i,k,j)*RTHRATEN(I,K,J)
      ENDDO
      ENDDO
      ENDDO

   ENDIF

! cumulus

   IF (config_flags%cu_physics .gt. 0) THEN

      DO J=jts,jtf
      DO I=its,itf
      DO K=kts,ktf
         RUCUTEN(I,K,J) =rho(i,k,j)*RUCUTEN(I,K,J)
         RVCUTEN(I,K,J) =rho(i,k,j)*RVCUTEN(I,K,J)
         RTHCUTEN(I,K,J)=rho(i,k,j)*RTHCUTEN(I,K,J)
         RQVCUTEN(I,K,J)=rho(i,k,j)*RQVCUTEN(I,K,J)
      ENDDO
      ENDDO
      ENDDO

      IF (P_QC .ge. PARAM_FIRST_SCALAR)THEN
         DO J=jts,jtf
         DO I=its,itf
         DO K=kts,ktf
            RQCCUTEN(I,K,J)=rho(i,k,j)*RQCCUTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF (P_QR .ge. PARAM_FIRST_SCALAR)THEN
         DO J=jts,jtf
         DO I=its,itf
         DO K=kts,ktf
            RQRCUTEN(I,K,J)=rho(i,k,j)*RQRCUTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF (P_QI .ge. PARAM_FIRST_SCALAR)THEN
         DO J=jts,jtf
         DO I=its,itf
         DO K=kts,ktf
            RQICUTEN(I,K,J)=rho(i,k,j)*RQICUTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF(P_QS .ge. PARAM_FIRST_SCALAR)THEN
         DO J=jts,jtf
         DO I=its,itf
         DO K=kts,ktf
            RQSCUTEN(I,K,J)=rho(i,k,j)*RQSCUTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

   ENDIF

! shallow cumulus

   IF (config_flags%shcu_physics .gt. 0) THEN

      DO J=jts,jtf
      DO I=its,itf
      DO K=kts,ktf
         RUSHTEN(I,K,J) =rho(i,k,j)*RUSHTEN(I,K,J)
         RVSHTEN(I,K,J) =rho(i,k,j)*RVSHTEN(I,K,J)
         RTHSHTEN(I,K,J)=rho(i,k,j)*RTHSHTEN(I,K,J)
         RQVSHTEN(I,K,J)=rho(i,k,j)*RQVSHTEN(I,K,J)
      ENDDO
      ENDDO
      ENDDO

      IF (P_QC .ge. PARAM_FIRST_SCALAR)THEN
         DO J=jts,jtf
         DO I=its,itf
         DO K=kts,ktf
            RQCSHTEN(I,K,J)=rho(i,k,j)*RQCSHTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF (P_QR .ge. PARAM_FIRST_SCALAR)THEN
         DO J=jts,jtf
         DO I=its,itf
         DO K=kts,ktf
            RQRSHTEN(I,K,J)=rho(i,k,j)*RQRSHTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF (P_QI .ge. PARAM_FIRST_SCALAR)THEN
         DO J=jts,jtf
         DO I=its,itf
         DO K=kts,ktf
            RQISHTEN(I,K,J)=rho(i,k,j)*RQISHTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF(P_QS .ge. PARAM_FIRST_SCALAR)THEN
         DO J=jts,jtf
         DO I=its,itf
         DO K=kts,ktf
            RQSSHTEN(I,K,J)=rho(i,k,j)*RQSSHTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF(P_QG .ge. PARAM_FIRST_SCALAR)THEN
         DO J=jts,jtf
         DO I=its,itf
         DO K=kts,ktf
            RQGSHTEN(I,K,J)=rho(i,k,j)*RQGSHTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

   ENDIF

! pbl

   IF (config_flags%bl_pbl_physics .gt. 0) THEN

      DO J=jts,jtf
      DO K=kts,ktf
      DO I=its,itf
         RUBLTEN(I,K,J) =rho(i,k,j)*RUBLTEN(I,K,J)
         RVBLTEN(I,K,J) =rho(i,k,j)*RVBLTEN(I,K,J)
         RTHBLTEN(I,K,J)=rho(i,k,j)*RTHBLTEN(I,K,J)
      ENDDO
      ENDDO
      ENDDO

      IF (P_QV .ge. PARAM_FIRST_SCALAR) THEN
         DO J=jts,jtf
         DO K=kts,ktf
         DO I=its,itf
            RQVBLTEN(I,K,J)=rho(i,k,j)*RQVBLTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF (P_QC .ge. PARAM_FIRST_SCALAR) THEN
         DO J=jts,jtf
         DO K=kts,ktf
         DO I=its,itf
           RQCBLTEN(I,K,J)=rho(i,k,j)*RQCBLTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF (P_QI .ge. PARAM_FIRST_SCALAR) THEN
         DO J=jts,jtf
         DO K=kts,ktf
         DO I=its,itf
            RQIBLTEN(I,K,J)=rho(i,k,j)*RQIBLTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

    ENDIF

! fdda
! note fdda u and v tendencies are staggered, also only interior points have muu/muv,
!   so only couple those

   IF (config_flags%grid_fdda .gt. 0) THEN

      DO J=jts,jtf
      DO K=kts,ktf
      DO I=itsu,itf
!     if( i == itf/2 .AND. j == jtf/2 .AND. k == ktf/2 ) &
!     write(*,'(a,3i6,e15.5)') 'u_ten before=',i,k,j, RUNDGDTEN(i,k,j)


   !      RUNDGDTEN(I,K,J) =(c1(k)*muu(I,J)+c2(k))*RUNDGDTEN(I,K,J)
		  RUNDGDTEN(I,K,J) = 0.5*(rho(i-1,k,j)+rho(i,k,j))*RUNDGDTEN(I,K,J)  ! hanshenyou 
		 
!        if( i == itf/2 .AND. j == jtf/2 .AND. k==ktf/2 ) &
!          write(*,'(a,2f15.5)') 'mu, muu=',rho(i,k,j), (c1(k)*muu(i,j)+c2(k))
!     if( i == itf/2 .AND. j == jtf/2 .AND. k == ktf/2 ) &
!     write(*,'(a,3i6,e15.5)') 'u_ten after=',i,k,j, RUNDGDTEN(i,k,j)
!     if( RUNDGDTEN(i,k,j) > 30.0 ) write(*,*) 'IKJ=',i,k,j
      ENDDO
      ENDDO
      ENDDO
!     write(*,'(a,e15.5)') 'u_ten MAXIMUM after=', maxval(RUNDGDTEN)
      DO J=jtsv,jtf
      DO K=kts,ktf
      DO I=its,itf
	  
     !    RVNDGDTEN(I,K,J) =(c1(k)*muv(I,J)+c2(k))*RVNDGDTEN(I,K,J)
		   RVNDGDTEN(I,K,J) = 0.5*(rho(i,k,j-1)+rho(i,k,j))*RVNDGDTEN(I,K,J) 
      ENDDO
      ENDDO
      ENDDO
	  
	  
      DO J=jts,jtf
      DO K=kts,ktf
      DO I=its,itf
!     if( i == itf/2 .AND. j == jtf/2 .AND. k == ktf/2 ) &
!     write(*,'(a,3i6,e15.5)') 'th before=',i,k,j, RTHNDGDTEN(I,K,J)
         RTHNDGDTEN(I,K,J)=rho(i,k,j)*RTHNDGDTEN(I,K,J)
!        RMUNDGDTEN(I,J) - no coupling
!     if( i == itf/2 .AND. j == jtf/2 .AND. k == ktf/2 ) &
!     write(*,'(a,3i6,e15.5)') 'th after=',i,k,j, RTHNDGDTEN(I,K,J)
      ENDDO
      ENDDO
      ENDDO

      IF (P_QV .ge. PARAM_FIRST_SCALAR) THEN
         DO J=jts,jtf
         DO K=kts,ktf
         DO I=its,itf
            RQVNDGDTEN(I,K,J)=rho(i,k,j)*RQVNDGDTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

    ENDIF

! 4d couple scalar tendencies that have only uncoupled physics tendencies at this point

   DO im = PARAM_FIRST_SCALAR,num_scalar
         DO J=jts,jtf
         DO K=kts,ktf
         DO I=its,itf
            scalar_tend(I,K,J,im)=rho(i,k,j)*scalar_tend(I,K,J,im)
         ENDDO
         ENDDO
         ENDDO
   ENDDO

   DO im = PARAM_FIRST_SCALAR,num_tracer
         DO J=jts,jtf
         DO K=kts,ktf
         DO I=its,itf
            tracer_tend(I,K,J,im)=rho(i,k,j)*tracer_tend(I,K,J,im)
         ENDDO
         ENDDO
         ENDDO
   ENDDO
 

END SUBROUTINE calculate_phy_tend_fm

!-----------------------------------------------------------------------

! only moist qv qc and so on 
SUBROUTINE calculate_phy_tend_fm2 (config_flags,        & 
                        rho,             &                  ! hanshenyou
             !        c1,c2,                     &
             !        mut,muu,muv,pi3d,                                 &
                     RTHRATEN,                                         &
                     RUBLTEN,RVBLTEN,RTHBLTEN,                         &
                     RQVBLTEN,RQCBLTEN,RQIBLTEN,                       &
                     RUCUTEN,RVCUTEN,RTHCUTEN,                         &
                     RQVCUTEN,RQCCUTEN,RQRCUTEN,                       &
                     RQICUTEN,RQSCUTEN,                                &
                     RUSHTEN,RVSHTEN,RTHSHTEN,                         &
                     RQVSHTEN,RQCSHTEN,RQRSHTEN,                       &
                     RQISHTEN,RQSSHTEN,RQGSHTEN,                       &
                     RUNDGDTEN,RVNDGDTEN,RTHNDGDTEN,RQVNDGDTEN,        &
                     RMUNDGDTEN,                                       &
                     scalar, scalar_tend, num_scalar,                  &
                     tracer, tracer_tend, num_tracer,                  &
                     ids,ide, jds,jde, kds,kde,                        &
                     ims,ime, jms,jme, kms,kme,                        &
                     its,ite, jts,jte, kts,kte                         )
!-----------------------------------------------------------------------
      IMPLICIT NONE

      TYPE(grid_config_rec_type), INTENT(IN)     ::      config_flags

      INTEGER,  INTENT(IN   )   ::          num_scalar, num_tracer
      INTEGER,  INTENT(IN   )   ::          ids,ide, jds,jde, kds,kde, &
                                            ims,ime, jms,jme, kms,kme, &
                                            its,ite, jts,jte, kts,kte

 !     REAL,     DIMENSION( ims:ime, kms:kme, jms:jme )               , &
 !               INTENT(IN   )   ::                               pi3d

  !    REAL,     DIMENSION( kms:kme )                                 , &
  !              INTENT(IN   )   ::                               c1, c2
                                                                 
  !    REAL,     DIMENSION( ims:ime, jms:jme )                        , &
  !              INTENT(IN   )   ::                                mut, &
  !                                                                muu, &
  !                                                                muv
        REAL,     DIMENSION( ims:ime, kms:kme, jms:jme ),                &
                INTENT(IN)   ::                           rho  
                                                           
! radiation

      REAL,     DIMENSION( ims:ime, kms:kme, jms:jme ),                &
                INTENT(INOUT)   ::                           RTHRATEN

! cumulus

      REAL,     DIMENSION( ims:ime , kms:kme , jms:jme ),              &
                INTENT(INOUT)   ::                                     &
                                                              RUCUTEN, &
                                                              RVCUTEN, &
                                                             RTHCUTEN, &
                                                             RQVCUTEN, &
                                                             RQCCUTEN, &
                                                             RQRCUTEN, &
                                                             RQICUTEN, &
                                                             RQSCUTEN, &
                                                              RUSHTEN, &
                                                              RVSHTEN, &
                                                             RTHSHTEN, &
                                                             RQVSHTEN, &
                                                             RQCSHTEN, &
                                                             RQRSHTEN, &
                                                             RQISHTEN, &
                                                             RQSSHTEN, &
                                                             RQGSHTEN

! pbl

      REAL,     DIMENSION( ims:ime, kms:kme, jms:jme )               , &
                INTENT(INOUT)   ::                            RUBLTEN, &
                                                              RVBLTEN, &
                                                             RTHBLTEN, &
                                                             RQVBLTEN, &
                                                             RQCBLTEN, &
                                                             RQIBLTEN

! fdda

      REAL,     DIMENSION( ims:ime, kms:kme, jms:jme )               , &
                INTENT(INOUT)   ::                            RUNDGDTEN, &
                                                              RVNDGDTEN, &
                                                             RTHNDGDTEN, &
                                                             RQVNDGDTEN
      REAL,     DIMENSION( ims:ime, jms:jme )               , &
                INTENT(INOUT)   ::                           RMUNDGDTEN

! 4d arrays

    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_tracer),INTENT(INOUT)   :: tracer
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_tracer),INTENT(INOUT)   :: tracer_tend
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_scalar),INTENT(INOUT)   :: scalar
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_scalar),INTENT(INOUT)   :: scalar_tend

      INTEGER :: i,k,j, im
      INTEGER :: itf,ktf,jtf,itsu,jtsv

!-----------------------------------------------------------------------

!<DESCRIPTION>
!
!  calculate_phy_tend couples the physics tendencies to the column mass (mu),
!  because prognostic equations are in flux form, but physics tendencies are
!  computed for uncoupled variables.
!
!</DESCRIPTION>

      itf=MIN(ite,ide-1)
      jtf=MIN(jte,jde-1)
      ktf=MIN(kte,kde-1)
      itsu=MAX(its,ids+1)
      jtsv=MAX(jts,jds+1)

! radiation

 

! cumulus

   IF (config_flags%cu_physics .gt. 0) THEN

      DO J=jts,jtf
      DO I=its,itf
      DO K=kts,ktf
 
         RQVCUTEN(I,K,J)=rho(i,k,j)*RQVCUTEN(I,K,J)
      ENDDO
      ENDDO
      ENDDO

      IF (P_QC .ge. PARAM_FIRST_SCALAR)THEN
         DO J=jts,jtf
         DO I=its,itf
         DO K=kts,ktf
            RQCCUTEN(I,K,J)=rho(i,k,j)*RQCCUTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF (P_QR .ge. PARAM_FIRST_SCALAR)THEN
         DO J=jts,jtf
         DO I=its,itf
         DO K=kts,ktf
            RQRCUTEN(I,K,J)=rho(i,k,j)*RQRCUTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF (P_QI .ge. PARAM_FIRST_SCALAR)THEN
         DO J=jts,jtf
         DO I=its,itf
         DO K=kts,ktf
            RQICUTEN(I,K,J)=rho(i,k,j)*RQICUTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF(P_QS .ge. PARAM_FIRST_SCALAR)THEN
         DO J=jts,jtf
         DO I=its,itf
         DO K=kts,ktf
            RQSCUTEN(I,K,J)=rho(i,k,j)*RQSCUTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

   ENDIF

! shallow cumulus

   IF (config_flags%shcu_physics .gt. 0) THEN

      DO J=jts,jtf
      DO I=its,itf
      DO K=kts,ktf
 
         RQVSHTEN(I,K,J)=rho(i,k,j)*RQVSHTEN(I,K,J)
      ENDDO
      ENDDO
      ENDDO

      IF (P_QC .ge. PARAM_FIRST_SCALAR)THEN
         DO J=jts,jtf
         DO I=its,itf
         DO K=kts,ktf
            RQCSHTEN(I,K,J)=rho(i,k,j)*RQCSHTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF (P_QR .ge. PARAM_FIRST_SCALAR)THEN
         DO J=jts,jtf
         DO I=its,itf
         DO K=kts,ktf
            RQRSHTEN(I,K,J)=rho(i,k,j)*RQRSHTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF (P_QI .ge. PARAM_FIRST_SCALAR)THEN
         DO J=jts,jtf
         DO I=its,itf
         DO K=kts,ktf
            RQISHTEN(I,K,J)=rho(i,k,j)*RQISHTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF(P_QS .ge. PARAM_FIRST_SCALAR)THEN
         DO J=jts,jtf
         DO I=its,itf
         DO K=kts,ktf
            RQSSHTEN(I,K,J)=rho(i,k,j)*RQSSHTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF(P_QG .ge. PARAM_FIRST_SCALAR)THEN
         DO J=jts,jtf
         DO I=its,itf
         DO K=kts,ktf
            RQGSHTEN(I,K,J)=rho(i,k,j)*RQGSHTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

   ENDIF

! pbl

   IF (config_flags%bl_pbl_physics .gt. 0) THEN

 

      IF (P_QV .ge. PARAM_FIRST_SCALAR) THEN
         DO J=jts,jtf
         DO K=kts,ktf
         DO I=its,itf
            RQVBLTEN(I,K,J)=rho(i,k,j)*RQVBLTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF (P_QC .ge. PARAM_FIRST_SCALAR) THEN
         DO J=jts,jtf
         DO K=kts,ktf
         DO I=its,itf
           RQCBLTEN(I,K,J)=rho(i,k,j)*RQCBLTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF (P_QI .ge. PARAM_FIRST_SCALAR) THEN
         DO J=jts,jtf
         DO K=kts,ktf
         DO I=its,itf
            RQIBLTEN(I,K,J)=rho(i,k,j)*RQIBLTEN(I,K,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

    ENDIF

! fdda
! note fdda u and v tendencies are staggered, also only interior points have muu/muv,
!   so only couple those

 

! 4d couple scalar tendencies that have only uncoupled physics tendencies at this point

   DO im = PARAM_FIRST_SCALAR,num_scalar
         DO J=jts,jtf
         DO K=kts,ktf
         DO I=its,itf
            scalar_tend(I,K,J,im)=rho(i,k,j)*scalar_tend(I,K,J,im)
         ENDDO
         ENDDO
         ENDDO
   ENDDO

   DO im = PARAM_FIRST_SCALAR,num_tracer
         DO J=jts,jtf
         DO K=kts,ktf
         DO I=its,itf
            tracer_tend(I,K,J,im)=rho(i,k,j)*tracer_tend(I,K,J,im)
         ENDDO
         ENDDO
         ENDDO
   ENDDO
 

END SUBROUTINE calculate_phy_tend_fm2
!-----------------------------------------------------------------------

SUBROUTINE positive_definite_filter ( a,                          &
                                      ids,ide, jds,jde, kds,kde,  &
                                      ims,ime, jms,jme, kms,kme,  &
                                      its,ite, jts,jte, kts,kte  )

  IMPLICIT NONE

  INTEGER,  INTENT(IN   )   ::          ids,ide, jds,jde, kds,kde, &
                                        ims,ime, jms,jme, kms,kme, &
                                        its,ite, jts,jte, kts,kte

  REAL, DIMENSION( ims:ime , kms:kme , jms:jme  ), INTENT(INOUT) :: a

  INTEGER :: i,k,j

!<DESCRIPTION>
!
! debug and testing code for bounding a variable
!
!</DESCRIPTION>

  DO j=jts,min(jte,jde-1)
  DO k=kts,kte-1
  DO i=its,min(ite,ide-1)
!    a(i,k,j) = max(a(i,k,j),0.)
    a(i,k,j) = min(1000.,max(a(i,k,j),0.))
  ENDDO
  ENDDO
  ENDDO

  END SUBROUTINE positive_definite_filter

!-----------------------------------------------------------------------

SUBROUTINE bound_tke ( tke, tke_upper_bound,       &
                       ids,ide, jds,jde, kds,kde,  &
                       ims,ime, jms,jme, kms,kme,  &
                       its,ite, jts,jte, kts,kte  )

  IMPLICIT NONE

  INTEGER,  INTENT(IN   )   ::          ids,ide, jds,jde, kds,kde, &
                                        ims,ime, jms,jme, kms,kme, &
                                        its,ite, jts,jte, kts,kte

  REAL, DIMENSION( ims:ime , kms:kme , jms:jme  ), INTENT(INOUT) :: tke
  REAL, INTENT(   IN) :: tke_upper_bound

  INTEGER :: i,k,j

!<DESCRIPTION>
!
! bounds tke between zero and tke_upper_bound.
!
!</DESCRIPTION>

  DO j=jts,min(jte,jde-1)
  DO k=kts,kte-1
  DO i=its,min(ite,ide-1)
    tke(i,k,j) = min(tke_upper_bound,max(tke(i,k,j),0.))
  ENDDO
  ENDDO
  ENDDO

  END SUBROUTINE bound_tke

!----------------------------------------------------------------------------------
!cyl: Implement the forward Lagrangian trajectory calculation in WRF
!Chiaying Lee RSMAS/UM
!----------------------------------------------------------------------------------
subroutine trajectory (     grid,config_flags,               &
                            dt,itimestep,ru_m, rv_m, ww_m,   &
                            mut,muu,muv,c1h,c2h,c1f,c2f,     &
                            rdx, rdy, rdn, rdnw,rdzw,        &
                            traj_i,traj_j,traj_k,            &
                            traj_long,traj_lat,              &
                            xlong,xlat,                      &
                            msft,msfu,msfv,                  &
                            ids, ide, jds, jde, kds, kde,    &
                            ims, ime, jms, jme, kms, kme,    &
                            its, ite, jts, jte, kts, kte    )

!---------------------------------------------------------------------------------------------------

  use module_date_time
  use module_utility
  use module_domain, only  : domain_clock_get
  use module_trajectory, only  : traject, traj_cnt

  implicit none
!---------------------------------------------------------------------------------------------------
! Subroutine trajectory calculates forward Lagrangian trajectory  
! of the selecting points. The trajectory of each point is subjected to u, v, and w.
! (Lee and Chen 2013, MWR).     
!
! The trajectories is initialized with given longitude (degree), latitude (degree), and height(eta level).
!
!--traj_i grid number in x direction (on mass grid), float
!--traj_j grid number in y direction (on mass grid), float
!--traj_k grid number in z direction (on mass grid), float
!--traj_long longitude of trajectories, float
!--traj_lat longitude of trajectories, float
!
!--------------------------------------------------------------------------------------------------
  TYPE(proj_info) :: proj
  TYPE(domain), INTENT(IN) :: grid
  TYPE (grid_config_rec_type) , INTENT(IN)          :: config_flags
  INTEGER ,                INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                             ims, ime, jms, jme, kms, kme, &
                                             its, ite, jts, jte, kts, kte
   INTEGER ,                                      INTENT(IN   ) :: itimestep
   REAL ,                                      INTENT(IN   ) :: rdx,     &
                                                                rdy,     &
                                                                dt
   REAL , DIMENSION( kms:kme ) ,               INTENT(IN   ) :: rdn,     &
                                                                rdnw

   REAL , DIMENSION( kms:kme ) ,               INTENT(IN   ) :: c1h, c2h, c1f, c2f

   REAL , DIMENSION( ims:ime , kms:kme, jms:jme  ) ,              &
                                        INTENT(IN   ) :: ru_m,     &
                                                         rv_m,     &
                                                         ww_m
   REAL , DIMENSION( ims:ime , jms:jme  ) ,              &
                                        INTENT(IN   ) ::xlong, xlat

  real, dimension(1:config_flags%num_traj), intent(inout) :: traj_i,traj_j,traj_k
  real, dimension(1:config_flags%num_traj), intent(inout) :: traj_long,traj_lat
  real, dimension(ims:ime,kms:kme,jms:jme),intent(in) :: rdzw
  real, dimension(ims:ime,kms:kme,jms:jme)::u,v,w
  real, dimension(ims:ime,jms:jme),intent(in)::msft,msfu,msfv
  real, dimension(ims:ime,jms:jme),intent(in)::muu,muv,mut
  integer :: j,k
  integer :: i_beg,i_end
  integer :: i_traj,j_traj,k_traj,tjk
  integer :: dm
  real :: traj_u,traj_v,traj_w
  real :: rdx_grid,rdy_grid,rdz_grid
  real :: deltx, delty, deltz,ax
  real :: const1
  real :: temp_i,temp_j
  integer :: i_u,j_v,k_w
  real :: u_a,u_b,u_c,u_d,v_a,v_b,v_c,v_d,w_a,w_b,w_c,w_d
  real :: d_a,d_b,d_c,d_d
  real :: u_temp_upper,u_temp_lower
  real :: v_temp_upper,v_temp_lower
  real :: w_temp_upper,w_temp_lower
  real :: eta_old, eta_new
  integer :: keta, keta_temp
  character(len=19)  :: current_timestr, next_timestr, wrk_timestr
  character(len=256) :: dbg_mes
  logical :: has_proj_map
! varalbe for map projectory
  real:: known_lat, known_lon
  TYPE (grid_config_rec_type) :: config_flags_temp
  TYPE(WRFU_Time) :: current_time, next_time
  TYPE(WRFU_Time) :: start_time, stop_time

  config_flags_temp = config_flags
  call trajmapproj(grid, config_flags_temp,proj)


! convert ru_m, rv_m and ww_m in u,v,w
   const1=1.0/2.0/sqrt(2.0)

   dm = grid%id

   write(dbg_mes,'(''trajectory('',i2.2,''): Entering trajectory'')') dm
   call wrf_debug(200,trim(dbg_mes) )

   i_beg = max( 1,its-1 )
   i_end = min( ite+2,ide-1 )
   do j=max(1,jts-1),min(jte+2,jde-1)
     do k=kms,kme-1
       u(its:i_end,k,j)=ru_m(its:i_end,k,j)/(c1h(k)*MUU(its:i_end,j)+c2h(k))*msfu(its:i_end,j)
       v(its:i_end,k,j)=rv_m(its:i_end,k,j)/(c1h(k)*MUV(its:i_end,j)+c2h(k))*msfv(its:i_end,j)
     enddo
     do k=kms,kme
       w(its:i_end,k,j)=ww_m(its:i_end,k,j)/(c1f(k)*MUT(its:i_end,j)+c2f(k))*msft(its:i_end,j)
     enddo
   enddo
          has_proj_map  = &
           ( proj%code .EQ. PROJ_LC           ) .OR. &
           ( proj%code .EQ. PROJ_PS_WGS84     ) .OR. &
           ( proj%code .EQ. PROJ_ALBERS_NAD83 ) .OR. &
           ( proj%code .EQ. PROJ_MERC         ) .OR. &
           ( proj%code .EQ. PROJ_LATLON       ) .OR. &
           ( proj%code .EQ. PROJ_CYL          ) .OR. &
           ( proj%code .EQ. PROJ_CASSINI      ) .OR. &
           ( proj%code .EQ. PROJ_GAUSS        ) .OR. &
           ( proj%code .EQ. PROJ_ROTLL        )

traj_loop: &
do tjk = 1,traj_cnt(dm)
!do tjk = 1,config_flags%num_traj
traj_is_active: &
    if (traj_i(tjk) .ne. -9999.0) then
      eta_old = 0.0
      eta_new = 0.0
      keta      = 0
      keta_temp = 0
      if( has_proj_map ) then
        call latlon_to_ij (proj, &
traj_lat(tjk),traj_long(tjk),traj_i(tjk),traj_j(tjk))
      end if

      i_traj = floor(traj_i(tjk)) ! find the lower_left_bottom corner for
!trajectory
      j_traj = floor(traj_j(tjk)) !
      k_traj = floor(traj_k(tjk)) !
traj_in_domain: &
      if ((i_traj .ge. its .and. i_traj .le. ite .and. i_traj .lt. ide) .and. &
          (j_traj .ge. jts .and. j_traj .le. jte .and. j_traj .lt. jde) .and. &
          (k_traj .le. kte .and. k_traj .lt. kde)) then
!-----------------------------------------------------------------------------
!  is trajectory in time interval?
!-----------------------------------------------------------------------------
         call domain_clock_get( grid, current_time=current_time, &
current_timestr=current_timestr)
         call geth_newdate( next_timestr, current_timestr, int(grid%dt) )
         call wrf_atotime( next_timestr, next_time )
         wrk_timestr(1:19) = traject(tjk,dm)%start_time(1:19)
!        write(*,*) ' '
!        write(*,*) traject(tjk,dm)
!        write(dbg_mes,'(''trajectory('',i2,2,''): tjk,start_time = '',i3,1x,a)') &
!                     dm,tjk,wrk_timestr
!        call wrf_debug( 200,trim(dbg_mes) )
         call wrf_atotime( wrk_timestr(1:19), start_time )
         wrk_timestr(1:19) = traject(tjk,dm)%stop_time(1:19)
         call wrf_atotime( wrk_timestr(1:19), stop_time )
is_in_time_interval: &
         if( next_time .ge. start_time .and. next_time .le. stop_time &
                                       .and. .not. traject(tjk,dm)%is_stationary ) then
! for u : check x stagger
        if (traj_i(tjk)-real(floor(traj_i(tjk))) .ge. 0.5 ) then
          i_u=floor(traj_i(tjk)) + 1
        else
          i_u=floor(traj_i(tjk))
        endif
! for layer k_traj
      if (k_traj .ge. 1 ) then
        u_a=u(i_u  ,k_traj,j_traj+1)
        u_b=u(i_u  ,k_traj,j_traj  )
        u_c=u(i_u+1,k_traj,j_traj  )
        u_d=u(i_u+1,k_traj,j_traj+1)
      else
        u_a=0.0
        u_b=0.0
        u_c=0.0
        u_d=0.0
      endif
        d_a=abs((real(i_u+1)-(traj_i(tjk)+0.5))*(real(j_traj  )-traj_j(tjk)))
        d_b=abs((real(i_u+1)-(traj_i(tjk)+0.5))*(real(j_traj+1)-traj_j(tjk)))
        d_c=abs((real(i_u  )-(traj_i(tjk)+0.5))*(real(j_traj+1)-traj_j(tjk)))
        d_d=abs((real(i_u  )-(traj_i(tjk)+0.5))*(real(j_traj  )-traj_j(tjk)))
        u_temp_lower=(u_a*d_a+u_b*d_b+u_c*d_c+u_d*d_d)/(d_a+d_b+d_c+d_d)

!for layer k_traj+1
        u_a=u(i_u  ,k_traj+1,j_traj+1)
        u_b=u(i_u  ,k_traj+1,j_traj  )
        u_c=u(i_u+1,k_traj+1,j_traj  )
        u_d=u(i_u+1,k_traj+1,j_traj+1)
        d_a=abs((real(i_u+1)-(traj_i(tjk)+0.5))*(real(j_traj  )-traj_j(tjk)))
        d_b=abs((real(i_u+1)-(traj_i(tjk)+0.5))*(real(j_traj+1)-traj_j(tjk)))
        d_c=abs((real(i_u  )-(traj_i(tjk)+0.5))*(real(j_traj+1)-traj_j(tjk)))
        d_d=abs((real(i_u  )-(traj_i(tjk)+0.5))*(real(j_traj  )-traj_j(tjk)))
        u_temp_upper=(u_a*d_a+u_b*d_b+u_c*d_c+u_d*d_d)/(d_a+d_b+d_c+d_d)

        traj_u=u_temp_upper*abs(real(k_traj)-traj_k(tjk))+u_temp_lower*abs(real(k_traj+1)-traj_k(tjk))

! for v: check y-stagger
        if (traj_j(tjk)-real(floor(traj_j(tjk))) .ge. 0.5 ) then
          j_v=floor(traj_j(tjk)) + 1
        else
          j_v=floor(traj_j(tjk))
        endif
! for layer k_traj
      if (k_traj .ge. 1 ) then
        v_a=v(i_traj  ,k_traj,j_v+1)
        v_b=v(i_traj  ,k_traj,j_v  )
        v_c=v(i_traj+1,k_traj,j_v  )
        v_d=v(i_traj+1,k_traj,j_v+1)
      else
        v_a=0.0
        v_b=0.0
        v_c=0.0
        v_d=0.0
      endif
        d_a=abs((real(i_traj+1)-traj_i(tjk))*(real(j_v  )-(traj_j(tjk)+0.5)))
        d_b=abs((real(i_traj+1)-traj_i(tjk))*(real(j_v+1)-(traj_j(tjk)+0.5)))
        d_c=abs((real(i_traj  )-traj_i(tjk))*(real(j_v+1)-(traj_j(tjk)+0.5)))
        d_d=abs((real(i_traj  )-traj_i(tjk))*(real(j_v  )-(traj_j(tjk)+0.5)))
        v_temp_lower=(v_a*d_a+v_b*d_b+v_c*d_c+v_d*d_d)/(d_a+d_b+d_c+d_d)
!for layer k_traj+1
        v_a=v(i_traj  ,k_traj+1,j_v+1)
        v_b=v(i_traj  ,k_traj+1,j_v  )
        v_c=v(i_traj+1,k_traj+1,j_v  )
        v_d=v(i_traj+1,k_traj+1,j_v+1)
        d_a=abs((real(i_traj+1)-traj_i(tjk))*(real(j_v  )-(traj_j(tjk)+0.5)))
        d_b=abs((real(i_traj+1)-traj_i(tjk))*(real(j_v+1)-(traj_j(tjk)+0.5)))
        d_c=abs((real(i_traj  )-traj_i(tjk))*(real(j_v+1)-(traj_j(tjk)+0.5)))
        d_d=abs((real(i_traj  )-traj_i(tjk))*(real(j_v  )-(traj_j(tjk)+0.5)))
        v_temp_upper=(v_a*d_a+v_b*d_b+v_c*d_c+v_d*d_d)/(d_a+d_b+d_c+d_d)

        traj_v=v_temp_upper*abs(real(k_traj)-traj_k(tjk))+v_temp_lower*abs(real(k_traj+1)-traj_k(tjk))


!for w: check for z-stagger
        if (traj_k(tjk)-real(floor(traj_k(tjk))) .ge. 0.5 ) then
          k_w=floor(traj_k(tjk)) + 1
        else
          k_w=floor(traj_k(tjk))
        endif
!for layer j_traj
        if (k_w .ge. 1) then
          w_b=w(i_traj  ,k_w  ,j_traj)
          w_c=w(i_traj+1,k_w  ,j_traj)
        else
          w_b=0.0
          w_c=0.0
        endif
        w_a=w(i_traj  ,k_w+1,j_traj)
        w_d=w(i_traj+1,k_w+1,j_traj)
        d_a=abs((real(i_traj+1)-traj_i(tjk))*(real(k_w  )-(traj_k(tjk)+0.5)))
        d_b=abs((real(i_traj+1)-traj_i(tjk))*(real(k_w+1)-(traj_k(tjk)+0.5)))
        d_c=abs((real(i_traj  )-traj_i(tjk))*(real(k_w+1)-(traj_k(tjk)+0.5)))
        d_d=abs((real(i_traj  )-traj_i(tjk))*(real(k_w  )-(traj_k(tjk)+0.5)))
        w_temp_lower=(w_a*d_a+w_b*d_b+w_c*d_c+w_d*d_d)/(d_a+d_b+d_c+d_d)
!for layer j_traj+1
        if (k_w .ge. 1) then
          w_b=w(i_traj  ,k_w  ,j_traj+1)
          w_c=w(i_traj+1,k_w  ,j_traj+1)
        else
          w_b=0.0
          w_c=0.0
        endif
        w_a=w(i_traj  ,k_w+1,j_traj+1)
        w_d=w(i_traj+1,k_w+1,j_traj+1)
        d_a=abs((real(i_traj+1)-traj_i(tjk))*(real(k_w  )-(traj_k(tjk)+0.5)))
        d_b=abs((real(i_traj+1)-traj_i(tjk))*(real(k_w+1)-(traj_k(tjk)+0.5)))
        d_c=abs((real(i_traj  )-traj_i(tjk))*(real(k_w+1)-(traj_k(tjk)+0.5)))
        d_d=abs((real(i_traj  )-traj_i(tjk))*(real(k_w  )-(traj_k(tjk)+0.5)))
        w_temp_upper=(w_a*d_a+w_b*d_b+w_c*d_c+w_d*d_d)/(d_a+d_b+d_c+d_d)

        traj_w=w_temp_upper*abs(real(j_traj)-traj_j(tjk))+w_temp_lower*abs(real(j_traj+1)-traj_j(tjk))
!get old eta
        eta_old=grid%znw(k_w+1)*abs(traj_k(tjk)+0.5-real(k_w))+grid%znw(k_w)*abs(traj_k(tjk)+0.5-real(k_w+1))

        rdx_grid=rdx*msft(i_traj,j_traj)! 1/(dx/msft)
        rdy_grid=rdy*msft(i_traj,j_traj)
        rdz_grid=rdnw(k_traj)
        deltx=traj_u*DT
        delty=traj_v*DT
        deltz=traj_w*DT
        traj_i(tjk)=traj_i(tjk)+deltx*rdx_grid
        traj_j(tjk)=traj_j(tjk)+delty*rdy_grid
        eta_new=eta_old+deltz
! get new traj_k(tjk)
        keta_temp = 0
        do keta=1, kme-1
          if (eta_new .le. grid%znw(keta) .and. eta_new .gt. grid%znw(keta+1)) then
           keta_temp=keta
          endif
        enddo
        if (keta_temp .eq. 0)  then
            traj_k(tjk) = traj_k(tjk)
        else
            traj_k(tjk) = (real(keta_temp)*abs(eta_new-grid%znw(keta_temp+1))+ &
                           real(keta_temp+1)*abs(eta_new-grid%znw(keta_temp))) &
                           /(grid%znw(keta_temp)-grid%znw(keta_temp+1))
            traj_k(tjk) = traj_k(tjk)-0.5
        endif
!! convert i,j,k into lon, lat
        if( has_proj_map ) then
          call ij_to_latlon (proj, traj_i(tjk), &
                             traj_j(tjk),traj_lat(tjk),traj_long(tjk))
        endif
      endif is_in_time_interval
     else traj_in_domain
        traj_i(tjk) = -9999.0
        traj_j(tjk) = -9999.0
        traj_k(tjk) = -9999.0
        traj_long(tjk) = -9999.0
        traj_lat(tjk) = -9999.0
     endif traj_in_domain
 endif traj_is_active
     traj_i(tjk) = wrf_dm_max_real(traj_i(tjk))
     traj_j(tjk) = wrf_dm_max_real(traj_j(tjk))
     traj_k(tjk) = wrf_dm_max_real(traj_k(tjk))
     traj_long(tjk) = wrf_dm_max_real(traj_long(tjk))
     traj_lat(tjk)  = wrf_dm_max_real(traj_lat(tjk))
enddo traj_loop
!end trajectory

END SUBROUTINE trajectory
!------------------------------------------------------------------------
!cyl:Implement the map prpjection code for trajectory calculation
!                                    Chiaying Lee  RSMAS/UM        
!------------------------------------------------------------------------
subroutine trajmapproj (grid,config_flags,ts_proj)
!------------------------------------------------------------------------
!!! This code calculates map-projection of trajectories. It is from share/wrf_timeseries.F
!------------------------------------------------------------------------

   IMPLICIT NONE


   ! Arguments
   TYPE (domain), INTENT(IN) :: grid
   TYPE (grid_config_rec_type) , INTENT(IN)  :: config_flags
   ! Externals
   LOGICAL, EXTERNAL :: wrf_dm_on_monitor
   INTEGER, EXTERNAL :: get_unused_unit


   ! Local variables
   INTEGER :: ntsloc_temp
   INTEGER :: i, k, iunit
   REAL :: ts_rx, ts_ry, ts_xlat, ts_xlong, ts_hgt
   REAL :: known_lat, known_lon
   CHARACTER (LEN=132) :: message
   TYPE (PROJ_INFO), INTENT(out) :: ts_proj


   INTEGER :: ids, ide, jds, jde, kds, kde,        &
              ims, ime, jms, jme, kms, kme,        &
              ips, ipe, jps, jpe, kps, kpe,        &
              imsx, imex, jmsx, jmex, kmsx, kmex,  &
              ipsx, ipex, jpsx, jpex, kpsx, kpex,  &
              imsy, imey, jmsy, jmey, kmsy, kmey,  &
              ipsy, ipey, jpsy, jpey, kpsy, kpey
   TYPE (grid_config_rec_type)               :: config_flags_temp


   config_flags_temp = config_flags

     CALL get_ijk_from_grid ( grid ,                               &
                               ids, ide, jds, jde, kds, kde,        &
                               ims, ime, jms, jme, kms, kme,        &
                               ips, ipe, jps, jpe, kps, kpe,        &
                               imsx, imex, jmsx, jmex, kmsx, kmex,  &
                               ipsx, ipex, jpsx, jpex, kpsx, kpex,  &
                               imsy, imey, jmsy, jmey, kmsy, kmey,  &
                               ipsy, ipey, jpsy, jpey, kpsy, kpey )

     CALL model_to_grid_config_rec ( grid%id , model_config_rec , config_flags_temp )


      ! Set up map transformation structure
      CALL map_init(ts_proj)


      IF (ips <= 1 .AND. 1 <= ipe .AND. &
          jps <= 1 .AND. 1 <= jpe) THEN
         known_lat = grid%xlat(1,1)
         known_lon = grid%xlong(1,1)
      ELSE
         known_lat = 9999.
         known_lon = 9999.
      END IF
      known_lat = wrf_dm_min_real(known_lat)
      known_lon = wrf_dm_min_real(known_lon)


      ! Mercator
      IF (config_flags%map_proj == PROJ_MERC) THEN
         CALL map_set(PROJ_MERC, ts_proj,               &
                      truelat1 = config_flags%truelat1, &
                      lat1     = known_lat,             &
                      lon1     = known_lon,             &
                      knowni   = 1.,                    &
                      knownj   = 1.,                    &
                      dx       = config_flags%dx)
      ! Lambert conformal
      ELSE IF (config_flags%map_proj == PROJ_LC) THEN
      CALL map_set(PROJ_LC, ts_proj,                  &
                      truelat1 = config_flags%truelat1,  &
                      truelat2 = config_flags%truelat2,  &
                      stdlon   = config_flags%stand_lon, &
                      lat1     = known_lat,              &
                      lon1     = known_lon,              &
                      knowni   = 1.,                     &
                      knownj   = 1.,                     &
                      dx       = config_flags%dx)
      ! Polar stereographic
      ELSE IF (config_flags%map_proj == PROJ_PS) THEN
         CALL map_set(PROJ_PS, ts_proj,                  &
                      truelat1 = config_flags%truelat1,  &
                      stdlon   = config_flags%stand_lon, &
                      lat1     = known_lat,              &
                      lon1     = known_lon,              &
                      knowni   = 1.,                     &
                      knownj   = 1.,                     &
                      dx       = config_flags%dx)


      ! Cassini (global ARW)
      ELSE IF (config_flags%map_proj == PROJ_CASSINI) THEN
         CALL map_set(PROJ_CASSINI, ts_proj,                            &
                      latinc   = grid%dy*360.0/(2.0*EARTH_RADIUS_M*PI), &
                      loninc   = grid%dx*360.0/(2.0*EARTH_RADIUS_M*PI), &
                      lat1     = known_lat,                             &
                      lon1     = known_lon,                             &
! We still need to get POLE_LAT and POLE_LON metadata variables before
!   this will work for rotated poles.
                      lat0     = 90.0,                                  &
                      lon0     = 0.0,                                   &
                      knowni   = 1.,                                    &
                      knownj   = 1.,                                    &
                      stdlon   = config_flags%stand_lon)


      ! Rotated latitude-longitude
      ELSE IF (config_flags%map_proj == PROJ_ROTLL) THEN
         CALL map_set(PROJ_ROTLL, ts_proj,                      &
! I have no idea how this should work for NMM nested domains
                      ixdim    = grid%e_we-1,                   &
                      jydim    = grid%e_sn-1,                   &
                      phi      = real(grid%e_sn-2)*grid%dy/2.0, &
                      lambda   = real(grid%e_we-2)*grid%dx,     &
                      lat1     = config_flags%cen_lat,          &
                      lon1     = config_flags%cen_lon,          &
                      latinc   = grid%dy,                       &
                      loninc   = grid%dx,                       &
                      stagger  = HH)


      END IF


end subroutine trajmapproj




SUBROUTINE diffusion_2ord_fm( tendency, config_flags,              &
                                     var, rdx,                           &
                                     ids, ide, jds, jde, kds, kde,        &
                                     ims, ime, jms, jme, kms, kme,        &
                                     its, ite, jts, jte, kts, kte        )
									 
									 
									 
									 
!-----------------------------------------------------------------------
! Begin declarations.

   IMPLICIT NONE

   TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

   INTEGER ,        INTENT(IN   ) ::        ids, ide, jds, jde, kds, kde, &
                                            ims, ime, jms, jme, kms, kme, &
                                            its, ite, jts, jte, kts, kte


   REAL , DIMENSION( ims:ime, kms:kme, jms:jme), INTENT(INOUT) ::tendency


   REAL , DIMENSION( ims:ime, kms:kme, jms:jme), INTENT(IN   ) ::   var 



   REAL ,                                        INTENT(IN   ) ::    rdx
! Local data

   INTEGER :: i, j, k, ktf

   INTEGER :: i_start, i_end, j_start, j_end, k_end 

   REAL :: fk 
   

! End declarations.
!-----------------------------------------------------------------------

   ktf=MIN(kte,kde-1)

!-----------------------------------------------------------------------
! u :   p (.), u(|), w(-)
!
!       p  u  p  u                                  u     u
!
! p  |  .  |  .  |  .  |   k+1                |  .  |  .  |  .  |   k+1
!
! w     - 13  -     -      k+1                     13               k+1
!
! p  |  11 O 11  |  .  |   k                  |  12 O 12  |  .  |   k
!
! w     - 13  -     -      k                       13               k
!
! p  |  .  |  .  |  .  |   k-1                |  .  |  .  |  .  |   k-1
!
!      i-1 i  i i+1                          j-1 j  j j+1 j+1
!

   i_start = its
   i_end   = ite
   j_start = jts
   j_end   = MIN(jte,jde-1)

   IF ( config_flags%open_xs .or. config_flags%specified .or. &
        config_flags%nested) i_start = MAX(ids+1,its)
   IF ( config_flags%open_xe .or. config_flags%specified .or. &
        config_flags%nested) i_end   = MIN(ide-1,ite)
   IF ( config_flags%open_ys .or. config_flags%specified .or. &
        config_flags%nested) j_start = MAX(jds+1,jts)
   IF ( config_flags%open_ye .or. config_flags%specified .or. &
        config_flags%nested) j_end   = MIN(jde-2,jte)
      IF ( config_flags%periodic_x ) i_start = its
      IF ( config_flags%periodic_x ) i_end = ite



!xkmh, rdx

fk = 75.0*rdx*rdx


  k_end = ktf - 1
 ! if(gridv == 'f' )  k_end = ktf 
  
     
   
   DO j = j_start, j_end
   DO k = 2,k_end
   DO i = i_start, i_end
   

   tendency(i, k, j) =  (var(i+1,k,j)+var(i-1,k,j)+var(i,k-1,j)+var(i,k+1,j)-var(i,k,j)*4.)*fk
   
   
   ENDDO
   
   
   
   ENDDO
   
     DO i = i_start, i_end 
         ! 下边界 
	     ! tendency(i, 1, j) =  (var(i+1,1,j)+var(i-1,1,j) +var(i,2,j)-var(i,1,j)*3.)*fk    ! 垂直1阶
		  
   	      tendency(i, 1, j) =  (var(i+1,1,j)+var(i-1,1,j)  -var(i,1,j)*2.)*fk               !  垂直2阶，下边界垂直扩散为0  
		  
     ENDDO
   
   
   ENDDO
   
   
	  
									 
END SUBROUTINE diffusion_2ord_fm	

								 
!only fore x-z 2D simulate ,and must dx=dz 
SUBROUTINE diffusion_2ord_fmh( name,tendency, config_flags,              &
                                     var, rdx,                           &
                                     ids, ide, jds, jde, kds, kde,        &
                                     ims, ime, jms, jme, kms, kme,        &
                                     its, ite, jts, jte, kts, kte        )
									 
									 
									 
									 
!-----------------------------------------------------------------------
! Begin declarations.

   IMPLICIT NONE


 CHARACTER(LEN=1) ,                          INTENT(IN   ) :: name
 
   TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

   INTEGER ,        INTENT(IN   ) ::        ids, ide, jds, jde, kds, kde, &
                                            ims, ime, jms, jme, kms, kme, &
                                            its, ite, jts, jte, kts, kte


   REAL , DIMENSION( ims:ime, kms:kme, jms:jme), INTENT(INOUT) ::tendency


   REAL , DIMENSION( ims:ime, kms:kme, jms:jme), INTENT(IN   ) ::   var 



   REAL ,                                        INTENT(IN   ) ::    rdx
! Local data

   INTEGER :: i, j, k, ktf

   INTEGER :: i_start, i_end, j_start, j_end, k_start,k_end 

   REAL :: fk 
   

! End declarations.
!-----------------------------------------------------------------------

   ktf=MIN(kte,kde-1)

!-----------------------------------------------------------------------
! u :   p (.), u(|), w(-)
!
!       p  u  p  u                                  u     u
!
! p  |  .  |  .  |  .  |   k+1                |  .  |  .  |  .  |   k+1
!
! w     - 13  -     -      k+1                     13               k+1
!
! p  |  11 O 11  |  .  |   k                  |  12 O 12  |  .  |   k
!
! w     - 13  -     -      k                       13               k
!
! p  |  .  |  .  |  .  |   k-1                |  .  |  .  |  .  |   k-1
!
!      i-1 i  i i+1                          j-1 j  j j+1 j+1
!

   i_start = its
   i_end   = ite
   j_start = jts
   j_end   = MIN(jte,jde-1)

   IF ( config_flags%open_xs .or. config_flags%specified .or. &
        config_flags%nested) i_start = MAX(ids+1,its)
   IF ( config_flags%open_xe .or. config_flags%specified .or. &
        config_flags%nested) i_end   = MIN(ide-1,ite)
   IF ( config_flags%open_ys .or. config_flags%specified .or. &
        config_flags%nested) j_start = MAX(jds+1,jts)
   IF ( config_flags%open_ye .or. config_flags%specified .or. &
        config_flags%nested) j_end   = MIN(jde-2,jte)
      IF ( config_flags%periodic_x ) i_start = its
      IF ( config_flags%periodic_x ) i_end = ite



!xkmh, rdx

fk = 75.0*rdx*rdx


  k_end = ktf - 1
   if(name == 'w' ) then 
	  k_start = 2 
   else
      k_start = 1 
   
   endif
  
     
   
   DO j = j_start, j_end
   DO k = k_start,k_end
   DO i = i_start, i_end
   

   tendency(i, k, j) =  tendency(i,k,j) + (var(i+1,k,j)+var(i-1,k,j)-var(i,k,j)*2.)*fk
   
   
   ENDDO
   
   
   
   ENDDO
   
 
   
   
   ENDDO
   
   
	  
									 
END SUBROUTINE diffusion_2ord_fmh	

!only fore x-z 2D simulate ,and must dx=dz 
SUBROUTINE diffusion_2ord_fmv( name,tendency, config_flags,              &
                                     field, rdx,                           &
                                     ids, ide, jds, jde, kds, kde,        &
                                     ims, ime, jms, jme, kms, kme,        &
                                     its, ite, jts, jte, kts, kte        )
									 
									 
									 
									 
!-----------------------------------------------------------------------
! Begin declarations.

   IMPLICIT NONE


   CHARACTER(LEN=1) ,                          INTENT(IN   ) :: name
 
   TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

   INTEGER ,        INTENT(IN   ) ::        ids, ide, jds, jde, kds, kde, &
                                            ims, ime, jms, jme, kms, kme, &
                                            its, ite, jts, jte, kts, kte


   REAL , DIMENSION( ims:ime, kms:kme, jms:jme), INTENT(INOUT) ::tendency


   REAL , DIMENSION( ims:ime, kms:kme, jms:jme), INTENT(IN   ) ::   field



   REAL ,                                        INTENT(IN   ) ::    rdx
! Local data

   INTEGER :: i, j, k, ktf

   INTEGER :: i_start, i_end, j_start, j_end, k_start,k_end 

   REAL :: fk 
   
      REAL , DIMENSION(its:ite, 0:kte+1) :: vflux

! End declarations.
!-----------------------------------------------------------------------


!xkmh, rdx

fk = 75.0*rdx*rdx


   ktf=MIN(kte,kde-1)
   
   IF (name .EQ. 'w')THEN

   
   i_start = its
   i_end   = MIN(ite,ide-1)
   j_start = jts
   j_end   = MIN(jte,jde-1)

j_loop_w : DO j = j_start, j_end

     DO k=kts,ktf-1
       DO i = i_start, i_end
	   
         
		   vflux(i,k)= 75.0*rdx*(field(i,k+1,j)-field(i,k,j))  
       ENDDO
     ENDDO

     DO i = i_start, i_end
       vflux(i,ktf)=0.
     ENDDO

     DO k=kts+1,ktf
       DO i = i_start, i_end
	   
 
            tendency(i,k,j)=tendency(i,k,j) + rdx*(vflux(i,k)-vflux(i,k-1))										 
										 
       ENDDO
     ENDDO

    ENDDO j_loop_w
     
  ELSE  

     i_start = its
     i_end   = MIN(ite,ide-1)
     j_start = jts
     j_end   = MIN(jte,jde-1)

j_loop_s : DO j = j_start, j_end

     DO k=kts,ktf-1
       DO i = i_start, i_end
				  
		 vflux(i,k) = 75.0*rdx *(field(i,k+1,j)-field(i,k,j))		  
       ENDDO
     ENDDO



     DO i = i_start, i_end
            vflux(i,0)=vflux(i,1)
     ENDDO

	IF(name .EQ. 'p') THEN
		DO i = i_start, i_end
			vflux(i,0) = 0.0 
		ENDDO
	ENDIF

     DO i = i_start, i_end
       vflux(i,ktf)=0.
     ENDDO

     DO k=kts,ktf
       DO i = i_start, i_end
	   
				
				
         tendency(i,k,j)=tendency(i,k,j) + rdx*(vflux(i,k)-vflux(i,k-1))				
				
       ENDDO
     ENDDO

 ENDDO j_loop_s

   ENDIF
   
   
	  
									 
END SUBROUTINE diffusion_2ord_fmv	

END MODULE module_em
