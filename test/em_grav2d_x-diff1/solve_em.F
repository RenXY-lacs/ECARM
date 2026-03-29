!WRF:MEDIATION_LAYER:SOLVER

SUBROUTINE solve_em ( grid , config_flags  &
! Arguments generated from Registry
#include "dummy_new_args.inc"
!
                    )
! Driver layer modules
   USE module_state_description
   USE module_domain, ONLY : &
                  domain, get_ijk_from_grid, get_ijk_from_subgrid                          &
                 ,domain_get_current_time, domain_get_start_time                           &
                 ,domain_get_sim_start_time, domain_clock_get,is_alarm_tstep
   USE module_domain_type, ONLY : history_alarm, restart_alarm, auxinput4_alarm            &
                 ,boundary_alarm
   USE module_configure, ONLY : grid_config_rec_type
   USE module_driver_constants
   USE module_machine
   USE module_tiles, ONLY : set_tiles
#ifdef DM_PARALLEL
   USE module_dm, ONLY : &
                  local_communicator, mytask, ntasks, ntasks_x, ntasks_y                   &
                 ,local_communicator_periodic, wrf_dm_maxval
   USE module_comm_dm, ONLY : &
                  halo_em_a_sub,halo_em_b_sub,halo_em_c2_sub,halo_em_chem_e_3_sub          &
                 ,halo_em_chem_e_5_sub,halo_em_chem_e_7_sub,halo_em_chem_old_e_5_sub       &
                 ,halo_em_chem_old_e_7_sub,halo_em_c_sub,halo_em_d2_3_sub                  &
                 ,halo_em_d2_5_sub,halo_em_d3_3_sub,halo_em_d3_5_sub,halo_em_d_sub         &
                 ,halo_em_e_3_sub,halo_em_e_5_sub,halo_em_hydro_uv_sub                     &
                 ,halo_em_moist_e_3_sub,halo_em_moist_e_5_sub,halo_em_moist_e_7_sub        &
                 ,halo_em_moist_old_e_5_sub,halo_em_moist_old_e_7_sub                      &
                 ,halo_em_scalar_e_3_sub,halo_em_scalar_e_5_sub,halo_em_scalar_e_7_sub     &
                 ,halo_em_scalar_old_e_5_sub,halo_em_scalar_old_e_7_sub,halo_em_tke_3_sub  &
                 ,halo_em_tke_5_sub,halo_em_tke_7_sub,halo_em_tke_advect_3_sub             &
                 ,halo_em_tke_advect_5_sub,halo_em_tke_old_e_5_sub                         &
                 ,halo_em_tke_old_e_7_sub,halo_em_tracer_e_3_sub,halo_em_tracer_e_5_sub    &
                 ,halo_em_tracer_e_7_sub,halo_em_tracer_old_e_5_sub                        &
                 ,halo_em_tracer_old_e_7_sub,halo_em_sbm_sub,period_bdy_em_a_sub                           &
                 ,period_bdy_em_b3_sub,period_bdy_em_b_sub,period_bdy_em_chem2_sub         &
                 ,period_bdy_em_chem_old_sub,period_bdy_em_chem_sub,period_bdy_em_d3_sub   &
                 ,period_bdy_em_d_sub,period_bdy_em_e_sub,period_bdy_em_moist2_sub         &
                 ,period_bdy_em_moist_old_sub,period_bdy_em_moist_sub                      &
                 ,period_bdy_em_scalar2_sub,period_bdy_em_scalar_old_sub                   &
                 ,period_bdy_em_scalar_sub,period_bdy_em_tke_old_sub, period_bdy_em_tke_sub  &
                 ,period_bdy_em_tracer2_sub,period_bdy_em_tracer_old_sub                   &
                 ,period_bdy_em_tracer_sub,period_em_da_sub,period_em_hydro_uv_sub         &
                 ,period_em_f_sub,period_em_g_sub                                          &
                 ,halo_em_f_1_sub,halo_em_init_4_sub,halo_em_thetam_sub,period_em_thetam_sub &
                 ,halo_em_d_pv_sub
#endif
   USE module_utility
! Mediation layer modules
! Model layer modules
   USE module_model_constants
   USE module_small_step_em
   USE module_em
   USE module_big_step_utilities_em
   USE module_bc
   USE module_bc_em
   USE module_solvedebug_em
   USE module_physics_addtendc
   USE module_diffusion_em
   USE module_polarfft
   USE module_microphysics_driver
   USE module_microphysics_zero_out
!   USE module_lightning_driver, ONLY : lightning_driver
   USE module_fddaobs_driver
!   USE module_diagnostics
#if (WRF_CHEM == 1)
   USE module_input_chem_data
   USE module_input_tracer
   USE module_chem_utilities
#endif
   USE module_dust_emis
   USE module_first_rk_step_part1
   USE module_first_rk_step_part2
   USE module_after_all_rk_steps
   USE module_llxy, ONLY : proj_cassini
   USE module_avgflx_em, ONLY : zero_avgflx, upd_avgflx
   USE module_cpl, ONLY : coupler_on, cpl_settime, cpl_store_input

   IMPLICIT NONE

   !  Input data.

   TYPE(domain) , TARGET          :: grid

   !  Definitions of dummy arguments to this routine (generated from Registry).
#include "dummy_new_decl.inc"

   !  Structure that contains run-time configuration (namelist) data for domain
   TYPE (grid_config_rec_type) , INTENT(IN)          :: config_flags

   ! Local data

   INTEGER                         :: k_start , k_end, its, ite, jts, jte
   INTEGER                         :: ids , ide , jds , jde , kds , kde , &
                                      ims , ime , jms , jme , kms , kme , &
                                      ips , ipe , jps , jpe , kps , kpe

   INTEGER                         :: sids , side , sjds , sjde , skds , skde , &
                                      sims , sime , sjms , sjme , skms , skme , &
                                      sips , sipe , sjps , sjpe , skps , skpe


   INTEGER ::              imsx, imex, jmsx, jmex, kmsx, kmex,    &
                           ipsx, ipex, jpsx, jpex, kpsx, kpex,    &
                           imsy, imey, jmsy, jmey, kmsy, kmey,    &
                           ipsy, ipey, jpsy, jpey, kpsy, kpey

   INTEGER                         :: ij , iteration
   INTEGER                         :: im , num_3d_m , ic , num_3d_c , is , num_3d_s
   INTEGER                         :: loop
   INTEGER                         :: sz
   INTEGER                         :: iswater

   LOGICAL                         :: specified_bdy, channel_bdy

   REAL                            :: t_new, time_duration_of_lbcs
   
   ! Changes in tendency at this timestep
   real ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33) :: h_tendency, &
                                                                                   z_tendency
                                                                                   
   ! Whether advection should produce decoupled horizontal and vertical advective tendency outputs
   LOGICAL                        :: tenddec
   
   ! Flag for producing diagnostic fields (e.g., radar reflectivity)
   LOGICAL                        :: diag_flag
   INTEGER                        :: ke_diag ! tells reflectivity calculation whether to do full depth or only k=1
   
#if (WRF_CHEM == 1)
   ! Index cross-referencing array for tendency accumulation
   INTEGER, DIMENSION( num_chem ) :: adv_ct_indices
#endif

! storage for tendencies and decoupled state (generated from Registry)

#include "i1_decl.inc"
! Previous time level of tracer arrays now defined as i1 variables;
! the state 4d arrays now redefined as 1-time level arrays in Registry.
! Benefit: save memory in nested runs, since only 1 domain is active at a
! time.  Potential problem on stack-limited architectures: increases
! amount of data on program stack by making these automatic arrays.

   INTEGER :: rc 
   INTEGER :: number_of_small_timesteps, rk_step
   INTEGER :: klevel,ijm,ijp,i,j,k,size1,size2    ! for prints/plots only
   INTEGER :: idum1, idum2, dynamics_option

   INTEGER :: rk_order, iwmax, jwmax, kwmax
   REAL :: dt_rk, dts_rk, dts, dtm, wmax
   REAL , ALLOCATABLE , DIMENSION(:)  :: max_vert_cfl_tmp, max_horiz_cfl_tmp
   LOGICAL :: leapfrog
   INTEGER :: l,kte,kk
   LOGICAL :: f_flux  ! flag for computing averaged fluxes in cu_gd
   REAL    :: curr_secs, curr_secs2
   INTEGER :: num_sound_steps
   INTEGER :: idex, jdex
   REAL    :: max_msft
   REAL    :: spacing

   INTEGER :: ii, jj !kk is above after l,kte
   REAL    :: dclat
   INTEGER :: debug_level

! urban related variables
   INTEGER :: NUM_ROOF_LAYERS, NUM_WALL_LAYERS, NUM_ROAD_LAYERS   ! urban

   TYPE(WRFU_TimeInterval)                    :: tmpTimeInterval, tmpTimeInterval2
   REAL                                       :: real_time
   LOGICAL                                    :: adapt_step_flag
   LOGICAL                                    :: fill_w_flag

! variables for flux-averaging code 20091223
   CHARACTER*256                              :: message, message2, message3
   REAL                                       :: old_dt
   TYPE(WRFU_Time)                            :: temp_time, CurrTime, restart_time
   INTEGER, PARAMETER                         :: precision = 100
   INTEGER                                    :: num, den
   TYPE(WRFU_TimeInterval)                    :: dtInterval, intervaltime,restartinterval

! Define benchmarking timers if -DBENCH is compiled
#include "bench_solve_em_def.h"

!----------------------
! Executable statements
!----------------------

!<DESCRIPTION>
!<pre>
! solve_em is the main driver for advancing a grid a single timestep.
! It is a mediation-layer routine -> DM and SM calls are made where 
! needed for parallel processing.  
!
! solve_em can integrate the equations using 3 time-integration methods
!      
!    - 3rd order Runge-Kutta time integration (recommended)
!      
!    - 2nd order Runge-Kutta time integration
!      
! The main sections of solve_em are
!     
! (1) Runge-Kutta (RK) loop
!     
! (2) Non-timesplit physics (i.e., tendencies computed for updating
!     model state variables during the first RK sub-step (loop)
!     
! (3) Small (acoustic, sound) timestep loop - within the RK sub-steps
!     
! (4) scalar advance for moist and chem scalar variables (and TKE)
!     within the RK sub-steps.
!     
! (5) time-split physics (after the RK step), currently this includes
!     only microphyics
!
! A more detailed description of these sections follows.
!</pre>
!</DESCRIPTION>

! Initialize timers if compiled with -DBENCH
#include "bench_solve_em_init.h"

!  set runge-kutta solver (2nd or 3rd order)

   dynamics_option = config_flags%rk_ord

!  Obtain dimension information stored in the grid data structure.

   CALL get_ijk_from_grid (  grid ,                   &
                             ids, ide, jds, jde, kds, kde,    &
                             ims, ime, jms, jme, kms, kme,    &
                             ips, ipe, jps, jpe, kps, kpe,    &
                             imsx, imex, jmsx, jmex, kmsx, kmex,    &
                             ipsx, ipex, jpsx, jpex, kpsx, kpex,    &
                             imsy, imey, jmsy, jmey, kmsy, kmey,    &
                             ipsy, ipey, jpsy, jpey, kpsy, kpey )
 
   CALL get_ijk_from_subgrid (  grid ,                   &
                             sids, side, sjds, sjde, skds, skde,    &
                             sims, sime, sjms, sjme, skms, skme,    &
                             sips, sipe, sjps, sjpe, skps, skpe    )
   k_start         = kps
   k_end           = kpe

   num_3d_m        = num_moist
   num_3d_c        = num_chem
   num_3d_s        = num_scalar

!  backward integration needs to advect only QV
   if (grid%dfi_stage .EQ. DFI_BCK) then
       num_3d_m    = P_QV
       num_3d_s    = PARAM_FIRST_SCALAR - 1
   endif

   f_flux = config_flags%do_avgflx_cugd .EQ. 1

!  Compute these starting and stopping locations for each tile and number of tiles.
!  See: http://www.mmm.ucar.edu/wrf/WG2/topics/settiles
   CALL set_tiles ( ZONE_SOLVE_EM, grid , ids , ide , jds , jde , ips , ipe , jps , jpe )
!   CALL set_tiles (  grid , ids , ide , jds , jde , ips , ipe , jps , jpe )

!  Max values of CFL for adaptive time step scheme

   ALLOCATE (max_vert_cfl_tmp(grid%num_tiles))
   ALLOCATE (max_horiz_cfl_tmp(grid%num_tiles))

  !
  ! Calculate current time in seconds since beginning of model run.
  !   Unfortunately, ESMF does not seem to have a way to return
  !   floating point seconds based on a TimeInterval.  So, we will
  !   calculate it here--but, this is not clean!!
  !
   tmpTimeInterval  = domain_get_current_time ( grid ) - domain_get_sim_start_time ( grid )
   tmpTimeInterval2 = domain_get_current_time ( grid ) - domain_get_start_time ( grid )
   curr_secs  = real_time(tmpTimeInterval)
   curr_secs2 = real_time(tmpTimeInterval2)

   old_dt = grid%dt   ! store old time step for flux averaging code at end of RK loop
!-----------------------------------------------------------------------------
! Adaptive time step: Added by T. Hutchinson, WSI  3/5/07
!   In this call, we do the time-step adaptation and set time-dependent lateral
!   boundary condition nudging weights.
!
   IF ( (config_flags%use_adaptive_time_step) .and. &
        ( (.not. grid%nested) .or. &
        ( (grid%nested) .and. (abs(grid%dtbc) < 0.0001) ) ) )THEN
      CALL adapt_timestep(grid, config_flags)
      adapt_step_flag = .TRUE.
   ELSE
      adapt_step_flag = .FALSE.
   ENDIF
! End of adaptive time step modifications
!-----------------------------------------------------------------------------
!
! Set diagnostic flag value history output time
!-----------------------------------------------------------------------------

   ke_diag = kms ! default to ke_diag=1 in case of nwp_diagnostics == 1
   diag_flag = .false.
   if ( Is_alarm_tstep(grid%domain_clock, grid%alarms(HISTORY_ALARM)) ) then
      diag_flag = .true.
      ke_diag = min(k_end,kde-1) ! set depth to full domain for reflectivity field
   endif
   IF (config_flags%nwp_diagnostics == 1) diag_flag = .true.

   grid%itimestep = grid%itimestep + 1
   grid%dtbc = grid%dtbc + grid%dt

   IF( coupler_on ) CALL cpl_store_input( grid, config_flags )

   IF (config_flags%polar) dclat = 90./REAL(jde-jds) !(0.5 * 180/ny)

#if (WRF_CHEM == 1)

   kte=min(k_end,kde-1)
# ifdef DM_PARALLEL
   if ( num_chem >= PARAM_FIRST_SCALAR ) then
!-----------------------------------------------------------------------
! see matching halo calls below for stencils
!--------------------------------------------------------------
     CALL wrf_debug ( 200 , ' call HALO_RK_CHEM' )
     IF      ( config_flags%h_mom_adv_order <= 4 ) THEN
#      include "HALO_EM_CHEM_E_3.inc"
       IF( config_flags%progn > 0 ) THEN
#         include "HALO_EM_SCALAR_E_3.inc"
       ENDIF
       IF( config_flags%cu_physics == CAMZMSCHEME ) THEN
#         include "HALO_EM_SCALAR_E_3.inc"
       ENDIF
     ELSE IF ( config_flags%h_mom_adv_order <= 6 ) THEN
#      include "HALO_EM_CHEM_E_5.inc"
       IF( config_flags%cu_physics == CAMZMSCHEME ) THEN
#         include "HALO_EM_SCALAR_E_5.inc"
       ENDIF
       IF( config_flags%progn > 0 ) THEN
#         include "HALO_EM_SCALAR_E_5.inc"
      ENDIF
     ELSE
       WRITE(wrf_err_message,*)'solve_em: invalid h_mom_adv_order = ',config_flags%h_mom_adv_order
       CALL wrf_error_fatal(TRIM(wrf_err_message))
     ENDIF
   ENDIF
   if ( num_tracer >= PARAM_FIRST_SCALAR ) then
!-----------------------------------------------------------------------
! see matching halo calls below for stencils
!--------------------------------------------------------------
     CALL wrf_debug ( 200 , ' call HALO_RK_tracer' )
     IF      ( config_flags%h_mom_adv_order <= 4 ) THEN
#      include "HALO_EM_TRACER_E_3.inc"
     ELSE IF ( config_flags%h_mom_adv_order <= 6 ) THEN
#      include "HALO_EM_TRACER_E_5.inc"
     ELSE
       WRITE(wrf_err_message,*)'solve_em: invalid h_mom_adv_order = ',config_flags%h_mom_adv_order
       CALL wrf_error_fatal(TRIM(wrf_err_message))
     ENDIF
   ENDIF
# endif
!--------------------------------------------------------------
   adv_ct_indices(   :  ) = 1
   IF ( config_flags%chemdiag == USECHEMDIAG ) THEN
   ! modify tendency list here
   ! note that the referencing direction here is opposite of that in chem_driver
       adv_ct_indices(p_co  ) = p_advh_co
       adv_ct_indices(p_o3  ) = p_advh_o3
       adv_ct_indices(p_no  ) = p_advh_no
       adv_ct_indices(p_no2 ) = p_advh_no2
       adv_ct_indices(p_hno3) = p_advh_hno3
       adv_ct_indices(p_iso ) = p_advh_iso
       adv_ct_indices(p_ho  ) = p_advh_ho
       adv_ct_indices(p_ho2 ) = p_advh_ho2
   END IF
#endif

   rk_order = config_flags%rk_ord

   IF ( grid%time_step_sound == 0 ) THEN
! This function will give 4 for 6*dx and 6 for 10*dx and returns even numbers only
     spacing = min(grid%dx, grid%dy)
     IF ( ( config_flags%use_adaptive_time_step ) .AND. ( config_flags%map_proj == PROJ_CASSINI ) ) THEN
       max_msft=MIN ( MAX(grid%max_msftx, grid%max_msfty) , &
                      1.0/COS(config_flags%fft_filter_lat*degrad) )
       num_sound_steps = max ( 2 * ( INT (300. * grid%dt / (spacing / max_msft) - 0.01 ) + 1 ), 4 )
     ELSE IF  ( config_flags%use_adaptive_time_step ) THEN
       max_msft= MAX(grid%max_msftx, grid%max_msfty)
       num_sound_steps = max ( 2 * ( INT (300. * grid%dt / (spacing / max_msft) - 0.01 ) + 1 ), 4 )
     ELSE
       num_sound_steps = max ( 2 * ( INT (300. * grid%dt /  spacing             - 0.01 ) + 1 ), 4 )
     END IF
     WRITE(wrf_err_message,*)'grid spacing, dt, time_step_sound=',spacing,grid%dt,num_sound_steps
     CALL wrf_debug ( 50 , wrf_err_message )
   ELSE
     num_sound_steps = grid%time_step_sound
   ENDIF

   dts = grid%dt/float(num_sound_steps)

   IF (config_flags%use_adaptive_time_step) THEN
  
     CALL get_wrf_debug_level( debug_level )
     IF ((config_flags%time_step < 0) .AND. (debug_level.GE.50)) THEN
#ifdef DM_PARALLEL
       CALL wrf_dm_maxval(grid%max_vert_cfl, idex, jdex)
#endif
       WRITE(wrf_err_message,*)'variable dt, max horiz cfl, max vert cfl: ',&
            grid%dt, grid%max_horiz_cfl, grid%max_vert_cfl
       CALL wrf_debug ( 0 , wrf_err_message )
     ENDIF

     grid%max_cfl_val = 0
     grid%max_horiz_cfl = 0
     grid%max_vert_cfl = 0
   ENDIF
   
   
   
		
		

! setting bdy tendencies to zero for DFI if constant_bc = true

     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles

!      IF( config_flags%specified .AND. grid%dfi_opt .NE. DFI_NODFI   &
!          .AND. config_flags%constant_bc .AND. (grid%dfi_stage .EQ. DFI_BCK .OR. grid%dfi_stage .EQ. DFI_FWD) ) THEN
       IF( config_flags%specified .AND. config_flags%constant_bc ) THEN

       CALL zero_bdytend (grid%u_btxs,grid%u_btxe,grid%u_btys,grid%u_btye,     &
                          grid%v_btxs,grid%v_btxe,grid%v_btys,grid%v_btye,     &
                          grid%pfm_btxs,grid%pfm_btxe,grid%pfm_btys,grid%pfm_btye, &
                          grid%rhofm_btxs,grid%rhofm_btxe,grid%rhofm_btys,grid%rhofm_btye,     &
                          grid%w_btxs,grid%w_btxe,grid%w_btys,grid%w_btye,     &
                   !       grid%mu_btxs,grid%mu_btxe,grid%mu_btys,grid%mu_btye, &
                          moist_btxs,moist_btxe,                               &
                          moist_btys,moist_btye,                               &
                          scalar_btxs,scalar_btxe,                              &
                          scalar_btys,scalar_btye,                              &
                          grid%spec_bdy_width,num_3d_m,num_3d_s,                &
                          ids,ide, jds,jde, kds,kde,                   &
                          ims,ime, jms,jme, kms,kme,                   &
                          ips,ipe, jps,jpe, kps,kpe,                   &
                          grid%i_start(ij), grid%i_end(ij),            &
                          grid%j_start(ij), grid%j_end(ij),            &
                          k_start, k_end                               )

       ENDIF

       !  If the user has requested to optionally select the moist theta (use_theta_m==1)
       !  switch, the first setting of the "old" value of theta_m uses the "old"
       !  value of Qv.  The moist_old variable does not exist until after the advection
       !  towards the end of the RK loop.  For the first time in the RK loop, we need
       !  a reasonable value for moist_old.

       CALL initialize_moist_old ( moist_old(:,:,:,P_Qv),              &
                          moist(:,:,:,P_Qv) ,                          &
                          ids,ide, jds,jde, kds,kde,                   &
                          ims,ime, jms,jme, kms,kme,                   &
                          grid%i_start(ij), grid%i_end(ij),            &
                          grid%j_start(ij), grid%j_end(ij),            &
                          k_start, k_end                               )
     ENDDO
     !$OMP END PARALLEL DO

     !  Now that we have initialized the moist_old values with P_Qv for
     !  computing a moist t_tendf after rk_step part2, fill in the halo 
     !  and period boundaries.

#if ( defined( DM_PARALLEL ) && ( ! defined( STUBMPI ) ) )
#  include "HALO_EM_MOIST_OLD_E_7.inc"
#  include "PERIOD_BDY_EM_MOIST_OLD.inc"
#endif
     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles
       im = P_Qv
       CALL set_physical_bc3d( moist_old(ims,kms,jms,im), 'p', config_flags,  &
                               ids, ide, jds, jde, kds, kde,                  &
                               ims, ime, jms, jme, kms, kme,                  &
                               ips, ipe, jps, jpe, kps, kpe,                  &
                               grid%i_start(ij), grid%i_end(ij),              &
                               grid%j_start(ij), grid%j_end(ij),              &
                               k_start    , k_end                             )
     END DO
     !$OMP END PARALLEL DO

!**********************************************************************
!
!  LET US BEGIN.......
!
!<DESCRIPTION>
!<pre>
! (1) RK integration loop is named the "Runge_Kutta_loop:"
!
!   Predictor-corrector type time integration.
!   Advection terms are evaluated at time t for the predictor step,
!   and advection is re-evaluated with the latest predicted value for
!   each succeeding time corrector step
!
!   2nd order Runge Kutta (rk_order = 2):
!   Step 1 is taken to the midpoint predictor, step 2 is the full step.
!
!   3rd order Runge Kutta (rk_order = 3):
!   Step 1 is taken to from t to dt/3, step 2 is from t to dt/2,
!   and step 3 is from t to dt.
!
!   non-timesplit physics are evaluated during first RK step and
!   these physics tendencies are stored for use in each RK pass.
!</pre>
!</DESCRIPTION>
!**********************************************************************

   Runge_Kutta_loop:  DO rk_step = 1, rk_order

   !  Set the step size and number of small timesteps for
   !  each part of the timestep
   
     write(6,*) 
     write(6,*)    '*******************    rk_step =    ********** ', rk_step   
	 write(6,*)
	
		 write(6,*) 
	 	 write(6,'(a,5(2x,5g12.5))') 'max uvwrp',  maxval(grid%u_2),maxval(grid%v_2),maxval(grid%w_2),maxval(grid%rhofm_2),maxval(grid%pfm_2)
		 write(6,'(a,5(2x,5g12.5))') 'min uvwrp',  minval(grid%u_2),minval(grid%v_2),minval(grid%w_2),minval(grid%rhofm_2),minval(grid%pfm_2) 	
		 write(6,*)
 		
     dtm = grid%dt
     IF ( rk_order == 1 ) THEN   

       write(wrf_err_message,*)' leapfrog removed, error exit for dynamics_option = ',dynamics_option
       CALL wrf_error_fatal( wrf_err_message )

     ELSE IF ( rk_order == 2 ) THEN   ! 2nd order Runge-Kutta timestep

       IF ( rk_step == 1) THEN
         dt_rk  = 0.5*grid%dt
         dts_rk = dts
         number_of_small_timesteps = num_sound_steps/2
       ELSE
         dt_rk = grid%dt
         dts_rk = dts
         number_of_small_timesteps = num_sound_steps
       ENDIF

     ELSE IF ( rk_order == 3 ) THEN ! third order Runge-Kutta

       IF ( rk_step == 1) THEN
         dt_rk = grid%dt/3.
         dts_rk = dt_rk
         number_of_small_timesteps = 1
       ELSE IF (rk_step == 2) THEN
         dt_rk  = 0.5*grid%dt
         dts_rk = dts
         number_of_small_timesteps = num_sound_steps/2
       ELSE
         dt_rk = grid%dt
         dts_rk = dts
         number_of_small_timesteps = num_sound_steps
       ENDIF

     ELSE

       write(wrf_err_message,*)' unknown solver, error exit for dynamics_option = ',dynamics_option
       CALL wrf_error_fatal( wrf_err_message )

     END IF


!
!  Time level t is in the *_2 variable in the first part 
!  of the step, and in the *_1 variable after the predictor.
!  the latest predicted values are stored in the *_2 variables.
!
     CALL wrf_debug ( 200 , ' call rk_step_prep ' )

BENCH_START(step_prep_tim)
     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )

     DO ij = 1 , grid%num_tiles
	             !   INTENT(  OUT) ::  ru,              rv,      rw,       ww,             php,      cqu,cqv,   cqw,       alt

      ! CALL rk_step_prep_fm( config_flags, grid%alt, grid%al, grid%alb,    &
      !                      ids, ide, jds, jde, kds, kde,     &
      !                      ims, ime, jms, jme, kms, kme,     &
      !                      grid%i_start(ij), grid%i_end(ij), &
       !                     grid%j_start(ij), grid%j_end(ij), &
      !                      k_start, k_end                   )
     

!       open(77,file='pfm.txt')
!       write(77,*) grid%pfm_2(grid%i_start(ij):grid%i_end(ij),k_start:k_end,2)
!       close(77)
!       !!!write(6,*) grid%i_start(ij), grid%i_end(ij),grid%j_start(ij),grid%j_end(ij),k_start,k_end
!          open(78,file='pfm-v.txt')       
!          write(78,*) grid%pfm_2(1,k_start:k_end,2)
!          close(78)
!stop

	 ! al is rhop, p is pp, alt is , rhob is 
     !!!write(6,*) 'max-pfm_base rhofm_base  = ', maxval(grid%pfm_base) ,maxval(grid%rhofm_base)

	 CALL calc_p_rho_fm( grid%al, grid%p, grid%pfm_2,grid%pfm_base,                 &
                        grid%alt, grid%rhofm_2, grid%rhofm_base, pm1,  &
                        grid%smdiv,            &
                        config_flags%non_hydrostatic, 0,            &
                        ids, ide, jds, jde, kds, kde,               &
                        ims, ime, jms, jme, kms, kme,               &
                        grid%i_start(ij), grid%i_end(ij),           &
                        grid%j_start(ij), grid%j_end(ij),           &
                        k_start    , k_end                          )
						
						
 
	


	 
 END DO
     !$OMP END PARALLEL DO
BENCH_END(step_prep_tim)

#ifdef DM_PARALLEL
!-----------------------------------------------------------------------
!  Stencils for patch communications  (WCS, 29 June 2001)
!  Note:  the small size of this halo exchange reflects the 
!         fact that we are carrying the uncoupled variables 
!         as state variables in the mass coordinate model, as
!         opposed to the coupled variables as in the height
!         coordinate model.
!
!                           * * * * *
!         *        * * *    * * * * *
!       * + *      * + *    * * + * * 
!         *        * * *    * * * * *
!                           * * * * *
!
!  3D variables - note staggering!  ru(X), rv(Y), ww(Z), php(Z)
!
!  ru     x
!  rv     x
!  ww     x
!  php    x
!  alt    x
!  ph_2   x
!  phb    x
!
!  the following are 2D (xy) variables
!
!  muu    x
!  muv    x
!  mut    x
!--------------------------------------------------------------
#    include "HALO_EM_A.inc"
#endif

! set boundary conditions on variables 
! from big_step_prep for use in big_step_proc

#ifdef DM_PARALLEL
#  include "PERIOD_BDY_EM_A.inc"
#endif

BENCH_START(set_phys_bc_tim)
     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij, ii, jj, kk )

     DO ij = 1 , grid%num_tiles

       CALL wrf_debug ( 200 , ' call rk_phys_bc_dry_1' )

!   grid%ww, ??
       CALL rk_phys_bc_dry_1_fm( config_flags, grid%u_2, grid%v_2, grid%w_2, grid%ww,      & 
                              grid%alt, grid%p,grid%msft,  &
								grid%ht,        &
                              ids, ide, jds, jde, kds, kde,      &
                              ims, ime, jms, jme, kms, kme,      &
                              ips, ipe, jps, jpe, kps, kpe,      &
                              grid%i_start(ij), grid%i_end(ij),  &
                              grid%j_start(ij), grid%j_end(ij),  &
                              k_start, k_end                )
							  
							  
							  
							  
       CALL set_physical_bc3d( grid%al, 'p', config_flags,            &
                              ids, ide, jds, jde, kds, kde,     &
                              ims, ime, jms, jme, kms, kme,     &
                              ips, ipe, jps, jpe, kps, kpe,     &
                              grid%i_start(ij), grid%i_end(ij), &
                              grid%j_start(ij), grid%j_end(ij), &
                              k_start    , k_end               )
							  


       CALL set_physical_bc3d( grid%rhofm_base, 'p', config_flags,            &      !hanshenyou   if need change to pfm_1 ? 
                               ids, ide, jds, jde, kds, kde,         &
                               ims, ime, jms, jme, kms, kme,         &
                               ips, ipe, jps, jpe, kps, kpe,         &
                               grid%i_start(ij), grid%i_end(ij),     &
                               grid%j_start(ij), grid%j_end(ij),     &
                               k_start    , k_end                    )
							   
       CALL set_physical_bc3d( grid%pfm_base, 'p', config_flags,            &      !hanshenyou   if need change to pfm_1 ? 
                               ids, ide, jds, jde, kds, kde,         &
                               ims, ime, jms, jme, kms, kme,         &
                               ips, ipe, jps, jpe, kps, kpe,         &
                               grid%i_start(ij), grid%i_end(ij),     &
                               grid%j_start(ij), grid%j_end(ij),     &
                               k_start    , k_end                    )

							   
							  
       CALL set_physical_bc3d( grid%rhofm_2, 'p', config_flags,            &
                              ids, ide, jds, jde, kds, kde,     &
                              ims, ime, jms, jme, kms, kme,     &
                              ips, ipe, jps, jpe, kps, kpe,     &
                              grid%i_start(ij), grid%i_end(ij), &
                              grid%j_start(ij), grid%j_end(ij), &
                              k_start    , k_end               )
       CALL set_physical_bc3d( grid%pfm_2, 'p', config_flags,            &
                              ids, ide, jds, jde, kds, kde,     &
                              ims, ime, jms, jme, kms, kme,     &
                              ips, ipe, jps, jpe, kps, kpe,     &
                              grid%i_start(ij), grid%i_end(ij), &
                              grid%j_start(ij), grid%j_end(ij), &
                              k_start    , k_end               )





     END DO
     !$OMP END PARALLEL DO
BENCH_END(set_phys_bc_tim)

     rk_step_is_one : IF (rk_step == 1) THEN ! only need to initialize diffusion tendencies

BENCH_START(calc_p_rho_tim)

!<DESCRIPTION>
!<pre>
!(2) The non-timesplit physics begins with a call to "phy_prep"
!    (which computes some diagnostic variables such as temperature,
!    pressure, u and v at p points, etc).  This is followed by
!    calls to the physics drivers:
!
!              radiation,
!              surface,
!              pbl,
!              cumulus,
!              fddagd,
!              3D TKE and mixing.
!<pre>
!</DESCRIPTION>

       IF (coupler_on) CALL cpl_settime( curr_secs2 )
	   
	   ! !!!!!!!write(6,*)  '===========  irst_rk_step_part1     ============  '

       CALL first_rk_step_part1 (    grid, config_flags         &
                             , moist , moist_tend               &
                             , chem  , chem_tend                &
                             , tracer, tracer_tend              &
                             , scalar , scalar_tend             &
                             , fdda3d, fdda2d                   &
                             , aerod                            &
                             , u_tendf, v_tendf               &
                             , w_tendf, pfm_tendf                &
                             , rhofm_tendf              &
                             , tke_tend                         &
                             , config_flags%use_adaptive_time_step &
                             , curr_secs                        &
                             , psim , psih , gz1oz0             &
                             , chklowq                          &
                             , cu_act_flag , hol , th_phy       &
                             , pi_phy , p_phy , grid%t_phy      &
                             , dz8w , p8w , t8w                 &
                             , ids, ide, jds, jde, kds, kde     &
                             , ims, ime, jms, jme, kms, kme     &
                             , ips, ipe, jps, jpe, kps, kpe     &
                             , imsx, imex, jmsx, jmex, kmsx, kmex    &
                             , ipsx, ipex, jpsx, jpex, kpsx, kpex    &
                             , imsy, imey, jmsy, jmey, kmsy, kmey    &
                             , ipsy, ipey, jpsy, jpey, kpsy, kpey    &
                             , k_start , k_end                  &
                             , f_flux                           &
                             , aerocu                           &
                            )
							
	   	!!!!!!!write(6,*) 'after CALL first_rk_step_part1,  rhofm_tend(i_end/2,ktf/2,1)',   rhofm_tendf(grid%i_end(ij)/2,kde/2,2)							

#ifdef DM_PARALLEL
       IF ( config_flags%bl_pbl_physics == MYNNPBLSCHEME2 .OR. &
            config_flags%bl_pbl_physics == MYNNPBLSCHEME3 ) THEN
#        include "HALO_EM_SCALAR_E_5.inc"
       ENDIF
#endif


!  hanshenyou
	if(.true.) then  ! 20220727 true for km kh  error 
       CALL first_rk_step_part2 (    grid, config_flags         &
                             , moist , moist_old , moist_tend   &
                             , chem  , chem_tend                &
                             , tracer, tracer_tend              &
                             , scalar , scalar_tend             &
                             , fdda3d, fdda2d                   &
                             , u_tendf, v_tendf               &
                             , w_tendf, pfm_tendf                &
                             , rhofm_tendf             &
                             ,  mu_tendf               &
                             , tke_tend                         &
                             , adapt_step_flag , curr_secs      &
                             , psim , psih , gz1oz0             &
                             , chklowq                          &
                             , cu_act_flag , hol , th_phy       &
                             , pi_phy , p_phy , grid%t_phy      &
                             , dz8w , p8w , t8w                 &
                             , nba_mij, num_nba_mij             & !JDM 
                             , nba_rij, num_nba_rij             & !JDM 
                             , ids, ide, jds, jde, kds, kde     &
                             , ims, ime, jms, jme, kms, kme     &
                             , ips, ipe, jps, jpe, kps, kpe     &
                             , imsx, imex, jmsx, jmex, kmsx, kmex    &
                             , ipsx, ipex, jpsx, jpex, kpsx, kpex    &
                             , imsy, imey, jmsy, jmey, kmsy, kmey    &
                             , ipsy, ipey, jpsy, jpey, kpsy, kpey    &
                             , k_start , k_end                  &
                            )

	     !write(6,*) 'after CALL first_rk_step_part2,  rhofm_tendf(grid%i_end(ij)/2,kde/2,2) =  ',   rhofm_tendf(grid%i_end(ij)/2,kde/2,2)
		 !write(6,*) 'after CALL first_rk_step_part2,    grid%xkmh(grid%i_end(ij)/2,kde/2,2) =  ',   grid%xkmh(grid%i_end(ij)/2,kde/2,2)
		 !write(6,*) 'grid%khdif  = ', grid%khdif 
		 write(6,*)  'after part2 max= ', maxval(abs(u_tendf)), maxval(abs(w_tendf)),maxval(abs(pfm_tendf)),maxval(abs(rhofm_tendf))
		 !stop
		
	endif
	
     END IF rk_step_is_one

BENCH_START(rk_tend_tim)
     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles
	 
	     !!!!!!!write(6,*) 'num_tiles : ij =  ' ,  grid%num_tiles 

       CALL wrf_debug ( 200 , ' call rk_tendency' )
	   
	! !!!write(6,*) '=============== before  rk_tendency  ================'
	



	   CALL rk_tendency ( config_flags, rk_step                                                                &
                         ,grid%u_tend, grid%v_tend, w_tend, pfm_tend, rhofm_tend                                 &  !hanshenyou
                         ,u_tendf, v_tendf, w_tendf, pfm_tendf, rhofm_tendf                                      &  !hanshenyou
                         ,grid%rthften                                                   &
                    !     ,grid%ru, grid%rv, grid%rw, grid%ww                                                   &    !
                         ,grid%u_2, grid%v_2, grid%w_2, grid%rhofm_2, grid%pfm_2                                    &   
						 ,grid%sigmdot                                                                            &   
                         ,grid%u_1, grid%v_1, grid%w_1, grid%rhofm_1, grid%pfm_1                &     
						 ,grid%rhofm_base,grid%pfm_base                                       & 
                         ,grid%h_diabatic, grid%phb, grid%t_init                                               &
                         ,grid%mu_2, grid%mut, grid%muu, grid%muv, grid%mub                                    &
                         ,grid%c1h, grid%c2h, grid%c1f, grid%c2f                                               &
                         ,grid%al, grid%alt, grid%p,                                                      &          ! hanshenyou 
						  grid%pb, grid%php, cqu, cqv, cqw                          &
                        ,grid%u_base, grid%v_base, grid%t_base, grid%qv_base, grid%z_base                     &
					 	 ,grid%msft, grid%ht,grid%zh                                          &    
						 ,grid%fzm,grid%fzp, grid%ztop                                    &    
                         ,grid%msfux,grid%msfuy, grid%msfvx, grid%msfvx_inv                                    &
                         ,grid%msfvy, grid%msftx,grid%msfty, grid%clat, grid%f, grid%e, grid%sina, grid%cosa   &
                         ,grid%fnm, grid%fnp, grid%rdn, grid%rdnw                                              &
                         ,grid%dt, grid%rdx, grid%rdy, grid%khdif, grid%kvdif, grid%xkmh, grid%xkhh            &
                         ,grid%diff_6th_opt, grid%diff_6th_factor                                              &
                         ,config_flags%momentum_adv_opt                                                        &
                         ,grid%dampcoef,grid%zdamp,config_flags%damp_opt,config_flags%rad_nudge                &
                         ,grid%cf1, grid%cf2, grid%cf3, grid%cfn, grid%cfn1, num_3d_m                          &
                         ,config_flags%non_hydrostatic, config_flags%top_lid                                   &
                      !   ,grid%u_frame, grid%v_frame                                                           &
                         ,ids, ide, jds, jde, kds, kde                                                         &
                         ,ims, ime, jms, jme, kms, kme                                                         &
                         ,grid%i_start(ij), grid%i_end(ij)                                                     &
                         ,grid%j_start(ij), grid%j_end(ij)                                                     &
                         ,k_start, k_end                                                                       &
                         ,max_vert_cfl_tmp(ij), max_horiz_cfl_tmp(ij)                                         )
						 
						 
			 
	
				 write(6,*)  'after rk_tendency tend = ', maxval(abs(grid%u_tend)), maxval(abs(w_tend)),maxval(abs(pfm_tend)),maxval(abs(rhofm_tend))			
				 write(6,*)  'after rk_tendency  max tendf = ', maxval(abs(u_tendf)), maxval(abs(w_tendf)),maxval(abs(pfm_tendf)),maxval(abs(rhofm_tendf))
		         write(6,*)  'after rk_tendency  maxloc grid%u_tend = ', maxloc(abs(grid%u_tend)) 
     END DO
     !$OMP END PARALLEL DO
BENCH_END(rk_tend_tim)

     IF (config_flags%use_adaptive_time_step) THEN
       DO ij = 1 , grid%num_tiles
         IF (max_horiz_cfl_tmp(ij) .GT. grid%max_horiz_cfl) THEN
           grid%max_horiz_cfl = max_horiz_cfl_tmp(ij)
         ENDIF
         IF (max_vert_cfl_tmp(ij) .GT. grid%max_vert_cfl) THEN
           grid%max_vert_cfl = max_vert_cfl_tmp(ij)
         ENDIF
       END DO
     
       IF (grid%max_horiz_cfl .GT. grid%max_cfl_val) THEN
         grid%max_cfl_val = grid%max_horiz_cfl
       ENDIF
       IF (grid%max_vert_cfl .GT. grid%max_cfl_val) THEN
         grid%max_cfl_val = grid%max_vert_cfl
       ENDIF
     ENDIF

BENCH_START(relax_bdy_dry_tim)
     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles

       IF ( (config_flags%specified .or. config_flags%nested) .and. ( rk_step == 1 ) ) THEN 

         CALL relax_bdy_dry ( config_flags,                                &
                              u_tendf, v_tendf, pfm_tendf, rhofm_tendf,                  &             ! hanshenyou 
                              w_save, mu_tend, grid%c1h, grid%c2h, grid%c1f, grid%c2f, &                     
                              grid%u_2, grid%v_2, grid%pfm_2, grid%rhofm_2,                           &    ! hanshenyou 
                              grid%w_2, grid%mu_2, grid%mut,                              &
                              grid%u_bxs,grid%u_bxe,grid%u_bys,grid%u_bye, &
                              grid%v_bxs,grid%v_bxe,grid%v_bys,grid%v_bye, &
                              grid%pfm_bxs,grid%pfm_bxe,grid%pfm_bys,grid%pfm_bye, &                  ! hanshenyou 
                              grid%rhofm_bxs,grid%rhofm_bxe,grid%rhofm_bys,grid%rhofm_bye, &          ! hanshenyou 
                              grid%w_bxs,grid%w_bxe,grid%w_bys,grid%w_bye, &
                              grid%mu_bxs,grid%mu_bxe,grid%mu_bys,grid%mu_bye, &
                              grid%u_btxs,grid%u_btxe,grid%u_btys,grid%u_btye, &
                              grid%v_btxs,grid%v_btxe,grid%v_btys,grid%v_btye, &
                              grid%pfm_btxs,grid%pfm_btxe,grid%pfm_btys,grid%pfm_btye,  &                ! hanshenyou 
                              grid%rhofm_btxs,grid%rhofm_btxe,grid%rhofm_btys,grid%rhofm_btye, &       ! hanshenyou     
                              grid%w_btxs,grid%w_btxe,grid%w_btys,grid%w_btye, &
                              grid%mu_btxs,grid%mu_btxe,grid%mu_btys,grid%mu_btye, &
                              config_flags%spec_bdy_width, grid%spec_zone, grid%relax_zone,       &
                              grid%dtbc, grid%fcx, grid%gcx,                              &
                              ids,ide, jds,jde, kds,kde,                   &
                              ims,ime, jms,jme, kms,kme,                   &
                              ips,ipe, jps,jpe, kps,kpe,                   &
                              grid%i_start(ij), grid%i_end(ij),            &
                              grid%j_start(ij), grid%j_end(ij),            &
                              k_start, k_end                              )

       ENDIF
	   
	  !  CALL rk_addtend_dry 
	  
	  		!!!!!!!write(6,*) 'before CALL rk_addtend_dry_fm,  rhofm_tend(i_end/2,ktf/2,1)',   rhofm_tend(grid%i_end(ij)/2,kde/2,2)	
		  		!!!!!!!write(6,*) 'before CALL rk_addtend_dry_fm,  rhofm_tendf(i_end/2,ktf/2,1)',   rhofm_tendf(grid%i_end(ij)/2,kde/2,2)		


        if(.true.) then 
       CALL rk_addtend_dry_fm( grid%u_tend,  grid%v_tend,  w_tend,  pfm_tend,  rhofm_tend,  &
                            u_tendf, v_tendf, w_tendf, pfm_tendf, rhofm_tendf, &
                            grid%u_save, grid%v_save, w_save, ph_save, grid%t_save, &
                            mu_tend, mu_tendf, rk_step,                      &
                            grid%c1h, grid%c2h,                              &
                            grid%h_diabatic, grid%mut, grid%msftx,           &
                            grid%msfty, grid%msfux,grid%msfuy,               &
                            grid%msfvx, grid%msfvx_inv, grid%msfvy,          &
                            ids,ide, jds,jde, kds,kde,                       &
                            ims,ime, jms,jme, kms,kme,                       &
                            ips,ipe, jps,jpe, kps,kpe,                       &
                            grid%i_start(ij), grid%i_end(ij),                &
                            grid%j_start(ij), grid%j_end(ij),                &
                            k_start, k_end                                  )
							
							
				 write(6,*)  'after rk_addtend_dry_fm max tend = ', maxval(abs(grid%u_tend)), maxval(abs(w_tend)),maxval(abs(pfm_tend)),maxval(abs(rhofm_tend))
		 !stop
							
							
		endif
		
							
							
		!!!!!!!write(6,*) 'after CALL rk_addtend_dry_fm,  rhofm_tend(i_end/2,ktf/2,1)',   rhofm_tend(grid%i_end(ij)/2,kde/2,2)		 			

       IF( config_flags%specified .or. config_flags%nested ) THEN 
         CALL spec_bdy_dry ( config_flags,                                    &
                             grid%u_tend, grid%v_tend, pfm_tend, rhofm_tend,               &
                             w_tend, mu_tend,                                &
                             grid%u_bxs,grid%u_bxe,grid%u_bys,grid%u_bye, &
                             grid%v_bxs,grid%v_bxe,grid%v_bys,grid%v_bye, &
                             grid%ph_bxs,grid%ph_bxe,grid%ph_bys,grid%ph_bye, &
                             grid%t_bxs,grid%t_bxe,grid%t_bys,grid%t_bye, &
                             grid%w_bxs,grid%w_bxe,grid%w_bys,grid%w_bye, &
                             grid%mu_bxs,grid%mu_bxe,grid%mu_bys,grid%mu_bye, &
                             grid%u_btxs,grid%u_btxe,grid%u_btys,grid%u_btye, &
                             grid%v_btxs,grid%v_btxe,grid%v_btys,grid%v_btye, &
                             grid%pfm_btxs,grid%pfm_btxe,grid%pfm_btys,grid%pfm_btye, &   ! hanshenyou 
                             grid%rhofm_btxs,grid%rhofm_btxe,grid%rhofm_btys,grid%rhofm_btye, &         ! hanshenyou
                             grid%w_btxs,grid%w_btxe,grid%w_btys,grid%w_btye, &
                             grid%mu_btxs,grid%mu_btxe,grid%mu_btys,grid%mu_btye, &
                             config_flags%spec_bdy_width, grid%spec_zone,                       &
                             ids,ide, jds,jde, kds,kde,  & ! domain dims
                             ims,ime, jms,jme, kms,kme,  & ! memory dims
                             ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                             grid%i_start(ij), grid%i_end(ij),                &
                             grid%j_start(ij), grid%j_end(ij),                &
                             k_start, k_end                                  )
							 
							 
			   	!!!!!!!write(6,*) 'after CALL spec_bdy_dry,  rhofm_tend(i_end/2,ktf/2,1)',   rhofm_tend(grid%i_end(ij)/2,kde/2,2)						 


       ENDIF

!---------------------------------------------------------------------------------------------
! KRS 9/12/2012: Inserted new IF block calls to spec_bdy_dry_perturb. If peturb_bdy=1, SKEBS 
! pattern passed in for perturbing the specified boundry conditions.  If peturb_bdy=2, user
! must provide pattern.  mu_2, mub, msf* also passed in for coupling needed for tendecies.
!---------------------------------------------------------------------------------------------

     if(.false.) then ! hanshenyou 
       IF( config_flags%specified .and. config_flags%perturb_bdy==1 ) THEN 
         CALL spec_bdy_dry_perturb ( config_flags,                                 &
                             grid%u_tend, grid%v_tend, rhofm_tend,                   &   !hanshenyou
                             grid%mu_2, grid%mub, grid%c1h, grid%c2h,              &
                             grid%msfux, grid%msfvx, grid%msft,	                   &
                             grid%ru_tendf_stoch, grid%rv_tendf_stoch, grid%rt_tendf_stoch, &
                             config_flags%spec_bdy_width, grid%spec_zone,                   &
                             grid%num_stoch_levels,      & ! stoch dims
                             ids,ide, jds,jde, kds,kde,  & ! domain dims
                             ims,ime, jms,jme, kms,kme,  & ! memory dims
                             ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                             grid%i_start(ij), grid%i_end(ij),                &
                             grid%j_start(ij), grid%j_end(ij),                &
                             k_start, k_end                                  )
     
       ENDIF

       IF( config_flags%specified .and. config_flags%perturb_bdy==2 ) THEN
         CALL spec_bdy_dry_perturb ( config_flags,                                 &
                             grid%u_tend, grid%v_tend, rhofm_tend,                   &   ! hanshenyou
                             grid%mu_2, grid%mub, grid%c1h, grid%c2h,              &
                             grid%msfux, grid%msfvx, grid%msft,                    &
                             grid%field_u_tend_perturb, grid%field_v_tend_perturb, grid%field_t_tend_perturb, &
                             config_flags%spec_bdy_width, grid%spec_zone,                   &
                             grid%num_stoch_levels,      & ! stoch dims
                             ids,ide, jds,jde, kds,kde,  & ! domain dims
                             ims,ime, jms,jme, kms,kme,  & ! memory dims
                             ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                             grid%i_start(ij), grid%i_end(ij),                &
                             grid%j_start(ij), grid%j_end(ij),                &
                             k_start, k_end                                  )
  
       ENDIF
	   
	  endif ! hanshenyou 

     END DO
     !$OMP END PARALLEL DO
BENCH_END(relax_bdy_dry_tim)

!<DESCRIPTION>
!<pre>
! (3) Small (acoustic,sound) steps.
!
!    Several acoustic steps are taken each RK pass.  A small step 
!    sequence begins with calculating perturbation variables 
!    and coupling them to the column dry-air-mass mu 
!    (call to small_step_prep).  This is followed by computing
!    coefficients for the vertically implicit part of the
!    small timestep (call to calc_coef_w).  
!
!    The small steps are taken
!    in the named loop "small_steps:".  In the small_steps loop, first 
!    the horizontal momentum (u and v) are advanced (call to advance_uv),
!    next mu and theta are advanced (call to advance_mu_t) followed by
!    advancing w and the geopotential (call to advance_w).  Diagnostic
!    values for pressure and inverse density are updated at the end of
!    each small_step.
!
!    The small-step section ends with the change of the perturbation variables
!    back to full variables (call to small_step_finish).
!</pre>
!</DESCRIPTION>

BENCH_START(small_step_prep_tim)
     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles

    ! Calculate coefficients for the vertically implicit acoustic/gravity wave
    ! integration.  We only need calculate these for the first pass through -
    ! the predictor step.  They are reused as is for the corrector step.
    ! For third-order RK, we need to recompute these after the first 
    ! predictor because we may have changed the small timestep -> grid%dts.

       CALL wrf_debug ( 200 , ' call small_step_prep ' )
	   
	!!!!write(6,*)  '===========  small_step_prep_fm     ============  '
	   
       CALL small_step_prep_fm( grid%u_1,grid%u_2,grid%v_1,grid%v_2,grid%w_1,grid%w_2,   &
                             grid%rhofm_1,grid%rhofm_2,grid%pfm_1,grid%pfm_2,                   &   ! hanshenyou
                             grid%u_save, grid%v_save, w_save,                        &
						     rhofm_save, pfm_save,                          &
                             rk_step,                             &
                             ids, ide, jds, jde, kds, kde,                            &
                             ims, ime, jms, jme, kms, kme,                            &
                             grid%i_start(ij), grid%i_end(ij),                        &
                             grid%j_start(ij), grid%j_end(ij),                        &
                             k_start    , k_end                                       )

	  
     
         write(6,*) '   ' 
	 	 write(6,'(a,i,5(2x,5f12.5))') 'after small_step_prep_fm max uvwrp',  rk_step,maxval(grid%u_2),maxval(grid%v_2),maxval(grid%w_2),maxval(grid%rhofm_2),maxval(grid%pfm_2)
		 write(6,'(a,i,5(2x,5f12.5))') 'after small_step_prep_fm min uvwrp',  rk_step,minval(grid%u_2),minval(grid%v_2),minval(grid%w_2),minval(grid%rhofm_2),minval(grid%pfm_2) 
         write(6,*) ' max pfm_2 :   ' 
		 write(6,*) maxval(abs(	grid%pfm_2) ) , maxloc(abs(	grid%pfm_2) )
		 write(6,*) maxval(abs(	grid%pfm_2(grid%i_start(ij):grid%i_end(ij),k_start:k_end,grid%j_start(ij):grid%j_end(ij))) ) ,  &
		            maxloc(abs(	grid%pfm_2(grid%i_start(ij):grid%i_end(ij),k_start:k_end,grid%j_start(ij):grid%j_end(ij))) )
						

       IF (config_flags%non_hydrostatic) THEN
         CALL calc_coef_w( a,alpha,gamma,                    &
                           grid%mut,                         &
                           grid%c1h, grid%c2h, grid%c1f, grid%c2f, &
                           grid%c3h, grid%c4h, grid%c3f, grid%c4f, &
                           cqw, grid%rdn, grid%rdnw, c2a,    &
                           dts_rk, g, grid%epssm,            &
                           config_flags%top_lid,             &
                           ids, ide, jds, jde, kds, kde,     &
                           ims, ime, jms, jme, kms, kme,     &
                           grid%i_start(ij), grid%i_end(ij), &
                           grid%j_start(ij), grid%j_end(ij), &
                           k_start    , k_end               )
       ENDIF

     ENDDO
     !$OMP END PARALLEL DO
BENCH_END(small_step_prep_tim)

#ifdef DM_PARALLEL
!-----------------------------------------------------------------------
!  Stencils for patch communications  (WCS, 29 June 2001)
!  Note:  the small size of this halo exchange reflects the 
!         fact that we are carrying the uncoupled variables 
!         as state variables in the mass coordinate model, as
!         opposed to the coupled variables as in the height
!         coordinate model.
!
!                              * * * * *
!            *        * * *    * * * * *
!          * + *      * + *    * * + * * 
!            *        * * *    * * * * *
!                              * * * * *
!
!  3D variables - note staggering!  ph_2(Z), u_save(X), v_save(Y)
!
!  ph_2      x
!  al        x
!  p         x
!  t_1       x
!  t_save    x
!  u_save    x
!  v_save    x
!
!  the following are 2D (xy) variables
!
!  mu_1      x
!  mu_2      x
!  mudf      x
!  php       x
!  alt       x
!  pb        x
!--------------------------------------------------------------
#      include "HALO_EM_B.inc"
#      include "PERIOD_BDY_EM_B.inc"
#endif

BENCH_START(set_phys_bc2_tim)
     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )

     DO ij = 1 , grid%num_tiles

       CALL set_physical_bc3d( grid%u_tend, 'u', config_flags,      &
                               ids, ide, jds, jde, kds, kde,         &
                               ims, ime, jms, jme, kms, kme,         &
                               ips, ipe, jps, jpe, kps, kpe,         &
                               grid%i_start(ij), grid%i_end(ij),     &
                               grid%j_start(ij), grid%j_end(ij),     &
                               k_start    , k_end                    )

       CALL set_physical_bc3d( grid%v_tend, 'v', config_flags,      &
                               ids, ide, jds, jde, kds, kde,         &
                               ims, ime, jms, jme, kms, kme,         &
                               ips, ipe, jps, jpe, kps, kpe,         &
                               grid%i_start(ij), grid%i_end(ij),     &
                               grid%j_start(ij), grid%j_end(ij),     &
                               k_start    , k_end                    )
							   
       CALL set_physical_bc3d(  rhofm_tend, 'p', config_flags,           &
                               ids, ide, jds, jde, kds, kde,         &
                               ims, ime, jms, jme, kms, kme,         &
                               ips, ipe, jps, jpe, kps, kpe,         &
                               grid%i_start(ij), grid%i_end(ij),     &
                               grid%j_start(ij), grid%j_end(ij),     &
                               k_start    , k_end                    )							   
							   
       CALL set_physical_bc3d(  pfm_tend, 'p', config_flags,           &
                               ids, ide, jds, jde, kds, kde,         &
                               ims, ime, jms, jme, kms, kme,         &
                               ips, ipe, jps, jpe, kps, kpe,         &
                               grid%i_start(ij), grid%i_end(ij),     &
                               grid%j_start(ij), grid%j_end(ij),     &
                               k_start    , k_end                    )
							   


				   



       CALL set_physical_bc3d( grid%rhofm_2, 'p', config_flags,          &    !hanshenyou
                               ids, ide, jds, jde, kds, kde,         &
                               ims, ime, jms, jme, kms, kme,         &
                               ips, ipe, jps, jpe, kps, kpe,         &
                               grid%i_start(ij), grid%i_end(ij),     &
                               grid%j_start(ij), grid%j_end(ij),     &
                               k_start    , k_end                    )
							   
							   
							   
       CALL set_physical_bc3d( grid%pfm_2, 'p', config_flags,          &    !hanshenyou
                               ids, ide, jds, jde, kds, kde,         &
                               ims, ime, jms, jme, kms, kme,         &
                               ips, ipe, jps, jpe, kps, kpe,         &
                               grid%i_start(ij), grid%i_end(ij),     &
                               grid%j_start(ij), grid%j_end(ij),     &
                               k_start    , k_end                    )
							   
							   
	    write(6,*)  'after set_physical_bc3d  grid%pfm_2  max min =  ', maxval(grid%pfm_2),minval(grid%pfm_2)						   
							   
				 write(6,*) maxval(abs(	grid%pfm_2(grid%i_start(ij):grid%i_end(ij),k_start:k_end,grid%j_start(ij):grid%j_end(ij))) ) ,  &
		            maxloc(abs(	grid%pfm_2(grid%i_start(ij):grid%i_end(ij),k_start:k_end,grid%j_start(ij):grid%j_end(ij))) )

					
       CALL set_physical_bc3d( pfm_save, 'p', config_flags,       &
                               ids, ide, jds, jde, kds, kde,         &
                               ims, ime, jms, jme, kms, kme,         &
                               ips, ipe, jps, jpe, kps, kpe,         &
                               grid%i_start(ij), grid%i_end(ij),     &
                               grid%j_start(ij), grid%j_end(ij),     &
                               k_start    , k_end                    )	

							   
     write(6,*)  'after set_physical_bc3d  pfm_save max min =  ', maxval(pfm_save),minval(pfm_save)	
	 
	 		 write(6,*) maxval(abs(	pfm_save(grid%i_start(ij):grid%i_end(ij),k_start:k_end,grid%j_start(ij):grid%j_end(ij))) ) ,  &
		            maxloc(abs(	pfm_save(grid%i_start(ij):grid%i_end(ij),k_start:k_end,grid%j_start(ij):grid%j_end(ij))) )
					

       CALL set_physical_bc3d( rhofm_save, 'p', config_flags,       &
                               ids, ide, jds, jde, kds, kde,         &
                               ims, ime, jms, jme, kms, kme,         &
                               ips, ipe, jps, jpe, kps, kpe,         &
                               grid%i_start(ij), grid%i_end(ij),     &
                               grid%j_start(ij), grid%j_end(ij),     &
                               k_start    , k_end                    )
							   
							   
			   
							   
							   

     END DO
     !$OMP END PARALLEL DO
BENCH_END(set_phys_bc2_tim)



       !$OMP PARALLEL DO   &
       !$OMP PRIVATE ( ij )
	   
       DO ij = 1 , grid%num_tiles

BENCH_START(calcu_abc_z)

         CALL calcu_abc_z( grid%u_2,grid%v_2,grid%w_2, grid%pfm_2, grid%rhofm_2,w_tend,  pfm_tend,  rhofm_tend,     &
                           grid%alt,grid%al,rhofm_save,pfm_save,grid%rhofm_base,                      &
						!   grid%a_fm,grid%b_fm,grid%c_fm,grid%wa_fm,grid%wb_fm,grid%ra,grid%rb,grid%pa,grid%pb,grid%pc,grid%pd,grid%pe,      
						   a_fm,b_fm,c_fm,wa_fm,wb_fm,ra,rb,pa,pb,pc,pd,pe,      &
						   config_flags,grid%spec_zone,                               &
                           grid%msft,grid%ht,                                      &
						   grid%zh,grid%fzm,grid%fzp,                              &      
                           grid%rdx, grid%rdy, dts_rk,grid%ztop,                   &				   
                           ids, ide, jds, jde, kds, kde,                          &
                           ims, ime, jms, jme, kms, kme,                          &
                           grid%i_start(ij), grid%i_end(ij),                      &
                           grid%j_start(ij), grid%j_end(ij),                      &
                           k_start    , k_end                                     )
						   

	 		 write(6,*) 'max pfm_save', maxval(abs(	pfm_save(grid%i_start(ij):grid%i_end(ij),k_start:k_end,grid%j_start(ij):grid%j_end(ij))) ) ,  &
		            maxloc(abs(	pfm_save(grid%i_start(ij):grid%i_end(ij),k_start:k_end,grid%j_start(ij):grid%j_end(ij))) )
					
						   
BENCH_END(calcu_abc_z)

       END DO
       !$OMP END PARALLEL DO
	   
	   


     small_steps : DO iteration = 1 , number_of_small_timesteps








	  !write(6,*) 
	  !write(6,*)  '******************  small_steps  iteration =  ',  iteration,number_of_small_timesteps 
	  !write(6,*) 
	  
       ! Boundary condition time (or communication time).  
#ifdef DM_PARALLEL
#      include "PERIOD_BDY_EM_B.inc"
#endif

       !$OMP PARALLEL DO   &
       !$OMP PRIVATE ( ij )

       DO ij = 1 , grid%num_tiles

BENCH_START(advance_uv_tim)

	   !!!!!!!write(6,*)  '===========  advance_uv_fm    ============  '	




		 
         CALL advance_uv_fm ( grid%u_2, grid%u_tend, grid%v_2, grid%v_tend,        &
                           grid%pfm_2,grid%alt,                                               &
						   config_flags,grid%spec_zone,                               &
                           grid%msft,grid%ht,                                      &
						   grid%zh,grid%fzm,grid%fzp,                              &      
                           grid%rdx, grid%rdy, dts_rk,grid%ztop,                   &				   
                           ids, ide, jds, jde, kds, kde,                          &
                           ims, ime, jms, jme, kms, kme,                          &
                           grid%i_start(ij), grid%i_end(ij),                      &
                           grid%j_start(ij), grid%j_end(ij),                      &
                           k_start    , k_end                                     )
						   
						   
			
			!write(6,*) 'grid%u_2(1,2,1) =  ', grid%u_2(1,2,1) 
				
           ! fore v = 0  template used 
            
           !!!write(6,*)  'after advance_uv_fm max min grid%u_2 =  ', maxval(grid%u_2),minval(grid%u_2)
		   grid%v_2 = 0.0 
	  
BENCH_END(advance_uv_tim)

       END DO
       !$OMP END PARALLEL DO





! hanshenyou 202112 22

BENCH_START(phys_bc_uv)
       !$OMP PARALLEL DO   &
       !$OMP PRIVATE ( ij )
       DO ij = 1 , grid%num_tiles
         CALL set_physical_bc3d( grid%u_2, 'u', config_flags,   &
                                 ids, ide, jds, jde, kds, kde,      &
                                 ims, ime, jms, jme, kms, kme,      &
                                 ips, ipe, jps, jpe, kps, kpe,      &
                                 grid%i_start(ij), grid%i_end(ij),  &
                                 grid%j_start(ij), grid%j_end(ij),  &
                                 k_start    , k_end                )

         CALL set_physical_bc3d( grid%v_2, 'v', config_flags,   &
                                 ids, ide, jds, jde, kds, kde,      &
                                 ims, ime, jms, jme, kms, kme,      &
                                 ips, ipe, jps, jpe, kps, kpe,      &
                                 grid%i_start(ij), grid%i_end(ij),  &
                                 grid%j_start(ij), grid%j_end(ij),  &
                                 k_start    , k_end                )
       END DO
       !$OMP END PARALLEL DO
BENCH_END(phys_bc_uv)
! hanshenyou 



      !!!write(6,*)  'after set_physical_bc3d max min grid%u_2 =  ', maxval(grid%u_2),minval(grid%u_2)













       !$OMP PARALLEL DO   &
       !$OMP PRIVATE ( ij )
       DO ij = 1 , grid%num_tiles

BENCH_START(spec_bdy_uv_tim)
         IF( config_flags%specified .or. config_flags%nested ) THEN
           CALL spec_bdyupdate(grid%u_2, grid%u_tend, dts_rk,      &
                               'u'         , config_flags, &
                                grid%spec_zone,                  &
                                ids,ide, jds,jde, kds,kde,  & ! domain dims
                                ims,ime, jms,jme, kms,kme,  & ! memory dims
                                ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                                grid%i_start(ij), grid%i_end(ij),         &
                                grid%j_start(ij), grid%j_end(ij),         &
                                k_start    , k_end             )

           CALL spec_bdyupdate(grid%v_2, grid%v_tend, dts_rk,      &
                                'v'         , config_flags, &
                                grid%spec_zone,                  &
                                ids,ide, jds,jde, kds,kde,  & ! domain dims
                                ims,ime, jms,jme, kms,kme,  & ! memory dims
                                ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                                grid%i_start(ij), grid%i_end(ij),         &
                                grid%j_start(ij), grid%j_end(ij),         &
                                k_start    , k_end             )

         ENDIF
BENCH_END(spec_bdy_uv_tim)

       END DO
       !$OMP END PARALLEL DO

#ifdef DM_PARALLEL
!
!  Stencils for patch communications  (WCS, 29 June 2001)
!
!         *                     *
!       * + *      * + *        +
!         *                     *
!
!  u_2               x
!  v_2                          x
!
#     include "HALO_EM_C.inc"
#endif

       !$OMP PARALLEL DO   &
       !$OMP PRIVATE ( ij )
       DO ij = 1 , grid%num_tiles

        !  advance the mass in the column, theta, and calculate ww
		
BENCH_START(advance_w_tim)

         write(6,*) '   ' 
	 	 write(6,'(a,i,5(2x,5f12.5))') 'before advance_w_p_rho_fm max uvwrp',  rk_step,maxval(grid%u_2),maxval(grid%v_2),maxval(grid%w_2),maxval(grid%rhofm_2),maxval(grid%pfm_2)
		 write(6,'(a,i,5(2x,5f12.5))') 'before advance_w_p_rho_fm min uvwrp',  rk_step,minval(grid%u_2),minval(grid%v_2),minval(grid%w_2),minval(grid%rhofm_2),minval(grid%pfm_2) 
         write(6,*) '   ' 
				


			   
         CALL advance_w_p_rho_fm_abc( grid%u_2,grid%v_2,grid%w_2, grid%pfm_2, grid%rhofm_2,w_tend,  pfm_tend,  rhofm_tend,     &
                           grid%alt,grid%al,rhofm_save,pfm_save,grid%rhofm_base,                      &
						!   grid%a_fm,grid%b_fm,grid%c_fm,grid%wa_fm,grid%wb_fm,grid%ra,grid%rb,grid%pa,grid%pb,grid%pc,grid%pd,grid%pe,      &
						      a_fm,b_fm,c_fm,wa_fm,wb_fm,ra,rb,pa,pb,pc,pd,pe,      &
						   config_flags,grid%spec_zone,                               &
                           grid%msft,grid%ht,                                      &
						   grid%zh,grid%fzm,grid%fzp,                              &      
                           grid%rdx, grid%rdy, dts_rk,grid%ztop,                   &				   
                           ids, ide, jds, jde, kds, kde,                          &
                           ims, ime, jms, jme, kms, kme,                          &
                           grid%i_start(ij), grid%i_end(ij),                      &
                           grid%j_start(ij), grid%j_end(ij),                      &
                           k_start    , k_end                                     )						   




         write(6,*) '   ' 
	 	 write(6,'(a,i,5(2x,5f12.5))') 'after advance_w_p_rho_fm max uvwrp',  rk_step,maxval(grid%u_2),maxval(grid%v_2),maxval(grid%w_2),maxval(grid%rhofm_2),maxval(grid%pfm_2)
		 write(6,'(a,i,5(2x,5f12.5))') 'after advance_w_p_rho_fm min uvwrp',  rk_step,minval(grid%u_2),minval(grid%v_2),minval(grid%w_2),minval(grid%rhofm_2),minval(grid%pfm_2) 
         write(6,*) '   ' 
				

        	 		 write(6,*) 'max pfm_save', maxval(abs(	pfm_save(grid%i_start(ij):grid%i_end(ij),k_start:k_end,grid%j_start(ij):grid%j_end(ij))) ) ,  &
		                                        maxloc(abs(	pfm_save(grid%i_start(ij):grid%i_end(ij),k_start:k_end,grid%j_start(ij):grid%j_end(ij))) )
					


BENCH_END(advance_w_tim)		
		


         CALL set_physical_bc3d( grid%w_2, 'w', config_flags,   &
                                 ids, ide, jds, jde, kds, kde,      &
                                 ims, ime, jms, jme, kms, kme,      &
                                 ips, ipe, jps, jpe, kps, kpe,      &
                                 grid%i_start(ij), grid%i_end(ij),  &
                                 grid%j_start(ij), grid%j_end(ij),  &
                                 k_start    , k_end                )







         CALL set_physical_bc3d( grid%rhofm_2, 'p', config_flags,             &
                                 ids, ide, jds, jde, kds, kde,     &
                                 ims, ime, jms, jme, kms, kme,     &
                                 ips, ipe, jps, jpe, kps, kpe,     &
                                 grid%i_start(ij), grid%i_end(ij), &
                                 grid%j_start(ij), grid%j_end(ij), &
                                 k_start    , k_end               )	  


       CALL set_physical_bc3d( grid%pfm_2, 'p', config_flags,            &
                              ids, ide, jds, jde, kds, kde,     &
                              ims, ime, jms, jme, kms, kme,     &
                              ips, ipe, jps, jpe, kps, kpe,     &
                              grid%i_start(ij), grid%i_end(ij), &
                              grid%j_start(ij), grid%j_end(ij), &
                              k_start    , k_end               )
							  
       !!!write(6,*)  'after  advance_w_p_rho_fm max min grid%rhofm_2 =  ', maxval(grid%rhofm_2),minval(grid%rhofm_2)
		
			    write(6,*)  'after set_physical_bc3d  grid%pfm_2  max min =  ', maxval(grid%pfm_2),minval(grid%pfm_2)						   
							   
				 write(6,*) maxval(abs(	grid%pfm_2(grid%i_start(ij):grid%i_end(ij),k_start:k_end,grid%j_start(ij):grid%j_end(ij))) ) ,  &
		            maxloc(abs(	grid%pfm_2(grid%i_start(ij):grid%i_end(ij),k_start:k_end,grid%j_start(ij):grid%j_end(ij))) )
				  

       ENDDO
       !$OMP END PARALLEL DO



BENCH_START(spec_bdy_t_tim)

       !$OMP PARALLEL DO   &
       !$OMP PRIVATE ( ij )
       DO ij = 1 , grid%num_tiles

         IF( config_flags%specified .or. config_flags%nested ) THEN

           CALL spec_bdyupdate(grid%rhofm_2, rhofm_tend, dts_rk,        &
                               't'         , config_flags,      &
                               grid%spec_zone,                  &
                               ids,ide, jds,jde, kds,kde,       &
                               ims,ime, jms,jme, kms,kme,       &
                               ips,ipe, jps,jpe, kps,kpe,       &
                               grid%i_start(ij), grid%i_end(ij),&
                               grid%j_start(ij), grid%j_end(ij),&
                               k_start    , k_end              )
							   
							   
           CALL spec_bdyupdate(grid%pfm_2, pfm_tend, dts_rk,        &
                               't'         , config_flags,      &
                               grid%spec_zone,                  &
                               ids,ide, jds,jde, kds,kde,       &
                               ims,ime, jms,jme, kms,kme,       &
                               ips,ipe, jps,jpe, kps,kpe,       &
                               grid%i_start(ij), grid%i_end(ij),&
                               grid%j_start(ij), grid%j_end(ij),&
                               k_start    , k_end              )

							   
						   
	    write(6,*)  'after spec_bdyupdate  grid%pfm_2  max min =  ', maxval(grid%pfm_2),minval(grid%pfm_2)						   
							   
				 write(6,*) maxval(abs(	grid%pfm_2(grid%i_start(ij):grid%i_end(ij),k_start:k_end,grid%j_start(ij):grid%j_end(ij))) ) ,  &
		            maxloc(abs(	grid%pfm_2(grid%i_start(ij):grid%i_end(ij),k_start:k_end,grid%j_start(ij):grid%j_end(ij))) )
	   

         
         ENDIF
BENCH_END(spec_bdy_t_tim)


        !!!write(6,*) 'before bc3d check _w %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    




 

! hanshenyou 





       ENDDO
       !$OMP END PARALLEL DO
	   
	   
	   
	   
	   
	 



       !$OMP PARALLEL DO   &
       !$OMP PRIVATE ( ij )
       DO ij = 1 , grid%num_tiles

BENCH_START(sumflux_tim)
         CALL sumflux ( grid%u_2, grid%v_2, grid%ww,          &
                        grid%u_save, grid%v_save, ww1,        &
                        grid%muu, grid%muv,                   &
                        grid%c1h, grid%c2h, grid%c1f, grid%c2f, &
                        grid%c3h, grid%c4h, grid%c3f, grid%c4f, &
                        grid%ru_m, grid%rv_m, grid%ww_m, grid%epssm,  &
                        grid%msfux, grid% msfuy, grid%msfvx,  &
                        grid%msfvx_inv, grid%msfvy,           &
                        iteration, number_of_small_timesteps, &
                        ids, ide, jds, jde, kds, kde,         &
                        ims, ime, jms, jme, kms, kme,         &
                        grid%i_start(ij), grid%i_end(ij),     &
                        grid%j_start(ij), grid%j_end(ij),     &
                        k_start    , k_end                   )
BENCH_END(sumflux_tim)

         IF( config_flags%specified .or. config_flags%nested ) THEN

BENCH_START(spec_bdynhyd_tim)
           IF (config_flags%non_hydrostatic)  THEN
             CALL spec_bdyupdate_ph( ph_save, grid%ph_2, ph_tend,     &
                                     mu_tend, grid%muts,              &
                                     grid%c1f, grid%c2f, dts_rk,      &
                                     'h'         , config_flags,      &
                                     grid%spec_zone,                  &
                                     ids,ide, jds,jde, kds,kde,       &
                                     ims,ime, jms,jme, kms,kme,       &
                                     ips,ipe, jps,jpe, kps,kpe,       &
                                     grid%i_start(ij), grid%i_end(ij),&
                                     grid%j_start(ij), grid%j_end(ij),&
                                     k_start    , k_end               )
             IF( config_flags%specified ) THEN
               CALL zero_grad_bdy ( grid%w_2,                         &
                                    'w'         , config_flags,       &
                                    grid%spec_zone,                   &
                                    ids,ide, jds,jde, kds,kde,        &
                                    ims,ime, jms,jme, kms,kme,        &
                                    ips,ipe, jps,jpe, kps,kpe,        &
                                    grid%i_start(ij), grid%i_end(ij), &
                                    grid%j_start(ij), grid%j_end(ij), &
                                    k_start    , k_end                )
             ELSE
               CALL spec_bdyupdate ( grid%w_2, w_tend, dts_rk,       &
                                     'h'         , config_flags,      &
                                     grid%spec_zone,                  &
                                     ids,ide, jds,jde, kds,kde,       &
                                     ims,ime, jms,jme, kms,kme,       &
                                     ips,ipe, jps,jpe, kps,kpe,       &
                                     grid%i_start(ij), grid%i_end(ij),&
                                     grid%j_start(ij), grid%j_end(ij),&
                                     k_start    , k_end               )
             ENDIF
           ENDIF
BENCH_END(spec_bdynhyd_tim)
         ENDIF



       ENDDO
       !$OMP END PARALLEL DO

#ifdef DM_PARALLEL
!
!  Stencils for patch communications  (WCS, 29 June 2001)
!
!         *                     *
!       * + *      * + *        +
!         *                     *
!
!  ph_2   x
!  al     x
!  p      x
!
!  2D variables (x,y)
!
!  mu_2   x
!  muts   x
!  mudf   x

#      include "HALO_EM_C2.inc"
#      include "PERIOD_BDY_EM_B3.inc"
#endif

BENCH_START(phys_bc_tim)
       !$OMP PARALLEL DO   &
       !$OMP PRIVATE ( ij )
       DO ij = 1 , grid%num_tiles

       ! boundary condition set for next small timestep


				 

								 
								 
! hanshenyou 20211215  20211222


	 		    !!!write(6,*) 'before set_physical_bc3d max pfm_2  :  '  , maxval(grid%pfm_2),maxloc(grid%pfm_2)  
		 	    !!!write(6,*) 'before set_physical_bc3d min pfm_2  :  '  , minval(grid%pfm_2),minloc(grid%pfm_2)
				 
         CALL set_physical_bc3d( grid%pfm_2, 'p', config_flags,            &
                                 ids, ide, jds, jde, kds, kde,     &
                                 ims, ime, jms, jme, kms, kme,     &
                                 ips, ipe, jps, jpe, kps, kpe,     &
                                 grid%i_start(ij), grid%i_end(ij), &
                                 grid%j_start(ij), grid%j_end(ij), &
                                 k_start    , k_end               )
								 
								 
								 
			    write(6,*)  'after set_physical_bc3d 1742  grid%pfm_2  max min =  ', maxval(grid%pfm_2),minval(grid%pfm_2)						   
							   
				 write(6,*) maxval(abs(	grid%pfm_2(grid%i_start(ij):grid%i_end(ij),k_start:k_end,grid%j_start(ij):grid%j_end(ij))) ) ,  &
		            maxloc(abs(	grid%pfm_2(grid%i_start(ij):grid%i_end(ij),k_start:k_end,grid%j_start(ij):grid%j_end(ij))) )
	    						 


         CALL set_physical_bc3d( grid%rhofm_2, 'p', config_flags,            &
                                 ids, ide, jds, jde, kds, kde,     &
                                 ims, ime, jms, jme, kms, kme,     &
                                 ips, ipe, jps, jpe, kps, kpe,     &
                                 grid%i_start(ij), grid%i_end(ij), &
                                 grid%j_start(ij), grid%j_end(ij), &
                                 k_start    , k_end               )

         write(6,*) '   ' 
	 	 write(6,'(a,i,5(2x,5f12.5))') 'before small_step_finish max uvwrp',  rk_step,maxval(grid%u_2),maxval(grid%v_2),maxval(grid%w_2),maxval(grid%rhofm_2),maxval(grid%pfm_2)
		 write(6,'(a,i,5(2x,5f12.5))') 'before small_step_finish min uvwrp',  rk_step,minval(grid%u_2),minval(grid%v_2),minval(grid%w_2),minval(grid%rhofm_2),minval(grid%pfm_2) 
         write(6,*) '   ' 
				
         write(6,*)  ' 	before small_step_finish max pfm_2  =  ',  maxval(abs(grid%pfm_2))  
								 

       END DO
       !$OMP END PARALLEL DO
BENCH_END(phys_bc_tim)




     END DO small_steps
	 
	 
	 


     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles

       CALL wrf_debug ( 200 , ' call rk_small_finish' )




BENCH_START(small_step_finish_tim)


	!!!!!!!write(6,*)  '===========  small_step_finish  ============  '	
	 
	 		!!!!!!!write(6,*) 'before  small_step_finish,  grid%pfm_2(i_end/2,ktf/2,1)',   grid%pfm_2(grid%i_end(ij)/2,kde/2,2)
	 
	  ! hanshenyou 20211205, still no used 
	     write(6,*)  ' 2 	before small_step_finish max pfm_2  =  ',  maxval(abs(grid%pfm_2))  
		    write(6,*)  ' 2 	before small_step_finish max pfm_save  =  ',  maxval(abs( pfm_save))  
		 
       CALL small_step_finish( grid%u_2, grid%u_1, grid%v_2, grid%v_1, grid%w_2, grid%w_1,     &
                               grid%rhofm_2, grid%rhofm_1, grid%pfm_2, grid%pfm_1, grid%sigmdot, ww1,    &  ! hanshenyou 
                               grid%mu_2, grid%mu_1,                       &
                               grid%mut, grid%muts, grid%muu, grid%muus, grid%muv, grid%muvs,  & 
                               grid%c1h, grid%c2h, grid%c1f, grid%c2f, &
                               grid%c3h, grid%c4h, grid%c3f, grid%c4f, &
                               grid%u_save, grid%v_save, w_save,           &
                               rhofm_save, pfm_save, mu_save,         &
                               grid%msfux,grid%msfuy, grid%msfvx,grid%msfvy, grid%msftx,grid%msfty, &
                               grid%h_diabatic,                       &
                               number_of_small_timesteps,dts_rk, &
                               rk_step, rk_order,                &
                               ids, ide, jds, jde, kds, kde,     &
                               ims, ime, jms, jme, kms, kme,     &
                               grid%i_start(ij), grid%i_end(ij), &
                               grid%j_start(ij), grid%j_end(ij), &
                               k_start    , k_end               )



	      write(6,*)  ' 	after small_step_finish max pfm_2  =  ',  maxval(abs(grid%pfm_2))   
		   write(6,*)  ' 	 grid%pfm_2(1,1,1,)   =  ',   grid%pfm_2(1,1,1)
		  

       IF (rk_step == rk_order) THEN

         CALL set_physical_bc3d( grid%ru_m, 'u', config_flags,   &
                                 ids, ide, jds, jde, kds, kde,      &
                                 ims, ime, jms, jme, kms, kme,      &
                                 ips, ipe, jps, jpe, kps, kpe,      &
                                 grid%i_start(ij), grid%i_end(ij),  &
                                 grid%j_start(ij), grid%j_end(ij),  &
                                 k_start    , k_end                )

         CALL set_physical_bc3d( grid%rv_m, 'v', config_flags,   &
                                 ids, ide, jds, jde, kds, kde,      &
                                 ims, ime, jms, jme, kms, kme,      &
                                 ips, ipe, jps, jpe, kps, kpe,      &
                                 grid%i_start(ij), grid%i_end(ij),  &
                                 grid%j_start(ij), grid%j_end(ij),  &
                                 k_start    , k_end                )

         CALL set_physical_bc3d( grid%ww_m, 'w', config_flags,   &
                                 ids, ide, jds, jde, kds, kde,      &
                                 ims, ime, jms, jme, kms, kme,      &
                                 ips, ipe, jps, jpe, kps, kpe,      &
                                 grid%i_start(ij), grid%i_end(ij),  &
                                 grid%j_start(ij), grid%j_end(ij),  &
                                 k_start    , k_end                )

         CALL set_physical_bc2d( grid%mut, 't', config_flags,   &
                                 ids, ide, jds, jde,               &
                                 ims, ime, jms, jme,                &
                                 ips, ipe, jps, jpe,                &
                                 grid%i_start(ij), grid%i_end(ij),  &
                                 grid%j_start(ij), grid%j_end(ij) )

         CALL set_physical_bc2d( grid%muts, 't', config_flags,   &
                                 ids, ide, jds, jde,               &
                                 ims, ime, jms, jme,                &
                                 ips, ipe, jps, jpe,                &
                                 grid%i_start(ij), grid%i_end(ij),  &
                                 grid%j_start(ij), grid%j_end(ij) )
								 
								 
! hanshenyou 20211215 
     !    CALL set_physical_bc3d( grid%ru_m, 'u', config_flags,   &

 
       END IF
	   
!=============================================================================================



!================================================================================================	   
	   
	   
	   
	   
	   
	   
	   

BENCH_END(small_step_finish_tim)

     END DO
     !$OMP END PARALLEL DO
	 
	 
	 


!-----------------------------------------------------------------------
!  add in physics tendency first if positive definite advection is used.
!  pd advection applies advective flux limiter on last runge-kutta step
!-----------------------------------------------------------------------
! first moisture

! tracer

      

#ifdef DM_PARALLEL
!
!  Stencils for patch communications  (WCS, 29 June 2001)
!
!          * * * * *            
!          * * * * *            
!          * * + * *            
!          * * * * *            
!          * * * * *            
!
! ru_m         x
! rv_m         x
! ww_m         x
! mut          x
!
!--------------------------------------------------------------

#  include "HALO_EM_D.inc"
! WCS addition 11/19/08
#  include "PERIOD_EM_DA.inc"
#endif

!<DESCRIPTION>
!<pre>
! (4) Still within the RK loop, the scalar variables are advanced.
!
!    For the moist and chem variables, each one is advanced
!    individually, using named loops "moist_variable_loop:"
!    and "chem_variable_loop:".  Each RK substep begins by
!    calculating the advective tendency, and, for the first RK step, 
!    3D mixing (calling rk_scalar_tend) followed by an update
!    of the scalar (calling rk_update_scalar).
!</pre>
!</DESCRIPTION>


       moist_scalar_advance: IF (num_3d_m >= PARAM_FIRST_SCALAR )  THEN

         moist_variable_loop: DO im = PARAM_FIRST_SCALAR, num_3d_m

! adv_moist_cond is set in module_physics_init based on mp_physics choice
!       true except for Ferrier scheme

           IF (grid%adv_moist_cond .or. im==p_qv ) THEN

             !$OMP PARALLEL DO   &
             !$OMP PRIVATE ( ij, tenddec )
             moist_tile_loop_1: DO ij = 1 , grid%num_tiles

               CALL wrf_debug ( 200 , ' call rk_scalar_tend' )
               tenddec = .false.

BENCH_START(rk_scalar_tend_tim)

        ! !!!!!!!write(6,*)     '===================   rk_scalar_tend im     =======================' 


               CALL rk_scalar_tend (  im, im, config_flags, tenddec,         & 
                           rk_step, dt_rk,                                   &
                           grid%ru_m, grid%rv_m, grid%ww_m,                  &
                           grid%muts, grid%mub, grid%mu_1,                   &
                           grid%c1h, grid%c2h,                               &
                           grid%alt,                                         &
                           moist_old(ims,kms,jms,im),                        &
                           moist(ims,kms,jms,im),                            &
                           moist_tend(ims,kms,jms,im),                       &
                           advect_tend,h_tendency,z_tendency,grid%rqvften,   & 
                           grid%qv_base, .true., grid%fnm, grid%fnp,         &
                           grid%msfux,grid%msfuy, grid%msfvx, grid%msfvx_inv,&
                           grid%msfvy, grid%msftx,grid%msfty,                & 
                           grid%rdx, grid%rdy, grid%rdn, grid%rdnw, grid%khdif, &
                           grid%kvdif, grid%xkhh,                            &
                           grid%diff_6th_opt, grid%diff_6th_factor,          &
                           config_flags%moist_adv_opt,                       &
                           grid%phb, grid%ph_2,                              &
                           config_flags%moist_mix2_off,                      &
                           config_flags%moist_mix6_off,                      &
                           ids, ide, jds, jde, kds, kde,     &
                           ims, ime, jms, jme, kms, kme,     &
                           grid%i_start(ij), grid%i_end(ij), &
                           grid%j_start(ij), grid%j_end(ij), &
                           k_start    , k_end               )

               IF( rk_step == 1 .AND. config_flags%use_q_diabatic == 1 )THEN
               IF( im.eq.p_qv .or. im.eq.p_qc )THEN
                   CALL q_diabatic_add  ( im, im,            & 
                           dt_rk, grid%mut,                  &
                           grid%c1h, grid%c2h,               &
                           grid%qv_diabatic,                 &
                           grid%qc_diabatic,                 &
                           moist_tend(ims,kms,jms,im),       &
                           ids, ide, jds, jde, kds, kde,     &
                           ims, ime, jms, jme, kms, kme,     &
                           grid%i_start(ij), grid%i_end(ij), &
                           grid%j_start(ij), grid%j_end(ij), &
                           k_start    , k_end               )
               ENDIF
               ENDIF

BENCH_END(rk_scalar_tend_tim)

BENCH_START(rlx_bdy_scalar_tim)
               IF( ( config_flags%specified .or. config_flags%nested ) .and. rk_step == 1 ) THEN 
                 IF ( im .EQ. P_QV .OR. config_flags%nested .OR. &
                    ( config_flags%specified .AND. config_flags%have_bcs_moist ) ) THEN
                   CALL relax_bdy_scalar ( moist_tend(ims,kms,jms,im),            & 
                                     moist(ims,kms,jms,im),  grid%mut,         &
                                     grid%c1h, grid%c2h,                       &
                                     moist_bxs(jms,kms,1,im),moist_bxe(jms,kms,1,im), &
                                     moist_bys(ims,kms,1,im),moist_bye(ims,kms,1,im), &
                                     moist_btxs(jms,kms,1,im),moist_btxe(jms,kms,1,im), &
                                     moist_btys(ims,kms,1,im),moist_btye(ims,kms,1,im), &
                                     config_flags%spec_bdy_width, grid%spec_zone, grid%relax_zone, &
                                     grid%dtbc, grid%fcx, grid%gcx,             &
                                     config_flags,               &
                                     ids,ide, jds,jde, kds,kde,  & ! domain dims
                                     ims,ime, jms,jme, kms,kme,  & ! memory dims
                                     ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                                     grid%i_start(ij), grid%i_end(ij),      &
                                     grid%j_start(ij), grid%j_end(ij),      &
                                     k_start, k_end                        )

                   CALL spec_bdy_scalar  ( moist_tend(ims,kms,jms,im),                &
                                     moist_bxs(jms,kms,1,im),moist_bxe(jms,kms,1,im), &
                                     moist_bys(ims,kms,1,im),moist_bye(ims,kms,1,im), &
                                     moist_btxs(jms,kms,1,im),moist_btxe(jms,kms,1,im), &
                                     moist_btys(ims,kms,1,im),moist_btye(ims,kms,1,im), &
                                     config_flags%spec_bdy_width, grid%spec_zone,                 &
                                     config_flags,               &
                                     ids,ide, jds,jde, kds,kde,  & ! domain dims
                                     ims,ime, jms,jme, kms,kme,  & ! memory dims
                                     ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                                     grid%i_start(ij), grid%i_end(ij),          &
                                     grid%j_start(ij), grid%j_end(ij),          &
                                     k_start, k_end                               )
                 ENDIF
               ENDIF
BENCH_END(rlx_bdy_scalar_tim)

             ENDDO moist_tile_loop_1
             !$OMP END PARALLEL DO

             !$OMP PARALLEL DO   &
             !$OMP PRIVATE ( ij, tenddec )
             moist_tile_loop_2: DO ij = 1 , grid%num_tiles

               CALL wrf_debug ( 200 , ' call rk_update_scalar' )
               tenddec = .false.

BENCH_START(update_scal_tim)
               CALL rk_update_scalar( scs=im, sce=im,                                  &
                               scalar_1=moist_old(ims,kms,jms,im),                     &
                               scalar_2=moist(ims,kms,jms,im),                         &
                               sc_tend=moist_tend(ims,kms,jms,im),                     &
                               advect_tend=advect_tend,                                &
                               h_tendency=h_tendency, z_tendency=z_tendency,           & 
                               msftx=grid%msftx,msfty=grid%msfty,                      &
                               c1=grid%c1h, c2=grid%c2h,                               &
                               mu_old=grid%mu_1, mu_new=grid%mu_2, mu_base=grid%mub,   &
                               rk_step=rk_step, dt=dt_rk, spec_zone=grid%spec_zone,    &
                               config_flags=config_flags, tenddec=tenddec,             & 
                               ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde,   &
                               ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme,   &
                               its=grid%i_start(ij), ite=grid%i_end(ij),               &
                               jts=grid%j_start(ij), jte=grid%j_end(ij),               &
                               kts=k_start    , kte=k_end                              )
               IF( rk_step == rk_order .AND. config_flags%use_q_diabatic == 1 )THEN
               IF( im.eq.p_qv .or. im.eq.p_qc )THEN
                   CALL q_diabatic_subtr( im, im,            & 
                           dt_rk,                            &
                           grid%qv_diabatic,                 &
                           grid%qc_diabatic,                 &
                           moist(ims,kms,jms,im),            &
                           ids, ide, jds, jde, kds, kde,     &
                           ims, ime, jms, jme, kms, kme,     &
                           grid%i_start(ij), grid%i_end(ij), &
                           grid%j_start(ij), grid%j_end(ij), &
                           k_start    , k_end               )
               ENDIF
               ENDIF
BENCH_END(update_scal_tim)

BENCH_START(flow_depbdy_tim)
               IF( config_flags%specified .AND. ( .NOT. config_flags%have_bcs_moist ) ) THEN
                 IF(im .ne. P_QV)THEN
                   CALL flow_dep_bdy  (  moist(ims,kms,jms,im),                 &
                                grid%ru_m, grid%rv_m, config_flags,             &
                                grid%spec_zone,                                 &
                                ids,ide, jds,jde, kds,kde,                      &
                                ims,ime, jms,jme, kms,kme,                      &
                                ips,ipe, jps,jpe, kps,kpe,                      &
                                grid%i_start(ij), grid%i_end(ij),               &
                                grid%j_start(ij), grid%j_end(ij),               &
                                k_start, k_end                               )
                 ENDIF
               ENDIF
BENCH_END(flow_depbdy_tim)

             ENDDO moist_tile_loop_2
             !$OMP END PARALLEL DO

           ENDIF  !-- if (grid%adv_moist_cond .or. im==p_qv ) then

         ENDDO moist_variable_loop

       ENDIF moist_scalar_advance

BENCH_START(tke_adv_tim)
       TKE_advance: IF (config_flags%km_opt .eq. 2.or.config_flags%km_opt.eq.5) then ! XZ
#ifdef DM_PARALLEL
         IF      ( config_flags%h_mom_adv_order <= 4 ) THEN
#       include "HALO_EM_TKE_ADVECT_3.inc"
         ELSE IF ( config_flags%h_mom_adv_order <= 6 ) THEN
#       include "HALO_EM_TKE_ADVECT_5.inc"
         ELSE
          WRITE(wrf_err_message,*)'solve_em: invalid h_mom_adv_order = ',config_flags%h_mom_adv_order
          CALL wrf_error_fatal(TRIM(wrf_err_message))
         ENDIF
#endif
         !$OMP PARALLEL DO   &
         !$OMP PRIVATE ( ij, tenddec )
         tke_tile_loop_1: DO ij = 1 , grid%num_tiles

           CALL wrf_debug ( 200 , ' call rk_scalar_tend for tke' )
           tenddec = .false.
           CALL rk_scalar_tend ( 1, 1, config_flags, tenddec,                      & 
                            rk_step, dt_rk,                                        &
                            grid%ru_m, grid%rv_m, grid%ww_m,                       &
                            grid%muts, grid%mub, grid%mu_1,                        &
                            grid%c1h, grid%c2h,                                    &
                            grid%alt,                                              &
                            grid%tke_1,                                            &
                            grid%tke_2,                                            &
                            tke_tend(ims,kms,jms),                                 &
                            advect_tend,h_tendency,z_tendency,grid%rqvften,        & 
                            grid%qv_base, .false., grid%fnm, grid%fnp,             &
                            grid%msfux,grid%msfuy, grid%msfvx, grid%msfvx_inv,     &
                            grid%msfvy, grid%msftx,grid%msfty,                     &
                            grid%rdx, grid%rdy, grid%rdn, grid%rdnw, grid%khdif,   &
                            grid%kvdif, grid%xkhh,                                 &
                            grid%diff_6th_opt, grid%diff_6th_factor,               &
                            config_flags%tke_adv_opt,                              &
                            grid%phb, grid%ph_2,                                   &
                            config_flags%tke_mix2_off,                             &
                            config_flags%tke_mix6_off,                             &
                            ids, ide, jds, jde, kds, kde,     &
                            ims, ime, jms, jme, kms, kme,     &
                            grid%i_start(ij), grid%i_end(ij), &
                            grid%j_start(ij), grid%j_end(ij), &
                            k_start    , k_end               )

         ENDDO tke_tile_loop_1
         !$OMP END PARALLEL DO

         !$OMP PARALLEL DO   &
         !$OMP PRIVATE ( ij, tenddec )
         tke_tile_loop_2: DO ij = 1 , grid%num_tiles

           CALL wrf_debug ( 200 , ' call rk_update_scalar' )
           tenddec = .false.
           CALL rk_update_scalar( scs=1,  sce=1,                                          &
                                  scalar_1=grid%tke_1,                                    &
                                  scalar_2=grid%tke_2,                                    &
                                  sc_tend=tke_tend(ims,kms,jms),                          &
                                  advect_tend=advect_tend,                                &
                                  h_tendency=h_tendency, z_tendency=z_tendency,           & 
                                  msftx=grid%msftx,msfty=grid%msfty,                      &
                                  c1=grid%c1h, c2=grid%c2h,                               &
                                  mu_old=grid%mu_1, mu_new=grid%mu_2, mu_base=grid%mub,   &
                                  rk_step=rk_step, dt=dt_rk, spec_zone=grid%spec_zone,    &
                                  config_flags=config_flags, tenddec=tenddec,             & 
                                  ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde,   &
                                  ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme,   &
                                  its=grid%i_start(ij), ite=grid%i_end(ij),               &
                                  jts=grid%j_start(ij), jte=grid%j_end(ij),               &
                                  kts=k_start    , kte=k_end                              )

! bound the tke (greater than 0, less than tke_upper_bound)

           CALL bound_tke( grid%tke_2, grid%tke_upper_bound,    &
                           ids, ide, jds, jde, kds, kde,        &
                           ims, ime, jms, jme, kms, kme,        &
                           grid%i_start(ij), grid%i_end(ij),    &
                           grid%j_start(ij), grid%j_end(ij),    &
                           k_start    , k_end                  )

           IF( config_flags%specified .or. config_flags%nested ) THEN
              CALL flow_dep_bdy (  grid%tke_2,                     &
                                   grid%ru_m, grid%rv_m, config_flags,               &
                                   grid%spec_zone,                              &
                                   ids,ide, jds,jde, kds,kde,  & ! domain dims
                                   ims,ime, jms,jme, kms,kme,  & ! memory dims
                                   ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                                   grid%i_start(ij), grid%i_end(ij),       &
                                   grid%j_start(ij), grid%j_end(ij),       &
                                   k_start, k_end                               )
           ENDIF
         ENDDO tke_tile_loop_2
         !$OMP END PARALLEL DO

       ENDIF TKE_advance
BENCH_END(tke_adv_tim)


!  next the chemical species

 !  update the pressure and density at the new time level

     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles

BENCH_START(calc_p_rho_tim)

  !     CALL calc_p_rho_phi( moist, num_3d_m, config_flags%hypsometric_opt,                    &
  !                          grid%al, grid%alb, grid%mu_2, grid%muts,                          &
  !                          grid%c1h, grid%c2h, grid%c3h, grid%c4h, grid%c3f, grid%c4f,       &
  !                          grid%ph_2, grid%phb, grid%p, grid%pb, grid%t_2,                   &
  !                          p0, t0, grid%p_top, grid%znu, grid%znw, grid%dnw, grid%rdnw,      &
  !                          grid%rdn, config_flags%non_hydrostatic,config_flags%use_theta_m,  &
   !                         ids, ide, jds, jde, kds, kde,                                     &
   !                         ims, ime, jms, jme, kms, kme,                                     &
   !                         grid%i_start(ij), grid%i_end(ij),                                 &
   !                         grid%j_start(ij), grid%j_end(ij),                                 &
   !                         k_start    , k_end                                                )

BENCH_END(calc_p_rho_tim)

     ENDDO
     !$OMP END PARALLEL DO

!  Reset the boundary conditions if there is another corrector step.
!  (rk_step < rk_order), else we'll handle it at the end of everything
!  (after the split physics, before exiting the timestep).

     rk_step_1_check: IF ( rk_step < rk_order ) THEN




BENCH_START(bc_end_tim)
       !$OMP PARALLEL DO   &
       !$OMP PRIVATE ( ij )
       tile_bc_loop_1: DO ij = 1 , grid%num_tiles
         CALL wrf_debug ( 200 , ' call rk_phys_bc_dry_2' )

         CALL rk_phys_bc_dry_2_fm( config_flags,                     &
                                grid%u_2, grid%v_2, grid%w_2,     &
                                grid%rhofm_2, grid%pfm_2,    &
                                ids, ide, jds, jde, kds, kde,     &
                                ims, ime, jms, jme, kms, kme,     &
                                ips, ipe, jps, jpe, kps, kpe,     &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij), &
                                k_start    , k_end               )

BENCH_START(diag_w_tim)
         IF (.not. config_flags%non_hydrostatic) THEN
           CALL diagnose_w( ph_tend, grid%ph_2, grid%ph_1, grid%w_2, grid%muts, &
                            grid%c1f, grid%c2f, dt_rk,              &
                            grid%u_2, grid%v_2, grid%ht,            &
                            grid%cf1, grid%cf2, grid%cf3, grid%rdx, grid%rdy, grid%msftx, grid%msfty, &
                            ids, ide, jds, jde, kds, kde,           &
                            ims, ime, jms, jme, kms, kme,           &
                            grid%i_start(ij), grid%i_end(ij),       &
                            grid%j_start(ij), grid%j_end(ij),       &
                            k_start    , k_end                     )
         ENDIF
BENCH_END(diag_w_tim)

         IF (num_3d_m >= PARAM_FIRST_SCALAR) THEN

           moisture_loop_bdy_1 : DO im = PARAM_FIRST_SCALAR , num_3d_m
  
             CALL set_physical_bc3d( moist(ims,kms,jms,im), 'p', config_flags,   &
                                     ids, ide, jds, jde, kds, kde,             &
                                     ims, ime, jms, jme, kms, kme,             &
                                     ips, ipe, jps, jpe, kps, kpe,             &
                                     grid%i_start(ij), grid%i_end(ij),                   &
                                     grid%j_start(ij), grid%j_end(ij),                   &
                                     k_start    , k_end                       )
           END DO moisture_loop_bdy_1

         ENDIF







         IF (config_flags%km_opt .eq. 2) THEN

           CALL set_physical_bc3d( grid%tke_2 , 'p', config_flags,  &
                                   ids, ide, jds, jde, kds, kde,            &
                                   ims, ime, jms, jme, kms, kme,            &
                                   ips, ipe, jps, jpe, kps, kpe,            &
                                   grid%i_start(ij), grid%i_end(ij),        &
                                   grid%j_start(ij), grid%j_end(ij),        &
                                   k_start    , k_end                      )
         END IF

       END DO tile_bc_loop_1
       !$OMP END PARALLEL DO
BENCH_END(bc_end_tim)





     ENDIF rk_step_1_check

   !!!!write(6,*)     '===================   Runge_Kutta_loop      =======================' 
		
		
!**********************************************************
!
!  end of RK predictor-corrector loop
!
!**********************************************************

BENCH_START(fm_test1)

       !$OMP PARALLEL DO   &
       !$OMP PRIVATE ( ij )
       DO ij = 1 , grid%num_tiles


			   

	  
       END DO 
       !$OMP END PARALLEL DO		



BENCH_END(fm_test1)


      write(6,*)  ' 	rkloop end  max pfm_2  =  ',  maxval(abs(grid%pfm_2))  , 'rk_step =  ', rk_step 
	  

		
   END DO Runge_Kutta_loop
   
!  grid%dmudt=grid%mu_2 - grid%mu_1


	 
     !!!write(6,*)  '  '
     !!!write(6,*)  '  '
	 write(6,*)  '#####################     after  Runge_Kutta_loop   max min uvwrp =     #####################' 

	 
	 	 write(6,*)  maxval(grid%u_2),maxval(grid%v_2),maxval(grid%w_2),maxval(grid%rhofm_2),maxval(grid%pfm_2)
		 write(6,*)  minval(grid%u_2),minval(grid%v_2),minval(grid%w_2),minval(grid%rhofm_2),minval(grid%pfm_2) 
		 
		 
		    
      write(6,*)  '  '
	
!         open(79,file='whsy.txt',form='binary',position='append')
!         write(79,*)  grid%pfm_2(grid%i_start(1):grid%i_end(1),k_start:k_end,2)	
!          close(79) 	
	
		
	!	stop


   IF (.not. config_flags%non_hydrostatic) THEN
#ifdef DM_PARALLEL
#    include "HALO_EM_HYDRO_UV.inc"
#    include "PERIOD_EM_HYDRO_UV.inc"
#endif
     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles
       CALL diagnose_w( ph_tend, grid%ph_2, grid%ph_1, grid%w_2, grid%muts, &
                       grid%c1f, grid%c2f, dt_rk,              &
                       grid%u_2, grid%v_2, grid%ht,            &
                       grid%cf1, grid%cf2, grid%cf3, grid%rdx, grid%rdy, grid%msftx, grid%msfty, &
                       ids, ide, jds, jde, kds, kde,           &
                       ims, ime, jms, jme, kms, kme,           &
                       grid%i_start(ij), grid%i_end(ij),       &
                       grid%j_start(ij), grid%j_end(ij),       &
                       k_start    , k_end                     )

     END DO
     !$OMP END PARALLEL DO

   END IF

   CALL wrf_debug ( 200 , ' call chem polar filter ' )



!-----------------------------------------------------------
!  end filter for chem and scalar variables at end of timestep
!-----------------------------------------------------------

   !  We're finished except for boundary condition (and patch) update

   ! Boundary condition time (or communication time).  At this time, we have
   ! implemented periodic and symmetric physical boundary conditions.

   ! b.c. routine for data within patch.

   ! we need to do both time levels of 
   ! data because the time filter only works in the physical solution space.

   ! First, do patch communications for boundary conditions (periodicity)

!-----------------------------------------------------------
!  Stencils for patch communications  (WCS, 29 June 2001)
!
!  here's where we need a wide comm stencil - these are the 
!  uncoupled variables so are used for high order calc in
!  advection and mixong routines.
!
!                              * * * * *
!            *        * * *    * * * * *
!          * + *      * + *    * * + * * 
!            *        * * *    * * * * *
!                              * * * * *
!
!   grid%u_1                            x
!   grid%u_2                            x
!   grid%v_1                            x
!   grid%v_2                            x
!   grid%w_1                            x
!   grid%w_2                            x
!   grid%t_1                            x
!   grid%t_2                            x
!  grid%ph_1                            x
!  grid%ph_2                            x
!  grid%tke_1                           x
!  grid%tke_2                           x
!
!    2D variables
!  grid%mu_1     x
!  grid%mu_2     x
!
!    4D variables
!  moist                         x
!   chem                         x
! scalar                         x
!----------------------------------------------------------


#ifdef DM_PARALLEL
   IF      ( config_flags%h_mom_adv_order <= 4 ) THEN
#    include "HALO_EM_D3_3.inc"
   ELSE IF ( config_flags%h_mom_adv_order <= 6 ) THEN
#    include "HALO_EM_D3_5.inc"
   ELSE 
      WRITE(wrf_err_message,*)'solve_em: invalid h_mom_adv_order = ',config_flags%h_mom_adv_order
      CALL wrf_error_fatal(TRIM(wrf_err_message))
   ENDIF
#  include "PERIOD_BDY_EM_D3.inc"
#  include "PERIOD_BDY_EM_MOIST.inc"
#  include "PERIOD_BDY_EM_CHEM.inc"
#  include "PERIOD_BDY_EM_TRACER.inc"
#  include "PERIOD_BDY_EM_SCALAR.inc"
#endif

!  now set physical b.c on a patch

BENCH_START(bc_2d_tim)
   !$OMP PARALLEL DO   &
   !$OMP PRIVATE ( ij )
   tile_bc_loop_2: DO ij = 1 , grid%num_tiles

     CALL wrf_debug ( 200 , ' call set_phys_bc_dry_2' )

     CALL set_phys_bc_dry_2( config_flags,                           &
                             grid%u_1, grid%u_2, grid%v_1, grid%v_2, grid%w_1, grid%w_2,           &
                             grid%t_1, grid%t_2, grid%ph_1, grid%ph_2, grid%mu_1, grid%mu_2,       &
                             ids, ide, jds, jde, kds, kde,           &
                             ims, ime, jms, jme, kms, kme,           &
                             ips, ipe, jps, jpe, kps, kpe,           &
                             grid%i_start(ij), grid%i_end(ij),       &
                             grid%j_start(ij), grid%j_end(ij),       &
                             k_start    , k_end                     )

     CALL set_physical_bc3d( grid%tke_1, 'p', config_flags,   &
                             ids, ide, jds, jde, kds, kde,            &
                             ims, ime, jms, jme, kms, kme,            &
                             ips, ipe, jps, jpe, kps, kpe,            &
                             grid%i_start(ij), grid%i_end(ij),        &
                             grid%j_start(ij), grid%j_end(ij),        &
                             k_start    , k_end-1                    )

     CALL set_physical_bc3d( grid%tke_2 , 'p', config_flags,  &
                             ids, ide, jds, jde, kds, kde,            &
                             ims, ime, jms, jme, kms, kme,            &
                             ips, ipe, jps, jpe, kps, kpe,            &
                             grid%i_start(ij), grid%i_end(ij),        &
                             grid%j_start(ij), grid%j_end(ij),        &
                             k_start    , k_end                      )

     moisture_loop_bdy_2 : DO im = PARAM_FIRST_SCALAR , num_3d_m

       CALL set_physical_bc3d( moist(ims,kms,jms,im), 'p',           &
                               config_flags,                           &
                               ids, ide, jds, jde, kds, kde,           &
                               ims, ime, jms, jme, kms, kme,           &
                               ips, ipe, jps, jpe, kps, kpe,           &
                               grid%i_start(ij), grid%i_end(ij),       &
                               grid%j_start(ij), grid%j_end(ij),       &
                               k_start    , k_end                     )

     END DO moisture_loop_bdy_2



   END DO tile_bc_loop_2
   !$OMP END PARALLEL DO
BENCH_END(bc_2d_tim)

!  this code forces boundary values to specified values to avoid drift

 
! reset surface w for consistency

#ifdef DM_PARALLEL
#  include "HALO_EM_C.inc"
#  include "PERIOD_BDY_EM_E.inc"
#endif

   CALL wrf_debug ( 10 , ' call set_w_surface' )
   fill_w_flag = .false.

   !$OMP PARALLEL DO   &
   !$OMP PRIVATE ( ij )
   DO ij = 1 , grid%num_tiles
   
   ! !!!!!!!write(6,*)  '=================  set_w_surface  ======================='
   
      CALL set_w_surface_fm( config_flags, grid%znw, fill_w_flag,              &
                           grid%w_2, grid%ht,  grid%u_2, grid%v_2,          &
                           grid%rdx, grid%rdy,                  &
                           grid%msft,                          &
                           ids, ide, jds, jde, kds, kde,                    &
                           ims, ime, jms, jme, kms, kme,                    &
                           grid%i_start(ij), grid%i_end(ij),                &
                           grid%j_start(ij), grid%j_end(ij),                &
                           k_start, k_end                                   )
!                          its, ite, jts, jte, k_start, min(k_end,kde-1),   &


							
							

   END DO
   !$OMP END PARALLEL DO

!-----------------------------------------------------------
!  After all of the RK steps, after the microphysics, after p-rho-phi,
!  after w, after filtering, we have data ready to use.
!-----------------------------------------------------------


     !!!write(6,*)  '  '
     !!!write(6,*)  '  '
	    !write(6,*)  '#####################     after  set_w_surface_fm    max min uvwrp =     #####################' 

	 	write(6,'(a,5(2x,5f15.5))') 'full max uvwrp  ',  maxval(grid%u_2),maxval(grid%v_2),maxval(grid%w_2),maxval(grid%rhofm_2),maxval(grid%pfm_2)
		write(6,'(a,5(2x,5f15.5))') 'full min uvwrp  ',  minval(grid%u_2),minval(grid%v_2),minval(grid%w_2),minval(grid%rhofm_2),minval(grid%pfm_2) 



#ifdef DM_PARALLEL
!-----------------------------------------------------------------------
! see above
!--------------------------------------------------------------
   CALL wrf_debug ( 200 , ' call HALO_RK_E' )
   IF      ( config_flags%h_mom_adv_order <= 4 ) THEN
#    include "HALO_EM_E_3.inc"
   ELSE IF ( config_flags%h_mom_adv_order <= 6 ) THEN
#    include "HALO_EM_E_5.inc"
   ELSE
     WRITE(wrf_err_message,*)'solve_em: invalid h_mom_adv_order = ',config_flags%h_mom_adv_order
     CALL wrf_error_fatal(TRIM(wrf_err_message))
   ENDIF
#endif

#ifdef DM_PARALLEL
   IF ( num_moist >= PARAM_FIRST_SCALAR  ) THEN
!-----------------------------------------------------------------------
! see above
!--------------------------------------------------------------
     CALL wrf_debug ( 200 , ' call HALO_RK_MOIST' )
     IF      ( config_flags%h_mom_adv_order <= 4 ) THEN
#      include "HALO_EM_MOIST_E_3.inc"
     ELSE IF ( config_flags%h_mom_adv_order <= 6 ) THEN
#      include "HALO_EM_MOIST_E_5.inc"
     ELSE
       WRITE(wrf_err_message,*)'solve_em: invalid h_mom_adv_order = ',config_flags%h_mom_adv_order
       CALL wrf_error_fatal(TRIM(wrf_err_message))
     ENDIF
   ENDIF

#endif

!  Max values of CFL for adaptive time step scheme

   DEALLOCATE(max_vert_cfl_tmp)
   DEALLOCATE(max_horiz_cfl_tmp)

   CALL wrf_debug ( 200 , ' call end of solve_em' )
   


 !   open(77,file='w.grd',form='binary',position='append')
  !  write(77)  grid%w_2(1:257,1:33,2)
 !   close(77)   
     !!!write(6,*) 'max-w  =  ', maxval(grid%w_2)
!     stop





     !  open(77,file='w4.grd',form='binary',position='append')
     !  write(77)  grid%u_2(1:ide,1:kde-1,1),grid%w_2(1:ide,1:kde,1),grid%rhofm_2(1:ide,1:kde-1,1),grid%pfm_2(1:ide,1:kde-1,1) ,grid%pfm_base(1:ide,1:kde-1,1) ,grid%rhofm_base(1:ide,1:kde-1,1)
     !  close(77)   
   
   
   
   
     !!!write(6,*) 'max-w  =  ', maxval(grid%w_2)
   ! stop 'stop at solve_em end '

		
	if( maxval(abs(grid%u_2)) >= 300.0) then 
		  write(6,*) 'max-u  =  ', maxval(grid%u_2)	
			
        stop
		
	 endif	

	 


  !!!write(6,*)  '=================  end of solve_em  ======================='
  !!!write(6,*) 
  !!!write(6,*) 
  !!!write(6,*) 
  !!!write(6,*) 
  !!!write(6,*) 
			
  ! stop  'stop at end of solve_em  '
   
   

!  Are we about to read SST input from the wrflowinput file?  That data is saved
!  for use in fractional merging of external/coupled SST and input SST. 
   IF ( coupler_on )   grid%just_read_auxinput4 = Is_alarm_tstep(grid%domain_clock, grid%alarms(AUXINPUT4_ALARM))

!  Are we about to read the lateral boundary file?  This is a domain one action only.
   IF ( grid%id .EQ. 1 ) grid%just_read_boundary = Is_alarm_tstep(grid%domain_clock, grid%alarms(BOUNDARY_ALARM))

! Finish timers if compiled with -DBENCH.
#include "bench_solve_em_end.h"

   RETURN

END SUBROUTINE solve_em
