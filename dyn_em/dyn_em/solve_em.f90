

SUBROUTINE solve_em ( grid , config_flags  &







,moist,moist_bxs,moist_bxe,moist_bys,moist_bye,moist_btxs,moist_btxe,moist_btys,moist_btye,dfi_moist,dfi_moist_bxs,dfi_moist_bxe, &
dfi_moist_bys,dfi_moist_bye,dfi_moist_btxs,dfi_moist_btxe,dfi_moist_btys,dfi_moist_btye,scalar,scalar_bxs,scalar_bxe,scalar_bys, &
scalar_bye,scalar_btxs,scalar_btxe,scalar_btys,scalar_btye,dfi_scalar,dfi_scalar_bxs,dfi_scalar_bxe,dfi_scalar_bys, &
dfi_scalar_bye,dfi_scalar_btxs,dfi_scalar_btxe,dfi_scalar_btys,dfi_scalar_btye,aerod,aerocu,ozmixm,aerosolc_1,aerosolc_2,fdda3d, &
fdda2d,advh_t,advz_t,tracer,tracer_bxs,tracer_bxe,tracer_bys,tracer_bye,tracer_btxs,tracer_btxe,tracer_btys,tracer_btye,nba_mij, &
nba_rij,sbmradar,chem &


                    )

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
   USE module_utility


   USE module_model_constants
   USE module_small_step_em
   USE module_em
   USE module_big_step_utilities_em
   USE module_advect_em
   USE module_bc
   USE module_bc_em
   USE module_solvedebug_em
   USE module_physics_addtendc
   USE module_diffusion_em
   USE module_polarfft
   USE module_microphysics_driver
   USE module_microphysics_zero_out

   USE module_fddaobs_driver

   USE module_dust_emis
   USE module_first_rk_step_part1
   USE module_first_rk_step_part2
   USE module_after_all_rk_steps
   USE module_llxy, ONLY : proj_cassini
   USE module_avgflx_em, ONLY : zero_avgflx, upd_avgflx
   USE module_cpl, ONLY : coupler_on, cpl_settime, cpl_store_input

   IMPLICIT NONE

   

   TYPE(domain) , TARGET          :: grid

   






real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_moist)           :: moist
real      ,DIMENSION(grid%sm33:grid%em33,grid%sm32:grid%em32,grid%spec_bdy_width,num_moist)           :: moist_bxs
real      ,DIMENSION(grid%sm33:grid%em33,grid%sm32:grid%em32,grid%spec_bdy_width,num_moist)           :: moist_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%spec_bdy_width,num_moist)           :: moist_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%spec_bdy_width,num_moist)           :: moist_bye
real      ,DIMENSION(grid%sm33:grid%em33,grid%sm32:grid%em32,grid%spec_bdy_width,num_moist)           :: moist_btxs
real      ,DIMENSION(grid%sm33:grid%em33,grid%sm32:grid%em32,grid%spec_bdy_width,num_moist)           :: moist_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%spec_bdy_width,num_moist)           :: moist_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%spec_bdy_width,num_moist)           :: moist_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_dfi_moist)           :: dfi_moist
real      ,DIMENSION(grid%sm33:grid%em33,grid%sm32:grid%em32,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bxs
real      ,DIMENSION(grid%sm33:grid%em33,grid%sm32:grid%em32,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_bye
real      ,DIMENSION(grid%sm33:grid%em33,grid%sm32:grid%em32,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btxs
real      ,DIMENSION(grid%sm33:grid%em33,grid%sm32:grid%em32,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%spec_bdy_width,num_dfi_moist)           :: dfi_moist_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_scalar)           :: scalar
real      ,DIMENSION(grid%sm33:grid%em33,grid%sm32:grid%em32,grid%spec_bdy_width,num_scalar)           :: scalar_bxs
real      ,DIMENSION(grid%sm33:grid%em33,grid%sm32:grid%em32,grid%spec_bdy_width,num_scalar)           :: scalar_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%spec_bdy_width,num_scalar)           :: scalar_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%spec_bdy_width,num_scalar)           :: scalar_bye
real      ,DIMENSION(grid%sm33:grid%em33,grid%sm32:grid%em32,grid%spec_bdy_width,num_scalar)           :: scalar_btxs
real      ,DIMENSION(grid%sm33:grid%em33,grid%sm32:grid%em32,grid%spec_bdy_width,num_scalar)           :: scalar_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%spec_bdy_width,num_scalar)           :: scalar_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%spec_bdy_width,num_scalar)           :: scalar_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_dfi_scalar)           :: dfi_scalar
real      ,DIMENSION(grid%sm33:grid%em33,grid%sm32:grid%em32,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bxs
real      ,DIMENSION(grid%sm33:grid%em33,grid%sm32:grid%em32,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_bye
real      ,DIMENSION(grid%sm33:grid%em33,grid%sm32:grid%em32,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btxs
real      ,DIMENSION(grid%sm33:grid%em33,grid%sm32:grid%em32,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%spec_bdy_width,num_dfi_scalar)           :: dfi_scalar_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_aerod)           :: aerod
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_aerocu)           :: aerocu
real      ,DIMENSION(grid%sm31:grid%em31,1:grid%levsiz,grid%sm33:grid%em33,num_ozmixm)           :: ozmixm
real      ,DIMENSION(grid%sm31:grid%em31,1:grid%paerlev,grid%sm33:grid%em33,num_aerosolc)           :: aerosolc_1
real      ,DIMENSION(grid%sm31:grid%em31,1:grid%paerlev,grid%sm33:grid%em33,num_aerosolc)           :: aerosolc_2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_fdda3d)           :: fdda3d
real      ,DIMENSION(grid%sm31:grid%em31,1:1,grid%sm33:grid%em33,num_fdda2d)           :: fdda2d
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_advh_t)           :: advh_t
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_advz_t)           :: advz_t
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_tracer)           :: tracer
real      ,DIMENSION(grid%sm33:grid%em33,grid%sm32:grid%em32,grid%spec_bdy_width,num_tracer)           :: tracer_bxs
real      ,DIMENSION(grid%sm33:grid%em33,grid%sm32:grid%em32,grid%spec_bdy_width,num_tracer)           :: tracer_bxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%spec_bdy_width,num_tracer)           :: tracer_bys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%spec_bdy_width,num_tracer)           :: tracer_bye
real      ,DIMENSION(grid%sm33:grid%em33,grid%sm32:grid%em32,grid%spec_bdy_width,num_tracer)           :: tracer_btxs
real      ,DIMENSION(grid%sm33:grid%em33,grid%sm32:grid%em32,grid%spec_bdy_width,num_tracer)           :: tracer_btxe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%spec_bdy_width,num_tracer)           :: tracer_btys
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%spec_bdy_width,num_tracer)           :: tracer_btye
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_nba_mij)           :: nba_mij
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_nba_rij)           :: nba_rij
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_sbmradar)           :: sbmradar
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_chem)           :: chem


   
   TYPE (grid_config_rec_type) , INTENT(IN)          :: config_flags

   

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
   
   
   real ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33) :: h_tendency, &
                                                                                   z_tendency
                                                                                   
   
   LOGICAL                        :: tenddec
   
   
   LOGICAL                        :: diag_flag
   INTEGER                        :: ke_diag 
   









real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: u_tendf
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: v_tendf
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: w_tend
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: w_tendf
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: pfm_tend
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: pfm_tendf
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: pfm_save
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: rhofm_tend
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: rhofm_tendf
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: rhofm_save
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: engk1
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: a_fm
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: b_fm
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: c_fm
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: ra
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: rb
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: pa
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: pb_fm
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: pc
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: pd
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: pe
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: wa_fm
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: wb_fm
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: ru_tendf
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: rv_tendf
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: ww1
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: wwp
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: rw_tend
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: rw_tendf
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: w_save
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: ph_tend
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: ph_tendf
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: ph_save
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: t_tend
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: t_tendf
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: t_2save
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: muave
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: mu_save
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: mu_tend
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: mu_tendf
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: tke_tend
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: advect_tend
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: alpha
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: a
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: gamma
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: c2a
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: phm
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: cqu
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: cqv
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: cqw
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: cqx
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: pm1
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: psim
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: psih
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: gz1oz0
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: chklowq
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: ht_loc
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: th_phy
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: pi_phy
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: p_phy
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: p8w
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: t8w
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: rho_phy
logical   ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: cu_act_flag
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: hol
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_moist)           :: moist_tend
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_moist)           :: moist_old
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_dfi_moist)           :: dfi_moist_tend
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_dfi_moist)           :: dfi_moist_old
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_scalar)           :: scalar_tend
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_scalar)           :: scalar_old
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_dfi_scalar)           :: dfi_scalar_tend
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_dfi_scalar)           :: dfi_scalar_old
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_tracer)           :: tracer_tend
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_tracer)           :: tracer_old
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_chem)           :: chem_tend
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_chem)           :: chem_old







   INTEGER :: rc 
   INTEGER :: number_of_small_timesteps, rk_step
   INTEGER :: klevel,ijm,ijp,i,j,k,size1,size2    
   INTEGER :: idum1, idum2, dynamics_option

   INTEGER :: rk_order, iwmax, jwmax, kwmax
   REAL :: dt_rk, dts_rk, dts, dtm, wmax
   REAL , ALLOCATABLE , DIMENSION(:)  :: max_vert_cfl_tmp, max_horiz_cfl_tmp
   LOGICAL :: leapfrog
   INTEGER :: l,kte,kk
   LOGICAL :: f_flux  
   REAL    :: curr_secs, curr_secs2
   INTEGER :: num_sound_steps
   INTEGER :: idex, jdex
   REAL    :: max_msft
   REAL    :: spacing

   INTEGER :: ii, jj 
   REAL    :: dclat
   INTEGER :: debug_level


   INTEGER :: NUM_ROOF_LAYERS, NUM_WALL_LAYERS, NUM_ROAD_LAYERS   

   TYPE(WRFU_TimeInterval)                    :: tmpTimeInterval, tmpTimeInterval2
   REAL                                       :: real_time
   LOGICAL                                    :: adapt_step_flag
   LOGICAL                                    :: fill_w_flag


   CHARACTER*256                              :: message, message2, message3
   REAL                                       :: old_dt
   TYPE(WRFU_Time)                            :: temp_time, CurrTime, restart_time
   INTEGER, PARAMETER                         :: precision = 100
   INTEGER                                    :: num, den
   TYPE(WRFU_TimeInterval)                    :: dtInterval, intervaltime,restartinterval




 REAL    :: tmass,teng,tengk,tie,tpe 

































































































   dynamics_option = config_flags%rk_ord



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


   if (grid%dfi_stage .EQ. DFI_BCK) then
       num_3d_m    = P_QV
       num_3d_s    = PARAM_FIRST_SCALAR - 1
   endif

   f_flux = config_flags%do_avgflx_cugd .EQ. 1



   CALL set_tiles ( ZONE_SOLVE_EM, grid , ids , ide , jds , jde , ips , ipe , jps , jpe )




   ALLOCATE (max_vert_cfl_tmp(grid%num_tiles))
   ALLOCATE (max_horiz_cfl_tmp(grid%num_tiles))

  
  
  
  
  
  
   tmpTimeInterval  = domain_get_current_time ( grid ) - domain_get_sim_start_time ( grid )
   tmpTimeInterval2 = domain_get_current_time ( grid ) - domain_get_start_time ( grid )
   curr_secs  = real_time(tmpTimeInterval)
   curr_secs2 = real_time(tmpTimeInterval2)

   old_dt = grid%dt   





   IF ( (config_flags%use_adaptive_time_step) .and. &
        ( (.not. grid%nested) .or. &
        ( (grid%nested) .and. (abs(grid%dtbc) < 0.0001) ) ) )THEN
      CALL adapt_timestep(grid, config_flags)
      adapt_step_flag = .TRUE.
   ELSE
      adapt_step_flag = .FALSE.
   ENDIF






   ke_diag = kms 
   diag_flag = .false.
   if ( Is_alarm_tstep(grid%domain_clock, grid%alarms(HISTORY_ALARM)) ) then
      diag_flag = .true.
      ke_diag = min(k_end,kde-1) 
   endif
   IF (config_flags%nwp_diagnostics == 1) diag_flag = .true.

   grid%itimestep = grid%itimestep + 1
   grid%dtbc = grid%dtbc + grid%dt

   IF( coupler_on ) CALL cpl_store_input( grid, config_flags )

   IF (config_flags%polar) dclat = 90./REAL(jde-jds) 


   rk_order = config_flags%rk_ord

   IF ( grid%time_step_sound == 0 ) THEN

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
       WRITE(wrf_err_message,*)'variable dt, max horiz cfl, max vert cfl: ',&
            grid%dt, grid%max_horiz_cfl, grid%max_vert_cfl
       CALL wrf_debug ( 0 , wrf_err_message )
     ENDIF

     grid%max_cfl_val = 0
     grid%max_horiz_cfl = 0
     grid%max_vert_cfl = 0
   ENDIF
   
   
   
		
		



     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles






       IF( config_flags%specified .AND. config_flags%constant_bc ) THEN

       CALL zero_bdytend (grid%u_btxs,grid%u_btxe,grid%u_btys,grid%u_btye,     &
                          grid%v_btxs,grid%v_btxe,grid%v_btys,grid%v_btye,     &
                          grid%pfm_btxs,grid%pfm_btxe,grid%pfm_btys,grid%pfm_btye, &
                          grid%rhofm_btxs,grid%rhofm_btxe,grid%rhofm_btys,grid%rhofm_btye,     &
                          grid%w_btxs,grid%w_btxe,grid%w_btys,grid%w_btye,     &
                   
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

       
       
       
       
       
	   


   
	
	


		
		 
		   
		 
       CALL initialize_moist_old ( moist_old(:,:,:,P_Qv),              &
                          moist(:,:,:,P_Qv) ,                          &
                          ids,ide, jds,jde, kds,kde,                   &
                          ims,ime, jms,jme, kms,kme,                   &
                          grid%i_start(ij), grid%i_end(ij),            &
                          grid%j_start(ij), grid%j_end(ij),            &
                          k_start, k_end                               )
     ENDDO
     !$OMP END PARALLEL DO

     
     
     

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
							   
							   
		
       CALL set_physical_bc3d( moist(ims,kms,jms,im), 'p', config_flags,  &
                               ids, ide, jds, jde, kds, kde,                  &
                               ims, ime, jms, jme, kms, kme,                  &
                               ips, ipe, jps, jpe, kps, kpe,                  &
                               grid%i_start(ij), grid%i_end(ij),              &
                               grid%j_start(ij), grid%j_end(ij),              &
                               k_start    , k_end                             )							   
     END DO
     !$OMP END PARALLEL DO



































   Runge_Kutta_loop:  DO rk_step = 1, rk_order

   
   
   
 
 		
     dtm = grid%dt
     IF ( rk_order == 1 ) THEN   

       write(wrf_err_message,*)' leapfrog removed, error exit for dynamics_option = ',dynamics_option
       CALL wrf_error_fatal3("<stdin>",673,&
wrf_err_message )

     ELSE IF ( rk_order == 2 ) THEN   

       IF ( rk_step == 1) THEN
         dt_rk  = 0.5*grid%dt
         dts_rk = dts
         number_of_small_timesteps = num_sound_steps/2
       ELSE
         dt_rk = grid%dt
         dts_rk = dts
         number_of_small_timesteps = num_sound_steps
       ENDIF

     ELSE IF ( rk_order == 3 ) THEN 

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
       CALL wrf_error_fatal3("<stdin>",707,&
wrf_err_message )

     END IF







     CALL wrf_debug ( 200 , ' call rk_step_prep ' )


     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )

     DO ij = 1 , grid%num_tiles



      CALL rk_step_prep_fm( config_flags, grid%alt, grid%al, grid%alb,    &
	                          cqu, cqv, cqw,cqx,                   & 
	                         moist,num_3d_m,       &
                            ids, ide, jds, jde, kds, kde,     &
                           ims, ime, jms, jme, kms, kme,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start, k_end                   )
     

 
 
	 CALL calc_p_rho_fm_new( moist, num_3d_m,                      &
	                    grid%al, grid%alt,grid%p, grid%rhofm_2, grid%rhofm_base,grid%pfm_2,grid%pfm_base,                 &
                        ids, ide, jds, jde, kds, kde,               &
                        ims, ime, jms, jme, kms, kme,               &
                        grid%i_start(ij), grid%i_end(ij),           &
                        grid%j_start(ij), grid%j_end(ij),           &
                        k_start    , k_end                          )	



 

	 	
		


          
		  
 	

          
		  
		
	 	
		
		
	 	
		 
		 
		 
		 
 END DO
     !$OMP END PARALLEL DO








     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij, ii, jj, kk )

     DO ij = 1 , grid%num_tiles

       CALL wrf_debug ( 200 , ' call rk_phys_bc_dry_1' )



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
							  


       CALL set_physical_bc3d( grid%rhofm_base, 'p', config_flags,            &      
                               ids, ide, jds, jde, kds, kde,         &
                               ims, ime, jms, jme, kms, kme,         &
                               ips, ipe, jps, jpe, kps, kpe,         &
                               grid%i_start(ij), grid%i_end(ij),     &
                               grid%j_start(ij), grid%j_end(ij),     &
                               k_start    , k_end                    )
							   
       CALL set_physical_bc3d( grid%pfm_base, 'p', config_flags,            &      
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


     rk_step_is_one : IF (rk_step == 1) THEN 



















       IF (coupler_on) CALL cpl_settime( curr_secs2 )
	   
	   	  
 

       CALL first_rk_step_part1 (    grid, config_flags         &
                             , moist , moist_tend               &
                             , chem  , chem_tend                &
                             , tracer, tracer_tend              &
                             , scalar , scalar_tend             &
                             , fdda3d, fdda2d                   &
                             , aerod                            &
                             , u_tendf, v_tendf               &
                             , w_tendf, pfm_tendf                &       
                             , rhofm_tendf,t_tendf               &
                             , tke_tend                         &
                             , config_flags%use_adaptive_time_step &
                             , curr_secs                        &
                             , psim , psih , gz1oz0             &
                             , chklowq                          &
                             , cu_act_flag , hol , th_phy       &
                             , pi_phy , p_phy , grid%t_phy      &
                             ,  p8w , t8w                 &
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
							
	  	 
		 
			  




  if(.true.) then 



 
 
	
	
		 
	
	
		 
	
	 

          CALL t_tend_to_p_tend (pfm_tendf, t_tendf, grid%rhofm_2, grid%pfm_2, th_phy,config_flags,   &
                                  moist,moist_tend,num_3d_m,                             &
                                  ids, ide, jds, jde, kds, kde,                          &
                                  ims, ime, jms, jme, kms, kme,                          &
                                  grid%i_start(ij), grid%i_end(ij),                      &
                                  grid%j_start(ij), grid%j_end(ij),                      &
                                  k_start, k_end                               )
								  
			

			
			
			
endif  




	if(.true.) then  
							
			
   
				
		       CALL first_rk_step_part2 (    grid, config_flags         &
                             , moist , moist_old , moist_tend   &
                             , chem  , chem_tend                &
                             , tracer, tracer_tend              &
                             , scalar , scalar_tend             &
                             , fdda3d, fdda2d                   &
                             , u_tendf, v_tendf               &
                             , w_tendf, rhofm_tendf, pfm_tendf                &
                             , ph_tendf, mu_tendf,t_tendf               &
                             , tke_tend                         &
                             , adapt_step_flag , curr_secs      &
                             , psim , psih , gz1oz0             &
                             , chklowq                          &
                             , cu_act_flag , hol , th_phy       &
                             , pi_phy , p_phy , grid%t_phy      &
                             , p8w , t8w                 &
                             , nba_mij, num_nba_mij             & 
                             , nba_rij, num_nba_rij             & 
                             , ids, ide, jds, jde, kds, kde     &
                             , ims, ime, jms, jme, kms, kme     &
                             , ips, ipe, jps, jpe, kps, kpe     &
                             , imsx, imex, jmsx, jmex, kmsx, kmex    &
                             , ipsx, ipex, jpsx, jpex, kpsx, kpex    &
                             , imsy, imey, jmsy, jmey, kmsy, kmey    &
                             , ipsy, ipey, jpsy, jpey, kpsy, kpey    &
                             , k_start , k_end                  &
                            )					
							
         
 
		
		  
 
	endif
	
     END IF rk_step_is_one


     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles
	 
	     

       CALL wrf_debug ( 200 , ' call rk_tendency' )
	   
	
	
        













 








			  
	   CALL rk_tendency ( config_flags, rk_step                                                                &
                         ,grid%u_tend, grid%v_tend, w_tend, pfm_tend, rhofm_tend                                 &  
                         ,u_tendf, v_tendf, w_tendf, pfm_tendf, rhofm_tendf                                      &  
						 ,grid%u_save, grid%v_save, w_save, pfm_save, rhofm_save        &    
                         ,grid%rthften                                                   &
                         ,grid%u_2, grid%v_2, grid%w_2, grid%rhofm_2, grid%pfm_2                  &   
						 ,grid%pfm_base,grid%rhofm_base           & 
						 ,grid%sigmdot                                                                            &   
                         ,grid%u_1, grid%v_1, grid%w_1, grid%rhofm_1, grid%pfm_1       &        
                         ,grid%h_diabatic, grid%phb, grid%t_init                                               &
                         ,grid%mu_2, grid%mut, grid%muu, grid%muv, grid%mub                                    &
                         ,grid%c1h, grid%c2h, grid%c1f, grid%c2f                                               &
                         ,grid%al, grid%alt, grid%p,                                          &   
						  grid%pb, grid%php, cqu, cqv, cqw,cqx                          &
                        ,grid%u_base, grid%v_base, grid%t_base, grid%qv_base, grid%z_base                     &
					 	 ,grid%msft, grid%ht,grid%zh                                          &    
						 ,grid%fzm,grid%fzp, grid%ztop                                    &    
                         ,grid%msfux,grid%msfuy, grid%msfvx, grid%msfvx_inv                                    &
                         ,grid%msfvy, grid%msftx,grid%msfty, grid%clat, grid%f, grid%e, grid%sina, grid%cosa   &
                         ,grid%fnm, grid%fnp, grid%rdn, grid%rdnw,grid%dnw                                              &
                         ,grid%dt, grid%rdx, grid%rdy, grid%khdif, grid%kvdif, grid%xkmh, grid%xkhh            &
                         ,grid%diff_6th_opt, grid%diff_6th_factor                                              &
                         ,config_flags%momentum_adv_opt                                                        &
                         ,grid%dampcoef,grid%zdamp,config_flags%damp_opt,config_flags%rad_nudge                &
                         ,grid%cf1, grid%cf2, grid%cf3, grid%cfn, grid%cfn1, num_3d_m                          &
                         ,config_flags%non_hydrostatic, config_flags%top_lid                                   &
                      
                         ,ids, ide, jds, jde, kds, kde                                                         &
                         ,ims, ime, jms, jme, kms, kme                                                         &
                         ,grid%i_start(ij), grid%i_end(ij)                                                     &
                         ,grid%j_start(ij), grid%j_end(ij)                                                     &
                         ,k_start, k_end                                                                       &
                         ,max_vert_cfl_tmp(ij), max_horiz_cfl_tmp(ij)                                         )
						 
						 
  	
  	
  	
	
	  
  	
  	
		   
     END DO
     !$OMP END PARALLEL DO





     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles

        IF ( (config_flags%specified .or. config_flags%nested) .and. ( rk_step == 1 ) ) THEN 



          if ( .true. ) then 
		  
      

 
         CALL relax_bdy_dry ( config_flags,                                &
                              grid%u_save, grid%v_save,  pfm_save, rhofm_save,                  &             
                              w_save, mu_tend, grid%c1h, grid%c2h, grid%c1f, grid%c2f, &                     
                           grid%u_2, grid%v_2, grid%pfm_2, grid%rhofm_2,                           &    
					      
                              grid%w_2, grid%mu_2, grid%mut,                              &
                              grid%u_bxs,grid%u_bxe,grid%u_bys,grid%u_bye, &
                              grid%v_bxs,grid%v_bxe,grid%v_bys,grid%v_bye, &
                              grid%pfm_bxs,grid%pfm_bxe,grid%pfm_bys,grid%pfm_bye, &                  
                              grid%rhofm_bxs,grid%rhofm_bxe,grid%rhofm_bys,grid%rhofm_bye, &          
                              grid%w_bxs,grid%w_bxe,grid%w_bys,grid%w_bye, &
                              grid%mu_bxs,grid%mu_bxe,grid%mu_bys,grid%mu_bye, &
                              grid%u_btxs,grid%u_btxe,grid%u_btys,grid%u_btye, &
                              grid%v_btxs,grid%v_btxe,grid%v_btys,grid%v_btye, &
                              grid%pfm_btxs,grid%pfm_btxe,grid%pfm_btys,grid%pfm_btye,  &                
                              grid%rhofm_btxs,grid%rhofm_btxe,grid%rhofm_btys,grid%rhofm_btye, &       
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

 
		   
		 
     
      endif 



       ENDIF
 
	
	  
	   
	   
	 

       CALL rk_addtend_dry_fm( grid%u_tend,  grid%v_tend,  w_tend,  pfm_tend,  rhofm_tend,  &
                            u_tendf, v_tendf, w_tendf, pfm_tendf, rhofm_tendf, &
                           grid%u_save, grid%v_save, w_save, pfm_save,  rhofm_save, &           
                                  rk_step,                      &
                         
                         
                         
                         
                            ids,ide, jds,jde, kds,kde,                       &
                            ims,ime, jms,jme, kms,kme,                       &
                            ips,ipe, jps,jpe, kps,kpe,                       &
                            grid%i_start(ij), grid%i_end(ij),                &
                            grid%j_start(ij), grid%j_end(ij),                &
                            k_start, k_end                                  )
							
			

 
		
							
							
		 
            
       IF( config_flags%specified .or. config_flags%nested ) THEN 
	   
	  

  
         CALL spec_bdy_dry ( config_flags,                                    &
                             grid%u_tend, grid%v_tend, pfm_tend, rhofm_tend,               &
                             w_tend, mu_tend,                                &
                             grid%u_bxs,grid%u_bxe,grid%u_bys,grid%u_bye, &
                             grid%v_bxs,grid%v_bxe,grid%v_bys,grid%v_bye, &
                             grid%pfm_bxs,grid%pfm_bxe,grid%pfm_bys,grid%pfm_bye, &
                             grid%rhofm_bxs,grid%rhofm_bxe,grid%rhofm_bys,grid%rhofm_bye, &
                             grid%w_bxs,grid%w_bxe,grid%w_bys,grid%w_bye, &
                             grid%mu_bxs,grid%mu_bxe,grid%mu_bys,grid%mu_bye, &
                             grid%u_btxs,grid%u_btxe,grid%u_btys,grid%u_btye, &
                             grid%v_btxs,grid%v_btxe,grid%v_btys,grid%v_btye, &
                             grid%pfm_btxs,grid%pfm_btxe,grid%pfm_btys,grid%pfm_btye, &   
                             grid%rhofm_btxs,grid%rhofm_btxe,grid%rhofm_btys,grid%rhofm_btye, &         
                             grid%w_btxs,grid%w_btxe,grid%w_btys,grid%w_btye, &
                             grid%mu_btxs,grid%mu_btxe,grid%mu_btys,grid%mu_btye, &
                             config_flags%spec_bdy_width, grid%spec_zone,                       &
                             ids,ide, jds,jde, kds,kde,  & 
                             ims,ime, jms,jme, kms,kme,  & 
                             ips,ipe, jps,jpe, kps,kpe,  & 
                             grid%i_start(ij), grid%i_end(ij),                &
                             grid%j_start(ij), grid%j_end(ij),                &
                             k_start, k_end                                  )
							 
							 
			   	
	 
             
		
		  
       ENDIF




     END DO
     !$OMP END PARALLEL DO



























     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles

    
    
    
    
    

       CALL wrf_debug ( 200 , ' call small_step_prep ' )
	   
	 
	 

	   
       CALL small_step_prep_fm( grid%u_1,grid%u_2,grid%v_1,grid%v_2,grid%w_1,grid%w_2,   &
                             grid%rhofm_1,grid%rhofm_2,grid%pfm_1,grid%pfm_2,                   &   
                             grid%u_save, grid%v_save, w_save,                        &
						     rhofm_save, pfm_save,                          &
                             rk_step,                             &
                             ids, ide, jds, jde, kds, kde,                            &
                             ims, ime, jms, jme, kms, kme,                            &
                             grid%i_start(ij), grid%i_end(ij),                        &
                             grid%j_start(ij), grid%j_end(ij),                        &
                             k_start    , k_end                                       )

	  	 
		         CALL calc_p_damp(  grid%pfm_2,                 &
                           pm1,  &
                           grid%smdiv,            &
                             0,    &                
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start    , k_end               )	
						  


     
       
     

     ENDDO
     !$OMP END PARALLEL DO




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
							   


				   
	 


       CALL set_physical_bc3d( grid%rhofm_2, 'p', config_flags,          &    
                               ids, ide, jds, jde, kds, kde,         &
                               ims, ime, jms, jme, kms, kme,         &
                               ips, ipe, jps, jpe, kps, kpe,         &
                               grid%i_start(ij), grid%i_end(ij),     &
                               grid%j_start(ij), grid%j_end(ij),     &
                               k_start    , k_end                    )
							   
							   
							   
       CALL set_physical_bc3d( grid%pfm_2, 'p', config_flags,          &    
                               ids, ide, jds, jde, kds, kde,         &
                               ims, ime, jms, jme, kms, kme,         &
                               ips, ipe, jps, jpe, kps, kpe,         &
                               grid%i_start(ij), grid%i_end(ij),     &
                               grid%j_start(ij), grid%j_end(ij),     &
                               k_start    , k_end                    )
							   
							   

					
       CALL set_physical_bc3d( pfm_save, 'p', config_flags,       &
                               ids, ide, jds, jde, kds, kde,         &
                               ims, ime, jms, jme, kms, kme,         &
                               ips, ipe, jps, jpe, kps, kpe,         &
                               grid%i_start(ij), grid%i_end(ij),     &
                               grid%j_start(ij), grid%j_end(ij),     &
                               k_start    , k_end                    )	

							   

	 


       CALL set_physical_bc3d( rhofm_save, 'p', config_flags,       &
                               ids, ide, jds, jde, kds, kde,         &
                               ims, ime, jms, jme, kms, kme,         &
                               ips, ipe, jps, jpe, kps, kpe,         &
                               grid%i_start(ij), grid%i_end(ij),     &
                               grid%j_start(ij), grid%j_end(ij),     &
                               k_start    , k_end                    )
							   
							   
			   
							   
							   

     END DO
     !$OMP END PARALLEL DO




       !$OMP PARALLEL DO   &
       !$OMP PRIVATE ( ij )
	   
       DO ij = 1 , grid%num_tiles



   

	         CALL calcu_abc_z( grid%u_2,grid%v_2,grid%w_2, grid%pfm_2, grid%rhofm_2,w_tend,  pfm_tend,  rhofm_tend,     &
                           grid%alt,cqw,cqx,rhofm_save,pfm_save,grid%rhofm_base, grid%pfm_base,      &
						
						   a_fm,b_fm,c_fm,wa_fm,wb_fm,ra,rb,pa,pb_fm,pc,pd,pe,      &
						   config_flags,grid%spec_zone,                               &
                           grid%msft,grid%ht,                                      &
						   grid%zh,grid%fzm,grid%fzp, grid%dnw,                             &      
                           grid%rdx, grid%rdy, dts_rk,grid%ztop,                   &				   
                           ids, ide, jds, jde, kds, kde,                          &
                           ims, ime, jms, jme, kms, kme,                          &
                           grid%i_start(ij), grid%i_end(ij),                      &
                           grid%j_start(ij), grid%j_end(ij),                      &
                           k_start    , k_end                                     )
	   









	 
	 
	 
	   


       END DO
       !$OMP END PARALLEL DO
	   
	   


     small_steps : DO iteration = 1 , number_of_small_timesteps

 
 
 

 



	  
	  
	  
	  
       

       !$OMP PARALLEL DO   &
       !$OMP PRIVATE ( ij )

       DO ij = 1 , grid%num_tiles



	   

		call calcu_engk_test( engk1,rhofm_save,grid%rhofm_2,grid%u_2+grid%u_save,grid%v_2+grid%v_save,grid%w_2+w_save,config_flags,     &
                           ids, ide, jds, jde, kds, kde,                          &
                           ims, ime, jms, jme, kms, kme,                          &
                           grid%i_start(ij), grid%i_end(ij),                      &
                           grid%j_start(ij), grid%j_end(ij),                      &
                           k_start    , k_end                                     )

         CALL set_physical_bc3d( engk1, 'p', config_flags,             &
                                 ids, ide, jds, jde, kds, kde,     &
                                 ims, ime, jms, jme, kms, kme,     &
                                 ips, ipe, jps, jpe, kps, kpe,     &
                                 grid%i_start(ij), grid%i_end(ij), &
                                 grid%j_start(ij), grid%j_end(ij), &
                                 k_start    , k_end               )	 

		 
         CALL advance_uv_fm ( grid%u_2, grid%u_tend, grid%v_2, grid%v_tend,        &
                           grid%pfm_2,grid%alt,cqu, cqv,                                               &
						   config_flags,grid%spec_zone,                               &
                           grid%msft,grid%ht,                                      &
						   grid%zh,grid%fzm,grid%fzp,                              &      
                           grid%rdx, grid%rdy, dts_rk,grid%ztop,                   &				   
                           ids, ide, jds, jde, kds, kde,                          &
                           ims, ime, jms, jme, kms, kme,                          &
                           grid%i_start(ij), grid%i_end(ij),                      &
                           grid%j_start(ij), grid%j_end(ij),                      &
                           k_start    , k_end                                     )
						   
						   
 

	  


       END DO
       !$OMP END PARALLEL DO








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



 






       !$OMP PARALLEL DO   &
       !$OMP PRIVATE ( ij )
       DO ij = 1 , grid%num_tiles


         IF( config_flags%specified .or. config_flags%nested ) THEN
           CALL spec_bdyupdate(grid%u_2, grid%u_tend, dts_rk,      &
                               'u'         , config_flags, &
                                grid%spec_zone,                  &
                                ids,ide, jds,jde, kds,kde,  & 
                                ims,ime, jms,jme, kms,kme,  & 
                                ips,ipe, jps,jpe, kps,kpe,  & 
                                grid%i_start(ij), grid%i_end(ij),         &
                                grid%j_start(ij), grid%j_end(ij),         &
                                k_start    , k_end             )

           CALL spec_bdyupdate(grid%v_2, grid%v_tend, dts_rk,      &
                                'v'         , config_flags, &
                                grid%spec_zone,                  &
                                ids,ide, jds,jde, kds,kde,  & 
                                ims,ime, jms,jme, kms,kme,  & 
                                ips,ipe, jps,jpe, kps,kpe,  & 
                                grid%i_start(ij), grid%i_end(ij),         &
                                grid%j_start(ij), grid%j_end(ij),         &
                                k_start    , k_end             )

         ENDIF


       END DO
       !$OMP END PARALLEL DO


       !$OMP PARALLEL DO   &
       !$OMP PRIVATE ( ij )
       DO ij = 1 , grid%num_tiles

        
		


  
	
         CALL advance_w_p_rho_fm_abc( grid%u_2,grid%v_2,grid%w_2, grid%pfm_2, grid%rhofm_2,w_tend,pfm_tend,rhofm_tend,      &
                           grid%alt,cqw,cqx,rhofm_save,pfm_save,grid%rhofm_base, grid%pfm_base,                &
						
						      a_fm,b_fm,c_fm,wa_fm,wb_fm,ra,rb,pa,pb_fm,pc,pd,pe,      &
                            grid%u_save,grid%v_save,w_save,grid%z,engk1,                                     & 
						   config_flags,grid%spec_zone,                               &
                           grid%msft,grid%ht,                                      &
						   grid%zh,grid%fzm,grid%fzp,grid%dnw,                              &      
                           grid%rdx, grid%rdy, dts_rk,grid%ztop, grid%z_at_w,                   &				   
                           ids, ide, jds, jde, kds, kde,                          &
                           ims, ime, jms, jme, kms, kme,                          &
                           grid%i_start(ij), grid%i_end(ij),                      &
                           grid%j_start(ij), grid%j_end(ij),                      &
                           k_start    , k_end                                     )						   




   
		
		


if(.false.) then 

         CALL calc_p_damp(  grid%pfm_2,                 &
                           pm1,  &
                           grid%smdiv,            &
                             iteration,    &
                          ids, ide, jds, jde, kds, kde,     &
                          ims, ime, jms, jme, kms, kme,     &
                          grid%i_start(ij), grid%i_end(ij), &
                          grid%j_start(ij), grid%j_end(ij), &
                          k_start    , k_end               )

endif



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
							  
 







       ENDDO
       !$OMP END PARALLEL DO





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

							   
						   

         
         ENDIF




    




 







       ENDDO
       !$OMP END PARALLEL DO
	   
	   
	   
	   
	   
	 



       !$OMP PARALLEL DO   &
       !$OMP PRIVATE ( ij )
       DO ij = 1 , grid%num_tiles


         CALL sumflux ( grid%u_2, grid%v_2, grid%w_2, grid%rhofm_2,          &
                        grid%u_save, grid%v_save, w_save,  rhofm_save,       &
                     
                     
                     
                        grid%ru_m, grid%rv_m, grid%ww_m, grid%rho,   grid%epssm,  &
                     
                     
                        iteration, number_of_small_timesteps, &
                        ids, ide, jds, jde, kds, kde,         &
                        ims, ime, jms, jme, kms, kme,         &
                        grid%i_start(ij), grid%i_end(ij),     &
                        grid%j_start(ij), grid%j_end(ij),     &
                        k_start    , k_end                   )


	 












         IF( config_flags%specified .or. config_flags%nested ) THEN


  

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



       ENDDO
       !$OMP END PARALLEL DO



       !$OMP PARALLEL DO   &
       !$OMP PRIVATE ( ij )
       DO ij = 1 , grid%num_tiles

       


				 

								 
								 



	 		    
		 	    
				 
         CALL set_physical_bc3d( grid%pfm_2, 'p', config_flags,            &
                                 ids, ide, jds, jde, kds, kde,     &
                                 ims, ime, jms, jme, kms, kme,     &
                                 ips, ipe, jps, jpe, kps, kpe,     &
                                 grid%i_start(ij), grid%i_end(ij), &
                                 grid%j_start(ij), grid%j_end(ij), &
                                 k_start    , k_end               )
								 
								 
	

         CALL set_physical_bc3d( grid%rhofm_2, 'p', config_flags,            &
                                 ids, ide, jds, jde, kds, kde,     &
                                 ims, ime, jms, jme, kms, kme,     &
                                 ips, ipe, jps, jpe, kps, kpe,     &
                                 grid%i_start(ij), grid%i_end(ij), &
                                 grid%j_start(ij), grid%j_end(ij), &
                                 k_start    , k_end               )


	 	 
		 

				

								 

       END DO
       !$OMP END PARALLEL DO





     END DO small_steps
	 
	 
	 


     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles

       CALL wrf_debug ( 200 , ' call rk_small_finish' )







  
		 
       CALL small_step_finish( grid%u_2, grid%u_1, grid%v_2, grid%v_1, grid%w_2, grid%w_1,     &
                               grid%rhofm_2, grid%rhofm_1, grid%pfm_2, grid%pfm_1, grid%sigmdot, ww1,    &  
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
		
       CALL set_physical_bc3d( grid%rho, 'p', config_flags,            &
                              ids, ide, jds, jde, kds, kde,     &
                              ims, ime, jms, jme, kms, kme,     &
                              ips, ipe, jps, jpe, kps, kpe,     &
                              grid%i_start(ij), grid%i_end(ij), &
                              grid%j_start(ij), grid%j_end(ij), &
                              k_start    , k_end               )
		
 

 
       END IF
	   






		   
	   
	   
	   
	   
	   
	   



     END DO
     !$OMP END PARALLEL DO
	 
	 
	 









     IF ((config_flags%moist_adv_opt /= ORIGINAL .and. config_flags%moist_adv_opt /= WENO_SCALAR) .and. &
         (rk_step == rk_order)) THEN

       !$OMP PARALLEL DO   &
       !$OMP PRIVATE ( ij )
       DO ij = 1 , grid%num_tiles
         CALL wrf_debug ( 200 , ' call rk_update_scalar_pd' )
         DO im = PARAM_FIRST_SCALAR, num_3d_m
		 
		   
		 
		 
		   
		   
           CALL rk_update_scalar_pd_fm( im, im,                                   &
                                     moist_old(ims,kms,jms,im),                &
                                     moist_tend(ims,kms,jms,im),               &
									 grid%rhofm_1, grid%rhofm_2,              &      
                                     rk_step, dt_rk, grid%spec_zone,           &
                                     config_flags,                             &
                                     ids, ide, jds, jde, kds, kde,             &
                                     ims, ime, jms, jme, kms, kme,             &
                                     grid%i_start(ij), grid%i_end(ij),         &
                                     grid%j_start(ij), grid%j_end(ij),         &
                                     k_start    , k_end                       )
									 
									 
 					 
									 
         ENDDO
       END DO
       !$OMP END PARALLEL DO



      



       !$OMP PARALLEL DO   &
       !$OMP PRIVATE ( ij )
       DO ij = 1 , grid%num_tiles
         IF (num_3d_m >= PARAM_FIRST_SCALAR) THEN
           DO im = PARAM_FIRST_SCALAR , num_3d_m
             CALL set_physical_bc3d( moist_old(ims,kms,jms,im), 'p', config_flags,   &
                                     ids, ide, jds, jde, kds, kde,                  &
                                     ims, ime, jms, jme, kms, kme,                  &
                                     ips, ipe, jps, jpe, kps, kpe,                  &
                                     grid%i_start(ij), grid%i_end(ij),              &
                                     grid%j_start(ij), grid%j_end(ij),              &
                                     k_start    , k_end                            )
           END DO
         ENDIF
       END DO
       !$OMP END PARALLEL DO

     END IF  
































       moist_scalar_advance: IF (num_3d_m >= PARAM_FIRST_SCALAR )  THEN

         moist_variable_loop: DO im = PARAM_FIRST_SCALAR, num_3d_m




           IF (grid%adv_moist_cond .or. im==p_qv ) THEN

             !$OMP PARALLEL DO   &
             !$OMP PRIVATE ( ij, tenddec )
             moist_tile_loop_1: DO ij = 1 , grid%num_tiles

               CALL wrf_debug ( 200 , ' call rk_scalar_tend' )
               tenddec = .false.




   
 	  

               CALL rk_scalar_tend_fm (  im, im, config_flags, tenddec,         & 
                           rk_step, dt_rk,                                   &
						    grid%rhofm_1,                                    &       
						 
                           grid%msft,grid%ht,                                      &  
						   grid%zh,grid%fzm,grid%fzp,                              &   
                           grid%ztop,                      &							    
                           grid%ru_m, grid%rv_m, grid%ww_m,                 &    
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
                   CALL q_diabatic_add_fm  ( im, im,            & 
                           dt_rk, grid%rhofm_1,     &
					
                    
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




               IF( ( config_flags%specified .or. config_flags%nested ) .and. rk_step == 1 ) THEN 
                 IF ( im .EQ. P_QV .OR. config_flags%nested .OR. &
                    ( config_flags%specified .AND. config_flags%have_bcs_moist ) ) THEN
					
 
		 
						  
                   CALL relax_bdy_scalar_fm ( moist_tend(ims,kms,jms,im),            & 
                                     moist(ims,kms,jms,im),  grid%rhofm_1,         &
                                     moist_bxs(jms,kms,1,im),moist_bxe(jms,kms,1,im), &
                                     moist_bys(ims,kms,1,im),moist_bye(ims,kms,1,im), &
                                     moist_btxs(jms,kms,1,im),moist_btxe(jms,kms,1,im), &
                                     moist_btys(ims,kms,1,im),moist_btye(ims,kms,1,im), &
                                     config_flags%spec_bdy_width, grid%spec_zone, grid%relax_zone, &
                                     grid%dtbc, grid%fcx, grid%gcx,             &
                                     config_flags,               &
                                     ids,ide, jds,jde, kds,kde,  & 
                                     ims,ime, jms,jme, kms,kme,  & 
                                     ips,ipe, jps,jpe, kps,kpe,  & 
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
                                     ids,ide, jds,jde, kds,kde,  & 
                                     ims,ime, jms,jme, kms,kme,  & 
                                     ips,ipe, jps,jpe, kps,kpe,  & 
                                     grid%i_start(ij), grid%i_end(ij),          &
                                     grid%j_start(ij), grid%j_end(ij),          &
                                     k_start, k_end                               )
									 


 
					 
 
			 
									 
                 ENDIF
               ENDIF


             ENDDO moist_tile_loop_1
             !$OMP END PARALLEL DO

             !$OMP PARALLEL DO   &
             !$OMP PRIVATE ( ij, tenddec )
             moist_tile_loop_2: DO ij = 1 , grid%num_tiles

               CALL wrf_debug ( 200 , ' call rk_update_scalar' )
               tenddec = .false.




  
	   
 
			  
               CALL rk_update_scalar( scs=im, sce=im,                                  &
                               scalar_1=moist_old(ims,kms,jms,im),                     &
                               scalar_2=moist(ims,kms,jms,im),                         &
                               sc_tend=moist_tend(ims,kms,jms,im),                     &
                               advect_tend=advect_tend,                                &
                               h_tendency=h_tendency, z_tendency=z_tendency,           & 
                             
                             
                             
							     rho_old=grid%rhofm_1, rho_new=grid%rhofm_2,                 &      
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


             ENDDO moist_tile_loop_2
             !$OMP END PARALLEL DO

           ENDIF  

         ENDDO moist_variable_loop

       ENDIF moist_scalar_advance



 
if(.false.) then 


       TKE_advance: IF (config_flags%km_opt .eq. 2.or.config_flags%km_opt.eq.5) then 
         !$OMP PARALLEL DO   &
         !$OMP PRIVATE ( ij, tenddec )
         tke_tile_loop_1: DO ij = 1 , grid%num_tiles

           CALL wrf_debug ( 200 , ' call rk_scalar_tend for tke' )
           tenddec = .false.
           CALL rk_scalar_tend_fm ( 1, 1, config_flags, tenddec,                      & 
                            rk_step, dt_rk,                                        &
			 	    grid%rhofm_1,                                    & 
			
                           grid%msft,grid%ht,                                      &  
						   grid%zh,grid%fzm,grid%fzp,                              &   
                           grid%ztop,                      &							    
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
                                 
                                 
                                 
								  rho_old=grid%rhofm_1, rho_new=grid%rhofm_2,                 &   
                                  rk_step=rk_step, dt=dt_rk, spec_zone=grid%spec_zone,    &
                                  config_flags=config_flags, tenddec=tenddec,             & 
                                  ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde,   &
                                  ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme,   &
                                  its=grid%i_start(ij), ite=grid%i_end(ij),               &
                                  jts=grid%j_start(ij), jte=grid%j_end(ij),               &
                                  kts=k_start    , kte=k_end                              )



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
                                   ids,ide, jds,jde, kds,kde,  & 
                                   ims,ime, jms,jme, kms,kme,  & 
                                   ips,ipe, jps,jpe, kps,kpe,  & 
                                   grid%i_start(ij), grid%i_end(ij),       &
                                   grid%j_start(ij), grid%j_end(ij),       &
                                   k_start, k_end                               )
           ENDIF
         ENDDO tke_tile_loop_2
         !$OMP END PARALLEL DO

       ENDIF TKE_advance



endif 


 

 

     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij )
     DO ij = 1 , grid%num_tiles



  
  
  
  
  
  
   
   
   
   
   



     ENDDO
     !$OMP END PARALLEL DO





     rk_step_1_check: IF ( rk_step < rk_order ) THEN





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






     ENDIF rk_step_1_check


		
		








       !$OMP PARALLEL DO   &
       !$OMP PRIVATE ( ij )
       DO ij = 1 , grid%num_tiles
	   
	 
	 
	  
       END DO 
       !$OMP END PARALLEL DO		







	  

		
   END DO Runge_Kutta_loop
   




 
   CALL wrf_debug ( 200 , ' call chem polar filter ' )









   
   

   

   
   

   







































 



   CALL wrf_debug ( 10 , ' call set_w_surface' )
   fill_w_flag = .false.

   !$OMP PARALLEL DO   &
   !$OMP PRIVATE ( ij )
   DO ij = 1 , grid%num_tiles
   

   
      CALL set_w_surface_fm( config_flags, grid%znw, grid%ztop, fill_w_flag,              &
                           grid%w_2, grid%ht,  grid%u_2, grid%v_2,          &
                           grid%rdx, grid%rdy,                  &
                           grid%msft,                          &
                           ids, ide, jds, jde, kds, kde,                    &
                           ims, ime, jms, jme, kms, kme,                    &
                           grid%i_start(ij), grid%i_end(ij),                &
                           grid%j_start(ij), grid%j_end(ij),                &
                           k_start, k_end                                   )



							


 								

   END DO
   !$OMP END PARALLEL DO









	    

     write(6,*)  '*********************************************  '
	 	  write(6,'(a,5(2x,5f15.5))') 'full max uvwrp  ',  maxval(grid%u_2),maxval(grid%v_2),maxval(grid%w_2),maxval(grid%rhofm_2),maxval(grid%pfm_2)
	 	  write(6,'(a,5(2x,5f15.5))') 'full min uvwrp  ',  minval(grid%u_2),minval(grid%v_2),minval(grid%w_2),minval(grid%rhofm_2),minval(grid%pfm_2) 
		 write(6,*) 'u-maxloc: ', maxloc(grid%u_2) 
		 write(6,*) 'u-minloc: ', minloc(grid%u_2) 
		 write(6,*) 'w-maxloc: ', maxloc(grid%w_2) 
		 write(6,*) 'w-minloc: ', minloc(grid%w_2) 
		 write(6,*) 'rol-maxloc: ', maxloc(grid%rhofm_2) 
		 write(6,*) 'rol-minloc: ', minloc(grid%rhofm_2) 
		 write(6,*) 'p-maxloc: ', maxloc(grid%pfm_2) 
		 write(6,*) 'p-minloc: ', minloc(grid%pfm_2) 
     write(6,*)  '*********************************************  '
  
   !$OMP PARALLEL DO   &
   !$OMP PRIVATE ( ij )
   DO ij = 1 , grid%num_tiles


CALL calcu_engmass(tmass,teng,tengk,tie,tpe,       &
                   config_flags,                                           &
                   grid%u_2, grid%v_2, grid%w_2, grid%pfm_2, grid%rhofm_2, &
                   grid%z_at_w,grid%rdx,grid%rdy,        &
                   ids, ide, jds, jde, kds, kde,                           &
                   ims, ime, jms, jme, kms, kme,                           &
                   grid%i_start(ij), grid%i_end(ij),                       &
		           grid%j_start(ij), grid%j_end(ij),                       &                              
                   k_start, k_end)
        write(6,*) ' tmass = ', tmass
    write(6,*) 'tekip =  ',  teng,tengk,tie,tpe





       write(6,*) 'changetekip =  ',teng-134534.9 ,tengk-114.9899, tie-101677.1 ,tpe-32742.83 
  END DO
   !$OMP END PARALLEL DO





   !$OMP PARALLEL DO   &
   !$OMP PRIVATE ( ij )
   DO ij = 1 , grid%num_tiles


     CALL wrf_debug ( 200 , ' call advance_ppt' )
     CALL advance_ppt(grid%rthcuten,grid%rqvcuten,grid%rqccuten,grid%rqrcuten, &
                      grid%cldfra_cup,                                         & 
                      grid%rqicuten,grid%rqscuten,           &
                      grid%rainc,grid%raincv,grid%rainsh,grid%pratec,grid%pratesh, &
                      grid%nca,grid%htop,grid%hbot,grid%cutop,grid%cubot,  &
                      grid%cuppt, grid%dt, config_flags,                   &
                      ids,ide, jds,jde, kds,kde,             &
                      ims,ime, jms,jme, kms,kme,             &
                      grid%i_start(ij), grid%i_end(ij),      &
                      grid%j_start(ij), grid%j_end(ij),      &
                      k_start    , k_end                    )


   ENDDO
  !$OMP END PARALLEL DO


if (.false. ) then  


   !$OMP PARALLEL DO   &
   !$OMP PRIVATE ( ij )
   DO ij = 1 , grid%num_tiles
        CALL wrf_debug ( 200 , ' call phy_prep_part2' )
        CALL phy_prep_part2_fm ( config_flags,                              &
		             grid%rhofm_2,           &  
                       
                       
                        grid%rthraten,                                   &
                        grid%rthblten, grid%rublten, grid%rvblten,       &
                        grid%rqvblten, grid%rqcblten, grid%rqiblten,     &
                        grid%rucuten,  grid%rvcuten,  grid%rthcuten,     &
                        grid%rqvcuten, grid%rqccuten, grid%rqrcuten,     &
                        grid%rqicuten, grid%rqscuten,                    &
                        grid%rushten,  grid%rvshten,  grid%rthshten,     &
                        grid%rqvshten, grid%rqcshten, grid%rqrshten,     &
                        grid%rqishten, grid%rqsshten, grid%rqgshten,     &
                        grid%rthften,  grid%rqvften,                     &
                        grid%RUNDGDTEN, grid%RVNDGDTEN, grid%RTHNDGDTEN, &
                        grid%RPHNDGDTEN,grid%RQVNDGDTEN, grid%RMUNDGDTEN,&
                        grid%t_2, th_phy, moist(ims,kms,jms,P_QV),       &
                        ids, ide, jds, jde, kds, kde,                    &
                        ims, ime, jms, jme, kms, kme,                    &
                        grid%i_start(ij), grid%i_end(ij),                &
                        grid%j_start(ij), grid%j_end(ij),                &
                        k_start, k_end                                   )
   ENDDO
  !$OMP END PARALLEL DO


endif 





















   IF( config_flags%specified .or. config_flags%nested ) THEN
     sz = grid%spec_zone
   ELSE
     sz = 0
   ENDIF

   IF (config_flags%mp_physics /= 0)  then

     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij, its, ite, jts, jte )

     scalar_tile_loop_1a: DO ij = 1 , grid%num_tiles

       IF ( config_flags%periodic_x ) THEN
         its = max(grid%i_start(ij),ids)
         ite = min(grid%i_end(ij),ide-1)
       ELSE
         its = max(grid%i_start(ij),ids+sz)
         ite = min(grid%i_end(ij),ide-1-sz)
       ENDIF
       jts = max(grid%j_start(ij),jds+sz)
       jte = min(grid%j_end(ij),jde-1-sz)

       CALL wrf_debug ( 200 , ' call moist_physics_prep' )




       CALL moist_physics_prep_em( grid%pfm_2, grid%rhofm_2,    grid%rho, &       
                                     grid%t_phy,p8w,t8w,           &
                                   th_phy, pi_phy , p_phy, &
                                   grid%z, grid%z_at_w, grid%dz8w,        &
                                   dtm, grid%h_diabatic,                       &
                                   moist(ims,kms,jms,P_QV),grid%qv_diabatic,   &
                                   moist(ims,kms,jms,P_QC),grid%qc_diabatic,   &
                                   config_flags,grid%fnm, grid%fnp,            &
                                   ids, ide, jds, jde, kds, kde,     &
                                   ims, ime, jms, jme, kms, kme,     &
                                   its, ite, jts, jte,               &
                                   k_start    , k_end               )
								   
								   

 
       
	   
	
       IF (config_flags%dust_emis.eq.1 .AND. config_flags%mp_physics.eq.thompsonaero)  then
         CALL wrf_debug ( 200 , ' call bulk_dust_emis' )
         CALL bulk_dust_emis (grid%itimestep,dtm,config_flags%num_soil_layers  &
               ,grid%u_phy,grid%v_phy,grid%rho,grid%alt                        &
               ,grid%u10,grid%v10,p8w,grid%dz8w,grid%smois,grid%erod                &
               ,grid%ivgtyp,grid%isltyp,grid%vegfra,grid%albbck,grid%xland     &
               ,grid%dx, g, grid%qnifa2d, ids,ide, jds,jde, kds,kde            &
               ,ims,ime, jms,jme, kms,kme                                      &
               ,its,ite, jts,jte, k_start,k_end )
       ENDIF


     END DO scalar_tile_loop_1a
     !$OMP END PARALLEL DO

     CALL wrf_debug ( 200 , ' call microphysics_driver' )

     grid%sr = 0.
     specified_bdy = config_flags%specified .OR. config_flags%nested
     channel_bdy = config_flags%specified .AND. config_flags%periodic_x






















 
   
	  

  
 
   
   
 
  
  
  if(.true.)then 
  
     CALL microphysics_driver(                                            &
      &         DT=dtm             ,DX=grid%dx              ,DY=grid%dy   &
      &        ,DZ8W=grid%dz8w          ,F_ICE_PHY=grid%f_ice_phy              &
      &        ,ITIMESTEP=grid%itimestep             ,LOWLYR=grid%lowlyr  &
      &        ,P8W=p8w            ,P=p_phy            ,PI_PHY=pi_phy     &
  
     &        ,RHO=grid%rhofm_2      &
      &  	  ,SPEC_ZONE=grid%spec_zone              &
      &        ,SR=grid%sr              ,TH=th_phy                        &
      &        ,refl_10cm=grid%refl_10cm                                  & 
      &        ,vmi3d=grid%vmi3d                                          & 
      &        ,di3d=grid%di3d                                            & 
      &        ,rhopo3d=grid%rhopo3d                                      & 
      &        ,phii3d=grid%phii3d                                        & 
      &        ,vmi3d_2=grid%vmi3d_2                                      & 
      &        ,di3d_2=grid%di3d_2                                        & 
      &        ,rhopo3d_2=grid%rhopo3d_2                                  & 
      &        ,phii3d_2=grid%phii3d_2                                    & 
      &        ,vmi3d_3=grid%vmi3d_3                                      & 
      &        ,di3d_3=grid%di3d_3                                        & 
      &        ,rhopo3d_3=grid%rhopo3d_3                                  & 
      &        ,phii3d_3=grid%phii3d_3                                    & 
      &        ,itype=grid%itype                                          & 
      &        ,itype_2=grid%itype_2                                      & 
      &        ,itype_3=grid%itype_3                                      & 
      &        ,WARM_RAIN=grid%warm_rain                                  &
      &        ,T8W=t8w                                                   &
      &        ,CLDFRA=grid%cldfra, EXCH_H=grid%exch_h &
      &        ,NSOURCE=grid%qndropsource                                 &
      &        ,XLAND=grid%xland,SNOWH=grid%SNOW                           &  
      &        ,SPECIFIED=specified_bdy, CHANNEL_SWITCH=channel_bdy       &
      &        ,F_RAIN_PHY=grid%f_rain_phy                                &
      &        ,F_RIMEF_PHY=grid%f_rimef_phy                              &
      &        ,MP_PHYSICS=config_flags%mp_physics                        &
      &        ,ID=grid%id                                                &
      &        ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde         &
      &        ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme         &
      &        ,IPS=ips,IPE=ipe, JPS=jps,JPE=jpe, KPS=kps,KPE=kpe         &
      &        ,I_START=grid%i_start,I_END=min(grid%i_end, ide-1)         &
      &        ,J_START=grid%j_start,J_END=min(grid%j_end, jde-1)         &
      &        ,KTS=k_start, KTE=min(k_end,kde-1)                         &
      &        ,NUM_TILES=grid%num_tiles                                  &
      &        ,NAER=grid%naer                                            &

      &        ,IRRIGATION=grid%irrigation                                &
      &        ,SF_SURF_IRR_SCHEME=config_flags%sf_surf_irr_scheme        &
      &        ,IRR_DAILY_AMOUNT=config_flags%irr_daily_amount            & 
      &        ,IRR_START_HOUR=config_flags%irr_start_hour                &
      &        ,IRR_NUM_HOURS=config_flags%irr_num_hours                  &
      &        ,JULIAN_IN=grid%julian                                     &
      &        ,IRR_START_JULIANDAY=config_flags%irr_start_julianday      &
      &        ,IRR_END_JULIANDAY=config_flags%irr_end_julianday          &
      &        ,IRR_FREQ=config_flags%irr_freq,IRR_PH=config_flags%irr_ph &
      &        ,IRR_RAND_FIELD=grid%irr_rand_field                        &
      &        ,GMT=grid%gmt,XTIME=grid%xtime                             & 

      
      &        ,DLF=grid%dlf,DLF2=grid%dlf2,T_PHY=grid%t_phy,P_HYD=grid%p_hyd  &
      &        ,P8W_HYD=grid%p_hyd_w,TKE_PBL=grid%tke_pbl                 &
      &        ,Z_AT_W=grid%z_at_w,QFX=grid%qfx,RLIQ=grid%rliq            &
      &        ,TURBTYPE3D=grid%turbtype3d,SMAW3D=grid%smaw3d             &
      &        ,WSEDL3D=grid%wsedl3d,CLDFRA_OLD_MP=grid%cldfra_old_mp     &
      &        ,CLDFRA_MP=grid%cldfra_mp,CLDFRA_MP_ALL=grid%cldfra_mp_ALL &
      &        ,LRADIUS=grid%LRADIUS, IRADIUS=grid%IRADIUS                & 
      &        ,CLDFRAI=grid%cldfrai                                      &
      &        ,CLDFRAL=grid%cldfral,CLDFRA_CONV=grid%CLDFRA_CONV         &
      &        ,ALT=grid%alt                                              &
      &        ,ACCUM_MODE=config_flags%accum_mode                        &
      &        ,AITKEN_MODE=config_flags%aitken_mode                      &
      &        ,COARSE_MODE=config_flags%coarse_mode                      &
      &        ,ICWMRSH3D=grid%icwmrsh,ICWMRDP3D=grid%icwmrdp3d           &
      &        ,SHFRC3D=grid%shfrc3d,CMFMC3D=grid%cmfmc                   &
      &        ,CMFMC2_3D=grid%cmfmc2,CONFIG_FLAGS=config_flags           &
      &        ,FNM=grid%fnm,FNP=grid%fnp,RH_OLD_MP=grid%rh_old_mp        &
      &        ,LCD_OLD_MP=grid%lcd_old_mp                                &

                 
      &        , RAINNC=grid%rainnc, RAINNCV=grid%rainncv                 &
      &        , SNOWNC=grid%snownc, SNOWNCV=grid%snowncv                 &
      &        , GRAUPELNC=grid%graupelnc, GRAUPELNCV=grid%graupelncv     & 
      &        , HAILNC=grid%hailnc, HAILNCV=grid%hailncv                 &
      &        , W=grid%w_2, Z=grid%z, HT=grid%ht                         &
      &        , MP_RESTART_STATE=grid%mp_restart_state                   &
      &        , TBPVS_STATE=grid%tbpvs_state                             & 
      &        , TBPVS0_STATE=grid%tbpvs0_state                           & 
      &        , QV_CURR=moist(ims,kms,jms,P_QV), F_QV=F_QV               &
      &        , QC_CURR=moist(ims,kms,jms,P_QC), F_QC=F_QC               &
      &        , QR_CURR=moist(ims,kms,jms,P_QR), F_QR=F_QR               &
      &        , QI_CURR=moist(ims,kms,jms,P_QI), F_QI=F_QI               &
      &        , QS_CURR=moist(ims,kms,jms,P_QS), F_QS=F_QS               &
      &        , QG_CURR=moist(ims,kms,jms,P_QG), F_QG=F_QG               &
      &        , QH_CURR=moist(ims,kms,jms,P_QH), F_QH=F_QH               & 
      &        , QIC_CURR=moist(ims,kms,jms,P_QIC), F_QIC=F_QIC               &
      &        , QIP_CURR=moist(ims,kms,jms,P_QIP), F_QIP=F_QIP               &
      &        , QID_CURR=moist(ims,kms,jms,P_QID), F_QID=F_QID               &
      &        , QNDROP_CURR=scalar(ims,kms,jms,P_QNDROP), F_QNDROP=F_QNDROP &
      &        , QT_CURR=scalar(ims,kms,jms,P_QT), F_QT=F_QT              &
      &        , QNN_CURR=scalar(ims,kms,jms,P_QNN), F_QNN=F_QNN          &
      &        , QNI_CURR=scalar(ims,kms,jms,P_QNI), F_QNI=F_QNI          &
      &        , QNC_CURR=scalar(ims,kms,jms,P_QNC), F_QNC=F_QNC          &
      &        , QNR_CURR=scalar(ims,kms,jms,P_QNR), F_QNR=F_QNR          &
      &        , QNS_CURR=scalar(ims,kms,jms,P_QNS), F_QNS=F_QNS          &
      &        , QNG_CURR=scalar(ims,kms,jms,P_QNG), F_QNG=F_QNG          &
      &        , QNWFA_CURR=scalar(ims,kms,jms,P_QNWFA), F_QNWFA=F_QNWFA  & 
      &        , QNIFA_CURR=scalar(ims,kms,jms,P_QNIFA), F_QNIFA=F_QNIFA  & 
      &        , QNH_CURR=scalar(ims,kms,jms,P_QNH), F_QNH=F_QNH          & 
      &        , QNIC_CURR=scalar(ims,kms,jms,P_QNIC), F_QNIC=F_QNIC          &
      &        , QNIP_CURR=scalar(ims,kms,jms,P_QNIP), F_QNIP=F_QNIP          &
      &        , QNID_CURR=scalar(ims,kms,jms,P_QNID), F_QNID=F_QNID          &
      &        , QIR_CURR=scalar(ims,kms,jms,P_QIR), F_QIR=F_QIR          & 
      &        , QIB_CURR=scalar(ims,kms,jms,P_QIB), F_QIB=F_QIB          & 
      &        , QVOLI_CURR=scalar(ims,kms,jms,P_QVOLI), F_QVOLI=F_QVOLI  & 
      &        , QAOLI_CURR=scalar(ims,kms,jms,P_QAOLI), F_QAOLI=F_QAOLI  & 
      &        , QI2_CURR=moist(ims,kms,jms,P_QI2), F_QI2=F_QI2           & 
      &        , QNI2_CURR=scalar(ims,kms,jms,P_QNI2), F_QNI2=F_QNI2      & 
      &        , QIR2_CURR=scalar(ims,kms,jms,P_QIR2), F_QIR2=F_QIR2      & 
      &        , QIB2_CURR=scalar(ims,kms,jms,P_QIB2), F_QIB2=F_QIB2      & 
      &        , QVOLI2_CURR=scalar(ims,kms,jms,P_QVOLI2), F_QVOLI2=F_QVOLI2  & 
      &        , QAOLI2_CURR=scalar(ims,kms,jms,P_QAOLI2), F_QAOLI2=F_QAOLI2  & 
      &        , QI3_CURR=moist(ims,kms,jms,P_QI3), F_QI3=F_QI3           & 
      &        , QNI3_CURR=scalar(ims,kms,jms,P_QNI3), F_QNI3=F_QNI3      & 
      &        , QVOLI3_CURR=scalar(ims,kms,jms,P_QVOLI3), F_QVOLI3=F_QVOLI3  & 
      &        , QAOLI3_CURR=scalar(ims,kms,jms,P_QAOLI3), F_QAOLI3=F_QAOLI3  & 





      &        , QVOLG_CURR=scalar(ims,kms,jms,P_QVOLG), F_QVOLG=F_QVOLG    & 
      &        , QVOLH_CURR=scalar(ims,kms,jms,P_QVOLH), F_QVOLH=F_QVOLH    & 
      &        , cu_used=config_flags%cu_used                             &
      &        , qrcuten=grid%rqrcuten, qscuten=grid%rqscuten             &
      &        , qicuten=grid%rqicuten, qccuten=grid%rqccuten             &
      &        , HAIL=config_flags%gsfcgce_hail                           & 
      &        , ICE2=config_flags%gsfcgce_2ice                           & 
      &        , PHYS_TOT=grid%phys_tot                                   & 
      &        , PHYSC=grid%physc                                         & 
      &        , PHYSE=grid%physe                                         & 
      &        , PHYSD=grid%physd                                         & 
      &        , PHYSS=grid%physs                                         & 
      &        , PHYSM=grid%physm                                         & 
      &        , PHYSF=grid%physf                                         & 

      &        , ACPHYS_TOT=grid%acphys_tot                               & 
      &        , ACPHYSC=grid%acphysc                                     & 
      &        , ACPHYSE=grid%acphyse                                     & 
      &        , ACPHYSD=grid%acphysd                                     & 
      &        , ACPHYSS=grid%acphyss                                     & 
      &        , ACPHYSM=grid%acphysm                                     & 
      &        , ACPHYSF=grid%acphysf                                     & 

      &        , RE_CLOUD_GSFC=grid%re_cloud_gsfc                         & 
      &        , RE_RAIN_GSFC=grid%re_rain_gsfc                           & 
      &        , RE_ICE_GSFC=grid%re_ice_gsfc                             & 
      &        , RE_SNOW_GSFC=grid%re_snow_gsfc                           & 
      &        , RE_GRAUPEL_GSFC=grid%re_graupel_gsfc                     & 
      &        , RE_HAIL_GSFC=grid%re_hail_gsfc                           & 
      &        , PRECR3D=grid%precr3d, PRECI3D=grid%preci3d, PRECS3D=grid%precs3d  &
      &        , PRECG3D=grid%precg3d, PRECH3D=grid%prech3d               &




      &        , RI_CURR=grid%rimi                                          &
      &        , re_cloud=grid%re_cloud, re_ice=grid%re_ice, re_snow=grid%re_snow & 
      &        , has_reqc=grid%has_reqc, has_reqi=grid%has_reqi, has_reqs=grid%has_reqs & 
      &        , qnwfa2d=grid%qnwfa2d, qnifa2d=grid%qnifa2d                             & 
      &        , diagflag=diag_flag, do_radar_ref=config_flags%do_radar_ref &
      &        , ke_diag=ke_diag                                           &
      &        ,u=grid%u_phy,v=grid%v_phy &
      &        ,scalar=scalar,num_scalar=num_scalar                             &
      &        ,TH_OLD=grid%th_old                                        &
      &        ,QV_OLD=grid%qv_old                                        &
      &        ,xlat=grid%xlat,xlong=grid%xlong,IVGTYP=grid%ivgtyp  &
      &        , EFFR_CURR=scalar(ims,kms,jms,P_EFFR), F_EFFR=F_EFFR          & 
      &        , ICE_EFFR_CURR=scalar(ims,kms,jms,P_ICE_EFFR), F_ICE_EFFR=F_ICE_EFFR          & 
      &        , TOT_EFFR_CURR=scalar(ims,kms,jms,P_TOT_EFFR), F_TOT_EFFR=F_TOT_EFFR          & 
      &        , QIC_EFFR_CURR=scalar(ims,kms,jms,P_QIC_EFFR), F_QIC_EFFR=F_QIC_EFFR          & 
      &        , QIP_EFFR_CURR=scalar(ims,kms,jms,P_QIP_EFFR), F_QIP_EFFR=F_QIP_EFFR          & 
      &        , QID_EFFR_CURR=scalar(ims,kms,jms,P_QID_EFFR), F_QID_EFFR=F_QID_EFFR          & 
      &        ,kext_ql=grid%kext_ql                                       &
      &        ,kext_qs=grid%kext_qs                                       &
      &        ,kext_qg=grid%kext_qg                                       &
      &        ,kext_qh=grid%kext_qh                                       &
      &        ,kext_qa=grid%kext_qa                                       &
      &        ,kext_qic=grid%kext_qic                                       &
      &        ,kext_qip=grid%kext_qip                                       &
      &        ,kext_qid=grid%kext_qid                                       &
      &        ,kext_ft_qic=grid%kext_ft_qic                                       &
      &        ,kext_ft_qip=grid%kext_ft_qip                                       &
      &        ,kext_ft_qid=grid%kext_ft_qid                                       &
      &        ,kext_ft_qs=grid%kext_ft_qs                                       &
      &        ,kext_ft_qg=grid%kext_ft_qg         &
      &        ,height=grid%height                                         &
      &        ,tempc=grid%tempc                                         &
      &        ,ccn_conc=grid%ccn_conc                                   & 
      &        ,sbmradar=sbmradar,num_sbmradar=num_sbmradar              & 
      &        ,sbm_diagnostics=config_flags%sbm_diagnostics             & 
      &        ,aerocu=aerocu                                            &
      &        ,aercu_fct=config_flags%aercu_fct                         &
      &        ,aercu_opt=config_flags%aercu_opt                         &
      &        ,no_src_types_cu=grid%no_src_types_cu                     &
      &        ,PBL=grid%bl_pbl_physics,EFCG=grid%EFCG,EFIG=grid%EFIG,EFSG=grid%EFSG &
      &        ,WACT=grid%WACT,CCN1_GS=grid%CCN1_GS,CCN2_GS=grid%CCN2_GS,CCN3_GS=grid%CCN3_GS  &
      &        ,CCN4_GS=grid%CCN4_GS,CCN5_GS=grid%CCN5_GS,CCN6_GS=grid%CCN6_GS &
      &        ,CCN7_GS=grid%CCN7_GS,NR_CU=grid%NR_CU,QR_CU=grid%QR_CU,NS_CU=grid%NS_CU &
      &        ,QS_CU=grid%QS_CU,CU_UAF=grid%CU_UAF,mskf_refl_10cm=grid%mskf_refl_10cm)
 



endif 


 





	


 
	
	
     CALL wrf_debug ( 200 , ' call moist_physics_finish' )


     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij, its, ite, jts, jte, im, ii, jj, kk )

     DO ij = 1 , grid%num_tiles

       its = max(grid%i_start(ij),ids)
       ite = min(grid%i_end(ij),ide-1)
       jts = max(grid%j_start(ij),jds)
       jte = min(grid%j_end(ij),jde-1)

       CALL microphysics_zero_outb (                                    &
                      moist , num_moist , config_flags ,                &
                      ids, ide, jds, jde, kds, kde,                     &
                      ims, ime, jms, jme, kms, kme,                     &
                      its, ite, jts, jte,                               &
                      k_start    , k_end                                )

       CALL microphysics_zero_outb (                                    &
                      scalar , num_scalar , config_flags ,              &
                      ids, ide, jds, jde, kds, kde,                     &
                      ims, ime, jms, jme, kms, kme,                     &
                      its, ite, jts, jte,                               &
                      k_start    , k_end                                )

       CALL microphysics_zero_outb (                                    &
                      chem , num_chem , config_flags ,              &
                      ids, ide, jds, jde, kds, kde,                     &
                      ims, ime, jms, jme, kms, kme,                     &
                      its, ite, jts, jte,                               &
                      k_start    , k_end                                )
       CALL microphysics_zero_outb (                                    &
                      tracer , num_tracer , config_flags ,              &
                      ids, ide, jds, jde, kds, kde,                     &
                      ims, ime, jms, jme, kms, kme,                     &
                      its, ite, jts, jte,                               &
                      k_start    , k_end                                )

       IF ( config_flags%periodic_x ) THEN
         its = max(grid%i_start(ij),ids)
         ite = min(grid%i_end(ij),ide-1)
       ELSE
         its = max(grid%i_start(ij),ids+sz)
         ite = min(grid%i_end(ij),ide-1-sz)
       ENDIF
       jts = max(grid%j_start(ij),jds+sz)
       jte = min(grid%j_end(ij),jde-1-sz)

       CALL microphysics_zero_outa (                                    &
                      moist , num_moist , config_flags ,                &
                      ids, ide, jds, jde, kds, kde,                     &
                      ims, ime, jms, jme, kms, kme,                     &
                      its, ite, jts, jte,                               &
                      k_start    , k_end                                )

       CALL microphysics_zero_outa (                                    &
                      scalar , num_scalar , config_flags ,              &
                      ids, ide, jds, jde, kds, kde,                     &
                      ims, ime, jms, jme, kms, kme,                     &
                      its, ite, jts, jte,                               &
                      k_start    , k_end                                )

       CALL microphysics_zero_outa (                                    &
                      chem , num_chem , config_flags ,                  &
                      ids, ide, jds, jde, kds, kde,                     &
                      ims, ime, jms, jme, kms, kme,                     &
                      its, ite, jts, jte,                               &
                      k_start    , k_end                                )

       CALL microphysics_zero_outa (                                    &
                      tracer , num_tracer , config_flags ,              &
                      ids, ide, jds, jde, kds, kde,                     &
                      ims, ime, jms, jme, kms, kme,                     &
                      its, ite, jts, jte,                               &
                      k_start    , k_end                                )

       
 
     
	 
 
 
			
       CALL moist_physics_finish_fm( grid%t_2, grid%t_1, t0,   &
	                                 grid%pfm_2,grid%rhofm_2,        & 
	                                  
									     th_phy,  &
                                     grid%h_diabatic, dtm,                       &
                                     moist(ims,kms,jms,P_QV),grid%qv_diabatic,   &
                                     moist(ims,kms,jms,P_QC),grid%qc_diabatic,   &
                                     grid%th_phy_m_t0,                           &
                                     config_flags,                               &
                                     ids, ide, jds, jde, kds, kde,               &
                                     ims, ime, jms, jme, kms, kme,               &
                                     its, ite, jts, jte,                         &
                                     k_start    , k_end                          )

 
    
 



	 
	 
	 

	   
 write(6,*)  'max grid%rainnc =  ', maxval(grid%rainnc) ,   maxloc(grid%rainnc)
 write(6,*)  'max grid%rainc =  ', maxval(grid%rainc) ,   maxloc(grid%rainc)

 


 	 


		
 END DO
     !$OMP END PARALLEL DO

        its=ips ; ite = ipe
        jts=jps ; jte = jpe
        CALL set_physical_bc3d( grid%h_diabatic, 'p', config_flags,      &
                                ids, ide, jds, jde, kds, kde,     &
                                ims, ime, jms, jme, kms, kme,     &
                                ips, ipe, jps, jpe, kps, kpe,     &
                                its, ite, jts, jte,               &
                                k_start    , k_end                )
   ENDIF  












 
	




		












   !$OMP PARALLEL DO   &
   !$OMP PRIVATE ( ij )
   tile_bc_loop_2: DO ij = 1 , grid%num_tiles

     CALL wrf_debug ( 200 , ' call set_phys_bc_dry_2' )

  
   
   
   
    
    
    
    
    
							 
							 
         CALL rk_phys_bc_dry_2_fm( config_flags,                     &
                                grid%u_2, grid%v_2, grid%w_2,     &
                                grid%rhofm_2, grid%pfm_2,    &
                                ids, ide, jds, jde, kds, kde,     &
                                ims, ime, jms, jme, kms, kme,     &
                                ips, ipe, jps, jpe, kps, kpe,     &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij), &
                                k_start    , k_end               )							 


         CALL rk_phys_bc_dry_2_fm( config_flags,                     &
                                grid%u_1, grid%v_1, grid%w_1,     &
                                grid%rhofm_1, grid%pfm_1,    &
                                ids, ide, jds, jde, kds, kde,     &
                                ims, ime, jms, jme, kms, kme,     &
                                ips, ipe, jps, jpe, kps, kpe,     &
                                grid%i_start(ij), grid%i_end(ij), &
                                grid%j_start(ij), grid%j_end(ij), &
                                k_start    , k_end               )

								
							 

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





   
   
   IF( config_flags%specified .or. config_flags%nested ) THEN 
   !$OMP PARALLEL DO   &
   !$OMP PRIVATE ( ij )
   tile_bc_loop_3: DO ij = 1 , grid%num_tiles

     CALL wrf_debug ( 200 , ' call spec_bdy_final' )
	 
	 
	 		
		    

     CALL spec_bdy_final_fm   ( grid%u_2,     &
	                 
                                grid%u_bxs, grid%u_bxe, grid%u_bys, grid%u_bye,  &
                                grid%u_btxs,grid%u_btxe,grid%u_btys,grid%u_btye, &
                                'u', config_flags,                               &
                                config_flags%spec_bdy_width, grid%spec_zone,     &
                                grid%dtbc,                                       &
                                ids,ide, jds,jde, kds,kde,  & 
                                ims,ime, jms,jme, kms,kme,  & 
                                ips,ipe, jps,jpe, kps,kpe,  & 
                                grid%i_start(ij), grid%i_end(ij),       &
                                grid%j_start(ij), grid%j_end(ij),       &
                                k_start    , k_end                     )
								
								
 						

     CALL spec_bdy_final_fm   ( grid%v_2,         &
	                         
                                grid%v_bxs, grid%v_bxe, grid%v_bys, grid%v_bye,  &
                                grid%v_btxs,grid%v_btxe,grid%v_btys,grid%v_btye, &
                                'v', config_flags,                               &
                                config_flags%spec_bdy_width, grid%spec_zone,     &
                                grid%dtbc,                                       &
                                ids,ide, jds,jde, kds,kde,  & 
                                ims,ime, jms,jme, kms,kme,  & 
                                ips,ipe, jps,jpe, kps,kpe,  & 
                                grid%i_start(ij), grid%i_end(ij),       &
                                grid%j_start(ij), grid%j_end(ij),       &
                                k_start    , k_end                     )

     IF( config_flags%nested) THEN
       CALL spec_bdy_final   ( grid%w_2, grid%muts, grid%c1f, grid%c2f, grid%msfty, &
                                grid%w_bxs, grid%w_bxe, grid%w_bys, grid%w_bye,  &
                                grid%w_btxs,grid%w_btxe,grid%w_btys,grid%w_btye, &
                                'w', config_flags,                               &
                                config_flags%spec_bdy_width, grid%spec_zone,     &
                                grid%dtbc,                                       &
                                ids,ide, jds,jde, kds,kde,  & 
                                ims,ime, jms,jme, kms,kme,  & 
                                ips,ipe, jps,jpe, kps,kpe,  & 
                                grid%i_start(ij), grid%i_end(ij),       &
                                grid%j_start(ij), grid%j_end(ij),       &
                                k_start    , k_end                     )
     ENDIF

     CALL spec_bdy_final_fm   ( grid%rhofm_2,    &
	                    
                                grid%rhofm_bxs, grid%rhofm_bxe, grid%rhofm_bys, grid%rhofm_bye,  &
                                grid%rhofm_btxs,grid%rhofm_btxe,grid%rhofm_btys,grid%rhofm_btye, &
                                't', config_flags,                               &
                                config_flags%spec_bdy_width, grid%spec_zone,     &
                                grid%dtbc,                                       &
                                ids,ide, jds,jde, kds,kde,  & 
                                ims,ime, jms,jme, kms,kme,  & 
                                ips,ipe, jps,jpe, kps,kpe,  & 
                                grid%i_start(ij), grid%i_end(ij),       &
                                grid%j_start(ij), grid%j_end(ij),       &
                                k_start    , k_end                     )

     CALL spec_bdy_final_fm   ( grid%pfm_2,               &
	                        
                                grid%pfm_bxs, grid%pfm_bxe, grid%pfm_bys, grid%pfm_bye,  &
                                grid%pfm_btxs,grid%pfm_btxe,grid%pfm_btys,grid%pfm_btye, &

								  't', config_flags,                               &
                                config_flags%spec_bdy_width, grid%spec_zone,     &
                                grid%dtbc,                                       &
                                ids,ide, jds,jde, kds,kde,  & 
                                ims,ime, jms,jme, kms,kme,  & 
                                ips,ipe, jps,jpe, kps,kpe,  & 
                                grid%i_start(ij), grid%i_end(ij),       &
                                grid%j_start(ij), grid%j_end(ij),       &
                                k_start    , k_end                     )
								
								
								
							

   
   
   
   

   
   
   
   
   
   
   
   
   

     moisture_loop_bdy_3 : DO im = PARAM_FIRST_SCALAR , num_3d_m

     IF ( im .EQ. P_QV .OR. config_flags%nested .OR. &
             ( config_flags%specified .AND. config_flags%have_bcs_moist ) ) THEN
        CALL spec_bdy_final_fm_moist   ( moist(ims,kms,jms,im),  grid%rhofm_2,    &
		                      
                               
                                moist_bxs(jms,kms,1,im),moist_bxe(jms,kms,1,im), &
                                moist_bys(ims,kms,1,im),moist_bye(ims,kms,1,im), &
                                moist_btxs(jms,kms,1,im),moist_btxe(jms,kms,1,im), &
                                moist_btys(ims,kms,1,im),moist_btye(ims,kms,1,im), &
                                't', config_flags,                               &
                                config_flags%spec_bdy_width, grid%spec_zone,     &
                                grid%dtbc,                                       &
                                ids,ide, jds,jde, kds,kde,  & 
                                ims,ime, jms,jme, kms,kme,  & 
                                ips,ipe, jps,jpe, kps,kpe,  & 
                                grid%i_start(ij), grid%i_end(ij),       &
                                grid%j_start(ij), grid%j_end(ij),       &
                                k_start    , k_end                     )
								
		    
			 
			  
			
		
			
     ENDIF

     END DO moisture_loop_bdy_3


     tracer_species_bdy_loop_3 : DO im = PARAM_FIRST_SCALAR , num_tracer

     IF( ( config_flags%nested ) ) THEN
        CALL spec_bdy_final_fm   ( tracer(ims,kms,jms,im),    & 
                    
                    
                                tracer_bxs(jms,kms,1,im),tracer_bxe(jms,kms,1,im), &
                                tracer_bys(ims,kms,1,im),tracer_bye(ims,kms,1,im), &
                                tracer_btxs(jms,kms,1,im),tracer_btxe(jms,kms,1,im), &
                                tracer_btys(ims,kms,1,im),tracer_btye(ims,kms,1,im), &
                                't', config_flags,                               &
                                config_flags%spec_bdy_width, grid%spec_zone,     &
                                grid%dtbc,                                       &
                                ids,ide, jds,jde, kds,kde,  & 
                                ims,ime, jms,jme, kms,kme,  & 
                                ips,ipe, jps,jpe, kps,kpe,  & 
                                grid%i_start(ij), grid%i_end(ij),       &
                                grid%j_start(ij), grid%j_end(ij),       &
                                k_start    , k_end                     )
     ENDIF

     END DO tracer_species_bdy_loop_3

     scalar_species_bdy_loop_3 : DO is = PARAM_FIRST_SCALAR , num_3d_s

     IF( ( config_flags%nested ) ) THEN
        CALL spec_bdy_final_fm   ( scalar(ims,kms,jms,is),      &
		                       
                               
                                scalar_bxs(jms,kms,1,is),scalar_bxe(jms,kms,1,is), &
                                scalar_bys(ims,kms,1,is),scalar_bye(ims,kms,1,is), &
                                scalar_btxs(jms,kms,1,is),scalar_btxe(jms,kms,1,is), &
                                scalar_btys(ims,kms,1,is),scalar_btye(ims,kms,1,is), &
                                't', config_flags,                               &
                                config_flags%spec_bdy_width, grid%spec_zone,     &
                                grid%dtbc,                                       &
                                ids,ide, jds,jde, kds,kde,  & 
                                ims,ime, jms,jme, kms,kme,  & 
                                ips,ipe, jps,jpe, kps,kpe,  & 
                                grid%i_start(ij), grid%i_end(ij),       &
                                grid%j_start(ij), grid%j_end(ij),       &
                                k_start    , k_end                     )
     ENDIF

     END DO scalar_species_bdy_loop_3

   END DO tile_bc_loop_3
   !$OMP END PARALLEL DO

   ENDIF




   CALL wrf_debug ( 10 , ' call set_w_surface' )
   fill_w_flag = .false.

   !$OMP PARALLEL DO   &
   !$OMP PRIVATE ( ij )
   DO ij = 1 , grid%num_tiles

      CALL set_w_surface_fm( config_flags, grid%znw, grid%ztop, fill_w_flag,              &
                           grid%w_2, grid%ht,  grid%u_2, grid%v_2,          &
                           grid%rdx, grid%rdy,                  &
                           grid%msft,                          &
                           ids, ide, jds, jde, kds, kde,                    &
                           ims, ime, jms, jme, kms, kme,                    &
                           grid%i_start(ij), grid%i_end(ij),                &
                           grid%j_start(ij), grid%j_end(ij),                &
                           k_start, k_end                                   )








 END DO
   !$OMP END PARALLEL DO






  CALL after_all_rk_steps ( grid, config_flags,                  &
                            moist, chem, tracer, scalar,         &
                            th_phy, pi_phy, p_phy, rho_phy,      &   
                            p8w, t8w, grid%dz8w,                      &
                            REAL(curr_secs,8), curr_secs2,       &
                            diag_flag,                           &
                            ids,  ide,  jds,  jde,  kds,  kde,   &
                            ims,  ime,  jms,  jme,  kms,  kme,   &
                            ips,  ipe,  jps,  jpe,  kps,  kpe,   &
                            imsx, imex, jmsx, jmex, kmsx, kmex,  &
                            ipsx, ipex, jpsx, jpex, kpsx, kpex,  &
                            imsy, imey, jmsy, jmey, kmsy, kmey,  &
                            ipsy, ipey, jpsy, jpey, kpsy, kpey   )








































   DEALLOCATE(max_vert_cfl_tmp)
   DEALLOCATE(max_horiz_cfl_tmp)

   CALL wrf_debug ( 200 , ' call end of solve_em' )
   



  
 







   


		
	if( maxval(abs(grid%w_2)) >= 100.0) then 
		 write(6,*) 'w is too big so stop for max-w  =  ', maxval(grid%w_2)	
		 write(6,*) 'w is too big so stop for min-w  =  ', minval(grid%w_2)	
		 write(6,*) 'maxloc: ', maxloc(grid%w_2) 
		 write(6,*) 'minloc: ', minloc(grid%w_2) 		
        stop
		
	 endif	

	if( maxval(abs(grid%u_2)) >= 110.0) then 
		 write(6,*) 'u is too big so stop for max-u  =  ', maxval(grid%u_2)	
		 write(6,*) 'u is too big so stop for min-u  =  ', minval(grid%u_2) 
		 write(6,*) 'maxloc: ', maxloc(grid%u_2) 
		 write(6,*) 'minloc: ', minloc(grid%u_2) 
        stop
		
	 endif	
	 
	if( maxval(abs(grid%v_2)) >= 110.0) then 
		 write(6,*) 'v is too big so stop for max-v  =  ', maxval(grid%v_2)	
		 write(6,*) 'v is too big so stop for min-v  =  ', minval(grid%v_2)		
		 write(6,*) 'maxloc: ', maxloc(grid%v_2) 
		 write(6,*) 'minloc: ', minloc(grid%v_2) 
        stop
		
	 endif		 












			

   
   



   IF ( coupler_on )   grid%just_read_auxinput4 = Is_alarm_tstep(grid%domain_clock, grid%alarms(AUXINPUT4_ALARM))


   IF ( grid%id .EQ. 1 ) grid%just_read_boundary = Is_alarm_tstep(grid%domain_clock, grid%alarms(BOUNDARY_ALARM))



























































   RETURN
   
   

END SUBROUTINE solve_em
