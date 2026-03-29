







MODULE module_small_step_em

   USE module_configure
   USE module_model_constants
   USE module_advect_em, only: calcu_gradx,calcu_grady,calcu_gradz,calcu_divg,calcu_divg_pres,calcu_divg_p,  &
                                  calculate_sgmdot,calculate_sgmdot2,calcu_engk,     &
								   average_vertical_w,average_vertical_pw        
   

CONTAINS




SUBROUTINE small_step_prep_fm( u_1, u_2, v_1, v_2, w_1, w_2, &          
                            rho_1, rho_2, p_1, p_2,         &  
                            u_save, v_save, w_save,       &
                            rho_save, p_save,          &    
                            rk_step,                      &
                            ids,ide, jds,jde, kds,kde,    &
                            ims,ime, jms,jme, kms,kme,    &
                            its,ite, jts,jte, kts,kte    )



  IMPLICIT NONE  



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





  INTEGER :: i, j, k
  INTEGER :: i_start, i_end, j_start, j_end, k_start, k_end
  INTEGER :: i_endu, j_endv












    i_start = its
    i_end   = min(ite,ide-1)
    j_start = jts
    j_end   = min(jte,jde-1)
    k_start = kts
    k_end = min(kte,kde-1)

    i_endu = ite
    j_endv = jte

    
    

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
	  

    ELSE  




    END IF

    

  
  
  
  
  
  



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


  IMPLICIT NONE  



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




  INTEGER         :: i,j,k
  INTEGER :: i_start, i_end, j_start, j_end, i_endu, j_endv








  i_start = its
  i_end   = min(ite,ide-1)
  j_start = jts
  j_end   = min(jte,jde-1)

  i_endu = ite
  j_endv = jte






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

   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: fzm,  &     
                                                                  fzp,  &     
                                                                  zh

   REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy,  &
																  dts,  &
																  ztop
																  
  
   
   
   
 
   
    REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) :: tendency 
	


 

      INTEGER :: i,j,k, i_start, i_end, j_start, j_end, k_start, k_end
      INTEGER :: i_endu, j_endv, k_endw
      INTEGER :: i_start_up, i_end_up, j_start_up, j_end_up
      INTEGER :: i_start_vp, i_end_vp, j_start_vp, j_end_vp
      INTEGER :: i_start_u_tend, i_end_u_tend, j_start_v_tend, j_end_v_tend









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
	
	
	
	


    
	
	
	
	
	
  u_outer_j_loop2: DO j = j_start, j_end
   DO k = k_start, k_end
     DO i = i_start_up, i_end_up
	 
	 
       u(i,k,j)=u(i,k,j) + dts*(alt(i-1,k,j)+alt(i,k,j))*cqu(i,k,j)*0.5*tendency(i,k,j)
	   
	   
	   
     ENDDO
   ENDDO
   ENDDO u_outer_j_loop2


		
    tendency = 0.0 
	
	
	
	call    calcu_grady( p, tendency,                 &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )   
   

 
	   


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







						
SUBROUTINE calc_p_rho_fm( al, p, pfm,  pb,                  &
                       alt, rho , rhob, pm1,      &
                       smdiv,             &
                       non_hydrostatic, step,        &
                       ids, ide, jds, jde, kds, kde, &
                       ims, ime, jms, jme, kms, kme, &
                       its,ite, jts,jte, kts,kte    )

  IMPLICIT NONE  



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



  INTEGER :: i, j, k
  INTEGER :: i_start, i_end, j_start, j_end, k_start, k_end
  REAL    :: ptmp



















   i_start = its
   i_end   = min(ite,ide-1)
   j_start = jts
   j_end   = min(jte,jde-1)
   k_start = kts
   k_end = min(kte,kde-1)


     DO j=j_start, j_end
     DO k=k_start, k_end
     DO i=i_start, i_end



      al(i,k,j) = rho(i,k,j) - rhob(i,k,j)



      p(i,k,j) = pfm(i,k,j) - pb(i,k,j)
	  alt(i,k,j) = 1.0/rho(i,k,j)  
	  

     ENDDO
     ENDDO
     ENDDO


END SUBROUTINE calc_p_rho_fm





SUBROUTINE calc_p_rho_fm_new ( moist, n_moist,      &
                            al, alt,p,rhofm,rhofm_base,pfm,pfm_base,    &
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            its, ite, jts, jte, kts, kte          )

  IMPLICIT NONE
   
   



  INTEGER ,          INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                      ims, ime, jms, jme, kms, kme, &
                                      its, ite, jts, jte, kts, kte

  INTEGER ,          INTENT(IN   ) :: n_moist
 

  REAL, DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(IN   ) :: rhofm,rhofm_base,pfm,pfm_base 


  REAL, DIMENSION( ims:ime , kms:kme , jms:jme, n_moist ), INTENT(IN   ) :: moist

  REAL, DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(  OUT) :: al, alt,p


  

  INTEGER :: i, j, k, kk, itf, jtf, ktf, ispe
  REAL    :: qtot, qf1, qf2, qvf
  
  real :: rhot
  
  REAL, DIMENSION( its:ite) :: temp,cpovcv_v
 














  itf=MIN(ite,ide-1)
  jtf=MIN(jte,jde-1)
  ktf=MIN(kte,kde-1)

  cpovcv_v = cpovcv


          al = 0.0 
		 





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
			
			
			rhot   = rhofm(i,k,j)*(1.0+qtot) 
		 
			
            alt(i,k,j) = 1.0/rhofm(i,k,j)  
			
			
			
          	
            
   
		
		
		
		     al(i,k,j) = rhot  - rhofm_base(i,k,j)  
			
			
	
     
	 
	 
	 
	 
	 
	 
	 
	 
			
			

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
      
               p(i,k,j) = pfm(i,k,j) - pfm_base(i,k,j)   
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
   
   
   

   
   
   TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags
   
   INTEGER,      INTENT(IN   )    :: spec_zone
   
   INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte


  
   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN) :: w
    REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN  ) :: rho,p
	
	
   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: u,v

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: w_tend, drho,dp

    
   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) ::  alt,cqw,cqx,rho_save,p_save, rhob,pfmb
   
  
   REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft, &
																	zs

   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: fzm,  &     
                                                                  fzp, dnw,  &     
                                                                  zh

   REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy,  &
																  dts,  &
																  ztop
																  
  
     REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  ::  sigmdot

     REAL , DIMENSION( ims:ime , jms:jme ) :: rmm,rmsft,zs1 

     REAL , DIMENSION( kms:kme ) :: sgmh,dsgm,ax,bx 

 

   
   
 
 

      INTEGER :: i,j,k, i_start, i_end, j_start, j_end, k_start, k_end
	  
	
	
	
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



      
       gama = rcv+1



       
       epssm = 1.0   
	 
       alfa = 0.5*(1.-epssm)
       beta = 0.5*(1.+epssm)







      i_start = its
      i_end   = min(ite,ide-1)
      j_start = jts
      j_end   = min(jte,jde-1)
	  
      k_start = kts
      k_end   = kte-1
	  


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


 
	 rhox = rho_save*(1.+cqx )  
	 
   if(.false.) then 


     




 
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
    
    
      k_start = kts
      k_end   = MIN(kte,kde)
   
   

    DO  j = j_start, j_end
    DO  k = kts+1,kde-1  
    DO  i = i_start, i_end
       wb(i,k,j) =   -1.*beta*dts*altw(i,k,j)*g
    enddo
    enddo
    enddo

 



     pc = 0.0
	 pd = 0.0 
	 
   pc = 0.5*beta*dts*(1.0-gama)*rhob*g
   

   
   
   
  
  
  
  
  
    DO  j = j_start, j_end
    DO  k = kts+1,kde-1  
    DO  i = i_start, i_end
	
       pd(i,k,j)= 0.5*beta*dts*(gama-1.0)*rho_save(i,k,j)* ( 1.+ cqx(i,k,j)      ) 

	 enddo
    enddo
    enddo
   
  
	 
       pe = 0.0 



						
						
   
    
	  
		
  
		
      DO  j = j_start, j_end
      DO  k = kts   , kte
      DO  i = i_start, i_end
            pe(i,k,j) =  -1.0 * altw(i,k,j)*pz(i,k,j)    
      ENDDO
      ENDDO 
      ENDDO


 
	   



  
  

 
  



      k_end = MIN(kte,kde-1)

   
      

	  ra = 0.0 
	  pa = 0.0
	  
    DO  j = j_start, j_end
    DO  k = k_start, k_end
    DO  i = i_start, i_end

     ra(i,k,j) = dts*msft(i,j)*msft(i,j)* beta/dsgm(k)
  
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
						
						

    DO  j = j_start, j_end
    DO  k = kts, kte
    DO  i = i_start, i_end
	 
	
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


     





		  


   





  
END SUBROUTINE calcu_abc_z


SUBROUTINE sumflux ( ru, rv, ww, rho,                            &
                     u_lin, v_lin, ww_lin,rho_lin,                   &
                   
                   
                   
                     ru_m, rv_m, ww_m, rho_m,epssm,                &
                  
                     iteration , number_of_small_timesteps,  &
                     ids,ide, jds,jde, kds,kde,              &
                     ims,ime, jms,jme, kms,kme,              &
                     its,ite, jts,jte, kts,kte              )


  IMPLICIT NONE  



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




 
 

  INTEGER :: mini, minj, mink


  REAL, INTENT(IN   )  ::  epssm
  INTEGER   :: i,j,k









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
  
				      + u_lin(i,k,j)
					  
      rv_m(i,k,j)  = rv_m(i,k,j) / number_of_small_timesteps   &
                     + v_lin(i,k,j)

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
                        u_save,v_save,w_save,z,engk1,                                             &   
                        config_flags,spec_zone,               &
                        msft,zs,                              &
                        zh,fzm, fzp,dnw,                          &
                        rdx, rdy, dts, ztop, z_at_w,                 &
                        ids, ide, jds, jde, kds, kde,         &
                        ims, ime, jms, jme, kms, kme,         &
                        its, ite, jts, jte, kts, kte  )

   IMPLICIT NONE
   
   
   

   
   
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

    
   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) ::  alt,cqw,cqx,rho_save,p_save, rhob,pfmb
   
  
   REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft, &
																	zs

   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: fzm,  &     
                                                                  fzp,dnw,  &     
                                                                  zh

   REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy,  &
																  dts,  &
																  ztop
																  
																  
	 REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) ::  z_at_w															  
																  
  
     REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  ::  sigmdot

     REAL , DIMENSION( ims:ime , jms:jme ) :: rmm,rmsft,zs1 

     REAL , DIMENSION( kms:kme ) :: sgmh,dsgm,ax,bx 



 

   
   
 
 

      INTEGER :: i,j,k, i_start, i_end, j_start, j_end, k_start, k_end
	  
	
	
	
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
	   
	   
	  
	  REAL, DIMENSION (kts:kte) :: dampwt
      real :: htop,hbot,hdepth,hk
      real    :: pi,dampmag
      REAL :: muthk, muthkm1


      REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: u_save,v_save,w_save,z,engk1
       REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  ::  dengk,engk,engk2,letot
	   REAL :: co1,co2,co3





      
       gama = rcv+1



       
       epssm = 1.0   
	 
       alfa = 0.5*(1.-epssm)
       beta = 0.5*(1.+epssm)







      i_start = its
      i_end   = min(ite,ide-1)
      j_start = jts
      j_end   = min(jte,jde-1)
	  
      k_start = kts
      k_end   = kte-1
	  
	  
	  
	  
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
						
						
 

 
    call calcu_divg(p_save,dpf,           & 
                       u, v, sigmdot,                       &
                   config_flags,                  &
                      msft,zh,zs,fzm,fzp,                    &
                       rdx,rdy,ztop,                  &
                       ids, ide, jds, jde, kds, kde,  &
                     ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )
						
 
 
 
 
  
  
   
   



	dpf1 = 0.0 
	
 	  call     calculate_ptend_part2_fm_small_test( p, dpf1,                 &

                        u, v, w2,sigmdot,                       &
						rho_save,rhob,         &
						alt,     &
                     config_flags,                  &
                        msft,zs,zh,fzm,fzp,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )



 

  
   
   
	 dpf= cpovcv*dpf + dpf1					
						
 	   
       


     dwsc =  w_tend  + alfa*dwf 
	 
	 
	 
	      
          
		 

		 
 
	    
	    
		
     	
	   
     dpsc = dp + dpf
	 
     drhosc = drho + drhof


	  

     ws = dts*dwsc + w   
     ps = dts*dpsc + p
     rs = dts*drhosc+rho
	 

   


   
						
	

		
      IF ( config_flags%periodic_y ) j_start = jts
      IF ( config_flags%periodic_y ) j_end = jte 
	
       


	

  call    calculate_sgmdot2( sigmdot,                 &     
                     u, v, w2,                       &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  ) 
						
if(.false.) then 
       CALL set_w_surface_fm_small( config_flags,              &
                           w, zs,  u, v,          &
                            rdx,  rdy,                  &
                            msft,                          &
                              ids, ide, jds, jde, kds, kde,                    &
                              ims, ime, jms, jme, kms, kme,                    &
                              its, ite, jts, jte, kts, kte                     )
    
    
    
    
endif
	  
    DO  j = j_start, j_end
 
    DO  i = i_start, i_end
  
	   
    
	   w(i,1,j) = -sigmdot(i,1,j)/(ztop/(ztop-zs(i,j)))/beta  
      w(i,kte,j) = 0.
    ENDDO
    ENDDO


 

   
    
      k_start = kts+1
      k_end   = MIN(kte,kde-1)
    
    
 

    
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


   if(.false.) then    
   
   
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
						
 
    drhof = 0.0  

		call calcu_divg( rho_save, drhof,                 &
                        u, v, sigmdot,                       &
                     config_flags,                  &
                        msft,zh,zs,fzm,fzp,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )
						
						
	
	  drhosc= drho+drhof 

   
	  
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

    	 engk2 =(rho_save+rho)*engk2 
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
   
   
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte
											  
											  
		REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: pp

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: u,     &
                                                                      v, w,&
                                                                      sgmdot
																	  
																	  
	    REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: rho,rhob,alt  
		

		 
		 
																		 
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT( INOUT) :: tendency

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,zs


         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh,fzm,fzp

         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   
   
   
   
         INTEGER :: i, j, k, ktf
         INTEGER :: i_start, i_end, j_start, j_end
      
      
       
       
         

	   
	   REAL , DIMENSION( ims:ime , kms:kme , jms:jme )   :: px,py,pz,altw

      
			

    LOGICAL :: specified
	
     
	 


 
   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.



 
   
 	   px = 0.0
	   py = 0.0 
	   pz = 0.0
	   


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
                      	 - rhob(i,k,j) * g * 0.5* (  w(i,k,j) + w(i,k+1,j)   )           &  				   
	                                          )

      ENDDO
      ENDDO

      ENDDO 

  


	
   


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
   
   
   
 
if(.false.) then 
   
   
  

   


   i_start = its
   i_end=MIN(ite,ide-1)

   

   IF ( (config_flags%open_ys) .and. jts == jds ) THEN

     j=jts

       DO k = kts, kte-1
     DO i = i_start, i_end
	   
	   
	   
     
      
      
      
							  
							  
							  
        tendency(i,k,j) =  tendency(i,k,j)  +      &
	 
	       rcv*rho(i,k,j)* 0.5* v(i,k,j+1) * ( alt(i,k,j+1) + alt(i,k,j) )   * py(i,k,j+1)    

 







	 		   
							  
       ENDDO
     ENDDO

   END IF

   


   

   IF ( (config_flags%open_ye) .and. jte == jde ) THEN

     j=jte-1

       DO k = kts, kte-1 
     DO i = i_start, i_end
	   
	   
       
        
			   
		
			   
       
		
        tendency(i,k,j) =  tendency(i,k,j)  +      &
	 	       rcv*rho(i,k,j)* 0.5*   v(i,k,j) *  ( alt(i,k,j)  + alt(i,k,j-1) ) * py(i,k,j) 
       ENDDO
     ENDDO

   END IF
   
   
   
   
   
      

   j_start = jts
   j_end=MIN(jte,jde-1)

   

   IF ( (config_flags%open_xs) .and. its == ids ) THEN

     i=its

	   DO j = j_start, j_end
     DO k = kts, kte-1 
	 
       
       

	  tendency(i,k,j) =  tendency(i,k,j)  +      &
	 	       rcv*rho(i,k,j)* 0.5* u(i+1,k,j) * ( alt(i+1,k,j) + alt(i,k,j) )   * px(i+1,k,j)   

       ENDDO

        
     ENDDO

   END IF

   

   IF ( (config_flags%open_xe) .and. ite == ide ) THEN

     i = ite-1

     DO j = j_start, j_end
     DO k = kts, kte-1 
	 
       
			   
      


	  tendency(i,k,j) =  tendency(i,k,j)  +      &
	 	       rcv*rho(i,k,j)* 0.5* u(i,k,j)  *  ( alt(i,k,j) + alt(i-1,k,j) )  * px(i,k,j)    

       ENDDO



     ENDDO

   END IF




 

 endif

 
 

END SUBROUTINE calculate_ptend_part2_fm_small_test




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
   
   
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte
											  
											  
		REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: pp

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: u,     &
                                                                      v, w,&
                                                                      sgmdot
																	  
																	  
	    REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: rho,rhob,alt  
		

		 
		 
																		 
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT( INOUT) :: tendency

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,zs


         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh,fzm,fzp

         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   
   
   
   
         INTEGER :: i, j, k, ktf
         INTEGER :: i_start, i_end, j_start, j_end
      
      
       
       
         

	   
	   REAL , DIMENSION( ims:ime , kms:kme , jms:jme )   :: px,py,pz,altw

      
			

    LOGICAL :: specified
	
     
	 


 
   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.



 
   
 	   px = 0.0
	   py = 0.0 
	   pz = 0.0
	   
	  
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
        
	                                          )

      ENDDO
      ENDDO

      ENDDO 

  


	
   


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


SUBROUTINE calc_p_damp( p, pm1,      &
                       smdiv,             &
                       step,        &
                       ids, ide, jds, jde, kds, kde, &
                       ims, ime, jms, jme, kms, kme, &
                       its,ite, jts,jte, kts,kte    )

  IMPLICIT NONE  



  INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
  INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
  INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte

  INTEGER,      INTENT(IN   )    :: step



  REAL, DIMENSION(ims:ime, kms:kme, jms:jme),INTENT(INOUT) :: p,pm1


  REAL,                                       INTENT(IN   ) :: smdiv





  INTEGER :: i, j, k
  INTEGER :: i_start, i_end, j_start, j_end, k_start, k_end
  REAL    :: ptmp



















   i_start = its
   i_end   = min(ite,ide-1)
   j_start = jts
   j_end   = min(jte,jde-1)
   k_start = kts
   k_end = min(kte,kde-1)

     

     IF (step == 0) then   
       DO j=j_start, j_end
       DO k=k_start, k_end
       DO i=i_start, i_end
         pm1(i,k,j)=p(i,k,j)
       ENDDO
       ENDDO
       ENDDO
     ELSE                     
       DO j=j_start, j_end    
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

     REAL ,       INTENT(IN   )  :: rdx, rdy   

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
















     jm1_limit = jds        
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
