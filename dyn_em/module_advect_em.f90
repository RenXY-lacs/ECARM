



MODULE module_advect_em

  USE module_bc
  USE module_model_constants
  USE module_wrf_error

CONTAINS



     




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
	  
	  
  
    i_start = MAX(its-1,ids)
    i_end   = MIN(ite+1,ide-1)
    j_start = MAX(jts-1,jds)
    j_end   = MIN(jte+1,jde-1) 
	  

	   
	   
	   
	
    

    
    
	  
	  
    
    
    
    
	  
       engk  = 0.0 
	   
	   
      DO      j = j_start, j_end
      DO  k = k_start, k_end
      DO  i = i_start, i_end
          
			  engk(i,k,j) = ( u(i,k,j)**2 + u(i+1,k,j)**2 + v(i,k,j)**2 + v(i,k,j+1)**2 + w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25 
      ENDDO
      ENDDO
      ENDDO
	  

if(.false.) then 


	  
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
	  
	
endif 
	
    END SUBROUTINE calcu_engk




	  
	  


	  
     SUBROUTINE calculate_sgmdot( sgmdot,                 &
                        u, v, w,                       &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) ::     u,     &
                                                                      v, &
                                                                      w
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) :: sgmdot

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
		 
	     INTEGER :: i_endu, j_endv
				 
				 
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         
         
         REAL , DIMENSION( kms:kme ) :: sgmh,dsgm,dsgmh
      
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: sx,sy 
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
	  
	  
	   sgmdot = 0.0 






 
		  
		  sx = 0.0
		  sy = 0.0 
	

      i_start = its  
      i_end   = MIN(ide-1,ite)
	  
	  j_start = jts
      j_end   = MIN(jte,jde-1)
	  
      DO  j = j_start, j_end
      DO  k = kts, ktf
      DO  i = i_start, i_endu 
	  
	  
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



   
   


      IF ( config_flags%periodic_x ) i_start = its
      IF ( config_flags%periodic_x ) i_end = ite


    
    

	  


      i_start = its
      i_end   = MIN(ide-1,ite)
	  
	  j_start = jts
      j_end   = MIN(jte,jde-1)	  
	  
	  
     DO  j = j_start, j_end
      DO  k = kts+1, ktf
      DO  i = i_start, i_end 
	  
    	
	  
	     
	  
          sxxs = 0.5*((sx(i+1,k-1,j) + sx(i,k-1,j))*dsgm(k) +(sx(i+1,k,j) + sx(i,k,j))*dsgm(k-1) ) /(dsgm(k-1) + dsgm(k) ) 
	  
	      syys = 0.5*((sy(i,k-1,j+1) + sy(i,k-1,j))*dsgm(k) +(sy(i,k,j+1) + sy(i,k,j))*dsgm(k-1) ) /(dsgm(k-1) + dsgm(k) ) 
	  
	  
	      sgmdot(i,k,j) = - sxxs - syys + ztop/(ztop-zs(i,j))*w(i,k,j)
	  

		
	  ENDDO
      ENDDO
      ENDDO

  
 
  
	  

      DO  j = j_start, j_end
      DO  i = i_start, i_end 	  
	  
	        sgmdot(i,kts,j) =  0. 
	       sgmdot(i,kte,j) =  0. 
		   
		   
      ENDDO
      ENDDO 


	

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
   
   
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) ::     u,     &
                                                                      v, &
                                                                      w
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) :: sgmdot

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         
         
         REAL , DIMENSION( kms:kme ) :: sgmh,dsgm,dsgmh
      
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: sx,sy 
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
   
   
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var,u,     &
                                                                      v, &
                                                                      sgmdot
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: divg

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh,fzm,fzp


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx
         
        
       
	   
	   
	   
	   
	   
		 
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


 



	  

   
      CALL set_physical_bc3d( var, 'p', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              its, ite, jts, jte, kts, kte,       &
                              its, ite, jts, jte, kts, kte )
	  
	  

     divg = 0.0 
	 

 

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
  
    
    

      	divg(i,k,j) =  -  (  0.5*((ztop-zs(i,j))*var(i,k,j)/msft(i,j)  + (ztop-zs(i+1,j))*var(i+1,k,j)/msft(i+1,j) ) *u(i+1,k,j) &
                     -   0.5*((ztop-zs(i-1,j))*var(i-1,k,j)/msft(i-1,j)  + (ztop-zs(i,j))*var(i,k,j)/msft(i,j) ) *u(i,k,j) )*rdx
      ENDDO
      ENDDO
      ENDDO
	  
	  
	IF ( config_flags%open_xs  .and.  its == ids ) THEN 
		
	
		
      DO  j = j_start, j_end
      DO k=kts,ktf
  
  
 
     
     


					 
	
	         ub = MIN( 0.5*(u(its,k,j)+u(its+1,k,j)), 0. )
		 
         divg(its,k,j) =                      &
               - (ztop-zs(its,j))/msft(its,j)*rdx*(                                             &
                       ub*(   var(its+1,k,j)                 &
                            - var(its  ,k,j)   ) +           &
                       var(its,k,j)*(u(its+1,k,j)-u(its,k,j))  &       
                                                                )				 
					 
      ENDDO
      ENDDO
 
	  
	ENDIF	  
	

	  
	    IF ( config_flags%open_xe  .and.  ite == ide ) THEN 
		
	 
		
		 i_end   = MIN(ite,ide-1)
		 
		
      DO  j = j_start, j_end
      DO k=kts,ktf
 
  
      
      
	
	     ub = MAX( 0.5*(u(ite-1,k,j)+u(ite,k,j)), 0. )	 
		 
         divg(i_end,k,j) =                  &
               - (ztop-zs(i_end,j))/msft(i_end,j)*rdx*(                                               &
                       ub*(  var(i_end  ,k,j)                  &
                           - var(i_end-1,k,j) ) +              &
                       var(i_end,k,j)*(u(ite,k,j)-u(ite-1,k,j))  &   
                                                                    )				 
					 
					 
      ENDDO
      ENDDO
  

	ENDIF	
	  
	  
	  
	  



      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
      
      IF ( config_flags%open_ys .or. specified .or. config_flags%polar ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified .or. config_flags%polar ) j_end   = MIN(jde-2,jte) 	  

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end

   
   
					 
 		 
    	divg(i,k,j) = divg(i,k,j) -   ( 0.5*((ztop-zs(i,j))*var(i,k,j)/msft(i,j)  + (ztop-zs(i,j+1))*var(i,k,j+1)/msft(i,j+1) ) *v(i,k,j+1)    &
                 -   0.5*((ztop-zs(i,j-1))*var(i,k,j-1)/msft(i,j-1)  + (ztop-zs(i,j))*var(i,k,j)/msft(i,j) ) *v(i,k,j) )*rdy					 

   

       ENDDO
      ENDDO
      ENDDO
	  
	  

	  
	    IF ( config_flags%open_ys  .and.  jts == jds ) THEN 
		
	
		
      DO k=kts,ktf
      DO  i = i_start, i_end		
		
    
    
	
	
	
	      vb = MIN( 0.5*(v(i,k,jts)+v(i,k,jts+1)), 0. )
		  
        divg(i,k,jts) =    divg(i,k,jts)                &
               - (ztop-zs(i,jts))/msft(i,jts)*rdy*(                                             &
                       vb*(  var(i,k,jts+1)                  &
                           - var(i,k,jts  ) ) +              &
                       var(i,k,jts)*(v(i,k,jts+1)-v(i,k,jts))  &     
                                                               )
					 
      ENDDO
      ENDDO
	    ENDIF	
		
 
	  
	    IF ( config_flags%open_ye  .and.  jte == jde ) THEN 
		
	 
		
		j_end   = MIN(jte,jde-1)
		
      DO k=kts,ktf
      DO  i = i_start, i_end		
		
	  
       
	
		  vb = MAX( 0.5*(v(i,k,jte-1)+v(i,k,jte)), 0. ) 
		 
         divg(i,k,j_end) =    divg(i,k,j_end)               &
               -  (ztop-zs(i,j_end))/msft(i,j_end)*rdy*(                                               &
                       vb*(   var(i,k,j_end  )                 &
                            - var(i,k,j_end-1) ) +             &
                       var(i,k,j_end)*(v(i,k,jte)-v(i,k,jte-1))  &    
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








			
   
      SUBROUTINE calcu_vorty( vorty,                 &
                        u, v, w,                       &
                     config_flags,                  &
                        msft,msfty,xlat,zh,zs,fzm,fzp,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: u,     &
                                                                      v,w
																	  
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: vorty

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,msfty,xlat,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh,fzm,fzp


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx
  
       
	   
       

       REAL , DIMENSION( kms:kme ) :: sgmh,dsgm,dsgmh
       REAL , DIMENSION( ims:ime , jms:jme )  ::  rmsft,rmm
	   
	   REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  ::  uw,uwy,vds,vss,udxzs,tmp,wvdsh 

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme )   :: vorty1,vorty2,vorty21, vorty22
		 
		 
       LOGICAL :: specified
	   
	   integer,dimension(3) ::  fdim(3)
	  
 
	   vorty = 0.0 
	   vorty1 = 0.0
	   vorty2 = 0.0
	   vorty21 = 0.0 
	   vorty22 = 0.0 


   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

 
	

	 
       ktf=MIN(kte,kde-1)

	  DO k=kts,kde-1
          sgmh(k) = 0.5*(zh(k)+zh(k+1))
          dsgm(k) = zh(k+1) - zh(k) 
      ENDDO

 	  DO k=2,kde-1
             dsgmh(k) = 0.5*(dsgm(k) + dsgm(k-1))   
			 
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

  




	uwy = 0. 
	
		
	
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
		  


      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = jte

    
     
     
     
	 
	 

	
	
   
    
	
       IF ( config_flags%periodic_y ) j_start = jts           
      IF ( config_flags%periodic_y ) j_end = jte	  
	  

	  
	  

      IF ( config_flags%open_ys .or. specified ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified ) j_end   = MIN(jde-1,jte)
	  
	
    
	  
 

      IF ( config_flags%open_xs  .or. specified ) i_start = MAX(ids+1,its)   
      IF ( config_flags%open_xe  .or. specified ) i_end   = MIN(ide-1,ite)
	  

	  
		  
      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
    
         
         
		 
		 
             vorty21(i,k,j) =   0.5*( (u(i,k,j)+u(i,k,j-1))*0.5 * (v(i,k,j) - v(i-1,k,j))*rdx   &
                                                   +  (u(i+1,k,j)+u(i+1,k,j-1))*0.5 * (v(i+1,k,j) - v(i,k,j))*rdx  )       

      ENDDO
      ENDDO
      ENDDO
	  
	  
    
     
	
	
     
       
		
		
	  
      
	  
	  
	  


	   
    
     
     
     
     

     
    
	  
     
    
    
    
		
    



    
    

    
    
	  


      DO  j = j_start, j_end
      DO k=2,ktf 
      DO  i = i_start, i_end  
    
      	vds(i,k,j) = ( v(i,k,j) - v(i,k-1,j) )/dsgmh(k) 
	
		
      ENDDO
      ENDDO
      ENDDO	 
	  
	  
	 
	  DO  j = j_start, j_end
      DO  i = i_start, i_end

    
      	vds(i,kte,j) = vds(i,kte-1 ,j) 
        vds(i,1,j) =  2.0*( v(i,1,j) - 0.0 )/dsgmh(2)   
	
		
      ENDDO
      ENDDO
	  
	       

	   
	   vss = 0.0 

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
     
      	
          vss(i,k,j) = 0.5*( vds(i,k,j) + vds(i,k+1,j) ) 
      ENDDO
      ENDDO
      ENDDO
	  

	  
 


 
	  
	udxzs = 0.0 
	
	
			j_end =  jte+1   
	  
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
							  
							  


				j_end =  jte   

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
     
      	vorty22(i,k,j)=    - 0.5*( udxzs(i,k,j)  + udxzs(i,k,j+1)  ) * vss(i,k,j) 

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
	  
	  
	  
	  
	  
	  



   
   
	  

	  
	  

     i_start = its
     i_end   =  MIN(ite,ide-1)
	 
	 
   
   

	  
    
    

	  
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



	  

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
 
      	vorty(i,k,j) = vorty1(i,k,j) + vorty2(i,k,j)  -   (ztop/(ztop-zs(i,j)) + ztop/(ztop-zs(i,j-1)))*0.5*wvdsh(i,k,j)

      ENDDO
      ENDDO
      ENDDO	  
	  


		  
      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end

      
      
      
      


       
                
                
                

  
        vorty(i,k,j) = vorty(i,k,j)  - 0.25* ( u(i+1,k,j)*u(i+1,k,j)   + u(i,k,j)*u(i,k,j)  + u(i,k,j-1)*u(i,k,j-1)  + u(i,k,j-1)*u(i,k,j-1)  )  & 
                      *(msft(i,j)-msft(i,j-1))*rdy    &
                                        + v(i,k,j) * ( u(i,k,j)*(msft(i,j)-msft(i-1,j))*rdx   &
                                                   +  u(i,k,j-1)*(msft(i,j-1)-msft(i-1,j-1))*rdx   &
                                                   +  u(i+1,k,j-1)*(msft(i+1,j-1)-msft(i,j-1))*rdx    &
                                                   +  u(i+1,k,j)*(msft(i+1,j)-msft(i,j))*rdx  ) * 0.25 
     



      ENDDO
      ENDDO
      ENDDO	  

  
      if(.true.)then   
	  
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






      

   
      SUBROUTINE calcu_vortz(  vortz,                 &
                        u, v, w,                      &
                     config_flags,                  &
                        msft,xlat,zh,zs,fzm,fzp,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: u,v,w
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: vortz

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,xlat,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh,fzm,fzp


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop
		 
		 
		    REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: ufull

   
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx
  


       REAL , DIMENSION( kms:kme ) :: sgmh,dsgm,dsgmh
       REAL , DIMENSION( ims:ime , jms:jme )  ::  rmsft,rmm
	   
       
	   REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  ::  uv,vortz1,vortz21,udxzs,temp,vortz22,vdyzs,vortz31,vortz32,vortz4


       LOGICAL :: specified


   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

 
      ktf=MIN(kte,kde-1)
	  
 

	  DO k=kts,kde-1
          sgmh(k) = 0.5*(zh(k)+zh(k+1))
          dsgm(k) = zh(k+1) - zh(k) 
      ENDDO

 	  DO k=2,kde-1
             dsgmh(k) = 0.5*(dsgm(k) + dsgm(k-1))   
			 
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
	  


	 

		vortz1 = 0.0 
		
        call  calcu_gradz( uv, vortz1,                    &
                      config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )





		 	call average_vertical_w_x( ufull, u,                 &
                     config_flags,                  &
                         fzm,fzp,                    &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )	

  
						
      CALL set_physical_bc3d( ufull, 'u', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              its, ite, jts, jte, kts, kte,       &
                              its, ite, jts, jte, kts, kte )	


		
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
        vortz22(i,kts,j) = vortz22(i,kts+1,j)      
        vortz22(i,kte,j) = vortz22(i,ktf,j) 
		
      ENDDO
		
      ENDDO


      

       
	   
     vortz31 = 0.0 
      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
    
            vortz31(i,k,j) =  0.5*((v(i,k,j)+v(i,k+1,j))*0.5 * (w(i,k,j)-w(i,k,j-1))*rdy                           &
                                                   + (v(i,k,j+1)+v(i,k+1,j+1))*0.5 * (w(i,k,j+1)-w(i,k,j))*rdy )
        

      ENDDO
      ENDDO
      ENDDO
	  


      
	  vdyzs = 0.0 
	  
      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
      

      
    
            vdyzs(i,k,j) =  (ztop-sgmh(k))/(ztop-zs(i,j))* ( v(i,k,j)* ( zs(i,j)- zs(i,j-1) )*rdy +  v(i,k,j+1)* ( zs(i,j+1)- zs(i,j) )*rdy)*0.5   &
                              * (w(i,k+1,j)-w(i,k,j))/dsgm(k)
							  
							  


      ENDDO
      ENDDO
      ENDDO
	  

			

			

      
	  
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
	
      

			
			
	  vortz4 = 0.0 		
			
      DO  j = j_start, j_end
      DO k=kts+1,ktf
      DO  i = i_start, i_end
            vortz4(i,k,j) =  fzm(k)*uv(i,k,j)*2.0  + fzp(k)*uv(i,k-1,j)*2.0   
      ENDDO
      ENDDO
	  
      DO  i = i_start, i_end
        vortz4(i,kts,j) =  vortz4(i,kts+1,j)
        vortz4(i,kte,j) =  vortz4(i,ktf,j) 
      ENDDO
	  
      ENDDO


	  
	  
	  
 	
       vortz = 0.0 	
	  
      DO  j = j_start, j_end
      DO k=kts,ktf+1 
      DO  i = i_start, i_end
 
      	vortz(i,k,j) =  -1.0*vortz1(i,k,j) -msft(i,j)*(vortz21(i,k,j)-vortz22(i,k,j) + vortz31(i,k,j) - vortz32(i,k,j)  )   &    
			  + reradius*vortz4(i,k,j)         

      ENDDO
      ENDDO
      ENDDO	  
	  



      END SUBROUTINE calcu_vortz







	
	

	  	  

   
      SUBROUTINE calcu_vortx( vortx,                 &
                        u, v, w,                    &
                     config_flags,                  &
                        msft,xlat,zh,zs,fzm,fzp,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) ::u,v,w
                                                                     
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT( OUT) :: vortx

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,xlat,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh,fzm,fzp


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx

         
        
       
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: vortx1,vortx21,vortx22 


      REAL , DIMENSION( kms:kme ) :: sgmh,dsgm,dsgmh
       REAL , DIMENSION( ims:ime , jms:jme )  ::  rmsft,rmm
	   
	   REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  ::  vw,us,tmp,vwx,uds,wudsh

    LOGICAL :: specified




   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

 
      ktf=MIN(kte,kde-1)
	  
 
	 
	  DO k=kts,kde-1
          sgmh(k) = 0.5*(zh(k)+zh(k+1))
          dsgm(k) = zh(k+1) - zh(k) 
      ENDDO

 	  DO k=2,kde-1
             dsgmh(k) = 0.5*(dsgm(k) + dsgm(k-1))   
			 
      ENDDO
	  
	 
       
      i_start = its
      i_end   =  MIN(ite,ide-1)

      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
	  
if(.false.) then   

      DO      j = j_start, j_end
      DO  k = kts, ktf
      DO  i = i_start, i_end

      
     
      
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
	
      call  calcu_gradx( vw, vortx1,                 &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )


	  



      i_start = its
      i_end   = ite  


      IF ( config_flags%open_xs ) i_start = MAX(ids+1,its)   
      IF ( config_flags%open_xe ) i_end   = MIN(ide-1,ite)
	  
	  
    
    
	
       IF ( specified ) i_start = MAX(ids+2,its)            
       IF ( specified ) i_end   = MIN(ide-2,ite)
	
	  
      IF ( config_flags%periodic_x ) i_start = its           
      IF ( config_flags%periodic_x ) i_end = ite


      IF ( config_flags%open_ys .or. specified ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified ) j_end   = MIN(jde-1,jte)
	  
	  



     vortx21  = 0.0 
         DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end

    
    
      	vortx21(i,k,j) =  0.5*( (v(i-1,k,j+1)+ v(i,k,j+1))*0.5 * (u(i,k,j+1) - u(i,k,j) )*rdy   &
                             +  (v(i-1,k,j)  + v(i,k,j) )*0.5  * (u(i,k,j)   - u(i,k,j-1) )*rdy )

        
        
        

        
        
        


        
         

         
      ENDDO
      ENDDO
      ENDDO







	   vortx22  = 0.0 
	  
      DO  j = j_start, j_end
	  
	  
      DO k=2,ktf-1
      DO  i = i_start, i_end

        
        

           


           vortx22(i,k,j) =  - 0.5*((ztop-sgmh(k))/(ztop-zs(i,j) ) *0.5*(v(i,k,j+1)  *(zs(i,j+1)   - zs(i,j)  )*rdy + v(i,k,j)  *(zs(i,j)   - zs(i,j-1)   )*rdy)   &
                              +    (ztop-sgmh(k))/(ztop-zs(i-1,j)) *0.5*(v(i-1,k,j+1)*(zs(i-1,j+1) - zs(i-1,j))*rdy + v(i-1,k,j)*(zs(i-1,j) - zs(i-1,j-1) )*rdy)   )   &
                   *  0.5*  ( (u(i,k,j)   - u(i,k-1,j) )/dsgmh(k)  & 
				            + (u(i,k+1,j) - u(i,k,j) )/dsgmh(k+1) )    



      ENDDO
      ENDDO
	
	   k = kts
	   DO  i = i_start, i_end  
	              vortx22(i,k,j) =  - 0.5*(   &
                  (ztop-sgmh(k))/(ztop-zs(i,j))   *0.5*(v(i,k,j+1)*(zs(i,j+1)     - zs(i,j))  *rdy + v(i,k,j)  *(zs(i,j)   - zs(i,j-1))*rdy)   &
             +    (ztop-sgmh(k))/(ztop-zs(i-1,j)) *0.5*(v(i-1,k,j+1)*(zs(i-1,j+1) - zs(i-1,j))*rdy + v(i-1,k,j)*(zs(i-1,j) - zs(i-1,j-1))*rdy) &
               )   &
                   *  0.5*  (2.0*(u(i,k,j) - 0.  )/dsgm(kts)    &
                               + (u(i,k+1,j) - u(i,k,j) )/dsgmh(k+1) )    

	    enddo 

	
	   k = ktf
	   DO  i = i_start, i_end  
           vortx22(i,k,j) =  - 0.5*((ztop-sgmh(k))/(ztop-zs(i,j) ) *0.5*(v(i,k,j+1)  *(zs(i,j+1)   - zs(i,j)  )*rdy + v(i,k,j)  *(zs(i,j)   - zs(i,j-1)   )*rdy)   &
                              +    (ztop-sgmh(k))/(ztop-zs(i-1,j)) *0.5*(v(i-1,k,j+1)*(zs(i-1,j+1) - zs(i-1,j))*rdy + v(i-1,k,j)*(zs(i-1,j) - zs(i-1,j-1) )*rdy)   )   &
                   *(u(i,k,j)-u(i,k-1,j))/dsgmh(k)  
				               

	    enddo 
      ENDDO


		
		vortx = 0.0 
		
      j_start = jts
      j_end   = MIN(jte,jde-1)

	  

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end

    
    
      	vortx(i,k,j) =  (-1.0*vortx1(i,k,j) -  vortx21(i,k,j) - vortx22(i,k,j) )  * (msft(i,j)+msft(i+1,j))*0.5   

      ENDDO
      ENDDO
      ENDDO


	  
	  
	  
	  
 
 
   j_start = jts
   j_end   = min(jte,jde-1)

   IF ( config_flags%open_xs ) i_start = MAX(ids+1,its)
   IF ( config_flags%open_xe ) i_end   = MIN(ide-1,ite)

   IF ( config_flags%open_ys .or. specified ) i_start = MAX(ids+1,its)
   IF ( config_flags%open_ye .or. specified ) i_end   = MIN(ide-1,ite)
      IF ( config_flags%periodic_x ) i_start = its
      IF ( config_flags%periodic_x ) i_end = ite
	  
	  
	  
      DO  j = j_start, j_end
      DO k=2,ktf      
      DO  i = i_start, i_end

    
    
      	uds(i,k,j) = ( u(i,k,j) - u(i,k-1,j) )/dsgmh(k) 

      ENDDO
      ENDDO
      ENDDO	 
	  
	  
	 
	  DO  j = j_start, j_end
      DO  i = i_start, i_end

    
      	uds(i,kte,j) = uds(i,kte-1,j) 
        uds(i,kts,j) = uds(i,kts+1,j)
      ENDDO
      ENDDO
	  
	  
	  DO  j = j_start, j_end
      DO k=kts,kte  
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
  
      
      
      
      

           vortx(i,k,j) = vortx(i,k,j)  - 0.25* ( v(i-1,k,j)*v(i-1,k,j)   + v(i,k,j)*v(i,k,j)  + v(i-1,k,j+1)*v(i-1,k,j+1)  + v(i,k,j+1)*v(i,k,j+1)  )  & 
                       *(msft(i,j)-msft(i-1,j))*rdx    &
                             + u(i,k,j) * ( v(i,k,j)*(msft(i,j)-msft(i,j-1))*rdy    &
                                                   +  v(i-1,k,j)*(msft(i-1,j)-msft(i-1,j-1))*rdy    &
                                                   +  v(i,k,j+1)*(msft(i,j+1)-msft(i,j))*rdy    &
                                                   +  v(i-1,k,j+1)*(msft(i-1,j+1)-msft(i-1,j))*rdy  ) * 0.25 
       
        
        

      ENDDO
      ENDDO
      ENDDO	  


if (.true.) then  
      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end
  

   
    
       
        
    vortx(i,k,j) = vortx(i,k,j)  +  u(i,k,j)*reradius*(0.25*( v(i-1,k,j)   + v(i,k,j)  + v(i-1,k,j+1)  + v(i,k,j+1)  )*tan(xlat(i,j)*degrad) &
                                 - 0.25*(w(i-1,k+1,j)+w(i-1,k,j)+w(i,k+1,j)+w(i,k,j)) )

	
		

      ENDDO
      ENDDO
      ENDDO	   
 endif 
 


      END SUBROUTINE calcu_vortx	
















	  
	  
	
     
      SUBROUTINE calcu_gradx( var, gradx,                 &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var    
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: gradx

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx            
       
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: varz 

        REAL , DIMENSION( kms:kme ) :: sgmh,dsgm,dsgmh
		
		 REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: vardx 
		
          REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  ::  varbc
		
        LOGICAL :: specified


   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

 
      ktf=MIN(kte,kde-1)

 	  DO k=kts,kde-1
          sgmh(k) = 0.5*(zh(k)+zh(k+1))
          dsgm(k) = zh(k+1) - zh(k) 
      ENDDO

 	  DO k=2,kde-1
             dsgmh(k) = 0.5*(dsgm(k) + dsgm(k-1))   
			 
      ENDDO
	  
	  

         varbc  = var
	  
      CALL set_physical_bc3d( varbc, 'p', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              its, ite, jts, jte, kts, kte,       &
                              its, ite, jts, jte, kts, kte )



	 
       
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
 
      DO  i = i_start, i_end
       varz(i,kts,j) =  varz(i,kts+1,j)  
      varz(i,kte,j) =  varz(i,ktf,j)  
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

	
	   
	
   
      
      i_start = MAX(its-1,ids)
      i_end   = MIN(ite+1, ide)
	  
	  
      j_start = jts
      j_end   = MIN(jte,jde-1)	

	
	  
	  


      IF ( config_flags%open_xs ) i_start = MAX(ids+1,its)   
      IF ( config_flags%open_xe ) i_end   = MIN(ide-1,ite)
      IF ( specified ) i_start = MAX(ids+1,its-1)            
      IF ( specified ) i_end   = MIN(ide-1,ite+1)
	  
      IF ( config_flags%periodic_x ) i_start = its           
      IF ( config_flags%periodic_x ) i_end = ite



 
 



 
 
 
 
 
 

 

 
 


		
  loop_j :     DO  j = j_start, j_end
	
		 
	 
      DO  k = kts, ktf
      DO  i = i_start, i_end 
	  
	    mrdx= 0.5* ( msft(i-1,j) + msft(i,j) ) *rdx
		
	    c1 = (msft(i,j)*(ztop-sgmh(k))/(ztop-zs(i,j)) + msft(i-1,j)*(ztop-sgmh(k))/(ztop-zs(i-1,j)))*0.5
	    c2 = rdx*(zs(i,j)-zs(i-1,j))*( (varz(i-1,k,j) + varz(i,k,j)) + (varz(i-1,k+1,j) + varz(i,k+1,j)) ) *0.25	  
	


    
	
	   gradx(i,k,j) =  - mrdx*(varbc(i,k,j)-varbc(i-1,k,j))  + c1*c2
	  

	  ENDDO
      ENDDO	 

      ENDDO    loop_j
	  


      END SUBROUTINE calcu_gradx



 
      SUBROUTINE calcu_grady( var0, grady,                 &
 
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   
   

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

   
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdy
         
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: varz  
		   REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: var 

       REAL , DIMENSION( kms:kme ) :: sgmh,dsgm,dsgmh
    LOGICAL :: specified

   
   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

 
      ktf=MIN(kte,kde-1)


 

	  DO k=kts,kde-1
          sgmh(k) = 0.5*(zh(k)+zh(k+1))
          dsgm(k) = zh(k+1) - zh(k) 
      ENDDO

 	  DO k=2,kde-1
             dsgmh(k) = 0.5*(dsgm(k) + dsgm(k-1))   
			 
      ENDDO



	 
       
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
							  
							  

      IF ( config_flags%periodic_y ) j_start = jts            
      IF ( config_flags%periodic_y ) j_end = jte 


       varz = 0.0 

      DO  j = j_start, j_end
      DO k=2,ktf
      DO  i = i_start, i_end
    
      	varz(i,k,j) =  (var(i,k,j) - var(i,k-1,j) ) / dsgmh(k)
      ENDDO
      ENDDO
      ENDDO
	  
	  


      DO      j = j_start, j_end
  
      DO  i = i_start, i_end
		varz(i,kts,j) =  varz(i,kts+1,j)  
	   
		varz(i,kte,j) =  varz(i,kte-1,j)  
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
	  



 
     
    
    
	
      i_start = its
      i_end   =  MIN(ite,ide-1)
   
 
	   
	   
	   j_start = MAX(jts-1,jds)
       j_end   = MIN(jte+1,jde-1)
	  
      IF ( config_flags%periodic_y ) j_start = jts          
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
	   
	   
	
	
	  grady(i,k,j) =       - mrdy*(var(i,k,j)-var(i,k,j-1)) + c1*c2
		 
	  ENDDO
      ENDDO
      ENDDO    
	  



      END SUBROUTINE calcu_grady


  
     SUBROUTINE calcu_gradz( var, gradz,                 &
 
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var
	
    
    
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: gradz

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   
   
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
 
 
	  DO k=kts,kde-1
          sgmh(k) = 0.5*(zh(k)+zh(k+1))
          dsgm(k) = zh(k+1) - zh(k) 
      ENDDO

 	  DO k=2,kde-1
             dsgmh(k) = 0.5*(dsgm(k) + dsgm(k-1))   
			 
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



      DO      j = j_start, j_end
  
      DO  i = i_start, i_end
	  
       varz(i,kts,j) =  varz(i,kts+1,j)  
	   
      varz(i,kte,j) =  varz(i,kte-1,j)  
       	  
      ENDDO
      ENDDO
	  
	  
	  
	  


	  

	  gradz = 0.0
	  
      DO  j = j_start, j_end
      DO  k = kts, ktf
      DO  i = i_start, i_end 
			gradz(i,k,j) =  - ztop/(ztop-zs(i,j))*varz(i,k,j)
	  ENDDO
      ENDDO
      ENDDO    
	  
	   
	  

      END SUBROUTINE calcu_gradz
	  
	  




 SUBROUTINE average_vertical_w( varw, var,                 &
                     config_flags,                  &
                         fzm,fzp,                    &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var

        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: varw
		 
		 	  REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) ::  fzm,fzp
			  
			  


   
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx
         
     




  
    LOGICAL :: specified


	specified = .false.
	if(config_flags%specified .or. config_flags%nested) specified = .true.
   
     ktf=MIN(kte,kde-1)
 
 


      i_start = its
      i_end   = MIN(ide-1,ite)



      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
	  

      varw  = 0. 
	  
	  


    
    

    
    
   

    
    


    
    
	  
	  

      DO  j = j_start, j_end
      DO k=kts+1,ktf
      DO  i = i_start, i_end
       
      	varw(i,k,j) =  fzm(k)*var(i,k,j)  + fzp(k)*var(i,k-1,j)     
      ENDDO
      ENDDO
      ENDDO


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
   
   
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var,pfm_base,rhofm_base

        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: varw
		 
		 	  REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) ::  fzm,fzp,dnw
			  
		real :: pi,temp,theta	  


   
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx
         
     




  
    LOGICAL :: specified


	specified = .false.
	if(config_flags%specified .or. config_flags%nested) specified = .true.
   
     ktf=MIN(kte,kde-1)
 
 


      i_start = its
      i_end   = MIN(ide-1,ite)



      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
	  

      varw  = 0. 
	  
	  


    
    

    
    
   

    
    


    
    
	  
	  

      DO  j = j_start, j_end
      DO k=kts+1,ktf
      DO  i = i_start, i_end
       
      	varw(i,k,j) =  exp(fzm(k)*log(var(i,k,j))  + fzp(k)*log(var(i,k-1,j)))     
      ENDDO
      ENDDO
      ENDDO


      DO      j = j_start, j_end
  
      DO  i = i_start, i_end
     
	   
      varw(i,kte,j) =  2.0* var(i,kte-1,j)  -  varw(i,kte-1,j)   
      ENDDO
      ENDDO
	   
      DO      j = j_start, j_end
  
      DO  i = i_start, i_end
	    temp =  pfm_base(i,kts,j)/rhofm_base(i,kts,j)/r_d
		pi    = (var(i,kts,j)/p0)**(r_d/Cp)
		theta = temp/pi-6.5/1000.0*0.5*dnw(kts) 
		varw(i,kts,j) = p0*(temp/theta)**(Cp/r_d) + var(i,kts,j) - pfm_base(i,kts,j)
		
   
	   
     
      ENDDO
      ENDDO	  
	 
 

      END SUBROUTINE average_vertical_pw
	  
	  

 SUBROUTINE average_vertical_pws( varw, var, rhofm_base,                 &
                     config_flags,                  &
                         fzm,fzp,                    &
						          dnw,  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var, rhofm_base

        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: varw
		 
		 	  REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) ::  fzm,fzp,dnw
			  
		real :: pi,temp,theta	  


   
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx
         
     




  
    LOGICAL :: specified


	specified = .false.
	if(config_flags%specified .or. config_flags%nested) specified = .true.
   
     ktf=MIN(kte,kde-1)
 
 


      i_start = its
      i_end   = MIN(ide-1,ite)



      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
	  

      varw  = 0. 
	  
	  


    
    

    
    
   

    
    


    
    
	  
	  

      DO  j = j_start, j_end
      DO k=kts+1,ktf
      DO  i = i_start, i_end
       
     
        varw(i,k,j) =  exp(fzm(k)*log(var(i,k,j))  + fzp(k)*log(var(i,k-1,j)))   		
      ENDDO
      ENDDO
      ENDDO


      DO      j = j_start, j_end
  
      DO  i = i_start, i_end
    
	   
      varw(i,kte,j) =  2.0* var(i,kte-1,j)  -  varw(i,kte-1,j)   
      ENDDO
      ENDDO
	   
      DO      j = j_start, j_end
  
      DO  i = i_start, i_end
	    temp =  var(i,kts,j)/rhofm_base(i,kts,j)/r_d
		pi    = (var(i,kts,j)/p0)**(r_d/Cp)
		theta = temp/pi-6.5/1000.0*0.5*dnw(kts) 
		varw(i,kts,j) = p0*(temp/theta)**(Cp/r_d)   
		
   
	   
     
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
   
   
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var

        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: varw
		 
		 	  REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) ::  fzm,fzp
			  
			  


   
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx
         
     




  
    LOGICAL :: specified


	specified = .false.
	if(config_flags%specified .or. config_flags%nested) specified = .true.
   
     ktf=MIN(kte,kde-1)
 
 


      i_start = its
      i_end   = MIN(ide,ite)



      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
	  

      varw  = 0. 
	  
	  



 



    
    
   

      IF ( config_flags%periodic_x ) i_start = its           
      IF ( config_flags%periodic_x ) i_end = ite


 
 
   	  
	  

      DO  j = j_start, j_end
      DO k=kts+1,ktf
      DO  i = i_start, i_end
       
      	varw(i,k,j) =  fzm(k)*var(i,k,j)  + fzp(k)*var(i,k-1,j)     
      ENDDO
      ENDDO
      ENDDO


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

   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: fzm,  &     
                                                                  fzp,  &     
                                                                  zh

   REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy,  &
																  ztop
																  
   INTEGER ,                                     INTENT(IN   ) :: time_step

   
  
   
   
   
   
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
			
	ENDIf 
	
	
	
	cbx = 25.0 


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
						
						



     tendency = tendency + 	 engk_gradx  +   vortx  




		
	 
		

      IF ( (config_flags%open_xs) .and. its == ids ) THEN

        j_start = jts
        j_end   = MIN(jte,jde-1)

        DO j = j_start, j_end
        DO k = kts, ktf
		
         
		
	     ub = MIN(u(its,k,j)-cbx, 0.)   

		
    	   tendency(its,k,j) =     - rdx*ub*(u(its+1,k,j) - u(its,k,j))			  
 
		  
        ENDDO
        ENDDO

      ENDIF
	  
	  
	  

      IF ( (config_flags%open_xe) .and. ite == ide ) THEN

        j_start = jts
        j_end   = MIN(jte,jde-1)

        DO j = j_start, j_end
        DO k = kts, ktf
       	
		  ub = MAX(u(ite,k,j)+cbx, 0.)
	 
          tendency(ite,k,j) =                                  &
                      - rdx*ub*(u(ite,k,j) - u(ite-1,k,j))
		
		  
        ENDDO
        ENDDO

      ENDIF





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

         mrdy=msfux(i,jts)*rdy       
         ip = MIN( imax, i   )
         im = MAX( imin, i-1 )

         DO k=kts,ktf

    
		
		  vw = 0.5*(v(ip,k,jts)+v(im,k,jts)) 
          vb = MIN( vw, 0. )
	  
		
          tendency(i,k,jts) =    tendency(i,k,jts) - mrdy*vb*(u(i,k,jts+1)-u(i,k,jts))   
		
	
						   
         ENDDO
      ENDDO

   ENDIF

   IF( (config_flags%open_ye) .and. (jte == jde)) THEN

      DO i = i_start, i_end

         mrdy=msfux(i,jte-1)*rdy     
         ip = MIN( imax, i   )
         im = MAX( imin, i-1 )

         DO k=kts,ktf

	 
	 
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

   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: fzm,  &     
                                                                  fzp,  &     
                                                                  zh

   REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy,  &
																  ztop
																  
   INTEGER ,                                     INTENT(IN   ) :: time_step

   

   
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
	  
	ENDIf 
	
	
	cbx = 25.0  

   specified = .false.

   if(config_flags%specified .or. config_flags%nested) specified = .true.

       ktf=MIN(kte,kde-1) 
	  










call    calcu_engk(engk,u,v,w,config_flags,                       &
                                  ids, ide, jds, jde, kds, kde,  &
                                  ims, ime, jms, jme, kms, kme,  &
                                  its, ite, jts, jte, kts, kte  )



	
		
	
call    calcu_grady( engk, engk_grady,                 &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )




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






      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = jte

      IF ( config_flags%open_ys ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye ) j_end   = MIN(jde-1,jte)
   
   
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
	
 
	else
	
	
   

      IF ( (config_flags%open_ys) .and. jts == jds ) THEN

        i_start = its
        i_end   = MIN(ite,ide-1)

        DO i = i_start, i_end
        DO k = kts, ktf
		
      
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
		
      
          vb = MAX(v(i,k,jte)+cbx, 0.)
		  
          tendency(i,k,jte) =                   &
                    - rdy*vb*(v(i,k,jte) - v(i,k,jte-1))
			  
        ENDDO
        ENDDO

      ENDIF





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



   IF( (config_flags%open_xs) .and. (its == ids)) THEN

      DO j = j_start, j_end

         mrdx=msfvy(its,j)*rdx       
         jp = MIN( jmax, j   )
         jm = MAX( jmin, j-1 )

         DO k=kts,ktf

  
          uw = 0.5*(u(its,k,jp)+u(its,k,jm))
          ub = MIN( uw, 0. )
       

	
	
		       
              tendency(its,k,j) =  tendency(its,k,j) -mrdx*ub*(v(its+1,k,j)-v(its,k,j))   

		   
         ENDDO
      ENDDO

   ENDIF

   IF( (config_flags%open_xe) .and. (ite == ide) ) THEN
      DO j = j_start, j_end

         mrdx=msfvy(ite-1,j)*rdx     
         jp = MIN( jmax, j   )
         jm = MAX( jmin, j-1 )

         DO k=kts,ktf

   
          uw = 0.5*(u(ite,k,jp)+u(ite,k,jm))
          ub = MAX( uw, 0. )
	   
		  

     
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

   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: fzm,  &     
                                                                  fzp,  &     
                                                                  zh

   REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy,  &
																  ztop
																  
   INTEGER ,                                     INTENT(IN   ) :: time_step



   
   
   INTEGER :: i, j, k, itf, jtf, ktf
   INTEGER :: i_start, i_end, j_start, j_end, k_start, k_end 

   
   
      INTEGER :: jmin, jmax, jp, jm, imin, imax
	  REAL    :: mrdx, mrdy, ub, vb, uw, vw
	  

       REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) :: engk, engk_gradz, vortz   
   
   
   

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

	


     ktf=MIN(kte,kde-1)


   
      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)


	  
	
	
	 
	     k_start = kts
         k_end   = kte

	DO j = j_start, j_end
	
	
	DO k = k_start, k_end
     DO i = i_start, i_end
	     tendency(i,k,j)  = tendency(i,k,j)  + engk_gradz(i,k,j)  + vortz(i,k,j)  
 
 

     ENDDO
     ENDDO
	 
     ENDDO 
	 



   
 



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












   IMPLICIT NONE
   
   
   
   TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

   LOGICAL ,                 INTENT(IN   ) :: tenddec  

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

   
   
   INTEGER :: i, j, k, itf, jtf, ktf
   INTEGER :: i_start, i_end, j_start, j_end
   INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
   INTEGER :: jmin, jmax, jp, jm, imin, imax

   REAL    :: mrdx, mrdy, ub, vb, uw, vw, mu



   REAL,  DIMENSION( its-1:ite+2, kts:kte, jts-1:jte+2  ) :: fqx, fqy, fqz

   REAL,  DIMENSION( its-1:ite+2, kts:kte, jts-1:jte+2  ) :: fqxl, fqyl, fqzl

   INTEGER :: horz_order, vert_order
   
   LOGICAL :: degrade_xs, degrade_ys
   LOGICAL :: degrade_xe, degrade_ye

   INTEGER :: jp1, jp0, jtmp

   REAL,DIMENSION( its-1:ite+2, kts:kte, jts-1:jte+2  ) :: flux_out, ph_low
   REAL :: scale
   
   REAL, PARAMETER :: eps=1.e-20




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





    REAL     :: dx,dy,dz

    LOGICAL, PARAMETER :: pd_limit = .true.





    
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








   degrade_xs = .true.
   degrade_xe = .true.
   degrade_ys = .true.
   degrade_ye = .true.





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



      ktf=MIN(kte,kde-1)
      i_start = its-1
      i_end   = MIN(ite,ide-1)+1
      j_start = jts-1
      j_end   = MIN(jte,jde-1)+1



      IF(degrade_xs) i_start = its
      IF(degrade_xe) i_end   = MIN(ite,ide-1)
      IF(degrade_ys) j_start = MAX(jts,jds+1)
      IF(degrade_ye) j_end = MIN(jte,jde-2)



      DO j = j_start, j_end+1
        DO k=kts,ktf
        DO i = i_start, i_end
 
		   
		   fqy(i,k,j) =  0.5*v(i,k,j)*          &
                  ((ztop-zs(i,j))*field(i,k,j)/msft(i,j) + (ztop-zs(i,j-1))*field(i,k,j-1)/msft(i,j-1))				  
				  
		   
        ENDDO
        ENDDO
      ENDDO






	
    
      IF(degrade_xs) i_start = its
      IF(degrade_xe) i_end   = MIN(ite,ide-1)
      IF(degrade_ys) j_start = jts  
      IF(degrade_ye) j_end = MIN(jte,jde-1)  
	  


      DO j = j_start, j_end
        DO k=kts,ktf
        DO i = i_start, i_end+1
          
          
           vel =  u(i,k,j)
            cr = vel*dt/dx 
          
		      fqxl(i,k,j) =  (dx/dt)*flux_upwind((ztop-zs(i-1,j))*field(i-1,k,j)/msft(i-1,j), (ztop-zs(i,j))*field(i,k,j)/msft(i,j), cr)
		  
          
           
				  
            fqx( i,k,j ) = 0.5*u(i,k,j)*          &
                  ((ztop-zs(i,j))*field(i,k,j)/msft(i,j) + (ztop-zs(i-1,j))*field(i-1,k,j)/msft(i-1,j))				  
				  
				  

           fqx(i,k,j) = fqx(i,k,j) - fqxl(i,k,j)
        ENDDO
        ENDDO
      ENDDO



   ELSE

      WRITE ( wrf_err_message , * ) 'module_advect: advect_scalar_pd, h_order not known ',horz_order
      CALL wrf_error_fatal3("<stdin>",4326,&
TRIM( wrf_err_message ) )

   ENDIF horizontal_order_test


























      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)



   IF( (config_flags%open_xs) .and. (its == ids) ) THEN

       DO j = j_start, j_end
       DO k = kts, ktf
         ub = MIN( 0.5*(u(its,k,j)+u(its+1,k,j)), 0. )
         tendency(its,k,j) = tendency(its,k,j)                     &
             
             
             
             
             
                 - (ztop-zs(its,j))/msft(its,j)*rdx*(                                             &
                       ub*(   field_old(its+1,k,j)                 &
                            - field_old(its  ,k,j)   ) +           &
                       field_old(its,k,j)*(u(its+1,k,j)-u(its,k,j))  &       
                                                                ) * msft(its,j)*msft(its,j)/(ztop-zs(its,j))  																	
																
																
       ENDDO
       ENDDO

   ENDIF

   IF( (config_flags%open_xe) .and. (ite == ide) ) THEN

       DO j = j_start, j_end
       DO k = kts, ktf
         ub = MAX( 0.5*(u(ite-1,k,j)+u(ite,k,j)), 0. )
         tendency(i_end,k,j) = tendency(i_end,k,j)                   &
        
         
         
         
         
																	
                - (ztop-zs(i_end,j))/msft(i_end,j)*rdx*(                                               &
                       ub*(  field_old(i_end  ,k,j)                  &
                           - field_old(i_end-1,k,j) ) +              &
                       field_old(i_end,k,j)*(u(ite,k,j)-u(ite-1,k,j))  &   
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
                       field_old(i,k,jts)*(v(i,k,jts+1)-v(i,k,jts))  &     
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
                       field_old(i,k,j_end)*(v(i,k,jte)-v(i,k,jte-1))  &    
                                                                    )* msft(i,j_end)*msft(i,j_end)/(ztop-zs(i,j_end)) 
       ENDDO
       ENDDO

   ENDIF

 

   

  

  
 

  























      i_start = its-1
      i_end   = MIN(ite,ide-1)+1
      j_start = jts-1
      j_end   = MIN(jte,jde-1)+1



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

         
         
          
         
         
        
         
		   
	   fqz(i,k,j) = sigmdot(i,k,j) *(fzm(k)*field(i,k,j)+fzp(k)*field(i,k-1,j))*(ztop-zs(i,j))/msft(i,j)    
		   

        ENDDO
        ENDDO

      ENDDO

   ELSE

      WRITE (wrf_err_message,*) ' advect_scalar_pd, v_order not known ',vert_order
      CALL wrf_error_fatal3("<stdin>",4524,&
wrf_err_message )

   ENDIF vert_order_test
   
   
   
   
   
   
   
   
   
   
   
   
   
   

   IF (pd_limit) THEN



   i_start = its-1
   i_end   = MIN(ite,ide-1)+1
   j_start = jts-1
   j_end   = MIN(jte,jde-1)+1



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
   


   
   



   DO j=j_start, j_end
   DO k=kts, ktf
!DIR$ vector always
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


       IF( fqz (i,k+1,j) .lt. 0.) fqz(i,k+1,j) = scale*fqz(i,k+1,j)
       IF( fqz (i,k  ,j) .gt. 0.) fqz(i,k  ,j) = scale*fqz(i,k  ,j)

     END IF

   ENDDO
   ENDDO
   ENDDO

   END IF
   
   




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
                                   
									  )

  ENDDO
  ENDDO
  ENDDO
  

  
  IF(tenddec) THEN
  DO j = j_start, j_end
  DO k = kts, ktf
  DO i = i_start, i_end

     z_tendency (i,k,j) = 0. -rdzw(k)* msft(i,j)*msft(i,j)/(ztop-zs(i,j))*( fqz (i,k+1,j)-fqz (i,k,j)  &
                                    
									  )

  ENDDO
  ENDDO
  ENDDO
  END IF



  IF(degrade_xs) i_start = MAX(its,ids+1)
  IF(degrade_xe) i_end   = MIN(ite,ide-2)

  DO j = j_start, j_end
  DO k = kts, ktf
!DEC$ vector always  
  DO i = i_start, i_end

     
     tendency (i,k,j) = tendency(i,k,j)                           &
               - msft(i,j)*msft(i,j)/(ztop-zs(i,j))*( rdx*( fqx (i+1,k,j)-fqx (i,k,j)     &
                               
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
                              
								   ))

  ENDDO
  ENDDO
  ENDDO
  END IF



  i_start = its
  i_end   = MIN(ite,ide-1)
  IF(degrade_ys) j_start = MAX(jts,jds+1)
  IF(degrade_ye) j_end   = MIN(jte,jde-2)

  DO j = j_start, j_end
  DO k = kts, ktf
!DEC$ vector always
  DO i = i_start, i_end

     
     
     tendency (i,k,j) = tendency(i,k,j)                           &
               -  msft(i,j)*msft(i,j)/(ztop-zs(i,j))*( rdy*( fqy (i,k,j+1)-fqy (i,k,j)     &
             
								   ))

  ENDDO
  ENDDO
  ENDDO
  

  

  IF(tenddec) THEN
  DO j = j_start, j_end
  DO k = kts, ktf
  DO i = i_start, i_end

     h_tendency (i,k,j) = h_tendency (i,k,j)                      &
               -  msft(i,j)*msft(i,j)/(ztop-zs(i,j))*( rdy*( fqy (i,k,j+1)-fqy (i,k,j)     &
                             
								  ) )

  ENDDO
  ENDDO
  ENDDO
  END IF

END SUBROUTINE advect_scalar_pd_fm_old



SUBROUTINE advect_scalar_fm2   ( field, field_old, tendency,    &
                             
                                u, v, sigmdot,                   &
                           
                            
                                time_step, config_flags,       &
								msft,   &
                                fzm, fzp, ztop,zs,                       &
                                rdx, rdy, rdzw, dt,            &
                                ids, ide, jds, jde, kds, kde,  &
                                ims, ime, jms, jme, kms, kme,  &
                                its, ite, jts, jte, kts, kte  )












   IMPLICIT NONE
   
   
   
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

   
   
   INTEGER :: i, j, k, itf, jtf, ktf
   INTEGER :: i_start, i_end, j_start, j_end
   INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
   INTEGER :: jmin, jmax, jp, jm, imin, imax

   REAL    :: mrdx, mrdy, ub, vb, uw, vw, mu



   REAL,  DIMENSION( its-1:ite+2, kts:kte, jts-1:jte+2  ) :: fqx, fqy, fqz

   REAL,  DIMENSION( its-1:ite+2, kts:kte, jts-1:jte+2  ) :: fqxl, fqyl, fqzl

   INTEGER :: horz_order, vert_order
   
   LOGICAL :: degrade_xs, degrade_ys
   LOGICAL :: degrade_xe, degrade_ye

   INTEGER :: jp1, jp0, jtmp

   REAL,DIMENSION( its-1:ite+2, kts:kte, jts-1:jte+2  ) :: flux_out, ph_low
   REAL :: scale
   
   REAL, PARAMETER :: eps=1.e-20




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





    REAL     :: dx,dy,dz

    LOGICAL, PARAMETER :: pd_limit = .false.





    
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








   degrade_xs = .true.
   degrade_xe = .true.
   degrade_ys = .true.
   degrade_ye = .true.





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



      ktf=MIN(kte,kde-1)
      i_start = its-1
      i_end   = MIN(ite,ide-1)+1
      j_start = jts-1
      j_end   = MIN(jte,jde-1)+1



      IF(degrade_xs) i_start = its
      IF(degrade_xe) i_end   = MIN(ite,ide-1)
      IF(degrade_ys) j_start = MAX(jts,jds+1)
      IF(degrade_ye) j_end = MIN(jte,jde-2)



      DO j = j_start, j_end+1
        DO k=kts,ktf
        DO i = i_start, i_end
 
		   
		   fqy(i,k,j) =  0.5*v(i,k,j)*          &
                  ((ztop-zs(i,j))*field(i,k,j)/msft(i,j) + (ztop-zs(i,j-1))*field(i,k,j-1)/msft(i,j-1))				  
				  
		   
        ENDDO
        ENDDO
      ENDDO



      DO j = j_start, j_end
        DO k=kts,ktf
        DO i = i_start, i_end+1
 
				  
            fqx( i,k,j ) = 0.5*u(i,k,j)*          &
                  ((ztop-zs(i,j))*field(i,k,j)/msft(i,j) + (ztop-zs(i-1,j))*field(i-1,k,j)/msft(i-1,j))				  
				  
				  
 
        ENDDO
        ENDDO
      ENDDO



   ELSE

      WRITE ( wrf_err_message , * ) 'module_advect: advect_scalar_pd, h_order not known ',horz_order
      CALL wrf_error_fatal3("<stdin>",4991,&
TRIM( wrf_err_message ) )

   ENDIF horizontal_order_test


























      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)



   IF( (config_flags%open_xs) .and. (its == ids) ) THEN

       DO j = j_start, j_end
       DO k = kts, ktf
         ub = MIN( 0.5*(u(its,k,j)+u(its+1,k,j)), 0. )
         tendency(its,k,j) = tendency(its,k,j)                     &
             
             
             
             
             
                 - (ztop-zs(its,j))/msft(its,j)*rdx*(                                             &
                       ub*(   field_old(its+1,k,j)                 &
                            - field_old(its  ,k,j)   ) +           &
                       field_old(its,k,j)*(u(its+1,k,j)-u(its,k,j))  &       
                                                                ) * msft(its,j)*msft(its,j)/(ztop-zs(its,j))  																	
																
																
       ENDDO
       ENDDO

   ENDIF

   IF( (config_flags%open_xe) .and. (ite == ide) ) THEN

       DO j = j_start, j_end
       DO k = kts, ktf
         ub = MAX( 0.5*(u(ite-1,k,j)+u(ite,k,j)), 0. )
         tendency(i_end,k,j) = tendency(i_end,k,j)                   &
        
         
         
         
         
																	
                - (ztop-zs(i_end,j))/msft(i_end,j)*rdx*(                                               &
                       ub*(  field_old(i_end  ,k,j)                  &
                           - field_old(i_end-1,k,j) ) +              &
                       field_old(i_end,k,j)*(u(ite,k,j)-u(ite-1,k,j))  &   
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
                       field_old(i,k,jts)*(v(i,k,jts+1)-v(i,k,jts))  &     
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
                       field_old(i,k,j_end)*(v(i,k,jte)-v(i,k,jte-1))  &    
                                                                    )* msft(i,j_end)*msft(i,j_end)/(ztop-zs(i,j_end)) 
       ENDDO
       ENDDO

   ENDIF

 

   

  

  
 

  























      i_start = its-1
      i_end   = MIN(ite,ide-1)+1
      j_start = jts-1
      j_end   = MIN(jte,jde-1)+1



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

         
         
          
         
         
        
         
		   
	   fqz(i,k,j) = sigmdot(i,k,j) *(fzm(k)*field(i,k,j)+fzp(k)*field(i,k-1,j))*(ztop-zs(i,j))/msft(i,j)    
		   

        ENDDO
        ENDDO

      ENDDO

   ELSE

      WRITE (wrf_err_message,*) ' advect_scalar_pd, v_order not known ',vert_order
      CALL wrf_error_fatal3("<stdin>",5189,&
wrf_err_message )

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
                                   
									  )

  ENDDO
  ENDDO
  ENDDO
  

  
 



  IF(degrade_xs) i_start = MAX(its,ids+1)
  IF(degrade_xe) i_end   = MIN(ite,ide-2)

  DO j = j_start, j_end
  DO k = kts, ktf
!DEC$ vector always  
  DO i = i_start, i_end

     
     tendency (i,k,j) = tendency(i,k,j)                           &
               - msft(i,j)*msft(i,j)/(ztop-zs(i,j))*( rdx*( fqx (i+1,k,j)-fqx (i,k,j)     &
                               
								   ) )

  ENDDO
  ENDDO
  ENDDO
  
 



  i_start = its
  i_end   = MIN(ite,ide-1)
  IF(degrade_ys) j_start = MAX(jts,jds+1)
  IF(degrade_ye) j_end   = MIN(jte,jde-2)

  DO j = j_start, j_end
  DO k = kts, ktf
!DEC$ vector always
  DO i = i_start, i_end

     
     
     tendency (i,k,j) = tendency(i,k,j)                           &
               -  msft(i,j)*msft(i,j)/(ztop-zs(i,j))*( rdy*( fqy (i,k,j+1)-fqy (i,k,j)     &
             
								   ))

  ENDDO
  ENDDO
  ENDDO
  

  

 

END SUBROUTINE advect_scalar_fm2




SUBROUTINE calcu_divg2( var, divg,                 &
                        u, v, sgmdot,                       &
                     config_flags,                  &
                        msft,zh,zs,fzm,fzp,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var,u,     &
                                                                      v, &
                                                                      sgmdot
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: divg

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh,fzm,fzp


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx
         
        
       
	   
	   
	   
	   
	   
		 
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


 



	  

   
      CALL set_physical_bc3d( var, 'p', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              its, ite, jts, jte, kts, kte,       &
                              its, ite, jts, jte, kts, kte )
	  
	  

     divg = 0.0 
	 
	 

	 

 

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)


      IF ( .NOT. config_flags%periodic_x ) THEN
        IF ( config_flags%open_xs .or. specified ) i_start = MAX(ids+1,its)
        IF ( config_flags%open_xe .or. specified ) i_end   = MIN(ide-2,ite)
      ENDIF

      
	  
	  
	IF ( config_flags%open_xs  .and.  its == ids ) THEN 
		
	
		
      DO  j = j_start, j_end
      DO k=kts,ktf
   
 					 
	
	         ub = MIN( 0.5*(u(its,k,j)+u(its+1,k,j)), 0. )
		 
         divg(its,k,j) =                      &
               - (ztop-zs(its,j))/msft(its,j)*rdx*(                                             &
                       ub*(   var(its+1,k,j)                 &
                            - var(its  ,k,j)   ) +           &
                       var(its,k,j)*(u(its+1,k,j)-u(its,k,j))  &       
                                                                )				 
					 
      ENDDO
      ENDDO
 
	  
	ENDIF	  
	

	  
	    IF ( config_flags%open_xe  .and.  ite == ide ) THEN 
		
	 
		
		 i_end   = MIN(ite,ide-1)
		 
		
      DO  j = j_start, j_end
      DO k=kts,ktf
 
  
     
	     ub = MAX( 0.5*(u(ite-1,k,j)+u(ite,k,j)), 0. )	 
		 
         divg(i_end,k,j) =                  &
               - (ztop-zs(i_end,j))/msft(i_end,j)*rdx*(                                               &
                       ub*(  var(i_end  ,k,j)                  &
                           - var(i_end-1,k,j) ) +              &
                       var(i_end,k,j)*(u(ite,k,j)-u(ite-1,k,j))  &   
                                                                    )				 
					 
					 
      ENDDO
      ENDDO
  

	ENDIF	
	  
	  
	  
	  



      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
      
      IF ( config_flags%open_ys .or. specified .or. config_flags%polar ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified .or. config_flags%polar ) j_end   = MIN(jde-2,jte) 	  

       
	  
	    IF ( config_flags%open_ys  .and.  jts == jds ) THEN 
		
	
		
      DO k=kts,ktf
      DO  i = i_start, i_end		

	
	      vb = MIN( 0.5*(v(i,k,jts)+v(i,k,jts+1)), 0. )
		  
        divg(i,k,jts) =    divg(i,k,jts)                &
               - (ztop-zs(i,jts))/msft(i,jts)*rdy*(                                             &
                       vb*(  var(i,k,jts+1)                  &
                           - var(i,k,jts  ) ) +              &
                       var(i,k,jts)*(v(i,k,jts+1)-v(i,k,jts))  &     
                                                               )
					 
      ENDDO
      ENDDO
	    ENDIF	
		
 
	  
	    IF ( config_flags%open_ye  .and.  jte == jde ) THEN 
		
	 
		
		j_end   = MIN(jte,jde-1)
		
      DO k=kts,ktf
      DO  i = i_start, i_end		
		
	 
	
		  vb = MAX( 0.5*(v(i,k,jte-1)+v(i,k,jte)), 0. ) 
		 
         divg(i,k,j_end) =    divg(i,k,j_end)               &
               -  (ztop-zs(i,j_end))/msft(i,j_end)*rdy*(                                               &
                       vb*(   var(i,k,j_end  )                 &
                            - var(i,k,j_end-1) ) +             &
                       var(i,k,j_end)*(v(i,k,jte)-v(i,k,jte-1))  &    
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
	  
	



      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
      
      IF ( config_flags%open_ys .or. specified .or. config_flags%polar ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified .or. config_flags%polar ) j_end   = MIN(jde-2,jte) 	  

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end

   
   
					 
 		 
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
   
   
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var,u,     &
                                                                      v, &
                                                                      sgmdot
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: divg

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh,fzm,fzp


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx
         
        
       
	   
	   
	   
	   
	   
		 
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
	  
	  
	IF ( config_flags%open_xs  .and.  its == ids ) THEN 
		
	
		
      DO  j = j_start, j_end
      DO k=kts,ktf
  
  
 
     
     


					 
	
	         ub = MIN( 0.5*(u(its,k,j)+u(its+1,k,j)), 0. )
		 
         divg(its,k,j) =   divg(its,k,j)                   &
               - (ztop-zs(its,j))/msft(its,j)*rdx*(                                             &
                       ub*(   var(its+1,k,j)                 &
                            - var(its  ,k,j)   ) +           &
                       var(its,k,j)*(u(its+1,k,j)-u(its,k,j))  &       
                                                                )				 
					 
      ENDDO
      ENDDO
 
	  
	ENDIF	  
	

	  
	    IF ( config_flags%open_xe  .and.  ite == ide ) THEN 
		
	 
		
		 i_end   = MIN(ite,ide-1)
		 
		
      DO  j = j_start, j_end
      DO k=kts,ktf
 
  
      
      
	
	     ub = MAX( 0.5*(u(ite-1,k,j)+u(ite,k,j)), 0. )	 
		 
         divg(i_end,k,j) =   divg(i_end,k,j)                &
               - (ztop-zs(i_end,j))/msft(i_end,j)*rdx*(                                               &
                       ub*(  var(i_end  ,k,j)                  &
                           - var(i_end-1,k,j) ) +              &
                       var(i_end,k,j)*(u(ite,k,j)-u(ite-1,k,j))  &   
                                                                    )				 
					 
					 
      ENDDO
      ENDDO
  

	ENDIF	
	  
	  
	  
	  



      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
      
      IF ( config_flags%open_ys .or. specified .or. config_flags%polar ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified .or. config_flags%polar ) j_end   = MIN(jde-2,jte) 	  

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end

   
   
					 
 		 
    	divg(i,k,j) = divg(i,k,j) -   ( 0.5*((ztop-zs(i,j))*var(i,k,j)/msft(i,j)  + (ztop-zs(i,j+1))*var(i,k,j+1)/msft(i,j+1) ) *v(i,k,j+1)    &
                 -   0.5*((ztop-zs(i,j-1))*var(i,k,j-1)/msft(i,j-1)  + (ztop-zs(i,j))*var(i,k,j)/msft(i,j) ) *v(i,k,j) )*rdy					 

   

       ENDDO
      ENDDO
      ENDDO
	  
	  

	  
	    IF ( config_flags%open_ys  .and.  jts == jds ) THEN 
		
	
		
      DO k=kts,ktf
      DO  i = i_start, i_end		
		
    
    
	
	
	
	      vb = MIN( 0.5*(v(i,k,jts)+v(i,k,jts+1)), 0. )
		  
        divg(i,k,jts) =    divg(i,k,jts)                &
               - (ztop-zs(i,jts))/msft(i,jts)*rdy*(                                             &
                       vb*(  var(i,k,jts+1)                  &
                           - var(i,k,jts  ) ) +              &
                       var(i,k,jts)*(v(i,k,jts+1)-v(i,k,jts))  &     
                                                               )
					 
      ENDDO
      ENDDO
	    ENDIF	
		
 
	  
	    IF ( config_flags%open_ye  .and.  jte == jde ) THEN 
		
	 
		
		j_end   = MIN(jte,jde-1)
		
      DO k=kts,ktf
      DO  i = i_start, i_end		
		
	  
       
	
		  vb = MAX( 0.5*(v(i,k,jte-1)+v(i,k,jte)), 0. ) 
		 
         divg(i,k,j_end) =    divg(i,k,j_end)               &
               -  (ztop-zs(i,j_end))/msft(i,j_end)*rdy*(                                               &
                       vb*(   var(i,k,j_end  )                 &
                            - var(i,k,j_end-1) ) +             &
                       var(i,k,j_end)*(v(i,k,jte)-v(i,k,jte-1))  &    
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
	  
	  
	  


SUBROUTINE calcu_divg_bak2024( var, divg,                 &
                        u, v, sgmdot,                       &
                     config_flags,                  &
                        msft,zh,zs,fzm,fzp,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var,u,     &
                                                                      v, &
                                                                      sgmdot
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: divg

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh,fzm,fzp


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx
         
        
       
	   
	   
	   
	   
	   
		 
		    REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) ::  varw
		 
		 
		 
		 


      REAL , DIMENSION( kms:kme ) :: sgmh,dsgm
       REAL , DIMENSION( ims:ime , jms:jme )  ::  rmm    
	   
	   REAL :: ub,vb 
				 
				 
		
				 
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


 



	  

   
      CALL set_physical_bc3d( var, 'p', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              its, ite, jts, jte, kts, kte,       &
                              its, ite, jts, jte, kts, kte )
	  
	  

     divg = 0.0 
	 
	 

	 

 

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)


      IF ( .NOT. config_flags%periodic_x ) THEN
        IF ( config_flags%open_xs .or. specified ) i_start = MAX(ids+1,its)
        IF ( config_flags%open_xe .or. specified ) i_end   = MIN(ide-2,ite)
      ENDIF

      
	  
	  
	IF ( config_flags%open_xs  .and.  its == ids ) THEN 
		
	
		
      DO  j = j_start, j_end
      DO k=kts,ktf
   
 					 
	
	         ub = MIN( 0.5*(u(its,k,j)+u(its+1,k,j)), 0. )
		 
         divg(its,k,j) =                      &
               - (ztop-zs(its,j))/msft(its,j)*rdx*(                                             &
                       ub*(   var(its+1,k,j)                 &
                            - var(its  ,k,j)   ) +           &
                       var(its,k,j)*(u(its+1,k,j)-u(its,k,j))  &       
                                                                )				 
					 
      ENDDO
      ENDDO
 
	  
	ENDIF	  


	  
	    IF ( config_flags%open_xe  .and.  ite == ide ) THEN 
		
	 
		
		 i_end   = MIN(ite,ide-1)
		 
		
      DO  j = j_start, j_end
      DO k=kts,ktf
 
  
     
	     ub = MAX( 0.5*(u(ite-1,k,j)+u(ite,k,j)), 0. )	 
		 
         divg(i_end,k,j) =                  &
               - (ztop-zs(i_end,j))/msft(i_end,j)*rdx*(                                               &
                       ub*(  var(i_end  ,k,j)                  &
                           - var(i_end-1,k,j) ) +              &
                       var(i_end,k,j)*(u(ite,k,j)-u(ite-1,k,j))  &   
                                                                    )				 
					 
					 
      ENDDO
      ENDDO
  

	ENDIF	
	  
	  
	
	  



      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
      
      IF ( config_flags%open_ys .or. specified .or. config_flags%polar ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified .or. config_flags%polar ) j_end   = MIN(jde-2,jte) 	  

       
	  
	    IF ( config_flags%open_ys  .and.  jts == jds ) THEN 
		
	
		
      DO k=kts,ktf
      DO  i = i_start, i_end		

	
	      vb = MIN( 0.5*(v(i,k,jts)+v(i,k,jts+1)), 0. )
		  
        divg(i,k,jts) =    divg(i,k,jts)                &
               - (ztop-zs(i,jts))/msft(i,jts)*rdy*(                                             &
                       vb*(  var(i,k,jts+1)                  &
                           - var(i,k,jts  ) ) +              &
                       var(i,k,jts)*(v(i,k,jts+1)-v(i,k,jts))  &     
                                                               )
					 
      ENDDO
      ENDDO
	    ENDIF	
		

	  
	    IF ( config_flags%open_ye  .and.  jte == jde ) THEN 
		
	 
		
		j_end   = MIN(jte,jde-1)
		
      DO k=kts,ktf
      DO  i = i_start, i_end		
		
	 
	
		  vb = MAX( 0.5*(v(i,k,jte-1)+v(i,k,jte)), 0. ) 
		 
         divg(i,k,j_end) =    divg(i,k,j_end)               &
               -  (ztop-zs(i,j_end))/msft(i,j_end)*rdy*(                                               &
                       vb*(   var(i,k,j_end  )                 &
                            - var(i,k,j_end-1) ) +             &
                       var(i,k,j_end)*(v(i,k,jte)-v(i,k,jte-1))  &    
                                                                    )
					 
      ENDDO
      ENDDO
	    ENDIF
		
		









	 
      field  = var  	 
	 
	 












	  


      i_start = its-1
      i_end   = MIN(ite,ide-1)+1
      j_start = jts-1
      j_end   = MIN(jte,jde-1)+1



    IF(degrade_xs) i_start = MAX(its-1,ids)
    IF(degrade_xe) i_end   = MIN(ite+1,ide-1)
    IF(degrade_ys) j_start = MAX(jts-1,jds)
    IF(degrade_ye) j_end   = MIN(jte+1,jde-1)

    

      DO j = j_start, j_end

         DO i = i_start, i_end
           fqz(i,1,j)  = 0.
        
           fqz(i,kde,j)  = 0.
          
         ENDDO

         DO k=kts+1,ktf
         DO i = i_start, i_end


		   
	   fqz(i,k,j) = sgmdot(i,k,j) *(fzm(k)*field(i,k,j)+fzp(k)*field(i,k-1,j))*(ztop-zs(i,j))/msft(i,j)    
		   

        ENDDO
        ENDDO

      ENDDO

    
   




  
	  
	  	call average_vertical_w( varw, var,                 &
                     config_flags,                  &
                         fzm,fzp,                    &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )	



  
  
  
        i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
	  
      DO  j = j_start, j_end
      DO  i = i_start, i_end
      
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



      ktf=MIN(kte,kde-1)
      i_start = its-1
      i_end   = MIN(ite,ide-1)+1
      j_start = jts-1
      j_end   = MIN(jte,jde-1)+1



      IF(degrade_xs) i_start = its
      IF(degrade_xe) i_end   = MIN(ite,ide-1)
      IF(degrade_ys) j_start = MAX(jts,jds+1)
      IF(degrade_ye) j_end = MIN(jte,jde-2)


      DO j = j_start, j_end+1
        DO k=kts,ktf
        DO i = i_start, i_end
 
		   
	
    

				  
	 
		   
		    fqy(i,k,j) =    0.5*((ztop-zs(i,j-1))*var(i,k,j-1)/msft(i,j-1)  + (ztop-zs(i,j))*var(i,k,j)/msft(i,j) ) *v(i,k,j) 	
        ENDDO
        ENDDO
      ENDDO
	  
	  
	  
	  j_start = jts 
      j_end   = MIN(jte,jde-1) 
	  
	  
      DO j = j_start, j_end
        DO k=kts,ktf
        DO i = i_start, i_end+1
 
				  
      
      

		 	 fqx( i,k,j ) =  	  0.5*((ztop-zs(i-1,j))*var(i-1,k,j)/msft(i-1,j)  + (ztop-zs(i,j))*var(i,k,j)/msft(i,j) ) *u(i,k,j)   
 
 
 
 
        ENDDO
        ENDDO
      ENDDO
	  
	  



	 
	 
	 

 

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
	  





      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
      
      IF ( config_flags%open_ys .or. specified .or. config_flags%polar ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified .or. config_flags%polar ) j_end   = MIN(jde-2,jte) 	  

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end

    
    
					 
 		 
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












   IMPLICIT NONE
   
   
   
   TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

   LOGICAL ,                 INTENT(IN   ) :: tenddec  

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

   
   
   INTEGER :: i, j, k, itf, jtf, ktf
   INTEGER :: i_start, i_end, j_start, j_end
   INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
   INTEGER :: jmin, jmax, jp, jm, imin, imax

   REAL    :: mrdx, mrdy, ub, vb, uw, vw, mu



   REAL,  DIMENSION( its-1:ite+2, kts:kte, jts-1:jte+2  ) :: fqx, fqy, fqz

   REAL,  DIMENSION( its-1:ite+2, kts:kte, jts-1:jte+2  ) :: fqxl, fqyl, fqzl

   INTEGER :: horz_order, vert_order
   
   LOGICAL :: degrade_xs, degrade_ys
   LOGICAL :: degrade_xe, degrade_ye

   INTEGER :: jp1, jp0, jtmp

   REAL,DIMENSION( its-1:ite+2, kts:kte, jts-1:jte+2  ) :: flux_out, ph_low
   REAL :: scale
   
   REAL, PARAMETER :: eps=1.e-20




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





    REAL     :: dx,dy,dz

    LOGICAL, PARAMETER :: pd_limit = .true.





    
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








   degrade_xs = .true.
   degrade_xe = .true.
   degrade_ys = .true.
   degrade_ye = .true.





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



      ktf=MIN(kte,kde-1)
      i_start = its-1
      i_end   = MIN(ite,ide-1)+1
      j_start = jts-1
      j_end   = MIN(jte,jde-1)+1



      IF(degrade_xs) i_start = its
      IF(degrade_xe) i_end   = MIN(ite,ide-1)
      IF(degrade_ys) j_start = MAX(jts,jds+1)
      IF(degrade_ye) j_end = MIN(jte,jde-2)



      DO j = j_start, j_end+1
        DO k=kts,ktf
        DO i = i_start, i_end
 
 
           dy = 2./(msftx(i,j)+msftx(i,j-1))/rdy  
        
          vel =  v(i,k,j)
          cr = vel*dt/dy       
        

 
		            fqyl(i,k,j) =  (dy/dt)*flux_upwind((ztop-zs(i,j-1))*field(i,k,j-1)/msft(i,j-1), (ztop-zs(i,j))*field(i,k,j)/msft(i,j), cr)
		   
		   fqy(i,k,j) =  0.5*v(i,k,j)*          &
                  ((ztop-zs(i,j))*field(i,k,j)/msft(i,j) + (ztop-zs(i,j-1))*field(i,k,j-1)/msft(i,j-1))				  
				  
		       fqy(i,k,j) = fqy(i,k,j) - fqyl(i,k,j)
			   
			   
			   
        ENDDO
        ENDDO
      ENDDO






	
    
      IF(degrade_xs) i_start = its
      IF(degrade_xe) i_end   = MIN(ite,ide-1)
      IF(degrade_ys) j_start = jts  
      IF(degrade_ye) j_end = MIN(jte,jde-1)  
	  


      DO j = j_start, j_end
        DO k=kts,ktf
        DO i = i_start, i_end+1
            dx = 2./(msfty(i,j)+msfty(i-1,j))/rdx 
          
           vel =  u(i,k,j)
         cr = vel*dt/dx 
          fqxl(i,k,j) =  (dx/dt)*flux_upwind(field_old(i-1,k,j), field_old(i,k,j  ), cr)
          
           
				  
            fqx( i,k,j ) = 0.5*u(i,k,j)*          &
                  ((ztop-zs(i,j))*field(i,k,j)/msft(i,j) + (ztop-zs(i-1,j))*field(i-1,k,j)/msft(i-1,j))				  
				  
				  

          fqx(i,k,j) = fqx(i,k,j) - fqxl(i,k,j)
        ENDDO
        ENDDO
      ENDDO



   ELSE

      WRITE ( wrf_err_message , * ) 'module_advect: advect_scalar_pd, h_order not known ',horz_order
      CALL wrf_error_fatal3("<stdin>",6805,&
TRIM( wrf_err_message ) )

   ENDIF horizontal_order_test


























      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)



   IF( (config_flags%open_xs) .and. (its == ids) ) THEN

       DO j = j_start, j_end
       DO k = kts, ktf
         ub = MIN( 0.5*(u(its,k,j)+u(its+1,k,j)), 0. )
         tendency(its,k,j) = tendency(its,k,j)                     &
             
             
             
             
             
                 - (ztop-zs(its,j))/msft(its,j)*rdx*(                                             &
                       ub*(   field_old(its+1,k,j)                 &
                            - field_old(its  ,k,j)   ) +           &
                       field_old(its,k,j)*(u(its+1,k,j)-u(its,k,j))  &       
                                                                ) * msft(its,j)*msft(its,j)/(ztop-zs(its,j))  																	
																
																
       ENDDO
       ENDDO

   ENDIF

   IF( (config_flags%open_xe) .and. (ite == ide) ) THEN

       DO j = j_start, j_end
       DO k = kts, ktf
         ub = MAX( 0.5*(u(ite-1,k,j)+u(ite,k,j)), 0. )
         tendency(i_end,k,j) = tendency(i_end,k,j)                   &
        
         
         
         
         
																	
                - (ztop-zs(i_end,j))/msft(i_end,j)*rdx*(                                               &
                       ub*(  field_old(i_end  ,k,j)                  &
                           - field_old(i_end-1,k,j) ) +              &
                       field_old(i_end,k,j)*(u(ite,k,j)-u(ite-1,k,j))  &   
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
                       field_old(i,k,jts)*(v(i,k,jts+1)-v(i,k,jts))  &     
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
                       field_old(i,k,j_end)*(v(i,k,jte)-v(i,k,jte-1))  &    
                                                                    )* msft(i,j_end)*msft(i,j_end)/(ztop-zs(i,j_end)) 
       ENDDO
       ENDDO

   ENDIF

 

   

  

  
 

  























      i_start = its-1
      i_end   = MIN(ite,ide-1)+1
      j_start = jts-1
      j_end   = MIN(jte,jde-1)+1



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
         
            vel = sigmdot(i,k,j)  
           cr = vel*dt/dz  
         
        
         
		  fqzl(i,k,j) = (dz/dt)*flux_upwind(field_old(i,k-1,j), field_old(i,k,j  ), cr)  *(ztop-zs(i,j))/msft(i,j)   

		
	   fqz(i,k,j) = sigmdot(i,k,j) *(fzm(k)*field(i,k,j)+fzp(k)*field(i,k-1,j))*(ztop-zs(i,j))/msft(i,j)    
		  fqz(i,k,j) = fqz(i,k,j) - fqzl(i,k,j)  

        ENDDO
        ENDDO

      ENDDO

   ELSE

      WRITE (wrf_err_message,*) ' advect_scalar_pd, v_order not known ',vert_order
      CALL wrf_error_fatal3("<stdin>",7005,&
wrf_err_message )

   ENDIF vert_order_test
   
   
   
   
   
   
   
   
   
   
   
   
   
   

   IF (pd_limit) THEN



   i_start = its-1
   i_end   = MIN(ite,ide-1)+1
   j_start = jts-1
   j_end   = MIN(jte,jde-1)+1



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
   


   
   



   DO j=j_start, j_end
   DO k=kts, ktf
!DIR$ vector always
   DO i=i_start, i_end

  
  
   
   
   
                 

   
   
   
  
   


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


       IF( fqz (i,k+1,j) .lt. 0.) fqz(i,k+1,j) = scale*fqz(i,k+1,j)
       IF( fqz (i,k  ,j) .gt. 0.) fqz(i,k  ,j) = scale*fqz(i,k  ,j)

     END IF

   ENDDO
   ENDDO
   ENDDO

   END IF
   
   




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



  IF(degrade_xs) i_start = MAX(its,ids+1)
  IF(degrade_xe) i_end   = MIN(ite,ide-2)

  DO j = j_start, j_end
  DO k = kts, ktf
!DEC$ vector always  
  DO i = i_start, i_end

     
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



  i_start = its
  i_end   = MIN(ite,ide-1)
  IF(degrade_ys) j_start = MAX(jts,jds+1)
  IF(degrade_ye) j_end   = MIN(jte,jde-2)

  DO j = j_start, j_end
  DO k = kts, ktf
!DEC$ vector always
  DO i = i_start, i_end

     
     
     tendency (i,k,j) = tendency(i,k,j)                           &
               -  msft(i,j)*msft(i,j)/(ztop-zs(i,j))*( rdy*( fqy (i,k,j+1)-fqy (i,k,j)     &
                                 +fqyl(i,k,j+1)-fqyl(i,k,j)       &
								   ))

  ENDDO
  ENDDO
  ENDDO
  

  

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












   IMPLICIT NONE
   
   
   
   TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

   LOGICAL ,                 INTENT(IN   ) :: tenddec  

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

   
   
   INTEGER :: i, j, k, itf, jtf, ktf
   INTEGER :: i_start, i_end, j_start, j_end
   INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
   INTEGER :: jmin, jmax, jp, jm, imin, imax

   REAL    :: mrdx, mrdy, ub, vb, uw, vw, mu



 
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





    REAL     :: dx,dy,dz
	
	

    LOGICAL, PARAMETER :: mono_limit = .true.




    


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
















   degrade_xs = .true.
   degrade_xe = .true.
   degrade_ys = .true.
   degrade_ye = .true.





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



      ktf=MIN(kte,kde-1)
      i_start = its-1
      i_end   = MIN(ite,ide-1)+1
      j_start = jts-1
      j_end   = MIN(jte,jde-1)+1



      IF(degrade_xs) i_start = its
      IF(degrade_xe) i_end   = MIN(ite,ide-1)
      IF(degrade_ys) j_start = MAX(jts,jds+1)
      IF(degrade_ye) j_end = MIN(jte,jde-2)



      DO j = j_start, j_end+1
        DO k=kts,ktf
        DO i = i_start, i_end
 
 
           dy = 2./(msftx(i,j)+msftx(i,j-1))/rdy  
        
          vel =  v(i,k,j)
          cr = vel*dt/dy       
        

 
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






	
    
      IF(degrade_xs) i_start = its
      IF(degrade_xe) i_end   = MIN(ite,ide-1)
      IF(degrade_ys) j_start = jts  
      IF(degrade_ye) j_end = MIN(jte,jde-1)  
	  


      DO j = j_start, j_end
        DO k=kts,ktf
        DO i = i_start, i_end+1
            dx = 2./(msfty(i,j)+msfty(i-1,j))/rdx 
          
           vel =  u(i,k,j)
         cr = vel*dt/dx 
          fqxl(i,k,j) =  (dx/dt)*flux_upwind((ztop-zs(i-1,j))*field_old(i-1,k,j)/msft(i-1,j), (ztop-zs(i,j))*field_old(i,k,j  )/msft(i,j), cr)
		  
		
		  
		  
          
           
				  
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



   ELSE

      WRITE ( wrf_err_message , * ) 'module_advect: advect_scalar_pd, h_order not known ',horz_order
      CALL wrf_error_fatal3("<stdin>",7578,&
TRIM( wrf_err_message ) )

   ENDIF horizontal_order_test


























      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)



   IF( (config_flags%open_xs) .and. (its == ids) ) THEN

       DO j = j_start, j_end
       DO k = kts, ktf
         ub = MIN( 0.5*(u(its,k,j)+u(its+1,k,j)), 0. )
         tendency(its,k,j) = tendency(its,k,j)                     &
             
             
             
             
             
                 - (ztop-zs(its,j))/msft(its,j)*rdx*(                                             &
                       ub*(   field_old(its+1,k,j)                 &
                            - field_old(its  ,k,j)   ) +           &
                       field_old(its,k,j)*(u(its+1,k,j)-u(its,k,j))  &       
                                                                ) * msft(its,j)*msft(its,j)/(ztop-zs(its,j))  																	
																
																
       ENDDO
       ENDDO

   ENDIF

   IF( (config_flags%open_xe) .and. (ite == ide) ) THEN

       DO j = j_start, j_end
       DO k = kts, ktf
         ub = MAX( 0.5*(u(ite-1,k,j)+u(ite,k,j)), 0. )
         tendency(i_end,k,j) = tendency(i_end,k,j)                   &
        
         
         
         
         
																	
                - (ztop-zs(i_end,j))/msft(i_end,j)*rdx*(                                               &
                       ub*(  field_old(i_end  ,k,j)                  &
                           - field_old(i_end-1,k,j) ) +              &
                       field_old(i_end,k,j)*(u(ite,k,j)-u(ite-1,k,j))  &   
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
                       field_old(i,k,jts)*(v(i,k,jts+1)-v(i,k,jts))  &     
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
                       field_old(i,k,j_end)*(v(i,k,jte)-v(i,k,jte-1))  &    
                                                                    )* msft(i,j_end)*msft(i,j_end)/(ztop-zs(i,j_end)) 
       ENDDO
       ENDDO

   ENDIF

 

   

  

  
 

  























      i_start = its-1
      i_end   = MIN(ite,ide-1)+1
      j_start = jts-1
      j_end   = MIN(jte,jde-1)+1



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
         
            vel = sigmdot(i,k,j)  
           cr = vel*dt/dz  
         
        
         
		  fqzl(i,k,j) = (dz/dt)*flux_upwind(field_old(i,k-1,j), field_old(i,k,j  ), cr)  *(ztop-zs(i,j))/msft(i,j)   

		
	   fqz(i,k,j) = sigmdot(i,k,j) *(fzm(k)*field(i,k,j)+fzp(k)*field(i,k-1,j))*(ztop-zs(i,j))/msft(i,j)    
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
      CALL wrf_error_fatal3("<stdin>",7790,&
wrf_err_message )

   ENDIF vert_order_test
   
   
   
   
   
   
   
   
   
   
   
   
   
   






   IF (mono_limit) THEN



   i_start = its-1
   i_end   = MIN(ite,ide-1)+1
   j_start = jts-1
   j_end   = MIN(jte,jde-1)+1



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



  IF(degrade_xs) i_start = MAX(its,ids+1)
  IF(degrade_xe) i_end   = MIN(ite,ide-2)

  DO j = j_start, j_end
  DO k = kts, ktf
!DEC$ vector always  
  DO i = i_start, i_end

     
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



  i_start = its
  i_end   = MIN(ite,ide-1)
  IF(degrade_ys) j_start = MAX(jts,jds+1)
  IF(degrade_ye) j_end   = MIN(jte,jde-2)

  DO j = j_start, j_end
  DO k = kts, ktf
!DEC$ vector always
  DO i = i_start, i_end

     
     
     tendency (i,k,j) = tendency(i,k,j)                           &
               -  msft(i,j)*msft(i,j)/(ztop-zs(i,j))*( rdy*( fqy (i,k,j+1)-fqy (i,k,j)     &
                                 +fqyl(i,k,j+1)-fqyl(i,k,j)       &
								   ))

  ENDDO
  ENDDO
  ENDDO
  

  

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
   
   
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var,u,v,    &
                                                                      pfm_base,rhofm_base, &
                                                                      sgmdot
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: divg

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh,fzm,fzp,dnw


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx
         
        
       
	   
	   
	   
	   
	   
		 
		    REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) ::  varw
		 
		 
		 
		 


      REAL , DIMENSION( kms:kme ) :: sgmh,dsgm
       REAL , DIMENSION( ims:ime , jms:jme )  ::  rmsft,rmm
	   
	   REAL :: ub,vb 
				 
				 
		
				 
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


 



	  

   
      CALL set_physical_bc3d( var, 'p', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              its, ite, jts, jte, kts, kte,       &
                              its, ite, jts, jte, kts, kte )
	  
	  

     divg = 0.0 
	 
	 

	 

 

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)


      IF ( .NOT. config_flags%periodic_x ) THEN
        IF ( config_flags%open_xs .or. specified ) i_start = MAX(ids+1,its)
        IF ( config_flags%open_xe .or. specified ) i_end   = MIN(ide-2,ite)
      ENDIF

      
	  
	  
	IF ( config_flags%open_xs  .and.  its == ids ) THEN 
		
	
		
      DO  j = j_start, j_end
      DO k=kts,ktf
   
 					 
	
	         ub = MIN( 0.5*(u(its,k,j)+u(its+1,k,j)), 0. )
		 
         divg(its,k,j) =                      &
               - (ztop-zs(its,j))/msft(its,j)*rdx*(                                             &
                       ub*(   var(its+1,k,j)                 &
                            - var(its  ,k,j)   ) +           &
                       var(its,k,j)*(u(its+1,k,j)-u(its,k,j))  &       
                                                                )				 
					 
      ENDDO
      ENDDO
 
	  
	ENDIF	  


	  
	    IF ( config_flags%open_xe  .and.  ite == ide ) THEN 
		
	 
		
		 i_end   = MIN(ite,ide-1)
		 
		
      DO  j = j_start, j_end
      DO k=kts,ktf
 
  
     
	     ub = MAX( 0.5*(u(ite-1,k,j)+u(ite,k,j)), 0. )	 
		 
         divg(i_end,k,j) =                  &
               - (ztop-zs(i_end,j))/msft(i_end,j)*rdx*(                                               &
                       ub*(  var(i_end  ,k,j)                  &
                           - var(i_end-1,k,j) ) +              &
                       var(i_end,k,j)*(u(ite,k,j)-u(ite-1,k,j))  &   
                                                                    )				 
					 
					 
      ENDDO
      ENDDO
  

	ENDIF	
	  
	  
	
	  



      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
      
      IF ( config_flags%open_ys .or. specified .or. config_flags%polar ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified .or. config_flags%polar ) j_end   = MIN(jde-2,jte) 	  

       
	  
	    IF ( config_flags%open_ys  .and.  jts == jds ) THEN 
		
	
		
      DO k=kts,ktf
      DO  i = i_start, i_end		

	
	      vb = MIN( 0.5*(v(i,k,jts)+v(i,k,jts+1)), 0. )
		  
        divg(i,k,jts) =    divg(i,k,jts)                &
               - (ztop-zs(i,jts))/msft(i,jts)*rdy*(                                             &
                       vb*(  var(i,k,jts+1)                  &
                           - var(i,k,jts  ) ) +              &
                       var(i,k,jts)*(v(i,k,jts+1)-v(i,k,jts))  &     
                                                               )
					 
      ENDDO
      ENDDO
	    ENDIF	
		
 
	  
	    IF ( config_flags%open_ye  .and.  jte == jde ) THEN 
		
	 
		
		j_end   = MIN(jte,jde-1)
		
      DO k=kts,ktf
      DO  i = i_start, i_end		
		
	 
	
		  vb = MAX( 0.5*(v(i,k,jte-1)+v(i,k,jte)), 0. ) 
		 
         divg(i,k,j_end) =    divg(i,k,j_end)               &
               -  (ztop-zs(i,j_end))/msft(i,j_end)*rdy*(                                               &
                       vb*(   var(i,k,j_end  )                 &
                            - var(i,k,j_end-1) ) +             &
                       var(i,k,j_end)*(v(i,k,jte)-v(i,k,jte-1))  &    
                                                                    )
					 
      ENDDO
      ENDDO
	    ENDIF
		
		









	 
      field  = var  	 
	 
	 












	  


      i_start = its-1
      i_end   = MIN(ite,ide-1)+1
      j_start = jts-1
      j_end   = MIN(jte,jde-1)+1



    IF(degrade_xs) i_start = MAX(its-1,ids)
    IF(degrade_xe) i_end   = MIN(ite+1,ide-1)
    IF(degrade_ys) j_start = MAX(jts-1,jds)
    IF(degrade_ye) j_end   = MIN(jte+1,jde-1)

    

      DO j = j_start, j_end

         DO i = i_start, i_end
           fqz(i,1,j)  = 0.
        
           fqz(i,kde,j)  = 0.
          
         ENDDO

         DO k=kts+1,ktf
         DO i = i_start, i_end


		   
	   fqz(i,k,j) = sgmdot(i,k,j) *(fzm(k)*field(i,k,j)+fzp(k)*field(i,k-1,j))*(ztop-zs(i,j))/msft(i,j)    
		   

        ENDDO
        ENDDO

      ENDDO

    
   




  
	  
	  	call average_vertical_pw( varw, var,pfm_base,rhofm_base,                 &
                     config_flags,                  &
                         fzm,fzp,dnw,                    &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )	



  
  
  
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



      ktf=MIN(kte,kde-1)
      i_start = its-1
      i_end   = MIN(ite,ide-1)+1
      j_start = jts-1
      j_end   = MIN(jte,jde-1)+1



      IF(degrade_xs) i_start = its
      IF(degrade_xe) i_end   = MIN(ite,ide-1)
      IF(degrade_ys) j_start = MAX(jts,jds+1)
      IF(degrade_ye) j_end = MIN(jte,jde-2)


      DO j = j_start, j_end+1
        DO k=kts,ktf
        DO i = i_start, i_end
 
		   
	
    

				  
	 	    fqy(i,k,j) =    0.5*((ztop-zs(i,j-1))*var(i,k,j-1)*rmsft(i,j-1)  + (ztop-zs(i,j))*var(i,k,j)*rmsft(i,j) ) *v(i,k,j) 	  
		   
        ENDDO
        ENDDO
      ENDDO
	  
	  
	  
	  j_start = jts 
      j_end   = MIN(jte,jde-1) 
	  
	  
      DO j = j_start, j_end
        DO k=kts,ktf
        DO i = i_start, i_end+1
 
				  
      
      

		 	 fqx( i,k,j ) =  	  0.5*((ztop-zs(i-1,j))*var(i-1,k,j)/msft(i-1,j)  + (ztop-zs(i,j))*var(i,k,j)/msft(i,j) ) *u(i,k,j)   
 
        ENDDO
        ENDDO
      ENDDO
	  
	  



	 
	 
	 

 

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
	  
	



      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
      
      IF ( config_flags%open_ys .or. specified .or. config_flags%polar ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified .or. config_flags%polar ) j_end   = MIN(jde-2,jte) 	  

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end

    
    
					 
 		 
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
   
   
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var,u,v,    &
                                                                       rhofm_save, &
                                                                      sgmdot
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: divg

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh,fzm,fzp,dnw


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx
         
        
       
	   
	   
	   
	   
	   
		 
		    REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) ::  varw
		 
		 
		 
		 


      REAL , DIMENSION( kms:kme ) :: sgmh,dsgm
       REAL , DIMENSION( ims:ime , jms:jme )  ::  rmsft,rmm
	   
	   REAL :: ub,vb 
				 
				 
		
				 
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


 



	  

   
      CALL set_physical_bc3d( var, 'p', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              its, ite, jts, jte, kts, kte,       &
                              its, ite, jts, jte, kts, kte )
	  
	  

     divg = 0.0 
	 
	 

	 

 

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)


      IF ( .NOT. config_flags%periodic_x ) THEN
        IF ( config_flags%open_xs .or. specified ) i_start = MAX(ids+1,its)
        IF ( config_flags%open_xe .or. specified ) i_end   = MIN(ide-2,ite)
      ENDIF

      
	  
	  
	IF ( config_flags%open_xs  .and.  its == ids ) THEN 
		
	
		
      DO  j = j_start, j_end
      DO k=kts,ktf
   
 					 
	
	         ub = MIN( 0.5*(u(its,k,j)+u(its+1,k,j)), 0. )
		 
         divg(its,k,j) =                      &
               - (ztop-zs(its,j))/msft(its,j)*rdx*(                                             &
                       ub*(   var(its+1,k,j)                 &
                            - var(its  ,k,j)   ) +           &
                       var(its,k,j)*(u(its+1,k,j)-u(its,k,j))  &       
                                                                )				 
					 
      ENDDO
      ENDDO
 
	  
	ENDIF	  


	  
	    IF ( config_flags%open_xe  .and.  ite == ide ) THEN 
		
	 
		
		 i_end   = MIN(ite,ide-1)
		 
		
      DO  j = j_start, j_end
      DO k=kts,ktf
 
  
     
	     ub = MAX( 0.5*(u(ite-1,k,j)+u(ite,k,j)), 0. )	 
		 
         divg(i_end,k,j) =                  &
               - (ztop-zs(i_end,j))/msft(i_end,j)*rdx*(                                               &
                       ub*(  var(i_end  ,k,j)                  &
                           - var(i_end-1,k,j) ) +              &
                       var(i_end,k,j)*(u(ite,k,j)-u(ite-1,k,j))  &   
                                                                    )				 
					 
					 
      ENDDO
      ENDDO
  

	ENDIF	
	  
	  
 


      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
      
      IF ( config_flags%open_ys .or. specified .or. config_flags%polar ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified .or. config_flags%polar ) j_end   = MIN(jde-2,jte) 	  

       
	  
	    IF ( config_flags%open_ys  .and.  jts == jds ) THEN 
		
	
		
      DO k=kts,ktf
      DO  i = i_start, i_end		

	
	      vb = MIN( 0.5*(v(i,k,jts)+v(i,k,jts+1)), 0. )
		  
        divg(i,k,jts) =    divg(i,k,jts)                &
               - (ztop-zs(i,jts))/msft(i,jts)*rdy*(                                             &
                       vb*(  var(i,k,jts+1)                  &
                           - var(i,k,jts  ) ) +              &
                       var(i,k,jts)*(v(i,k,jts+1)-v(i,k,jts))  &     
                                                               )
					 
      ENDDO
      ENDDO
	    ENDIF	
		

	  
	    IF ( config_flags%open_ye  .and.  jte == jde ) THEN 
		
	 
		
		j_end   = MIN(jte,jde-1)
		
      DO k=kts,ktf
      DO  i = i_start, i_end		
		
	 
	
		  vb = MAX( 0.5*(v(i,k,jte-1)+v(i,k,jte)), 0. ) 
		 
         divg(i,k,j_end) =    divg(i,k,j_end)               &
               -  (ztop-zs(i,j_end))/msft(i,j_end)*rdy*(                                               &
                       vb*(   var(i,k,j_end  )                 &
                            - var(i,k,j_end-1) ) +             &
                       var(i,k,j_end)*(v(i,k,jte)-v(i,k,jte-1))  &    
                                                                    )
					 
      ENDDO
      ENDDO
	    ENDIF
		
		
 






	 
      field  = var  	 
	 
	 












	  


      i_start = its-1
      i_end   = MIN(ite,ide-1)+1
      j_start = jts-1
      j_end   = MIN(jte,jde-1)+1



    IF(degrade_xs) i_start = MAX(its-1,ids)
    IF(degrade_xe) i_end   = MIN(ite+1,ide-1)
    IF(degrade_ys) j_start = MAX(jts-1,jds)
    IF(degrade_ye) j_end   = MIN(jte+1,jde-1)

    

      DO j = j_start, j_end

         DO i = i_start, i_end
           fqz(i,1,j)  = 0.
        
           fqz(i,kde,j)  = 0.
          
         ENDDO

         DO k=kts+1,ktf
         DO i = i_start, i_end


		   
	   fqz(i,k,j) = sgmdot(i,k,j) *(fzm(k)*field(i,k,j)+fzp(k)*field(i,k-1,j))*(ztop-zs(i,j))/msft(i,j)    
		   

        ENDDO
        ENDDO

      ENDDO

    
   




  
	  
	  	call average_vertical_pws( varw, var, rhofm_save,                 &
                     config_flags,                  &
                         fzm,fzp,dnw,                    &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )	



  
  
  
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



      ktf=MIN(kte,kde-1)
      i_start = its-1
      i_end   = MIN(ite,ide-1)+1
      j_start = jts-1
      j_end   = MIN(jte,jde-1)+1



      IF(degrade_xs) i_start = its
      IF(degrade_xe) i_end   = MIN(ite,ide-1)
      IF(degrade_ys) j_start = MAX(jts,jds+1)
      IF(degrade_ye) j_end = MIN(jte,jde-2)


      DO j = j_start, j_end+1
        DO k=kts,ktf
        DO i = i_start, i_end
 
		   
	
    

				  
	 	    fqy(i,k,j) =    0.5*((ztop-zs(i,j-1))*var(i,k,j-1)*rmsft(i,j-1)  + (ztop-zs(i,j))*var(i,k,j)*rmsft(i,j) ) *v(i,k,j) 	  
		   
        ENDDO
        ENDDO
      ENDDO
	  
	  
	  
	  j_start = jts 
      j_end   = MIN(jte,jde-1) 
	  
	  
      DO j = j_start, j_end
        DO k=kts,ktf
        DO i = i_start, i_end+1
 
				  
      
      

		 	 fqx( i,k,j ) =  	  0.5*((ztop-zs(i-1,j))*var(i-1,k,j)/msft(i-1,j)  + (ztop-zs(i,j))*var(i,k,j)/msft(i,j) ) *u(i,k,j)   
 
        ENDDO
        ENDDO
      ENDDO
	  
	  



	 
	 
	 

 

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
	  
	



      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
      
      IF ( config_flags%open_ys .or. specified .or. config_flags%polar ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified .or. config_flags%polar ) j_end   = MIN(jde-2,jte) 	  

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end

    
    
					 
 		 
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
   
   
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var,u,     &
                                                                      v, &
                                                                      sgmdot
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: divg

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh,fzm,fzp


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx
         
        
       
	   
	   
	   
	   
	   
		 
		    REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) ::  varw
		 
		 
		 
		 


      REAL , DIMENSION( kms:kme ) :: sgmh,dsgm
       REAL , DIMENSION( ims:ime , jms:jme )  ::  rmsft,rmm
	   
	   REAL :: ub,vb 
				 
				 
		
				 
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


 
   fqx = 0.0 
   fqy = 0.0 
   fqz = 0.0 
	  

   
      CALL set_physical_bc3d( var, 'p', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              its, ite, jts, jte, kts, kte,       &
                              its, ite, jts, jte, kts, kte )
	  
	  

     divg = 0.0 
	 
	 

	 

 

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)


      IF ( .NOT. config_flags%periodic_x ) THEN
        IF ( config_flags%open_xs .or. specified ) i_start = MAX(ids+1,its)
        IF ( config_flags%open_xe .or. specified ) i_end   = MIN(ide-2,ite)
      ENDIF

      
	  
	  
	IF ( config_flags%open_xs  .and.  its == ids ) THEN 
		
	
		
      DO  j = j_start, j_end
      DO k=kts,ktf
   
 					 
	
	         ub = MIN( 0.5*(u(its,k,j)+u(its+1,k,j)), 0. )
		 
         divg(its,k,j) =                      &
               - (1.0-zs(its,j)/ztop)/msft(its,j)*rdx*(                                             &
                       ub*(   var(its+1,k,j)                 &
                            - var(its  ,k,j)   ) +           &
                       var(its,k,j)*(u(its+1,k,j)-u(its,k,j))  &       
                                                                )				 
					 
      ENDDO
      ENDDO
 
	  
	ENDIF	  


	  
	    IF ( config_flags%open_xe  .and.  ite == ide ) THEN 
		
	 
		
		 i_end   = MIN(ite,ide-1)
		 
		
      DO  j = j_start, j_end
      DO k=kts,ktf
 
  
     
	     ub = MAX( 0.5*(u(ite-1,k,j)+u(ite,k,j)), 0. )	 
		 
         divg(i_end,k,j) =                  &
               - (1.0-zs(i_end,j)/ztop)/msft(i_end,j)*rdx*(                                               &
                       ub*(  var(i_end  ,k,j)                  &
                           - var(i_end-1,k,j) ) +              &
                       var(i_end,k,j)*(u(ite,k,j)-u(ite-1,k,j))  &   
                                                                    )				 
					 
					 
      ENDDO
      ENDDO
  

	ENDIF	
	  
	  
 
	  



      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
      
      IF ( config_flags%open_ys .or. specified .or. config_flags%polar ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified .or. config_flags%polar ) j_end   = MIN(jde-2,jte) 	  

       
	  
	    IF ( config_flags%open_ys  .and.  jts == jds ) THEN 
		
	
		
      DO k=kts,ktf
      DO  i = i_start, i_end		

	
	      vb = MIN( 0.5*(v(i,k,jts)+v(i,k,jts+1)), 0. )
		  
        divg(i,k,jts) =    divg(i,k,jts)                &
               - (1.0-zs(i,jts)/ztop)/msft(i,jts)*rdy*(                                             &
                       vb*(  var(i,k,jts+1)                  &
                           - var(i,k,jts  ) ) +              &
                       var(i,k,jts)*(v(i,k,jts+1)-v(i,k,jts))  &     
                                                               )
					 
      ENDDO
      ENDDO
	    ENDIF	
		

	  
	    IF ( config_flags%open_ye  .and.  jte == jde ) THEN 
		
	 
		
		j_end   = MIN(jte,jde-1)
		
      DO k=kts,ktf
      DO  i = i_start, i_end		
		
	 
	
		  vb = MAX( 0.5*(v(i,k,jte-1)+v(i,k,jte)), 0. ) 
		 
         divg(i,k,j_end) =    divg(i,k,j_end)               &
               -  (1.0-zs(i,j_end)/ztop)/msft(i,j_end)*rdy*(                                               &
                       vb*(   var(i,k,j_end  )                 &
                            - var(i,k,j_end-1) ) +             &
                       var(i,k,j_end)*(v(i,k,jte)-v(i,k,jte-1))  &    
                                                                    )
					 
      ENDDO
      ENDDO
	    ENDIF
		
 




	 
      field  = var  	 
	 
	 












	  


      i_start = its-1
      i_end   = MIN(ite,ide-1)+1
      j_start = jts-1
      j_end   = MIN(jte,jde-1)+1



    IF(degrade_xs) i_start = MAX(its-1,ids)
    IF(degrade_xe) i_end   = MIN(ite+1,ide-1)
    IF(degrade_ys) j_start = MAX(jts-1,jds)
    IF(degrade_ye) j_end   = MIN(jte+1,jde-1)

    

      DO j = j_start, j_end
      DO i = i_start, i_end
		 
           fqz(i,1,j)  = 0.
        
           fqz(i,kde,j)  = 0.
          
         ENDDO

         DO k=kts+1,ktf
         DO i = i_start, i_end


		   
	   fqz(i,k,j) = sgmdot(i,k,j) *(fzm(k)*field(i,k,j)+fzp(k)*field(i,k-1,j))*(1.0-zs(i,j)/ztop)/msft(i,j)    
		   

        ENDDO
        ENDDO

      ENDDO

    
   




  
	  
	   	call average_vertical_w( varw, var,                 &
                      config_flags,                  &
                         fzm,fzp,                    &
                         ids, ide, jds, jde, kds, kde,  &
                         ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )	

	  
      
      
       
      
       

  
  
  
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


           
             divg(i,k,j) = divg(i,k,j) -    ( varw(i,k+1,j)* rmm(i,j) * sgmdot(i,k+1,j) - varw(i,k,j)* rmm(i,j) * sgmdot(i,k,j)   ) /dsgm(k) * (1.0-zs(i,j)/ztop)

            
									  
			 

  

			 
									  
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



      ktf=MIN(kte,kde-1)
      i_start = its-1
      i_end   = MIN(ite,ide-1)+1
      j_start = jts-1
      j_end   = MIN(jte,jde-1)+1



      IF(degrade_xs) i_start = its
      IF(degrade_xe) i_end   = MIN(ite,ide-1)
      IF(degrade_ys) j_start = MAX(jts,jds+1)
      IF(degrade_ye) j_end = MIN(jte,jde-1) 


      DO j = j_start, j_end 
        DO k=kts,ktf
        DO i = i_start, i_end
 
		   
	 

				  
	 	    fqy(i,k,j) =    0.5*((1.0-zs(i,j-1)/ztop)*var(i,k,j-1)*rmsft(i,j-1)  + (1.0-zs(i,j)/ztop)*var(i,k,j)*rmsft(i,j) ) *v(i,k,j) 	  
		   
        ENDDO
        ENDDO
      ENDDO
	  

 
  

  
	  
	  j_start = jts 
      j_end   = MIN(jte,jde-1) 
	  
	  
      DO j = j_start, j_end
  
 
      DO k=kts,ktf
        DO i = i_start, i_end
 
	 	

		 	 fqx( i,k,j ) =  	  0.5*((1.0-zs(i-1,j)/ztop)*var(i-1,k,j)/msft(i-1,j)  + (1.0-zs(i,j)/ztop)*var(i,k,j)/msft(i,j) ) *u(i,k,j)   
 
  
 
        ENDDO
        ENDDO
      ENDDO
	  
 
 
	 
	 
	 

 

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
	  


 
 
	



      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
      
      IF ( config_flags%open_ys .or. specified .or. config_flags%polar ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified .or. config_flags%polar ) j_end   = MIN(jde-2,jte) 	  

      DO  j = j_start, j_end
 
      DO k=kts,ktf
      DO  i = i_start, i_end

 
 		 
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
         	divg(i,k,j) = divg(i,k,j)  * msft(i,j)*msft(i,j)/(1.0-zs(i,j)/ztop) 
      ENDDO
      ENDDO
      ENDDO
	  
 
 
	  


 
 
 

      END SUBROUTINE calcu_divg_pres  
	  

 SUBROUTINE average_vertical_w_log( varw, var,                 &
                     config_flags,                  &
                         fzm,fzp,                    &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var

        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: varw
		 
		 	  REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) ::  fzm,fzp
			  
			  


   
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx
         
     




  
    LOGICAL :: specified


	specified = .false.
	if(config_flags%specified .or. config_flags%nested) specified = .true.
   
     ktf=MIN(kte,kde-1)
 
 


      i_start = its
      i_end   = MIN(ide-1,ite)



      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
	  

      varw  = 0. 
	  
	  


    
    

    
    
   

    
    


    
    
	  
	  

      DO  j = j_start, j_end
      DO k=kts+1,ktf
      DO  i = i_start, i_end
       
      	varw(i,k,j) = exp(  fzm(k)*log(var(i,k,j) ) + fzp(k)*log(var(i,k-1,j) ) )     
      ENDDO
      ENDDO
      ENDDO


      DO      j = j_start, j_end
  
      DO  i = i_start, i_end
       varw(i,kts,j) =   2.0* var(i,kts,j)  -  varw(i,kts+1,j)  
	   
      varw(i,kte,j) =  2.0* var(i,kte-1,j)  -  varw(i,kte-1,j)   
      ENDDO
      ENDDO
	   
	  

      END SUBROUTINE average_vertical_w_log 	







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
       
  
    i_start = MAX(its-1,ids)
    i_end   = MIN(ite+1,ide-1)
    j_start = MAX(jts-1,jds)
    j_end   = MIN(jte+1,jde-1) 
	  

	   
	   
	   









	  
       engk  = 0.0  
	  

      DO      j = j_start, j_end
      DO  k = k_start, k_end
      DO  i = i_start, i_end
         
			  engk(i,k,j) = ( v(i,k,j)**2 + v(i,k,j+1)**2 + w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25 
      ENDDO
      ENDDO
      ENDDO
	  
	  
if(.false.) then 
	  
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
	  

	
      IF ( specified .AND. its .LE. ids+1 .AND. .NOT. config_flags%periodic_x .AND. jts .LE. jds+1  ) THEN
	 	 
	    i = ids 
        j = jds 
		
        DO k=kts,ktf
 

   
           vb = v(i,k,j)
	
           IF (v(i,k,j+1) .LT. 0.) vb = v(i,k,j+1)   

              engk(i,k,j) = ( vb**2 + v(i,k,j+1)**2 + w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25  

 
        ENDDO
      ENDIF	  
	  

    IF ( specified .AND. ite .GE. ide-1 .AND. .NOT. config_flags%periodic_x .AND. jts .LE. jds+1  ) THEN
	
	
       j = jds 
        DO k=kts,ktf
           i =  ide-1  
         
     
		     vb = v(i,k,j)
			 
   
		   IF (v(i,k,j+1) .LT. 0.) vb = v(i,k,j+1)    
		  
		  
            engk(i,k,j) = (    vb**2 + v(i,k,j+1)**2 + w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25 
        ENDDO
      
		
      ENDIF		  
	  
	  
	  

	
      IF ( specified .AND. its .LE. ids+1 .AND. .NOT. config_flags%periodic_x .AND. jte .GE. jde-1  ) THEN
	 	 
	    i = ids 
        j = jde-1 
		
        DO k=kts,ktf
 

        
		   
		   
		   vb = v(i,k,j+1)
           IF (v(i,k,j) .GT. 0.) vb = v(i,k,j)
		   
		   

              engk(i,k,j) = ( v(i,k,j)**2 + vb**2 +  w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25  

 
        ENDDO
      ENDIF	  	  



    IF ( specified .AND. ite .GE. ide-1 .AND. .NOT. config_flags%periodic_x .AND. jte .GE. jde-1  ) THEN
	
	
       j = jde-1 
        DO k=kts,ktf
           i =  ide-1  
         
          
 


		   vb = v(i,k,j+1)
           IF (v(i,k,j) .GT. 0.) vb = v(i,k,j) 		   
		  
		  
            engk(i,k,j) = (   v(i,k,j+1)**2 + w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25 
        ENDDO
    
		
      ENDIF		  	 
	  
endif 
	  
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
       
  
    i_start = MAX(its-1,ids)
    i_end   = MIN(ite+1,ide-1)
    j_start = MAX(jts-1,jds)
    j_end   = MIN(jte+1,jde-1) 
	  

	   
	   
	   






	  


	  
       engk  = 0.0  
	  

      DO      j = j_start, j_end
      DO  k = k_start, k_end
      DO  i = i_start, i_end
  
			  engk(i,k,j) = ( u(i,k,j)**2 + u(i+1,k,j)**2   + w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25 
      ENDDO
      ENDDO
      ENDDO
	  
	  
	  
if(.false.) then 
	  
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
	  

	
      IF ( specified .AND. its .LE. ids+1 .AND. .NOT. config_flags%periodic_x .AND. jts .LE. jds+1  ) THEN
	 	 
	    i = ids 
        j = jds 
		
        DO k=kts,ktf
 

           ub = u(i ,k,j)
          
		   IF (u(i+1,k,j) .LT. 0.) ub = u(i+1,k,j)
         

              engk(i,k,j) = ( ub**2 + u(i+1,k,j)**2 +  w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25  

 
        ENDDO
      ENDIF	  
	  

    IF ( specified .AND. ite .GE. ide-1 .AND. .NOT. config_flags%periodic_x .AND. jts .LE. jds+1  ) THEN
	
	
       j = jds 
        DO k=kts,ktf
           i =  ide-1  
         
           ub = u(i+1 ,k,j)
	 
			 
          IF (u(i,k,j) .GT. 0.) ub = u(i,k,j) 
  
		  
		  
            engk(i,k,j) = (  u(i,k,j)**2 + ub**2 +    w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25 
        ENDDO
     
		
      ENDIF		  
	  
	  
	  

	
      IF ( specified .AND. its .LE. ids+1 .AND. .NOT. config_flags%periodic_x .AND. jte .GE. jde-1  ) THEN
	 	 
	    i = ids 
        j = jde-1 
		
        DO k=kts,ktf
 
           ub = u(i ,k,j)
       
		   IF (u(i+1,k,j) .LT. 0.) ub = u(i+1,k,j)
     
              engk(i,k,j) = ( ub**2 + u(i+1,k,j)**2 +   w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25  
 
        ENDDO
      ENDIF	  	  



    IF ( specified .AND. ite .GE. ide-1 .AND. .NOT. config_flags%periodic_x .AND. jte .GE. jde-1  ) THEN
	
	
       j = jde-1 
        DO k=kts,ktf
           i =  ide-1  
         
           ub = u(i+1 ,k,j)
 	 
          IF (u(i,k,j) .GT. 0.) ub = u(i,k,j) 
 
            engk(i,k,j) = (  u(i,k,j)**2 + ub**2 +    w(i,k,j)**2 + w(i,k+1,j)**2 ) * 0.25 
        ENDDO
    
		
      ENDIF		  	 
	  

endif 
	  
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
       
  
    i_start = MAX(its-1,ids)
    i_end   = MIN(ite+1,ide-1)
    j_start = MAX(jts-1,jds)
    j_end   = MIN(jte+1,jde-1) 
	  

	   
	   
	   






	  


	  
       engk  = 0.0 

      DO      j = j_start, j_end
      DO  k = k_start, k_end
      DO  i = i_start, i_end
    
			  engk(i,k,j) = ( u(i,k,j)**2 + u(i+1,k,j)**2 + v(i,k,j)**2 + v(i,k,j+1)**2   ) * 0.25 
      ENDDO
      ENDDO
      ENDDO
	  
	  
	  
if(.false.) then 
	  
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

	  
endif 
	  
	  
    END SUBROUTINE calcu_engk_uv




     SUBROUTINE calcu_grady_new( var0, grady,                 &
 
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   
   

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

   
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdy
         
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: varz  
		   REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: var 

       REAL , DIMENSION( kms:kme ) :: sgmh,dsgm,dsgmh
    LOGICAL :: specified

   
   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

 
      ktf=MIN(kte,kde-1)


 

	  DO k=kts,kde-1
          sgmh(k) = 0.5*(zh(k)+zh(k+1))
          dsgm(k) = zh(k+1) - zh(k) 
      ENDDO

 	  DO k=2,kde-1
             dsgmh(k) = 0.5*(dsgm(k) + dsgm(k-1))   
			 
      ENDDO



	 
       
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
							  
							  

      IF ( config_flags%periodic_y ) j_start = jts            
      IF ( config_flags%periodic_y ) j_end = jte 


       varz = 0.0 

      DO  j = j_start, j_end
      DO k=2,ktf
      DO  i = i_start, i_end
    
      	varz(i,k,j) =  (var(i,k,j) - var(i,k-1,j) ) / dsgmh(k)
      ENDDO
      ENDDO
      ENDDO
	  
	  


      DO      j = j_start, j_end
  
      DO  i = i_start, i_end
		varz(i,kts,j) =  varz(i,kts+1,j)  
	   
		varz(i,kte,j) =  varz(i,kte-1,j)  
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
	  



 
     
    
    
	
      i_start = its
      i_end   =  MIN(ite,ide-1)
      j_start = jts
      j_end   = jte   
	  
      IF ( config_flags%periodic_y ) j_start = jts          
      IF ( config_flags%periodic_y ) j_end = jte



    IF ( config_flags%open_ys .or. specified ) j_start = MAX(jds+1,jts)
    IF ( config_flags%open_ye .or. specified ) j_end   = MIN(jde-1,jte)
	  

      grady = 0.0 
  
      DO  j = j_start, j_end
      DO  k = kts, ktf
      DO  i = i_start, i_end 
	  
	  
	     mrdy= 0.5* ( msft(i,j-1) + msft(i,j) ) *rdy
	   
	   
	   
		 
	    
		 
		    c2  =   (   msft(i,j)  *(ztop-sgmh(k))/(ztop-zs(i,j))   )  *   (   varz(i,k,j)  + varz(i,k+1,j)  )               &  
		         +  (   msft(i,j-1)*(ztop-sgmh(k))/(ztop-zs(i,j-1 ))  )*   (   varz(i,k,j-1)  + varz(i,k+1,j-1)  )  
			   
			    c1 = rdy*(zs(i,j)-zs(i,j-1)) 
				
				
		 
	  
	  
	  
	    grady(i,k,j) =       - mrdy*(var(i,k,j)-var(i,k,j-1)) +  c1*c2*0.25
		


		 
	  ENDDO
      ENDDO
      ENDDO    
	  


  
	  

      END SUBROUTINE calcu_grady_new






	
     
      SUBROUTINE calcu_gradx_new( var, gradx,                 &
                     config_flags,                  &
                        msft,zh,zs,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var    
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: gradx

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx            
       
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: varz 

        REAL , DIMENSION( kms:kme ) :: sgmh,dsgm,dsgmh
		
		 REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  :: vardx 
		
          REAL , DIMENSION( ims:ime , kms:kme , jms:jme )  ::  varbc
		
        LOGICAL :: specified


   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

 
      ktf=MIN(kte,kde-1)

 	  DO k=kts,kde-1
          sgmh(k) = 0.5*(zh(k)+zh(k+1))
          dsgm(k) = zh(k+1) - zh(k) 
      ENDDO

 	  DO k=2,kde-1
             dsgmh(k) = 0.5*(dsgm(k) + dsgm(k-1))   
			 
      ENDDO
	  
	  

         varbc  = var
	  
      CALL set_physical_bc3d( varbc, 'p', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              its, ite, jts, jte, kts, kte,       &
                              its, ite, jts, jte, kts, kte )




 
      
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
 
      DO  i = i_start, i_end
       varz(i,kts,j) =  varz(i,kts+1,j)  
      varz(i,kte,j) =  varz(i,ktf,j)  
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

	

      
      i_start = its
      i_end   = ite 
	  
	  
      j_start = jts
      j_end   = MIN(jte,jde-1)	

	
	  
	  


      IF ( config_flags%open_xs ) i_start = MAX(ids+1,its)   
      IF ( config_flags%open_xe ) i_end   = MIN(ide-1,ite)
      IF ( specified ) i_start = MAX(ids+1,its)            
      IF ( specified ) i_end   = MIN(ide-1,ite)
	  
      IF ( config_flags%periodic_x ) i_start = its           
      IF ( config_flags%periodic_x ) i_end = ite



 
 



 
 
 
 
 
 

 

 
 


		
  loop_j :     DO  j = j_start, j_end
	
		 
	 
      DO  k = kts, ktf
      DO  i = i_start, i_end 
	  
	    mrdx= 0.5* ( msft(i-1,j) + msft(i,j) ) *rdx
		
	 
	 
	


		    c2  =   (   msft(i,j)  *(ztop-sgmh(k))/(ztop-zs(i,j))   )  *   (   varz(i,k,j)  + varz(i,k+1,j)  )               &  
		         +  (   msft(i-1,j)*(ztop-sgmh(k))/(ztop-zs(i-1,j ))  )*   (   varz(i-1,k,j)  + varz(i-1,k+1,j)  )  
			   
			    c1 = rdx*(zs(i,j)-zs(i-1,j)) 
				

    
	
	   gradx(i,k,j) =  - mrdx*(varbc(i,k,j)-varbc(i-1,k,j))  + c1*c2
	  

	  ENDDO
      ENDDO	 

      ENDDO    loop_j
	  


      END SUBROUTINE calcu_gradx_new





SUBROUTINE calcu_divg( var, divg,                 &
                        u, v, sgmdot,                       &
                     config_flags,                  &
                        msft,zh,zs,fzm,fzp,                    &
                        rdx,rdy,ztop,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE
   
   
   

            INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: var,u,     &
                                                                      v, &
                                                                      sgmdot
        TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

      
         REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: divg

         REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msft,zs

         REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: zh,fzm,fzp


         REAL ,                                        INTENT(IN   ) :: rdx,  &
                                                                  rdy

         REAL ,                                        INTENT(IN   ) :: ztop

   
   
         INTEGER :: i, j, k, itf, jtf, ktf
         INTEGER :: i_start, i_end, j_start, j_end
         INTEGER :: i_start_f, i_end_f, j_start_f, j_end_f
         INTEGER :: jmin, jmax, jp, jm, imin, imax, im, ip
         INTEGER :: jp1, jp0, jtmp
         REAL :: c1,c2,mrdx
         
        
       
	   
	   
	   
	   
	   
		 
		    REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) ::  varw
		 
		 
		 
		 


      REAL , DIMENSION( kms:kme ) :: sgmh,dsgm
       REAL , DIMENSION( ims:ime , jms:jme )  ::  rmm    
	   
	   REAL :: ub,vb 
				 
				 
		
				 
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


 



	  

   
      CALL set_physical_bc3d( var, 'p', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              its, ite, jts, jte, kts, kte,       &
                              its, ite, jts, jte, kts, kte )
	  
	  

     divg = 0.0 
	 
	 

	 

 

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)


      IF ( .NOT. config_flags%periodic_x ) THEN
        IF ( config_flags%open_xs .or. specified ) i_start = MAX(ids+1,its)
        IF ( config_flags%open_xe .or. specified ) i_end   = MIN(ide-2,ite)
      ENDIF

      
	  
	  
	IF ( config_flags%open_xs  .and.  its == ids ) THEN 
		
	
		
      DO  j = j_start, j_end
      DO k=kts,ktf
   
 					 
	
	         ub = MIN( 0.5*(u(its,k,j)+u(its+1,k,j)), 0. )
		 
         divg(its,k,j) =                      &
               - (ztop-zs(its,j))/msft(its,j)*rdx*(                                             &
                       ub*(   var(its+1,k,j)                 &
                            - var(its  ,k,j)   ) +           &
                       var(its,k,j)*(u(its+1,k,j)-u(its,k,j))  &       
                                                                )				 
					 
      ENDDO
      ENDDO
 
	  
	ENDIF	  


	  
	    IF ( config_flags%open_xe  .and.  ite == ide ) THEN 
		
	 
		
		 i_end   = MIN(ite,ide-1)
		 
		
      DO  j = j_start, j_end
      DO k=kts,ktf
 
  
     
	     ub = MAX( 0.5*(u(ite-1,k,j)+u(ite,k,j)), 0. )	 
		 
         divg(i_end,k,j) =                  &
               - (ztop-zs(i_end,j))/msft(i_end,j)*rdx*(                                               &
                       ub*(  var(i_end  ,k,j)                  &
                           - var(i_end-1,k,j) ) +              &
                       var(i_end,k,j)*(u(ite,k,j)-u(ite-1,k,j))  &   
                                                                    )				 
					 
					 
      ENDDO
      ENDDO
  

	ENDIF	
	  
	  
	
	  



      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
      
      IF ( config_flags%open_ys .or. specified .or. config_flags%polar ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified .or. config_flags%polar ) j_end   = MIN(jde-2,jte) 	  

       
	  
	    IF ( config_flags%open_ys  .and.  jts == jds ) THEN 
		
	
		
      DO k=kts,ktf
      DO  i = i_start, i_end		

	
	      vb = MIN( 0.5*(v(i,k,jts)+v(i,k,jts+1)), 0. )
		  
        divg(i,k,jts) =    divg(i,k,jts)                &
               - (ztop-zs(i,jts))/msft(i,jts)*rdy*(                                             &
                       vb*(  var(i,k,jts+1)                  &
                           - var(i,k,jts  ) ) +              &
                       var(i,k,jts)*(v(i,k,jts+1)-v(i,k,jts))  &     
                                                               )
					 
      ENDDO
      ENDDO
	    ENDIF	
		

	  
	    IF ( config_flags%open_ye  .and.  jte == jde ) THEN 
		
	 
		
		j_end   = MIN(jte,jde-1)
		
      DO k=kts,ktf
      DO  i = i_start, i_end		
		
	 
	
		  vb = MAX( 0.5*(v(i,k,jte-1)+v(i,k,jte)), 0. ) 
		 
         divg(i,k,j_end) =    divg(i,k,j_end)               &
               -  (ztop-zs(i,j_end))/msft(i,j_end)*rdy*(                                               &
                       vb*(   var(i,k,j_end  )                 &
                            - var(i,k,j_end-1) ) +             &
                       var(i,k,j_end)*(v(i,k,jte)-v(i,k,jte-1))  &    
                                                                    )
					 
      ENDDO
      ENDDO
	    ENDIF
		
		









	 
      field  = var  	 
	 
	 












	  


      i_start = its 
      i_end   = MIN(ite,ide-1) 
      j_start = jts 
      j_end   = MIN(jte,jde-1) 



    IF(degrade_xs) i_start = MAX(its-1,ids)
    IF(degrade_xe) i_end   = MIN(ite+1,ide-1)
    IF(degrade_ys) j_start = MAX(jts-1,jds)
    IF(degrade_ye) j_end   = MIN(jte+1,jde-1)

    

      DO j = j_start, j_end

         DO i = i_start, i_end
           fqz(i,1,j)  = 0.
        
           fqz(i,kde,j)  = 0.
          
         ENDDO

         DO k=kts+1,ktf
         DO i = i_start, i_end


		   
	   fqz(i,k,j) = sgmdot(i,k,j) *(fzm(k)*field(i,k,j)+fzp(k)*field(i,k-1,j))*(ztop-zs(i,j))/msft(i,j)    
		   

        ENDDO
        ENDDO

      ENDDO

    
   




  
	  
	  	call average_vertical_w( varw, var,                 &
                     config_flags,                  &
                         fzm,fzp,                    &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte  )	



  
  
  
        i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
	  
      DO  j = j_start, j_end
      DO  i = i_start, i_end
      
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



      ktf=MIN(kte,kde-1)
	  
	  
	    
		
      i_start = its
      i_end   = MIN(ite,ide-1) 
      j_start = jts
      j_end   = MIN(jte+1,jde) 



  
  
  
  


      DO j = j_start, j_end 
        DO k=kts,ktf
        DO i = i_start, i_end
 
		   
	
    

				  
	 
		   
		    fqy(i,k,j) =    0.5*((ztop-zs(i,j-1))*var(i,k,j-1)/msft(i,j-1)  + (ztop-zs(i,j))*var(i,k,j)/msft(i,j) ) *v(i,k,j) 	
        ENDDO
        ENDDO
      ENDDO
	  
	  
	  
      i_start = its
      i_end   = MIN(ite+1,ide) 
      j_start = jts
      j_end   = MIN(jte,jde-1) 
	  
	  
      DO j = j_start, j_end
        DO k=kts,ktf
        DO i = i_start, i_end 
 
				  
      
      

		 	 fqx( i,k,j ) =  	  0.5*((ztop-zs(i-1,j))*var(i-1,k,j)/msft(i-1,j)  + (ztop-zs(i,j))*var(i,k,j)/msft(i,j) ) *u(i,k,j)   
 
 
 
 
        ENDDO
        ENDDO
      ENDDO
	  
	  



	 
	 
	 

 

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
	  



 
	 
	 


      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)
	  
      
      IF ( config_flags%open_ys .or. specified .or. config_flags%polar ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified .or. config_flags%polar ) j_end   = MIN(jde-2,jte) 	  

      DO  j = j_start, j_end
      DO k=kts,ktf
      DO  i = i_start, i_end

    
    
					 
 		 
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

j_end = jts

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





  rhoengk = rho*engk
  potential_energy  = rho*z_at_mass*g

  internal_energy = p*717.5/287.0















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




      SUBROUTINE calcu_engk_test(engk,rho,rhop,u,v,w,config_flags,                       &
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

