!> @file
!> @brief Subroutine that computes the velocity potential and
!> streamfunction from isobaric winds.
!>
!><pre>
!> This routine uses discrete numerical analysis to compute the
!> velocity potential and streamfunction from the wind velocities
!> at the given P-point. The process is based on the CALVOR 
!> subroutine in UPP_PHYSICS, authored by J. Meng. The routine
!> returns: 
!> 	  PSI: the streamfunction for P-point
!> 	  CHI: the velocity potential for P-point
!></pre>
!>   
!> @param[in] UP real U-wind component (m/s) at P-point.
!> @param[in] VP real V-wind component (m/s) at P-point.
!> @param[out] CHI real velocity potential (m^2/s) at P-point.
!> @param[out] PSI real streamfunction (m^2/s) at P-point.
!>
!> ### Program history log:
!> Date       | Programmer   | Comments
!> -----------|--------------|---------
!> 2024-08-27 | Karina Asmar | Initial   
!>
!> @author Karina Asmar EMC/VPPPG @date 2024-08-27
!-----------------------------------------------------------------------
!> @brief Subroutine that computes velocity potential and streamfunction
!> from isobaric winds.  
!>   
!> @param[in] UP real U-wind component (m/s) at P-point.
!> @param[in] VP real V-wind component (m/s) at P-point.
!> @param[out] CHI real velocity potential (m^2/s) at P-point.
!> @param[out] PSI real streamfunction (m^2/s) at P-point.
!-----------------------------------------------------------------------
      SUBROUTINE CALCHIPSI(UP,VP,CHI,PSI)
!
!     INCLUDE ETA GRID DIMENSIONS.  SET/DERIVE OTHER PARAMETERS.
      use masks,        only: gdlat, gdlon
      use params_mod,   only: dtr, small, erad
      use ctlblk_mod,   only: jsta_2l, jend_2u, spval, &
                              jsta, jend, im, jm, &
                              ista, iend, ista_2l, iend_2u
!
      implicit none
!
!     DECLARE VARIABLES.
!     
      REAL, dimension(ista_2l:iend_2u,jsta_2l:jend_2u), intent(in)    :: UP, VP
      REAL, dimension(ista_2l:iend_2u,jsta_2l:jend_2u), intent(out) :: CHI, PSI
      REAL, dimension(IM,2) :: GLATPOLES, COSLPOLES, UPOLES, VPOLES, PSIPOLES, CHIPOLES
      REAL, dimension(IM,JSTA:JEND) :: COSLTEMP, PSITEMP, CHITEMP
!
      real,    allocatable ::  wrk1(:,:), wrk3(:,:), cosl(:,:)
      REAL, dimension(ista_2l:iend_2u,jsta_2l:jend_2u) :: DCHI, DPSI
!
      integer, parameter :: npass2=2, npass3=3
      integer I,J,ii,iir,iil,jj,JMT2,imb2, npass, nn, jtem
      real    tx1(im+2), tx2(im+2), tx3(im+2), tx4(im+2)
!     
!***************************************************************************
!     START CALCHIPSI HERE.
!     
!     LOOP TO COMPUTE STREAMFUNCTION AND VELOCITY POTENTIAL FROM WINDS.
   
!$omp  parallel do private(i,j)
        DO J=JSTA_2L,JEND_2U
          DO I=ISTA_2L,IEND_2U
            DPSI(I,J) = SPVAL
            DCHI(I,J) = SPVAL
          ENDDO
        ENDDO

      CALL EXCH(UP)
      CALL EXCH(VP)
!
        CALL EXCH(GDLAT(ISTA_2L,JSTA_2L))
        CALL EXCH(GDLON(ISTA_2L,JSTA_2L))

        allocate (wrk1(ista:iend,jsta:jend),        &
     &            wrk3(ista:iend,jsta:jend), cosl(ista_2l:iend_2u,jsta_2l:jend_2u))

        imb2 = im/2

!$omp  parallel do private(i,j)
        DO J=JSTA,JEND 
          do i=ista,iend
            cosl(i,j) = cos(gdlat(i,j)*dtr)
            IF(cosl(i,j) >= SMALL) then
              wrk1(i,j) = ERAD*cosl(i,j)
            else
              wrk1(i,j) = 0.
            end if    
          enddo
        enddo
        CALL EXCH(cosl)

        call fullpole( cosl(ista_2l:iend_2u,jsta_2l:jend_2u),coslpoles)
        call fullpole(gdlat(ista_2l:iend_2u,jsta_2l:jend_2u),glatpoles)

!$omp  parallel do private(i,j,ii)
        DO J=JSTA,JEND
          if (j == 1) then
           if(gdlat(ista,j) > 0.) then ! count from north to south
              do i=ista,iend
                ii = i + imb2
                if (ii > im) ii = ii - im
                wrk3(i,j) = (180.-GDLAT(i,J+1)-GLATPOLES(ii,1))*DTR !1/dphi
              enddo
            else ! count from south to north
              do i=ista,iend
                ii = i + imb2
                if (ii > im) ii = ii - im
                wrk3(i,j) = (180.+GDLAT(i,J+1)+GLATPOLES(ii,1))*DTR !1/dphi
              enddo
            end if      
          elseif (j == JM) then
            if(gdlat(ista,j) < 0.) then ! count from north to south
              do i=ista,iend
                ii = i + imb2
                if (ii > im) ii = ii - im
                wrk3(i,j) = (180.+GDLAT(i,J-1)+GLATPOLES(ii,2))*DTR
              enddo
            else ! count from south to north
              do i=ista,iend
                ii = i + imb2
                if (ii > im) ii = ii - im
                wrk3(i,j) = (180.-GDLAT(i,J-1)-GLATPOLES(ii,2))*DTR
              enddo
            end if  
          else
            do i=ista,iend
              wrk3(i,j) = (GDLAT(I,J-1)-GDLAT(I,J+1))*DTR !1/dphi
            enddo
          endif
        enddo  

        npass = 0

        jtem = jm / 18 + 1
      
        call fullpole(UP(ista_2l:iend_2u,jsta_2l:jend_2u),upoles)
        call fullpole(VP(ista_2l:iend_2u,jsta_2l:jend_2u),vpoles)

!$omp  parallel do private(i,j,ii)
        DO J=JSTA,JEND
          if (j == 1) then
           if(gdlat(ista,j) > 0.) then ! count from north to south
              do i=ista,iend
                  if(UP(I,J)==SPVAL) cycle
                  DPSI(I,J) = UP(I,J)*wrk3(i,j) * wrk1(i,j)   
                  if(VP(I,J)==SPVAL) cycle
                  DCHI(I,J) = -1.0*VP(I,J)*wrk3(i,j) * wrk1(i,j)  
              enddo
            else ! count from south to north
              jj=2
              do i=ista,iend
                  if(UP(I,J)==SPVAL) cycle
                  DPSI(I,J) = UP(I,jj)*wrk3(i,jj) * wrk1(i,jj)   
                  if(VP(I,J)==SPVAL) cycle                  
                  DCHI(I,J) = -1.0*VP(I,jj)*wrk3(i,jj) * wrk1(i,jj)  
              enddo
            end if      
          elseif (j == JM) then
            if(gdlat(ista,j) < 0.) then ! count from north to south
              do i=ista,iend
                  if(UP(I,J)==SPVAL) cycle              
                  DPSI(I,J) = UP(I,J)*wrk3(i,j) * wrk1(i,j)   
                  if(VP(I,J)==SPVAL) cycle            
                  DCHI(I,J) = -1.0*VP(I,J)*wrk3(i,j) * wrk1(i,j)
              enddo
            else ! count from south to north
              jj = jm-1
              do i=ista,iend
                  if(UP(I,J)==SPVAL) cycle              
                  DPSI(I,J) = UP(I,jj)*wrk3(i,jj) * wrk1(i,jj)  
                  if(VP(I,J)==SPVAL) cycle                  
                  DCHI(I,J) = -1.0*VP(I,jj)*wrk3(i,jj) * wrk1(i,jj)
              enddo
            end if  
          else
            do i=ista,iend
              if(UP(I,J)==SPVAL) cycle            
              DPSI(I,J)   = UP(I,J)*wrk3(i,j) * wrk1(i,j)  
              if(VP(I,J)==SPVAL) cycle              
              DCHI(I,J)   = -1.0*VP(I,J)*wrk3(i,j) * wrk1(i,j)  
            enddo
          endif
          if (npass > 0) then
            do i=ista,iend
              tx1(i) = dpsi(i,j)
              tx3(i) = dchi(i,j)
            enddo
            do nn=1,npass
              do i=ista,iend
                tx2(i+1) = tx1(i)
                tx4(i+1) = tx3(i)
              enddo
              tx2(1)    = tx2(im+1)
              tx2(im+2) = tx2(2)
              tx4(1)    = tx4(im+1)
              tx4(im+2) = tx4(2)
              do i=2,im+1
                tx1(i-1) = 0.25 * (tx2(i-1) + tx2(i+1)) + 0.5*tx2(i)
                tx3(i-1) = 0.25 * (tx4(i-1) + tx4(i+1)) + 0.5*tx4(i)
              enddo
            enddo
            do i=ista,iend
              dpsi(i,j) = tx1(i)
              dchi(i,j) = tx4(i)
            enddo
          endif
        END DO                               ! end of J loop

!        DO J=JSTA,JEND
!          DO I=ISTA,IEND
!            PSI(I,J)=DPSI(I,J)
!            CHI(I,J)=DCHI(I,J)
!          ENDDO
!        ENDDO

        call fullpole(psi(ista_2l:iend_2u,jsta_2l:jend_2u),psipoles)
        call fullpole(chi(ista_2l:iend_2u,jsta_2l:jend_2u),chipoles)

!$omp  parallel do private(i,j,ii)
        DO J=JSTA,JEND
          if (j == 1) then
           if(gdlat(ista,j) > 0.) then ! count from north to south
              do i=ista,iend
                ii = i + imb2
                if (ii > im) ii = ii - im
                psi(i,J+1) = dpsi(i,j) + psipoles(ii,1)
                chi(i,J+1) = dchi(i,j) + chipoles(ii,1)
                print*,'j=1',psi(i,j-1),dpsi(i,j)
              enddo
            else ! count from south to north
              do i=ista,iend
                ii = i + imb2
                if (ii > im) ii = ii - im
                psi(i,J+1) = dpsi(i,j) - psipoles(ii,1)
                chi(i,J+1) = dchi(i,j) - chipoles(ii,1)      
              enddo
            end if      
          elseif (j == JM) then
            if(gdlat(ista,j) < 0.) then ! count from north to south
              do i=ista,iend
                ii = i + imb2
                if (ii > im) ii = ii - im
                psi(i,J-1) = dpsi(i,j) - psipoles(ii,2)
                chi(i,J-1) = dchi(i,j) - chipoles(ii,2)
                print*,'j=jm',psi(i,j-1),dpsi(i,j)
              enddo
            else ! count from south to north
              do i=ista,iend
                ii = i + imb2
                if (ii > im) ii = ii - im
                psi(i,J-1) = dpsi(i,j) + psipoles(ii,2)
                chi(i,J-1) = dchi(i,j) + chipoles(ii,2)
              enddo
            end if  
          else
            do i=ista,iend
                ii = i + imb2
                if (ii > im) ii = ii - im
              PSI(I,J) = sum(dpsi(i,1:j)) + psipoles(ii,1)
              CHI(I,J) = DCHI(I,J-1) + DCHI(I,J+1) 

              print*,'jsta,jend',jsta,jend
              print*,'j', j
              print*,'dpsi add',dpsi(i,1)+dpsi(1,2)+dpsi(1,3)
              print*,'psi j', psi(i,j)
              print*,'dpsi shape :j', j, shape(dpsi(i,1:j))
              print*,'dpsi :j', dpsi(i,1:j)
              print*,'sum', sum(dpsi(i,1:j))
              print*,'sum2', sum(dpsi(i,1:3))
            enddo
          endif
        enddo              ! end of J loop


! GFS use lon avg as one scaler value for pole point

        call exch(psi(ista_2l:iend_2u,jsta_2l:jend_2u))
        call exch(chi(ista_2l:iend_2u,jsta_2l:jend_2u))
        !call fullpole(psi(ista_2l:iend_2u,jsta_2l:jend_2u),psipoles)     
        !call fullpole(chi(ista_2l:iend_2u,jsta_2l:jend_2u),chipoles)     

        cosltemp=spval
        if(jsta== 1) cosltemp(1:im, 1)=coslpoles(1:im,1)
        if(jend==jm) cosltemp(1:im,jm)=coslpoles(1:im,2)
        psitemp=spval
        if(jsta== 1) psitemp(1:im, 1)=psipoles(1:im,1)
        if(jend==jm) psitemp(1:im,jm)=psipoles(1:im,2)
        chitemp=spval
        if(jsta== 1) chitemp(1:im, 1)=chipoles(1:im,1)
        if(jend==jm) chitemp(1:im,jm)=chipoles(1:im,2)
        
        call poleavg(IM,JM,JSTA,JEND,SMALL,cosltemp(1,jsta),SPVAL,psitemp(1,jsta))
        call poleavg(IM,JM,JSTA,JEND,SMALL,cosltemp(1,jsta),SPVAL,chitemp(1,jsta))

        if(jsta== 1) psi(ista:iend, 1)=psitemp(ista:iend, 1)
        if(jend==jm) psi(ista:iend,jm)=psitemp(ista:iend,jm)
        if(jsta== 1) chi(ista:iend, 1)=chitemp(ista:iend, 1)
        if(jend==jm) chi(ista:iend,jm)=chitemp(ista:iend,jm)
    
        deallocate (wrk1, wrk3, cosl)

!     
!     END OF ROUTINE.
!     
      RETURN
      END
