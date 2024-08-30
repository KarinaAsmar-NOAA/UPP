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
      use vrbls2d,      only: f
      use masks,        only: gdlat, gdlon, dx, dy
      use params_mod,   only: d00, dtr, small, erad
      use ctlblk_mod,   only: jsta_2l, jend_2u, spval, modelname, global, &
                              jsta, jend, im, jm, jsta_m, jend_m, gdsdegr,&
                              ista, iend, ista_m, iend_m, ista_2l, iend_2u, me, num_procs
      use gridspec_mod, only: gridtype, dyval
      use upp_math,     only: DVDXDUDY, DDVDX, DDUDY, UUAVG

      implicit none
!
!     DECLARE VARIABLES.
!     
      REAL, dimension(ista_2l:iend_2u,jsta_2l:jend_2u), intent(in)    :: UP, VP
      REAL, dimension(ista_2l:iend_2u,jsta_2l:jend_2u), intent(out) :: CHI, PSI
      REAL, dimension(IM,2) :: GLATPOLES, COSLPOLES, UPOLES, PSIPOLES
      REAL, dimension(IM,JSTA:JEND) :: COSLTEMP, PSITEMP
!
      real,    allocatable ::  wrk1(:,:), wrk2(:,:), wrk3(:,:), cosl(:,:)
      INTEGER, allocatable ::  IHE(:),IHW(:), IE(:),IW(:)
!
      integer, parameter :: npass2=2, npass3=3
      integer I,J,ip1,im1,ii,iir,iil,jj,JMT2,imb2, npass, nn, jtem
      real    R2DX,R2DY,DVDX,DUDY,UAVG,TPH1,TPHI, tx1(im+2), tx2(im+2)
!     
!***************************************************************************
!     START CALCHIPSI HERE.
!     
!     LOOP TO COMPUTE STREAMFUNCTION AND VELOCITY POTENTIAL FROM WINDS.
!     
!$omp  parallel do private(i,j)
        DO J=JSTA_2L,JEND_2U
          DO I=ISTA_2L,IEND_2U
            PSI(I,J) = SPVAL
            CHI(I,J) = SPVAL
          ENDDO
        ENDDO

      CALL EXCH(UP)
      CALL EXCH(VP)
!
        CALL EXCH(GDLAT(ISTA_2L,JSTA_2L))
        CALL EXCH(GDLON(ISTA_2L,JSTA_2L))

        allocate (wrk1(ista:iend,jsta:jend), wrk2(ista:iend,jsta:jend),          &
     &            wrk3(ista:iend,jsta:jend), cosl(ista_2l:iend_2u,jsta_2l:jend_2u))
        allocate(iw(im),ie(im))

        imb2 = im/2
!$omp  parallel do private(i)
      do i=ista,iend
        ie(i) = i+1
        iw(i) = i-1
      enddo

!$omp  parallel do private(i,j,ip1,im1)
        DO J=JSTA,JEND
          do i=ista,iend
            ip1 = ie(i)
            im1 = iw(i)
            cosl(i,j) = cos(gdlat(i,j)*dtr)
            IF(cosl(i,j) >= SMALL) then
              wrk1(i,j) = 1.0 / (ERAD*cosl(i,j))
            else
              wrk1(i,j) = 0.
            end if    
            if(i == im .or. i == 1) then
              wrk2(i,j) = 1.0 / ((360.+GDLON(ip1,J)-GDLON(im1,J))*DTR) !1/dlam
            else
              wrk2(i,j) = 1.0 / ((GDLON(ip1,J)-GDLON(im1,J))*DTR)      !1/dlam
            end if
          enddo
        enddo
!       CALL EXCH(cosl(1,JSTA_2L))
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
          !     wrk3(i,j) = 1.0 / ((180.-GDLAT(i,J+1)-GDLAT(II,J))*DTR) !1/dphi
                wrk3(i,j) = 1.0 / ((180.-GDLAT(i,J+1)-GLATPOLES(ii,1))*DTR) !1/dphi
              enddo
            else ! count from south to north
              do i=ista,iend
                ii = i + imb2
                if (ii > im) ii = ii - im
          !     wrk3(i,j) = 1.0 / ((180.+GDLAT(i,J+1)+GDLAT(II,J))*DTR) !1/dphi
                wrk3(i,j) = 1.0 / ((180.+GDLAT(i,J+1)+GLATPOLES(ii,1))*DTR) !1/dphi
!
              enddo
            end if      
          elseif (j == JM) then
            if(gdlat(ista,j) < 0.) then ! count from north to south
              do i=ista,iend
                ii = i + imb2
                if (ii > im) ii = ii - im
          !      wrk3(i,j) = 1.0 / ((180.+GDLAT(i,J-1)+GDLAT(II,J))*DTR)
                wrk3(i,j) = 1.0 / ((180.+GDLAT(i,J-1)+GLATPOLES(ii,2))*DTR)
              enddo
            else ! count from south to north
              do i=ista,iend
                ii = i + imb2
                if (ii > im) ii = ii - im
          !     wrk3(i,j) = 1.0 / ((180.-GDLAT(i,J-1)-GDLAT(II,J))*DTR)
                wrk3(i,j) = 1.0 / ((180.-GDLAT(i,J-1)-GLATPOLES(ii,2))*DTR)
              enddo
            end if  
          else
            do i=ista,iend
              wrk3(i,j) = 1.0 / ((GDLAT(I,J-1)-GDLAT(I,J+1))*DTR) !1/dphi
            enddo
          endif
        enddo  

        npass = 0

        jtem = jm / 18 + 1
      
        call fullpole(UP(ista_2l:iend_2u,jsta_2l:jend_2u),upoles)

!$omp  parallel do private(i,j,ip1,im1,ii,jj,tx1,tx2)
        DO J=JSTA,JEND
          IF(J == 1) then                            ! Near North or South pole
            if(gdlat(ista,j) > 0.) then ! count from north to south
              IF(cosl(ista,j) >= SMALL) THEN            !not a pole point
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  if(VP(ip1,J)==SPVAL .or. VP(im1,J)==SPVAL .or. &
                     UP(II,J)==SPVAL .or. UP(I,J+1)==SPVAL) cycle
                     UPOLES(II,1)==SPVAL .or. UP(I,J+1)==SPVAL cycle
                  PSI(I,J) = ((UP(ip1,J)-UP(im1,J))*wrk3(i,j) +  (upoles(II,1)*coslpoles(II,1)))  
                  CHI(I,J) = ((UP(ip1,J)-UP(im1,J))*wrk2(i,j) +  (upoles(II,1)*coslpoles(II,1)))  
                enddo
              ELSE                                   !pole point, compute at j=2
                jj = 2
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  if(VP(ip1,JJ)==SPVAL .or. VP(im1,JJ)==SPVAL .or. &
                     UP(I,J)==SPVAL .or. UP(I,jj+1)==SPVAL) cycle
                  PSI(I,J) = ((UP(ip1,JJ)-UP(im1,JJ))*wrk3(i,jj) -  (UP(I,J)*COSL(I,J)))
                  CHI(I,J) = ((UP(ip1,J)-UP(im1,J))*wrk2(i,j) -  (UP(I,J)*COSL(I,J))))  
                enddo
              ENDIF
            else
              IF(cosl(ista,j) >= SMALL) THEN            !not a pole point
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  if(VP(ip1,J)==SPVAL .or. VP(im1,J)==SPVAL .or. &
                     UPOLES(II,1)==SPVAL .or. UP(I,J+1)==SPVAL) cycle
                  PSI(I,J) = ((UP(ip1,J)-UP(im1,J))*wrk3(i,j))              
                  CHI(I,J) = ((UP(ip1,J)-UP(im1,J))*wrk2(i,j))
                enddo
              ELSE                                   !pole point, compute at j=2
                jj = 2
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  if(VP(ip1,JJ)==SPVAL .or. VP(im1,JJ)==SPVAL .or. &
                     UP(I,J)==SPVAL .or. UP(I,jj+1)==SPVAL) cycle
                  PSI(I,J) = ((UP(ip1,JJ)-UP(im1,JJ))*wrk3(i,jj))
                  CHI(I,J) = ((UP(ip1,JJ)-UP(im1,JJ))*wrk2(i,jj))
                enddo
              ENDIF
            endif
          ELSE IF(J == JM) THEN                      ! Near North or South Pole
            if(gdlat(ista,j) < 0.) then ! count from north to south
              IF(cosl(ista,j) >= SMALL) THEN            !not a pole point
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  if(VP(ip1,J)==SPVAL .or. VP(im1,J)==SPVAL .or. &
!                    UP(I,J-1)==SPVAL .or. UP(II,J)==SPVAL) cycle
                     UP(I,J-1)==SPVAL .or. UPOLES(II,2)==SPVAL) cycle
                  PSI(I,J) = ((UP(ip1,J)-UP(im1,J))*wrk3(i,j) -  (UP(I,J-1)*COSL(I,J-1)))      
                  CHI(I,J) = ((UP(ip1,J)-UP(im1,J))*wrk2(i,j) -  (UP(I,J-1)*COSL(I,J-1)))      
                enddo
              ELSE                                   !pole point,compute at jm-1
                jj = jm-1
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  if(VP(ip1,JJ)==SPVAL .or. VP(im1,JJ)==SPVAL .or. &
                     UP(I,jj-1)==SPVAL .or. UP(I,J)==SPVAL) cycle
                  PSI(I,J) = ((UP(ip1,JJ)-UP(im1,JJ))*wrk3(i,jj) -  (UP(I,jj-1)*COSL(I,Jj-1))) 
                  PSI(I,J) = ((UP(ip1,JJ)-UP(im1,JJ))*wrk2(i,jj) -  (UP(I,jj-1)*COSL(I,Jj-1))) 
                enddo
              ENDIF
            else
              IF(cosl(ista,j) >= SMALL) THEN            !not a pole point
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  if(VP(ip1,J)==SPVAL .or. VP(im1,J)==SPVAL .or. &
!                    UP(I,J-1)==SPVAL .or. UP(II,J)==SPVAL) cycle
                     UP(I,J-1)==SPVAL .or. UP(II,2)==SPVAL) cycle
                  PSI(I,J) = ((UP(ip1,J)-UP(im1,J))*wrk3(i,j) +  (UP(I,J-1)*COSL(I,J-1)))         
                  PSI(I,J) = ((UP(ip1,J)-UP(im1,J))*wrk2(i,j) +  (UP(I,J-1)*COSL(I,J-1)))         
                enddo
              ELSE                                   !pole point,compute at jm-1
                jj = jm-1
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  if(VP(ip1,JJ)==SPVAL .or. VP(im1,JJ)==SPVAL .or. &
                     UP(I,jj-1)==SPVAL .or. UP(I,J)==SPVAL) cycle
                  PSI(I,J) = ((UP(ip1,JJ)-UP(im1,JJ))*wrk3(i,jj) +  (UP(I,jj-1)*COSL(I,Jj-1))) 
                  PSI(I,J) = ((UP(ip1,JJ)-UP(im1,JJ))*wrk2(i,jj) +  (UP(I,jj-1)*COSL(I,Jj-1))) 
                enddo
              ENDIF
            endif
          ELSE
            DO I=ISTA,IEND
              ip1 = ie(i)
              im1 = iw(i)
              if(VP(ip1,J)==SPVAL .or. VP(im1,J)==SPVAL .or. &
                 UP(I,J-1)==SPVAL .or. UP(I,J+1)==SPVAL) cycle
              PSI(I,J)   = ((UP(ip1,J)-UP(im1,J))*wrk3(i,j) -  (UP(I,J-1)*COSL(I,J-1)))
              PSI(I,J)   = ((UP(ip1,J)-UP(im1,J))*wrk2(i,j) -  (UP(I,J-1)*COSL(I,J-1)))
            ENDDO
          END IF
          if (npass > 0) then
            do i=ista,iend
              tx1(i) = psi(i,j)
            enddo
            do nn=1,npass
              do i=ista,iend
                tx2(i+1) = tx1(i)
              enddo
              tx2(1)    = tx2(im+1)
              tx2(im+2) = tx2(2)
              do i=2,im+1
                tx1(i-1) = 0.25 * (tx2(i-1) + tx2(i+1)) + 0.5*tx2(i)
              enddo
            enddo
            do i=ista,iend
              psi(i,j) = tx1(i)
            enddo
          endif
        END DO                               ! end of J loop

!       deallocate (wrk1, wrk2, wrk3, cosl)
! GFS use lon avg as one scaler value for pole point

      ! call poleavg(IM,JM,JSTA,JEND,SMALL,COSL(1,jsta),SPVAL,psi(1,jsta))

        call exch(psi(ista_2l:iend_2u,jsta_2l:jend_2u))
        call fullpole(psi(ista_2l:iend_2u,jsta_2l:jend_2u),psipoles)     

        cosltemp=spval
        if(jsta== 1) cosltemp(1:im, 1)=coslpoles(1:im,1)
        if(jend==jm) cosltemp(1:im,jm)=coslpoles(1:im,2)
        psitemp=spval
        if(jsta== 1) psitemp(1:im, 1)=psipoles(1:im,1)
        if(jend==jm) psitemp(1:im,jm)=psipoles(1:im,2)
        
        call poleavg(IM,JM,JSTA,JEND,SMALL,cosltemp(1,jsta),SPVAL,psitemp(1,jsta))

        if(jsta== 1) psi(ista:iend, 1)=psitemp(ista:iend, 1)
        if(jend==jm) psi(ista:iend,jm)=psitemp(ista:iend,jm)
    
        deallocate (wrk1, wrk2, wrk3, cosl, iw, ie)
!     
!     END OF ROUTINE.
!     
      RETURN
      END