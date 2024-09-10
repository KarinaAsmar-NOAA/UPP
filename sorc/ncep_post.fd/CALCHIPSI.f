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
      REAL, dimension(IM,2) :: GLATPOLES, COSLPOLES, UPOLES, VPOLES, PSIPOLES, CHIPOLES
      REAL, dimension(IM,JSTA:JEND) :: COSLTEMP, PSITEMP, CHITEMP
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

! CALCULATE DISTANCES
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
            cosl(i,j) = cos(gdlat(i,j)*dtr)  ! cosine of latitude in radians
            IF(cosl(i,j) >= SMALL) then
              wrk1(i,j) = ERAD*cosl(i,j)  ! earth radius in meters times cosine of latitude in radians
            else
              wrk1(i,j) = 0.
            end if    
            if(i == im .or. i == 1) then
              wrk2(i,j) = (360.+GDLON(ip1,J)-GDLON(im1,J))*DTR !dlam
            else
              wrk2(i,j) = (GDLON(ip1,J)-GDLON(im1,J))*DTR      !dlam
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
                wrk3(i,j) = (180.-GDLAT(i,J+1)-GLATPOLES(ii,1))*DTR !dphi
              enddo
            else ! count from south to north
              do i=ista,iend
                ii = i + imb2
                if (ii > im) ii = ii - im
                wrk3(i,j) = (180.+GDLAT(i,J+1)+GLATPOLES(ii,1))*DTR !dphi
!
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
              wrk3(i,j) = (GDLAT(I,J-1)-GDLAT(I,J+1))*DTR !dphi
            enddo
          endif
        enddo  

! CALCULATE STRM AND VPOT
        npass = 0

        jtem = jm / 18 + 1
      
        call fullpole(UP(ista_2l:iend_2u,jsta_2l:jend_2u),upoles)
        call fullpole(VP(ista_2l:iend_2u,jsta_2l:jend_2u),vpoles)
        
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
                     UPOLES(II,1)==SPVAL .or. UP(I,J+1)==SPVAL) cycle
                  PSI(I,J) = (UPOLES(II,1)*COSLPOLES(II,1) + UP(I,J+1)*COSL(I,J+1))*wrk3(i,j)
                  CHI(I,J) = (VPOLES(II,1)*COSLPOLES(II,1) + VP(I,J+1)*COSL(I,J+1))*wrk3(i,j)
                enddo
              ELSE                                   !pole point, compute at j=2
                jj = 2
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  if(UP(ip1,JJ)==SPVAL .or. VP(im1,JJ)==SPVAL .or. &
                     UP(I,J)==SPVAL .or. UP(I,jj+1)==SPVAL) cycle
                  PSI(I,J) = (UP(I,J)*COSL(I,J)-UP(I,jj+1)*COSL(I,Jj+1))*wrk3(i,jj) 
                  CHI(I,J) = (VP(I,J)*COSL(I,J)-VP(I,jj+1)*COSL(I,Jj+1))*wrk3(i,jj) 
                enddo
              ENDIF                              ! end pole point if block
            else                              ! gdlat(ista,j) not positive
              IF(cosl(ista,j) >= SMALL) THEN            !not a pole point
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  if(VP(ip1,J)==SPVAL .or. VP(im1,J)==SPVAL .or. &
                     UPOLES(II,1)==SPVAL .or. UP(I,J+1)==SPVAL) cycle
                  PSI(I,J) = (upoles(II,1)*coslpoles(II,1)+UP(I,J+1)*COSL(I,J+1))*wrk3(i,j)
                  CHI(I,J) = (vpoles(II,1)*coslpoles(II,1)+VP(I,J+1)*COSL(I,J+1))*wrk3(i,j)
                enddo
              ELSE                                   !pole point, compute at j=2
                jj = 2
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  if(VP(ip1,JJ)==SPVAL .or. VP(im1,JJ)==SPVAL .or. &
                     UP(I,J)==SPVAL .or. UP(I,jj+1)==SPVAL) cycle
                  PSI(I,J) = (UP(I,J)*COSL(I,J)-UP(I,jj+1)*COSL(I,Jj+1))*wrk3(i,jj)
                  CHI(I,J) = (VP(I,J)*COSL(I,J)-VP(I,jj+1)*COSL(I,Jj+1))*wrk3(i,jj)
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
                     UP(I,J-1)==SPVAL .or. UPOLES(II,2)==SPVAL) cycle
                  PSI(I,J) = (UP(I,J-1)*COSL(I,J-1)+UP(II,J)*COSL(II,J))*wrk3(i,j)   &
     &                      +   upoles(II,2)*coslpoles(II,2)*wrk3(i,j) 
                  CHI(I,J) = (VP(I,J-1)*COSL(I,J-1)+VP(II,J)*COSL(II,J))*wrk3(i,j)   &
     &                      +   vpoles(II,2)*coslpoles(II,2)*wrk3(i,j)     
                enddo
              ELSE                                   !pole point,compute at jm-1
                jj = jm-1
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  if(VP(ip1,JJ)==SPVAL .or. VP(im1,JJ)==SPVAL .or. &
                     UP(I,jj-1)==SPVAL .or. UP(I,J)==SPVAL) cycle
                  PSI(I,J) = (UP(I,jj-1)*COSL(I,Jj-1)-UP(I,J)*COSL(I,J))*wrk3(i,jj)
                  CHI(I,J) = (VP(I,jj-1)*COSL(I,Jj-1)-VP(I,J)*COSL(I,J))*wrk3(i,jj)
                enddo
              ENDIF                        ! end pole point if block
            else                        ! gdlat(ista,j) not negative
              IF(cosl(ista,j) >= SMALL) THEN            !not a pole point
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  if(VP(ip1,J)==SPVAL .or. VP(im1,J)==SPVAL .or. &
                     UP(I,J-1)==SPVAL .or. UPOLES(II,2)==SPVAL) cycle
                  PSI(I,J) = (UP(I,J-1)*COSL(I,J-1)+UP(II,J)*COSL(II,J))*wrk3(i,j)   &
     &                      +   upoles(II,2)*coslpoles(II,2)*wrk3(i,j)
                  CHI(I,J) = (VP(I,J-1)*COSL(I,J-1)+VP(II,J)*COSL(II,J))*wrk3(i,j)   &
     &                      +   vpoles(II,2)*coslpoles(II,2)*wrk3(i,j)
                enddo
              ELSE                                   !pole point,compute at jm-1
                jj = jm-1
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  if(VP(ip1,JJ)==SPVAL .or. VP(im1,JJ)==SPVAL .or. &
                     UP(I,jj-1)==SPVAL .or. UP(I,J)==SPVAL) cycle
                  PSI(I,J) = (UP(I,jj-1)*COSL(I,Jj-1)-UP(I,J)*COSL(I,J))*wrk3(i,jj)
                  CHI(I,J) = (VP(I,jj-1)*COSL(I,Jj-1)-VP(I,J)*COSL(I,J))*wrk3(i,jj)
                enddo
              ENDIF
            endif
          ELSE                                    ! J NOT 1 OR JM
            DO I=ISTA,IEND
              ip1 = ie(i)
              im1 = iw(i)
              if(VP(ip1,J)==SPVAL .or. VP(im1,J)==SPVAL .or. &
                 UP(I,J-1)==SPVAL .or. UP(I,J+1)==SPVAL) cycle
              CHI(I,J) = (UP(I,J-1)*COSL(I,J-1)-UP(I,J+1)*COSL(I,J+1))*wrk3(i,j)  !! TEST PSI
              print*,'wrong psi',chi(i,j)
              PSI(I,J)   = (-1*(VP(ip1,J)+VP(im1,J))*wrk2(i,j)               &
     &                    +  (UP(I,J-1)*COSL(I,J-1)                          &
                          +   UP(I,J+1)*COSL(I,J+1))*wrk3(i,j)) * wrk1(i,j)*0.5  
 !             PSI(I,J) = 0.5*(UP(I,J-1)*COSL(I,J-1) + UP(I,J+1)*COSL(I,J+1))*wrk3(i,j)      &
 !     &           - 0.5*(VP(im1,J)+VP(ip1,J))*wrk2(i,j)
              print*,'psi',PSI(I,J)
      !        CHI(I,J) = (VP(I,J-1)*COSL(I,J-1)-VP(I,J+1)*COSL(I,J+1))*wrk3(i,j)
            ENDDO
          END IF                              ! END J IF BLOCK
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
              chi(i,j) = tx1(i)
            enddo
          endif                        ! end npass>0 if block
        END DO                               ! end of J loop

!       deallocate (wrk1, wrk2, wrk3, cosl)
! GFS use lon avg as one scaler value for pole point

      ! call poleavg(IM,JM,JSTA,JEND,SMALL,COSL(1,jsta),SPVAL,psi(1,jsta))

        call exch(psi(ista_2l:iend_2u,jsta_2l:jend_2u))
        call fullpole(psi(ista_2l:iend_2u,jsta_2l:jend_2u),psipoles)  
        call fullpole(chi(ista_2l:iend_2u,jsta_2l:jend_2u),chipoles)

        cosltemp=spval
        if(jsta== 1) cosltemp(1:im, 1)=coslpoles(1:im,1)
        if(jend==jm) cosltemp(1:im,jm)=coslpoles(1:im,2)
        psitemp=spval
        chitemp=spval
        if(jsta== 1) psitemp(1:im, 1)=psipoles(1:im,1)
        if(jend==jm) psitemp(1:im,jm)=psipoles(1:im,2)
        if(jsta== 1) chitemp(1:im, 1)=chipoles(1:im,1)
        if(jend==jm) chitemp(1:im,jm)=chipoles(1:im,2)
        
        call poleavg(IM,JM,JSTA,JEND,SMALL,cosltemp(1,jsta),SPVAL,psitemp(1,jsta))
        call poleavg(IM,JM,JSTA,JEND,SMALL,cosltemp(1,jsta),SPVAL,chitemp(1,jsta))

        if(jsta== 1) psi(ista:iend, 1)=psitemp(ista:iend, 1)
        if(jend==jm) psi(ista:iend,jm)=psitemp(ista:iend,jm)
        if(jsta== 1) chi(ista:iend, 1)=chitemp(ista:iend, 1)
        if(jend==jm) chi(ista:iend,jm)=chitemp(ista:iend,jm)
    
        deallocate (wrk1, wrk2, wrk3, cosl, iw, ie)
!     
!     END OF ROUTINE.
!     
      RETURN
      END
