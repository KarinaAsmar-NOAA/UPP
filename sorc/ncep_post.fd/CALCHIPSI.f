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
      SUBROUTINE CALCHIPSI(UP,VP,CHI,DPSI)
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
      REAL, dimension(ista_2l:iend_2u,jsta_2l:jend_2u), intent(out) :: CHI, DPSI
      REAL, dimension(IM,2) :: GLATPOLES, COSLPOLES, UPOLES, VPOLES, PSIPOLES, CHIPOLES
      REAL, dimension(IM,JSTA:JEND) :: COSLTEMP, PSITEMP, CHITEMP
!
      real,    allocatable ::  wrk1(:,:), wrk2(:,:), wrk3(:,:), cosl(:,:)
      INTEGER, allocatable ::  IHE(:),IHW(:), IE(:),IW(:)
      REAL, dimension(ista_2l:iend_2u,jsta_2l:jend_2u) :: DCHI !, DPSI
!
      integer, parameter :: npass2=2, npass3=3
      integer I,J,ip1,im1,ii,iir,iil,jj,JMT2,imb2, npass, nn, jtem
      real    R2DX,R2DY,DVDX,DUDY,UAVG,TPH1,TPHI, tx1(im+2), tx2(im+2)
!     
!***************************************************************************
!     START CALCHIPSI HERE.
!     
!     LOOP TO COMPUTE STREAMFUNCTION AND VELOCITY POTENTIAL FROM WINDS.
   
!$omp  parallel do private(i,j)
        DO J=JSTA_2L,JEND_2U
          DO I=ISTA_2L,IEND_2U
            DPSI(I,J) = SPVAL
            CHI(I,J) = SPVAL
          ENDDO
        ENDDO

      CALL EXCH(UP)
      CALL EXCH(VP)
!
      IF (MODELNAME == 'GFS' .or. global) THEN
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
              wrk1(i,j) = ERAD*cosl(i,j)
            else
              wrk1(i,j) = 0.
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
!                    UP(II,J)==SPVAL .or. UP(I,J+1)==SPVAL) cycle
                     UPOLES(II,1)==SPVAL .or. UP(I,J+1)==SPVAL) cycle
                  DPSI(I,J) = UP(I,J+1)*wrk3(i,j) * wrk1(i,j)  
                enddo
              ELSE                                   !pole point, compute at j=2
                jj = 2
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  if(VP(ip1,JJ)==SPVAL .or. VP(im1,JJ)==SPVAL .or. &
                     UP(I,J)==SPVAL .or. UP(I,jj+1)==SPVAL) cycle
                  DPSI(I,J) = UP(I,jj+1)*wrk3(i,jj) * wrk1(i,jj) 
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
!                    UP(II,J)==SPVAL .or. UP(I,J+1)==SPVAL) cycle
                     UPOLES(II,1)==SPVAL .or. UP(I,J+1)==SPVAL) cycle
                  DPSI(I,J) = UP(I,J+1)*wrk3(i,j) * wrk1(i,j)  
                enddo
              ELSE                                   !pole point, compute at j=2
                jj = 2
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  if(VP(ip1,JJ)==SPVAL .or. VP(im1,JJ)==SPVAL .or. &
                     UP(I,J)==SPVAL .or. UP(I,jj+1)==SPVAL) cycle
                  DPSI(I,J) = UP(I,jj+1)*wrk3(i,jj) * wrk1(i,jj) 
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
                  DPSI(I,J) = upoles(II,2)*wrk3(i,j) * wrk1(i,j)   
                enddo
              ELSE                                   !pole point,compute at jm-1
                jj = jm-1
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  if(VP(ip1,JJ)==SPVAL .or. VP(im1,JJ)==SPVAL .or. &
                     UP(I,jj-1)==SPVAL .or. UP(I,J)==SPVAL) cycle
                  DPSI(I,J) = UP(I,J)*wrk3(i,jj) * wrk1(i,jj) 
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
                     UP(I,J-1)==SPVAL .or. UPOLES(II,2)==SPVAL) cycle
                  DPSI(I,J) = upoles(II,2)*wrk3(i,j) * wrk1(i,j)   
                enddo
              ELSE                                   !pole point,compute at jm-1
                jj = jm-1
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  if(VP(ip1,JJ)==SPVAL .or. VP(im1,JJ)==SPVAL .or. &
                     UP(I,jj-1)==SPVAL .or. UP(I,J)==SPVAL) cycle
                  DPSI(I,J) = UP(I,J)*wrk3(i,jj) * wrk1(i,jj) 
                enddo
              ENDIF
            endif
          ELSE
            DO I=ISTA,IEND
              ip1 = ie(i)
              im1 = iw(i)
              if(VP(ip1,J)==SPVAL .or. VP(im1,J)==SPVAL .or. &
                 UP(I,J-1)==SPVAL .or. UP(I,J+1)==SPVAL) cycle
              DPSI(I,J)   = UP(I,J+1)*wrk3(i,j) * wrk1(i,j)  
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
        call fullpole(psi(ista_2l:iend_2u,jsta_2l:jend_2u),avpoles)     

        cosltemp=spval
        if(jsta== 1) cosltemp(1:im, 1)=coslpoles(1:im,1)
        if(jend==jm) cosltemp(1:im,jm)=coslpoles(1:im,2)
        avtemp=spval
        if(jsta== 1) avtemp(1:im, 1)=avpoles(1:im,1)
        if(jend==jm) avtemp(1:im,jm)=avpoles(1:im,2)
        
        call poleavg(IM,JM,JSTA,JEND,SMALL,cosltemp(1,jsta),SPVAL,avtemp(1,jsta))

        if(jsta== 1) psi(ista:iend, 1)=avtemp(ista:iend, 1)
        if(jend==jm) psi(ista:iend,jm)=avtemp(ista:iend,jm)
    
        deallocate (wrk1, wrk2, wrk3, cosl, iw, ie)

      ELSE !(MODELNAME == 'GFS' .or. global)

      IF (GRIDTYPE == 'B')THEN
        CALL EXCH(VP)
        CALL EXCH(UP)
      ENDIF
     
      CALL DVDXDUDY(UP,VP)

      IF(GRIDTYPE == 'A')THEN
!$omp parallel do  private(i,j,jmt2,tphi,r2dx,r2dy,dvdx,dudy,uavg)
        DO J=JSTA_M,JEND_M
          JMT2 = JM/2+1
          TPHI = (J-JMT2)*(DYVAL/gdsdegr)*DTR
          DO I=ISTA_M,IEND_M
            IF(DDVDX(I,J)<SPVAL.AND.DDUDY(I,J)<SPVAL.AND. &
               UUAVG(I,J)<SPVAL.AND.UP(I,J)<SPVAL.AND.  &
     &         UP(I,J+1)<SPVAL.AND.UP(I,J-1)<SPVAL) THEN
              DVDX   = DDVDX(I,J)
              DUDY   = DDUDY(I,J)
              UAVG   = UUAVG(I,J) 
!  is there a (f+tan(phi)/erad)*u term?
              IF(MODELNAME  == 'RAPR' .OR. MODELNAME  == 'FV3R') then
                 psi(I,J) = DVDX - DUDY + F(I,J)   ! for run RAP over north pole      
              else
                 psi(I,J) = DVDX - DUDY + F(I,J) + UAVG*TAN(GDLAT(I,J)*DTR)/ERAD  ! not sure about this???
              endif
            END IF
          END DO
        END DO

      ELSE IF (GRIDTYPE == 'E')THEN
       allocate(ihw(JSTA_2L:JEND_2U), IHE(JSTA_2L:JEND_2U))
!$omp  parallel do private(j)
        DO J=JSTA_2L,JEND_2U
          IHW(J) = -MOD(J,2)
          IHE(J) = IHW(J)+1
        ENDDO
!$omp parallel do  private(i,j,jmt2,tphi,r2dx,r2dy,dvdx,dudy,uavg)
        DO J=JSTA_M,JEND_M
          JMT2 = JM/2+1
          TPHI = (J-JMT2)*(DYVAL/1000.)*DTR
          TPHI = (J-JMT2)*(DYVAL/gdsdegr)*DTR
          DO I=ISTA_M,IEND_M
            IF(VP(I+IHE(J),J) < SPVAL.AND.VP(I+IHW(J),J) < SPVAL .AND.   &
     &         UP(I,J+1) < SPVAL     .AND.UP(I,J-1) < SPVAL) THEN
              DVDX   = DDVDX(I,J)
              DUDY   = DDUDY(I,J)
              UAVG   = UUAVG(I,J)
!  is there a (f+tan(phi)/erad)*u term?
              psi(I,J) = DVDX - DUDY + F(I,J) + UAVG*TAN(TPHI)/ERAD 
            END IF
          END DO
        END DO
       deallocate(ihw, IHE)
      ELSE IF (GRIDTYPE == 'B')THEN
!        CALL EXCH(VP)      !done before dvdxdudy() Jesse 20200520
        DO J=JSTA_M,JEND_M
          JMT2 = JM/2+1
          TPHI = (J-JMT2)*(DYVAL/gdsdegr)*DTR
          DO I=ISTA_M,IEND_M         
            if(VP(I,  J)==SPVAL .or. VP(I,  J-1)==SPVAL .or. &
               VP(I-1,J)==SPVAL .or. VP(I-1,J-1)==SPVAL .or. &
               UP(I,  J)==SPVAL .or. UP(I-1,J)==SPVAL .or. &
               UP(I,J-1)==SPVAL .or. UP(I-1,J-1)==SPVAL) cycle
              DVDX   = DDVDX(I,J)
              DUDY   = DDUDY(I,J)
              UAVG   = UUAVG(I,J)
!  is there a (f+tan(phi)/erad)*u term?
           psi(I,J) = DVDX - DUDY + F(I,J) + UAVG*TAN(TPHI)/ERAD 
          END DO
        END DO 
      END IF 
      END IF
!     
!     END OF ROUTINE.
!     
      RETURN
      END
