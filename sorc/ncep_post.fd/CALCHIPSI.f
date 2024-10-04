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

      SUBROUTINE CALCHIPSI(USPL,VSPL,CHI,PSI)
!
!     INCLUDE ETA GRID DIMENSIONS.  SET/DERIVE OTHER PARAMETERS.
!     
      use ctlblk_mod, only: ISTA, IEND, JSTA, JEND, IM, JM, LSM, ME, SPVAL, MPI_COMM_COMP,&
                            num_procs, icnt, idsp, isxa, iexa, jsxa, jexa, jsta_2l, jend_2u,&
			    ista_2l,iend_2u
      use masks,        only: gdlat, gdlon
      use params_mod,   only: dtr, small, erad
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none

      include 'mpif.h'
!
!     DECLARE VARIABLES.
!     
      integer :: I, J, L, IERR
      REAL, dimension(ista_2l:iend_2u,jsta_2l:jend_2u), intent(in)    :: USPL, VSPL
      REAL, dimension(ista_2l:iend_2u,jsta_2l:jend_2u), intent(out) :: CHI, PSI
!
      integer k, m
      INTEGER, allocatable ::  IHE(:),IHW(:), IE(:),IW(:)
      integer ip1,im1,ii,iir,iil,jj,JMT2,imb2, npass, nn, jtem
      REAL, dimension(ista_2l:iend_2u,jsta_2l:jend_2u) :: DCHI, DPSI
      real, allocatable :: CHI1(:),CHISUB(:),PSI1(:),PSISUB(:),DCHI_FULL(:,:),DPSI_FULL(:,:),      &
                              CHI_OUT(:,:),PSI_OUT(:,:),cosl(:,:)
!     
!***************************************************************************
!     START CALCHIPSI HERE.
!    
!     COLLECT ALL DELTA PSI/CHI VALUES FOR ENTIRE IM,JM DOMAIN  

      ALLOCATE(DCHI_FULL(IM,JM))
      ALLOCATE(DPSI_FULL(IM,JM))
      ALLOCATE(CHI_OUT(IM,JM))
      ALLOCATE(PSI_OUT(IM,JM))

      CALL CALDELTACHISPI(USPL,VSPL,DCHI,DPSI)
      
      CALL COLLECT_ALL(DCHI(ISTA:IEND,JSTA:JEND),DCHI_FULL)
      CALL COLLECT_ALL(DPSI(ISTA:IEND,JSTA:JEND),DPSI_FULL)
       
      IF (ME==0) THEN 

        allocate(cosl(im,jm))
        allocate(iw(im),ie(im))

        imb2 = im/2
        do i=1,im
          ie(i) = i+1
          iw(i) = i-1
        enddo

        DO J=1,JM
          IF(J == 1) then                            ! Near North or South pole
            if(gdlat(ista,j) > 0.) then ! count from north to south
              IF(cosl(ista,j) >= SMALL) THEN            !not a pole point
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  PSI(I,J) = 0.0 
                  CHI(I,J) = 0.0
                enddo
              ELSE                                   !pole point, compute at j=2
                jj = 2
                DO I=1,IM
                  ip1 = ie(i)
                  im1 = iw(i)
                  PSI(ip1,jj+1) = DPSI(I,J) + PSI(im1,jj)
                  CHI(ip1,jj+1) = DCHI(I,J) + CHI(im1,jj)
                enddo
              ENDIF
            else
              IF(cosl(ista,j) >= SMALL) THEN            !not a pole point
                DO I=1,IM
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  PSI(ip1,J+1) = DPSI(I,J) 
                  CHI(ip1,J+1) = DCHI(I,J)    
                enddo
              ELSE                                   !pole point, compute at j=2
                jj = 2
                DO I=1,IM
                  ip1 = ie(i)
                  im1 = iw(i)
                  PSI(ip1,jj+1) = DPSI(I,J) + PSI(im1,jj)             
                enddo
              ENDIF
            endif
          ELSE IF(J == JM) THEN                      ! Near North or South Pole
            if(gdlat(ista,j) < 0.) then ! count from north to south
              IF(cosl(ista,j) >= SMALL) THEN            !not a pole point
                DO I=1,IM
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  PSI(ip1,2) = DPSI(I,J) + PSI(im1,J-1) 
                  CHI(ip1,2) = DCHI(I,J) + CHI(im1,J-1)                
                enddo
              ELSE                                   !pole point,compute at jm-1
                jj = jm-1
                DO I=1,IM
                  ip1 = ie(i)
                  im1 = iw(i)
                  PSI(ip1,jj) = DPSI(I,J) + PSI(im1,jj-1)
                  CHI(ip1,jj) = DCHI(I,J) + CHI(im1,jj-1)      
                enddo
              ENDIF
            else
              IF(cosl(ista,j) >= SMALL) THEN            !not a pole point
                DO I=1,IM
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  PSI(ip1,2) = DPSI(I,J) + PSI(im1,J-1)
                  CHI(ip1,2) = DCHI(I,J) + CHI(im1,J-1)             
                enddo
              ELSE                                   !pole point,compute at jm-1
                jj = jm-1
                DO I=1,IM
                  ip1 = ie(i)
                  im1 = iw(i)
                  PSI(ip1,jj) = DPSI(I,J) + PSI(im1,jj-1)
                  CHI(ip1,jj) = DPSI(I,J) + CHI(im1,jj-1)              
                enddo
              ENDIF
            endif
          ELSE
            DO I=1,IM
              ip1 = ie(i)
              im1 = iw(i)
              PSI(ip1,J-1) = DPSI(I,J) + PSI(im1,J+1)
              CHI(ip1,J+1) = DCHI(I,J) + CHI(im1,J+1)                 
            ENDDO
          END IF
	  ENDDO

      ENDIF                             ! END OF ME=0 BLOCK

      ALLOCATE(CHI1(im*jm))
      ALLOCATE(CHISUB(icnt(me)))
      ALLOCATE(PSI1(im*jm))
      ALLOCATE(PSISUB(icnt(me)))

      IF (ME==0) THEN
        k=0
        DO m=0,num_procs-1
        DO J=jsxa(m),jexa(m)
        DO I=isxa(m),iexa(m)
          k=k+1
            CHI1(k)=CHI_OUT(I,J)
            PSI1(k)=PSI_OUT(I,J)
        ENDDO
        ENDDO
	ENDDO
      ENDIF
      
      CALL MPI_SCATTERV(CHI1,icnt,idsp,MPI_REAL, &
                           CHISUB,icnt(me),MPI_REAL,0,MPI_COMM_WORLD,IERR)
      CALL MPI_SCATTERV(PSI1,icnt,idsp,MPI_REAL, &
                           PSISUB,icnt(me),MPI_REAL,0,MPI_COMM_WORLD,IERR)

 
       k=0
       DO J=JSTA,JEND
       DO I=ISTA,IEND
         k=k+1
           CHI(I,J)=CHISUB(k)
           PSI(I,J)=PSISUB(k)
       ENDDO
       ENDDO
  
      IF(ALLOCATED(CHI1)) DEALLOCATE(CHI1)
      IF(ALLOCATED(CHISUB)) DEALLOCATE(CHISUB)
      IF(ALLOCATED(PSI1)) DEALLOCATE(PSI1)
      IF(ALLOCATED(PSISUB)) DEALLOCATE(PSISUB)
      IF(ALLOCATED(CHI_OUT)) DEALLOCATE(CHI_OUT)
      IF(ALLOCATED(PSI_OUT)) DEALLOCATE(PSI_OUT)
!     
!     END OF ROUTINE.
      RETURN
      END

!-----------------------------------------------------------------------
      SUBROUTINE CALDELTACHISPI(UP,VP,DCHI,DPSI)
!
!     INCLUDE ETA GRID DIMENSIONS.  SET/DERIVE OTHER PARAMETERS.
      use masks,        only: gdlat, gdlon
      use params_mod,   only: dtr, small, erad
      use ctlblk_mod,   only: jsta_2l, jend_2u, spval, &
                              jsta, jend, im, jm, &
                              ista, iend, ista_2l, iend_2u, me
!
      implicit none
!
!     DECLARE VARIABLES.
!     
      REAL, dimension(ista_2l:iend_2u,jsta_2l:jend_2u), intent(in)    :: UP, VP
      REAL, dimension(ista_2l:iend_2u,jsta_2l:jend_2u), intent(out) :: DCHI, DPSI
      REAL, dimension(IM,2) :: GLATPOLES, COSLPOLES, UPOLES, VPOLES, DPSIPOLES, DCHIPOLES
      REAL, dimension(IM,JSTA:JEND) :: COSLTEMP, DPSITEMP, DCHITEMP
!
      real,    allocatable ::  wrk1(:,:), wrk2(:,:), wrk3(:,:), cosl(:,:)
      INTEGER, allocatable ::  IHE(:),IHW(:), IE(:),IW(:)
!
      integer, parameter :: npass2=2, npass3=3
      integer I,J,ip1,im1,ii,iir,iil,jj,JMT2,imb2, npass, nn, jtem
      real    tx1(im+2), tx2(im+2), tx3(im+2), tx4(im+2)
!     
!***************************************************************************
!     START CALDELTACHIPSI HERE.
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
            if(i == im .or. i == 1) then
              wrk2(i,j) = (360.+GDLON(ip1,J)-GDLON(im1,J))*DTR !1/dlam
            else
              wrk2(i,j) = (GDLON(ip1,J)-GDLON(im1,J))*DTR      !1/dlam
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
              wrk3(i,j) = (GDLAT(I,J-1)-GDLAT(I,J+1))*DTR !1/dphi
            enddo
          endif
        enddo  

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
                     UP(I,J)==SPVAL .or. UP(I,J+1)==SPVAL) cycle
                  DPSI(I,J) = (0.5*((UP(I,J+1)*COSL(I,J+1)+upoles(I,1)*coslpoles(I,1)))*wrk3(i,j) &
		  & - 0.5*(VP(ip1,J)+VP(im1,J))*wrk2(i,j))*wrk1(i,j)
                  DCHI(I,J) = (-0.5*((UP(I,J+1)+UP(I,J)))*wrk2(i,j) - 0.5*(VP(ip1,J)+VP(im1,J+1))*wrk3(i,j))*wrk1(i,j)
                enddo
              ELSE                                   !pole point, compute at j=2
                jj = 2
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  if(VP(im1,jj)==SPVAL .or. VP(ip1,jj)==SPVAL .or. &
                     UP(I,jj+1)==SPVAL .or. UP(I,J)==SPVAL) cycle
                  DPSI(I,J) = (0.5*((UP(I,jj+1)*COSL(I,jj+1)+UP(I,J)*COSL(I,J)))*wrk3(i,jj) &
		  & - 0.5*(VP(ip1,jj)+VP(im1,jj))*wrk2(i,jj))*wrk1(i,jj)
                  DCHI(I,J) = (-0.5*((UP(im1,jj+1)+UP(ip1,J)))*wrk2(i,jj) - 0.5*(VP(ip1,J)+VP(im1,jj+1))*wrk3(i,jj))*wrk1(i,jj)
                enddo
              ENDIF
            else
              IF(cosl(ista,j) >= SMALL) THEN            !not a pole point
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  if(VP(im1,J)==SPVAL .or. VP(ip1,J)==SPVAL .or. &
                     UP(I,1+1)==SPVAL .or. UP(I,J)==SPVAL) cycle
                  DPSI(I,J) = (0.5*((UP(I,J+1)*COSL(I,J+1)+upoles(I,1)*coslpoles(I,1)))*wrk3(i,j) &
		  & - 0.5*(VP(ip1,J)+VP(im1,J))*wrk2(i,j))*wrk1(i,j)
                  DCHI(I,J) = (-0.5*((UP(im1,J+1)+UP(ip1,J)))*wrk2(i,j) - 0.5*(VP(ip1,J)+VP(im1,J+1))*wrk3(i,j))*wrk1(i,j)              
                enddo
              ELSE                                   !pole point, compute at j=2
                jj = 2
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  if(VP(im1,jj)==SPVAL .or. VP(ip1,jj)==SPVAL .or. &
                     UP(I,jj+1)==SPVAL .or. UP(I,J)==SPVAL) cycle
                  DPSI(I,J) = (0.5*((UP(I,jj+1)*COSL(I,jj+1)+UP(I,J)*COSL(I,J)))*wrk3(i,jj) &
		  & - 0.5*(VP(ip1,jj)+VP(im1,jj))*wrk2(i,jj))*wrk1(i,jj)
                  DCHI(I,J) = (-0.5*((UP(im1,jj+1)+UP(ip1,J)))*wrk2(i,jj) - 0.5*(VP(ip1,J)+VP(im1,jj+1))*wrk3(i,jj))*wrk1(i,jj)                
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
                  if(VP(im1,J)==SPVAL .or. VP(ip1,J)==SPVAL .or. &
                     UP(I,J)==SPVAL .or. UP(I,J-1)==SPVAL) cycle
                  DPSI(I,J) = (0.5*((upoles(I,2)*coslpoles(I,2)+UP(I,J-1)*COSL(I,J-1)))*wrk3(i,j)&
		  & - 0.5*(VP(ip1,J)+VP(im1,J))*wrk2(i,j))*wrk1(i,j)
                  DCHI(I,J) = (-0.5*((UP(im1,J)+UP(ip1,J-1)))*wrk2(i,j) - 0.5*(VP(ip1,J-1)+VP(im1,J))*wrk3(i,j))*wrk1(i,j)                
                enddo
              ELSE                                   !pole point,compute at jm-1
                jj = jm-1
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  if(VP(im1,jj)==SPVAL .or. VP(ip1,jj)==SPVAL .or. &
                     UP(I,J)==SPVAL .or. UP(I,jj-1)==SPVAL) cycle
                  DPSI(I,J) = (0.5*((UP(I,J)*COSL(I,J)+UP(I,jj-1)*COSL(I,jj-1)))*wrk3(i,jj) &
		  & - 0.5*(VP(ip1,jj)+VP(im1,jj))*wrk2(i,jj))*wrk1(i,jj)
                  DCHI(I,J) = (-0.5*((UP(im1,J)+UP(ip1,jj-1)))*wrk2(i,jj) - 0.5*(VP(ip1,jj-1)+VP(im1,J))*wrk3(i,jj))*wrk1(i,jj)          
                enddo
              ENDIF
            else
              IF(cosl(ista,j) >= SMALL) THEN            !not a pole point
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  ii = i + imb2
                  if (ii > im) ii = ii - im
                  if(VP(im1,J)==SPVAL .or. VP(ip1,J)==SPVAL .or. &
                     UP(I,J)==SPVAL .or. UPOLES(I,J-1)==SPVAL) cycle
                  DPSI(I,J) = (0.5*((upoles(I,2)*coslpoles(I,2)+UP(I,J-1)*COSL(I,J-1)))*wrk3(i,j) &
		  & - 0.5*(VP(ip1,J)+VP(im1,J))*wrk2(i,j))*wrk1(i,j)
                  DCHI(I,J) = (-0.5*((UP(im1,J)+UP(ip1,J-1)))*wrk2(i,j) - 0.5*(VP(ip1,J-1)+VP(im1,2))*wrk3(i,j))*wrk1(i,j)              
                enddo
              ELSE                                   !pole point,compute at jm-1
                jj = jm-1
                DO I=ISTA,IEND
                  ip1 = ie(i)
                  im1 = iw(i)
                  if(VP(im1,jj)==SPVAL .or. VP(ip1,jj)==SPVAL .or. &
                     UP(I,J)==SPVAL .or. UP(I,jj-1)==SPVAL) cycle
                  DPSI(I,J) = (0.5*((UP(I,J)*COSL(I,J)+UP(I,jj-1)*COSL(I,jj-1)))*wrk3(i,jj) &
		  & - 0.5*(VP(ip1,jj)+VP(im1,jj))*wrk2(i,jj))*wrk1(i,jj)
                  DCHI(I,J) = (-0.5*((UP(im1,J)+UP(ip1,jj-1)))*wrk2(i,jj) - 0.5*(VP(ip1,jj-1)+VP(im1,J))*wrk3(i,jj))*wrk1(i,jj)                
                enddo
              ENDIF
            endif
          ELSE
            DO I=ISTA,IEND
              ip1 = ie(i)
              im1 = iw(i)
              if(VP(im1,J)==SPVAL .or. VP(ip1,J)==SPVAL .or. &
                 UP(I,J+1)==SPVAL .or. UP(I,J-1)==SPVAL) cycle
                  DPSI(I,J) = (0.5*(UP(I,J+1)*COSL(J,J+1)+UP(I,J-1)*COSL(I,J-1))*wrk3(i,j) &
		  & - 0.5*(VP(ip1,J)+VP(im1,J))*wrk2(i,j))*wrk1(i,j)
                  DCHI(I,J) = (-0.5*((UP(im1,J+1)+UP(ip1,J-1)))*wrk2(i,j) - 0.5*(VP(ip1,J+1)+VP(im1,J))*wrk3(i,j))*wrk1(i,j)                   
            ENDDO
          END IF
          if (npass > 0) then
!$omp  parallel do private(i,j,tx1,tx2,tx3,tx4)
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

! GFS use lon avg as one scaler value for pole point
!        call exch(dpsi(ista_2l:iend_2u,jsta_2l:jend_2u))
!        call exch(dchi(ista_2l:iend_2u,jsta_2l:jend_2u))
!        call fullpole(dpsi(ista_2l:iend_2u,jsta_2l:jend_2u),dpsipoles)     
!        call fullpole(dchi(ista_2l:iend_2u,jsta_2l:jend_2u),dchipoles)    

!        cosltemp=spval
!        if(jsta== 1) cosltemp(1:im, 1)=coslpoles(1:im,1)
!        if(jend==jm) cosltemp(1:im,jm)=coslpoles(1:im,2)
!        dpsitemp=spval
!        if(jsta== 1) dpsitemp(1:im, 1)=dpsipoles(1:im,1)
!        if(jend==jm) dpsitemp(1:im,jm)=dpsipoles(1:im,2)
!        dchitemp=spval
!        if(jsta== 1) dchitemp(1:im, 1)=dchipoles(1:im,1)
!        if(jend==jm) dchitemp(1:im,jm)=dchipoles(1:im,2)
        
!        call poleavg(IM,JM,JSTA,JEND,SMALL,cosltemp(1,jsta),SPVAL,dpsitemp(1,jsta))
!        call poleavg(IM,JM,JSTA,JEND,SMALL,cosltemp(1,jsta),SPVAL,dchitemp(1,jsta))

!        if(jsta== 1) dpsi(ista:iend, 1)=dpsitemp(ista:iend, 1)
!       if(jend==jm) dpsi(ista:iend,jm)=dpsitemp(ista:iend,jm)
!        if(jsta== 1) dchi(ista:iend, 1)=dchitemp(ista:iend, 1)
!        if(jend==jm) dchi(ista:iend,jm)=dchitemp(ista:iend,jm)
    
        deallocate (wrk1, wrk2, wrk3, cosl)
!     
!     END OF ROUTINE.
!     
      RETURN
      END
