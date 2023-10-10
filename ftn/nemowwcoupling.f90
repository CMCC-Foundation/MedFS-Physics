MODULE NemoWwCoupling

INTEGER, PARAMETER, PUBLIC :: NemoWw_NEMO=1
INTEGER, PARAMETER, PUBLIC :: NemoWw_WW=0
CHARACTER(len=200), PUBLIC :: NemoWwLog

LOGICAL, PRIVATE :: IsMaster = .FALSE.
INTEGER, PRIVATE :: ModelKind 
INTEGER, PRIVATE :: mpprank_master_global, mpprank_other_master_global
INTEGER, PRIVATE :: CommGlob, CommModel
CHARACTER(LEN=5), PRIVATE :: myname

CONTAINS


SUBROUTINE NemoWwFinalize()

   IMPLICIT NONE

   INCLUDE 'mpif.h'

   INTEGER :: ierr

   CALL MPI_BARRIER ( CommGlob, ierr ) 

   END SUBROUTINE NemoWwFinalize


SUBROUTINE NemoWwInit(MK,Tag,CommGlobModel)

   IMPLICIT NONE

   INCLUDE 'mpif.h'

   INTEGER,INTENT(IN) :: MK
   INTEGER,INTENT(IN) :: Tag
   INTEGER,INTENT(INOUT) :: CommGlobModel

   LOGICAL :: mpi_was_called = .FALSE.
   INTEGER :: otherMK
   INTEGER :: mppsize_glob , mppsize_model 
   INTEGER :: mpprank_glob , mpprank_model
   INTEGER :: ierr, i
   !CHARACTER(LEN=5) :: myname
   INTEGER, ALLOCATABLE :: modrank(:)

   ModelKind=MK
   CommGlob=CommGlobModel

   IF ( ModelKind == NemoWw_NEMO) THEN
      myname="NEMO"
      ELSE
      myname="WW"
      END IF
!   print *, "NemoWwInit : it is ",myname

   ! Is MPI already initialized?
   CALL MPI_INITIALIZED(mpi_was_called, ierr)
   CALL ERRCHECK(ierr, 'ERROR: unable to check for MPI!')

   ! Get the total number of processes (nemo+ww)
   CALL MPI_COMM_SIZE(CommGlob, mppsize_glob, ierr)
   CALL MPI_ERRCHECK(ierr)

   ! Get my global rank
   CALL MPI_COMM_RANK(CommGlob, mpprank_glob, ierr)
   CALL MPI_ERRCHECK(ierr)

   ! Allocate and initialize models name array
   IF (.NOT. ALLOCATED(modrank)) ALLOCATE(modrank(mppsize_glob),STAT=ierr)
!   i = ierr
!   IF (.NOT. ALLOCATED(modname)) ALLOCATE(modname(mppsize_glob),STAT=ierr)
!   ierr = ierr + i
   CALL ERRCHECK(ierr, 'ERROR: failed memory allocation!')
!   modrank(mpprank_glob+1) = mpprank_glob
!   modname(mpprank_glob+1) = TRIM(myname)

!   ! Exchange partecipating models rank
!   CALL MPI_ALLGATHER((/mpprank_model/), 1, MPI_INTEGER,     &
!       modrank, 1, MPI_INTEGER, CommGlob, ierr)
!   CALL MPI_ERRCHECK(ierr)

   ! Split the global communicator into local communicators
   ! The result is the local communicator for the processes
   ! belonging to the same model (nemo here)
   CALL MPI_COMM_SPLIT(CommGlob, ModelKind, mpprank_glob, CommModel, ierr)
   CALL MPI_ERRCHECK(ierr)

   CALL MPI_COMM_SIZE(CommModel, mppsize_model, ierr)
   CALL MPI_ERRCHECK(ierr)

   IF ( mppsize_model == mppsize_glob ) THEN

      mpprank_other_master_global=-1 ! to manage the case only one model is running (no coupling)
      mpprank_master_global=0

      ELSE

      ! Get my local rank
      CALL MPI_COMM_RANK(CommModel, mpprank_model, ierr)
      CALL MPI_ERRCHECK(ierr)

      ! Exchange partecipating models rank
      CALL MPI_ALLGATHER((/mpprank_model/), 1, MPI_INTEGER,     &
          modrank, 1, MPI_INTEGER, CommGlob, ierr)
      CALL MPI_ERRCHECK(ierr)   

      mpprank_other_master_global=-2 ! conventional value in case not master and coupling 

      IF ( mpprank_model == 0 ) THEN
         IsMaster=.TRUE.
         mpprank_master_global=mpprank_glob
         DO i=1,mppsize_glob
!         PRINT *, i, modrank(i)
            IF ( modrank(i) == 0 .and. mpprank_glob /= i-1 ) THEN
               mpprank_other_master_global=i-1
               END IF 
            END DO
         END IF

      CommGlobModel=CommModel

      IF ( IsMaster ) THEN
         WRITE(NemoWwLog,FMT="(a,i3,a,i3,a,a,a,i3,a,i3)")                        &
            "NemoWwCoupling(",Tag,                                         &
            ") : globsz ",mppsize_glob,                               &
            ", Master ",TRIM(myname),                                      &
            " globrk ",mpprank_master_global,                         &
            ", coupled globrk ",mpprank_other_master_global
         ELSE
         WRITE(NemoWwLog,FMT="(a,i3,a,i3,a,a,a,i3,a,i3)")                        &
            "NemoWwCoupling(",Tag,                                         &
            ") : globsz ",mppsize_glob,                               &
            ", NoMast ",TRIM(myname),                                      &
            " globrk ",mpprank_glob,                         &
            ", modrk ",mpprank_model
         END IF

      NemoWwLog=TRIM(NemoWwLog)
      PRINT *,NemoWwLog

      END IF

   END SUBROUTINE NemoWwInit


SUBROUTINE NemoWwHandshake(sdt_n,sdu_n,sdv_n,cds_n,sdt_w,sdu_w,sdv_w,cds_w,gdxd,gdyd,nsb,sdxd,sdyd,ldi,ldj,nsea,iaproc,naproc,mapsf)

   IMPLICIT NONE

   INCLUDE 'mpif.h'

   INTEGER :: mpprank_glob, mpprank_model
   INTEGER :: mppsize_model
   INTEGER :: ierr

   INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(12,307)

   REAL(wp), DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: sdu_n,sdv_n
   REAL(wp), DIMENSION(:,:), INTENT(IN), OPTIONAL :: sdt_n
   REAL(wp), DIMENSION(:,:), INTENT(OUT), OPTIONAL :: cds_n

   REAL(wp), DIMENSION(:), INTENT(IN), OPTIONAL :: cds_w
   REAL(wp), DIMENSION(:), INTENT(OUT), OPTIONAL :: sdt_w,sdu_w,sdv_w

   INTEGER, INTENT(IN) :: gdxd, gdyd                           !! global domain x and y -dimensions

   INTEGER, INTENT(IN), OPTIONAL :: nsb                        !! number of subdomains
   INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: sdxd, sdyd   !! subdomain x and y -dimensions of all tn, un, vn
   INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: ldi, ldj     !! local domain i and j -index of all tn, un, vn in global domain
!   INTEGER, PARAMETER :: dim=1 !! DRU da rimuovere

   INTEGER, INTENT(IN), OPTIONAL :: nsea                       !! number of sea point in global domain and in local subdomain
   INTEGER, INTENT(IN), OPTIONAL :: iaproc, naproc             !! number of actual process and number of processes involved in computation
   INTEGER, INTENT(IN), DIMENSION(:,:), OPTIONAL :: mapsf      !! map from vector to matrix

   !! ARRAY AT INTERFACE BETWEEN TWO MODELS
   REAL(wp), DIMENSION(:), ALLOCATABLE :: gdsu       !! vect of surface for global u field
   !REAL(wp), DIMENSION(:), ALLOCATABLE :: gdsd !!DRU da eliminare
   REAL(wp), DIMENSION(:), ALLOCATABLE :: lts        !! vect of surface for local field

   INTEGER :: i, j, p
   INTEGER :: nseal_c

   INTEGER :: istatus(MPI_STATUS_SIZE)
   REAL(wp) :: test(1)

   CALL MPI_COMM_RANK(CommGlob, mpprank_glob, ierr)
   CALL MPI_ERRCHECK(ierr)

   CALL MPI_COMM_RANK(CommModel, mpprank_model, ierr)
   CALL MPI_ERRCHECK(ierr)

   CALL MPI_COMM_SIZE(CommModel, mppsize_model, ierr)
   CALL MPI_ERRCHECK(ierr)

   !ALLOCATE(gdsu(gdxd*gdyd))

!!--------------------------------------
!! NOT MASTER NEMO ---------------------
!!--------------------------------------

   IF ( (.not. mpprank_other_master_global == -1) .and. (.not. isMaster) .and. ModelKind == NemoWw_NEMO ) THEN
      PRINT *, "NemoWwCoupling : NEMO not master ",mpprank_glob,mpprank_model,sdxd(mpprank_model+1),sdyd(mpprank_model+1)

      ALLOCATE(lts(sdxd(mpprank_model+1)*sdyd(mpprank_model+1)))

      CALL MPI_RECV(lts, sdxd(mpprank_model+1)*sdyd(mpprank_model+1), MPI_DOUBLE_PRECISION, 0,mpprank_model,CommModel,istatus,ierr)
      DO i=1,sdxd(mpprank_model+1)
         DO j=1,sdyd(mpprank_model+1)
            cds_n(i,j)=lts( IndexMat2Vec(i,j,sdxd(mpprank_model+1)) )
            !ierr=ierr
            END DO
         END DO

      !! 1. T - STORE IN ARRAY THE LOCAL SURFACE DATA

      CALL MatrixToVect(sdt_n(:,:),sdxd(mpprank_model+1),sdyd(mpprank_model+1),lts)

      !! 2. SEND ARRAY FROM LOCAL TO MASTER
      PRINT *,"NemoWwCoupling : send T",mpprank_model,sdxd(mpprank_model+1),sdyd(mpprank_model+1),CommModel
      CALL MPI_SEND(lts, sdxd(mpprank_model+1)*sdyd(mpprank_model+1), MPI_DOUBLE_PRECISION, 0, mpprank_model, CommModel, ierr)


      !! 3. U - STORE IN ARRAY THE LOCAL SURFACE DATA

      CALL MatrixToVect(sdu_n(:,:,1),sdxd(mpprank_model+1),sdyd(mpprank_model+1),lts)

      !! 4. SEND ARRAY FROM LOCAL TO MASTER
      PRINT *,"NemoWwCoupling : send U",mpprank_model,sdxd(mpprank_model+1),sdyd(mpprank_model+1),CommModel
      CALL MPI_SEND(lts, sdxd(mpprank_model+1)*sdyd(mpprank_model+1), MPI_DOUBLE_PRECISION, 0, mpprank_model, CommModel, ierr)


      !! 5. V - STORE IN ARRAY THE LOCAL SURFACE DATA

      CALL MatrixToVect(sdv_n(:,:,1),sdxd(mpprank_model+1),sdyd(mpprank_model+1),lts)

      !! 6. SEND ARRAY FROM LOCAL TO MASTER
      PRINT *,"NemoWwCoupling : send V",mpprank_model,sdxd(mpprank_model+1),sdyd(mpprank_model+1),CommModel
      CALL MPI_SEND(lts, sdxd(mpprank_model+1)*sdyd(mpprank_model+1), MPI_DOUBLE_PRECISION, 0, mpprank_model, CommModel, ierr)



      !! 3. END

      DEALLOCATE(lts)

      END IF

!!--------------------------------------
!! NOT MASTER WW -----------------------
!!--------------------------------------

   IF (  (.not. mpprank_other_master_global == -1) .and. (.not. isMaster) .and. ModelKind == NemoWw_WW ) THEN
      nseal_c=1+(nsea-iaproc)/naproc
      PRINT *, "NemoWwCoupling : WW_not_master ",mpprank_glob,"mpprank_model",mpprank_model,iaproc,mppsize_model,naproc,nseal_c 
      PRINT *, "NemoWwCoupling : WW_2not_master ",mpprank_glob,minval(cds_w(1:nsea)),maxval(cds_w(1:nsea))
!      PRINT *, "DRU_TEST_CDS123_nmww mpprank_model",mpprank_model," j1 j2 j3 i1 i2 i3 cds_w i1 i2 i3 ",(iaproc),   &
!                  (iaproc+naproc),(iaproc+2*naproc),cds_w(iaproc),cds_w(iaproc+naproc),cds_w(iaproc+2*naproc)
      !! 1. STORE IN ARRAY THE LOCAL SURFACE DATA
      ALLOCATE(lts(nseal_c))
      DO j=1,nseal_c
         i=iaproc+(j-1)*naproc
         lts(j)=cds_w(i)
         END DO
      PRINT *, "NemoWwCoupling : WW_3not_master ",minval(lts),maxval(lts)
!      PRINT *, "DRU_TEST_CDS123_nmww mpprank_model",mpprank_model," j1 j2 j3 lts j1 j2 j3 ",lts(1),lts(2),lts(3)
      !! 2. SEND ARRAY FROM LOCAL TO MASTER
      CALL MPI_SEND(lts, nseal_c, MPI_DOUBLE_PRECISION, 0, mpprank_model, CommModel, ierr)
      CALL MPI_ERRCHECK(ierr)
      !! 3. END
      DEALLOCATE(lts)

      CALL MPI_BCAST(sdt_w, nsea, MPI_DOUBLE_PRECISION, 0, CommModel, ierr)
      CALL MPI_ERRCHECK(ierr)

      CALL MPI_BCAST(sdu_w, nsea, MPI_DOUBLE_PRECISION, 0, CommModel, ierr)
      CALL MPI_ERRCHECK(ierr)

      CALL MPI_BCAST(sdv_w, nsea, MPI_DOUBLE_PRECISION, 0, CommModel, ierr)
      CALL MPI_ERRCHECK(ierr)

      END IF   

!!--------------------------------------
!! MASTER NEMO -------------------------
!!--------------------------------------

   IF (  (.not. mpprank_other_master_global == -1) .and. isMaster .and. ModelKind == NemoWw_NEMO ) THEN
      PRINT *, "NemoWwCoupling : NEMO master ",mpprank_glob,mpprank_model
      PRINT *, "NemoWwCoupling : Master for ",myname," has global rank ",mpprank_master_global, &
        " and remote master global rank is ",mpprank_other_master_global

!      ALLOCATE(gdsu(gdxd*gdyd))
!      IF ( .not. ALLOCATED(gdsu) ) PRINT *, "NemoWwCoupling : failed allocation gdsu"
!      CALL GatherBuild(gdxd,gdyd,sdxd,sdyd,ldi,ldj,sdu_n,gdsu)

      ALLOCATE(gdsu(gdxd*gdyd))
      IF ( .not. ALLOCATED(gdsu) ) PRINT *, "NemoWwCoupling : failed allocation gdsu"

!      IF ( .not. ALLOCATED(gdsd) ) PRINT *, "NemoWwCoupling : NEMO failed allocation gdsd",size(gdsd)
!      IF ( ALLOCATED(gdsd) ) PRINT *, "NemoWwCoupling : NEMO ok allocation gdsd",size(gdsd)
      PRINT *, "NemoWwCoupling : NEMO rec gdsu start",mpprank_other_master_global,CommGlob,size(gdsu)
      CALL MPI_RECV (gdsu, gdxd*gdyd, MPI_DOUBLE_PRECISION, mpprank_other_master_global, 3, CommGlob, istatus, ierr)
      CALL MPI_ERRCHECK(ierr)
      CALL Scatter(gdxd,gdyd,sdxd,sdyd,ldi,ldj,gdsu,cds_n)
      PRINT *, "NemoWwCoupling : NEMO rec gdsu done",mpprank_other_master_global

      PRINT *,"NemoWwCoupling : NEMO ric T"
      CALL GatherBuild(gdxd,gdyd,sdxd,sdyd,ldi,ldj,sdt_n(:,:),gdsu)
      PRINT *, "NemoWwCoupling : NEMO send A start ", mpprank_master_global
      CALL MPI_SEND (gdsu, gdxd*gdyd, MPI_DOUBLE_PRECISION, mpprank_other_master_global, 4, CommGlob, ierr)
      CALL MPI_ERRCHECK(ierr)
      PRINT *, "NemoWwCoupling : NEMO send A done ", mpprank_master_global

      PRINT *,"NemoWwCoupling : NEMO ric U"
      CALL GatherBuild(gdxd,gdyd,sdxd,sdyd,ldi,ldj,sdu_n(:,:,1),gdsu)
      PRINT *, "NemoWwCoupling : NEMO send B start ", mpprank_master_global
      CALL MPI_SEND (gdsu, gdxd*gdyd, MPI_DOUBLE_PRECISION, mpprank_other_master_global, 5, CommGlob, ierr)
      CALL MPI_ERRCHECK(ierr)
      PRINT *, "NemoWwCoupling : NEMO send B done ", mpprank_master_global

      PRINT *,"NemoWwCoupling : NEMO ric V"
      CALL GatherBuild(gdxd,gdyd,sdxd,sdyd,ldi,ldj,sdv_n(:,:,1),gdsu)
      PRINT *, "NemoWwCoupling : NEMO send C start ", mpprank_master_global
      CALL MPI_SEND (gdsu, gdxd*gdyd, MPI_DOUBLE_PRECISION, mpprank_other_master_global, 6, CommGlob, ierr)
      CALL MPI_ERRCHECK(ierr)
      PRINT *, "NemoWwCoupling : NEMO send C done ", mpprank_master_global

      DEALLOCATE(gdsu)

      END IF

!!------------------------------------
!! MASTER WW -------------------------
!!------------------------------------

   IF (  (.not. mpprank_other_master_global == -1) .and. isMaster .and. ModelKind == NemoWw_WW ) THEN
      PRINT *, "NemoWwCoupling : WW master ",mpprank_glob,mpprank_model,iaproc,mppsize_model,naproc
      PRINT *, "NemoWwCoupling : Master for ",myname," has global rank ",mpprank_master_global, &
         " and remote master global rank is ",mpprank_other_master_global

      ALLOCATE(gdsu(gdxd*gdyd))
      gdsu=0.0

      !! 1. RECEIVE ARRAY FROM ALL OTHER PROCESSES
      DO p=1,mppsize_model-1
         !iaproc=p+1
         nseal_c=1+(nsea-(p+1))/naproc
         ALLOCATE(lts(nseal_c))
         CALL MPI_RECV(lts, nseal_c, MPI_DOUBLE_PRECISION, p, p, CommModel, istatus, ierr)
         CALL MPI_ERRCHECK(ierr)
         PRINT *,"NemoWwCoupling : WW_master_ric",mpprank_glob,mpprank_model,iaproc,mppsize_model,naproc,p+1,nseal_c
         PRINT *,"NemoWwCoupling : WW_2master_ric",p+1,minval(lts),maxval(lts)
!         PRINT *, "DRU_TEST_CDS123_mww mpprank_model",mpprank_model,"p",p," j1 j2 j3 lts j1 j2 j3 ",lts(1),lts(2),lts(3)
!         PRINT *, "DRU_TEST_CDS123_mww mpprank_model",mpprank_model,"p",p," j1 j2 j3 i1 i2 i3",p+1,(p+1)+naproc,(p+1)+2*naproc
      !! 2. BUILD GLOBAL DOMAIN
         DO j=1,nseal_c
!!            i=iaproc+(j-1)*naproc  !!BUG
            i=(p+1)+(j-1)*naproc
            gdsu(IndexMat2Vec(mapsf(i,1),mapsf(i,2),gdxd))=lts(j)
            !gdsd(iaproc+(j-1)*naproc)=lts(j)
            END DO
         DEALLOCATE(lts)
         END DO
      p=0
      nseal_c=1+(nsea-(p+1))/naproc
      DO j=1,nseal_c
         i=(p+1)+(j-1)*naproc
         gdsu(IndexMat2Vec(mapsf(i,1),mapsf(i,2),gdxd))=cds_w(i)
         END DO 
      !! 3. SEND DATA TO OTHER MODEL
         PRINT *,"NemoWwCoupling : WW_3master_ric",minval(cds_w),maxval(cds_w)
      PRINT *, "NemoWwCoupling : WW_send_CDS_start",mpprank_other_master_global,CommGlob,size(gdsu),minval(gdsu),maxval(gdsu)
      CALL MPI_SEND (gdsu, gdxd*gdyd, MPI_DOUBLE_PRECISION, mpprank_other_master_global, 3, CommGlob, ierr)
      CALL MPI_ERRCHECK(ierr)
      PRINT *, "NemoWwCoupling : WW send CDS done",mpprank_other_master_global,CommGlob,size(gdsu)

      PRINT *, "NemoWwCoupling : WW rec T start"
      CALL MPI_RECV (gdsu, gdxd*gdyd, MPI_DOUBLE_PRECISION, mpprank_other_master_global, 4, CommGlob, istatus, ierr)
      CALL MPI_ERRCHECK(ierr)
      CALL MatrixToSea(gdsu,gdxd,nsea,naproc,mapsf,sdt_w)
      PRINT *, "NemoWwCoupling : WW rec T done"
      CALL MPI_BCAST(sdt_w, nsea, MPI_DOUBLE_PRECISION, 0, CommModel, ierr)
      CALL MPI_ERRCHECK(ierr)
      PRINT *, "NemoWwCoupling : WW bcast T done"

      PRINT *, "NemoWwCoupling : WW rec U start"
      CALL MPI_RECV (gdsu, gdxd*gdyd, MPI_DOUBLE_PRECISION, mpprank_other_master_global, 5, CommGlob, istatus, ierr)
      CALL MPI_ERRCHECK(ierr)
      CALL MatrixToSea(gdsu,gdxd,nsea,naproc,mapsf,sdu_w)
      PRINT *, "NemoWwCoupling : WW rec U done"
      CALL MPI_BCAST(sdu_w, nsea, MPI_DOUBLE_PRECISION, 0, CommModel, ierr)
      CALL MPI_ERRCHECK(ierr)
      PRINT *, "NemoWwCoupling : WW bcast U done"

      PRINT *, "NemoWwCoupling : WW rec V start"
      CALL MPI_RECV (gdsu, gdxd*gdyd, MPI_DOUBLE_PRECISION, mpprank_other_master_global, 6, CommGlob, istatus, ierr)
      CALL MPI_ERRCHECK(ierr)
      CALL MatrixToSea(gdsu,gdxd,nsea,naproc,mapsf,sdv_w)
      PRINT *, "NemoWwCoupling : WW rec V done"
      CALL MPI_BCAST(sdv_w, nsea, MPI_DOUBLE_PRECISION, 0, CommModel, ierr)
      CALL MPI_ERRCHECK(ierr)
      PRINT *, "NemoWwCoupling : WW bcast V done"

      DEALLOCATE(gdsu)

      END IF


   END SUBROUTINE NemoWwHandshake


SUBROUTINE MatrixToSea(gdsu,gdxd,nsea,naproc,mapsf,sdu_w)

   INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(12,307)

   REAL(wp), DIMENSION(:), INTENT(IN) :: gdsu
   INTEGER, INTENT(IN) :: gdxd
   INTEGER, INTENT(IN) :: nsea
   INTEGER, INTENT(IN) :: naproc
   INTEGER, DIMENSION(:,:), INTENT(IN) :: mapsf
   REAL(wp), DIMENSION(:), INTENT(OUT) :: sdu_w
   
   INTEGER :: i !, j  !! not required

   DO i=1,nsea
      ! j=1+(i-1)/naproc !! not required
      sdu_w(i)=gdsu(IndexMat2Vec(mapsf(i,1),mapsf(i,2),gdxd))
      END DO

   END SUBROUTINE

SUBROUTINE MatrixToVect(sdu_n,sdxd,sdyd,lts)

   INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(12,307)

   REAL(wp), DIMENSION(:,:), INTENT(IN) :: sdu_n
   INTEGER, INTENT(IN) :: sdxd, sdyd   !! subdomain x and y -dimensions of all tn, un, vn
   REAL(wp), DIMENSION(:), INTENT(OUT) :: lts

   INTEGER :: i, j

   DO i=1,sdxd
      DO j=1,sdyd
         !lts(  (j-1)*sdxd(mpprank_model+1) + i )=sdu_n(i,j,1)
         lts( IndexMat2Vec(i,j,sdxd) )=sdu_n(i,j)
         END DO
      END DO

   END SUBROUTINE


SUBROUTINE GatherBuild(gdxd,gdyd,sdxd,sdyd,ldi,ldj,sdu_n,gdsu)

   IMPLICIT NONE

   INCLUDE 'mpif.h'

   INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(12,307)

   INTEGER, INTENT(IN) :: gdxd, gdyd                 !! global domain x and y -dimensions
   INTEGER, INTENT(IN), DIMENSION(:) :: sdxd, sdyd   !! subdomain x and y -dimensions of all tn, un, vn
   INTEGER, INTENT(IN), DIMENSION(:) :: ldi, ldj     !! local domain i and j -index of all tn, un, vn in global domain
   REAL(wp), DIMENSION(:,:), INTENT(IN) :: sdu_n
   REAL(wp), INTENT(OUT), DIMENSION(:) :: gdsu

   REAL(wp), DIMENSION(:), ALLOCATABLE :: lts        !! vect of surface for local field
   INTEGER :: mppsize_model, ierr
   INTEGER :: istatus(MPI_STATUS_SIZE)
   INTEGER :: p, i, j

   CALL MPI_COMM_SIZE(CommModel, mppsize_model, ierr)
   CALL MPI_ERRCHECK(ierr)

      DO p=1,mppsize_model-1
         PRINT *,"NemoWwCoupling : NEMO ric",p,sdxd(p+1),sdyd(p+1),CommModel,ldi(p+1),ldj(p+1),gdxd,gdyd
         ALLOCATE(lts(sdxd(p+1)*sdyd(p+1)))
         CALL MPI_RECV(lts, sdxd(p+1)*sdyd(p+1), MPI_DOUBLE_PRECISION, p, p, CommModel, istatus, ierr)
         CALL MPI_ERRCHECK(ierr)
         !! 2. BUILD GLOBAL DOMAIN
         DO i=1,sdxd(p+1)
            DO j=1,sdyd(p+1)
               gdsu( IndexMat2Vec(i+ldi(p+1)-1,j+ldj(p+1)-1,gdxd)  )=lts( IndexMat2Vec(i,j,sdxd(p+1)) )
               !ierr=ierr
               END DO
            END DO
         DEALLOCATE(lts)
         END DO
      p=0
      DO i=1,sdxd(p+1)
         DO j=1,sdyd(p+1)
            gdsu( IndexMat2Vec(i+ldi(p+1)-1,j+ldj(p+1)-1,gdxd)  )=sdu_n(i,j)
            !ierr=ierr
            END DO
         END DO
   END SUBROUTINE

SUBROUTINE Scatter(gdxd,gdyd,sdxd,sdyd,ldi,ldj,gdsu,sdu_n)
   INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(12,307)
   INTEGER, INTENT(IN) :: gdxd, gdyd                 !! global domain x and y -dimensions
   INTEGER, INTENT(IN), DIMENSION(:) :: sdxd, sdyd   !! subdomain x and y -dimensions of all tn, un, vn
   INTEGER, INTENT(IN), DIMENSION(:) :: ldi, ldj     !! local domain i and j -index of all tn, un, vn in global domain
   REAL(wp), DIMENSION(:,:), INTENT(OUT) :: sdu_n
   REAL(wp), INTENT(IN), DIMENSION(:) :: gdsu

   REAL(wp), DIMENSION(:), ALLOCATABLE :: lts        !! vect of surface for local field
   INTEGER :: mppsize_model, ierr
   !INTEGER :: istatus(MPI_STATUS_SIZE)
   INTEGER :: p, i, j

   CALL MPI_COMM_SIZE(CommModel, mppsize_model, ierr)
   CALL MPI_ERRCHECK(ierr)

      DO p=1,mppsize_model-1
         PRINT *,"NemoWwCoupling : NEMO send",p,sdxd(p+1),sdyd(p+1),CommModel,ldi(p+1),ldj(p+1),gdxd,gdyd
         ALLOCATE(lts(sdxd(p+1)*sdyd(p+1)))
         !! 2. BUILD SUB DOMAIN
         DO i=1,sdxd(p+1)
            DO j=1,sdyd(p+1)
               lts( IndexMat2Vec(i,j,sdxd(p+1)) )=gdsu( IndexMat2Vec(i+ldi(p+1)-1,j+ldj(p+1)-1,gdxd)  )
               !ierr=ierr
               END DO
            END DO
         CALL MPI_SEND(lts, sdxd(p+1)*sdyd(p+1), MPI_DOUBLE_PRECISION, p, p, CommModel, ierr)
         CALL MPI_ERRCHECK(ierr)
         DEALLOCATE(lts)
         END DO
      p=0
      DO i=1,sdxd(p+1)
         DO j=1,sdyd(p+1)
            sdu_n(i,j)=gdsu( IndexMat2Vec(i+ldi(p+1)-1,j+ldj(p+1)-1,gdxd)  )
            !ierr=ierr
            END DO
         END DO

   END SUBROUTINE

FUNCTION IndexMat2Vec(x,y,nx) RESULT(i)
   IMPLICIT NONE
   INTEGER :: nx,x,y,i
   i=(y-1)*nx+x
   END FUNCTION IndexMat2Vec


SUBROUTINE ERRCHECK(ierr, str)

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: ierr
   CHARACTER(*), INTENT(IN) :: str

   IF (ierr /= 0) THEN
     WRITE(*,FMT='(A)') TRIM(str)
     STOP
   END IF

END SUBROUTINE ERRCHECK



  SUBROUTINE MPI_ERRCHECK(ierr)

   IMPLICIT NONE

   INCLUDE 'mpif.h'

   INTEGER, INTENT(IN) :: ierr

   CHARACTER(LEN=MPI_MAX_ERROR_STRING) :: str
   INTEGER :: str_len, lerr

   IF (ierr /= MPI_SUCCESS) THEN
     CALL MPI_ERROR_STRING(ierr, str, str_len, ierr)
     WRITE(*,FMT='(A)') TRIM(str)
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr, lerr)
   END IF

  END SUBROUTINE MPI_ERRCHECK


END MODULE NemoWwCoupling

