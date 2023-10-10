       PROGRAM FORMAT_WIND
!
!      2012/07 :: M. Adani
!
!-----------------------------------------------------
!       usage ::: FORMAT_WIND start_hour ndays
!-----------------------------------------------------

       IMPLICIT NONE

       INTEGER , PARAMETER       :: jpis = 245, jpjs = 73
       INTEGER                   :: i,j,error,n,t
       INTEGER                   :: iyr,imo,ida
       INTEGER                   :: start_hour
       INTEGER                   :: ndays
       INTEGER                   :: count_write,tot_write
       INTEGER                   :: iargc
       CHARACTER(len=2),DIMENSION (4)  :: c_hour
       DATA                      c_hour/'00','06','12','18'/

       REAL,DIMENSION(jpis,jpjs) :: mslec,clcec,u10ec,v10ec,t2mec,rhmec
       CHARACTER(len=34)         :: ecmwf_file
       CHARACTER(len=12)         :: pfile
       CHARACTER(len=4)          :: cyr
       CHARACTER(len=2)          :: cmo,cda
       CHARACTER(len=2)          :: cstart_hour
       CHARACTER(len=3)          :: cndays
       CHARACTER(len=8)          :: ymd
       CHARACTER(len=6)          :: hms

!---------------------begin program-------------------------------------

        count_write=0

!Get argument
        if(iargc().ne.2)stop 'Stop wrong number of arguments'
        call getarg(1,cstart_hour)
        call getarg(2,cndays)

        read(cstart_hour,"(i2)")start_hour
        read(cndays,"(i3)")ndays

        tot_write=ndays*4+1
 
!Open file with list of ECMWF files to be processed       
       open(20,file='fnames',form='formatted')


!       if (start_hour .ne. 0) ndays=ndays+1
       ndays=ndays+1

       do n=1,ndays

20     read(20,'(12A)',end=999) pfile

!Find year month date
       cyr=pfile(1:4)
       read(cyr(1:4),'(i4)') iyr

       cmo=pfile(5:6)
       if (cmo(1:1) .eq. '0' ) then
       read(cmo(2:2),'(i1)') imo
       else
       read(cmo,'(i2)') imo
       endif

       cda=pfile(7:8)
       if (cda(1:1) .eq. '0' ) then
       read(cda(2:2),'(i1)') ida
       else
       read(cda,'(i2)') ida
       endif

!Open Output file
       if (n .eq. 1) then
!       if (ndays .lt. 10) then
!          open(unit=50,file='ww3_wind'//pfile(1:8)//'_'//cndays(1:1)//'.ascii',       &
!               form='FORMATTED',status='new',iostat=error)
!       elseif (ndays .lt. 100) then
!           open(unit=50,file='ww3_wind'//pfile(1:8)//'_'//cndays(1:2)//'.ascii',      &
!               form='FORMATTED',status='new',iostat=error)
!       else
!           open(unit=50,file='ww3_wind'//pfile(1:8)//'_'//cndays//'.ascii',           &
!                form='FORMATTED',status='new',iostat=error)
!       endif
           open(unit=50,file='ww3_wind.ascii',           &
                form='FORMATTED',status='new',iostat=error)
       endif
          ymd=pfile(1:8)

!Open Input file
          open(40,file=pfile,form='unformatted',status='old')
          do t=1,4 !time records in a file

! Read  ECMWF data (00:00 UTC)
! --------------------------------------
             read(40) ((mslec(i,j),i=1,jpis),j=jpjs,1,-1)
             read(40) ((clcec(i,j),i=1,jpis),j=jpjs,1,-1)
             read(40) ((u10ec(i,j),i=1,jpis),j=jpjs,1,-1)
             read(40) ((v10ec(i,j),i=1,jpis),j=jpjs,1,-1)
             read(40) ((t2mec(i,j),i=1,jpis),j=jpjs,1,-1)
             read(40) ((rhmec(i,j),i=1,jpis),j=jpjs,1,-1)
 
             if ( (count_write .eq. 0) .and.  (n .eq. 1) .and. (c_hour(t) .ne. cstart_hour) ) goto 100

             hms=c_hour(t)//'0000'

             write (50,13) ymd,hms
             write (50,14) ((u10ec(i,j),i=1,jpis),j=1,jpjs)
             write (50,14) ((v10ec(i,j),i=1,jpis),j=1,jpjs)
             count_write=count_write+1
             if ( count_write .eq. tot_write ) goto 200

100  continue
             enddo  !number of time records
          close(40)
       enddo  !ndays

200  continue
       close(50)
999    close(20)

 13    format (a8,1x,a6)
 14    format (8f10.4)
 
       END PROGRAM

