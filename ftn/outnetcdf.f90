MODULE OutNetCDF

IMPLICIT NONE

PUBLIC WriteNc


CONTAINS



SUBROUTINE WriteNc(NX,NY,OutField,UNDEF,SUFIX,FLAGS,NOGRD,TIME,minLon,maxLon,minLat,maxLat)

   USE NETCDF
   IMPLICIT NONE
   
   INTEGER                            :: NX,NY,NOGRD,i,j,numTimes,count,delta,vartimeID,hh,ss
   INTEGER,DIMENSION(:),ALLOCATABLE   :: prevHr
   INTEGER                            :: ncid,status,lonID,latID,timeID,varID,varlonID,varlatID
   REAL,DIMENSION(NX,NY),INTENT(in)   :: OutField
   REAL,DIMENSION(NX,NY)              :: lon,lat
   INTEGER,DIMENSION(2)               :: TIME
   REAL                               :: UNDEF,minLon,maxLon,minLat,maxLat,x,y,dx,dy
   CHARACTER                          :: SUFIX*4, ENAME*4,name*4,day*8,hour*6
   CHARACTER, DIMENSION (NOGRD)       :: VarName*4
   LOGICAL,DIMENSION(NOGRD),INTENT(in):: FLAGS

 
    count=0
    DO i=1,NOGRD
        IF ( FLAGS(i) ) then       
            count=count+1 
            VarName(count)=ID2name(i)
        END IF
    END DO
    !    maxLon=minLon+(1/16)*(NX-1)
    !    maxLat=minLat+(1/16)*(NY-1)
    
    !Populate lon-lat grid
       x = minLon
       dx = (maxLon-x)/(NX-1)
       DO j=1,NY
        DO i=1,NX
            lon(i,j) = x 
            x = x + dx
         END DO
          x = minLon
       END DO
       y = minLat
       dy = (maxLat-y)/(NY-1)
       DO j=1,NY
         DO i=1,NX
          lat(i,j) = y
         END DO
         y = y + dy
       END DO 

  write (day, "(I8)") TIME(1)
  if ( TIME(2) < 100000 ) then
    write (hour, "(A1,I5)") '0',TIME(2)
  else
    write (hour, "(I6)") TIME(2)
  end if
      read(hour(1:2),"(I2)") hh
      read(hour(3:4),"(I2)") ss 
    
 status = nf90_create('out_grd.nc', nf90_noclobber+nf90_netcdf4, ncid) 
   if(status==NF90_EEXIST) then
!!      print *, "File gia' esistente: verranno aggiornati i campi dati",SUFIX(2:4)
      status=nf90_open('out_grd.nc', nf90_WRITE, ncid) 
      if (status /= nf90_noerr) call handle_err(status) 
      
!---  Get the ID of the variable  ----
                 status = nf90_inq_varid(ncid,SUFIX(2:4), varID)     
             if (status /= nf90_noerr) call handle_err(status)
             status = nf90_inq_varid(ncid,'time_counter', vartimeID)     
             if (status /= nf90_noerr) call handle_err(status)
                     
!---  Get the ID of the dimension time ----
                 status = nf90_inq_dimid(ncid,"time_counter", timeID)     
             if (status /= nf90_noerr)  call handle_err(status)
                
!---  Get the length for the unlimited dimension = number of records written so far ----
                 status = nf90_inquire_dimension(ncid,timeID,name,numTimes)     
             if (status /= nf90_noerr) call handle_err(status) 
             
             ALLOCATE(prevHr(numTimes))
               
      !!  print*,'before update: ',name,'=',numTimes
      !!  print*, minLon,maxLon,NX,minLat,maxLat,NY,UNDEF,SUFIX(1:4),' '
      !!  print*,VarName
      !!  print*, trim(SUFIX(1:4)) ,' =? ',trim(VarName(1))
       if ( trim(SUFIX(1:4))  .ne. trim(VarName(1)) ) then  
       
           delta = 0
       else 
           delta = 1
            status = nf90_get_var(ncid, vartimeID, prevHr)    
            if (status /= nf90_noerr) call handle_err(status)   
           status = nf90_put_var(ncid, vartimeID, prevHr(numTimes)+3600 , &
                              start = (/ numTimes+delta /) )! , &
                             ! count = (/  1 /)      )      
       end if 
         
       status = nf90_put_var(ncid, varID, OutField  , &
                              start = (/ 1, 1,numTimes+delta /) , &
                              count = (/ NX, NY, 1 /)      )   
                              
             
       if (status /= nf90_noerr) call handle_err(status)

!--- Close netCDF file   ----
       status = nf90_close(ncid) 
       if (status /= nf90_noerr) call handle_err(status)         
      
      DEALLOCATE(prevHr)
    else
    
             
     print*,'     NX,      minLat,                maxLat,                  NY,              UNDEF,               VarName:'
     print*, NX, minLat,maxLat,NY, UNDEF,' ',VarName(1:count) 
    status = nf90_def_dim(ncid, 'x', NX, lonID)       
    if (status /= nf90_noerr) call handle_err(status) 
    status = nf90_def_dim(ncid, 'y', NY, latID)       
    if (status /= nf90_noerr) call handle_err(status)  
!---  Create unlimited (record) dimension 'time'   ----
    status = nf90_def_dim(ncid,'time_counter',nf90_unlimited, timeID)
    if (status /= nf90_noerr) call handle_err(status) 
    
    status = nf90_def_var(ncid,'nav_lon',NF90_FLOAT, &         
                                    (/lonID,latID/),varlonID)
    status = nf90_def_var(ncid,'nav_lat',NF90_FLOAT, &         
                                    (/lonID,latID/),varlatID)
    status = nf90_def_var(ncid,'time_counter',NF90_FLOAT, &         
                                    (/timeID/),vartimeID)        
    
    DO i=1,count
    Name=VarName(i)
    status = nf90_def_var(ncid,trim((Name(2:4))),NF90_FLOAT, &           
                                    (/lonID,latID,timeID/),varID)  
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_def_var_deflate(ncid,varID,0,1,1)
    if (status /= nf90_noerr) call handle_err(status)
    !status =nf90_put_att(ncid,varID,'units',trim(UNITS))
    !if (status /= nf90_noerr) call handle_err(status)  
    status =nf90_put_att(ncid,varID,'missing_value',UNDEF)
    if (status /= nf90_noerr) call handle_err(status)
    
    
    END DO
    
    status =nf90_put_att(ncid,varlonID,'valid_min',minLon)
    if (status /= nf90_noerr) call handle_err(status)
    status =nf90_put_att(ncid,varlonID,'valid_max',maxLon)
    if (status /= nf90_noerr) call handle_err(status)
    status =nf90_put_att(ncid,varlatID,'valid_min',minLat)
    if (status /= nf90_noerr) call handle_err(status)
    status =nf90_put_att(ncid,varlatID,'valid_max',maxLat)
    if (status /= nf90_noerr) call handle_err(status)
    status =nf90_put_att(ncid,vartimeID,'units', &
    'seconds since '//day(1:4)//'-'//day(5:6)//'-'//day(7:8)//' 00:00:00' )
                                                      !//hour(1:2)//':'//hour(3:4)//':00')
     
!---  End definitions -> leave define mode  ----
    status = nf90_enddef(ncid)  
    if (status /= nf90_noerr) call handle_err(status)
    
    
 !---  Get the ID of the variable  ----
     status = nf90_inq_varid(ncid,SUFIX(2:4), varID)     
     if (status /= nf90_noerr) call handle_err(status)   
!---  Write variable in netCDF file  ----
    status = nf90_put_var(ncid, varlonID, lon) 
    if (status /= nf90_noerr) call handle_err(status) 
    status = nf90_put_var(ncid, varlatID, lat) 
    if (status /= nf90_noerr) call handle_err(status)
    
    status = nf90_put_var(ncid, vartimeID, hh*3600+ss) 
    if (status /= nf90_noerr) call handle_err(status)
    
    status = nf90_put_var(ncid, varID, OutField) 
    if (status /= nf90_noerr) call handle_err(status) 
!---  Close netCDF file  ----
    status = nf90_close(ncid) 
    if (status /= nf90_noerr) call handle_err(status)
    
    
    endif
   

  end SUBROUTINE
  

SUBROUTINE handle_err(status)
    use netcdf
    implicit none 
    integer, intent (in) :: status    
    character (len = 80) :: nf_90_strerror
    if (status /= nf90_noerr) then
       write (*,*) nf90_strerror(status) 
       stop 'Stopped'
      end if
    end subroutine handle_err
    
    
    
    
FUNCTION ID2name (ID) RESULT (ENAME) 
    implicit none 
    integer, intent (in)        :: ID 
    CHARACTER                   :: ENAME*4
 
            IF ( ID .EQ.  1 ) THEN
       
                ENAME  = '.dpt'
              ELSE IF ( ID .EQ.  2 ) THEN
                ENAME  = '.cur'
              ELSE IF ( ID .EQ.  3 ) THEN
                ENAME  = '.wnd'
              ELSE IF ( ID .EQ.  4 ) THEN
                ENAME  = '.dt'
              ELSE IF ( ID .EQ.  5 ) THEN
                ENAME  = '.ust'
              ELSE IF ( ID .EQ.  6 ) THEN
                ENAME  = '.hs'
              ELSE IF ( ID .EQ.  7 ) THEN
                ENAME  = '.l'
              ELSE IF ( ID .EQ.  8 ) THEN
                ENAME  = '.t'
              ELSE IF ( ID .EQ.  9 ) THEN
                ENAME  = '.dir'
              ELSE IF ( ID .EQ. 10 ) THEN
                ENAME  = '.spr'
              ELSE IF ( ID .EQ. 11 ) THEN
                ENAME  = '.fp'
              ELSE IF ( ID .EQ. 12 ) THEN
                ENAME  = '.dp'
              ELSE IF ( ID .EQ. 13 ) THEN
                ENAME  = '.fpl'
              ELSE IF ( ID .EQ. 14 ) THEN
                ENAME  = '.dpl'
              ELSE IF ( ID .EQ. 15 ) THEN
                ENAME  = '.phs'
              ELSE IF ( ID .EQ. 16 ) THEN
                ENAME  = '.ptp'
              ELSE IF ( ID .EQ. 17 ) THEN
                ENAME  = '.plp'
              ELSE IF ( ID .EQ. 18 ) THEN
                ENAME  = '.pth'
              ELSE IF ( ID .EQ. 19 ) THEN
                ENAME  = '.psi'
              ELSE IF ( ID .EQ. 20 ) THEN
                ENAME  = '.pws'
              ELSE IF ( ID .EQ. 21 ) THEN
                ENAME  = '.wsf'
              ELSE IF ( ID .EQ. 22 ) THEN
                ENAME  = '.pnr'
              ELSE IF ( ID .EQ. 23 ) THEN
                ENAME  = '.dtd'
              ELSE IF ( ID .EQ. 24 ) THEN
                ENAME  = '.fc'
              ELSE IF ( ID .EQ. 25 ) THEN
                ENAME  = '.ice'
              ELSE IF ( ID .EQ. 26 ) THEN
                ENAME  = '.wlv'
              ELSE IF ( ID .EQ. 27 ) THEN
                ENAME  = '.abr'
              ELSE IF ( ID .EQ. 28 ) THEN
                ENAME  = '.ubr'
              ELSE IF ( ID .EQ. 29 ) THEN
                ENAME  = '.Sxy'
              ELSE IF ( ID .EQ. 30 ) THEN
                ENAME  = '.Cdg'
              ELSE IF ( ID .EQ. 31 ) THEN
                ENAME  = '.SDx'
              ELSE IF ( ID .EQ. 32 ) THEN
                ENAME  = '.SDy'
              ELSE IF ( ID .EQ. 33 ) THEN
                ENAME  = '.nws'
          ELSE
                ENAME = 'MISS'  
       END IF
 

 END FUNCTION  

END MODULE OutNetCDF
