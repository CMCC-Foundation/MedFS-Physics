MODULE UV2Tgrid

IMPLICIT NONE


PRIVATE
INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(12,307)

PUBLIC InterpU
PUBLIC InterpV


CONTAINS


SUBROUTINE InterpU(rField,rFieldOut,nNumX,nNumY,nNumZd,nTime)
   implicit none
   integer nNumX,nNumY,nNumZd ,nTime
   real(wp),dimension(nNumX,nNumY,nNumZd,nTime),intent(in) :: rField
   real(wp),dimension(nNumX,nNumY,nNumZd,nTime) :: rField0
   real(wp),dimension(nNumX,nNumY,nNumZd,nTime),intent(out) :: rFieldOut
   integer i,j,k,t
   
!   character(len=8) :: fmt,x1,y1 ! format descriptor

!   fmt = '(I5.5)' ! an integer of width 5 with zeros at the left
!   write (x1,fmt) nNumX ! converting integer to string using a 'internal file'
!   write (y1,fmt) nNumY ! converting integer to string using a 'internal file'


   call mv2zero(rField,rField0,nNumX,nNumY,nNumZd,nTime)
   

!   open(unit=50,file='velocityU_input'//trim(x1)//'_'//trim(y1)//'.dat',form='FORMATTED',status='REPLACE')
!   write (50,14) ((rField(i,j,1,1),i=1,nNumX),j=1,nNumY)
!   close(50)
!   open(unit=50,file='velocityU_missingToZero'//trim(x1)//'_'//trim(y1)//'.dat',form='FORMATTED',status='REPLACE')
!   write (50,14) ((rField0(i,j,1,1),i=1,nNumX),j=1,nNumY)
!   close(50)
 
   
   
   rFieldOut(:,:,:,:)=0
   do t=1,nTime
      do k=1,nNumZd
         do j=1,nNumY
            do i=1,(nNumX-1)
               rFieldOut(i+1,j,k,t)=(rField0(i+1,j,k,t)+rField0(i,j,k,t))/2
               end do
            end do
          end do
      end do
!   open(unit=90,file='velocityU_output'//trim(x1)//'_'//trim(y1)//'.dat',form='FORMATTED',status='REPLACE')
!   write (90,14) ((rFieldout(i,j,1,1),i=1,nNumX),j=1,nNumY)
!   close(90)
 14    format (8f10.4)  
  ! print *,'InterpU after',rFieldOut(30,70,1,1)
  ! do t=1,nTime
  !    do k=1,nNumZd
  !       do j=1,nNumY
  !          rFieldOut(1,j,k,t)=0
  !          end do
  !        end do
  !    end do
   end subroutine

SUBROUTINE InterpV(rField,rFieldOut,nNumX,nNumY,nNumZd,nTime)
   implicit none
   integer nNumX,nNumY,nNumZd ,nTime
   real(wp),dimension(nNumX,nNumY,nNumZd,nTime),intent(in) :: rField
   real(wp),dimension(nNumX,nNumY,nNumZd,nTime) :: rField0
   real(wp),dimension(nNumX,nNumY,nNumZd,nTime),intent(out) :: rFieldOut
   integer i,j,k,t


   call mv2zero(rField,rField0,nNumX,nNumY,nNumZd,nTime)

   rFieldOut(:,:,:,:)=0
   do t=1,nTime
      do k=1,nNumZd
         do j=1,(nNumY-1)
            do i=1,nNumX
               rFieldOut(i,j+1,k,t)=(rField0(i,j+1,k,t)+rField0(i,j,k,t))/2
               end do
            end do
          end do
      end do
   end subroutine
   
 
SUBROUTINE mv2zero(rField,rField0,nNumX,nNumY,nNumZd,nTime)
   implicit none
   integer nNumX,nNumY,nNumZd ,nTime
   real(wp),dimension(nNumX,nNumY,nNumZd,nTime),intent(in) :: rField
   real(wp),dimension(nNumX,nNumY,nNumZd,nTime),intent(out) :: rField0

   rField0(:,:,:,:)=rField(:,:,:,:)
   WHERE( ABS( rField0(:,:,:,:) ) >= 999 ) rField0(:,:,:,:) = 0.0

   end subroutine 
 
   
END MODULE UV2Tgrid
