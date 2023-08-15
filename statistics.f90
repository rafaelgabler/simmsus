!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: statistics				  !					         
!Last update: 16/07/2023			  !
!*************************************************!

subroutine statistics
use variables
use subroutines

write(*,*) 'Initializing statistics module...'
write(*,*) ''
! Allocating some variables
npast=npast/n2
ncor=((1+(N-1))*(N-1))/2
allocate(U(rea,N,3))
allocate(V(rea,N,6))
allocate(flut(N,npast,rea,3))
allocate(velmedia(npast,rea,3))
allocate(vmedia(npast,3))
allocate(errovmedia(npast,3))
allocate(auxcor(rea,npast,3))
allocate(auxiliarcor(rea,3))
allocate(errocor(npast,3))
allocate(funcaor(npast,3))
allocate(dif(npast,3))
allocate(aux_erro_vel(npast,rea,3))
allocate(aux_erro_var(npast,rea,6))
allocate(variancia(npast,rea,6))
allocate(var(npast,6))
allocate(errovar(npast,6))

! Defining some formats

509 FORMAT(F30.4,F30.4,F30.4)
510 FORMAT(E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x)
513 FORMAT(E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2)
514 FORMAT(E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x)
1012 FORMAT(E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2)

!********************************** AVERAGE VELOCITY ***************************************!

! The average velocity of the suspension in a given time-step is considered here as an 
! ensemble average based on all particles in all simultaneous numerical experiments.

do i=1,rea
write(rea_char, '(I3)') i
open (rea+i,file='velocity'//rea_char//'.plt',STATUS='OLD')
end do

do j=1,rea
read (rea+j,'(A)') linha1
do k=1,npast-1
read(rea+j,'(A)') linha2
do i=1,N
read(rea+j,509)U(j,i,1),U(j,i,2),U(j,i,3)
end do
velmedia(k,j,1)= sum(U(j,:,1))/N
velmedia(k,j,2)= sum(U(j,:,2))/N
velmedia(k,j,3)= sum(U(j,:,3))/N
end do
end do

do k=1,npast-1
vmedia(k,1)=sum(velmedia(k,:,1))/rea
vmedia(k,2)=sum(velmedia(k,:,2))/rea
vmedia(k,3)=sum(velmedia(k,:,3))/rea
end do

do k=1,npast-1
do j=1,rea
aux_erro_vel(k,j,1)=(velmedia(k,j,1)-vmedia(k,1))**2.0
aux_erro_vel(k,j,2)=(velmedia(k,j,2)-vmedia(k,2))**2.0
aux_erro_vel(k,j,3)=(velmedia(k,j,3)-vmedia(k,3))**2.0
end do
end do

do k=1,npast-1
errovmedia(k,1)=((1.0/(rea))*sum(aux_erro_vel(k,:,1)))**0.5
errovmedia(k,2)=((1.0/(rea))*sum(aux_erro_vel(k,:,2)))**0.5
errovmedia(k,3)=((1.0/(rea))*sum(aux_erro_vel(k,:,3)))**0.5
end do

do j=1,rea
 close(rea+j)
end do

i=2*rea+1
open (i,file='average_velocity.plt')
write(i,*) 'Variables="U","V","W","UAVERAGE","DU","DV","DW","T"'

do k=1,npast-1
write(i,510)vmedia(k,1),vmedia(k,2),vmedia(k,3),   &
(((vmedia(k,1)**2.0)+(vmedia(k,2)**2.0)+(vmedia(k,3)**2.0))**0.5),   &
errovmedia(k,1),errovmedia(k,2),errovmedia(k,3),k*dt*n2
end do


write(*,*) 'Statistical analysis over the average velocity - OK'

!****************************************** VARIANCE ***************************************!

do i=1,rea
write(rea_char, '(I3)') i
open (rea+i,file='velocity'//rea_char//'.plt',STATUS='OLD')
end do


do j=1,rea
read (rea+j,'(A)') linha1
do k=1,npast-1
read(rea+j,'(A)') linha2
do i=1,N
read(rea+j,509)U(j,i,1),U(j,i,2),U(j,i,3)

if(mistura)then

if(i.le.(N*percentual)) then
V(j,i,1) = (U(j,i,1)-velmedia(k,j,1))**2.0
V(j,i,2) = (U(j,i,1)-velmedia(k,j,1))*(U(j,i,2)-velmedia(k,j,2))
V(j,i,3) = (U(j,i,1)-velmedia(k,j,1))*(U(j,i,3)-velmedia(k,j,3))
V(j,i,4) = (U(j,i,2)-velmedia(k,j,2))**2.0
V(j,i,5) = (U(j,i,2)-velmedia(k,j,2))*(U(j,i,3)-velmedia(k,j,3))
V(j,i,6) = (U(j,i,3)-velmedia(k,j,3))**2.0
else
V(j,i,1) = 0.0
V(j,i,2) = 0.0
V(j,i,3) = 0.0
V(j,i,4) = 0.0
V(j,i,5) = 0.0
V(j,i,6) = 0.0
end if

else
V(j,i,1) = (U(j,i,1)-velmedia(k,j,1))**2.0
V(j,i,2) = (U(j,i,1)-velmedia(k,j,1))*(U(j,i,2)-velmedia(k,j,2))
V(j,i,3) = (U(j,i,1)-velmedia(k,j,1))*(U(j,i,3)-velmedia(k,j,3))
V(j,i,4) = (U(j,i,2)-velmedia(k,j,2))**2.0
V(j,i,5) = (U(j,i,2)-velmedia(k,j,2))*(U(j,i,3)-velmedia(k,j,3))
V(j,i,6) = (U(j,i,3)-velmedia(k,j,3))**2.0
end if
end do

if(mistura) then
variancia(k,j,1)=sum(V(j,:,1))/(N*percentual)
variancia(k,j,2)=sum(V(j,:,2))/(N*percentual)
variancia(k,j,3)=sum(V(j,:,3))/(N*percentual)
variancia(k,j,4)=sum(V(j,:,4))/(N*percentual)
variancia(k,j,5)=sum(V(j,:,5))/(N*percentual)
variancia(k,j,6)=sum(V(j,:,6))/(N*percentual)
else
variancia(k,j,1)=sum(V(j,:,1))/N
variancia(k,j,2)=sum(V(j,:,2))/N
variancia(k,j,3)=sum(V(j,:,3))/N
variancia(k,j,4)=sum(V(j,:,4))/N
variancia(k,j,5)=sum(V(j,:,5))/N
variancia(k,j,6)=sum(V(j,:,6))/N
end if
V=0.0
end do
end do

do k=1,npast-1
var(k,1)=sum(variancia(k,:,1))/rea
var(k,2)=sum(variancia(k,:,2))/rea
var(k,3)=sum(variancia(k,:,3))/rea
var(k,4)=sum(variancia(k,:,4))/rea
var(k,5)=sum(variancia(k,:,5))/rea
var(k,6)=sum(variancia(k,:,6))/rea
end do

do k=1,npast-1
do j=1,rea
aux_erro_var(k,j,1)=(variancia(k,j,1)-var(k,1))**2.0
aux_erro_var(k,j,2)=(variancia(k,j,2)-var(k,2))**2.0
aux_erro_var(k,j,3)=(variancia(k,j,3)-var(k,3))**2.0
aux_erro_var(k,j,4)=(variancia(k,j,4)-var(k,4))**2.0
aux_erro_var(k,j,5)=(variancia(k,j,5)-var(k,5))**2.0
aux_erro_var(k,j,6)=(variancia(k,j,6)-var(k,6))**2.0
end do
end do

do k=1,npast-1
errovar(k,1)=((1.0/(rea))*sum(aux_erro_var(k,:,1)))**0.5
errovar(k,2)=((1.0/(rea))*sum(aux_erro_var(k,:,2)))**0.5
errovar(k,3)=((1.0/(rea))*sum(aux_erro_var(k,:,3)))**0.5
errovar(k,4)=((1.0/(rea))*sum(aux_erro_var(k,:,4)))**0.5
errovar(k,5)=((1.0/(rea))*sum(aux_erro_var(k,:,5)))**0.5
errovar(k,6)=((1.0/(rea))*sum(aux_erro_var(k,:,6)))**0.5
end do

do j=1,rea
 close(rea+j)
end do

i=2*rea+2
open (i,file='variance.plt')
write(i,*) 'Variables="V11","V12","V13","V22","V23","V33","DV11","DV12","DV13","DV22","DV23","DV33","T"'

do k=2,npast-1
write(i,510)var(k,1),var(k,2),var(k,3),   &
var(k,4),var(k,5),var(k,6),   &
errovar(k,1),errovar(k,2),errovar(k,3),errovar(k,4),errovar(k,5),errovar(k,6),k*dt*n2
end do

write(*,*) 'Statistical analysis over the variance - OK'

!********************************* SELF-CORRELATION FUNCTION *******************************!

do i=1,rea
write(rea_char, '(I3)') i
open (rea+i,file='velocity'//rea_char//'.plt',STATUS='OLD')
end do

do j=1,rea
read (rea+j,'(A)') linha1
do k=1,npast-1
read(rea+j,'(A)') linha2

do i=1,N
read(rea+j,509)U(j,i,1),U(j,i,2),U(j,i,3)

if(mistura)then
if(i.le.(N*percentual))then
flut(i,k,j,1)= U(j,i,1)-velmedia(k,j,1)
flut(i,k,j,2)= U(j,i,2)-velmedia(k,j,2)
flut(i,k,j,3)= U(j,i,3)-velmedia(k,j,3)
else
flut(i,k,j,1)= 0.0
flut(i,k,j,2)= 0.0
flut(i,k,j,3)= 0.0
end if
else
flut(i,k,j,1)= U(j,i,1)-velmedia(k,j,1)
flut(i,k,j,2)= U(j,i,2)-velmedia(k,j,2)
flut(i,k,j,3)= U(j,i,3)-velmedia(k,j,3)
end if
end do
end do
end do

do j=1,rea
 close(rea+j)
end do

do j=1,rea
do i=1,npast
do k=1,npast-i
auxcor(j,k,1)=sum((flut(:,k,j,1)*flut(:,k+i-1,j,1)))/sum((flut(:,k,j,1)**2.0))
auxcor(j,k,2)=sum((flut(:,k,j,2)*flut(:,k+i-1,j,2)))/sum((flut(:,k,j,2)**2.0))
auxcor(j,k,3)=sum((flut(:,k,j,3)*flut(:,k+i-1,j,3)))/sum((flut(:,k,j,3)**2.0))
end do

funcaor(i,1)=sum(auxcor(:,:,1))/((npast-i+1))
funcaor(i,2)=sum(auxcor(:,:,2))/((npast-i+1))
funcaor(i,3)=sum(auxcor(:,:,3))/((npast-i+1))

auxiliarcor(j,1)=(funcaor(i,1)-(auxcor(j,i,1)/(npast-i+1)))**2.0
auxiliarcor(j,2)=(funcaor(i,2)-(auxcor(j,i,2)/(npast-i+1)))**2.0
auxiliarcor(j,3)=(funcaor(i,3)-(auxcor(j,i,3)/(npast-i+1)))**2.0

errocor(i,1)=((1.0/rea)*sum(auxiliarcor(:,1)))**0.5
errocor(i,2)=((1.0/rea)*sum(auxiliarcor(:,2)))**0.5
errocor(i,3)=((1.0/rea)*sum(auxiliarcor(:,3)))**0.5

auxcor=0.0
end do
end do

i=2*rea+3
open (i,file='self_correlation.plt')
write(i,*) 'Variables="R1","R2","R3","DR1","DR2","DR3","T"'

do k=1,npast-1
write(i,510)funcaor(k,1),funcaor(k,2),funcaor(k,3),errocor(k,1),errocor(k,2),errocor(k,3),k*dt*n2
end do

write(*,*) 'Statistical analysis over the self-correlation function- OK'

!********************************* SYSTEM CORRELATION TIME **********************************!

dif(1,1)=0.0
dif(1,2)=0.0
dif(1,3)=0.0

do k=2,npast-3
dif(k,1)=dif(k-1,1)+((funcaor(k,1)+funcaor(k-1,1))*(n2*dt)/2.0)
dif(k,2)=dif(k-1,2)+((funcaor(k,2)+funcaor(k-1,2))*(n2*dt)/2.0)
dif(k,3)=dif(k-1,3)+((funcaor(k,3)+funcaor(k-1,3))*(n2*dt)/2.0)
end do

i=2*rea+4
open (i,file='correlation_time.plt')
write(i,*) 'Variables="T1","T2","T3","T"'

do k=1,npast-1
write(i,510)dif(k,1),dif(k,2),dif(k,3),k*dt*n2
end do

write(*,*) 'Statistical analysis over the calculation time - OK'

!*******************************************************************************************!

write(*,*) 'Generation of output files - OK'

if(grafmag)then
 CALL SYSTEM('gnuplot -persist "script4.gnu"')
else
 CALL SYSTEM('gnuplot -persist "script5.gnu"')
end if

! Deallocating variables 

write(*,*) ''
write(*,*) 'Deallocating the variables used in the statistics module...'
write(*,*) ''

deallocate(U)
write(*,*) 'Deallocating matrix 1 - OK'
deallocate(V)
write(*,*) 'Deallocating matrix 2 - OK'
deallocate(flut)
write(*,*) 'Deallocating matrix 3 - OK'
deallocate(velmedia)
write(*,*) 'Deallocating matrix 4 - OK'
deallocate(vmedia)
write(*,*) 'Deallocating matrix 5 - OK'
deallocate(errovmedia)
write(*,*) 'Deallocating matrix 6 - OK'
deallocate(auxcor)
write(*,*) 'Deallocating matrix 7 - OK'
deallocate(funcaor)
write(*,*) 'Deallocating matrix 8 - OK'
deallocate(dif)
write(*,*) 'Deallocating matrix 9 - OK'
deallocate(aux_erro_vel)
write(*,*) 'Deallocating matrix 10- OK'
deallocate(aux_erro_var)
write(*,*) 'Deallocating matrix 11- OK'
deallocate(variancia)
write(*,*) 'Deallocating matrix 12- OK'
deallocate(var)
write(*,*) 'Deallocating matrix 13- OK'
deallocate(errovar)
write(*,*) 'Deallocating matrix 14- OK'

end subroutine statistics
