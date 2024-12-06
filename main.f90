!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: main				  !					         
!Last update: 16/07/2023			  !
!*************************************************!

subroutine main

use variables
use subroutines

! Pre-calibrated numerical variables

pi = acos(-1.0)
nb=125 ! number of phiscal cells for periodic interactions
nbr=27 ! number of reciprocal cells for periodic interactions
Str=1.0E-01 ! rotational Stokes number
shearratei=shearrate ! shear-rate
Per=(4.0/3.0)*Pe ! rotational Peclet number

! "estatica" = TRUE: Monte carlo simulation (number of time-steps = 2)
! "estatica" = FALSE: Dynamic simulation (number of time-steps = time/dt)

if(estatica) then
npast=2
else
npast=tempo/dt
end if

! Number of random numbers used in each time-step

nnr=3*N*rea

! Allocating variables in the memory

 call allocatevariables
 
! Defining number formats for the output files

1012 FORMAT(E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2)
2024 FORMAT(E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2)

! Cleaning some important variables

X=0.0
U=0.0
aux1=0.0
aux2=0.0
aux3=0.0
Di=0.0
nr=0.0
FORCAS=0.0
FT=0.0
hidrodinamica_aux1=0.0
hidrodinamica_aux2=0.0
hidro1=0.0
hidro2=0.0
contribuicao_self=0.0
contribuicao_fisico=0.0
contribuicao_reciproco=0.0
torquereal=0.0
torquereciproco=0.0
auxt=0.0

! Defining the size of the particles

 call particle_distribution

! Calculating the size of the simulation box based on the 
! number of particles and on the volume fraction defined
! by the user in the simconfig.dat

 call box_size

! Creating the initial particle distribution

 call condicao_inicial

! Creating the files to write the results

if(.not.continua) then
 call gera_arquivos(posicao,velocidade,rea)
end if

! Defining local simulation parameters

qsi=1.0*((pi**0.5)/((l*l*h)**(1.0/3.0))) ! convergence parameter for the periodic sums

! Building a table with all the Green-functions
! And building the periodic structure to compute
! the Ewald summations

if(periodicidade) then
 call tabelagreen(qsi,l,nb,nbr,h)
 call estrutura_periodica
end if

! Creating the initial distribution of the particles dipole moments

if(magpart) then
if(.not.continua) then
 call distribui_dipolo(Di,rea,N)
end if
end if


print *,'******************************************************************************'
print *,'*                                                                            *'
print *,'*                   INITIAL CONDITIONS SUCCESSFULLY GENERATED                *'
print *,'*                                                                            *'
print *,'******************************************************************************'
print *,''

 call field_excitations

! Checking if we are continuing an old simulation or starting a new one

if(continua)then
iter=iter
auxiliar_continua=npast
else
iter=1
auxiliar_continua=npast-1
end if
aux_real=auxiliar_continua

! Printing initial fields of the local volume fraction

if(printphi)then
call campo_phi(rea,k)   
end if

! Beggining of the main loop

print *,'******************************************************************************'
print *,'*                                                                            *'
print *,'*                                SIMULATING                                  *'
print *,'*                                                                            *'
print *,'******************************************************************************'
print *,''

do k=iter, auxiliar_continua

k_real=k

! Updating the instant shear-rate 

if(shear) then
if(oscillatory) then
if(bifshear) then
shearrate=gpvetor(k)*sin(freq*dt*k)
else
shearrate=shearratei*sin(freq*dt*k)
end if
end if
end if

! Calculating Brownian forces and torques 

if(browniano)then
 call brownian
end if

! Calculating gravitational forces

if(gravidade)then
 call gravity
end if

! Calculating repulsion forces between all particles and
! contact forces between overlapped particles

 call repulsion

! Calculating non-periodic long-range dipolar forces between
! the particles

if(magpart)then
if(.not.fmagper) then
 call forca_magnetica
end if
end if

! Calculating magnetic forces due to an external magnetic field
! At this point we have been considering that magnetic field
! gradients are null within the suspension space due to the
! scale of the simulation box, which is a continuum volume point.
! This way the magnetic force due to an external magnetic field
! should be null, but if the user wants to simulate the effects
! of a non-null magnetic field gradient within the suspension
! space, then the commented command "call campo_externo" should
! be enabled.

!if(externo) then
! call campo_externo
!else
FORCAS(5,:,:,:)=0.0
!end if

! Calculating the total force acting on each particle

 !$OMP PARALLEL DO
do j=1,rea
do i=1,N
FT(j,i,1)=FORCAS(1,j,i,1)+FORCAS(2,j,i,1)+FORCAS(3,j,i,1)+FORCAS(4,j,i,1)+FORCAS(5,j,i,1)+FORCAS(6,j,i,1)
FT(j,i,2)=FORCAS(1,j,i,2)+FORCAS(2,j,i,2)+FORCAS(3,j,i,2)+FORCAS(4,j,i,2)+FORCAS(5,j,i,2)+FORCAS(6,j,i,2)
FT(j,i,3)=FORCAS(1,j,i,3)+FORCAS(2,j,i,3)+FORCAS(3,j,i,3)+FORCAS(4,j,i,3)+FORCAS(5,j,i,3)+FORCAS(6,j,i,3)
end do
end do
 !$OMP END PARALLEL DO

! Calculating periodic interactions

if(periodicidade) then
call periodic_interactions
end if

! Calculating the particles velocities
if(.not.ligaih) then
if(inertia) then
 !$OMP PARALLEL DO
do j=1,rea
do i=1,N
 call resvel(U(j,i,1),dt,St,FT(j,i,1))
 call resvel(U(j,i,2),dt,St,FT(j,i,2))
 call resvel(U(j,i,3),dt,St,FT(j,i,3))
end do
end do
 !$OMP END PARALLEL DO
else
U=FT
endif
end if

! Adding the shear contribution to the velocity of each particle

if(shear)then
 !$OMP PARALLEL DO
do j=1,rea
do i=1,N
U(j,i,2) = U(j,i,2) + shearrate*X(j,i,3)
end do
end do
 !$OMP END PARALLEL DO
end if


!***************************** FLUCTUATION MODE OPTION ***********************************!

! Here the velocity of the particles is subtracted from the average velocity of the system
! and the user will observe the behavior of the particles moving with a coordinate system
! at the average speed of the suspension and will only see flucuations induced by long-range
! interactions.

if(leito)then
do i=1,rea
usistema(i,1)=sum(U(i,:,1))/N
usistema(i,2)=sum(U(i,:,2))/N
usistema(i,3)=sum(U(i,:,3))/N
end do

 !$OMP PARALLEL DO
do q=1,rea
do i=1,N
U(q,i,1)=U(q,i,1)-usistema(q,1)
U(q,i,2)=U(q,i,2)-usistema(q,2)
U(q,i,3)=U(q,i,3)-usistema(q,3)
end do
end do
 !$OMP END PARALLEL DO
end if
!*****************************************************************************************!

! Calculating the current position of the particles using a straighfoward Euler integration
 !$OMP PARALLEL DO
do j=1,rea
 do i=1,N
 call respos(X(j,i,1),dt,U(j,i,1))
 call respos(X(j,i,2),dt,U(j,i,2))
 call respos(X(j,i,3),dt,U(j,i,3))
 end do
end do
 !$OMP END PARALLEL DO

!*****************************************************************************************!
! Imposing periodic boundary conditions to avoid any dependence of the system with respect
! to physical walls

 !$OMP PARALLEL DO
do j=1,rea
do i=1,N
if(X(j,i,1).gt.l) then
X(j,i,1)=X(j,i,1)-l
end if
if(X(j,i,1).lt.0.0)then
X(j,i,1)=l-X(j,i,1)
end if
if(X(j,i,2).gt.l) then
X(j,i,2)=X(j,i,2)-l
end if
if(X(j,i,2).lt.0.0)then
X(j,i,2)=l-X(j,i,2)
end if
if(X(j,i,3).gt.h) then
X(j,i,3)=X(j,i,3)-h
end if
if(X(j,i,3).lt.0.0)then
X(j,i,3)=h-X(j,i,3)
end if
end do
end do
 !$OMP END PARALLEL DO


!********************************************************************************************************!

! Writting in a data file the positions, velocities and orientations of the particles in the current time
! step

 call writting_files(k,k_real)
  

! SOLUTION OF THE ROTATIONAL MOTION OF THE PARTICLES

if(torque) then

! Calculating long-range magnetic torques on each particle due to particle interaction

 if(magpart)then
 if(.not.tmagper)then
 call torque_magnetico
 end if
 
! Calculating the magnetic torques induced by an external field 

 if(externo) then
 if(rotating)then
 call rotating_field(alpha, freqcampo*k*dt)
 else
 if(oscilacampo) then
 call torque_externo(alpha*campo(k))
 else
 call torque_externo(alpha)
 end if
 end if
 end if
 end if

! Brownian torques have already been computed in subroutine Browian

! Calculating the total torques acting on the particles

do j=1,rea
do i=1,N
if(browniano) then
Tt(j,i,1)= TORQUES(1,j,i,1) + TORQUES(2,j,i,1) + TORQUES(3,j,i,1)
Tt(j,i,2)= TORQUES(1,j,i,2) + TORQUES(2,j,i,2) + TORQUES(3,j,i,2)
Tt(j,i,3)= TORQUES(1,j,i,3) + TORQUES(2,j,i,3) + TORQUES(3,j,i,3)
else
Tt(j,i,1)= TORQUES(1,j,i,1) + TORQUES(2,j,i,1) 
Tt(j,i,2)= TORQUES(1,j,i,2) + TORQUES(2,j,i,2) 
Tt(j,i,3)= TORQUES(1,j,i,3) + TORQUES(2,j,i,3) 
end if
end do
end do


! Solving the angular velocity of the particles

if(mistura)then
do j=1,rea
do i=(percentual*N)+1,N
 call resomega(W(j,i,1),dt,Str,Tt(j,i,1))
 call resomega(W(j,i,2),dt,Str,Tt(j,i,2))
 call resomega(W(j,i,3),dt,Str,Tt(j,i,3))
end do
end do
else
do j=1,rea
do i=1,N
 call resomega(W(j,i,1),dt,Str,Tt(j,i,1))
 call resomega(W(j,i,2),dt,Str,Tt(j,i,2))
 call resomega(W(j,i,3),dt,Str,Tt(j,i,3))
end do
end do
end if

if(shear)then
do j=1,rea
do i=1,N
W(j,i,1) = W(j,i,1) - shearrate*0.5
end do
end do
end if

! Evolving the dipole moment of the particles with their angular velocities

if(mistura)then
do j=1,rea
do i=(N*percentual)+1,N
 call evoldip(Di(j,i,1),Di(j,i,2),Di(j,i,3),W(j,i,2),W(j,i,3),dt)
 call evoldip(Di(j,i,2),Di(j,i,3),Di(j,i,1),W(j,i,3),W(j,i,1),dt)
 call evoldip(Di(j,i,3),Di(j,i,1),Di(j,i,2),W(j,i,1),W(j,i,2),dt)
end do
end do
else
do j=1,rea
do i=1,N
 call evoldip(Di(j,i,1),Di(j,i,2),Di(j,i,3),W(j,i,2),W(j,i,3),dt)
 call evoldip(Di(j,i,2),Di(j,i,3),Di(j,i,1),W(j,i,3),W(j,i,1),dt)
 call evoldip(Di(j,i,3),Di(j,i,1),Di(j,i,2),W(j,i,1),W(j,i,2),dt)
end do
end do
end if

! Normalizing the dipole moments

do j=1,rea
if(mistura)then
do i=1,(percentual*N)
Di(j,i,1)=0.0
Di(j,i,2)=0.0
Di(j,i,3)=0.0
end do
do i=(percentual*N)+1,N
modulodipolo=(Di(j,i,1)**2.0 + Di(j,i,2)**2.0 + Di(j,i,3)**2.0)**0.5
Di(j,i,1)=Di(j,i,1)/modulodipolo
Di(j,i,2)=Di(j,i,2)/modulodipolo
Di(j,i,3)=Di(j,i,3)/modulodipolo
end do
else
do i=1,N
modulodipolo=(Di(j,i,1)**2.0 + Di(j,i,2)**2.0 + Di(j,i,3)**2.0)**0.5
Di(j,i,1)=Di(j,i,1)/modulodipolo
Di(j,i,2)=Di(j,i,2)/modulodipolo
Di(j,i,3)=Di(j,i,3)/modulodipolo
end do
end if
end do

! Calculating the magnetization of the suspension

if(grafmag) then
 call media_dipolo(Di,N,rea,magtempo(1,k),1)
 call media_dipolo(Di,N,rea,magtempo(2,k),2)
 call media_dipolo(Di,N,rea,magtempo(3,k),3)

! Calculating the magnetization errorbar

do j=1,rea
do i=1,N
flutmag(i,j)=(Di(j,i,3)-magtempo(3,k))**2.0
end do
end do
erromag=((1.0/(N*rea))*sum(flutmag))**0.5

! Calculating the magnetization derivatives in each direction

derivada1=(magtempo(1,k)-magtempo(1,k-1))/dt
derivada2=(magtempo(2,k)-magtempo(2,k-1))/dt
derivada3=(magtempo(3,k)-magtempo(3,k-1))/dt

! Writting the current magnetization components and their derivatives in a data file

write(5*rea,2024) campo(k),y(k), magtempo(1,k),magtempo(2,k),magtempo(3,k), derivada1, derivada2, derivada3, k*dt

! Writting aditional magnetization datafiles related to dynamical increase of the field's frequency

contfreqinteiro1= ((k-1)*dt)/intervalo
contfreqinteiro2= (k*dt)/intervalo

if(contfreqinteiro1.ne.contfreqinteiro2) then
multiplofreq=multiplofreq+1
end if

if(bifurcation)then 
write(400*rea+multiplofreq+1,2024) campo(k),y(k), magtempo(1,k),magtempo(2,k),magtempo(3,k), derivada1, derivada2, derivada3, k*dt
end if
end if

tempototal(k)=k*dt
end if
end do

! Printing local volume fraction maps inside the simulation box
if(printphi) then
call campo_phi(rea,k)
end if

! Calculating the final structure factor of the suspension
if(fator)then
 call fator_estrutura(X,N,l,h,dt,rea)
end if

! Closing files
do j=1,2*rea
 close(j)
end do

 close(100*rea)
 close(300*rea)


! Deallocating all matrices and vectors

deallocate(X, STAT = DeAllocateStatus)
deallocate(U, STAT = DeAllocateStatus)
deallocate(FORCAS, STAT = DeAllocateStatus)
deallocate(FT, STAT = DeAllocateStatus)
deallocate(nr, STAT = DeAllocateStatus)
deallocate(hidrodinamica_aux1, STAT = DeAllocateStatus)
deallocate(hidrodinamica_aux2, STAT = DeAllocateStatus)
deallocate(hidro1, STAT = DeAllocateStatus)
deallocate(hidro2, STAT = DeAllocateStatus)
deallocate(ILF, STAT = DeAllocateStatus)
deallocate(ILR, STAT = DeAllocateStatus)
deallocate(XI, STAT = DeAllocateStatus)
deallocate(Tt, STAT = DeAllocateStatus)
deallocate(Di, STAT = DeAllocateStatus)
deallocate(aux1, STAT = DeAllocateStatus)
deallocate(aux2, STAT = DeAllocateStatus)
deallocate(aux3, STAT = DeAllocateStatus)
deallocate(aux4, STAT = DeAllocateStatus)
deallocate(contribuicao_self, STAT = DeAllocateStatus)
deallocate(contribuicao_fisico, STAT = DeAllocateStatus)
deallocate(contribuicao_reciproco, STAT = DeAllocateStatus)
if(tmagper)then
deallocate(auxt, STAT = DeAllocateStatus)
deallocate(torquereal, STAT = DeAllocateStatus)
deallocate(torquereciproco, STAT = DeAllocateStatus)
deallocate(cof4, STAT = DeAllocateStatus)
deallocate(cof5, STAT = DeAllocateStatus)
deallocate(cof7, STAT = DeAllocateStatus)
end if
if(fmagper) then
deallocate(cof6, STAT = DeAllocateStatus)
deallocate(cof8, STAT = DeAllocateStatus)
deallocate(auxf, STAT = DeAllocateStatus)
deallocate(forcareal, STAT = DeAllocateStatus)
deallocate(forcareciproca, STAT = DeAllocateStatus)
end if
deallocate(ILF, STAT = DeAllocateStatus)
deallocate(ILR, STAT = DeAllocateStatus)
deallocate(XI, STAT = DeAllocateStatus)
deallocate(cof1, STAT = DeAllocateStatus)
deallocate(cof2, STAT = DeAllocateStatus)
deallocate(cof3, STAT = DeAllocateStatus)
if(leito)then
deallocate(usistema, STAT = DeAllocateStatus)
end if
if(grafmag)then
deallocate(magtempo, STAT = DeAllocateStatus)
deallocate(flutmag, STAT = DeAllocateStatus)
end if
deallocate(tempototal, STAT = DeAllocateStatus)
if(agregado_inicial) then
deallocate(centro_massa, STAT = DeAllocateStatus)
end if
deallocate(DIAM, STAT = DeAllocateStatus)
deallocate(beta, STAT = DeAllocateStatus)
deallocate(diarand, STAT = DeAllocateStatus)

write(*,*) ''
write(*,*) 'End of the processing module...'
write(*,*) ''
end subroutine main
