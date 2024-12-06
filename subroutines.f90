!*************************************************!
! 		     SIMMSUS			  ! 
!MODULE: subroutines				  !					         
!Last update: 16/07/2023			  !
!*************************************************!

module subroutines
use variables
contains

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: allocatevariables			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

subroutine allocatevariables
 
if(bifshear) then
allocate(gpvetor(npast))
end if
allocate(trap(npast))
allocate(X(rea,N,3))
allocate(U(rea,N,3))
allocate(W(rea,N,3))
if(browniano)then
allocate(FORCAS(6,rea,N,3))
allocate(TORQUES(3,rea,N,3))
else
allocate(FORCAS(5,rea,N,3))
allocate(TORQUES(2,rea,N,3))
end if
allocate(FT(rea,N,3))
allocate(Tt(rea,N,3))
allocate(Di(rea,N,3))
allocate(aux1(rea,N))
allocate(aux2(rea,N))
allocate(aux3(rea,N))
allocate(aux4(rea,N))
allocate(nr(nnr))
allocate(hidrodinamica_aux1(N,3))
allocate(hidrodinamica_aux2(N,3))
allocate(contribuicao_self(rea,N))
allocate(contribuicao_fisico(rea,N))
allocate(contribuicao_reciproco(rea,N))
allocate(hidro1(nb,3))
allocate(hidro2(nbr,3))
if(tmagper)then
allocate(auxt(N,3))
allocate(torquereal(nb,3))
allocate(torquereciproco(nbr,3))
allocate(cof4(2,10000))
allocate(cof5(2,10000))
allocate(cof7(2,10000))
end if
if(fmagper) then
allocate(cof6(2,10000))
allocate(cof8(2,10000))
allocate(auxf(N,3))
allocate(forcareal(nb,3))
allocate(forcareciproca(nbr,3))
end if
allocate(ILF(nb,3))
allocate(ILR(nbr,3))
allocate(XI(nb,rea,N,3))
allocate(cof1(2,10000))
allocate(cof2(2,10000))
allocate(cof3(2,10000))
if(leito)then
allocate(usistema(rea,3))
end if
if(grafmag)then
allocate(magtempo(3,npast))
allocate(flutmag(N,rea))
end if
allocate(tempototal(npast))
if(agregado_inicial) then
allocate(centro_massa(rea,3))
end if
allocate(DIAM(rea,N))
allocate(beta(rea,N))
allocate(diarand(rea*N))
allocate(campo(npast))
allocate(y(npast))

end subroutine allocatevariables

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: particle_distribution		  !					         
!Last update: 16/07/2023			  !
!*************************************************!

subroutine particle_distribution

if(polidispersidade)then

call randomica(1.5,2.5,diarand,(N*rea),8)

do j=1,rea
do i=1,N
DIAM(j,i)=diarand((j-1)*N + i)
end do
end do

else

do j=1,rea
do i=1,N
DIAM(j,i)=2.0
end do
end do

end if

do j=1,rea
do i=1,n
beta(j,i) = DIAM(j,i)/DIAM(j,1)
end do
end do

end subroutine particle_distribution

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: box_size				  !					         
!Last update: 16/07/2023			  !
!*************************************************!

subroutine box_size

if(agregado_inicial) then

ragreg=(N/phi)**(1.0/3.0)
l=100.0*ragreg
h=l

else
 l=((N/(razao*phi))*(4.0*3.1416)/(3.0))**(1.0/3.0)
 razao2=razao
 h=razao2*l
end if

end subroutine box_size

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: gera_arquivos			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine responsible for creating all the fi- !
! les used in the simulation 			  !
!*************************************************!
subroutine gera_arquivos(a,b,c)

logical a,b
integer c,i

! Subtitle of the units used for file storage

! From 1 to rea -> posicao.plt
! From rea+1 to 2*rea -> velocidade.plt
! From 2*rea+1 to 3*rea -> dipolo.plt
! From 3*rea+1 to 4*rea -> dip_parcial.plt
! 5*rea -> mag_tempo.plt
! 100*rea-> particula_teste.plt
! 200*rea-> distribuicao_agregados.plt
! 300*rea -> duffing_field.plt
! 1093*j para j=1,rea -> phi_local.plt
! 2012*j para j=1,rea -> energia.plt
! From 400*rea +1 to 400*rea + nfreq


! Creating a file that will be used to plot the trajectory of an arbitrary test particle

open(100*c,file='test_particle.plt')
write(100*c,*) 'Variables="X","Y","Z","Dx","Dy","Dz","t"'


if(bifurcation.or.bifshear) then
do i=1,nfreq
write(rea_char, '(I3)') i
open (400*c+i,file='magnetization'//rea_char//'.plt')
write(400*c+i,*) 'Variables="H","dH/dt","Mx","My","Mz","dMx","dMy","dMz","t"'
end do
end if

if(duffing) then
open(300*c,file='duffing_field.plt')
write(300*c,*)'Variables="t","H(t)","dH/dt"'
end if

if(fator)then
open(6*rea,file='structure_factor.plt')
end if

! Creating the files for each realization

if(a)then
do i=1,c
write(rea_char, '(I3)') i
if(.not.ovito)then
open (i,file='position'//rea_char//'.plt')
else
open (i,file='position'//rea_char//'.xyz')
end if
if(.not.ovito)then
if(gravadipolo) then
write(i,*) 'Variables="X","Y","Z","Dx","Dy","Dz","D"'
else
write(i,*) 'Variables="X","Y","Z","D"'
end if
end if
end do
end if

if(b)then
do i=1,c
write(rea_char, '(I3)') i
open (c+i,file='velocity'//rea_char//'.plt')
write(c+i,*) 'Variables="U","V","W"'
end do
end if

if(grafmag) then
  open(5*c,file='magnetization.plt')
  write(5*c,*)'Variables="H","dH/dt","Mx","My","Mz","dMx","dMy","dMz","t"'
end if

if(printphi) then
do i=1,rea
write(rea_char, '(I3)') i
open (1093*i,file='local_phi'//rea_char//'.plt')
write(1093*i,*) 'Variables="X","Y","Z","Phi"'
end do
end if

! Right now all the necessary files for storaging the position, 
! velocity, magnetization and other important physical variables 
! have been created!

write(*,*) 'All output files have been created'
write(*,*) ''

end subroutine gera_arquivos

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: tabelagreen			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!************************************************!
! Subroutine resposible for pre-calculating the  !
! Green functions used in the calculation of pe- !
! riodic interactions, due to hydrodynamic for-  !
! ces and magnetic forces and torques.           ! 
! This procedure is done to decrease the compu-  !
! tational cost, since theses functions are pre- !
! calculated only once and are simply called any-!
! time a force or torque on a particle is calcu- !
! lated through the code in future loops.        !
!************************************************!

subroutine tabelagreen(a,e,f,g,y)

real a   ! parameter "qsi" of the sum convergence
real RR1,RR2,RR3,RR4,SIG,SIG2,SIG3,SIG4,SIG5,SIG7 ! auxiliary parameters used to compute these functions
real b(2,10000), c(2,10000), d(2,10000), tb(2,10000), tc(2,10000),td(2,10000), tf(2,10000) ! periodic torques coefficient
integer f,g,i ! f = nb e g = nbr
real e,y ! e=l, y=h
real pi, RR7, RR5, te(2,10000)

pi=acos(-1.0)
	
! Making a table with all the coefficients

	do i=1,10000

	 b(1,i)=2.0+(i-1)*(((3.0**0.5)*e*(f**(1.0/3.0)))/9999.0)
	 c(1,i)=b(1,i)
 
 	if(tmagper)then
 	tb(1,i)=b(1,i)
 	tc(1,i)=b(1,i)
 	end if
 	if(fmagper)then
 	td(1,i)=b(1,i)
 	end if

! Calculating the powers of all the possible distances between the particles

	RR1=b(1,i)
	RR2=b(1,i)**2.0
	RR3=b(1,i)**3.0
	RR4=b(1,i)**4.0
	RR5=b(1,i)**5.0
	RR7=b(1,i)**7.0
	
	SIG=qsi
	
	SIG2=qsi**2.0
	SIG3=qsi**3.0
	SIG4=qsi**4.0
	SIG5=qsi**5.0
	SIG7=qsi**7.0

! Computing the functions used to calculate hydrodynamic interactions (in real space)
 	if(ligaih) then
 	b(2,i) = (0.75/RR1 + 0.5/RR3)*ERFC(SIG*RR1) + (4.0*SIG7*RR4 + 3.0*SIG3*RR2 - &
        	 & 20.0*SIG5*RR2 - 4.50*SIG + 14.0*SIG3 + SIG/RR2)*EXP(-SIG2*RR2)/SQRT(pi) 
      	c(2,i) = (0.75/RR1 - 1.5/RR3)*ERFC(SIG*RR1) + (-4.0*SIG7*RR4 - 3.0*SIG3*RR2 + &
        	 & 16.0*SIG5*RR2 + 1.50*SIG - 2.0*SIG3 - 3.0*SIG/RR2)*EXP(-SIG2*RR2)/SQRT(pi)
	end if
 ! Computing the functions used to calculate periodic magnetic torques (in real space)	
	if(tmagper)then
 	tb(2,i)=(ERFC(SIG*RR1)+ ((2.0*SIG*RR1)/(SQRT(pi)))*EXP(-SIG2*RR2))/RR3 
	tc(2,i)=(3.0*ERFC(SIG*RR1) + (((2.0*SIG*RR1)/(SQRT(pi)))*(3.0+(2.0*SIG2*RR2))*EXP(-SIG2*RR2)))/RR5
	end if
 ! Computing the functions used to calculate periodic magnetic forces (in real space)	
	if(fmagper)then
	td(2,i)=(15.0*ERFC(SIG*RR1)+ ((2.0*SIG*RR1)/(SQRT(pi)))*(15.0+(10.0*SIG2*RR2+4.0*SIG4*RR4))*EXP(-SIG2*RR2))/RR7 
	end if

! Same procedure done now for the reciprocal space

! Calculating all possible wave numbers

	if(g**(1.0/3.0).eq.5.0)then
	d(1,i)= (2.0*pi/e)+ ((i-1)*(((4.0*(3.0**0.5)*pi/e)-(2.0*pi/e))/9999.0))
	end if

	if(g**(1.0/3.0).eq.3.0)then
	d(1,i)= (2.0*pi/e)+ ((i-1)*(((2.0*(3.0**0.5)*pi/e)-(2.0*pi/e))/9999.0))
	end if
! Determining the 2th and 4th powers of these possible wave numbers

	RR2=d(1,i)**2.0
	RR4=d(1,i)**4.0

! Computing the functions used to calculate hydrodynamic interactions (in reciprocal space)	
	if(ligaih) then
	d(2,i)=(1.0/(e*e*y))*(1.0 - 1.0*RR2/3.0)*(1.0 + 0.25*RR2/SIG2 + 0.125*RR4/SIG4)*&
             & 6.0*pi*EXP(-0.25*RR2/SIG2)/RR2
	end if
 ! Computing the functions used to calculate periodic magnetic torques (in reciprocal space)	
	if(tmagper)then
	te(1,i)=d(1,i)
	te(2,i)=-(1.0/(e*e*y))*(4.0*pi*(EXP(-((pi/e)**2.0)*RR2/SIG2)))/1.0 !RR2
	end if
 ! Computing the functions used to calculate periodic magnetic forces (in reciprocal space)
	if(fmagper)then
	tf(1,i)=d(1,i)	
	tf(2,i)=((8.0*(pi**2.0))/((e**4.0)*(RR2)))*EXP(-((pi**2.0)*RR2)/(SIG2*(e**2.0)))/RR2 
	end if

	end do

	cof1=b
 	cof2=c
 	cof3=d

	 if(tmagper)then
	 cof4=tb
	 cof5=tc
	 cof7=te
	 end if

	 if(fmagper) then
	 cof6=td
	 cof8=tf
	 end if

end subroutine tabelagreen

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: estrutura_periodica		  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Suroutine responsible for creating a periodic   !
! strucutre (lattice) containing copies of a phy- !
! sical cell with the interacting particles repli-!
! cated through the real and reciprocal spaces    !
! using the Ewald summation technique (1921) des- !
! cribed in details in the works of Beenakker     !
! (1986), Abade (2002), Gontijo (2013)            !
!*************************************************!

subroutine estrutura_periodica

integer auxper(5), auxper2(5), auxper3(5)

pi=acos(-1.0)

! Creating the lattices indeces in the real space

auxper=0

if(nb**(1.0/3.0).eq.5.0)then
auxper(1)=1
auxper(2)=2
auxper(3)=3
auxper(4)=4
auxper(5)=5

auxper2(1)=0
auxper2(2)=5
auxper2(3)=10
auxper2(4)=15
auxper2(5)=20

auxper3(1)=0
auxper3(2)=25
auxper3(3)=50
auxper3(4)=75
auxper3(5)=100
end if

if(nb**(1.0/3.0).eq.3.0)then
auxper(1)=1
auxper(2)=2
auxper(3)=3

auxper2(1)=0
auxper2(2)=3
auxper2(3)=6

auxper3(1)=0
auxper3(2)=9
auxper3(3)=18
end if

do a=1,(nb**(1.0/3.0))
do b=1,(nb**(1.0/3.0))
do c=1,(nb**(1.0/3.0))

! Number of physical boxes

s=auxper3(a)+ auxper2(b) +auxper(c)

if(nb**(1.0/3.0).eq.5.0) then
ILF(s,1)=a-3
ILF(s,2)=b-3
ILF(s,3)=c-3
end if

if(nb**(1.0/3.0).eq.3.0) then
ILF(s,1)=a-2
ILF(s,2)=b-2
ILF(s,3)=c-2
end if

end do
end do
end do

! Creating the initial configuration of all the physical lattices

do a=1,nb
do b=1,rea
do i=1,N

XI(a,b,i,1)= X(b,i,1) + ILF(a,1)*l
XI(a,b,i,2)= X(b,i,2) + ILF(a,2)*l
XI(a,b,i,3)= X(b,i,3) + ILF(a,3)*h

end do
end do
end do

if(k.eq.1)then
open(872,file='condicao_inicial_periodica.plt')
write(872,*)'Variables= "X1","X2","X3"'
do a=1,nb
do b=1,N
write(872,*) XI(a,1,b,1),XI(a,1,b,2),XI(a,1,b,3)
end do
end do
end if

! Creating the lattice's indeces in the reciprocal space

if(nbr**(1.0/3.0).eq.5.0)then
auxper(1)=1
auxper(2)=2
auxper(3)=3
auxper(4)=4
auxper(5)=5

auxper2(1)=0
auxper2(2)=5
auxper2(3)=10
auxper2(4)=15
auxper2(5)=20

auxper3(1)=0
auxper3(2)=25
auxper3(3)=50
auxper3(4)=75
auxper3(5)=100
end if

if(nbr**(1.0/3.0).eq.3.0)then
auxper(1)=1
auxper(2)=2
auxper(3)=3

auxper2(1)=0
auxper2(2)=3
auxper2(3)=6

auxper3(1)=0
auxper3(2)=9
auxper3(3)=18
end if

do a=1,(nbr**(1.0/3.0))
do b=1,(nbr**(1.0/3.0))
do c=1,(nbr**(1.0/3.0))

! Number of reciprocal lattices

s=auxper3(a)+ auxper2(b) +auxper(c)

if(nbr**(1.0/3.0).eq.5.0) then
ILR(s,1)=a-3
ILR(s,2)=b-3
ILR(s,3)=c-3
end if

if(nbr**(1.0/3.0).eq.3.0) then
ILR(s,1)=a-2
ILR(s,2)=b-2
ILR(s,3)=c-2
end if

end do
end do
end do

end subroutine estrutura_periodica

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: distribui_dipolo			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine responsible for creating an initial  !
! distribution of the dipole moments of all parti-!
! cles in all realizations			  !
!*************************************************!

subroutine distribui_dipolo(a,b,c)

integer b,c ! b= number of realizations, c= number of particles
real a(b,c,3) ! Dipoles
real d,e ! 
integer f 

! If dipoles are distributed in an ordered way

if(dipolo_ordenado) then

do j=1,rea
do i=1,N
a(j,i,1)=0.0
a(j,i,2)=1.0
a(j,i,3)=0.0
end do
end do

else
! For a random dipole distribution
e=percentual*c
f=e

call randomica(-1.0,1.0,nr,(3*N*rea),3)

! If we are mixing magnetic particles with non magnetic ones...
if(mistura)then
 !$OMP PARALLEL DO
do j=1,b
do i=1,f

a(j,i,1)= 0.0
a(j,i,2)= 0.0
a(j,i,3)= 0.0

end do
end do
 !$OMP END PARALLEL DO

 !$OMP PARALLEL DO
do j=1,b
do i=f+1,c
a(j,i,1)= nr((i*2+(i-2)+(c*3*(j-1))))
a(j,i,2)= nr((i*2+(i-1)+(c*3*(j-1))))
a(j,i,3)= nr((i*2+(i)+(c*3*(j-1))))
end do
end do
 !$OMP END PARALLEL DO

! Normalizing the vectors

 !$OMP PARALLEL DO
do j=1,b
do i=1,f

a(j,i,1)= 0.0
a(j,i,2)= 0.0
a(j,i,3)= 0.0

end do
end do
 !$OMP END PARALLEL DO

 !$OMP PARALLEL DO
do j=1,b
do i=f+1,N
d=((a(j,i,1)**2.0)+(a(j,i,2)**2.0)+(a(j,i,3)**2.0))**0.5
a(j,i,1)=a(j,i,1)/d
a(j,i,2)=a(j,i,2)/d
a(j,i,3)=a(j,i,3)/d
end do
end do
 !$OMP END PARALLEL DO
Di=a

else
! If all the particles are magnetic particles, then...

 !$OMP PARALLEL DO
do j=1,b
do i=1,c
a(j,i,1)= nr((i*2+(i-2)+(c*3*(j-1))))
a(j,i,2)= nr((i*2+(i-1)+(c*3*(j-1))))
a(j,i,3)= nr((i*2+(i)+(c*3*(j-1))))
end do
end do
 !$OMP END PARALLEL DO

! Normalizing the vectors

 !$OMP PARALLEL DO
do j=1,rea
do i=1,N
d=((a(j,i,1)**2.0)+(a(j,i,2)**2.0)+(a(j,i,3)**2.0))**0.5
a(j,i,1)=a(j,i,1)/d
a(j,i,2)=a(j,i,2)/d
a(j,i,3)=a(j,i,3)/d
end do
end do
 !$OMP END PARALLEL DO
Di=a
end if
end if

end subroutine distribui_dipolo

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: condicao_inicial			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine responsible for creating the initial !
! particle distribution for all simultaneous 	  !
! numerical experiments.			  !
!						  !
! Several possible initial configurations are im- !
! plemented here. These include:		  !
!						  !
! 1 - An initial random distribution		  !
! 2 - An ordered NxNxN array of particles 	  !
! 3 - An initial spherical aggregate		  !
!*************************************************!

subroutine condicao_inicial

integer auxiliar1, loop, loop2

! If you have to continue a previous simulation, then get into this "if clause"
if(continua) then

! Opening the "particula_teste.plt" file in which we print the trajectory of an
! arbitrary test particle

open(100*rea,file='particula_teste.plt', STATUS='OLD')

! Defining important file's format
509 FORMAT(F30.4,F30.4,F30.4)
666 FORMAT(F30.4,F30.4,F30.4,F30.4,F30.4,F30.4,F30.4)

! Opening velocity files 

do i=1,rea
write(rea_char, '(I3)') i
open (rea+i,file='velocidade'//rea_char//'.plt', STATUS='OLD')
end do

! Reading the velocities of all the particles in several files (one file per realization)

do j=1,rea
read (rea+j,'(A)') linha1
do k=1,(iter/n2) 
read(rea+j,'(A)') linha2
do i=1,N
read(rea+j,509)U(j,i,1),U(j,i,2),U(j,i,3)
end do
end do
end do

! Opening position files

do i=1,rea
write(rea_char, '(I3)') i
open (i,file='posicao'//rea_char//'.plt', STATUS='OLD')
end do

! Reading the position and dipole moments of all the particles (in order to define a configuration from which the
! simulation must continue, we need to know not only the final position of the particles, but also their orientation

do j=1,rea
read (j,'(A)') linha1
do k=1,(iter/n2) 
read(j,'(A)') linha2
do i=1,N
read(j,666)X(j,i,1),X(j,i,2),X(j,i,3),Di(j,i,1),Di(j,i,2),Di(j,i,3),DIAM(j,i)
end do
end do
end do

else

! If we are starting a new simulation, then we will check if the initial condition is a spherical aggregate

if(agregado_inicial) then

open(666,file='vel_tempo.plt')
write(666,*)'Variables="t","Ux","Uy","Uz"'

! Calculating the aggregate's radius

ragreg=(N/phi)**(1.0/3.0)

! Inserting the aggregate in a giant lattice to avoid periodicity condition 
l=100.0*ragreg
h=l

 call randomica(-1.0,1.0,nr,(3*N*rea),5)

! Defining the center of the aggregate

xcentro=l/2.0
ycentro=l/2.0
zcentro=h - 2.0*ragreg

do j=1,rea
do i=1,N

! Range of X position of the particles in the aggregate

xmin= l/2.0 - ragreg
xmax= xmin + 2.0*ragreg

 X(j,i,1)= xcentro + (xmax-xmin)*0.5*nr((i*2+(i-2)+(N*3*(j-1))))

! Using the circle equation to distribute the particles in the y direction (sphere projection on a xy plane)

ymin=ycentro-(ragreg**2.0 -(X(j,i,1)-xcentro)**2.0)**0.5
ymax=ycentro+(ragreg**2.0 -(X(j,i,1)-xcentro)**2.0)**0.5

 X(j,i,2)= ycentro +(ymax-ymin)*0.5*nr((i*2+(i-1)+(N*3*(j-1))))

! Using the circle equation to define the z positions considering now the sphere projection on planes zy and zx

zmin=zcentro-(ragreg**2.0 -(X(j,i,1)-xcentro)**2.0 -(X(j,i,2)-ycentro)**2.0)**0.5
zmax=zcentro+(ragreg**2.0 -(X(j,i,1)-xcentro)**2.0 -(X(j,i,2)-ycentro)**2.0)**0.5

 X(j,i,3)= zcentro +(zmax-zmin)*0.5*nr((i*2+(i)+(N*3*(j-1))))

end do
end do

write(*,*) 'O raio do agregado e:',ragreg
write(*,*) 'A largura do box e:', l
write(*,*) 'A altura do box e:',h

! Checking for particle overlaping in the initial condition

118 call randomica(-1.0,1.0,nr,(3*N*rea),6)

! Calculating the center of mass of the aggregate

 centro_massa(k,1)=sum(X(k,:,1))/N
 centro_massa(k,2)=sum(X(k,:,2))/N
 centro_massa(k,3)=sum(X(k,:,3))/N

else

! If you want to create an ordered distribution then...

if(ordenado) then

 !$OMP PARALLEL DO
do j=1,rea
loop=0
loop2=0
do i=1,N

reale= (i-1.0)/(n**(1.0/3.0))
inteiro= (i-1)/(n**(1.0/3.0))

reale2= (loop+1.0)/(n**(1.0/3.0))
inteiro2= (loop+1)/(n**(1.0/3.0))

if(reale.ne.0.0)then
if(reale.eq.inteiro) then
loop=loop+1
end if
end if

if(reale2.ne.1.0) then
if(reale2.eq.inteiro2) then
loop2=loop2+1
end if
end if

auxiliar1=i/(N**(2.0/3.0)) 

 a=i-loop*(N**(1.0/3.0))
 b=1+loop-auxiliar1*(N**(1.0/3.0))
 c=i/(N**(2.0/3.0)) 

X(j,i,1)= l/(2.0*(N**(1.0/3.0))) + (a)*(l-(l/(N**(1.0/3.0))))/((N**(1.0/3.0))-1.0)
X(j,i,2)= l/(2.0*(N**(1.0/3.0))) + (b)*(l-(l/(N**(1.0/3.0))))/((N**(1.0/3.0))-1.0)
X(j,i,3)= h/(2.0*(N**(1.0/3.0))) + (c)*(h-(h/(N**(1.0/3.0))))/((N**(1.0/3.0))-1.0)


end do
end do
 !$OMP END PARALLEL DO

else

! If you want to generate a random initial distribution, then you must call the random number generation routine...

 call randomica(0.0,1.0,nr,(3*N*rea),1)

 !$OMP PARALLEL DO
do j=1,rea
do i=1,N
 X(j,i,1)=0.0+(l-0.0)*nr((i*2+(i-2)+(N*3*(j-1))))
 X(j,i,2)=0.0+(l-0.0)*nr((i*2+(i-1)+(N*3*(j-1))))
 X(j,i,3)=0.0+(h-0.0)*nr((i*2+(i)+(N*3*(j-1))))
end do
end do
 !$OMP END PARALLEL DO

! Now that we already created an initial random configuration, we must check for particle overlaps

! We call the random number generation subroutine

 !$OMP PARALLEL DO
do k=1,rea
do i=1,N
do j=1,N

103 call randomica(-1.0,1.0,nr,3,2)

if(i.ne.j)then
! Calculate the distance between all the pair of particles in the suspension
r=((X(k,i,1)-X(k,j,1))**2.0+(X(k,i,2)-X(k,j,2))**2.0+(X(k,i,3)-X(k,j,3))**2.0)**0.5
! If at any point the distance is smaller then 0.01*a, then we give a Brownian kick in the overlaped particles
if(r.le.2.01)then

! Creating a small random vector used to "kick" one of the overlaped particles

nr1=nr(1)
nr2=nr(2)
nr3=nr(3)

modrand=((nr1**2.0)+(nr2**2.0)+(nr3**2.0))**0.5

! Normalizing this vector

nr1=nr1/modrand
nr2=nr2/modrand
nr3=nr3/modrand

nr1=nr1*0.25
nr2=nr2*0.25
nr3=nr3*0.25

! Changing the position of one of the overlaped particles with a Brownian "kick"

X(k,i,1)=X(k,i,1)+nr1
X(k,i,2)=X(k,i,2)+nr2
X(k,i,3)=X(k,i,3)+nr3

end if
! Calculate the new distance between the "problematic" particles
r=((X(k,i,1)-X(k,j,1))**2.0+(X(k,i,2)-X(k,j,2))**2.0+(X(k,i,3)-X(k,j,3))**2.0)**0.5

! If at any point the particles are still overlaped we give new 
! Brownian kicks in the particles that present this problem until the
! situation is solved

if(r.le.2.01)then
go to 103
end if

! Imposing conditions to avoid particles outside the box

if(X(k,i,1).lt.0) then
X(k,i,1)=abs(X(k,i,1))
go to 103
end if

if(X(k,i,2).lt.0) then
X(k,i,2)=abs(X(k,i,2))
go to 103
end if

if(X(k,i,3).lt.0) then
X(k,i,3)=abs(X(k,i,3))
go to 103
end if

if(X(k,i,1).gt.(l)) then
X(k,i,1)=X(k,i,1) - l
go to 103
end if

if(X(k,i,2).gt.(l)) then
X(k,i,2)=X(k,i,2) - l
go to 103
end if

if(X(k,i,3).gt.(h)) then
X(k,i,3)=X(k,i,3) - h
go to 103
end if

end if
end do
end do
end do
 !$OMP END PARALLEL DO
end if

! Lets define the initial velocity of the particles as the Stokes velocity (mobility problem)

 !$OMP PARALLEL DO
do j=1,rea
do i=1,N
U(j,i,1)=0.0
U(j,i,2)=0.0
U(j,i,3)=-1.0
end do
end do
 !$OMP END PARALLEL DO
end if

end if

end subroutine condicao_inicial

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: field_excitations			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine resonsible for pre-calculating the	  !
! magnetic field excitation and shear-rate ramp in!
! case the user wants to build a bifurcation dia- ! 
! gram of the magnetization response		  !
!						  !
! The following possibilities are considered for  !
! the magnetic field H(t)			  !	
!						  !	
! 1 - H(t) as a solution of the Duffing oscilator !
! 2 - H(t) = cos(w1*t) + cos(w2*t) 		  !
! 3 - H(t) = sin(w*t) [oscillatory field]	  !
! 4 - H(t) = sin(i*w*t) [i varying periodically]  !
! 5 - H(t) = sin(w*t)ex + cos(w*t)ey [rotating]   !
!					          !
! Field excitation number 4 implies that the field!
! frequency will periodically change as the simula!
! tion evolves. This is useful to build bifurca-  !
! tion diagrams. 				  !
! 						  !
! This subroutine also produces an oscillatory    !
! shear ramp.				   	  !
!*************************************************!
 
subroutine field_excitations

510 FORMAT(F30.4,F30.4,F30.4)

if(duffing) then

! Condições iniciais para o campo H e sua derivada H'=y

y(1)=0.0
campo(1)=1.0
tempototal(1)=0.0

write(300*rea,510) tempototal(1),campo(1),y(1)

! Solução por meio do método de Ruge-Kutta de 4 ordem da excitação H(t) e de sua derivada H'(t)

do k=2,npast

k1=C4*cos(freqcampo*tempototal(k-1)) - (C1*y(k-1) + C2*campo(k-1) + C3*(campo(k-1)**3.0))
k2=C4*cos(freqcampo*(tempototal(k-1)+ 0.5*dt))-(C1*(y(k-1)+ 0.5*dt*k1)+C2*campo(k-1)+C3*(campo(k-1)**3.0))
k3=C4*cos(freqcampo*(tempototal(k-1)+ 0.5*dt))-(C1*(y(k-1)+ 0.5*dt*k2)+C2*campo(k-1)+C3*(campo(k-1)**3.0))
k4=C4*cos(freqcampo*(tempototal(k-1)+dt))-(C1*(y(k-1)+ dt*k3)+C2*campo(k-1)+C3*(campo(k-1)**3.0))

g1=y(k-1)
g2=y(k-1) + dt*0.5*g1
g3=y(k-1) + dt*0.5*g2
g4=y(k-1) + dt*g3

y(k) = y(k-1) + (dt/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)
campo(k) = campo(k-1) + (dt/6.0)*(g1 + 2.0*g2 + 2.0*g3 + g4)
tempototal(k) = tempototal(k-1) + dt 

write(300*rea,510) tempototal(k),campo(k),y(k)

end do

print *,'******************************************************************************'
print *,'*                                                                            *'
print *,'*             DUFFING HARMONIC EXCITATION SUCCESSFULLY PRE-CALCULATED        *'
print *,'*                                                                            *'
print *,'******************************************************************************'
print *,''

end if

! Caso estejamos trabalhando com uma excitação de campo do tipo batimento:

if(beating) then

! Condições iniciais para o campo H e sua derivada H'=y

y(1)=0.0
campo(1)=2.0
tempototal(1)=0.0

do k=2,npast
tempototal(k) = tempototal(k-1) + dt 
campo(k) = cos(freqcampo*tempototal(k)) + cos(freqbeat*tempototal(k))
y(k) = -freqcampo*sin(freqcampo*tempototal(k)) - freqbeat*sin(freqbeat*tempototal(k))
end do

print *,'******************************************************************************'
print *,'*                                                                            *'
print *,'*           BEATING PATTERN FIELD EXCITATION SUCCESSFULLY PRE-CALCULATED     *'
print *,'*                                                                            *'
print *,'******************************************************************************'
print *,''

end if

! No caso em que temos um campo oscilatorio e não estamos lidando com excitações do tipo
! duffing ou por padrões de batimento, estamos então no contexto de uma excitação harmônica do tipo seno

if(oscilacampo) then
if(.not.duffing) then
if(.not.beating) then

! Condições iniciais para o campo H e sua derivada H'=y

y(1)=freqcampo
campo(1)=0.0
tempototal(1)=0.0

! Calculando o intervalo de tempo de atuação de cada frequência (para construção dos diagramas de bifurcação)

intervalo=npast*dt/nfreq
multiplofreq=0

do k=2,npast
! MONTANDO O CAMPO MAGNÉTICO APLICADO CONSIDERANDO A MUDANÇA DA FREQUÊNCIA A CADA PERÍODO PRÉ-DETERMINADO PELO USUÁRIO
if(bifurcation) then
tempototal(k) = tempototal(k-1) + dt 
contfreqinteiro1=tempototal(k-1)/intervalo
contfreqinteiro2=tempototal(k)/intervalo
if(contfreqinteiro1.ne.contfreqinteiro2) then
multiplofreq=multiplofreq+1 
end if
campo(k)=sin((freqcampo+(bifmax-freqcampo)*multiplofreq)*tempototal(k)) 
y(k)=(freqcampo + (bifmax-freqcampo)*multiplofreq)*cos((freqcampo + (bifmax-freqcampo)*multiplofreq)*tempototal(k))
else
tempototal(k) = tempototal(k-1) + dt 
campo(k) = sin(freqcampo*tempototal(k)) 
y(k) = freqcampo*cos(freqcampo*tempototal(k))
end if
end do
multiplofreq=0


print *,'******************************************************************************'
print *,'*                                                                            *'
print *,'*                HARMONIC FIELD EXCITATION SUCCESSFULLY PRE-CALCULATED       *'
print *,'*                                                                            *'
print *,'******************************************************************************'
print *,''

end if
end if
end if

! MONTANDO A RAMPA DE SHEAR DO CISALHAMENTO OSCILATÓRIO PARA O CASO DE DIAGRAMA DE BIFURCAÇÃO PARA GAMMA PONTO
do k=2,npast
if(bifshear) then
tempototal(k) = tempototal(k-1) + dt 
contfreqinteiro1=tempototal(k-1)/intervalo
contfreqinteiro2=tempototal(k)/intervalo
if(contfreqinteiro1.ne.contfreqinteiro2) then
multiplofreq=multiplofreq+1 
end if
gpvetor(k)= shearrate + multiplofreq*((bifmax-shearrate)/nfreq) 
end if
end do
multiplofreq=0

end subroutine field_excitations

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: campo_phi				  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine responsible for creating a tecplot   !	
! format file with a 3D field of the local volume !
! fraction of particles inside the simulation box !
! domain. This subroutine divides the simulation  !
! box space into 10x10x10 cells, counts the number!
! of particles in each cell and divide by the num-!
! ber of particles within the suspension space.   !
!*************************************************!

subroutine campo_phi(a,s)
real dx,dy,dz
integer nc, ncubo ! number of sub-divisions in the central cell
integer i,j,k,o,p,a 
integer auxiliar, s
integer, allocatable :: auxper(:), auxper2(:), auxper3(:)
real, allocatable :: philocal(:,:)
real reali, realj, realk, realn

nc=10
ncubo=nc**3.0
realn=nc

allocate(philocal(a,ncubo))
allocate(auxper(nc))
allocate(auxper2(nc))
allocate(auxper3(nc))

dx = l/realn
dy = l/realn
dz = h/realn

philocal=0.0

do i=1,nc
auxper(i)=i
auxper2(i)= (i-1)*nc
auxper3(i)= (i-1)*(nc**2.0)
end do

do i=1,nc
do j=1,nc
do k=1,nc

do p=1,rea
do o=1,N

if(X(p,o,1).ge.(i-1)*dx)then
if(X(p,o,1).le.(i*dx))then

if(X(p,o,2).ge.(j-1)*dy)then
if(X(p,o,2).le.(j*dy))then

if(X(p,o,3).ge.(k-1)*dz)then
if(X(p,o,3).le.(k*dz))then

auxiliar = auxper3(k)+auxper2(j)+auxper(i)

philocal(p,auxiliar)=philocal(p,auxiliar) + 1

end if
end if
end if
end if
end if
end if

end do
end do

end do
end do
end do

do p=1,rea

! If it works I should change this "10" by "nc" in text mode
write(1093*p,*) 'zone F=POINT,I=10,J=10,K=10'

do i=1,nc
do j=1,nc
do k=1,nc
auxiliar = auxper3(k)+auxper2(j)+auxper(i)

reali=i
realj=j
realk=k

!This expression (1.0 + 1.2/realn) is just a numerical artificie used for visualization purposes in tecplot

write(1093*p,*) (reali-1.0)*dx*(1.0 + (1.2/realn)), (realj-1.0)*dy*(1.0 + (1.2/realn)), (realk-1.0)*(1.0 + (1.2/realn))*dz, philocal(p,auxiliar)


end do
end do
end do

end do

deallocate(philocal)
end subroutine campo_phi


!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: brownian				  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine responsible for implementing Brownian!
! forces and torques				  !
!*************************************************!

subroutine brownian
 call randomica(-1.0,1.0,nr,(3*N*rea),3+k)
 !$OMP PARALLEL DO
do j=1,rea
do i=1,N

nr1= nr((i*2+(i-2)+(N*3*(j-1))))
nr2= nr((i*2+(i-1)+(N*3*(j-1))))
nr3= nr((i*2+(i)+(N*3*(j-1))))

modrand=((nr1**2.0)+(nr2**2.0)+(nr3**2.0))**0.5
nr1=nr1/modrand
nr2=nr2/modrand
nr3=nr3/modrand

if(gravidade) then
FORCAS(6,j,i,1)= (beta(j,i)**(-2.0))*((6.0/(Pe*dt))**0.5)*nr1
FORCAS(6,j,i,2)= (beta(j,i)**(-2.0))*((6.0/(Pe*dt))**0.5)*nr2
FORCAS(6,j,i,3)= (beta(j,i)**(-2.0))*((6.0/(Pe*dt))**0.5)*nr3
else
if(shear) then
FORCAS(6,j,i,1)= (beta(j,i)**(-2.0))*((6.0/(Pe*dt))**0.5)*nr1
FORCAS(6,j,i,2)= (beta(j,i)**(-2.0))*((6.0/(Pe*dt))**0.5)*nr2
FORCAS(6,j,i,3)= (beta(j,i)**(-2.0))*((6.0/(Pe*dt))**0.5)*nr3
else
FORCAS(6,j,i,1)= (beta(j,i)**(-2.0))*((6.0/(dt))**0.5)*nr1
FORCAS(6,j,i,2)= (beta(j,i)**(-2.0))*((6.0/(dt))**0.5)*nr2
FORCAS(6,j,i,3)= (beta(j,i)**(-2.0))*((6.0/(dt))**0.5)*nr3
end if
end if
end do
end do
 !$OMP END PARALLEL DO

if(torque)then

 call randomica(-1.0,1.0,nr,(3*N*rea),npast+k)
 !$OMP PARALLEL DO
do j=1,rea
do i=1,N
nr1= nr((i*2+(i-2)+(N*3*(j-1))))
nr2= nr((i*2+(i-1)+(N*3*(j-1))))
nr3= nr((i*2+(i)+(N*3*(j-1))))

modrand=((nr1**2.0)+(nr2**2.0)+(nr3**2.0))**0.5
nr1=nr1/modrand
nr2=nr2/modrand
nr3=nr3/modrand

if(gravidade)then
TORQUES(3,j,i,1)= (beta(j,i)**(-2.0))*((9.0/(2.0*Pe*dt))**0.5)*nr1
TORQUES(3,j,i,2)= (beta(j,i)**(-2.0))*((9.0/(2.0*Pe*dt))**0.5)*nr2
TORQUES(3,j,i,3)= (beta(j,i)**(-2.0))*((9.0/(2.0*Pe*dt))**0.5)*nr3
else
if(shear) then
TORQUES(3,j,i,1)= (beta(j,i)**(-2.0))*((9.0/(2.0*Pe*dt))**0.5)*nr1
TORQUES(3,j,i,2)= (beta(j,i)**(-2.0))*((9.0/(2.0*Pe*dt))**0.5)*nr2
TORQUES(3,j,i,3)= (beta(j,i)**(-2.0))*((9.0/(2.0*Pe*dt))**0.5)*nr3
else
TORQUES(3,j,i,1)= (beta(j,i)**(-2.0))*((9.0/(2.0*dt))**0.5)*nr1
TORQUES(3,j,i,2)= (beta(j,i)**(-2.0))*((9.0/(2.0*dt))**0.5)*nr2
TORQUES(3,j,i,3)= (beta(j,i)**(-2.0))*((9.0/(2.0*dt))**0.5)*nr3
end if
end if

end do
end do
 !$OMP END PARALLEL DO

end if

end subroutine brownian

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: gravity				  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine responsible for implementing gravita-!
! tional forces					  !
!*************************************************!
subroutine gravity
 !$OMP PARALLEL DO
do j=1,rea
do i=1,N
if(gravidade)then
! Gravitational forces
FORCAS(3,j,i,1)=0.0
FORCAS(3,j,i,2)=0.0
FORCAS(3,j,i,3)=-beta(j,i)**3.0
else
! Gravitational forces
FORCAS(3,j,i,1)=0.0
FORCAS(3,j,i,2)=0.0
FORCAS(3,j,i,3)=0.0
end if
end do
end do
 !$OMP END PARALLEL DO
end subroutine gravity

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: repulsion				  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine responsible for computing repulsive  !
! forces (lubrication) and contact forces for  	  !
! overlapped particles using the Hertz model	  !
!*************************************************!

subroutine repulsion

 !$OMP PARALLEL DO
do j=1,rea
do i=1,N
do q=1,N
if(i.ne.q) then
! Calculating the distance between all the particles in the suspension (pair by pair)
r=(((X(j,i,1)-X(j,q,1))**2.0)+((X(j,i,2)-X(j,q,2))**2.0)+((X(j,i,3)-X(j,q,3))**2.0))**0.5

! If the distance is less than 2.2, i.e. less then 0.2*a, where a denotes the particle radius, 
! then the repulsive force is turned on. It is important to notice that this force is turned off 
! in case of particle overlap (this condition is really difficult to occur for particles 
! wihtout inertia (mobility problem), because in this context we use a contact (Hertz) force.

if(r.le.2.2) then
if(r.ge.2.0) then
FORCAS(1,j,i,1)= 10.0*exp(-r/0.01)*(X(j,i,1)-X(j,q,1)/r)
FORCAS(1,j,i,2)= 10.0*exp(-r/0.01)*(X(j,i,2)-X(j,q,2)/r)
FORCAS(1,j,i,3)= 10.0*exp(-r/0.01)*(X(j,i,3)-X(j,q,3)/r)
end if
end if
if(r.ge.2.2) then
if(r.le.2.0) then
FORCAS(1,j,i,1)=0.0
FORCAS(1,j,i,2)=0.0
FORCAS(1,j,i,3)=0.0
end if
end if

! Here we calculate contact forces for overlapped particles

if(r.le.(beta(j,i)+beta(j,q))) then
Eij=abs(2.0-r)
FORCAS(2,j,i,1)= (100.0*Eij**(3.0/2.0))*(X(j,i,1)-X(j,q,1))/Eij
FORCAS(2,j,i,2)= (100.0*Eij**(3.0/2.0))*(X(j,i,2)-X(j,q,2))/Eij
FORCAS(2,j,i,3)= (100.0*Eij**(3.0/2.0))*(X(j,i,3)-X(j,q,3))/Eij
else
FORCAS(2,j,i,1)=0.0
FORCAS(2,j,i,2)=0.0
FORCAS(2,j,i,3)=0.0
end if

end if
end do
end do
end do
 !$OMP END PARALLEL DO
end subroutine repulsion

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: forca_magnetica			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine responsible for computing long-range !
! non-periodic dipolar interactions between the	  !
! particles					  !
!*************************************************!

subroutine forca_magnetica

integer auxiliary

if(mistura)then
 auxiliary=(percentual*N)+1
else
 auxiliary=1
end if

if(gravidade)then
lambda=alpha2*6.0/Pe
else
if(shear) then
lambda=alpha2*6.0/Pe
else
lambda=alpha2*6.0
end if
end if

do j=1,rea
do i=auxiliary,N
do q=auxiliary,N
if(i.ne.q) then
! Checking the distance between the particles
r=(((X(j,i,1)-X(j,q,1))**2.0)+((X(j,i,2)-X(j,q,2))**2.0)+((X(j,i,3)-X(j,q,3))**2.0))**0.5
! If particles are to close we turn off magnetic interactions to avoid overlap
if(r.le.2.2) then
aux1(j,q)=0.0
aux2(j,q)=0.0
aux3(j,q)=0.0
else
! Calculating the vector R_{ij} which connects a particle i to a particle j
rij(1)=X(j,i,1)-X(j,q,1)
rij(2)=X(j,i,2)-X(j,q,2)
rij(3)=X(j,i,3)-X(j,q,3)
! Normalizing vector R_{ij}
modrij=((rij(1)**2.0)+(rij(2)**2.0)+(rij(3)**2.0))**0.5
rij(1)=rij(1)/modrij
rij(2)=rij(2)/modrij
rij(3)=rij(3)/modrij

termo1=((Di(j,i,1)*Di(j,q,1))+(Di(j,i,2)*Di(j,q,2))+(Di(j,i,3)*Di(j,q,3)))*rij(1)
termo2=((Di(j,i,1)*rij(1))+(Di(j,i,2)*rij(2))+(Di(j,i,3)*rij(3)))*Di(j,q,1)
termo3=((Di(j,q,1)*rij(1))+(Di(j,q,2)*rij(2))+(Di(j,q,3)*rij(3)))*Di(j,i,1)
termo4=(((Di(j,i,1)*rij(1))+(Di(j,i,2)*rij(2))   &
+(Di(j,i,3)*rij(3)))*((Di(j,q,1)*rij(1))   &
+(Di(j,q,2)*rij(2))+(Di(j,q,3)*rij(3))))*rij(1)

aux1(j,q)=(lambda/(1.0*(r**4.0)))*(termo1+termo2+termo3-5.0*termo4)

termo1=((Di(j,i,1)*Di(j,q,1))+(Di(j,i,2)*Di(j,q,2))+(Di(j,i,3)*Di(j,q,3)))*rij(2)
termo2=((Di(j,i,1)*rij(1))+(Di(j,i,2)*rij(2))+(Di(j,i,3)*rij(3)))*Di(j,q,2)
termo3=((Di(j,q,1)*rij(1))+(Di(j,q,2)*rij(2))+(Di(j,q,3)*rij(3)))*Di(j,i,2)
termo4=(((Di(j,i,1)*rij(1))+(Di(j,i,2)*rij(2))   &
+(Di(j,i,3)*rij(3)))*((Di(j,q,1)*rij(1))   &
+(Di(j,q,2)*rij(2))+(Di(j,q,3)*rij(3))))*rij(2)

aux2(j,q)=(lambda/(1.0*(r**4.0)))*(termo1+termo2+termo3-5*termo4)

termo1=((Di(j,i,1)*Di(j,q,1))+(Di(j,i,2)*Di(j,q,2))+(Di(j,i,3)*Di(j,q,3)))*rij(3)
termo2=((Di(j,i,1)*rij(1))+(Di(j,i,2)*rij(2))+(Di(j,i,3)*rij(3)))*Di(j,q,3)
termo3=((Di(j,q,1)*rij(1))+(Di(j,q,2)*rij(2))+(Di(j,q,3)*rij(3)))*Di(j,i,3)
termo4=(((Di(j,i,1)*rij(1))+(Di(j,i,2)*rij(2))+   &
(Di(j,i,3)*rij(3)))*((Di(j,q,1)*rij(1))   &
+(Di(j,q,2)*rij(2))+(Di(j,q,3)*rij(3))))*rij(3)

aux3(j,q)=(lambda/(1.0*(r**4.0)))*(termo1+termo2+termo3-5.0*termo4)

end if
end if
end do

FORCAS(4,j,i,1)=sum(aux1(j,:))
FORCAS(4,j,i,2)=sum(aux2(j,:))
FORCAS(4,j,i,3)=sum(aux3(j,:))

aux1=0.0
aux2=0.0
aux3=0.0
end do
end do

end subroutine forca_magnetica

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: campo_externo			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine responsible for computing magnetic	  !
! forces between an external magnetic field and	  !
! each particle in suspension.  		  !
!*************************************************!

subroutine campo_externo
integer auxiliary

! These subroutine may consider a magnetic fluidized bed mixing magnetic and non-magnetic particles, the logical 
! variable "mistura" activates this possibility. If mistura = TRUE then we have to set zero dipoles for a certain 
! percentage "percentual" of the particles.

if(mistura)then
 auxiliary=(percentual*N)+1
else
 auxiliary=1
end if

 !$OMP PARALLEL DO
do j=1,rea
do i=auxiliary,N
! posicao_campo = 1 -> Applied field in the lower wall
if(posicao_campo.eq.1)then
if(Pe.eq.0.0) then
if(X(j,i,3).ge.(2.0)) then
FORCAS(5,j,i,3)= 2.0*(alpha)*(Di(j,i,3))/(1.0*((X(j,i,3))**3.0))
else
FORCAS(5,j,i,3)= 2.0*(alpha)*(Di(j,i,3))/(1.0*(2.0**3.0))
end if
else
if(X(j,i,3).ge.(2.0)) then
FORCAS(5,j,i,3)= 2.0*(alpha)*(Di(j,i,3))/(Pe*((X(j,i,3))**3.0))
else
FORCAS(5,j,i,3)= 2.0*(alpha)*(Di(j,i,3))/(Pe*(2.0**3.0))
end if
end if
end if

! posicao_campo = 2 -> Applied field on the upper wall

if(posicao_campo.eq.2) then
if(Pe.eq.0.0) then
if(X(j,i,3).le.(h-2.0)) then
FORCAS(5,j,i,3)= 2.0*(alpha)*(Di(j,i,3))/(1.0*((h-X(j,i,3))**3.0))
else
FORCAS(5,j,i,3)= 2.0*(alpha)*(Di(j,i,3))/(1.0*((h-2.0)**3.0))
end if
else
if(X(j,i,3).le.(h-2.0)) then
FORCAS(5,j,i,3)= 2.0*(alpha)*(Di(j,i,3))/(Pe*((h-X(j,i,3))**3.0))
else
FORCAS(5,j,i,3)= 2.0*(alpha)*(Di(j,i,3))/(Pe*((h-2.0)**3.0))
end if
end if
end if

! posicao_campo = 3 -> Applied field on the right side

if(posicao_campo.eq.3) then
if(Pe.eq.0.0) then
if(X(j,i,1).le.(l-2.0)) then
FORCAS(5,j,i,1)= 2.0*(alpha)*(Di(j,i,1))/(1.0*((l-X(j,i,1))**3.0))
else
FORCAS(5,j,i,1)= 2.0*(alpha)*(Di(j,i,1))/(1.0*((l-2.0)**3.0))
end if
else
if(X(j,i,1).le.(l-2.0)) then
FORCAS(5,j,i,1)= 2.0*(alpha)*(Di(j,i,1))/(Pe*((l-X(j,i,1))**3.0))
else
FORCAS(5,j,i,1)= 2.0*(alpha)*(Di(j,i,1))/(Pe*((l-2.0)**3.0))
end if
end if
end if

! posicao_campo = 4 -> Applied field on the left side

if(posicao_campo.eq.4) then
if(Pe.eq.0.0)then
if(X(j,i,1).le.(l-2.0)) then
FORCAS(5,j,i,1)= -2.0*(alpha)*(Di(j,i,1))/(1.0*((l-X(j,i,1))**3.0))
else
FORCAS(5,j,i,1)= -2.0*(alpha)*(Di(j,i,1))/(1.0*((l-2.0)**3.0))
end if
else
if(X(j,i,1).le.(l-2.0)) then
FORCAS(5,j,i,1)= -2.0*(alpha)*(Di(j,i,1))/(Pe*((l-X(j,i,1))**3.0))
else
FORCAS(5,j,i,1)= -2.0*(alpha)*(Di(j,i,1))/(Pe*((l-2.0)**3.0))
end if
end if
end if

end do
end do
 !$OMP END PARALLEL DO

end subroutine campo_externo

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: periodic_interactions		  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Suroutine responsible for computing properly all!
! periodic interactions. 			  !
!						  !
! This subroutine is called whenever the program  !
! identifies the presence of any periodic interac-!
! tion. These periodic interactions can be of the !
! following kinds:				  !
!						  !
! 1 - Long range hydrodynamic interactions	  !
! 2 - Long range periodic magnetic torques        !
! 3 - Long range periodic magnetic forces  	  !
!						  !
! This is probably the most expensive and complex !
! subroutine in this code.			  !
!*************************************************!

subroutine periodic_interactions

real mobilidade_self(3,3), mobilidade1(3,3), mobilidade2(3,3) 
real coeficiente1, coeficiente2,coeficiente3, coeficiente4
real coeficiente5, coeficiente6, coeficiente7, coeficiente8

! We start by computing first the sums in the real space 
! Here, nb denotes the number of physical boxes. This
! number can be different from nbr, which denotes the
! nubmer of reciprocal boxes. In simconfig.dat we 
! recommend the user to set nb=125 and nbr=27. These
! values ensure a good precision when computing the 
! average sedimentation velocity of a suspension of
! spheres in Creeping-flow

! This loop works like this: we fix a realization and
! a given particle and the make a sweep computing the
! long-range interactions between this particles and
! all the other particles in the surrounding boxes.
! We then, extend this procedure to all real particles
! in all the simultaneous numerical experiments.

do q=1,rea
do i=1,N
do s=1,nb
do j=1,N

aux_periodico= abs(ILF(s,1))+abs(ILF(s,2))+abs(ILF(s,3))

! Checking the distance between a particle "i"  and other particle "j" in the real boxes (physical and images)  

! Calculating the "r" vector (distance)

if(.not.shear)then
rij(1)= X(q,i,1)-(X(q,j,1)+(ILF(s,1)*l))
rij(2)= X(q,i,2)-(X(q,j,2)+(ILF(s,2)*l))
rij(3)= X(q,i,3)-(X(q,j,3)+(ILF(s,3)*h))
else
rij(1)= X(q,i,1)-(X(q,j,1)+(ILF(s,1)*l))
rij(2)= X(q,i,2)-(X(q,j,2)+((ILF(s,2)*l)+(k*ILF(s,3)*shearrate*dt)))
rij(3)= X(q,i,3)-(X(q,j,3)+(ILF(s,3)*h))
end if

! Normalizing the "r" vector -> r/|r|

modrij=(rij(1)**2.0 +rij(2)**2.0 + rij(3)**2.0)**0.5

rn(1)=rij(1)/modrij
rn(2)=rij(2)/modrij
rn(3)=rij(3)/modrij

d=1.0 + (9999.0*(modrij-2))/((3.0**0.5)*l*(nb**(1.0/3.0)))
diferenca_interpol1=(modrij-cof1(1,d))

!************ SUMS IN THE PHYSICAL SPACE IN THE REAL BOX (WHERE THE REAL PARTICLES ARE) AND "IMAGE" BOXES***********!


!******************************* HYDRODYNAMIC INTERACTIONS *********************************************************!

! Building the self-mobility matrix

if(aux_periodico.eq.0.0) then
if(ligaih) then
if(i.eq.j) then

do a=1,3
mobilidade_self(a,a)= 1.0-(6.0*(pi**(-0.5))*qsi)+((40.0/3.0)*(pi**(-0.5))*(qsi**3.0))
end do

end if
end if
end if


! Building the called mobility matrix 1 (lattice sum in the physical space)

if(ligaih) then
if(modrij.gt.(2.0))then

! Interpolating the pre-calculated Green functions

if(diferenca_interpol1.gt.0.0)then
coeficiente1=cof1(2,d) + ((modrij-cof1(1,d))/(cof1(1,d+1)-cof1(1,d)))*(cof1(2,d+1)-cof1(2,d))
coeficiente2=cof2(2,d) + ((modrij-cof2(1,d))/(cof2(1,d+1)-cof2(1,d)))*(cof2(2,d+1)-cof2(2,d))
else
coeficiente1=cof1(2,d-1) + ((modrij-cof1(1,d-1))/(cof1(1,d)-cof1(1,d-1)))*(cof1(2,d)-cof1(2,d-1))
coeficiente2=cof2(2,d-1) + ((modrij-cof2(1,d-1))/(cof2(1,d)-cof2(1,d-1)))*(cof2(2,d)-cof2(2,d-1))
end if

mobilidade1(1,1)=coeficiente1 + coeficiente2*(rn(1)*rn(1))
mobilidade1(1,2)=coeficiente2*(rn(1)*rn(2))
mobilidade1(1,3)=coeficiente2*(rn(1)*rn(3))

mobilidade1(2,1)=coeficiente2*(rn(2)*rn(1))
mobilidade1(2,2)=coeficiente1 + coeficiente2*(rn(2)*rn(2)) 
mobilidade1(2,3)=coeficiente2*(rn(2)*rn(3))
 
mobilidade1(3,1)=coeficiente2*(rn(3)*rn(1))
mobilidade1(3,2)=coeficiente2*(rn(3)*rn(2))
mobilidade1(3,3)=coeficiente1 +coeficiente2*(rn(3)*rn(3))

hidrodinamica_aux1(j,1)=mobilidade1(1,1)*FT(q,j,1)+mobilidade1(1,2)*FT(q,j,2)+mobilidade1(1,3)*FT(q,j,3)
hidrodinamica_aux1(j,2)=mobilidade1(2,1)*FT(q,j,1)+mobilidade1(2,2)*FT(q,j,2)+mobilidade1(2,3)*FT(q,j,3)
hidrodinamica_aux1(j,3)=mobilidade1(3,1)*FT(q,j,1)+mobilidade1(3,2)*FT(q,j,2)+mobilidade1(3,3)*FT(q,j,3)
end if ! ending the "if" of the condition (modrij.gt.2.0)
end if ! ending the "if" of the condition (ligaih)

!****************** MAGNETIC TORQUES **************************************!

if(tmagper) then
if(modrij.gt.2.0)then

! Interpolating the pre-calculated Green functions

if(diferenca_interpol1.gt.0.0)then
coeficiente4=cof4(2,d) + ((modrij-cof4(1,d))/(cof4(1,d+1)-cof4(1,d)))*(cof4(2,d+1)-cof4(2,d))
coeficiente5=cof5(2,d) + ((modrij-cof5(1,d))/(cof5(1,d+1)-cof5(1,d)))*(cof5(2,d+1)-cof5(2,d))
else
coeficiente4=cof4(2,d-1) + ((modrij-cof4(1,d-1))/(cof4(1,d)-cof4(1,d-1)))*(cof4(2,d)-cof4(2,d-1))
coeficiente5=cof5(2,d-1) + ((modrij-cof5(1,d-1))/(cof5(1,d)-cof5(1,d-1)))*(cof5(2,d)-cof5(2,d-1))
end if

! Computing periodic torques due to magnetic interactions in the real space
!eps=10.0

!if(aux_periodico.eq.0.0)then
!termo5= 4.0*pi/((2.0+eps)*(l**3.0)) ! Surface term (only present in the physical lattice)
!else
!termo5=0.0
!end if

termo5=0.0

termo2=((Di(q,j,1)*rij(1))+(Di(q,j,2)*rij(2))+(Di(q,j,3)*rij(3)))

if(gravidade)then
lambda=alpha2*6.0*pi/Pe
else
lambda=alpha2*6.0*pi
end if

termo1=(Di(q,i,3)*Di(q,j,2))-(Di(q,i,2)*Di(q,j,3))
termo3=(Di(q,i,2)*rij(3))-(Di(q,i,3)*rij(2))
auxt(j,1)=lambda*((termo1*coeficiente4)+(termo2*termo3*coeficiente5) + termo5*termo1)  

termo1=(Di(q,i,1)*Di(q,j,3))-(Di(q,i,3)*Di(q,j,1))
termo3=(Di(q,i,3)*rij(1))-(Di(q,i,1)*rij(3))
auxt(j,2)=lambda*((termo1*coeficiente4)+(termo2*termo3*coeficiente5) + termo5*termo1) 

termo1=(Di(q,i,2)*Di(q,j,1))-(Di(q,i,1)*Di(q,j,2))
termo3=(Di(q,i,1)*rij(2))-(Di(q,i,2)*rij(1))
auxt(j,3)=lambda*((termo1*coeficiente4)+(termo2*termo3*coeficiente5) + termo5*termo1) 

end if ! ending the "if" of the condition (modrij.gt.2.0)
end if ! ending the "if" of the condition (tmagper)

!*********************************** MAGNETIC FORCES *****************************************!
if(fmagper) then
if(modrij.gt.(2.0))then

! Interpolating the pre-calculated Green functions

if(diferenca_interpol1.gt.0.0)then
coeficiente6=cof6(2,d) + ((modrij-cof6(1,d))/(cof6(1,d+1)-cof6(1,d)))*(cof6(2,d+1)-cof6(2,d))
else
coeficiente6=cof6(2,d-1) + ((modrij-cof6(1,d-1))/(cof6(1,d)-cof6(1,d-1)))*(cof6(2,d)-cof6(2,d-1))
end if

! Computing periodic forces due to magnetic interactions in the real space

if(gravidade)then
lambda=alpha2*8.0*pi/Pe
else
lambda=alpha2*8.0*pi
end if

! (di.dj)rij
termo1=((Di(q,i,1)*Di(q,j,1))+(Di(q,i,2)*Di(q,j,2))+(Di(q,i,3)*Di(q,j,3)))*rij(1)
! (dj.rij)di
termo3=((Di(q,j,1)*rij(1))+(Di(q,j,2)*rij(2))+(Di(q,j,3)*rij(3)))*Di(q,i,1)
! (di.rij)dj
termo2=((Di(q,i,1)*rij(1))+(Di(q,i,2)*rij(2))+(Di(q,i,3)*rij(3)))*Di(q,j,1)
!((di.rij)(dij.rij))rij
termo4=(((Di(1,i,1)*rij(1))+(Di(1,i,2)*rij(2))   &
+(Di(1,i,3)*rij(3)))*((Di(q,j,1)*rij(1))   &
+(Di(q,j,2)*rij(2))+(Di(q,j,3)*rij(3))))*rij(1)

auxf(j,1)=lambda*(((termo1+termo2+termo3)*coeficiente5) - termo4*coeficiente6)

! (di.dj)rij
termo1=((Di(q,i,1)*Di(q,j,1))+(Di(q,i,2)*Di(q,j,2))+(Di(q,i,3)*Di(q,j,3)))*rij(2)
! (dj.rij)di
termo3=((Di(q,j,1)*rij(1))+(Di(q,j,2)*rij(2))+(Di(q,j,3)*rij(3)))*Di(q,i,2)
! (di.rij)dj
termo2=((Di(q,i,1)*rij(1))+(Di(q,i,2)*rij(2))+(Di(q,i,3)*rij(3)))*Di(q,j,2)
!((di.rij)(dij.rij))rij
termo4=(((Di(1,i,1)*rij(1))+(Di(1,i,2)*rij(2))   &
+(Di(1,i,3)*rij(3)))*((Di(q,j,1)*rij(1))   &
+(Di(q,j,2)*rij(2))+(Di(q,j,3)*rij(3))))*rij(2)

auxf(j,2)=lambda*(((termo1+termo2+termo3)*coeficiente5) - termo4*coeficiente6)

! (di.dj)rij
termo1=((Di(q,i,1)*Di(q,j,1))+(Di(q,i,2)*Di(q,j,2))+(Di(q,i,3)*Di(q,j,3)))*rij(3)
! (dj.rij)di
termo3=((Di(q,j,1)*rij(1))+(Di(q,j,2)*rij(2))+(Di(q,j,3)*rij(3)))*Di(q,i,3)
! (di.rij)dj
termo2=((Di(q,i,1)*rij(1))+(Di(q,i,2)*rij(2))+(Di(q,i,3)*rij(3)))*Di(q,j,3)
!((di.rij)(dij.rij))rij
termo4=(((Di(1,i,1)*rij(1))+(Di(1,i,2)*rij(2))   &
+(Di(1,i,3)*rij(3)))*((Di(q,j,1)*rij(1))   &
+(Di(q,j,2)*rij(2))+(Di(q,j,3)*rij(3))))*rij(3)

auxf(j,3)=lambda*(((termo1+termo2+termo3)*coeficiente5) - termo4*coeficiente6)
end if ! ending the "if" of the condition  (modrij.gt.2.0) 
end if ! ending the "if" of the condition (fmagper)

end do
!************************** SUMING EVERYTHING *********************************!

! Computing the sum of the previous vectors and matrices

if(ligaih) then
hidro1(s,1)=sum(hidrodinamica_aux1(:,1))
hidro1(s,2)=sum(hidrodinamica_aux1(:,2))
hidro1(s,3)=sum(hidrodinamica_aux1(:,3))
end if

if(fmagper) then
forcareal(s,1)=sum(auxf(:,1))
forcareal(s,2)=sum(auxf(:,2))
forcareal(s,3)=sum(auxf(:,3))
end if

if(tmagper) then
torquereal(s,1)=sum(auxt(:,1))
torquereal(s,2)=sum(auxt(:,2))
torquereal(s,3)=sum(auxt(:,3))
auxt=0.0
end if

end do

! Computing the contributions of the real space interactions in the velocities of each particle

if(ligaih) then
U(q,i,1)=mobilidade_self(1,1)*FT(q,i,1) + sum(hidro1(:,1))
U(q,i,2)=mobilidade_self(2,2)*FT(q,i,2) + sum(hidro1(:,2))
U(q,i,3)=mobilidade_self(3,3)*FT(q,i,3) + sum(hidro1(:,3))
end if

! Computing the contributions of the real space interactions in the magnetic forces and torques acting on each particle

if(fmagper) then
FORCAS(4,q,i,1)=sum(forcareal(:,1))
FORCAS(4,q,i,2)=sum(forcareal(:,2))
FORCAS(4,q,i,3)=sum(forcareal(:,3))
end if

if(tmagper) then
TORQUES(1,q,i,1)=sum(torquereal(:,1))
TORQUES(1,q,i,2)=sum(torquereal(:,2))
TORQUES(1,q,i,3)=sum(torquereal(:,3))
end if

contribuicao_self(q,i)=mobilidade_self(3,3)*FT(q,i,3)
contribuicao_fisico(q,i)=sum(hidro1(:,3))

end do
end do

!************************************* SUMS IN THE RECIPROCAL SPACE ***********************************!
do q=1,rea
do i=1,N
do s=1,nbr
do j=1,N

aux_periodico= abs(ILR(s,1))+abs(ILR(s,2))+abs(ILR(s,3))

if(aux_periodico.ne.0.0) then

rij(1)= X(q,i,1)-X(q,j,1)
rij(2)= X(q,i,2)-X(q,j,2)
rij(3)= X(q,i,3)-X(q,j,3)

! Calculating the normalized "r" vector = r/|r|

modrij=(rij(1)**2.0 +rij(2)**2.0 + rij(3)**2.0)**0.5

rn(1)=rij(1)/modrij
rn(2)=rij(2)/modrij
rn(3)=rij(3)/modrij

if(modrij.ne.0) then

! Computing the wave number vector based on the lattice index 

konda(1)=ILR(s,1)*2.0*pi/l
konda(2)=ILR(s,2)*2.0*pi/l
konda(3)=ILR(s,3)*2.0*pi/h

! Calculating the wave number module

modk=((konda(1)**2.0)+(konda(2)**2.0)+(konda(3)**2.0))**0.5

! Normalized wave number vector

knormal(1)=konda(1)/modk
knormal(2)=konda(2)/modk
knormal(3)=konda(3)/modk


!********************** HYDRODYNAMIC INTERACTIONS *********************************************************!

if(ligaih) then

! Interpolating the pre-calculated Green function

 call interpola_reciproco(nbr,cof3,cof3,coeficiente3,modk,l)

! Calculating the mobility matrix for the reciprocal space contribution

mobilidade2(1,1)=coeficiente3*(1.0 - knormal(1)*knormal(1))
mobilidade2(1,2)=coeficiente3*(-knormal(1)*knormal(2))
mobilidade2(1,3)=coeficiente3*(-knormal(1)*knormal(3))
mobilidade2(2,1)=coeficiente3*(-knormal(2)*knormal(1))
mobilidade2(2,2)=coeficiente3*(1.0 - knormal(2)*knormal(2))
mobilidade2(2,3)=coeficiente3*(- knormal(2)*knormal(3))
mobilidade2(3,1)=coeficiente3*(-knormal(3)*knormal(1))
mobilidade2(3,2)=coeficiente3*(-knormal(3)*knormal(2))
mobilidade2(3,3)=coeficiente3*(1.0 - knormal(3)*knormal(3))
kr=(knormal(1)*rn(1))+(knormal(2)*rn(2))+(knormal(3)*rn(3))

kr2=cos((konda(1)*rij(1))+(konda(2)*rij(2))+(konda(3)*rij(3)))

hidrodinamica_aux2(j,1)=kr*(mobilidade2(1,1)*FT(q,j,1)+mobilidade2(1,2)*FT(q,j,2)+mobilidade2(1,3)*FT(q,j,3))
hidrodinamica_aux2(j,2)=kr*(mobilidade2(2,1)*FT(q,j,1)+mobilidade2(2,2)*FT(q,j,2)+mobilidade2(2,3)*FT(q,j,3))
hidrodinamica_aux2(j,3)=kr*(mobilidade2(3,1)*FT(q,j,1)+mobilidade2(3,2)*FT(q,j,2)+mobilidade2(3,3)*FT(q,j,3))
end if !closing the "if" of the condition (ligaih)

!**************************** MAGNETIC TORQUES ***********************************************!

if(tmagper) then

kr2=1.0*((konda(1)*rij(1))+(konda(2)*rij(2))+(konda(3)*rij(3)))

! Interpolating the pre-calculated Green function

 call interpola_reciproco(nbr,cof3,cof7,coeficiente7,modk,l)

! Computing the reciprocal space contribution
termo2=((Di(q,j,1)*konda(1))+(Di(q,j,2)*konda(2))+(Di(q,j,3)*konda(3)))

!if(Pe.eq.0.0)then
!lambda=alpha2*8.0
!else
!lambda=alpha2*8.0/Per
!end if

if(gravidade)then
lambda=alpha2*6.0*pi/Pe
else
lambda=alpha2*6.0*pi
end if


termo1=(Di(q,i,2)*konda(3))-(Di(q,i,3)*konda(2))
auxt(j,1)=lambda*(coeficiente7*(termo1*termo2)*cos(kr2))   
termo1=(Di(q,i,3)*konda(1))-(Di(q,i,1)*konda(3))
auxt(j,2)=lambda*(coeficiente7*(termo1*termo2)*cos(kr2))   
termo1=(Di(q,i,1)*konda(2))-(Di(q,i,2)*konda(1))
auxt(j,3)=lambda*(coeficiente7*(termo1*termo2)*cos(kr2))  
end if !closing the "if" of the condition (tmagper)

!********************************* MAGNETIC FORCES **********************************************!
if(fmagper) then

! Interpolating the pre-calculated Green function

 call interpola_reciproco(nbr,cof3,cof8,coeficiente8,modk,l)

! Computing the reciprocal space contribution

kr=(konda(1)*rij(1))+(konda(2)*rij(2))+(konda(3)*rij(3))
termo1=((Di(q,i,1)*konda(1))+(Di(q,i,2)*konda(2))+(Di(q,i,3)*konda(3)))
termo2=((Di(q,j,1)*konda(1))+(Di(q,j,2)*konda(2))+(Di(q,j,3)*konda(3)))

if(Pe.eq.0.0)then
lambda=alpha2*8.0
else
lambda=alpha2*8.0/Pe
end if
auxf(j,1)=lambda*(coeficiente8*termo1*termo2*sin(2.0*pi*kr/l))*konda(1)
auxf(j,2)=lambda*(coeficiente8*termo1*termo2*sin(2.0*pi*kr/l))*konda(2) 
auxf(j,3)=lambda*(coeficiente8*termo1*termo2*sin(2.0*pi*kr/l))*konda(3)   
end if !closing the "if" of the condition (fmagper)

end if !closing the "if" of the condition (modrij.ne.0)
end if !closing the "if" of the condition (auxperiodico.ne.0)

end do

! Making the final sum for the physical and reciprocal contribution for all possible interactions

if(ligaih) then
hidro2(s,1)=sum(hidrodinamica_aux2(:,1))
hidro2(s,2)=sum(hidrodinamica_aux2(:,2))
hidro2(s,3)=sum(hidrodinamica_aux2(:,3))
end if

if(fmagper) then
forcareciproca(s,1)=sum(auxf(:,1))
forcareciproca(s,2)=sum(auxf(:,2))
forcareciproca(s,3)=sum(auxf(:,3))
end if

if(tmagper) then
torquereciproco(s,1)=sum(auxt(:,1))
torquereciproco(s,2)=sum(auxt(:,2))
torquereciproco(s,3)=sum(auxt(:,3))
end if

end do

if(ligaih)then
! Computing the velocities of the particles, now with the reciprocal space sum contribution
U(q,i,1)= U(q,i,1) + sum(hidro2(:,1))
U(q,i,2)= U(q,i,2) + sum(hidro2(:,2))
U(q,i,3)= U(q,i,3) + sum(hidro2(:,3))
contribuicao_reciproco(q,i)=sum(hidro2(:,3))
end if

if(fmagper)then
! Computing the magnetic forces acting on the particles, now with the reciprocal space sum contribution
FORCAS(4,q,i,1)=FORCAS(4,q,i,1)+sum(forcareciproca(:,1))
FORCAS(4,q,i,2)=FORCAS(4,q,i,2)+sum(forcareciproca(:,2))
FORCAS(4,q,i,3)=FORCAS(4,q,i,3)+sum(forcareciproca(:,3))
end if

if(tmagper)then
! Computing the magnetic torques acting on the particles, now with the reciprocal space sum contribution
TORQUES(1,q,i,1)=TORQUES(1,q,i,1)+sum(torquereciproco(:,1))
TORQUES(1,q,i,2)=TORQUES(1,q,i,2)+sum(torquereciproco(:,2))
TORQUES(1,q,i,3)=TORQUES(1,q,i,3)+sum(torquereciproco(:,3))
end if

end do
end do

end subroutine periodic_interactions

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: interpola_reciproco		  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine resposible for interpolating the va- !
! lues of the pre-calculated the Green functions  !
! used in the calculation of periodic interactions!
! due to hydrodynamic forces, magnetic forces and !
! torques (only in the reciprocal space)          !
!*************************************************!

subroutine interpola_reciproco(a,b,c,d,e,f)

! a = number of reciprocal lattices
! b = cof3
! c = cofx
! d = coeficientex
! e = norm of k  (modk)
! f = l (lenght of the box)

integer a
real b(2,10000), c(2,10000), d
real g, modk, pi, l, diferenca_interpol2

pi=acos(-1.0)

if(a**(1.0/3.0).eq.5.0)then
g=1+9999*((e-(2*pi/f))/(((4*pi/f)*(3**0.5))-(2*pi/f)))
end if
if(nbr**(1.0/3.0).eq.3.0)then
g=1+9999*((e-(2*pi/f))/(((2*pi/f)*(3**0.5))-(2*pi/f)))
end if
diferenca_interpol2=(e-b(1,g))
if(diferenca_interpol2.gt.0.0)then
d=c(2,g) + ((e-c(1,g))/(c(1,g+1)-c(1,g)))*(c(2,g+1)-c(2,g))
else
d=c(2,g-1) + ((e-c(1,g-1))/(c(1,g)-c(1,g-1)))*(c(2,g)-c(2,g-1))
end if

end subroutine interpola_reciproco

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: resvel				  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine resposible for calculating the parti-!
! cles velocity using a 4th order Runge-Kutta.    !
!*************************************************!

! IT IS IMPORTANT TO NOTICE THAT THIS SUBROUTINE 
! DOES NOT APPLY FOR A MOBILITY PROBLEM. IT ONLY 
! MAKES SENSE WHEN WE CONSIDER PARTICLE INERTIA 
! (RESISTANCE FORMULATION)

subroutine resvel(a,b,c,d)
real a                      ! velocity component
real b			    ! time step
real c                      ! Stokes number (null for zero inertia)
real d                      ! sum of forces in a given direction
real k1,k2,k3,k4            ! internal variables

k1=b*(-a+d)/c
k2=b*((-a-0.5*k1)+d)/c
k3=b*((-a-0.5*k2)+d)/c
k4=b*((-a-k3)+d)/c

a=a+(1.0/6.0)*(k1+2.0*k2+2.0*k3+k4)

end subroutine resvel

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: respos				  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine resposible for calculating the parti-!
! cles position.				  !
!*************************************************!

subroutine respos(a,b,c)
real a                      ! position
real b			    ! time-step
real c                      ! velocity component
a=a+b*c
end subroutine respos

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: writting_files			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine resposible for writting output files !
! containing the positions, velocities and orien- !
! tations of all particles in all realiations     !
!*************************************************!

subroutine writting_files(k,k_real)
integer k, teste2
real k_real, teste1

509 FORMAT(F30.4,F30.4,F30.4,F30.4)
666 FORMAT(F30.4,F30.4,F30.4,F30.4,F30.4,F30.4,F30.4)
if(continua)then
teste1=k/n3
teste2=k/n2

if(teste1.eq.teste2) then

if(posicao)then
do j=1,rea
if(.not.ovito) then
write(j,'(A12,I6,A1)') 'zone t="',k,'"'
else
write(j,*) N
write(j,*) 'Time =',k*dt
end if
do i=1,N
write(j,666)X(j,i,1),X(j,i,2),   &
X(j,i,3),Di(j,i,1),Di(j,i,2),Di(j,i,3),DIAM(j,i)
end do
end do
end if

write(100*rea,*) X(1,1,1),X(1,1,2),X(1,1,3), Di(1,1,1), Di(1,1,2), Di(1,1,3), k*dt
write(*,*) 'SIMULATION PROGRESS:', (k_real/aux_real)*100, '%'

if(agregado_inicial) then
write(666,*) k*dt, sum(U(:,:,1))/(N*rea),sum(U(:,:,2))/(N*rea),sum(U(:,:,3))/(N*rea)
end if

if(velocidade)then
do j=1,rea
write(rea+j,'(A12,I6,A1)') 'zone t="',k,'"'
do i=1,N
write(rea+j,509)U(j,i,1),U(j,i,2),   &
U(j,i,3)
end do
end do
end if

end if

else

if(k.eq.1) then
if(posicao)then
do j=1,rea
if(.not.ovito) then
write(j,'(A12,I6,A1)') 'zone t="',k,'"'
else
write(j,*) N
write(j,*) 'Time =',k*dt
end if
do i=1,N
if(gravadipolo)then
write(j,666)X(j,i,1),X(j,i,2),   &
X(j,i,3),Di(j,i,1),Di(j,i,2),Di(j,i,3),DIAM(j,i)
else
write(j,509)X(j,i,1),X(j,i,2),X(j,i,3),DIAM(j,i)
end if
end do
end do
end if

write(100*rea,*) X(1,1,1),X(1,1,2),X(1,1,3), Di(1,1,1), Di(1,1,2), Di(1,1,3), k*dt
write(*,*) 'SIMULATION PROGRESS:', (k_real/aux_real)*100, '%'

if(agregado_inicial) then
write(666,*) k*dt, sum(U(:,:,1))/(N*rea),sum(U(:,:,2))/(N*rea),sum(U(:,:,3))/(N*rea)
end if

if(velocidade)then
do j=1,rea
write(rea+j,'(A12,I6,A1)') 'zone t="',k,'"'
do i=1,N
write(rea+j,509)U(j,i,1),U(j,i,2),   &
U(j,i,3)
end do
end do
end if

end if

teste1=k/n3
teste2=k/n2

if(k.ne.1) then

if(teste1.eq.teste2) then
if(posicao)then
do j=1,rea
if(.not.ovito) then
write(j,*) 'zone t="',k,'"'
else
write(j,*) N
write(j,*) 'Time =',k*dt
end if
do i=1,N
if(gravadipolo)then
write(j,666)X(j,i,1),X(j,i,2),   &
X(j,i,3),Di(j,i,1),Di(j,i,2),Di(j,i,3),DIAM(j,i)
else
write(j,509)X(j,i,1),X(j,i,2),X(j,i,3),DIAM(j,i)
end if
end do
end do
end if

write(100*rea,*) X(1,1,1),X(1,1,2),X(1,1,3), Di(1,1,1), Di(1,1,2), Di(1,1,3), k*dt

write(*,*) 'SIMULATION PROGRESS:', (k_real/aux_real)*100, '%'  

if(agregado_inicial) then
write(666,*) k*dt, sum(U(:,:,1))/(N*rea),sum(U(:,:,2))/(N*rea),sum(U(:,:,3))/(N*rea)
end if

if(velocidade)then
do j=1,rea
write(rea+j,'(A12,I6,A1)') 'zone t="',k,'"'
do i=1,N
write(rea+j,509)U(j,i,1),U(j,i,2),   &
U(j,i,3)
end do
end do
end if
end if
end if
end if

end subroutine writting_files

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: torque_magnetico			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine resposible for computing magnetic    !
! torques due to dipolar interactions when this   !
! computation is performed in a non-periodic way  !
!*************************************************!

subroutine torque_magnetico
integer auxiliary

if(gravidade)then
lambda=alpha2*24.0
else
if(shear) then
lambda=alpha2*3.0/Pe
else
lambda=alpha2*9.0/2.0
end if
end if

! These subroutine may consider a magnetic fluidized 
! bed mixing magnetic and non-magnetic particles, 
! the logical variable "mistura" activates this po-
! ssibility. If mistura = TRUE then we have to set 
! zero dipoles for a certain percentage "percentual" 
! of the particles.

if(mistura)then
 auxiliary=(percentual*N)+1
else
 auxiliary=1
end if

do j=1,rea
do i=auxiliary,N
do q=auxiliary,N
if(i.ne.q) then
! Calculating the distance between all the pairs of particles
r=(((X(j,i,1)-X(j,q,1))**2.0)+((X(j,i,2)-X(j,q,2))**2.0)+((X(j,i,3)-X(j,q,3))**2.0))**0.5
if(r.le.2.0) then
!if(r.le.2.0) then
aux1(j,q)=0.0
aux2(j,q)=0.0
aux3(j,q)=0.0
else

! Calculating the vector that connects a particle i to a particle j
rij(1)=X(j,i,1)-X(j,q,1)
rij(2)=X(j,i,2)-X(j,q,2)
rij(3)=X(j,i,3)-X(j,q,3)
! Normalizing this vector
modrij=((rij(1)**2.0)+(rij(2)**2.0)+(rij(3)**2.0))**0.5
rij(1)=rij(1)/modrij
rij(2)=rij(2)/modrij
rij(3)=rij(3)/modrij

termo2=((Di(j,q,1)*rij(1))+(Di(j,q,2)*rij(2))+(Di(j,q,3)*rij(3)))

termo1=((Di(j,i,2)*Di(j,q,3))-(Di(j,i,3)*Di(j,q,2)))
termo3=((Di(j,i,2)*rij(3))-(Di(j,i,3)*rij(2)))

aux1(j,q)=(lambda/(1.0*(r**3.0)))*(((-1.0/3.0)*termo1)+(termo2*termo3))

termo1=((Di(j,i,3)*Di(j,q,1))-(Di(j,i,1)*Di(j,q,3)))
termo3=((Di(j,i,3)*rij(1))-(Di(j,i,1)*rij(3)))

aux2(j,q)=(lambda/(1.0*(r**3.0)))*(((-1.0/3.0)*termo1)+(termo2*termo3))

termo1=((Di(j,i,1)*Di(j,q,2))-(Di(j,i,2)*Di(j,q,1)))
termo3=((Di(j,i,1)*rij(2))-(Di(j,i,2)*rij(1)))

aux3(j,q)=(lambda/(1.0*(r**3.0)))*(((-1.0/3.0)*termo1)+(termo2*termo3))

end if
end if
end do
TORQUES(1,j,i,1)=sum(aux1(j,:))
TORQUES(1,j,i,2)=sum(aux2(j,:))
TORQUES(1,j,i,3)=sum(aux3(j,:))

aux1=0.0
aux2=0.0
aux3=0.0
r=0.0
modrij=0.0

end do
end do

end subroutine torque_magnetico

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: torque_externo			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine resposible for computing magnetic    !
! torques due to an external field		  !
!*************************************************!

subroutine torque_externo(alfa)
integer auxiliary
real alfa

! These subroutine may consider a magnetic fluidized 
! bed mixing magnetic and non-magnetic particles, 
! the logical variable "mistura" activates this po-
! ssibility. If mistura = TRUE then we have to set 
! zero dipoles for a certain percentage "percentual" 
! of the particles.

if(mistura)then
 auxiliary=(percentual*N)+1
else
 auxiliary=1
end if


 !$OMP PARALLEL DO
do j=1,rea
do i=auxiliary,N

! posicao_campo = 1 -> Applied field on the lower wall

if(posicao_campo.eq.1)then
if(gravidade) then
TORQUES(2,j,i,1)=-(3.0*alfa*Di(j,i,2)/4.0*Pe)
TORQUES(2,j,i,2)=(3.0*alfa*Di(j,i,1)/4.0*Pe)
TORQUES(2,j,i,3)=0.0
else
if(shear) then
TORQUES(2,j,i,1)=-(3.0*alfa*Di(j,i,2)/4.0*Pe)
TORQUES(2,j,i,2)=(3.0*alfa*Di(j,i,1)/4.0*Pe)
TORQUES(2,j,i,3)=0.0
else
TORQUES(2,j,i,1)=-(3.0*alfa*Di(j,i,2)/4.0)
TORQUES(2,j,i,2)=(3.0*alfa*Di(j,i,1)/4.0)
TORQUES(2,j,i,3)=0.0
end if
end if
end if

! posicao_campo = 2 ->  Applied field on the upper wall

if(posicao_campo.eq.2)then
if(gravidade) then
TORQUES(2,j,i,1)=(3.0*alfa*Di(j,i,2)/4.0*Pe)
TORQUES(2,j,i,2)=-(3.0*alfa*Di(j,i,1)/4.0*Pe)
TORQUES(2,j,i,3)=0.0
else
if(shear) then
TORQUES(2,j,i,1)=(3.0*alfa*Di(j,i,2)/4.0*Pe)
TORQUES(2,j,i,2)=-(3.0*alfa*Di(j,i,1)/4.0*Pe)
TORQUES(2,j,i,3)=0.0
else
TORQUES(2,j,i,1)=(3.0*alfa*Di(j,i,2)/4.0)
TORQUES(2,j,i,2)=-(3.0*alfa*Di(j,i,1)/4.0)
TORQUES(2,j,i,3)=0.0
end if
end if
end if

! posicao_campo = 3 ->  Applied field on the right side

if(posicao_campo.eq.3) then
if(gravidade) then
TORQUES(2,j,i,1)=0.0
TORQUES(2,j,i,2)=-(3.0*alfa*Di(j,i,3)/4.0*Pe)
TORQUES(2,j,i,3)=(3.0*alfa*Di(j,i,2)/4.0*Pe)
else
if(shear) then
TORQUES(2,j,i,1)=0.0
TORQUES(2,j,i,2)=-(3.0*alfa*Di(j,i,3)/4.0*Pe)
TORQUES(2,j,i,3)=(3.0*alfa*Di(j,i,2)/4.0*Pe)
else
TORQUES(2,j,i,1)=0.0
TORQUES(2,j,i,2)=-(3.0*alfa*Di(j,i,3)/4.0)
TORQUES(2,j,i,3)=(3.0*alfa*Di(j,i,2)/4.0)
end if
end if
end if

! posicao_campo = 4 ->  Applied field on the left side

if(posicao_campo.eq.4) then
if(gravidade) then
TORQUES(2,j,i,1)=0.0
TORQUES(2,j,i,2)=(3.0*alfa*Di(j,i,3)/4.0*Pe)
TORQUES(2,j,i,3)=-(3.0*alfa*Di(j,i,2)/4.0*Pe)
else
if(shear) then
TORQUES(2,j,i,1)=0.0
TORQUES(2,j,i,2)=(3.0*alfa*Di(j,i,3)/4.0*Pe)
TORQUES(2,j,i,3)=-(3.0*alfa*Di(j,i,2)/4.0*Pe)
else
TORQUES(2,j,i,1)=0.0
TORQUES(2,j,i,2)=(3.0*alfa*Di(j,i,3)/4.0)
TORQUES(2,j,i,3)=-(3.0*alfa*Di(j,i,2)/4.0)
end if
end if
end if

end do
end do
 !$OMP END PARALLEL DO

end subroutine torque_externo

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: rotating_field			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine resposible for computing magnetic    !
! torques due to an external rotating field	  !
!*************************************************!

subroutine rotating_field(alfa, wt)
integer auxiliary
real alfa, wt

! These subroutine may consider a magnetic fluidized 
! bed mixing magnetic and non-magnetic particles, 
! the logical variable "mistura" activates this po-
! ssibility. If mistura = TRUE then we have to set 
! zero dipoles for a certain percentage "percentual" 
! of the particles.

if(mistura)then
 auxiliary=(percentual*N)+1
else
 auxiliary=1
end if

 !$OMP PARALLEL DO
do j=1,rea
do i=auxiliary,N

! Rotating field on the upper side

if(gravidade) then
TORQUES(2,j,i,1)=(3.0*alfa/4.0*Pe)*(Di(j,i,2)*cos(wt) - (Di(j,i,3)*sin(wt)))
TORQUES(2,j,i,2)=(3.0*alfa/4.0*Pe)*(-(Di(j,i,1)*cos(wt)))
TORQUES(2,j,i,3)=(3.0*alfa/4.0*Pe)*((Di(j,i,1)*sin(wt)))
else
if(shear) then
TORQUES(2,j,i,1)=(3.0*alfa/4.0*Pe)*(Di(j,i,2)*cos(wt) - (Di(j,i,3)*sin(wt)))
TORQUES(2,j,i,2)=(3.0*alfa/4.0*Pe)*(-(Di(j,i,1)*cos(wt)))
TORQUES(2,j,i,3)=(3.0*alfa/4.0*Pe)*((Di(j,i,1)*sin(wt)))
else
TORQUES(2,j,i,1)=(3.0*alfa/4.0)*(Di(j,i,2)*cos(wt) - (Di(j,i,3)*sin(wt)))
TORQUES(2,j,i,2)=(3.0*alfa/4.0)*(-(Di(j,i,1)*cos(wt)))
TORQUES(2,j,i,3)=(3.0*alfa/4.0)*((Di(j,i,1)*sin(wt)))
end if
end if
end do
end do
 !$OMP END PARALLEL DO

end subroutine rotating_field


!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: resvel_sem_inercia			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine resposible for calculating the velo- !
! city of the particles in a non inertial context !
!*************************************************!

subroutine resvel_sem_inercia(a,b,c)
real a   ! velocity component
real b	 ! sum of forces in a given direction
real c   ! 0 or 1, 0 = no gravity, 1 = gravity

a=-c+b
end subroutine resvel_sem_inercia

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: resomega				  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine resposible for calculating the angu- !
! lar velocity of the particles 		  !
!*************************************************!

subroutine resomega(a,b,c,d)
real a  ! rotational velocitty in a given direction
real b	! time step
real c  ! Rotational Stokes number
real d  ! Sum of torques in a given direction
real k1,k2,k3,k4  ! Internal variables

k1=b*(-a+d)/c
k2=b*((-a-0.5*k1)+d)/c
k3=b*((-a-0.5*k2)+d)/c
k4=b*((-a-k3)+d)/c

a=a+(1.0/6.0)*(k1+2.0*k2+2.0*k3+k4)

end subroutine resomega

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: evoldip				  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine resposible for evolving the dipole	  !
! orientation for each particle			  !
!*************************************************!

subroutine evoldip(a,b,c,d,e,f)
real a	! Dipole moment orientation to be calculated
real b	! Particle dipole moment in the j direction
real c	! Particle dipole moment in the k direction
real d  ! Angular velocity in the j direction
real e  ! Angular velocity in the k direction
real f  ! Time-step

! d(di)/dt = omega_i x di for particle i

a=a+(d*c -e*b)*f

end subroutine evoldip

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: randomica				  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine resposible for calculating all random!
! numbers used in all simulations		  !
!*************************************************!

subroutine randomica(a,b,c,n,d)
real a,b       ! a,b = random number range
integer n, m 
real c(n)      ! c = generated random sequence
integer d,i,e
integer f(8)
integer, allocatable :: seed(:)

 call random_seed(size = m)
allocate (seed(m))

 CALL DATE_AND_TIME(values=f)
 CALL SYSTEM_CLOCK(count=e)

do i = 1,m
seed(i) =  47*d + f(8)*i*d*12000 + e*(3*d+i)
end do

 call random_seed(put = seed)

 call random_number(c)

 c = a+(b-a)*c

end subroutine randomica

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: fator_estrutura			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine resposible for calculating the suspen!
! sion structure factor: NEEDS MORE TESTS	  !
!*************************************************!

subroutine fator_estrutura(X,N,l,h,dt,rea)
integer N,i,j, q, contador,it,r,p, loop,rea,mediai, nk
real X(rea,N,3),dt, mediar
real l,h, rij1,rij2,rij3
real lenghx, lenghz, pi, DKR
complex*16 contagem, IM
integer a,b,c, nn

integer, allocatable :: kx(:)
integer, allocatable :: ky(:)
integer, allocatable :: kz(:)
real, allocatable :: S(:,:)
real, allocatable :: Saux(:)
real, allocatable :: Sk(:)
real, allocatable :: modk(:)
real, allocatable :: Sf(:,:)
complex*16, allocatable :: SI(:,:)
complex*16, allocatable :: SII(:,:)
integer, allocatable :: auxper(:)
integer, allocatable :: auxper2(:)
integer, allocatable :: auxper3(:)

! Defining the number of wave vectors directions

nn=55

nk=nn**3.D0

allocate(kx(nk))
allocate(ky(nk))
allocate(kz(nk))
allocate(S(nk,rea))
allocate(SI(nk,rea))
allocate(SII(nk,rea))
allocate(Saux(nk))
allocate(Sk(nk))
allocate(modk(nk))
allocate(auxper(nk))
allocate(auxper2(nk))
allocate(auxper3(nk))

! Creating the nk vectors for several possible directions of the wave number vector

auxper(1)=1
auxper2(1)=0
auxper3(1)=0

do i=2,nn
auxper(i)=i
auxper2(i)=nn*(i-1)
auxper3(i)=(nn*nn)*(i-1)
end do

do a=1,nn
do b=1,nn
do c=1,nn
q=auxper3(a)+ auxper2(b) +auxper(c)
kx(q)=a- (nn -(nn-1))
ky(q)=b- (nn -(nn-1))
kz(q)=c- (nn- (nn-1))
end do
end do
end do

! Defining pi, I and other parameters

pi = 3.141592653589793D0
IM = CMPLX(0.D0,1.D0)
SI = CMPLX(0.D0,0.D0)
DKR=2.D0*pi/l

! Calculating the Structure fator for several wave number vectors

do it=1,rea
SI = CMPLX(0.D0,0.D0)
do q=1,nk
do i=1,N
rij1=X(it,i,1)
rij2=X(it,i,2)
rij3=X(it,i,3)
SI(q,it) = SI(q,it) + EXP(IM*DKR*(kx(q)*X(it,i,1) + ky(q)*X(it,i,2) + kz(q)*X(it,i,3)))
end do
SI(q,it) = SI(q,it)*CONJG(SI(q,it))
S(q,it) = (1.D0/N)*REAL(SI(q,it))
SI(q,it) = CMPLX(0.D0,0.D0)
end do
end do

!************************************* Ploting the curve S(k) x k for isotropic configurations ****************************************!

! We already calculated the structure factor for several directions. The problem is that many "direction vectors" can produce the same
! wave number magnitude. Lets organize this...

! Calculating the wave numbers corresponding for each direction and making an essemble avarage of the structure factor over
! all the realizations.

do q=1,nk
modk(q) = SQRT(DKR*(kx(q)**2.D0 + ky(q)**2.D0 + kz(q)**2.D0))
Sk(q)=(1.D0/rea)*sum(S(q,:))
end do

! Organizing the table S x k in ascending order

 call ordem_crescente(modk,Sk,nk)

! Now we need to take an average over different directions that originate the same wave number "k"

! Let's first count how many wave numbers we have

i=1

do q=2,nk
if(modk(q)-modk(q-1).ne.0.D0) then
i=i+1
end if
end do

! We have "i" values of wave number. Let's allocate this final array with "i" terms

allocate(Sf(i,2))

! Setting the initial terms

Sf(1,1)=modk(1)
Sf(1,2)=Sk(1)

j=1
a=1

! Taking an average over different directions that originate the same wave number "k"

do q=3,nk
if(modk(q)-modk(q-1).eq.0.D0) then
j=j+1
else
Sf(a+1,2)=sum(Sk(a+1:a+j))/j
Sf(a+1,1)=sum(modk(a+1:a+j))/j

j=1
a=a+1
end if
end do

! Writing in a file
write(6*rea,*) 'variables="k","S"'
write(6*rea,*) 'zone t="Sxk(isotropic)"'
do q=2,i-1
if(Sf(q,1).ne.0.D0)then
write(6*rea,'(F20.12,F20.12)') Sf(q,1),Sf(q,2)
end if
end do

deallocate(Sf)

!****************************** Ploting the curve S(k) x k for directions kx, ky and kz (non-isotropic) *******************************!


!************************************************************* kx *********************************************************************!

do q=1,nk
modk(q) = DKR*kx(q)
end do

! Organizing the table S x k_x in ascending order

 call ordem_crescente(modk,Sk,nk)

! Let's first count how many wave numbers in the x directions we have

i=1

do q=2,nk
if(modk(q)-modk(q-1).ne.0.D0) then
i=i+1
end if
end do

! We have "i" values of wave number. Let's allocate this final array with "i" terms

allocate(Sf(i,2))

! Setting the initial terms

Sf(1,1)=modk(1)
Sf(1,2)=Sk(1)

j=1
a=1

! Taking an average over different kx that originate different values of S (due to the other directions)

do q=3,nk
if(modk(q)-modk(q-1).eq.0.D0) then
Sf(a+1,1)=modk(q)
j=j+1
else
Sf(a+1,2)=sum(Sk(a+1:a+j))/j

j=1
a=a+1
end if
end do

! Writing in a file
write(6*rea,*) 'zone t="Sxk_x"'
do q=2,i-1
if(Sf(q,1).ne.0.D0)then
write(6*rea,'(F20.12,F20.12)') Sf(q,1),Sf(q,2)
end if
end do

deallocate(Sf)

!************************************************************* ky *********************************************************************!

do q=1,nk
modk(q) = DKR*ky(q)
end do

! Organizing the table S x k_x in ascending order

 call ordem_crescente(modk,Sk,nk)

! Let's first count how many wave numbers in the x directions we have

i=1

do q=2,nk
if(modk(q)-modk(q-1).ne.0.D0) then
i=i+1
end if
end do

! We have "i" values of wave number. Let's allocate this final array with "i" terms

allocate(Sf(i,2))

! Setting the initial terms

Sf(1,1)=modk(1)
Sf(1,2)=Sk(1)

j=1
a=1

! Taking an average over different kx that originate different values of S (due to the other directions)

do q=3,nk
if(modk(q)-modk(q-1).eq.0.D0) then
Sf(a+1,1)=modk(q)
j=j+1
else
Sf(a+1,2)=sum(Sk(a+1:a+j))/j
j=1
a=a+1
end if
end do

! Writing in a file
write(6*rea,*) 'zone t="Sxk_y"'
do q=2,i-1
if(Sf(q,1).ne.0.D0)then
write(6*rea,'(F20.12,F20.12)') Sf(q,1),Sf(q,2)
end if
end do

deallocate(Sf)

!************************************************************* kz *********************************************************************!

do q=1,nk
modk(q) = DKR*kz(q)
end do

! Organizing the table S x k_x in ascending order

 call ordem_crescente(modk,Sk,nk)

! Let's first count how many wave numbers in the x directions we have

i=1

do q=2,nk
if(modk(q)-modk(q-1).ne.0.D0) then
i=i+1
end if
end do

! We have "i" values of wave number. Let's allocate this final array with "i" terms

allocate(Sf(i,2))

! Setting the initial terms

Sf(1,1)=modk(1)
Sf(1,2)=Sk(1)

j=1
a=1

! Taking an average over different kx that originate different values of S (due to the other directions)

do q=3,nk
if(modk(q)-modk(q-1).eq.0.D0) then
Sf(a+1,1)=modk(q)
j=j+1
else
Sf(a+1,2)=sum(Sk(a+1:a+j))/j
j=1
a=a+1
end if
end do

! Writing in a file
write(6*rea,*) 'zone t="Sxk_z"'
do q=2,i-1
if(Sf(q,1).ne.0.D0)then
write(6*rea,'(F20.12,F20.12)') Sf(q,1),Sf(q,2)
end if
end do

deallocate(Sf)


deallocate(kx)
deallocate(ky)
deallocate(kz)
deallocate(S)
deallocate(SI)
deallocate(Saux)
deallocate(Sk)
deallocate(modk)
deallocate(auxper)
deallocate(auxper2)
deallocate(auxper3)

end subroutine fator_estrutura

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: ordem_crescente			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine resposible for allocating numbers in !
! ascending order				  !
!*************************************************!

subroutine ordem_crescente(k,s,n)

real k(n), s(n)
integer i,j, q
real temp1, temp2
real temp3, temp4


do j=1,n-1
do q=j+1,n

aux=k(q)-k(j)

if(aux.lt.0.0) then

temp1= k(q)
temp2= k(j)

temp3= s(q)
temp4= s(j)

k(j)=temp1
k(q)=temp2
s(j)=temp3
s(q)=temp4

end if

end do
end do
end subroutine ordem_crescente

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: media_dipolo			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine resposible for calculating an ensem- !
! ble average of the particles orientation to cal-!
! culate the suspension's magnetization		  !
!*************************************************!

subroutine media_dipolo(U,N,rea,media,direcao)
integer N,rea
integer i,j,k
integer direcao
real U(rea,N,3)
real Di(rea,N,3)
real media

Di=U

media=sum(Di(:,:,direcao))/(N*rea)

end subroutine media_dipolo

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: parte_ativa			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine resposible for calculating the number!
! of particles that belongs to an active part of  !
! the suspension far from the system's boundaries.!
!
! This subroutine has not been used in recent ver-!
! sion due to the abscence of physical walls on   !
! the system. But in future versions, if anyone   !
! decides to work of wall bounded systems, this   !
! subroutine could be useful to the statistics.   !
!*************************************************!

subroutine parte_ativa(X,N,rea,l,h)
integer N,rea
integer i,j,k
real l,h,D
real X(rea,N,3)
real distancia
integer contador(rea)
integer partat(rea)

distancia=tam_ativa
contador=0
 !$OMP PARALLEL DO
do k=1,rea
do i=1,N

if(X(k,i,1).ge.distancia) then
if(X(k,i,1).le.(l-distancia)) then
if(X(k,i,2).ge.distancia) then
if(X(k,i,2).le.(l-distancia)) then
if(X(k,i,3).ge.distancia) then
if(X(k,i,1).le.(h-distancia)) then

contador(k)=contador(k)+1

end if
end if
end if
end if
end if
end if

end do
end do
 !$OMP END PARALLEL DO

 !$OMP PARALLEL DO
do k=1,rea
partat(k)=contador(k)
end do
 !$OMP END PARALLEL DO

write(*,*) 'O numero de particulas na parte ativa e de:',sum(partat)

end subroutine parte_ativa

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: acopla_particulas			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine resposible for sticking particles to !
! each other as soon as they form an aggregate.   !
! This subroutine prevents aggregate breakups.    !
! WARNING: NEEDS MORE TESTS			  !
!*************************************************!

subroutine acopla_particulas(X,N,D,rea,Di)

integer N,i,j,k,s, iter, iter2
integer LOOP, AGRE, PARTAGRE
integer AGREGPART(N/2)
integer PARTICULATESTE
integer PARTPORAGRE
integer D, rea
integer agregtotal(N)
integer parteste(N)
integer agreteste((N/2),N)
real X(rea,N,3),dt
integer G(N)
integer A(N/2,N)
real aux1,aux2,aux3,aux4,aux5
real Di(rea,N,3)

G=0

if(D.eq.1) then
agregtotal=0
end if

LOOP=0
AGRE=0
PARTAGRE=0
PARTICULATESTE=0
AGREGPART=0
parteste=0
agreteste=0

i=1
j=2

10 if(parteste(i).ne.0) then
LOOP=1
end if

if(parteste(j).ne.0) then
LOOP=1
end if


if(i.ne.j) then
r=(((X(D,i,1)-X(D,j,1))**2.0)+((X(D,i,2)-X(D,j,2))**2.0)+((X(D,i,3)-X(D,j,3))**2.0))**0.5

if(r.le.2.06) then

if(LOOP.eq.0.0) then
AGRE=AGRE+1
PARTAGRE=PARTAGRE+1
parteste(i)=i
parteste(j)=j

agreteste(AGRE,i)=i
agreteste(AGRE,j)=j

AGREGPART(AGRE)=PARTAGRE
else
do k=1,N/2
if(agreteste(k,i).eq.i) then
AGRE=k
end if
end do
parteste(i)=i
parteste(j)=j
PARTAGRE=AGREGPART(AGRE) + 1
AGREGPART(AGRE)=PARTAGRE
agreteste(AGRE,j)=j
end if


end if
end if

j=j+1

if(j.eq.N) then
i=i+1
j=i+1
LOOP=0
PARTAGRE=1
end if

if(i.ne.N-1) then
LOOP=0
PARTAGRE=1
go to 10
end if


do k=1,AGRE
do i=2,N
A(k,i)=agreteste(k,i)
end do
end do

! Coupling a variable of particles that belong to the same aggregate
! This quantity can be force, torque or dipole orientation of the 
! particles. This guarantee that the particles inside the same 
! aggregate structure will move and rotate together.

do k=1,AGRE

iter=1

do i=2,N
if(A(k,i).ne.0) then
if(iter.eq.1)then
aux3=Di(D,i,1)
aux4=Di(D,i,2)
aux5=Di(D,i,3)
iter=iter+1
else
Di(D,i,1)=aux3
Di(D,i,2)=aux4
Di(D,i,3)=aux5
end if
end if

end do
end do

end subroutine acopla_particulas

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: analise_de_agregados		  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Subroutine resposible for analyzing the number  !
! of aggregates and couting the particle distri-  !
! bution inside aggregates. It works together with!
! subroutine fdp_agregados			  ! 
! WARNING: NEEDS MORE TESTS			  !
!*************************************************!

subroutine analise_de_agregados(X,N,D,rea,s)

integer N, i,j,k,s
integer LOOP, AGRE, PARTAGRE
integer AGREGPART(N/2)
integer PARTICULATESTE
integer PARTPORAGRE
integer D, rea
integer agregtotal(N)
integer parteste(N)
integer agreteste((N/2),N)
real X(rea,N,3),dt
integer G(N)
real aux1,aux2

G=0

if(D.eq.1) then
agregtotal=0
end if

LOOP=0
AGRE=0
PARTAGRE=0
PARTICULATESTE=0
AGREGPART=0
parteste=0
agreteste=0

i=1
j=2

10 if(parteste(i).ne.0) then
LOOP=1
end if

if(parteste(j).ne.0) then
LOOP=1
end if


if(i.ne.j) then
r=(((X(D,i,1)-X(D,j,1))**2.0)+((X(D,i,2)-X(D,j,2))**2.0)+((X(D,i,3)-X(D,j,3))**2.0))**0.5

if(r.le.2.06) then

if(LOOP.eq.0.0) then
AGRE=AGRE+1
PARTAGRE=PARTAGRE+1
parteste(i)=i
parteste(j)=j

agreteste(AGRE,i)=i
agreteste(AGRE,j)=j

AGREGPART(AGRE)=PARTAGRE
else
do k=1,N/2
if(agreteste(k,i).eq.i) then
AGRE=k
end if
end do
parteste(i)=i
parteste(j)=j
PARTAGRE=AGREGPART(AGRE) + 1
AGREGPART(AGRE)=PARTAGRE
agreteste(AGRE,j)=j
end if

end if
end if

j=j+1

if(j.eq.N) then
i=i+1
j=i+1
LOOP=0
PARTAGRE=1
end if

if(i.ne.N-1) then
LOOP=0
PARTAGRE=1
go to 10
end if


! Counting the number of aggregates with N particles

do i=2,N
do k=1,AGRE
if(AGREGPART(k).eq.i) then
G(i)=G(i)+1
end if
end do
end do

do i=2,N
agregtotal(i)=agregtotal(i)+G(i)
end do


if(D.eq.rea) then

 call fdp_agregados(agregtotal,N,s)

end if

end subroutine analise_de_agregados

!*************************************************!
! 		     SIMMSUS			  ! 
!SUBROUTINE: fdp_agregados			  !					         
!Last update: 16/07/2023			  !
!*************************************************!

!*************************************************!
! Acessory subroutine used in subroutine:	  !
! analise_de_agregados				  ! 
! WARNING: NEEDS MORE TESTS			  !
!*************************************************!

subroutine fdp_agregados(A,N,k)
integer N,i,k
integer A(N)
integer numero_agregados
real aux2,aux1
real f(N)

numero_agregados=sum(A)
aux1=numero_agregados

do i=2,N
f(i)=A(i)!/aux1
end do

write(200*rea,*) f(2),f(3),f(4),f(5),f(6),f(7),f(8),f(9),f(10),k*dt

end subroutine fdp_agregados

!************************************************
end module subroutines
