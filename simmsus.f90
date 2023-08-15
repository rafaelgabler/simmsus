program simmsus
use variables
implicit none

print *,'******************************************************************************'
print *,'*                                                                            *'
print *,'*                SIMMSUS - SIMULATION OF MAGNETIC SUSPENSIONS	              *'
print *,'*____________________________________________________________________________*'
print *,'*                                                                            *'
print *,'*                     PROF. RAFAEL GABLER GONTIJO, PhD                       *'
print *,'*____________________________________________________________________________*'
print *,'*                                                                            *'
print *,'*                        IN DEVELOPMENT SINCE 2009                           *'
print *,'*____________________________________________________________________________*'
print *,'*                                                                            *'
print *,'*                         LAST UPDATE: 16/07/2023                            *'
print *,'*                                                                            *'
print *,'******************************************************************************'
print *,''
print *,'******************************************************************************'
print *,'*                                                                            *'
print *,'*       Numerical simulation of magnetic suspensions of hard spheres         *'  
print *,'*____________________________________________________________________________*'
print *,'*                                                                            *'
print *,'*                       Langevin and Stokesian Dynamics                      *'
print *,'*                                                                            *'
print *,'******************************************************************************'
print *,''

! Collecting simulation data
 call input
 
! Start to count the simulation time
 call cpu_time(ti)
 
! Calling the main subroutine responsible for executing the simulation
 call main
 
! Calling the statistical analysis subroutine
if(estatistica) then
 call statistics
end if

! Stops counting the simulation time
 call cpu_time(tf)
 tpros=tf-ti
 
print *, 'TOTAL SIMULATION TIME:',tpros,'SECONDS'
write(*,*) ''
end program simmsus
