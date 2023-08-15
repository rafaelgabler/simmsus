subroutine input

use variables

! Reads the configuration file "simconfig.dat" containing the necessary data to run the simulations
! To be used with the code SIMSUS
! Author: Rafael Gabler Gontijo

505   FORMAT(1X,A40,1X,E11.4E2)
507   FORMAT(1X,A40,1X,I6)
508   FORMAT(1X,A40,L10)

open(3,file='simconfig.dat')
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ (3,508) texto,torque
      READ (3,508) texto,magpart
      READ (3,508) texto,estatica
      READ (3,508) texto,gravidade
      READ (3,508) texto,leito
      READ (3,508) texto,browniano
      READ (3,508) texto,ligaih
      READ (3,508) texto,continua
      READ (3,508) texto,inertia
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ (3,508) texto,tmagper
      READ (3,508) texto,fmagper
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ (3,508) texto,posicao
      READ (3,508) texto,velocidade
      READ (3,508) texto,gravadipolo
      READ (3,508) texto,ovito
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ(3,507) texto,N
      READ(3,507) texto,rea
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ (3,508) texto,ordenado
      READ (3,508) texto,mistura
      READ (3,508) texto,agregado_inicial
      READ (3,508) texto,dipolo_ordenado
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ (3,508) texto,polidispersidade
      READ(3,505) texto,phi
      READ(3,507) texto,razao
      READ(3,505) texto,percentual
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ (3,508) texto,externo
      READ(3,507) texto,posicao_campo
      READ (3,508) texto,oscilacampo
      READ (3,508) texto,rotating
      READ (3,508) texto,duffing
      READ (3,508) texto,beating
      READ (3,508) texto,bifurcation
      READ(3,505) texto, freqcampo
      READ(3,505) texto, freqbeat
      READ(3,505) texto, C1
      READ(3,505) texto, C2
      READ(3,505) texto, C3
      READ(3,505) texto, C4
      READ(3,505) texto, bifmax
      READ(3,507) texto,nfreq
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ (3,508) texto,shear
      READ (3,508) texto,oscillatory
      READ (3,508) texto,bifshear
      READ(3,505) texto, shearrate
      READ(3,505) texto, freq
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ(3,505) texto,alpha2
      READ(3,505) texto,alpha
      READ(3,505) texto,Pe
      READ(3,505) texto,St
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ(3,505) texto,tempo
      READ(3,505) texto,dt
      READ(3,507) texto,n2
      READ(3,507) texto,iter
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ (3,'(A)') texto
      READ (3,508) texto,estatistica
      READ (3,508) texto,fator
      READ (3,508) texto,printphi
n3=n2
 close(3)

if(ligaih.or.fmagper.or.tmagper) then
periodicidade=.TRUE.
else
periodicidade=.FALSE.
end if

if(externo) then
grafmag=.TRUE.
end if

end subroutine input
