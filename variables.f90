module variables

! Declaring all variabler used in the code


! INTEGERS (NON ALLOCATABLE)
integer i,j,k,q,s,a,b,c,d,g
integer N, auxper(5), auxper2(5), auxper3(5)
integer razao,iter,aux_periodico,posicao_campo
integer nb,nbr,npast,ncor,rea
integer auxint1,auxint2,auxint3,auxint4
integer nnr,semente(1),auxiliar1
integer inteiro, inteiro2, loop, loop2
integer n2,auxiliar_continua,DeAllocateStatus
integer teste2, nfreq, contfreqinteiro1, contfreqinteiro2, multiplofreq

! REAL (NON ALLOCATABLE)
real razao2, diferenca_interpol1, diferenca_interpol2
real diferenca, ajuda, lambda, kreal, aux_real, bifmax
real derivada1, derivada2, derivada3, trapezio, intervalo
real coeficiente1, coeficiente2,coeficiente3, coeficiente4
real coeficiente5, coeficiente6, coeficiente7, coeficiente8
real modulodipolo,kr2,dt,percentual,qsi, modk, kr
real mobilidade(3,3),resistencia(3,3), teste_identidade(3,3)
real mobilidade_self(3,3), mobilidade1(3,3), mobilidade2(3,3) 
real resistencia_self(3,3), resistencia1(3,3), resistencia2(3,3), pi
real contri1,contri2,contri3,contri4,contri5,contri6
real rij(3), rn(3),h2, freqmax,C1, C2, C3, C4, k_real
real l,h,phi,nr1, nr2, nr3, Pe, Per, ragreg
real r, alpha2,alpha,xcentro, ycentro, zcentro
real xmin, xmax, ymin, ymax, zmin, zmax
real modrij, modrand,auxcont, freq, shearratei
real termo1,termo2,termo3,termo4,termo5,eps
real tempo,ti, tf, tpros,dist,difusao
real reale, reale2, Str, freqcampo, freqbeat
real n3, modip, shearrate, St
real teste1, frequencia, umsobreene
real UMEDIA(3),SIGMA(3), konda(3), knormal(3)
real k1,k2,k3,k4,g1,g2,g3,g4

! REAL (ALLOCATABLES)
real, allocatable :: gpvetor(:)
real, allocatable :: campo(:)
real, allocatable :: y(:)
real, allocatable :: trap(:)
real, allocatable :: cof1(:,:)
real, allocatable :: cof2(:,:)
real, allocatable :: cof3(:,:)
real, allocatable :: cof4(:,:)
real, allocatable :: cof5(:,:)
real, allocatable :: cof6(:,:)
real, allocatable :: cof7(:,:)
real, allocatable :: cof8(:,:)
real, allocatable :: Di(:,:,:)
real, allocatable :: centro_massa(:,:)
real, allocatable :: aux1(:,:),aux2(:,:),aux3(:,:), aux4(:,:)
real, allocatable :: auxt(:,:),auxf(:,:)
real, allocatable :: torquereal(:,:),torquereciproco(:,:)
real, allocatable :: forcareal(:,:),forcareciproca(:,:)
real, allocatable :: contribuicao_self(:,:)
real, allocatable :: contribuicao_fisico(:,:)
real, allocatable :: contribuicao_reciproco(:,:)
real, allocatable :: contribuicoes(:,:)
real, allocatable :: FORCAS(:,:,:,:)
real, allocatable :: TORQUES(:,:,:,:)
real, allocatable :: FT(:,:,:)
real, allocatable :: Tt(:,:,:)
real, allocatable :: X(:,:,:)
real, allocatable :: W(:,:,:)
real, allocatable :: XI(:,:,:,:)
real, allocatable :: U(:,:,:)
real, allocatable :: V(:,:,:)
real, allocatable :: var(:,:)
real, allocatable :: errovar(:,:)
real, allocatable :: flutmag(:,:)
real, allocatable :: auxiliarcor(:,:)
real, allocatable :: nr(:)
real, allocatable :: tempototal(:)
real, allocatable :: aux_erro_vel(:,:,:)
real, allocatable :: aux_erro_var(:,:,:)
real, allocatable :: variancia(:,:,:)
real, allocatable :: flut(:,:,:,:)
real, allocatable :: auxcor(:,:,:)
real, allocatable :: errocor(:,:)
real, allocatable :: funcaor(:,:)
real, allocatable :: dif(:,:)
real, allocatable :: velmedia(:,:,:)
real, allocatable :: vmedia(:,:)
real, allocatable :: errovmedia(:,:)
real, allocatable :: fmedia(:,:)
real, allocatable :: errofmedia(:,:)
real, allocatable :: flut_tempo(:,:)
real, allocatable :: usistema(:,:)
real, allocatable :: ILF(:,:)
real, allocatable :: ILR(:,:)
real, allocatable :: DIAM(:,:)
real, allocatable :: auxflut(:,:,:)
real, allocatable :: diarand(:)
real, allocatable :: hidrodinamica_aux1(:,:)
real, allocatable :: hidrodinamica_aux2(:,:)
real, allocatable :: hidro1(:,:)
real, allocatable :: hidro2(:,:)
real, allocatable :: magtempo(:,:)
real, allocatable :: beta(:,:)
real, allocatable :: difusaoaux(:,:)

! LOGICAL
logical beating,posicao,fator,velocidade,estatistica, dipolo_ordenado
logical estatica,ordenado,torque,externo,gravidade,leito,mistura
logical gravadipolo,grafmag,browniano,tmagper,fmagper,ligaih
logical periodicidade,shear,oscillatory,agregado_inicial,continua
logical polidispersidade,printphi,oscilacampo,inertia,rotating
logical duffing,bifurcation, bifshear, ovito,magpart

! STRINGS
 character(3) rea_char
 character(23) linha1
 character(20) linha2
 character(25) linha3
 
end module variables
