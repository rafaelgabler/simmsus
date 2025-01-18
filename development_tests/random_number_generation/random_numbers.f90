program random_test

integer n, n_inter, i,j
real minimo, maximo, range_inter
real, allocatable :: sequencia_gerada(:)
real, allocatable :: sequencia_leandro(:)
integer, allocatable :: numeros_por_intervalo(:)
integer, allocatable :: numeros_por_intervalo2(:)

n=1000000
n_inter=100
minimo = -1
maximo = 1
range_inter = (maximo-minimo)/n_inter

allocate(sequencia_gerada(n))
allocate(sequencia_leandro(n))
allocate(numeros_por_intervalo(n_inter))
allocate(numeros_por_intervalo2(n_inter))

call randomica(minimo,maximo,sequencia_gerada,n,12)

do i=1,n
do j=1,n_inter
if(sequencia_gerada(i).ge.minimo+(j-1)*range_inter.and.sequencia_gerada(i).le.minimo+(j*range_inter)) then
numeros_por_intervalo(j)=numeros_por_intervalo(j)+1
end if
end do
end do


open(1,file='rnd1M.out')
do j=1,n
read(1,*) sequencia_leandro(j)
end do

do i=1,n
do j=1,n_inter
if(sequencia_leandro(i).ge.minimo+(j-1)*range_inter.and.sequencia_leandro(i).le.minimo+(j*range_inter)) then
numeros_por_intervalo2(j)=numeros_por_intervalo2(j)+1
end if
end do
end do

open(2,file='comparacao.plt')
write(2,*) 'variables="N","freq_simmsus","freq_simmsus_cpp"'
do j=1,n_inter
write(2,*) j, numeros_por_intervalo(j), numeros_por_intervalo2(j)
end do


end




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
