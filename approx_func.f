!  firstchange.f90 
!
!  FUNCTIONS:
!  firstchange - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: firstchange
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

 
 program firstchange

    implicit none
    integer :: N, n0, i ,k , m0 = 4   ! N - разбиение 
    real :: c, d, h, sigma, gr, eps = 0,0001 ! c,d - левая/правая граница, h = (d-c)/N
    real,allocatable,dimension (:) :: X, Y, polinom, mas,gradient,V
    logical :: flag_of_iteration ! для отслеживания итерации ( если начало метода - .True. )
    
    write (*,*) 'c: input  '
    read *,c
    write (*,*) 'd: input  '
    read *,d
    write (*,*) 'N: input  '
    read *,N
    write (*,*) 'n0 input  '
    read *,n0
    write (*,*) 'hh'
    h = (d - c)/N
    !print *,h
     allocate (X(N+1))
    do i=1,N+1
        X(i) = c + (i - 1)*h
    end do
    do i=1,N+1
        print *,X(i)
    end do
    ! вычисляем таблицу значений функции на сетке
     allocate(Y(N+1))
     call function_Y(X,Y,N+1)
     ! модифицированный метод скорейшего спуска
    allocate (polinom(n0+m))
    allocate (V(n0+m))
    allocate (mas(n0+m))
    do i=1,n0+m
        polinom(i) = 0
    end do
    do i=1,n0+m
        mas(i) = 0
    end do
    do i=1,n0+m
        V(i) = 0
    end do
    do i=1,n0+m
        gradient(i) = 0
    end do
    do k=n0,m
        i = 0 
        flag_of_iteration =  .true.
        sigma = 1000
       ! sigma = 1
       ! write(*,*) i
        ! метод минимазиции , sigma - мера приближения
        do while (sigma>=eps)
            mas = polinom
            i = i + 1
            call count_V()
            polinom = polinom + lambda()*V
            if (flag_of_iteration .eqv. .True.) then 
                flag_of_iteration = .False.
            else
                sigma = count_sigma()
            end if
            gr = sqrt(dot_product(gradient,gradient))
            write(*,*)'номер итерации',i,'|grad(F)|=',gr
        end do
            
        
    deallocate(polinom)
    deallocate(Y)
    deallocate (X)
    deallocate (V)
    deallocate (mas)
    deallocate(gradient)
    
  contains
  real function lambda()
    
        real:: m1,m2,m3
        ! m2 = F(xk)
        m1 = count_F(polinom - V)
        m2 = count_F(polinom)
        m3 = count_F(polinom + V)
        lambda = (m1 - m3)/(2*(m1 - 2*m2 + m3))
        
    end function lamba
    
    
    subroutine count_V()
        real::mas
        
        if (flag_of_iteration .eqv. .True.) then ! начало v = - gradF
            call count_gradient()
            V = (-1)*gradient
        else    
            mas = dot_product(gradient,gradient) ! перемножение векторов
            call count_gradient()
            ! b = dot_product(gradient,gradient)/mas
            ! V = (-1)*gradient + b*V
            V = (-1)*gradient + (dot_product(gradient,gradient)/mas)*V
        end if
          
    end subroutine getV
    
  
    
    real function count_F(array)
    
        real,dimension::array
        integer::i
        count_F = 0
        do i=1,N+1
            count_F=count_F + (Y(i) - P(X(i),array))**2
        end do
        
    end function count_F
    
    real function P(xx,array)
        
        real::xx
        integer::i
        real,dimension::array
        
        P = 0
        do i =0,k
            P = array(k + 1 - i) + P*xx
        end do
        
    end function P
    
    subroutine  count_gradient()
    
      integer::i,j
      do j=1,k+1
        gradient(j) = 0
        do i=1,N+1
            gradient(j) = gradient(j) + 2*(X(i)**(j-1))*(P(X(i),polinom)-Y(i))
        end do
      end do
      
    end subroutine  count_gradient
    
    real function count_sigma()
        
        integer::i
        real::q
        
        mas=mas - polinom
        count_sigma=1000
        
        do i=1,k+1
            if (polinom(i)/=0) then
                q =abs(mas(i)/polinom(i))
                if (q<count_sigma) then
                    count_sigma = q
                end if
            end if
        end do
    end function count_sigma
    
     
    end program firstchange
    
       
    subroutine function_Y(X,Y,N) 
        integer:: N, i
        real, dimension (N) :: X, Y
        
        do i=1,N
            Y(i) = (X(i) + 1)*exp(-2*X(i))
            !print *,Y(i)
        end do
    
    end  subroutine function_Y
