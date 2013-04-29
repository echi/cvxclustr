module constants
  implicit none
  integer, parameter :: dble_prec = kind(0.0d0)
  real(kind=dble_prec), parameter :: zero  = 0.0_dble_prec
  real(kind=dble_prec), parameter :: one   = 1.0_dble_prec
  real(kind=dble_prec), parameter :: two   = 2.0_dble_prec
  real(kind=dble_prec), parameter :: three = 3.0_dble_prec
  real(kind=dble_prec), parameter :: four  = 4.0_dble_prec
  real(kind=dble_prec), parameter :: five  = 5.0_dble_prec
  real(kind=dble_prec), parameter :: six   = 6.0_dble_prec
  real(kind=dble_prec), parameter :: seven = 7.0_dble_prec
  real(kind=dble_prec), parameter :: eight = 8.0_dble_prec
  real(kind=dble_prec), parameter :: nine  = 9.0_dble_prec
  real(kind=dble_prec), parameter :: ten  = 10.0_dble_prec
  real(kind=dble_prec), parameter :: half  = 0.5_dble_prec
end module constants

module hilbert
!
!    This module provides basic operations with vectors: inner products and common norms.
!
  use constants
  implicit none
  
contains

  function norm(vector,L)
!
!     This function returns the l_1 norm (L=1), the l_2 norm (L=2), 
!     or the l_infinity norm (L=3) of a vector.
!
    implicit none
    integer :: L, n
    real(kind=dble_prec) :: norm
    real(kind=dble_prec) :: vector(:)

    n = size(vector)

    select case(L)
    case(1)
       norm = sum(abs(vector))
    case(2)
       norm = sqrt(dot_product(vector,vector))
    case(3)
       norm = maxval(abs(vector))
    end select
  end function norm

  function inner_product(a,b)
    implicit none
    integer :: n
    real(kind=dble_prec) :: inner_product
    real(kind=dble_prec) :: a(:), b(:)

    n = size(a)
    inner_product = dot_product(a,b)

  end function inner_product

end module hilbert

module sort
  use constants

contains
  recursive subroutine quick_sort(list)
    !
    !     This subroutine sorts the real array LIST into nonincreasing
    !     order by a combination of quick sort and exchange sort.
    !
    !     Taken from Ken's code.
    !
    implicit none
    integer :: i,j,n
    real(kind=dble_prec) :: split,uniform
    real(kind=dble_prec) :: list(:)
    !
    n = size(list)
    !
    !     For small lists do an exchange sort.
    !
    if (n<=6) then
       do i = 1,n
          do j = i+1,n
             if (list(i)<list(j)) then
                call swap(list(i),list(j))
             end if
          end do
       end do
       return
    end if
    !
    !     Otherwise carry out a quick sort.           
    !
    call random_number(uniform)
    i = int(n*uniform)+1
    split = list(i)
    list(i) = list(1)
    list(1) = split
    i = 1
    do j = 2,n
       if (list(j)>=split) then
          i = i+1
          call swap(list(i),list(j))
       end if
    end do
    list(1) = list(i)
    list(i) = split
    if (i>2) call quick_sort(list(1:i-1))
    if (i+1<n) call quick_sort(list(i+1:n))
  end subroutine quick_sort
  
  subroutine swap(a,b)
    !
    !     This subroutine swaps two reals.
    !
    implicit none
    real(kind=dble_prec) :: a,b,temp
    !
    temp = a
    a = b
    b = temp
  end subroutine swap
end module sort

module proximal
  use constants
  use hilbert
  
contains
  
  subroutine proj_L2(x,n,px,tau)
    implicit none
    integer :: n
    real(kind=dble_prec) :: lv, px(n), tau, x(n)
    lv = norm(x,2)
    px = x
    if (lv > tau) then
       px = (tau/lv)*x
    end if
  end subroutine proj_L2

  subroutine proj_Linf(x,n,px,tau)
    implicit none
    integer :: i,n
    real(kind=dble_prec) :: ax(n), px(n), tau, x(n)

    ax = abs(x)
    do i=1,n
       px(i) = min(ax(i),tau)
    end do
    px = sign(px,x)
  end subroutine proj_Linf
  
  subroutine prox_L1(x,n,px,tau)
    implicit none
    integer :: i,n
    real(kind=dble_prec) :: px(n), tau, x(n)
      
    px = abs(x) - tau
    do i=1,n
       px(i) = max(px(i),zero)
    end do
    px = sign(px,x)
    
  end subroutine prox_L1
  
  subroutine prox_L2(x,n,px,tau)
    implicit none
    integer :: n
    real(kind=dble_prec) :: lv, px(n), tau, x(n)
    px = zero
 !   lv = dsqrt(dot_product(x,x))                                                                                                                                                                     
    lv = norm(x,2)
    
    if (lv.eq.zero) then
       px = x
    else
       px = max(zero,one-(tau/lv))*x
    end if
  end subroutine prox_L2

  subroutine prox_Linf(x,n,px,tau)
!
! This function performs the proximal mapping of tau * L-infinity norm.
! It is computed via Moreau decomposition and Projection onto
! an L1-ball of radius tau.
!
!   px = x - project_L1(x,tau)
!
    implicit none
    integer :: n
    real(kind=dble_prec) :: lv, px(n), tau, x(n)
    px = (one/tau)*x
    call proj_L1(x,n,px,tau)
    px = x - px
  end subroutine prox_Linf

  subroutine project_to_simplex(x,n,z)
    use sort
    implicit none
    integer :: j, n, rho
    real(kind=dble_prec) :: x(n), z
    real(kind=dble_prec) :: cumsum, mu(n), theta
    mu = x
    call quick_sort(mu)
    cumsum = mu(1)
    do j = 2,n
       cumsum = cumsum + mu(j)
       if (dble(j)*mu(j) - cumsum + z .LE. zero) then
          rho = j-1
          exit
       end if
       rho = j
    end do
    theta = (sum(mu(1:rho)) - z) / dble(rho)
    do j = 1,n
       x(j) = max(x(j) - theta, zero)
    end do
  end subroutine project_to_simplex

  subroutine proj_L1(x,n,px,tau)
    implicit none
    integer :: n
    real(kind=dble_prec) :: px(n), tau, x(n)
    px = abs(x)
    call project_to_simplex(px,n,tau)
    px = sign(px,x)
  end subroutine proj_L1

  subroutine prox(vector_in,n,vector_out,tau,type)
!
! This function is a wrapper for performing various proximal and projection mappings.
!
! Menu of types:
!  type = 1 : proximal : L1
!  type = 2 : proximal : L2
!  type = 3 : proximal : L-infinity
!  type = 4 : project  : L-infinity
!  type = 5 : project  : L2
!  type = 6 : project  : L1
!  type = 7 : project  : L12 *** Not implemented
!  type = 8 : proximal : L12 *** Not implememted
!
! Notes:
!  March 13, 2013: Added proximal mapping of L-infinity norm
!  March 13, 2013: Swapped order 1->4, 2->5, 3->6, 4->1, 5->2, 6->3
!
  implicit none
  integer :: n, type
  real(kind=dble_prec) :: vector_in(n), vector_out(n), tau

  select case(type)
    case(1)
       call prox_L1(vector_in,n,vector_out,tau)
    case(2)
       call prox_L2(vector_in,n,vector_out,tau)
    case(3)
       call prox_Linf(vector_in,n,vector_out,tau)
    case(4)
       call proj_Linf(vector_in,n,vector_out,tau)
    case(5)
       call proj_L2(vector_in,n,vector_out,tau)
    case(6)
       call proj_L1(vector_in,n,vector_out,tau)
    end select

  end subroutine prox

  subroutine project(vector_in,n,vector_out,tau,type)
!
! This function is a wrapper for performing various projection mappings.
! Note that the types are chosen to complement the order of the prox calls.
!
! Menu of types:
!  type = 1 : L-infinity
!  type = 2 : L2
!  type = 3 : L1
!
  implicit none
  integer :: n, type
  real(kind=dble_prec) :: vector_in(n), vector_out(n), tau

  select case(type)
    case(1)
       call proj_Linf(vector_in,n,vector_out,tau)
    case(2)
       call proj_L2(vector_in,n,vector_out,tau)
    case(3)
       call proj_L1(vector_in,n,vector_out,tau)
    end select

  end subroutine project

end module proximal

subroutine proxA(vector_in,n,vector_out,tau,type)
  use proximal
  implicit none
  integer :: n, type
  real(kind=dble_prec) :: vector_in(n), vector_out(n), tau

  call prox(vector_in,n,vector_out,tau,type)

end subroutine proxA

subroutine projectA(vector_in,n,vector_out,tau,type)
  use proximal
  implicit none
  integer :: n, type
  real(kind=dble_prec) :: vector_in(n), vector_out(n), tau
  
  call project(vector_in,n,vector_out,tau,type)
end subroutine projectA

module admm
  use constants
  implicit none
  
contains
  
  subroutine get_xbar(X,xbar,p,q)
    implicit none
    integer :: p,q
    real(kind=dble_prec) :: X(q,p), xbar(q)
    xbar = (one/p)*sum(X,2)
  end subroutine get_xbar

  subroutine residual_primal(U,V,ix,p,q,nK,residual)
!
! This subroutine computes the L2 norm of the primal residual.
!
    implicit none
    integer :: ix(nK,2),p,q,nK
    real(kind=dble_prec) :: U(q,p), V(q,nK), residual
 
    residual = sqrt(sum( (U(:,ix(:,1)) - U(:,ix(:,2)) - V)**2))

  end subroutine residual_primal
  
  subroutine tri2vecA(i,j,p,k,n)
    implicit none
    integer :: j,p,n
    integer :: i(n), k(n)
    k = p*(i-1) - i*(i-1)/2 + j - i
  end subroutine tri2vecA
  
  subroutine tri2vecB(i,j,p,k,n)
    implicit none
    integer :: i,p,n
    integer :: j(n), k(n)
    k = p*(i-1) - i*(i-1)/2 + j - i
  end subroutine tri2vecB

  subroutine residual_dual(V,V_old,p,q,nK,nu,residual)
!
! This subroutine computes the L2 norm of the dual residual
!
    use hilbert
    implicit none
    integer :: ii, p, q, nK
    integer :: ix(p), seq(p)
    real(kind=dble_prec) :: nu, residual, V(q,nK), V_old(q,nK)
    real(kind=dble_prec) :: L1(q), L2(q)
    
    seq = (/(ii, ii=1,p, 1)/)
! Loop unrolling: ii=1
    ii = 1
    call tri2vecB(ii,seq((ii+1):p),p,ix(1:(p-ii)),p-ii)
    L2 = sum(V(:,ix(1:(p-ii))),2)
    L2 = L2 - sum(V_old(:,ix(1:(p-ii))),2)
    residual = norm(L2,2)**2
! Loop over ii=2:p-1
    do ii=2,p-1
       call tri2vecA(seq(1:(ii-1)),ii,p,ix(1:(ii-1)),ii-1)
       L1 = sum(V(:,ix(1:(ii-1))),2)
       L1 = L1 - sum(V_old(:,ix(1:(ii-1))),2)
       call tri2vecB(ii,seq((ii+1):p),p,ix(1:(p-ii)),p-ii)
       L2 = sum(V(:,ix(1:(p-ii))),2)
       L2 = L2 - sum(V_old(:,ix(1:(p-ii))),2)
       residual = residual + norm(L1-L2,2)**2
    end do
! Loop unrolling: ii=p
    ii = p
    call tri2vecA(seq(1:(ii-1)),ii,p,ix(1:(ii-1)),ii-1)
    L1 = sum(V(:,ix(1:(ii-1))),2)
    L1 = L1 - sum(V_old(:,ix(1:(ii-1))),2)
    residual = nu*(residual + norm(L1,2)**2)
  end subroutine residual_dual

  subroutine update_U(X,Lambda,U,V,xbar,p,q,nK,nu)
    implicit none
    integer :: ii,p,q,nK
    integer :: ix(p), seq(p)
    real(kind=dble_prec) :: X(q,p), U(q,p), V(q,nK), Lambda(q,nK), xbar(q)
    real(kind=dble_prec) :: L1(q), L2(q), nu, omega
    seq = (/(ii, ii=1,p, 1)/)
    
    omega = one/(one + p*nu)
! Loop unrolling: ii=1
    ii = 1
    call tri2vecB(ii,seq((ii+1):p),p,ix(1:(p-ii)),p-ii)
    L2 = sum(Lambda(:,ix(1:(p-ii))),2)
    L2 = L2 + nu*sum(V(:,ix(1:(p-ii))),2)
    U(:,ii) = omega*(X(:,ii) + L2) + (one-omega)*xbar
! Loop over ii=2:p-1
    do ii=2,p-1
       call tri2vecA(seq(1:(ii-1)),ii,p,ix(1:(ii-1)),ii-1)
       L1 = sum(Lambda(:,ix(1:(ii-1))),2)
       L1 = L1 + nu*sum(V(:,ix(1:(ii-1))),2)
       call tri2vecB(ii,seq((ii+1):p),p,ix(1:(p-ii)),p-ii)
       L2 = sum(Lambda(:,ix(1:(p-ii))),2)
       L2 = L2 + nu*sum(V(:,ix(1:(p-ii))),2)
       U(:,ii) = omega*(X(:,ii) - L1 + L2) + (one-omega)*xbar
    end do
! Loop unrolling: ii = p
    ii = p
    call tri2vecA(seq(1:(ii-1)),ii,p,ix(1:(ii-1)),ii-1)
    L1 = sum(Lambda(:,ix(1:(ii-1))),2)
    L1 = L1 + nu*sum(V(:,ix(1:(ii-1))),2)
    U(:,ii) = omega*(X(:,ii) - L1) + (one-omega)*xbar
    
  end subroutine update_U

  subroutine update_V(U,Lambda,V,w,gamma,nu,ix,q,p,nK,type)
    use proximal
    implicit none
    integer :: kk,nK,p,q
    integer :: i,j,ix(nK,2),type
    real(kind=dble_prec) :: gamma, nu
    real(kind=dble_prec) :: U(q,p), Lambda(q,nK), V(q,nK), w(nK)
    real(kind=dble_prec) :: z(q)
  
    do kk=1,nK
       i = ix(kk,1)
       j = ix(kk,2)
       z = U(:,i) - U(:,j) - (one/nu)*Lambda(:,kk)
       call prox(z,q,V(:,kk),w(kk)*gamma/nu,type)
    end do
  
  end subroutine update_V

  subroutine update_Lambda(Lambda,U,V,nu,ix,q,p,nK)
    implicit none
    integer :: nK, p, q
    integer :: ix(nK,2)
    real(kind=dble_prec) :: nu, Lambda(q,nK), U(q,p), V(q,nK)
    Lambda = Lambda - nu*(U(:,ix(:,1)) - U(:,ix(:,2)) - V)
  end subroutine update_Lambda

end module admm

module ama
  use constants
  implicit none
  
contains
  
  subroutine loss_primal(X,U,gamma,ix,p,q,nK,w,output,type)
    use hilbert
    implicit none
    integer :: k,p,q,nK,type
    integer :: ix(nK,2)
    real(kind=dble_prec) :: DU(q), X(q,p), U(q,p), gamma, w(nK), output, penalty
    
    penalty = zero
    do k=1,nK
       DU = U(:,ix(k,1)) - U(:,ix(k,2))
       penalty = penalty + w(k)*norm(DU,type)
    end do
    output = half*sum((X - U)**two) + gamma*penalty
  end subroutine loss_primal

  subroutine loss_dual(X,Lambda,ix,p,q,nK,s1,s2,M1,M2,mix1,mix2,output)
    use hilbert
    implicit none
    integer :: ii, mix1, mix2, nK, p, q
    integer :: ix(nK,2), M1(mix1,p), M2(mix2,p), s1(p), s2(p)
    real(kind=dble_prec) :: output, Lambda(q,nK), X(q,p)
    real(kind=dble_prec) :: first_term, second_term, L1(q), L2(q)
    
    first_term = zero
    do ii=1,p
       L1 = zero
       L2 = zero
       if (s1(ii) > 0) then
          L1 = sum(Lambda(:,M1(1:s1(ii),ii)),2) 
       end if
       if (s2(ii) > 0) then
          L2 = sum(Lambda(:,M2(1:s2(ii),ii)),2)
       end if
       first_term = first_term + sum( (L1-L2)**two )
    end do
    second_term = zero
    do ii=1,nK
       second_term = second_term + inner_product(X(:,ix(ii,1))-X(:,ix(ii,2)),Lambda(:,ii))
    end do
    output = -half*first_term - second_term
    
  end subroutine loss_dual

  subroutine update_U(X,Lambda,U,M1,M2,s1,s2,mix1,mix2,p,q,nK)
    implicit none
    integer :: ii,mix1,mix2,nK,p,q
    integer :: s1(p), s2(p), M1(mix1,p), M2(mix2,p)
    real(kind=dble_prec) :: Lambda(q,nK), U(q,p), X(q,p)
    
    do ii=1,p
       U(:,ii) = X(:,ii)
       if (s1(ii) > 0) then
          U(:,ii) = U(:,ii) + sum(Lambda(:,M1(1:s1(ii),ii)),2)
       end if
       if (s2(ii) > 0) then
          U(:,ii) = U(:,ii) - sum(Lambda(:,M2(1:s2(ii),ii)),2)
       end if
    end do
    
  end subroutine update_U
  
  subroutine update_V(U,Lambda,V,w,gamma,nu,ix,q,p,nK,type)
    use proximal
    implicit none
    integer :: kk,nK,p,q
    integer :: i,j,ix(nK,2),type
    real(kind=dble_prec) :: gamma, nu
    real(kind=dble_prec) :: U(q,p), Lambda(q,nK), V(q,nK), w(nK)
    real(kind=dble_prec) :: z(q)
    
    do kk=1,nK
       i = ix(kk,1)
       j = ix(kk,2)
       z = U(:,i) - U(:,j) - (one/nu)*Lambda(:,kk)
       call prox(z,q,V(:,kk),w(kk)*gamma/nu,type)
    end do
  end subroutine update_V
  
  subroutine update_Lambda(Lambda,U,nu,gamma,ix,q,p,nK,w,type)
    use proximal
    implicit none
    integer :: ii, nK, p, q, type
    integer :: ix(nK,2)
    real(kind=dble_prec) :: gamma, Lambda(q,nK), nu, U(q,p), w(nK), Z(q,nK), buffer(q)
    Z = Lambda - nu*(U(:,ix(:,1)) - U(:,ix(:,2)))
    do ii=1,nK
       call project(Z(:,ii),q,Lambda(:,ii),gamma*w(ii),type)
    end do
  end subroutine update_Lambda

end module ama

subroutine convex_cluster_ama(X,Lambda,U,V,q,p,nK,ix,w,gamma,nu,eta,s1,s2,M1,M2,&
     mix1,mix2,primal,dual,max_iter,iter,tol,type)
  use ama
  implicit none
  integer :: iter, max_iter, mix1, mix2, nK, p, q, type
  integer :: ix(nK,2), M1(mix1,p), M2(mix2,p), s1(p), s2(p)
  real(kind=dble_prec) :: X(q,p), Lambda(q,nK), U(q,p), V(q,nK), w(nK), gamma, nu, eta
  real(kind=dble_prec) :: primal(max_iter), dual(max_iter), tol
  real(kind=dble_prec) :: Lambda_old(q,nK), Ld(q,nK)
  real(kind=dble_prec) :: fp, fd

  do iter=1,max_iter
     Lambda_old = Lambda
     call update_U(X,Lambda_old,U,M1,M2,s1,s2,mix1,mix2,p,q,nK)
     call update_Lambda(Lambda,U,nu,gamma,ix,q,p,nK,w,type)
     call loss_primal(X,U,gamma,ix,p,q,nK,w,fp,type)
     primal(iter) = fp
     call loss_dual(X,Lambda,ix,p,q,nK,s1,s2,M1,M2,mix1,mix2,fd)
     dual(iter) = fd
     if (fp-fd < tol) exit
  end do
  call update_V(U,Lambda,V,w,gamma,nu,ix,q,p,nK,type)
end subroutine convex_cluster_ama

subroutine convex_cluster_ama_backtrack(X,Lambda,U,V,q,p,nK,ix,w,gamma,nu,eta,s1,s2,M1,M2,&
     mix1,mix2,primal,dual,max_iter,iter,tol,type)
  use ama
  use constants
  implicit none
  integer :: iter, max_iter, mix1, mix2, nK, p, q, type
  integer :: ix(nK,2), M1(mix1,p), M2(mix2,p), s1(p), s2(p)
  real(kind=dble_prec) :: Lambda(q,nK), U(q,p), V(q,nK), X(q,p), w(nK)
  real(kind=dble_prec) :: eta, gamma, nu
  real(kind=dble_prec) :: primal(max_iter), dual(max_iter), tol
  real(kind=dble_prec) :: Lambda_old(q,nK), Ld(q,nK)
  real(kind=dble_prec) :: f, f_last, fp, fd, qloss

  do iter=1,max_iter
     Lambda_old = Lambda
     call update_U(X,Lambda_old,U,M1,M2,s1,s2,mix1,mix2,p,q,nK)
     call loss_dual(X,Lambda_old,ix,p,q,nK,s1,s2,M1,M2,mix1,mix2,f_last)
     do
        Lambda = Lambda_old
        call update_Lambda(Lambda,U,nu,gamma,ix,q,p,nK,w,type)
        Ld = Lambda - Lambda_old
!        qloss = sum(Ld*(U(:,ix(:,1)) - U(:,ix(:,2)))) + (half/nu)*sum(Ld**2) - f_last
        qloss = -(half/nu)*sum(Ld*Ld) - f_last
        call loss_dual(X,Lambda,ix,p,q,nK,s1,s2,M1,M2,mix1,mix2,f)
        if (qloss+f.ge.zero) then
           exit
        end if
        nu = nu/eta
     end do
     call loss_primal(X,U,gamma,ix,p,q,nK,w,fp,type)
     primal(iter) = fp
     call loss_dual(X,Lambda,ix,p,q,nK,s1,s2,M1,M2,mix1,mix2,fd)
     dual(iter) = fd
!     if (fp-fd < tol) exit
     if (fp-fd < tol*(one + half*(fp+fd))) exit
  end do
  iter = min(iter,max_iter)
  call update_V(U,Lambda,V,w,gamma,nu,ix,q,p,nK,type)
end subroutine convex_cluster_ama_backtrack

subroutine convex_cluster_ama_fista(X,Lambda,U,V,q,p,nK,ix,w,gamma,nu,eta,s1,s2,M1,M2,&
     mix1,mix2,primal,dual,max_iter,iter,tol,type)
!
! This subroutine solves the convex clustering problem with a fixed step size using FISTA extrapolation steps.
!
  use ama
  use constants
  implicit none
  integer :: iter, max_iter, mix1, mix2, nK, p, q, type
  integer :: ix(nK,2), M1(mix1,p), M2(mix2,p), s1(p), s2(p)
  real(kind=dble_prec) :: X(q,p), Lambda(q,nK), U(q,p), V(q,nK), w(nK), gamma, nu, eta
  real(kind=dble_prec) :: primal(max_iter), dual(max_iter), tol
  real(kind=dble_prec) :: Lambda_old(q,nK), S(q,nK)
  real(kind=dble_prec) :: fp, fd

  do iter=1,2
     Lambda_old = Lambda
     call update_U(X,Lambda_old,U,M1,M2,s1,s2,mix1,mix2,p,q,nK)
     call update_Lambda(Lambda,U,nu,gamma,ix,q,p,nK,w,type)
     call loss_primal(X,U,gamma,ix,p,q,nK,w,fp,type)
     primal(iter) = fp
     call loss_dual(X,Lambda,ix,p,q,nK,s1,s2,M1,M2,mix1,mix2,fd)
     dual(iter) = fd
  end do
  do iter=3,max_iter
     S = Lambda + (dble(iter-2)/dble(iter+1))*(Lambda - Lambda_old)     
     call update_U(X,S,U,M1,M2,s1,s2,mix1,mix2,p,q,nK)
     call update_Lambda(S,U,nu,gamma,ix,q,p,nK,w,type)
     Lambda_old = Lambda
     Lambda = S
     call loss_primal(X,U,gamma,ix,p,q,nK,w,fp,type)
     primal(iter) = fp
     call loss_dual(X,Lambda,ix,p,q,nK,s1,s2,M1,M2,mix1,mix2,fd)
     dual(iter) = fd
     if (fp-fd < tol*(one+half*(fp+fd))) exit
!     if (fp-fd < tol) exit
  end do
  call update_V(U,Lambda,V,w,gamma,nu,ix,q,p,nK,type)
end subroutine convex_cluster_ama_fista

subroutine convex_cluster_ama_fista_backtrack(X,Lambda,U,V,q,p,nK,ix,w,gamma,nu,eta,s1,s2,M1,M2,&
mix1,mix2,primal,dual,max_iter,iter,tol,type)
  use ama
  use constants
  implicit none
  integer :: iter, max_iter, mix1, mix2, nK, p, q, type
  integer :: ix(nK,2), M1(mix1,p), M2(mix2,p), s1(p), s2(p)
  real(kind=dble_prec) :: Lambda(q,nK), U(q,p), V(q,nK), X(q,p), w(nK)
  real(kind=dble_prec) :: eta, gamma, nu
  real(kind=dble_prec) :: primal(max_iter), dual(max_iter), tol
  real(kind=dble_prec) :: Lambda_old(q,nK), G(q,nK)
  real(kind=dble_prec) :: f, f_last, fp, fd, qloss
  real(kind=dble_prec) :: alpha, alpha_old, S_old(q,nK), S(q,nK)

  alpha_old = one
  S_old = Lambda
  Lambda_old = S_old
  do iter=1,max_iter
     call update_U(X,Lambda_old,U,M1,M2,s1,s2,mix1,mix2,p,q,nK)
     call loss_dual(X,Lambda_old,ix,p,q,nK,s1,s2,M1,M2,mix1,mix2,f_last)
     do
        Lambda = Lambda_old
        call update_Lambda(Lambda,U,nu,gamma,ix,q,p,nK,w,type)
        G = ten**(log10(Lambda)-log10(nu)) - ten**(log10(Lambda_old) - log10(nu))
        qloss = sum(G*(U(:,ix(:,1)) - U(:,ix(:,2)))) + (half*nu)*sum(G*G) - f_last
        call loss_dual(X,Lambda,ix,p,q,nK,s1,s2,M1,M2,mix1,mix2,f)
        if (qloss+f.ge.zero) then
           exit
        end if
        nu = nu/eta
     end do
     call loss_primal(X,U,gamma,ix,p,q,nK,w,fp,type)
     primal(iter) = fp
     call loss_dual(X,Lambda,ix,p,q,nK,s1,s2,M1,M2,mix1,mix2,fd)
     dual(iter) = fd
     if (fp-fd < tol*(one + half*(fp+fd))) exit
     S = Lambda
     alpha = half*(one + sqrt(one + four*(alpha_old*alpha_old)))
     Lambda = S + ((alpha_old-one)/alpha)*(S-S_old)
     S_old = S
     alpha_old = alpha
     Lambda_old = Lambda
  end do
  iter = min(iter,max_iter)
  call update_V(U,S_old,V,w,gamma,nu,ix,q,p,nK,type)
end subroutine convex_cluster_ama_fista_backtrack

subroutine convex_cluster_admm(X,Lambda,U,V,q,p,nK,ix,w,gamma,nu,primal,dual,max_iter,iter,tol,type)
  use admm
  use hilbert
  implicit none
  integer :: iter, max_iter, nK, p, q, type
  integer :: ix(nK,2)
  real(kind=dble_prec) :: gamma, Lambda(q,nK), nu, U(q,p), V(q,nK), V_old(q,nK), w(nK), X(q,p)
  real(kind=dble_prec) :: primal(max_iter), dual(max_iter), tol
  real(kind=dble_prec) :: xbar(q), rp, rd, upsilon, LambdaNorm
  
  call get_xbar(X,xbar,p,q)
  
  upsilon = zero
  do iter=1,p
     upsilon = upsilon + norm(X(:,iter) - xbar,2)**2
  end do
  upsilon = sqrt(upsilon)
  
  do iter=1,max_iter
     V_old = V
     call update_U(X,Lambda,U,V,xbar,p,q,nK,nu)
     call update_V(U,Lambda,V,w,gamma,nu,ix,q,p,nK,type)
     call update_Lambda(Lambda,U,V,nu,ix,q,p,nK)
     LambdaNorm = sqrt(sum(Lambda**2))
     call residual_primal(U,V,ix,p,q,nK,rp)
     primal(iter) = rp
     call residual_dual(V,V_old,p,q,nK,nu,rd)
     dual(iter) = rd
     !     if (max(rp,rd) < tol) exit
     if (rp*LambdaNorm + upsilon*rd < tol) exit
  end do
  iter = min(iter,max_iter)
end subroutine convex_cluster_admm

subroutine convex_cluster_admm_acc(X,Lambda,U,V,q,p,nK,ix,w,gamma,nu,primal,dual,max_iter,iter,tol,type)
  use admm
  use hilbert
  implicit none
  integer :: iter,q,p,nK,max_iter,type
  integer :: ix(nK,2)
  real(kind=dble_prec) :: X(q,p), Lambda(q,nK), U(q,p), V(q,nK), V_old(q,nK), w(nK), gamma, nu
  real(kind=dble_prec) :: primal(max_iter), dual(max_iter), tol
  real(kind=dble_prec) :: xbar(q), rp, rd
  real(kind=dble_prec) :: alpha, alpha_old, S_old(q,nK), Lambda_old(q,nK), Ld(q,nK), r_last
  real(kind=dble_prec) :: upsilon, LambdaNorm
  
  call get_xbar(X,xbar,p,q)
  upsilon = zero
  do iter=1,p
     upsilon = upsilon + norm(X(:,iter) - xbar,2)**2
  end do
  upsilon = sqrt(upsilon)
  
  V_old = V
  Lambda_old = Lambda
  alpha_old = one
  r_last = huge(one)
  do iter=1,max_iter
     call update_U(X,Lambda,U,V,xbar,p,q,nK,nu)
     call update_V(U,Lambda,V,w,gamma,nu,ix,q,p,nK,type)
     call update_Lambda(Lambda,U,V,nu,ix,q,p,nK)
     LambdaNorm = sqrt(sum(Lambda**2))
     call residual_primal(U,V,ix,p,q,nK,rp)
     primal(iter) = rp
     call residual_dual(V,V_old,p,q,nK,nu,rd)
     dual(iter) = rd
     !     if (max(rp,rd) < tol) exit
     if (rp*LambdaNorm + upsilon*rd < tol) exit
     if (max(rp,rd) > r_last) then
        alpha_old = one
        r_last = huge(one)
     else
        alpha = half*(one + sqrt(1 + four*(alpha_old**2)))
        V = V + ((alpha_old-one)/alpha)*(V-V_old)
        Lambda = Lambda + ((alpha_old-one)/alpha)*(Lambda-Lambda_old)
        r_last = max(rp,rd)
        alpha_old = alpha
     end if
     V_old = V
     Lambda_old = Lambda
  end do
end subroutine convex_cluster_admm_acc

