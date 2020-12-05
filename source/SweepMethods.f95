module SweepMethods
    implicit none
    contains

    subroutine progonka(A,B,C,D,Im,U)
      !A*U(i-1)+ B*U(i)+C*U(i+1) = D
      implicit none
      integer::i,Im
      real:: A(Im), B(Im),D(Im),C(Im),U(Im),F(Im),G(Im)
      real:: t
      F(1)=-C(1)/B(1)
      G(1)=D(1)/B(1)
      do i=2, Im-1
          t=B(i)+A(i)*F(i-1)
          F(i)=-C(i) / t
          G(i)=(D(i)-A(i)*G(i-1))/t
      enddo
      U(Im)=(D(Im)-A(Im)*G(Im-1))/(B(Im)+A(Im)*F(Im-1))
      do i=Im-1,1,-1
          U(i)=F(i)*U(i+1)+G(i)
      enddo
      return
    end subroutine progonka

end module
