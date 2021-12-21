

        SUBROUTINE splicoff2(x1a,x2a,ya,m,n,y2a)
        INTEGER m,n,NN
        REAL*8 x1a(m),x2a(n),y2a(m,n),ya(m,n),ypp1
        PARAMETER (NN=100)
! Derived from the method described in the Numerical Recipes for Fortran 77 textbook
! and adapted to the context.

        INTEGER j,k
        REAL*8 y2tmp(NN),ytmp(NN)
        do j=1,m
            do k=1,n
				ytmp(k)=ya(j,k)
            enddo
            ypp1=1.e30
            call SPLINE_CUBIC_SET(n,x2a,ytmp,3,ypp1,3,ypp1,y2tmp)


            do k=1,n
				y2a(j,k)=y2tmp(k)
            enddo
        enddo
        return
        END


        SUBROUTINE spli2d(x1a,x2a,ya,y2a,m,n,x1,x2,y)
        INTEGER m,n,NN
        REAL*8 x1,x2,y,x1a(m),x2a(n),y2a(m,n),ya(m,n),ypp1
        REAL*8 rjl(1), yrjl(1),yy1(1), yy2(1)
        PARAMETER (NN=100) 

        INTEGER j,k
        REAL*8 y2tmp(NN),ytmp(NN),yytmp(NN)
		
        do j=1,m 
            do k=1,n
                ytmp(k)=ya(j,k)
                y2tmp(k)=y2a(j,k)
            enddo
            rjl(1)=x2 
            call SPLINE_CUBIC_VAL(n,x2a,ytmp,y2tmp,rjl(1),yrjl(1),
     &  yy1,yy2)
            yytmp(j)=yrjl(1)
        enddo
        ypp1=1.30e30
        call SPLINE_CUBIC_SET(m,x1a,yytmp,3,ypp1,3,ypp1,y2tmp)
        
        rjl(1)=x1
        call SPLINE_CUBIC_VAL(m,x1a,yytmp,y2tmp,rjl(1),yrjl(1),
     &  yy1,yy2)
        y=yrjl(1)
        return
        END



