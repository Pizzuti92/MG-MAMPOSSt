      SUBROUTINE powell(p,xi,n,np,ftol,iter,fret)
      INTEGER*4 iter,n,np,NMAX,ITMAX
      DOUBLE PRECISION fret,ftol,p(np),xi(np,np),func
      EXTERNAL func
      PARAMETER (NMAX=20,ITMAX=500)
CU    USES func,linmin
      INTEGER i,ibig,j
      DOUBLE PRECISION del,fp,fptt,t,pt(NMAX),ptt(NMAX),xit(NMAX)

      fret=func(p)

      do 11 j=1,n
        pt(j)=p(j)
11    continue
      iter=0
1     iter=iter+1
      fp=fret
      ibig=0
      del=0.d0
      do 13 i=1,n
        do 12 j=1,n
          xit(j)=xi(j,i)
12      continue
        fptt=fret

ccc      write(*,*) ' in Powell: ', fret,p

        call linmin(p,xit,n,fret)
        if(abs(fptt-fret).gt.del)then
          del=abs(fptt-fret)
          ibig=i
        endif
13    continue
      if(2.d0*abs(fp-fret).le.ftol*(abs(fp)+abs(fret)))return
c      if(iter.eq.ITMAX) pause 'powell exceeding maximum iterations'
      if(iter.eq.ITMAX) then
         write(*,*) 'powell exceeding maximum iterations'
         return
      endif
      do 14 j=1,n
        ptt(j)=2.d0*p(j)-pt(j)
        xit(j)=p(j)-pt(j)
        pt(j)=p(j)
14    continue
      fptt=func(ptt)
      if(fptt.ge.fp)goto 1
      t=2.d0*(fp-2.d0*fret+fptt)*(fp-fret-del)**2-del*(fp-fptt)**2
      if(t.ge.0.d0)goto 1
      call linmin(p,xit,n,fret)
      do 15 j=1,n
        xi(j,ibig)=xi(j,n)
        xi(j,n)=xit(j)
15    continue
      goto 1
      END
C  (C) Copr. 1986-92 Numerical Recipes Software )!0.d0
