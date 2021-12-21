c----------------------------------------------------------------------------- 
        SUBROUTINE XMIDMEAN(XDATA,N,XMID)
c----------------------------------------------------------------------------- 
 
c--- This subroutine calculates the MIDMEAN for a set of N ORDERED
c    statistics stored in XDATA. The MIDMEAN is defined to be the
c    mean of the central 50% of the data set. This corresponds to 
c    a 25% Trimmed Mean. The value of the MIDMEAN is returned as 
c    XMID and is defined:
c
c                   XMID = TRIM(.25) 
c
c    where TRIM(.25) is the 25% trimmed mean as defined above. For
c    more information on the MIDMEAN or TRIMMED MEAN see pages 312,
c    313 and pages 307, 308 in UREDA.
 
c**************************************************************************** 
 
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
        dimension xdata(n)
        data n1,n2,d1,d2,df,zero/1,2,1.0,2.0,0.25,0.0/

        ig = int(df*n)
        r = (df*n) - dfloat(ig)
        sum1 = zero
        do 11 i=ig + n2,n-ig-n1
        sum1 = sum1 + xdata(i)
11      continue
        sum3 = (d1-r)*(xdata(ig+n1) + xdata(n-ig))
        xmid = (d1/(dfloat(n)*(d1-(d2*df))))*(sum3+sum1)

        return
        end
