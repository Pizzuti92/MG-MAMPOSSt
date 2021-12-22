	PROGRAM MRSWFP
C
C       ==============================================================
C       Purpose: This program computes the radial prolate spheriodal 
C                functions of the first and second kinds, and their 
C                derivatives using subroutine RSWFP
C       Input :  m  --- Mode parameter, m = 0,1,2,...
C                n  --- Mode parameter, n = m,m+1,m+2,...
C                c  --- Spheroidal parameter
C                cv --- Characteristic value
C                x  --- Argument of radial function ( x > 1.0 )
C                KF --- Function code
C                       KF=1 for the first kind
C                       KF=2 for the second kind
C                       KF=3 for both the first and second kinds
C       Output:  R1F --- Radial function of the first kind
C                R1D --- Derivative of the radial function of
C                        the first kind
C                R2F --- Radial function of the second kind
C                R2D --- Derivative of the radial function of
C                        the second kind
C       Example:
C                KD= 1, m = 2, n = 3, c = 5.0 and cv =19.1359819110
C
C    x      R23(1)(c,x)   R23(1)'(c,x)    R23(2)(c,x)   R23(2)'(c,x)
C -------------------------------------------------------------------
C 1.(7).1 (-8) 1.6240735      1.6240735 ( 6)-3.0786785 (14) 3.0786784
C  1.005  (-3) 8.0600009      1.5998506     -6.2737713 ( 3) 1.2299041
C   1.1   (-1) 1.3578875      1.0702727 (-1)-4.3693218      3.5698419
C   1.5   (-1) 1.2002958 (-1)-8.0657929 (-1) 1.2623910 (-1) 4.8469847
C   5.0   (-2) 3.9455888 (-2) 4.2556949 (-2)-1.0099734 (-1) 2.0031280
C       ==============================================================
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION EG(200)
	WRITE(*,*)'Please enter: KF --- Kind choice code'
	WRITE(*,*)'KF = ?'
	READ(*,*)KF
	WRITE(*,10)KF
	WRITE(*,*)'Please enter m, n, c and x ( x > 1 )'
	READ(*,*)M,N,C,X
	CALL SEGV(M,N,C,1,CV,EG)
	WRITE(*,20)M,N,C,CV,X
	WRITE(*,*)
	CALL RSWFP(M,N,C,X,CV,KF,R1F,R1D,R2F,R2D)
	WRITE(*,*)'   x      Rmn(1)(c,x)      Rmn''(1)(c,x)     ',
     &            'Rmn(2)(c,x)      Rmn''(2)(c,x)'
	WRITE(*,*)'-----------------------------------------------',
     &            '---------------------------'
	IF (KF.EQ.1) THEN
	   WRITE(*,30)X,R1F,R1D
	ELSE IF (KF.EQ.2) THEN
	   WRITE(*,40)X,R2F,R2D
	ELSE IF (KF.EQ.3) THEN
	   WRITE(*,30)X,R1F,R1D,R2F,R2D
	ENDIF
	IF (KF.EQ.3) THEN
	   WRITE(*,50)R1F*R2D-R2F*R1D,1.0D0/(C*(X*X-1.0D0))
	   WRITE(*,60)
	ENDIF
10      FORMAT(1X,3HKF=,I3)
20      FORMAT(1X,2Hm=,I2,',   ',2Hn=,I2,',   ',2Hc=,F5.1,
     &         ',   ',4Hcv =,D20.12,',  ',2Hx=,F5.2)
30      FORMAT(1X,F5.2,4D17.8)
40      FORMAT(1X,F5.2,34X,4D17.8)
50      FORMAT(1X,/1X,'Wronskian check:',/1X,'Computed value =',
     &         D17.8,5X,'Exact value =',D17.8)
60      FORMAT(1X,/1X,'Caution: This check is not accurate if it ',
     &        'involves',/1X,'         the subtraction of two ',
     &        'similar numbers')
	END


	SUBROUTINE RSWFP(M,N,C,X,CV,KF,R1F,R1D,R2F,R2D)
C
C       ==============================================================
C       Purpose: Compute prolate spheriodal radial functions of the
C                first and second kinds, and their derivatives
C       Input :  m  --- Mode parameter, m = 0,1,2,...
C                n  --- Mode parameter, n = m,m+1,m+2,...
C                c  --- Spheroidal parameter
C                x  --- Argument of radial function ( x > 1.0 )
C                cv --- Characteristic value
C                KF --- Function code
C                       KF=1 for the first kind
C                       KF=2 for the second kind
C                       KF=3 for both the first and second kinds
C       Output:  R1F --- Radial function of the first kind
C                R1D --- Derivative of the radial function of
C                        the first kind
C                R2F --- Radial function of the second kind
C                R2D --- Derivative of the radial function of
C                        the second kind
C       Routines called:
C            (1) SDMN for computing expansion coefficients dk
C            (2) RMN1 for computing prolate and oblate radial 
C                functions of the first kind
C            (3) RMN2L for computing prolate and oblate radial
C                functions of the second kind for a large argument
C            (4) RMN2SP for computing the prolate radial function 
C                of the second kind for a small argument
C       ==============================================================
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION DF(200)
	KD=1
	CALL SDMN(M,N,C,CV,KD,DF)
	IF (KF.NE.2) THEN
	   CALL RMN1(M,N,C,X,DF,KD,R1F,R1D)
	ENDIF
	IF (KF.GT.1) THEN
	   CALL RMN2L(M,N,C,X,DF,KD,R2F,R2D,ID)
	   IF (ID.GT.-8) THEN
	      CALL RMN2SP(M,N,C,X,CV,DF,KD,R2F,R2D)
	   ENDIF
	ENDIF
	RETURN
	END


	SUBROUTINE SDMN(M,N,C,CV,KD,DF)
C
C       =====================================================
C       Purpose: Compute the expansion coefficients of the
C                prolate and oblate spheroidal functions, dk
C       Input :  m  --- Mode parameter
C                n  --- Mode parameter
C                c  --- Spheroidal parameter
C                cv --- Characteristic value
C                KD --- Function code
C                       KD=1 for prolate; KD=-1 for oblate
C       Output:  DF(k) --- Expansion coefficients dk;
C                          DF(1), DF(2), ... correspond to
C                          d0, d2, ... for even n-m and d1,
C                          d3, ... for odd n-m
C       =====================================================
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION A(200),D(200),G(200),DF(200)
	NM=25+INT(0.5*(N-M)+C)
	IF (C.LT.1.0D-10) THEN
	   DO 5 I=1,NM
5             DF(I)=0D0
	   DF((N-M)/2+1)=1.0D0
	   RETURN
	ENDIF   
	CS=C*C*KD
	IP=1
	IF (N-M.EQ.2*INT((N-M)/2)) IP=0
	DO 10 I=1,NM+2
	   IF (IP.EQ.0) K=2*(I-1)
	   IF (IP.EQ.1) K=2*I-1
	   DK0=M+K
	   DK1=M+K+1
	   DK2=2*(M+K)
	   D2K=2*M+K
	   A(I)=(D2K+2.0)*(D2K+1.0)/((DK2+3.0)*(DK2+5.0))*CS
	   D(I)=DK0*DK1+(2.0*DK0*DK1-2.0*M*M-1.0)/((DK2-1.0)
     &          *(DK2+3.0))*CS
	   G(I)=K*(K-1.0)/((DK2-3.0)*(DK2-1.0))*CS
10      CONTINUE
	FS=1.0D0
	F1=0.0D0
	F0=1.0D-100
	KB=0
	DF(NM+1)=0.0D0
	DO 30 K=NM,1,-1
	   F=-((D(K+1)-CV)*F0+A(K+1)*F1)/G(K+1)
	   IF (DABS(F).GT.DABS(DF(K+1))) THEN
	      DF(K)=F
	      F1=F0
	      F0=F
	      IF (DABS(F).GT.1.0D+100) THEN
		 DO 12 K1=K,NM
12                  DF(K1)=DF(K1)*1.0D-100
		 F1=F1*1.0D-100
		 F0=F0*1.0D-100
	      ENDIF  
	   ELSE
	      KB=K
	      FL=DF(K+1)
	      F1=1.0D-100
	      F2=-(D(1)-CV)/A(1)*F1
	      DF(1)=F1
	      IF (KB.EQ.1) THEN
		 FS=F2
	      ELSE IF (KB.EQ.2) THEN
		 DF(2)=F2
		 FS=-((D(2)-CV)*F2+G(2)*F1)/A(2)
	      ELSE 
		 DF(2)=F2
		 DO 20 J=3,KB+1
		    F=-((D(J-1)-CV)*F2+G(J-1)*F1)/A(J-1)
		    IF (J.LE.KB) DF(J)=F
		    IF (DABS(F).GT.1.0D+100) THEN
		       DO 15 K1=1,J
15                        DF(K1)=DF(K1)*1.0D-100
		       F=F*1.0D-100
		       F2=F2*1.0D-100
		    ENDIF  
		    F1=F2
20                  F2=F
		 FS=F
	      ENDIF
	      GO TO 35
	   ENDIF
30      CONTINUE
35      SU1=0.0D0
	R1=1.0D0
	DO 40 J=M+IP+1,2*(M+IP)
40         R1=R1*J
	SU1=DF(1)*R1
	DO 45 K=2,KB
	   R1=-R1*(K+M+IP-1.5D0)/(K-1.0D0)
45           SU1=SU1+R1*DF(K)
	SU2=0.0D0
	DO 50 K=KB+1,NM
	   IF (K.NE.1) R1=-R1*(K+M+IP-1.5D0)/(K-1.0D0)
	   SU2=SU2+R1*DF(K)
	   IF (DABS(SW-SU2).LT.DABS(SU2)*1.0D-14) GOTO 55
50         SW=SU2
55      R3=1.0D0
	DO 60 J=1,(M+N+IP)/2
60         R3=R3*(J+0.5D0*(N+M+IP))
	R4=1.0D0
	DO 65 J=1,(N-M-IP)/2
65         R4=-4.0D0*R4*J
	S0=R3/(FL*(SU1/FS)+SU2)/R4
	DO 70 K=1,KB
70         DF(K)=FL/FS*S0*DF(K)
	DO 75 K=KB+1,NM
75         DF(K)=S0*DF(K)
	RETURN
	END


	SUBROUTINE RMN1(M,N,C,X,DF,KD,R1F,R1D)
C
C       =======================================================
C       Purpose: Compute prolate and oblate spheroidal radial
C                functions of the first kind for given m, n,
C                c and x
C       Routines called:
C            (1) SCKB for computing expansion coefficients c2k
C            (2) SPHJ for computing the spherical Bessel
C                functions of the first kind     
C       =======================================================
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION CK(200),DF(200),SJ(0:251),DJ(0:251)
	EPS=1.0D-14
	IP=1
	NM1=INT((N-M)/2)
	IF (N-M.EQ.2*NM1) IP=0
	NM=25+NM1+INT(C)
	REG=1.0D0
	IF (M+NM.GT.80) REG=1.0D-200
	R0=REG
	DO 10 J=1,2*M+IP
10         R0=R0*J
	R=R0    
	SUC=R*DF(1)
	DO 15 K=2,NM
	   R=R*(M+K-1.0)*(M+K+IP-1.5D0)/(K-1.0D0)/(K+IP-1.5D0)
	   SUC=SUC+R*DF(K)
	   IF (K.GT.NM1.AND.DABS(SUC-SW).LT.DABS(SUC)*EPS) GO TO 20
15         SW=SUC
20      CONTINUE
	IF (X.EQ.0.0) THEN
	   CALL SCKB(M,N,C,DF,CK)
	   SUM=0.0D0
	   DO 25 J=1,NM
	      SUM=SUM+CK(J)
	      IF (DABS(SUM-SW1).LT.DABS(SUM)*EPS) GO TO 30
25            SW1=SUM
30         R1=1.0D0
	   DO 35 J=1,(N+M+IP)/2
35            R1=R1*(J+0.5D0*(N+M+IP))
	   R2=1.0D0
	   DO 40 J=1,M
40            R2=2.0D0*C*R2*J
	   R3=1.0D0
	   DO 45 J=1,(N-M-IP)/2
45            R3=R3*J
	   SA0=(2.0*(M+IP)+1.0)*R1/(2.0**N*C**IP*R2*R3)
	   IF (IP.EQ.0) THEN
	      R1F=SUM/(SA0*SUC)*DF(1)*REG
	      R1D=0.0D0
	   ELSE IF (IP.EQ.1) THEN
	      R1F=0.0D0
	      R1D=SUM/(SA0*SUC)*DF(1)*REG
	   ENDIF
	   RETURN
	ENDIF
	CX=C*X
	NM2=2*NM+M
	CALL SPHJ(NM2,CX,NM2,SJ,DJ)
	A0=(1.0D0-KD/(X*X))**(0.5D0*M)/SUC  
	R1F=0.0D0
	DO 50 K=1,NM
	   L=2*K+M-N-2+IP
	   IF (L.EQ.4*INT(L/4)) LG=1
	   IF (L.NE.4*INT(L/4)) LG=-1
	   IF (K.EQ.1) THEN
	      R=R0
	   ELSE
	      R=R*(M+K-1.0)*(M+K+IP-1.5D0)/(K-1.0D0)/(K+IP-1.5D0)
	   ENDIF
	   NP=M+2*K-2+IP
	   R1F=R1F+LG*R*DF(K)*SJ(NP)
	   IF (K.GT.NM1.AND.DABS(R1F-SW).LT.DABS(R1F)*EPS) GO TO 55
50         SW=R1F
55      R1F=R1F*A0
	B0=KD*M/X**3.0D0/(1.0-KD/(X*X))*R1F    
	SUD=0.0D0
	DO 60 K=1,NM
	   L=2*K+M-N-2+IP
	   IF (L.EQ.4*INT(L/4)) LG=1
	   IF (L.NE.4*INT(L/4)) LG=-1
	   IF (K.EQ.1) THEN
	      R=R0
	   ELSE
	      R=R*(M+K-1.0)*(M+K+IP-1.5D0)/(K-1.0D0)/(K+IP-1.5D0)
	   ENDIF
	   NP=M+2*K-2+IP
	   SUD=SUD+LG*R*DF(K)*DJ(NP)
	   IF (K.GT.NM1.AND.DABS(SUD-SW).LT.DABS(SUD)*EPS) GO TO 65
60         SW=SUD
65      R1D=B0+A0*C*SUD
	RETURN
	END


	SUBROUTINE SCKB(M,N,C,DF,CK)
C
C       ======================================================
C       Purpose: Compute the expansion coefficients of the
C                prolate and oblate spheroidal functions, c2k
C       Input :  m  --- Mode parameter
C                n  --- Mode parameter
C                c  --- Spheroidal parameter
C                DF(k) --- Expansion coefficients dk
C       Output:  CK(k) --- Expansion coefficients ck;
C                          CK(1), CK(2), ... correspond to
C                          c0, c2, ...
C       ======================================================
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION DF(200),CK(200)
	IF (C.LE.1.0D-10) C=1.0D-10
	NM=25+INT(0.5*(N-M)+C)
	IP=1
	IF (N-M.EQ.2*INT((N-M)/2)) IP=0
	REG=1.0D0
	IF (M+NM.GT.80) REG=1.0D-200
	FAC=-0.5D0**M
	DO 35 K=0,NM-1
	   FAC=-FAC
	   I1=2*K+IP+1
	   R=REG
	   DO 10 I=I1,I1+2*M-1
10            R=R*I
	   I2=K+M+IP
	   DO 15 I=I2,I2+K-1
15            R=R*(I+0.5D0)
	   SUM=R*DF(K+1)
	   DO 20 I=K+1,NM
	      D1=2.0D0*I+IP
	      D2=2.0D0*M+D1
	      D3=I+M+IP-0.5D0
	      R=R*D2*(D2-1.0D0)*I*(D3+K)/(D1*(D1-1.0D0)*(I-K)*D3)
	      SUM=SUM+R*DF(I+1)
	      IF (DABS(SW-SUM).LT.DABS(SUM)*1.0D-14) GOTO 25
20            SW=SUM
25         R1=REG
	   DO 30 I=2,M+K
30            R1=R1*I
35         CK(K+1)=FAC*SUM/R1
	RETURN
	END 


	SUBROUTINE SPHJ(N,X,NM,SJ,DJ)
C
C       =======================================================
C       Purpose: Compute spherical Bessel functions jn(x) and
C                their derivatives
C       Input :  x --- Argument of jn(x)
C                n --- Order of jn(x)  ( n = 0,1,��� )
C       Output:  SJ(n) --- jn(x)
C                DJ(n) --- jn'(x)
C                NM --- Highest order computed
C       Routines called:
C                MSTA1 and MSTA2 for computing the starting
C                point for backward recurrence
C       =======================================================
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION SJ(0:N),DJ(0:N)
	NM=N
	IF (DABS(X).EQ.1.0D-100) THEN
	   DO 10 K=0,N
	      SJ(K)=0.0D0
10            DJ(K)=0.0D0
	   SJ(0)=1.0D0
	   DJ(1)=.3333333333333333D0
	   RETURN
	ENDIF
	SJ(0)=DSIN(X)/X
	SJ(1)=(SJ(0)-DCOS(X))/X
	IF (N.GE.2) THEN
	   SA=SJ(0)
	   SB=SJ(1)
	   M=MSTA1(X,200)
	   IF (M.LT.N) THEN
	      NM=M
	   ELSE
	      M=MSTA2(X,N,15)
	   ENDIF
	   F0=0.0D0
	   F1=1.0D0-100
	   DO 15 K=M,0,-1
	      F=(2.0D0*K+3.0D0)*F1/X-F0
	      IF (K.LE.NM) SJ(K)=F
	      F0=F1
15            F1=F
	   IF (DABS(SA).GT.DABS(SB)) CS=SA/F
	   IF (DABS(SA).LE.DABS(SB)) CS=SB/F0
	   DO 20 K=0,NM
20            SJ(K)=CS*SJ(K)
	ENDIF      
	DJ(0)=(DCOS(X)-DSIN(X)/X)/X
	DO 25 K=1,NM
25         DJ(K)=SJ(K-1)-(K+1.0D0)*SJ(K)/X
	RETURN
	END


	INTEGER FUNCTION MSTA1(X,MP)
C
C       ===================================================
C       Purpose: Determine the starting point for backward  
C                recurrence such that the magnitude of    
C                Jn(x) at that point is about 10^(-MP)
C       Input :  x     --- Argument of Jn(x)
C                MP    --- Value of magnitude
C       Output:  MSTA1 --- Starting point   
C       ===================================================
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	A0=DABS(X)
	N0=INT(1.1*A0)+1
	F0=ENVJ(N0,A0)-MP
	N1=N0+5
	F1=ENVJ(N1,A0)-MP
	DO 10 IT=1,20             
	   NN=N1-(N1-N0)/(1.0D0-F0/F1)                  
	   F=ENVJ(NN,A0)-MP
	   IF(ABS(NN-N1).LT.1) GO TO 20
	   N0=N1
	   F0=F1
	   N1=NN
 10        F1=F
 20     MSTA1=NN
	RETURN
	END


	INTEGER FUNCTION MSTA2(X,N,MP)
C
C       ===================================================
C       Purpose: Determine the starting point for backward
C                recurrence such that all Jn(x) has MP
C                significant digits
C       Input :  x  --- Argument of Jn(x)
C                n  --- Order of Jn(x)
C                MP --- Significant digit
C       Output:  MSTA2 --- Starting point
C       ===================================================
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	A0=DABS(X)
	HMP=0.5D0*MP
	EJN=ENVJ(N,A0)
	IF (EJN.LE.HMP) THEN
	   OBJ=MP
	   N0=INT(1.1*A0)
	ELSE
	   OBJ=HMP+EJN
	   N0=N
	ENDIF
	F0=ENVJ(N0,A0)-OBJ
	N1=N0+5
	F1=ENVJ(N1,A0)-OBJ
	DO 10 IT=1,20
	   NN=N1-(N1-N0)/(1.0D0-F0/F1)
	   F=ENVJ(NN,A0)-OBJ
	   IF (ABS(NN-N1).LT.1) GO TO 20
	   N0=N1
	   F0=F1
	   N1=NN
10         F1=F
20      MSTA2=NN+10
	RETURN
	END

	REAL*8 FUNCTION ENVJ(N,X)
	DOUBLE PRECISION X
	ENVJ=0.5D0*DLOG10(6.28D0*N)-N*DLOG10(1.36D0*X/N)
	RETURN
	END


	SUBROUTINE RMN2L(M,N,C,X,DF,KD,R2F,R2D,ID)
C
C       ========================================================
C       Purpose: Compute prolate and oblate spheroidal radial
C                functions of the second kind for given m, n,
C                c and a large cx
C       Routine called:
C                SPHY for computing the spherical Bessel
C                functions of the second kind      
C       ========================================================
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION DF(200),SY(0:251),DY(0:251)
	EPS=1.0D-14
	IP=1
	NM1=INT((N-M)/2)
	IF (N-M.EQ.2*NM1) IP=0
	NM=25+NM1+INT(C)
	REG=1.0D0
	IF (M+NM.GT.80) REG=1.0D-200
	NM2=2*NM+M
	CX=C*X
	CALL SPHY(NM2,CX,NM2,SY,DY)
	R0=REG
	DO 10 J=1,2*M+IP
10         R0=R0*J
	R=R0    
	SUC=R*DF(1)
	DO 15 K=2,NM
	   R=R*(M+K-1.0)*(M+K+IP-1.5D0)/(K-1.0D0)/(K+IP-1.5D0)
	   SUC=SUC+R*DF(K)
	   IF (K.GT.NM1.AND.DABS(SUC-SW).LT.DABS(SUC)*EPS) GO TO 20
15         SW=SUC
20      A0=(1.0D0-KD/(X*X))**(0.5D0*M)/SUC
	R2F=0.0
	DO 50 K=1,NM
	   L=2*K+M-N-2+IP
	   IF (L.EQ.4*INT(L/4)) LG=1
	   IF (L.NE.4*INT(L/4)) LG=-1
	   IF (K.EQ.1) THEN
	      R=R0
	   ELSE
	      R=R*(M+K-1.0)*(M+K+IP-1.5D0)/(K-1.0D0)/(K+IP-1.5D0)
	   ENDIF
	   NP=M+2*K-2+IP
	   R2F=R2F+LG*R*(DF(K)*SY(NP))
	   EPS1=DABS(R2F-SW)
	   IF (K.GT.NM1.AND.EPS1.LT.DABS(R2F)*EPS) GO TO 55
50         SW=R2F
55      ID1=INT(LOG10(EPS1/DABS(R2F)+EPS))
	R2F=R2F*A0
	IF (NP.GE.NM2) THEN
	   ID=10
	   RETURN
	ENDIF
	B0=KD*M/X**3.0D0/(1.0-KD/(X*X))*R2F                
	SUD=0.0D0
	DO 60 K=1,NM
	   L=2*K+M-N-2+IP
	   IF (L.EQ.4*INT(L/4)) LG=1
	   IF (L.NE.4*INT(L/4)) LG=-1
	   IF (K.EQ.1) THEN
	      R=R0
	   ELSE
	      R=R*(M+K-1.0)*(M+K+IP-1.5D0)/(K-1.0D0)/(K+IP-1.5D0)
	   ENDIF
	   NP=M+2*K-2+IP
	   SUD=SUD+LG*R*(DF(K)*DY(NP))
	   EPS2=DABS(SUD-SW)
	   IF (K.GT.NM1.AND.EPS2.LT.DABS(SUD)*EPS) GO TO 65
60         SW=SUD
65      R2D=B0+A0*C*SUD       
	ID2=INT(LOG10(EPS2/DABS(SUD)+EPS))
	ID=MAX(ID1,ID2)
	RETURN
	END


	SUBROUTINE SPHY(N,X,NM,SY,DY)
C
C       ======================================================
C       Purpose: Compute spherical Bessel functions yn(x) and
C                their derivatives
C       Input :  x --- Argument of yn(x) ( x � 0 )
C                n --- Order of yn(x) ( n = 0,1,��� )
C       Output:  SY(n) --- yn(x)
C                DY(n) --- yn'(x)
C                NM --- Highest order computed
C       ======================================================
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION SY(0:N),DY(0:N)
	NM=N
	IF (X.LT.1.0D-60) THEN
	   DO 10 K=0,N
	      SY(K)=-1.0D+300
10            DY(K)=1.0D+300
	   RETURN
	ENDIF
	SY(0)=-DCOS(X)/X
	SY(1)=(SY(0)-DSIN(X))/X
	F0=SY(0)
	F1=SY(1)
	DO 15 K=2,N
	   F=(2.0D0*K-1.0D0)*F1/X-F0
	   SY(K)=F
	   IF (DABS(F).GE.1.0D+300) GO TO 20              
	   F0=F1
15         F1=F
20      NM=K-1
	   DY(0)=(DSIN(X)+DCOS(X)/X)/X
	   DO 25 K=1,NM
25            DY(K)=SY(K-1)-(K+1.0D0)*SY(K)/X
	RETURN
	END


	SUBROUTINE RMN2SP(M,N,C,X,CV,DF,KD,R2F,R2D)
C
C       ======================================================
C       Purpose: Compute prolate spheroidal radial function
C                of the second kind with a small argument
C       Routines called:
C            (1) LPMNS for computing the associated Legendre
C                functions of the first kind    
C            (2) LQMNS for computing the associated Legendre
C                functions of the second kind  
C            (3) KMN for computing expansion coefficients
C                and joining factors
C       ======================================================
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION PM(0:251),PD(0:251),QM(0:251),QD(0:251),
     &            DN(200),DF(200)
	IF (DABS(DF(1)).LT.1.0D-280) THEN
	   R2F=1.0D+300
	   R2D=1.0D+300
	   RETURN
	ENDIF
	EPS=1.0D-14
	IP=1
	NM1=INT((N-M)/2)
	IF (N-M.EQ.2*NM1) IP=0
	NM=25+NM1+INT(C)
	NM2=2*NM+M
	CALL KMN(M,N,C,CV,KD,DF,DN,CK1,CK2)
	CALL LPMNS(M,NM2,X,PM,PD)
	CALL LQMNS(M,NM2,X,QM,QD)
	SU0=0.0D0
	DO 10 K=1,NM
	  J=2*K-2+M+IP
	  SU0=SU0+DF(K)*QM(J)
	  IF (K.GT.NM1.AND.DABS(SU0-SW).LT.DABS(SU0)*EPS) GO TO 15                                                                   
10        SW=SU0
15      SD0=0.0D0
	DO 20 K=1,NM
	  J=2*K-2+M+IP
	  SD0=SD0+DF(K)*QD(J)
	  IF (K.GT.NM1.AND.DABS(SD0-SW).LT.DABS(SD0)*EPS) GO TO 25
20        SW=SD0
25        SU1=0.0D0
	  SD1=0.0D0
	  DO 30 K=1,M
	     J=M-2*K+IP
	     IF (J.LT.0) J=-J-1
	     SU1=SU1+DN(K)*QM(J)
30           SD1=SD1+DN(K)*QD(J)
	  GA=((X-1.0D0)/(X+1.0D0))**(0.5D0*M)
	  DO 55 K=1,M
	     J=M-2*K+IP
	     IF (J.GE.0) GO TO 55
	     IF (J.LT.0) J=-J-1
	     R1=1.0D0
	     DO 35 J1=1,J
35              R1=(M+J1)*R1
	     R2=1.0D0
	     DO 40 J2=1,M-J-2
40              R2=J2*R2
	     R3=1.0D0
	     SF=1.0D0
	     DO 45 L1=1,J
		R3=0.5D0*R3*(-J+L1-1.0)*(J+L1)/((M+L1)*L1)*(1.0-X)
45              SF=SF+R3
	     IF (M-J.GE.2) GB=(M-J-1.0D0)*R2
	     IF (M-J.LE.1) GB=1.0D0
	     SPL=R1*GA*GB*SF
	     SU1=SU1+(-1)**(J+M)*DN(K)*SPL
	     SPD1=M/(X*X-1.0D0)*SPL
	     GC=0.5D0*J*(J+1.0)/(M+1.0)
	     SD=1.0D0
	     R4=1.0D0
	     DO 50 L1=1,J-1
		R4=0.5D0*R4*(-J+L1)*(J+L1+1.0)/((M+L1+1.0)*L1)
     &             *(1.0-X)
50              SD=SD+R4
	     SPD2=R1*GA*GB*GC*SD
	     SD1=SD1+(-1)**(J+M)*DN(K)*(SPD1+SPD2)
55        CONTINUE
	  SU2=0.0D0
	  KI=(2*M+1+IP)/2
	  NM3=NM+KI
	  DO 60 K=KI,NM3
	     J=2*K-1-M-IP
	     SU2=SU2+DN(K)*PM(J)
	     IF (J.GT.M.AND.DABS(SU2-SW).LT.DABS(SU2)*EPS) GO TO 65
60           SW=SU2
65        SD2=0.0D0
	  DO 70 K=KI,NM3
	     J=2*K-1-M-IP
	     SD2=SD2+DN(K)*PD(J)
	     IF (J.GT.M.AND.DABS(SD2-SW).LT.DABS(SD2)*EPS) GO TO 75
70           SW=SD2
75      SUM=SU0+SU1+SU2
	SDM=SD0+SD1+SD2
	R2F=SUM/CK2
	R2D=SDM/CK2
	RETURN
	END


	SUBROUTINE LPMNS(M,N,X,PM,PD)
C
C       ========================================================
C       Purpose: Compute associated Legendre functions Pmn(x)
C                and Pmn'(x) for a given order
C       Input :  x --- Argument of Pmn(x)
C                m --- Order of Pmn(x),  m = 0,1,2,...,n
C                n --- Degree of Pmn(x), n = 0,1,2,...,N
C       Output:  PM(n) --- Pmn(x)
C                PD(n) --- Pmn'(x)
C       ========================================================
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION PM(0:N),PD(0:N)
	DO 10 K=0,N
	   PM(K)=0.0D0
10         PD(K)=0.0D0
	IF (DABS(X).EQ.1.0D0) THEN
	   DO 15 K=0,N
	      IF (M.EQ.0) THEN
		 PM(K)=1.0D0
		 PD(K)=0.5D0*K*(K+1.0)
		 IF (X.LT.0.0) THEN
		    PM(K)=(-1)**K*PM(K)
		    PD(K)=(-1)**(K+1)*PD(K)
		 ENDIF
	      ELSE IF (M.EQ.1) THEN
		 PD(K)=1.0D+300
	      ELSE IF (M.EQ.2) THEN
		 PD(K)=-0.25D0*(K+2.0)*(K+1.0)*K*(K-1.0)
		 IF (X.LT.0.0) PD(K)=(-1)**(K+1)*PD(K)
	      ENDIF
15         CONTINUE
	   RETURN
	ENDIF
	X0=DABS(1.0D0-X*X)
	PM0=1.0D0
	PMK=PM0
	DO 20 K=1,M
	   PMK=(2.0D0*K-1.0D0)*DSQRT(X0)*PM0
20         PM0=PMK
	PM1=(2.0D0*M+1.0D0)*X*PM0
	PM(M)=PMK
	PM(M+1)=PM1
	DO 25 K=M+2,N
	   PM2=((2.0D0*K-1.0D0)*X*PM1-(K+M-1.0D0)*PMK)/(K-M)
	   PM(K)=PM2
	   PMK=PM1
25         PM1=PM2
	PD(0)=((1.0D0-M)*PM(1)-X*PM(0))/(X*X-1.0)  
	DO 30 K=1,N
30          PD(K)=(K*X*PM(K)-(K+M)*PM(K-1))/(X*X-1.0D0)
	RETURN
	END


	SUBROUTINE LQMNS(M,N,X,QM,QD)
C
C       ========================================================
C       Purpose: Compute associated Legendre functions Qmn(x)
C                and Qmn'(x) for a given order
C       Input :  x --- Argument of Qmn(x)
C                m --- Order of Qmn(x),  m = 0,1,2,...
C                n --- Degree of Qmn(x), n = 0,1,2,...
C       Output:  QM(n) --- Qmn(x)
C                QD(n) --- Qmn'(x)
C       ========================================================
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION QM(0:N),QD(0:N)
	DO 10 K=0,N
	   QM(K)=0.0D0
10         QD(K)=0.0D0
	IF (DABS(X).EQ.1.0D0) THEN
	   DO 15 K=0,N
	      QM(K)=1.0D+300
15            QD(K)=1.0D+300
	   RETURN
	ENDIF
	LS=1
	IF (DABS(X).GT.1.0D0) LS=-1
	XQ=DSQRT(LS*(1.0D0-X*X))
	Q0=0.5D0*DLOG(DABS((X+1.0)/(X-1.0)))
	Q00=Q0
	Q10=-1.0D0/XQ
	Q01=X*Q0-1.0D0
	Q11=-LS*XQ*(Q0+X/(1.0D0-X*X))
	QF0=Q00
	QF1=Q10
	DO 20 K=2,M
	   QM0=-2.0D0*(K-1.0)/XQ*X*QF1-LS*(K-1.0)*(2.0-K)*QF0
	   QF0=QF1
20         QF1=QM0
	IF (M.EQ.0) QM0=Q00
	IF (M.EQ.1) QM0=Q10
	QM(0)=QM0
	IF (DABS(X).LT.1.0001D0) THEN
	   IF (M.EQ.0.AND.N.GT.0) THEN
	      QF0=Q00
	      QF1=Q01
	      DO 25 K=2,N
		 QF2=((2.0*K-1.0D0)*X*QF1-(K-1.0)*QF0)/K
		 QM(K)=QF2
		 QF0=QF1
25               QF1=QF2
	   ENDIF
	   QG0=Q01
	   QG1=Q11
	   DO 30 K=2,M
	      QM1=-2.0D0*(K-1.0)/XQ*X*QG1-LS*K*(3.0-K)*QG0
	      QG0=QG1
30            QG1=QM1
	   IF (M.EQ.0) QM1=Q01
	   IF (M.EQ.1) QM1=Q11
	   QM(1)=QM1
	   IF (M.EQ.1.AND.N.GT.1) THEN
	      QH0=Q10
	      QH1=Q11
	      DO 35 K=2,N
		 QH2=((2.0*K-1.0D0)*X*QH1-K*QH0)/(K-1.0)
		 QM(K)=QH2
		 QH0=QH1
35               QH1=QH2
	   ELSE IF (M.GE.2) THEN
	      QG0=Q00
	      QG1=Q01
	      QH0=Q10
	      QH1=Q11
	      DO 45 L=2,N
		 Q0L=((2.0D0*L-1.0D0)*X*QG1-(L-1.0D0)*QG0)/L
		 Q1L=((2.0*L-1.0D0)*X*QH1-L*QH0)/(L-1.0D0)
		 QF0=Q0L
		 QF1=Q1L
		 DO 40 K=2,M
		    QMK=-2.0D0*(K-1.0)/XQ*X*QF1-LS*(K+L-1.0)*
     &                  (L+2.0-K)*QF0
		    QF0=QF1
40                  QF1=QMK
		 QM(L)=QMK
		 QG0=QG1
		 QG1=Q0L
		 QH0=QH1
45               QH1=Q1L
	   ENDIF
	ELSE
	   IF (DABS(X).GT.1.1) THEN
	      KM=40+M+N
	   ELSE
	      KM=(40+M+N)*INT(-1.0-1.8*LOG(X-1.0))
	   ENDIF
	   QF2=0.0D0
	   QF1=1.0D0
	   DO 50 K=KM,0,-1
	      QF0=((2.0*K+3.0D0)*X*QF1-(K+2.0-M)*QF2)/(K+M+1.0)
	      IF (K.LE.N) QM(K)=QF0
	      QF2=QF1
50            QF1=QF0
	   DO 55 K=0,N
55            QM(K)=QM(K)*QM0/QF0
	ENDIF
	IF (DABS(X).LT.1.0D0) THEN
	   DO 60 K=0,N
60            QM(K)=(-1)**M*QM(K)
	ENDIF
	QD(0)=((1.0D0-M)*QM(1)-X*QM(0))/(X*X-1.0)
	DO 65 K=1,N
65         QD(K)=(K*X*QM(K)-(K+M)*QM(K-1))/(X*X-1.0)
	RETURN
	END


	SUBROUTINE KMN(M,N,C,CV,KD,DF,DN,CK1,CK2)
C
C       ===================================================
C       Purpose: Compute the expansion coefficients of the
C                prolate and oblate spheroidal functions
C                and joining factors
C       ===================================================
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION U(200),V(200),W(200),DF(200),DN(200),
     &            TP(200),RK(200)
	NM=25+INT(0.5*(N-M)+C)
	NN=NM+M
	CS=C*C*KD
	IP=1
	IF (N-M.EQ.2*INT((N-M)/2)) IP=0
	DO 10 I=1,NN+3     
	   IF (IP.EQ.0) K=-2*(I-1)
	   IF (IP.EQ.1) K=-(2*I-3)
	   GK0=2.0D0*M+K
	   GK1=(M+K)*(M+K+1.0D0)
	   GK2=2.0D0*(M+K)-1.0D0
	   GK3=2.0D0*(M+K)+3.0D0
	   U(I)=GK0*(GK0-1.0D0)*CS/(GK2*(GK2+2.0D0))
	   V(I)=GK1-CV+(2.0D0*(GK1-M*M)-1.0D0)*CS/(GK2*GK3)
10         W(I)=(K+1.0D0)*(K+2.0D0)*CS/((GK2+2.0D0)*GK3)
	DO 20 K=1,M
	   T=V(M+1)
	   DO 15 L=0,M-K-1
15            T=V(M-L)-W(M-L+1)*U(M-L)/T
20         RK(K)=-U(K)/T
	R=1.0D0
	DO 25 K=1,M
	   R=R*RK(K)
25         DN(K)=DF(1)*R
	TP(NN)=V(NN+1)
	DO 30 K=NN-1,M+1,-1
	   TP(K)=V(K+1)-W(K+2)*U(K+1)/TP(K+1)
	   IF (K.GT.M+1) RK(K)=-U(K)/TP(K)
30      CONTINUE
	IF (M.EQ.0) DNP=DF(1)
	IF (M.NE.0) DNP=DN(M)
	DN(M+1)=(-1)**IP*DNP*CS/((2.0*M-1.0)*(2.0*M+1.0-4.0*IP)
     &          *TP(M+1))
	DO 35 K=M+2,NN
35         DN(K)=RK(K)*DN(K-1)
	R1=1.0D0
	DO 40 J=1,(N+M+IP)/2
40         R1=R1*(J+0.5D0*(N+M+IP))
	NM1=(N-M)/2
	R=1.0D0
	DO 45 J=1,2*M+IP
45         R=R*J
	SU0=R*DF(1)
	DO 50 K=2,NM
	   R=R*(M+K-1.0)*(M+K+IP-1.5D0)/(K-1.0D0)/(K+IP-1.5D0)
	   SU0=SU0+R*DF(K)
	   IF (K.GT.NM1.AND.DABS((SU0-SW)/SU0).LT.1.0D-14) GO TO 55
50         SW=SU0
55      IF (KD.EQ.1) GOTO 70
	R2=1.0D0
	DO 60 J=1,M
60         R2=2.0D0*C*R2*J
	R3=1.0D0
	DO 65 J=1,(N-M-IP)/2
65         R3=R3*J
	SA0=(2.0*(M+IP)+1.0)*R1/(2.0**N*C**IP*R2*R3*DF(1))
	CK1=SA0*SU0
	IF (KD.EQ.-1) RETURN
70      R4=1.0D0
	DO 75 J=1,(N-M-IP)/2
75         R4=4.0D0*R4*J
	R5=1.0D0
	DO 80 J=1,M
80         R5=R5*(J+M)/C
	G0=DN(M)
	IF (M.EQ.0) G0=DF(1)
	SB0=(IP+1.0)*C**(IP+1)/(2.0*IP*(M-2.0)+1.0)/(2.0*M-1.0)
	CK2=(-1)**IP*SB0*R4*R5*G0/R1*SU0
	RETURN
	END


	SUBROUTINE SEGV(M,N,C,KD,CV,EG)
C
C       =========================================================
C       Purpose: Compute the characteristic values of spheroidal
C                wave functions
C       Input :  m  --- Mode parameter
C                n  --- Mode parameter
C                c  --- Spheroidal parameter
C                KD --- Function code
C                       KD=1 for Prolate; KD=-1 for Oblate
C       Output:  CV --- Characteristic value for given m, n and c
C                EG(L) --- Characteristic value for mode m and n'
C                          ( L = n' - m + 1 )
C       =========================================================
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION B(100),H(100),D(300),E(300),F(300),CV0(100),
     &            A(300),G(300),EG(200)
	IF (C.LT.1.0D-10) THEN
	   DO 5 I=1,N
5             EG(I)=(I+M)*(I+M-1.0D0)
	   GO TO 70
	ENDIF                                           
	ICM=(N-M+2)/2
	NM=10+INT(0.5*(N-M)+C)
	CS=C*C*KD
	DO 60 L=0,1
	   DO 10 I=1,NM
	      IF (L.EQ.0) K=2*(I-1)
	      IF (L.EQ.1) K=2*I-1
	      DK0=M+K
	      DK1=M+K+1
	      DK2=2*(M+K)
	      D2K=2*M+K
	      A(I)=(D2K+2.0)*(D2K+1.0)/((DK2+3.0)*(DK2+5.0))*CS
	      D(I)=DK0*DK1+(2.0*DK0*DK1-2.0*M*M-1.0)/((DK2-1.0)
     &             *(DK2+3.0))*CS
10            G(I)=K*(K-1.0)/((DK2-3.0)*(DK2-1.0))*CS
	   DO 15 K=2,NM
	      E(K)=DSQRT(A(K-1)*G(K))
15            F(K)=E(K)*E(K)
	   F(1)=0.0D0
	   E(1)=0.0D0
	   XA=D(NM)+DABS(E(NM))
	   XB=D(NM)-DABS(E(NM))
	   NM1=NM-1
	   DO 20 I=1,NM1
	      T=DABS(E(I))+DABS(E(I+1))
	      T1=D(I)+T
	      IF (XA.LT.T1) XA=T1
	      T1=D(I)-T
	      IF (T1.LT.XB) XB=T1
20         CONTINUE
	   DO 25 I=1,ICM
	      B(I)=XA
25            H(I)=XB
	   DO 55 K=1,ICM
	      DO 30 K1=K,ICM
		 IF (B(K1).LT.B(K)) THEN
		    B(K)=B(K1)
		    GO TO 35
		 ENDIF
30            CONTINUE
35            IF (K.NE.1.AND.H(K).LT.H(K-1)) H(K)=H(K-1)
40            X1=(B(K)+H(K))/2.0D0
	      CV0(K)=X1
	      IF (DABS((B(K)-H(K))/X1).LT.1.0D-14) GO TO 50
	      J=0
	      S=1.0D0
	      DO 45 I=1,NM
		 IF (S.EQ.0.0D0) S=S+1.0D-30
		 T=F(I)/S
		 S=D(I)-T-X1
		 IF (S.LT.0.0D0) J=J+1
45            CONTINUE
	      IF (J.LT.K) THEN
		 H(K)=X1
	      ELSE
		 B(K)=X1
		 IF (J.GE.ICM) THEN
		    B(ICM)=X1
		 ELSE
		    IF (H(J+1).LT.X1) H(J+1)=X1
		    IF (X1.LT.B(J)) B(J)=X1
		 ENDIF
	      ENDIF
	      GO TO 40
50            CV0(K)=X1
	      IF (L.EQ.0) EG(2*K-1)=CV0(K)
	      IF (L.EQ.1) EG(2*K)=CV0(K)
55         CONTINUE
60      CONTINUE
70      CV=EG(N-M+1)
	RETURN
	END
