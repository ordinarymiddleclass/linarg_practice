*DECK DBESJ
      FUNCTION DLNGAM(xx)
      REAL*8 DLNGAM,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      DLNGAM=tmp+log(stp*ser/x)
      RETURN
      END
*DECK DBESJ
      SUBROUTINE DBESJ (X, ALPHA, N, Y, NZ)
C***BEGIN PROLOGUE  DBESJ
C***PURPOSE  Compute an N member sequence of J Bessel functions
C            J_{ALPHA+K-1}(X), K=1,...,N for non-negative ALPHA
C            and X.
C***LIBRARY   SLATEC
C***CATEGORY  C10A3
C***TYPE      DOUBLE PRECISION (BESJ-S, DBESJ-D)
C***KEYWORDS  J BESSEL FUNCTION, SPECIAL FUNCTIONS
C***AUTHOR  Amos, D. E., (SNLA)
C           Daniel, S. L., (SNLA)
C           Weston, M. K., (SNLA)
C***DESCRIPTION
C
C     Abstract  **** a double precision routine ****
C         DBESJ computes an N member sequence of J Bessel functions
C         J_{ALPHA+K-1}(X), K=1,...,N for non-negative ALPHA and X.
C         A combination of the power series, the asymptotic expansion
C         for X to infinity and the uniform asymptotic expansion for
C         NU to infinity are applied over subdivisions of the (NU,X)
C         plane.  For values of (NU,X) not covered by one of these
C         formulae, the order is incremented or decremented by integer
C         values into a region where one of the formulae apply. Backward
C         recursion is applied to reduce orders by integer values except
C         where the entire sequence lies in the oscillatory region.  In
C         this case forward recursion is stable and values from the
C         asymptotic expansion for X to infinity start the recursion
C         when it is efficient to do so. Leading terms of the series and
C         uniform expansion are tested for underflow.  If a sequence is
C         requested and the last member would underflow, the result is
C         set to zero and the next lower order tried, etc., until a
C         member comes on scale or all members are set to zero.
C         Overflow cannot occur.
C
C         The maximum number of significant digits obtainable
C         is the smaller of 14 and the number of digits carried in
C         double precision arithmetic.
C
C     Description of Arguments
C
C         Input      X,ALPHA are double precision
C           X      - X .GE. 0.0D0
C           ALPHA  - order of first member of the sequence,
C                    ALPHA .GE. 0.0D0
C           N      - number of members in the sequence, N .GE. 1
C
C         Output     Y is double precision
C           Y      - a vector whose first N components contain
C                    values for J_{ALPHA+K-1}(X), K=1,...,N
C           NZ     - number of components of Y set to zero due to
C                    underflow,
C                    NZ=0   , normal return, computation completed
C                    NZ .NE. 0, last NZ components of Y set to zero,
C                             Y(K)=0.0D0, K=N-NZ+1,...,N.
C
C     Error Conditions
C         Improper input arguments - a fatal error
C         Underflow  - a non-fatal error (NZ .NE. 0)
C
C***REFERENCES  D. E. Amos, S. L. Daniel and M. K. Weston, CDC 6600
C                 subroutines IBESS and JBESS for Bessel functions
C                 I(NU,X) and J(NU,X), X .GE. 0, NU .GE. 0, ACM
C                 Transactions on Mathematical Software 3, (1977),
C                 pp. 76-92.
C               F. W. J. Olver, Tables of Bessel Functions of Moderate
C                 or Large Orders, NPL Mathematical Tables 6, Her
C                 Majesty's Stationery Office, London, 1962.
C***ROUTINES CALLED  D1MACH, DASYJY, DJAIRY, DLNGAM, I1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   750101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DBESJ
      EXTERNAL DJAIRY
      INTEGER I,IALP,IDALP,IFLW,IN,INLIM,IS,I1,I2,K,KK,KM,KT,N,NN,
     1        NS,NZ
      INTEGER I1MACH
      DOUBLE PRECISION AK,AKM,ALPHA,ANS,AP,ARG,COEF,DALPHA,DFN,DTM,
     1           EARG,ELIM1,ETX,FIDAL,FLGJY,FN,FNF,FNI,FNP1,FNU,
     2           FNULIM,GLN,PDF,PIDT,PP,RDEN,RELB,RTTP,RTWO,RTX,RZDEN,
     3           S,SA,SB,SXO2,S1,S2,T,TA,TAU,TB,TEMP,TFN,TM,TOL,
     4           TOLLN,TRX,TX,T1,T2,WK,X,XO2,XO2L,Y,SLIM,RTOL
      SAVE RTWO, PDF, RTTP, PIDT, PP, INLIM, FNULIM
      DOUBLE PRECISION D1MACH, DLNGAM
      DIMENSION Y(*), TEMP(3), FNULIM(2), PP(4), WK(7)
      DATA RTWO,PDF,RTTP,PIDT                    / 1.34839972492648D+00,
     1 7.85398163397448D-01, 7.97884560802865D-01, 1.57079632679490D+00/
      DATA  PP(1),  PP(2),  PP(3),  PP(4)        / 8.72909153935547D+00,
     1 2.65693932265030D-01, 1.24578576865586D-01, 7.70133747430388D-04/
      DATA INLIM           /      150            /
      DATA FNULIM(1), FNULIM(2) /      100.0D0,     60.0D0     /
C***FIRST EXECUTABLE STATEMENT  DBESJ
      NZ = 0
      KT = 1
      NS=0
C     I1MACH(14) REPLACES I1MACH(11) IN A DOUBLE PRECISION CODE
C     I1MACH(15) REPLACES I1MACH(12) IN A DOUBLE PRECISION CODE
      TA = D1MACH(3)
      TOL = MAX(TA,1.0D-15)
      I1 = I1MACH(14) + 1
      I2 = I1MACH(15)
      TB = D1MACH(5)
      ELIM1 = -2.303D0*(I2*TB+3.0D0)
      RTOL=1.0D0/TOL
      SLIM=D1MACH(1)*RTOL*1.0D+3
C     TOLLN = -LN(TOL)
      TOLLN = 2.303D0*TB*I1
      TOLLN = MIN(TOLLN,34.5388D0)
      IF (N-1) 720, 10, 20
   10 KT = 2
   20 NN = N
      IF (X) 730, 30, 80
   30 IF (ALPHA) 710, 40, 50
   40 Y(1) = 1.0D0
      IF (N.EQ.1) RETURN
      I1 = 2
      GO TO 60
   50 I1 = 1
   60 DO 70 I=I1,N
        Y(I) = 0.0D0
   70 CONTINUE
      RETURN
   80 CONTINUE
      IF (ALPHA.LT.0.0D0) GO TO 710
C
      IALP = INT(ALPHA)
      FNI = IALP + N - 1
      FNF = ALPHA - IALP
      DFN = FNI + FNF
      FNU = DFN
      XO2 = X*0.5D0
      SXO2 = XO2*XO2
C
C     DECISION TREE FOR REGION WHERE SERIES, ASYMPTOTIC EXPANSION FOR X
C     TO INFINITY AND ASYMPTOTIC EXPANSION FOR NU TO INFINITY ARE
C     APPLIED.
C
      IF (SXO2.LE.(FNU+1.0D0)) GO TO 90
      TA = MAX(20.0D0,FNU)
      IF (X.GT.TA) GO TO 120
      IF (X.GT.12.0D0) GO TO 110
      XO2L = LOG(XO2)
      NS = INT(SXO2-FNU) + 1
      GO TO 100
   90 FN = FNU
      FNP1 = FN + 1.0D0
      XO2L = LOG(XO2)
      IS = KT
      IF (X.LE.0.50D0) GO TO 330
      NS = 0
  100 FNI = FNI + NS
      DFN = FNI + FNF
      FN = DFN
      FNP1 = FN + 1.0D0
      IS = KT
      IF (N-1+NS.GT.0) IS = 3
      GO TO 330
  110 ANS = MAX(36.0D0-FNU,0.0D0)
      NS = INT(ANS)
      FNI = FNI + NS
      DFN = FNI + FNF
      FN = DFN
      IS = KT
      IF (N-1+NS.GT.0) IS = 3
      GO TO 130
  120 CONTINUE
      RTX = SQRT(X)
      TAU = RTWO*RTX
      TA = TAU + FNULIM(KT)
      IF (FNU.LE.TA) GO TO 480
      FN = FNU
      IS = KT
C
C     UNIFORM ASYMPTOTIC EXPANSION FOR NU TO INFINITY
C
  130 CONTINUE
      I1 = ABS(3-IS)
      I1 = MAX(I1,1)
      FLGJY = 1.0D0
      CALL DASYJY(DJAIRY,X,FN,FLGJY,I1,TEMP(IS),WK,IFLW)
      IF(IFLW.NE.0) GO TO 380
      GO TO (320, 450, 620), IS
  310 TEMP(1) = TEMP(3)
      KT = 1
  320 IS = 2
      FNI = FNI - 1.0D0
      DFN = FNI + FNF
      FN = DFN
      IF(I1.EQ.2) GO TO 450
      GO TO 130
C
C     SERIES FOR (X/2)**2.LE.NU+1
C
  330 CONTINUE
      GLN = DLNGAM(FNP1)
      ARG = FN*XO2L - GLN
      IF (ARG.LT.(-ELIM1)) GO TO 400
      EARG = EXP(ARG)
  340 CONTINUE
      S = 1.0D0
      IF (X.LT.TOL) GO TO 360
      AK = 3.0D0
      T2 = 1.0D0
      T = 1.0D0
      S1 = FN
      DO 350 K=1,17
        S2 = T2 + S1
        T = -T*SXO2/S2
        S = S + T
        IF (ABS(T).LT.TOL) GO TO 360
        T2 = T2 + AK
        AK = AK + 2.0D0
        S1 = S1 + FN
  350 CONTINUE
  360 CONTINUE
      TEMP(IS) = S*EARG
      GO TO (370, 450, 610), IS
  370 EARG = EARG*FN/XO2
      FNI = FNI - 1.0D0
      DFN = FNI + FNF
      FN = DFN
      IS = 2
      GO TO 340
C
C     SET UNDERFLOW VALUE AND UPDATE PARAMETERS
C     UNDERFLOW CAN ONLY OCCUR FOR NS=0 SINCE THE ORDER MUST BE LARGER
C     THAN 36. THEREFORE, NS NEE NOT BE TESTED.
C
  380 Y(NN) = 0.0D0
      NN = NN - 1
      FNI = FNI - 1.0D0
      DFN = FNI + FNF
      FN = DFN
      IF (NN-1) 440, 390, 130
  390 KT = 2
      IS = 2
      GO TO 130
  400 Y(NN) = 0.0D0
      NN = NN - 1
      FNP1 = FN
      FNI = FNI - 1.0D0
      DFN = FNI + FNF
      FN = DFN
      IF (NN-1) 440, 410, 420
  410 KT = 2
      IS = 2
  420 IF (SXO2.LE.FNP1) GO TO 430
      GO TO 130
  430 ARG = ARG - XO2L + LOG(FNP1)
      IF (ARG.LT.(-ELIM1)) GO TO 400
      GO TO 330
  440 NZ = N - NN
      RETURN
C
C     BACKWARD RECURSION SECTION
C
  450 CONTINUE
      IF(NS.NE.0) GO TO 451
      NZ = N - NN
      IF (KT.EQ.2) GO TO 470
C     BACKWARD RECUR FROM INDEX ALPHA+NN-1 TO ALPHA
      Y(NN) = TEMP(1)
      Y(NN-1) = TEMP(2)
      IF (NN.EQ.2) RETURN
  451 CONTINUE
      TRX = 2.0D0/X
      DTM = FNI
      TM = (DTM+FNF)*TRX
      AK=1.0D0
      TA=TEMP(1)
      TB=TEMP(2)
      IF(ABS(TA).GT.SLIM) GO TO 455
      TA=TA*RTOL
      TB=TB*RTOL
      AK=TOL
  455 CONTINUE
      KK=2
      IN=NS-1
      IF(IN.EQ.0) GO TO 690
      IF(NS.NE.0) GO TO 670
      K=NN-2
      DO 460 I=3,NN
        S=TB
        TB = TM*TB - TA
        TA=S
        Y(K)=TB*AK
        DTM = DTM - 1.0D0
        TM = (DTM+FNF)*TRX
        K = K - 1
  460 CONTINUE
      RETURN
  470 Y(1) = TEMP(2)
      RETURN
C
C     ASYMPTOTIC EXPANSION FOR X TO INFINITY WITH FORWARD RECURSION IN
C     OSCILLATORY REGION X.GT.MAX(20, NU), PROVIDED THE LAST MEMBER
C     OF THE SEQUENCE IS ALSO IN THE REGION.
C
  480 CONTINUE
      IN = INT(ALPHA-TAU+2.0D0)
      IF (IN.LE.0) GO TO 490
      IDALP = IALP - IN - 1
      KT = 1
      GO TO 500
  490 CONTINUE
      IDALP = IALP
      IN = 0
  500 IS = KT
      FIDAL = IDALP
      DALPHA = FIDAL + FNF
      ARG = X - PIDT*DALPHA - PDF
      SA = SIN(ARG)
      SB = COS(ARG)
      COEF = RTTP/RTX
      ETX = 8.0D0*X
  510 CONTINUE
      DTM = FIDAL + FIDAL
      DTM = DTM*DTM
      TM = 0.0D0
      IF (FIDAL.EQ.0.0D0 .AND. ABS(FNF).LT.TOL) GO TO 520
      TM = 4.0D0*FNF*(FIDAL+FIDAL+FNF)
  520 CONTINUE
      TRX = DTM - 1.0D0
      T2 = (TRX+TM)/ETX
      S2 = T2
      RELB = TOL*ABS(T2)
      T1 = ETX
      S1 = 1.0D0
      FN = 1.0D0
      AK = 8.0D0
      DO 530 K=1,13
        T1 = T1 + ETX
        FN = FN + AK
        TRX = DTM - FN
        AP = TRX + TM
        T2 = -T2*AP/T1
        S1 = S1 + T2
        T1 = T1 + ETX
        AK = AK + 8.0D0
        FN = FN + AK
        TRX = DTM - FN
        AP = TRX + TM
        T2 = T2*AP/T1
        S2 = S2 + T2
        IF (ABS(T2).LE.RELB) GO TO 540
        AK = AK + 8.0D0
  530 CONTINUE
  540 TEMP(IS) = COEF*(S1*SB-S2*SA)
      IF(IS.EQ.2) GO TO 560
      FIDAL = FIDAL + 1.0D0
      DALPHA = FIDAL + FNF
      IS = 2
      TB = SA
      SA = -SB
      SB = TB
      GO TO 510
C
C     FORWARD RECURSION SECTION
C
  560 IF (KT.EQ.2) GO TO 470
      S1 = TEMP(1)
      S2 = TEMP(2)
      TX = 2.0D0/X
      TM = DALPHA*TX
      IF (IN.EQ.0) GO TO 580
C
C     FORWARD RECUR TO INDEX ALPHA
C
      DO 570 I=1,IN
        S = S2
        S2 = TM*S2 - S1
        TM = TM + TX
        S1 = S
  570 CONTINUE
      IF (NN.EQ.1) GO TO 600
      S = S2
      S2 = TM*S2 - S1
      TM = TM + TX
      S1 = S
  580 CONTINUE
C
C     FORWARD RECUR FROM INDEX ALPHA TO ALPHA+N-1
C
      Y(1) = S1
      Y(2) = S2
      IF (NN.EQ.2) RETURN
      DO 590 I=3,NN
        Y(I) = TM*Y(I-1) - Y(I-2)
        TM = TM + TX
  590 CONTINUE
      RETURN
  600 Y(1) = S2
      RETURN
C
C     BACKWARD RECURSION WITH NORMALIZATION BY
C     ASYMPTOTIC EXPANSION FOR NU TO INFINITY OR POWER SERIES.
C
  610 CONTINUE
C     COMPUTATION OF LAST ORDER FOR SERIES NORMALIZATION
      AKM = MAX(3.0D0-FN,0.0D0)
      KM = INT(AKM)
      TFN = FN + KM
      TA = (GLN+TFN-0.9189385332D0-0.0833333333D0/TFN)/(TFN+0.5D0)
      TA = XO2L - TA
      TB = -(1.0D0-1.5D0/TFN)/TFN
      AKM = TOLLN/(-TA+SQRT(TA*TA-TOLLN*TB)) + 1.5D0
      IN = KM + INT(AKM)
      GO TO 660
  620 CONTINUE
C     COMPUTATION OF LAST ORDER FOR ASYMPTOTIC EXPANSION NORMALIZATION
      GLN = WK(3) + WK(2)
      IF (WK(6).GT.30.0D0) GO TO 640
      RDEN = (PP(4)*WK(6)+PP(3))*WK(6) + 1.0D0
      RZDEN = PP(1) + PP(2)*WK(6)
      TA = RZDEN/RDEN
      IF (WK(1).LT.0.10D0) GO TO 630
      TB = GLN/WK(5)
      GO TO 650
  630 TB=(1.259921049D0+(0.1679894730D0+0.0887944358D0*WK(1))*WK(1))
     1 /WK(7)
      GO TO 650
  640 CONTINUE
      TA = 0.5D0*TOLLN/WK(4)
      TA=((0.0493827160D0*TA-0.1111111111D0)*TA+0.6666666667D0)*TA*WK(6)
      IF (WK(1).LT.0.10D0) GO TO 630
      TB = GLN/WK(5)
  650 IN = INT(TA/TB+1.5D0)
      IF (IN.GT.INLIM) GO TO 310
  660 CONTINUE
      DTM = FNI + IN
      TRX = 2.0D0/X
      TM = (DTM+FNF)*TRX
      TA = 0.0D0
      TB = TOL
      KK = 1
      AK=1.0D0
  670 CONTINUE
C
C     BACKWARD RECUR UNINDEXED
C
      DO 680 I=1,IN
        S = TB
        TB = TM*TB - TA
        TA = S
        DTM = DTM - 1.0D0
        TM = (DTM+FNF)*TRX
  680 CONTINUE
C     NORMALIZATION
      IF (KK.NE.1) GO TO 690
      S=TEMP(3)
      SA=TA/TB
      TA=S
      TB=S
      IF(ABS(S).GT.SLIM) GO TO 685
      TA=TA*RTOL
      TB=TB*RTOL
      AK=TOL
  685 CONTINUE
      TA=TA*SA
      KK = 2
      IN = NS
      IF (NS.NE.0) GO TO 670
  690 Y(NN) = TB*AK
      NZ = N - NN
      IF (NN.EQ.1) RETURN
      K = NN - 1
      S=TB
      TB = TM*TB - TA
      TA=S
      Y(K)=TB*AK
      IF (NN.EQ.2) RETURN
      DTM = DTM - 1.0D0
      TM = (DTM+FNF)*TRX
      K=NN-2
C
C     BACKWARD RECUR INDEXED
C
      DO 700 I=3,NN
        S=TB
        TB = TM*TB - TA
        TA=S
        Y(K)=TB*AK
        DTM = DTM - 1.0D0
        TM = (DTM+FNF)*TRX
        K = K - 1
  700 CONTINUE
      RETURN
C
C
C
  710 CONTINUE
      CALL XERMSG ('SLATEC', 'DBESJ', 'ORDER, ALPHA, LESS THAN ZERO.',
     +   2, 1)
      RETURN
  720 CONTINUE
      CALL XERMSG ('SLATEC', 'DBESJ', 'N LESS THAN ONE.', 2, 1)
      RETURN
  730 CONTINUE
      CALL XERMSG ('SLATEC', 'DBESJ', 'X LESS THAN ZERO.', 2, 1)
      RETURN
      END

*DECK DJAIRY
      SUBROUTINE DJAIRY (X, RX, C, AI, DAI)
C***BEGIN PROLOGUE  DJAIRY
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DBESJ and DBESY
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (JAIRY-S, DJAIRY-D)
C***AUTHOR  Amos, D. E., (SNLA)
C           Daniel, S. L., (SNLA)
C           Weston, M. K., (SNLA)
C***DESCRIPTION
C
C                  DJAIRY computes the Airy function AI(X)
C                   and its derivative DAI(X) for DASYJY
C
C                                   INPUT
C
C         X - Argument, computed by DASYJY, X unrestricted
C        RX - RX=SQRT(ABS(X)), computed by DASYJY
C         C - C=2.*(ABS(X)**1.5)/3., computed by DASYJY
C
C                                  OUTPUT
C
C        AI - Value of function AI(X)
C       DAI - Value of the derivative DAI(X)
C
C***SEE ALSO  DBESJ, DBESY
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   750101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891009  Removed unreferenced variable.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910408  Updated the AUTHOR section.  (WRB)
C***END PROLOGUE  DJAIRY
C
      INTEGER I, J, M1, M1D, M2, M2D, M3, M3D, M4, M4D, N1, N1D, N2,
     1 N2D, N3, N3D, N4, N4D
      DOUBLE PRECISION A,AI,AJN,AJP,AK1,AK2,AK3,B,C,CCV,CON2,
     1 CON3, CON4, CON5, CV, DA, DAI, DAJN, DAJP, DAK1, DAK2, DAK3,
     2 DB, EC, E1, E2, FPI12, F1, F2, RTRX, RX, SCV, T, TEMP1, TEMP2,
     3 TT, X
      DIMENSION AJP(19), AJN(19), A(15), B(15)
      DIMENSION AK1(14), AK2(23), AK3(14)
      DIMENSION DAJP(19), DAJN(19), DA(15), DB(15)
      DIMENSION DAK1(14), DAK2(24), DAK3(14)
      SAVE N1, N2, N3, N4, M1, M2, M3, M4, FPI12, CON2, CON3,
     1 CON4, CON5, AK1, AK2, AK3, AJP, AJN, A, B,
     2 N1D, N2D, N3D, N4D, M1D, M2D, M3D, M4D, DAK1, DAK2, DAK3,
     3 DAJP, DAJN, DA, DB
      DATA N1,N2,N3,N4/14,23,19,15/
      DATA M1,M2,M3,M4/12,21,17,13/
      DATA FPI12,CON2,CON3,CON4,CON5/
     1 1.30899693899575D+00, 5.03154716196777D+00, 3.80004589867293D-01,
     2 8.33333333333333D-01, 8.66025403784439D-01/
      DATA AK1(1), AK1(2), AK1(3), AK1(4), AK1(5), AK1(6), AK1(7),
     1     AK1(8), AK1(9), AK1(10),AK1(11),AK1(12),AK1(13),
     2     AK1(14)         / 2.20423090987793D-01,-1.25290242787700D-01,
     3 1.03881163359194D-02, 8.22844152006343D-04,-2.34614345891226D-04,
     4 1.63824280172116D-05, 3.06902589573189D-07,-1.29621999359332D-07,
     5 8.22908158823668D-09, 1.53963968623298D-11,-3.39165465615682D-11,
     6 2.03253257423626D-12,-1.10679546097884D-14,-5.16169497785080D-15/
      DATA AK2(1), AK2(2), AK2(3), AK2(4), AK2(5), AK2(6), AK2(7),
     1     AK2(8), AK2(9), AK2(10),AK2(11),AK2(12),AK2(13),AK2(14),
     2     AK2(15),AK2(16),AK2(17),AK2(18),AK2(19),AK2(20),AK2(21),
     3     AK2(22),AK2(23) / 2.74366150869598D-01, 5.39790969736903D-03,
     4-1.57339220621190D-03, 4.27427528248750D-04,-1.12124917399925D-04,
     5 2.88763171318904D-05,-7.36804225370554D-06, 1.87290209741024D-06,
     6-4.75892793962291D-07, 1.21130416955909D-07,-3.09245374270614D-08,
     7 7.92454705282654D-09,-2.03902447167914D-09, 5.26863056595742D-10,
     8-1.36704767639569D-10, 3.56141039013708D-11,-9.31388296548430D-12,
     9 2.44464450473635D-12,-6.43840261990955D-13, 1.70106030559349D-13,
     1-4.50760104503281D-14, 1.19774799164811D-14,-3.19077040865066D-15/
      DATA AK3(1), AK3(2), AK3(3), AK3(4), AK3(5), AK3(6), AK3(7),
     1     AK3(8), AK3(9), AK3(10),AK3(11),AK3(12),AK3(13),
     2     AK3(14)         / 2.80271447340791D-01,-1.78127042844379D-03,
     3 4.03422579628999D-05,-1.63249965269003D-06, 9.21181482476768D-08,
     4-6.52294330229155D-09, 5.47138404576546D-10,-5.24408251800260D-11,
     5 5.60477904117209D-12,-6.56375244639313D-13, 8.31285761966247D-14,
     6-1.12705134691063D-14, 1.62267976598129D-15,-2.46480324312426D-16/
      DATA AJP(1), AJP(2), AJP(3), AJP(4), AJP(5), AJP(6), AJP(7),
     1     AJP(8), AJP(9), AJP(10),AJP(11),AJP(12),AJP(13),AJP(14),
     2     AJP(15),AJP(16),AJP(17),AJP(18),
     3     AJP(19)         / 7.78952966437581D-02,-1.84356363456801D-01,
     4 3.01412605216174D-02, 3.05342724277608D-02,-4.95424702513079D-03,
     5-1.72749552563952D-03, 2.43137637839190D-04, 5.04564777517082D-05,
     6-6.16316582695208D-06,-9.03986745510768D-07, 9.70243778355884D-08,
     7 1.09639453305205D-08,-1.04716330588766D-09,-9.60359441344646D-11,
     8 8.25358789454134D-12, 6.36123439018768D-13,-4.96629614116015D-14,
     9-3.29810288929615D-15, 2.35798252031104D-16/
      DATA AJN(1), AJN(2), AJN(3), AJN(4), AJN(5), AJN(6), AJN(7),
     1     AJN(8), AJN(9), AJN(10),AJN(11),AJN(12),AJN(13),AJN(14),
     2     AJN(15),AJN(16),AJN(17),AJN(18),
     3     AJN(19)         / 3.80497887617242D-02,-2.45319541845546D-01,
     4 1.65820623702696D-01, 7.49330045818789D-02,-2.63476288106641D-02,
     5-5.92535597304981D-03, 1.44744409589804D-03, 2.18311831322215D-04,
     6-4.10662077680304D-05,-4.66874994171766D-06, 7.15218807277160D-07,
     7 6.52964770854633D-08,-8.44284027565946D-09,-6.44186158976978D-10,
     8 7.20802286505285D-11, 4.72465431717846D-12,-4.66022632547045D-13,
     9-2.67762710389189D-14, 2.36161316570019D-15/
      DATA A(1),   A(2),   A(3),   A(4),   A(5),   A(6),   A(7),
     1     A(8),   A(9),   A(10),  A(11),  A(12),  A(13),  A(14),
     2     A(15)           / 4.90275424742791D-01, 1.57647277946204D-03,
     3-9.66195963140306D-05, 1.35916080268815D-07, 2.98157342654859D-07,
     4-1.86824767559979D-08,-1.03685737667141D-09, 3.28660818434328D-10,
     5-2.57091410632780D-11,-2.32357655300677D-12, 9.57523279048255D-13,
     6-1.20340828049719D-13,-2.90907716770715D-15, 4.55656454580149D-15,
     7-9.99003874810259D-16/
      DATA B(1),   B(2),   B(3),   B(4),   B(5),   B(6),   B(7),
     1     B(8),   B(9),   B(10),  B(11),  B(12),  B(13),  B(14),
     2     B(15)           / 2.78593552803079D-01,-3.52915691882584D-03,
     3-2.31149677384994D-05, 4.71317842263560D-06,-1.12415907931333D-07,
     4-2.00100301184339D-08, 2.60948075302193D-09,-3.55098136101216D-11,
     5-3.50849978423875D-11, 5.83007187954202D-12,-2.04644828753326D-13,
     6-1.10529179476742D-13, 2.87724778038775D-14,-2.88205111009939D-15,
     7-3.32656311696166D-16/
      DATA N1D,N2D,N3D,N4D/14,24,19,15/
      DATA M1D,M2D,M3D,M4D/12,22,17,13/
      DATA DAK1(1), DAK1(2), DAK1(3), DAK1(4), DAK1(5), DAK1(6),
     1     DAK1(7), DAK1(8), DAK1(9), DAK1(10),DAK1(11),DAK1(12),
     2    DAK1(13),DAK1(14)/ 2.04567842307887D-01,-6.61322739905664D-02,
     3-8.49845800989287D-03, 3.12183491556289D-03,-2.70016489829432D-04,
     4-6.35636298679387D-06, 3.02397712409509D-06,-2.18311195330088D-07,
     5-5.36194289332826D-10, 1.13098035622310D-09,-7.43023834629073D-11,
     6 4.28804170826891D-13, 2.23810925754539D-13,-1.39140135641182D-14/
      DATA DAK2(1), DAK2(2), DAK2(3), DAK2(4), DAK2(5), DAK2(6),
     1     DAK2(7), DAK2(8), DAK2(9), DAK2(10),DAK2(11),DAK2(12),
     2     DAK2(13),DAK2(14),DAK2(15),DAK2(16),DAK2(17),DAK2(18),
     3     DAK2(19),DAK2(20),DAK2(21),DAK2(22),DAK2(23),
     4     DAK2(24)        / 2.93332343883230D-01,-8.06196784743112D-03,
     5 2.42540172333140D-03,-6.82297548850235D-04, 1.85786427751181D-04,
     6-4.97457447684059D-05, 1.32090681239497D-05,-3.49528240444943D-06,
     7 9.24362451078835D-07,-2.44732671521867D-07, 6.49307837648910D-08,
     8-1.72717621501538D-08, 4.60725763604656D-09,-1.23249055291550D-09,
     9 3.30620409488102D-10,-8.89252099772401D-11, 2.39773319878298D-11,
     1-6.48013921153450D-12, 1.75510132023731D-12,-4.76303829833637D-13,
     2 1.29498241100810D-13,-3.52679622210430D-14, 9.62005151585923D-15,
     3-2.62786914342292D-15/
      DATA DAK3(1), DAK3(2), DAK3(3), DAK3(4), DAK3(5), DAK3(6),
     1     DAK3(7), DAK3(8), DAK3(9), DAK3(10),DAK3(11),DAK3(12),
     2    DAK3(13),DAK3(14)/ 2.84675828811349D-01, 2.53073072619080D-03,
     3-4.83481130337976D-05, 1.84907283946343D-06,-1.01418491178576D-07,
     4 7.05925634457153D-09,-5.85325291400382D-10, 5.56357688831339D-11,
     5-5.90889094779500D-12, 6.88574353784436D-13,-8.68588256452194D-14,
     6 1.17374762617213D-14,-1.68523146510923D-15, 2.55374773097056D-16/
      DATA DAJP(1), DAJP(2), DAJP(3), DAJP(4), DAJP(5), DAJP(6),
     1     DAJP(7), DAJP(8), DAJP(9), DAJP(10),DAJP(11),DAJP(12),
     2     DAJP(13),DAJP(14),DAJP(15),DAJP(16),DAJP(17),DAJP(18),
     3     DAJP(19)        / 6.53219131311457D-02,-1.20262933688823D-01,
     4 9.78010236263823D-03, 1.67948429230505D-02,-1.97146140182132D-03,
     5-8.45560295098867D-04, 9.42889620701976D-05, 2.25827860945475D-05,
     6-2.29067870915987D-06,-3.76343991136919D-07, 3.45663933559565D-08,
     7 4.29611332003007D-09,-3.58673691214989D-10,-3.57245881361895D-11,
     8 2.72696091066336D-12, 2.26120653095771D-13,-1.58763205238303D-14,
     9-1.12604374485125D-15, 7.31327529515367D-17/
      DATA DAJN(1), DAJN(2), DAJN(3), DAJN(4), DAJN(5), DAJN(6),
     1     DAJN(7), DAJN(8), DAJN(9), DAJN(10),DAJN(11),DAJN(12),
     2     DAJN(13),DAJN(14),DAJN(15),DAJN(16),DAJN(17),DAJN(18),
     3     DAJN(19)        / 1.08594539632967D-02, 8.53313194857091D-02,
     4-3.15277068113058D-01,-8.78420725294257D-02, 5.53251906976048D-02,
     5 9.41674060503241D-03,-3.32187026018996D-03,-4.11157343156826D-04,
     6 1.01297326891346D-04, 9.87633682208396D-06,-1.87312969812393D-06,
     7-1.50798500131468D-07, 2.32687669525394D-08, 1.59599917419225D-09,
     8-2.07665922668385D-10,-1.24103350500302D-11, 1.39631765331043D-12,
     9 7.39400971155740D-14,-7.32887475627500D-15/
      DATA DA(1),  DA(2),  DA(3),  DA(4),  DA(5),  DA(6),  DA(7),
     1     DA(8),  DA(9),  DA(10), DA(11), DA(12), DA(13), DA(14),
     2     DA(15)          / 4.91627321104601D-01, 3.11164930427489D-03,
     3 8.23140762854081D-05,-4.61769776172142D-06,-6.13158880534626D-08,
     4 2.87295804656520D-08,-1.81959715372117D-09,-1.44752826642035D-10,
     5 4.53724043420422D-11,-3.99655065847223D-12,-3.24089119830323D-13,
     6 1.62098952568741D-13,-2.40765247974057D-14, 1.69384811284491D-16,
     7 8.17900786477396D-16/
      DATA DB(1),  DB(2),  DB(3),  DB(4),  DB(5),  DB(6),  DB(7),
     1     DB(8),  DB(9),  DB(10), DB(11), DB(12), DB(13), DB(14),
     2     DB(15)          /-2.77571356944231D-01, 4.44212833419920D-03,
     3-8.42328522190089D-05,-2.58040318418710D-06, 3.42389720217621D-07,
     4-6.24286894709776D-09,-2.36377836844577D-09, 3.16991042656673D-10,
     5-4.40995691658191D-12,-5.18674221093575D-12, 9.64874015137022D-13,
     6-4.90190576608710D-14,-1.77253430678112D-14, 5.55950610442662D-15,
     7-7.11793337579530D-16/
C***FIRST EXECUTABLE STATEMENT  DJAIRY
      IF (X.LT.0.0D0) GO TO 90
      IF (C.GT.5.0D0) GO TO 60
      IF (X.GT.1.20D0) GO TO 30
      T = (X+X-1.2D0)*CON4
      TT = T + T
      J = N1
      F1 = AK1(J)
      F2 = 0.0D0
      DO 10 I=1,M1
        J = J - 1
        TEMP1 = F1
        F1 = TT*F1 - F2 + AK1(J)
        F2 = TEMP1
   10 CONTINUE
      AI = T*F1 - F2 + AK1(1)
C
      J = N1D
      F1 = DAK1(J)
      F2 = 0.0D0
      DO 20 I=1,M1D
        J = J - 1
        TEMP1 = F1
        F1 = TT*F1 - F2 + DAK1(J)
        F2 = TEMP1
   20 CONTINUE
      DAI = -(T*F1-F2+DAK1(1))
      RETURN
C
   30 CONTINUE
      T = (X+X-CON2)*CON3
      TT = T + T
      J = N2
      F1 = AK2(J)
      F2 = 0.0D0
      DO 40 I=1,M2
        J = J - 1
        TEMP1 = F1
        F1 = TT*F1 - F2 + AK2(J)
        F2 = TEMP1
   40 CONTINUE
      RTRX = SQRT(RX)
      EC = EXP(-C)
      AI = EC*(T*F1-F2+AK2(1))/RTRX
      J = N2D
      F1 = DAK2(J)
      F2 = 0.0D0
      DO 50 I=1,M2D
        J = J - 1
        TEMP1 = F1
        F1 = TT*F1 - F2 + DAK2(J)
        F2 = TEMP1
   50 CONTINUE
      DAI = -EC*(T*F1-F2+DAK2(1))*RTRX
      RETURN
C
   60 CONTINUE
      T = 10.0D0/C - 1.0D0
      TT = T + T
      J = N1
      F1 = AK3(J)
      F2 = 0.0D0
      DO 70 I=1,M1
        J = J - 1
        TEMP1 = F1
        F1 = TT*F1 - F2 + AK3(J)
        F2 = TEMP1
   70 CONTINUE
      RTRX = SQRT(RX)
      EC = EXP(-C)
      AI = EC*(T*F1-F2+AK3(1))/RTRX
      J = N1D
      F1 = DAK3(J)
      F2 = 0.0D0
      DO 80 I=1,M1D
        J = J - 1
        TEMP1 = F1
        F1 = TT*F1 - F2 + DAK3(J)
        F2 = TEMP1
   80 CONTINUE
      DAI = -RTRX*EC*(T*F1-F2+DAK3(1))
      RETURN
C
   90 CONTINUE
      IF (C.GT.5.0D0) GO TO 120
      T = 0.4D0*C - 1.0D0
      TT = T + T
      J = N3
      F1 = AJP(J)
      E1 = AJN(J)
      F2 = 0.0D0
      E2 = 0.0D0
      DO 100 I=1,M3
        J = J - 1
        TEMP1 = F1
        TEMP2 = E1
        F1 = TT*F1 - F2 + AJP(J)
        E1 = TT*E1 - E2 + AJN(J)
        F2 = TEMP1
        E2 = TEMP2
  100 CONTINUE
      AI = (T*E1-E2+AJN(1)) - X*(T*F1-F2+AJP(1))
      J = N3D
      F1 = DAJP(J)
      E1 = DAJN(J)
      F2 = 0.0D0
      E2 = 0.0D0
      DO 110 I=1,M3D
        J = J - 1
        TEMP1 = F1
        TEMP2 = E1
        F1 = TT*F1 - F2 + DAJP(J)
        E1 = TT*E1 - E2 + DAJN(J)
        F2 = TEMP1
        E2 = TEMP2
  110 CONTINUE
      DAI = X*X*(T*F1-F2+DAJP(1)) + (T*E1-E2+DAJN(1))
      RETURN
C
  120 CONTINUE
      T = 10.0D0/C - 1.0D0
      TT = T + T
      J = N4
      F1 = A(J)
      E1 = B(J)
      F2 = 0.0D0
      E2 = 0.0D0
      DO 130 I=1,M4
        J = J - 1
        TEMP1 = F1
        TEMP2 = E1
        F1 = TT*F1 - F2 + A(J)
        E1 = TT*E1 - E2 + B(J)
        F2 = TEMP1
        E2 = TEMP2
  130 CONTINUE
      TEMP1 = T*F1 - F2 + A(1)
      TEMP2 = T*E1 - E2 + B(1)
      RTRX = SQRT(RX)
      CV = C - FPI12
      CCV = COS(CV)
      SCV = SIN(CV)
      AI = (TEMP1*CCV-TEMP2*SCV)/RTRX
      J = N4D
      F1 = DA(J)
      E1 = DB(J)
      F2 = 0.0D0
      E2 = 0.0D0
      DO 140 I=1,M4D
        J = J - 1
        TEMP1 = F1
        TEMP2 = E1
        F1 = TT*F1 - F2 + DA(J)
        E1 = TT*E1 - E2 + DB(J)
        F2 = TEMP1
        E2 = TEMP2
  140 CONTINUE
      TEMP1 = T*F1 - F2 + DA(1)
      TEMP2 = T*E1 - E2 + DB(1)
      E1 = CCV*CON5 + 0.5D0*SCV
      E2 = SCV*CON5 - 0.5D0*CCV
      DAI = (TEMP1*E1-TEMP2*E2)*RTRX
      RETURN
      END

*DECK DASYJY
      SUBROUTINE DASYJY (FUNJY, X, FNU, FLGJY, IN, Y, WK, IFLW)
C***BEGIN PROLOGUE  DASYJY
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DBESJ and DBESY
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (ASYJY-S, DASYJY-D)
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C                 DASYJY computes Bessel functions J and Y
C               for arguments X.GT.0.0 and orders FNU .GE. 35.0
C               on FLGJY = 1 and FLGJY = -1 respectively
C
C                                  INPUT
C
C      FUNJY - External subroutine JAIRY or YAIRY
C          X - Argument, X.GT.0.0D0
C        FNU - Order of the first Bessel function
C      FLGJY - Selection flag
C              FLGJY =  1.0D0 gives the J function
C              FLGJY = -1.0D0 gives the Y function
C         IN - Number of functions desired, IN = 1 or 2
C
C                                  OUTPUT
C
C         Y  - A vector whose first IN components contain the sequence
C       IFLW - A flag indicating underflow or overflow
C                    return variables for BESJ only
C      WK(1) = 1 - (X/FNU)**2 = W**2
C      WK(2) = SQRT(ABS(WK(1)))
C      WK(3) = ABS(WK(2) - ATAN(WK(2)))  or
C              ABS(LN((1 + WK(2))/(X/FNU)) - WK(2))
C            = ABS((2/3)*ZETA**(3/2))
C      WK(4) = FNU*WK(3)
C      WK(5) = (1.5*WK(3)*FNU)**(1/3) = SQRT(ZETA)*FNU**(1/3)
C      WK(6) = SIGN(1.,W**2)*WK(5)**2 = SIGN(1.,W**2)*ZETA*FNU**(2/3)
C      WK(7) = FNU**(1/3)
C
C     Abstract   **** A Double Precision Routine ****
C         DASYJY implements the uniform asymptotic expansion of
C         the J and Y Bessel functions for FNU.GE.35 and real
C         X.GT.0.0D0. The forms are identical except for a change
C         in sign of some of the terms. This change in sign is
C         accomplished by means of the flag FLGJY = 1 or -1. On
C         FLGJY = 1 the Airy functions AI(X) and DAI(X) are
C         supplied by the external function JAIRY, and on
C         FLGJY = -1 the Airy functions BI(X) and DBI(X) are
C         supplied by the external function YAIRY.
C
C***SEE ALSO  DBESJ, DBESY
C***ROUTINES CALLED  D1MACH, I1MACH
C***REVISION HISTORY  (YYMMDD)
C   750101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   891004  Correction computation of ELIM.  (WRB)
C   891009  Removed unreferenced variable.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910408  Updated the AUTHOR section.  (WRB)
C***END PROLOGUE  DASYJY
      INTEGER I, IFLW, IN, J, JN,JR,JU,K, KB,KLAST,KMAX,KP1, KS, KSP1,
     * KSTEMP, L, LR, LRP1, ISETA, ISETB
      INTEGER I1MACH
      DOUBLE PRECISION ABW2, AKM, ALFA, ALFA1, ALFA2, AP, AR, ASUM, AZ,
     * BETA, BETA1, BETA2, BETA3, BR, BSUM, C, CON1, CON2,
     * CON548,CR,CRZ32, DFI,ELIM, DR,FI, FLGJY, FN, FNU,
     * FN2, GAMA, PHI,  RCZ, RDEN, RELB, RFN2,  RTZ, RZDEN,
     * SA, SB, SUMA, SUMB, S1, TA, TAU, TB, TFN, TOL, TOLS, T2, UPOL,
     *  WK, X, XX, Y, Z, Z32
      DOUBLE PRECISION D1MACH
      DIMENSION Y(*), WK(*), C(65)
      DIMENSION ALFA(26,4), BETA(26,5)
      DIMENSION ALFA1(26,2), ALFA2(26,2)
      DIMENSION BETA1(26,2), BETA2(26,2), BETA3(26,1)
      DIMENSION GAMA(26), KMAX(5), AR(8), BR(10), UPOL(10)
      DIMENSION CR(10), DR(10)
      EQUIVALENCE (ALFA(1,1),ALFA1(1,1))
      EQUIVALENCE (ALFA(1,3),ALFA2(1,1))
      EQUIVALENCE (BETA(1,1),BETA1(1,1))
      EQUIVALENCE (BETA(1,3),BETA2(1,1))
      EQUIVALENCE (BETA(1,5),BETA3(1,1))
      SAVE TOLS, CON1, CON2, CON548, AR, BR, C,
     1 ALFA1, ALFA2, BETA1, BETA2, BETA3, GAMA
      DATA TOLS            /-6.90775527898214D+00/
      DATA CON1,CON2,CON548/
     1 6.66666666666667D-01, 3.33333333333333D-01, 1.04166666666667D-01/
      DATA  AR(1),  AR(2),  AR(3),  AR(4),  AR(5),  AR(6),  AR(7),
     A      AR(8)          / 8.35503472222222D-02, 1.28226574556327D-01,
     1 2.91849026464140D-01, 8.81627267443758D-01, 3.32140828186277D+00,
     2 1.49957629868626D+01, 7.89230130115865D+01, 4.74451538868264D+02/
      DATA  BR(1), BR(2), BR(3), BR(4), BR(5), BR(6), BR(7), BR(8),
     A      BR(9), BR(10)  /-1.45833333333333D-01,-9.87413194444444D-02,
     1-1.43312053915895D-01,-3.17227202678414D-01,-9.42429147957120D-01,
     2-3.51120304082635D+00,-1.57272636203680D+01,-8.22814390971859D+01,
     3-4.92355370523671D+02,-3.31621856854797D+03/
      DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10),
     1     C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18),
     2     C(19), C(20), C(21), C(22), C(23), C(24)/
     3       -2.08333333333333D-01,        1.25000000000000D-01,
     4        3.34201388888889D-01,       -4.01041666666667D-01,
     5        7.03125000000000D-02,       -1.02581259645062D+00,
     6        1.84646267361111D+00,       -8.91210937500000D-01,
     7        7.32421875000000D-02,        4.66958442342625D+00,
     8       -1.12070026162230D+01,        8.78912353515625D+00,
     9       -2.36408691406250D+00,        1.12152099609375D-01,
     A       -2.82120725582002D+01,        8.46362176746007D+01,
     B       -9.18182415432400D+01,        4.25349987453885D+01,
     C       -7.36879435947963D+00,        2.27108001708984D-01,
     D        2.12570130039217D+02,       -7.65252468141182D+02,
     E        1.05999045252800D+03,       -6.99579627376133D+02/
      DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32),
     1     C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40),
     2     C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48)/
     3        2.18190511744212D+02,       -2.64914304869516D+01,
     4        5.72501420974731D-01,       -1.91945766231841D+03,
     5        8.06172218173731D+03,       -1.35865500064341D+04,
     6        1.16553933368645D+04,       -5.30564697861340D+03,
     7        1.20090291321635D+03,       -1.08090919788395D+02,
     8        1.72772750258446D+00,        2.02042913309661D+04,
     9       -9.69805983886375D+04,        1.92547001232532D+05,
     A       -2.03400177280416D+05,        1.22200464983017D+05,
     B       -4.11926549688976D+04,        7.10951430248936D+03,
     C       -4.93915304773088D+02,        6.07404200127348D+00,
     D       -2.42919187900551D+05,        1.31176361466298D+06,
     E       -2.99801591853811D+06,        3.76327129765640D+06/
      DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56),
     1     C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64),
     2     C(65)/
     3       -2.81356322658653D+06,        1.26836527332162D+06,
     4       -3.31645172484564D+05,        4.52187689813627D+04,
     5       -2.49983048181121D+03,        2.43805296995561D+01,
     6        3.28446985307204D+06,       -1.97068191184322D+07,
     7        5.09526024926646D+07,       -7.41051482115327D+07,
     8        6.63445122747290D+07,       -3.75671766607634D+07,
     9        1.32887671664218D+07,       -2.78561812808645D+06,
     A        3.08186404612662D+05,       -1.38860897537170D+04,
     B        1.10017140269247D+02/
      DATA ALFA1(1,1), ALFA1(2,1), ALFA1(3,1), ALFA1(4,1), ALFA1(5,1),
     1     ALFA1(6,1), ALFA1(7,1), ALFA1(8,1), ALFA1(9,1), ALFA1(10,1),
     2     ALFA1(11,1),ALFA1(12,1),ALFA1(13,1),ALFA1(14,1),ALFA1(15,1),
     3     ALFA1(16,1),ALFA1(17,1),ALFA1(18,1),ALFA1(19,1),ALFA1(20,1),
     4     ALFA1(21,1),ALFA1(22,1),ALFA1(23,1),ALFA1(24,1),ALFA1(25,1),
     5     ALFA1(26,1)     /-4.44444444444444D-03,-9.22077922077922D-04,
     6-8.84892884892885D-05, 1.65927687832450D-04, 2.46691372741793D-04,
     7 2.65995589346255D-04, 2.61824297061501D-04, 2.48730437344656D-04,
     8 2.32721040083232D-04, 2.16362485712365D-04, 2.00738858762752D-04,
     9 1.86267636637545D-04, 1.73060775917876D-04, 1.61091705929016D-04,
     1 1.50274774160908D-04, 1.40503497391270D-04, 1.31668816545923D-04,
     2 1.23667445598253D-04, 1.16405271474738D-04, 1.09798298372713D-04,
     3 1.03772410422993D-04, 9.82626078369363D-05, 9.32120517249503D-05,
     4 8.85710852478712D-05, 8.42963105715700D-05, 8.03497548407791D-05/
      DATA ALFA1(1,2), ALFA1(2,2), ALFA1(3,2), ALFA1(4,2), ALFA1(5,2),
     1     ALFA1(6,2), ALFA1(7,2), ALFA1(8,2), ALFA1(9,2), ALFA1(10,2),
     2     ALFA1(11,2),ALFA1(12,2),ALFA1(13,2),ALFA1(14,2),ALFA1(15,2),
     3     ALFA1(16,2),ALFA1(17,2),ALFA1(18,2),ALFA1(19,2),ALFA1(20,2),
     4     ALFA1(21,2),ALFA1(22,2),ALFA1(23,2),ALFA1(24,2),ALFA1(25,2),
     5     ALFA1(26,2)     / 6.93735541354589D-04, 2.32241745182922D-04,
     6-1.41986273556691D-05,-1.16444931672049D-04,-1.50803558053049D-04,
     7-1.55121924918096D-04,-1.46809756646466D-04,-1.33815503867491D-04,
     8-1.19744975684254D-04,-1.06184319207974D-04,-9.37699549891194D-05,
     9-8.26923045588193D-05,-7.29374348155221D-05,-6.44042357721016D-05,
     1-5.69611566009369D-05,-5.04731044303562D-05,-4.48134868008883D-05,
     2-3.98688727717599D-05,-3.55400532972042D-05,-3.17414256609022D-05,
     3-2.83996793904175D-05,-2.54522720634871D-05,-2.28459297164725D-05,
     4-2.05352753106481D-05,-1.84816217627666D-05,-1.66519330021394D-05/
      DATA ALFA2(1,1), ALFA2(2,1), ALFA2(3,1), ALFA2(4,1), ALFA2(5,1),
     1     ALFA2(6,1), ALFA2(7,1), ALFA2(8,1), ALFA2(9,1), ALFA2(10,1),
     2     ALFA2(11,1),ALFA2(12,1),ALFA2(13,1),ALFA2(14,1),ALFA2(15,1),
     3     ALFA2(16,1),ALFA2(17,1),ALFA2(18,1),ALFA2(19,1),ALFA2(20,1),
     4     ALFA2(21,1),ALFA2(22,1),ALFA2(23,1),ALFA2(24,1),ALFA2(25,1),
     5     ALFA2(26,1)     /-3.54211971457744D-04,-1.56161263945159D-04,
     6 3.04465503594936D-05, 1.30198655773243D-04, 1.67471106699712D-04,
     7 1.70222587683593D-04, 1.56501427608595D-04, 1.36339170977445D-04,
     8 1.14886692029825D-04, 9.45869093034688D-05, 7.64498419250898D-05,
     9 6.07570334965197D-05, 4.74394299290509D-05, 3.62757512005344D-05,
     1 2.69939714979225D-05, 1.93210938247939D-05, 1.30056674793963D-05,
     2 7.82620866744497D-06, 3.59257485819352D-06, 1.44040049814252D-07,
     3-2.65396769697939D-06,-4.91346867098486D-06,-6.72739296091248D-06,
     4-8.17269379678658D-06,-9.31304715093561D-06,-1.02011418798016D-05/
      DATA ALFA2(1,2), ALFA2(2,2), ALFA2(3,2), ALFA2(4,2), ALFA2(5,2),
     1     ALFA2(6,2), ALFA2(7,2), ALFA2(8,2), ALFA2(9,2), ALFA2(10,2),
     2     ALFA2(11,2),ALFA2(12,2),ALFA2(13,2),ALFA2(14,2),ALFA2(15,2),
     3     ALFA2(16,2),ALFA2(17,2),ALFA2(18,2),ALFA2(19,2),ALFA2(20,2),
     4     ALFA2(21,2),ALFA2(22,2),ALFA2(23,2),ALFA2(24,2),ALFA2(25,2),
     5     ALFA2(26,2)     / 3.78194199201773D-04, 2.02471952761816D-04,
     6-6.37938506318862D-05,-2.38598230603006D-04,-3.10916256027362D-04,
     7-3.13680115247576D-04,-2.78950273791323D-04,-2.28564082619141D-04,
     8-1.75245280340847D-04,-1.25544063060690D-04,-8.22982872820208D-05,
     9-4.62860730588116D-05,-1.72334302366962D-05, 5.60690482304602D-06,
     1 2.31395443148287D-05, 3.62642745856794D-05, 4.58006124490189D-05,
     2 5.24595294959114D-05, 5.68396208545815D-05, 5.94349820393104D-05,
     3 6.06478527578422D-05, 6.08023907788436D-05, 6.01577894539460D-05,
     4 5.89199657344698D-05, 5.72515823777593D-05, 5.52804375585853D-05/
      DATA BETA1(1,1), BETA1(2,1), BETA1(3,1), BETA1(4,1), BETA1(5,1),
     1     BETA1(6,1), BETA1(7,1), BETA1(8,1), BETA1(9,1), BETA1(10,1),
     2     BETA1(11,1),BETA1(12,1),BETA1(13,1),BETA1(14,1),BETA1(15,1),
     3     BETA1(16,1),BETA1(17,1),BETA1(18,1),BETA1(19,1),BETA1(20,1),
     4     BETA1(21,1),BETA1(22,1),BETA1(23,1),BETA1(24,1),BETA1(25,1),
     5     BETA1(26,1)     / 1.79988721413553D-02, 5.59964911064388D-03,
     6 2.88501402231133D-03, 1.80096606761054D-03, 1.24753110589199D-03,
     7 9.22878876572938D-04, 7.14430421727287D-04, 5.71787281789705D-04,
     8 4.69431007606482D-04, 3.93232835462917D-04, 3.34818889318298D-04,
     9 2.88952148495752D-04, 2.52211615549573D-04, 2.22280580798883D-04,
     1 1.97541838033063D-04, 1.76836855019718D-04, 1.59316899661821D-04,
     2 1.44347930197334D-04, 1.31448068119965D-04, 1.20245444949303D-04,
     3 1.10449144504599D-04, 1.01828770740567D-04, 9.41998224204238D-05,
     4 8.74130545753834D-05, 8.13466262162801D-05, 7.59002269646219D-05/
      DATA BETA1(1,2), BETA1(2,2), BETA1(3,2), BETA1(4,2), BETA1(5,2),
     1     BETA1(6,2), BETA1(7,2), BETA1(8,2), BETA1(9,2), BETA1(10,2),
     2     BETA1(11,2),BETA1(12,2),BETA1(13,2),BETA1(14,2),BETA1(15,2),
     3     BETA1(16,2),BETA1(17,2),BETA1(18,2),BETA1(19,2),BETA1(20,2),
     4     BETA1(21,2),BETA1(22,2),BETA1(23,2),BETA1(24,2),BETA1(25,2),
     5     BETA1(26,2)     /-1.49282953213429D-03,-8.78204709546389D-04,
     6-5.02916549572035D-04,-2.94822138512746D-04,-1.75463996970783D-04,
     7-1.04008550460816D-04,-5.96141953046458D-05,-3.12038929076098D-05,
     8-1.26089735980230D-05,-2.42892608575730D-07, 8.05996165414274D-06,
     9 1.36507009262147D-05, 1.73964125472926D-05, 1.98672978842134D-05,
     1 2.14463263790823D-05, 2.23954659232457D-05, 2.28967783814713D-05,
     2 2.30785389811178D-05, 2.30321976080909D-05, 2.28236073720349D-05,
     3 2.25005881105292D-05, 2.20981015361991D-05, 2.16418427448104D-05,
     4 2.11507649256221D-05, 2.06388749782171D-05, 2.01165241997082D-05/
      DATA BETA2(1,1), BETA2(2,1), BETA2(3,1), BETA2(4,1), BETA2(5,1),
     1     BETA2(6,1), BETA2(7,1), BETA2(8,1), BETA2(9,1), BETA2(10,1),
     2     BETA2(11,1),BETA2(12,1),BETA2(13,1),BETA2(14,1),BETA2(15,1),
     3     BETA2(16,1),BETA2(17,1),BETA2(18,1),BETA2(19,1),BETA2(20,1),
     4     BETA2(21,1),BETA2(22,1),BETA2(23,1),BETA2(24,1),BETA2(25,1),
     5     BETA2(26,1)     / 5.52213076721293D-04, 4.47932581552385D-04,
     6 2.79520653992021D-04, 1.52468156198447D-04, 6.93271105657044D-05,
     7 1.76258683069991D-05,-1.35744996343269D-05,-3.17972413350427D-05,
     8-4.18861861696693D-05,-4.69004889379141D-05,-4.87665447413787D-05,
     9-4.87010031186735D-05,-4.74755620890087D-05,-4.55813058138628D-05,
     1-4.33309644511266D-05,-4.09230193157750D-05,-3.84822638603221D-05,
     2-3.60857167535411D-05,-3.37793306123367D-05,-3.15888560772110D-05,
     3-2.95269561750807D-05,-2.75978914828336D-05,-2.58006174666884D-05,
     4-2.41308356761280D-05,-2.25823509518346D-05,-2.11479656768913D-05/
      DATA BETA2(1,2), BETA2(2,2), BETA2(3,2), BETA2(4,2), BETA2(5,2),
     1     BETA2(6,2), BETA2(7,2), BETA2(8,2), BETA2(9,2), BETA2(10,2),
     2     BETA2(11,2),BETA2(12,2),BETA2(13,2),BETA2(14,2),BETA2(15,2),
     3     BETA2(16,2),BETA2(17,2),BETA2(18,2),BETA2(19,2),BETA2(20,2),
     4     BETA2(21,2),BETA2(22,2),BETA2(23,2),BETA2(24,2),BETA2(25,2),
     5     BETA2(26,2)     /-4.74617796559960D-04,-4.77864567147321D-04,
     6-3.20390228067038D-04,-1.61105016119962D-04,-4.25778101285435D-05,
     7 3.44571294294968D-05, 7.97092684075675D-05, 1.03138236708272D-04,
     8 1.12466775262204D-04, 1.13103642108481D-04, 1.08651634848774D-04,
     9 1.01437951597662D-04, 9.29298396593364D-05, 8.40293133016090D-05,
     1 7.52727991349134D-05, 6.69632521975731D-05, 5.92564547323195D-05,
     2 5.22169308826976D-05, 4.58539485165361D-05, 4.01445513891487D-05,
     3 3.50481730031328D-05, 3.05157995034347D-05, 2.64956119950516D-05,
     4 2.29363633690998D-05, 1.97893056664022D-05, 1.70091984636413D-05/
      DATA BETA3(1,1), BETA3(2,1), BETA3(3,1), BETA3(4,1), BETA3(5,1),
     1     BETA3(6,1), BETA3(7,1), BETA3(8,1), BETA3(9,1), BETA3(10,1),
     2     BETA3(11,1),BETA3(12,1),BETA3(13,1),BETA3(14,1),BETA3(15,1),
     3     BETA3(16,1),BETA3(17,1),BETA3(18,1),BETA3(19,1),BETA3(20,1),
     4     BETA3(21,1),BETA3(22,1),BETA3(23,1),BETA3(24,1),BETA3(25,1),
     5     BETA3(26,1)     / 7.36465810572578D-04, 8.72790805146194D-04,
     6 6.22614862573135D-04, 2.85998154194304D-04, 3.84737672879366D-06,
     7-1.87906003636972D-04,-2.97603646594555D-04,-3.45998126832656D-04,
     8-3.53382470916038D-04,-3.35715635775049D-04,-3.04321124789040D-04,
     9-2.66722723047613D-04,-2.27654214122820D-04,-1.89922611854562D-04,
     1-1.55058918599094D-04,-1.23778240761874D-04,-9.62926147717644D-05,
     2-7.25178327714425D-05,-5.22070028895634D-05,-3.50347750511901D-05,
     3-2.06489761035552D-05,-8.70106096849767D-06, 1.13698686675100D-06,
     4 9.16426474122779D-06, 1.56477785428873D-05, 2.08223629482467D-05/
      DATA GAMA(1),   GAMA(2),   GAMA(3),   GAMA(4),   GAMA(5),
     1     GAMA(6),   GAMA(7),   GAMA(8),   GAMA(9),   GAMA(10),
     2     GAMA(11),  GAMA(12),  GAMA(13),  GAMA(14),  GAMA(15),
     3     GAMA(16),  GAMA(17),  GAMA(18),  GAMA(19),  GAMA(20),
     4     GAMA(21),  GAMA(22),  GAMA(23),  GAMA(24),  GAMA(25),
     5     GAMA(26)        / 6.29960524947437D-01, 2.51984209978975D-01,
     6 1.54790300415656D-01, 1.10713062416159D-01, 8.57309395527395D-02,
     7 6.97161316958684D-02, 5.86085671893714D-02, 5.04698873536311D-02,
     8 4.42600580689155D-02, 3.93720661543510D-02, 3.54283195924455D-02,
     9 3.21818857502098D-02, 2.94646240791158D-02, 2.71581677112934D-02,
     1 2.51768272973862D-02, 2.34570755306079D-02, 2.19508390134907D-02,
     2 2.06210828235646D-02, 1.94388240897881D-02, 1.83810633800683D-02,
     3 1.74293213231963D-02, 1.65685837786612D-02, 1.57865285987918D-02,
     4 1.50729501494096D-02, 1.44193250839955D-02, 1.38184805735342D-02/
C***FIRST EXECUTABLE STATEMENT  DASYJY
      TA = D1MACH(3)
      TOL = MAX(TA,1.0D-15)
      TB = D1MACH(5)
      JU = I1MACH(15)
      IF(FLGJY.EQ.1.0D0) GO TO 6
      JR = I1MACH(14)
      ELIM = -2.303D0*TB*(JU+JR)
      GO TO 7
    6 CONTINUE
      ELIM = -2.303D0*(TB*JU+3.0D0)
    7 CONTINUE
      FN = FNU
      IFLW = 0
      DO 170 JN=1,IN
        XX = X/FN
        WK(1) = 1.0D0 - XX*XX
        ABW2 = ABS(WK(1))
        WK(2) = SQRT(ABW2)
        WK(7) = FN**CON2
        IF (ABW2.GT.0.27750D0) GO TO 80
C
C     ASYMPTOTIC EXPANSION
C     CASES NEAR X=FN, ABS(1.-(X/FN)**2).LE.0.2775
C     COEFFICIENTS OF ASYMPTOTIC EXPANSION BY SERIES
C
C     ZETA AND TRUNCATION FOR A(ZETA) AND B(ZETA) SERIES
C
C     KMAX IS TRUNCATION INDEX FOR A(ZETA) AND B(ZETA) SERIES=MAX(2,SA)
C
        SA = 0.0D0
        IF (ABW2.EQ.0.0D0) GO TO 10
        SA = TOLS/LOG(ABW2)
   10   SB = SA
        DO 20 I=1,5
          AKM = MAX(SA,2.0D0)
          KMAX(I) = INT(AKM)
          SA = SA + SB
   20   CONTINUE
        KB = KMAX(5)
        KLAST = KB - 1
        SA = GAMA(KB)
        DO 30 K=1,KLAST
          KB = KB - 1
          SA = SA*WK(1) + GAMA(KB)
   30   CONTINUE
        Z = WK(1)*SA
        AZ = ABS(Z)
        RTZ = SQRT(AZ)
        WK(3) = CON1*AZ*RTZ
        WK(4) = WK(3)*FN
        WK(5) = RTZ*WK(7)
        WK(6) = -WK(5)*WK(5)
        IF(Z.LE.0.0D0) GO TO 35
        IF(WK(4).GT.ELIM) GO TO 75
        WK(6) = -WK(6)
   35   CONTINUE
        PHI = SQRT(SQRT(SA+SA+SA+SA))
C
C     B(ZETA) FOR S=0
C
        KB = KMAX(5)
        KLAST = KB - 1
        SB = BETA(KB,1)
        DO 40 K=1,KLAST
          KB = KB - 1
          SB = SB*WK(1) + BETA(KB,1)
   40   CONTINUE
        KSP1 = 1
        FN2 = FN*FN
        RFN2 = 1.0D0/FN2
        RDEN = 1.0D0
        ASUM = 1.0D0
        RELB = TOL*ABS(SB)
        BSUM = SB
        DO 60 KS=1,4
          KSP1 = KSP1 + 1
          RDEN = RDEN*RFN2
C
C     A(ZETA) AND B(ZETA) FOR S=1,2,3,4
C
          KSTEMP = 5 - KS
          KB = KMAX(KSTEMP)
          KLAST = KB - 1
          SA = ALFA(KB,KS)
          SB = BETA(KB,KSP1)
          DO 50 K=1,KLAST
            KB = KB - 1
            SA = SA*WK(1) + ALFA(KB,KS)
            SB = SB*WK(1) + BETA(KB,KSP1)
   50     CONTINUE
          TA = SA*RDEN
          TB = SB*RDEN
          ASUM = ASUM + TA
          BSUM = BSUM + TB
          IF (ABS(TA).LE.TOL .AND. ABS(TB).LE.RELB) GO TO 70
   60   CONTINUE
   70   CONTINUE
        BSUM = BSUM/(FN*WK(7))
        GO TO 160
C
   75   CONTINUE
        IFLW = 1
        RETURN
C
   80   CONTINUE
        UPOL(1) = 1.0D0
        TAU = 1.0D0/WK(2)
        T2 = 1.0D0/WK(1)
        IF (WK(1).GE.0.0D0) GO TO 90
C
C     CASES FOR (X/FN).GT.SQRT(1.2775)
C
        WK(3) = ABS(WK(2)-ATAN(WK(2)))
        WK(4) = WK(3)*FN
        RCZ = -CON1/WK(4)
        Z32 = 1.5D0*WK(3)
        RTZ = Z32**CON2
        WK(5) = RTZ*WK(7)
        WK(6) = -WK(5)*WK(5)
        GO TO 100
   90   CONTINUE
C
C     CASES FOR (X/FN).LT.SQRT(0.7225)
C
        WK(3) = ABS(LOG((1.0D0+WK(2))/XX)-WK(2))
        WK(4) = WK(3)*FN
        RCZ = CON1/WK(4)
        IF(WK(4).GT.ELIM) GO TO 75
        Z32 = 1.5D0*WK(3)
        RTZ = Z32**CON2
        WK(7) = FN**CON2
        WK(5) = RTZ*WK(7)
        WK(6) = WK(5)*WK(5)
  100   CONTINUE
        PHI = SQRT((RTZ+RTZ)*TAU)
        TB = 1.0D0
        ASUM = 1.0D0
        TFN = TAU/FN
        RDEN=1.0D0/FN
        RFN2=RDEN*RDEN
        RDEN=1.0D0
        UPOL(2) = (C(1)*T2+C(2))*TFN
        CRZ32 = CON548*RCZ
        BSUM = UPOL(2) + CRZ32
        RELB = TOL*ABS(BSUM)
        AP = TFN
        KS = 0
        KP1 = 2
        RZDEN = RCZ
        L = 2
        ISETA=0
        ISETB=0
        DO 140 LR=2,8,2
C
C     COMPUTE TWO U POLYNOMIALS FOR NEXT A(ZETA) AND B(ZETA)
C
          LRP1 = LR + 1
          DO 120 K=LR,LRP1
            KS = KS + 1
            KP1 = KP1 + 1
            L = L + 1
            S1 = C(L)
            DO 110 J=2,KP1
              L = L + 1
              S1 = S1*T2 + C(L)
  110       CONTINUE
            AP = AP*TFN
            UPOL(KP1) = AP*S1
            CR(KS) = BR(KS)*RZDEN
            RZDEN = RZDEN*RCZ
            DR(KS) = AR(KS)*RZDEN
  120     CONTINUE
          SUMA = UPOL(LRP1)
          SUMB = UPOL(LR+2) + UPOL(LRP1)*CRZ32
          JU = LRP1
          DO 130 JR=1,LR
            JU = JU - 1
            SUMA = SUMA + CR(JR)*UPOL(JU)
            SUMB = SUMB + DR(JR)*UPOL(JU)
  130     CONTINUE
          RDEN=RDEN*RFN2
          TB = -TB
          IF (WK(1).GT.0.0D0) TB = ABS(TB)
          IF(RDEN.LT.TOL) GO TO 131
          ASUM = ASUM + SUMA*TB
          BSUM = BSUM + SUMB*TB
          GO TO 140
  131     IF(ISETA.EQ.1) GO TO 132
          IF(ABS(SUMA).LT.TOL) ISETA=1
          ASUM=ASUM+SUMA*TB
  132     IF(ISETB.EQ.1) GO TO 133
          IF(ABS(SUMB).LT.RELB) ISETB=1
          BSUM=BSUM+SUMB*TB
  133     IF(ISETA.EQ.1 .AND. ISETB.EQ.1) GO TO 150
  140   CONTINUE
  150   TB = WK(5)
        IF (WK(1).GT.0.0D0) TB = -TB
        BSUM = BSUM/TB
C
  160   CONTINUE
        CALL FUNJY(WK(6), WK(5), WK(4), FI, DFI)
        TA=1.0D0/TOL
        TB=D1MACH(1)*TA*1.0D+3
        IF(ABS(FI).GT.TB) GO TO 165
        FI=FI*TA
        DFI=DFI*TA
        PHI=PHI*TOL
  165   CONTINUE
        Y(JN) = FLGJY*PHI*(FI*ASUM+DFI*BSUM)/WK(7)
        FN = FN - FLGJY
  170 CONTINUE
      RETURN
      END

*DECK XERMSG
      SUBROUTINE XERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL)
C***BEGIN PROLOGUE  XERMSG
C***PURPOSE  Process error messages for SLATEC and other libraries.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERMSG-A)
C***KEYWORDS  ERROR MESSAGE, XERROR
C***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
C***DESCRIPTION
C
C   XERMSG processes a diagnostic message in a manner determined by the
C   value of LEVEL and the current value of the library error control
C   flag, KONTRL.  See subroutine XSETF for details.
C
C    LIBRAR   A character constant (or character variable) with the name
C             of the library.  This will be 'SLATEC' for the SLATEC
C             Common Math Library.  The error handling package is
C             general enough to be used by many libraries
C             simultaneously, so it is desirable for the routine that
C             detects and reports an error to identify the library name
C             as well as the routine name.
C
C    SUBROU   A character constant (or character variable) with the name
C             of the routine that detected the error.  Usually it is the
C             name of the routine that is calling XERMSG.  There are
C             some instances where a user callable library routine calls
C             lower level subsidiary routines where the error is
C             detected.  In such cases it may be more informative to
C             supply the name of the routine the user called rather than
C             the name of the subsidiary routine that detected the
C             error.
C
C    MESSG    A character constant (or character variable) with the text
C             of the error or warning message.  In the example below,
C             the message is a character constant that contains a
C             generic message.
C
C                   CALL XERMSG ('SLATEC', 'MMPY',
C                  *'THE ORDER OF THE MATRIX EXCEEDS THE ROW DIMENSION',
C                  *3, 1)
C
C             It is possible (and is sometimes desirable) to generate a
C             specific message--e.g., one that contains actual numeric
C             values.  Specific numeric values can be converted into
C             character strings using formatted WRITE statements into
C             character variables.  This is called standard Fortran
C             internal file I/O and is exemplified in the first three
C             lines of the following example.  You can also catenate
C             substrings of characters to construct the error message.
C             Here is an example showing the use of both writing to
C             an internal file and catenating character strings.
C
C                   CHARACTER*5 CHARN, CHARL
C                   WRITE (CHARN,10) N
C                   WRITE (CHARL,10) LDA
C                10 FORMAT(I5)
C                   CALL XERMSG ('SLATEC', 'MMPY', 'THE ORDER'//CHARN//
C                  *   ' OF THE MATRIX EXCEEDS ITS ROW DIMENSION OF'//
C                  *   CHARL, 3, 1)
C
C             There are two subtleties worth mentioning.  One is that
C             the // for character catenation is used to construct the
C             error message so that no single character constant is
C             continued to the next line.  This avoids confusion as to
C             whether there are trailing blanks at the end of the line.
C             The second is that by catenating the parts of the message
C             as an actual argument rather than encoding the entire
C             message into one large character variable, we avoid
C             having to know how long the message will be in order to
C             declare an adequate length for that large character
C             variable.  XERMSG calls XERPRN to print the message using
C             multiple lines if necessary.  If the message is very long,
C             XERPRN will break it into pieces of 72 characters (as
C             requested by XERMSG) for printing on multiple lines.
C             Also, XERMSG asks XERPRN to prefix each line with ' *  '
C             so that the total line length could be 76 characters.
C             Note also that XERPRN scans the error message backwards
C             to ignore trailing blanks.  Another feature is that
C             the substring '$$' is treated as a new line sentinel
C             by XERPRN.  If you want to construct a multiline
C             message without having to count out multiples of 72
C             characters, just use '$$' as a separator.  '$$'
C             obviously must occur within 72 characters of the
C             start of each line to have its intended effect since
C             XERPRN is asked to wrap around at 72 characters in
C             addition to looking for '$$'.
C
C    NERR     An integer value that is chosen by the library routine's
C             author.  It must be in the range -99 to 999 (three
C             printable digits).  Each distinct error should have its
C             own error number.  These error numbers should be described
C             in the machine readable documentation for the routine.
C             The error numbers need be unique only within each routine,
C             so it is reasonable for each routine to start enumerating
C             errors from 1 and proceeding to the next integer.
C
C    LEVEL    An integer value in the range 0 to 2 that indicates the
C             level (severity) of the error.  Their meanings are
C
C            -1  A warning message.  This is used if it is not clear
C                that there really is an error, but the user's attention
C                may be needed.  An attempt is made to only print this
C                message once.
C
C             0  A warning message.  This is used if it is not clear
C                that there really is an error, but the user's attention
C                may be needed.
C
C             1  A recoverable error.  This is used even if the error is
C                so serious that the routine cannot return any useful
C                answer.  If the user has told the error package to
C                return after recoverable errors, then XERMSG will
C                return to the Library routine which can then return to
C                the user's routine.  The user may also permit the error
C                package to terminate the program upon encountering a
C                recoverable error.
C
C             2  A fatal error.  XERMSG will not return to its caller
C                after it receives a fatal error.  This level should
C                hardly ever be used; it is much better to allow the
C                user a chance to recover.  An example of one of the few
C                cases in which it is permissible to declare a level 2
C                error is a reverse communication Library routine that
C                is likely to be called repeatedly until it integrates
C                across some interval.  If there is a serious error in
C                the input such that another step cannot be taken and
C                the Library routine is called again without the input
C                error having been corrected by the caller, the Library
C                routine will probably be called forever with improper
C                input.  In this case, it is reasonable to declare the
C                error to be fatal.
C
C    Each of the arguments to XERMSG is input; none will be modified by
C    XERMSG.  A routine may make multiple calls to XERMSG with warning
C    level messages; however, after a call to XERMSG with a recoverable
C    error, the routine should return to the user.  Do not try to call
C    XERMSG with a second recoverable error after the first recoverable
C    error because the error package saves the error number.  The user
C    can retrieve this error number by calling another entry point in
C    the error handling package and then clear the error number when
C    recovering from the error.  Calling XERMSG in succession causes the
C    old error number to be overwritten by the latest error number.
C    This is considered harmless for error numbers associated with
C    warning messages but must not be done for error numbers of serious
C    errors.  After a call to XERMSG with a recoverable error, the user
C    must be given a chance to call NUMXER or XERCLR to retrieve or
C    clear the error number.
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  FDUMP, J4SAVE, XERCNT, XERHLT, XERPRN, XERSVE
C***REVISION HISTORY  (YYMMDD)
C   880101  DATE WRITTEN
C   880621  REVISED AS DIRECTED AT SLATEC CML MEETING OF FEBRUARY 1988.
C           THERE ARE TWO BASIC CHANGES.
C           1.  A NEW ROUTINE, XERPRN, IS USED INSTEAD OF XERPRT TO
C               PRINT MESSAGES.  THIS ROUTINE WILL BREAK LONG MESSAGES
C               INTO PIECES FOR PRINTING ON MULTIPLE LINES.  '$$' IS
C               ACCEPTED AS A NEW LINE SENTINEL.  A PREFIX CAN BE
C               ADDED TO EACH LINE TO BE PRINTED.  XERMSG USES EITHER
C               ' ***' OR ' *  ' AND LONG MESSAGES ARE BROKEN EVERY
C               72 CHARACTERS (AT MOST) SO THAT THE MAXIMUM LINE
C               LENGTH OUTPUT CAN NOW BE AS GREAT AS 76.
C           2.  THE TEXT OF ALL MESSAGES IS NOW IN UPPER CASE SINCE THE
C               FORTRAN STANDARD DOCUMENT DOES NOT ADMIT THE EXISTENCE
C               OF LOWER CASE.
C   880708  REVISED AFTER THE SLATEC CML MEETING OF JUNE 29 AND 30.
C           THE PRINCIPAL CHANGES ARE
C           1.  CLARIFY COMMENTS IN THE PROLOGUES
C           2.  RENAME XRPRNT TO XERPRN
C           3.  REWORK HANDLING OF '$$' IN XERPRN TO HANDLE BLANK LINES
C               SIMILAR TO THE WAY FORMAT STATEMENTS HANDLE THE /
C               CHARACTER FOR NEW RECORDS.
C   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
C           CLEAN UP THE CODING.
C   890721  REVISED TO USE NEW FEATURE IN XERPRN TO COUNT CHARACTERS IN
C           PREFIX.
C   891013  REVISED TO CORRECT COMMENTS.
C   891214  Prologue converted to Version 4.0 format.  (WRB)
C   900510  Changed test on NERR to be -9999999 < NERR < 99999999, but
C           NERR .ne. 0, and on LEVEL to be -2 < LEVEL < 3.  Added
C           LEVEL=-1 logic, changed calls to XERSAV to XERSVE, and
C           XERCTL to XERCNT.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERMSG
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      CHARACTER*8 XLIBR, XSUBR
      CHARACTER*72  TEMP
      CHARACTER*20  LFIRST
      INTEGER NERR,LEVEL
      RETURN
      END
*

      SUBROUTINE coulfg(xx,eta1,xlmin,xlmax,fc,gc,fcp,gcp,mode1,kfn,
     1                  ifail)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c  revised coulomb wavefunction program using steed's method           c
c                                                                      c
c  a. r. barnett           manchester  march   1981                    c
c                                                                      c
c  original program 'rcwfn'      in    cpc  8 (1974) 377-395           c
c                 + 'rcwff'      in    cpc 11 (1976) 141-142           c
c  full description of algorithm in    cpc 21 (1981) 297-314           c
c  this version written up       in    cpc xx (1982) yyy-zzz           c
c                                                                      c
c  coulfg returns f,g,f',g', for real xx.gt.0,real eta1 (including 0), c
c   and real lambda(xlmin) .gt. -1 for integer-spaced lambda values    c
c   thus giving positive-energy solutions to the coulomb schrodinger   c
c   equation,to the klein-gordon equation and to suitable forms of     c
c   the dirac equation ,also spherical + cylindrical bessel equations  c
c                                                                      c
c  for a range of lambda values (xlmax - xlmin) must be an integer,    c
c  starting array element is m1 = max0(  int(xlmin+accur),0) + 1       c
c      see text for modifications for integer l-values                 c
c                                                                      c
c  if 'mode' = 1  get f,g,f',g'   for integer-spaced lambda values     c
c            = 2      f,g      unused arrays must be dimensioned in    c
c            = 3      f               call to at least length (1)      c
c  if 'kfn'  = 0 real        coulomb functions are returned            c
c            = 1 spherical   bessel                                    c
c            = 2 cylindrical bessel                                    c
c  the use of 'mode' and 'kfn' is independent                          c
c                                                                      c
c  precision&  results to within 2-3 decimals of 'machine accuracy'    c
c   in oscillating region x .ge. eta1 + sqrt(eta1**2 + xlm(xlm+1))     c
c   coulfg is coded for real*8 on ibm or equivalent  accur = 10**-16   c
c   use autodbl + extended precision on hx compiler  accur = 10**-33   c
c   for mantissas of 56 + 112 bits. for single precision cdc (48 bits) c
c   reassign dsqrt=sqrt etc.  see text for complex arithmetic version  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     .. scalar arguments ..
      DOUBLE PRECISION eta1,xlmax,xlmin,xx
      INTEGER ifail,kfn,mode1
c     ..
c     .. array arguments ..
      DOUBLE PRECISION fc(1),fcp(1),gc(1),gcp(1)
c     ..
c     .. scalars in common ..
      DOUBLE PRECISION paccq
      INTEGER iexp,m1,nfp,npq
c     ..
c     .. local scalars ..
      DOUBLE PRECISION a,abort,acc,acc4,acch,accur,ai,alpha,ar,b,beta,
     1                 bi,br,c,d,dell,df,di,dp,dq,dr,e2mm1,ek,el,eta,f,
     2                 fcl,fcl1,fcm,fjwkb,fpl,gam,gcl,gcl1,gjwkb,gpl,
     3                 half,one,p,pk,pk1,px,q,rl,rt2epi,sl,ta,ten2,tk,
     4                 tm30,two,w,wi,x,xi,xl,xll,xlm,zero
      INTEGER l,l1,lp,lxtra,maxl,mode
      LOGICAL etane0,xlturn
c     ..
c     .. external subroutines ..
      EXTERNAL jwkb
c     ..
c     .. intrinsic functions ..
      INTRINSIC abs,dble,int,max,max0,min,mod,sign,sqrt
c     ..
c     .. common blocks ..
      COMMON /steed/paccq,nfp,npq,iexp,m1
c     ..
c     .. data statements ..
c***  common block is for information only.  not required in code
c***  coulfg has calls to&  sqrt, abs, amod ,  int, sign, float, min1
      DATA zero,one,two,ten2,abort/0.0D0,1.0D0,2.0D0,1.0D2,2.0D4/
      DATA half,tm30/0.5D0,1.0D-30/
      DATA rt2epi/0.797884560802865D0/
c     ..
      accur = 2.d-15
c ***            change accur to suit machine and precision required
      mode = 1
      IF (mode1.EQ.2 .OR. mode1.EQ.3) mode = mode1
      ifail = 0
      iexp = 1
      npq = 0
      eta = eta1
      gjwkb = zero
      paccq = one
      IF (kfn.NE.0) eta = zero
      etane0 = eta .NE. zero
      acc = accur
      acc4 = acc*ten2*ten2
      acch = sqrt(acc)
c ***    test range of xx, exit if.le. sqrt(accur) or if negative
c
      IF (xx.LE.acch) GO TO 130
      x = xx
      xlm = xlmin
      IF (kfn.EQ.2) xlm = xlm - half
      IF (xlm.LE.-one .OR. xlmax.LT.xlmin) GO TO 140
      e2mm1 = eta*eta + xlm*xlm + xlm
      xlturn = x* (x-two*eta) .LT. xlm*xlm + xlm
      dell = xlmax - xlmin + acc
      IF (abs(mod(dell,one)).GT.acc) WRITE (6,FMT=9060) xlmax,xlmin,
     1    dell
      lxtra = int(dell)
      xll = xlm + dble(lxtra)
c ***       lxtra is number of additional lambda values to be computed
c ***       xll  is max lambda value, or 0.5 smaller for j,y bessels
c ***         determine starting array element (m1) from xlmin
      m1 = max0(int(xlmin+acc),0) + 1
      l1 = m1 + lxtra
c
c ***    evaluate cf1  =  f   =  fprime(xl,eta,x)/f(xl,eta,x)
c
      xi = one/x
      fcl = one
      pk = xll + one
      px = pk + abort
   10 ek = eta/pk
      f = (ek+pk*xi)*fcl + (fcl-one)*xi
      pk1 = pk + one
c ***   test ensures b1 .ne. zero for negative eta; fixup is exact.
      IF (abs(eta*x+pk*pk1).GT.acc) GO TO 20
      fcl = (one+ek*ek)/ (one+ (eta/pk1)**2)
      pk = two + pk
      GO TO 10
C
   20 d = one/ ((pk+pk1)* (xi+ek/pk1))
      df = -fcl* (one+ek*ek)*d
      IF (fcl.NE.one) fcl = -one
      IF (d.LT.zero) fcl = -fcl
      f = f + df
c
c ***   begin cf1 loop on pk = k = lambda + 1
c
      p = one
   30 pk = pk1
      pk1 = pk1 + one
      ek = eta/pk
      tk = (pk+pk1)* (xi+ek/pk1)
      d = tk - d* (one+ek*ek)
      IF (abs(d).GT.acch) GO TO 40
      WRITE (6,FMT=9000) d,df,acch,pk,ek,eta,x
      p = p + one
      IF (p.GT.two) GO TO 150
   40 d = one/d
      IF (d.LT.zero) fcl = -fcl
      df = df* (d*tk-one)
      f = f + df
      IF (pk.GT.px) GO TO 150
      IF (abs(df).GE.abs(f)*acc) GO TO 30
      nfp = pk - xll - 1
      IF (lxtra.EQ.0) GO TO 60
c
c *** downward recurrence to lambda = xlm. array gc,if present,stores rl
c
      fcl = fcl*tm30
      fpl = fcl*f
      IF (mode.EQ.1) fcp(l1) = fpl
      fc(l1) = fcl
      xl = xll
      rl = one
      el = zero
      DO 50 lp = 1,lxtra
          IF (etane0) el = eta/xl
          IF (etane0) rl = sqrt(one+el*el)
          sl = el + xl*xi
          l = l1 - lp
          fcl1 = (fcl*sl+fpl)/rl
          fpl = fcl1*sl - fcl*rl
          fcl = fcl1
          fc(l) = fcl
          IF (mode.EQ.1) fcp(l) = fpl
          IF (mode.NE.3 .AND. etane0) gc(l+1) = rl
          xl = xl - one
   50 CONTINUE
      IF (fcl.EQ.zero) fcl = acc
      f = fpl/fcl
c ***    now we have reached lambda = xlmin = xlm
c ***    evaluate cf2 = p + i.q  again using steed's algorithm
c ***    see text for compact complex code for sp cdc or non-ansi ibm
c
   60 IF (xlturn) CALL jwkb(x,eta,max(xlm,zero),fjwkb,gjwkb,iexp)
      IF (iexp.GT.1 .OR. gjwkb.GT.one/ (acch*ten2)) GO TO 80
      xlturn = .FALSE.
      ta = two*abort
      pk = zero
      wi = eta + eta
      p = zero
      q = one - eta*xi
      ar = -e2mm1
      ai = eta
      br = two* (x-eta)
      bi = two
      dr = br/ (br*br+bi*bi)
      di = -bi/ (br*br+bi*bi)
      dp = -xi* (ar*di+ai*dr)
      dq = xi* (ar*dr-ai*di)
   70 p = p + dp
      q = q + dq
      pk = pk + two
      ar = ar + pk
      ai = ai + wi
      bi = bi + two
      d = ar*dr - ai*di + br
      di = ai*dr + ar*di + bi
      c = one/ (d*d+di*di)
      dr = c*d
      di = -c*di
      a = br*dr - bi*di - one
      b = bi*dr + br*di
      c = dp*a - dq*b
      dq = dp*b + dq*a
      dp = c
      IF (pk.GT.ta) GO TO 160
      IF (abs(dp)+abs(dq).GE. (abs(p)+abs(q))*acc) GO TO 70
      npq = pk/two
      paccq = half*acc/min(abs(q),one)
      IF (abs(p).GT.abs(q)) paccq = paccq*abs(p)
c
c *** solve for fcm = f at lambda = xlm,then find norm factor w=w/fcm
c
      gam = (f-p)/q
      IF (q.LE.acc4*abs(p)) GO TO 170
      w = one/sqrt((f-p)*gam+q)
      GO TO 90
c *** arrive here if g(xlm) .gt. 10**6 or iexp .gt. 250 + xlturn = .true.
   80 w = fjwkb
      gam = gjwkb*w
      p = f
      q = one
c
c *** normalise for spherical or cylindrical bessel functions
c
   90 alpha = zero
      IF (kfn.EQ.1) alpha = xi
      IF (kfn.EQ.2) alpha = xi*half
      beta = one
      IF (kfn.EQ.1) beta = xi
      IF (kfn.EQ.2) beta = sqrt(xi)*rt2epi
      fcm = sign(w,fcl)*beta
      fc(m1) = fcm
      IF (mode.EQ.3) GO TO 100
      IF (.NOT.xlturn) gcl = fcm*gam
      IF (xlturn) gcl = gjwkb*beta
      IF (kfn.NE.0) gcl = -gcl
      gc(m1) = gcl
      gpl = gcl* (p-q/gam) - alpha*gcl
      IF (mode.EQ.2) GO TO 100
      gcp(m1) = gpl
      fcp(m1) = fcm* (f-alpha)
  100 IF (lxtra.EQ.0) RETURN
c *** upward recurrence from gc(m1),gcp(m1)  stored value is rl
c *** renormalise fc,fcp at each lambda and correct regular derivative
c ***    xl   = xlm here  and rl = one , el = zero for bessels
      w = beta*w/abs(fcl)
      maxl = l1 - 1
      DO 120 l = m1,maxl
          IF (mode.EQ.3) GO TO 110
          xl = xl + one
          IF (etane0) el = eta/xl
          IF (etane0) rl = gc(l+1)
          sl = el + xl*xi
          gcl1 = ((sl-alpha)*gcl-gpl)/rl
          gpl = rl*gcl - (sl+alpha)*gcl1
          gcl = gcl1
          gc(l+1) = gcl1
          IF (mode.EQ.2) GO TO 110
          gcp(l+1) = gpl
          fcp(l+1) = w* (fcp(l+1)-alpha*fc(l+1))
  110     fc(l+1) = w*fc(l+1)
  120 CONTINUE
      RETURN
c
c ***    error messages
c
  130 ifail = -1
      WRITE (6,FMT=9010) xx,acch
      RETURN
C
  140 ifail = -2
      WRITE (6,FMT=9020) xlmax,xlmin,xlm
      RETURN
C
  150 ifail = 1
      WRITE (6,FMT=9030) abort,f,df,pk,px,acc
      RETURN
C
  160 ifail = 2
      WRITE (6,FMT=9040) abort,p,q,dp,dq,acc
      RETURN
C
  170 ifail = 3
      WRITE (6,FMT=9050) p,q,acc,dell,lxtra,m1
      RETURN
C
 9000 FORMAT (/' cf1 accuracy loss& d,df,acch,k,eta/k,eta,x = ',1p,
     1       7d9.2,/)
 9010 FORMAT (' for xx = ',1p,d12.3,' try small-x  solutions',' or x n',
     1       'egative',/' square root accuracy parameter =  ',d12.3,/)
 9020 FORMAT (/' problem with input order values&xlmax,xlmin,xlm = ',1p,
     1       3d15.6,/)
 9030 FORMAT (' cf1 has failed to converge after ',f10.0,' iterations',
     1       /' f,df,pk,px,accur =  ',1p,5d12.3,//)
 9040 FORMAT (' cf2 has failed to converge after ',f7.0,' iterations',
     1       /' p,q,dp,dq,accur =  ',1p,4d17.7,d12.3,//)
 9050 FORMAT (' final q.le. abs(p)*acc*10**4 , p,q,acc = ',1p,3d12.3,4x,
     1       ' dell,lxtra,m1 = ',d12.3,2i5,/)
 9060 FORMAT (' xlmax - xlmin = dell not an integer ',1p,3d20.10,/)
      END
      SUBROUTINE jwkb(xx,eta1,xl,fjwkb,gjwkb,iexp)
c     .. scalar arguments ..
      DOUBLE PRECISION eta1,fjwkb,gjwkb,xl,xx
      INTEGER iexp
c     ..
c     .. local scalars ..
      DOUBLE PRECISION aloge,dzero,eta,gh,gh2,half,hl,hll,one,phi,phi10,
     1                 rl2,rl35,six,sl,ten,x,xll1,zero
c     ..
c     .. intrinsic functions ..
      INTRINSIC atan2,dble,exp,int,log,max,sqrt
c     ..
c     .. data statements ..
c *** computes jwkb approximations to coulomb functions    for xl.ge. 0
c *** as modified by biedenharn et al. phys rev 97 (1955) 542-554
c *** calls amax1,sqrt,alog,exp,atan2,float,int        barnett feb 1981
      DATA zero,half,one,six,ten/0.0D0,0.5D0,1.0D0,6.0D0,10.0D0/
      DATA dzero,rl35,aloge/0.0D0,35.0D0,0.4342945D0/
c     ..
      x = xx
      eta = eta1
      gh2 = x* (eta+eta-x)
      xll1 = max(xl*xl+xl,dzero)
      IF (gh2+xll1.LE.zero) RETURN
      hll = xll1 + six/rl35
      hl = sqrt(hll)
      sl = eta/hl + hl/x
      rl2 = one + eta*eta/hll
      gh = sqrt(gh2+hll)/x
      phi = x*gh - half* (hl*log((gh+sl)**2/rl2)-log(gh))
      IF (eta.NE.zero) phi = phi - eta*atan2(x*gh,x-eta)
      phi10 = -phi*aloge
      iexp = int(phi10)
      IF (iexp.GT.250) gjwkb = ten** (phi10-dble(iexp))
      IF (iexp.LE.250) gjwkb = exp(-phi)
      IF (iexp.LE.250) iexp = 0
      fjwkb = half/ (gh*gjwkb)
      RETURN
C
      END

