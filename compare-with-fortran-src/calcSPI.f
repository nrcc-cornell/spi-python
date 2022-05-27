C     CALCULATE SPI
C
C     Using the routines from Colorado Climate Center
C
      PARAMETER (MAXYRS=100, N=1200, NRUN=1)

      REAL GAMSHAPE(12), GAMSCALE(12), PZERO(12), PROBNE(MAXYRS*12),
     .     PCPACC(MAXYRS*12), INDEXV(MAXYRS*12), PRECIP(N)

C      OPEN(UNIT=10,FILE='precip_data_month_045115.txt',STATUS='old')
      OPEN(UNIT=10,FILE='precip_data_month.txt',STATUS='old')

C     Read monthly precip data from file
      DO 100 I=1,N
          READ(10,*) PRECIP(I)
 100  CONTINUE

      CALL spigam(NRUN,PRECIP,GAMSCALE,GAMSHAPE,PZERO,INDEXV,PROBNE,
     .     PCPACC)

C      write(*,*) (GAMSHAPE(I), DO I=1,12)
C      write(*,*) (GAMSCALE(I), DO I=1,12)

      NMON = MAXYRS*12
      DO 150 I=1,NMON
          im = mod (I-1,12) + 1
          IF (im .ne. 8) THEN
              GOTO 150
          ENDIF
          WRITE(*,*) GAMSHAPE(im),GAMSCALE(im),PZERO(im),PCPACC(I),
     .               PROBNE(I), INDEXV(I)
 150  CONTINUE

      STOP
      END
