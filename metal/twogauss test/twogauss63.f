      PROGRAM TWOGAUSS5
*-
* after GENFIT by P.Francis, modified by Eileen
*
*  Reads in a list of spectra, sends them to their rest frame,
* rebins part of them, normalised by its mean. This then has
* the mean spectrum subtracted from it.
*
* A 2-G model is then fit to the spectra, and the parameters output
*
* Input:- List of spectra and redshifts.
*
* Calls:- RRFIG, AMOEBA
* A file listing the fit parameters
*
* Modified for WiggleZ, Michael Drinkwater, 2015 January
*     v3 specifically fits [OII] and [OIII] doublets separately, plus 2-compt H-beta
*     v4 -this version fixes NLINES=2-4, NBIN=1
*        ALSO fits linked 2-components to OIII doublet using NLINES=4
*     v5 calculates Bayesian Information Criterion for fits
*     v6 revises fitting parameters, fixing constant continuum among other changes
*     v6.1 adds text output of selected fits
*     v6.3 adds uncertainties
*-

      IMPLICIT NONE

      INTEGER I, J, COLMAX, RMAX, ISPEC, IOERR, IOERR2
      INTEGER SCREEN, KEYBD, INFILE, I4363
      INTEGER OUTFIL1, OUTFIL2, OUTFIL3, OUTFIL0, OUTFIL4, OUTFIL5
      INTEGER SPFIL1, SPFIL2, SPFIL3
      INTEGER COUNT, MP, NP, NLINES,OFLAG
      INTEGER MAXSIZ, NBIN, NUM, NPAR, NPARD

      PARAMETER (COLMAX=800, RMAX=600, SCREEN=6,
     :     KEYBD=5, INFILE=7, OUTFIL1=8, OUTFIL2=10, OUTFIL3=12,  
     :     OUTFIL0=20, OUTFIL4=22, OUTFIL5=24, SPFIL1=14, SPFIL2=16, 
     :     SPFIL3=18, MAXSIZ=8000, MP=40, NP=40)

      REAL WAV(MAXSIZ), SIG(MAXSIZ), ERR(MAXSIZ)
      REAL BINSIZ, RBDAT(MAXSIZ),FITH3,FITH4,CENTRE
      REAL RBWAV(MAXSIZ), RBSIG(MAXSIZ), TRED,SFACT
      REAL FIT1(MAXSIZ), FIT2(MAXSIZ), FIT21(MAXSIZ), FIT22(MAXSIZ)
      REAL INT1,INT1E,FW1,CEN1,INT2,FW2,CEN2,CONTIN
      REAL CHI2,BIC,CORR2,ADCORR2,WAVMIN,WAVMAX
      REAL INT1D,INT1DE,FW1D,CEN1D,INT2D,FW2D,CEN2D
      REAL CHI2D,BICD,CORR2D,ADCORR2D
      REAL OIIF, OIIE, OIIIF, OIIIE, HBF, HBE, HGF, HGAE, OWF, OWE
      REAL EW1, EWB
      CHARACTER*40 ONAME, TNAME, FILNAM, SPECNAME
      CHARACTER*100 CHLINE
      COMMON /SPECTRA/RBWAV,RBDAT,RBSIG,NBIN,NLINES,CENTRE,OFLAG


c     Input/output files
      IOERR=1
      DO WHILE (IOERR.NE.0)
         WRITE (SCREEN,'(A$)')' List of spectra: '
         FILNAM='toto1'
         READ (KEYBD,'(A)') FILNAM
         OPEN (INFILE, FILE=FILNAM, STATUS='OLD',IOSTAT=IOERR)
         IF (IOERR.EQ.0) READ(INFILE,'(A)',IOSTAT=IOERR) CHLINE
      END DO


      WRITE(SCREEN,'(A$)')' Base output file name ('//
     +     FILNAM(1:LNBLNK(FILNAM))//'): '
      READ(KEYBD,'(A)')  ONAME
      IF (ONAME.EQ.'') ONAME=FILNAM

c     Summary text output files
      FILNAM=ONAME(1:LNBLNK(ONAME))//'_4363.txt'
      OPEN(OUTFIL0, FILE=FILNAM)
      WRITE(OUTFIL0,'(A)')
     :     '# OII          NB     Flux  '//
     :     '  FluxE     FWHM     Wave    Sfact     NP    Chi2 '//
     :     '      BIC       R2      AR2    Flux1    FWHM1    Wave1'

      FILNAM=ONAME(1:LNBLNK(ONAME))//'_O2.txt'
      OPEN(OUTFIL1, FILE=FILNAM)
      WRITE(OUTFIL1,'(A)')
     :     '# OII          NB     Flux  '//
     :     '  FluxE     FWHM     Wave    Sfact     NP    Chi2 '//
     :     '      BIC       R2      AR2    Flux1    FWHM1    Wave1'//
     :     '    Flux2    FWHM2    Wave2    Ratio    NP2'//
     :     '     Chi2      BIC       R2      AR2    dBIC      dAR2'
      FILNAM=ONAME(1:LNBLNK(ONAME))//'_O3.txt'
      OPEN(OUTFIL2, FILE=FILNAM)
      WRITE(OUTFIL2,'(A)')
     :     '# OIII         NB     Flux  '//
     :     '  FluxE     FWHM     Wave     Sfact     NP    Chi2 '//
     :     '      BIC       R2      AR2    Flux1    FWHM1    Wave1'//
     :     '    Flux2    FWHM2    Wave2    Ratio    NP2'//
     :     '     Chi2      BIC       R2      AR2    dBIC      dAR2'
      FILNAM=ONAME(1:LNBLNK(ONAME))//'_Hb.txt'
      OPEN(OUTFIL3, FILE=FILNAM)
      WRITE(OUTFIL3,'(A)')
     :     '# Hbeta        NB     Flux  '//
     :     '  FluxE     FWHM     Wave     Sfact     NP    Chi2 '//
     :     '      BIC       R2      AR2    Flux1    FWHM1    Wave1'//
     :     '    Flux2    FWHM2    Wave2    Ratio    NP2'//
     :     '     Chi2      BIC       R2      AR2    dBIC      dAR2'
     
      FILNAM=ONAME(1:LNBLNK(ONAME))//'_Hg.txt'
      OPEN(OUTFIL4, FILE=FILNAM)
      WRITE(OUTFIL4,'(A)')
     :     '# Hgamma        NB     Flux  '//
     :     '  FluxE     FWHM     Wave     Sfact     NP    Chi2 '//
     :     '      BIC       R2      AR2    Flux1    FWHM1    Wave1'//
     :     '    Flux2    FWHM2    Wave2    Ratio    NP2'//
     :     '     Chi2      BIC       R2      AR2    dBIC      dAR2'
     
      FILNAM=ONAME(1:LNBLNK(ONAME))//'_summary.txt'
      OPEN(OUTFIL5, FILE=FILNAM)
      WRITE(OUTFIL5,'(A)')
     :     '# Spectrum     4363 flux       4363 error'//
     :     '    3727 flux     3727 error    4959 flux    4959 error '//
     :     '      Hb flux     hb error   EWhb   hg flux    hg error'

c     spectrum (data + components + fit) output files
      FILNAM=ONAME(1:LNBLNK(ONAME))//'_O2.dat'
      OPEN(SPFIL1, FILE=FILNAM)
      WRITE(SPFIL1,'(A)')
     :     '# Wave Data Sigma Fit1 Fit2a Fit2b Fit2'
      FILNAM=ONAME(1:LNBLNK(ONAME))//'_O3.dat'
      OPEN(SPFIL2, FILE=FILNAM)
      WRITE(SPFIL2,'(A)')
     :     '# Wave Data Sigma Fit1 Fit2a Fit2b Fit2'
      FILNAM=ONAME(1:LNBLNK(ONAME))//'_Hb.dat'
      OPEN(SPFIL3, FILE=FILNAM)
      WRITE(SPFIL3,'(A)')
     :     '# Wave Data Sigma Fit1 Fit2a Fit2b Fit2'


* setup plots
      CALL PGBEGIN(0,'?',1,1)
c     use this to prevent pauses before each plot - for fast processing
      WRITE(SCREEN,'(A$)')' Inspect every spectrum (y)? '
      ONAME='toto3'
      READ(KEYBD,'(A)')  ONAME
      IF (ONAME.EQ.'N'.OR.ONAME.EQ.'n') CALL PGASK(.FALSE.)
c     select spectrum data to save
      WRITE(SCREEN,'(A$)')' Save fits for which spectrum (e.g. 1)? '
      READ(KEYBD,*)  ISPEC
      


*  Set the rebinning parameters. Negative binsize for no rebinning.
      BINSIZ=-1
C      WRITE(SCREEN,'(A$)')' Binsize to use (A, negative for none)? '
C      READ(KEYBD,*) BINSIZ
C      WRITE(SCREEN,'(A$)')' No. of lines to fit? '
C      READ(KEYBD,*) NLINES
C      WRITE(SCREEN,'(A)')' '
      WRITE(SCREEN,'(A,F4.1,A$)')' Continuum to add (',CONTIN,')? '
      READ(KEYBD,*) CONTIN
      WRITE(SCREEN,'(A)')' '

* ask if fitting OII 4363
      WRITE(SCREEN,'(A$)')' Fit OIII 4363 (y)? '
      READ(KEYBD,'(A)')  ONAME
      I4363=1
      IF (ONAME.EQ.'N'.OR.ONAME.EQ.'n') I4363=0

**********************************************************************************
*  Loop through all the spectra in the list:
      COUNT=0
      DO WHILE (IOERR.EQ.0)

*  Read in spectrum name:
*  -change for WiggleZ: SET TRED=0 IF NOT GIVEN; all are shifted to rest
         IOERR2=0
         READ(CHLINE,*,IOSTAT=IOERR2) TNAME,TRED
         IF (IOERR2.NE.0) THEN
            READ(CHLINE,*) TNAME
            TRED=0
         END IF
         WRITE(SCREEN,'(A)') 
     +        '################################################'
         WRITE(SCREEN,'(A,x,F9.3)') TNAME,TRED

*  Read in the spectrum, bring the spectrum to its rest frame and optionally add continuum
         CALL RRFIG(TNAME, NUM, WAV, SIG, ERR, MAXSIZ)
	 DO J=1,NUM
	    WAV(J) = WAV(J)/(1.0+TRED)
            SIG(J) = SIG(J) + CONTIN
	 END DO
         WAVMIN=WAV(1)
         WAVMAX=WAV(NUM)
         WRITE(SCREEN,*) NUM,WAV(1),SIG(1),WAVMIN,WAVMAX

C Line pairs
C [OII] 3727 = 3726 + 3728.8
C [OIII] xxxx = 4958.9 + 5006.8
C Hbeta = 4861.3 + 4861.3

         IF (I4363.EQ.1) THEN
            CENTRE=4363.0
            NLINES=1
            OFLAG=0
            CALL TWOFIT(TNAME,NUM,WAV,SIG,ERR,BINSIZ,'[4363]',SFACT
     :        ,INT1,INT1E,FW1,CEN1,INT2,FW2,CEN2
     :        ,NPAR,CHI2,BIC,CORR2,ADCORR2
     :        ,FIT1,FIT21,FIT22, EW1)
	    OWF = INT1
	    OWE = INT1E
            WRITE(OUTFIL0,'(A10,I7,5F9.2,I7,2F9.2,2F9.5)')
     :        TNAME,NBIN

     :        ,INT1,INT1E,FW1,CEN1,SFACT
     :        ,NPAR,CHI2,BIC,CORR2,ADCORR2
            WRITE(SCREEN,'(A10,I7,5F9.2,I7,2F9.2,2F9.5)')
     :        TNAME,NBIN
     :        ,INT1,INT1E,FW1,CEN1,SFACT
     :        ,NPAR,CHI2,BIC,CORR2,ADCORR2
         END IF

         CENTRE=3727.0
         NLINES=1
         OFLAG=0
         CALL TWOFIT(TNAME,NUM,WAV,SIG,ERR,BINSIZ,'[OII]',SFACT
     :        ,INT1,INT1E,FW1,CEN1,INT2,FW2,CEN2
     :        ,NPAR,CHI2,BIC,CORR2,ADCORR2
     :        ,FIT1,FIT21,FIT22, EW1)
     	 OIIF = INT1
	 OIIE = INT1E
         CENTRE=3726.0
         NLINES=1
         OFLAG=2
         CALL TWOFIT(TNAME,NUM,WAV,SIG,ERR,BINSIZ,'[OII]-2',SFACT
     :        ,INT1D,INT1DE,FW1D,CEN1D,INT2D,FW2D,CEN2D
     :        ,NPARD,CHI2D,BICD,CORR2D,ADCORR2D
     :        ,FIT2,FIT21,FIT22, EW1)
	WRITE(OUTFIL1,'(A10,I7,5F9.2,I7,2F9.2,2F9.5,7F9.2
     :,I7,2F9.3,2F9.5,F8.2,F10.5)')
     :        TNAME,NBIN
     :        ,INT1,INT1E,FW1,CEN1,SFACT
     :        ,NPAR,CHI2,BIC,CORR2,ADCORR2
     :        ,INT1D,FW1D,CEN1D,INT2D,FW2D,CEN2D,INT2D/INT1D
     :        ,NPARD,CHI2D,BICD,CORR2D,ADCORR2D
     :        ,BIC-BICD,ADCORR2D-ADCORR2
        IF (I.EQ.ISPEC) THEN
         DO J=1,NBIN
            WRITE(SPFIL1,'(7F9.4)')
     :      RBWAV(J),RBDAT(J),RBSIG(J),FIT1(J),FIT21(J),FIT22(J),FIT2(J)
         END DO
        END IF

         CENTRE =4958.9
         NLINES=1
         OFLAG=1
         CALL TWOFIT(TNAME,NUM,WAV,SIG,ERR,BINSIZ,'[OIII]',SFACT
     :        ,INT1,INT1E,FW1,CEN1,INT2,FW2,CEN2
     :        ,NPAR,CHI2,BIC,CORR2,ADCORR2
     :        ,FIT1,FIT21,FIT22, EW1)
     	 OIIIF = INT1
	 OIIIE = INT1E
         CENTRE =4958.9
         NLINES=2
         OFLAG=1
         CALL TWOFIT(TNAME,NUM,WAV,SIG,ERR,BINSIZ,'[OIII]-2',SFACT
     :        ,INT1D,INT1DE,FW1D,CEN1D,INT2D,FW2D,CEN2D
     :        ,NPARD,CHI2D,BICD,CORR2D,ADCORR2D
     :        ,FIT2,FIT21,FIT22, EW1)
	WRITE(OUTFIL2,'(A10,I7,5F9.2,I7,2F9.2,2F9.5,7F9.2
     :,I7,2F9.3,2F9.5,F8.2,F10.5)')
     :        TNAME,NBIN
     :        ,INT1,INT1E,FW1,CEN1,SFACT
     :        ,NPAR,CHI2,BIC,CORR2,ADCORR2
     :        ,INT1D,FW1D,CEN1D,INT2D,FW2D,CEN2D,INT2D/INT1D
     :        ,NPARD,CHI2D,BICD,CORR2D,ADCORR2D
     :        ,BIC-BICD,ADCORR2D-ADCORR2
        IF (I.EQ.ISPEC) THEN
         DO J=1,NBIN
            WRITE(SPFIL2,'(7F9.4)')
     :      RBWAV(J),RBDAT(J),RBSIG(J),FIT1(J),FIT21(J),FIT22(J),FIT2(J)
         END DO
         END IF

         CENTRE =4861.3
         NLINES=1
         OFLAG=0
         CALL TWOFIT(TNAME,NUM,WAV,SIG,ERR,BINSIZ,'Hbeta',SFACT
     :        ,INT1,INT1E,FW1,CEN1,INT2,FW2,CEN2
     :        ,NPAR,CHI2,BIC,CORR2,ADCORR2
     :        ,FIT1,FIT21,FIT22, EW1)
         CENTRE =4861.3
         NLINES=2
         OFLAG=0
	 HBF = INT1
	 HBE = INT1E
	 EWB = EW1
*	 WRITE(SCREEN, 'A'
         CALL TWOFIT(TNAME,NUM,WAV,SIG,ERR,BINSIZ,'Hbeta-2',SFACT
     :        ,INT1D,INT1DE,FW1D,CEN1D,INT2D,FW2D,CEN2D
     :        ,NPARD,CHI2D,BICD,CORR2D,ADCORR2D
     :        ,FIT2,FIT21,FIT22, EW1)
         WRITE(OUTFIL3,'(A10,I7,5F9.2,I7,2F9.2,2F9.5,7F9.2
     :,I7,2F9.3,2F9.5,F8.2,F10.5)')
     :        TNAME,NBIN
     :        ,INT1,INT1E,FW1,CEN1,SFACT
     :        ,NPAR,CHI2,BIC,CORR2,ADCORR2
     :        ,INT1D,FW1D,CEN1D,INT2D,FW2D,CEN2D,INT2D/INT1D
     :        ,NPARD,CHI2D,BICD,CORR2D,ADCORR2D
     :        ,BIC-BICD,ADCORR2D-ADCORR2
        IF (I.EQ.ISPEC) THEN
         DO J=1,NBIN
            WRITE(SPFIL3,'(7F9.4)')
     :      RBWAV(J),RBDAT(J),RBSIG(J),FIT1(J),FIT21(J),FIT22(J),FIT2(J)
         END DO
         END IF



         CENTRE =4341.0
         NLINES=1
         OFLAG=0
         CALL TWOFIT(TNAME,NUM,WAV,SIG,ERR,BINSIZ,'Hgamma',SFACT
     :        ,INT1,INT1E,FW1,CEN1,INT2,FW2,CEN2
     :        ,NPAR,CHI2,BIC,CORR2,ADCORR2
     :        ,FIT1,FIT21,FIT22, EW1)
         CENTRE =4341.0
         NLINES=2
         OFLAG=0
	 HGF = INT1
	 HGAE = INT1E
*	 TNAME = SPECNAME
*         CALL TWOFIT(TNAME,NUM,WAV,SIG,ERR,BINSIZ,'Hgamma-2',SFACT
*    :        ,INT1D,INT1DE,FW1D,CEN1D,INT2D,FW2D,CEN2D
*    :        ,NPARD,CHI2D,BICD,CORR2D,ADCORR2D
*    :        ,FIT2,FIT21,FIT22, EW1)
         WRITE(OUTFIL4,'(A10,I7,5F9.2,I7,2F9.2,2F9.5,7F9.2
     :,I7,2F9.3,2F9.5,F8.2,F10.5)')
     :        TNAME,NBIN
     :        ,INT1,INT1E,FW1,CEN1,SFACT
     :        ,NPAR,CHI2,BIC,CORR2,ADCORR2
     :        ,INT1D,FW1D,CEN1D,INT2D,FW2D,CEN2D,INT2D/INT1D
     :        ,NPARD,CHI2D,BICD,CORR2D,ADCORR2D
     :        ,BIC-BICD,ADCORR2D-ADCORR2
        IF (I.EQ.ISPEC) THEN
         DO J=1,NBIN
            WRITE(SPFIL3,'(7F9.4)')
     :      RBWAV(J),RBDAT(J),RBSIG(J),FIT1(J),FIT21(J),FIT22(J),FIT2(J)
         END DO
         END IF
	 WRITE(OUTFIL5,'(A10,11F9.2)')
     :        TNAME,OWF
     :        ,OWE,OIIF,OIIE,OIIIF,OIIIE
     :        ,HBF,HBE,EWB,HGF,HGAE


         READ(INFILE,'(A)',IOSTAT=IOERR) CHLINE
      ENDDO

      CLOSE(INFILE)
      CLOSE(OUTFIL0)
      CLOSE(OUTFIL1)
      CLOSE(OUTFIL2)
      CLOSE(OUTFIL3)
      CLOSE(OUTFIL4)
      CLOSE(OUTFIL5)

      CLOSE(SPFIL1)
      CLOSE(SPFIL2)
      CLOSE(SPFIL3)

      CALL PGEND
 
      END


      SUBROUTINE TWOFIT(TNAME,NUM,WAV,DAT,SIG,BINSIZ
     :     ,LNAME,SFACTOR,INT1,INT1E,FW1,CEN1,INT2,FW2,CEN2
     :     ,NPAR,CHI2,BIC,CORR2,ADCORR2
     :     ,MODEL,MODEL1,MODEL2,EW1)
*-
* after GENFIT by P.Francis, modified by Eileen
*
*  Reads in a list of spectra, sends them to their rest frame,
* rebins part of them, normalised by its mean. This then has
* the mean spectrum subtracted from it.
*
* A 2-G model is then fit to the spectra, and the parameters output
*
* Input:- List of spectra and redshifts.
*
* Calls:- RRFIG, AMOEBA
* A file listing the fit parameters
*
* Paul Francis, Melbourne, 4th Jan 1994. Modified from the PCA
* program spca.f
*-

      IMPLICIT NONE

      INTEGER J, COLMAX, RMAX, K, NLINES, STORENL
      INTEGER SCREEN, KEYBD
      INTEGER NPAR, MP, NP, MAXROOT
      INTEGER MAXSIZ, NBIN, NUM, ITER, N, ITMAX

c      PARAMETER (COLMAX=800, RMAX=600, SCREEN=6, ITMAX=1500,
c     : KEYBD=5, MAXSIZ=8000, MP=40, NP=40)
c      PARAMETER (COLMAX=800, RMAX=600, SCREEN=6, ITMAX=6500,
      PARAMETER (COLMAX=800, RMAX=600, SCREEN=6, ITMAX=20000,
     : KEYBD=5, MAXSIZ=8000, MP=40, NP=40, MAXROOT=500)

      REAL FUNK, FUNK2, FUNC1, RTBIS
      REAL P(MP,NP), PAR(MP), TEMPAR(MP), Y(MP)
      REAL WAV(MAXSIZ),   DAT(MAXSIZ),   SIG(MAXSIZ), FACTOR(MP)
      REAL RBWAV(MAXSIZ), RBDAT(MAXSIZ), RBSIG(MAXSIZ)
      REAL START, END, BINSIZ, TEMP, HIY
      REAL XARRAY(MAXSIZ),XMEAN,XRMS,RMS,SFACTOR
      REAL MEAN, FTOL, FMIN, MODEL(MAXSIZ)
      REAL PMIN(MP), MODEL1(MAXSIZ), MODEL2(MAXSIZ)
      REAL ALFA,WAVE,FLUX,PI,SCALE,RANGE,CENTRE
      REAL FW1,CEN1,FW2,CEN2,INT1,INT1E,INT2,FLUX0,CHI2,BIC
      REAL YMEAN,YRMS,SUMXY,XYMEAN,CORR2,ADCORR2,EW1
      REAL CMIN,EMIN,ETEST,ELOWER,EUPPER,EPLUS,EMINUS,ECHANGE
      INTEGER ISTART,IEND,IX1,IX2,NX,OUNIT,OFLAG

      CHARACTER  LNAME*(*),TNAME*(*),WARNING*(12)
      character*80 TITLE

      COMMON /SPECTRA/RBWAV,RBDAT,RBSIG,NBIN,NLINES,CENTRE,OFLAG
      COMMON /FITTING/PAR

      data pi/3.141592/

      EXTERNAL FUNK
      EXTERNAL FUNK2
      EXTERNAL FUNC1
      EXTERNAL RTBIS

* setup fractional tolorance

      FTOL = 5.0E-7
      FTOL = 5.0E-6
      IF (LNAME(1:5).eq.'Hbeta') FTOL  = 5.0E-5

      RANGE = 205.0
      RANGE = 55.0
c      RANGE = 36.0

*  SET the wavelength range in WAVELENGTH units (a bit more for H-beta to pass OIII)
	START=CENTRE-RANGE
c        IF (OFLAG.EQ.1) START=START+(5006.8 - 4958.9)/2.0
        IF (OFLAG.EQ.1) START=START+((5006.8 - 4958.9)/2.0) -8.0
	IF(ABS(CENTRE-4363.).LT.10.) START=4348.0

	END=START+RANGE*2.0
c	IF(ABS(CENTRE-4861.).LT.20.) END=END+RANGE/2.0

C     TEST IF REGION OF INTEREST IS IN THE PROVIDED SPECTRUM
	IF(WAV(1).LE.START.AND.WAV(NUM).GE.END) THEN


c     Rebin if BINSIZ > 0, otherwise just copy data directly
c     in the latter case, also apply rough scaling to noise array
           IF (BINSIZ.GT.0.0) THEN
*     Prepare the rebinned wavelength scale
              NBIN = INT((END-START)/BINSIZ)
              TEMP=START
              DO J=1,NBIN
                 RBWAV(J)=TEMP
                 TEMP=TEMP+BINSIZ
              END DO
*     Rebin data and sigma
              CALL REBIN(NUM,NBIN,WAV,DAT,RBWAV,RBDAT)
              CALL REBIN(NUM,NBIN,WAV,SIG,RBWAV,RBSIG)
           ELSE
*     Find start and end points in original data
              ISTART=1
              DO WHILE (WAV(ISTART).LT.START.AND.ISTART.LT.NUM)
                 ISTART=ISTART+1
              END DO
              IF (ISTART.GE.NUM) STOP 'Cannot find start index'
              IEND=NUM
              DO WHILE (WAV(IEND).GT.END.AND.IEND.GT.0)
                 IEND=IEND-1
              END DO
              IF (IEND.LE.0) STOP 'Cannot find end index'
*     Copy selected data from original arrays
              NBIN=0
              DO J=ISTART,IEND
                 NBIN=NBIN+1
                 RBWAV(NBIN)=WAV(J)
                 RBDAT(NBIN)=DAT(J)
                 RBSIG(NBIN)=SIG(J)
              END DO
c     Scale noise array to match approximate RMS of nearby data
c     -normally above selected region but below for H-beta
              IF(ABS(CENTRE-4861.).LT.20.) THEN
                 IX1=ISTART-INT((IEND-ISTART)/2)
                 IX2=ISTART
              ELSE
                 IX1=IEND
                 IX2=IEND+INT((IEND-ISTART)/2)
              END IF

              IF (IX1.GE.1.AND.IX2.LE.NUM) THEN
                 NX=0
                 DO J=IX1,IX2
                    NX=NX+1
                    XARRAY(NX) = DAT(J)
                 END DO
                 CALL MEAN_RMS(XARRAY,NX,XMEAN,XRMS)
                 RMS = XRMS
                 NX=0
                 DO J=IX1,IX2
                    NX=NX+1
                    XARRAY(NX) = SIG(J)
                 END DO
                 CALL MEAN_RMS(XARRAY,NX,XMEAN,XRMS)
                 SFACTOR=RMS/XMEAN
                 DO J=1, NBIN
                    RBSIG(J) = RBSIG(J)*SFACTOR
                 END DO
c                 WRITE(SCREEN,*) 'SCALE:',ISTART,IEND,IX1,IX2,RMS,XMEAN
              ELSE
                 STOP 'SCALE: LIMITS OUT OF RANGE'
              END IF

           END IF



*  Find the mean flux, scale data AND uncertainty by it; 
c
c     DO NOT RESCALE so we can do line ratios (M.Drinkwater 21/12/15)
c
           MEAN = 0.0
           HIY = -1000.0
           DO J=1, NBIN
              MEAN = MEAN + RBDAT(J)
	   END DO
           MEAN = MEAN/REAL(NBIN)
           DO J=1, NBIN
c              RBDAT(J) = RBDAT(J)/MEAN
c              IF (OUNIT.GT.0) WRITE(OUNIT-2,*) RBDAT(J)
c              RBSIG(J) = RBSIG(J)/MEAN
              HIY = MAX(HIY,RBDAT(J))
	   END DO
c     PLOT data and uncertainty
	   HIY = HIY * 1.5
	   CALL PGSLS(1)
	   CALL PGENV(RBWAV(1),RBWAV(NBIN),0.0,HIY,0,0)
           TITLE=LNAME//' '//TNAME
	   CALL PGLABEL('Rest Wavelength','Flux',TITLE)
c     default colour
           CALL PGSLS(1)
	   CALL PGLINE(NBIN,RBWAV,RBDAT)
           CALL PGSLS(3)
	   CALL PGLINE(NBIN,RBWAV,RBSIG)
           CALL PGSLS(1)



           WRITE(SCREEN,'(A,2A12,A)')
     :         '------------- Fitting ',TNAME
     :          ,LNAME,'  ----------------'
           WRITE(SCREEN,*) 'SCALE:',ISTART,IEND,IX1,IX2,RMS,XMEAN

c  OLD fitting parameters:
c	PAR(1) is the power law scale (wavelength)
c	PAR(2) is the power law alpha
c	PAR(3) is central wavelength for the lines
c	PAR(4) is FWHM for the broad gaussians
c	PAR(5) is EW for the broad gaussian lines
c	PAR(6) is FWHM for the narrow gaussians
c	PAR(7) is EW for the narrow gaussian lines
c	PAR(8) (NLINES) =1 for broad only =2 for broad + narrow

c  NEW fitting parameters:
c	PAR(1) is the inverse continuum slope (wavelength)
c	PAR(2) is the continuum level
c	PAR(3) is FWHM for the broad gaussians (km/s)
c	PAR(4) is EW for the broad gaussian lines
c	PAR(5) is FWHM for the narrow gaussians  (km/s)
c	PAR(6) is EW for the narrow gaussian lines
c -optional for NLINES>2 fitting 2 components to OIII doublet = 4 lines in total
c	PAR(7) is centre for the broad  gaussian lines
c	PAR(8) is centre for the narrow gaussian lines

c  NEWER fitting parameters: SINGLE LINE (+flag to add OIII doublet)
c	PAR(1) is the continuum level
c	PAR(2) is EW for the narrow gaussians (km/s)
c	PAR(3) is FWHM for the narrow gaussian lines
c	PAR(4) is centre for the narrow  gaussian lines

c  NEWER fitting parameters: DOUBLE LINE (+flag to add OIII doublet)
c	PAR(1) is the continuum level
c	PAR(2) is EW for the narrow gaussians (km/s)
c	PAR(3) is FWHM for the narrow gaussian lines
c	PAR(4) is centre for the narrow  gaussian lines
c	PAR(5) is EW for the broad gaussians (km/s)
c	PAR(6) is FWHM for the broad gaussian lines
c	PAR(7) is centre for the broad  gaussian lines

c Set up starting simplex
c NOTE NLINES,CENTRE,OFLAG ALREADY SET BEFORE CALLING THIS ROUTINE
c OII:   NL=1, NL=1 +OFLAG=2 (fixing separation, FW only)
c Hbeta: NL=1, NL=2
c OIII:  NL=1, NL=2 +OFLAG=1 (fixing separation, FW and EW)

           PAR(1) = 0.75
           PAR(1) = 0.8
           PAR(1) = MEAN
           PAR(2) = 10.0
           PAR(3) = 250.0
           PAR(4) = CENTRE
	   PAR(5) = 10.0
	   PAR(6) = 250.0
           PAR(7) = CENTRE
	   NPAR=4
	   IF (NLINES.EQ.2) NPAR=7
	   IF (OFLAG.EQ.2) NPAR=5

           IF (LNAME(1:5).eq.'Hbeta'.AND.NLINES.EQ.2) THEN
              PAR(1) = MEAN
              PAR(2) = 10.0
              PAR(3) = 250.0
              PAR(5) = -4.0
              PAR(6) = 900.0
           END IF

	   IF (LNAME(1:5).eq.'[OIII'.AND.NLINES.EQ.2) THEN
              PAR(2) = 7.0
              PAR(3) = 250.0
              PAR(5) = 1.5
              PAR(6) = 900.0
              PAR(7) = CENTRE-2.0
	   END IF 

* Plot starting model
  	  ALFA = PAR(1)
c	  IF (NLINES.GT.1) THEN
  	  DO N = 1, NBIN
	    WAVE = RBWAV(N)
	    FLUX = ALFA
	    CALL LINEPROF(PAR,FLUX,WAVE)
            MODEL(N) = FLUX
c            WRITE(SCREEN,*) N,WAVE,ALFA,FLUX
	  END DO
c     blue
          CALL PGSCI(4)
	  CALL PGSLS(2)
	  CALL PGLINE(NBIN,RBWAV,MODEL)
          CALL PGSCI(1)
c	  END IF

* Use variable starting factors in each dimension
           FACTOR(1) = 0.02
           FACTOR(1) = 0.2
           FACTOR(1) = 0.2
	   FACTOR(2) = 2.0 
           FACTOR(3) = 10.0
           FACTOR(4) = 0.2
	   FACTOR(5) = 2.0
	   FACTOR(6) = 10.0
	   FACTOR(7) = 0.1

           DO K=1, (NPAR+1)
               DO J=1, NPAR
                  P(K,J) = PAR(J)
                  IF (K .EQ. J) THEN
                     P(K,J) = P(K,J) + FACTOR(J)
                  END IF
                  TEMPAR(J) = P(K,J)
	       END DO
               Y(K) = FUNK(TEMPAR)
	   END DO
           CALL AMOEBA(P,Y,MP,NP,NPAR,FTOL,FUNK,ITER,ITMAX)

* If converged chose output, compute the model
           WRITE(SCREEN,'(A,I9)') 'Number of iterations = ',ITER
	   IF (ITER.LT.ITMAX) THEN
           FMIN = 1.0E15
           DO J=1,NPAR+1
               DO K=1, NPAR
                  PAR(K) = P(1,K)
	       END DO
               IF (FUNK(PAR) .LT. FMIN) THEN
                  FMIN = FUNK(PAR)
                  DO K=1, NPAR
                     PMIN(K) = PAR(K)
	       END DO
               END IF
	   END DO
           DO N=1,NPAR
              PAR(N) = PMIN(N)
           END DO
           ALFA = PAR(1)
           DO N = 1, NBIN
              WAVE = RBWAV(N)
              FLUX = ALFA
              CALL LINEPROF(PAR,FLUX,WAVE)
              MODEL(N) = FLUX
c              IF (OUNIT.GT.0) WRITE(OUNIT,*) FLUX
           END DO

c     calculate statistics on model
           CHI2 = FUNK2(PAR)
           BIC = CHI2 + REAL(NPAR)*LOG(REAL(NBIN))

           CALL MEAN_RMS(RBDAT,NBIN,XMEAN,XRMS)
           CALL MEAN_RMS(MODEL,NBIN,YMEAN,YRMS)
           SUMXY=0.0
           DO N = 1, NBIN
              SUMXY = SUMXY + MODEL(N)*RBDAT(N)
           END DO
           XYMEAN=SUMXY/REAL(NBIN)
           CORR2 = ((XYMEAN - XMEAN*YMEAN)/(XRMS*YRMS))**2
           ADCORR2=CORR2-((1-CORR2)*(REAL(NPAR)/(REAL(NBIN-NPAR)-1)))


c     plot combined model
c     green
           CALL PGSCI(3)
           CALL PGLINE(NBIN,RBWAV,MODEL)
           CALL PGSCI(1)
           WRITE(SCREEN,'(A)') 'plotted model '


c     if two-line fit then plot components
c     make individual components by repeating with 1 line and taking difference. Crude!
           IF (NLINES.GT.1) THEN
              STORENL=NLINES
              NLINES=NLINES-1
              DO N = 1, NBIN
                 WAVE = RBWAV(N)
                 FLUX = ALFA
                 FLUX0 = FLUX
                 CALL LINEPROF(PAR,FLUX,WAVE)
                 MODEL1(N) = FLUX
                 MODEL2(N) = MODEL(N) - FLUX + FLUX0
              END DO
              NLINES=STORENL
              CALL PGSLS(2)
c     red
              CALL PGSCI(2)
              CALL PGLINE(NBIN,RBWAV,MODEL1)
c     magenta
              CALL PGSCI(6)
              CALL PGLINE(NBIN,RBWAV,MODEL2)
              CALL PGSCI(1)
           END IF


c     return selected parameters for lines
          FLUX=  ALFA
c     -line 1
          INT1 = FLUX * PAR(2)
	  FW1=PAR(3)
	  CEN1=PAR(4)
	  EW1 = PAR(2)
c     -line 2
          INT2 = FLUX * PAR(5)
	  FW2=PAR(6)
	  CEN2=PAR(7)
          WARNING=' '
c     find uncertainty range in line flux based on 1-sigma range of EW=PAR(2)
c     this is found by varying the parameter up and down until the Chi-2 difference
c     is 3.5 above the best fit value. This only applies to the single-line fits when
c     we are only fitting 3 parameters. The target chi-2 value is stored in PAR(40).
c     See: Wall & Jenkins 2008 Practical Statistics for Astronomers, Table 6.1
           WRITE(SCREEN,'(A)') 'searching range '
           CMIN = FUNK2(PAR)
           PAR(40) = CMIN + 3.5
           EMIN = PAR(2)
           ETEST = EMIN
           ECHANGE = ABS(ETEST/2)
           WRITE(SCREEN,'(A)') 'searching UP '
           K=0
           DO WHILE (FUNC1(ETEST).LE.0.AND.K.LT.MAXROOT)
              K=K+1
              IF (K.LT.20) WRITE(SCREEN,'(4F14.5)')
     +             CMIN,PAR(40),ETEST,FUNC1(ETEST)
              ETEST = ETEST + ECHANGE
           END DO
           IF (K.GE.MAXROOT)WRITE(SCREEN,'(A)') '**NO ROOT**'
           EUPPER = ETEST
           ETEST = EMIN
           ECHANGE = ABS(ETEST/2)
           WRITE(SCREEN,'(A)') 'searching DOWN '
           K=0
           DO WHILE (FUNC1(ETEST).LE.0.AND.K.LT.MAXROOT)
              K=K+1
              IF (K.LT.20) WRITE(SCREEN,'(4F14.5)')
     +             CMIN,PAR(40),ETEST,FUNC1(ETEST)
              ETEST = ETEST - ECHANGE
           END DO
           IF (K.GE.MAXROOT)WRITE(SCREEN,'(A)') '**NO ROOT**'
           ELOWER = ETEST
           WRITE(SCREEN,'(A,F12.3)') 'Original min chi2: ', CMIN
           WRITE(SCREEN,'(A,2F12.3)') 'Lower: ', ELOWER, FUNC1(ELOWER)
           WRITE(SCREEN,'(A,2F12.3)') 'Min:   ', EMIN, FUNC1(EMIN)
           WRITE(SCREEN,'(A,2F12.3)') 'Upper: ', EUPPER, FUNC1(EUPPER)
           FTOL=1E-5
           EPLUS = RTBIS(FUNC1,EMIN,EUPPER,FTOL)
           EMINUS = RTBIS(FUNC1,ELOWER,EMIN,FTOL)
           WRITE(SCREEN,'(A)')
     +          '      EW          Chi2-3.5    Flux      Flux-best'
           WRITE(SCREEN,'(4F12.3)')
     +          EPLUS, FUNC1(EPLUS),FLUX*EPLUS,FLUX*(EPLUS-EMIN)
           WRITE(SCREEN,'(3F12.3)') EMIN, FUNC1(EMIN),FLUX*EMIN
           WRITE(SCREEN,'(4F12.3)') 
     +          EMINUS, FUNC1(EMINUS),FLUX*EMINUS,FLUX*(EMINUS-EMIN)
           INT1E = FLUX*(EPLUS-EMINUS)/2.0
           WRITE(SCREEN,'(A,F12.3)') 'Flux error= ',INT1E

	ELSE
	  FW1=0.0
	  CEN1=0.0
	  FW2=0.0
	  CEN2=0.0
          INT1=-99
          INT2=-99
          WARNING=' ***FAILED'
	END IF

        WRITE(SCREEN,'(A)')
     :       '    ITER   NBIN   NPAR'//
     :       '      PAR(1),      PAR(2),      PAR(3),'//
     :       '      PAR(4),        INT1,      INT1E,      PAR(5),'//
     :       '      PAR(6),      PAR(7),        INT2,'//
     :       '   INT2/INT1,        CHI2,          BIC'//
     :       '       CORR2,     ADCORR2'
        WRITE(SCREEN,'(X,3I7,15F13.6,A)')
     :       ITER,NBIN,NPAR,PAR(1),PAR(2),PAR(3),PAR(4)
     :       ,INT1,INT1E,PAR(5),PAR(6),PAR(7),INT2,INT2/INT1
     :       ,CHI2,BIC,CORR2,ADCORR2,WARNING

      ELSE
	  FW1=0.0
	  CEN1=0.0
	  FW2=0.0
	  CEN2=0.0
          INT1=-99
          INT2=-99
          WARNING=' ***OUT OF RANGE'
          WRITE(SCREEN,'(30X,A)') ' ***OUT OF RANGE'
      END IF

      END

 
	REAL FUNCTION FUNK(PAR)
c     calculates RMS of (model spectrum - data)

	IMPLICIT NONE
	INTEGER N
	INTEGER MAXSIZ
        INTEGER I, NBIN, NLINES,OFLAG
	PARAMETER (MAXSIZ=8000)
	REAL PAR(*)
	REAL WAVE
        REAL FMAX,CENTRE
	REAL PI, SCALE, DAT(MAXSIZ),RBSIG(MAXSIZ)
	REAL FLUX, ALFA, WAVARRAY(MAXSIZ), DIFF
 	COMMON /SPECTRA/WAVARRAY,DAT,RBSIG,NBIN,NLINES,CENTRE,OFLAG
	DATA PI/3.14159265/

        FMAX = -1.
 	ALFA = PAR(1)
        I = 0
        DIFF=0.0

	DO 80 N = 1, NBIN
	  WAVE = WAVARRAY(N)
c	    FLUX = (WAVE/SCALE)**(-(2+ALFA))
	    FLUX =  ALFA
	  CALL LINEPROF(PAR,FLUX,WAVE)
          DIFF = DIFF + (FLUX - DAT(N))**2
80	CONTINUE
        FUNK = SQRT(DIFF/REAL(NBIN))

	END


	REAL FUNCTION FUNK2(PAR)
c     calculates Chi^2 of (model spectrum - data)

	IMPLICIT NONE
	INTEGER N
	INTEGER MAXSIZ
        INTEGER I, NBIN, NLINES,OFLAG
	PARAMETER (MAXSIZ=8000)
	REAL PAR(*)
	REAL WAVE
        REAL FMAX,CENTRE
	REAL PI, SCALE, DAT(MAXSIZ),RBSIG(MAXSIZ)
	REAL FLUX, ALFA, WAVARRAY(MAXSIZ), DIFF
 	COMMON /SPECTRA/WAVARRAY,DAT,RBSIG,NBIN,NLINES,CENTRE,OFLAG
	DATA PI/3.14159265/

        FMAX = -1.
 	ALFA = PAR(1)
        I = 0
        DIFF=0.0

	DO 80 N = 1, NBIN
	  WAVE = WAVARRAY(N)
c	    FLUX = (WAVE/SCALE)**(-(2+ALFA))
	    FLUX =  ALFA
	  CALL LINEPROF(PAR,FLUX,WAVE)
c        write(6,*) 'FUNK2: ',N,WAVE,FLUX,DAT(N),RBSIG(N),DIFF
          DIFF = DIFF + ((FLUX - DAT(N))/RBSIG(N))**2
80	CONTINUE
        FUNK2 = DIFF
c        write(6,*) 'FUNK2: ',DIFF,NBIN

	END



	REAL FUNCTION FUNC1(XTEST)
c     calculates Chi^2(model spectrum - data) - PAR(40)
c     BUT only varies one parameter, X=PAR(2)
c     other values of PAR from common block
c     used to find value of PAR(2) where Chi^2 = PAR(40)

	IMPLICIT NONE
	INTEGER N
	INTEGER MAXSIZ
        INTEGER I, NBIN, NLINES,OFLAG
	PARAMETER (MAXSIZ=8000)
	REAL PAR(40)
	REAL WAVE,XTEST
        REAL FMAX,CENTRE
	REAL PI, SCALE, DAT(MAXSIZ),RBSIG(MAXSIZ)
	REAL FLUX, ALFA, WAVARRAY(MAXSIZ), DIFF
 	COMMON /SPECTRA/WAVARRAY,DAT,RBSIG,NBIN,NLINES,CENTRE,OFLAG
        COMMON /FITTING/PAR
	DATA PI/3.14159265/

        FMAX = -1.
 	ALFA = PAR(1)
        PAR(2) = XTEST
        I = 0
        DIFF=0.0

	DO N = 1, NBIN
	   WAVE = WAVARRAY(N)
c	   FLUX = (WAVE/SCALE)**(-(2+ALFA))
	   FLUX =  ALFA
           CALL LINEPROF(PAR,FLUX,WAVE)
c     write(6,*) 'FUNC1: ',N,WAVE,FLUX,DAT(N),RBSIG(N),DIFF
           DIFF = DIFF + ((FLUX - DAT(N))/RBSIG(N))**2
        ENDDO
        FUNC1 = DIFF - PAR(40)
c        write(6,*) 'FUNC1: ',DIFF,NBIN

	END




      SUBROUTINE AMOEBA(P,Y,MP,NP,NDIM,FTOL,FUNK,ITER,ITMAX)
*-
* Downhill simplex minimisation routine
* From Numerical Recipies, p292
*-

      PARAMETER (NMAX=20,ALPHA=1.0,BETA=0.5,GAMMA=2.0)
      DIMENSION P(MP,NP),Y(MP),PR(NMAX),PRR(NMAX),PBAR(NMAX)
      MPTS=NDIM+1
      ITER=0

* Find the highest, second highest and lowest (best) point

1     ILO=1
      IF (Y(1) .GT. Y(2)) THEN
         IHI=1
         INHI=2
      ELSE
         IHI=2
         INHI=1
      END IF
      DO 11 I=1, MPTS
         IF(Y(I).LT.Y(ILO)) ILO=I
         IF(Y(I).GT.Y(IHI)) THEN
            INHI=IHI
            IHI=I
         ELSE IF(Y(I).GT.Y(INHI))THEN
            IF(I.NE.IHI) INHI=I
         ENDIF
11    CONTINUE

* Compute the fractional range from highest to lowest

      RTOL=2.*ABS(Y(IHI)-Y(ILO))/(ABS(Y(IHI))+ABS(Y(ILO)))
      IF(RTOL.LT.FTOL)RETURN
      IF(ITER.EQ.ITMAX) THEN
C       PAUSE 'Amoeba exceeding maximum iterations'
	Print *,'Amoeba exceeding maximum iterations; NO FIT'
	RETURN
      END IF
      ITER=ITER+1

* Start the new iteration

      DO 12 J=1,NDIM
         PBAR(J) = 0.
12    CONTINUE

* Find the centre of the face

      DO 14 I=1,MPTS
         IF (I.NE.IHI)THEN
            DO 13 J=1,NDIM
               PBAR(J)=PBAR(J)+P(I,J)
13          CONTINUE
         END IF
14    CONTINUE

* Reflect the simplex by factor ALPHA

      DO 15 J=1,NDIM
         PBAR(J) = PBAR(J)/NDIM
         PR(J) = (1.+ALPHA)*PBAR(J)-ALPHA*P(IHI,J)
15    CONTINUE

* Evaluate the new vertex

      YPR = FUNK(PR)

* Case 1: the new point is best of all, so try going GAMMA further!

      IF(YPR.LE.Y(ILO))THEN
         DO 16 J=1,NDIM
            PRR(J)=GAMMA*PR(J)+(1.-GAMMA)*PBAR(J)
16       CONTINUE
         YPRR=FUNK(PRR)   
* If it worked...
         IF(YPRR.LT.Y(ILO))THEN
            DO 17 J=1,NDIM
               P(IHI,J)=PRR(J)
17          CONTINUE
         Y(IHI)=YPRR
* If it didn't, go back to the original reflected point
         ELSE
            DO 18 J=1,NDIM
               P(IHI,J)=PR(J)
18          CONTINUE
            Y(IHI)=YPR
         END IF

* Case 2: the reflected point is worse than the second-highest

      ELSE IF(YPR.GE.Y(INHI))THEN
* If it is better than the worst point, the replace the worst
         IF(YPR.LT.Y(IHI))THEN
            DO 19 J=1,NDIM
               P(IHI,J)=PR(J)
19          CONTINUE
            Y(IHI)=YPR
         ENDIF
* But look for a lower intermediate point
         DO 21 J=1,NDIM
            PRR(J)=BETA*P(IHI,J)+(1.-BETA)*PBAR(J)
21       CONTINUE
         YPRR=FUNK(PRR)
* If it helps
         IF(YPRR.LT.Y(IHI))THEN
         DO 22 J=1,NDIM
            P(IHI,J)=PRR(J)
22       CONTINUE
         Y(IHI)=YPRR
* All else fails, so shrink the amoeba around its best point
         ELSE
            DO 24 I=1, MPTS
               IF(I.NE.ILO) THEN
                  DO 23 J=1,NDIM
                     PR(J)=0.5*(P(I,J)+P(ILO,J))
                     P(I,J)=PR(J)
23                CONTINUE
                  Y(I)=FUNK(PR)
               ENDIF
24          CONTINUE
         ENDIF

* Case 3: the reflection gave a middling point - use it

      ELSE
         DO 25 J=1,NDIM
            P(IHI,J)=PR(J)
25       CONTINUE
         Y(IHI)=YPR
      END IF

      GOTO 1

      END




      SUBROUTINE REBIN(NUM,NUMB,WAVI,DATI,WAVB,RBDAT)
*-
* REBIN
*   Rebins the spectrum in WAVI, DATI, to the wavelength scale
* of WAVB, putting the result in RBDAT.
*
* Modified from REBLIM to use single precision, and not to bother with
* starts and finishes of the data region, and errors. Feb 1990.
* Modified from REBIN Sept 1989. Modified to scale errors by the change
* in binsize, Feb 1990.
*-
      INTEGER NUM,NUMB,K,I,J,LOLIM,HILIM
      REAL WAVI(NUM),DATI(NUM),
     :     WAVB(NUMB),RBDAT(NUMB),A,B,C,D,
     :     INI,INB,STI,STB

* Find start and interval of both spectra.

      STI=WAVI(1)
      INI=(WAVI(NUM)-WAVI(1))/(NUM-1)
      STB=WAVB(1)
      INB=(WAVB(NUMB)-WAVB(1))/(NUMB-1)

* Prepare array to hold rebinned first spectrum.

      DO 5 K=1,NUMB
        RBDAT(K)=0.0
5     CONTINUE

* Perform rebinning.

      DO 7 I=1,NUM

        LOLIM=INT((WAVI(I)-0.5*(INI+INB)-STB)/INB)
        HILIM=INT((WAVI(I)+0.5*(INI+INB)-STB)/INB)+2

        DO 6 J=LOLIM,HILIM
          IF ((J.GE.1).AND.(J.LE.NUMB)) THEN
            A=WAVI(I)-0.5*INI
            B=WAVI(I)+0.5*INI
            C=WAVB(J)-0.5*INB
            D=WAVB(J)+0.5*INB

            IF ((C.LE.B).AND.(A.LE.D)) THEN
              IF (C.LT.A) THEN
                IF (D.LT.B) THEN
                  RBDAT(J)=RBDAT(J)+(D-A)*DATI(I)/(B-A)
                ELSE
                  RBDAT(J)=RBDAT(J)+DATI(I)
                ENDIF
              ELSE
                IF (D.LT.B) THEN
                  RBDAT(J)=RBDAT(J)+(D-C)*DATI(I)/(B-A)
                ELSE
                  RBDAT(J)=RBDAT(J)+(B-C)*DATI(I)/(B-A)
                ENDIF
              ENDIF
            ENDIF
          ENDIF
6       CONTINUE
7     CONTINUE

* Eliminate edge effects.

      DO 20 I=1,NUMB
         IF (((WAVB(I)-0.5*INB).LT.(STI+0.5*INI)).OR.
     :   ((WAVB(I)+0.5*INB).GT.(STI+NUM*INI-0.5*INI))) THEN
            RBDAT(I)=0.0
         ENDIF
20    CONTINUE

      END




       SUBROUTINE RRFIG(TNAME, NUM, WAV, SIG, ERRARR, MAXSIZ)
*-
* RRFIG
*   Given a character string TNAME, this subroutine reads into the
* arrays WAV, SIG and ERR the IUE ascii file data of the spectrum
* TNAME.mrg
*
* Paul Francis, University of Melbourne, Dec 1993
*-

      INTEGER I, FILNUM, NUM
      PARAMETER (FILNUM=9)
      REAL WAV(MAXSIZ), SIG(MAXSIZ), ERRARR(MAXSIZ)
      REAL WAVE, FLUX, SIGMA
      Character  TNAME*20

* Open the file for read access only.

      OPEN(FILNUM, FILE=TNAME, STATUS='OLD')

* Read in the spectrum, discarding zero wavelength rows.

      NUM = 0
      DO 10 I=1,10000

         READ(FILNUM,*,END=20) WAVE, FLUX, SIGMA

         IF (WAVE .GT. 1000.0) THEN

            NUM = NUM + 1
            WAV(NUM) = WAVE
            SIG(NUM) = FLUX
            ERRARR(NUM) = SIGMA

         END IF

10    CONTINUE

* Close the file.

20    CLOSE(FILNUM)

      END


	SUBROUTINE LINEPROF(PAR,flux,wave)

c       CALCULATING THE EMISSION LINE PROFILES
c       FOR DOUBLE GAUSSIANS              
 
	implicit none
	real linen, integ, fwhm, sigma, h, flux
	real line, wave, integn, fwhmn, sigman, hn
	real wavebar, wavebarn, wavebar2, wavebarn2, doublet
        real PAR(*), pi, wid, widn, wmid, wmidn
	INTEGER NLINES
	data pi/3.14159265/
c
        INTEGER MAXSIZ, NBIN,OFLAG
        PARAMETER (MAXSIZ=8000)
        REAL CENTRE
	REAL DAT(MAXSIZ)
	REAL WAVARRAY(MAXSIZ),RBSIG(MAXSIZ)
 	COMMON /SPECTRA/WAVARRAY,DAT,RBSIG,NBIN,NLINES,CENTRE,OFLAG

c	CALCULATING THE EMISSION LINE PROFILES
c	(THEY ARE 2-GAUSSIANS AT THE MOMENT)
c
	line = 0.0
	linen = 0.0
	wid = 65.0
	widn = 25.0

        wavebarn = PAR(4)
        wavebar  = PAR(7)
        wmidn = wavebarn
        wmid = wavebar

c     extra settings to fit both OIII doublet lines at once
        IF (OFLAG.EQ.1) THEN
           wid = 55.0
           widn = 55.0
           doublet = 5006.8 - 4958.9
           wavebar2 = wavebar + doublet
           wavebarn2 = wavebarn + doublet
           wmidn = wavebarn + doublet/2.0
           wmid = wmidn
        END IF 
c     extra settings to fit both OII doublet lines at once
c     BUT we allow EW to vary in the second line
        IF (OFLAG.EQ.2) THEN
           wid = 55.0
           widn = 55.0
           doublet = 3728.8 - 3726.0
           wavebar2 = wavebar + doublet
           wavebarn2 = wavebarn + doublet
           wmidn = wavebarn + doublet/2.0
           wmid = wmidn
        END IF 

c  THIS IS THE NARROW GAUSSIAN - try to force positive EW:
        if(wave.ge.wmidn-widn.and.wave.lt.wmidn+widn)then
c           integn = flux*abs(PAR(2))
c removed abs - MJD 8/12/15
           integn = flux*(PAR(2))
           fwhmn = (wavebarn*abs(PAR(3)))/(300000.)
           sigman = fwhmn/(2.0*sqrt(2*log(2.0)))
           hn = integn/(sigman*sqrt(2.0*pi))
           linen = hn*exp(-0.5*((wave-wavebarn)/sigman)**2)
           if(OFLAG.EQ.1)
     +          linen=linen+3*hn*exp(-0.5*((wave-wavebarn2)/sigman)**2)
           if(OFLAG.EQ.2)then
c removed abs - MJD 8/12/15
c              integn = flux*abs(PAR(5))
              integn = flux*(PAR(5))
              fwhmn = (wavebarn*abs(PAR(3)))/(300000.)
              sigman = fwhmn/(2.0*sqrt(2*log(2.0)))
              hn = integn/(sigman*sqrt(2.0*pi))
              linen=linen+hn*exp(-0.5*((wave-wavebarn2)/sigman)**2)
           end if
        endif

c  THIS IS THE BROAD GAUSSIAN
        IF(NLINES.EQ.2.OR.NLINES.EQ.4) THEN
           if(wave.ge.wmid-wid.and.wave.lt.wmid+wid) then
              integ = flux*PAR(5)
              fwhm = (wavebar*PAR(6))/(300000.)
              sigma = fwhm/(2.0*sqrt(2*log(2.0)))
              h = integ/(sigma*sqrt(2.0*pi))
              line = h*exp(-0.5*((wave-wavebar)/sigma)**2)
              if(OFLAG.EQ.1)
     +             line=line+3*h*exp(-0.5*((wave-wavebar2)/sigma)**2)
           end if
        end if

c     ADD THE EMISSION LINES TO THE CONTINUUM
        flux = flux + line + linen
 
	return
	end

 
        subroutine mean_rms(x,nx,xmean,xrms)
* Simple routine to find mean and standard deviation of an
* array of values.
        implicit none
        real sumx,sumxx,rn,xmean,xrms,x(nx)
        integer nx,i

        if (nx.le.1) then
                stop 'mean_rms: error with N < 0'
                return
        end if

        rn = real(nx)
        sumx = 0.0
        sumxx = 0.0
        do i = 1, nx
                sumx = sumx + x(i)
                sumxx = sumxx + x(i)*x(i)
        end do
        xmean = sumx / rn
        xrms = sqrt(sumxx/rn - (sumx/rn)**2 )

        return
        end

      FUNCTION RTBIS(FUNC1,X1,X2,XACC)
      PARAMETER (JMAX=40)
      FMID=FUNC1(X2)
      F=FUNC1(X1)
C      IF(F*FMID.GE.0.) PAUSE 'Root must be bracketed for bisection.'
      IF(F*FMID.GE.0.) PRINT *,'Root must be bracketed for bisection.'
      IF(F.LT.0.)THEN
        RTBIS=X1
        DX=X2-X1
      ELSE
        RTBIS=X2
        DX=X1-X2
      ENDIF
      DO 11 J=1,JMAX
        DX=DX*.5
        XMID=RTBIS+DX
        FMID=FUNC1(XMID)
        IF(FMID.LE.0.)RTBIS=XMID
        IF(ABS(DX).LT.XACC .OR. FMID.EQ.0.) RETURN
11    CONTINUE
c      PAUSE 'too many bisections'
      PRINT *,'too many bisections'
      END


