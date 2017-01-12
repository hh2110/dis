!-----------------------------------------------------------------------
! SUBROUTINE TO DEFINE THE GAMMA SURFACES
!-----------------------------------------------------------------------
      SUBROUTINE GSF
      USE VAR 
      IMPLICIT NONE
! DEFINE VARIABLES
      REAL :: NORM     !! NORMALISATION FACTOR
      INTEGER :: I     !! GSF COEFFICIENT INDEX
!
!$OMP PARALLEL DO PRIVATE(PHID,I)        
      DO PHID=1,NPHIS      
!-----------------------------------------------------------------------
! GAMMA PHASE
!-----------------------------------------------------------------------
      IF(GSURF(PHID) == 10) THEN !! SINGLE LAYER LI-SHEN GSF
        GSF_EQ(PHID)=10
        GC(PHID,1)  = 376.2336646
        GC(PHID,2)  = 0.0 
        GC(PHID,3)  = 0.0 
        GC(PHID,4)  =-107.8224430 
        GC(PHID,5)  = 0.0 
        GC(PHID,6)  = 0.0
        GC(PHID,7)  =-15.30011049
        GC(PHID,8)  =-4.4112215 
        GC(PHID,9)  = 1.061276727 
        GC(PHID,10) = 0.0 
        GC(PHID,11) =-102.1717455
        GC(PHID,12) = 0.0
        GC(PHID,13) = 0.0
        GC(PHID,14) = 35.86964829
        GC(PHID,15) = 8.68759370       
      ELSE IF(GSURF(PHID) == 11) THEN !! SINGLE LAYER VOSKOBOINIKOV-VORONTSOV GSF
        GSF_EQ(PHID)=10
        GC(PHID,1)  = 709.4070619583151
        GC(PHID,2)  = 0.0
        GC(PHID,3)  = 0.0
        GC(PHID,4)  =-212.06195869445665
        GC(PHID,5)  = 0.0
        GC(PHID,6)  = 0.0
        GC(PHID,7)  =-35.60558246205095
        GC(PHID,8)  = 8.530979347228337
        GC(PHID,9)  = 1.333770578253783
        GC(PHID,10) = 0.0
        GC(PHID,11) =-320.00000000000017
        GC(PHID,12) = 0.0
        GC(PHID,13) = 0.0
        GC(PHID,14) =-21.461829674214595
        GC(PHID,15) = 0.7658514012257662        
      ELSE IF(GSURF(PHID) == 12) THEN !! DOUBLE-LAYER VOSKOBOINIKOV-VORONTSOV GSF (solver-4x)
        GSF_EQ(PHID)=20
        GC(PHID,1)  = 325.4444444444444
        GC(PHID,2)  = 0.0
        GC(PHID,3)  = 0.0
        GC(PHID,4)  =-11.833333333333348
        GC(PHID,5)  =-77.70370370370368
        GC(PHID,6)  =-0.9259259259259454
        GC(PHID,7)  = 0.9259259259259487
        GC(PHID,8)  =-18.944444444444432
        GC(PHID,9)  = 0.0
        GC(PHID,10) = 0.0
        GC(PHID,11) =-37.81644263192049
        GC(PHID,12) = 0.0
        GC(PHID,13) =-37.23909236273086
        GC(PHID,14) = 0.0
        GC(PHID,15) = 0.0          
!-----------------------------------------------------------------------
! GAMMA PRIME PHASE
!-----------------------------------------------------------------------        
      ELSE IF(GSURF(PHID) == 20) THEN !! SINGLE LAYER SCHOECK-SHEN GSF
        GSF_EQ(PHID)=11
        GC(PHID,1)  = 728.73356
        GC(PHID,2)  =-127.927598
        GC(PHID,3)  = 60.5122581
        GC(PHID,4)  =-164.700179
        GC(PHID,5)  = 7.23846759
        GC(PHID,6)  =-3.15972992
        GC(PHID,7)  = 0.0
        GC(PHID,8)  =-22.1128735
        GC(PHID,9)  = 0.0
        GC(PHID,10) = 207.178649
        GC(PHID,11) =-284.966972
        GC(PHID,12) = 2.23832338 
        GC(PHID,13) = 6.59233717 
        GC(PHID,14) = 18.8309473
        GC(PHID,15) = 0.0          
      ELSE IF(GSURF(PHID) == 21) THEN !! SINGLE LAYER VOSKOBOINIKOV-VORONTSOV GSF
        GSF_EQ(PHID)=12
        GC(PHID,1)  = 851.0672863405989
        GC(PHID,2)  =-51.961524227066356
        GC(PHID,3)  =-10.71152422706633
        GC(PHID,4)  =-193.60257121980004
        GC(PHID,5)  =-0.26923788646691355
        GC(PHID,6)  =-4.288475772933729
        GC(PHID,7)  =-32.0224287801998
        GC(PHID,8)  = 9.435904553133536
        GC(PHID,9)  = 0.0
        GC(PHID,10) = 70.00000000000003
        GC(PHID,11) =-342.44273848271786
        GC(PHID,12) =-1.0806716851095137
        GC(PHID,13) =-4.999999999999988
        GC(PHID,14) =-19.60771296085209
        GC(PHID,15) = 0.0          
      ELSE IF(GSURF(PHID) == 22) THEN !! DOUBLE-LAYER VOSKOBOINIKOV 1 (solver-4x)
        GSF_EQ(PHID)=21
        GC(PHID,1)  = 433.4444444444445
        GC(PHID,2)  =-22.27307312721769
        GC(PHID,3)  =-53.00000000000001
        GC(PHID,4)  = 31.166666666666675
        GC(PHID,5)  =-76.7037037037037
        GC(PHID,6)  = 5.657407407407426
        GC(PHID,7)  =-5.324074074074083
        GC(PHID,8)  =-29.444444444444443
        GC(PHID,9)  = 5.4397397938843515
        GC(PHID,10) =-0.5534180126147842
        GC(PHID,11) =-55.71430097679889
        GC(PHID,12) = 4.416666666666689
        GC(PHID,13) =-47.63139720814412
        GC(PHID,14) = 5.220084679281477
        GC(PHID,15) = 0.0          
      ELSE IF(GSURF(PHID) == 23) THEN !! SINGLE LAYER HH - 100mJ SISF 
        GSF_EQ(PHID)=12
        GC(PHID,1)  =848.983953007 
        GC(PHID,2)  =-51.9615242271
        GC(PHID,3)  =-6.5448575604
        GC(PHID,4)  =-195.685904553
        GC(PHID,5)  =-2.3525712198
        GC(PHID,6)  =-4.28847577293
        GC(PHID,7)  =-29.9390954469
        GC(PHID,8)  =10.1303489976
        GC(PHID,9)  = 0.0
        GC(PHID,10) =70.0
        GC(PHID,11) =-338.8342993
        GC(PHID,12) =-4.68911086754
        GC(PHID,13) =-5.0
        GC(PHID,14) =-18.4048999
        GC(PHID,15) =0.0
      ELSE IF(GSURF(PHID) == 24) THEN !! SINGLE LAYER HH - 50mJ SISF
        GSF_EQ(PHID)=12
        GC(PHID,1)  =851.067286341
        GC(PHID,2)  =-51.9615242271
        GC(PHID,3)  =-10.7115242271
        GC(PHID,4)  =-193.60257122
        GC(PHID,5)  =-0.269237886467
        GC(PHID,6)  =-4.28847577293
        GC(PHID,7)  =-32.0224287802
        GC(PHID,8)  =9.43590455313
        GC(PHID,9)  =0.0
        GC(PHID,10) =70.0
        GC(PHID,11) =-342.442738483
        GC(PHID,12) =-1.08067168511
        GC(PHID,13) =-5.0
        GC(PHID,14) =-19.6077129609
        GC(PHID,15) = 0.0
      ELSE IF(GSURF(PHID) == 25) THEN !! SINGLE LAYER HH - 10mJ SISF
        GSF_EQ(PHID)=12
        GC(PHID,1)  =852.733953007
        GC(PHID,2)  =-51.9615242271
        GC(PHID,3)  =-14.0448575604
        GC(PHID,4)  =-191.935904553
        GC(PHID,5)  =1.3974287802
        GC(PHID,6)  =-4.28847577293
        GC(PHID,7)  =-33.6890954469
        GC(PHID,8)  =8.88034899758
        GC(PHID,9)  =0.0
        GC(PHID,10) =70.0
        GC(PHID,11) =-345.329489829
        GC(PHID,12) =1.80607966084
        GC(PHID,13) =-5.0
        GC(PHID,14) =-20.5699634095
        GC(PHID,15) =0.0
      ELSE
        STOP '>>ERROR - INVALID GSF FOR GAMMA PRIME PHASE!'
      END IF
! PRINT GSF COEFFICIENTS TO SCREEN
        PRINT *, ''
        PRINT *, '>>NOTE - GSF FOURIER COEFFICIENTS'
        PRINT '(" PHASE = PHI(",I1,")")',PHID    
        DO I=1,15
          PRINT '(" GC(",I1,I2,")=",F8.2)',PHID,I,(GC(PHID,I))
        END DO   
! NORMALISE FOURIER COEFFICEINTS TO GET GSF IN REDUCED UNITS
        NORM =AA/SQRT(3.0)*C44X/10.0*1.5*1000.0
        DO I=1,15
          GC(PHID,I) =GC(PHID,I)/NORM
        END DO       
!        
      END DO
!$OMP END PARALLEL DO             
!      
      END SUBROUTINE GSF
