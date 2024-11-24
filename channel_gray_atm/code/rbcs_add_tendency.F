#include "RBCS_OPTIONS.h"

CBOP
C !ROUTINE: RBCS_ADD_TENDENCY

C !INTERFACE: ==========================================================
      SUBROUTINE RBCS_ADD_TENDENCY(
     U                    gTendency,
     I                    k, bi, bj, tracerNum,
     I                    myTime, myIter, myThid )

C !DESCRIPTION:
C     Add to tendency array the contribution from 3-D field relaxation

C !USES: ===============================================================
      use  simple_Sat_Vapor_Pres_Mod, ONLY: ESCOMP
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#ifdef ALLOW_PTRACERS
#include "PTRACERS_SIZE.h"
#include "PTRACERS_FIELDS.h"
#endif
#include "RBCS_SIZE.h"
#include "RBCS_PARAMS.h"
#include "RBCS_FIELDS.h"

C !INPUT/OUTPUT PARAMETERS: ============================================
C  gTendency      :: the tendency array
C  k              :: vertical level index
C  bi,bj          :: tile indices
C  tracerNum      :: tracer number (1=Temp, 2=Salt, >2 : ptracer)
C  myTime         :: current time
C  myIter         :: current timestep
C  myThid         :: my Thread Id number
      _RL gTendency(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER k, bi, bj
      INTEGER tracerNum
      _RL myTime
      INTEGER myIter
      INTEGER myThid

C !LOCAL VARIABLES: ====================================================
C  i,j            :: loop indices
      INTEGER i,j
      INTEGER irbc
      _RL rbcsVanishingFac
      _RL rec_tauRlx
#ifdef ALLOW_PTRACERS
      INTEGER iTracer
#endif
#ifdef ALLOW_ATM_PHYS
      _RL conv_theta2T
      _RL temp(0:sNx+1,0:sNy+1)
      _RL qSat(0:sNx+1,0:sNy+1)
#endif
CEOP

#ifdef ALLOW_RBCS

      IF ( rbcsVanishingTime.GT.0. _d 0 ) THEN
        rbcsVanishingFac =
     &      MAX( 0. _d 0 , 1. _d 0 - myTime / rbcsVanishingTime )
      ELSE
        rbcsVanishingFac = 1. _d 0
      ENDIF

#ifndef DISABLE_RBCS_MOM
      IF ( tracerNum.EQ.-1 .AND. useRBCuVel ) THEN
        rec_tauRlx = rbcsVanishingFac/tauRelaxU
        DO j=0,sNy+1
         DO i=0,sNx+1
          gTendency(i,j) = gTendency(i,j)
     &       - RBC_maskU(i,j,k,bi,bj)*rec_tauRlx
     &        *( uVel(i,j,k,bi,bj)- RBCuVel(i,j,k,bi,bj) )
         ENDDO
        ENDDO
      ENDIF
      IF ( tracerNum.EQ.-2 .AND. useRBCvVel ) THEN
        rec_tauRlx = rbcsVanishingFac/tauRelaxV
        DO j=0,sNy+1
         DO i=0,sNx+1
          gTendency(i,j) = gTendency(i,j)
     &       - RBC_maskV(i,j,k,bi,bj)*rec_tauRlx
     &        *( vVel(i,j,k,bi,bj)- RBCvVel(i,j,k,bi,bj) )
         ENDDO
        ENDDO
      ENDIF
#endif /* DISABLE_RBCS_MOM */

      IF ( tracerNum.EQ.1 .AND. useRBCtemp ) THEN
        irbc = MIN(maskLEN,tracerNum)
        rec_tauRlx = rbcsVanishingFac/tauRelaxT
        DO j=0,sNy+1
         DO i=0,sNx+1
          gTendency(i,j) = gTendency(i,j)
     &       - RBC_mask(i,j,k,bi,bj,irbc)*rec_tauRlx
     &        *( theta(i,j,k,bi,bj)- RBCtemp(i,j,k,bi,bj) )

         ENDDO
        ENDDO
      ENDIF

      IF ( tracerNum.EQ.2 .AND. useRBCsalt ) THEN
        irbc = MIN(maskLEN,tracerNum)
        rec_tauRlx = rbcsVanishingFac/tauRelaxS
#ifdef ALLOW_ATM_PHYS
C-- addition starts
       IF ( useAtm_Phys ) THEN
        conv_theta2T = (rC(k)/atm_po)**atm_kappa
        DO j=0,sNy+1
         DO i=0,sNx+1
c         Temp(i,j) = theta(i,j,k,bi,bj)*conv_theta2T
          Temp(i,j) = ( tRef(k)
     &      + maskC(i,j,k,bi,bj)*( theta(i,j,k,bi,bj) - tRef(k) )
     &                )*conv_theta2T
         ENDDO
        ENDDO
        CALL ESCOMP( Temp, qSat )
        DO j=0,sNy+1
         DO i=0,sNx+1
           RBCsalt(i,j,k,bi,bj) = 0.8 _d 0*qSat(i,j)
         ENDDO
        ENDDO
       ENDIF
C-- addition ends
#endif
        DO j=0,sNy+1
         DO i=0,sNx+1
          gTendency(i,j) = gTendency(i,j)
     &       - RBC_mask(i,j,k,bi,bj,irbc)*rec_tauRlx
     &        *( salt(i,j,k,bi,bj)- RBCsalt(i,j,k,bi,bj) )
         ENDDO
        ENDDO
      ENDIF

#ifdef ALLOW_PTRACERS
      IF ( usePTRACERS .AND. tracerNum.GT.2 ) THEN
       iTracer = tracerNum-2
       IF ( useRBCpTrNum(iTracer) ) THEN
        irbc = MIN(maskLEN,tracerNum)
        rec_tauRlx = rbcsVanishingFac/tauRelaxPTR(iTracer)
        DO j=0,sNy+1
         DO i=0,sNx+1
          gTendency(i,j) = gTendency(i,j)
     &       - RBC_mask(i,j,k,bi,bj,irbc)*rec_tauRlx
     &        *( pTracer(i,j,k,bi,bj,iTracer)
     &           - RBC_ptracers(i,j,k,bi,bj,iTracer) )
         ENDDO
        ENDDO
       ENDIF
      ENDIF
#endif /* ALLOW_PTRACERS */

#endif /* ALLOW_RBCS */

      RETURN
      END
