
#include "ICEPLUME_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: ICEPLUME_TENDENCY_APPLY_T

C !INTERFACE:
      SUBROUTINE ICEPLUME_TENDENCY_APPLY_T(
     U                    gT_arr,
     I                    iMin,iMax,jMin,jMax,
     I                    k, bi, bj, myTime, myIter, myThid )

C     !DESCRIPTION:
C     As for S/R ICEPLUME_TENDENCY_APPLY_T
C     Add iceplume tendency terms to T (theta) tendency.
C     Routine works for one level at a time.

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "GRID.h"
#include "EEPARAMS.h"
#include "DYNVARS.h"
#include "ICEPLUME.h"


C     !INPUT/OUTPUT PARAMETERS:
C     gT_arr    :: the tendency array
C     k         :: Current vertical level index
C     bi,bj     :: Current tile indices
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number
C     myThid    :: my Thread Id number
      _RL     gT_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER iMin, iMax, jMin, jMax
      INTEGER k, bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid

C     !LOCAL VARIABLES:
      INTEGER i, j

#ifdef ALLOW_ICEPLUME

catn: this is the background melt tendency
      IF ( applyIcePlumeBGTendT ) THEN
       DO j=jMin,jMax
        DO i=iMin,iMax
          gT_arr(i,j) = gT_arr(i,j)
     &         + iceplumeBG_TendT(i,j,k,bi,bj)
        ENDDO
       ENDDO
      ENDIF

#endif /* ALLOW_ICEPLUME */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: ICEPLUME_TENDENCY_APPLY_S

C !INTERFACE:
      SUBROUTINE ICEPLUME_TENDENCY_APPLY_S(
     U                    gS_arr,
     I                    iMin,iMax,jMin,jMax,
     I                    k, bi, bj, myTime, myIter, myThid )

C     !DESCRIPTION:
C     As for S/R ICEPLUME_TENDENCY_APPLY_S
C     Add iceplume tendency terms to S tendency.
C     Routine works for one level at a time.

C     !INPUT PARAMETERS:
      IMPLICIT NONE
#include "SIZE.h"
#include "GRID.h"
#include "EEPARAMS.h"
#include "DYNVARS.h"
#include "ICEPLUME.h"

C     !INPUT/OUTPUT PARAMETERS:
C     gS_arr    :: the tendency array
C     k         :: Current vertical level index
C     bi,bj     :: Current tile indices
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number
C     myThid    :: my Thread Id number
      _RL     gS_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER iMin, iMax, jMin, jMax
      INTEGER k, bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES:
      INTEGER i, j

#ifdef ALLOW_ICEPLUME

catn: this is the background melt tendency
      IF ( applyIcePlumeBGTendS ) THEN
       DO j=jMin,jMax
        DO i=iMin,iMax
          gS_arr(i,j) = gS_arr(i,j)
     &         + iceplumeBG_TendS(i,j,k,bi,bj)
        ENDDO
       ENDDO
      ENDIF

#endif /* ALLOW_ICEPLUME */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0

#ifdef ALLOW_PTRACERS
      SUBROUTINE ICEPLUME_TENDENCY_APPLY_PTRACERS(
     U                    gPtracer,
     I                    k, bi, bj, iTracer, myTime, myIter, myThid )

      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "ICEPLUME.h"
#include "ICEPLUME_OPTIONS.h"
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"

      _RL     gPtracer(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER k, bi, bj, iTracer
      _RL     myTime
      INTEGER myIter
      INTEGER myThid

      INTEGER i, j

      DO j=1,sNy
       DO i=1,sNx
        IF ( hFacC(i,j,k,bi,bj) .NE. 0. _d 0 ) THEN
         gPtracer(i,j) = gPtracer(i,j)
     &    + addMass3Dplume(i,j,k,bi,bj) * mass2rUnit
     &    * ( ptr_addMass3D(i,j,k,bi,bj,iTracer)
     &        - pTracer(i,j,k,bi,bj,iTracer) )
     &    * recip_rA(i,j,bi,bj)
     &    * recip_drF(k)
     &    * _recip_hFacC(i,j,k,bi,bj)
        ENDIF
       ENDDO
      ENDDO

      RETURN
      END
#endif
