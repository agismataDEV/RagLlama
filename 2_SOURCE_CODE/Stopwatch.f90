MODULE StopWatch

! =============================================================================
!
!   Programmer: Jasper de Koning
!
!   Language: Fortran 90
!
!   Version Date    Comment
!   ------- -----   -------
!   1.0     060222  Original code (JdK)
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The module Stopwatch contains subroutines and functions that help to 
!   keep track of the amount of time that separate parts of the program use.
!
! *****************************************************************************
!
!   MODULES USED
!
	USE Types
!#ifndef RUNESTIMATOR
!    USE ParnacMPIDef
!#endif

! *****************************************************************************
!
!   GLOBAL DECLARATIONS
!
!   cStopWatch  type   Contains all information on a specific Stopwatch
!       iStartTime  i8b    start time in 100s of seconds of a Stopwatch run
!       iTotalTime  i8b    total time in 100s of seconds counted since start 
!                          of program
!       iLastTime   i8b    time in 100s of seconds elapsed since iStartTime at
!                          the stopping of a Stopwatch run
!       iCallCount  i8b    number of calls to the routine SWStartAndCount since
!                          since start of program.
!       bRunning    int    whether the Stopwatch is running (1) or not (0)
!       iTimeF      i8b    unit number of the Stopwatch log file
!       pszText     char   identification text in the Stopwatch log file
!
	IMPLICIT NONE

	SAVE

	type cStopWatch
		integer(i8b)::			iStartTime, iTotalTime, iLastTime;
		integer(i8b)::			iCallCount;
		integer::				bRunning;
		integer(i8b)::			iTimeF;
		character(LEN=80)::		pszText;
	end type cStopWatch

! *****************************************************************************
!
!   CONTAINED SUBROUTINES AND FUNCTIONS
!
!   SWInit           sub   Initializes the Stopwatch
!   SWStart          sub   Starts a Stopwatch timer session
!   SWStartAndCount  sub   Same as SWStart, but also increases iCallCount
!   SWStop           sub   Stops a timer session and adds runtime to iTotalTime
!   SWStopAndWrite   sub   Same as SWStop, but also echoes runtime to SW log
!   SWReset          sub   Resets the Stopwatch
!   SWOutputToFile   sub   Echoes total runtime to SW log
!   GetIntTime       fun   Translates (hrs,mins,secs,100s) format to 100s    
!   GetIntTimeNOW    fun   Obtain the current time from the system
!   TimeElapsed      fun   Calculates the time elapsed since iStartTime
!   OpenTimeFile     sub   Opens the SW log file
!   WriteToTimeFile  sub   Writes text to SW log file
!   CloseTimeFile    sub   Closes the SW log file
!
! =============================================================================
!
	CONTAINS

	SUBROUTINE SWInit(cswIn, iUnit, pszText)

! =============================================================================
!
!   Programmer: Jasper de Koning
!
!   Language: Fortran 90
!
!   Version Date    Comment
!   ------- -----   -------
!   1.0     060222  Original code (JdK)
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The subroutine SWInit initializes a specific stopwatch
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cswIn   io  type(cStopWatch) Stopwatch to be initialized
!   iUnit   i   int              Unit number of SW log file
!   pszText i   char             Identification text for this Stopwatch
!   si  o   dp    value of Si(x)
!
		type(cStopWatch), INTENT(inout)::	cswIn;
		INTEGER, INTENT(IN)::			iUnit;
		character(LEN=*), INTENT(IN)::		pszText;

! *****************************************************************************
! 
!   LOCAL PARAMETERS      
!
!   none   
!
! *****************************************************************************
!
!   I/O
!
!   none
!   
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   none
!   
! =============================================================================

		cswIn%iTimeF		= iUnit;	
		cswIn%pszText		= pszText;
		cswIn%iLastTime		= 0;
		cswIn%iTotalTime	= 0;
		cswIn%iCallCount	= 0;

	END SUBROUTINE SWInit
    
    SUBROUTINE SWStart(swIn)

! =============================================================================
!
!   Programmer: Jasper de Koning
!
!   Language: Fortran 90
!
!   Version Date    Comment
!   ------- -----   -------
!   1.0     060222  Original code (JdK)
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The subroutine SWStart starts a Stopwatch timing session
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   swIn   io  type(cStopWatch) Stopwatch to be started
!
		type(cStopWatch), INTENT(inout)::		swIn;

! *****************************************************************************
! 
!   LOCAL PARAMETERS      
!
!   none   
!
! *****************************************************************************
!
!   I/O
!
!   none
!   
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   GetIntTimeNow
!   
! =============================================================================

		swIn%bRunning		= 1;
		swIn%iStartTime		= GetIntTimeNow();

	END SUBROUTINE SWStart

	SUBROUTINE SWStartAndCount(swIn)

! =============================================================================
!
!   Programmer: Jasper de Koning
!
!   Language: Fortran 90
!
!   Version Date    Comment
!   ------- -----   -------
!   1.0     060222  Original code (JdK)
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The subroutine SWStartAndCount starts a Stopwatch timing session, 
!   and it also increases the iCallCount counter - may be used to 
!   count the number of calls of a specific subroutine, function, loop...
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   swIn   io  type(cStopWatch) Stopwatch to be started
!
		type(cStopWatch), INTENT(inout)::		swIn;

! *****************************************************************************
! 
!   LOCAL PARAMETERS      
!
!   none   
!
! *****************************************************************************
!
!   I/O
!
!   none
!   
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   GetIntTimeNow
!   
! =============================================================================

		swIn%bRunning		= 1;
		swIn%iStartTime		= GetIntTimeNow();
		swIn%iCallCount		= swIn%iCallCount + 1;

	END SUBROUTINE SWStartAndCount

    SUBROUTINE SWStop(swIn)

! =============================================================================
!
!   Programmer: Jasper de Koning
!
!   Language: Fortran 90
!
!   Version Date    Comment
!   ------- -----   -------
!   1.0     060222  Original code (JdK)
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The subroutine SWStop ends a Stopwatch timing session. iLastTime
!   now contains the time elapsed since the last start time, and iTotalTime 
!   is increased with the elapsed time.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   swIn   io  type(cStopWatch) Stopwatch to be stopped
!
		type(cStopWatch), INTENT(inout)::		swIn;

! *****************************************************************************
! 
!   LOCAL PARAMETERS      
!
!   none   
!
! *****************************************************************************
!
!   I/O
!
!   none
!   
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   TimeElapsed
!   GetIntTimeNow
!   
! =============================================================================

		swIn%bRunning		= 0;
		swIn%iLastTime		= TimeElapsed(swIn%iStartTime);
		swIn%iTotalTime		= swIn%iTotalTime + swIn%iLastTime;

		swIn%iStartTime		= GetIntTimeNow();

	END SUBROUTINE SWStop
    
	SUBROUTINE SWStopAndWrite(cswIn)

! =============================================================================
!
!   Programmer: Jasper de Koning
!
!   Language: Fortran 90
!
!   Version Date    Comment
!   ------- -----   -------
!   1.0     060222  Original code (JdK)
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The subroutine SWStopAndWrite ends a Stopwatch timing session. iLastTime
!   now contains the time elapsed since the last start time, and iTotalTime 
!   is increased with the elapsed time. The elapsed time is output to the SW
!   log file.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   swIn   io  type(cStopWatch) Stopwatch to be stopped
!
		type(cStopWatch), INTENT(inout)::		cswIn;

! *****************************************************************************
! 
!   LOCAL PARAMETERS      
!
!   none   
!
! *****************************************************************************
!
!   I/O
!
!   Current time to the SW log file
!   
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   SWStop
!   
! =============================================================================

		call SWStop(cswIn);

		write (UNIT = cswIn%iTimeF, FMT='(A40, I8, I10)') TRIM(cswIn%pszText), cswIn%iCallCount, cswIn%iLastTime;

	END SUBROUTINE SWStopAndWrite

	SUBROUTINE SWReset(swIn)

! =============================================================================
!
!   Programmer: Jasper de Koning
!
!   Language: Fortran 90
!
!   Version Date    Comment
!   ------- -----   -------
!   1.0     060222  Original code (JdK)
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The subroutine SWReset resets the Stopwatch to zero time and count values.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   swIn   io  type(cStopWatch) Stopwatch to be reset
!
		type(cStopWatch), INTENT(inout)::		swIn;

! *****************************************************************************
! 
!   LOCAL PARAMETERS      
!
!   none   
!
! *****************************************************************************
!
!   I/O
!
!   none
!   
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   none
!   
! =============================================================================

		swIn%bRunning		= 0;
		swIn%iStartTime		= 0;
		swIn%iTotalTime		= 0;
		swIn%iCallCount		= 0;

	END SUBROUTINE SWReset

	SUBROUTINE SWOutputToFile(cswIn)

! =============================================================================
!
!   Programmer: Jasper de Koning
!
!   Language: Fortran 90
!
!   Version Date    Comment
!   ------- -----   -------
!   1.0     060222  Original code (JdK)
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The subroutine SWOutputToFile echoes the total run time of this stopwatch
!   to the logfile.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cswIn   io  type(cStopWatch) Stopwatch to be echoed
!
		type(cStopWatch), INTENT(IN)::			cswIn;

! *****************************************************************************
! 
!   LOCAL PARAMETERS      
!
!   none   
!
! *****************************************************************************
!
!   I/O
!
!   Total time to the SW log file
!   
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   none
!   
! =============================================================================

		write (UNIT = cswIn%iTimeF, FMT='(A40, I8, I10)') TRIM(cswIn%pszText), cswIn%iCallCount, cswIn%iTotalTime;

	END SUBROUTINE SWOutputToFile

	INTEGER(i8b) FUNCTION GetIntTime(iHr, iMin, iSec, iHunSec)

! =============================================================================
!
!   Programmer: Jasper de Koning
!
!   Language: Fortran 90
!
!   Version Date    Comment
!   ------- -----   -------
!   1.0     060222  Original code (JdK)
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The function GetIntTime translates the (Hr, Min, Sec, HunSec) format
!   of a system time to a single number containing the total 100s of secs.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   GetIntTime  o   i8b  total time in 100s of seconds
!   iHr         i   i8b  part of the time containing the hours
!   iMin        i   i8b  part of the time containing the minutes
!   iSec        i   i8b  part of the time containing the seconds
!   iHunSec     i   i8b  part of the time containing the 100s of seconds
!
		INTEGER(i8b), INTENT(IN)::		iHr, iMin, iSec, iHunSec

! *****************************************************************************
! 
!   LOCAL PARAMETERS      
!
!   none   
!
! *****************************************************************************
!
!   I/O
!
!   none
!   
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   none
!   
! =============================================================================

		GetIntTime	= ((iHr * 60_i8b + iMin) * 60_i8b + iSec) * 100_i8b + iHunSec;

	END FUNCTION GetIntTime

	integer(i8b) FUNCTION GetIntTimeNOW()

! =============================================================================
!
!   Programmer: Jasper de Koning
!
!   Language: Fortran 90
!
!   Version Date    Comment
!   ------- -----   -------
!   1.0     060222  Original code (JdK)
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The function GetIntTimeNOW returns the current system time in 100s of 
!   seconds.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   GetIntTimeNOW   o   i8b   current system time in 100s of seconds
!
! *****************************************************************************
! 
!   LOCAL PARAMETERS      
!
!   iHr      i   i8b  part of the time containing the hours
!   iMin     i   i8b  part of the time containing the minutes
!   iSec     i   i8b  part of the time containing the seconds
!   iHunSec  i   i8b  part of the time containing the 100s of seconds
!
		integer(i8b)::				iHr, iMin, iSec, iHunSec;
!
! *****************************************************************************
!
!   I/O
!
!   none
!   
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   GetIntTimeNow
!   GETTIM or MPI_WTIME, depending on the system
!   
! =============================================================================

!		call GETTIM(iHr, iMin, iSec, iHunSec);
!       GetIntTimeNOW = GetIntTime(iHr, iMin, iSec, iHunSec)
#ifndef RUNESTIMATOR
                GetIntTimeNOW = 100*MPI_Wtime();		! This is multiplied by 100, in order to measure 1/100ths of secs, rather than whole seconds
!#else
!                GetIntTimeNOW  = 0                      ! Function disabled for RUNESTIMATOR
#endif
	END FUNCTION GetIntTimeNOW

    integer(i8b) FUNCTION TimeElapsed(iStartTime)

! =============================================================================
!
!   Programmer: Jasper de Koning
!
!   Language: Fortran 90
!
!   Version Date    Comment
!   ------- -----   -------
!   1.0     060222  Original code (JdK)
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The function TimeElapsed returns the current system time in 100s of 
!   seconds relative to a certain start time.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   TimeElapsed     o   i8b   Elapsed time since the given Start time
!   iStartTime      i   i8b   Start time from which the elapsed time is given.
!
		integer(i8b), INTENT(in)::		iStartTime
!
! *****************************************************************************
! 
!   LOCAL PARAMETERS      
!
!   none
!
! *****************************************************************************
!
!   I/O
!
!   none
!   
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   GetIntTimeNow
!   
! =============================================================================

		TimeElapsed = GetIntTimeNow() - iStartTime;

	END FUNCTION TimeElapsed
    
	SUBROUTINE OpenTimeFile(iUnit, sFileName)

! =============================================================================
!
!   Programmer: Jasper de Koning
!
!   Language: Fortran 90
!
!   Version Date    Comment
!   ------- -----   -------
!   1.0     060222  Original code (JdK)
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The subroutine OpenTimeFile opens the Stopwatch log file. It is rather 
!   simple, no error checking whatsoever.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   iUnit      i   int   Unit number to be assigned to the SW logfile
!   sFileName  i   char  Name of the SW log file
!
		INTEGER, INTENT(IN)::			iUnit;
		CHARACTER(LEN=*), INTENT(IN)::		sFileName;
!
! *****************************************************************************
! 
!   LOCAL PARAMETERS      
!
!   none
!
! *****************************************************************************
!
!   I/O
!
!   Opens the Stopwatch log file
!   
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   none
!   
! =============================================================================

		OPEN(UNIT = iUnit, FILE = sFileName, STATUS="REPLACE",form="FORMATTED");

	END SUBROUTINE OpenTimeFile

	SUBROUTINE WriteToTimeFile(iUnit, pszText)

! =============================================================================
!
!   Programmer: Jasper de Koning
!
!   Language: Fortran 90
!
!   Version Date    Comment
!   ------- -----   -------
!   1.0     060222  Original code (JdK)
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The subroutine WriteToTimeFile writes text to the Stopwatch log file.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   iUnit    i   int   Unit number of the SW log file
!   pszText  i   char  Text to be output to the SW log file
!
		integer, INTENT(IN)::			iUnit;
		character(LEN=*), INTENT(IN)::		pszText;
!
! *****************************************************************************
! 
!   LOCAL PARAMETERS      
!
!   none
!
! *****************************************************************************
!
!   I/O
!
!   Writes text to the Stopwatch log file
!   
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   none
!   
! =============================================================================

		write (UNIT = iUnit, FMT=*) pszText

	END SUBROUTINE WriteToTimeFile

	SUBROUTINE CloseTimeFile(iUnit)

! =============================================================================
!
!   Programmer: Jasper de Koning
!
!   Language: Fortran 90
!
!   Version Date    Comment
!   ------- -----   -------
!   1.0     060222  Original code (JdK)
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The subroutine CloseTimeFile closes the Stopwatch log file.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   iUnit    i   int   Unit number of the SW log file
!
		INTEGER, INTENT(IN)::			iUnit;
!
! *****************************************************************************
! 
!   LOCAL PARAMETERS      
!
!   none
!
! *****************************************************************************
!
!   I/O
!
!   Closes the Stopwatch log file
!   
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   none
!   
! =============================================================================
		
		close(iUnit);

	END SUBROUTINE
	
END MODULE STOPWATCH
