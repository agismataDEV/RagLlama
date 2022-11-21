MODULE ParnacGeneral

! =============================================================================
!
!   Programmer: Jasper de Koning
!
!   Language: Fortran 90
!
!   Version Date    Comment
!   ------- -----   -------
!   1.0     090505  Original code (JdK)
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The module ParnacGeneral contains global declarations, subroutines and
!   functions that are used on many occasions within the main program and
!   modules of the Parnacprogram..
!
! *****************************************************************************
!
!   MODULES USED
!
    USE Types
    USE StopWatch

! *****************************************************************************
!
!   GLOBAL DECLARATIONS
!
!   sSystemSlash    char   either '/' or '\', depending on the system being a
!                          Unix or a Windows system
!   sInputDir       char   Input directory for all config and input files
!   sOutputDir      char   Output directory for all output and log files
!   sTempDir        char   Temporary directory for all temporary storage
!   sDefInputDir    char   Default input directory
!   sDefOutputDir   char   Default output directory
!   sDefTempDir     char   Default temporary directory
!   sInputName      char   Name of the input configuration file
!   sOutputRoot     char   Root name of all output files
!   iImportUNIT     i4b    Unit number for all import activities
!   iExportUNIT     i4b    Unit number for all export activities
!   iLogFileUNIT    i4b    Unit number for the general log file
!   iLogSWFileUNIT  i4b    Unit number for the Stopwatch log file
!   iDebugLvl       i8b    Overall debug output message level
!   iTimeStamp      i8b    Time stamp at the beginning of the program
!   iMemAllocated   i8b    Amount of memory used - counting only large arrays
!   cswXXXX         type(cStopwatch)   numerous definitions of various
!                          stopwatches that keep track of different parts of
!                          the program
!
    IMPLICIT NONE
    SAVE

    !*************************************************************************************
    !
    !        Global file and directory variables
    !
    !*************************************************************************************

    character:: sSystemSlash = '/' ! Is this unix (/) or windows (\)?
    ! windows also understands '/'

    character(LEN=128) ::         sInputDir = ''
    character(LEN=128) ::         sOutputDir = ''
    character(LEN=128) ::         sTempDir = ''
    character(LEN=128) ::         sDefInputDir = '.'
    character(LEN=128) ::         sDefOutputDir = '.'
    character(LEN=128) ::         sDefTempDir = '.'

    character(LEN=128) ::         sInputName = ''
    character(LEN=256)::        sOutputRoot = ''

    integer(i4b), parameter ::  iImportUNIT = 42
    integer(i4b), parameter ::  iExportUNIT = 43

    !****************************************************************************************
    !
    !        Global variables for the log and stopwatch files
    !
    !****************************************************************************************

    ! debug parameter
    !  0: The biggest steps in the algorithm
    !  1: The start and end of steps like: Compute the convolution between two grids
    !  2: Steps inside of functions of level 1
    !  3: Large loops / functions inside of level 2 functions
    !  4: etc.
    integer(i8b) :: iDebugLvl

    integer(i4b) :: iLogFileUNIT = 30
    integer(i4b), parameter :: iLogSWFileUNIT = 36        ! The output unit for the stopwatches

    ! For debugging
    integer(i8b) :: iTimeStamp
    integer(i8b) ::        iMemAllocated

    ! ********************************************************************************************************
    ! Variables to compute the amount of time spent on each specific part of the code. This is done using a
    ! a (class) Stopwatch. See the class itself for more information
    ! ********************************************************************************************************

    ! First we write down the variables, which measure parts of the iteration scheme
    type(cStopwatch)::                                        cswProgram
    type(cStopwatch)::                                        cswIniti
    type(cStopwatch)::                                        cswLinStep
    type(cStopwatch)::                                        cswCorrStep
    type(cStopwatch)::                                        cswNonLinContr

    ! Secondly, we write down the variables, which measure parts of the code
    ! Green's function
    type(cStopwatch)::                                        cswGreen
    ! The transforms
    type(cStopwatch)::                                        cswTrans
    type(cStopwatch)::                                        cswTransXY
    type(cStopwatch)::                                        cswTransXYInv
    type(cStopwatch)::                                        cswTransXYZ
    type(cStopwatch)::                                        cswTransXYZInv
    type(cStopwatch)::                                        cswTransGreenXY
    type(cStopwatch)::                                        cswTransGreenXYZ
    type(cStopwatch)::                                        cswTransT
    type(cStopwatch)::                                        cswTransTInv
    ! Differentiation
    type(cStopwatch)::                                        cswDeriv
    ! Block-related
    type(cStopwatch)::                                        cswBlock
    type(cStopwatch)::                                        cswBlockRedist0_1
    type(cStopwatch)::                                        cswBlockRedist1_2
    type(cStopwatch)::                                        cswBlockRedist2_1
    type(cStopwatch)::                                        cswBlockRedist1_0
    type(cStopwatch)::                                        cswBlockObtain
    type(cStopwatch)::                                        cswBlockPut
    ! Maps
    type(cStopwatch)::                                        cswMaps
    type(cStopwatch)::                                        cswMapVt
    type(cStopwatch)::                                        cswMapVtInv
    type(cStopwatch)::                                        cswMapVxy
    type(cStopwatch)::                                        cswMapVxyInv
    ! Miscellaneous
    type(cStopwatch)::                                        cswMiscl
    type(cStopwatch)::                                        cswLinConv
    type(cStopwatch)::                                        cswNLinConv
    type(cStopwatch)::                                        cswDiskAcces
    ! Writing/Reading from disk
    type(cStopwatch)::                                        cswDisk
    type(cStopwatch)::                  cswDBG1
    type(cStopwatch)::                  cswDBG2
    type(cStopwatch)::                  cswDBG3
    type(cStopwatch)::                  cswDBG4
    type(cStopwatch)::                  cswDBG5

! *****************************************************************************
!
!   CONTAINED SUBROUTINES AND FUNCTIONS
!
!   StopWatchesInit  sub   Initializes all stopwatches to be used
!   StopWatchesStop  sub   Stop all stopwatches and write result to SW logfile
!   PrintToLog       sub   Print a message to the general log file
!   OpenLogFile      sub   Open the general log file
!   CloseLogFile     sub   Close the general log file
!   int2str          fun   translate an integer to a string in a standard
!                          format
!   alloctest        sub   Test whether a memory allocation has been performed
!                          well
!
! =============================================================================
!
CONTAINS

    SUBROUTINE StopWatchesInit(iProcID)

! =============================================================================
!
!   Programmer: Jasper de Koning
!
!   Language: Fortran 90
!
!   Version Date    Comment
!   ------- -----   -------
!   1.0     090505  Original code (JdK)
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The subroutine StopWatchesInit initializes the stopwatches, which are used
!   to measure the performance of the program. This function is called at the
!   start of the program. The Stopwatch log file is opened and the various
!   stopwatches are initialized
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   iProcID     i8b     processor ID
!
        integer(i8b), intent(in)::                iProcId; 
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
!   Stopwatch log file is opened and first sentence in SW log file is written.
!
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   OpenTimeFile
!   WriteToTimeFile
!   SWInit
!   SWStartAndCount
!
! =============================================================================

        ! Open Stopwatch log file and write first line
        call OpenTimeFile(iLogSWFileUNIT, trim(sOutputDir)//trim(sOutputRoot)//trim(int2str(iProcID))//".sw"); 
        call WriteToTimeFile(iLogSWFileUNIT, "Number of calls and execution time stopwatch results: ")

        ! Initialize stopwatches that measure parts of the iteration scheme
        call SWInit(cswProgram, iLogSWFileUNIT, "Total program"); 
        call SWInit(cswIniti, iLogSWFileUNIT, " Initialization")
        call SWInit(cswLinStep, iLogSWFileUNIT, " Linear step")
        call SWInit(cswCorrStep, iLogSWFileUNIT, " Correction step")
        call SWInit(cswNonLinContr, iLogSWFileUNIT, " Computation of the nonlinear correction")

        ! Initialize stopwatches that measure parts of the code
        ! Green's function
        call SWInit(cswGreen, iLogSWFileUNIT, "Green's function")
        ! The transforms
        call SWInit(cswTrans, iLogSWFileUNIT, "Total of Fourier Transforms")
        call SWInit(cswTransXY, iLogSWFileUNIT, " Transform XY")
        call SWInit(cswTransXYInv, iLogSWFileUNIT, " Transform XY Inverse")
        call SWInit(cswTransXYZ, iLogSWFileUNIT, " Transform XYZ")
        call SWInit(cswTransXYZInv, iLogSWFileUNIT, " Transform XYZ Inverse")
        call SWInit(cswTransGreenXY, iLogSWFileUNIT, " Transform Green XY")
        call SWInit(cswTransGreenXYZ, iLogSWFileUNIT, " Transform Green XYZ")
        call SWInit(cswTransT, iLogSWFileUNIT, " Transform T")
        call SWInit(cswTransTInv, iLogSWFileUNIT, " Transform T Inverse")
        ! Differentiation
        call SWInit(cswDeriv, iLogSWFileUNIT, "Computing the derivative")
        ! Block-related
        call SWInit(cswBlock, iLogSWFileUNIT, "Total of block-related costs")
        call SWInit(cswBlockRedist0_1, iLogSWFileUNIT, " Redistribute from 0 to 1")
        call SWInit(cswBlockRedist1_2, iLogSWFileUNIT, " Redistribute from 1 to 2")
        call SWInit(cswBlockRedist2_1, iLogSWFileUNIT, " Redistribute from 2 to 1")
        call SWInit(cswBlockRedist1_0, iLogSWFileUNIT, " Redistribute from 1 to 0")
        call SWInit(cswBlockObtain, iLogSWFileUNIT, " Obtain an xyz block")
        call SWInit(cswBlockPut, iLogSWFileUNIT, " Put an xyz block")
        ! Maps
        call SWInit(cswMaps, iLogSWFileUNIT, "Total of map-related costs")
        call SWInit(cswMapVt, iLogSWFileUNIT, " MapVt")
        call SWInit(cswMapVtInv, iLogSWFileUNIT, " MapVtInv")
        call SWInit(cswMapVxy, iLogSWFileUNIT, " MapVxy")
        call SWInit(cswMapVxyInv, iLogSWFileUNIT, " MapVxyInv")
        ! Miscellaneous
        call SWInit(cswDiskAcces, iLogSWFileUNIT, "Time spent on disk-access"); 
        call SWInit(cswLinConv, iLogSWFileUNIT, " Convolution in the linear step")
        call SWInit(cswNLinConv, iLogSWFileUNIT, " Convolution in the non-linear step")
        call SWInit(cswDBG1, iLogSWFileUNIT, " DBG 1")
        call SWInit(cswDBG2, iLogSWFileUNIT, " DBG 2")
        call SWInit(cswDBG3, iLogSWFileUNIT, " DBG 3")
        call SWInit(cswDBG4, iLogSWFileUNIT, " DBG 4")
        call SWInit(cswDBG5, iLogSWFileUNIT, " DBG 5")

        ! Start overhead cswProgram Stopwatch
        call SwStartAndCount(cswProgram); 
    END SUBROUTINE StopWatchesInit

    SUBROUTINE StopWatchesStop()

! =============================================================================
!
!   Programmer: Jasper de Koning
!
!   Language: Fortran 90
!
!   Version Date    Comment
!   ------- -----   -------
!   1.0     090505  Original code (JdK)
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The subroutine StopWatchesStop stops all stopwatches and outputs the
!   results to the Stopwatch log file.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   none
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
!   Total times and call counts recorded by all stopwatches is output to the
!   Stopwatch log file, and it is closed.
!
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   SWStop
!   SWOutputToFile
!   CloseTimeFile
!
! =============================================================================

        ! First we write down the variables, which measure parts of the
        ! iteration scheme
        call SWStop(cswProgram); 
        call SWOutputToFile(cswProgram); 
        call SWOutputToFile(cswLinStep); 
        call SWOutputToFile(cswCorrStep); 
        call SWOutputToFile(cswNonLinContr); 
        ! Secondly, we write down the variables, which measure parts of the code
        ! Green's function
        call SWOutputToFile(cswGreen); 
        ! The transforms
        call SWOutputToFile(cswTrans); 
        call SWOutputToFile(cswTransXY); 
        call SWOutputToFile(cswTransXYInv); 
        call SWOutputToFile(cswTransXYZ); 
        call SWOutputToFile(cswTransXYZInv); 
        call SWOutputToFile(cswTransGreenXY); 
        call SWOutputToFile(cswTransGreenXYZ); 
        call SWOutputToFile(cswTransT); 
        call SWOutputToFile(cswTransTInv); 
        ! Differentiation
        call SWOutputToFile(cswDeriv); 
        ! Block-related
        call SWOutputToFile(cswBlock); 
        call SWOutputToFile(cswBlockRedist0_1); 
        call SWOutputToFile(cswBlockRedist1_2); 
        call SWOutputToFile(cswBlockRedist2_1); 
        call SWOutputToFile(cswBlockRedist1_0); 
        call SWOutputToFile(cswBlockObtain); 
        call SWOutputToFile(cswBlockPut); 
        ! Maps
        call SWOutputToFile(cswMapVt); 
        call SWOutputToFile(cswMapVtInv); 
        call SWOutputToFile(cswMapVxy); 
        call SWOutputToFile(cswMapVxyInv); 
        ! Miscellaneous
        call SWOutputToFile(cswLinConv)
        call SWOutputToFile(cswNLinConv)
        call SWOutputToFile(cswDiskAcces); 
        call SWOutputToFile(cswDBG1)
        call SWOutputToFile(cswDBG2)
        call SWOutputToFile(cswDBG3)
        call SWOutputToFile(cswDBG4)
        call SWOutputToFile(cswDBG5)

        call CloseTimeFile(iLogSWFileUNIT); 
    END SUBROUTINE StopWatchesStop

    SUBROUTINE PrintToLog(cString, iDebugRequest)

! =============================================================================
!
!   Programmer: Jasper de Koning
!
!   Language: Fortran 90
!
!   Version Date    Comment
!   ------- -----   -------
!   1.0     090505  Original code (JdK)
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The function PrintToLog echoes a certain string to the general log file.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cString         char   Text to be output
!   iDebugRequest   i4b    Debuglevel of the text, determining the indentation
!                          of the text in the file
!
        character(LEN=*), intent(in)::                        cString
        integer(i4b), intent(in)::                iDebugRequest
!
! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   iStartmessage  i8b   start index of the string in the text line
!   iEndmessage    i8b   end index of the string in the text line
!   iCtr           i8b   loop counter
!   rTimeStap      dp    time stamp to be echoed to the log file
!   rMemMB         dp    amount of memory in MB to be echoed to the log file
!   cMessage       char  total text line to the file
!
        integer(i8b):: iStartmessage, iEndmessage, iCtr
        real(dp) :: rTimeStamp, rMemMB
        character(LEN=120) :: cMessage

! *****************************************************************************
!
!   I/O
!
!   If the DebugRequest>Debuglevel, the string is output to the general log
!   file, together with the time stamp and the current memory usage.
!
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   TimeElapsed
!
! =============================================================================

!#ifndef RUNESTIMATOR
        if (iDebugRequest .le. iDebugLvl) then

            iStartmessage = max(4, 4 + 4*iDebugRequest)
            iEndmessage = min(iStartmessage + len_trim(cString) - 1, 120_i8b)
            do iCtr = 1, iStartmessage - 1
                cMessage(iCtr:iCtr) = char(32)
            end do
            cMessage(iStartmessage:iEndmessage) = cString(1:(iEndmessage - iStartmessage + 1))
            if (iEndmessage < 120) then
                do iCtr = iEndmessage + 1, 120
                    cMessage(iCtr:iCtr) = char(32)
                end do
            end if

            rTimeStamp = real(TimeElapsed(iTimeStamp))/100.0
            rMemMB = real(iMemAllocated)/(1024**2)

            write (iLogFileUNIT, '(A5, F10.2, A6, F10.2,A1, A120)') '|Mem ', rMemMB, ' |Time', rTimeStamp, '|', cMessage; 
!#ifndef NOFLUSH_LOG
            if (iLogFileUNIT /= 6) then
                flush (iLogFileUNIT)
            end if
!#endif

        end if
!#endif
    END SUBROUTINE

    SUBROUTINE OpenLogFile(iProcID)

! =============================================================================
!
!   Programmer: Jasper de Koning
!
!   Language: Fortran 90
!
!   Version Date    Comment
!   ------- -----   -------
!   1.0     090505  Original code (JdK)
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The subroutine OpenLogFile opens the general log file. No error checking
!   included.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   iProcID    i8b  Processor Identification number
!
        integer(i8b), intent(in) :: iProcID

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   sSuffix     char    suffix at the end of the log filename identifying the
!                       processor for which the log file is for.
!
        character(len=4)::                sSuffix; 
! *****************************************************************************
!
!   I/O
!
!   Opens the general log file; if iProcID==0, output log file to screen
!
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   GetIntTimeNow
!   int2str
!
! =============================================================================

        iTimeStamp = GetIntTimeNow(); 
        sSuffix = int2str(iProcID)

        if (iProcID > 0) then
            open (UNIT=iLogFileUNIT, STATUS='REPLACE', FILE=trim(sOutputDir)//trim(sOutputRoot)//trim(sSuffix)//".log"); 
        else
            iLogFileUNIT = 6 !Output to screen
        end if

    END SUBROUTINE

    SUBROUTINE CloseLogFile

! =============================================================================
!
!   Programmer: Jasper de Koning
!
!   Language: Fortran 90
!
!   Version Date    Comment
!   ------- -----   -------
!   1.0     090505  Original code (JdK)
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The subroutine CloseLogFile closes the general log file. No error checking
!   included.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   none
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
!   Closes general log file
!
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   none
!
! =============================================================================

        close (iLogFileUNIT)

    END SUBROUTINE

    character(LEN=4) FUNCTION int2str(val)

! =============================================================================
!
!   Programmer: Koos Huijssen
!
!   Language: Fortran 90
!
!   Version Date    Comment
!   ------- -----   -------
!   1.0     090505  Original code (KH)
!
!
!   Copyright (C)   Laboratory of Electromagnetic Research
!                   Delft University of Technology, Delft, The Netherlands
!
!   Permission to copy or distribute this software or documentation
!   in hard copy or soft copy granted only by written license
!   obtained from the Laboratory of Electromagnetic Research.
!   All rights reserved. No parts of this publication may be
!   reproduced, stored in a retrieval system (e.g. in memory, disk,
!   or core) or be transmitted by any means, electronic, mechanical,
!   photocopy, recording, or otherwise, without written permission
!   from the publisher.
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The function int2str returns the string representation of an integer,
!   always using three digits (zeros or numbers) and preceded by an '-'
!   or '_' depending on whether the number is positive or negative.
!   If abs(val)>999, then only the smallest three digits are used.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   val  i8b   Value to be translated
!
        integer(i8b), INTENT(IN) :: val
!
! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   absval      i8b    absolute of val
!   i1          i8b    4th digit of val
!   i2          i8b    3rd digit of val
!   i3          i8b    2nd digit of val
!   i4          i8b    1st digit of val
!   plusminus   char   either '-' or '_' depending on the sign of val
!
        integer(i8b) :: absval, i1, i2, i3, i4
        character(1) :: plusminus

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

        if (val < 0) then
            plusminus = '-'
        else
            plusminus = '_'
        end if

        absval = abs(val)

        i1 = absval/1000
        i2 = (absval - i1*1000)/100
        i3 = (absval - i1*1000 - i2*100)/10
        i4 = absval - i1*1000 - i2*100 - i3*10

        int2str = trim(plusminus)//char(48 + i2)//char(48 + i3)//char(48 + i4)

    END FUNCTION int2str

    SUBROUTINE alloctest(error, errtext)

! =============================================================================
!
!   Programmer: Koos Huijssen
!
!   Language: Fortran 90
!
!   Version Date    Comment
!   ------- -----   -------
!   1.0     090505  Original code (KH)
!
!
!   Copyright (C)   Laboratory of Electromagnetic Research
!                   Delft University of Technology, Delft, The Netherlands
!
!   Permission to copy or distribute this software or documentation
!   in hard copy or soft copy granted only by written license
!   obtained from the Laboratory of Electromagnetic Research.
!   All rights reserved. No parts of this publication may be
!   reproduced, stored in a retrieval system (e.g. in memory, disk,
!   or core) or be transmitted by any means, electronic, mechanical,
!   photocopy, recording, or otherwise, without written permission
!   from the publisher.
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The subroutine alloctest tests whether an allocation has been successfully
!   performed. If not, it echoes an error text and makes the program to exit.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   error     i8b    error value that has been returned by an ALLOCATE command
!   errtext   char   text to be echoed if an error has occurred
!
        integer(i8b), INTENT(in) :: error
        character(LEN=*), OPTIONAL, INTENT(in) :: errtext

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   acTemp   char    text to be echoed to the log file
!
        character(LEN=120) :: acTemp

! *****************************************************************************
!
!   I/O
!
!   If error /=0, Echo error number and errtext to log file
!
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   PrintToLog
!
! =============================================================================

        if (error /= 0) then
            write (acTemp, "('error ',I5,': memory allocation failed')") error
            call PrintToLog(acTemp, -1)
            if (present(errtext)) then
                call PrintToLog(errtext, -1); 
            end if
            stop
        end if

    END SUBROUTINE alloctest

END MODULE ParnacGeneral
