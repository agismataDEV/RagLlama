MODULE ParnacDataRedist

! =============================================================================
!
!   Programmer: Jasper de Koning / Koos Huijssen
!
!   Language: Fortran 90
!
!   Version Date    Comment
!   ------- -----   -------
!   1.0     090505  Original code (KH)
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
!   The module ParnacDataRedist contains subroutines that redistribute the data
!   in the grid structure to distribution 0, 1 and 2 and vice versa. The most
!   complex one is the transition from distr. 1 to 2 and vv, as it concerns a
!   communication of the data between all processors. Also contained are a
!   number of subroutines that copy a block of data for a single frequency
!   from a Grid in distribution 2.
!
! *****************************************************************************
!
!   MODULES USED
!
	USE Types
    USE Constants
	USE ParnacGeneral
	USE ParnacDataDef
!    USE ParnacMPIDef

	USE ParnacDataSpaceInit, ONLY : &
        GridDistr0CreateEmpty, GridDistr0Allocate, GridDistr0Deallocate, &
        GridDistr1CreateEmpty, GridDistr1Allocate, GridDistr1Deallocate, &
        GridDistr2Allocate, GridDistr2Deallocate

! *****************************************************************************
!
!   GLOBAL DECLARATIONS
!
!   none
!
IMPLICIT NONE

! *****************************************************************************
!
!   CONTAINED SUBROUTINES AND FUNCTIONS
!
!   ReorderDistr1ToDistr2   sub   reorder a grid from distribution 1 to 2
!   ReorderDistr2ToDistr1   sub   reorder a grid from distribution 2 to 1
!   ReorderDistr0ToDistr1   sub   reorder a grid from distribution 0 to 1
!   ReorderDistr1ToDistr0   sub   reorder a grid from distribution 1 to 0
!   Distr2ObtainXYZBlock    sub   obtain a single frequency XYZ block from a 
!                                 grid
!   Distr2ObtainYMirroredXYZBlock  sub  obtain a single frequency XYZ block
!                                 from a grid in mirrored form
!   Distr2PutXYZBlock       sub   put a single frequency XYZ block in a grid
!   Distr2AddToXYZBlock     sub   add a single frequency XYZ block to a block
!                                 in a grid
!
! =============================================================================
!
CONTAINS

SUBROUTINE ReorderDistr1ToDistr2(cGrid)
	
! =============================================================================
!
!   Programmer: Jasper de Koning / Koos Huijssen
!
!   Language: Fortran 90
!
!   Version Date    Comment
!   ------- -----   -------
!   1.0     090505  Original code (KH)
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The subroutine ReorderDistr1ToDistr2 reorders the data in a grid from 
!   distribution 1 to distr. 2, i.e. from a T-local distribution to an XYZ-
!   local distribution. This involves a communication of the data between all
!   processors.
!
!
!   Memory usage: Currently this implementation uses quite a lot of memory, as
!   it involves a buffer which has the size of the original grid. If possible,
!   improving this will remove a bottleneck in the memory usage. The only
!   way to do this would be to use asymmetrical data transfer and complex inter-
!   processor communication.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cGrid   io   type(Grid)   The grid that contains the data array which needs
!                             to be redistributed.
!
	type(Grid), INTENT(inout) ::		cGrid;

! *****************************************************************************
! 
!   LOCAL PARAMETERS      
!
!   TotalProc   i8b   Total number of processors
!   Local1T     i8b   Total T-positions on each processor in Distr. 1
!   Local2T     i8b   Total T-positions on each processor in Distr. 2
!   Local1XYZ   i8b   Total XYZ-positions on each processor in Distr. 1
!   Local2XYZ   i8b   Total XYZ-positions on each processor in Distr. 2
!   iProcL      i4b   Loop counter over the processors
!   iT2L        i8b   Loop counter over the T-positions in Distr. 2
!   iXYZ1L      i8b   Loop counter over the XYZ-positions in Distr. 1
!   iSrcInd     i8b   Index in the source data array
!   iDestInd    i8b   Index in the destination data array
!   iErr        i4b   Error number
!   acTemp      char  Temporary char array for output log messages
!
	integer(i8b) ::				TotalProc, Local1T, Local2T, Local1XYZ, Local2XYZ
	integer(i8b) ::				iProcL, iT2L, iXYZ1L
	integer(i8b) ::				iSrcInd, iDestInd
	integer(i4b) ::				iErr
	character(LEN = 2048) ::		acTemp;
!	integer(i4b)::			aiStatus(MPI_STATUS_SIZE) 
!	integer(i4b)::			aiReqSent(cGrid.iProcN), aiReqRecv(cGrid.iProcN);
	
! *****************************************************************************
!
!   I/O
!
!   Log file entries 
!   
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   SWstart
!   SWstartandcount
!   SWstop
!   MPI_Barrier
!   MPI_Alltoall
!   GridDistr2Allocate
!   GridDistr1Deallocate
!   PrintToLog
!   
! *****************************************************************************
!
!   PSEUDOCODE
!
!   Allocate the Distr. 2 array as a buffer
!   Pack the data in Distr. 1 such that we have data blocks of size 
!     local2T*local1XYZ, to facilitate the redistribution of contiguous blocks.
!     the stride in T is still 1...
!   Redistribute the data across all processors to change from the T-local
!     distribution 1 to the XYZ-local distribution 2
!   Unpack the data in Distr. 2 such that the stride in XYZ becomes 1
!   Deallocate the Distr. 1 array
!   
! =============================================================================
	
	write (acTemp, '("ReorderDistr1ToDistr2")');	call PrintToLog(acTemp, 3);

	call SWStart(cswDBG1);
!KH Is this one necessary since the MPI_WAIT is also included???
!KH It does make it easy to check the waiting time and therefore
!KH the processor load imbalance
    call MPI_Barrier(MPI_COMM_WORLD, iErr);
	call SWStop(cswDBG1);

	call SwStartAndCount(cswBlock);	call SwStartAndCount(cswBlockRedist1_2)
	
	! If it is not in the correct distribution, we can do nothing and return an error
	if (cGrid.iDistr /= 1) then
			write (acTemp, '("ReorderDistr1ToDistr2: the grid provided as input is not in distribution 1, but in distribution ", I3, "")') cGrid%iDistr;
			call PrintToLog(acTemp, -1);
		stop
	end if

    ! If iProcN = 1, then the grid is in T-local as well as XYZ-local distribution,
    ! no need to redistribute
    if (cGrid.iProcN>1) then

	    write (acTemp, '("Allocate mem for redistribution")');	call PrintToLog(acTemp, 4);
	    ! First we allocate cGrid.pacD2, we use this as a buffer
	    call GridDistr2Allocate(cGrid);

	    !Initialize redistribution variables (for an easier understanding of the redistribution process)
	    TotalProc = cGrid.iProcN
	    Local1T   = cGrid.iD2GlobN	! = cGrid.iD1TL = cSpace%iDimT+1
	    Local2XYZ = cGrid.iD1GlobN	! = cGrid.iD1GlobN = cSpace%iDimX * cSpace%iDimY * cSpace%iDimZ
	    Local1XYZ = cGrid.iD1LocN	! = Local2XYZ/cGrid.iProcN
	    Local2T   = cGrid.iD2LocN	! = Local1T/cGrid.iProcN
	
	    write (acTemp, '(" Pack")');	call PrintToLog(acTemp, 4);
	    ! We continue by packing the data into the second buffer
	    do iProcl = 0, TotalProc-1
		    do iXYZ1l = 0, Local1XYZ-1
			    ! the following most inner loop is not needed, since we can use vectorized copying
			    ! but it is added for clarity
			    !do iT2l = 0, Local2T-1

			    iSrcInd		= 1 + iXYZ1l  * Local1T			 + iProcl * Local2T  ! + iT2l
			    iDestInd	= 1 + (iProcl * Local1XYZ + iXYZ1l) * Local2T		 ! + iT2l	
			    
			    cGrid.pacD2(iDestInd:iDestInd+Local2T-1) 	= cGrid.pacD1(iSrcInd:iSrcInd+Local2T-1);

			    !end do
		    end do
	    end do
		    
	    ! Now we have allocated the buffer space and have packed the data, it is time to start with 
	    ! the distributing itself. We have to send one block of data from each processor, to each other processor. 
        ! This is done by the MPI_Alltoall primitive.
	    write (acTemp, '(" Communicate")');	call PrintToLog(acTemp, 4);

        call MPI_Alltoall(cGrid.pacD2,              &
                          int(Local1XYZ*Local2T),   &
                          MPI_DOUBLE_COMPLEX,       &
                          cGrid.pacD1,              &
                          int(Local1XYZ*Local2T),   &
                          MPI_DOUBLE_COMPLEX,       &
                          MPI_COMM_WORLD,           &
                          iErr)
	    
! MPI_Alltoall is an out-of-place command, which means that the redistribution 
! of the data needs two buffers to be allocated simultaneously. This makes it 
! expensive in terms of memory usage. To improve this, in the future we may 
! adopt an alternative implementation; maybe the implementation using non-
! blocking Sends and Receives that is included below (and that is also out-
! of-place) aids in the implementation.
!
!	    do iProcl = 0, TotalProc-1
!		    call MPI_ISend(cGrid.pacD2(1 + iProcl * Local1XYZ*Local2T),  &
!				    int(Local1XYZ*Local2T), &
!				    MPI_DOUBLE_COMPLEX, &
!				    int(iProcl),	&
!				    int(cGrid.iProcID),	&
!				    MPI_COMM_WORLD, &
!				    aiReqSent(iProcl+1), &
!				    iErr);
!           end do
!            do iProcl = 0, TotalProc-1
!		    call MPI_IRecv(cGrid.pacD1(1 + iProcl * Local1XYZ*Local2T), &
!				    int(Local1XYZ*Local2T), &
!				    MPI_DOUBLE_COMPLEX, &
!				    int(iProcl), &
!				    int(iProcl), &
!				    MPI_COMM_WORLD, &
!				    aiReqRecv(iProcl+1), &
!				    iErr);
!	    end do
!	    
!	    write (acTemp, '(" Synchronize")');	call PrintToLog(acTemp, 4);
!	    ! We wait untill we have obtained and sent all the blocks
!	    do iProcl = 0, TotalProc-1
!		    call MPI_Wait(aiReqRecv(iProcl+1), aiStatus, iErr);
!		    call MPI_Wait(aiReqSent(iProcl+1), aiStatus, iErr);
!	    end do

	    write (acTemp, '(" Unpack")');	call PrintToLog(acTemp, 4);
	    ! Now we have distributed all the data, it will still be a mess, it is not yet in distribution 2. At this
	    !  moment the data is still distributed such that the stride in T-Direction is 1. So, to finish this job,
	    !  we have to redistribute this data again, but now locally.
	    !
	    ! The index of the destination can be written as:
	    !  A(x, y, z, i) = Av(x + y * Sy + z * Sz + i * Si), the data as we get it is not in a nice order
	    !  , we can find the value of i easily (this is T), it is hard however to obtain the 
	    !  x, y, z values. However, we can find x + y * Sy + z * Sz = iBl + iPl * cGrid.iD0LocN,
	    !  because the t-blocks were block distributed over the processes in distribution 0,1.
	    do iProcl = 0, TotalProc-1
		    do iXYZ1l = 0, Local1XYZ-1
			    do iT2l = 0, Local2T-1
				    
				    iSrcInd			= 1 + (iProcl * Local1XYZ + iXYZ1l) * Local2T + iT2l;	
				    iDestInd		= 1 + iT2l * Local2XYZ + (iProcl * Local1XYZ + iXYZ1l) 

				    cGrid.pacD2(iDestInd)	= cGrid.pacD1(iSrcInd);
				    
			    end do
		    end do
	    end do

	    ! We finish by deallocating the space that was used by distribution 1
	    call GridDistr1DeAllocate(cGrid);
	
    end if

	! And finally, we can set the flag to 2!
	cGrid.iDistr		= 2;

!	call SWStart(cswDBG1);
!!KH Is this one necessary since the MPI_WAIT is also included???
!    call MPI_Barrier(MPI_COMM_WORLD, iErr);
!	call SWStop(cswDBG1);
	call SWStop(cswBlock);	call SWStop(cswBlockRedist1_2)
END SUBROUTINE ReorderDistr1ToDistr2

SUBROUTINE ReorderDistr2ToDistr1(cGrid)
	
! =============================================================================
!
!   Programmer: Jasper de Koning / Koos Huijssen
!
!   Language: Fortran 90
!
!   Version Date    Comment
!   ------- -----   -------
!   1.0     090505  Original code (KH)
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The subroutine ReorderDistr2ToDistr1 performs the inverse operation of 
!   ReorderDistr1ToDistr2. See its header for more info.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cGrid   io   type(Grid)   The grid that contains the data array which needs
!                             to be redistributed.
!
	type(Grid), INTENT(inout) ::		cGrid;

! *****************************************************************************
! 
!   LOCAL PARAMETERS      
!
!   TotalProc   i8b   Total number of processors
!   Local1T     i8b   Total T-positions on each processor in Distr. 1
!   Local2T     i8b   Total T-positions on each processor in Distr. 2
!   Local1XYZ   i8b   Total XYZ-positions on each processor in Distr. 1
!   Local2XYZ   i8b   Total XYZ-positions on each processor in Distr. 2
!   iProcL      i4b   Loop counter over the processors
!   iT2L        i8b   Loop counter over the T-positions in Distr. 2
!   iXYZ1L      i8b   Loop counter over the XYZ-positions in Distr. 1
!   iSrcInd     i8b   Index in the source data array
!   iDestInd    i8b   Index in the destination data array
!   iErr        i4b   Error number
!   acTemp      char  Temporary char array for output log messages
!
	integer(i8b) ::				TotalProc, Local1T, Local2T, Local1XYZ, Local2XYZ
	integer(i8b) ::				iProcl,    iT2l, iXYZ1L
	integer(i8b) ::				iSrcInd, iDestInd
	integer(i4b) ::				iErr;
	character(LEN = 2048) ::		acTemp;
!	integer(i4b)::			aiStatus(MPI_STATUS_SIZE) 
!	integer(i4b)::			aiReqSent(cGrid.iProcN), aiReqRecv(cGrid.iProcN);
	
! *****************************************************************************
!
!   I/O
!
!   Log file entries 
!   
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   SWstart
!   SWstartandcount
!   SWstop
!   MPI_Barrier
!   MPI_Alltoall
!   GridDistr2Allocate
!   GridDistr1Deallocate
!   PrintToLog
!   
! =============================================================================
	
	write (acTemp, '("ReorderDistr2ToDistr1")');	call PrintToLog(acTemp, 3);

	call SWStart(cswDBG1);
!KH Is this one necessary since the MPI_WAIT is also included???
    call MPI_Barrier(MPI_COMM_WORLD, iErr);
	call SWStop(cswDBG1);

	call SwStartAndCount(cswBlock);	call SwStartAndCount(cswBlockRedist2_1)

	! If it is not in the correct distribution, we can do nothing and return an error
	if (cGrid.iDistr /= 2) then
			write (acTemp, '("ReorderDistr2ToDistr1: the grid provided as input is not in distribution 2, but in distribution ", I3, "")') cGrid%iDistr;
			call PrintToLog(acTemp, -1);
		stop
	end if

    ! If iProcN = 1, then the grid is in T-local as well as XYZ-local distribution,
    ! no need to redistribute
    if (cGrid.iProcN>1) then

	    ! First we allocate cGrid.acD1, we use this as a buffer
	    write (acTemp, '("Allocate mem for redistribution")');	call PrintToLog(acTemp, 4);
	    call GridDistr1Allocate(cGrid);

	    !Initialize redistribution variables (for better understanding of the redistribution)
	    TotalProc = cGrid.iProcN
	    Local1T   = cGrid.iD1TL						! = cGrid.iD2GlobN
	    Local2T   = cGrid.iD1TL/cGrid.iProcN		! = cGrid.iD2LocN
	    Local1XYZ = cGrid.iD1GlobN/cGrid.iProcN		! = cGrid.iD1LocN
	    Local2XYZ = cGrid.iD1GlobN

	    ! We continue by packing the data into the second buffer
	    write (acTemp, '(" Pack")');	call PrintToLog(acTemp, 4);
	    do iProcl = 0, TotalProc-1
		    do iXYZ1l = 0, Local1XYZ-1
			    do iT2l = 0, Local2T-1
				    
				    iSrcInd		= 1 + iT2l * Local2XYZ + (iProcl * Local1XYZ + iXYZ1l) 
				    iDestInd	= 1 + (iProcl * Local1XYZ + iXYZ1l) * Local2T + iT2l;	

				    cGrid.pacD1(iDestInd)	= cGrid.pacD2(iSrcInd);
				    
			    end do
		    end do
	    end do

	    ! Now we have allocated the buffer space and have packed the data, it is time to start with 
	    ! the distributing itself. We have to send one block of data from each processor, to each other processor. 
	    ! We will do this using synchronous sends (ISend). We did not opt for the usage of MPI_Alltoall, because
	    ! this operation is not in-place. This command is not inplace either, but we could fiddle around with it for a bit,
	    ! such that it might BECOME inplace. Therefore, to keep options open we use this method.
	    write (acTemp, '(" Communicate")');	call PrintToLog(acTemp, 4);

        call MPI_Alltoall(cGrid.pacD1,              &
                          int(Local1XYZ*Local2T),   &
                          MPI_DOUBLE_COMPLEX,       &
                          cGrid.pacD2,              &
                          int(Local1XYZ*Local2T),   &
                          MPI_DOUBLE_COMPLEX,       &
                          MPI_COMM_WORLD,           &
                          iErr)

    !	if (cGrid.iProcN > 1) then
    !		call SWStart(cswDBG1);
    !!KH Is this one necessary since the MPI_WAIT is also included???
    !        call MPI_Barrier(MPI_COMM_WORLD, iErr);
    !		call SWStop(cswDBG1);
    !	end if
    !	do iProcl = 0, TotalProc-1
    !		call MPI_ISend(cGrid.pacD1(1 + iProcl * Local1XYZ*Local2T),  &
    !				int(Local1XYZ*Local2T), &
    !				MPI_DOUBLE_COMPLEX, &
    !				int(iProcl),	&
    !				int(cGrid.iProcID),	&
    !				MPI_COMM_WORLD, &
    !				aiReqSent(iProcl+1), &
    !				iErr);
    !		call MPI_IRecv(cGrid.pacD2(1 + iProcl * Local1XYZ*Local2T), &
    !				int(Local1XYZ*Local2T), &
    !				MPI_DOUBLE_COMPLEX, &
    !				int(iProcl), &
    !				int(iProcl), &
    !				MPI_COMM_WORLD, &
    !				aiReqRecv(iProcl+1), &
    !				iErr);
    !	end do
    !	
    !	! We wait untill we have obtained and sent all the blocks
    !	write (acTemp, '(" Synchronize")');	call PrintToLog(acTemp, 4);
    !	do iProcl = 0, TotalProc-1
    !		call MPI_Wait(aiReqRecv(iProcl+1), aiStatus, iErr);
    !		call MPI_Wait(aiReqSent(iProcl+1), aiStatus, iErr);
    !	end do
	    
	    ! Now we have distributed all the data, it will still be a mess, it is not yet in distribution 1. At this
	    !  moment the data is still distributed such that the stride in T-Direction is big. So, to finish this job
	    !  we have to redistribute this data again, but now locally.
	    !
	    ! The index of the destination can be written as: A(t, i) = Av(t + i * Si), the data as we get it is not 
	    !  in a nice order. It is stored as: A_org(t, i) = A_org(p*cGrid.iD1TL/cGrid.iProcN + t, i), with p
	    !  the ID of the processor that sent the data. So, we get:	!
	    write (acTemp, '(" Unpack")');	call PrintToLog(acTemp, 4);
	    do iProcl = 0, TotalProc-1
		    do iXYZ1l = 0, Local1XYZ-1
			    ! the following most inner loop is not needed, since we can use vectorized copying
			    ! but it is added for clarity
			    !do iT2l = 0, Local2T-1

			    iSrcInd		= 1 + (iProcl * Local1XYZ + iXYZ1l) * Local2T		 ! + iT2l	
			    iDestInd	= 1 + iXYZ1l  * Local1T			 + iProcl * Local2T  ! + iT2l
			    
			    cGrid.pacD1(iDestInd:iDestInd+Local2T-1) 	= cGrid.pacD2(iSrcInd:iSrcInd+Local2T-1);

			    !end do
		    end do
	    end do

	    ! We finish by deallocating the space that was used by distribution 2
	    call GridDistr2DeAllocate(cGrid);
	
    end if

	! And finally, we can set the flag to 1!
	cGrid.iDistr		= 1;
!call SWStart(cswDBG1);
!!KH Is this one necessary since the MPI_WAIT is also included???
!    call MPI_Barrier(MPI_COMM_WORLD, iErr);
!call SWStop(cswDBG1);
	call SWStop(cswBlock);	call SWStop(cswBlockRedist2_1)
END SUBROUTINE ReorderDistr2ToDistr1

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %
! % 		ReorderDistr0ToDistr1
! %
! % Transform a cGrid from Distribution 0 to Distribution 1. We do not deallocate the data in
! % distribution 0. 
! % Input:
! % 	cGrid: this grid should obviously be in distribution 0. 
! % Output 
! % 	cGrid: The same data, but now stored in distribution 1. If pacD1 was already allocated
! %		the data stored in it is overwritten and is lost for all eternity. The data
! %		in parD0 will not be deallocated however.
! %
! %	We assume that cGrid.iD1T = 1 <--- Should be better
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE ReorderDistr0ToDistr1(cGrid)
	
! =============================================================================
!
!   Programmer: Jasper de Koning / Koos Huijssen
!
!   Language: Fortran 90
!
!   Version Date    Comment
!   ------- -----   -------
!   1.0     090505  Original code (KH)
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The subroutine ReorderDistr0ToDistr1 reorders the data in a grid from 
!   distribution 0 to distr. 1, i.e. from a T-local distribution, real numbers,
!   no wraparound region in T, to a T-local distribution, real numbers stored 
!   in complex data type, including wraparound region in T. It is memory 
!   intensive as it has two buffers in memory simultaneously, however it is 
!   never a bottleneck as ReorderDistr1ToDistr2 takes much more memory.
!
!   One unnecessary assumption in the current implementation is that 
!   cGrid.iD1T = 1. This could be generalized but this is not very urgent.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cGrid   io   type(Grid)   The grid that contains the data array which needs
!                             to be redistributed.
!
	type(Grid), intent(inout)::	cGrid

! *****************************************************************************
! 
!   LOCAL PARAMETERS      
!
!   iLI      i8b   Loop counter over all XYZ-positions
!   acTemp   char  Temporary char array for output log messages
!
	integer(i8b) ::			iLT, iLI
	character(LEN = 2048) ::		acTemp;
	
! *****************************************************************************
!
!   I/O
!
!   Log file entries 
!   
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   SWstartandcount
!   SWstop
!   GridDistr1CreateEmpty
!   GridDistr0Deallocate
!   PrintToLog
!   
! =============================================================================
	
	call SwStartAndCount(cswBlock);	call SwStartAndCount(cswBlockRedist0_1)
	write (acTemp, '("ReorderDistr0ToDistr1")');	call PrintToLog(acTemp, 3);
	! First we will check whether we are in distribution 0
	if (cGrid.iDistr /= 0) then
		write (acTemp, '("ReorderDistr0ToDistr1: the grid provided as input is not in distribution 0, but in distribution ", I3, "")') cGrid%iDistr;
		call PrintToLog(acTemp,0)
		stop
	end if
	
	! Allocate the data for distribution 1
	call GridDistr1CreateEmpty(cGrid)
	if(mod(cGrid.iD0TL,2_i8b)==1) then
		do iLI = 0, cGrid.iD0LocN-1
			cGrid.pacD1(1 + iLi*cGrid.iD1IS : 1 + iLi*cGrid.iD1IS + (cGrid.iD0TL/2-1)*cGrid.iD1TS : cGrid.iD1TS) =&
					cGrid.parD0(1 + iLi*cGrid.iD0IS : 1 + iLi*cGrid.iD0IS + (cGrid.iD0TL-2)*cGrid.iD0TS : 2*cGrid.iD0TS) + &
				im* cGrid.parD0(2 + iLi*cGrid.iD0IS : 1 + iLi*cGrid.iD0IS + (cGrid.iD0TL-2)*cGrid.iD0TS : 2*cGrid.iD0TS)
			cGrid.pacD1(1 + iLi*cGrid.iD1IS + ((cGrid.iD0TL+1)/2-1)*cGrid.iD1TS) =&
					cGrid.parD0(1 + iLi*cGrid.iD0IS + (cGrid.iD0TL-1)*cGrid.iD0TS);
		end do
	else
		do iLI = 0, cGrid.iD0LocN-1
			cGrid.pacD1(1 + iLi*cGrid.iD1IS : 1 + iLi*cGrid.iD1IS + (cGrid.iD0TL/2-1)*cGrid.iD1TS : cGrid.iD1TS) =&
					cGrid.parD0(1 + iLi*cGrid.iD0IS : 1 + iLi*cGrid.iD0IS + (cGrid.iD0TL-1)*cGrid.iD0TS : 2*cGrid.iD0TS) + &
				im* cGrid.parD0(2 + iLi*cGrid.iD0IS : 1 + iLi*cGrid.iD0IS + (cGrid.iD0TL-1)*cGrid.iD0TS : 2*cGrid.iD0TS)
		end do
	end if
    call GridDistr0DeAllocate(cGrid);
    
	! Set the flag to 1
	cGrid.iDistr		= 1;

	call SWStop(cswBlock);	call SWStop(cswBlockRedist0_1)
END SUBROUTINE ReorderDistr0ToDistr1

SUBROUTINE ReorderDistr1ToDistr0(cGrid)
	
! =============================================================================
!
!   Programmer: Jasper de Koning / Koos Huijssen
!
!   Language: Fortran 90
!
!   Version Date    Comment
!   ------- -----   -------
!   1.0     090505  Original code (KH)
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The subroutine ReorderDistr1ToDistr0 performs the inverse operation of 
!   ReorderDistr0ToDistr1. See its header for more info.
!
!   One unnecessary assumption in the current implementation is that 
!   cGrid.iD1T = 1. This could be generalized but this is not very urgent.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cGrid   io   type(Grid)   The grid that contains the data array which needs
!                             to be redistributed.
!
	type(Grid), intent(inout)::	cGrid

! *****************************************************************************
! 
!   LOCAL PARAMETERS      
!
!   iLI      i8b   Loop counter over all XYZ-positions
!   acTemp   char  Temporary char array for output log messages
!
	integer(i8b) ::			iLT, iLI
	character(LEN = 2048) ::		acTemp;
	
! *****************************************************************************
!
!   I/O
!
!   Log file entries 
!   
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   SWstartandcount
!   SWstop
!   GridDistr0CreateEmpty
!   GridDistr1Deallocate
!   PrintToLog
!   
! =============================================================================

	write (acTemp, '("ReorderDistr1ToDistr0")');	call PrintToLog(acTemp, 3);
	
	call SwStartAndCount(cswBlock);	call SwStartAndCount(cswBlockRedist1_0)
	! First we will check whehter we are in distribution 1
	if (cGrid.iDistr /= 1) then
		write (acTemp, '("ReorderDistr1ToDistr0: the grid provided as input is not in distribution 1, but in distribution ", I3, "")') cGrid%iDistr;
		call PrintToLog(acTemp,0)
		stop
	end if

	! Allocate the data for distribution 0
	call GridDistr0CreateEmpty(cGrid)

	if(mod(cGrid.iD0TL,2_i8b)==1) then
		do iLi = 0, cGrid.iD0LocN-1
			cGrid.parD0(1 + iLi*cGrid.iD0IS : 1 + iLi*cGrid.iD0IS + (cGrid.iD0TL-1)*cGrid.iD0TS: 2*cGrid.iD0TS) = &
			dreal(cGrid.pacD1(1 + iLi*cGrid.iD1IS : 1 + iLi*cGrid.iD1IS + ((cGrid.iD0TL+1)/2-1)*cGrid.iD1TS : cGrid.iD1TS))
			cGrid.parD0(2 + iLi*cGrid.iD0IS : 1 + iLi*cGrid.iD0IS + (cGrid.iD0TL-1)*cGrid.iD0TS: 2*cGrid.iD0TS) = &
			dimag(cGrid.pacD1(1 + iLi*cGrid.iD1IS : 1 + iLi*cGrid.iD1IS + (cGrid.iD0TL/2-1)*cGrid.iD1TS : cGrid.iD1TS))
		end do
	else
		do iLi = 0, cGrid.iD0LocN-1
			cGrid.parD0(1 + iLi*cGrid.iD0IS : 1 + iLi*cGrid.iD0IS + (cGrid.iD0TL-1)*cGrid.iD0TS: 2*cGrid.iD0TS) = &
			dreal(cGrid.pacD1(1 + iLi*cGrid.iD1IS : 1 + iLi*cGrid.iD1IS + (cGrid.iD0TL/2-1)*cGrid.iD1TS : cGrid.iD1TS))
			cGrid.parD0(2 + iLi*cGrid.iD0IS : 1 + iLi*cGrid.iD0IS + (cGrid.iD0TL-1)*cGrid.iD0TS: 2*cGrid.iD0TS) = &
			dimag(cGrid.pacD1(1 + iLi*cGrid.iD1IS : 1 + iLi*cGrid.iD1IS + (cGrid.iD0TL/2-1)*cGrid.iD1TS : cGrid.iD1TS))
		end do
	end if	
	
    call GridDistr1DeAllocate(cGrid);

	! Set the flag to 0
	cGrid.iDistr		= 0;

	call SWStop(cswBlock);	call SWStop(cswBlockRedist1_0)
END SUBROUTINE ReorderDistr1ToDistr0

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %
! % 		Distr2ObtainXYZBlock
! %
! %	This functions puts the value of the XYZ-block with index iIndex from distribution 2 in 
! %	Cgrid, into the array pacBlock. pacBlock is a 3-dimensional array with dimensions:
! %	(2*cGrid.iD2XL, 2*cGrid.iD2YL, 2*cGrid.iD2ZL);
! %
! %	input:
! %		cGrid:		a grid, stored in distribution 2.
! %		iIndex:		an index of the xyz-block we will obtain
! %		pacBlock:	an array of complex values, with dimensions 
! %				(2*cGrid.iD2XL, 2*cGrid.iD2YL, 2*cGrid.iD2ZL)
! %	output:
! %		pacBlock:	an array containing the values of the grid, with  dimensions 
! %				(2*cGrid.iD2XL, 2*cGrid.iD2YL, 2*cGrid.iD2ZL) and strides
! %				(cGrid.iD2XS, cGrid.iD2YS, cGrid.iD2ZS) 
! %
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE Distr2ObtainXYZBlock(cGridS, cGridD, iOmega)
	
! =============================================================================
!
!   Programmer: Jasper de Koning / Koos Huijssen
!
!   Language: Fortran 90
!
!   Version Date    Comment
!   ------- -----   -------
!   1.0     090505  Original code (KH)
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The subroutine Distr2ObtainXYZBlock extracts an XYZ block at a specific
!   frequency from a grid which has its data stored in Distr. 2. The XYZ block
!   is stored in the first frequency in cGridD.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cGridS   i   type(Grid)   The grid that contains the source data array
!   cGridD   io  type(Grid)   The grid that will contain the XYZ block at exit
!   iOmega   i    i8b         The desired frequency of the XYZ block
!
	type(Grid), intent(in)::		cGridS
	type(Grid), intent(inout)::		cGridD
	integer(i8b), intent(in) ::		iOmega;

! *****************************************************************************
! 
!   LOCAL PARAMETERS      
!
!   iL       i8b   Loop counter over all XYZ-positions
!   iLX      i8b   Loop counter over all X-positions
!   iLY      i8b   Loop counter over all Y-positions
!   iLZ      i8b   Loop counter over all Z-positions
!   iIndS    i8b   Index to a position in the source array
!   iIndD    i8b   Index to a position in the destination array
!   acTemp   char  Temporary char array for output log messages
!
	integer(i8b) ::				iL, iLX, iLY, iLZ, iIndS, iIndD
	character(LEN = 2048) ::		acTemp;
	
! *****************************************************************************
!
!   I/O
!
!   Log file entries 
!   
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   SWstartandcount
!   SWstop
!   PrintToLog
!   
! =============================================================================

	write (acTemp, '("Distr2ObtainXYZBlock")');	call PrintToLog(acTemp, 4);
!	print *, " Punto di Interesse Odierno 0 New"
	call SwStartAndCount(cswBlock); call SwStartAndCount(cswBlockObtain);
!    print *, " Punto di Interesse Odierno 1 New"
	! First we initialize the block to zero
	do iL = 0, cGridD%ID2XL*cGridD%ID2YL*cGridD%ID2ZL-1
		cGridD%pacD2(1+iL)	= 0.0_dp;
	end do
!	print *, " Punto di Interesse Odierno 2 New"
	! Copy cGridS to cGridD on the positions starting from X/Y/Z beam index (0,0,0)
!	print *, "size(cGridD%pacD2)"
!	print *, size(cGridD%pacD2)
!	print *, "size(cGridS%pacD2)"
!	print *, size(cGridS%pacD2)
	do iLX = 0, cGridS%iD2XL-1
		do iLY = 0, cGridS%iD2YL-1
			do iLZ = 0, cGridS%iD2ZL-1
				iIndS= 1 + iOmega * cGridS%iD2IS + iLX * cGridS%iD2XS + iLY * cGridS%iD2YS + iLZ * cGridS%iD2ZS
				iIndD= 1 + iLX * cGridD%iD2XS + iLY * cGridD%iD2YS + iLZ * cGridD%iD2ZS
				!print *, "iIndD"
				!print *, iIndD
				!print *, "iIndS"
				!print *, iIndS
				
				cGridD%pacD2(iIndD) = cGridS%pacD2(iIndS);
			end do
		end do
	end do
!	print *, " Punto di Interesse Odierno 3 New"
	call SWStop(cswBlock); call SWStop(cswBlockObtain);

END SUBROUTINE Distr2ObtainXYZBlock

SUBROUTINE Distr2ObtainYMirroredXYZBlock(cGridS, cGridD, iOmega)

! =============================================================================
!
!   Programmer: Jasper de Koning / Koos Huijssen
!
!   Language: Fortran 90
!
!   Version Date    Comment
!   ------- -----   -------
!   1.0     090505  Original code (KH)
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The subroutine Distr2ObtainXYZBlock extracts an XYZ block at a specific
!   frequency from a grid which has its data stored in Distr. 2. The XYZ block
!   is stored in the first frequency in cGridD, and it is mirrored in the
!   Y-dimension. The line Y=0 is also excluded from the copy. This subroutine
!   is used to represent the contrast sources at the negative y-axis for a 
!   problem that is symmetric in Y, and it is assumed that cGridS only stores
!   values on the positive Y-axis starting at Y=0.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cGridS   i   type(Grid)   The grid that contains the source data array
!   cGridD   io  type(Grid)   The grid that will contain the XYZ block at exit
!   iOmega   i    i8b         The frequency of the XYZ block
!
	type(Grid), intent(in)::		cGridS
	type(Grid), intent(inout)::		cGridD
	integer(i8b), intent(in) ::		iOmega;

! *****************************************************************************
! 
!   LOCAL PARAMETERS      
!
!   iL       i8b   Loop counter over all XYZ-positions
!   iLX      i8b   Loop counter over all X-positions
!   iLY      i8b   Loop counter over all Y-positions
!   iLZ      i8b   Loop counter over all Z-positions
!   iIndS    i8b   Index to a position in the source array
!   iIndD    i8b   Index to a position in the destination array
!   acTemp   char  Temporary char array for output log messages
!
	integer(i8b) ::				iL, iLX, iLY, iLZ, iIndS, iIndD
	character(LEN = 2048) ::		acTemp;
	
! *****************************************************************************
!
!   I/O
!
!   Log file entries 
!   
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   SWstartandcount
!   SWstop
!   PrintToLog
!   
! =============================================================================

	write (acTemp, '("Distr2ObtainYMirroredXYZBlock")');	call PrintToLog(acTemp, 4);
	
	call SwStartAndCount(cswBlock); call SwStartAndCount(cswBlockObtain);
	! First we initialize the block to zero
	do iL = 0, cGridD%ID2XL*cGridD%ID2YL*cGridD%ID2ZL-1
		cGridD%pacD2(1+iL)	= 0;
	end do
	
	! Copy cSpaceS to cSpaceD on the positions starting from X/Y/Z beam index (0,0,0)
	! But mirror the source indices in the Y dimension
	! and exclude y=0 from the copy (since we have already calculated the influence
	! of the contrast source at y=0 in the un-mirrored block...
	do iLX = 0, cGridS%iD2XL-1
		do iLY = 0, cGridS%iD2YL-2
			do iLZ = 0, cGridS%iD2ZL-1
				iIndS= 1 + iOmega * cGridS%iD2IS + iLX * cGridS%iD2XS + (cGridS%iD2YL-1-iLY) * cGridS%iD2YS + iLZ * cGridS%iD2ZS
				iIndD= 1 + iLX * cGridD%iD2XS + iLY * cGridD%iD2YS + iLZ * cGridD%iD2ZS
				cGridD%pacD2(iIndD) = cGridS%pacD2(iIndS);
			end do
		end do
	end do
	
	call SWStop(cswBlock); call SWStop(cswBlockObtain);

END SUBROUTINE Distr2ObtainYMirroredXYZBlock

SUBROUTINE Distr2PutXYZBlock(cGridS, cGridD, iOmega)
	
! =============================================================================
!
!   Programmer: Jasper de Koning / Koos Huijssen
!
!   Language: Fortran 90
!
!   Version Date    Comment
!   ------- -----   -------
!   1.0     090505  Original code (KH)
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The subroutine Distr2ObtainXYZBlock puts an XYZ block at a specific
!   frequency in a Distr. 2 grid, i.e. it performs the inverse operation of 
!   Distr2ObtainXYZBlock.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cGridS   i   type(Grid)   The grid that contains the XYZ block
!   cGridD   io  type(Grid)   The grid that from which the XYZ block at
!                             frequency iOmega is updated at exit
!   iOmega   i    i8b         The frequency of the XYZ block
!
	type(Grid), intent(in)::		cGridS
	type(Grid), intent(inout)::		cGridD
	integer(i8b), intent(in) ::		iOmega;

! *****************************************************************************
! 
!   LOCAL PARAMETERS      
!
!   iL       i8b   Loop counter over all XYZ-positions
!   iLX      i8b   Loop counter over all X-positions
!   iLY      i8b   Loop counter over all Y-positions
!   iLZ      i8b   Loop counter over all Z-positions
!   iIndS    i8b   Index to a position in the source array
!   iIndD    i8b   Index to a position in the destination array
!   acTemp   char  Temporary char array for output log messages
!
	integer(i8b) ::				iL, iLX, iLY, iLZ, iIndS, iIndD
	character(LEN = 2048) ::		acTemp;
	
! *****************************************************************************
!
!   I/O
!
!   Log file entries 
!   
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   SWstartandcount
!   SWstop
!   PrintToLog
!   
! =============================================================================

	write (acTemp, '("Distr2PutXYZBlock")');	call PrintToLog(acTemp, 4);

	call SwStartAndCount(cswBlock);	call SwStartAndCount(cswBlockPut);
	! Copy cSpaceS starting from X/Y/Z beam index (0,0,0) to cSpaceD
	do iLX = 0, cGridD%iD2XL-1
		do iLY = 0, cGridD%iD2YL-1
			do iLZ = 0, cGridD%iD2ZL-1
				iIndS= 1 + iLX * cGridS%iD2XS + iLY * cGridS%iD2YS + iLZ * cGridS%iD2ZS
				iIndD= 1 + iOmega * cGridD%iD2IS + iLX * cGridD%iD2XS + iLY * cGridD%iD2YS + iLZ * cGridD%iD2ZS
				cGridD%pacD2(iIndD) = cGridS%pacD2(iIndS);
			end do
		end do
	end do
	
	call SWStop(cswBlock); call SWStop(cswBlockPut);
END SUBROUTINE Distr2PutXYZBlock

SUBROUTINE Distr2AddToXYZBlock(cGridS, cGridD, iOmega)

! =============================================================================
!
!   Programmer: Jasper de Koning / Koos Huijssen
!
!   Language: Fortran 90
!
!   Version Date    Comment
!   ------- -----   -------
!   1.0     090505  Original code (KH)
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The subroutine Distr2AddToXYZBlock adds an XYZ block to the values already
!   present at a specific frequency in a Distr. 2 grid.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cGridS   i   type(Grid)   The grid that contains the XYZ block
!   cGridD   io  type(Grid)   The grid that from which the XYZ block at
!                             frequency iOmega is updated at exit
!   iOmega   i    i8b         The frequency of the XYZ block
!
	type(Grid), intent(in)::		cGridS
	type(Grid), intent(inout)::		cGridD
	integer(i8b), intent(in) ::		iOmega;

! *****************************************************************************
! 
!   LOCAL PARAMETERS      
!
!   iL       i8b   Loop counter over all XYZ-positions
!   iLX      i8b   Loop counter over all X-positions
!   iLY      i8b   Loop counter over all Y-positions
!   iLZ      i8b   Loop counter over all Z-positions
!   iIndS    i8b   Index to a position in the source array
!   iIndD    i8b   Index to a position in the destination array
!   acTemp   char  Temporary char array for output log messages
!
	integer(i8b) ::				iL, iLX, iLY, iLZ, iIndS, iIndD
	character(LEN = 2048) ::		acTemp;
	
! *****************************************************************************
!
!   I/O
!
!   Log file entries 
!   
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   SWstartandcount
!   SWstop
!   PrintToLog
!   
! =============================================================================

	write (acTemp, '("Distr2AddToXYZBlock")');	call PrintToLog(acTemp, 4);

	call SwStartAndCount(cswBlock);	call SwStartAndCount(cswBlockPut);
	! Copy cSpaceS starting from X/Y/Z beam index (0,0,0) to cSpaceD
	do iLX = 0, cGridD%iD2XL-1
		do iLY = 0, cGridD%iD2YL-1
			do iLZ = 0, cGridD%iD2ZL-1
				iIndS= 1 + iLX * cGridS%iD2XS + iLY * cGridS%iD2YS + iLZ * cGridS%iD2ZS
				iIndD= 1 + iOmega * cGridD%iD2IS + iLX * cGridD%iD2XS + iLY * cGridD%iD2YS + iLZ * cGridD%iD2ZS
				cGridD%pacD2(iIndD) = cGridD%pacD2(iIndD) + cGridS%pacD2(iIndS);
			end do
		end do
	end do
	
	call SWStop(cswBlock); call SWStop(cswBlockPut);

END SUBROUTINE Distr2AddToXYZBlock

END MODULE ParnacDataRedist
