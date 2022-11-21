MODULE ParnacOutput

! =============================================================================
!
!   Programmer: Koos Huijssen / Jasper de Koning
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
!   The module ParnacOutput contains subroutines that facilitate an easy
!   export of data from Parnac to output files. The files are in the HDF5
!   format, which is a generic data format which becomes more and more popular.
!
!   more info: www.hdfgroup.org/hdf5
!
! *****************************************************************************
!
!   MODULES USED
!
    USE Types
    USE ParnacGeneral
    USE HDF5
    USE ParnacDataDef
    USE ParnacParamDef

    USE ParnacDataSpaceInit, ONLY: &
        iBeamOffsetT, iBeamOffsetX, iBeamOffsetY

! *****************************************************************************
!
!   GLOBAL DECLARATIONS
!
!   Exportarray   Generic interface to all specific ExportXDArray and
!                 ExportXDComplexArray subroutines
!
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: ExportSlice, ExportArray, Print_Error

    INTERFACE ExportArray
        MODULE PROCEDURE Export1DArray
        MODULE PROCEDURE Export1DComplexArray
        MODULE PROCEDURE Export2DArray
        MODULE PROCEDURE Export2DComplexArray
        MODULE PROCEDURE Export3DArray
        MODULE PROCEDURE Export3DComplexArray
    END INTERFACE

! *****************************************************************************
!
!   CONTAINED SUBROUTINES AND FUNCTIONS
!
!   ExportSlice        sub   Exports a subset from a data space to HDF5 file
!   Export1DArray      sub   Exports a 1D real array to HDF5 file
!   Export1DComplexArray   sub   Exports a 1D complex array to HDF5 file
!   Expord2DArray      sub   Exports a 1D real array to HDF5 file
!   Export2DComplexArray   sub   Exports a 1D complex array to HDF5 file
!   Export3DArray      sub   Exports a 1D real array to HDF5 file
!   Export3DComplexArray   sub   Exports a 1D complex array to HDF5 file
!   WriteAttributeDP   sub   Writes a dp attribute to a HDF5 file
!   WriteAttributei8b  sub   Writes an i8b attribute to a HDF5 file
!
! =============================================================================
!
CONTAINS

    SUBROUTINE ExportSlice(pszFileName, pszVarName, cSpace, aiStart, aiLength, SliceIndex, aiLengthSlice, bIsReal)

! =============================================================================
!
!   Programmer: Koos Huijssen / Jasper de Koning
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
!   The subroutine ExportSlice exports a subset from a data space to a HDF5
!   output file. This subset is generally smaller than the complete data space
!   and is defined by the ranges in all dimensions. It is always the space in
!   terms of the beam indices that is saved. The distribution of the space
!   may be 0, 1 or 2, and the data in it may be real or complex.
!
!   Each processor writes the part of the data that is contained on that
!   processor, and the complete data space has to be reconstructed from those
!   parts. The space is basically stored as a 1D block of data, where the order
!   of the elements and the striding is depending on the distribution. To
!   characterize the data space and to enable its reconstruction in the global
!   grid, the following attributes are saved along with it:
!
!   iStart:                    global index of the first position in T/X/Y/Z of the subset
!   iT/X/Yoffset: offset tables in the T,X and Y dimension with length DimZ,
!                 giving the number of grid points that needs to be shifted
!                 as a function of z to obtain the correct skewness of the
!                 beam.
!   rTheta:                    the angles of the beam in T/X/Y, in radians
!   rStepsize:          the step size in each dimension, in the order of T/X/Y/Z,
!                 in seconds and meters
!   iSavelength:  The length of each dimension; the dimensions are ordered
!                 in respecting the permutation order given in iPermute.
!                 E.g. if iPermute(1)=4, then iSavelength(4) gives the length
!                 of the temporal dimension (which is normally dimension 1)
!   iPermute:          the permutation order of the dimensions in the saved
!                 subset. The 'normal' order would be T=1, X=2, Y=3, Z=4.
!                 E.g. an iPermute of (/ 4, 1, 2, 3 /) means that the
!                 temporal dimension is the fourth dimension, X is the first,
!                 Y is the second and Z is the third. The first dimension
!                 has a stride of 1, the last dimension varies slowest.
!   iSliceIndex:  Slice index                                                   ! added by A.M
!   iIsReal:          1 if the data in the data array is real, 0 if it is complex
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   pszFileName  i   char   Filename of the export file
!   pszVarName   i   char   Name of the variable in the export file
!   cSpace       i  type(Space)   Space from which the subset needs to be taken
!   aiStart      i   i8b    Starting beam indices of the subset within the beam
!   aiLength     i   i8b    Length of the subset in each dimension
!   bIsReal      i   lgt    Whether the data in the slice has to be treated as
!                           real or complex - this is not trivial for all
!                           situations
!
        character(len=*), intent(in)::        pszFileName, pszVarName
        type(Space), intent(in)::                cSpace
        integer(i8b), intent(in)::                aiStart(4), aiLength(4), aiLengthSlice; 
        integer(i8b), intent(in)::                SliceIndex(cModelParams%Numslices)
        logical, intent(in) ::                        bIsReal

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   iErr       i4b    Error number
!   iFile     HID_T   HDF5 File ID
!   iSpaceID  HID_T   HDF5 Space ID
!   iSetID    HID_T   HDF5 Set ID
!   iLi        i8b    Index over all XYZ positions
!   iInd       i8b    Index inside the temporary buffer array
!   iIndSrc    i8b    Index inside the data space array
!   iCount     i8b    Number of elements on this processor to be exported
!   iXi        i8b    Beam index in x
!   iYi        i8b    Beam index in y
!   iZi        i8b    Beam index in z
!   acTemp     char   Temporary char array for output log messages
!   acTempFile char   Filename including directory path
!   acTempVar  char   Variable name
!   arBuffer   dp     Data buffer of the points to be exported to file
!   dX/Y/Z/T   dp     Step size in X/Y/Z/T in meters and seconds
!   aiTempStart   i8b  Starting beam index of the subset in T/X/Y/Z
!   aiTempStart2  i8b  Starting global index of the subset in T/X/Y/Z
!   aiTempLength  i8b  Adjusted length of the subset in T/X/Y/Z
!   aiSaveLength  i8b  Lengths as saved to output file: the order is not
!                      necessarily T/X/Y/Z
!   aiPermute     i8b  Permutation order of the T/X/Y/Z dimensions in the
!                      data array - e.g. if aiPermute(1)=4, then the temporal
!                      dimension is in the fourth, i.e. slowest varying
!                      dimension (with the largest stride)
!   aiTempOffset  i8b  array containing the offsets in the T/X/Y dimensions
!                      as a function of Z
!
        integer(i4b)::              iErr
        integer(HID_T)::            iFile, iSpaceID, iSetID
        integer(i8b)::                                iLi, iInd, iCount, iXi, iYi, iZi, iIndSrc
        character(len=1024)::                acTemp, acTempFile, acTempVar; 
        real(dp), allocatable::                arBuffer(:)
        real(dp)::                                        dX, dY, dZ, dT; 
        integer(i8b)::                                aiTempStart(4), aiTempStart2(4), &
                       aiTempLength(4), aiSaveLength(4), aiPermute(4), &
                       aiTempOffset(aiLength(4), 3)

! *****************************************************************************
!
!   I/O
!
!   log file entries and export to HDF5 output file
!
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   PrintToLog
!   h5open_f
!   h5fcreate_f
!   h5screate_simple_f
!   h5dwrite_f
!   h5sclose_f
!   WriteAttributei8b
!   WriteAttributeDP
!   h5dclose_f
!   h5fclose_f
!   h5close_f
!
! *****************************************************************************
!
!   PSEUDOCODE
!
!   Adjust the requested lengths of the subset in all dimensions to the
!     possible
!   If the distribution of the data space is Dist. 0
!       Count the number of elements in the data space that are in the
!         requested subset and on the current processor
!       If the number of elements is larger than zero
!          Allocate the data buffer
!           Fill the data buffer with the elements on the current processor that
!             fall within the requested subset from the data space from Dist. 0
!           Set the aiSavelength and aiPermute vectors
!   If the distribution is Dist. 1
!       Perform the same operations for Dist. 1
!   If the distribution is Dist. 2
!       Perform the same operations for Dist. 2
!   If the number of elements is larger than zero
!       Open the output file and create a HDF structure
!       Export the data buffer to the output file
!       Export the global starting indices as the attribute iStart
!       Export the beam offsets in T/X/Y as a function of Z as the attributes
!         iT/X/Yoffset
!       Export the beam angles in T/X/Y as the attribute rTheta
!       Export the step sizes in the T/X/Y/Z dimensions as the attribute
!         rStepsize
!       Export the lengths of the dimensions of the data array as the attribute
!         aiSavelength
!       Export the order of the dimensions in the data array as the attribute
!         aiPermute
!       Export the total number of processors as the attribute iNumProc
!       Export the flag whether the data in the data array is real as the
!         attribute iIsReal
!       Close the HDF structure and the output file
!       Deallocate the data buffer
!   Else
!       Don't save anything for this processor, write message to log
!
! =============================================================================
        write (acTemp, '("Exporting ", A, " to file ", A)') trim(pszVarName), trim(pszFilename)//trim(int2str(cSpace.cGrid.iProcID))//".h5"
        call PrintToLog(acTemp, 1); 
        !Possible length may be less than desired because aiStart+aiLength may be larger
        ! than iDim
        aiTempStart = aiStart
        if (cSpace.cGrid.iDistr == 0) then
            aiTempLength(1) = min(aiLength(1), cSpace.iDimT - aiTempStart(1))
        elseif (bIsReal) then
            aiTempLength(1) = min(aiLength(1), 2*(cSpace.iDimT + 1) - aiTempStart(1))
        else
            aiTempLength(1) = min(aiLength(1), cSpace.iDimT + 1 - aiTempStart(1))
        end if
        aiTempLength(2) = min(aiLength(2), cSpace.iDimX - aiTempStart(2))
        aiTempLength(3) = min(aiLength(3), cSpace.iDimY - aiTempStart(3))
        aiTempLength(4) = min(aiLength(4), cSpace.iDimZ - aiTempStart(4))
!        open (8,  file=trim(trim(sOutputDir)//'Bubble_indexes_output'//int2str(cSpace%cGrid%iProcID)),  position='append')
!        write(8,*) iLi
!        close(8,status='delete')
        !Handle different distributions
        ! In Distr. 1 and 2, data can be real as well as complex; in distr. 0, it can only be real
        if (cSpace.cGrid.iDistr == 0) then
            ! First we count how many points in the slice are stored on this processor
            iCount = 0; 
            do iLi = 0, cSpace.cGrid.iD0LocN - 1
                if (cSpace.cGrid.aiD0Loc(1 + iLi, 1) >= aiTempStart(2) .and. cSpace.cGrid.aiD0Loc(1 + iLi, 1) < aiTempStart(2) + aiTempLength(2) .and. &
                    cSpace.cGrid.aiD0Loc(1 + iLi, 2) >= aiTempStart(3) .and. cSpace.cGrid.aiD0Loc(1 + iLi, 2) < aiTempStart(3) + aiTempLength(3) .and. &
                    cSpace.cGrid.aiD0Loc(1 + iLi, 3) >= aiTempStart(4) .and. cSpace.cGrid.aiD0Loc(1 + iLi, 3) < aiTempStart(4) + aiTempLength(4)) then
                    iCount = iCount + aiTempLength(1); 
                end if
            end do

            if (iCount > 0) then
                ! Allocate the buffer and fill it
                allocate (arBuffer(iCount)); 
                iInd = 1; 
                do iLi = 0, cSpace.cGrid.iD0LocN - 1
                    if (cSpace.cGrid.aiD0Loc(1 + iLi, 1) >= aiTempStart(2) .and. cSpace.cGrid.aiD0Loc(1 + iLi, 1) < aiTempStart(2) + aiTempLength(2) .and. &
                        cSpace.cGrid.aiD0Loc(1 + iLi, 2) >= aiTempStart(3) .and. cSpace.cGrid.aiD0Loc(1 + iLi, 2) < aiTempStart(3) + aiTempLength(3) .and. &
                        cSpace.cGrid.aiD0Loc(1 + iLi, 3) >= aiTempStart(4) .and. cSpace.cGrid.aiD0Loc(1 + iLi, 3) < aiTempStart(4) + aiTempLength(4)) then
                        arBuffer(iInd:iInd + aiTempLength(1) - 1) = &
                            cSpace.cGrid.parD0(1 + iLi*cSpace.cGrid.iD0IS + aiTempStart(1): &
                                               1 + iLi*cSpace.cGrid.iD0IS + aiTempStart(1) + aiTempLength(1) - 1); 
                        iInd = iInd + aiTempLength(1); 
                    end if
                end do

                !Set the lengths of the axes in the order of storage in the 1-D array
                !Order of storage is TXYZ
                !If another order was chosen for XYZ in InitGridDistrLookupTables, change
                ! aiSaveLength and aiPermute as well!
                aiSaveLength(1) = aiTempLength(1)
                aiSaveLength(2) = aiTempLength(2)
                aiSaveLength(3) = aiTempLength(3)
                aiSaveLength(4) = aiTempLength(4)

                !Set the permute dimensions array such as to reconfigure an array in the order (T,X,Y,Z)
                !aiPermute(1) is the dimension that has to be placed first, i.e. if X is first in the
                ! required order, then aiPermute(1) gets number 2 here...
                aiPermute(1) = 1
                aiPermute(2) = 2
                aiPermute(3) = 3
                aiPermute(4) = 4
            end if

        else if (cSpace.cGrid.iDistr == 1) then

            !since the real values are stored as though they were complex,
            ! we need to halve aiTempStart(1), and aiTempLength(1), rounded up
            if (bIsReal) then
                aiTempStart(1) = aiTempStart(1)/2
                aiTempLength(1) = (aiTempLength(1) + 1)/2
            end if

            ! First we count how many points in the slice are stored on this processor
            iCount = 0; 
            do iLi = 0, cSpace.cGrid.iD1LocN - 1
                if (cSpace.cGrid.aiD1Loc(1 + iLi, 1) >= aiTempStart(2) .and. cSpace.cGrid.aiD1Loc(1 + iLi, 1) < aiTempStart(2) + aiTempLength(2) .and. &
                    cSpace.cGrid.aiD1Loc(1 + iLi, 2) >= aiTempStart(3) .and. cSpace.cGrid.aiD1Loc(1 + iLi, 2) < aiTempStart(3) + aiTempLength(3) .and. &
                    cSpace.cGrid.aiD1Loc(1 + iLi, 3) >= aiTempStart(4) .and. cSpace.cGrid.aiD1Loc(1 + iLi, 3) < aiTempStart(4) + aiTempLength(4)) then
                    iCount = iCount + aiTempLength(1)*2; ! times 2 for complex
                end if
            end do

            if (iCount > 0) then
                ! Allocate the buffer and fill it
                allocate (arBuffer(iCount)); 
                iInd = 1; 
                do iLi = 0, cSpace.cGrid.iD1LocN - 1
                    if (cSpace.cGrid.aiD1Loc(1 + iLi, 1) >= aiTempStart(2) .and. cSpace.cGrid.aiD1Loc(1 + iLi, 1) < aiTempStart(2) + aiTempLength(2) .and. &
                        cSpace.cGrid.aiD1Loc(1 + iLi, 2) >= aiTempStart(3) .and. cSpace.cGrid.aiD1Loc(1 + iLi, 2) < aiTempStart(3) + aiTempLength(3) .and. &
                        cSpace.cGrid.aiD1Loc(1 + iLi, 3) >= aiTempStart(4) .and. cSpace.cGrid.aiD1Loc(1 + iLi, 3) < aiTempStart(4) + aiTempLength(4)) then
                        arBuffer(iInd:(iInd + 2*(aiTempLength(1) - 1)):2) = &
                            REAL(cSpace.cGrid.pacD1(1 + iLi*cSpace.cGrid.iD1IS + aiTempStart(1): &
                                                    1 + iLi*cSpace.cGrid.iD1IS + aiTempStart(1) + aiTempLength(1) - 1), dp); 
                        arBuffer((iInd + 1):(iInd + 1 + 2*(aiTempLength(1) - 1)):2) = &
                            DIMAG(cSpace.cGrid.pacD1(1 + iLi*cSpace.cGrid.iD1IS + aiTempStart(1): &
                                                     1 + iLi*cSpace.cGrid.iD1IS + aiTempStart(1) + aiTempLength(1) - 1)); 
                        iInd = iInd + aiTempLength(1)*2; ! times 2 for complex
                    end if
                end do

                !Reset aiTempStart(1) and aiTempLength(1) for the file if necessary
                if (bIsReal) then
                    aiTempStart(1) = aiTempStart(1)*2
                    aiTempLength(1) = aiTempLength(1)*2
                end if

                !Set the lengths of the axes in the order of storage in the 1-D array
                !Order of storage is TXYZ
                !If another order was chosen for XYZ in InitGridDistrLookupTables, change
                ! aiSaveLength and aiPermute as well!
                aiSaveLength(1) = aiTempLength(1)
                aiSaveLength(2) = aiTempLength(2)
                aiSaveLength(3) = aiTempLength(3)
                aiSaveLength(4) = aiTempLength(4)

                !Set the permute dimensions array such as to reconfigure an array in the order (T,X,Y,Z)
                !aiPermute(1) is the dimension that has to be placed first, i.e. if X is first in the
                ! required order, then aiPermute(1) gets number 2 here...
                aiPermute(1) = 1
                aiPermute(2) = 2
                aiPermute(3) = 3
                aiPermute(4) = 4
            end if

        else if (cSpace.cGrid.iDistr == 2) then

            !since the real values are stored as though they were complex,
            ! we need to halve the start and the length, rounded up
            if (bIsReal) then
                aiTempStart(1) = aiTempStart(1)/2
                aiTempLength(1) = (aiTempLength(1) + 1)/2
            end if

            ! First we count how many points in the slice are stored on this processor
            iCount = 0; 
            do iLi = 0, cSpace.cGrid.iD2LocN - 1
                if (cSpace.cGrid.aiD2Loc(1 + iLi) >= aiTempStart(1) .and. cSpace.cGrid.aiD2Loc(1 + iLi) < aiTempStart(1) + aiTempLength(1)) then
                    iCount = iCount + aiTempLength(2)*aiTempLength(3)*aiTempLength(4)*2; ! times 2 for complex, and for two reals
                end if
            end do

            if (iCount > 0) then
                ! Allocate the buffer and fill it
                allocate (arBuffer(iCount)); 
                iInd = 1; 
                do iLi = 0, cSpace.cGrid.iD2LocN - 1
                    if (cSpace.cGrid.aiD2Loc(1 + iLi) >= aiTempStart(1) .and. cSpace.cGrid.aiD2Loc(1 + iLi) < aiTempStart(1) + aiTempLength(1)) then

                        !KH We have the order i = x + y * DimX + z*DimX*DimY + omega * DimX*DimY*DimZ
                        !If another order was chosen in InitGridDistrLookupTables, change the order here!!
                        do iZi = aiTempStart(4), aiTempStart(4) + aiTempLength(4) - 1
                            do iYi = aiTempStart(3), aiTempStart(3) + aiTempLength(3) - 1

                                iIndSrc = 1 + iLi*cSpace.cGrid.iD2IS + iZi*cSpace.cGrid.iD2ZS + &
                                          iYi*cSpace.cGrid.iD2YS + aiTempStart(2)

                                !take care! If bIsReal, then values are stored rather peculiar in Dist. 2,
                                ! namely in couples t0,t1 for all XYZ traces as though they are real and imag parts,
                                ! then t2,t3 for all XYZ traces etc...
                                if (bIsReal) then

                                    arBuffer(iInd:(iInd + aiTempLength(2) - 1)) = &
                                        REAL(cSpace.cGrid.pacD2(iIndSrc:(iIndSrc + aiTempLength(2) - 1)), dp); 
                                    arBuffer((iInd + aiTempLength(2)*aiTempLength(3)*aiTempLength(4)): &
                                             (iInd + aiTempLength(2)*aiTempLength(3)*aiTempLength(4) + aiTempLength(2) - 1)) = &
                                        DIMAG(cSpace.cGrid.pacD2(iIndSrc:iIndSrc + aiTempLength(2) - 1)); 
                                    iInd = iInd + aiTempLength(2); 
                                else

                                    arBuffer(iInd:(iInd + 2*(aiTempLength(2) - 1)):2) = &
                                        REAL(cSpace.cGrid.pacD2(iIndSrc:iIndSrc + aiTempLength(2) - 1), dp); 
                                    arBuffer((iInd + 1):(iInd + 1 + 2*(aiTempLength(2) - 1)):2) = &
                                        DIMAG(cSpace.cGrid.pacD2(iIndSrc:iIndSrc + aiTempLength(2) - 1)); 
                                    iInd = iInd + aiTempLength(2)*2; ! times 2 for complex

                                end if
                            end do
                        end do

                        !after the XYZ block, correct for the second block stored in the loop...
                        if (bIsReal) then

                            iInd = iInd + aiTempLength(2)*aiTempLength(3)*aiTempLength(4); 
                        end if

                    end if
                end do

                !Reset aiTempStart(1) and aiTempLength(1) for the file if necessary
                if (bIsReal) then
                    aiTempStart(1) = aiTempStart(1)*2
                    aiTempLength(1) = aiTempLength(1)*2
                end if

                !Set the lengths of the axes in the order of storage in the 1-D array
                !Order of storage is XYZT
                !If another order was chosen for XYZ in InitGridDistrLookupTables, change
                ! aiSaveLength and aiPermute as well!
                aiSaveLength(1) = aiTempLength(2)
                aiSaveLength(2) = aiTempLength(3)
                aiSaveLength(3) = aiTempLength(4)
                aiSaveLength(4) = aiTempLength(1)

                !Set the permute dimensions array such as to reconfigure an array in the order (T,X,Y,Z)
                !aiPermute(1) is the dimension that has to be placed first, i.e. if X is first in the
                ! required order, then aiPermute(1) gets number 2 here...
                aiPermute(1) = 4
                aiPermute(2) = 1
                aiPermute(3) = 2
                aiPermute(4) = 3
            end if

        end if

        if (iCount > 0) then

            call PrintToLog("Start writing file", 3); 
            ! Open the file and initialize the necessary data
            acTempFile = trim(pszFileName)//trim(int2str(cSpace.cGrid.iProcID))//".h5"; 
            acTempVar = trim(pszVarName); 
            call h5open_f(iErr); 
            call h5fcreate_f(trim(acTempFile), H5F_ACC_TRUNC_F, iFile, iErr)
            ! Export the slice
            ! Export the data
            call h5screate_simple_f(1, (/int(iCount, HSIZE_T)/), iSpaceID, iErr)
            call h5dcreate_f(iFile, trim(acTempVar), H5T_NATIVE_DOUBLE, iSpaceID, iSetID, iErr); 
            call h5dwrite_f(iSetID, H5T_NATIVE_DOUBLE, arBuffer, (/int(iCount, HSIZE_T)/), iErr)
            call h5sclose_f(iSpaceID, iErr)

            !         Export attributes: Integer start position and offsets, and real angles and step sizes
            !        Also export dimensions, permutation order, number of processors and the IsReal argument for reconstruction of the array
            !Store global start index
            aiTempStart2(1) = cSpace.iStartT + iBeamOffsetT(cSpace, aiTempStart(1)); 
            aiTempStart2(2) = cSpace.iStartX + iBeamOffsetX(cSpace, aiTempStart(1)); 
            aiTempStart2(3) = cSpace.iStartY + iBeamOffsetY(cSpace, aiTempStart(1)); 
!                aiTempStart2(4) = aiTempStart(4) + cSpace.iStartZ;  ! If something is wrong add the first term  Added By A.M
            aiTempStart2(4) = cSpace.iStartZ; 
            call WriteAttributei8b(iFile, iSetID, "iStart", 4_i8b, aiTempStart2)
            !Store offset tables
            do iLi = 0, aiTempLength(4) - 1
                !In saving a slice, account for the possibility that iBeamStartZ/=0, since
                !it might be an iSI_GREENSFUNCTION...
                iZi = mod(iLi + aiTempStart(4) - cSpace%iBeamIndexStartZ, cSpace%iDimZ) + cSpace%iBeamIndexStartZ
                aiTempOffset(1 + iLi, 1) = iBeamOffsetT(cSpace, iZi)
                aiTempOffset(1 + iLi, 2) = iBeamOffsetX(cSpace, iZi)
                aiTempOffset(1 + iLi, 3) = iBeamOffsetY(cSpace, iZi)
            end do
            call WriteAttributei8b(iFile, iSetID, "iToffset", aiTempLength(4), aiTempOffset(:, 1)); 
            call WriteAttributei8b(iFile, iSetID, "iXoffset", aiTempLength(4), aiTempOffset(:, 2)); 
            call WriteAttributei8b(iFile, iSetID, "iYoffset", aiTempLength(4), aiTempOffset(:, 3)); 
            call WriteAttributeDP(iFile, iSetID, "rTheta", 3_i8b, atan(1/(/cSpace.dTanT, cSpace.dTanX, cSpace.dTanY/)))
            dX = cSpace.dDx*cMediumParams%c0/cModelParams%freq0; 
            dT = cSpace.dDt/cModelParams%freq0; 
            ! STATIC VALUES
!        dX                = 110e-6
!        dT                = 0.074224021e-6

            call WriteAttributeDP(iFile, iSetID, "rStepsize", 4_i8b, (/dT, dX, dX, dX/))
            call WriteAttributei8b(iFile, iSetID, "iSaveLength", 5_i8b, (/aiSaveLength, aiLengthSlice/)) ! Added By Agis
            call WriteAttributei8b(iFile, iSetID, "iPermute", 4_i8b, aiPermute)
            iLi = cSpace.cGrid.iProcN; 
            call WriteAttributei8b(iFile, iSetID, "iNumProc", 1_i8b, (/iLi/))
            if (bIsReal .or. (cSpace.cGrid.iDistr == 0)) then
                iLi = 1
            else
                iLi = 0
            end if
            call WriteAttributei8b(iFile, iSetID, "iSliceIndex", size(SliceIndex, 1)*1_i8b, (/SliceIndex/))  ! Added By Agis
            call WriteAttributei8b(iFile, iSetID, "iIsReal", 1_i8b, (/iLi/))

            ! Deinitialize
            call h5dclose_f(iSetID, iErr)
            call h5fclose_f(iFile, iErr); 
            call h5close_f(iErr)

            call PrintToLog("End writing file", 3); 
            deallocate (arBuffer); 
        else

            write (acTemp, '("No export; requested slice had no data on this processor.")')
            call PrintToLog(acTemp, 3); 
        end if

    END SUBROUTINE ExportSlice

    SUBROUTINE Export1DArray(pszFileName, pszVarName, variable, iProcID)

! =============================================================================
!
!   Programmer: Koos Huijssen / Jasper de Koning
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
!   The subroutine Export1DArray exports a 1D real array to a HDF5 output file.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   pszFileName  i   char   Filename of the export file
!   pszVarName   i   char   Name of the variable in the export file
!   variable     i    dp    Variable to be stored
!   iProcID      i   i8b    ID of the current processor
!
        character(len=*), intent(in)::                pszFileName, pszVarName
        real(dp), intent(in) ::                                variable(:)
        integer(i8b), intent(in) ::                        iProcID

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   iErr       i4b    Error number
!   iFile     HID_T   HDF5 File ID
!   iSpaceID  HID_T   HDF5 Space ID
!   iSetID    HID_T   HDF5 Set ID
!   iDims    HSIZE_T  Length of the dimension(s) of the array
!   acTemp     char   Temporary char array for output log messages
!   acTempFile char   Filename including directory path
!   acTempVar  char   Variable name
!
        integer(i4b)::              iErr
        integer(HID_T)::            iFile, iSpaceID, iSetID
        integer(HSIZE_T) ::                        iDims(1)
        character(len=1024)::                acTemp, acTempFile, acTempVar; 
! *****************************************************************************
!
!   I/O
!
!   log file entries and export to HDF5 output file
!
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   PrintToLog
!   h5open_f
!   h5fcreate_f
!   h5screate_simple_f
!   h5dcreate_f
!   h5dwrite_f
!   h5sclose_f
!   h5dclose_f
!   h5fclose_f
!   h5close_f
!
! =============================================================================

        write (acTemp, '("Exporting ", A, " to file ", A)') trim(pszVarName), trim(pszFilename)//trim(int2str(iProcID))//".h5"
        call PrintToLog(acTemp, 1); 
        ! Allocate the buffer and fill it
        iDims = (/size(variable, 1)/); 
        ! Open the file and initialize the necessary data
        acTempFile = trim(pszFileName)//trim(int2str(iProcID))//".h5"; 
        acTempVar = trim(pszVarName); 
        call h5open_f(iErr); 
        CALL h5fcreate_f(trim(acTempFile), H5F_ACC_TRUNC_F, iFile, iErr)

        ! Export the slice
        !         Export the data
        call h5screate_simple_f(1, iDims, iSpaceID, iErr)
        call h5dcreate_f(iFile, trim(acTempVar), H5T_NATIVE_DOUBLE, iSpaceID, iSetID, iErr); 
        call h5dwrite_f(iSetID, H5T_NATIVE_DOUBLE, variable, iDims, iErr)
        call h5sclose_f(iSpaceID, iErr)

        ! Deinitialize
        call h5dclose_f(iSetID, iErr)
        call h5fclose_f(iFile, iErr); 
        call h5close_f(iErr)

    END SUBROUTINE Export1DArray

    SUBROUTINE Export1DComplexArray(pszFileName, pszVarName, variable, iProcID)

! =============================================================================
!
!   Programmer: Koos Huijssen / Jasper de Koning
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
!   The subroutine Export1DComplexArray exports a 1D complex array to a HDF5
!   output file.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   pszFileName  i   char   Filename of the export file
!   pszVarName   i   char   Name of the variable in the export file
!   variable     i   dpc    Variable to be stored
!   iProcID      i   i8b    ID of the current processor
!
        character(len=*), intent(in)::                pszFileName, pszVarName
        complex(dpc), intent(in) ::                        variable(:)
        integer(i8b), intent(in) ::                        iProcID

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   iErr       i4b    Error number
!   iFile     HID_T   HDF5 File ID
!   iSpaceID  HID_T   HDF5 Space ID
!   iSetID    HID_T   HDF5 Set ID
!   iDims    HSIZE_T  Length of the dimension(s) of the array
!   acTemp     char   Temporary char array for output log messages
!   acTempFile char   Filename including directory path
!   acTempVar  char   Variable name
!
        integer(i4b)::              iErr
        integer(HID_T)::            iFile, iSpaceID, iSetID
        integer(HSIZE_T) ::                        iDims(1)
        character(len=1024)::                acTemp, acTempFile, acTempVar; 
! *****************************************************************************
!
!   I/O
!
!   log file entries and export to HDF5 output file
!
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   PrintToLog
!   h5open_f
!   h5fcreate_f
!   h5screate_simple_f
!   h5dcreate_f
!   h5dwrite_f
!   h5sclose_f
!   h5dclose_f
!   h5fclose_f
!   h5close_f
!
! =============================================================================

        write (acTemp, '("Exporting ", A, " to file ", A)') trim(pszVarName), trim(pszFilename)//trim(int2str(iProcID))//".h5"
        call PrintToLog(acTemp, 1); 
        ! Allocate the buffer and fill it
        iDims = (/size(variable, 1)/); 
        ! Open the file and initialize the necessary data
        acTempFile = trim(pszFileName)//trim(int2str(iProcID))//".h5"; 
        acTempVar = trim(pszVarName); 
        call h5open_f(iErr); 
        CALL h5fcreate_f(trim(acTempFile), H5F_ACC_TRUNC_F, iFile, iErr)

        ! Export the slice
        !         Export the data
        call h5screate_simple_f(1, iDims, iSpaceID, iErr)
        call h5dcreate_f(iFile, trim(acTempVar), H5T_NATIVE_DOUBLE, iSpaceID, iSetID, iErr); 
        call h5dwrite_f(iSetID, H5T_NATIVE_DOUBLE, real(variable, dp), iDims, iErr)
        call h5dcreate_f(iFile, trim(acTempVar//"_cplx"), H5T_NATIVE_DOUBLE, iSpaceID, iSetID, iErr); 
        call h5dwrite_f(iSetID, H5T_NATIVE_DOUBLE, dimag(variable), iDims, iErr)
        call h5sclose_f(iSpaceID, iErr)

        ! Deinitialize
        call h5dclose_f(iSetID, iErr)
        call h5fclose_f(iFile, iErr); 
        call h5close_f(iErr)

    END SUBROUTINE Export1DComplexArray

    SUBROUTINE Export2DArray(pszFileName, pszVarName, variable, iProcID)

! =============================================================================
!
!   Programmer: Koos Huijssen / Jasper de Koning
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
!   The subroutine Export2DArray exports a 2D real array to a HDF5 output file.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   pszFileName  i   char   Filename of the export file
!   pszVarName   i   char   Name of the variable in the export file
!   variable     i    dp    Variable to be stored
!   iProcID      i   i8b    ID of the current processor
!
        character(len=*), intent(in)::                pszFileName, pszVarName
        real(dp), intent(in) ::                                variable(:, :)
        integer(i8b), intent(in) ::                        iProcID

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   iErr       i4b    Error number
!   iFile     HID_T   HDF5 File ID
!   iSpaceID  HID_T   HDF5 Space ID
!   iSetID    HID_T   HDF5 Set ID
!   iDims    HSIZE_T  Length of the dimension(s) of the array
!   acTemp     char   Temporary char array for output log messages
!   acTempFile char   Filename including directory path
!   acTempVar  char   Variable name
!
        integer(i4b)::              iErr
        integer(HID_T)::            iFile, iSpaceID, iSetID
        integer(HSIZE_T) ::                        iDims(2)
        character(len=1024)::                acTemp, acTempFile, acTempVar; 
! *****************************************************************************
!
!   I/O
!
!   log file entries and export to HDF5 output file
!
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   PrintToLog
!   h5open_f
!   h5fcreate_f
!   h5screate_simple_f
!   h5dcreate_f
!   h5dwrite_f
!   h5sclose_f
!   h5dclose_f
!   h5fclose_f
!   h5close_f
!
! =============================================================================

        write (acTemp, '("Exporting ", A, " to file ", A)') trim(pszVarName), trim(pszFilename)//trim(int2str(iProcID))//".h5"
        call PrintToLog(acTemp, 1); 
        ! Allocate the buffer and fill it
        iDims = (/size(variable, 1), size(variable, 2)/); 
        ! Open the file and initialize the necessary data
        acTempFile = trim(pszFileName)//trim(int2str(iProcID))//".h5"; 
        acTempVar = trim(pszVarName); 
        call h5open_f(iErr); 
        CALL h5fcreate_f(trim(acTempFile), H5F_ACC_TRUNC_F, iFile, iErr)

        ! Export the slice
        !         Export the data
        call h5screate_simple_f(2, iDims, iSpaceID, iErr)
        call h5dcreate_f(iFile, trim(acTempVar), H5T_NATIVE_DOUBLE, iSpaceID, iSetID, iErr); 
        call h5dwrite_f(iSetID, H5T_NATIVE_DOUBLE, variable, iDims, iErr)
        call h5sclose_f(iSpaceID, iErr)
        ! Deinitialize
        call h5dclose_f(iSetID, iErr)
        call h5fclose_f(iFile, iErr); 
        call h5close_f(iErr)

    END SUBROUTINE Export2DArray

    SUBROUTINE Export2DComplexArray(pszFileName, pszVarName, variable, iProcID)

! =============================================================================
!
!   Programmer: Koos Huijssen / Jasper de Koning
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
!   The subroutine Export2DComplexArray exports a 2D complex array to a HDF5
!   output file.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   pszFileName  i   char   Filename of the export file
!   pszVarName   i   char   Name of the variable in the export file
!   variable     i   dpc    Variable to be stored
!   iProcID      i   i8b    ID of the current processor
!
        character(len=*), intent(in)::                pszFileName, pszVarName
        complex(dpc), intent(in) ::                        variable(:, :)
        integer(i8b), intent(in) ::                        iProcID

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   iErr       i4b    Error number
!   iFile     HID_T   HDF5 File ID
!   iSpaceID  HID_T   HDF5 Space ID
!   iSetID    HID_T   HDF5 Set ID
!   iDims    HSIZE_T  Length of the dimension(s) of the array
!   acTemp     char   Temporary char array for output log messages
!   acTempFile char   Filename including directory path
!   acTempVar  char   Variable name
!
        integer(i4b)::              iErr
        integer(HID_T)::            iFile, iSpaceID, iSetID
        integer(HSIZE_T) ::                        iDims(2)
        character(len=1024)::                acTemp, acTempFile, acTempVar; 
! *****************************************************************************
!
!   I/O
!
!   log file entries and export to HDF5 output file
!
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   PrintToLog
!   h5open_f
!   h5fcreate_f
!   h5screate_simple_f
!   h5dcreate_f
!   h5dwrite_f
!   h5sclose_f
!   h5dclose_f
!   h5fclose_f
!   h5close_f
!
! =============================================================================

        write (acTemp, '("Exporting ", A, " to file ", A)') trim(pszVarName), trim(pszFilename)//trim(int2str(iProcID))//".h5"
        call PrintToLog(acTemp, 1); 
        ! Allocate the buffer and fill it
        iDims = (/size(variable, 1), size(variable, 2)/); 
        ! Open the file and initialize the necessary data
        acTempFile = trim(pszFileName)//trim(int2str(iProcID))//".h5"; 
        acTempVar = trim(pszVarName); 
        call h5open_f(iErr); 
        CALL h5fcreate_f(trim(acTempFile), H5F_ACC_TRUNC_F, iFile, iErr)

        ! Export the slice
        !         Export the data
        call h5screate_simple_f(2, iDims, iSpaceID, iErr)
        call h5dcreate_f(iFile, trim(acTempVar), H5T_NATIVE_DOUBLE, iSpaceID, iSetID, iErr); 
        call h5dwrite_f(iSetID, H5T_NATIVE_DOUBLE, real(variable, dp), iDims, iErr)
        call h5dcreate_f(iFile, trim(acTempVar//"_cplx"), H5T_NATIVE_DOUBLE, iSpaceID, iSetID, iErr); 
        call h5dwrite_f(iSetID, H5T_NATIVE_DOUBLE, dimag(variable), iDims, iErr)
        call h5sclose_f(iSpaceID, iErr)

        ! Deinitialize
        call h5dclose_f(iSetID, iErr)
        call h5fclose_f(iFile, iErr); 
        call h5close_f(iErr)

    END SUBROUTINE Export2DComplexArray

    SUBROUTINE Export3DArray(pszFileName, pszVarName, variable, iProcID)

! =============================================================================
!
!   Programmer: Koos Huijssen / Jasper de Koning
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
!   The subroutine Export3DArray exports a 3D real array to a HDF5 output file.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   pszFileName  i   char   Filename of the export file
!   pszVarName   i   char   Name of the variable in the export file
!   variable     i    dp    Variable to be stored
!   iProcID      i   i8b    ID of the current processor
!
        character(len=*), intent(in)::                pszFileName, pszVarName
        real(dp), intent(in) ::                                variable(:, :, :)
        integer(i8b), intent(in) ::                        iProcID

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   iErr       i4b    Error number
!   iFile     HID_T   HDF5 File ID
!   iSpaceID  HID_T   HDF5 Space ID
!   iSetID    HID_T   HDF5 Set ID
!   iDims    HSIZE_T  Length of the dimension(s) of the array
!   acTemp     char   Temporary char array for output log messages
!   acTempFile char   Filename including directory path
!   acTempVar  char   Variable name
!
        integer(i4b)::              iErr
        integer(HID_T)::            iFile, iSpaceID, iSetID
        integer(HSIZE_T) ::                        iDims(3)
        character(len=1024)::                acTemp, acTempFile, acTempVar; 
! *****************************************************************************
!
!   I/O
!
!   log file entries and export to HDF5 output file
!
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   PrintToLog
!   h5open_f
!   h5fcreate_f
!   h5screate_simple_f
!   h5dcreate_f
!   h5dwrite_f
!   h5sclose_f
!   h5dclose_f
!   h5fclose_f
!   h5close_f
!
! =============================================================================

        write (acTemp, '("Exporting ", A, " to file ", A)') trim(pszVarName), trim(pszFilename)//trim(int2str(iProcID))//".h5"
        call PrintToLog(acTemp, 1); 
        ! Allocate the buffer and fill it
        iDims = (/size(variable, 1), size(variable, 2), size(variable, 3)/); 
        ! Open the file and initialize the necessary data
        acTempFile = trim(pszFileName)//trim(int2str(iProcID))//".h5"; 
        acTempVar = trim(pszVarName); 
        call h5open_f(iErr); 
        CALL h5fcreate_f(trim(acTempFile), H5F_ACC_TRUNC_F, iFile, iErr)

        ! Export the slice
        !         Export the data
        call h5screate_simple_f(3, iDims, iSpaceID, iErr)
        call h5dcreate_f(iFile, trim(acTempVar), H5T_NATIVE_DOUBLE, iSpaceID, iSetID, iErr); 
        call h5dwrite_f(iSetID, H5T_NATIVE_DOUBLE, variable, iDims, iErr)
        call h5sclose_f(iSpaceID, iErr)

        ! Deinitialize
        call h5dclose_f(iSetID, iErr)
        call h5fclose_f(iFile, iErr); 
        call h5close_f(iErr)

    END SUBROUTINE Export3DArray

    SUBROUTINE Export3DComplexArray(pszFileName, pszVarName, variable, iProcID)

! =============================================================================
!
!   Programmer: Koos Huijssen / Jasper de Koning
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
!   The subroutine Export3DComplexArray exports a 3D complex array to a HDF5
!   output file.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   pszFileName  i   char   Filename of the export file
!   pszVarName   i   char   Name of the variable in the export file
!   variable     i   dpc    Variable to be stored
!   iProcID      i   i8b    ID of the current processor
!
        character(len=*), intent(in)::                pszFileName, pszVarName
        complex(dpc), intent(in) ::                                variable(:, :, :)
        integer(i8b), intent(in) ::                        iProcID

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   iErr       i4b    Error number
!   iFile     HID_T   HDF5 File ID
!   iSpaceID  HID_T   HDF5 Space ID
!   iSetID    HID_T   HDF5 Set ID
!   iDims    HSIZE_T  Length of the dimension(s) of the array
!   acTemp     char   Temporary char array for output log messages
!   acTempFile char   Filename including directory path
!   acTempVar  char   Variable name
!
        integer(i4b)::              iErr
        integer(HID_T)::            iFile, iSpaceID, iSetID
        integer(HSIZE_T) ::                        iDims(3)
        character(len=1024)::                acTemp, acTempFile, acTempVar; 
! *****************************************************************************
!
!   I/O
!
!   log file entries and export to HDF5 output file
!
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   PrintToLog
!   h5open_f
!   h5fcreate_f
!   h5screate_simple_f
!   h5dcreate_f
!   h5dwrite_f
!   h5sclose_f
!   h5dclose_f
!   h5fclose_f
!   h5close_f
!
! =============================================================================

        write (acTemp, '("Exporting ", A, " to file ", A)') trim(pszVarName), trim(pszFilename)//trim(int2str(iProcID))//".h5"
        call PrintToLog(acTemp, 1); 
        ! Allocate the buffer and fill it
        iDims = (/size(variable, 1), size(variable, 2), size(variable, 3)/); 
        ! Open the file and initialize the necessary data
        acTempFile = trim(pszFileName)//trim(int2str(iProcID))//".h5"; 
        acTempVar = trim(pszVarName); 
        call h5open_f(iErr); 
        CALL h5fcreate_f(trim(acTempFile), H5F_ACC_TRUNC_F, iFile, iErr)

        ! Export the slice
        !         Export the data
        call h5screate_simple_f(3, iDims, iSpaceID, iErr)
        call h5dcreate_f(iFile, trim(acTempVar), H5T_NATIVE_DOUBLE, iSpaceID, iSetID, iErr); 
        call h5dwrite_f(iSetID, H5T_NATIVE_DOUBLE, real(variable, dp), iDims, iErr)
        call h5dcreate_f(iFile, trim(acTempVar//"_cplx"), H5T_NATIVE_DOUBLE, iSpaceID, iSetID, iErr); 
        call h5dwrite_f(iSetID, H5T_NATIVE_DOUBLE, dimag(variable), iDims, iErr)
        call h5sclose_f(iSpaceID, iErr)

        ! Deinitialize
        call h5dclose_f(iSetID, iErr)
        call h5fclose_f(iFile, iErr); 
        call h5close_f(iErr)

    END SUBROUTINE Export3DComplexArray

    SUBROUTINE WriteAttributeDP(iFileID, iSetID, pszAttrName, iLen, arData)

! =============================================================================
!
!   Programmer: Koos Huijssen / Jasper de Koning
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
!   The subroutine WriteAttributeDP exports a dp 1D array as an attribute to
!   a HDF5 output file.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   iFileID      i   HID_T  HDF5 File ID
!   iSetID       i   HID_T  HDF5 Set ID
!   pszAttrName  i   char   Name of the attribute in the export file
!   iLen         i   i8b    Number of elements in the attribute variable
!   arData       i   dp     Attribute variable to be stored
!
        integer(HID_T), intent(in)::            iFileID, iSetID
        character(len=*), intent(in)::        pszAttrName
        integer(i8b), intent(in)::              iLen
        real(dp), intent(in)::                  arData(:)

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   iSpaceID  HID_T   HDF5 Space ID
!   iTypeID   HID_T   HDF5 Type ID
!   iAttrID    HID_T  HDF5 Attribute ID
!   iErr       i4b    Error number
!
        integer(HID_T)::                        iSpaceID, iTypeID, iAttrID
        integer(i4b)::                          iErr

! *****************************************************************************
!
!   I/O
!
!   export to HDF5 output file
!
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   h5tcopy_f
!   h5screate_simple_f
!   h5acreate_f
!   h5awrite_f
!   h5aclose_f
!   h5sclose_f
!   h5tclose_f
!
! =============================================================================

        call h5tcopy_f(H5T_NATIVE_DOUBLE, iTypeID, iErr); 
        call h5screate_simple_f(1, (/int(iLen, HSIZE_T)/), iSpaceID, iErr); 
        call h5acreate_f(iSetID, pszAttrName, iTypeID, iSpaceID, iAttrID, iErr)
        call h5awrite_f(iAttrID, iTypeID, arData, (/int(iLen, HSIZE_T)/), iErr); 
        call h5aclose_f(iAttrID, iErr)
        call h5sclose_f(iSpaceID, iErr); 
        call h5tclose_f(iTypeID, iErr); 
    END SUBROUTINE

    SUBROUTINE WriteAttributei8b(iFileID, iSetID, pszAttrName, iLen, aiData)

! =============================================================================
!
!   Programmer: Koos Huijssen / Jasper de Koning
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
!   The subroutine WriteAttributei8b exports an i8b 1D array as an attribute to
!   a HDF5 output file.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   iFileID      i   HID_T  HDF5 File ID
!   iSetID       i   HID_T  HDF5 Set ID
!   pszAttrName  i   char   Name of the attribute in the export file
!   iLen         i   i8b    Number of elements in the attribute variable
!   arData       i   i8b    Attribute variable to be stored
!
        integer(HID_T), intent(in)::            iFileID, iSetID
        character(len=*), intent(in)::        pszAttrName
        integer(i8b), intent(in)::              iLen
        integer(i8b), intent(in)::              aiData(:)

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   iSpaceID  HID_T   HDF5 Space ID
!   iTypeID   HID_T   HDF5 Type ID
!   iAttrID    HID_T  HDF5 Attribute ID
!   iErr       i4b    Error number
!
        integer(HID_T)::                        iSpaceID, iTypeID, iAttrID
        integer(i4b)::                          iErr

! *****************************************************************************
!
!   I/O
!
!   export to HDF5 output file
!
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   h5tcopy_f
!   h5screate_simple_f
!   h5acreate_f
!   h5awrite_f
!   h5aclose_f
!   h5sclose_f
!   h5tclose_f
!
! =============================================================================

        call h5tcopy_f(H5T_NATIVE_INTEGER, iTypeID, iErr); 
        call h5screate_simple_f(1, (/int(iLen, HSIZE_T)/), iSpaceID, iErr); 
        call h5acreate_f(iSetID, pszAttrName, iTypeID, iSpaceID, iAttrID, iErr)
        call h5awrite_f(iAttrID, iTypeID, int(aiData), (/int(iLen, HSIZE_T)/), iErr); 
        call h5aclose_f(iAttrID, iErr)
        call h5sclose_f(iSpaceID, iErr); 
        call h5tclose_f(iTypeID, iErr); 
    END SUBROUTINE

    SUBROUTINE Print_Error(error, err_unit_number, error_file)
!
!      Use parameters, Only: path
!     _________________________________________________________________________
        Implicit None
!
        Integer(4), Intent(In) :: err_unit_number
        Real(8), Intent(In) :: error
        Character(*), Intent(In) :: error_file

        Open (Unit=err_unit_number, File=trim(error_file), Status='old', &
       & Position='append')
        Write (err_unit_number, *) error
        Close (Unit=err_unit_number)

    END SUBROUTINE Print_Error

END MODULE ParnacOutput
