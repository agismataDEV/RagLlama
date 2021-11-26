MODULE ParnacTransformDef

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
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The module TransformDef contains the definition of a structure that
!   contains the plans of various FFT transforms. The structure forms a
!   part of the Grid structure, which forms a part of the Space structure.
!   These are defined in ParnacDataDef.
!
! *****************************************************************************
!
!   MODULES USED
!	
    USE Types
!
! *****************************************************************************
!
!   GLOBAL DECLARATIONS
!
!   Transforms    type   Definition of structure containing the plans for the 
!                        FFT transforms for a space. The plans are set in 
!                        GridInitTransforms.
!
! =============================================================================
!	
    IMPLICIT NONE

	type Transforms
		!Standard transforms, used in the convolution operator
		integer(i8b) ::			iPlanTransform1D, iPlanTransform1D_inv
		integer(i8b) ::			iPlanTransformT, iPlanTransformT_inv		
		integer(i8b) ::			iPlanTransformXY, iPlanTransformXY_inv
		integer(i8b) ::			iPlanTransformXYZ, iPlanTransformXYZ_inv

		!KH these transforms are sufficient for the convolution. However, for efficiency we also implement
		integer(i8b) ::			iPlanTransformXY_relevantZ, iPlanTransformZ 
		integer(i8b) ::			iPlanTransformZ_inv, iPlanTransformXY_relevantZ_inv 
		!KH since for the contrast source a large upper z part is zero (lower part forward XY transform, then full Z),
		!KH and for the field result, a large upper z part is unimportant (full inverse Z, then lower part inverse XY transform)
		!KH Sizes are determined by iDimZ and iBeamindexstartZ

		!Small T transforms, used in the nonlinear operator
		integer(i8b) ::			iPlanTransformT_sml, iPlanTransformT_sml_inv
		!Large T transforms, used in the nonlinear operator
		integer(i8b) ::			iPlanTransformT_lrg, iPlanTransformT_lrg_inv 		

		!Small XYZ transforms, used in the linear contrast operator
		integer(i8b) ::			iPlanTransformXYZ_sml, iPlanTransformXYZ_sml_inv 		

	end type Transforms
	
END MODULE ParnacTransformDef
