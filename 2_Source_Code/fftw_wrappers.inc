! include file for wrappers to use fftw libaries from erwin .dll files
! lines added to use .dll instead of .lib

!BEWARE: not all wrappers are here, there are also 4 wrapper definitions in ParnacConvInit.f90
! This is a workaround for a problem with a compiler error

!DEC$ ATTRIBUTES ALIAS: "dfftw_plan_dft_1d_"        :: dfftw_plan_dft_1d
!DEC$ ATTRIBUTES ALIAS: "dfftw_plan_dft_2d_"        :: dfftw_plan_dft_2d
!DEC$ ATTRIBUTES ALIAS: "dfftw_plan_dft_3d_"        :: dfftw_plan_dft_3d
!DEC$ ATTRIBUTES ALIAS: "dfftw_plan_dft_r2c_1d_"    :: dfftw_plan_dft_r2c_1d
!DEC$ ATTRIBUTES ALIAS: "dfftw_plan_dft_c2r_1d_"    :: dfftw_plan_dft_c2r_1d
!DEC$ ATTRIBUTES ALIAS: "dfftw_plan_dft_r2c_2d_"    :: dfftw_plan_dft_r2c_2d
!DEC$ ATTRIBUTES ALIAS: "dfftw_plan_dft_c2r_2d_"    :: dfftw_plan_dft_c2r_2d

!DEC$ ATTRIBUTES ALIAS: "dfftw_plan_many_dft_"      :: dfftw_plan_many_dft
!DEC$ ATTRIBUTES ALIAS: "dfftw_plan_many_dft_r2c_"  :: dfftw_plan_many_dft_r2c
!DEC$ ATTRIBUTES ALIAS: "dfftw_plan_many_dft_c2r_"  :: dfftw_plan_many_dft_c2r

!DEC$ ATTRIBUTES ALIAS: "dfftw_execute_"            :: dfftw_execute
!DEC$ ATTRIBUTES ALIAS: "dfftw_execute_dft_"        :: dfftw_execute_dft
!DEC$ ATTRIBUTES ALIAS: "dfftw_execute_dft_r2c_"    :: dfftw_execute_dft_r2c
!DEC$ ATTRIBUTES ALIAS: "dfftw_execute_dft_c2r_"    :: dfftw_execute_dft_c2r

!DEC$ ATTRIBUTES ALIAS: "dfftw_destroy_plan_"       :: dfftw_destroy_plan

!DEC$ ATTRIBUTES ALIAS: "dfftw_plan_guru_dft_"       :: dfftw_plan_guru_dft
!DEC$ ATTRIBUTES ALIAS: "dfftw_plan_guru_dft_c2r_"       :: dfftw_plan_guru_dft_c2r
!DEC$ ATTRIBUTES ALIAS: "dfftw_plan_guru_dft_r2c_"       :: dfftw_plan_guru_dft_r2c





