module ModLumping
    
use Mod_kind_param, only : wp
    
implicit none
!public module variables (many are defined in SETTINGS_2DLumping.txt input file):
character(len=5),public :: elcpname
integer,public :: cluster_num, nxintervals, nyintervals, YaxisChoice, vpMethod 
logical,public :: change_nd, massweightedKmeans, onlySurrogatesOutput
real(wp),public :: eladdMassConc, SetVolatThresholdVal

!--
!Explicit interfaces to subroutines in submodules SubModactcoeffratio and SubModLumpingSchemes:
interface
    module subroutine ActCoeffRatio_Volatility(nn, TKelvin, cpnameinp, filenameInp, folderpathout)
        integer,intent(in) :: nn                                    !nn is original number of neutral components
        real(wp),intent(in) :: TKelvin                              !the input temperature [K]
        character(len=*),dimension(nn),intent(in) :: cpnameinp      !list of assigned component names (here for neutral components only)
        character(len=*),intent(in) :: filenameInp                  !input file name, used to construct file names for input/output
        character(len=*),intent(in) :: folderpathout
    end subroutine ActCoeffRatio_Volatility

    module subroutine Lumping_schemes(nn, norg1, psatmethod, TK, psat, OtoC, HtoC, meanOS_C, &
                        & actcoeff_ratio, MolarMass, MassConc1, EVAP_paramA, EVAP_paramB, SMILES)
    integer,intent(in) :: norg1, nn  !nn is the number of neutral components
    character(len=*),intent(in) :: psatmethod
    real(wp),intent(in) :: TK
    real(wp),dimension(:),intent(in) :: psat, OtoC, HtoC, meanOS_C, actcoeff_ratio, &
                                        & MolarMass, MassConc1, EVAP_paramA, EVAP_paramB
    character(len=*),dimension(:),intent(in) :: SMILES
    end subroutine Lumping_schemes 
end interface
    
end module ModLumping