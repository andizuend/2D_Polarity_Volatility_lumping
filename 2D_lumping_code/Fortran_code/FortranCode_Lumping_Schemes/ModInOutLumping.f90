!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Module containing subroutines to read or write ASCII text files for input/output   *
!*   of component properties use with or produced via the product lumping schemes.      *
!*                                                                                      *
!*   :: Authors & Copyright ::                                                          *
!*   Ampritta Amaladhasan, Andi Zuend                                                   *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2018                                                            *
!*   -> latest changes: 2025-07-07                                                      *
!*                                                                                      *
!*   :: License ::                                                                      *
!*   This program is free software: you can redistribute it and/or modify it under the  *
!*   terms of the GNU General Public License as published by the Free Software          *
!*   Foundation, either version 3 of the License, or (at your option) any later         *
!*   version.                                                                           *
!*   The AIOMFAC model code is distributed in the hope that it will be useful, but      *
!*   WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or      *
!*   FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more      *
!*   details.                                                                           *
!*   You should have received a copy of the GNU General Public License along with this  *
!*   program. If not, see <http://www.gnu.org/licenses/>.                               *
!*                                                                                      *
!*   :: List of subroutines and functions contained in this module:                     *
!*   --------------------------------------------------------------                     *
!*   -  subroutine ReadCompVaporPressures                                               *
!*   -  subroutine ReadMolarComposition                                                 *
!*   -  subroutine OutputVaporPressureParam                                             *
!*   -  subroutine OutputComponentProperties                                            *
!*   -  subroutine OutputLumpedConcentrations                                           *
!*   -  subroutine OutputLumpedSMILES                                                   *
!*   -  subroutine Output_aw_levels                                                     *
!*   -  subroutine OutputConc_and_T                                                     *
!*   -  subroutine OutputGridlines                                                      *
!*   -  subroutine SMILES_to_AIOMFAC_inputFile                                          *
!*   -  subroutine Add_Electrolyte_Component                                            *
!*                                                                                      *
!****************************************************************************************
module ModInOutLumping

use Mod_kind_param, only : wp
use ModLumping, only : change_nd, onlySurrogatesOutput, YaxisChoice

implicit none
!public module variables:
character(len=900),public :: fpathin, fpathout
character(len=30),public :: chyaxis
character(len=4),public :: inpnum, outpnum
character(len=:),allocatable,public :: parentInpFile
!--
integer,public :: ndiout, linemax
!--
real(wp),parameter,public :: Pascal_per_atm = 101325.0E0_wp    !for vapour pressure unit conversion
real(wp),dimension(:),allocatable,public :: awlevels

contains  !module subroutines...

 
!==================================================================================================================================   
    subroutine ReadCompVaporPressures(nn, TKelvin, refcp1, vpMethod, EVAP_paramA, EVAP_paramB, EVAP_psat, SMILES)

    implicit none
    !interface variables:
    integer,intent(in) :: nn
    real(wp),intent(in) :: TKelvin
    integer,intent(in) :: refcp1, vpMethod
    real(wp),dimension(:),intent(out) :: EVAP_paramA, EVAP_paramB, EVAP_psat
    character(len=200),dimension(:),intent(out) :: SMILES
    !local variables:
    character(len=50) :: txtcheck
    character(len=20) :: dummy
    character(len=200) :: filename
    character(len=1000) :: fname
    integer :: i, istat, norg1, unitx, fileiostat, cpno1
    logical :: fileexists, filevalid
    real(wp),dimension(18) :: AllParams
    !...........................................
    
    !check if input vaporpressure file exists and read its content if true:
    filename = 'vaporpressure_'//inpnum//'.txt'
    fname = trim(fpathin)//trim(filename)
    inquire(FILE = fname, EXIST = fileexists) !inpfilesize is the file size in [bytes]
    !read file
    if (fileexists) then
        open(NEWUNIT = unitx, FILE = fname, IOSTAT=fileiostat, ACTION='read', STATUS='OLD')
        read(unitx,*) dummy, txtcheck           !read only first 2 words on this line for subsequent check
        if (txtcheck(1:6) == "SMILES") then     !indicates header of a correct input file
            filevalid = .true.
            !read input data from file:
            backspace unitx                     !jump back to beginning of record (to the beginning of the line)
            !read first line of headers
            read(unitx,*) dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy
        else !file not valid; it will be closed and deleted below
            filevalid = .false.
            write(*,*) ""
            write(*,*) "ERROR from ReadCompVaporPressures: Vapour pressure input file", filename," contains invalid data structure."
            read(*,*)
        endif
        !loop over smiles mixture components with log10 (pure component vapor pressures / atm) from Nanoolal method, EVAPORATION and the Parameters A, B from EVAPORATION model:
        if (filevalid) then
            do i = 1,nn !or until exit
                read(unitx,*,IOSTAT=istat) cpno1, SMILES(cpno1), AllParams(1:18)  !, EVAP_paramA(cpno1), EVAP_paramB(cpno1)
                !assign only the parameters for the selected pure-component vapor pressure calculation method
                EVAP_paramA(cpno1) = AllParams(vpMethod*2 - 1)
                EVAP_paramB(cpno1) = AllParams(vpMethod*2)
                if (istat /= 0) then    !end of file reached or read error
                    exit
                endif
            enddo
        endif !filevalid
        close(unitx)
    endif !fileexists 
    
    !add pure-component, liquid-state vapour pressure of water at given temperature:
    if (refcp1 == 1) then !water present
        norg1 = 2
        EVAP_psat(1) = exp(54.842763_wp - 6763.22_wp / TKelvin - 4.210_wp*log(TKelvin) + 0.000367_wp*TKelvin &
            & + tanh(0.0415_wp*(TKelvin - 218.8_wp))*(53.878_wp - 1331.22_wp / TKelvin &
            & - 9.44523_wp*log(TKelvin) + 0.014025_wp*TKelvin)) !psat of water in [Pa]
    else
        norg1 = 1
    endif    
    do i = norg1,nn
        EVAP_psat(i) = Pascal_per_atm * 10.0_wp**( EVAP_paramA(i) + (EVAP_paramB(i) / TKelvin**1.5_wp) ) ![Pa]; 10^(A+(B/T^1.5)) for vapor pressure calculation in atm converted to Pa                   
    enddo
    
    end subroutine ReadCompVaporPressures
!==================================================================================================================================


!!==================================================================================================================================
!    subroutine ReadCompVaporPressures_OLD(nn, TKelvin, refcp1, log10atmNanoolal, EVAP_paramA, EVAP_paramB, EVAP_psat, SMILES)
!
!    implicit none
!    !interface variables:
!    integer,intent(in) :: nn
!    real(wp),intent(in) :: TKelvin
!    integer,intent(in) :: refcp1
!    real(wp),dimension(:),intent(out) :: log10atmNanoolal, EVAP_paramA, EVAP_paramB, EVAP_psat
!    character(len=200),dimension(:),intent(out) :: SMILES
!    !local variables:
!    character(len=50) :: txtcheck
!    character(len=20) :: dummy
!    character(len=200) :: filename
!    character(len=1000) :: fname
!    integer :: i, istat, norg1, unitx, fileiostat, cpno1
!    real(wp) :: logpEVAP
!    logical :: fileexists, filevalid
!    !...........................................
!    
!    !check if input vaporpressure file exists and read its content if true:
!    filename = 'vaporpressure_'//inpnum//'.txt'
!    fname = trim(fpathin)//trim(filename)
!    inquire(FILE = fname, EXIST = fileexists) !inpfilesize is the file size in [bytes]
!    !read file
!    if (fileexists) then
!        open (NEWUNIT = unitx, FILE = fname, IOSTAT=fileiostat, ACTION='read', STATUS='OLD')
!        read(unitx,*) dummy, txtcheck           !read only first 2 words on this line for subsequent check
!        if (txtcheck(1:6) == "SMILES") then     !indicates header of a correct input file
!            filevalid = .true.
!            !read input data from file:
!            backspace unitx !jump back to beginning of record (to the beginning of the line)
!            read(unitx,*) dummy, dummy, dummy, dummy, dummy !read first line of headers
!        else !file not valid. It will be closed and deleted below
!            filevalid = .false.
!            write(*,*) ""
!            write(*,*) "ERROR from ReadCompVaporPressures: Vapour pressure input file", filename," contains invalid data structure."
!            read(*,*)
!        endif
!        !loop over smiles mixture components with log10 (pure component vapor pressures / atm) from Nanoolal method, EVAPORATION and the Parameters A, B from EVAPORATION model:
!        if (filevalid) then
!            do i = 1,nn !or until exit
!                read(unitx, *, IOSTAT=istat) cpno1, SMILES(cpno1), log10atmNanoolal(cpno1), logpEVAP, EVAP_paramA(cpno1), EVAP_paramB(cpno1)
!                if (istat /= 0) then !end of file reached or read error
!                    exit
!                endif
!            enddo
!        endif !filevalid
!        close(unitx)
!    endif !fileexists 
!    
!    !add pure-component, liquid-state vapour pressure of water at given temperature:
!    if (refcp1 == 1) then !water present
!        norg1 = 2
!        EVAP_psat(1) = exp(54.842763_wp - 6763.22_wp/TKelvin - 4.210_wp*log(TKelvin) + 0.000367_wp*TKelvin &
!            & + tanh(0.0415_wp*(TKelvin - 218.8_wp))*(53.878_wp - 1331.22_wp / TKelvin &
!            & - 9.44523_wp*log(TKelvin) + 0.014025_wp*TKelvin)) !psat of water in [Pa]
!        !add water psat also to Nanoolal psat on logscale:
!        log10atmNanoolal(1) = log10(EVAP_psat(1) / Pascal_per_atm)
!    else
!        norg1 = 1
!    endif    
!    do i = norg1,nn
!        EVAP_psat(i) = Pascal_per_atm * 10.0E0_wp**( EVAP_paramA(i) +(EVAP_paramB(i)/TKelvin**1.5_wp) ) ![Pa]; 10^(A+(B/T^1.5)) for vapor pressure calculation in atm converted to Pa
!    enddo
!    
!    end subroutine ReadCompVaporPressures_OLD
!!==================================================================================================================================

    
!==================================================================================================================================
    subroutine ReadMolarComposition(nn, ntmoles)

    implicit none
    !interface variables:
    integer,intent(in):: nn
    real(wp),dimension(:),intent(out) :: ntmoles 
    !local variables:
    character(len=200) :: filename
    character(len=1000) :: fname
    character(len=20) :: dummy, concunit
    integer :: istat, ncp, nn2, unitx
    real(wp),parameter :: NAvo = 6.02214076E23_wp       ![molec/mol] Avogadro's constant
    logical :: fileexists
    !.........................
    
    !Check for external total number of moles data file for components
    write(*,'(A,/)') "read molar amounts from external composition file"
    filename = 'input_concentrations_'//inpnum//'.txt'
    !-----
    fname = trim(fpathin)//trim(filename)
    inquire(FILE = fname, EXIST = fileexists)
    !!attempt to read a supposedly valid file
    if (fileexists) then
        nn2 = nn
        open (NEWUNIT = unitx, FILE = fname, ACTION='read', STATUS='OLD')
        read(unitx,*) dummy             !read just first word of first line, then go to next line
        read(unitx,*) dummy             !read second line
        read(unitx,*) concunit 
        read(unitx,*)                   !read empty line 
        do ncp = 1,nn
            read(unitx,*,IOSTAT=istat) ntmoles(ncp)  !save data read
            if (istat /= 0) then 
                nn2 = ncp -1
                exit                    !end of file reached or reading error.
            endif
        enddo
        close(unitx)
    else
        write(*,*) ""
        write(*,*) "ERROR in ReadMolarComposition: total moles for components input file not found; file path: ", fname
        read(*,*)
    endif !fileexists
    !check whether the last component's concentration is missing (presumably because water was not added to the concentration file):
    if (nn2 < nn) then
        if (nn - nn2 == 1) then                         !assume water concentration was missing in input file, so shift index by one and add a dummy amount for water:
            block 
                real(wp),dimension(nn) :: ntmoles_copy
                ntmoles_copy(1:nn2) = ntmoles(1:nn2)
                ntmoles(1) = 5.0E-01_wp                 !dummy amount for component 1 = water;
                ntmoles(2:nn) = ntmoles_copy(1:nn2)
            end block
        else
            write(*,*) ""
            write(*,*) "ERROR in ReadMolarComposition: concentration input rows from file differ in number of components from expected number by: ", nn - nn2
            read(*,*)
        endif
    endif
    !check whether unit conversion is necessary:
    if (concunit(1:6) == '[molec') then
        ntmoles = (ntmoles / (NAvo*1.0E-6_wp))            ![molec/cm^3] to [mol/m^3] conversion
    endif
    
    end subroutine ReadMolarComposition
!==================================================================================================================================

    
!==================================================================================================================================
    subroutine OutputVaporPressureParam(nn, norg1, filename, Lmethod, EVAP_paramA, EVAP_paramB, SMILES, MassConc)
    
    implicit none
    !interface variables:
    integer,intent(in):: nn, norg1
    character(len=*),intent(in) :: filename, Lmethod
    real(wp),dimension(:),intent(in) :: EVAP_paramA, EVAP_paramB
    character(len=200),dimension(:),intent(in) :: SMILES
    real(wp),dimension(:),intent(in) :: MassConc
    !local variables:
    character(len=1) :: kn
    character(len=4) :: Iformat
    character(len=60) :: tformat
    character(len=300) :: title
    character(len=1000) :: fname
    integer :: i, k, unitx
    logical :: lumpedsyst
    !.........................
    
    if (Lmethod(1:4) == 'Full') then
        lumpedsyst = .false.
    else
        lumpedsyst = .true.
    endif
    k = max(2, ceiling(log10(real(nn))) )       !determined size of integer format specifier
    write(kn,'(I0)') k
    Iformat = trim('I'//kn//'.'//kn)
    
    fname = trim(fpathout)//trim(filename)
    open (NEWUNIT = unitx, FILE = fname, STATUS = "UNKNOWN")
    title = "Lumping method ** "//trim(Lmethod)//" ** | Pure organic component vapor pressure, param. A and B computed using EVAPORATION method | YaxisChoice = "//trim(chyaxis)
    write(tformat,'(I0)') len_trim(title)
    tformat = '(A'//trim(tformat)//')'
    write(unitx,'(A)') "==================================================================================================================="
    write(unitx,tformat) adjustl(title)
    write(unitx,'(A)') "==================================================================================================================="
    write(unitx,*) ""
    k = count(MassConc(norg1:nn) > 0.0_wp)
    write(unitx,'(A,I0,A,I0)') "nneutral (full system) : ", nn, "; # of organic surrogate components: ", k
    write(unitx,*) ""
    title = "  EVAPParam_A,      EVAPParam_B,       SMILES"
    write(unitx,'(A)') title
    do i = norg1,nn
        if (onlySurrogatesOutput .and. lumpedsyst) then !filter out the surrogate components only for output
            if (MassConc(i) > 0.0_wp) then
                write(unitx,'(2(ES13.6,",",4X), A,2X)') EVAP_paramA(i), EVAP_paramB(i), SMILES(i)
            endif
        else
            write(unitx,'(2(ES13.6,",",4X), A,2X)') EVAP_paramA(i), EVAP_paramB(i), SMILES(i)
        endif
    enddo !i
    close(unitx)

    end subroutine OutputVaporPressureParam
!==================================================================================================================================

    
!==================================================================================================================================
    subroutine OutputComponentProperties(nn, norg1, TKelvin, filename, Lmethod, psatmethod, resol, psat, OtoC, HtoC, &
                                       & meanOS_C, actcoeff_ratio, MolarMass, MassConc)
    
    implicit none
    !interface variables:
    integer,intent(in) :: nn, norg1
    real(wp),intent(in) :: TKelvin 
    character(len=*),intent(in) :: filename, Lmethod, psatmethod, resol
    real(wp),dimension(:),intent(in) :: psat, OtoC, HtoC, meanOS_C, actcoeff_ratio, MolarMass, MassConc
    !local variables:
    real(wp),parameter :: Rgas = 8.314462618_wp
    character(len=1) :: kn
    character(len=4) :: Iformat
    character(len=60) :: tformat
    character(len=300) :: title
    character(len=1000) :: fname
    integer :: i, k, unitx
    logical :: lumpedsyst
    !......................................
    
    if (Lmethod(1:4) == 'Full') then
        lumpedsyst = .false.
    else
        lumpedsyst = .true.
    endif
    k = max(2, ceiling(log10(real(nn))) )       !determined size of integer format specifier
    write(kn,'(I0)') k
    Iformat = trim('I'//kn//'.'//kn)

    fname = trim(fpathout)//trim(filename)
    open (NEWUNIT = unitx, FILE = fname, STATUS = "UNKNOWN")
    title = "Lumping method ** "//trim(Lmethod)//" ** | Pure organic component properties | psat calculated by "//trim(psatmethod)//" method"
    write(unitx,'(A)') "==================================================================================================================="
    write(unitx,'(A)') adjustl(title)
    title = "YaxisChoice = "//trim(chyaxis)
    write(unitx,'(A)') adjustl(title)
    title = "Grid resolution = "//trim(resol)
    write(unitx,'(A)') adjustl(title)
    write(unitx,'(A)') "==================================================================================================================="
    write(unitx,*) ""
    
    k = count(MassConc(norg1:nn) > 0.0_wp)
    write(unitx,'(A,I0,A,I0)') "nneutral (full system) : ", nn, "; # of organic surrogate components: ", k
    write(unitx,*) ""
    title = "Cpno.,   OtoC_ratio,    HtoC_ratio,      meanOS_C,    Pure_Comp_VP_[Pa], actcoeff_ratio, Molar_mass_[kg/mol], C0_[ug/m^3]"
    write(unitx,'(A)') adjustl(title)
    tformat = '('//Iformat//',",",A,*(ES13.6,:,",",A))'
    do i = norg1,nn
        if (onlySurrogatesOutput .and. lumpedsyst) then !filter out the surrogate components only for output
            if (MassConc(i) > 0.0_wp) then
                write(unitx, tformat) i, char(9), OtoC(i), char(9), HtoC(i), char(9), meanOS_C(i), char(9), psat(i), char(9), &
                    & actcoeff_ratio(i), char(9), MolarMass(i), char(9), psat(i)*MolarMass(i)/(Rgas*TKelvin)*1.0E9_wp
            endif
        else
            write(unitx, tformat) i, char(9), OtoC(i), char(9), HtoC(i), char(9), meanOS_C(i), char(9), psat(i), char(9), &
                & actcoeff_ratio(i), char(9), MolarMass(i), char(9), psat(i)*MolarMass(i)/(Rgas*TKelvin)*1.0E9_wp
        endif
    enddo
    close(unitx)
    
    end subroutine OutputComponentProperties
!==================================================================================================================================


!==================================================================================================================================
    subroutine OutputLumpedConcentrations(nn, norg1, filename, Lmethod, psatmethod, resol, MassConcOrig, &
                                        & MassConcLumped, MolarMass)
    
    implicit none
    !interface variables:
    integer,intent(in) :: nn, norg1
    character(len=*),intent(in) :: filename, Lmethod, psatmethod, resol
    real(wp),dimension(:),intent(in) :: MassConcOrig, MassConcLumped, MolarMass
    !local variables:
    character(len=1) :: kn
    character(len=4) :: Iformat
    character(len=100) :: tformat
    character(len=1000) :: fname
    integer :: i, k, unitx
    logical :: lumpedsyst
    !......................................
    
    if (Lmethod(1:4) == 'Full') then
        lumpedsyst = .false.
    else
        lumpedsyst = .true.
    endif
    k = max(2, ceiling(log10(real(nn))) )       !determined size of integer format specifier
    write(kn,'(I0)') k
    Iformat = trim('I'//kn//'.'//kn)

    fname = trim(fpathout)//trim(filename)
    open (NEWUNIT = unitx, FILE = fname, STATUS = "UNKNOWN")
    write(unitx,'(A)') "==================================================================================================================="
    write(unitx,'(A)') "Lumping method ** "//trim(Lmethod)//" ** | Organic surrogate components total (gas + condensed) composition | psat calculated by "//trim(psatmethod)//" method"
    write(unitx,'(A)') 'YaxisChoice = '//trim(chyaxis)
    write(unitx,'(A)') "Grid resolution = "//trim(resol)
    write(unitx,'(A)') "==================================================================================================================="
    write(unitx,*) ""
    k = count(MassConcLumped(norg1:nn) > 0.0_wp)
    write(unitx,'(A,I0.4,A,I0.4)') "nneutral (full system) : ", nn, "; # of organic surrogate components: ", k
    write(unitx,*) ""
    write(unitx,'(A)') "Comp_ID, TotMassConc_lumped_[ug/m^3], ntot_lumped_[mol/m^3], TotMassConc_original_[ug/m^3], ntot_original_[mol/m^3], Molar_mass_[kg/mol]"
    tformat = '('//Iformat//',",", *(7X,ES13.6,:,","))'
    do i = norg1,nn
        if (onlySurrogatesOutput .and. lumpedsyst) then     !filter out the surrogate components only for output
            if (MassConcLumped(i) > 0.0_wp) then
                write(unitx, tformat) i, MassConcLumped(i), MassConcLumped(i)/(MolarMass(i)*1.0E9_wp), MassConcOrig(i), MassConcOrig(i)/(MolarMass(i)*1.0E9_wp), MolarMass(i)
            endif
        else
            write(unitx, tformat) i, MassConcLumped(i), MassConcLumped(i)/(MolarMass(i)*1.0E9_wp), MassConcOrig(i), MassConcOrig(i)/(MolarMass(i)*1.0E9_wp), MolarMass(i)
        endif
    enddo
    close(unitx)
        
    end subroutine OutputLumpedConcentrations
!==================================================================================================================================

    
!==================================================================================================================================
    !write files serving as input for SMILES processing to create AIOMFAC input files for lumped surrogate systems
    subroutine OutputLumpedSMILES(nn, norg1, filename, MassConcLumped, SMILES)
    
    implicit none
    !interface variables:
    integer,intent(in) :: nn, norg1
    character(len=*),intent(in) :: filename
    real(wp),dimension(:),intent(in) :: MassConcLumped
    character(len=200),dimension(:),intent(in) :: SMILES
    !local variables:
    character(len=1000) :: fname
    integer :: i, unitx
    !......................................

    fname = trim(fpathout)//trim(filename)
    open (NEWUNIT = unitx, FILE = fname, STATUS = "UNKNOWN")
    do i = norg1,nn
        if (MassConcLumped(i) > 0.0_wp) then
            write(unitx, '(A200)') adjustl(SMILES(i))
        endif
    enddo
    close(unitx)
    
    end subroutine OutputLumpedSMILES
!==================================================================================================================================
    

!==================================================================================================================================
    !write files serving as input for SMILES processing to create AIOMFAC input files for lumped surrogate systems
    subroutine Output_aw_levels(filename, awlevels)
    
    implicit none
    !interface variables:
    character(len=*),intent(in) :: filename
    real(wp),dimension(:),intent(in) :: awlevels
    !local variables:
    character(len=1000) :: fname
    integer :: k, kmax, unitx
    !......................................
    
    fname = trim(fpathout)//trim(filename)
    open (NEWUNIT = unitx, FILE = fname, STATUS = "UNKNOWN")
    write(unitx, '(A)') "point,   aw_level (bulk RH)"
    kmax = size(awlevels)
    do k = 1,kmax
        write(unitx, '(I2.2,",",3X,F8.5)') k, awlevels(k)
    enddo
    wait(unitx)
    close(unitx)
    
    end subroutine Output_aw_levels
!==================================================================================================================================
    
    
!==================================================================================================================================
    !write files serving as input for AIOMFAC gas-particle partitioning calculations
    subroutine OutputConc_and_T(nn, norg1, linemax, filename, Lmethod, TK, MassConc, MolarMass, elcpname, eladdMassConc)
    
    implicit none
    !interface variables:
    integer,intent(in) :: nn, norg1, linemax
    character(len=*),intent(in) :: filename, Lmethod 
    real(wp),intent(in) :: TK
    real(wp),dimension(:),intent(in) :: MassConc, MolarMass
    character(len=*),intent(in) :: elcpname
    real(wp),intent(in) :: eladdMassConc
    !local variables:
    character(len=10) :: cpnstring
    character(len=60) :: tformat
    character(len=300) :: textline
    character(len=1000) :: fname
    integer :: i, m, line, nlast, unitx
    real(wp) :: elMM
    real(wp),dimension(size(MassConc)+1) :: ntlist   !total molar concentration [mol/m^3] of air, for output;
    logical :: lumpedsyst
    !......................................

    if (Lmethod(1:4) == 'Full') then
        lumpedsyst = .false.
    else
        lumpedsyst = .true.
    endif
    
    fname = trim(fpathout)//trim(filename)
    open (NEWUNIT = unitx, FILE = fname, STATUS = "UNKNOWN")
    write(unitx, '(A)') "mixture composition and temperature:                                                                "
    write(unitx, '(A)') "mass fraction? 0  !here composition is given as nt (total number of moles in system per m^3)        "
    write(unitx, '(A)') "mole fraction? 2  !set to 2 to indicate moles per m^3 (nt) as input                                 "
    write(unitx, '(A)') "----"
    if (onlySurrogatesOutput .and. lumpedsyst) then
        m = norg1-1
        do i = norg1,nn
            if (MassConc(i) > 0.0_wp) then
                m = m+1
                ntlist(m) = MassConc(i) / (MolarMass(i)*1.0E9_wp)    ![mol/m^3]  only save the surrogates
            endif
        enddo
        nlast = m
    else
        nlast = nn
        ntlist(1:nn) = MassConc(1:nn) / (MolarMass(1:nn)*1.0E9_wp)   ![mol/m^3]  save all
    endif
    !add potential electrolyte component amount:
    if (trim(elcpname) /= 'none') then
        nlast = nlast+1
        !select molar mass of added electrolyte component:
        select case(trim(elcpname))
        case('AS')
            elMM = 0.132134_wp    !molar mass in [kg/mol]
        case('abs')
            elMM = 0.115103_wp   
        case default
            elMM = 0.0_wp        !this will cause division by zero below and alert that an uncovered option was chosen;  
        end select
        ntlist(nlast) = eladdMassConc / (elMM*1.0E9_wp)     ![mol/m^3]
    endif
    
    write(cpnstring,'(I0)') nlast
    tformat = '(I3.3,",",2X,F6.2,",",2X, *(ES13.6,:,",",1X))'
    textline = "point,  T_K,  cp02, cp03, cp04, ..., cp"//trim(cpnstring)
    write(unitx, '(A)') adjustl(textline)
    do line = 1,linemax !linemax here indicates 'linemax' water activity (RH) levels used in the gas-particle partitioning calculations with AIOMFAC
        !the output lines are repeated kmax times.
        write(unitx, tformat) line, TK, (ntlist(i), i = norg1,nlast)
    enddo        
    wait(unitx)
    close(unitx)
    
    end subroutine OutputConc_and_T
!==================================================================================================================================
    
    
!==================================================================================================================================
    !write files serving as input for AIOMFAC gas-particle partitioning calculations
    subroutine OutputGridlines(xgridl, ygridl, filename)
    
    implicit none
    !interface variables:
    real(wp),dimension(:),intent(in) :: xgridl, ygridl
    character(len=*),intent(in) :: filename
    !local variables:
    integer :: nxi, nyi, un
    character(len=30) :: tformat
    !.....................................  
    
    nxi = size(xgridl)
    nyi = size(ygridl)
    open (NEWUNIT = un, FILE = trim(fpathout)//filename, STATUS='UNKNOWN')
    write(un,'(A)') "Grid coordinates for the axis system used. Values may be log10 depending on choice of axis system."
    write(un,*) ""
    write(un,'(A)') 'YaxisChoice = '//trim(chyaxis)
    write(un,*) ""
    write(un,'(A,I3.3)') "#-x-Gridlines: ", nxi
    write(un,'(A,I3.3)') "#-y-Gridlines: ", nyi
    write(un,*) ""
    tformat = '(*(2X,ES13.6,:,","))'
    write(un, '(A)') "x-grid line coordinates:"
    write(un, tformat) xgridl(1:nxi)  !the log10 values here..., needs to be considerd in custom plotting
    tformat = '(*(2X,ES13.6,:,","))'
    write(un, '(A)') "y-grid line coordinates:"
    write(un, tformat) ygridl(1:nyi)
    close(un)
    
    end subroutine OutputGridlines
!==================================================================================================================================
    
!==================================================================================================================================
!write files serving as input for AIOMFAC gas-particle partitioning calculations
    subroutine OutputKmeansCluster(nn, norg1, dim_num, cluster_num, cluster, cluster_center, cluster_population, cluster_surrogate, filename)
    
    implicit none
    !interface variables:
    integer, intent(in) :: nn, norg1, dim_num, cluster_num
    integer,dimension(:),allocatable,intent(in) :: cluster                               !dimension is (norg1:nn)
    real(wp),dimension(:,:),allocatable,intent(in)  :: cluster_center                        !dimensions are (1:dim_num, 1:cluster_num)
    integer,dimension(:),allocatable,intent(in) :: cluster_population, cluster_surrogate !dimension is (1:cluster_num)
    character(len=*),intent(in) :: filename
    !local variables:
    integer :: j, un
    !.....................................  
    
    open (NEWUNIT = un, FILE = trim(fpathout)//filename, STATUS='UNKNOWN')
    write(un,'(A)') 'Selected properties of the different identified k-means clusters.'
    write(un,*) ''
    write(un,'(A)') 'YaxisChoice = '//trim(chyaxis)
    write(un,'(A,1X,I0)') 'Number of clusters = ', cluster_num
    write(un,*) ''
    write(un,'(A)') 'legend:  cc-x, cc-y = cluster center x and y coordinates (log-scale where applicable); icsc = index of cluster surrogate component; cpop = cluster population number;'
    write(un,*) ''
    write(un,'(A)') 'No.        cc-x           cc-y       icsc    cpop'
    write(un,'(A)') '-------------------------------------------------'
    do j = 1,cluster_num
        write(un,'(I0.3,2X,2(ES13.5,2X),2X,I0.4,4X,I0.4)') j, cluster_center(1:dim_num,j), cluster_surrogate(j), cluster_population(j) 
    enddo
    write(un,*) ''
    write(un,*) ''
    write(un,'(A)') 'List of components and to which cluster they belong:'
    write(un,'(A)') 'Comp_no.,  cluster_no'
    write(un,'(A)') '---------------------'
    do j = norg1,nn
        write(un,'(I0.5,",",5X,I0.5)') j, cluster(j)        
    enddo
    
    close(un)
    
    end subroutine OutputKmeansCluster
!==================================================================================================================================


!==================================================================================================================================
    !Generate an AIOMFAC-web style input file using only the SMILES of the determined surrogate components.
    !This subroutine will do this via access to our Python SMILES conversion code (from the folder SMILES_to_AIOMFAC_inp).
    !For this subroutine to work properly, Python needs to be installed as well as other packages used in 
    !SMILES_to_AIOMFAC_inp, such as the epam Indigo API.
    subroutine SMILES_to_AIOMFAC_inputFile(Smilesfname, outfilename)
    
    use ModStringFunctions, only : Replace_Text
    use ModOScommands, only : isWindowsPlatform, copy_file
    
    implicit none
    !interface variables:
    character(len=*),intent(in) :: Smilesfname, outfilename
    !local variables:
    character(len=200) :: inpPath, pathTools
    character(len=300) :: command, command2, f1pathName, f2pathName
    integer :: Estat, Cstat
    logical :: isWindowsOS
    !.................
    
    !define commands to call Python script for SMILES--SMARTS--AIOMFAC subgroup conversion and input file writing:
    inpPath = '../Output_lumping/'                  !relative path for accessing smiles_XXXX.txt file.
    pathTools = '../../SMILES_to_AIOMFAC_inp/'      !relative path from 2D lumping framework project directory to the SMILES_to_AIOMFAC_inp tool
    
    !(1) For ease of use of local folder structure, the SMILES file for conversion is first copied to the InputFiles 
    !folder of "SMILES_to_AIOMFAC_inp";
    f1pathName = trim(inpPath)//Smilesfname
    f2pathName = trim(pathTools)//'InputFiles/'//Smilesfname
    call copy_file(f1pathName, f2pathName)
    
    !(2) run the Python script for SMILES conversion to input file:    
    isWindowsOS = isWindowsPlatform() 
    command = 'python '//trim(pathTools)//'SMILES_to_AIOMFAC_input.py '//trim(pathTools)//'InputFiles/'//trim(Smilesfname)
    if (isWindowsOS) then
        command2 = Replace_Text(command, "/", "\") !replace forward by backslashes for Windows commands
    else !Linux, MacOS, ...
        command2 = command
    endif
    call execute_command_line(trim(command2), EXITSTAT=Estat, CMDSTAT=Cstat)
    
    !(3) copy the produced output file ('input_XXXX.txt") back to the 2D-Lumping framework output folder:
    write(*,*) '... copying the AIOMFAC input file to the local "Output_lumping" directory;'
    f1pathName = trim(pathTools)//'OutputFiles/'//trim(outfilename)
    f2pathName = trim(inpPath)//trim(outfilename)
    call copy_file(f1pathName, f2pathName)
    
    end subroutine SMILES_to_AIOMFAC_inputFile
!==================================================================================================================================

    
!==================================================================================================================================
    !Modify the generated AIOMFAC-web style input file (see SMILES_to_AIOMFAC_inputFile) to include an electrolyte component
    !after the last organic (neutral) component in the input file.
    subroutine Add_Electrolyte_Component(outfilename, eltype, elMassC)
    
    use ModStringFunctions, only : Replace_Text
    use ModOScommands, only : isWindowsPlatform, copy_file
    
    implicit none
    !interface variables:
    character(len=*),intent(in) :: outfilename, eltype 
    real(wp) :: elMassC
    !local variables:
    character(len=200) :: inpPath
    character(len=300) :: f1pathName
    character(len=9)   :: cpn
    character(len=25)  :: txtcheck, elname
    character(len=:),dimension(:),allocatable :: datarows
    integer :: i, ilen, iplusrow, ilast, istat, maxlen, nncp, trimdatlen, un
    logical :: save_rows
    !............................
    
    if (trim(eltype) == 'none') then
        return      !leave subroutine since no electrolyte was selected to be added    
    endif
    
    inpPath = '../Output_lumping/'                      !relative path for accessing an input_XXXX.txt file.
    f1pathName = trim(inpPath)//trim(outfilename)
    
    !(1) open file, read row data until finding the ++++ indicator of having reached past the last componen,
    !    then save the data in the last rows containing composition info for appending later;
    open (NEWUNIT=un, FILE=f1pathName, IOSTAT=istat, ACTION='readwrite', STATUS='OLD')
    save_rows = .false.
    ilen = 300
    maxlen = ilen
    do while (maxlen > ilen-12)                         !this outer loop will repeat the below procedure in cases where the data rows near the end are very long;
        ilen = maxlen + maxlen
        if (save_rows) then
            deallocate(datarows)
            rewind(un)
            save_rows = .false.
        endif
        i = 0
        do                                              !loop over rows in input file
            i = i +1
            if (save_rows) then                         !record the data in the rows as part of a character array
                read(un,'(A)',IOSTAT=istat) datarows(i) !read whole row
            else
                read(un,*,IOSTAT=istat) txtcheck        !read only first word/argument for subsequent check
            endif
            if (save_rows) then
                trimdatlen = len_trim(datarows(i))      !evaluate the character data row length;
                maxlen = max(trimdatlen, maxlen)
                if (datarows(i)(1:4) == '====') then    !'====' indicates the last row in the input file
                    ilast = i
                    exit
                endif
            else if (txtcheck(1:4) == '++++') then      !"++++" indicates end of components definition part
                iplusrow = i
                allocate( character(len=ilen) :: datarows(iplusrow+1:iplusrow+8) )
                datarows(:) = ''
                txtcheck = ''
                save_rows = .true.
            else if (istat /= 0) then                   
                write(*,'(A)') "WARNING from Add_Electrolyte_Component: istat /= 0 was detected while reading input file!"
                read(*,*)
                exit
            endif
        enddo !i
    enddo !maxlen
    
    !(2) backspace to file row position of last component and add the additional electrolyte component:
    do i = 1,(ilast+1 -iplusrow)
        backspace(un)
    enddo
    !determine last neutral component number:
    i = index(datarows(iplusrow+5), 'cp', BACK = .TRUE.)
    read(datarows(iplusrow+5)(i+2:i+12), *) nncp
    write(un,'(A,A,I0.2)')  'component no.:', achar(9), nncp+1    !achar(9) is a tab character for tab-delimited text file;
    select case(trim(eltype))
    case('AS')
        elname = "'"//"(NH4)2SO4"//"'"
        write(un,'(A,A,A)')     'component name:', achar(9), trim(elname)
        write(un,'(A,A,I0.3,A,A,I0.2)') 'subgroup no., qty:', achar(9), 204,',', achar(9), 02
        write(un,'(A,A,I0.3,A,A,I0.2)') 'subgroup no., qty:', achar(9), 261,',', achar(9), 01
    case('abs')
        elname = "'"//"NH4HSO4"//"'"
        write(un,'(A,A,A)')     'component name:', achar(9), trim(elname)
        write(un,'(A,A,I0.3,A,A,I0.2)') 'subgroup no., qty:', achar(9), 204,',', achar(9), 01
        write(un,'(A,A,I0.3,A,A,I0.2)') 'subgroup no., qty:', achar(9), 248,',', achar(9), 01
    case default         
        write(*,*) "WARNING from Add_Electrolyte_Component: no valid electrolyte for addition identified"
        read(*,*)
    end select
    write(un,'(A)')  '----'
    write(un,'(A)')  '++++'
    !modify the relevant datarows text where the new component info needs to be added too:
    write(cpn, '(I0)') nncp+1
    datarows(iplusrow+5) = trim(datarows(iplusrow+5))//', el_cp'//trim(cpn)
    if (elMassC <= 0.0_wp) then
        datarows(iplusrow+6) = trim(datarows(iplusrow+6))//', 0.0e+00'
    else
        datarows(iplusrow+6) = trim(datarows(iplusrow+6))//', 1.0e-12'
    endif
    do i = iplusrow+1,ilast
        write(un,'(A)') trim(datarows(i))
    enddo
    close(un)
    
    end subroutine Add_Electrolyte_Component
!==================================================================================================================================
    
    
end module ModInOutLumping