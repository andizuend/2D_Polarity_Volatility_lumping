submodule (ModLumping) SubModactcoeffratio

contains
    
!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Subroutine to compute component-specific activity coeff. ratios with respect to    * 
!*   separate, binary dissolution in two distinct solvents (water or hexanediol)        *
!*                                                                                      *
!*   :: Authors & Copyright ::                                                          *
!*   Dalrin Ampritta Amaladhasan, Andi Zuend, Dan Hassan-Barthaux                       *
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
!****************************************************************************************
module subroutine ActCoeffRatio_Volatility(nn, TKelvin, cpnameinp, filenameInp, folderpathout)

!used module procedures and variables:
use Mod_kind_param, only : wp
use ModSystemProp, only : cpname, compN, definemixtures, nneutral, Mmass, waterpresent, ITAB
use ModSubgroupProp, only : O2C_H2C_component, topsubno
use ModAIOMFACvar, only : wtf, actcoeff_n
use ModCalcActCoeff, only : AIOMFAC_calc
use ModInOutLumping, only : ReadCompVaporPressures, ReadMolarComposition, Add_Electrolyte_Component, &
    & OutputComponentProperties, OutputLumpedConcentrations, OutputLumpedSMILES, Output_aw_levels,   &
    & OutputConc_and_T, OutputVaporPressureParam, awlevels, chYaxis, linemax, fpathin, fpathout,     &
    & folderpathout, inpnum, parentInpFile, ndiout, outpnum
use ModOScommands, only : copy_file
use qsort_c_module, only : QsortC

implicit none
!-- interface variables:
integer,intent(in) :: nn                                    !nn is original number of neutral (non-electrolyte) components
real(wp),intent(in) :: TKelvin                              !the input temperature [K]
character(len=*),dimension(nn),intent(in) :: cpnameinp      !list of assigned component names (here for neutral components only)
character(len=*),intent(in) :: filenameInp                  !input file name, used to construct file names for input/output
character(len=*),intent(in) :: folderpathout
!-- local variables:
character(len=15) :: psatmethod, Lmethod, resol
character(len=200) :: filename, filein, fileout
character(len=200),dimension(:),allocatable :: SMILES
character(len=200),dimension(2) :: cpnameRefcp 
character(len=100),dimension(9) :: methodNames 
!--
integer :: i, k, ndi, ncp, addi, addn, nc, ncomp, refcp1, refcp2, norg1, nsubset, nstart, nend, nTemp
integer,dimension(:),allocatable :: compIDdat, compIDdatTemp
integer,dimension(:,:),allocatable :: cpsubgdat, cpsubgdatTemp
!--
real(wp),parameter :: Tref = 298.15_wp, ntiny = sqrt(tiny(1.0_wp))
real(wp) :: cpCarbon, cpHydrogen, cpOxygen, cpNitrogen, cpSulfur, wtftestcp
real(wp),dimension(:),allocatable :: actcoeff_ratio, actcoeffRef1, actcoeffRef2, EVAP_paramA, EVAP_paramB, &
    & EVAP_psat, MolarMass, ntmoles, OtoC, HtoC, NtoC, StoC, meanOS_C, psat, psatCopy, TotalMassConc
!--
logical :: duplicates
!........................................................

!######## Set model parameters for this run (other parameters read from SETTINGS input file in Main_IO_Lumping) #########

!set water activity (RH) levels
awlevels = [0.99_wp, 0.98_wp, 0.96_wp, 0.95_wp, 0.94_wp, 0.92_wp, 0.90_wp, 0.85_wp, 0.80_wp, 0.75_wp, &
          & 0.70_wp, 0.65_wp, 0.60_wp, 0.50_wp, 0.40_wp, 0.30_wp, 0.20_wp, 0.15_wp, 0.10_wp, 0.01_wp, 0.001_wp]
linemax = size(awlevels)        !number of water activity (RH) levels to be used in AIOMFAC partitioning calculations
methodNames(1:9) = ["VP_N_BP_N", "VP_N_BP_SB", "VP_N_BP_JR", "VP_MY_BP_N", &
                    & "VP_MY_BP_SB", "VP_MY_BP_JR", "EVAP", "EVAP2", "SIMPOL"]

!####################################################

!initialize:
nneutral = 0
fpathout = folderpathout
k = len_trim(filenameInp)
allocate(character(len=k) :: parentInpFile)
parentInpFile = trim(filenameInp)

allocate( actcoeff_ratio(nn), actcoeffRef1(nn), actcoeffRef2(nn), MolarMass(nn), ntmoles(nn), &
    & OtoC(nn), HtoC(nn), NtoC(nn), StoC(nn), meanOS_C(nn), TotalMassConc(nn) )

OtoC = -7.777777_wp      !set to an impossible value at initialization
HtoC = -7.777777_wp
NtoC = -7.777777_wp
StoC = -7.777777_wp
meanOS_C = -7.777777_wp
!
!====  Activity coefficient ratio calculation  =======================================================
!
!(1) Check whether water is already present (then use as reference component 1)
refcp1 = 0
refcp2 = 0
if (waterpresent) then
    refcp1 = findloc(ITAB(1:min(10,nn),16), VALUE=1, dim=1) !usually = 1
    addn = 1
    norg1 = 2
    if (refcp1 /= 1) then
        write(*,'(A,/)') "WARNING: water is present but not the first component (which it should be)! &
            &Incorrect numbering for first organic comp. norg1 will occur!"
        read(*,*)
    endif
else
    addn = 2
    norg1 = 1
endif

!(2) Set or introduce the two reference solvents to be used as neutral system components 
!    (at the end of current neutral components, but prior to electrolyte components of the input system).
!    In this case, re-package data for transfer, since the number of input component and array sizes will have changed:
allocate( compIDdat(nn+addn), cpsubgdat(nn+addn,topsubno), cpname(nn+addn) )
compIDdat = 0
cpsubgdat = 0
cpname(1:nn) = cpnameinp(1:nn)
compIDdat(1:nn) = CompN(1:nn)
do i = 1,200
    cpsubgdat(1:nn,i) = ITAB(1:nn,i)
enddo
if (.not. waterpresent) then
    refcp1 = nn+1
    cpname(nn+1) = "H2O"
    compIDdat(nn+1) = 401  !add water as new neutral component after present neutral components
    cpsubgdat(nn+1,16) = 1 !add water in form of its subgroup
endif
refcp2 = nn+addn
cpname(nn+addn) = "1,2-Hexanediol"
compIDdat(nn+addn) = 618   !1,2-Hexanediol as ref. component 2
cpsubgdat(nn+addn,145) = 1 !add 1,2-Hexanediol subgroups
cpsubgdat(nn+addn,146) = 3
cpsubgdat(nn+addn,150) = 1
cpsubgdat(nn+addn,151) = 1
cpsubgdat(nn+addn,153) = 2

cpnameRefcp(1) = cpname(refcp1)
cpnameRefcp(2) = cpname(refcp2)

write(*,*) ""
write(*,'(A,/)') "running AIOMFAC_calc to determine the activity coeff. ratios"
!(3) Re-define mixture properties and component properties of this system, including mapping of cpsubgdat to ITAB.
!    First brake the system into chunks of no more than 'nsubset' components.
nsubset = 5     !max. number of components in a subsystem calculation for activity coeff. ratios; (small nsubset < 10 seems best).
                !This is important for speeding up calculations with systems of thousands of components:
deallocate(cpname)
nend = 0
do !nsubset loop
    
    nstart = nend + 1
    nend = min(nend + nsubset, nn)
    nTemp = nend -nstart + 1
    if (allocated(compIDdatTemp)) then
        deallocate(compIDdatTemp, cpsubgdatTemp, cpname)
    endif
    if (nstart > nn) then                   !test for exit condition
        allocate(cpname(nn+addn))           !re-allocate this array as it had been changed
        cpname(1:nn) = cpnameinp(1:nn)
        exit                                !done with nsubset loop
    endif
    !transfer component data for this chunk of system components:
    allocate( compIDdatTemp(nTemp+2), cpsubgdatTemp(nTemp+2,topsubno), cpname(nTemp+2) ) !the +2 are for the two reference components
    compIDdatTemp = 0
    cpsubgdatTemp = 0
    cpname(1:nTemp) = cpnameinp(nstart:nend)
    cpname(nTemp+1) = cpnameRefcp(1)
    cpname(nTemp+2) = cpnameRefcp(2)
    compIDdatTemp(1:nTemp) = compIDdat(nstart:nend)
    compIDdatTemp(nTemp+1) = 401            !water
    compIDdatTemp(nTemp+2) = compIDdat(nn+addn)
    do i = 1,200
        cpsubgdatTemp(1:nTemp,i) = cpsubgdat(nstart:nend,i)
    enddo
    cpsubgdatTemp(nTemp+1,1:200) = cpsubgdat(refcp1,1:200)
    cpsubgdatTemp(nTemp+2,1:200) = cpsubgdat(refcp2,1:200)
    !now re-define subsystem composition for this data chunk:
    call definemixtures(1, nTemp+2, compIDdatTemp, cpsubgdatTemp) 

    !(4) Loop over all system components, each individually dissolved in water as reference solvent and component 1.
    wtftestcp = 0.01_wp
    wtf = 0.0_wp
    do nc = 1,nTemp
        if (nc > 1) then
            wtf(nc-1) = 0.0_wp                      !reset prior component's value
        endif
        wtf(nTemp+1) = 1.0_wp - wtftestcp           !which is the refcp1
        wtf(nc) = wtf(nc) + wtftestcp
        call AIOMFAC_calc(wtf, Tref)                !calculate at given mass fraction and reference temperature
        actcoeffRef1(nstart+nc-1) = actcoeff_n(nc)  !save value at 1:nn index position
    enddo
    
    !(5) Loop over all system components, each individually dissolved in reference solvent 2.
    wtftestcp = 0.01_wp
    wtf = 0.0_wp
    do nc = 1,nTemp
        if (nc > 1) then
            wtf(nc-1) = 0.0_wp                      !reset prior component's value
        endif
        wtf(nTemp+2) = 1.0_wp - wtftestcp           !which is the refcp2
        wtf(nc) = wtf(nc) + wtftestcp
        call AIOMFAC_calc(wtf, Tref)                !calculate at given mass fraction and reference temperature
        ncomp = nstart+nc-1
        actcoeffRef2(ncomp) = actcoeff_n(nc)
        ! also retrieve O:C, H:C, N:C ratios of components for output:
        call O2C_H2C_component(nc, cpCarbon, cpHydrogen, cpOxygen, cpNitrogen, cpSulfur, OtoC(ncomp), HtoC(ncomp), NtoC(ncomp), StoC(ncomp))
        MolarMass(ncomp) = Mmass(nc)
    enddo
    
enddo !nsubset loop

!(6) Compute activity coefficient ratio for all neutral system components as (value in ref. cp. 2) / (value in ref.cp 1):
actcoeff_ratio = actcoeffRef2 / actcoeffRef1
write(*,'(A,/)') "done with activity coeff. ratio calculation"
deallocate(actcoeffRef1, actcoeffRef2)

!
!====  Volatility axis parameters and input/output of other properties  =====================
!
allocate( EVAP_paramA(nn), EVAP_paramB(nn), EVAP_psat(nn), psat(nn), SMILES(nn) ) 
!Set file path and name properties:
k = index(fpathout, '/Output')
fpathin = fpathout(1:k-1)//'/Input'//trim(fpathout(k+7:))
filename = trim(filenameInp)
i = index(filename, ".txt")                         !returns starting position ofstring input within string filename
inpnum = trim(filename(i-4:i))
read(inpnum,*) ndi                                  !save inpnum also as integer ndi
if (change_nd) then 
    ndiout = 1260                                   !set numbering of output files starting with 1260 regardless of current input file number (for use in AIOMFAC G2P calc)
else
    ndiout = ndi
endif
write(outpnum,'(I4.4)') ndiout                      !convert integer to character string
!--
!Reading pure-component vapor pressure input file for getting EVAPORATION parameters A, B or equivalent parameters for the selected VP method:
call ReadCompVaporPressures(nn, TKelvin, refcp1, vpMethod, EVAP_paramA, EVAP_paramB, EVAP_psat, SMILES)
!
!====  Read total molar input concentrations [mol/m^3 of air] of the components ======
!
call ReadMolarComposition(nn, ntmoles) 

do ncp = 1,nn
    !write(*,*) "Molar mass of component ", ncp, MolarMass(ncp)     !molar mass here in [kg/mol]
    TotalMassConc(ncp) = (ntmoles(ncp)*1.0E9_wp) *MolarMass(ncp)    ![mol/m^3] to [ug/m^3]
enddo

psat = EVAP_psat
psatmethod = trim(methodNames(vpMethod))

deallocate( EVAP_psat ) 

!estimate the mean oxidation state of carbon; see supplement of Kroll et al. (2011, doi:10.1038/NCHEM.948)
meanOS_C(norg1:nn) = 2*OtoC(norg1:nn) - HtoC(norg1:nn) - 5*NtoC(norg1:nn) 

!## Pre-lumping of identical components ##
!Determine if two or more components have the exact same SMILES and vapour pressure and are thus identical on that level of structure resolution.
!If it is the case, lump those component mass concentrations together into the component of lower component number;
!This helps with correct displaying of the actual mass concentrations of distinct components in plots.
!--
write(*,'(A,/)') "optional pre-lumping of duplicate SMILES"
!first sort vapour pressure array:
allocate (psatCopy(nn))
psatCopy = psat
call QsortC(psatCopy)                                                   !sort array

!determine if there are identical psat values in array, indicating the same component:
addi = 0
duplicates = .true.
do while(duplicates)
    duplicates = .false.
    do ncp = norg1,nn-1
        if (abs(psatCopy(ncp) - psatCopy(ncp+1)) < ntiny) then          !found the same component
            nc = findloc(psat(:), VALUE = psatCopy(ncp), dim=1)         !find array index of unsorted array
            k =  findloc(psat(:), VALUE = psatCopy(ncp+1), BACK = .true., dim=1)
            !check to confirm that they are duplicates indeed:
            if (trim(SMILES(nc)) == trim(SMILES(k))) then               !confirmed; transfer mass conc. and moles to the first component
                if (ntmoles(k) > 0.0_wp) then
                    duplicates = .true.
                    ntmoles(nc) = ntmoles(nc) + ntmoles(k)
                    ntmoles(k) = 0.0_wp
                    TotalMassConc(nc) = TotalMassConc(nc) + TotalMassConc(k)
                    TotalMassConc(k) = 0.0_wp
                    addi = addi + 1                                     !count number of duplicates
                endif
            endif
        endif
    enddo
enddo
write(*,'(A,/)') "done with pre-lumping of duplicate SMILES"

deallocate(ntmoles)
!
!==== OUTPUT data from benchmark system to files ========================
!
filename = "input_"//outpnum//".txt"
filein = trim(fpathin)//trim(parentInpFile)
fileout = trim(fpathout)//trim(filename)
call copy_file(filein, fileout)

!potentially add electrolyte component information to generated input file:
call Add_Electrolyte_Component(filename, elcpname, eladdMassConc)

select case(YaxisChoice) 
case(1)
    chyaxis = 'log10(act.coeff_ratio)'
case(2)
    chyaxis = 'O:C_ratio'
case(3)  
    chyaxis = 'log10(O:C_ratio)'
case(4)
    chyaxis = 'meanOS_C'
end select

Lmethod = "FullSystem"
resol = "fullRes"

filename = "input_"//outpnum//"_EVAP_AB.txt"                            !e.g. input_1260_EVAP_AB
call OutputVaporPressureParam(nn, norg1, filename, Lmethod, EVAP_paramA, EVAP_paramB, SMILES, TotalMassConc)

filename = "SystemCompProp_"//outpnum//'_'//trim(Lmethod)//'.txt'
call OutputComponentProperties(nn, norg1, TKelvin, filename, Lmethod, psatmethod, resol, psat, OtoC, HtoC, &
                              & meanOS_C, actcoeff_ratio, MolarMass, TotalMassConc)  

filename = "LumpedConc_"//outpnum//'_'//trim(Lmethod)//'.txt'
call OutputLumpedConcentrations(nn, norg1, filename, Lmethod, psatmethod, resol, TotalMassConc, TotalMassConc, MolarMass)

filename = "input_"//outpnum//'_aw_inp.txt'
call Output_aw_levels(filename, awlevels)

filename = "input_"//outpnum//'_Conc_and_T.txt'
call OutputConc_and_T(nn, norg1, linemax, filename, Lmethod, TKelvin, TotalMassConc, MolarMass, elcpname, eladdMassConc)

filename = "input_"//outpnum//'_SMILES.txt'                             !save used SMILES file of the full system also with the input_ prefix
call OutputLumpedSMILES(nn, norg1, filename, TotalMassConc, SMILES)

!
!====  Run the different lumping schemes with selected 2D axes ===========================
!
write(*,*)
write(*,'(A)') "running the various gridded lumping schemes"
call Lumping_schemes(nn, norg1, psatmethod, TKelvin, psat, OtoC, HtoC, meanOS_C, &
        actcoeff_ratio, MolarMass, TotalMassConc, EVAP_paramA, EVAP_paramB, SMILES)

deallocate( cpname, parentInpFile, actcoeff_ratio, EVAP_paramA, EVAP_paramB, OtoC, HtoC, psat, SMILES, awlevels, MolarMass, TotalMassConc ) 

write(*,*)
write(*,'(A,/)') "completed lumping system components and mass concentrations to sets of surrogate components"
write(*,'(A,/)') ">> see output in folder 'Output_lumping' <<"
write(*,*)
write(*,'(A)') "press Enter to continue"
read(*,*)       !wait for user action

end subroutine ActCoeffRatio_Volatility

end submodule SubModactcoeffratio
! ======================= end =======================================================