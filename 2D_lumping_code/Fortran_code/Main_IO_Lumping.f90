!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Main program for the 2D polarity--volatility lumping framework based on AIOMFAC    *
!*   activity coefficient ratios or other polarity proxies like O:C ratios.             *
!*   This main program interfaces with the AIOMFAC core model with its input from an    * 
!*   AIOMFAC-web style input file provided in folder Input_lumping.                     *
!*   The SETTINGS_2DLumping.txt file contains all main parameters for a specific        *
!*   lumping framework run.                                                             *
!*   Output data files are directed to the Output_lumping folder.                       *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Andi Zuend, Dalrin Ampritta Amaladhasan, Dan Hassan-Barthaux                       *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University (2013 - present)         *
!*                                                                                      *
!*   -> created:        2019  (this file)                                               *
!*   -> latest changes: 2025-09-18                                                      *
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
    
program Main_IO_Lumping

!module variables:
use Mod_kind_param, only : wp
use ModSystemProp, only : errorflag_clist, errorflagmix, nneutral, SetSystem, topsubno
use ModSubgroupProp, only : SubgroupAtoms, SubgroupNames
use ModMRpart, only : MRdata
use ModSRparam, only : SRdata
use Mod_InputOutput, only : RepErrorWarning, ReadInputFile
use ModLumping, only : ActCoeffRatio_Volatility, change_nd, cluster_num, eladdMassConc, elcpname, &
    & massweightedKmeans, nxintervals, nyintervals, onlySurrogatesOutput, SetVolatThresholdVal, vpMethod, YaxisChoice

implicit none
!set preliminary AIOMFAC-input-related parameters:
integer,parameter :: maxpoints = 2                          !limit maximum number of composition points; here set to only 2 for activity ratio computation
integer,parameter :: ninpmax = int(1E+5)                    !set the maximum number of mixture components allowed (preliminary parameter)
!local variables:
character(len=4) :: version_2DPVLF, version_AIOMFAC
character(len=200) :: filename, text
character(len=3000) :: filepath, folderpathout, txtfilein  
character(len=200),dimension(:),allocatable :: cpnameinp    !list of assigned component names (from input file)
integer :: allocstat, errorind, i, istat, ncp, nn, npoints, &
    & nspecmax, unito, warningflag, warningind, unt
integer,dimension(:,:),allocatable :: cpsubg                !list of input component subgroups and corresponding subgroup quantities
real(wp) :: TKelvin
real(wp),dimension(:),allocatable :: T_K
real(wp),dimension(:,:),allocatable :: composition
logical :: filevalid, verbose, xinputtype
logical,dimension(size(errorflag_clist)) :: errflag_list
!...................................................................................

!
!==== INITIALIZATION section =======================================================
!
version_2DPVLF  = "1.0"     !version of the 2D Polarity-Volatility lumping framework program
version_AIOMFAC = "3.12"    !AIOMFAC-web version of core AIOMFAC code used in this project
verbose = .true.            !if true, some debugging information will be printed to the unit "unito" (errorlog file)
nspecmax = 0
errorind = 0                !0 means no error found
warningind = 0              !0 means no warnings found
txtfilein = ''
!
!==== INPUT data section ===========================================================
!
txtfilein = '.././Input_lumping/SETTINGS_2DLumping.txt'
!---
filepath = adjustl(trim(txtfilein))
write(*,*) ""
write(*,'(A)')   "MESSAGE from 2D_Polarity-Volatility_lumping: program started"
write(*,'(A,A)') "parameter settings file = ", trim(filepath)
write(*,*) ""
!---
!read SETTINGS_2DLumping file and set related model parameters:
open (NEWUNIT = unt, FILE = filepath, IOSTAT=istat, ACTION='read', STATUS='OLD')
read(unt,*) text !read over irrelevant text...
read(unt,*) text
read(unt,*) 
read(unt,*) text
read(unt,*) text
read(unt,*) filename
read(unt,*) 
read(unt,*) TKelvin
read(unt,*) elcpname
read(unt,*) eladdMassConc
read(unt,*) 
read(unt,*) change_nd
read(unt,*) vpMethod
read(unt,*) onlySurrogatesOutput
read(unt,*) massweightedKmeans
read(unt,*) 
read(unt,*) YaxisChoice
read(unt,*) nxintervals
read(unt,*) nyintervals
read(unt,*) cluster_num
read(unt,*) 
read(unt,*) SetVolatThresholdVal
close(unt)
!---

allocate(cpsubg(ninpmax,topsubno), cpnameinp(ninpmax), composition(maxpoints,ninpmax), T_K(maxpoints), STAT=allocstat)

!-- Read AIOMFAC-web input file to define system components:
filepath = '.././Input_lumping/'//trim(filename)
call ReadInputFile(filepath, folderpathout, filename, ninpmax, maxpoints, unito, verbose, ncp, npoints, &
    & warningind, errorind, filevalid, cpnameinp, cpsubg, T_K, composition, xinputtype)
!--
if (filevalid) then
    !
    !==== AIOMFAC initialization section ===============================
    !
    if (verbose) then
        write(unito,*) ""
        write(unito,'(A,/)') "MESSAGE from AIOMFAC: input file read, starting AIOMFAC mixture definitions and initialization... "
    endif
    
    call MRdata()           !initialize the MR data for the interaction coeff. calculations
    call SRdata()           !initialize data for the SR part coefficient calculations
    call SubgroupNames()    !initialize the subgroup names for the construction of component subgroup strings
    call SubgroupAtoms()
    !--
    !set mixture properties based on the data from the input file:
    call SetSystem(1, .true., ncp, cpnameinp(1:ncp), cpsubg(1:ncp,1:topsubno) )

    !transfer composition data to adequate array size:
    deallocate(cpsubg, composition, STAT=allocstat)
    
    if (errorflagmix /= 0) then
        call RepErrorWarning(unito, errorflagmix, warningflag, errflag_list, i, errorind, warningind)
    endif

    if (errorind == 0) then 
        !
        !@@@@ special AIOMFAC calls to compute activity coefficient ratios of components 
        !     in two reference solvents and subsequently the 2D lumping methods
        !
        nn = nneutral
        call ActCoeffRatio_Volatility(nn, TKelvin, cpnameinp(1:nn), filename, folderpathout)
        !!@@@@
    endif !errorind
else
    write(unito,'(A)')  "ERROR: No valid input file found!"    
    write(*,'(A)')      "ERROR: No valid input file found!"
    read(*,*)
endif !file valid

write(unito,*) "+-+-+-+-+"
write(unito,*) "Final warning indicator (an entry '00' means no warnings found):"
write(unito,'(I2.2)') warningind
write(unito,*) "+-+-+-+-+"
write(unito,*) ""
write(unito,*) "########"
write(unito,*) "Final error indicator (an entry '00' means no errors found):"
write(unito,'(I2.2)') errorind
write(unito,*) "########"
close(unito)    !close the error log-file

!write(*,*) ""
!write(*,*) "MESSAGE from AIOMFAC: End of program; final error indicator: ", errorind
!write(*,*) ""
!read(*,*)  !Pause; just for debugging and testing.
!
!==== The end ======================================================================
!
end program Main_IO_Lumping