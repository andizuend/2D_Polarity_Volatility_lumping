submodule (ModLumping) SubModLumpingSchemes
    
contains
    
!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Main subroutine to test different surrogate component selection approaches by      *
!*   of which involve a form of lumping/clustering of component mass concentrations.    *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Andreas Zuend, Dalrin Ampritta Amaladhasan                                         *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2019-03-07                                                      *
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
module subroutine Lumping_schemes(nn, norg1, psatmethod, TK, psat, OtoC, HtoC, meanOS_C, &
                    & actcoeff_ratio, MolarMass, MassConc1, EVAP_paramA, EVAP_paramB, SMILES)

use Mod_kind_param, only : wp
use ModInOutLumping, only : ReadCompVaporPressures, ReadMolarComposition, Add_Electrolyte_Component, &
    & OutputComponentProperties, OutputLumpedConcentrations, OutputLumpedSMILES, Output_aw_levels,   &
    & OutputConc_and_T, OutputVaporPressureParam, OutputGridlines, OutputKmeansCluster,              &
    & SMILES_to_AIOMFAC_inputFile, awlevels, ndiout, linemax  
use ModGridCellSurrogate, only : aspectRatio_XtoY, magnX, magnY, KmeansClusteringSurrogates, &
    & LocGridCellSurrogate_Medoid, LocGridCellSurrogate_MassWeighted_Medoid, LocGridCellSurrogate_Midpoint

implicit none
!-- interface variables:
integer,intent(in) :: norg1, nn  !nn is the number of neutral components
character(len=*),intent(in) :: psatmethod
real(wp),intent(in) :: TK
real(wp),dimension(:),intent(in) :: psat, OtoC, HtoC, meanOS_C, actcoeff_ratio, &
                                & MolarMass, MassConc1, EVAP_paramA, EVAP_paramB
character(len=*),dimension(:),intent(in) :: SMILES
!-- local variables:
character(len=2) :: chx, chy
character(len=4) :: chfnum
character(len=30) :: Lmethod, resol
character(len=1000) :: filename, fname2
!--
integer :: allocstat, clust_no, dim_num, i, j, k, m, nc, ns, nschemes, fnum, scID, scIDHighVolat, schemeID, seed
integer,dimension(nn) :: gridcomp
integer,dimension(:,:,:),allocatable :: grid
integer,dimension(:,:),allocatable :: midcomp, counter
integer,dimension(:),allocatable :: index_k, cluster, cluster_population, cluster_surrogate, &
    & cluster_k, cluster_population_k, cluster_surrogate_k
!--
real(wp),parameter :: tolerance = 1.0E-2_wp*sqrt(epsilon(1.0_wp))
real(wp) :: intervalx, intervaly, midpx, midpy, lylim, uylim, uxlim, lxlim, rno, rno2, &
    & MassConcHighVolatSurrogate, sumx
real(wp),dimension(nn) :: MassConc, MassConc2, MassConc3, MassConc4, MassConc5, MassConcLumped
real(wp),dimension(nn) :: Xaxis, Yaxis, lgpsat
real(wp),dimension(:,:),allocatable :: cluster_center, cluster_center_k
real(wp),dimension(:),allocatable :: xgridl, ygridl, Xax_k, Yax_k, MassConc_k, MassConc5_k
!--
logical :: modifyClusterNum
logical,dimension(:),allocatable :: is_highVolatComp
!...........................................

!--- parameter settings (now defined in SETTINGS input file)
if (cluster_num < 0) then
    modifyClusterNum = .true.
    cluster_num = nxintervals*nyintervals   !set as an input for the k-means method to yield similar numbers of surrogates as the other methods or set further below.
else
    modifyClusterNum = .false.
endif
!define the aspect ratio of normalized x-axis range relative to y-axis range:
aspectRatio_XtoY = 2.0_wp / (real(nxintervals, kind=wp) / real(nyintervals, kind=wp))    
    !That is, the (prescribed) multiplying factor of the normalized, dimensionless value range along the x-axis that is regarded as equal to 
    !the normalized, dimensionless range along the y-axis. For normalization, divide the values along the x-axis by the x-axis range of interest 
    !(max(x) - min(x)) and analogously for the y-axis when distance measures are involved, e.g. in medoid lumping, but we also multiply the x-axis 
    !distance values by aspectRatio_XtoY to achieve selected weighting regardless of the absolute value ranges. 
!---

!Set axis and ranges of 2-D coordinate grid for lumping schemes:
lgpsat = log10(psat)
Xaxis = lgpsat
select case(YaxisChoice) 
case(1) !log10(activity coeff. ratio)
    Yaxis = log10(actcoeff_ratio)
case(2) !(O:C ratio)
    Yaxis = OtoC   
case(3) !log10(O:C ratio)
    Yaxis = log10(max(OtoC, 1.0E-2_wp))                 !if O:C is zero, it will use 0.01
case(4) !mean oxidation state of carbon (estimated)
    Yaxis = meanOS_C 
end select

MassConc = MassConc1

!Lump conc. of all components whose vp is > a threshold log10(psat) value to one special surrogate component.
allocate(is_highVolatComp(nn))
is_highVolatComp = .false.
gridcomp = 0
m = 0
do i = norg1,nn
    if (lgpsat(i) > SetVolatThresholdVal) then  !SetVolatThresholdVal is set in Settings_2DLumping.txt file
        is_highVolatComp(i) = .true.
        m = m + 1
        gridcomp(m) = i
    endif
enddo

!find the extreme axis values among the components of non-zero mass conc.:
uxlim = maxval(Xaxis(norg1:), mask = MassConc(norg1:) > 0.0_wp .and. .not. is_highVolatComp(norg1:) )
lxlim = minval(Xaxis(norg1:), mask = MassConc(norg1:) > 0.0_wp .and. .not. is_highVolatComp(norg1:) )
uylim = maxval(Yaxis(norg1:), mask = MassConc(norg1:) > 0.0_wp .and. .not. is_highVolatComp(norg1:) )
lylim = minval(Yaxis(norg1:), mask = MassConc(norg1:) > 0.0_wp .and. .not. is_highVolatComp(norg1:) )

magnX = abs(uxlim - lxlim)  !save magnitude of axis value ranges for normalization in ModGridCellSurrogate
magnY = abs(uylim - lylim)

!Determine highly volatile components and lump them into a special surrogate based on the weighted medoid method
if (m > 0) then
    call LocGridCellSurrogate_MassWeighted_Medoid(m, Xaxis, Yaxis, scIDHighVolat, gridcomp, MassConc)
    MassConcHighVolatSurrogate = sum(MassConc(gridcomp(1:m)))
    MassConc(gridcomp(1:m)) = 0.0_wp
    MassConc(scIDHighVolat) = MassConc1(scIDHighVolat)      !keep original value for the designated high-vol surrogate
else
    scIDHighVolat = -1                                      !indicating that no high-volatility surrogate exists;
    MassConcHighVolatSurrogate = 0.0_wp
endif

!check whether set grid resolution is possible:
if (nxintervals*nyintervals > nn) then 
    write(*,*) ""
    write(*,'(A)') "ERROR detected in Lumping_schemes: the chosen grid resolution is higher than the number of system components! & 
                    & This is incorrect since it is not possible to assign more surrogate components than components. &
                    & Therefore, check the set 'nxintervals' and 'nyintervals' values in the 'SETTINGS_2DLumping.txt' file."  
    read(*,*) 
endif

allocate( counter(nxintervals,nyintervals), midcomp(nxintervals,nyintervals), grid(nxintervals,nyintervals,nn), &
    & xgridl(nxintervals+1), ygridl(nyintervals+1), STAT=allocstat)

grid = 0
intervalx = 1.001_wp*magnX / real(nxintervals, kind=wp)
intervaly = 1.001_wp*magnY / real(nyintervals, kind=wp)
counter = 0
gridcomp = 0

!assign x-axis and y-axis grid line coordinates as vectors with regular (equal) grid line spacings:
xgridl = [(lxlim +intervalx*j, j = 0,nxintervals)]
ygridl = [(lylim +intervaly*j, j = 0,nyintervals)]
!assign the components to grid cells within the selected coordinate system:
do i = norg1,nn
    if (.not. is_highVolatComp(i)) then
        j = minloc(Xaxis(i) - xgridl(:), mask = Xaxis(i) - xgridl(:) > 0.0_wp, dim=1) 
        k = minloc(Yaxis(i) - ygridl(:), mask = Yaxis(i) - ygridl(:) > 0.0_wp, dim=1)
        if (MassConc(i) <= 0.0_wp) then  !the grid indices might be out of bounds, so check:
            j = max(min(j, nxintervals), 1) 
            k = max(min(k, nyintervals), 1) 
        endif
        grid(j,k,i) = i
        counter(j,k) = counter(j,k)+1
    endif
enddo

!use the three different approaches to determine surrogates in each grid cell:
MassConc2 = 0.0_wp
MassConc3 = 0.0_wp
MassConc4 = 0.0_wp
do j = 1, nxintervals
    do k = 1,nyintervals
        m = 0
        gridcomp = 0
        do i = norg1, nn
            if (grid(j,k,i) > 0 .and. counter(j,k) > 1) then
                m = m+1
                gridcomp(m) = i
            else if (grid(j,k,i) > 0 .and. counter(j,k) == 1) then
                scID = grid(j,k,i)
            endif
        enddo
        if (m > 1) then
            do ns = 1,3
                select case(ns)
                case(1)             !medoid (unweighted) output of lumped surrogate system
                    call LocGridCellSurrogate_medoid(m, Xaxis, Yaxis, scID, gridcomp)
                    MassConc2(scID) = sum( MassConc(gridcomp(1:m)) )
                case(2)             !midpoint (unweighted) output of lumped surrogate system
                    midpx = 0.5_wp*( xgridl(j) + xgridl(j+1) )
                    midpy = 0.5_wp*( ygridl(k) + ygridl(k+1) )
                    call LocGridCellSurrogate_Midpoint(m, Xaxis, Yaxis, midpx, midpy, scID, gridcomp)
                    MassConc3(scID) = sum( MassConc(gridcomp(1:m)) )
                case(3)             !mass-weighted medoid
                    call LocGridCellSurrogate_MassWeighted_Medoid(m, Xaxis, Yaxis, scID, gridcomp, MassConc)
                    MassConc4(scID) = sum( MassConc(gridcomp(1:m)) )
                end select
            enddo
        else if (counter(j,k) == 1) then
            MassConc2(scID) = MassConc(scID)
            MassConc3(scID) = MassConc(scID)
            MassConc4(scID) = MassConc(scID)
        endif
    enddo
enddo

if (scIDHighVolat > 0) then
    MassConc2(scIDHighVolat) = MassConcHighVolatSurrogate 
    MassConc3(scIDHighVolat) = MassConcHighVolatSurrogate 
    MassConc4(scIDHighVolat) = MassConcHighVolatSurrogate
endif

!
!--- k-means clustering method:
!
if (modifyClusterNum) then !set cluster_num to the same value as the number of surrogate components used with the gridded methods
    cluster_num = count(MassConc4(:) > 0.0_wp)
    if (scIDHighVolat > 0) then
        cluster_num = cluster_num - 1
    endif
endif
if (nn < cluster_num) then
    write(*,'(A)') "ERROR from Lumping_schemes: the set number of clusters (cluster_num) > number of system components!"
    read(*,*)
endif
dim_num = 2
allocate(cluster_center_k(dim_num, cluster_num))
!initialize the cluster_centers with the scaled aspect ratio and axis magnitudes from grid cell (midpoint) coordinates (where possible):
j = 0
do i = 1,nxintervals
    do k = 1,nyintervals
        j = j+1
        if (j > cluster_num) then
            exit
        endif
        cluster_center_k(1,j) = 0.5_wp*(xgridl(i) + xgridl(i+1))*aspectRatio_XtoY / magnX       !scaled x-coord
        cluster_center_k(2,j) = 0.5_wp*(ygridl(k) + ygridl(k+1)) / magnY                        !scaled y-coord
    enddo
enddo
if (j < cluster_num) then               !if needed, add extra cluster centers
    seed = 11314177
    call random_seed(seed)
    do j = j+1,cluster_num
        call random_number(rno)
        rno = 0.05_wp + rno*0.9_wp        !limit the interval to [0.05, 0.95]
        rno2 = 1.0_wp - rno
        cluster_center_k(1,j) = ( rno*xgridl(1) + rno2*xgridl(nxintervals) )*aspectRatio_XtoY / magnX      !scaled random x-coord
        cluster_center_k(2,j) = ( rno*ygridl(1) + rno2*ygridl(nyintervals) ) / magnY                       !scaled random y-coord
    enddo
endif

!filter out components that are already associated with the special high-volatility surrogate "cluster".
!those should not be part of the kmeans clustering:
nc = count(.not. is_highVolatComp(norg1:))                  !the number of organic components to be clustered
allocate(Xax_k(nc), Yax_k(nc), MassConc_k(nc), MassConc5_k(nc), index_k(nc), stat = allocstat)
Xax_k = huge(1.0_wp)
Yax_k = huge(1.0_wp)
MassConc_k = 0.0_wp
MassConc5 = 0.0_wp
j = 0
do i = norg1,nn
    if (.not. is_highVolatComp(i)) then
        j = j + 1
        index_k(j) = i                                  !mapping index to associate cluster set component number to the component number in the total set
        Xax_k(j) = Xaxis(i)
        Yax_k(j) = Yaxis(i)
        MassConc_k(j) = MassConc(i)
    endif
enddo

!now call the KmeansClusteringSurrogates subroutine from ModGridCellSurrogate:
call KmeansClusteringSurrogates(nc, 1, cluster_num, dim_num, Xax_k, Yax_k, MassConc_k, cluster_k, cluster_center_k, &
                             & cluster_population_k, cluster_surrogate_k, MassConc5_k)

!map back from cluster index to the full-set component numbering:
if (scIDHighVolat > 0) then
    clust_no = cluster_num + 1
else
    clust_no = cluster_num
endif
allocate( cluster(norg1:nn), cluster_population(clust_no), cluster_surrogate(clust_no), &
    & cluster_center(dim_num, clust_no), STAT=allocstat )   !the +1 is to accommondate the high-volatility surrogate as extra cluster

cluster = 0

!also add high-volatility compound as a surrogate and special cluster population (for plots)
if (scIDHighVolat > 0) then
    MassConc5(scIDHighVolat) = MassConcHighVolatSurrogate
    cluster_population(clust_no) = count(is_highVolatComp(:))
    cluster_surrogate(clust_no) = scIDHighVolat
    cluster(:) = clust_no                                !assign all components the high volatility cluster by default (changed for non-high-volatility components below)
    cluster_center(1:2,clust_no) = [Xaxis(scIDHighVolat), Yaxis(scIDHighVolat)]
endif
!map from kmeans arrays to full set arrays:
MassConc5(index_k(1:nc)) = MassConc5_k(1:nc)
cluster(index_k(1:nc)) = cluster_k(1:nc)
cluster_center(1:2,1:cluster_num) = cluster_center_k(1:2,1:cluster_num)
cluster_population(1:cluster_num) = index_k(cluster_population_k(1:cluster_num))
cluster_surrogate(1:cluster_num) = index_k(cluster_surrogate_k(1:cluster_num))

!check mass conservation:
sumx = sum(MassConc1(norg1:nn))
if (abs(sumx - sum(MassConc5(norg1:nn))) / sqrt(sumx) > tolerance) then
    write(*,*) ""
    write(*,'(A,*(ES12.4,:,2X))') "Total new comp mass conc from kmeans approach is NOT equal to total old comp composition", sum(MassConc5(norg1:nn)), sum(MassConc(norg1:nn))
    read(*,*)
endif

deallocate(cluster_surrogate_k, cluster_k, cluster_center_k, cluster_population_k, STAT=allocstat)


!------ OUTPUT data to files ------------------------------------------------------------
!
!save grid resolution used:
write(chx,'(I2.2)') nxintervals
write(chy,'(I2.2)') nyintervals
resol = trim(chx)//"x"//trim(chy)
!output gridline coordinates file:
write(chfnum,'(I4.4)') ndiout
filename = "GridlineCoord_"//chfnum//'_'//trim(resol)//'.txt'
call OutputGridlines(xgridl, ygridl, filename)

schemeID = 0 !initialize
nschemes = 4  !number of different schemes evaluated above (for now = 3)
do ns = 1,nschemes
    select case(ns)
    case(1) !Medoid (unweighted) output of lumped surrogate system
        Lmethod = "Medoid"
        schemeID = 1
        MassConcLumped = MassConc2
    case(2) !Midpoint (unweighted) output of lumped surrogate system
        Lmethod = "Midpoint"
        schemeID = 2
        MassConcLumped = MassConc3
    case(3)
        Lmethod = "Weighted_Medoid"
        schemeID = 3
        MassConcLumped = MassConc4
    case(4)
        Lmethod = "Kmeans"
        schemeID = 4
        MassConcLumped = MassConc5
    case default
        write(*,*) "WARNING from Lumping_schemes: no scheme case defined! ns = ", ns
        read(*,*)
    end select
    
    !== call the output subroutines to generate files in 'Output_lumping' folder ==
    if (change_nd) then
        fnum = ndiout + schemeID        !e.g. 1261
    else
        fnum = ndiout                   !e.g. 1403 
    endif
    write(chfnum,'(I4.4)') fnum
    filename = "SystemCompProp_"//chfnum//'_'//trim(Lmethod)//'_'//trim(resol)//'.txt'
    call OutputComponentProperties(nn, norg1, TK, filename, Lmethod, psatmethod, resol, psat, OtoC, HtoC, &
                                   & meanOS_C, actcoeff_ratio, MolarMass, MassConcLumped)
    !output lumped and original (total) mass concentrations and molar amounts:
    filename = "LumpedConc_"//chfnum//'_'//trim(Lmethod)//'_'//trim(resol)//'.txt'
    call OutputLumpedConcentrations(nn, norg1, filename, Lmethod, psatmethod, resol, MassConc, MassConcLumped, MolarMass) 
    !special output file with additional data for k-means cluster populations:
    if (schemeID == 4) then
        filename = "KmeansClusters_"//chfnum//'_'//trim(resol)//'.txt'
        call OutputKmeansCluster(nn, norg1, dim_num, clust_no, cluster, cluster_center, cluster_population, cluster_surrogate, filename)
    endif
    !----
    !output files used as inputs directly with AIOMFAC gas-particle partitioning calculation:
    if (change_nd) then                                                         !output files used for input in Gas2Part calc:
        filename = "input_"//chfnum//'_EVAP_AB.txt'  
        call OutputVaporPressureParam(nn, norg1, filename, Lmethod, EVAP_paramA, EVAP_paramB, SMILES, MassConcLumped)
        filename = "input_"//chfnum//'_SMILES.txt'  !!"smiles_"//chfnum//'.txt'
        call OutputLumpedSMILES(nn, norg1, filename, MassConcLumped, SMILES)
        !generate an AIOMFAC-web style input file using the SMILES list of surrogate components only (for use with AIOMFAC gas--particle partitioning calc.):
        fname2 = filename                                                       !name of generated SMILES file of surrogate compounds;
        filename = "input_"//chfnum//".txt"
        call SMILES_to_AIOMFAC_inputFile(fname2, filename)                      !needs access and uses external Python code (Tools_for_SMILES_conversion);
        call Add_Electrolyte_Component(filename, elcpname, eladdMassConc)   !potentially add an electrolyte component to the generated input file;
        filename = "input_"//chfnum//'_aw_inp.txt'
        call Output_aw_levels(filename, awlevels)
        filename = "input_"//chfnum//'_Conc_and_T.txt'
        call OutputConc_and_T(nn, norg1, linemax, filename, Lmethod, TK, MassConcLumped, MolarMass, elcpname, eladdMassConc)
        call SLEEP(1)
    endif
    !----
enddo !nschemes
!----------- end of output ------------------------------------------------------------------------------


deallocate(grid, counter, midcomp, xgridl, ygridl, STAT=allocstat)
deallocate(cluster_surrogate, cluster, cluster_center, cluster_population, STAT=allocstat)

end subroutine Lumping_schemes
!====================================================================================================================
    
end submodule SubModLumpingSchemes