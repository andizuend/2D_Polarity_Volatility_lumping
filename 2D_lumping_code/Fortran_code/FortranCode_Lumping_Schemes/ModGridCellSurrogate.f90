!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Module containing subroutines for computing the surrogate component of a selected  *
!*   grid cell by means of different methods.                                           * 
!*                                                                                      *
!*   :: Authors & Copyright ::                                                          *
!*   Ampritta Amaladhasan, Andi Zuend                                                   *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2018                                                            *
!*   -> latest changes: 2025-07-17                                                      *
!*                                                                                      *
!*   :: List of subroutines and functions contained in this module:                     *
!*   --------------------------------------------------------------                     *
!*   -  subroutine LocGridCellSurrogate_Medoid                                          *
!*   -  subroutine LocGridCellSurrogate_MassWeighted_Medoid                             *
!*   -  subroutine LocGridCellSurrogate_Midpoint                                        *
!*   -  subroutine KmeansClusteringSurrogates                                           *
!*                                                                                      *
!****************************************************************************************
module ModGridCellSurrogate

use Mod_kind_param, only : wp

implicit none
!public module variables:
real(wp),public :: aspectRatio_XtoY
real(wp),public :: magnX, magnY                 !the magnitudes of axis value ranges for X and Y axis; used for distance normalization.

!===================
    contains
!===================

    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Subroutine to determine a surrogate component within a selected grid cell based on *
    !*   the principle of finding the 2-D point (component) with minimum squared Eucledian  *
    !*   distance between itself and all other points (components) of the cell. This        *
    !*   provides a reasonable option to deal with any kind of clustering, including cells  *
    !*   in which points are clustered unevenly (e.g. all located near one corner).         *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Ampritta Amaladhasan, Andreas Zuend                                                *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2019-03-07                                                      *
    !*   -> latest changes: 2019-10-23                                                      *
    !*                                                                                      *
    !****************************************************************************************
    pure subroutine LocGridCellSurrogate_Medoid(m, xp, yp, scID, gridcomp)

    implicit none
    !-- interface variables:
    integer,intent(in) :: m                             !number of components within this grid cell
    real(wp),dimension(:),intent(in) :: xp, yp          !x-axis and y-axis coordinates of a set of points (components)
    integer,dimension(:),intent(in) :: gridcomp         !component index number inside the grid (for range 1:m)
    integer,intent(out) :: scID                         !index value of the determined surrogate component for this grid cell, 
                                                        !such that [xp(scID), yp(scID)] are the 2-D coordinates of the surrogate component.
    !-- local variables:
    integer :: i, n
    real(wp),dimension(m) :: sqdist                     !cumulative squared distance of point i from all other points in grid cell
    !.................................................
    
    sqdist = 0.0_wp
    do i = 1,m
        do n = 1,m
            if (i /= n) then    !the calculated squared distances are normalized by the magnitudes of X and Y axis ranges and the given aspect ratio.
                sqdist(i) = sqdist(i) + ( (xp(gridcomp(n)) - xp(gridcomp(i)))*aspectRatio_XtoY / magnX )**2 &
                            & + ( (yp(gridcomp(n)) - yp(gridcomp(i))) / magnY )**2
            endif
        enddo
    enddo
    scID = gridcomp(minloc(sqdist(1:m), dim=1))

    end subroutine LocGridCellSurrogate_Medoid
    !==================================================================================================================================

    

    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Subroutine to determine a surrogate component within a selected grid cell based on *
    !*   the principle of finding the 2-D point (component) with minimum squared Eucledian  *
    !*   distance between itself and the center of mass (concentration) within the cell.    *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andreas Zuend, Ampritta Amaladhasan                                                *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2019-10-29                                                      *
    !*   -> latest changes: 2019-10-29                                                      *
    !*                                                                                      *
    !****************************************************************************************
    pure subroutine LocGridCellSurrogate_MassWeighted_Medoid(m, xp, yp, scID, gridcomp, MassConc) 

    implicit none
    !-- interface variables:
    integer,intent(in) :: m                                 !number of components within this grid cell
    real(wp),dimension(:),intent(in)  :: MassConc           ![ug/m^3] component mass conc.
    real(wp),dimension(:),intent(in) :: xp, yp              !x-axis and y-axis coordinates of a set of points (components)
    integer,dimension(:),intent(in) :: gridcomp             !component index number inside the grid cell
    integer,intent(out) :: scID                             !index value of the determined surrogate component for this grid cell, 
                                                            !such that [xp(scID), yp(scID)] are the coordinates of the surrogate component.
    !-- local variables:
    real(wp),parameter :: hugen = huge(1.0_wp)
    real(wp) :: sumMassConc, massCenterX, massCenterY
    real(wp),dimension(m) :: sqdist                         !squared distance of point i from the center of mass in this grid cell
    !.................................................
    
    !compute coordinates of center of mass:
    sumMassConc = sum( MassConc(gridcomp(1:m)) )
    if (sumMassConc > 0.0_wp) then
        massCenterX = sum( MassConc(gridcomp(1:m))*xp(gridcomp(1:m)) ) / sumMassConc
        massCenterY = sum( MassConc(gridcomp(1:m))*yp(gridcomp(1:m)) ) / sumMassConc
    else !special case
        sqdist = 0.0_wp
        scID = gridcomp(1)                                  !since all components are zero mass, choice of surrogate does not matter.
        return
    endif
    
    !compute squared distance of components (of non-zero mass conc.) from center of mass (including x--y magnitude normalization):
    sqdist = hugen !assign a huge default value to components, so that components of zero mass conc. will not be selectd below
    where (MassConc(gridcomp(1:m)) > 0.0_wp)
        sqdist(1:m) = ( (massCenterX - xp(gridcomp(1:m)))*aspectRatio_XtoY / magnX )**2 + &
                    & ( (massCenterY - yp(gridcomp(1:m))) / magnY )**2
    endwhere
    !select the component ID corresponding to the minimum distance value from center of mass:
    scID = gridcomp(minloc(sqdist(1:m), dim=1))

    end subroutine LocGridCellSurrogate_MassWeighted_Medoid
    !==================================================================================================================================
  
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Subroutine to determine a surrogate component within a selected grid cell based on *
    !*   the principle of finding the component with minimum squared Eucledian              *
    !*   distance between itself the mid point of the grid cell.                            *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Ampritta Amaladhasan, Andreas Zuend                                                *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2019-03-07                                                      *
    !*   -> latest changes: 2020-04-02                                                      *
    !*                                                                                      *
    !****************************************************************************************
    pure subroutine LocGridCellSurrogate_Midpoint(m, xp, yp, midpx, midpy, scID, gridcomp)

    implicit none
    !-- interface arguments:
    integer,intent(in) :: m                             !number of components within this grid cell
    real(wp),dimension(:),intent(in) :: xp, yp          !x-axis and y-axis coordinates of a set of points (components)
    real(wp),intent(in) :: midpx, midpy                 !coordinates of the midpoint of the current grid cell
    integer,dimension(:),intent(in) :: gridcomp         !component index number inside the grid cell
    integer,intent(out) :: scID                         !index value of the determined surrogate component for this grid cell, 
                                                        !such that [xp(scID), yp(scID)] are the coordinates of the surrogate component.
    !-- local variables:
    integer :: i
    real(wp),dimension(m) :: sqdist                     !squared distance of point i from grid cell midpoint
    !.................................................
    
    sqdist = 0.0_wp
    do i = 1,m      !the calculated squared distances of points from the grid cell midpoint are normalized by the magnitudes of X and Y axis ranges and the given aspect ratio.
        sqdist(i) = ( (midpx - xp(gridcomp(i)))*aspectRatio_XtoY / magnX )**2 + ( ( midpy - yp(gridcomp(i))) / magnY )**2
    enddo
    scID = gridcomp(minloc(sqdist(1:m), dim=1))

    end subroutine LocGridCellSurrogate_Midpoint
    !==================================================================================================================================
   
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Subroutine to determine the surrogate components for the given system using a      *
    !*   variant of the k-means clustering algorithm.                                       *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Ampritta Amaladhasan, Andi Zuend                                                   *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2019-03-07                                                      *
    !*   -> latest changes: 2025-07-05                                                      *
    !*                                                                                      *
    !****************************************************************************************   
    subroutine KmeansClusteringSurrogates(nc, n1, cluster_num, dim_num, Xaxis, Yaxis, MassConc_k, cluster_k, cluster_center_k, &
                                        & cluster_population_k, cluster_surrogate_k, MassConcLumped_k)
    use ModLumping, only : massweightedKmeans

    implicit none
    !-- interface arguments:
    integer,intent(in) :: nc                                            !nc is the number of neutral components for clustering
    integer,intent(in) :: n1                                            !n1 is the number of the first organic component in the clustering set (usually 1)
    integer,intent(in) :: cluster_num, dim_num
    real(wp),dimension(:),intent(in) :: Xaxis, Yaxis, MassConc_k        !coordinates and mass conc. of the components for clustering
    integer,dimension(:),allocatable,intent(out) :: cluster_k
    real(wp),dimension(:,:),intent(inout) :: cluster_center_k
    integer,dimension(:),allocatable,intent(out) :: cluster_population_k
    integer,dimension(:),allocatable,intent(out) :: cluster_surrogate_k
    real(wp),dimension(:),intent(out) :: MassConcLumped_k
    !-- local variables:
    integer :: allocstat, i, j, it_max, it_num, npoints !, seed
    real(wp),parameter :: numtiny = tiny(1.0_wp)
    real(wp) :: sqdist_i
    real(wp),allocatable,dimension(:,:) :: point
    real(wp),allocatable,dimension(:) :: cluster_energy, cluster_variance, weights
    real(wp),dimension(cluster_num) :: minSqdistCluster, clusterMassConc
    !....................................

    write(*,* ) ' '
    write(*,'(A,/)') 'running the k-means mass-weighted algorithm'
    
    !parameters:
    npoints = nc - n1 + 1
    it_max = 60
    
    allocate( point(dim_num, n1:nc), weights(npoints) )     !note: second dimension does start at n1, not 1;
    !assign the different properly scaled x- and y-axis coordinates, accounting for aspect ratio and axis scaling:
    point(1,n1:nc) = Xaxis(n1:nc)*aspectRatio_XtoY / magnX
    point(2,n1:nc) = Yaxis(n1:nc) / magnY
    allocate( cluster_k(n1:nc), cluster_energy(cluster_num), cluster_population_k(cluster_num), cluster_variance(cluster_num), &
        & cluster_surrogate_k(cluster_num), STAT=allocstat )

    write(*,'(A,I0,/)') 'number of iterations allowed: ', it_max
    
    if (massweightedKmeans) then                            !the weights of the data used are their mass concentrations
        weights(1:npoints) = max(MassConc_k(n1:nc), numtiny)
    else
        weights(1:npoints) = 1.0_wp                         !use uniform weighting
    endif
    !now call a specific version of the K-means algorithm with weighted points:
    call kmeans_w_01(dim_num, npoints, cluster_num, it_max, it_num, point, weights, cluster_k, cluster_center_k, cluster_population_k, cluster_energy)
    
    write(*,'(A,I0,/)') 'number of k-means iterations taken with 1st kmeans_w_01: ', it_num
    call cluster_variance_compute (dim_num, npoints, cluster_num, point, cluster_k(1:cluster_num), cluster_center_k(1:2,1:cluster_num), cluster_variance(1:cluster_num))
    call cluster_print_summary (npoints, cluster_num, cluster_population_k(1:cluster_num), cluster_energy(1:cluster_num), cluster_variance(1:cluster_num))
    
    !!@@@ for tests:
    !seed = 123456789
    !call cluster_initialize_5 ( dim_num, npoints, cluster_num, point, seed, cluster_center_k )
    !weights(1:npoints) = MassConc_k(n1:nc)
    !call kmeans_w_03(dim_num, npoints, cluster_num, it_max, it_num, point, weights, cluster_k, cluster_center_k, cluster_population_k, cluster_energy)
    !
    !write(*,'(A)' ) ' '
    !write(*,'(A,I0)' ) ' Number of k-means iterations taken with 2nd kmeans_w_03: ', it_num
    !call cluster_variance_compute (dim_num, npoints, cluster_num, point, cluster_k, cluster_center_k, cluster_variance)
    !call cluster_print_summary (npoints, cluster_num, cluster_population_k, cluster_energy, cluster_variance)
    !!@@@
    
    cluster_surrogate_k = 0
    minSqdistCluster = huge(1.0_wp)
    clusterMassConc = 0.0_wp
    do i = n1,nc
        if (MassConc_k(i) > 0.0_wp) then                    !only compare components of non-zero mass concentration;
            j = cluster_k(i)
            sqdist_i = sum((cluster_center_k(1:dim_num,j) - point(1:dim_num,i))**2)
            clusterMassConc(j) = clusterMassConc(j) + MassConc_k(i)
            if (sqdist_i < minSqdistCluster(j)) then        !save the new cluster_k surrogate component and the minimal value for this cluster_k
                cluster_surrogate_k(j) = i
                minSqdistCluster(j) = sqdist_i
            endif
        endif
    enddo
    !now that the surrogate components have been determined, save the cluster_k mass conc. under their component index and set all others to zero.
    MassConcLumped_k = 0.0_wp
    do j = 1,cluster_num
        MassConcLumped_k(cluster_surrogate_k(j)) = clusterMassConc(j)
    enddo
    
    !map the cluster_k center coordinates back to regular, unscaled x, y coordinates:
    do j = 1,cluster_num
        cluster_center_k(1,j) = cluster_center_k(1,j)*magnX / aspectRatio_XtoY
        cluster_center_k(2,j) = cluster_center_k(2,j)*magnY
    enddo

    deallocate (cluster_energy, cluster_variance, point, weights)

    end subroutine KmeansClusteringSurrogates
    !==================================================================================================================================
    
end module ModGridCellSurrogate