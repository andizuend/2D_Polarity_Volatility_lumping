subroutine ch_cap ( c )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
use Mod_kind_param, only : wp
implicit none

character c
integer itemp

itemp = ichar ( c )

if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
end if

return
end
function ch_eqi ( c1, c2 )

!*****************************************************************************80
!
!! CH_EQI is a case insensitive comparison of two characters for equality.
!
!  Example:
!
!    CH_EQI ( 'A', 'a' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
use Mod_kind_param, only : wp
implicit none

logical ch_eqi
character c1
character c1_cap
character c2
character c2_cap

c1_cap = c1
c2_cap = c2

call ch_cap ( c1_cap )
call ch_cap ( c2_cap )

if ( c1_cap == c2_cap ) then
    ch_eqi = .true.
else
    ch_eqi = .false.
end if

return
end


subroutine ch_to_digit ( c, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
!
!  Example:
!
!     C   DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer DIGIT, the corresponding value.  If C was
!    'illegal', then DIGIT is -1.
!
use Mod_kind_param, only : wp
implicit none

character c
integer digit

if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

else if ( c == ' ' ) then

    digit = 0

else

    digit = -1

end if

return
end


subroutine cluster_energy_compute ( dim_num, point_num, cluster_num, point, &
    cluster, cluster_center, cluster_energy )

!*****************************************************************************80
!
!! CLUSTER_ENERGY_COMPUTE computes the energy of the clusters.
!
!  Discussion:
!
!    The cluster energy is defined as the sum of the distance
!    squared from each point to its cluster center.  It is the goal
!    of the H-means and K-means algorithms to find, for a fixed number
!    of clusters, a clustering that minimizes this energy
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the number of spatial dimensions.
!
!    Input, integer POINT_NUM, the number of data points.
!
!    Input, integer CLUSTER_NUM, the number of clusters.
!
!    Input, real(wp) POINT(DIM_NUM,POINT_NUM), the data points.
!
!    Input, integer CLUSTER(POINT_NUM), the cluster to which each
!    data point belongs.
!
!    Input, real(wp) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), the
!    centers associated with the minimal energy clustering.
!
!    Output, real(wp) CLUSTER_ENERGY(CLUSTER_NUM), the energy
!    associated with each cluster.
!
use Mod_kind_param, only : wp
implicit none

integer cluster_num
integer dim_num
integer point_num

integer, dimension ( point_num ) :: cluster
real(wp), dimension ( dim_num, cluster_num ) :: cluster_center
real(wp), dimension ( cluster_num ) :: cluster_energy
integer i
integer j
real(wp), dimension ( dim_num, point_num ) :: point
real(wp) point_energy

cluster_energy(1:cluster_num) = 0.0_wp

do i = 1, point_num

    j = cluster(i)

    point_energy = sum ( &
        ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )

    cluster_energy(j) = cluster_energy(j) + point_energy

end do

return
end


subroutine cluster_initialize_1 ( dim_num, point_num, cluster_num, point, &
    cluster_center )

!*****************************************************************************80
!
!! CLUSTER_INITIALIZE_1 initializes the clusters to data points.
!
!  Discussion:
!
!    The cluster centers are simply chosen to be the first data points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the number of spatial dimensions.
!
!    Input, integer POINT_NUM, the number of points.
!
!    Input, integer CLUSTER_NUM, the number of clusters.
!
!    Input, real(wp) POINT(DIM_NUM,POINT_NUM), the coordinates
!    of the points.
!
!    Output, real(wp) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
!    the coordinates of the cluster centers.
!
use Mod_kind_param, only : wp
implicit none

integer cluster_num
integer dim_num
integer point_num

real(wp) cluster_center(dim_num,cluster_num)
integer i
real(wp) point(dim_num,point_num)

do i = 1, cluster_num
    cluster_center(1:dim_num,i) = point(1:dim_num,i)
end do

return
end


subroutine cluster_initialize_2 ( dim_num, point_num, cluster_num, point, &
    seed, cluster_center )

!*****************************************************************************80
!
!! CLUSTER_INITIALIZE_2 initializes the cluster centers to random values.
!
!  Discussion:
!
!    In this case, the hyperbox containing the data is computed.
!
!    Then the cluster centers are chosen uniformly at random within
!    this hyperbox.
!
!    Of course, if the data is not smoothly distributed throughout
!    the box, many cluster centers will be isolated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the number of spatial dimensions.
!
!    Input, integer POINT_NUM, the number of points.
!
!    Input, integer CLUSTER_NUM, the number of clusters.
!
!    Input, real(wp) POINT(DIM_NUM,POINT_NUM), the coordinates
!    of the points.
!
!    Input/output, integer SEED, a seed for the random
!    number generator.
!
!    Output, real(wp) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
!    the coordinates of the cluster centers.
!
use Mod_kind_param, only : wp
implicit none

integer cluster_num
integer dim_num
integer point_num

real(wp) cluster_center(dim_num,cluster_num)
integer i
real(wp) point(dim_num,point_num)
real(wp) r(dim_num)
real(wp) r_max(dim_num)
real(wp) r_min(dim_num)
integer seed

r_min = minval ( point, 2 )
r_max = maxval ( point, 2 )

do i = 1, cluster_num

    call r8vec_uniform_01 ( dim_num, seed, r )

    cluster_center(1:dim_num,i) = &
        ( 1.0_wp - r(1:dim_num) ) * r_min(1:dim_num) &
        + r(1:dim_num)   * r_max(1:dim_num)
end do

return
end


subroutine cluster_initialize_3 ( dim_num, point_num, cluster_num, point, &
    seed, cluster_center )

!*****************************************************************************80
!
!! CLUSTER_INITIALIZE_3 initializes the cluster centers to random values.
!
!  Discussion:
!
!    In this case, each point is randomly assigned to a cluster, and
!    the cluster centers are then computed as the centroids of the points
!    in the cluster.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the number of spatial dimensions.
!
!    Input, integer POINT_NUM, the number of points.
!
!    Input, integer CLUSTER_NUM, the number of clusters.
!
!    Input, real(wp) POINT(DIM_NUM,POINT_NUM), the coordinates
!    of the points.
!
!    Input/output, integer SEED, a seed for the random
!    number generator.
!
!    Output, real(wp) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
!    the coordinates of the cluster centers.
!
use Mod_kind_param, only : wp
implicit none

integer cluster_num
integer dim_num
integer point_num

real(wp) cluster_center(dim_num,cluster_num)
integer cluster_population(cluster_num)
integer i
integer i4_uniform
integer j
real(wp) point(dim_num,point_num)
integer seed
!
!  Assign one point to each cluster center.
!
do i = 1, cluster_num
    cluster_center(1:dim_num,i) = point(1:dim_num,i)
end do

cluster_population(1:cluster_num) = 1
!
!  The rest of the points get assigned randomly.
!
do i = cluster_num+1, point_num
    j = i4_uniform ( 1, cluster_num, seed )
    cluster_center(1:dim_num,j) = cluster_center(1:dim_num,j) + &
        point(1:dim_num,i)
    cluster_population(j) = cluster_population(j) + 1
end do
!
!  Now average the points to get the centroid.
!
do i = 1, cluster_num
    cluster_center(1:dim_num,i) = cluster_center(1:dim_num,i) / &
        real ( cluster_population(i), kind=wp )
end do

return
end


subroutine cluster_initialize_4 ( dim_num, point_num, cluster_num, point, &
    seed, cluster_center )

!*****************************************************************************80
!
!! CLUSTER_INITIALIZE_4 initializes the cluster centers to random values.
!
!  Discussion:
!
!    In this case, each data point is divided randomly among the
!    the cluster centers.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the number of spatial dimensions.
!
!    Input, integer POINT_NUM, the number of points.
!
!    Input, integer CLUSTER_NUM, the number of clusters.
!
!    Input, real(wp) POINT(DIM_NUM,POINT_NUM), the coordinates
!    of the points.
!
!    Input/output, integer SEED, a seed for the random
!    number generator.
!
!    Output, real(wp) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
!    the coordinates of the cluster centers.
!
use Mod_kind_param, only : wp
implicit none

integer cluster_num
integer dim_num
integer point_num

real(wp) cluster_center(dim_num,cluster_num)
real(wp) cluster_factor(cluster_num)
real(wp) cluster_weight(cluster_num)
real(wp) divisor
integer i
integer j
real(wp) point(dim_num,point_num)
integer seed

cluster_center(1:dim_num,1:cluster_num) = 0.0_wp
cluster_weight(1:cluster_num) = 0.0_wp

do i = 1, point_num

    call r8vec_uniform_01 ( cluster_num, seed, cluster_factor )

    divisor = sum ( cluster_factor(1:cluster_num) )
    cluster_factor(1:cluster_num) = cluster_factor(1:cluster_num) / divisor

    do j = 1, cluster_num
        cluster_center(1:dim_num,j) = cluster_center(1:dim_num,j) &
            + cluster_factor(j) * point(1:dim_num,i)
    end do

    cluster_weight(1:cluster_num) = cluster_weight(1:cluster_num) &
        + cluster_factor(1:cluster_num)

end do
!
!  Now normalize,  so that each cluster center is now a convex
!  combination of the points.
!
do i = 1, cluster_num
    cluster_center(1:dim_num,i) = cluster_center(1:dim_num,i) &
        / cluster_weight(i)
end do

return
end


subroutine cluster_initialize_5 ( dim_num, point_num, cluster_num, point, &
    seed, cluster_center )

!*****************************************************************************80
!
!! CLUSTER_INITIALIZE_5 initializes the cluster centers to random values.
!
!  Discussion:
!
!    In this case, each cluster center is a random convex combination
!    of the data points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the number of spatial dimensions.
!
!    Input, integer POINT_NUM, the number of points.
!
!    Input, integer CLUSTER_NUM, the number of clusters.
!
!    Input, real(wp) POINT(DIM_NUM,POINT_NUM), the coordinates
!    of the points.
!
!    Input/output, integer SEED, a seed for the random
!    number generator.
!
!    Output, real(wp) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
!    the coordinates of the cluster centers.
!
use Mod_kind_param, only : wp
implicit none

integer cluster_num
integer dim_num
integer point_num

real(wp) cluster_center(dim_num,cluster_num)
real(wp) column_sum
real(wp) factor(point_num,cluster_num)
integer j
real(wp) point(dim_num,point_num)
integer seed
!
!  Get a PxC block of random factors.
!
call r8mat_uniform_01 ( point_num, cluster_num, seed, factor )
!
!  Make each column of factors have unit sum.
!
do j = 1, cluster_num
    column_sum = sum ( factor(1:point_num,j) )
    factor(1:point_num,j) = factor(1:point_num,j) / column_sum
end do
!
!  Set centers = points * factors.
!
cluster_center(1:dim_num,1:cluster_num) = &
    matmul ( point(1:dim_num,1:point_num), factor(1:point_num,1:cluster_num) )

end subroutine cluster_initialize_5


subroutine cluster_print_summary ( point_num, cluster_num, &
    cluster_population, cluster_energy, cluster_variance )

!*****************************************************************************80
!
!! CLUSTER_PRINT_SUMMARY prints a summary of data about a clustering.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer POINT_NUM, the number of points.
!
!    Input, integer CLUSTER_NUM, the number of clusters.
!
!    Input, integer CLUSTER_POPULATION(CLUSTER_NUM), the number of
!    points assigned to each cluster.
!
!    Input, real(wp) CLUSTER_ENERGY(CLUSTER_NUM), the energy of
!    the clusters.
!
!    Input, real(wp) CLUSTER_VARIANCE(CLUSTER_NUM), the variance of
!    the clusters.
!
use Mod_kind_param, only : wp
implicit none

integer cluster_num

real(wp) ce
integer cep
real(wp) ce_total
real(wp) cluster_energy(cluster_num)
integer cluster_population(cluster_num)
real(wp) cluster_variance(cluster_num)
integer cp
integer cpp
real(wp) cv
integer i
integer point_num

ce_total = sum ( cluster_energy(1:cluster_num) )

write ( *, '(a)' ) ' '
write ( *, '(a)' ) '  Clustering statistics:'
write ( *, '(a)' ) ' '
write ( *, '(a,i8)' ) '    Number of clusters is ', cluster_num
write ( *, '(a,i8)' ) '    Number of points is   ', point_num
write ( *, '(a,g14.6)' ) '    Total energy is       ', ce_total
write ( *, '(a)' ) ' '
write ( *, '(a)' ) '    Cluster   Population        Energy          Variance'
write ( *, '(a)' ) &
    '    -------  -----------  -----------------  --------------'
write ( *, '(a)' ) '                  #    %     value        %'
write ( *, '(a)' ) ' '

do i = 1, cluster_num
    cp = cluster_population(i)
    cpp = int ( real ( 100 * cp, kind=wp ) / real ( point_num, kind=wp ) )
    ce = cluster_energy(i)
    cep = int ( ( ce * 100.0_wp ) / ce_total )
    cv = cluster_variance(i)
    write ( *, '(6x,i3,2x,i8,2x,i3,g14.6,2x,i3,2x,g14.6)' ) &
        i, cp, cpp, ce, cep, cv
end do

cp = sum ( cluster_population(1:cluster_num) )
cpp = 100
ce = sum ( cluster_energy(1:cluster_num) )
cep = 100
cv = sum ( cluster_population(1:cluster_num) &
    * cluster_variance(1:cluster_num) ) / cp

write ( *, '(a)' ) ' '
write ( *, '(a9,2x,i8,2x,i3,g14.6,2x,i3,2x,g14.6)' ) &
    '    Total', cp, cpp, ce, cep, cv

return
end


subroutine cluster_variance_compute ( dim_num, point_num, cluster_num, point, &
    cluster, cluster_center, cluster_variance )

!*****************************************************************************80
!
!! CLUSTER_VARIANCE_COMPUTE computes the variance of the clusters.
!
!
!  Discussion:
!
!    The cluster variance (from the cluster center) is the average of the
!    sum of the squares of the distances of each point in the cluster to the
!    cluster center.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the number of spatial dimensions.
!
!    Input, integer POINT_NUM, the number of data points.
!
!    Input, integer CLUSTER_NUM, the number of clusters.
!
!    Input, real(wp) POINT(DIM_NUM,POINT_NUM), the data points.
!
!    Input, integer CLUSTER(POINT_NUM), the cluster to which each
!    data point belongs.
!
!    Input, real(wp) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), the
!    centers associated with the minimal energy clustering.
!
!    Output, real(wp) CLUSTER_VARIANCE(CLUSTER_NUM), the variance
!    associated with each cluster.
!
use Mod_kind_param, only : wp
implicit none

integer cluster_num
integer dim_num
integer point_num

integer cluster(point_num)
real(wp) cluster_center(dim_num,cluster_num)
integer cluster_population(cluster_num)
real(wp) cluster_variance(cluster_num)
integer i
integer j
real(wp) point(dim_num,point_num)
real(wp) point_variance

cluster_population(1:cluster_num) = 0
cluster_variance(1:cluster_num) = 0.0_wp

do i = 1, point_num

    j = cluster(i)

    point_variance = sum ( &
        ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )

    cluster_variance(j) = cluster_variance(j) + point_variance
    cluster_population(j) = cluster_population(j) + 1

end do

cluster_variance(1:cluster_num) = cluster_variance(1:cluster_num) &
    / cluster_population(1:cluster_num)

return
end


subroutine file_column_count ( input_filename, column_num )

!*****************************************************************************80
!
!! FILE_COLUMN_COUNT counts the number of columns in the first line of a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!    Most lines of the file is presumed to consist of COLUMN_NUM words,
!    separated by spaces.  There may also be some blank lines, and some
!    comment lines,
!    which have a "#" in column 1.
!
!    The routine tries to find the first non-comment non-blank line and
!    counts the number of words in that line.
!
!    If all lines are blanks or comments, it goes back and tries to analyze
!    a comment line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FIlenAME, the name of the file.
!
!    Output, integer COLUMN_NUM, the number of columns in the file.
!
use Mod_kind_param, only : wp
implicit none

integer column_num
logical got_one
character ( len = * ) input_filename
integer input_status
integer input_unit
character ( len = 255 ) line
!
!  Open the file.
!
call get_unit ( input_unit )

open ( unit = input_unit, file = input_filename, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = input_status )

if ( input_status /= 0 ) then
    column_num = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' &
        // trim ( input_filename ) // '" on unit ', input_unit
    return
end if
!
!  Read one line, but skip blank lines and comment lines.
!
got_one = .false.

do

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
        exit
    end if

    if ( len_trim ( line ) == 0 ) then
        cycle
    end if

    if ( line(1:1) == '#' ) then
        cycle
    end if

    got_one = .true.
    exit

end do

if ( .not. got_one ) then

    rewind ( input_unit )

    do

        read ( input_unit, '(a)', iostat = input_status ) line

        if ( input_status /= 0 ) then
            exit
        end if

        if ( len_trim ( line ) == 0 ) then
            cycle
        end if

        got_one = .true.
        exit

    end do

end if

close ( unit = input_unit )

if ( .not. got_one ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Warning!'
    write ( *, '(a)' ) '  The file does not seem to contain any data.'
    column_num = -1
    return
end if

call s_word_count ( line, column_num )

return
end


subroutine file_row_count ( input_filename, row_num )

!*****************************************************************************80
!
!! FILE_ROW_COUNT counts the number of row records in a file.
!
!  Discussion:
!
!    It does not count lines that are blank, or that begin with a
!    comment symbol '#'.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FIlenAME, the name of the input file.
!
!    Output, integer ROW_NUM, the number of rows found.
!
use Mod_kind_param, only : wp
implicit none

integer bad_num
integer comment_num
integer ierror
character ( len = * ) input_filename
integer input_status
integer input_unit
character ( len = 255 ) line
integer record_num
integer row_num

call get_unit ( input_unit )

open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = input_status )

if ( input_status /= 0 ) then
    row_num = -1;
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_ROW_COUNT - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
        trim ( input_filename ) // '" on unit ', input_unit
    stop 1
end if

comment_num = 0
row_num = 0
record_num = 0
bad_num = 0

do

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
        ierror = record_num
        exit
    end if

    record_num = record_num + 1

    if ( line(1:1) == '#' ) then
        comment_num = comment_num + 1
        cycle
    end if

    if ( len_trim ( line ) == 0 ) then
        comment_num = comment_num + 1
        cycle
    end if

    row_num = row_num + 1

end do

close ( unit = input_unit )

return
end


subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the open command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IUNIT, the free unit number.
!
use Mod_kind_param, only : wp
implicit none

integer i
integer ios
integer iunit
logical lopen

iunit = 0

do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

        inquire ( unit = i, opened = lopen, iostat = ios )

        if ( ios == 0 ) then
            if ( .not. lopen ) then
                iunit = i
                return
            end if
        end if

    end if

end do

return
end


!subroutine hmeans_w_01 ( dim_num, point_num, cluster_num, it_max, it_num, &
!  point, weight, cluster, cluster_center, cluster_population, cluster_energy )
!use Mod_kind_param, only : wp
!
!!*****************************************************************************80
!!
!!! HMEANS_W_01 applies the weighted H-Means algorithm.
!!
!!  Discussion:
!!
!!    The input data for the weight H-Means problem includes:
!!    * a set of N data points X in M dimensions,
!!    * a set of N nonnegative weights W,
!!    * a desired number of clusters K.
!!    * an initial set of cluster centers Z,
!!    * an (optional) initial set of cluster assignments.
!!
!!    The goal is to determine K points Z, called cluster centers, and
!!    to assign each point X(I) to some cluster Z(J), so that we minimize
!!    the weighted standard deviation of the distance of each data point
!!    to the center of its cluster.  Writing J = CLUSTER(I) to
!!    indicate the index of the nearest cluster center Z(J) to the
!!    point X(I), the quantity we are trying to minimize is the sum
!!    of the weighted cluster energies E(J), where:
!!
!!      E(J) = Sum ( 1 <= I <= N ) W(I) * || X(I) - Z(J) ||**2
!!
!!    Here, we assume that we are using the Euclidean norm, so that
!!
!!      || X(I) - Z(J) ||**2 = Sum ( 1 <= K <= M )
!!        ( X(I)(K) - Z(J)(K) )**2
!!
!!    In this notation, X(I)(K) is the K-th spatial component of the
!!    I-th data point.
!!
!!    Note that this routine should give the same results as HMEANS_01
!!    in any case in which all the entries of the WEIGHT vector are equal.
!!
!!  Licensing:
!!
!!    This code is distributed under the GNU LGPL license.
!!
!!  Modified:
!!
!!    29 June 2006
!!
!!  Author:
!!
!!    John Burkardt
!!
!!  Reference:
!!
!!    Wendy Martinez, Angel Martinez,
!!    Computational Statistics Handbook with MATLAB,
!!    pages 373-376,
!!    Chapman and Hall / CRC, 2002.
!!
!!  Parameters:
!!
!!    Input, integer DIM_NUM, the number of spatial dimensions.
!!
!!    Input, integer POINT_NUM, the number of data points.
!!
!!    Input, integer CLUSTER_NUM, the number of clusters.
!!
!!    Input, integer IT_MAX, the maximum number of iterations.
!!
!!    Output, integer IT_NUM, the number of iterations taken.
!!
!!    Input, real(wp) POINT(DIM_NUM,POINT_NUM), the data points.
!!
!!    Input, real(wp) WEIGHT(POINT_NUM), the weights
!!    assigned to the data points.  These must be nonnegative, and
!!    at least one must be strictly positive.
!!
!!    Input/output, integer CLUSTER(POINT_NUM).  On input, the user
!!    may specify an initial cluster for each point, or leave all entrie of
!!    CLUSTER set to 0.  On output, CLUSTER contains the index of the
!!    cluster to which each data point belongs.
!!
!!    Input/output, real(wp) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), the
!!    centers associated with the minimal energy clustering.
!!
!!    Output, integer CLUSTER_POPULATION(CLUSTER_NUM),
!!    the populuation of each cluster.
!!
!!    Output, real(wp) CLUSTER_ENERGY(CLUSTER_NUM), the energy
!!    associated with each cluster.
!!
!use Mod_kind_param, only : wp
!  implicit none
!
!  integer cluster_num
!  integer dim_num
!  integer point_num
!
!  integer c
!  real(wp), dimension ( dim_num, cluster_num ) :: centroid
!  integer, dimension ( point_num ) :: cluster
!  real(wp), dimension ( dim_num, cluster_num ) :: cluster_center
!  real(wp), dimension ( cluster_num ) :: cluster_energy
!  integer, dimension ( cluster_num ) :: cluster_population
!  real(wp), dimension ( cluster_num ) :: cluster_weight
!  logical, parameter :: debug = .true.
!  real(wp) energy
!  integer i
!  integer it_max
!  integer it_num
!  integer j
!  integer missed
!  real(wp), dimension ( dim_num, point_num ) :: point
!  real(wp) point_energy
!  real(wp) point_energy_min
!  integer swap
!  real(wp) weight(point_num)
!!
!!  Data checks.
!!
!  if ( cluster_num < 1 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'HMEANS_W_01 - Fatal error!'
!    write ( *, '(a)' ) '  CLUSTER_NUM < 1.'
!    stop 1
!  end if
!
!  if ( dim_num < 1 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'HMEANS_W_01 - Fatal error!'
!    write ( *, '(a)' ) '  DIM_NUM < 1.'
!    stop 1
!  end if
!
!  if ( point_num < 1 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'HMEANS_W_01 - Fatal error!'
!    write ( *, '(a)' ) '  POINT_NUM < 1.'
!    stop 1
!  end if
!
!  if ( it_max < 0 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'HMEANS_W_01 - Fatal error!'
!    write ( *, '(a)' ) '  IT_MAX < 0.'
!    stop 1
!  end if
!
!  if ( any ( weight(1:point_num) < 0.0_wp ) ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'HMEANS_W_01 - Fatal error!'
!    write ( *, '(a)' ) '  Some weight entry is negative.'
!    stop 1
!  end if
!
!  if ( all ( weight(1:point_num) <= 0.0_wp ) ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'HMEANS_W_01 - Fatal error!'
!    write ( *, '(a)' ) '  No weight entry is positive.'
!    stop 1
!  end if
!!
!!  On input, legal entries in CLUSTER are preserved, but
!!  otherwise, each point is assigned to its nearest cluster.
!!
!  do i = 1, point_num
!    if ( cluster(i) <= 0 .or. cluster_num < cluster(i) ) then
!
!      point_energy_min = huge ( point_energy_min )
!
!      do j = 1, cluster_num
!
!        point_energy = sum ( &
!          ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )
!
!        if ( point_energy < point_energy_min ) then
!          point_energy_min = point_energy
!          cluster(i) = j
!        end if
!
!      end do
!
!    end if
!  end do
!
!  it_num = 0
!
!  do while ( it_num < it_max )
!
!    it_num = it_num + 1
!!
!!  #1:
!!  Reassign points to clusters:
!!  Assign each point to the cluster whose center is nearest;
!!  Count the number of points whose cluster assignment is changed.
!!
!    swap = 0
!
!    do i = 1, point_num
!
!      point_energy_min = huge ( point_energy_min )
!      c = cluster(i)
!
!      do j = 1, cluster_num
!
!        point_energy = sum ( &
!          ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )
!
!        if ( point_energy < point_energy_min ) then
!          point_energy_min = point_energy
!          cluster(i) = j
!        end if
!
!      end do
!
!      if ( c /= cluster(i) ) then
!        swap = swap + 1
!      end if
!
!    end do
!!
!!  If no point changed its cluster assignment, the algorithm can make no
!!  more improvements, so terminate.
!!
!    if ( 1 < it_num ) then
!      if ( swap == 0 ) then
!        exit
!      end if
!    end if
!!
!!  Determine the current energy.
!!
!    energy = 0.0
!    do i = 1, point_num
!      energy = energy + weight(i) * sum ( &
!          ( point(1:dim_num,i) - cluster_center(1:dim_num,cluster(i)) )**2 )
!    end do
!    write ( *, * ) it_num, energy
!!
!!  #2:
!!  Determine the centroids of the clusters, and set the
!!  cluster center to the cluster centroid.
!!
!    centroid(1:dim_num,1:cluster_num) = 0.0_wp
!    cluster_population(1:cluster_num) = 0
!    cluster_weight(1:cluster_num) = 0.0_wp
!
!    do i = 1, point_num
!      c = cluster(i)
!      cluster_population(c) = cluster_population(c) + 1
!      cluster_weight(c) = cluster_weight(c) + weight(i)
!      centroid(1:dim_num,c) = centroid(1:dim_num,c) &
!        + weight(i) * point(1:dim_num,i)
!    end do
!
!    missed = 0
!
!    do c = 1, cluster_num
!
!      if ( cluster_weight(c) /= 0.0_wp ) then
!        centroid(1:dim_num,c) = centroid(1:dim_num,c) / cluster_weight(c)
!      else
!        missed = missed + 1
!        centroid(1:dim_num,c) = point(1:dim_num,missed)
!      end if
!
!    end do
!
!    cluster_center(1:dim_num,1:cluster_num) = centroid(1:dim_num,1:cluster_num)
!
!  end do
!!
!!  Compute the energy based on the final value of the cluster centers.
!!
!  cluster_energy(1:cluster_num) = 0.0_wp
!
!  do i = 1, point_num
!
!    c = cluster(i)
!
!    cluster_energy(c) = cluster_energy(c) + weight(i) * sum ( &
!      ( point(1:dim_num,i) - cluster_center(1:dim_num,c) )**2 )
!
!  end do
!
!  return
!end


!subroutine hmeans_w_02 ( dim_num, point_num, cluster_num, it_max, it_num, &
!  point, weight, cluster, cluster_center, cluster_population, &
!  cluster_energy, seed )
!
!!*****************************************************************************80
!!
!!! HMEANS_W_02 applies the weighted H-Means algorithm.
!!
!!  Discussion:
!!
!!    The input data for the weight H-Means problem includes:
!!    * a set of N data points X in M dimensions,
!!    * a set of N nonnegative weights W,
!!    * a desired number of clusters K.
!!    * an initial set of cluster centers Z,
!!    * an (optional) initial set of cluster assignments.
!!
!!    The goal is to determine K points Z, called cluster centers, and
!!    to assign each point X(I) to some cluster Z(J), so that we minimize
!!    the weighted standard deviation of the distance of each data point
!!    to the center of its cluster.  Writing J = CLUSTER(I) to
!!    indicate the index of the nearest cluster center Z(J) to the
!!    point X(I), the quantity we are trying to minimize is the sum
!!    of the weighted cluster energies E(J), where:
!!
!!      E(J) = Sum ( 1 <= I <= N ) W(I) * || X(I) - Z(J) ||**2
!!
!!    Here, we assume that we are using the Euclidean norm, so that
!!
!!      || X(I) - Z(J) ||**2 = Sum ( 1 <= K <= M )
!!        ( X(I)(K) - Z(J)(K) )**2
!!
!!    In this notation, X(I)(K) is the K-th spatial component of the
!!    I-th data point.
!!
!!    Note that this routine should give the same results as HMEANS_02
!!    in any case in which all the entries of the WEIGHT vector are equal.
!!
!!  Licensing:
!!
!!    This code is distributed under the GNU LGPL license.
!!
!!  Modified:
!!
!!    29 June 2006
!!
!!  Author:
!!
!!    John Burkardt
!!
!!  Parameters:
!!
!!    Input, integer DIM_NUM, the number of spatial dimensions.
!!
!!    Input, integer POINT_NUM, the number of points.
!!
!!    Input, integer CLUSTER_NUM, the number of clusters.
!!
!!    Input, integer IT_MAX, the maximum number of iterations.
!!
!!    Output, integer IT_NUM, the number of iterations taken.
!!
!!    Input, real(wp) POINT(DIM_NUM,POINT_NUM), the coordinates
!!    of the points.
!!
!!    Input/output, integer CLUSTER(POINT_NUM).  On input, the user
!!    may specify an initial cluster for each point, or leave all entrie of
!!    CLUSTER set to 0.  On output, CLUSTER contains the index of the
!!    cluster to which each data point belongs.
!!
!!    Input/output, real(wp) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
!!    the coordinates of the cluster centers.
!!
!!    Output, integer CLUSTER_POPULATION(CLUSTER_NUM), the number of
!!    points assigned to each cluster.
!!
!!    Output, real(wp) CLUSTER_ENERGY(CLUSTER_NUM), the energy of
!!    the clusters.
!!
!!    Input/output, integer SEED, a seed for the random
!!    number generator.
!!
!use Mod_kind_param, only : wp
!  implicit none
!
!  integer cluster_num
!  integer dim_num
!  integer point_num
!
!  integer c
!  integer cluster(point_num)
!  real(wp) cluster_center(dim_num,cluster_num)
!  real(wp) cluster_energy(cluster_num)
!  integer cluster_population(cluster_num)
!  real(wp) cluster_weight(cluster_num)
!  logical, parameter :: debug = .false.
!  real(wp) energy(cluster_num)
!  integer i
!  integer i4_uniform
!  integer it_max
!  integer it_num
!  integer j
!  integer list(1)
!  real(wp) point(dim_num,point_num)
!  real(wp) point_energy
!  real(wp) point_energy_min
!  integer seed
!  integer swap
!  real(wp) weight(point_num)
!!
!!  Data checks.
!!
!  if ( cluster_num < 1 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'HMEANS_W_02 - Fatal error!'
!    write ( *, '(a)' ) '  CLUSTER_NUM < 1.'
!    stop 1
!  end if
!
!  if ( dim_num < 1 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'HMEANS_W_02 - Fatal error!'
!    write ( *, '(a)' ) '  DIM_NUM < 1.'
!    stop 1
!  end if
!
!  if ( point_num < 1 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'HMEANS_W_02 - Fatal error!'
!    write ( *, '(a)' ) '  POINT_NUM < 1.'
!    stop 1
!  end if
!
!  if ( it_max < 0 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'HMEANS_W_02 - Fatal error!'
!    write ( *, '(a)' ) '  IT_MAX < 0.'
!    stop 1
!  end if
!
!  if ( any ( weight(1:point_num) < 0.0_wp ) ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'HMEANS_W_02 - Fatal error!'
!    write ( *, '(a)' ) '  Some weight entry is negative.'
!    stop 1
!  end if
!
!  if ( all ( weight(1:point_num) <= 0.0_wp ) ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'HMEANS_W_02 - Fatal error!'
!    write ( *, '(a)' ) '  No weight entry is positive.'
!    stop 1
!  end if
!!
!!  On input, legal entries in CLUSTER are preserved, but
!!  otherwise, each point is assigned to its nearest cluster.
!!
!  do i = 1, point_num
!    if ( cluster(i) <= 0 .or. cluster_num < cluster(i) ) then
!
!      point_energy_min = huge ( point_energy_min )
!
!      do j = 1, cluster_num
!
!        point_energy = sum ( &
!          ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )
!
!        if ( point_energy < point_energy_min ) then
!          point_energy_min = point_energy
!          cluster(i) = j
!        end if
!
!      end do
!
!    end if
!  end do
!
!  it_num = 0
!
!  do
!!
!!  Given centers, assign points to nearest center.
!!
!    cluster_population(1:cluster_num) = 0
!    cluster_weight(1:cluster_num) = 0.0_wp
!    cluster_energy(1:cluster_num) = 0.0_wp
!
!    swap = 0
!
!    do i = 1, point_num
!
!      do j = 1, cluster_num
!        energy(j) = sum ( &
!          ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )
!      end do
!
!      list = minloc ( energy(1:cluster_num) )
!
!      c = list(1)
!
!      if ( c /= cluster(i) ) then
!        swap = swap + 1
!      end if
!
!      cluster(i) = c
!      cluster_energy(c) = cluster_energy(c) + weight(i) * energy(c)
!      cluster_population(c) = cluster_population(c) + 1
!      cluster_weight(c) = cluster_weight(c) + weight(i)
!
!    end do
!
!    if ( debug ) then
!      write ( *, '(i3,g14.6)' ) it_num, sum ( cluster_energy(1:cluster_num) )
!    end if
!
!    if ( 0 < it_num ) then
!      if ( swap == 0 ) then
!        exit
!      end if
!    end if
!
!    if ( it_max <= it_num ) then
!      exit
!    end if
!
!    it_num = it_num + 1
!!
!!  Given points in cluster, replace center by weighted centroid.
!!
!    cluster_center(1:dim_num,1:cluster_num) = 0.0_wp
!
!    do i = 1, point_num
!      c = cluster(i)
!      cluster_center(1:dim_num,c) = cluster_center(1:dim_num,c) &
!        + weight(i) * point(1:dim_num,i)
!    end do
!
!    do i = 1, cluster_num
!      if ( cluster_weight(i) /= 0 ) then
!        cluster_center(1:dim_num,i) = cluster_center(1:dim_num,i) / &
!          cluster_weight(i)
!      else
!        j = i4_uniform ( 1, point_num, seed )
!        cluster_center(1:dim_num,i) = point(1:dim_num,j)
!      end if
!    end do
!
!  end do
!!
!!  Compute the energy based on the final value of the cluster centers.
!!
!  cluster_energy(1:cluster_num) = 0.0_wp
!
!  do i = 1, point_num
!
!    c = cluster(i)
!
!    cluster_energy(c) = cluster_energy(c) + weight(i) * sum ( &
!      ( point(1:dim_num,i) - cluster_center(1:dim_num,c) )**2 )
!
!  end do
!
!  return
!end
function i4_uniform ( a, b, seed )

!*****************************************************************************80
!
!! I4_UNifORM returns a scaled pseudorandom I4.
!
!  Discussion:
!
!    An I4 is an integer value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer A, B, the limits of the interval.
!
!    Input/output, integer SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer I4_UNifORM, a number between A and B.
!
use Mod_kind_param, only : wp
implicit none

integer a
integer b
integer i4_uniform
integer k
real ( kind = 4 ) r
integer seed
integer value

if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_UNifORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
end if

k = seed / 127773

seed = 16807 * ( seed - k * 127773 ) - k * 2836

if ( seed < 0 ) then
    seed = seed + 2147483647
end if

r = real ( seed, kind = 4 ) * 4.656612875E-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
r = ( 1.0E+00 - r ) * ( real ( min ( a, b ), kind = 4 ) - 0.5E+00 ) &
    +             r   * ( real ( max ( a, b ), kind = 4 ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
value = nint ( r, kind = 4 )

value = max ( value, min ( a, b ) )
value = min ( value, max ( a, b ) )

i4_uniform = value

return
end


subroutine i4mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! I4MAT_write writes an I4MAT file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FIlenAME, the output file name.
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Input, integer TABLE(M,N), the table data.
!
use Mod_kind_param, only : wp
implicit none

integer m
integer n

integer j
character ( len = * ) output_filename
integer output_status
integer output_unit
character ( len = 30 ) string
integer table(m,n)
!
!  Open the file.
!
call get_unit ( output_unit )

open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_write - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
        trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop 1
end if
!
!  Create a format string.
!
write ( string, '(a1,i8,a4)' ) '(', m, 'i10)'
!
!  Write the data.
!
do j = 1, n
    write ( output_unit, string ) table(1:m,j)
end do
!
!  Close the file.
!
close ( unit = output_unit )

return
end


subroutine kmeans_w_01 ( dim_num, point_num, cluster_num, it_max, it_num, &
    point, weight, cluster, cluster_center, cluster_population, cluster_energy )

!*****************************************************************************80
!
!! KMEANS_W_01 applies the weighted K-Means algorithm.
!
!  Discussion:
!
!    The input data for the weight K-Means problem includes:
!    * a set of N data points X in M dimensions,
!    * a set of N nonnegative weights W,
!    * a desired number of clusters K.
!    * an initial set of cluster centers Z,
!    * an (optional) initial set of cluster assignments.
!
!    The goal is to determine K points Z, called cluster centers, and
!    to assign each point X(I) to some cluster Z(J), so that we minimize
!    the weighted standard deviation of the distance of each data point
!    to the center of its cluster.  Writing J = CLUSTER(I) to
!    indicate the index of the nearest cluster center Z(J) to the
!    point X(I), the quantity we are trying to minimize is the sum
!    of the weighted cluster energies E(J), where:
!
!      E(J) = Sum ( 1 <= I <= N ) W(I) * || X(I) - Z(J) ||^2
!
!    Here, we assume that we are using the Euclidean norm, so that
!
!      || X(I) - Z(J) ||^2 = Sum ( 1 <= K <= M )
!         ( X(I)(K) - Z(J)(K) )^2
!
!    In this notation, X(I)(K) is the K-th spatial component of the
!    I-th data point.
!
!    Note that this routine should give the same results as KMEANS_01
!    in any case in which all the entries of the WEIGHT vector are equal.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!    10 October 2011
!
!  Author:
!    John Burkardt
!
!  Reference:
!    David Sparks,
!    Algorithm AS 58: Euclidean Cluster Analysis,
!    Applied Statistics,
!    Volume 22, Number 1, 1973, pages 126-130.
!
!  Parameters:
!
!    Input, integer DIM_NUM, the number of spatial dimensions.
!    Input, integer POINT_NUM, the number of points.
!    Input, integer CLUSTER_NUM, the number of clusters.
!    Input, integer IT_MAX, the maximum number of iterations.
!    Output, integer IT_NUM, the number of iterations taken.
!    Input, real(wp) POINT(DIM_NUM,POINT_NUM), the points.
!    Input, real(wp) WEIGHT(POINT_NUM), the weights.
!    Output, integer CLUSTER(POINT_NUM), indicates which cluster
!    each point belongs to.
!    Input/output, real(wp) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
!    the cluster centers.
!    Output, integer CLUSTER_POPULATION(CLUSTER_NUM), the number
!    of points in each cluster.
!    Output, real(wp) CLUSTER_ENERGY(CLUSTER_NUM), the
!    cluster energies.
!
use Mod_kind_param, only : wp
implicit none

integer,intent(in) :: cluster_num
integer,intent(in) :: dim_num
integer,intent(in) :: point_num
integer,intent(in) :: it_max
integer,intent(out) :: it_num
real(wp),intent(in) :: point(dim_num, point_num)
real(wp),intent(in) :: weight(point_num)
integer,intent(out) :: cluster(point_num)
real(wp),intent(inout) :: cluster_center(dim_num, cluster_num)
real(wp),intent(out) :: cluster_energy(cluster_num)
integer,intent(out) :: cluster_population(cluster_num)
!local variables:
integer :: c
real(wp) :: cluster_weight(cluster_num)
real(wp) :: dc
real(wp) :: de
real(wp) :: f(point_num)
integer :: i
integer :: il
integer :: ir
integer :: j
integer :: k
integer :: list(1)
integer :: swap

it_num = 0
!
!  Idiot checks.
!
if ( cluster_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_W_01 - Fatal error!'
    write ( *, '(a)' ) '  CLUSTER_NUM < 1.'
    stop 1
end if

if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_W_01 - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    stop 1
end if

if ( point_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_W_01 - Fatal error!'
    write ( *, '(a)' ) '  POINT_NUM < 1.'
    stop 1
end if

if ( it_max < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_W_01 - Fatal error!'
    write ( *, '(a)' ) '  IT_MAX < 0.'
    stop 1
end if

if ( any ( weight(1:point_num) < 0.0_wp ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_W_01 - Fatal error!'
    write ( *, '(a)' ) '  Some weight entry is negative.'
    stop 1
end if

if ( all ( weight(1:point_num) <= 0.0_wp ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_W_01 - Fatal error!'
    write ( *, '(a)' ) '  No weight entry is positive.'
    stop 1
end if
!
!  Assign each point to the nearest cluster.
!
do i = 1, point_num
    do j = 1, cluster_num
        cluster_energy(j) = sum( ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )
    end do
    list = minloc(cluster_energy(1:cluster_num))
    cluster(i) = list(1)
end do
!
!  Determine the cluster populations and weights.
!
cluster_population(1:cluster_num) = 0
cluster_weight(1:cluster_num) = 0.0_wp

do i = 1, point_num
    j = cluster(i)
    cluster_population(j) = cluster_population(j) + 1
    cluster_weight(j) = cluster_weight(j) + weight(i)
end do
!
!  Calculate the mean and sum of squares for each cluster.
!
cluster_center(1:dim_num,1:cluster_num) = 0.0_wp

do i = 1, point_num
    j = cluster(i)
    cluster_center(1:dim_num,j) = cluster_center(1:dim_num,j) &
        + weight(i) * point(1:dim_num,i)
end do

do i = 1, cluster_num
    if ( 0.0_wp < cluster_weight(i) ) then
        cluster_center(1:dim_num,i) = cluster_center(1:dim_num,i) &
            / cluster_weight(i)
    end if
end do
!
!  Set the point energies.
!
f(1:point_num) = 0.0_wp

do i = 1, point_num
    j = cluster(i)
    f(i) = sum ( ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )
end do
!
!  Set the cluster energies.
!
cluster_energy(1:cluster_num) = 0.0_wp

do i = 1, point_num
    j = cluster(i)
    cluster_energy(j) = cluster_energy(j) + weight(i) * f(i)
end do
!
!  Adjust the point energies by a weight factor.
!
do i = 1, point_num
    j = cluster(i)
    if ( weight(i) < cluster_weight(j) ) then
        f(i) = f(i) * cluster_weight(j) / ( cluster_weight(j) - weight(i) )
    end if
end do
!
!  Examine each observation in turn to see if it should be
!  reassigned to a different cluster.
!
it_num = 0

do while ( it_num < it_max )

    it_num = it_num + 1

    swap = 0

    do i = 1, point_num

        il = cluster(i)
        ir = il

        if ( cluster_population(il) <= 1 ) then
            cycle
        end if

        dc = f(i)

        do j = 1, cluster_num

            if ( j /= il ) then

                de = sum ( &
                    ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 ) &
                    * cluster_weight(j) / ( cluster_weight(j) + weight(i) )

                if ( de < dc ) then
                    dc = de
                    ir = j
                end if

            end if

        end do
        !
        !  If the lowest value was obtained by staying in the current cluster,
        !  then cycle.
        !
        if ( ir == il ) then
            cycle
        end if
        !
        !  Reassign the point from cluster IL to cluster IR.
        !
        cluster_center(1:dim_num,il) = &
            ( cluster_weight(il) * cluster_center(1:dim_num,il) &
            - weight(i) * point(1:dim_num,i) ) &
            / ( cluster_weight(il) - weight(i) )

        cluster_center(1:dim_num,ir) = &
            ( cluster_weight(ir) * cluster_center(1:dim_num,ir) &
            + weight(i) * point(1:dim_num,i) ) &
            / ( cluster_weight(ir) + weight(i) )

        cluster_weight(il) = cluster_weight(il) - weight(i)
        cluster_weight(ir) = cluster_weight(ir) + weight(i)

        cluster_energy(il) = cluster_energy(il) - f(i)
        cluster_energy(ir) = cluster_energy(ir) + dc

        cluster_population(ir) = cluster_population(ir) + 1
        cluster_population(il) = cluster_population(il) - 1

        cluster(i) = ir
        !
        !  Adjust the value of F for all points in clusters IL and IR.
        !
        do j = 1, point_num
            k = cluster(j)
            if ( k == il .or. k == ir ) then
                f(j) = sum ( &
                    ( point(1:dim_num,j) - cluster_center(1:dim_num,k) )**2 )
                if ( weight(j) < cluster_weight(k) ) then
                    f(j) = f(j) * cluster_weight(k) &
                        / ( cluster_weight(k) - weight(j) )
                end if
            end if
        end do
        swap = swap + 1
    end do
    !
    !  Exit if no reassignments were made during this iteration.
    !
    if ( swap == 0 ) then
        exit
    end if

end do
!
!  Compute the energy based on the final value of the cluster centers.
!
cluster_energy(1:cluster_num) = 0.0_wp

do i = 1, point_num
    c = cluster(i)
    cluster_energy(c) = cluster_energy(c) + weight(i) * sum ( &
        ( point(1:dim_num,i) - cluster_center(1:dim_num,c) )**2 )
end do

end subroutine kmeans_w_01


subroutine kmeans_w_03 ( dim_num, point_num, cluster_num, it_max, it_num, &
    point, weight, cluster, cluster_center, cluster_population, cluster_energy )

!*****************************************************************************80
!
!! KMEANS_W_03 applies the weighted K-Means algorithm.
!
!  Discussion:
!
!    The input data for the weight K-Means problem includes:
!    * a set of N data points X in M dimensions,
!    * a set of N nonnegative weights W,
!    * a desired number of clusters K.
!    * an initial set of cluster centers Z,
!    * an (optional) initial set of cluster assignments.
!
!    The goal is to determine K points Z, called cluster centers, and
!    to assign each point X(I) to some cluster Z(J), so that we minimize
!    the weighted standard deviation of the distance of each data point
!    to the center of its cluster.  Writing J = CLUSTER(I) to
!    indicate the index of the nearest cluster center Z(J) to the
!    point X(I), the quantity we are trying to minimize is the sum
!    of the weighted cluster energies E(J), where:
!
!      E(J) = Sum ( 1 <= I <= N ) W(I) * || X(I) - Z(J) ||**2
!
!    Here, we assume that we are using the Euclidean norm, so that
!
!      || X(I) - Z(J) ||**2 = Sum ( 1 <= K <= M )
!        ( X(I)(K) - Z(J)(K) )**2
!
!    In this notation, X(I)(K) is the K-th spatial component of the
!    I-th data point.
!
!    Note that this routine should give the same results as KMEANS_01
!    in any case in which all the entries of the WEIGHT vector are equal.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Wendy Martinez, Angel Martinez,
!    Computational Statistics Handbook with MATLAB,
!    pages 373-376,
!    Chapman and Hall / CRC, 2002.
!
!  Parameters:
!
!    Input, integer DIM_NUM, the number of spatial dimensions.
!
!    Input, integer POINT_NUM, the number of data points.
!
!    Input, integer CLUSTER_NUM, the number of clusters.
!
!    Input, integer IT_MAX, the maximum number of iterations.
!
!    Output, integer IT_NUM, the number of iterations taken.
!
!    Input, real(wp) POINT(DIM_NUM,POINT_NUM), the data points.
!
!    Input, real(wp) WEIGHT(POINT_NUM), the weights.
!
!    Input/output, integer CLUSTER(POINT_NUM), the cluster
!    to which each point belongs.  On output, these may have been altered.
!
!    Input/output, real(wp) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM), the
!    centers associated with the clustering.  On output, these may
!    have been altered.
!
!    Output, integer CLUSTER_POPULATION(CLUSTER_NUM), the number
!    of points in each cluster.
!
!    Output, real(wp) CLUSTER_ENERGY(CLUSTER_NUM), the energy of
!    the clusters.
!
use Mod_kind_param, only : wp
implicit none

integer cluster_num
integer dim_num
integer point_num

integer ci
integer cj
integer, dimension ( point_num ) :: cluster
real(wp), dimension ( dim_num, cluster_num ) :: cluster_center
real(wp), dimension ( cluster_num ) :: cluster_energy
integer, dimension ( cluster_num ) :: cluster_population
real(wp), dimension ( cluster_num ) :: cluster_weight
logical, parameter :: debug = .true.
real(wp), dimension ( cluster_num ) :: distsq
integer i
integer it_max
integer it_num
integer j
integer list(1)
real(wp), dimension ( dim_num, point_num ) :: point
integer swap
real(wp) weight(point_num)
!
!  Check the input.
!
if ( cluster_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_W_03 - Fatal error!'
    write ( *, '(a)' ) '  CLUSTER_NUM < 1.'
    stop 1
end if

if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_W_03 - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    stop 1
end if

if ( point_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_W_03 - Fatal error!'
    write ( *, '(a)' ) '  POINT_NUM < 1.'
    stop 1
end if

if ( it_max < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_W_03 - Fatal error!'
    write ( *, '(a)' ) '  IT_MAX < 0.'
    stop 1
end if

if ( any ( weight(1:point_num) < 0.0_wp ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_W_03 - Fatal error!'
    write ( *, '(a)' ) '  Some weight entry is negative.'
    stop 1
end if

if ( all ( weight(1:point_num) <= 0.0_wp ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_W_03 - Fatal error!'
    write ( *, '(a)' ) '  No weight entry is positive.'
    stop 1
end if
!
!  Assign each observation to the nearest cluster center.
!
do i = 1, point_num

    do j = 1, cluster_num
        cluster_energy(j) = sum ( &
            ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )
    end do

    list = minloc ( cluster_energy(1:cluster_num) )
    cluster(i) = list(1)

end do
!
!  Determine the cluster populations and weights.
!
cluster_population(1:cluster_num) = 0
cluster_weight(1:cluster_num) = 0.0_wp

do i = 1, point_num
    ci = cluster(i)
    cluster_population(ci) = cluster_population(ci) + 1
    cluster_weight(ci) = cluster_weight(ci) + weight(i)
end do
!
!  Average the points in each cluster to get a new cluster center.
!
cluster_center(1:dim_num,1:cluster_num) = 0.0_wp

do i = 1, point_num
    j = cluster(i)
    cluster_center(1:dim_num,j) = cluster_center(1:dim_num,j) &
        + weight(i) * point(1:dim_num,i)
end do

do i = 1, cluster_num
    if ( 0.0_wp < cluster_weight(i) ) then
        cluster_center(1:dim_num,i) = cluster_center(1:dim_num,i) &
            / cluster_weight(i)
    end if
end do
!
!  Carry out the iteration.
!
it_num = 0

do while ( it_num < it_max )

    it_num = it_num + 1

    swap = 0

    do i = 1, point_num

        ci = cluster(i)

        if ( cluster_population(ci) <= 1 ) then
            cycle
        end if

        do cj = 1, cluster_num

            if ( cj == ci ) then

                distsq(cj) = sum ( &
                    ( point(1:dim_num,i) - cluster_center(1:dim_num,cj) )**2 ) &
                    * cluster_weight(cj)  &
                    / ( cluster_weight(cj) - weight(i) )

            else if ( cluster_population(cj) == 0 ) then

                cluster_center(1:dim_num,cj) = point(1:dim_num,i)
                distsq(cj) = 0.0_wp

            else

                distsq(cj) = sum ( &
                    ( point(1:dim_num,i) - cluster_center(1:dim_num,cj) )**2 ) &
                    * cluster_weight(cj) &
                    / ( cluster_weight(cj) + weight(i) )

            end if

        end do
        !
        !  Find the index of the minimum value of DISTSQ.
        !
        list = minloc ( distsq(1:cluster_num) )
        !
        !  If that's not the cluster to which point I now belongs, move it there.
        !
        if ( list(1) == ci ) then
            cycle
        end if

        cj = list(1)

        cluster_center(1:dim_num,ci) = &
            ( cluster_weight(ci) * cluster_center(1:dim_num,ci) &
            - weight(i) * point(1:dim_num,i) ) &
            / ( cluster_weight(ci) - weight(i) )

        cluster_center(1:dim_num,cj) = &
            ( cluster_weight(cj) * cluster_center(1:dim_num,cj) &
            + weight(i) * point(1:dim_num,i) ) &
            / ( cluster_weight(cj) + weight(i) )

        cluster_population(ci) = cluster_population(ci) - 1
        cluster_population(cj) = cluster_population(cj) + 1

        cluster_weight(ci) = cluster_weight(ci) - weight(i)
        cluster_weight(cj) = cluster_weight(cj) + weight(i)

        cluster(i) = cj

        swap = swap + 1

    end do
    !
    !  Exit if no reassignments were made during this iteration.
    !
    if ( swap == 0 ) then
        exit
    end if

end do
!
!  Compute the energy based on the final value of the cluster centers.
!
cluster_energy(1:cluster_num) = 0.0_wp

do i = 1, point_num

    ci = cluster(i)

    cluster_energy(ci) = cluster_energy(ci) + weight(i) * sum ( &
        ( point(1:dim_num,i) - cluster_center(1:dim_num,ci) )**2 )

end do

end subroutine kmeans_w_03


function r4_uniform_01 ( seed )

!*****************************************************************************80
!
!! R4_UNifORM_01 returns a unit pseudorandom R4.
!
!  Discussion:
!
!    An R4 is a real ( kind = 4 ) value.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r4_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R4_UNifORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer SEED, the "seed" value,
!    which should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 4 ) R4_UNifORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
use Mod_kind_param, only : wp
implicit none

integer k
integer seed
real ( kind = 4 ) r4_uniform_01

if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_UNifORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
end if

k = seed / 127773

seed = 16807 * ( seed - k * 127773 ) - k * 2836

if ( seed < 0 ) then
    seed = seed + 2147483647
end if

r4_uniform_01 = real ( seed, kind = 4 ) * 4.656612875E-10

return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNifORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real(wp) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNifORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real(wp) R8_UNifORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
use Mod_kind_param, only : wp
implicit none

integer k
real(wp) r8_uniform_01
integer seed

if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNifORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
end if

k = seed / 127773

seed = 16807 * ( seed - k * 127773 ) - k * 2836

if ( seed < 0 ) then
    seed = seed + 2147483647
end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
r8_uniform_01 = real ( seed, kind=wp ) * 4.656612875E-10_wp

return
end


subroutine r8mat_data_read ( input_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_DATA_read reads data from an R8MAT file.
!
!  Discussion:
!
!    The file may contain more than N points, but this routine will
!    return after reading N of them.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FIlenAME, the name of the input file.
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Output, real(wp) TABLE(M,N), the table data.
!
use Mod_kind_param, only : wp
implicit none

integer m
integer n

integer ierror
character ( len = * ) input_filename
integer input_status
integer input_unit
integer j
character ( len = 255 ) line
real(wp) table(m,n)
real(wp) x(m)

ierror = 0

call get_unit ( input_unit )

open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = input_status )

if ( input_status /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_DATA_read - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
        trim ( input_filename ) // '" on unit ', input_unit
    stop 1
end if

j = 0

do while ( j < n )

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8MAT_DATA_read - Fatal error!'
        write ( *, '(a)' ) '  Error while reading lines of data.'
        write ( *, '(a,i8)' ) '  Number of values expected per line M = ', m
        write ( *, '(a,i8)' ) '  Number of data lines read, J =         ', j
        write ( *, '(a,i8)' ) '  Number of data lines needed, N =       ', n
        stop 1
    end if

    if ( line(1:1) == '#' .or. len_trim ( line ) == 0 ) then
        cycle
    end if

    call s_to_r8vec ( line, m, x, ierror )

    if ( ierror /= 0 ) then
        cycle
    end if

    j = j + 1

    table(1:m,j) = x(1:m)

end do

close ( unit = input_unit )

return
end


subroutine r8mat_header_read ( input_filename, m, n )

!*****************************************************************************80
!
!! R8MAT_HEADER_read reads the header from an R8MAT file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FIlenAME, the name of the input file.
!
!    Output, integer M, spatial dimension.
!
!    Output, integer N, the number of points.
!
use Mod_kind_param, only : wp
implicit none

character ( len = * ) input_filename
integer m
integer n

call file_column_count ( input_filename, m )

if ( m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_HEADER_read - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data columns in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop 1
end if

call file_row_count ( input_filename, n )

if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_HEADER_read - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data rows in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop 1
end if

return
end


subroutine r8mat_uniform_01 ( m, n, seed, r )

!*****************************************************************************80
!
!! R8MAT_UNifORM_01 returns a unit pseudorandom R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of real(wp) values.
!
!    For now, the input quantity SEED is an integer variable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns
!    in the array.
!
!    Input/output, integer SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real(wp) R(M,N), the array of pseudorandom values.
!
use Mod_kind_param, only : wp
implicit none

integer m
integer n

integer i
integer j
integer k
integer seed
real(wp) r(m,n)

if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_UNifORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
end if

do j = 1, n

    do i = 1, m

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed < 0 ) then
            seed = seed + 2147483647
        end if

        r(i,j) = real ( seed, kind=wp ) * 4.656612875E-10_wp

    end do
end do

return
end


subroutine r8mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_write writes an R8MAT file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FIlenAME, the output file name.
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points.
!
!    Input, real(wp) TABLE(M,N), the table data.
!
use Mod_kind_param, only : wp
implicit none

integer m
integer n

integer j
character ( len = * ) output_filename
integer output_status
integer output_unit
character ( len = 30 ) string
real(wp) table(m,n)
!
!  Open the file.
!
call get_unit ( output_unit )

open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_write - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
        trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop 1
end if
!
!  Create a format string.
!
!  For greater precision in the output file, try:
!
!                                            '(', m, 'g', 24, '.', 16, ')'
!
write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 14, '.', 6, ')'
!
!  Write the data.
!
do j = 1, n
    write ( output_unit, string ) table(1:m,j)
end do
!
!  Close the file.
!
close ( unit = output_unit )

return
end


subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNifORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of real(wp) values.
!
!    For now, the input quantity SEED is an integer variable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer N, the number of entries in
!    the vector.
!
!    Input/output, integer SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real(wp) R(N), the vector of pseudorandom values.
!
use Mod_kind_param, only : wp
implicit none

integer n

integer i
integer k
integer seed
real(wp) r(n)

if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNifORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
end if

do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
        seed = seed + 2147483647
    end if

    r(i) = real ( seed, kind=wp ) * 4.656612875E-10_wp

end do

return
end


subroutine random_initialize ( seed )

!*****************************************************************************80
!
!! RANdoM_INITIALIZE initializes the FORTRAN90 random number seed.
!
!  Discussion:
!
!    If you don't initialize the random number generator, its behavior
!    is not specified.  If you initialize it simply by:
!
!      call random_seed ( )
!
!    its behavior is not specified.  On the DEC ALPHA, if that's all you
!    do, the same random number sequence is returned.  In order to actually
!    try to scramble up the random number generator a bit, this routine
!    goes through the tedious process of getting the size of the random
!    number seed, making up values based on the current time, and setting
!    the random number seed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer SEED.
!    If SEED is zero on input, then you're asking this routine to come up
!    with a seed value, which is returned as output.
!    If SEED is nonzero on input, then you're asking this routine to
!    use the input value of SEED to initialize the random number generator,
!    and SEED is not changed on output.
!
use Mod_kind_param, only : wp
implicit none

integer count
integer count_max
integer count_rate
logical, parameter :: debug = .false.
integer i
integer seed
integer, allocatable :: seed_vector(:)
integer seed_size
real(wp) t
!
!  Initialize the random number seed.
!
call random_seed ( )
!
!  Determine the size of the random number seed.
!
call random_seed ( size = seed_size )
!
!  Allocate a seed of the right size.
!
allocate ( seed_vector(seed_size) )

if ( seed /= 0 ) then

    if ( debug ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RANdoM_INITIALIZE'
        write ( *, '(a,i20)' ) '  Initialize RANdoM_NUMBER, user SEED = ', seed
    end if

else

    call system_clock ( count, count_rate, count_max )

    seed = count

    if ( debug ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RANdoM_INITIALIZE'
        write ( *, '(a,i20)' ) '  Initialize RANdoM_NUMBER, arbitrary SEED = ', &
            seed
    end if

end if
!
!  Now set the seed.
!
seed_vector(1:seed_size) = seed

call random_seed ( put = seed_vector(1:seed_size) )
!
!  Free up the seed space.
!
deallocate ( seed_vector )
!
!  Call the random number routine a bunch of times.
!
do i = 1, 100
    call random_number ( harvest = t )
end do

return
end


subroutine s_to_r8 ( s, r, ierror, lchar )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 from a string.
!
!  Discussion:
!
!    This routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the real number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 spaces
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon.
!
!    with most quantities optional.
!
!  Example:
!
!    S                 R
!
!    '1'               1.0
!    '     1   '       1.0
!    '1A'              1.0
!    '12,34,56'        12.0
!    '  34 7'          34.0
!    '-1E2ABCD'        -100.0
!    '-1X2ABCD'        -1.0
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0
!    '17d2'            1700.0
!    '-14e-2'         -0.14
!    'e2'              100.0
!    '-12.73e-9.23'   -12.73 * 10.0^(-9.23)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal real.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real(wp) R, the real value that was read from the string.
!
!    Output, integer IERROR, error flag.
!
!    0, no errors occurred.
!
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer LCHAR, the number of characters read from
!    the string to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
use Mod_kind_param, only : wp
implicit none

logical ch_eqi
character c
integer ierror
integer ihave
integer isgn
integer iterm
integer jbot
integer jsgn
integer jtop
integer lchar
integer nchar
integer ndig
real(wp) r
real(wp) rbot
real(wp) rexp
real(wp) rtop
character ( len = * ) s
character, parameter :: TAB = char ( 9 )

nchar = len_trim ( s )
ierror = 0
r = 0.0_wp
lchar = - 1
isgn = 1
rtop = 0.0_wp
rbot = 1.0_wp
jsgn = 1
jtop = 0
jbot = 1
ihave = 1
iterm = 0

do

    lchar = lchar + 1
    c = s(lchar+1:lchar+1)
    !
    !  Blank or TAB character.
    !
    if ( c == ' ' .or. c == TAB ) then

        if ( ihave == 2 ) then

        else if ( ihave == 6 .or. ihave == 7 ) then
            iterm = 1
        else if ( 1 < ihave ) then
            ihave = 11
        end if
        !
        !  Comma.
        !
    else if ( c == ',' .or. c == ';' ) then

        if ( ihave /= 1 ) then
            iterm = 1
            ihave = 12
            lchar = lchar + 1
        end if
        !
        !  Minus sign.
        !
    else if ( c == '-' ) then

        if ( ihave == 1 ) then
            ihave = 2
            isgn = - 1
        else if ( ihave == 6 ) then
            ihave = 7
            jsgn = - 1
        else
            iterm = 1
        end if
        !
        !  Plus sign.
        !
    else if ( c == '+' ) then

        if ( ihave == 1 ) then
            ihave = 2
        else if ( ihave == 6 ) then
            ihave = 7
        else
            iterm = 1
        end if
        !
        !  Decimal point.
        !
    else if ( c == '.' ) then

        if ( ihave < 4 ) then
            ihave = 4
        else if ( 6 <= ihave .and. ihave <= 8 ) then
            ihave = 9
        else
            iterm = 1
        end if
        !
        !  Exponent marker.
        !
    else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

        if ( ihave < 6 ) then
            ihave = 6
        else
            iterm = 1
        end if
        !
        !  Digit.
        !
    else if ( ihave < 11 .and. lge ( c, '0' ) .and. lle ( c, '9' ) ) then

        if ( ihave <= 2 ) then
            ihave = 3
        else if ( ihave == 4 ) then
            ihave = 5
        else if ( ihave == 6 .or. ihave == 7 ) then
            ihave = 8
        else if ( ihave == 9 ) then
            ihave = 10
        end if

        call ch_to_digit ( c, ndig )

        if ( ihave == 3 ) then
            rtop = 10.0_wp * rtop + real ( ndig, kind=wp )
        else if ( ihave == 5 ) then
            rtop = 10.0_wp * rtop + real ( ndig, kind=wp )
            rbot = 10.0_wp * rbot
        else if ( ihave == 8 ) then
            jtop = 10 * jtop + ndig
        else if ( ihave == 10 ) then
            jtop = 10 * jtop + ndig
            jbot = 10 * jbot
        end if
        !
        !  Anything else is regarded as a terminator.
        !
    else
        iterm = 1
    end if
    !
    !  If we haven't seen a terminator, and we haven't examined the
    !  entire string, go get the next character.
    !
    if ( iterm == 1 .or. nchar <= lchar+1 ) then
        exit
    end if

end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LCHAR is equal to NCHAR.
!
if ( iterm /= 1 .and. lchar+1 == nchar ) then
    lchar = nchar
end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then

    ierror = ihave

    return
end if
!
!  Number seems OK.  Form it.
!
if ( jtop == 0 ) then
    rexp = 1.0_wp
else

    if ( jbot == 1 ) then
        rexp = 10.0_wp**( jsgn * jtop )
    else
        rexp = jsgn * jtop
        rexp = rexp / jbot
        rexp = 10.0_wp**rexp
    end if

end if

r = isgn * rexp * rtop / rbot

return
end


subroutine s_to_r8vec ( s, n, rvec, ierror )

!*****************************************************************************80
!
!! S_TO_R8VEC reads an R8VEC from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be read.
!
!    Input, integer N, the number of values expected.
!
!    Output, real(wp) RVEC(N), the values read from the string.
!
!    Output, integer IERROR, error flag.
!    0, no errors occurred.
!    -K, could not read data for entries -K through N.
!
use Mod_kind_param, only : wp
implicit none

integer n

integer i
integer ierror
integer ilo
integer lchar
real(wp) rvec(n)
character ( len = * ) s

i = 0

ilo = 1

do while ( i < n )

    i = i + 1

    call s_to_r8 ( s(ilo:), rvec(i), ierror, lchar )

    if ( ierror /= 0 ) then
        ierror = -i
        exit
    end if

    ilo = ilo + lchar

end do

return
end


subroutine s_word_count ( s, nword )

!*****************************************************************************80
!
!! S_WORD_COUNT counts the number of "words" in a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be examined.
!
!    Output, integer NWORD, the number of "words" in the string.
!    Words are presumed to be separated by one or more blanks.
!
use Mod_kind_param, only : wp
implicit none

logical blank
integer i
integer lens
integer nword
character ( len = * ) s

nword = 0
lens = len ( s )

if ( lens <= 0 ) then
    return
end if

blank = .true.

do i = 1, lens

    if ( s(i:i) == ' ' ) then
        blank = .true.
    else if ( blank ) then
        nword = nword + 1
        blank = .false.
    end if

end do

return
end


subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
use Mod_kind_param, only : wp
implicit none

character ( len = 8 ) ampm
integer d
integer h
integer m
integer mm
character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
integer n
integer s
integer values(8)
integer y

call date_and_time ( values = values )

y = values(1)
m = values(2)
d = values(3)
h = values(5)
n = values(6)
s = values(7)
mm = values(8)

if ( h < 12 ) then
    ampm = 'AM'
else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
        ampm = 'Noon'
    else
        ampm = 'PM'
    end if
else
    h = h - 12
    if ( h < 12 ) then
        ampm = 'PM'
    else if ( h == 12 ) then
        if ( n == 0 .and. s == 0 ) then
            ampm = 'Midnight'
        else
            ampm = 'AM'
        end if
    end if
end if

write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

return
end
