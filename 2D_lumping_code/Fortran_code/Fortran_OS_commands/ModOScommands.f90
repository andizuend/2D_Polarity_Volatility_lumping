!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Module with subroutines and functions to carry out different operating system      *
!*   commands working both on Windows and Linux platforms. Main uses include copying of *
!*   files via commands executed in a terminal command line.                            * 
!*   Some functions require Python to be installed and available from command line on   * 
!*   the local system.                                                                  *
!*                                                                                      *
!*   :: Authors & Copyright ::                                                          *
!*   Andi Zuend                                                                         *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2020                                                            *
!*   -> latest changes: 2020-05-20                                                      *
!*                                                                                      *
!*   :: List of subroutines and functions contained in this module:                     *
!*   --------------------------------------------------------------                     *
!*   -  subroutine copy_file                                                            *
!*   -  function   isWindowsPlatform                                                    *
!*                                                                                      *
!****************************************************************************************
module ModOScommands

use Mod_kind_param, only : wp

implicit none
public

    contains
    
    !--------------------------------------------------------------------------------------                                                                                                                           
    !   This subroutine copies a file 'file_name' to 'file_name_new' via a system command.
    !   The file_name string can include a relative path on the sytem.
    !--------------------------------------------------------------------------------------
    subroutine copy_file(file_name, file_name_new)

    use ModStringFunctions

    implicit none

    character(len=*),intent(in) :: file_name, file_name_new
    !local:
    character(len = 50 +max(len(file_name_new), len(file_name))) :: command, command2
    integer :: Estat, Cstat
    logical :: isWindowsOS
    !.................................
    
    isWindowsOS = isWindowsPlatform()  !function call
    if (isWindowsOS) then
        !call execute_command_line('cd', EXITSTAT=Estat, CMDSTAT=Cstat) !display current path
        command = 'copy "'//trim(file_name) //'" "'//trim(file_name_new)//'" > NUL'
        command2 = Replace_Text(command, "/", "\") !replace forward- by backslashes for Windows commands
    else !on a LINUX OS?
        !call execute_command_line('pwd', EXITSTAT=Estat, CMDSTAT=Cstat) !display current path
        command2 = 'cp "'//trim(file_name) //'" "'//trim(file_name_new)//'" > NUL'
    endif
    call execute_command_line(trim(command2), EXITSTAT=Estat, CMDSTAT=Cstat)

    end subroutine copy_file
    !------------------------------------------------------------------------------------
    
    
    !------------------------------------------------------------------------------------
    !A logical querry function that returns 'true' if the operating system / platform is a version of Windows.
    !It makes use of a Python command and therefore requires Python to be installed; this makes the OS commmand platform-independent.
    function isWindowsPlatform()
    
    implicit none
    
    logical :: isWindowsPlatform
    character(len=20) :: cplatform
    character(len=80) :: command
    integer :: Estat, Cstat, un1
    !...............................
    
    !(1) Generate a temporary file to store command terminal data:
    open(newunit = un1, FILE = 'tempF1.txt', STATUS = 'UNKNOWN')
    close(un1)
    !(2) use a Python command on the terminal to determine the platform/OS in use:
    command = 'python -c "import platform; print(platform.system());" > tempF1.txt'
    call execute_command_line(command, EXITSTAT=Estat, CMDSTAT=Cstat)
    open(newunit = un1, FILE = 'tempF1.txt', STATUS = 'UNKNOWN')
    read(un1,*) cplatform
    close(un1, STATUS = 'DELETE')
    
    select case(cplatform(1:3))
    case('Win', 'win', 'WIN')
        isWindowsPlatform = .true.
    case default
        isWindowsPlatform = .false.
    end select
    
    end function isWindowsPlatform
    !------------------------------------------------------------------------------------

end module ModOScommands