MODULE COMMANDS
    implicit none

CONTAINS

SUBROUTINE check_and_create_directory(dir_name, dir_exists)
    character(len=*), intent(in) :: dir_name
    logical, intent(out) :: dir_exists

    call check_directory_exists(dir_name, dir_exists)
    if (.not. dir_exists) then
        call create_directory(dir_name)
    end if
END SUBROUTINE check_and_create_directory

SUBROUTINE check_directory_exists(dir_name, exists)
    character(len=*), intent(in) :: dir_name
    logical, intent(out) :: exists
    inquire(file=trim(dir_name), exist=exists)
END SUBROUTINE check_directory_exists

SUBROUTINE create_directory(dir_name)
    character(len=*), intent(in) :: dir_name
    integer :: ios
    character(len=255) :: command

    ! Construct the system command to create a directory
    command = 'mkdir -p ' // trim(dir_name)

    ! Execute the system command
    call execute_command(command, ios)

    if (ios /= 0) then
        print *, "Error creating directory:", dir_name
        stop 1
    end if
END SUBROUTINE create_directory

SUBROUTINE execute_command(command, ios)
    character(len=*), intent(in) :: command
    integer, intent(out) :: ios
    call system(command)
    ios = 0  ! assuming system call was successful

END SUBROUTINE execute_command


END MODULE COMMANDS
