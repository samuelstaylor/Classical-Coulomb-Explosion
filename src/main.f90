! Program:      Classical Coulomb Explosion
! Description:  Simulates molecles of given charges exploding from eachother
! Author:       Samuel S. Taylor, 2024

PROGRAM MAIN
    USE COULOMB_EXPLOSION
    implicit none

    call initialize
    call calculate
    call cleanup
    
END PROGRAM MAIN