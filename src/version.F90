SUBROUTINE VERSION(version_string)
      ! Retrives git information and writes to stdout
      !
      ! See here: https://git-scm.com/docs/git-describe#_examples
      !
      ! Paul Gierz
      ! AWI Bremerhaven
      IMPLICIT NONE
      CHARACTER(len=100), intent(out) :: version_string

      CALL execute_command_line("git describe --dirty  > gitversion.tmp")
      OPEN(unit=99, file="gitversion.tmp")
      READ(99,*) version_string
      WRITE(*,*) "You are using this version of FESOM"
      WRITE(*,*) "Based upon git describe --dirty): ", len_trim(version_string)

END SUBROUTINE
