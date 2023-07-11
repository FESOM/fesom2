module fesom_version_info_module
  private
  public fesom_git_sha
  
contains

  
  pure function fesom_git_sha() result(x)
    character(:), allocatable :: x
    ! EO args
    x = "@FESOM_GIT_SHA@"
  end function

end module
