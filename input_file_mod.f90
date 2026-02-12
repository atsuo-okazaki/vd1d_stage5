module input_file_mod
  use kind_params, only : i4b
  implicit none
  
  ! Name of main parameter file (default: ad1d.in)
  character(len=256) :: infile = 'ad1d.in'
  
end module input_file_mod
