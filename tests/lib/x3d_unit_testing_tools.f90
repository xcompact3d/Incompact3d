module x3d_unit_testing_tools

   use MPI
   use iso_c_binding

   implicit none

   logical, save :: check_error

   private
   public :: log, check_error

contains

   subroutine log(filename, l1, l2, linf, time)

      use decomp_2d, only : mytype, nrank, nproc, real_type, decomp_2d_abort

      implicit none

      character(len=*), intent(in) :: filename
      real(mytype), intent(in) :: l1, l2, linf, time

      logical :: file_found
      integer :: iounit, code
      real(mytype), dimension(nproc) :: alltime

      alltime(:) = 0.d0
      alltime(nrank+1) = time
      call MPI_ALLREDUCE(MPI_IN_PLACE,alltime,nproc,real_type,MPI_SUM,MPI_COMM_WORLD,code)
      if (code /= 0) call decomp_2d_abort(code, "MPI_ALLREDUCE")

      if (nrank == 0) then
         write(*,*) trim(filename)//" : ", l1, l2, linf, alltime
         inquire(file=filename//".dat", exist=file_found)
         if (file_found) then
            open(newunit=iounit, file=trim(filename)//".dat", form='formatted', status='old', access='append')
         else
            open(newunit=iounit, file=trim(filename)//".dat", form='formatted', status='new')
         endif
         write(iounit, *) l1, l2, linf, alltime
         close(iounit)
      endif

   end subroutine log

end module x3d_unit_testing_tools
