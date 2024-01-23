!! SPDX-License-Identifier: BSD-3-Clause

!!! A few utility routines to find factors of integer numbers
module factor

   implicit none

   private

   public :: findfactor
   public :: primefactors

contains

   subroutine findfactor(num, factors, nfact)

      implicit none

      integer, intent(IN) :: num
      integer, intent(OUT), dimension(*) :: factors
      integer, intent(OUT) :: nfact
      integer :: i, m

      ! find the factors <= sqrt(num)
      ! Cast the int as double to make sure of the correct result of sqrt
      ! IntelLLVM got an issue with 1.0 but not with 1.d0
      m = int(sqrt(num * 1.d0))
      nfact = 1
      do i = 1, m
         if (num / i * i == num) then
            factors(nfact) = i
            nfact = nfact + 1
         end if
      end do
      nfact = nfact - 1

      ! derive those > sqrt(num)
      if (factors(nfact)**2 /= num) then
         do i = nfact + 1, 2 * nfact
            factors(i) = num / factors(2 * nfact - i + 1)
         end do
         nfact = nfact * 2
      else
         do i = nfact + 1, 2 * nfact - 1
            factors(i) = num / factors(2 * nfact - i)
         end do
         nfact = nfact * 2 - 1
      end if

      return

   end subroutine findfactor

   subroutine primefactors(num, factors, nfact)

      implicit none

      integer, intent(IN) :: num
      integer, intent(OUT), dimension(*) :: factors
      integer, intent(INOUT) :: nfact

      integer :: i, n

      i = 2
      nfact = 1
      n = num
      do
         if (mod(n, i) == 0) then
            factors(nfact) = i
            nfact = nfact + 1
            n = n / i
         else
            i = i + 1
         end if
         if (n == 1) then
            nfact = nfact - 1
            exit
         end if
      end do

      return

   end subroutine primefactors
end module
