MODULE mo_mrm_MDF_tools
!< date: October 2019
!< author: Maren Kaluza

  IMPLICIT NONE

  public :: get_decomposition_input

  private

CONTAINS
subroutine get_decomposition_input(iDomain, nNodes, nLinks, iStart11)
  use mo_kind, only : i4
  use mo_mrm_global_variables, only : &
          level11,      & ! IN: for number of nCells
          L11_nOutlets    ! IN: number of nodes minus number of outlets is the
                          !        number of links
  integer(i4), intent(in)  :: iDomain
  integer(i4), intent(out) :: nNodes
  integer(i4), intent(out) :: nLinks
  integer(i4), intent(out) :: iStart11

  nNodes = level11(iDomain)%nCells
  nLinks = nNodes - L11_nOutlets(iDomain)
  iStart11 = level11(iDomain)%iStart
end subroutine get_decomposition_input

END MODULE mo_mrm_MDF_tools
