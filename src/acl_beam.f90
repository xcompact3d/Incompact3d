module actuator_line_beam_model

    use decomp_2d, only: mytype, nrank
    use actuator_line_element
    use xbeam_shared

type BeamType
    character(len=100) :: name                                    ! Beam model name (the same as the turbine name)
    real(mytype), allocatable :: pos(:,:)                         ! Positions
    real(mytype), allocatable :: StructuralTwist(:)               ! Structural twist
    real(mytype), allocatable :: frame_of_reference_delta(:,:,:)  ! Frame of reference delta
    integer, allocatable :: Conn(:,:)                             ! Connectivity
    integer :: Nnodes                                             ! Number of nodes
    integer :: Ndofs                                              ! Number of degrees of freedom=6*(num_nodes-3) (blades are clamped at the root)
    integer :: NElems                                             ! Number of elements
    type(xbelem), allocatable :: elem(:)                          ! Element information
    type(xbnode), allocatable :: node(:)                          ! Nodal information
end type BeamType

contains

    !*******************************************************************************
    !
    subroutine actuator_line_beam_model_init(beam,acl,Nblades)
    !
    !*******************************************************************************

      implicit none
      type(BeamType) :: beam
      integer, intent(in) :: Nblades
      type(ActuatorLineType), intent(in), dimension(3) :: acl(3)
      integer :: i,j,k

      ! Degrees of freedom for the beam. It should be equal to the number of blades times twice the number of elements plus one (midpoints + edges)
      beam%Nnodes=Nblades*(2*acl(1)%Nelem+1)
      beam%NElems=Nblades*acl(1)%Nelem
      beam%Ndofs=6*(beam%Nnodes-Nblades)

      ! First allocate
      allocate(beam%pos(beam%Nnodes,3))
      allocate(beam%StructuralTwist(beam%Nnodes))
      allocate(beam%frame_of_reference_delta(Nblades*acl(1)%Nelem,3,3))
      allocate(beam%Conn(acl(1)%Nelem,3))
      allocate(beam%elem(beam%NElems))
      allocate(beam%node(beam%Nnodes))

      ! Init the coordinates
      do i=1,Nblades
         do j=1,acl(i)%Nelem
            beam%pos((i-1)*(2*acl(i)%NElem+1)+2*j-1,1)=acl(i)%QCx(j)    ! First-point of the element
            beam%pos((i-1)*(2*acl(i)%NElem+1)+2*j,1)=acl(i)%PEx(j)      ! Mid-point of the element
            beam%pos((i-1)*(2*acl(i)%NElem+1)+2*j+1,1)=acl(i)%QCx(j+1)  ! Last-point of the element
            beam%pos((i-1)*(2*acl(i)%NElem+1)+2*j-1,2)=acl(i)%QCy(j)    ! First-point of the element
            beam%pos((i-1)*(2*acl(i)%NElem+1)+2*j,2)=acl(i)%PEy(j)      ! Mid-point of the element
            beam%pos((i-1)*(2*acl(i)%NElem+1)+2*j+1,2)=acl(i)%QCy(j+1)  ! Last-point of the element
            beam%pos((i-1)*(2*acl(i)%NElem+1)+2*j-1,3)=acl(i)%QCz(j)    ! First-point of the element
            beam%pos((i-1)*(2*acl(i)%NElem+1)+2*j,3)=acl(i)%PEz(j)      ! Mid-point of the element
            beam%pos((i-1)*(2*acl(i)%NElem+1)+2*j+1,3)=acl(i)%QCz(j+1)  ! Last-point of the element
         enddo
      enddo
      beam%StructuralTwist=0.
          
      return

    end subroutine actuator_line_beam_model_init

    !*******************************************************************************
    !
    subroutine actuator_line_beam_solve(beam,dt)
    !
    !*******************************************************************************

      use cbeam3_solv ! This loads the module from xbeam (WInc3D needs to be compiled against the xbeam library)

      implicit none
      type(BeamType), intent(inout) :: beam
      real(mytype), intent(in) :: dt

      ! Update location of the beams

      ! Update the velocities

      ! Call the cbeam solver (timestep variation -- not sure if that is the correct thing to do)
      ! call cbeam3_solv_nlndyn_step(beam%Ndofs, &
      !                              beam%NElems,&
      !                              beam%Nnodes,&
      !                              dt,&
      !                              elem,&
      !                              node,&
      !                              static_forces,&
      !                              dynamic_forces,&
      !                              gravity_forces,&
      !                              quat,&
      !                              for_vel,&
      !                              for_acc,&
      !                              pos_ini,&
      !                              psi_ini,&
      !                              pos_def,&
      !                              psi_def,&
      !                              pos_dot_def,&
      !                              psi_dot_def,&
      !                              options)

      return

    end subroutine actuator_line_beam_solve

end module actuator_line_beam_model
