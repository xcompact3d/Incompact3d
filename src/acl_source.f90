
module actuator_line_source

    use decomp_2d, only: mytype, decomp_2d_abort
    use actuator_line_model_utils
    use actuator_line_model

    implicit none
    real(mytype),save :: constant_epsilon, meshFactor, thicknessFactor,chordFactor
    real(mytype),save, allocatable :: Sx(:),Sy(:),Sz(:),Sc(:),Se(:),Sh(:),Su(:),Sv(:),Sw(:),SFX(:),SFY(:),SFZ(:)
    real(mytype),save, allocatable :: Su_part(:),Sv_part(:),Sw_part(:)
    real(mytype),save, allocatable :: Snx(:),Sny(:),Snz(:),Stx(:),Sty(:),Stz(:),Ssx(:),Ssy(:),Ssz(:),Ssegm(:)
    real(mytype),save, allocatable :: A(:,:)
    logical, allocatable :: inside_the_domain(:)
    integer :: NSource
    logical, save :: rbf_interpolation=.false.
    logical, save :: pointwise_interpolation=.false.
    logical, save :: anisotropic_projection=.false.
    logical, save :: has_mesh_based_epsilon=.false.
    logical, save :: has_constant_epsilon=.false.
    public get_locations, get_forces, set_vel, initialize_actuator_source

contains

    subroutine initialize_actuator_source

    implicit none
    integer :: counter,itur,iblade,ielem,ial

    counter=0
    if (Ntur>0) then
        do itur=1,Ntur
            ! Blades
            do iblade=1,Turbine(itur)%Nblades
                do ielem=1,Turbine(itur)%Blade(iblade)%Nelem
                counter=counter+1
                end do
            end do
            ! Tower
            if(turbine(itur)%has_tower) then
                do ielem=1,Turbine(itur)%Tower%Nelem
                counter=counter+1
                end do
            endif
        end do
    endif

    if (Nal>0) then
        do ial=1,Nal
            do ielem=1,actuatorline(ial)%NElem
                counter=counter+1
            end do
        end do
    endif
    NSource=counter
    allocate(Sx(NSource),Sy(NSource),Sz(NSource),Sc(Nsource),Su(NSource),Sv(NSource),Sw(NSource),Se(NSource),Sh(NSource),Sfx(NSource),Sfy(NSource),Sfz(NSource))
    allocate(Su_part(NSource),Sv_part(NSource),Sw_part(NSource))
    allocate(Snx(NSource),Sny(NSource),Snz(NSource),Stx(Nsource),Sty(NSource),Stz(NSource),Ssx(NSource),Ssy(NSource),Ssz(NSource),Ssegm(NSource))
    allocate(A(NSource,NSource))
    allocate(inside_the_domain(NSource))

    end subroutine initialize_actuator_source

    subroutine get_locations

    implicit none
    integer :: counter,itur,iblade,ielem,ial

    counter=0
    if (Ntur>0) then
        do itur=1,Ntur
            do iblade=1,Turbine(itur)%Nblades
                do ielem=1,Turbine(itur)%Blade(iblade)%Nelem
                counter=counter+1
                Sx(counter)=Turbine(itur)%Blade(iblade)%PEX(ielem)
                Sy(counter)=Turbine(itur)%Blade(iblade)%PEY(ielem)
                Sz(counter)=Turbine(itur)%Blade(iblade)%PEZ(ielem)
                Sc(counter)=Turbine(itur)%Blade(iblade)%EC(ielem)
                Snx(counter)=Turbine(itur)%Blade(iblade)%nEx(ielem)
                Sny(counter)=Turbine(itur)%Blade(iblade)%nEy(ielem)
                Snz(counter)=Turbine(itur)%Blade(iblade)%nEz(ielem)
                Stx(counter)=Turbine(itur)%Blade(iblade)%tEx(ielem)
                Sty(counter)=Turbine(itur)%Blade(iblade)%tEy(ielem)
                Stz(counter)=Turbine(itur)%Blade(iblade)%tEz(ielem)
                Ssx(counter)=Turbine(itur)%Blade(iblade)%sEx(ielem)
                Ssy(counter)=Turbine(itur)%Blade(iblade)%sEy(ielem)
                Ssz(counter)=Turbine(itur)%Blade(iblade)%sEz(ielem)
                Ssegm(counter)=Turbine(itur)%Blade(iblade)%EDS(ielem)

                end do
            end do
                !Tower
                if(turbine(itur)%has_tower) then
                do ielem=1,Turbine(itur)%Tower%Nelem
                counter=counter+1
                Sx(counter)=Turbine(itur)%Tower%PEX(ielem)
                Sy(counter)=Turbine(itur)%Tower%PEY(ielem)
                Sz(counter)=Turbine(itur)%Tower%PEZ(ielem)
                Sc(counter)=Turbine(itur)%Tower%EC(ielem)
                Snx(counter)=Turbine(itur)%Tower%nEx(ielem)
                Sny(counter)=Turbine(itur)%Tower%nEy(ielem)
                Snz(counter)=Turbine(itur)%Tower%nEz(ielem)
                Stx(counter)=Turbine(itur)%Tower%tEx(ielem)
                Sty(counter)=Turbine(itur)%Tower%tEy(ielem)
                Stz(counter)=Turbine(itur)%Tower%tEz(ielem)
                Ssx(counter)=Turbine(itur)%Tower%sEx(ielem)
                Ssy(counter)=Turbine(itur)%Tower%sEy(ielem)
                Ssz(counter)=Turbine(itur)%Tower%sEz(ielem)
                Ssegm(counter)=Turbine(itur)%Tower%EDS(ielem)
                end do
                endif
        end do
    endif

    if (Nal>0) then
        do ial=1,Nal
            do ielem=1,actuatorline(ial)%NElem
                counter=counter+1
                Sx(counter)=actuatorline(ial)%PEX(ielem)
                Sy(counter)=actuatorline(ial)%PEY(ielem)
                Sz(counter)=actuatorline(ial)%PEZ(ielem)
                Sc(counter)=actuatorline(ial)%EC(ielem)
                Snx(counter)=actuatorline(ial)%nEx(ielem)
                Sny(counter)=actuatorline(ial)%nEy(ielem)
                Snz(counter)=actuatorline(ial)%nEz(ielem)
                Stx(counter)=actuatorline(ial)%tEx(ielem)
                Sty(counter)=actuatorline(ial)%tEy(ielem)
                Stz(counter)=actuatorline(ial)%tEz(ielem)
                Ssx(counter)=actuatorline(ial)%sEx(ielem)
                Ssy(counter)=actuatorline(ial)%sEy(ielem)
                Ssz(counter)=actuatorline(ial)%sEz(ielem)
                Ssegm(counter)=actuatorline(ial)%EDS(ielem)
            end do
        end do
    endif

    end subroutine get_locations

    subroutine set_vel

    implicit none
    integer :: counter,itur,iblade,ielem,ial

    counter=0
    if (Ntur>0) then
        do itur=1,Ntur
            ! Blades
            do iblade=1,Turbine(itur)%Nblades
                do ielem=1,Turbine(itur)%Blade(iblade)%Nelem
                counter=counter+1
                Turbine(itur)%Blade(iblade)%EVX(ielem)=Su(counter)
                Turbine(itur)%Blade(iblade)%EVY(ielem)=Sv(counter)
                Turbine(itur)%Blade(iblade)%EVZ(ielem)=Sw(counter)
                Turbine(itur)%Blade(iblade)%Eepsilon(ielem)=Se(counter)
                end do
            end do
            ! Tower
            if(turbine(itur)%has_tower) then
                do ielem=1,Turbine(itur)%Tower%Nelem
                counter=counter+1
                Turbine(itur)%Tower%EVX(ielem)=Su(counter)
                Turbine(itur)%Tower%EVY(ielem)=Sv(counter)
                Turbine(itur)%Tower%EVZ(ielem)=Sw(counter)
                Turbine(itur)%Tower%Eepsilon(ielem)=Se(counter)
                end do
            endif
        end do
    endif

    if (Nal>0) then
        do ial=1,Nal
            do ielem=1,actuatorline(ial)%NElem
                counter=counter+1
                actuatorline(ial)%EVX(ielem)=Su(counter)
                actuatorline(ial)%EVY(ielem)=Sv(counter)
                actuatorline(ial)%EVZ(ielem)=Sw(counter)
                actuatorline(ial)%Eepsilon(ielem)=Se(counter)
            end do
        end do
    endif


    end subroutine set_vel

    subroutine get_forces

    implicit none
    integer :: counter,itur,iblade,ielem,ial

    counter=0
    if (Ntur>0) then
        do itur=1,Ntur
            !Blade
            do iblade=1,Turbine(itur)%Nblades
                do ielem=1,Turbine(itur)%Blade(iblade)%Nelem
                counter=counter+1
                Sfx(counter)=Turbine(itur)%Blade(iblade)%EFX(ielem)
                Sfy(counter)=Turbine(itur)%Blade(iblade)%EFY(ielem)
                Sfz(counter)=Turbine(itur)%Blade(iblade)%EFZ(ielem)
                end do
            end do

            !Tower
            if(turbine(itur)%has_tower) then
                do ielem=1,Turbine(itur)%Tower%Nelem
                counter=counter+1
                Sfx(counter)=Turbine(itur)%Tower%EFX(ielem)
                Sfy(counter)=Turbine(itur)%Tower%EFY(ielem)
                Sfz(counter)=Turbine(itur)%Tower%EFZ(ielem)
                end do
            endif
        end do
    endif

    if (Nal>0) then
        do ial=1,Nal
            do ielem=1,actuatorline(ial)%NElem
                counter=counter+1
                Sfx(counter)=actuatorline(ial)%EFX(ielem)
                Sfy(counter)=actuatorline(ial)%EFY(ielem)
                Sfz(counter)=actuatorline(ial)%EFZ(ielem)
            end do
        end do
    endif

    end subroutine get_forces

    subroutine Compute_Momentum_Source_Term_pointwise

        use decomp_2d, only: mytype, nproc, xstart, xend, xsize, update_halo
        use MPI
        use param, only: dx,dy,dz,eps_factor,xnu,istret,xlx,yly,zlz
        use var, only: ux1, uy1, uz1, FTx, FTy, FTz, yp

        implicit none
        real(mytype), allocatable, dimension(:,:,:) :: ux1_halo, uy1_halo, uz1_halo
        real(mytype) :: xmesh, ymesh,zmesh
        real(mytype) :: dist, epsilon, Kernel
        real(mytype) :: min_dist, ymax,ymin,zmin,zmax
        real(mytype) :: x0,y0,z0,x1,y1,z1,x,y,z,u000,u100,u001,u101,u010,u110,u011,u111
        real(mytype) :: t1,t2, alm_proj_time
        integer :: min_i,min_j,min_k
        integer :: i_lower, j_lower, k_lower, i_upper, j_upper, k_upper
        integer :: i,j,k, isource, code

        ! First we need to compute the locations
        call get_locations

        ! Zero the velocities
        Su(:)=0.0
        Sv(:)=0.0
        Sw(:)=0.0
        ! This is not optimum but works
        ! Define the domain

        if (istret.eq.0) then
        ymin=(xstart(2)-1)*dy-dy/2. ! Add -dy/2.0 overlap
        ymax=(xend(2)-1)*dy+dy/2.   ! Add +dy/2.0 overlap
        else
        ymin=yp(xstart(2))
        ymax=yp(xend(2))
        endif

        zmin=(xstart(3)-1)*dz-dz/2. ! Add a -dz/2.0 overlap
        zmax=(xend(3)-1)*dz+dz/2.   ! Add a +dz/2.0 overlap

        ! Check if the points lie outside the fluid domain
        do isource=1,Nsource
        if((Sx(isource)>xlx).or.(Sx(isource)<0).or.(Sy(isource)>yly).or.(Sy(isource)<0).or.(Sz(isource)>zlz).or.(Sz(isource)<0)) then
            print *, 'Point outside the fluid domain'
            stop
        endif
        enddo
        !write(*,*) 'Rank=', nrank, 'X index Limits=', xstart(1), xend(1), 'X lims=', (xstart(1)-1)*dx, (xend(1)-1)*dx
        !write(*,*) 'Rank=', nrank, 'Y index Limits=', xstart(2), xend(2), 'Y lims=', ymin, ymax
        !write(*,*) 'Rank=', nrank, 'Z index Limits=', xstart(3), xend(3), 'Z lims=', zmin, zmax
        call update_halo(ux1,ux1_halo,1,opt_global=.true.)
        call update_halo(uy1,uy1_halo,1,opt_global=.true.)
        call update_halo(uz1,uz1_halo,1,opt_global=.true.)
        !print *,  nrank, shape(ux1), shape(ux1_halo)
        !do j=xstart(2),xend(2)
        !print *, nrank, ux1(xstart(1),j,xstart(3)), ux1_halo(xstart(1),j,xstart(3))
        !enddo

        do isource=1,NSource

        min_dist=1e6
        if((Sy(isource)>=ymin).and.(Sy(isource)<ymax).and.(Sz(isource)>=zmin).and.(Sz(isource)<zmax)) then
            !write(*,*) 'nrank= ',nrank, 'owns this node'
            do k=xstart(3),xend(3)
            zmesh=(k-1)*dz
            do j=xstart(2),xend(2)
            if (istret.eq.0) ymesh=(j-1)*dy
            if (istret.ne.0) ymesh=yp(j)
            do i=xstart(1),xend(1)
            xmesh=(i-1)*dx
            dist = sqrt((Sx(isource)-xmesh)**2.+(Sy(isource)-ymesh)**2.+(Sz(isource)-zmesh)**2.)

            if (dist<min_dist) then
                min_dist=dist
                min_i=i
                min_j=j
                min_k=k
            endif

            enddo
            enddo
            enddo

            if(Sy(isource)>ymax.or.Sy(isource)<ymin) then
            write(*,*) 'In processor ', nrank
            write(*,*) 'Sy =', Sy(isource),'is not within the', ymin, ymax, 'limits'
            stop
            endif
            if(Sz(isource)>zmax.or.Sz(isource)<zmin) then
            write(*,*) 'In processor ', nrank
            write(*,*) 'Sz =', Sz(isource),'is not within the', zmin, zmax, 'limits'
            stop
            endif

            if(Sx(isource)>(min_i-1)*dx) then
                i_lower=min_i
                i_upper=min_i+1
            else if(Sx(isource)<(min_i-1)*dx) then
                i_lower=min_i-1
                i_upper=min_i
            else if(Sx(isource)==(min_i-1)*dx) then
                i_lower=min_i
                i_upper=min_i
            endif

            if (istret.eq.0) then
            if(Sy(isource)>(min_j-1)*dy.and.Sy(isource)<(xend(2)-1)*dy) then
                j_lower=min_j
                j_upper=min_j+1
            else if(Sy(isource)>(min_j-1)*dy.and.Sy(isource)>(xend(2)-1)*dy) then
                j_lower=min_j
                j_upper=min_j+1 ! THis is in the Halo domain
            else if(Sy(isource)<(min_j-1)*dy.and.Sy(isource)>(xstart(2)-1)*dy) then
                j_lower=min_j-1
                j_upper=min_j
            else if(Sy(isource)<(min_j-1)*dy.and.Sy(isource)<(xstart(2)-1)*dy) then
                j_lower=min_j-1 ! THis is in the halo domain
                j_upper=min_j
            else if (Sy(isource)==(min_j-1)*dy) then
                j_lower=min_j
                j_upper=min_j
            endif
            else

            if(Sy(isource)>yp(min_j)) then
                j_lower=min_j
                j_upper=min_j+1
            else if(Sy(isource)<yp(min_j)) then
                j_lower=min_j-1
                j_upper=min_j
            else if (Sy(isource)==yp(min_j)) then
                j_lower=min_j
                j_upper=min_j
            endif
            endif

            if(Sz(isource)>(min_k-1)*dz.and.Sz(isource)<(xend(3)-1)*dz) then
                k_lower=min_k
                k_upper=min_k+1
            elseif(Sz(isource)>(min_k-1)*dz.and.Sz(isource)>(xend(3)-1)*dz) then
                k_lower=min_k
                k_upper=min_k+1 ! This in the halo doamin
            else if(Sz(isource)<(min_k-1)*dz.and.Sz(isource)>(xstart(3)-1)*dz) then
                k_lower=min_k-1
                k_upper=min_k
            else if(Sz(isource)<(min_k-1)*dz.and.Sz(isource)<(xstart(3)-1)*dz) then
                k_lower=min_k-1 ! This is in the halo domain
                k_upper=min_k
            else if (Sz(isource)==(min_k-1)*dz) then
                k_lower=min_k
                k_upper=min_k
            endif

            ! Prepare for interpolation
            x0=(i_lower-1)*dx
            x1=(i_upper-1)*dx
            if (istret.eq.0) then
            y0=(j_lower-1)*dy
            y1=(j_upper-1)*dy
            else
            y0=yp(j_lower)
            y1=yp(j_upper)
            endif
            z0=(k_lower-1)*dz
            z1=(k_upper-1)*dz

            x=Sx(isource)
            y=Sy(isource)
            z=Sz(isource)

            !if(x>x1.or.x<x0.or.y>y1.or.y<y0.or.z>z1.or.z<z0) then
            !    write(*,*) 'x0, x1, x', x0, x1, x
            !    write(*,*) 'y0, y1, y', y0, y1, y
            !    write(*,*) 'z0, z1, z', z0, z1, z
            !    write(*,*) 'Problem with the trilinear interpolation';
            !    stop
            !endif
            ! Apply interpolation kernels from 8 neighboring nodes

            Su_part(isource)= trilinear_interpolation(x0,y0,z0, &
                                                  x1,y1,z1, &
                                                  x,y,z, &
                                                  ux1_halo(i_lower,j_lower,k_lower), &
                                                  ux1_halo(i_upper,j_lower,k_lower), &
                                                  ux1_halo(i_lower,j_lower,k_upper), &
                                                  ux1_halo(i_upper,j_lower,k_upper), &
                                                  ux1_halo(i_lower,j_upper,k_lower), &
                                                  ux1_halo(i_upper,j_upper,k_lower), &
                                                  ux1_halo(i_lower,j_upper,k_upper), &
                                                  ux1_halo(i_upper,j_upper,k_upper))

             Sv_part(isource)= trilinear_interpolation(x0,y0,z0, &
                                                  x1,y1,z1, &
                                                  x,y,z, &
                                                  uy1_halo(i_lower,j_lower,k_lower), &
                                                  uy1_halo(i_upper,j_lower,k_lower), &
                                                  uy1_halo(i_lower,j_lower,k_upper), &
                                                  uy1_halo(i_upper,j_lower,k_upper), &
                                                  uy1_halo(i_lower,j_upper,k_lower), &
                                                  uy1_halo(i_upper,j_upper,k_lower), &
                                                  uy1_halo(i_lower,j_upper,k_upper), &
                                                  uy1_halo(i_upper,j_upper,k_upper))

            Sw_part(isource)= trilinear_interpolation(x0,y0,z0, &
                                                  x1,y1,z1, &
                                                  x,y,z, &
                                                  uz1_halo(i_lower,j_lower,k_lower), &
                                                  uz1_halo(i_upper,j_lower,k_lower), &
                                                  uz1_halo(i_lower,j_lower,k_upper), &
                                                  uz1_halo(i_upper,j_lower,k_upper), &
                                                  uz1_halo(i_lower,j_upper,k_lower), &
                                                  uz1_halo(i_upper,j_upper,k_lower), &
                                                  uz1_halo(i_lower,j_upper,k_upper), &
                                                  uz1_halo(i_upper,j_upper,k_upper))

        else
            Su_part(isource)=0.0
            Sv_part(isource)=0.0
            Sw_part(isource)=0.0
            !write(*,*) 'Warning: I do not own this node'
        endif
        enddo
        !$OMP END PARALLEL DO

        call MPI_ALLREDUCE(Su_part,Su,Nsource,MPI_REAL8,MPI_SUM, &
            MPI_COMM_WORLD,code)
        if (code.ne.0) call decomp_2d_abort(code, "MPI_ALLREDUCE")
        call MPI_ALLREDUCE(Sv_part,Sv,Nsource,MPI_REAL8,MPI_SUM, &
            MPI_COMM_WORLD,code)
        if (code.ne.0) call decomp_2d_abort(code, "MPI_ALLREDUCE")
        call MPI_ALLREDUCE(Sw_part,Sw,Nsource,MPI_REAL8,MPI_SUM, &
            MPI_COMM_WORLD,code)
        if (code.ne.0) call decomp_2d_abort(code, "MPI_ALLREDUCE")

        ! Zero the Source term at each time step
        FTx(:,:,:)=0.0
        FTy(:,:,:)=0.0
        FTz(:,:,:)=0.0
        Visc=xnu
        !## Send the velocities to the
        call set_vel
        !## Compute forces
        call actuator_line_model_compute_forces
        !## Get Forces
        call get_forces

        if(nrank==0) then
            write(*,*) 'Projecting the AL Momentum Source term ... '
        endif
        t1 = MPI_WTIME()


            !## Add the source term
            do k=1,xsize(3)
            zmesh=(k+xstart(3)-1-1)*dz
            do j=1,xsize(2)
            if (istret.eq.0) ymesh=(xstart(2)+j-1-1)*dy
            if (istret.ne.0) ymesh=yp(xstart(2)+j-1)
            do i=1,xsize(1)
            xmesh=(i-1)*dx

            do isource=1,NSource

            dist = sqrt((Sx(isource)-xmesh)**2+(Sy(isource)-ymesh)**2+(Sz(isource)-zmesh)**2)
            epsilon=eps_factor*(dx*dy*dz)**(1.0/3.0)
            if (dist<10.0*epsilon) then
                Kernel= 1.0/(epsilon**3.0*pi**1.5)*dexp(-(dist/epsilon)**2.0)
            else
                Kernel=0.0
            endif
            ! First apply a constant lift to induce the
            FTx(i,j,k)=FTx(i,j,k)-SFx(isource)*Kernel
            FTy(i,j,k)=FTy(i,j,k)-SFy(isource)*Kernel
            FTz(i,j,k)=FTz(i,j,k)-SFz(isource)*Kernel

            enddo

            enddo
            enddo
            enddo

        alm_proj_time=MPI_WTIME()-t1
        call MPI_ALLREDUCE(alm_proj_time,t1,1,MPI_REAL8,MPI_SUM, &
                   MPI_COMM_WORLD,code)
        if (code.ne.0) call decomp_2d_abort(code, "MPI_ALLREDUCE")

        if(nrank==0) then
            alm_proj_time=alm_proj_time/float(nproc)
            write(*,*) 'AL Momentum Source term projection completed in :', alm_proj_time ,'seconds'
        endif

    end subroutine Compute_Momentum_Source_Term_pointwise


end module actuator_line_source
