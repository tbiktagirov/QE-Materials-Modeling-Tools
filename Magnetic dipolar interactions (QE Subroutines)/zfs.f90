!
! Copyright (C) 2001-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE zfs
  !-----------------------------------------------------------------------
  !
  !  
  USE kinds,                  ONLY : dp 
  USE io_global,              ONLY : stdout
  USE constants,              ONLY : fpi, angstrom_au, rytoev, bohr_radius_si
  USE fft_base,               ONLY : dffts
  USE symme,                  ONLY : symmatrix
  USE klist,                  ONLY : nks
  USE mp,                     ONLY : mp_sum
  USE gipaw_module,           ONLY : nbnd_occ, iverbosity

  !-- constants ----------------------------------------------------------
  IMPLICIT NONE
  real(dp), parameter :: mu0_by_fpi  = 1e-7
  real(dp), parameter :: mu_Bohr     = 9.27401e-24_dp
  real(dp), parameter :: bohr_radius = bohr_radius_si
  real(dp), parameter :: g_factor    = 2.002319_dp
  real(dp), parameter :: J_to_cm     = 5.03445e+22_dp
  !-- local variables ----------------------------------------------------
  real(dp) :: zfs_bare(3,3), zfs_bare_i(3,3)
  real(dp) :: zfs_gipaw_ae(3,3), zfs_gipaw_ps(3,3)
  real(dp) :: zfs_ss_tot(3,3)
  real(dp) :: v(3), axis(3,3), Dval, Eval
  integer :: alpha, beta, gamma
  real :: spin, factor
  logical :: ifound

  call start_clock('zfs')

  ! calculating a common factor, chcking spin
  spin =  0.5 * ( maxval(nbnd_occ(1:nks)) - minval(nbnd_occ(1:nks)) )
  factor =  mu0_by_fpi * g_factor**2 * mu_bohr**2 / bohr_radius ** 3 * J_to_cm
  write(stdout,*)
  write(stdout,1003) 'Effective spin is S = ', spin
  write(stdout,*)
1003 FORMAT(5X,A,F6.1)

  ! calculate the bare contribution
  write(stdout,'(5X,''Calculating the plane-wave contribution'')')
  call zfs_ss_bare(zfs_bare, zfs_bare_i)

  ! print results
  if (iverbosity > 20) then
    write(stdout,*)
    write(stdout,'(5X,''Testing/debugging only: plane-wave MAE (not symmetrized)'')')
      do beta = 1, 3
        write(stdout,1000) (factor * zfs_bare(alpha,beta) / 4, alpha=1,3)
      enddo
    write(stdout,*)
  endif

  write(stdout, '(5X,''Done with plane-wave term - calculating PAW reconstruction'')')

  ! calculate the reconstruction part 
  call zfs_correction("AE", zfs_gipaw_ae)
  call zfs_correction("PS", zfs_gipaw_ps)

  if (iverbosity > 20) then
    write(stdout,*)
    write(stdout,'(5X,''Testing/debugging only: all-electrom PAW term (not symmetrized)'')')
      do beta = 1, 3
        write(stdout,1000) (factor * zfs_gipaw_ae(alpha,beta) /4, alpha=1,3)
      enddo
    write(stdout,*)
  endif

  if (iverbosity > 20) then
    write(stdout,'(5X,''Testing/debugging only: pseudo-density PAW term (not symmetrized)'')')
      do beta = 1, 3
        write(stdout,1000) (factor * zfs_gipaw_ps(alpha,beta) /4, alpha=1,3)
      enddo
    write(stdout,*)
    write(stdout,*)
  endif

  !sum of all contributions
  zfs_ss_tot(:,:) = zfs_bare(:,:) + zfs_gipaw_ae(:,:) - zfs_gipaw_ps(:,:)

! TB: use symmatrix - not symtensor
  call symmatrix(zfs_ss_tot)
  call symmatrix(zfs_bare)
  call symmatrix(zfs_gipaw_ae)
  call symmatrix(zfs_gipaw_ae)

  if (iverbosity > 10) then
    write(stdout,*)
    write(stdout,'(5X,''Spin-spin contribution to magnetic anisotropy energy (MAE):'')')
    write(stdout,'(5X,''*** CAUTION: MAE is not zero-field splitting: ZFS = MAE/[S(S-1/2)] ***'')')
    write(stdout,*)
  endif

  if (iverbosity > 10) then
    write(stdout,'(5X,''----- MAE: plane-wave term, cm-1 -----'')')
      do beta = 1, 3
        write(stdout,1000) (factor * zfs_bare(alpha,beta) /4, alpha=1,3)
      enddo
    write(stdout,*)
  endif

  if (iverbosity > 10) then
    write(stdout,'(5X,''----- MAE: PAW reconstruction, cm-1 -----'')')
       do beta = 1, 3
         write(stdout,1000) (factor * (zfs_gipaw_ae(alpha,beta) &
                           - zfs_gipaw_ps(alpha,beta)) /4, alpha=1,3)
       enddo
    write(stdout,*)
  endif

  if (iverbosity > 10) then
    write(stdout,'(5X,''----- Total spin-spin MAE, cm-1 (it is not ZFS yet!) -----'')')
      do beta = 1, 3
        write(stdout,1000) (factor * zfs_ss_tot(alpha,beta) /4, alpha=1,3)
      enddo
    write(stdout,*)
  endif


1000 FORMAT(5X,3(ES14.7,2X))

  ! print ZFS values
  if (spin .gt. 0.5) then 
    write(stdout,*)
    write(stdout,'(5X,''Spin-spin Zero-Field Splitting with PAW reconstruction:'')')
    write(stdout,'(5X,''(Reference: T. Biktagirov, W.G. Schmidt, U. Gerstmann, PRB 97, 115135 (2018))'')')
    write(stdout,*)

    write(stdout,'(5X,''----- Spin-spin ZFS tensor, cm-1 (symmetrized) -----'')')
      do beta = 1, 3
        write(stdout,1000) (factor * zfs_ss_tot(alpha,beta)/(2*spin*(2*spin-1)), alpha=1,3)
      enddo
    write(stdout,*)

    write(stdout, '(5X,''Principal values and principal directions of the spin-spin D tensor, cm-1'')')   
    call principal_axis(zfs_ss_tot(:,:)*factor/(2*spin*(2*spin-1)), v, axis)
    write(stdout,1001) 'Dxx=', v(1), 'axis=(', axis(1:3,1), ')'
    write(stdout,1001) 'Dyy=', v(2), 'axis=(', axis(1:3,2), ')'
    write(stdout,1001) 'Dzz=', v(3), 'axis=(', axis(1:3,3), ')'
    write(stdout,*)

    ifound = .false.
    do alpha = 1, 3
      do beta = 1, 3
        if (beta .eq. alpha) cycle
        do gamma = 1, 3
          if ( (gamma .eq. beta) .or. (gamma .eq. alpha) ) cycle
          Dval = v(alpha) - 0.5*(v(beta)+v(gamma))
          Eval = 0.5*(v(gamma)-v(beta))
          if ( (Eval/Dval .lt. 1.0_dp/3) .and. (Eval/Dval .ge. 0.0_dp) .and. .not.ifound) then
            write(stdout, '(5X,''ZFS expressed as D and E parameters (0 =< E/D =< 1/3):'')')
            write(stdout,1002) 'D =', Dval, ' cm-1' 
            write(stdout,1002) 'E =', Eval, ' cm-1'
            write(stdout,1002) 'Rhombicity (E/D):', Eval/Dval, ' '
            write(stdout,*)
            write(stdout,*)
            ifound = .true.
          endif
        enddo
      enddo
    enddo
    if (.not.ifound) then
      write(stdout, '(5X,''Rhombic case (E/D = 1/3): D = 1/3Dzz is not uniquely defined.'')')
      write(stdout,*)
    endif


  else
    write(stdout, '(5X,''ZFS is not defined for S < 1'')')
    write(stdout,*)
  endif

1001 FORMAT(5X,A,F10.6,4X,A,3F10.6,A)
1002 FORMAT(5X,A,F10.6,A)

  call stop_clock('zfs')

END SUBROUTINE zfs



!-----------------------------------------------------------------------
SUBROUTINE zfs_ss_bare(zfs_bare, zfs_bare_i)
  !-----------------------------------------------------------------------
  !
  ! ... Get the charge density on the smooth grid
  !  
  USE kinds,                  ONLY : dp 
  USE mp_wave
  USE io_global,              ONLY : stdout,ionode_id
  USE mp_world,               ONLY : mpime,nproc
  USE mp,                     ONLY : mp_sum
  USE mp_global,              ONLY : inter_pool_comm, intra_pool_comm
  USE lsda_mod,               ONLY : current_spin, isk, nspin
  USE wvfct,                  ONLY : nbnd, npwx, npw, wg, g2kin, &
                                     current_k
  USE constants,              ONLY : tpi, fpi
  USE gvecw,                  ONLY : gcutw
  USE klist,                  ONLY : nks, xk, igk_k, ngk
  USE gvect,                  ONLY : ngm, ngm_g, g, gg, gstart, ig_l2g
  USE gvecs,                  ONLY : ngms
  USE wavefunctions_module,   ONLY : evc
  USE cell_base,              ONLY : tpiba2, omega
  USE io_files,               ONLY : nwordwfc, iunwfc
  USE buffers,                ONLY : get_buffer
  USE fft_base,               ONLY : dffts, dfftp
  USE fft_interfaces,         ONLY : invfft, fwfft
  USE gipaw_module,           ONLY : job, nbnd_occ, iverbosity
  !-- parameters ---------------------------------------------------------
  IMPLICIT NONE
  real(dp), intent(out) :: zfs_bare(3,3), zfs_bare_i(3,3)  
  !-- local variables ----------------------------------------------------
  complex(dp) :: psic_i(dffts%nnr), psic_j(dffts%nnr) 
  complex(dp) :: rho1(dffts%nnr), rho2(dffts%nnr),rho3(dffts%nnr)
  complex(dp), allocatable :: zfs_g(:,:,:)
  integer :: ibnd, jbnd, nk, ik, jk, alpha, beta, ig, jg, ir
  integer :: current_spin_i, current_spin_j, s_weight
  real(dp) :: fac
  integer :: gknum, igm_l2g(ngm)
  real(dp) :: tmp_g(3,ngm_g) 
  real(dp) :: trace 
  logical ::  ifound
  complex(dp):: tmp
  complex(dp) :: rho1_g(ngms),rho2_g(ngms),rho3_g(ngms) !for future self: don't do the size = npw
  complex(dp) :: rho2_g_tot(ngm_g)


  zfs_bare(:,:) = (0.d0,0.d0)
  zfs_bare_i(:,:) = (0.d0,0.d0)

  allocate( zfs_g(ngms,3,3) )
  zfs_g(:,:,:) = (0.d0,0.d0)

  !!!building a list of minus_g indices 
  tmp_g(:,:) = (0.d0,0.d0)
  igm_l2g(:) = 0.d0  
  do ig = 1, ngm
      tmp_g(:,ig_l2g(ig)) = g(:,ig)
  enddo
  call mp_sum( tmp_g , intra_pool_comm )

  gknum = 1
  do ig = 1, ngm
    ifound = .false.
    do jg = 1, ngm_g
      if ( all((g(:,ig) + tmp_g(:,jg)) == 0.d0) ) then
        igm_l2g(ig) = jg
        ifound = .true.
        exit
      endif
    enddo
    if (.not.ifound) write(stdout,*) "WARNING: there is no g minus for ig =", ig
  enddo


do nk = 1, nks/2

  do ik = nk, nk+nks/2, nks/2
     current_spin_i = isk(ik)
     do jk = nk, nk+nks/2, nks/2
       current_spin_j = isk(jk)

       if (iverbosity > 10) then
       write(stdout,1005) 'for k-points', ik, ' and', jk
       endif

       ! loop over bands
       do ibnd = 1, nbnd 
         do jbnd = ibnd, nbnd

           if ((current_spin_i .eq. current_spin_j) .and. (ibnd .eq. jbnd)) CYCLE
           if ((ibnd > nbnd_occ(ik)) .or. (jbnd > nbnd_occ(jk))) CYCLE

           zfs_g(:,:,:) = (0.d0,0.d0)
           psic_i(:)    = (0.d0,0.d0)
           psic_j(:)    = (0.d0,0.d0)
           rho1(:)      = (0.d0,0.d0)
           rho2(:)      = (0.d0,0.d0)
           rho3(:)      = (0.d0,0.d0)
           rho1_g(:)    = (0.d0,0.d0)
           rho2_g(:)    = (0.d0,0.d0)
           rho3_g(:)    = (0.d0,0.d0)

           call gk_sort(xk(1,ik), ngm, g, gcutw, npw, igk_k(1,ik), g2kin)
           call get_buffer (evc, nwordwfc, iunwfc, ik)
           psic_i(dffts%nl(igk_k(1:npw,ik))) = evc(1:npw,ibnd)
           call invfft ('Wave', psic_i, dffts)

           call gk_sort(xk(1,jk), ngm, g, gcutw, npw, igk_k(1,jk), g2kin)
           call get_buffer (evc, nwordwfc, iunwfc, jk)
           psic_j(dffts%nl(igk_k(1:npw,jk))) = evc(1:npw,jbnd)
           call invfft ('Wave', psic_j, dffts)

           ! Calculate densities
           rho1(:) = conjg(psic_i(:)) * psic_i(:) / omega
           rho2(:) = conjg(psic_j(:)) * psic_j(:) / omega
           rho3(:) = conjg(psic_i(:)) * psic_j(:) / omega

           ! transform to reciprocal space
           CALL fwfft ('Rho', rho1, dffts)
           CALL fwfft ('Rho', rho2, dffts)
           CALL fwfft ('Rho', rho3, dffts)

           rho1_g(1:npw) = rho1(dffts%nl(1:npw))
           rho2_g(1:npw) = rho2(dffts%nl(1:npw))
           rho3_g(1:npw) = rho3(dffts%nl(1:npw))

           ! building rho2 as a function of -G
           rho2_g_tot(:)=(0.d0,0.d0)
           call mergewf(rho2_g,rho2_g_tot,npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
           call splitwf(rho2_g,rho2_g_tot,npw,igm_l2g,mpime,nproc,ionode_id,intra_pool_comm)

           !different signs for Coulomb-like and exchange-like terms
           s_weight = -1
           if ( current_spin_j .eq. current_spin_i) s_weight = +1

           !common factor; also takes into account double counting
           fac = fpi * omega
           if (ibnd .eq. jbnd) fac = tpi * omega

           do ig = gstart, npw
             trace = 1.d0/3.d0 * gg(ig)
             do alpha = 1, 3
                zfs_g(ig,alpha,alpha) = - trace
                do beta = 1, 3
                      zfs_g(ig,alpha,beta) = (zfs_g(ig,alpha,beta) &
                                + g(alpha,ig)*g(beta,ig))/gg(ig) &
                                * sqrt(wg(ibnd,ik)*wg(jbnd,jk)) &
                               *( rho1_g(ig) * rho2_g(ig) &   
                                - conjg(rho3_g(ig))*rho3_g(ig)  )
                enddo
             enddo
           enddo

            do ig = gstart, ngms
               zfs_bare_i(:,:) = zfs_bare_i(:,:) + s_weight * fac * aimag(zfs_g(ig,:,:))
               zfs_bare(:,:) = zfs_bare(:,:) + s_weight * fac * real(zfs_g(ig,:,:), kind=dp)
           enddo

         enddo
       enddo

     enddo         
  enddo

enddo
deallocate( zfs_g )

1005 FORMAT(5X,A,I6,A,I6)

#ifdef __MPI
  call mp_sum( zfs_bare, intra_pool_comm ) 
#endif


END SUBROUTINE zfs_ss_bare 




!-----------------------------------------------------------------------
SUBROUTINE zfs_correction(contr, zfs_corr_tens)
  !-----------------------------------------------------------------------
  !
  !The method developed in: T. Biktagirov, W.G. Schmidt, U. Gerstmann, PRB 97, 115135 (2018)
  !Ref[1] G. Kresse and D. Joubert, Phys. Rev. B (1999)
  !Ref[2] J. Paier, R. Hirschl, M. Marsman, and G. Kresse, J. Chem. Phys (2005)
  !
  USE io_files,              ONLY : nwordwfc, iunwfc
  USE kinds,                 ONLY : dp
  USE uspp,                  ONLY : ap
  USE atom,                  ONLY : rgrid
  USE gvect,                 ONLY : g, ngm
  USE klist,                 ONLY : nks, xk, wk, igk_k
  USE cell_base,             ONLY : tpiba2
  USE ions_base,             ONLY : nat, ityp, ntyp => nsp
  USE wvfct,                 ONLY : npwx, nbnd, npw, g2kin, &
                                    current_k, wg
  USE lsda_mod,              ONLY : current_spin, nspin, isk
  USE wavefunctions_module,  ONLY : evc
  USE becmod,                ONLY : calbec
  USE constants,             ONLY : pi, fpi
  USE buffers
  USE gvecw,                 ONLY : gcutw
  USE scf,                   ONLY : rho
  USE io_global,             ONLY : stdout
  USE paw_gipaw,             ONLY : paw_recon, paw_nkb, paw_vkb, paw_becp
  USE gipaw_module,          ONLY : job, nbnd_occ
  USE mp_global,             ONLY : intra_pool_comm, inter_pool_comm
  USE mp,                    ONLY : mp_sum
  !-- parameters ---------------------------------------------------------
  IMPLICIT NONE
  real(dp), intent(out) :: zfs_corr_tens(3,3)
  character(len=2),intent(in) :: contr ! "AE"= all-electron or "PS"=pseudo
  !-- local variables ----------------------------------------------------
  integer :: j, nt, ibnd, jbnd, il1, il2, il3, il4, nk, ik, jk, nbs1, nbs2, nbs3, nbs4 
  integer :: lm, m1, m2, m3, m4, lm1, lm2, lm3, lm4, l1, l2, l3, l4 
  integer :: nrc, kkpsi, nbeta
  integer :: ijkb0, ih, jh, ph, qh, na, ikb, jkb, pkb, qkb, r_first
  integer :: current_spin_i, current_spin_j
  integer :: s_weight, delta_ij, delta_pq
  real(dp) :: al(2), ql(2), al_0(2), ql_0(2)
  complex(dp) :: bec_product 
  real(dp) :: fac
  real(dp) :: tmp1, tmp2
  real(dp), allocatable ::  wee1(:,:,:,:,:),wee2(:,:,:,:,:)
  complex(dp), allocatable :: zfs_corr(:), paw_becp_ik(:,:), paw_becp_jk(:,:)
  
  allocate( zfs_corr(9) )
  zfs_corr = 0.0_dp

  fac = 0.5_dp  ! accounts for double counting

! calculating densities
  nbeta = maxval(paw_recon(1:ntyp)%paw_nbeta)
  allocate( wee1(ntyp,nbeta,nbeta,nbeta,nbeta) )
  allocate( wee2(ntyp,nbeta,nbeta,nbeta,nbeta) )

  wee1 = 0.0_dp
  wee2 = 0.0_dp

  do nt = 1, ntyp
     do il1 = 1, paw_recon(nt)%paw_nbeta
       do il2 = 1, paw_recon(nt)%paw_nbeta
         do il3 = 1, paw_recon(nt)%paw_nbeta
           do il4 = 1, paw_recon(nt)%paw_nbeta
           call ZFS_wee(nt, il1, il2, il3, il4, contr, 1, tmp1)
           wee1(nt, il1, il2, il3, il4) = tmp1 
           call ZFS_wee(nt, il1, il2, il3, il4, contr, 2, tmp2)
           wee2(nt, il1, il2, il3, il4) = tmp2
           enddo
         enddo
       enddo
     enddo
  enddo

 !  calculate the reconstruction part
 ! loop over k-points
 do nk = 1, nks/2
  do ik = nk, nk+nks/2, nks/2
   current_spin_i = isk(ik)
   call gk_sort ( xk(1,ik), ngm, g, gcutw, npw, igk_k(:,ik), g2kin )
   call get_buffer ( evc, nwordwfc, iunwfc, ik)
   call init_gipaw_2 ( npw, igk_k(1,ik), xk(1,ik), paw_vkb )
   call calbec ( npw, paw_vkb, evc, paw_becp )
   allocate(paw_becp_ik(size(paw_becp,1),size(paw_becp,2)))
   paw_becp_ik = paw_becp

   do jk = nk, nk+nks/2, nks/2
     current_spin_j = isk(jk)

     if ( current_spin_j .eq.  current_spin_i) then
        s_weight = +1
     else
        s_weight = -1
     endif

     call gk_sort ( xk(1,jk), ngm, g, gcutw, npw, igk_k(:,jk), g2kin )
     call get_buffer ( evc, nwordwfc, iunwfc, jk)
     call init_gipaw_2 ( npw, igk_k(1,jk), xk(1,jk), paw_vkb )
     call calbec ( npw, paw_vkb, evc, paw_becp )
     allocate(paw_becp_jk(size(paw_becp,1),size(paw_becp,2)))
     paw_becp_jk = paw_becp
     
     do ibnd = 1,nbnd
      do jbnd = 1,nbnd

        if ((current_spin_i .eq. current_spin_j) .and. (ibnd .eq. jbnd)) CYCLE
        if ((ibnd > nbnd_occ(ik)) .or. (jbnd > nbnd_occ(jk))) CYCLE

        ijkb0 = 0
        do nt = 1, ntyp
           do na = 1, nat
            if ( ityp(na) == nt ) then

               do ih = 1, paw_recon(nt)%paw_nh
                ikb = ijkb0 + ih
                nbs1 = paw_recon(nt)%paw_indv(ih)
                l1 = paw_recon(nt)%paw_nhtol(ih)
                m1 = paw_recon(nt)%paw_nhtom(ih)
                lm1 = m1 + l1**2
                  do ph = 1, paw_recon(nt)%paw_nh
                   pkb = ijkb0 + ph
                   nbs3 = paw_recon(nt)%paw_indv(ph)
                   l3 = paw_recon(nt)%paw_nhtol(ph)
                   m3 = paw_recon(nt)%paw_nhtom(ph)
                   lm3 = m3 + l3**2
                   do qh = 1, paw_recon(nt)%paw_nh
                    qkb = ijkb0 + qh
                    nbs4 = paw_recon(nt)%paw_indv(qh)
                    l4 = paw_recon(nt)%paw_nhtol(qh)
                    m4 = paw_recon(nt)%paw_nhtom(qh)
                    lm4 = m4 + l4**2

                    delta_pq = 0
                    if (ph .eq. qh) delta_pq = 1

                    bec_product= paw_becp_ik(ikb,ibnd) * conjg( paw_becp_ik(ikb,ibnd) ) &
                               * paw_becp_jk(pkb,jbnd) * conjg( paw_becp_jk(qkb,jbnd) ) &
                               -(paw_becp_jk(ikb,jbnd) * conjg( paw_becp_ik(ikb,ibnd) ) &
                               * paw_becp_ik(pkb,ibnd) * conjg( paw_becp_jk(qkb,jbnd) ) )

                    do lm = 5, 9
                      zfs_corr(lm) = zfs_corr(lm) &
                                   + fac * s_weight * sqrt(wg(ibnd,ik)*wg(jbnd,jk)) &
                                   * bec_product &
                                   * ( ap(lm,lm1,lm1) * delta_pq &
                                   * wee1(nt,nbs1,nbs1,nbs3,nbs4)  &  !density for r'> r
                                   +   ap(lm,lm3,lm4)  &  
                                   * wee2(nt,nbs1,nbs1,nbs3,nbs4) ) !density for r'< r 
                    enddo
                   enddo
                  enddo
               enddo

               ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
            endif

           enddo
        enddo
        
      enddo
     enddo

     deallocate(paw_becp_jk)
   enddo  

   deallocate(paw_becp_ik)
  enddo
enddo  


  !  transform in cartesian coordinates
  zfs_corr_tens(1,1) =  sqrt(3.0_dp) * zfs_corr(8) - zfs_corr(5)
  zfs_corr_tens(2,2) = -sqrt(3.0_dp) * zfs_corr(8) - zfs_corr(5)
  zfs_corr_tens(3,3) = 2.0_dp * zfs_corr(5)
  zfs_corr_tens(1,2) = sqrt(3.0_dp) * zfs_corr(9)
  zfs_corr_tens(2,1) = zfs_corr_tens(1,2)
  zfs_corr_tens(1,3) = -zfs_corr(6) * sqrt(3.0_dp)
  zfs_corr_tens(3,1) = zfs_corr_tens(1,3)
  zfs_corr_tens(2,3) = -zfs_corr(7) * sqrt(3.0_dp)
  zfs_corr_tens(3,2) = zfs_corr_tens(2,3)
  
  zfs_corr_tens = -sqrt(4.0_dp*pi/5.0_dp)*zfs_corr_tens

  deallocate ( zfs_corr )
  deallocate ( wee1 )
  deallocate ( wee2 )
  
END SUBROUTINE zfs_correction


!----------------------------------------------------------------------------
SUBROUTINE ZFS_wee(nt, l1, l2, l3, l4, contr,flag,wee)
!----------------------------------------------------------------------------
  !Calculating the two-electron four-partial-wave integral
  !
  USE kinds,                 ONLY : dp
  USE paw_gipaw,             ONLY : paw_recon
  USE atom,                  ONLY : rgrid
  USE constants,             ONLY : fpi
  USE io_global,             ONLY : stdout
  !
  !
  IMPLICIT NONE
  !
  integer, intent(in) :: nt
  integer, intent(in) :: l1, l2, l3, l4
  character(len=2),intent(in) :: contr ! "AE"= all-electron or "PS"=pseudo 
  integer, intent(in) :: flag
  real(dp), intent(out) :: wee
  !
  integer :: r_first, nrc, i, j
  real(dp) :: w1, w2
  real(dp), allocatable :: aux1(:), aux2(:)
  real(dp) :: al(2), ql(2), al_0(2), ql_0(2)
  !
  !
   r_first = 1
   if ( abs ( rgrid(nt)%r(1) ) < 1d-8 ) r_first = 2
   nrc = paw_recon(nt)%psphi(l1)%label%nrc
   allocate(aux1(nrc))
   allocate(aux2(nrc))
   aux1 = 0.0_dp

   if (flag .eq. 1) then
      do i = r_first, nrc !loop for r
          aux2 = 0.0_dp
          do j = r_first, i !loop for r'
               if (contr == 'AE') then
               aux2(j) = paw_recon(nt)%aephi(l3)%psi(j) &
                       * paw_recon(nt)%aephi(l4)%psi(j)
               else
               aux2(j) = paw_recon(nt)%psphi(l3)%psi(j) &
                       * paw_recon(nt)%psphi(l4)%psi(j) 
               endif
          enddo
          call simpson(nrc,aux2,rgrid(nt)%rab,w1)
          if (contr == 'AE') then
          aux1(i) = paw_recon(nt)%aephi(l1)%psi(i) &
                  * paw_recon(nt)%aephi(l2)%psi(i) * w1  / rgrid(nt)%r(i)**3 
          else
          aux1(i) = paw_recon(nt)%psphi(l1)%psi(i) &
                  * paw_recon(nt)%psphi(l2)%psi(i) * w1  / rgrid(nt)%r(i)**3 
          endif
      enddo
   else
      do i = r_first, nrc   !loop for r' 
          aux2 = 0.0_dp
          do j = r_first, i   !loop for r 
               if (contr == 'AE') then
               aux2(j) = paw_recon(nt)%aephi(l1)%psi(j) &
                       * paw_recon(nt)%aephi(l2)%psi(j)
               else
               aux2(j) = paw_recon(nt)%psphi(l1)%psi(j) &
                       * paw_recon(nt)%psphi(l2)%psi(j) 
               endif
          enddo 
          call simpson(nrc,aux2,rgrid(nt)%rab,w2)
          if (contr == 'AE') then
          aux1(i) = paw_recon(nt)%aephi(l3)%psi(i) &
                  * paw_recon(nt)%aephi(l4)%psi(i) / rgrid(nt)%r(i)**3 * w2 
          else
          aux1(i) = paw_recon(nt)%psphi(l3)%psi(i) &
                  * paw_recon(nt)%psphi(l4)%psi(i) / rgrid(nt)%r(i)**3 * w2 
          endif
      enddo
   endif
   call simpson(nrc,aux1,rgrid(nt)%rab,wee)
   deallocate(aux1)
   deallocate(aux2)

END SUBROUTINE ZFS_wee

