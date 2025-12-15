!
! Copyright (c) 1989-2014 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
!
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
 subroutine run_diag_r(lmax,np,fp,ep,debl,lloc,irc,npx,lpx,fpx,nxtra, &
&                    vkb,evkb,nproj,rr,vfull,vp,zz,mmax)

!diagnostics for semi-local and Vanderbilt-Kleinman-bylander pseudopotentials
!fully-relativistic case

!lmax  maximum angular momentum
!np  principal quantum number for corresponding all-electron state
!ep  bound-state of scattering state reference energies for vkb potentials
!debl  energy sift for second projector
!fp  occupancy
!lloc  l for local potential
!irc  core radii indices
!npx  principal quantum number for corresponding all-electron state
!lpx  l's for extra valence states
!fpx  occupancies for extra valence states
!nxtra  number of extra valence states
!vkb  VKB projectors
!evkb  coefficients of VKB projectors
!nproj  number of vkb projectors for each l
!rr  log radial grid
!vfull  all-electron potential
!vp  semi-local pseudopotentials (vp(:,5) is local potential if linear comb.)
!zz  atomic number
!mmax  size of radial grid

 implicit none
 integer, parameter :: dp=kind(1.0d0)

!Input variables
 integer :: lmax,lloc,mmax,nlim,nxtra
 integer :: np(6),npx(6),lpx(6),irc(6),nproj(6)
 real(dp) :: zz
 real(dp) :: rr(mmax),vp(mmax,5,2),vfull(mmax),vkb(mmax,2,4,2)
 real(dp):: ep(6,2),debl(6,2),fp(6),fpx(6),evkb(2,4,2)

!Output variables - printing only

!Local variables
 integer :: ll,l1,ii,jj,ierr,mch,mchf
 integer :: ikap,mkap,kap
 real(dp) :: al,emax,emin,etest,sls,umch,upmch,cnorm,uldf,gam,gpr
 real(dp), allocatable :: uu(:),up(:),ur(:,:),urp(:,:)

 allocate(uu(mmax),up(mmax),ur(mmax,2),urp(mmax,2))

 al = 0.01d0 * dlog(rr(101)/rr(1))

! loop for diagnostic output

 write(6,'(/a)') 'Diagnostic tests using semi-local pseudopotentials'
 write(6,'(/2a)')'  l kap rcore       rmatch      e in        ', &
&   'e test    norm test  slope test'
!
 do l1 = 1, lmax + 1
  ll = l1 - 1
  if(l1==1) then
   mkap=1
  else
   mkap=2
  end if
! loop on J = ll +/- 1/2
  do ikap=1,mkap
   if(ikap==1) kap=-(ll+1)
   if(ikap==2) kap=  ll
   mchf=max(irc(l1),irc(lloc))+5
   etest=ep(l1,ikap)
   if(fp(l1)/=0.0d0) then
     call ldiracfb(np(l1),ll,kap,ierr,etest, &
&                  rr,zz,vfull,ur,urp,mmax,mch)
   else
     call ldiracfs(ll,kap,ierr,etest, &
&                  rr,zz,vfull,ur,urp,mmax,mchf)
   end if
   call renorm_r(ur,rr,ll,kap,zz,mmax,cnorm)
   umch=ur(mchf,1)
   upmch=cnorm*urp(mchf,1)
   uldf=upmch/umch
   etest=ep(l1,ikap)
   if(fp(l1)/=0.0d0) then
     call lschpb(ll+1,ll,ierr,etest, &
&                rr,vp(1,l1,ikap),uu,up,mmax,mch)
     if(ierr/=0) write(6,*) 'l1,ierr,etest',l1,ierr,etest
   else
     call lschpse(ll+1,ll,ierr,etest,uldf, &
&                 rr,vp(1,l1,ikap),uu,up,mmax,mchf)
     if(ierr/=0) write(6,*) 'l1,ierr,etest',l1,ierr,etest
   end if
   gam=dabs(umch/uu(mchf))
   gpr=dabs(upmch/up(mchf))
!
   write(6,'(2i3,6f12.7)') ll,kap,rr(irc(l1)),rr(mchf),&
&        ep(l1,ikap),etest,gam,gpr

  end do !ikap
 end do !l1
!
! diagnostic output for extra valence states in the presence of
! shallow core treated as valence
!
 if(nxtra .gt. 0) then
  write(6,'(a)') ''
  do jj = 1, nxtra
    ll = lpx(jj)
    l1 = ll + 1
   if(l1==1) then
    mkap=1
   else
    mkap=2
   end if
! loop on J = ll +/- 1/2
   do ikap=1,mkap
     if(ikap==1) kap=-(ll+1)
     if(ikap==2) kap=  ll
     mchf=max(irc(l1),irc(lloc))+5
     etest=ep(l1,ikap)+debl(l1,ikap)
     call ldiracfb(npx(jj),ll,kap,ierr,etest, &
&                  rr,zz,vfull,ur,urp,mmax,mch)
     umch=ur(mchf,1)
     upmch=urp(mchf,1)
     call lschpb(ll+2,ll,ierr,etest, &
&                rr,vp(1,l1,ikap),uu,up,mmax,mch)
     gam=dabs(umch/uu(mchf))
     gpr=dabs(upmch/up(mchf))
!
     write(6,'(2i3,6f12.7)') ll,kap,rr(irc(l1)),rr(mchf),&
&          ep(l1,ikap)+debl(l1,ikap),etest,gam,gpr
   end do !ikap
  end do !jj (nxtra)
 end if
!
! loop for diagnostic output using Vanderbilt Kleinman-Bylander projectors
!
 write(6,'(/a/a)') &
& 'Diagnostic tests using Vanderbilt-Kleinman-Bylander pseudopotentials',&
& '  1 or 2 projectors used as specified by nproj input data'
 write(6,'(/2a)')'  l kap rcore       rmatch      e in        ', &
&   'e test    norm test  slope test'
!
 do l1 = 1, lmax + 1
  ll = l1 - 1
  if(l1==1) then
   mkap=1
  else
   mkap=2
  end if
! loop on J = ll +/- 1/2
  do ikap=1,mkap
   if(ikap==1) kap=-(ll+1)
   if(ikap==2) kap=  ll
   if(ll/=lloc) then
! find cutoff radius for projectors
    mchf=max(irc(l1),irc(lloc))+5
    etest=ep(l1,ikap)
    if(fp(l1)/=0.0d0 .or. ep(l1,ikap)<0.d0) then
      call ldiracfb(np(l1),ll,kap,ierr,etest, &
&                   rr,zz,vfull,ur,urp,mmax,mch)
    else
      call ldiracfs(ll,kap,ierr,etest, &
&                  rr,zz,vfull,ur,urp,mmax,mchf)
    end if
    umch=ur(mchf,1)
    upmch=urp(mchf,1)
    uldf=upmch/umch
    etest=ep(l1,ikap)

    if(fp(l1)/=0.0d0 .or. ep(l1,ikap)<0.0d0) then
      sls=(l1-1)*l1
      emax=vp(mmax,l1,ikap)+0.5d0*sls/rr(mmax)**2
      emin=emax
      do jj=1,mmax
        emin=dmin1(emin,vp(jj,l1,ikap)+0.5d0*sls/rr(jj)**2)
      end do
      call lschvkbb(ll+1,ll,nproj(l1),ierr,etest,emin,emax, &
&                   rr,vp(1,lloc+1,1),vkb(1,1,l1,ikap),evkb(1,l1,ikap), &
&                   uu,up,mmax,mch)
      if(ierr/=0) write(6,*) 'l1,ierr,etest',l1,ierr,etest
    else
      if(etest>0.0d0) then
       emin=0.0d0
       emax=2.d0*etest
      else
       emin=2.d0*etest
       emax=0.0d0
      end if
      call lschvkbbe(ll+1,ll,nproj(l1),ierr,etest,uldf,emin,emax, &
&                  rr,vp(1,lloc+1,1),vkb(1,1,l1,ikap),evkb(1,l1,ikap), &
&                  uu,up,mmax,mchf)
      if(ierr/=0) write(6,*) 'l1,ierr,etest',l1,ierr,etest
    end if
    gam=dabs(umch/uu(mchf))
    gpr=dabs(upmch/up(mchf))
!
    write(6,'(2i3,6f12.7)') ll,kap,rr(irc(l1)),rr(mchf),&
&         ep(l1,ikap),etest,gam,gpr

   end if
  end do !ikap
 end do !l1
!
! diagnostic output for extra valence states in the presence of
! shallow core treated as valence
!
 if(nxtra .gt. 0) then
  write(6,'(a)') ''
  do jj = 1, nxtra
    ll = lpx(jj)
    l1 = ll + 1
   if(ll/=lloc .and. nproj(l1)==1) then
    if(l1==1) then
     mkap=1
    else
     mkap=2
    end if
! loop on J = ll +/- 1/2
    do ikap=1,mkap
      if(ikap==1) kap=-(ll+1)
      if(ikap==2) kap=  ll
      mchf=max(irc(l1),irc(lloc))+5
      etest=ep(l1,ikap)+debl(l1,ikap)
      call ldiracfb(npx(jj),ll,kap,ierr,etest, &
&                   rr,zz,vfull,ur,urp,mmax,mch)
      umch=ur(mchf,1)
      upmch=urp(mchf,1)
      sls=(l1-1)*l1
      emax=vp(mmax,l1,ikap)+0.5d0*sls/rr(mmax)**2
      emin=emax
      do ii=1,mmax
        emin=dmin1(emin,vp(ii,l1,ikap)+0.5d0*sls/rr(ii)**2)
      end do
      call lschvkbb(ll+2,ll,nproj(l1),ierr,etest,emin,emax, &
&                   rr,vp(1,lloc+1,1),vkb(1,1,l1,ikap),evkb(1,l1,ikap), &
&                   uu,up,mmax,mch)
      if(ierr/=0) write(6,*) 'l1,ierr,etest',l1,ierr,etest
      gam=dabs(umch/uu(mchf))
      gpr=dabs(upmch/up(mchf))
      write(6,'(2i3,6f12.7)') ll,kap,rr(irc(l1)),rr(mchf),&
&           ep(l1,ikap)+debl(l1,ikap),etest,gam,gpr
    end do !ikap
   end if !lloc, nproj==1
  end do !jj (nxtra)
 end if
!
! diagnostic output for 2nd projectors based on scattering states
!
 write(6,'(a)') ''
 do l1 = 1, lmax + 1
  ll = l1 - 1
  if(ll==lloc) cycle
  if(nproj(l1)==2) then
   if(l1==1) then
    mkap=1
   else
    mkap=2
   end if
! loop on J = ll +/- 1/2
   do ikap=1,mkap
    if(ikap==1) kap=-(ll+1)
    if(ikap==2) kap=  ll

! find cutoff radius for projectors
    mchf=max(irc(l1),irc(lloc))+5
    etest=ep(l1,ikap)+debl(l1,ikap)
    if(etest<0.0d0) then
      call ldiracfb(np(l1)+1,ll,kap,ierr,etest, &
&                    rr,zz,vfull,ur,urp,mmax,mch)
     else
       call ldiracfs(ll,kap,ierr,etest, &
&                  rr,zz,vfull,ur,urp,mmax,mchf)
    end if
    umch=ur(mchf,1)
    upmch=urp(mchf,1)
    uldf=upmch/umch
    etest=ep(l1,ikap)+debl(l1,ikap)
    if(etest<0.0d0) then
      sls=(l1-1)*l1
      emax=vp(mmax,l1,ikap)+0.5d0*sls/rr(mmax)**2
      emin=emax
      do ii=1,mmax
        emin=dmin1(emin,vp(ii,l1,ikap)+0.5d0*sls/rr(ii)**2)
      end do
      call lschvkbb(ll+2,ll,nproj(l1),ierr,etest,emin,emax, &
&                   rr,vp(1,lloc+1,1),vkb(1,1,l1,ikap),evkb(1,l1,ikap), &
&                   uu,up,mmax,mch)
      if(ierr/=0) write(6,*) 'l1,ierr,etest',l1,ierr,etest
    else
      etest=ep(l1,ikap)+debl(l1,ikap)
      emin=0.0d0
      emax=2.d0*etest
       call lschvkbbe(ll+2,ll,nproj(l1),ierr,etest,uldf,emin,emax, &
 &                  rr,vp(1,lloc+1,1),vkb(1,1,l1,ikap),evkb(1,l1,ikap), &
 &                  uu,up,mmax,mchf)
      if(ierr/=0) then
       etest=ep(l1,ikap)+debl(l1,ikap)
       emin=0.0d0
       emax=2.d0*etest
        call lschvkbbe(ll+1,ll,nproj(l1),ierr,etest,uldf,emin,emax, &
 &                   rr,vp(1,lloc+1,1),vkb(1,1,l1,ikap),evkb(1,l1,ikap), &
 &                   uu,up,mmax,mchf)
       if(ierr/=0) write(6,*) 'l1,ierr,etest',l1,ierr,etest
      end if
    end if
    gam=dabs(umch/uu(mchf))
    gpr=dabs(upmch/up(mchf))
 !
    write(6,'(2i3,6f12.7)') ll,kap,rr(irc(l1)),rr(mchf),&
 &         ep(l1,ikap)+debl(l1,ikap),etest,gam,gpr
   end do !ikap
  end if !nproj==2
 end do !l1
 
 deallocate(uu,up,ur,urp)
 return
 end subroutine run_diag_r
