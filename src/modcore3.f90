!
! Copyright (c) 1989-2019 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
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
! Creates monotonic polynomial model core charge matching all-electron
! core charge and 4 derivatives at "crossover" radius.
! Polynomial is 8th-order with no linear term.

! Performs analysis and based on "hardness" criterion described in
! Teter, Phys. Rev. B 48, 5031 (1993) , Appendix, as

subroutine modcore3(icmod,rhops,rhotps,rhoc,rhoae,rhotae,rhomod, &
&                   fcfact,rcfact,mmax,rr,nc,nv,la,zion,iexc)

!icmod  3 coefficient optimizaion, 4 for specivied fcfact and rfact
!rhops  state-by-state pseudocharge density
!rhotps  total pseudocharge density
!rhoc  core-charge density
!rhoae  state-by-state all-electron valence charge density
!rhotae  total all-electron valence charge density
!rhomod  model core density and 4 derivatives
!fcfact  prefactor for model amplitude (multiplies crossover value)
!rcfact  prefactor for model scale (multiplies crossover radius)
!irps  rr index of maximum rc
!mmax  dimension of log grid
!rr log radial grid
!nc  number of core states
!nv  number of valence states
!la  angular-momenta
!zion  ion charge
!iexc  exchange-correlation function to be used

   implicit none
   integer, parameter :: dp=kind(1.0d0)

!Input variables
   integer :: icmod,nv,nc,iexc,mmax
   integer :: la(30)
   real(dp) :: rhoae(mmax,nv),rhops(mmax,nv),rhotae(mmax)
   real(dp) :: rhotps(mmax),rhoc(mmax),rr(mmax)
   real(dp) :: zion,fcfact,rcfact

!Output variables
   real(dp) :: rhomod(mmax,5)

!convergence criterion
   real(dp), parameter :: eps=1.0d-7
   real(dp), parameter :: blend=2.0d0

!Local variables
   real(dp) :: al
   real(dp) :: d2mdiff,rmatch,rhocmatch,r0,rcross
   real(dp) :: gg,tt,yy
   real(dp) :: drint,rtst,rint(20),fint(20)  !ad-hoc smoothing variables
   real(dp), allocatable :: vxcae(:),vxcpsp(:),vo(:),d2excae(:,:),d2excps(:,:)
   real(dp), allocatable :: dvxcae(:,:),dvxcps(:,:),vxct(:)
   integer :: ii,ircc,ircross,irmod,jj,kk
   integer :: iint  !ad-hoc smoothing variables

!2-dimensional Nelder-Mead variables
   real(dp), parameter :: alpha_nm=1.0d0
   real(dp), parameter :: gamma_nm=2.0d0
   real(dp), parameter :: rho_nm=-0.5d0
   real(dp), parameter :: sigma_nm=0.5d0
   real(dp) :: xx(2,3),ff(3),xt(2),ft,x0(2),xr(2),fr,xe(2),fe,xc(2),fc
   real(dp) :: d2ref,fta(10),d2min


   allocate(vxcae(mmax),vxcpsp(mmax),vo(mmax))
   allocate(dvxcae(mmax,nv),dvxcps(mmax,nv),vxct(mmax))
   allocate(d2excae(nv,nv),d2excps(nv,nv))

   d2excae(:,:)=0.0d0
   d2excps(:,:)=0.0d0

!set limit for d2Exc calculation to radius beyond which rhomod=rhoc
   irmod=mmax

   write(6,'(/a/a)') 'Model core correction analysis',&
   &                  '  based on d2Exc/dn_idn_j'

   call der2exc(rhotae,rhoc,rhoae,rr,d2excae,d2excps,d2mdiff, &
   &                   zion,iexc,nc,nv,la,irmod,mmax)

   write(6,'(/a/)') 'd2excae - all-electron derivatives'
   do kk=1,nv
      write(6,'(1p,4d16.6)') (d2excae(kk,jj),jj=1,nv)
   end do


! set model charge to zero
   rhomod(:,:)=0.0d0

! compute d2excps with no core correction
   rhomod(:,:)=0.0d0

   call der2exc(rhotps,rhomod(1,1),rhops,rr,d2excps,d2excae,d2mdiff, &
   &                   zion,iexc,nc,nv,la,irmod,mmax)

   write(6,'(/a/)') 'd2excps - pseudofunction derivatives with no core correction'
   do kk=1,nv
      write(6,'(1p,4d16.6)') (d2excps(kk,jj),jj=1,nv)
   end do
   write(6,'(/a,1p,e16.6)') 'rms 2nd derivative error',d2mdiff

   d2ref=d2mdiff

! find valence pseudocharge - core charge crossover
   ircc = 0
   do ii = mmax,1,-1
      if(rhoc(ii) > rhotps(ii)) then
         ircc=ii
         rmatch=rr(ircc)
         rhocmatch=rhoc(ircc)
         exit
      end if
   end do

   if(ircc == 0) then
      write(6,'(/a)') 'modcore3: ERROR ircc (core-valence charge crossover) &
      &        not found'
      stop
   else
      write(6,'(/a,2f10.4)') 'rmatch, rhocmatch',rmatch,rhocmatch
   end if

   gg=0.d0
   yy=0.d0

!option for optimization
   if(icmod==4) then

      fcfact=1.0d0

!Coarse-grid search for minimum

      write(6,'(/a)') 'Coarse scan for minimum error'
      write(6,'(a)') '  matrix elements: rms 2nd-derivative errors (mHa)'
      write(6,'(a)') '  column index : amplitude prefactor to rhocmatch'
      write(6,'(a)') '  row index : scale prefactor to rmatch'

      d2min=10.0d0
      write(6,'(/7x,10f7.3/)') (1.5d0+0.5d0*(jj-1), jj=1,10)
      do  kk=1,10
         xt(2)=(1.0d0+0.1d0*(kk-1))*rmatch
         do jj=1,10
            xt(1)=(1.5d0+0.5d0*(jj-1))*rhocmatch

            r0=1.5d0*xt(2)
            do ii=mmax,1,-1
               if(rr(ii)<r0) then
                  call gg1cc(gg,yy)
                  if(xt(1)*gg<rhoc(ii)) then
                     rcross=rr(ii)
                     exit
                  end if
               end if
            end do

            do ii=1,mmax
               yy=rr(ii)/xt(2)
               call gg1cc(gg,yy)
               rhomod(ii,1)= xt(1)*gg
            end do

            call der2exc(rhotps,rhomod(1,1),rhops,rr,d2excps,d2excae,d2mdiff, &
            &                     zion,iexc,nc,nv,la,irmod,mmax)
            fta(jj)=1.0d3*d2mdiff

            if(d2mdiff<d2min) then
               xx(:,1)=xt(:)
               d2min=d2mdiff
            end if

         end do
         write(6,'(f5.1,f9.3,9f7.3)') 1.0d0+0.1d0*(kk-1),(fta(jj),jj=1,10)
      end do

!Initial Nelder-Mead simplex from coarse-search minimum

      xx(1,2)=xx(1,1)+0.25d0*rhocmatch
      xx(2,2)=xx(2,1)
      xx(1,3)=xx(1,1)
      xx(2,3)=xx(2,1)+0.05d0*rmatch
      xx(1,1)=xx(1,1)-0.125d0*rhocmatch
      xx(2,1)=xx(2,1)-0.025d0*rmatch

!Fill function values for initial simplex

      do kk=1,3
         xt(:)=xx(:,kk)


         r0=1.5d0*xt(2)
         do ii=mmax,1,-1
            if(rr(ii)<r0) then
               call gg1cc(gg,yy)
               if(xt(1)*gg<rhoc(ii)) then
                  rcross=rr(ii)
                  exit
               end if
            end if
         end do

         do ii=1,mmax
            yy=rr(ii)/xt(2)
            call gg1cc(gg,yy)
            tt=(rr(ii)-r0-2.0d0*rcross)/(r0-rcross)
            rhomod(ii,1)= xt(1)*gg
         end do

         call der2exc(rhotps,rhomod(1,1),rhops,rr,d2excps,d2excae,d2mdiff, &
         &                    zion,iexc,nc,nv,la,irmod,mmax)
         ff(kk)=d2mdiff

!  write(6,'(i4,a,2f10.4,1p,e14.4)') kk,'  xt',xt(1),xt(2),ff(kk)
      end do

!Nelder-Mead iteration loop
      write(6,'(/a)') 'Nelder-Mead iteration'
      do kk=1,101

!(1) Order (dumb bubble sort)
         do jj=1,3
            if(ff(3)<ff(2)) then
               ft=ff(2) ; xt(:)=xx(:,2)
               ff(2)=ff(3) ; xx(:,2)=xx(:,3)
               ff(3)=ft ; xx(:,3)=xt(:)
            end if
            if(ff(2)<ff(1)) then
               ft=ff(1) ; xt(:)=xx(:,1)
               ff(1)=ff(2) ; xx(:,1)=xx(:,2)
               ff(2)=ft ; xx(:,2)=xt(:)
            end if
         end do

!stopping criterion
         if(ff(3)-ff(1)<1.0d-4*d2ref) then
            write(6,'(a,i4,a)') ' converged in',kk-1,' steps'
            write(6,'(a)') 'amplitude prefactor, scale prefactor'
            write(6,'(2f10.4)') xx(1,1)/rhocmatch,xx(2,1)/rmatch
            exit
         else if(kk==101) then
            write(6,'(a,i4,a)') ' WARNING: not fully converged in 100 steps'
            write(6,'(a)') 'amplitude prefactor, scale prefactor'
            write(6,'(2f10.4)') xx(1,1)/rhocmatch,xx(2,1)/rmatch
            exit
         end if

!(2) Centroid
         x0(:)=0.5d0*(xx(:,1)+xx(:,2))

!(3) Reflection
         xr(:)=x0(:)+alpha_nm*(x0(:)-xx(:,3))


         r0=1.5d0*xr(2)
         do ii=mmax,1,-1
            if(rr(ii)<r0) then
               call gg1cc(gg,yy)
               if(xr(1)*gg<rhoc(ii)) then
                  rcross=rr(ii)
                  exit
               end if
            end if
         end do

         do ii=1,mmax
            yy=rr(ii)/xr(2)
            call gg1cc(gg,yy)
            tt=(rr(ii)-r0-2.0d0*rcross)/(r0-rcross)
            rhomod(ii,1)= xr(1)*gg
         end do

         call der2exc(rhotps,rhomod(1,1),rhops,rr,d2excps,d2excae,d2mdiff, &
         &                    zion,iexc,nc,nv,la,irmod,mmax)
         fr=d2mdiff
! write(6,'(i4,a,2f10.4,1p,e14.4)') kk,'  xr',xr(1),xr(2),fr

         if(ff(1)<= fr .and. fr<ff(2)) then
            ff(3)=fr ; xx(:,3)=xr(:)
            cycle  !kk loop
         end if

!(4) Expansion
         if(fr<ff(1)) then
            xe(:)=x0(:)+gamma_nm*(x0(:)-xx(:,3))

            r0=1.5d0*xe(2)
            do ii=mmax,1,-1
               if(rr(ii)<r0) then
                  call gg1cc(gg,yy)
                  if(xe(1)*gg<rhoc(ii)) then
                     rcross=rr(ii)
                     exit
                  end if
               end if
            end do

            do ii=1,mmax
               yy=rr(ii)/xe(2)
               call gg1cc(gg,yy)
               tt=(rr(ii)-r0-2.0d0*rcross)/(r0-rcross)
               rhomod(ii,1)= xe(1)*gg
            end do

            call der2exc(rhotps,rhomod(1,1),rhops,rr,d2excps,d2excae,d2mdiff, &
            &                     zion,iexc,nc,nv,la,irmod,mmax)
            fe=d2mdiff
! write(6,'(i4,a,2f10.4,1p,e14.4)') kk,'  xe',xe(1),xe(2),fe

            if(fe<fr) then
               ff(3)=fe ; xx(:,3)=xe(:)
               cycle  !kk
            else
               ff(3)=fr ; xx(:,3)=xr(:)
               cycle  !kk
            end if
         end if

!(5) Contraction
         xc(:)=x0(:)+rho_nm*(x0(:)-xx(:,3))


         r0=1.5d0*xc(2)
         do ii=mmax,1,-1
            if(rr(ii)<r0) then
               call gg1cc(gg,yy)
               if(xc(1)*gg<rhoc(ii)) then
                  rcross=rr(ii)
                  exit
               end if
            end if
         end do

         do ii=1,mmax
            yy=rr(ii)/xc(2)
            call gg1cc(gg,yy)
            tt=(rr(ii)-r0-2.0d0*rcross)/(r0-rcross)
            rhomod(ii,1)= xc(1)*gg
         end do

         call der2exc(rhotps,rhomod(1,1),rhops,rr,d2excps,d2excae,d2mdiff, &
         &                    zion,iexc,nc,nv,la,irmod,mmax)
         fc=d2mdiff
!  write(6,'(i4,a,2f10.4,1p,e14.4)') kk,'  xc',xc(1),xc(2),fc
         if(fc<ff(3)) then
            ff(3)=fc ; xx(:,3)=xc(:)
            cycle  !kk
         end if

!(6) Reduction
         do jj=2,3
            xx(:,jj)=xx(:,1)+sigma_nm*(xx(:,jj)-xx(:,1))

            r0=1.5d0*xx(2,jj)
            do ii=mmax,1,-1
               if(rr(ii)<r0) then
                  call gg1cc(gg,yy)
                  if(xx(1,jj)*gg<rhoc(ii)) then
                     rcross=rr(ii)
                     exit
                  end if
               end if
            end do

            do ii=1,mmax
               yy=rr(ii)/xx(2,jj)
               call gg1cc(gg,yy)
               tt=(rr(ii)-r0-2.0d0*rcross)/(r0-rcross)
               rhomod(ii,1)= xx(1,jj)*gg
            end do

            call der2exc(rhotps,rhomod(1,1),rhops,rr,d2excps,d2excae,d2mdiff, &
            &                     zion,iexc,nc,nv,la,irmod,mmax)
            ff(jj)=d2mdiff
!  write(6,'(i4,a,2f10.4,1p,e14.4)') kk,' xrd',xx(1,jj),xx(2,jj),ff(jj)
         end do  !jj

      end do  !kk

      write(6,'(/a)') 'Optimized Teter model core charge'

!option for specifying prefactors as input
   else if(icmod==3) then
      xx(1,1)=fcfact*rhocmatch
      xx(2,1)=rcfact*rmatch

      write(6,'(/a/)') &
      &       'Teter function model core charge with specified parameters'
   end if

   r0=1.5d0*xx(2,1)
   do ii=mmax,1,-1
      if(rr(ii)<r0) then
         call gg1cc(gg,yy)
         if( xx(1,1)*gg<rhoc(ii)) then
            rcross=rr(ii)
            ircross=ii
            exit
         end if
      end if
   end do

! blend the Teter function tail into the all-electron rhoc
! first two derivatives are filled in analytically before the blend starts
   rhomod(:,:)=0.0d0
   do ii=1,mmax
      yy=rr(ii)/xx(2,1)
      call gg1cc(gg,yy)
      tt=(rr(ii)-r0-2.0d0*rr(ircross))/(r0-rr(ircross))

      rhomod(ii,1)=xx(1,1)*gg
      if(ii<ircross) then
         call gp1cc(gg,yy)
         rhomod(ii,2)=xx(1,1)*gg/xx(2,1)
         call gpp1cc(gg,yy)
         rhomod(ii,3)=xx(1,1)*gg/xx(2,1)**2
      end if
   end do

   call der2exc(rhotps,rhomod(1,1),rhops,rr,d2excps,d2excae,d2mdiff, &
   &                   zion,iexc,nc,nv,la,irmod,mmax)

   write(6,'(/a/)') 'd2excps - pseudofunction derivatives with core correction'
   do kk=1,nv
      write(6,'(1p,4d16.6)') (d2excps(kk,jj),jj=1,nv)
   end do
   write(6,'(/a,1p,e16.6)') 'rms 2nd derivative error',d2mdiff

!7-point numerical first derivatives applied successively
!skip non-blended section for 1st and 2nd derivatives

   al = 0.01d0 * dlog(rr(101) / rr(1))

   do jj=2,3
      do ii=ircross-6,mmax-3
         rhomod(ii,jj)=(-rhomod(ii-3,jj-1)+ 9.d0*rhomod(ii-2,jj-1)&
         &     -45.d0*rhomod(ii-1,jj-1)+45.d0*rhomod(ii+1,jj-1)&
         &     -9.d0*rhomod(ii+2,jj-1)+rhomod(ii+3,jj-1))&
         &     /(60.d0*al*rr(ii))
      end do
   end do

!ad-hoc treatment of numerical noise near origin
!set up a mesh on which 2nd derivative will have a stable
!polynomial representation
!assumes dpnint remains 7th order
   drint=0.02d0*rr(ircross)
   rtst=0.5d0*drint
   do jj=1,4
      do ii=1,ircross
         if(rr(ii)>rtst) then
            rint(4+jj)=rr(ii)
            rint(5-jj)=-rr(ii)
            fint(4+jj)=rhomod(ii,3)
            fint(5-jj)=rhomod(ii,3)
            iint=ii
            rtst=rr(ii)+drint
            exit
         end if
      end do
   end do

   call dpnint(rint,fint,8,rr,rhomod(1,3),iint-1)

   jj=4
   do ii=3*jj-8,mmax-3
      rhomod(ii,jj)=(-rhomod(ii-3,jj-1)+ 9.d0*rhomod(ii-2,jj-1)&
      &     -45.d0*rhomod(ii-1,jj-1)+45.d0*rhomod(ii+1,jj-1)&
      &     -9.d0*rhomod(ii+2,jj-1)+rhomod(ii+3,jj-1))&
      &     /(60.d0*al*rr(ii))
   end do

!set up a mesh on which 3rd derivative will have a stable
   drint=0.95d0*drint
   rtst=0.5d0*drint
   rtst=0.5d0*drint
   do jj=1,4
      do ii=1,ircross
         if(rr(ii)>rtst) then
            rint(4+jj)=rr(ii)
            rint(5-jj)=-rr(ii)
            fint(4+jj)=rhomod(ii,4)
            fint(5-jj)=-rhomod(ii,4)
            iint=ii
            rtst=rr(ii)+drint
            exit
         end if
      end do
   end do

   call dpnint(rint,fint,8,rr,rhomod(1,4),iint-1)

   jj=5
   do ii=3*jj-8,mmax-3
      rhomod(ii,jj)=(-rhomod(ii-3,jj-1)+ 9.d0*rhomod(ii-2,jj-1)&
      &     -45.d0*rhomod(ii-1,jj-1)+45.d0*rhomod(ii+1,jj-1)&
      &     -9.d0*rhomod(ii+2,jj-1)+rhomod(ii+3,jj-1))&
      &     /(60.d0*al*rr(ii))
   end do

!set up a mesh on which 4th derivative will have a stable
   drint=0.95d0*drint
   rtst=0.5d0*drint
   do jj=1,4
      do ii=1,ircross
         if(rr(ii)>rtst) then
            rint(4+jj)=rr(ii)
            rint(5-jj)=-rr(ii)
            fint(4+jj)=rhomod(ii,5)
            fint(5-jj)=rhomod(ii,5)
            iint=ii
            rtst=rr(ii)+drint
            exit
         end if
      end do
   end do

   call dpnint(rint,fint,8,rr,rhomod(1,5),iint-1)

   deallocate(vxcae,vxcpsp,vo)
   deallocate(dvxcae,dvxcps)
   deallocate(d2excae,d2excps)

   return
end subroutine modcore3
