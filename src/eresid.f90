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
! calculates residual energy operator matrix elements as a function of q

subroutine eresid(ll,irc,nnull,nbas,mmax,rr,qmax,qroot, &
&                  uu,pswf0_sb,pswfnull_sb, nqout,qout, &
&                  eresid0,eresiddot,eresidmat)

!ll  anglar momentum
!irc  log mesh index of rc
!nnull  number of unconstrained basis vectors for residual minimization
!nbas  number of sbf basis functions
!mmax  index of largest radius for non-zero value of wave function
!rr(mmax)  log radial mesh
!dr  linear radial mesh spacing for Fourier transforms
!dq  linear wave vector mesh spacing for residual kinetic energy integration
!qmax  largest wave vector for E_res integration
!qroot(nbas)  q values for sbfs in basis
!uu(mmax)  r*all-electron wave function
!pswf0_sb(nbas)  sbf basis  coefficients for constrained part of pseudo
!                    wave function
!pswfnull_sb(nbas,nnull) sbf coefficients for unconstrained basis of pseudo
!                       wave function
!nqout  number of wave vector lower cutoffs
!qout(nqout)  wave vector lower cutoffs for E_res for which output is wanted,
!       rounded on output to nearest multiple of dq
!nqout  number of q values in [0,qmax] for which results are to be saved
!qout(nqout)  set of these q values
!eresid0(nqout)  set of <pswf0| E_resid |pswf0> matrix elements
!eresiddot(nnull,nqout) set of <pswfnull| E_resid |pswf0> matrix elements
!eresidevec(nnull,null,nqout)  set of <pswfnull| E_res |pswfnull'> matrices

   implicit none
   integer, parameter :: dp=kind(1.0d0)
   real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp

!Input variables
   integer :: irc,ll,nnull,nbas,mmax,nqout
   real(dp) :: rr(mmax),uu(mmax),pswf0_sb(nbas),pswfnull_sb(nbas,nnull)
   real(dp) :: qroot(nbas)
!real(dp) :: dr,dq,qmax
   real(dp) :: qmax

!Output variables
   real(dp) :: eresid0(nqout),eresiddot(nnull,nqout)
   real(dp) :: eresidmat(nnull,nnull,nqout)

!Input/output variables
   real(dp) :: qout(nqout)

!Local variables
   real(dp), allocatable :: rlin(:),uulin(:),sbf(:),ps0(:)
   real(dp), allocatable :: psnull(:,:),tpsnull(:),psnullq(:)
   real(dp), allocatable :: dac0(:),dacdot(:,:),dacmat(:,:,:)
   real(dp), allocatable :: acdot(:),acmat(:,:)
   real(dp) :: sb_out(5)
   real(dp) :: ac0,rc,ps0q,tps0,xx,qq,qq4
   real(dp) :: dr,dq

   integer :: nq,irclin,mmaxlin
   integer :: ii,iq,jj,ll1


! new calculations of dq, dr
! larger denominators yield better convergence at the cost of more time
! dq=2.0d0*pi/(300.0d0*rr(irc))
   dq=2.0d0*pi/(150.0d0*rr(irc))
! dq=2.0d0*pi/(75.0d0*rr(irc))
! dr=2.0d0*pi/(150.0d0*qmax)
   dr=2.0d0*pi/(75.0d0*qmax)

!set up linear radial mesh
   rc=rr(irc)
   irclin=int(rc/dr + 1.5d0)
   dr=rc/(irclin-1)

   do ii=mmax,1,-1
      if(abs(uu(ii))>1.0d-10) then
         mmaxlin=int(rr(ii)/dr - 2.0d0)
         exit
      end if
   end do
!mmaxlin=rr(mmax)/dr - 1.5d0 !leave non-zero room for cubic interpolation

   allocate(rlin(mmaxlin),uulin(mmaxlin),ps0(mmaxlin))
   allocate(psnull(irclin,nnull),tpsnull(nnull),psnullq(nnull))
   allocate(sbf(mmaxlin))
   allocate(dac0(4),dacdot(nnull,4),dacmat(nnull,nnull,4))
   allocate(acdot(nnull),acmat(nnull,nnull))


   do ii=1,mmaxlin
      rlin(ii)=(ii-1)*dr
   end do

!interpolate all-electron wave function onto linear mesh, even though we
!will only use the tail region

   call dpnint(rr,uu,mmax,rlin,uulin,mmaxlin)

   do ii=irclin,mmaxlin
      ps0(ii)=uulin(ii)
   end do

!set up constrained basis functions on linear grid inside rc (multiplied by r)
   ll1=ll+1
   do ii=1,irclin
      tps0=0.0d0
      tpsnull(:)=0.0d0
      do jj=1,nbas
         xx=rlin(ii)*qroot(jj)
         call sbf8(ll1,xx,sb_out)
         tps0      =tps0      +pswf0_sb(jj)     *sb_out(ll1)
         tpsnull(:)=tpsnull(:)+pswfnull_sb(jj,:)*sb_out(ll1)
      end do
      ps0(ii)=rlin(ii)*tps0
      psnull(ii,:)=rlin(ii)*tpsnull(:)
   end do

!big q integration loop

!q integrals are performed on all matrix elements simultaneously starting
!inward from qmax (considered infinity) to qq=0, saving snapshots along the
!way

   nq=int(qmax/dq)+1

!null accumulators and running registers

   ac0=0.0d0
   acdot(:)=0.0d0
   acmat(:,:)=0.0d0

   dac0(:)=0.0d0
   dacdot(:,:)=0.0d0
   dacmat(:,:,:)=0.0d0

   do iq=nq,1,-1
      qq=dq*(iq-1)

!sbfs on full rlin mesh for this qq
!multiply by rlin for integration
      do ii=1,mmaxlin
         xx=qq*rlin(ii)
         call sbf8(ll1,xx,sb_out)
         sbf(ii)=rlin(ii)*sb_out(ll1)
      end do

!Fourier transforms of basis functions
!null functions exist for r<rc

      ps0q=( 9.0d0*ps0(1       )*sbf(1       ) &
      &      +28.0d0*ps0(2       )*sbf(2       ) &
      &      +23.0d0*ps0(3       )*sbf(3    )    &
      &      +23.0d0*ps0(irclin-2)*sbf(irclin-2) &
      &      +28.0d0*ps0(irclin-1)*sbf(irclin-1) &
      &      + 9.0d0*ps0(irclin  )*sbf(irclin  ))/24.0d0

      psnullq(:)=( 9.0d0*psnull(1       ,:)*sbf(1       ) &
      &            +28.0d0*psnull(2       ,:)*sbf(2       ) &
      &            +23.0d0*psnull(3       ,:)*sbf(3    )    &
      &            +23.0d0*psnull(irclin-2,:)*sbf(irclin-2) &
      &            +28.0d0*psnull(irclin-1,:)*sbf(irclin-1) &
      &            + 9.0d0*psnull(irclin  ,:)*sbf(irclin  ))/24.0d0

      do ii=4,irclin-3
         ps0q=ps0q+ps0(ii)*sbf(ii)
         psnullq(:)=psnullq(:)+psnull(ii,:)*sbf(ii)
      end do

      ps0q=ps0q +( 9.0d0*ps0(irclin   )*sbf(irclin   ) &
      &            +28.0d0*ps0(irclin+1 )*sbf(irclin+1 ) &
      &            +23.0d0*ps0(irclin+2 )*sbf(irclin+2 )    &
      &            +23.0d0*ps0(mmaxlin-2)*sbf(mmaxlin-2) &
      &            +28.0d0*ps0(mmaxlin-1)*sbf(mmaxlin-1) &
      &            + 9.0d0*ps0(mmaxlin  )*sbf(mmaxlin  ))/24.0d0

      do ii=irclin+3,mmaxlin-3
         ps0q=ps0q+ps0(ii)*sbf(ii)
      end do

      ps0q=dr*ps0q
      psnullq(:)=dr*psnullq(:)

!stuff integrands for E_resid q integration
      qq4=dq*qq**4
      dac0(1)=qq4*ps0q**2
      dacdot(:,1)=qq4*ps0q*psnullq(:)
      do jj=1,nnull
         dacmat(:,jj,1)=qq4*psnullq(:)*psnullq(jj)
      end do

!inward integration as soon as dac arrays ("derivative accumulator arrays")
!are full
      if(iq<nq-3) then
         ac0=ac0+(9.0d0*dac0(1)+19.0d0*dac0(2) &
         &      -5.0d0*dac0(3)+dac0(4))/24.0d0

         acdot(:)=acdot(:)+(9.0d0*dacdot(:,1)+19.0d0*dacdot(:,2) &
         &      -5.0d0*dacdot(:,3)+dacdot(:,4))/24.0d0

         acmat(:,:)=acmat(:,:)+(9.0d0*dacmat(:,:,1)+19.0d0*dacmat(:,:,2) &
         &      -5.0d0*dacmat(:,:,3)+dacmat(:,:,4))/24.0d0
      end if

!shift dac* arrays one step inward

      do ii=4,2,-1
         dac0(ii)=dac0(ii-1)
         dacdot(:,ii)=dacdot(:,ii-1)
         dacmat(:,:,ii)=dacmat(:,:,ii-1)
      end do
      dac0(1)=0.0d0
      dacdot(:,1)=0.0d0
      dacmat(:,:,1)=0.0d0

!test if present qq value is a "save" point
!1/pi factor from implicit sbf*sbf integral
      do ii=1,nqout
         if(abs(qout(ii)-qq)<0.5d0*dq) then
            eresid0(ii)=ac0/pi
            eresiddot(:,ii)=acdot(:)/pi
            eresidmat(:,:,ii)=acmat(:,:)/pi
            qout(ii)=qq
         end if
      end do

   end do  !big iq loop


   deallocate(rlin,uulin,ps0)
   deallocate(psnull,tpsnull,psnullq)
   deallocate(sbf)
   deallocate(dac0,dacdot,dacmat)
   deallocate(acdot,acmat)

   return
end subroutine eresid
