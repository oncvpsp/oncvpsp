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
 program oncvpsp
!
! Creates and tests optimized norm-conserving Vanderbilt or Kleinman-Bylander
! pseudopotentials based on D. R. Hamann, Phys. Rev. B 88, 085117 (2013)
! and references therein.
!
!   D. R. Hamann
!   Mat-Sim Research LLC
!   P.O. Box 742
!   Murray Hill, NJ 07974
!   USA
!
!   Developed from original "gncpp" code of March 8,1987
!
!   Output format for ABINIT pspcod=8 and upf format for quantumespresso
!
 implicit none
 integer, parameter :: dp=kind(1.0d0)

!
 integer :: ii,ierr,iexc,iexct,ios,iprint,irps,it,icmod,lpopt
 integer :: jj,kk,ll,l1,lloc,lmax,lt,inline
 integer :: mch,mchf,mmax,n1,n2,nc,nlim,nlloc,nlmax,irpsh,nrl
 integer :: nv,irct,ncnf,nvt
 integer :: iprj,mxprj
 integer,allocatable :: ae_bound_well_n_qn(:,:)
!
 integer :: dtime(8),na(30),la(30),np(6)
 integer :: nacnf(30,5),lacnf(30,5),nvcnf(5)
 integer :: nat(30),lat(30),indxr(30),indxe(30)
 integer :: irc(6),nodes(4)
 integer :: nproj(6),npx(6),lpx(6)
 integer :: ncon(6),nbas(6)

 real(dp) :: al,amesh,csc,csc1,deblt,depsh,depsht,drl,eeel
 real(dp) :: eeig,eexc
 real(dp) :: emax,epsh1,epsh1t,epsh2,epsh2t
 real(dp) :: et,etest,emin,sls
 real(dp) :: fcfact,rcfact,dvloc0
 real(dp) :: rr1,rcion,rciont,rcmax,rct,rlmax,rpkt
 real(dp) :: sf,zz,zion,zval,etot
 real(dp) :: xdummy
!
 real(dp) :: cl(6),debl(6),ea(30),ep(6),fa(30),facnf(30,5)
 real(dp) :: fat(30,2)
 real(dp) :: fnp(6),fp(6)
 real(dp) :: qcut(6),qmsbf(6),rc(6),rc0(6)
 real(dp) :: rpk(30)
 real(dp) :: epx(6),fpx(6)
 real(dp) :: epstot
 real(dp), parameter :: eps=1.0d-8

 real(dp), allocatable :: vkb_coef(:,:),cvgplt(:,:,:,:),ae_bound_well_overlap(:,:)
 real(dp), allocatable :: rr(:)
 real(dp), allocatable :: ps_rho_val(:),rhoc(:),rhot(:)
 real(dp), allocatable :: uu(:),up(:)
 real(dp), allocatable :: v_ps_sl(:,:),vfull(:),vkb_proj(:,:,:),ps_rpsi(:,:,:)
 real(dp), allocatable :: vwell(:)
 real(dp), allocatable :: vpuns(:,:)
 real(dp), allocatable :: vo(:),vxc(:)
 real(dp), allocatable :: rhomod(:,:),ae_psi2_val(:,:),ps_psi2_val(:,:),ae_rho_val(:)
 real(dp), allocatable :: uupsa(:,:) !pseudo-atomic orbitals array
 real(dp), allocatable :: ae_bound_well_eig(:,:),fpa(:,:)
 real(dp), allocatable :: ae_bound_well_rpsi(:,:),ae_bound_well_drpsi_dr(:,:)
 real(dp), allocatable :: vr(:,:,:)

 character*2 :: atsym
 character*4 :: psfile

 logical :: srel,cset

 write(6,'(a/a//)') &
&      'ONCVPSP  (Optimized Norm-Conservinng Vanderbilt PSeudopotential)', &
&      'scalar-relativistic version 4.0.1 03/01/2019'

 write(6,'(a/a/a//)') &
&      'While it is not required under the terms of the GNU GPL, it is',&
&      'suggested that you cite D. R. Hamann, Phys. Rev. B 88, 085117 (2013)', &
&      'in any publication utilizing these pseudopotentials.'

 srel=.true.
!srel=.false.

 nproj(:)=0
 fcfact=0.0d0
 rcfact=0.0d0
 rc(:)=0.0d0
 ep(:)=0.0d0

! read input data
 inline=0

! atom and reference configuration
 call cmtskp(inline)
 read(5,*,iostat=ios) atsym,zz,nc,nv,iexc,psfile
 call read_error(ios,inline)

 call cmtskp(inline)
 do ii=1,nc+nv
   read(5,*,iostat=ios) na(ii),la(ii),fa(ii)
   call read_error(ios,inline)
 end do

! pseudopotential and optimization
 call cmtskp(inline)
 read(5,*,iostat=ios) lmax
 call read_error(ios,inline)

 call cmtskp(inline)
 do l1=1,lmax+1
   read(5,*,iostat=ios) lt,rc(l1),ep(l1),ncon(l1),nbas(l1),qcut(l1)
   if(lt/=l1-1) ios=999
   call read_error(ios,inline)
 end do

! local potential
 call cmtskp(inline)
 read(5,*, iostat=ios) lloc,lpopt,rc(5),dvloc0
 call read_error(ios,inline)

! Vanderbilt-Kleinman-Bylander projectors
 call cmtskp(inline)
 do l1=1,lmax+1
   read(5,*,iostat=ios) lt,nproj(l1),debl(l1)
   if(lt/=l1-1) ios=999
   call read_error(ios,inline)
 end do

! model core charge
 call cmtskp(inline)
 read(5,*, iostat=ios) icmod, fcfact
 if(ios==0 .and. icmod==3) then
  backspace(5)
  read(5,*, iostat=ios) icmod, fcfact, rcfact
 end if
 call read_error(ios,inline)

! log derivative analysis
 call cmtskp(inline)
 read(5,*,iostat=ios) epsh1, epsh2, depsh
 call read_error(ios,inline)

! output grid
 call cmtskp(inline)
 read(5,*,iostat=ios) rlmax,drl
 call read_error(ios,inline)

! test configurations
 call cmtskp(inline)
 read(5,*,iostat=ios) ncnf
 call read_error(ios,inline)

 do jj=2,ncnf+1
   call cmtskp(inline)
   read(5,*,iostat=ios) nvcnf(jj)
   call read_error(ios,inline)
   do ii=nc+1,nc+nvcnf(jj)
     call cmtskp(inline)
     read(5,*,iostat=ios) nacnf(ii,jj),lacnf(ii,jj),facnf(ii,jj)
     call read_error(ios,inline)
   end do
 end do

! end of reading input data

 nvcnf(1)=nv
 do ii=1,nc+nv
   nacnf(ii,1)=na(ii)
   lacnf(ii,1)=la(ii)
   facnf(ii,1)=fa(ii)
 end do

 do jj=2,ncnf+1
   do ii=1,nc
     nacnf(ii,jj)=na(ii)
     lacnf(ii,jj)=la(ii)
     facnf(ii,jj)=fa(ii)
   end do
 end do

 if(icmod==0) then
   fcfact=0.0d0
 end if

 call check_data(atsym,zz,fcfact,rcfact,epsh1,epsh2,depsh,rlmax,drl,fa,facnf, &
&                rc,ep,qcut,debl,nc,nv,iexc,lmax,lloc,lpopt,icmod, &
&                ncnf,na,la,nvcnf,nacnf,lacnf,ncon,nbas,nproj,psfile)

 nrl=int((rlmax/drl)-0.5d0)+1

!PWSCF wants an even number of mesh pointe
!if(trim(psfile)=='upf') then
  if(mod(nrl,2)/=0) nrl=nrl+1
!end if

!amesh=1.012d0
 amesh=1.006d0
!amesh=1.003d0
!amesh=1.0015d0

 mxprj=5

 al=dlog(amesh)
 rr1=0.0005d0/zz
 rr1=dmin1(rr1,0.0005d0/10)
 mmax=dlog(45.0d0 /rr1)/al

!calculate zion for output
 zion=zz
 do ii=1,nc
  zion=zion-fa(ii)
 end do


 allocate(rr(mmax))
 allocate(ps_rho_val(mmax),rhoc(mmax),rhot(mmax))
 allocate(uu(mmax),up(mmax),uupsa(mmax,30))
 allocate(vkb_coef(mxprj,4), cvgplt(2,7,mxprj,4),ae_bound_well_overlap(mxprj,mxprj))
 allocate(v_ps_sl(mmax,5),vfull(mmax),vkb_proj(mmax,mxprj,4),ps_rpsi(mmax,mxprj,4))
 allocate(vwell(mmax))
 allocate(vpuns(mmax,5))
 allocate(vo(mmax),vxc(mmax))
 allocate(ae_psi2_val(mmax,nv),ps_psi2_val(mmax,nv),ae_rho_val(mmax))
 allocate(ae_bound_well_n_qn(mxprj,6))
 allocate(ae_bound_well_eig(mxprj,6),fpa(mxprj,6))
 allocate(ae_bound_well_rpsi(mmax,mxprj),ae_bound_well_drpsi_dr(mmax,mxprj))
 allocate(vr(mmax,mxprj,6))

 vr(:,:,:)=0.0d0
 v_ps_sl(:,:)=0.0d0
 vkb_proj(:,:,:)=0.0d0
 ae_bound_well_eig(:,:)=0.0d0

 do ii=1,mmax
  rr(ii)=rr1*exp(al*(ii-1))
 end do

!
! full potential atom solution
!
   call sratom(na,la,ea,fa,rpk,nc,nc+nv,it,rhoc,ps_rho_val, &
   &           rr,vfull,vxc,zz,mmax,iexc,etot,ierr,srel)
!
!

! Drop digits beyond 5 decimals for input rcs before making any use of them
 do l1=1,max(lmax+1,lloc+1)
    jj=int(rc(l1)*10.0d5)
    rc(l1)=jj/10.0d5
 end do

 rcmax=0.0d0
 do l1=1,lmax+1
   rcmax=dmax1(rcmax,rc(l1))
 end do
 do l1=1,lmax+1
   if(rc(l1)==0.0d0) then
     rc(l1)=rcmax
   end if
 end do

 nproj(lloc+1)=0
 rc0(:)=rc(:)

! output printing (echos input data, with all-electron eigenvalues added)

 write(6,'(a)') '# ATOM AND REFERENCE CONFIGURATION'
 write(6,'(a)') '# atsym  z   nc   nv     iexc    psfile'
 write(6,'(a,a,f6.2,2i5,i8,2a)') '  ',trim(atsym),zz,nc,nv,iexc, &
&      '      ',psfile
 write(6,'(a/a)') '#','#   n    l    f        energy (Ha)'
 do ii=1,nc+nv
   write(6,'(2i5,f8.2,1pe18.7)') na(ii),la(ii),fa(ii),ea(ii)
 end do

 write(6,'(a/a/a)') '#','# PSEUDOPOTENTIAL AND OPTIMIZATION','# lmax'
 write(6,'(i5)')  lmax
 write(6,'(a/a)') '#','#   l,   rc,      ep,       ncon, nbas, qcut'
 do l1=1,lmax+1
   write(6,'(i5,2f10.5,2i5,f10.5)') l1-1,rc(l1),ep(l1),ncon(l1),&
&        nbas(l1),qcut(l1)
 end do

 write(6,'(a/a/a,a)') '#','# LOCAL POTENTIAL','# lloc, lpopt,  rc(5),', &
&      '   dvloc0'
 write(6,'(2i5,f10.5,a,f10.5)') lloc,lpopt,rc(5),'   ',dvloc0

 write(6,'(a/a/a)') '#','# VANDERBILT-KLEINMAN-BYLANDER PROJECTORs', &
&      '# l, nproj, debl'
 do l1=1,lmax+1
   write(6,'(2i5,f10.5)') l1-1,nproj(l1),debl(l1)
 end do

 write(6,'(a/a/a)') '#','# MODEL CORE CHARGE', &
&      '# icmod, fcfact, rcfact'
 write(6,'(i5,2f10.5)') icmod,fcfact,rcfact

 write(6,'(a/a/a)') '#','# LOG DERIVATIVE ANALYSIS', &
&      '# epsh1, epsh2, depsh'
 write(6,'(3f8.2)') epsh1,epsh2,depsh

 write(6,'(a/a/a)') '#','# OUTPUT GRID','# rlmax, drl'
 write(6,'(2f8.4)') rlmax,drl

 write(6,'(a/a/a)') '#','# TEST CONFIGURATIONS','# ncnf'
 write(6,'(i5)') ncnf
 write(6,'(a/a)') '# nvcnf','#   n    l    f'
 do jj=2,ncnf+1
   write(6,'(i5)') nvcnf(jj)
   do ii=nc+1,nc+nvcnf(jj)
     write(6,'(2i5,f8.2)') nacnf(ii,jj),lacnf(ii,jj),facnf(ii,jj)
   end do
   write(6,'(a)') '#'
 end do

 write(6,'(//a)') 'Reference configufation results'
 write(6,'(a,i6)') '  iterations',it
 if(it .ge. 100) then
  write(6,'(a)') 'oncvpsp: ERROR all-electron reference atom not converged'
  stop
 end if
 write(6,'(a,1p,d18.8)') '  all-electron total energy (Ha)',etot

!find log mesh point nearest input rc
 rcmax=0.0d0
 irc(:)=0
 do l1=1,max(lmax+1,lloc+1)
  rct=rc(l1)
  irc(l1)=0
  do ii=2,mmax
   if(rr(ii)>rct) then
    irc(l1)=ii
    rc(l1)=rr(ii)
    exit
   end if
  end do
  rcmax=dmax1(rcmax,rc(l1))
 end do

!
 cvgplt(:,:,:,:)=0.0d0
!
! loop to construct pseudopotentials for all angular momenta
!
 write(6,'(/a/a)') 'Begin loop to  construct optimized pseudo wave functions',&
&      'and semi-local pseudopoentials for all angular momenta'

!temporarily set this to 1 so that the pseudo wave function needed for the
!local potential will be generated.  Reset after run_vkb.
 nproj(lloc+1)=1
 do l1=1,lmax+1
   ll=l1-1
   uu(:)=0.0d0; ae_bound_well_overlap(:,:)=0.0d0
   iprj=0

!get principal quantum number for the highest core state for this l
   ae_bound_well_n_qn(1,l1)=l1
   do kk=1,nc
     if(la(kk)==l1-1) ae_bound_well_n_qn(1,l1)=na(kk)+1
   end do !kk

!get all-electron bound states for projectors
   if(nv/=0) then
     do kk=nc+1,nc+nv
       if(la(kk)==l1-1) then
         iprj=iprj+1
         et=ea(kk)
         call lschfb(na(kk),la(kk),ierr,et, &
&                      rr,vfull,uu,up,zz,mmax,mch,srel)
         if(ierr /= 0) then
           write(6,'(/a,3i4)') 'oncvpsp-387: lschfb convergence ERROR n,l,iter=', &
&           na(ii),la(ii),it
           stop
         end if
         ae_bound_well_eig(iprj,l1)=ea(kk)
         ae_bound_well_n_qn(iprj,l1)=na(kk)
         ae_bound_well_rpsi(:,iprj)=uu(:)
         ae_bound_well_drpsi_dr(:,iprj)=up(:)
       end if !la(kk)==l1-1
       if(iprj==nproj(l1)) exit
     end do !kk
   end if !nv/=0

!get all-electron well states for projectors
!if there were no valence states, use ep from input data for 1st well state
!otherwise shift up by input debl
   if(iprj==0) ae_bound_well_eig(1,l1)=ep(l1)
   if(iprj<nproj(l1))then
     do kk=1,nproj(l1)-iprj
       iprj=iprj+1
       if(iprj>1 .and. debl(l1)<=0.0d0) then
         write(6,'(a,f8.3,a/a)') 'oncvpsp: ERROR debl =',debl, 'for l=', &
&              ' ERROR not allowed with 2 or more scattering states', &
&              'program will stop'
         stop
       end if
       if(iprj>1) then
         ae_bound_well_eig(iprj,l1)=ae_bound_well_eig(iprj-1,l1)+debl(l1)
         ae_bound_well_n_qn(iprj,l1)=ae_bound_well_n_qn(iprj-1,l1)+1
       end if

       call wellstate(ae_bound_well_n_qn(iprj,l1),ll,irc(l1),ae_bound_well_eig(iprj,l1),rr, &
&                     vfull,uu,up,zz,mmax,mch,srel)
       ae_bound_well_rpsi(:,iprj)=uu(:)
       ae_bound_well_drpsi_dr(:,iprj)=up(:)
     end do !kk
   end if !iprj<nproj(l1)

   do iprj=1,nproj(l1)

!calculate relativistic correction to potential to force projectors to 0 at rc
     call vrel(ll,ae_bound_well_eig(iprj,l1),rr,vfull,vr(1,iprj,l1),ae_bound_well_rpsi(1,iprj),ae_bound_well_drpsi_dr(1,iprj), &
&              zz,mmax,irc(l1),srel)

   end do

!get all-electron overlap matrix
   do jj=1,nproj(l1)
     do ii=1,jj
       call fpovlp(ae_bound_well_rpsi(1,ii),ae_bound_well_rpsi(1,jj),irc(l1),ll,zz,ae_bound_well_overlap(ii,jj),rr,srel)
       ae_bound_well_overlap(jj,ii)=ae_bound_well_overlap(ii,jj)
     end do
   end do

   call run_optimize(ae_bound_well_eig(1,l1),ll,mmax,mxprj,rr,ae_bound_well_rpsi,ae_bound_well_overlap, &
&                    irc(l1),qcut(l1),qmsbf(l1),ncon(l1),nbas(l1),nproj(l1), &
&                    ps_rpsi(1,1,l1),v_ps_sl(1,l1),vkb_proj(1,1,l1),vfull,cvgplt(1,1,1,l1))

 end do !l1

! construct Vanderbilt / Kleinman-Bylander projectors

 write(6,'(/a,a)') 'Construct Vanderbilt / Kleinmman-Bylander projectors'

 call run_vkb(lmax,lloc,lpopt,dvloc0,irc,nproj,rr,mmax,mxprj,ps_rpsi,vfull,v_ps_sl, &
&             vkb_coef,vkb_proj,nlim,vr)

!restore this to its proper value
 nproj(lloc+1)=0

 deallocate(ae_bound_well_rpsi,ae_bound_well_drpsi_dr)

! accumulate charge and eigenvalues
! pseudo wave functions are calculated with VKB projectors for
! maximum consistency of unscreening
! get all-electron and pseudopotential valence-state by valence-state
! charge densities

! null charge and eigenvalue accumulators
 uupsa(:,:)=0.0d0
 eeig=0.0d0
 zval=0.0d0
 ps_rho_val(:)=0.0d0
 nodes(:)=0
 ae_rho_val(:)=0.0d0
 irps=0
 do kk=1,nv

   et=ea(nc+kk)
   ll=la(nc+kk)
   l1=ll+1
   call lschfb(na(nc+kk),ll,ierr,et, &
&                rr,vfull,uu,up,zz,mmax,mch,srel)
   if(ierr /= 0) then
     write(6,'(/a,3i4)') 'oncvpsp-461: lschfb convergence ERROR n,l,iter=', &
&     na(ii),la(ii),it
     stop
   end if

   ae_psi2_val(:,kk)=(uu(:)/rr(:))**2

   ae_rho_val(:)=ae_rho_val(:) + fa(nc+kk)*ae_psi2_val(:,kk)

   emax=0.75d0*et
   emin=1.25d0*et

   call lschvkbb(ll+nodes(l1)+1,ll,nproj(l1),ierr,et,emin,emax, &
&                rr,v_ps_sl(1,lloc+1),vkb_proj(1,1,l1),vkb_coef(1,l1), &
&                uu,up,mmax,mch)

   if(ierr/=0) then
     write(6,'(a,3i4)') 'oncvpsp: lschvkbb ERROR',ll+nodes(l1)+1,ll,ierr
     flush(6)
     stop
   end if

! save valence pseudo wave functions for upfout
   uupsa(:,kk)=uu(:)

   ps_psi2_val(:,kk)=(uu(:)/rr(:))**2
   ps_rho_val(:)=ps_rho_val(:)+fa(nc+kk)*ps_psi2_val(:,kk)
   eeig=eeig+fa(nc+kk)*et

   zval=zval+fa(nc+kk)
   nodes(l1)=nodes(l1)+1
   irps=max(irps,irc(l1))
 end do !kk

 allocate(rhomod(mmax,5))

 rhomod(:,:)=0.0d0

! construct model core charge based on monotonic polynomial fit
! or Teter function fit

 if(icmod==1) then
   call modcore(icmod,ps_psi2_val,ps_rho_val,rhoc,ae_psi2_val,ae_rho_val,rhomod, &
&               fcfact,rcfact,irps,mmax,rr,nc,nv,la,zion,iexc)

 else if(icmod==2) then
   call modcore2(icmod,ps_psi2_val,ps_rho_val,rhoc,ae_psi2_val,ae_rho_val,rhomod, &
&               fcfact,rcfact,irps,mmax,rr,nc,nv,la,zion,iexc)

 else if(icmod>=3) then
   call modcore3(icmod,ps_psi2_val,ps_rho_val,rhoc,ae_psi2_val,ae_rho_val,rhomod, &
&               fcfact,rcfact,irps,mmax,rr,nc,nv,la,zion,iexc)

 end if

! screening potential for pseudocharge

 call vout(1,ps_rho_val,rhomod(1,1),vo,vxc,zval,eeel,eexc,rr,mmax,iexc)

! total energy output

 epstot= eeig + eexc - 0.5d0*eeel
 write(6,'(/a,f12.6/)') 'Pseudoatom total energy', epstot

 call run_diag(lmax,ae_bound_well_n_qn,ae_bound_well_eig,lloc,irc, &
&                    vkb_proj,vkb_coef,nproj,rr,vfull,v_ps_sl,zz,mmax,mxprj,srel)

 call run_ghosts(lmax,la,ea,nc,nv,lloc,irc,qmsbf, &
&                    vkb_proj,vkb_coef,nproj,rr,v_ps_sl,mmax,mxprj)

! unscreen semi-local potentials

 do l1=1,max(lmax+1,lloc+1)
   vpuns(:,l1)=v_ps_sl(:,l1)-vo(:)
 end do

!fix unscreening error due to greater range of all-electron charge
 do ii=mmax,1,-1
   if(ps_rho_val(ii)==0.0d0) then
     do l1=1,max(lmax+1,lloc+1)
       vpuns(ii,l1)=-zion/rr(ii)
     end do
   else
     exit
   end if
 end do


! loop over reference plus test atom configurations
!if(.false.) then

!ncnf=0
 rhot(:)=ps_rho_val(:)
 do jj=1,ncnf+1

   write(6,'(/a,i2)') 'Test configuration',jj-1

! charge density is initialized to that of reference configuration

   rhot(:)=ps_rho_val(:)

   call run_config(jj,nacnf,lacnf,facnf,nc,nvcnf,rhot,rhomod,rr,zz, &
&                  rcmax,mmax,mxprj,iexc,ea,etot,epstot,nproj,vpuns, &
&                  lloc,vkb_proj,vkb_coef,srel)

 end do !jj

 call run_plot(lmax,ae_bound_well_n_qn,ae_bound_well_eig,lloc,irc, &
&                    vkb_proj,vkb_coef,nproj,rr,vfull,v_ps_sl,vpuns,zz,mmax,mxprj,drl,nrl, &
&                    ps_rho_val,rhoc,rhomod,srel,cvgplt)



 call run_phsft(lmax,lloc,nproj,ae_bound_well_eig,epsh1,epsh2,depsh,vkb_proj,vkb_coef, &
&               rr,vfull,v_ps_sl,zz,mmax,mxprj,irc,srel)

 call gnu_script(ae_bound_well_eig,vkb_coef,lmax,lloc,mxprj,nproj)

 if(trim(psfile)=='psp8' .or. trim(psfile)=='both') then


  call linout(lmax,lloc,rc,vkb_proj,vkb_coef,nproj,rr,vpuns,ps_rho_val,rhomod, &
&             ae_rho_val,rhoc,zz,zion,mmax,mxprj,iexc,icmod,nrl,drl,atsym, &
&             na,la,ncon,nbas,nvcnf,nacnf,lacnf,nc,nv,lpopt,ncnf, &
&             fa,rc0,ep,qcut,debl,facnf,dvloc0,fcfact,rcfact, &
&             epsh1,epsh2,depsh,rlmax,psfile)
 end if

 if(trim(psfile)=='upf' .or. trim(psfile)=='both') then
  call upfout(lmax,lloc,rc,vkb_proj,vkb_coef,nproj,rr,vpuns,ps_rho_val,rhomod, &
&             zz,zion,mmax,mxprj,iexc,icmod,nrl,drl,atsym,epstot, &
&             na,la,ncon,nbas,nvcnf,nacnf,lacnf,nc,nv,lpopt,ncnf, &
&             fa,rc0,ep,qcut,debl,facnf,dvloc0,fcfact,rcfact, &
&             epsh1,epsh2,depsh,rlmax,psfile,uupsa,ea)
 end if

 stop
 end program oncvpsp


 subroutine cmtskp(inline)
! skips lines of standard input (file 5) whose first character is #
 implicit none

!In/Out variable
 integer :: inline

!Local variable
 character*1 tst

 tst='#'
 do while (tst=='#')
  read(5,*) tst
  inline=inline+1
 end do
 backspace(5)
 inline=inline-1

 return
 end subroutine cmtskp

 subroutine read_error(ios,inline)
! report data read error and stop
 implicit none

!Input variables
 integer :: ios,inline

 inline=inline+1
 if(ios/=0) then
  write(6,'(a,i4)') 'Read ERROR, input data file line',inline
  write(6,'(a)') 'Program will stop'
  stop
 end if

 return
 end subroutine read_error
