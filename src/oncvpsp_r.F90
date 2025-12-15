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
 program oncvpsp_r
!
! Creates and tests optimized norm-conserving Vanderbilt or Kleinman-Bylander
! pseudopotentials based on D. R. Hamann, Phys. Rev. B 88, 085117 (2013)
! and references therein.
! This version is for the fully-relativistic case
!
!   D. R. Hamann
!   Mat-Sim Research LLC
!   P.O. Box 742
!   Murray Hill, NJ 07974
!   USA
!
!   Developed from original "gncpp" code of March 8,1987
!
!   Output format for ABINIT pspcod=8
!
 implicit none
 integer, parameter :: dp=kind(1.0d0)

!
 integer :: ii,ierr,iexc,iexct,ios,iprint,it,icmod,lpopt
 integer :: jj,ll,l1,lloc,lmax,lt,inline
 integer :: mch,mchf,mmax,n1,n2,nc,nlim,nlloc,nlmax,irpsh,nrl
 integer :: nv,irct,ncnf
 integer :: ikap,kap,mkap
!
 integer :: dtime(8),na(30),la(30),np(6)
 integer :: nacnf(30,5),lacnf(30,5),nvcnf(5)
 integer :: irc(6)
 integer :: ncut,nxtra,nproj(6),npx(6),lpx(6)
 integer :: ncon(6),nbas(6)

 real(dp) :: al,amesh,csc,csc1,deblt,depsh,depsht,drl,eeel,fj
 real(dp) :: eeig,eexc
 real(dp) :: emax,epsh1,epsh1t,epsh2,epsh2t
 real(dp) :: et,etest,emin,sls
 real(dp) :: fcfact,dvloc0,soscale
 real(dp) :: rr1,rcion,rciont,rcmax,rct,rlmax,rpkt
 real(dp) :: sf,vsc,vsc1,zz,zion,zval,etot
 real(dp) :: cnorm,ebar,eprmin
 real(dp), parameter :: eps=1.0d-8
!
 real(dp) :: cl(6),evkb(2,4,2),debl(6,2),ea(30,2),ep(6,2),fa(30),facnf(30,5)
 real(dp) ::fnp(6),fp(6)
 real(dp) :: qcut(6),rc(6),rc0(6)
 real(dp) :: rpk(30,2),cvgplt(2,4,2,4,2)
 real(dp) :: epx(6,2),fpx(6)
 real(dp) :: qq(2,2)
 real(dp) :: esr(4,4),eso(4,4)
 real(dp) :: epstot

 real(dp), allocatable :: rr(:)
 real(dp), allocatable :: rho(:),rhoc(:),rhot(:)
 real(dp), allocatable :: uu(:,:),uu2(:,:),up(:,:),up2(:,:)
 real(dp), allocatable :: vp(:,:,:),vfull(:),vkb(:,:,:,:),pswf(:,:,:,:)
 real(dp), allocatable :: vwell(:)
 real(dp), allocatable :: vpuns(:,:)
 real(dp), allocatable :: vo(:),vxc(:)
 real(dp), allocatable :: rhomod(:,:)
 real(dp), allocatable :: vsr(:,:,:),vso(:,:,:)


 character*2 :: atsym
 character*4 :: psfile

 write(6,'(a/a//)') &
&      'ONCVPSP  (Optimized Norm-Conservinng Vanderbilt PSeudopotential)', &
&      'fully-relativistic version 2.1.1, 03/26/2014'

 write(6,'(a,a,a,a//)') &
&      'While it is not required under the terms of the GNU GPL, it is',&
&      'suggested that you cite the arXiv preprint or Phys. Rev. B paper',&
&      'through which you found this code in any publication using these',&
&      'pseudopotentials.  The reference will be here when available.'

 nproj(:)=0

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
   read(5,*,iostat=ios) lt,rc(l1),ep(l1,1),ncon(l1),nbas(l1),qcut(l1)
   if(lt/=l1-1) ios=999
   call read_error(ios,inline)
   ep(l1,2)=ep(l1,1)
 end do

! local potential
 call cmtskp(inline)
 read(5,*, iostat=ios) lloc,lpopt,rc(5),dvloc0
 call read_error(ios,inline)

! Vanderbilt-Kleinman-Bylander projectors
 call cmtskp(inline)
 do l1=1,lmax+1
   read(5,*,iostat=ios) lt,nproj(l1),debl(l1,1)
   if(lt/=l1-1) ios=999
   call read_error(ios,inline)
   debl(l1,2)=debl(l1,1)
 end do

! model core charge
 call cmtskp(inline)
 read(5,*, iostat=ios) icmod, fcfact
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

 call check_data(atsym,zz,fcfact,epsh1,epsh2,depsh,rlmax,drl,fa,facnf, &
&                rc,ep(1,1),qcut,debl(1,1),nc,nv,iexc,lmax,lloc,lpopt,icmod, &
&                ncnf,na,la,nvcnf,nacnf,lacnf,ncon,nbas,nproj,psfile)

 nrl=int(rlmax/drl)+1

!PWSCF wants an even number of mesh pointe
 if(trim(psfile)=='upf') then
  if(mod(nrl,2)/=0) nrl=nrl+1
 end if


 amesh=1.012d0
!amesh=1.006d0
 al=dlog(amesh)
 rr1=.0005d0/zz
 mmax=dlog(45.0d0 /rr1)/al

!calculate zion for output
 zion=zz
 do ii=1,nc
  zion=zion-fa(ii)
 end do

 allocate(rr(mmax))
 allocate(rho(mmax),rhoc(mmax),rhot(mmax))
 allocate(uu(mmax,2),uu2(mmax,2),up(mmax,2),up2(mmax,2))
 allocate(vp(mmax,5,2),vfull(mmax),vkb(mmax,2,4,2),pswf(mmax,2,4,2))
 allocate(vwell(mmax))
 allocate(vpuns(mmax,5))
 allocate(vo(mmax),vxc(mmax))
 allocate(vsr(mmax,4,4),vso(mmax,4,4))

 do ii=1,mmax
  rr(ii)=rr1*exp(al*(ii-1))
 end do

! full potential atom solution

   call relatom(na,la,ea,fa,rpk,nc,nc+nv,it,rhoc,rho, &
&              rr,vfull,zz,mmax,iexc,etot,ierr)

! output printing (echos input data, with all-electron eigenvalues added)

 write(6,'(a)') '# ATOM AND REFERENCE CONFIGURATION'
 write(6,'(a)') '# atsym  z    nc    nv    iexc   psfile'
 write(6,'(a,a,f6.2,3i6,2a)') '  ',trim(atsym),zz,nc,nv,iexc, &
&      '      ',psfile
 write(6,'(a/a)') '#','#   n    l    f              energy (Ha)'
 write(6,'(a)') '#   n    l    f        l+1/2             l-1/2'
 do ii=1,nc+nv
  if(la(ii)==0) then
   write(6,'(2i5,f8.2,1p,d18.7)') na(ii),la(ii),fa(ii),ea(ii,1)
  else
   write(6,'(2i5,f8.2,1p,2d18.7)') na(ii),la(ii),fa(ii),ea(ii,1),ea(ii,2)
  end if
 end do

 write(6,'(a/a/a)') '#','# PSEUDOPOTENTIAL AND OPTIMIZATION','# lmax'
 write(6,'(i5)')  lmax
 write(6,'(a/a)') '#','#   l,   rc,     ep,   ncon, nbas, qcut'
 do l1=1,lmax+1
   write(6,'(i5,2f8.2,2i5,f8.2)') l1-1,rc(l1),ep(l1,1),ncon(l1),nbas(l1),qcut(l1)
 end do

 write(6,'(a/a/a,a)') '#','# LOCAL POTENTIAL','# lloc, lpopt,  rc(5),', &
&      '   dvloc0'
 write(6,'(2i5,f10.2,a,f8.2)') lloc,lpopt,rc(5),'   ',dvloc0

 write(6,'(a/a/a)') '#','# VANDERBILT-KLEINMAN-BYLANDER PROJECTORs', &
&      '# l, nproj, debl'
 do l1=1,lmax+1
   write(6,'(2i5,f10.4)') l1-1,nproj(l1),debl(l1,1)
 end do

 write(6,'(a/a/a)') '#','# MODEL CORE CHARGE', &
&      '# icmod, fcfact'
 write(6,'(i5,f8.2,2i5,f8.2)') icmod,fcfact

 write(6,'(a/a/a)') '#','# LOG DERIVATIVE ANALYSIS', &
&      '# epsh1, epsh2, depsh'
 write(6,'(3f8.2)') epsh1,epsh2,depsh

 write(6,'(a/a/a)') '#','# OUTPUT GRID','# rlmax, drl'
 write(6,'(2f8.2)') rlmax,drl

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
  write(6,'(a)') 'oncvpsp: all-electron reference atom not converged'
  stop
 end if
 write(6,'(a,1p,d18.8)') '  all-electron total energy (Ha)',etot

! collect pseudopotential information from all-electron configuration
! and eigenvalues
! nxtra counts valence states where semi-core of same l have
! been included as "valence"

 vsc=1.2d0 !scaling factor to get an rc value for "extra" states

 np(:)=0
 fp(:)=0.0d0

 nxtra=0
 npx(:)=0
 lpx(:)=0
 fpx(:)=0.0d0
 epx(:,:)=0.0d0


 if(nv/=0) then
   do ii=nc+1,nc+nv
     l1=la(ii)+1
     if(np(l1)==0) then
       np(l1)=na(ii)
       fp(l1)=fa(ii)
       ep(l1,:)=ea(ii,:)
!test to put 2nd projector scattering states at equal energies
       ebar=(l1*ep(l1,1)+(l1-1)*ep(l1,2))/dble(2*l1-1)
       debl(l1,:)=debl(l1,:)+ebar-ep(l1,:)
     else if (fa(ii) .gt. 0.0d0) then
       nxtra = nxtra + 1
       npx(nxtra)=na(ii)
       lpx(nxtra)=la(ii)
       fpx(nxtra)=fa(ii)
       epx(nxtra,:)=ea(ii,:)
       debl(l1,:)=ea(ii,:)-ep(l1,:)
     end if
   end do
 end if

 rcmax=0.0d0
 do l1=1,lmax+1
   rcmax=dmax1(rcmax,rc(l1))
 end do
 do l1=1,lmax+1
   if(rc(l1)==0.0d0) then
     rc(l1)=rcmax
   end if
 end do

! principal quantum number for all-electron well states
 do ii=nc,1,-1
   l1=la(ii)+1
   if(np(l1)==0) then
    np(l1)=na(ii)+1
   end if
 end do

 do l1=1,lmax+1
  if(np(l1)==0) then
   np(l1)=l1
  end if
 end do

 nproj(lloc+1)=0
 rc0(:)=rc(:)

!find log mesh point nearest input rc
 rcmax=0.0d0
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
! null charge and eigenvalue accumulators
 eeig=0.0d0
 zval=0.0d0
 rho(:)=0.0d0
 cvgplt(:,:,:,:,:)=0.0d0
!
! loop to construct pseudopotentials for all angular momenta
!
 write(6,'(/a/a)') 'Begin loop to  construct optimized pseudo wave functions',&
&      'and semi-local pseudopoentials for all angular momenta'

 do l1=1,lmax+1
  ll=l1-1
  uu(:,:)=0.0d0; uu2(:,:)=0.0d0; qq(:,:)=0.0d0

  if(l1==1) then
   mkap=1
  else
   mkap=2
  end if
! loop on J = ll +/- 1/2
  do ikap=1,mkap
   if(ikap==1) then
     kap=-(ll+1)
     fj=dble(ll+1)*fp(l1)/dble(2*ll+1)
   else
     kap=  ll
     fj=dble(ll)*fp(l1)/dble(2*ll+1)
   end if
! all-electron wave function for first projector
   if(fp(l1) .gt. 0.0d0) then
     call ldiracfb(np(l1),ll,kap,ierr,ep(l1,ikap), &
&                  rr,zz,vfull,uu,up,mmax,mch)
     if(ierr .ne. 0) then
       write(6,'(a,i2,a,i2)') 'oncvpsp 409: lschfb convergence error, n=',&
&       np(l1),'  l=',ll
       stop
     end if
   else
     call wellstate_r(np(l1),ll,kap,irc(l1),ep(l1,ikap), &
&                     rr,vfull,uu,zz,mmax)
   end if

! find outermost peak of wavefunction
   do jj=mch-1,1,-1
     if(up(jj,1)*up(jj+1,1)<0.0d0) then
       rpkt=rr(jj)
       exit
     end if
   end do

   call renorm_r(uu,rr,ll,kap,zz,mmax,cnorm)

   write(6,'(/a/a,f12.8,a,i2,a,i2,a,i2)') &
&    '  All-electron reference function large-component-only renormalization',&
&    '    cnorm= ',cnorm,'  l= ',ll,'  n= ',np(l1),'  kap= ',kap

   call fpovlp_r(uu,uu,irc(l1),ll,kap,zz,qq(1,1),rr)

! all-electron wave function for second projector
   if(nproj(l1)==2) then
     et=ep(l1,ikap)+debl(l1,ikap)
     if(et<0.0d0) then
       call ldiracfb(np(l1)+1,ll,kap,ierr,et, &
&                    rr,zz,vfull,uu2,up2,mmax,mch)
       if(ierr .ne. 0) then
         write(6,'(a,i2,a,i2)') 'oncvpsp 441: lschfb convergence error, n=',&
&         np(l1),'  l=',ll
         stop
       end if
       if(abs(et-ep(l1,ikap)-debl(l1,ikap))>eps) then
         write(6,'(/a,f10.6,a/a,f10.6,a/a)') &
&              'WARNING: negative energy of ',ep(l1,ikap)+debl(l1,ikap), &
&              ' specified for second projector', &
&              ' Actual bound-state energy of ',et, ' will be used', &
&              ' Action: debl should probably be increased'
         debl(l1,ikap)=et-ep(l1,ikap)
       end if
     else
       call wellstate_r(np(l1)+1,ll,kap,irc(l1),et, &
&                       rr,vfull,uu2,zz,mmax)
     end if

     call renorm_r(uu2,rr,ll,kap,zz,mmax,cnorm)

     write(6,'(/a/a,f12.8,a,i2,a,i2,a,i2)') &
&    '  All-electron reference function large-component-only renormalization',&
&    '    cnorm= ',cnorm,'  ll= ',ll,',  n= ',np(l1)+1,'  kap= ',kap

     call fpovlp_r(uu,uu2,irc(l1),ll,kap,zz,qq(1,2),rr)
     call fpovlp_r(uu2,uu2,irc(l1),ll,kap,zz,qq(2,2),rr)
   end if

   write(6,'(/a,f6.3,a,i2)') &
&        'First projector wave function outermost peak radius = ',&
&        rpkt,', l=',ll

   write(6,'(/a,i2/a,f12.8,a,f12.8,a,f12.8)') &
&        '  All-electron projector functions [0,rc] norms and overlap, l=',&
&        ll,'    qq(1,1)=',qq(1,1),'  qq(1,2)=',qq(1,2),'  qq(2,2)=',qq(2,2)

! calculate optimized pseudo wave functions and related potentials
   call run_optimize(ep(l1,ikap),ep(l1,ikap)+debl(l1,ikap),ll,mmax,rr, &
&                    uu(1,1),uu2(1,1),qq, &
&                    irc(l1),qcut(l1),ncon(l1),nbas(l1),nproj(l1), &
&                    pswf(1,1,l1,ikap),vp(1,l1,ikap),vkb(1,1,l1,ikap),vfull, &
&                    cvgplt(1,1,1,l1,ikap))

  end do !ikap
 end do !l1

! construct Vanderbilt / Kleinman-Bylander projectors

 write(6,'(/a,a)') 'Construct Vanderbilt / Kleinmman-Bylander projectors'

 call run_vkb_r(lmax,lloc,lpopt,dvloc0,irc,nproj,rr,mmax,pswf,vfull,vp, &
&             evkb,vkb,nlim)


! accumulate charge and eigenvalues
! pseudo wave functions are calculated with VKB projectors for
! maximum consistency of unscreening
 do l1=1,lmax+1
  ll=l1-1
  if(ll==0) then
    mkap=1
  else
   mkap=2
  end if
  do ikap=1,mkap
   if(ikap==1) then
     kap=-(ll+1)
     fj=dble(ll+1)*fp(l1)/dble(2*ll+1)
   else
     kap=  ll
     fj=dble(ll)*fp(l1)/dble(2*ll+1)
   end if
   if(fp(l1) .gt. 0.0d0) then
     etest=ep(l1,ikap)
     sls=(l1-1)*l1
     emax=vp(mmax,l1,ikap)+0.5d0*sls/rr(mmax)**2
     emin=emax
     do jj=1,mmax
       emin=dmin1(emin,vp(jj,l1,ikap)+0.5d0*sls/rr(jj)**2)
     end do
     call lschvkbb(ll+1,ll,nproj(l1),ierr,etest,emin,emax, &
&                  rr,vp(1,lloc+1,ikap),vkb(1,1,l1,ikap),evkb(1,l1,ikap), &
&                  uu(1,1),up(1,1),mmax,mch)
     if(ierr .ne. 0) then
       write(6,'(a,i2,a,i2)') 'oncvpsp: lschvkbb convergence error, n=',&
&       ll+1,'  l=',ll
       stop
     end if
     eeig=eeig+fj*etest
     zval=zval+fj
     rho(:)=rho(:)+fj*(uu(:,1)/rr(:))**2
   end if
  end do !ikap
 end do !l

! "extra" valence states with single-node pseudowavefunctions add
! screening charge
!
 if(nxtra>0) then
   do ii=1,nxtra
     ll=lpx(ii)
     l1=ll+1
     if(ll==0) then
      mkap=1
     else
      mkap=2
     end if
! loop on J = ll +/- 1/2
     do ikap=1,mkap
      if(ikap==1) then
        kap=-(ll+1)
        fj=dble(ll+1)*fpx(ii)/dble(2*ll+1)
      else
        kap=  ll
        fj=dble(ll)*fpx(ii)/dble(2*ll+1)
      end if
      etest=epx(l1,ikap)
      sls=(l1-1)*l1
      emax=vp(mmax,l1,ikap)+0.5d0*sls/rr(mmax)**2
      emin=emax
      do jj=1,mmax
        emin=dmin1(emin,vp(jj,l1,ikap)+0.5d0*sls/rr(jj)**2)
      end do
      call lschvkbb(ll+2,ll,nproj(l1),ierr,etest,emin,emax, &
&                   rr,vp(1,lloc+1,ikap),vkb(1,1,l1,ikap),evkb(1,l1,ikap), &
&                   uu(1,1),up(1,1),mmax,mch)
      if(ierr .ne. 0) then
        write(6,'(a,i2,a,i2)') 'oncvpsp: lschvkbb convergence error, n=',&
&        ll+2,'  l=',ll
        stop
      end if
      eeig=eeig+fj*etest
      rho(:)=rho(:)+fj*(uu(:,1)/rr(:))**2
    end do !ikap
    zval=zval+fpx(ii)
   end do !ii
 end if
!
 allocate(rhomod(mmax,5))

 rhomod(:,:)=0.0d0

! construct model core charge based on monotonic polynomial fit
 if(icmod==1 .or. icmod==2) then
   call modcore(rho,rhoc,rhomod,fcfact,irct,mmax,rr)
 end if

! screening potential for pseudocharge

 call vout(1,rho,rhomod(1,1),vo,vxc,zval,eeel,eexc,rr,mmax,iexc)

! total energy output

 epstot= eeig + eexc - 0.5d0*eeel
 write(6,'(/a,f12.6/)') 'Pseudoatom total energy', epstot

 call run_diag_r(lmax,np,fp,ep,debl,lloc,irc,npx,lpx,fpx,nxtra, &
&              vkb,evkb,nproj,rr,vfull,vp,zz,mmax)

! for abinit output, need decomposition into scalar-relativistic and
! spin-orbit projectors

 if(trim(psfile)=='psp8') then
  call sr_so_r(lmax,irc,nproj,rr,mmax,evkb,vkb, &
&                vsr,esr,vso,eso)

! drop sr, so orthonormal projectors with neglibible coefficients
! moldify cutoff if desired

 eprmin=2.0d-5
 write(6,'(/a,1p,e10.2,a)') 'Orthonormal projectors with coefficients <', &
&      eprmin,' Ha will be dropped'

 do l1=1,lmax+1
  if(abs(esr(3,l1))<eprmin) esr(3,l1)=0.0d0
  if(abs(esr(4,l1))<eprmin) esr(4,l1)=0.0d0
  if(abs(eso(3,l1))<eprmin) eso(3,l1)=0.0d0
  if(abs(eso(4,l1))<eprmin) eso(4,l1)=0.0d0
 end do

! diagnostics are re-run for this case for error checking and to check
! the errors introduced by truncation, if any from 4 sr and 4 so projectors

   call run_diag_sr_so_r(lmax,np,fp,ep,debl,lloc,irc,npx,lpx,fpx,nxtra, &
&                      vsr,esr,vso,eso,nproj,rr,vfull,vp,zz,mmax)

 end if

! unscreen semi-local potentials taking scalar-relativistic average

 do l1=1,max(lmax+1,lloc+1)
   ll=l1-1
   if(ll==0  .or. ll==4) then
     vpuns(:,l1)=vp(:,l1,1)-vo(:)
   else
     vpuns(:,l1)=(dble(ll+1)*vp(:,l1,1)+ dble(ll)*vp(:,l1,2))/dble(2*ll+1) &
&                -vo(:)
   end if
 end do

! loop over reference plus test atom configurations
 do jj=1,ncnf+1
   do ii=nv+1,nvcnf(jj)
     ea(ii,:)=ea(nv,:)
   end do

   write(6,'(/a,i2)') 'Test configuration',jj-1
   call run_config_r(nacnf(1,jj),lacnf(1,jj),ea,facnf(1,jj),nc,nvcnf(jj), &
&                  rhomod,rr,zz,rcmax,mmax,iexc,etot,epstot,nproj,vpuns, &
&                  lloc,vkb,evkb)
 end do

!create potential and wave function data for plotting
 call run_plot_r(lmax,np,fp,ep,lloc,rc,npx,lpx,fpx,epx,nxtra, &
&              vkb,evkb,nproj,rr,vfull,vp,vpuns,zz,mmax,drl,nrl, &
&              rho,rhoc,rhomod,pswf,cvgplt)

!set radius for log derivatives, max(rc) if no cores are included,
! max(rc)+1.0 a_b if they are.
!
 if(nxtra==0) then
  irpsh=nlim
 else
  irpsh=nint(dlog((rr(nlim)+1.0d0)/rr(1))/al)
 end if

!create arctan(log derivative) data for plotting
 call run_phsft_r(lmax,lloc,nproj,ep,epsh1,epsh2,depsh,vkb,evkb, &
&               rr,vfull,vp,zz,mmax,irpsh)

!create script for automated plotting
 call gnu_script_r(lmax,lloc,nproj,nxtra,lpx)


!psp output for Abinit or Pwscf
 if(trim(psfile)=='psp8') then
  soscale=1.0d0  !change to rescale spin-orbit strength
  call linout_r(lmax,lloc,rc,vsr,esr,vso,eso,rr,vpuns,rho,rhomod, &
&                  zz,zion,mmax,iexc,fcfact,nrl,drl,atsym,soscale)

 else if(trim(psfile)=='upf') then
  call upfout_r(lmax,lloc,rc,vkb,evkb,nproj,rr,vpuns,rho,rhomod, &
&             zz,zion,mmax,iexc,icmod,nrl,drl,atsym,epstot, &
&             na,la,ncon,nbas,nvcnf,nacnf,lacnf,nc,nv,lpopt,ncnf, &
&             fa,rc0,ep,qcut,debl,facnf,dvloc0,fcfact, &
&             epsh1,epsh2,depsh,rlmax,psfile)
 end if

 stop
 end program oncvpsp_r


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
  write(6,'(a,i4)') 'Read error, input data file line',inline
  write(6,'(a)') 'Program will stop'
  stop
 end if

 return
 end subroutine read_error
