program Ionboost

    use iso_fortran_env, only: int32, real32, int64, real64
    use mod_math, only: erfc, asinh, Thot, Tcold

    implicit none
    integer(int32), parameter :: nmax=150100, njmax=200
    integer(int32) :: ntotal, nsort, iline, nstep, iline_vide, itime
    integer(kind=4) ::  idebut(0:nmax),irang(0:nmax),m,djout
 integer(int32) :: i, j, idico, iter
   
    real(real64) :: time, dt, cs, Th, Tc, cs2, cs2old, PI
    real(real64) :: rgauss, alpha, Whot1, Whot2, Wcold
    real(real64) :: length,ni0,n0cold,n0hot,neh,nec,neold,kvide
    real(real64) :: nti,nthot,ntcold,nte,lmax,nte1,nte2
    real(real64) :: nthot0,LSS,nombre,xvide
    real(real64) :: nhm3,nhm2,nhm1,nu,nnhot,nncold
    real(real64) :: lDebye,lgrade,lgradi,nLSS,tb,vmax0,nii,vv,qq,k1,k2
    real(real64) :: k3,k4,nhot1,vmax,Emax,C1,C2,D1,D2,a0,a1,a2,vfinal,vfinal2
    real(real64) :: borne2,borne1,Integ,a3,a4,xr,temp,ekmoy,Inv2,vintr,E1s2r
    real(real64) :: p1,p2,p3,p4,Hold,Hnew,venew,phi0,phiconv,slop,Et1,Et2,t1,t2
    real(real64) :: phimoy,path,xbar,Rbar,phibar,Integ3,cs3,cscin,csmoyen
    real(real64) :: fe_min, prog_v
    real(real64) ::   x0 (0:nmax), xt (0:nmax), xtold(0:nmax)
    real(real64) ::    v  (0:nmax), vint (0:nmax), vemoy(0:nmax)
    real(real64) ::    v2(0:nmax), vintold(0:nmax)
    real(real64) ::    E  (0:nmax),  Eold (0:nmax), E2(0:nmax)
    real(real64) ::    phi (0:nmax),  ni (0:nmax)
    real(real64) ::    nilisse (0:nmax),  nelisse (0:nmax)
    real(real64) ::    phiold (0:nmax), phi2old(0:nmax)
    real(real64) ::    phi1m(0:nmax), phi2m(0:nmax), phi3m(0:nmax)
    real(real64) ::    ne (0:nmax),  rho (0:nmax)
    real(real64) ::    nhot (0:nmax),  ncold (0:nmax)
    real(real64) ::    nhotold (0:nmax),  ncoldold (0:nmax)
    real(real64) ::    gradnh (0:nmax),  gradnc (0:nmax)
    real(real64) ::    ghold (0:nmax),  gcold (0:nmax)
    real(real64) ::    pe (0:nmax)
    real(real64) ::    nss (0:nmax)
    real(real64) ::    dWhot (0:nmax), SWhot(0:nmax)
    real(real64) ::    phi2(0:nmax), xt2(0:nmax), Enhot(0:nmax), csc(0:nmax)
    real(real64) ::    En(0:nmax), Encold(0:nmax)

    real(real64) ::   vmoy (1:nmax), Emoy (1:nmax), dndv (1:nmax), dndE (1:nmax)
    real(real64) ::   dx0  (1:nmax),  dxt (1:nmax), dxtold(1:nmax)
    real(real64) ::   ni1s2 (1:nmax), E1s2 (1:nmax), E1s2old(1:nmax)
    real(real64) ::   niSS (0:nmax), qiSS(0:nmax), dxt2(1:nmax)
      
    real(real64) :: a (1:nmax), b(0:nmax), c(0:nmax-1), f(0:nmax),bx(0:nmax), fx(0:nmax)

!real(real64) :: dve2s2(0:njmax), dve2s2_0(0:njmax)
    real(real64) :: ve (0:njmax), fe0(0:njmax), fe(0:njmax)
    real(real64) :: Te(0:njmax), Te_0(0:njmax)
    real(real64) :: veold(0:njmax), ue(0:njmax)
    real(real64) :: ve1m(0:njmax), ve2m(0:njmax),ve3m(0:njmax)
    real(real64) :: Inv0(0:njmax)

    real(real64) :: fe0_cold(0:njmax), fe_cold(0:njmax), Te_cold(0:njmax)
    real(real64) :: fe0_hot(0:njmax), fe_hot(0:njmax),Te_hot(0:njmax),dve2s2_0(0:njmax)
    real(real64) :: dve, nephi, gnephi, pephi2, pephi, gpephi, H, period, Spe, Spe0, deltaE, poly, disc, Spe1
    real(real64) :: nehot, necold
     
    logical lfini,nb_cons,En_cons,EOS,VTS,laststep,break,bitemp

    character*10 profil
!*               1 - demarrage et conditions initiales
    
!*         1.1 lecture des donnees
      
    integer(int32) :: ncell, nvide, nlisse,Nbe, itmax, jcoldmax
    real(real64)   :: ve_max, prog, iter0, iter1, fiter2, fiter3, dti,tmax
    real(real64)   :: Tcmax,Thmax, T_MeV
real(real64)    :: omegadt

    open(unit=10, status='unknown', file='boost72.par')


    read(10,*) ncell
    read(10,*) nvide
    read(10,*) nlisse
    read(10,*) ve_max
    read(10,*) Nbe
    read(10,*) prog
    read(10,*) itmax
    read(10,*) tmax
    read(10,*) iter0
    read(10,*) iter1
    read(10,*) fiter2
    read(10,*) fiter3
    read(10,*) dti
    read(10,*) n0cold
    read(10,*) n0hot
    read(10,*)  Thmax
    read(10,*)  Tcmax
    read(10,*)  lfini
    read(10,*)  lmax
    read(10,*)  nb_cons
    read(10,*)  En_cons
    read(10,*)  EOS
    read(10,*)  T_MeV
    read(10,*)  LSS
    read(10,*)  nLSS
    read(10,*)  nu
    read(10,*)  profil
    read(10,*)  VTS
    read(10,*)  djout
    close(13)

    if(prog.lt.0.9)then
       prog=10.0**(-5.0/ncell)
    endif
    
     if(profil.ne.'sack'.and.profil.ne.'step'.and.profil.ne.'expo'.and.profil.ne.'gaus')then
       write (*,*) 'What is the initial density profile'
       stop
     endif

    if((.not.lfini).and.(nb_cons.or.En_cons))then
       nb_cons=.false.
       En_cons=.false.
       print*,    'lfini=.false. therefore En_cons=.false. and nb_cons=.false.'
     endif
     
     if((.not.nb_cons).and.En_cons) then
       nb_cons=.true.
       print*,    'En_cons=.true.  therefore  nb_cons=.true.'
     endif
     
     if(nb_cons.and.(n0cold.lt.1.d-04))then
       nb_cons=.false.
       print*,    'n0cold.lt.1.d-04  therefore  nb_cons=.false.'
     endif
     
     if((.not.lfini).and.VTS) then
       VTS=.false.
       print*,    'lfini=.false. therefore VTT=.false.'
     endif
     
     ntotal=ncell+nvide
     
     if(ntotal.gt.nmax)then
       print*, 'ntotal=ncell+nvide is bigger than nmax'
       stop
     endif


    omegadt=sqrt(n0cold+n0hot)*dti
    if(omegadt.gt.2.) then
      print*, 'No reduced time because of stability issue'
      dti=2./sqrt(n0cold+n0hot)
      write(*,*)'dti=',dti
    endif

    bitemp=.false.
    if(Tcmax.GT.0.and.n0cold.gt.0.) bitemp=.true.

    nstep=(tmax/dti+0.5d0)
    nsort=1
    if(nstep.gt.2000)nsort=nstep/2000
    laststep=.false.

    iline=1
    iline_vide=1
    if(ncell.gt.3999)iline=ncell/2000
    if(nvide.gt.799)iline_vide=nvide/400


!*         1.2 construction du maillage
      
!*     1.2.1 determination de la longueur du plasma

      if(.not.lfini) then
        time=0.0
        dt=dti
        cs2=(n0cold+n0hot)/(n0cold/Tcmax+n0hot/Tcmax)
        cs=sqrt(cs2)
        if(profil.eq.'sack')length=(10.d0*LSS+25.d0)*cs
        if(profil.eq.'expo')length=(LSS+25.d0)*cs
        if(profil.eq.'step')length=25.d0*cs

        if(itmax.eq.1) goto 91
!
        Th=Thmax*Thot(0.d0)
        Tc=Tcmax*Tcold(0.d0)
!
        do 90 itime=1,itmax
!
          if(time.ge.(tmax-1.d-5)) goto 91
          if(itime.eq.itmax) goto 91
          if((time+dt).gt.(tmax-1.d-5)) dt=tmax-time

          cs2old=cs2
!
          time=time+dt
!
          Th=Thmax*Thot(time)
          Tc=Tcmax*Tcold(time)
          cs2=(n0cold+n0hot)/(n0cold/Tc+n0hot/Th)

          cs=0.5d0*(sqrt(cs2old)+sqrt(cs2))
          length=length+cs*dt
!
 90     continue
 91     continue
      endif
    
    if(lfini) length=lmax


!*     1.2.2 maillage spatial

    PI=4.D0*DATAN(1.D0)
    rgauss=2.d0*length/sqrt(PI)
    ni0=n0cold+n0hot

    if(prog.eq.1.d0) then
        dx0(1)=length/ncell
    else
        dx0(1)=length*(1.-prog)/(1.-prog**ncell)
    endif
        x0(0)=-length

    x0(1)=x0(0)+dx0(1)

    do i = 0,1
      niSS(i)=ni0
      if(profil.eq.'expo'.and.x0(i).gt.-LSS)  niSS(i)=ni0*exp(-(x0(i)+LSS)/LSS)
      if(profil.eq.'sack') niSS(i)=ni0*(2.d0/pi)*DATAN(exp(-x0(i)/LSS))
      if(profil.eq.'gaus') niSS(i)=ni0*exp(-((x0(i)+length)/rgauss)**2)
    end do

    do i = 2,ncell
      if(i.eq.2) then
        dx0(2)=(1.d0+prog)*dx0(1)*niSS(0)/niSS(1)-dx0(1)
      else
        dx0(i)=(dx0(i-2)+dx0(i-1))*prog*niSS(i-2)/niSS(i-1)-dx0(i-1)
      endif
      if((profil.eq.'expo'.or.profil.eq.'sack') .and.dx0(i).gt.dx0(i-1).and.dx0(i).gt.LSS/5.d0)  dx0(i)=LSS/5.d0
      if(profil.eq.'gaus'.and.dx0(i).gt.dx0(i-1).and.dx0(i).gt.lmax/5.d0)  dx0(i)=lmax/5.d0

      x0(i)=x0(i-1)+dx0(i)
      niSS(i)=ni0
      if(profil.eq.'expo'.and.x0(i).gt.-LSS) niSS(i)=ni0*exp(-(x0(i)+LSS)/LSS)
      if(profil.eq.'sack') niSS(i)=ni0*(2.d0/pi)*DATAN(exp(-x0(i)/LSS))
      if(profil.eq.'gaus') niSS(i)=ni0*exp(-((x0(i)+length)/rgauss)**2)
    end do


     qiSS(0)=niSS(0)*dx0(1)
     do i = 1,ncell-1
        qiSS(i)=niSS(i)*(dx0(i)+dx0(i+1))/2.d0
    end do
    qiSS(ncell)=niSS(ncell)*dx0(ncell)*(1.d0+prog)/2.
!
!    *     1.2.3 mise en memoire de l'ordre initial des ions
!
    do i = 0,ncell
       idebut(i)=i
       irang(i)=i
    end do



!*         1.3 conditions initiales et estimation du
!*             potentiel phi initial
!
        break=.false. !pas de deferlement au debut

    Th=Thmax*Thot(0.d0)
    Tc=Tcmax*Tcold(0.d0)
      if(Th.eq.0.d0) then
         write(*,*) 'The hot temperature must be > 0.'
         stop
      endif


    phi(ncell)=(n0cold*Tc+n0hot*Th)/(n0cold+n0hot)
    E(ncell)=sqrt(2.*(n0cold*Tc*exp(-phi(ncell)/Tc)+n0hot*Th*exp(-phi(ncell)/Th)))
    alpha=E(ncell)/phi(ncell)

    do i = 0,ncell
       xt(i)=x0(i)
       v(i)=0.d0
           vemoy(i)=0.d0
       phi(i)=phi(ncell)*exp(alpha*(x0(i)-x0(ncell)))
    end do


    do i =ncell+1,ntotal
       v(i)=0.d0
           vemoy(i)=0.d0
           v2(i)=0.d0
       ni(i)=1.d-12
    enddo

    E(0)=0.d0
    gradnh(0)=0.d0
    gradnc(0)=0.d0
    time=0.d0
    dt=dti
    nhm2=n0hot
    nhm1=n0hot
    Whot1=0.d0
    Whot2=0.d0
    Wcold=0.d0
    idico=0

!*     1.3.1 initialisation de la fonction de distribution des electrons
!*             et des temperatures locales

    xvide= ve_max**2/sqrt(8.d0)+0.3d0
    ve_max=ve_max*sqrt(Thmax)
    fe_min=n0hot/sqrt(Thmax)*exp(-ve_max**2/Thmax/2.d0)

    if(Thmax.eq.Tcmax.or.(.not.bitemp)) then
       prog_v=1.d0
       dve=ve_max/(Nbe-1)
    else
       prog_v=(Thmax/Tcmax)**(0.5d0/(Nbe-2))
       dve=ve_max*(prog_v-1.d0)/(prog_v**(Nbe-1)-1.d0)
    endif


    jcoldmax=Nbe
       do j=0,Nbe-1
      if(j.eq.0)then
        ve(0)=0.d0
      else
        ve(j)=ve(j-1)+dve
        dve=dve*prog_v
      endif
          fe0_hot(j)=n0hot/sqrt(Thmax)*exp(-(ve(j)**2)/Thmax/2.d0)
      if(bitemp) then
        fe0_cold(j)=n0cold/sqrt(Tcmax)*exp(-(ve(j)**2)/Tcmax/2.d0)
          else
            fe0_cold(j)=0.d0
      endif
      if(fe0_cold(j).lt.fe_min.and.j.lt.jcoldmax)then
        jcoldmax=j
      endif
          fe0(j)=fe0_hot(j)+fe0_cold(j)
      fe_hot(j)=fe0_hot(j)
      fe_cold(j)=fe0_cold(j)
      fe(j)=fe0(j)
        enddo


    do j=0,Nbe-2
      Te(j)=-0.5d0*(ve(j+1)**2-ve(j)**2)/(log(fe(j+1))-log(fe(j)))
      Te_0(j)=Te(j)
      Te_hot(j)= Thmax
      Te_cold(j)= Tcmax
    enddo
    Te(Nbe-1)=Thmax
    Te_hot(Nbe-1)=Thmax
    Te_cold(Nbe-1)=Tcmax

    do j=0,Nbe-1
      write(*,*) ve(j)
      write(*,*) fe0(j), fe0_hot(j), fe0_cold(j)
      write(*,*) Te(j), Te_hot(j), Te_cold(j)
      write(*,*)
    enddo

    do j=0,Nbe-1
       ve2m(j)=ve(j)
       ve1m(j)=ve(j)
    enddo




end program Ionboost
