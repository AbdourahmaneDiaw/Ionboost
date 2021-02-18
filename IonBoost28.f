c***5***10***15***20***25***30***35***40***45***50***55***60***65***70**
*
*		Program IonBoost28.f  -  3 aout 2009
*
*       Calcul de l'expansion  dans le vide d'un plasma de temperature
*         homogene
*       Les ions sont traites comme des plans charges
*       Les electrons sont en equilibre de Boltzmann
*	La separation de charge est prise en compte par la resolution
*	    de l'equation de Poisson
*
*	La temperature est fonction du temps
*	  et on calcule la gaine electronique dans la partie vide d'ions
*
*	La fonction de distribution des electrons
*	    est bi-maxwellienne.
*
*       On tient compte de la conservation
*         - du nombre de particules de chaque espece
*         - de l'energie
*
*       On traite le deferlement en reordonnant les particules
*	  - idebut correspond au rang initial d'une particule en
*		fonction de son rang final
*	  - irang correspond au contraire au rang final d'une
*		particule en fonction de son rang initial
*
*	On inclut la possibilite d'avoir une viscosite "nu"
*
*	Differentes formes de profil initial sont possibles.
*	Le maillage est redefini.
*
*	On inclut des ions test de Z/m differents
*        
*       Chaque ion a maintenant une charge Z/m
*
*       On inclut une premiere phase, correspondant au temps de 
*           montee du laser
*
*	On inclut des ions négatifs infiniment lourds
*
*
c***5***10***15***20***25***30***35***40***45***50***55***60***65***70**

      parameter (nmax=180000)
      implicit double precision (a-h,o-z)
      
      integer idebut(0:nmax),irang(0:nmax),itest(0:nmax)

      double precision length,ni0,n0cold,n0hot,neh,nec,neold,kvide
     *                 ,nti,nthot,ntcold,nte,lmax
     *                 ,nthot0
     *                 ,nhm3,nhm2,nhm1,LSS,nu,nnhot,nncold
     *		       ,lDebye,lgrade,lgradi,lay1,lay2,mix,mixp	
    
      double precision   x0 (0:nmax), xt (0:nmax), xtold(0:nmax)
     *         , v  (0:nmax),  vint (0:nmax)
     *         , E  (0:nmax),  Eold (0:nmax)
     *         , phi (0:nmax),  ni (0:nmax)
     *         , nilisse (0:nmax),  nelisse (0:nmax)
     *         , phiold (0:nmax)
     *         , ne (0:nmax),  rho (0:nmax)
     *         , nhot (0:nmax),  ncold (0:nmax)
     *         , nhotold (0:nmax),  ncoldold (0:nmax)
     *         , gradnh (0:nmax),  gradnc (0:nmax)
     *         , ghold (0:nmax),  gcold (0:nmax)
     *         , phot (0:nmax),  pcold (0:nmax)
     *         , pe (0:nmax)
     *         , nss (0:nmax),  difth (0:nmax)
     *         , dWhot (0:nmax), SWhot(0:nmax)
     *         , x0test (0:nmax), xttest (0:nmax), xttestold(0:nmax)
     *         , vtest (0:nmax), vtestint (0:nmax)
     *         , Etest  (0:nmax),  Etestold (0:nmax)
     
      double precision vmoy (1:nmax), Emoy (1:nmax)
     *         ,dndv (1:nmax), dndE (1:nmax)

      double precision dx0  (1:nmax),  dxt (1:nmax), dxtold(1:nmax)
     *         ,ni1s2 (1:nmax), E1s2 (1:nmax)
     *         , niSS (0:nmax), qiSS(0:nmax)
     *         , charge(0:nmax)
      
      double precision a (1:nmax), b(0:nmax), c(0:nmax-1), f(0:nmax)
     *         ,bx(0:nmax), fx(0:nmax)
     
      logical lfini,nb_cons,En_cons,EOS,VTT,laststep,multicouche
     *        ,ions_negatifs,multiphase

      character*10 profil
     
      external Thot,Tcold

      namelist/par/ncell,nvide,nlisse,prog,itmax,tmax
     *            ,iter0,iter1,iter2,iter3
     *            ,dti,n0cold,n0hot,Thmax,Tcmax,lfini,lmax,nb_cons
     *            ,En_cons,EOS,T_MeV,LSS,nLSS,nu,profil,VTT,Ztest
     *            ,multicouche,lay1,lay2,mix,charge2
     *            ,ions_negatifs,p_negatif
     *            ,multiphase,trise
      
c***5***10***15***20***25***30***35***40***45***50***55***60***65***70**
      
	
*               1 - demarrage et conditions initiales
	
*         1.1 lecture des donnees
      
      open (13,file='boost.par',status='unknown',form='formatted')
      read (13,par)
      close(13)
      
      write(*,par)
      
      if(profil.ne.'sack'.and.profil.ne.'step'.and.profil.ne.'expo'
     *   .and.profil.ne.'gaus')
     *   then
         write (*,*) 'c''est quoi le profil?'
         stop
      endif
      
      if((.not.lfini).and.(nb_cons.or.En_cons))then
         nb_cons=.false.
         En_cons=.false.
         write(*,*) 
     *    'lfini=.false. et donc En_cons=.false. et nb_cons=.false.'
      endif
      
      if((.not.nb_cons).and.En_cons) then
         nb_cons=.true.
         write(*,*) 'En_cons=.true.  et donc  nb_cons=.true.'
      endif
      
      if(nb_cons.and.(n0cold.lt.1.d-04))then
         nb_cons=.false.
         write(*,*) 'n0cold.lt.1.d-04  et donc  nb_cons=.false.'
      endif
      
      if((.not.lfini).and.VTT) then
         VTT=.false.
         write(*,*) 'lfini=.false. et donc VTT=.false.'
      endif
      
      ntotal=ncell+nvide
      
      if(ntotal.gt.nmax)then
         write(*,*) 'ntotal=ncell+nvide est plus grand que nmax'
         stop
      endif
	
	nstep=(tmax/dti+0.5d0)
	nsort=1
	if(nstep.gt.2000)nsort=nstep/2000
	laststep=.false.
	
	iline=1
	if(ncell.gt.3999)iline=ncell/2000
	
*         1.2 construction du maillage
      
*     1.2.1 determination de la longueur du plasma

      if(.not.lfini) then
         time=0.d0
         dt=dti
         Th=Thmax
         Tc=Tcmax
         cs2=(n0cold+n0hot)/(n0cold/Tc+n0hot/Th)
         cs=sqrt(cs2)
         if(profil.eq.'sack')length=10.d0*LSS+25.d0*cs
         if(profil.eq.'expo')length=LSS+25.d0*cs
         if(profil.eq.'step')length=25.d0*cs
         
         if(itmax.eq.1) goto 91
	
         Th=Thmax*Thot(multiphase,trise,0.)
         Tc=Tcmax*Tcold(0.)

         do 90 itime=1,itmax
	
            if(time.ge.(tmax-1.d-5)) goto 91
            if(itime.eq.itmax) goto 91
            if((time+dt).gt.(tmax-1.d-5)) dt=tmax-time
	
            cs2old=cs2
	
            time=time+dt
	
            Th=Thmax*Thot(multiphase,trise,time)
            Tc=Tcmax*Tcold(time)
            cs2=(n0cold+n0hot)/(n0cold/Tc+n0hot/Th)
	
            cs=0.5d0*(sqrt(cs2old)+sqrt(cs2))
            length=length+cs*dt

 90      continue
 91      continue
      endif
	
      if(ions_negatifs)length=length*sqrt(1.d0/(1.d0-p_negatif))
      if(lfini)length=lmax
      
*     1.2.2 maillage spatial

      PI=4.D0*DATAN(1.D0)
      rgauss=2.d0*length/sqrt(PI)
	ni0=n0cold+n0hot
	      
	if(prog.eq.1.d0)then
         dx0(1)=length/ncell
	else
         dx0(1)=length*(1.-prog)/(1.-prog**ncell)
	endif
      x0(0)=-length
	x0(1)=x0(0)+dx0(1)
      
	do i = 0,1
	   niSS(i)=ni0
	   if(profil.eq.'expo'.and.x0(i).gt.-LSS) 
     *      niSS(i)=ni0*exp(-(x0(i)+LSS)/LSS)
	   if(profil.eq.'sack') 
     *      niSS(i)=ni0*(2.d0/pi)*DATAN(exp(-x0(i)/LSS))
	   if(profil.eq.'gaus') 
     *      niSS(i)=ni0*exp(-((x0(i)+length)/rgauss)**2)
	end do

	do i = 2,ncell
	   if(i.eq.2)then
	     dx0(2)=(1.d0+prog)*dx0(1)*niSS(0)/niSS(1)-dx0(1)
	   else
	     dx0(i)=(dx0(i-2)+dx0(i-1))*prog*niSS(i-2)/niSS(i-1)-dx0(i-1)
	   endif
	   if((profil.eq.'expo'.or.profil.eq.'sack')
     *      .and.dx0(i).gt.dx0(i-1)
     *      .and.dx0(i).gt.LSS/5.d0) dx0(i)=LSS/5.d0
	   if(profil.eq.'gaus'
     *      .and.dx0(i).gt.dx0(i-1)
     *      .and.dx0(i).gt.lmax/5.d0) dx0(i)=lmax/5.d0
	   x0(i)=x0(i-1)+dx0(i)
	   niSS(i)=ni0
	   if(profil.eq.'expo'.and.x0(i).gt.-LSS) 
     *      niSS(i)=ni0*exp(-(x0(i)+LSS)/LSS)
	   if(profil.eq.'sack') 
     *      niSS(i)=ni0*(2.d0/pi)*DATAN(exp(-x0(i)/LSS))
	   if(profil.eq.'gaus') 
     *      niSS(i)=ni0*exp(-((x0(i)+length)/rgauss)**2)
	end do
			
	qiSS(0)=niSS(0)*dx0(1)
	do i = 1,ncell-1
	   qiSS(i)=niSS(i)*(dx0(i)+dx0(i+1))/2.d0
	end do
	qiSS(ncell)=niSS(ncell)*dx0(ncell)*(1.d0+prog)/2.
      
*     1.2.3 mise en memoire de l'ordre initial des ions

	do i = 0,ncell
	   idebut(i)=i
	   irang(i)=i
	end do

	
*         1.3 conditions initiales et estimation du
*             potentiel phi initial
	
	Tnorm=T_MeV/.511d0
	e0=0.5*Tnorm*(1.+9.*Tnorm/4.+3.*Tnorm*Tnorm/4.)
     *             /(1.+3.*Tnorm/2.+3.*Tnorm*Tnorm/8.)
      if(EOS) then
         Tn0=2.*e0*(1.+9.*e0/2.+3.*e0*e0/2.)/(1.+6.*e0+3.*e0*e0)
      else
         Tn0=2.*e0
      endif

	Th=Thmax*Thot(multiphase,trise,0.)
	Tc=Tcmax*Tcold(0.)
      if(Th.eq.0.d0)then
         write(*,*) 'la temperature chaude ne doit pas s''annuler'
         stop
      endif
	
	phi(ncell)=(n0cold*Tc+n0hot*Th)/(n0cold+n0hot)
	E(ncell)=sqrt(2.*(n0cold*Tc*exp(-phi(ncell)/Tc)
     *                 +n0hot*Th*exp(-phi(ncell)/Th)))
      alpha=E(ncell)/phi(ncell)

          if(ions_negatifs)then
             mixp=-p_negatif/(1.d0-p_negatif)
          else
             mixp=0.d0
          endif

	do i = 0,ncell
	   xt(i)=x0(i)
	   v(i)=0.d0
	   phi(i)=phi(ncell)*exp(alpha*(x0(i)-x0(ncell)))
	   charge(i)=1.d0
          if(multicouche)then
             if(x0(i).gt.lay1.and.x0(i).lt.lay2) then
                if((i/2)*2.eq.i) then
                   qiSS(i)=qiSS(i)*2.d0*(1.d0-mix)
                else
                   charge(i)=charge2
                   qiSS(i)=qiSS(i)*2.d0*mix
                endif
             else
                if((i/2)*2.eq.i) then
                   qiSS(i)=qiSS(i)*(2.d0-1.d-8)
                else
                   charge(i)=charge2
                   qiSS(i)=qiSS(i)*1.d-8
                endif
             endif
          endif

          if(ions_negatifs)then
             qiSS(i)=qiSS(i)*(1.d0-mixp)
          endif

	end do
	
	do i = 0,ncell
	   x0test(i)=x0(i)
	   xttest(i)=x0(i)
	   vtest(i)=0.d0
	   itest(i)=i
	end do
	
	do i =ncell+1,ntotal
	   v(i)=0.
	   ni(i)=1.d-12
	enddo
	
	E(0)=0.d0
	Etest(0)=0.d0
	gradnh(0)=0.d0
	gradnc(0)=0.d0
	time=0.d0
	dt=dti
        if(multiphase.and.dt.gt.trise/100.d0)dt=trise/100.d0
	nhm2=n0hot
	nhm1=n0hot
	Thm2=Thmax*Thot(multiphase,trise,trise)
	Thm1=Thm2
	Tcm2=Tc
	Tcm1=Tc
	Whot1=0.d0
	Whot2=0.d0
	Wcold=0.d0
	idico=0
	
	open(9,file='conservation',status='unknown')
	write(9,'(37a)')'time',char(9), 'nti',char(9)
     *   ,'nthot', char(9), 'ntcold', char(9), 'nte'
     *   , char(9), 'n0hot', char(9), 'n0cold'
     *   , char(9), 'En_ion', char(9), 'Whot1'
     *   , char(9), 'Whot2', char(9), 'Whot', char(9), 'Wcold'
     *   , char(9), 'Th', char(9), 'Tc', char(9), 'En_elec'
     *   , char(9), 'En_delta'
     *   , char(9), 'vmax', char(9), 'vfinal'

	open(10,file='historique',status='unknown')
	write(10,'(29a)')'time',char(9), 'xi',char(9)
     *   ,'v(ivmax)', char(9), 'Energy_max', char(9)
     *   ,'E(ivmax)'
     *   , char(9), 'ne(ivmax)', char(9),'ni(ivmax)'
     *   , char(9),'ni(0)'
     *   , char(9),'nhot(0)',char(9),'ncold(0)'
     *   , char(9),'lDebye', char(9),'lgrade', char(9),'lgradi'
     *   , char(9),'ivmax',char(9),'idebut(ivmax)'
	
	open(14,file='histobis',status='unknown')
	write(14,'(9a)')'time',char(9), 't/R0',char(9), 'xi',char(9)
     *   ,'R/R0', char(9), 'Eq42'
	
*               2 - boucle temporelle
	
	do 100 itime=1,itmax
	
	if(itime.le.3) then
	   iter=iter0
	else
	   iter=iter1
	endif


*         2.1 determination ou estimation des temperatures et de n0hot
*             et calcul de la nouvelle densite ionique

*	2.1.1 - determination des temperatures dans le cas ou leur
*             dependance est fixee par les fonctions Thot et Tcold

	Thold=Th
	Tcoold=Tc
	
	if((multiphase.and.time.le.trise).or.(.not.En_cons))then
	   Th=Thmax*Thot(multiphase,trise,time)
	   Tc=Tcmax*Tcold(time)
	   idico=0
         if(Th.eq.0.d0)then
            write(*,*) 'la temperature chaude ne doit pas s''annuler'
            stop
         endif
	   if((Th-Thold)/Thold.gt.0.05d0) then
*	      Th=Thold*(1.05d0)
	      iter=iter0
*	      idico=1
	   endif
	   if((Th-Thold)/Thold.lt.-0.05d0) then
*	      Th=Thold*(0.95d0)
	      iter=iter0
*	      idico=1
   	   endif
	endif
		
*	2.1.2 - estimation de la modification de n0hot liee a la 
*             conservation du nombre d'electrons chauds

	if(nb_cons)then
	   nhm3=nhm2
	   nhm2=nhm1
	   nhm1=n0hot
	   n0hot=3.d0*nhm1-3.d0*nhm2+nhm3
	endif

*	2.1.3 - estimation de la modification de Th et Tc liee a la 
*              conservation de l'energie

	if(En_cons.and.(.not.multiphase.or.time.gt.trise))then
	   Thm3=Thm2
	   Thm2=Thm1
	   Thm1=Th
	   Th=3.d0*Thm1-3.d0*Thm2+Thm3
	   Tcm3=Tcm2
	   Tcm2=Tcm1
	   Tcm1=Tc
	   Tc=3.d0*Tcm1-3.d0*Tcm2+Tcm3
	endif
	
*****************************
	if(itime.eq.(itime/nsort)*nsort) then
	   write(*,*)
	   write(*,*)time
	endif
*	write(*,*)Th, Tc, n0hot
*****************************

*	2.1.4 - mise en memoire des anciennes valeurs de la densite
*             electronique, du potentiel, du champ electrique,
*             des gradients de densite electronique,
*             et de la position dans la gaine vide d'ions

      if(itime.gt.1) then
         do i=0,ntotal
            nhotold(i)=nhot(i)
            ncoldold(i)=ncold(i)
            phiold(i)=phi(i)
            Eold(i)=E(i)
            ghold(i)=gradnh(i)
            gcold(i)=gradnc(i)
         enddo
         do i=0,ncell
           Etestold(i)=Etest(i)
         enddo
         do i=ncell+1,ntotal
            xtold(i)=xt(i)
         enddo
         do i=1,ntotal
            dxtold(i)=dxt(i)
         enddo
         Whot1old=Whot1
         Whot2old=Whot2
         Wcoldold=Wcold
      endif
	
*     2.1.5 calcul de la nouvelle densite ionique
	
	do i = 1,ncell
	   dxt(i)=xt(i)-xt(i-1)
*	   ni1s2(i)=(qiSS(i-1)+qiSS(i))/dxt(i)/2.
	end do
	
	ni(0)=qiSS(0)/dxt(1)
	
*	do i = 1,ncell-1
*	   ni(i)=(dxt(i+1)*ni1s2(i)+dxt(i)*ni1s2(i+1))
*     *         /(dxt(i)+dxt(i+1))
*	end do
	
	do i = 1,ncell-1
	   ni(i)=2.d0*qiSS(i)/(dxt(i)+dxt(i+1))
	end do
	
*	ni(ncell)=1.5d0*ni1s2(ncell)-0.5d0*ni1s2(ncell-1)
* 	ni(ncell)=( (dxt(ncell-1)+2.d0*dxt(ncell))*ni1s2(ncell)
*     *             -dxt(ncell)*ni1s2(ncell-1) )
*     *		 /(dxt(ncell-1)+dxt(ncell))
*      if(ni(ncell).lt.0.d0)ni(ncell)=1.d-10
	
 	ni(ncell)=2.d0*qiSS(ncell)/dxt(ncell)
     * /(1.d0+2.d0*dxt(ncell)/dxt(ncell-1)-dxt(ncell-1)/dxt(ncell-2))

	do i = 1,ncell
	   ni1s2(i)=(ni(i-1)+ni(i))/2.d0
	end do

*	   2.2 - iteration pour le calcul de phi, Th, Tc

	do 10 iphi=1,iter

*	2.2.1 - construction des tableaux a,b,c, et f

	neh=n0hot*exp(-phi(0)/Th)
	nec=n0cold*exp(-phi(0)/Tc)
	b(0)=-2.d0/dxt(1)**2-neh/Th-nec/Tc
	c(0)=2.d0/dxt(1)**2
	f(0)=ni(0)-neh*(1.d0+phi(0)/Th)-nec*(1.d0+phi(0)/Tc)
	
	do i=1,ncell-1
	   neh=n0hot*exp(-phi(i)/Th)
	   nec=n0cold*exp(-phi(i)/Tc)
	   a(i)=2.d0/dxt(i)/(dxt(i)+dxt(i+1))
	   b(i)=-2.d0/(dxt(i)*dxt(i+1))-neh/Th-nec/Tc
	   c(i)=2.d0/dxt(i+1)/(dxt(i)+dxt(i+1))
	   f(i)=ni(i)-neh*(1.d0+phi(i)/Th)-nec*(1.d0+phi(i)/Tc)
	enddo
	
	neh=n0hot*exp(-phi(ncell)/Th)
	nec=n0cold*exp(-phi(ncell)/Tc)
	neold=neh+nec
	peold=neh*Th+nec*Tc
	a(ncell)=2.d0/dxt(ncell)**2
	b(ncell)=-a(ncell)
     *         -sqrt(2.d0/peold)*neold/dxt(ncell)
     *         -neh/Th-nec/Tc
        f(ncell)=ni(ncell)-neh*(1.d0+phi(ncell)/Th)
     *                  -nec*(1.d0+phi(ncell)/Tc)
     *         -sqrt(8.d0*peold)/dxt(ncell)
     *                   *(1.d0+0.5d0*neold*phi(ncell)/peold)

	do i=0,ncell
           if(ions_negatifs.and.xt(i).le.0.)then
             f(i)=f(i)+mixp
           endif
	enddo
	
     
*	2.2.2 - inversion de la matrice tridiagonale

	call gauss(nmax,ncell,a,b,c,f,bx,fx,phi)
	
      if(iphi.ne.iter.and.itime.eq.1) goto 10
      if(iphi.ne.iter.and.((multiphase.and.time.le.trise).or.
     *   ((.not.nb_cons).and.(.not.En_cons)))) goto 10


*     2.2.3 - calcul de la densite electronique dans le plasma

	do i=0,ncell
	   nhot(i)=n0hot*exp(-phi(i)/Th)
  	   ncold(i)=n0cold*exp(-phi(i)/Tc)
	   ne(i)=nhot(i)+ncold(i)
	   rho(i)=ni(i)-ne(i)
           if(ions_negatifs.and.xt(i).le.0.)then
             rho(i)= rho(i)+mixp
           endif
	enddo
      
*	2.2.4 - calcul de la pression electronique
*		  et des gradients de densite
*		  electronique et ionique au bord
	
	phot(ncell)=nhot(ncell)*Th
	pcold(ncell)=ncold(ncell)*Tc
	pe(ncell)=phot(ncell)+pcold(ncell)
	grade=(log(ne(ncell))-log(ne(ncell-1)))
     *           /(xt(ncell)-xt(ncell-1))
	gradi=(log(ni(ncell))-log(ni(ncell-1)))
     *           /(xt(ncell)-xt(ncell-1))
	
*     2.2.5 calcul de la densite electronique dans la gaine vide d'ions

      kvide=sqrt(n0hot/2.d0/Th)
      uphi=0.d0
      
      do i=ncell+1,ntotal
         if(i.eq.(ncell+1))then
            dxt(ncell+1)=0.
         else
            dxt(i)=0.0250*sqrt(Th/nhot(i-1))
         endif
         xt(i)=xt(i-1)+dxt(i)
         uphi=uphi+kvide*sqrt(1.d0+pcold(i-1)/phot(i-1))*dxt(i)
         phi(i)=phi(ncell)
     *          +2.d0*Th*log(1.d0+uphi*exp(-0.5d0*phi(ncell)/Th))
	   nhot(i)=n0hot*exp(-phi(i)/Th)
  	   ncold(i)=n0cold*exp(-phi(i)/Tc)
	   ne(i)=nhot(i)+ncold(i)
	   rho(i)=-ne(i)
	   phot(i)=nhot(i)*Th
	   pcold(i)=ncold(i)*Tc
         pe(i)=phot(i)+pcold(i)
         E(i)=sqrt(2.d0*pe(i))
      enddo
	
*     2.2.6 calcul du nombre total d'electrons et d'ions

      nti=0.d0
      do i=1,ncell
         nti=nti+ni1s2(i)*dxt(i)
      enddo

      nthot=0.d0
      ntcold=0.d0
      do i=1,ntotal
         nthot=nthot+0.5d0*(nhot(i-1)+nhot(i))*dxt(i)
         ntcold=ntcold+0.5d0*(ncold(i-1)+ncold(i))*dxt(i)
      enddo
      nthot=nthot+E(ntotal)
      nte=nthot+ntcold
      
      if(itime.eq.1) En_hot0=nthot*(e0/Tn0)*Th
      if(multiphase.and.time.le.trise) En_hot0=nthot*(e0/Tn0)*Th

      En_cold=ntcold*0.5d0*Tc

*     2.2.7 ajustement de n0hot pour conserver le nombre
*           total d'electrons chauds (alors le nombre d'electrons
*           froids est automatiquement conserve, car le potentiel
*           en 0 s'ajuste en consequence)

      if(nb_cons)then
         if(itime.eq.1)then
            nthot0=nthot
         else
            if(iphi.gt.1) n0hot=n0hot*(nthot0/nthot)
         endif
      endif
      
	if(En_cons.and.itime.eq.1)then
	   En_cold0=En_cold
	   write(*,*)
	   write(*,*)'En_hot0=', En_hot0
	   write(*,*)'En_cold0=', En_cold0
	   write(*,*)
	endif
	
*     2.2.7 calcul de la vitesse d'interface des cellules vides d'ions

      if(itime.gt.1)then
         do i=ncell+1,ntotal
            vint(i)=(xt(i)-xtold(i))/dt
         enddo
      endif

*     2.2.8 calcul des gradients de densite electronique

      do i=1,ntotal-1
         gradnh(i)=(nhot(i+1)-nhot(i-1))/(dxt(i)+dxt(i+1))
         gradnc(i)=(ncold(i+1)-ncold(i-1))/(dxt(i)+dxt(i+1))
      enddo
      
      gradnh(ncell)=gradnh(ncell-1)*nhot(ncell)/nhot(ncell-1)
      gradnh(ntotal)=gradnh(ntotal-1)*nhot(ntotal)/nhot(ntotal-1)
      gradnh(ncell+1)=gradnh(ncell)
      
      gradnc(ncell)=gradnc(ncell-1)
      gradnc(ntotal)=gradnc(ntotal-1)
      if(n0cold.gt.1.d-04)then
        if(ncold(ncell-1).gt.1.d-30)
     *   gradnc(ncell)=gradnc(ncell-1)*ncold(ncell)/ncold(ncell-1)
        if(ncold(ntotal-1).gt.1.d-30)
     *   gradnc(ntotal)=gradnc(ntotal-1)*ncold(ntotal)/ncold(ntotal-1)
      endif
      gradnc(ncell+1)=gradnc(ncell)

*     2.2.9 calcul de l'energie fournie par les electrons chauds [et
*           reajustement de Th (ancienne methode)]

      if(itime.gt.1.and.iphi.gt.1)then
         dWhot(0)=0.5d0*(phiold(0)+phi(0))
     *           *(nhot(0)-nhotold(0)-vint(0)*dt
     *           *0.5d0*(gradnh(0)+ghold(0)))
         Whot1=Whot1old+0.25d0*dWhot(0)*(dxtold(1)+dxt(1))
         SWhot(0)=0.25d0*dWhot(0)*(dxtold(1)+dxt(1))
         do i=1,ncell
            dWhot(i)=0.5d0*(phiold(i)+phi(i))
     *           *(nhot(i)-nhotold(i)
     *            -vint(i)*dt*0.5d0*(gradnh(i)+ghold(i)))
            Whot1=Whot1+0.25d0*dWhot(i)
     *           *(dxtold(i)+dxt(i)+dxtold(i+1)+dxt(i+1))
            SWhot(i)=SWhot(i-1)+0.25d0*dWhot(i)
     *           *(dxtold(i)+dxt(i)+dxtold(i+1)+dxt(i+1))
         enddo
         Whot2=Whot2old
         do i=ncell+1,ntotal-1
            dWhot(i)=0.5d0*(phiold(i)+phi(i))
     *           *(nhot(i)-nhotold(i)
     *            -vint(i)*dt*0.5d0*(gradnh(i)+ghold(i)))
            Whot2=Whot2+0.25d0*dWhot(i)
     *           *(dxtold(i)+dxt(i)+dxtold(i+1)+dxt(i+1))
            SWhot(i)=SWhot(i-1)+0.25d0*dWhot(i)
     *           *(dxtold(i)+dxt(i)+dxtold(i+1)+dxt(i+1))
         enddo
         dWhot(ntotal)=0.5d0*(phiold(ntotal)+phi(ntotal))
     *           *(nhot(ntotal)-nhotold(ntotal)-vint(ntotal)
     *           *dt*0.5d0*(gradnh(ntotal)+ghold(ntotal)))
         Whot2=Whot2+0.25d0*dWhot(ntotal)
     *           *(dxtold(ntotal)+dxt(ntotal))
         SWhot(ntotal)=SWhot(ntotal-1)+0.25d0*dWhot(ntotal)
     *           *(dxtold(ntotal)+dxt(ntotal))
         Whot=Whot1+Whot2

*         if(En_cons)then
*            if(iphi.ne.iter)then
*               Th=(5.d0*Th+2.d0*(En_hot0-Whot)/nthot)/6.d0
*               Th=2.d0*(En_hot0-Whot)/nthot
*            else
*               Th=2.d0*(En_hot0-Whot)/nthot
*            endif
*         endif

      endif
      
*     2.2.10 calcul de l'energie fournie par les electrons froids et
*           reajustement de Tc

      if(itime.gt.1.and.iphi.gt.1.and.
     *      (.not.multiphase.or.time.gt.trise))then
         Wcold=Wcoldold+0.125d0*(phiold(0)+phi(0))
     *           *(ncold(0)-ncoldold(0)-vint(0)*dt
     *           *0.5d0*(gradnc(0)+gcold(0)))*(dxtold(1)+dxt(1))
         do i=1,ntotal-1
            Wcold=Wcold+0.125d0*(phiold(i)+phi(i))
     *           *(ncold(i)-ncoldold(i)
     *            -vint(i)*dt*0.5d0*(gradnc(i)+gcold(i)))
     *           *(dxtold(i)+dxt(i)+dxtold(i+1)+dxt(i+1))
         enddo
         Wcold=Wcold+0.125d0*(phiold(ntotal)+phi(ntotal))
     *           *(ncold(ntotal)-ncoldold(ntotal)-vint(ntotal)
     *           *dt*0.5d0*(gradnc(ntotal)+gcold(ntotal)))
     *           *(dxtold(ntotal)+dxt(ntotal))
         if(En_cons.and.nb_cons)then
            if(iphi.ne.iter)then
               Tc=(3.d0*Tc+2.d0*(En_cold0-Wcold)/ntcold)/4.d0
*               Tc=2.d0*(En_cold0-Wcold)/ntcold
            else
               Tc=(1.d0*Tc+2.d0*(En_cold0-Wcold)/ntcold)/2.d0
*               Tc=2.d0*(En_cold0-Wcold)/ntcold
            endif
         endif
      endif
      
*     2.2.11 calcul du champ electrique dans toute la detente

	do i = 1,ncell
	   E1s2(i)=(phi(i)-phi(i-1))/dxt(i)
	end do
		
	do i = 1,ncell-1
	   E(i)=(dxt(i+1)*E1s2(i)+dxt(i)*E1s2(i+1))/(dxt(i)+dxt(i+1))
	end do
	
*	do i = 1,ncell-1
*	   E(i)=E1s2(i)+qiSS(i)/(1.d0+prog)
*     *        -(3.d0*ne(i)+ne(i-1))/8.d0*dxt(i)
*	end do
		
 	E(ncell)=sqrt(2.d0*pe(ncell))
 	
	do i = 1,ncell
	   if(xttest(i).ge.xt(itest(i)))then
	     do j=itest(i),ntotal-1
	       if(xttest(i).le.xt(j+1))then
	         itest(i)=j
	         goto 301
	       endif
	     enddo
	     itest(i)=ntotal-1
301      continue
	   else
	     do j=itest(i)-1,0,-1
	       if(xttest(i).ge.xt(j))then
	         itest(i)=j
	         goto 302
	       endif
	     enddo
	     itest(i)=0
302      continue
         endif
	   
	   if(itest(i).ne.ncell) then
	     Etest(i)=(E(itest(i))*(xt(itest(i)+1)-xttest(i))
     *             +E(itest(i)+1)*(xttest(i)-xt(itest(i))))
     *            /dxt(itest(i)+1)
         else
	     Etest(i)=E(ncell)
	   endif
	end do
 	

*     2.2.12 ajustement de la vitesse aux temps 'entiers'
	
	if(itime.gt.1) then
	   do i = 1,ncell
	      v(i)=vint(i)+dt*(3.d0*E(i)+Eold(i))*charge(i)/8.d0
	   end do	
	endif

	if(itime.gt.1) then
	   do i = 1,ncell
	      vtest(i)=vtestint(i)
     *               +dt*(3.d0*Etest(i)+Etestold(i))*Ztest/8.d0
	   end do	
	endif

*     2.2.13 calcul de l'energie cinetique des ions

      En_ion=0.
	do i = 1,ncell-1
	   En_ion=En_ion+qiSS(i)*v(i)**2/2.d0/charge(i)
	end do	
	if(profil.eq.'step') En_ion=
     *              En_ion+qiSS(ncell)*v(ncell)**2/4.d0/charge(i)
	if(profil.ne.'step') En_ion=
     *              En_ion+qiSS(ncell)*v(ncell)**2/2.d0/charge(i)

*     2.2.14 calcul de l'energie electrostatique

      En_elec=0.
	do i = 1,ntotal-1
	   En_elec=En_elec+0.25d0*(E(i)**2)*(dxt(i)+dxt(i+1))
	end do
	En_elec=En_elec+0.25d0*(E(ntotal)**2)*dxt(ntotal)
	En_elec=En_elec+E(ntotal)*Th
			
*     2.2.15 calcul de l'energie totale

	if(itime.eq.1.or.(multiphase.and.time.le.trise))then
*	   En_totale=En_ion+En_elec+0.5d0*(nthot*Th+ntcold*Tc)
	   En_totale=En_ion+En_elec+0.5d0*ntcold*Tc+En_hot0
	endif

*	2.2.16 reajustement de Th 
*	       (remplace l'ajustement precedent)

      if(En_cons.and.itime.gt.1.and.iphi.gt.1.and.
     *    (.not.multiphase.or.time.gt.trise))then

	   enew=(En_totale-En_ion-En_elec-0.5d0*ntcold*Tc)
     *         *Tn0/nthot/Thmax
         if(EOS) then
	     Thnew=2.*Thmax*enew*(1.+9.*enew/2.+3.*enew*enew/2.)
     *                /(1.+6.*enew+3.*enew*enew)/Tn0
         else
	     Thnew=2.*Thmax*enew/Tn0
         endif
         if(iphi.ne.iter)then
	     Th=(iter2*Th+Thnew)/(iter2+1)
         else
	     Th=(iter3*Th+Thnew)/(iter3+1)
         endif
	endif
	
10    continue


*         2.3 champs 'theoriques'

*	2.3.1 normalisation au champ theorique du modele self-similaire

 	if(itime.eq.1)then
 	   Enorm=0.d0
 	else
 	   Ess=1.d0*sqrt(Th)/time
 	   Enorm=E(ncell)/Ess
 	endif

*	2.3.2 comparaison au champ 'theorique' suppose 'fite' le
*	      champ au bord a tout instant

 	Ebord0=sqrt(2.d0*n0hot)/exp(0.5d0)
 	Ebordth=2.d0*Ebord0*sqrt(Th)/sqrt(4.d0+(Ebord0*time)**2)

	
*         2.4 test d'arret sur le nombre d'iterations
	
	if(itmax.eq.1) goto 101


*         2.5 sorties de fonctions dependant du temps

	vmax=0.d0
	ivmax=ncell
	do i=1,ncell
	   if(v(i).gt.vmax) then
	     vmax=v(i)
	     ivmax=i
	   endif
	enddo
	vfinal=v(ivmax)+E(ivmax)*charge(ivmax)*time

*	vfinal=v(ncell)+E(ncell)*time
	lDebye=sqrt(Th/ne(ncell))
	lgrade=-1./grade
	if(itime.eq.1)then
	   lgradi=0.
	else
	   lgradi=-1./gradi
	endif

	nstep=(tmax/dt)
	nsort=1
	if(nstep.gt.2000)nsort=nstep/2000

	if(itime-1.ne.((itime-1)/nsort)*nsort.and..not.laststep)
     *   go to 102
	write(9,
     *'(f8.2,17(a,f10.5))')
     *   time,char(9),nti,char(9),nthot,char(9),ntcold
     *   ,char(9),nte,char(9),n0hot,char(9),n0cold
     *   ,char(9),En_ion,char(9),Whot1
     *   ,char(9),Whot2,char(9),Whot,char(9),Wcold,char(9),Th
     *   ,char(9),Tc,char(9),En_elec,char(9),En_elec+En_ion-Whot
     *   ,char(9),vmax,char(9),vfinal

	write(10,
     *'(f8.2,a,f12.5,a,f7.4,a,f8.4,a,f7.5
     *  ,a,f10.7,a,f10.7,3(a,f8.6),a,f10.5,a,f10.5,a,f10.5
     *  ,a,i6,a,i6)')
     *   time,char(9),xt(ivmax),char(9),v(ivmax),char(9),.5*v(ivmax)**2
     *   ,char(9),E(ivmax)
     *   ,char(9),ne(ivmax),char(9),ni(ivmax)
     *   ,char(9),ni(0)
     *   ,char(9),nhot(0),char(9),ncold(0)
     *   ,char(9),lDebye
     *   ,char(9),lgrade,char(9),lgradi
     *   ,char(9),ivmax,char(9),idebut(ivmax)
	
	Rs=rgauss/ni(0)
	xi=(xt(ivmax)+lmax)/Rs
	tsR0=time/rgauss
	RsR0=Rs/rgauss
	Eq42=2.d0*xi*exp(xi**2/2.d0)/rgauss/sqrt(RsR0)

	write(14,
     *'(f8.2,a,f7.3,a,f7.4,a,f8.3,a,f8.4,a,f7.4)')
     *   time,char(9),tsR0,char(9),xi,char(9),RsR0,char(9),Eq42

102	continue

*         2.6 test d'arret sur le temps (ici s'arrete normalement
*		  un calcul standard avec itmax.gt.1 et tmax.gt.0 apres
*		  un nombre fini d'iterations)
	
	if(time.ge.(tmax-1.d-5)) goto 101
	if(itime.eq.itmax) goto 101


*         2.7 incrementation du temps
	
	dtold=dt
	dt=dti
	if(VTT.and..not.multicouche.and.
     *     (.not.multiphase.or.time.gt.trise)) dt=dti/sqrt(ni(0))
	if(VTT.and.multicouche.and.
     *     (.not.multiphase.or.time.gt.trise)) 
     *       dt=dti/sqrt(5.d0*qiSS(0)/(xt(10)-xt(0)))
*	if(idico.eq.1)dt=dti/100.d0
        if(multiphase.and.dt.gt.trise/100.d0.and.time.lt.trise)
     *      dt=trise/100.d0
        if(multiphase.and.time.lt.trise.and.(time+dt).gt.trise)
     *      dt=trise-time
	if((time+dt).gt.(tmax-1.d-5)) then
	   dt=tmax-time
	   laststep=.true.
	endif
	time=time+dt

	
*         2.8 modification de la vitesse aux temps intermediaires
*             'demi-entiers'
	
	if(itime.eq.1) then

	   do i = 1,ncell
	      vint(i)=v(i)+0.5d0*dt*E(i)*charge(i)
	   end do
	
	   do i = 1,ncell
	      vtestint(i)=vtest(i)+0.5d0*dt*Etest(i)*Ztest
	   end do
	
	else

*	   do i = 1,ncell
*	      vint(i)=vint(i)+0.5*(dtold+dt)*E(i)*charge(i)
*	   end do

	   do i = 1,ncell
	      vtestint(i)=vtestint(i)+0.5*(dtold+dt)*Etest(i)*Ztest
	   end do


*            2.8.1 - construction des tableaux a,b,c, et f pour 
*			  le calcul implicite de la vitesse avec viscosite

	   delt=0.5*(dtold+dt)
	   b(0)=1.d0/delt
	   c(0)=0.d0
	   f(0)=0.d0
	
	   do i=1,ncell-1
	      a(i)=-2.d0*nu/dxt(i)/(dxt(i)+dxt(i+1))
	      b(i)=1.d0/delt+2.d0*nu/(dxt(i)*dxt(i+1))
	      c(i)=-2.d0*nu/dxt(i+1)/(dxt(i)+dxt(i+1))
	      f(i)=E(i)*charge(i)+vint(i)/delt
	   enddo
	
	   a(ncell)=0.d0
	   b(ncell)=1.d0/delt
         f(ncell)=E(ncell)*charge(ncell)+vint(ncell)/delt
     
*	      2.8.2 - inversion de la matrice tridiagonale

	   call gauss(nmax,ncell,a,b,c,f,bx,fx,vint)
	
	endif
	

*         2.9 modification de la position
	
	do i = 1,ncell
	   xt(i)=xt(i)+dt*vint(i)
           if(xt(i).lt.x0(0))then
              xt(i)=2.d0*x0(0)-xt(i)
              vint(i)=-vint(i)
           endif
	end do
	
	do i = 1,ncell
	   xttest(i)=xttest(i)+dt*vtestint(i)
           if(xttest(i).lt.x0(0))then
              xttest(i)=2.d0*x0(0)-xttest(i)
              vtestint(i)=-vtestint(i)
           endif
	end do
	

*         2.10 rearrangement des numeros des ions

	
	do j=2,ncell
	   xxt=xt(j)
	   xx0=x0(j)
	   xxtold=xtold(j)
	   vv=v(j)
	   vvint=vint(j)
	   iidebut=idebut(j)
	   qqiSS=qiSS(j)
	   nnhot=nhot(j)
	   nncold=ncold(j)
	   pphi=phi(j)
	   EE=E(j)
	   ggradnh=gradnh(j)
	   ggradnc=gradnc(j)
	   ddxt=dxt(j)
          ccharge=charge(j)
	   do i=j-1,1,-1
	      if(xt(i).le.xxt) goto 95
	      xt(i+1)=xt(i)
	      x0(i+1)=x0(i)
	      xtold(i+1)=xtold(i)
	      v(i+1)=v(i)
	      vint(i+1)=vint(i)
	      idebut(i+1)=idebut(i)
		qiSS(i+1)=qiSS(i)
		nhot(i+1)=nhot(i)
		ncold(i+1)=ncold(i)
		phi(i+1)=phi(i)
		E(i+1)=E(i)
		gradnh(i+1)=gradnh(i)
		gradnc(i+1)=gradnc(i)
		dxt(i+1)=dxt(i)
              charge(i+1)=charge(i)
	   enddo
	   i=0
 95	   xt(i+1)=xxt
	   x0(i+1)=xx0
	   xtold(i+1)=xxtold
	   v(i+1)=vv
	   vint(i+1)=vvint
	   idebut(i+1)=iidebut
	   qiSS(i+1)=qqiSS
	   nhot(i+1)=nnhot
	   ncold(i+1)=nncold
	   phi(i+1)=pphi
	   E(i+1)=EE
	   gradnh(i+1)=ggradnh
	   gradnc(i+1)=ggradnc
	   dxt(i+1)=ddxt
          charge(i+1)=ccharge
	enddo
	
	do j=1,ncell
	   irang(idebut(j))=j
	enddo
	
100   continue 
101	continue
	close(9)
	close(10)
	close(14)

	write(*,*) 'OK etape2'

	
*               3 - calcul de la densite theorique (modele
*		    self-similaire)
	

      if(itime.eq.1)then
         do i=ncell+1,ntotal
            x0(i)=xt(i)
         enddo
      endif
	
	if(itmax.eq.1.or.time.le.1.d-04) goto 201

	do i=0,ntotal
	   xi=xt(i)/sqrt(Th)/time
	   if(xi.lt.-1.d0) then
	      nss(i)=n0hot
	   else
	      nss(i)=n0hot*exp(-(xi+1.d0))
	   endif
	enddo
	
	do i=0,ntotal
	   difth(i)=(ni(i)-nss(i))/nss(i)
	enddo

201	continue

	write(*,*) 'OK etape3'


*			4 - sorties graphiques

*		4.1 sorties des profils 

	open(11,file='profil',status='unknown')
		
	write(11,'(a,15(a,a))')'x0',char(9),'xt',char(9),'v'
     *      ,char(9),'phi',char(9),'E',char(9),'ni'
     *      ,char(9),'ne',char(9),'rho'
     *      ,char(9),'nhot',char(9),'ncold'
     *      ,char(9),'nss',char(9),'difnor',char(9),'dWhot'
     *      ,char(9),'SWhot',char(9),'charge'
	
	do i=0,ncell,iline
       write(11,'(f11.4,a,f12.4,a,f7.4,2(a,f9.6),6(a,f12.8),a,f9.6
     *              ,2(a,e11.4),a,f5.2)') 
     *	   x0(i),char(9),xt(i),char(9),v(i)
     *         ,char(9),phi(i),char(9),E(i),char(9),ni(i)
     *         ,char(9),ne(i),char(9),rho(i)
     *         ,char(9),nhot(i),char(9),ncold(i)
     *         ,char(9),nss(i),char(9),rho(i)/ni(i)
     *         ,char(9),dWhot(i),char(9),SWhot(i),char(9),charge(i)
	enddo
	
	do i=ncell+1,ntotal
	   write(11,'(f11.4,a,f12.4,a,f7.4,2(a,f9.6),6(a,f12.8),a,f9.6
     *              ,2(a,e11.4))') 
     *	   x0(i),char(9),xt(i),char(9),v(i)
     *         ,char(9),phi(i),char(9),E(i),char(9),ni(i)
     *         ,char(9),ne(i),char(9),rho(i)
     *         ,char(9),nhot(i),char(9),ncold(i)
     *         ,char(9),nss(i),char(9),ni(i)
     *         ,char(9),dWhot(i),char(9),SWhot(i)
	enddo
	
	close(11)
	
	open(11,file='profil_test',status='unknown')
		
	write(11,'(a,4(a,a))')'x0test',char(9),'xttest',char(9),'vtest'
     *             ,char(9),'itest',char(9),'E'
	
	do i=0,ncell,iline
	   write(11,'(f11.4,a,f12.4,a,f7.4,a,i4,a,f9.6)') 
     *	   x0test(i),char(9),xttest(i),char(9),vtest(i),char(9),itest(i)
     *     ,char(9),Etest(i)
	enddo
		
	close(11)
	
	open(11,file='profil_lisse',status='unknown')
	
	write(11,'(a,4(a,a))')'x0',char(9),'xt',char(9),'xt+lmax'
     *      ,char(9),'nilisse',char(9),'nelisse'
     
      do i=0,ncell,iline
         if(i.le.nlisse.or.i.ge.ncell-nlisse) then
            nilisse(i)=ni(i)
            nelisse(i)=ne(i)
         else
            qq=qiSS(i)
            div=dxt(i)+dxt(i+1)
            nelisse(i)=ne(i)
            do j=1,nlisse
               qq=qq+qiSS(i-j)+qiSS(i+j)
               div=div+dxt(i-j)+dxt(i-j+1)+dxt(i+j)+dxt(i+j+1)
               nelisse(i)=nelisse(i)+ne(i-j)+ne(i+j)
            enddo
            nilisse(i)=2.d0*qq/div
            nelisse(i)=nelisse(i)/(2*nlisse+1)
         endif
	   write(11,'(f11.4,a,f12.4,a,f12.4,2(a,f12.8))') 
     *	   x0(i),char(9),xt(i),char(9),xt(i)+lmax
     *         ,char(9),nilisse(i),char(9),nelisse(i)
      enddo
      
	do i=ncell+1,ntotal
	   write(11,'(f11.4,a,f12.4,a,f12.4,2(a,f12.8))') 
     *	   x0(i),char(9),xt(i),char(9),xt(i)+lmax
     *         ,char(9),ni(i),char(9),ne(i)
      enddo
	
      close(11)
	
	open(11,file='profil2',status='unknown')
		
	write(11,'(a,13(a,a))')'x0',char(9),'xt',char(9),'v'
     *      ,char(9),'phi',char(9),'E',char(9),'ni'
     *      ,char(9),'ne',char(9),'rho'
     *      ,char(9),'nhot',char(9),'ncold'
     *      ,char(9),'nss',char(9),'difnor',char(9),'dWhot'
     *      ,char(9),'SWhot'
		
	do i = 1,ncell
	   dxt(i)=xt(i)-xt(irang(idebut(i)-1))
*	   dxt(i)=abs(dxt(i))
	end do
	
	do i = 1,ncell-1
	   if(idebut(i).ne.ncell)then
	      ni(i)=2.d0*qiSS(i)/(dxt(i)+dxt(irang(idebut(i)+1)))
	   else
	      ni(i)=2.d0*qiSS(i)/dxt(i)/(1.+prog)
	   endif
	   ni(i)=abs(ni(i))
	end do
	
	do i=0,ncell,iline
	   write(11,'(f11.4,a,f12.4,a,f7.4,2(a,f9.6),6(a,f12.8),a,f9.6
     *              ,2(a,e11.4))') 
     *	   x0(i),char(9),xt(i),char(9),v(i)
     *         ,char(9),phi(i),char(9),E(i),char(9),ni(i)
     *         ,char(9),ne(i),char(9),rho(i)
     *         ,char(9),nhot(i),char(9),ncold(i)
     *         ,char(9),nss(i),char(9),rho(i)/ni(i)
     *         ,char(9),dWhot(i),char(9),SWhot(i)
	enddo
	
	do i=ncell+1,ntotal
	   write(11,'(f11.4,a,f12.4,a,f7.4,2(a,f9.6),6(a,f12.8),a,f9.6
     *              ,2(a,e11.4))') 
     *	   x0(i),char(9),xt(i),char(9),v(i)
     *         ,char(9),phi(i),char(9),E(i),char(9),ni(i)
     *         ,char(9),ne(i),char(9),rho(i)
     *         ,char(9),nhot(i),char(9),ncold(i)
     *         ,char(9),nss(i),char(9),ni(i)
     *         ,char(9),dWhot(i),char(9),SWhot(i)
	enddo
	
	close(11)
	
		
*		4.2 spectres en vitesse et en energie
	
	do i=1,ncell
	   vmoy(i)=0.5*(v(i-1)+v(i))
	   Emoy(i)=0.25*(v(i-1)**2+v(i)**2)
	   dndv(i)=(niSS(i-1)+niSS(i))*dx0(i)/(v(i)-v(i-1))/2.d0
	   dndv(i)=abs(dndv(i))
	   dndE(i)=(niSS(i-1)+niSS(i))*dx0(i)/(v(i)**2-v(i-1)**2)
	   dndE(i)=abs(dndE(i))
	enddo
	
	open(12,file='spectres', status='unknown')
	
	write(12,'(a,3(a,a))')'vmoy', char(9), 'Emoy'
     *      , char(9), 'dndv', char(9), 'dndE'
      do i=1,ncell,iline
         if(dndE(i).lt.1.d03)then
	      write(12,'(f7.4,a,f8.4,a,f10.6,a,f10.6)')
     *            vmoy(i), char(9), Emoy(i)
     *            ,char(9), dndv(i), char(9), dndE(i)
         endif
      enddo
      
      close(12)
	
	write(*,*) 
	write(*,*) time
	write(*,*) 
	stop
	end      

c***5***10***15***20***25***30***35***40***45***50***55***60***65***70**
c***5***10***15***20***25***30***35***40***45***50***55***60***65***70**


      Subroutine Gauss(nmax,ncell,A,B,C,F,bx,fx,phi)

c   tridiagonal matrix reduction   M*PHI = F
c
c      B0 C0	           PHI0     F0
c      A1 B1 C1              PHI1     F1
c          .  .  .        *   .    =  .
c            Ai Bi Ci        PHIi     Fi
c                .  .  .      .       .
c                  An Bn     PHIn     Fn
c

      double precision phi (0:nmax)
	double precision a (1:nmax),  b(0:nmax), c(0:nmax-1), f(0:nmax)
     *                ,bx(0:nmax), fx(0:nmax)

C... ascending recursion for Gauss elimination into bidiagonal system
C... first elements in j=0 and C in general are unchanged

      j=0
      Bx(j) = B(j)
      Fx(j) = F(j)
	
      do j=1,ncell
         Bx(j) = B(j) - A(j)*C(j-1) /Bx(j-1) 
         Fx(j) = F(j) - A(j)*Fx(j-1)/Bx(j-1) 
      end do
	
C... descending recursion for solution

      j = ncell
      phi(j) = Fx(j)/Bx(j)
	
      do j=ncell-1,0,-1
         phi(j) = ( Fx(j)-C(j)*phi(j+1) )/Bx(j)
      end do

      Return
      End
	
c***5***10***15***20***25***30***35***40***45***50***55***60***65***70**
c***5***10***15***20***25***30***35***40***45***50***55***60***65***70**

      function Thot(multiphase,trise,time)
      
      implicit double precision (a-h,o-z)
      logical multiphase

      if(multiphase.and.time.le.trise) then
	   Thot=0.01d0+0.99d0*time/trise
	else
	   Thot=1.0d0
	endif
      
      Return
      End
	
c***5***10***15***20***25***30***35***40***45***50***55***60***65***70**
c***5***10***15***20***25***30***35***40***45***50***55***60***65***70**

      function Tcold(time)
      
      implicit double precision (a-h,o-z)
      if(time.le.9.9d0) then
	   Tcold=1.d0
	else
	   Tcold=1.d0
	endif
      
      Return
      End
	
c***5***10***15***20***25***30***35***40***45***50***55***60***65***70**
c***5***10***15***20***25***30***35***40***45***50***55***60***65***70**

