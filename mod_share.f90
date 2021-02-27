module mod_share

     use iso_fortran_env, only: int32, real32,real64

     implicit none

     integer(int32), parameter :: nmax=180000, njmax=200
     integer(int32)  :: ntotal, nsort, iline, nstep, iline_vide, itime,ivmax
     integer(kind=4) :: idebut(0:nmax),irang(0:nmax),itest(0:nmax), m,djout
     integer(int32)  :: i, j, idico, iter, iphi, iidebut
     real(real64) :: length,ni0,n0cold,n0hot,neh,nec,neold,kvide
     real(real64) :: nti,nthot,ntcold,nte,lmax
     real(real64) :: nthot0nhm3,nhm2,nhm1,LSS,nu,nnhot,nncold, nhm3, Tcm3, Thm3
     real(real64) :: lDebye,lgrade,lgradi,lay1,lay2,mix,mixp
     real(real64) :: x0 (0:nmax), xt (0:nmax), xtold(0:nmax)
     real(real64) :: v  (0:nmax),  vint (0:nmax)
     real(real64) :: E  (0:nmax),  Eold (0:nmax)
     real(real64) :: phi (0:nmax),  ni (0:nmax)
     real(real64) :: nilisse (0:nmax),  nelisse (0:nmax)
     real(real64) :: phiold (0:nmax)
     real(real64) :: ne (0:nmax),  rho (0:nmax)
     real(real64) :: nhot (0:nmax),  ncold (0:nmax)
     real(real64) :: nhotold (0:nmax),  ncoldold (0:nmax)
     real(real64) :: gradnh (0:nmax),  gradnc (0:nmax)
     real(real64) :: ghold (0:nmax),  gcold (0:nmax)
     real(real64) :: phot (0:nmax),  pcold (0:nmax)
     real(real64) :: pe (0:nmax)
     real(real64) :: nss (0:nmax),  difth (0:nmax)
     real(real64) :: dWhot (0:nmax), SWhot(0:nmax)
     real(real64) :: vmoy (1:nmax), Emoy (1:nmax), dndv (1:nmax), dndE (1:nmax)
     real(real64) :: dx0  (1:nmax),  dxt (1:nmax), dxtold(1:nmax)
     real(real64) :: ni1s2 (1:nmax), E1s2 (1:nmax), niSS (0:nmax), qiSS(0:nmax), charge(0:nmax)

     real(real64) :: a (1:nmax), b(0:nmax), c(0:nmax-1), f(0:nmax), bx(0:nmax), fx(0:nmax)
    
     logical lfini,nb_cons,En_cons,laststep
     integer(int32) :: ncell, nvide,Nbe, itmax, jcoldmax
     real(real64)   :: ve_max, prog, iter0, iter1, iter2, iter3, dti,tmax
     real(real64)   :: Tcmax,Thmax, T_MeV, nLSS, charge2, p_negatif, trise
     real(real64)    :: omegadt, Th, Tc, dt, cs, cs2, time, cs2old, rgauss
     real(real64) :: PI, e0, Tnorm, Tn0, alpha, whot
     real(real64) :: Thm1, Thm2, Whot2, Whot1, Tcm1,Tcm2,Wcold, Tcold, Thot, Thold,Tcoold, En_ion, En_totale, En_elec
     real(real64) :: enew, thnew
     real(real64) :: whot2old, whot1old, wcoldold, peold, grade, gradi, uphi, en_hot0, En_cold, En_cold0,nthot0

    real(real64)  :: xxtold, xxt, xx0, vvint, vv, qqiss, pphi, ggradnh, ggradnc, EE, ccharge, xi
    real(real64)  :: Ess,Enorm, Ebord0, Ebordth, vfinal, vmax, dtold, delt, ddxt
    character*10     profil
    integer       :: fileunit



end module mod_share

