      subroutine fwleaps 
     1    (nv,it,kx,nf,no,ib,fd,mb,rt,nd,nc,iw,nw,rw,nr,tl,s2,ne,iv,
     2     RSS,Subss,nret)

      double precision rt,rw,RSS,Subss 
      dimension rt(nd,nc),iw(nw),rw(nr),RSS(nret),Subss(nret)

c ******************************************************************** *
c Adaptation of Furnival and Wilson leaps and bounds algorithm, adapted
c by Ian Painter, 2004. Modification consists of commenting out all
c write statements, and returning results in the vectors RSS and Subss.
c Code kindly provided by Sanford Weisberg.
c                                                                      *
c                   Regressions by leaps and bounds                    *
c                     G.M.Furnival and R.W.Wilson                      *
c                           Version 8/10/81                            *
c                                                                      *
c call screen(nv,it,kx,nf,no,ib,fd,mb,rt,nd,nc,iw,nw,rw,nr,tl,s2,ne,iv)*
c                                                                      *
c  nv = number of rows and columns in input matrix (nv.le.m where m    *
c       is the number of bits in the fractional part of the real       *
c       numbers employed in the computations)                          *
c  it = type of input matrix (0=products about 0, 1=products about     *
c       means, 2=correlation matrix)                                   *
c  kx = number of independent variables (kx.lt.nv)                     *
c  nf = number of indepent variables to be forced i.e. included in     *
c       every subset (nf.lt.kx)                                        *
c  no = number of observations (no.gt.kx when it=0, no.gt kx+1         *
c       otherwise)                                                     *
c  ib = selection criterion code (1=r**2, 2=adjusted r**2, 3=mallows   *
c       c(p) with frane's variable penalty)                            *
c  fd = variable penalty or f-to-delete for c(p). a value of 2.0       *
c       gives mallows original c(p). (a value of fd .le.0 is replaced  *
c       by 2.0)                                                        *
c  mb = number of best regressions desired-for each size subset when   *
c       ib=1, in total when ib.gt.1                                    *
c  rt = real two-dimensional array. first nv columns contain input     *
c       matrix of products or correlations. rest of array is working   *
c       storage.                                                       *
c  nd = first dimension of rt(nd.ge.nv) 
c  nc = second dimension of rt (nc.ge.3*nd)                            *
c  iw = integer one-dimensional array. first kx+1 elements contain     *
c       subscripts of variables to be used in problem with subscripts  *
c       of forced variables in first nf elements, subscripts of free   *
c       variables in next kx-nf elements and subscript of dependent    *
c       variable in element kx+1. rest of array is working storage.    *
c  nw = dimension of iw (nw.ge.4*nd)                                   *
c  rw = real one dimensional array. when it=1,first nv elements must   *
c       contain vector of means. rest of array is working storage.     *
c  nr = dimension of rw (nr.ge.2*mb*kx+7*nd if ib=1, nr.ge.2*mb+7*nd   *
c       if ib.gt.1)                                                    *
c  tl = minimum allowable tolerance (tl.lt.kx*kx/2**(m/2) or.ge.1      *
c       reset to kx*kx/2**(m/2))                                       *
c  s2 = estimated variance of residuals for c(p) (value.le.0 replaced  *
c       by mean square residual from full model)                       *
c  ne = error code (1=input parameter out of bounds, 2=invalid list    *
c       of subscripts in iw, 3=zero or negative diagonal in input      *
c       matrix rt, 4=set of forced variables nearly collinear,         *
c       5=full model nearly collinear-s2 must be supplied for c(p),    *
c       6=r**2 so close to 1 (r**2.gt.1-2**(-m/2)) that selection of   *
c       best subsets is uncertain.                                     *
c  iv = version of algorithm (0=partial reordering, 2=full reordering) *
c 
c
c             RSS, Subss and nret added by ISP 5-2004
c RSS = real one dimensional array that returns the R**2 values of the *
c       selected subsets.                                              *
c Subss = real one dimensional array that returns the subsets selected.*
c nret = length of RSS and Subss, equal to ib (see below)              *
c                                                                      *
c ******************************************************************** *
      ne=0
      if(ib.eq.1) ns=mb*kx
      if(ib.ne.1) ns=mb 
      if(
     1   kx+min0(1,it).lt.no.and.nc.ge.3*nd.and.nw.ge.4*nd.and.
     2   nr.ge.2*ns+7*nd.and.ib.ge.1.and.ib.le.3.and.it.ge.0. 
     3   and.it.le.2.and.mb.gt.0.and.iv.ge.0.and.iv.le.1) 
     4   go to 10
      ne=1
      go to 11
   10 nrm=7*nd+ns+1 
      call forwd(rt(1,1),rt(1,nd+1),rt(1,2*nd+1),rw(1),rw(nd+1),
     1   rw(2*nd+1),rw(3*nd+1),rw(4*nd+1),rw(5*nd+1),rw(6*nd+1),
     2   rw(7*nd+1),rw(nrm),iw(1),iw(nd+1),iw(2*nd+1),iw(3*nd+1),ns,
     3   nd,kx,no,it,nf,ib,mb,s2,tl,nv,fd,ne,iv,RSS(1),Subss(1))
c   11 if(ne.ne.0) write(*,12) ne
c   12 format(' screen error code is',i3)
   11  return
      end 
c                      bounds from forward pivots 
      subroutine forwd(rr,xi,xn,xm,dd,d,dt,co,sc,sq,cl,rm,inn,ipp,ind,
     1   int,ns,nd,kx,no,it,nf,ib,mb,s2,tl,nv,fd,ne,iv,rss,sbs)
      double precision rr,xi,xm,xn,dd,d,dt,co,sc,sq,cl,rm,ci,cn,zip
      double precision tol,tal,tul,rt,rs,sig,bit,df,dint,pen,wt,bnd
      double precision wst,rss,sbs,rsret,tn,sel,beta
      double precision tel,f,yint
      real dx,di
      dimension rr(nd,nd),xi(nd,nd),xn(nd,nd),xm(nd),dd(nd),d(nd),
     1   dt(nd),co(nd),sc(nd),sq(nd),cl(ns),rm(ns),inn(nd),ind(nd), 
     2   int(nd),ipp(nd),rss(ns),sbs(ns)
c    multreg changes,next 8 lines 
      common/tab/lent,labl(34)
      common/qscren/iscrn 
      bit=2.0 
      do 10 nbit=1,1000 
      if(bit.eq.bit+1.0) go to 11 
      bit=bit+bit 
   10 continue
   11 if(nv.le.nbit) go to 13 
      ne=1
      return
   13 continue
      do 14 l=1,nv
      ind(l)=l
      int(l)=0
      co(l)=2.0**(nv-l) 
   14 continue
      zip=0.0 
      call copy(rr,rr,ind,ind,nv,0,nv,nd,sc,cn,zip,0) 
c                          check subscript list 
      kz=kx+1 
      do 17 l=1,kz
      m=inn(l)
      ind(l)=m
      if(m.ge.1.and.m.le.nv.and.int(m).eq.0) go to 15 
      ne=2
      return
   15    if(rr(m,m).gt.zip) go to 16
      ne=3
      return
   16    tn=rr(m,m) 
c                         scale by powers of 16 
      sc(m)=16.0**ifix(alog(real(tn))/2.0/alog(16.0)) 
      sq(m)=rr(m,m)/(sc(m)*sc(m)) 
      int(m)=1
   17 continue
      dx=real(zip)
      pen=fd
      if(pen.le.zip.and.ib.eq.3) pen=2.0
      df=no-min0(it,1)
      sel=1.0 
      ky=ind(kz)
c                            set tolerances 
      tul=sq(ky)/2.0**(nbit/2)
      tol=float(kx*kx)/2.0**(nbit/2)
      tal=tl
      if(tal.lt.tol.or.tal.ge.1.0) tal=tol
      call copy(xn,rr,ind,ind,kz,1,kz,nd,sc,cn,zip,1) 
c          forced variables and mean square residual for c(p) 
      sig=s2/(sc(ky)*sc(ky))
      do 20 l=1,kx
      if(l.eq.nf+1) call copy(xi,xn,ind,ind,kz,1,kz,nd,sc,ci,cn,0)
      if(l.gt.nf.and.(ib.ne.3.or.sig.gt.zip)) go to 21
      n=ind(l)
      ipp(l)=l+1
      call test(wt,n,xn,int,d,dd,sq,dt,l-1,l,nd)
      if(wt.ge.tal.or.l.gt.nf) go to 18 
      ne=4
      return
   18    if(wt.ge.tol) go to 19 
      ne=5
      return
   19    call pivot(n,l-1,l+1,
     1      kx,xn,int,ind,cn,co(n),dt,1,nd,nn,ky,l+1,kx)
   20 continue
      sig=xn(ky,ky)/(df-float(kx))
c                load transform storage with worst value
   21 call trans(wst,mr,sq(ky)+sq(ky),df,kx,sig+sig,pen,ib,mb)
      do 22 l=1,ns
      cl(l)=zip 
      rm(l)=wst 
   22 continue
      mn=nf+1 
      ipp(mn)=mn
c                           forward algorithm 
   23 ls=ipp(mn)
      mh=mn 
      mt=mn 
      if(ls.lt.kx-2) call copy
     1      (xn,xi,int,ind,mn-1,ls,mn+kz-ls,nd,sc,cn,ci,0)
      do 28 lb=ls,kx
      mn=mh+lb-ls 
      ipp(mn)=lb+1
c                             select pivot
      nb=lb 
      do 24 l=lb,kx 
      j=ind(l)
      if(mn.eq.1) dd(l)=-xi(j,ky)*xi(j,ky)/xi(j,j)
      if(mn.gt.1) dd(l)=xi(j,j)/sq(j) 
      if(dd(l).lt.dd(nb)) nb=l
   24       continue
      n=ind(nb) 
      ind(nb)=ind(lb) 
      ind(lb)=n 
      call test(wt,n,xi,int,d,dd,sq,dt,mn-1,mn,nd)
      if(wt.lt.tol) go to 29
c                       avoid pivot for last x
      if(lb.lt.kx) go to 25 
      rs=xi(ky,ky)-xi(n,ky)*xi(n,ky)/xi(n,n)
      ci=ci+co(n) 
      go to 26
   25       call pivot(n,mn-1,lb+1,kx,
     1         xi,int,ind,ci,co(n),dt,1,nd,nn,ky,lb+1,kx) 
      rs=xi(ky,ky)
   26       if(wt.lt.tal) go to 28
c suspending rs check, as this apparently d
c     if(rs.ge.tul) go to 27
c     ne=6
c      return
c     27          mt=mn+1
      mt=mn+1
      call qstore(rs,ci,cl,rm,mb,df,mn,sig,pen,ib,ns)
   28    continue 
      ci=ci-co(n) 
   29    if(mt.le.mn) sel=tal 
      if(ipp(mn).le.kx) go to 32
      if(wt.lt.tal.or.ls.ge.kx-2) go to 30
      call pivot(n,0,kz,kx,xi,int,ind,ci,co(n),dt,1,nd,nn,ky,kz,0)
c                        begin branch and bound 
      mn=mh 
      call landb(xi,xn,d,dd,ci,cn,co,cl,rm,ind,ns,nd,df,
     1         ib,mb,ls,mn,ky,sig,pen,zip,kx,kz,int,iv) 
      di=real((rs-xi(ky,ky))/sq(ky))
      dx=amax1(dx,abs(di))
      call copy(xi,xn,int,ind,mn-1,ls,mn+kz-ls,nd,sc,ci,cn,0) 
c                               backtrack 
   30    if(mn.eq.1) go to 33 
      mn=mn-1 
      n=int(mn) 
      call pivot(n,mn-1,ipp(mn),
     1         kx,xi,int,ind,ci,-co(n),dt,0,nd,nn,ky,ipp(mn),kx)
      if(mn.le.nf.or.mn.gt.mt) go to 31 
      if(wt.lt.tol.or.mn.lt.mh) go to 32
      call trans(bnd,mr,rs,df,mn,sig,pen,ib,mb) 
      if(bnd.lt.rm(mr)) go to 32
   31    go to 30 
   32 go to 23
c                                output
   33 icur=1  
      do 42 lb=1,ns 
      if(cl(lb).eq.zip) go to 42
c                             decode labels 
      sbs(icur)=cl(lb)
      m=0 
      do 34 l=1,nv
      if(cl(lb).lt.co(l)) go to 34
      m=m+1 
      ind(m)=l
      cl(lb)=cl(lb)-co(l) 
   34    continue 
      ind(m+1)=ky 
c    
c      if(mod(lb-1,mb).eq.0.and.ib.eq.1) 
c     1  call intpr("n variables = ",-1,m,1)
c	call dblepr("regid = ",-1,sbs(icur),1)
c	call intpr("subset = ",-1,ind,m)
      call crit(rs,rm(lb),sq(ky),sig,df,pen,it,ib,m,rsret)
      rss(icur)=rsret
      icur=icur+1
      call copy(xn,rr,ind,ind,m+1,0,m+1,nd,sc,cn,zip,1) 
c                           invert submatrix
      do 37 l=1,m 
      n=ind(l)
      d(n)=xn(n,n)
      call pivot(n,l-1,l+1,m, 
     1         xn,ind,ind,cn,zip,dt,0,nd,nn,ky,l+1,m+1) 
   37    continue 
c                         regression statistics 
      rt=xn(ky,ky)/(df-float(m))
      do 39 l=1,m 
      n=ind(l)
      dd(l)=sc(ky)*xn(n,ky)/(d(n)*sc(n))
      beta=dd(l)
      tel=-d(n)*d(n)/sq(n)/xn(n,n)
      f=-xn(n,ky)*xn(n,ky)/(xn(n,n)*rt) 
   39    continue 
      di=real((rs-xn(ky,ky))/sq(ky))
      dx=amax1(dx,abs(di))
      if(it.ne.1) go to 42
c                               intercept 
      dint=xm(ky) 
      do 40 l=1,m 
      n=ind(l)
      dint=dint-dd(l)*xm(n) 
   40       continue
      yint=dint 
   42 continue
      di=real((sq(ky)-xi(ky,ky))/sq(ky))
      dx=amax1(dx,abs(di))
c  multreg next two lines 
c      if(dx.gt.tol)
c     1 call dblepr('Largest discrepancy observed for R**2 is',-1,dx,1)
      lb=1
      return
      end 
c                           leaps and bounds
      subroutine landb(xi,xn,d,yi,ci,cn,co,cl,rm,ind,ns,nd,df,ib, 
     1   mb,iz,mn,ky,sig,pen,zip,kx,kz,ni,iv) 
      double precision xi,xn,d,yi,cn,ci,co,cl,rm,zip,rs,sig,df,bnd,pen 
      dimension xn(nd,nd),xi(nd,nd),d(nd),yi(nd),co(nd),cl(ns),rm(ns),
     1   ni(nd),ind(nd) 
      ni(iz)=kx 
      ip=iz 
      yi(iz)=xi(ky,ky)
      if(iv.eq.0) go to 11
c                  initial p+1 variable regressions 
      do 10 l=iz,kx 
      n=ind(l)
      rs=xn(ky,ky)-xn(n,ky)*xn(n,ky)/xn(n,n)
      call qstore(rs,cn+co(n),cl,rm,mb,df,mn,sig,pen,ib,ns)
   10    continue 
c                            first fathoming
   11 ls=ip+1+iv
      do 12 lb=ls,kx
      call trans(bnd,mr,yi(ip),df,mn+kx-lb+iv,sig,pen,ib,mb)
      if(bnd.lt.rm(mr)) go to 13
   12    continue 
      ip=ip-1 
      go to 27
   13    ld=kx+ls-lb
      if(ip.gt.iz) go to 15 
c                           finish inversion
      do 14 l=iz,kx 
      n=ind(l)
      xi(n,n)=-d(n) 
      call pivot(n,0,iz,l-1,
     1            xi,ind,ind,ci,zip,d,3,nd,nn,ky,iz,0)
   14       continue
      lo=kx 
      go to 18
c                          extend prior pivots
   15    if(iv.eq.1.or.ld.le.ni(ip-1)) go to 17 
      il=ip-2 
      do 16 l=iz,il 
      if(ld.le.ni(l+1)) go to 16
      m=ind(l)
      if(xn(m,m).lt.zip) call pivot(m,0,l+1,ld,xn,
     1            ind,ind,cn,zip,d,3,nd,ni(l+1),ky,ni(l+1)+1,0) 
      if(xn(m,m).gt.zip) call pivot(m,0,l+1,ld,xi,
     1            ind,ind,ci,zip,d,3,nd,ni(l+1),ky,ni(l+1)+1,0) 
   16       continue
   17    lo=ld+iv*(kx-ld) 
      call pivot(n,0,ip,lo, 
     1      xi,ind,ind,ci,-co(n),d,0,nd,ni(ip),ky,ip,0) 
c                           bound regressions 
   18    do 21 lb=ip,lo 
      n=ind(lb) 
      rs=xi(ky,ky)-xi(n,ky)*xi(n,ky)/xi(n,n)
      if(ld.eq.kx) call qstore 
     1         (rs,ci-co(n),cl,rm,mb,df,kx-ip+mn-1,sig,pen,ib,ns) 
c                          re-order variables 
      m=lb
   19       if(m.eq.ip.or.rs.le.yi(m)) go to 20 
      yi(m+1)=yi(m) 
      ind(m)=ind(m-1) 
      m=m-1 
      go to 19
   20       ind(m)=n
      yi(m+1)=rs
   21    continue 
      if(ls.ge.kx) go to 27 
      ld=min0(ld,kx-1)
c                           second fathoming
      do 22 lb=ls,ld
      lc=ld+ls-lb 
      call trans(bnd,mr,yi(lc+1),df,mn+ld-lb+iv,sig,pen,ib,mb)
      if(bnd.lt.rm(mr)) go to 23
   22       continue
      go to 27
c                 fixed and p+1 variable regressions
   23       do 26 lb=ls,lc
      ip=lb-iv
      n=ind(ip-1) 
      call pivot(n,0,ip,lc+iv*(kz-lc)-1,
     1            xn,ind,ind,cn,co(n),d,0,nd,ni(ip),ky,ip,0)
      if(iv.eq.1) go to 24
      call qstore(xn(ky,ky),cn,cl,rm,mb,df,mn,sig,pen,ib,ns) 
      mn=mn+1 
      go to  26 
   24          mn=mn+1
      do 25 l=ip,kx 
      m=ind(l)
      rs=xn(ky,ky)-xn(m,ky)*xn(m,ky)/xn(m,m)
      call qstore(rs,cn+co(m),cl,rm,mb,df,mn,sig,pen,ib,ns)
   25          continue 
   26       continue
c                               backtrack 
   27    if(ip.le.iz) return
      n=ind(ip-1) 
      if(xn(n,n).lt.zip) go to 28 
      call pivot(n,0,ip,ni(ip),xi,
     1         ind,ind,ci,co(n),d,0,nd,nn,ky,ip,0)
      ip=ip-1 
      go to 27
   28    mn=mn-1
      call pivot(n,0,ip,ni(ip), 
     1      xn,ind,ind,cn,-co(n),d,0,nd,nn,ky,ip,0) 
      go to 11
      end 
c                           check tolerances
      subroutine test(wt,n,xi,int,d,dd,sq,dt,ml,mn,nd)

      double precision wt,xi,d,dd,sq,dt
      dimension xi(nd,nd),int(nd),d(nd),dd(nd),sq(nd),dt(nd)
      d(n)=xi(n,n)
      wt=d(n)/sq(n) 
      dd(mn)=wt*d(n)
      int(mn)=n 
      if(ml.eq.0) go to 11
      do 10 l=1,ml
      j=int(l)
c*****the next two lines of code are modified for cdc (?)
      dt(l)=xi(j,j) 
      if(xi(n,n).ne.0.)dt(l)=dt(l)-xi(n,j)*xi(n,j)/xi(n,n)
      if(-dd(l)/dt(l).lt.wt) wt=-dd(l)/dt(l)
   10    continue 
   11 return
      end 
c                             partial sweep 
      subroutine pivot(n,ml,ls,lx,xn,int,ind,cn,co,dt,iq,nd,nn,ky,ll,ly)
      double precision xn,b,dt,cn,co 
      dimension xn(nd,nd),dt(nd),ind(nd),int(nd)
      nn=lx 
      if(iq.eq.3) go to 10
      cn=cn+co
      xn(n,n)=-xn(n,n)
      xn(ky,ky)=xn(ky,ky)+xn(n,ky)*xn(n,ky)/xn(n,n) 
   10 if(ml.eq.0) go to 13
c                            gauss-jordan 
      do 12 l=1,ml
      j=int(l)
      b=xn(n,j)/xn(n,n) 
      if(iq.eq.0) xn(j,j)=xn(j,j)+b*xn(n,j) 
      if(iq.eq.1) xn(j,j)=dt(l) 
      do 11 m=ll,ly 
      k=ind(m)
      xn(j,k)=xn(j,k)+b*xn(n,k) 
      xn(k,j)=xn(j,k) 
   11       continue
   12    continue 
   13 if(ll.gt.lx) go to 16 
c                                gauss
      do 15 l=ll,lx 
      j=ind(l)
      b=xn(n,j)/xn(n,n) 
      do 14 m=ls,l
      k=ind(m)
      xn(j,k)=xn(j,k)+b*xn(n,k) 
      xn(k,j)=xn(j,k) 
   14       continue
      xn(j,ky)=xn(j,ky)+b*xn(n,ky)
      xn(ky,j)=xn(j,ky) 
   15    continue 
   16 return
      end 
c                   labels and saves best regressions 
      subroutine qstore(rs,cn,cl,rm,mb,df,mn,sig,pen,ib,ns) 
      double precision rs,cn,cl,rm,df,sig,pen,tr
      dimension cl(ns),rm(ns) 
      call trans(tr,mr,rs,df,mn,sig,pen,ib,mb)
      if(rm(mr).le.tr) return
      l=0
      do 10 m=1,mb
      l=mr-mb+m 
      if(cl(l).eq.cn) return
   10 continue
      if(mb.eq.1) go to 12
      do 11 m=2,mb
      if(tr.ge.rm(l-1)) go to 12
      rm(l)=rm(l-1) 
      cl(l)=cl(l-1) 
      l=l-1 
   11    continue 
   12 rm(l)=tr
      cl(l)=cn
      return
      end 
c                              matrix copy
      subroutine copy(xn,xi,int,ind,ml,mt,mp,nd,sc,cn,ci,iq)
      double precision xn,xi,sc,cn,ci
      dimension int(nd),ind(nd),xn(nd,nd),xi(nd,nd),sc(nd)
      do 12 l=1,mp
      if(l.le.ml) go to 10
      m=mt+l-ml-1 
      int(l)=ind(m) 
   10    k=int(l) 
      do 11 m=1,l 
      j=int(m)
      if(iq.eq.0) xn(j,k)=xi(j,k) 
      if(iq.eq.1) xn(j,k)=xi(j,k)/(sc(j)*sc(k)) 
      xn(k,j)=xn(j,k) 
   11    continue 
   12 continue
      cn=ci 
      return
      end 
c             transform of criterion-smaller must be better 
      subroutine trans(tr,mr,rs,df,mn,sig,pen,ib,mb)
      double precision rs,df,sig,pen,tr
      mr=mb 
c     go to (10,11,12),ib
      select case (ib)
      case (1)
         tr=rs 
         mr=mn*mb
         return
      case (2)
         tr=rs/(df-float(mn))
         return
      case (3)
         tr=rs+pen*float(mn)*sig 
         return
      end select
      end 
c                   criterion and r**2 from transform 
      subroutine crit(rs,tr,ss,sig,df,pen,it,ib,mn,rsrt) 
      double precision rs,tr,ss,sig,df,pen,cr,rsrt
c     go to (11,14,16),ib
      cr = 0.0
      select case (ib)
      case (1)
         rs=tr 
         cr=1.0-rs/ss
c         rsrt=cr
c	call dblepr('R^2 = ',-1,cr,1)
c         return
      case (2)
         cr=1.0-df*tr/ss 
         rs=tr*(df-float(mn))
c         go to 12
      case (3)
         cr=tr/sig-df+float(min0(it,1))*(pen-1.0)
         rs=tr-pen*float(mn)*sig 
c         go to 12
      end select
      rsrt=cr
      end 
