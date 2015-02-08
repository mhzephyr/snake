        subroutine descin
c**************************************************************
c                                                             *
c                          d e s c i n                        *
c                                                             *
c**************************************************************
c
c                read the description file
c                fill lhe commons /relfra/ , /fidef/ and /alter/
c
        include 'snake.inc'
c                next line : begening of the /plot/ package
        logical lbox,luser,bare
        character unit*3,name*2,rep*1
        common /cplot/unit(15),name(15)
        common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @              ,nplvold,h1(2,maxq),v1(2,maxq),nq
     @              ,ip,ih,iv,kp,bplot(3,maxplv)
     @              ,plobox(5,4,2,maxre),nface(maxre),lbox
     @              ,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @    ,iplpi(maxplv),iplpf(maxplv)
     @    ,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @		     ,luser,nclic,lclic(maxre),clight(4,maxre),
     @    epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3),nrs(maxre)
c                previous line : end of the /plot/ package
c                next line : begening of the /relfra/ package
        character epsw*4,rname*8,rnamer*8,title*80,cerfiln*20
        common /crelfr/epsw(maxre),rname(maxre),rnamer(maxre)
     @                 ,title(maxre),cerfiln(maxepr,maxre)
        common /relfra/ ver(maxve,12),xyz0(maxre,3)
     @                 ,frbox(2,maxre,3)
     @                 ,zang(maxre),xang(maxre),yang(maxre)
     @                 ,indve,indre,indep,ilive,adata(maxdat,maxre)
     @                 ,indrer(maxre),tept(maxepr,maxre)
     @                 ,yepmin(maxre),yepmax(maxre),kept(maxepr,maxre)
     @                 ,yepstp(maxre),yept(maxepr,maxre),nep(maxre)
     @                 ,colli(maxepr,maxre,6)
c                previous line : end of the /relfra/ package
        character*80 line*80,ch*1
c                next line : begening of the /fidef/ package
      logical cyl,lfil
        character*80 fifi
c
c  add in for adding two map files -JJL
c        
        character*80 file(2)
        real factor(2)
c  end add in  -JJL        
        common /fidef/method(maxre),itype(maxre),indfi(maxre)
     @               ,fact(maxre),cyl(maxre),lfil
        common /cfidef/ fifi(maxre)
c                previous line : end of the /fidef/ package
        common/alter/jflag,prec,hstmax,hstmin,sinc
c  lsynchro : compute average energy losses by synchrotron radiation?
c  lspin : compute spin precession?
      logical lsynchro,lspin
      common/logic/lsynchro,lspin
c  x_mag_mon: particle magnetic momentum in unit of ehbar/2m
c  e_over_m: particle cyclotron pulsation in rad s-1 T-1
c  twice_spin: 1. for spin 1/2 particle, 2. for spin 1 particle...
      common/particle/z_particle,xm
      common/spin/x_mag_mon,e_over_m,twice_spin,xmref
      dimension rel(3)
c
c
c                         min,max and increment of the step change
c		and mass of the particle (in GeV) for average synchrotron
c		radiation losses:
c
        read(2,*)hstmin,hstmax,sinc
      read(2,*)lspin,lsynchro
      read(2,*)z_particle,xm
      read(2,*)x_mag_mon,e_over_m,twice_spin,xmref
c
        indre=0
      iepz=0
c
c**************************new region:*************************
c
c       define the magnetic field in this region
c
 200     continue
        read(2,'(a)',end=300)line
        if(line(1:4).eq.'end*')goto300
        indre=indre+1
        if(indre.gt.maxre)stop'snake:max.number of regions surpassed'
        title(indre)=line
        write(6,*)title(indre)
c                name of this new region :
        read(2,'(a)')rname(indre)
c                name of a previous region ( or 'absolute') for incremental
c                         frame definition :
        read(2,'(a)')rnamer(indre)
c                absol. coord. of the relat. frame origine point:
        read(2,*)(xyz0(indre,i),i=1,3)
c                rotation angle/z axis (posit.def. when one goes from abs...
c                         ...to relat frame)
        read(2,*)zangd,xangd,yangd
        zang(indre)=zangd*dtor
c                rotation angle/x axis (posit.def. when one goes from abs...
c                         ...to relat frame)
        xang(indre)=xangd*dtor
c                rotation angle/y axis (posit.def. when one goes from abs...
c                         ...to relat frame)
        yang(indre)=yangd*dtor
c                perform , if necessary , the incrementation on
c                 xyz0,zang,yang and xang  of the new region :
        write(6,*)'    ',rname(indre),' refer.= ',rnamer(indre)
        if(rnamer(indre).ne.'absolute')then
                 call itrans
        endif
c                relat. coord. of the free-box
        read(2,*)(frbox(1,indre,i),i=1,3)
        read(2,*)(frbox(2,indre,i),i=1,3)
c		free-box : test max > min ?
      do 1 i=1,3
         if(frbox(1,indre,i).ge.frbox(2,indre,i))then
	        	write(6,*)' descin:zero or negative box volume'
			stop
         endif
 1    continue
c                method of transport:1:no field,2:matrix,3:rung-khuta
        read(2,*)method(indre)
c              type of matrix/field:1:--------,1st order,constant
c                                     2:--------,2nd order,analytic
c                                      3:--------,---------,read on a file
        read(2,*)itype(indre)
c    index of the matrix/field in this type(over read if field on a file):
        read(2,*)indfi(indre)
c Adding two 3-D maps -JJL
      if((itype(indre).eq.3).and.(indfi(indre).eq.7))then
         read(2,*)file(1),factor(1)
c         write(6,*)file(1),factor(1)
         read(2,*)file(2),factor(2)
         read(2,*)fifi(indre)
         call mapadd(file,factor,fifi(indre))
         go to 747
       endif 
c end of adding maps -JJL                
c                field file name:
        read(2,'(a)')line
      if(line(1:1).eq.'#')then
		fifi(indre)=line
      else
		i1=1
		i2=i1
        call nextcoma(i2,line)
        read(line(i1:i2-1),'(a)')fifi(indre)
      endif
c		cartesian or cylindrical coors.syst for box and map?
c				(but not for end-plans !)
 747	read(2,'(a)')rep
      cyl(indre)=rep.eq.'c'
c		test radius > 0. and convert degrees into radians :
      if(cyl(indre))then
        if(frbox(1,indre,1).le.0.)then
         write(6,*)' descin: zero or negative minimum',
     @				' radius for a cylindrical region'
         stop
      endif
		frbox(1,indre,2)=frbox(1,indre,2)*dtor
		frbox(2,indre,2)=frbox(2,indre,2)*dtor
      endif
c   multiplicative factor of the field(result in tesla) if method=3...
c   ....or reference momentum if method=2:
        read(2,*)fact(indre)
c     number of additional data to describe the region:
        read(2,*)nad
        if(nad.gt.maxdat)stop ' descin: nb.of additional data >maxdat'
c     additional data
        do  8 i=1,maxdat
                 adata(i,indre)=0.
                 if(i.le.nad)read(2,*)adata(i,indre)
  8   continue
c        define the end-planes of this mew region(yep refers to relat.coor.)
      do 421 iep=1,maxepr
 421		cerfiln(iep,indre)='none'
        read(2,'(a)')epsw(indre)
        if( epsw(indre).eq.'loop')then
             read(2,*)yepmin(indre),yepmax(indre),yepstp(indre)
             nep(indre)=((yepmax(indre)-yepmin(indre))/yepstp(indre))+1
	     if(nep(indre).gt.maxepr)stop
     s			' descin:too many end-planes in this region !!!'
        else
             nep(indre)=0
             do201 i=1,maxepr+1
                          read(2,'(a)')line
                          if(line(1:3).eq.'eol')goto202
                          nep(indre)=i
			  i1=1
			  i2=i1
			  call nextcoma(i2,line)
                          read(line(i1:i2-1),'(a)')ch
			  i1=i2+1
			  i2=i1
			  call nextcoma(i2,line)
                          read(line(i1:i2-1),*)yept(i,indre)
			  i1=i2+1
			  i2=i1
			  call nextcoma(i2,line)
                          read(line(i1:i2-1),*)tept(i,indre)
			  i1=i2+1
			  i2=i1
			  call nextcoma(i2,line)
                          read(line(i1:i2-1),'(a)')cerfiln(i,indre)
             do icol=1,6
			    i1=i2+1
			    i2=i1
			    call nextcoma(i2,line)
                            read(line(i1:i2-1),*)colli(i,indre,icol)
             enddo
              if(ch.eq.'x')then
				kept(i,indre)=1
               elseif(ch.eq.'y')then
				kept(i,indre)=2
               else
				kept(i,indre)=3
              endif
 201      tept(i,indre)=tept(i,indre)*dtor
             write(6,*)' descin:too many end-planes in this region !!!'
             stop
 202          if(nep(indre).le.0)then
                     write(6,*)' descin:no end plane in this region !!!'
                     stop
             endif
             yepmin(indre)=yept(1,indre)
             yepmax(indre)=yept(nep(indre),indre)
        endif
c        if(yepmin(indre).gt.yepmax(indre))
c     @  write(6,*)' warning : y end-plane max < min'
c        if(yepmin(indre).lt.frbox(1,indre,2)
c     @  .or.yepmax(indre).gt.frbox(2,indre,2))
c     @       write(6,*)' warning : free-box < end-plane'
c         store the abs. coor. of the submits of the free-box:
        ir=indre
        it=0
      do 5 ix=1,2
          do 5 iy=1,2
            do 5 iz=1,2
           if(cyl(ir))then
			frx=frbox(ix,ir,1)*cos(frbox(iy,ir,2))
			fry=frbox(ix,ir,1)*sin(frbox(iy,ir,2))
			frz=frbox(iz,ir,3)
           else
			frx=frbox(ix,ir,1)
			fry=frbox(iy,ir,2)
			frz=frbox(iz,ir,3)
           endif
              it=it+1
 5             call plorot(frx,fry,frz
     s        ,boxp(it,ir,1),boxp(it,ir,2),boxp(it,ir,3),xb,yb,zb)
c region frames:
	
		nrs(ir)=0
c plot of region frames not yet implemented in cylindrical regions:
      if(cyl(ir))goto101
      do iq=1,3
			iqp=mod(iq,3)+1
			iqm=mod(iq+1,3)+1
            if(frbox(1,ir,iqp)*frbox(2,ir,iqp).lt.0.
     s			.and.frbox(1,ir,iqm)*frbox(2,ir,iqm).lt.0.
     s			.and.frbox(2,ir,iq).gt.0.)then
				nrs(ir)=nrs(ir)+1
c   starting point:
				rel(iq)=max(0.,frbox(1,ir,iq))
				rel(iqp)=0.
				rel(iqm)=0.
               call plorot(rel(1),rel(2),rel(3),
     s				regp(1,nrs(ir),ir,1),
     s				regp(1,nrs(ir),ir,2),
     s				regp(1,nrs(ir),ir,3),
     s                           xxb,yyb,zzb)
c   ending point:
				rel(iq)=frbox(2,ir,iq)
				rel(iqp)=0.
				rel(iqm)=0.
              call plorot(rel(1),rel(2),rel(3),
     s				regp(2,nrs(ir),ir,1),
     s				regp(2,nrs(ir),ir,2),
     s				regp(2,nrs(ir),ir,3),
     s                           xxb,yyb,zzb)
              endif
         enddo

c end planes:
      do 100 iepr=1,nep(ir)
		iepz=iepz+1
        if(iepz.gt.maxep)stop'descin:too many end planes!'
		neps(iepz)=0
c plot of end planes not yet implemented in cylindrical regions:
        if(cyl(ir))goto 100
        if(epsw(ir).eq.'loop')then
c
c	case 'loop':
c
			  yep=yepmin(ir)+yepstp(ir)*(iepr-1)
			  do iq=1,3,2
			  iqp=mod(iq,3)+1
			  iqm=mod(iq+1,3)+1
c   starting point:
			    rel(1)=0.
			    rel(3)=0.
			    rel(iq)=frbox(1,ir,iq)
			    rel(2)=yep
			    if(rel(iqp).gt.frbox(1,ir,iqp)
     s			    .and.rel(iqp).lt.frbox(2,ir,iqp))then
			        neps(iepz)=neps(iepz)+1
			        call plorot(rel(1),rel(2),rel(3),
     s				epp(1,neps(iepz),iepz,1),
     s				epp(1,neps(iepz),iepz,2),
     s				epp(1,neps(iepz),iepz,3),
     s                           xxb,yyb,zzb)
c   ending point:
			        rel(iq)=frbox(2,ir,iq)
			        call plorot(rel(1),rel(2),rel(3),
     s				epp(2,neps(iepz),iepz,1),
     s				epp(2,neps(iepz),iepz,2),
     s				epp(2,neps(iepz),iepz,3),
     s                           xxb,yyb,zzb)
			     endif
               enddo
      else
c
c	case 'list':
c  1/segment along "iqp":
			yep=yept(iepr,ir)
			iq=kept(iepr,ir)
			iqp=mod(iq,3)+1
			iqm=mod(iq+1,3)+1
			tep=tept(iepr,ir)
c   starting point:
			rel(iq)=yep
			rel(iqp)=frbox(1,ir,iqp)
			rel(iqm)=0.
            if(rel(iq).gt.frbox(1,ir,iq)
     s			.and.rel(iq).lt.frbox(2,ir,iq)
     s			.and.rel(iqm).gt.frbox(1,ir,iqm)
     s			.and.rel(iqm).lt.frbox(2,ir,iqm))then
			        neps(iepz)=neps(iepz)+1
              call plorot(rel(1),rel(2),rel(3),
     s				epp(1,neps(iepz),iepz,1),
     s				epp(1,neps(iepz),iepz,2),
     s				epp(1,neps(iepz),iepz,3),
     s                           xxb,yyb,zzb)
c   ending point:
			        rel(iq)=yep
			        rel(iqp)=frbox(2,ir,iqp)
			        rel(iqm)=0.
                  call plorot(rel(1),rel(2),rel(3),
     s				epp(2,neps(iepz),iepz,1),
     s				epp(2,neps(iepz),iepz,2),
     s				epp(2,neps(iepz),iepz,3),
     s                           xxb,yyb,zzb)
                endif
c  2/segment along "iqm/iq":
			sit=sin(tep)
           if(sit.ne.0.)then
				xl1=(frbox(1,ir,iq)-yep)/sit
				xl2=(frbox(2,ir,iq)-yep)/sit
            else
				xl1=-1.e30
				xl2=1.e30
           endif
			cot=cos(tep)
            if(cot.ne.0.)then
				xl3=frbox(1,ir,iqm)/cot
				xl4=frbox(2,ir,iqm)/cot
            else
				xl1=1.e30
				xl2=1.e30
            endif
			neps(iepz)=neps(iepz)+1
c   starting point:
			xl=max(xl1,xl2)
			xlm=max(xl3,xl4)
			xlmax=min(xl,xlm)
			rel(iq)=yep+xlmax*sit
			rel(iqp)=0.
			rel(iqm)=xlmax*cot
            call plorot(rel(1),rel(2),rel(3),
     s				epp(1,neps(iepz),iepz,1),
     s				epp(1,neps(iepz),iepz,2),
     s				epp(1,neps(iepz),iepz,3),
     s                           xxb,yyb,zzb)
c   ending point:
			xl=min(xl1,xl2)
			xlm=min(xl3,xl4)
			xlmin=max(xl,xlm)
			rel(iq)=yep+xlmin*sit
			rel(iqp)=0.
			rel(iqm)=xlmin*cot
            call plorot(rel(1),rel(2),rel(3),
     s				epp(2,neps(iepz),iepz,1),
     s				epp(2,neps(iepz),iepz,2),
     s				epp(2,neps(iepz),iepz,3),
     s                           xxb,yyb,zzb)
        endif
 100	continue
 101	continue
        goto200
c
 300     continue
        write(6,*)
        write(6,*)
c      save /cplot/,/plot/,/crelfr/,/relfra/,/fidef/,/cfidef/,/alter/,
c     &     /logic/,/particle/,/spin/	  
        return
        end
      subroutine nextcoma(i,line)
c*************************************************************************
c                                                                        *
c                       n e x t c o m a                                  *
c                                                                        *
c*************************************************************************
c
c	Updates i to the location of the first occurence in 'line' 
c		of a coma at right of the ith. character 
	character line * (*)
	j=i
      do 1 i=j+1,80
		if(line(i:i).eq.',')return
 1      continue
      i=i+1
      return
      end
      subroutine plotin(i,x,y,z,bx,by,bz,px,py,pz,flag)
c*************************************************************************
c                                                                        *
c                       p l o t i n                                      *
c                                                                        *
c*************************************************************************
c
c		add a plot point to the ith trajectory wile keeping
c			the old trajectories
c		if(flag) : create a new plot trajectory and fill
c			its first point
c
          include 'snake.inc'
c                next line : begening of the /plot/ package
        logical lbox,luser,bare
        character unit*3,name*2
        common /cplot/unit(15),name(15)
        common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @              ,nplvold,h1(2,maxq),v1(2,maxq),nq
     @              ,ip,ih,iv,kp,bplot(3,maxplv)
     @              ,plobox(5,4,2,maxre),nface(maxre),lbox
     @              ,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @    ,iplpi(maxplv),iplpf(maxplv)
     @    ,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @		     ,luser,nclic,lclic(maxre),clight(4,maxre),
     @    epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3),nrs(maxre)
c                previous line : end of the /plot/ package
c
      logical flag
        if(i+nplvold.gt.maxplv)return
      iplv=i+nplvold
      if(flag)then
		nplv=iplv
		nplp(iplv)=0
      endif
        if(nplp(iplv).lt.maxplp)then
                nplp(iplv)=nplp(iplv)+1
		iplp=nplp(iplv)
                plott(iplp,iplv,1)=x
                plott(iplp,iplv,2)=y
                plott(iplp,iplv,3)=z
                plott(iplp,iplv,4)=bx
                plott(iplp,iplv,5)=by
                plott(iplp,iplv,6)=bz
                plott(iplp,iplv,7)=px
                plott(iplp,iplv,8)=py
                plott(iplp,iplv,9)=pz
        endif
c      save /cplot/,/plot/
      return
      end
        subroutine box
c*************************************************************************
c                                                                        *
c                       b o x                                            *
c                                                                        *
c*************************************************************************
c
c                pre-plot the free boxes
c
c       input:
c             frbox(maxi or mini,region index,x y or z):relat.coord.
c                of the free-box
c             ip:index of the absol. axis(1:x,2:y,3:z)which is orthog./screen
c             ih:       "       "           "           "      horizontal "
c             iv:       "       "           "           "      vertical   "
c             kp= 1:observer is at ip>0.
c                -1:   "     "  "  ip<0.
c       output:
c             plobox(submit index,horiz. or vertic.,region index):abs.coord
c                of submit of free-box(special submit list for 3d projection)
          include 'snake.inc'
c                next line : begening of the /plot/ package
        logical lbox,luser,bare
        character unit*3,name*2
        common /cplot/unit(15),name(15)
        common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @              ,nplvold,h1(2,maxq),v1(2,maxq),nq
     @              ,ip,ih,iv,kp,bplot(3,maxplv)
     @              ,plobox(5,4,2,maxre),nface(maxre),lbox
     @              ,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @    ,iplpi(maxplv),iplpf(maxplv)
     @    ,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @		     ,luser,nclic,lclic(maxre),clight(4,maxre),
     @    epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3),nrs(maxre)
c                previous line : end of the /plot/ package
        common /absfra/ vea(maxve,0:maxep,21),nuve,nuep,nure
c                next line : begening of the /relfra/ package
        character epsw*4,rname*8,rnamer*8,title*80,cerfiln*20
        common /crelfr/epsw(maxre),rname(maxre),rnamer(maxre)
     @                 ,title(maxre),cerfiln(maxepr,maxre)
        common /relfra/ ver(maxve,12),xyz0(maxre,3)
     @                 ,frbox(2,maxre,3)
     @                 ,zang(maxre),xang(maxre),yang(maxre)
     @                 ,indve,indre,indep,ilive,adata(maxdat,maxre)
     @                 ,indrer(maxre),tept(maxepr,maxre)
     @                 ,yepmin(maxre),yepmax(maxre),kept(maxepr,maxre)
     @                 ,yepstp(maxre),yept(maxepr,maxre),nep(maxre)
     @                 ,colli(maxepr,maxre,6)
c                previous line : end of the /relfra/ package
      common /comlight/cl(3)
        dimension ls(5,6)
        dimension nh(3),nv(3),np(3)
        data ls/1,3,4,2,1,
     @          5,6,8,7,5,
     @          1,2,6,5,1,
     @		3,7,8,4,3,
     @		1,5,7,3,1,
     @		2,4,8,6,2/
c
c
      ul2=0.
         do 5 i=1,3
                 nh(i)=0.
                 nv(i)=0.
 5        ul2=ul2+cl(i)*cl(i)
      ul=sqrt(ul2)
        nh(ih)=1
        nv(iv)=1
        np(1)=nh(2)*nv(3)-nh(3)*nv(2)
        np(2)=nh(3)*nv(1)-nh(1)*nv(3)
        np(3)=nh(1)*nv(2)-nh(2)*nv(1)
        do 6 i=1,3
                 if(np(i).ne.0)then
                          ip=i
                          kp=np(i)
                 endif
 6       continue
       do 4 ir=1,nure
          indre=ir
	  nface(ir)=0 
	  do8if=1,6
c		  put this face in the screen (h*v*p) coord. syst.:
		hl1=boxp(ls(1,if),indre,ih)
		vl1=boxp(ls(1,if),indre,iv)
		pl1=boxp(ls(1,if),indre,ip)*float(kp)
		hl2=boxp(ls(2,if),indre,ih)
		vl2=boxp(ls(2,if),indre,iv)
		pl2=boxp(ls(2,if),indre,ip)*float(kp)
		hl3=boxp(ls(3,if),indre,ih)
		vl3=boxp(ls(3,if),indre,iv)
		pl3=boxp(ls(3,if),indre,ip)*float(kp)
c		  is this face seen direct or inverse ?
		dh2=hl2-hl1
		dv2=vl2-vl1
		dp2=pl2-pl1
		dh3=hl3-hl1
		dv3=vl3-vl1
		dp3=pl3-pl1
		sh=dv2*dp3-dp2*dv3
		sv=dp2*dh3-dh2*dp3
		sp=dh2*dv3-dv2*dh3
        if(sp.lt.0.)then
c			   this face has to be plotted :
			nface(ir)=nface(ir)+1
            do 7 is=1,5
			  plobox(is,nface(ir),1,indre)=
     @			    boxp(ls(is,if),indre,ih)
 7			  plobox(is,nface(ir),2,indre)=
     @			    boxp(ls(is,if),indre,iv)
			arg=sh*sh+sv*sv+sp*sp
           if(arg.ne.0.)then
				coef=1./(sqrt(arg)*ul)
           else
				coef=1.
           endif
			clight(nface(ir),indre)=coef*(sh*cl(ih)+
     s			sv*cl(iv)+sp*cl(ip)*float(kp))
           if(nface(ir).ge.4)goto 4
           endif
 8	  continue
 4       continue
c      save /cplot/,/plot/,/absfra/,/crelfr/,/relfra/,/comlight/
        return
        end
        subroutine catch(ah,av,nreg)
c*************************************************************************
c                                                                        *
c                       c a t c h                                        *
c                                                                        *
c*************************************************************************
c
c                catch the free boxes whose a plotted submit is the 
c		narest from cusor location
c
c       input:
c		ah,av(horizontal and vertical cursor lcation 
c			in physical gks coordinates)
c       output:
c		nreg(number of the region whose free box has been caugth)
          include 'snake.inc'
c                next line : begening of the /plot/ package
        logical lbox,luser,bare
        character unit*3,name*2
        common /cplot/unit(15),name(15)
        common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @              ,nplvold,h1(2,maxq),v1(2,maxq),nq
     @              ,ip,ih,iv,kp,bplot(3,maxplv)
     @              ,plobox(5,4,2,maxre),nface(maxre),lbox
     @              ,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @    ,iplpi(maxplv),iplpf(maxplv)
     @    ,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @		    ,luser,nclic,lclic(maxre),clight(4,maxre),
     @    epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3),nrs(maxre)
c                previous line : end of the /plot/ package
        common /absfra/ vea(maxve,0:maxep,21),nuve,nuep,nure
      logical flag
      d2min=1.e10
      do 1 i=1,nure
      do 2 j=1,5
        do 2 k=1,nface(i)
		d2=(ah-plobox(j,k,1,i))**2+(av-plobox(j,k,2,i))**2
        if(d2.lt.d2min)then
			d2min=d2
			nreg=i
         endif
 2	   continue
 1      continue
c	  save the number of clicked boxes and their list in chronlogical
c		order :
c  supress redundency:
      flag=.false.
      do 3 ic=1,nclic
        if(lclic(ic).eq.nreg)then
			nclic=nclic-1
			flag=.true.
        endif
        if(flag)then
			if(ic.le.nclic)lclic(ic)=lclic(ic+1)
        endif
 3      continue
      nclic=nclic+1
      lclic(nclic)=nreg
c      save /cplot/,/plot/,/absfra/ 
      return
      end
        subroutine unclip
c*************************************************************************
c                                                                        *
c                       u n c l i p                                      *
c                                                                        *
c*************************************************************************
c
c                unclip the plot vectors :
c
          include 'snake.inc'
        common /absfra/ vea(maxve,0:maxep,21),nuve,nuep,nure
c                next line : begening of the /plot/ package
        logical lbox,luser,bare
        character unit*3,name*2
        common /cplot/unit(15),name(15)
        common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @              ,nplvold,h1(2,maxq),v1(2,maxq),nq
     @              ,ip,ih,iv,kp,bplot(3,maxplv)
     @              ,plobox(5,4,2,maxre),nface(maxre),lbox
     @              ,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @    ,iplpi(maxplv),iplpf(maxplv)
     @    ,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @		     ,luser,nclic,lclic(maxre),clight(4,maxre),
     @    epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3),nrs(maxre)
c                previous line : end of the /plot/ package
c
      do 1 i=1,nplv
		iplpi(i)=1
 1		iplpf(i)=nplp(i)
c      save /cplot/,/plot/
      return
      end
        subroutine clip(sup,iqb,qb)
c*************************************************************************
c                                                                        *
c                           c l i p                                      *
c                                                                        *
c*************************************************************************
c
c                clip the plot vectors :
c
          include 'snake.inc'
        common /absfra/ vea(maxve,0:maxep,21),nuve,nuep,nure
c                next line : begening of the /plot/ package
        logical lbox,luser,bare
        character unit*3,name*2
        common /cplot/unit(15),name(15)
        common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @              ,nplvold,h1(2,maxq),v1(2,maxq),nq
     @              ,ip,ih,iv,kp,bplot(3,maxplv)
     @              ,plobox(5,4,2,maxre),nface(maxre),lbox
     @              ,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @    ,iplpi(maxplv),iplpf(maxplv)
     @    ,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @		     ,luser,nclic,lclic(maxre),clight(4,maxre),
     @    epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3),nrs(maxre)
c                previous line : end of the /plot/ package
c
      do 1 i=1,nplv
		js=1
        do 2 j=1,nplp(i)
          if((qb-plott(j,i,iqb))*sup.gt.0.)goto3
 2         js=j
 3      iplpi(i)=js
		js=nplp(i)
       do 4 j=nplp(i),iplpi(i),-1
            if((qb-plott(j,i,iqb))*sup.gt.0.)goto 5
 4			js=j
 5		iplpf(i)=js
 1      continue			
c      save /absfra/,/cplot/,/plot/
      return
      end
        subroutine pretra(*)
c*************************************************************************
c                                                                        *
c                       p r e t r a                                      *
c                                                                        *
c*************************************************************************
c
c                find the initial window
c       input:
c             nplv<=maxplv:number of trajectories to be plotted
c             nplp(traj.index)<=maxplp:number of submits in the plot-line
c             plott(submit index,traj.index,x y z bx by bz...
c                    ...px py pz):absol.coord.of submit,
c				  absol. field components and
c				  absol. spin components at this point
c             ip:index of the absol. axis(1:x,2:y,3:z)which is orthog./screen
c             ih:       "       "           "           "      horizontal "
c             iv:       "       "           "           "      vertical   "
c       output:
          include 'snake.inc'
        common /absfra/ vea(maxve,0:maxep,21),nuve,nuep,nure
c                next line : begening of the /plot/ package
        logical lbox,luser,bare
        character unit*3,name*2
        common /cplot/unit(15),name(15)
        common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @              ,nplvold,h1(2,maxq),v1(2,maxq),nq
     @              ,ip,ih,iv,kp,bplot(3,maxplv)
     @              ,plobox(5,4,2,maxre),nface(maxre),lbox
     @              ,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @    ,iplpi(maxplv),iplpf(maxplv)
     @    ,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @		     ,luser,nclic,lclic(maxre),clight(4,maxre),
     @    epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3),nrs(maxre)
c                previous line : end of the /plot/ package
c transformation world-->cm initiale
      common /initiale/ coxi,shxi,coyi,shyi
c transformation world-->cm courante
      common /courante/ cox,shx,coy,shy,zmarge
c transformation world<-->normalisee courante
      common/ntw/xmax1,xmin1,ymax1,ymin1,xnmin,xnmax,ynmin,ynmax
c menu:
      character*10 menu(6)
      data menu/'initial','zoom','unzoom','hard copy',
     s   'paint','exit'/
        logical hsvu,hsvz
        dimension work(maxplp)
        parameter (fracti=10.)
        character formh*15,formv*15
        common/ctext/formh,formv
      DIMENSION XZ(8),YZ(8)
c fleche diagonale en cm:
      DATA XZ/.1,0.,.5,0.,21.,20.5,21.,20.9/
      DATA YZ/.5,0.,.1,0.,29.7,29.6,29.7,29.2/
c
c                find the plot-window according to ih and iv:
c
c      write(6,*)'In pretra unit=',unit
c      write(10,*)'In pretra name=',name
c      write(10,*)(plott(5,5,l),l=1,9)
        hsvu=unit(iv).eq.unit(ih)
        xmin= 1.e10
        xmax=-1.e10
        ymin= 1.e10
        ymax=-1.e10
        if(lbox)then
c             then scale also on the free-boxes (11 submits * nure regions)
         do 4 ire=1,nure
         do 4 is=1,5
         do 4 if=1,nface(ire)
              if(plobox(is,if,1,ire).lt.xmin)xmin=plobox(is,if,1,ire)
              if(plobox(is,if,1,ire).gt.xmax)xmax=plobox(is,if,1,ire)
              if(plobox(is,if,2,ire).lt.ymin)ymin=plobox(is,if,2,ire)
4             if(plobox(is,if,2,ire).gt.ymax)ymax=plobox(is,if,2,ire)
        endif
c                then scale on the plot vectors:
        do 1 indve=1,nplv
            np=nplp(indve)
            if(ih.eq.10)then
                 call plotl(indve,work)
            elseif(ih.eq.11.or.ih.eq.12.or.ih.eq.15)then
                 call plotm(indve,ih,work)
            elseif(ih.eq.13)then
                 call plota(indve,work)
            elseif(ih.eq.14)then
                 call plotra(indve,work)
            else
                 call plotg(indve,ih,work)
            endif
            do 2 n=1,np
              if(work(n).lt.xmin)xmin=work(n)
 2             if(work(n).gt.xmax)xmax=work(n)
            if(iv.eq.10)then
                   call plotl(indve,work)
            elseif(iv.eq.11.or.iv.eq.12.or.iv.eq.15)then
                   call plotm(indve,iv,work)
            elseif(iv.eq.13)then
                 call plota(indve,work)
            elseif(iv.eq.14)then
                 call plotra(indve,work)
            else
                 call plotg(indve,iv,work)
            endif
            do 3 n=1,np
              if(work(n).lt.ymin)ymin=work(n)
 3             if(work(n).gt.ymax)ymax=work(n)
 1         continue
c add 1% margins:
      hm=1.e-2*(xmax-xmin)
      xmin=xmin-hm
      xmax=xmax+hm
      vm=1.e-2*(ymax-ymin)
      ymin=ymin-vm
      ymax=ymax+vm
c
            hsvz=hsvu
c                open gks and graphic screen
      xmini=xmin
      xmaxi=xmax
      ymini=ymin
      ymaxi=ymax
      iwr=1
      zmarge=2.
      call map1(xmin,ymin,xmax,ymax,10,21.,28.5,zmarge,hsvz,
     s		cox,shx,coy,shy)
      coxi=cox
      shxi=shx
      coyi=coy
      shyi=shy
 29      CALL IGSSE(6,1)
      call map1(xmin,ymin,xmax,ymax,10,21.,28.5,zmarge,hsvz,
     s		cox,shx,coy,shy)
c  higz: the size of the physical window on the wk is fixed by the ith. line
c of the file 'higz_windows.dat' (if it exists ) , where 'i' is the second
c argument of igsse. in that line data are: h(def=0),v(def=0),dh(def=600) and
c  dv(def=600). Unit is pixel. (h,v) is the position of the upper left corner,
c  (dh,dv) is the size of the window. for sun h<1052, v<913. Inside this 
c window, higz plots on the biggest viewport it can draw, keeping its dimensions
c  h' and v' in the same ratio as d2 and d3 .(Keeping also the aspect ratio)
        if(xmin.eq.xmax.or.ymin.eq.ymax)then
      call igend
c      		CALL IGSA (1)   :facultatif en multi-fenetrage
                 write(6,*)' pretra:can''t scale'
                 return1
        endif
 30      continue
        CALL ICLRWK(0,1)
      CALL IGRNG(21.,29.7)
c plot en coord. world:
      call iselnt(2)
      call isclip(1)
        call isln(3)
      if(.not.bare)then
c	write(6,*)'calling ipl'
            do 32 i=1,nq
c          write(10,*)h1(1,i),v1(1,i)
 32         call ipl(2,h1(1,i),v1(1,i))
      endif
        call isln(1)
      call iselnt(0)
      if(.not.bare)call grad
      call iselnt(2)
c gestion d'eventuels overflow (+-30 normalise):
		alp=(xmax1-xmin1)/(xnmax-xnmin)
		xp=xmin1+(30.-xnmin)*alp
		xm=xmin1-(30.+xnmin)*alp
		alp=(ymax1-ymin1)/(ynmax-ynmin)
		yp=ymin1+(30.-ynmin)*alp
		ym=ymin1-(30.+ynmin)*alp
c
c
      call traj(xp,xm,yp,ym)
c
c
c plot en cm:
      call iselnt(1)
c      CALL IPL(3,XZ,YZ)
c      CALL IPL(3,XZ(6),YZ(6))
      if(iwr.eq.2)then
		iwr=1
      call igend
      close(10)
      goto 29
      endif
c plot du menu en cm:
c      call iselnt(1)
      call ischh(1.)
c      call istxci(1)
      dx=20.9/6.
      dz=dx/10.
      y1=28.5+.1
      y2=29.7-dz
      x1=.1
      do i=1,6
        x2=x1+dx-dz
        call igpave(x1,x2,y1,y2,dz,0.,0.,'TR')
      write(6,*)'calling itx ', menu(i)
      write(6,*)'params=',x1+dz,y1+dz/2.
        call itx(x1+dz,y1+dz/2.,menu(i))
        x1=x1+dx
      enddo
 31      call igloc(1,nt,ibn,xndc,yndc,xwc,ywc)
 33      xw=xndc*29.7
      yw=yndc*29.7
      if(yw.ge.y1.and.yw.le.y2)then	
		x1=.1
      do i=1,6
			x2=x1+dx-dz
			if(xw.gt.x1.and.xw.le.x2)
     s				goto(101,102,103,104,105,106),i
		x1=x1+dx
      enddo
      write(6,*)'mal vise en x'
      goto 31
      endif
      write(6,*)'mal vise en y'
      goto 31
c initial:
 101      write(6,*)menu(1)
      hsvz=hsvu
      xmin=xmini
      xmax=xmaxi
      ymin=ymini
      ymax=ymaxi
      cox=coxi
      shx=shxi
      coy=coyi
      shy=shyi
      call map1(xmini,ymini,xmaxi,ymaxi,10,21.,28.5,
     s		zmarge,hsvz,coxi,shxi,coyi,shyi)
      goto 30
c zoom:
 102      write(6,*)menu(2)
      call igloc(1,nt,ibn,xndc1,yndc1,xwc,ywc)
		x1=xndc1*29.7
		y1=yndc1*29.7
  	        call cross(x1,y1)
      call igloc(1,nt,ibn,xndc2,yndc2,xwc,ywc)
		x2=xndc2*29.7
		y2=yndc2*29.7
      if(abs(x2-x1).lt..1.and.abs(y2-y1).lt..1)then
			hsvz=.false.
			call igloc(1,nt,ibn,xndc2,yndc2,xwc,ywc)
      else
			hsvz=hsvu
      endif
		xmin=(min(xndc1,xndc2)-zmarge/29.7)/cox+shx
		ymin=(min(yndc1,yndc2)-zmarge/29.7)/coy+shy
		xmax=(max(xndc1,xndc2)-zmarge/29.7)/cox+shx
		ymax=(max(yndc1,yndc2)-zmarge/29.7)/coy+shy
		call map1(xmin,ymin,xmax,ymax,10,21.,28.5,zmarge,hsvz,
     s		cox,shx,coy,shy)
      goto 30
c unzoom:
 103      write(6,*)menu(3)
      call igloc(1,nt,ibn,xndc,yndc,xwc,ywc)
		x1=(xndc-zmarge/29.7)/cox+shx
		y1=(yndc-zmarge/29.7)/coy+shy
		call igloc(1,nt,ibn,xndc,yndc,xwc,ywc)
		x2=(xndc-zmarge/29.7)/cox+shx
		y2=(yndc-zmarge/29.7)/coy+shy
		alp=(xmax-xmin)/(max(x1,x2)-min(x1,x2))
		bet=xmin-alp*min(x1,x2)
		xmin=alp*xmin+bet
		xmax=alp*xmax+bet
		alp=(ymax-ymin)/(max(y1,y2)-min(y1,y2))
		bet=ymin-alp*min(y1,y2)
		ymin=alp*ymin+bet
		ymax=alp*ymax+bet
      call map1(xmin,ymin,xmax,ymax,10,21.,28.5,zmarge,hsvz,
     s		cox,shx,coy,shy)
      goto 30
c hard copy:
 104      write(6,*)menu(4)
      iwr=2
      call igend
      open(10,file='higz.ps',status='unknown')
      call iopks(6)
c  		 PostScript metafile  A4 Portrait (-111):
      call iopwk(1,10,-111)
      call iacwk(1)
      call map1(xmin,ymin,xmax,ymax,10,21.,28.5,zmarge,hsvz,
     s		cox,shx,coy,shy)
      go to 30
c paint:
 105      write(6,*)menu(5)
      call igloc(1,nt,ibn,xndc,yndc,xwc,ywc)
      if(yndc*29.7.ge.y1)goto 33
		ah=(xndc-zmarge/29.7)/cox+shx
		av=(yndc-zmarge/29.7)/coy+shy
                call catch(ah,av,nreg)
      call iselnt(2)
      call paint(nreg)
      goto 105
c exit:
 106      write(6,*)menu(6)
      call igend
c      		CALL IGSA (1)   :facultatif en multi-fenetrage
c      save /absfra/,/cplot/,/plot/,/initiale/,/courante/,/ntw/,/ctext/
      return
      END
      subroutine map1(xmin,ymin,xmax,ymax,n,xsize1,ysize1,zmarge,hsvz,
     s		cox,shx,coy,shy)
c*************************************************************************
c                                                                        *
c                       m a p 1                                          *
c                                                                        *
c*************************************************************************
          include 'snake.inc'
        logical hsvz
        parameter (fracti=10.)
c                next line : begening of the /plot/ package
        logical lbox,luser,bare
        character unit*3,name*2
        common /cplot/unit(15),name(15)
        common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @              ,nplvold,h1(2,maxq),v1(2,maxq),nq
     @              ,ip,ih,iv,kp,bplot(3,maxplv)
     @              ,plobox(5,4,2,maxre),nface(maxre),lbox
     @              ,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @    ,iplpi(maxplv),iplpf(maxplv)
     @    ,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @		     ,luser,nclic,lclic(maxre),clight(4,maxre),
     @    epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3),nrs(maxre)
c                previous line : end of the /plot/ package
        character formh*15,formv*15
        real a
      common/ctext/formh,formv
      common/text/hldep,bh,vldep,bv,iexph,iexpv,ith,itv
      common/ntw/xmax1,xmin1,ymax1,ymin1,xnmin,xnmax,ynmin,ynmax
      external scale
c      write(10,*)'in map1'
c      write(10,*)xmin,ymin,xmax,ymax,n,xsize1,ysize1,zmarge
      xnmin=zmarge/29.7
      xnmax=(xsize1-zmarge)/29.7
      ynmin=zmarge/29.7
      ynmax=(ysize1-zmarge)/29.7
      hsvzn=(xnmax-xnmin)/(ynmax-ynmin)
      if(hsvz)then
		dx=xmax-xmin
		dy=ymax-ymin
        if(dx/dy.gt.hsvzn)then
			ymoy=(ymax+ymin)*.5
			dys2=.5*dx/hsvzn
			xmin1=xmin
			xmax1=xmax
			ymin1=ymoy-dys2
			ymax1=ymoy+dys2
         else
			xmoy=(xmax+xmin)*.5
			dxs2=.5*dy*hsvzn
			xmin1=xmoy-dxs2
			xmax1=xmoy+dxs2
			ymin1=ymin
			ymax1=ymax
         endif
      else
			xmin1=xmin
			xmax1=xmax
			ymin1=ymin
			ymax1=ymax
      endif
c      write(10,*)xmin1,xmax1,ymin1,ymax1
c      write(11,*)xnmin,xnmax,ynmin,ynmax
      call iswn(2,xmin1,xmax1,ymin1,ymax1)
      call isvp(2,xnmin,xnmax,ynmin,ynmax)
      cox=(xnmax-xnmin)/(xmax1-xmin1)
      shx=xmin1
      coy=(ynmax-ynmin)/(ymax1-ymin1)
      shy=ymin1
c       compute the scale and store the tic marks:
      dtic=(xmax1-xmin1)/fracti
      nq=0
c      write(6,*)'calling scale',xmin1,xmax1,dtic,a,bh,formh,ith,iexph	
      call scale(xmin1,xmax1,dtic,a,bh,formh,ith,iexph)
c      write(6,*)'return from scale',xmin1,xmax1,dtic,a,bh,formh,ith,iexph
      hldep=a+bh
      do 3 hl=hldep,xmax1,bh
                 if(nq.ge.maxq)goto5
                 nq=nq+1
                 h1(1,nq)=hl
                 v1(1,nq)=ymin1
                 h1(2,nq)=hl
 3                v1(2,nq)=ymax1
      dtic=(ymax1-ymin1)*hsvzn/fracti
        call scale(ymin1,ymax1,dtic,a,bv,formv,itv,iexpv)
      vldep=a+bv
      do 4 vl=vldep,ymax1,bv
                 if(nq.ge.maxq)goto 5
                 nq=nq+1
                 h1(1,nq)=xmin1
                 v1(1,nq)=vl
                 h1(2,nq)=xmax1
 4               v1(2,nq)=vl
c 5    save /cplot/,/plot/,/ctext/,/text/,/ntw/  
 5    return
      end

        subroutine paint(nreg)
c*************************************************************************
c                                                                        *
c                       p a i n t                                        *
c                                                                        *
c*************************************************************************
c
           include 'snake.inc'
c                next line : begening of the /plot/ package
        logical lbox,luser,bare
        character unit*3,name*2
        common /cplot/unit(15),name(15)
        common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @              ,nplvold,h1(2,maxq),v1(2,maxq),nq
     @              ,ip,ih,iv,kp,bplot(3,maxplv)
     @              ,plobox(5,4,2,maxre),nface(maxre),lbox
     @              ,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @    ,iplpi(maxplv),iplpf(maxplv)
     @    ,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @		     ,luser,nclic,lclic(maxre),clight(4,maxre),
     @    epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3),nrs(maxre)
c                previous line : end of the /plot/ package
      dimension icoll(3)
c
      call iselnt(2)
c  set:  style=solid
      call isln(1)
      call isfais(1)
c paint the 3 faces by filling them with background or foreground
c  color depending on the lightning:
      do 6 if=1,nface(nreg)
      if(clight(if,nreg).lt.0.)then
			icols=0
			icoll(if)=1
      else
			icols=1
			icoll(if)=0
      endif
      call isfaci(icols)
 6      call ifa(5,plobox(1,if,1,nreg),plobox(1,if,2,nreg))
c draw again the ridges of the 3 faces (erased by gfa) :
      do 7 if=1,nface(nreg)
        call isplci(icoll(if))
 7      call ipl(5,plobox(1,if,1,nreg),plobox(1,if,2,nreg))
      call isplci(1)
c?	call igterm
c      save /cplot/,/plot/
      return
      end
c
c
        subroutine cross(h,v)
c*************************************************************************
c                                                                        *
c                       c r o s s                                        *
c                                                                        *
c*************************************************************************
c         plot a cross at the position h*v (normalized coordinates)
c
        dimension hch(2),hcv(2),vch(2),vcv(2)
        hch(1)=h
        hch(2)=h
        hcv(1)=0.
        hcv(2)=21.
        vch(1)=0.
        vch(2)=28.5
        vcv(1)=v
        vcv(2)=v
        call isln(2)
        call ipl(2,hch(1),vch(1))
        call ipl(2,hcv(1),vcv(1))
        return
        end
c
        subroutine plotl(indve,work)
c*************************************************************************
c                                                                        *
c                       p l o t l                                        *
c                                                                        *
c*************************************************************************
c
c         for the trajectory of index=indve :
c         integrate the track-length from the data of plott
c        and put it in work(maxplp)
c
      include 'snake.inc'
c                next line : begening of the /plot/ package
      logical lbox,luser,bare
      character unit*3,name*2
      common /cplot/unit(15),name(15)
      common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @              ,nplvold,h1(2,maxq),v1(2,maxq),nq
     @              ,ip,ih,iv,kp,bplot(3,maxplv)
     @              ,plobox(5,4,2,maxre),nface(maxre),lbox
     @              ,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @    ,iplpi(maxplv),iplpf(maxplv)
     @    ,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu 
     @		     ,luser,nclic,lclic(maxre),clight(4,maxre),
     @    epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3),nrs(maxre)
c                previous line : end of the /plot/ package
      dimension work(maxplp)
c
      work(1)=0.
      do 1 i=2,nplp(indve)
            dl2=0.
            do 2 j=1,3
 2           dl2=dl2+(plott(i,indve,j)-plott(i-1,indve,j))**2
 1         work(i)=work(i-1)+sqrt(dl2)
c      save /cplot/,/plot/
      return
      end
      subroutine plota(indve,work)
c*************************************************************************
c                                                                        *
c                       p l o t a                                        *
c                                                                        *
c*************************************************************************
c
c         for the trajectory of index=indve :
c         compute the angular component (t) of cylindrical coordinates (r*t*z)
c        and put it in work(maxplp)
c        axis of the cylindrical coord. : z
c        origine of the angles on axis x    (t=arctan(y/x))
c        warning!: t is given in deg. , the starting point is taken
c        between -180. and 180. but the other points of the trajectory
c        can have their t outside this range .
c
      include 'snake.inc'
c                next line : begening of the /plot/ package
      logical lbox,luser,bare
      character unit*3,name*2
      common /cplot/unit(15),name(15)
      common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @              ,nplvold,h1(2,maxq),v1(2,maxq),nq
     @              ,ip,ih,iv,kp,bplot(3,maxplv)
     @              ,plobox(5,4,2,maxre),nface(maxre),lbox
     @              ,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @    ,iplpi(maxplv),iplpf(maxplv)
     @    ,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @		     ,luser,nclic,lclic(maxre),clight(4,maxre),
     @    epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3),nrs(maxre)
c                previous line : end of the /plot/ package
      dimension work(maxplp)
c
c        do1 i=1,nplp(indve)
c            arg=atan2(plott(i,indve,2),plott(i,indve,1))
c1           work(i)=arg*rtod
c        return
        scal(xn,yn,xa,ya)=xn*xa+yn*ya
        vect(xn,yn,xa,ya)=xn*ya-yn*xa
c
        x1=plott(1,indve,1)
        y1=plott(1,indve,2)
c         define starting angle in  ) -pi , pi ) :
        if(x1.eq.0..and.y1.eq.0.)then
                 work(1)=0.
        else
                 work(1)=atan2(y1,x1)*rtod
        endif
c        next code allows the angle to go outside ) -pi , pi ) :
       do 1 i=2,nplp(indve)
               x=plott(i,indve,1)
               y=plott(i,indve,2)
               if((x1.eq.0..and.y1.eq.0.).or.(x.eq.0..and.y.eq.0.))then
                 da=0.
               else
                 da=atan2(vect(x1,y1,x,y),scal(x1,y1,x,y))
               endif
               work(i)=work(i-1)+da*rtod
               x1=x
 1        y1=y
c      save /cplot/,/plot/
      return
      end
c
      subroutine plotra(indve,work)
c*************************************************************************
c                                                                        *
c                       p l o t r a                                      *
c                                                                        *
c*************************************************************************
c
c         for the trajectory of index=indve :
c         compute the angular component (t) by calling plota,
c		  the radius (r) 
c        and put their product in work(maxplp)
c        axis of the cylindrical coord. : z
c        origine of the angles on axis x    (t=arctan(y/x))
c        warning!: t is given in deg. , the starting point is taken
c        between -180. and 180. but the other points of the trajectory
c        can have their t outside this range .
c
      include 'snake.inc'
c                next line : begening of the /plot/ package
      logical lbox,luser,bare
      character unit*3,name*2
      common /cplot/unit(15),name(15)
      common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @              ,nplvold,h1(2,maxq),v1(2,maxq),nq
     @              ,ip,ih,iv,kp,bplot(3,maxplv)
     @              ,plobox(5,4,2,maxre),nface(maxre),lbox
     @              ,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @    ,iplpi(maxplv),iplpf(maxplv)
     @    ,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @		     ,luser,nclic,lclic(maxre),clight(4,maxre),
     @    epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3),nrs(maxre)
c                previous line : end of the /plot/ package
      dimension work(maxplp)
      call plota(indve,work)
      do 1 i=1,nplp(indve)
                 b2=0.
            do 2 j=1,2
 2             b2=b2+plott(i,indve,j)**2
 1      work(i)=work(i)*dtor*sqrt(b2)
c      save /cplot/,/plot/
      return
      end
      subroutine plotg(indve,iq,work)
c*************************************************************************
c                                                                        *
c                       p l o t g                                        *
c                                                                        *
c*************************************************************************
c
c         for the trajectory of index=indve :
c         copy plott(i,indve,iq) in work(i)
c
      include 'snake.inc'
c                next line : begening of the /plot/ package
      logical lbox,luser,bare
      character unit*3,name*2
      common /cplot/unit(15),name(15)
      common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @              ,nplvold,h1(2,maxq),v1(2,maxq),nq
     @              ,ip,ih,iv,kp,bplot(3,maxplv)
     @              ,plobox(5,4,2,maxre),nface(maxre),lbox
     @              ,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @    ,iplpi(maxplv),iplpf(maxplv)
     @    ,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @		     ,luser,nclic,lclic(maxre),clight(4,maxre),
     @    epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3),nrs(maxre)
c                previous line : end of the /plot/ package
      dimension work(maxplp)
c
      do 1 i=1,nplp(indve)
 1      work(i)=plott(i,indve,iq)
c      save /cplot/,/plot/
      return
      end
c
      subroutine plotm(indve,iq,work)
c*************************************************************************
c                                                                        *
c                       p l o t m                                        *
c                                                                        *
c*************************************************************************
c
c         for the trajectory of index=indve :
c         compute the module of the field or the radius from the data of plott
c        and put it in work(maxplp)
c
      include 'snake.inc'
c                next line : begening of the /plot/ package
      logical lbox,luser,bare
      character unit*3,name*2
      common /cplot/unit(15),name(15)
      common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @              ,nplvold,h1(2,maxq),v1(2,maxq),nq
     @              ,ip,ih,iv,kp,bplot(3,maxplv)
     @              ,plobox(5,4,2,maxre),nface(maxre),lbox
     @              ,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @    ,iplpi(maxplv),iplpf(maxplv)
     @    ,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @		     ,luser,nclic,lclic(maxre),clight(4,maxre),
     @    epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3),nrs(maxre)
c                previous line : end of the /plot/ package
      dimension work(maxplp)
c
      if(iq.eq.11)then
                j1=4
                j2=6
      else if(iq.eq.12)then
                j1=1
                j2=2
      else if(iq.eq.15)then
		j1=7
		j2=9
      else
		stop'plotm: invalid "iq" argument'
        endif
       do 1 i=1,nplp(indve)
                 b2=0.
                 do2 j=j1,j2
 2            b2=b2+plott(i,indve,j)**2
 1       work(i)=sqrt(b2)
c      save /cplot/,/plot/
      return
      end
c
      subroutine traj(xp,xm,yp,ym)
c*************************************************************************
c                                                                        *
c                       t r a j                                          *
c                                                                        *
c*************************************************************************
c
c                plot trajectories,tic marks ( and free-boxes if necessary)
c
c       input:
c             nure:number of region(s)  defined
c             nplv<=maxplv:number of trajectories to be plotted
c             nplp(traj.index)<=maxplp:number of submits in the plot-line
c             plott(submit index,traj.index,x y z bx by bz...
c                    ...px py pz):absol.coord.of submit,
c				  absol. field components and
c				  absol. spin components at this point
c             nq<=maxq:number of tic marks to be plotted
c             h1(submit index,tic mark index):absol.horiz.coord of the submit
c                of the tic mark
c             v1(submit index,tic mark index):absol.vertic.coord of the submit
c                of the tic mark
c             ip:index of the absol. axis(1:x,2:y,3:z)which is orthog./screen
c             ih:       "       "           "           "      horizontal "
c             iv:       "       "           "           "      vertical   "
c             plobox(submit index,horiz. or vertic.,region index):abs.coord
c                of submit of free-box(special submit list for 3d projection)
      include 'snake.inc'
c                next line : begening of the /plot/ package
      logical lbox,luser,bare
      character unit*3,name*2
      common /cplot/unit(15),name(15)
      common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @              ,nplvold,h1(2,maxq),v1(2,maxq),nq
     @              ,ip,ih,iv,kp,bplot(3,maxplv)
     @              ,plobox(5,4,2,maxre),nface(maxre),lbox
     @              ,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @		    ,iplpi(maxplv),iplpf(maxplv)
     @    ,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @		     ,luser,nclic,lclic(maxre),clight(4,maxre),
     @    epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3),nrs(maxre)
c                previous line : end of the /plot/ package
      common /absfra/ vea(maxve,0:maxep,21),nuve,nuep,nure
      dimension workh(maxplp),workv(maxplp)
c	dimension dist(2)
c
c                trace vectors:
c
      call isln(1)
      do 1 i=1,nplv
            if(ih.eq.10)then
            call plotl(i,workh)
            elseif(ih.eq.11.or.ih.eq.12.or.ih.eq.15)then
            call plotm(i,ih,workh)
            elseif(ih.eq.13)then
            call plota(i,workh)
            elseif(ih.eq.14)then
            call plotra(i,workh)
            else
            call plotg(i,ih,workh)
            endif
      if(iv.eq.10)then
      call plotl(i,workv)
      elseif(iv.eq.11.or.iv.eq.12.or.iv.eq.15)then
      call plotm(i,iv,workv)
      elseif(iv.eq.13)then
      call plota(i,workv)
      elseif(iv.eq.14)then
      call plotra(i,workv)
      else
      call plotg(i,iv,workv)
      endif
c next 6 temp:
c	eps=1.e-3
c	dist(1)=workh(iplpf(i))-workh(iplpi(i))
c	dist(2)=workv(iplpf(i))-workv(iplpi(i))
c	if(abs(dist(1)).gt.eps.or.abs(dist(2)).gt.eps)then
c		write(6,*)'vec#',i,' dist h/v=',dist
c	endif
      iplp=iplpf(i)-iplpi(i)+1
      if(iplp.gt.1)then
c	call ipl(iplp,workh,workv)
c  HIGZ est limite a 1000 points:
	    iplpc=iplp
	    iplpini=iplpi(i)
 200	    if(iplpc.gt.1000)then
      call ipl1(1000,workh(iplpini),workv(iplpini),
     s			xp,xm,yp,ym)
		iplpini=iplpini+999
		iplpc=iplpc-999
      goto 200
      endif
      call ipl1(iplpc,workh(iplpini),workv(iplpini),xp,xm,yp,ym)
c next temp:
cq
      call igterm
      endif
 1      continue
c
c                trace user plot:
c
      call isln(1)
      if(lbox.and.luser)then
         do i=1,nulu
         call ipl1(nuseu(i),plotus(1,i,ih),plotus(1,i,iv),
     s			xp,xm,yp,ym)
c next temp:
c	  call igterm
      end do
      endif
c
c                trace free-boxes:
c			first : not clickd boxes
      call isln(1)
      if(lbox.and..not.bare)then
        do 3 i=1,nure
        do 4 ic=1,nclic
 4      if(i.eq.lclic(ic))goto 3
      do 6 if=1,nface(i)	
         call ipl1(5,plobox(1,if,1,i),plobox(1,if,2,i),
     s			xp,xm,yp,ym)
c next temp:
c	  call igterm
 6      continue
 3      continue
c			second : clicked boxes,in chronological order

      do 5 i=1,nclic
      nreg=lclic(i)
      call paint(nreg)
c next temp:
c	  call igterm
 5       continue
c
c	plot region frames:
c
c isplci(color): 1=background, 2=red, 3=green, 4=dark blue, 5=yellow,
c		6=magenta (red-purple), 7=cyan (light blue)
c isln(type): 1=solid, 2=dashed, 3=dotted, 4=dashed-dotted
      call isplci(2)
      do ir=1,nure
       do irs=1,nrs(ir)
       call ipl1(2,regp(1,irs,ir,ih),regp(1,irs,ir,iv),
     s			xp,xm,yp,ym)
       enddo
      enddo
      call isplci(1)
c
c	plot end-planes:
c
c isplci(color): 1=background, 2=red, 3=green, 4=dark blue, 5=yellow,
c		6=magenta (red-purple), 7=cyan (light blue)
c isln(type): 1=solid, 2=dashed, 3=dotted, 4=dashed-dotted
      call isplci(4)
      do iep=1,nuep
        do ieps=1,neps(iep)
        call ipl1(2,epp(1,ieps,iep,ih),epp(1,ieps,iep,iv),
     s			xp,xm,yp,ym)
        enddo
      enddo
      call isplci(1)
       endif
c      save /cplot/,/plot/,/absfra/
      return
      end
c
      subroutine grad
c*************************************************************************
c                                                                        *
c                       g r a d                                          *
c                                                                        *
c*************************************************************************
c         plot numerical values of cordinates corresponding to the tic marks
c           for simplicity the plot is done in gks's normalized coordinates
c       input:
c             /text/:data to plot strings
c
      include 'snake.inc'
c                next line : begening of the /plot/ package
      logical lbox,luser,bare
      character unit*3,name*2
      common /cplot/unit(15),name(15)
      common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @              ,nplvold,h1(2,maxq),v1(2,maxq),nq
     @              ,ip,ih,iv,kp,bplot(3,maxplv)
     @              ,plobox(5,4,2,maxre),nface(maxre),lbox
     @              ,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @		    ,iplpi(maxplv),iplpf(maxplv)
     @    ,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @		     ,luser,nclic,lclic(maxre),clight(4,maxre),
     @    epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3),nrs(maxre)
c                previous line : end of the /plot/ package
      character formh*15,formv*15,string*20
      common/ctext/formh,formv
      common/text/hldep,bh,vldep,bv,iexph,iexpv,ith,itv
      common/ntw/xmax1,xmin1,ymax1,ymin1,xnmin,xnmax,ynmin,ynmax
      call ischh(1./120.)
      call ischup(-1.,0.)
c      call istxci(1)
      coef=10.**(-iexph)
      alp=(xnmax-xnmin)/(xmax1-xmin1)
        do 13 hl=hldep,xmax1,bh
                          hlco=hl*coef
                          write(string,formh)hlco
			  hnp=hn
			  hn=xnmin+hl*alp-xmin1*alp
           write(6,*)'call itx 1',hn,ynmin,string
                          call itx(hn,ynmin,string)
 13      continue
      call ischup(0.,1.)
      coef=10.**(-iexpv)
      alp=(ynmax-ynmin)/(ymax1-ymin1)
      do 16 vl=vldep,ymax1,bv
                          vlco=vl*coef
                          write(string,formv)vlco
			  vnp=vn
			  vn=ynmin+vl*alp-ymin1*alp
           write(6,*)'call itx 2',xnmin,vn,string
                call itx(xnmin,vn,string)
 16      continue
      call ischh(1./60.)
      vn=(vnp*5.+vn*3.)/8.
      call itx(xnmin,vn,name(iv)//'-'//unit(iv))
      call ischup(-1.,0.)
      hn=(hnp*3.+hn*5.)/8.
      call itx(hn,ynmin,name(ih)//'-'//unit(ih))
      call ischup(0.,1.)
      call ischh(1./3.)
c      save /cplot/,/plot/,/ctext/,/text/,/ntw/
      return
      end
c
c
c
      subroutine ipl1(n,x,y,xp,xm,yp,ym)
c*************************************************************************
c                                                                        *
c                       i p l 1                                          *
c                                                                        *
c*************************************************************************
c 
c    La routine ipl de Higz genere un trace incorrect lorsque la polyline
c  a l'un de ses sommets situe a grande distance de la fenetre effectivement
c  plottee (environ 38 en coordonnee normalisee pour la version 93c de cernlib).
c    ipl1(n,x,y,xp,xm,yp,ym) se substitue a ipl(n,x,y) en assurant un
c  ecretage sur le rectangle (xp,xm,yp,ym) (coord. world) suppose 
c  correspondre a x=+-30, y=+-30 (coord. normalisees).
c    La polyline initiale est decoupee en sous-polylines dont seuls le point
c  origine et le point extremite sortent eventuellement du cadre. Ces 2 points
c  sont alors rabatus sur le cadre par interpolation lineaire, vers le 
c  point suivant pour le premier, vers le point precedent pour le dernier.
c  Les parties de polylines entierement a l'exterieur du cadre sont converties
c  en points de plot (polyline a 2 sommets coincidents) et ne genent par
c  le plot.
c modif P.V. 9/10/95: Les parties de polylines entierement a l'exterieur
c  du cadre ne donnent lieu a aucun plot
      parameter (maxpoly=1000)
      dimension x(*),y(*),xt(maxpoly),yt(maxpoly)
      logical previous_in
      if(n.gt.maxpoly)then
      write(6,*)'ipl1: polyline size limited to ',maxpoly,'!'
      stop
      endif
      do ii=1,n
      if(x(ii).lt.xp.and.x(ii).gt.xm.and.
     s		   y(ii).lt.yp.and.y(ii).gt.ym)then
            if(ii.eq.1)then
				if=1
				xt(if)=x(ii)
				yt(if)=y(ii)
            elseif(previous_in)then
				if=if+1
				xt(if)=x(ii)
				yt(if)=y(ii)
				if(ii.eq.n)then
				call ipl(if,xt,yt)
c               write(10,*)if,xt,yt
            endif
            else
				if=1
                call rabat(ii-1,ii,if,x,y,xt,yt,
     s					xp,xm,yp,ym)
				if=if+1
				xt(if)=x(ii)
				yt(if)=y(ii)
                if(ii.eq.n)then
                call ipl(if,xt,yt)
              endif
            endif
            previous_in=.true.
      else
      if(ii.ne.1.and.previous_in)then
				if=if+1
            call rabat(ii,ii-1,if,x,y,xt,yt,
     s					xp,xm,yp,ym)
            call ipl(if,xt,yt)
            elseif(ii.ne.1)then
				den=x(ii)-x(ii-1)
				xl1=2.
                if(den.ne.0.)xl1=(xp-x(ii-1))/den
				y1=xl1*y(ii)+(1.-xl1)*y(ii-1)
				xl2=2.
                if(den.ne.0.)xl2=(xm-x(ii-1))/den
				y2=xl2*y(ii)+(1.-xl2)*y(ii-1)
				den=y(ii)-y(ii-1)
				xl3=2.
                if(den.ne.0.)xl3=(yp-y(ii-1))/den
				x1=xl3*x(ii)+(1.-xl3)*x(ii-1)
c si les points ii et ii-1 sont en dehors, test si le segment (ii->ii-1)
c intersect la fenetre en calculant ses intersections eventuelles aves
c 3 des 4 cotes de la fenetre:
                if((xl1.gt.0..and.xl1.lt.1..and.
     s				   y1.gt.ym.and.y1.lt.yp).or.
     s				   (xl2.gt.0..and.xl2.lt.1..and.
     s				   y2.gt.ym.and.y2.lt.yp).or.
     s				   (xl3.gt.0..and.xl3.lt.1..and.
     s				   x1.gt.xm.and.x1.lt.xp))then
                  if=1
                  call rabat(ii-1,ii,if,x,y,xt,yt,
     s					xp,xm,yp,ym)
                        if=if+1
                  call rabat(ii,ii-1,if,x,y,xt,yt,
     s					xp,xm,yp,ym)
                  call ipl(if,xt,yt)
                  endif
            endif
            previous_in=.false.
            endif
      enddo
      return
      end
c
c
	subroutine ipl1old(n,x,y,xp,xm,yp,ym)
c*************************************************************************
c                                                                        *
c                       i p l 1                                          *
c                                                                        *
c*************************************************************************
c 
c    La routine ipl de Higz genere un trace incorrect lorsque la polyline
c  a l'un de ses sommets situe a grande distance de la fenetre effectivement
c  plottee (environ 38 en coordonnee normalisee pour la version 93c de cernlib).
c    ipl1(n,x,y,xp,xm,yp,ym) se substitue a ipl(n,x,y) en assurant un
c  ecretage sur le rectangle (xp,xm,yp,ym) (coord. world) suppose 
c  correspondre a x=+-30, y=+-30 (coord. normalisees).
c    La polyline initiale est decoupee en sous-polylines dont seuls le point
c  origine et le point extremite sortent eventuellement du cadre. Ces 2 points
c  sont alors rabatus sur le cadre par interpolation lineaire, vers le 
c  point suivant pour le premier, vers le point precedent pour le dernier.
c  Les parties de polylines entierement a l'exterieur du cadre sont converties
c  en points de plot (polyline a 2 sommets coincidents) et ne genent par
c  le plot.
      parameter (maxpoly=1000)
      dimension x(*),y(*),xt(maxpoly),yt(maxpoly)
      if(n.gt.maxpoly)then
            write(6,*)'ipl1: polyline size limited to ',maxpoly,'!'
            stop
      endif
      idep=0
 2      if(x(idep+1).gt.xp.or.x(idep+1).lt.xm.or.
     s	y(idep+1).gt.yp.or.y(idep+1).lt.ym)then
      call rabat(idep+1,idep+2,1,x,y,xt,yt,xp,xm,yp,ym)
      else
        xt(1)=x(idep+1)
        yt(1)=y(idep+1)
      endif
      do 1 i=idep+2,n
		ic=i
      if(x(i).gt.xp.or.x(i).lt.xm.or.
     s		y(i).gt.yp.or.y(i).lt.ym)then
            call rabat(i,i-1,i-idep,
     s  		x,y,xt,yt,xp,xm,yp,ym)
            call ipl(i-idep,xt,yt)
            if(i.eq.n)return
			idep=i-1
            goto 2
        else
			xt(i-idep)=x(i)
			yt(i-idep)=y(i)
      endif
 1      continue
      call ipl(ic-idep,xt,yt)
      return
      end
c
c*************************************************************************
c                                                                        *
c                       r a b a t                                        *
c                                                                        *
c*************************************************************************
c 
c  copie dans (xt,yt) a l'adresse j le point (x,y), adresse i, ecrete par
c  le rectangle (xp,xm,yp,ym). L'interpolation est faite de (x,y)(i) vers
c  (x,y)(ip) avec ip=i+-1. Le point (x,y)(ip) est suppose situe de sorte 
c  que le segment i-->ip intercepte le rectangle (xp,xm,yp,ym)
c
      subroutine rabat(i,ip,j,x,y,xt,yt,xp,xm,yp,ym)
      dimension x(*),y(*),xt(*),yt(*)
      xc=x(i)
      yc=y(i)
      if(xc.gt.xp)then
		yc=(y(ip)*(xc-xp)+yc*(xp-x(ip)))/(xc-x(ip))
		xc=xp
      elseif(xc.lt.xm)then
		yc=(y(ip)*(xc-xm)+yc*(xm-x(ip)))/(xc-x(ip))
		xc=xm
      endif
      if(yc.gt.yp)then
		xc=(x(ip)*(yc-yp)+xc*(yp-y(ip)))/(yc-y(ip))
		yc=yp
      elseif(yc.lt.ym)then
		xc=(x(ip)*(yc-ym)+xc*(ym-y(ip)))/(yc-y(ip))
		yc=ym
      endif
      xt(j)=xc
      yt(j)=yc
      return
      end
c
      subroutine scale(ba,bb,fn,a,fnsr,form,it,iexp)
c*************************************************************************
c                                                                        *
c                       s c a l e                                        *
c                                                                        *
c*************************************************************************
c
c                find the most economic way to mark from ba to bb
c       with an increment aproximately=fn
c
c       input:
c          ba,bb and fn
c       output:
c          a(near ba) and fnsr(near fn)
c          form:most economical write format
c          it:number of chatacter(s) in form
c          iexp:value of the decimal exponent (if any)
c
      character form*15,formm*8,forme*7,formf*17
       b1=ba
       b2=bb
      if(b1.eq.b2)then
c            in this case:a and b arbitrary:
      a=b1
      b=1.
      return
      endif
      call c125(abs(fn),b,imin)
       fnsr=sign(b,b2-b1)
       ap=abs(b1)
       ap=ap-amod(ap,b)
       a=sign(ap,b1)
c  find first and last usefull digit:
      if(a.ne.0.)then
        i1max=ifix(alog10(abs(a))+.0001)
        if(abs(a).lt.1.)i1max=i1max-1
        else
        i1max=-1
        endif
      if(b2.ne.0.)then
        i2max=ifix(alog10(abs(b2))+.0001)
        if(abs(b2).lt.1.)i2max=i2max-1
        else
        i2max=-1
      endif
       imax=max(i1max,i2max)
        imoy=(imin+imax)/2+1
        iexp=3*(imoy/3)
      if(imin.le.2.and.imax.ge.-3)iexp=0
       imin=imin-iexp
       imax=imax-iexp
       ifmin=min(imin,0)
       ifmax=max(imax,-1)
       it=ifmax-ifmin+2
      if(b1.lt.0..or.b2.lt.0.)it=it+1
       if=-ifmin
       write(formm,100)it,if
 100      format('(f',i2,'.',i2)
        forme=')'
      if(iexp.ne.0)then
        nexp=ifix(alog10(abs(float(iexp))+.0001))+1
        if(iexp.lt.0)nexp=nexp+1
        write(formf,200)nexp
 200     format('('',''''e'',i',i1,','''''')'')')
        write(forme,formf)iexp
        it=it+1+nexp
        endif
        form=formm//forme
       return
       end
       subroutine c125(x,xp,imin)
c*************************************************************************
c                                                                        *
c                       c 1 2 5                                          *
c                                                                        *
c*************************************************************************
c
c                for x>0. find xp as close as possible from x
c       such that: xp=(1. , 2. , 2.5  or 5.)*10.**n
c
c       input: x
c       output:
c              xp
c              imin:distance in digit between the last stginficant digit
c                         and the decimal point( >0 at right)
c
       dimension fl(5),xl(5),id(5)
       data fl/1.,2.,2.5,5.,10./
       data id/0 ,0 ,-1 ,0 , 1 /
       do 18 i=1,5
         xl(i)=alog10(fl(i))
 18   continue
       xlc=alog10(x)
       ixl=int(xlc+.0001)
       if(x.lt.1.)ixl=ixl-1
       dixl=10.**ixl
       xl0=xlc-ixl
       do 21 i=2,5
         if(xl0.lt.xl(i))then
                 xp=fl(i)*dixl
                 goto1
         endif
 21   continue
1      xm=fl(i-1)*dixl
       imin=id(i)+ixl
       if(xp/x.gt.x/xm)then
        xp=xm
        imin=id(i-1)+ixl
       endif
       return
       end
c
c*************************************************************************
c                                                                        *
c                       i t r a n s                                      *
c                                                                        *
c*************************************************************************
c
        subroutine itrans
      include 'snake.inc'
c                next line : begening of the /relfra/ package
      character epsw*4,rname*8,rnamer*8,title*80,cerfiln*20
      common /crelfr/epsw(maxre),rname(maxre),rnamer(maxre)
     @                 ,title(maxre),cerfiln(maxepr,maxre)
      common /relfra/ ver(maxve,12),xyz0(maxre,3)
     @                 ,frbox(2,maxre,3)
     @                 ,zang(maxre),xang(maxre),yang(maxre)
     @                 ,indve,indre,indep,ilive,adata(maxdat,maxre)
     @                 ,indrer(maxre),tept(maxepr,maxre)
     @                 ,yepmin(maxre),yepmax(maxre),kept(maxepr,maxre)
     @                 ,yepstp(maxre),yept(maxepr,maxre),nep(maxre)
     @                 ,colli(maxepr,maxre,6)
c                previous line : end of the /relfra/ package

      call plorot(0.,0.,0.,x0,y0,z0,xb,yb,zb)
      call plorot(1.,0.,0.,xb,yb,zb,xx,xy,xz)
      call plorot(0.,1.,0.,xb,yb,zb,yx,yy,yz)
      call plorot(0.,0.,1.,xb,yb,zb,zx,zy,zz)
        irsave=indre
      do 3 i=1,indre-1
                 if(rnamer(indre).eq.rname(i))then
                          indrer(indre)=i
                          goto4
                 endif
 3      continue
        write(6,*)' can''t find this region :',rnamer(indre)
        stop
 4      indre=indrer(indre)
        call plorot(x0,y0,z0,x0,y0,z0,xb,yb,zb)
        call plorot(xx,xy,xz,xb,yb,zb,xx,xy,xz)
        call plorot(yx,yy,yz,xb,yb,zb,yx,yy,yz)
        call plorot(zx,zy,zz,xb,yb,zb,zx,zy,zz)
        indre=irsave
        xyz0(indre,1)=x0
        xyz0(indre,2)=y0
        xyz0(indre,3)=z0
        xang(indre)=asin(yz)
        if((xz.eq.0.and.zz.eq.0.).or.(yx.eq.0..and.yy.eq.0.))then
c         polar ambiguity:
                 yang(indre)=0.
                 zang(indre)=atan2(xy,xx)
        else
                 yang(indre)=atan2(-xz,zz)
                 zang(indre)=atan2(-yx,yy)
        endif
        zad=zang(indre)*rtod
        xad=xang(indre)*rtod
        yad=yang(indre)*rtod
c      save /crelfr/,/relfra/
		return
      end
      subroutine rot
c*************************************************************************
c                                                                        *
c                       r o t                                            *
c                                                                        *
c*************************************************************************
c
c                at each region convert absolute coord. in relative one
c                to fill the current relative coord. buffer
c       input:
c             nuve:number of input vector(s) defined
c             vea(iv,iep,ic):description of the particules at each end-plane
c                iv<=nuve<=maxve:vector index
c                iep<=nuep<=maxep:end-plane index
c                ic=1:x coord. in absolute frame
c                   2:y            "
c                   3:z            "
c                   4:x component of the normalized mommentum in absol.frame
c                   5:y            "
c                   6:z            "
c                   7:module of the mommentum
c                   8:trace length measured from the "spring point"
c                   9:life-flag(-1.:buried,0.:just dead,1.:alive)
c            10 to 15:the same as 1 to 6 but in relative frame coordinates
c		   16:x component of the spin
c		   17:y          "
c		   18:z          "
c	     19 to 21: the same as 16 to 18 but in relative frame coord.
c             indre=current region   index
c             indep=   "    end-plane "
c             xyz0(region index,x y or z):absolute coord. of the rel.origine
c             zang(region index),xang(region index) and yang(region index):
c                rotation angles to define the relat.frame in the absol.frame
c             adata(data index < maxdat,region index):additional data to
c                describe a region
c       output:
c             ver(iv,ic):description of the particules in the current region
c                iv<=nuve<=maxve:vector index
c                ic=1:x coord.relative to the current region(in relative frame)
c                   2:y            "
c                   3:z            "
c                   4:x component of the normalized mommentum(in relat.frame)
c                   5:y            "
c                   6:z            "
c                   7:module of the mommentum
c                   8:trace length measured from the "spring point"
c                   9:life-flag(-1.:buried,0.:just dead,1.:alive)
c		   10:x component of the spin (in relat.frame)
c		   11:y          "
c		   12:z          "
c       output for subroutine rot and plorot:
c             x0,y0,z0,stx,ctx,stz,ctz,sty,cty=rotation matrix elements
c
      include 'snake.inc'
c                next line : begening of the /plot/ package
      logical lbox,luser,bare
      character unit*3,name*2
      common /crot/x0,y0,z0,stx,ctx,stz,ctz,sty,cty
      common /cplot/unit(15),name(15)
      common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @              ,nplvold,h1(2,maxq),v1(2,maxq),nq
     @              ,ip,ih,iv,kp,bplot(3,maxplv)
     @              ,plobox(5,4,2,maxre),nface(maxre),lbox
     @              ,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @    ,iplpi(maxplv),iplpf(maxplv)
     @    ,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @		     ,luser,nclic,lclic(maxre),clight(4,maxre),
     @    epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3),nrs(maxre)
c                previous line : end of the /plot/ package
      common /absfra/ vea(maxve,0:maxep,21),nuve,nuep,nure
      common/regin/veri(maxve,maxre,12)
c                next line : begening of the /relfra/ package
      character epsw*4,rname*8,rnamer*8,title*80,cerfiln*20
      common /crelfr/epsw(maxre),rname(maxre),rnamer(maxre)
     @                 ,title(maxre),cerfiln(maxepr,maxre)
      common /relfra/ ver(maxve,12),xyz0(maxre,3)
     @                 ,frbox(2,maxre,3)
     @                 ,zang(maxre),xang(maxre),yang(maxre)
     @                 ,indve,indre,indep,ilive,adata(maxdat,maxre)
     @                 ,indrer(maxre),tept(maxepr,maxre)
     @                 ,yepmin(maxre),yepmax(maxre),kept(maxepr,maxre)
     @                 ,yepstp(maxre),yept(maxepr,maxre),nep(maxre)
     @                 ,colli(maxepr,maxre,6)
c                previous line : end of the /relfra/ package
c                next line : begening of the /fidef/ package
      logical cyl,lfil
      character*80 fifi
      common /fidef/method(maxre),itype(maxre),indfi(maxre),fact(maxre)
     @		,cyl(maxre),lfil
      common /cfidef/ fifi(maxre)
c                previous line : end of the /fidef/ package
c
        x0=xyz0(indre,1)
        y0=xyz0(indre,2)
        z0=xyz0(indre,3)
          stz=sin(zang(indre))
          ctz=cos(zang(indre))
          stx=sin(xang(indre))
          ctx=cos(xang(indre))
          sty=sin(yang(indre))
          cty=cos(yang(indre))
c
c
        do 1 i=1,nuve
c                if this vector is dead or buried, no rotation:
                 if(vea(i,indep,9).le.0.)then
                          ver(i,9)=vea(i,indep,9)
		 	  do 3 iq=1,12
3				veri(i,indre,iq)=0.
			  veri(i,indre,9)=vea(i,indep,9)
                          goto1
                 endif
                 if(method(indre).le.2)then
c                 no field in the new region:
c			plot these points to zero the field:
c 19/11/92 next not usefull for the first region: this plot has been done 
c		by inject 
		  if(indre.ne.1)call plotin(i,vea(i,indep,1),
     @				vea(i,indep,2),
     @				vea(i,indep,3),
     @				0.,0.,0.,vea(i,indep,16),
     @  			vea(i,indep,17),
     @				vea(i,indep,18),.false.)
                 endif
c                first:translation (the new origine is the point "0"):
                 x=vea(i,indep,1)-x0
                 y=vea(i,indep,2)-y0
                 z=vea(i,indep,3)-z0
                 cx=vea(i,indep,4)
                 cy=vea(i,indep,5)
                 cz=vea(i,indep,6)
                 px=vea(i,indep,16)
                 py=vea(i,indep,17)
                 pz=vea(i,indep,18)
c                second:rotations if necessary:
                 if(zang(indre).ne.0.)then
c                         in this case:rotation/z axis
                          xt=x*ctz+y*stz
                          y=y*ctz-x*stz
                          x=xt
                          cxt=cx*ctz+cy*stz
                          cy=cy*ctz-cx*stz
                          cx=cxt
                          pxt=px*ctz+py*stz
                          py=py*ctz-px*stz
                          px=pxt
                 endif
                 if(xang(indre).ne.0.)then
c                         in this case:rotation/x axis
                          yt=y*ctx+z*stx
                          z=z*ctx-y*stx
                          y=yt
                          cyt=cy*ctx+cz*stx
                          cz=cz*ctx-cy*stx
                          cy=cyt
                          pyt=py*ctz+pz*stz
                          pz=pz*ctz-py*stz
                          py=pyt
                 endif
                 if(yang(indre).ne.0.)then
c                         in this case:rotation/y axis
                          zt=z*cty+x*sty
                          x=x*cty-z*sty
                          z=zt
                          czt=cz*cty+cx*sty
                          cx=cx*cty-cz*sty
                          cz=czt
                          pzt=pz*ctz+px*stz
                          px=px*ctz-pz*stz
                          pz=pzt
                 endif
                 ver(i,1)=x
                 ver(i,2)=y
                 ver(i,3)=z
                 ver(i,4)=cx
                 ver(i,5)=cy
                 ver(i,6)=cz
                 ver(i,7)=vea(i,indep,7)
                 ver(i,8)=vea(i,indep,8)
                 ver(i,9)=vea(i,indep,9)
                 ver(i,10)=px
                 ver(i,11)=py
                 ver(i,12)=pz
c  fill  common/regin/veri(maxve,maxre,12)= rel. coord. at region entrance
      do 2 iq=1,12
 2      veri(i,indre,iq)=ver(i,iq)
 1      continue
c      save /crot/,/cplot/,/plot/,/absfra/,/regin/,/crelfr/,/relfra/,
c     &     /fidef/,/cfidef/
      return
      end
c
c
      subroutine unrot
c*************************************************************************
c                                                                        *
c                       u n r o t                                        *
c                                                                        *
c*************************************************************************
c
c                at each end-plane convert relative coord. in absolute one
c                to fill the absolute coord. buffer and the plot buffer
c       input from subroutine rot:
c             x0,y0,z0,stx,ctx,stz,ctz,sty,cty=rotation matrix elements
c       input:
c             nuve:number of input vector(s) defined
c             indre=current region   index
c             indep=   "    end-plane "
c             ver(iv,ic):description of the particules in the current region
c                iv<=nuve<=maxve:vector index
c                ic=1:x coord.relative to the current region(in relative frame)
c                   2:y            "
c                   3:z            "
c                   4:x component of the normalized mommentum(in relat.frame)
c                   5:y            "
c                   6:z            "
c                   7:module of the mommentum
c                   8:trace length measured from the "spring point"
c                   9:life-flag(-1.:buried,0.:just dead,1.:alive)
c		   10:x component of the spin (in relat.frame)
c		   11:y          "
c		   12:z          "
c             xyz0(region index,x y or z):absolute coord. of the rel.origine
c             zang(region index),xang(region index) and yang(region index):
c                rotation angles to define the relat.frame in the absol.frame
c             adata(data index < maxdat,region index):additional data to
c                describe a region
c       output:
c             vea(iv,iep,ic):description of the particules at each end-plane
c                iv<=nuve<=maxve:vector index
c                iep<=nuep<=maxep:end-plane index
c                ic=1:x coord. in absolute frame
c                   2:y            "
c                   3:z            "
c                   4:x component of the normalized mommentum in absol.frame
c                   5:y            "
c                   6:z            "
c                   7:module of the mommentum
c                   8:trace length measured from the "spring point"
c                   9:life-flag(-1.:buried,0.:just dead,1.:alive)
c            10 to 15:the same as 1 to 6 but in relative frame coordinates
c		   16:x component of the spin
c		   17:y          "
c		   18:z          "
c	     19 to 21: the same as 16 to 18 but in relative frame coord.
c             ilive=   "    number of alive vector(s)
c             nplv<=maxplv:number of trajectories to be plotted
c             nplp(traj.index)<=maxplp:number of submits in the plot-line
c             plott(submit index,traj.index,x y z bx by bz...
c                    ...px py pz):absol.coord.of submit,
c				  absol. field components and
c				  absol. spin components at this point
c
      include 'snake.inc'
c                next line : begening of the /plot/ package
      logical lbox,luser,bare
      character unit*3,name*2
      common /crot/x0,y0,z0,stx,ctx,stz,ctz,sty,cty
      common /cplot/unit(15),name(15)
      common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @              ,nplvold,h1(2,maxq),v1(2,maxq),nq
     @              ,ip,ih,iv,kp,bplot(3,maxplv)
     @              ,plobox(5,4,2,maxre),nface(maxre),lbox
     @              ,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @    ,iplpi(maxplv),iplpf(maxplv)
     @    ,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @		     ,luser,nclic,lclic(maxre),clight(4,maxre),
     @    epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3),nrs(maxre)
c                previous line : end of the /plot/ package
      common /absfra/ vea(maxve,0:maxep,21),nuve,nuep,nure
c                next line : begening of the /relfra/ package
      character epsw*4,rname*8,rnamer*8,title*80,cerfiln*20
      common /crelfr/epsw(maxre),rname(maxre),rnamer(maxre)
     @                 ,title(maxre),cerfiln(maxepr,maxre)
      common /relfra/ ver(maxve,12),xyz0(maxre,3)
     @                 ,frbox(2,maxre,3)
     @                 ,zang(maxre),xang(maxre),yang(maxre)
     @                 ,indve,indre,indep,ilive,adata(maxdat,maxre)
     @                 ,indrer(maxre),tept(maxepr,maxre)
     @                 ,yepmin(maxre),yepmax(maxre),kept(maxepr,maxre)
     @                 ,yepstp(maxre),yept(maxepr,maxre),nep(maxre)
     @                 ,colli(maxepr,maxre,6)
c                previous line : end of the /relfra/ package
c                next line : begening of the /fidef/ package
      logical cyl,lfil
      character*80 fifi
      common /fidef/method(maxre),itype(maxre),indfi(maxre),fact(maxre)
     @			,cyl(maxre),lfil
      common /cfidef/ fifi(maxre)
c                previous line : end of the /fidef/ package
       indep=indep+1
        if(indep.gt.maxep)stop'unrot:max num. of end-planes surpas.'
c
c
c
        ilive=0
        do 2 i=1,nuve
c                if the vector is buried:no rotation
                 if(ver(i,9).le.-1.)then
            do 5 j=1,21
5                            vea(i,indep,j)=0.
                          vea(i,indep,9)=-1.
                          goto 2
                 endif
                 x=ver(i,1)
                 y=ver(i,2)
                 z=ver(i,3)
                 cx=ver(i,4)
                 cy=ver(i,5)
                 cz=ver(i,6)
                 px=ver(i,10)
                 py=ver(i,11)
                 pz=ver(i,12)
c                first:rotations if necessary:
                 if(yang(indre).ne.0.)then
c                         in this case:rotation/y axis
                          zt=z*cty-x*sty
                          x=x*cty+z*sty
                          z=zt
                          czt=cz*cty-cx*sty
                          cx=cx*cty+cz*sty
                          cz=czt
                          pzt=pz*cty-px*sty
                          px=px*cty+pz*sty
                          pz=pzt
                 endif
                 if(xang(indre).ne.0.)then
c                         in this case:rotation/x axis
                          yt=y*ctx-z*stx
                          z=z*ctx+y*stx
                          y=yt
                          cyt=cy*ctx-cz*stx
                          cz=cz*ctx+cy*stx
                          cy=cyt
                          pyt=py*ctx-pz*stx
                          pz=pz*ctx+py*stx
                          py=pyt
                 endif
                 if(zang(indre).ne.0.)then
c                         in this case:rotation/z axis
                          xt=x*ctz-y*stz
                          y=y*ctz+x*stz
                          x=xt
                          cxt=cx*ctz-cy*stz
                          cy=cy*ctz+cx*stz
                          cx=cxt
                          pxt=px*ctz-py*stz
                          py=py*ctz+px*stz
                          px=pxt
                 endif
c                second:translation:
                 vea(i,indep,1)=x+x0
                 vea(i,indep,2)=y+y0
                 vea(i,indep,3)=z+z0
                 vea(i,indep,4)=cx
                 vea(i,indep,5)=cy
                 vea(i,indep,6)=cz
                 vea(i,indep,7)=ver(i,7)
                 vea(i,indep,8)=ver(i,8)
                 vea(i,indep,9)=ver(i,9)
                 if(ver(i,9).gt.0.)ilive=ilive+1
                 vea(i,indep,10)=ver(i,1)
                 vea(i,indep,11)=ver(i,2)
                 vea(i,indep,12)=ver(i,3)
                 vea(i,indep,13)=ver(i,4)
                 vea(i,indep,14)=ver(i,5)
                 vea(i,indep,15)=ver(i,6)
                 vea(i,indep,16)=px
                 vea(i,indep,17)=py
                 vea(i,indep,18)=pz
                 vea(i,indep,19)=ver(i,10)
                 vea(i,indep,20)=ver(i,11)
                 vea(i,indep,21)=ver(i,12)
                 if(i.le.maxplv)then
c       plot these points to end properly the previous plot:
c	(the last field seen in rungk (at landing time)...
c			..has been saved in bplot(j,i))
           	    call plotin(i,vea(i,indep,1),
     @				vea(i,indep,2),
     @				vea(i,indep,3),
     @    		        bplot(1,i),
     @				bplot(2,i),
     @				bplot(3,i),
     @				vea(i,indep,16),
     @				vea(i,indep,17),
     @				vea(i,indep,18),
     @				.false.)
c	reset bplot:
                     do4j=1,3
 4                    bplot(j,i)=0.
        endif
 2       continue
c      save /crot/,/cplot/,/plot/,/absfra/,/crelfr/,/relfra/,
c	 &     /fidef/,/cfidef/
      return
      end
c
      subroutine plorot(xr,yr,zr,xa,ya,za,xb,yb,zb)
c*************************************************************************
c                                                                        *
c                       p l o r o t                                      *
c                                                                        *
c*************************************************************************
c
c       convert the relative coortdinates of a point (xr,yr,zr)
c            in the absolutes one (xa,ya,za)
c
c       input:
c             xr,yr,zr
c             xyz0(region index,x y or z):absolute coord. of the rel.origine
c             zang(region index),xang(region index) and yang(region index):
c                rotation angles to define the relat.frame in the absol.frame
c       output:
c              xa,ya,za
c
      include 'snake.inc'
c                next line : begening of the /relfra/ package
      character epsw*4,rname*8,rnamer*8,title*80,cerfiln*20
      common /crelfr/epsw(maxre),rname(maxre),rnamer(maxre)
     @                 ,title(maxre),cerfiln(maxepr,maxre)
      common /relfra/ ver(maxve,12),xyz0(maxre,3)
     @                 ,frbox(2,maxre,3)
     @                 ,zang(maxre),xang(maxre),yang(maxre)
     @                 ,indve,indre,indep,ilive,adata(maxdat,maxre)
     @                 ,indrer(maxre),tept(maxepr,maxre)
     @                 ,yepmin(maxre),yepmax(maxre),kept(maxepr,maxre)
     @                 ,yepstp(maxre),yept(maxepr,maxre),nep(maxre)
     @                 ,colli(maxepr,maxre,6)
c                previous line : end of the /relfra/ package
      data indrep/0/
      if(indre.ne.indrep)then
                 indrep=indre
                 x0=xyz0(indre,1)
                 y0=xyz0(indre,2)
                 z0=xyz0(indre,3)
                 stz=sin(zang(indre))
                 ctz=cos(zang(indre))
                 stx=sin(xang(indre))
                 ctx=cos(xang(indre))
                 sty=sin(yang(indre))
                 cty=cos(yang(indre))
        endif
                 x=xr
                 y=yr
                 z=zr
                 if(yang(indre).ne.0.)then
c                         in this case:rotation/y axis
                          zt=z*cty-x*sty
                          x=x*cty+z*sty
                          z=zt
                 endif
                 if(xang(indre).ne.0.)then
c                         in this case:rotation/x axis
                          yt=y*ctx-z*stx
                          z=z*ctx+y*stx
                          y=yt
                 endif
                 if(zang(indre).ne.0.)then
c                         in this case:rotation/z axis
                          xt=x*ctz-y*stz
                          y=y*ctz+x*stz
                          x=xt
                 endif
                 xa=x+x0
                 ya=y+y0
                 za=z+z0
                 xb=x
                 yb=y
                 zb=z
c      save /crelfr/,/relfra/
      return
      end
c
      subroutine unplorot(xa,ya,za,xr,yr,zr,xb,yb,zb)
c*************************************************************************
c                                                                        *
c                       u n p l o r o t                                  *
c                                                                        *
c*************************************************************************
c
c       convert the absolute coortdinates of a point (xr,yr,zr)
c            into the relative one (xa,ya,za)
c
c       input:
c              xa,ya,za
c             xyz0(region index,x y or z):absolute coord. of the rel.origine
c             zang(region index),xang(region index) and yang(region index):
c                rotation angles to define the relat.frame in the absol.frame
c       output:
c             xr,yr,zr
c
      include 'snake.inc'
c                next line : begening of the /relfra/ package
      character epsw*4,rname*8,rnamer*8,title*80,cerfiln*20
      common /crelfr/epsw(maxre),rname(maxre),rnamer(maxre)
     @                 ,title(maxre),cerfiln(maxepr,maxre)
      common /relfra/ ver(maxve,12),xyz0(maxre,3)
     @                 ,frbox(2,maxre,3)
     @                 ,zang(maxre),xang(maxre),yang(maxre)
     @                 ,indve,indre,indep,ilive,adata(maxdat,maxre)
     @                 ,indrer(maxre),tept(maxepr,maxre)
     @                 ,yepmin(maxre),yepmax(maxre),kept(maxepr,maxre)
     @                 ,yepstp(maxre),yept(maxepr,maxre),nep(maxre)
     @                 ,colli(maxepr,maxre,6)
c                previous line : end of the /relfra/ package
      data indrep/0/
        if(indre.ne.indrep)then
                 indrep=indre
                 x0=xyz0(indre,1)
                 y0=xyz0(indre,2)
                 z0=xyz0(indre,3)
                 stz=sin(zang(indre))
                 ctz=cos(zang(indre))
                 stx=sin(xang(indre))
                 ctx=cos(xang(indre))
                 sty=sin(yang(indre))
                 cty=cos(yang(indre))
        endif
                 x=xa-x0
                 y=ya-y0
                 z=za-z0
                 if(zang(indre).ne.0.)then
c                         in this case:rotation/z axis
                          xt=x*ctz+y*stz
                          y=y*ctz-x*stz
                          x=xt
                 endif
                 if(xang(indre).ne.0.)then
c                         in this case:rotation/x axis
                          yt=y*ctx+z*stx
                          z=z*ctx-y*stx
                          y=yt
                 endif

                 if(yang(indre).ne.0.)then
c                         in this case:rotation/y axis
                          zt=z*cty+x*sty
                          x=x*cty-z*sty
                          z=zt
                 endif
                 xr=x
                 yr=y
                 zr=z
                 x=xa
                 y=ya
                 z=za
                 if(zang(indre).ne.0.)then
c                         in this case:rotation/z axis
                          xt=x*ctz+y*stz
                          y=y*ctz-x*stz
                          x=xt
                 endif 
                if(xang(indre).ne.0.)then
c                         in this case:rotation/x axis
                          yt=y*ctx+z*stx
                          z=z*ctx-y*stx
                          y=yt
                 endif
                 if(yang(indre).ne.0.)then
c                         in this case:rotation/y axis
                          zt=z*cty+x*sty
                          x=x*cty-z*sty
                          z=zt
                 endif
                 xb=x
                 yb=y
                 zb=z
c      save /crelfr/,/relfra/
      return
      end
c
      subroutine arot(xu,yu,zu,xv,yv,zv,ltr)
c*************************************************************************
c                                                                        *
c                       a r o t                                          *
c                                                                        *
c*************************************************************************
c
c             move the plot data point (xu,yu,zu)
c            by translation (no translation if ltr=.false.)
c                and rotation / z , x and y axis.
c
c       input:
c             xu,yu,zu
c             xtp,ytp,ztp=translation vector
c             zap,xap,yap=set of angles of rotation in radian
c       output:
c              xv,yv,zv
c
      include 'snake.inc'
      logical ltr
c                next line : begening of the /plot/ package
      logical lbox,luser,bare
      character unit*3,name*2
      common /cplot/unit(15),name(15)
      common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @              ,nplvold,h1(2,maxq),v1(2,maxq),nq
     @              ,ip,ih,iv,kp,bplot(3,maxplv)
     @              ,plobox(5,4,2,maxre),nface(maxre),lbox
     @              ,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @    ,iplpi(maxplv),iplpf(maxplv)
     @    ,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @		     ,luser,nclic,lclic(maxre),clight(4,maxre),
     @    epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3),nrs(maxre)
c                previous line : end of the /plot/ package
        if(ltr)then
                 x=xu+xtp
                 y=yu+ytp
                 z=zu+ztp
        else
                 x=xu
                 y=yu
                 z=zu
        endif
        if(zap.ne.0.)then
c                in this case:rotation/z axis
                 stz=sin(zap)
                 ctz=cos(zap)
                 xt=x*ctz-y*stz
                 y=y*ctz+x*stz
                 x=xt
        endif
        if(xap.ne.0.)then
c                in this case:rotation/x axis
                 stx=sin(xap)
                 ctx=cos(xap)
                 yt=y*ctx-z*stx
                 z=z*ctx+y*stx
                 y=yt
        endif
        if(yap.ne.0.)then
c                in this case:rotation/y axis
                 sty=sin(yap)
                 cty=cos(yap)
                 zt=z*cty-x*sty
                 x=x*cty+z*sty
                 z=zt
        endif
        xv=x
        yv=y
        zv=z
c      save /cplot/,/plot/
      return
      end
c
c
      subroutine add_field
c*************************************************************************
c                                                                        *
c                       a d d _ f i e l d                                *
c                                                                        *
c*************************************************************************
c
      include 'snake.inc'
      common /absfra/ vea(maxve,0:maxep,21),nuve,nuep,nure
c                next line : begening of the /relfra/ package
      character epsw*4,rname*8,rnamer*8,title*80,cerfiln*20
      common /crelfr/epsw(maxre),rname(maxre),rnamer(maxre)
     @                 ,title(maxre),cerfiln(maxepr,maxre)
      common /relfra/ ver(maxve,12),xyz0(maxre,3)
     @                 ,frbox(2,maxre,3)
     @                 ,zang(maxre),xang(maxre),yang(maxre)
     @                 ,indve,indre,indep,ilive,adata(maxdat,maxre)
     @                 ,indrer(maxre),tept(maxepr,maxre)
     @                 ,yepmin(maxre),yepmax(maxre),kept(maxepr,maxre)
     @                 ,yepstp(maxre),yept(maxepr,maxre),nep(maxre)
     @                 ,colli(maxepr,maxre,6)
c                previous line : end of the /relfra/ package
c                next line : begening of the /mapp/ package
      character*50 form
      character*1 sym
      common /chfield/sym(3),form
      common /mapp/fc(4),b(maxx,maxy,maxz,3),db(maxy),ylis(maxy)
     @  ,ndx,xmin,xmax,xstep,ndy,ymin,ymax,ystep,ndz,zmin,zmax,zstep
     @  ,jmem
c                previous line : end of the /mapp/ package
c                next line : begening of the /fidef/ package
      logical cyl,lfil
      character*80 fifi
      common /fidef/method(maxre),itype(maxre),indfi(maxre),fact(maxre)
     @			,cyl(maxre),lfil
      common /cfidef/ fifi(maxre)
c                previous line : end of the /fidef/ package
      dimension indre_arg(10)
      character nomreg*8,sym_cour*1
      dimension nomreg(10),sym_cour(3),xyz(3),f(3)
      dimension b1(maxx,maxy,maxz,3),ylis_cour(maxy)
      logical lmap,ladd,eject
c
c
c  extract the name and index of the argument regions:
      i1=2
      i2=i1
      call nextcoma(i2,fifi(indre))
      read(fifi(indre)(i1:i2-1),*)ibarg
      if(ibarg.gt.10)stop'add_field: too many regions'
      do ib=1,ibarg
		i1=i2+1
		i2=i1
        call nextcoma(i2,fifi(indre))
                read(fifi(indre)(i1:i2-1),'(a)')nomreg(ib)
            do 5 i=1,nure
                 	if(nomreg(ib).eq.rname(i))then
                          isave=i
                          goto6
                 	endif
 5            continue
       write(6,*)' can''t find this region :',nomreg(ib)
       stop
 6      indre_arg(ib)=isave
      if(fifi(indre_arg(ib))(1:1).eq.'#')then
         write(6,*)'add_field#1: region ',nomreg(ib),
     @ 			'is already an addition type region!'
         stop
        endif
      if(method(indre_arg(ib)).ne.3)then
         write(6,*)'add_field: region ',nomreg(ib),
     @ 			'is not a field distribution type region!'
         stop
        endif
      enddo
      i1=i2+1
      i2=i1
      call nextcoma(i2,fifi(indre))
c extract the symmmetry word and save the current region data:
      read(fifi(indre)(i1:i2-1),'(3a1)')sym_cour
      do l=1,3
        if(sym_cour(l).ne.'s'.and.sym_cour(l).ne.'a'.and.
     @		sym_cour(l).ne.'n'.and.
     @		sym_cour(l).ne.'S'.and.sym_cour(l).ne.'A'.and.
     @		sym_cour(l).ne.'N')
     @		stop'add_field: invalid symmetry word'
      enddo
      if(cyl(indre))
     @		stop'add_field: cylindrical map not yet implemented'
      indre_cour=indre
      ndx_cour=1
      xmin_cour=0.
      xmax_cour=0.
      xstep_cour=0.
      ndy_cour=1
      ymin_cour=0.
      ymax_cour=0.
      ystep_cour=0.
      ndz_cour=1
      zmin_cour=0.
      zmax_cour=0.
      zstep_cour=0.
      goto(801,802,803,804,805,806)indfi(indre)
c
c		1 dim. Poisson type map
c
 801      ndy_cour=nint(adata(1,indre))
      ymin_cour=adata(2,indre)
      ymax_cour=adata(3,indre)
      ystep_cour=(ymax_cour-ymin_cour)/real(ndy_cour-1)
      goto 100
c
c		2 dim. measured type map, uniform y steps
c
 802      ndx_cour=nint(adata(1,indre))
      xmin_cour=adata(2,indre)
      xmax_cour=adata(3,indre)
      ndy_cour=nint(adata(4,indre))
      ymin_cour=adata(5,indre)
      ymax_cour=adata(6,indre)
      xstep_cour=(xmax_cour-xmin_cour)/real(ndx_cour-1)
      ystep_cour=(ymax_cour-ymin_cour)/real(ndy_cour-1)
      goto 100
c
c		2 dim. measured type map, variable y steps
c
 803      ndx_cour=nint(adata(1,indre))
      xmin_cour=adata(2,indre)
      xmax_cour=adata(3,indre)
      ndy_cour=nint(adata(4,indre))
      do j=1,ndy_cour
         ylis_cour(j)=adata(j+4,indre)
      enddo
      xstep_cour=(xmax_cour-xmin_cour)/real(ndx_cour-1)
      goto 100
c
c		2 dim. Poisson type map
c
 804      ndy_cour=nint(adata(1,indre))
      ymin_cour=adata(2,indre)
      ymax_cour=adata(3,indre)
      ndz_cour=nint(adata(4,indre))
      zmin_cour=adata(5,indre)
      zmax_cour=adata(6,indre)
      ystep_cour=(ymax_cour-ymin_cour)/real(ndy_cour-1)
      zstep_cour=(zmax_cour-zmin_cour)/real(ndz_cour-1)
      goto 100
c
c		3 dim. Poisson type map
c
 805      continue
c index=6: same as 5
 806      ndx_cour=nint(adata(1,indre))
      xmin_cour=adata(2,indre)
      xmax_cour=adata(3,indre)
      ndy_cour=nint(adata(4,indre))
      ymin_cour=adata(5,indre)
      ymax_cour=adata(6,indre)
      ndz_cour=nint(adata(7,indre))
      zmin_cour=adata(8,indre)
      zmax_cour=adata(9,indre)
      xstep_cour=(xmax_cour-xmin_cour)/real(ndx_cour-1)
      ystep_cour=(ymax_cour-ymin_cour)/real(ndy_cour-1)
      zstep_cour=(zmax_cour-zmin_cour)/real(ndz_cour-1)
      goto 100
c
 100       if(ndx_cour.gt.maxx)
     @		stop'ifield:max. field points in x surpassed'
        if(ndy_cour.gt.maxy)
     @		stop'ifield:max. field points in y surpassed'
        if(ndz_cour.gt.maxz)
     @		stop'ifield:max. field points in z surpassed'
c
c		built the field map of the current region
c	by adding the field distribution of the two argument regions
c
      do i=1,ndx_cour
            do j=1,ndy_cour
                  do k=1,ndz_cour
                    do l=1,3
					b1(i,j,k,l)=0.
                    enddo
                  enddo
            enddo
      enddo
c in case of diag. in field:
      indve=1
c
      do ib=1,ibarg
	   indre=indre_arg(ib)
	   lmap=.true.
	   call ifield(lmap,ladd)
       if(ladd)then
         write(6,*)'add_field#2: region ',nomreg(ib),
     @ 			'is already an addition type region!'
            stop
         endif
      do i=1,ndx_cour
	      x_cour=xmin_cour+xstep_cour*real(i-1)
	      do j=1,ndy_cour
		 y_cour=ymin_cour+ystep_cour*real(j-1)
	         do k=1,ndz_cour
		    z_cour=zmin_cour+zstep_cour*real(k-1)
		    indre=indre_cour
		    call plorot(x_cour,y_cour,z_cour,
     @			        x_abs,y_abs,z_abs,
     @			        xb,yb,zb)
		    indre=indre_arg(ib)
		    call unplorot(x_abs,y_abs,z_abs,
     @			        xyz(1),xyz(2),xyz(3),
     @			        xb,yb,zb)
c don't eject if outside the free box:
		    call field(eject,xyz,f,
     @				.false.,.false.,.false.,.false.)
c no contribution to the final field in case of eject:
		    if(.not.eject)then
		        call plorot(f(1),f(2),f(3),
     @			            xb,yb,zb,			
     @			            bx_abs,by_abs,bz_abs)
		        indre=indre_cour
		        call unplorot(bx_abs,by_abs,bz_abs,
     @			            xb,yb,zb,
     @			            bx_cour,by_cour,bz_cour)
			b1(i,j,k,1)=b1(i,j,k,1)+bx_cour
			b1(i,j,k,2)=b1(i,j,k,2)+by_cour
			b1(i,j,k,3)=b1(i,j,k,3)+bz_cour
		    endif
         enddo
       enddo
      enddo
      enddo
c
c  initialize the current region field map:
      indre=indre_cour
      do l=1,3
		sym(l)=sym_cour(l)
      enddo
      ndx=ndx_cour
      xmin=xmin_cour
      xmax=xmax_cour
      xstep=xstep_cour
      ndy=ndy_cour
      ymin=ymin_cour
      ymax=ymax_cour
      ystep=ystep_cour
      ndz=ndz_cour
      zmin=zmin_cour
      zmax=zmax_cour
      zstep=zstep_cour
      jmem=2
      do i=1,ndx
        do j=1,ndy
		        if(indfi(indre).eq.3)ylis(j)=ylis_cour(j)
            do k=1,ndz
                do l=1,3
					b(i,j,k,l)=b1(i,j,k,l)
                enddo
            enddo
        enddo
      enddo
c
c		1 dim. Poisson type map: compute the gradient
c
      if(indfi(indre).eq.1)then
	     if(ndy.lt.2)stop'add_field: not enough y in a index=1 map'
	     db(1)=(b(1,2,1,3)-b(1,1,1,3))/ystep
	     do j=2,ndy-1
		db(j)=(b(1,j+1,1,3)-b(1,j-1,1,3))/(2*ystep)
	     enddo
	     db(ndy)=(b(1,ndy,1,3)-b(1,ndy-1,1,3))/ystep
      endif
c      save /absfra/,/crelfr/,/relfra/,/chfield/,/mapp/,/fidef/,/cfidef/
      return
      end
c
c
      subroutine ifield(lmap,ladd)
c*************************************************************************
c                                                                        *
c                       i f i e l d                                      *
c                                                                        *
c*************************************************************************
c
c  initialize the field in the new region for the subroutine send and/or field
c        (input and output list depend on methode,type and field index)
c
c       parameters:
c             maxx:maximum number of x values for the field
c             maxy:   "      "    "  y  "      "   "   "
c             maxz:   "      "    "  z  "      "   "   "
c             maxdat: "      "    "  additional data to describe a region
c       input:
c             indre=current region   index
c             method(region index):methode to use to carry particules
c                in this region
c             itype(region index):type in this method
c             indfi(region index):index in this type and in this method
c             fifi(region index):field file name
c             adata(data index < maxdat,region index):additional data to
c                describe a region
c       output for send:
c             a(6,6):1st order transport type matrix
c       output for subroutine field:
c             indfi(region index): over-read on the field-file
c             sym(x y or z)='s' for symmetry,='a' for antisymmetry,else='n'
c             ndx,xmin,xmax:control param. for x-mesh
c             ndy,ymin,ymax:       "           y  "
c             ndz,zmin,zmax:       "           z  "
c             form:read format for field and/or y values
c             b(xindex,yindex,zindex,bx by or bz):components of the magn.field
c             db(yindex)=y component of the gradient of the field
c             ylis(yindex)=y value for yindex in case of irregular mesh
c		In case of field map to be read on a disk file:
c   lmap=.t. : read the header, the field map and the user plot (if any)
c   lmap=.f. : read only the header because the map and the user plot will be 
c             overwriten by fil.
c
      include 'snake.inc'
c                next line : begening of the /plot/ package
      logical lbox,luser,bare
      character unit*3,name*2
      common /cplot/unit(15),name(15)
      common /plot/plott(maxplp,maxplv,9),nplv,nplp(maxplv)
     @              ,nplvold,h1(2,maxq),v1(2,maxq),nq
     @              ,ip,ih,iv,kp,bplot(3,maxplv)
     @              ,plobox(5,4,2,maxre),nface(maxre),lbox
     @              ,boxp(8,maxre,3),xtp,ytp,ztp,xap,yap,zap,bare
     @		    ,iplpi(maxplv),iplpf(maxplv)
     @		    ,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @		     ,luser,nclic,lclic(maxre),clight(4,maxre),
     @    epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3),nrs(maxre)
c                previous line : end of the /plot/ package
c                next line : begening of the /relfra/ package
      character epsw*4,rname*8,rnamer*8,title*80,cerfiln*20
      common /crelfr/epsw(maxre),rname(maxre),rnamer(maxre)
     @                 ,title(maxre),cerfiln(maxepr,maxre)
      common /relfra/ ver(maxve,12),xyz0(maxre,3)
     @                 ,frbox(2,maxre,3)
     @                 ,zang(maxre),xang(maxre),yang(maxre)
     @                 ,indve,indre,indep,ilive,adata(maxdat,maxre)
     @                 ,indrer(maxre),tept(maxepr,maxre)
     @                 ,yepmin(maxre),yepmax(maxre),kept(maxepr,maxre)
     @                 ,yepstp(maxre),yept(maxepr,maxre),nep(maxre)
     @                 ,colli(maxepr,maxre,6)
c                previous line : end of the /relfra/ package
c                next line : begening of the /mapp/ package
      character*50 form
      character*1 sym
      common /chfield/sym(3),form
      common /mapp/fc(4),b(maxx,maxy,maxz,3),db(maxy),ylis(maxy)
     @  ,ndx,xmin,xmax,xstep,ndy,ymin,ymax,ystep,ndz,zmin,zmax,zstep
     @  ,jmem
c                previous line : end of the /mapp/ package
c   for meth=3 , type=2 , index=1 (dipole with grad. ):
      common/dipgra/rm,alp,bet
c   for meth=3 , type=2 , index=16+17 (arbitrary field numeric.deriv.)
      common/dipff/xtra,ytra,atra,dref,rms,alps,bets
     @           ,xc(2),yc(2),r0(2),e0(2),s0(2),s1(2),s2(2),s3(2)
     @           ,tbound,dst
c   for meth=3 , type=2 , index=2 to 14 and 18 (raytrace library ):
      common/raytra/fact1(5),qrad1,fact2(5),qrad2,dl
      common/diprt/dipdata(75)
c   for meth=3 , type=2 , index=15 (clam):
      common/clam1/d0,dx,dy,bdref,ta
c                next line : begening of the /fidef/ package
	logical cyl,lfil
      character*80 fifi
      common /fidef/method(maxre),itype(maxre),indfi(maxre),fact(maxre)
     @			,cyl(maxre),lfil
      common /cfidef/ fifi(maxre)
c                previous line : end of the /fidef/ package
      dimension ad(3),bd(3),cd(3),plotusc(3)
      logical l,lmap,ladd
      common/matrix/a(6,6),dty,s(3,3)
c  2 next=temporary:
c	dimension xyz(3),f(3)
c	logical eject
c
      ladd=.false.
c	write(6,*)"in ifield"
      go to(101,102,103),method(indre)
c
c*********************no field**************************
c
c       nothing to do:
c 101   save /cplot/,/plot/,/crelfr/,/relfra/,/chfield/,/mapp/,/dipgra/,
c     &       /dipff/,/raytra/,/diprt/,/clam1/,/fidef/,/cfidef/,
c     &       /matrix/	
 101   return
c
c**********************matrix:**************************
c
 102     goto(301,302),itype(indre)
c
c**********************1st order matrix*
c
c       fill the current matrix array acordind to the index value:
 301     goto(201,202),indfi(indre)
c
c       matrix for 1000.mm long drif,parallele faces:
 201     dty=1000.
        goto 5
c
c       matrix for 500.mm long drif,parallele faces:
 202     dty=500.
        goto 5
c
 5       continue
c       matrix for drift:
        do 1 i=1,6
                 do 1 j=1,6
                          a(i,j)=0.
1                         s(i,j)=0.
c       x:
        a(1,1)=1.
        s(1,1)=1.
        a(2,1)=dty
c       xp:
        a(2,2)=1.
        s(2,2)=1.
c       z:
        a(3,3)=1.
        s(3,3)=1.
        a(4,3)=dty
c       zp:
        a(4,4)=1.
c       l:
        a(5,5)=1.
c       d:
        a(6,6)=1.
c      save /cplot/,/plot/,/crelfr/,/relfra/,/chfield/,/mapp/,/dipgra/,
c     &       /dipff/,/raytra/,/diprt/,/clam1/,/fidef/,/cfidef/,
c     &       /matrix/
      return
c
c******************2nd order matrix*
c
 302     continue
        stop'ifield:not yet implemented'
c
c************************rung-khuta************************
c
 103     continue
c       initialize the subroutine "field"
        goto(501,502,503),itype(indre)
c
c       itype=1:homogenous field
c
 501     goto(401,402,403),indfi(indre)
c       indfi=1:pure dipolar field given in data:
 401     fc(1)=adata(1,indre)
        fc(2)=adata(2,indre)
        fc(3)=adata(3,indre)
c      save /cplot/,/plot/,/crelfr/,/relfra/,/chfield/,/mapp/,/dipgra/,
c     &       /dipff/,/raytra/,/diprt/,/clam1/,/fidef/,/cfidef/,
c     &       /matrix/
        return
c       indfi=2:pure quadrupolar field given in data:
c  ( fc(1)>0. for the field of a qpole having x as axis,focusing positive
c  particles in y )
c  ( fc(2)>0. for the field of a qpole having y as axis,focusing positive
c  particles in z )
c  ( fc(3)>0. for the field of a qpole having x as axis,focusing positive
c  particles in x )
 402     fc(1)=adata(1,indre)
        fc(2)=adata(2,indre)
        fc(3)=adata(3,indre)
c      save /cplot/,/plot/,/crelfr/,/relfra/,/chfield/,/mapp/,/dipgra/,
c     &       /dipff/,/raytra/,/diprt/,/clam1/,/fidef/,/cfidef/,
c     &       /matrix/
        return
c       indfi=3: dipolar field along z, qpole+sextupole having y as axis:
c  ( fc(1)=dipolar field, fc(2)>0.=qpolar gradiant focusing in z,
c    fc(3)>0.=sextupolar grad. defocusing in x, no effect in x,
c    fc(4)>0.=sextupolar grad. focusing in z, no effect in x )
 403     fc(1)=adata(1,indre)
        fc(2)=adata(2,indre)
        fc(3)=adata(3,indre)
        fc(4)=adata(4,indre)
c      save /cplot/,/plot/,/crelfr/,/relfra/,/chfield/,/mapp/,/dipgra/,
c     &       /dipff/,/raytra/,/diprt/,/clam1/,/fidef/,/cfidef/,
c     &       /matrix/
        return
c
c       itype=2:analytical field
c
 502     goto(901,902,903,904,905,906,907,908,
     s       909,910,911,912,913,914,915,916,917,918,919,920),
     s       indfi(indre)
c
c         indfi=1;field with gradient:
c
 901     rm=adata(1,indre)
        alp=adata(2,indre)
        bet=adata(3,indre)
c      save /cplot/,/plot/,/crelfr/,/relfra/,/chfield/,/mapp/,/dipgra/,
c     &       /dipff/,/raytra/,/diprt/,/clam1/,/fidef/,/cfidef/,
c     &       /matrix/
        return
c
c         indfi=2:unit gradient (1 t/mm ) along x and z (qpole axis = y )
c                nothing to do
c      save /cplot/,/plot/,/crelfr/,/relfra/,/chfield/,/mapp/,/dipgra/,
c     &       /dipff/,/raytra/,/diprt/,/clam1/,/fidef/,/cfidef/,
c     &       /matrix/
 902  return
c
c         indfi=3:quadrupole
c
 903     fact1(1)=adata(1,indre)
        qrad1=adata(2,indre)
c      save /cplot/,/plot/,/crelfr/,/relfra/,/chfield/,/mapp/,/dipgra/,
c     &       /dipff/,/raytra/,/diprt/,/clam1/,/fidef/,/cfidef/,
c     &       /matrix/
        return
c
c         indfi=4:dipole entrance fringing field:
c
 904     fact1(1)=adata(1,indre)
        qrad1=adata(2,indre)
c      save /cplot/,/plot/,/crelfr/,/relfra/,/chfield/,/mapp/,/dipgra/,
c     &       /dipff/,/raytra/,/diprt/,/clam1/,/fidef/,/cfidef/,
c     &       /matrix/
        return
c
c         indfi=5:dipole exit fringing field:
c
 905     fact1(1)=adata(1,indre)
        qrad1=adata(2,indre)
c      save /cplot/,/plot/,/crelfr/,/relfra/,/chfield/,/mapp/,/dipgra/,
c     &       /dipff/,/raytra/,/diprt/,/clam1/,/fidef/,/cfidef/,
c     &       /matrix/
       return
c
c         indfi=6:uniform dipole :
c
 906     fact1(1)=adata(1,indre)
        fact1(2)=adata(2,indre)
        fact1(3)=adata(3,indre)

c      save /cplot/,/plot/,/crelfr/,/relfra/,/chfield/,/mapp/,/dipgra/,
c     &       /dipff/,/raytra/,/diprt/,/clam1/,/fidef/,/cfidef/,
c     &       /matrix/
        return
c
c         indfi=7:multipole entrance fringing field:
c
 907     fact1(1)=adata(1,indre)
        fact1(2)=-adata(2,indre)
        fact1(3)=adata(3,indre)
        fact1(4)=-adata(4,indre)
        fact1(5)=adata(5,indre)
        qrad1=adata(6,indre)
c			write(6,*)'multipole entrance fringing field '
c      save /cplot/,/plot/,/crelfr/,/relfra/,/chfield/,/mapp/,/dipgra/,
c     &       /dipff/,/raytra/,/diprt/,/clam1/,/fidef/,/cfidef/,
c     &       /matrix/
        return
c
c         indfi=8:uniform multipole :
c
 908     fact1(1)=adata(1,indre)
        fact1(2)=-adata(2,indre)
        fact1(3)=adata(3,indre)
        fact1(4)=-adata(4,indre)
        fact1(5)=adata(5,indre)
        qrad1=adata(6,indre)
c      save /cplot/,/plot/,/crelfr/,/relfra/,/chfield/,/mapp/,/dipgra/,
c     &       /dipff/,/raytra/,/diprt/,/clam1/,/fidef/,/cfidef/,
c     &       /matrix/
        return
c
c         indfi=9:multipole exit fringing field:
c
 909     fact1(1)=adata(1,indre)
        fact1(2)=-adata(2,indre)
        fact1(3)=adata(3,indre)
        fact1(4)=-adata(4,indre)
        fact1(5)=adata(5,indre)
        qrad1=adata(6,indre)
c      save /cplot/,/plot/,/crelfr/,/relfra/,/chfield/,/mapp/,/dipgra/,
c     &       /dipff/,/raytra/,/diprt/,/clam1/,/fidef/,/cfidef/,
c     &       /matrix/
        return
c
c         indfi=10:qq overlaping fields:
c
 910     fact1(1)=adata(1,indre)
        fact1(2)=-adata(2,indre)
        fact1(3)=adata(3,indre)
        fact1(4)=-adata(4,indre)
        fact1(5)=adata(5,indre)
        qrad1=adata(6,indre)
        fact2(1)=adata(7,indre)
        fact2(2)=-adata(8,indre)
        fact2(3)=adata(9,indre)
        fact2(4)=-adata(10,indre)
        fact2(5)=adata(11,indre)
        qrad2=adata(12,indre)
        dl=adata(13,indre)
c      save /cplot/,/plot/,/crelfr/,/relfra/,/chfield/,/mapp/,/dipgra/,
c     &       /dipff/,/raytra/,/diprt/,/clam1/,/fidef/,/cfidef/,
c     &       /matrix/
        return
c
c         indfi=11:qd overlaping fields:
c
 911     fact1(1)=adata(1,indre)
        fact1(2)=-adata(2,indre)
        fact1(3)=adata(3,indre)
        fact1(4)=-adata(4,indre)
        fact1(5)=adata(5,indre)
        qrad1=adata(6,indre)
        fact2(1)=adata(7,indre)
        fact2(2)=adata(8,indre)
	qrad2=adata(9,indre)
        dl=adata(10,indre)
c      save /cplot/,/plot/,/crelfr/,/relfra/,/chfield/,/mapp/,/dipgra/,
c     &       /dipff/,/raytra/,/diprt/,/clam1/,/fidef/,/cfidef/,
c     &       /matrix/
        return
c
c         indfi=12:dd overlaping fields:
c
 912     fact1(1)=adata(1,indre)
        fact1(2)=adata(2,indre)
	qrad1=adata(3,indre)
        fact2(1)=adata(4,indre)
        fact2(2)=adata(5,indre)
	qrad2=adata(6,indre)
        dl=adata(7,indre)
c      save /cplot/,/plot/,/crelfr/,/relfra/,/chfield/,/mapp/,/dipgra/,
c     &       /dipff/,/raytra/,/diprt/,/clam1/,/fidef/,/cfidef/,
c     &       /matrix/
        return
c
c         indfi=13:dq overlaping fields (old version, PV 4/9/97):
c
 913     fact1(1)=adata(1,indre)
        fact1(2)=adata(2,indre)
	qrad1=adata(3,indre)
        fact2(1)=adata(4,indre)
        fact2(2)=-adata(5,indre)
        fact2(3)=adata(6,indre)
        fact2(4)=-adata(7,indre)
        fact2(5)=adata(8,indre)
        qrad2=adata(9,indre)
        dl=adata(10,indre)
c      save /cplot/,/plot/,/crelfr/,/relfra/,/chfield/,/mapp/,/dipgra/,
c     &       /dipff/,/raytra/,/diprt/,/clam1/,/fidef/,/cfidef/,
c     &       /matrix/
        return
c
c         indfi=14:dq overlaping fields (new version, PV 4/9/97):
c
 914     fact1(1)=adata(1,indre)
        fact1(2)=adata(2,indre)
	qrad1=adata(3,indre)
        fact2(1)=adata(4,indre)
        fact2(2)=-adata(5,indre)
        fact2(3)=adata(6,indre)
        fact2(4)=-adata(7,indre)
        fact2(5)=adata(8,indre)
        qrad2=adata(9,indre)
        dl=adata(10,indre)
        open(3,file=fifi(indre),status='old')
        read(3,107)(dipdata(j),j=1,5),(dipdata(j),j=11,22 ),
     s          ( dipdata(j ) , j=25,64)
        close(3)
c      save /cplot/,/plot/,/crelfr/,/relfra/,/chfield/,/mapp/,/dipgra/,
c     &       /dipff/,/raytra/,/diprt/,/clam1/,/fidef/,/cfidef/,
c     &       /matrix/
        return
c
c         indfi=15:clam central field:
c
 915     bref=adata(1,indre)
        dref=adata(2,indre)
        bdref=bref*dref
        beta=adata(3,indre)*dtor
c           % degre assumed!!!
        alpha=adata(4,indre)
c           % radian assumed!!
        ta=tan(alpha)
        dx=cos(angt)
        dy=sin(angt)
        do 23 i=1,3
                 ad(i)=xyz0(indre,i)
                 bd(i)=0.
 23   continue
        cd(1)=-sin(beta)
        cd(2)= cos(beta)
        cd(3)=0.
        call dist(ad,bd,cd,d02,d0)
c      save /cplot/,/plot/,/crelfr/,/relfra/,/chfield/,/mapp/,/dipgra/,
c     &       /dipff/,/raytra/,/diprt/,/clam1/,/fidef/,/cfidef/,
c     &       /matrix/
        return
c
c         indfi=16;field with gradient , numeric. derivation:
c
 916     rms=adata(1,indre)
        alps=adata(2,indre)
        bets=adata(3,indre)
        xc(1)=adata(4,indre)
        yc(1)=adata(5,indre)
        r0(1)=adata(6,indre)
        e0(1)=adata(7,indre)
        s0(1)=adata(8,indre)
        s1(1)=adata(9,indre)
        s2(1)=adata(10,indre)
        s3(1)=adata(11,indre)
        xc(2)=adata(12,indre)
        yc(2)=adata(13,indre)
        r0(2)=adata(14,indre)
        e0(2)=adata(15,indre)
        s0(2)=adata(16,indre)
        s1(2)=adata(17,indre)
        s2(2)=adata(18,indre)
        s3(2)=adata(19,indre)
        tbound=adata(20,indre)*dtor
        dst=adata(21,indre)
	luser=.true.
c      save /cplot/,/plot/,/crelfr/,/relfra/,/chfield/,/mapp/,/dipgra/,
c     &       /dipff/,/raytra/,/diprt/,/clam1/,/fidef/,/cfidef/,
c     &       /matrix/
        return
c
c         indfi=17;dipole field from raytrace, old fashion:
c
 917   open(3,file=fifi(indre),status='old')
      read (3,107) ( dipdata( j ) , j=1,5 ), ( dipdata( j ), j=11,22 ),
     s          ( dipdata( j ) , j=25,64)
 107   format( 5f10.5/ 5f10.5/3f10.5/4f10.5/ 4f10.5/ 6f10.5/ 6f10.5/
     s        6f10.5/ 4f10.5/ 7f10.5/ 7f10.5                           )
      close(3)
c      save /cplot/,/plot/,/crelfr/,/relfra/,/chfield/,/mapp/,/dipgra/,
c     &       /dipff/,/raytra/,/diprt/,/clam1/,/fidef/,/cfidef/,
c     &       /matrix/
      return
c  
c          indfi=18; dipole field from raytrace, new fashion:
c
 918   open(3,file=fifi(indre),status='old')
      read (3,107) ( dipdata( j ) , j=1,5 ), ( dipdata( j ), j=11,22 ),
     s          ( dipdata( j ) , j=25,64)
      close(3)
c      save /cplot/,/plot/,/crelfr/,/relfra/,/chfield/,/mapp/,/dipgra/,
c     &       /dipff/,/raytra/,/diprt/,/clam1/,/fidef/,/cfidef/,
c     &       /matrix/
      return
c
c         indfi=19:helmotz coil
c
 919     fact1(1)=adata(1,indre)
        qrad1=adata(2,indre)
c      save /cplot/,/plot/,/crelfr/,/relfra/,/chfield/,/mapp/,/dipgra/,
c     &       /dipff/,/raytra/,/diprt/,/clam1/,/fidef/,/cfidef/,
c     &       /matrix/
        return
c
c        indfi=20: Solenoid
c
 920    fact1(1)=adata(1,indre)
        fact1(2)=adata(2,indre)
        fact1(3)=adata(3,indre)
c        write(6,*)"collected soln data",fact1
c      save /cplot/,/plot/,/crelfr/,/relfra/,/chfield/,/mapp/,/dipgra/,
c     &       /dipff/,/raytra/,/diprt/,/clam1/,/fidef/,/cfidef/,
c     &       /matrix/
        return
c
c       itype=3:read the file:fifi(indre)
c
 503	continue
c  built a new map by adding two existing field distributions?
      if(fifi(indre)(1:1).eq.'#')then
		ladd=.true.
c      save /cplot/,/plot/,/crelfr/,/relfra/,/chfield/,/mapp/,/dipgra/,
c     &       /dipff/,/raytra/,/diprt/,/clam1/,/fidef/,/cfidef/,
c     &       /matrix/
      return
      endif
        open(3,file=fifi(indre),status='old')
      write(6,'(2a)')'reading map file: ',fifi(indre)
c   overread the interpolation method type to be used:
        read(3,*)indfi(indre)
c   read the symetries of the field:
        read(3,'(3a1)')sym
c   next reads depend on indfi:
        goto(801,802,803,804,805,806),indfi(indre)
c
c       indfi=1:1 dim. poisson type table(bz(0,y,0) and dbz/dy(0,y,0))
c         assume:bz(x,y,z)=bz(0,y,z) and bz(x,y,-z)=bz(x,y,z)
 801     read(3,*)ndy,ymin,ymax
        if(ndy.gt.maxy)stop'ifield:max. field points in y surpassed'
c 		if cylindrical , convert degree to radian :
	if(cyl(indre))then
		ymin=ymin*dtor
		ymax=ymax*dtor
	endif
        ystep=(ymax-ymin)/real(ndy-1)
        q1=ymin+ystep
       q2=ymax
       f1=frbox(1,indre,2)
       a1=abs(f1)
        f2=frbox(2,indre,2)
       a2=abs(f2)
       l=sym(2).eq.'N'.or.sym(2).eq.'n'
        if((l.and.(f1.lt.q1.or.f2.gt.q2)).or.(.not.l.and.
     @     (a1.lt.q1.or.a2.gt.q2.or.a1.gt.q2.or.a2.lt.q1)))
     @     write(6,*)' warning : field < free-box in y'
c   read the format:
        read(3,'(a)')form
c   read the field and gradian:
      if(lfil)then
	  ndy=ndy+2
	  ymin=ymin-ystep
 	  ymax=ymax+ystep
      endif
 	  ndx=1
 	  xmin=0.
 	  xmax=0.
 	  xstep=0.
 	  ndz=1
 	  zmin=0.
 	  zmax=0.
 	  zstep=0.
        if(.not.lmap)then
		close(3)
c      save /cplot/,/plot/,/crelfr/,/relfra/,/chfield/,/mapp/,/dipgra/,
c     &       /dipff/,/raytra/,/diprt/,/clam1/,/fidef/,/cfidef/,
c     &       /matrix/
      return
      endif
		read(3,form)(b(1,i,1,3),db(i),i=1,ndy)
c 		if cylindrical , convert degree to radian :
		if(cyl(indre))then
		do 811 i=1,ndy
 811			db(i)=db(i)*rtod
		endif
      goto 100
c
c       indfi=2:2 dim. measured type table(bz(x,y,0))constant y increment
c         assume: bz(x,y,-z)=bz(x,y,z)
 802     read(3,*)ndx,xmin,xmax
        if(ndx.gt.maxx)stop'ifield:max. field points in x surpassed'
        xstep=(xmax-xmin)/real(ndx-1)
        q1=xmin+xstep
       q2=xmax
       f1=frbox(1,indre,1)
       a1=abs(f1)
        f2=frbox(2,indre,1)
       a2=abs(f2)
       l=sym(1).eq.'N'.or.sym(1).eq.'n'
        if((l.and.(f1.lt.q1.or.f2.gt.q2)).or.(.not.l.and.
     @     (a1.lt.q1.or.a2.gt.q2.or.a1.gt.q2.or.a2.lt.q1)))
     @     write(6,*)' warning : field < free-box in x'
        read(3,*)ndy,ymin,ymax
        if(ndy.gt.maxy)stop'ifield:max. field points in y surpassed'
c 		if cylindrical , convert degree to radian :
      if(cyl(indre))then
		ymin=ymin*dtor
		ymax=ymax*dtor
      endif
        ystep=(ymax-ymin)/real(ndy-1)
        q1=ymin+ystep
       q2=ymax
       f1=frbox(1,indre,2)
       a1=abs(f1)
        f2=frbox(2,indre,2)
       a2=abs(f2)
       l=sym(2).eq.'N'.or.sym(2).eq.'n'
        if((l.and.(f1.lt.q1.or.f2.gt.q2)).or.(.not.l.and.
     @     (a1.lt.q1.or.a2.gt.q2.or.a1.gt.q2.or.a2.lt.q1)))
     @     write(6,*)' warning : field < free-box in y'
      ndz=1
      zmin=0.
      zmax=0.
      zstep=0.
c   read the format:
        read(3,'(a)')form
c   read the field
        if(.not.lmap)then
		close(3)
c      save /cplot/,/plot/,/crelfr/,/relfra/,/chfield/,/mapp/,/dipgra/,
c     &       /dipff/,/raytra/,/diprt/,/clam1/,/fidef/,/cfidef/,
c     &       /matrix/
        return
        endif
          read(3,form)((b(i,j,1,3),i=1,ndx),j=1,ndy)
      goto 100
c
c       indfi=3:2 dim. measured type table(bz(x,y,0))variable y increment
c         assume: bz(x,y,-z)=bz(x,y,z)
 803     read(3,*)ndx,xmin,xmax
        if(ndx.gt.maxx)stop'ifield:max. field points in x surpassed'
        xstep=(xmax-xmin)/real(ndx-1)
        q1=xmin+xstep
       q2=xmax
       f1=frbox(1,indre,1)
       a1=abs(f1)
        f2=frbox(2,indre,1)
       a2=abs(f2)
       l=sym(1).eq.'N'.or.sym(1).eq.'n'
        if((l.and.(f1.lt.q1.or.f2.gt.q2)).or.(.not.l.and.
     @     (a1.lt.q1.or.a2.gt.q2.or.a1.gt.q2.or.a2.lt.q1)))
     @     write(6,*)' warning : field < free-box in x'
        read(3,*)ndy
        if(ndy.gt.maxy)stop'ifield:max. field points in y surpassed'
      ndz=1
      zmin=0.
      zmax=0.
      zstep=0.
c   read the format:
        read(3,'(a)')form
c   read the field
        read(3,form)(ylis(j),j=1,ndy)
c 		if cylindrical , convert degree to radian :
      if(cyl(indre))then
      do 813 j=1,ndy 
 813		ylis(j)=ylis(j)*dtor
      endif
        if(.not.lmap)then
      close(3)
c      save /cplot/,/plot/,/crelfr/,/relfra/,/chfield/,/mapp/,/dipgra/,
c     &       /dipff/,/raytra/,/diprt/,/clam1/,/fidef/,/cfidef/,
c     &       /matrix/
      return
      endif
        read(3,form)((b(i,j,1,3),i=1,ndx),j=1,ndy)
        q1=ylis(2)
       q2=ylis(ndy)
       f1=frbox(1,indre,2)
       a1=abs(f1)
        f2=frbox(2,indre,2)
       a2=abs(f2)
       l=sym(2).eq.'N'.or.sym(2).eq.'n'
        if((l.and.(f1.lt.q1.or.f2.gt.q2)).or.(.not.l.and.
     @     (a1.lt.q1.or.a2.gt.q2.or.a1.gt.q2.or.a2.lt.q1)))
     @     write(6,*)' warning : field < free-box in y'
cprov deb.:
c        open(3,file='new-'//fifi(indre))
c        write(3,*)indfi(indre)
c        write(3,'(3a1)')sym
c        write(3,*)ndx,xmin,xmax
c        write(3,*)ndy
c        write(3,'(a)')form
c        write(3,form)(ylis(j),j=1,ndy)
c        write(3,form)((b(i,j,1,3),i=1,ndx),j=1,ndy)
c        close(3)
cprov fin:
        jmem=2
      goto 100
c
c       indfi=4:2 dim. poisson type table(by(0,y,z) and bz(0,y,z))
c         assume: bz(x,y,z)=bz(0,y,z)
c
804     read(3,*)ndy,ymin,ymax
        if(ndy.gt.maxy)stop'ifield:max. field points in y surpassed'
c 		if cylindrical , convert degree to radian :
      if(cyl(indre))then
		ymin=ymin*dtor
		ymax=ymax*dtor
      endif
        ystep=(ymax-ymin)/real(ndy-1)
        q1=ymin+ystep
       q2=ymax
       f1=frbox(1,indre,2)
       a1=abs(f1)
        f2=frbox(2,indre,2)
       a2=abs(f2)
       l=sym(2).eq.'N'.or.sym(2).eq.'n'
        if((l.and.(f1.lt.q1.or.f2.gt.q2)).or.(.not.l.and.
     @     (a1.lt.q1.or.a2.gt.q2.or.a1.gt.q2.or.a2.lt.q1)))
     @     write(6,*)' warning : field < free-box in y'
        read(3,*)ndz,zmin,zmax
        if(ndz.gt.maxz)stop'ifield:max. field points in z surpassed'
        zstep=(zmax-zmin)/real(ndz-1)
        q1=zmin+zstep
       q2=zmax
       f1=frbox(1,indre,3)
       a1=abs(f1)
        f2=frbox(2,indre,3)
       a2=abs(f2)
       l=sym(3).eq.'N'.or.sym(3).eq.'n'
        if((l.and.(f1.lt.q1.or.f2.gt.q2)).or.(.not.l.and.
     @     (a1.lt.q1.or.a2.gt.q2.or.a1.gt.q2.or.a2.lt.q1)))
     @     write(6,*)' warning : field < free-box in z'
      ndx=1
      xmin=0.
      xmax=0.
      xstep=0.
c   read the format:
        read(3,'(a)')form
c   read the field
        if(.not.lmap)then
        close(3)
c      save /cplot/,/plot/,/crelfr/,/relfra/,/chfield/,/mapp/,/dipgra/,
c     &       /dipff/,/raytra/,/diprt/,/clam1/,/fidef/,/cfidef/,
c     &       /matrix/
        return
        endif
        read(3,form)((b(1,j,k,2),b(1,j,k,3),j=1,ndy),k=1,ndz)
      goto 100
c
c       indfi=5:3 dim table
c
c
 805	continue
c index 6: same initialization as index 5:
 806     read(3,*)ndx,xmin,xmax
        if(ndx.gt.maxx)stop'ifield:max. field points in x surpassed'
        xstep=(xmax-xmin)/real(ndx-1)
        if(indfi(indre).eq.5)then
        	q1=xmin+xstep
        else
        	q1=xmin
        endif
       q2=xmax
       f1=frbox(1,indre,1)
       a1=abs(f1)
        f2=frbox(2,indre,1)
       a2=abs(f2)
       l=sym(1).eq.'N'.or.sym(1).eq.'n'
        if((l.and.(f1.lt.q1.or.f2.gt.q2)).or.(.not.l.and.
     @     (a1.lt.q1.or.a2.gt.q2.or.a1.gt.q2.or.a2.lt.q1)))
     @     write(6,*)' warning : field < free-box in x'
        read(3,*)ndy,ymin,ymax
        if(ndy.gt.maxy)stop'ifield:max. field points in y surpassed'
c 		if cylindrical , convert degree to radian :
      if(cyl(indre))then
		ymin=ymin*dtor
		ymax=ymax*dtor
      endif
        ystep=(ymax-ymin)/real(ndy-1)
        if(indfi(indre).eq.5)then
        	q1=ymin+ystep
        else
        	q1=ymin
        endif
        q2=ymax
       f1=frbox(1,indre,2)
       a1=abs(f1)
        f2=frbox(2,indre,2)
       a2=abs(f2)
       l=sym(2).eq.'N'.or.sym(2).eq.'n'
        if((l.and.(f1.lt.q1.or.f2.gt.q2)).or.(.not.l.and.
     @     (a1.lt.q1.or.a2.gt.q2.or.a1.gt.q2.or.a2.lt.q1)))
     @     write(6,*)' warning : field < free-box in y'
        read(3,*)ndz,zmin,zmax
        if(ndz.gt.maxz)stop'ifield:max. field points in z surpassed'
        zstep=(zmax-zmin)/real(ndz-1)
        if(indfi(indre).eq.5)then
        	q1=zmin+zstep
        else
        	q1=zmin
        endif
        q2=zmax
       f1=frbox(1,indre,3)
       a1=abs(f1)
        f2=frbox(2,indre,3)
       a2=abs(f2)
       l=sym(3).eq.'N'.or.sym(3).eq.'n'
        if((l.and.(f1.lt.q1.or.f2.gt.q2)).or.(.not.l.and.
     @     (a1.lt.q1.or.a2.gt.q2.or.a1.gt.q2.or.a2.lt.q1)))
     @     write(6,*)' warning : field < free-box in z'
c   read the format:
       		 read(3,'(a)')form
        if(.not.lmap)then
      close(3)
c      save /cplot/,/plot/,/crelfr/,/relfra/,/chfield/,/mapp/,/dipgra/,
c     &       /dipff/,/raytra/,/diprt/,/clam1/,/fidef/,/cfidef/,
c     &       /matrix/
      return
      endif
c   read the field
       		 if(form(1:3).eq.'bin')then
       		     read(3)(((b(i,j,k,1),b(i,j,k,2),b(i,j,k,3)
     @                   ,i=1,ndx),j=1,ndy),k=1,ndz)
       		 else
        	    read(3,form)(((b(i,j,k,1),b(i,j,k,2),b(i,j,k,3)
     @                   ,i=1,ndx),j=1,ndy),k=1,ndz)
       		 endif
	goto100
c
c	read the user plot asuming that it is in relat. coord.,
c   convert it in absol. coord. and add it to previous data in /plot/:
 100      read(3,*,end=200)inulu
      do 300 i=1,inulu
		nulu=nulu+1
        if(nulu.gt.nulumax)then
	  	  nulu=nulumax
		  write(6,*)' warning ifield: too many user line'
          goto 200
        endif
		luser=.true.
		read(3,*)nuseuc
		nuseu(nulu)=nuseuc
		if(nuseu(nulu).gt.nuseumax)nuseu(nulu)=nuseumax
      do 400 is=1,nuseuc
         read(3,'(3e12.5)')plotusc
         if(is.le.nuseumax)then
         call plorot(plotusc(1),plotusc(2),plotusc(3),
     @		plotus(is,nulu,1),plotus(is,nulu,2),plotus(is,nulu,3),
     @		xb,yb,zb)
         endif
 400      continue
 300      continue
 200      close(3)
      write(6,'(a)')'done'
c      save /cplot/,/plot/,/crelfr/,/relfra/,/chfield/,/mapp/,/dipgra/,
c     &       /dipff/,/raytra/,/diprt/,/clam1/,/fidef/,/cfidef/,
c     &       /matrix/
      return
      end
      subroutine field(eject,xyz,f,logx,logy,logz,lfreeb)
c*************************************************************************
c                                                                        *
c                       f i e l d                                        *
c                                                                        *
c*************************************************************************
c
c                find the 3 components of the field f(3) at the point xyz(3)
c         field is called by rungk  .true. output value for eject causes
c         the death of the particle.
c
c       input:
c             xyz(x y or z):relative coord. of the point
c       input for subroutine ifield:
c             fact(region index):multipl.factor of the field or ref.mommentum
c             indfi(region index): over-read on the field-file(if any)
c             sym(x y or z)='s' for symmetry,='a' for antisymmetry,else='n'
c             ndx,xmin,xmax:control param. for x-mesh (if any)
c             ndy,ymin,ymax:       "           y  "      "
c             ndz,zmin,zmax:       "           z  "      "
c             form:read format for field (b),gradient (db) and y values (ylis)
c             b(xindex,yindex,zindex,bx by or bz):components of the magn.field
c             db(yindex)=y component of the gradient of the field
c             ylis(yindex)=y value for yindex in case of irregular mesh
c       output:
c             f(x y or z):field component along x y or z
c             eject:value .true. cause the death of the particule
c
      include 'snake.inc'
c                next line : begening of the /diag/ package
      character diag*100,diagt*100
      common/cdiag/diag(maxve),diagt
      common/ndiag/iepdead(maxve),iredead(maxve)
c                previous line : end of the /diag/ package
      character mdiag*3,qdiag(6)*1
c                next line : begening of the /relfra/ package
      character epsw*4,rname*8,rnamer*8,title*80,cerfiln*20
      common /crelfr/epsw(maxre),rname(maxre),rnamer(maxre)
     @                 ,title(maxre),cerfiln(maxepr,maxre)
      common /relfra/ ver(maxve,12),xyz0(maxre,3)
     @                 ,frbox(2,maxre,3)
     @                 ,zang(maxre),xang(maxre),yang(maxre)
     @                 ,indve,indre,indep,ilive,adata(maxdat,maxre)
     @                 ,indrer(maxre),tept(maxepr,maxre)
     @                 ,yepmin(maxre),yepmax(maxre),kept(maxepr,maxre)
     @                 ,yepstp(maxre),yept(maxepr,maxre),nep(maxre)
     @                 ,colli(maxepr,maxre,6)
c                previous line : end of the /relfra/ package
c                next line : begening of the /mapp/ package
      character*50 form
      character*1 sym
      common /chfield/sym(3),form
      common /mapp/fc(4),b(maxx,maxy,maxz,3),db(maxy),ylis(maxy)
     @  ,ndx,xmin,xmax,xstep,ndy,ymin,ymax,ystep,ndz,zmin,zmax,zstep
     @  ,jmem
c                previous line : end of the /mapp/ package
c                next line : begening of the /fidef/ package
      logical cyl,lfil
      character*80 fifi
      common /fidef/method(maxre),itype(maxre),indfi(maxre),fact(maxre)
     @			,cyl(maxre),lfil
      common /cfidef/ fifi(maxre)
c                previous line : end of the /fidef/ package
      dimension f(3),xyz(3),xyzc(3)
      logical eject,logx,logy,logz,ok,lfreeb
      data qdiag/'x','y','z','r','t','z'/
c	write(6,*)'in field'
c
c
c
c		convert to cylindrical coord. ?
c
      if(cyl(indre))then
c		save the cartesian coord. in xyzc :
      do 13 i=1,3
 13		xyzc(i)=xyz(i)
		xyz(1)=sqrt(xyzc(1)*xyzc(1)+xyzc(2)*xyzc(2))
		xyz(2)=atan2(xyzc(2),xyzc(1))
      endif
c
c                in the free-box?
      eject=.false.
c no free box test in case of field addition:
      if(lfreeb)then
        	do 11 k=1,3
         if(xyz(k).lt.frbox(1,indre,k))then
			idiag=k
			mdiag='min'
			eject=.true.
          endif
         if(xyz(k).gt.frbox(2,indre,k))then
			idiag=k
			mdiag='max'
			eject=.true.
         endif
 11      continue
      if(cyl(indre))idiag=idiag+3
c		built the diag:
      if(eject)then
          write(diagt,'(4a)')'exits free box in :'
     s		   ,qdiag(idiag),mdiag,' during step by step ray tracing'
c      save /cdiag/,/ndiag/,/crelfr/,/relfra/,/chfield/,/mapp/,/fidef/,
c     &     /cfidef/	  
      return
      endif
      endif
c
c                get the field acording to its type:
c
c	write(6,*)'in field itype=',itype(indre)
        goto(601,602,603)itype(indre)
c
c                itype=1:homogenous field
c
 601      goto(611,612,613)indfi(indre)
 611     f(1)=fc(1)
        f(2)=fc(2)
        f(3)=fc(3)
        goto 1000
 612     f(1)=fc(2)*xyz(3)+fc(3)*xyz(2)
        f(2)=fc(1)*xyz(3)+fc(3)*xyz(1)
        f(3)=fc(2)*xyz(1)+fc(1)*xyz(2)
        goto 1000
613     f(1)=fc(2)*xyz(3)+fc(3)*xyz(1)*xyz(3)
     s      +fc(4)*(xyz(3)**2-xyz(1)**2)*.5
        f(2)=0.
        f(3)=fc(1)+fc(2)*xyz(1)+fc(3)*(xyz(1)**2-xyz(3)**2)*.5
     s      +fc(4)*xyz(1)*xyz(3)
        goto 1000
c
c                itype=2:analytical field
c
c602     write(6,*)"xyz=",xyz
 602     call analyf(eject,xyz,f,indfi(indre),indve)
c       write(6,*)'return from analyf'
        goto 1000
c
c                itype=3:possible symetry and interpolation
c
 603     continue
c       1/range reduction due to symetry
        x=xyz(1)
        y=xyz(2)
        z=xyz(3)
        if(sym(1).ne.'N'.and.sym(1).ne.'n'.and.x.lt.0.)x=-xyz(1)
        if(sym(2).ne.'N'.and.sym(2).ne.'n'.and.y.lt.0.)y=-xyz(2)
        if(sym(3).ne.'N'.and.sym(3).ne.'n'.and.z.lt.0.)z=-xyz(3)
c       2/interpolation
c
        goto(701,702,703,704,705,706)indfi(indre)
c
c
c       indfi=1:1 dim. poisson type table(bz(0,y,0) and dbz/dy(0,y,0))
c         assume:bz(x,y,z)=bz(0,y,z) and bz(x,y,-z)=bz(x,y,z)
c
701     continue
        stop'field:not yet implemented'
c
c       indfi=2:2 dim. measured type table(bz(x,y,0))constant y increment
c         assume: bz(x,y,-z)=bz(x,y,z)
c
c
702     continue
c        call geo(xyz,eject)
c	if(eject)return
c             finds the mesh
        xr=(x-xmin)/xstep
        ixm=int(xr)
      if(logx)then
		if(mod(ixm,2).eq.0)ixm=ixm-1
		ix=ixm+2
		ixp=ix+2
		p=.5*(xr-real(ix-1))
      else
        	ix=ixm+1
       		ixp=ix+1
	        p=xr-real(ixm)
      endif
        eject=eject.or.ix.lt.0
        eject=eject.or.ix.gt.ndx
        yr=(y-ymin)/ystep
        iym=int(yr)
      if(logy)then
        if(mod(iym,2).eq.0)iym=iym-1
		iy=iym+2
		iyp=iy+2
		q=.5*(yr-real(iy-1))
      else
        	iy=iym+1
       		iyp=iy+1
  	        q=yr-real(iym)
      endif
        eject=eject.or.iy.lt.0
        eject=eject.or.iy.gt.ndy
        if(eject)then
c		built the diag:
        write(diagt,'(5(a,i3))')'exits field map :x index='
     s           ,ix,'[1,',ndx,'] y index=',iy,'[1,',ndy,']'
c      save /cdiag/,/ndiag/,/crelfr/,/relfra/,/chfield/,/mapp/,/fidef/,
c     &     /cfidef/	  
        return
        endif
        f00= b(ix,iy,1,3)
        f0m= b(ix,iym,1,3)
        f0p= b(ix,iyp,1,3)
        fm0= b(ixm,iy,1,3)
        fp0= b(ixp,iy,1,3)
        fpp= b(ixp,iyp,1,3)
      ok=f00.ne.1.e30.and.f0m.ne.1.e30.and.f0p.ne.1.e30
       ok=ok.and.fm0.ne.1.e30.and.fp0.ne.1.e30.and.fpp.ne.1.e30
      if(ok)then
		ixs=1
		iys=1
		goto322
      elseif(logx.or.logy)then
		eject=.true.
c      save /cdiag/,/ndiag/,/crelfr/,/relfra/,/chfield/,/mapp/,/fidef/,
c     &     /cfidef/	  
       return
       endif
c            finds the 6-point formula center
        do 331 jj=1,4
        jx=jj/3
        jy=1-mod(jj,2)
        iwx=1+2*jx
        iwy=1+2*jy
            do 332 iax=ix-jx,ix-jx+iwx,iwx
            do 332 iay=iy-jy,iy-jy+iwy,iwy
            ixcen=iax
            iycen=iay
            if((iax-1).lt.1.or.(iax+1).gt.ndx)goto332
            if((iay-1).lt.1.or.(iay+1).gt.ndy)goto332
                do 303 lx=ixcen-1,ixcen+1
                 do 303 ly=iycen-1,iycen+1
                if(b(lx,ly,1,3).eq.1.e30)goto332
 303             continue
        goto311
 332    continue
 331    continue
 304    eject=.true.
c		built the diag:
      write(diagt,'(a,5(a,i3))')
     s		'can''t find 6 non zero field points'
     s           ,':x index=',ix,'[1,',ndx,'] y index=',iy,'[1,',ndy,']'
c      save /cdiag/,/ndiag/,/crelfr/,/relfra/,/chfield/,/mapp/,/fidef/,
c     &     /cfidef/	  
        return
c             checks the 6th point
 311    ixd=ixcen-ix
        iyd=iycen-iy
        ix0=isign(1,-ixd)
        iy0=isign(1,-iyd)
        do 312 jx=ix0,-ix0,-2*ix0
        do 312 jy=iy0,-iy0,-2*iy0
        ixs=jx
        iys=jy
        if((ixcen+jx).lt.1.or.(ixcen+jx).gt.ndx)goto312
        if((iycen+jy).lt.1.or.(iycen+jy).gt.ndy)goto312
        if(b(ixcen+jx,iycen+jy,1,3).eq.1.e30)goto312
        goto 321
 312    continue
        goto 304
c          2-dim 2nd order inter/extra-polation for bz (x,y,0)
 321    p=(p-real(ixd))*real(ixs)
        q=(q-real(iyd))*real(iys)
        f00= b(ix+ixd     ,iy+iyd     ,1,3)
        f0m= b(ix+ixd     ,iy+iyd-iys ,1,3)
        f0p= b(ix+ixd     ,iy+iyd+iys ,1,3)
        fm0= b(ix+ixd-ixs ,iy+iyd     ,1,3)
        fp0= b(ix+ixd+ixs ,iy+iyd     ,1,3)
        fpp= b(ix+ixd+ixs ,iy+iyd+iys ,1,3)
 322    b30= f0m*q*(q-1.) + fm0*p*(p-1.) + fp0*p*(p-2.*q+1.)
     +        + f0p*q*(q-2.*p+1.)
        b30=.5*b30 + f00*( 1. +p*q - p**2 - q**2) + fpp*p*q
c         1st & 2nd field derivatives  in the (x,y) plane
        ffp1 =  fm0 + fp0 -2.*f00
        ffp2 = -fp0 - f0p + f00 + fpp
        ffq1 =  f0m + f0p -2.*f00
        grx1 = (p*ffp1+q*ffp2+.5*(-fm0+fp0))*real(ixs)/xstep
        gry1 = (p*ffp2+q*ffq1+.5*(-f0m+f0p))*real(iys)/ystep
        grx2 = ffp1/xstep/xstep
        gry2 = ffq1/ystep/ystep
      if(cyl(indre))then
		gry1=gry1/xyz(1)
		xlaplace=grx2+grx1/xyz(1)+gry2/(xyz(1)*xyz(1))
      else
		xlaplace=grx2+gry2
      endif
c                    taylor's  serie for z.ne.0.
        f(1) = z*grx1
        f(2) = z*gry1
        f(3) = b30 -.5*z*z*(xlaplace)
        goto 12
c
c
c       indfi=3:2 dim. measured type table(bz(x,y,0))variable y increment
c         assume: bz(x,y,-z)=bz(x,y,z)
c
 703      continue
      if(cyl(indre))stop' not yet impl. with cylindr. coord. option'
c                find the mesh:
        rim=(x-xmin)/xstep
        im=int(rim)
        ex=rim-real(im)
        eject=eject.or.im.lt.1
        i=im+1
        ip=i+1
        eject=eject.or.i.gt.(ndx-1)
c        use again the previous value of j?:
      j=jmem
c       if(j.lt.2.or.j.ge.ndy)j=2
        jm=j-1
        jp=j+1
        if(y.lt.ylis(j))then
c                search backward:
                  do 3 j=jm,2,-1
                          if(ylis(j).lt.y)then
                                   jp=j+1
                                   goto 4
                          endif
3                continue
                 eject=.true.
                 bid=0.
        endif
        j1=jp+1
        if(y.gt.ylis(jp))then
c                search forward:
                  do 6 jp=j1,ndy
                          j=jp-1
                          if(ylis(jp).gt.y)goto4
6                continue
                 eject=.true.
        endif
4       if(eject)then
                 jmem=2
         write(diagt,'(5(a,i3))')
     s		'exits field map :x index='
     s           ,i,'[1,',ndx,'] y index=',j,'[1,',ndy,']'
c      save /cdiag/,/ndiag/,/crelfr/,/relfra/,/chfield/,/mapp/,/fidef/,
c     &     /cfidef/	  
        return
        endif
        jm=j-1
        yl=ylis(j)
        ylp=ylis(jp)
        ey=(y-yl)/(ylp-yl)
        py=ylp-yl
        ph=yl-ylis(jm)
        px=xstep
        zh=z
c                2dim. interpolation/extrapolation at 2nd order:
      c0=b(i,j,1,3)
      c1=b(im,j,1,3)
      c2=b(ip,j,1,3)
      c3=b(i,jm,1,3)
      c4=b(i,jp,1,3)
      c5=((c4*ph**2-c3*py**2)/(ph+py)-(ph-py)*c0)/ph
      c6=c2+c1-2.*c0
      c7=(c4*ph+c3*py)/(ph+py)-c0
      c8=b(ip,jp,1,3)+c0-c2-c4
      f(1)=(0.5*(c2-c1)+ex*c6+ey*c8)/px*zh
      f(2)=((c5+ex*c8)/py+ey*2.*c7/ph)*zh
      f(3)=c0+0.5*ex*(c2-c1)+ey*c5+0.5*(ex**2-(zh/px)**2)*c6+ex*ey*c8+
     @ (ey**2*py/ph-zh**2/(py*ph))*c7
        goto 12
c
c       indfi=4:2 dim. poisson type table(bz(0,y,z) and by(0,y,z))
c         assume: bz(x,y,z)=bz(0,y,z)
c
704     stop'field:not yet implemented'
c
c       indfi=5:3 dim. table
c
705     continue
c         find the cube:
        xr=((x-xmin)/xstep)+.5
c        xr=(x-xmin)/xstep
        ixm=int(xr)
      if(logx)then
		if(mod(ixm,2).eq.0)ixm=ixm-1
		ix=ixm+2
		ixp=ix+2
		p=.5*(xr-real(ix-1))
      else
        	ix=ixm+1
       		ixp=ix+1
	        p=xr-real(ixm)
      endif
        eject=eject.or.ixm.lt.1
        eject=eject.or.ixp.gt.ndx
c
        yr=((y-ymin)/ystep)+.5
c        yr=(y-ymin)/ystep
        iym=int(yr)
      if(logy)then
        if(mod(iym,2).eq.0)iym=iym-1
		iy=iym+2
		iyp=iy+2
		q=.5*(yr-real(iy-1))
      else
        	iy=iym+1
       		iyp=iy+1
  	        q=yr-real(iym)
      endif   
        eject=eject.or.iym.lt.1
        eject=eject.or.iyp.gt.ndy
c
        zr=((z-zmin)/zstep)+.5
c        zr=(z-zmin)/zstep
        izm=int(zr)
      if(logz)then
        if(mod(izm,2).eq.0)izm=izm-1
		iz=izm+2
		izp=iz+2
		r=.5*(zr-real(iz-1))
      else
        	iz=izm+1
       		izp=iz+1
 	        r=zr-real(izm)
      endif
        eject=eject.or.izm.lt.1
        eject=eject.or.izp.gt.ndz
        if(eject)then
         write(diagt,'(7(a,i3))')
     s		'exits field map :x index='
     s           ,ix,'[1,',ndx,'] y index=',iy,'[1,',ndy
     s           ,'] z index=',iz,'[1,',ndz,']'
c      save /cdiag/,/ndiag/,/crelfr/,/relfra/,/chfield/,/mapp/,/fidef/,
c     &     /cfidef/	  
        return
        endif
c                test the geometry:
c        call geo(xyz,eject)
      if(eject)then
c       save /cdiag/,/ndiag/,/crelfr/,/relfra/,/chfield/,/mapp/,/fidef/,
c     &     /cfidef/	  
       return
      endif
c         3 dim. 2nd order interpolation inside the cube for each
c                component of the field.
c         use 11 values of the field:
        do 7 i=1,3
          f000=b(ix,iy,iz,i)
	  eject=eject.or.(f000.gt.50.)
          f00p=b(ix,iy,izp,i)
	  eject=eject.or.(f00p.gt.50.)
          f00m=b(ix,iy,izm,i)
	  eject=eject.or.(f00m.gt.50.)
          f0p0=b(ix,iyp,iz,i)
	  eject=eject.or.(f0p0.gt.50.)
          f0m0=b(ix,iym,iz,i)
	  eject=eject.or.(f0m0.gt.50.)
          fp00=b(ixp,iy,iz,i)
	  eject=eject.or.(fp00.gt.50.)
          fm00=b(ixm,iy,iz,i)
	  eject=eject.or.(fm00.gt.50.)
          f0pp=b(ix,iyp,izp,i)
	  eject=eject.or.(f0pp.gt.50.)
          fp0p=b(ixp,iy,izp,i)
	  eject=eject.or.(fp0p.gt.50.)
          fpp0=b(ixp,iyp,iz,i)
	  eject=eject.or.(fpp0.gt.50.)
          fppp=b(ixp,iyp,izp,i)
	  eject=eject.or.(fppp.gt.50.)
          if(eject)then
		write(diagt,'(a,7(a,i3))')'try to use undefined '
     s		,'field values :x index='
     s           ,ix,'[1,',ndx,'] y index=',iy,'[1,',ndy
     s           ,'] z index=',iz,'[1,',ndz,']'
c      save /cdiag/,/ndiag/,/crelfr/,/relfra/,/chfield/,/mapp/,/fidef/,
c     &     /cfidef/	  
        return
        endif
          c=2.*f000
          cp=fp00-c+fm00
          cq=f0p0-c+f0m0
          cr=f00p-c+f00m
          dp=f000-fp00+fpp0-f0p0
          dq=f000-f0p0+f0pp-f00p
          dr=f000-f00p+fp0p-fp00
          e=-dp-f0pp+f00p-fp0p+fppp
c                compute the taylor's serie:
          pq=p*q
       pqr=pq*r
       qr=q*r
       pr=p*r
          f(i)=f000+.5*(p*(fp00-fm00)+p*p*cp
     @                 +q*(f0p0-f0m0)+q*q*cq
     @                 +r*(f00p-f00m)+r*r*cr)
     @             +pq*dp+qr*dq+pr*dr+pqr*e
7       continue
        goto 12
c
c
c       indfi=6:3 dim. table, no extra mesh needed for interpolation
c
 706      continue
c         find the cube:
        xr=(x-xmin)/xstep
        xrlim=real(ndx-1)*.5
       if(xr.ge.xrlim)then
      if(logx)then
        	ix=int(xr/2.)*2+1
		ixm=ix-2
		ixp=ix+2
		p=.5*(xr-real(ix-1))
        	eject=eject.or.ixm.lt.1
        	eject=eject.or.ixp.gt.ndx
      	else
        	ixm=int(xr)
        	ix=ixm+1
       		ixp=ix+1
	        p=xr-real(ixm)
        	eject=eject.or.ixm.lt.1
        	eject=eject.or.ixp.gt.ndx
       	endif
       else
      if(logx)then
        	ix=int(xr/2.)*2+3
		ixm=ix+2
		ixp=ix-2
		p=-.5*(xr-real(ix-1))
        	eject=eject.or.ixm.gt.ndx
        	eject=eject.or.ixp.lt.1
      	else
        	ix=int(xr)+2
        	ixp=ix-1
       		ixm=ix+1
	        p=real(ix-1)-xr
        	eject=eject.or.ixm.gt.ndx
        	eject=eject.or.ixp.lt.1
       	endif
       endif
c
        yr=(y-ymin)/ystep
        yrlim=real(ndy-1)*.5
       if(yr.ge.yrlim)then
      if(logy)then
        	iy=int(yr/2.)*2+1
		iym=iy-2
		iyp=iy+2
		q=.5*(yr-real(iy-1))
        	eject=eject.or.iym.lt.1
        	eject=eject.or.iyp.gt.ndy
      	else
        	iym=int(yr)
        	iy=iym+1
       		iyp=iy+1
	        q=yr-real(iym)
        	eject=eject.or.iym.lt.1
        	eject=eject.or.iyp.gt.ndy
       	endif
       else
      if(logy)then
        	iy=int(yr/2.)*2+3
		iym=iy+2
		iyp=iy-2
		q=-.5*(yr-real(iy-1))
        	eject=eject.or.iym.gt.ndy
        	eject=eject.or.iyp.lt.1
      	else
        	iy=int(yr)+2
        	iyp=iy-1
       		iym=iy+1
	        q=real(iy-1)-yr
        	eject=eject.or.iym.gt.ndy
        	eject=eject.or.iyp.lt.1
       	endif
       endif
c
        zr=(z-zmin)/zstep
        zrlim=real(ndz-1)*.5
       if(zr.ge.zrlim)then
      if(logz)then
        	iz=int(zr/2.)*2+1
		izm=iz-2
		izp=iz+2
		r=.5*(zr-real(iz-1))
        	eject=eject.or.izm.lt.1
        	eject=eject.or.izp.gt.ndz
      	else
        	izm=int(zr)
        	iz=izm+1
       		izp=iz+1
	        r=zr-real(izm)
        	eject=eject.or.izm.lt.1
        	eject=eject.or.izp.gt.ndz
       	endif
       else
      if(logz)then
        	iz=int(zr/2.)*2+3
		izm=iz+2
		izp=iz-2
		r=-.5*(zr-real(iz-1))
        	eject=eject.or.izm.gt.ndz
        	eject=eject.or.izp.lt.1
      	else
        	iz=int(zr)+2
        	izp=iz-1
       		izm=iz+1
	        r=real(iz-1)-zr
        	eject=eject.or.izm.gt.ndz
        	eject=eject.or.izp.lt.1
       	endif
       endif 
                if(eject)then
		write(diagt,'(7(a,i3))')
     s		'exits field map :x index='
     s           ,ix,'[1,',ndx,'] y index=',iy,'[1,',ndy
     s           ,'] z index=',iz,'[1,',ndz,']'
c      save /cdiag/,/ndiag/,/crelfr/,/relfra/,/chfield/,/mapp/,/fidef/,
c     &     /cfidef/	  
       return
       endif
c                test the geometry:
c        call geo(xyz,eject)
      if(eject)then
c       save /cdiag/,/ndiag/,/crelfr/,/relfra/,/chfield/,/mapp/,/fidef/,
c     &     /cfidef/	  
       return
      endif
c         3 dim. 2nd order interpolation inside the cube for each
c                component of the field.
c         use 11 values of the field:
        do  i=1,3
          f000=b(ix,iy,iz,i)
	  eject=eject.or.(abs(f000).gt.50.)
          f00p=b(ix,iy,izp,i)
	  eject=eject.or.(abs(f00p).gt.50.)
          f00m=b(ix,iy,izm,i)
	  eject=eject.or.(abs(f00m).gt.50.)
          f0p0=b(ix,iyp,iz,i)
	  eject=eject.or.(abs(f0p0).gt.50.)
          f0m0=b(ix,iym,iz,i)
	  eject=eject.or.(abs(f0m0).gt.50.)
          fp00=b(ixp,iy,iz,i)
	  eject=eject.or.(abs(fp00).gt.50.)
          fm00=b(ixm,iy,iz,i)
	  eject=eject.or.(abs(fm00).gt.50.)
          f0pp=b(ix,iyp,izp,i)
	  eject=eject.or.(abs(f0pp).gt.50.)
          fp0p=b(ixp,iy,izp,i)
	  eject=eject.or.(abs(fp0p).gt.50.)
          fpp0=b(ixp,iyp,iz,i)
	  eject=eject.or.(abs(fpp0).gt.50.)
          fppp=b(ixp,iyp,izp,i)
	  eject=eject.or.(abs(fppp).gt.50.)
          if(eject)then
           write(diagt,'(a,7(a,i3))')'try to use undefined '
     s		,'field values :x index='
     s           ,ix,'[1,',ndx,'] y index=',iy,'[1,',ndy
     s           ,'] z index=',iz,'[1,',ndz,']'
c      save /cdiag/,/ndiag/,/crelfr/,/relfra/,/chfield/,/mapp/,/fidef/,
c     &     /cfidef/	  
          return
          endif
          c=2.*f000
          cp=fp00-c+fm00
          cq=f0p0-c+f0m0
          cr=f00p-c+f00m
          dp=f000-fp00+fpp0-f0p0
          dq=f000-f0p0+f0pp-f00p
          dr=f000-f00p+fp0p-fp00
          e=-dp-f0pp+f00p-fp0p+fppp
c                compute the taylor's serie:
          pq=p*q
       pqr=pq*r
       qr=q*r
       pr=p*r
          f(i)=f000+.5*(p*(fp00-fm00)+p*p*cp
     @                 +q*(f0p0-f0m0)+q*q*cq
     @                 +r*(f00p-f00m)+r*r*cr)
     @             +pq*dp+qr*dq+pr*dr+pqr*e
      enddo
        goto 12
c
c       3/range restitution due to symetry
 12      if(x.ne.xyz(1))then
                 if(sym(1).eq.'S'.or.sym(1).eq.'s')then
                          f(2)=-f(2)
                          f(3)=-f(3)
                 else
                          f(1)=-f(1)
                 endif
        endif
        if(y.ne.xyz(2))then
                 if(sym(2).eq.'S'.or.sym(2).eq.'s')then
                          f(3)=-f(3)
                          f(1)=-f(1)
                 else
                          f(2)=-f(2)
                 endif
        endif
        if(z.ne.xyz(3))then
                 if(sym(3).eq.'S'.or.sym(3).eq.'s')then
                          f(1)=-f(1)
                          f(2)=-f(2)
                 else
                          f(3)=-f(3)
                 endif
        endif
c       4/normalize the field and restore cartesian coord.
c		and field components :
 1000    continue
        if(eject)then
c        save /cdiag/,/ndiag/,/crelfr/,/relfra/,/chfield/,/mapp/,/fidef/,
c     &     /cfidef/	  
        return
        endif
      if(cyl(indre))then
		cst=cos(xyz(2))
		sst=sin(xyz(2))
		bx=f(1)*cst-f(2)*sst
c next line bug corrected 19/06/91: ...-f(2)...
		by=f(1)*sst+f(2)*cst
		f(1)=bx
		f(2)=by
        	do 14 i=1,3
		   xyz(i)=xyzc(i)
14	        continue
      endif
        do 2 i=1,3
                 f(i)=f(i)*fact(indre)
2       continue
c      save /cdiag/,/ndiag/,/crelfr/,/relfra/,/chfield/,/mapp/,/fidef/,
c     &     /cfidef/	  
      return
      end
c
      subroutine dist(a,b,c,d2,d)
c*************************************************************************
c                                                                        *
c                       d i s t                                          *
c                                                                        *
c*************************************************************************
c
c       compute d=distance between a point a and a straight line
c of equation b + d*c (a,b and c = vectors , c assumed of unit module !!)
      dimension a(3),b(3),c(3),bma(3)
        bmasc=0.
      do 29 i=1,3
                 bmac=b(i)-a(i)
                 bma(i)=bmac
                 bmasc=bmasc+bmac*c(i)
 29   continue
        d2=0.
        do 32 i=1,3
                 dc=bma(i)-bmasc*c(i)
                 d2=d2+dc*dc
 32   continue
        d=sqrt(d2)
      return
      end
        
      subroutine mapadd(file,factor,outfile)        
C subroutine that takes 2 3-d maps for SNAKE (output of mapwrt_3d) and adds them 
c together having multiplied each by a scaling factor -JJL 9/24/2004
c
c file(1)=name of 1st map file
c factor(1)=multiplicative factor for 1st map
c file(2)=name of 2nd map file
c factor(2)=multiplicative factor for 2nd map
c outfile=name of output map file
c
      Parameter (nn=900000)
      integer n(2),ndx(2),ndy(2),ndz(2),index
      character*1 typ(2,3)
      real factor(2)
      real xmax(2),ymax(2),zmax(2)
      real xmin(2),ymin(2),zmin(2),b(3,nn)
      character*80 file(2),outfile
      character*10 frmt(2)
      
      write(6,*)file(1),file(2)
      write(6,*)factor(1),factor(2)
      
      
      open(15,file=file(1),status='old',form='formatted')  !,err=97)
      open(16,file=file(2),status='old',form='formatted',err=98)
      
      do i=1,2
       k=i+14
       read(k,'(i2)')n(i)
       read(k,'(3a1)') (typ(i,j),j=1,3)
       read(k,*) ndx(i), xmin(i), xmax(i)
       read(k,*) ndy(i), ymin(i), ymax(i)
       read(k,*) ndz(i), zmin(i), zmax(i)
       read(k,*) frmt(i)
       write(6,'(i2)')n(i)
       write(6,'(3a1)') (typ(i,j),j=1,3)
       write(6,*) ndx(i), xmin(i), xmax(i)
       write(6,*) ndy(i), ymin(i), ymax(i)
       write(6,*) ndz(i), zmin(i), zmax(i)
       write(6,*) frmt(i)
      enddo

c check that files match

      if((n(1).ne.n(2)).or.
     &    (ndx(1).ne.ndx(2)).or.
     &    (ndy(1).ne.ndy(2)).or.
     &    (ndz(1).ne.ndz(2)).or.
     &    (xmax(1).ne.xmax(2)).or.
     &    (ymax(1).ne.ymax(2)).or.
     &    (zmax(1).ne.zmax(2)).or.
     &    (xmin(1).ne.xmin(2)).or.
     &    (ymin(1).ne.ymin(2)).or.
     &    (zmin(1).ne.zmin(2)).or.
     &    (typ(1,1).ne.typ(2,1)).or.
     &    (typ(1,2).ne.typ(2,2)).or.
     &    (typ(1,3).ne.typ(2,3))) then
       go to 96
       else
       write(6,*) 'map files match'
       endif

c open output file

      open(3,file=outfile,status='unknown',form='formatted',err=95)
      write(3,'(i2)')n(1)
      write(3,'(3a1)') (typ(1,j),j=1,3)
      write(3,*) ndx(1), xmin(1), xmax(1)
      write(3,*) ndy(1), ymin(1), ymax(1)
      write(3,*) ndz(1), zmin(1), zmax(1)
      write(3,*) frmt(1) 
      
      index=ndx(1)*ndy(1)*ndz(1)*3
      do i=1,2
       k=i+14
       read(k,(frmt(i)))(b(i,j),j=1,index)
      enddo
      close(15)
      close(16)
      do i=1,index
      b(3,i)=(b(1,i)*factor(1))+(b(2,i)*factor(2))
      enddo
      write(3,(frmt(1)))(b(3,i),i=1,index)
      close(3)    
                
      return 
 95   stop 'error opening outfile'
 96   stop 'map files do not match'
c 97   stop 'open error on file 1'
 98   stop 'open error on file 2'   
      end
      
      
