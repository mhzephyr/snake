
c*************************************************************************
c                                                                        *
c                       s n a k e                                        *
c                                                                        *
c*************************************************************************
c
c        general ray tracing program for multi regions problems
c                       p.vernin 21/07/86
c                version 2 16/07/88
c                version 3 31/05/91
c
c       units:mm,gev/c,t
c       input files:
c             -geometry and field are defined region after region in
c                the "description file"
c             -the description file may itself refer to "field files"
c       output file:
c             -vectors for mudifi (cern program for multi dimensional fit)
c                may be stored on a file
c       parameters:
c             maxve:maximum number of input vector(s) that can be defined
c             maxre:   "      "    "  region(s)       "
c             maxep:   "      "    "  end-plane(s)    "
c             maxepr:  "      "    "  end-plane(s) by region
c             maxplv:  "      "    "  trajectories   that can be plotted
c             maxplp:  "      "    "  point per traj.    "
c             maxq:    "      "    "  tic marks          "
c             maxdat: "      "    "  additional data to describe a region
c       common /absfra/:
c             nuve:number of input vector(s) defined
c             nuep:number of end-plane(s) defined
c             nure:number of region(s)  defined
c             vea(iv,iep,ic):description of the particules at each end-plane
c                iv<=nuve<=maxve:vector index
c                iep<=nuep<=maxep:end-plane index
c                ic=1:x coord. in absolute frame
c                   2:y            "
c                   3:z            "
c                   4:x component of the normalized mommentum in absol.frame
c                   5:y            "
c                   6:z            "
c                   7:p module of the mommentum
c                   8:trace length measured from the "spring point"
c                   9:life-flag(-1.:buried,0.:just dead,1.:alive)
c            10 to 15:the same as 1 to 6 but in relative frame coordinates
c		   16:x component of the spin
c		   17:y          "
c		   18:z          "
c	     19 to 21: the same as 16 to 18 but in relative frame coord.
c       common /relfra/:
c             indre=current region   index
c             indep=   "    end-plane "
c             ilive=   "    number of alive vector(s)
c             ver(iv,ic):description of the particules in the current region
c                iv<=nuve<=maxve:vector index
c                ic=1:x coord.relative to the current region(in relative frame)
c                   2:y            "
c                   3:z            "
c                   4:x component of the normalized mommentum(in relat.frame)
c                   5:y            "
c                   6:z            "
c                   7:p:module of the mommentum
c                   8:trace length measured from the "spring point"
c                   9:life-flag(-1.:buried,0.:just dead,1.:alive)
c		   10:x component of the spin (in relat.frame)
c		   11:y          "
c		   12:z          "
c             xyz0(region index,x y or z):absolute coord. of the rel.origine
c             frbox(maxi or mini,region index,x y or z):relat.coord.
c                of the free-box
c             zang(region index),xang(region index) and yang(region index):
c                rotation angles to define the relat.frame in the absol.frame
c             adata(data index < maxdat,region index):additional data to
c                describe a region
c       common /fidef/:
c             method(region index):methode to use to carry particles
c                in this region
c             itype(region index):type in this method
c             indfi(region index):index in this type and in this method
c             fact(region index):multipl.factor of the field or ref.mommentum
c             fifi(region index):field file name
c       common /plot/:
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
c             kp= 1:observer is at ip>0.
c                -1:   "     "  "  IP<0.
c             plobox(submit index,horiz. or vertic.,region index):abs.coord
c                of submit of free-box(special submit list for 3d projection)
        program snake
c#include </opt/SUNWspro/SC4.2/include/f77/f77_floatingpoint.h>
        include 'snake.inc'
        logical firout
c                next line : begening of the /diag/ package
      character diag*100,diagt*100
      common/cdiag/diag(maxve),diagt
      common/ndiag/iepdead(maxve),iredead(maxve)
c                previous line : end of the /diag/ package
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
      common /comlight/cl(3)
        character*20 fidi,fibl,fidit,fide,fimu,fimut,firay
        character*80 line
        character*1 rep,rip
        character rep2*2
c	character*8 tyme
      character*24 daet 
c                next line : begening of the /fidef/ package
      logical cyl,lfil,lmap,ladd
        character*80 fifi
       common /fidef/method(maxre),itype(maxre),indfi(maxre),fact(maxre)
     @			,cyl(maxre),lfil
       common /cfidef/ fifi(maxre)
c                previous line : end of the /fidef/ package
         common/alter/jflag,prec,hstmax,hstmin,sinc
      dimension xmax(8),xmin(8),plottn(2,maxplv,6)
        data firout/.true./
       data fide/'DIR.DAT'/
c        data fide/'solenoid_dir.dat'/
c       data fide/'comm_8_dir.dat'/
c       data fide/'sn-test3-dir.dat'/
c       data fide/'sn-900-dir.dat'/
c        data fide/'smapl_dir.dat'/
c       data fide/'sn-smc-dir.dat'/
c       data fide/'rtmapb_dir.dat'/
c	data fide/'sn-test1-dir.dat'/
c	data fide/'sn-ber-dir.dat'/
c	data fide/'cone11fb-dir.dat'/
c	data fide/'snk-dir.dat'/
        data fibl/'                    '/
        data rip/'r'/
        data fimu/'mud.mud'/
c        data fimu/'sn.out'/
        data firay/'temp.dat'/
c        data tyme/'00:00:00'/
c        data daet/' unknown '/
c
       COMMON/PAWC/H(20000)
       CALL MZEBRA(-3)
       CALL MZPAW(20000,' ')
       CALL IGINIT(0)

c
      data cl/1.,1.,1./
c	istat=ieee_handler('set','common',SIGFPE_ABORT)
c	if(istat.ne.0)write(6,*)' IEEE_HANDLER "set" action failed'
      lfil=.false.
        name(1)='x'
        name(2)='y'
        name(3)='z'
        name(4)='bx'
        name(5)='by'
        name(6)='bz'
        name(7)='sx'
        name(8)='sy'
        name(9)='sz'
        name(10)='l'
        name(11)='b'
        name(12)='r'
        name(13)='t'
        name(14)='rt'
        name(15)='s'
        unit(1)='mm'
        unit(2)='mm'
        unit(3)='mm'
        unit(4)='T'
        unit(5)='T'
        unit(6)='T'
        unit(7)=' '
        unit(8)=' '
        unit(9)=' '
        unit(10)='mm'
        unit(11)='T'
        unit(12)='mm'
        unit(13)='deg'
        unit(14)='mm'
        unit(15)=' '
c
c
        fidi=fide
      nplv=0
      nplvold=0
100     continue
        write(6,'(/,''$ directive-file name??(def='',a,'')'')')fidi
        read(5,'(a)')fidit
        if(fidit.ne.fibl)fidi=fidit
        write(6,*)fidi
        write(6,*)
        open(2,file=fidi,err=428,status='old')
        goto429
428     continue
        write(6,'('' can''''t open this file'')')
        goto100
c
c                read directives :
c
429     call descin
        close(2)
        nure=indre
c
c                fill the absolute coord. buff.,set region index to 0:
c
528     call spring
c
c                set end-plane index to 0:
c
        indep=0
        indved=0
      nulu=0
      luser=.false.
c**************************new region:*************************
c
        do 200 indrel=1,nure
        indre=indrel
        write(6,*)title(indre)
	if(ilive.le.0)then
		write(6,'(''  alive:'',i4,'' vector, region skipped'')')
     s			ilive
		goto 200
	endif
c
c       fill the relat. coodr. buff. with the abs. one(last e-p index)
c                accordind to new geometrical specif.
        call rot
c                initialize the field routine according to new field specif.
      lmap=.true.
      call ifield(lmap,ladd)
      if(ladd)call add_field
c       loop on the end-planes of this region:
        do 1 iep=1,nep(indre)
                 if(epsw(indre).eq.'loop')then
		     kep=2			
                     yep=yepmin(indre)+yepstp(indre)*(iep-1)
                     tep=0.
                 else
		     kep=kept(iep,indre)
                     yep=yept(iep,indre)
                     tep=tept(iep,indre)
                 endif
c               send the particle to the next end-plane...
c                ...by working on the relat.coord.buffer
c      write(6,*)'calling send'
                 call send(kep,yep,tep,iep)
c                increment the end-plane index and fill the abs.coo.buff.:
                 call unrot
                 write(6,'(''  alive:'',i4,'' vector(s)'')')ilive
1       continue
200     continue
c
c*************************end of the ray-tracing****************
c
c       save the number of end-planes:
        nuep=indep
c		unclip the plot:
      call unclip
        write(6,*)
600     write(6,601)
601     format(/,' p(lot)',
     @          ',d(isplay)',
     @          ',o(utput for mudifi)',
     @          ',u(ser defined output)',/,
     @          ',r(aytrace style output)',
     @          ',n(ew problem),',/,
     @          ' m(ove the plot)',
     @          ',c(lip the plot)',
     @          ',l(ightning)',
     @          ',e(nvelope)',/,
     @          's(ection)',
     @          ',b(are plot)',
     @          ' or q(uit)? (def=p)')
        read(5,'(a)')rep
        if(rep.eq.' ')rep='p'
        if(rep.eq.'P'.or.rep.eq.'p')goto500
        if(rep.eq.'D'.or.rep.eq.'d')goto510
        if(rep.eq.'O'.or.rep.eq.'o')goto520
        if(rep.eq.'U'.or.rep.eq.'u')goto520
        if(rep.eq.'R'.or.rep.eq.'r')goto525
        if(rep.eq.'N'.or.rep.eq.'n')then
    	  write(6,'(''$keep the old plot trajectories (y/n def=n) ?'')')
	        read(5,'(a)')rep
	        if(rep.ne.'Y'.and.rep.ne.'y')then
			nplvold=0
		else
			nplvold=nplv
		endif
		goto400
      endif
        if(rep.eq.'M'.or.rep.eq.'m')goto550
        if(rep.eq.'C'.or.rep.eq.'c')goto560
        if(rep.eq.'L'.or.rep.eq.'l')goto570
        if(rep.eq.'E'.or.rep.eq.'e')goto580
        if(rep.eq.'S'.or.rep.eq.'s')goto590
        if(rep.eq.'B'.or.rep.eq.'b')goto800
        if(rep.eq.'Q'.or.rep.eq.'q')stop
                 goto600
c
c		lightning:
c570	write(6,'(''enter the x,y and z components of a vector '',
c     s	''of arbitrary module'',/,
c     s  '' in the direction of the source of light (def='',
c     s  3f6.2,'')'')')cl
570	write(6,'(2a,/,a,3f6.2,a)')
     s		'enter the x,y and z components of a vector ',
     s		'of arbitrary module',
     s		' in the direction of the source of light (def=',
     s		cl,')'
	read(5,'(a)')line
	if(line(1:1).ne.' ')read(line,*)cl
	goto600
c
c            plot:
500     if(nplv.le.0)then
                 write(6,*)' snake:nothing to plot'
                 goto600
        endif
c                horizontal axis:
        write(6,'(''$horizontal axis:(x,y,z,bx,by,bz,sx,sy,sz,'' 
     s      ,''l,B,r,t,rt or S def=y)'')')
        read(5,'(a)')rep2
        if(rep2.eq.' ')rep2='y'
        if(rep2.eq.'X' .or.rep2.eq.'x' )ih=1
        if(rep2.eq.'Y' .or.rep2.eq.'y' )ih=2
        if(rep2.eq.'Z' .or.rep2.eq.'z' )ih=3
        if(rep2.eq.'BX'.or.rep2.eq.'bx')ih=4
        if(rep2.eq.'BY'.or.rep2.eq.'by')ih=5
        if(rep2.eq.'BZ'.or.rep2.eq.'bz')ih=6
        if(rep2.eq.'SX'.or.rep2.eq.'sx')ih=7
        if(rep2.eq.'SY'.or.rep2.eq.'sy')ih=8
        if(rep2.eq.'SZ'.or.rep2.eq.'sz')ih=9
        if(rep2.eq.'L' .or.rep2.eq.'l' )ih=10
        if(rep2.eq.'B' .or.rep2.eq.'b' )ih=11
        if(rep2.eq.'R' .or.rep2.eq.'r' )ih=12
        if(rep2.eq.'T' .or.rep2.eq.'t' )ih=13
        if(rep2.eq.'RT'.or.rep2.eq.'rt')ih=14
        if(rep2.eq.'S' .or.rep2.eq.'s' )ih=15
c                vertical axis:
        write(6,'(''$vertical axis:(x,y,z,bx,by,bz,sx,sy,sz,'' 
     s      ,''l,B,r,t,rt or S def=x)'')')
        read(5,'(a)')rep2
        if(rep2.eq.' ')rep2='x'
        if(rep2.eq.'X' .or.rep2.eq.'x' )iv=1
        if(rep2.eq.'Y' .or.rep2.eq.'y' )iv=2
        if(rep2.eq.'Z' .or.rep2.eq.'z' )iv=3
        if(rep2.eq.'BX'.or.rep2.eq.'bx')iv=4
        if(rep2.eq.'BY'.or.rep2.eq.'by')iv=5
        if(rep2.eq.'BZ'.or.rep2.eq.'bz')iv=6
        if(rep2.eq.'SX'.or.rep2.eq.'sx')iv=7
        if(rep2.eq.'SY'.or.rep2.eq.'sy')iv=8
        if(rep2.eq.'SZ'.or.rep2.eq.'sz')iv=9
        if(rep2.eq.'L' .or.rep2.eq.'l' )iv=10
        if(rep2.eq.'B' .or.rep2.eq.'b' )iv=11
        if(rep2.eq.'R' .or.rep2.eq.'r' )iv=12
        if(rep2.eq.'T' .or.rep2.eq.'t' )iv=13
        if(rep2.eq.'RT'.or.rep2.eq.'rt')iv=14
        if(rep2.eq.'S' .or.rep2.eq.'s' )iv=15
        lbox=ih.le.3.and.iv.le.3
        if(lbox)call box
	nclic=0
        call pretra(*600)
	goto600
c
c                display a vector:
c
510     if(nuve.le.0)then
                 write(6,*)' snake:nothing to display'
                 goto600
        endif
        if(indved.lt.nuve)indved=indved+1
515     write(6,'(''$vector index (def='',i4,'')?'')')indved
        read(5,'(i3)')indve
        if(indve.le.0)then
                 indve=indved
        else
                 indved=indve
        endif
        if(indve.gt.nuve)then
                 write(6,*)' snake:vector index must be <',nuve+1
                 write(6,*)
                 indved=nuve
                 goto515
        endif
516       write(6,'(''$output end-plane index(def='',i4,'')?'')')nuep
          read(5,'(i3)')indep
          if(indep.le.0)indep=nuep
          if(indep.gt.nuep)then
                 write(6,*)' snake:end-plane index must be <',nuep+1
                 write(6,*)
                 goto516
        endif
c		find the region index (indre) of this e-p:
      indt=0
      do i=1,nure
        do j=1,nep(i)
		indt=indt+1
         if(indt.eq.indep)then
			indre=i
			jep=j
         endif
        enddo
      enddo
      inep=indep-jep
      write(6,'(4(a,i2),3a,/)')' e-pl #',indep,' is the e-pl #',jep,
     s     ' (over ',nep(indre),') of region #',indre,
     s     ' "',rname(indre),'"'
      write(6,'(a,i2,a,i2,a)')
     s  '            ep# 0              ep#',inep,
     s  '                     ep#',indep,'          unit'
      write(6,'(a,i2,a,i2)')
     s'            abs.         abs.      reg.#',indre,
     s '         abs.      reg.#',indre
      write(6,'(a5,5f13.6,a6)')
     s   'x',vea(indve,0,1),vea(indve,inep,1),veri(indve,indre,1),
     s		vea(indve,indep,1),vea(indve,indep,10),'mm',
     s   'y',vea(indve,0,2),vea(indve,inep,2),veri(indve,indre,2),
     s		vea(indve,indep,2),vea(indve,indep,11),'mm',
     s   'z',vea(indve,0,3),vea(indve,inep,3),veri(indve,indre,3),
     s		vea(indve,indep,3),vea(indve,indep,12),'mm',
     s   'cx',vea(indve,0,4),vea(indve,inep,4),veri(indve,indre,4),
     s		vea(indve,indep,4),vea(indve,indep,13),'-',
     s   'cy',vea(indve,0,5),vea(indve,inep,5),veri(indve,indre,5),
     s		vea(indve,indep,5),vea(indve,indep,14),'-',
     s   'cz',vea(indve,0,6),vea(indve,inep,6),veri(indve,indre,6),
     s		vea(indve,indep,6),vea(indve,indep,15),'-',
     s   'p',vea(indve,0,7),vea(indve,inep,7),veri(indve,indre,7),
     s		vea(indve,indep,7),vea(indve,indep,7),'Gev/c',
     s   'l',vea(indve,0,8),vea(indve,inep,8),veri(indve,indre,8),
     s		vea(indve,indep,8),vea(indve,indep,8),'mm',
     s   'live?',vea(indve,0,9),vea(indve,inep,9),veri(indve,indre,9),
     s		vea(indve,indep,9),vea(indve,indep,9),'-',
     s   'sx',vea(indve,0,16),vea(indve,inep,16),veri(indve,indre,10),
     s		vea(indve,indep,16),vea(indve,indep,19),'-',
     s   'sy',vea(indve,0,17),vea(indve,inep,17),veri(indve,indre,11),
     s		vea(indve,indep,17),vea(indve,indep,20),'-',
     s   'sz',vea(indve,0,18),vea(indve,inep,18),veri(indve,indre,12),
     s		vea(indve,indep,18),vea(indve,indep,21),'-'
      if(vea(indve,indep,9).ne.1)write(6,'(a,i3,2a,i3,a,/,1x,a)')
     s	' death in region no.',iredead(indve)
     s ,' when trying to reach the end'
     s ,' plane no.',iepdead(indve),' . diagnostic :',diag(indve)
        write(6,*)
        goto600
c
c                compute the envelope:
c
580     if(nuve.le.0)then
                 write(6,*)' snake:nothing to display'
                 goto600
        endif
586       write(6,'('' end-plane index(def='',i4,'')?'')')nuep
          read(5,'(i3)')indep
          if(indep.le.0)indep=nuep
          if(indep.gt.nuep)then
                 write(6,*)' snake:end-plane index must be <',nuep+1
                 write(6,*)
                 goto586
        endif
c		find the region index (indre) of this e-p:
      indt=0
      do i=1,nure
        do j=1,nep(i)
        indt=indt+1
        if(indt.eq.indep)then
			indre=i
			jep=j
        endif
        enddo
       enddo
      write(6,'(4(a,i2),3a)')' e-pl #',indep,' is the e-pl #',jep,
     s     ' (over ',nep(indre),') of region #',indre,
     s     ' "',rname(indre),'"'
      do j=1,8
		nlive=0
		xmax(j)=-1.e30
		xmin(j)=1.e30
        do i=1,nuve
            if(vea(i,indep,9).gt.0.)then
			nlive=nlive+1
			xval=vea(i,indep,j)
             if(xval.gt.xmax(j))xmax(j)=xval
             if(xval.lt.xmin(j))xmin(j)=xval
            endif
        enddo
       enddo
      write(6,*)nlive,' traj. still alive on this e-p'
      if(nlive.le.0)goto600
        write(6,581)
        write(6,584)'max',(xmax(i),i=1,8)
        write(6,584)'min',(xmin(i),i=1,8)
        write(6,582)
      do j=1,8
		xmax(j)=-1.e30
		xmin(j)=1.e30
        do i=1,nuve
            if(vea(i,indep,9).gt.0.)then
			xval=vea(i,indep,j+9)
            if(xval.gt.xmax(j))xmax(j)=xval
            if(xval.lt.xmin(j))xmin(j)=xval
            endif
         enddo
       enddo
        write(6,584)'max',(xmax(i),i=1,6)
        write(6,584)'min',(xmin(i),i=1,6)
        write(6,*)
        goto600
581     format(/,'absolute:',
     @     'x(mm)',5x,'y(mm)',5x,'z(mm)',4x,'cx',8x,'cy',8x,'cz'
     @   ,4x,'p(gev)   l(mm) ')
582     format(/,'relative:',
     @     'x(mm)',5x,'y(mm)',5x,'z(mm)',4x,'cx',8x,'cy',8x,'cz')
584     format(a3,3f10.3,3f10.6,f7.4,f9.1,f5.0)
c
c                replace the actual plot by its section:
c		(the section plane is the y abs.=0. plane)
c
590     if(nplv.le.0)then
                 write(6,*)' snake:nothing to cut'
                 goto600
        endif
      write(6,*)' >0. (>) , <0. (<)  or both (b) crossing (def=b)?'
      read(5,'(a)')rep
      scros=0.
      if(rep.eq.'>')scros=1.
      if(rep.eq.'<')scros=-1.
      write(6,*)' needle total size (def=1.mm)?'
      read(5,'(a)')line
      sneedle=1.
      if(line(1:1).ne.' ')read(line,*)sneedle
c  store the scatter plot in a temporary arrey named "plottn(2,maxplv,6)"
      sneedle=.5*sneedle
      nplvn=0
      iep=0
      do i=1,nure
		iep=iep+nep(i)
      enddo
c limit the search to the non old plot:
      do 591 iv=nplvold+1,nplv
c limit the search to vectors still alive at the last e-p:
		ivec=iv+nplvold
		if(vea(ivec,iep,9).ne.1)goto591
c limit the search to the non clip part:
      do 592 ip=iplpi(iv),iplpf(iv)-1
		pini=plott(ip,iv,2)
		pfin=plott(ip+1,iv,2)
        if((pini*pfin).le.0.)then
            if((pfin-pini)*scros.ge.0.)then
               if(nplvn.ge.maxplv)then
			    write(6,*)' no room for the whole plot'
			    goto593
                endif
			  nplvn=nplvn+1
			  den=1./(pfin-pini)
			  cini=pfin*den
			  cfin=-pini*den
		pintx=plott(ip,iv,1)*cini+plott(ip+1,iv,1)*cfin
		pintz=plott(ip,iv,3)*cini+plott(ip+1,iv,3)*cfin
			  plottn(1,nplvn,1)=pintx+sneedle
			  plottn(1,nplvn,2)=sneedle*sign(1.,pini)
			  plottn(1,nplvn,3)=pintz+sneedle
		dpintx=(plott(ip+1,iv,1)-plott(ip,iv,1))*den
		dpintz=(plott(ip+1,iv,3)-plott(ip,iv,3))*den
			  plottn(1,nplvn,4)=dpintx
	 	 	  plottn(1,nplvn,5)=0.
			  plottn(1,nplvn,6)=dpintz
			  plottn(2,nplvn,1)=pintx-sneedle
		  	  plottn(2,nplvn,2)=sneedle*sign(1.,pfin)
			  plottn(2,nplvn,3)=pintz-sneedle
			  plottn(2,nplvn,4)=dpintx
		 	  plottn(2,nplvn,5)=0.
			  plottn(2,nplvn,6)=dpintz
            endif
         endif
 592		continue
 591      continue
 593     if(nplvn.le.0)then
                 write(6,*)' snake:no cut, the initial plot is kept'
                 goto600
      else
		nplvold=0
		doiv=1,nplvn
		  call plotin(iv,plottn(1,iv,1),plottn(1,iv,2),
     s   plottn(1,iv,3),plottn(1,iv,4),plottn(1,iv,5),
     s   plottn(1,iv,6),0.,0.,0.,.true.)
		  call plotin(iv,plottn(2,iv,1),plottn(2,iv,2),
     s   plottn(2,iv,3),plottn(2,iv,4),plottn(2,iv,5),
     s   plottn(2,iv,6),0.,0.,0.,.false.)
		enddo
		call unclip
        endif
        goto600
c
c                output for mudifi:
c
520     write(6,'(''$input,output e-p index(def= 0,'',i3,'')?'')')nuep
        read(5,'(a)')line
        if(line(1:1).eq.' ')then
                 in1=0
                 in2=nuep
        else
c                read(line,'(2i3)')in1,in2
                  read(line,*)in1,in2
        endif
        write(6,'(''$ output-file name?(def='',a,'')'')')fimu
        read(5,'(a)')fimut
        if(fimut.ne.fibl)fimu=fimut
        write(6,*)fimu
        write(6,*)
        iw=0
        open(4,file=fimu,status='unknown')
        if(rep.eq.'U'.or.rep.eq.'u')then
        call userout(in1,in2)
        close(4)
        goto 600
      endif
c        write(4,'(a9,3x,a8)')daet,tyme
      call fdate(daet)
      write(4,'(a24)')daet
        do 521 i=1,nuve
                 alive=vea(i,in1,9)
                 if(alive.le.0.)goto521
                 alive=vea(i,in2,9)
                 if(alive.le.0.)goto521
                 iw=iw+1
                 tl1=vea(i,in1,8)
                 tl2=vea(i,in2,8)
c       build a set of independent variables:
                 xr1=      vea(i,in1,10)
                 axr1=asin(vea(i,in1,13))
                 zr1=      vea(i,in1,12)
                 azr1=asin(vea(i,in1,15))
		 pxr1=    (vea(i,in1,19))
		 pyr1=    (vea(i,in1,20))
		 pzr1=    (vea(i,in1,21))
                 xr2=      vea(i,in2,10)
                 axr2=asin(vea(i,in2,13))
                 zr2=      vea(i,in2,12)
                 azr2=asin(vea(i,in2,15))
		 pxr2=    (vea(i,in2,19))
		 pyr2=    (vea(i,in2,20))
		 pzr2=    (vea(i,in2,21))
                 dl=tl2-tl1
                 psq=vea(i,in1,7)
cc  temporary: convert to transport's conventions (cm , mr proj. , %):
C OFFERMANN PATCH (M, TAN(ANG), FRACTION) 11/30/95
C was really giving M, ANG, FRACTION corrected 11/12/01  -jjl
                 psqref=vea(1,0,7)                 

c                 psq=((psq-psqref)/psqref)*100.  

                 psq=((psq-psqref)/psqref)    !fractional value  
                 cy1r=vea(1,in1,14)
                 axr1r=atan2(vea(1,in1,13),cy1r)
                 cy1=vea(i,in1,14)

c                 xr1=(xr1-vea(1,in1,10))/10.
c
c lines commented out with c2 are needed for output of trajectory parameters
c relative to a reference.
c When removing the c2's comment out the line after. JJL 2/14/2008
c
c lines commented out with c1 are needed for output of trajectory parameters
c not relative to a reference.
c When removing the c1's comment out the line before with a c2. JJL 4/5/2007
c 
c2                 xr1=(xr1-vea(1,in1,10))/1000.  !meters
                   xr1=xr1/1000.  ! w/o reference traj meters

c2                 axr1=(atan2(vea(i,in1,13),cy1)-axr1r)   !theta-theta_ref
                  axr1=atan2(vea(i,in1,13),cy1)   !theta
                 axr1=tan(axr1)              ! 11/13/01 Tan(theta-theta_ref)

c                zr1=(zr1-vea(1,in1,12))/10.

c2                zr1=(zr1-vea(1,in1,12))/1000.   !meters
                  zr1=zr1/1000.   !w/o reference traj meters

cc        assume azr1r=axr2r=0. and azr1r and azr2r neglectable: not anymore! 11/12/01
c                 azr1=atan2(vea(i,in1,15),cy1)*1000.

                 azr1r=atan2(vea(1,in1,15),cy1r) !phi reference
c2                 azr1=atan2(vea(i,in1,15),cy1)-azr1r   !phi-phi reference
                   azr1=atan2(vea(i,in1,15),cy1)   !phi
                 azr1=tan(azr1)     !tan(phi-phi_ref)
                 cy2=vea(i,in2,14)
c                 xr2=(xr2-vea(1,in2,10))/10.
c2                 xr2=(xr2-vea(1,in2,10))/1000.  !meters
                   xr2=xr2/1000.  !meters w/o reference traj
                 cy2r=vea(1,in2,14)
                 axr2r=atan2(vea(1,in2,13),cy2r)  !theta_ref

c                 axr2=(atan2(vea(i,in2,13),cy2)-axr2r)*1000.

c2                 axr2=(atan2(vea(i,in2,13),cy2)-axr2r)   !theta-theta_ref
                   axr2=atan2(vea(i,in2,13),cy2)   !theta w/o reference traj
                 axr2=tan(axr2)

c                 zr2=(zr2-vea(1,in2,12))/10.

c2                 zr2=(zr2-vea(1,in2,12))/1000.    !meters
                   zr2=zr2/1000.    !meters w/o reference traj

c                 azr2=atan2(vea(i,in2,15),cy2)*1000.

                 azr2r=atan2(vea(1,in2,15),cy2r)   !phi_ref
c2                 azr2=atan2(vea(i,in2,15),cy2)-azr2r     !phi-phi_ref
                   azr2=atan2(vea(i,in2,15),cy2)     !phi w/o reference tra 
                 azr2=tan(azr2)                          !tan(phi-phi_ref)
                 tl1r=vea(1,in1,8)
                 tl2r=vea(1,in2,8)
                 dlr=tl2r-tl1r
c                 dl=-(dl-dlr)/10.
                 dl=-(dl-dlr)/1000.          !meters

523           write(4,'(10e14.6)')xr1,axr1,zr1,azr1
     s			         ,xr2,axr2,zr2,azr2
     s		                 ,dl,psq
c  temporary end

c523           write(4,'(16e14.6)')xr1,axr1,zr1,azr1,pxr1,pyr1,pzr1
c     s			         ,xr2,axr2,zr2,azr2,pxr2,pyr2,pzr2
c     s		                 ,dl,psq
521     continue
        close(4)
        write(6,522)iw
c522     format(
c     @  i4,' non dead vector(s) of 16 components each have been stored'
c     @  ,' ',/
c     @  ,' component   content   unit',/
c     @  ,'     1         xin       mm',/
c     @  ,'     2         x''in     rad',/
c     @  ,'     3         xin       mm',/
c     @  ,'     4         z''in     rad',/
c     @  ,'     5         pxin      - ',/
c     @  ,'     6         pyin      - ',/
c     @  ,'     7         pzin      - ',/
c     @  ,'     8         xout      mm',/
c     @  ,'     9         x''out    rad',/
c     @  ,'     10        xout      mm',/
c     @  ,'     11        z''out    rad',/
c     @  ,'     12        pxout     - ',/
c     @  ,'     13        pyout     - ',/
c     @  ,'     14        pzout     - ',/
c     @  ,'     15      lin-lout    mm',/
c     @  ,'     16         p     Gev/c',/
c     @  ,' x''=angle of the momentum with the plane x=0.',/
c     @  ,' z''=angle of the momentum with the plane z=0.',/
c     @  ,' coordinates are in relative frame depending on the region.')
 522     format(
     @  i4,' non dead vector(s) of 10 components each have been stored'
     @  ,' ',/
     @  ,' component   content   unit',/
     @  ,'     1         xin       m',/
     @  ,'     2      Tan(x'') in',/
     @  ,'     3         xin       m',/
     @  ,'     4      Tan(z'') in',/
     @  ,'     5         xout      m',/
     @  ,'     6      Tan(x'') out ',/
     @  ,'     7         xout      m',/
     @  ,'     8       Tan(z'') out',/
     @  ,'     9       lin-lout    m',/
     @  ,'     10         p       fraction',/
     @  ,' coordinates are in relative frame depending on the region.')

        goto 600
c
c                binary output for raytrace:
c
525     write(6,'(''$input,output e-p index(def= 0,'',i3,'')?'')')nuep
        read(5,'(a)')line
        if(line(1:1).eq.' ')then
                 in1=0
                 in2=nuep
        else
c                read(line,'(2i3)')in1,in2
                  read(line,*)in1,in2
        endif
        write(6,'(''$ output-file name?(def='',a,'')'')')firay
        read(5,'(a)')fimut
        if(fimut.ne.fibl)firay=fimut
        write(6,*)firay
        write(6,*)
        open(4,file=firay,status='unknown'
     s          ,form='unformatted')
        IW=0
        DO 524I=1,NUVE
            ALIVE=VEA(I,IN1,9)
            IF(ALIVE.LT.0)GOTO524
            ALIVE=VEA(I,IN2,9)
            IF(ALIVE.LT.0)GOTO524
            IW=IW+1
 524     CONTINUE
        WRITE(4)IW,0.
      write(6,*)iw,' traj.'
        do 526 i=1,nuve
                 alive=vea(i,in1,9)
                 if(alive.le.0.)goto526
                 alive=vea(i,in2,9)
                 if(alive.le.0.)goto526
C   Raytrace type output data (MeV,cm,mrad,angles between the ray 
C	projected on a plane and an axis):
            eneray=vea(i,in1,7)*1000.
            xray=vea(i,in1,10)*.1
            tray=atan2(vea(i,in1,13),vea(i,in1,14))*1000.
            yray=-vea(i,in1,12)*.1
            pray=atan2(-vea(i,in1,15),vea(i,in1,14))*1000.
            zray=vea(i,in1,11)*.1
            x0ray=vea(i,in2,10)*.1
            xsray=atan2(vea(i,in2,13),vea(i,in2,14))*1000.
            y0ray=-vea(i,in2,12)*.1
            ysray=atan2(-vea(i,in2,15),vea(i,in2,14))*1000.
         write(4)eneray,xray,tray,yray,pray,zray,0.,
     s           x0ray,xsray,y0ray,ysray,0.,0.
 526      continue
      close(4)
      goto 600	
c
c                clip the plot:
c
 560    write(6,'(''$x y or z clipping (def=x) ?'')')
        read(5,'(a)')rep
      iqb=1
      if(rep.eq.'Y'.or.rep.eq.'y')iqb=2
      if(rep.eq.'Z'.or.rep.eq.'z')iqb=3
        write(6,'(''$clipping value ?'')')
        read(5,*)qb
        write(6,'(''$ <'',f10.3,'' or >'',f10.3,''(def=<) ?'')')
     s       qb,qb
        read(5,'(a)')rep
      sup=-1.
      if(rep.eq.'>')sup=1.
      call clip(sup,iqb,qb)
      goto 600
c
c                move the plot data:
c
 550      write(6,'(a,i3,3a)')' defauts are those of region #(def=',
     s	indre,' "',rname(indre),'")  ?'
        read(5,'(a)')line
        if(line(1:1).ne.' ')read(line,*)indre
        write(6,'(''$translate (y/n def=y) ?'')')
        read(5,'(a)')rep
        if(rep.ne.'N'.and.rep.ne.'n')then
              write(6,'(4(a,f10.3))')' x,y and z components (def=',
     s        -xyz0(indre,1),' , ',-xyz0(indre,2),' , ',
     s	      -xyz0(indre,3),') ?'
        	read(5,'(a)')line
        	if(line(1:1).ne.' ')then
			read(line,*)xtp,ytp,ztp
        else
			xtp=-xyz0(indre,1)
			ytp=-xyz0(indre,2)
			ztp=-xyz0(indre,3)
        endif
        else
              xtp=0.
              ytp=0.
              ztp=0.
        endif
        write(6,'(''$rotation (y/n def=y) ?'')')
        read(5,'(a)')rep
        if(rep.ne.'N'.and.rep.ne.'n')then
              write(6,'(4(a,f8.3))')
     s		' rotation angles/ z,x and y in deg.(def=',
     s		-zang(indre)*rtod,' , ',-xang(indre)*rtod,' , ',
     s		-yang(indre)*rtod,') ?'
        	read(5,'(a)')line
        	if(line(1:1).ne.' ')then
                read(line,*)zad,xad,yad
           	   	xap=xad*dtor
             	 	yap=yad*dtor
              		zap=zad*dtor
        else
			zap=-zang(indre)
			xap=-xang(indre)
			yap=-yang(indre)
        endif
        else
              xap=0.
              yap=0.
              zap=0.
        endif
c
c  move trajectories:
c
      do 6 indve=1,nplv
          np=nplp(indve)
         do 6 n=1,np
           call arot(plott(n,indve,1),plott(n,indve,2),plott(n,indve,3)
     s              ,plott(n,indve,1),plott(n,indve,2),plott(n,indve,3)
     s              ,.true.)
           call arot(plott(n,indve,4),plott(n,indve,5),plott(n,indve,6)
     s              ,plott(n,indve,4),plott(n,indve,5),plott(n,indve,6)
     s              ,.false.)
6          call arot(plott(n,indve,7),plott(n,indve,8),plott(n,indve,9)
     s              ,plott(n,indve,7),plott(n,indve,8),plott(n,indve,9)
     s              ,.false.)
c
c  move free boxes:
c
        do 7 ir=1,nure
          do 7 it=1,8
7          call arot(boxp(it,ir,1),boxp(it,ir,2),boxp(it,ir,3)
     s              ,boxp(it,ir,1),boxp(it,ir,2),boxp(it,ir,3)
     s              ,.true.)
c
c  move region frames:
c
      do ir=1,nure
        do irs=1,nrs(ir)
        do iend=1,2
	   	  call arot(regp(iend,irs,ir,1),regp(iend,irs,ir,2),
     s		  regp(iend,irs,ir,3),regp(iend,irs,ir,1),
     s		  regp(iend,irs,ir,2),regp(iend,irs,ir,3),.true.)
        enddo
        enddo
       enddo
c
c  move end-planes:
c
      do iep=1,nuep
        do ieps=1,neps(iep)
        do iend=1,2
	   	  call arot(epp(iend,ieps,iep,1),epp(iend,ieps,iep,2),
     s		  epp(iend,ieps,iep,3),epp(iend,ieps,iep,1),
     s		  epp(iend,ieps,iep,2),epp(iend,ieps,iep,3),.true.)
        enddo
        enddo
      enddo
c
c  move user plot:
c
           do i=1,nulu
          	 do j=1,nuseu(i)
           call arot(plotus(j,i,1),plotus(j,i,2),plotus(j,i,3)
     s              ,plotus(j,i,1),plotus(j,i,2),plotus(j,i,3)
     s              ,.true.)
         	  end do
           end do
        call arot(cl(1),cl(2),cl(3),cl(1),cl(2),cl(3),.false.)
        goto 600
c
c                Bare plot (remove in the plot everything except the filaments)
c
 800     bare=.not.bare
      goto 600
c                to the beginning:
c
 400     continue
        goto 100
        end
        subroutine userout(in1,in2)
c*************************************************************************
c                                                                        *
c                       u s e r o u t                                    *
c                                                                        *
c*************************************************************************
c
c uere defined output of ray tracing data from end planes in1 and in2
        include 'snake.inc'
        common /absfra/ vea(maxve,0:maxep,21),nuve,nuep,nure
      iw=0
      do 521 i=1,nuve
                 alive=vea(i,in1,9)
                 if(alive.le.0.)goto521
                 alive=vea(i,in2,9)
                 if(alive.le.0.)goto521
                 iw=iw+1
                 tl1=vea(i,in1,8)
                 tl2=vea(i,in2,8)
c       built a set of independent variables:
                 xr1=      vea(i,in1,10)
                 axr1=asin(vea(i,in1,13))
                 zr1=      vea(i,in1,12)
                 azr1=asin(vea(i,in1,15))
		 pxr1=    (vea(i,in1,19))
		 pyr1=    (vea(i,in1,20))
		 pzr1=    (vea(i,in1,21))
                 xr2=      vea(i,in2,10)
                 axr2=asin(vea(i,in2,13))
                 zr2=      vea(i,in2,12)
                 azr2=asin(vea(i,in2,15))
		 pxr2=    (vea(i,in2,19))
		 pyr2=    (vea(i,in2,20))
		 pzr2=    (vea(i,in2,21))
                 dl=tl2-tl1
                 psq=vea(i,in1,7)
c  temporary: convert to transport's conventions (cm , mr proj. , %):
                 psqref=vea(1,0,7)
                 psq=((psq-psqref)/psqref)*100.
                 cy1r=vea(1,in1,14)
                 axr1r=atan2(vea(1,in1,13),cy1r)
                 cy1=vea(i,in1,14)
                 xr1=(xr1-vea(1,in1,10))/10.
                 axr1=(atan2(vea(i,in1,13),cy1)-axr1r)*1000.
                 zr1=(zr1-vea(1,in1,12))/10.
c               assume azr1r=axr2r=0. and azr1r and azr2r negligable:
                 azr1=atan2(vea(i,in1,15),cy1)*1000.
                 cy2=vea(i,in2,14)
                 xr2=(xr2-vea(1,in2,10))/10.
                 cy2r=vea(1,in2,14)
                 axr2r=atan2(vea(1,in2,13),cy2r)
                 axr2=(atan2(vea(i,in2,13),cy2)-axr2r)*1000.
                 zr2=(zr2-vea(1,in2,12))/10.
                 azr2=atan2(vea(i,in2,15),cy2)*1000.
                 tl1r=vea(1,in1,8)
                 tl2r=vea(1,in2,8)
                 dlr=tl2r-tl1r
                 dl=-(dl-dlr)/10.
c  temporary end
523           write(4,'(10e14.6)')xr1,axr1,zr1,azr1
     s			         ,xr2,axr2,zr2,azr2
     s		                 ,dl,psq
521     continue
c      save /absfra/
      return
      end
c
c
c
        subroutine spring
c*************************************************************************
c                                                                        *
c                       s p r i n g                                      *
c                                                                        *
c*************************************************************************
c
c                create the initial vector(s) set (end-plane index=0)
c     use 6 nested loops whose parameters (initial,maximum and increment value)
c         are given by the user
c       output:
c             nuve:number of input vector(s) defined
c             vea(iv,0,ic):description of the particules at each end-plane
c                iv<=nuve<=maxve:vector index
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
c             nplv<=maxplv:number of trajectories to be plotted
c             nplp(traj.index)<=maxplp:number of submits in the plot-line
c             plott(1,traj.index,x y or z):absol.coord.of submit
c
        include 'snake.inc'
        character rep*1,filtra*20,fitr*20,fibl*20,line*60
        logical lrand,ldisk,lref
c                next line : begening of the /plot/ package
        logical lbox,luser,bare
      integer ul
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
c        data fitr/'sn-test4-traj.dat'/
c       data fitr/'sn-900-traj.dat'/
c        data fitr/'sol_traj.dat'/
c        data fitr/'td_traj2.dat'/
        data fitr/'TRAJ.DAT'/
        data fibl/' '/
 400     nuve=0
 300     format(3f13.6/)
 330     format(3f13.6)
c
c                data from terminal or disk ?
c
        write(6,'(a)')'$ data from t(erminal) or d(isk) (def.=t):'
        read(5, '(a)')rep
        ldisk=rep.eq.'d'.or.rep.eq.'D'
        if(ldisk)then
c
          write(6,'(''$ trajectory-file name?(def='',a,'')'')')fitr
          read(5,'(a)')filtra
          if(filtra.ne.fibl)fitr=filtra
          write(6,*)fitr
        ul=7	
          open(ul,file=fitr,status='unknown')
      else
	  ul=5
        endif
c
c
c                loop random or individual ?
c
        write(6,'(a)')'$ l(oops) r(andom) or i(ndividual) (def.=l):'
        read(ul, '(a)')rep
        if(ldisk)write(6,*)rep
      if(rep.eq.' ')rep='l'
        lrand=rep.eq.'R'.or.rep.eq.'r'
        if(lrand)then
                 write(6,'(a)')'$ number of trajectories :'
                 read(ul,*)istat
                if(ldisk)write(6,*)istat
               open(8,file='sn-rand.dat',status='old')
                 read(8,*)inidef
                 close(8)
                 write(6,'(a,i15,a)')
     s           '$ initial value of random generator (def from disk='
     s           ,inidef,'):'
                 read(ul,'(i15)')iniran
                 if(iniran.eq.0)iniran=inidef
                 write(6,*)iniran
                 call ranset(iniran)
        endif
c
c******************************
c			case rep.e.'i' : loop or random :
	if(rep.ne.'I'.and.rep.ne.'i')then
c
c       input min,max and step from terminal:
c                x:
        write(6,'(''$ x(def=0.,0.,1.mm):min,max,step?'')')
        read(ul,'(a)')line
      if(line(1:1).ne.' ')then
		read(line,*)xmin,xmax,xstep
      else
		xmin=0.
		xmax=0.
		xstep=0.
      endif
        if(xstep.eq.0.)xstep=1.
        if(((xmax-xmin)*xstep).lt.0.)xmax=xmin
        write(6,330)xmin,xmax,xstep
c
c                y:
        write(6,'(''$ y(def=0.,0.,1.mm):min,max,step?'')')
        read(ul,'(a)')line
      if(line(1:1).ne.' ')then
		read(line,*)ymin,ymax,ystep
      else
		ymin=0.
		ymax=0.
		ystep=0.
      endif
        if(ystep.eq.0.)ystep=1.
        if(((ymax-ymin)*ystep).lt.0.)ymax=ymin
        write(6,330)ymin,ymax,ystep
c
c                z:
        write(6,'(''$ z(def=0.,0.,1.mm):min,max,step?'')')
        read(ul,'(a)')line
      if(line(1:1).ne.' ')then
		read(line,*)zmin,zmax,zstep
      else
		zmin=0.
		zmax=0.
		zstep=0.
      endif
        if(zstep.eq.0.)zstep=1.
        if(((zmax-zmin)*zstep).lt.0.)zmax=zmin
        write(6,330)zmin,zmax,zstep
c
c                theta x:
        write(6,'(''$ theta x (proj.angle)''
     s   ,''(def=0.,0.,1.rad):min,max,step?'')')
        read(ul,'(a)')line
      if(line(1:1).ne.' ')then
		read(line,*)txmin,txmax,txstep
      else
		txmin=0.
		txmax=0.
		txstep=0.
      endif
        if(txstep.eq.0.)txstep=1.
        if(((txmax-txmin)*txstep).lt.0.)txmax=txmin
        write(6,330)txmin,txmax,txstep
c
c                theta z:
        write(6,'(''$ theta z (space angle)''
     s   ,''(def=0.,0.,1.rad):min,max,step?'')')
        read(ul,'(a)')line
      if(line(1:1).ne.' ')then
		read(line,*)tzmin,tzmax,tzstep
      else
		tzmin=0.
		tzmax=0.
		tzstep=0.
      endif
        if(tzstep.eq.0.)tzstep=1.
        if(((tzmax-tzmin)*tzstep).lt.0.)tzmax=tzmin
        write(6,330)tzmin,tzmax,tzstep
c
c                p/q:
        write(6,'(''$ p/q(def=1.,1.,1.gev/c):min,max,step?'')')
        read(ul,'(a)')line
      if(line(1:1).ne.' ')then
		read(line,*)psqmin,psqmax,psqstep
      else
		psqmin=1.
		psqmax=1.
		psqstep=1.
      endif
        if(psqstep.eq.0.)psqstep=1.
        if(((psqmax-psqmin)*psqstep).lt.0.)psqmax=psqmin
        write(6,330)psqmin,psqmax,psqstep

c
c                spin x:
        write(6,'(''$ spin x''
     s   ,''(def=0.,0.,1.):min,max,step?'')')
        read(ul,'(a)')line
      if(line(1:1).ne.' ')then
		read(line,*)spxmin,spxmax,spxstep
      else
		spxmin=0.
		spxmax=0.
		spxstep=0.
      endif
        if(spxstep.eq.0.)spxstep=1.
        if(((spxmax-spxmin)*spxstep).lt.0.)spxmax=spxmin
        write(6,330)spxmin,spxmax,spxstep
c
c                spin y:
        write(6,'(''$ spin y''
     s   ,''(def=0.,0.,1.):min,max,step?'')')
        read(ul,'(a)')line
      if(line(1:1).ne.' ')then
		read(line,*)spymin,spymax,spystep
      else
		spymin=0.
		spymax=0.
		spystep=0.
      endif
        if(spystep.eq.0.)spystep=1.
        if(((spymax-spymin)*spystep).lt.0.)spymax=spymin
        write(6,330)spymin,spymax,spystep
c
c                spin z:
        write(6,'(''$ spin z''
     s   ,''(def=0.,0.,1.):min,max,step?'')')
        read(ul,'(a)')line
      if(line(1:1).ne.' ')then
		read(line,*)spzmin,spzmax,spzstep
      else
		spzmin=0.
		spzmax=0.
		spzstep=0.
      endif
        if(spzstep.eq.0.)spzstep=1.
        if(((spzmax-spzmin)*spzstep).lt.0.)spzmax=spzmin
        write(6,330)spzmin,spzmax,spzstep
c
c                         insert ref. traj. ?
        write(6,'(a)')'$ insert ref. traj. ? (y or n , def.=y)'
        read(ul,'(a)')rep
        if(ldisk)write(6,*)rep
        lref=rep.ne.'N'.and.rep.ne.'n'
c
c
c                         o k ?
        write(6,'(a)')'$ o.k. ? (y or n , def.=y)'
        read(5,'(a)')rep
        if(rep.eq.'N'.or.rep.eq.'n')goto400
c
c  if random , then compute the loop param. to obtain istat*1*1*1*1*1=istat
c          traj. :
c
        if(lrand)then
                 ixmax=istat
                 iymax=1
                 izmax=1
                 itxmax=1
                 itzmax=1
                 ipsqmax=1
                 ispxmax=1
                 ispymax=1
                 ispzmax=1
c
c   else compute the loop parameter :
c
        else
           ixmax=(xmax-xmin)/xstep+1.5
           iymax=(ymax-ymin)/ystep+1.5
           izmax=(zmax-zmin)/zstep+1.5
           itxmax=(txmax-txmin)/txstep+1.5
           itzmax=(tzmax-tzmin)/tzstep+1.5
           ipsqmax=(psqmax-psqmin)/psqstep+1.5
           ispxmax=(spxmax-spxmin)/spxstep+1.5
           ispymax=(spymax-spymin)/spystep+1.5
           ispzmax=(spzmax-spzmin)/spzstep+1.5
        endif
c
c
c                creation of the input vectors:
c
c
c               insert reference trajectory in 1st position ?
c
        if(lref)then
	     call inject(0.,0.,0.,0.,1.,0.,(psqmin+psqmax)*.5,
     s			0.,1.,0.)
        endif
c
c          loop on the 9 indep. components(x,y,z,tx,tz,p,px,py,pz):
c
c                x:
        do 1 ix=1,ixmax
          x=(ix-1)*xstep+xmin
c                  y:
          do 2 iy=1,iymax
            y=(iy-1)*ystep+ymin
c                    z:
            do 3 iz=1,izmax
              z=(iz-1)*zstep+zmin
c                      theta x:
              do 4 itx=1,itxmax
                tx=(itx-1)*txstep+txmin
c                         theta z:
                do 5 itz=1,itzmax
                  tz=(itz-1)*tzstep+tzmin
c                          p:
                     do 6 ipsq=1,ipsqmax
                       psq=(ipsq-1)*psqstep+psqmin
c
                    	 do 7 ispx=1,ispxmax
                     	  px=(ispx-1)*spxstep+spxmin
c
                    	    do 8 ispy=1,ispymax
                     	     py=(ispy-1)*spystep+spymin
c
                              do 9 ispz=1,ispzmax
                     	       pz=(ispz-1)*spzstep+spzmin
                       if(lrand)then
                         x=aphasa(xmin,xmax)
                         y=aphasa(ymin,ymax)
                         z=aphasa(zmin,zmax)
                         tx=aphasa(txmin,txmax)
                         tz=aphasa(tzmin,tzmax)
                         psq=aphasa(psqmin,psqmax)
			 px=aphasa(spxmin,spxmax)
			 py=aphasa(spymin,spymax)
			 pz=aphasa(spzmin,spzmax)
                       endif
                       cz=sin(tz)
                       cx=sin(tx)*cos(tz)
                       cy2=1.-cx*cx-cz*cz
                       if(cy2.lt.0.)stop' spring:bad angles'
c                  take the positive solution:
                       cy=sqrt(cy2)
c                         index of the new vector:
                          if(nuve.ge.maxve)then
        write(6,*)' room for only',maxve,' vector(s)'
        write(6,*)'  first ignored vector:'
        write(6,100)x,y,z,tx,tz,psq,px,py,pz
100     format(' x=',f10.3,' y=',f10.3,' z=',f10.3,' tx=',f10.6,
     @     ' tz=',f10.6,' p=',f8.3,' sx=',f8.3,
     @     ' sy=',f8.3,' sz=',f8.3)
                                   goto10
			   endif
      call inject(x,y,z,cx,cy,cz,psq,px,py,pz)
     @		    
9			   continue
8             continue
7		      continue
6                  continue
5               continue
4             continue
3           continue
2         continue
1       continue
10       continue
        if(lrand)then
          call ranget(iniran)
          write(6,'(a,i15/a)')' final value of the random generator :'
     s   ,iniran,' ( saved on disk file ''snake-rand.dat'')'
          open(8,file='sn-rand.dat',status='old')
          write(8,*)iniran
          close(8)
      endif
c***********************
      else
c**********************
c			case rep.eq.'i' : individual
                 write(6,'(a)')'$ number of trajectories :'
                 read(ul,*)istat
                if(ldisk)write(6,*)istat
       do 108,ix=1,istat
       if(.not. ldisk)write(6,*)' type in x, y, z, tx, tz, psq',
     &                         'px, py, pz    for ray #',ix
       read(ul,*)x,y,z,tx,tz,psq
       if(ldisk)write(6,100)x,y,z,tx,tz,psq
                       cz=sin(tz)
                       cx=sin(tx)*cos(tz)
                       cy2=1.-cx*cx-cz*cz
                       if(cy2.lt.0.)stop' spring:bad angles'
c                  take the positive solution:
                       cy=sqrt(cy2)
c                         index of the new vector:
                          if(nuve.ge.maxve)then
        write(6,*)' room for only',maxve,' vector(s)'
        write(6,*)'  first ignored vector:'
        write(6,100)x,y,z,tx,tz,psq,px,py,pz
                                   goto107
                          endif
      call inject(x,y,z,cx,cy,cz,psq,px,py,pz)
 108    continue
      endif
c**************************
c
107     if(ldisk)close(ul)
c      save /cplot/,/plot/,/absfra/,/crelfr/,/relfra/
        return
        end
c
        subroutine inject(x,y,z,cx,cy,cz,psq,px,py,pz)
c*************************************************************************
c                                                                        *
c                       i n j e c t                                      *
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
c
                          nuve=nuve+1
                          indve=nuve
			  ilive=nuve
c       fill the abs. coodr. buff. with this vector,region index=0
                          indre=0
                          vea(indve,0,1)=x
                          vea(indve,0,2)=y
                          vea(indve,0,3)=z
                          vea(indve,0,4)=cx
                          vea(indve,0,5)=cy
                          vea(indve,0,6)=cz
                          vea(indve,0,7)=psq
                          vea(indve,0,8)=0.
                          vea(indve,0,9)=1.
                          vea(indve,0,10)=x
                          vea(indve,0,11)=y
                          vea(indve,0,12)=z
                          vea(indve,0,13)=cx
                          vea(indve,0,14)=cy
                          vea(indve,0,15)=cz
                          vea(indve,0,16)=px
                          vea(indve,0,17)=py
                          vea(indve,0,18)=pz
                          vea(indve,0,19)=px
                          vea(indve,0,20)=py
                          vea(indve,0,21)=pz
c         fill the plot buffer with absolute data if room:
      call plotin(indve,x,y,z,0.,0.,0.,
     s				px,py,pz,.true.)
c      save /absfra/,/crelfr/,/relfra/
      return
      end
c
        subroutine send(kep,yep,tep,iep)
c*************************************************************************
c                                                                        *
c                       s e n d                                          *
c                                                                        *
c*************************************************************************
c                send the vectors from one end-plane to the next one
c       according to method.control the death of the vectors
c       input:
c             nuve:number of input vector(s) defined
c             method(region index):methode to use to carry particules
c                in this region
c             itype(region index):type in this method
c             fact(region index):multipl.factor of the field or ref.mommentum
c       input/output:
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
c
        include 'snake.inc'
	parameter (maxmir=10)
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
c                next line : begening of the /cer/ package
      common/cer/nmir,rmirin,rmirout,dzmir,dxmir,qmir(3,maxmir)
     s		,tiltz(maxmir),tiltx(maxmir),axp,axm,azp,azm
      common/mir/tv(3,maxmir),sv(3,maxmir),cv(3,maxmir)
c                previous line : end of the /cer/ package
c                next line : begening of the /diag/ package
      character diag*100,diagt*100
      common/cdiag/diag(maxve),diagt
      common/ndiag/iepdead(maxve),iredead(maxve)
c                previous line : end of the /diag/ package
      character mdiag*3,qdiag(6)*1
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
         logical eject,lcer
         common/alter/jflag,prec,hstmax,hstmin,sinc
         common/params/pin(10),fin,ctep,step,idir,epsil,pout(10)
     s   ,h,ifail,pal
         common/matrix/a(6,6),dy,s(3,3)
        dimension q(6)
        equivalence (x,q(1)),(y,q(2)),(z,q(3)),(xp,q(4)),(yp,q(5))
     @             ,(zp,q(6))
      data qdiag/'x','y','z','r','t','z'/
      dimension xyz(3),pm(5,3)
c                data for rungk:
c       end-plane orthogonal to axis y (1:x,2:y,3:z):
        idir=kep
c       accuracy of the final value of the coord. idir:(mm)
        epsil=.001
c       variable step (jflag=1):
        jflag=1
c   data for step change:(adapted for spectro.900:r=1800.,half gap=60.)
c       threshold value of relative field change for step change:
        prec=.02
        ctep=cos(tep)
        step=sin(tep)
        fin=yep*ctep
c			cerenkov mirrors?:
c
      lcer=cerfiln(iep,indre).ne.'none'
      if(.not.lcer)goto 197
      open(9,file=cerfiln(iep,indre),status='old')
      read(9,*)nmir,rmirin,rmirout,dzmir,dxmir
c next temp:
c	write(6,*)nmir,rmirin,rmirout,dzmir,dxmir
      if(nmir.gt.maxmir)then
		nmir=maxmir
		write(6,*)'too many mirrors,max=',maxmir
      endif
      read(9,*)axp,axm,azp,azm
c next temp:
      write(6,*)axp,axm,azp,azm
      axp=axp*dtor
      axm=axm*dtor
      azp=azp*dtor
      azm=azm*dtor
      do 198 imir=1,nmir
        read(9,*)(qmir(iq,imir),iq=1,3),tiltz(imir),tiltx(imir)
c next temp:
c		write(6,*)(qmir(iq,imir),iq=1,3),tiltz(imir),tiltx(imir)
      tiltz(imir)=tiltz(imir)*dtor
      tiltx(imir)=tiltx(imir)*dtor
c   I call "middle" of a mirror the barycenter of its 4 corners
c  distance (qmir="middle" of the mirror)_(cv=center of the sphere):
      xm=sqrt(rmirin**2-.25*(dxmir**2+dzmir**2))
c sv and tv are two unit vectors, perpendicular between them. s is normal to 
c the planes limiting the mirror in z'',t is normal to the planes limiting the 
c mirrors in x''.(x'',z'') coincidate with (xrel.,zrel.) if
c  tiltx=tiltz=0. To go from rel. frame to '' frame one rotates first /z, 
c  second /x
      cosz=cos(tiltz(imir))
      sinz=sin(tiltz(imir))
      cosx=cos(tiltx(imir))
      sinx=sin(tiltx(imir))
      tv(1,imir)=cosz
      tv(2,imir)=sinz*cosx
      tv(3,imir)=sinz*sinx
      sv(1,imir)=0.
      sv(2,imir)=-sinx
      sv(3,imir)=cosx
c  unit vector between "middle" and center:
      cv(1,imir)=tv(2,imir)*sv(3,imir)-tv(3,imir)*sv(2,imir)
      cv(2,imir)=tv(3,imir)*sv(1,imir)-tv(1,imir)*sv(3,imir)
      cv(3,imir)=tv(1,imir)*sv(2,imir)-tv(2,imir)*sv(1,imir)
      do 932 iq=1,3
c  compute the position ( cv ) of the center of this mirror:
      cv(iq,imir)=qmir(iq,imir)+cv(iq,imir)*xm
c  plot this mirror as a rectangle whose submits are those of the inner
c  (reflecting) surface of the mirror:
      pm(1,iq)=qmir(iq,imir)+.5*(tv(iq,imir)*dxmir+sv(iq,imir)*dzmir)
      pm(2,iq)=qmir(iq,imir)+.5*(tv(iq,imir)*dxmir-sv(iq,imir)*dzmir)
      pm(3,iq)=qmir(iq,imir)+.5*(-tv(iq,imir)*dxmir-sv(iq,imir)*dzmir)
      pm(4,iq)=qmir(iq,imir)+.5*(-tv(iq,imir)*dxmir+sv(iq,imir)*dzmir)
c    s   1         2         3         4         5         6         7
      pm(5,iq)=pm(1,iq)
 932      continue	
		nulu=nulu+1
      if(nulu.gt.nulumax)then
	  	  nulu=nulumax
		  write(6,*)' warning send: too many user line'
		  goto198
      endif
		luser=.true.
      nuseu(nulu)=5
      if(nuseu(nulu).gt.nuseumax)nuseu(nulu)=nuseumax
      do 400 is=1,nuseu(nulu)
 400      call plorot(pm(is,1),pm(is,2),pm(is,3),
     @		plotus(is,nulu,1),plotus(is,nulu,2),plotus(is,nulu,3),
     @		xb,yb,zb)
 198      continue
      close(9)
 197      continue
c
c                loop over the nuve vectors:
c      write(6,*)'in send'
c
        do 1 indve=1,nuve
      i=indve
c       dead?
        if(ver(i,9).le.0.)then
c       in this case:fill vect. with zero value:
                 do 2 j=1,8
                          ver(i,j)=0.
2                continue
                do 3 j=10,12
                          ver(i,j)=0.
3                continue
c                if it was dead or buried,now it is buried:
                 ver(i,9)=-1.
                 goto1
        endif
      diag(indve)=' '
        x=ver(i,1)
        y=ver(i,2)
        z=ver(i,3)
        cx=ver(i,4)
        cy=ver(i,5)
        cz=ver(i,6)
        psq=ver(i,7)
      px=ver(i,10)
      py=ver(i,11)
      pz=ver(i,12)
c               input point in the free-box?
      if(cyl(indre))then
		xyz(1)=sqrt(x*x+y*y)
		xyz(2)=atan2(y,x)			
      else
		xyz(1)=x
		xyz(2)=y
      endif
		xyz(3)=z
        eject=.false.
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
      write(diag(indve),'(3a)')'injected outside the free box in :'
     s			,qdiag(idiag),mdiag
      goto 200
      endif
c
c	cerenkov? If yes, emitt light from the current point with the
c  current direction, i.e. the radiator is the previous end-plane:
c
      if(lcer)then
c   x,y,z: point de depart de la lumiere cerenkov
c   cx,cy,cz: direction de depart de la lumiere cerenkov
c    x1,y1,z1: point de reflexion de la lumiere par le miroir "imir"
c   (pas de reflexion: imir=0, x1...=x...)
c    x2,y2,z2: point d'interception de la lumiere par le miroir "imir1"
c   (pas d'interception: imir1=0, x2...=x...)
      call mirror(x,y,z,cx,cy,cz,x1,y1,z1,x2,y2,z2,imir,imir1)
          call plorot(x,y,z,xa,ya,za,xb,yb,zb)
          call plotin(indve,xa,ya,za,0.,0.,0.,0.,0.,0.,.false.)
          call plorot(x1,y1,z1,xa1,ya1,za1,xb,yb,zb)
          call plotin(indve,xa1,ya1,za1,0.,0.,0.,0.,0.,0.,.false.)
          call plorot(x2,y2,z2,xa2,ya2,za2,xb,yb,zb)
          call plotin(indve,xa2,ya2,za2,0.,0.,0.,0.,0.,0.,.false.)
c  now, plot back to the starting point:
          call plotin(indve,xa1,ya1,za1,0.,0.,0.,0.,0.,0.,.false.)
          call plotin(indve,xa,ya,za,0.,0.,0.,0.,0.,0.,.false.)
      endif
c
c       no field,matrix or rung-khuta type?
c
        goto(101,102,103),method(indre)
c
c                no field:
c
 101     continue
c        write(6,*)'no field'
		yp=-x*step+y*ctep
        cyp=-cx*step+cy*ctep
        if(cyp.le.0.)goto15
        gly=(fin-yp)/cyp
        if(gly.lt.0.)goto15
        xt=x+gly*cx
        yt=y+gly*cy
        zt=z+gly*cz
c                        output point in the free box?
        if(xt.gt.frbox(1,indre,1).and.xt.lt.frbox(2,indre,1).and.
     @     yt.gt.frbox(1,indre,2).and.yt.lt.frbox(2,indre,2).and.
     @     zt.gt.frbox(1,indre,3).and.zt.lt.frbox(2,indre,3))then
                 x=xt
                 y=yt
                 z=zt
                 pal=gly
c test for collimator clearance
c      write(6,*)' calling collitest ',x,z
       if(colli(iep,indre,1).ne.0.) call collitest(x,z,iep,i,eject)
c      write(6,*)' return from collitest eject=',eject
                 goto199
        endif
c                find the escape point:
cbug 11/06/92:        glxz=(frbox(2,indre,2)-y)/cy
      glxz=1.e30
        if(xt.le.frbox(1,indre,1))then
		glxt=(frbox(1,indre,1)-x)/cx
        if(glxt.lt.glxz)then
			glxz=glxt
			idiag=1
			mdiag='min'
        endif
        elseif(xt.ge.frbox(2,indre,1))then
		glxt=(frbox(2,indre,1)-x)/cx
		if(glxt.lt.glxz)then
			glxz=glxt
			idiag=1
			mdiag='max'
        endif
        endif
        if(yt.le.frbox(1,indre,2))then
		glxt=(frbox(1,indre,2)-y)/cy
		if(glxt.lt.glxz)then
			glxz=glxt
			idiag=2
			mdiag='min'
        endif
        elseif(yt.ge.frbox(2,indre,2))then
		glxt=(frbox(2,indre,2)-y)/cy
        if(glxt.lt.glxz)then
			glxz=glxt
			idiag=2
			mdiag='max'
        endif
        endif
        if(zt.le.frbox(1,indre,3))then
		glxt=(frbox(1,indre,3)-z)/cz
        if(glxt.lt.glxz)then
			glxz=glxt
			idiag=3
			mdiag='min'
        endif
        elseif(zt.ge.frbox(2,indre,3))then
		glxt=(frbox(2,indre,3)-z)/cz
        if(glxt.lt.glxz)then
			glxz=glxt
			idiag=3
			mdiag='max'
        endif
        endif
        x=x+glxz*cx
        y=y+glxz*cy
        z=z+glxz*cz
        pal=glxz
        eject=.true.
c		built the diag:
	write(diag(indve),'(4a)')'exits free box in :'
     s			,qdiag(idiag),mdiag,' in a no field region'
        goto 199
c                can't reach the end-plane:
15      pal=0.
        eject=.true.
c		built the diag:
      write(diag(indve),'(2a)')'wrong direction to reach the end plane'
     s           , ' in a no field region'
      goto 200
c
c			cerenkov mirrors?:
c
 199      continue
      if(cerfiln(iep,indre).eq.'none')goto 200
      goto 200
c
c                 matrix:
c
 102     continue
        tx=asin(cx)
        tz=asin(cz)
        if(fact(indre).eq.0.)stop'send/1stmatrix:ref. momentum=0.'
        del=psq/fact(indre)-1.
c  length of the reference trajectory is taken in the directive file
c                (first auxiliary data) :
        pal=adata(1,indre)
        goto(301,302),itype(indre)
c                1st order matrix:
 301     xt=x*a(1,1)+tx*a(2,1)+z*a(3,1)+tz*a(4,1)+del*a(6,1)
        txt=x*a(1,2)+tx*a(2,2)+z*a(3,2)+tz*a(4,2)+del*a(6,2)
        zt=x*a(1,3)+tx*a(2,3)+z*a(3,3)+tz*a(4,3)+del*a(6,3)
        tzt=x*a(1,4)+tx*a(2,4)+z*a(3,4)+tz*a(4,4)+del*a(6,4)
        pl=x*a(1,5)+tx*a(2,5)+z*a(3,5)+tz*a(4,5)+del*a(6,5)
        pal=pal+pl
        cx=sin(txt)
        cz=sin(tzt)
        x=xt
        z=zt
        y=yep
      pxt=px*s(1,1)+py*s(2,1)+pz*s(3,1)
      pyt=px*s(1,2)+py*s(2,2)+pz*s(3,2)
      pzt=px*s(1,3)+py*s(2,3)+pz*s(3,3)
      px=pxt
      py=pyt
      pz=pzt
        cy2=1.-cx*cx-cz*cz
        if(cy2.lt.0.)then
                 eject=.true.
c		built the diag:
      write(diag(indve),'(2a)')'crazy direction after a 1st order'
     s           , ' matrix'
         goto 200
        endif
c                  take the positive solution:
        cy=sqrt(cy2)
c              output point in the free-box?
        do 12 k=1,3
      if(q(k).lt.frbox(1,indre,k))then
		idiag=k
		mdiag='min'
		eject=.true.
      endif
      if(q(k).gt.frbox(2,indre,k))then
		idiag=k
		mdiag='max'
		eject=.true.
      endif
 12      continue
c		built the diag:
      if(eject)then
      write(diag(indve),'(4a)')'exit free box in :'
     s			,qdiag(idiag),mdiag,' after 1st order matrix'
      goto 200
      endif
c test for collimator clearance
c      write(6,*)' calling collitest ',x,z
       if(colli(iep,indre,1).ne.0.) call collitest(x,z,iep,i,eject)
c      write(6,*)' return from collitest eject=',eject
      goto 200
c
c                2nd order matrix
c
302     stop'send/2nd matrix:not yet implemented'
c
c                rung-khuta:
c


c       fill the standard rungk buffer:
103     continue
        pin(1)=x
        pin(2)=y
        pin(3)=z
        pin(4)=cx
        pin(5)=cy
        pin(6)=cz
        pin(7)=px
        pin(8)=py
        pin(9)=pz
        pin(10)=psq
        h=hstmin
c		write(6,*)'calling rungk'
        call rungk(eject,i)
        x=pout(1)
        y=pout(2)
        z=pout(3)
        cx=pout(4)
        cy=pout(5)
        cz=pout(6)
        px=pout(7)
        py=pout(8)
        pz=pout(9)
c  next:sync.rad.aver.ener.loss:
      psq=pout(10)
      if(eject)goto 200
c test for collimator clearance
c      write(6,*)' calling collitest ',x,z
       if(colli(iep,indre,1).ne.0.) call collitest(x,z,iep,i,eject)
c      write(6,*)' return from collitest eject=',eject
c
c
 200     continue
        ver(i,1)=x
        ver(i,2)=y
        ver(i,3)=z
        ver(i,4)=cx
        ver(i,5)=cy
        ver(i,6)=cz
c       (no change in ver(i,7)=psq)
c  next:sync.rad.aver.ener.loss:
      ver(i,7)=psq
        ver(i,8)=ver(i,8)+pal
        if(eject)then
		ver(i,9)=0.
c      		 (if not:no change)
c at this point indep is the index of the previous end plane :
		iepdead(i)=indep+1
		iredead(i)=indre
		diag(i)=diagt
        endif
        ver(i,10)=px
        ver(i,11)=py
        ver(i,12)=pz
c
c                to next vector
c
 1       continue
c
c      save /cplot/,/plot/,/cer/,/mir/,/cdiag/,/ndiag/,/absfra/,/crelfr/,
c     &     /relfra/,/fidef/,/cfidef/,/alter/,/params/,/matrix/
        return
        end
c
c*************************************************************************
c                                                                        *
c                       c o l l i t e s t                                *
c                                                                        *
c*************************************************************************
c
      subroutine collitest(x,z,iep,iv,eject)
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
c                next line : begening of the /diag/ package
      character diag*100,diagt*100
      common/cdiag/diag(maxve),diagt
      common/ndiag/iepdead(maxve),iredead(maxve)
c                previous line : end of the /diag/ package
        logical eject
        real k,m1,m2
c 
c  there are two kinds of collimators; trapazoid and ellipse
c  parameters in colli are as follows
c                           trapezoid      ellipse
c       colli(iep,indre,1)    xmin          xmax
c       colli(iep,indre,2)    xmax          zmax
c       colli(iep,indre,3)    zmax1          h
c       colli(iep,indre,4)    zmax2          k
c       colli(iep,indre,5)    zmin1          -
c       colli(iep,indre,6)    zmin2          -
c
c      write(6,100)(colli(iep,indre,j),j=1,6)
c 100  format(1x,'colli(i)=',6(1pg10.2))
      if(colli(iep,indre,1).eq.0.) go to 1
      if(colli(iep,indre,2).eq.0.) go to 1
c
c    test for trapezoid or ellipse
c
      if((colli(iep,indre,5).eq.0.).and.(colli(iep,indre,6).eq.0.))then
c  ellipse
      xmax=colli(iep,indre,1)
      zmax=colli(iep,indre,2)
      h=colli(iep,indre,3)
      k=colli(iep,indre,4)
      ans=((x-h)/xmax)**2+((z-k)/zmax)**2
       if(ans.gt.1.)then
        eject=.true.
c        write(6,110)iv,iep
c		built the diag:
      write(diag(iv),'(2a,i3,a)')'discarded by elliptical collimator'
     s			,' at end-plane #',iep,' of current region'
c 2 next done in send:
c	iepdead(iv)=iep
c	iredead(iv)=indre
       endif
      go to 1
      else
c  trapezoid
      xmin=colli(iep,indre,1)
      xmax=colli(iep,indre,2)
      zmax1=colli(iep,indre,3)
      zmax2=colli(iep,indre,4)
      zmin1=colli(iep,indre,5)
      zmin2=colli(iep,indre,6)
      if(xmax.eq.xmin)then
       write(6,*)' ************colli error xmax=xmin**********'
       go to 1
      endif
c      write(6,*)'  xmax=',xmax, '    xmin=',xmin
c      write(6,*)' zmax1=',zmax1,'   zmax2=',zmax2
c      write(6,*)' zmin1=',zmin1,'   zmin2=',zmin2
      m1=(zmax2-zmax1)/(xmax-xmin)
      b1=zmax1-(m1*xmin)
      m2=(zmin2-zmin1)/(xmax-xmin)
      b2=zmin1-(m2*xmin)
c      write(6,*)' m1=',m1,'  b1=',b1
c      write(6,*)' m2=',m2,'  b2=',b2      
c      write(6,*)' x=',x,'  z=',z
       if((z.lt.(x*m2)+b2).or.(z.gt.(x*m1)+b1)
     &  .or.(x.lt.xmin).or.(x.gt.xmax))then
        eject=.true.
c		built the diag:
      write(diag(iv),'(2a,i3,a)')'discarded by trapezoidal collimator'
     s			,' at end-plane #',iep,' of current region'
c 2 next done in send:
c	iepdead(iv)=iep
c	iredead(iv)=indre
        write(6,110)iv,iep
       endif
       endif
 110   format(1x,'vector #',i4,'   intercepted by collimator #',i3)
c      save /crelfra/,/relfra/,/cdiag/,/ndiag/
  1    return
      end
c
c
      subroutine mirror(x,y,z,cx,cy,cz,x1,y1,z1,x2,y2,z2,imir,imir1)
c*************************************************************************
c                                                                        *
c                       m i r o r                                        *
c                                                                        *
c*************************************************************************
c
	parameter (maxmir=10)
c                next line : begening of the /cer/ package
	common/cer/nmir,rmirin,rmirout,dzmir,dxmir,qmir(3,maxmir)
     s		,tiltz(maxmir),tiltx(maxmir),axp,axm,azp,azm
	 common/mir/tv(3,maxmir),sv(3,maxmir),cv(3,maxmir)
c                previous line : end of the /cer/ package
	dimension av(3),uv(3),bv(3),vv(3),wv(3)
	dimension bvold(3),vvold(3),fvold(3)
c
      av(1)=x
      av(2)=y
      av(3)=z
      uv(1)=cx
      uv(2)=cy
      uv(3)=cz
      xlold=1.e30
      iold=0
      imir=0
      imir1=0
      rmir=rmirin
      imirn=0
      sign=1.
      call sphere(av,uv,xlold,iold,imirn,rmir,bvold,vvold,sign)
      if(iold.gt.0)then
		imir=iold
        do 10 i=1,3
			bv(i)=bvold(i)
 10			vv(i)=vvold(i)
		x1=bv(1)
		y1=bv(2)
		z1=bv(3)
      else
c pas d'intersection du rayon incident:
		x1=x
		y1=y
		z1=z
		x2=x
		y2=y
		z2=z
		return
      endif
c  calcul du rayon reflechi:
      vsu=0.
      do 4 i=1,3
 4		vsu=vsu+vv(i)*uv(i)
      do 5 i=1,3
 5		wv(i)=uv(i)-2.*vsu*vv(i)
c rayon reflechi intercepte par un autre miroir?
      xlold=.5*rmirin
      iold=0
c next temp:	rmir=rmirout
      rmir=rmirin
      imirn=imir
      sign=-1.
      call sphere(bv,wv,xlold,iold,imirn,rmir,fvold,vvold,sign)
      if(iold.gt.0)then
c  point 2 = point d'interception:
		imir1=iold
		x2=fvold(1)
		y2=fvold(2)
		z2=fvold(3)
      else
c  pas d'interception du rayon reflechi:
c  point 2 = rayon reflechi prolonge sur .5*rmirin:
		x2=bv(1)+.5*rmirin*wv(1)
		y2=bv(2)+.5*rmirin*wv(2)
		z2=bv(3)+.5*rmirin*wv(3)
      endif
c      save /cer/,/mir/
      return
      end
c
      subroutine sphere(av,uv,xlold,iold,imirn,rmir,bvold,vvold,sign)
c*************************************************************************
c                                                                        *
c                       s p h e r e                                      *
c                                                                        *
c*************************************************************************
c
c  computes the intersection between a ray and a set of "spherical rectangles"
c  called "mirrors".
c  If there is an intersection, and if it occurs after a flight distance
c  smaller than "xlold", updates "xlold", returns the intersection point "bv",
c  the unit vector "vv" (from "bv" to the center "cv" of the concerned mirror),
c  the number "iold" of this mirror. Dont test mirror#"imirn". Dont keep
c  solutions for which xl<0. If "sign"=+1. choose the solution for which 
c  xl is the greatest (snallest for "sign"=-1). Returns 0 for "iold" if
c  no intersection. If sevral intersection, keep that for which xl is 
c  the smallest.
      parameter (maxmir=10)
c                next line : begening of the /cer/ package
      common/cer/nmir,rmirin,rmirout,dzmir,dxmir,qmir(3,maxmir)
     s		,tiltz(maxmir),tiltx(maxmir),axp,axm,azp,azm
      common/mir/tv(3,maxmir),sv(3,maxmir),cv(3,maxmir)
c                previous line : end of the /cer/ package
      dimension av(3),uv(3),dg(3),bv(3),dv(3),vv(3),ts(3),vn(3),gv(3)
      dimension bvold(3),vvold(3)
c  boucle sur les miroirs:
      do 2 jmir=1,nmir
        if(jmir.eq.imirn)goto2
c  unit vector between "middle" and center:
      ts(1)=tv(2,jmir)*sv(3,jmir)-tv(3,jmir)*sv(2,jmir)
      ts(2)=tv(3,jmir)*sv(1,jmir)-tv(1,jmir)*sv(3,jmir)
      ts(3)=tv(1,jmir)*sv(2,jmir)-tv(2,jmir)*sv(1,jmir)
c  calcul du libre parcours xl:
		xl0=0.
		dg2=0.
      do 1 i=1,3
			dg(i)=cv(i,jmir)-av(i)
			xl0=xl0+dg(i)*uv(i)
 1			dg2=dg2+dg(i)*dg(i)
		delta=xl0*xl0-dg2+rmir*rmir
      if(delta.le.0.)goto2
c  solution "+" = celle a xl grand:
		xl=xl0+sign*sqrt(delta)
      if(xl.lt.0.)goto 2
      if(xl.gt.xlold)goto 2
c dans les limites du miroir?
		det=0.
		des=0.
      do 3 i=1,3
			bv(i)=av(i)+xl*uv(i)
			dv(i)=bv(i)-cv(i,jmir)
 3			vv(i)=dv(i)/rmir
		proj=0.
      do 4 i=1,3
			vn(i)=tv(i,jmir)*cos(axp)-ts(i)*sin(axp)
			gv(i)=qmir(i,jmir)+.5*tv(i,jmir)*dxmir-bv(i)
 4			proj=proj+vn(i)*gv(i)
        if(proj.lt.0.)goto 2
		proj=0.
        do 5 i=1,3
			vn(i)=tv(i,jmir)*cos(axm)-ts(i)*sin(axm)
			gv(i)=qmir(i,jmir)-.5*tv(i,jmir)*dxmir-bv(i)
 5			proj=proj+vn(i)*gv(i)
      if(proj.gt.0.)goto2
		proj=0.
      do 6 i=1,3
			vn(i)=sv(i,jmir)*cos(azp)+ts(i)*sin(azp)
			gv(i)=qmir(i,jmir)+.5*sv(i,jmir)*dzmir-bv(i)
 6			proj=proj+vn(i)*gv(i)
        if(proj.lt.0.)goto2
		proj=0.
      do 7 i=1,3
			vn(i)=sv(i,jmir)*cos(azm)+ts(i)*sin(azm)
			gv(i)=qmir(i,jmir)-.5*sv(i,jmir)*dzmir-bv(i)
 7			proj=proj+vn(i)*gv(i)
      if(proj.gt.0.)goto2
        do 9 i=1,3
			bvold(i)=bv(i)
 9			vvold(i)=vv(i)
		xlold=xl
		iold=jmir
 2      continue
c      save /cer/,/mir/
      return
      end
c
        subroutine rungk(eject,indve)
c*************************************************************************
c                                                                        *
c                       r u n g k                                        *
c                                                                        *
c*************************************************************************
c
c             carry a charged particule through the magnetic field
c                given by the subroutine "field" by integrating
c                the equation of motion by the rung-khuta method.
c                the integration stops on a "end-plane".
c
c             input parameters:
c                pin:input vector(x,y,z,cx,cy,cz,px,py,pz,p)
c					(mm,gev/c)
c                jflag=1:autom.adjust. of the step when field change
c                     =2:step fixed to "hstmin"(mm)
c                prec:threshold value of rel. field change for step adj.
c                hstmax,hstmin:max and min value of step
c                sinc:abs.val. of the relative step change
c                fin:coord. of the end-plane
c                idir=1:the end-plane is orthogonal to x axis
c                    =2: "   "    "    "     "       " y  "
c                    =3: "   "    "    "     "       " z  "
c                epsil:acuracy of the end-plane reachig
c                h:current value of the step length(mm)
c             output parameters:
c                pout:output vector(x,y,z,cx,cy,cz,px,py,pz,p)
c					(mm,gev/c)
c             output flag:
c                ifail=0:o.k.
c                     =1:can't reach the end-plane
c                eject:the subr. field has been called outside its
c                  definition region, or the partic. can't reach the end-plane
c
        include 'snake.inc'
        parameter (tmax=.017)
        parameter (turn=2*pi)
c  electron mass (GeV) and classical radius (mm) for sync.rad.aver.ener.losses
c  computation:
      parameter (xme=.510999e-3,re=2.817940e-12)
      parameter (coef=(2./3.)*re*xme)
c
         logical eject
         logical lpplot
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
c                next line : begening of the /diag/ package
      character diag*100,diagt*100
      common/cdiag/diag(maxve),diagt
      common/ndiag/iepdead(maxve),iredead(maxve)
c                previous line : end of the /diag/ package
c                next line : begening of the /plot/ package
        logical lbox,luser,bare
         logical logx,logy,logz
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
         common/params/pin(10),fin,ctep,step,idir,epsil,pout(10)
     s   ,h,ifail,pal
         dimension f(3)
         dimension secxs(4),secys(4),seczs(4),cold(3),csync(3)
         dimension spinxs(4),spinys(4),spinzs(4)
         dimension hsave(3),hnext(3),hmean(3),hgrad(3)
        dimension xyz(3),xyzt(3),abc(3),xz(9),xx(9),pxyz(3)
        equivalence (x,xyz(1)),(y,xyz(2)),(z,xyz(3)),(xt,xyzt(1))
     @,(yt,xyzt(2)),(zt,xyzt(3)),(a,abc(1)),(b,abc(2)),(c,abc(3)
     @),(xz(1),xyz(1)),(xz(4),abc(1)),(px,pxyz(1)),(py,pxyz(2)),
     @(pz,pxyz(3)),(xz(7),pxyz(1))
      logical start
        parameter (thrd=1./3.)
      data xlight_velocity/.299792458e-3/
      data logx/.false./,logy/.false./,logz/.false./
c	write(6,*)'in rungk'
      xm2=xm*xm
c in case lsynchro false:
      dp=0.
                  call ucopy(pin ,xx,9)
c	write(6,*)'return from ucopy'
c  convention: zero momentum means infinite momentum --> no curvature:
c   zero mom. rays are plotted at each step
      p7=pin(10)
      if(p7.eq.0.)then
		pinv=0.
		gama=1.
      else
c          	this constant for units mm,gev/c,t:
c     (06/11/94: on remplace c=299 792 51 par c=299 792 458)
        	pinv=xlight_velocity/p7
		gama=sqrt((p7*p7+xm2)/xm2)
      endif
*********
* computation of the coefficients for the spin precession:
        if (lspin) then
	  gsur2=x_mag_mon*xm/(twice_spin*xmref)
          beta1 = 1. + gama * ( gsur2 - 1.)
          alpha1 = (gsur2 - beta1)
        endif
*********
        ctmax=cos(tmax)
        tarn=0.
        pal=0.
c tilted end-plane: the tilt axis is idir+1 . so the e-p stands
c         between the plane xx(idir)=0. (for tilt angle=0.)
c             and the plane xx(jdir)=0. (for tilt angle=90.deg)
c               where jdir=idir-1
        jdir=mod(idir+1,3)+1
        vtil=-xx(jdir)*step+xx(idir)*ctep
        do=vtil-fin
        sdo=sign(1.,do)
        iflag=1
        call ucopy(pin(4),csync,3)
        call ucopy(pin(4),cold,3)
c        write(6,*)"xx=",xx
c	write(6,*)'calling field'
        call field(eject,xx,f,logx,logy,logz,.true.)
         if(eject)then
                 call ucopy(pin,xz,9)
                 goto300
         endif
        if(jflag.eq.1) call ucopy(f,hsave,3)
      start=.true.
 10      vtil=-xx(jdir)*step+xx(idir)*ctep
        vpil=-xx(jdir+3)*step+xx(idir+3)*ctep
c  kill the track if injected after the e.p:
      if(start.and.vtil.gt.fin)then
		eject=.true.
		call ucopy(pin,pout,10)
        return
      endif
c plot the first point (to check field continuity):
      if(start)then
          call plorot(xx(1),xx(2),xx(3),xa,ya,za,xb,yb,zb)
          call plorot(f(1),f(2),f(3),xc,yc,zc,bx,by,bz)
          call plorot(px,py,pz,xc,yc,zc,pxa,pya,pza)
          call plotin(indve,xa,ya,za,bx,by,bz,pxa,pya,pza,.false.)
      endif
      start=.false.
        s=h*vpil+vtil
        d2=s-fin
c          snake:don't test the end if y<y end or backw. direct.:
        if(d2.lt.0..or.vpil.lt.0.)goto 30
 20      h=h-d2/vpil
        iflag=iflag+1
 30      continue
        call ucopy(xx,xz,9)
        h2=0.5*h
        h4=0.25*h
        ph=pinv*h
        ph2=0.5*ph
        secxs(1)=(b*f(3)-c*f(2))*ph2*z_particle
        secys(1)=(c*f(1)-a*f(3))*ph2*z_particle
        seczs(1)=(a*f(2)-b*f(1))*ph2*z_particle
c spin deb.
      if(lspin)then
          alpha2 = alpha1 * (a*f(1) + b*f(2) + c*f(3))
          spinxs(1) = (py*(alpha2*c + beta1*f(3))
     s              -  pz*(alpha2*b + beta1*f(2))) * ph2
          spinys(1) = (pz*(alpha2*a + beta1*f(1))
     s              -  px*(alpha2*c + beta1*f(3))) * ph2
          spinzs(1) = (px*(alpha2*b + beta1*f(2))
     s              -  py*(alpha2*a + beta1*f(1))) * ph2
          pxt=px+spinxs(1)
          pyt=py+spinys(1)
          pzt=pz+spinzs(1)
      endif
c spin fin
        xt=x+h2*a+h4*secxs(1)
        yt=y+h2*b+h4*secys(1)
        zt=z+h2*c+h4*seczs(1)
        at=a+secxs(1)
        bt=b+secys(1)
        ct=c+seczs(1)
c        write(6,*)"xyzt=",xyzt
        if(iflag.le.2)call field(eject,xyzt,f,logx,logy,logz,
     @    			.true.)
         if(eject)goto300
        secxs(2)=(bt*f(3)-ct*f(2))*ph2*z_particle
        secys(2)=(ct*f(1)-at*f(3))*ph2*z_particle
        seczs(2)=(at*f(2)-bt*f(1))*ph2*z_particle
c spin deb.
	if(lspin)then
          alpha2 = alpha1 * (at*f(1) + bt*f(2) + ct*f(3))
          spinxs(2) = (pyt*(alpha2*ct + beta1*f(3))
     s              -  pzt*(alpha2*bt + beta1*f(2))) * ph2
          spinys(2) = (pzt*(alpha2*at + beta1*f(1))
     s              -  pxt*(alpha2*ct + beta1*f(3))) * ph2
          spinzs(2) = (pxt*(alpha2*bt + beta1*f(2))
     s              -  pyt*(alpha2*at + beta1*f(1))) * ph2
          pxt=px+spinxs(2)
          pyt=py+spinys(2)
          pzt=pz+spinzs(2)
      endif
c spin fin
        at=a+secxs(2)
        bt=b+secys(2)
        ct=c+seczs(2)
         secxs(3)=(bt*f(3)-ct*f(2))*ph2*z_particle
         secys(3)=(ct*f(1)-at*f(3))*ph2*z_particle
         seczs(3)=(at*f(2)-bt*f(1))*ph2*z_particle
c spin deb.
      if(lspin)then
          alpha2 = alpha1 * (at*f(1) + bt*f(2) + ct*f(3))
          spinxs(3) = (pyt*(alpha2*ct + beta1*f(3))
     s              -  pzt*(alpha2*bt + beta1*f(2))) * ph2
          spinys(3) = (pzt*(alpha2*at + beta1*f(1))
     s              -  pxt*(alpha2*ct + beta1*f(3))) * ph2
          spinzs(3) = (pxt*(alpha2*bt + beta1*f(2))
     s              -  pyt*(alpha2*at + beta1*f(1))) * ph2
          pxt=px+spinxs(3)+spinxs(3)
          pyt=py+spinys(3)+spinys(3)
          pzt=pz+spinzs(3)+spinzs(3)
      endif
c spin fin
         xt=x+h*a+h*secxs(3)
         yt=y+h*b+h*secys(3)
         zt=z+h*c+h*seczs(3)
         at=a+secxs(3)+secxs(3)
         bt=b+secys(3)+secys(3)
         ct=c+seczs(3)+seczs(3)
c         write(6,*)"xyzt=",xyzt
         if(jflag.le.2) call field(eject,xyzt,f,logx,logy,logz,
     @				.true.)
         if(eject)goto300
         if(iflag.eq.1) call ucopy(f,hnext,3)
         secxs(4)=(bt*f(3)-ct*f(2))*ph2*z_particle
         secys(4)=(ct*f(1)-at*f(3))*ph2*z_particle
         seczs(4)=(at*f(2)-bt*f(1))*ph2*z_particle
c spin deb.
      if(lspin)then
          alpha2 = alpha1 * (at*f(1) + bt*f(2) + ct*f(3))
          spinxs(4) = (pyt*(alpha2*ct + beta1*f(3))
     s              -  pzt*(alpha2*bt + beta1*f(2))) * ph2
          spinys(4) = (pzt*(alpha2*at + beta1*f(1))
     s              -  pxt*(alpha2*ct + beta1*f(3))) * ph2
          spinzs(4) = (pxt*(alpha2*bt + beta1*f(2))
     s              -  pyt*(alpha2*at + beta1*f(1))) * ph2
          px=px+(spinxs(1)+spinxs(4)+2.0*(spinxs(2)+spinxs(3)))*thrd
          py=py+(spinys(1)+spinys(4)+2.0*(spinys(2)+spinys(3)))*thrd
          pz=pz+(spinzs(1)+spinzs(4)+2.0*(spinzs(2)+spinzs(3)))*thrd
      endif
c spin fin
         z=z+(c+(seczs(1)+seczs(2)+seczs(3))*thrd)*h
         y=y+(b+(secys(1)+secys(2)+secys(3))*thrd)*h
         x=x+(a+(secxs(1)+secxs(2)+secxs(3))*thrd)*h
         a=a+(secxs(1)+secxs(4)+2.0*(secxs(2)+secxs(3)))*thrd
         b=b+(secys(1)+secys(4)+2.0*(secys(2)+secys(3)))*thrd
         c=c+(seczs(1)+seczs(4)+2.0*(seczs(2)+seczs(3)))*thrd
        pal=pal+h
c    plot if there is enough change in direction:
        ctot=0.
        lpplot=.false.
      du2=0.
      do 1 i=1,3
	  du=abc(i)-csync(i)
	  du2=du2+du*du
 1         ctot=ctot+cold(i)*abc(i)
c sync.rad.aver.ener.loss:
      if(lsynchro.and.p7.ne.0.)then
		p2=p7*p7
		p3=p2*p7
		dp=coef*du2*p3*gama/(h*xm*xm2)
		p7=p7-dp
c          	this constant for units mm,gev/c,t:
        	pinv=xlight_velocity/p7
		gama=sqrt((p7*p7+xm2)/xm2)
	        call ucopy(abc,csync,3)
      endif
        if(ctot.lt.ctmax.or.p7.eq.0.)then
c                periodic check-up:
          tarn=tarn+acos(ctot)
          if(tarn.ge.turn)then
c           this trajectory is probably trapped in the field:
            eject=.true.
c		built the diag:
		write(diag(indve),'(3(a,e10.3))')
     s		 'total rotation angle exceeds max.:'
     s           ,tarn, ' >',turn,' rad.'
            goto300
          endif
          call ucopy(abc,cold,3)
c         plot this point if room:
          lpplot=.true.
          call plorot(x,y,z,xa,ya,za,xb,yb,zb)
          call plorot(f(1),f(2),f(3),xc,yc,zc,bx,by,bz)
          call plorot(px,py,pz,xc,yc,zc,pxa,pya,pza)
	  call plotin(indve,xa,ya,za,bx,by,bz,pxa,pya,pza,.false.)
        endif
c       write(6,200)x,y,z,h,f,iflag
200      format(8f10.4)
         if(iflag.gt.1)goto90
c snake: if(sdo.eq.sign(1.,abc(idir))) goto 120
         if(jflag.eq.1)  goto40
c
c                fixed step:
c
          call ucopy(xz,xx,9)
c                to next step:
       goto 10
c
c                variable step:
c
40       do 50 i=1,3
         hmean(i)=(hnext(i)+hsave(i))*0.5
50       hgrad(i)=hnext(i)-hsave(i)
          ht=vmodul(hmean)
         hratio=0
         sh=sign(1.,h)
c            if the field is smaller than 10**-2 tesla (100 gauss)...
c       ...the step can increase or stay constant,else it can ....
c       ...also decrease
         if(ht.le.1.e-2) goto 60
         hg=vmodul(hgrad)
         hratio=hg/ht
         sh=sign(1.,h)
         if(abs(h).le.hstmin) goto 60
         if(hratio.gt.prec) goto80
 60       call ucopy(xz,xx,9)
c                step const. or increase?
         if(hratio.ge.0.1*prec) goto 70
c                step increase:
         h=h+sinc*sh
         h=amin1(abs(h),hstmax)*sh
 70       call ucopy(hnext,hsave,3)
         goto 10
c                step decrease and go back to previouv step:
 80       pal=pal-h
      if(p7.ne.0.)p7=p7+dp
         if(lpplot.and.(indve+nplvold).le.maxplv)then
            if(nplp(indve+nplvold).lt.maxplp)
     s      nplp(indve+nplvold)=nplp(indve+nplvold)-1
	 endif
         h=h-sinc*sh
         h=amax1(abs(h),hstmin)*sh
         call ucopy(hsave,f,3)
c                to next step
         goto 10
c                acurate landing:
90      vtil=-xz(jdir)*step+xz(idir)*ctep
         d1=vtil-fin
         hold=abs(h)
         h=-d1/vpil
c                one can finish by linear approx.?
         if(abs(h).lt.epsil)goto100
c                is there a problem?
         if(abs(h).gt.hold)goto120
         call ucopy(xz,xx,9)
         iflag=iflag+1
         goto 30
c                linear approx.:all is o.k
 100      ifail=0
         do 110 mm=1,3
          pout(mm)=xyz(mm)+abc(mm)*h
         pout(mm+3)=abc(mm)
 110      pout(mm+6)=pxyz(mm)
         pout(10)=p7
         pal=pal+h
          if(indve.le.maxplv)then
            call plorot(f(1),f(2),f(3),xa,ya,za,bx,by,bz)
            bplot(1,indve)=bx
            bplot(2,indve)=by
            bplot(3,indve)=bz
          endif
         goto 130
c                there is a problem about landing:
 120      ifail=2
c		built the diag:
        write(diag(indve),'(a)')
     s		 'unsuccessfull landing on the end plane'
          eject=.true.
         goto 300
130      return
300      continue
         do 301 mm=1,9
301              pout(mm)=xz(mm)
         pout(10)=p7
         return
         end
          function vmodul(a)
          dimension a(3)
          vmodul=sqrt(a(1)**2+a(2)**2+a(3)**2)
c      save /alter/,/logic/,/particle/,/spin/,/cdiag/,/ndiag/,/cplot/,
c     &     /plot/,/params/
          return
          end
          subroutine ucopy(a,b,n)
          dimension a(*),b(*)
          if(n.eq.0) return
           do 2 ii=1,n
2         b(ii)=a(ii)
          return
          end
        subroutine circle(xc,yc,zc,r,nuc)
c*************************************************************************
c                                                                        *
c                       c i r c l e                                      *
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
     @		    ,iplpi(maxplv),iplpf(maxplv)
     @    ,plotus(nuseumax,nulumax,3),nuseu(nulumax),nulu
     @		     ,luser,nclic,lclic(maxre),clight(4,maxre),
     @    epp(2,3,maxep,3) ,neps(maxep),regp(2,3,maxre,3),nrs(maxre)
c                previous line : end of the /plot/ package
      if(nulu.eq.nulumax)then
		write(6,*)' warning circle: too many user line'
      return
      endif
      nulu=nulu+1
      dt=2.*pi/real(nuseumax-1)
      t=0.
      do 1 i=1,nuseumax
		t=t+dt
		x=xc+r*cos(t)
		y=yc+r*sin(t)
		z=zc
                call plorot(x,y,z,xa,ya,za,xb,yb,zb)
		plotus(i,nuc,1)=xa
 		plotus(i,nuc,2)=ya
 1		plotus(i,nuc,3)=za
      nuseu(nuc)=nuseumax
c      save /cplot/,/plot/
      return
      end
        subroutine analyf(eject,xyz,f,index,indve)
c*************************************************************************
c                                                                        *
c                       a n a l y f                                      *
c                                                                        *
c*************************************************************************
c
c                find the 3 components of the field f(3) at the point xyz(3)
c         analyf is called by field in the case of analytic field definition
c         .true.output value for eject causes the death of the particule.
c
c       input:
c             xyz(x y or z):relative coord. of the point
c             index:field index of the region
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
        common/aject/ieject
c   for meth=3 , type=2 , index=1 (dipole with grad. ):
        common/dipgra/rm,alp,bet
c   for meth=3 , type=2 , index=16+17 (arbitrary field numeric.deriv.)
        common/dipff/xtra,ytra,atra,dref,rms,alps,bets
     @           ,xc(2),yc(2),r0(2),e0(2),s0(2),s1(2),s2(2),s3(2)
     @           ,tbound,dst
c   for meth=3 , type=2 , index=2 to 14 (raytrace library ):
        common/raytra/fact1(5),qrad1,fact2(5),qrad2,dl
c   for meth=3 , type=2 , index=15 (clam):
        common/clam1/d0,dx,dy,bdref,ta
c   
        dimension f(3),xyz(3)
         logical eject
      real x,y,z
        x=xyz(1)
        y=xyz(2)
        z=xyz(3)
c        write(6,*)"xyz=",xyz
        goto(101,102,103,104,105,106,107,108,109,
     s       110,111,112,113,114,115,116,117,118,119,120),index
c
c                index=1:field with gradient (alpha and beta)
c                        (x=0.,y=0.,z=0.)=center of symmetry
c                         z=0.:symmetry plane
c
c       rectangular-->cylindrical coord.:
 101     r=sqrt(x*x+y*y)
        rn=r/rm
        te=atan2(y,x)
        zn=z/rm
        dr=rn-1.
c       z component of the field in sym.pl.and its derivatives/rn:
        b0=1.-alp*dr+bet*dr*dr
        db=  -alp+2.*bet*dr
        d2b=     +2.*bet
c       extrapolation from the symmetry plane:
        cr=zn*db
        f(3)=b0-.5*zn*zn*(d2b+db/rn)
c       cylindrical--->rectangular field components:
        f(1)=cr*cos(te)
        f(2)=cr*sin(te)
c      save /cdiag/,/ndiag/,/aject/,/dipgra/,/dipff/,/raytra/,/clam1/
        return
c
c         index=2:unit gradient (1.t/mm ) along x and z (qpole axis = y )
c
 102     f(1)=z
        f(2)=0.
        f(3)=x
c      save /cdiag/,/ndiag/,/aject/,/dipgra/,/dipff/,/raytra/,/clam1/
        return
c
c               index=3:quadrupole
c
 103	continue
        call quad(xyz(1),xyz(3),fact1(1),qrad1,f(1),f(3))
        f(2)=0.
c      save /cdiag/,/ndiag/,/aject/,/dipgra/,/dipff/,/raytra/,/clam1/
        return
c
c              index=4:dipole entrance fringe
c
 104     ydum=-xyz(2)
        xdum=xyz(1)
        zdum=xyz(3)
        call dfringe(xdum,ydum,zdum,fact1(1),qrad1,f(1),f(2),f(3))
        f(2)=-f(2)
c      save /cdiag/,/ndiag/,/aject/,/dipgra/,/dipff/,/raytra/,/clam1/
        return
c
c              index=5:dipole exit fringe
c
 105     ydum=xyz(2)
        xdum=xyz(1)
        zdum=xyz(3)
        call dfringe(xdum,ydum,zdum,fact1(1),qrad1,f(1),f(2),f(3))
c      save /cdiag/,/ndiag/,/aject/,/dipgra/,/dipff/,/raytra/,/clam1/
        return
c
c               index=6:uniform dipole
c
 106     f(1)=fact1(1)
        f(2)=fact1(2)
        f(3)=fact1(3)
        return
c
c           index=7:entrance multipole
c
 107     xdum=xyz(1)
        ydum=-xyz(2)
        zdum=xyz(3)
        call mpoles(1,xdum,zdum,ydum,fact1,qrad1,f(1),f(3),f(2))
        f(2)=-f(2)
c      save /cdiag/,/ndiag/,/aject/,/dipgra/,/dipff/,/raytra/,/clam1/
        return
c
c           index=8:uniform field multipole
c
 108     xdum=xyz(1)
        ydum=xyz(2)
        zdum=xyz(3)
        call mpoles(2,xdum,zdum,ydum,fact1,qrad1,f(1),f(3),f(2))
        return
c
c           index=9:exit multipole
c
 109     xdum=xyz(1)
        ydum=xyz(2)
        zdum=xyz(3)
c      write(6,*)'calling mpoles, exit'
        call mpoles(3,xdum,zdum,ydum,fact1,qrad1,f(1),f(3),f(2))
c      write(6,*)'exit mpoles'
c      save /cdiag/,/ndiag/,/aject/,/dipgra/,/dipff/,/raytra/,/clam1/
        return
c
c           index=10 qq overlapping fringe fields
c
 110     xdum=xyz(1)
        ydum=xyz(2)
        zdum=xyz(3)
        call qqfringe(xdum,ydum,zdum,dl,f(1),f(2),f(3),fact1,fact2,
     s               qrad1,qrad2)
        return
c
c           index=11 qd overlapping fringe fields
c
 111     xdum=xyz(1)
        ydum=xyz(2)
        zdum=xyz(3)
c        write(9,936)xdum,ydum,zdum,dl,(fact1(i),i=1,5),
        call qdfringe(xdum,ydum,zdum,dl,f(1),f(2),f(3),fact1,fact2,
     &                qrad1,qrad2)
c      save /cdiag/,/ndiag/,/aject/,/dipgra/,/dipff/,/raytra/,/clam1/
        return
c
c           index=12 dd overlapping fringe fields
c
 112     xdum=xyz(1)
        ydum=xyz(2)
        zdum=xyz(3)
        call ddfringe(xdum,ydum,zdum,dl,f(1),f(2),f(3),fact1,fact2,
     &                qrad1,qrad2)
        return
c
c           index=13 dq overlapping fringe fields (old version,PV 4/9/96)
c
 113     xdum=xyz(1)
        ydum=xyz(2)
        zdum=xyz(3)
        call dqfringe(xdum,ydum,zdum,dl,f(1),f(2),f(3),fact1,fact2,
     &                qrad1,qrad2)
c      save /cdiag/,/ndiag/,/aject/,/dipgra/,/dipff/,/raytra/,/clam1/
        return
c
c           index=14 dq overlapping fringe fields (new version,PV 4/9/96)
c
 114     xdum=xyz(1)
        ydum=xyz(2)
        zdum=xyz(3)
        call dqfringe2(xdum,ydum,zdum,dl,f(1),f(2),f(3),fact1,fact2,
     &                qrad1,qrad2)
        return
c
c                index=15:central field of a clam:
c
 115     d1=abs(d0+x*dx+y*dy)
        d2=d1*d1
        eject=d1.le.0..or.abs(z).gt.(d1*ta)
      if(eject)then
c		built the diag:
        write(diagt,'(a)')
     s		'try to go in a forbidden region of a clam'
c      save /cdiag/,/ndiag/,/aject/,/dipgra/,/dipff/,/raytra/,/clam1/
		return
      endif
        ccf=bdref/(d2+z*z)
        f(1)=-ccf*z*dx
        f(2)=-ccf*z*dy
        f(3)=ccf*d1
c      save /cdiag/,/ndiag/,/aject/,/dipgra/,/dipff/,/raytra/,/clam1/
        return
c
c                index=16:dipole with gradient numerical deriv:
c
 116     call dipole(x,y,b00)
        call dipole(x+dst,y,bp0)
        call dipole(x-dst,y,bm0)
        call dipole(x,y+dst,b0p)
        call dipole(x,y-dst,b0m)
        div1=1./(2.*dst)
        div2=1./(dst*dst)
        grx1=(bp0-bm0)*div1
        gry1=(b0p-b0m)*div1
        b002=2.*b00
        grx2=bp0+bm0-b002
        gry2=b0p+b0m-b002
        f(1)=z*grx1
        f(2)=z*gry1
        f(3)=b00-.5*z*z*(grx2+gry2)*div2
c      save /cdiag/,/ndiag/,/aject/,/dipgra/,/dipff/,/raytra/,/clam1/
        return
c
c                index=17:raytrace style dipole,old version:
c
  117	xdum=xyz(1)/10.
        ydum=xyz(2)/10.
        zdum=xyz(3)/10.
c        write(7,*)xdum,ydum,zdum
        call dipolert(xdum,zdum,ydum,f(1),f(3),f(2))
c        write(7,*)xdum,ydum,zdum,(f(i),i=1,3)
c  get the right sign on the field
       do 12 i=1,3
  12   f(i)=-f(i)
c      save /cdiag/,/ndiag/,/aject/,/dipgra/,/dipff/,/raytra/,/clam1/
        return
c
c    		index=18:raytrace style dipole, new version
c
  118   xdum=xyz(1)/10.
        ydum=xyz(2)/10.
        zdum=xyz(3)/10.
c        write(7,*)xdum,ydum,zdum
         call dipolert2(xdum,zdum,ydum,f(1),f(3),f(2))
c        write(7,*)xdum,ydum,zdum,(f(i),i=1,3)
c  get the right sign on the field
       do 13 i=1,3
  13   f(i)=-f(i)
        return
c
c              index=19:polarized target
c
 119    ydum=xyz(2)/1000.
        xdum=xyz(1)/1000.
        zdum=xyz(3)/1000.
        call hlmhltz(xdum,ydum,zdum,qrad1,fact1(1),f(1),f(2),f(3))
c      save /cdiag/,/ndiag/,/aject/,/dipgra/,/dipff/,/raytra/,/clam1/
        return
c
c        index=20: Solenoid
c
 120    ydum=xyz(2)/1000.
        xdum=xyz(1)/1000.
        zdum=xyz(3)/1000.
c        write(6,*)"Calling solenoid",xyz(1)
        call solenoid(xdum,ydum,zdum,fact1(1),fact1(2),fact1(3),
     &  bx,by,bz)
c        write(6,*)"return from solenoid",bx,by,bz
        f(1)=bx
        f(2)=by
        f(3)=bz   
c      save /cdiag/,/ndiag/,/aject/,/dipgra/,/dipff/,/raytra/,/clam1/
      return
        end
c
      subroutine mpoles(in,nx,ny,nz,ngrad,nrad,nbx,nby,nbz)
c stolen from raytrace jjl 10/8/86
c****                                                                   
c**** calculation of multipole(poles) field components                  
c****                                                                   
c****                                                                   
c****                                                                   
c**** 2 - quadrupole  (grad1)                                           
c**** 3 - hexapole    (grad2)                                           
c**** 4 - octapole    (grad3)                                           
c**** 5 - decapole    (grad4)                                           
c**** 6 - dodecapole  (grad5)                                           
c****                                                                   
c****                                                                   
      implicit real*8(a-h,o-z)                                          
      real nx,ny,nz,ngrad(5),nrad,nbx,nby,nbz
      common  /blck91/  c0, c1, c2, c3, c4, c5                          
c      data c0, c1, c2, c3, c4, c5/0.1122d0,8.5d0,-1.4982d0,3.5882d0,
c     &     -2.1209d0,1.7230d0/
c  above are not quite right -JJL 9/4/98
      data c0, c1, c2, c3, c4, c5/0.1039d0,6.27108d0,-1.51247d0,
     &     3.59946d0,-2.1323d0,1.7230d0/
      if (nrad.eq.0.)then
      write(6,*)' error in mpoles,  nrad= 0.'
      call exit(0)
      endif
      grad1 = -ngrad(1)/nrad
      grad2 =  ngrad(2)/nrad**2
      grad3 = -ngrad(3)/nrad**3
      grad4 =  ngrad(4)/nrad**4
      grad5 = -ngrad(5)/nrad**5
      rad=nrad
      d = 2. * rad                                                      
      frh  = 1.d0
      fro  = 1.d0
      frd  = 1.d0
      frdd = 1.d0
      dh  = frh *d
      do  = fro *d
      dd  = frd *d
      ddd = frdd*d
      x = nx
      y = ny
      z = nz
      x2 = x*x                                                          
      x3 = x2*x                                                         
      x4 = x3*x                                                         
      x5 = x4*x                                                         
      x6 = x5*x
      x7 = x6*x
      y2 = y*y                                                          
      y3 = y2*y                                                         
      y4 = y3*y                                                         
      y5 = y4*y                                                         
      y6 = y5*y
      y7 = y6*y
      go to ( 2, 1, 2 ) , in                                            
      print 3, in                                                       
    3 format( '  error in bpoles in= ', i5 ///)                          
      call exit(0)   
    1 continue                                                          
      b2x = grad1*y                                                     
      b2y = grad1*x                                                     
      b3x = grad2*2.*x*y                                                
      b3y = grad2*(x2-y2)                                               
      b4x = grad3*(3.*x2*y-y3)                                          
      b4y = grad3*(x3-3.*x*y2)                                          
      b5x = grad4*4.*(x3*y-x*y3)                                        
      b5y = grad4*(x4-6.*x2*y2+y4)                                      
      b6x = grad5*(5.*x4*y-10.*x2*y3+y5)                                
      b6y = grad5*(x5-10.*x3*y2+5.*x*y4)                                
      bx = b2x + b3x + b4x + b5x + b6x                                  
      by = b2y + b3y + b4y + b5y + b6y                                  
      bz = 0.                                                           
      bt =   sqrt( bx*bx + by*by )                                     
      nbx=bx
      nby=by
      nbz=bz
      return                                                            
c****
c****
c****
    2 s = z/d                                                           
      call bpls( 2, d, s, re, g1, g2, g3, g4, g5, g6 )
      b2x = grad1*( re*y - (g2/12.)*(3.*x2*y + y3) +                    
     1   (g4/384.)*(5.*x4*y + 6.*x2*y3 + y5 ) -                         
     2   (g6/23040.)*(7.*x6*y + 15.*x4*y3 + 9.*x2*y5 + y7)  )
      b2y = grad1*( re*x - (g2/12.)*(x3 + 3.*x*y2) +                    
     1   (g4/384.)*(x5 + 6.*x3*y2 + 5.*x*y4 ) -                         
     2   (g6/23040.)*(x7 + 9.*x5*y2 + 15.*x3*y4 + 7.*x*y6) )
      b2z = grad1*( g1*x*y - (g3/12.)*(x3*y + x*y3 ) +                  
     1   (g5/384.)*(x5*y +2.*x3*y3 + x*y5)  )
c****
c****
      ss = z/dh  + dsh
      call bpls( 3, dh, ss, re, g1, g2, g3, g4, g5, g6 )
      b3x = grad2*( re*2.*x*y - (g2/48.)*(12.*x3*y + 4.*x*y3 ) )        
      b3y = grad2*( re*(x2-y2) - (g2/48.)*(3.*x4 + 6.*x2*y2 - 5.*y4 ) ) 
      b3z = grad2*( g1*(x2*y - y3/3.) - (g3/48.)*(3.*x4*y+2.*x2*y3-y5)) 
c****
c****
      ss = z/do  + dso
      call bpls( 4, do, ss, re, g1, g2, g3, g4, g5, g6 )
      b4x = grad3*( re*(3.*x2*y - y3) - (g4/80.)*(20.*x4*y - 4.*y5 ) )  
      b4y = grad3*( re*(x3 - 3.*x*y2) - (g4/80.)*(4.*x5-20.*x*y4 ) )    
      b4z = grad3*g1*(x3*y - x*y3 )                                     
c****
c****
      ss = z/dd  + dsd
      call bpls( 5, dd, ss, re, g1, g2, g3, g4, g5, g6 )
      b5x = grad4*re*(4.*x3*y - 4.*x*y3)                                
      b5y = grad4*re*(x4 - 6.*x2*y2 + y4 )                              
      b5z = grad4*g1*(x4*y - 2.*x2*y3 + y5/5. )                         
c****
c****
      ss = z/ddd + dsdd
      call bpls( 6, ddd,ss, re, g1, g2, g3, g4, g5, g6 )
      b6x = grad5*re*(5.*x4*y - 10.*x2*y3 + y5 )                        
      b6y = grad5*re*(x5 - 10.*x3*y2 + 5.*x*y4 )                        
      b6z = 0.                                                          
c****
c****
      bx = b2x + b3x + b4x + b5x + b6x                                  
      by = b2y + b3y + b4y + b5y + b6y                                  
      bz = b2z + b3z + b4z + b5z + b6z                                  
      bt =   sqrt( bx*bx + by*by + bz*bz )                             
      nbx=bx
      nby=by
      nbz=bz
c      save /blck91/
      return                                                            
      end                                                               
      subroutine bpls ( igp, d, s, re, g1, g2, g3, g4, g5, g6 )
c****
c****
c****
      implicit real*8 (a-h,o-z)
c****
c****
      common  /blck91/  c0, c1, c2, c3, c4, c5                          
c****
c****
      s2 = s*s
      s3 = s2*s
      s4 = s2*s2
      s5 = s4*s
      cs = c0 + c1*s + c2*s2 + c3*s3 + c4*s4 + c5*s5                    
      cp1 =(c1 + 2.*c2*s + 3.*c3*s2 + 4.*c4*s3 + 5.*c5*s4) / d          
      cp2 = (2.*c2 + 6.*c3*s + 12.*c4*s2 + 20.*c5*s3  ) / (d*d)         
      cp3 = ( 6.*c3 + 24.*c4*s + 60.*c5*s2 ) / (d**3)                   
      cp4 = ( 24.*c4 + 120.*c5*s ) / (d**4)                             
c****
      cp5 = 120.*c5/(d**5)
c****
c****
c****
      if( abs(cs) .gt. 70. )  cs = sign(70.d0, cs )                   
      e = exp(cs)                                                      
      re = 1./(1. + e)                                                  
      ere = e*re                                                        
      ere1= ere*re
      ere2= ere*ere1                                                    
      ere3= ere*ere2                                                    
      ere4= ere*ere3                                                    
c****
      ere5= ere*ere4
      ere6= ere*ere5
c****
c****
      cp12 = cp1*cp1                                                    
      cp13 = cp1*cp12                                                   
      cp14 = cp12*cp12                                                  
      cp22 = cp2*cp2                                                    
c****
      cp15 = cp12*cp13
      cp16 = cp13*cp13
      cp23 = cp2*cp22
      cp32 = cp3*cp3
c****
c****
      if( igp .eq. 6 ) return
      g1 = -cp1*ere1                                                    
c****
c****
      if( igp .eq. 5 ) return
      if( igp .eq. 4 ) go to 1
      g2 =-( cp2+cp12   )*ere1    + 2.*cp12 * ere2                      
      g3 =-(cp3 + 3.*cp1*cp2 + cp13  ) * ere1      +                    
     1   6.*(cp1*cp2 + cp13)*ere2 - 6.*cp13*ere3                        
c****
c****
      if( igp .eq. 3 ) return
1     g4 = -(cp4 + 4.*cp1*cp3 + 3.*cp22 + 6.*cp12*cp2 + cp14)*ere1  +   
     1   (8.*cp1*cp3 + 36.*cp12*cp2 + 6.*cp22 + 14.*cp14)*ere2    -     
     2   36.*(cp12*cp2 + cp14)*ere3       + 24.*cp14*ere4               
c****
c****
      if( igp .ne. 2 ) return
      g5 = (-cp5 - 5.*cp1*cp14 - 10.*cp2*cp3 - 10.*cp12*cp3 -
     1     15.*cp1*cp22 - 10.*cp13*cp2 - cp15)*ere1 +
     2     (10.*cp1*cp4 +20.*cp2*cp3 +60.*cp12*cp3 + 90.*cp1*cp22 +
     3     140.*cp13*cp2 +30.*cp15)*ere2 + (-60.*cp12*cp3 -
     4     90.*cp1*cp22 - 360.*cp13*cp2 - 150.*cp15)*ere3 +
     5     (240.*cp13*cp2 +240.*cp15)*ere4 + (-120.*cp15)*ere5
      g6 = (-6.*cp1*cp5 - 15.*cp2*cp4 - 15.*cp12*cp4 - 10.*cp32 -
     1     60.*cp1*cp2*cp3 - 20.*cp13*cp3 - 15.*cp23 - 45.*cp12*cp22 -
     2     15.*cp14*cp2 - cp16)*ere1 + (12.*cp1*cp5 + 30.*cp2*cp4 +
     3     90.*cp12*cp4 +20.*cp32 + 360.*cp1*cp2*cp3 +280.*cp13*cp3 +
     4     90.*cp23 + 630.*cp12*cp22 + 450.*cp14*cp2 + 62.*cp16)*ere2 +
     5     (-90.*cp12*cp4 - 360.*cp1*cp2*cp3 -720.*cp13*cp3 -90.*cp23 -
     6     1620.*cp12*cp22 -2250.*cp14*cp2 - 540.*cp16)*ere3 +
     7     (480.*cp13*cp3 + 1080.*cp12*cp22 + 3600.*cp14*cp2 +
     8     1560.*cp16)*ere4 + (-1800.*cp14*cp2 - 1800.*cp16)*ere5 +
     9     720.*cp16*ere6
c****
c      save /blck91/
      return
      end
c
c
      subroutine quad(x,y,fact,a,fx,fy)
      real x,y,fx,fy,fact,a
c
c    calculates quadrupole field components fx and fy for given
c    input coordinates x,y     jjl 10/2/86
c
c          fact = max field strength at r=a
c             a = aperture radius 
      if(x**2+y**2.gt.a**2)go to 900
      fx=-y*fact/a
      fy=-x*fact/a
      return
  900 write(6,*)'point outside quadrupole aperture'
      return
      end
      subroutine dipolert (x,y,z,bbx,bby,bbz)
c*************************************************************************
c                                                                        *
c                       d i p o l e r t                                  *
c                                                                        *
c*************************************************************************
c
c  dipolert- analytic calculation of dipole magnetic fields
c            stolen from raytrace adapted for use in snake
c                                         -jjl  8/21/89
c    data(i) are read in snake in ifield
c    dipolert is called from snake in analyf
c                snake provides x,y,z in the c-axis system
c                dipolert returns bx,by,bz in the same system
c
      implicit real*8(a-h,o-z)
      real x,y,z,bbx,bby,bbz,data(75),pii
      real*8  ndx, k
      common  /diprt/data
      common  /block10old/  bx, by, bz, k, tc, dtc
      common  /blck20/  ndx,bet1,gama,delt,csc
      common  /blck21/  rca,dels,br,s2,s3,s4,s5,s6,s7,s8,scor
      common  /blck22/  d, dg, s, bf, bt
      common  /blck23/  c0, c1, c2, c3, c4, c5
      common  /blck24/  rb, xc, zc
      common  /blck25/  in, mtyp
      dimension tc(6), dtc(6)
C pii = pi/180
	  data pii/0.0174532925/
      lf1  = dble(data(  1  ))
      lu1  = dble(data(  2  ))
      lf2  = dble(data(  3  ))
      dg   = dble(data(  4  ))
      mtyp = dble(data(  5  ))
      a    = dble(data( 11  ))
      b    = dble(data( 12  ))
      d    = dble(data( 13  ))
      rb   = dble(data( 14  ))
      bf   = dble(data( 15  ))
      phi  = dble(data( 16  ))
      alpha= dble(data( 17  ))
      beta = dble(data( 18  ))
      ndx  = dble(data( 19  ))
      bet1 = dble(data( 20  ))
      gama = dble(data( 21  ))
      delt = dble(data( 22  ))
      z11  = dble(data( 25  ))
      z12  = dble(data( 26  ))
      z21  = dble(data( 27  ))
      z22  = dble(data( 28  ))
      br1  = dble(data( 41  ))
      br2  = dble(data( 42  ))
      xcr1 = dble(data( 43  ))
      xcr2 = dble(data( 44  ))
c      write(8,*)' in dipolert 1',x,y,z
c**** transform from c-axis to b-axis
c****
      xpric=x+(2.*rb*sin(pii*phi/2.)*sin(pii*((phi/2.)-beta)))
      zpric=z+(2.*rb*sin(pii*phi/2.)*cos(pii*((phi/2.)-beta)))
      cosa=cos(pii*(phi-alpha-beta))
      sina=sin(pii*(phi-alpha-beta))
      tc(1)=(-xpric*cosa)+(zpric*sina)
      tc(2)=dble(y)
      tc(3)=(-xpric*sina)-(zpric*cosa)
c      write(8,*)' in dipolert 1',(tc(i),i=1,3)
c
c  test for region
c
      if(tc(3).gt.z11)then
       bx=0.
       by=0.
       bz=0.
       go to 10
      endif
      if(tc(3).ge.z12) go to 1
      if(tc(3).lt.z12) go to 2
c      go to 15
c****
c**** in designates magnet regions for bfun
c****
  1   in = 1
      xc= rb*cos( alpha/ 57.29578 )
      zc=-rb*sin( alpha/ 57.29578 )
c****
      c0   = dble(data( 29  ))
      c1   = dble(data( 30  ))
      c2   = dble(data( 31  ))
      c3   = dble(data( 32  ))
      c4   = dble(data( 33  ))
      c5   = dble(data( 34  ))
      dels = dble(data( 45  ))
      rca  = dble(data( 47  ))
      csc = cos( alpha/57.29578 )
      scor = dble(data(49 ))
      s2   = dble(data( 51  )) / rb    + rca/2.d0
      s3   = dble(data( 52  )) / rb**2
      s4   = dble(data( 53  )) / rb**3 + rca**3/8.d0
      s5   = dble(data( 54  )) / rb**4
      s6   = dble(data( 55  )) / rb**5 + rca**5/16.d0
      s7   = dble(data( 56  )) / rb**6
      s8   = dble(data( 57  )) / rb**7 + rca**7/25.6d0
c
      call ndipold
c      write(8,*)' return from ndipold',bx,by,bz
c**** transform to c-axis system
c***
      bxdum=bx
      bzdum=bz
      bx=(-bxdum*cosa)-(bzdum*sina)
      bz=(+bxdum*sina)-(bzdum*cosa)
      go to 10
c**** transform to second vfb coord system
c***
  2   copab =cos( (phi-alpha-beta)/57.29578)
      sipab =sin( (phi-alpha-beta)/57.29578)
      cospb =cos( (phi/2.-beta)/57.29578 )
      sinpb =sin( (phi/2.-beta)/57.29578 )
      sip2 =sin( (phi/2.)/57.29578 )
      xt = tc(1)
      zt = tc(3)
      vxt = tc(4)
      vzt = tc(6)
      tc(3) = - zt  *copab +  xt  *sipab -2.*rb*sip2*cospb
      tc(1) = - zt  *sipab -  xt  *copab -2.*rb*sip2*sinpb
      tc(6) = - vzt *copab +  vxt *sipab
      tc(4) = - vzt *sipab -  vxt *copab
      if(tc(3).lt.z21)go to 3
      if(tc(3).gt.z22)then
       bx=0.
       by=0.
       bz=0.
       go to 10
      endif
      if(tc(3).le.z22)go to 4
c****
c****
c**** uniform field integration region
c****
c****
  3   in = 2
      xc=-rb*cos( beta / 57.29578 )
      zc=-rb*sin( beta / 57.29578 )
c
      call ndipold
c      btot=btot+(by*rho/ 57.29587)
c      if(theta.eq.22.5d+00)br0=by
c      write(1,*)' ',theta,by
      go to 10
c***
c***
c**** setup for second fringe field and integration
c****
c****
  4   br   = br2
      c0   = dble(data( 35  ))
      c1   = dble(data( 36  ))
      c2   = dble(data( 37  ))
      c3   = dble(data( 38  ))
      c4   = dble(data( 39  ))
      c5   = dble(data( 40  ))
      dels = dble(data( 46  ))
      rca  = dble(data( 48  ))
      scor = dble(data(50 ))
      csc = cos( beta /57.29578 )
      s2   = dble(data( 58  )) / rb    + rca/2.d0
      s3   = dble(data( 59  )) / rb**2
      s4   = dble(data( 60  )) / rb**3 + rca**3/8.d0
      s5   = dble(data( 61  )) / rb**4
      s6   = dble(data( 62  )) / rb**5 + rca**5/16.d0
      s7   = dble(data( 63  )) / rb**6
      s8   = dble(data( 64  )) / rb**7 + rca**7/25.6d0
      in = 3
      xc=-rb*cos( beta / 57.29578 )
      zc=-rb*sin( beta / 57.29578 )
c
      call ndipold
c      btot=btot+(by*rho/ 57.29587)
c      write(1,*)'  ',theta,by
  10  continue
      bbx=bx
      bby=by
      bbz=bz
c      write(8,*)' returning ',bx,by,bz,bbx,bby,bbzc
c      save /biprt/,/block10old/,/blck20/,/blck21/,/blck22/,/blck23/,
c     &     /blck24/,/blck25/
      return
      end
c
      subroutine ndipold
c****
c****
c**** mtyp = 3 or 4
c**** this version of bfun is mainly for nonuniform field magnets
c**** the central field region is represented to 3'rd order on-and-
c**** off the midplane by analytic expressions. see slac no. 75
c**** fringe field regions represented by fermi type fall-off
c**** along with radial fall-off
c**** components of 'b' in fringe region evaluated by numerical methods
c****
c****
c**** the relationship between b0, ......... b12 and b(i,j) relative to
c**** axes (z,x) is given by
c****
c****
c**** b0  = b( 0, 0 )
c**** b1  = b( 1, 0 )
c**** b2  = b( 2, 0 )
c**** b3  = b( 1, 1 )
c**** b4  = b( 1,-1 )
c**** b5  = b( 0, 1 )
c**** b6  = b( 0, 2 )
c**** b7  = b( 0,-1 )
c**** b8  = b( 0,-2 )
c**** b9  = b(-1, 0 )
c**** b10 = b(-2, 0 )
c**** b11 = b(-1, 1 )
c**** b12 = b(-1,-1 )
c****
c****
      implicit real*8(a-h,o-z)
      real*8  ndx, k
      common  /block10old/  bx, by, bz, k, tc, dtc
      common  /blck20/  ndx,bet1,gama,delt,csc
      common  /blck21/  rca,dels,br,s2,s3,s4,s5,s6,s7,s8,scor
      common  /blck22/  d, dg, s, bf, bt
      common  /blck23/  c0, c1, c2, c3, c4, c5
      common  /blck24/  rb, xc, zc
      common  /blck25/  in, mtyp
      dimension tc(6), dtc(6)
      x = tc(1)
      y = tc(2)
      z = tc(3)
      dx = x - xc
      dz = z - zc
      rp =sqrt( dx**2 + dz**2 )
      dr = rp - rb
      go to ( 1, 2, 3, 14 ), in
    7 print 8, in, mtyp
      call exit(0)
    8 format (    '0 error -go to -  in bfun   in=', i3, '   mtyp=',i4 )
    2 drr1 = dr/rb
      drr2 = drr1*drr1
      drr3 = drr2*drr1
      drr4 = drr3*drr1
      if( y .ne. 0. )  go to 4
c****
c**** mid-plane uniform field region
c****
      bx = 0.
      by = 0.
      if( mtyp .eq. 3) by=
     1     bf* ( 1. - ndx*drr1 + bet1*drr2 + gama*drr3 + delt*drr4 )
      if( mtyp .eq. 4) by= bf/ (1. + ndx*drr1 )
      bz = 0.
      bt = by
      return
c****
c**** non mid-plane uniform field region
c****
    4 yr1 = y/rb
      yr2 = yr1*yr1
      yr3 = yr2*yr1
      yr4 = yr3*yr1
      rr1 = rb/rp
      rr2 = rr1*rr1
      rr3 = rr2*rr1
      if( mtyp .eq. 3 ) go to 11
      if( mtyp .eq. 4 ) go to 12
      go to 7
c****
c**** mtyp = 3
c****
   11 brr = bf*( ( -ndx + 2.*bet1*drr1 + 3.*gama*drr2 + 4.*delt*drr3 )
     1   *yr1 - (ndx*rr2 + 2.*bet1*rr1*(1.-rr1*drr1) +
     2   3.*gama*( 2. + 2.*rr1*drr1 - rr2*drr2 ) +
     3   4.*delt*( 6.*drr1 + 3.*rr1*drr2 - rr2*drr3 ))*yr3/6. )
      by = bf* ( 1. - ndx*drr1 + bet1*drr2 + gama*drr3 + delt*drr4 -
     1   .5*yr2*( -ndx*rr1 + 2.*bet1*( 1. + rr1*drr1) +
     2   3.*gama*drr1*( 2. + rr1*drr1) + 4.*delt*drr2*(3. + rr1*drr1) )
     3   + yr4*( -ndx*rr3 + 2.*bet1*( rr3*drr1 - rr2) +
     4   3.*gama*( 4.*rr1 - 2.*rr2*drr1 + rr3*drr2 ) +
     5   4.*delt*( 6. + 12.*rr1*drr1 - 3.*rr2*drr2 + rr3*drr3 ) )/24. )
      go to 13
c****
c**** mtyp = 4
c****
   12 dnr1 = 1. + ndx*drr1
      dnr2 = dnr1*dnr1
      dnr3 = dnr2*dnr1
      dnr4 = dnr3*dnr1
      dnr5 = dnr4*dnr1
      brr = bf*ndx*( -yr1/dnr2 + yr3*( 6.*ndx*ndx/dnr4 -
     1   2.*ndx*rr1/dnr3 - rr2/dnr2 ) /6.  )
      by = bf*( 1./dnr1 + .5*yr2*ndx*( -2.*ndx/dnr3 + rr1/dnr2) +
     2   yr4*ndx*( 24.*ndx**3 /dnr5 - 12.*ndx*ndx*rr1/dnr4 -
     3   2.*ndx*rr2/dnr3 - rr3/dnr2 ) /24.  )
c****
c****
   13 bx = brr*dx/rp
      bz = brr*dz/rp
      bt  =sqrt(bx*bx + by*by + bz*bz)
      return
c****
c****
    1 sine = -1.
      go to 5
    3 sine = 1.
    5 if( z  .gt. 0. ) dr = x * sine*csc
      call ndppold( b0, z, x, y, dr      )
      if( y  .ne. 0. )  go to 6
c****
c**** mid-plane fringing field region
c****
      bx = 0.
      by = b0
      bz = 0.
      bt   = b0
      return
c****
c**** non mid-plane fringing field region
c****
    6 if( z .gt. 0. )  go to 9
      dr1  =       (sqrt( dx**2 + (dz+dg)**2 ) - rb )
      dr2  =       (sqrt( dx**2 + (dz+2.*dg)**2 ) - rb )
      dr3  =       (sqrt( (dx+dg)**2 + (dz+dg)**2 )  - rb )
      dr4  =       (sqrt( (dx-dg)**2 + (dz+dg)**2 )  - rb )
      dr5  =       (sqrt( (dx+dg)**2 + dz**2 ) - rb )
      dr6  =       (sqrt( (dx+ 2.*dg)**2 + dz**2 ) - rb )
      dr7  =       (sqrt( (dx-dg)**2 + dz**2 ) - rb )
      dr8  =       (sqrt( (dx- 2.*dg)**2 + dz**2 ) - rb )
      dr9  =       (sqrt( dx**2 + (dz-dg)**2 ) - rb )
      dr10 =       (sqrt( dx**2 + (dz-2.*dg)**2 ) - rb )
      dr11 =       (sqrt( (dx+dg)**2 + (dz-dg)**2 )  - rb )
      dr12 =       (sqrt( (dx-dg)**2 + (dz-dg)**2 )  - rb )
      go to 10
    9 dr1  = sine* x*csc
      dr2  = dr1
      dr9  = dr1
      dr10 = dr1
      dr3  = sine* ( x + dg )*csc
      dr5  = dr3
      dr11 = dr3
      dr4  = sine*( x - dg )*csc
      dr7  = dr4
      dr12 = dr4
      dr6  = sine* ( x + 2.*dg )*csc
      dr8  = sine* ( x - 2.*dg )*csc
c****
c****
   10 call ndppold ( b1 , z + dg, x , y , dr1 )
      call ndppold ( b2 , z + 2.*dg, x , y , dr2 )
      call ndppold ( b3 , z + dg, x + dg , y , dr3 )
      call ndppold ( b4 , z + dg, x - dg , y , dr4 )
      call ndppold ( b5 , z , x + dg , y, dr5 )
      call ndppold ( b6 , z , x + 2.*dg , y , dr6 )
      call ndppold ( b7 , z , x - dg , y, dr7 )
      call ndppold ( b8 , z , x - 2.*dg , y , dr8 )
      call ndppold ( b9 , z - dg, x , y , dr9 )
      call ndppold ( b10, z - 2.*dg, x, y, dr10 )
      call ndppold ( b11, z - dg, x + dg , y , dr11 )
      call ndppold ( b12, z - dg, x - dg , y , dr12 )
      yg1 = y/dg
      yg2 = yg1**2
      yg3 = yg1**3
      yg4 = yg1**4
      bx = yg1 * ( (b5-b7)*2./3. - (b6-b8)/12. )  +
     1     yg3*( (b5-b7)/6. - (b6-b8)/12. -
     2     (b3 + b11 - b4 - b12 - 2.*b5 + 2.*b7 ) / 12. )
      by = b0 - yg2*( ( b1 + b9 + b5 + b7 - 4.*b0 ) *2./3. -
     1     ( b2 + b10 + b6 + b8 - 4.*b0 ) / 24. ) +
     2     yg4* (-( b1 + b9 + b5 + b7 - 4.*b0 ) / 6. +
     3     ( b2 + b10 + b6 + b8 - 4.*b0 ) / 24. +
     4     ( b3 + b11 + b4 + b12 - 2.*b1 - 2.*b9 -
     5     2.*b5 - 2.*b7 + 4.*b0 ) / 12. )
      bz = yg1*( (b1 - b9 ) *2./3. - ( b2 - b10 ) /12. ) +
     1     yg3*( ( b1 - b9 ) / 6. - ( b2 - b10 ) / 12. -
     2     ( b3 + b4 - b11 - b12 - 2.*b1 + 2.*b9 ) / 12.  )
      bt  =sqrt(bx*bx + by*by + bz*bz)
      return
   14 bx = 0.
      by = br
      bz = 0.
      bt = br
c      save /block10old/,/blck20/,/blck21/,/blck22/,/blck23/,
c    &     /blck24/,/blck25/
      return
      end
      subroutine  ndppold ( bfld, z, x, y , dr )
c****
c****
c****
c****
      implicit real*8(a-h,o-z)
      real*8  ndx, k
      common  /blck10old/  bx, by, bz, k, tc, dtc
      common  /blck20/  ndx,bet1,gama,delt,csc
      common  /blck21/  rca,dels,br,s2,s3,s4,s5,s6,s7,s8,scor
      common  /blck22/  d, dg, s, bf, bt
      common  /blck23/  c0, c1, c2, c3, c4, c5
      common  /blck24/  rb, xc, zc
      common  /blck25/  in, mtyp
      dimension tc(6), dtc(6)
      drr1 = dr/rb
      drr2 = drr1*drr1
      drr3 = drr2*drr1
      drr4 = drr3*drr1
c****
c**** mtyp    :                  modified iterative procedure
c****
      xp = x
      xp2 = xp*xp
      xp3 = xp2*xp
      xp4 = xp3 * xp
      zp = -(s2*xp2 + s3*xp3 + s4*xp4 + s5*xp4*xp + s6*xp4*xp2 +
     1       s7*xp4*xp3 + s8*xp4*xp4 )
      az = (z-zp)/10.d0
      azmax = sqrt(  x*x + z*z  )
      if( az  .gt.  azmax  )  az = azmax
      zsign = z-zp
      rinv4 = 0.
      do 11 i=1,21
      xp   = x + az*(i-11)
      xp2 = xp*xp
      xp3 = xp2*xp
      xp4 = xp3*xp
      zp = -(s2*xp2 + s3*xp3 + s4*xp4 + s5*xp4*xp + s6*xp4*xp2 +
     1       s7*xp4*xp3 + s8*xp4*xp4 )
      xxp = x-xp
      zzp = z-zp
      dd =            xxp*xxp + zzp*zzp
      if( dd  .lt.  1.d-15 )  dd = 1.d-15
      if( dd  .gt.  1.d15  )  dd = 1.d15
      rinv4 = rinv4 + 1.0d0 / (dd*dd )
   11 continue
      dp = sqrt( 1.d0/rinv4 )
      dp = sqrt( dp )
      s = 1.9023d0* sign( 1.d0, zsign ) * dp/d + dels
c****
c**** first guess for closest point is
c****
c*    xp = x
c*    xp2 = xp*xp
c*    xp3 = xp2*xp
c*    xp4 = xp3*xp
c****
c**** calculate zp on curve for corresponding xp
c****
c*    zp = -( s2*xp2 + s3*xp3 + s4*xp4 + s5*xp4*xp + s6*xp4*xp2 +
c*   1   s7*xp4*xp3 + s8*xp4*xp4 )
c*    zsign = z-zp
c****
c**** slope of curve at xp, zp
c****
c*    do 4 i=1,3
c*    dzdxc = -(2.*s2*xp + 3.*s3*xp2+ 4.*s4*xp3 + 5.*s5*xp4 +
c*   1   6.*s6*xp4*xp + 7.*s7*xp4*xp2 + 8.*s8*xp4*xp3 )
c****
c**** next approximation to closest point is
c****
c*    xp = ( dzdxc*(z-zp)  +  dzdxc*dzdxc*xp + x ) / (1.+dzdxc*dzdxc)
c*    if( i  .eq.  1 )  xp = (3.*xp +  x ) / 4.
c*    xp2 = xp*xp
c*    xp3 = xp2*xp
c*    xp4 = xp3*xp
c*    zp = -( s2*xp2 + s3*xp3 + s4*xp4 + s5*xp4*xp + s6*xp4*xp2 +
c*   1   s7*xp4*xp3 + s8*xp4*xp4 )
c*  4 continue
c*    xxp = x-xp
c*    zzp = z-zp
c*    s = sign( 1.d0,zsign) * sqrt( xxp*xxp + zzp*zzp) / d + dels
c****
c****
c****
c****
      cs=c0+s*(c1+s*(c2+s*(c3+s*(c4+s*c5))))
      if( abs(cs)  .gt.  70.  )  cs =sign( 70.d0 ,cs  )
      e=exp(cs)
      p0 = 1.0 + e
      db=bf-br
      bfld = 0.
      if( mtyp .eq. 3 ) bfld =
     1       br +( 1. - ndx*drr1 + bet1*drr2+gama*drr3+delt*drr4)*db/p0
      if( mtyp .eq. 4 ) bfld = br + ( 1./(1. +ndx*drr1) )*db/p0
c****
c**** print 100, x, y, z,  dr, s, bfld
c*100 format( 1p6d15.4 )
c****
c      save /block10old/,/blck20/,/blck21/,/blck22/,/blck23/,
c     &     /blck24/,/blck25/
      return
      end
      subroutine qqfringe(x,y,z,dy,bx,by,bz,fact1,fact2,qrad1,qrad2)
c*************************************************************************
c  qqfringe -                                                            *
c      calculates combined fringe fields of two adjacent quads, q1 & q2. *
c      the exit of q1 is located at (x, y, z) and the entrance of q2 is  *
c      dy downstream from the exit of q1.                                *
c      facti describes qi - bquad, bhex, boct, bdec, bddec               *
c      qradi radius of qi                                                *
c*************************************************************************
      real fact1(5),fact2(5),b(3)
      bx=0.
      by=0.
      bz=0.
c  exit of q1 field
      xdum=x
      ydum=y
      zdum=z
      call mpoles(3,xdum,zdum,ydum,fact1,qrad1,b(1),b(3),b(2))
      bx=b(1)
      by=b(2)
      bz=b(3)
c      write(9,100)xdum,ydum,zdum,bx,by,bz
c 100  format(1x,'x,y,z=',3f10.3,' bi=',3f11.6)
c 101  format(1x,'      ',3f10.3,'    ',3f11.6)
c  entrance of q2 field
      ydum=dy-ydum
c      write(9,102)(fact2(i),i=1,5),qrad2
c 102  format(1x,'calling mpoles  ',6f8.2)
      call mpoles(1,xdum,zdum,ydum,fact2,qrad2,b(1),b(3),b(2))
      bx=bx+b(1)
      by=by-b(2)
      bz=bz+b(3)
c      write(9,101)xdum,ydum,zdum,bx,by,bz
      return
      end
      subroutine qdfringe(x,y,z,dy,bx,by,bz,fact1,fact2,qrad1,qrad2)
c*************************************************************************
c  qdfringe -                                                            *
c      calculates combined fringe fields of one quad, q1, followed by a  *
c      dipole, d.                                                        *
c      the exit of q1 is located at x=y=z=0 and the entrance of d is     *
c      dy downstream from the exit of q1 and rotated by angle th.        *
c      fact1 describes q1 - bquad, bhex, boct, bdec, bddec               *
c      fact2 describes dipole - bdip, theta
c*************************************************************************
      real fact1(5),fact2(5),b(3),bpri(3)
c  pifac=pi/180
      data pifac/1.745329e-02/
      th=fact2(2)*pifac
      bx=0.
      by=0.
      bz=0.
c  exit of q1 field
      xdum=x
      ydum=y
      zdum=z
      call mpoles(3,xdum,zdum,ydum,fact1,qrad1,b(1),b(3),b(2))
      bx=b(1)
      by=-b(2)
      bz=b(3)
c  entrance dipole field
c   translate and reflect y axis
      ypri=y-dy
      xpri=x
      zpri=z
c   rotate by th
      xdpri=(xpri*cos(th))+(ypri*sin(th))
      ydpri=(-xpri*sin(th))+(ypri*cos(th))
      zdpri=z
      ydpri=-ydpri
      call dfringe(xdpri,ydpri,zdpri,fact2(1),qrad2,b(1),b(2),b(3))
      b(2)=-b(2)
      bpri(1)=(b(1)*cos(-th))+(b(2)*sin(-th))
      bpri(2)=(-b(1)*sin(-th))+(b(2)*cos(-th))
      bpri(3)=b(3)
      bx=bx+bpri(1)
      by=by+bpri(2)
      bz=bz+bpri(3)
      return
      end
      
      subroutine ddfringe(x,y,z,dl,bx,by,bz,fact1,fact2,qrad1,qrad2)
c*************************************************************************
c  ddfringe -                                                            *
c      calculates combined fringe fields of two dipoles                  *
c      both d1 and d2 pole faces are rotated by thi with respect to the  *
c      optic axis making a total bend from the center of d1 to the center*
c      of d2 of th=th1+th2. (see diagrams in notebook 29 april 1987)     *
c      the exit of d1 is located at x=y=z=0 and the entrance of d2 is    *
c      dl downstream from the exit of d1 and rotated by angle th.        *
c      fact1 describes d1 - bdip, theta1                                 *
c      fact2 describes d2 - bdip, theta2                                 *
c*************************************************************************
      real fact1(5),fact2(5),b(3)
c  pifac=pi/180
      data pifac/1.745329e-02/
      th2=fact2(2)*pifac
      th1=fact1(2)*pifac
      bx=0.
      by=0.
      bz=0.
c  exit of d1 field
      xdum=x
      ydum=y
      zdum=z
      call dfringe(xdum,ydum,zdum,fact1(1),qrad,b(1),b(2),b(3))
      by=-by
c  entrance dipole field
c   translate x & y axes
      dy=dl*sin(th1)
      dx=dl*cos(th1)
      ypri=y-dy
      xpri=x+dx
      zpri=z
c   rotate by th
      th=th1+th2
      xdpri=(xpri*cos(th))+(ypri*sin(th))
      ydpri=(-xpri*sin(th))+(ypri*cos(th))
      zdpri=z
      call dfringe(xdpri,ydpri,zdpri,fact2(1),qrad,b(1),b(2),b(3))
      bx=bx+b(1)
      by=by+b(2)
      bz=bz+b(3)
      return
      end
      
      subroutine dqfringe(x,y,z,dl,bx,by,bz,fact1,fact2,qrad1,qrad2)
c***************************************************************************
c    dqfringe - computes overlap fringe fields for a dipole quadrupole     *
c    pair.  the dipole center exit is located at x=y=z=0. the optic axis   *
c    makes an angle th with the dipole face.  the quad is located a        *
c    distance dl along the optic axis from the dipole.                     *
c    diagrams in notebook 29 april 1987                                    *
c                                                                          *
c      fact1 describes dipole - bdip, th                                   *
c      fact2 describes quad - bquad, bhex, boct, bdec, bddec               *
c***************************************************************************
      real fact1(5),fact2(5),b(3)
c  pifac=pi/180
      data pifac/1.745329e-02/
      th=fact1(2)*pifac
      bx=0.
      by=0.
      bz=0.
c  exit of d1 field
      xdum=x
      ydum=y
      zdum=z
c next: bug corrected (?) P.V 4/10/97
c      call dfringe(xdum,ydum,zdum,fact1(1),qrad,bx,by,bz)
      call dfringe(xdum,ydum,zdum,fact1(1),qrad1,bx,by,bz)
      by=-by
c  entrance quad field
c   translate x & y axes
      dy=dl*sin(th)
      dx=dl*cos(th)
      ypri=y-dy
      xpri=x+dx
      zpri=z
c   rotate by th
      xdpri=(xpri*cos(th))+(ypri*sin(th))
      ydpri=(-xpri*sin(th))+(ypri*cos(th))
      zdpri=z
      ydpri=-ydpri
      call mpoles(1,xdpri,zdpri,ydpri,fact2,qrad2,b(1),b(3),b(2))
      bx=bx+b(1)
      by=by+b(2)
      bz=bz+b(3)
      return
      end
      subroutine dqfringe2(x,y,z,dl,bx,by,bz,fact1,fact2,qrad1,qrad2)
c***************************************************************************
c    dqfringe - computes overlap fringe fields for a dipole quadrupole     *
c    pair.  the dipole center exit is located at x=y=z=0. the optic axis   *
c    makes an angle th with the dipole face.  the quad is located a        *
c    distance dl along the optic axis from the dipole.                     *
c    diagrams in notebook 29 april 1987                                    *
c    modified version of dqfringe. altered to include dipolert 8/23/89 jjl *
c                                                                          *
c      fact1 describes dipole - bdip, th                                   *
c      fact2 describes quad - bquad, bhex, boct, bdec, bddec               *
c***************************************************************************
      real fact1(5),fact2(5),b(3),bpri(3)
c  pifac=pi/180
      data pifac/1.745329e-02/,ten/10./
      th=fact1(2)*pifac
c      write(7,*)' dl, th' ,dl,th
      bx=0.
      by=0.
      bz=0.
c  exit of d1 field
        xdum=x/ten
        ydum=y/ten
        zdum=z/ten
c      write(7,*)' xdum,ydum,zdum ',xdum,ydum,zdum
        call dipolert2(xdum,zdum,ydum,bx,bz,by)
c  get the right sign on the field
        bx=-bx
        by=-by
        bz=-bz
c      write(7,*)' dipole',bx,by,bz
c  entrance quad field
c   translate x & y axes
      dx=dl*sin(th)
      dy=dl*cos(th)
      ypri=y-dy
      xpri=x+dx
      zpri=z
c      write(7,*)' primed coord ',xpri,ypri,zpri
c   rotate by th
      xdpri=(xpri*cos(th))+(ypri*sin(th))
      ydpri=(-xpri*sin(th))+(ypri*cos(th))
      zdpri=z
c      write(7,*)' dprimed coord ',xdpri,ydpri,zdpri
      ydpri=-ydpri
      call mpoles(1,xdpri,zdpri,ydpri,fact2,qrad2,b(1),b(3),b(2))
      b(2)=-b(2)
c      write(7,*)' mpoles',(b(i),i=1,3)
      bpri(1)=(b(1)*cos(-th))+(b(2)*sin(-th))
      bpri(2)=(-b(1)*sin(-th))+(b(2)*cos(-th))
      bpri(3)=b(3)
      bx=bx+bpri(1)
      by=by+bpri(2)
      bz=bz+bpri(3)
      return
      end
      subroutine qdfringe2(x,y,z,dy,bx,by,bz,fact1,fact2,qrad1,qrad2)
c*************************************************************************
c  qdfringe -                                                            *
c      calculates combined fringe fields of one quad, q1, followed by a  *
c      dipole, d.                                                        *
c      the exit of q1 is located at x=y=z=0 and the entrance of d is     *
c      dy downstream from the exit of q1 and rotated by angle th.        *
c      fact1 describes q1 - bquad, bhex, boct, bdec, bddec               *
c      fact2(1)=dx
c      fact2(2)=dy
c      fact2(3)=th
c*************************************************************************
      real fact1(5),fact2(5),b(3),bpri(3)
c  pifac=pi/180
      data pifac/1.745329e-02/,ten/10./
      th=fact2(3)*pifac
      dx=fact2(1)
      dy=fact2(2)
      bx=0.
      by=0.
      bz=0.
c  exit of q1 field
      xdum=x
      ydum=y
      zdum=z
      call mpoles(3,xdum,zdum,ydum,fact1,qrad1,b(1),b(3),b(2))
      bx=b(1)
      by=b(2)
      bz=b(3)
c  entrance dipole field
c   translate and reflect y axis
      ypri=y-dy
      xpri=x-dx
      zpri=z
c   rotate by th
      xdpri=((xpri*cos(th))+(ypri*sin(th)))/ten
      ydpri=((-xpri*sin(th))+(ypri*cos(th)))/ten
      zdpri=z/ten
c      call dfringe(xdpri,ydpri,zdpri,fact2(1),qrad2,b(1),b(2),b(3))
        call dipolert2(xdpri,zdpri,ydpri,b(1),b(3),b(2))
c  get the right sign on the field
      do 10 i=1,3
  10   b(i)=-b(i)
      bpri(1)=(b(1)*cos(-th))+(b(2)*sin(-th))
      bpri(2)=(-b(1)*sin(-th))+(b(2)*cos(-th))
      bpri(3)=b(3)
      bx=bx+bpri(1)
      by=by+bpri(2)
      bz=bz+bpri(3)
      return
      end

c
c  dipolert2- analytic calculation of dipole magnetic fields
c            stolen from raytrace adapted for use in snake
c                                         -jjl  8/21/89
c    data(i) are read in snake in ifield
c    dipolert2 is called from snake in analyf
c                snake provides x,y,z in the c-axis system in cm
c                dipolert2 returns bx,by,bz in the same system
c
      subroutine dipolert2(x,y,z,bbx,bby,bbz)
      implicit real*8(a-h,o-z)
      real x,y,z,bbx,bby,bbz,data(75),pii
      real*8  ndx, k
      common  /diprt/data
      common  /blck10/  bx, by, bz, k, tc, dtc
      common  /blck20/  ndx,bet1,gama,delt,csc
      common  /blck21/  rca,dels,br,s2,s3,s4,s5,s6,s7,s8,scor
      common  /blck22/  d, dg, s, bf, bt
      common  /blck23/  c0, c1, c2, c3, c4, c5
      common  /blck24/  rb, xc, zc
      common  /blck25/  in, mtyp
      dimension tc(2,6), dtc(6), bb(2,0:12)
	  data pii/0.0174532925/
      lf1  = dble(data(  1  ))
      lu1  = dble(data(  2  ))
      lf2  = dble(data(  3  ))
      dg   = dble(data(  4  ))
      mtyp = dble(data(  5  ))
      a    = dble(data(  11  ))
      b    = dble(data(  12  ))
      d    = dble(data(  13  ))
      rb   = dble(data(  14  ))
      bf   = dble(data(  15  ))
      phi  = dble(data(  16  ))
      alpha= dble(data(  17  ))
      beta = dble(data(  18  ))
      ndx  = dble(data(  19  ))
      bet1 = dble(data(  20  ))
      gama = dble(data(  21  ))
      delt = dble(data(  22  ))
      z11  = dble(data(  25  ))
      z12  = dble(data(  26  ))
      z21  = dble(data(  27  ))
      z22  = dble(data(  28  ))
      br1  = dble(data(  41  ))
      br2  = dble(data(  42  ))
      xcr1 = dble(data(  43  ))
      xcr2 = dble(data(  44  ))
c      write(8,*)' in dipolert2 1',x,y,z
c**** transform from c-axis to b-axis
c****
      xpric=x+(2.*rb*dsin(pii*phi/2.)*dsin(pii*((phi/2.)-beta)))
      zpric=z+(2.*rb*dsin(pii*phi/2.)*dcos(pii*((phi/2.)-beta)))
c bug check -jjl 12/14/95
c      xpric=x+(2.*rb*dsind(phi/2.)*dsind((phi/2.)-alpha))
c      zpric=z+(2.*rb*dsind(phi/2.)*dcosd((phi/2.)-alpha))
c end bug check
      cosa=dcos(pii*(phi-alpha-beta))
      sina=dsin(pii*(phi-alpha-beta))
      tc(1,1)=(-xpric*cosa)+(zpric*sina)
      tc(1,2)=dble(y)
      tc(1,3)=(-xpric*sina)-(zpric*cosa)
c**** transform to second vfb coord system
c***
      copab =dcos( (phi-alpha-beta)/57.29578)
      sipab =dsin( (phi-alpha-beta)/57.29578)
      cospb =dcos( (phi/2.-beta)/57.29578 )
      sinpb =dsin( (phi/2.-beta)/57.29578 )
      sip2 =dsin( (phi/2.)/57.29578 )
      xt = tc(1,1)
      zt = tc(1,3)
      vxt = tc(1,4)
      vzt = tc(1,6)
      tc(2,3) = - zt  *copab +  xt  *sipab -2.*rb*sip2*cospb
      tc(2,1) = - zt  *sipab -  xt  *copab -2.*rb*sip2*sinpb
      tc(2,6) = - vzt *copab +  vxt *sipab
      tc(2,4) = - vzt *sipab -  vxt *copab
      tc(2,2)=dble(y)
c      write(8,*)' in dipolert2 1',(tc(i),i=1,3)
c
c  test for region
c
      if(tc(1,3).gt.z11)then
       bx=0.
       by=0.
       bz=0.
       go to 10
      endif
      if((tc(1,3).ge.z12).and.(tc(2,3).lt.z21)) go to 1  !entrance fringe
      if((tc(1,3).ge.z12).and.(tc(2,3).ge.z21)) go to 25 !overlapping entr&exit
      if(tc(1,3).lt.z12) go to 2                         !uniform field or exit
c      go to 15
c****
c**** in designates magnet regions for bfun
c****
  1   in = 1
      xc= rb*dcos( alpha/ 57.29578 )
      zc=-rb*dsin( alpha/ 57.29578 )
c****
      c0   = dble(data(  29  ))
      c1   = dble(data(  30  ))
      c2   = dble(data(  31  ))
      c3   = dble(data(  32  ))
      c4   = dble(data(  33  ))
      c5   = dble(data(  34  ))
      dels = dble(data(  45  ))
      rca  = dble(data(  47  ))
      csc = dcos( alpha/57.29578 )
      scor = dble(data( 49 ))
      s2   = dble(data(  51  )) / rb    + rca/2.d0
      s3   = dble(data(  52  )) / rb**2
      s4   = dble(data(  53  )) / rb**3 + rca**3/8.d0
      s5   = dble(data(  54  )) / rb**4
      s6   = dble(data(  55  )) / rb**5 + rca**5/16.d0
      s7   = dble(data(  56  )) / rb**6
      s8   = dble(data(  57  )) / rb**7 + rca**7/25.6d0
c
      call ndip(1)
c      write(8,*)' return from ndip',bx,by,bz
c**** transform to c-axis system
c***
      bxdum=bx
      bzdum=bz
      bx=(-bxdum*cosa)-(bzdum*sina)
      bz=(+bxdum*sina)-(bzdum*cosa)
      go to 10
  2   continue
      if(tc(2,3).lt.z21)go to 3      !uniform field
      if(tc(2,3).gt.z22)then         !outside of fringe region
       bx=0.
       by=0.
       bz=0.
       go to 10
      endif
      if(tc(2,3).le.z22)go to 4      !exit fringe
c****
c****
c**** uniform field integration region
c****
c****
  3   in = 2
ctest      xc=-rb*dcos( beta / 57.29578 )
      xc= rb*dcos( beta / 57.29578 )
      zc=-rb*dsin( beta / 57.29578 )
c
      call ndip(1)
c      btot=btot+(by*rho/ 57.29587)
c      if(theta.eq.22.5d+00)br0=by
c      write(1,*)' ',theta,by
      go to 10
c***
c***
c**** setup for second fringe field and integration
c****
c****
  4   br   = br2
      c0   = dble(data(  35  ))
      c1   = dble(data(  36  ))
      c2   = dble(data(  37  ))
      c3   = dble(data(  38  ))
      c4   = dble(data(  39  ))
      c5   = dble(data(  40  ))
      dels = dble(data(  46  ))
      rca  = dble(data(  48  ))
      scor = dble(data( 50 ))
      csc = dcos( beta /57.29578 )
      s2   = dble(data(  58  )) / rb    + rca/2.d0
      s3   = dble(data(  59  )) / rb**2
      s4   = dble(data(  60  )) / rb**3 + rca**3/8.d0
      s5   = dble(data(  61  )) / rb**4
      s6   = dble(data(  62  )) / rb**5 + rca**5/16.d0
      s7   = dble(data(  63  )) / rb**6
      s8   = dble(data(  64  )) / rb**7 + rca**7/25.6d0
      in = 3
      xc=-rb*dcos( beta / 57.29578 )
      zc=-rb*dsin( beta / 57.29578 )
c
      call ndip(2)
c      btot=btot+(by*rho/ 57.29587)
c      write(1,*)'  ',theta,by
  10  continue
      bbx=bx
      bby=by
      bbz=bz
c      write(8,*)' returning ',bx,by,bz,bbx,bby,bbz
c      save /biprt/,/block10old/,/blck20/,/blck21/,/blck22/,/blck23/,
c     &     /blck24/,/blck25/
      return
c
c overlapping entrance and exit fringe fields
c
  25  in = 1                         !setup for entrance field
      xc= rb*dcos( alpha/ 57.29578 )
      zc=-rb*dsin( alpha/ 57.29578 )
c****
      c0   = dble(data(  29  ))
      c1   = dble(data(  30  ))
      c2   = dble(data(  31  ))
      c3   = dble(data(  32  ))
      c4   = dble(data(  33  ))
      c5   = dble(data(  34  ))
      dels = dble(data(  45  ))
      rca  = dble(data(  47  ))
      csc = dcos( alpha/57.29578 )
      scor = dble(data(49 ))
      s2   = dble(data(  51  )) / rb    + rca/2.d0
      s3   = dble(data(  52  )) / rb**2
      s4   = dble(data(  53  )) / rb**3 + rca**3/8.d0
      s5   = dble(data(  54  )) / rb**4
      s6   = dble(data(  55  )) / rb**5 + rca**5/16.d0
      s7   = dble(data(  56  )) / rb**6
      s8   = dble(data(  57  )) / rb**7 + rca**7/25.6d0
c
      call nndip(1,bb)
c***
c***
c**** setup for exit fringe field
c****
c****
      br   = br2
      c0   = dble(data(  35  ))
      c1   = dble(data(  36  ))
      c2   = dble(data(  37  ))
      c3   = dble(data(  38  ))
      c4   = dble(data(  39  ))
      c5   = dble(data(  40  ))
      dels = dble(data(  46  ))
      rca  = dble(data(  48  ))
      scor = dble(data(50 ))
      csc = dcos( beta /57.29578 )
      s2   = dble(data(  58  )) / rb    + rca/2.d0
      s3   = dble(data(  59  )) / rb**2
      s4   = dble(data(  60  )) / rb**3 + rca**3/8.d0
      s5   = dble(data(  61  )) / rb**4
      s6   = dble(data(  62  )) / rb**5 + rca**5/16.d0
      s7   = dble(data(  63  )) / rb**6
      s8   = dble(data(  64  )) / rb**7 + rca**7/25.6d0
      in = 3
      xc=-rb*dcos( beta / 57.29578 )
      zc=-rb*dsin( beta / 57.29578 )
      call nndip(2,bb)
      if(y.eq.0.)then
       bx=0.
       bz=0.
       by=bb(1,0)+bb(2,0)-bf
      else
      do 26 i=0,12
       bb(1,i)=bb(1,i)+bb(2,i)-bf
  26  continue
      yg1 = y/dg
      yg2 = yg1**2
      yg3 = yg1**3
      yg4 = yg1**4
c        1         2         3         4         5         6         7
      bx = yg1 * ( (bb(1,5)-bb(1,7))*2./3. - (bb(1,6)-bb(1,8))/12. ) +
     1     yg3*( (bb(1,5)-bb(1,7))/6. - (bb(1,6)-bb(1,8))/12. -
     2     (bb(1,3) + bb(1,11) - bb(1,4) - bb(1,12)
     3      - 2.*bb(1,5) + 2.*bb(1,7) ) / 12. )
      by = bb(1,0) - yg2*( ( bb(1,1) + bb(1,9) + bb(1,5) + bb(1,7)
     1     - 4.*bb(1,0) ) *2./3. -
     2     ( bb(1,2) + bb(1,10) + bb(1,6) + bb(1,8)-4.*bb(1,0))/24.) +
     3     yg4*(-(bb(1,1)+ bb(1,9)+ bb(1,5)+ bb(1,7)- 4.*bb(1,0) )/6.+
     4     ( bb(1,2)+ bb(1,10)+ bb(1,6) + bb(1,8) - 4.*bb(1,0))/24. +
     5     ( bb(1,3)+bb(1,11)+bb(1,4) +bb(1,12)-2.*bb(1,1)-2.*bb(1,9)-
     6     2.*bb(1,5) - 2.*bb(1,7) + 4.*bb(1,0) ) / 12. )
      bz = yg1*((bb(1,1)-bb(1,9))*2./3. - (bb(1,2)-bb(1,10) ) /12. ) +
     1     yg3*( (bb(1,1)- bb(1,9))/6. - (bb(1,2) - bb(1,10) ) / 12. -
     2     ( bb(1,3) + bb(1,4) - bb(1,11) - bb(1,12) 
     3      - 2.*bb(1,1) + 2.*bb(1,9) ) / 12.  )
      bt  =dsqrt(bx*bx + by*by + bz*bz)
      endif
      go to 10
      end
c
      subroutine ndip(l)
c****
c****
c**** mtyp = 3 or 4
c**** this version of bfun is mainly for nonuniform field magnets
c**** the central field region is represented to 3'rd order on-and-
c**** off the midplane by analytic expressions. see slac no. 75
c**** fringe field regions represented by fermi type fall-off
c**** along with radial fall-off
c**** components of 'b' in fringe region evaluated by numerical methods
c****
c****
c**** the relationship between b0, ......... b12 and b(i,j) relative to
c**** axes (z,x) is given by
c****
c****
c**** b0  = b( 0, 0 )
c**** b1  = b( 1, 0 )
c**** b2  = b( 2, 0 )
c**** b3  = b( 1, 1 )
c**** b4  = b( 1,-1 )
c**** b5  = b( 0, 1 )
c**** b6  = b( 0, 2 )
c**** b7  = b( 0,-1 )
c**** b8  = b( 0,-2 )
c**** b9  = b(-1, 0 )
c**** b10 = b(-2, 0 )
c**** b11 = b(-1, 1 )
c**** b12 = b(-1,-1 )
c****
c****
      implicit real*8(a-h,o-z)
      real*8  ndx, k
      common  /blck10/  bx, by, bz, k, tc, dtc
      common  /blck20/  ndx,bet1,gama,delt,csc
      common  /blck21/  rca,dels,br,s2,s3,s4,s5,s6,s7,s8,scor
      common  /blck22/  d, dg, s, bf, bt
      common  /blck23/  c0, c1, c2, c3, c4, c5
      common  /blck24/  rb, xc, zc
      common  /blck25/  in, mtyp
      dimension tc(2,6), dtc(6)
      x = tc(l,1)
      y = tc(l,2)
      z = tc(l,3)
      dx = x - xc
      dz = z - zc
      rp =dsqrt( dx**2 + dz**2 )
      dr = rp - rb
      go to ( 1, 2, 3, 14 ), in
    7 print 8, in, mtyp
      call exit(0)
    8 format (    '0 error -go to -  in bfun   in=', i3, '   mtyp=',i4 )
    2 drr1 = dr/rb
      drr2 = drr1*drr1
      drr3 = drr2*drr1
      drr4 = drr3*drr1
      if( y .ne. 0. )  go to 4
c****
c**** mid-plane uniform field region
c****
      bx = 0.
      by = 0.
      if( mtyp .eq. 3) by=
     1     bf* ( 1. - ndx*drr1 + bet1*drr2 + gama*drr3 + delt*drr4 )
      if( mtyp .eq. 4) by= bf/ (1. + ndx*drr1 )
      bz = 0.
      bt = by
      return
c****
c**** non mid-plane uniform field region
c****
    4 yr1 = y/rb
      yr2 = yr1*yr1
      yr3 = yr2*yr1
      yr4 = yr3*yr1
      rr1 = rb/rp
      rr2 = rr1*rr1
      rr3 = rr2*rr1
      if( mtyp .eq. 3 ) go to 11
      if( mtyp .eq. 4 ) go to 12
      go to 7
c****
c**** mtyp = 3
c****
   11 brr = bf*( ( -ndx + 2.*bet1*drr1 + 3.*gama*drr2 + 4.*delt*drr3 )
     1   *yr1 - (ndx*rr2 + 2.*bet1*rr1*(1.-rr1*drr1) +
     2   3.*gama*( 2. + 2.*rr1*drr1 - rr2*drr2 ) +
     3   4.*delt*( 6.*drr1 + 3.*rr1*drr2 - rr2*drr3 ))*yr3/6. )
      by = bf* ( 1. - ndx*drr1 + bet1*drr2 + gama*drr3 + delt*drr4 -
     1   .5*yr2*( -ndx*rr1 + 2.*bet1*( 1. + rr1*drr1) +
     2   3.*gama*drr1*( 2. + rr1*drr1) + 4.*delt*drr2*(3. + rr1*drr1) )
     3   + yr4*( -ndx*rr3 + 2.*bet1*( rr3*drr1 - rr2) +
     4   3.*gama*( 4.*rr1 - 2.*rr2*drr1 + rr3*drr2 ) +
     5   4.*delt*( 6. + 12.*rr1*drr1 - 3.*rr2*drr2 + rr3*drr3 ) )/24. )
      go to 13
c****
c**** mtyp = 4
c****
   12 dnr1 = 1. + ndx*drr1
      dnr2 = dnr1*dnr1
      dnr3 = dnr2*dnr1
      dnr4 = dnr3*dnr1
      dnr5 = dnr4*dnr1
      brr = bf*ndx*( -yr1/dnr2 + yr3*( 6.*ndx*ndx/dnr4 -
     1   2.*ndx*rr1/dnr3 - rr2/dnr2 ) /6.  )
      by = bf*( 1./dnr1 + .5*yr2*ndx*( -2.*ndx/dnr3 + rr1/dnr2) +
     2   yr4*ndx*( 24.*ndx**3 /dnr5 - 12.*ndx*ndx*rr1/dnr4 -
     3   2.*ndx*rr2/dnr3 - rr3/dnr2 ) /24.  )
c****
c****
   13 bx = brr*dx/rp
      bz = brr*dz/rp
      bt  =dsqrt(bx*bx + by*by + bz*bz)
      return
c****
c****
    1 sine = -1.
      go to 5
    3 sine = 1.
    5 if( z  .gt. 0. ) dr = x * sine*csc
      call ndpp( b0, z, x, y, dr      )
      if( y  .ne. 0. )  go to 6
c****
c**** mid-plane fringing field region
c****
      bx = 0.
      by = b0
      bz = 0.
      bt   = b0
      return
c****
c**** non mid-plane fringing field region
c****
    6 if( z .gt. 0. )  go to 9
      dr1  =       (dsqrt( dx**2 + (dz+dg)**2 ) - rb )
      dr2  =       (dsqrt( dx**2 + (dz+2.*dg)**2 ) - rb )
      dr3  =       (dsqrt( (dx+dg)**2 + (dz+dg)**2 )  - rb )
      dr4  =       (dsqrt( (dx-dg)**2 + (dz+dg)**2 )  - rb )
      dr5  =       (dsqrt( (dx+dg)**2 + dz**2 ) - rb )
      dr6  =       (dsqrt( (dx+ 2.*dg)**2 + dz**2 ) - rb )
      dr7  =       (dsqrt( (dx-dg)**2 + dz**2 ) - rb )
      dr8  =       (dsqrt( (dx- 2.*dg)**2 + dz**2 ) - rb )
      dr9  =       (dsqrt( dx**2 + (dz-dg)**2 ) - rb )
      dr10 =       (dsqrt( dx**2 + (dz-2.*dg)**2 ) - rb )
      dr11 =       (dsqrt( (dx+dg)**2 + (dz-dg)**2 )  - rb )
      dr12 =       (dsqrt( (dx-dg)**2 + (dz-dg)**2 )  - rb )
      go to 10
    9 dr1  = sine* x*csc
      dr2  = dr1
      dr9  = dr1
      dr10 = dr1
      dr3  = sine* ( x + dg )*csc
      dr5  = dr3
      dr11 = dr3
      dr4  = sine*( x - dg )*csc
      dr7  = dr4
      dr12 = dr4
      dr6  = sine* ( x + 2.*dg )*csc
      dr8  = sine* ( x - 2.*dg )*csc
c****
c****
   10 call ndpp ( b1 , z + dg, x , y , dr1 )
      call ndpp ( b2 , z + 2.*dg, x , y , dr2 )
      call ndpp ( b3 , z + dg, x + dg , y , dr3 )
      call ndpp ( b4 , z + dg, x - dg , y , dr4 )
      call ndpp ( b5 , z , x + dg , y, dr5 )
      call ndpp ( b6 , z , x + 2.*dg , y , dr6 )
      call ndpp ( b7 , z , x - dg , y, dr7 )
      call ndpp ( b8 , z , x - 2.*dg , y , dr8 )
      call ndpp ( b9 , z - dg, x , y , dr9 )
      call ndpp ( b10, z - 2.*dg, x, y, dr10 )
      call ndpp ( b11, z - dg, x + dg , y , dr11 )
      call ndpp ( b12, z - dg, x - dg , y , dr12 )
      yg1 = y/dg
      yg2 = yg1**2
      yg3 = yg1**3
      yg4 = yg1**4
c
c fudge to fix dipole fringe field -4/8/02  jjl
c
      fudge=1.30986      
      bx = fudge*yg1 * ( (b5-b7)*2./3. - (b6-b8)/12. )  +
     1     yg3*( (b5-b7)/6. - (b6-b8)/12. -
     2     (b3 + b11 - b4 - b12 - 2.*b5 + 2.*b7 ) / 12. )
      by = b0 - yg2*( ( b1 + b9 + b5 + b7 - 4.*b0 ) *2./3. -
     1     ( b2 + b10 + b6 + b8 - 4.*b0 ) / 24. ) +
     2     yg4* (-( b1 + b9 + b5 + b7 - 4.*b0 ) / 6. +
     3     ( b2 + b10 + b6 + b8 - 4.*b0 ) / 24. +
     4     ( b3 + b11 + b4 + b12 - 2.*b1 - 2.*b9 -
     5     2.*b5 - 2.*b7 + 4.*b0 ) / 12. )
      bz = fudge*yg1*( (b1 - b9 ) *2./3. - ( b2 - b10 ) /12. ) +
     1     yg3*( ( b1 - b9 ) / 6. - ( b2 - b10 ) / 12. -
     2     ( b3 + b4 - b11 - b12 - 2.*b1 + 2.*b9 ) / 12.  )
      bt  =dsqrt(bx*bx + by*by + bz*bz)
      return
   14 bx = 0.
      by = br
      bz = 0.
      bt = br
c      save /blck10/,/blck20/,/blck21/,/blck22/,/blck23/,/blck24/,
c     $     /blck25/
      return
      end
c
      subroutine nndip(l,bb)
c****
c****
c**** mtyp = 3 or 4
c**** this version of bfun is mainly for nonuniform field magnets
c**** the central field region is represented to 3'rd order on-and-
c**** off the midplane by analytic expressions. see slac no. 75
c**** fringe field regions represented by fermi type fall-off
c**** along with radial fall-off
c**** components of 'b' in fringe region evaluated by numerical methods
c****
c****
c**** the relationship between b0, ......... b12 and b(i,j) relative to
c**** axes (z,x) is given by
c****
c****
c**** bb(l,0)  = b( 0, 0 )
c**** bb(l,1)  = b( 1, 0 )
c**** bb(l,2)  = b( 2, 0 )
c**** bb(l,3)  = b( 1, 1 )
c**** bb(l,4)  = b( 1,-1 )
c**** bb(l,5)  = b( 0, 1 )
c**** bb(l,6)  = b( 0, 2 )
c**** bb(l,7)  = b( 0,-1 )
c**** bb(l,8)  = b( 0,-2 )
c**** bb(l,9)  = b(-1, 0 )
c**** bb(l,10) = b(-2, 0 )
c**** bb(l,11) = b(-1, 1 )
c**** bb(l,12) = b(-1,-1 )
c****
c****
c modified for use with snake in order to deal with everlapping entrance and
c exit fringing fields - jjl 5/21/91
      implicit real*8(a-h,o-z)
      real*8  ndx, k
      common  /blck10/  bx, by, bz, k, tc, dtc
      common  /blck20/  ndx,bet1,gama,delt,csc
      common  /blck21/  rca,dels,br,s2,s3,s4,s5,s6,s7,s8,scor
      common  /blck22/  d, dg, s, bf, bt
      common  /blck23/  c0, c1, c2, c3, c4, c5
      common  /blck24/  rb, xc, zc
      common  /blck25/  in, mtyp
      dimension tc(2,6), dtc(6),bb(2,0:12)
      x = tc(l,1)
      y = tc(l,2)
      z = tc(l,3)
      dx = x - xc
      dz = z - zc
      rp =dsqrt( dx**2 + dz**2 )
      dr = rp - rb
      go to ( 1, 7, 3, 7 ), in
    7 print 8, in, mtyp
      call exit(0)
    8 format ('0 error -go to -  in nndip   in=', i3, '   mtyp=',i4)
c****
c****
    1 sine = -1.
      go to 5
    3 sine = 1.
    5 if( z  .gt. 0. ) dr = x * sine*csc
      call ndpp( bb(l,0), z, x, y, dr      )
      if( y  .ne. 0. )  go to 6
c****
c**** mid-plane fringing field region
c****
      return
c****
c**** non mid-plane fringing field region
c****
    6 if( z .gt. 0. )  go to 9
      dr1  =       (dsqrt( dx**2 + (dz+dg)**2 ) - rb )
      dr2  =       (dsqrt( dx**2 + (dz+2.*dg)**2 ) - rb )
      dr3  =       (dsqrt( (dx+dg)**2 + (dz+dg)**2 )  - rb )
      dr4  =       (dsqrt( (dx-dg)**2 + (dz+dg)**2 )  - rb )
      dr5  =       (dsqrt( (dx+dg)**2 + dz**2 ) - rb )
      dr6  =       (dsqrt( (dx+ 2.*dg)**2 + dz**2 ) - rb )
      dr7  =       (dsqrt( (dx-dg)**2 + dz**2 ) - rb )
      dr8  =       (dsqrt( (dx- 2.*dg)**2 + dz**2 ) - rb )
      dr9  =       (dsqrt( dx**2 + (dz-dg)**2 ) - rb )
      dr10 =       (dsqrt( dx**2 + (dz-2.*dg)**2 ) - rb )
      dr11 =       (dsqrt( (dx+dg)**2 + (dz-dg)**2 )  - rb )
      dr12 =       (dsqrt( (dx-dg)**2 + (dz-dg)**2 )  - rb )
      go to 10
    9 dr1  = sine* x*csc
      dr2  = dr1
      dr9  = dr1
      dr10 = dr1
      dr3  = sine* ( x + dg )*csc
      dr5  = dr3
      dr11 = dr3
      dr4  = sine*( x - dg )*csc
      dr7  = dr4
      dr12 = dr4
      dr6  = sine* ( x + 2.*dg )*csc
      dr8  = sine* ( x - 2.*dg )*csc
c****
c****
   10 call ndpp ( bb(l,1) , z + dg, x , y , dr1 )
      call ndpp ( bb(l,2) , z + 2.*dg, x , y , dr2 )
      call ndpp ( bb(l,3) , z + dg, x + dg , y , dr3 )
      call ndpp ( bb(l,4) , z + dg, x - dg , y , dr4 )
      call ndpp ( bb(l,5) , z , x + dg , y, dr5 )
      call ndpp ( bb(l,6) , z , x + 2.*dg , y , dr6 )
      call ndpp ( bb(l,7) , z , x - dg , y, dr7 )
      call ndpp ( bb(l,8) , z , x - 2.*dg , y , dr8 )
      call ndpp ( bb(l,9) , z - dg, x , y , dr9 )
      call ndpp ( bb(l,10), z - 2.*dg, x, y, dr10 )
      call ndpp ( bb(l,11), z - dg, x + dg , y , dr11 )
      call ndpp ( bb(l,12), z - dg, x - dg , y , dr12 )
c      save /blck10/,/blck20/,/blck21/,/blck22/,/blck23/,/blck24/,
c     $     /blck25/
      return
      end
c
      subroutine  ndpp ( bfld, z, x, y , dr )
c****
c****
c****
c****
      implicit real*8(a-h,o-z)
      real*8  ndx, k
      common  /blck10/  bx, by, bz, k, tc, dtc
      common  /blck20/  ndx,bet1,gama,delt,csc
      common  /blck21/  rca,dels,br,s2,s3,s4,s5,s6,s7,s8,scor
      common  /blck22/  d, dg, s, bf, bt
      common  /blck23/  c0, c1, c2, c3, c4, c5
      common  /blck24/  rb, xc, zc
      common  /blck25/  in, mtyp
      dimension tc(2,6), dtc(6)
      drr1 = dr/rb
      drr2 = drr1*drr1
      drr3 = drr2*drr1
      drr4 = drr3*drr1
c****
c**** mtyp    :                  modified iterative procedure
c****
      xp = x
      xp2 = xp*xp
      xp3 = xp2*xp
      xp4 = xp3 * xp
      zp = -(s2*xp2 + s3*xp3 + s4*xp4 + s5*xp4*xp + s6*xp4*xp2 +
     1       s7*xp4*xp3 + s8*xp4*xp4 )
      az = (z-zp)/10.d0
      azmax = dsqrt(  x*x + z*z  )
      if( az  .gt.  azmax  )  az = azmax
      zsign = z-zp
      rinv4 = 0.
      do 11 i=1,21
      xp   = x + az*(i-11)
      xp2 = xp*xp
      xp3 = xp2*xp
      xp4 = xp3*xp
      zp = -(s2*xp2 + s3*xp3 + s4*xp4 + s5*xp4*xp + s6*xp4*xp2 +
     1       s7*xp4*xp3 + s8*xp4*xp4 )
      xxp = x-xp
      zzp = z-zp
      dd =            xxp*xxp + zzp*zzp
      if( dd  .lt.  1.d-15 )  dd = 1.d-15
      if( dd  .gt.  1.d15  )  dd = 1.d15
      rinv4 = rinv4 + 1.0d0 / (dd*dd )
   11 continue
      dp = dsqrt( 1.d0/rinv4 )
      dp = dsqrt( dp )
      s = 1.9023d0* dsign( 1.d0, zsign ) * dp/d + dels
c****
c**** first guess for closest point is
c****
c*    xp = x
c*    xp2 = xp*xp
c*    xp3 = xp2*xp
c*    xp4 = xp3*xp
c****
c**** calculate zp on curve for corresponding xp
c****
c*    zp = -( s2*xp2 + s3*xp3 + s4*xp4 + s5*xp4*xp + s6*xp4*xp2 +
c*   1   s7*xp4*xp3 + s8*xp4*xp4 )
c*    zsign = z-zp
c****
c**** slope of curve at xp, zp
c****
c*    do 4 i=1,3
c*    dzdxc = -(2.*s2*xp + 3.*s3*xp2+ 4.*s4*xp3 + 5.*s5*xp4 +
c*   1   6.*s6*xp4*xp + 7.*s7*xp4*xp2 + 8.*s8*xp4*xp3 )
c****
c**** next approximation to closest point is
c****
c*    xp = ( dzdxc*(z-zp)  +  dzdxc*dzdxc*xp + x ) / (1.+dzdxc*dzdxc)
c*    if( i  .eq.  1 )  xp = (3.*xp +  x ) / 4.
c*    xp2 = xp*xp
c*    xp3 = xp2*xp
c*    xp4 = xp3*xp
c*    zp = -( s2*xp2 + s3*xp3 + s4*xp4 + s5*xp4*xp + s6*xp4*xp2 +
c*   1   s7*xp4*xp3 + s8*xp4*xp4 )
c*  4 continue
c*    xxp = x-xp
c*    zzp = z-zp
c*    s = dsign( 1.d0,zsign) * dsqrt( xxp*xxp + zzp*zzp) / d + dels
c****
c****
c****
c****
      cs=c0+s*(c1+s*(c2+s*(c3+s*(c4+s*c5))))
      if( dabs(cs)  .gt.  70.  )  cs =dsign( 70.d0 ,cs  )
      e=dexp(cs)
      p0 = 1.0 + e
      db=bf-br
      bfld = 0.
      if( mtyp .eq. 3 ) bfld =
     1       br +( 1. - ndx*drr1 + bet1*drr2+gama*drr3+delt*drr4)*db/p0
      if( mtyp .eq. 4 ) bfld = br + ( 1./(1. +ndx*drr1) )*db/p0
c****
c**** print 100, x, y, z,  dr, s, bfld
c*100 format( 1p6d15.4 )
c****
c      save /blck10/,/blck20/,/blck21/,/blck22/,/blck23/,/blck24/,
c     $     /blck25/
      return
      end
      subroutine hlmhltz(x,y,z,r,amp,bx,by,bz)
c************************************************************************
c hlmhltz - calculates field for a pair of helmholtz coils
c           x,y,z = field point coords, origin exactly between coils (meters)
c               r = coil radius and separation
c             amp = current in coils
c              bi = x,y,z components of field (tesla)
c 7/8/87
c*************************************************************************
      real x,y,z,r,amp,bx,by,bz,rho,br1,br2,by1,by2,theta,y1,y2
c      write(6,*)"in hlmhltz",r, amp
      y1=(r/2)+y
      y2=y-(r/2)
      rho=sqrt(x**2+z**2)
      if (rho.eq.0.)then
      theta=0.
      go to 1
      endif
      theta=atan2(z,x)
 1    call loop(r,amp,rho,y1,br1,by1)
      call loop(r,amp,rho,y2,br2,by2)
      by=by1+by2
      br=br1+br2
      bx=br*cos(theta)
      bz=br*sin(theta)
      return
      end
      subroutine solenoid(x,y,z,r,rl,amp,bx,by,bz)
c************************************************************************
c solenoid - calculates approximate field for a solenoid (sum over many current loops)
c           x,y,z = field point coords, origin center of coil (meters)
c               r = coil radius
c              rl = length
c             amp = current in coils
c              bi = x,y,z components of field (tesla)
c 1/7/05
c*************************************************************************
      parameter (loops=1000) !even number of loops
      real x,y,z,r,rl,amp,bx,by,bz,rho,br1,br2,by1,by2,theta,yy(loops)
c      write(6,*)"in solenoid ",x,y,z,r,rl,amp
      dl=rl/(loops-1) !space between loops
      loops2 = loops/2
      bx = 0.
      by = 0.
      bz = 0.
      br = 0.
      do 10 i=1,loops
      j=i-1
      yy(i)=(dl*(j-loops2))+y
      rho=sqrt(x**2+z**2)
      if (rho.eq.0.)then
      theta=0.
      go to 1
      endif
      theta=atan2(z,x)
 1    call loop(r,amp,rho,yy(i),br1,by1)
      by=by+by1
      bx=bx+br1*cos(theta)
      bz=Bz+br1*sin(theta)
   10 continue
c      write(6,*)bx,by,bz  
      return
      end
      
      subroutine loop(a,amp,rho,y,br,by)
c************************************************************************
c  loop - calculates the magnetic field for a current loop
c         a = radius of loop
c       amp = current in loop
c       rho = distance from loop axis (meters)
c         y = distance from loop along axis (meters)
c        br = radial component of field (tesla)
c        by = axial component of field
c
c  reference: garrett, journal of applied physics, 34 sept 1963, p. 2567
c  7/8/87
c******************************************************************
      real a,amp,rho,y,br,by,r1sq,r2sq,ksq,kpsq
c	real vk,ve
      real*8 dksq,dvk,dve
c evaluate parameters
      r1sq=(a+rho)**2+y**2
      r2sq=(a-rho)**2+y**2
      ksq=4.*a*rho/r1sq
      kpsq=1-ksq
c evaluate elliptical integrals
      dksq=ksq
      call fb01ad(dksq,dvk,dve)
c      ve=dve
c      vk=dvk
c      write(6,*)ksq,vk,ve
c calculate fields
      if(r2sq.eq.0.)then
      write(6,*)' bull''s-eye!!! you just hit a coil.'
      by=1.69e+38
      go to 15
      endif
      by=(amp*2.0e-07/sqrt(r1sq))*((dvk-dve)+(a-rho)*2*a*dve/r2sq)
 15   if(rho.eq.0.)then
      br=0.
      return
      endif
      if(kpsq.eq.0.)kpsq=0.3e-36
      br=-(amp*1.e-07*y)/(rho*sqrt(r1sq)*kpsq)*((2-ksq)*(dvk-dve)
     &     -ksq*dvk)
      return
      end      
      subroutine fb01ad(c,  vk,ve)
      implicit real*8(a-h,o-z)
c*ibm real*8 xlg/  z7fffffffffffffff /
      real * 8 xlg/'7ffffffff'x/
      d=1d0-c
      if(d .gt. 0d0)e=-log(d)
c**** harwell version of fb01ad
      if(c .ge. 1d0)go to 2
           ve=e*((((((((((
     a     3.18591956555015718d-5*d  +.989833284622538479d-3)*d
     b    +.643214658643830177d-2)*d +.16804023346363385d-1)*d
     c    +.261450147003138789d-1)*d +.334789436657616262d-1)*d
     d    +.427178905473830956d-1)*d +.585936612555314917d-1)*d
     e    +.937499997212031407d-1)*d +.249999999999901772d0)*d)
     f    +(((((((((
     g     .149466217571813268d-3*d  +.246850333046072273d-2)*d
     h    +.863844217360407443d-2)*d+.107706350398664555d-1)*d
     i    +.782040406095955417d-2)*d +.759509342255943228d-2)*d
     j    +.115695957452954022d-1)*d +.218318116761304816d-1)*d
     k    +.568051945675591566d-1)*d +.443147180560889526d0)*d
     l    +1d0
c****
c**** routine modified to calculate vk and ve always
c****
c****
           vk=e*((((((((((
     a     .297002809665556121d-4*d   +.921554634963249846d-3)*d
     b    +.597390429915542916d-2)*d  +.155309416319772039d-1)*d
     c    +.239319133231107901d-1)*d  +.301248490128989303d-1)*d
     d    +.373777397586236041d-1)*d  +.48828041906862398d-1)*d
     e    +.703124997390383521d-1)*d  +.124999999999908081d0)*d
     f    +.5d0)+(((((((((
     g     .139308785700664673d-3*d   +.229663489839695869d-2)*d
     h    +.800300398064998537d-2)*d  +.984892932217689377d-2)*d
     i    +.684790928262450512d-2)*d  +.617962744605331761d-2)*d
     j    +.878980187455506468d-2)*d  +.149380135326871652d-1)*d
     k    +.308851462713051899d-1)*d  +.965735902808562554d-1)*d
     l    +1.38629436111989062d0
      return
    2 ve=1d0
      vk=xlg
      return
      end
      subroutine dfringe(x,y,z,fact,a,fx,fy,fz)
      real x,y,z,fx,fy,fz,fact,a,delta
      real c0(0:5),s,b(-2:2,-2:2)
      data c0/0.2383,1.7395,-0.4768,0.5288,-0.1299,0.0222/
      data delta/80./
c
c   calculates  entrance fringe field of dipole magnet
c   in the same manner as raytrace, for use in snake. jjl 2/17/87
c   n.b. can only be used for bending in the x-y plane.
c
c          fact = central field 
c             a = gap size (radius)
c         x,y,z = coordinates relative to magnet entrance
c                  (set entrance=0 in the relative ref frame)
c            fz = returned field value
c         c0(i) = fringe field coefficients
c     
c            fz = fact/(1+exp(s))
c             s = c0(0) + c0(1)*(y/a) + c0(2)*((y/a)**2) ....etc.
c****************************************************************
c
c  set up grid of b's for expansion (see raytrace manual p. 11-12)
c
      do 20 j=-2,2
      do 20 k=-2,2
      if(abs(j)+abs(k).ge.3) go to 20
      s=c0(0)
      do 10 i=1,5
c
c  when curved efb's are introduced y will depend on x in the following
c
 10    s=s+(c0(i)*(((y+(j*delta))/a)**i))
 20   b(j,k)=fact/(1+exp(s))
c
c  calculate fields
c        1         2         3         4         5         6         7
      fz=b(0,0)
     &   -((z**2/delta**2)*
     &      (((2./3.)*(b(1,0)+b(-1,0)+b(0,1)+b(0,-1)-(4.*b(0,0))))
     &     -((1./24.)*(b(2,0)+b(-2,0)+b(0,2)+b(0,-2)-(4.*b(0,0)))))
     &   +((z**4/delta**4)*
     &      (((-1./6.)*(b(1,0)+b(-1,0)+b(0,1)+b(0,-1)-(4.*b(0,0))))
     &      +((1./24.)*(b(2,0)+b(-2,0)+b(0,2)+b(0,-2)-(4.*b(0,0))))
     &      +((1./12.)*(b(1,1)+b(-1,1)+b(1,-1)+b(-1,-1)))
     &      -((1./6.0)*(b(1,0)+b(-1,0)+b(0,1)+b(0,-1)-(2.*b(0,0)))))))
c        1         2         3         4         5         6         7
      fy=((z/delta)*
     &      (((2./3.)*(b(1,0)-b(-1,0)))
     &       -((1./12.)*(b(2,0)-b(-2,0)))))
     &   +(((z**3)/(delta**3))*
     &      (((1./6.)*(b(1,0)-b(-1,0)))
     &       -((1./12.)*(b(2,0)-b(-2,0)))
     &       -((1./12.)*(b(1,1)+b(1,-1)-b(-1,1)-b(-1,-1)
     &           -(2.*b(1,0))+(2.*b(-1,0))))))
c        1         2         3         4         5         6         7
      fx=((z/delta)*
     &      (((2./3.)*(b(0,1)-b(0,-1)))
     &       -((1./12.)*(b(0,2)-b(0,-2)))))
     &   +(((z**3)/(delta**3))*
     &      (((1./6.)*(b(0,1)-b(0,-1)))
     &       -((1./12.)*(b(0,2)-b(0,-2)))
     &       -((1./12.)*(b(1,1)+b(-1,1)-b(1,-1)-b(-1,-1)
     &           -(2.*b(0,1))+(2.*b(0,-1))))))
      return
      end

       subroutine dipole(x,y,b0)
c*************************************************************************
c                                                                        *
c                       d i p o l e                                      *
c                                                                        *
c*************************************************************************
c   for meth=3 , type=2 , index=16+17 (arbitrary field numeric.deriv.)
        common/dipff/xtra,ytra,atra,dref,rms,alps,bets
     @           ,xc(2),yc(2),r0(2),e0(2),s0(2),s1(2),s2(2),s3(2)
     @           ,tbound,dst
c
c                field with gradient (alpha and beta)
c                        (x=0.,y=0.,z=0.)=center of symmetry
c                         z=0.:symmetry plane
c                and curvated faces
c      this routine compute the field only in the symm. plane
c      the fringing field of the faces is computed
c      out of plane field needs numerical derivation.
c
c       rectangular-->cylindrical coord.:
101     r=sqrt(x*x+y*y)
        rn=r/rms
        te=atan2(y,x)
        dr=rn-1.
c       z component of the field in sym.pl.:
        b0=1.-alps*dr+bets*dr*dr
c               input or exit fringing field:
        if(te.ge.tbound)then
                i=2
        else
                i=1
        endif
                ar0=abs(r0(i))
                sign=r0(i)/ar0
                er=e0(i)/b0
                rr=sqrt((x-xc(i))**2+(y-yc(i))**2)
                drf=(rr-ar0)*sign/er
c               if(i.eq.2.and.drf.gt.-1.)then
c                 i=2
c               endif
c               if(drf.le.0.)goto1
                dr2=drf*drf
                sr=s0(i)+s1(i)*drf+s2(i)*dr2+s3(i)*dr2*drf
                ch=1./(1.+exp(sr))
                b0=b0*ch
c 1    save /dipff/       
  1    return
        end
        subroutine geo(xyz,eject)
c*************************************************************************
c                                                                        *
c                       g e o                                            *
c                                                                        *
c*************************************************************************
c
c       eject=.true. if the point xyz is not in a free space.
c        coordinates are relative to the frame of the region indre.
c
        include 'snake.inc'
c                         declarations for tor : begening
      parameter (dpi=2.*pi)
      parameter (ngalette=8)
c           % nombre de galettes
      parameter (alphs2=pi/ngalette)
      parameter (alpha=2.*alphs2)
c           %angle entre les plans de 2 galet. voisines
      parameter (r1=200.)
c           %dist. mini. axe du tore--->centre du conducteur
      parameter (r2=1200.)
c           %dist. maxi               "
      parameter (r3=(r2-r1)/2.)
c           %petit rayon du tore
      parameter (dr=50.)
c           % demi-extension du conducteur le long du rayon
      parameter (dz=10.)
c           % demi-extension du conducteur perpend. au rayon
      parameter (rmoy=(r1+r2)/2.)
c           %grand rayon du tore
      parameter (margr=50.)
      parameter (margz=50.)
      parameter (drt=dr+margr)
      parameter (dzt=dz+margz)
      parameter (rmaxt=r3+drt)
      parameter (rmaxt2=rmaxt**2)
c                         declarations for tor : end
c
        logical eject,first
        dimension xyz(3)
c                next line : begening of the /diag/ package
	character diag*100,diagt*100
	common/cdiag/diag(maxve),diagt
	common/ndiag/iepdead(maxve),iredead(maxve)
c                previous line : end of the /diag/ package
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
c                         declarations for tordu : begening
      common/geoh/zcol,xcol,rcol,yaa,yaaa,zlim,zmin,ay2max,zmax
      common/geod/xtab(11),ytab(11),ztab(11),ttheta,zv,zh
     s ,xa,ya,za,xb,yb,zb,xc,yc,zc,ya1,za1,yc1,zc1,yc2,zc2
c                         declarations for tordu : end
c
        data first/.true./
c
      if(rname(indre).eq.'tor')goto2000
      eject=.false.
c      save /cdiag/,/ndiag/,/crelfr/,/relfra/,/geoh/,/geod/
      return
c
c                         case of the tore :
c
2000     if(first)then
                 ta=tan(alphs2)
                 cas2=cos(alphs2)
                 sas2=sin(alphs2)
                 first=.false.
        endif
c        on passe dans les coordonnees de la region 2:
        eject=.false.
                 ax=abs(xyz(1))
                 y=xyz(2)
                 az=abs(xyz(3))
        if(ax.eq.0.)goto10
        eject=az/ax.gt.ta
        if(eject)then
c		built the diag:
		write(diagt,'(2a,2(e10.3,a))')
     s		' try to cross the'
     s           , ' plane of a coil ; angle ='
     s		,az/ax,' max =',ta,' (rad)'
c      save /cdiag/,/ndiag/,/crelfr/,/relfra/,/geoh/,/geod/
      return
      endif
 10   continue
c
c        region du champ central du tore : encombrement des 2 galettes
c       superieures , dans les coordonnees relatives a cette region :
c
c        distance au plan de la galette:
        z1=az*cas2-ax*sas2
c           %dans le rep. ou z=axe galette
        x1=ax*cas2+az*sas2
	if((x1+2.67*y+4200.).lt.0.)return
	if((x1-.49*y-620.).lt.0.)return
        if(abs(z1).gt.dzt)return
c  3 next lines : temporary for infinite pancake 
	eject=.true.
	az1=abs(z1)
c		built the diag:
	write(diagt,'(2a,2(e10.3,a))')
     s		'try to go inside the cryostat'
     s           , ' ; dist./coil plane ='
     s		,az1,' min =',dzt
c      save /cdiag/,/ndiag/,/crelfr/,/relfra/,/geoh/,/geod/
      return
c		end of temporary
   20     x1=ax*cas2+az*sas2-rmoy
c        distance**2 a l'axe de la galette:
        ray2=y*y+x1*x1
        eject=ray2.lt.rmaxt2
      if(eject)write(diagt,'(2a,2(e10.3,a))')
     s		'try to go inside the cyl.'
     s           , ';dist./coil plane ='
     s		,az1,' min =',dzt
c      save /cdiag/,/ndiag/,/crelfr/,/relfra/,/geoh/,/geod/
        return
        end
c
c*********************************************************************
c                                                                    *
c         a p h a s a , r a n f , r a n g e t  and  r a n s e t      *
c                                                                    *
c*********************************************************************
       function aphasa(a1,a2)
       common /rand/ix
c vax:       (ran is a vax system subroutine)
c vax:       aphasa=a1+ran(ix)*(a2-a1)
c sun:
        aphasa=a1+ranf(ibidon)*(a2-a1)
c
       return
       end
c
       function ranf(ibidon)
       common /rand/ix
       equivalence (r,ir)
       iz=ix*899
c       if(iz)8,9,9
c 8     iz=iz+2147483647+1
c 9     r=iz
      if(iz.le.0)then
	  iz=iz+2147483647+1
	  endif
	  r=iz
	  
       ranf=r/2147483647.	
       ix=iz
c      save /rand/
       return
       end
c
       subroutine ranget(irand)
       common /rand/ix
       irand=ix
c      save /rand/
       return
       end
c
       subroutine ranset(irand)
       common /rand/ix
       ix=irand
c      save /rand/
       return
       end
c
c  end of snake:f
