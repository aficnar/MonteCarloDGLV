C23456789012345678901234567890123456789012345678901234567890123456789012
	program frontend
	implicit none
	integer iseed,idum,n,nn,log,i,j,contseed,contplot,contkmax,
     *		kpoints,contmg,timer1,timer2
	real rand
	double precision E,x,k,mu,mg,Mq,L,lg,CR,as,mcarlo,p,dN0,xmin,
     *				 xmax,xstep,kmin,kmax,kstep,oldmcarlo(1000),pi,hc,
     *				 cum(1000),cumsc(1000),karray(1000),mcarlo2,
     *				 mcarloerr(1000),oldmcarloerr(1000),errortemp,
     *			     mcarloerrsc(1000)
	character cx*4,cj*2,ce*5,name*48,date*6,datef*24,textseed*14,
     *		  textplot*10,textkmax*14,textmg*14,output*8
	parameter (pi=3.1415927d0,hc=0.197d0)

	call fdate(datef)
	date = datef(5:10)

C	*** Reading input from the input file	
	open(unit=3,file='input.dat',status='old')
	read(3,*)
	read(3,*)
	read(3,*)
	read(3,*)
	read(3,*)
	read(3,*)
	read(3,*)
	read(3,*) x,E,mu,Mq,mg,L,lg,CR,as
	read(3,*)
	read(3,*)
	read(3,*)
	read(3,*)
	read(3,*) n,nn,iseed,kmin,kmax,kpoints
	read(3,*)
	read(3,*)
	read(3,*)
	read(3,*)
	read(3,*) contseed,contplot,contkmax,contmg
	close(3)
C	*** Utilizing control options from the input file
	if (contseed.EQ.1) iseed=TIME()
	if (contkmax.EQ.1) kmax=2*E*x*(1-x)
	if (contkmax.EQ.2) then
	  if (x .LT. 0.5) then
	    kmax=2*E*x
		textkmax='     2Ex      '
	  else
	    kmax=2*E*(1-x)
		textkmax='   2E(1-x)    '
	  endif
	endif
	if (contmg.EQ.1) then
	  Mq=0.0
	  mg=0.0
	endif
	if (contmg.EQ.2) then
	  Mq=0.25
	  mg=0.354
	endif
	if (contseed .EQ. 0) textseed='  From input  '
	if (contseed .EQ. 1) textseed='    TIME()    '
	if (contplot .EQ. 0) textplot='  Linear  '
	if (contplot .EQ. 1) textplot='   Log    '
	if (contkmax .EQ. 0) textkmax='  From input  '
	if (contkmax .EQ. 1) textkmax='   2Ex(1-x)   '
	if (contmg .EQ. 0) textmg='  From input  '
	if (contmg .EQ. 1) textmg='   Massless   '
	if (contmg .EQ. 2) textmg=' Thermal Mass '

	idum=-iseed
	p=sqrt(E*E-Mq*Mq)
	do 03 i=1,kpoints
	  oldmcarlo(i)=0.D0
	  oldmcarloerr(i)=0.D0
03	continue

C	*** Loop 02 goes over all orders to be computed
	do 02 j=1,n
C	*** Create the name of the output file for order j and open it
	  write(cx,901) x
901	  format (F4.2)
	  write(cj,902) j
902	  format (I2)
	  if (j.LT.10) cj(1:1)='0'
	  write(ce,903) E
903	  format (F5.1)
	  if (E.LT.100.D0) ce(1:1)='0'
	  name='output_x='//cx//'_E='//ce//'_N='//cj//'_'//date//'_Z.dat'
	  open(unit=4,file=name,status='unknown')
C	*** Write the preamble of the output file
	  write(4,*) '--------------------------------'
	  write(4,*) '----- DGLV Output File Z0.25 ---'
	  write(4,*) '--------------------------------'
	  write(4,*) datef
	  write(4,*) ''
	  write(4,*) 'Physical parameters used:'
	  write(4,*) ' ----------------------------------------------------------------------------------------'
	  write(4,*) '|    x    |    E	  |   mu    |    M     |   mg    |   L    | lambda |   C_R   | alpha_s |'
	  write(4,*) '|---------|---------|---------|----------|---------|--------|--------|---------|---------|'
	  write(4,904) '|',x,'|',E,'|',mu,'|',Mq,'|',mg,'|',L,'|',lg,'|',CR,'|',as,'|'
904	  format(1X,A1,2X,F5.3,2X,A1,2X,F5.1,2X,A1,2X,F5.3,2X,A1,2X,F6.3,2X,A1,2X,F5.3,2X,A1,2X,F4.1,2X,A1,2X,F4.1,2X,A1,2X,F5.3,2X,A1,2X,F5.3,2X,A1)
	  write(4,*) ' ----------------------------------------------------------------------------------------'
	  write(4,*) ''
	  write(4,*) 'Other parameters used:'
	  write(4,*) ' -----------------------------------------------------------------'
	  write(4,*) '| Order | Sampl. pts |     Seed     |  Min. k  |  Max. k  | k-pts |'
	  write(4,*) '|-------|------------|--------------|----------|----------|-------|'
	  write(4,905) '|',j,'|',nn,'|',iseed,'|',kmin,'|',kmax,'|',kpoints,'|'
905	  format(1X,A1,1X,I3,3X,A1,1X,I8,3X,A1,2X,I10,2X,A1,1X,F6.2,3X,A1,2X,F6.2,2X,A1,2X,I3,2X,A1)
	  write(4,*) ' -----------------------------------------------------------------'
	  write(4,*) ''
	  write(4,*) 'Controls used:'
	  write(4,*) ' -------------------------------------------------------'
	  write(4,*) '|     Seed     | Distrib. |    Max. k    |     Mass     |'
	  write(4,*) '|--------------|----------|--------------|--------------|'
	  write(4,*) '|',textseed,'|',textplot,'|',textkmax,'|',textmg,'|'
	  write(4,*) ' -------------------------------------------------------'
C	*** Loop 01 goes over all k points to be computed and records the results
	  timer1=time()
	  do 01 i=1,kpoints
	    if (contplot .EQ. 0) k=kmin+(i-1)*(kmax-kmin)/(kpoints-1)
	    if (contplot .EQ. 1) k=kmin*(kmax/kmin)**(float(i-1)/float(kpoints-1))
	    call INTEGRATION(j,nn,E,x,k,mu,mg,Mq,L,lg,CR,as,mcarlo,mcarlo2,idum)
		karray(i)=k
		cum(i)=mcarlo+oldmcarlo(i)
		oldmcarlo(i)=cum(i)
	    dN0=CR*as*(1-x+x*x/2)/(pi*pi*k*k)
		cumsc(i)=cum(i)/dN0
	    errortemp=sqrt((mcarlo2-mcarlo*mcarlo)/nn)
	    mcarloerr(i)=sqrt(oldmcarloerr(i)**2+errortemp**2)
	    oldmcarloerr(i)=mcarloerr(i)
	    mcarloerrsc(i)=mcarloerr(i)/dN0
01	  continue	
	  timer2=time()
	  call timer(timer2-timer1,output)
	  write(4,*) ''
	  write(4,*) 'Time it took for this order: ',output
	  write(4,*) ''
	  write(4,*) 'Output: k, dN(n), Error[dN(n)], dN(n)/dN(0), Error[dN(n)/dN(0)]'
	  write(4,*) '---------------------------------------------------------------'
C	*** Loop 05 writes the output to the output file
	  do 05 i=1,kpoints
	    write(4,*) karray(i),cum(i),mcarloerr(i),cumsc(i),mcarloerrsc(i)    
05	  continue
	  close(4)
02	continue

	stop
	end



C ----------------------------------------------------------------------------------------------------------
C	*** Converts interval interv in seconds and gives output string in form of hh:mm:ss
	subroutine timer(interv,output)
	integer interv,min,sec,hour
	character chour*2,cmin*2,csec*2,output*8
      hour=interv/3600
	min=(interv-hour*3600)/60
	sec=interv-hour*3600-min*60
	write(chour,04) hour
	write(cmin,04) min
	write(csec,04) sec
04	format (I2)
	if (hour.LT.10) chour(1:1)='0'
	if (hour.LT.1) chour(2:2)='0'
	if (min.LT.10) cmin(1:1)='0'
	if (min.LT.1) cmin(2:2)='0'
	if (sec.LT.10) csec(1:1)='0'
	if (sec.LT.1) csec(2:2)='0'
	output=chour//':'//cmin//':'//csec
	return
	end
