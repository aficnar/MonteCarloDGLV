C23456789012345678901234567890123456789012345678901234567890123456789012

	SUBROUTINE INTEGRATION(n,nn,E,x,k,mu,mg,Mq,L,lg,CR,as,mcarlo,mcarlo2,idum)
	implicit none
	integer n,nn,delta(4096,0:12),l1,u1,l2,u2,u,v,w,i,j,m,s,sign,
     *		factorial,idum
	double precision E,x,k,mu,mg,Mq,L,Le,lg,CR,as,qq(0:12),phi(0:12),
     *			     dotmatrix(0:12,0:12),dotarray(2,13),dot,
     *			     cc,icsq,icsz,icsf,factor,partial,total,mcarlo,
     *				 Ex2inv,beta2,c1,c2,pi,hc,mcarlo2
	complex complexarray(13),phase,ph1,ph2
	common /array/ delta
	data qq,phi/13*0.d0,13*0.d0/
	parameter (pi=3.1415927d0,hc=0.197d0)

C	*** Precomputation of various factors
	Ex2inv=1./(2*E*x*hc)
	Le=L/(n+1)
	beta2=mg*mg*(1-x)+Mq*Mq*x*x
	factor=2*CR*as*(L/lg)**n/(pi*pi)/factorial(n)/nn

C	*** The first element in the arrays corresponds to k
	qq(0)=k
	phi(0)=0.0d0

C	*** All the possible combinations of 0 and 1 are generated, with the corresponding sign in position '0'
	call srdelta(delta,n)

C	*** This is the main loop that does the Monte Carlo evaluation, nn are the sampling points, n is the order
	mcarlo=0.d0
	mcarlo2=0.d0
	do 10 w=1,nn
C	  *** Generate the random variables and store them into correspoding arrays
	  do 11 i=1,n
	    qq(i)=icsq(mu,idum)
	    phi(i)=icsf(idum)
11	  continue
C	  *** All the possible scalar products of q's are computed and stored in the matrix dotmatrix
	  do 12 i=0,n
	    dotmatrix(i,i)=qq(i)*qq(i)
	    do 112 j=i+1,n
	      dotmatrix(i,j)=qq(i)*qq(j)*dcos(phi(i)-phi(j))
	      dotmatrix(j,i)=dotmatrix(i,j)
112	    continue
12	  continue
C	  *** Here starts the loop over all the terms from prod_i(v^2(q_i)-delta(q_i))
	  total=0.0d0
	  do 13 s=1,(2**n-1)
	    sign=delta(s,0)
C		*** All the relevant dot products are precomputed and stored in dotarray
	    do 113 i=1,n+1
	      dotarray(1,i)=dot(1,n,i,n,dotmatrix,s)
	      dotarray(2,i)=dot(i,n,i,n,dotmatrix,s)+beta2
	      complexarray(i)=cmplx(1.d0, -dotarray(2,i)*Ex2inv*Le)
113	    continue
C		*** Here starts the main loop for the summation in the formula
	    partial=0.0d0
	    c1=cc(1,dotarray)
	    do 313 m=1,n
	      c2=cc(m+1,dotarray)
	      ph2= phase(m,dotarray,complexarray,Ex2inv,Le)
	      ph1= ph2* complexarray(1)/
     *       (1+(dotarray(2,1)*Ex2inv*dotarray(2,1)*Ex2inv*Le*Le))
	      partial=partial+(c1-c2)*(REAL(ph2)-REAL(ph1))
	      c1=c2
313	    continue
	    total=total+sign*partial
13	  continue
	  mcarlo=mcarlo+total
	  mcarlo2=mcarlo2+total*total
10	continue
	mcarlo=mcarlo*factor
	mcarlo2=mcarlo2*factor*factor*nn
	return
	end

C	--------------------------------------------------------------------------------------------------------------

C	*** This is the function that contructs the dot products (k - q_l1 - ... - q_u1) .dot. (k - q_l2 - ... - q_u2)
	double precision function dot(l1,u1,l2,u2,dotmatrix,s)
	  double precision dotmatrix(0:12,0:12)
	  integer l1,u1,l2,u2,i,j,s
	  integer delta(4096,0:12)
	  common /array/ delta
	  dot=dotmatrix(0,0)
	  do 50 i=l1,u1
	    dot=dot-dotmatrix(0,i)*delta(s,i)
50	  continue
	  do 51 j=l2,u2
	    dot=dot-dotmatrix(0,j)*delta(s,j)
51	  continue
	  do 52 i=l1,u1
	    do 152 j=l2,u2
	        dot=dot+dotmatrix(i,j)*delta(s,i)*delta(s,j)
152	    continue
52	  continue
	return
	end

C	--------------------------------------------------------------------------------------------------------------

C	*** Here we construct the terms cc(m,...) = C(1,...,n) .dot. C(m,...,n)
	double precision function cc(m,dotarray)
	  integer m
	  double precision dotarray(2,13)
	  cc=dotarray(1,m)/(dotarray(2,1)*dotarray(2,m))
	return
	end

C	--------------------------------------------------------------------------------------------------------------

C	*** This is the term inside the COS functions
	complex function phase(m,dotarray,complexarray,Ex2inv,Le)
	  integer m,i
	  double precision dotarray(2,13),Ex2inv,Le,temp
	  complex complexarray(13)
	  phase=(1.d0,0.d0)
	  do 53 i=2,m
	    temp= 1+(dotarray(2,i)*Ex2inv*dotarray(2,i)*Ex2inv*Le*Le)
	    phase=phase*complexarray(i)/temp
53	  continue
	return
	end

C	--------------------------------------------------------------------------------------------------------------

C	*** All the possible combinations of 0 and 1 are stored in the array delta. The sign is stored at position 0
	subroutine srdelta(delta,n)
	  integer delta(4096,0:12),n
	  integer switch,sign,i,j
	  do 70 i=1,2**n
	    switch=i
	    sign=0
	    do 170 j=1,n
	      delta(i,j)=switch-(switch/2)*2
	      switch=switch/2
	      sign=sign+delta(i,j)
170	    continue
	  delta(i,0)=(-1)**(n-sign)
70	  continue
	return
	end