C23456789012345678901234567890123456789012345678901234567890123456789012

C	*** Inverse cumulative distribution function of the Z variable (exponential)
	double precision function icsz(L,n,idum)
	  integer n,idum
	  double precision L
	  icsz=-L/(n+1.)*log(1.-ran2(idum))
	return
	end

C	*** Inverse cumulative distribution function of the Q variable (yukawa)
	double precision function icsq(mu,idum)
	  integer idum
	  double precision mu
	  double precision temp1
	  temp1=ran2(idum)
	  icsq=mu*sqrt(temp1/(1.-temp1))
	return
	end

C	*** Inverse cumulative distribution function of the phi variable (flat)
	double precision function icsf(idum)
	  integer idum
	  icsf=6.2831853d0*ran2(IDUM)
	return
	end

C	*** Factorial function (used to compute the factor)
	integer function factorial(a)
	  integer temp2,a,i
	  temp2=1
	  do 71 i=2,a
	    temp2=temp2*i
71	  continue
	  factorial=temp2
	return
	end

C	*** RAN2
	function ran2(idum)
	  integer idum,IM1,IM2,INM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
	  real ran2,AM,EPS,RNMX
	  parameter(IM1=2147483563,IM2=2147483399,AM=1.d0/IM1,INM1=IM1-1,
     *   IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *   NTAB=32,NDIV=1+INM1/NTAB,EPS=1.2d-7,RNMX=1.-EPS)
C     *** Long period (>2x10^{18}) random number generator of L'Ecuyer
C     *** with Bays-Durham shuffle and added safeguards. Returns a uniform
C     *** random deviate between 0.0 and 1.0 (exclusive of the endpoint
C     *** values). Call with idum a negative integer to initialise;
C     *** thereafter, do not alter idum between successive deviate in a
C     *** sequence. RNMX should approximate the largest floating value
C     *** that is less than 1.
	  integer idum2,j,k,iv(NTAB),iy
	  save iv,iy,idum2
	  data idum2/123456789/,iv/NTAB*0/,iy/0/

	  if(idum .le. 0) then
	    idum = max(-idum,1)
	    idum2 = idum
	    do 72, j=NTAB+8,1,-1
	      k=idum/IQ1
	      idum = IA1*(idum-k*IQ1)-k*IR1
	      if(idum .lt. 0) idum = idum + IM1
	      if(j .le. NTAB) iv(j) = idum
72	    continue
	    iy = iv(1)
	  end if
	  k = idum/IQ1
	  idum = IA1*(idum-k*IQ1)-k*IR1
	  if(idum .lt. 0) idum = idum + IM1
	  k = idum2/IQ2
	  idum2 = IA2*(idum2-k*IQ2)-k*IR2
	  if(idum2 .lt. 0) idum2 = idum2 + IM2
	  j = 1 + iy/NDIV
	  iy = iv(j) - idum2
	  iv(j) = idum
	  if(iy .lt. 1) iy = iy + INM1
	  ran2 = min(AM*iy,RNMX)
	return
	end
