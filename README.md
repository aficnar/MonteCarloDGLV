# MonteCarloDGLV
Monte Carlo code for calculating the DGLV radiative energy loss at arbitrary orders in opacity.

* At the beginning of my PhD, together with A. Buzzatti, I have developed a flexible and efficient Fortran code 
for numerical evaluation of the induced gluon distribution in the [DGLV model](http://arxiv.org/abs/nucl-th/0310076) 
at an arbitrary order in opacity. This distribution tells us the details of the energy loss that
a high energy particle traveling through a quark-gluon plasma (created in heavy ion collisions at RHIC and LHC) suffers. 
For more information, see Section 2.3 in my [PhD thesis](http://academiccommons.columbia.edu/catalog/ac%3A178237).

* At the heart of our code is the importance sampling Monte Carlo integration algorithm,
which is built from first principles. In general, importance sampling is one of the best
variance reduction techniques for estimating an integral using the Monte Carlo integration
method. 

* The relevant formula for gluon distribution is notoriously complicated (a 3n-dimensional
integral at opacity order n), but we were able to construct a relatively fast automated numerical integration code, 
that takes about 3^n seconds at a given order, on an average personal computer and with average settings.

* This code was one of first steps in the development of the [CUJET model](http://arxiv.org/abs/1106.3061), 
a state-of-the-art implementation of the DGLV model.
