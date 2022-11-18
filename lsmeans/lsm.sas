data lsm;
   input A $ cov dg B;
   datalines;
   a 350  970  1
   a 400 1000  1
   a 360  980  1
   a 350  980  2
   a 340  970  2
   b 390  990  1
   b 340  950  1
   b 410  980  2
   b 430  990  2
   b 390  980  2
   c 400  990  1
   c 320  940  1
   c 330  930  2
   c 390 1000  2
   c 420 1000  2
;

proc glm data=lsm;
   class a b;
   model dg = A B A*B cov;
   lsmeans A B A*B / stderr cl;
run;
