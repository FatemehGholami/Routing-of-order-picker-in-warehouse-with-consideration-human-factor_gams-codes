$title Improved Version of eps-Constraint Method for Multiobjective Optimization (EPSCMMIP,SEQ=384)

$onText
The Augmented Epsilon Constraint Method version 2 (AUGMECON2)
The method is applied to a Multi-Objective Integer Programming problem
(specifically a Multi-Objective Multi-Dimensional Knapsack Problem)
with 50 binary variables X, 2 objective functions and 2 constraints

The AUGMECON2 can be used to generate the exact Pareto set (all the Pareto
optimal solutions) if the step size (i.e.the interval between the grid points
of the objective functions that are used as constraints) is appropriately
chosen. For problems with integer objective function coeffcients the step size
should be at most equal to unity.

The exact Pareto set of the specific problem consists of 35 Pareto optimal
solutions. The solution time is approximately 4 secs on a 2012 machine. The
gridpoints are set to 491 = the second objective function range. The output file
2kp50_augmecon2_results.txt contains the payoff table, the gridpoints and the
Pareto optimal solutions. The indication "jump" is used to flag when one or
more grid points are skipped.

The model is separable. The first part of the model (till line 129) is the
problem description and the second part (from line 130) is the implementation
of the method

Additional information can be found at:

http://www.gams.com/modlib/adddocs/epscmmip.pdf


Mavrotas, G, Effective implementation of the eps-constraint method in
Multi-Objective Mathematical Programming problems.
Applied Mathematics and Computation 213, 2 (2009), 455-465.

Mavrotas, G, and Florios, K, An improved version of the augmented
eps-constraint method (AUGMECON2) for finding the exact Pareto set in
Multi-Objective Integer Programming problems.
Applied Mathematics and Computation 219, 18 (2013), 9652-9669

Keywords: mixed integer linear programming, multiobjective optimization,
          eps-constraint method, exact Pareto set, mathematics
$offText

$eolCom //

$sTitle Example Model Definitions
Set
  i /0*6/
  i2(i) index of products(vertex) /1*6/
  K 'objective functions' / k1* k2 /;

   alias(i,j);
   alias(i2,j2);

parameters
   dir(k) 'direction of the objective functions 1 for max and -1 for min' / k1 -1, k2 -1 /
   w(i) weighte of product i   /0 8,1 22,2 18,3 14,4 14,5 4,6 10/
   ww /-0.1/
   ;

scalar
   BW operator body weight /90/
   a coefficient friction of the surface /0.8/
   s sex /1/
   vmax maximum speed /0.9/
   v1 /0.008/
   v2 /0.4/
   M /9999999999999999/;


table d(i,j) Distance Matrix

      0   1    2    3    4    5    6

   0  0   5    15   12   19   44   35
   1  5   0    8    7    14   39   30
   2  15  8    0    17   24   29   40
   3  12  7    17   0    15   32   29
   4  19  14   24   15   0    25   14
   5  44  39   29   32   25   0    23
   6  35  30   40   29   14   23   0
   ;

variables

  Z(K) 'objective function variables'
   ;

positive variable
   u(i)
   t(i,j)
   deltaE(i,j)
   Eprim(i,j)
   E(i,j)
   wprim(j)
   b(i,j)
   v(i,j)

   ;
binary variable
   x(i,j);

equations
  ObjectiveFunction1
  ObjectiveFunction2
  co1(i2)
  co2(j2)
  co3
  co4
  co5(i,j)
  co6(j2)
*co7(i,j)
 co8(i,j)
 co9(i,j)
*co99(i,j)
 co10(i,j)
 co11(i,j)
*co19(i,j)
*co20(i,j)
*co21(i,j)

 co22(i,j)
 co23(i,j)

 co31(i)
 co41(i,j)
 co61
 ;

ObjectiveFunction1                                         ..  Z('k1')=e=sum((i,j)$ (ord(i)<>ord(j)),t(i,j)*x(i,j));
ObjectiveFunction2                                         ..  Z('k2')=e=sum((i,j)$ (ord(i)<>ord(j)),Eprim(i,j)*x(i,j));

co1(i2)                                                    ..  sum((j),x(i2,j))=e=1;
co2(j2)                                                    ..  sum((i),x(i,j2))=e=1;
co3                                                        ..  sum((j),x('0',j))=e=1;
co4                                                        ..  sum((i),x(i,'0'))=e=1;
co5(i,j)$((ord(i)> 2)and (ord(j)> 2)and (ord(i) <>ord(j))) ..  u(i)-u(j)=l= card(j) - card(j)*x(i,j) - 1;
co6(j2)                                                    ..  wprim(j2)=e= w("0")*x("0",j2)+sum(i2,wprim(i2)*x(i2,j2))+w(j2);
*co7(i2,j2)                                                ..  wprim(j2)=g=b(i2,j2)+w(j2);    w("0")*x("0",j2)+  sum(i2,b(i2,j2))
co8(i,j) $( ord(i)<>ord(j))                                ..  deltaE(i,j)=e= 0.01*d(i,j)*(10.08+(1.15*wprim(i)*a)+(0.505*s*wprim(i)*a));
co9(i,j) $( ord(i)<>ord(j))                                ..  Eprim(i,j)=e= (0.036*t(i,j))+deltaE(i,j);

*co99(i,j)$( ord(i)<>ord(j))                                ..  E(i,j)*t(i,j)=e= Eprim(i,j);

co10(i,j)$( ord(i)<>ord(j))                                ..  v(i,j) =e=vmax*(1-v1*wprim(i));
co11(i,j)$( ord(i)<>ord(j))                                ..  t(i,j)*v(i,j)=e=d(i,j);


*co19(i2,j)                                                  ..  b(i2,j)-M*x(i2,j)=l=0;
*co20(i2,j)                                                  ..  b(i2,j)-wprim(i2)=l=0;
*co21(i2,j)                                                  ..  b(i2,j)=g=wprim(i2)- M*(1-x(i2,j));

co22(i,j)                                                  .. t(i,j)=g=0;
co23(i,j)                                                  .. Eprim(i,j)=g=0;

co31(i)                                                    ..  x(i,i)=e=0;
co41(i,j)                                                  ..  x(i,j)+x(j,i)=l=1;
co61                                                       ..  wprim("0")=e=w("0");

Model example / all /;
option limcol=100;
option limrow=30 , reslim=2000;
option optca=0 , optcr=0;
Option MINLP=lindo;

$sTitle eps-constraint Method
Set
   k1(k)  'the first element of k'
   km1(k) 'all but the first elements of k'
   kk(k)  'active objective function in constraint allobj';

k1(k)$(ord(k) = 1) = yes;
km1(k)  = yes;
km1(k1) =  no;

Parameter
   rhs(k)    'right hand side of the constrained obj functions in eps-constraint'
   maxobj(k) 'maximum value from the payoff table'
   minobj(k) 'minimum value from the payoff table'
   numk(k)   'ordinal value of k starting with 1';

Scalar
   iter         'total number of iterations'
   infeas       'total number of infeasibilities'
   elapsed_time 'elapsed time for payoff and e-sonstraint'
   start        'start time'
   finish       'finish time';

Variable
   a_objval 'auxiliary variable for the objective function'
   obj      'auxiliary variable during the construction of the payoff table'
   sl(k)    'slack or surplus variables for the eps-constraints';

Positive Variable sl;

Equation
   con_obj(k) 'constrained objective functions'
   augm_obj   'augmented objective function to avoid weakly efficient solutions'
   allobj     'all the objective functions in one expression';

con_obj(km1).. z(km1) - dir(km1)*sl(km1) =e= rhs(km1);

* We optimize the first objective function and put the others as constraints
* the second term is for avoiding weakly efficient points

augm_obj..
   a_objval =e= sum(k1,dir(k1)*z(k1))
         + 1e-3*sum(km1,power(10,-(numk(km1) - 1))*sl(km1)/(maxobj(km1) - minobj(km1)));

allobj.. sum(kk, dir(kk)*z(kk)) =e= obj;

Model
   mod_payoff    / example, allobj            /
   mod_epsmethod / example, con_obj, augm_obj /;

Parameter payoff(k,k) 'payoff tables entries';

Alias (k,kp);

option optCr = 0, limRow = 0, limCol = 0, solPrint = off, solveLink = %solveLink.LoadLibrary%;

* Generate payoff table applying lexicographic optimization
loop(kp,
   kk(kp) = yes;
   repeat
      solve mod_payoff using minlp maximizing obj;
      payoff(kp,kk) = z.l(kk);
      z.fx(kk) = z.l(kk); // freeze the value of the last objective optimized
      kk(k++1) = kk(k);   // cycle through the objective functions
   until kk(kp);
   kk(kp) = no;
*  release the fixed values of the objective functions for the new iteration
   z.up(k) =  inf;
   z.lo(k) = -inf;
);
if(mod_payoff.modelStat <> %modelStat.optimal% and
   mod_payoff.modelStat <> %modelStat.integer Solution%,
   abort 'no optimal solution for mod_payoff';);

File fx / 2kp50_augmecon2_results.txt /;
put  fx ' PAYOFF TABLE'/;
loop(kp,
   loop(k, put payoff(kp,k):12:2;);
   put /;
);

minobj(k) = smin(kp,payoff(kp,k));
maxobj(k) = smax(kp,payoff(kp,k));

* gridpoints are calculated as the range (difference between max and min) of
* the 2nd objective function from the payoff table
$if not set gridpoints $set gridpoints 20
Set
   g         'grid points' / g0*g15 /
   grid(k,g) 'grid';

Parameter
   gridrhs(k,g) 'RHS of eps-constraint at grid point'
   maxg(k)      'maximum point in grid for objective'
   posg(k)      'grid position of objective'
   firstOffMax  'some counters'
   lastZero     'some counters'
*  numk(k) 'ordinal value of k starting with 1'
   numg(g)      'ordinal value of g starting with 0'
   step(k)      'step of grid points in objective functions'
   jump(k)      'jumps in the grid points traversing';

lastZero = 1;
loop(km1,
   numk(km1) = lastZero;
   lastZero  = lastZero + 1;
);
numg(g) = ord(g) - 1;

grid(km1,g) = yes; // Here we could define different grid intervals for different objectives
maxg(km1)   = smax(grid(km1,g), numg(g));
step(km1)   = (maxobj(km1) - minobj(km1))/maxg(km1);
gridrhs(grid(km1,g))$(dir(km1) = -1) = maxobj(km1) - numg(g)/maxg(km1)*(maxobj(km1) - minobj(km1));
gridrhs(grid(km1,g))$(dir(km1) =  1) = minobj(km1) + numg(g)/maxg(km1)*(maxobj(km1) - minobj(km1));

put / ' Grid points' /;
loop(g,
   loop(km1, put gridrhs(km1,g):12:2;);
   put /;
);
put / 'Efficient solutions' /;

* Walk the grid points and take shortcuts if the model becomes infeasible or
* if the calculated slack variables are greater than the step size
posg(km1) = 0;
iter   = 0;
infeas = 0;
start  = jnow;

repeat
   rhs(km1) = sum(grid(km1,g)$(numg(g) = posg(km1)), gridrhs(km1,g));
   solve mod_epsmethod maximizing a_objval using minlp;
   iter = iter + 1;
   if(mod_epsmethod.modelStat<>%modelStat.optimal% and
      mod_epsmethod.modelStat<>%modelStat.integer Solution%,
      infeas = infeas + 1; // not optimal is in this case infeasible
      put iter:5:0, '  infeasible' /;
      lastZero = 0;
      loop(km1$(posg(km1)  > 0 and lastZero = 0), lastZero = numk(km1));
      posg(km1)$(numk(km1) <= lastZero) = maxg(km1); // skip all solves for more demanding values of rhs(km1)
   else
      put iter:5:0;
      loop(k, put z.l(k):12:2;);
      jump(km1) = 1;
*     find the first off max (obj function that hasn't reach the final grid point).
*     If this obj.fun is k then assign jump for the 1..k-th objective functions
*     The jump is calculated for the innermost objective function (km=1)
      jump(km1)$(numk(km1) = 1) = 1 + floor(sl.L(km1)/step(km1));
      loop(km1$(jump(km1)  > 1), put '   jump';);
      put /;
   );
*  Proceed forward in the grid
   firstOffMax = 0;
   loop(km1$(posg(km1) < maxg(km1) and firstOffMax = 0),
      posg(km1)   = min((posg(km1) + jump(km1)),maxg(km1));
      firstOffMax = numk(km1);
   );
   posg(km1)$(numk(km1) < firstOffMax) = 0;
   abort$(iter > 1000) 'more than 1000 iterations, something seems to go wrong';
until sum(km1$(posg(km1) = maxg(km1)),1) = card(km1) and firstOffMax = 0;

finish = jnow;
elapsed_time = (finish - start)*60*60*24;

put /;
put 'Infeasibilities = ', infeas:5:0 /;
put 'Elapsed time: ',elapsed_time:10:2, ' seconds' /;
