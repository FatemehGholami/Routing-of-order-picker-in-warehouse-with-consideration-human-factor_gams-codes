sets
i /0*6/
i2(i) index of products(vertex) /1*6/
alias(i,j);
alias(i2,j2);
sets
iter iteration /1*11/;

parameters
w(i) weighte of product i   /0 8,1 20,2 10,3 24,4 10,5 8,6 6/
ww /-0.1/
result(iter,*)

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

   0  0   9    21   14   28   23   33
   1  9   0    12   17   29   26   30
   2  21  12   0    29   17   28   18
   3  14  17   29   0    14   13   21
   4  28  29   17   14   0    17   7
   5  23  26   28   13   17   0    10
   6  33  30   18   21   7    10   0
   ;




variables
z1
z2
zf
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
ObjectiveFunction3
co1(i)
co2(j)
co3
co4
co5(i,j)
co6(j)
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

ObjectiveFunction1                                         ..  z1=e=sum((i,j)$ (ord(i)<>ord(j)),t(i,j)*x(i,j));
ObjectiveFunction2                                         ..  z2=e=sum((i,j)$ (ord(i)<>ord(j)),Eprim(i,j)*x(i,j));
ObjectiveFunction3                                         ..  zf=e=ww*((sum((i,j)$ (ord(i)<>ord(j)),t(i,j)*x(i,j))-z1.l)/z1.l)+(1-ww)*((sum((i,j)$ (ord(i)<>ord(j)), E(i,j)*x(i,j))-z2.l)/z2.l);

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

*model ssf /ObjectiveFunction3,co1,co2,co3,co4,co5,co8,co9,co91,co10,co11,co12,co19,co20,co21/;
*u.fx(j)$(ord(j) = 1) = 0;
*u.up(j) = card(j) - 1;


model ss1 /ObjectiveFunction1,co1,co2,co3,co4,co5,co6,co8,co9,co10,co11,CO22,co23,co31,co41,co61/
option optca=0,optcr=0,minlp=baron,reslim=2500;

model ss2 /ObjectiveFunction2,co1,co2,co3,co4,co5,co6,co8,co9,co10,co11,co22,CO23,co31,co41,co61/
option optca=0,optcr=0,minlp=baron,reslim=2500;

model ssf /ObjectiveFunction3,co1,co2,co3,co4,co5,co6,co8,co9,co10,co11,CO22,CO23,co31,co41,co61/
option optca=0,optcr=0,minlp=baron,reslim=2500;

*solve ss1 using minlp min z1;
*solve ss2 using minlp min z2;

loop(iter,
ww=ww+.1;

solve ss1 using minlp min z1;
solve ss2 using minlp min z2;
solve ssf using minlp min zf;


result(iter,'x1')=sum((i,j)$ (ord(i)<>ord(j)),t.l(i,j)*x.l(i,j));

result(iter,'x2')=sum((i,j)$ (ord(i)<>ord(j)),Eprim.l(i,j)*x.l(i,j));

result(iter,'w1')=ww;
result(iter,'w2')=1-ww;

);
display result;






*model ss2 /all/;

*model ss2 /all/;
*option limcol=100;
*option limrow=30 , reslim=2500;
*option optca=0 , optcr=0;
*Option MINLP=baron;


*solve ss2 using minlp min z1;
*solve ss2 using minlp min z1;
*display x.l,z1.l,t.l;
