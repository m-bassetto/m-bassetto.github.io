function [cts,ctseq,dcts,dctseq]=nashconendo(temp,alph,bet,hbet,gam,n,phi,sig,cfnetrk,...
   cfpio,cfpig,uyth,uoth,k,vvec,mink,maxk,cprod);
% Copyright by Marco Bassetto, 1998-2006. This code can be freely
% distributed and modified for research purposes only, 
% provided this copyright notice is included in the modified code. 
% Proper credit should be given in all publications arising from
% modifications of this code; this should include a citation of 
% "Political Economy of Taxation in an Overlapping-Generations Economy,"
% by Marco Bassetto
%This function computes the difference of the utilities
%for the two players, government spending, and the labor supply, 
%to be returned as constraints (all of these must be nonnegative).
%It also computes the appropriate gradient.

kp=temp(1);
co=temp(2);
ty=temp(3);

temp1=cos(vvec*acos(2*(kp-mink)/(maxk-mink)-1));
netrex=(temp1*cfnetrk)^(alph-1);
toex=temp1*cfpio;
pigp=temp1*cfpig;

coex=netrex*kp+toex;

netw=gam*coex/(bet*netrex);
cy=netw/gam;
l=(cy+kp-ty)/netw;

dtemp1=((maxk-mink)^2/4-(kp-(mink+maxk)/2)^2)^(-1/2)*vvec.*sin(...
   vvec*acos(2*(kp-mink)/(maxk-mink)-1));
dnetrex=(alph-1)*[dtemp1*cfnetrk; 0; 0]*(temp1*cfnetrk)^(alph-2);
dtoex=[dtemp1*cfpio; 0; 0];
dcoex=dnetrex*kp+netrex*[1; 0; 0]+dtoex;
dpigp=[dtemp1*cfpig; 0; 0];

dnetw=gam*dcoex/(bet*netrex)-gam*coex*dnetrex/(bet*netrex^2);
dcy=dnetw/gam;
dl=(dcy+[1; 0; 0]-[0; 0; 1])/netw-(cy+kp-ty)*dnetw/netw^2;

ln=l;
dln=dl;
if l<0,
   l=0;
   dl=0;
end;

g=(cprod*k^alph*(l*n)^(1-alph)-(cy+kp)*n-co)/(1+n);
dg=((1-alph)*cprod*(k/l)^alph*n^(1-alph)*dl-(dcy+[1; 0; 0])*n-[0; 1; 0])/(1+n);
gn=g;
dgn=dg;

if g<0,
   g=0;
   dg=0;
end;

if g==0,
   gs=0;
else;
   gs=g^(-sig);
end;

uy=log(cy)+bet*log(coex)-gam*l+phi*g^(1-sig)/(1-sig)+...
   hbet*phi*pigp^(1-sig)/(1-sig);
duy=dcy/cy+bet*dcoex/coex-gam*dl+phi*gs*dg+hbet*phi*...
   pigp^(-sig)*dpigp;

diffuy=uy-uyth;

diffuo=log(co)+(hbet/bet)*phi*g^(1-sig)/(1-sig)-uoth;
duo=[0; 1; 0]/co+(hbet/bet)*phi*gs*dg;

cts=[-diffuy; -diffuo; -gn; -ln];
dcts=-[duy duo dgn dln];

ctseq=[];
dctseq=[];
