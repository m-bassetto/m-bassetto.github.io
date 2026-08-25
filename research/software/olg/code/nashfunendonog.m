function [prdiff,dprdiff]=nashfunendonog(temp,alph,bet,gam,n,cfnetrk,...
   cfpio,uyth,uoth,k,vvec,mink,maxk,cprod);
% Copyright by Marco Bassetto, 1998-2007. This code can be freely
% distributed and modified for research purposes only, 
% provided this copyright notice is included in the modified code. 
% Proper credit should be given in all publications arising from
% modifications of this code; this should include a citation of 
% "Political Economy of Taxation in an Overlapping-Generations Economy,"
% by Marco Bassetto
%This function computes the product of the difference
%in the utilities to find the Nash bargaining solution in the case of no
% government
% It also computes the relevant gradient
kp=temp(1);
ty=temp(2);

temp1=cos(vvec*acos(2*(kp-mink)/(maxk-mink)-1));
netrex=(temp1*cfnetrk)^(alph-1);
toex=temp1*cfpio;

coex=netrex*kp+toex;

netw=gam*coex/(bet*netrex);
cy=netw/gam;
l=(cy+kp-ty)/netw;

dtemp1=((maxk-mink)^2/4-(kp-(mink+maxk)/2)^2)^(-1/2)*vvec.*sin(...
   vvec*acos(2*(kp-mink)/(maxk-mink)-1));
dnetrex=(alph-1)*[dtemp1*cfnetrk; 0]*(temp1*cfnetrk)^(alph-2);
dtoex=[dtemp1*cfpio; 0];
dcoex=dnetrex*kp+netrex*[1; 0]+dtoex;

dnetw=gam*dcoex/(bet*netrex)-gam*coex*dnetrex/(bet*netrex^2);
dcy=dnetw/gam;
dl=(dcy+[1; 0]-[0; 1])/netw-(cy+kp-ty)*dnetw/netw^2;

ln=l;
dln=dl;
if l<0,
   l=0;
   dl=0;
   co=0;
   dco=zeros(2,1);
else,
    co=cprod*k^alph*(l*n)^(1-alph)-(cy+kp)*n;
    dco=(1-alph)*cprod*(k/l)^alph*n^(1-alph)*dl-(dcy+[1; 0])*n;
end;

con=co;
dcon=dco;
if co<.0001,
   co=0.0001;
   dco=zeros(2,1);
end;

uy=log(cy)+bet*log(coex)-gam*l;
duy=dcy/cy+bet*dcoex/coex-gam*dl;

diffuy=uy-uyth;

diffuo=log(co)-uoth;
duo=dco/co;

prdiff=0;
dprdiff=zeros(2,1);
if diffuy>=0 & diffuo>=0 & ln>=0 & co>=0,
   prdiff=-diffuy^n*diffuo;
   dprdiff=-(diffuy^n*duo+n*diffuy^(n-1)*duy*diffuo);
%   disp('ok');
else,
   if diffuy<0,
%       disp('Young');
      prdiff=prdiff+1e+7*diffuy^2;
      dprdiff=dprdiff+2e+7*diffuy*duy;
   end;
   if diffuo<0,
%      disp('Old');
      prdiff=prdiff+1e+7*diffuo^2;
      dprdiff=dprdiff+2e+7*diffuo*duo;
   end;
   if ln<0,
%      disp('l');
      prdiff=prdiff+1e+7*ln^2;
      dprdiff=dprdiff+2e+7*ln*dln;
   end;
   if con<.0001,
%      disp('co');
       prdiff=prdiff+1e+7*(con-.0001)^2;
       dprdiff=dprdiff+2e+7*(con-.0001)*dcon;
   end;
end;
