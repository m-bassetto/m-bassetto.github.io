function [cts,ctseq,dcts,dctseq]=nashconexo(temp,bet,hbet,gam,R,n,phi,sig,cfpik,...
   cfpio,cfpig,uyth,k,vvec,mink,maxk);
% Copyright by Marco Bassetto, 1998-2006. This code can be freely
% distributed and modified for research purposes only, 
% provided this copyright notice is included in the modified code. 
% Proper credit should be given in all publications arising from
% modifications of this code; this should include a citation of 
% "Political Economy of Taxation in an Overlapping-Generations Economy,"
% by Marco Bassetto
% This function computes the difference of the utilities
%for the two players, to be returned as constraints
%(both must be nonnegative), and government spending
%(for the same reason). It also computes the appropriate gradient.

kp=temp(1);
co=temp(2);
ty=temp(3);

temp1=cos(vvec*acos(2*(kp-mink)/(maxk-mink)-1));
netrex=R*(1-temp1*cfpik);
toex=temp1*cfpio;
pigp=temp1*cfpig;

coex=netrex*kp+toex;

g=n/(gam*(1+n))+bet*n*netrex*(kp-ty)/(gam*coex*(1+n))-...
   n*coex/(bet*(1+n)*netrex)-kp*n/(1+n)+(R*k-co)/(1+n);

dtemp1=((maxk-mink)^2/4-(kp-(mink+maxk)/2)^2)^(-1/2)*vvec.*sin(...
   vvec*acos(2*(kp-mink)/(maxk-mink)-1));
dnetrex=[-R*dtemp1*cfpik; 0; 0];
dtoex=[dtemp1*cfpio; 0; 0];
dcoex=dnetrex*kp+netrex*[1; 0; 0]+dtoex;
dpigp=[dtemp1*cfpig; 0; 0];

temp2=dnetrex*(kp-ty)+netrex*[1; 0; -1];

gn=g;

dg=(bet*n/(1+n))*(temp2/(gam*coex)-netrex*(kp-ty)*dcoex/(gam*coex^2))-...
   (n/(bet*(1+n)))*(dcoex/netrex-coex*dnetrex/netrex^2)-...
   (n/(1+n))*[1; 0; 0]+[0; -1; 0]/(1+n);
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

uy=(1+bet)*log(coex)-log(bet)-log(netrex)-1-...
   bet*netrex*(kp-ty)/coex+phi*g^(1-sig)/(1-sig)+...
   hbet*phi*pigp^(1-sig)/(1-sig);
duy=(1+bet)*dcoex/coex-dnetrex/netrex-bet*(temp2/coex-...
   netrex*(kp-ty)*dcoex/coex^2)+...
   phi*gs*dg+hbet*phi*pigp^(-sig)*dpigp;

diffuy=uy-uyth;
diffuo=log(co)-log(R*k)+(hbet/bet)*phi*g^(1-sig)/(1-sig);
duo=[0; 1; 0]/co+(hbet/bet)*phi*gs*dg;

dprdiff=zeros(3,1);
if diffuy>0 & diffuo>0 & gn>=0,
   dprdiff=-(diffuy^n*duo+n*diffuy^(n-1)*duy*diffuo);
else,
   if diffuy<0,
      dprdiff=dprdiff+2e+7*diffuy*duy;
   end;
   if diffuo<0,
      dprdiff=dprdiff+2e+7*diffuo*duo;
   end;
   if gn<0,
      dprdiff=dprdiff+2e+7*gn*dgn;
   end;
end;

cts=[-diffuy; -diffuo; -gn];
dcts=-[duy duo dgn];
ctseq=[];
dctseq=[];
