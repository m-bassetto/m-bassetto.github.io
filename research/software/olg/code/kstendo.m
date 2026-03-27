function [z,dz]=kstendo(kl,alph,bet,gam,n,cfpio,cfnetrk,vvec,mink,...
   maxk,k,cprod);
% Copyright by Marco Bassetto, 1998-2006. This code can be freely
% distributed and modified for research purposes only, 
% provided this copyright notice is included in the modified code. 
% Proper credit should be given in all publications arising from
% modifications of this code; this should include a citation of 
% "Political Economy of Taxation in an Overlapping-Generations Economy,"
% by Marco Bassetto
%This function computes the residual for solving for
%k^* and w^* i.e. the capital saved at the threat point 
%and the labor supply with endogenous prices

kstar=kl(1);
lstar=kl(2);

temp=cos(vvec*acos(2*(kstar-mink)/(maxk-mink)-1));
dtemp=sin(vvec*acos(2*(kstar-mink)/(maxk-mink)-1)).*vvec/...
   sqrt((kstar-mink)*(maxk-kstar));

wstar=(1-alph)*cprod*(k/(n*lstar))^alph;
dwstar=[0; -alph*(1-alph)*cprod*(k/n)^alph*lstar^(-alph-1)];

z1=kstar-bet*wstar/gam+temp*cfpio/(temp*cfnetrk)^(alph-1);
dz1=[1; 0]-bet*dwstar/gam+[dtemp*cfpio/(temp*cfnetrk); 0]+(1-alph)*...
   [dtemp*cfnetrk*(temp*cfpio); 0]*(temp*cfnetrk)^(-alph);

z2=lstar-1/gam-kstar/wstar;
dz2=[0; 1]-[1/wstar; 0]+kstar*dwstar/wstar^2;  

z=z1^2+z2^2;
dz=2*z1*dz1+2*z2*dz2;
