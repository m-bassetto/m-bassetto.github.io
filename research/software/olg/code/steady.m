function y=steady(k,mink,maxk,vvec,cfkp);

% Copyright by Marco Bassetto, 1998-2006. This code can be freely
% distributed and modified for research purposes only, 
% provided this copyright notice is included in the modified code. 
% Proper credit should be given in all publications arising from
% modifications of this code; this should include a citation of 
% "Political Economy of Taxation in an Overlapping-Generations Economy,"
% by Marco Bassetto
% This function computes the difference between future and current capital
% to compute a steady state

temp=cos(acos(2*(k-mink)/(maxk-mink)-1)*vvec);
kp=temp*cfkp;
y=k-kp;