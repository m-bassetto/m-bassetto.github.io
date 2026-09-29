function z=kstexo(kstar,bet,gam,R,cfpio,cfpik,vvec,mink,maxk);
% Copyright by Marco Bassetto, 1998-2006. This code can be freely
% distributed and modified for research purposes only, 
% provided this copyright notice is included in the modified code. 
% Proper credit should be given in all publications arising from
% modifications of this code; this should include a citation of 
% "Political Economy of Taxation in an Overlapping-Generations Economy,"
% by Marco Bassetto
%This function computes the residual for solving for
%k^*, i.e. the capital saved at the threat point

temp=cos(vvec*acos(2*(kstar-mink)/(maxk-mink)-1));
z=kstar-(bet/gam)+temp*cfpio/(R*(1-temp*cfpik));

