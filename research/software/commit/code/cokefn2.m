function y=cokefn2(x,alph1,alph2,alph3,gam,qdl,qdh,L,B);

%This function is the one that has to be 0 at lhat
%L is the labor supplied by everybody else

% Copyright by Marco Bassetto
% This code can be freely distributed and modified for research purposes only, 
% provided this copyright notice is included in the modified code. 
% Proper credit should be given in all publications arising from
% modifications of this code; this should include a citation of 
% "Equilibrium and Government Commitment" by Marco Bassetto

% We use a trapezoid approximation to the integral
% Notice that the integral is not smooth, so
% Gauss-Chebyshev might not be such a good idea.

%Number of points
n=200;
q=linspace(-qdl,qdh,n);
%[qdl x B]

%At each point, I need to evaluate whether the government
%defaults
if B<L-qdl,
    prodfactor=max(alph1,(alph1^(1/gam)+alph2^(1/gam))^gam*...
        (1-B./(L+q)).^(1-gam));
else,
    prodfactor=alph1*ones(1,n);
end;

y=alph3-sum(prodfactor(2:n-1).*(x+q(2:n-1)).^(-gam))/n-prodfactor(1)*(x-qdl)/(2*n)-...
   -prodfactor(n)*(x+qdh)/(2*n);