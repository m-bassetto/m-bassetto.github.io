% Main code
% Copyright by Marco Bassetto, 1998-2006. This code can be freely
% distributed and modified for research purposes only, 
% provided this copyright notice is included in the modified code. 
% Proper credit should be given in all publications arising from
% modifications of this code; this should include a citation of 
% "Political Economy of Taxation in an Overlapping-Generations Economy,"
% by Marco Bassetto

warning off optim:fmincon:SwitchingToMediumScale
% Part I: Solution of the stationary policy with exogenous prices
%Initialise the parameters
gam=1;
bet=.2181;
hbet=.1119;
phi=.497;
%load phi;
sig=.6977;
R=(1.06^30);
n=2;
gridk=15;
grph=0;
%gridpi must be equal to gridk in this version.
%We can put gridpi less than gridk using regression
%instead of interpolation, but the code needs to be changed
%because val is no longer square and cannot be inverted
gridpi=gridk;
mink=bet/(3*gam);
maxk=bet/gam+.1;
relax=.3;

% We use Chebyshev interpolation on [mink,maxk]
rts=cos((2*(1:gridk)-1)*pi/(2*gridk))';
rts=sort(rts);
% the value of each polynomial is on a column
vvec=0:gridpi-1;
val=cos(acos(rts)*vvec);
k=((maxk-mink)*rts+maxk+mink)/2;
ival=inv(val);

%Initial guess for the policy
cfpik=zeros(gridpi,1);
cfpio=zeros(gridpi,1);
cfpig=[.001; zeros(gridpi-1,1)];

%Lower and upper bounds (should not play any role,except
%for Ty);
vlb=[mink+.0001; 1e-15; 0];
vub=[maxk-.0001; 200; 200];

%Initialisation
%init=[bet*ones(1,gridk)/gam+.00001; R*k'-.002; .001*ones(1,gridk)];

kpopt=zeros(gridk,1);
coopt=zeros(gridk,1);
tyopt=zeros(gridk,1);
taulopt=zeros(gridk,1);
taukopt=zeros(gridk,1);
taukoptold=zeros(gridk,1);
toopt=zeros(gridk,1);
tooptold=zeros(gridk,1);
gopt=zeros(gridk,1);
goptold=zeros(gridk,1);
zx=linspace(mink,maxk);
%zx=k';

epsi=100;

tic;
ite=1;
while epsi>1e-08, 
  
   %the relevant threat point depends very little on k_t
   %compute most quantities outside the k_t loop
   %first, solve for k^*
   options=optimset('TolX',1e-16,'Display','off');
   [kstar,fval,exitflag]=fzero(@(kstarsolve) kstexo(kstarsolve,bet,gam,R,cfpio,cfpik,vvec,...
      mink,maxk), [mink; maxk],options);
  if exitflag<0,
      exitflag
      error('Failure to converge');
  end;
   temp=cos(vvec*acos(2*(kstar-mink)/(maxk-mink)-1));
   pikstar=temp*cfpik;
   piostar=temp*cfpio;
   pigstar=temp*cfpig;
   init=[kstar*ones(1,gridk); R*k'-.02; .001*ones(1,gridk)];
   uyth=-(1+bet)*log(gam)+bet*log(bet)+bet*log(R)+bet*...
      log(1-pikstar)-(1+bet)+gam*piostar/(R*(1-pikstar))+hbet*phi*...
      pigstar^(1-sig)/(1-sig);

   %We now recompute the optimal policy at the grid points
   options=optimset('GradConstr','on','Gradobj','on','TolX',1e-16,'Display','off');
   for i=1:gridk,
      [temp,fval,exitflag]=fmincon(@(tempsol) nashfunexo(tempsol,bet,hbet,gam,R,n,phi,sig,cfpik,...
          cfpio,cfpig,uyth,k(i),vvec,mink,maxk),init(:,i),[],[],[],[],vlb,vub,...
          @(tempsol) nashconexo(tempsol,bet,hbet,gam,R,n,phi,sig,cfpik,...
          cfpio,cfpig,uyth,k(i),vvec,mink,maxk),options);
      if exitflag<1,
          exitflag
          error('Failure to converge');
      end;
      %temp
      kpopt(i)=temp(1);
      vf(i)=nashfunexo(temp,bet,hbet,gam,R,n,phi,sig,cfpik,cfpio,cfpig,uyth,...
         k(i),vvec,mink,maxk);
      temp2=cos(vvec*acos(2*(kpopt(i)-mink)/(maxk-mink)-1));
      pikopt=temp2*cfpik;
      pioopt=temp2*cfpio;
      pigopt=temp2*cfpig;

      coopt(i)=temp(2);
      tyopt(i)=temp(3);
      netrex=R*(1-pikopt);
      copopt=netrex*kpopt(i)+pioopt;
      gopt(i)=n/(gam*(1+n))+bet*n*netrex*(kpopt(i)-tyopt(i))/...
         (gam*copopt*(1+n))-n*copopt/(bet*(1+n)*netrex)-n*kpopt(i)/...
         (1+n)+(R*k(i)-coopt(i))/(1+n);
      diffuo=log(coopt(i))+phi*(hbet/bet)*gopt(i)^(1-sig)/(1-sig)...
         -log(R*k(i));
      taukopt(i)=1-coopt(i)/(R*k(i)*(1+diffuo));
      toopt(i)=coopt(i)*diffuo/(1+diffuo);
      %init(:,i)=temp;
   end;
   cfpik2=ival*taukopt;
   cfpig2=ival*gopt;
   cfpio2=ival*toopt;
   epsi2=max(max(abs([cfpik cfpig cfpio]-[cfpik2 cfpig2 cfpio2])))
%   epsi2=max(max(abs(([taukopt gopt toopt]-[taukoptold goptold tooptold])./...
%       [taukopt gopt toopt])))
   ite
   if epsi2<epsi,
       taukoptold=taukopt;
       goptold=gopt;
       tooptold=toopt;
       cfpikold=cfpik;
       cfpigold=cfpig;
       cfpioold=cfpio;
       cfpik=relax*cfpik2+(1-relax)*cfpik;
       cfpig=relax*cfpig2+(1-relax)*cfpig;
       cfpio=relax*cfpio2+(1-relax)*cfpio;
       epsi=epsi2;
       ite=ite+1;
       relax=.2;
   else,
       relax=relax/2;
       cfpik=relax*cfpik2+(1-relax)*cfpikold;
       cfpig=relax*cfpig2+(1-relax)*cfpigold;
       cfpio=relax*cfpio2+(1-relax)*cfpioold;
       if relax<.0001,
           error('Failure to converge');
       end;
   end;
   if grph==1,
      subplot(2,2,1);
      temp=cos(acos(2*(zx-mink)/(maxk-mink)-1)'*vvec);
      plot(zx,temp*cfpik,'-',zx,temp*cfpik2,':');
      title('taukopt');
      subplot(2,2,2);
      plot(zx,temp*cfpig,zx,temp*cfpig2,':');
      title('gopt');
      subplot(2,2,3);
      plot(zx,temp*cfpio,zx,temp*cfpio2,':');
      title('toopt');
      subplot(2,2,4);
      plot(0,10000*epsi,'+');
      pause;
   end;

   toc;
end;

%Compute variables not required in the iterations
temp=cos(acos(2*(kpopt-mink)/(maxk-mink)-1)*vvec);
pikopt=temp*cfpik;
pioopt=temp*cfpio;
pigopt=temp*cfpig;
taulopt=1-gam*kpopt/bet-gam*pioopt./(bet*R*(1-pikopt));
lopt=(1+bet)/gam-tyopt./(1-taulopt)-pioopt./(R*(1-pikopt).*(1-taulopt));
cyopt=(1-taulopt)/gam;
copopt=R*(1-pikopt).*kpopt+pioopt;
uyopt=log(cyopt)+bet*log(copopt)-gam*lopt+phi*gopt.^(1-sig)/(1-sig)+...
   hbet*phi*pigopt.^(1-sig)/(1-sig);
uoopt=log(coopt)+phi*(hbet/bet)*gopt.^(1-sig)/(1-sig);
uoth=log(R*k);
gnp=cyopt*n/(1+n)+coopt/(1+n)+gopt+kpopt*n/(1+n);
gnp2=R*k/(1+n)+lopt*n/(1+n);
g2gnp=gopt./gnp;
cy2gnp=cyopt*n./(gnp*(1+n));
co2gnp=coopt./(gnp*(1+n));
to2gnp=toopt./(gnp*(1+n));

%Find the steady state
options=optimset('FunValCheck','on','TolX',1e-16);
cfkp=ival*kpopt;
kss=fzero(@(ksssolve) steady(ksssolve,mink,maxk,vvec,cfkp),...
    [mink+.0001; maxk-.0001],options);
if exitflag<0,
    exitflag
    error('Failure to converge');
end;
temp=cos(acos(2*(kss-mink)/(maxk-mink)-1)*vvec);
taukss=temp*cfpik;
toss=temp*cfpio;
gss=temp*cfpig;
temp2=temp*ival;
taulss=temp2*taulopt;
lss=temp2*lopt;
cyss=temp2*cyopt;
coss=temp2*coopt;
coss=temp2*coopt;
gnpss=temp2*gnp;
g2gnpss=gss/gnpss;
cy2gnpss=cyss*n/(gnpss*(1+n));
co2gnpss=coss/(gnpss*(1+n));
to2gnpss=toss/(gnpss*(1+n));
taukhat=1-((R*(1-taukss))^(1/30)-1)/(R^(1/30)-1);
oldss=R*kss-coss;
youngss=lss-cyss-kss;
ratss=oldss/(gss*(1+n));

subplot(2,2,1);
plot(k,g2gnp);
set(gca,'FontSize',10);
title(['Provision of the public good as a fraction of GNP']);
   xlabel('Capital per old person');
   axis([0.05 0.35 0.14 0.24]);

subplot(2,2,2);
plot(k,taulopt);
set(gca,'FontSize',10);
title(['Tax rate on labor income']);
   xlabel('Capital per old person');
   axis([0.05 0.35 0.2 0.3]);

subplot(2,2,3);
plot(k,1-((R*(1-taukopt)).^(1/30)-1)./(R^(1/30)-1));
set(gca,'FontSize',10);
   title(['Tax rate on capital income']);
   xlabel('Capital per old person');
   axis([0.05 0.35 0.15 0.3]);
   
subplot(2,2,4);
plot(k,to2gnp);
set(gca,'FontSize',10);
   title(['Transfers to the old as a fraction of GNP']);
   xlabel('Capital per old person');
   axis([0.05 0.35 0.04 0.14]);
   
print -deps fig1.eps
print -dps plots.ps

figure;
plot(k,taulopt.*lopt,'-',k,taukopt*R.*k-toopt,':');
xlabel('Capital per old person');
text(k(10),taulopt(10)*lopt(10),'Young',...
   'HorizontalAlignment','right','VerticalAlignment',...
   'bottom');
text(k(2),taukopt(2)*R*k(2)-toopt(2),'Old',...
   'HorizontalAlignment','left','VerticalAlignment',...
   'top');

print -deps fig2.eps
print -dps -append plots.ps
save olgpexo;

clear;


% Part II computation of the transition
% (This duplicates most of part I)

%This code computes the effects of a
%demographic transition
%Initialise the parameters
gam=1;
bet=.2181;
hbet=.1119;
phi=.497;
sig=.6977;
R=1.06^30;
% nn=[2; 11/9]; % Unanticipated experiment
nn=[2; 2; 11/9]; % Anticipated experiment
%N is the number of policy rules we have to keep track of
N=size(nn,1)-1;
%T is the number of transition periods I consider
T=8; 
gridk=15;
%gridpi must be equal to gridk in this version.
%We can put gridpi less than gridk using regression
%instead of interpolation, but the code needs to be changed
%because val is no longer square and cannot be inverted
gridpi=gridk;
mink=0.1;
maxk=bet/gam+.1;
relax=.3;

% We use Chebyshev interpolation on [mink,maxk]
rts=cos((2*(1:gridk)-1)*pi/(2*gridk))';
rts=sort(rts);
% the value of each polynomial is on a column
vvec=0:gridpi-1;
val=cos(acos(rts)*vvec);
k=((maxk-mink)*rts+maxk+mink)/2;
ival=inv(val);

%Initial guess for the policy
cfpik=zeros(gridpi,1);
cfpio=zeros(gridpi,1);
cfpig=[.001; zeros(gridpi-1,1)];

%Lower and upper bounds (should not play any role,except
%for Ty);
vlb=[mink+.0001; 1e-15; 0];
vub=[maxk-.0001; 200; 200];

%Initialisation
%init=[bet*ones(1,gridk)/gam+.00001; R*k'-.002; .001*ones(1,gridk)];

kpopt=zeros(gridk,1);
coopt=zeros(gridk,1);
tyopt=zeros(gridk,1);
taulopt=zeros(gridk,1);
taukopt=zeros(gridk,1);
taukoptold=zeros(gridk,1);
toopt=zeros(gridk,1);
tooptold=zeros(gridk,1);
gopt=zeros(gridk,1);
goptold=zeros(gridk,1);
zx=linspace(mink,maxk);
%zx=k';

ktr=zeros(T,1);
tauktr=zeros(T,1);
totr=zeros(T,1);
gtr=zeros(T,1);
taultr=zeros(T,1);
ltr=zeros(T,1);
cytr=zeros(T,1);
cotr=zeros(T,1);
gnptr=zeros(T,1);
uytr=zeros(T,1);
uotr=zeros(T,1);
uythtr=zeros(T,1);
uothtr=zeros(T,1);

kpl=zeros(gridpi,N);
taukpl=zeros(gridpi,N);
topl=zeros(gridpi,N);
gpl=zeros(gridpi,N);
taulpl=zeros(gridpi,N);
lpl=zeros(gridpi,N);
cypl=zeros(gridpi,N);
copl=zeros(gridpi,N);
gnppl=zeros(gridpi,N);
uypl=zeros(gridpi,N);
uopl=zeros(gridpi,N);

%First, we compute the initial steady-state
n=nn(1);
epsi=1;
ite=1;
tic;
while epsi>1e-08,
   %the relevant threat point depends very little on k_t
   %compute most quantities outside the k_t loop
   %first, solve for k^*
   options=optimset('TolX',1e-16,'Display','off');
   [kstar,fval,exitflag]=fzero(@(kstarsolve) kstexo(kstarsolve,bet,gam,R,cfpio,cfpik,vvec,...
      mink,maxk), [mink; maxk],options);
  if exitflag<0,
      exitflag
      error('Failure to converge');
  end;
   temp=cos(vvec*acos(2*(kstar-mink)/(maxk-mink)-1));
   pikstar=temp*cfpik;
   piostar=temp*cfpio;
   pigstar=temp*cfpig;
   init=[kstar*ones(1,gridk)+.00001; R*k'-.02; .001*ones(1,gridk)];
   uyth=-(1+bet)*log(gam)+bet*log(bet)+bet*log(R)+bet*...
      log(1-pikstar)-(1+bet)+gam*piostar/(R*(1-pikstar))+hbet*phi*...
      pigstar^(1-sig)/(1-sig);

   %We now recompute the optimal policy at the grid points
   options=optimset('GradConstr','on','Gradobj','on','TolX',1e-16,'Display','off');
   for i=1:gridk,
      [temp,fval,exitflag]=fmincon(@(tempsol) nashfunexo(tempsol,bet,hbet,gam,R,n,phi,sig,cfpik,...
          cfpio,cfpig,uyth,k(i),vvec,mink,maxk),init(:,i),[],[],[],[],vlb,vub,...
          @(tempsol) nashconexo(tempsol,bet,hbet,gam,R,n,phi,sig,cfpik,...
          cfpio,cfpig,uyth,k(i),vvec,mink,maxk),options);
      if exitflag<1,
          exitflag
          error('Failure to converge');
      end;
      kpopt(i)=temp(1);
      temp2=cos(vvec*acos(2*(kpopt(i)-mink)/(maxk-mink)-1));
      pikopt=temp2*cfpik;
      pioopt=temp2*cfpio;
      pigopt=temp2*cfpig;

      coopt(i)=temp(2);
      tyopt(i)=temp(3);
      netrex=R*(1-pikopt);
      copopt=netrex*kpopt(i)+pioopt;
      gopt(i)=n/(gam*(1+n))+bet*n*netrex*(kpopt(i)-tyopt(i))/...
         (gam*copopt*(1+n))-n*copopt/(bet*(1+n)*netrex)-n*kpopt(i)/...
         (1+n)+(R*k(i)-coopt(i))/(1+n);
      diffuo=log(coopt(i))+phi*(hbet/bet)*gopt(i)^(1-sig)/(1-sig)...
         -log(R*k(i));
      taukopt(i)=1-coopt(i)/(R*k(i)*(1+diffuo));
      toopt(i)=coopt(i)*diffuo/(1+diffuo);
   end;
   cfpik2=ival*taukopt;
   cfpig2=ival*gopt;
   cfpio2=ival*toopt;
      epsi2=max(max(abs([cfpik cfpig cfpio]-[cfpik2 cfpig2 cfpio2])))
%   epsi2=max(max(abs(([taukopt gopt toopt]-[taukoptold goptold tooptold])./...
%       [taukopt gopt toopt])))
   ite
   if epsi2<epsi,
       taukoptold=taukopt;
       goptold=gopt;
       tooptold=toopt;
       cfpikold=cfpik;
       cfpigold=cfpig;
       cfpioold=cfpio;
       cfpik=relax*cfpik2+(1-relax)*cfpik;
       cfpig=relax*cfpig2+(1-relax)*cfpig;
       cfpio=relax*cfpio2+(1-relax)*cfpio;
       epsi=epsi2;
       ite=ite+1;
       relax=.2;
   else,
       relax=relax/2;
       cfpik=relax*cfpik2+(1-relax)*cfpikold;
       cfpig=relax*cfpig2+(1-relax)*cfpigold;
       cfpio=relax*cfpio2+(1-relax)*cfpioold;
       if relax<.0001,
           error('Failure to converge');
       end;
   end;
   toc;
end;

%Compute variables not required in the iterations
temp=cos(acos(2*(kpopt-mink)/(maxk-mink)-1)*vvec);
pikopt=temp*cfpik;
pioopt=temp*cfpio;
pigopt=temp*cfpig;
taulopt=1-gam*kpopt/bet-gam*pioopt./(bet*R*(1-pikopt));
lopt=(1+bet)/gam-tyopt./(1-taulopt)-pioopt./(R*(1-pikopt).*(1-taulopt));
cyopt=(1-taulopt)/gam;
uyopt=log(cyopt)+bet*log(copopt)-gam*lopt+phi*gopt.^(1-sig)/(1-sig)+...
   hbet*phi*pigopt.^(1-sig)/(1-sig);
uoopt=log(coopt)+phi*(hbet/bet)*gopt.^(1-sig)/(1-sig);
uoth=log(R*k);
gnp=cyopt*n/(1+n)+coopt/(1+n)+gopt+kpopt*n/(1+n);

%Find the steady state
options=optimset('FunValCheck','on','TolX',1e-16);
cfkp=ival*kpopt;
ktr(1)=fzero(@(ksssolve) steady(ksssolve,mink,maxk,vvec,cfkp),...
    [mink+.0001; maxk-.0001],options);
if exitflag<0,
    exitflag
    error('Failure to converge');
end;
ktr(2)=ktr(1); 
temp=cos(acos(2*(ktr(1)-mink)/(maxk-mink)-1)*vvec);
tauktr(1)=temp*cfpik;
totr(1)=temp*cfpio;
gtr(1)=temp*cfpig;
temp2=temp*ival;
taultr(1)=temp2*taulopt;
ltr(1)=temp2*lopt;
cytr(1)=temp2*cyopt;
cotr(1)=temp2*coopt;
gnptr(1)=temp2*gnp;
uytr(1)=temp2*uyopt;
uotr(1)=temp2*uoopt;
uythtr(1)=uyth;
uothtr(1)=log(R*ktr(1));

%We now compute the policy rules in the final steady state

n=nn(N+1);
epsi=1;
ite=1;
while epsi>1e-08,
   %the relevant threat point depends very little on k_t
   %compute most quantities outside the k_t loop
   %first, solve for k^*
   options=optimset('TolX',1e-16,'Display','off');
   [kstar,fval,exitflag]=fzero(@(kstarsolve) kstexo(kstarsolve,bet,gam,R,cfpio,cfpik,vvec,...
      mink,maxk), [mink; maxk],options);
  if exitflag<0,
      exitflag
      error('Failure to converge');
  end;
   temp=cos(vvec*acos(2*(kstar-mink)/(maxk-mink)-1));
   pikstar=temp*cfpik;
   piostar=temp*cfpio;
   pigstar=temp*cfpig;
   init=[kstar*ones(1,gridk)+.00001; R*k'-.02; .001*ones(1,gridk)];
   uyth=-(1+bet)*log(gam)+bet*log(bet)+bet*log(R)+bet*...
      log(1-pikstar)-(1+bet)+gam*piostar/(R*(1-pikstar))+hbet*phi*...
      pigstar^(1-sig)/(1-sig);

   %We now recompute the optimal policy at the grid points
   options=optimset('GradConstr','on','Gradobj','on','TolX',1e-16,'Display','off');
   for i=1:gridk,
       [temp,fval,exitflag]=fmincon(@(tempsol) nashfunexo(tempsol,bet,hbet,gam,R,n,phi,sig,cfpik,...
          cfpio,cfpig,uyth,k(i),vvec,mink,maxk),init(:,i),[],[],[],[],vlb,vub,...
          @(tempsol) nashconexo(tempsol,bet,hbet,gam,R,n,phi,sig,cfpik,...
          cfpio,cfpig,uyth,k(i),vvec,mink,maxk),options);
      if exitflag<1,
          exitflag
          error('Failure to converge');
      end;
      kpopt(i)=temp(1);
      temp2=cos(vvec*acos(2*(kpopt(i)-mink)/(maxk-mink)-1));
      pikopt=temp2*cfpik;
      pioopt=temp2*cfpio;
      pigopt=temp2*cfpig;

      coopt(i)=temp(2);
      tyopt(i)=temp(3);
      netrex=R*(1-pikopt);
      copopt=netrex*kpopt(i)+pioopt;
      gopt(i)=n/(gam*(1+n))+bet*n*netrex*(kpopt(i)-tyopt(i))/...
         (gam*copopt*(1+n))-n*copopt/(bet*(1+n)*netrex)-n*kpopt(i)/...
         (1+n)+(R*k(i)-coopt(i))/(1+n);
      diffuo=log(coopt(i))+phi*(hbet/bet)*gopt(i)^(1-sig)/(1-sig)...
         -log(R*k(i));
      taukopt(i)=1-coopt(i)/(R*k(i)*(1+diffuo));
      toopt(i)=coopt(i)*diffuo/(1+diffuo);
   end;
   cfpik2=ival*taukopt;
   cfpig2=ival*gopt;
   cfpio2=ival*toopt;
      epsi2=max(max(abs([cfpik cfpig cfpio]-[cfpik2 cfpig2 cfpio2])))
%   epsi2=max(max(abs(([taukopt gopt toopt]-[taukoptold goptold tooptold])./...
%       [taukopt gopt toopt])))
   ite
   if epsi2<epsi,
       taukoptold=taukopt;
       goptold=gopt;
       tooptold=toopt;
       cfpikold=cfpik;
       cfpigold=cfpig;
       cfpioold=cfpio;
       cfpik=relax*cfpik2+(1-relax)*cfpik;
       cfpig=relax*cfpig2+(1-relax)*cfpig;
       cfpio=relax*cfpio2+(1-relax)*cfpio;
       epsi=epsi2;
       ite=ite+1;
       relax=.2;
   else,
       relax=relax/2;
       cfpik=relax*cfpik2+(1-relax)*cfpikold;
       cfpig=relax*cfpig2+(1-relax)*cfpigold;
       cfpio=relax*cfpio2+(1-relax)*cfpioold;
       if relax<.0001,
           error('Failure to converge');
       end;
   end;
   toc;
end;

%Compute variables not required in the iterations
temp=cos(acos(2*(kpopt-mink)/(maxk-mink)-1)*vvec);
pikopt=temp*cfpik;
pioopt=temp*cfpio;
pigopt=temp*cfpig;
taulopt=1-gam*kpopt/bet-gam*pioopt./(bet*R*(1-pikopt));
lopt=(1+bet)/gam-tyopt./(1-taulopt)-pioopt./(R*(1-pikopt).*(1-taulopt));
cyopt=(1-taulopt)/gam;
uyopt=log(cyopt)+bet*log(copopt)-gam*lopt+phi*gopt.^(1-sig)/(1-sig)+...
   hbet*phi*pigopt.^(1-sig)/(1-sig);
uoopt=log(coopt)+phi*(hbet/bet)*gopt.^(1-sig)/(1-sig);
uoth=log(R*k);
gnp=cyopt*n/(1+n)+coopt/(1+n)+gopt+kpopt*n/(1+n);

kpl(:,N)=ival*kpopt;
taukpl(:,N)=cfpik;
topl(:,N)=cfpio;
gpl(:,N)=cfpig;
taulpl(:,N)=ival*taulopt;
lpl(:,N)=ival*lopt;
cypl(:,N)=ival*cyopt;
copl(:,N)=ival*coopt;
gnppl(:,N)=ival*gnp;
uypl(:,N)=ival*uyopt;
uopl(:,N)=ival*uoopt;
uythtr(N+1:T)=uyth;

%Now we find the transition policy rules by backward induction

for inpl=N-1:-1:1,
   n=nn(inpl+1);
   %the relevant threat point depends very little on k_t
   %compute most quantities outside the k_t loop
   %first, solve for k^*
   options=optimset('TolX',1e-16,'Display','off');
   [kstar,fval,exitflag]=fzero(@(kstarsolve) kstexo(kstarsolve,bet,gam,R,cfpio,cfpik,vvec,...
      mink,maxk), [mink; maxk],options);
  if exitflag<0,
      exitflag
      error('Failure to converge');
  end;
   temp=cos(vvec*acos(2*(kstar-mink)/(maxk-mink)-1));
   pikstar=temp*cfpik;
   piostar=temp*cfpio;
   pigstar=temp*cfpig;
   init=[kstar*ones(1,gridk)+.00001; R*k'-.02; .001*ones(1,gridk)];
   uyth=-(1+bet)*log(gam)+bet*log(bet)+bet*log(R)+bet*...
      log(1-pikstar)-(1+bet)+gam*piostar/(R*(1-pikstar))+hbet*phi*...
      pigstar^(1-sig)/(1-sig);

   %We now recompute the optimal policy at the grid points
   options=optimset('GradConstr','on','Gradobj','on','TolX',1e-16,'Display','off');
   for i=1:gridk,
       [temp,fval,exitflag]=fmincon(@(tempsol) nashfunexo(tempsol,bet,hbet,gam,R,n,phi,sig,cfpik,...
          cfpio,cfpig,uyth,k(i),vvec,mink,maxk),init(:,i),[],[],[],[],vlb,vub,...
          @(tempsol) nashconexo(tempsol,bet,hbet,gam,R,n,phi,sig,cfpik,...
          cfpio,cfpig,uyth,k(i),vvec,mink,maxk),options);
      if exitflag<1,
          exitflag
          error('Failure to converge');
      end;
      kpopt(i)=temp(1);
      temp2=cos(vvec*acos(2*(kpopt(i)-mink)/(maxk-mink)-1));
      pikopt=temp2*cfpik;
      pioopt=temp2*cfpio;
      pigopt=temp2*cfpig;
      coopt(i)=temp(2);
      tyopt(i)=temp(3);
      netrex=R*(1-pikopt);
      copopt=netrex*kpopt(i)+pioopt;
      gopt(i)=n/(gam*(1+n))+bet*n*netrex*(kpopt(i)-tyopt(i))/...
         (gam*copopt*(1+n))-n*copopt/(bet*(1+n)*netrex)-n*kpopt(i)/...
         (1+n)+(R*k(i)-coopt(i))/(1+n);
      diffuo=log(coopt(i))+phi*(hbet/bet)*gopt(i)^(1-sig)/(1-sig)...
         -log(R*k(i));
      taukopt(i)=1-coopt(i)/(R*k(i)*(1+diffuo));
      toopt(i)=coopt(i)*diffuo/(1+diffuo);
   end;
   cfpik=ival*taukopt;
   cfpig=ival*gopt;
   cfpio=ival*toopt;
   inpl
   toc;

   %Compute variables not required in the iterations
   temp=cos(acos(2*(kpopt-mink)/(maxk-mink)-1)*vvec);
   pikopt=temp*taukpl(:,inpl+1);
   pioopt=temp*topl(:,inpl+1);
   pigopt=temp*gpl(:,inpl+1);

   taulopt=1-gam*kpopt/bet-gam*pioopt./(bet*R*(1-pikopt));
   lopt=(1+bet)/gam-tyopt./(1-taulopt)-pioopt./(R*(1-pikopt).*...
      (1-taulopt));
   cyopt=(1-taulopt)/gam;
   uyopt=log(cyopt)+bet*log(copopt)-gam*lopt+phi*gopt.^(1-sig)/(1-sig)+...
      hbet*phi*pigopt.^(1-sig)/(1-sig);
   uoopt=log(coopt)+phi*(hbet/bet)*gopt.^(1-sig)/(1-sig);
   uoth=log(R*k);
   gnp=cyopt*n/(1+n)+coopt/(1+n)+gopt+kpopt*n/(1+n);

   kpl(:,inpl)=ival*kpopt;
   taukpl(:,inpl)=cfpik;
   topl(:,inpl)=cfpio;
   gpl(:,inpl)=cfpig;
   taulpl(:,inpl)=ival*taulopt;
   lpl(:,inpl)=ival*lopt;
   cypl(:,inpl)=ival*cyopt;
   copl(:,inpl)=ival*coopt;
   gnppl(:,inpl)=ival*gnp;
   uypl(:,inpl)=ival*uyopt;
   uopl(:,inpl)=ival*uoopt;
   uythtr(inpl+1)=uyth;
end;

%We compute the transition path
for intr=2:N,
   temp=cos(acos(2*(ktr(intr)-mink)/(maxk-mink)-1)*vvec);
   ktr(intr+1)=temp*kpl(:,intr-1);
   tauktr(intr)=temp*taukpl(:,intr-1);
   totr(intr)=temp*topl(:,intr-1);
   gtr(intr)=temp*gpl(:,intr-1);
   taultr(intr)=temp*taulpl(:,intr-1);
   ltr(intr)=temp*lpl(:,intr-1);
   cytr(intr)=temp*cypl(:,intr-1);
   cotr(intr)=temp*copl(:,intr-1);
   gnptr(intr)=temp*gnppl(:,intr-1);
   uytr(intr)=temp*uypl(:,intr-1);
   uotr(intr)=temp*uopl(:,intr-1);
   uothtr(intr)=log(R*ktr(intr));
end;
   
for intr=N+1:T-1,
   temp=cos(acos(2*(ktr(intr)-mink)/(maxk-mink)-1)*vvec);
   ktr(intr+1)=temp*kpl(:,N);
   tauktr(intr)=temp*taukpl(:,N);
   totr(intr)=temp*topl(:,N);
   gtr(intr)=temp*gpl(:,N);
   taultr(intr)=temp*taulpl(:,N);
   ltr(intr)=temp*lpl(:,N);
   cytr(intr)=temp*cypl(:,N);
   cotr(intr)=temp*copl(:,N);
   gnptr(intr)=temp*gnppl(:,N);
   uytr(intr)=temp*uypl(:,N);
   uotr(intr)=temp*uopl(:,N);
   uothtr(intr)=log(R*ktr(intr));
end;

temp=cos(acos(2*(ktr(T)-mink)/(maxk-mink)-1)*vvec);
tauktr(T)=temp*taukpl(:,N);
totr(T)=temp*topl(:,N);
gtr(T)=temp*gpl(:,N);
taultr(T)=temp*taulpl(:,N);
ltr(T)=temp*lpl(:,N);
cytr(T)=temp*cypl(:,N);
cotr(T)=temp*copl(:,N);
gnptr(T)=temp*gnppl(:,N);
uytr(T)=temp*uypl(:,N);
uotr(T)=temp*uopl(:,N);
uothtr(T)=log(R*ktr(T));

n=[nn; nn(N+1)*ones(T-N-1,1)];
g2gnp=gtr./gnptr;
cy2gnp=cytr.*n./(gnptr.*(1+n));
co2gnp=cotr./(gnptr.*(1+n));
to2gnp=totr./(gnptr.*(1+n));
% Comment out the next block for the anticipated experiment
% save tnaexo;
% time=0:6;
% figure;
% subplot(2,2,1);
% plot(time,[gtr(1); gtr(1:6)]);
% set(gca,'FontSize',10);
% title(['Provision of the public good']);
%    xlabel('Time');
% 
% subplot(2,2,2);
% plot(time,[taultr(1); taultr(1:6)]);
% set(gca,'FontSize',10);
% title(['Tax rate on labor income']);
%    xlabel('Time');
% 
% subplot(2,2,3);
% plot(time,1-((R*(1-[tauktr(1); tauktr(1:6)])).^(1/30)-1)/(R^(1/30)-1));
% set(gca,'FontSize',10);
%    title(['Tax rate on capital income']);
%    xlabel('Time');
%    
% subplot(2,2,4);
% plot(time,[totr(1); totr(1:6)]);
% set(gca,'FontSize',10);
%    title(['Transfers to the old']);
%    xlabel('Time');
%    
% print -deps fig3.eps
% print -dps -append plots.ps
% 
% figure;
% %Present value of net payments
% pvp=[taultr(1)*ltr(1)+tauktr(1)*ktr(1)-totr(1)/R;...
%       taultr(1:6).*ltr(1:6)+tauktr(2:7).*ktr(2:7)-totr(2:7)/R];
% pvp=[pvp(1); pvp(1:6)];
% %Present value of the public good being provided
% pvg=[gtr(1)*(1+1/R); gtr(1:6)+gtr(2:7)/R];
% pvg=[pvg(1); pvg(1:6)];
% 
% plot(time,pvp-pvg);
% xlabel('Period in which the generation dies');
% ylabel('PV of payments less PV of benefits'); 
% 
% print -deps fig4.eps
% print -dps -append plots.ps

% Uncomment this part for the anticipated experiment
save taexo
time=0:6;

figure;
subplot(2,2,1);
plot(time,gtr(1:7));
set(gca,'FontSize',10);
title(['Provision of the public good']);
   xlabel('Time');

subplot(2,2,2);
plot(time,taultr(1:7));
set(gca,'FontSize',10);
title(['Tax rate on labor income']);
   xlabel('Time');

subplot(2,2,3);
plot(time,1-((R*(1-tauktr(1:7))).^(1/30)-1)/(R^(1/30)-1));
set(gca,'FontSize',10);
   title(['Tax rate on capital income']);
   xlabel('Time');
   
subplot(2,2,4);
plot(time,totr(1:7));
set(gca,'FontSize',10);
   title(['Transfers to the old']);
   xlabel('Time');
   
print -deps fig5.eps
print -dps -append plots.ps

figure;
%Present value of net payments
pvp=[taultr(1)*ltr(1)+tauktr(1)*ktr(1)-totr(1)/R;...
      taultr(1:7).*ltr(1:7)+tauktr(2:8).*ktr(2:8)-totr(2:8)/R];
%Present value of the public good being provided
pvg=[gtr(1)*(1+1/R); gtr(1:7)+gtr(2:8)/R];
   
plot(time,pvp(1:7)-pvg(1:7));
xlabel('Period in which the generation dies');
ylabel('PV of payments less PV of benefits'); 

print -deps fig6.eps
print -dps -append plots.ps

clear;

% Part III Solution of the stationary policy with endogenous prices
%Initialise the parameters
alph=1/3;
gam=1;
bet=.4362;
hbet=.1738;
phi=1.084;
%load phi;
sig=.6622;
n=2;
cprod=.8504;

gridk=15;
%grph controls whether graphs are displayed at each round
grph=0;
%gridpi must be equal to gridk in this version.
%We can put gridpi less than gridk using regression
%instead of interpolation, but the code needs to be changed
%because val is no longer square and cannot be inverted
gridpi=gridk;
%mink=bet/(3*gam);
mink=.01;
maxk=.1;
relaxbase=.7; 
relax=relaxbase;

% We use Chebyshev interpolation on [mink,maxk]
rts=cos((2*(1:gridk)-1)*pi/(2*gridk))';
rts=sort(rts);
% the value of each polynomial is on a column
vvec=0:gridpi-1;
val=cos(acos(rts)*vvec);
k=((maxk-mink)*rts+maxk+mink)/2;
ival=inv(val);

%Initial guess for the policy
R=(1.03^30)*ones(gridk,1);
cfnetrk=[R(1)^(1/(alph-1)); zeros(gridpi-1,1)];
cfpio=zeros(gridpi,1);
cfpig=[.001; zeros(gridpi-1,1)];

%Lower and upper bounds for the kst function
vlk=[mink+.0001; 0.0001];
vuk=[maxk-.001; 1000];

%Lower and upper bounds (should not play any role,except
%for Ty);
vlb=[mink+.001; 1e-15; 0];
vub=[maxk-.001; 200; 200];

%Initialisation
%init=[bet*ones(1,gridk)/gam+.00001; ; .001*ones(1,gridk)];
init=zeros(3,gridk);
initstar=[bet/gam; (1+bet)/gam]*ones(1,gridk);

kpopt=zeros(gridk,1);
coopt=zeros(gridk,1);
tyopt=zeros(gridk,1);
cyopt=zeros(gridk,1);
netwopt=zeros(gridk,1);
taulopt=zeros(gridk,1);
taukopt=zeros(gridk,1);
toopt=zeros(gridk,1);
gopt=zeros(gridk,1);
lopt=zeros(gridk,1);
%zx=linspace(mink,maxk);
zx=k';

epsi=1;
tic;
ite=1;
while epsi>1e-08,
   
   %We now recompute the optimal policy at the grid points
   for i=1:gridk,
      %compute the threat point
      options=optimset('Gradobj','on','TolX',1e-6,'TolFun',1e-7,'Display','off');
      [temp,fval,exitflag]=fmincon(@(tempsol) kstendo(tempsol,alph,bet,gam,n,...
          cfpio,cfnetrk,vvec,mink,maxk,k(i),cprod),initstar(:,i),[],[],...
          [],[],mink,maxk,[],options);
      if exitflag<1 | fval>1e-06,
          exitflag
          fval
          error('Failure to converge');
      end;
      initstar(:,i)=temp;
      
      if min(abs(temp-vlk))<.001,
         disp('We might have a problem here with the min');
         kstar=mink;
         lstar=1/gam;
      end;
      if max(abs(temp-vuk))<.001,
         disp('We might have a problem here with the max');
      end;
       
      kstar=temp(1);
      lstar=temp(2);
      temp=cos(vvec*acos(2*(kstar-mink)/(maxk-mink)-1));
      netrkstar=(temp*cfnetrk)^(alph-1);
      piostar=temp*cfpio;
      pigstar=temp*cfpig;
      cystar=(1-alph)*cprod*(k(i)/n)^alph*lstar^(1-alph)-kstar;
      coexstar=kstar*netrkstar+piostar;
      uyth=log(cystar)+bet*log(coexstar)-gam*lstar+hbet*phi*...
         pigstar^(1-sig)/(1-sig);
      uoth=log(cprod)+log(alph)+alph*log(k(i))+(1-alph)*log(lstar*n);
      init(:,i)=[kstar-.001; alph*cprod*k(i)^alph*(n*lstar)^(1-alph)-.002; .001];
     
      options=optimset('Gradobj','on','TolX',1e-6,'TolFun',1e-7,'Display',...
          'off','GradConstr','on');
     [temp,fval,exitflag]=fmincon(@(tempsol) nashfunendo(tempsol,alph,bet,...
         hbet,gam,n,phi,sig,cfnetrk,cfpio,cfpig,uyth,uoth,k(i),vvec,mink,...
         maxk,cprod),init(:,i),[],[],[],[],vlb,vub,@(tempsol) nashconendo(...
         tempsol,alph,bet,hbet,gam,n,phi,sig,cfnetrk,cfpio,cfpig,uyth,...
         uoth,k(i),vvec,mink,maxk,cprod),options);
     if exitflag<1,
          exitflag
          error('Failure to converge');
      end;

     kpopt(i)=temp(1);
      vf(i)=fval;
      temp2=cos(vvec*acos(2*(kpopt(i)-mink)/(maxk-mink)-1));
      netrex=(temp2*cfnetrk)^(alph-1);
      pioopt=temp2*cfpio;
      pigopt=temp2*cfpig;
      
      coopt(i)=temp(2);
      tyopt(i)=temp(3);
      copopt=netrex*kpopt(i)+pioopt;
      netwopt(i)=gam*copopt/(bet*netrex);
      cyopt(i)=netwopt(i)/gam;
      lopt(i)=(cyopt(i)+kpopt(i)-tyopt(i))/netwopt(i);
      gopt(i)=(cprod*k(i)^alph*(lopt(i)*n)^(1-alph)-(cyopt(i)+kpopt(i))*n-coopt(i))/(1+n);
      
      diffuo=log(coopt(i))+phi*(hbet/bet)*gopt(i)^(1-sig)/(1-sig)-uoth;
      taukopt(i)=1-coopt(i)/(alph*cprod*k(i)^alph*(lopt(i)*n)^(1-alph)*(1+diffuo));
      toopt(i)=coopt(i)*diffuo/(1+diffuo);
      R(i)=alph*cprod*k(i)^(alph-1)*(lopt(i)*n).^(1-alph);
      init(:,i)=temp;
   end;
   cfnetrk2=ival*((1-taukopt).*R).^(1/(alph-1));
   cfpig2=ival*gopt;
   cfpio2=ival*toopt;
   cfR=ival*R;
   epsi2=max(max(abs([cfnetrk cfpig cfpio]-[cfnetrk2 cfpig2 cfpio2])))
   ite
   if epsi2<epsi,
       taukoptold=taukopt;
       goptold=gopt;
       tooptold=toopt;
       cfnetrkold=cfnetrk;
       cfpigold=cfpig;
       cfpioold=cfpio;
       cfnetrk=relax*cfnetrk2+(1-relax)*cfnetrk;
       cfpig=relax*cfpig2+(1-relax)*cfpig;
       cfpio=relax*cfpio2+(1-relax)*cfpio;
       epsi=epsi2;
       ite=ite+1;
       relax=relaxbase;
   else,
       relax=relax/2;
       cfnetrk=relax*cfnetrk2+(1-relax)*cfnetrkold;
       cfpig=relax*cfpig2+(1-relax)*cfpigold;
       cfpio=relax*cfpio2+(1-relax)*cfpioold;
       if relax<.0001,
           error('Failure to converge');
       end;
   end;
   if grph==1,
      subplot(2,2,1);
      temp=cos(acos(2*(zx-mink)/(maxk-mink)-1)'*vvec);
      plot(zx,(temp*cfnetrk).^(alph-1),'-',zx,(temp*cfnetrk2).^(alph-1),':');
      title('netrex');
      subplot(2,2,2);
      plot(zx,temp*cfpig,zx,temp*cfpig2,':');
      title('gopt');
      subplot(2,2,3);
      plot(zx,temp*cfpio,zx,temp*cfpio2,':');
      title('toopt');
      subplot(2,2,4);
      %      plot(0,10000*epsi,'+');
      %plot(zx,temp*cfR,'-');
      plot(zx,temp*cfnetrk,zx,temp*cfnetrk2,':');
      pause;
   end;
   toc;
end;

%Compute variables not required in the iterations
temp=cos(acos(2*(kpopt-mink)/(maxk-mink)-1)*vvec);
netrkopt=(temp*cfnetrk).^(alph-1);
pioopt=temp*cfpio;
pigopt=temp*cfpig;
taulopt=1-netwopt./((1-alph)*cprod*(k./(n*lopt)).^alph);
copopt=netrkopt.*kpopt+pioopt;
uyopt=log(cyopt)+bet*log(copopt)-gam*lopt+phi*gopt.^(1-sig)/(1-sig)+...
   hbet*phi*pigopt.^(1-sig)/(1-sig);
uoopt=log(coopt)+phi*(hbet/bet)*gopt.^(1-sig)/(1-sig);
gnp=cyopt*n/(1+n)+coopt/(1+n)+gopt+kpopt*n/(1+n);
gnp2=cprod*k.^alph.*(lopt*n).^(1-alph)/(1+n);
Ropt=(alph*cprod*k.^(alph-1).*(lopt*n).^(1-alph));
wopt=(1-alph)*cprod*k.^(alph).*(lopt*n).^(-alph);
g2gnp=gopt./gnp;
cy2gnp=cyopt*n./(gnp*(1+n));
co2gnp=coopt./(gnp*(1+n));
to2gnp=toopt./(gnp*(1+n));

%Find the steady state
options=optimset('FunValCheck','on','TolX',1e-16);
cfkp=ival*kpopt;
kss=fzero(@(ksssolve) steady(ksssolve,mink,maxk,vvec,cfkp),[mink+.0001;...
    maxk-.0001],options);
temp=cos(acos(2*(kss-mink)/(maxk-mink)-1)*vvec);
netrkss=(temp*cfnetrk)^(alph-1);
toss=temp*cfpio;
gss=temp*cfpig;
temp2=temp*ival;
taukss=temp2*taukopt;
taulss=temp2*taulopt;
lss=temp2*lopt;
cyss=temp2*cyopt;
coss=temp2*coopt;
gnpss=temp2*gnp;
g2gnpss=gss/gnpss;
cy2gnpss=cyss*n/(gnpss*(1+n));
co2gnpss=coss/(gnpss*(1+n));
to2gnpss=toss/(gnpss*(1+n));
Rss=alph*cprod*kss^(alph-1)*(lss*n)^(1-alph);
taukhat=1-((Rss*(1-taukss))^(1/30)-1)/(Rss^(1/30)-1);
oldss=Rss*kss-coss;
youngss=(1-alph)*cprod*kss^alph*n^(-alph)*lss^(1-alph)-cyss-kss;
ratss=oldss/(gss*(1+n));
save olgpendo;
subplot(2,2,1);
plot(k,g2gnp);
set(gca,'FontSize',10);
title(['Provision of the public good as a fraction of GNP']);
   xlabel('Capital per old person');

subplot(2,2,2);
plot(k,taulopt);
set(gca,'FontSize',10);
title(['Tax rate on labor income']);
   xlabel('Capital per old person');

subplot(2,2,3);
plot(k,1-((Ropt.*(1-taukopt)).^(1/30)-1)./(Ropt.^(1/30)-1));
set(gca,'FontSize',10);
   title(['Tax rate on capital income']);
   xlabel('Capital per old person');
   
subplot(2,2,4);
plot(k,to2gnp);
set(gca,'FontSize',10);
   title(['Transfers to the old as a fraction of GNP']);
   xlabel('Capital per old person');
      
print -deps fig7.eps;
print -dps -append plots.ps;

save olgpendo;
clear;

% Part IV: computation of the transition with endogenous prices. This
% replicates most of part III

%This code computes the effects of a
%demographic transition with endogenous factor prices
%Initialise the parameters
gam=1;
alph=1/3;
bet=.4362;
hbet=.1738;
phi=1.084;
sig=.6622;
cprod=.8504;
nn=[2; 11/9]; % Unanticipated experiment
% nn=[2; 2; 11/9]; % Anticipated experiment
%N is the number of policy rules we have to keep track of
N=size(nn,1)-1;
%T is the number of transition periods I consider
T=7;
gridk=15;
%grph controls whether graphs are displayed at each round
grph=0;
%gridpi must be equal to gridk in this version.
%We can put gridpi less than gridk using regression
%instead of interpolation, but the code needs to be changed
%because val is no longer square and cannot be inverted
gridpi=gridk;
mink=0.01;
maxk=.1;
relaxbase=.7; 
relax=relaxbase;

% We use Chebyshev interpolation on [mink,maxk]
rts=cos((2*(1:gridk)-1)*pi/(2*gridk))';
rts=sort(rts);
% the value of each polynomial is on a column
vvec=0:gridpi-1;
val=cos(acos(rts)*vvec);
k=((maxk-mink)*rts+maxk+mink)/2;
ival=inv(val);

%Initial guess for the policy
R=(1.03^30)*ones(gridk,1);
cfnetrk=[R(1)^(1/(alph-1)); zeros(gridpi-1,1)];
cfpio=zeros(gridpi,1);
cfpig=[.001; zeros(gridpi-1,1)];

%Lower and upper bounds for the kst function
vlk=[mink+.0001; 0.0001];
vuk=[maxk-.001; 1000];

%Lower and upper bounds (should not play any role,except
%for Ty);
vlb=[mink+.001; 1e-15; 0];
vub=[maxk-.001; 200; 200];

%Initialisation
%init=[bet*ones(1,gridk)/gam+.00001; ; .001*ones(1,gridk)];
init=zeros(3,gridk);
initstar=[bet/gam; (1+bet)/gam]*ones(1,gridk);

kpopt=zeros(gridk,1);
coopt=zeros(gridk,1);
tyopt=zeros(gridk,1);
cyopt=zeros(gridk,1);
netwopt=zeros(gridk,1);
taulopt=zeros(gridk,1);
taukopt=zeros(gridk,1);
toopt=zeros(gridk,1);
gopt=zeros(gridk,1);
lopt=zeros(gridk,1);
%zx=linspace(mink,maxk);
zx=k';

ktr=zeros(T,1);
tauktr=zeros(T,1);
totr=zeros(T,1);
gtr=zeros(T,1);
taultr=zeros(T,1);
ltr=zeros(T,1);
Rtr=zeros(T,1);
cytr=zeros(T,1);
cotr=zeros(T,1);
gnptr=zeros(T,1);
uytr=zeros(T,1);
uotr=zeros(T,1);
uythtr=zeros(T,1);
uothtr=zeros(T,1);


cfnetrkpl=zeros(gridpi,N);
kpl=zeros(gridpi,N);
taukpl=zeros(gridpi,N);
topl=zeros(gridpi,N);
gpl=zeros(gridpi,N);
taulpl=zeros(gridpi,N);
lpl=zeros(gridpi,N);
cypl=zeros(gridpi,N);
copl=zeros(gridpi,N);
gnppl=zeros(gridpi,N);
uypl=zeros(gridpi,N);
uopl=zeros(gridpi,N);

%First, we compute the initial steady-state
n=nn(1);
epsi=1;
tic;
ite=1;
while epsi>1e-06,
   %We now recompute the optimal policy at the grid points
   for i=1:gridk,
      %compute the threat point
      options=optimset('Gradobj','on','TolX',1e-6,'TolFun',1e-7,'Display','off');
      [temp,fval,exitflag]=fmincon(@(tempsol) kstendo(tempsol,alph,bet,gam,n,...
          cfpio,cfnetrk,vvec,mink,maxk,k(i),cprod),initstar(:,i),[],[],...
          [],[],mink,maxk,[],options);
      if exitflag<1 | fval>1e-06,
          exitflag
          fval
          error('Failure to converge');
      end;
      initstar(:,i)=temp;

      if min(abs(temp-vlk))<.001,
         disp('We might have a problem here with the min');
         kstar=mink;
         lstar=1/gam;
      end;
      if max(abs(temp-vuk))<.001,
         disp('We might have a problem here with the max');
      end;
      
      kstar=temp(1);
      lstar=temp(2);
      
      temp=cos(vvec*acos(2*(kstar-mink)/(maxk-mink)-1));
      netrkstar=(temp*cfnetrk)^(alph-1);
      piostar=temp*cfpio;
      pigstar=temp*cfpig;
      cystar=(1-alph)*cprod*(k(i)/n)^alph*lstar^(1-alph)-kstar;
      coexstar=kstar*netrkstar+piostar;
      uyth=log(cystar)+bet*log(coexstar)-gam*lstar+hbet*phi*...
         pigstar^(1-sig)/(1-sig);
      uoth=log(cprod)+log(alph)+alph*log(k(i))+(1-alph)*log(lstar*n);
      %if init(:,i)==zeros(3,1),
         init(:,i)=[kstar-.001; alph*cprod*k(i)^alph*(n*lstar)^(1-alph)-.002; .001];
      %end;
      
     options=optimset('Gradobj','on','TolX',1e-6,'TolFun',1e-7,'Display',...
          'off','GradConstr','on');
     [temp,fval,exitflag]=fmincon(@(tempsol) nashfunendo(tempsol,alph,bet,...
         hbet,gam,n,phi,sig,cfnetrk,cfpio,cfpig,uyth,uoth,k(i),vvec,mink,...
         maxk,cprod),init(:,i),[],[],[],[],vlb,vub,@(tempsol) nashconendo(...
         tempsol,alph,bet,hbet,gam,n,phi,sig,cfnetrk,cfpio,cfpig,uyth,...
         uoth,k(i),vvec,mink,maxk,cprod),options);
     if exitflag<1,
          exitflag
          error('Failure to converge');
      end;

      kpopt(i)=temp(1);
      vf(i)=fval;
      temp2=cos(vvec*acos(2*(kpopt(i)-mink)/(maxk-mink)-1));
      netrex=(temp2*cfnetrk)^(alph-1);
      pioopt=temp2*cfpio;
      pigopt=temp2*cfpig;
      
      %i
      %pause;
      coopt(i)=temp(2);
      tyopt(i)=temp(3);
      copopt=netrex*kpopt(i)+pioopt;
      netwopt(i)=gam*copopt/(bet*netrex);
      cyopt(i)=netwopt(i)/gam;
      lopt(i)=(cyopt(i)+kpopt(i)-tyopt(i))/netwopt(i);
      gopt(i)=(cprod*k(i)^alph*(lopt(i)*n)^(1-alph)-(cyopt(i)+kpopt(i))*n-coopt(i))/(1+n);
      
      diffuo=log(coopt(i))+phi*(hbet/bet)*gopt(i)^(1-sig)/(1-sig)-uoth;
      taukopt(i)=1-coopt(i)/(alph*cprod*k(i)^alph*(lopt(i)*n)^(1-alph)*(1+diffuo));
      toopt(i)=coopt(i)*diffuo/(1+diffuo);
      R(i)=alph*cprod*k(i)^(alph-1)*(lopt(i)*n).^(1-alph);
      init(:,i)=temp;
   end;
   cfnetrk2=ival*((1-taukopt).*R).^(1/(alph-1));
   cfpig2=ival*gopt;
   cfpio2=ival*toopt;
   cfR=ival*R;
   epsi2=max(max(abs([cfnetrk cfpig cfpio]-[cfnetrk2 cfpig2 cfpio2])))
   ite
   if epsi2<epsi,
       taukoptold=taukopt;
       goptold=gopt;
       tooptold=toopt;
       cfnetrkold=cfnetrk;
       cfpigold=cfpig;
       cfpioold=cfpio;
       cfnetrk=relax*cfnetrk2+(1-relax)*cfnetrk;
       cfpig=relax*cfpig2+(1-relax)*cfpig;
       cfpio=relax*cfpio2+(1-relax)*cfpio;
       epsi=epsi2;
       ite=ite+1;
       relax=relaxbase;
   else,
       relax=relax/2;
       cfnetrk=relax*cfnetrk2+(1-relax)*cfnetrkold;
       cfpig=relax*cfpig2+(1-relax)*cfpigold;
       cfpio=relax*cfpio2+(1-relax)*cfpioold;
       if relax<.0001,
           error('Failure to converge');
       end;
   end;
   if grph==1,
      subplot(2,2,1);
      temp=cos(acos(2*(zx-mink)/(maxk-mink)-1)'*vvec);
      plot(zx,(temp*cfnetrk).^(alph-1),'-',zx,(temp*cfnetrk2).^(alph-1),':');
      title('netrex');
      subplot(2,2,2);
      plot(zx,temp*cfpig,zx,temp*cfpig2,':');
      title('gopt');
      subplot(2,2,3);
      plot(zx,temp*cfpio,zx,temp*cfpio2,':');
      title('toopt');
      subplot(2,2,4);
      %      plot(0,10000*epsi,'+');
      %plot(zx,temp*cfR,'-');
      plot(zx,temp*cfnetrk,zx,temp*cfnetrk2,':');
      pause;
   end;
   %toc;
end;

%Compute variables not required in the iterations
temp=cos(acos(2*(kpopt-mink)/(maxk-mink)-1)*vvec);
netrkopt=(temp*cfnetrk).^(alph-1);
pioopt=temp*cfpio;
pigopt=temp*cfpig;
taulopt=1-netwopt./((1-alph)*cprod*(k./(n*lopt)).^alph);
copopt=netrkopt.*kpopt+pioopt;
uyopt=log(cyopt)+bet*log(copopt)-gam*lopt+phi*gopt.^(1-sig)/(1-sig)+...
   hbet*phi*pigopt.^(1-sig)/(1-sig);
uoopt=log(coopt)+phi*(hbet/bet)*gopt.^(1-sig)/(1-sig);
gnp=cyopt*n/(1+n)+coopt/(1+n)+gopt+kpopt*n/(1+n);
gnp2=cprod*k.^alph.*(lopt*n).^(1-alph)/(1+n);

%Find the steady state
options=optimset('FunValCheck','on','TolX',1e-16);
cfkp=ival*kpopt;
ktr(1)=fzero(@(ksssolve) steady(ksssolve,mink,maxk,vvec,cfkp),[mink+.0001;...
    maxk-.0001],options);
temp=cos(acos(2*(ktr(1)-mink)/(maxk-mink)-1)*vvec);
ktr(2)=ktr(1);
totr(1)=temp*cfpio;
gtr(1)=temp*cfpig;
temp2=temp*ival;
tauktr(1)=temp2*taukopt;
taultr(1)=temp2*taulopt;
ltr(1)=temp2*lopt;
cytr(1)=temp2*cyopt;
cotr(1)=temp2*coopt;
gnptr(1)=temp2*gnp;

%We now compute the policy rules in the final steady state

n=nn(N+1);
epsi=1;
tic;
ite=1;
while epsi>1e-06,
   
   %We now recompute the optimal policy at the grid points
   for i=1:gridk,
      %compute the threat point
      options=optimset('Gradobj','on','TolX',1e-6,'TolFun',1e-7,'Display','off');
      [temp,fval,exitflag]=fmincon(@(tempsol) kstendo(tempsol,alph,bet,gam,n,...
          cfpio,cfnetrk,vvec,mink,maxk,k(i),cprod),initstar(:,i),[],[],...
          [],[],mink,maxk,[],options);
      if exitflag<1 | fval>1e-06,
          exitflag
          fval
          error('Failure to converge');
      end;
      initstar(:,i)=temp;
      
       if min(abs(temp-vlk))<.001,
         disp('We might have a problem here with the min');
         kstar=mink;
         lstar=1/gam;
      end;
      if max(abs(temp-vuk))<.001,
         disp('We might have a problem here with the max');
      end;
      
      kstar=temp(1);
      lstar=temp(2);
  
      temp=cos(vvec*acos(2*(kstar-mink)/(maxk-mink)-1));
      netrkstar=(temp*cfnetrk)^(alph-1);
      piostar=temp*cfpio;
      pigstar=temp*cfpig;
      cystar=(1-alph)*cprod*(k(i)/n)^alph*lstar^(1-alph)-kstar;
      coexstar=kstar*netrkstar+piostar;
      uyth=log(cystar)+bet*log(coexstar)-gam*lstar+hbet*phi*...
         pigstar^(1-sig)/(1-sig);
      uoth=log(cprod)+log(alph)+alph*log(k(i))+(1-alph)*log(lstar*n);
      %if init(:,i)==zeros(3,1),
         init(:,i)=[kstar-.001; alph*cprod*k(i)^alph*(n*lstar)^(1-alph)-.002; .001];
      %end;
      
     options=optimset('Gradobj','on','TolX',1e-6,'TolFun',1e-7,'Display',...
          'off','GradConstr','on');
     [temp,fval,exitflag]=fmincon(@(tempsol) nashfunendo(tempsol,alph,bet,...
         hbet,gam,n,phi,sig,cfnetrk,cfpio,cfpig,uyth,uoth,k(i),vvec,mink,...
         maxk,cprod),init(:,i),[],[],[],[],vlb,vub,@(tempsol) nashconendo(...
         tempsol,alph,bet,hbet,gam,n,phi,sig,cfnetrk,cfpio,cfpig,uyth,...
         uoth,k(i),vvec,mink,maxk,cprod),options);
     if exitflag<1,
          exitflag
          error('Failure to converge');
      end;

      kpopt(i)=temp(1);
      vf(i)=fval;
      temp2=cos(vvec*acos(2*(kpopt(i)-mink)/(maxk-mink)-1));
      netrex=(temp2*cfnetrk)^(alph-1);
      pioopt=temp2*cfpio;
      pigopt=temp2*cfpig;
      
      %i
      %pause;
      coopt(i)=temp(2);
      tyopt(i)=temp(3);
      copopt=netrex*kpopt(i)+pioopt;
      netwopt(i)=gam*copopt/(bet*netrex);
      cyopt(i)=netwopt(i)/gam;
      lopt(i)=(cyopt(i)+kpopt(i)-tyopt(i))/netwopt(i);
      gopt(i)=(cprod*k(i)^alph*(lopt(i)*n)^(1-alph)-(cyopt(i)+kpopt(i))*n-coopt(i))/(1+n);
      
      diffuo=log(coopt(i))+phi*(hbet/bet)*gopt(i)^(1-sig)/(1-sig)-uoth;
      taukopt(i)=1-coopt(i)/(alph*cprod*k(i)^alph*(lopt(i)*n)^(1-alph)*(1+diffuo));
      toopt(i)=coopt(i)*diffuo/(1+diffuo);
      R(i)=alph*cprod*k(i)^(alph-1)*(lopt(i)*n).^(1-alph);
      init(:,i)=temp;
   end;
   cfnetrk2=ival*((1-taukopt).*R).^(1/(alph-1));
   cfpig2=ival*gopt;
   cfpio2=ival*toopt;
   cfR=ival*R;
   epsi2=max(max(abs([cfnetrk cfpig cfpio]-[cfnetrk2 cfpig2 cfpio2])))
   ite
   if epsi2<epsi,
       taukoptold=taukopt;
       goptold=gopt;
       tooptold=toopt;
       cfnetrkold=cfnetrk;
       cfpigold=cfpig;
       cfpioold=cfpio;
       cfnetrk=relax*cfnetrk2+(1-relax)*cfnetrk;
       cfpig=relax*cfpig2+(1-relax)*cfpig;
       cfpio=relax*cfpio2+(1-relax)*cfpio;
       epsi=epsi2;
       ite=ite+1;
       relax=relaxbase;
   else,
       relax=relax/2;
       cfpik=relax*cfnetrk2+(1-relax)*cfnetrkold;
       cfpig=relax*cfpig2+(1-relax)*cfpigold;
       cfpio=relax*cfpio2+(1-relax)*cfpioold;
       if relax<.0001,
           error('Failure to converge');
       end;
   end;
   if grph==1,
      subplot(2,2,1);
      temp=cos(acos(2*(zx-mink)/(maxk-mink)-1)'*vvec);
      plot(zx,(temp*cfnetrk).^(alph-1),'-',zx,(temp*cfnetrk2).^(alph-1),':');
      title('netrex');
      subplot(2,2,2);
      plot(zx,temp*cfpig,zx,temp*cfpig2,':');
      title('gopt');
      subplot(2,2,3);
      plot(zx,temp*cfpio,zx,temp*cfpio2,':');
      title('toopt');
      subplot(2,2,4);
      %      plot(0,10000*epsi,'+');
      %plot(zx,temp*cfR,'-');
      plot(zx,temp*cfnetrk,zx,temp*cfnetrk2,':');
      pause;
   end;
   %toc;
end;

%Compute variables not required in the iterations
temp=cos(acos(2*(kpopt-mink)/(maxk-mink)-1)*vvec);
netrkopt=(temp*cfnetrk).^(alph-1);
pioopt=temp*cfpio;
pigopt=temp*cfpig;
taulopt=1-netwopt./((1-alph)*cprod*(k./(n*lopt)).^alph);
copopt=netrkopt.*kpopt+pioopt;
uyopt=log(cyopt)+bet*log(copopt)-gam*lopt+phi*gopt.^(1-sig)/(1-sig)+...
   hbet*phi*pigopt.^(1-sig)/(1-sig);
uoopt=log(coopt)+phi*(hbet/bet)*gopt.^(1-sig)/(1-sig);
gnp=cyopt*n/(1+n)+coopt/(1+n)+gopt+kpopt*n/(1+n);
gnp2=cprod*k.^alph.*(lopt*n).^(1-alph)/(1+n);

cfnetrkpl(:,N)=cfnetrk;
kpl(:,N)=ival*kpopt;
taukpl(:,N)=ival*taukopt;
topl(:,N)=cfpio;
gpl(:,N)=cfpig;
taulpl(:,N)=ival*taulopt;
lpl(:,N)=ival*lopt;
cypl(:,N)=ival*cyopt;
copl(:,N)=ival*coopt;
gnppl(:,N)=ival*gnp;
uypl(:,N)=ival*uyopt;
uopl(:,N)=ival*uoopt;

%Now we find the transition policy rules by backward induction

for inpl=N-1:-1:1,
   n=nn(inpl+1);
   
   %We now recompute the optimal policy at the grid points
   for i=1:gridk,
      %compute the threat point
      options=optimset('Gradobj','on','TolX',1e-6,'TolFun',1e-7,'Display','off');
      [temp,fval,exitflag]=fmincon(@(tempsol) kstendo(tempsol,alph,bet,gam,n,...
          cfpio,cfnetrk,vvec,mink,maxk,k(i),cprod),initstar(:,i),[],[],...
          [],[],mink,maxk,[],options);
      if exitflag<1 | fval>1e-06,
          exitflag
          fval
          error('Failure to converge');
      end;
      initstar(:,i)=temp;
      
      if min(abs(temp-vlk))<.001,
         disp('We might have a problem here with the min');
         kstar=mink;
         lstar=1/gam;
      end;
      if max(abs(temp-vuk))<.001,
         disp('We might have a problem here with the max');
      end;
      
      kstar=temp(1);
      lstar=temp(2);
      
      temp=cos(vvec*acos(2*(kstar-mink)/(maxk-mink)-1));
      netrkstar=(temp*cfnetrk)^(alph-1);
      piostar=temp*cfpio;
      pigstar=temp*cfpig;
      cystar=(1-alph)*cprod*(k(i)/n)^alph*lstar^(1-alph)-kstar;
      coexstar=kstar*netrkstar+piostar;
      uyth=log(cystar)+bet*log(coexstar)-gam*lstar+hbet*phi*...
         pigstar^(1-sig)/(1-sig);
      uoth=log(cprod)+log(alph)+alph*log(k(i))+(1-alph)*log(lstar*n);
      %if init(:,i)==zeros(3,1),
         init(:,i)=[kstar-.001; alph*cprod*k(i)^alph*(n*lstar)^(1-alph)-.002; .001];
      %end;
      
     options=optimset('Gradobj','on','TolX',1e-6,'TolFun',1e-7,'Display',...
          'off','GradConstr','on');
     [temp,fval,exitflag]=fmincon(@(tempsol) nashfunendo(tempsol,alph,bet,...
         hbet,gam,n,phi,sig,cfnetrk,cfpio,cfpig,uyth,uoth,k(i),vvec,mink,...
         maxk,cprod),init(:,i),[],[],[],[],vlb,vub,@(tempsol) nashconendo(...
         tempsol,alph,bet,hbet,gam,n,phi,sig,cfnetrk,cfpio,cfpig,uyth,...
         uoth,k(i),vvec,mink,maxk,cprod),options);
     if exitflag<1,
          exitflag
          error('Failure to converge');
      end;

      kpopt(i)=temp(1);
      vf(i)=fval;
      temp2=cos(vvec*acos(2*(kpopt(i)-mink)/(maxk-mink)-1));
      netrex=(temp2*cfnetrk)^(alph-1);
      pioopt=temp2*cfpio;
      pigopt=temp2*cfpig;
      
      %i
      %pause;
      coopt(i)=temp(2);
      tyopt(i)=temp(3);
      copopt=netrex*kpopt(i)+pioopt;
      netwopt(i)=gam*copopt/(bet*netrex);
      cyopt(i)=netwopt(i)/gam;
      lopt(i)=(cyopt(i)+kpopt(i)-tyopt(i))/netwopt(i);
      gopt(i)=(cprod*k(i)^alph*(lopt(i)*n)^(1-alph)-(cyopt(i)+kpopt(i))*n-coopt(i))/(1+n);
      
      diffuo=log(coopt(i))+phi*(hbet/bet)*gopt(i)^(1-sig)/(1-sig)-uoth;
      taukopt(i)=1-coopt(i)/(alph*cprod*k(i)^alph*(lopt(i)*n)^(1-alph)*(1+diffuo));
      toopt(i)=coopt(i)*diffuo/(1+diffuo);
      R(i)=alph*cprod*k(i)^(alph-1)*(lopt(i)*n).^(1-alph);
      init(:,i)=temp;
   end;
   cfnetrk=ival*((1-taukopt).*R).^(1/(alph-1));
   cfpig=ival*gopt;
   cfpio=ival*toopt;

%   inpl
   toc;
   
   %Compute variables not required in the iterations
   temp=cos(acos(2*(kpopt-mink)/(maxk-mink)-1)*vvec);
   netrkopt=(temp*cfnetrk).^(alph-1);
   pioopt=temp*cfpio;
   pigopt=temp*cfpig;
   taulopt=1-netwopt./((1-alph)*cprod*(k./(n*lopt)).^alph);
   copopt=netrkopt.*kpopt+pioopt;
   uyopt=log(cyopt)+bet*log(copopt)-gam*lopt+phi*gopt.^(1-sig)/(1-sig)+...
      hbet*phi*pigopt.^(1-sig)/(1-sig);
   uoopt=log(coopt)+phi*(hbet/bet)*gopt.^(1-sig)/(1-sig);
   gnp=cyopt*n/(1+n)+coopt/(1+n)+gopt+kpopt*n/(1+n);
   gnp2=cprod*k.^alph.*(lopt*n).^(1-alph)/(1+n);
   
   cfnetrkpl(:,inpl)=cfnetrk;
   kpl(:,inpl)=ival*kpopt;
   taukpl(:,inpl)=ival*taukopt;
   topl(:,inpl)=cfpio;
   gpl(:,inpl)=cfpig;
   taulpl(:,inpl)=ival*taulopt;
   lpl(:,inpl)=ival*lopt;
   cypl(:,inpl)=ival*cyopt;
   copl(:,inpl)=ival*coopt;
   gnppl(:,inpl)=ival*gnp;
   uypl(:,inpl)=ival*uyopt;
   uopl(:,inpl)=ival*uoopt;
end;

%We compute the transition path
for intr=2:N,
   temp=cos(acos(2*(ktr(intr)-mink)/(maxk-mink)-1)*vvec);
   ktr(intr+1)=temp*kpl(:,intr-1);
   tauktr(intr)=temp*taukpl(:,intr-1);
   totr(intr)=temp*topl(:,intr-1);
   gtr(intr)=temp*gpl(:,intr-1);
   taultr(intr)=temp*taulpl(:,intr-1);
   ltr(intr)=temp*lpl(:,intr-1);
   cytr(intr)=temp*cypl(:,intr-1);
   cotr(intr)=temp*copl(:,intr-1);
   gnptr(intr)=temp*gnppl(:,intr-1);
   uytr(intr)=temp*uypl(:,intr-1);
   uotr(intr)=temp*uopl(:,intr-1);
end;

   
for intr=N+1:T-1,
   temp=cos(acos(2*(ktr(intr)-mink)/(maxk-mink)-1)*vvec);
   ktr(intr+1)=temp*kpl(:,N);
   tauktr(intr)=temp*taukpl(:,N);
   totr(intr)=temp*topl(:,N);
   gtr(intr)=temp*gpl(:,N);
   taultr(intr)=temp*taulpl(:,N);
   ltr(intr)=temp*lpl(:,N);
   cytr(intr)=temp*cypl(:,N);
   cotr(intr)=temp*copl(:,N);
   gnptr(intr)=temp*gnppl(:,N);
   uytr(intr)=temp*uypl(:,N);
   uotr(intr)=temp*uopl(:,N);
end;

temp=cos(acos(2*(ktr(T)-mink)/(maxk-mink)-1)*vvec);
tauktr(T)=temp*taukpl(:,N);
totr(T)=temp*topl(:,N);
gtr(T)=temp*gpl(:,N);
taultr(T)=temp*taulpl(:,N);
ltr(T)=temp*lpl(:,N);
cytr(T)=temp*cypl(:,N);
cotr(T)=temp*copl(:,N);
gnptr(T)=temp*gnppl(:,N);
uytr(T)=temp*uypl(:,N);
uotr(T)=temp*uopl(:,N);

n=[nn; nn(N+1)*ones(T-N-1,1)];
g2gnptr=gtr./gnptr;
cy2gnptr=cytr.*n./(gnptr.*(1+n));
co2gnptr=cotr./(gnptr.*(1+n));
to2gnptr=totr./(gnptr.*(1+n));
Rtr=alph*cprod*ktr.^(alph-1).*(ltr.*n).^(1-alph);
wtr=(1-alph)*cprod*ktr.^(alph).*(ltr.*n).^(-alph);
taukhat=1-((Rtr.*(1-tauktr)).^(1/30)-1)./(Rtr.^(1/30)-1);
oldtr=Rtr.*ktr-cotr;
youngtr=(1-alph)*cprod*ktr.^alph.*n.^(-alph).*ltr.^(1-alph)-cytr-ktr;
rattr=oldtr./(gtr.*(1+n));

% This part should be commented for the anticipated transition
time=0:6;

figure;
subplot(2,2,1);
plot(time,[gtr(1); gtr(1:6)]);
set(gca,'FontSize',10);
title(['Provision of the public good']);
   xlabel('Time');

subplot(2,2,2);
plot(time,[taultr(1); taultr(1:6)]);
set(gca,'FontSize',10);
title(['Tax rate on labor income']);
   xlabel('Time');

subplot(2,2,3);
plot(time,1-(([Rtr(1); Rtr(1:6)].*(1-[tauktr(1); tauktr(1:6)])).^...
    (1/30)-1)./([Rtr(1); Rtr(1:6)].^(1/30)-1));
set(gca,'FontSize',10);
   title(['Tax rate on capital income']);
   xlabel('Time');
   
subplot(2,2,4);
plot(time,[totr(1); totr(1:6)]);
set(gca,'FontSize',10);
   title(['Transfers to the old']);
   xlabel('Time');
   
print -deps fig8.eps;
print -dps -append plots.ps;

save tnaendo;

% This part should be uncommented for the anticipated transition
% time=0:6;
% 
% figure;
% subplot(2,2,1);
% plot(time,gtr);
% set(gca,'FontSize',10);
% title(['Provision of the public good']);
%    xlabel('Time');
% 
% subplot(2,2,2);
% plot(time,taultr);
% set(gca,'FontSize',10);
% title(['Tax rate on labor income']);
%    xlabel('Time');
% 
% subplot(2,2,3);
% plot(time,1-((Rtr.*(1-tauktr)).^(1/30)-1)./(Rtr.^(1/30)-1));
% set(gca,'FontSize',10);
%    title(['Tax rate on capital income']);
%    xlabel('Time');
%    
% subplot(2,2,4);
% plot(time,totr);
% set(gca,'FontSize',10);
%    title(['Transfers to the old']);
%    xlabel('Time');
%    
% print -deps fig9.eps;
% print -dps -append plots.ps;
% 
% save taendo;


