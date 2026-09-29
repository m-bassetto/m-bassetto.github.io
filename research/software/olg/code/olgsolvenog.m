% This code computes the solution of the general equilibrium economy
% with no government spending
% Copyright by Marco Bassetto, 1998-2007. This code can be freely
% distributed and modified for research purposes only, 
% provided this copyright notice is included in the modified code. 
% Proper credit should be given in all publications arising from
% modifications of this code; this should include a citation of 
% "Political Economy of Taxation in an Overlapping-Generations Economy,"
% by Marco Bassetto

warning off optim:fmincon:SwitchingToMediumScale

% Part III Solution of the stationary policy with endogenous prices
%Initialise the parameters
alph=1/3;
gam=1;
bet=.4362;
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
% Since general equilibrium effects are very important
% a constant guess is a bad idea, as it leads to no government
% in the first round, with a number of binding constraints that
% damage numerical precision
% We start instead from the laissez faire equilibrium
R=alph*cprod*k.^(alph-1)*n^(1-alph)*((1+bet)/gam)^(1-alph);
cfnetrk=ival*R.^(1/(alph-1));
cfpio=zeros(gridpi,1);

%Lower and upper bounds for the kst function
vlk=[mink+.0001; 0.0001];
vuk=[maxk-.001; 1000];

%Lower and upper bounds (should not play any role,except
%for Ty);
vlb=[mink+.001; 0];
vub=[maxk-.001; 200];

%Initialisation
%init=[bet*ones(1,gridk)/gam+.00001; ; .001*ones(1,gridk)];
init=zeros(2,gridk);
initstar=[bet/gam; (1+bet)/gam]*ones(1,gridk);

kpopt=zeros(gridk,1);
coopt=zeros(gridk,1);
tyopt=zeros(gridk,1);
cyopt=zeros(gridk,1);
netwopt=zeros(gridk,1);
taulopt=zeros(gridk,1);
taukopt=zeros(gridk,1);
toopt=zeros(gridk,1);
lopt=zeros(gridk,1);
%zx=linspace(mink,maxk);
zx=k';

epsi=1000000;
tic;
ite=1;
while epsi>1e-08,
   
   %We now recompute the optimal policy at the grid points
   for i=1:gridk,
      %compute the threat point
      options=optimset('Gradobj','on','TolX',1e-8,'TolFun',1e-7,'Display','off');
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
      cystar=(1-alph)*cprod*(k(i)/n)^alph*lstar^(1-alph)-kstar;
      coexstar=kstar*netrkstar+piostar;
      uyth=log(cystar)+bet*log(coexstar)-gam*lstar;
      uoth=log(cprod)+log(alph)+alph*log(k(i))+(1-alph)*log(lstar*n);
      init(:,i)=[kstar-.01; .01];
     
      options=optimset('Gradobj','on','TolX',1e-10,'TolFun',1e-12,'Display',...
          'off','GradConstr','on');
     [temp,fval,exitflag]=fmincon(@(tempsol) nashfunendonog(tempsol,alph,bet,...
        gam,n,cfnetrk,cfpio,uyth,uoth,k(i),vvec,mink,...
         maxk,cprod),init(:,i),[],[],[],[],vlb,vub,@(tempsol) nashconendonog(...
         tempsol,alph,bet,gam,n,cfnetrk,cfpio,uyth,...
         uoth,k(i),vvec,mink,maxk,cprod),options);
     % We repeat the optimization as sometimes the algorithm gets
     % stuck in the first round because it carries over improper Hessian
     % information (the Hessian behaves badly close to the boundaries of
     % the constraints). The fresh start solves this problem.
     [temp,fval,exitflag]=fmincon(@(tempsol) nashfunendonog(tempsol,alph,bet,...
        gam,n,cfnetrk,cfpio,uyth,uoth,k(i),vvec,mink,...
         maxk,cprod),temp,[],[],[],[],vlb,vub,@(tempsol) nashconendonog(...
         tempsol,alph,bet,gam,n,cfnetrk,cfpio,uyth,...
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
      
      tyopt(i)=temp(2);

      copopt=netrex*kpopt(i)+pioopt;
      netwopt(i)=gam*copopt/(bet*netrex);
      cyopt(i)=netwopt(i)/gam;
      lopt(i)=(cyopt(i)+kpopt(i)-tyopt(i))/netwopt(i);
      coopt(i)=cprod*k(i)^alph*(lopt(i)*n)^(1-alph)-(cyopt(i)+kpopt(i))*n;
      diffuo=log(coopt(i))-uoth;
      taukopt(i)=1-coopt(i)/(alph*cprod*k(i)^alph*(lopt(i)*n)^(1-alph)*(1+diffuo));
      toopt(i)=coopt(i)*diffuo/(1+diffuo);
      R(i)=alph*cprod*k(i)^(alph-1)*(lopt(i)*n).^(1-alph);
      init(:,i)=temp;
   end;
   cfnetrk2=ival*((1-taukopt).*R).^(1/(alph-1));
   cfpio2=ival*toopt;
   cfR=ival*R;
   epsi2=max(max(abs([cfnetrk cfpio]-[cfnetrk2 cfpio2])))
   ite
   if epsi2<epsi,
       taukoptold=taukopt;
       tooptold=toopt;
       cfnetrkold=cfnetrk;
       cfpioold=cfpio;
       cfnetrk=relax*cfnetrk2+(1-relax)*cfnetrk;
       cfpio=relax*cfpio2+(1-relax)*cfpio;
       epsi=epsi2;
       ite=ite+1;
       relax=relaxbase;
       toc;
   else,
       relax=relax/2;
       cfnetrk=relax*cfnetrk2+(1-relax)*cfnetrkold;
       cfpio=relax*cfpio2+(1-relax)*cfpioold;
       if relax<1e-20,
           cfnetrk=cfnetrkold;
           cfpio=cfpioold;
           warning('Failure to converge');
           epsi=0;
       end;
   end;
   if grph==1,
      subplot(2,2,1);
      temp=cos(acos(2*(zx-mink)/(maxk-mink)-1)'*vvec);
      plot(zx,(temp*cfnetrk).^(alph-1),'-',zx,(temp*cfnetrk2).^(alph-1),':');
      title('netrex');
      subplot(2,2,2);
      plot(zx,temp*cfpio,zx,temp*cfpio2,':');
      title('toopt');
      subplot(2,2,3);
      plot(zx,1-netwopt./((1-alph)*cprod*(k./(n*lopt)).^alph));
      subplot(2,2,4);
      %      plot(0,10000*epsi,'+');
      %plot(zx,temp*cfR,'-');
      plot(zx,temp*cfnetrk,zx,temp*cfnetrk2,':');
      pause;
   end;
end;

%Compute variables not required in the iterations
temp=cos(acos(2*(kpopt-mink)/(maxk-mink)-1)*vvec);
netrkopt=(temp*cfnetrk).^(alph-1);
pioopt=temp*cfpio;
taulopt=1-netwopt./((1-alph)*cprod*(k./(n*lopt)).^alph);
copopt=netrkopt.*kpopt+pioopt;
uyopt=log(cyopt)+bet*log(copopt)-gam*lopt;
uoopt=log(coopt);
gnp=cyopt*n/(1+n)+coopt/(1+n)+kpopt*n/(1+n);
gnp2=cprod*k.^alph.*(lopt*n).^(1-alph)/(1+n);
Ropt=(alph*cprod*k.^(alph-1).*(lopt*n).^(1-alph));
wopt=(1-alph)*cprod*k.^(alph).*(lopt*n).^(-alph);
cy2gnp=cyopt*n./(gnp*(1+n));
co2gnp=coopt./(gnp*(1+n));
to2gnp=toopt./(gnp*(1+n));
ty2gnp=tyopt*n./(gnp*(1+n));

subplot(2,2,1);
plot(k,ty2gnp);
set(gca,'FontSize',10);
title(['Transfer to the young as a fraction of GNP']);
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