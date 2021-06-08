clear; clc;
lambda=1;
m=1;
dn = [1e-9,1e-9];  % differential of the variables A,b

N=2^5;      % 2*N is the number of points taken for each dimension of the grid
L=12;       % Window length (in the range -L:L)
dimtr=2;    % length of the third dimension of space
dimskxy=3;  % length of the third dimension of frequency
[ dx,TH,RHO,skx,sky ] = gridfft2( L,N,dimskxy,dimtr );
ct=0;       % counter

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%----------------------------- SOS algorithm -----------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

npob=25;        % size population
nvar=2;         % number of optiimized parameters
ittot=200;       % maximum number of iterations
lb=[0.65,1];      % lower bound of the parameters
ub=[5*lambda,10];       % upper bound of the parameters

% Variable initialization
graf=10;
fita=zeros(ittot/graf,1);
ind=1;
fit=zeros(npob,1);

% Initial ecosystem
UL=repmat(ub-lb,npob,1);
LB=repmat(lb,npob,1);
xt=rand(npob,nvar).*(UL) + LB;

tic
% Measures the fitness of each organism
for itin=1:npob
    fit(itin)=fitness_sat(lambda,m,dx,TH,RHO,skx,sky,xt(itin,:),dn);
end

% Find the organism with best fitness
ibest=find(fit==min(fit));
ibest=ibest(1);


iter = 0;
while min(fit)>1e-5 && iter<ittot+1
    ct1=0;
for it=1:npob
%------------------------------ Mutualism -------------------------------
% choosing the random organism
jt=randi(npob-1);
jt=jt+(jt>=it);

% Mutual vector
mv=(xt(it,:)+xt(jt,:))/2;

% benefit factor vector
bf=randi(2,[2,1]);

% new organisms
xin=xt(it,:) + rand(1,nvar).*(xt(ibest,:)-mv*bf(1));
xjn=xt(jt,:) + rand(1,nvar).*(xt(ibest,:)-mv*bf(2));
ut=rand(1,nvar).*(ub-lb) + lb;
xin(xin>ub)=ut(xin>ub); xin(xin<lb)=ut(xin<lb);
xjn(xjn>ub)=ut(xjn>ub); xjn(xjn<lb)=ut(xjn<lb);

% objective function value of old organisms 
fiti=fitness_sat(lambda,m,dx,TH,RHO,skx,sky,xt(it,:),dn);
fitj=fitness_sat(lambda,m,dx,TH,RHO,skx,sky,xt(jt,:),dn);

% objective function value of new organisms
fitin=fitness_sat(lambda,m,dx,TH,RHO,skx,sky,xin,dn);
fitjn=fitness_sat(lambda,m,dx,TH,RHO,skx,sky,xjn,dn);

% Replace old organism if the new organisms have a better (lower) fitness
if fitin<fiti
    xt(it,:)=xin;
end
if fitjn<fitj
    xt(jt,:)=xjn;
    ct1=1;
end

%------------------------- Comensalism --------------------------------
% choosing the random organism
jt=randi(npob-1);
jt=jt+(jt>=it);

% new organism
xin=xt(it,:) + (2*rand(1,nvar)-1).*(xt(ibest,:)-xt(jt,:));
ut=rand(1,nvar).*(ub-lb) + lb;
xin(xin>ub)=ut(xin>ub); xin(xin<lb)=ut(xin<lb);

% objective function value of old organisms
fiti=fitness_sat(lambda,m,dx,TH,RHO,skx,sky,xt(it,:),dn);

% objective function value of new organisms
fitin=fitness_sat(lambda,m,dx,TH,RHO,skx,sky,xin,dn);

% Replace old organism if the new organisms have a better (lower) fitness
if fitin<fiti
    xt(it,:)=xin;
    ct1=1;
end

%---------------------------- Parasitism ----------------------------------
% choosing the random organism
jt=randi(npob-1);
jt=jt+(jt>=it);

% number of dimensions that are going to be modified
nd=randi(nvar);
% dimensions that are going to be modified
cd=randperm(nvar);
cd=cd(1:nd);

% parasite vector
xp=xt(it,:);
xp(cd)=rand(1,nd).*(ub(cd)-lb(cd))+lb(cd);

% fitness of host organism
fitj=fitness_sat(lambda,m,dx,TH,RHO,skx,sky,xt(jt,:),dn);

% fitness of parasite organism
fitp=fitness_sat(lambda,m,dx,TH,RHO,skx,sky,xp,dn);

% Replace host organism if parasite organism has a better fitness
if fitp<fitj
    xt(jt,:)=xp;
    ct1=1;
end
end
% Enter his loop if at least one organism changed
if ct1==1
    % Evaluate fitness
    for itin=1:npob
        fit(itin)=fitness_sat(lambda,m,dx,TH,RHO,skx,sky,xt(itin,:),dn);
    end

    % Find best fitness
    ibest=find(fit==min(fit));
    ibest=ibest(1);
end
if min(fit)<1e-4
    ct=ct+1;
    if ct==1
        N=2^5;
        [ dx,TH,RHO,skx,sky ] = gridfft2( L,N,dimskxy,dimtr );
        for itin=1:npob
            fit(itin)=fitness_sat(lambda,m,dx,TH,RHO,skx,sky,xt(itin,:),dn);
        end
        % Find best fitness
        ibest=find(fit==min(fit));
        ibest=ibest(1);
    end
end

if mod(iter,graf)==0
    % Acumulated fitness
    fita(ind)=min(fit);
    ind=ind+1;
end
iter = iter + 1;
end
toc