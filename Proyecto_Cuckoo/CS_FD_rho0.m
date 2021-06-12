function [bestnest]=CS_FD_rho0(a,lambda,DA)
if nargin<1
   lambda=1;a=15;DA=1; 
elseif nargin<2
   lambda=1;DA=1;
elseif nargin<3
   DA=1;
end
% a is half of the window width
% Lambda is a parameter of the differential equation
%% Cuckoo search parameter
pa=0.5;% Discovery rate of alien eggs/solutions
%% Change this if you want to get better results
N_IterTotal=8000;
n=40; %Number of solutions proposed
du=1e-10;%Step to calculate the derivative
%% Simple bounds of the search domain
nd=6;%Dimensionality of the problem
Lb=0.85.*ones(1,nd);%Lower bounds
Ub=10.*ones(1,nd);%.*ones(1,nd);%Upper bounds
% Random initial solutions
Lb=repmat(Lb,n,1);Ub=repmat(Ub,n,1);
nest=Lb+(Ub-Lb).*rand(size(Lb));
%% Function Parameters
    N=128;%Number of points
    x=linspace(-a,a,N+1);
    x=x(1:end-1);
    [X,Y]=meshgrid(x);
%%    
    rho2(:,:,1)=(X-10).^2+(Y).^2;
    rho2(:,:,2)=(X+10).^2+(Y).^2;
    rho2(:,:,3)=(X).^2+(Y-10).^2;
    %rho2(:,:,4)=(X).^2+(Y+5).^2;
    theta(:,:,1)=atan2(Y,X-10);
    theta(:,:,2)=atan2(Y,X+10);
    theta(:,:,3)=atan2(Y-10,X);
    %theta(:,:,4)=atan2(Y+5,X);
    rho=sqrt(rho2);
   
%% Fourier Parameters
    dx=x(2)-x(1);
    dk=(2*pi/(2*a));
    kx=X.*1/dx*dk;
    ky=Y.*1/dx*dk;
    kx=fftshift(kx);
    ky=fftshift(ky);
%% Get the current best
fitness=10^10*ones(n,1); %Initial Fitness Values
[fmin,bestnest,nest,fitness]=get_best_nest(nest,nest,fitness,du,kx,ky,rho2,rho,theta,lambda,dx,N,DA);
N_iter=0;
iter=1;
%% Starting iterations
while (iter<N_IterTotal) && fmin>0.1
    % Generate new solutions (but keep the current best)
     new_nest=get_cuckoos(nest,bestnest,Lb,Ub);   
     [fnew,best,nest,fitness]=get_best_nest(nest,new_nest,fitness,du,kx,ky,rho2,rho,theta,lambda,dx,N,DA);
    % Update the counter
      N_iter=N_iter+n; 
    % Discovery and randomization
      new_nest=empty_nests(nest,Lb,Ub,pa) ;
    
    % Evaluate this set of solutions
      [fnew,best,nest,fitness]=get_best_nest(nest,new_nest,fitness,du,kx,ky,rho2,rho,theta,lambda,dx,N,DA);
    % Update the counter again
      N_iter=N_iter+n;
    % Find the best objective so far  
        if fnew<fmin
            fmin=fnew;
            bestnest=best;
        end
        if fmin<du
           du=fmin/2;
        end
    iter=iter+1;
%     if mod(iter,10)==0
%     disp([du,fmin,bestnest,iter])
%     end
end %% End of iterations
%% Post-optimization processing
%% Display all the nests
% bestnest=[bestnest,fmin];
disp(fmin)
end
%% --------------- All subfunctions
%% Get cuckoos by ramdom walk
function nest=get_cuckoos(nest,best,Lb,Ub)
% Levy flights
%n=size(nest,1);
% Levy exponent and coefficient
% For details, see equation (2.21), Page 16 (chapter 2) of the book
% X. S. Yang, Nature-Inspired Metaheuristic Algorithms, 2nd Edition, Luniver Press, (2010).
beta=3/2;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    s=nest;
    % This is a simple way of implementing Levy flights
    % For standard random walks, use step=1;
    %% Levy flights by Mantegna's algorithm
    u=randn(size(s)).*sigma;
    v=randn(size(s));
    step=u./abs(v).^(1/beta);
  
    % In the next equation, the difference factor (s-best) means that 
    % when the solution is the best solution, it remains unchanged.     
    stepsize=0.01*step.*(s-best);
    % Here the factor 0.01 comes from the fact that L/100 should the typical
    % step size of walks/flights where L is the typical lenghtscale; 
    % otherwise, Levy flights may become too aggresive/efficient, 
    % which makes new solutions (even) jump out side of the design domain 
    % (and thus wasting evaluations).
    % Now the actual random walks or flights
    s=s+stepsize.*randn(size(s));
   % Apply simple bounds/limits
   nest=simplebounds(s,Lb,Ub);
end
%% Find the current best nest
function [fmin,best,nest,fitness]=get_best_nest(nest,newnest,fitness,du,kx,ky,rho2,rho,theta,lambda,dx,N,DA)
% Evaluating all new solutions
    newnest4=permute(newnest,[4,2,3,1]);
    fnew=fobj(newnest4,du,kx,ky,rho2,rho,theta,lambda,dx,N,DA);
    fnew=permute(fnew,[2,1]);
    I=fnew<fitness;
    fitness(I)=fnew(I);
    nest(I,:)=newnest(I,:);
[fmin,K]=min(fitness) ;
best=nest(K,:);
end
%% Replace some nests by constructing new solutions/nests

function new_nest=empty_nests(nest,Lb,Ub,pa)
% A fraction of worse nests are discovered with a probability pa
n=size(nest,1);
% Discovered or not -- a status vector
K=rand(size(nest))>pa;
% In the real world, if a cuckoo's egg is very similar to a host's eggs, then 
% this cuckoo's egg is less likely to be discovered, thus the fitness should 
% be related to the difference in solutions.  Therefore, it is a good idea 
% to do a random walk in a biased way with some random step sizes.  
%% New solution by biased/selective random walks
stepsize=rand*(nest(randperm(n),:)-nest(randperm(n),:));
new_nest=nest+stepsize.*K;
  s=new_nest;
  new_nest=simplebounds(s,Lb,Ub);
end
% Application of simple constraints
function s=simplebounds(s,Lb,Ub)
  % Apply the lower bound
  ns_tmp=s;
  I=ns_tmp<Lb;
  ns_tmp(I)=Lb(I)+(Ub(I)-Lb(I)).*rand(size(Lb(I)));
  
  % Apply the upper bounds 
  J=ns_tmp>Ub;
  ns_tmp(J)=Ub(J)+(Ub(J)-Ub(J)).*rand(size(Ub(J)));
  % Update this new move 
  s=ns_tmp;
end
%% Fitness Funtion
function z=fobj(u,du,kx,ky,rho2,rho,theta,lambda,dx,n,DA)
%% Creating the array u+du
N=size(u,2);
v=eye(N,N);
v=permute(v,[3,2,1]);
if DA==2
    v=cat(3,2.*v,v,-v,-2.*v);
else
    v=cat(3,v,-v);
end
u=u+v.*du;

%% Proposed Ansatz
 U1=u(1,1,:,:).*exp(-rho2(:,:,1)./u(1,2,:,:).^2).*rho2(:,:,1).*(cosd(20).*cos(2.*theta(:,:,1))+1i.*sind(20).*sin(2.*theta(:,:,1)));
 U2=u(1,3,:,:).*exp(-rho2(:,:,2)./u(1,4,:,:).^2).*rho2(:,:,2).*(cosd(30).*cos(2.*theta(:,:,2))+1i.*sind(30).*sin(2.*theta(:,:,2)));
 U3=u(1,5,:,:).*exp(-rho2(:,:,3)./u(1,6,:,:).^2).*rho2(:,:,3).*(cosd(40).*cos(2.*theta(:,:,3))+1i.*sind(40).*sin(2.*theta(:,:,3)));
% A0=0.9195/1; B0=2.5608;
% A1=0.9401/1; B1=2.5947;
% A2=0.9651/1; B2=2.618;

%S=0.05
U=U1+U2+U3;
%U=U_1+U_2+U_3;%+U_4;
%2^(-(m+1)/2).*
%% Derivatives and other calculations
%% with respect to X and Y
FU=fft2(U);
FUx=ifft2(1i.*kx.*FU); 
FUy=ifft2(1i.*ky.*FU);
FUx=abs(FUx).^2;FUy=abs(FUy).^2;
U2=abs(U).^2;
%U4=U2.^2;

%% Medium
%NL=-1/2.*U4; %Kerr
s=0.05;
NL=(log(s.*U2+1)-s.*U2)./s^2; %Saturation

%% Differential Equationsize()
G=(lambda.*U2+FUx+FUy+NL);

%% Simpsons 1/3 integration Rule
D=(2*dx./3).*(G(end,1,:,:)+G(1,1,:,:)+2.*sum(G(2:2:(n-1),:,:,:))+4.*sum(G(3:2:n,:,:,:)));
D=(2*dx./3).*(D(1,end,:,:)+D(1,1,:,:)+2.*sum(D(1,2:2:(n-1),:,:))+4.*sum(D(1,3:2:n,:,:)));
D=permute(D,[3,2,1,4]);

%% Finite difference Derivative for A, B, C ... Accuracy 4
if DA==2
    D=(-D(1:N,:)+8.*D((N+1):2*N,:)-8.*D((2*N+1):3*N,:)+D((N+1):2*N,:))./(12*du);
else
    D=(D(1:N,:)-D((N+1):end,:))./(2*du);
end
%   4.3851    1.4519    4.1880    2.0114    0.7800    0.7800
% 
%    4.2476    1.4541    2.4690    2.5310    1.0137    2.8330
%   2.5409    1.8755    2.4743    1.9465    3.0709    1.8156
%% Gradient
z=sqrt(sum(D.*D));
end