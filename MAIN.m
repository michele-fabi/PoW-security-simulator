%% CODE FOR ``Blockchain Reskilling''
% Nakamoto Consensus
% Author: Michele Fabi
% PRELIMINARY VERSION!
%% INITIALIZATION
close all
clear 
clc
rng('default'); % sets simulation seed

% find path of current .m file and sets scripts folder as main. 
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile);
addpath(fullfile(pathstr,'/Scripts'));
userpath(pathstr);
%% CONFIGURATION

n = 100;                     % # participating miners
k = 6;                       % confirmation depth;
rho = 1/3;                   % fraction of malicious miners (forming ''the attacker'' as a colaition) 
T = 60*24;                   % total attack time window in minutes (e.g. 60*24*7 = one week)
lambda = 1/(10*n);           % block rate:  lambda = 6/n ->  one block each 10 seconds
                             %              lambda = 1/(10*n) -> one block
                             %              each 10 minutes 
Delta = 1/6;                   % transmission delay in minutes (Delta = 1/6 -> 10 seconds delay)

%%%%% used directly within the simulation:
%%%   alpha = n*(1-rho)*lambda;    % honest nodes block rate
%%%   beta = n*rho*lambda;         % adversarial block rate
%%%%%

%% SET SIMULATION SIZE 

%%% Increases the simulation length until simulated mining rates convergence
%%% to the theoretical ones 
%%% tic/toc measure time to convergence

eps = 10^(-4); % convergence tollerance 
stepsim = 100; % simulation size increase in each iteration

format long % display more decimals

 tic 
 [simSize, interBlockTime] = Compute_Simulation_Size(T,n,lambda,eps,stepsim);
 toc
 %histogram(interBlockTime) % just a sanity check: they have to be exponentially distributed

format short % restore short number format
 

%% SINGLE SIMULATION ROUND 

simsize = 10^4;  % uncomment to set simulation size manually 

rng(1) % reset the seed at same value to have the same random draws. 
tic
[pfailD, chaingrateD, AvForkRate, AvForkSize, ForkedChain] = Simulate_with_delay_VECTORIZED(T,k,rho,n,lambda,Delta,simsize);
toc

rng(1)
tic
[pfail, chaingrate]= Simulate_no_delay_VECTORIZED(T,k,rho,n,lambda,simsize);
toc
%% GRAPHS 
stepsim = 100; % simulation size increase in each iteration
eps = 10^(-5); % convergence tollerance 

%%% show that attack probability falls exponentially in confirmation lag
kvals = 1:25;
pfailvec = nan(size(kvals));
for k=kvals
    [simSize, ~] = Compute_Simulation_Size(T,n,lambda,eps,stepsim);
    [pfailD, ~, ~, ~, ~] = Simulate_with_delay_VECTORIZED(T,k,rho,n,lambda,Delta,simsize);
pfailvec(kvals == k) = pfailD;
k
end    
figure()
plot(pfailvec);
ylabel('attack probability');
xlabel('confirmation depth');


%%% show that fork rate is a convex fn. of the block rate and fork size grows
%%% linearly with block rate
k = 6;
T= 60; % one hour
lamvec = 1/n.*[1/10:1/10:2];
forkvec = nan(length(lamvec),2);
for lambda=lamvec
    [simSize, ~] = Compute_Simulation_Size(T,n,lambda,eps,stepsim);
    [~, ~, AvForkRate, AvForkSize, ~] = Simulate_with_delay_VECTORIZED(T,k,rho,n,lambda,Delta,simsize);
forkvec(lamvec == lambda,1) = AvForkRate;
forkvec(lamvec == lambda,2) = AvForkSize;
lambda
end    

figure()
plot(lamvec,forkvec(:,1));
ylabel('fork rate');
xlabel('block rate');

figure()
plot(lamvec,forkvec(:,2));
ylabel('av. fork size');
xlabel('block rate');


%%% show that growth rate of logest chain reduces as tramission delay
%%% increases
k = 6;
lambda = 1/(10*n);
Deltavec = 1/2:1/2:20;
gratevec = nan(length(Deltavec),1);
for Delta=Deltavec
    [simSize, ~] = Compute_Simulation_Size(T,n,lambda,eps,stepsim);
    [~, grate, ~, ~, ~] = Simulate_with_delay_VECTORIZED(T,k,rho,n,lambda,Delta,simsize);
gratevec(Deltavec == Delta,1) = grate;
Delta
end    


figure()
plot(Deltavec,gratevec);
ylabel('LC growth rate');
xlabel('transmission delay');

