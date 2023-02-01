function [simsize, interblockTime] = Compute_Simulation_Size(T,n,lambda,eps,stepsim)
% Author: Michele Fabi


Lambda = lambda*n; % aggregate mining rate
Eblocks = Lambda*T; % expected number of blocks mined in the reference time window

% initialization 
simrate = 1000*Lambda; % pick arbitrarily large simulated rate to start (factor *1000 is arbitrary)
simsize = 0; 

while abs(simrate - Lambda) > eps % convergence criterion 

simsize = simsize + stepsim;
completed = 0; 
stepblocks = 1;
while completed == 0
    % simulates inter-arrival block times 
    interblockTime = -log(rand(stepblocks*ceil(Eblocks),simsize))/Lambda; 
    % the previous step uses the inverse-CDF technique: 
    % quantiles are uniformely distributed over [0,1].
    % since waiting times are exponentially distributed with CDF 
    % F = 1 - exp(- Lambda*t),
    % we can invert F directly and solve for t: t = - ln(1-F)/Lambda
    % since F is uniform[0,1], 1-F is also uniform[0,1] 

    last_block_time = sum(interblockTime,1); %last block mined within the time horizon 

   if sum(last_block_time  < T) > 0
       stepblocks = stepblocks +1; % need more time 
   else 
       stepblocks = 0;
       completed = 1;
       blockTime = cumsum(interblockTime,1);
       interblockTime(blockTime > T) = nan;
       blockTime(blockTime > T) = nan;       
   end
end

blockRateVec = sum(~isnan(blockTime),1)/T;
simrate = sum(blockRateVec)/simsize;
[simrate Lambda]
end
