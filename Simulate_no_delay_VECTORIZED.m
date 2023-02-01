function [pfail, chaingrate] = Simulate_no_delay_VECTORIZED(T,k,rho,n,lambda,simsize)
% Author: Michele Fabi

alpha = n*(1-rho)*lambda;    % honest nodes block rate
beta = n*rho*lambda;         % adversarial block rate


Eblocks = (alpha + beta)*T;
completed = 0;
step = 1;
while completed == 0
   interblockTime = -log(rand(step*ceil(Eblocks),simsize))/(alpha + beta);
   last_block_time = sum(interblockTime,1);
   if sum(last_block_time  < T) > 0
       step = step +1;
   else
       step = 0;
       completed = 1;
       blockTime = cumsum(interblockTime,1);
       interblockTime(blockTime > T) = nan;
       blockTime(blockTime > T) = nan;       
   end
end

Fdraws = rand(size(blockTime)); % random draws
Ablocks = nan(size(blockTime));
Hblocks = nan(size(blockTime));

Ablocks(Fdraws <=rho) = blockTime(Fdraws <=rho);
Hblocks(Fdraws > rho) = blockTime(Fdraws > rho);

%initial block 


 AblocksSorted = sort(Ablocks);
 HblocksSorted = sort(Hblocks);

AblocksSorted(isnan(AblocksSorted)) = Inf;
HblocksSorted(isnan(HblocksSorted)) = Inf;

AblockCount = sum(~isinf(AblocksSorted),1); 
HblockCount = sum(~isinf(HblocksSorted),1); 


A_kth_block = AblocksSorted(k,:);
H_kth_block = HblocksSorted(k,:);


racekresult = H_kth_block - A_kth_block;
runlostFasterk = racekresult > 0; % if honest nodes are slow, the attacker wins

no_pre_mine = find(runlostFasterk == 0);

Diff = nan(size(HblocksSorted(k+1:end,:)));
Diff(:,no_pre_mine) = HblocksSorted(k+1:end,no_pre_mine) - AblocksSorted(k+1:end,no_pre_mine);
runlostFaster = sum(Diff > 0,1) > 0;


lost =  runlostFasterk + runlostFaster;
pfail =  sum(lost)/simsize;

%% COMPUTES HONEST CHAIN GROWTH RATE 
chaingrate = sum(HblockCount/T)/simsize;
end