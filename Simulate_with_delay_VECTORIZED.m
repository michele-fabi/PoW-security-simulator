function [pfail, chaingrate, AvForkRate, AvForkSize, ForkedChain] = Simulate_with_delay_VECTORIZED(T,k,rho,n,lambda,Delta,simsize)
% Author: Michele Fabi

alpha = n*(1-rho)*lambda;    % honest nodes block rate
beta = n*rho*lambda;         % adversarial block rate

%% SIMULATE BLOCK TIMES 
Eblocks = (alpha + beta)*T; %expected number of blocks
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

%% DETERMINE HONEST AND ADVERSARIAL BLOCKS 

Fdraws = rand(size(blockTime)); % random draws
Ablocks = nan(size(blockTime));
Hblocks = nan(size(blockTime));

Ablocks(Fdraws <=rho) = blockTime(Fdraws <=rho);
Hblocks(Fdraws > rho) = blockTime(Fdraws > rho);

%% COMPUTES BLOCKCHAIN GROWTH 

% inizializes blocks adding the first 
Hblocks = [zeros(1,simsize); Hblocks]; 
Ablocks = [zeros(1,simsize); Ablocks]; 

% initialize the variable that sotres summary stats of the blockchain:
ForkedChain = nan(size(Hblocks,1),2,simsize);    
%dim1 = chain length, %dim2 = fork time, %dim3 = simulation index
% NB: heightindex = 1 corresponds to the initial block, already existing at t=0 (with height 0).

% initialization 
    heightindex = ones(simsize,1); % vector with creation times of reference honest blocks
    Hheightindex = ones(simsize,1); % vector with creation times reference honest blocks on the longest chain 
    chainlength = size(Hblocks,1);
    completedindex = zeros(simsize,1); %index that keeps track of which simulations have been fully analyzed
    incompleindexpos = 1:simsize; % at start no simulation has been scanned
    simsizeids = (1:simsize)'; % simulation id's 
    
% to deal with all simulations simultaneously we loop over a vector index
% completedindex indicates which simulations have been fully scanned by the
% loop
while sum(completedindex) < simsize

    % finds linear index of next blocks on the longest chain in all
    % simulations have not been scanned fully
 linearheightindex = sub2ind(size(Hblocks),heightindex(incompleindexpos),simsizeids(incompleindexpos)); 
 refblocks = Hblocks(linearheightindex); % gets the creation time of blocks on the longest chain 

 refsimblocks = -1*ones(simsize,1);
 refsimblocks(incompleindexpos) = refblocks;

 % get nan and non-nan values
 refHblocks = refsimblocks(~isnan(refsimblocks) & refsimblocks ~= -1);


% finds indices: 
 Hblockssimi = find(~isnan(refsimblocks) & refsimblocks ~= -1); % honest blocks in simulations that are not fully scanned yet
 Ablockssimi = find(isnan(refsimblocks)); % adversarial blocks
 
 % this pice of codes computes forks  
 refHblocksvec = repmat(refHblocks',chainlength,1);     
 % the number Nchains of chains that extend a reference block in the longest chain 
 % equals number of blocks that are less than Delta away from it
 Nchains = sum(Hblocks(:,Hblockssimi) - refHblocksvec <= Delta ...
            & Hblocks(:,Hblockssimi) - refHblocksvec >= 0,1); 
 % get positions in ForkedChain where to store block time of new blocks on
 % the logest chain
 timeids = sub2ind(size(ForkedChain),Hheightindex(Hblockssimi),1*ones(length(Hblockssimi),1),Hblockssimi); 
 %%get positions in ForkedChain where to store the number of forks
 %%extending new blocks on the longest chain 
 Nchainids = sub2ind(size(ForkedChain),Hheightindex(Hblockssimi),2*ones(length(Hblockssimi),1),Hblockssimi);
 
 % store creation times of new blocks on the longest chain and number of
 % forks 
 ForkedChain(timeids) = refHblocks;
 ForkedChain(Nchainids) = Nchains;

 % update indexes
 heightindex(Ablockssimi) = heightindex(Ablockssimi) + 1;           % moves to next block for the adversary
 heightindex(Hblockssimi) = heightindex(Hblockssimi) + Nchains';    % moves after forked blocks for the honest miners
 
 % non-nan ref blocks
 Hheightindex(Hblockssimi) = Hheightindex(Hblockssimi) + ones(length((Hblockssimi)),1); %updates lenth of the logest chain

%%% check which simulations have been fully scanned
completedindex = heightindex >= (chainlength+1)*ones(size(Hblocks,2),1); %simulation is scanned if all blocks have been analyzed 
incompleindexpos =  find(~completedindex);
end
%% COMPUTE FAIL PROBABILITY

HblocksNet = squeeze(ForkedChain(:,1,:)); % obtain times of block on the longest chain 

AblocksSorted = [sort(Ablocks(2:end,:))];
HblocksSorted = [sort(HblocksNet(2:end,:))];

% replace blocks that never arrive with infinite block time 
AblocksSorted(isnan(AblocksSorted)) = Inf;
HblocksSorted(isnan(HblocksSorted)) = Inf;

% count # blocks 
AblockCount = sum(~isinf(AblocksSorted),1); 
HblockCount = sum(~isinf(HblocksSorted),1); 

% block time of the k-th block 
A_kth_block = AblocksSorted(k,:);
H_kth_block = HblocksSorted(k,:);

racekresult = H_kth_block - A_kth_block;
runlostFasterk = racekresult > 0; % attacker catches up releasing a longer chain once the 
% confirmation depth is reached

no_pre_mine = find(runlostFasterk == 0); % in these simulation the attacker won by releasing a longer chain 
% at the time the confirmation lag is reached 

Diff = nan(size(HblocksSorted(k+1:end,:)));
Diff(:,no_pre_mine) = HblocksSorted(k+1:end,no_pre_mine) - AblocksSorted(k+1:end,no_pre_mine);
runlostFaster = sum(Diff > 0,1) > 0; % in these simulation the attacker won by racing with honest miners after the confirmation depth
% has been reached 

lost =  runlostFasterk + runlostFaster; % the attacker wins in both of the previous events
pfail =  sum(lost)/simsize; % probability that the protocol is attacked

%% COMPUTES FORK RATE 
Chains = squeeze(ForkedChain(:,2,:));  %obtian number of branches at each height of the longest chain 
NForks =  sum(Chains > 1,1); % the longest chain is forked with more than one branch 
Forkrate = NForks/T;
AvForkRate = sum(Forkrate)/simsize;

AvForkSize =  sum(Chains(Chains > 1))/sum(NForks);

%% COMPUTES GROWTH RATE OF THE LONGEST CHAIN
chaingrate = sum(HblockCount/T)/simsize;
end