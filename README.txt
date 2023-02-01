CODE FOR ''BLOCKCHAIN RESKILLING'' 
% AUTHOR: Michele Fabi
% Matlab version used: R2021b



This archive contains the Matlab code to analyze a mining race between honest and 
mailicious miners. The script MAIN.m will initialize the code and plot security metrics 
as a function of the protocol's parameters. 

Apart from MAIN.m, this archive contains the following functions: 

FUNCTION								|PURPOSE

Simulate_no_delay_VECTORIZED					Outputs the probability that the attack succeeds and the growth rate of the honest chain.	
Simulate_with_delay_VECTORIZED				Outputs the probability that the attack succeeds, growth rate of the honest chain, and several statistics on the effects of forking 
Compute_Simulation_Size						Compute the simulation size as a fn of a precision threshold
rowwiseLast								Finds the locations of the final non-zero value in each row of a matrix**

All codes are preliminary and in testing phase. 
For feedback, please contact me at my institutional email (michele.fabi@ensae.fr). 


**I thank the online Mathworks community, in particular user 'Matt J', for making rowwiseLast publicly available.
