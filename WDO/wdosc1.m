%-------------------------------------------------------------------------
% Sample Matlab / Octave Code for the Wind Driven Optimization.
% Optimization of the Sphere Function in the range of [-5, 5].
% by Dr. Zikri Bayraktar - thewdoalgorithm@gmail.com
%
% DISCLAIMER: This code is provided for educational purposes
% only. Use at own your risk!
%-------------------------------------------------------------------------
%
% Please refer to the following journal article in your research papers:
% Z. Bayraktar, M. Komurcu, J. A. Bossard and D. H. Werner, "The Wind 
% Driven Optimization Technique and its Application in Electromagnetics," 
% IEEE Transactions on Antennas and Propagation, Volume 61, Issue 5, 
% pages 2745 - 2757, May 2013.
%-------------------------------------------------------------------------
function [outfreq] = wdosc1(num,hi,lamda,gamma,dh,sm)
tic; 
format long g;
% delete('WDOoutput.txt');  
% delete('WDOpressure.txt');  
% delete('WDOposition.txt');
% fid=fopen('WDOoutput.txt','a');
%--------------------------------------------------------------

% User defined WDO parameters:
param.popsize = num;		% population size.
param.npar = 256;			% Dimension of the problem.
param.maxit = 4000;		% Maximum number of iterations.
param.RT = 3;			% RT coefficient.
param.g = 0.2;			% gravitational constant.
param.alp = 0.4;		% constants in the update eq.
param.c = 0.4;			% coriolis effect.
maxV = 0.3;			% maximum allowed speed.
dimMin = zeros(1,param.npar);			% Lower dimension boundary.
dimMax= 255.*ones(1,param.npar);			% Upper dimension boundary.
%---------------------------------------------------------------

% Initialize WDO population, position and velocity:
% Randomize population in the range of [-1, 1]:
pos = 2*(rand(param.popsize,param.npar)-0.5);
% Randomize velocity:
vel = maxV * 2 * (rand(param.popsize,param.npar)-0.5);  
	
%---------------------------------------------------------------

% Evaluate initial population: (Sphere Function)
for K=1:param.popsize,
	x(K,:) = (dimMax - dimMin).*((pos(K,:)+1)./2) + dimMin;
    	pres(K,1) = fitnessfuncsc1(x(K,:),hi,lamda,gamma,dh,sm);
end
%----------------------------------------------------------------

% Finding best air parcel in the initial population :
[globalpres,indx] = min(pres);
globalpos = pos(indx,:);
minpres(1) = min(pres);			% minimum pressure
%-----------------------------------------------------------------

% Rank the air parcels:
[sorted_pres rank_ind] = sort(pres);
% Sort the air parcels:
pos = pos(rank_ind,:);
keepglob(1) = globalpres;
%-----------------------------------------------------------------

% Start iterations :
iter = 1;   % iteration counter
for ij = 2:param.maxit,
    	% Update the velocity:
    	for i=1:param.popsize
		% choose random dimensions:
		a = randperm(param.npar);        			
		% choose velocity based on random dimension:
    		velot(i,:) = vel(i,a);				
        	vel(i,:) = (1-param.alp)*vel(i,:)-(param.g*pos(i,:))+ ...
				    abs(1-1/i)*((globalpos-pos(i,:)).*param.RT)+ ...
				    (param.c*velot(i,:)/i);
    	end
    
        	% Check velocity:
        	vel = min(vel, maxV);
        	vel = max(vel, -maxV);
		% Update air parcel positions:
    		pos = pos + vel;
        	pos = min(pos, 1.0);
        	pos = max(pos, -1.0); 
		% Evaluate population: (Pressure)
		for K=1:param.popsize,
			x(K,:) = (dimMax - dimMin).*((pos(K,:)+1)./2) + dimMin;
    			pres(K,1) = fitnessfuncsc1(x(K,:),hi,lamda,gamma,dh,sm);
		end

    	%----------------------------------------------------
    	% Finding best particle in population
    	[minpres,indx] = min(pres);
    	minpos = pos(indx,:);           	% min location for this iteration
    	%----------------------------------------------------
    	% Rank the air parcels:
    	[sorted_pres rank_ind] = sort(pres);
    	% Sort the air parcels position, velocity and pressure:
    	pos = pos(rank_ind,:);
    	vel = vel(rank_ind,:);
    	pres = sorted_pres;  
    
    	% Updating the global best:
    	better = minpres < globalpres;
    	if better
        		globalpres = minpres;             % initialize global minimum
        		globalpos = minpos;
   	end
	% Keep a record of the progress:
    	keepglob(ij) = globalpres;
    	save WDOposition.txt pos -ascii -tabs;
        
        minval = rank_ind(1);
   ij
end
minval = rank_ind(1);
bestvalue = x(minval,:)
outfreq = bestvalue;
         
% 	%Save values to the final file.
%     	pressure = transpose(keepglob);
%     	save WDOpressure.txt pressure -ascii -tabs;
%         
%         % Plot the pressure function progress over iterations:
%         semilogy(keepglob, 'k' ,'LineWidth',2)
%         title(['Global Best Pressure is " ',num2str(keepglob(1,param.maxit)),' ".'])
%         xlabel('Number of Iterations')
%         ylabel('Global Pressure (i.e. fitness) in log scale')
%         grid on
%         xlim([0, param.maxit])
        
        
    	%END
%-----------------------------------------------------
