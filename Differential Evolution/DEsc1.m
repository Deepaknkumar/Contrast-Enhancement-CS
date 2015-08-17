%% MODE
% Multi-objective Evolutionary Algorithm (MOEA) based on Differential
% Evolution (DE).
% It implements a greedy selection based on pure dominance.
% DE algorithm has been introduced in:
%
% Storn, R., Price, K., 1997. Differential evolution: A simple and 
% efficient heuristic for global optimization over continuous spaces. 
% Journal of Global Optimization 11, 341 – 359.
%%
%% Beta version 
% Copyright 2006 - 2012 - CPOH  
% Predictive Control and Heuristic Optimization Research Group
%      http://cpoh.upv.es
% ai2 Institute
%      http://www.ai2.upv.es
% Universitat Politècnica de València - Spain.
%      http://www.upv.es
%%
%% Author
% Gilberto Reynoso Meza
% gilreyme@upv.es
% http://cpoh.upv.es/en/gilberto-reynoso-meza.html
% http://www.mathworks.es/matlabcentral/fileexchange/authors/289050
%%
%% For new releases and bug fixing of this Tool Set please visit:
% http://cpoh.upv.es/en/research/software.html
% Matlab Central File Exchange
%%
%% Overall Description
% This code implements a basic multi-objective optimization algorithm based
% on Diferential Evolution (DE) algorithm.
%
% When one objective is optimized, the standard DE runs; if two or more
% objectives are optimized, the greedy selection step in DE algorithm is 
% performed using a dominance relation.
%%
%
%% 

function [outfreq]=DEsc1(num,hi,lamda,gamma,dh,sm)

NOBJ = 1;                          % Number of objectives
MODEDat.NRES = 0;                          % Number of constraints
Nvar   = 256;                       % Numero of decision variables
MODEDat.mop = str2func('fitnessfunc');    % Cost function

MODEDat.FieldD = [zeros(256,1),255*ones(256,1)]; % Initialization bounds
MODEDat.Initial= [zeros(256,1),255*ones(256,1)]; % Optimization bounds

%MODEDat.XPOP = n*NOBJ;             % Population size
% MODEDat.Esc = 0.5;                         % Scaling factor
% MODEDat.Pm= 0.2;                           % Croosover Probability

%MODEDat.MAXGEN =10000;                     % Generation bound
MODEDat.MAXFUNEVALS = 150*Nvar...  % Function evaluations bound
    *NOBJ;   

MODEDat.CounterGEN=0;
MODEDat.CounterFES=0;
%% Reading parameters from MODEDat
Generaciones  = 4000; %MODEDat.MAXGEN;    % Maximum number of generations.
Xpop          =   num*NOBJ; %MODEDat.XPOP;      % Population size.
%Nvar;          %= Nvar;      % Number of decision variables.
Nobj          = 1; % MODEDat.NOBJ;      % Number of objectives.
Bounds        = MODEDat.FieldD;    % Optimization bounds.
Initial       = MODEDat.Initial;   % Initialization bounds.
ScalingFactor = 0.5;  %MODEDat.Esc;       % Scaling fator in DE algorithm.
CrossOverP    = 0.2;  %MODEDat.Pm;        % Crossover probability in DE algorithm.
mop           = MODEDat.mop;       % Cost function.
MODEDat.InitialPop = [];


%% Initial random population
Parent = zeros(Xpop,Nvar);  % Parent population.
Mutant = zeros(Xpop,Nvar);  % Mutant population.
Child  = zeros(Xpop,Nvar);  % Child population.
FES    = 0;                 % Function Evaluation.

for xpop=1:Xpop
    for nvar=1:Nvar
        Parent(xpop,nvar) = Initial(nvar,1)+(Initial(nvar,2)...
                            - Initial(nvar,1))*rand();
    end
end

if size(MODEDat.InitialPop,1)>=1
    Parent(1:size(MODEDat.InitialPop,1),:)=MODEDat.InitialPop;
end

sof = size(Parent);
for i=1:sof(1)
    JxParent(i,1) = fitnessfuncsc1(Parent(i,:),hi,lamda,gamma,dh,sm);
end
FES = FES+Xpop;   

%% Evolution process

for n=1:Generaciones 
    
    for xpop=1:Xpop
        rev=randperm(Xpop);
        
        %% Mutant vector calculation
        Mutant(xpop,:)= Parent(rev(1,1),:)+ScalingFactor*...
                       (Parent(rev(1,2),:)-Parent(rev(1,3),:));
        
        for nvar=1:Nvar %Bounds are always verified
            if Mutant(xpop,nvar)<Bounds(nvar,1)
                Mutant(xpop,nvar) = Bounds(nvar,1);
            elseif Mutant(xpop,nvar)>Bounds(nvar,2)
                Mutant(xpop,nvar)=Bounds(nvar,1);
            end
        end
        
        %% Crossover operator
        for nvar=1:Nvar
            if rand() > CrossOverP
                Child(xpop,nvar) = Parent(xpop,nvar);
            else
                Child(xpop,nvar) = Mutant(xpop,nvar);
            end
        end

    end

    sofc = size(Child);
for i=1:sofc(1)
    JxChild(i,1) = fitnessfuncsc1(Child(i,:),hi,lamda,gamma,dh,sm);
end

    FES=FES+Xpop;

    %% Selection
    for xpop=1:Xpop
        if JxChild(xpop,:) <= JxParent(xpop,:) 
            Parent(xpop,:) = Child(xpop,:);
            JxParent(xpop,:) = JxChild(xpop,:);
        end
    end
    
	PFront=JxParent;
	PSet=Parent;

    OUT.Xpop           = Parent;   % Population
    OUT.Jpop           = JxParent; % Poopulation's Objective Vector
    OUT.PSet           = PSet;     % Pareto Set
    OUT.PFront         = PFront;   % Pareto Front
    OUT.Param          = MODEDat;  % MODE Parameters
    MODEDat.CounterGEN = n;
    MODEDat.CounterFES = FES;
    
%     [OUT MODEDat]=PrinterDisplay(OUT,MODEDat); % To print results on screen
%     
%     if FES>MODEDat.MAXFUNEVALS || n>MODEDat.MAXGEN
%         disp('Termination criteria reached.')
%         break;
[OUT.PFront, OUT.PSet]=DominanceFilter(PFront,PSet);
n
end


OUT.Xpop=PSet;
OUT.Jpop=PFront;
[OUT.PFront, OUT.PSet]=DominanceFilter(PFront,PSet); %A Dominance Filter
outfreq = OUT.PSet;

%% Dominance Filter
%
% A filter based on dominance criteria
%
function [Frente Conjunto]=DominanceFilter(F,C)

Xpop=size(F,1);
Nobj=size(F,2);
Nvar=size(C,2);
Frente=zeros(Xpop,Nobj);
Conjunto=zeros(Xpop,Nvar);
k=0;

for xpop=1:Xpop
    Dominado=0;
    
    for compara=1:Xpop
        if F(xpop,:)==F(compara,:)
            if xpop > compara
                Dominado=1;
                break;
            end
        else
            if F(xpop,:)>=F(compara,:)
                Dominado=1;
                break;
            end
        end
    end
    
    if Dominado==0
        k=k+1;
        Frente(k,:)=F(xpop,:);
        Conjunto(k,:)=C(xpop,:);
    end
end
Frente=Frente(1:k,:);
Conjunto=Conjunto(1:k,:);
