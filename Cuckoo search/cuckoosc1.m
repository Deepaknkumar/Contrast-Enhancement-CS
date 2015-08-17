%sir,in this code,changing N_iterTotal at line 30, we will get respective
%graphs,
%for PSNR Plot to gradually increase,uncommment the lines from 20 to 23
%and comment the lines, 16 to 18

function [bestnest,fmin,psnrc,N_iterTotal] = cuckoosc1(n,hi,lamda,gamma,dh,sm)

if nargin<1
    n = 25;
end

number_of_solution = 256;
Lb = zeros(1,number_of_solution);
Ub = 255.*ones(1,number_of_solution);

for i=1:n
    nest(i,:) = Lb + (Ub - Lb).*rand(size(Lb));
end

% for i=1:n-1
%     nest(i,:) = Lb + (Ub - Lb).*rand(size(Lb));
% end
% nest(i,:) = hi;

pa = .25;           %discovery rate

fitness = 10^10.*ones(n,1);
[fmin,bestnest,nest,fitness] = get_best_nest(nest,nest,fitness,hi,lamda,gamma,dh,sm);

N_iterTotal = 5000; %total number of iterations

N_iter = 0;

for iter=1:N_iterTotal
    iter
    new_nest = get_cuckoos(nest,bestnest,Lb,Ub);
    [fnew,best,nest,fitness] = get_best_nest(nest,new_nest,fitness,hi,lamda,gamma,dh,sm);
    
    N_iter = N_iter + n;
    
    new_nest = empty_nests(nest,Lb,Ub,pa);      %discovery and Randomization
    [fnew,best,nest,fitness] = get_best_nest(nest,new_nest,fitness,hi,lamda,gamma,dh,sm);
    
    N_iter = N_iter + n;
    
    if fnew < fmin
        fmin = fnew;
        bestnest = best;
    end
fminval(iter) = fmin;
psnrc(iter,:) = bestnest;
end

    x = 1:1:N_iterTotal;
    figure,plot(x,fminval,'LineWidth',2);
    title('Convergence plot ( fmin vs no. of iterations)');
    xlabel('Number Of Iterations');
    ylabel('Value of objective function (fmin)');
    hold on
fmin
bestnest

%subfunctions

%cuckoos using levy flight
function nest = get_cuckoos(nest,best,Lb,Ub)
n = size(nest,1);
beta=3/2;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);

for j=1:n
    s=nest(j,:);
    u =randn(size(s))*sigma;
    v =randn(size(s));
    step = u./abs(v).^(1/beta);
    
    stepsize = .01*step.*(s - best);
    
    s = s + stepsize.*randn(size(s));

    nest(j,:) = SimpleBounds(s,Lb,Ub);
end

%current best nest

function [fmin,best,nest,fitness] = get_best_nest(nest,newnest,fitness,hi,lamda,gamma,dh,sm)

for j = 1:size(nest,1)
    fnew = fobj(newnest(j,:),hi,lamda,gamma,dh,sm);
    if fnew < fitness(j)
        fitness(j) = fnew;
        nest(j,:) = newnest(j,:);
    end
end
%current bests
[fmin,K] = min(fitness);
best = nest(K,:);

%replacing some nests by constructing new solutions/nests
function new_nest = empty_nests(nest,Lb,Ub,pa)
n = size(nest,1);

k = rand(size(nest)) > pa;
stepsize = rand*(nest(randperm(n),:) - nest(randperm(n),:));
 new_nest = nest + stepsize.* k ;
 
 for j=1:size(new_nest)
     s = new_nest(j,:);
     new_nest(j,:) = SimpleBounds(s,Lb,Ub);
 end
 
 %function for simple bounds
    function s  = SimpleBounds(s,Lb,Ub)
        ns_temp = s;
        I = s<Lb;
        ns_temp(I) = Lb(I);
        
        J = s>Ub;
        ns_temp(J) = Ub(J);
        
        s = ns_temp;
        
 %objective Function
 
        function z = fobj(u,hi,lamda,gamma,dh,sm)
            dh = dh*transpose(u);
            dh = sum(dh);
            s1 = sum((transpose(u)-hi).^2);
            s2 = sum((transpose(u)-sm).^2);
            s3 = sum(dh.^2);
            z = s1 + lamda*s2 + gamma*s3;