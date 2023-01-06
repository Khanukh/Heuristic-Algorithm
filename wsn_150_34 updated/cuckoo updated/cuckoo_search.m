% -----------------------------------------------------------------
% Cuckoo Search (CS) algorithm by Xin-She Yang and Suash Deb      %
% Programmed by Xin-She Yang at Cambridge University              %
% Programming dates: Nov 2008 to June 2009                        %
% Last revised: Dec  2009   (simplified version for demo only)    %
% -----------------------------------------------------------------
% Papers -- Citation Details:
% 1) X.-S. Yang, S. Deb, Cuckoo search via Levy flights,
% in: Proc. of World Congress on Nature & Biologically Inspired
% Computing (NaBIC 2009), December 2009, India,
% IEEE Publications, USA,  pp. 210-214 (2009).
% http://arxiv.org/PS_cache/arxiv/pdf/1003/1003.1594v1.pdf 
% 2) X.-S. Yang, S. Deb, Engineering optimization by cuckoo search,
% Int. J. Mathematical Modelling and Numerical Optimisation, 
% Vol. 1, No. 4, 330-343 (2010). 
% http://arxiv.org/PS_cache/arxiv/pdf/1005/1005.2908v2.pdf
% ----------------------------------------------------------------%
% This demo program only implements a standard version of         %
% Cuckoo Search (CS), as the Levy flights and generation of       %
% new solutions may use slightly different methods.               %
% The pseudo code was given sequentially (select a cuckoo etc),   %
% but the implementation here uses Matlab's vector capability,    %
% which results in neater/better codes and shorter running time.  % 
% This implementation is different and more efficient than the    %
% the demo code provided in the book by 
%    "Yang X. S., Nature-Inspired Metaheuristic Algoirthms,       % 
%     2nd Edition, Luniver Press, (2010).                 "       %
% --------------------------------------------------------------- %

% =============================================================== %
% Notes:                                                          %
% Different implementations may lead to slightly different        %
% behavour and/or results, but there is nothing wrong with it,    %
% as this is the nature of random walks and all metaheuristics.   %
% -----------------------------------------------------------------

function [bestnest,fmin]=cuckoo_search(n)

if nargin<1
% Number of nests (or different solutions)
n=30;
end

% Discovery rate of alien eggs/solutions
pa=0.25;

%% Change this if you want to get better results
% Tolerance
Tol=3300;
maximum_itter =2000;
beta = 10000000;
number_of_test = 30;
%% Simple bounds of the search domain
% Lower bounds
global  excel_string_ouput;
global ssinks;
nd=length(ssinks); 
Lb=1*ones(1,nd); 
% Upper bounds
Ub=4*ones(1,nd);

excel_ouput = zeros(4,number_of_test);

for iii = 1 : number_of_test
tic
% Random initial solutions
for i=1:n,
nest(i,:)=Lb*(randi(4));
end

% Get the current best
fitness=10^10*ones(n,1);
[fmin,bestnest,nest,fitness]=get_best_nest(nest,nest,fitness);

N_iter=0;
%% Starting iterations
while (fmin>Tol && N_iter < maximum_itter),

    % Generate new solutions (but keep the current best)
     new_nest=get_cuckoos(nest,bestnest,Lb,Ub);   
     [fnew,best,nest,fitness]=get_best_nest(nest,new_nest,fitness);
    % Update the counter
      N_iter=N_iter+n; 
    % Discovery and randomization
      new_nest=empty_nests(nest,Lb,Ub,pa) ;
    
    % Evaluate this set of solutions
      [fnew,best,nest,fitness]=get_best_nest(nest,new_nest,fitness);
    % Update the counter again
      N_iter=N_iter+n;
    % Find the best objective so far  
    
    if fnew<fmin,
        fmin=fnew;
        bestnest=round(best);
    end
end %% End of iterations
tend = toc;
excel_ouput(1,iii) = iii;
excel_ouput(2,iii) = tend;
excel_ouput(3,iii) = sum(bestnest>1.5);

excel_ouput(4,iii) = N_iter;
excel_ouput(5,iii) = fmin;


%% output generation 


 

end
format long
numbers_of_failed_divergance  = 0 ;
for i = 1 : number_of_test
    if excel_ouput(5,i) > beta
        numbers_of_failed_divergance = numbers_of_failed_divergance +1 ;
        excel_ouput(5,i) = -1;
    end
end
for i = 1 : number_of_test
    if  excel_ouput(5,i) == -1
                excel_ouput(5,i) = max(excel_ouput(5,:));
                
               max(excel_ouput(5,:))
    end
end
 
excel_ouput(:,number_of_test+1) = mean(excel_ouput,2);
excel_ouput(:,number_of_test+2)=std(excel_ouput,0,2);
excel_ouput(2,number_of_test+3) =numbers_of_failed_divergance;
excel_string_ouput = string(excel_ouput);
excel_string_ouput(1,number_of_test+1) ="mean";
excel_string_ouput(1,number_of_test+2) ="standard deviation";
excel_string_ouput(1,number_of_test+3) ="failed divergance ";
excel_string_ouput(3,number_of_test+3) =" mean time for mean itteration ";
excel_string_ouput(5,number_of_test+3) ="";
 
excel_string_ouput(:,(2 : number_of_test+4)) = excel_string_ouput;
excel_string_ouput(1,1) = "test";
excel_string_ouput(2,1) = "time";
excel_string_ouput(3,1) = "sinks used";
excel_string_ouput(4,1) = "divergance";
excel_string_ouput(5,1) = "energy";
%% Post-optimization processing
%% Display all the nests

disp(strcat('Total number of iterations=',num2str(N_iter)));
fmin
bestnest

%% --------------- All subfunctions are list below ------------------
%% Get cuckoos by ramdom walk
function nest=get_cuckoos(nest,best,Lb,Ub)
% Levy flights
n=size(nest,1);
% Levy exponent and coefficient
% For details, see equation (2.21), Page 16 (chapter 2) of the book
% X. S. Yang, Nature-Inspired Metaheuristic Algorithms, 2nd Edition, Luniver Press, (2010).
beta=3/2;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);

for j=1:n,
    s=nest(j,:);
    % This is a simple way of implementing Levy flights
    % For standard random walks, use step=1;
    %% Levy flights by Mantegna's algorithm
    u=round(randn(size(s))*sigma);
    v=round(randn(size(s)))+0.000001;
    step=u./abs(v).^(1/beta);
  
    % In the next equation, the difference factor (s-best) means that 
    % when the solution is the best solution, it remains unchanged.     
    stepsize=step.*(s-best);
    % Here the factor 0.01 comes from the fact that L/100 should the typical
    % step size of walks/flights where L is the typical lenghtscale; 
    % otherwise, Levy flights may become too aggresive/efficient, 
    % which makes new solutions (even) jump out side of the design domain 
    % (and thus wasting evaluations).
    % Now the actual random walks or flights
    s=round(s+stepsize.*randn(size(s)));
   
   % Apply simple bounds/limits
   nest(j,:)=simplebounds(s,Lb,Ub);
end

%% Find the current best nest
function [fmin,best,nest,fitness]=get_best_nest(nest,newnest,fitness)
% Evaluating all new solutions
for j=1:size(nest,1),
    fnew=fobj(newnest(j,:));
    if fnew<=fitness(j),
       fitness(j)=fnew;
       nest(j,:)=newnest(j,:);
    end
end
% Find the current best
[fmin,K]=min(fitness) ;
best=nest(K,:);

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
for j=1:size(new_nest,1)
    s=new_nest(j,:);
  new_nest(j,:)=simplebounds(s,Lb,Ub);  
end

% Application of simple constraints
function s=simplebounds(s,Lb,Ub)
  % Apply the lower bound
  ns_tmp=s;
  I=ns_tmp<Lb;
  ns_tmp(I)=Lb(I);
  
  % Apply the upper bounds 
  J=ns_tmp>Ub;
  ns_tmp(J)=Ub(J);
  % Update this new move 
  s=ns_tmp;

%% You can replace the following by your own functions
% A d-dimensional objective function
function z=fobj(x)
%% d-dimensional sphere function sum_j=1^d (u_j-1)^2. 
%  with a minimum at (1,1, ...., 1); 
x = round(x);
%% function 
global ssensors;
global ssink_types;
global  ssinks;
global wwsn;
global YY;
global Y_tmp_sensor_total_path_to_sinks;
global Y_Total_communication_energy;
% sensor data   vertex name kc maxl cpu ram bw
sensors = ssensors;
% sink_type data   cpu ram bw cost
sink_types = ssink_types;
 
% sinks   vertex x y 
sinks= ssinks;
sink_vertex = min(sinks)+1;
 
 
 
sink_initial_type = x;
 
 
% wsn
wsn =wwsn;
wsn(:,1:2)= wsn(:,1:2) +1;
 
 
 
nodes = max( max(wsn(:,1)) ,max(wsn(:,2)));
% creating names
str_names = string; 
    for i = 1 : nodes
        if(i <sink_vertex ) 
            str_names(i) = "sensor" + i ;
        end
        if(i >=sink_vertex ) 
            str_names(i) = "sink" + (i - sink_vertex +1 ) ;
        end
    end
G = graph(wsn(:,1),wsn(:,2),wsn(:,3),str_names) ; 
 
 
 
%% Y matrix 
%%% updated
 
% y matrix shortest path for each sink 
Y = YY;
for i = 1 :  sink_vertex -1
    for j = 1 : nodes - sink_vertex + 1  
                    if sink_initial_type(j) == 1
                        Y(i,j) = 0;
                    end         
    end  
end
%% p1 reliability constraint
 
p1 = 0 ; 
for i = 1 : sink_vertex -1
    p1 = p1 +  max(0 , sensors(i,3) - sum(Y(i,:) ) );
end
 
%% p2 max cpu workload  for each sink sum of sensors connected to sink's workload / number of sinks sensor is connected 2
 
p2 = 0 ; 
sink_sum_load = 0;
for i = 1 : nodes - sink_vertex + 1 
    for j =1 : sink_vertex -1
        if (Y(j,i) == 1)
            sink_sum_load = sink_sum_load + sensors(j,5) / sum(Y(j,:));
        end
        
    end
    
    p2 = p2 +  max(0 , sink_sum_load - sink_types(sink_initial_type(i) ,1) ) ;
    p2_data(i,1) = sink_sum_load;
    p2_data(i,2) = sink_types(sink_initial_type(i) ,1);
    sink_sum_load = 0;
end
 
%% p3 max bandwith workload
 
p3 = 0 ; 
sink_sum_load = 0;
for i = 1 : nodes - sink_vertex + 1 
    for j =1 : sink_vertex -1
        if (Y(j,i) == 1 )
            sink_sum_load = sink_sum_load + sensors(j,6) / sum(Y(j,:));
        end
        
    end
    
    p3 = p3 +  max(0 , sink_sum_load - sink_types(sink_initial_type(i) ,2) ) ;
    p3_data(i,1) = sink_sum_load;
    p3_data(i,2) = sink_types(sink_initial_type(i) ,2);
    sink_sum_load = 0;
end
 
%% p4 max memory workload
 
p4 = 0 ; 
sink_sum_load = 0;
for i = 1 : nodes - sink_vertex + 1 
    for j =1 : sink_vertex -1
        if (Y(j,i) == 1 )
            sink_sum_load = sink_sum_load + sensors(j,7) / sum(Y(j,:));
        end
        
    end
    
    p4 = p4 +  max(0 , sink_sum_load - sink_types(sink_initial_type(i) ,3) ) ;
    p4_data(i,1) = sink_sum_load;
    p4_data(i,2) = sink_types(sink_initial_type(i) ,3);
    sink_sum_load = 0;
end
 
%% ps power consumption
 
ps_communication = zeros(1,length(sinks)) ;
 
E_elec = 10; % nj/Bit
K_E_elec = 2;% in Kbit
K_total = length(sensors) * 2;
 
%%%%% clalculating each sensor max path %%% updated
 
YY_tmp_sensor_total_path_to_sinks = Y_tmp_sensor_total_path_to_sinks.*Y;
H =sum(YY_tmp_sensor_total_path_to_sinks,2)';
 
%%%%% clalculating each sensor energy %%% updated
    YY_Total_communication_energy = Y_Total_communication_energy.*Y;
   for i = 1 : length(sensors) 
      for j = 1 : length(sinks)
        if Y(i,j) == 1 
            if (H(i) == 0 )
                YY_Total_communication_energy(i,j) = 0;
                continue;
            end
             YY_Total_communication_energy(i,j) =   (K_E_elec*  E_elec ) * exp( YY_tmp_sensor_total_path_to_sinks(i,j) /H(i));             
        end       
      end
   end
 ps_communication = sum(YY_Total_communication_energy,2)';
%     display ("end sensor "  +  i + " total energy of " + Total_communication_energy )
 
 
%%%%%%%%%%%%%%%%%%%%%%%    process      %%%%%%%%%%%%%%%%%%%%%%%%%%%%
ps_process = zeros(1,length(sinks)) ;
 
sink_sum_load = 0;
for i = 1 : nodes - sink_vertex + 1 
    p_idle = sink_types(sink_initial_type(i) ,6) * sink_types( sink_initial_type(i) ,5 );
    p_residual = sink_types(sink_initial_type(i) ,6) - p_idle ;
    for j =1 : sink_vertex -1
        if (Y(j,i) == 1)
            sink_sum_load = sink_sum_load + sensors(j,5) / sum(Y(j,:));
     
        end
        
    end
        if(sink_types(sink_initial_type(i) ,1) ~= 0)
            
        ps_process(i) = sink_sum_load/sink_types(sink_initial_type(i) ,1)*p_residual + p_idle;
        else
            ps_process(i) = 0  ;
        end
       sink_sum_load = 0;     
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%    end      %%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
PS = sum(ps_process) + sum(ps_communication)/(E_elec + K_E_elec); % dimension lessed of energy
 
 
 
 
%% Fitness Function 
alpha = 1 ;
beta = 10000000 ;
z = alpha*(PS   )  + beta*(p1 + p2 + p3 + p4 )  ;


