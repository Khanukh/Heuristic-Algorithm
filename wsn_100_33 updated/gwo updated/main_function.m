clear all 
clc



SearchAgents_no=30; % Number of search agents
Function_name = 'F10';
  lb=-100;
  ub=100;
  step = 1;
  dim= 4; % sink numbers
  Max_iteration=100;
%% fitness function
fobj = @F10;


%% algorithm

[Best_score,Best_pos,GWO_cg_curve]=GWO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);

%% plot 
figure('Position',[500 500 660 290])
%Draw search space
subplot(1,2,1);
x=lb:step:ub; y=x;%[-500,500]
L=length(x);
f=[];
for i=1:L
    for j=1:L
          f(i,j)=fobj([x(i),y(j)]);
    end
end
surfc(x,y,f,'LineStyle','none');
title('Parameter space')
xlabel('x_1');
ylabel('x_2');
zlabel([Function_name,'( x_1 , x_2 )'])

%Draw objective space
subplot(1,2,2);
semilogy(GWO_cg_curve,'Color','r')
title('Objective space')
xlabel('Iteration');
ylabel('Best score obtained so far');

axis tight
grid on
box on
legend('GWO')

 display(['The best solution obtained by GWO is : ', num2str(Best_pos)]);
 display(['The best optimal value of the objective funciton found by GWO is : ', num2str(Best_score)]);
% 
%         


%% function 
function fobj = F10(x)

% sensor data   vertex name kc maxl cpu ram bw
sensors = xlsread('sensors.csv');

% sink_type data   cpu ram bw cost
sink_types = xlsread('sink_types.csv');


% sinks   vertex x y 
sinks= xlsread('sinks.csv');
sink_vertex = min(sinks)+1;

for i = 1 : length(sinks)
sink_initial_type(i) = [x(i)];
end

% wsn
wsn = xlsread('graph.csv');
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


% y matrix shortest path for each sink 
Y = zeros(sink_vertex-1, length(sinks) );
for i = 1 :  sink_vertex -1
    for j = 1 : nodes - sink_vertex + 1  
            if (   (length( shortestpath(G,i,j + sink_vertex-1) )  - 1 )  <= sensors(i,4) )
                    Y(i,j) = 1 ;
                    if sink_initial_type(j) == 1
                        Y(i,j) = 0;
                    end
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
E_elec = 0;
K_E_elec = 10;
E_amp = 1;
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

%%%%%%%%%%%%%%%%%%%%%%%    communnication      %%%%%%%%%%%%%%%%%%%%%%%%%%%%
ps_communication = zeros(1,length(sinks)) ;


%%%%%%%%%%%%%%%%%%%%%%%    end      %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PS = sum(ps_process) + sum(ps_communication);

%% cost 
cost = 0 ;    
cost = sum(sink_types(sink_initial_type,4));


%% Fitness Function 
alpha = 1 ;
beta = 1000 ;
fobj = alpha*(cost + PS )  + beta*(p1 + p2 + p3 + p4 )  ;
end