
%%
tic
SearchAgents_no=4;
Function_name='F10';
Max_iteration=60; 
 lb=1;
 ub=4;
 dim=length(ssinks);
step = 1;
  format longG 
 %% fitness function
fobj = @F10;


%% algorithm
 [Best_score,Best_pos,WOA_cg_curve]=WOA(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
toc 
 %% plot 
figure('Position',[500 500 660 290])
%Draw search space
subplot(1,2,1);
x=lb:step:ub; y=x;%[-500,500]
L=length(x);
f=[];

% for i=1:L
%     for j=1:L
%            f(i,j)=fobj([x(i),y(j)]);
%      end
%       
%  end          
% 
% surfc(x,y,f,'LineStyle','none');
% 
% 
% title('Parameter space')
% xlabel('x_1');
% ylabel('x_2');
% zlabel([Function_name,'( x_1 , x_2 )'])
% %Draw objective space
subplot(1,2,2);
semilogy(WOA_cg_curve,'Color','r')
title('Objective space')
xlabel('Iteration');
ylabel('Best score obtained so far');

axis tight
grid on
box on
legend('WOA')

display(['The best solution obtained by WOA is : ', num2str(Best_pos)]);
display(['The best optimal value of the objective funciton found by WOA is : ', num2str(Best_score)]);

        

%% function 
function fobj = F10(x)
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
beta = 1000000 ;
fobj = alpha*(PS   )  + beta*(p1 + p2 + p3 + p4 )  ;
end