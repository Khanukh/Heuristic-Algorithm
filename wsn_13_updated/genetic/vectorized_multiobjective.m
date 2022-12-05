function scores = vectorized_multiobjective(t)
IntCon = [1, 1,1,1]; 
x = t ;

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
beta = 1000000 ;
fobj = alpha*(PS   )  + beta*(p1 + p2 + p3 + p4 )  ;


   %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%
  global Gvalue_num;
  global G_value;
   G_value(Gvalue_num)= fobj;
    Gvalue_num = Gvalue_num + 1;
    
   
     popSize = size(t,1); % Population size
     numObj = 1;  % Number of objectives
     % initialize scores
     scores = zeros(popSize, numObj);
     % Compute wn
%      scores(:,1) = abs(((sqrt(KKeq/M)/(2*pi)-WWn)/WWn))*100 ;
     % Compute sd
     scores(:,1) =fobj;
     % Compute BNEA
%      scores(:,3) = abs((sqrt(4*1.38e-23*300/(M * 1e3 *Q*(Wn*2* pi)^3)) / (Sd / 1e3) - BBNEA)/BBNEA)*130;
   
  