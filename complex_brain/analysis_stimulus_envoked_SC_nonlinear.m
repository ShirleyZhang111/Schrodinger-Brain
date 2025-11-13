%% Nonlinear structure-informed model
SC_path = 'your_path_of_SC'; 
load(SC_path);
d = sum(SC,2);
SC = diag(d.^(-0.5))*SC*diag(d.^(-0.5)); % normalization of SC
[node_items,~] = size(SC);

% complex-valued model parameter (individual)
dave_path = 'your_save_path_of_SC_nonlinear_model'; % save_path from calc_coefficients_SC_nonlinear.m
file_name = 'your_file_name';
load([dave_path file_name(4:end-4) '_SC_nonlinear_fitting_complex.mat']); % complex-valued

% stimulus = 10 on V1_L, region-th, (HCP360 atlas: region = 181)
state0 = zeros(node_items,1);
state0(region,1) = 10;
z0 = [real(state0) ; imag(state0)];
t0 = 0:0.01:9.99; % 10TRs
[T,Z] = ode89(@(t,z) complexODE(t,z,coe,SC), t0, z0);
S_construct = Z(:,1:node_items) + Z(:,node_items+1:node_items*2).*1i;
S_construct = S_construct.';

% time: response time; 
S0 = abs(S_construct);
for id_cortex = 1:node_items
    if ismember(id_cortex,region) == 0
        [pks1,locs1] = findpeaks(S0(id_cortex,:));
        if(isempty(locs1)==0)
            time_complex(id_cortex,1) = locs1(1);
        else
            time_complex(id_cortex,1) = 1000;
        end
    end
end
time_complex(region,1) = 0;

% value: response amplitude; 
value_complex = zeros(node_items,1);
for id_cortex = 1:node_items
    if(ismember(id_cortex,region) == 0)
        value_complex(id_cortex,1) = S0(id_cortex,time_complex(id_cortex));
    end
end
value_complex(region,1) = 10;

label = sortrows(cell2mat(LabelID(:,[2,5])),1);
[~,G0,~,~,~] = calc_network_measures(SC,label);
coup_strength = G0(:,region);  % Communicability of specific region 


% Response amplitude - Communicability
Amp_com = [coup_strength,value_complex];
Amp_com = sortrows(Amp_com,1);
figure;
loglog(Amp_com(:,1),Amp_com(:,2),'LineWidth',0.7,'Color',[150  59  121]/255);hold on;
xlim([min(coup_strength),max(coup_strength)])

% Response time - Communicability
Time_com = [coup_strength,time_complex/100];
Time_com = sortrows(Time_com,1);
figure;
semilogx(Time_com(:,1),Time_com(:,2),'LineWidth',1,'Color',[150  59  121]/255);hold on;
xlim([min(coup_strength),max(coup_strength)])
