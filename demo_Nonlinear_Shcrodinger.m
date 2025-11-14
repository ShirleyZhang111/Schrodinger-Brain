% Nonlinear Schrodinger-like model
% This code implements the Nonlinear comple-valued model for fMRI data analysis

% Load rs-fMRI signals (available datasets: HCP, HCPex, HCP-voxel, UKB)
load('TC_HCP_379.mat');
% load('TC_HCPex_426.mat');
% load('TC_UKB_426.mat');
cfg.is_largescale = 0;

% For large-scale voxel-level analysis:
% load('TC_HCP_9982.mat');
% cfg.is_largescale = 1;

% Set configuration parameters for nonlinear model
cfg.TC = TC;
cfg.maxit = 2000;
cfg.tol = 1e-12;
cfg.field = 'complex';
cfg.T = 101:1100;
cfg.mu = 0.1*length(cfg.T);

% Estimate Model Parameters from rs-fMRI data using nonlinear model
[coeff, coupling_mat, res, ~, ~, energy,Cor] = calc_coefficients_nonlinear(cfg);



% Figure 1: Nonzeros pattren of sparse coupling matrix
figure; spy(coupling_mat);

% Figure 2: Distribution of sparse coupling matrix angles
% Shows the phase distribution of non-zero coupling elements in polar coordinates
figure;
polarhistogram(angle(coupling_mat(coupling_mat~=0)), 100,'FaceColor',[237 141 90]/255,'FaceAlpha',1,'EdgeColor',[1.00,0.80,0.16]); 
set(gcf,'unit','centimeters','position',[0.1 0.1 6 3.7])
set(gca,'FontUnits','points','FontName','Arial',"FontSize",7);

% Figure 3: Histogram of model fitting accuracy (correlation coefficients)
% Displays the distribution of correlation between model predictions and actual data
figure;
edges = [0 0.5:0.03:1 1];
histogram(Cor,edges,'FaceColor',[189 119 149]/255,'FaceAlpha',1); hold on;
xline(mean(mean(Cor)),'Color',[150 59 121]/255,'LineWidth',1.5,'Alpha',1);
set(gca, 'Box', 'off');
set(gcf,'unit','centimeters','position',[0.1 5 7 4.8]);
set(gca,'FontUnits','points','FontName','Arial',"FontSize",7);
xlabel('Correlation','FontUnits','points','FontName','Arial',"FontSize",7);
ylabel('Frequency','FontUnits','points','FontName','Arial',"FontSize",7)

% Figure 4: Distribution of Hamiltonian energy angles
% Shows the phase distribution of Hamiltonian energy values
figure;
polarhistogram(angle(energy), 100,'FaceColor',[128 197 162]/255,'FaceAlpha',1,'EdgeColor',[128 197 162]/255); 
set(gcf,'unit','centimeters','position',[0.1 0.1 3.8 3.8])
set(gca,'FontUnits','points','FontName','Arial',"FontSize",7);

%% SC-informed Nonlinear Model and stimulus-envoked dynamics

% Load structural connectivity matrix
load('SC_379.mat');
% Normalize SC matrix
d = sum(SC_379,2);
SC = diag(d.^(-0.5))*SC_379*diag(d.^(-0.5)); 

% Load fMRI data
load('TC_HCP_379.mat');

% Set configuration parameters
cfg.TC = TC;
cfg.SC = SC;
cfg.field = 'complex';
cfg.T = 101:1100;

% Estimate model parameters with SC constraint
[coefficient, g] = calc_coefficients_SC_nonlinear(cfg);

% Simulate stimulus-evoked dynamics
% Apply stimulus to V1_L region (region 181 in HCP360 atlas)
node_items = 379; 
region = 181;

% Set initial state with stimulus
state0 = zeros(node_items,1);
state0(region,1) = 10;          % Stimulus amplitude = 10
z0 = [real(state0) ; imag(state0)];
t0 = 0:0.01:9.99;               % Time vector for simulation (10 TRs)

% Solve nonlinear ODE system
[T,Z] = ode89(@(t,z) complexODE(z,[coefficient;g],SC), t0, z0);
S_construct = Z(:,1:node_items) + Z(:,node_items+1:node_items*2).*1i;
S_construct = S_construct.';


% Calculate response characteristics

% Response time calculation
S0 = abs(S_construct);          % Signal amplitude
for id_cortex = 1:node_items
    if ismember(id_cortex,region) == 0
        [pks1,locs1] = findpeaks(S0(id_cortex,:));
        if(isempty(locs1)==0)
            time_complex(id_cortex,1) = locs1(1); % First peak time
        else
            time_complex(id_cortex,1) = 1000;     % No peak found
        end
    end
end
time_complex(region,1) = 0;     % Stimulated region time = 0


% Response amplitude calculation  
value_complex = zeros(node_items,1);
for id_cortex = 1:node_items
    if(ismember(id_cortex,region) == 0)
        value_complex(id_cortex,1) = S0(id_cortex,time_complex(id_cortex));
    end
end
value_complex(region,1) = 10;   % Stimulated region amplitude = 10

% Communicability analysis
% Compute communicability matrix
d = sum(abs(SC),2);
W = diag(d.^(-1/2))*SC*diag(d.^(-1/2));
[E0,D0] = eig(W);
G0 = E0*diag(exp(diag(D0)))*E0';

% Get communicability from stimulated region
coup_strength = G0(:,region);  

% Plot: Response amplitude vs Communicability
Amp_com = [coup_strength,value_complex];
Amp_com = sortrows(Amp_com,1);
figure;
loglog(Amp_com(:,1),Amp_com(:,2),'LineWidth',0.7,'Color',[150  59  121]/255);hold on;
xlim([min(coup_strength),max(coup_strength)])

% Plot: Response time vs Communicability  
Time_com = [coup_strength,time_complex/100];
Time_com = sortrows(Time_com,1);
figure;
semilogx(Time_com(:,1),Time_com(:,2),'LineWidth',1,'Color',[150  59  121]/255);hold on;

xlim([min(coup_strength),max(coup_strength)])
