data_path = 'your_save_path_of_nonlinear_model'; %  save_path in calc_coefficients_nonlinear.m
items = dir(data_path);
items = items(3:end);
% Group parameters
subjects_number = length(items);
for id = 1:subjects_number
    item = items(id).name;
    hopf = load(strcat(data_path, '/', item));    
    Coefficient(:,:,id) = hopf.coefficient; 
    Coupling_mat(:,:,id) = hopf.coupling_mat;
    Energy(:,:,id) = hopf.energy;
    MSE(:,:,id) = hopf.mse_real;
end

% Nonzero-pattern of Coupling Matrix (Individual)
figure; spy(Coupling_mat(:,:,1));

% Angles of Nonzeros of Complex-valued Coupling Matrix (Group)
figure;
polarhistogram(angle(Coupling_mat(Coupling_mat~=0)), 100,'FaceColor',[237 141 90]/255,'FaceAlpha',1,'EdgeColor',[1.00,0.80,0.16]); 
set(gcf,'unit','centimeters','position',[0.1 0.1 6 3.7])
set(gca,'FontUnits','points','FontName','Arial',"FontSize",7);

% Fitting accuracy (MSE) of nonlinear model (Group)
figure;
edges = [0 0.5:0.03:1 1];
histogram(MSE,edges,'FaceColor',[189 119 149]/255,'FaceAlpha',1); hold on;
xline(mean(mean(MSE)),'Color',[150 59 121]/255,'LineWidth',1.5,'Alpha',1);
set(gca, 'Box', 'off');
set(gcf,'unit','centimeters','position',[0.1 5 7 4.8]);
set(gca,'FontUnits','points','FontName','Arial',"FontSize",7);
xlabel('Correlation','FontUnits','points','FontName','Arial',"FontSize",7);
ylabel('Frequency','FontUnits','points','FontName','Arial',"FontSize",7)


% Angles of Eigenvalues of Complex-valued Coupling Matrix (Group)
lambda = zeros(subjects_number,1);
for k = 1:subjects_number
    [E,D] = eig(Coupling_mat(:,:,k));
    lambda(:,k) = diag(D);
end
figure;
polarhistogram(angle(lambda(lambda~=0)), 100,'FaceColor',[237 141 90]/255,'FaceAlpha',1,'EdgeColor',[1.00,0.80,0.16]); 
set(gcf,'unit','centimeters','position',[0.1 0.1 6 3.7])
set(gca,'FontUnits','points','FontName','Arial',"FontSize",7);

% Angles of Energy (Group)
figure;
polarhistogram(angle(Energy), 100,'FaceColor',[128 197 162]/255,'FaceAlpha',1,'EdgeColor',[128 197 162]/255); 
set(gcf,'unit','centimeters','position',[0.1 0.1 3.8 3.8])
set(gca,'FontUnits','points','FontName','Arial',"FontSize",7);

% Amplitude of nonlinear polynomial coefficient (Group)
Amp = abs(Coefficient);
figure;
data  = Amp((7:10),:)';
hb = boxplot(data,...              
                    'Color','k',...                                   
                    'symbol',' ',...                                  
                    'Notch','on',...                                 
                    'OutlierSize',4);                              
set(hb,'LineWidth',1)                                               
colors = [151 208 197;151 208 197;45 136 117;151 208 197]/255;
h = findobj(gca,'Tag','Box');
for j = 1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.9);
end
ylim([0 0.0015]);box on;set(gca,'FontSize',10);set(gca, 'LineWidth', 0.5);

% Angles of global coupling parameter 'g' in SC-informed nonlinear model (Group)
data_path_sc = 'your_save_path_of_SC_nonlinear_model'; %  save_path in calc_coefficients_SC_nonlinear.m
items = dir(data_path_sc);
items = items(3:end);
for id = 1:subjects_number
    item = items(id).name;
    hopf = load(strcat(data_path_sc, '/', item));    
    g(:,id) = hopf.g;
end
figure;
histogram(angle(g),20,'FaceAlpha',1,'FaceColor',[237 141 90]/255,'EdgeColor','none');hold on;
xline(mean(angle(g)),'LineWidth',2,'Color',[0.84,0.34,0.07],'Alpha',1);
