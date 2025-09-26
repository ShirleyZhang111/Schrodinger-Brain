function corr_mode = calc_recon_with_modes(inCfg)
% Data reconstruction using schrodinger-derived/SC-derived modes
% Input: 
%     mode : structural-eigenmode or spatiotemporal mode
%     TC : space * time
% Output:
%     corr_mode : the peason's correlation of the reconstructed FC (calculated by the
%     reconstructed data) and the empirical FC 

% input configuration
mode = inCfg.mode; % structural-eigenmode or spatiotemporal mode
data_path = inCfg.data_path;
file_name = inCfg.file_name;

% start processing
load([data_path file_name],'TC');% N*T
FC = corrcoef(TC');
corr_mode = zeros(length(TC(:,1)),1);
% data reconstruction with len modes
for len = 1:length(TC(:,1))
    U = mode(:,1:len);     
    C = U'* TC;
    signal_recon = U*C;
    FC2 = corrcoef(signal_recon');
    [r2,~] = corr(FC(:),FC2(:));
    corr_mode(len,1) = r2;
end