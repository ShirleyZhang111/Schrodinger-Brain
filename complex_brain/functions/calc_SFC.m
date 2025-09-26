function SFC = calc_SFC(H,SC)
% Calculate structure-functional connectivity
% Input: 
%     H : functional coupling matrix
%     TC : structural connectivity matrix
% Output:
%     SFC : the structure-functional connectivity of H

% input configuration
% group functional connectivity matrix: FC/H_complex/H_real
H = H - diag(diag(H));
% group structural connectivity matrix
SC = SC - diag(diag(SC));

% compute structural-functional coupling
SFC = zeros(length(SC(:,1)),1);
for row = 1:length(SC(:,1))
    [r1,~] = corr(SC(row,:)',H(row,:)','type','pearson');
    SFC(row,1) = r1; 
end