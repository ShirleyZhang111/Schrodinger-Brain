function [mode,frequency] = calc_spatiotemporal_modes_dynamics(inCfg)
% Input: 
%     idx : the idx-th slowest mode     
%     t : temporal duration of mode dynamics
%     lambda : the eigenvalue of Q
%     G : all the modes
% Output: 
%     mode: the idx-th slowest spatiotemporal mode
%     frequency: the frequency of this mode

% input configuration
if isfield(inCfg,'idx');idx = inCfg.idx;else;idx = 1;end % The idx-th slowest mode
if isfield(inCfg,'t');t = inCfg.t;else;t = 100;end % mode dynamics temporal
lambda = inCfg.lambda;
G = inCfg.G;

phase_lambda = angle(lambda)/(2*pi*0.72);
phase_lambda = sortrows([(1:length(lambda))',phase_lambda],2);
frequency = phase_lambda(idx,2);

mode = zeros(length(lambda),t);
mode(:,1) = real(G(:,phase_lambda(idx,1)));
for k = 1:t-1
    mode(:,k+1) = real((lambda(idx,1).^k).*G(:,phase_lambda(idx,1)));
end