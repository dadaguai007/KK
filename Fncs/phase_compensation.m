

% W=400;
% pilotIndex=1:1:W;
phi_mean=angle(mean(data_kk(pilotIndex,:)./...
    qam_signal_mat(pilotIndex,:)));
data_kk=data_kk.*...
    repmat(exp(-1j.*phi_mean),size(data_kk,1),1);