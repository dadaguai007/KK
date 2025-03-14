function [sync_bitsequence,shift] = sync_prbs_seq(bitdet,bitref)
% bitdet: the decoded bit sequence
% bitref: the reference prbs bit sequence
% if size(bitdet) ~= size(bitref)
%     error('Size mismatch!');
% end
if isrow(bitdet)
    bitdet = bitdet';
    bitref = bitref';
end
Cor = xcorr(bitdet-mean(bitdet),bitref-mean(bitref));
% plot(Cor);
[~,peakindex] = max(abs(Cor));
shift = peakindex-length(bitdet);
sync_bitsequence = circshift(bitdet,-shift);
% if max(Cor)~=max(abs(Cor)) % to check whether it is data bar
%     sync_bitsequence = 1-sync_bitsequence;
% end