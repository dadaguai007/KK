% OFDM_tran
sigRxo=ssfm(sigTxo,param);
power2=signalpower(sigRxo);
fprintf(' after ssfm signal power: %.2f dBm\n', 10 * log10(power2 / 1e-3));
