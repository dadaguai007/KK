function out=Filter_Work(input_date,fc,fs,fil_type)
% 使用matlab自带的各种类型的滤波器设计和作用函数

% fc 为 向量 或 数 格式，
if strcmp(fil_type,'band')
    [out,d] = bandpass(input_date,fc,fs);
    out(1:length(d.Coefficients)+1) = [];
    out(length(out)-length(d.Coefficients):length(out))=[];

elseif strcmp(fil_type,'low')
    [out,d] = lowpass(input_date,fc,fs);
    out(1:length(d.Coefficients)+1) = [];
    out(length(out)-length(d.Coefficients):length(out))=[];
elseif strcmp(fil_type,'high')
    [out,d] = highpass(input_date,fc,fs);
    out(1:length(d.Coefficients)+1) = [];
    out(length(out)-length(d.Coefficients):length(out))=[];

elseif strcmp(fil_type,'notch')
    [out,d] = bandstop(input_date,fc,fs);
    out(1:length(d.Coefficients)+1) = [];
    out(length(out)-length(d.Coefficients):length(out))=[];

end

end