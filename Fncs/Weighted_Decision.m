% Weight Decsition ,代替硬判决；
function S_wd=Weighted_Decision(s)
M=4;
s_i=real(s);
s_q=imag(s);

s_i_hd=hard_decision_pam(M,s_i);
E_i=s_i_hd-s_i;

s_q_hd=hard_decision_pam(M,s_q);
E_q=s_q_hd-s_q;


gamma_i=reliability(s_i);
gamma_q=reliability(s_q);

a=5;
b=0.5;
f_i=reliability_nonlinear(gamma_i,a,b);
f_q=reliability_nonlinear(gamma_q,a,b);


%WD
s_i_wd=s_i+(f_i.').*(E_i);

s_q_wd=s_q+(f_q.').*(E_q);


S_wd=s_i_wd+1j*s_q_wd;

end
