function snr_est = estimate_snr(R, p)
[L, N]=size(R);     
r_m = mean(R,2);
R_m = repmat(r_m,[1 N]); % mean of each band
R_o = R - R_m;           % data with zero-mean
[Ud,~,~] = svds(R_o*R_o'/N,p);  % computes the p-projection matrix
x_p =  Ud' * R_o;                 % project the zero-mean data onto p-subspace
P_y = sum(R(:).^2)/N;
P_x = sum(x_p(:).^2)/N + r_m'*r_m;
snr_est = 10*log10( (P_x - p/L*P_y)/(P_y- P_x) );
end