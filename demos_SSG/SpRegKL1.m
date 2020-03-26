function S = SpRegKL1(AtX, AtA, opt)
% kernel sparse regression 
% date: 2013-5-22

% global set for stopping condition
global itol; global istop;
itol = 1; istop = 0;

[J, I] = size(AtX);

% param set
if ~isfield(opt,'lam'), lam = 1e-4; else lam = opt.lam; end
if ~isfield(opt,'mu'), mu = 1e-3; else mu = opt.mu; end
if ~isfield(opt,'tol'), tol = 1e-3; else tol = opt.tol; end
if ~isfield(opt,'maxit'), maxit = 500; else maxit = opt.maxit; end
if ~isfield(opt, 'verb'), verb = 1; else verb = opt.verb; end
if ~isfield(opt, 'F'), F = (AtA+mu*eye(J)) \ eye(J); else F = opt.F; end

% inital
s = zeros(J, I); d = s;
S = F * AtX;
% algorithm
for iter = 1 : maxit,
    S0 = S;
    S = F * (AtX + mu*(s+d));
    s = soft(S-d, lam/mu);
    d = d - (S - s);
    stop = print_stop(iter, S, S0, tol, verb);
    if stop, break; end
end
% show and exit
if verb > 0, disp(['SpRegKL1 Iterations:', num2str(iter)]); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = soft(u, a)
z = sign(u) .* max(abs(u)-a, 0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stop = print_stop(iter, S, S0, tol, verb)
global itol; global istop;
stop = 0;
if iter < 2, return; end
dtol = sqrt(sum(sum((S-S0).^2)) / sum(sum(S0.^2)));
ddtol = abs((dtol-itol) / itol);
if verb == 2, disp([dtol ddtol]); end
if dtol < tol, stop = 1; end
if ddtol < 0.009,
    istop = istop + 1;
    if istop > 3, stop = 2; end
else istop = 0; 
end
itol = dtol;
if rem(iter,10)==0, if verb == 1, fprintf('.'); end; end
end