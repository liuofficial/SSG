function P = SpDenNLH1(P, opt, NL)
% Author:  JianJun Liu
% Date:     2013-06-15

global itol; global istop;

if ~isfield(opt, 'beta'), beta = 20; else beta = opt.beta; end
if ~isfield(opt, 'maxit'), maxit = 500; else maxit = opt.maxit; end
if ~isfield(opt, 'tol'), tol = 1e-3; else tol = opt.tol; end
if ~isfield(opt, 'verb'), verb = 1; else verb = opt.verb; end

imsize = opt.imsize; rows = imsize(1); cols = imsize(2);
idx = opt.idx; Pt = opt.Pt;

nlh1_pars = [rows; cols; NL.neigh; beta; 1];

itol = 1; istop = 0;
for iter = 1 : maxit,
    P0 = P;
    P = mexp_NLH1(single(P)', single(P)', single(NL.W), int32(NL.Y), single(nlh1_pars), 1)';
    P(:,idx) = Pt;
    if print_stop(iter, P, P0, tol, verb), break; end
end
if verb > 0, disp(strcat(['SpDenNLH1 Iterations:', num2str(iter)])); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stop = print_stop(iter, Q, Q0, tol, verb)
global itol; global istop;
stop = 0;
if iter < 2, return; end
dtol = sqrt(sum(sum((Q-Q0).^2)) / sum(sum(Q0.^2)));
if verb == 2, disp([dtol abs((dtol-itol) / itol)]); end
if dtol < tol, stop = 1; end
if abs((dtol-itol) / itol) < 0.009,
    istop = istop + 1;
    if istop > 3, stop = 1; end
else istop = 0;
end
itol = dtol;
if rem(iter,10)==0, if verb == 1, fprintf('.'); end; end
end

function z = soft(u, a)
z = sign(u) .* max(abs(u)-a, 0);
end