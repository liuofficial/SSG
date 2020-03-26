function demos
n = 1;
[img img_gt rows cols bands] = load_data(n);
nClass = length(unique(img_gt)) - 1;
gams = 2.^0;  mus = 10.^(-3); lams = 10.^(-4); sigs = sqrt(0.5 ./ gams);

[train_idx, test_idx] = load_train_test(n, 2, 3, 1);
[Train, Test] = set_train_test(train_idx, test_idx, img, img_gt);
sig = sigs; mu = mus; lam = lams;

opt = [];
opt.mu = mu; opt.lam = lam;
[S, pred, t_alg, t_prd] = KSR(Train.dat, img, sig, 10000, opt, img_gt(:)', nClass, Train.lab);
acc = class_eval(pred(test_idx), Test.lab);
disp(acc.OA);

beta = 20;
NL = [];
NL.nwin = 2; NL.nbloc = 9;
NL.nbest = 6; NL.nnear = 8;
NL.a = 2;
NL.X = img;
NL.rows = rows; NL.cols = cols;
NL.dist = 0;
snr = estimate_snr(img, nClass);
NL.h = sqrt(norm(img,'fro')^2/numel(img)/10^(snr/10)) * 2 * pi;
[NL.W, NL.Y, NL.neigh, NL.time] = compute_nl_weights(NL);

P = S;
sci_idx = SCI(P, 0.7, Train.lab, nClass);
opt = [];
opt.beta = beta; opt.idx = sci_idx; opt.Pt = P(:,sci_idx); opt.imsize = [rows cols]; opt.maxit = 50;

t0 = tic;
S = SpDenNLH1(P, opt, NL);
t_alg = toc(t0);
[pred, t_prd] = CLS_PRED(Train.dat, img, S, sig, 10000, img_gt(:)', nClass, Train.lab);
acc = class_eval(pred(test_idx), Test.lab);
disp(acc.OA);
end

function sci_idx = SCI(S, tol, train_lab, nClass)
T = coef_pred(train_lab);
sci = max(T * abs(S));
sci = bsxfun(@times, sci, 1./sum(abs(S)));
sci = (nClass * sci - 1) / (nClass - 1);
sci_idx = sci > tol;
end

function [pred, t_prd] = CLS_PRED(A, X, S, sig, bloc, gt, nClass, lab)
I = size(X,2);
pred = zeros(nClass, I);
AtA = RBF(A, A, sig, 10000);
nbloc = ceil(I / bloc);
t_prd = 0;
for i = 1 : nbloc,
    idx = (i-1)*bloc+1 : min(i*bloc, I);
    AtX = RBF(A, X(:,idx), sig, 10000);
    t0 = tic;
    pred(idx) = class_ker_pred(AtX, AtA, S(:,idx), gt(idx), lab);
    t_prd = t_prd + toc(t0);
end
end

function [S, pred, t_alg, t_prd] = KSR(A, X, sig, bloc, opt, gt, nClass, lab)
J = size(A,2); I = size(X,2);
S = zeros(J, I);
pred = zeros(nClass, I);
t0 = tic;
AtA = RBF(A, A, sig, 10000);
opt.F = (AtA+opt.mu*eye(J)) \ eye(J);
t_alg = toc(t0);
t_prd = 0;
nbloc = ceil(I / bloc);
for i = 1 : nbloc,
    idx = (i-1)*bloc+1 : min(i*bloc, I);
    t0 = tic;
    AtX = RBF(A, X(:,idx), sig, 10000);
    S(:,idx) = SpRegKL1(AtX, AtA, opt);
    t_alg = t_alg + toc(t0);
    t0 = tic;
    pred(idx) = class_ker_pred(AtX, AtA, S(:,idx), gt(idx), lab);
    t_prd = t_prd + toc(t0);
end
end