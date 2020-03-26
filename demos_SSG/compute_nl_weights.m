function [W, Y, neigh, t] = compute_nl_weights(pt)
% 2013-5-29
h = pt.h;
nwin = pt.nwin; nbloc = pt.nbloc;
nbest = pt.nbest; nnear = pt.nnear;
a = pt.a;
rows = pt.rows; cols = pt.cols;
X = pt.X; 
dist = pt.dist; % 0 SAD, 1 ED

m = 2*nwin+1; w = 2*nbloc+1;
G = fspecial('gaussian', [m m], a);
weight_pars = [rows; cols; m; w; nnear; nbest; h;];

if dist == 0,
    d = 1./sqrt((sum(X.^2)+eps));
    X = bsxfun(@times, X, d);
    t0 = tic;
    [W,Y] = mex_NLW_SAD(single(X), single(G), single(weight_pars));
    t = toc(t0);
else
    t0 = tic;
     [W,Y] = mex_NLW_ED(single(X), single(G), single(weight_pars));
     t = toc(t0);
end


neigh = nnear + nbest;
end