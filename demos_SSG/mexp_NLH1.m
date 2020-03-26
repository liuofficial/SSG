
% mexp_NLH1
% ||u-im0|| + mu||HW(u)||
% u0 = mexp_NLH1(single(u0), single(im0), single(W), int32(Y), single(params), nPal)
% params = [rows; cols; iNbNeigh; mu; inner]
% u0 and im0 : rows*cols # 0
% [rows cols] = size(image)
% W: nonlocal weights
% Y: nonlocal loaction
% iNbNeigh: neighbors 4 + (>=6)
% inner : iteration number of gauss-seidel, 1 or 2 is good in iteration
% nPal : 1 parallel, 0 not
% NOte: mu = R+ is better
% Author: JianJun Liu
% Date: 2013-3-26