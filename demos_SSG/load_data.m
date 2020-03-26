function [img img_gt rows cols bands] = load_data(n)
% 20130603
switch n
    case 1,
        load data\indian\Indian_gt.mat;
        img_gt = indian_pines_gt;
        load data\indian\AVIRIS_Indiana_imgreal.mat;
        [img, rows, cols, bands] = d3_to_d2(img, 1);
        img_gt = img_gt(:);
        rdown = 0.001; rup = 0.999;
        img = line_dat(img, rdown, rup);
        img([104:108 150:163 220],:) = [];
        return;
    case 2,
        load data\pavia\PaviaU_gt.mat;
        img_gt = double(paviaU_gt);
        load data\pavia\PaviaU.mat;
        img = paviaU;
        if m > 0,
            load data\pavia\PaviaU_im.mat;
            img = im;
            [rows cols] = size(img_gt);
            bands = size(img,1);
            img_gt = img_gt(:);
            rdown = 0.001; rup = 0.999;
            img = line_dat(img, rdown, rup);
            return;
        end
    case 3,
        load data\ksc\KSC.mat;
        img_gt = KSC_gt; img = KSC;
        img_gt = img_gt(53:end,155:end); img = img(53:end,155:end,:);
    case 4,
        load data\center\Pavia_gt.mat;
        img_gt = double(pavia_gt); img_gt = img_gt(:,224:end);
        load data\center\Pavia.mat;
        img = pavia; img = img(:,224:end,:);
end

[img, rows, cols, bands] = reshape3d_2d(img);
img_gt = img_gt(:);
rdown = 0.001; rup = 0.999;
img = line_dat(img, rdown, rup);
end