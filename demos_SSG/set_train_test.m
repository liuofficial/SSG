function [Train, Test, Back, Ground] = set_train_test(train_idx, test_idx, img, img_gt)
% 20130609
img_size = length(img_gt); img_idx = 1 : img_size; 

Train.idx = train_idx; Test.idx = test_idx;
Train.dat = img(:, train_idx); Test.dat = img(:, test_idx);
Train.lab = img_gt(train_idx)'; Test.lab = img_gt(test_idx)';
train_size = length(train_idx); test_size = length(test_idx);

back_idx = true(img_size,1); back_idx(train_idx) = false; back_idx(test_idx) = false;
back_idx = img_idx(back_idx);
Back.idx = back_idx;
Back.dat = img(:, back_idx);
Back.lab = 0 .* back_idx;
back_size = length(back_idx);

ground_idx = true(img_size,1); ground_idx(train_idx) = false;
ground_idx = img_idx(ground_idx);
Ground.idx = ground_idx;
Ground.dat = img(:, ground_idx);
Ground.lab = 0 .* ground_idx;
ground_size = length(ground_idx);

Train.size = train_size; Test.size = test_size;
Back.size = back_size; Ground.size = ground_size;
end