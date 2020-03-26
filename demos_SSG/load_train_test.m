function [trainidx testidx] = load_train_test(dat, type, i, j)
switch dat,
    case 1,
        if type == 1,
            load(['data\indian\train_test\Indian_n' num2str(i) '_' num2str(j) '.mat']);
        end
        if type == 2,
            load(['data\indian\train_test\Indian_k' num2str(i) '_' num2str(j) '.mat']);
        end
        if type == 3, % percent
            load(['data\indian\train_test\Indian_r' num2str(i) '_' num2str(j) '.mat']);
        end
    case 2,
        if type == 1,
            load(['data\pavia\train_test\Pavia_n' num2str(i) '_' num2str(j) '.mat']);
        end
        if type == 2, % contain 40002, training sample, random sampling
            load(['data\pavia\train_test\Pavia_g' num2str(i) '_' num2str(j) '.mat']);
        end
        if type == 3, % contain 43923, random sampling
            load(['data\pavia\train_test\Pavia_t' num2str(i) '_' num2str(j) '.mat']);
        end
    case 3,
        if type == 1,
            load(['data\ksc\train_test\KSC_n' num2str(i) '_' num2str(j) '.mat']);
        end
        if type == 2,
            load(['data\ksc\train_test\KSC_k' num2str(i) '_' num2str(j) '.mat']);
        end
         if type == 3,
            load(['data\ksc\train_test\KSC_m' num2str(i) '_' num2str(j) '.mat']);
         end
    case 4,
        if type == 1,
            load(['data\center\train_test\Pavia_n' num2str(i) '_' num2str(j) '.mat']);
        end
        if type == 2,
        end
end
trainidx = train_idx; testidx = test_idx;
end