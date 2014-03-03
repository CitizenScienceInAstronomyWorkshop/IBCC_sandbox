vbData = dlmread('/homes/49/edwin/matlab/combination/data/vbIbccPaper/exp2_inv/baseData/nIts_vb');
emData = dlmread('/homes/49/edwin/matlab/combination/data/vbIbccPaper/exp2_inv/baseData/nIts_em');
gibbsData = dlmread('/homes/49/edwin/matlab/combination/data/vbIbccPaper/exp2/baseData/nIts_gibbs');

data = [vbData' emData' gibbsData'];
figure;
hist(data);

means = mean(data, 1);
vars = var(data, 0, 1);