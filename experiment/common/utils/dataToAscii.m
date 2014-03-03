
dataDirectory = '/homes/49/edwin/matlab/combination/gzsnBundle/data';

%filename for the dataset should be here. Requires the following fields:
%1. classification ID
%2. agent ID
%3. asset ID
%4. PTF type
%5. PTF class
filename = [dataDirectory '/gzsnData.mat'];

if loadDataFromFile && (~exist('snRawData', 'var') || isempty(snRawData))
    load(filename);
    snRawData = snData;
end

fid = fopen([dataDirectory '/gzsnData.csv'], 'w');
for n=1:length(snRawData{1})
    fprintf(fid, '%d, %d, %d, %s, %s, %d\n', snRawData{2}(n), snRawData{3}(n), ...
        snRawData{4}(n), snRawData{7}{n}, snRawData{8}{n}, snRawData{9}(n) );
    display(num2str(n));
end
fclose(fid);

clear snData;
