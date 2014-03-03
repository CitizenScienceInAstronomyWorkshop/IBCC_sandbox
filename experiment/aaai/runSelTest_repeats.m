nRepeats = 1;
aucs_set = cell(1, nRepeats);
lam_set = cell(1, nRepeats);

for repeat=1:nRepeats
    display(['repeats: ' num2str(repeat)]);

    runSelTest
    
    aucs_set{repeat} = aucs;
    lam_set{repeat} = lam;    
end

fclose('all');