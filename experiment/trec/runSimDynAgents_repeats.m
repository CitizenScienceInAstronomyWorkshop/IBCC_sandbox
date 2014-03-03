nRepeats = 3;
nWorkers = 20; 

x_0 = cell(1, nRepeats);
aucs_dyn0 = cell(1, nRepeats);
lam_dyn0 = cell(1, nRepeats);

x_10 = cell(1, nRepeats);
aucs_dyn10 = cell(1, nRepeats);
lam_dyn10 = cell(1, nRepeats);

x_25 = cell(1, nRepeats);
aucs_dyn25 = cell(1, nRepeats);
lam_dyn25 = cell(1, nRepeats);

x_50 = cell(1, nRepeats);
aucs_dyn50 = cell(1, nRepeats);
lam_dyn50 = cell(1, nRepeats);

x_75 = cell(1, nRepeats);
aucs_dyn75 = cell(1, nRepeats);
lam_dyn75 = cell(1, nRepeats);

for repeat=1:nRepeats

%     display(['repeats: ' num2str(repeat)]);
%     nDegrade = 0;
% 
%     runSimDynAgents
%     
%     x_0{repeat} = x;
%     aucs_dyn0{repeat} = aucs;
%     lam_dyn0{repeat} = lam;    
% 
%     nDegrade = 2;
% 
%     display(['repeats: ' num2str(repeat)]);
% 
%     runSimDynAgents
%     
%     x_10{repeat} = x;    
%     aucs_dyn10{repeat} = aucs;
%     lam_dyn10{repeat} = lam;    

    nDegrade = 5;

    display(['repeats: ' num2str(repeat)]);

    runSimDynAgents
    
    x_25{repeat} = x;    
    aucs_dyn25{repeat} = aucs;
    lam_dyn25{repeat} = lam;    

    nDegrade = 10;

    display(['repeats: ' num2str(repeat)]);

    runSimDynAgents
    
    x_50{repeat} = x;    
    aucs_dyn50{repeat} = aucs;
    lam_dyn50{repeat} = lam;        

    nDegrade = 15;

    display(['repeats: ' num2str(repeat)]);

    runSimDynAgents
    
    x_75{repeat} = x;    
    aucs_dyn75{repeat} = aucs;
    lam_dyn75{repeat} = lam;    
end

fclose('all');