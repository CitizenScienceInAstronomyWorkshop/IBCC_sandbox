settingsFuncs = {@settings.thesis.exp1, ...
    @settings.thesis.exp2, ...
    @settings.thesis.exp3, ...
    @settings.thesis.exp4, ...
    @settings.thesis.exp6, ...
    @settings.thesis.exp7
    };

for expNo=1:numel(settingsFunc)
    
    varySettingFunc = settingsFuncs{expNo};
    [expSet, synthSet, bccSet, ~, range, varName] = varySettingFunc(1);
    
    for i=1:numel(range)+1
        filename = [expSet.outputDir '/' num2str(i) '.fig'];
        hgload()
    end
end