function dirName = checkDataDir( expSettings )
    dirName = expSettings.getDataDir();

    while expSettings.noOverwriteData && exist(dirName, 'file')==7
        dirName = sprintf('%s0', dirName);
    end

    mkdir(dirName);
end

