function [numClassData, numTypeData] = translateGZSNLabels(objIds, labelText, screenText)

%Translate the PTF labels into scores - only -1, 1 and 3 are allowed.
classCounts = zeros(length(objIds), 3); %record the number of PTF entries supporting each classification
typeCounts = zeros(length(objIds), 2);
numTypeData = zeros(length(objIds), 1); %count the number of PTF entries saying an asset is a transient
numClassData = zeros(length(objIds), 1); %count the number of PTF entries saying an asset is a SN
nClassOnly = 0; %number for which we have no PTF type 

for i=1:length(objIds)
    
    % TYPE (SCREEN LABELS) ------------------------------------------------
    type = screenText(i);
    type = type{1};
    tokens = textscan(type, '%s');
    if strcmp(type, '""')
        tokens = {[]};
    end
    for t=1:length(tokens{1})
        token = tokens{1}{t};
        
        if strcmp(token, '"')
            continue
        end
        
        if ~isempty(strfind(token, 'Transient'))
            typeCounts(i, 2) = typeCounts(i, 2) + 1;
        else
%             display(['token: ' token])
            typeCounts(i, 1) = typeCounts(i, 1) + 1;
        end
    end
    
    %CLASS (TRUE LABELS) --------------------------------------------------
    
    %ptfclass (e.g. supernova type ii/a)
    class = labelText(i);
    class = class{1};

    tokens = textscan(class, '%s');
    for t=1:length(tokens{1})
        token = tokens{1}{t};       
        if strcmp(token, '""')
            %do nothing wid dis            
        elseif ~isempty(strfind(token, 'SN')) || ...
                ~isempty(strfind(token, 'SN?')) || ...
                ~isempty(strfind(token, 'nova')) || ...
                ~isempty(strfind(token, 'SNI')) || ...
                ~isempty(strfind(token, 'Ia')) || ...
                ~isempty(strfind(token, 'Ic')) || ...
                ~isempty(strfind(token, 'II')) || ...
                ~isempty(strfind(token, 'IIb')) %|| ...
            %Positive Examples
            classCounts(i, 2) = classCounts(i, 2) + 1;  
        elseif ~isempty(strfind(token, 'AGN')) ||...%AGN = active galactic nucleus
                ~isempty(strfind(token, 'CV')) || ...%CV = cataclysmic variable star
                ~isempty(strfind(token, 'galaxy')) || ... 
                ~isempty(strfind(token, 'Junk')) || ... 
                ~isempty(strfind(token, 'CV?')) || ... 
                ~isempty(strfind(token, 'Rock')) ||...
                ~isempty(strfind(token, 'rock')) ||...%might need to treat as unknown as rock classification is unreliable
                ~isempty(strfind(token, 'type')) 
            %Negative Examples
            classCounts(i, 1) = classCounts(i, 1) + 1;
        elseif ~isempty(strfind(token, 'unknown'))||...
            ~isempty(strfind(token, 'Maybe')) 
            %Maybe?
            classCounts(i, 3) = classCounts(i, 3) + 1;
        elseif ~isempty(strfind(token, 'Transient')) 
            %Maybe - add to type category
            typeCounts(i,2) = typeCounts(i,2) + 1;
        elseif  ~isempty(strfind(token, 'varstar')) ||...
                ~isempty(strfind(token, 'VarStar')) 
            %Probably not - put this in the type category
            typeCounts(i,1) = typeCounts(i,1) + 1;
        else
            %ignore other stuff but print it in case we missed something.
%             display(['Other string found in PTF class: ' token]);
        end
    end
    % PROCESS TYPE LABELS -------------------------------------------------
    [maxScore, maxIdx] = max(typeCounts(i,:));
    maxFlags =  maxScore>0 & typeCounts(i, :) == maxScore;
    % Non-transient -> -1, transient -> 1, >1 transient labellings -> 3
    if sum(maxFlags) == 1 %One answer was most popular
        if maxIdx==1
            numTypeData(i) = -1;
        elseif maxScore>1
            numTypeData(i) = 3;
        end
    else %equal number of results for both
        if maxScore > 0
%             display(['altering ptftype interpretation: ' num2str(typeCounts(i,1)) ', ' num2str(typeCounts(i,2))]);
            numTypeData(i) = 1;
        else %no labels
            numTypeData(i) = 0;
        end
    end    
    % PROCESS CLASS LABELS ------------------------------------------------
    [maxScore, maxIdx] = max(classCounts(i,1:2));
    maxFlags = classCounts(i, :) == maxScore;
    if maxScore == 0
        %no positive/negative confirmations
        if classCounts(i,3)>0
            numClassData(i) = -1;
        else
            numClassData(i) = 0;
        end
    elseif sum(maxFlags) == 1
%         display(maxScore);
        if maxIdx==2
            numClassData(i) = 2;
        elseif maxIdx==1
            numClassData(i) = 1;
        end
    else
        numClassData(i) = -1;
    end
    
    if numTypeData(i) == 0 && maxScore > 0
        nClassOnly = nClassOnly + 1;
    end
end