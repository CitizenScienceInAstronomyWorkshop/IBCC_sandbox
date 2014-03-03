alpha = zeros(6, 24);
acc = zeros(10, 5);

nAgents = max(trec7_topicIdx(:,6));
userAlpha = zeros(nAgents, 5);

nRows = size(trec7_topicIdx,1);

respClasses = unique(trec7_topicIdx(:,3));

for i=1:nRows
    row = trec7_topicIdx(i,:);
    trueClass = row(1);
    trueClass = find(respClasses==trueClass)-1;
    testClass = trueClass;
    if row(2)==0
        trueClass = 6;       
    end
    response = row(4);
    
    respClass = row(3);
    if respClass == 0 
        respClass = 6;
    else
        respClass = find(respClasses==respClass)-1;
    end        
  
    %for original dataset
%     if response == 3
%         alpha(trueClass, respClass+12) = alpha(trueClass,respClass+12) + 1;
%     elseif response == 1
%         alpha(trueClass,respClass) = alpha(trueClass,respClass) + 1;
%     else
%         alpha(trueClass, 6+respClass) = alpha(trueClass, 6+respClass)+1;
%     end

    conf = row(5);

    if response == 3 && conf
        col = respClass+12;
    elseif response==3 && ~conf
        col = respClass + 18;
    elseif response ==1 && conf
        col = respClass;
    elseif response ==1 && ~conf
        col = respClass + 6;
    else
        display('error');
    end
    alpha(trueClass, col) = alpha(trueClass, col) + 1;
    
    if trueClass==6
        testClass = testClass + 5;
    end
    if (respClass==testClass && response==1) || (respClass==6 && testClass>5) ...
            || (response==3 && testClass>5 && respClass==testClass-5)
        if conf
            col = 1;
        else
            col = 2;
        end
    elseif response==1
        if conf
            col = 3;
        else
            col = 4;
        end
    else
        col = 5;
    end
        
    acc(testClass, col) = acc(testClass, col) + 1;
    userAlpha(row(6), col) = userAlpha(row(6), col) + 1;
end

alpha
