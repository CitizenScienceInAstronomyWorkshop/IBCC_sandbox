function AUC=singlePointAUC(fpr, tpr)

topLeft = fpr.*(1-tpr);
bottomLeft = fpr.*tpr ./ 2;
topRight = (1-fpr).*(1-tpr)./2;

AUC = 1 - topLeft - bottomLeft - topRight;