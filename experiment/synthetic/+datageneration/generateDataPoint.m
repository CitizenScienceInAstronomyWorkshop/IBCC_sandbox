function [ row, t ] = generateDataPoint( nSensors, nInfSensors, means, deviations, p_c1, pMode, class1ErrorOnly )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
   
   if ~exist('p_c1','var')
       p_c1 = 0.5;
   end

   row = randn(1, nSensors);

   c = rand(1);
   m = rand(1, nInfSensors);
   
   if ~exist('pMode','var')
        pMode = 0.8; 
   end
%    meanDiff = means(1, 1:nInfSensors) - means(2, 1:nInfSensors);  

   if c > p_c1
       %Class 0 object
        if exist('class1ErrorOnly','var') && class1ErrorOnly==true
            pMode = 0;
        elseif exist('class1ErrorOnly','var') && class1ErrorOnly==2
            pMode = 1 - pMode;
        end
        muWrong = means(1, 1:nInfSensors); %class 1 means
        muRight = means(2, 1:nInfSensors); % class 0 means
        mu = muWrong.*(m<pMode) + muRight.*(m>=pMode);
%         row(1:nInfSensors) = row(1:nInfSensors) .* deviations(2, 1:nInfSensors) + mu;
        row(1:nInfSensors) = normrnd(mu, abs(deviations(2, 1:nInfSensors)));

        t = 0;
   else
       %Class 1 object
        muWrong = means(2, 1:nInfSensors); 
        muRight = means(1, 1:nInfSensors);
        mu = muWrong.*(m<pMode) + muRight.*(m>=pMode);
%         row(1:nInfSensors) = row(1:nInfSensors) .* deviations(1, 1:nInfSensors) + mu;
        row(1:nInfSensors) = normrnd(mu, abs(deviations(1, 1:nInfSensors)));

        t = 1;
   end   
   
end

