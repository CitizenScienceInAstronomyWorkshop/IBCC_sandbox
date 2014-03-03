function [ rating, classRating, igRating ] = piRating( Alpha, Nu, hd )
%piRating Score the agents according to the informativeness of their
%confusion matrices. Uses sum over pairs of rows Hellinger distance between dirichlet-multinonmial
%distributions, normalised to between 0 and 1.

nAgents = size(Alpha,3);
nClasses = size(Alpha,1);

classRating = sparse(nClasses,nAgents);
rating = sparse(1,nAgents);
nonZero = find(sum(Alpha(1,:,:),2)>0);

nClasses = size(Alpha,1);
% if nargin<3 || hd
%     for i=1:nClasses
%         for j=1:nClasses
% 
%             if i==j
%                 %rating = rating + 0;
%             else
%                 %Hellinger distance between rows i and j
% 
%                 gammasumi = gammaln(sum(Alpha(i,:,nonZero),2));
%                 gammasumj = gammaln(sum(Alpha(j,:,nonZero),2));
% 
%                 gammasumin = gammaln(sum(Alpha(i,:,:),2)+1);
%                 gammasumjn = gammaln(sum(Alpha(j,:,:),2)+1);               
% 
%                 ratingij = sparse(1,nAgents);
%                 for c=1:size(Alpha,2)            
%                     fx = gammasumi-gammasumin  + gammaln(Alpha(i,c,nonZero)+1) - gammaln(Alpha(i,c,nonZero)); 
%                     gx = gammasumj-gammasumjn  + gammaln(Alpha(j,c,nonZero)+1) - gammaln(Alpha(j,c,nonZero));
% 
%     %                 ratingij(nonZero) = sqrt(1 - sqrt(exp(reshape(fx,1,nAgents) + reshape(gx,1,nAgents))));
%                     ratingijc = (sqrt(exp(reshape(fx,1,nAgents))) - sqrt(exp(reshape(gx,1,nAgents))) ).^2;
%                     ratingij(nonZero) = ratingij(nonZero) + ratingijc;
% 
%                 end
%                 rating(nonZero) = rating(nonZero) + (1./sqrt(2) .* ratingij(nonZero).^0.5);
%                 classRating(i, nonZero) = classRating(i, nonZero) + (1./sqrt(2) .* ratingij(nonZero).^0.5);
%             end
%         end
%     end
%     rating = rating ./ (nClasses*nClasses);
%     classRating = classRating ./ nClasses;
% else
%     rating = [];
%     classRating = [];
% end
if ndims(Nu)==3
    Kappa = Nu ./ repmat(sum(Nu), size(Nu,1),size(Nu,2));
    state = {Alpha, Nu, Kappa, []};
else
    Kappa = Nu ./ sum(Nu);    
    state = {Alpha, Nu, Kappa', []};
end

igRating = infgain.batchInfoGain(1, state);

end

