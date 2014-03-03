function [IG Hci Hci_t] = evaluatePair2(i, P, Pi, lnPi, Alpha0, p_casej, p_totalc, p_total)
%UNTITLED3 Evaluate the IG for a worker and document pair. Use the second
%Calculation variant
%   Should give same ordering of pairs as the first calculation but may be more
%   efficient for some cases due to rearranging the equations and avoiding
%   calculating some costly constant terms. The Values produced should
%   differ by constant amounts. 
    
    pti = repmat(P(i,:)', 1, size(Pi,2));

    Pi_t = sum( pti .* Pi, 1 );
    lnPi_t = sum( pti .* lnPi, 1 );
    Hci = -sum( Pi_t .* lnPi_t, 2);
    
    Hci_t = 0;
    
    nScores = size(Pi,2);
    
    Alpha0_total = repmat(sum(Alpha0,2), 1,nScores);       
                
    for n=1:size(p_casej,3)
        pat = p_casej(:,:,n); %prob of alpha and total
        
        if pat==0
            break;
        end	
	
        alpha = p_totalc(:,:,n) + Alpha0;
        total = p_total(:,:,n) + Alpha0_total;
        
        Pi_c = alpha ./ total;
        lnPi_c = psi(alpha) - psi(total);
        Hci_t = Hci_t - sum(sum( pti .* pat .* Pi_c .* lnPi_c, 1));
    end
    IG = Hci - Hci_t;
end

