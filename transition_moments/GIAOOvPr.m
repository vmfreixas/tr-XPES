function S = GIAOOvPr(a,b,A,B,alpha,beta,ka,kb)
    LA = a(1)+a(2)+a(3);
    LB = b(1)+b(2)+b(3);

    xi = alpha*beta/(alpha+beta);
    zeta = alpha+beta;
    P = (alpha*A+beta*B)/(alpha+beta);
    S = 0;
    
    if (LA+LB==0)
        S = ((pi/zeta)^(3/2))*exp(-xi*dot(A-B,A-B));
        S = S*exp(-dot(ka+kb,ka+kb)/(4*zeta));
        S = S*exp(j*(dot(ka,P-A)+dot(kb,P-B)));
        return;
        
%     elseif (LA<LB)
%         if b(1)>0
%             ii = 1;
%         elseif b(2)>0
%             ii = 2;
%         elseif b(3)
%             ii = 3;
%         end
%         S = GIAOOvPr(pli(a,ii),mi(b,ii),A,B,alpha,beta,ka,kb) ...
%             +(A(ii)-B(ii))*GIAOOvPr(a,mi(b,ii),A,B,alpha,beta,ka,kb);
%         
    elseif ((LA>0)&&(LB==0))
        if a(1)>0
            ii = 1;
        elseif a(2)>0
            ii = 2;
        elseif a(3)>0
            ii = 3;
        end
        S = (P(ii)-A(ii)+j*(ka(ii)+kb(ii))/(2*zeta))*GIAOOvPr(mi(a,ii),b,A,B,alpha,beta,ka,kb);
        if a(ii)>1
            S = S+1/(2*zeta)*(a(ii)-1)*GIAOOvPr(mi(mi(a,ii),ii),b,A,B,alpha,beta,ka,kb);
        end
    elseif ((LB>0)&&(LA==0))
        if b(1)>0
            ii = 1;
        elseif b(2)>0
            ii = 2;
        elseif b(3)>0
            ii = 3;
        end
        S = (P(ii)-B(ii)+j*(ka(ii)+kb(ii))/(2*zeta))*GIAOOvPr(a,mi(b,ii),A,B,alpha,beta,ka,kb);
        if b(ii)>1
            S = S+1/(2*zeta)*(b(ii)-1)*GIAOOvPr(a,mi(mi(b,ii),ii),A,B,alpha,beta,ka,kb);
        end
    elseif ((LA>0)&&(LB>0))
        if b(1)>0
            ii = 1;
        elseif b(2)>0
            ii = 2;
        elseif b(3)>0
            ii = 3;
        end
%         vertical recursion
%         S = (P(ii)-B(ii)+j*(ka(ii)+kb(ii))/(2*zeta))*GIAOOvPr(a,mi(b,ii),A,B,alpha,beta,ka,kb);
%         if b(ii)>1
%             S = S + 1/(2*zeta)*(b(ii)-1)*GIAOOvPr(a,mi(mi(b,ii),ii),A,B,alpha,beta,ka,kb);
%         end 
%         if a(ii)>0
%             S = S + 1/(2*zeta)*a(ii)*GIAOOvPr(mi(a,ii),mi(b,ii),A,B,alpha,beta,ka,kb);
%         end
        
        S = GIAOOvPr(pli(a,ii),mi(b,ii),A,B,alpha,beta,ka,kb) ...
            +(A(ii)-B(ii))*GIAOOvPr(a,mi(b,ii),A,B,alpha,beta,ka,kb);
        
    end
%     return;
    
        