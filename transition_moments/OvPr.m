function S = OvPr(a,b,A,B,alpha,beta)
    LA = a(1)+a(2)+a(3);
    LB = b(1)+b(2)+b(3);

    xi = alpha*beta/(alpha+beta);
    zeta = alpha+beta;
    P = (alpha*A+beta*B)/(alpha+beta);
    S = 0;
    
    if (LA+LB==0)
        S = ((pi/zeta)^(3/2))*exp(-xi*dot(A-B,A-B));
        return;
    elseif ((LA>0)&&(LB==0))
        if a(1)>0
            ii = 1;
        elseif a(2)>0
            ii = 2;
        elseif a(3)>0
            ii = 3;
        end
        S = (P(ii)-A(ii))*OvPr(mi(a,ii),b,A,B,alpha,beta);
        if a(ii)>1
            S = S+1/(2*zeta)*(a(ii)-1)*OvPr(mi(mi(a,ii),ii),b,A,B,alpha,beta);
        end
    elseif ((LB>0)&&(LA==0))
        if b(1)>0
            ii = 1;
        elseif b(2)>0
            ii = 2;
        elseif b(3)>0
            ii = 3;
        end
        S = (P(ii)-B(ii))*OvPr(a,mi(b,ii),A,B,alpha,beta);
        if b(ii)>1
            S = S+1/(2*zeta)*(b(ii)-1)*OvPr(a,mi(mi(b,ii),ii),A,B,alpha,beta);
        end
    elseif ((LA>0)&&(LB>0))
        if b(1)>0
            ii = 1;
        elseif b(2)>0
            ii = 2;
        elseif b(3)>0
            ii = 3;
        end
        S = (P(ii)-B(ii))*OvPr(a,mi(b,ii),A,B,alpha,beta);
        if b(ii)>1
            S = S + 1/(2*zeta)*(b(ii)-1)*OvPr(a,mi(mi(b,ii),ii),A,B,alpha,beta);
        end 
        if a(ii)>0
            S = S + 1/(2*zeta)*a(ii)*OvPr(mi(a,ii),mi(b,ii),A,B,alpha,beta);
        end
    end
    return;
    
        