function S = compute_dysonint(Shells,ShellPairs,usecart,wavevec,mu)

nShell = length(Shells)-1;

% %print
% nshellindyson=length(Shells)

nShellPairs = length(ShellPairs);
% BH = [0 0 0.1];

for nSP = 1:nShellPairs
    index_i = ShellPairs(nSP).iShell.index;
    index_j = ShellPairs(nSP).jShell.index;
%     temp = bottomup_overlap(Shells(index_i),Shells(index_j),ShellPairs(nSP));
%     temp = bottomup_overlap_os(Shells(index_i),Shells(index_j),ShellPairs(nSP));
    if (usecart)
            IS = Shells(index_i).shellcartstart;
            IE = Shells(index_i).shellcartstart+Shells(index_i).sizecart-1;
            JS = Shells(index_j).shellcartstart;
            JE = Shells(index_j).shellcartstart+Shells(index_j).sizecart-1;
            JS = 1;
            JE = 1;
    else
        IS = Shells(index_i).shellsphstart;
        IE = Shells(index_i).shellsphstart + 2*(Shells(index_i).l); 
        JS = Shells(index_j).shellsphstart; 
        JE = Shells(index_j).shellsphstart + 2*(Shells(index_j).l);
        JS = 1;
        JE = 1;
    end

    % start using recursion
    
    nContraction_i = length(Shells(index_i).zeta);
    nContraction_j = length(Shells(index_j).zeta);
    nContraction_j = 1;

    S_cart = [];  % overlap integral in cartesian gaussian 
    
    for car1 = 1:Shells(index_i).sizecart
        for car2 = 1:Shells(index_j).sizecart
            lA = Shells(index_i).cartesian(car1).xyz;
            lB = Shells(index_j).cartesian(car2).xyz;
            A = Shells(index_i).O;
            B = Shells(index_j).O;
            
%             ka = 1/2*cross(A,BH);
            ka = wavevec;
            onei = [0 0 0];
            onei(mu)=1;
            kb = [0.0 0.0 0.0];
            
            Selement = 0.0;
            for k =1:nContraction_i
                for l = 1:nContraction_j
                    
                    Stemp = GIAOOvPr(lA+onei,lB,A,B,Shells(index_i).zeta(k),Shells(index_j).zeta(l),ka,kb)...
                        +A(mu)*GIAOOvPr(lA,lB,A,B,Shells(index_i).zeta(k),Shells(index_j).zeta(l),ka,kb); % GIAO
                    
                    Selement = Selement+Stemp*Shells(index_i).normcoeff(k) ...
                        *Shells(index_j).normcoeff(l);
           
                end
            end
          %  S(IS+car1-1,JS+car2-1) = Selement;
            S_cart(car1,car2) = Selement;
        end
    end
    
    if ( usecart )
        S(IS:IE,JS:JE) = S_cart;
    else 
        S_sph = [];
        for isph = 1:(2*Shells(index_i).l+1) 
            for jsph = 1:(2*Shells(index_j).l+1)
                transtmp = 0;
                for icart = 1:Shells(index_i).sizecart
                    for jcart = 1:Shells(index_j).sizecart
                        transtmp = transtmp+Shells(index_i).car2sph_matrix(isph,icart) ...
                            *Shells(index_j).car2sph_matrix(jsph,jcart)*S_cart(icart,jcart);
                    end
                end
                S_sph(isph,jsph) = transtmp;
            end
        end
        %copy S_sph to S
%         size(S_sph)
        S(IS:IE,JS:JE) = S_sph;
    end
    
    % end using recursion
    
%     temp=0;
%     S(IS:IE,JS:JE)=temp;
end


