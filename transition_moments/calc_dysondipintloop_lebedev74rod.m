%scho%calculate dyson integral by contracting MO coeff 
clear 

for igridpoint = 1:74
    stringofnum=string(igridpoint);
    
    file1 = "dysonorbint/total74rodMTSOdysonxreal"+stringofnum+".txt"
    file2 = "dysonorbint/total74rodMTSOdysonximag"+stringofnum+".txt"
    
    intRe = load(file1);
    intIm = load(file2);
    
    %dysonxint = load("dysonorbint/MTSOdysonxreal"+stringofnum+".txt")+i*load("dysonorbint/MTSOdysonximag"+stringofnum+".txt");
    dysonxint = intRe+i*intIm;
    
    
    file1 = "dysonorbint/total74rodMTSOdysonyreal"+stringofnum+".txt"
    file2 = "dysonorbint/total74rodMTSOdysonyimag"+stringofnum+".txt"
    
    intRe = load(file1);
    intIm = load(file2);
    
    %dysonyint = load("dysonorbint/MTSOdysonyreal"+stringofnum+".txt")+i*load("dysonorbint/MTSOdysonyimag"+stringofnum+".txt");
    dysonyint = intRe+i*intIm;
    
    file1 = "dysonorbint/total74rodMTSOdysonzreal"+stringofnum+".txt"
    file2 = "dysonorbint/total74rodMTSOdysonzimag"+stringofnum+".txt"
    
    intRe = load(file1);
    intIm = load(file2);
    
    %dysonzint = load("dysonorbint/MTSOdysonzreal"+stringofnum+".txt")+i*load("dysonorbint/MTSOdysonzimag"+stringofnum+".txt");
    dysonzint = intRe+i*intIm;
    
    
    %%
    % read k vector value 
    evperhartree = 27.211386245988;
    % this part is MTSO
    wavevec=[0.0 0.0 0.0];
    % k from 164.8 eV to 167 eV
    icounter = 0;
    ecounter = 0;
    for Eelectron = 0.1:0.1:16.0 
        %momentum
        ecounter = ecounter+1;
        Elist(ecounter)=Eelectron;
        kvalue = sqrt(Eelectron*2);
        kvalueEh = kvalue/evperhartree;
        
        for theta = 0:pi/12:2*pi
            zcomp = kvalueEh*cos(theta);
            
            for phi = 0:2*pi/36:35/36*2*pi
                icounter = icounter + 1;
                xcomp = kvalueEh*sin(theta)*cos(phi);
                ycomp = kvalueEh*sin(theta)*sin(phi);
                
                klist(icounter,1)=xcomp;
                klist(icounter,2)=ycomp;
                klist(icounter,3)=zcomp;
            end
        end
    end
    
    % kvec = load('../../../pyrrole/dyson/brainhole/dysonorbint/MTSOLedgeElist.txt');
    numkvec = icounter;
    angle = "60";
    dis = "17";
    for ii = 1:1
    
        inputnamereal = "/Users/victormanuelfreixaslemus/Desktop/Projects/Photoelectron_spectroscopy/tr_XPES_code/tr-XPES/DO_AO_" + string(ii) + ".dat";
        
        realtemp = importdata(inputnamereal); 
       
        MOneu1cat2 = realtemp;
    
        alphatemp1 = transpose(MOneu1cat2);
        NBasis = length(alphatemp1);
        alphatemp = alphatemp1(1:NBasis);
       
        intneu1cat2(:,1) = transpose((alphatemp)*dysonxint);
        intneu1cat2(:,2) = transpose((alphatemp)*dysonyint);
        intneu1cat2(:,3) = transpose((alphatemp)*dysonzint);
        
        outputnamereal = "New_dipoles_0/dipolesRe" + stringofnum + "_" + string(ii) + ".dat";
        outputnameimag = "New_dipoles_0/dipolesIm" + stringofnum + "_" + string(ii) + ".dat";
    
        writematrix(real(intneu1cat2),outputnamereal);
        writematrix(imag(intneu1cat2),outputnameimag);
    end

end 
