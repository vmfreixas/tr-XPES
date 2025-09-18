%scho%calculate dyson integral by contracting MO coeff 
%clear 

function calc_dysondipintloop_lebedev74rod(orbIntDir, DOFileName, nDO, outputDir, e0, ef, de, r)

for igridpoint = 1:74
    stringofr=string(r);
    stringofnum=string(igridpoint);
    
    file1 = orbIntDir + "/total74rodMTSOdysonxreal"+stringofnum+"_"+stringofr+".txt"
    file2 = orbIntDir + "/total74rodMTSOdysonximag"+stringofnum+"_"+stringofr+".txt"
    
    intRe = load(file1);
    intIm = load(file2);
    
    %dysonxint = load("dysonorbint/MTSOdysonxreal"+stringofnum+".txt")+i*load("dysonorbint/MTSOdysonximag"+stringofnum+".txt");
    dysonxint = intRe+i*intIm;
    
    
    file1 = orbIntDir + "/total74rodMTSOdysonyreal"+stringofnum+"_"+stringofr+".txt"
    file2 = orbIntDir + "/total74rodMTSOdysonyimag"+stringofnum+"_"+stringofr+".txt"
    
    intRe = load(file1);
    intIm = load(file2);
    
    %dysonyint = load("dysonorbint/MTSOdysonyreal"+stringofnum+".txt")+i*load("dysonorbint/MTSOdysonyimag"+stringofnum+".txt");
    dysonyint = intRe+i*intIm;
    
    file1 = orbIntDir + "/total74rodMTSOdysonzreal"+stringofnum+"_"+stringofr+".txt"
    file2 = orbIntDir + "/total74rodMTSOdysonzimag"+stringofnum+"_"+stringofr+".txt"
    
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
    for Eelectron = e0:de:ef 
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
    for ii = 1:nDO
    
        inputnamereal = DOFileName + string(ii) + ".dat"
        
        realtemp = importdata(inputnamereal); 
       
        MOneu1cat2 = realtemp;
    
        alphatemp1 = transpose(MOneu1cat2);
        NBasis = length(alphatemp1);
        alphatemp = alphatemp1(1:NBasis);
       
        intneu1cat2(:,1) = transpose((alphatemp)*dysonxint);
        intneu1cat2(:,2) = transpose((alphatemp)*dysonyint);
        intneu1cat2(:,3) = transpose((alphatemp)*dysonzint);
        
        outputnamereal = outputDir + "/dipolesRe" + stringofnum + "_" + string(ii) + "_" + stringofr + ".dat";
        outputnameimag = outputDir + "/dipolesIm" + stringofnum + "_" + string(ii) + "_" + stringofr + ".dat";
    
        writematrix(real(intneu1cat2),outputnamereal);
        writematrix(imag(intneu1cat2),outputnameimag);
    end

end

end
