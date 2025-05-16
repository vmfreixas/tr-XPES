function  out = calcdysonint_up( xyz,Atoms,TotalCharge,Options,integralDir,e0,ef,de,r )

Set = Options.BasisSet;  % name of basis set, e.g. sto-3g

Basis = basisread(Set);

NA = length(Atoms);   % number of atoms
nShell = 0;           % total number of shells
nShell_atom = 0;      % number of shells in an atom, SP shell count as one shell here
Shells = ([]);        % Shells store all the shells

% generate shells start
for k = 1:NA   % loop over atoms 
    nShell_atom = length(Basis{Atoms(k)});
    for l = 1:nShell_atom  % loop over shells in an atom
        if Basis{Atoms(k)}(l).shelltype == 'S'
            nShell = nShell+1 ;
            Shells(nShell).l = 0;  % angular momentum of the shell
            Shells(nShell).O = xyz(k,:);   % origin of the shell
            Shells(nShell).zeta = Basis{Atoms(k)}(l).exponents;
            Shells(nShell).coeffs = Basis{Atoms(k)}(l).coeffs;
            ncontra = length(Shells(nShell).coeffs);
            for nm = 1:ncontra
                Shells(nShell).norm(nm) = (2/pi)^(3/4)*(Shells(nShell).zeta(nm))^(3/4);
            end
        
        elseif Basis{Atoms(k)}(l).shelltype == 'SP'
            nShell = nShell+1 ;
            Shells(nShell).l = 0;  % angular momentum of the shell
            Shells(nShell).O = xyz(k,:);   % origin of the shell
            Shells(nShell).zeta = Basis{Atoms(k)}(l).exponents;
            Shells(nShell).coeffs = Basis{Atoms(k)}(l).coeffs(1,:);
            ncontra = length(Shells(nShell).coeffs);
            for nm = 1:ncontra
                Shells(nShell).norm(nm) = (2/pi)^(3/4)*(Shells(nShell).zeta(nm))^(3/4);
            end
            % so far add the S shell
            
            nShell = nShell+1 ;
            Shells(nShell).l = 1;  % angular momentum of the shell
            Shells(nShell).O = xyz(k,:);   % origin of the shell
            Shells(nShell).zeta = Basis{Atoms(k)}(l).exponents;
            Shells(nShell).coeffs = Basis{Atoms(k)}(l).coeffs(2,:);
            ncontra = length(Shells(nShell).coeffs);
            for nm = 1:ncontra
                Shells(nShell).norm(nm) = (2/pi)^(3/4)*2*(Shells(nShell).zeta(nm))^(5/4);
            end
            % add the P shell
        
        elseif Basis{Atoms(k)}(l).shelltype == 'P'
            nShell = nShell+1 ;
            Shells(nShell).l = 1;  % angular momentum of the shell
            Shells(nShell).O = xyz(k,:);   % origin of the shell
            Shells(nShell).zeta = Basis{Atoms(k)}(l).exponents;
            Shells(nShell).coeffs = Basis{Atoms(k)}(l).coeffs;
            ncontra = length(Shells(nShell).coeffs);
            for nm = 1:ncontra
                Shells(nShell).norm(nm) = (2/pi)^(3/4)*2*(Shells(nShell).zeta(nm))^(5/4);
            end
            
        elseif Basis{Atoms(k)}(l).shelltype == 'D'
            nShell = nShell+1 ;
            Shells(nShell).l = 2;  % angular momentum of the shell
            Shells(nShell).O = xyz(k,:);   % origin of the shell
            Shells(nShell).zeta = Basis{Atoms(k)}(l).exponents;
            Shells(nShell).coeffs = Basis{Atoms(k)}(l).coeffs;
            ncontra = length(Shells(nShell).coeffs);
            for nm = 1:ncontra
                Shells(nShell).norm(nm) = (2/pi)^(3/4)*(2^2)*(Shells(nShell).zeta(nm))^(7/4)/sqrt(3);
            end
            
        elseif Basis{Atoms(k)}(l).shelltype == 'F'
            nShell = nShell+1 ;
            Shells(nShell).l = 3;  % angular momentum of the shell
            Shells(nShell).O = xyz(k,:);   % origin of the shell
            Shells(nShell).zeta = Basis{Atoms(k)}(l).exponents;
            Shells(nShell).coeffs = Basis{Atoms(k)}(l).coeffs;
            ncontra = length(Shells(nShell).coeffs);
            for nm = 1:ncontra
                Shells(nShell).norm(nm) = (2/pi)^(3/4)*(2^3)*(Shells(nShell).zeta(nm))^(9/4)/sqrt(5*3);
            end
            
        elseif Basis{Atoms(k)}(l).shelltype == 'G'
            nShell = nShell+1 ;
            Shells(nShell).l = 4;  % angular momentum of the shell
            Shells(nShell).O = xyz(k,:);   % origin of the shell
            Shells(nShell).zeta = Basis{Atoms(k)}(l).exponents;
            Shells(nShell).coeffs = Basis{Atoms(k)}(l).coeffs;
            ncontra = length(Shells(nShell).coeffs);
            for nm = 1:ncontra
                Shells(nShell).norm(nm) = (2/pi)^(3/4)*(2^4)*(Shells(nShell).zeta(nm))^(11/4)/sqrt(7*5*3);
            end
            
        elseif Basis{Atoms(k)}(l).shelltype == 'H'
            nShell = nShell+1 ;
            Shells(nShell).l = 5;  % angular momentum of the shell
            Shells(nShell).O = xyz(k,:);   % origin of the shell
            Shells(nShell).zeta = Basis{Atoms(k)}(l).exponents;
            Shells(nShell).coeffs = Basis{Atoms(k)}(l).coeffs;
            ncontra = length(Shells(nShell).coeffs);
            for nm = 1:ncontra
                Shells(nShell).norm(nm) = (2/pi)^(3/4)*(2^5)*(Shells(nShell).zeta(nm))^(13/4)/sqrt(9*7*5*3);
            end
            
        else
            1000000000000000
        end

    end
end

% % print!!!
% trueshells=size(Shells)

% extra shell is 
Shells(nShell+1).l=0;
Shells(nShell+1).O=[0.0 0.0 0.0];
Shells(nShell+1).zeta=0.0;
Shells(nShell+1).coeffs=1.0;
Shells(nShell+1).norm=1.0;

% generate shells end

% start of the shell
cartsize_so_far = 0;
sphsize_so_far = 0;
for k = 1:nShell+1
    Shells(k).shellcartstart = cartsize_so_far+1;
    Shells(k).shellsphstart = sphsize_so_far+1;
    langular = Shells(k).l;
    Shells(k).sizecart = (langular+1)*(langular+2)/2;
    cartsize_so_far = cartsize_so_far + Shells(k).sizecart;
    sphsize_so_far = sphsize_so_far + langular*2+1;
    Shells(k).normcoeff = (Shells(k).norm).*Shells(k).coeffs;
            % list of cartesian angular momentum index
            fuck = 0;
            for xx = 0:Shells(k).l
                x = Shells(k).l-xx;
                for yy = 0:xx
                    y = xx-yy;
                    z = Shells(k).l-x-y;
                    fuck = fuck +1;
                    Shells(k).cartesian(fuck).xyz=[x,y,z];
                end
            end
    % here generate cart to sph transform coeff
    carsize =  (Shells(k).l+1)*(Shells(k).l+2)/2;
    sphsize =  2*Shells(k).l+1;
    % Shells(k).car2sph_matrix(sph,car) first index is spherical
    % index, second is cartesian index
    Shells(k).car2sph_matrix = zeros(sphsize,carsize);
    if (Shells(k).l == 0)
        Shells(k).car2sph_matrix(1,1)=1;
    elseif ( Shells(k).l == 1 )
        Shells(k).car2sph_matrix(1,1)=1;
        Shells(k).car2sph_matrix(2,2)=1;
        Shells(k).car2sph_matrix(3,3)=1;
    else
        for p = 0:Shells(k).l
            for q = 1:carsize
                ltemp = Shells(k).cartesian(q).xyz;
                m = p-Shells(k).l;
                scalenorm = sqrt( doubfac(Shells(k).l)/(doubfac(ltemp(1))*doubfac(ltemp(2))*doubfac(ltemp(3))) );
                if ( m<0 )
                    cplxcoeff = car2sph(Shells(k).l,m,ltemp(1),ltemp(2),ltemp(3));
                    Shells(k).car2sph_matrix(p+1,q)=scalenorm*sqrt(2)*(-imag(cplxcoeff));
                    Shells(k).car2sph_matrix(2*Shells(k).l+1-p,q)=scalenorm*sqrt(2)*real(cplxcoeff);
                elseif (m==0)
                    Shells(k).car2sph_matrix(p+1,q)=scalenorm*real(car2sph(Shells(k).l,m,ltemp(1),ltemp(2),ltemp(3)));
                end
            end
        end
    end
            
end
    
% generate 2e specific shell pairs start
nShellpairs = nShell * (nShell+1)/2; % number of shell pairs
nS = 0;   % loop over shellpairs
for k = 1:nShell
    l=nShell+1;
    nS = nS+1;
    % make sure angular momentum of bra is greater or equal to ket
        
    ShellPairs2e(nS).iShell.index = k;
    ShellPairs2e(nS).jShell.index = l;
    index_i = k;
    index_j = l;
            
    A = Shells(index_i).O;
    B = Shells(index_j).O;
    ShellPairs2e(nS).A = A;
    ShellPairs2e(nS).B = B;
    ShellPairs2e(nS).AB = A-B;
        
        nContraction_i = length(Shells(index_i).zeta);
        nContraction_j = length(Shells(index_j).zeta);
        nPGTOPair = 0; % number of primitive gaussian type orbital pair
        for m = 1:nContraction_i
            for n = 1:nContraction_j
                nPGTOPair = nPGTOPair+1;
                ShellPairs2e(nS).idx_a(nPGTOPair)= m;
                ShellPairs2e(nS).idx_b(nPGTOPair)= n;
                ShellPairs2e(nS).Zeta(nPGTOPair) = Shells(index_i).zeta(m)+Shells(index_j).zeta(n);
                ShellPairs2e(nS).xi(nPGTOPair)   = Shells(index_i).zeta(m)*Shells(index_j).zeta(n)/ShellPairs2e(nS).Zeta(nPGTOPair);
                ShellPairs2e(nS).P(nPGTOPair,:)  = (Shells(index_i).zeta(m)*A+Shells(index_j).zeta(n)*B)/ShellPairs2e(nS).Zeta(nPGTOPair);
                ShellPairs2e(nS).K(nPGTOPair)    = (pi/ShellPairs2e(nS).Zeta(nPGTOPair))^(3/2)*exp(-ShellPairs2e(nS).xi(nPGTOPair)*dot(A-B,A-B));
                ShellPairs2e(nS).Klib(nPGTOPair) = 1/ShellPairs2e(nS).Zeta(nPGTOPair)*exp(-ShellPairs2e(nS).xi(nPGTOPair)*dot(A-B,A-B));
                ShellPairs2e(nS).tmpval(nPGTOPair) = (pi^(3/2))*ShellPairs2e(nS).Klib(nPGTOPair)/sqrt(ShellPairs2e(nS).Zeta(nPGTOPair)) ...
                    *Shells(index_i).coeffs(m)*Shells(index_j).coeffs(n) ...
                    *Shells(index_i).norm(m)*Shells(index_j).norm(n);
                ShellPairs2e(nS).K(nPGTOPair) = ShellPairs2e(nS).K(nPGTOPair)*Shells(index_i).coeffs(m)*Shells(index_j).coeffs(n);
                ShellPairs2e(nS).K(nPGTOPair) = ShellPairs2e(nS).K(nPGTOPair)*Shells(index_i).norm(m)*Shells(index_j).norm(n);
            end
        end
        ShellPairs2e(nS).nPGTOPair = nPGTOPair;
   
end

% generate 2e specific shell pairs end


nS = 0;
% generate shellpairs
for k = 1:nShell
    for l = k:nShell
        nS = nS+1;
%         if (Shells(k).l>=Shells(l).l)
%             ShellPairs(nS).iShell.index = k;
%             ShellPairs(nS).jShell.index = l;
%             index_i = k;
%             index_j = l;
%         else
%             ShellPairs(nS).iShell.index = l;
%             ShellPairs(nS).jShell.index = k;
%             index_i = l;
%             index_j = k;
%         end
        % make sure angular momentum of bra is greater or equal to ket
        
            ShellPairs(nS).iShell.index = k;
            ShellPairs(nS).jShell.index = l;
            index_i = k;
            index_j = l;
            
        A = Shells(index_i).O;
        B = Shells(index_j).O;
        ShellPairs(nS).A = A;
        ShellPairs(nS).B = B;
        ShellPairs(nS).AB = A-B;
        
        nContraction_i = length(Shells(index_i).zeta);
        nContraction_j = length(Shells(index_j).zeta);
        nPGTOPair = 0; % number of primitive gaussian type orbital pair
        for m = 1:nContraction_i
            for n = 1:nContraction_j
                nPGTOPair = nPGTOPair+1;
%                 ShellPairs(nS).idx_a(nPGTOPair)= m;
%                 ShellPairs(nS).idx_b(nPGTOPair)= n;
                ShellPairs(nS).Zeta(nPGTOPair) = Shells(index_i).zeta(m)+Shells(index_j).zeta(n);
                ShellPairs(nS).xi(nPGTOPair)   = Shells(index_i).zeta(m)*Shells(index_j).zeta(n)/ShellPairs(nS).Zeta(nPGTOPair);
                ShellPairs(nS).P(nPGTOPair,:)  = (Shells(index_i).zeta(m)*A+Shells(index_j).zeta(n)*B)/ShellPairs(nS).Zeta(nPGTOPair);
                ShellPairs(nS).K(nPGTOPair)    = (pi/ShellPairs(nS).Zeta(nPGTOPair))^(3/2)*exp(-ShellPairs(nS).xi(nPGTOPair)*dot(A-B,A-B));
                ShellPairs(nS).Klib(nPGTOPair) = 1/ShellPairs(nS).Zeta(nPGTOPair)*exp(-ShellPairs(nS).xi(nPGTOPair)*dot(A-B,A-B));
                ShellPairs(nS).tmpval(nPGTOPair) = (pi^(3/2))*ShellPairs(nS).Klib(nPGTOPair)/sqrt(ShellPairs(nS).Zeta(nPGTOPair)) ...
                    *Shells(index_i).coeffs(m)*Shells(index_j).coeffs(n) ...
                    *Shells(index_i).norm(m)*Shells(index_j).norm(n);
                ShellPairs(nS).K(nPGTOPair) = ShellPairs(nS).K(nPGTOPair)*Shells(index_i).coeffs(m)*Shells(index_j).coeffs(n);
                ShellPairs(nS).K(nPGTOPair) = ShellPairs(nS).K(nPGTOPair)*Shells(index_i).norm(m)*Shells(index_j).norm(n);
            end
        end
        ShellPairs(nS).nPGTOPair = nPGTOPair;
    end
end
% generate shellpaie end
% % print!!!
% shellpairsize = size(ShellPairs)

% whether or not using cartesian gaussian
usecart = false; 
% useGIAO = true;
useGIAO = false;

% define magnetic field BH

% Bint = [0.0 0.0 0.02];
Bint = [0.0 0.0 0.0];

% bottom up calculate overlap S

S = compute_overlapS(Shells(1:nShell),ShellPairs,usecart,useGIAO,Bint);   % unnormalized upper triangle
[nbf1,nbf2]= size(S);

if (useGIAO) 
    for k = 2:nbf1
       for l = 1:k-1
          S(k,l) = conj(S(l,k));
          
       end
    end
else
    for k = 2:nbf1
      for l = 1:k-1
         S(k,l) = S(l,k);
         
      end
    end
end


for k = 1:nbf1
    NormConst(k)=1/sqrt(S(k,k));
end

% NormConst

% this part is pyrrole molecule 
% wavevec=[0.0 0.0 0.0];
% % load('k.mat');
% 
% 
% part1=load('k_z.mat');
% part2=load('k2.mat');
% 
% kvec=[part1.k,part2.k];
% 
% writematrix(kvec','dysonorbint/kvectorlist.txt');
% % kvec
% cosk = cos(pi/4)*kvec;
% kvec = cosk;
% 
% mu = 1; % x component of dipole
% for ii = 1:length(kvec)
%     wavevec(1)=kvec(ii);
%     wavevec(2)=kvec(ii);
%     dysoninttemp = compute_dysonint(Shells,ShellPairs2e,usecart,wavevec,mu);
%     dysonint(:,ii) = (NormConst').*dysoninttemp;
% end
% 
% mu = 2; % y component of dipole
% for ii = 1:length(kvec)
%     wavevec(1)=kvec(ii);
%     wavevec(2)=kvec(ii);
%     dysoninttemp = compute_dysonint(Shells,ShellPairs2e,usecart,wavevec,mu);
%     dysoninty(:,ii) = (NormConst').*dysoninttemp;
% end
% 
% mu = 3; % z component of dipole
% for ii = 1:length(kvec)
%     wavevec(1)=kvec(ii);
%     wavevec(2)=kvec(ii);
%     dysoninttemp = compute_dysonint(Shells,ShellPairs2e,usecart,wavevec,mu);
%     dysonintz(:,ii) = (NormConst').*dysoninttemp;
% end

% read lebedev grid
lebedevthetaphilist= readmatrix('lebedevthetaphilist74.txt');

evperhartree = 27.211386245988;
% this part is MTSO
wavevec=[0.0 0.0 0.0];
% k from 164.8 eV to 167 eV


for igrid = 1:74
    
    igrid
    
    lebtheta = lebedevthetaphilist(igrid,1);
    lebphi   = lebedevthetaphilist(igrid,2);
    
    xyzofgrid = [sin(lebtheta)*cos(lebphi),sin(lebtheta)*sin(lebphi),cos(lebtheta)];
    
    ecounter = 0;
    icounter = 0;
    for Eelectron = e0:de:ef  % (eV) 
        
        icounter = icounter+1;
        %momentum
        ecounter = ecounter+1;
        Elist(ecounter)=Eelectron;
        kvalue = sqrt(Eelectron/evperhartree*2);
        kvalueEh = kvalue;
        
        krotated(1) = xyzofgrid(1)*kvalueEh;
        krotated(2) = xyzofgrid(2)*kvalueEh;
        krotated(3) = xyzofgrid(3)*kvalueEh;

        klist(icounter,1)=krotated(1);
        klist(icounter,2)=krotated(2);
        klist(icounter,3)=krotated(3);
    end
    
    % compute integral 
    mu = 1; % x component of dipole
    for ii = 1:length(klist)
    %     ii
        kofuse=klist(ii,:);
        dysoninttemp = compute_dysonint(Shells,ShellPairs2e,usecart,kofuse,mu);
        dysonintx(:,ii) = (NormConst').*dysoninttemp;
    end

    mu = 2; % y component of dipole
    for ii = 1:length(klist)
        kofuse=klist(ii,:);
        dysoninttemp = compute_dysonint(Shells,ShellPairs2e,usecart,kofuse,mu);
        dysoninty(:,ii) = (NormConst').*dysoninttemp;
    end

    mu = 3; % z component of dipole
    for ii = 1:length(klist)
        kofuse=klist(ii,:);
        dysoninttemp = compute_dysonint(Shells,ShellPairs2e,usecart,kofuse,mu);
        dysonintz(:,ii) = (NormConst').*dysoninttemp;
    end
    
    clear klist;
    % write matrix
    writematrix(real(dysonintx),integralDir+"/total74rodMTSOdysonxreal"+string(igrid)+"_"+string(r)+".txt");
    writematrix(imag(dysonintx),integralDir+"/total74rodMTSOdysonximag"+string(igrid)+"_"+string(r)+".txt");

    writematrix(real(dysoninty),integralDir+"/total74rodMTSOdysonyreal"+string(igrid)+"_"+string(r)+".txt");
    writematrix(imag(dysoninty),integralDir+"/total74rodMTSOdysonyimag"+string(igrid)+"_"+string(r)+".txt");
    
    writematrix(real(dysonintz),integralDir+"/total74rodMTSOdysonzreal"+string(igrid)+"_"+string(r)+".txt");
    writematrix(imag(dysonintz),integralDir+"/total74rodMTSOdysonzimag"+string(igrid)+"_"+string(r)+".txt");
    clear dysontx;
    clear dysonty;
    clear dysontz;
    
end
%writematrix(Elist,'MTSOLedgeElistnew.txt');





out.shells = Shells;
out.shellpairs = ShellPairs;
out.S = S;
out.norm=NormConst;
out.dysonx=dysonintx;
out.dysony=dysoninty;
out.dysonz=dysonintz;

end


