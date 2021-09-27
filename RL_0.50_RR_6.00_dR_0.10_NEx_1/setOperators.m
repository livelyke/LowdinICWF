
fprintf('Setting Operators\n')
valsDiag = [-2*ones(nDim,1)];
vals = [ones(nDim,1)];
nLapl = spdiags(valsDiag,0,nDim,nDim) ...
            + spdiags(vals, 1, nDim,nDim)     + spdiags(vals, -1, nDim,nDim);

nLapl = (1/dxn^2)*nLapl; 
Tn = (-1/(2*mu_n))*nLapl;

Vnn = 1./nAxis;

[chi, dum] = eigs(nLapl,NnucOrb);
for nu=1:NnucOrb
    chi(:,nu) = chi(:,nu)./sqrt(chi(:,nu)'*chi(:,nu)*dxn);
end

%% Kernel matrices

KTn = chi'*Tn*chi*dxn; KTn(abs(KTn)<1e-13) = 0;
KVnn = chi'*(repmat(Vnn,1,NnucOrb).*chi)*dxn; 
In = eye(NnucOrb);

% Read the file as cell string line by line:
file = 'Te_Vee/Te';
fid = fopen(file, 'r');
if fid < 0, error('Cannot open file: %s', file); end
DataTe = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
fclose(fid);

file = strcat('Te_Vee/Vee');
% Read the file as cell string line by line:
fid = fopen(file, 'r');
if fid < 0, error('Cannot open file: %s', file); end
DataVee = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
fclose(fid);

KTe = sparse(NeleSpatialOrb,NeleSpatialOrb);
for i=1:size(DataTe{1}) %Read in excess for now, fix later to only read in what's needed
        a = strsplit(DataTe{1}{i},' ');
        i1 = a(2); i1 = str2num(i1{1});
        i2 = a(3); i2 = str2num(i2{1});
        MTeval = a(end); MTeval = strsplit(MTeval{1},',');
        MTeval = MTeval(1); MTeval = MTeval{1}; MTeval = str2double(MTeval(2:end));
        KTe(i1,i2) = MTeval;
end


%% KVee
% Collective index n1n2, n3n4 as
%   n1n2 = n1 + NeleSpatialOrb*(n2-1);
KVee = zeros(NeleSpatialOrb^2,NeleSpatialOrb^2);
for i=2:size(DataVee{1})
        a = strsplit(DataVee{1}{i},' ');
        n1 = a(2);    n1 = str2num(n1{1});
        n2 = a(4);    n2 = str2num(n2{1});
        n3 = a(6);    n3 = str2num(n3{1});
        n4 = a(8);    n4 = str2num(n4{1});
        MTeval = a(10); MTeval = str2double(MTeval{1});
        n1n2 = n1 + NeleSpatialOrb*(n2-1);
        n3n4 = n3 + NeleSpatialOrb*(n4-1);
        KVee(n1n2,n3n4) = MTeval;
end


%% Assemble Hamiltonian, H_{IJ,I'J'} flattened with nuclear index running
%fastest
%i.e. index IJ = J + N_J*(I-1)
%Taking only the upper diagonal part

NIJ = NSlater*NnucOrb;

%Could be made faster by getting nonzero indices in advance and
%   filling with dummy values
MTn = sparse(NIJ,NIJ);
MTe = sparse(NIJ,NIJ);
MVee = sparse(NIJ,NIJ);
MVen = sparse(NIJ,NIJ);
MVnn = sparse(NIJ,NIJ);


for I=1:NSlater
    for Ip=I:NSlater
        IJ = (1:NnucOrb) + NnucOrb*(I-1);
        IpJp = (1:NnucOrb) + NnucOrb*(Ip-1);

        %MTn
        if(Ip==I)
            MTn(IJ,IpJp) = triu(KTn);
        end
        
        if(CITruncation=='doubles')
            %In this case there are no slater determinants which vary by
            %only one spin orbital
            if(Ip==I)
                %Get spin orbital indeces associated to I. 
                spinInd =    slaterIndices(I).spinOrbitals;
                spaceInd = slaterIndices(I).spatialOrbitals;
                
                %MTe 
                MTeval=0;
                MVeeval =0;
                for i=1:Nele
                    MTeval = MTeval + KTe(spaceInd(i),spaceInd(i));
                    for j=(i+1):Nele
                        ii = spaceInd(i) + NeleSpatialOrb*(spaceInd(i)-1);
                        jj = spaceInd(j) + NeleSpatialOrb*(spaceInd(j)-1);
                        ij = spaceInd(i) + NeleSpatialOrb*(spaceInd(j)-1);
                        ji = spaceInd(j) + NeleSpatialOrb*(spaceInd(i)-1);
                        MVeeval = MVeeval  +  KVee(ii,jj) - KVee(ij,ji);
                    end
                end
                MTe(IJ,IpJp) = MTeval*In;
                
                MVee(IJ,IpJp) = MVeeval*In;
                
            end
        end
        
        
    end
end
            
        

