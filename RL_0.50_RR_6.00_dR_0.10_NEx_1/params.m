
dxn = 0.1;

nBoxL = 0.5; nBoxR = 10; nAxis = (nBoxL:dxn:nBoxR).'; nDim = max(size(nAxis));

m1 = 1836; m2 = 1836; mu_n = m1*m2/(m1+m2);

NnucOrb = 2;

%% Determine number of unique slater determinants with the same spin
Nele = 2;
NeleSpatialOrb = 2; 
NeleSpinOrb = 2*NeleSpatialOrb;

%multiplicity
if(mod(Nele,2)==0)
    spinState='singlet';
end

%Number of total possible slater determinants
NPossibleSlater = nchoosek(NeleSpinOrb,Nele);

CITruncation = 'doubles';  %Where to truncate the possible combinations 
if CITruncation=='singles'
        NSlater = 1 + Nele*(NeleSpinOrb - Nele)/2; % for singlet restricted moleclue.
elseif CITruncation=='doubles'
        NSlater = 1 + (Nele/2 * (NeleSpinOrb-Nele)/2)^2; % for singlet restricted moleclue.
elseif CITruncation == 'singles-doubles'
        NSlater = 1 + Nele*(NeleSpinOrb - Nele)/2 + ...
                                  (Nele/2 * (NeleSpinOrb-Nele)/2)^2; % for singlet restricted moleclue.
end

spatialOrbIndex = 1:NeleSpatialOrb;
spinOrbitalIndex = 1:NeleSpinOrb;

%Construct SlaterIndex Objects
unExcited = 1:Nele;

slaterIndices = [SlaterIndex([1,2])]; %Unexcited state
if(CITruncation == 'doubles') % Only implemented for even number Ne, single inital state
    for s1=1:2:Nele
        set = unExcited;
        for s1_promotion=(Nele+1):2:NeleSpinOrb
            set(s1) = s1_promotion;
            for s2=2:2:Nele
                for s2_promotion=(Nele+2):2:NeleSpinOrb
                    set(s2) = s2_promotion;
                    slaterIndices = [slaterIndices SlaterIndex(set)];
                end
            end
        end
    end
end
            

%% Time parameters
dt = 0.01; dtImag = 0.1;


%% Save Parameters
savePath = strcat('dump_NeSp_',num2str(NeleSpatialOrb),...
                                    '_boxL_',num2str(nBoxL),'_boxR_',num2str(nBoxR),...
                                    '_dxn_',num2str(dxn));
mkdir(savePath);

