%% Debug

%% Check the order of kron to manual construction

NSlater = 2;
NnucOrb=2;

KTn = rand(NnucOrb); KTn = (KTn + KTn.')/2;
Se1 = eye(NSlater);

A1 = kron(Se1,KTn);

A2 = zeros(NSlater*NnucOrb);

for I=1:NSlater
    for Ip = I:NSlater
        if(I==Ip)
            IJ = (1:NnucOrb) + NnucOrb*(I-1);
            IJp = (1:NnucOrb) + NnucOrb*(Ip-1);
            A2(IJ,IJp) = triu(KTn);
        end
    end
end

