
%Slater Determinant class
%Has the spin Orbital indices
%TODO: Can perform maximal alignment proecedure
classdef SlaterIndex
        properties
                spinOrbitals
                spatialOrbitals
        end
        methods
            function obj = SlaterIndex(spinOrbitals)
                obj.spinOrbitals = spinOrbitals;
                spatialOrbitals = [];
                for i=1:max(size(spinOrbitals))
                    spatialOrbitals = [spatialOrbitals ceil(spinOrbitals(i)/2) ];
                end
                obj.spatialOrbitals = spatialOrbitals;
            end
        end
end
