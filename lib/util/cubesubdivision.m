% Gerhard Kurz, Florian Pfaff, Uwe D. Hanebeck,
% Discretization of SO(3) Using Recursive Tesseract Subdivision
% Proceedings of the 2017 IEEE International Conference on Multisensor Fusion and Integration for Intelligent Systems (MFI 2017)
% Daegu, Korea, November 2017.

function result = cubesubdivision(n, normalize)
    if nargin < 2
        normalize = true;
    end
    
    %number of points is 6*4^n+2
    quads = generateCube();

    while length(unique(cell2mat(quads),'rows')) < n
        newQuads = cell(4*length(quads),1);
        for i=1:length(quads)
            newQuads(4*i-3:4*i,1) = subdividequad(quads{i});
        end
        quads = newQuads;
    end
    
    result = unique(cell2mat(quads),'rows');

    if normalize    
        result = result./repmat(sqrt(sum(result.^2,2)), 1, 3);
    end
end

function quads = generateCube()
    quad2d = (dec2bin(0:3)-'0')*2-1;
    
    quads = cell(6,1);
    for i=1:3
        %insert column with +1 or -1 at location i
        quads{2*i-1} = [quad2d(:,1:i-1), ones(4,1), quad2d(:,i:end)];
        quads{2*i} = [quad2d(:,1:i-1), -ones(4,1), quad2d(:,i:end)];
    end
end

function quads = subdividequad(quad)
    %divides given cube into 8 cubes
    assert(all(size(quad) == [4 3]));
    center = mean(quad);
    quads = cell(4,1);
    for i=1:4 %create 4 new quads
        p = quad(i,:); 
        %create cube with corners p and center        
        newCube = generateQuad(p,center);
        quads {i} = newCube;
    end
end

function newQuads = generateQuad(a,b)
    %generate axis aligned cube with opposite corners a and b
    x = zeros(8,3); %corner candidates
    %try all combinations of min/max in each dimension
    %todo optimize performance
    for i=1:8
        for j=1:3
            if bitget(i,j) == 1
                x(i,j) = min(a(j),b(j));
            else
                x(i,j) = max(a(j),b(j));
            end
        end
    end
    %throw away duplicates
    newQuads = unique(x,'rows');
end