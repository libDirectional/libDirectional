% Gerhard Kurz, Florian Pfaff, Uwe D. Hanebeck,
% Discretization of SO(3) Using Recursive Tesseract Subdivision
% Proceedings of the 2017 IEEE International Conference on Multisensor Fusion and Integration for Intelligent Systems (MFI 2017)
% Daegu, Korea, November 2017.

function result = tesseractsubdivision(n, normalize)
    if nargin < 2
        normalize = true;
    end
    
    %numer of points is f(2^m) where m=0,1, ... and f(n)=(n+1).^4-(n-1).^4
    %f is https://oeis.org/A008511
    cubes = generateTesseract();

    while length(unique(cell2mat(cubes),'rows')) < n
        newCubes = cell(8*length(cubes),1);
        for i=1:length(cubes)
            newCubes(8*i-7:8*i,1) = subdividecube(cubes{i});
        end
        cubes = newCubes;
    end
    
    result = unique(cell2mat(cubes),'rows');
    
    if normalize
        result = result./repmat(sqrt(sum(result.^2,2)), 1, 4);
    end
end

function cubes = generateTesseract()
    cube3d = (dec2bin(0:7)-'0')*2-1;
    
    cubes = cell(8,1);
    for i=1:4
        %insert column with +1 or -1 at location i
        cubes{2*i-1} = [cube3d(:,1:i-1), ones(8,1), cube3d(:,i:end)];
        cubes{2*i} = [cube3d(:,1:i-1), -ones(8,1), cube3d(:,i:end)];
    end
end

function cubes = subdividecube(cube)
    %divides given cube into 8 cubes
    assert(all(size(cube) == [8 4]));
    center = mean(cube);
    cubes = cell(8,1);
    for i=1:8 %create 8 new cubes
        p = cube(i,:); 
        %create cube with corners p and center        
        newCube = generateCube(p,center);
        cubes {i} = newCube;
    end
end

function newCube = generateCube(a,b)
    %generate axis aligned cube with opposite corners a and b
    x = zeros(16,4); %corner candidates
    %try all combinations of min/max in each dimension
    %todo optimize performance
    for i=1:16
        for j=1:4
            if bitget(i,j) == 1
                x(i,j) = min(a(j),b(j));
            else
                x(i,j) = max(a(j),b(j));
            end
        end
    end
    %throw away duplicates
    newCube = unique(x,'rows');
end