%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLinMap.m                                                               %
%    Turn piece-wise constant realization into piece-wise linear.         %
%    realization                                                          %
% Author: Geir-Arne Fuglstad <geirarne.fuglstad@gmail.com> (2019)         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Slin] = PLinMap(tv, vLoc, S, loc)
    % Create coordinates of centroids of triangles
    cVec = zeros(size(tv,1), 3);
    for i = 1:size(tv,1)
        for j = 1:3
            cVec(i,:) = cVec(i,:) + vLoc(tv(i,j),:);
        end
    end
    cVec = cVec/3;
    
    % Create vt map
    tmpI = zeros(prod(size(tv)), 1);
    tmpJ = tmpI;
    cnt = 1;
    for i = 1:size(tv,1)
        for j = 1:size(tv,2)
            tmpI(cnt) = tv(i,j);
            tmpJ(cnt) = i;
            cnt = cnt + 1;
        end
    end
    
    % Sort
    [vIdx, Idx] = sort(tmpI);
    tIdx = tmpJ(Idx);
    
    % Make indexes to look up triangles for each vertex
    [a,~] = histc(vIdx, 1:size(vLoc,1));
    a = [1;a];
    sIdx = cumsum(a);
    eIdx = sIdx(2:end)-1;
    sIdx(end) = [];
    
    % Calculate weights on each vertex
    vWs = cell(size(vLoc,1),1);
    for i = 1:size(vLoc,1)
        ws = zeros(eIdx(i)-sIdx(i)+1, 1);
        for j = 1:size(ws,1)
           ws(j) = norm(cVec(tIdx(sIdx(i)+j-1),:)-vLoc(i,:),2);
        end
        ws = ws/sum(ws);
        vWs{i} = [tIdx(sIdx(i):eIdx(i))';ws'];
    end
    
    % Calculate piece-wise linear S
    nLoc = size(S,1);
    tmpI = [];
    tmpJ = [];
    tmpV = [];
  
    for i = 1:nLoc
        % Which cell
        cIdx = find(S(i,:) == 1);
        
        % Get vertexes
        v1 = vLoc(tv(cIdx,1),:);
        v2 = vLoc(tv(cIdx,2),:);
        v3 = vLoc(tv(cIdx,3),:);
        
        % Shift coordinate system
        vec1 = (v2-v1);
        vec2 = (v3-v1);
        
        % Calculate barycentric coordinates
        cord = abs([vec1',vec2']\(loc(i,:)'-v1'));
        cord = [1-sum(cord);cord];
        
        % create linear combination needed
        ttJ = [];
        ttV = [];
        for j = 1:3
            ttJ = [ttJ; vWs{tv(cIdx,j)}(1,:)'];
            ttV = [ttV; cord(j)*vWs{tv(cIdx,j)}(2,:)'];
        end
        ttI = ones(size(ttJ,1),1)*i;
        
        % Update full vectors
        tmpI = [tmpI;ttI];
        tmpJ = [tmpJ;ttJ];
        tmpV = [tmpV;ttV];
    end
    
    % Assemble matrix
    Slin = sparse(tmpI,tmpJ,tmpV, size(S,1), size(S,2));
end