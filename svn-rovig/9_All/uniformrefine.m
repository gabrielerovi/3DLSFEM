function [node,elem,T_edge_bcF,boundaryF] = uniformrefine(node,elem,T_edge_bcC)
%% UNIFORMREFINE  uniformly refine a 2-D triangulation.
%
% [node,elem] = uniformrefine(node,elem) divides each triangle into 4 small
% similar triangles.
%
% [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag) also update boundary
% conditions represented by bdFlag.
%
% [node,elem,~,HB] = uniformrefine(node,elem) outpus HB array which is
% useful for nodal interpolation.
%
% Example
%
%      node = [0,0; 1,0; 1,1; 0,1];
%      elem = [2,3,1; 4,1,3];
%      figure(1); subplot(1,3,1); showmesh(node,elem);
%      [node,elem] = uniformrefine(node,elem);
%      figure(1); subplot(1,3,2); showmesh(node,elem);
%      bdFlag = setboundary(node,elem,'Dirichlet','all','Neumann','y==1');
%      [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
%      figure(1); subplot(1,3,3); showmesh(node,elem);
%
% See also uniformrefine3, uniformbisect, bisect, bisect3
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Construct data structure

elem=sort(elem,2);
totalEdge = uint32(sort([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])],2));
matlabversion = version;
if str2double(matlabversion(end-5:end-2)) > 2012
    [edge, tempvar, j] = unique(totalEdge,'rows','legacy'); %#ok<*ASGLU>
else
    [edge, tempvar, j] = unique(totalEdge,'rows'); %#ok<*ASGLU>
end
N = size(node,1); NT = size(elem,1); NE = size(edge,1);
elem2edge = uint32(reshape(j,NT,3));

%% Add new nodes: middle points of all edges
node(N+1:N+NE,:) = (node(edge(:,1),:)+node(edge(:,2),:))/2; 


HB(:,[1 2 3]) = [(N+1:N+NE)', edge(:,1:2)]; 
edge2newNode = uint32((N+1:N+NE)');

%% Refine each triangle into four triangles as follows
% 3
% | \
% 5- 4
% |\ |\
% 1- 6- 2
t = 1:NT;
 p(t,1:3) = elem(t,1:3);
 p(t,4:6) = edge2newNode(elem2edge(t,1:3));
 boundaryF=[];
% elem(t,:) = [p(t,1), p(t,6), p(t,5)];
% elem(NT+1:2*NT,:) = [p(t,6), p(t,2), p(t,4)];
% elem(2*NT+1:3*NT,:) = [p(t,5), p(t,4), p(t,3)];
% elem(3*NT+1:4*NT,:) = [p(t,4), p(t,5), p(t,6)];


%check if there are internal elements
T_edge_bcF=zeros(4*NT,3);
for t=1:NT
    
elem(4*t-3,:) = [p(t,1), p(t,6), p(t,5)];
elem(4*t-2,:) = [p(t,6), p(t,2), p(t,4)];
elem(4*t-1,:) = [p(t,5), p(t,4), p(t,3)];
elem(4*t,:) = [p(t,4), p(t,5), p(t,6)];
T_edge_bcF(4*t-3,:)=[T_edge_bcC(t,1),T_edge_bcC(t,2),0];
T_edge_bcF(4*t-2,:)=[T_edge_bcC(t,1),T_edge_bcC(t,3), 0];
T_edge_bcF(4*t-1,:)=[T_edge_bcC(t,2),T_edge_bcC(t,3), 0];
T_edge_bcF(4*t,:)=[0,0,0];



        if(T_edge_bcC(t,1) >0)
         boundaryF(end+1,:)=[p(t,1), p(t,6),T_edge_bcC(t,1)];
         boundaryF(end+1,:)=[p(t,2), p(t,6),T_edge_bcC(t,1)];
        end
        if(T_edge_bcC(t,2) >0)
         boundaryF(end+1,:)=[p(t,1), p(t,5),T_edge_bcC(t,2)];
         boundaryF(end+1,:)=[p(t,3), p(t,5),T_edge_bcC(t,2)];
        end
        if(T_edge_bcC(t,3) >0)
         boundaryF(end+1,:)=[p(t,2), p(t,4),T_edge_bcC(t,3)];
         boundaryF(end+1,:)=[p(t,3), p(t,4),T_edge_bcC(t,3)];
        end        

end



elem=sort(elem,2);


% %% Update boundary edges
% if exist('bdFlag','var') && ~isempty(bdFlag)
%     bdFlag(NT+1:2*NT,[1 3]) = bdFlag(t,[1 3]); 
%     bdFlag(2*NT+1:3*NT,[1 2]) = bdFlag(t,[1 2]); 
%     bdFlag(3*NT+1:4*NT,1) = 0;
%     bdFlag(t,1) = 0;
% else
%     bdFlag = [];
% end