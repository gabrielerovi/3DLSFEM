function [RTCtoRTF,P1CtoP1F,NDCtoNDF,NDtoRT]=projections_levels2D(mesh,RTCtoRTF,P1CtoP1F,NDCtoNDF,NDtoRT)

L=size(mesh);
L=L(1);

% edge to remove because we impose Dirichlet
E_remove=mesh{L}.E_remove;
% node to remove because we impose Dirichlet
N_remove=mesh{L}.N_remove;
% node to remove because we impose Dirichlet on the relative edge
type_of_dof=2;


NDtoRT{L}(E_remove,:)=[];
%NDtoRT{L}(:,Ned_remove)=[];

for lev=1:L-1
    NC_remove=mesh{lev}.N_remove; 
    EC_remove=mesh{lev}.E_remove;
    
    NF_remove=mesh{lev+1}.N_remove;
    EF_remove=mesh{lev+1}.E_remove;
    
    RTCtoRTF{lev}(EF_remove,:)=[];
    RTCtoRTF{lev}(:,EC_remove)=[];
    
    P1CtoP1F{lev}(NF_remove,:)=[];
    P1CtoP1F{lev}(:,NC_remove)=[];
    
%    NDCtoNDF{lev}(Ned_remove,:)=[];
%    NDCtoNDF{lev}(:,Ned_remove)=[];
    
    NDtoRT{lev}(EC_remove,:)=[];
%    NDtoRT{lev}(:,Ned_remove)=[];
    
end



end