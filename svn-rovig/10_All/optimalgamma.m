function gamma = optimalgamma(mesh,parameters)
% loop on all faces of GammaC
% then compute for each 
% 1) normal stress: gamma=1/6 * Area
% 2) normal displacement: gamma=1/6 * (sum Area which the vertex belong)

L=length(mesh);


for lev=1:L
N3=parameters.dim * mesh{lev}.N;
NF3=parameters.dim * mesh{lev}.NF;

gamma{lev}=zeros(N3+NF3,1);

%param=parameters.C_contact * 2/3.0;

param=parameters.C_contact * 1/3.0;

for ii = mesh{lev}.F_contact
    
    face=mesh{lev}.face(ii,:);
    nodes=mesh{lev}.node(face,:);
    A=AreaTriangle(nodes) ;
%     % first component of stress is the normal one
%     gamma{lev}(ii)=gamma{lev}(ii)+param*A;
%     % first component of displacement is the normal one
%     gamma{lev}(face+NF3)=gamma{lev}(face+NF3)+param*A;

    % first component of stress is the normal one
    gamma{lev}(ii)=gamma{lev}(ii)+3*param*A;
    % first component of displacement is the normal one
    gamma{lev}(face+NF3)=gamma{lev}(face+NF3)+param*A;

end

end


end