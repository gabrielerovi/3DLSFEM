function [A,F,C]=system_divdiv2D(mesh,f1,f2,qrule)

NE=mesh.NE;
N=mesh.N;
NT=mesh.NT;
A=zeros(NE,NE);
C=zeros(N,N);
F=zeros(NE,1);


    
for t=1:NT
    
    nodet=mesh.elem(t,:);
    edget=mesh.elemE(t,:);
    node=mesh.node(nodet,:);
    
    Massloc=local_mass_matrix(node);
    Divloc=local_div_matrix(node);
    M_loc=Massloc+Divloc;
    
    
    A(edget,edget)=A(edget,edget)+M_loc;
    
     F_loc = rhs_raviart(node,qrule,f1,f2);
     
     
     F(edget)=F(edget)+F_loc;

    Curl_loc = curl_curl_lagrangian  (node);
    C(nodet,nodet)=C(nodet,nodet)+Curl_loc;

end


end