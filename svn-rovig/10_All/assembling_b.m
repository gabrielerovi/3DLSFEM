function F=assembling_b(qrule,node,face_per_elem,fx,fy,fz,coeff_equilibrium)


[q_point,weights,Volume]=quadrature_points_3D(qrule,node);
[RT0,RT0_div] = phiRT3Dcell(q_point,node);
number_of_qp=length(q_point(:,1));
value=zeros(number_of_qp,1);
F=cell(3,1);

   for kk=1:3
       F{kk}=zeros(face_per_elem,1);
   end

   for nn=1:face_per_elem 

            for qp=1:number_of_qp
                F{1}(nn,1)= F{1}(nn,1)+ weights(qp) * (-1.0)*coeff_equilibrium * RT0_div(nn) *fx(q_point(qp,1),q_point(qp,2),q_point(qp,3));
                F{2}(nn,1)= F{2}(nn,1)+ weights(qp) * (-1.0)*coeff_equilibrium * RT0_div(nn) *fy(q_point(qp,1),q_point(qp,2),q_point(qp,3));
                F{3}(nn,1)= F{3}(nn,1)+ weights(qp) * (-1.0)*coeff_equilibrium * RT0_div(nn) *fz(q_point(qp,1),q_point(qp,2),q_point(qp,3));
            end           
   end

    F{1}= Volume * F{1};              
    F{2}= Volume * F{2};   
    F{3}= Volume * F{3};   

    
end

    