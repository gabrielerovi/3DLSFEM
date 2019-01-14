function [dirichlet,stronger_dirichlet,bool_bc]= boundary_value_bool(type_of_dof)
% dirichlet contains the value of the bc
% n_and_or_t say if only the direction of the bc is fixed or also the orthogonal
% bool_bc defines which border is a real bc
contactYES=1;
contactNO=0;

if(type_of_dof==1)
    disp=1:12;
    disp=zeros(12,1);
  
value=-0.01;
    
%%%%%%%%%% SPHERE NO NEUMANN %%%%%%%%%%  
% dirichlet=[     0 0 0     contactYES;
%                 0 0 value     contactNO;];
% % %%%%%%%%%% SPHERE NEUMANN %%%%%%%%%%  
% dirichlet=[     0 0 value     contactNO;
%                 0 0 0     contactNO;
%                 0 0 0     contactYES];
%%%%%%%%%% CONTACT CUBE %%%%%%%%%%
dirichlet=[     0 0 0     contactYES;
                0 0 0     contactNO;
                0 0 0     contactNO;
                0 0 0     contactNO;
                0 0 0     contactNO;
                0 0 value contactNO;];
% %%%%%%%%%% CUBE DIRICHLET %%%%%%%%%%
% dirichlet=[     0 0 0.     contactNO;
%                 0 0 0     contactNO;
%                 0 0 0     contactNO;
%                 0 0 0     contactNO;
%                 0 0 0     contactNO;
%                 0 0 0.     contactNO;];
%             
            
% contact no neumann sphere
% bool_bc=[0 1];
% contact neumann sphere
bool_bc=[1 0 0];

% CONTACT CUBE cube
bool_bc=[0 0 0 0 0 1];
% linear dirichlet cube
% bool_bc=[1 1 1 1 1 1];

% contact no neumann sphere
% stronger_dirichlet=[0; 1; ]; 

% contact neumann sphere
% stronger_dirichlet=[1; 0; 0];

% CONTACT CUBE 
stronger_dirichlet=[0;
            0;
            0;
            0;
            0;
            1;]; 
% stronger_dirichlet=[1;
%             0;
%             0;
%             0;
%             0;
%             1;];       
else if(type_of_dof==3) 
% %%%%%%%%%% SPHERE NO NEUMANN %%%%%%%%%%  
% dirichlet=[     0.0 0.0 0.0  contactYES;
%                 0.0 0.0 0.0  contactNO;];
% %%%%%%%%%% SPHERE NEUMANN %%%%%%%%%%  
% dirichlet=[     0.0 0.0 0.0  contactNO;
%                 0.0 0.0 0.0  contactNO;
%                 0.0 0.0 0.0  contactYES;];
            
% %%%%%%%%%% CONTACT CUBE  %%%%%%%%%%
 dirichlet=[    0.0 0.0 0.0  contactYES;
                0.0 0.0 0.0  contactNO;
                0.0 0.0 0.0  contactNO;
                0.0 0.0 0.0  contactNO;
                0.0 0.0 0.0  contactNO;
                0.0 0.0 0.0  contactNO;];  

% %%%%%%%%%% LINEAR CUBE DIRICHLET %%%%%%%%%%
%  dirichlet=[     0.0 0.0 0.0  contactNO;
%                 0.0 0.0 0.0  contactNO;
%                 0.0 0.0 0.0  contactNO;
%                 0.0 0.0 0.0  contactNO;
%                 0.0 0.0 0.0  contactNO;
%                 0.0 0.0 0.0  contactNO;];   
% contact no neumann sphere
% bool_bc=[0 0];

% contact neumann sphere
% bool_bc=[0 1 0 ];

% contact cube
bool_bc=[0 1 1 1 1 0];
% linear dirichlet cube
% bool_bc=[0 0 0 0 0 0]; 
  
  
stronger_dirichlet=[     0 0 0;
           0 0 0;
           0 0 0;
           0 0 0;
           0 0 0;
           0 0 0;];
    end

end
end