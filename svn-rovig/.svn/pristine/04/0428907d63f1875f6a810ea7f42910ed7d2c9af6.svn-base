function [dirichlet,n_and_or_t,bool_bc]= boundary_value_bool(type_of_dof)
% dirichlet contains the value of the bc
% n_and_or_t say if only the direction of the bc is fixed or also the orthogonal
% bool_bc defines which border is a real bc
contactYES=1;
contactNO=0;

if(type_of_dof==1)
    disp=1:12;
    disp=zeros(12,1);
  
value=-0.01;
    
    
disp1x=0.;
disp1y=0;
disp2x=0.0;%-0.05;
disp2y=0.;



disp3x=0.;
disp3y=value;
%disp3y=-0.05;

disp4x=0.0;
disp4y=0.;
    
dirichlet=[     disp1x disp1y contactYES;
                disp2x disp2y contactNO;
                disp3x disp3y contactNO;
                disp4x disp4y contactNO;
                0.0 0.0       contactYES;
                0.0 0.0       contactYES;];


% beam
%bool_bc=[0 0 0 1 0 0 0];

% L-beam
bool_bc=[0 0 0 0 0 1];

% dirichlet
%bool_bc=[1 1 1 1 0 0 0];

% contact cube
bool_bc=[0 0 1 0 0 0 0];
% contact circle
bool_bc=[0 0 1 0 0 0 0];



%bool_bc=[0 1 0 1 0 0 0];
%bool_bc=[0 0 0 0 0 0 0];

%bool_bc=[1 1 1 1 0 0 0];

%bool_bc=[0 0 0 1 0 0 0];
% bool_bc=[1 1 1 1 1 1];

%bool_bc=[0 0 0 1 1 1];

n_and_or_t=[1 1;
            1 1;
            1 1;
            1 1;
            1 1;
            1 1;]; 
       
else if(type_of_dof==2) 
        
value=-0.000005;
dirichlet=[     0.0 0.0 contactYES;
                0.0 0.0 contactNO;
                0.0 value contactNO;
                0.0 0.0 contactNO;
                0.0 value contactYES;
                0.0 0.0 contactYES;];



%beam
%bool_bc=[1 1 1 0 0 0 0];

% L-beam
bool_bc=[1 1 1 1 1 0];


% contact cube
bool_bc=[0 1 0 1 0 0 0];

% contact circle
%bool_bc=[0 1 0 0 0 0 0];




%dirichlet
%bool_bc=[ 0 0 0 0 0 0];



%bool_bc=[0 0 0 0 0 0];



% bool_bc=[1 1 1 1 0 0];

n_and_or_t=[0 0;
            0 0;
            0 0;
            0 0;
            0 0;
            0 0;];
 else
dirichlet=[     0.0 0.0 0.0;
                0.0 0.0 0.0;
                0.0 0.0 0.0;
                0.0 0.0 0.0;
                0.0 0.0 0.0;
                0.0 0.0 0.0;];
            
bool_bc=[0 0 0 0 0 0];

n_and_or_t=[     0 0 0;
           0 0 0;
           0 0 0;
           0 0 0;
           0 0 0;
           0 0 0;];
    end

end
end