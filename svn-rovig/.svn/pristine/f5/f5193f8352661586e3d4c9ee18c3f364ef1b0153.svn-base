function parameters =parameters()


% LSstressblock, LSpoisson,
% LSelasticity,LSelasticityAsymmetric,DispElasticity
parameters.input_name='LSelasticity';
%parameters.input_mesh='Cookmesh.msh';
% hfine = 0.01;
% hcoarse=0.1;
parameters.input_mesh='sphere2.msh';
parameters.create_only_square=true;
parameters.a=1;
parameters.b=1;
parameters.COARSE=1;
parameters.FINE=5;
parameters.number_of_levels=1+parameters.FINE-parameters.COARSE;
parameters.L=parameters.number_of_levels;

parameters.components=3;
parameters.N_components=3;
parameters.E_components=3;
parameters.F_components=3;
parameters.smoothing_steps=5;
parameters.dim=3;
parameters.qrule=2;


parameters.E=20;%70
parameters.nu=0.248;%499;%0.33;%0.499;
parameters.K=parameters.E/(3 * (1-2*parameters.nu) );

parameters.max_iter=150;
parameters.toll_loop=10^(-10);
parameters.lambda = 10^50;%4.038;  %1;%10^50;%parameters.E * parameters.nu / ( (1+parameters.nu) * (1-2* parameters.nu ) ) ;
parameters.mu = 1; %3/2 * (parameters.K-parameters.lambda) ;


parameters.beta = 1.0 / ( 2.0 * parameters.mu);
parameters.alpha = - parameters.beta*parameters.lambda/(parameters.dim *parameters.lambda + 2.0 * parameters.mu);
parameters.gamma=parameters.lambda+2*parameters.mu;
 mu= parameters.mu;
 lambda=parameters.lambda;
 parameters.force1=@(x,y,z) (5*pi^2*sin(pi*x).*sin(pi*y).*sin(pi*z) - 2*pi^2*cos(pi*x).*cos(pi*z).*sin(pi*y) - 2*pi^2*cos(pi*x).*cos(pi*y).*sin(pi*z));
 parameters.force2=@(x,y,z) (5*pi^2*sin(pi*x).*sin(pi*y).*sin(pi*z) - 2*pi^2*cos(pi*y).*cos(pi*z).*sin(pi*x) - 2*pi^2*cos(pi*x).*cos(pi*y).*sin(pi*z));
 parameters.force3=@(x,y,z) (5*pi^2*sin(pi*x).*sin(pi*y).*sin(pi*z) - 2*pi^2*cos(pi*y).*cos(pi*z).*sin(pi*x) - 2*pi^2*cos(pi*x).*cos(pi*z).*sin(pi*y));
% parameters.force1=@(x,y,z) (x*0.0);
% parameters.force2=@(x,y,z) (y*0.0);
% parameters.force3=@(x,y,z) (z*0.0);
parameters.istheexternalforcenonzero=false;

% gap function for sphere (body normal)
%parameters.Radius=0.4;
% parameters.gap=@(x,y) (parameters.Radius)/(cos(atan(abs(x/(parameters.Radius-y))))+10^(-30))-(parameters.Radius);

% gap function for sphere (obstacle normal)
parameters.gap=@(x,y) parameters.Radius*( 1- sin(acos (x/parameters.Radius) ) );
% gap function for square
%parameters.gap=@(x,y) 0.0;

parameters.hertz=1;
parameters.frictionless=1;
parameters.contact=1;
parameters.penalty=0;
parameters.body_normal_bool=false;
parameters.body_normal_x=@(x,y,z) 0.0;
parameters.body_normal_y=@(x,y,z) 0.0;
parameters.body_normal_z=@(x,y,z) -1.0;


parameters.GS_is_symmetric=false;
parameters.MGiter=100;
parameters.MGiterSS=20;
parameters.MGiterUU=20;
parameters.toll=10^(-9);
parameters.LineSearch=false;
parameters.h_dependence=false;

parameters.C_contact=1 * 10^1;
parameters.C_eq=1.0 * 10^2;
parameters.C_const=1.0;
parameters.C_asym=0.0 * 10^(0.0);

% Benchmarl dirichlet solution
parameters.disp1=@(x,y,z) (sin(pi * x) .* sin(pi * y) .* sin(pi * z));
parameters.disp2=@(x,y,z) (sin(pi * x) .* sin(pi * y) .* sin(pi * z));
parameters.disp3=@(x,y,z) (sin(pi * x) .* sin(pi * y) .* sin(pi * z));


parameters.dispgradx = @(x,y,z) (pi*cos(pi*x).*sin(pi*y).*sin(pi*z));
parameters.dispgrady = @(x,y,z) (pi*cos(pi*y).*sin(pi*x).*sin(pi*z));
parameters.dispgradz = @(x,y,z) (pi*cos(pi*z).*sin(pi*x).*sin(pi*y));
 
parameters.stressxx=@(x,y,z)(  lambda*(pi*cos(pi*x).*sin(pi*y).*sin(pi*z) + pi*cos(pi*y).*sin(pi*x).*sin(pi*z) + pi*cos(pi*z).*sin(pi*x).*sin(pi*y)) + 2*mu*pi*cos(pi*x).*sin(pi*y).*sin(pi*z)  );
parameters.stressxy=@(x,y,z)(  2*mu*((pi*cos(pi*x).*sin(pi*y).*sin(pi*z))/2 + (pi*cos(pi*y).*sin(pi*x).*sin(pi*z))/2) );
parameters.stressxz=@(x,y,z)(  2*mu*((pi*cos(pi*x).*sin(pi*y).*sin(pi*z))/2 + (pi*cos(pi*z).*sin(pi*x).*sin(pi*y))/2) );
parameters.stressyx=@(x,y,z)(  2*mu*((pi*cos(pi*x).*sin(pi*y).*sin(pi*z))/2 + (pi*cos(pi*y).*sin(pi*x).*sin(pi*z))/2) );
parameters.stressyy=@(x,y,z)(  lambda*(pi*cos(pi*x).*sin(pi*y).*sin(pi*z) + pi*cos(pi*y).*sin(pi*x).*sin(pi*z) + pi*cos(pi*z).*sin(pi*x).*sin(pi*y)) + 2*mu*pi*cos(pi*y).*sin(pi*x).*sin(pi*z) );
parameters.stressyz=@(x,y,z)(  2*mu*((pi*cos(pi*y).*sin(pi*x).*sin(pi*z))/2 + (pi*cos(pi*z).*sin(pi*x).*sin(pi*y))/2)  );
parameters.stresszx=@(x,y,z)(  2*mu*((pi*cos(pi*x).*sin(pi*y).*sin(pi*z))/2 + (pi*cos(pi*z).*sin(pi*x).*sin(pi*y))/2)  );
parameters.stresszy=@(x,y,z)(  2*mu*((pi*cos(pi*y).*sin(pi*x).*sin(pi*z))/2 + (pi*cos(pi*z).*sin(pi*x).*sin(pi*y))/2)  );
parameters.stresszz=@(x,y,z)(  lambda*(pi*cos(pi*x).*sin(pi*y).*sin(pi*z) + pi*cos(pi*y).*sin(pi*x).*sin(pi*z) + pi*cos(pi*z).*sin(pi*x).*sin(pi*y)) + 2*mu*pi*cos(pi*z).*sin(pi*x).*sin(pi*y)  );

parameters.divstress1=@(x,y,z) (lambda*(pi^2*cos(pi*x).*cos(pi*y).*sin(pi*z) + pi^2*cos(pi*x).*cos(pi*z).*sin(pi*y) - pi^2*sin(pi*x).*sin(pi*y).*sin(pi*z)) + 2*mu*((pi^2*cos(pi*x).*cos(pi*y).*sin(pi*z))/2 - (pi^2*sin(pi*x).*sin(pi*y).*sin(pi*z))/2) + 2*mu*((pi^2*cos(pi*x).*cos(pi*z).*sin(pi*y))/2 - (pi^2*sin(pi*x).*sin(pi*y).*sin(pi*z))/2) - 2*mu*pi^2*sin(pi*x).*sin(pi*y).*sin(pi*z));
parameters.divstress2=@(x,y,z) (lambda*(pi^2*cos(pi*x).*cos(pi*y).*sin(pi*z) + pi^2*cos(pi*y).*cos(pi*z).*sin(pi*x) - pi^2*sin(pi*x).*sin(pi*y).*sin(pi*z)) + 2*mu*((pi^2*cos(pi*x).*cos(pi*y).*sin(pi*z))/2 - (pi^2*sin(pi*x).*sin(pi*y).*sin(pi*z))/2) + 2*mu*((pi^2*cos(pi*y).*cos(pi*z).*sin(pi*x))/2 - (pi^2*sin(pi*x).*sin(pi*y).*sin(pi*z))/2) - 2*mu*pi^2*sin(pi*x).*sin(pi*y).*sin(pi*z));
parameters.divstress3=@(x,y,z) (lambda*(pi^2*cos(pi*x).*cos(pi*z).*sin(pi*y) + pi^2*cos(pi*y).*cos(pi*z).*sin(pi*x) - pi^2*sin(pi*x).*sin(pi*y).*sin(pi*z)) + 2*mu*((pi^2*cos(pi*x).*cos(pi*z).*sin(pi*y))/2 - (pi^2*sin(pi*x).*sin(pi*y).*sin(pi*z))/2) + 2*mu*((pi^2*cos(pi*y).*cos(pi*z).*sin(pi*x))/2 - (pi^2*sin(pi*x).*sin(pi*y).*sin(pi*z))/2) - 2*mu*pi^2*sin(pi*x).*sin(pi*y).*sin(pi*z));


end

