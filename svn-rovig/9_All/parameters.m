function parameters =parameters()


% LSstressblock, LSpoisson,
% LSelasticity,LSelasticityAsymmetric,DispElasticity
parameters.input_name='LSelasticity';
%parameters.input_mesh='Cookmesh.msh';
parameters.input_mesh='SemiCircular.msh';
parameters.create_only_square=true;
parameters.a=1;
parameters.b=1;
parameters.COARSE=2;
parameters.FINE=3;
parameters.number_of_levels=1+parameters.FINE-parameters.COARSE;
parameters.L=parameters.number_of_levels;

parameters.components=2;
parameters.N_components=2;
parameters.E_components=2;
parameters.F_components=0;
parameters.smoothing_steps=5;
parameters.dim=2;
parameters.qrule=3;


parameters.E=20;%70
parameters.nu=0.248;%499;%0.33;%0.499;
parameters.K=parameters.E/(3 * (1-2*parameters.nu) );

parameters.max_iter=50;
parameters.toll_loop=10^(-10);
parameters.lambda = 1;%10^50;%parameters.E * parameters.nu / ( (1+parameters.nu) * (1-2* parameters.nu ) ) ;
parameters.mu = 1; %3/2 * (parameters.K-parameters.lambda) ;


parameters.beta = 1.0 / ( 2.0 * parameters.mu);
parameters.alpha = - parameters.beta*parameters.lambda/(parameters.dim *parameters.lambda + 2.0 * parameters.mu);
parameters.gamma=parameters.lambda+2*parameters.mu;
 
%  parameters.force1=@(x,y) ((1.0)*( 4.0 * pi*pi * sin(pi * x) * sin(pi * y) - 2.0 * pi*pi *cos( pi * x) * cos( pi * y) ));
%  parameters.force2=@(x,y) ((1.0)*( 4.0 * pi*pi * sin(pi * x) * sin(pi * y) - 2.0 * pi*pi *cos( pi * x) * cos( pi * y) ));
parameters.force1=@(x,y) (0.0);
parameters.force2=@(x,y) (0.0);
parameters.fzero=@(x,y) (0.0);

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
parameters.body_normal_x=@(x,y) 0.0;
parameters.body_normal_y=@(x,y) -1.0;


parameters.GS_is_symmetric=false;
parameters.MGiter=100;
parameters.MGiterSS=20;
parameters.MGiterUU=20;
parameters.toll=10^(-9);
parameters.LineSearch=false;
parameters.h_dependence=false;

parameters.C_contact=1.0 * 10^1;
parameters.C_eq=1.0 * 10^2;
parameters.C_const=1.0;
parameters.C_asym=0.0 * 10^(1.0);

% Benchmarl dirichlet solution
parameters.disp1=@(x,y) (sin(pi*x)*sin(pi*y));
parameters.disp2=@(x,y) (sin(pi*x)*sin(pi*y));
phi1=@(x,y) pi .* cos(pi*x)*sin(pi*y);
phi2=@(x,y) pi .* sin(pi*x)*cos(pi*y);
phi3=@(x,y) pi .* sin(pi*x)*cos(pi*y);
phi4=@(x,y) pi .* cos(pi*x)*sin(pi*y);

parameters.stressxx=@(x,y)(parameters.gamma * pi .* cos(pi*x)*sin(pi*y)  + parameters.lambda * pi .* sin(pi*x)*cos(pi*y));
parameters.stressxy=@(x,y)(parameters.mu.*phi3(x,y)     + parameters.mu     .*phi4(x,y));
parameters.stressyx=@(x,y)(parameters.mu.*phi3(x,y)     + parameters.mu     .*phi4(x,y));
parameters.stressyy=@(x,y)(parameters.lambda.*phi1(x,y) + parameters.gamma .*phi2(x,y));
end
