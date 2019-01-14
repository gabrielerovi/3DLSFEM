function [RT0_mass,RT0_div,P1_Grad,RT0_Grad,RT0_divergence,signRT0_divergence,q_point,weights]=assemblingreference()

qrule=4;
node=[0     0     0;
      1     0     0;
      0     1     0;
      0     0     1];


[RT0_XX,RT0_XY,RT0_XZ,RT0_YX,RT0_YY,RT0_YZ,RT0_ZX,RT0_ZY,RT0_ZZ,RT0_div,RT0_divergence,signRT0_divergence,...
RT0_XGradX,RT0_XGradY,RT0_XGradZ,RT0_YGradX,RT0_YGradY,RT0_YGradZ,RT0_ZGradX,RT0_ZGradY,RT0_ZGradZ,...
GradXGradX,GradXGradY,GradXGradZ,GradYGradX,GradYGradY,GradYGradZ,GradZGradX,GradZGradY,GradZGradZ,q_point,weights] = phiRT3D(qrule,node);


q_point=q_point';

RT0_mass=[RT0_XX',RT0_XY',RT0_XZ',RT0_YX',RT0_YY',RT0_YZ',RT0_ZX',RT0_ZY',RT0_ZZ'];

P1_Grad=[GradXGradX,GradXGradY,GradXGradZ,GradYGradX,GradYGradY,GradYGradZ,GradZGradX,GradZGradY,GradZGradZ];

RT0_Grad=[RT0_XGradX,RT0_XGradY,RT0_XGradZ,RT0_YGradX,RT0_YGradY,RT0_YGradZ,RT0_ZGradX,RT0_ZGradY,RT0_ZGradZ];



end
