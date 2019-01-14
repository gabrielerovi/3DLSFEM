


function [u,lambda,WorkingSet] = ArnoldActivesetAugmented(A,B,b,c,WorkingSet,gamma)

% we solve for min H, with H=0.5 x' A x - x' f - lambda (c-B x)
% structure of the problem
% |A B'| |x     |= |b|
% |B 0 | |lambda|  |c|
% in particular ST= matrix for quality constraints, CT= corresponding rhs
% in particular BT= matrix for inequality constraints, C= corresponding rhs


toll=10^(-9);
n=length(b);
nconstraint=length(c);


continueplease=true;


%WorkingSet=[];
constraint_tot=[];
% [LA,UA] = lu(full(A));
% toll=10^(-12);
% omegafixed=0.5;
Apenalty=sparse(n,n);
Rhspenaltytmp=zeros(n,1);
gamma_zero=zeros(n,1);
 while(continueplease)   
B_k = B(WorkingSet,:);
c_k = c(WorkingSet);

LW=length(WorkingSet);


if(isempty(WorkingSet))
    M_k=[A    B_k';
         B_k  sparse(LW,LW);];
 F_k=[b;c_k];
else
    
    Apenalty(:,:)=0;
for ii=1:length(WorkingSet)
Apenalty(WorkingSet(ii),WorkingSet(ii))=gamma(WorkingSet(ii));
end


Rhspenalty=Rhspenaltytmp+Apenalty*B_k'*c_k;
M_k=[A+Apenalty    B_k';
     B_k           sparse(LW,LW);];
  F_k=[b+Rhspenalty;c_k];
end




% use iterative refinement to compute the solution of the system:
% M_k * correction= F_k
% by doing 
% 1) M_k * xx= F_k
% 2) M_k * cc= F_k - M_k * xx
% 3) correction=xx+cc


[L,U] = lu(full(M_k));
yy=L\F_k;
correction=U\yy;

for hh=1:1
res=F_k-M_k*correction;
yy=L\res;
cc=U\yy;
correction=correction+cc;
end
u_k=correction(1:n);
lambda_k=correction(1+n:end);




% sizeB=size(B_k);
% if(sizeB(1)==0)
% y1=LA\b;
% x1=UA\y1;   
% lambda_new=[];
% else
% lambda_new=zeros(sizeB(1),1);
% res=toll+1;
% resold=res;
% res1=toll+1;
% omega=omegafixed;
% while(res>toll && norm(res1)>toll)
% F=b-B_k'*lambda_new;
% y1=LA\F;
% z1=UA\y1;
% res1=B_k*z1-c_k;
% lambda_old=lambda_new+omega*(res1);
% res2=b-A*z1-B_k'*lambda_old;
% res=norm([res1;res2]);
% 
% if(res>resold && omega>0.01)
% res=resold;
% omega=omega/2;
% else
% omega=omega;
% x1=z1; 
% lambda_new=lambda_old;
% end
% 
% end
% %correction=M_k\F_k;
% end
% norm(u_k-x1)
% norm(lambda_k-lambda_new)



WorkingSetOld=WorkingSet;
check=c-B*(u_k);
lambda_neg=find(lambda_k <-toll);
check_neg=find(check<-toll);

% remove negative multipliers
NotWorkingSet=WorkingSet(lambda_neg);
WorkingSet(lambda_neg)=[];
% add non feasible constraint
WorkingSet=unique([WorkingSet;check_neg]);
WorkingSet=setdiff(WorkingSet,NotWorkingSet);


if( isempty(check_neg) && isempty(lambda_neg)) 
    continueplease=false;
end
end

 u=u_k;
 lambda=lambda_k;
 end





