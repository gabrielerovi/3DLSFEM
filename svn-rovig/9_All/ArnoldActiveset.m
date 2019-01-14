


function [u,lambda,WorkingSet] = ArnoldActiveset(A,B,b,c,WorkingSet)

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

 while(continueplease)   
B_k = B(WorkingSet,:);
c_k = c(WorkingSet);

LW=length(WorkingSet);

M_k=[A    B_k';
     B_k  sparse(LW,LW);];


F_k=[b;c_k];


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

%correction=M_k\F_k;

u_k=correction(1:n);

lambda_k=correction(1+n:end);
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





