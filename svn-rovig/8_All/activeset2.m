function [u,lambda,WorkingSet] = activeset2(A,STT,BT,f,ST,CT,B,C,xold)

% we solve for min H, with H=0.5 x' A x - x' f - lambda (B x -c)
% structure of the problem
% |A -B'| |x     |= |f|
% |B  0 | |lambda|  |c|
% in particular ST= matrix for quality constraints, CT= corresponding rhs
% in particular BT= matrix for inequality constraints, C= corresponding rhs


toll=10^(-12);
n=length(f);
nconstraint=length(C);


continueplease=true;

check=B*xold-C;
cont=0;
WorkingSet=[];
constraint_tot=[];




 while(continueplease)
%         
BT_k = BT(:,WorkingSet);   
B_k = B(WorkingSet,:);
c_k = C(WorkingSet);
% 
% % LW number of active inequality constraints
% % LS number of equality constraints (always active)
% 
LW=length(WorkingSet);
LS=length(CT);
LWS=LW+LS;
% 
M_k=[A   -STT -BT_k;
    ST   sparse(LS,LWS)
    B_k   sparse(LW,LWS);];


F_k=[f;CT;sparse(length(WorkingSet),1)];

correction=M_k\F_k;

xold=correction(1:n);

c_lambda_EQ=correction(1+n:n+LS);

c_lambda_k=correction(1+n+LS:end);


[minlambda, pos_minlambda]=min(c_lambda_k);

% check if, after solving the problem, all the constraints are satisfied
% if not, then remove the constraint that has the most negative multiplier
if(minlambda<-toll)
   WorkingSet( pos_minlambda ) = [ ];
else
% if all the multipliers are satisfied, then the correction is meaningful,
% but we have to find a proper alpha such that this correction can satisfy
% also the other constraints not belonging to the working set

% if check(i)<0, then the i constraint is not fullfilled
% since minlambda>0, the check corresponding to the WorkingSet should be
% zero
check=B*(xold)-C;
blocking_constraint=[];
NotWorkingSet=setdiff(1:nconstraint,WorkingSet);
alpha=1;
WorkingSetOld=WorkingSet;
WorkingSet=find(check<-toll);

  if(isempty(WorkingSet) )
      continueplease=false;
  end
  
    
    
end
    

    

end

WorkingSet=WorkingSetOld;

J=[f;CT;sparse(length(WorkingSet),1)]

sol_tot=[xold;c_lambda_k];
% norm(J-M_k*sol_tot)

 u=xold;
 lambda=c_lambda_k;
 end

