function hfeas=LMI_TAC18a_th2(A,C,L,h,delta)
% This MATLAB program checks the feasibility of LMIs from Theorem 1 of the paper 
% A. Selivanov and E. Fridman, "Boundary observers for a reaction-diffusion system under time-delayed and sampled-data measurements," IEEE Transactions on Automatic Control, 2018.

% The program uses YALMIP parser (http://users.isy.liu.se/johanl/yalmip/)

% Input: 
% A, C, L  - parameters from (23)
% h        - upper sampling bound from (47)
% delta    - decay rate 

% Output: 
% hfeas =h if feasible, =0 if not feasible

N=size(A,1); 
%% Decision variables 
P=sdpvar(N); 
W=sdpvar(N); 
P2=sdpvar(N,N,'f'); 
P3=sdpvar(N,N,'f'); 
%% LMIs 
Ups=blkvar; 
Ups(1,1)=(A-L*C)'*P2+P2'*(A-L*C)+2*delta*P; 
Ups(1,2)=P-P2'+(A-L*C)'*P3; 
Ups(1,3)=-P2'*L*C; 
Ups(2,2)=-P3-P3'+h^2*exp(2*delta*h)*W; 
Ups(2,3)=-P3'*L*C; 
Ups(3,3)=-pi^2/4*W; 
Ups=sdpvar(Ups); 
%% Solution of LMIs
LMIs=[Ups<=0,P>=0,W>=0]; 

options=sdpsettings('solver','lmilab','verbose',0); 
sol=optimize(LMIs,[],options); 

hfeas=0; 
if sol.problem==0
    [primal,~]=check(LMIs); 
    if min(primal)>=0
        hfeas=h; 
    end
else
    yalmiperror(sol.problem) 
end