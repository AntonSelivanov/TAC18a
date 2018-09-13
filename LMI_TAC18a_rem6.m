function [hfeas,L]=LMI_TAC18a_rem6(A,C,h,delta,epsilon)
% This MATLAB program checks the feasibility of LMIs from Remark 6 of the paper 
% A. Selivanov and E. Fridman, "Boundary observers for a reaction-diffusion system under time-delayed and sampled-data measurements," IEEE Transactions on Automatic Control, 2018.

% The program uses YALMIP parser (http://users.isy.liu.se/johanl/yalmip/)

% Input: 
% A, C      - parameters from (23)
% h         - upper sampling bound from (47)
% delta     - decay rate 
% epsilon   - tuning parameter

% Output: 
% hfeas     =h if feasible, =0 if not feasible
% L         - injection gain 

N=size(A,1); % the number of modes
%% Decision variables 
P=sdpvar(N); 
W=sdpvar(N); 
P2=sdpvar(N,N,'f'); 
Y=sdpvar(N,1,'f'); 
%% LMIs 
Ups=blkvar; 
Ups(1,1)=A'*P2+P2'*A-C'*Y'-Y*C+2*delta*P; 
Ups(1,2)=P-P2'+epsilon*A'*P2-epsilon*C'*Y'; 
Ups(1,3)=-Y*C; 
Ups(2,2)=-epsilon*P2-epsilon*P2'+h^2*exp(2*delta*h)*W; 
Ups(2,3)=-epsilon*Y*C; 
Ups(3,3)=-pi^2/4*W; 
Ups=sdpvar(Ups); 
%% Solution of LMIs
LMIs=[P>=0,W>=0,Ups<=0]; 

options=sdpsettings('solver','lmilab','verbose',0); 
sol=optimize(LMIs,[],options); 

hfeas=0; L=[]; 
if sol.problem==0
    [primal,~]=check(LMIs); 
    if min(primal)>0
        hfeas=h; 
        P2v=double(P2); 
        L=(P2v')^(-1)*double(Y);         
    end
else
    yalmiperror(sol.problem) 
end