function [taufeas,L]=LMI_TAC18a_rem4(A,C,tauM,delta,epsilon)
% This MATLAB program checks the feasibility of LMIs from Remark 4 of the paper 
% A. Selivanov and E. Fridman, "Boundary observers for a reaction-diffusion system under time-delayed and sampled-data measurements," IEEE Transactions on Automatic Control, 2018.

% The program uses YALMIP parser (http://users.isy.liu.se/johanl/yalmip/)

% Input: 
% A, C      - parameters from (23)
% tauM      - upper delay bound for (4)
% delta     - decay rate 
% epsilon   - tuning parameter

% Output: 
% taufeas   =tauM if feasible, =0 if not feasible
% L         - injection gain 

N=size(A,1); % the number of modes
%% Decision variables 
P=sdpvar(N); 
S=sdpvar(N); 
R=sdpvar(N); 
P2=sdpvar(N,N,'f'); 
Y=sdpvar(N,1,'f'); 
G=sdpvar(N,N,'f'); 
%% LMIs 
Phi=blkvar; 
Phi(1,1)=A'*P2+P2'*A+2*delta*P+S-exp(-2*delta*tauM)*R; 
Phi(1,2)=P-P2'+epsilon*A'*P2; 
Phi(1,3)=exp(-2*delta*tauM)*(R-G)-Y*C; 
Phi(1,4)=exp(-2*delta*tauM)*G; 
Phi(2,2)=-epsilon*P2-epsilon*P2'+tauM^2*R; 
Phi(2,3)=-epsilon*Y*C; 
Phi(3,3)=-exp(-2*delta*tauM)*(2*R-G-G'); 
Phi(3,4)=exp(-2*delta*tauM)*(R-G); 
Phi(4,4)=-exp(-2*delta*tauM)*(S+R); 
Phi=sdpvar(Phi); 

Park=[R G; G' R]; 
%% Solution of LMIs
LMIs=[Phi<=0,Park>=0,P>=0,S>=0,R>=0]; 

options=sdpsettings('solver','lmilab','verbose',0); 
sol=optimize(LMIs,[],options); 

taufeas=0; L=[]; 
if sol.problem==0
    [primal,~]=check(LMIs); 
    if min(primal)>=0
        taufeas=tauM; 
        P2val=double(P2);         
        L=(P2val')^(-1)*double(Y); 
    end
else
    yalmiperror(sol.problem) 
end