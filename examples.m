a=25;       % reaction coefficient 
delta=1;	% decay rate 
%% Injection gain 
N=ceil(sqrt(delta+a)/pi-1/2);           % N is such that (16) is true 
lambda=((2*(1:N)-1)*pi/2).^2;           % =(lambda_1,...,lambda_N)
A=diag(-lambda+a); C=sqrt(2)*ones(1,N); % parameters from (23)
%% Observer with delayed measurements (Theorem 1 and Remark 4)
epsilon=.5; % tuning parameter for the design
tau0=.02;   % initial approximation of tauM
if LMI_TAC18a_rem4(A,C,tau0,delta,epsilon)==0 
    disp('Remark 4: the LMIs are not feasible'); 
else
    tauM=fminsearch(@(x) -LMI_TAC18a_rem4(A,C,x,delta,epsilon),tau0);   % maximum tauM so that the design LMIs are feasible
    [~,L]=LMI_TAC18a_rem4(A,C,tauM,delta,epsilon);                      % calculates the injection gain using the design LMIs 
    tauM=fminsearch(@(x) -LMI_TAC18a_th1(A,C,L,x,delta),tauM);          % solve the LMIs of Theorem 1 to find the maximum tauM
    disp(['tauM=' num2str(tauM)]); 
end
%% Observer with sampled in time measurements (Theorem 1 and Remark 4)
epsilon=.5; % tuning parameter for the design
h0=.04;     % initial approximation of h
if LMI_TAC18a_rem6(A,C,h0,delta,epsilon)==0
    disp('Remark 6: the LMIs are not feasible'); 
else
    h=fminsearch(@(x) -LMI_TAC18a_rem6(A,C,x,delta,epsilon),h0);    % maximum h so that the design LMIs are feasible
    [~,L]=LMI_TAC18a_rem6(A,C,h,delta,epsilon);                     % calculates the injection gain using the design LMIs 
    h=fminsearch(@(x) -LMI_TAC18a_th2(A,C,L,h,delta),h);            % solve the LMIs of Theorem 2 to find the maximum h 
    disp(['h=' num2str(h)]); 
end