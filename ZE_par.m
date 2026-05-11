addpath(genpath('./../QETLAB-1.0'))

global D d n
d=2; 
n=4; 
D=d^n;

max=[];
for i=1:d
    for j=1:d
        if j==i
        max(i,j)=1;
        end
    end
end
max=reshape(max,[],1);
keti=kron(max,transpose(max)); %identity choi operator

PS=SymmetricProjection(d,2);
dsym=trace(PS);
PA=eye(d^2)-PS;
danti=trace(PA);
PSPS=PermuteSystems(kron(PS,PS),[1,3,2,4],[d,d,d,d]);
PAPA=PermuteSystems(kron(PA,PA),[1,3,2,4],[d,d,d,d]);

A1=kron(keti,kron(keti,(PSPS/dsym+PAPA/danti)));
A2=PermuteSystems(A1,[7,8,1,2,3,4,5,6]);
A3=PermuteSystems(A1,[5,6,7,8,1,2,3,4]);
A4=PermuteSystems(A1,[3,4,5,6,7,8,1,2]);
A5=PermuteSystems(A4,[3,4,1,2,5,6,7,8]);
A6=PermuteSystems(A4,[1,2,3,4,7,8,5,6]);

A1=PermuteSystems(A1,[1,3,5,7,2,4,6,8]); % out out out out in in in in 
A2=PermuteSystems(A2,[1,3,5,7,2,4,6,8]);
A3=PermuteSystems(A3,[1,3,5,7,2,4,6,8]);
A4=PermuteSystems(A4,[1,3,5,7,2,4,6,8]);
A5=PermuteSystems(A5,[1,3,5,7,2,4,6,8]);
A6=PermuteSystems(A6,[1,3,5,7,2,4,6,8]);

A1=sparse(A1);
A2=sparse(A2);
A3=sparse(A3);
A4=sparse(A4);
A5=sparse(A5);
A6=sparse(A6);

cvx_clear
cvx_solver mosek

    cvx_begin SDP
        variable T1(D^2,D^2) semidefinite;
        variable T2(D^2,D^2) semidefinite;
        variable T3(D^2,D^2) semidefinite;
        variable T4(D^2,D^2) semidefinite;
        variable T5(D^2,D^2) semidefinite;
        variable T6(D^2,D^2) semidefinite;
        variable r(D,D) semidefinite;

        maximize (1/6*real(trace(A1*T1)+trace(A2*T2)+trace(A3*T3)+trace(A4*T4)+trace(A5*T5)+trace(A6*T6)))

        subject to
            -T1-T2-T3-T4-T5-T6+kron(eye(D),r) ==  semidefinite(D^2)
            
            trace(T1*A2)==0
            trace(T1*A3)==0
            trace(T1*A4)==0
            trace(T1*A5)==0
            trace(T1*A6)==0

            trace(T2*A1)==0
            trace(T2*A3)==0
            trace(T2*A4)==0
            trace(T2*A5)==0
            trace(T2*A6)==0

            trace(T3*A2)==0
            trace(T3*A1)==0
            trace(T3*A4)==0
            trace(T3*A5)==0
            trace(T3*A6)==0

            trace(T4*A2)==0
            trace(T4*A3)==0
            trace(T4*A1)==0
            trace(T4*A5)==0
            trace(T4*A6)==0

            trace(T5*A2)==0
            trace(T5*A3)==0
            trace(T5*A4)==0
            trace(T5*A1)==0
            trace(T5*A6)==0

            trace(T6*A2)==0
            trace(T6*A3)==0
            trace(T6*A4)==0
            trace(T6*A5)==0
            trace(T6*A1)==0

            trace(r)==1
           
    cvx_end;