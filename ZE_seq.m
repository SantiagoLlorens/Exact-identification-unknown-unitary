global D d n
d=2; 
n=4;
D=d^n;

%definition of tester normalization conditions for arbitrary n
%(convention out x in)
MC=cell(d^(4*n-1),1);

% conditions that 1xR is a tensor with the identity (n-th degree of the condition is different from the rest)
it2=0;
for k=1:(d^(2*n-1))
    for l=1:(d^(2*n-1))
        it2=it2+1;
        MC{it2}=sparse(d^(2*n-1)+l,k,1,D^2,D^2);
    end;
end;

for k=1:(d^(2*n-1))
    for l=k:(d^(2*n-1))                            
        it2=it2+1;
        MC{it2}=sparse([l,d^(2*n-1)+l],[k,d^(2*n-1)+k],[1,-1],D^2,D^2);                       
    end;
end;

for a=1:(n-1)
    for k=1:(d^(2*a-1))
        for l=1:(d^(2*a-1))
                pomv1=zeros(2^(n-a),1);
                pomv2=zeros(2^(n-a),1);
                pomv3=zeros(2^(n-a),1);
                for b=1:(2^(n-a))
                    pomv1(b)=odskok(b,a)+d^(2*a-1)+l;
                    pomv2(b)=odskok(b,a)+k;
                    pomv3(b)=1;    
                end;
                it2=it2+1;
                MC{it2}=sparse(pomv1,pomv2,pomv3,D^2,D^2); 
        end;
    end;

    for k=1:(d^(2*a-1))
        for l=k:(d^(2*a-1))                    
                pomv1=zeros(2^(n-a+1),1);
                pomv2=zeros(2^(n-a+1),1);
                pomv3=zeros(2^(n-a+1),1);
                for b=1:(2^(n-a))
                    pomv1(b)=odskok(b,a)+l;
                    pomv2(b)=odskok(b,a)+k;
                    pomv3(b)=1;    
                    
                    pomv1(2^(n-a)+b)=d^(2*a-1)+odskok(b,a)+l;
                    pomv2(2^(n-a)+b)=d^(2*a-1)+odskok(b,a)+k;
                    pomv3(2^(n-a)+b)=-1;                        
                end;                      
                it2=it2+1;                
                MC{it2}=sparse(pomv1,pomv2,pomv3,D^2,D^2);                 
        end;
    end;
end;
MCsize=it2

max=[];
for i=1:d
    for j=1:d
        if j==i
        max(i,j)=1;
        end
    end
end
max=reshape(max,[],1);
keti=kron(max,transpose(max));%identity choi operator

PS=SymmetricProjection(d,2);
dsym=trace(PS);
PA=eye(d^2)-PS;
danti=trace(PA);
PSPS=PermuteSystems(kron(PS,PS),[1,3,2,4],[d,d,d,d]);
PAPA=PermuteSystems(kron(PA,PA),[1,3,2,4],[d,d,d,d]);

A1=kron(keti,kron(keti,(PSPS/dsym+PAPA/danti))); % out in out in out in out in 
A2=PermuteSystems(A1,[7,8,1,2,3,4,5,6]);
A3=PermuteSystems(A1,[5,6,7,8,1,2,3,4]);
A4=PermuteSystems(A1,[3,4,5,6,7,8,1,2]);
A5=PermuteSystems(A4,[3,4,1,2,5,6,7,8]);
A6=PermuteSystems(A4,[1,2,3,4,7,8,5,6]);

cvx_clear
cvx_solver mosek

    cvx_begin SDP
        variable T1(D^2,D^2) semidefinite;
        variable T2(D^2,D^2) semidefinite;
        variable T3(D^2,D^2) semidefinite;
        variable T4(D^2,D^2) semidefinite;
        variable T5(D^2,D^2) semidefinite;
        variable T6(D^2,D^2) semidefinite;
        variable R(D^2,D^2) semidefinite;

        maximize (real(trace(A1*T1)+trace(A2*T2)+trace(A3*T3)+trace(A4*T4)+trace(A5*T5)+trace(A6*T6))/6)

        subject to
            T1+T2+T3+T4+T5+T6<=R
            
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
            
            trace(T3*A1)==0
            trace(T3*A2)==0
            trace(T3*A4)==0
            trace(T3*A5)==0
            trace(T3*A6)==0

            trace(T4*A1)==0
            trace(T4*A2)==0
            trace(T4*A3)==0
            trace(T4*A5)==0
            trace(T4*A6)==0

            trace(T5*A1)==0
            trace(T5*A2)==0
            trace(T5*A3)==0
            trace(T5*A4)==0
            trace(T5*A6)==0

            trace(T6*A1)==0
            trace(T6*A2)==0
            trace(T6*A3)==0
            trace(T6*A4)==0
            trace(T6*A5)==0
            
            trace(R)==D

            for k = 1:MCsize
                trace(MC{k}*R)==0;
            end;
           
    cvx_end;
