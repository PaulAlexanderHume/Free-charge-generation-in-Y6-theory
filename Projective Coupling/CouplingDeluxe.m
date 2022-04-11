clear all
close all

%% Input Section
cd tda1
ReferenceState1=1; ReferenceState2=4; 
NStates=20; dS=1:NStates;
LStates=14; dL=1:LStates;
Ref=importdata('Ref_norm_allstates.mat');
Ed_eV=importdata('E_ES.txt'); Ed=(Ed_eV/27.2114); 
Sr=importdata('S.dat'); C=importdata('C.txt')'; 
NE=importdata('number_electrons.dat');
nbf=importdata('number_basis_functions.dat');
CI=importdata('CI.mat');
CI=CI*sqrt(2);
nso=length(CI(:,:,1));   C=C(:,1:nso);
Ref1=Ref(:,ReferenceState1);  Ref2=Ref(:,ReferenceState2); 
CI1=importdata('CI1raw.mat'); CI1=CI1*sqrt(2);  
cd ..

cd r1
S1=importdata('S.dat'); C1=importdata('C.txt')'; NE1=importdata('number_electrons.dat'); 
nso1=length(CI1(:,:,1));  C1=C1(:,1:nso1); 
cd ..

%% GS Coupling Section
cd tda1
GS=importdata('D0_Check.mat')^2;
RefGSDimerES=importdata('RefGSDimerES_Check.mat');
RefESDimerGS=importdata('RefESDimerGS_Check.mat');
RefGS=[GS;RefGSDimerES];
RefES=[RefESDimerGS;Ref];                                                    
Ref=[RefGS RefES]; Ed=[0;Ed];
ReferenceState1=ReferenceState1+1; ReferenceState2=ReferenceState2+1;
Ref1=Ref(:,ReferenceState1);       Ref2=Ref(:,ReferenceState2);
dS=(1:NStates+1);                  dL=(1:LStates+1);
CInew =zeros(length(CI), length(CI), NStates+1);
CI1new=zeros(length(CI1),length(CI1),LStates+1);
for k=dS
    if k==1
       CInew(:,:,k) =zeros(length(CI) ,length(CI));
    else
       CInew(:,:,k) =CI(:,:,k-1);
    end
end
for k=dL
    if k==1
       CI1new(:,:,k)=zeros(length(CI1),length(CI1));
    else
       CI1new(:,:,k)=CI1(:,:,k-1);
    end
end
CI=CInew; CI1=CI1new;

%% Diabatisation of Reference States
P11R=C1*CI1(:,:,ReferenceState1)*CI1(:,:,ReferenceState1)'*C1';    
P12R=C1*CI1(:,:,ReferenceState1)*CI1(:,:,ReferenceState2)'*C1';
P21R=C1*CI1(:,:,ReferenceState2)*CI1(:,:,ReferenceState1)'*C1';
P22R=C1*CI1(:,:,ReferenceState2)*CI1(:,:,ReferenceState2)'*C1';
PS11R=P11R*S1;         PS12R=P12R*S1;  
PS21R=P21R*S1;         PS22R=P22R*S1;
A11R=diag(PS11R);      A12R=diag(PS12R);
A21R=diag(PS21R);      A22R=diag(PS22R);
x11_DR = sum(A11R(1:(nbf/2)));  x11_AR = sum(A11R(((nbf/2)+1):nbf));
dx11R  = x11_DR-x11_AR;
x22_DR = sum(A22R(1:(nbf/2)));  x22_AR = sum(A22R(((nbf/2)+1):nbf));
dx22R  = x22_DR-x22_AR;
dx12R =(sum(A12R(1:(nbf/2))) - sum(A12R(((nbf/2)+1):nbf)) + sum(A21R(1:(nbf/2)))- sum(A21R(((nbf/2)+1):nbf)))/2;

wR=atan(2*dx12R/(dx11R-dx22R))/2;
Ref1= cos(wR)*Ref(:,ReferenceState1)+sin(wR)*Ref(:,ReferenceState2);
Ref2=-sin(wR)*Ref(:,ReferenceState1)+cos(wR)*Ref(:,ReferenceState2);
Ref(:,ReferenceState1)=Ref1; Ref(:,ReferenceState2)=Ref2;
% overlap=c(1,1)*c(1,2) + c(1,2)*c(2,2)
% save overlap overlap
[a b]=max(max(abs(CI1(:,:,ReferenceState1))));
[c d]=max(abs(CI1(:,b,ReferenceState1)));
Dominant_Excitation_Ref_State1=[d b CI1(d,b,ReferenceState1)]

[a b]=max(max(abs(CI1(:,:,ReferenceState2))));
[c d]=max(abs(CI1(:,b,ReferenceState2)));
Dominant_Excitation_Ref_State2=[d b CI1(d,b,ReferenceState2)]

%% Form Excitonic States from Adiabatic States (Symmetric Systems Only)
% Exc1=Ref(:,ReferenceState1);   
% Exc2=Ref(:,ReferenceState2);    
% Ref(:,ReferenceState1)=1/sqrt(2)*(Exc1+Exc2);
% Ref(:,ReferenceState2)=1/sqrt(2)*(Exc1-Exc2);
% Ref1=1/sqrt(2)*(Exc1+Exc2);
% Ref2=1/sqrt(2)*(Exc1-Exc2);

%% Two State Coupling
%   Overlap correction based on assumption of only two interacting states
S12=(Ref1'*Ref2)/(norm(Ref1)*norm(Ref2));

%   Two-State Electronic Coupling
H11=sum((Ref1(dS)).*Ed(dS).*(Ref1(dS))/(norm(Ref1)*norm(Ref1)));
H12=sum((Ref1(dS)).*Ed(dS).*(Ref2(dS))/(norm(Ref1)*norm(Ref2)));
H22=sum((Ref2(dS)).*Ed(dS).*(Ref2(dS))/(norm(Ref2)*norm(Ref2)));
H=[H11 H12;H12 H22];

E1=0.5*(1/(1-S12^2))*(H11+H22-2*H12*S12*(H11-H22)*sqrt(1-S12^2));
V=(H12-1/2*(H11+H22)*S12)/(1-S12^2);   
E2=0.5*(1/(1-S12^2))*(H11+H22-2*H12*S12*(H22-H11)*sqrt(1-S12^2));
E_Check=(Ed(2)-Ed(1))/2;

E1eV=27.2114*H11;
VeV=27.2114*V
E2eV=27.2114*H22;

%% Coupling Calculated Using Lowdin Orthogonalisation
S=zeros(LStates,LStates);
H=zeros(LStates,LStates);
for i=1:LStates
    for j=1:LStates
        S(i,j)=Ref(:,i)'*Ref(:,j)/(norm(Ref(:,i))*norm(Ref(:,j)));
        H(i,j)=sum((Ref(dS,i)).*Ed(dS).*(Ref(dS,j))/(norm(Ref(dS,i))*norm(Ref(dS,j))));
    end
end
Sinv=sqrtm(inv(S));

%sinv=inv(sqrtm(s));
[U,D]=eig(S); D12=sqrt(inv(D)); %D12(1,:)=-D12(1,:);
Sinv=U*D12*U';

% I=eye(size(s));
% sinv=I-1/2*(s-I)+3/8*(s-I)^2-5/16*(s-I)^3;


Heff=Sinv*H*Sinv; Heff_eV=Heff.*27.2114;
Heff(ReferenceState1,ReferenceState2)*27.2114
Heff(ReferenceState1,ReferenceState1)*27.2114;
Heff(ReferenceState2,ReferenceState2)*27.2114;
H_eV=27.2114*H;

%% Print Key Couplings - Heff_eV
Ex1=1;     Ex2=4;     CT1=2;     CT2=7;    G=0;

Ex1=Ex1+1; Ex2=Ex2+1; CT1=CT1+1; CT2=CT2+1; G=1;

E_Ex1_Ex2_CT1_CT2=[Heff_eV(Ex1,Ex1) Heff_eV(Ex2,Ex2) Heff_eV(CT1,CT1) Heff_eV(CT2,CT2)]'

H_Ex1=[Heff_eV(Ex1,Ex2) Heff_eV(Ex1,CT1) Heff_eV(Ex1,CT2) Heff_eV(Ex1,G)]';

H_Ex2=[Heff_eV(Ex2,Ex1) Heff_eV(Ex2,CT1) Heff_eV(Ex2,CT2) Heff_eV(Ex2,G)]';

H_CT1=[Heff_eV(CT1,Ex1) Heff_eV(CT1,Ex2) Heff_eV(CT1,CT2) Heff_eV(CT1,G)]';

H_CT2=[Heff_eV(CT2,Ex1) Heff_eV(CT2,Ex2) Heff_eV(CT2,CT1) Heff_eV(CT2,G)]';

tda_table=[E_Ex1_Ex2_CT1_CT2 [H_Ex1 H_Ex2 H_CT1 H_CT2]'];

save E_Ex1_Ex2_CT1_CT2 E_Ex1_Ex2_CT1_CT2 
save tda_table tda_table

cd ..

%% FED Coupling
% V=0;m=2; n=6;
% P11=C*CI(:,:,m)*CI(:,:,m)'*C';    
% P12=C*CI(:,:,m)*CI(:,:,n)'*C';
% P21=C*CI(:,:,n)*CI(:,:,m)'*C';
% P22=C*CI(:,:,n)*CI(:,:,n)'*C';
% PS11=P11*Sr;         PS12=P12*Sr;  
% PS21=P21*Sr;         PS22=P22*Sr;
% A11=diag(PS11);     A12=diag(PS12);
% A21=diag(PS21);     A22=diag(PS22);
% x11_D = sum(A11(1:(nbf/2)));  x11_A = sum(A11(((nbf/2)+1):nbf));
% dx11  = x11_D-x11_A;
% x22_D = sum(A22(1:(nbf/2)));  x22_A = sum(A22(((nbf/2)+1):nbf));
% dx22  = x22_D-x22_A;
% dx12 =(sum(A12(1:(nbf/2))) - sum(A12(((nbf/2)+1):nbf)) + sum(A21(1:(nbf/2)))- sum(A21(((nbf/2)+1):nbf)))/2;
% 
% % w=atan(2*dx12/(dx11-dx22))/2
% % 
% % overlap=c(1,1)*c(1,2) + c(1,2)*c(2,2)
% % save overlap overlap
% 
% V=(Ed(n)-Ed(m))*abs(dx12)/sqrt(((dx11-dx22)^2)+4*dx12^2)*27.2114
