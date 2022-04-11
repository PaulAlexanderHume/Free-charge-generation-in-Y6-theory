clear all
format long
'Importing Data for Reference State'
cd r1
C1=importdata('C.txt'); ES_ref=importdata('ES.mat');
nbf1=importdata('number_basis_functions.dat');                           
NE1=importdata('number_electrons.dat');                                 
C1=C1'; %C1=[C1 zeros(length(C1),1)];
[sizeC1a sizeC1b]=size(C1); 							  
Occupied1=1:NE1/2;                                                     
NRefStates=20; 
cd ..

'Importing Data from Dimer'
cd tda1
Cd=importdata('C.txt'); Sd=importdata('S.dat'); 
nbf=importdata('number_basis_functions.dat');
NE=importdata('number_electrons.dat');                             
ES=importdata('ES.mat'); 				                                                                                                         
Cd=Cd'; %Cd=[Cd zeros(length(Cd),2)];
[sizeCda sizeCdb]=size(Cd);								 
OccupiedDimer=1:NE/2;                                              
NStates=20; 

nso=min(sizeC1b,sizeCdb);					 			 %new, number of spin orbitals
Unoccupied1=(NE1/2+1):nso; 					 			 
UnoccupiedDimer=(NE/2+1):nso;				 			 
C1=C1(:,1:nso);								 			 
Cd=Cd(:,1:nso);								 			 

'Constructing CI Matrices'
CI1=zeros(nso,nso,NRefStates);	 						 
for j=1:NRefStates
    for i=1:length(ES_ref)
        if ES_ref(i,1,j)~=0 
        CI1(ES_ref(i,1,j),ES_ref(i,2,j),j)=ES_ref(i,3,j);                            
        end 
    end
end
save CI1raw CI1

CI=zeros(nso,nso,NStates);								 
for k=1:NStates
    for i=1:length(ES)
        if ES(i,1,k)~=0 
        CI(ES(i,1,k),ES(i,2,k),k)=ES(i,3,k);                            
        end
    end
end
save CI CI

'Defining Which States Contribute'
for j=1:NRefStates
   Occupied1_nonzero_temp=find(~all(CI1(:,:,j)==0,2));
   a1(j)=min(Occupied1_nonzero_temp);
   b1(j)=max(Occupied1_nonzero_temp);
   Unoccupied1_nonzero_temp=find(~all(CI1(:,:,j)==0,1));
   c1(j)=max(Unoccupied1_nonzero_temp);
end
Occupied1_nonzero=min(a1):max(b1);          
Unoccupied1=(NE/2+1):max(c1);
C1=C1(:,[Occupied1 Unoccupied1]);		
CI1=CI1([Occupied1 Unoccupied1],[Occupied1 Unoccupied1],:);	

for k=1:NStates
   OccupiedDimer_nonzero_temp=find(~all(CI(:,:,k)==0,2));
   ad(k)=min(OccupiedDimer_nonzero_temp);
   bd(k)=max(OccupiedDimer_nonzero_temp);
   UnoccupiedDimer_nonzero_temp=find(~all(CI(:,:,k)==0,1));
   c1(k)=max(UnoccupiedDimer_nonzero_temp);
end

OccupiedDimer_nonzero=min(ad):max(bd);
UnoccupiedDimer=(NE/2+1):max(c1);
Cd=Cd(:,[OccupiedDimer UnoccupiedDimer]);	
CI=CI([OccupiedDimer UnoccupiedDimer],[OccupiedDimer UnoccupiedDimer],:);		

'Constructing Dimer Overlap Matrices'
S=C1'*Sd*Cd;                                                                                  
S0=S(Occupied1,OccupiedDimer);                                                                                                       
D0=det(S0); 

'Renormalising Determinants'
N=zeros(NE1/2,nbf1); Sm=C1'*Sd*C1;
Sm0=Sm(Occupied1,Occupied1);
Dm0=det(Sm0);
for a=Occupied1
    for r=Unoccupied1
        Dm=Sm0;
        Dm(a,:)=Sm(r,Occupied1);
        Dm(:,a)=Sm(Occupied1,r);
        Dm(a,a)=Sm(r,r);
        N(a,r)=sqrt(det(Dm));
    end
end    
save N N
N=importdata('N.mat');                                                       

'Pre-Calculating Row Sums'
X=zeros(max(Unoccupied1),max(UnoccupiedDimer),NRefStates);
for j=1:NRefStates
    CI1(Occupied1,Unoccupied1,j)=CI1(Occupied1,Unoccupied1,j)./N(:,Unoccupied1);
    X(:,:,j)=CI1(:,:,j)*S;
end
save RefCI1 CI1

'Pre-Calculating Column Sums'
Y=zeros(max(Unoccupied1),max(UnoccupiedDimer),NStates);
for k=1:NStates
    Y(:,:,k)=S*CI(:,:,k)';
end

%'Pre-Calculating Cross-Terms'
%Z=zeros(nbf,nbf,NRefStates,NStates);
%for j=1:NRefStates
%   for k=1:NStates
%       Z(:,:,j,k)=CI1(:,:,j)*S*CI(:,:,k)';
%    end
%end		

%'Calculating Cofactor Matrices'
%Cf=zeros(NE/2,NE/2,NE/2,NRefStates);                    %Cf is the cofactor matrix
%for j=1:NRefStates
%    for a=Occupied1
%        M=S0;									                          %M is the working matrix
%        M(a,:)=X(a,OccupiedDimer,j);                    %Replacing row a with sum over r
%        for b=OccupiedDimer
%            for i=Occupied1
%                redM=M;
%                redM(i,:)=[];
%                redM(:,b)=[];
%                Cf(i,b,a,j)=((-1)^(i+b))*det(redM);    %i and b are the row/column you delete. a is the row X is inserted. j is the reference state
%            end
%        end
%    end
%end

'Calculating Overlaps of Monomer-Dimer Electronic States'
S_ES=zeros(NStates,NRefStates);

for j=1:NRefStates
    for k=1:NStates
    if S_ES(k,j)==0      
        D=zeros(NE/2,NE/2);
        Zjk=CI1(:,:,j)*S*CI(:,:,k)';
        for a=Occupied1
            for b=OccupiedDimer
             %Direct Algorithm           
               M=S0;									                 %M is the working matrix
               M(a,:)=X(a,OccupiedDimer,j);            %Replacing row a with sum over r
               M(:,b)=Y(Occupied1,b,k);			           %Replacing column b with sum over s
               M(a,b)=Zjk(a,b);						             %Replacing entry (a,b) with sum over r and s
               D(a,b)=2*det(M)*(D0/sqrt(Dm0));			     %Calculating determinant, multiply by 2 because cofactor normalisation. Multiply by beta spin block
               
              %Cofactor Algorithm
%              v=Y(Occupied1,b,k);                   %Replacing column b with sum over s
%              v(a)=Zjk(a,b);						            %Replacing entry (a,b) with sum over r and s
%              D(a,b)=2*Cf(:,b,a,j)'*v*(D0/sqrt(Dm0));
            end										    
        end
        S_ES(k,j)=sum(sum(D));							%Sum terms to give projection of states.
        save Ref_norm_allstates S_ES 
    end
    end
end
cd ..

