% 20-05-14 08:23 Added By Xingyu Zhu
% This file calculates the current/velocity of each component


% % -- main program begin to set the matrix elements --
k=sqrt(kz^2+kx^2);
bsab=kx*rhocsab;
bsab(abs(bsab)<1e-50)=1e-50;  % to avoid singular when k_perp=0

bsab2=bsab.^2;
M=sparse(NN,NN);

% initialize
b11_s=zeros(S); b12_s=zeros(S); b13_s=zeros(S); 
b21_s=zeros(S); b22_s=zeros(S); b23_s=zeros(S); 
b31_s=zeros(S); b32_s=zeros(S); b33_s=zeros(S);
cnj_s=zeros(3*SNJ,S);
b11nj_s=cnj_s.*0; b12nj_s=cnj_s.*0; b13nj_s=cnj_s.*0;
b21nj_s=cnj_s.*0; b22nj_s=cnj_s.*0; b23nj_s=cnj_s.*0;
b31nj_s=cnj_s.*0; b32nj_s=cnj_s.*0; b33nj_s=cnj_s.*0;

for s=1:S % species
  nj=0;
  for n=-N:N % Bessel function
    %             Gamn=exp(-bs2(s))*besseli(n,bs2(s)); % large k_perp will NaN
    %             Gamnp=exp(-bs2(s))*(besseli(n+1,bs2(s))+...
    %                 besseli(n-1,bs2(s))-2*besseli(n,bs2(s)))/2;
    
    for j=1:J % poles of Z(zeta)
      nj=nj+1;
      
      for iab=1:2 % 2018-10-13 10:58
      
      Gamn=besseli(n,bsab2(iab,s),1); % 2014-10-13 12:51
      Gamnp=(besseli(n+1,bsab2(iab,s),1)+...
        besseli(n-1,bsab2(iab,s),1)-2*besseli(n,bsab2(iab,s),1))/2;
  
      cnj_s(nj,s)=czj(j)*kz*vtzs(s)+kz*vds(s)+n*wcs(s); %

      cnj=cnj_s(nj,s);
      bj0ab=vds(s)+(1-1/lmdTab(iab,s))*czj(j)*vtzs(s);

      %
      if(n==0)  % for A_nj
        bnj1=bj0ab/(czj(j)*vtzs(s)+vds(s)); % avoid cnj=0
      else
        bnj1=kz*bj0ab/cnj;
      end
      bnj2=1-bnj1;

      tmp=wps2(s)*bzj(j);

      b11nj_s(nj,s)=b11nj_s(nj,s)+rsab(iab,s)*tmp*bnj2*n^2*Gamn/bsab2(iab,s);
      b11_s(s)=b11_s(s)+rsab(iab,s)*tmp*bnj1*n^2*Gamn/bsab2(iab,s);

      b12nj_s(nj,s)=b12nj_s(nj,s)+rsab(iab,s)*tmp*bnj2*1i*n*Gamnp;
      b12_s(s)=b12_s(s)+rsab(iab,s)*tmp*bnj1*1i*n*Gamnp;
      b21nj_s(nj,s)=-b12nj_s(nj,s);
      b21_s(s)=-b12_s(s);

      b22nj_s(nj,s)=b22nj_s(nj,s)+rsab(iab,s)*tmp*bnj2*(n^2*Gamn/bsab2(iab,s)...
          -2*bsab2(iab,s)*Gamnp);
      b22_s(s)=b22_s(s)+rsab(iab,s)*tmp*bnj1*(n^2*Gamn/bsab2(iab,s)...
          -2*bsab2(iab,s)*Gamnp);

      %
      if(n==0)  % for eta_n*A_nj
        bnj1=0; % avoid cnj=0 when kz=0
      else
        bnj1=n*wcs(s)*bj0ab/cnj/vtzs(s);
      end
      bnj2=czj(j)/lmdTab(iab,s)+bnj1; %

      b13nj_s(nj,s)=b13nj_s(nj,s)+rsab(iab,s)*tmp*bnj2*n*...
          sqrt(2*lmdTab(iab,s))*Gamn/bsab(iab,s);
      b13_s(s)=b13_s(s)-rsab(iab,s)*tmp*bnj1*n*sqrt(2*lmdTab(iab,s))*Gamn/bsab(iab,s); %
      b31nj_s(nj,s)=b13nj_s(nj,s);
      b31_s(s)=b13_s(s);

      b23nj_s(nj,s)=b23nj_s(nj,s)-rsab(iab,s)*1i*tmp*bnj2*...
          sqrt(2*lmdTab(iab,s))*Gamnp*bsab(iab,s);
      b23_s(s)=b23_s(s)+rsab(iab,s)*1i*tmp*bnj1*sqrt(2*lmdTab(iab,s))*Gamnp*bsab(iab,s); %
      b32nj_s(nj,s)=-b23nj_s(nj,s);
      b32_s(s)=-b23_s(s);

      %
      if(bj0ab==0 || kz==0)  % for eta_n^2*A_nj
        bnj1=0;
        bnj2=czj(j)*czj(j);
      else
        % bnj1=n^2*bj0ab/cnj/vtzs(s)^2/kz;
        bnj1=n^2*wcs(s)^2*bj0ab/cnj/vtzs(s)^2/kz; % !!fixed bug of missed wcs^2, 18-10-01 10:19
        bnj2=(vds(s)/vtzs(s)+czj(j))*czj(j)/lmdTab(iab,s)+...
          n*wcs(s)*bj0ab*(1-n*wcs(s)/cnj)/vtzs(s)^2/kz; %
      end

      b33nj_s(nj,s)=b33nj_s(nj,s)+rsab(iab,s)*tmp*bnj2*2*lmdTab(iab,s)*Gamn;
      b33_s(s)=b33_s(s)+rsab(iab,s)*tmp*bnj1*2*lmdTab(iab,s)*Gamn;
      end
    end
  end
end


wtmp = wws2(jpa,1,jpl)*wcs1;
dEx_tmp=Pola_SI(jpa,1,jpl,1); dEy_tmp=Pola_SI(jpa,1,jpl,2); dEz_tmp=Pola_SI(jpa,1,jpl,3);
dBx_tmp=Pola_SI(jpa,1,jpl,4); dBy_tmp=Pola_SI(jpa,1,jpl,5); dBz_tmp=Pola_SI(jpa,1,jpl,6);
for s=1:S
    Js(jpa,1,s,jpl)=0.0; Js(jpa,2,s,jpl)=0.0; Js(jpa,3,s,jpl)=0.0;
    Jxs=0.0; Jys=0.0; Jzs=0.0;
    nj=0;
    for n=-N:N
        for j=1:J
            nj=nj+1;
            Jxs=Jxs+b11nj_s(nj,s)/(wtmp-cnj_s(nj,s))*dEx_tmp+...
                b12nj_s(nj,s)/(wtmp-cnj_s(nj,s))*dEy_tmp+...
                b13nj_s(nj,s)/(wtmp-cnj_s(nj,s))*dEz_tmp;
            Jys=Jys+b21nj_s(nj,s)/(wtmp-cnj_s(nj,s))*dEx_tmp+...
                b22nj_s(nj,s)/(wtmp-cnj_s(nj,s))*dEy_tmp+...
                b23nj_s(nj,s)/(wtmp-cnj_s(nj,s))*dEz_tmp;
            Jzs=Jzs+b31nj_s(nj,s)/(wtmp-cnj_s(nj,s))*dEx_tmp+...
                b32nj_s(nj,s)/(wtmp-cnj_s(nj,s))*dEy_tmp+...
                b33nj_s(nj,s)/(wtmp-cnj_s(nj,s))*dEz_tmp;
        end
    end
    Js(jpa,1,s,jpl)=-1i*epsilon0*(Jxs+b11_s(s)/wtmp*dEx_tmp+b12_s(s)/wtmp*dEy_tmp+b13_s(s)/wtmp*dEz_tmp);
    Js(jpa,2,s,jpl)=-1i*epsilon0*(Jys+b21_s(s)/wtmp*dEx_tmp+b22_s(s)/wtmp*dEy_tmp+b23_s(s)/wtmp*dEz_tmp);
    Js(jpa,3,s,jpl)=-1i*epsilon0*(Jzs+b31_s(s)/wtmp*dEx_tmp+b32_s(s)/wtmp*dEy_tmp+b33_s(s)/wtmp*dEz_tmp);
    dV(jpa,1,s,jpl)=Js(jpa,1,s,jpl)/(qs(s)*ns0(s));
    dV(jpa,2,s,jpl)=Js(jpa,2,s,jpl)/(qs(s)*ns0(s));
    dV(jpa,3,s,jpl)=Js(jpa,3,s,jpl)/(qs(s)*ns0(s));
    dVnorm(jpa,1,s,jpl)=dV(jpa,1,s,jpl)/vA;
    dVnorm(jpa,2,s,jpl)=dV(jpa,2,s,jpl)/vA;
    dVnorm(jpa,3,s,jpl)=dV(jpa,3,s,jpl)/vA;
    %density fluctuation
    xinorm(jpa,s,jpl) = (dV(jpa,1,s,jpl)*kx+dV(jpa,3,s,jpl)*kz)/wtmp;
    % dissipation
    JE(jpa,1,s,jpl)=-(Js(jpa,1,s,jpl)*conj(dEx_tmp)+conj(Js(jpa,1,s,jpl))*dEx_tmp)/4;
    JE(jpa,2,s,jpl)=-(Js(jpa,2,s,jpl)*conj(dEy_tmp)+conj(Js(jpa,2,s,jpl))*dEy_tmp)/4;
    JE(jpa,3,s,jpl)=-(Js(jpa,3,s,jpl)*conj(dEz_tmp)+conj(Js(jpa,3,s,jpl))*dEz_tmp)/4;
    EBenergy=(dEx_tmp*conj(dEx_tmp)+dEy_tmp*conj(dEy_tmp)+dEz_tmp*conj(dEz_tmp))*epsilon0*0.5+...
        (dBx_tmp*conj(dBx_tmp)+dBy_tmp*conj(dBy_tmp)+dBz_tmp*conj(dBz_tmp))/mu0*0.5;
    JE(jpa,1,s,jpl)=JE(jpa,1,s,jpl)/EBenergy/2;
    JE(jpa,2,s,jpl)=JE(jpa,2,s,jpl)/EBenergy/2;
    JE(jpa,3,s,jpl)=JE(jpa,3,s,jpl)/EBenergy/2;
end
