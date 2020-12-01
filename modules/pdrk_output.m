% 18-10-06 07:05 Hua-sheng XIE, huashengxie@gmail.com, FRI-ENN, China
% Ackn.: Richard Denton (Dartmouth), Xin Tao (USTC), Jin-song Zhao (PMO), 
% etc ...
% Use this file to output the results.

% iout=2;
% npl1=4;
npl1=npl;
if(iout==2)
  % calculate the polarizations for wws
  sp=1; nw=1;
  Pola=zeros(npa,npb,npl,20);
  wws2=wws;
  icalp=1; % tell kernel.m which can calculate polarization now
  for jpl=1:npl1
    % in eigs(), use wws(jpa,jpb,jpl) as intitial guess for each (pa,pb)
    run ../modules/pdrk_kernel; % to update
  end
end

%% added by Xingyu Zhu
%       FileName = '/output/velocity.dat';
%       matrix=[pas',real(wws(:,1,1)),imag(wws(:,1,1)),real(Pola_norm(:,1,1,1)),imag(Pola_norm(:,1,1,1)),...
%           real(Pola_norm(:,1,1,2)),imag(Pola_norm(:,1,1,2)),real(Pola_norm(:,1,1,3)),imag(Pola_norm(:,1,1,3)),...
%           real(Pola_norm(:,1,1,4)),imag(Pola_norm(:,1,1,4)),real(Pola_norm(:,1,1,5)),imag(Pola_norm(:,1,1,5)),...
%           real(Pola_norm(:,1,1,6)),imag(Pola_norm(:,1,1,6)),real(dVnorm(:,1,3,1)),imag(dVnorm(:,1,3,1)),...
%           real(dVnorm(:,2,3,1)),imag(dVnorm(:,2,3,1)),real(dVnorm(:,3,3,1)),imag(dVnorm(:,3,3,1)),...
%           real(dVnorm(:,1,1,1)),imag(dVnorm(:,1,1,1)),real(dVnorm(:,2,1,1)),imag(dVnorm(:,2,1,1)),...
%           real(dVnorm(:,3,1,1)),imag(dVnorm(:,3,1,1)),real(dVnorm(:,1,2,1)),imag(dVnorm(:,1,2,1)),...
%           real(dVnorm(:,2,2,1)),imag(dVnorm(:,2,2,1)),real(dVnorm(:,3,2,1)),imag(dVnorm(:,3,2,1))];
%       save(FileName,'matrix','-ascii');


% % To do: add group velocity, etc, 18-10-18 15:00
%%
%close all;
if(iout==2)
% h=figure('unit','normalized','Position',[0.01 0.25 0.5 0.6],...
%   'DefaultAxesFontSize',15);
for jpl=1:npl1
    if(ipa==ipb) % plot 1D polarizations, to update. 2018-10-19 18:10
        if iloga == 1
            pas = 10.^pas;
            pa = pas;
        end
      h1 = figure('unit','normalized','Position',[0.1 0.1 0.8 0.8],...
      'DefaultAxesFontSize',8);;
      subplot(451); hold on; box on;
      plot(pas,real(wws2(:,1,jpl)),'-','Color',pltc(jpl,:),'linewidth',2);
      xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
      xlabel(strpa);ylabel('\omega_r/\omega_{c1}');
      
      subplot(452); hold on; box on;
      plot(pas,imag(wws2(:,1,jpl)),'-','Color',pltc(jpl,:),'linewidth',2);
      xlim([min(pa),max(pa)]);
      xlabel(strpa);ylabel('\omega_i/\omega_{c1}');
      
      subplot(453); hold on; box on;
      plot(pas,Pola_SI(:,1,jpl,7)./Pola_SI(:,1,jpl,8),'-','Color',...
          pltc(jpl,:),'linewidth',2);
      xlim([min(pa),max(pa)]);
      xlabel(strpa);ylabel('Energy E/Energy B');
      
      subplot(456); hold on; box on;
      plot(pas,real(Pola_norm(:,1,jpl,2)./(1i*Pola_norm(:,1,jpl,1))),'-','Color',...
          pltc(jpl,:),'linewidth',2);
      plot(pas,real(Pola_norm(:,1,jpl,2)./(1i*Pola_norm(:,1,jpl,1))),'-','Color',...
          pltc(jpl,:),'linewidth',2);
      xlim([min(pa),max(pa)]);
      xlabel(strpa);ylabel('E_y/(iE_x)');
      
      subplot(457); hold on; box on;
      plot(pas,abs(Pola_norm(:,1,jpl,3))./abs(Pola_norm(:,1,jpl,1)),'-','Color',...
          pltc(jpl,:),'linewidth',2);
      plot(pas,abs(Pola_norm(:,1,jpl,3))./abs(Pola_norm(:,1,jpl,1)),'-','Color',...
          pltc(jpl,:),'linewidth',2);
      xlim([min(pa),max(pa)]);
      xlabel(strpa);ylabel('|E_z|/|E_x|');
      
      subplot(458); hold on; box on;
      plot(pas,real(Pola_norm(:,1,jpl,1)),'-','Color',pltc(jpl,:),'linewidth',2);
      hold on;
      plot(pas,imag(Pola_norm(:,1,jpl,1)),'--','Color',pltc(jpl,:),'linewidth',2);
      xlim([min(pa),max(pa)]);
      xlabel(strpa);ylabel('E_x/B_0V_A');
      legend('Re','Im'); legend('boxoff');
      
      subplot(459); hold on; box on;
      plot(pas,real(Pola_norm(:,1,jpl,2)),'-','Color',pltc(jpl,:),'linewidth',2);
      hold on;
      plot(pas,imag(Pola_norm(:,1,jpl,2)),'--','Color',pltc(jpl,:),'linewidth',2);
      xlim([min(pa),max(pa)]);
      xlabel(strpa);ylabel('E_y/B_0V_A');
      
      subplot(4,5,10); hold on; box on;
      plot(pas,real(Pola_norm(:,1,jpl,3)),'-','Color',pltc(jpl,:),'linewidth',2);
      hold on;
      plot(pas,imag(Pola_norm(:,1,jpl,3)),'--','Color',pltc(jpl,:),'linewidth',2);
      xlim([min(pa),max(pa)]);
      xlabel(strpa);ylabel('E_z/B_0V_A');
      
      subplot(4,5,11); hold on; box on;
      plot(pas,real(Pola_norm(:,1,jpl,5)./(1i*Pola_norm(:,1,jpl,4))),'-','Color',...
          pltc(jpl,:),'linewidth',2);
      plot(pas,real(Pola_norm(:,1,jpl,5)./(1i*Pola_norm(:,1,jpl,4))),'-','Color',...
          pltc(jpl,:),'linewidth',2);
      xlim([min(pa),max(pa)]);
      xlabel(strpa);ylabel('B_y/(iB_x)');
      
      subplot(4,5,12); hold on; box on;
      plot(pas,abs(Pola_norm(:,1,jpl,6))./abs(Pola_norm(:,1,jpl,4)),'-','Color',...
          pltc(jpl,:),'linewidth',2);
      plot(pas,abs(Pola_norm(:,1,jpl,6))./abs(Pola_norm(:,1,jpl,4)),'-','Color',...
          pltc(jpl,:),'linewidth',2);
      xlim([min(pa),max(pa)]);
      xlabel(strpa);ylabel('|B_z|/|B_x|');
      ytickformat('%.1f');
      
      subplot(4,5,13); hold on; box on;
      plot(pas,real(Pola_norm(:,1,jpl,4)),'-','Color',pltc(jpl,:),'linewidth',2);
      hold on;
      plot(pas,imag(Pola_norm(:,1,jpl,4)),'--','Color',pltc(jpl,:),'linewidth',2);
      xlim([min(pa),max(pa)]);
      xlabel(strpa);ylabel('B_x/B_0');
      legend('Re','Im'); legend('boxoff');
      
      subplot(4,5,14); hold on; box on;
      plot(pas,real(Pola_norm(:,1,jpl,5)),'-','Color',pltc(jpl,:),'linewidth',2);
      hold on;
      plot(pas,imag(Pola_norm(:,1,jpl,5)),'--','Color',pltc(jpl,:),'linewidth',2);
      xlim([min(pa),max(pa)]);
      xlabel(strpa);ylabel('B_y/B_0');
      
      subplot(4,5,15); hold on; box on;
      plot(pas,real(Pola_norm(:,1,jpl,6)),'-','Color',pltc(jpl,:),'linewidth',2);
      hold on;
      plot(pas,imag(Pola_norm(:,1,jpl,6)),'--','Color',pltc(jpl,:),'linewidth',2);
      xlim([min(pa),max(pa)]);
      xlabel(strpa);ylabel('B_z/B_0');
      
      subplot(4,5,16); hold on; box on;
      plot(pas,real(Pola_norm(:,1,jpl,10)./(1i*Pola_norm(:,1,jpl,9))),'-','Color',...
          pltc(jpl,:),'linewidth',2);
      plot(pas,real(Pola_norm(:,1,jpl,10)./(1i*Pola_norm(:,1,jpl,9))),'-','Color',...
          pltc(jpl,:),'linewidth',2);
      xlim([min(pa),max(pa)]);
      xlabel(strpa);ylabel('J_y/(iJ_x)');
      
      subplot(4,5,17); hold on; box on;
      plot(pas,abs(Pola_norm(:,1,jpl,11))./abs(Pola_norm(:,1,jpl,9)),'-','Color',...
          pltc(jpl,:),'linewidth',2);
      plot(pas,abs(Pola_norm(:,1,jpl,11))./abs(Pola_norm(:,1,jpl,9)),'-','Color',...
          pltc(jpl,:),'linewidth',2);
      xlim([min(pa),max(pa)]);
      xlabel(strpa);ylabel('|J_z|/|J_x|');
      ytickformat('%.1f');
      
      subplot(4,5,18); hold on; box on;
      plot(pas,real(Pola_norm(:,1,jpl,9)),'-','Color',pltc(jpl,:),'linewidth',2);
      hold on;
      plot(pas,imag(Pola_norm(:,1,jpl,9)),'--','Color',pltc(jpl,:),'linewidth',2);
      xlim([min(pa),max(pa)]);
      xlabel(strpa);ylabel('J_x/(q_in_iV_A)');
      legend('Re','Im'); legend('boxoff');
      
      subplot(4,5,19); hold on; box on;
      plot(pas,real(Pola_norm(:,1,jpl,10)),'-','Color',pltc(jpl,:),'linewidth',2);
      hold on;
      plot(pas,imag(Pola_norm(:,1,jpl,10)),'--','Color',pltc(jpl,:),'linewidth',2);
      xlim([min(pa),max(pa)]);
      xlabel(strpa);ylabel('J_y/(q_in_iV_A)');
      
      subplot(4,5,20); hold on; box on;
      plot(pas,real(Pola_norm(:,1,jpl,11)),'-','Color',pltc(jpl,:),'linewidth',2);
      hold on;
      plot(pas,imag(Pola_norm(:,1,jpl,11)),'--','Color',pltc(jpl,:),'linewidth',2);
      xlim([min(pa),max(pa)]);
      xlabel(strpa);ylabel('J_z/(q_in_iV_A)');
      
      
      
      %plot polarisation
      print(h1,'-dpng',[savepath,'fig_pdrk_',figstr,'_pola.png']);
      savefig([savepath,'fig_pdrk_',figstr,'_pola.fig']);
      
       %Edited by Duan 19/12/01 for magnetic helicity and compressibility
     
      tempBx = Pola(:,1,jpl,4);
      tempBy = Pola(:,1,jpl,5);
      tempBz = Pola(:,1,jpl,6);
      
      if ipa+ipb-2 == 0 %(scan k, fixed theta)
        kdirection = [sin(theta/180*pi) 0 cos(theta/180*pi)];
        Bpara = tempBx * kdirection(1) + tempBy * kdirection(2) + tempBz * kdirection(3);
        Bperp1 = tempBy;
        kperp2 = cross(kdirection,[0 1 0]);
        kperp2 = kperp2/norm(kperp2);
        Bperp2 = tempBx * kperp2(1) + tempBy * kperp2(2) + tempBz * kperp2(3);
      end
      
      Mcompress = abs(tempBz).^2 ./ (abs(Bpara).^2 + abs(Bperp1).^2 + abs(Bperp2).^2);
     % Mhelicity = 2*imag(Bperp1.*conj(Bperp2))./(abs(Bpara).^2 + abs(Bperp1).^2 + abs(Bperp2).^2);
     % Mhelicity = 2*cos(theta*pi/180)*1i*tempBx./tempBy./(cos(theta*pi/180)^2+abs(tempBx).^2./abs(tempBy).^2);
      Mhelicity = 2/cos(theta*pi/180)*real(1i*tempBx.*conj(tempBy))./(abs(tempBx).^2 + abs(tempBy).^2 + abs(tempBz).^2);
       h2 = figure;
      subplot(121)
       hold on; box on;
      plot(pas,Mcompress,'-',...
          'Color',pltc(jpl,:),'linewidth',2);
      
      xlim([min(pa),max(pa)]);
      xlabel(strpa);ylabel('dB_{//}^2/|dB|^2');
      
       subplot(122)
       hold on; box on;
      plot(pas,Mhelicity,'-',...
          'Color',pltc(jpl,:),'linewidth',2);
      
      xlim([min(pa),max(pa)]);
      xlabel(strpa);ylabel('\sigma_m');
      
      print(h2,'-dpng',[savepath,'fig_pdrk_',figstr,'_helicity.png']);
      savefig([savepath,'fig_pdrk_',figstr,'_helicity.fig']);
       %Edited by Duan 19/12/01 for magnetic helicity and compressibility
      
       %Added by Xingyu Zhu 2020-05-14
      run ../modules/pkues_plot_comp_velocity
      run ../modules/pkues_plot_growing_rate
      run ../modules/pkues_add_polarization_1
      run ../modules/pkues_add_polarization_2
      
    else % plot 2D polarizations, to update
      
      wwjp=squeeze(wws2(:,:,jpl));
      Polajp=squeeze(Pola(:,:,jpl,:));
      
      subplot(221);
      surf(ppa,ppb,real(wwjp)); hold on; box on; set(gca,'BoxStyle','full');
      xlabel([strpa,',ilogx=',num2str(iloga)]);
      ylabel([strpb,',ilogy=',num2str(ilogb)]);axis tight;
      zlabel('\omega_r/\omega_{c1}');
      
      subplot(222);
      surf(ppa,ppb,imag(wwjp)); hold on; box on;set(gca,'BoxStyle','full');
      xlabel(strpa); ylabel(strpb); axis tight;
      zlabel('\omega_i/\omega_{c1}');
    
      EyoveriEx=Polajp(:,:,2)./(1i*Polajp(:,:,1));
      
      subplot(223);
      surf(ppa,ppb,real(EyoveriEx));
      hold on; box on;set(gca,'BoxStyle','full');
      xlabel(strpa); ylabel(strpb); axis tight;
      zlabel('Re[E_y/(iE_x)]');
      
      subplot(223);
      surf(ppa,ppb,imag(EyoveriEx));
      hold on; box on;set(gca,'BoxStyle','full');
      xlabel(strpa); ylabel(strpb); axis tight;
      zlabel('Im[E_y/(iE_x)]');
      
    end
end

end

filename=['out_pdrk_S=',num2str(S),'_J=',num2str(J),...
    '_N=',num2str(N),'_B0=',num2str(B0),'.mat'];

allvars=whos;
tosave=cellfun(@isempty, regexp({allvars.class}, ...
    '^matlab\.(ui|graphics)\.')); % save workspace exclude figure(s)
save([savepath,filename], allvars(tosave).name); 
