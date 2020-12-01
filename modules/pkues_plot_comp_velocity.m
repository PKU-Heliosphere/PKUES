% 20-05-14 added by Xingyu Zhu and Jiansen He
% plot velocity Polarizations Comp1: Vy/iVx; Comp1/Comp2: Vx1/Vx2;
if (S==3) %for two ion components (core+beam)
h3 = figure('unit','normalized','Position',[0.01 0.1 0.9 0.8],...
  'DefaultAxesFontSize',10);
subplot(361); hold on; box on;
plot(pas,real(dVnorm(:,2,1,jpl)./(1i*dVnorm(:,1,1,jpl))),'-','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('Viy_{core}/(iVix_{core})');

subplot(362); hold on; box on;
plot(pas,real(dVnorm(:,3,1,jpl)./(1i*dVnorm(:,1,1,jpl))),'-','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('Viz_{core}/(iVix_{core})');

subplot(363); hold on; box on;
plot(pas,abs(dVnorm(:,1,1,jpl))./abs(dVnorm(:,1,2,jpl)),'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,zeros(npa)+1,'--','Color','b');
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('|Vix_{core}|/|Vix_{beam}|');

subplot(364); hold on; box on;
plot(pas,real(dVnorm(:,1,1,jpl)),'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,imag(dVnorm(:,1,1,jpl)),'--','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]);
xlabel(strpa);ylabel('Vx_{core}');
legend('Re','Im'); legend('boxoff');

subplot(365); hold on; box on;
plot(pas,real(dVnorm(:,2,1,jpl)),'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,imag(dVnorm(:,2,1,jpl)),'--','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]);
xlabel(strpa);ylabel('Vy_{core}');
legend('Re','Im'); legend('boxoff');

subplot(366); hold on; box on;
plot(pas,real(dVnorm(:,3,1,jpl)),'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,imag(dVnorm(:,3,1,jpl)),'--','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]);
xlabel(strpa);ylabel('Vz_{core}');
legend('Re','Im'); legend('boxoff');

      
subplot(367); hold on; box on;
plot(pas,real(dVnorm(:,2,2,jpl)./(1i*dVnorm(:,1,2,jpl))),'-','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('Viy_{beam}/(iVix_{beam})');

subplot(368); hold on; box on;
plot(pas,real(dVnorm(:,3,2,jpl)./(1i*dVnorm(:,1,2,jpl))),'-','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('Viz_{beam}/(iVix_{beam})');

subplot(369); hold on; box on;
plot(pas,abs(dVnorm(:,1,2,jpl))./abs(dVnorm(:,1,3,jpl)),'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,zeros(npa)+1,'--','Color','b');
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('|Vix_{beam}|/|Vex|');

subplot(3,6,10); hold on; box on;
plot(pas,real(dVnorm(:,1,2,jpl)),'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,imag(dVnorm(:,1,2,jpl)),'--','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]);
xlabel(strpa);ylabel('Vx_{beam}');
legend('Re','Im'); legend('boxoff');

subplot(3,6,11); hold on; box on;
plot(pas,real(dVnorm(:,2,2,jpl)),'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,imag(dVnorm(:,2,2,jpl)),'--','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]);
xlabel(strpa);ylabel('Vy_{beam}');
legend('Re','Im'); legend('boxoff');
subplot(3,6,11); hold on; box on;
plot(pas,real(dVnorm(:,2,2,jpl)),'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,imag(dVnorm(:,2,2,jpl)),'--','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]);
xlabel(strpa);ylabel('Vy_{beam}');
legend('Re','Im'); legend('boxoff');

subplot(3,6,12); hold on; box on;
plot(pas,real(dVnorm(:,3,2,jpl)),'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,imag(dVnorm(:,3,2,jpl)),'--','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]);
xlabel(strpa);ylabel('Vz_{beam}');
legend('Re','Im'); legend('boxoff');

subplot(3,6,13); hold on; box on;
plot(pas,real(dVnorm(:,2,3,jpl)./(1i*dVnorm(:,1,3,jpl))),'-','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('Vey/(iVex)');

subplot(3,6,14); hold on; box on;
plot(pas,real(dVnorm(:,3,3,jpl)./(1i*dVnorm(:,1,3,jpl))),'-','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('Vez/(iVex)');

subplot(3,6,15); hold on; box on;
plot(pas,abs(dVnorm(:,1,3,jpl))./abs(dVnorm(:,1,1,jpl)),'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,zeros(npa)+1,'--','Color','b');
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('|Vex|/|Vix_{core}|');

subplot(3,6,16); hold on; box on;
plot(pas,real(dVnorm(:,1,3,jpl)),'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,imag(dVnorm(:,1,3,jpl)),'--','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]);
xlabel(strpa);ylabel('Vx_e');
legend('Re','Im'); legend('boxoff');

subplot(3,6,17); hold on; box on;
plot(pas,real(dVnorm(:,2,3,jpl)),'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,imag(dVnorm(:,2,3,jpl)),'--','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]);
xlabel(strpa);ylabel('Vy_e');
legend('Re','Im'); legend('boxoff');

subplot(3,6,18); hold on; box on;
plot(pas,real(dVnorm(:,3,3,jpl)),'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,imag(dVnorm(:,3,3,jpl)),'--','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]);
xlabel(strpa);ylabel('Vz_e');
legend('Re','Im'); legend('boxoff');

end

if(S==2) %for one ion component 
h3 = figure('unit','normalized','Position',[0.01 0.1 0.9 0.6],...
  'DefaultAxesFontSize',15);
subplot(261); hold on; box on;
plot(pas,real(dVnorm(:,2,1,jpl)./(1i*dVnorm(:,1,1,jpl))),'-','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('V_{i,y}/(iV_{i,x})');

subplot(262); hold on; box on;
plot(pas,real(dVnorm(:,3,1,jpl)./(1i*dVnorm(:,1,1,jpl))),'-','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('V_{i,z}/(iV_{i,x})');
subplot(263); hold on; box on;
plot(pas,abs(dVnorm(:,1,1,jpl))./abs(dVnorm(:,1,2,jpl)),'-','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('|V_{i,x}|/|V_{e,x}|');

subplot(264); hold on; box on;
plot(pas,real(dVnorm(:,1,1,jpl)),'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,imag(dVnorm(:,1,1,jpl)),'--','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]);
xlabel(strpa);ylabel('V_{i,x}');
legend('Re','Im'); legend('boxoff');

subplot(265); hold on; box on;
plot(pas,real(dVnorm(:,2,1,jpl)),'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,imag(dVnorm(:,2,1,jpl)),'--','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]);
xlabel(strpa);ylabel('V_{i,y}');
legend('Re','Im'); legend('boxoff');

subplot(266); hold on; box on;
plot(pas,real(dVnorm(:,3,1,jpl)),'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,imag(dVnorm(:,3,1,jpl)),'--','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]);
xlabel(strpa);ylabel('V_{i,z}');
legend('Re','Im'); legend('boxoff');

      
subplot(267); hold on; box on;
plot(pas,real(dVnorm(:,2,2,jpl)./(1i*dVnorm(:,1,2,jpl))),'-','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('V_{e,y}/(iV_{e,x})');

subplot(268); hold on; box on;
plot(pas,real(dVnorm(:,3,2,jpl)./(1i*dVnorm(:,1,2,jpl))),'-','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('V_{e,z}/(iV_{e,x})');

subplot(269); hold on; box on;
plot(pas,abs(dVnorm(:,1,1,jpl))./abs(dVnorm(:,1,2,jpl)),'-','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('|V_{i,y}|/|V_{e,y}|');

subplot(2,6,10); hold on; box on;
plot(pas,real(dVnorm(:,1,2,jpl)),'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,imag(dVnorm(:,1,2,jpl)),'--','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]);
xlabel(strpa);ylabel('V_{e,x}');
legend('Re','Im'); legend('boxoff');

subplot(2,6,11); hold on; box on;
plot(pas,real(dVnorm(:,2,2,jpl)),'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,imag(dVnorm(:,2,2,jpl)),'--','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]);
xlabel(strpa);ylabel('V_{e,y}');
legend('Re','Im'); legend('boxoff');

subplot(2,6,12); hold on; box on;
plot(pas,real(dVnorm(:,3,2,jpl)),'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,imag(dVnorm(:,3,2,jpl)),'--','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]);
xlabel(strpa);ylabel('V_{e,z}');
legend('Re','Im'); legend('boxoff');

end

print(h3,'-dpng',[savepath,'fig_pdrk_',figstr,'_velocity.png']);
savefig([savepath,'fig_pdrk_',figstr,'_velocity.fig']);