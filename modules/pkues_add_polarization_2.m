% 20-05-20 added by Xingyu Zhu and Jiansen He
% plot growing rate according to Eq.(3) of He et al. (2019)
if(S==3) %for core ions
%%%%%%%
phi_Bx_Vix = mod((angle(Pola_norm(:,1,1,4))-angle(dV(:,1,1,1)))*180/pi+360,360);
% phi_Bx_Vix = (angle(Pola_norm(:,1,1,4))-angle(dV(:,1,1,1)))*180/pi;
phi_By_Viy = mod((angle(Pola_norm(:,1,1,5))-angle(dV(:,2,1,1)))*180/pi+360,360);
phi_Bz_Viz = mod((angle(Pola_norm(:,1,1,6))-angle(dV(:,3,1,1)))*180/pi+360,360);
phi_Ex_Vix = mod((angle(Pola_norm(:,1,1,1))-angle(dV(:,1,1,1)))*180/pi+360,360);
phi_Ey_Viy = mod((angle(Pola_norm(:,1,1,2))-angle(dV(:,2,1,1)))*180/pi+360,360);
phi_Ez_Viz = mod((angle(Pola_norm(:,1,1,3))-angle(dV(:,3,1,1)))*180/pi+360,360);
phi_Bx_Ex = mod((angle(Pola_norm(:,1,1,4))-angle(Pola_norm(:,1,1,1)))*180/pi+360,360);
phi_By_Ey = mod((angle(Pola_norm(:,1,1,5))-angle(Pola_norm(:,1,1,2)))*180/pi+360,360);
phi_Bz_Ez = mod((angle(Pola_norm(:,1,1,6))-angle(Pola_norm(:,1,1,3)))*180/pi+360,360);
h3 = figure('unit','normalized','Position',[0.01 0.1 0.7 0.8],...
  'DefaultAxesFontSize',15);
subplot(331); hold on; box on;
plot(pas,phi_Bx_Vix,'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,zeros(npa),'--','Color','b');hold on;
plot(pas,zeros(npa)+90,'--','Color','b');hold on;
plot(pas,zeros(npa)+180,'--','Color','b');hold on;
plot(pas,zeros(npa)+270,'--','Color','b');hold on;
plot(pas,zeros(npa)+360,'--','Color','b');hold on;
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('$\phi(Bx-Vix) (\circ)$','Interpreter','Latex');

subplot(332); hold on; box on;
plot(pas,phi_By_Viy,'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,zeros(npa),'--','Color','b');hold on;
plot(pas,zeros(npa)+90,'--','Color','b');hold on;
plot(pas,zeros(npa)+180,'--','Color','b');hold on;
plot(pas,zeros(npa)+270,'--','Color','b');hold on;
plot(pas,zeros(npa)+360,'--','Color','b');hold on;
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('$\phi(By-Viy) (\circ)$','Interpreter','Latex');

subplot(333); hold on; box on;
plot(pas,phi_Bz_Viz,'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,zeros(npa),'--','Color','b');hold on;
plot(pas,zeros(npa)+90,'--','Color','b');hold on;
plot(pas,zeros(npa)+180,'--','Color','b');hold on;
plot(pas,zeros(npa)+270,'--','Color','b');hold on;
plot(pas,zeros(npa)+360,'--','Color','b');hold on;
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('$\phi(Bz-Viz) (\circ)$','Interpreter','Latex');

subplot(334); hold on; box on;
plot(pas,phi_Ex_Vix,'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,zeros(npa),'--','Color','b');hold on;
plot(pas,zeros(npa)+90,'--','Color','b');hold on;
plot(pas,zeros(npa)+180,'--','Color','b');hold on;
plot(pas,zeros(npa)+270,'--','Color','b');hold on;
plot(pas,zeros(npa)+360,'--','Color','b');hold on;
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('$\phi(Ex-Vix) (\circ)$','Interpreter','Latex');

subplot(335); hold on; box on;
plot(pas,phi_Ey_Viy,'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,zeros(npa),'--','Color','b');hold on;
plot(pas,zeros(npa)+90,'--','Color','b');hold on;
plot(pas,zeros(npa)+180,'--','Color','b');hold on;
plot(pas,zeros(npa)+270,'--','Color','b');hold on;
plot(pas,zeros(npa)+360,'--','Color','b');hold on;
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('$\phi(Ey-Viy) (\circ)$','Interpreter','Latex');

subplot(336); hold on; box on;
plot(pas,phi_Ez_Viz,'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,zeros(npa),'--','Color','b');hold on;
plot(pas,zeros(npa)+90,'--','Color','b');hold on;
plot(pas,zeros(npa)+180,'--','Color','b');hold on;
plot(pas,zeros(npa)+270,'--','Color','b');hold on;
plot(pas,zeros(npa)+360,'--','Color','b');hold on;
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('$\phi(Ez-Viz) (\circ)$','Interpreter','Latex');

subplot(337); hold on; box on;
plot(pas,phi_Bx_Ex,'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,zeros(npa),'--','Color','b');hold on;
plot(pas,zeros(npa)+90,'--','Color','b');hold on;
plot(pas,zeros(npa)+180,'--','Color','b');hold on;
plot(pas,zeros(npa)+270,'--','Color','b');hold on;
plot(pas,zeros(npa)+360,'--','Color','b');hold on;
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('$\phi(Bx-Ex) (\circ)$','Interpreter','Latex');

subplot(338); hold on; box on;
plot(pas,phi_By_Ey,'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,zeros(npa),'--','Color','b');hold on;
plot(pas,zeros(npa)+90,'--','Color','b');hold on;
plot(pas,zeros(npa)+180,'--','Color','b');hold on;
plot(pas,zeros(npa)+270,'--','Color','b');hold on;
plot(pas,zeros(npa)+360,'--','Color','b');hold on;
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('$\phi(By-Ey) (\circ)$','Interpreter','Latex');

subplot(339); hold on; box on;
plot(pas,phi_Bz_Ez,'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,zeros(npa),'--','Color','b');hold on;
plot(pas,zeros(npa)+90,'--','Color','b');hold on;
plot(pas,zeros(npa)+180,'--','Color','b');hold on;
plot(pas,zeros(npa)+270,'--','Color','b');hold on;
plot(pas,zeros(npa)+360,'--','Color','b');hold on;
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('$\phi(Bz-Ez) (\circ)$','Interpreter','Latex');

end
if(S==2)
phi_Bx_Vix = mod((angle(Pola_norm(:,1,1,4))-angle(dV(:,1,1,1)))*180/pi+360,360);
phi_By_Viy = mod((angle(Pola_norm(:,1,1,5))-angle(dV(:,2,1,1)))*180/pi+360,360);
phi_Bz_Viz = mod((angle(Pola_norm(:,1,1,6))-angle(dV(:,3,1,1)))*180/pi+360,360);
phi_Ex_Vix = mod((angle(Pola_norm(:,1,1,1))-angle(dV(:,1,1,1)))*180/pi+360,360);
phi_Ey_Viy = mod((angle(Pola_norm(:,1,1,2))-angle(dV(:,2,1,1)))*180/pi+360,360);
phi_Ez_Viz = mod((angle(Pola_norm(:,1,1,3))-angle(dV(:,3,1,1)))*180/pi+360,360);
phi_Bx_Ex = mod((angle(Pola_norm(:,1,1,4))-angle(Pola_norm(:,1,1,1)))*180/pi+360,360);
phi_By_Ey = mod((angle(Pola_norm(:,1,1,5))-angle(Pola_norm(:,1,1,2)))*180/pi+360,360);
phi_Bz_Ez = mod((angle(Pola_norm(:,1,1,6))-angle(Pola_norm(:,1,1,3)))*180/pi+360,360);    
h3 = figure('unit','normalized','Position',[0.01 0.1 0.7 0.8],...
  'DefaultAxesFontSize',15);
subplot(331); hold on; box on;
plot(pas,phi_Bx_Vix,'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,zeros(npa),'--','Color','b');hold on;
plot(pas,zeros(npa)+90,'--','Color','b');hold on;
plot(pas,zeros(npa)+180,'--','Color','b');hold on;
plot(pas,zeros(npa)+270,'--','Color','b');hold on;
plot(pas,zeros(npa)+360,'--','Color','b');hold on;
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('$\phi(Bx-Vix) (\circ)$','Interpreter','Latex');

subplot(332); hold on; box on;
plot(pas,phi_By_Viy,'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,zeros(npa),'--','Color','b');hold on;
plot(pas,zeros(npa)+90,'--','Color','b');hold on;
plot(pas,zeros(npa)+180,'--','Color','b');hold on;
plot(pas,zeros(npa)+270,'--','Color','b');hold on;
plot(pas,zeros(npa)+360,'--','Color','b');hold on;
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('$\phi(By-Viy) (\circ)$','Interpreter','Latex');

subplot(333); hold on; box on;
plot(pas,phi_Bz_Viz,'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,zeros(npa),'--','Color','b');hold on;
plot(pas,zeros(npa)+90,'--','Color','b');hold on;
plot(pas,zeros(npa)+180,'--','Color','b');hold on;
plot(pas,zeros(npa)+270,'--','Color','b');hold on;
plot(pas,zeros(npa)+360,'--','Color','b');hold on;
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('$\phi(Bz-Viz) (\circ)$','Interpreter','Latex');

subplot(334); hold on; box on;
plot(pas,phi_Ex_Vix,'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,zeros(npa),'--','Color','b');hold on;
plot(pas,zeros(npa)+90,'--','Color','b');hold on;
plot(pas,zeros(npa)+180,'--','Color','b');hold on;
plot(pas,zeros(npa)+270,'--','Color','b');hold on;
plot(pas,zeros(npa)+360,'--','Color','b');hold on;
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('$\phi(Ex-Vix) (\circ)$','Interpreter','Latex');

subplot(335); hold on; box on;
plot(pas,phi_Ey_Viy,'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,zeros(npa),'--','Color','b');hold on;
plot(pas,zeros(npa)+90,'--','Color','b');hold on;
plot(pas,zeros(npa)+180,'--','Color','b');hold on;
plot(pas,zeros(npa)+270,'--','Color','b');hold on;
plot(pas,zeros(npa)+360,'--','Color','b');hold on;
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('$\phi(Ey-Viy) (\circ)$','Interpreter','Latex');

subplot(336); hold on; box on;
plot(pas,phi_Ez_Viz,'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,zeros(npa),'--','Color','b');hold on;
plot(pas,zeros(npa)+90,'--','Color','b');hold on;
plot(pas,zeros(npa)+180,'--','Color','b');hold on;
plot(pas,zeros(npa)+270,'--','Color','b');hold on;
plot(pas,zeros(npa)+360,'--','Color','b');hold on;
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('$\phi(Ez-Viz) (\circ)$','Interpreter','Latex');

subplot(337); hold on; box on;
plot(pas,phi_Bx_Ex,'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,zeros(npa),'--','Color','b');hold on;
plot(pas,zeros(npa)+90,'--','Color','b');hold on;
plot(pas,zeros(npa)+180,'--','Color','b');hold on;
plot(pas,zeros(npa)+270,'--','Color','b');hold on;
plot(pas,zeros(npa)+360,'--','Color','b');hold on;
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('$\phi(Bx-Ex) (\circ)$','Interpreter','Latex');

subplot(338); hold on; box on;
plot(pas,phi_By_Ey,'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,zeros(npa),'--','Color','b');hold on;
plot(pas,zeros(npa)+90,'--','Color','b');hold on;
plot(pas,zeros(npa)+180,'--','Color','b');hold on;
plot(pas,zeros(npa)+270,'--','Color','b');hold on;
plot(pas,zeros(npa)+360,'--','Color','b');hold on;
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('$\phi(By-Ey) (\circ)$','Interpreter','Latex');

subplot(339); hold on; box on;
plot(pas,phi_Bz_Ez,'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,zeros(npa),'--','Color','b');hold on;
plot(pas,zeros(npa)+90,'--','Color','b');hold on;
plot(pas,zeros(npa)+180,'--','Color','b');hold on;
plot(pas,zeros(npa)+270,'--','Color','b');hold on;
plot(pas,zeros(npa)+360,'--','Color','b');hold on;
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('$\phi(Bz-Ez) (\circ)$','Interpreter','Latex');

end

print(h3,'-dpng',[savepath,'fig_pdrk_',figstr,'_BEVipola.png']);
savefig([savepath,'fig_pdrk_',figstr,'_BEVipola.fig']);