% 20-05-20 Coded by Xingyu Zhu and Jiansen He
% plot pseudo growing rate according to Eq.(3) of He et al. (2019)
if (S==3) 
%%%%%%%
phi_Bx_Vex = mod((angle(Pola_norm(:,1,1,4))-angle(dV(:,1,3,1)))*180/pi+360,360);
phi_By_Vey = mod((angle(Pola_norm(:,1,1,5))-angle(dV(:,2,3,1)))*180/pi+360,360);
phi_Bz_Vez = mod((angle(Pola_norm(:,1,1,6))-angle(dV(:,3,3,1)))*180/pi+360,360);
phi_Ex_Vex = mod((angle(Pola_norm(:,1,1,1))-angle(dV(:,1,3,1)))*180/pi+360,360);
phi_Ey_Vey = mod((angle(Pola_norm(:,1,1,2))-angle(dV(:,2,3,1)))*180/pi+360,360);
phi_Ez_Vez = mod((angle(Pola_norm(:,1,1,3))-angle(dV(:,3,3,1)))*180/pi+360,360);
phi_Bx_Ex = mod((angle(Pola_norm(:,1,1,4))-angle(Pola_norm(:,1,1,1)))*180/pi+360,360);
phi_By_Ey = mod((angle(Pola_norm(:,1,1,5))-angle(Pola_norm(:,1,1,2)))*180/pi+360,360);
phi_Bz_Ez = mod((angle(Pola_norm(:,1,1,6))-angle(Pola_norm(:,1,1,3)))*180/pi+360,360);
h3 = figure('unit','normalized','Position',[0.01 0.1 0.7 0.8],...
  'DefaultAxesFontSize',15);
subplot(331); hold on; box on;
plot(pas,phi_Bx_Vex,'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,zeros(npa),'--','Color','b');hold on;
plot(pas,zeros(npa)+90,'--','Color','b');hold on;
plot(pas,zeros(npa)+180,'--','Color','b');hold on;
plot(pas,zeros(npa)+270,'--','Color','b');hold on;
plot(pas,zeros(npa)+360,'--','Color','b');hold on;
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('$\phi(Bx-Vex) (\circ)$','Interpreter','Latex');

subplot(332); hold on; box on;
plot(pas,phi_By_Vey,'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,zeros(npa),'--','Color','b');hold on;
plot(pas,zeros(npa)+90,'--','Color','b');hold on;
plot(pas,zeros(npa)+180,'--','Color','b');hold on;
plot(pas,zeros(npa)+270,'--','Color','b');hold on;
plot(pas,zeros(npa)+360,'--','Color','b');hold on;
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('$\phi(By-Vey) (\circ)$','Interpreter','Latex');

subplot(333); hold on; box on;
plot(pas,phi_Bz_Vez,'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,zeros(npa),'--','Color','b');hold on;
plot(pas,zeros(npa)+90,'--','Color','b');hold on;
plot(pas,zeros(npa)+180,'--','Color','b');hold on;
plot(pas,zeros(npa)+270,'--','Color','b');hold on;
plot(pas,zeros(npa)+360,'--','Color','b');hold on;
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('$\phi(Bz-Vez) (\circ)$','Interpreter','Latex');

subplot(334); hold on; box on;
plot(pas,phi_Ex_Vex,'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,zeros(npa),'--','Color','b');hold on;
plot(pas,zeros(npa)+90,'--','Color','b');hold on;
plot(pas,zeros(npa)+180,'--','Color','b');hold on;
plot(pas,zeros(npa)+270,'--','Color','b');hold on;
plot(pas,zeros(npa)+360,'--','Color','b');hold on;
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('$\phi(Ex-Vex) (\circ)$','Interpreter','Latex');
   
subplot(335); hold on; box on;
plot(pas,phi_Ey_Vey,'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,zeros(npa),'--','Color','b');hold on;
plot(pas,zeros(npa)+90,'--','Color','b');hold on;
plot(pas,zeros(npa)+180,'--','Color','b');hold on;
plot(pas,zeros(npa)+270,'--','Color','b');hold on;
plot(pas,zeros(npa)+360,'--','Color','b');hold on;
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('$\phi(Ey-Vey) (\circ)$','Interpreter','Latex');

subplot(336); hold on; box on;
plot(pas,phi_Ez_Vez,'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,zeros(npa),'--','Color','b');hold on;
plot(pas,zeros(npa)+90,'--','Color','b');hold on;
plot(pas,zeros(npa)+180,'--','Color','b');hold on;
plot(pas,zeros(npa)+270,'--','Color','b');hold on;
plot(pas,zeros(npa)+360,'--','Color','b');hold on;
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('$\phi(Ez-Vez) (\circ)$','Interpreter','Latex');

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
phi_Bx_Vex = mod((angle(Pola_norm(:,1,1,4))-angle(dV(:,1,2,1)))*180/pi+360,360);
phi_By_Vey = mod((angle(Pola_norm(:,1,1,5))-angle(dV(:,2,2,1)))*180/pi+360,360);
phi_Bz_Vez = mod((angle(Pola_norm(:,1,1,6))-angle(dV(:,3,2,1)))*180/pi+360,360);
phi_Ex_Vex = mod((angle(Pola_norm(:,1,1,1))-angle(dV(:,1,2,1)))*180/pi+360,360);
phi_Ey_Vey = mod((angle(Pola_norm(:,1,1,2))-angle(dV(:,2,2,1)))*180/pi+360,360);
phi_Ez_Vez = mod((angle(Pola_norm(:,1,1,3))-angle(dV(:,3,2,1)))*180/pi+360,360);
phi_Bx_Ex = mod((angle(Pola_norm(:,1,1,4))-angle(Pola_norm(:,1,1,1)))*180/pi+360,360);
phi_By_Ey = mod((angle(Pola_norm(:,1,1,5))-angle(Pola_norm(:,1,1,2)))*180/pi+360,360);
phi_Bz_Ez = mod((angle(Pola_norm(:,1,1,6))-angle(Pola_norm(:,1,1,3)))*180/pi+360,360);        
h3 = figure('unit','normalized','Position',[0.01 0.1 0.7 0.8],...
  'DefaultAxesFontSize',15);
subplot(331); hold on; box on;
plot(pas,phi_Bx_Vex,'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,zeros(npa),'--','Color','b');hold on;
plot(pas,zeros(npa)+90,'--','Color','b');hold on;
plot(pas,zeros(npa)+180,'--','Color','b');hold on;
plot(pas,zeros(npa)+270,'--','Color','b');hold on;
plot(pas,zeros(npa)+360,'--','Color','b');hold on;
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('$\phi(Bx-Vex) (\circ)$','Interpreter','Latex');

subplot(332); hold on; box on;
plot(pas,phi_By_Vey,'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,zeros(npa),'--','Color','b');hold on;
plot(pas,zeros(npa)+90,'--','Color','b');hold on;
plot(pas,zeros(npa)+180,'--','Color','b');hold on;
plot(pas,zeros(npa)+270,'--','Color','b');hold on;
plot(pas,zeros(npa)+360,'--','Color','b');hold on;
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('$\phi(By-Vey) (\circ)$','Interpreter','Latex');

subplot(333); hold on; box on;
plot(pas,phi_Bz_Vez,'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,zeros(npa),'--','Color','b');hold on;
plot(pas,zeros(npa)+90,'--','Color','b');hold on;
plot(pas,zeros(npa)+180,'--','Color','b');hold on;
plot(pas,zeros(npa)+270,'--','Color','b');hold on;
plot(pas,zeros(npa)+360,'--','Color','b');hold on;
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('$\phi(Bz-Vez) (\circ)$','Interpreter','Latex');

subplot(334); hold on; box on;
plot(pas,phi_Ex_Vex,'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,zeros(npa),'--','Color','b');hold on;
plot(pas,zeros(npa)+90,'--','Color','b');hold on;
plot(pas,zeros(npa)+180,'--','Color','b');hold on;
plot(pas,zeros(npa)+270,'--','Color','b');hold on;
plot(pas,zeros(npa)+360,'--','Color','b');hold on;
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('$\phi(Ex-Vex) (\circ)$','Interpreter','Latex');

subplot(335); hold on; box on;
plot(pas,phi_Ey_Vey,'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,zeros(npa),'--','Color','b');hold on;
plot(pas,zeros(npa)+90,'--','Color','b');hold on;
plot(pas,zeros(npa)+180,'--','Color','b');hold on;
plot(pas,zeros(npa)+270,'--','Color','b');hold on;
plot(pas,zeros(npa)+360,'--','Color','b');hold on;
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('$\phi(Ey-Vey) (\circ)$','Interpreter','Latex');

subplot(336); hold on; box on;
plot(pas,phi_Ez_Vez,'-','Color',pltc(jpl,:),'linewidth',2);
hold on;
plot(pas,zeros(npa),'--','Color','b');hold on;
plot(pas,zeros(npa)+90,'--','Color','b');hold on;
plot(pas,zeros(npa)+180,'--','Color','b');hold on;
plot(pas,zeros(npa)+270,'--','Color','b');hold on;
plot(pas,zeros(npa)+360,'--','Color','b');hold on;
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('$\phi(Ez-Vez) (\circ)$','Interpreter','Latex');

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

print(h3,'-dpng',[savepath,'fig_pdrk_',figstr,'_BEVepola.png']);
savefig([savepath,'fig_pdrk_',figstr,'_BEVepola.fig']);