% 20-05-20 added by Xingyu Zhu
% plot growing rate according to Eq.(3) of He et al. (2019)
if (S==3) %for two ion components (core+beam)
h3 = figure('unit','normalized','Position',[0.01 0.1 0.7 0.8],...
  'DefaultAxesFontSize',15);
subplot(341); hold on; box on;
plot(pas,JE(:,1,1,jpl),'-','Color',pltc(jpl,:),'linewidth',2);
symlog('y');
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('\gamma_{x,core} [s^{-1}]');
title('-J\cdotE/(dB^2+dE^2)/2')

subplot(342); hold on; box on;
plot(pas,JE(:,2,1,jpl),'-','Color',pltc(jpl,:),'linewidth',2);
symlog('y');
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('\gamma_{y,core}');

subplot(343); hold on; box on;
plot(pas,JE(:,3,1,jpl),'-','Color',pltc(jpl,:),'linewidth',2);
symlog('y');
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('\gamma_{z,core}');

subplot(344); hold on; box on;
plot(pas,(JE(:,1,1,jpl)+JE(:,2,1,jpl)+JE(:,3,1,jpl))./(JE(:,1,2,jpl)+JE(:,2,2,jpl)+JE(:,3,2,jpl)),...
    '-','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('\gamma_{trace,core}/\gamma_{trace,beam}');
   
subplot(345); hold on; box on;
plot(pas,JE(:,1,2,jpl),'-','Color',pltc(jpl,:),'linewidth',2);
symlog('y');
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('\gamma_{x,beam}');

subplot(346); hold on; box on;
plot(pas,JE(:,2,2,jpl),'-','Color',pltc(jpl,:),'linewidth',2);
symlog('y');
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('\gamma_{y,beam}');

subplot(347); hold on; box on;
plot(pas,JE(:,3,2,jpl),'-','Color',pltc(jpl,:),'linewidth',2);
symlog('y');
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('\gamma_{z,beam}');

subplot(348); hold on; box on;
plot(pas,(JE(:,1,2,jpl)+JE(:,2,2,jpl)+JE(:,3,2,jpl))./(JE(:,1,3,jpl)+JE(:,2,3,jpl)+JE(:,3,3,jpl)),...
    '-','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('\gamma_{trace,beam}/\gamma_{trace,e}');

subplot(349); hold on; box on;
plot(pas,JE(:,1,3,jpl),'-','Color',pltc(jpl,:),'linewidth',2);
symlog('y');
xlim([min(pa),max(pa)]);
xlabel(strpa);ylabel('\gamma_{x,electron}');

subplot(3,4,10); hold on; box on;
plot(pas,JE(:,2,3,jpl),'-','Color',pltc(jpl,:),'linewidth',2);
symlog('y');
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('\gamma_{y,electron}');

subplot(3,4,11); hold on; box on;
plot(pas,JE(:,3,3,jpl),'-','Color',pltc(jpl,:),'linewidth',2);
symlog('y');
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('\gamma_{z,electron}');

subplot(3,4,12); hold on; box on;
plot(pas,(JE(:,1,1,jpl)+JE(:,2,1,jpl)+JE(:,3,1,jpl))./(JE(:,1,3,jpl)+JE(:,2,3,jpl)+JE(:,3,3,jpl)),...
    '-','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('\gamma_{trace,core}/\gamma_{trace,e}');
end

if(S==2) %for one ion component
h3 = figure('unit','normalized','Position',[0.01 0.1 0.7 0.6],...
  'DefaultAxesFontSize',15);
subplot(241); hold on; box on;
plot(pas,JE(:,1,1,jpl),'-','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('\gamma_{x,proton} [s^{-1}]');
title('-J\cdotE/(dB^2+dE^2)/2')

subplot(242); hold on; box on;
plot(pas,JE(:,2,1,jpl),'-','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('\gamma_{y,proton}');

subplot(243); hold on; box on;
plot(pas,JE(:,3,1,jpl),'-','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('\gamma_{z,proton}');

subplot(244); hold on; box on;
plot(pas,(JE(:,1,1,jpl)+JE(:,2,1,jpl)+JE(:,3,1,jpl))./(JE(:,1,2,jpl)+JE(:,2,2,jpl)+JE(:,3,2,jpl)),...
    '-','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('\gamma_{trace,i}/\gamma_{trace,e}');

subplot(245); hold on; box on;
plot(pas,JE(:,1,2,jpl),'-','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('\gamma_{x,electron}');

subplot(246); hold on; box on;
plot(pas,JE(:,2,2,jpl),'-','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('\gamma_{y,electron}');

subplot(247); hold on; box on;
plot(pas,JE(:,3,2,jpl),'-','Color',pltc(jpl,:),'linewidth',2);
xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
xlabel(strpa);ylabel('\gamma_{z,electron}');

end

print(h3,'-dpng',[savepath,'fig_pdrk_',figstr,'_growingrate.png']);
savefig([savepath,'fig_pdrk_',figstr,'_growingrate.fig']);