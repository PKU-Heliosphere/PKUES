% 18-10-05 08:00 Hua-sheng XIE, huashengxie@gmail.com, FRI-ENN, China
% Ackn.: Richard Denton (Dartmouth), Xin Tao (USTC), Jin-song Zhao (PMO),
% etc ...
% This file plot all the solution of pdrk-em3d results, which can be used
% to select the initial data for separating different dispersion surface
%=================================================================================
% edited by Duan 2019/11/28
% draw the 1-D dispersion relation in different color
%
close all;

h=figure('unit','normalized','Position',[0.01 0.45 0.5 0.3],...
  'DefaultAxesFontSize',15);

% normalized omega and k
wwn=ww/wcs1;
kkn=kk*cwp;
kxxn=kxx*cwp;
kzzn=kzz*cwp;

f1 = figure(1);

if(ipa==ipb) % 1D plot

  %kkk=reshape(repmat(kkn,1,nw0),1,[]);
  papp=reshape(repmat(ppa,1,nw0),1,[]);
  www=reshape(wwn,1,[]);
  npp=length(papp);

  sp1 = subplot(121);
  % Edited by Duan 19/11/28
  color_arr = hsv(nw0);
  wiplots = cell(nw0,1);
  for iplot = 1:nw0
  if(iloga==0)
    %plot(papp,real(www),'g+','LineWidth',2); hold on;
    plot(ppa,real(wwn(:,1,iplot)),'Color',color_arr(iplot,:),'Marker','+','LineStyle','None','LineWidth',2); hold on;
  else
    %semilogx(10.^papp,real(www),'g+','LineWidth',2); hold on;
    wiplots{iplot} = loglog(10.^ppa,real(wwn(:,1,iplot)),'Color',color_arr(iplot,:),'Marker','+','LineStyle','None','LineWidth',2); hold on;
  end  
  end
  ylim([-1 10])
  % Edited by Duan 19/28/11
  xlabel([strpa,', runtime=',num2str(runtime),'s']); 
  ylabel(['\omega_r/\omega_{c1}, npa=',num2str(npa),',npb=',num2str(npb)]);
  title(['(a) \beta_{||}=',num2str(betasz,3),...
      ', \beta_\perp=',num2str(betasp,3)]);
  %xlim([min(10.^pa) max(10.^pa)]); 
  box on; %ylim([-2.5,2.5]);
  
  sp2 = subplot(122);
   % Edited by Duan 19/11/28
  for iplot = 1:nw0
  if(iloga==0)
    %plot(papp,imag(www),'g+','LineWidth',2); hold on;
    plot(ppa,imag(wwn(:,1,iplot)),'Color',color_arr(iplot,:),'Marker','+','LineStyle','None','LineWidth',2); hold on;
  else
    %semilogx(10.^papp,imag(www),'g+','LineWidth',2); hold on;
    loglog(10.^ppa,imag(wwn(:,1,iplot)),'Color',color_arr(iplot,:),'Marker','+','LineStyle','None','LineWidth',2); hold on;
  end
  end
   % Edited by Duan 19/11/28
  ylim([-2 0.5])
  xlabel([strpa,', (S=',num2str(S),',N=',num2str(N),',J=',num2str(J),')']); 
  ylabel('\omega_i/\omega_{c1}');
  title(['(b) v_A/c=',num2str(vA/sqrt(c2),2),', ',strpb,'=',...
    num2str(par(ipbtmp))]);
  %xlim([min(10.^pa),max(10.^pa)]); 
  box on; %ylim([-1.0,0.1]);
  
else % 2D plot
  for jp=1:25 % change here to plot more surfaces
    subplot(121);
    wwjp=squeeze(wwn(:,:,jp));
    surf(ppa,ppb,real(wwjp)); hold on; box on;
    xlabel([strpa,',ilogx=',num2str(iloga)]);
    ylabel([strpb,',ilogy=',num2str(ilogb)]);
    zlabel(['\omega_r/\omega_{c1},npa=',num2str(npa),',npb=',num2str(npb)]);
    title(['(a) \beta_{||}=',num2str(betasz,3),...
      ', \beta_\perp=',num2str(betasp,3)]);
    subplot(122);
    surf(ppa,ppb,imag(wwjp)); hold on; box on;
    xlabel(strpa); ylabel(strpb);
    zlabel(['\omega_i/\omega_{c1},N=',num2str(N),',J=',num2str(J)]);
    title(['(b) runtime=',num2str(runtime),'s']);
    %%
    % zoom in the figure to find dispersion surface data for plot_select.m
    %zlim([-0.5e0,0.1]);
    %%
  end
  
end

figstr=['S=',num2str(S),'_J=',num2str(J),'_N=',num2str(N),...
    '_npa=',num2str(npa),'_npb=',num2str(npb)];
print(gcf,'-dpng',[savepath,'fig_pdrk_',figstr,'_all.png']);
savefig([savepath,'fig_pdrk_',figstr,'_all.fig']);

%Added by Die Duan 2019/12/13
prompt = 'Please select the wave, for wr input 1, for wi input 2. Input others will end the selection.';
computeflag = input(prompt);
while computeflag == 1 || computeflag == 2
    dcm_obj = datacursormode(f1);
    datatipinfo = getCursorInfo(dcm_obj);
    wvalue = datatipinfo.Position(2);
    kindex = datatipinfo.DataIndex;        
    if computeflag == 1
        [~,wwwindex] = min(abs(real(www)-wvalue));
        plotindex = floor(wwwindex/length(ppa))+1;
        subplot(122);
        tempk = papp(wwwindex);
        tempwi = imag(www(wwwindex));
        disp(['kdi=',num2str(tempk),'wi=',num2str(tempwi)]);
        if iloga == 0
            plot(tempk,tempwi,'kh','MarkerSize',3);
            plot(tempk,tempwi,'kh','MarkerSize',10);
            plot(tempk,tempwi,'kh','MarkerSize',20);
            xlim([tempk-0.5 tempk+0.5])
            wpdat=[tempk,0,wvalue;];
        else
        loglog(10.^tempk,tempwi,'kh','MarkerSize',3);
        loglog(10.^tempk,tempwi,'kh','MarkerSize',10);
        loglog(10.^tempk,tempwi,'kh','MarkerSize',20);
        xlim([10^(tempk-2),10^(tempk+2)])
        if tempwi > 0
            ylim([tempwi/100,tempwi*100])
        else
            ylim([tempwi*100,tempwi/100])
        end
        wpdat=[10.^tempk,0,wvalue;];
        end
        
        rex=1;
        rey=1;
        rez=1;
        
        
    else
        if computeflag == 2
        [~,wwwindex] = min(abs(imag(www)-wvalue));
        plotindex = floor(wwwindex/length(ppa))+1;
        subplot(121);
        tempk = papp(wwwindex);
        tempwi = real(www(wwwindex));
        disp(['kdi=',num2str(tempk),'wr=',num2str(tempwi)]);
        if iloga == 0
            plot(tempk,tempwi,'kh','MarkerSize',3);
            plot(tempk,tempwi,'kh','MarkerSize',10);
            plot(tempk,tempwi,'kh','MarkerSize',20);
            xlim([tempk-0.5 tempk+0.5])
            wpdat=[tempk,0,wvalue*1i;];
        else
            loglog(10.^tempk,tempwi,'kh','MarkerSize',3);
            loglog(10.^tempk,tempwi,'kh','MarkerSize',10);
            loglog(10.^tempk,tempwi,'kh','MarkerSize',20);
            xlim([10^(tempk-2),10^(tempk+2)])
            if tempwi > 0
                ylim([tempwi/100,tempwi*100])
            else
                ylim([tempwi*100,tempwi/100])
            end
        wpdat=[10.^tempk,0,wvalue*1i;];
        end
        
        rex=1;
        rey=1;
        rez=1;
        end
    end
    run ../modules/pdrk_plot_select;
    run ../modules/pdrk_output;
    prompt = 'Please select the next wave, for wr input 1, for wi input 2. Input others will end the selection.';
    computeflag = input(prompt);
end