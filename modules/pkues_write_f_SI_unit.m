% written by Xingyu Zhu 2020-06-01
% last modified: 2020-06-19
%%%%%%%%%%%%%%%%%%%%%%%
VXarray = zeros(vxsteps+1,vysteps+1,vzsteps+1);
VYarray = zeros(vxsteps+1,vysteps+1,vzsteps+1);                                                                                                                                                                                                                     rray = zeros(vxsteps+1,vysteps+1,vzsteps+1);
VZarray = zeros(vxsteps+1,vysteps+1,vzsteps+1);
f0array = zeros(vxsteps+1,vysteps+1,vzsteps+1);
farray  = zeros(vxsteps+1,vysteps+1,vzsteps+1);
imfarray= zeros(vxsteps+1,vysteps+1,vzsteps+1);
dfarray = zeros(vxsteps+1,vysteps+1,vzsteps+1);
% dfqEarray = zeros(vxsteps+1,vysteps+1,vzsteps+1);
% dfqvBarray = zeros(vxsteps+1,vysteps+1,vzsteps+1);

%%%%%%%%%%%%%%%%%%%%%%%
x = wws(jpa_df,jpb_df,jpl_df)*wcs1;
dEx = Pola_SI(jpa_df,jpb_df,jpl_df,1)/scaling_factor;
dEy = Pola_SI(jpa_df,jpb_df,jpl_df,2)/scaling_factor;
dEz = Pola_SI(jpa_df,jpb_df,jpl_df,3)/scaling_factor;
dBx = Pola_SI(jpa_df,jpb_df,jpl_df,4)/scaling_factor;
dBy = Pola_SI(jpa_df,jpb_df,jpl_df,5)/scaling_factor;
dBz = Pola_SI(jpa_df,jpb_df,jpl_df,6)/scaling_factor;

% Start Time series:
maxf=-999.;
minf=999.;
maxdfim=-999.;
mindfim=999.;

j=s_df;

for timerun=0:num_periods*timesteps

    if (periods)
        time=2.d0*pi*(1.d0*timerun)/(abs(real(x))*1.d0*timesteps);
    else
        time=2.d0*pi*(1.d0*timerun)/(real(1.d0)*1.d0*timesteps);
    end
	disp(['time check',num2str(time)])

    filename=['pdrk_dist_deltaf_(kdi=' num2str(k*cwp,'%.3f') ...
        ')(ampl=' num2str(ampl,'%.2f') ')_' num2str(timerun,'%03d') '.vtk'];
    filename2=['pdrk_normEB_(kdi=' num2str(k*cwp,'%.3f') ...
        ')(ampl=' num2str(ampl,'%.2f') ')_' num2str(timerun,'%03d') '.vtk'];
    disp(['kdi=' num2str(k*cwp,'%.3f') ';' ...
        'wr/wci=' num2str(real(x/wcs1),'%.3f') ';' ...
        'wi/wci=' num2str(imag(x/wcs1),'%.3f')]);

    % Here should the loop over all v start:
    for ii=0:vxsteps
    for kk=0:vysteps
    for ll=0:vzsteps

    vx=vxrange(1)+(1.d0*ii)*(vxrange(2)-vxrange(1))/(1.d0*vxsteps); % ¡¾vA¡¿
    vy=vyrange(1)+(1.d0*kk)*(vyrange(2)-vyrange(1))/(1.d0*vysteps); % ¡¾vA¡¿
    vz=vzrange(1)+(1.d0*ll)*(vzrange(2)-vzrange(1))/(1.d0*vzsteps); % ¡¾vA¡¿
    vx=vx*vA; vy=vy*vA; vz=vz*vA; % ¡¾m¡¿

    vpar=vz; % ¡¾m¡¿
    vperp=sqrt(vx*vx+vy*vy); % ¡¾m¡¿
    phi=acos(vx/vperp);
    if vy<0.d0
        phi=2.d0*pi-phi; end
    if vperp==0.d0
        phi=0.d0; end

    dfdvpar=-2.d0*(vpar-vds(j))*exp(-vperp^2/vthpS(j)^2);
    dfdvpar=dfdvpar*exp(-(vpar-vds(j))^2/vthzS(j)^2);
    dfdvpar=dfdvpar/(pi^1.5*vthpS(j)^2*vthzS(j)^3)*ns0(j);
                     
    dfdvperp=-2.d0*vperp*exp(-vperp^2/vthpS(j)^2);
    dfdvperp=dfdvperp*exp(-(vpar-vds(j))^2/vthzS(j)^2);
    dfdvperp=dfdvperp/(pi^1.5*vthpS(j)^4*vthzS(j))*ns0(j);

    fnull=1.d0/(vthpS(j)^2*vthzS(j)*pi^1.5);
    fnull=fnull*exp(-(vpar-vds(j))^2/vthzS(j)^2)*exp(-vperp^2/vthpS(j)^2)*ns0(j);

    z=kx*vperp/wcs(j);

    %Summation over all Bessel functions:
    deltaf=0.d0;% dfqE=0.d0;dfqvB=0.d0;
    for m=-N:N

        a=(1.d0*m)*wcs(j)-x+kz*vpar;

        UStix=dfdvperp+kz*(vperp*dfdvpar-vpar*dfdvperp)/x;
        VStix=kx*(vperp*dfdvpar-vpar*dfdvperp)/x;
%         qEStix = dfdvperp;
%         qvBStix1 = kz*(vperp*dfdvpar-vpar*dfdvperp)/x;

        comp1=dEx*UStix*(a*1i*cos(phi)-wcs(j)*sin(phi))/(wcs(j)^2-a^2);
        comp2=dEy*UStix*(a*1i*sin(phi)+wcs(j)*cos(phi))/(wcs(j)^2-a^2);
        comp3=-1i*dEz*dfdvpar/a;
        comp4=-dEz*VStix*(a*1i*cos(phi)-wcs(j)*sin(phi))/(wcs(j)^2-a^2);
        
%         comp11=dEx*qEStix*(a*1i*cos(phi)-wcs(j)*sin(phi))/(wcs(j)^2-a^2);
%         comp12=dEy*qEStix*(a*1i*sin(phi)+wcs(j)*cos(phi))/(wcs(j)^2-a^2);
%         comp13=comp3;
%         
%         comp21=dEx*qvBStix1*(a*1i*cos(phi)-wcs(j)*sin(phi))/(wcs(j)^2-a^2);
%         comp22=dEy*qvBStix1*(a*1i*sin(phi)+wcs(j)*cos(phi))/(wcs(j)^2-a^2);
%         comp23=comp4;
        

        deltaf = deltaf - qs(j)/ms(j)*exp(1i*z*sin(phi))*exp(-1i*(1.d0*m)*phi)*besselj(m,z)*ampl*(comp1+comp2+comp3+comp4);
%         dfqE = dfqE - qs(j)/ms(j)*exp(1i*z*sin(phi))*exp(-1i*(1.d0*m)*phi)*besselj(m,z)*ampl*(comp11+comp12+comp13);
%         dfqvB = dfqvB - qs(j)/ms(j)*exp(1i*z*sin(phi))*exp(-1i*(1.d0*m)*phi)*besselj(m,z)*ampl*(comp21+comp22+comp23);
    end % end m loop

    if const_r
        if damping
            deltaf=deltaf*exp(-1i*time*x);
%             dfqE=dfqE*exp(-1i*time*x);
%             dfqvB=dfqvB*exp(-1i*time*x);
        else
            deltaf=deltaf*exp(-1i*time*real(x));
%             dfqE=dfqE*exp(-1i*time*real(x));
%             dfqvB=dfqvB*exp(-1i*time*real(x));
        end
    else 
        deltaf=deltaf*exp(-1i*2.d0*pi*(1.d0*timerun)/(1.d0*timesteps));
%         dfqE=dfqE*exp(-1i*2.d0*pi*(1.d0*timerun)/(1.d0*timesteps));
%         dfqvB=dfqvB*exp(-1i*2.d0*pi*(1.d0*timerun)/(1.d0*timesteps));
    end

    maxf=max((maxf),real(fnull+deltaf));
    minf=min((minf),real(fnull+deltaf));
    maxdfim=max((maxdfim),imag(deltaf));
    mindfim=min((mindfim),imag(deltaf));

	VXarray(ii+1,kk+1,ll+1)=real(vx);
	VYarray(ii+1,kk+1,ll+1)=real(vy);
	VZarray(ii+1,kk+1,ll+1)=real(vz);
    f0array(ii+1,kk+1,ll+1)=real(fnull);
	farray(ii+1,kk+1,ll+1)=real(fnull+deltaf);
    imfarray(ii+1,kk+1,ll+1)=imag(fnull+deltaf);
    dfarray(ii+1,kk+1,ll+1)=real(deltaf);
%     dfqEarray(ii+1,kk+1,ll+1)=real(dfqE);
%     dfqvBarray(ii+1,kk+1,ll+1)=real(dfqvB);
    end % end vz loop
    end % end vy loop
    end % end vx loop
    
    disp(['# Maximum value of f:',num2str(maxf)]);
    disp(['# Minimum value of f:',num2str(minf)]);
    disp(['# Maximum value of Im(delta f):',num2str(maxdfim)]);
    disp(['# Minimum value of Im(delta f):',num2str(mindfim)]);

    title1 = 'f0+deltaf';
    title2 = 'deltaf';
%     title3 = 'dfqE';
%     title4 = 'dfqvB';

%% save to vtk file to be visulized by Paraview
    E_direc = vxrange(2)*real([dEx,dEy,dEz]*exp(-1i*time*real(x))/norm(real([dEx,dEy,dEz])));
    B_direc = vxrange(2)*real([dBx,dBy,dBz]*exp(-1i*time*real(x))/norm(real([dBx,dBy,dBz])));
    filename = [savepath filename];
    filename2 = [savepath filename2];
    vtkwrite(filename,'structured_grid',VXarray/vA,VYarray/vA,VZarray/vA,...
        'scalars',title1,farray,'scalars',title2,dfarray,'Precision',20,'BINARY');
    %%% write E field and B field
    fid = fopen(filename2, 'w'); 
    % VTK files contain five major parts
    % 1. VTK DataFile Version
    fprintf(fid, '# vtk DataFile Version 2.0\n');
    % 2. Title
    fprintf(fid, 'VTK from Matlab\n');
    fprintf(fid, 'ASCII\n');
    fprintf(fid, 'DATASET POLYDATA\n');
    fprintf(fid, ['POINTS ' num2str(3) ' float\n']);
    spec = [repmat(['%0.', '15', 'f '], 1, 3), '\n'];   
    fprintf(fid, spec, [0.0,0.0,0.0]);
    fprintf(fid, spec, E_direc);
    fprintf(fid, spec, B_direc);
    fprintf(fid,'\nLINES %d %d\n',2,6);
    fprintf(fid,'2 %d %d\n',[0;1]);
    fprintf(fid,'2 %d %d\n',[0;2]);
    fprintf(fid,'\nCELL_DATA 2\n');
    fprintf(fid,'SCALARS cell_scalars float 1\n');
    fprintf(fid,'LOOKUP_TABLE my_table\n');
    fprintf(fid,'0.0\n');
    fprintf(fid,'1.0\n');
    fprintf(fid,'\nLOOKUP_TABLE my_table 2\n');
    fprintf(fid,'0.0 1.0 0.0 1.0\n');
    fprintf(fid,'1.0 1.0 0.0 1.0');
    fclose(fid);

end