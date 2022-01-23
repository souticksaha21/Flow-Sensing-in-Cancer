clear all;
clc;


%parameters for best possible fits
eta=3.0039e+03;
phi1=0.0018;
phi2=1.6612;
beta=0.9073;

%reading experimental data
flm=xlsread('Counter_combined');
flp=xlsread('Parallel_combined');
flT=xlsread('dT');
flF=xlsread('Flow');
flma=xlsread('Counter_above');
flpa=xlsread('Parallel_above');
flmb=xlsread('Counter_below');
flpb=xlsread('Parallel_below');



%extract angles
%Experimental data - TGF-\beta gradient only
for i=1:size(flT,1)
   if flT(i)>0 
   flT(i)=flT(i)*2*pi/360; %convert angles to radians
   else
   flT(i)=(flT(i)+360)*2*pi/360;     
   end
   CIT(i)=cos(flT(i)); %get CI values from angles
end

%Experimental data - flow only
for i=1:size(flF,1)
   if flF(i)>0 
   flF(i)=flF(i)*2*pi/360; 
   else
   flF(i)=(flF(i)+360)*2*pi/360;     
   end
   CIF(i)=cos(flF(i));
end

%Experimental data - effect of TGF-\beta gradient and flow in opposing
%directions
for i=1:size(flm,1)
   if flm(i)>0 
   flm(i)=flm(i)*2*pi/360; 
   else
   flm(i)=(flm(i)+360)*2*pi/360;     
   end
   CIm(i)=cos(flm(i));
end


%Experimental data - effect of TGF-\beta gradient and flow in same
%directions
for i=1:size(flp,1)
   if flp(i)>0 
   flp(i)=flp(i)*2*pi/360; 
   else
   flp(i)=(flp(i)+360)*2*pi/360;     
   end
   CIp(i)=cos(flp(i));
end

%Detemining \alpha from the highest mean CI from the experimental datasets
%used to fit
alpha=max([mean(CIT),mean(CIF),mean(CIm),mean(CIp)]);

%Experimental data - effect of TGF-\beta gradient and flow in opposing
%directions above detection limit

for i=1:size(flma,1)
   if flma(i)>0 
   flma(i)=flma(i)*2*pi/360; 
   else
   flma(i)=(flma(i)+360)*2*pi/360;     
   end
   CIma(i)=cos(flma(i));
end

%Experimental data - effect of TGF-\beta gradient and flow in same
%directions above detection limit

for i=1:size(flpa,1)
   if flpa(i)>0 
   flpa(i)=flpa(i)*2*pi/360; 
   else
   flpa(i)=(flpa(i)+360)*2*pi/360;     
   end
   CIpa(i)=cos(flpa(i));
end

%Experimental data - effect of TGF-\beta gradient and flow in opposing
%directions below detection limit

for i=1:size(flmb,1)
   if flmb(i)>0 
   flmb(i)=flmb(i)*2*pi/360; 
   else
   flmb(i)=(flmb(i)+360)*2*pi/360;     
   end
   CImb(i)=cos(flmb(i));
end

%Experimental data - effect of TGF-\beta gradient and flow in same
%directions below detection limit

for i=1:size(flpb,1)
   if flpb(i)>0 
   flpb(i)=flpb(i)*2*pi/360; 
   else
   flpb(i)=(flpb(i)+360)*2*pi/360;     
   end
   CIpb(i)=cos(flpb(i));
end

%Different explisons and c0's based on experimental condition

%epsilon=ga'/c0
epsilonp=0.01;
cp=9.5;

epsilonm=0.015;
cm=0.2;

epsilonmb=0.004; 
cmb=0.05; 

epsilonma = 0.03;
cma = 0.4;

epsilonpa = 0.02;
cpa = 9;

epsilonpb = 0.005;
cpb = 10;

cT=5;
epsilonT=0.08;


%Getting the \Delta m's (which is called \kappa here) for different conditions
kappaT=eta*epsilonT*cT*beta/(1+cT*beta)^2; %only TGF-\beta gradient present
kappaF=eta*phi1/(1+phi2)^2; %only flow present
kappap=eta*(epsilonp*cp*beta+phi1)/(1+cp*beta+phi2)^2; %Additive effect of flow and TGF-\beta gradient
kappam=eta*(epsilonm*cm*beta-phi1)/(1+cm*beta+phi2)^2; %Opposing effect of flow and TGF-\beta gradient

kappapa=eta*(epsilonpa*cpa*beta+phi1)/(1+cpa*beta+phi2)^2; %Additive effect of flow and TGF-\beta gradient 
%above limits
kappama=eta*(epsilonma*cma*beta-phi1)/(1+cma*beta+phi2)^2; %Opposing effect of flow and TGF-\beta gradient 
%above limit

kappapb=eta*(epsilonpb*cpb*beta+phi1)/(1+cpb*beta+phi2)^2; %Additive effect of flow and TGF-\beta gradient 
%below limits
kappamb=eta*(epsilonmb*cmb*beta-phi1)/(1+cmb*beta+phi2)^2; %Opposing effect of flow and TGF-\beta gradient 
%below limit

nbins=1000;
dtheta=2*pi/nbins;

%Model data extraction
%Flow only
for j=1:nbins
    theta(j)=2*pi*j/nbins; %Uniformly distributed angle
    %Defining probability distribution numerically
    panaF(j)=(1-alpha)/(2*pi)+alpha*exp(kappaF*cos(theta(j)))/(2*pi*besseli(0,kappaF)); 
    %Defining the cumulative of the distrbution
    cumanaF(j)=0;
    
    if j==1
       cumanaF(j)=panaF(j)*dtheta;
    else
       cumanaF(j)=cumanaF(j-1)+panaF(j)*dtheta;
    end

end

npts=1e5;

for k=1:npts
   z=rand;
   for j=1:nbins
      if j==1
         if z>0 && z<cumanaF(j)
            thetFana(k)=theta(j); %Extracting angles from the distribution
            CIanaF(k)=cos(thetFana(k)); %Extracting corresponding CI values
         end
      else
         if z>cumanaF(j-1) && z<cumanaF(j)
            thetFana(k)=theta(j);
            CIanaF(k)=cos(thetFana(k));
         end

         
      end
       
   end
    
end


sF=strcat("DAI(IF)=",num2str(median(CIanaF)));
disp(sF) %Median of CI in model from flow only case

%TGF-\beta grad only
for j=1:nbins
    theta(j)=2*pi*j/nbins;
    panaT(j)=(1-alpha)/(2*pi)+alpha*exp(kappaT*cos(theta(j)))/(2*pi*besseli(0,kappaT));
    cumanaT(j)=0;
    
    if j==1
       cumanaT(j)=panaT(j)*dtheta;
    else
       cumanaT(j)=cumanaT(j-1)+panaT(j)*dtheta;
    end

end

npts=1e5;

for k=1:npts
   z=rand;
   for j=1:nbins
      if j==1
         if z>0 && z<cumanaT(j)
            thetTana(k)=theta(j);
            CIanaT(k)=cos(thetTana(k));
         end
      else
         if z>cumanaT(j-1) && z<cumanaT(j)
            thetTana(k)=theta(j);
            CIanaT(k)=cos(thetTana(k));
         end

         
      end
       
   end
    
end

sT=strcat("DAI(\nabla T)=",num2str(median(CIanaT)));
disp(sT) %Median of CI in model from TGF-\beta gradient only case


%TGF-\beta grad+Flow only
for j=1:nbins
    theta(j)=2*pi*j/nbins;
    panap(j)=(1-alpha)/(2*pi)+alpha*exp(kappap*cos(theta(j)))/(2*pi*besseli(0,kappap));
    cumanap(j)=0;
    
    if j==1
       cumanap(j)=panap(j)*dtheta;
    else
       cumanap(j)=cumanap(j-1)+panap(j)*dtheta;
    end

end

npts=1e5;

for k=1:npts
   z=rand;
   for j=1:nbins
      if j==1
         if z>0 && z<cumanap(j)
            thetpana(k)=theta(j);
            CIanap(k)=cos(thetpana(k));
         end
      else
         if z>cumanap(j-1) && z<cumanap(j)
            thetpana(k)=theta(j);
            CIanap(k)=cos(thetpana(k));
         end

         
      end
       
   end
    
end

sp=strcat("DAI(\nabla T+IF)=",num2str(median(CIanap)));
disp(sp) %Median of CI in model from additive effect of flow and TGF-\beta gradient case


%TGF-\beta grad-Flow only
for j=1:nbins
    theta(j)=2*pi*j/nbins;
    panam(j)=(1-alpha)/(2*pi)+alpha*exp(kappam*cos(theta(j)))/(2*pi*besseli(0,kappam));
    cumanam(j)=0;
    
    if j==1
       cumanam(j)=panam(j)*dtheta;
    else
       cumanam(j)=cumanam(j-1)+panam(j)*dtheta;
    end

end

npts=1e5;

for k=1:npts
   z=rand;
   for j=1:nbins
      if j==1
         if z>0 && z<cumanam(j)
            thetmana(k)=theta(j);
            CIanam(k)=cos(thetmana(k));
         end
      else
         if z>cumanam(j-1) && z<cumanam(j)
            thetmana(k)=theta(j);
            CIanam(k)=cos(thetmana(k));
         end

         
      end
       
   end
    
end

sm=strcat("DAI(\nabla T-IF)=",num2str(median(CIanam)));
disp(sm) %Median of CI in model from opposing effect of flow and TGF-\beta gradient case





%TGF-\beta grad-Flow only above limit
for j=1:nbins
    theta(j)=2*pi*j/nbins;
    panama(j)=(1-alpha)/(2*pi)+alpha*exp(kappama*cos(theta(j)))/(2*pi*besseli(0,kappama));
    cumanama(j)=0;
    
    if j==1
       cumanama(j)=panama(j)*dtheta;
    else
       cumanama(j)=cumanama(j-1)+panama(j)*dtheta;
    end

end

npts=1e5;

for k=1:npts
   z=rand;
   for j=1:nbins
      if j==1
         if z>0 && z<cumanama(j)
            thetmaana(k)=theta(j);
            CIanama(k)=cos(thetmaana(k));
         end
      else
         if z>cumanama(j-1) && z<cumanama(j)
            thetmaana(k)=theta(j);
            CIanama(k)=cos(thetmaana(k));
         end

         
      end
       
   end
    
end

sma=strcat("DAI(+/-)=",num2str(median(CIanama)));
disp(sma) %Median of CI in model from opposing effect of flow and TGF-\beta gradient case above limit


%TGF-\beta grad-Flow only below limit
for j=1:nbins
    theta(j)=2*pi*j/nbins;
    panamb(j)=(1-alpha)/(2*pi)+alpha*exp(kappamb*cos(theta(j)))/(2*pi*besseli(0,kappamb));
    cumanamb(j)=0;
    
    if j==1
       cumanamb(j)=panamb(j)*dtheta;
    else
       cumanamb(j)=cumanamb(j-1)+panamb(j)*dtheta;
    end

end

npts=1e5;

for k=1:npts
   z=rand;
   for j=1:nbins
      if j==1
         if z>0 && z<cumanamb(j)
            thetmbana(k)=theta(j);
            CIanamb(k)=cos(thetmbana(k));
         end
      else
         if z>cumanamb(j-1) && z<cumanamb(j)
            thetmbana(k)=theta(j);
            CIanamb(k)=cos(thetmbana(k));
         end

         
      end
       
   end
    
end

smb=strcat("DAI(0/-)=",num2str(median(CIanamb)));
disp(smb) %Median of CI in model from opposing effect of flow and TGF-\beta gradient case below limit



%TGF-\beta grad+Flow only above limit
for j=1:nbins
    theta(j)=2*pi*j/nbins;
    panapa(j)=(1-alpha)/(2*pi)+alpha*exp(kappapa*cos(theta(j)))/(2*pi*besseli(0,kappapa));
    cumanapa(j)=0;
    
    if j==1
       cumanapa(j)=panapa(j)*dtheta;
    else
       cumanapa(j)=cumanapa(j-1)+panapa(j)*dtheta;
    end

end

npts=1e5;

for k=1:npts
   z=rand;
   for j=1:nbins
      if j==1
         if z>0 && z<cumanapa(j)
            thetpaana(k)=theta(j);
            CIanapa(k)=cos(thetpaana(k));
         end
      else
         if z>cumanapa(j-1) && z<cumanapa(j)
            thetpaana(k)=theta(j);
            CIanapa(k)=cos(thetpaana(k));
         end

         
      end
       
   end
    
end

spa=strcat("DAI(+/+)=",num2str(median(CIanapa)));
disp(spa) %Median of CI in model from additive effect of flow and TGF-\beta gradient case above limit



%TGF-\beta grad+Flow only below limit
for j=1:nbins
    theta(j)=2*pi*j/nbins;
    panapb(j)=(1-alpha)/(2*pi)+alpha*exp(kappapb*cos(theta(j)))/(2*pi*besseli(0,kappapb));
    cumanapb(j)=0;
    
    if j==1
       cumanapb(j)=panapb(j)*dtheta;
    else
       cumanapb(j)=cumanapb(j-1)+panapb(j)*dtheta;
    end

end

npts=1e5;

for k=1:npts
   z=rand;
   for j=1:nbins
      if j==1
         if z>0 && z<cumanapb(j)
            thetpbana(k)=theta(j);
            CIanapb(k)=cos(thetpbana(k));
         end
      else
         if z>cumanapb(j-1) && z<cumanapb(j)
            thetpbana(k)=theta(j);
            CIanapb(k)=cos(thetpbana(k));
         end

         
      end
       
   end
    
end

spb=strcat("DAI(0/+)=",num2str(median(CIanapb)));
disp(spb) %Median of CI in model from additive effect of flow and TGF-\beta gradient case below limit


%%
prp=[0.7 0 0.7];
or=[1 0.5 0];
%Median plot for fit
figure(1); clf
h1=bar(0,median(CIT),'BarWidth', 0.3,'FaceColor','r', 'LineWidth', 3.5); hold on;
h2=bar(0.4,median(CIanaT),'BarWidth', 0.3,'FaceColor','w','EdgeColor','r','LineWidth', 3.5);
h3=bar(1,median(CIF),'BarWidth', 0.3,'FaceColor','b', 'LineWidth', 3.5); hold on;
h4=bar(1.4,median(CIanaF),'BarWidth', 0.3,'FaceColor','w','EdgeColor','b','LineWidth', 3.5);
h5=bar(2,median(CIp),'BarWidth', 0.3,'FaceColor',prp, 'LineWidth', 3.5); hold on;
h6=bar(2.4,median(CIanap),'BarWidth', 0.3,'FaceColor','w','EdgeColor',prp,'LineWidth', 3.5);
h7=bar(3,median(CIm),'BarWidth', 0.3,'FaceColor',or, 'LineWidth', 3.5); hold on;
h8=bar(3.4,median(CIanam),'BarWidth', 0.3,'FaceColor','w','EdgeColor',or,'LineWidth', 3.5);
title('Data fit for model')
ylabel('DAI')
legend({'Experiment','Model'}, 'fontsize', 20,'Location','southwest')
   set(gca, 'LineWidth', 3.5)
set(gca,'fontsize',20);
ylim([-0.35 0.7])
xticks([0.2 1.2 2.2 3.2])
%xticklabels({'\nabla T','V_f','Parallel V_f','Counter V_f'})
xticklabels({'','','',''})
xtickangle(45)
yticks([-0.2 0 0.2 0.4 0.6])
yticklabels({'-0.2','0.0','0.2','0.4','0.6'})
set(gca, 'fontweight', 'bold', 'fontsize', 25)
xlim([-0.5 4]) 
set(gcf, 'PaperSize', [4 2]);
saveas(h1,'Fig_5C_part1.png')

%Median plot for prediction
figure(2); clf
h1=bar(2,median(CIpa),'BarWidth', 0.3,'FaceColor',prp, 'LineWidth', 3.5); hold on;
h6=bar(2.4,median(CIanapa),'BarWidth', 0.3,'FaceColor','w','EdgeColor',prp,'LineWidth', 3.5);
h7=bar(3,median(CIpb),'BarWidth', 0.3,'FaceColor',prp, 'LineWidth', 3.5); hold on;
h8=bar(3.4,median(CIanapb),'BarWidth', 0.3,'FaceColor','w','EdgeColor',prp,'LineWidth', 3.5);
h9=bar(4,median(CIma),'BarWidth', 0.3,'FaceColor',or, 'LineWidth', 3.5); hold on;
h10=bar(4.4,median(CIanama),'BarWidth', 0.3,'FaceColor','w','EdgeColor',or,'LineWidth', 3.5);
h11=bar(5,median(CImb),'BarWidth', 0.3,'FaceColor',or, 'LineWidth', 3.5); hold on;
h12=bar(5.4,median(CIanamb),'BarWidth', 0.3,'FaceColor','w','EdgeColor',or,'LineWidth', 3.5);
title('Prediction', 'fontsize', 20)
ylabel('DAI', 'fontsize', 25)
legend({'Experiment','Model'}, 'fontsize', 20,'Location','southwest')
   set(gca, 'LineWidth', 3.5)
ylim([-0.35 0.7])
xticks([0.2 1.2 2.2 3.2 4.2 5.2])

%xticklabels({'\nabla T','V_f','Parallel V_f+>limit','Parallel V_f+<limit','Counter V_f+>limit','Counter V_f+<limit'})
xticklabels({'','','','','',''})
yticks([-0.2 0 0.2 0.4 0.6])
yticklabels({'-0.2','0.0','0.2','0.4','0.6'})
set(gca, 'fontweight', 'bold', 'fontsize', 25)
xtickangle(45)
%set(gca,'FontName', 'Arial');
xlim([1.5 6])
saveas(h1,'Fig_5C_part2.png')

