clear all;
clc;



%extracting data from experiments
flm=xlsread('Counter_combined'); %opposing effect of flow and TGF-\beta gradient
flp=xlsread('Parallel_combined'); %additive effect of flow and TGF-\beta gradient
flT=xlsread('dT'); %TGF-\beta gradient only
flF=xlsread('Flow'); %Flow only
flma=xlsread('Counter_above'); %opposing effect of flow and TGF-\beta gradient above detection limit
flpa=xlsread('Parallel_above'); %additive effect of flow and TGF-\beta gradient above detection limit
flmb=xlsread('Counter_below'); %opposing effect of flow and TGF-\beta gradient below detection limit
flpb=xlsread('Parallel_below'); %additive effect of flow and TGF-\beta gradient below detection limit

%extract angles
%Experimental data - TGF-\beta gradient only
for i=1:size(flT,1)
   if flT(i)>0 
   flT(i)=flT(i)*2*pi/360; 
   else
   flT(i)=(flT(i)+360)*2*pi/360;     
   end
   CIT(i)=cos(flT(i));
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

epsilonmb=0.004; %epsilonT for region 1
cmb=0.05; %TGF conc for region 1

epsilonma = 0.03;
cma = 0.4;

epsilonpa = 0.02;
cpa = 9;

epsilonpb = 0.005;
cpb = 10;

cT=5;
epsilonT=0.08;


ntrials=1e4; %Number of parameter sets over which search is made

for k1=1:ntrials

eta=4000*rand; %uniformly searching for \eta between 0-4000
phi1=0.01*rand;  %uniformly searching for \phi_1 between 0-0.01
phi2=4*rand; %uniformly searching for \phi_2 between 0-4
beta=rand; %uniformly searching for \beta between 0-1

%These ranges have been decided after doing randomized searches over larger
%ranges

%Saving the parameters as arrays to extract later

xe(k1)=eta;
xp1(k1)=phi1;
xp2(k1)=phi2;
xb(k1)=beta;
    

%calculation of \Delta m or \kappa which decides sharpness of p(\theta)
kappaT=eta*epsilonT*cT*beta/(1+cT*beta)^2;
kappaF=eta*phi1/(1+phi2)^2;
kappap=eta*(epsilonp*cp*beta+phi1)/(1+cp*beta+phi2)^2;
kappam=eta*(epsilonm*cm*beta-phi1)/(1+cm*beta+phi2)^2;

kappapa=eta*(epsilonpa*cpa*beta+phi1)/(1+cpa*beta+phi2)^2;
kappama=eta*(epsilonma*cma*beta-phi1)/(1+cma*beta+phi2)^2;

kappapb=eta*(epsilonpb*cpb*beta+phi1)/(1+cpb*beta+phi2)^2;
kappamb=eta*(epsilonmb*cmb*beta-phi1)/(1+cmb*beta+phi2)^2;

nbins=100;
dtheta=2*pi/nbins;

%Flow only
for j=1:nbins
    theta(j)=2*pi*j/nbins;
    panaF(j)=(1-alpha)/(2*pi)+alpha*exp(kappaF*cos(theta(j)))/(2*pi*besseli(0,kappaF));
    cumanaF(j)=0;
    
    if j==1
       cumanaF(j)=panaF(j)*dtheta;
    else
       cumanaF(j)=cumanaF(j-1)+panaF(j)*dtheta;
    end

end

npts=1e4;

for k=1:npts
   z=rand;
   for j=1:nbins
      if j==1
         if z>0 && z<cumanaF(j)
            thetFana(k)=theta(j);
            CIanaF(k)=cos(thetFana(k));
         end
      else
         if z>cumanaF(j-1) && z<cumanaF(j)
            thetFana(k)=theta(j);
            CIanaF(k)=cos(thetFana(k));
         end

         
      end
       
   end
    
end


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

npts=1e4;

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



%TGF-\beta grad+Flow (additive effect)
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

npts=1e4;

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



%TGF-\beta grad-Flow (opposing effect)
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

npts=1e4;

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





%TGF-\beta grad-Flow above limit
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

npts=1e4;

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



%TGF-\beta grad-Flow below limit
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

npts=1e4;

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

%TGF-\beta grad+Flow above limit
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

npts=1e4;

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



%TGF-\beta grad+Flow below limit
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

npts=1e4;

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

disp(k1)

error(k1)=abs(median(CIT)-median(CIanaT))+abs(median(CIF)-median(CIanaF))+abs(median(CIp)-median(CIanap))+abs(median(CIm)-median(CIanam));

%storing median values in arrays for plot
CIpltT(k1)=median(CIanaT);
CIpltF(k1)=median(CIanaF);
CIpltp(k1)=median(CIanap);
CIpltm(k1)=median(CIanam);

CIpltpa(k1)=median(CIanapa);
CIpltma(k1)=median(CIanama);

CIpltpb(k1)=median(CIanapb);
CIpltmb(k1)=median(CIanamb);



end

%%
prp=[0.7 0 0.7];
or=[1 0.5 0];

[M,A]=min(error); %extracting minimal error
%Writing best parameters
s=strcat("Best parameters : \alpha = ",num2str(alpha),", \eta = ",num2str(xe(A)),", \phi_1 = ",num2str(xp1(A)),...
    ", \phi_2 = ",num2str(xp2(A)),", \beta = ",num2str(xb(A)));
disp(s)
%The best parameters will be slightly different each time because of the
%randomized search. But it will still be very close to the values used for
%plotting the figure in the paper

%plotting the best fitted data
figure(1); clf
h1=bar(0,median(CIT),'BarWidth', 0.3,'FaceColor','r', 'LineWidth', 3.5); hold on;
h2=bar(0.4,CIpltT(A),'BarWidth', 0.3,'FaceColor','w','EdgeColor','r','LineWidth', 3.5);
h3=bar(1,median(CIF),'BarWidth', 0.3,'FaceColor','b', 'LineWidth', 3.5); hold on;
h4=bar(1.4,CIpltF(A),'BarWidth', 0.3,'FaceColor','w','EdgeColor','b','LineWidth', 3.5);
h5=bar(2,median(CIp),'BarWidth', 0.3,'FaceColor',prp, 'LineWidth', 3.5); hold on;
h6=bar(2.4,CIpltp(A),'BarWidth', 0.3,'FaceColor','w','EdgeColor',prp,'LineWidth', 3.5);
h7=bar(3,median(CIm),'BarWidth', 0.3,'FaceColor',or, 'LineWidth', 3.5); hold on;
h8=bar(3.4,CIpltm(A),'BarWidth', 0.3,'FaceColor','w','EdgeColor',or,'LineWidth', 3.5);
%xlabel('condition')
%s=strcat('\Delta m_F^*=',num2str(DelmF),', \mu^*=',num2str(mu),', \nu^*=',num2str(nu),', \phi^*=',num2str(phi));
title('Best fit')
ylabel('DAI')
legend({'Experiment','Model'},'Location','southwest')
   set(gca, 'LineWidth', 3.5)
set(gca,'fontsize',20);
xticks([0.2 1.2 2.2 3.2])
xticklabels({'\nabla T','IF','\nabla T+IF','\nabla T-IF'})
xtickangle(0)
ylim([-0.3 0.65]) 
xlim([-0.5 4]) 
set(gcf, 'PaperSize', [4 2]);
saveas(h1,'Fig_5C_best_fit.png')

%plotting the corresponding predictions

figure(2); clf
h1=bar(2,median(CIpa),'BarWidth', 0.3,'FaceColor',prp, 'LineWidth', 3.5); hold on;
h6=bar(2.4,CIpltpa(A),'BarWidth', 0.3,'FaceColor','w','EdgeColor',prp,'LineWidth', 3.5);
h7=bar(3,median(CIpb),'BarWidth', 0.3,'FaceColor',prp, 'LineWidth', 3.5); hold on;
h8=bar(3.4,CIpltpb(A),'BarWidth', 0.3,'FaceColor','w','EdgeColor',prp,'LineWidth', 3.5);
h9=bar(4,median(CIma),'BarWidth', 0.3,'FaceColor',or, 'LineWidth', 3.5); hold on;
h10=bar(4.4,CIpltma(A),'BarWidth', 0.3,'FaceColor','w','EdgeColor',or,'LineWidth', 3.5);
h11=bar(5,median(CImb),'BarWidth', 0.3,'FaceColor',or, 'LineWidth', 3.5); hold on;
h12=bar(5.4,CIpltmb(A),'BarWidth', 0.3,'FaceColor','w','EdgeColor',or,'LineWidth', 3.5);
%xlabel('condition')
%s=strcat('\Delta m_F^*=',num2str(DelmF),', \mu^*=',num2str(mu),', \nu^*=',num2str(nu),', \phi^*=',num2str(phi));
title('Prediction')
ylabel('DAI')
xlabel('Signal State : \nabla T/IF')
legend({'Experiment','Model'},'Location','southwest')
   set(gca, 'LineWidth', 3.5)
set(gca,'fontsize',20);
ylim([-0.3 0.65])
xticks([0.2 1.2 2.2 3.2 4.2 5.2])
xticklabels({'\nabla T','V_f','+/+','0/+','+/-','0/-'})
xtickangle(0)
xlim([1.5 6])
saveas(h1,'Fig_5C_best_prediction.png')





