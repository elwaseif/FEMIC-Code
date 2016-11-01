%% Function to support FEMIC inversion parameter display
% this function may be called on each iteration of the model to display the
% model error, model vs. observation comparison, and model parameters
% profile at a selected station location.
%
% G. Schultz 2008

function [hplot] = FEMIC_plotINVinterim(gmodel,obs,f,iter,max_iter,p,mdlD,err,p0,pmin,pmax)

%% Postion the figure window in the screen
p0=1./10.^p0;p=1./10.^p;
Screen.size = get(0, 'ScreenSize');  %Determines the Size of the Screen
%Position                                  %Set for a screen of 1024x768 and to set above windows taskbar
    Screen.Left = Screen.size(1);
    Screen.Bottom = Screen.size(2);
    Screen.Width = Screen.size(3);
    Screen.Height = Screen.size(4);

hplot=figure(12);

set(gcf,'NumberTitle','off',...
            'Name','FEMIC Frequency-domain EMI Inversion',...
            'Resize','on');
%% Plots the forward model values and observations at each iteration
hold on
set(gcf,'Position',[floor((Screen.Width-610)/1.5) floor((Screen.Height-610)/2) 600 600]);
subplot(211); 
plot(f,obs,'k*',f,gmodel,'rd-','MarkerSize',6+1.5*iter);
df = (max(f)-min(f))/10;
set(gca,'XLim',[min(f)-df max(f)+df]);
set(gca,'YLim',[min([gmodel;obs])-0.2*min([gmodel;obs]) max([gmodel;obs])+0.2*max([gmodel;obs])]);
xlabel('Frequency (Hz)');ylabel('Quadrature'); legend('Observed','Predicted'); 
title(['Iteration ',num2str(iter),' - RMSerr= ',num2str(err)]);
set(gca,'FontSize',12)
%% Plots the RMS residual errors between model and observations
hold on
subplot(223);
plot(iter,log10(err),'*')
set(gca,'XLim',[0 max_iter]);
set(gca,'YLim',[-2 2]);
xlabel('Iteration'); ylabel('log_1_0(RMS Error)');
set(gca,'FontSize',12)

%% Sorts the model parameters vector for plotting
D=[mdlD];
D(end) = D(end-1);
d(1,:) = [0 1]; 
eca(1,:) = [p(1) p(1)];
ip(1,:) =[p0(1) p0(1)];
for i=2:length(D), 
    d(i,:) = [d(i-1,2) d(i-1,2)+D(i)];
    eca(i,:)=[p(i) p(i)]; 
    ip(i,:) = [p0(i) p0(i)];
end
dd = reshape(d',2*length(d),1);
eeca = reshape(eca',2*length(d),1);
oeca = reshape(ip',2*length(d),1);

%% Plots the current model
hold on;
subplot(224);
plot(oeca,-dd,'r.-','LineWidth',2);
hold on
plot(eeca,-dd,'k-','LineWidth',iter);
hold on
plot([pmin;pmin],[0;min(-dd)],'g-',[pmax;pmax],[0;min(-dd)],'g-',...
    'LineWidth',2);
set(gca,'YLim',[min(-dd) 0]);
set(gca,'XLim',[-2 1100]);
xlabel('Resist. (Ohm m)'); ylabel('Depth (m)');
set(gca,'FontSize',12)
shg