% Function program plot_secondary_yield.m
% ---------------------------------------------------------------------
% Plots the used secondary yield curve.
%
% ---------------------------------------------------------------------
% CALLS TO : error_message.m, clear_window.m
% 06/05/00 : Pasi Yla-Oijala, Rolf Nevanlinna Institute
% ---------------------------------------------------------------------

function plot_secondary_yield

fileis1 = exist('secy1');
fileis2 = exist('secy2');
fileis3 = exist('secy3');
fileis4 = exist('secy4');

cl = 0;
fileis0 = exist('param');
if fileis0 == 0
  cl = error_message('Parameter file param missing.');
else  
  if fileis1 == 0
    cl = error_message('Secondary yield file secy1 missing.');
  else
    cl = error_message('Plotting the secondary yield curve.');

    load secy1

    n = length(secy1(:,1));
    energy   = secy1(1:n-1,1);
    secyield = secy1(1:n-1,2);

    load param
    emin = param(8);
    emax = param(9);
    
    Window = figure(3);
    clf;
    set(Window,'name','MULTIPAC - Input Window II');

    subplot(1,1,1)
    plot(energy,secyield,'b',[0,max(energy)],[1,1],'-r')
    ax = axis;
    hold on
    if emax > ax(2)
      plot([emin,emin],[0,ax(4)],'-k')
    else
      plot([emin,emin],[0,ax(4)],'-k',[emax,emax],[0,ax(4)],'-k')
    end      
    grid
    xlabel('Impact energy [eV]')
    ylabel('Secondary yield')
    title(['MultiPac 2.0                   Secondary yield           ' date ])
  
    if fileis2 > 0
      load secy2
      n = length(secy2(:,1));
      energy   = secy2(1:n-1,1);
      secyield = secy2(1:n-1,2);
      plot(energy,secyield,'-k')
    end

    if fileis3 > 0
      load secy3
      n = length(secy3(:,1));
      energy   = secy3(1:n-1,1);
      secyield = secy3(1:n-1,2);
      plot(energy,secyield,'-g')
    end

    if fileis4 > 0
      load secy4
      n = length(secy4(:,1));
      energy   = secy4(1:n-1,1);
      secyield = secy4(1:n-1,2);
      plot(energy,secyield,'-m')
    end
    hold off
  end
end  

if cl == 1
  clear_window;
end  
% ---------------------------------------------------------------------
