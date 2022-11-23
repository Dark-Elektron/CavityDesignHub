% Function program plot_initials(s)
% ---------------------------------------------------------------------
% Plots the initial points and the geometry.
%
% ---------------------------------------------------------------------
% CALLS TO : check_fieldparam.m, error_message.m
% 06/03/00 : Pasi Yla-Oijala, Rolf Nevanlinna Institute
% 14/03/00 :                 - check for inputs added
% ---------------------------------------------------------------------

function plot_initials(s)

cl = 0;
fileis = check_fieldparam;
if fileis > 0
  load fieldparam
  gtype = fieldparam(1);  
  
  ok = 1;
  if s == 1
    if gtype <= 2
      fileis = exist('counter_initials.mat');
    else
      fileis1 = exist('counter_initialsl.mat');
      fileis2 = exist('counter_initialsr.mat');
      fileis  = fileis1*fileis2;
    end
    if fileis == 0
      cl = error_message(['Initials file is missing. Choose Create Inputs '... 
			  'in menu Inputs.']);
      ok = 0;
    else
      if gtype <= 2
        load counter_initials 
      else 
        load counter_initialsl
        load counter_initialsr
        initials = [initialsl;initialsr];
      end
    end
  else
    fileis = exist('initials_temp.mat');
    if fileis == 0
      cl = error_message(['Initials file is missing. Choose Create Inputs '...
			  'in menu Inputs.']);
    ok = 0;
    else
      load initials_temp
    end
  end  

  if ok == 1
    if s > 0
      cl = error_message('Plotting the initials.');
    end
    
    Window = figure(3);
    clf;
    set(Window,'name','MULTIPAC - Input Window II');
    subplot(1,1,1)

    if gtype <= 2
      load geodata.n
      n  = length(geodata(:,1));
      gr = geodata(4:n,1);
      gz = geodata(4:n,2);
      ir = initials(:,1);
      iz = initials(:,2);
      m  = length(ir);
    
      plot(gz,gr,'-r',iz,ir,'bo')
    else
      load geodatal.n
      load geodataw.n
      load geodatar.n

      nl  = length(geodatal(:,1));
      grl = geodatal(4:nl,1);
      gzl = geodatal(4:nl,2);

      nw  = length(geodataw(:,1));
      grw = geodataw(4:nw,1);
      gzw = geodataw(4:nw,2);

      nr  = length(geodatar(:,1));
      grr = geodatar(4:nr,1);
      gzr = geodatar(4:nr,2);

      m  = length(initials(:,1));
      ir = initials(:,1);
      iz = initials(:,2);

      plot(gzl,grl,'-r',gzw,grw,'-r',gzr,grr,'-r',iz,ir,'bo')
      hold on
      fill(gzw,grw,'y')
      hold off
    end
    axis equal
    grid
    xlabel('z axis [m]')
    ylabel('r axis [m]')
    title(['MultiPac 2.0        Initial Points         number of points ' ...
	    num2str(m) '        ' date ])
  end  
end
if cl == 1
  clear_window;
end  
% ---------------------------------------------------------------------
