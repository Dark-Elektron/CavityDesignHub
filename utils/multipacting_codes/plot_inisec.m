% Function program plot_inisec(s)
% ---------------------------------------------------------------------
% Plots the initial points, the geometry and the secondary yield curve.
%
% ---------------------------------------------------------------------
% CALLS TO : check_fieldparam.m, check_inputs.m, check_inputs_win.m, 
%            error_message.m, clear_window.m
% 14/05/00 : Pasi Yla-Oijala, Rolf Nevanlinna Institute
% ---------------------------------------------------------------------

function plot_inisec(s)

cl = 0;
ok = check_fieldparam;
if ok > 0
  load fieldparam
  gtype = fieldparam(1);
else
  return;
end  

if s == 1
  if gtype <= 2
    ok = check_inputs;
    if ok > 0
      load counter_initials  
    end
  else
    ok = check_inputs_win;
    if ok > 0
      load counter_initialsl
      load counter_initialsr
    end
  end  
else
  if gtype <= 2
    ok = exist('initials_temp.mat');
    if ok == 0
      cl = error_message(['Initials file is missing. Choose Create Inputs ' ...
			  'in menu Inputs.']);
      return;
    else
      load initials_temp
    end    
  else 
    ok1 = exist('initials_templ.mat');
    ok2 = exist('initials_tempr.mat');
    ok  = ok1*ok2;
    if ok0 == 0
      cl = error_message(['Initials file is missing. Choose Create Inputs ' ...
			  'in menu Inputs.']);
      return;
    else
      load initials_templ
      load initials_tempr
    end
  end
end  

if ok > 0
  cl = error_message('Plotting the geometry and the secondary yield.');

  Window = figure(3);
  clf;
  set(Window,'name','MULTIPAC - Input Window II');

  subplot(2,1,1)
  
  if gtype <= 2
    load geodata.n

    n  = length(geodata(:,1));
    gr = geodata(4:n,1);
    gz = geodata(4:n,2);
    m  = length(initials(:,1));
    ir = initials(:,1);
    iz = initials(:,2);

    plot(gz,gr,'-r',iz,ir,'bo')
  else
    load geodatal.n
    nl  = length(geodatal(:,1));
    grl = geodatal(4:nl,1);
    gzl = geodatal(4:nl,2);

    load geodataw.n
    nw  = length(geodataw(:,1));
    grw = geodataw(4:nw,1);
    gzw = geodataw(4:nw,2);

    load geodatar.n
    nr  = length(geodatar(:,1));
    grr = geodatar(4:nr,1);
    gzr = geodatar(4:nr,2);

    ml  = length(initialsl(:,1));
    irl = initialsl(:,1);
    izl = initialsl(:,2);
    mr  = length(initialsr(:,1));
    irr = initialsr(:,1);
    izr = initialsr(:,2);

    plot(gzl,grl,'-r',gzw,grw,'-r',gzr,grr,'-r',izl,irl,'bo',izr,irr,'bo')
    hold on
    fill(gzw,grw,'y')
    hold off
    m  = mr+ml;
  end

  axis equal
  grid
  xlabel('z axis [m]')
  ylabel('r axis [m]')
  title(['Initial Points         number of points ' num2str(m) ])
  
  fileis1 = exist('secy1');
  fileis2 = exist('secy2');
  fileis3 = exist('secy3');
  fileis4 = exist('secy4');
  
  subplot(2,1,2)
  if fileis1 > 0
    load secy1
    n        = length(secy1(:,1));
    energy   = secy1(1:n-1,1);
    secyield = secy1(1:n-1,2);
    plot(energy,secyield,'b',[0,max(energy)],[1,1],'-r')
    hold on
  else 
    error_message('Secondary yield curve secy1 missing.');
    return;
  end  
  if fileis2 > 0
    load secy2
    n        = length(secy2(:,1));
    energy   = secy2(1:n-1,1);
    secyield = secy2(1:n-1,2);
    plot(energy,secyield,'b')
  end  
  if fileis3 > 0
    load secy3
    n        = length(secy3(:,1));
    energy   = secy3(1:n-1,1);
    secyield = secy3(1:n-1,2);
    plot(energy,secyield,'b')
  end  
  if fileis4 > 0
    load secy4
    n        = length(secy4(:,1));
    energy   = secy4(1:n-1,1);
    secyield = secy4(1:n-1,2);
    plot(energy,secyield,'b')
  end  
  hold off
  grid
  xlabel('Impact energy [eV]')
  ylabel('Secondary yield')
  title('Secondary yield')
end  

if cl == 1
  clear_window;
end  
% ---------------------------------------------------------------------

