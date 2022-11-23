% Function program [ind,k,u] = search(korig,maara,ind,raja)
% ------------------------------------------------------------------
% Search resonances by streching the geometry.
%
% ------------------------------------------------------------------
% CALLS TO : eee.m
% 12/04/00 : Seppo J‰rvemp‰‰ - RNI
% ------------------------------------------------------------------

function [ind,k,u] = search(korig,maara,ind,raja)

offset1 = -0.9; offset2 = 0; offset3 = 0.9;

laskuri = 1;
k1 = eee(offset1, maara, 0);
[k,u,siirto] = eee(offset2, maara, 0);
k3 = eee(offset3, maara, 0);

if (ind == 0)
  [ero2 , ind] = min(abs(k-korig));
  if (abs(k(ind)) < 1e-6)
    ind = ind+1;
  end
  ind;
end
ero1 = korig-k1(ind); s1 = sign(ero1);
ero2 = korig-k(ind); s2 = sign(ero2);
ero3 = korig-k3(ind); s3 = sign(ero3);
%disp([s1 s2 s3])
if ((s1 == 1) & (s2 == 1) & (s3 == 1))
  disp('?????????????????????????????')
  disp([k1(ind) k(ind) k3(ind)])
  error('Can''t compress enough')
end
if ((s1 == -1) & (s2 == -1) & (s3 == -1))
  jatka = 1;
  while ((jatka > 0) & (jatka < 4))
    disp('problem, have not stretched enough, trying to fix')
    offset1 = offset1-0.9
    
    k3 = k; offset3 = offset2; s3 = s2;
    k = k1; offset2 = offset1; s2 = s1;
    k1 = eee(offset1, maara, 0);
    ero1 = korig-k1(ind); s1 = sign(ero3);
    if (s3 == 1)
      jatka = 0;
      disp('succeeded')
    else
      jatka = jatka+1;
    end
  end
  if (jatka ~= 0)
    error('FAILED')
  end  
end
while (abs(ero2) > raja)
  if (s1 ~= s2)
    s3 = s2; offset3 = offset2;
    offset2 = 0.5*(offset1+offset3);
    s3 = s2;
  elseif (s2 ~= s3)
    s1 = s2; offset1 = offset2;
    offset2 = 0.5*(offset2+offset3);
    s1 = s2;
  end
  [k,u,siirto] = eee(offset2, maara, 0);
  ero2 = korig-k(ind);
  s2 = sign(ero2);
%  disp([k(ind), ero2, offset2, siirto])
  disp([k(ind), ero2, siirto])
end
k = k(ind); u = u(:,ind);
% -----------------------------------------------------------------------
