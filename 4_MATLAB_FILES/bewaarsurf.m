function bewaarsurf(handvat,naam,fontsize,res)

% bewaarsurf(handvat,naam,fontsize)
% handvat = fig nummer, naam = naam file, fontsize = fontsize,
% res=resolution
% OBSOLETE, USE STORESURF() INSTEAD

if(nargin<4)
    res=150;
end

leg = findobj(handvat,'Tag','legend');
set(leg,'FontSize',fontsize)

C = findobj(handvat,'Tag','Colorbar');
set(C,'FontSize',fontsize)

lines = findobj(handvat,'Type','Line');

for i = 1:length(lines)
    if (get(lines(i),'linewidth') <2);
        set(lines(i),'linewidth',2);
    end
end

as = findobj(handvat,'Type','axes');

set(as,'FontSize',fontsize)

for i = 1:length(as)
    % if (get(as(i),'Tag') == '')
    X = get(as(i),'Xlabel');
    set(X,'FontSize',fontsize);
%    set(X,'FontSize',fontsize,'FontWeight','Bold');
    Y = get(as(i),'Ylabel');
    set(Y,'FontSize',fontsize);
%    set(Y,'FontSize',fontsize,'FontWeight','Bold');
    T = get(as(i),'Title');
    set(T,'FontSize',fontsize);
%    set(T,'FontSize',fontsize,'FontWeight','Bold');
    % end
end

% saveas(handvat,naam)
print(['-f' int2str(handvat)],'-djpeg100','-zbuffer',['-r' int2str(res)],naam)
% eval(['!jpeg2ps -q ' naam '.jpg >' naam '.eps'])
% print(['-f' int2str(handvat)],'-depsc','-zbuffer',['-r' int2str(res)],naam)
