function spy(sys,name)
    subplot(1,3,1);spy(sys.M); title('spy(M)');
    subplot(1,3,2);spy(sys.D); title('spy(D)');
    subplot(1,3,3);spy(sys.K); title('spy(K)');
    
    if nargin > 1
        onetitle(name);
    elseif ~isempty(sys.Name)
        onetitle(sys.Name);
    end
end

function onetitle(str)
    %   Create one common title for different subplots
    set(gcf,'NextPlot','add');
    ha = axes; h = title(str,'Interpreter','none');
    set(ha,'Visible','off');
    set(h,'Visible','on'); 
end