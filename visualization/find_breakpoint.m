function breakPoint = find_breakpoint(pathY)
    stayInPlace = pathY(1:end-1)-pathY(2:end);
    stayInPlace = stayInPlace==0;
    aboveBase = pathY>0.6;
    belowBase = pathY<0.3;
    outOfBase = aboveBase+belowBase;
    
%     firstHalfTop = pathY>=0.9;
%     firstHalfBot = pathY<=0.1;
%     
%     breakTop = find(firstHalfTop, 1, 'first');
%     breakBot = find(firstHalfBot, 1, 'first');
%     if isempty(breakTop)
%         breakTop = 1;
%     end
%     if isempty(breakBot)
%         breakBot = 1;
%     end
    breakPoint = find(stayInPlace.*outOfBase(1:end-1), 1, 'first');
    if isempty(breakPoint)
        breakPoint = length(pathY);
    end
end