%function gcs = calcGCS(results)
% Calculate GCS according to the yield space area

function gcs = calcGCS(results, results_WT)

%% Calculate wild type yield space area
yRange      = results_WT.yieldR;
mRange      = results_WT.muR;
pRange      = results_WT.prodR;

[muMaxWT,posMuMaxWT]    = max(mRange);

intPart     = 0;
if yRange(posMuMaxWT) == 0
    % calculate Ober/Untersumme
    for i=1:(length(mRange)-1)
        sumO        = (mRange(i)-mRange(i+1))*yRange(i);
        sumU        = (mRange(i)-mRange(i+1))*yRange(i+1);
        intPart     = ((sumO+sumU)/2)+intPart;
    end
    areaWT  = intPart;  
else
    % if GC occurs in WT
    for i=1:(posMuMaxWT-1)
        sumO        = (-mRange(i)+mRange(i+1))*yRange(i);
        sumU        = (-mRange(i)+mRange(i+1))*yRange(i+1);
        intPart     = ((sumO+sumU)/2)+intPart;
    end
    ysLow       = intPart;
    
    intPart     = 0;
    for i=posMuMaxWT:(length(mRange)-1)
        sumO        = (mRange(i)-mRange(i+1))*yRange(i);
        sumU        = (mRange(i)-mRange(i+1))*yRange(i+1);
        intPart     = ((sumO+sumU)/2)+intPart;
    end
%     areaWT  = intPart-ysLow;
    areaWT  = intPart;
end
        

% %% Calculate mutant yield space area
% yRange      = results.printData.yieldR;
% mRange      = results.printData.muR;
% pRange      = results.printData.prodR;
% 
% [muMax,posMuMax]    = max(mRange);
% 
% intPart     = 0;
% for i=1:(posMuMax-1)
%     sumO        = (mRange(i+1)-mRange(i))*yRange(i);
%     sumU        = (mRange(i+1)-mRange(i))*yRange(i+1);
%     intPart     = ((sumO+sumU)/2)+intPart;
% end
% areaMutLow      = intPart;
% 
% intPart     = 0;
% for i=posMuMax:(length(mRange)-1)
%     sumO        = (mRange(i)-mRange(i+1))*yRange(i);
%     sumU        = (mRange(i)-mRange(i+1))*yRange(i+1);
%     intPart     = ((sumO+sumU)/2)+intPart;
% end
% 
% areaMut     = intPart-areaMutLow;
% 
% gcs     = -areaMut/areaWT;



%% Alternative Calculate GCS for mutant
% Determine maximal growth rate for mutant
yRangeMut       = results.printData.yieldR;
mRangeMut       = results.printData.muR;
pRangeMut       = results.printData.prodR;

% yRangeMut       = results.yieldR;
% mRangeMut       = results.muR;
% pRangeMut       = results.prodR;

yRangeWT        = results_WT.yieldR;
mRangeWT        = results_WT.muR;
pRangeWT        = results_WT.prodR;


[muMax,posMuMax]    = max(mRangeMut);
yMuMax              = yRangeMut(posMuMax);

% area below mutant yield space
intPart     = 0;
for i=1:(posMuMax-1)
    sumO        = (mRangeMut(i+1)-mRangeMut(i))*yRangeMut(i);
    sumU        = (mRangeMut(i+1)-mRangeMut(i))*yRangeMut(i+1);
    intPart     = ((sumO+sumU)/2)+intPart;
end
areaMutLow      = intPart;

% area above mutant production envelop
% determine growth rate in wild type data
[muMax_wt,posMuMax_wt]    = max(mRangeWT);
for i=posMuMax_wt(end):length(mRangeWT)
    if mRangeWT(i)<=muMax
        posMuMax_mut_wt     = i;
        break
    end
end

% recalculate position f maximal growth (adapt to closest to wild type)
prev_muDiff     = abs(mRangeMut(posMuMax)-mRangeWT(posMuMax_mut_wt));
posMuMax_recalc     = [];
for i=posMuMax+1:length(mRangeMut)
    muDiff  = abs(mRangeMut(i)-mRangeWT(posMuMax_mut_wt));
    if muDiff > prev_muDiff
        posMuMax_recalc     = i-1;
        break
    else
        prev_muDiff     = muDiff;
    end
end
if isempty(posMuMax_recalc)
    posMuMax_recalc     = posMuMax;
end

intPart     = 0;
for i=(posMuMax_recalc):(length(mRangeMut)-1)
    sumO        = (mRangeMut(i)-mRangeMut(i+1))*yRangeMut(i);
    sumU        = (mRangeMut(i)-mRangeMut(i+1))*yRangeMut(i+1);
    intPart     = ((sumO+sumU)/2)+intPart;
end
areaMutHigh      = intPart;

intPart     = 0;
for i=posMuMax_mut_wt:(length(mRangeWT)-1)
    sumO        = (mRangeWT(i)-mRangeWT(i+1))*yRangeWT(i);
    sumU        = (mRangeWT(i)-mRangeWT(i+1))*yRangeWT(i+1);
    intPart     = ((sumO+sumU)/2)+intPart;
end
areaWtHigh      = intPart;

areaDiff_mut_wt     = areaWtHigh-areaMutHigh;


% area wild type yield space beyond maximal growth rate of mutant
intPart_low     = 0;
intPart_high    = 0;
pos             = 1;
for i=1:(length(mRange)-1)
    if mRangeWT(i) >= muMax
        if mRangeWT(i) == muMaxWT
            pos     = 0;
        end
            
        if pos
            sumO        = (-mRangeWT(i)+mRangeWT(i+1))*yRangeWT(i);
            sumU        = (-mRangeWT(i)+mRangeWT(i+1))*yRangeWT(i+1);
            intPart_low = ((sumO+sumU)/2)+intPart_low;
        else
            sumO        = (mRangeWT(i)-mRangeWT(i+1))*yRangeWT(i);
            sumU        = (mRangeWT(i)-mRangeWT(i+1))*yRangeWT(i+1);
            intPart_high = ((sumO+sumU)/2)+intPart_high;
        end
    end
end
areaWTMuMax     = intPart_high-intPart_low;


% for i=1:(length(mRange)-1)
%     if yRangeWT(i) > yMuMax
%         % if yield at maximal growth of the mutant is not on wild type
%         % yield space hull curve
%         intPart     = intPart+(mRangeWT(i-1)-muMax)*((yRangeWT(i-1)+yMuMax)/2)
%         break;
%     else
%         sumO        = (mRangeWT(i)-mRangeWT(i+1))*yRangeWT(i);
%         sumU        = (mRangeWT(i)-mRangeWT(i+1))*yRangeWT(i+1);
%         intPart     = ((sumO+sumU)/2)+intPart;
%     end
% end
% areaWTMuMax  = intPart;

%
if 0
% if yRange(posMuMax) ~= 0
    intPart     = 0;
    for i=1:posMuMaxWT
        sumO        = (mRangeWT(i+1)-mRangeWT(i))*yRangeWT(i);
        sumU        = (mRangeWT(i+1)-mRangeWT(i))*yRangeWT(i+1);
        intPart     = ((sumO+sumU)/2)+intPart;
    end
    areaWTLow   = intPart;
else
    areaWTLow   = 0;
end

% % consider WT area beyond mu max of the mutant
% areaLow     = areaWTMuMax+areaMutLow-areaWTLow;

% DO NOT consider WT area beyond mu max of the mutant
areaLow     = areaMutLow-areaWTLow;
% areaWT
% areaDiff_mut_wt
% areaWTMuMax

[yMin,~]    = min(yRangeMut);   % minimal accessible yield
muMaxYzero  = sum(mRangeMut(yRangeMut==0));


% Account for guaranteed yield at maximal growth
[minGY,~]   = min(yRangeMut(mRangeMut==muMax));
[maxY,~]    = max(yRange);  


% Calculate GCS (case analysis)
if abs(areaWT-areaWTMuMax) < 1e-4
    gcs     = -2;
else
    if yMin > 0
    %     gcs         = ((areaLow-areaDiff_mut_wt)/(areaWT-areaWTMuMax))*(minGY/maxY);
        gcs         = ((areaLow)/(areaWT-areaWTMuMax))*(minGY/maxY);
    elseif yMin == 0 && muMaxYzero == 0
    %     gcs         = -1+((areaLow/(areaWT+areaDiff_mut_wt-areaWTMuMax))*(minGY/maxY));
        gcs         = -1+(((areaLow)/(areaWT-areaWTMuMax))*(minGY/maxY));
    elseif muMaxYzero > 0
    %     gcs         = -2+((areaLow/(areaWT+areaDiff_mut_wt-areaWTMuMax))*(minGY/maxY)); 
        gcs         = -2+(((areaLow)/(areaWT-areaWTMuMax))*(minGY/maxY));
    end
end

end
