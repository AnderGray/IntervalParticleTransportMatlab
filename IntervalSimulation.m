%%%
%   A prototype for propagating bounded hulls though a particle transport
%   monte carlo
%%%

withXSUncertainty = 1;
totalMC = 1;
numFe = 613;

if withXSUncertainty
    
    XSDIR = '/Users/akgray/Documents/neutronics/juliaNeutronics/tendlHDF5/Fe056';

    dir = strcat(XSDIR,'/Fe056_0000.h5');
    energy = h5read(dir,"/energy");

    scatFe = zeros(numFe, numel(energy));
    abspFe = zeros(numFe, numel(energy));
    totalFe = zeros(numFe, numel(energy));


    for i =1:numFe
        if i <10
            nums = strcat("000",string(i));
        elseif i<100
            nums = strcat("00",string(i));
        else
            nums = strcat("0",string(i));
        end
        files = strcat("/Fe056_",nums,".h5");
        scatFe(i,:) = h5read(strcat(XSDIR,files),"/elastic");
        abspFe(i,:) = h5read(strcat(XSDIR,files),"/absorption");
        totalFe(i,:) = h5read(strcat(XSDIR,files),"/total");
    end

    totalMaxFe = zeros(numel(energy),1);
    totalMinFe = zeros(numel(energy),1);
    abspMaxFe  = zeros(numel(energy),1);
    abspMinFe  = zeros(numel(energy),1);
    scatMaxFe  = zeros(numel(energy),1);
    scatMinFe  = zeros(numel(energy),1);

    %abspFe = abspFe .*0.1;

    for i = 1:numel(energy)
        totalMaxFe(i) = max(totalFe(:,i));
        totalMinFe(i) = min(totalFe(:,i));

        abspMaxFe(i)  = max(abspFe(:,i));
        abspMinFe(i)  = min(abspFe(:,i));

        scatMaxFe(i)  = max(scatFe(:,i));
        scatMinFe(i)  = min(scatFe(:,i));
    end

    %scatFe  = 0;
    %abspFe  = 0;
    %totalFe = 0;
    
    %%%% CONTINUE HERE, WAS GOING TO CREATE INTERVAL XS
    
    
    total = [totalMinFe, totalMaxFe];
    scat  = [scatMinFe, scatMaxFe];
    absp  = [abspMinFe, abspMaxFe];
    
else

    filename = '/Users/akgray/Documents/neutronics/juliaNeutronics/tendlHDF5/Fe056/Fe056_0000.h5';

    energy  = h5read(filename,"/energy");
    total   = h5read(filename, "/total");
    scat    = h5read(filename, "/elastic");
    absp    = h5read(filename, "/absorption");
    
    
    total   = [total, total];
    scat    = [scat, scat];
    absp    = [absp, absp];
    

end

disp("Nuclear Data read")

%% Interval Monte Carlo
numBatch        = 10;
numParticles    = 100000;
energyTally     = 1e4:5e4:2e7;
tallymin        = zeros(numel(energyTally),numBatch);
tallymax        = zeros(numel(energyTally),numBatch);


tic;
parfor  j =1:numBatch
    localTallyMax = zeros(numel(energyTally),1);
    localTallyMin = zeros(numel(energyTally),1);
    for i = 1:numParticles
        current = [0,0];
        particleEnergy = [14.2e6,14.2e6];
        positionUncertainty = 0;
        current(1,2) = positionUncertainty + current(1,1);  
        Alive = 1;
        while Alive == 1
            randNums = rand(3,1);
            
            [particleEnergyOut, currentOut, Alive] = transport(particleEnergy, current, total, scat, absp, energy, randNums);

            [dis, Idx] = min( abs( energyTally' - particleEnergy ));
            dmin = currentOut(1) - current(1);
            dmax = currentOut(2) - current(2);

            Indexes = Idx(1):Idx(2);

            localTallyMin(Indexes) = localTallyMin(Indexes) + dmin;
            localTallyMax(Indexes) = localTallyMax(Indexes) + dmax;

            particleEnergy = particleEnergyOut;
            current = currentOut;
        end
    end
    tallymin(:,j) = localTallyMin;
    tallymax(:,j) = localTallyMax;
end
toc

tallyMean = zeros(numel(energyTally),2);
tallyVar = zeros(numel(energyTally),2);

for k = 1:numel(energyTally)
    
    [tallyMean(k,:), tallyVar(k,:)] = IntervalStatistics([tallymin(k,:);tallymax(k,:)]');

end 

figure;
plot(energyTally,tallyMean(:,1));
hold on
plot(energyTally,tallyMean(:,2));


figure;
loglog(energyTally(10:end),tallyMean(10:end,1));
hold on
loglog(energyTally(10:end),tallyMean(10:end,2));

figure;
plot(energyTally(10:end),tallyVar(10:end,1));
hold on
plot(energyTally(10:end),tallyVar(10:end,2));

disp("Interval MC finished")

%% Total Monte Carlo

if totalMC
    tic;
    numSimulations = 15;
    numBatch = 10;
    tendl                       = unidrnd( numFe, [numSimulations, 1]);
    current = [0,0];
    particleEnergy = [14.2e6,14.2e6];
    positionUncertainty(1) = 0;
    current(1,2) = positionUncertainty(1) + current(1,1);
    tmcTallyMean = zeros(numel(energyTally),numSimulations);
    tmcTallyVar = zeros(numel(energyTally),numSimulations);
    
    for k = 1:numSimulations
        disp("Starting simulation:")
        k
        tmcTallySingle = zeros(numBatch,numel(energyTally));
        parfor j = 1: numBatch
            
             totes = repmat(totalFe(tendl(k), :)',[1,2]); 
             scats = repmat(scatFe(tendl(k), :)',[1,2]);
             abps  = repmat(abspFe(tendl(k), :)',[1,2]);
            
            
            disp("Starting batch:")
            j
            tmcTallyBatch = zeros(numel(energyTally),1);
            for i = 1:numParticles
                Alive = 1; 
                pointParticleEnergy         = particleEnergy; %repmat((particleEnergy(1,1) + (particleEnergy(1,2) - particleEnergy(1,1)) .* rand(1,1))', [1,2]);
                pointParticlePosition       = current; %repmat((current(1,1)        + (current(1,2)        - current(1,1))        .* rand(1,1))', [1,2]);
                while Alive == 1
                    randNums = rand(3,1);
                    if withXSUncertainty
                        [pointParticleEnergyOut, pointParticlePositionOut, Alive] = transport( pointParticleEnergy, pointParticlePosition, ...
                            totes, scats,...
                            abps, energy, randNums); 
                    else
                        [pointParticleEnergyOut, pointParticlePositionOut, Alive] = transport( pointParticleEnergy,...
                            pointParticlePosition, total, scat,absp, energy, randNums); 
                    end
                    [~, Idx] = min( abs( energyTally' - pointParticleEnergy(1) ));
                    d = pointParticlePositionOut(1) - pointParticlePosition(1);

                    tmcTallyBatch(Idx) = tmcTallyBatch(Idx) + d;
                    pointParticleEnergy = pointParticleEnergyOut;
                    pointParticlePosition = pointParticlePositionOut;
                end
            end 
            tmcTallySingle(j,:) = tmcTallyBatch(:);
        end
        for l = 1:numel(energyTally)
            tmcTallyMean(l,k) = mean(tmcTallySingle(:,l));
            tmcTallyVar(l,k) = var(tmcTallySingle(:,l));
        end
    end
    toc
end

save('AnderGrandPlan')

%% Plotting

energyLocals = [0.5,1,5,10,14.2].*1e6;
[~, Idx] = min( abs( energyTally' - energyLocals ));

figure;
loglog(energyTally(10:end),tallyMean(10:end,1),'LineWidth',2);
hold on
loglog(energyTally(10:end),tallyMean(10:end,2),'LineWidth',2);

for k = 1:numSimulations
    loglog(energyTally(10:end),tmcTallyMean(10:end,k),'r-');
end

y = ylim;

for i =1:numel(energyLocals)
   plot([energyLocals(i),energyLocals(i)],[y(1),y(2)]) 
end

for i = 1:numel(energyLocals)
    
    figure;
    
    range = tallyMean(Idx(i),1)-4*sqrt(tallyVar(Idx(i),2)):0.01:tallyMean(Idx(i),2)+4*sqrt(tallyVar(Idx(i),2));

    cdfupper = zeros(numel(range),1);
    cdflower = zeros(numel(range),1);

    for j=1:numel(range)
        upperVal1 = normcdf(range(j),tallyMean(Idx(i),1),sqrt(tallyVar(Idx(i),1)));
        upperVal2 = normcdf(range(j),tallyMean(Idx(i),1),sqrt(tallyVar(Idx(i),2)));
        cdfupper(j) = max(upperVal1,upperVal2);

        lowerVal1 = normcdf(range(j),tallyMean(Idx(i),2),sqrt(tallyVar(Idx(i),1)));
        lowerVal2 = normcdf(range(j),tallyMean(Idx(i),2),sqrt(tallyVar(Idx(i),2)));
        cdflower(j) = min(lowerVal1,lowerVal2);
    end
    plot(range,cdflower,'LineWidth',2);
    hold on;
    plot(range,cdfupper,'LineWidth',2);
    
    for k = 1:numSimulations
        plot(range, normcdf(range,tmcTallyMean(Idx(i),k),sqrt(tmcTallyVar(Idx(i),k))), 'r-' );
    end
    ecdf(tmcTallyMean(Idx(i),:))
    
    oo=['Distributions at: ',num2str(energyLocals(i))];
    clear title;
    title(oo);
    
end





function [particleEnergyOut, particlePositionOut, alive] = transport(particleEnergyIn, particlePositionIn, total, scat, absp, energy, randSeq)

    totalIntersTemp = sort( interp1( energy, total, particleEnergyIn, 'linear'));
    scatIntersTemp  = sort( interp1( energy, scat,  particleEnergyIn, 'linear'));
    abspIntersTemp  = sort( interp1( energy, absp,  particleEnergyIn, 'linear'));

    totalInters = [min(totalIntersTemp(:,1)), max(totalIntersTemp(:,2))];
    scatInters = [min(scatIntersTemp(:,1)), max(scatIntersTemp(:,2))];
    abspInters = [min(abspIntersTemp(:,1)), max(abspIntersTemp(:,2))];

    [dis, Idx] = min( abs( energy - particleEnergyIn ));

    if energy(Idx(1)) < particleEnergyIn(1)
       Idx(1)=Idx(1)+1;
    end
    if energy(Idx(2)) > particleEnergyIn(2)
        Idx(end)=Idx(2)-1;
    end

    totalxs = [min( [ total([Idx(1):Idx(2)], 1)', totalInters ] ), max( [ total([Idx(1):Idx(2)], 2)', totalInters ] )];
    scatxs  = [min( [ scat([Idx(1):Idx(2)], 1)',  scatInters ] ),  max( [ scat( [Idx(1):Idx(2)], 2)', scatInters ] )];
    abspxs  = [min( [ absp([Idx(1):Idx(2)], 1)',  abspInters] ),   max( [ absp( [Idx(1):Idx(2)], 2)', abspInters ] )];

    jump = - log(randSeq(1))./totalxs;
    particlePositionOut = particlePositionIn + sort(jump);

    totalSummed = sort(scatxs) + sort(abspxs);
    probAbsp = [abspxs(1)/totalSummed(2), abspxs(2)/totalSummed(1)];

    %positionUncertainty = particlePositionOut(2) - particlePositionOut(1);

    if randSeq(2) < probAbsp(1)
        %j
        alive = 0;
        particleEnergyOut = [0,0];
    else
        awr = 55.454479886265396;
        e_loss = ((awr-1)/(awr+1))^2;
        particleEnergyOut = particleEnergyIn *randSeq(3)*e_loss;
        alive = 1;
    end
    if particleEnergyOut(1) <= 10000
        alive = 0;
    end

end

