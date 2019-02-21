%%%
%   A prototype for propagating bounded hulls though a particle transport
%   monte carlo
%%%

withXSUncertainty = 1;
withParticles = 1;
numFe = 613;

if withXSUncertainty
    

    XSDIR = '/home/angray/Documents/neutronics/juliaNeutronics/tendlHDF5/Fe056';

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

    filename = '/home/angray/Documents/neutronics/juliaNeutronics/tendlHDF5/Fe056/Fe056_0000.h5';

    energy = h5read(filename,"/energy");
    total = h5read(filename, "/total");
    scat = h5read(filename, "/elastic");
    absp = h5read(filename, "/absorption");
    
    
    total = [total, total];
    scat  = [scat, scat];
    absp  = [absp, absp];
    

end

%%

steps = 10;
particleEnergy = zeros(steps,2);
current = zeros(steps,2);

particleEnergy(1,:) = [14.2e6,14.2e6];

positionUncertainty(1) = 0.2;

current(1,2) = positionUncertainty(1) + current(1,1);  


randNums = rand(steps, 2);

tic;
[particleEnergy, current, positionUncertainty] = transport(particleEnergy, current, total, scat, absp, energy, randNums);
toc

% 
% figure(1);
% for k = 1:steps
%     rectangle('Position',[current(k,1),particleEnergy(k,1),current(k,2)-current(k,1),particleEnergy(k,2)-particleEnergy(k,1)])
%     hold all
%     plot((current(k,1)+current(k,2))/2,(particleEnergy(k,1)+particleEnergy(k,2))/2,'*')
%     drawnow
% end
% 
% title("Neutron phase space");
% ylabel("energy (eV)");
% xlabel("position (arb)");

positionUncertainty


if withParticles
    
    numParticles = 1000;
    tendl        = unidrnd( numFe, [numParticles, 1]);
    pointParticleEnergy         = zeros(steps, numParticles *2);
    pointParticlePosition       = zeros(steps, numParticles *2);
    pointParticleEnergy(1,:)    = repmat((particleEnergy(1,1) + (particleEnergy(1,2) - particleEnergy(1,1)) .* rand(numParticles,1))', [1,2]);
    pointParticlePosition(1,:)  = repmat((current(1,1)        + (current(1,2)        - current(1,1))        .* rand(numParticles,1))', [1,2]);
    
    tic;
    for i = 1:numParticles
        
        if withXSUncertainty
            [pointParticleEnergy(:,[i,i+numParticles]), pointParticlePosition(:,[i,i+numParticles]), ~] = transport( pointParticleEnergy(:,[i,i+numParticles])...
                , pointParticlePosition(:,[i,i+numParticles]), ...
                repmat(totalFe(tendl(i), :)',[1,2]), repmat(scatFe(tendl(i), :)',[1,2]), repmat(abspFe(tendl(i), :)',[1,2]), energy, randNums); 
        else
            [pointParticleEnergy(:,[i,i+numParticles]), pointParticlePosition(:,[i,i+numParticles]), ~] = transport( pointParticleEnergy(:,[i,i+numParticles])...
                , pointParticlePosition(:,[i,i+numParticles]), ...
                total, scat,absp, energy, randNums); 
        end
       % plot(pointParticlePosition(:,i),pointParticleEnergy(:,i) , '*')
       % drawnow
        
    end 
    toc
end


figure;
for k = 1:steps
    rectangle('Position',[current(k,1),particleEnergy(k,1),current(k,2)-current(k,1),particleEnergy(k,2)-particleEnergy(k,1)])
    hold on
    plot((current(k,1)+current(k,2))/2,(particleEnergy(k,1)+particleEnergy(k,2))/2,'*')
end

title("Neutron phase space");
ylabel("energy (eV)");
xlabel("position (arb)");


if withParticles

    for i = 1:numParticles
        plot(pointParticlePosition(:,i),pointParticleEnergy(:,i) , '*')
    end
end



function [particleEnergy, particlePosition, positionUncertainty] = transport(particleEnergy, particlePosition, total, scat, absp, energy, randSeq)


steps = length(randSeq);
positionUncertainty = zeros(steps,1);

    for i = 2:steps
        j = i-1;
        
        
        totalIntersTemp = sort( interp1( energy, total, particleEnergy(j,:), 'linear'));
        scatIntersTemp  = sort( interp1( energy, scat,  particleEnergy(j,:), 'linear'));
        abspIntersTemp  = sort( interp1( energy, absp,  particleEnergy(j,:), 'linear'));
        
        totalInters = [min(totalIntersTemp(:,1)), max(totalIntersTemp(:,2))];
        scatInters = [min(scatIntersTemp(:,1)), max(scatIntersTemp(:,2))];
        abspInters = [min(abspIntersTemp(:,1)), max(abspIntersTemp(:,2))];

        [dis, Idx] = min( abs( energy - particleEnergy(j,:) ));

        if energy(Idx(1)) < particleEnergy(j,1)
           Idx(1)=Idx(1)+1;
        end
        if energy(Idx(2)) > particleEnergy(j,2)
            Idx(end)=Idx(2)-1;
        end

        totalxs = [min( [ total([Idx(1):Idx(2)], 1)', totalInters ] ), max( [ total( [Idx(1):Idx(2)], 2)', totalInters ] )];
        scatxs  = [min( [ scat([Idx(1):Idx(2)], 1)', scatInters ] ),  max( [ scat( [Idx(1):Idx(2)], 2)', scatInters ] )];
        abspxs  = [min( [ absp([Idx(1):Idx(2)], 1)', abspInters] ),  max( [ absp( [Idx(1):Idx(2)], 2)', abspInters ] )];

        jump = - log(randSeq(j,1))./totalxs;
        particlePosition(i,:) = particlePosition(j,:) + sort(jump);

        totalSummed = sort(scatxs) + sort(abspxs);
        %probAbsp = [abspxs(1)/totalSummed(2), abspxs(2)/totalSummed(1)];
        probAbsp = [abspxs(1)/totalSummed(2), abspxs(2)/totalSummed(1)];
        
        positionUncertainty(i) = particlePosition(i,2) - particlePosition(i,1);

        if randSeq(j,2) < probAbsp(1)
            %j
            break;
        else
            particleEnergy(i,:) = particleEnergy(j,:) .*0.9;
        end
    end
end


