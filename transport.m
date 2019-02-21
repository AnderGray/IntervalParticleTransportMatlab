
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
