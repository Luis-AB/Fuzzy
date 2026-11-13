
SpeciesCount =  (popsize-1):-1:0; % SpeciesCount is also known as ranking

switch MigrateModel
                
        % Compute migration rates, assuming the population is sorted from most fit to least fit
        case 1  
            lambda = I*0.5*ones(1,popsize); % immigration rate
            parmu = E*SpeciesCount/popsize; % emigration rate
            
        case 2 
            lambda = I*(1-SpeciesCount/popsize); % immigration rate
            parmu = E*0.5*ones(1,popsize);  % emigration rate

        case 3  
            lambda = I*(1-SpeciesCount/popsize); % immigration rate
            parmu = E*SpeciesCount/popsize; % emigration rate

        case 4 
            ip = ceil((popsize+1)*0.5);
            lambda(ip+1:popsize) = I;
            lambda(1:ip) = I*2*(1-SpeciesCount(1:ip)/popsize); % immigration rate
            parmu(1:ip) = E;
            parmu(ip+1:popsize) = E*2*SpeciesCount(ip+1:popsize)/popsize; % emigration rate 
            parmu(parmu>E) = E; lambda(lambda>I) = I;
            
         case 5 
            lambda = I*(1-SpeciesCount/popsize).^2;  % immigration rate
            parmu = E*(SpeciesCount/popsize).^2;  % emigration rate

        case 6
            lambda = I*0.5*(cos(SpeciesCount*pi/popsize)+1); % immigration rate
            parmu = E*0.5*(-cos(SpeciesCount*pi/popsize)+1); % emigration rate
            
end   
 