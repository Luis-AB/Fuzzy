function [bestFit,bestInd,BestFitIter,BestIndIter,FAccessIter] = PSOGSAFUZZY(objectiveFcn,iterations,dims,bounds)

    fis = readfis("PSOC_FUZZY.fis");
    Div_max = 0;
    
    VarMin = bounds(1,1);
    VarMax = bounds(1,2);
    nVar = dims;
    VelMax = (VarMax - VarMin)*0.2;
    VelMin = -VelMax;
    
    MaxIt = iterations;      % Maximum Number of Iterations
    nPop = 50;        % Population Size (Swarm Size)
    
    current_fitness = zeros(nPop,1);
    bestInd = zeros(1,nVar);
    bestFit = inf;
    pBestScore = ones(nPop,1) * inf;
    pBest = zeros(nPop,nVar);
    BestIndIter=zeros(MaxIt,nVar);
    BestFitIter = ones(MaxIt,1)*Inf;
    FAccessIter=zeros(MaxIt,1);
    FAccess=0;
    
    G0=1; % gravitational constant

    current_position = rand(nPop,nVar).*(VarMax-VarMin)+VarMin; %initial positions in the problem's boundary
    mass(nPop) = 0;
    velocity = .3*randn(nPop,nVar) ;

    C1=0.5; %C1 in Equation (9)
    C2=1.5; %C2 in Equation (9)


    %Parametros nuevos %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Cstuck = 0;
    maxDiv = 0;
    limitDiv = 1;
    divChangeCount = 0;
    divIter=zeros(MaxIt,1);
    c1Iter=zeros(MaxIt,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    %% Diversity initialize

    temp = zeros(1,nVar);
    for n=1:nVar
        temp2 = zeros( nPop,1);
        for h=1: nPop
            %temp2(h) = particle(h).Position(n);
            temp2(h) = current_position(h,n);
        end
        temp(n)=median(temp2);
    end
    meandim=temp;
    temp = zeros(1,nVar);
    for n=1:nVar
        temp2 = 0;
        for h=1: nPop
            %temp2 = temp2 + abs(particle(h).Position(n) - meandim(n));
            temp2 = temp2 + abs(current_position(h,n) - meandim(n));
        end
        temp(n)=temp2/ nPop;
    end
    divj=sum(temp(:));
    
    div = divj/nVar;
    
    if div > maxDiv
        maxDiv = div;
    end
    
    divPercent = div/maxDiv;
    
    %% Diversity initialize
    
    %% Historico de diversidad y valor inicial de W
    iwt = divPercent;
    divIter(1) = iwt;
    
    c1Iter(1)=C1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    %%main loop
    iter = 0 ;                  % Iterations’ counter
    while( iter < MaxIt )
        G=G0*exp(-23*iter/MaxIt); %Equation (4)
        iter = iter + 1;
        force = zeros(nPop,nVar);
        mass(nPop) = 0;
        acceleration = zeros(nPop,nVar);

        for i = 1:nPop
            %///Bound the search Space///
            %Tp=current_position(i,:)>up;Tm=current_position(i,:)<low;current_position(i,:)=(current_position(i,:).*(~(Tp+Tm)))+up.*Tp+low.*Tm;                     
    
            Flag4up = current_position(i,:)>VarMax';
            Flag4low = current_position(i,:)<VarMin';
            current_position(i,:) = (current_position(i,:).*(~(Flag4up+Flag4low)))+VarMax'.*Flag4up+VarMin'.*Flag4low; 
            %////////////////////////////
    
                                 %-------------------------------------------%
                                 %         Evaluate the population           %           
                                 %-------------------------------------------%      
            fitness = objectiveFcn(current_position(i,:));
            FAccess = FAccess + 1;
            current_fitness(i)=fitness;    
        
            if(pBestScore(i)>fitness)
                pBestScore(i)=fitness;
                pBest(i,:)=current_fitness(i,:);
            end
            if(bestFit>fitness)
                bestFit=fitness;
                bestInd=current_position(i,:);
            end
        end

        best=min(current_fitness);
        worst=max(current_fitness);

        BestFitIter(iter) = bestFit;
        FAccessIter(iter) = FAccess;
        BestIndIter(iter,:) = bestInd;
         
        for pp=1:nPop
            if current_fitness(pp)==best
                break;
            end
        end
    
        bestIndex=pp;
            
        for pp=1:nVar
            best_fit_position(iter,1)=best;
            best_fit_position(iter,pp+1)=current_position(bestIndex,pp);   
        end
                                               %-------------------%
                                               %   Calculate Mass  %
                                               %-------------------%
        for i=1:nPop
            mass(i)=(current_fitness(i)-0.99*worst)/(best-worst);    
        end

        for i=1:nPop
            mass(i)=mass(i)*5/sum(mass);
        end

                                               %-------------------%
                                               %  Force    update  %
                                               %-------------------%
        for i=1:nPop
            for j=1:nVar
                for k=1:nPop
                    if(current_position(k,j)~=current_position(i,j))
                        % Equation (3)
                        force(i,j)=force(i,j)+ rand()*G*mass(k)*mass(i)*(current_position(k,j)-current_position(i,j))/abs(current_position(k,j)-current_position(i,j));
                    end
                end
            end
        end
                                               %------------------------------------%
                                               %  Accelations $ Velocities  UPDATE  %
                                               %------------------------------------%
        for i=1:nPop
            for j=1:nVar
                if(mass(i)~=0)
                	%Equation (6)
                    acceleration(i,j)=force(i,j)/mass(i);
                end
            end
        end   

        for i=1:nPop
            for j=1:nVar
                %Equation(9)
                velocity(i,j)=iwt*velocity(i,j)+C1*rand()*acceleration(i,j) + C2*rand()*(bestInd(j)-current_position(i,j));
                velocity(i,j) = min(velocity(i,j),VelMax(1));
                velocity(i,j) = max(velocity(i,j),VelMin(1));
            end
        end
                                               %--------------------------%
                                               %   positions   UPDATE     %
                                               %--------------------------%                                               
        %Equation (10) 
        current_position = current_position + velocity ;
        current_position = min(current_position,VarMax);
        current_position = max(current_position,VarMin);

        if iter >= 2 && BestFitIter(iter) == BestFitIter(iter-1)
            Cstuck = Cstuck+1;
        else
            Cstuck = 0;
        end


        %% Diversity
    if  limitDiv ~= -1
        temp = zeros(1,nVar);
        for n=1:nVar
            temp2 = zeros(nPop,1);
            for h=1:nPop
                %temp2(h) = particle(h).Position(l);
                temp2(h) = current_position(h,n);
            
            end
            temp(n)=median(temp2);
        end
        meandim=temp;
        temp = zeros(1,nVar);
        for n=1:nVar
            temp2 = 0;
            for h=1:nPop
                %temp2 = temp2 + abs(particle(h).Position(l) - meandim(l));
                temp2 = temp2 + abs(current_position(h,n) - meandim(n));
            end
            temp(n)=temp2/nPop;
        end
        divj=sum(temp(:));
        
        div = divj/nVar;
        
        if div > maxDiv
            maxDiv = div;
        end
        
        %divPercent = div/maxDiv*100;
        divPercent = div/maxDiv;
    end
    
    %% Ajuste de W
    if divPercent <= 0.01
        VelMax=(VarMax-VarMin);
        VelMin=-VelMax;
    end
    %else
    
    %inputs = [divPercent*100 FAccess];
    inputs = [min(Cstuck,100) divPercent*100];
    output = evalfis(fis,inputs);
    
    %w = divPercent;
    iwt = output(1);
    %c1=output(2);
    %end
    divIter(iter+1) = iwt;
    c1Iter(iter+1)=C1;
    

    end
   
end
