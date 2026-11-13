function [bestFit,bestInd,BestFitIter,BestIndIter,FAccessIter] = BLPSOFUZZY(objectiveFcn,iterations,dims,bounds)

    fis = readfis("PSOC_FUZZY.fis");
    Div_max = 0;
    
    D = dims;            % Number of Decision Variables
    Xmin = bounds(1,1);             % Lower Bound of Decision Variables
    Xmax = bounds(1,2);            % Upper Bound of Decision Variables
    maxGEN = iterations;      % Maximum Number of Iterations
    popsize = 50;        % Population Size (Swarm Size)

    
    %New Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Cstuck = 0;
    maxDiv = 0;
    limitDiv = 1;
    divChangeCount = 0;
    divIter=zeros(maxGEN,1);
    c1Iter=zeros(maxGEN,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    BestIndIter = zeros(maxGEN,D);
    BestFitIter = ones(maxGEN,1)*Inf;
    FAccessIter = zeros(maxGEN,1);
    FAccess = 0;

    pBest = zeros(popsize,D);
    val_pBest = ones(popsize,1) * inf;
    pBest_ind = zeros(popsize,D);
    pBest_f = zeros(popsize,D);
    
    %iwt = 0.9 - (1 : maxGEN) * (0.7 / maxGEN);
    c = 1.49445;
    I = 1; % max immigration rate 
    E = 1; % max emigration rate 

    % Compute migration rates, assuming the population is sorted from most fit to least fit
    MigrateModel = 5;
    migration_models;
  
    % Initialize the main population
	X = repmat(Xmin, popsize, 1) + rand(popsize,D) .* (repmat(Xmax-Xmin, popsize, 1));
    
    for i=1:popsize
        val_X = objectiveFcn(X(i,:));
        FAccess = FAccess + 1;
        pBest(i,:) = X(i,:); 
        val_pBest(i) = val_X;
    end
    
   	[~,indexG] = min(val_pBest);
   	bestInd = pBest(indexG,:); bestFit = val_pBest(indexG,:);    
   	Vmax = (Xmax - Xmin)*0.2;  Vmin = -Vmax;    
   	V = repmat(Vmin,popsize,1) + rand(popsize,D).*repmat(Vmax-Vmin,popsize,1);
    
    gapm = 5; stay_num = zeros(1,popsize);  
    cc = 0;
    
  	for i=1:popsize
        pBest_ind(i,:) = LearnIndex_BLPSO(val_pBest,popsize,D,i,parmu,lambda);
    end 

    %% Diversity initialize

    temp = zeros(1,D);
    for n=1:D
        temp2 = zeros( popsize,1);
        for h=1: popsize
            %temp2(h) = particle(h).Position(n);
            temp2(h) = X(h,n);
        end
        temp(n)=median(temp2);
    end
    meandim=temp;
    temp = zeros(1,D);
    for n=1:D
        temp2 = 0;
        for h=1: popsize
            %temp2 = temp2 + abs(particle(h).Position(n) - meandim(n));
            temp2 = temp2 + abs(X(h,n) - meandim(n));
        end
        temp(n)=temp2/ popsize;
    end
    divj=sum(temp(:));
    
    div = divj/D;
    
    if div > maxDiv
        maxDiv = div;
    end
    
    divPercent = div/maxDiv;
    
    %% Diversity initialize
    
    %% Historico de diversidad y valor inicial de W
    iwt = divPercent;
    divIter(1) = iwt;
    
    c1Iter(1)=c;

                 
 	GEN = 1;  
    
   	while (GEN <= maxGEN)
        for i = 1:popsize   
            %  Biogeography-based exemplar selection method
            if stay_num(i)>gapm 
                pBest_ind(i,:) = LearnIndex_BLPSO(val_pBest,popsize,D,i,parmu,lambda);
               	cc = cc+1; 
            	stay_num(i) = 0;
            end
            
          	% Biogeography-based learning updating    
          	for  j=1:D
            	pBest_f(i,j) = pBest(pBest_ind(i,j),j);  
            end  
          	V(i,:) = iwt*V(i,:) + c*rand(1,D).*(pBest_f(i,:)-X(i,:));  % update velocity
           	V(i,:) = boundary_repair(V(i,:),Vmin,Vmax,'absorb');   
           	X(i,:) = X(i,:)+V(i,:);    % update position
            X(i,:) = min(X(i,:),Xmax);
            X(i,:) = max(X(i,:),Xmin);


           	if all(X(i,:)<=Xmax) && all(X(i,:)>=Xmin)  % X(i,:) is feasible
                val = objectiveFcn(X(i,:));
                FAccess = FAccess + 1;
                if val<val_pBest(i)    % update pBest
                    pBest(i,:) = X(i,:);  val_pBest(i) = val;
                  	stay_num(i) = 0;
                else
                  	stay_num(i) = stay_num(i)+1;
                end 
            end 
  
        end  
        [~,indexG] = min(val_pBest);
     	bestInd = pBest(indexG,:); bestFit = val_pBest(indexG,:); 
 
        BestFitIter(GEN) = bestFit;
        BestIndIter(GEN,:) = bestInd;
        FAccessIter(GEN) = FAccess;

        if GEN >= 2 && BestFitIter(GEN) == BestFitIter(GEN-1)
            Cstuck = Cstuck+1;
        else
            Cstuck = 0;
        end

        %% Diversity
        if  limitDiv ~= -1
            temp = zeros(1,D);
            for l=1:D
                temp2 = zeros(popsize,1);
                for h=1:popsize
                    %temp2(h) = particle(h).Position(l);
                    temp2(h) = X(h,l);
                
                end
                temp(l)=median(temp2);
            end
            meandim=temp;
            temp = zeros(1,D);
            for l=1:D
                temp2 = 0;
                for h=1:popsize
                    %temp2 = temp2 + abs(particle(h).Position(l) - meandim(l));
                    temp2 = temp2 + abs(X(h,l) - meandim(l));
                end
                temp(l)=temp2/popsize;
            end
            divj=sum(temp(:));
            
            div = divj/D;
            
            if div > maxDiv
                maxDiv = div;
            end
            
            %divPercent = div/maxDiv*100;
            divPercent = div/maxDiv;
        end

        %% Ajuste de W
        if divPercent <= 0.01
            Vmax = (Xmax - Xmin);  
            Vmin = -Vmax;       
        end
        %else
        
        %inputs = [divPercent*100 FAccess];
        inputs = [min(Cstuck,100) divPercent*100];
        output = evalfis(fis,inputs);
        

        iwt = output(1);
        c=output(2);
       
      
        [Xpl(GEN), Xpt(GEN), Div_max] = EEMeasure(X, dims, popsize, Div_max);

        
        GEN = GEN+1;
    end


    
end
