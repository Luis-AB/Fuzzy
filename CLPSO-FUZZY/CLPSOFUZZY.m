function [bestFit,bestInd,BestFitIter,BestIndIter,FAccessIter] = CLPSOFUZZY(objectiveFcn,iterations,dims,bounds)
    fis = readfis("PSOC_FUZZY.fis");
    Div_max = 0;
    
    %counter1 = 0;
    
    D = dims;            % Number of Decision Variables
    Xmin = bounds(1,1);             % Lower Bound of Decision Variables
    Xmax = bounds(1,2);            % Upper Bound of Decision Variables
    maxGEN = iterations;      % Maximum Number of Iterations
    popsize = 50;        % Population Size (Swarm Size)
    
    BestIndIter = zeros(maxGEN,D);
    BestFitIter = ones(maxGEN,1)*Inf;
    FAccessIter = zeros(maxGEN,1);
    FAccess = 0;

    pBest = zeros(popsize,D);
    val_pBest = ones(popsize,1) * inf;
    pBest_ind = zeros(popsize,D);
    pBest_f = zeros(popsize,D);
    
    %iwt = 0.9 - (1 : maxGEN) * (0.5 / maxGEN);
    c = 1.49445;
    
    %Parametros nuevos
    Cstuck = 0;
    maxDiv = 0;
    limitDiv = 1;
    divChangeCount = 0;
    divIter=zeros(maxGEN,1);
    c1Iter=zeros(maxGEN,1);
    velOLD = zeros(popsize,D);
    velNEW = zeros(popsize,D);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   

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
    %Velocity limits
   	Vmax = (Xmax - Xmin)*0.2;  
    Vmin = -Vmax;    
   	V = repmat(Vmin,popsize,1) + rand(popsize,D).*repmat(Vmax-Vmin,popsize,1);
        
   	% Learning Probability Pc
   	t = 0 : 1/(popsize-1):1;
   	t = 5.*t;
 	Pc = 0.0+(0.5-0.0).*(exp(t)-exp(t(1)))./(exp(t(popsize))-exp(t(1)));
        
  	% Refreshing Gap
 	gapm = 5; stay_num = zeros(1,popsize); 
  	for i=1:popsize
        pBest_ind(i,:) = LearnIndex_CLPSO(val_pBest,popsize,D,i,Pc(i));
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
            % update exemplar index
            if stay_num(i)>gapm   
                pBest_ind(i,:) = LearnIndex_CLPSO(val_pBest,popsize,D,i,Pc(i));
               	stay_num(i) = 0;
            end 

          	%  Compehensive Learning Strategy       
          	for  j=1:D
            	pBest_f(i,j) = pBest(pBest_ind(i,j),j);  
            end  
          	V(i,:) = iwt*V(i,:) + c*rand(1,D).*(pBest_f(i,:)-X(i,:));  % update velocity
            velOLD(i,:) = V(i,:);
           	V(i,:) = boundary_repair(V(i,:),Vmin,Vmax,'absorb');
            velNEW(i,:) = V(i,:);
           	X(i,:) = X(i,:)+V(i,:);    % update position
            Comp(i,:) = velOLD(i,:) == velNEW(i,:);

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
        %Diversity(GEN) = divPercent;

        %% Ajuste de W
        if divPercent <= 0.01

            Vmax = (Xmax - Xmin);  
            Vmin = -Vmax;
            %counter1 = counter1 + 1;
            
        end
        %else
        
        %inputs = [divPercent*100 FAccess];
        inputs = [min(Cstuck,100) divPercent*100];
        output = evalfis(fis,inputs);
        
        
        iwt = output(1);
        c = output(2);
        
        %end
       
 
      GEN = GEN+1;
    end
    
end

