//Change line 26 of the "quickrun" function to define y how you want in order to solve arbitrary problems
//Not currently programmed for a poisson distribution--change the "EvaluatePoint" function to change how the points are defined

//Globals used in the calculation
global gb_rowops
global flipped_rows
global active_set
global gb_used_indices
gb_rowops = eye(1,1);
flipped_rows = []
active_set = []
gb_used_indices = [] //like the active set but elements are never removed

// Globals used as options
global rng_points //0 = 1s and -1s, 1 = 0.9s and -0.9s
rng_points = 1

// Globals for tracking performance and statistics
global gb_stats
global gb_timestats
global gb_costcoefs // notes the reduced cost value of the pivot column in each simplex iteration
global gb_us //counts the progression of the number of columns yet unsolved
global gb_flipcounter
gb_timestats = zeros(1,7); //time via timer() function
gb_stats = zeros(1:7);
gb_costcoefs = []
gb_flipcounter = []
stacksize('max')

//gb_stats(1) = time for make_phase_one_indir
//gb_stats(2) = time for simplex_indir
//gb_stats(3) = time for find_min_cost_coef
//gb_stats(4) = Total iterations of advanced flipping algorithm
//gb_stats(5) = Greatest depth of advanced flipping algorithm
//gb_stats(6) = number of simplex iterations
//gb_stats(7) = total sum of cost coefficients

//gb_timestats(1) = gen_rand_in_point_fix time (Y)
//gb_timestats(2) = make_phase_one_indir time (Y)
//gb_timestats(3) = simplex_indir time
//gb_timestats(4) = pivot-related calculations within simplex indir
//gb_timestats(5) = flipping_alg or fakingAlg
//gb_timestats(6) = flipping_alg_advanced
//gb_timestats(7) = manual_active_set_finder time (Y)

// ### Initializing Functions ###

//Generates and solves a random test problem with a specified N value
function[success,y,time] = quickrun(N)
    //note: time does not include time required to generate the test problems
    global gb_rowops
    global gb_stats
    global active_set
    global gb_used_indices
    global gb_costcoefs
    global gb_us
    global gb_flipcounter
    global gb_timestats
    gb_stats = zeros(1:7);
    gb_timestats = zeros(1:7)
    gb_rowops = 0;
    active_set = [];
    gb_costcoefs = []
    gb_us = []
    gb_flipcounter = []
    gb_used_indices = [1:((N-1)*N/2+1)]'
    timer()
    // y = target point
    y = gen_rand_in_point_fix(N); // <-- CHANGE THIS TO CHANGE THE OBJECTIVE
    Asize = (N-1)*(N)/2 + 1;
    Asizemax = 2^(N-1);
    Asize = min(Asizemax,Asize);
    // A = initial tableau
    A = evaluatePoints(1:Asize,N) // <-- Change evaluatePoints to change the A-matrix
    b = y;
    temp = size(A); temp = temp(1);
    A(temp+1,:) = ones(A(1,:));
    b(temp+1) = 1;
    tic()
    [A2,b2,tab] = make_phase_one_indir(A,b,N,%f);
    if tab == -1 
        success = %f
        time = -1
        return
    end
    [A3,b3,rez] = process_results_indir(A2,b2,N,tab);
    gb_used_indices = unique(gb_used_indices)
    if (abs(rez - y) < 10^(-4))
        //occasionally the basis does end up in the tableau by coincidence
        output = "Solution matches the original problem!";
        success = %t
        disp(output);
    else
        //Manually find active set
        aset = manual_active_set_finder(N,%f);
        asettab = evaluatePoints(aset,N);
        asettab($+1,:) = 1; //for convexity constraints
        for i = 1:(N*(N-1))/2
            asettab(i,:) = asettab(i,:)*flipped_rows(i)
        end
        asettab($+1,:) = 0; //for c row of tableau
        asettab = gb_rowops*asettab; // apply all previous row operations to the newly generated column
        asettab = round_fp(asettab); //fix near-zeros and near-ones caused by floating point
        A2 = asettab(1:$-1,:);
        asettab($+1,:) = aset';
        asettab(:,$+1) = tab(:,$);
        b2 = asettab(1:$-2,$);
        timer();
        [A3,b3,rez] = process_results_indir(A2,b2,N,asettab);
        if (abs(rez - y) < 10^(-4))
            output = "Solution matches the original problem!";
            success = %t
        else
            output = "Solution does not match the original problem, something went wrong!";
            success = %f
            disp(rez)
            disp(y)
            disp(rez - y)
        end
        disp(output)
    end
    time = toc()
endfunction

// Generate a random target point
function[rand_in_pt] = gen_rand_in_point_fix(N)
    global gb_timestats
    timer()
    m = 2^(N-1)
    n = min((N)^2, m)
    k = N*(N-1)/2 //size of points
    //choose the extreme points by index
    if floor(N/32) == 0
        extreme_points = grand(1,n,"uin",1,m)
        s = size(extreme_points)
        su = size(unique(extreme_points))
        while (s(2) <> su(2))
            // This is a slightly lazy implementation of ensuring uniqueness
            // but it works fine
            extreme_points = grand(1,n,"uin",1,m)
            s = size(extreme_points)
            su = size(unique(extreme_points))
        end
    else 
        //Too large to generate uniformly in one go
        //"unf" generation leads to all numbers being 2 or 4 or 8 apart, depending on size of numbers
        //offset by fudge to ensure good distribution of indices
        fudge = 2^(N-32) //this has the same problems past N =64, but this implementation cannot solve past N=52 anyway
        s = [1,4]
        su = [2,3]
        while (s(2) <> su(2))
            extreme_points = grand(1,n,"unf",1+fudge,m-fudge)
            extreme_points = extreme_points + grand(1,n,"uin",-1*fudge,fudge)
            extreme_points = ceil(extreme_points)
            s = size(extreme_points)
            su = size(unique(extreme_points))
        end
    end
    //generate those vectors
    extreme_points = evaluatePoints(extreme_points,N)
    for i = 1:n
        if %t //(rand() > 0.5) //make this lower to put points closer to the center
            rng(i) = rand();
        else
            rng(i) = 0;
        end
    end
    rng = rng/sum(rng);

    rand_in_pt = zeros(k,1);
    for i = 1:n
        rand_in_pt = rand_in_pt + rng(i)*extreme_points(:,i)
        //rng contains a set of numbers with a sum of 1
    end
    gb_timestats(1) = gb_timestats(1) + timer()
endfunction

// Determines if a column contains a single 1 value and all other values are 0
function[result] = isIdentityColumn(col)
    s = size(col); s = s(1);
    result = %t;
    if sum(col) <> 1
        result = %f;
    end
    for i = 1:s
        if col(i) <> 0 & col(i) <> 1
            result = %f;
        end
        if result == %f;
            break
        end
    end
endfunction

// Set up the simplex tableau
// "phase one" is the only phase for this sort of problem 
// so the name is somewhat vestigial
function [A,b,tab] = make_phase_one_indir(A,b,N,debugFlag)
    global gb_rowops
    global gb_stats
    global gb_used_indices
    global gb_timestats
    global flipped_rows
    timer()
    [m,n]              = size(A);
    tab                = zeros(m+1,n+1);
    tab(1:m,1:n)       = A;
    tab($,n+1:$-1)     = 1;
    tab(1:m,$)         = b(:);
    gb_rowops = eye(m+1,m+1)
    gb_stats = zeros(1:7);

    // Account for negative b-values
    flipped_rows = []
    for i = 1:m
        if b(i) < 0
            tab(i,:) = -1*tab(i,:)
            flipped_rows($+1) = -1
        else
            flipped_rows($+1) = 1
        end
    end
    
    //calculate reduced costs and setup gb_rowops
    for i = 1:m 
        tab($,:) = tab($,:)-tab(i,:);
        gb_rowops($,:) = gb_rowops($,:)-gb_rowops(i,:);
    end
    gb_rowops($,$) = 1;

    
    gb_timestats(2) = gb_timestats(2) + timer()
    tab = simplex_indir(tab(1:m,1:n),tab(1:m,$),tab($,1:n),N,debugFlag,'phase one'); //TODO calcfix (medium)
    if tab == -1
        A = -1
        b = -1
        return
    end

    A = tab(1:m,1:$-1);
    b = tab(1:m,$);
    gb_stats(1) = gb_stats(1) + toc()
endfunction

/// ### Main Loop Functions ###
function [tab] = simplex_indir(A,b,c,N,debugFlag,phase_name)
    global gb_rowops
    global gb_stats
    global gb_used_indices
    global gb_costcoefs
    global gb_us
    global gb_timestats
    timer()
    [m,n]        = size(A);
    tab          = zeros(m+1,n+1);
    tab(1:m,1:n) = A;
    tab(m+1,1:n) = c(:)';
    tab(1:m,$)   = b(:);
    tab(m+2,1:n) = zeros(tab(1,1:n))
    tab(m+2,1:n) = 1:(n); //this bottom row tracks the indices of tab entries

    keep_running = %t;
    while keep_running
        gb_stats(6) = gb_stats(6) + 1 //count the number of iterations
        if modulo(gb_stats(6),100) == 0
            unsolved = 0
            for i = 1:((N*(N-1))/2)
                //reduced cost of 0 is a good estimate of which columns are "solved"
                // ie: a column of 0s with a single 1
                if tab($-1,i) <> 0 
                    unsolved = unsolved + 1
                end
            end
            gb_us($+1) = unsolved;
        end
        ////adjust modulo to adjust frequency of updates
        if modulo(gb_stats(6),100) == 0 
            disp(gb_stats(6))
            disp(unsolved,"Unsolved:")
        end

        gb_timestats(3) = gb_timestats(3) + timer()
        [trueindex,tabindex,tab] = approx_min_cost_coef(tab, N)
        if trueindex == -1 | tabindex == -1
            disp("Failure in approx_min_cost_coef")
            tab = -1
            return
        end
        timer()
        if ~and(abs(gb_rowops($,1:$-1)) < 10^-8) //when all bottom-row entries of gb_rowops are 0 (except the farthest right) we're done
            [dummy,J] = min(tab($-1,1:n)); //find the most negative cost coefficient
            //overwrite J for new algorithm
            J = tabindex
            
            gb_stats(7) = gb_stats(7) + tab($,tabindex)
            gb_costcoefs($+1) = dummy
            gb_stats(2) = gb_stats(2) + toc()
            
            if and(tab(1:m,J)<=0)
                disp('ERROR: problem unbounded. All entries <= 0 in column:');
                disp(J);
            else
                pivot_row = 0;
                min_found = %inf;
                for i = 1:m //
                    if tab(i,J)>0
                        tmp = tab(i,$)/tab(i,J); 
                        if tmp < min_found
                            min_found = tmp;
                            pivot_row = i;
                        end
                    end
                end
                
                // normalize the pivot row(s)
                gb_rowops(pivot_row,:) = gb_rowops(pivot_row,:)/tab(pivot_row,J); // apply the constant division to our rowops matrix
                tab(pivot_row,:) = tab(pivot_row,:)/tab(pivot_row,J);

                gb_timestats(3) = gb_timestats(3) + timer()
                for i=1:m+1
                    if i ~= pivot_row
                        gb_rowops(i,:) = gb_rowops(i,:)-sign(tab(i,J))*abs(tab(i,J))*gb_rowops(pivot_row,:);
                        tab(i,:)=tab(i,:)-sign(tab(i,J))*abs(tab(i,J))*tab(pivot_row,:);
                    end
                end
                temp = timer()
                gb_timestats(3) = gb_timestats(3) + temp
                gb_timestats(4) = gb_timestats(4) + temp
            end
        else
            disp("Simplex has completed")
            keep_running=%f;
        end
    end
    
    gb_stats(2) = gb_stats(2) + toc()
    gb_timestats(3) = gb_timestats(3) + timer()
endfunction

// Setup for the flipping algorithm and update the tableau
function [trueindex,tabindex,tab] = approx_min_cost_coef(tab, N)
    global gb_rowops
    global gb_stats
    global gb_timestats
    global active_set
    timer()
    row = gb_rowops($,1:($-2))
    [trueindex,result] = flippingAlg(row,N)
    //[trueindex,result] = fakingAlg(row,N)
    if trueindex == -1 | result == -1
        disp("Failure in flipping algorithm")
        trueindex = -1
        tabindex = -1
        return
    end
    //now to insert this true-index value into the tableau,
    //store its tabindex

    [tab,tabindex] = tabinsert(trueindex,tab,N)
    gb_stats(3) = gb_stats(3) + toc()
endfunction

// For comparison to the baseline revised simplex method
// which simply picks at random until a negative reduced cost is found
// Runs in polynomial time, but slower than flipping_alg
function[index,result] = fakingAlg(x,N)
    global gb_timestats
    global flipped_rows
    timer()
    index = 1
    for i = 1:(N*(N-1)/2) //fix row-flips
        x(i) = x(i)*flipped_rows(i)
    end
    m = 2^(N-1)
    if N >= 32
        fudge = 2^(N-32)
    else
        fudge = 0
    end
    count = 0
    tic()
    while %t
        s= ceil(1000000*rand())//grand's seed gets reset at each evaluatePoints call
        grand("setsd",s)
        count = count + 1
        index = ceil(grand(1,1,"unf",1+fudge,m-fudge-1))
        if fudge <> 0
            index = index + grand(1,1,"uin",-1*fudge,fudge)
        end
        testvec = evaluatePoints(index,N)
        reducedcost = sum(testvec.*x')
        if reducedcost < 0
            //success!
            result = reducedcost
            //index = index
            break
        end
        if modulo(count,100) == 0
            disp("faking alg count",count)
        end
    end
    gb_timestats(5) = gb_timestats(5) + timer()
endfunction

// Finds a resultant vector with a strong (but not necessarily optimal) negative
// reduced cost.
// extremely quick in practice compared to the rest of the simplex iteration
function[index,result] = flippingAlg(x,N)
    global active_set
    global gb_flipcounter
    global flipped_rows
    global rng_points
    global gb_timestats //#5
    timer()
    gb_flipcounter($+1,1) = 0
    gb_flipcounter($,2)   = 0
    //disp("Flipping Alg begun")
    //the name "gameboard" is from the idea of making it a game:
    //press the right buttons to minimize the sum in play
    //not allowed to return an index that is already in the active_set
    //if x = -1, randomize for testing
    index = 0
    if x == -1
        x = grand(1,N*(N-1)/2,"uin",-20,20)'
    end
    s = size(x)
    if s(2) > s(1) //if there are more columns than rows
        x = x' //flip x upright
    end

    for i = 1:(N*(N-1)/2)
        x(i) = x(i)*flipped_rows(i)
    end
    
    gameboard = zeros(N-1,N-1)
    
    //choose a starting point generating vector:         
    //greedy method
    for i = 1:N-1
        if x(i) > 0 //if x is positive, we want negative genvec
            genvec(i) = -1
        else
            genvec(i) = 1
        end
    end
    
    //populate the game board
    for i = 1:(N-1)
        gameboard(i,i) = x(i)
    end
    
    count = N
    for j = 1:(N-2)
        for i = j+1:(N-1)
            gameboard(j,i) = x(count)
            count = count + 1;
        end
    end
    
    //Change the game board based on our starting point (genvec)
    for i = 1:N-1
        if genvec(i) == -1
            gameboard(i,:) = -1*gameboard(i,:)
            gameboard(:,i) = -1*gameboard(:,i)
            gameboard(i,i) = -1*gameboard(i,i)
        end
    end
    
    inp = 1
    oldchange = -1
    while %t
        maindiagonalsum = 0
        for i = 1:(N-1)
            maindiagonalsum = maindiagonalsum + gameboard(i,i)
        end
        //Determine the change caused by each flip
        for i = 1:(N-1)
            changefrom(i) = (-2)*(sum(gameboard(i,:)) + sum(gameboard(:,i)) - gameboard(i,i)) //avoid double-counting i
            changeopposite(i) = changefrom(i) - 2*(maindiagonalsum - 2*gameboard(i,i)) //result of swapping everything except i
        end
        gb_flipcounter($,1) = gb_flipcounter($,1) + 1
        //The double-min in this conditional finds the minimal value between the two vectors
        if (and(changefrom > 0) & and(changeopposite > 0)) | (min(min(changefrom,changeopposite)) == 0 & oldchange == 0) 
            //If no flips are non-positive, prepare to end the loop
            finalsum = sum(gameboard)
            if finalsum < 0 //If we've created a game state that sums to negative, then that's good enough
                break;
            else //Time to get complicated
                //disp("This problem is more complicated--fixes coming soon!")
                gb_timestats(5) = gb_timestats(5) + timer()
                [index,finalsum] = flippingAlgAdvanced(x,N,gameboard,genvec,'p') //issue = positive sum and all swaps are positive
                timer()
                break;
            end
        end
        [mini,inp] = min(changefrom)
        [miniflip,inp2] = min(changeopposite)
        //whichever minimum is smaller is used
        if mini <= miniflip
            oldchange = round_fp(mini); //if oldchange is zero, it affects our choices in the next loop
            //Implement the flip in our game state
            gameboard(inp,:) = -1*gameboard(inp,:)
            gameboard(:,inp) = -1*gameboard(:,inp)
            gameboard(inp,inp) = -1*gameboard(inp, inp)
            genvec(inp) = -1*genvec(inp)
        else //mini > miniflip
            oldchange = round_fp(miniflip); //if oldchange is zero, it affects our choices in the next loop
            inp = inp2
            //Implement the flip in our game state
            gameboard(inp,:) = -1*gameboard(inp,:)
            gameboard(:,inp) = -1*gameboard(:,inp)
            gameboard(inp,inp) = -1*gameboard(inp, inp)
            for i = 1:(N-1)
                gameboard(i,i) = -1*gameboard(i,i)
            end
            genvec = -1*genvec
            genvec(inp) = -1*genvec(inp)
        end
    end
    result = finalsum

    index = genvecTranslate(genvec,N) //need this in both version
    if rng_points == 1
        point = evaluatePoints(index,N)
        s(1,:) = size(point)
        s(2,:) = size(x)
        if s(1,:) == s(2,:)
            result = sum(point.*x)
        else
            result = sum(point.*x')
        end
    end
    gb_flipcounter($,2) = result
    //    disp(output)
    gb_timestats(5) = gb_timestats(5) + timer()
    if or(index==active_set)
        [index,result] = flippingAlgAdvanced(x,N,gameboard,genvec,'d')
    end
endfunction

// When the flipping algorithm fails, engage the advanced flipping algorithm
// the flipping algorithm very rarely fails, so this uses a more robust but more expensive calculation
// This should always terminate when an answer exists, however it may take exponential time in extremely rare cases
// less efficient than gameboard flipping, but it is otherwise difficult to account for an arbitrary number of flips
function[index,result] = flippingAlgAdvanced(x,N,gameboard,genvec,issue)
    global gb_timestats
    global gb_stats
    timer()
    
    loopcount = 0
    //start at 1, but if we're here, it will likely ramp up
    concurrent_flips = 1 
    flipsset = zeros(1,1)
    count = 1
    for i = 1:(N-1)
        flipsset(count,:) = i
    end
    gb_stats(5) = 1
    while %t //check each candidate ("cand") set of tier 0 flips
        gb_stats(4) = gb_stats(4) + 1
        temp = genvecTranslate(genvec)
        temp = evaluatePoints(temp,N)
        minresult = sum(temp.*x)
        minflip = -1
        s = size(flipsset)
        for i = 1:s(1) //iterate through each set of possible flips
            loopcount = loopcount + 1
            candGenvec = genvec
            for j = 1:s(2)
                candGenvec(flipsset(i,j)) = (-1)*candGenvec(flipsset(i,j))
            end
            candIndex = genvecTranslate(candGenvec)
            candPoint = evaluatePoints(candIndex,N)
            candResult = sum(candPoint.*x)
            if candResult < minresult 
                minresult = candResult
                minflipset = i
                for j = 1:s(2)
                    genvec(flipsset(i,j)) = (-1)*genvec(flipsset(i,j))


                end
                // Once an improvement has been made, we can return to 1 concurrent 
                // flips until we reach another impasse
                count = 1
                newflipsset = zeros(1,s(2)+1)
                for k = 1:(N-1)
                    for l = (k+1):(N-1)
                        flipsset(count,:) = [k,l]
                        count = count + 1
                    end
                end
                flipsset = newflipsset
                break
            end
        end
        if minflip == -1 //no more improvements are possible at this depth
            if minresult >= 0
                //Need to add another level to flipsset
                gb_stats(5) = gb_stats(5) + 1
                s = size(flipsset)
                count = 1
                newflipsset = zeros(1,s(2)+1)
                for i = 1:s(1)
                    for j = (flipsset(i,2)+1):(N-1)
                        newflipsset(count,1:s(2)) = flipsset(i,:)
                        newflipsset(count,$) = j
                        count = count + 1
                    end
                end
                flipsset = newflipsset
            else
                break
            end
        end
    end

    //Now fix everything up for the new calculations
    index = genvecTranslate(genvec,N)
    point = evaluatePoints(index,N)
    s(1,:) = size(point)
    s(2,:) = size(x)
    if s(1,:) == s(2,:)
        result = sum(point.*x)
    else
        result = sum(point.*x')
    end
    gb_timestats(6) = gb_timestats(6) + timer()
endfunction


//This replaces an old column of the tableau with the new pivot column
function[tab,tabindex] = tabinsert(trueindex, tab, N)
    global active_set
    global gb_timestats
    global gb_rowops
    global gb_used_indices
    global flipped_rows
    if ~or(find(gb_used_indices==trueindex)) //returns true if trueindex is already in the used indices list
        gb_used_indices($+1) = trueindex
    end
    //if the chosen index is already in the tableau, find it and return its index
    if or(tab($,:) == trueindex)
        s = size(tab)
        for i = 1:s(2)
            if tab($,i) == trueindex
                tabindex = i;
                break
                if ~or(active_set == trueindex) 
                    //if the active set does not contain the true index add it
                    //but don't replace one
                    active_set($+1) = trueindex
                end
            end
        end
        return
    end
    //Takes a tableau and a true-index, multiplies the true-index's vector by
    //gb_rowops and replaces an old tableau item with it
    dim = (N*(N-1)/2)+1 //size of the active set
    vec = evaluatePoint(trueindex,N)
    //Apply flipped rows
    vec = vec.*flipped_rows(1:$-1) //apply row-flips to the generated vector
    vec($+1,:) = 1; //for convexity constraints

    vec($+1,:) = 0; //for c row of tableau
    vec = gb_rowops*vec // apply all previous row operations to the newly generated column
    vec = round_fp(vec) //remove near-zeros caused by floating point
    vec($+1,:) = trueindex'; //for indirect indexing
    //Now to determine which column of the tab to replace

    //Simple removal: max cost coefficient is removed
    //find the max cost coefficient
    [maxv,tabindex] = max(tab($-1,1:dim))
    //remove this cost coefficient from the active set if it is present
    removal = tab($,tabindex)
    s = size(active_set)
    if or(removal == active_set) //if the index for removal is in the active set... remove it
        for i = 1:s(1)
            if active_set(i) == removal
                active_set(i) = trueindex
                output = "updated active set, " + string(removal) + " replaced by " + string(trueindex) + "."
                //                disp(active_set,output)
                break
            end
        end
    else //if the tab to be removed is not present, just add the new one
        active_set($+1) = trueindex
        output = "updated active set, " + string(trueindex) + " added."
        //        disp(active_set,"updated active set (addition):")
    end
    //replace the chosen column in the tableau
    tab(:,tabindex) = vec;
endfunction

// ################################
// ### Post-processing Functions###
// ################################

// Constructs the final basis, which does not always align with the final tableau
function[aset] = manual_active_set_finder(N, fullsearch)
    global gb_rowops
    global gb_used_indices
    global gb_timestats
    timer()
    aset = []
    if fullsearch
        for i = 1:2^(N-1)
            temp = evaluatePoint(i,N)
            temp($+1) = 1 //for convexity constraint
            temp = temp.*flipped_rows
            temp($+1) = 0 //for bottom row of tableau, cost coefficient
            potential = gb_rowops*temp
            potential = round_fp(potential)
            if ((sum(potential) == 1) & (min(potential) == 0) & (max(potential) == 1))
                aset($+1) = i
            end
        end
    else
        s = size(gb_used_indices)
        for j = s(1):-1:1 //check each index that was ever in the active set
            //start at the end
            i = gb_used_indices(j)
            temp = evaluatePoint(i,N)
            temp($+1) = 1 //for convexity constraint
            temp = temp.*flipped_rows
            temp($+1) = 0 //for bottom row of tableau, cost coefficient
            potential = gb_rowops*temp
            potential = round_fp(potential)
            if ((sum(potential) == 1) & (min(potential) == 0) & (max(potential) == 1))
                aset($+1) = i
            end
            s2 = size(aset)
            //if s2(1) == ??
            //break
            //end
        end
    end
    gb_timestats(7) = gb_timestats(7) + timer()
endfunction

// Confirms the validity of the result by calculating the convex combination
// represented by the tableau
// b3 gives the weights, A3 gives the resultant vectors
// such that b3(i) is the weight for A3(:,i)
function[A3, b3, combination] = process_results_indir(A2, b2, N, tab)
    combination = 0;
    s = size(A2);
    m = s(1); n = s(2); // m = rows; n = columns
    count = 1;
    disp("A2:")
    disp(A2)
    for i = 1:n
       if isIdentityColumn(A2(:,i))
           if tab($,i) == 0
               disp("Error, incomplete operation [process_results_indir]")
           end
           bindex = find(A2(:,i)) //gives the index of the "1"
           colindex = tab($,i)
           b3(count) = b2(bindex);
           if colindex == 0
               disp("ERROR SOLUTION NOT VALID")
               colindex = 1
           end
           A3(:,count) = evaluatePoint(colindex,N)
           combination = combination + b2(bindex)*A3(:,count)
           count = count + 1;
       end
    end
endfunction

// ########################
// ### Helper Functions ###
// ########################

// Transform an index into the resultant vector it represents
// This problem remains solvable for a wide variety of point evaluation algorithms
// so long as the high and low values are sufficiently distinguished
// and the calculation is deterministic
function [point] = evaluatePoint(index, N)
    global rng_points
    //    returns the value of the point at the given index
    //    for a given N, the size of the determining vectors
    //determining vector at index:
    fudgefactor = 0.2 // difference between results and 1/-1 setup
    //default = 0.2 fits the data
    detvec(1) = 1
    detvec(2:N) = decToBinVec(index-1,N-1)
    if rng_points == 0
        point = zeros((N*(N-1))/2,1)
        count = 0;
        for i = 1:(N-1)
            for j = (i+1):N
                count = count + 1;
                if detvec(i) == detvec(j)
                    point(count) = 1;
                else
                    point(count) = -1;
                end
            end
        end
    elseif rng_points == 1
        // set a unique, but unchanging seed for each index/size
        grand('setsd',modulo(index,N))
        point = zeros((N*(N-1))/2,1)
        count = 0;
        for i = 1:(N-1)
            for j = (i+1):N
                count = count + 1;
                if detvec(i) == detvec(j)
                    point(count) = (1 - fudgefactor) + grand(1, 1, "unf", 0, fudgefactor);
                else
                    point(count) = -(1 - fudgefactor) - grand(1, 1, "unf", 0, fudgefactor);
                end
            end
        end
    end
endfunction

// Runs evaluatePoint on an arbitrary number of indices, creating a matrix
// with each index's resultant vector as sequential columns
function[points] = evaluatePoints(indices,N)
    s = size(indices); s = prod(s);
    for i = 1:s
        points(:,i) = evaluatePoint(indices(i),N)
    end
endfunction

function[binVec] = decToBinVec(x,y)
    // Returns a y-bit binary, little-endian column vector 
    // containing the binary representation of x
    binVec = zeros(y,1);
    if(x < 0) //or (log2(x) > y)
        return -1;
    end
    remain = x;
    place = 2^y;
    for i = 1:y
        place = place/2;
        if (remain >= place)
            remain = remain - place;
            binVec(i) = 1;
        end
    end
endfunction

// Translates a generating vector to the resultant vector index it represents
function[index] = genvecTranslate(genvec,N)
    genvec = flipdim(genvec,1)
    index = 0
    for i = 1:N-1
        if genvec(i) == 1
            index = index + 2^(i-1)
        end
    end
    index = index + 1
endfunction

// sets near-0 and near-1 values to 0 and 1 to make sure that 
// floating points don't misbehave
function[newarray] = round_fp(array)
    newarray = array;
    s = size(array)
    for i = 1:prod(s)
        if abs(array(i)) < 10^(-6)
            newarray(i) = 0;
        elseif array(i) < 1 + 10^(-6) & array(i) > 1 - 10^(-6) 
            newarray(i) = 1;
        end
    end
endfunction
