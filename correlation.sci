// Zoë MacDonald - May 2018
// University of Toronto Master of Computer Science

//Change line 26 of the "quickrun" function to define y to the objective correlation
//Replace the "EvaluatePoint" function to change how the points are defined

//Globals used in the calculation
global gb_rowops
global flipped_rows
gb_rowops = eye(1,1);
flipped_rows = []

// Globals used as options
global rng_points //0 = 1s and -1s, 1 = 0.9s and -0.9s
rng_points = 1
global flipping_mode //decides which flipping mode to follow for iterations
//1 = exact calculation (too expensive past N = 9)
//2 = flipping alg, 3 = advanced flipping alg, 4 = random selection
//5,6,7 are for comparisons between methods
flipping_mode = 2

// Globals for tracking performance and statistics
global gb_stats
global gb_timestats
global gb_us //counts the progression of the number of columns yet unsolved
global gb_flipcounter
global reducedcosts //holds all reduced costs for each iteration (very expensive, only for efficiency testing)
global rccount //indexes the next reducedcosts matrix
global fliptests
gb_timestats = zeros(1,6); //time via timer() function
gb_stats = zeros(1:7);
gb_flipcounter = []
stacksize('max') //for scilab versions before 6.0's dynamic stack size

//Some stats may be deprecated

//gb_stats(1) = time for make_phase_one_indir
//gb_stats(2) = time for simplex_indir
//gb_stats(3) = time for find_min_cost_coef
//gb_stats(4) = Total iterations of advanced flipping algorithm
//gb_stats(5) = Greatest depth of advanced flipping algorithm
//gb_stats(6) = number of simplex iterations
//gb_stats(7) = total sum of cost coefficients

//gb_timestats(1) = gen_rand_in_point_fix time 
//gb_timestats(2) = make_phase_one_indir time
//gb_timestats(3) = simplex_indir time
//gb_timestats(4) = pivot-related calculations within simplex indir
//gb_timestats(5) = flipping_alg or fakingAlg
//gb_timestats(6) = flipping_alg_advanced

// ### Initializing Functions ###

//Generates and solves a random test problem with a specified N value
function[success,y,time] = quickrun(N)
    //note: time does not include time required to generate random test problems
    global gb_rowops
    global gb_stats
    global gb_us
    global gb_flipcounter
    global gb_timestats
    global reducedcosts //holds all reduced costs for each iteration (very expensive, only for efficiency testing)
    global rccount //indexes the next reducedcosts matrix
    global fliptests
    gb_stats = zeros(1:7);
    gb_timestats = zeros(1:7)
    gb_rowops = 0;
    gb_us = []
    gb_flipcounter = []
    if N < 10 //for flipping_mode
        reducedcosts = zeros(2^(N-1),2,2);
        rccount = 1;
    end
    fliptests = zeros(10,3)
    timer()
    // y = vectorized objective correlation matrix
    y = gen_rand_in_point_fix(N); // <-- CHANGE THIS TO CHANGE THE OBJECTIVE
    Asize = (N-1)*(N)/2 + 1;
    Asizemax = 2^(N-1);
    Asize = min(Asizemax,Asize);
    // E = initial feasible basis of artificial variables
    E = eye(zeros(Asize,Asize)) //(this includes the row which would be convexity constraints in A)
    b = y;
    //temp = size(A);
    temp = size(E);
    temp = temp(1);
    b(temp) = 1; //E-version
    tic()
    //[A2,b2,tab] = make_phase_one_indir(A,b,N,%f);
    [A2,b2,tab] = make_phase_one_indir(E,b,N,%f);
    if tab == -1 
        success = %f
        time = -1
        return
    end
    [A3,b3,rez] = process_results_indir(A2,b2,N,tab);
    if (abs(rez - y) < 10^(-4))
        //the final tableau represents the basis
        output = "This solution matches the original problem!";
        success = %t
        disp(output);
    else
        output = "Unexpected Error, final basis is incorrect."
        disp(output)
    end
    time = toc()
endfunction

// Generate a random target point
function[rand_in_pt] = gen_rand_in_point_fix(N)
    global gb_timestats
    timer()
    m = 2^(N-1) //range of indices of possible points
    n = N^2 //amount we want to select to randomly generate the interior point
    //n = min((N)^2, m)
    k = N*(N-1)/2 //size of points
    //choose the extreme points by index
    if m <= n //there are not enough points for our usual scheme, so just use them all
        extreme_points = 1:m
        n = m
    elseif N <= 15
        x = 1:m //this gets expensive at higher m-values
        extreme_points = samwr(n,1,x)
    elseif floor(N/32) == 0
        extreme_points = grand(1,n,"uin",1,m)
        s = size(extreme_points)
        su = size(unique(extreme_points))
        while (s(2) <> su(2))
            // This is a slightly lazy implementation of ensuring uniqueness
            // but it works fine for these larger N values
            extreme_points = grand(1,n,"uin",1,m)
            s = size(extreme_points)
            su = size(unique(extreme_points))
        end
    else 
        //Too large to generate uniformly in one go
        //"unf" generation at this size leads to all numbers being 2 or 4 or 8 apart, depending on size of numbers
        //due to floating point precision
        //must offset by fudge-factor to ensure good distribution of indices
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
        //tab($,:) = tab($,:)-tab(i,:); //enable for version A, disable for E
        gb_rowops($,:) = gb_rowops($,:)-gb_rowops(i,:);
    end
    gb_rowops($,$) = 1;
    tab($,:) = 0 //enable for E, disable for A
    
    gb_timestats(2) = gb_timestats(2) + timer()
    tab = simplex_indir(tab(1:m,1:n),tab(1:m,$),tab($,1:n),N,debugFlag,'phase one');
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
    global gb_us
    global gb_timestats
    timer()
    [m,n]        = size(A);
    tab          = zeros(m+1,n+1);
    tab(1:m,1:n) = A;
    tab(m+1,1:n) = c(:)';
    tab(1:m,$)   = b(:);
    tab(m+2,1:n) = zeros(tab(1,1:n))
    tab(m+2,1:n) = 0 //A version: 1:(n); //this bottom row tracks the indices of tab entries 
    //(E matrix doesn't need this row, so initially set to zero)
    keep_running = %t;
    continueflag = %t;
    while keep_running
        gb_stats(6) = gb_stats(6) + 1 //count the number of iterations
        
        if modulo(gb_stats(6),100) == 0 
            disp(gb_stats(6))
            unsolved = 0
            for i = 1:n
                if tab($,i) == 0 //if the index is zero, this is still a column of E and not A
                    unsolved = unsolved + 1
                end
            end
            if unsolved == 0 //should have exited, what is going on?
                disp("error in gb_rowops, alternative completion conditions activated")
                continueflag = %f
            end
            disp(unsolved,"Unsolved:")
        end
        
        if (~and(tab($,1:$-1) <> 0)) //when no indices from E remain, we're done
            gb_timestats(3) = gb_timestats(3) + timer()
            [trueindex,tabindex,tab] = approx_min_cost_coef(tab, N)
            if trueindex == -1 | tabindex == -1
                disp("Failure in approx_min_cost_coef")
                tab = -1
                return
            end
            timer()

            //overwrite J for new algorithm
            J = tabindex

            gb_stats(7) = gb_stats(7) + tab($,tabindex)
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
                //tab(pivot_row,:) = tab(pivot_row,:)/tab(pivot_row,J); //update whole column
                tab(pivot_row,$) = tab(pivot_row,$)/tab(pivot_row,J); //update b
                tab(pivot_row,J) = 1; //update pivot column

                gb_timestats(3) = gb_timestats(3) + timer()
                for i=1:m+1
                    if i ~= pivot_row
                        gb_rowops(i,:) = gb_rowops(i,:)-tab(i,J)*gb_rowops(pivot_row,:);
                        //tab(i,:)=tab(i,:)-sign(tab(i,J))*abs(tab(i,J))*tab(pivot_row,:); //update entire tableau row
                        tab(i,$)=tab(i,$)-tab(i,J)*tab(pivot_row,$); //update b
                        tab(i,J)=tab(i,J)-tab(i,J)*tab(pivot_row,J); //update pivot column
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
    global reducedcosts
    global rccount
    global fliptests
    global flipping_mode
    timer()
    row = gb_rowops($,1:($-2))
    //trueindex = zeros(4,1)
    result = zeros(4,1)
    s = size(row)
    if s(2) > s(1) //if there are more columns than rows
        row = row' //flip x upright
    end
    //Compare different flipping algorithm variants
    //default (and most efficient) is flipping_mode = 2
    if flipping_mode == 1
        [trueindex,result(1)] = flippingExact(row,N)
    elseif flipping_mode == 2
        [trueindex,result(2)] = flippingAlg(row,N)
    elseif flipping_mode == 3
        [trueindex,result(3)] = flippingAlgAdvanced(row,N)
    elseif flipping_mode == 4
        [trueindex,result(4)] = fakingAlg(row,N)
    elseif flipping_mode == 5 //compare flipping v. adv flipping
        [trueindex,result(3)] = flippingAlgAdvanced(row,N)
        [trueindex,result(2)] = flippingAlg(row,N)
    elseif flipping_mode == 6 //compare flipping, adv flipping, and random
        [trueindex,result(3)] = flippingAlgAdvanced(row,N)
        [trueindex,result(4)] = fakingAlg(row,N)
        [trueindex,result(2)] = flippingAlg(row,N)
    elseif flipping_mode == 7 //compare to exact, print out ratings for flip and advflip
        [trueindex(2),result(2)] = flippingAlg(row,N)
        [trueindex(3),result(3)] = flippingAlgAdvanced(row,N)
        [trueindex(1),result(1)] = flippingExact(row,N)
        truevalues = reducedcosts(:,:,rccount-1)
        testresults = zeros(1,3)
        for i = 1:3
            testresults(i) = find(truevalues(:,2)==trueindex(i))
        end
        fliptests(rccount-1,:) = testresults
        disp(fliptests)
    end
    trueindex = trueindex(1)
    if trueindex == -1 | result(2) == -1
        disp("Failure in flipping algorithm")
        trueindex = -1
        tabindex = -1
        return
    end
    //now to insert this true-index value into the tableau,
    //store its tabindex
    //set trueindex index to change algorithm in use
    [tab,tabindex] = tabinsert(trueindex,tab,N)
    gb_stats(3) = gb_stats(3) + toc()
endfunction

// checks AND SORTS *EVERY* A-column to find the minimum reduced cost
// only for testing, runs *VERY* slowly
function[index,result] = flippingExact(x,N)
    global reducedcosts
    global rccount
    global flipped_rows
    
    s = size(x)
    if s(2) > s(1) //if there are more columns than rows
        x = x' //flip x upright
    end
    for i = 1:(N*(N-1)/2)
        x(i) = x(i)*flipped_rows(i)
    end
    result = %inf
    rvec = zeros(2^(N-1),2)
    rvec(:,2) = [1:2^(N-1)]'
    for i = 1:2^(N-1)
        candpoint = evaluatePoints(i,N)
        candresult = sum(candpoint.*x)
        rvec(i,1) = candresult
        if candresult < result
            result = candresult
            index = i
        end
    end
    reducedcosts(:,:,rccount) = gsort(rvec,'lr','i')
    rccount = rccount + 1
endfunction

// For comparison to the baseline revised simplex method
// which simply picks at random until a negative reduced cost is found
// Runs in polynomial time, but slower than flipping_alg
function[index,result] = fakingAlg(x,N)
    global gb_timestats
    global flipped_rows
    timer()
    x = x'
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
    while %t
        s= ceil(1000000*rand())//grand's seed gets reset at each evaluatePoints call
        grand("setsd",s)
        count = count + 1
        index = ceil(grand(1,1,"unf",0+fudge,m-fudge))
        if fudge <> 0
            index = index + grand(1,1,"uin",-1*fudge,fudge)
        end
        if count > 1000
            disp(index)
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
    global gb_flipcounter
    global flipped_rows
    global rng_points
    global gb_timestats
    timer()
    gb_flipcounter($+1,1) = 0
    gb_flipcounter($,2)   = 0
    
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
            changeopposite(i) = changefrom(i) - 2*(maindiagonalsum - 2*gameboard(i,i)) //bad luck protection--remove to run slightly faster
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
                [index,finalsum] = flippingAlgAdvanced(x,N) //issue = positive sum and all swaps are positive
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
endfunction

// When the flipping algorithm fails, engage the advanced flipping algorithm
// the flipping algorithm very rarely fails, so this uses a more robust but more expensive calculation
// This should always terminate when an answer exists, however it may take exponential time in extremely rare cases
// less efficient than gameboard flipping, but it is otherwise difficult to account for an arbitrary number of flips
function[index,result] = flippingAlgAdvanced(x,N)
    global gb_timestats
    global gb_stats
    global flipped_rows
    timer()

    for i = 1:(N*(N-1)/2)
        x(i) = x(i)*flipped_rows(i)
    end

    for i = 1:N-1
        if x(i) > 0 //if x is positive, we want negative genvec
            genvec(i) = -1
        else
            genvec(i) = 1
        end
    end

    loopcount = 0
    //start at 1, but if we're here, it will likely ramp up
    concurrent_flips = 1 
    flipsset = zeros(1,1)
    count = 1
    for i = 1:(N-1)
        flipsset(i,1) = i
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
            candgenvec = genvec
            for j = 1:s(2)
                candgenvec(flipsset(i,j)) = (-1)*candgenvec(flipsset(i,j))
                if s(2) > 1
                    disp("candgenvec,flipsset,i,j:")
                    disp(candgenvec)
                    disp(flipsset)
                    disp(i)
                    disp(j)
                end
            end
            candindex = genvecTranslate(candgenvec)
            candpoint = evaluatePoints(candindex,N)
            candresult = sum(candpoint.*x)
            if candresult < minresult
                minresult = candresult
                minflip = i
            end
        end
        if minflip > 0
            for j = 1:s(2)
                genvec(flipsset(minflip,j)) = (-1)*genvec(flipsset(minflip,j))
            end
            flipsset = [1:(N-1)]'
            s = size(flipsset)
        end
        if minflip == -1 //no more improvements are possible at this depth
            if minresult >= 0
                //Need to add another level to flipsset
                gb_stats(5) = gb_stats(5) + 1
                s = size(flipsset)
                count = 1
                newflipsset = zeros(1,s(2)+1)
                for i = 1:s(1)
//                    disp("flipsset:")
//                    disp(flipsset)
//                    disp(i,"i:")
//                    disp(minresult,"minresult:")
//                    disp(s,"s:")
//                    disp(count,"count:")
                    for j = (flipsset(i,$)+1):(N-1)
                        newflipsset(count,1:s(2)) = flipsset(i,:)
                        newflipsset(count,$) = j
                        count = count + 1
                    end
                end
                flipsset = newflipsset
                s = size(flipsset)
                disp(flipsset,"new flipsset:")
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
    global gb_timestats
    global gb_rowops
    global flipped_rows
    s = size(tab)
    m = s(1)-2 //unsure
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
    //Pivot removal: determine the pivot row and replace the column which has that as its pivot row
    pivot_row = 0;
    min_found = %inf;
    for i = 1:m //
        if vec(i)>0
            tmp = tab(i,$)/vec(i); 
            if tmp < min_found
                min_found = tmp;
                pivot_row = i;
            end
        end
    end
    tabindex = pivot_row
    tab(:,tabindex) = vec;
endfunction

// ################################
// ### Post-processing Functions###
// ################################

// Confirms the validity of the result by calculating the convex combination
// represented by the tableau
// b3 gives the weights, A3 gives the resultant vectors
// such that b3(i) is the weight for A3(:,i)
function[A3, b3, combination] = process_results_indir(A2, b2, N, tab)
    combination = 0;
    s = size(A2);
    m = s(1); n = s(2); // m = rows; n = columns
    count = 1;
    //disp("A2:")
    //disp(A2)
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
