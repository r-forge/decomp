rank.wt <-
function(x,wt)
    {
    # calculate the length and reserve an empty rank vector of the right length
    n <- length(x)
    r <- vector(length = n)
    
    # order the variable
    myOrder <- order(x)
    
    # rescale sampling weights such that sum(weights) = 1 (Lerman & Yitzhaki (1989) par. 2)
    # this is done already when rank.wt is called from RCI, but does not change anything
    wt <- wt / sum(wt)
    
    # calculate the fractional rank and return it in the original order of the variable
    # in the first term, we use a modified vector of weights starting with 0 and the last value chopped off
    # see Lerman & Yitzhaki eq. (2)
    # a loop is avoided by using the cumulative sum
    r[myOrder] <- (c(0,cumsum(wt[myOrder][-n])) + wt[myOrder]/2)
    return(r)
    }

