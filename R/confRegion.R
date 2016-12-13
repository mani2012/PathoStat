# Function to calculate confidence interval and returns the margin of error
confinterval <- function(x, size, delta=0.0001) {
    x[which(x <= 0)] <- delta
    information <- matrix(c((1/x[1])+(1/x[3]), 1/x[3], 1/x[3], 
        (1/x[2])+(1/x[3])),2,2)
    information <- size*information
    variance <- solve(information)
    marginerror <- c(1.96*sqrt(variance[1,1]), 1.96*sqrt(variance[2,2]))
    return (marginerror)
} 

# Function to check whether the given point is with in confidence region
checkregion <- function(x1, x2, chisqval,x,information) {
    d1 <- x1-x[1]
    d2 <- x2-x[2]
    chi <- (d1*information[1,1]+d2*information[2,1])*d1 + 
        (d1*information[1,2]+d2*information[2,2])*d2 
    check <- (chi < chisqval)
    return (check)
}

# Function to calculate confidence interval and returns the margin of error
logitconfinterval <- function(x, size, delta=0.0001) {
    x[which(x <= 0)] <- delta
    information <- matrix(c((1/x[1])+(1/x[3]), 1/x[3], 1/x[3], 
        (1/x[2])+(1/x[3])),2,2)
    information <- size*information
    variance <- solve(information)
    marginerror <- c(1.96*sqrt(variance[1,1])/(x[1]*(1-x[1])), 
        1.96*sqrt(variance[2,2])/(x[2]*(1-x[2])))
    return (marginerror)
} 

# Function to check whether the given point is with in logit confidence region
logitcheckregion <- function(x1, x2, chisqval,x,information) {
    lx <- logit(x)
    d1 <- (x1-lx[1])*x[1]*(1-x[1])
    d2 <- (x2-lx[2])*x[2]*(1-x[2])
    chi <- (d1*information[1,1]+d2*information[2,1])*d1 + 
        (d1*information[1,2]+d2*information[2,2])*d2 
    check <- (chi < chisqval)
    return (check)
}

#' Compute the confidence region for the given proportions
#'
#' @param p1 Read counts for first taxon
#' @param p2 Read counts for second taxon
#' @param size Total read counts in the sample
#' @param uselogit Use logit transformation to compute confidence region
#' @param n Total number of simulation points to generate
#' @param seed Seed to use in random simulation
#' @param jit jitter option (FALSe by default) for the plot
#' @return Confidence region plot
#' @import stats graphics
#' @importFrom gtools logit
#' @export
#' @examples
#' p1 <- 20
#' p2 <- 25
#' size <- 200
#' plotConfRegion(p1, p2, size, uselogit=FALSE)
plotConfRegion <- 
    function(p1, p2, size=100, uselogit=TRUE, n=10000, seed=1000, jit=FALSE)
{
    if (p1 <= 0) p1 <- 1
    if (p2 <= 0) p2 <- 1
    if (size <= (p1+p2)) size <- (p1+p2)+1
    actualprop = c(p1,p2,(size-(p1+p2)))
    set.seed(seed, kind = NULL, normal.kind = NULL)
    A = rmultinom(n,size,actualprop)
    A = A/size
    amount <- 1/(5*size)
    if (jit)  {
        jitA <- jitter(A,1,amount)
        jitA <- abs(jitA)
        jitA <- apply(jitA, c(1, 2), function(x) {if(x > 1) 2-x else x})
    } else  {
        jitA <- A
        jitA <- apply(jitA, c(1, 2), function(x, amount) 
            {if(x <= 0) amount else x}, amount)
    }
    
    # Calculating interval based on averages from all points here
    # me <- apply(A, 2, confinterval, size)
    # p1up <- sum(A[1,] + me[,1])/n
    # p1low <- sum(A[1,] - me[,1])/n
    # p2up <- sum(A[2,] + me[,2])/n
    # p2low <- sum(A[2,] - me[,2])/n
    
    # Calculating interval based on the exact value
    x <- actualprop/size
    me <- confinterval(x, size)
    p1up <- x[1] + me[1]
    p1low <- x[1] - me[1]
    p2up <- x[2] + me[2]
    p2low <- x[2] - me[2]
    
    chisqval <- qchisq(0.95, 2, ncp = 0, lower.tail = TRUE, log.p = FALSE)
    # x <- rep(0,3)
    # x[1] <- (p1up+p1low)/2
    # x[2] <- (p2up+p2low)/2
    # x[3] <- 1-x[1]-x[2]
    information <- matrix(c((1/x[1])+(1/x[3]), 1/x[3], 1/x[3], 
                            (1/x[2])+(1/x[3])),2,2)
    information <- size*information
    
    if (!uselogit)  {
        plot(jitA[1,],jitA[2,], xlab="Proportion 1", ylab="Proportion 2")
        #plot(A[1,],A[2,], xlab="Proportion 1", ylab="Proportion 2")
        for (i in seq_len(n)) {
            if (checkregion(jitA[1,i],jitA[2,i],chisqval,x,information))  {
                points(jitA[1,i],jitA[2,i], col="green")
            }
        }
        
        abline(v=p1up, col="red")
        abline(v=p1low, col="red")
        abline(h=p2up, col="red")
        abline(h=p2low, col="red")
        
        # Create a title with a red, bold/italic font
        title(main="Two Proportions with 95% Confidence Interval", 
            col.main="red", font.main=4)
    } else {
        B <- jitA
        B <- logit(B)
        
        # Calculating interval based on averages from all points here
        # lme <- apply(A, 2, logitconfinterval, size)
        # lp1up <- sum(B[1,] + lme[,1])/n
        # lp1low <- sum(B[1,] - lme[,1])/n
        # lp2up <- sum(B[2,] + lme[,2])/n
        # lp2low <- sum(B[2,] - lme[,2])/n
        
        # Calculating interval based on the exact value
        lx <- logit(x)
        lme <- logitconfinterval(x, size)
        lp1up <- lx[1] + lme[1]
        lp1low <- lx[1] - lme[1]
        lp2up <- lx[2] + lme[2]
        lp2low <- lx[2] - lme[2]
        
        plot(jitA[1,],jitA[2,], xlab="Proportion 1", ylab="Proportion 2")
        for (i in seq_len(n)) {
            if (logitcheckregion(B[1,i],B[2,i],chisqval,x,information))  {
                points(jitA[1,i],jitA[2,i], col="green")
            }
        }
        ilp1up <- exp(lp1up)/(1+exp(lp1up))
        ilp1low <- exp(lp1low)/(1+exp(lp1low))
        ilp2up <- exp(lp2up)/(1+exp(lp2up))
        ilp2low <- exp(lp2low)/(1+exp(lp2low))
        abline(v=ilp1up, col="red")
        abline(v=ilp1low, col="red")
        abline(h=ilp2up, col="red")
        abline(h=ilp2low, col="red")
        # Create a title with a red, bold/italic font
        title(main="Inverse Logit Two Proportions with 95% Confidence Interval",
            col.main="red", font.main=4)
    }

}
