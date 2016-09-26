library(gtools)

# Function to calculate confidence interval and returns the margin of error
confinterval <- function(x, size) {
    delta <- 0.0001
    for (i in 1:3)  {
        if (x[i] <= 0)  {
            x[i] = delta
        }
    }
    information <- matrix(c((1/x[1])+(1/x[3]), 1/x[3], 1/x[3], 
        (1/x[2])+(1/x[3])),2,2)
    information <- size*information
    variance <- solve(information)
    marginerror <- c(1.96*sqrt(variance[1,1]), 1.96*sqrt(variance[2,2]))
    #  marginerror <- c(1.645*sqrt(variance[1,1]), 1.645*sqrt(variance[2,2]))
    return (marginerror)
} 

# Function to check whether the given point is with in confidence region
checkregion <- function(x1, x2, chisqval,x,information) {
    d1 <- x1-x[1]
    d2 <- x2-x[2]
    chi <- (d1*information[1,1]+d2*information[2,1])*d1 + 
        (d1*information[1,2]+d2*information[2,2])*d2 
    check <- FALSE
    if (chi < chisqval) {
        check <- TRUE
    }
    return (check)
} 

# Function to calculate confidence interval and returns the margin of error
logitconfinterval <- function(x, size) {
    delta <- 0.0001
    for (i in 1:3)  {
        if (x[i] <= 0)  {
            x[i] = delta
        }
    }
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
    check <- FALSE
    if (chi < chisqval) {
        check <- TRUE
    }
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
    function(p1, p2, size=100, uselogit=TRUE, n=10000, seed=1000)
{

    actualprop = c(p1,p2,(size-(p1+p2)))
    #actualprop = c(1,1,18)
    set.seed(seed, kind = NULL, normal.kind = NULL)
    A = rmultinom(n,size,actualprop)
    A = A/size
    jitA <- jitter(A)
    p1up <- 0
    p1low <- 0
    p2up <- 0
    p2low <- 0
    for (i in 1:n) {
        me <- confinterval(A[,i], size)
        p1up <- p1up + A[,i][1]+me[1]
        p1low <- p1low + A[,i][1]-me[1]
        p2up <- p2up + A[,i][2]+me[2]
        p2low <- p2low + A[,i][2]-me[2]
    }
    p1up <- p1up/n
    p1low <- p1low/n
    p2up <- p2up/n
    p2low <- p2low/n
    chisqval <- qchisq(0.95, 2, ncp = 0, lower.tail = TRUE, log.p = FALSE)
    x <- rep(0,3)
    x[1] <- (p1up+p1low)/2
    x[2] <- (p2up+p2low)/2
    x[3] <- 1-x[1]-x[2]
    information <- matrix(c((1/x[1])+(1/x[3]), 1/x[3], 1/x[3], 
                            (1/x[2])+(1/x[3])),2,2)
    information <- size*information
    
    if (!uselogit)  {
        plot(jitA[1,],jitA[2,], xlab="Proportion 1", ylab="Proportion 2")
        #plot(A[1,],A[2,], xlab="Proportion 1", ylab="Proportion 2")
        for (i in 1:n) {
            if (checkregion(A[1,i],A[2,i],chisqval,x,information))  {
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
        #quantile (A[1,],.95)
        #quantile (A[1,],.05)
        #quantile (A[2,],.95)
        #quantile (A[2,],.05)
        
        delta = 0.001
        B = logit(A+delta)
        #B = logit(A)
        
        lp1up <- 0
        lp1low <- 0
        lp2up <- 0
        lp2low <- 0
        for (i in 1:n) {
            lme <- logitconfinterval(A[,i], size)
            lp1up <- lp1up + B[,i][1]+lme[1]
            lp1low <- lp1low + B[,i][1]-lme[1]
            lp2up <- lp2up + B[,i][2]+lme[2]
            lp2low <- lp2low + B[,i][2]-lme[2]
        }
        lp1up <- lp1up/n
        lp1low <- lp1low/n
        lp2up <- lp2up/n
        lp2low <- lp2low/n
        jitB <- jitter(B)
        
        #plot(jitB[1,],jitB[2,], xlab="logit(Proportion 1)", 
        #    ylab="logit(Proportion 2)")

        #for (i in 1:n) {
        #    if (logitcheckregion(B[1,i],B[2,i],chisqval,x,information))  {
        #        points(jitB[1,i],jitB[2,i], col="green")
        #    }
        #}
        
        #abline(v=lp1up, col="red")
        #abline(v=lp1low, col="red")
        #abline(h=lp2up, col="red")
        #abline(h=lp2low, col="red")
        
        # Create a title with a red, bold/italic font
        #title(main="Logit Two Proportions with 95% Confidence Interval", 
        #    col.main="red", font.main=4)
        
        plot(jitA[1,],jitA[2,], xlab="Proportion 1", ylab="Proportion 2")
        for (i in 1:n) {
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
