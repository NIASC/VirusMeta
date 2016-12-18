###################################################################################
#OR calculation function
oddsratioWald.proc <- function(n00, n01, n10, n11, alpha = 0.05){
    #
    #  Compute the odds ratio between two binary variables:

    #M2 <- matrix(c(1438, 1427, 18496, 18515), nrow = 2)
    #rownames(M2) <- c("Exposure Yes", "Exposure No")
    #colnames(M2) <- c("Disease Yes", "Disease No")
    #M2
    #      Disease Yes Disease No
    #Exposure Yes       n00     n01
    #Exposure No        n10     n11

    #oddsratioWald.proc(M2[1,1],M2[1,2],M2[2,1],M2[2,2])

    OR <- (n00 * n11)/(n01 * n10)
    #
    #  Compute the Wald confidence intervals:
    #
    siglog <- sqrt((1/n00) + (1/n01) + (1/n10) + (1/n11))
    zalph <- qnorm(1 - alpha/2)
    logOR <- log(OR)
    loglo <- logOR - zalph * siglog
    loghi <- logOR + zalph * siglog
    #
    ORlo <- exp(loglo)
    ORhi <- exp(loghi)
    #
    oframe <- data.frame(LowerCI = ORlo, OR = OR, UpperCI = ORhi, alpha = alpha)
    oframe
}
####################################################################################

