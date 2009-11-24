library(ggplot2)

plot.gentimes.x.popsize = function(gentimes) {
    p = ggplot(gentimes) + geom_line(aes(x=Gen, y=Hours, group=PopSize))
    pvalues = unique(gentimes$PopSize)
    plabels = c()
    x = c()
    y = c()
    for (pvalue in pvalues) {
        plabels = append(plabels, paste("N = ", pvalue))
        x = append(x, max(subset(gentimes, PopSize==pvalue, Gen)))
        y = append(y, max(subset(gentimes, PopSize==pvalue, Hours)))
    }

    p = p + annotate("text", x=x, y=y, label=plabels, size=3, hjust=0.5, vjust=-0.5)

    p = p + xlab("Number of Generations") + ylab("Number of Hours")

    return(p)
}

plot.genpopsize = function(gentimes, gen=1000) {
    p = ggplot(subset(gentimes, Gen==gen)) + geom_line(aes(x=PopSize, y=Hours))
    p = p + xlab("Total Number of Organisms") + ylab("Number of Hours to Complete 1000 Generations")
    return(p)
}


plot.memusage.x.popsize = function(memusage) {
    p = ggplot(memusage) + geom_line(aes(x=Hour, y=Memory, group=PopSize))
    p = p + xlab("Number of Hours") + ylab("Memory Usage")
    pvalues = unique(memusage$PopSize)
    plabels = c()
    x = c()
    y = c()
    for (pvalue in pvalues) {
        plabels = append(plabels, paste("N = ", pvalue))
        x = append(x, max(subset(memusage, PopSize==pvalue, Hour)))
        y = append(y, max(subset(memusage, PopSize==pvalue, Memory)))
    }

    p = p + annotate("text", x=x, y=y, label=plabels, size=3, hjust=0.5, vjust=-0.5)
    return(p)
}
