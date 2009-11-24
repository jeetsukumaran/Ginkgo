library(ggplot2)
abcdplot.points.ppdivtime = function(gentimes) {
    p = ggplot(gentimes)
    p = p + geom_line(aes(x=Gen, y=Hours, group=PopSize))
    pvalues = unique(gentimes$PopSize)
    plabels = c()
    x = c()
    y = c()
    for (pvalue in pvalues) {
        plabels = append(plabels, paste("N = ", pvalue))
        x = append(x, max(subset(gentimes, PopSize==pvalue, Gen)))
        y = append(y, max(subset(gentimes, PopSize==pvalue, Hours)))
    }

    p = p + annotate("text", x=x-5, y=y+0.5, label=plabels)
    return(p)
}
