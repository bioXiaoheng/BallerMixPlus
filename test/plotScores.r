#plot the -logA and the score of each replicate
library(ggplot2)
#library(ggpubr)
args <- commandArgs(trailingOnly=TRUE)

filename <- args[1] 

figname <- args[2]


df = read.table(filename, sep='\t', header = TRUE )

p <- ggplot(df) + geom_line(aes(x = physPos, y = CLR, color = log10(s_hat)) ) + 
	xlab('Physical Position') + ylab('Composite Likelihood Ratio') +
	scale_color_gradient2(name = expression('log'[10]~hat(italic(a))), limits = c(-3, 9), midpoint = 0, high="red", mid = 'white', low="blue")

p <- p + theme_minimal()


png(figname, width=960, height=400)
print(p)
dev.off()