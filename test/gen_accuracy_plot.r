x <- read.table("all_results.txt")

#rename columns
attr(x, "names") <- c("mtf50", "angle", "dev")

angles <- unique(sort(x$angle))
mtfs <- unique(sort(x$mtf50))

y <- list()
max_dev <- 0
min_dev <- 1
for (a in angles) {
    y[[a]] <- subset(x, angle == a)
    max_dev <- max(max_dev, y[[a]]$dev)
    min_dev <- min(min_dev, y[[a]]$dev)
}

pdf("accuracy_plot.pdf", width=(200/25.4), height=(120/25.4), pointsize=8)

plot(mtfs, aggregate(y[[angles[1]]]$dev, list(Sig=y[[angles[1]]]$mtf50), median)$x, 
  ylim=c(min_dev, max_dev), type="l", 
  xlab="MTF50 (c/p)", ylab="Error (c/p)",
  main="MTF50 error by slant angle"
)
for (a in 2:length(angles)) {
    lines(mtfs, aggregate(y[[angles[a]]]$dev, list(Sig=y[[angles[a]]]$mtf50), median)$x, ylim=c(min_dev, max_dev), col=a)
}

legend("topleft", legend=c("angle=4 deg", "angle=10 deg", "angle=30 deg"), fill=c(1,2,3), inset=0.05)
dev.off()
