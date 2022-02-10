# Plot up the c_curve and lc_extrap files as sequenced reads against
# distinct reads.
# The script will rely on all of the data files being in the
# same directory as the script.

library(stringr)
library(tidyr)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
setwd(args[1])

# Need to read in multiple files to multiple dfs
# Then need to rename the expected_distinct column in each to contain
c_curve.file.list = list.files(path =".", pattern="*.c_curve.txt", full.names=T)
lc_extrap.file.list = list.files(path =".", pattern="*.lc_extrap.txt", full.names=T)
c_curve.file.list.short = list.files(path =".", pattern="*.c_curve.txt", full.names=F)
lc_extrap.file.list.short = list.files(path =".", pattern="*.lc_extrap.txt", full.names=F)

# Get the sample names and use these as the column names
c_curve.sample.names = lapply(X=c_curve.file.list.short, FUN=function(x){str_replace(x, "\\.c_curve\\.txt", "")})
lc_extrap.sample.names = lapply(X=lc_extrap.file.list.short, FUN=function(x){str_replace(x, "\\.lc_extrap\\.txt", "")})

# Read in the datafiles as a list of dfs
c_curve.df.list = lapply(c_curve.file.list, FUN=read.table, header=TRUE)
lc_extrap.df.list = lapply(lc_extrap.file.list, FUN=read.table, header=TRUE)

# Name the dfs according to the samples they represent
# We will use these names to label the columns of the concatenated df
names(c_curve.df.list) = c_curve.sample.names
names(lc_extrap.df.list) = lc_extrap.sample.names

c_curve.new.df.list = lapply(1:length(c_curve.df.list), FUN=function(i) {
  df = c_curve.df.list[[i]]
  sample.name = names(c_curve.df.list)[[i]]
  names(df)[names(df) == "distinct_reads"] = paste0(sample.name, "_distinct_reads")
  return(df)
})

lc_extrap.new.df.list = lapply(1:length(lc_extrap.df.list), FUN=function(i) {
  df = subset(lc_extrap.df.list[[i]], select= -c(LOWER_0.95CI, UPPER_0.95CI))
  sample.name = names(lc_extrap.df.list)[[i]]
  names(df)[names(df) == "EXPECTED_DISTINCT"] = paste0(sample.name, "_EXPECTED_DISTINCT")
  return(df)
})

# Get a list of the column names
c_curve.list.col.names = vapply(1:length(c_curve.new.df.list), FUN=function(x){names(c_curve.new.df.list[[x]])[[2]]}, FUN.VALUE="character")
lc_extrap.list.col.names = vapply(1:length(lc_extrap.new.df.list), FUN=function(x){names(lc_extrap.new.df.list[[x]])[[2]]}, FUN.VALUE="character")

c_curve.reduced.df = Reduce(f=function(x,y){merge(x, y, by="total_reads")}, x=c_curve.new.df.list)
lc_extrap.reduced.df = Reduce(f=function(x,y){merge(x, y, by="TOTAL_READS")}, x=lc_extrap.new.df.list)

# Combine into a single dataframe with the 

# Gather all of the different distinct_reads and EXPECTED_DISTINCT columns for 
# each df list into two key-value pairs so that they can be plotted as multiple
# lines on the plot.
c_curve.gathered.df = c_curve.reduced.df %>% pivot_longer(c_curve.list.col.names, names_to = "sample", values_to = "distinct_reads")
c_curve.gathered.df$prediction_type = "sampled"

lc_extrap.gathered.df = lc_extrap.reduced.df %>% pivot_longer(lc_extrap.list.col.names, names_to = "sample", values_to = "EXPECTED_DISTINCT")
lc_extrap.gathered.df$prediction_type = "extrapolated"

# Make the col names the same for merging
names(lc_extrap.gathered.df) = names(c_curve.gathered.df)
# merge
plotting_df = rbind(lc_extrap.gathered.df, c_curve.gathered.df)
# Finally order the sample factor so that the extrapolated samples come first in
# the plot
plotting_df$sample = factor(plotting_df$sample, levels=c(lc_extrap.list.col.names, c_curve.list.col.names))

# Finally plot the sample
complexity.plot = ggplot() + geom_line(data=plotting_df, aes(x=total_reads, y=distinct_reads, group=sample, color=prediction_type)) +
  xlim(0, 1.5e+09) + xlab("total_sequenced_reads")

ggsave("library.complexity.png", plot=complexity.plot)