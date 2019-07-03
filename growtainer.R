
library(tidyverse)
newdata = tibble(channel = dir("growtainerspectra/")) %>%
      mutate(dat = map(paste0("growtainerspectra/", channel), read.delim,
                            col.names=c("WL", "dark", "scope", "norm", "Relative"), header=F),
           fit=map(dat, function (df) {
              newdata = data.frame(WL=400:750)
              rel = predict(loess(Relative~WL, data=df, span=0.05), newdata=newdata)
              newdata$relintensity = pmax(rel/sum(rel), 0) * sum(df$norm) # Clamp and re-normalise
              newdata
           })) %>%
      unnest(fit) %>%
      mutate(relintensity=round(relintensity/max(relintensity), 4))
wide = newdata %>%
  spread(channel, relintensity)
write_tsv(wide, "data/growtainerspectra.tsv")


# Uncomment for plot
ggplot(newdata, aes(WL, relintensity, colour=channel)) +
  geom_line() + theme_bw()
