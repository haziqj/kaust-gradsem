invisible(suppressMessages(capture.output(source("R/kidney.R"))))
mod <- "
  GlyCon =~ y1 + y2 + y3
  KdnHlt =~ y4 + y5 + y6
  KdnHlt ~ GlyCon
"
code <- 'library(blavaan)
fit_blav <- bsem(mod, Kidney, meanstructure = TRUE, std.lv = TRUE)
print(fit_blav)'
source(textConnection(code), echo = TRUE, keep.source = TRUE, spaced = FALSE,
       prompt.echo = "> ", continue.echo = "+ ", max.deparse.length = Inf)
cat("> ")
