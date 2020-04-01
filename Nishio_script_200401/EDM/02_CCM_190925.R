
# Create output folder
load("01_Embedding_Out/01_Embedding_Out_temp190925.RData")
dir.create("02_CCM_Out")

# Load library
library(rEDM); packageVersion("rEDM") # v0.6.9
library(ggplot2)
library(reshape2)

# Load function
source("functions/S01_HelperFunctions.R")

# CCM
# vs. RNA
K4I.RNA.res <- CCMwSurrogate("K4_I", "RNA", Ex = Ek4I, Ey = Erna, lag = -8:4)
CCMplot(K4I.RNA.res, "K4_I", "RNA")
CCMplotMagnif(K4I.RNA.res, "K4_I", "RNA")

K4II.RNA.res <- CCMwSurrogate("K4_II", "RNA", Ex = Ek4II, Ey = Erna, lag = -8:4)
CCMplot(K4II.RNA.res, "K4_II", "RNA")
CCMplotMagnif(K4II.RNA.res, "K4_II", "RNA")

K4VI.RNA.res <- CCMwSurrogate("K4_VI", "RNA", Ex = Ek4VI, Ey = Erna, lag = -8:4)
CCMplot(K4VI.RNA.res, "K4_VI", "RNA")
CCMplotMagnif(K4VI.RNA.res, "K4_VI", "RNA")

K4VII.RNA.res <- CCMwSurrogate("K4_VII", "RNA", Ex = Ek4VII, Ey = Erna, lag = -8:4)
CCMplot(K4VII.RNA.res, "K4_VII", "RNA")
CCMplotMagnif(K4VII.RNA.res, "K4_VII", "RNA")

K4VIII.RNA.res <- CCMwSurrogate("K4_VIII", "RNA", Ex = Ek4VIII, Ey = Erna, lag = -8:4)
CCMplot(K4VIII.RNA.res, "K4_VIII", "RNA")
CCMplotMagnif(K4VIII.RNA.res, "K4_VIII", "RNA")



K27I.RNA.res.CCM <- CCMwSurrogate("K27_I", "RNA", Ex = Ek27I, Ey = Erna, lag = -8:4)
CCMplot(K27I.RNA.res.CCM, "K27_I", "RNA")
CCMplotMagnif(K27I.RNA.res.CCM, "K27_I", "RNA")

K27II.RNA.res <- CCMwSurrogate("K27_II", "RNA", Ex = Ek27II, Ey = Erna, lag = -8:4)
CCMplot(K27II.RNA.res, "K27_II", "RNA")
CCMplotMagnif(K27II.RNA.res, "K27_II", "RNA")

K27III.RNA.res <- CCMwSurrogate("K27_III", "RNA", Ex = Ek27III, Ey = Erna, lag = -8:4)
CCMplot(K27III.RNA.res, "K27_III", "RNA")
CCMplotMagnif(K27III.RNA.res, "K27_III", "RNA")

K27IV.RNA.res <- CCMwSurrogate("K27_IV", "RNA", Ex = Ek27IV, Ey = Erna, lag = -8:4)
CCMplot(K27IV.RNA.res, "K27_IV", "RNA")
CCMplotMagnif(K27IV.RNA.res, "K27_IV", "RNA")

K27V.RNA.res <- CCMwSurrogate("K27_V", "RNA", Ex = Ek27V, Ey = Erna, lag = -8:4)
CCMplot(K27V.RNA.res, "K27_V", "RNA")
CCMplotMagnif(K27V.RNA.res, "K27_V", "RNA")

K27VI.RNA.res <- CCMwSurrogate("K27_VI", "RNA", Ex = Ek27VI, Ey = Erna, lag = -8:4)
CCMplot(K27VI.RNA.res, "K27_VI", "RNA")
CCMplotMagnif(K27VI.RNA.res, "K27_VI", "RNA")

K27VII.RNA.res <- CCMwSurrogate("K27_VII", "RNA", Ex = Ek27VII, Ey = Erna, lag = -8:4)
CCMplot(K27VII.RNA.res, "K27_VII", "RNA")
CCMplotMagnif(K27VII.RNA.res, "K27_VII", "RNA")

K27VIII.RNA.res <- CCMwSurrogate("K27_VIII", "RNA", Ex = Ek27VIII, Ey = Erna, lag = -8:4)
CCMplot(K27VIII.RNA.res, "K27_VIII", "RNA")
CCMplotMagnif(K27VIII.RNA.res, "K27_VIII", "RNA")



# K27 vs. K4
K27I.K4I.res.CCM <- CCMwSurrogate("K27_I", "K4_I", Ex = Ek27I, Ey = Ek4I, lag = -8:4)
CCMplot(K27I.K4I.res.CCM, "K27_I", "K4_I")
CCMplotMagnif(K27I.K4I.res.CCM, "K27_I", "K4_I")

K27II.K4II.res.CCM <- CCMwSurrogate("K27_II", "K4_II", Ex = Ek27II, Ey = Ek4II, lag = -8:4)
CCMplot(K27II.K4II.res.CCM, "K27_II", "K4_II")
CCMplotMagnif(K27II.K4II.res.CCM, "K27_II", "K4_II")

K27VI.K4VI.res.CCM <- CCMwSurrogate("K27_VI", "K4_VI", Ex = Ek27VI, Ey = Ek4VI, lag = -8:4)
CCMplot(K27VI.K4VI.res.CCM, "K27_VI", "K4_VI")
CCMplotMagnif(K27VI.K4VI.res.CCM, "K27_VI", "K4_VI")

K27VII.K4VII.res.CCM <- CCMwSurrogate("K27_VII", "K4_VII", Ex = Ek27VII, Ey = Ek4VII, lag = -12:4)
CCMplot(K27VII.K4VII.res.CCM, "K27_VII", "K4_VII")
CCMplotMagnif(K27VII.K4VII.res.CCM, "K27_VII", "K4_VII")

K27VIII.K4VIII.res.CCM <- CCMwSurrogate("K27_VIII", "K4_VIII", Ex = Ek27VIII, Ey = Ek4VIII, lag = -12:4)
CCMplot(K27VIII.K4VIII.res.CCM, "K27_VIII", "K4_VIII")
CCMplotMagnif(K27VIII.K4VIII.res.CCM, "K27_VIII", "K4_VIII")



# K27 vs. K27
K27I.K27II.res <- CCMwSurrogate("K27_I", "K27_II", Ex = Ek27I, Ey = Ek27II, lag = -8:4)
CCMplot(K27I.K27II.res, "K27_I", "K27_II")
CCMplotMagnif(K27I.K27II.res, "K27_I", "K27_II")

K27II.K27III.res <- CCMwSurrogate("K27_II", "K27_III", Ex = Ek27II, Ey = Ek27III, lag = -8:4)
CCMplot(K27II.K27III.res, "K27_II", "K27_III")
CCMplotMagnif(K27II.K27III.res, "K27_II", "K27_III")

K27III.K27IV.res <- CCMwSurrogate("K27_III", "K27_IV", Ex = Ek27III, Ey = Ek27IV, lag = -8:4)
CCMplot(K27III.K27IV.res, "K27_III", "K27_IV")
CCMplotMagnif(K27III.K27IV.res, "K27_III", "K27_IV")

K27IV.K27V.res <- CCMwSurrogate("K27_IV", "K27_V", Ex = Ek27IV, Ey = Ek27V, lag = -8:4)
CCMplot(K27IV.K27V.res, "K27_IV", "K27_V")
CCMplotMagnif(K27IV.K27V.res, "K27_IV", "K27_V")

K27V.K27VI.res <- CCMwSurrogate("K27_V", "K27_VI", Ex = Ek27V, Ey = Ek27VI, lag = -8:4)
CCMplot(K27V.K27VI.res, "K27_V", "K27_VI")
CCMplotMagnif(K27V.K27VI.res, "K27_V", "K27_VI")

K27VI.K27VII.res <- CCMwSurrogate("K27_VI", "K27_VII", Ex = Ek27VI, Ey = Ek27VII, lag = -8:4)
CCMplot(K27VI.K27VII.res, "K27_VI", "K27_VII")
CCMplotMagnif(K27VI.K27VII.res, "K27_VI", "K27_VII")

K27VII.K27VIII.res <- CCMwSurrogate("K27_VII", "K27_VIII", Ex = Ek27VII, Ey = Ek27VIII, lag = -8:4)
CCMplot(K27VII.K27VIII.res, "K27_VII", "K27_VIII")
CCMplotMagnif(K27VII.K27VIII.res, "K27_VII", "K27_VIII")



# K27_III-V vs. K27_I
K27I.K27III.res <- CCMwSurrogate("K27_I", "K27_III", Ex = Ek27I, Ey = Ek27III, lag = -8:4)
CCMplot(K27I.K27III.res, "K27_I", "K27_III")
CCMplotMagnif(K27I.K27III.res, "K27_I", "K27_III")

K27I.K27IV.res <- CCMwSurrogate("K27_I", "K27_IV", Ex = Ek27I, Ey = Ek27IV, lag = -8:4)
CCMplot(K27I.K27IV.res, "K27_I", "K27_IV")
CCMplotMagnif(K27I.K27IV.res, "K27_I", "K27_IV")

K27I.K27V.res <- CCMwSurrogate("K27_I", "K27_V", Ex = Ek27I, Ey = Ek27V, lag = -8:4)
CCMplot(K27I.K27V.res, "K27_I", "K27_V")
CCMplotMagnif(K27I.K27V.res, "K27_I", "K27_V")



# K4 vs. K4
K4I.K4II.res.CCM <- CCMwSurrogate("K4_I", "K4_II", Ex = Ek4I, Ey = Ek4II, lag = -8:4)
CCMplot(K4I.K4II.res.CCM, "K4_I", "K4_II")
CCMplotMagnif(K4I.K4II.res.CCM, "K4_I", "K4_II")

K4VI.K4VII.res.CCM <- CCMwSurrogate("K4_VI", "K4_VII", Ex = Ek4VI, Ey = Ek4VII, lag = -8:4)
CCMplot(K4VI.K4VII.res.CCM, "K4_VI", "K4_VII")
CCMplotMagnif(K4VI.K4VII.res.CCM, "K4_VI", "K4_VII")

K4VII.K4VIII.res.CCM <- CCMwSurrogate("K4_VII", "K4_VIII", Ex = Ek4VII, Ey = Ek4VIII, lag = -8:4)
CCMplot(K4VII.K4VIII.res.CCM, "K4_VII", "K4_VIII")
CCMplotMagnif(K4VII.K4VIII.res.CCM, "K4_VII", "K4_VIII")



# K4_I vs. K4_VI-VIII
K4I.K4VI.res.CCM <- CCMwSurrogate("K4_I", "K4_VI", Ex = Ek4I, Ey = Ek4VI, lag = -8:4)
CCMplot(K4I.K4VI.res.CCM, "K4_I", "K4_VI")
CCMplotMagnif(K4I.K4VI.res.CCM, "K4_I", "K4_VI")

K4I.K4VII.res.CCM <- CCMwSurrogate("K4_I", "K4_VII", Ex = Ek4I, Ey = Ek4VII, lag = -8:4)
CCMplot(K4I.K4VII.res.CCM, "K4_I", "K4_VII")
CCMplotMagnif(K4I.K4VII.res.CCM, "K4_I", "K4_VII")

K4I.K4VIII.res.CCM <- CCMwSurrogate("K4_I", "K4_VIII", Ex = Ek4I, Ey = Ek4VIII, lag = -8:4)
CCMplot(K4I.K4VIII.res.CCM, "K4_I", "K4_VIII")
CCMplotMagnif(K4I.K4VIII.res.CCM, "K4_I", "K4_VIII")


# K27_I vs. K4_VI-VIII
K27I.K4VI.res.CCM <- CCMwSurrogate("K27_I", "K4_VI", Ex = Ek27I, Ey = Ek4VI, lag = -8:4)
CCMplot(K27I.K4VI.res.CCM, "K27_I", "K4_VI")
CCMplotMagnif(K27I.K4VI.res.CCM, "K27_I", "K4_VI")

K27I.K4VII.res.CCM <- CCMwSurrogate("K27_I", "K4_VII", Ex = Ek27I, Ey = Ek4VII, lag = -8:4)
CCMplot(K27I.K4VII.res.CCM, "K27_I", "K4_VII")
CCMplotMagnif(K27I.K4VII.res.CCM, "K27_I", "K4_VII")

K27I.K4VIII.res.CCM <- CCMwSurrogate("K27_I", "K4_VIII", Ex = Ek27I, Ey = Ek4VIII, lag = -8:4)
CCMplot(K27I.K4VIII.res.CCM, "K27_I", "K4_VIII")
CCMplotMagnif(K27I.K4VIII.res.CCM, "K27_I", "K4_VIII")




# vs. Temp
K4I.Temp.res.CCM <- CCMwSurrogateTemp("K4_I", "Temp", Ex = Ek4I, lag = -8:4)
CCMplotTemp(K4I.Temp.res.CCM, "K4_I", "Temp")
CCMplotMagnifTemp(K4I.Temp.res.CCM, "K4_I", "Temp")

K4II.Temp.res.CCM <- CCMwSurrogateTemp("K4_II", "Temp", Ex = Ek4II, lag = -8:4)
CCMplotTemp(K4II.Temp.res.CCM, "K4_II", "Temp")
CCMplotMagnifTemp(K4II.Temp.res.CCM, "K4_II", "Temp")

K4VI.Temp.res.CCM <- CCMwSurrogateTemp("K4_VI", "Temp", Ex = Ek4VI, lag = -8:4)
CCMplotTemp(K4VI.Temp.res.CCM, "K4_VI", "Temp")
CCMplotMagnifTemp(K4VI.Temp.res.CCM, "K4_VI", "Temp")

K4VII.Temp.res.CCM <- CCMwSurrogateTemp("K4_VII", "Temp", Ex = Ek4VII, lag = -8:4)
CCMplotTemp(K4VII.Temp.res.CCM, "K4_VII", "Temp")
CCMplotMagnifTemp(K4VII.Temp.res.CCM, "K4_VII", "Temp")

K4VIII.Temp.res.CCM <- CCMwSurrogateTemp("K4_VIII", "Temp", Ex = Ek4VIII, lag = -8:4)
CCMplotTemp(K4VIII.Temp.res.CCM, "K4_VIII", "Temp")
CCMplotMagnifTemp(K4VIII.Temp.res.CCM, "K4_VIII", "Temp")

K27I.Temp.res.CCM <- CCMwSurrogateTemp("K27_I", "Temp", Ex = Ek27I, lag = -8:4)
CCMplotTemp(K27I.Temp.res.CCM, "K27_I", "Temp")
CCMplotMagnifTemp(K27I.Temp.res.CCM, "K27_I", "Temp")

K27II.Temp.res.CCM <- CCMwSurrogateTemp("K27_II", "Temp", Ex = Ek27II, lag = -8:4)
CCMplotTemp(K27II.Temp.res.CCM, "K27_II", "Temp")
CCMplotMagnifTemp(K27II.Temp.res.CCM, "K27_II", "Temp")

K27III.Temp.res.CCM <- CCMwSurrogateTemp("K27_III", "Temp", Ex = Ek27III, lag = -8:4)
CCMplotTemp(K27III.Temp.res.CCM, "K27_III", "Temp")
CCMplotMagnifTemp(K27III.Temp.res.CCM, "K27_III", "Temp")

K27IV.Temp.res.CCM <- CCMwSurrogateTemp("K27_IV", "Temp", Ex = Ek27IV, lag = -8:4)
CCMplotTemp(K27IV.Temp.res.CCM, "K27_IV", "Temp")
CCMplotMagnifTemp(K27IV.Temp.res.CCM, "K27_IV", "Temp")

K27V.Temp.res.CCM <- CCMwSurrogateTemp("K27_V", "Temp", Ex = Ek27V, lag = -8:4)
CCMplotTemp(K27V.Temp.res.CCM, "K27_V", "Temp")
CCMplotMagnifTemp(K27V.Temp.res.CCM, "K27_V", "Temp")

K27VI.Temp.res.CCM <- CCMwSurrogateTemp("K27_VI", "Temp", Ex = Ek27VI, lag = -8:4)
CCMplotTemp(K27VI.Temp.res.CCM, "K27_VI", "Temp")
CCMplotMagnifTemp(K27VI.Temp.res.CCM, "K27_VI", "Temp")

K27VII.Temp.res.CCM <- CCMwSurrogateTemp("K27_VII", "Temp", Ex = Ek27VII, lag = -8:4)
CCMplotTemp(K27VII.Temp.res.CCM, "K27_VII", "Temp")
CCMplotMagnifTemp(K27VII.Temp.res.CCM, "K27_VII", "Temp")

K27VIII.Temp.res.CCM <- CCMwSurrogateTemp("K27_VIII", "Temp", Ex = Ek27VIII, lag = -8:4)
CCMplotTemp(K27VIII.Temp.res.CCM, "K27_VIII", "Temp")
CCMplotMagnifTemp(K27VIII.Temp.res.CCM, "K27_VIII", "Temp")

RNA.Temp.res.CCM <- CCMwSurrogateTemp("RNA", "Temp", Ex = Erna, lag = -8:4)
CCMplotTemp(RNA.Temp.res.CCM, "RNA", "Temp")
CCMplotMagnifTemp(RNA.Temp.res.CCM, "RNA", "Temp")



# Save output
save.image("02_CCM_Out/02_CCM_Out_temp190925.RData")
