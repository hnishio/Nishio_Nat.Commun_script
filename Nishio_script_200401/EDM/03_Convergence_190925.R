
# Create output folder
load("01_Embedding_Out/01_Embedding_Out.RData")
dir.create("05_Convergence_Out")

# Load library
library(rEDM); packageVersion("rEDM") # v0.6.9, 2018.4.22
library(ggplot2)
library(reshape2)

# Load function
source("functions/S01_HelperFunctions.R")

# Convergence

### vs. RNA
# K27_I vs RNA
K27I.RNA.res <- Convergence(d, x="K27_I", y="RNA", Ex=Ek27I, Ey=Erna, tpx=-1, tpy=-1)
Convergence_plot(K27I.RNA.res, x="K27_I", y="RNA", tpx=-1, tpy=-1)

# K27_II vs RNA
K27II.RNA.res <- Convergence(d, x="K27_II", y="RNA", Ex= Ek27II, Ey=Erna, tpx=-1, tpy=-1)
Convergence_plot(K27II.RNA.res, x="K27_II", y="RNA", tpx=-1, tpy=-1)

# K27_III vs RNA
K27III.RNA.res <- Convergence(d, x="K27_III", y="RNA", Ex=Ek27III, Ey=Erna, tpx=-3, tpy=-3)
Convergence_plot(K27III.RNA.res, x="K27_III", y="RNA", tpx=-3, tpy=-3)

# K27_IV vs RNA
K27IV.RNA.res <- Convergence(d, x="K27_IV", y="RNA", Ex=Ek27IV, Ey=Erna, tpx=-4, tpy=-4)
Convergence_plot(K27IV.RNA.res, x="K27_IV", y="RNA", tpx=-4, tpy=-4)

# K27_V vs RNA
K27V.RNA.res <- Convergence(d, x="K27_V", y="RNA", Ex=Ek27V, Ey=Erna, tpx=-4, tpy=-4)
Convergence_plot(K27V.RNA.res, x="K27_V", y="RNA", tpx=-4, tpy=-4)

# K27_VI vs RNA
K27VI.RNA.res <- Convergence(d, x="K27_VI", y="RNA", Ex=Ek27VI, Ey=Erna, tpx=-4, tpy=-4)
Convergence_plot(K27VI.RNA.res, x="K27_VI", y="RNA", tpx=-4, tpy=-4)

# K27_VII vs RNA
K27VII.RNA.res <- Convergence(d, x="K27_VII", y="RNA", Ex=Ek27VII, Ey=Erna, tpx=-4, tpy=-4)
Convergence_plot(K27VII.RNA.res, x="K27_VII", y="RNA", tpx=-4, tpy=-4)

# K27_VIII vs RNA
K27VIII.RNA.res <- Convergence(d, x="K27_VIII", y="RNA", Ex=Ek27VIII, Ey=Erna, tpx=-6, tpy=-6)
Convergence_plot(K27VIII.RNA.res, x="K27_VIII", y="RNA", tpx=-6, tpy=-6)

# K4_I vs RNA
K4I.RNA.res <- Convergence(d, x="K4_I", y="RNA", Ex=Ek4I, Ey=Erna, tpx=0, tpy=-1)
Convergence_plot(K4I.RNA.res, x="K4_I", y="RNA", tpx=0, tpy=-1)



### K27 vs. K4
# K27_I vs K4_I
K27I.K4I.res <- Convergence(d, x="K27_I", y="K4_I", Ex=Ek27I, Ey=Ek4I, tpx=-1, tpy=-1)
Convergence_plot(K27I.K4I.res, x="K27_I", y="K4_I", tpx=-1, tpy=-1)

# K27_II vs K4_II
K27II.K4II.res <- Convergence(d, x="K27_II", y="K4_II", Ex=Ek27II, Ey=Ek4II, tpx=-1, tpy=-1)
Convergence_plot(K27II.K4II.res, x="K27_II", y="K4_II", tpx=-1, tpy=-1)

# K27_VI vs K4_VI
K27VI.K4VI.res <- Convergence(d, x="K27_VI", y="K4_VI", Ex=Ek27VI, Ey=Ek4VI, tpx=-4, tpy=-4)
Convergence_plot(K27VI.K4VI.res, x="K27_VI", y="K4_VI", tpx=-4, tpy=-4)

# K27_VII vs K4_VII
K27VII.K4VII.res <- Convergence(d, x="K27_VII", y="K4_VII", Ex=Ek27VII, Ey=Ek4VII, tpx=0, tpy=0)
Convergence_plot(K27VII.K4VII.res, x="K27_VII", y="K4_VII", tpx=0, tpy=0)

# K27_VIII vs K4_VIII
K27VIII.K4VIII.res <- Convergence(d, x="K27_VIII", y="K4_VIII", Ex=Ek27VIII, Ey=Ek4VIII, tpx=-9, tpy=-9)
Convergence_plot(K27VIII.K4VIII.res, x="K27_VIII", y="K4_VIII", tpx=-9, tpy=-9)



### K27 vs. K27
# K27_I vs K27_II
K27I.K27II.res <- Convergence(d, x="K27_I", y="K27_II", Ex=Ek27I, Ey=Ek27II, tpx=0, tpy=-1)
Convergence_plot(K27I.K27II.res, x="K27_I", y="K27_II", tpx=0, tpy=-1)

# K27_I vs K27_III
K27I.K27III.res <- Convergence(d, x="K27_I", y="K27_III", Ex=Ek27I, Ey=Ek27III, tpx=1, tpy=-2)
Convergence_plot(K27I.K27III.res, x="K27_I", y="K27_III", tpx=1, tpy=-2)

# K27_I vs K27_IV
K27I.K27IV.res <- Convergence(d, x="K27_I", y="K27_IV", Ex=Ek27I, Ey=Ek27IV, tpx=1, tpy=-3)
Convergence_plot(K27I.K27IV.res, x="K27_I", y="K27_IV", tpx=1, tpy=-3)

# K27_I vs K27_V
K27I.K27V.res <- Convergence(d, x="K27_I", y="K27_V", Ex=Ek27I, Ey=Ek27V, tpx=1, tpy=-3)
Convergence_plot(K27I.K27V.res, x="K27_I", y="K27_V", tpx=1, tpy=-3)

# K27_II vs K27_III
K27II.K27III.res <- Convergence(d, x="K27_II", y="K27_III", Ex=Ek27II, Ey=Ek27III, tpx=1, tpy=-2)
Convergence_plot(K27II.K27III.res, x="K27_II", y="K27_III", tpx=1, tpy=-2)

# K27_III vs K27_IV
K27III.K27IV.res <- Convergence(d, x="K27_III", y="K27_IV", Ex=Ek27III, Ey=Ek27IV, tpx=0, tpy=-1)
Convergence_plot(K27III.K27IV.res, x="K27_III", y="K27_IV", tpx=0, tpy=-1)

# K27_IV vs K27_V
K27IV.K27V.res <- Convergence(d, x="K27_IV", y="K27_V", Ex=Ek27IV, Ey=Ek27V, tpx=-1, tpy=-1)
Convergence_plot(K27IV.K27V.res, x="K27_IV", y="K27_V", tpx=-1, tpy=-1)

# K27_V vs K27_VI
K27V.K27VI.res <- Convergence(d, x="K27_V", y="K27_VI", Ex=Ek27V, Ey=Ek27VI, tpx=-1, tpy=-1)
Convergence_plot(K27V.K27VI.res, x="K27_V", y="K27_VI", tpx=-1, tpy=-1)

# K27_VI vs K27_VII
K27VI.K27VII.res <- Convergence(d, x="K27_VI", y="K27_VII", Ex=Ek27VI, Ey=Ek27VII, tpx=-1, tpy=-2)
Convergence_plot(K27VI.K27VII.res, x="K27_VI", y="K27_VII", tpx=-1, tpy=-2)

# K27_VII vs K27_VIII
K27VII.K27VIII.res <- Convergence(d, x="K27_VII", y="K27_VIII", Ex=Ek27VII, Ey=Ek27VIII, tpx=0, tpy=-4)
Convergence_plot(K27VII.K27VIII.res, x="K27_VII", y="K27_VIII", tpx=0, tpy=-4)



### K4 vs. K4
# K4_I vs K4_II
K4I.K4II.res <- Convergence(d, x="K4_I", y="K4_II", Ex=Ek4I, Ey=Ek4II, tpx=-2, tpy=-1)
Convergence_plot(K4I.K4II.res, x="K4_I", y="K4_II", tpx=-2, tpy=-1)

# K4_VI vs K4_VII
K4VI.K4VII.res <- Convergence(d, x="K4_VI", y="K4_VII", Ex=Ek4VI, Ey=Ek4VII, tpx=0, tpy=0)
Convergence_plot(K4VI.K4VII.res, x="K4_VI", y="K4_VII", tpx=0, tpy=0)

# K4_VII vs K4_VIII
K4VII.K4VIII.res <- Convergence(d, x="K4_VII", y="K4_VIII", Ex=Ek4VII, Ey=Ek4VIII, tpx=0, tpy=0)
Convergence_plot(K4VII.K4VIII.res, x="K4_VII", y="K4_VIII", tpx=0, tpy=0)



# K4_I vs K4_VI
K4I.K4VI.res <- Convergence(d, x="K4_I", y="K4_VI", Ex=Ek4I, Ey=Ek4VI, tpx=-3, tpy=0)
Convergence_plot(K4I.K4VI.res, x="K4_I", y="K4_VI", tpx=-3, tpy=0)

# K4_I vs K4_VII
K4I.K4VII.res <- Convergence(d, x="K4_I", y="K4_VII", Ex=Ek4I, Ey=Ek4VII, tpx=-3, tpy=0)
Convergence_plot(K4I.K4VII.res, x="K4_I", y="K4_VII", tpx=-3, tpy=0)

# K4_I vs K4_VIII
K4I.K4VIII.res <- Convergence(d, x="K4_I", y="K4_VIII", Ex=Ek4I, Ey=Ek4VIII, tpx=-3, tpy=-1)
Convergence_plot(K4I.K4VIII.res, x="K4_I", y="K4_VIII", tpx=-3, tpy=-1)




# K27_I vs K4_VII
K27I.K4VII.res <- Convergence(d, x="K27_I", y="K4_VII", Ex=Ek27I, Ey=Ek4VII, tpx=-3, tpy=0)
Convergence_plot(K27I.K4VII.res, x="K27_I", y="K4_VII", tpx=-3, tpy=0)

# K27_I vs K4_VIII
K27I.K4VIII.res <- Convergence(d, x="K27_I", y="K4_VIII", Ex=Ek27I, Ey=Ek4VIII, tpx=-5, tpy=1)
Convergence_plot(K27I.K4VIII.res, x="K27_I", y="K4_VIII", tpx=-5, tpy=1)








source("functions/S01_HelperFunctions.R")
### vs. temperature
K4I.Temp.res <- Convergence_Temp(d, x="K4_I", y="Temp", Ex=Ek4I, tpx=-2)
Convergence_plot_Temp(K4I.Temp.res, x="K4_I", y="Temp", tpx=-2)

K4II.Temp.res <- Convergence_Temp(d, x="K4_II", y="Temp", Ex=Ek4II, tpx=-3)
Convergence_plot_Temp(K4II.Temp.res, x="K4_II", y="Temp", tpx=-3)

K4VI.Temp.res <- Convergence_Temp(d, x="K4_VI", y="Temp", Ex=Ek4VI, tpx=0)
Convergence_plot_Temp(K4VI.Temp.res, x="K4_VI", y="Temp", tpx=0)

K4VII.Temp.res <- Convergence_Temp(d, x="K4_VII", y="Temp", Ex=Ek4VII, tpx=-1)
Convergence_plot_Temp(K4VII.Temp.res, x="K4_VII", y="Temp", tpx=-1)

K4VIII.Temp.res <- Convergence_Temp(d, x="K4_VIII", y="Temp", Ex=Ek4VIII, tpx=0)
Convergence_plot_Temp(K4VIII.Temp.res, x="K4_VIII", y="Temp", tpx=0)


K27I.Temp.res <- Convergence_Temp(d, x="K27_I", y="Temp", Ex=Ek27I, tpx=-2)
Convergence_plot_Temp(K27I.Temp.res, x="K27_I", y="Temp", tpx=-2)

K27II.Temp.res <- Convergence_Temp(d, x="K27_II", y="Temp", Ex=Ek27II, tpx=-3)
Convergence_plot_Temp(K27II.Temp.res, x="K27_II", y="Temp", tpx=-3)

K27III.Temp.res <- Convergence_Temp(d, x="K27_III", y="Temp", Ex=Ek27III, tpx=-4)
Convergence_plot_Temp(K27III.Temp.res, x="K27_III", y="Temp", tpx=-4)

K27IV.Temp.res <- Convergence_Temp(d, x="K27_IV", y="Temp", Ex=Ek27IV, tpx=-5)
Convergence_plot_Temp(K27IV.Temp.res, x="K27_IV", y="Temp", tpx=-5)

K27V.Temp.res <- Convergence_Temp(d, x="K27_V", y="Temp", Ex=Ek27V, tpx=-6)
Convergence_plot_Temp(K27V.Temp.res, x="K27_V", y="Temp", tpx=-6)


K27VI.Temp.res <- Convergence_Temp(d, x="K27_VI", y="Temp", Ex=Ek27VI, tpx=-6)
Convergence_plot_Temp(K27VI.Temp.res, x="K27_VI", y="Temp", tpx=-6)

K27VII.Temp.res <- Convergence_Temp(d, x="K27_VII", y="Temp", Ex=Ek27VII, tpx=-7)
Convergence_plot_Temp(K27VII.Temp.res, x="K27_VII", y="Temp", tpx=-7)

K27VIII.Temp.res <- Convergence_Temp(d, x="K27_VIII", y="Temp", Ex=Ek27VIII, tpx=-6)
Convergence_plot_Temp(K27VIII.Temp.res, x="K27_VIII", y="Temp", tpx=-6)

RNA.Temp.res <- Convergence_Temp(d, x="RNA", y="Temp", Ex=Erna, tpx=-3)
Convergence_plot_Temp(RNA.Temp.res, x="RNA", y="Temp", tpx=-3)



# Save output
save.image("03_Convergence_Out/03_Convergence_temp190925.RData")
