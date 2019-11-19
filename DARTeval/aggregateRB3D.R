library(daRt)
library(reshape2)
library(dplyr)
library(ggplot2)
library(raster)
library(doParallel)
library(QOLfunctions)
library(DBI)
library(gganimate)

plotDir <- "C:/Users/micromet/Desktop/"
dartSimulationDir <- "V:/Tier_processing/hv867657/DART_5-7-5_1126/user_data/simulations"
#the DEM that is exactly centered on the DART simulation domain.
DEMfileName <- "D:/gitProjects/GoogleEarthModelCreator/simulations/city2km/out/DEM_merged.tif"
dartBaseSimName <- "SPARTACUS_SW"
sequenceID <- "SW_perType"
###simulations to read
iterVals <- c("ILLUDIR", "ITER5")
variablesRB3DVals <- c("Intercepted", "Scattered",
                       "Absorbed", "+ZFaceExit", "+ZFaceEntry")
#how many cells to keep that are "underground" i.e. how much vertical variation underground
#(e.g. exposed foundations) should be resolved. relax this parameter at expense of array sizes
#and memory footprint.
maxUndergroundCellsVal <- 5
#max memory to use for raster (before cache to HDD)
rasterOptions(maxmemory = 1.6e+10)
#DART parameter: height above ground that was set in DART for the 3D model(s)
DARTmodelElevationParam <- 5
#the "typeNums" in DART - and what they should refer to in real terms e.g. 103_TypeNum = "Wall".
typeNumsAll <- c(     "",     paste0(c("103", "104",  "105",  "106",  "107",  "108",    "109"), "_TypeNum"))
typeNumsAll_labs <- c("Unclassified", "Wall", "Wall", "Wall", "Roof", "Wall", "Ground", "Wall")

#convert raw dart type numbers to a labelled factor
typeNumConvert <- function(x) {
  factor(x = x,
         levels = typeNumsAll,
         labels = typeNumsAll_labs)
}

dartSimBaseDir <- file.path(dartSimulationDir, dartBaseSimName, "daRtinput", sequenceID)
setwd(dartSimBaseDir)
con <- dbConnect(RSQLite::SQLite(), dbname = paste0(dartBaseSimName, "_JTT.db"))
dbDF <- dbReadTable(con, dartBaseSimName)
simDir <- gsub(" ", "", dbDF$scriptDirectory[dbDF$exitStatus == 0])
#hack if usign windows mounted drive
simDir <- gsub("/storage/basic/micromet/", "V:/", simDir)
message(paste("Numbers of unfinished sims:", length(which(dbDF$exitStatus != 0))))

DEM <- raster::raster(x = DEMfileName)
zStatsList <- list()

sF_all <- daRt::simulationFilter(
  variables = "RADIATIVE_BUDGET",
  bands = 0L,
  product = "rb3D", iters = iterVals,
  variablesRB3D = variablesRB3DVals,
  typeNums = typeNumsAll)


allFiles <- daRt::getFiles(simDir, sF = sF_all)
daRt::deleteFiles(x = allFiles, trianglesInput = TRUE, maketOutput = TRUE, #trianglesOutput = TRUE,
                  deleteSimulationFiles = FALSE)
seqParams <- daRt::sequenceParameters(allFiles)

iterNo <- 1 #preallocate a tracker for the iteration(s)
for (i in 1:length(simDir)) {
  print(paste(iterNo, "/", length(simDir)))
  sF <- daRt::simulationFilter(
    variables = "RADIATIVE_BUDGET",
    product = "rb3D", iters = iterVals,
    variablesRB3D = variablesRB3DVals,
    bands = 0L,
    typeNums = typeNumsAll)
  d <- daRt::getData(x = simDir[i], sF = sF, nCores = 4)
  d1 <- daRt::rb3DtoNc(x = d)
  rm(d); gc()
  d1 <- daRt::removeRelief(x = d1, DEM = DEM, DARTmodelElevation = DARTmodelElevationParam,
                           maxUndergroundCells = maxUndergroundCellsVal)
  zStatsList[[iterNo]] <- as.data.frame(d1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(typeNum_raw = typeNum) %>%
    dplyr::mutate(typeNum = typeNumConvert(typeNum)) %>%
    dplyr::group_by(Z, variableRB3D, band, iter, typeNum, simName) %>%
    dplyr::summarise(nTypeNos = length(unique(typeNum_raw)),
                     #as the sample sizes are the same, the below equates to the mean of the sum
                     #of any grouped typeNos
                     meanVal = mean(value) * nTypeNos)
  iterNo <- iterNo + 1
}

zStats <- dplyr::bind_rows(zStatsList)

zStats_params <- zStats %>%
  dplyr::left_join(seqParams)

zStats_params_groupedAzimuth <- zStats_params %>%
  dplyr::group_by(Z, variableRB3D, band, iter, typeNum, parametre4, parametre11) %>%
  dplyr::summarise(mean_meanVal = mean(meanVal),
                   p5_meanVal = quantile(meanVal, 0.05),
                   p95_meanVal = quantile(meanVal, 0.95),
                   min_meanVal = min(meanVal),
                   max_meanVal = max(meanVal),
                   n_meanVal = length(meanVal))

pltOut <- ggplot(zStats_params_groupedAzimuth) +
  geom_ribbon(aes(x = Z, ymin = min_meanVal, ymax = max_meanVal,
                  fill = typeNum,
                  colour = typeNum),
              lwd = 0.4, alpha = 0.8) +
  facet_grid(parametre4 + iter ~ parametre11 + variableRB3D, scales = "free_x") +
  coord_flip() +
  xlab("Z (m)") +
  theme_gray() +
  labs(fill = "typeNum") +
  guides(colour = FALSE) +
  scale_fill_brewer(palette = 3, type = "qual") +
  scale_color_brewer(palette = 3, type = "qual")

plotName <- file.path(plotDir, paste0(sequenceID, jYHHMMSS(), ".pdf"))
ggsave(filename = plotName, plot = pltOut,
       height = 24, width = 24)

toSave <- zStats_params_groupedAzimuth %>%
  dplyr::rename(SZA = parametre4, albedo = parametre11)
textName <- file.path(plotDir, paste0(sequenceID, jYHHMMSS(), ".txt"))

write.table(toSave, textName, sep = ",", row.names = FALSE)
