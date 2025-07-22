setwd("/Users/vincedanniballe/Library/Mobile Documents/com~apple~CloudDocs/Documents/02 - Labs/Goldstein Lab/AD/AD Manuscript/Cell_chat")
library(anndata)
library(Matrix)
library(CellChat)
library(patchwork)
library(future)

ad <- read_h5ad("250711_AD_19_samples_cell_chat_granular_final.h5ad")

counts <- Matrix::t( Matrix::Matrix(ad$X, sparse = TRUE) )        # genes × cells
rownames(counts) <- ad$var_names                                  # gene symbols
colnames(counts) <- rownames(ad$obs)                              # barcodes

lib.size   <- Matrix::colSums(counts)
data.input <- as( log1p( Matrix::t( Matrix::t(counts) / lib.size ) * 1e4 ),
                  "dgCMatrix" )

meta               <- ad$obs
meta$samples       <- meta$orig_patients          # sample IDs CellChat needs
meta$labels        <- meta$TCAT_all # cell-type labels
status_levels      <- unique(meta$Alz_status)

## ── 2 · run the FULL CellChat pipeline for each status ───────────────────
CellChatDB.use <- subsetDB(CellChatDB.human)       # drops “Non-protein Signaling”
plan(multisession, workers = 4)                    # parallel for heavy steps

object.list <- lapply(status_levels, function(st) {
  message("→  processing Alz_status: ", st)
  sel <- meta$Alz_status == st
  
  # ----- create CellChat object ------------------------------------------
  cc  <- createCellChat(object = data.input[, sel],
                        meta   = meta[sel, ],
                        group.by = "labels")
  setIdent(cc, ident.use = "labels")
  cc@DB <- CellChatDB.use
  
  # ----- per-tutorial preprocessing & inference --------------------------
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, population.size = FALSE)    # triMean default
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc)
  
  return(cc)
})

warnings()

names(object.list) <- status_levels
saveRDS(object.list, "cellchat_object.list_three_statuses_done.rds")
object.list <- readRDS("cellchat_object.list_three_statuses_done.rds")



## ── interactions and interaction strength  -------------------------
object.list <- object.list[c("Control", "Pre-Clinical", "Clinical")]
cellchat    <- mergeCellChat(object.list, add.names = names(object.list),
                             cell.prefix = TRUE)

gg1 <- compareInteractions(cellchat, group = 1:3, show.legend = FALSE)         # #links
gg2 <- compareInteractions(cellchat, group = 1:3, measure = "weight",
                           show.legend = FALSE)                                # strength
gg1 + gg2


pair_CP <- mergeCellChat(
  object.list[c("Control", "Pre-Clinical")],
  add.names   = c("Control", "PreClinical"),
  cell.prefix = TRUE
)
pair_CC <- mergeCellChat(
  object.list[c("Control", "Clinical")],
  add.names   = c("Control", "Clinical"),
  cell.prefix = TRUE
)
#### 1 · Differential expression on merged pair ──────────────────────────
features.name <- "PreClinical.merged"

pair_CP <- identifyOverExpressedGenes(
  object        = pair_CP,
  group.dataset = "datasets",
  pos.dataset   = "PreClinical",      # the condition with ↑ fold‐change
  features.name = features.name,
  only.pos      = FALSE,
  thresh.pc     = 0.1,                # percent‐expressed cutoff
  thresh.fc     = 0.05,               # logFC cutoff
  thresh.p      = 0.05,
  group.DE.combined = FALSE
)

# map DE results onto your communication network
net <- netMappingDEG(pair_CP, features.name = features.name, variable.all = TRUE)

# pull out only the ligand‐receptor pairs up in PreClinical vs Control
net.up.Pre <- subsetCommunication(
  object   = pair_CP,
  net      = net,
  datasets = "PreClinical",
  ligand.logFC    = 0.05,
  receptor.logFC  = NULL
)

# rename for clarity
cc_Pre <- object.list[["Pre-Clinical"]]

pair_CC <- identifyOverExpressedGenes(
  object        = pair_CC,           # your Control+Clinical merge
  group.dataset = "datasets",
  pos.dataset   = "Clinical",        # condition with ↑ fold-change
  features.name = features.name,
  only.pos      = FALSE,
  thresh.pc     = 0.1,
  thresh.fc     = 0.05,
  thresh.p      = 0.05,
  group.DE.combined = FALSE
)

# map DE onto the communication network
net.Cli    <- netMappingDEG(pair_CC, features.name, variable.all = TRUE)

# extract only the LR pairs up in Clinical vs Control
net.up.Clin <- subsetCommunication(
  object   = pair_CC,
  net      = net.Cli,
  datasets = "Clinical",
  ligand.logFC   = 0.05,
  receptor.logFC = NULL
)


#### 2 · Gene-level chord: mOSN + iOSN → all other types ─────────────────────
# single-condition CellChat object for Clinical
cc_Clin <- object.list[["Clinical"]]

# define sources & targets on pair_CP itself
cp_levels <- levels(pair_CP@idents)

pair_CP <- mergeCellChat(
  object.list[c("Control","Pre-Clinical")],
  add.names   = c("Control","PreClinical"),
  cell.prefix = FALSE     
)

library(circlize)

# 1) prep
cc_Pre <- object.list[["Pre-Clinical"]]
circos.clear()

cp_levels <- levels(cc_Pre@idents)
#OSNs as targets
tgt       <- which(cp_levels %in% c("mOSN","iOSN"))
src       <- which(! cp_levels %in% c("mOSN","iOSN","INP"))

# 3) chord plot of DE’d LR pairs
circos.clear()
netVisual_chord_gene(
  object      = cc_Pre,
  sources.use = src,
  targets.use = tgt,
  slot.name   = "net",
  net         = net.up.Pre,
  lab.cex     = 0.8,
  small.gap   = 1,       # much smaller per-sector gap
  title.name  = "Up in Pre-Clinical\n(all → mOSN+iOSN )"
)

#Repeat for OSNs as sources, and for Control vs. Clin 
