gs2eg <- function(gs, chr=NULL, strand=NULL, startp=NULL, endp=NULL)
  # This function maps gene symbol to entrez ID
  # The mapping takes several steps
  # 1. map the gene symbol as official gene symbol to entrez ID
  # 2. for entenz IDs that are not mapped in the first step, map them as alias to entrez ID
  # In ether step if there are multiple matches, use chr, strand, and position to refine the mapping.
  # If after refining, there are still multiple matches, then use '_&_' to separate the matched entrez IDs.
  
  # ****** Need to change MT in chr to M before matching
  # ****** Need to make sure that endp >= startp before using this function.  
  
{
  ori_gs = gs
  
  #load("~/Desktop/Projects/Pilot1/CTRP_Prediction/ProcessedData/eg2CHR_singleLocation.rda")
  
  library(org.Hs.eg.db)
  offGS2EG <- org.Hs.egSYMBOL2EG
  mapped_genes <- mappedkeys(offGS2EG)
  offGS2EG <- as.list(offGS2EG[mapped_genes])
  offGS2EG[["HBD"]] = "3045"
  offGS2EG[["RNR1"]] = "6052"
  offGS2EG[["RNR2"]] = "6053"
  offGS2EG[["TEC"]] = "7006"  
  offGS2EG[["MEMO1"]] = "51072"
  offGS2EG[["MMD2"]] = "221938"  
  offGS2EG[["DEL11P13"]] = "100528024"
  offGS2EG[["DEL1P36"]] = "100240737"
  offGS2EG[["SMIM44"]] = "122405565"
  offGS2EG[["TRNAK-CUU"]] = "107984908"  
  offGS2EG[["TRNAG-CCC"]] = "107985600"
  offGS2EG[["TRNAN-GUU"]] = "107985601"
  offGS2EG[["TRNAE-UUC"]] = "107985602"
  offGS2EG[["TRNAC-GCA"]] = "107985606"
  offGS2EG[["TRNAV-CAC"]] = "107985608"
  offGS2EG[["TRNAQ-CUG"]] = "107985616"
  offGS2EG[["TRNAG-GCC"]] = "107985618"
  offGS2EG[["TRNAI-AAU"]] = "107986679"
  offGS2EG[["TRNAM-CAU"]] = "107986680"
  offGS2EG[["TRNAL-CAA"]] = "107986681"
  offGS2EG[["RNA5-8SP"]] = "110467518"
  offGS2EG[["RNA28SP"]] = "110467525"
  offGS2EG[["RNA18SP"]] = "110467532"
  
  offGS2EG[["AZI1"]] = "22994"
  offGS2EG[["C1ORF220"]] = "400798"
  
  if (sum(duplicated(names(offGS2EG))) > 0)
  {
    stop("offGS2EG has duplicated names")
  }
  if (length(unique(sapply(offGS2EG, length))) > 1)
  {
    stop("some gene maps to more than one entrez ID in offGS2EG")
  }  
  alias2EG <- org.Hs.egALIAS2EG
  mapped_genes <- mappedkeys(alias2EG)
  alias2EG <- as.list(alias2EG[mapped_genes])
  if (sum(duplicated(names(alias2EG))) > 0)
  {
    stop("alias2EG has duplicated names")
  }
  
  if ((!is.null(chr)) & (is.null(strand)))
  {
    ref = paste(eg2CHR[, "entrezID"], eg2CHR[, "chr"], sep = "|")
  }  
  if ((!is.null(chr)) & (!is.null(strand)))
  {
    ref = paste(eg2CHR[, "entrezID"], eg2CHR[, "chr"], eg2CHR[, "strand"], sep = "|")
  }
  
  entrezID = rep("", length(gs))
  for (i in 1:length(gs))
  {
    idi = which(names(offGS2EG) == gs[i])
    if (length(idi) == 0)
    {
      idii = which(names(alias2EG) == gs[i])
      if (length(idii) == 0)
      {
        next()
      }
      
      mapped = alias2EG[[gs[i]]]
      if (length(mapped) == 1)
      {
        entrezID[i] = mapped   
      } else {
        if (is.null(chr))
        {
          entrezID[i] = paste(mapped, collapse = "_&_")
          next()
        }
        if ((!is.null(chr)) & (is.null(strand)))
        {
          toMatch = paste(mapped, rep(chr[i], length(mapped)), sep = "|")
        }  
        if ((!is.null(chr)) & (!is.null(strand)))
        {
          toMatch = paste(mapped, rep(chr[i], length(mapped)), rep(strand[i], length(mapped)), sep = "|")
        }        
        idj = which(ref %in% toMatch)
        if (length(idj) == 0)
        {
          next()
        }
        if (length(idj) == 1)
        {
          entrezID[i] = eg2CHR[idj, "entrezID"]
        } else {
          if ((!is.null(startp)) & (!is.null(endp)))
          {
            # p = (as.mumeric(startp[i])+as.numeric(endp[i]))/2
            # idj = idj[order(abs((as.numeric(eg2CHR[idj, "startp"]) + as.numeric(eg2CHR[idj, "endp"]))/2 - rep(p, length(idj))))[1]]
            dis = rep(NA, length(idj))
            for (j in 1:length(idj))
            {
              if (as.numeric(startp[i]) > as.numeric(eg2CHR[idj[j], "endp"]))
              {
                dis[j] = as.numeric(startp[i]) - as.numeric(eg2CHR[idj[j], "endp"])
              } else if (as.numeric(endp[i]) < as.numeric(eg2CHR[idj[j], "startp"]))
              {
                dis[j] = as.numeric(eg2CHR[idj[j], "startp"]) - as.numeric(endp[i])
              } else {
                dis[j] = 0
              }
            }
            idj = idj[which(dis == min(dis))]
          } 
          entrezID[i] = paste(eg2CHR[idj, "entrezID"], collapse = "_&_")
        }
      }    
    } else {
      entrezID[i] = offGS2EG[[gs[i]]]
    }
  }
  
  # "DEC1"      "8553_&_50514" 
  id = intersect(which(ori_gs == "DEC1"), which(entrezID == "8553_&_50514")) 
  entrezID[id] = "50514"
  
  # "DUSP27"    "92235_&_338599" 
  id = intersect(which(ori_gs == "DUSP27"), which(entrezID == "92235_&_338599")) 
  entrezID[id] = "92235"
  
  # "GARS"      "2617_&_2618"       
  id = intersect(which(ori_gs == "GARS"), which(entrezID == "2617_&_2618")) 
  entrezID[id] = "2617"
  
  # "GIF"       "2694_&_4282_&_4504"
  id = intersect(which(ori_gs == "GIF"), which(entrezID == "2694_&_4282_&_4504")) 
  entrezID[id] = "2694"
  
  # "HIST1H2BC" "3017_&_8347" 
  id = intersect(which(ori_gs == "HIST1H2BC"), which(entrezID == "3017_&_8347")) 
  entrezID[id] = "8347"
  
  # "LOR"       "4014_&_4017"    
  id = intersect(which(ori_gs == "LOR"), which(entrezID == "4014_&_4017")) 
  entrezID[id] = "4014"
  
  # "MARS"      "4141_&_84174" 
  id = intersect(which(ori_gs == "MARS"), which(entrezID == "4141_&_84174")) 
  entrezID[id] = "4141"
  
  # "MUM1"      "3662_&_84939" 
  id = intersect(which(ori_gs == "MUM1"), which(entrezID == "3662_&_84939")) 
  entrezID[id] = "84939"
  
  # "NAT6"      "3373_&_24142"      
  id = intersect(which(ori_gs == "NAT6"), which(entrezID == "3373_&_24142")) 
  entrezID[id] = "24142"
  
  # "NOV"       "4856_&_5361_&_6134"
  id = intersect(which(ori_gs == "NOV"), which(entrezID == "4856_&_5361_&_6134")) 
  entrezID[id] = "4856"
  
  # "QARS"      "2058_&_5859"    
  id = intersect(which(ori_gs == "QARS"), which(entrezID == "2058_&_5859")) 
  entrezID[id] = "5859"
  
  # "SARS"      "6301_&_54938"      
  id = intersect(which(ori_gs == "SARS"), which(entrezID == "6301_&_54938")) 
  entrezID[id] = "6301"
  
  # "SEPT2"     "4735_&_23157"      
  id = intersect(which(ori_gs == "SEPT2"), which(entrezID == "4735_&_23157")) 
  entrezID[id] = "4735"
  
  # "SLC35E2"   "9906_&_728661" 
  id = intersect(which(ori_gs == "SLC35E2"), which(entrezID == "9906_&_728661")) 
  entrezID[id] = "9906"
  
  # "FAM21A"    "253725_&_387680"               
  id = intersect(which(ori_gs == "FAM21A"), which(entrezID == "253725_&_387680")) 
  entrezID[id] = "387680"
  
  # "C10orf2"   "56652_&_57132"                 
  id = intersect(which(ori_gs == "C10orf2"), which(entrezID == "56652_&_57132")) 
  entrezID[id] = "56652"
  
  # "AIM1"      "202_&_9212_&_51151"            
  id = intersect(which(ori_gs == "AIM1"), which(entrezID == "202_&_9212_&_51151")) 
  entrezID[id] = "202"
  
  # "CD97"      "976_&_30817"                   
  id = intersect(which(ori_gs == "CD97"), which(entrezID == "976_&_30817")) 
  entrezID[id] = "976"
  
  # "SHFM1"     "1749_&_7979"                   
  id = intersect(which(ori_gs == "SHFM1"), which(entrezID == "1749_&_7979")) 
  entrezID[id] = "7979"
  
  # "AGPAT9"    "79888_&_84803"                 
  id = intersect(which(ori_gs == "AGPAT9"), which(entrezID == "79888_&_84803")) 
  entrezID[id] = "84803"
  
  # "EFTUD1"    "1997_&_79631"                  
  id = intersect(which(ori_gs == "EFTUD1"), which(entrezID == "1997_&_79631")) 
  entrezID[id] = "79631"
  
  # "ADC"       "113451_&_339896"               
  id = intersect(which(ori_gs == "ADC"), which(entrezID == "113451_&_339896")) 
  entrezID[id] = "113451"
  
  # "CSRP2BP"   "57325_&_100303755"             
  id = intersect(which(ori_gs == "CSRP2BP"), which(entrezID == "57325_&_100303755")) 
  entrezID[id] = "57325"
  
  # "C11orf48"  "79081_&_102288414"             
  id = intersect(which(ori_gs == "C11orf48"), which(entrezID == "79081_&_102288414")) 
  entrezID[id] = "79081"
  
  # "C2orf47"   "51614_&_79568"                 
  id = intersect(which(ori_gs == "C2orf47"), which(entrezID == "51614_&_79568")) 
  entrezID[id] = "79568"
  
  # "C7orf55"   "154791_&_100996928"            
  id = intersect(which(ori_gs == "C7orf55"), which(entrezID == "154791_&_100996928")) 
  entrezID[id] = "154791"
  
  # "HIST1H2BI" "3017_&_8346"                   
  id = intersect(which(ori_gs == "HIST1H2BI"), which(entrezID == "3017_&_8346")) 
  entrezID[id] = "8346"
  
  # "STRA13"    "8553_&_201254"                 
  id = intersect(which(ori_gs == "STRA13"), which(entrezID == "8553_&_201254")) 
  entrezID[id] = "201254"
  
  # "B3GNT1"    "10678_&_11041"                 
  id = intersect(which(ori_gs == "B3GNT1"), which(entrezID == "10678_&_11041")) 
  entrezID[id] = "11041"
  
  # "APITD1"    "378708_&_100526739"            
  id = intersect(which(ori_gs == "APITD1"), which(entrezID == "378708_&_100526739")) 
  entrezID[id] = "378708"
  
  # "DDC8"      "7077_&_100653515"              
  id = intersect(which(ori_gs == "DDC8"), which(entrezID == "7077_&_100653515")) 
  entrezID[id] = "100653515"
  
  # "TCEB3C"    "162699_&_728929"               
  id = intersect(which(ori_gs == "TCEB3C"), which(entrezID == "162699_&_728929")) 
  entrezID[id] = "162699"
  
  # "C2orf27B"  "29798_&_408029"                
  id = intersect(which(ori_gs == "C2orf27B"), which(entrezID == "29798_&_408029")) 
  entrezID[id] = "408029"
  
  # "HIST1H2BG" "3017_&_8339"                   
  id = intersect(which(ori_gs == "HIST1H2BG"), which(entrezID == "3017_&_8339")) 
  entrezID[id] = "8339"
  
  # "HN1"       "51155_&_100462977"             
  id = intersect(which(ori_gs == "HN1"), which(entrezID == "51155_&_100462977")) 
  entrezID[id] = "51155"
  
  # "HIST1H2BE" "3017_&_8344"                   
  id = intersect(which(ori_gs == "HIST1H2BE"), which(entrezID == "3017_&_8344")) 
  entrezID[id] = "8344"
  
  # "HIST1H2BF" "3017_&_8343"                   
  id = intersect(which(ori_gs == "HIST1H2BF"), which(entrezID == "3017_&_8343")) 
  entrezID[id] = "8343"
  
  # "HN1L"      "57585_&_90861"                 
  id = intersect(which(ori_gs == "HN1L"), which(entrezID == "57585_&_90861")) 
  entrezID[id] = "90861"
  
  # "NOTCH2NL"  "388677_&_100996717_&_100996763"
  id = intersect(which(ori_gs == "NOTCH2NL"), which(entrezID == "388677_&_100996717_&_100996763")) 
  entrezID[id] = "388677"
  
  # "TCEB3CL"   "728929_&_107983955"            
  id = intersect(which(ori_gs == "TCEB3CL"), which(entrezID == "728929_&_107983955")) 
  entrezID[id] = "728929"
  
  # "C15orf38"  "348110_&_100526783"            
  id = intersect(which(ori_gs == "C15orf38"), which(entrezID == "348110_&_100526783")) 
  entrezID[id] = "348110"
  
  # "TRAPPC2P1" "6399_&_10597_&_284306"         
  id = intersect(which(ori_gs == "TRAPPC2P1"), which(entrezID == "6399_&_10597_&_284306")) 
  entrezID[id] = "10597"
  
  # "LIMS3L"    "96626_&_100288695"             
  id = intersect(which(ori_gs == "LIMS3L"), which(entrezID == "96626_&_100288695")) 
  entrezID[id] = "100288695"
  
  # "ATP6C"     "527_&_528"   
  id = intersect(which(ori_gs == "ATP6C"), which(entrezID == "527_&_528")) 
  entrezID[id] = "527"
  
  # "ABP1"     "26_&_28988"   
  id = intersect(which(ori_gs == "ABP1"), which(entrezID == "26_&_28988")) 
  entrezID[id] = "26"
  
  # "DBC1"     "1620_&_57805"
  id = intersect(which(ori_gs == "DBC1"), which(entrezID == "1620_&_57805")) 
  entrezID[id] = "1620"
  
  # "CCRL1"    "1524_&_51554"
  id = intersect(which(ori_gs == "CCRL1"), which(entrezID == "1524_&_51554")) 
  entrezID[id] = "51554"
  
  # "MLL2"     "8085_&_9757" 
  id = intersect(which(ori_gs == "MLL2"), which(entrezID == "8085_&_9757")) 
  entrezID[id] = "8085"
  
  # "CXXC11"   "54620_&_285093"
  id = intersect(which(ori_gs == "CXXC11"), which(entrezID == "54620_&_285093")) 
  entrezID[id] = "285093"
  
  # "GGTA1P"   "2681_&_121328"  
  id = intersect(which(ori_gs == "GGTA1P"), which(entrezID == "2681_&_121328")) 
  entrezID[id] = "2681"
  
  # "PPYR1"    "5540_&_100996758"  
  id = intersect(which(ori_gs == "PPYR1"), which(entrezID == "5540_&_100996758")) 
  entrezID[id] = "5540"
  
  id = intersect(which(ori_gs == "MEGT1"), which(entrezID == "58530_&_110599563")) 
  entrezID[id] = "110599563"
  
  #id = intersect(which(ori_gs == "LST3"), which(entrezID == "28234_&_84173_&_338821")) 
  #entrezID[id] = "110599563"
  id = intersect(which(ori_gs == "MPP6"), which(entrezID == "10200_&_51678")) 
  entrezID[id] = "51678"
  
  id = intersect(which(ori_gs == "GGT2"), which(entrezID == "728441_&_102724197")) 
  entrezID[id] = "728441"
  
  id = intersect(which(ori_gs == "TAZ"), which(entrezID == "6901_&_25937")) 
  entrezID[id] = "6901"
  
  id = intersect(which(ori_gs == "CT45A4"), which(entrezID == "441519_&_441521")) 
  entrezID[id] = "441519"
  
  id = intersect(which(ori_gs == "SMIM11B"), which(entrezID == "54065_&_102723553")) 
  entrezID[id] = "102723553"
  
  id = intersect(which(ori_gs == "CBSL"), which(entrezID == "875_&_102724560")) 
  entrezID[id] = "102724560"
  
  id = grep("&", entrezID)
  mapping = cbind(GeneSymbol = ori_gs[id], EntrezID = entrezID[id])
  print(mapping)
  
  entrezID
}

# Function to convert entrz id back to gene symbol 
eg2gs <- function(eg)
{
  library(org.Hs.eg.db)
  eg2GS <- org.Hs.egSYMBOL
  mapped_genes <- mappedkeys(eg2GS)
  eg2GS <- as.list(eg2GS[mapped_genes])
  eg2GS[["4550"]] = "MT-RNR2"
  eg2GS[["4549"]] = "MT-RNR1"
  
  gs = rep("", length(eg))
  for (i in seq(1, length(eg)))
  {
    if ((i %% ceiling(length(eg)/20) == 0) & (length(eg) > 100000))
    {
      print(paste(round(i/length(eg), 2)*100, " percents finished.", sep = ""))
    }
    
    if (length(grep("_&_", eg[i])) > 0)
    {
      eg_s = strsplit(eg[i], split = "_&_")[[1]]
      eg_s = eg_s[which(eg_s %in% names(eg2GS))]
      if (length(eg_s) > 0)
      {
        gs[i] = paste(as.character(eg2GS[eg_s]), collapse = "_&_")  
      }
    } else if (eg[i] %in% names(eg2GS))
    {
      #gs[i] = eg2GS[[eg[i]]]  
      gs[i] = eg2GS[toString(eg[i])]
    }
  }
  gs
}


# Function to convert Ensembl ID to Entrez ID
ensembl2eg <- function(ensembl_list)
{
  library(org.Hs.eg.db)
  x = org.Hs.egENSEMBL
  mapped_genes <- mappedkeys(x)
  xx <- as.list(x[mapped_genes])
  ensembl2entrez = c()
  ensembl2entrez_name = c()
  for (i in seq(1, length(xx)))
  {
    ensembl2entrez = c(ensembl2entrez, rep(names(xx)[i], length(xx[[i]])))
    ensembl2entrez_name = c(ensembl2entrez_name, xx[[i]])
  }
  
  ensembl = unique(ensembl2entrez_name)
  entrezID = rep("", length(ensembl))
  for (i in seq(1, length(ensembl)))
  {
    idi = which(ensembl2entrez_name == ensembl[i])
    if (length(idi) > 1)
    {
      entrezID[i] = paste(ensembl2entrez[idi], collapse = "_&_")
    } else {
      entrezID[i] = ensembl2entrez[idi]
    }
  }
  mapping = cbind(entrezID = entrezID, ensembl = ensembl)
  rownames(mapping) = mapping[, "ensembl"]
  
  entrezID = rep("", length(ensembl_list))
  id = which(ensembl_list %in% rownames(mapping))
  entrezID[id] = mapping[ensembl_list[id], "entrezID"]
  
  # id = grep("&", entrezID)
  # mapping = cbind(EnsembleID = ensembl_list[id], EntrezID = entrezID[id])
  # print(mapping)
  
  # [1,] "ENSG00000004866" "7982_&_93655"   
  id = intersect(which(ensembl_list == "ENSG00000004866"), which(entrezID == "7982_&_93655")) 
  entrezID[id] = "7982"
  
  # [2,] "ENSG00000048545" "118142757_&_2978"  
  id = intersect(which(ensembl_list == "ENSG00000048545"), which(entrezID == "118142757_&_2978")) 
  entrezID[id] = "2978"
  
  # [3,] "ENSG00000063587" "105373378_&_10838"                                    
  id = intersect(which(ensembl_list == "ENSG00000063587"), which(entrezID == "105373378_&_10838")) 
  entrezID[id] = "10838"
  
  # [4,] "ENSG00000075790" "115253422_&_55973"                                    
  id = intersect(which(ensembl_list == "ENSG00000075790"), which(entrezID == "115253422_&_55973")) 
  entrezID[id] = "55973"
  
  # [5,] "ENSG00000076928" "100505585_&_9138"                                     
  id = intersect(which(ensembl_list == "ENSG00000076928"), which(entrezID == "100505585_&_9138")) 
  entrezID[id] = "9138"
  
  # [6,] "ENSG00000086200" "101180901_&_51194"                                    
  id = intersect(which(ensembl_list == "ENSG00000086200"), which(entrezID == "101180901_&_51194")) 
  entrezID[id] = "51194"
  
  # [7,] "ENSG00000088298" "111089941_&_55741"                                    
  id = intersect(which(ensembl_list == "ENSG00000088298"), which(entrezID == "111089941_&_55741")) 
  entrezID[id] = "55741"
  
  # [8,] "ENSG00000099977" "100037417_&_1652"                                     
  id = intersect(which(ensembl_list == "ENSG00000099977"), which(entrezID == "100037417_&_1652")) 
  entrezID[id] = "1652"
  
  # [9,] "ENSG00000100197" "107987478_&_107987479_&_1565"                         
  id = intersect(which(ensembl_list == "ENSG00000100197"), which(entrezID == "107987478_&_107987479_&_1565")) 
  entrezID[id] = "1565"
  
  # [10,] "ENSG00000105793" "101927446_&_85865"                                    
  id = intersect(which(ensembl_list == "ENSG00000105793"), which(entrezID == "101927446_&_85865")) 
  entrezID[id] = "85865"
  
  # [11,] "ENSG00000105887" "136319_&_767558"                                      
  id = intersect(which(ensembl_list == "ENSG00000105887"), which(entrezID == "136319_&_767558")) 
  entrezID[id] = "767558"
  
  # [12,] "ENSG00000105889" "256227_&_401312"                                      
  id = intersect(which(ensembl_list == "ENSG00000105889"), which(entrezID == "256227_&_401312")) 
  entrezID[id] = "256227"
  
  # [13,] "ENSG00000107014" "6013_&_6019"                                          
  id = intersect(which(ensembl_list == "ENSG00000107014"), which(entrezID == "6013_&_6019")) 
  entrezID[id] = "6019"
  
  # [14,] "ENSG00000108379" "101929777_&_7473"                                     
  id = intersect(which(ensembl_list == "ENSG00000108379"), which(entrezID == "101929777_&_7473")) 
  entrezID[id] = "7473"
  
  # [15,] "ENSG00000111215" "11272_&_5554"                                         
  id = intersect(which(ensembl_list == "ENSG00000111215"), which(entrezID == "11272_&_5554")) 
  entrezID[id] = "11272"
  
  # [16,] "ENSG00000111671" "105369632_&_84727"                                    
  id = intersect(which(ensembl_list == "ENSG00000111671"), which(entrezID == "105369632_&_84727")) 
  entrezID[id] = "84727"
  
  # [17,] "ENSG00000111850" "57150_&_63914"                                        
  id = intersect(which(ensembl_list == "ENSG00000111850"), which(entrezID == "57150_&_63914")) 
  entrezID[id] = "57150"
  
  # [18,] "ENSG00000112541" "10846_&_90632"                                        
  id = intersect(which(ensembl_list == "ENSG00000112541"), which(entrezID == "10846_&_90632")) 
  entrezID[id] = "10846"
  
  # [19,] "ENSG00000114316" "107986084_&_7375"                                     
  id = intersect(which(ensembl_list == "ENSG00000114316"), which(entrezID == "107986084_&_7375")) 
  entrezID[id] = "7375"
  
  # [20,] "ENSG00000114374" "64595_&_8287"                                         
  id = intersect(which(ensembl_list == "ENSG00000114374"), which(entrezID == "64595_&_8287")) 
  entrezID[id] = "8287"
  
  # [21,] "ENSG00000115239" "100302652_&_51130"                                    
  id = intersect(which(ensembl_list == "ENSG00000115239"), which(entrezID == "100302652_&_51130")) 
  entrezID[id] = "51130"
  
  # [22,] "ENSG00000120088" "104909134_&_1394"                                     
  id = intersect(which(ensembl_list == "ENSG00000120088"), which(entrezID == "104909134_&_1394")) 
  entrezID[id] = "1394"
  
  # [23,] "ENSG00000120341" "111240474_&_89866"                                    
  id = intersect(which(ensembl_list == "ENSG00000120341"), which(entrezID == "111240474_&_89866")) 
  entrezID[id] = "89866"
  
  # [24,] "ENSG00000120709" "100128966_&_51307"                                    
  id = intersect(which(ensembl_list == "ENSG00000120709"), which(entrezID == "100128966_&_51307")) 
  entrezID[id] = "51307"
  
  # [25,] "ENSG00000124343" "100132596_&_7499"                                     
  id = intersect(which(ensembl_list == "ENSG00000124343"), which(entrezID == "100132596_&_7499")) 
  entrezID[id] = "7499"
  
  # [26,] "ENSG00000124713" "107080644_&_27232"                                    
  id = intersect(which(ensembl_list == "ENSG00000124713"), which(entrezID == "107080644_&_27232")) 
  entrezID[id] = "27232"
  
  # [27,] "ENSG00000125498" "112267881_&_3802"                                     
  id = intersect(which(ensembl_list == "ENSG00000125498"), which(entrezID == "112267881_&_3802")) 
  entrezID[id] = "3802"
  
  # [28,] "ENSG00000126861" "101927057_&_4974"                                     
  id = intersect(which(ensembl_list == "ENSG00000126861"), which(entrezID == "101927057_&_4974")) 
  entrezID[id] = "4974"
  
  # [29,] "ENSG00000128383" "100913187_&_200315"                                   
  id = intersect(which(ensembl_list == "ENSG00000128383"), which(entrezID == "100913187_&_200315")) 
  entrezID[id] = "200315"
  
  # [30,] "ENSG00000130703" "105369209_&_9885"                                     
  id = intersect(which(ensembl_list == "ENSG00000130703"), which(entrezID == "105369209_&_9885")) 
  entrezID[id] = "9885"
  
  # [31,] "ENSG00000131737" "100653049_&_3885"                                     
  id = intersect(which(ensembl_list == "ENSG00000131737"), which(entrezID == "100653049_&_3885")) 
  entrezID[id] = "3885"
  
  # [32,] "ENSG00000133816" "84953_&_9645"                                         
  id = intersect(which(ensembl_list == "ENSG00000133816"), which(entrezID == "84953_&_9645")) 
  entrezID[id] = "9645"
  
  # [33,] "ENSG00000137843" "106821730_&_56924"                                    
  id = intersect(which(ensembl_list == "ENSG00000137843"), which(entrezID == "106821730_&_56924")) 
  entrezID[id] = "56924"
  
  # [34,] "ENSG00000137936" "723788_&_8412"                                        
  id = intersect(which(ensembl_list == "ENSG00000137936"), which(entrezID == "723788_&_8412")) 
  entrezID[id] = "8412"
  
  # [35,] "ENSG00000138594" "112268148_&_29766"                                    
  id = intersect(which(ensembl_list == "ENSG00000138594"), which(entrezID == "112268148_&_29766")) 
  entrezID[id] = "29766"
  
  # [36,] "ENSG00000143702" "645455_&_9859"                                        
  id = intersect(which(ensembl_list == "ENSG00000143702"), which(entrezID == "645455_&_9859")) 
  entrezID[id] = "9859"
  
  # [37,] "ENSG00000144455" "100130207_&_285362"                                   
  id = intersect(which(ensembl_list == "ENSG00000144455"), which(entrezID == "100130207_&_285362")) 
  entrezID[id] = "285362"
  
  # [38,] "ENSG00000145736" "2966_&_728340"                                        
  id = intersect(which(ensembl_list == "ENSG00000145736"), which(entrezID == "2966_&_728340")) 
  entrezID[id] = "2966"
  
  # [39,] "ENSG00000145979" "107080638_&_51256"                                    
  id = intersect(which(ensembl_list == "ENSG00000145979"), which(entrezID == "107080638_&_51256")) 
  entrezID[id] = "51256"
  
  # [40,] "ENSG00000146112" "107987457_&_170954"                                   
  id = intersect(which(ensembl_list == "ENSG00000146112"), which(entrezID == "107987457_&_170954")) 
  entrezID[id] = "170954"
  
  # [41,] "ENSG00000149557" "105369550_&_9638"                                     
  id = intersect(which(ensembl_list == "ENSG00000149557"), which(entrezID == "105369550_&_9638")) 
  entrezID[id] = "9638"
  
  # [42,] "ENSG00000151327" "101927178_&_283635"                                   
  id = intersect(which(ensembl_list == "ENSG00000151327"), which(entrezID == "101927178_&_283635")) 
  entrezID[id] = "283635"
  
  # [43,] "ENSG00000152926" "109504726_&_51351"                                    
  id = intersect(which(ensembl_list == "ENSG00000152926"), which(entrezID == "109504726_&_51351")) 
  entrezID[id] = "51351"
  
  # [44,] "ENSG00000153071" "112267931_&_1601"                                     
  id = intersect(which(ensembl_list == "ENSG00000153071"), which(entrezID == "112267931_&_1601")) 
  entrezID[id] = "1601"
  
  # [45,] "ENSG00000155307" "388813_&_64092"                                       
  id = intersect(which(ensembl_list == "ENSG00000155307"), which(entrezID == "388813_&_64092")) 
  entrezID[id] = "64092"
  
  # [46,] "ENSG00000156273" "100379661_&_571"                                      
  id = intersect(which(ensembl_list == "ENSG00000156273"), which(entrezID == "100379661_&_571")) 
  entrezID[id] = "571"
  
  # [47,] "ENSG00000157224" "102723899_&_9069"                                     
  id = intersect(which(ensembl_list == "ENSG00000157224"), which(entrezID == "102723899_&_9069")) 
  entrezID[id] = "9069"
  
  # [48,] "ENSG00000158301" "100528062_&_114928"                                   
  id = intersect(which(ensembl_list == "ENSG00000158301"), which(entrezID == "100528062_&_114928")) 
  entrezID[id] = "114928"
  
  # [49,] "ENSG00000158747" "100532736_&_4681"                                     
  id = intersect(which(ensembl_list == "ENSG00000158747"), which(entrezID == "100532736_&_4681")) 
  entrezID[id] = "4681"
  
  # [50,] "ENSG00000159216" "100506403_&_861"                                      
  id = intersect(which(ensembl_list == "ENSG00000159216"), which(entrezID == "100506403_&_861")) 
  entrezID[id] = "861"
  
  # [51,] "ENSG00000160209" "105372824_&_8566"                                     
  id = intersect(which(ensembl_list == "ENSG00000160209"), which(entrezID == "105372824_&_8566")) 
  entrezID[id] = "8566"
  
  # [52,] "ENSG00000161270" "107985317_&_4868"                                     
  id = intersect(which(ensembl_list == "ENSG00000161270"), which(entrezID == "107985317_&_4868")) 
  entrezID[id] = "4868"
  
  # [53,] "ENSG00000163156" "100534012_&_79005"                                    
  id = intersect(which(ensembl_list == "ENSG00000163156"), which(entrezID == "100534012_&_79005")) 
  entrezID[id] = "79005"
  
  # [54,] "ENSG00000164053" "111822955_&_84126"                                    
  id = intersect(which(ensembl_list == "ENSG00000164053"), which(entrezID == "111822955_&_84126")) 
  entrezID[id] = "84126"
  
  # [55,] "ENSG00000165269" "112267859_&_364"                                      
  id = intersect(which(ensembl_list == "ENSG00000165269"), which(entrezID == "112267859_&_364")) 
  entrezID[id] = "364"
  
  # [56,] "ENSG00000166272" "102724307_&_54838"                                    
  id = intersect(which(ensembl_list == "ENSG00000166272"), which(entrezID == "102724307_&_54838")) 
  entrezID[id] = "54838"
  
  # [57,] "ENSG00000166848" "105371348_&_54386"                                    
  id = intersect(which(ensembl_list == "ENSG00000166848"), which(entrezID == "105371348_&_54386")) 
  entrezID[id] = "54386"
  
  # [58,] "ENSG00000167476" "105372240_&_126306"                                   
  id = intersect(which(ensembl_list == "ENSG00000167476"), which(entrezID == "105372240_&_126306")) 
  entrezID[id] = "126306"
  
  # [59,] "ENSG00000168746" "101927242_&_140834"                                   
  id = intersect(which(ensembl_list == "ENSG00000168746"), which(entrezID == "101927242_&_140834")) 
  entrezID[id] = "140834"
  
  # [60,] "ENSG00000168802" "113455421_&_54921"                                    
  id = intersect(which(ensembl_list == "ENSG00000168802"), which(entrezID == "113455421_&_54921")) 
  entrezID[id] = "54921"
  
  # [61,] "ENSG00000169627" "552900_&_654483"                                      
  id = intersect(which(ensembl_list == "ENSG00000169627"), which(entrezID == "552900_&_654483")) 
  entrezID[id] = "654483"
  
  # [62,] "ENSG00000171711" "100289462_&_1673"                                     
  id = intersect(which(ensembl_list == "ENSG00000171711"), which(entrezID == "100289462_&_1673")) 
  entrezID[id] = "1673"
  
  # [63,] "ENSG00000176124" "100874074_&_10301"                                    
  id = intersect(which(ensembl_list == "ENSG00000176124"), which(entrezID == "100874074_&_10301")) 
  entrezID[id] = "10301"
  
  # [64,] "ENSG00000176797" "414325_&_55894"                                       
  id = intersect(which(ensembl_list == "ENSG00000176797"), which(entrezID == "414325_&_55894")) 
  entrezID[id] = "414325"
  
  # [65,] "ENSG00000177257" "100289462_&_1673"                                     
  id = intersect(which(ensembl_list == "ENSG00000177257"), which(entrezID == "100289462_&_1673")) 
  entrezID[id] = "100289462"
  
  # [66,] "ENSG00000178115" "727909_&_728047"                                      
  id = intersect(which(ensembl_list == "ENSG00000178115"), which(entrezID == "727909_&_728047")) 
  entrezID[id] = "727909"
  
  # [67,] "ENSG00000178882" "100533183_&_144347"                                   
  id = intersect(which(ensembl_list == "ENSG00000178882"), which(entrezID == "100533183_&_144347")) 
  entrezID[id] = "144347"
  
  # [68,] "ENSG00000178934" "3963_&_653499"                                        
  id = intersect(which(ensembl_list == "ENSG00000178934"), which(entrezID == "3963_&_653499")) 
  entrezID[id] = "653499"
  
  # [69,] "ENSG00000181135" "101928160_&_286075"                                   
  id = intersect(which(ensembl_list == "ENSG00000181135"), which(entrezID == "101928160_&_286075")) 
  entrezID[id] = "286075"
  
  # [70,] "ENSG00000183426" "642799_&_9284"                                        
  id = intersect(which(ensembl_list == "ENSG00000183426"), which(entrezID == "642799_&_9284")) 
  entrezID[id] = "9284"
  
  # [71,] "ENSG00000183632" "102723655_&_102723713_&_102724101_&_102724127_&_24150"
  id = intersect(which(ensembl_list == "ENSG00000183632"), which(entrezID == "102723655_&_102723713_&_102724101_&_102724127_&_24150")) 
  entrezID[id] = "24150"
  
  # [72,] "ENSG00000185164" "102723728_&_283820"                                   
  id = intersect(which(ensembl_list == "ENSG00000185164"), which(entrezID == "102723728_&_283820")) 
  entrezID[id] = "283820"
  
  # [73,] "ENSG00000186047" "103689915_&_220107"                                   
  id = intersect(which(ensembl_list == "ENSG00000186047"), which(entrezID == "103689915_&_220107")) 
  entrezID[id] = "220107"
  
  # [74,] "ENSG00000186448" "10168_&_110354863"                                    
  id = intersect(which(ensembl_list == "ENSG00000186448"), which(entrezID == "10168_&_110354863")) 
  entrezID[id] = "10168"
  
  # [75,] "ENSG00000186572" "245910_&_503614"                                      
  id = intersect(which(ensembl_list == "ENSG00000186572"), which(entrezID == "245910_&_503614")) 
  entrezID[id] = "245910"
  
  # [76,] "ENSG00000186599" "245908_&_504180"                                      
  id = intersect(which(ensembl_list == "ENSG00000186599"), which(entrezID == "245908_&_504180")) 
  entrezID[id] = "504180"
  
  # [77,] "ENSG00000186854" "105374836_&_129293"                                   
  id = intersect(which(ensembl_list == "ENSG00000186854"), which(entrezID == "105374836_&_129293")) 
  entrezID[id] = "129293"
  
  # [78,] "ENSG00000189064" "100101629_&_645037_&_729447"                          
  id = intersect(which(ensembl_list == "ENSG00000189064"), which(entrezID == "100101629_&_645037_&_729447")) 
  entrezID[id] = "729447"
  
  # [1,] "ENSG00000100181" "112268292_&_387590"                                    
  id = intersect(which(ensembl_list == "ENSG00000100181"), which(entrezID == "112268292_&_387590")) 
  entrezID[id] = "387590"
  
  # [2,] "ENSG00000106610" "101929736_&_64940"                                     
  id = intersect(which(ensembl_list == "ENSG00000106610"), which(entrezID == "101929736_&_64940")) 
  entrezID[id] = "64940"
  
  # [3,] "ENSG00000160336" "110116772_&_388561"   
  id = intersect(which(ensembl_list == "ENSG00000160336"), which(entrezID == "110116772_&_388561")) 
  entrezID[id] = "388561"
  
  # [4,] "ENSG00000183292" "150527_&_646743"                                       
  id = intersect(which(ensembl_list == "ENSG00000183292"), which(entrezID == "150527_&_646743")) 
  entrezID[id] = "150527"
  
  # [5,] "ENSG00000185040" "102723555_&_285955"                                    
  id = intersect(which(ensembl_list == "ENSG00000185040"), which(entrezID == "102723555_&_285955")) 
  entrezID[id] = "102723555"
  
  # [6,] "ENSG00000188681" "100132288_&_100233156"                                 
  id = intersect(which(ensembl_list == "ENSG00000188681"), which(entrezID == "100132288_&_100233156")) 
  entrezID[id] = "100132288"
  
  # [7,] "ENSG00000196605" "100505555_&_162993"                                    
  id = intersect(which(ensembl_list == "ENSG00000196605"), which(entrezID == "100505555_&_162993")) 
  entrezID[id] = "162993"
  
  # [8,] "ENSG00000196993" "100507607_&_105379464"                                 
  id = intersect(which(ensembl_list == "ENSG00000196993"), which(entrezID == "100507607_&_105379464")) 
  entrezID[id] = "100507607"
  
  # [9,] "ENSG00000197302" "107983990_&_124411"  
  id = intersect(which(ensembl_list == "ENSG00000197302"), which(entrezID == "107983990_&_124411")) 
  entrezID[id] = "124411"
  
  # [10,] "ENSG00000197912" "101930112_&_6687"                                      
  id = intersect(which(ensembl_list == "ENSG00000197912"), which(entrezID == "101930112_&_6687")) 
  entrezID[id] = "6687"
  
  # [11,] "ENSG00000198129" "245910_&_503614"                                       
  id = intersect(which(ensembl_list == "ENSG00000198129"), which(entrezID == "245910_&_503614")) 
  entrezID[id] = "503614"
  
  # [12,] "ENSG00000203668" "1122_&_23596"                                          
  id = intersect(which(ensembl_list == "ENSG00000203668"), which(entrezID == "1122_&_23596")) 
  entrezID[id] = "1122"
  
  # [13,] "ENSG00000204131" "340527_&_392490"                                       
  id = intersect(which(ensembl_list == "ENSG00000204131"), which(entrezID == "340527_&_392490")) 
  entrezID[id] = "340527"
  
  # [14,] "ENSG00000204209" "113523636_&_1616"                                      
  id = intersect(which(ensembl_list == "ENSG00000204209"), which(entrezID == "113523636_&_1616")) 
  entrezID[id] = "1616"
  
  # [15,] "ENSG00000204314" "100507547_&_80863"                                     
  id = intersect(which(ensembl_list == "ENSG00000204314"), which(entrezID == "100507547_&_80863")) 
  entrezID[id] = "80863"
  
  # [16,] "ENSG00000204577" "102725035_&_107987425_&_107987462_&_11025"             
  id = intersect(which(ensembl_list == "ENSG00000204577"), which(entrezID == "102725035_&_107987425_&_107987462_&_11025")) 
  entrezID[id] = "11025"
  
  # [17,] "ENSG00000205076" "3963_&_653499"                                         
  id = intersect(which(ensembl_list == "ENSG00000205076"), which(entrezID == "3963_&_653499")) 
  entrezID[id] = "3963"
  
  # [18,] "ENSG00000205456" "102723655_&_102723713_&_102724101_&_102724127_&_729264"
  id = intersect(which(ensembl_list == "ENSG00000205456"), which(entrezID == "102723655_&_102723713_&_102724101_&_102724127_&_729264")) 
  entrezID[id] = "729264"
  
  # [19,] "ENSG00000205457" "102723655_&_102723713_&_102724101_&_102724127_&_653550"
  id = intersect(which(ensembl_list == "ENSG00000205457"), which(entrezID == "102723655_&_102723713_&_102724101_&_102724127_&_653550")) 
  entrezID[id] = "653550"
  
  # [20,] "ENSG00000205571" "6606_&_6607"                                           
  id = intersect(which(ensembl_list == "ENSG00000205571"), which(entrezID == "6606_&_6607")) 
  entrezID[id] = "6607"
  
  # [21,] "ENSG00000205572" "728492_&_8293"                                         
  id = intersect(which(ensembl_list == "ENSG00000205572"), which(entrezID == "728492_&_8293")) 
  entrezID[id] = "728492"
  
  # [22,] "ENSG00000211689" "445347_&_6966"                                         
  id = intersect(which(ensembl_list == "ENSG00000211689"), which(entrezID == "445347_&_6966")) 
  entrezID[id] = "6966"
  
  # [23,] "ENSG00000212127" "106707243_&_50840"                                     
  id = intersect(which(ensembl_list == "ENSG00000212127"), which(entrezID == "106707243_&_50840")) 
  entrezID[id] = "50840"
  
  # [24,] "ENSG00000213160" "100526832_&_151230"                                    
  id = intersect(which(ensembl_list == "ENSG00000213160"), which(entrezID == "100526832_&_151230")) 
  entrezID[id] = "151230"
  
  # [25,] "ENSG00000213694" "1903_&_286223"                                         
  id = intersect(which(ensembl_list == "ENSG00000213694"), which(entrezID == "1903_&_286223")) 
  entrezID[id] = "1903"
  
  # [26,] "ENSG00000213753" "100271626_&_65996"                                     
  id = intersect(which(ensembl_list == "ENSG00000213753"), which(entrezID == "100271626_&_65996")) 
  entrezID[id] = "65996"
  
  # [27,] "ENSG00000213809" "100528032_&_22914"                                     
  id = intersect(which(ensembl_list == "ENSG00000213809"), which(entrezID == "100528032_&_22914")) 
  entrezID[id] = "22914"
  
  # [28,] "ENSG00000213999" "100271849_&_4207"                                      
  id = intersect(which(ensembl_list == "ENSG00000213999"), which(entrezID == "100271849_&_4207")) 
  entrezID[id] = "100271849"
  
  # [29,] "ENSG00000214026" "107987373_&_6150"                                      
  id = intersect(which(ensembl_list == "ENSG00000214026"), which(entrezID == "107987373_&_6150")) 
  entrezID[id] = "6150"
  
  # [30,] "ENSG00000215190" "106660612_&_728411"                                    
  id = intersect(which(ensembl_list == "ENSG00000215190"), which(entrezID == "106660612_&_728411")) 
  entrezID[id] = "106660612"
  
  # [31,] "ENSG00000215269" "2576_&_2578_&_2579_&_26748_&_645073"                   
  id = intersect(which(ensembl_list == "ENSG00000215269"), which(entrezID == "2576_&_2578_&_2579_&_26748_&_645073")) 
  entrezID[id] = "645073"
  
  # [32,] "ENSG00000223443" "377630_&_401447"                                       
  id = intersect(which(ensembl_list == "ENSG00000223443"), which(entrezID == "377630_&_401447")) 
  entrezID[id] = "377630"
  
  # [33,] "ENSG00000224659" "26749_&_729396"                                        
  id = intersect(which(ensembl_list == "ENSG00000224659"), which(entrezID == "26749_&_729396")) 
  entrezID[id] = "729396"
  
  # [34,] "ENSG00000226023" "255313_&_728062"                                       
  id = intersect(which(ensembl_list == "ENSG00000226023"), which(entrezID == "255313_&_728062")) 
  entrezID[id] = "728062"
  
  # [35,] "ENSG00000228695" "107987423_&_51716"                                     
  id = intersect(which(ensembl_list == "ENSG00000228695"), which(entrezID == "107987423_&_51716")) 
  entrezID[id] = "51716"
  
  # [36,] "ENSG00000228696" "100506084_&_100996709_&_51326"                         
  id = intersect(which(ensembl_list == "ENSG00000228696"), which(entrezID == "100506084_&_100996709_&_51326")) 
  entrezID[id] = "100506084"
  
  # [37,] "ENSG00000229571" "441873_&_645359"                                       
  id = intersect(which(ensembl_list == "ENSG00000229571"), which(entrezID == "441873_&_645359")) 
  entrezID[id] = "441873"
  
  # [38,] "ENSG00000230124" "100527964_&_84320"                                     
  id = intersect(which(ensembl_list == "ENSG00000230124"), which(entrezID == "100527964_&_84320")) 
  entrezID[id] = "84320"
  
  # [39,] "ENSG00000231431" "100420005_&_440910"                                    
  id = intersect(which(ensembl_list == "ENSG00000231431"), which(entrezID == "100420005_&_440910")) 
  entrezID[id] = "100420005"
  
  # [40,] "ENSG00000236125" "392188_&_645402"                                       
  id = intersect(which(ensembl_list == "ENSG00000236125"), which(entrezID == "392188_&_645402")) 
  entrezID[id] = "645402"
  
  # [41,] "ENSG00000236362" "100008586_&_2574_&_2576_&_2578_&_2579_&_26748_&_729408"
  id = intersect(which(ensembl_list == "ENSG00000236362"), which(entrezID == "100008586_&_2574_&_2576_&_2578_&_2579_&_26748_&_729408")) 
  entrezID[id] = "100008586"
  
  # [42,] "ENSG00000239533" "401634_&_84559"                                        
  id = intersect(which(ensembl_list == "ENSG00000239533"), which(entrezID == "401634_&_84559")) 
  entrezID[id] = "84559"
  
  # [43,] "ENSG00000241370" "202658_&_79897"                                        
  id = intersect(which(ensembl_list == "ENSG00000241370"), which(entrezID == "202658_&_79897")) 
  entrezID[id] = "79897"
  
  # [44,] "ENSG00000244355" "110599563_&_58530"                                     
  id = intersect(which(ensembl_list == "ENSG00000244355"), which(entrezID == "110599563_&_58530")) 
  entrezID[id] = "58530"
  
  # [45,] "ENSG00000244731" "110384692_&_720"                                       
  id = intersect(which(ensembl_list == "ENSG00000244731"), which(entrezID == "110384692_&_720")) 
  entrezID[id] = "720"
  
  # [46,] "ENSG00000248472" "100287596_&_100288486_&_727856"                        
  id = intersect(which(ensembl_list == "ENSG00000248472"), which(entrezID == "100287596_&_100288486_&_727856")) 
  entrezID[id] = "100288486"
  
  # [47,] "ENSG00000248767" "100288524_&_116033993"                                 
  id = intersect(which(ensembl_list == "ENSG00000248767"), which(entrezID == "100288524_&_116033993")) 
  entrezID[id] = "116033993"
  
  # [48,] "ENSG00000255251" "100131608_&_100133251"                                 
  id = intersect(which(ensembl_list == "ENSG00000255251"), which(entrezID == "100131608_&_100133251")) 
  entrezID[id] = "100131608"
  
  # [49,] "ENSG00000255374" "259289_&_259291"                                       
  id = intersect(which(ensembl_list == "ENSG00000255374"), which(entrezID == "259289_&_259291")) 
  entrezID[id] = "259289"
  
  # [50,] "ENSG00000257046" "115072896_&_338821"                                    
  id = intersect(which(ensembl_list == "ENSG00000257046"), which(entrezID == "115072896_&_338821")) 
  entrezID[id] = "115072896"
  
  # [51,] "ENSG00000257365" "100529261_&_2342"                                      
  id = intersect(which(ensembl_list == "ENSG00000257365"), which(entrezID == "100529261_&_2342")) 
  entrezID[id] = "2342"
  
  # [52,] "ENSG00000259295" "440300_&_728121"                                       
  id = intersect(which(ensembl_list == "ENSG00000259295"), which(entrezID == "440300_&_728121")) 
  entrezID[id] = "728121"
  
  # [53,] "ENSG00000260128" "100288380_&_89838"                                     
  id = intersect(which(ensembl_list == "ENSG00000260128"), which(entrezID == "100288380_&_89838")) 
  entrezID[id] = "100288380"
  
  # [54,] "ENSG00000261247" "653075_&_653125"                                       
  id = intersect(which(ensembl_list == "ENSG00000261247"), which(entrezID == "653075_&_653125")) 
  entrezID[id] = "653075"
  
  # [55,] "ENSG00000261509" "102723655_&_102723713_&_102724101_&_102724127_&_729355"
  id = intersect(which(ensembl_list == "ENSG00000261509"), which(entrezID == "102723655_&_102723713_&_102724101_&_102724127_&_729355")) 
  entrezID[id] = "729355"
  
  # [56,] "ENSG00000261556" "100506060_&_107984138"                                 
  id = intersect(which(ensembl_list == "ENSG00000261556"), which(entrezID == "100506060_&_107984138")) 
  entrezID[id] = "100506060"
  
  # [57,] "ENSG00000262628" "653166_&_8386"                                         
  id = intersect(which(ensembl_list == "ENSG00000262628"), which(entrezID == "653166_&_8386")) 
  entrezID[id] = "8386"
  
  # [58,] "ENSG00000263715" "104909134_&_1394"                                      
  id = intersect(which(ensembl_list == "ENSG00000263715"), which(entrezID == "104909134_&_1394")) 
  entrezID[id] = "104909134"
  
  # [59,] "ENSG00000269028" "100462981_&_100463498"                                 
  id = intersect(which(ensembl_list == "ENSG00000269028"), which(entrezID == "100462981_&_100463498")) 
  entrezID[id] = "100463498"
  
  # [60,] "ENSG00000269136" "100129976_&_641367"                                    
  id = intersect(which(ensembl_list == "ENSG00000269136"), which(entrezID == "100129976_&_641367")) 
  entrezID[id] = "100129976"
  
  # [61,] "ENSG00000273513" "101060351_&_102723859"                                 
  id = intersect(which(ensembl_list == "ENSG00000273513"), which(entrezID == "101060351_&_102723859")) 
  entrezID[id] = "101060351"
  
  # [62,] "ENSG00000274512" "101060376_&_729873"                                    
  id = intersect(which(ensembl_list == "ENSG00000274512"), which(entrezID == "101060376_&_729873")) 
  entrezID[id] = "101060376"
  
  # [63,] "ENSG00000274808" "414059_&_414060"                                       
  id = intersect(which(ensembl_list == "ENSG00000274808"), which(entrezID == "414059_&_414060")) 
  entrezID[id] = "414059"
  
  # [64,] "ENSG00000276070" "388372_&_9560"                                         
  id = intersect(which(ensembl_list == "ENSG00000276070"), which(entrezID == "388372_&_9560")) 
  entrezID[id] = "9560"
  
  # [65,] "ENSG00000276085" "414062_&_6349"                                         
  id = intersect(which(ensembl_list == "ENSG00000276085"), which(entrezID == "414062_&_6349")) 
  entrezID[id] = "414062"
  
  # [66,] "ENSG00000277147" "388692_&_57234"                                        
  id = intersect(which(ensembl_list == "ENSG00000277147"), which(entrezID == "388692_&_57234")) 
  entrezID[id] = "57234"
  
  # [67,] "ENSG00000278662" "643707_&_647042"                                       
  id = intersect(which(ensembl_list == "ENSG00000278662"), which(entrezID == "643707_&_647042")) 
  entrezID[id] = "647042"
  
  # [68,] "ENSG00000283154" "100505385_&_29970"                                     
  id = intersect(which(ensembl_list == "ENSG00000283154"), which(entrezID == "100505385_&_29970")) 
  entrezID[id] = "100505385"
  
  # [69,] "ENSG00000283196" "107985911_&_645166_&_654342"
  id = intersect(which(ensembl_list == "ENSG00000283196"), which(entrezID == "107985911_&_645166_&_654342")) 
  entrezID[id] = "645166"
  
  ### Priyanka's modifications
  # [70,] "ENSG00000236172" "102800314_&_102800315" 
  id = intersect(which(ensembl_list == "ENSG00000236172"), which(entrezID == "102800314_&_102800315")) 
  entrezID[id] = "102800315"
  
  # [71,] "ENSG00000267452" "105371828_&_440446"  
  id = intersect(which(ensembl_list == "ENSG00000267452"), which(entrezID == "105371828_&_440446")) 
  entrezID[id] = "440446"
  
  # [72,] "ENSG00000214401" "107984142_&_644246"  
  id = intersect(which(ensembl_list == "ENSG00000214401"), which(entrezID == "107984142_&_644246")) 
  entrezID[id] = "644246"
  
  # [73,] "ENSG00000231304" "100874028_&_101927829" 
  id = intersect(which(ensembl_list == "ENSG00000231304"), which(entrezID == "100874028_&_101927829")) 
  entrezID[id] = "100874028"
  
  # [74,] "ENSG00000229425" "101927745_&_105369302" 
  id = intersect(which(ensembl_list == "ENSG00000229425"), which(entrezID == "101927745_&_105369302")) 
  entrezID[id] = "101927745"
  
  # [75,] "ENSG00000225914" "101929163_&_414764" 
  id = intersect(which(ensembl_list == "ENSG00000225914"), which(entrezID == "101929163_&_414764")) 
  entrezID[id] = "414764"
  
  # [76,] "ENSG00000214900" "196913_&_283551"
  id = intersect(which(ensembl_list == "ENSG00000214900"), which(entrezID == "196913_&_283551")) 
  entrezID[id] = "196913"
  
  # [77,] "ENSG00000259905" "101928840_&_791114"
  id = intersect(which(ensembl_list == "ENSG00000259905"), which(entrezID == "101928840_&_791114")) 
  entrezID[id] = "791114"
  
  # [78,] "ENSG00000151067" "100874369_&_775"  
  id = intersect(which(ensembl_list == "ENSG00000151067"), which(entrezID == "100874369_&_775")) 
  entrezID[id] = "775"
  
  # [79,] "ENSG00000224533" "100507404_&_101927830"  
  id = intersect(which(ensembl_list == "ENSG00000224533"), which(entrezID == "100507404_&_101927830")) 
  entrezID[id] = "100507404"
  
  # [80,] "ENSG00000145391" "105377622_&_80854" 
  id = intersect(which(ensembl_list == "ENSG00000145391"), which(entrezID == "105377622_&_80854")) 
  entrezID[id] = "80854"
  
  # [81,] "ENSG00000245750" "104472713_&_145837_&_414926" 
  id = intersect(which(ensembl_list == "ENSG00000245750"), which(entrezID == "104472713_&_145837_&_414926")) 
  entrezID[id] = "414926"
  
  # [82,] "ENSG00000226833" "100505774_&_112267877"
  id = intersect(which(ensembl_list == "ENSG00000226833"), which(entrezID == "100505774_&_112267877")) 
  entrezID[id] = "112267877"
  
  # [83,] "ENSG00000255760" "105369723_&_105369725"  
  id = intersect(which(ensembl_list == "ENSG00000255760"), which(entrezID == "105369723_&_105369725")) 
  entrezID[id] = "105369723"
  
  # [84,] "ENSG00000277610" "101954264_&_101954268"  
  id = intersect(which(ensembl_list == "ENSG00000277610"), which(entrezID == "101954264_&_101954268")) 
  entrezID[id] = "101954264"
  
  # [85,] "ENSG00000224897" "101928283_&_401398"  
  id = intersect(which(ensembl_list == "ENSG00000224897"), which(entrezID == "101928283_&_401398")) 
  entrezID[id] = "401398"
  
  # [86,] "ENSG00000228262" "104355287_&_104355288" 
  id = intersect(which(ensembl_list == "ENSG00000228262"), which(entrezID == "104355287_&_104355288")) 
  entrezID[id] = "104355287"
  
  # [87,] "ENSG00000215483" "400123_&_646982" 
  id = intersect(which(ensembl_list == "ENSG00000215483"), which(entrezID == "400123_&_646982")) 
  entrezID[id] = "646982"
  
  # [88,] "ENSG00000266256" "284276_&_400660"  
  id = intersect(which(ensembl_list == "ENSG00000266256"), which(entrezID == "284276_&_400660")) 
  entrezID[id] = "284276"
  
  # [89,] "ENSG00000177910" "441452_&_645961"  
  id = intersect(which(ensembl_list == "ENSG00000177910"), which(entrezID == "441452_&_645961")) 
  entrezID[id] = "645961"
  
  # [90,] "ENSG00000227036" "100499467_&_400619"  
  id = intersect(which(ensembl_list == "ENSG00000227036"), which(entrezID == "100499467_&_400619")) 
  entrezID[id] = "400619"
  
  # [91,] "ENSG00000199337" "100169754_&_100169763"
  id = intersect(which(ensembl_list == "ENSG00000199337"), which(entrezID == "100169754_&_100169763")) 
  entrezID[id] = "100169754"
  
  # [92,] "ENSG00000275239" "105376064_&_105379447"
  id = intersect(which(ensembl_list == "ENSG00000275239"), which(entrezID == "105376064_&_105379447")) 
  entrezID[id] = "105376064"
  
  # [93,] "ENSG00000184154" "120356739_&_220074" 
  id = intersect(which(ensembl_list == "ENSG00000184154"), which(entrezID == "120356739_&_220074")) 
  entrezID[id] = "120356739"
  
  # [94,] "ENSG00000231131" "101928687_&_105378305" 
  id = intersect(which(ensembl_list == "ENSG00000231131"), which(entrezID == "101928687_&_105378305")) 
  entrezID[id] = "101928687"
  
  # [95,] ""ENSG00000279561" "100132249_&_105376060 
  id = intersect(which(ensembl_list == "ENSG00000279561"), which(entrezID == "100132249_&_105376060")) 
  entrezID[id] = "100132249"
  
  # [96,] "ENSG00000236829" "100134368_&_729480"
  id = intersect(which(ensembl_list == "ENSG00000236829"), which(entrezID == "100134368_&_729480")) 
  entrezID[id] = "100134368"
  
  # [97,] "ENSG00000232527" "100996732_&_107985200"  
  id = intersect(which(ensembl_list == "ENSG00000232527"), which(entrezID == "100996732_&_107985200")) 
  entrezID[id] = "100996732"
  
  # [98,] "ENSG00000184319" "118433_&_284942"   
  id = intersect(which(ensembl_list == "ENSG00000184319"), which(entrezID == "118433_&_284942")) 
  entrezID[id] = "284942"
  
  # [99,] "ENSG00000237298" "100506866_&_101927055"   
  id = intersect(which(ensembl_list == "ENSG00000237298"), which(entrezID == "100506866_&_101927055")) 
  entrezID[id] = "100506866"
  
  # [100,] "ENSG00000267322" "103091864_&_677769"  
  id = intersect(which(ensembl_list == "ENSG00000267322"), which(entrezID == "103091864_&_677769")) 
  entrezID[id] = "103091864"
  
  # [101,] "ENSG00000248489" "100289230_&_102724810"   
  id = intersect(which(ensembl_list == "ENSG00000248489"), which(entrezID == "100289230_&_102724810")) 
  entrezID[id] = "102724810"
  
  # [102,] "ENSG00000236107" "101929680_&_102724058"   
  id = intersect(which(ensembl_list == "ENSG00000236107"), which(entrezID == "101929680_&_102724058")) 
  entrezID[id] = "101929680"
  
  # [103,] "ENSG00000233614" "100287029_&_727856"   
  id = intersect(which(ensembl_list == "ENSG00000233614"), which(entrezID == "100287029_&_727856")) 
  entrezID[id] = "100287029"
  
  # [104,] "ENSG00000179818" "151516_&_400960"    
  id = intersect(which(ensembl_list == "ENSG00000179818"), which(entrezID == "151516_&_400960")) 
  entrezID[id] = "400960"
  
  # [105,] "ENSG00000270011" "100529215_&_7730"    
  id = intersect(which(ensembl_list == "ENSG00000270011"), which(entrezID == "100529215_&_7730")) 
  entrezID[id] = "100529215"
  
  # [106,] "ENSG00000236663" "106481742_&_106780825"    
  id = intersect(which(ensembl_list == "ENSG00000236663"), which(entrezID == "106481742_&_106780825")) 
  entrezID[id] = "106481742"
  
  # [107,] "ENSG00000230417" "100132987_&_414243"   
  id = intersect(which(ensembl_list == "ENSG00000230417"), which(entrezID == "100132987_&_414243")) 
  entrezID[id] = "100132987"
  
  # [108,] "ENSG00000251655" "5542_&_653247"  
  id = intersect(which(ensembl_list == "ENSG00000251655"), which(entrezID == "5542_&_653247")) 
  entrezID[id] = "5542"
  
  # [109,] "ENSG00000236362" "100008586_&_2574_&_2576_&_2578_&_2579_&_26748"  
  id = intersect(which(ensembl_list == "ENSG00000236362"), which(entrezID == "100008586_&_2574_&_2576_&_2578_&_2579_&_26748")) 
  entrezID[id] = "100008586"
  
  # [110,] "ENSG00000245694" "101927480_&_643911"   
  id = intersect(which(ensembl_list == "ENSG00000245694"), which(entrezID == "101927480_&_643911")) 
  entrezID[id] = "643911"
  
  # [111,] "ENSG00000237159" "105376023_&_415056"   
  id = intersect(which(ensembl_list == "ENSG00000237159"), which(entrezID == "105376023_&_415056")) 
  entrezID[id] = "415056"
  
  # [112,] "ENSG00000223972" "100287102_&_100287596_&_102725121_&_727856_&_84771"   
  id = intersect(which(ensembl_list == "ENSG00000223972"), which(entrezID == "100287102_&_100287596_&_102725121_&_727856_&_84771")) 
  entrezID[id] = "100287102"
  
  # [113,] "ENSG00000226808" "100506835_&_414260"    
  id = intersect(which(ensembl_list == "ENSG00000226808"), which(entrezID == "100506835_&_414260")) 
  entrezID[id] = "414260"
  
  # [114,] "ENSG00000183625" "105377067_&_1232"     
  id = intersect(which(ensembl_list == "ENSG00000183625"), which(entrezID == "105377067_&_1232")) 
  entrezID[id] = "1232"
  
  # [115,] "ENSG00000248538" "101929128_&_157273"      
  id = intersect(which(ensembl_list == "ENSG00000248538"), which(entrezID == "101929128_&_157273")) 
  entrezID[id] = "102724880"
  
  # [116,] "ENSG00000242153" "386690_&_83864"      
  id = intersect(which(ensembl_list == "ENSG00000242153"), which(entrezID == "386690_&_83864")) 
  entrezID[id] = "386690"
  
  # [117,] "ENSG00000229140" "106144608_&_137196_&_728724"        
  id = intersect(which(ensembl_list == "ENSG00000229140"), which(entrezID == "106144608_&_137196_&_728724")) 
  entrezID[id] = "137196"
  
  # [118,] "ENSG00000274808" "101060321_&_414059"         
  id = intersect(which(ensembl_list == "ENSG00000274808"), which(entrezID == "101060321_&_414059")) 
  entrezID[id] = "414059"
  
  # [119,] "ENSG00000198948" "101928198_&_9848"          
  id = intersect(which(ensembl_list == "ENSG00000198948"), which(entrezID == "101928198_&_9848")) 
  entrezID[id] = "9848"
  
  # [120,] "ENSG00000248461" "101929745_&_105374731"           
  id = intersect(which(ensembl_list == "ENSG00000248461"), which(entrezID == "101929745_&_105374731")) 
  entrezID[id] = "101929745"
  
  # [121,] "ENSG00000213904" "100996307_&_101930071"           
  id = intersect(which(ensembl_list == "ENSG00000213904"), which(entrezID == "100996307_&_101930071")) 
  entrezID[id] = "100996307"
  
  # [122,] "ENSG00000231527" "100132948_&_105379444"          
  id = intersect(which(ensembl_list == "ENSG00000231527"), which(entrezID == "100132948_&_105379444")) 
  entrezID[id] = "100132948"
  
  # [123,] "ENSG00000286105" "121725057_&_221468"           
  id = intersect(which(ensembl_list == "ENSG00000286105"), which(entrezID == "121725057_&_221468")) 
  entrezID[id] = "121725057"
  
  # [124,] "ENSG00000130600" "102724852_&_283120"           
  id = intersect(which(ensembl_list == "ENSG00000130600"), which(entrezID == "102724852_&_283120")) 
  entrezID[id] = "283120"
  
  # [125,] "ENSG00000233098" "339260_&_440416"           
  id = intersect(which(ensembl_list == "ENSG00000233098"), which(entrezID == "339260_&_440416")) 
  entrezID[id] = "440416"
  
  # [126,] "ENSG00000215378" "100287083_&_170949"          
  id = intersect(which(ensembl_list == "ENSG00000215378"), which(entrezID == "100287083_&_170949")) 
  entrezID[id] = "170949"
  
  # [127,] "ENSG00000143384" "4170_&_574406"           
  id = intersect(which(ensembl_list == "ENSG00000143384"), which(entrezID == "4170_&_574406")) 
  entrezID[id] = "4170"
  
  # [128,] "ENSG00000196242" "148824_&_81472"           
  id = intersect(which(ensembl_list == "ENSG00000196242"), which(entrezID == "148824_&_81472")) 
  entrezID[id] = "81472"
  
  # [129,] "ENSG00000171714" "102723370_&_203859"            
  id = intersect(which(ensembl_list == "ENSG00000171714"), which(entrezID == "102723370_&_203859")) 
  entrezID[id] = "203859"
  
  # [130,] "ENSG00000239149" "677882_&_677885"            
  id = intersect(which(ensembl_list == "ENSG00000239149"), which(entrezID == "677882_&_677885")) 
  entrezID[id] = "677885"
  
  # [131,] "ENSG00000237541" "3117_&_3118"            
  id = intersect(which(ensembl_list == "ENSG00000237541"), which(entrezID == "3117_&_3118")) 
  entrezID[id] = "3118"
  
  # [132,] "ENSG00000189064" "100101629_&_26749_&_645037_&_729408_&_729447"          
  id = intersect(which(ensembl_list == "ENSG00000189064"), which(entrezID == "100101629_&_26749_&_645037_&_729408_&_729447")) 
  entrezID[id] = "729447"
  
  # [133,] "ENSG00000140650" "100130283_&_5373"         
  id = intersect(which(ensembl_list == "ENSG00000140650"), which(entrezID == "100130283_&_5373")) 
  entrezID[id] = "5373"
  
  # [134,] "ENSG00000228340" "284757_&_729296"        
  id = intersect(which(ensembl_list == "ENSG00000228340"), which(entrezID == "284757_&_729296")) 
  entrezID[id] = "284757"
  
  # [135,] "ENSG00000275496" "102723360_&_102724701"       
  id = intersect(which(ensembl_list == "ENSG00000275496"), which(entrezID == "102723360_&_102724701")) 
  entrezID[id] = "102723360"
  
  # [136,] "ENSG00000227110" "100288428_&_101927394"      
  id = intersect(which(ensembl_list == "ENSG00000227110"), which(entrezID == "100288428_&_101927394")) 
  entrezID[id] = "100288428"
  
  # [137,] "ENSG00000249981" "107987420_&_107987434"      
  id = intersect(which(ensembl_list == "ENSG00000249981"), which(entrezID == "107987420_&_107987434")) 
  entrezID[id] = "107987420"
  
  # [138,] "ENSG00000109927" "116804918_&_7007"       
  id = intersect(which(ensembl_list == "ENSG00000109927"), which(entrezID == "116804918_&_7007")) 
  entrezID[id] = "7007"
  
  # [139,] "ENSG00000233008" "101927560_&_101927587"      
  id = intersect(which(ensembl_list == "ENSG00000233008"), which(entrezID == "101927560_&_101927587")) 
  entrezID[id] = "101927587"
  
  # [140,] "ENSG00000242147" "102723629_&_105376382"      
  id = intersect(which(ensembl_list == "ENSG00000242147"), which(entrezID == "102723629_&_105376382")) 
  entrezID[id] = "105376382"
  
  # [141,] "ENSG00000269226" "122394733_&_286527"      
  id = intersect(which(ensembl_list == "ENSG00000269226"), which(entrezID == "122394733_&_286527")) 
  entrezID[id] = "122394733"
  
  # [142,] "ENSG00000250302" "107986281_&_152578"      
  id = intersect(which(ensembl_list == "ENSG00000250302"), which(entrezID == "107986281_&_152578")) 
  entrezID[id] = "152578"
  
  # [143,] "ENSG00000245293" "101929595_&_107986298"       
  id = intersect(which(ensembl_list == "ENSG00000245293"), which(entrezID == "101929595_&_107986298")) 
  entrezID[id] = "101929595"
  
  # [144,] "ENSG00000143429" "645166_&_654342"       
  id = intersect(which(ensembl_list == "ENSG00000143429"), which(entrezID == "645166_&_654342")) 
  entrezID[id] = "654342"
  
  # [145,] "ENSG00000254911" "100158262_&_619383"       
  id = intersect(which(ensembl_list == "ENSG00000254911"), which(entrezID == "100158262_&_619383")) 
  entrezID[id] = "619383"
  
  # [146,] "ENSG00000119812" "105374454_&_25940"      
  id = intersect(which(ensembl_list == "ENSG00000119812"), which(entrezID == "105374454_&_25940")) 
  entrezID[id] = "25940"
  
  # [147,] "ENSG00000281880" "103157000_&_440034"      
  id = intersect(which(ensembl_list == "ENSG00000281880"), which(entrezID == "103157000_&_440034")) 
  entrezID[id] = "103157000"
  
  # [148,] "ENSG00000249642" "105377603_&_107983963"       
  id = intersect(which(ensembl_list == "ENSG00000249642"), which(entrezID == "105377603_&_107983963")) 
  entrezID[id] = "105377603"
  entrezID
}

# Function to convert entrz id to ensembl ID 
eg2ensembl <- function(eg)
{
  library(org.Hs.eg.db)
  
  x <- org.Hs.egENSEMBL
  mapped_genes <- mappedkeys(x)
  ori_eg2ensembl <- as.list(x[mapped_genes])
  eg2ensembl = rep("", length = length(ori_eg2ensembl))
  names(eg2ensembl) = names(ori_eg2ensembl)
  for (i in 1:length(eg2ensembl))
  {
    if (i == 1)
    {
      eg2ensembl[1] = ori_eg2ensembl[[1]][1]
    } else {
      idi = which(!(ori_eg2ensembl[[i]] %in% eg2ensembl[1:(i - 1)]))
      if (length(idi) == 0)
      {
        eg2ensembl[i] = ori_eg2ensembl[[i]][1]
      } else {
        eg2ensembl[i] = ori_eg2ensembl[[i]][idi[1]]
      }
    }
  }
  
  ensembl = rep("", length(eg))
  for (i in seq(1, length(eg)))
  {
    if ((i %% ceiling(length(eg)/20) == 0) & (length(eg) > 100000))
    {
      print(paste(round(i/length(eg), 2)*100, " percents finished.", sep = ""))
    }
    
    if (length(grep("_&_", eg[i])) > 0)
    {
      eg_s = strsplit(eg[i], split = "_&_")[[1]]
      eg_s = eg_s[which(eg_s %in% names(eg2ensembl))]
      if (length(eg_s) > 0)
      {
        ensembl[i] = paste(as.character(eg2ensembl[eg_s]), collapse = "_&_")  
      }
    } else if (eg[i] %in% names(eg2ensembl))
    {
      ensembl[i] = eg2ensembl[as.character(eg[i])]
    }
  }
  #Manual changes
  # [1,]       
  id = intersect(which(ensembl == "ENSG00000179818"), which(eg == "151516")) 
  ensembl[id] = "ENSG00000244617"
  
  # [2,]       
  id = intersect(which(ensembl == "ENSG00000105887"), which(eg == "767558")) 
  ensembl[id] = ''
  
  # [3,]       
  id = intersect(which(ensembl == "ENSG00000105889"), which(eg == "401312")) 
  ensembl[id] = ''
  
  # [4,]       
  id = intersect(which(ensembl == "ENSG00000145979"), which(eg == "107080638")) 
  ensembl[id] = ''
  
  # [5,]       
  id = intersect(which(ensembl == "ENSG00000213160"), which(eg == "100526832")) 
  ensembl[id] = ''
  
  # [6,]       
  id = intersect(which(ensembl == "ENSG00000213999"), which(eg == "4207")) 
  ensembl[id] = 'ENSG00000064489'
  
  # [7,]       
  id = intersect(which(ensembl == "ENSG00000239149"), which(eg == "677882")) 
  ensembl[id] = 'ENSG00000266079'
  
  # [8,]       
  id = intersect(which(ensembl == "ENSG00000288170"), which(eg == "100287596")) 
  ensembl[id] = 'ENSG00000236875'
  
  # [9,]       
  id = intersect(which(ensembl == "ENSG00000288170"), which(eg == "727856")) 
  ensembl[id] = 'ENSG00000227159'
  
  # [10,]       
  id = intersect(which(ensembl == "ENSG00000158747"), which(eg == "100532736")) 
  ensembl[id] = 'ENSG00000270136'
  
  ensembl
}


### Mutation file
mut_file = read.csv(file = '../CCLE_Multiomics_Data/CCLE_mutations.csv')
gs <- mut_file[which(mut_file['Entrez_Gene_Id']!=0), 1]
gs <- gs[! gs %in% 'LST3'] # Remove LST3 as, 28234-SLCO1B3, 84173-ELMOD3, 338821-SLCO1B7 are suitable mappings.These three genes are already mapped
#Convert entrenz id from mutation file to gene symbol
en_id <- mut_file[which(mut_file['Entrez_Gene_Id']!=0), 2]
gs_new2 <- eg2gs(en_id)
ens_id = eg2ensembl(en_id)
#map1 <- data.frame(orig_GS = gs, entrez_id = en_id_new, new_GS = gs_new)
map_new <- data.frame(orig_GS = gs,old_entrez_id = en_id, en2gs_GS = unlist(gs_new2), ens_id = unlist(ens_id))
write.csv(map_new, '../Data_Curation_final/Maps/eg_gs_map_mut.csv')

### Gene expression and gene copy number
es_gs_file = read.csv(file = '../Maps/es_gs_orig.csv')
es<-es_gs_file['entrenz_id']
gs_org <- es_gs_file['gene_symbol']
g_s <- eg2gs(unlist(es))
ens <- eg2ensembl(unlist(es))
map_new <- data.frame(orig_GS = unlist(gs_org),old_entrez_id = unlist(es), en2gs_GS = unlist(g_s), ens_id = unlist(ens))
write.csv(map_new, '../Maps/eg_gs_map_cn_ge.csv')

### CCLE RRBS file 
rrbs_file <- read.delim('../CCLE_Multiomics_Data/CCLE_RRBS_TSS_1kb_20180614.txt')
gs_orig <- rrbs_file['gene']
e_id <- gs2eg(unlist(gs_orig))
# Convert back the entrenz id to gene symbol
gs_new <- eg2gs(e_id)
ens_id = eg2ensembl(e_id)
map_rrbs <- data.frame(orig_GS = unlist(gs_orig), new_entrenz_id = unlist(e_id), new_GS = unlist(gs_new), ens_id = unlist(ens_id))
write.csv(map_rrbs, '../Maps/gs_eg_gs_map_rrbs.csv')

#CCLE_gene_expression_full
gs_full_file = read.csv(file = '../Maps/ens_gs_gef_orig.csv')
ens<-gs_full_file['ensembl_id']
gs_org <- gs_full_file['gene_symbol']
e_s <- ensembl2eg(unlist(ens)) # Convert Ensembl ID to Entrenz ID
gs_new <- eg2gs(e_s) # Convert Entrenz ID to Gene Symbol
map_new <- data.frame(orig_GS = unlist(gs_org),orig_ensembl_id = unlist(ens), entrenzID = unlist(as.character(e_s)), new_GS = unlist(gs_new))
write.csv(map_new, '../Maps/ens_eg_gs_map_gef.csv')

