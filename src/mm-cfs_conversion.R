#conversion factors
#area_ft2_ORO = area_mi2_ORO*5280^2 # 1 mile = 5280ft
#fnf_mm_day_ORO = fnf_cfs/area_ft2_ORO*3600*24*304.8 # 1ft = 304.8mm 
mm_to_cfs = 5280^2/3600/24/304.8

area_mi2_MHB = 536
area_mi2_FOL = 1885
area_mi2_MIL = 1675
area_mi2_GDW = 986
area_mi2_ISB = 2074
area_mi2_MRC = 1061
area_mi2_MKM = 544
area_mi2_NHG = 363
area_mi2_NML = 900
area_mi2_ORO = 3607
area_mi2_PNF = 1545
area_mi2_BND = 8900
area_mi2_SHA = 6665
area_mi2_SCC = 393
area_mi2_TRM = 561
area_mi2_CLE = 688
area_mi2_TLG = 1538
area_mi2_WHI = 201
area_mi2_YRS = 1108

area_mi2 = c(area_mi2_SHA,area_mi2_ORO,area_mi2_YRS,area_mi2_FOL,area_mi2_MKM,area_mi2_NHG,area_mi2_NML,area_mi2_TLG,area_mi2_MRC,
                 area_mi2_MIL,area_mi2_PNF,area_mi2_TRM,area_mi2_SCC,area_mi2_ISB)
area_lab = c('SHA', 'ORO','YRS','FOL','MKM','NHG','NML','TLG','MRC','MIL','PNF','TRM','SCC','ISB')
area_calc = cbind(area_lab,area_mi2)
rownames(area_calc)<-area_lab
