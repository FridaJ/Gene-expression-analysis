##### qPCR assignment for gene expr analysis course #####

setwd("Documents/R for gene expression analysis/")
install.packages("pcr")
install.packages("readr")
library(pcr)
library(readr)

fl <- system.file("extdata", "ct1.csv", package = "pcr") #filepath to example file in pcr package
ct1 <- readr::read_csv(fl) #read csv example file

group_var <- rep(c("brain", "kidney"), each = 6) #Grouping variable. First 6 are brain, the rest kidney

##### Investigating structure of data #####

#> ct1
## A tibble: 12 x 2
#c_myc GAPDH
#  <dbl> <dbl>
#1   30.7  23.7
#2   30.3  23.6
#3   30.6  23.5
#4   30.3  23.6
#5   30.5  23.7
#6   30.4  23.7
#7   27.1  22.8
#8   27.0  22.6
#9   27.0  22.6
#10  27.1  22.6
#11  27.0  22.6
#12  26.9  22.8

##### Task 1, calculate the relative expression between brain/kidneys #####

pcr_analyze(ct1,
            group_var = group_var,
            reference_gene = 'GAPDH',
            reference_group = 'brain',
            method = 'delta_delta_ct',
            plot = FALSE)

#   group  gene normalized calibrated relative_expression      error    lower    upper
#1  brain c_myc      6.860      0.000            1.000000 0.17395402 0.886410 1.128146
#2 kidney c_myc      4.365     -2.495            5.637283 0.09544632 5.276399 6.022850

#The relative expression of c_myc in kidney compared to brain is 5.64.

#The following can be found under "Details" of the pcr_ddct function:
#"The comparative C_T methods assume that the cDNA templates of the gene/s of interest as well as the 
#control/reference gene have similar amplification efficiency. And that this amplification efficiency is near 
#perfect. Meaning, at a certain threshold during the linear portion of the PCR reaction, the amount of the gene 
#of the interest and the control double each cycle. Another assumptions is that, the expression difference 
#between two genes or two samples can be captured by subtracting one (gene or sample of interest) from another 
#(reference). This final assumption requires also that these references don't change with the treatment or the 
#course in question."


##### Task 2, Calculating effiency #####

fl <- system.file('extdata', 'ct3.csv', package = 'pcr')
ct3 <- readr::read_csv(fl)

# make amount/dilution variable:
amount <- rep(c(1, .5, .2, .1, .05, .02, .01), each = 3)

# calculate the standard curve
pcr_assess(ct3,
           amount = amount,
           method = 'standard_curve',
           plot = FALSE)

#  gene  intercept slope r_squared
#1 c_myc  25.69669 -3.388095 0.9965504
#2 GAPDH  22.68221 -3.414551 0.9990278

#As you remember from the lectures the value of the amplification efficiency is given by E=10^(−1/slope) -1. 
#What is the slope and the efficiency for c_myc? Include the values in your report.

slope1 = -3.388
E1 = 10^(-1/slope1)-1 #0.973

res_lm <- lm(ct3$c_myc ~log10(amount)) #Check that slope is almost the same. (-3.388...)
summary(res_lm)

##### Task 3, plotting efficiency #####

pcr_assess(ct3,
           amount = amount,
           reference_gene = "GAPDH",
           method = "efficiency",
           plot = TRUE)

#   gene intercept      slope  r_squared
#1 c_myc   3.01448 0.02645619 0.02070273

slope2 = 0.0265 #Should be as close to 0 as possible (if slopes of target and reference are the same). <0.01 ok
#(E2 = 10^(-1/slope2) # efficiency is -1 (1.8e-38 - 1)

# For report: Look at the slope of the trend line. What does the slope indicate? 
# What are the values of the slope and R2? 


##### Task 4, Standard curve quantification #####

res_ct3 <- pcr_assess(ct3, amount = amount, method = 'standard_curve', plot = FALSE)
cmyc_col <- which(res_ct3$gene == "c_myc") # get index for c_myc
slope_ct3 <- res_ct3$slope[cmyc_col] #use index found for c_myc
intercept_ct3 <- res_ct3$intercept[cmyc_col]

#Redo task 1 with more accurate values on slope and intercept
pcr_analyze(ct1,
            group_var = group_var,
            reference_gene = 'GAPDH',
            reference_group = 'brain',
            intercept = intercept_ct3,
            slope = slope_ct3,
            method = 'relative_curve',
            plot = FALSE)

#   group  gene  normalized calibrated       error     lower    upper
#1  brain c_myc 0.009470115   1.000000 0.001108129 0.8829867 1.117013
#2 kidney c_myc 0.051454535   5.433359 0.003304818 5.3691310 5.497587

## Q. Does the result differ from Exercise 1? If so how?


##### Task 5, Testing for statistical significance #####

fl <- system.file('extdata', 'ct4.csv', package = 'pcr')
ct4 <- readr::read_csv(fl) #Dataframe 24x2
group <- rep(c("ctrl", "treatment"), each = 12) #Are the first 12 controls? I will assume so.

#Do t.test for difference between ref and target gene in ctrl/treatment groups
pcr_test(ct4, group_var = group, reference_gene = "ref", reference_group = "ctrl", test = "t.test")
#Do wilcox.test for difference between ref and target gene in ctrl/treatment groups
pcr_test(ct4, group_var = group, reference_gene = "ref", reference_group = "ctrl", test = "wilcox.test")

#           gene  estimate      p_value     lower     upper
#t-ctrl   target -0.684825 3.429195e-05 -0.955952 -0.413698
#w-ctrl   target  -0.6354  2.958409e-06 -0.8805   -0.4227
 
#Is reference_control really an existing argument?? I did not use it, instead reference gene/group. 

#The "estimate" is the calculated delta-delta C_T for the genes. Since it is a negative value, I assume that the 
#target gene is down_regulated, by a factor of 0.685 and 0.635 when using the t.test and wilcox.test, respectively.
#The p-values given show that the wilcox test gives a slightly higher significance, by an order of ten.

#The simple t-test can be used to test the significance of the difference between two conditions Δ C_T. 
#t-test assumes in addition, that the input C_T values are normally distributed and the variance between 
#conditions are comparable. Wilcoxon test can be used when sample size is small and those two last assumptions 
#are hard to achieve.





