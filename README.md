# Analysis code for Keung, et al. 2019 Regulation of evidence accumulation by pupil-linked processes

`SuppSoftware_pupil.m` runs all the analyses described in [Keung, et al. 2019 Nat. Hum. Beh.][Keung], and reproduces all the 
main figures in the paper. Running each cell consecutively will generate each analysis/figure in the paper.

`regs_test.mat` contains test data for two participants. Data are organized into a n element struct array 
(n is the number of participants). Full dataset can be found at https://osf.io/37yk8/

```
Each element has the following fields:
      clicksLR:         click stimulus,                         # clicks x # trials matrix,  coded as +1 for left click and -1 for right click
      choiceLR:         participant's choices,                  # trials vector,             coded as +1 for left and -1 right 
      dN_all:           net difference in left vs right clicks, # trials vector
      RT:               reaction time in seconds,               # trials vector
      clicks_on_pd_z.D: z-scored pupil diameter data (30 Hz), centered on stimulus onset, each data point is averaged across both left and right eyes.
```

[Keung]:http://nature.com/articles/s41562-019-0551-4
