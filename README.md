# forecast-Long_Term_Projection_SARIMA
Long term projection of short term SARIMA estimate of raw liquidity need

1. cover1.csv data file, time series of raw liquidity need: quasi-stationary time series with periodic spikes of random magnitude
2. true_roberto.R R script
3. Liquidity_Model_WhitePaper.docx model White Paper
4. WhitePaper_Improvement_Links.txt links to further model improvement
5. OU_Calibration.py module to calibrate a 1-DIM Ornstein-Uhlenbeck process on a time series via MLE (see last TODO bullet)

TODO:
* Improve the model by computing the distribution of the Max 
* Extend to numerical + theoretical analysis of the periodic jumps of random magnitude
* Alternative approach based on calibrating 1-DIM OU on the time series directly
