# README #

This code is originally developed by Riad Akhundov (Riad.Akhundov@uon.edu.au) and adapted by Miel Willems and van der Have Arthur in order to fit the requirements of the Human Movement Biomechanics Research Group Leuven. 

This code automaticly classifies EMG (Noise, Unusable, Usable, Good) using a pretrained deep learning neural network (AlexNet), filters the EMG, stores the classification results in a excel file and in specific folders. At last it stores all the data in a 
matlab strucuture which can be loaded back in to further postproces. You can find the published paper explaining the classification methods in this folder. 

#How to run the code: 
1) click run. 
2) Select the folder where you stored all the c3d files of the participant you want to analyze
3) Define your low pass and bass pass filter in the dialog window. You could manually crop the window if you want by defining the start and end timestamp in a later pop up window.
4) Indictate which EMG system you used
5) Define which muscles you measured on which channel with the specific muscles

# You should end up with: (everything is stored witin the participants folder)
1) an excel file with the classification results
2) pictures of all signals in classified folders 
3) a matlab structure called EMG.mat containing all EMG data of that participant

