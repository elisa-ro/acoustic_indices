## Acoustic Indices

This folder is dedicated to calculated acoustic indices from .wav files. Many of the files are used for exploratory purposes (converting audio files into .wav format and coverting from stereo to mono audio files). 

### To calculate bioacoustic indicators:
Given a folder of .wav files who's bioacoustic indices need to be calculated, the Report.ipynb should be run. The folder path should be stated at the beginning of the file.
Lastly run all of the cells.
The csv file with all of the bioacoustic indicators for each file will be found in the path folder.

Layman instructions for a beginner:
1. Download the Report.ipynb file that I've attached (it's the same one on my github with modifications if you don't have Jupyter Notebooks)
2. If one doesn't have an app to read the .ipynb files, I'd recommend downloading VSC (https://code.visualstudio.com/)
3. Open then Report.ipynb file on VSC. If  there are pop-ups because you don't have python installed, download it
4. In your computer, you should have a folder with audio files, copy the path to that folder and replace the variable directory = "<your audio file folder>" instead of  "C:\\Users\\elisarog\\Desktop\\arbimon_labels"
5. For each 'cell' in the Report ipynb file, press play or 'run it' 
6. Once you've run through all of them, in your folder with the audio files, you should find a csv file called "bioacoustic indicators" that has the ecoindicators
