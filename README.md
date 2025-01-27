# VWstudyPaper
A repository containing analysis code pertaining to my performance study of the CLS vibrating wire setup for magnetic measurements. This repository is directly related to the journal article in progress, "Performance Study of the Vibrating Wire Technique to Determine Longitudinal Magnetic Field Profile Using Scans to High Wire Harmonic."


--------------------
Date: 2025-01-08
Cameron Baribeau
ReadMe for shared GitRepo for vibrating wire performance study (2024)
--------------------

--------------------
Files in this directory:
- Code: data analysis macros to process vibrating wire scans written for Igor Pro v6.37; plain text Mathematica notebook to convert the B field of a Hall probe scan into DST Sine series coefficients; auxiliary settings file for vibrating wire macros


--------------------
File format of vibrating wire scans:
- (Important to understand the scan file format in parsing the code.)
- Comma separated variables in raw TXT file
- One file per frequency sweep
- Scan date listed at start of file name
- Wire drive frequency is the swept variable ("fscan")
- Wire harmonic listed at end of file name (..._nNN.text)
- Data columns listed below:
- (1) Temp: ambient room temperature
- (2) Frequency: frequency set at signal generator
- (3) FCX: position of horizontal (X) motorized stage relative to wire (or coil) origin
- (4) FCY: position of vertical (Y) motorized stage relative to wire (or coil) origin
- (5) MeanX1: Mean from (N) recorded samples, X channel, lock-in amplifier #1
- (6) MeanY1: Mean from (N) recorded samples, Y channel, lock-in amplifier #1
- (7) MeanMag1: Mean from (N) recorded samples, magnitude of X+Y channels, amp #1
- (8) MeanPhase1: Mean from (N) recorded samples, phase of X+Y channels, amp #1
- (9) stdX1: Standard deviation from (N) recorded samples, X channel, amp #1
- (10) stdY1: Standard deviation from (N) recorded samples, Y channel, amp #1
- (11) stdMag1: Standard deviation from (N) recorded samples, magnitude of X+Y, amp #1
- (12) stdPhase1: Standard deviation from (N) recorded samples, phase of X+Y, amp #1
- (13...20) Same data types as columns 5...12 now for lock-in amplifier #2

--------------------
Vibrating wire scans, additional notes:
- Igor code is written to automatically detect the wire drive current IF it is listed in the file name
- We varied the wire drive current and did not fit it into the file name
- The wire drive current must be set manually in the Igor wave "vwScanCurrents" after import of scan files; one can refer to "vwScanCurrents" in the ProcessedData to see the appropriate wire currents (in Amps) for each scan 


--------------------
File format of Hall probe scans:
- Comma separated variables in raw text file
- ".hp" file extension: hall probe; a text editor should be able to view the file
- ".hpd" file extension: hall probe data; a text editor should be able to view the file
- Two files included, one .hp (raw scan voltage) and one .hpd (converted to magnetic units)
- Five header rows preceded with # character
- Data columns listed below for hp scan file:
- (1) X Position: Position of horizontal (X) motorized stage relative to HP origin
- (2) Y Position: Position of vertical (Y) motorized stage relative to HP origin
- (3) Z Position: Position of longitudinal (Z) motorized stage relative to HP origin
- (4) X Data: Voltage output of horizontal (X) Hall sensor
- (5) Y Data: Voltage output of vertical (Y) Hall sensor
- (6) Z Data: Voltage output of longitudinal (Z) Hall sensor
- (7) TM box: Temperature output of Hall probe electronics box
- (8) TM box: Temperature output of Hall probe
- (9) Time: RELATIVE time (seconds) counting from start of scan
- HPD processed scan file contains X/Y/Z position and BX/BY/BZ field data
  
- EOF
