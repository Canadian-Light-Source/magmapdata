#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// Import, Process, Graph Vibrating Wire data chapter of CLS Magnet Mapping Facility Analysis Macros
//
// Revision Date: January 2025
// Author: Cameron Baribeau
//
// Overview: The macros in this module handle vibrating wire frequency-scan data, including import,
//	processing, and B-field reconstruction from a set of multi-harmonic frequency sweeps.
//
// Vibrating wire scan file name structure: <deviceName>_<###mA>_<n##>_<xy##>.txt
//	The order of the drive current, wire harmonic, and x/y scan tags is not critical. However, it is
//	it is assumed that a device name is given at the start of the file name and that the wire drive
//	drive current, wire harmonic, and (for quad centering scans) an X/Y scan tag are all included.
//	Underscore "_" delimiters are necessary for the scripts to parse information from the scan name.
//	Multi-harmonic frequency sweeps do not require the X/Y scan tag.
//
// Future work:
//	- Improve usability (e.g. menu interface for global parameters), error handling, polish etc.
//
// List of macros in this chapter:
//	cma_vw_main
//	cma_vw_init
//	cma_vw_impDefaults
//	cma_vw_impFolder
//	cma_vw_impFile
//	cma_vw_params
//	cma_vw_export
//	cma_vwa_getSin
//	cma_vwa_checkSin
//	cma_vwa_getFw
//	cma_vwa_fitFw
//	cma_vwa_fitFw2
// 	cma_vwa_fitFxn
//	cma_vwa_fwFitMask
//	cma_vwa_refitFw
//	cma_vwa_getBn
//	cma_vwa_makeB
//	cma_vw_show
//	cma_vw_showPos
//	cma_vw_showFw
//	cma_vw_layout
//	cma_vw_results
//	cma_vw_plotBn
//	cma_vw_plotFund
//	cma_vw_plotN
//	cma_vw_showAP
//	cma_vw_modAP
//	cma_vw_showAA
//	cma_vw_modAA

//------------------------------------
// Main function that scripts calls to all key subroutines
Function cma_vw_main()
	String queryList = "No;Yes", commandStr
	Variable i, iStart, iEnd, keepImporting
	NVAR flag_skipInitImport
	Silent 1
	if((flag_skipInitImport == 0) || (exists("flag_skipInitImport") ==0)) 
		// Initialize vibrating wire parameters
		commandStr = "cma_vw_init()"
		Execute/Q commandStr
		Print "Vibrating wire globals initialized. Settings listed in cma_vw_init can be modified and the analysis re-ran..."
		// Declare / establish connection to global waves created in cma_vw_init
		Wave vwScanCurrents
		Wave/T vwLoadedScans
		// Begin file import loop
		do
			// Prompt user for import
			sprintf commandStr, "cma_vw_impFolder()"
			Execute/Q commandStr
			// Poll user: are there more folders to import?
			Prompt keepImporting, "More Folders to Import?", popup, queryList
			DoPrompt "Continue / Halt Data Import",keepImporting
		while(keepImporting == 2)
		Print "Data Import Complete. Analyzing data..."
	endif
	// Establish connection to global variables created in cma_vw_init
	NVAR sort_ByN, sort_ByStr, flag_printOut, flag_skipInitImport, flag_export
	// Also establish connection to global waves if this is not our first time through _main...
	if(flag_skipInitImport == 1)
		Wave vwScanCurrents
		Wave/T vwLoadedScans
	endif
	// Set flag to skip initialization and import on follow-up runs
	flag_skipInitImport = 1
	// -------------------
	// Process imported scans
	commandStr = "cma_vwa_getSin()"
	Execute commandStr
	commandStr = "cma_vwa_getFw()"
	Execute commandStr
	// Sort data
	Wave vwScanModes, vwScan_MaxAmps, vw_fundFreq, vwScan_SinNs, vw_Tavg
	if(sort_ByN == 1)
		Print "Sorting Data by mode number..."
		Sort vwScanModes, vwLoadedScans, vwScan_MaxAmps, vwScanCurrents, vwScanModes, vw_fundFreq, vwScan_SinNs, vw_Tavg
	endif
	if(sort_ByStr == 1)
		Print "Sorting Data by mode strength..."
		Sort/R vwScan_MaxAmps, vwLoadedScans, vwScan_MaxAmps, vwScanCurrents, vwScanModes, vw_fundFreq, vwScan_SinNs, vw_Tavg
	endif
	// Further analysis (everything else...)
	commandStr = "cma_vwa_fitFw()"
	Execute commandStr
	Print "Building Bn coefficients from F_w fit results..."
	commandStr = "cma_vwa_getBn()"
	Execute commandStr
	Print "Reconstructing Field..."
	commandStr = "cma_vwa_makeB(0, DimSize(vwScanModes,0))"
	Execute commandStr
	// Make summary table
	commandStr = "cma_vw_results()"
	Execute commandStr
	// Export results, if requirested
	if(flag_export == 1)
		commandStr = "cma_vw_export()"
		Execute commandStr
	endif
	Print "Data Analysis Complete! Results ready for graphing via cma_vw_show()!"
EndMacro

//------------------------------------
// Initialize global wave + string for vibrating wire data
Proc cma_vw_init()
	// Import default parameter values into vwDefaults
	cma_vw_impDefaults()
	Make/T/O/N=1 vwLoadedScans=""			// Wave containing names of all loaded vibrating wire scans
	Make/O/N=1 vwScanCurrents=0			// Wave containing drive current [A] for each imported wire scan
	Make/O/N=1 vwScanModes=0				// Wave containing harmonic number for each imported wire scan
	Make/O/N=1 vwScan_SinNs				// Wave containing relative amplitudes at position detector for each wire scan
	Make/O/N=1 vw_fundFreq					// Wave containing calculated frequency of n=1 harmoninc for each wire scan
	Make/O/N=1 vw_Tavg						// Wave containing average temperature for each wire scan
	DeletePoints 0,1, vwLoadedScans, vwScanCurrents, vwScanModes, vw_fundFreq, vwScan_SinNs, vw_Tavg
	// Declaration of global strings and variables...
	String /G deviceName = ""
	Variable /G photogateNodeCutoff = str2num(vwDefaults[8])	// Ignore modes where SIN( xSensor * n * Pi / wireLength) falls below this threshold
	Variable /G photgateCal1 = str2num(vwDefaults[9]) 		// Wire position detector amplitude [V] to wire displacement [um] (HRZ position detector)
	Variable /G photgateCal2 = str2num(vwDefaults[10]) 		// Wire position (VRT position detector)
	Variable /G wireLength = str2num(vwDefaults[11])		// [m]
	Variable /G xSensor_ch1 = str2num(vwDefaults[12])	// [m]
	Variable /G xSensor_ch2 = str2num(vwDefaults[13])	// [m]
	Variable /G wireMassPerL = str2num(vwDefaults[14])	// Wire mass density [kg/m]
	Variable /G phCleaner1 = str2num(vwDefaults[15])	// Threshold under which adjacent phase data are considered flat
	Variable /G phCleaner2 = str2num(vwDefaults[16])	// Threshold over which an isolated phase datum is considered an outlier, hence omitted from Fw fitting
	Variable /G phCleaner3 = str2num(vwDefaults[17])	// Threshold over which the SDEV of a phase measurement is considered too large, hence omitted from Fw fitting
	Variable /G tempResult = 0				// Global variable used to clumsily pass results between subroutines
	// Declaration of global user parameters (can be set via Prompt at execution time)
	Variable /G flag_skipInitImport = 0			// Flag to skip initialization and import (i.e. for re-running cma_vw_main in the same experiment)
	Variable /G flag_export = 1					// Flag to export results at end of cma_vw_main
	Variable /G flag_printOut = 1				// Flag to show some results at end of cma_vw_main
	Variable /G flag_old_files = str2num(vwDefaults[21])	// Flag to expect old (pre-2023) 1channel scan file formats
	Variable /G sort_ByStr = 0					// Flag to sort data according to mode contribution strength
	Variable /G sort_ByN = 1					// Flag to sort data according to mode number
	// F_w curve fitting mask settings
	Make/O/N=6 vw_fwFitParams
	vw_fwFitParams[0] = 0		// Binary flag to fit to (0) Temnykh's F_w or (1) wire amplitude function; *always* 0 for this work; references to this setting trimmed from code
	vw_fwFitParams[1] = 1		// Binary flag to (0) do nothing, (1) mask data far from Fw//Aw resonance from curve fitting
	vw_fwFitParams[2] = 0		// Binary flag to (0) do nothing; (1) mask data close to Fw/Aw resonance from curve fitting
	vw_fwFitParams[3] = 55		// Integer flag to (0) do nothing; (N) *stop* rejecting data too close to Fw/Aw resonance above N
	vw_fwFitParams[4] = str2num(vwDefaults[18])	// Reject points whose indices are beyond this threshold from the apparent wire resonance
	vw_fwFitParams[5] = str2num(vwDefaults[19])	// Reject points whose wire phase data are within this amount of 90deg
	// Hard coded error cutoff thresholds for [0] amplitude; [1,2] wire damping; and [3] frequency outlier
	Make/O/N=4 vw_fitErrCutoffs
	vw_fitErrCutoffs[0] = 50
	vw_fitErrCutoffs[1] = 0.1
	vw_fitErrCutoffs[2] = 100
	vw_fitErrCutoffs[3] = 0.15
	// Clean up vwDefaults wave, no longer needed
	KillWaves/Z vwDefaults
EndMacro

//------------------------------------
// Load default parameters from external file
// 2025-01-23: Rewrote function to prompt user for directory. 
// This file *can* be written with a hard-coded pathname if and when the init file is placed in a fixed directory
Proc cma_vw_impDefaults()
	String fileName = "CMA_v3_initVWdefaults.txt"
	// Create new path, prompting user for full directory
	Print "Provide directory to vibrating wire initialization file:"
	NewPath /Q/O pathName
	KillWaves/Z wave0, wave1					// Clear waves to be used for import
	Print "Attempting to load VW default parameters from selected directory..."
	// Load data
	LoadWave /J/N/K=2/Q  /P=pathName fileName
	// Test that incoming file has expected length
	if(DimSize(wave0, 0) != 24)
		Print "Warning! Imported default parameters file has an unexpected file length."
	endif
	// Clean up
	Make/O/T/N=(numPnts(wave1)) vwDefaults = wave1
	Print "VW default parameters loaded successfully. See wave vwDefaults."
	KillWaves/Z wave0, wave1
EndMacro

//------------------------------------
// Import and name all TXT files found inside a given directory, treating them as vibrating wire data
// 	N.B. Incoming data files have 20 columns and use/overwrite wave0,wave1...wave19
Proc cma_vw_impFolder()
	String fileName
	Variable i = 0
	Silent 1															  // Importing data files...
	// Create new path, prompting user for full directory
	NewPath /Q/O pathName
	// Make temporary wave for storing wavenames
	Make /O/T/N=1 tempName=""
	// Clear waves to be used for import
	KillWaves/Z wave0,wave1,wave2,wave3,wave4,wave5,wave6,wave7,wave8
	KillWaves/Z wave9,wave10,wave11,wave12,wave13,wave14,wave15,wave16,wave17,wave18,wave19
	do
		// Get next file name
		fileName = IndexedFile(pathName, i, "????")
		// Break loop when we find an empty file (or run out of files)
		if (strlen(fileName) == 0)
			break
		endif
		// Try to import every TXT file
		if (StringMatch(fileName, "*.txt"))
			cma_vw_impFile(fileName)
		endif				
		i += 1
	while(1)
	// Kill temporary wave
	KillWaves  tempName	
EndMacro

//------------------------------------
// Import vibrating wire file and rename waves, etc.
Proc cma_vw_impFile(fileName)
	String fileName, tstS
	String baseName, trialName
	String nextSrc, nextDest
	Variable i, loopCount, skipTrial = 0, skipRefMag = 0, nCols = 19
	// Record global device name (assumption is that this doesn't change in a single run of data)
	deviceName = StringFromList(0, fileName, "_")
	baseName = fileName
	// Expecting device name (magnet or insertion device) to be at front of file name before "_" delimiter
	if(StringMatch(baseName, "*_*") == 1)
		baseName = RemoveListItem(0, baseName, "_") 	// discard device name
	endif
	baseName = RemoveEnding(baseName, ".txt") 		// discard file extension
	// Load file
	LoadWave /J/N/Q /L={0,1,0,0,0} /P=pathName fileName
	if(skipTrial ==0)
		// Prep wave names (hard coded)
		Make/T/O/N=20 columnNames
		columnNames[0] = baseName + "_temp"
		columnNames[1] = baseName + "_freq"
		columnNames[2] = baseName + "_xpos"
		columnNames[3] = baseName + "_ypos"
		columnNames[4] = baseName + "_ch1_Xavg"
		columnNames[5] = baseName + "_ch1_Yavg"
		columnNames[6] = baseName + "_ch1_Ravg"
		columnNames[7] = baseName + "_ch1_Pavg"
		columnNames[8] = baseName + "_ch1_Xsdev"
		columnNames[9] = baseName + "_ch1_Ysdev"
		columnNames[10] = baseName + "_ch1_Rsdev"
		columnNames[11] = baseName + "_ch1_Psdev"
		columnNames[12] = baseName + "_ch2_Xavg"
		columnNames[13] = baseName + "_ch2_Yavg"
		columnNames[14] = baseName + "_ch2_Ravg"
		columnNames[15] = baseName + "_ch2_Pavg"
		columnNames[16] = baseName + "_ch2_Xsdev"
		columnNames[17] = baseName + "_ch2_Ysdev"
		columnNames[18] = baseName + "_ch2_Rsdev"
		columnNames[19] = baseName + "_ch2_Psdev"
		// Alternate column names for old file format (pre- second lock-in amplifier)
		if(flag_old_files == 1)
			nCols = 11
			columnNames[0] = baseName + "_temp"
			columnNames[1] = baseName + "_freq"
			columnNames[2] = baseName + "_ch1_Xavg"
			columnNames[3] = baseName + "_ch1_Yavg"
			columnNames[4] = baseName + "_ch1_Ravg"
			columnNames[5] = baseName + "_ch1_Pavg"
			columnNames[6] = baseName + "_ch1_Xsdev"
			columnNames[7] = baseName + "_ch1_Ysdev"
			columnNames[8] = baseName + "_ch1_Rsdev"
			columnNames[9] = baseName + "_ch1_Psdev"
			columnNames[10] = baseName + "_meter1avg"
			columnNames[11] = baseName + "_meter1sdev"
		endif
		// Rename scan data columns
		i = 0
		do
			nextSrc = "wave" + num2str(i)
			nextDest = columnNames[i]
			Rename $nextSrc, $nextDest
			i += 1
		while(i <= nCols)
		// Append to waves tracking name, drive current, harmonic number, average temperature of each imported scan
		Redimension /N=(DimSize(vwLoadedScans,0)+1)  vwLoadedScans, vwScanCurrents, vwScanModes, vw_fundFreq, vwScan_SinNs, vw_Tavg
		vwLoadedScans[DimSize(vwLoadedScans,0)-1] = baseName
		// Parse scan name for wire drive current
		if (StringMatch(baseName, "*mA_*"))
			i = 0
			// Loop through "_" delimiters in file name until the "###mA"
			do
				tstS = StringFromList(i, baseName, "_")
				if (StringMatch(tstS, "*mA*"))
					break
				endif
				i += 1
			while(1)
			// tstS should now be "###mA". Trim "mA" and save result
			vwScanCurrents[DimSize(vwLoadedScans,0)-1] = str2num(ReplaceString("mA", tstS, ""))
		endif
		// If wire current not found in scan name, put a default value
		if(StringMatch(baseName, "*mA_*") == 0)
			Print "Warning! No wire drive current found in scan name during import. File: " + fileName
			vwScanCurrents[DimSize(vwLoadedScans,0)-1] = 0.999
		endif
		// Parse scan name for wire harmonic number
		if(StringMatch(baseName, "*_n*"))
			i = 0
			// Loop through "_" delimiters in file name until the "n###"
			do
				tstS = StringFromList(i, baseName, "_")
				if (StringMatch(tstS, "n*"))
					break
				endif
				i += 1
			while(1)
			// tstS should now be "n###". Trim "n" and save result
			vwScanModes[DimSize(vwLoadedScans,0)-1] = str2num(ReplaceString("n", tstS, ""))
		endif
		// If mode number not found in scan name, put a default value
		if(StringMatch(baseName, "*_n*") == 0)
			Print "Warning! No wire harmonic number found in scan name during import. File: " + fileName
			vwScanModes[DimSize(vwLoadedScans,0)-1] = 99
		endif
		// Distill temperature wave down to an average value
		WaveStats /Q $(columnNames[0])
		vw_Tavg[DimSize(vwLoadedScans,0)-1] = V_avg
	endif
EndMacro

//------------------------------------
Proc cma_vw_export(exportDir)
	String exportDir
	Prompt exportDir, "Full Path Name for Export Directory of Fw Fitting Results"
	String vwTableName
	// Table name
	vwTableName = deviceName + "_FwFitResults"
	vwTableName = ReplaceString("-", vwTableName, "_")	// just in case
	// Build tables with given wave names
	Edit/N=$vwTableName vwLoadedScans, vwScanModes, vw_fundFreq, fitDat_a, fitDat_b, fitDat_c, vw_SinSeriesBn
	// Export path
	NewPath /Q/O exportPath, exportDir
	// Hardwire header info for export
	Make/T/O/N=1 VW_header = {"ScanName, Mode Num, n=1 Freq, Fw_Coef_a, Fw_Coef_b, Fw_Coef_c,SineSeries_Bn"}
	Save /O/G /M="\r\n" /P=exportPath VW_header as vwTableName + ".csv"
	// Append tables to headers
	SaveTableCopy /Z/N=0 /A=2 /T=2 /P=exportPath /W=$vwTableName as vwTableName + ".csv"
	// Cleanup
	KillWindow $vwTableName
	KillWaves/Z VW_header
EndMacro 

//------------------------------------
// Update and bring up key experiment parameters
// Note! vw_ParamV is only updated at runtime. If a user changes any global variable (e.g. wireLength) then vw_ParamV can fall out of sync
Proc cma_vw_params()
	PauseUpdate; Silent 1
	Make/O/T /N=5 vw_ParamT = {"Wire Length [m]","Wire Mass/L [m/kg]"," Detector Pos1 [m]","Detector Pos2 [m]", "Detector Gain1 [V/um]", "Detector Gain2 [V/um]"}
	Make /O/N=5 vw_ParamV = {wireLength, wireMassPerL, xSensor_ch1, xSensor_ch2, photGateCal1, photGateCal2}
	Edit/W=(24,382.25,281.25,557) vw_ParamT,vw_ParamV
	ModifyTable format(Point)=1,width(Point)=35,width(vw_ParamT)=110,width(vw_ParamV)=65
EndMacro

//------------------------------------
// Set wire harmonic sensitivity coefficients from position of detector (xSensor) and length of wire
Proc cma_vwa_getSin(chNum)
	Variable chNum
	Variable useXS = $("xSensor_ch" + num2str(chNum))
	vwScan_SinNs = SIN((Pi*vwScanModes/wireLength)*useXS)
EndMacro
//------------------------------------
// Check for modes where the photogate is sitting too close to a node
Proc cma_vwa_checkSin()
	Variable i, tst_badscan
	PauseUpdate; Silent 1
	// Prep node cutoff tracker
	Make/O/N=1 logNodeCutoffs
	DeletePoints 0, 1, logNodeCutoffs
	// Loop through wire harmonics to log bad fits
	i = 0
	do
		FindValue/V=(vwScanModes[i]) logBadScans
		tst_badscan = V_value	
		// Avoid double counting bad/empty scans
		if( !(Abs(vwScan_SinNs[i]) >= photogateNodeCutoff) && (tst_badscan == -1) )
			Redimension/N=(DimSize(logNodeCutoffs,0)+1) logNodeCutoffs
			logNodeCutoffs[DimSize(logNodeCutoffs,0)-1] = vwScanModes[i]
		endif
		i+=1
	while(i< DimSize(vwLoadedScans,0))
EndMacro

//------------------------------------
// Inputs: vibrating wire scan files
// Outputs: wire position, max normalized wire amplitude, adjusted wire phase all imported scans, F(w) waves
	// Channel 1: horizontal position detector
	// Channel 2: vertical position detector
Proc cma_vwa_getFw(chNum)
	Variable chNum
	PauseUpdate; Silent 1
	String baseName, nextMag, nextDMag, nextPhase, nextDPhase, nextFreq
	String dWire, dDWire, adjPhase, dAdjPhase, nextFw
	Variable i = 0, loopCount = DimSize(vwLoadedScans,0)
	Variable i2, loopCount2, tstA, tstB, boolTst, iTst, isnt_badscan
	Variable useCal = $("photgateCal" + num2str(chNum))
	// Prep wave for results
	Make/O/N=(loopCount) vwScan_MaxAmps = 0
	Make/O/N=1 tempWW
	// Prep log waves for rejected phase and Fw data
	Make/O/N=1 logBadPh_i, logBadFw_i
	Make/O/N=1/T logBadPh_scan, logBadFw_scan
	DeletePoints 0, 1, logBadPh_i, logBadPh_scan, logBadFw_i, logBadFw_scan
	// Loop through measurements
	do
		baseName = vwLoadedScans[i] + "_ch" + num2str(chNum)
		nextFreq = vwLoadedScans[i] + "_freq"
		// Get wave names
		nextMag = baseName + "_Ravg"
		nextDMag = baseName + "_Rsdev"
		nextPhase = baseName + "_Pavg"
		nextDPhase = baseName + "_Psdev"
		dWire = baseName + "_Uavg"
		dDWire = baseName +  "_Usdev"
		adjPhase = baseName + "_PavgAdj"
		dAdjPhase = baseName + "_PsdevAdj"
		nextFw = baseName + "_Fw"
		// Initialize results waves		
		Make/O/N=(DimSize($nextMag,0)) $dWire, $dDWire
		Make/O/N=(DimSize($nextPhase,0)) $adjPhase, $dAdjPhase
		Make/O/N=(DimSize($nextMag,0)) $nextFw
		Make/O/N=(DimSize($nextMag,0)) tempDeriv
		// ------
		// Convert photogate voltage [V] to wire displacement [um]
		$dWire = $nextMag / useCal
		$dDWire = $nextDMag / useCal
		// Convert from lock-in RMS reading to Pk reading
		$dWire *= Sqrt(2)
		$dDWire *= Sqrt(2)
		// Normalize wire displacement to drive current and sin(xSensor) (this isn't really used)
		Redimension /N=(DimSize($dWire,0)) tempWW
		tempWW = Abs($dWire) / (vwScanCurrents[i]* Abs(vwScan_SinNs[i]))
		WaveStats/Q tempWW
		vwScan_MaxAmps[i] = V_max
		// Where response is near a node, filter data to a small nonzero value (keep something to sort by)
		if(Abs(vwScan_SinNs[i]) < photogateNodeCutoff)
			vwScan_MaxAmps[i] /= 1000000
		endif
		// ------
		// Duplicate phase wave and adjust the twin while keeping the original for reference
		$adjPhase = $nextPhase
		$dAdjPhase = $nextDPhase
		loopCount2 = DimSize($adjPhase,0)
		// Toss suspicious phase data (spiky outlier points and/or massive SDEV)
		i2 = 1
		do
			// boolTst: test for outlier scenario
			boolTst = (Abs($adjPhase[i2-1] - $adjPhase[i2+1]) < phCleaner1) & (Abs($adjPhase[i2] - $adjPhase[i2+1]) > phCleaner2)
			if( boolTst || ($dadjPhase[i2] > phCleaner3) )
				$adjPhase[i2] = NaN
				$dAdjPhase[i2] = NaN
				// Increment and populate log waves for rejected phase data
				Redimension/N=(DimSize(logBadPh_i,0)+1) logBadPh_scan, logBadPh_i
				logBadPh_scan[DimSize(logBadPh_scan,0)-1] = baseName
				logBadPh_i[DimSize(logBadPh_i,0)-1] = i2
			endif
			i2 += 1
		while(i2 < loopCount2-1)
		// Check SDEV of first and final points
		i2 = 0
		do
			if( $dadjPhase[i2] > phCleaner3 )
				$adjPhase[i2] = NaN
				$dAdjPhase[i2] = NaN
				// Increment and populate log waves for rejected phase data
				Redimension/N=(DimSize(logBadPh_i,0)+1) logBadPh_scan, logBadPh_i
				logBadPh_scan[DimSize(logBadPh_scan,0)-1] = baseName
				logBadPh_i[DimSize(logBadPh_i,0)-1] = i2
			endif
			i2 += (loopCount2-1)
		while(i2 < loopCount2)
		// Unwrap phase to avoid discontinuities when "rolling over" +/- 180 phase (this is purely cosmetic)
		i2 = 2
		do
			iTst = 1
			// Compare present phase against last non-NaN data (or at worst the 1st index)
			tstA = $adjPhase[i2]
			do
				if (NumType($nextFw[i2-iTst]) != 2)
					break
				endif	
				iTst += 1
			while( (i2 - iTst) > 0)
			tstB = $adjPhase[i2 - iTst]
			// Offset by 2Pi any data point that has a step change, above threshold, from last non-NaN neighbour
			if(Abs(tstA - tstB) > 110)
				$adjPhase[i2] += 360*Sign($adjPhase[i2-iTst])
			endif
			i2 += 1
		while(i2 < loopCount2)	
		// Sometimes a scan begins near -180 and flows further negative (or vice versa). Offset whole scan if data >>360deg
		WaveStats/Q $adjPhase
		if(V_min < -540)
			$adjPhase += 360
		endif
		if(V_max > 540)
			$adjPhase -= 360
		endif
		// ------
		// F_w = Amplitude[m] * DriveCurrent[A]/2 * cos(phi[rad])  >>  (note, dWire is in [um])
		$nextFw = 0.5*$dWire*vwScanCurrents[i]*COS($adjPhase*Pi/180) * (0.000001)//[um]>>[m]
		// ------
		// Set scale/indexing of amplitude, phase, and Fw waves based on frequency sweep
		SetScale/P x $nextFreq[0], Abs($nextFreq[1]-$nextFreq[0]), "", $nextMag, $nextDMag, $nextPhase, $nextDPhase
		SetScale/P x $nextFreq[0], Abs($nextFreq[1]-$nextFreq[0]), "", $dWire, $dDwire, $adjPhase, $dAdjPhase, $nextFw
		i += 1
	while(i < loopCount)
	// Fw data cleaning, read flags and assign NaNs
	i = 0
	if( numPnts(logBadFw_i) > 0)
		do
			nextFw = logBadFw_scan[i] + "_Fw"
			i2 = logBadFw_i[i]
			$nextFw[i2] = NaN
			i += 1	
		while(i < DimSize(logBadFw_i,0))
	endif
	// One last pass to tag Fw containing overwhelming amount of NaNs
	i = 0
	Make/O/N=1 logBadScans
	DeletePoints 0, 1, logBadScans
	do
		// Wave names
		baseName = vwLoadedScans[i] + "_ch" + num2str(chNum)
		nextFw = baseName + "_Fw"
		// Check point count
		WaveStats/Q $nextFw
		isnt_badscan = V_npnts > 5
		if( !isnt_badscan)
			Redimension/N=(DimSize(logBadScans,0)+1) logBadScans
			logBadScans[DimSize(logBadScans,0)-1] = vwScanModes[i]
		endif
		i += 1
	while(i < loopCount)
	// Clean Up
	if( numPnts(logBadPh_i) > 0)
		Print "Warning: " + num2str(numPnts(logBadPh_i)) + " phase data points rejected for isolated spike; see logBadPh wave."
	endif
	if( numPnts(logBadFw_i) > 0)
		Print "Warning: " + num2str(numPnts(logBadFw_i)) + " F(w) data points rejected for isolated sign change; see logBadFw wave."
	endif
	if( numPnts(logBadScans) > 0)
		Print "Warning: " + num2str(numPnts(logBadScans)) + " F(w) data waves tagged as empty static; see logBadScans wave."
	endif	
	KillWaves/Z tempWW
	ResumeUpdate
EndMacro

//------------------------------------
// Perform curve fitting on all Fw waves of the selected signal channel (1 or 2)
// Save fitting coefficient results to waves fitDat_a, fitDat_b, fitDat_c, fitDat_d
Proc cma_vwa_fitFw(chNum)
	Variable chNum
	PauseUpdate; Silent 1
	Variable i=0, nScans = DimSize(vwLoadedScans,0)
	Variable tst_badscan, fwRMS, fwFitRMS, tst_badAmp, tst_badDamp
	String baseName, fwName, maskName, fName, phName, fitName
	// Prep coef summary waves
	Make/N=(nScans)/O fitDat_a, fitDat_b, fitDat_c, fitDat_d, fitDat_chiSq
	Make/N=(nScans)/O fitSigma_a, fitSigma_b, fitSigma_c, fitSigma_d
	Make/O tempFWmask
	// Hard coded fitting constraints for to "fitDat_c" aka K2, the damping term, and "fitDat_d" aka K3, the DC offset term
	Make/O/T/N=4 T_Constraints
	T_Constraints[0] = {"K2 > 0","K2 < 2.5","K3 > -1","K3 < 1"}
	Print "Fitting Fw curves..."
	do
		// Get wave names
		baseName = vwLoadedScans[i] + "_ch" + num2str(chNum)
		phName = baseName + "_PavgAdj"
		fwName = baseName + "_Fw"
		maskName = baseName + "_fitMask"
		fName = vwLoadedScans[i] + "_freq"
		fitName = "fit_" + fwName
		// Skip bad/empty scans
		FindValue/V=(vwScanModes[i]) logBadScans
		tst_badscan = V_value
		if( tst_badscan != -1 )
			// Make dummy/placeholder fit wave and write NaNs to results waves
			Make/O /N=2 $fitName=0
			cma_vwa_wipeFit(i)
		endif
		// Proceed with non bad scans
		if( tst_badscan == -1 )
			// Set fitting mask
			cma_vwa_fwFitMask(i, fwName, maskName, fName, phName)
			// Execute Fit
			cma_vwa_fitFw2($fwName, $fName, fitName)
			// Save fit and coefficients
			$fitName = tempFit
			fitDat_a[i] = W_coef[0]
			fitDat_b[i] = W_coef[1]
			fitDat_c[i] = W_coef[2]
			fitDat_d[i] = W_coef[3]
			// Save errors and chi square
			fitSigma_a[i] = W_sigma[0]
			fitSigma_b[i] = W_sigma[1]
			fitSigma_c[i] = W_sigma[2]
			fitSigma_d[i] = W_sigma[3]
			fitDat_chiSq[i] = tempResult
		endif // endif: empty scan filter
		i += 1
	while(i < nScans)
	// ------
	// Prep bad fit tracker
	Make/O/N=1 logBadFits
	DeletePoints 0, 1, logBadFits
	// Loop through wire harmonics
	i = 0
	do
		// Check if scan is in bad list
		FindValue/V=(vwScanModes[i]) logBadScans
		tst_badscan = V_value
		// Get wave names
		baseName = vwLoadedScans[i] + "_ch" + num2str(chNum)
		fwName = baseName + "_Fw"
		fitName = "fit_" + fwName
		// Establish basic logic blocks
		// Compare RMS of Fw and Fw_fit. If they are not in the same order of magnitude, the fit went wild
		WaveStats/Q $fwName
		fwRMS = V_rms
		WaveStats/Q $fitName
		fwFitRMS = V_rms
		// Boolean 1 indicates bad scan (either based on amplitude or damping)
		tst_badAmp = (Max(fwRMS, fwFitRMS) / Min(fwRMS, fwFitRMS)) > vw_fitErrCutoffs[0]
		tst_badDamp = !((Abs(fitDat_c[i]) < vw_fitErrCutoffs[2]) && (Abs(fitDat_c[i]) > vw_fitErrCutoffs[1]))
		// Log bad fits, but don't double-count bad/empty scans
		if( (tst_badscan == -1) && (tst_badAmp || tst_badDamp) )
			Redimension/N=(DimSize(logBadFits,0)+1) logBadFits
			logBadFits[DimSize(logBadFits,0)-1] = vwScanModes[i]
			cma_vwa_wipeFit(i)
		endif
		i+=1
	while(i< nScans)
	// Check flagged bad fits. Some can be restored by nudging the amplitude or frequency fit coeff's
	cma_vwa_refitFw(chNum)
	// Clean up
	KillWaves/Z tempFit, tempFWmask
	ResumeUpdate
EndMacro
//------------------------------------
// Helper function, depends upon cma_vwa_fitFw to set up the global fit mask and constraints
// Calling this function directly (without sitting up mask and constraints manually) will fail
Function cma_vwa_fitFw2(wavFw, wavFreq, sFit)
	Wave wavFw, wavFreq
	String sFit
	Wave tempFWmask, T_Constraints, tempDeriv, tempFit, vw_fwFitParams
	Variable initA, initB, initC, initD
	Variable V_FitError=0, V_FitQuitReason=0
	Variable fwRMS, fwFitRMS, fwSkew, fwFitTest
	NVAR tempResult
	// Establish baseline skewness of Fw near L/R curve edges. Used later to assess fits.
	Variable minFreq, maxFreq, freqL, freqR, fwSkewL, fwSkewR, fwFitSkewL, fwFitSkewR
	WaveStats/Q wavFw
	fwRMS = V_rms
	WaveStats/Q wavFreq
	minFreq = V_min
	maxFreq = V_max
	freqL = V_min + Abs(V_max-V_min)*0.3
	freqR = V_max - Abs(V_max-V_min)*0.3
	WaveStats/Q/R=(minFreq,freqL) wavFw
	fwSkewL = V_skew
	WaveStats/Q/R=(freqR,maxFreq) wavFw
	fwSkewR = V_skew
	// a: amplitude
	// b: resonant frequency
	// c: damping term ( ~1-2 E+00)
	// Initial guesses:
	initA = 0.0001
	initC = 1
	initD = 0
	Make/O/N=(DimSize(wavFw,0)) tempDeriv
	// Initial guess for frequency is location of largest derivative in F_w
	Differentiate wavFw /X = wavFreq /D = tempDeriv
	WaveStats/Q tempDeriv
	initB = wavFreq[V_maxRowLoc]
	if((Abs(V_min) > Abs(V_max)) || (V_maxRowLoc == 0) || (V_maxRowLoc == DimSize(wavFw, 0)))
		initB = wavFreq[V_minRowLoc]
	endif
	// Execute fit
	Make/D/N=4/O W_coef
	W_coef[0] = {initA, initB, initC, initD}
	V_FitError = 0
	V_FitQuitReason = 0
	FuncFit/Q/G/N/NTHR=0/W=2 cma_vwa_fitFxn W_coef wavFw /M = tempFWmask /X=wavFreq /D /C=T_Constraints
	// Pass local result into global variable for cma_vwa_fitFw
	tempResult = V_chisq
	// Connect string reference to fit wave, now that it exists
	Wave wavFit = $sFit
	Duplicate/O wavFit, tempFit
	// -----
	// Refine Fit. Initial fit can get stuck in a local minimum, so we try to nudge it free
	Wave W_sigma
	Variable i2, i3, tryA = 4, tryB = 6
	// Snapshot starting values
	Variable snapA, snapB, snapC, snapD, snapChi, snapFitTst, sigA, sigB, sigC, sigD
	snapA = W_coef[0]
	snapB = W_coef[1]
	snapC = W_coef[2]
	snapD = W_coef[3]
	sigA = W_sigma[0]
	sigB = W_sigma[1]
	sigC = W_sigma[2]
	sigD = W_sigma[3]
	snapChi = V_chisq
	WaveStats/Q wavFw
	fwFitRMS = V_rms
	WaveStats/Q/R=(minFreq,freqL) wavFw
	fwFitSkewL = V_skew
	WaveStats/Q/R=(freqR,maxFreq) wavFw
	fwFitSkewR = V_skew
	snapFitTst = 0.2*Abs(fwRMS - fwFitRMS)/fwRMS + 0.4*Abs(fwSkewL - fwFitSkewL)/Abs(fwSkewL) + 0.4*Abs(fwSkewR - fwFitSkewR)/Abs(fwSkewR)
	i2=0
	i3=0
	do
		i3=0
		do
			// Refine initial fitting guesses
			initA = snapA*(1 + 0.66*Ceil((tryA-1)/2)*((-1)^(tryA-1)))
			initB = snapB + 0.01*Ceil((tryB-1)/2)*((-1)^(tryB-1))
			initC = snapC
			initD = snapD
			W_coef[0] = {initA, initB, initC, initD}
			// Retry Fit
			FuncFit/Q/G/N/NTHR=0/W=2 cma_vwa_fitFxn W_coef wavFw /M = tempFWmask /X=wavFreq /D /C=T_Constraints
			tempResult = V_chisq
			// Take statistics to assess quality of fit
			WaveStats/Q wavFw
			fwFitRMS = V_rms
			WaveStats/Q/R=(minFreq,freqL) wavFw
			fwFitSkewL = V_skew
			WaveStats/Q/R=(freqR,maxFreq) wavFw
			fwFitSkewR = V_skew
			fwFitTest = 0.2*Abs(fwRMS - fwFitRMS)/fwRMS + 0.4*Abs(fwSkewL - fwFitSkewL)/Abs(fwSkewL) + 0.4*Abs(fwSkewR - fwFitSkewR)/Abs(fwSkewR)
			// Overwrite results if the fit quality is improved
			if(fwFitTest < snapFitTst)
				snapA = W_coef[0]
				snapB = W_coef[1]
				snapC = W_coef[2]
				snapD = W_coef[3]
				// Save errors and chi square
				sigA = W_sigma[0]
				sigB = W_sigma[1]
				sigC = W_sigma[2]
				sigD = W_sigma[3]
				snapChi = tempResult
				// Save the fit itself
				tempFit = wavFit
				snapFitTst = fwFitTest
			endif
			i3 += 1
		while(i3 < tryB)
		i2 +=1
	while(i2 < tryA)
	// Write best coefs/snaps back to W_ waves
	W_coef[0] = snapA
	W_coef[1] = snapB
	W_coef[2] = snapC
	W_coef[3] = snapD
	W_sigma[0] = sigA
	W_sigma[1] = sigB
	W_sigma[2] = sigC
	W_sigma[3] = sigD
	tempResult = snapChi	
	wavFit = tempFit
	// -----
	// Cleanup
	KillWaves/Z tempDeriv
EndMacro
//------------------------------------
Function cma_vwa_fitFxn(w0,w) : FitFunc
	Wave w0
	Variable w
	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(w) = a*(w-b)/(4*w*(w-b)^2+(w*c^2)) + d
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ w
	//CurveFitDialog/ Coefficients 4
	//CurveFitDialog/ w0[0] = a
	//CurveFitDialog/ w0[1] = b
	//CurveFitDialog/ w0[2] = c
	//CurveFitDialog/ w0[3] = d
	// f(w) = a*(w-b)/(4*w*(w-b)^2+(w*c^2)) + d
	//
	return w0[0]*(w-w0[1])/(4*w*(w-w0[1])^2+(w*w0[2]^2))+w0[3]
End
//------------------------------------
// Compare global settings stored in vw_fwFitParams against wire harmonic (harNum)
// Set fit mask for given F(w) curve
Proc cma_vwa_fwFitMask(hNum, fwName, fitMaskName, fName, phName)
	Variable hNum
	String fwName, fitMaskName, fName, phName
	Variable iMaskEdge1, iMaskEdge2
	// Initialize fitting mask
	Make/O /N=(DimSize($fwName,0)) $fitMaskName
	Redimension /N=(DimSize($fwName,0)) tempFWmask
	// Initialize temporary waves
	Make/O tempFWdat, tempMask1, tempMask2
	Redimension /N=(DimSize($fwName,0)) tempMask1, tempMask2
	tempMask1 = 1
	tempMask2 = 1
	// Choose fit mask methodology based on global settings and given wire harmonic 
	if(vw_fwFitParams[3] != 0)	
		vw_fwFitParams[2] = 1
		if(hNum >= vw_fwFitParams[3])
			vw_fwFitParams[2] = 0
		endif
	endif
	// Apply chosen fit mask methodology (if any)
	// Identify landmark to build fit filter around, using location of sharpest slope (we assume the derivative peaks on resonance)
	Duplicate/O $fwName, tempFWdat
	Differentiate $fwName /X = $fName /D = tempFWdat
	tempFWdat = Abs(tempFWdat)
	WaveStats/Q tempFWdat
	// Build fitting mask, Step 1 (reject data too far from resonance)
	if(vw_fwFitParams[1] == 1)
		iMaskEdge1 = V_maxRowLoc - vw_fwFitParams[4]
		iMaskEdge2 = V_maxRowLoc + vw_fwFitParams[4]
		tempMask1 = 0
		tempMask1[iMaskEdge1, iMaskEdge2] = 1
	endif
	// Build fitting mask, Step 2 (reject data too close to resonance)
	if(vw_fwFitParams[2] == 1)
		// Reject points within a threshold of 90deg (i.e. too close to resonance)
		tempMask2 = Limit(Floor(Abs(90 - Abs($phName))/vw_fwFitParams[5]),0,1)
	endif
	// Combine waves
	$fitMaskName = tempMask1*tempMask2
	// Send mask to global wave visible to cma_vwa_fitFw2
	tempFWmask = $fitMaskName
	// Clean up
	KillWaves/Z tempFWdat, tempMask1, tempMask2
End
//------------------------------------
// Wipe coefficients from a Fw fit identified as bad (NaN's make it obvious)
Proc cma_vwa_wipeFit(i)
	Variable i	
	fitDat_a[i] = NaN
	fitDat_b[i] = NaN
	fitDat_c[i] = NaN
	fitDat_d[i] = NaN
	fitSigma_a[i] = NaN
	fitSigma_b[i] = NaN
	fitSigma_c[i] = NaN
	fitSigma_d[i] = NaN
	fitDat_chiSq[i] = NaN
	vw_fundFreq[i] = NaN
EndMacro
//------------------------------------
// Review list of bad fits and attempt re-fitting using different initial guesses for amplitude and frequency
// Note: this function assumes cma_vwa_fitFw has already been called to set up T_Constraints
Proc cma_vwa_refitFw(chNum)
	Variable chNum
	PauseUpdate; Silent 1
	Variable i=0, iBad, keepTrying, iTryA, iTryB, maxTryA=7, maxTryB=10
	Variable freqGuess, initA, initB, initC, initD, fit_test
	Variable nBadFitsI, nBadFitsF, compHi, compLo
	Variable fwRMS, fwSkew, fwFitRMS, fwFitSkew, fwFitTest, fwFitBest, tst_badAmp, tst_badDamp
	Variable minFreq, maxFreq, freqL, freqR, fwSkewL, fwSkewR, fwFitSkewL, fwFitSkewR
	Variable V_FitError=0, V_FitQuitReason=0
	String baseName, fwName, maskName, fName, phName, fitName
	// Loop through flagged bad fits
	nBadFitsI = DimSize(logBadFits,0)
	nBadFitsF = nBadFitsI
	// Update wire frequencies
	vw_fundFreq = fitDat_b / vwScanModes
	// In the event that there are no bad fits: Great! Do nothing
	if(nBadFitsI >0)
		do
			// Get scan index for next bad fit
			FindValue/Z/V=(logBadFits[i]) vwScanModes
			iBad = V_value
			// Get wave names
			baseName = vwLoadedScans[iBad] + "_ch" + num2str(chNum)
			fwName = baseName + "_Fw"
			fName = vwLoadedScans[iBad] + "_freq"
			fitName = "fit_" + fwName
			maskName = baseName + "_fitMask"
			// All rejected/tagged scans are NaN'd, so the average fundFreq (un-polluted by wonky fit results) is a safe guess
			WaveStats/Q vw_fundFreq
			freqGuess = V_avg
			// Skip scans where too much Fw data has been rejected as NaN (bad scan)
			WaveStats/Q $fwName
			if(V_npnts > 5)
				// Establish baseline stats of Fw near L/R curve edges. Used to assess fits.
				WaveStats/Q $fName
				minFreq = V_min
				maxFreq = V_max
				freqL = V_min + Abs(V_max-V_min)*0.3
				freqR = V_max - Abs(V_max-V_min)*0.3
				WaveStats/Q/R=(minFreq,freqL) $fwName
				fwSkewL = V_skew
				WaveStats/Q/R=(freqR,maxFreq) $fwName
				fwSkewR = V_skew
				WaveStats/Q $fwName
				fwRMS = V_rms
				// Resize fitting mask to current scan
				Redimension /N=(DimSize($fwName,0)) tempFWmask
				// Assign fitting mask
				tempFWmask = $maskName
				// Enter Try Loop A (sweep amplitude)
				keepTrying = 1
				iTryA = 1
				do
					// Entry Try Loop B (sweep frequency)
					iTryB = 1
					do
						// Initial fitting guesses. Frequency (+)(-)(+)(-)... staircases away from resonant frequency
						initA = 0.0001*10^(iTryA-1)
						initB = (freqGuess + 0.02*Ceil((iTryB-1)/2)*((-1)^(iTryB-1)))*logBadFits[i]
						initC = 1
						initD = 0
						Make/D/N=4/O W_coef
						W_coef[0] = {initA, initB, initC, initD}
						// Retry fit (and reset error flags)
						V_FitError=0
						V_FitQuitReason=0
						FuncFit/Q/G/N/NTHR=0/W=2 cma_vwa_fitFxn W_coef $fwName /M = tempFWmask /X=$fName /D /C=T_Constraints
						tempResult = V_chisq
						// Take statistics to assess quality of fit
						WaveStats/Q $fitName
						fwFitRMS = V_rms
						WaveStats/Q/R=(minFreq,freqL) $fitName
						fwFitSkewL = V_skew
						WaveStats/Q/R=(freqR,maxFreq) $fitName
						fwFitSkewR = V_skew
						fwFitTest= 0.2*Abs(fwRMS - fwFitRMS)/fwRMS + 0.4*Abs(fwSkewL - fwFitSkewL)/Abs(fwSkewL) + 0.4*Abs(fwSkewR - fwFitSkewR)/Abs(fwSkewR)
						// Check whether this fit clears error thresholds -- all logic signals must be zero
						tst_badAmp = (Max(fwRMS, fwFitRMS) / Min(fwRMS, fwFitRMS)) > vw_fitErrCutoffs[0]
						tst_badDamp = !((Abs(W_coef[2]) < vw_fitErrCutoffs[2]) && (Abs(W_coef[2]) > vw_fitErrCutoffs[1]))
						//if((tst_badAmp == 0) && (tst_badDamp == 0) && (V_FitError == 0) && (V_FitQuitReason == 0))
						if( (tst_badAmp + tst_badDamp + V_fitError + V_FitQuitReason) == 0)
							// On first successful fit:
							if(keepTrying == 1)
								fitDat_a[iBad] = W_coef[0]
								fitDat_b[iBad] = W_coef[1]
								fitDat_c[iBad] = W_coef[2]
								fitDat_d[iBad] = W_coef[3]
								// Save errors and chi square
								fitSigma_a[iBad] = W_sigma[0]
								fitSigma_b[iBad] = W_sigma[1]
								fitSigma_c[iBad] = W_sigma[2]
								fitSigma_d[iBad] = W_sigma[3]
								fitDat_chiSq[iBad] = tempResult
								// Save the fit itself
								Duplicate/O $fitName, tempFit2
								fwFitBest = fwFitTest
								keepTrying = 0
							endif
							// Finish Try Loop B, peeking for better fits
							if(keepTrying == 0)
								// Overwrite results if the fit quality is improved
								if(fwFitTest < fwFitBest)
									fitDat_a[iBad] = W_coef[0]
									fitDat_b[iBad] = W_coef[1]
									fitDat_c[iBad] = W_coef[2]
									fitDat_d[iBad] = W_coef[3]
									// Save errors and chi square
									fitSigma_a[iBad] = W_sigma[0]
									fitSigma_b[iBad] = W_sigma[1]
									fitSigma_c[iBad] = W_sigma[2]
									fitSigma_d[iBad] = W_sigma[3]
									fitDat_chiSq[iBad] = tempResult
									// Save the fit itself
									tempFit2 = $fitName
									fwFitBest = fwFitTest
								endif
							endif
						endif
						iTryB += 1
					while(iTryB < maxTryB)
					iTryA += 1
				//while((keepTrying == 1) & (iTryA < maxTryA))
				while(iTryA < maxTryA)
				// After all attempts,if we found a good fit: restore best fit, remove scan from list of bad fits, decrement accordingly
				if(keepTrying == 0)
					DeletePoints i, 1, logBadFits
					i -= 1
					nBadFitsF -= 1
					$fitName = tempFit2
				endif
			endif
			i += 1
		while(i < nBadFitsF)
		// Bad fits still occasionally slip through. Final test: reject Fw fits whose calculated fund. freq stand out beyond a threshold
		vw_fundFreq = fitDat_b / vwScanModes
		// Loop until no hi/lo outliers remain
		keepTrying = 1
		do
			WaveStats/Q vw_fundFreq
			compHi = V_max - V_avg
			compLo = V_avg - V_min
			// Reject too large max
			if( compHi > vw_fitErrCutoffs[3] )
				Redimension/N=(DimSize(logBadFits,0)+1) logBadFits
				logBadFits[DimSize(logBadFits,0)-1] = vwScanModes[V_maxloc]
				cma_vwa_wipeFit(V_maxloc)
				nBadFitsF += 1
			endif
			// Reject too low min
			if( compLo > vw_fitErrCutoffs[3] )
				Redimension/N=(DimSize(logBadFits,0)+1) logBadFits
				logBadFits[DimSize(logBadFits,0)-1] = vwScanModes[V_minloc]
				cma_vwa_wipeFit(V_minloc)
				nBadFitsF += 1
			endif
			// Break condition
			if ((compHi < vw_fitErrCutoffs[3]) & (compLo < vw_fitErrCutoffs[3]))
				keepTrying = 0
			endif
		while(keepTrying)
		// Log results
		Print "Retrying fits on N = " + num2str(nBadFitsI) + " Fw curves identified as bad...N = " + num2str(nBadFitsF) + " Fw curves remain bad after retries. See logBadFits."
	endif
	KillWaves/Z tempFit2
EndMacro

//------------------------------------
// Calculate B_n cofficients of sine series representation of B-field across vibrating wire, applying some 
//	filtering for measurements that may produce bad/overwhelming contributions to the reconstructed field
Proc cma_vwa_getBn()
	Variable i = 0
	Variable testA, testB, testC
	PauseUpdate; Silent 1
	// Prep results wave
	Make/N=(DimSize(vwLoadedScans,0))/O vw_SinSeriesBn = 0
	// Re-check for any modes sitting on a wire node
	cma_vwa_checkSin()
	// Loop through wire harmonics
	do
		// Check if scan index is tagged as bad
		FindValue/V=(vwScanModes[i]) logBadFits
		testA = V_value
		FindValue/V=(vwScanModes[i]) logNodeCutoffs
		testB = V_value
		FindValue/V=(vwScanModes[i]) logBadScans
		testC = V_value
		// Reject Bn's listed in any of the above
		if( (testA == -1) && (testB == -1) && (testC == -1) )
			vw_SinSeriesBn[i] = fitDat_a[i] / vwScan_SinNs[i] * (2*wireMassPerL / vwScanCurrents[i]^2)
			// Extra normalization
			vw_SinSeriesBn[i] *= sqrt(pi)*wireLength
		endif
		i+=1
	while(i< DimSize(vwLoadedScans,0))
EndMacro

//------------------------------------
// Reconstruct B field across vibrating wire by summing calculated sin series B_n coefficients
Proc cma_vwa_makeB(iStart,iEnd)
	Variable iStart, iEnd
	PauseUpdate; Silent 1
	Variable nPts = wireLength*1000, i
	// Initialize intermediate and result waves
	Make/O/N=(nPts) rebuiltField = 0, tempSin, tempSinScaled
	SetScale /I x 0, wireLength, tempSin, tempSinScaled, rebuiltField
	i = iStart
	do
		// Calculate next sine contribution
		tempSin = sin(Pi*x * vwScanModes[i] / wireLength)
		tempSinScaled = tempSin * vw_SinSeriesBn[i]
		// Add next contribution to total
		rebuiltField += tempSinScaled
		i+=1
	while(i < iEnd)
	// Normalize field
	rebuiltField *= sqrt(2*Pi/wireLength)
	// Cleanup
	KillWaves/Z tempSin, tempSinScaled
EndMacro

//------------------------------------
// Build sequential plots showing measured wire position and calculated F_w curves for all wire harmonics
// between nStart and nEnd.
// 	chNum: lock-in amplifier channel number (1 or 2)
// 	nStart: First harmonic number to include in graphs
//	nEnd: Last harmonic number to include in graphs
//	gPerLayout: Number of plots on each 1-page Layout (e.g. 2, 3, 4, 6, 8...)
// 	savePic: [1/0] argument to export Layouts as PNG files
Proc cma_vw_show(chNum, nStart, nEnd, gPerLayout, savePic)
	Variable chNum nStart, nEnd, gPerLayout, savePic
	Prompt chNum, "Select Lock-In Amplifier Channel [1/2]"
	Prompt nStart, "First wire harmonic in sequence of graphs to be built" 
	Prompt nEnd, "Last wire harmonic in sequence of graphs to be built"
	Prompt gPerLayout, "Number of graphs to include in each 1-page layout"
	Prompt savePic, "Binary [0/1] option to export layout as a PNG file"
	// Declarations
	Variable i = 0, gStart, gEnd
	Variable gsRemaining = nEnd
	String readName1, readName2, writeName1, writeName2
	// If exporting then create new path, prompting user for full directory
	if(savePic == 1)
		NewPath /Q/O picPathName
	endif
	// Start from scratch by killing all existing Graphs and Layouts (yup, sorry)
	cma_clear_windows()
	// Make fresh set of graphs
	cma_vw_showPos(chNum, nStart, nEnd)
	cma_vw_showFw(chNum, nStart, nEnd)
	// Begin loop to build final layouts
	do
		gStart = 0 + i*(gPerLayout)
		if(gsRemaining > gPerLayout)
			gEnd = gStart + gPerLayout -1
		endif
		if(gsRemaining <= gPerLayout)
			gEnd = gStart + gsRemaining
		endif
		cma_vw_layout("WireAndPhase", gStart, gEnd)
		cma_vw_layout("Graph_Fw", gStart, gEnd)
		gsRemaining -= gPerLayout
		// Possibly export graphs
		if(savePic == 1)
			// Get layout names
			readName1 = "Layout" + num2str(2*i)
			readName2 = "Layout" + num2str(1 + 2*i)
			// Make file names
			writeName1 = "WireDispPlot_n" + num2str(gStart + nStart-1) + "to" + num2str(gEnd + nStart-1) + ".png"
			writeName2 = "FwPlot_n" + num2str(gStart + nStart-1) + "to" + num2str(gEnd + nStart-1) + ".png"
			// Execute file export
			SavePICT/P=picPathName/E=-5/B=288 /WIN=$readName1 as writeName1
			SavePICT/P=picPathName/E=-5/B=288 /WIN=$readName2 as writeName2
		endif		
		i += 1
	while(gsRemaining > 0)
EndMacro

//------------------------------------
// Create many graphs, each showing w. displacement and unwrapped phase data versus freq for one mode
Proc cma_vw_showPos(chNum, modeStart, modeEnd)
	Variable chNum, modeStart, modeEnd
	String baseName, fName, maskName, phiName, dphiName, wName, dwName, Cmd
	Variable i = modeStart, loopCount = modeEnd
	Variable axisHi, axisLo, useI
	PauseUpdate; Silent 1
	// Cheat to get graph names to work: create initial graph (which has no idex) as a blank slate, never used
	Display /N=WireAndPhase
	SetWindow kwTopWin hide=1
	Make/O/N=(DimSize(vwScanModes,0)) tempI
	do
		// Find index for next desired mode number
		tempI = Abs(vwScanModes - i)
		WaveStats/Q tempI
		useI = V_minloc
		baseName = vwLoadedScans[useI] + "_ch" + num2str(chNum)
		// Create graph, get names
		Display /N=WireAndPhase
		fName = vwLoadedScans[useI] + "_freq"
		maskName = baseName + "_fitMask"
		phiName = baseName + "_PavgAdj"
		dphiName = baseName + "_PsdevAdj"
		wName = baseName + "_Uavg"
		dwName = baseName + "_Usdev"
		// Append and format wire pos
		AppendToGraph $wName vs $fName
		ErrorBars $wName Y,wave=($dwName,$dwName)
		ModifyGraph rgb($wName) = (65280,0,0)
		// Append and format phase difference between wire pos and drive current
		AppendToGraph/R $phiName vs $fName
		ErrorBars $phiName Y,wave=($dphiName,$dphiName)
		ModifyGraph rgb($phiName) = (24576,24576,65280)
		// Textbox
		SprintF Cmd, "n = %s", num2str(vwScanModes[useI])
		TextBox/C/N=text0/A=LT/X=0.00/Y=0.00 Cmd
		// Formatting
		ModifyGraph mode=4,marker=8,lstyle=1,grid(left)=2
		Label left, "Wire Amplitude (um)"
		Label right "(xWire, I0) Phase Offset (deg)"
		Label bottom "Frequency (Hz)"
		// Hide Graph
		SetWindow kwTopWin hide=1
		i += 1
	while(i <= loopCount)
	// Clean Up
	KillWaves/Z tempI
EndMacro

//------------------------------------
// Create many graphs, each showing F_w data and curve fit versus freq for one mode
Proc cma_vw_showFw(chNum, modeStart, modeEnd)
	Variable chNum, modeStart, modeEnd
	String baseName, fwName, fName, fitName, maskName, Cmd
	Variable i = modeStart, loopCount = modeEnd, useI
	PauseUpdate; Silent 1
	// Cheat to get graph names to work: create initial graph (which has no idex) as a blank slate, never used
	Display /N=Graph_Fw
	SetWindow kwTopWin hide=1
	Make/O/N=(DimSize(vwScanModes,0)) tempI
	do
		// Find index for next desired mode number
		tempI = Abs(vwScanModes - i)
		WaveStats/Q tempI
		useI = V_minloc
		baseName = vwLoadedScans[useI] + "_ch" + num2str(chNum)
		// Create graph, get names
		Display /N=Graph_Fw
		fwName = baseName + "_Fw"
		fName = vwLoadedScans[useI] + "_freq"
		maskName = baseName + "_fitMask"
		fitName = "fit_" + fwName
		// Append and format Fw
		AppendToGraph $fwName vs $fName
		ModifyGraph mode($fwName) = 3, marker($fwName) = 8, msize($fwName)=2, rgb($fwName) = (0,0,0)
		AppendToGraph $fitName
		AppendToGraph/R $maskName vs $fName
		ModifyGraph lstyle($maskName)=2,lsize($maskName)=1.5, rgb($maskName)=(39168,39168,39168)
		// Textbox
		SprintF Cmd, "n = %s", num2str(vwScanModes[useI])
		TextBox/C/N=text0/A=LT/X=0.00/Y=0.00 Cmd
		// Graph formatting
		ModifyGraph grid(left)=2
		Label left "F_w"
		Label bottom "Frequency (Hz)"
		// Hide Graph
		SetWindow kwTopWin hide=1
		i += 1
	while(i <= loopCount)
	// Clean Up
	KillWaves/Z tempI
EndMacro

//------------------------------------
// Build a layout of already-created graphs 
//	(requires that they were created by the earlier cma_vw_show functions, which have procedurally determinable names)
Proc cma_vw_layout(baseName, iStart, iEnd)
	String baseName
	String graphName
	Variable iStart, iEnd
	Variable i = iStart
	PauseUpdate; Silent 1
	Layout/C=1/W=(5.25,43.25,373.5,475.25)
	do
		graphName = baseName + num2str(i)
		AppendToLayout $graphName
		i += 1
	while(i <= iEnd)
	Tile
EndMacro

//------------------------------------
Proc cma_vw_results()
	// Update n=1 frequencies
	vw_fundFreq = fitDat_b / vwScanModes
	// Make results table	
	Edit/W=(561.75,53,1419.75,247.25) vwLoadedScans,vwScanCurrents,vwScanModes,vwScan_MaxAmps
	AppendToTable vwScan_SinNs,vw_fundFreq, fitDat_a,fitDat_b,fitDat_c,fitDat_d,vw_SinSeriesBn
	ModifyTable format(Point)=1,width(Point)=30,width(vwScanCurrents)=82,width(vwScanModes)=68
	ModifyTable width(vw_fundFreq)=64, width(fitDat_a)=68,width(fitDat_b)=60,width(fitDat_c)=66,width(fitDat_d)=68
	// Tabulate logged snarls in data
	Edit/W=(824,78.25,1397,209.5) logBadPh_scan,logBadPh_i,logBadFw_scan,logBadFw_i
	AppendToTable logBadScans, logBadFits,logNodeCutoffs
	ModifyTable format(Point)=1,width(Point)=35,width(logBadPh_scan)=86,width(logBadPh_i)=62, width(logBadFw_scan)=81
	ModifyTable width(logBadFw_scan)=81,width(logBadFw_i)=62,width(logBadFits)=59,width(logNodeCutoffs)=77
	// Make results graphs
	cma_vw_plotBn()
	cma_vw_plotFund()
EndMacro

//------------------------------------
Window cma_vw_plotBn() : Graph
	PauseUpdate; Silent 1
	Display /W=(1026.75,278,1421.25,486.5) vw_SinSeriesBn vs vwScanModes
	ModifyGraph mode=8, marker=8, rgb=(0,0,65280), msize=2, grid=2
	Label left "Sine Series Stregnth (arb.) \\u#2"
	Label bottom "Wire Harmonic (n)"
	TextBox/C/N=text0/X=1.00/Y=1.00 "\\s(vw_SinSeriesBn) VW Scan"
EndMacro

//------------------------------------
Window cma_vw_plotFund() : Graph
	PauseUpdate; Silent 1
	Display /W=(1028.25,518,1422.75,701.75) vw_fundFreq vs vwScanModes
	ModifyGraph rgb=(39168,39168,39168), grid(left)=2
	Label left "Frequency of Wire 1st Harmonic (Hz)"
	Label bottom "Wire Harmonic (n)"
	WaveStats/Q vw_fundFreq
	SetAxis left Floor(V_min),Ceil(V_max)
EndMacro

// ---------
// Create detailed graphs for wire one harmonic
Proc cma_vw_plotN(chNum, useN)
	Variable chNum, useN
	Variable useI
	String baseName, fName, maskName, phiName, dphiName, wName, dwName, Cmd
	String fwName, fitName
	PauseUpdate; Silent 1
	// Find index for desired mode number
	Make/O/N=(DimSize(vwScanModes,0)) tempI
	tempI = Abs(vwScanModes - useN)
	WaveStats/Q tempI
	useI = V_minloc
	// Get wave names
	baseName = vwLoadedScans[useI] + "_ch" + num2str(chNum)
	fName = vwLoadedScans[useI] + "_freq"
	maskName = baseName + "_fitMask"
	phiName = baseName + "_PavgAdj"
	dphiName = baseName + "_PsdevAdj"
	wName = baseName + "_Uavg"
	dwName = baseName + "_Usdev"	
	fwName = baseName + "_Fw"
	fitName = "fit_" + fwName
	// First plot: wire displacement
	Display $wName vs $fName
	ErrorBars $wName Y,wave=($dwName,$dwName)
	ModifyGraph rgb($wName) = (65280,0,0)
	// Append and format phase difference between wire pos and drive current
	AppendToGraph/R $phiName vs $fName
	ErrorBars $phiName Y,wave=($dphiName,$dphiName)
	ModifyGraph rgb($phiName) = (24576,24576,65280)
	SprintF Cmd, "n = %s", num2str(vwScanModes[useI])
	TextBox/C/N=text0/A=LT/X=0.00/Y=0.00 Cmd
	// Formatting
	ModifyGraph mode=4,marker=8,lstyle=1,grid(left)=2
	Label left, "Wire Amplitude (um)"
	Label right "(xWire, I0) Phase Offset (deg)"
	Label bottom "Frequency (Hz)"
	// ---
	// Second plot: Fw
	Display $fwName vs $fName
	ModifyGraph mode($fwName) = 3, marker($fwName) = 8, msize($fwName)=2, rgb($fwName) = (0,0,0)
	AppendToGraph $fitName
	AppendToGraph/R $maskName vs $fName
	ModifyGraph lstyle($maskName)=2,lsize($maskName)=1.5, rgb($maskName)=(39168,39168,39168)
	// Textbox
	SprintF Cmd, "n = %s", num2str(vwScanModes[useI])
	TextBox/C/N=text0/A=LT/X=0.00/Y=0.00 Cmd
	// Graph formatting
	ModifyGraph grid(left)=2
	Label left "F_w"
	Label bottom "Frequency (Hz)"
	// Clean Up
	KillWaves/Z tempI
EndMacro


// ---------
// Plot wire amplitude and phase from one X/Y channel across a list of harmonics in the wave named by given string.
// The named wave (e.g. useList = "wave0") must contain a list of integers corresponding to the wire
// 	harmonics to be plotted, e.g. wave0 = [1,2,3,4,5,11,101]
Proc cma_vw_showAP(chNum, useList)
	Variable chNum
	String useList
	PauseUpdate; Silent 1		// building window...
	Display
	Variable i = 0, nPts = DimSize($useList,0)
	String nextWire, nextF, nextN, nextPh
	do
		nextN = num2str($useList[i])
		if($useList[i] < 10)
			nextN = "0" + nextN
		endif
		nextWire = "fscan_n" + nextN + "_ch" + num2str(chNum) + "_Uavg"
		nextPh =  "fscan_n" + nextN + "_ch" + num2str(chNum) + "_PavgAdj"
		nextF = "fscan_n" + nextN + "_freq"
		AppendToGraph $nextWire vs $nextF
		AppendToGraph/R $nextPh vs $nextF
		ModifyGraph rgb($nextPh) = (24576,24576,65280)
		i += 1
	while(i<nPts)
	ModifyGraph mode=4, marker=8, lStyle=1, grid(left)=2, grid(bottom)=2
	// Call subroutine for further modifications to graph (trace offsets etc)
	cma_vw_modAP(chNum, useList)
EndMacro
// ---------
Proc cma_vw_modAP(chNum, useList)
	Variable chNum
	String useList
	PauseUpdate; Silent 1		// building window...
	String nextN, nextScan, nextWire, nextPh
	Variable i = 0, nPts = DimSize($useList,0), nextFreq
	Variable ampInit, ampPk, ampMul, ampOffset
	// Loop through waves
	do
		nextN = num2str($useList[i])
		if($useList[i] < 10)
			nextN = "0" + nextN
		endif
		nextScan = "fscan_n" + nextN
		nextWire = "fscan_n" + nextN + "_ch" + num2str(chNum) + "_Uavg"
		nextPh =  "fscan_n" + nextN + "_ch" + num2str(chNum) + "_PavgAdj"
		FindValue /TEXT=nextScan vwLoadedScans
		if(V_value == -1)
			Print "Warning! Expected scan not found in vwLoadedScans!"
		endif
		nextFreq = fitDat_b[V_value]		
		// Offset and scale traces to put the minimum value at 0.1 and the max at 1.0 arbitrary units
		WaveStats/Q $nextWire
		ampPk = V_max
		ampInit = V_min
		ampOffset = (ampPk - 10*ampInit)/9
		ampMul = 1/(ampPk + ampOffset)
		ampOffset *= ampMul
		ModifyGraph muloffset($nextWire)={0,ampMul}
		ModifyGraph offset($nextWire)={-nextFreq,ampOffset}
		ModifyGraph offset($nextPh)={-nextFreq,0}
		i += 1
	while(i < nPts)
	Label left "Wire Amplitude (r.u.)";DelayUpdate
	Label bottom "Frequency about Resonance (Hz)";DelayUpdate
	Label right "Wire Phase (deg)"
	ModifyGraph mode=0
	SetAxis left 0,1
EndMacro


// ---------
// Plot wire amplitudes from both X/Y channels across a list of harmonics in the wave named by given string.
// The named wave (e.g. useList = "wave0") must contain a list of integers corresponding to the wire
// 	harmonics to be plotted, e.g. wave0 = [1,2,3,4,5,11,101]
Proc cma_vw_showAA(useList)
	String useList
	PauseUpdate; Silent 1		// building window...
	Display
	Variable i = 0, nPts = DimSize($useList,0)
	String nextWire, nextF, nextN, nextPh, nextDWire
	do
		nextN = num2str($useList[i])
		if($useList[i] < 10)
			nextN = "0" + nextN
		endif
		// Channel 1
		nextWire = "fscan_n" + nextN + "_ch1_Uavg"
		nextPh =  "fscan_n" + nextN + "_ch1_PavgAdj"
		nextF = "fscan_n" + nextN + "_freq"
		nextDWire = "fscan_n" + nextN + "_ch1_Usdev"
		AppendToGraph $nextWire vs $nextF
		ErrorBars $nextWire Y,wave=($nextDWire,$nextDWire)
		ModifyGraph lstyle($nextWire)=1
		// ---
		// Channel 2
		nextWire = "fscan_n" + nextN + "_ch2_Uavg"
		nextPh =  "fscan_n" + nextN + "_ch2_PavgAdj"
		nextF = "fscan_n" + nextN + "_freq"
		nextDWire = "fscan_n" + nextN + "_ch2_Usdev"
		AppendToGraph/R $nextWire vs $nextF
		ErrorBars $nextWire Y,wave=($nextDWire,$nextDWire)
		ModifyGraph rgb($nextWire) = (19456,39168,0)
		ModifyGraph lstyle($nextWire)=3
		i += 1
	while(i<nPts)
	ModifyGraph mode=4, marker=8, grid(left)=2, grid(bottom)=2
	// Call subroutine for further modifications to graph (trace offsets etc)
	cma_vw_modAA(useList)
EndMacro
// ---------
// Companion macro to cma_showAA: applies offsets to renormalize curves
Proc cma_vw_modAA(useList)
	String useList
	PauseUpdate; Silent 1		// building window...
	String nextN, nextScan, nextWire, nextPh
	Variable i = 0, nPts = DimSize($useList,0), nextFreq
	Variable ampInit, ampPk, ampMul, ampOffset
	// Loop through waves
	do
		nextN = num2str($useList[i])
		if($useList[i] < 10)
			nextN = "0" + nextN
		endif
		// Channel 1
		nextScan = "fscan_n" + nextN
		nextWire = "fscan_n" + nextN + "_ch1_Uavg"
		nextPh =  "fscan_n" + nextN + "_ch1_PavgAdj"
		FindValue /TEXT=nextScan vwLoadedScans
		if(V_value == -1)
			Print "Warning! Expected scan not found in vwLoadedScans!"
		endif
		nextFreq = fitDat_b[V_value]
		// Offsets
		WaveStats/Q $nextWire
		ampPk = V_max
		ampInit = V_min
		ampOffset = (ampPk - 10*ampInit)/9
		ampMul = 1/(ampPk + ampOffset)
		ampOffset *= ampMul
		ModifyGraph muloffset($nextWire)={0,ampMul}
		ModifyGraph offset($nextWire)={-nextFreq,ampOffset}
		// Channel 2
		nextWire = "fscan_n" + nextN + "_ch2_Uavg"
		nextPh =  "fscan_n" + nextN + "_ch2_PavgAdj"
		FindValue /TEXT=nextScan vwLoadedScans
		if(V_value == -1)
			Print "Warning! Expected scan not found in vwLoadedScans!"
		endif
		nextFreq = fitDat_b[V_value]
		// Offsets
		WaveStats/Q $nextWire
		ampPk = V_max
		ampInit = V_min
		ampOffset = (ampPk - 10*ampInit)/9
		ampMul = 1/(ampPk + ampOffset)
		ampOffset *= ampMul
		ModifyGraph muloffset($nextWire)={0,ampMul}
		ModifyGraph offset($nextWire)={-nextFreq,ampOffset}
		i += 1
	while(i < nPts)
	Label left "Wire Amplitude (r.u.)"
	Label bottom "Frequency about Resonance (Hz)"
	ModifyGraph mode=0
	SetAxis left 0,1
EndMacro


//------------------------------------
// EOF
