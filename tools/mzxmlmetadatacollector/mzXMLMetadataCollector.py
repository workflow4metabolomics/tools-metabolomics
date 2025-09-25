# -*- coding: utf-8 -*-
"""
Created on Dec 04 2023
@author: quruin

Last update on Nov 14 2024

@author: quruin
"""
# Packages
import os
import re
import sys

import numpy as np


def get_info(file, format1):
    # mzml files #
    if (format1.casefold() == 'mzml'):
        # Saving all info of first file until first scan
        text1 = read_mzml_file(file)

        # Getting version format
        version1 = get_version(text1)

        # Getting ko size
        taillek1 = get_kb_size(file)

        # Getting Mo size
        taillem1 = get_mb_size(file)

        # Getting Go size
        tailleg1 = get_gb_size(file)

        # Getting MS level
        mslevel1 = get_ms_level(text1)

        # Getting spectrum type
        spectrum1 = get_spectrum_type(text1)

        # Getting sourc file
        source1 = get_source_file(text1)

        # Getting acquisition date
        date1 = get_acquisition_date(text1)

        # Getting software used
        softwaresList1 = get_softwares_list(text1)

        # Getting processing methods
        ProcessList1 = get_processing_list(text1)

        # Getting machine model
        modele1 = get_instrument_model(text1)

        # Getting points and scans number
        nbscans1 = count_scans(text1)
        nbpoints1 = count_points(file)

        # Encoding
        encodage1 = get_encoding(text1)

        return format1, version1, taillek1, taillem1, tailleg1, mslevel1, \
            spectrum1, source1, date1, softwaresList1, ProcessList1, \
            modele1, nbscans1, nbpoints1, encodage1
    
    # mzXml files #
    if (format1.casefold() == 'mzxml'):
        # Saving all info of first file until first scan
        text1 = read_mzxml_file(file)

        # Getting format version
        version1 = get_mzxml_version(text1)

        # Getting ko size
        taillek1 = get_kb_size(file)

        # Getting Mo size
        taillem1 = get_mb_size(file)

        # Getting Go size
        tailleg1 = get_gb_size(file)

        # Getting MS level
        mslevel1 = "Not available in mzXML"

        # Getting spectrum type
        spectrum1 = get_mzxml_spectrum_type(text1)

        # Getting source file
        source1 = get_source_mzxml(text1)

        # Getting acquisition date
        date1 = "Not available in mzXML"

        # Getting used softwares
        softwaresList1 = get_mzxml_software_list(text1)

        # Processing
        ProcessList1 = "Not available in mzXML"

        # Getting Machine model
        modele1 = get_mzxml_instrument_model(text1)

        # Getting scans and points number
        nbscans1 = count_scans_mzxml(text1)
        nbpoints1 = count_points_mzxml(file)

        # Encoding
        encodage1 = get_encoding_mzxml(text1)

        return format1, version1, taillek1, taillem1, tailleg1, \
            mslevel1, spectrum1, source1, date1, softwaresList1, \
            ProcessList1, \
            modele1, nbscans1, nbpoints1, encodage1
    
    # Other types of files
    ######################################################
    if ((format1.casefold() != 'mzml') & (format1.casefold() != 'mzxml')):
        return ['' for i in range(0, 15)]

def get_encoding_mzxml(text1):
    try:
        encodage1 = text1.split('precision="')[1].split('"')[0] + "-bit"
    except Exception:
        encodage1 = 'Not available'
    return encodage1

def count_points_mzxml(file):
    try:
        nbpoints1 = 0
        with open(file, 'r') as f:
            line = 'start'
            while (line != ''):
                line = f.readline()
                if ('peaksCount="' in line):
                    nbpoints1 = nbpoints1 + \
                            int(line.split('peaksCount="')[1].split('"')[0])
    except Exception:
        nbpoints1 = "Calculation failed"
    return nbpoints1

def count_scans_mzxml(text1):
    try:
        nbscans1 = text1.split('msRun scanCount="')[1].split('"')[0]
    except Exception:
        nbscans1 = "Not found"
    return nbscans1

def get_mzxml_spectrum_type(text1):
    try:
        spectrum1 = text1.split('centroided="')[1].split('"')[0]
        if spectrum1 == "0":
            spectrum1 = "profile"
        else:
            if spectrum1 == "1":
                spectrum1 = "centroid"
            else:
                spectrum1 = "Not found"
    except Exception:
        spectrum1 = "Not found"
    return spectrum1

def get_source_mzxml(text1):
    try:
        source1 = text1.split('parentFile fileName="')[1].split('"')[0]
    except Exception:
        source1 = "Not found"
    return source1

def get_mzxml_software_list(text1):
    try:
        softwaresList1 = ''
        for i in range(text1.count("<software")):
            if i != 0:
                softwaresList1 = softwaresList1 + ' + '
            softwaresList1 = softwaresList1 + \
                    text1.split("<software")[i+1].split('type="')[1]\
                    .split('"')[0] + ': ' +\
                    text1.split("<software")[i+1].split('name="')[1]\
                    .split('"')[0] + ' ' +\
                    text1.split("<software")[i+1]\
                    .split('version="')[1]\
                    .split('"')[0]
    except Exception:
        softwaresList1 = "Not found"
    return softwaresList1

def get_mzxml_instrument_model(text1):
    try:
        if ('<msModel category="msModel" value="' in text1):
            modele1 = text1.split('<msModel category="msModel" \
                    value="')[1].split('"')[0]
        else:
            if ('accession="MS:1001547"' in text1):
                modele1 = text1.split('<cvParam cvRef="MS" \
                        accession="MS:1001547" name="')[1].split('"')[0]
            else:
                if ('accession="MS:1003123"' in text1):
                    modele1 = text1.split('<cvParam cvRef="MS" \
                            accession="MS:1003123" name="')[1].split('"')[0]
                else:
                    if ('accession="MS:1000495"' in text1):
                        modele1 = text1.split('<cvParam cvRef="MS" \
                                accession="MS:1000495" name="')[1]\
                                    .split('"')[0]
                    else:
                        modele1 = "Not found"
    except Exception:
        modele1 = "Not found"
    return modele1

def get_mzxml_version(text1):
    try:
        version1 = text1.split('<mzXML')[1].split('mzXML_')[1]\
                .split('"')[0].split(' ')[0]
    except Exception:
        version1 = "Not found"
    return version1

def read_mzxml_file(file):
    f = open(file, 'r+', encoding="utf-8")
    with f:
        text1 = ''
        while ('</scan>' not in text1):
            text1 = text1 + f.readline()
    f.close()
    return text1

def read_mzml_file(file):
    with open(file, 'r+', encoding="utf-8") as f:
        text1 = ''
        while ('<binary>' not in text1):
            text1 = text1 + f.readline()
    return text1

def get_gb_size(file):
    try:
        tailleg1 = np.round(os.stat(file).st_size/1024**3, 2)
    except Exception:
        tailleg1 = "Calculation failed"
    return tailleg1

def get_mb_size(file):
    try:
        taillem1 = np.round(os.stat(file).st_size/1024**2, 2)
    except Exception:
        taillem1 = "Calculation failed"
    return taillem1

def get_kb_size(file):
    try:
        taillek1 = np.round(os.stat(file).st_size/1024, 2)
    except Exception:
        taillek1 = "Calculation failed"
    return taillek1

def get_version(text1):
    try:
        version1 = text1.split('<mzML')[1] \
                .split('version="')[1].split('"')[0]
    except Exception:
        version1 = "Not found"
    return version1

def get_spectrum_type(text1):
    try:
        if ('<cvParam cvRef="MS" accession="MS:1000127" \
                    name="' in text1):
            spectrum1 = text1.split('<cvParam cvRef="MS"  \
                    accession="MS:1000127" name="')[1].split('"')[0]
        else:
            if ('<cvParam cvRef="MS" accession="MS:1000127" \
                        value="" name="' in text1):
                spectrum1 = text1.split('<cvParam cvRef="MS" \
                        accession="MS:1000127" value="" name="')[1]\
                            .split('"')[0]
            else:
                if ('<cvParam accession="MS:1000127" \
                            cvRef="MS" name="' in text1):
                    spectrum1 = text1.split('<cvParam \
                            accession="MS:1000127" \
                                cvRef="MS" name="')[1].split('"')[0]
                else:
                    spectrum1 = "Not found"
    except Exception:
        spectrum1 = "Not found"
    return spectrum1

def get_source_file(text1):
    try:
        source1 = text1.split('<sourceFile ')[1].split('">')[0]\
                .split('location="')[1].split('"')[0]
    except Exception:
        try:
            source1 = text1.split('<sourceFile id="WIFF" name="')[1]\
                    .split('"')[0]
        except Exception:
            try:
                source1 = text1.split('<sourceFile id=" ')[1].split('"')[0]
            except Exception:
                source1 = "Not found"
    return source1

def get_acquisition_date(text1):
    try:
        date1 = text1.split('startTimeStamp="')[1].split('"')[0]\
                .split('T')[0] + " " + text1.split('startTimeStamp="')[1]\
                .split('"')[0].split('T')[1]
    except Exception:
        date1 = "Not found"
    return date1

def get_softwares_list(text1):
    try:
        softwaresList1 = ''
        subtext1 = text1.split('softwareList count')[1] \
                .split('</softwareList>')[0]
        for i in range(subtext1.count("<software")):
            if i != 0:
                softwaresList1 = softwaresList1 + ' + '
            softwaresversions1 = subtext1.split("<software")[i+1]\
                    .split('version="')[1].split('"')[0]
            if ('cvRef="MS" name="' in subtext1):
                softwares1 = subtext1.split('cvRef="MS" name="')[i+1]\
                        .split('"')[0]
            else:
                softwares1 = subtext1.split('<software id="')[i+1]\
                        .split('"')[0]
            softwaresList1 = softwaresList1 + \
                    softwares1 + ' ' + softwaresversions1
    except Exception:
        softwaresList1 = "Not found"
    return softwaresList1

def get_processing_list(text1):
    try:
        ProcessList1 = ''
        subtext1 = text1.split('<dataProcessingList count')[1] \
                .split('</dataProcessingList>')[0]
        for i in range(subtext1.count("<processingMethod ")):
            if i != 0:
                ProcessList1 = ProcessList1 + ' + '
            methods1 = subtext1.split("<processingMethod ")[i+1] \
                    .split('name="')[1].split('"')[0]
            softwares1 = subtext1.split('softwareRef="')[1] \
                    .split('"')[0]
            ProcessList1 = ProcessList1 + softwares1 + ' ' + methods1
    except Exception:
        ProcessList1 = "Not found"
    return ProcessList1

def get_encoding(text1):
    encodage1 = re.search('name="(.*)-bit float"', text1).group(0).replace('name="', "").replace('"', '')
    return encodage1

def count_points(file):
    try:
        nbpoints1 = 0
        with open(file, 'r') as f:
            line = 'start'
            while (line != ''):
                line = f.readline()
                if 'defaultArrayLength=' in line:
                    nbpoints1 = nbpoints1 + int(
                            line.split('defaultArrayLength="')[1]
                            .split('"')[0])
        f.close()
    except Exception:
        nbpoints1 = "Calculation failed"
    return nbpoints1

def count_scans(text1):
    try:
        nbscans1 = text1.split('<spectrumList count="')[1].split('"')[0]
    except Exception:
        nbscans1 = "Not found"
    return nbscans1

def get_instrument_model(text1):
    try:
        if ('<cvParam cvRef="MS" accession="MS:1000703" name="' in text1):
            modele1 = text1.split('<cvParam cvRef="MS" \
                    accession="MS:1000703" name="')[1].split('"')[0]
        else:
            if ('accession="MS:1001547"' in text1):
                modele1 = text1.split('<cvParam cvRef="MS" \
                        accession="MS:1001547" name="')[1].split('"')[0]
            else:
                if ('<cvParam cvRef="MS" \
                            accession="MS:1003094" name="' in text1):
                    modele1 = text1.split('<cvParam cvRef="MS" \
                            accession="MS:1003094" name="')[1].split('"')[0]
                else:
                    if ('accession="MS:1003123" name="' in text1):
                        modele1 = text1.split('<cvParam cvRef="MS" \
                                accession="MS:1003123" name="')[1] \
                                    .split('"')[0]
                    else:
                        if ('<cvParam cvRef="MS" \
                                    accession="MS:1000495"' in text1):
                            modele1 = text1.split('<cvParam cvRef="MS" \
                                    accession="MS:1000495" name="')[1]\
                                        .split('value="')[1].split('"')[0]
                        else:
                            if ('<cvParam accession="MS:1000703" \
                                        cvRef="MS" name="' in text1):
                                modele1 = text1.split('<cvParam \
                                        accession="MS:1000703" \
                                            cvRef="MS" name="')[1] \
                                                .split('"')[0]
                            else:
                                if ('<cvParam cvRef="MS" accession= \
                                        "MS:1000483" \
                                            value="" name="' in text1):
                                    modele1 = text1.split('<cvParam \
                                            cvRef="MS" accession="MS:1000483" \
                                                value="" name="')[1] \
                                                    .split('"')[0]
                                else:
                                    if ('<cvParam accession="MS:1000703" \
                                                cvRef="MS" name="' in text1):
                                        modele1 = text1.split('<cvParam \
                                                accession="MS:1000703" \
                                                    cvRef="MS" \
                                                        name="')[1] \
                                                            .split('"')[0]
                                    else:
                                        modele1 = 'Not available'
        if (modele1 == ''):
            modele1 = 'Not available'
    except Exception:
        modele1 = "Not found"
    return modele1

def get_ms_level(text1):
    try:
        if ('<cvParam cvRef="MS" accession="MS:1000580" name="' in text1):
            return text1.split('<cvParam cvRef="MS" accession="MS:1000580" name="')[1].split('"')[0]
        
        if ('<cvParam cvRef="MS" accession="MS:1000580" value="" name="' in text1):
            return text1.split('<cvParam cvRef="MS" accession="MS:1000580" value="" name="')[1].split('"')[0]

        if ('accession="MS:1000579" cvRef="MS" name="' in text1):
            return text1.split('accession="MS:1000579" cvRef="MS" name="')[1].split('"')[0]

        if ('cvRef="MS" accession="MS:1000579" name="' in text1):
            return text1.split('cvRef="MS" accession="MS:1000579" name="')[1].split('"')[0]
        if (' <cvParam cvRef="MS" accession="MS:1000580" name="' in text1):
            return text1.split(' <cvParam cvRef="MS" accession="MS:1000580" name="')[1].split('"')[0]
        return "Not found"
    except Exception:
        mslevel1 = 'Error'
    return mslevel1


def update_entries(infile, singleFile):
    if (singleFile == "SINGLE"):
        info = ''
        info += 'Name\tformat\tversion\tsize(ko)\tsize(Mo)\tsize(Go)\
        \tMSlevel\tSpectrum type\tSource file\tAcquisition \
        date\tSoftware(s) used\tProcessing method(s)\tMachine\t \
        Number of scans\tNumber of points\tEncoding\n'
        x = get_info(infile, sys.argv[4].split('.')[-1])
        info += sys.argv[4] + '\t'
        for i in range(len(x)):
            info += str(x[i]) + '\t'
        info += '\n'
        return info
    else:
        info = ''
        info += 'Name\tformat\tversion\tsize(ko)\tsize(Mo)\tsize(Go)\t \
        MSlevel\tSpectrum type\tSource file\tAcquisition \
        date\tSoftware(s) used\tProcessing method(s)\tMachine\t \
        Number of scans\tNumber of points\tEncoding\n'
        ii = 4
        for f in infile.split(','):
            x = get_info(f, sys.argv[ii].split('.')[-1])
            info += sys.argv[ii] + '\t'
            for i in range(len(x)):
                info += str(x[i]) + '\t'
            info += '\n'
            ii += 1
        return info


print([i for i in sys.argv])
outfile = sys.argv[1]
singleFile = sys.argv[2]
infile = sys.argv[3]
res = open(outfile, "w")
res.write(update_entries(infile, singleFile))
res.close()
